#!/usr/bin/env Rscript

# Usage (examples):
#   Rscript transform_chrUnk_base_to_insilico.R \
#     --insilico_vcf insilicovcf.vcf \
#     --base_vcf base_vcf \
#     --Outputfilename out_vcf \
#
# Notes:
#     --insilico_vcf : The final reference anchored insilico references vcf which has gone through all iterations of insilico
#     --base_vcf : The final reference anchored vcf for all samples.
#     --Outputfilename : The name of the transfored reference anchored vcf for all samples \



args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (grepl("^--", a)) {
      key <- sub("^--", "", a)
      if (grepl("=", key, fixed = TRUE)) {
        kv <- strsplit(key, "=", fixed = TRUE)[[1]]
        out[[kv[1]]] <- kv[2]
      } else {
        if (i < length(args) && !grepl("^--", args[i + 1])) {
          out[[key]] <- args[i + 1]
          i <- i + 1L
        } else {
          out[[key]] <- TRUE
        }
      }
    }
    i <- i + 1L
  }
  out
}

opts <- parse_args(args)

get_opt <- function(x, key, default = NULL) {
  if (!is.null(x[[key]])) x[[key]] else default
}

insilico_vcf <- get_opt(opts, "insilico_vcf", "")
base_vcf     <- get_opt(opts, "base_vcf", "")
out_vcf      <- get_opt(opts, "Outputfilename", "")

if (!nzchar(insilico_vcf) || !nzchar(base_vcf) || !nzchar(out_vcf)) {
  stop(
    "Usage:\n",
    "  Rscript transform_chrUnk_base_to_insilico_streaming.R \\\n",
    "    --insilico_vcf insilico.vcf[.gz] \\\n",
    "    --base_vcf base.vcf[.gz] \\\n",
    "    --Outputfilename vcftransformed.vcf[.gz] \n"
  )
}

# -----------------------------
# File helpers
# -----------------------------
open_text_in <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }
}

open_text_out <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, "wt")
  } else {
    file(path, "wt")
  }
}



# DNA helpers (derived from main anchoring script)
revcomp_iupac <- function(x) {
  map <- c(
    A="T", C="G", G="C", T="A",
    M="K", R="Y", W="W", S="S", Y="R", K="M",
    V="B", H="D", D="H", B="V", N="N"
  )

  rc_one <- function(s) {
    if (is.na(s) || s == "") return(NA_character_)
    s <- toupper(trimws(as.character(s)))
    if (!nzchar(s)) return(NA_character_)
    ch <- strsplit(s, "", fixed = TRUE)[[1]]
    rc <- rev(ifelse(ch %in% names(map), map[ch], ch))
    paste0(rc, collapse = "")
  }

  vapply(x, rc_one, character(1))
}

norm_base <- function(x) {
  x <- toupper(trimws(as.character(x)))
  if (length(x) == 0L) return(NA_character_)
  if (is.na(x) || x == "") return(NA_character_)
  x
}

same_base <- function(a, b) {
  a <- norm_base(a)
  b <- norm_base(b)
  !is.na(a) && !is.na(b) && a == b
}

# Flip GT 0<->1 only in the GT portion, preserving anything after ":"
# Examples:
#   0/0:35   -> 1/1:35
#   1/1:35   -> 0/0:35
#   0/1:35   -> 0/1:35    (normalized for unphased hets)
#   1/0:35   -> 0/1:35    (normalized for unphased hets)
#   0|1:35   -> 1|0:35    (phasing preserved)
flip_gt_01 <- function(vec, normalize_unphased_het = TRUE) {
  x <- as.character(vec)
  x[is.na(x)] <- ""

  gt   <- sub("^([^:]+).*", "\\1", x)
  rest <- sub("^[^:]*", "", x)

  flip_one <- function(g) {
    if (g %in% c("", ".", "./.", ".|.")) {
      return(if (g == "") "" else "./.")
    }

    phased <- grepl("\\|", g)
    sep <- if (phased) "|" else "/"
    alleles <- strsplit(g, "[/|]", perl = TRUE)[[1]]

    out <- vapply(alleles, function(a) {
      if (a == ".") return(".")
      if (a == "0") return("1")
      if (a == "1") return("0")
      # leave other allele indices unchanged
      a
    }, character(1))

    # normalize unphased het 1/0 -> 0/1
    if (!phased && normalize_unphased_het && length(out) == 2L) {
      if (setequal(out, c("0", "1"))) out <- c("0", "1")
    }

    paste(out, collapse = sep)
  }

  gt2 <- vapply(gt, flip_one, character(1))
  paste0(gt2, rest)
}

# -----------------------------
# VCF header reader
# -----------------------------
read_vcf_header <- function(path) {
  con <- open_text_in(path)
  on.exit(close(con), add = TRUE)

  meta <- character(0)
  header <- NULL

  repeat {
    ln <- readLines(con, n = 1L, warn = FALSE)
    if (!length(ln)) break
    ln <- sub("^\ufeff", "", ln, perl = TRUE)

    if (startsWith(ln, "##")) {
      meta <- c(meta, ln)
    } else if (startsWith(ln, "#")) {
      header <- ln
      break
    } else {
      stop("Malformed VCF: missing #CHROM header line")
    }
  }

  if (is.null(header)) stop("Malformed VCF: missing #CHROM header line")

  list(meta = meta, header = header)
}

# -----------------------------
# Build insilico chrUnk lookup by ID
# Stores only first biallelic chrUnk record per ID
# -----------------------------
build_insilico_lookup <- function(path) {
  message("Building insilico chrUnk lookup from: ", path)

  H <- read_vcf_header(path)
  cols <- strsplit(sub("^#", "", H$header), "\t", fixed = TRUE)[[1]]

  required <- c("CHROM", "ID", "REF", "ALT")
  if (!all(required %in% cols)) {
    stop("Insilico VCF missing required columns: ", paste(setdiff(required, cols), collapse = ", "))
  }

  ix_chrom <- match("CHROM", cols)
  ix_id    <- match("ID", cols)
  ix_ref   <- match("REF", cols)
  ix_alt   <- match("ALT", cols)

  con <- open_text_in(path)
  on.exit(close(con), add = TRUE)

  # skip header
  repeat {
    ln <- readLines(con, n = 1L, warn = FALSE)
    if (!length(ln)) stop("Malformed VCF while re-reading insilico header")
    ln2 <- sub("^\ufeff", "", ln, perl = TRUE)
    if (startsWith(ln2, "#") && !startsWith(ln2, "##")) break
  }

  lookup <- new.env(hash = TRUE, parent = emptyenv())
  n_chrunk <- 0L
  n_stored <- 0L

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    if (!nzchar(line)) next

    f <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(f) < max(ix_chrom, ix_id, ix_ref, ix_alt)) next

    chrom <- f[ix_chrom]
    id    <- f[ix_id]
    ref   <- f[ix_ref]
    alt   <- f[ix_alt]

    if (!same_base(chrom, "chrUnk")) next
    n_chrunk <- n_chrunk + 1L

    if (is.na(id) || id == "" || id == ".") next
    if (grepl(",", alt, fixed = TRUE)) next  # biallelic only
    if (exists(id, envir = lookup, inherits = FALSE)) next

    lookup[[id]] <- c(REF = norm_base(ref), ALT = norm_base(alt))
    n_stored <- n_stored + 1L
  }

  message("  chrUnk records seen in insilico: ", n_chrunk)
  message("  biallelic chrUnk IDs stored:     ", n_stored)

  list(
    header = H,
    cols = cols,
    lookup = lookup,
    n_stored = n_stored
  )
}

# -----------------------------
# Main transform
# -----------------------------
transform_base_vcf <- function(ins_lookup, base_vcf, out_vcf) {
  message("Reading base header: ", base_vcf)

  H <- read_vcf_header(base_vcf)
  cols <- strsplit(sub("^#", "", H$header), "\t", fixed = TRUE)[[1]]

  required <- c("CHROM", "POS", "ID", "REF", "ALT")
  if (!all(required %in% cols)) {
    stop("Base VCF missing required columns: ", paste(setdiff(required, cols), collapse = ", "))
  }

  ix_chrom <- match("CHROM", cols)
  ix_pos   <- match("POS", cols)
  ix_id    <- match("ID", cols)
  ix_ref   <- match("REF", cols)
  ix_alt   <- match("ALT", cols)
  ix_fmt   <- match("FORMAT", cols)

  sample_idx <- if (!is.na(ix_fmt) && ix_fmt < length(cols)) {
    (ix_fmt + 1L):length(cols)
  } else {
    integer(0)
  }


  message("Writing output to: ", out_vcf)

  in_con  <- open_text_in(base_vcf)
  out_con <- open_text_out(out_vcf)
  on.exit(try(close(in_con), silent = TRUE), add = TRUE)
  on.exit(try(close(out_con), silent = TRUE), add = TRUE)

  # Copy base header exactly as-is
  repeat {
    ln <- readLines(in_con, n = 1L, warn = FALSE)
    if (!length(ln)) stop("Unexpected EOF while copying base header")
    writeLines(ln, out_con)

    ln2 <- sub("^\ufeff", "", ln, perl = TRUE)
    if (startsWith(ln2, "#") && !startsWith(ln2, "##")) break
  }

  stats <- c(
    exact_match = 0L,
    swapped_flip = 0L,
    revcomp_ref_no_flip = 0L,
    revcomp_alt_flip = 0L,
    not_in_insilico = 0L,
    unresolved = 0L,
    unchanged_non_chrunk = 0L,
    unchanged_multiallelic = 0L,
    total_records = 0L
  )

  audit_path <- paste0(out_vcf, ".audit")
  audit_con <- file(audit_path, open = "wt")
  on.exit(try(close(audit_con), silent = TRUE), add = TRUE)

  writeLines(
    paste(
      c(
        "ID",
        "base_REF_before",
        "base_ALT_before",
        "insilico_REF",
        "insilico_ALT",
        "rc_base_REF",
        "rc_base_ALT",
        "action"
      ),
      collapse = "\t"
    ),
    audit_con
  )

  repeat {
    line <- readLines(in_con, n = 1L, warn = FALSE)
    if (!length(line)) break

    stats["total_records"] <- stats["total_records"] + 1L

    if (!nzchar(line)) {
      writeLines(line, out_con)
      next
    }

    f <- strsplit(line, "\t", fixed = TRUE)[[1]]

    # malformed line -> pass through unchanged
    if (length(f) < max(ix_chrom, ix_id, ix_ref, ix_alt)) {
      writeLines(line, out_con)
      next
    }

    chrom <- f[ix_chrom]
    id    <- f[ix_id]
    bref  <- f[ix_ref]
    balt  <- f[ix_alt]

    bref_n <- norm_base(bref)
    balt_n <- norm_base(balt)

    # Only target chrUnk rows in base
    if (!same_base(chrom, "chrUnk")) {
      stats["unchanged_non_chrunk"] <- stats["unchanged_non_chrunk"] + 1L
      writeLines(line, out_con)
      next
    }

    # Only biallelic
    if (grepl(",", balt_n, fixed = TRUE)) {
      stats["unchanged_multiallelic"] <- stats["unchanged_multiallelic"] + 1L
      writeLines(line, out_con)
      next
    }

    # Must exist in insilico lookup
    if (is.na(id) || id == "" || id == "." || !exists(id, envir = ins_lookup, inherits = FALSE)) {
      stats["not_in_insilico"] <- stats["not_in_insilico"] + 1L
      writeLines(line, out_con)
      next
    }

    ins <- ins_lookup[[id]]
    iref <- norm_base(ins[["REF"]])
    ialt <- norm_base(ins[["ALT"]])

    rc_bref <- revcomp_iupac(bref_n)
    rc_balt <- revcomp_iupac(balt_n)

    action <- NULL
    flip_needed <- FALSE

    # 1) exact match
    if (same_base(bref_n, iref) && same_base(balt_n, ialt)) {
      stats["exact_match"] <- stats["exact_match"] + 1L
      action <- "exact_match"

      writeLines(line, out_con)
      writeLines(
        paste(id, bref_n, balt_n, iref, ialt, rc_bref, rc_balt, action, sep = "\t"),
        audit_con
      )
      next
    }

    # 2) swapped REF/ALT -> flip GT
    if (same_base(bref_n, ialt) && same_base(balt_n, iref)) {
      action <- "swapped_flip"
      flip_needed <- TRUE
      stats["swapped_flip"] <- stats["swapped_flip"] + 1L

    # 3) insilico REF == reverse complement of base REF -> no GT flip
    } else if (same_base(rc_bref, iref)) {
      action <- "revcomp_ref_no_flip"
      flip_needed <- FALSE
      stats["revcomp_ref_no_flip"] <- stats["revcomp_ref_no_flip"] + 1L

    # 4) insilico REF == reverse complement of base ALT -> flip GT
    } else if (same_base(rc_balt, iref)) {
      action <- "revcomp_alt_flip"
      flip_needed <- TRUE
      stats["revcomp_alt_flip"] <- stats["revcomp_alt_flip"] + 1L

    # unresolved
    } else {
      stats["unresolved"] <- stats["unresolved"] + 1L
      action <- "unresolved"

      writeLines(line, out_con)
      writeLines(
        paste(id, bref_n, balt_n, iref, ialt, rc_bref, rc_balt, action, sep = "\t"),
        audit_con
      )
      next
    }

    # Apply allele changes to match insilico
    f[ix_ref] <- iref
    f[ix_alt] <- ialt

    # Flip genotypes if required
    if (flip_needed && length(sample_idx)) {
      f[sample_idx] <- flip_gt_01(f[sample_idx], normalize_unphased_het = TRUE)
    }

    new_line <- paste(f, collapse = "\t")
    writeLines(new_line, out_con)

    writeLines(
      paste(id, bref_n, balt_n, iref, ialt, rc_bref, rc_balt, action, sep = "\t"),
      audit_con
    )
  }

  close(in_con)
  close(out_con)
  close(audit_con)

  summary_path <- paste0(out_vcf, ".summary")
  summary_df <- data.frame(
    metric = names(stats),
    count = as.integer(stats),
    stringsAsFactors = FALSE
  )

  utils::write.table(
    summary_df,
    file = summary_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  message("Done.")
  message("Summary:")
  for (nm in names(stats)) {
    message(sprintf("  %-24s %d", nm, stats[[nm]]))
  }
  message("Audit file:   ", audit_path)
  message("Summary file: ", summary_path)

  invisible(list(
    output = out_vcf,
    audit = audit_path,
    summary = summary_path,
    stats = stats
  ))
}

# -----------------------------
# Run
# -----------------------------
ins <- build_insilico_lookup(insilico_vcf)
transform_base_vcf(ins$lookup, base_vcf, out_vcf)
