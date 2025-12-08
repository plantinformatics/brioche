#!/usr/bin/env Rscript
# 
# Convert unanchored genotypes to an anchored VCF
#
# Usage (examples):
#   Rscript Anchoring_script.R \
#     --Rawgenotypes Chickpea_ALL_RawGenotypes_filtered.namesTrimmed.tsv \
#     --Briochemappings AVR_Ta_Hv_40K_v1_2_with_cicar.CDCFrontier.gnm3.QT0PBrioche_all_markers1to1stagingforvcf.csv \
#     --IsVCFraws false \
#     --IsSNPchip true \
#     --isDArTfile false \
#     --dobiallelic false \
#     --Outputfilename my_output.vcf \
#     --Accession "" \
#     --genomename ""
#
# Notes:
#   - --Rawgenotypes : path to the raw genotype table/vcf file. Must have minimum values of Name REF, ALT, and a genotype matrix of data 
#   - --Briochemappings : path to Brioche mappings CSV *Brioche_all_markers1to1stagingforvcf
#   - --IsVCFraws : if TRUE, sample raw genotype file is already in a vcf format
#   - --IsSNPchip : if TRUE, data being input is from an ilumina SNP chip
#   - --isDArTfile false : If TRUE, sample raw genotype file is in DArT SNP format. File must be 1 row format from Dart e.g., 1 allele with genotypes coded as 0,1,2 (as opposed to two row where there are two rows for the same allele with one representing the ref and one the alt and genotypes are coded as 1 0 binary
#   - --dobiallelic false : If TRUE brioche format REF and ALT have been transformed to match only the REF/ALT of the targets file, no RC required in anchoring script.
#   - --Outputfilename : path/filename to write VCF (default: output.vcf)
#   - --Accession : Accession of the reference genome used for mapping in Brioche e.g., GCF_000331145.2 (don't do GCA)  If provided, metadata from NCBI will be incorporated into vcf header info
#   - --genomename :genome name (title) of the reference genome used for mapping in Brioche e.g., Cicar.CDCFrontier_v2.0  If provided, metadata from NCBI will be incorporated into vcf header info. If Accession and genomename are provided, Accession is used first then genomename afterwards if no metadata was found using accession



# install packages if they are missing
if (!requireNamespace("dplyr", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("dplyr", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
if (!requireNamespace("getopt", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("getopt", repos = "https://cloud.r-project.org", quiet = TRUE)))
}

if (!requireNamespace("stringr", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("stringr", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
if (!requireNamespace("vcfR", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("vcfR", repos = "https://cloud.r-project.org", quiet = TRUE)))
}

if (!requireNamespace("rentrez", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("rentrez", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
if (!requireNamespace("readr", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("readr", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("readr", repos = "https://cloud.r-project.org", quiet = TRUE)))
}


suppressPackageStartupMessages({
  library(getopt)
  library(dplyr)
  library(stringr)
  library(vcfR)
  library(rentrez)
  library(readr)
  library(data.table)
})


args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (grepl("^--", a)) {
      a2 <- sub("^--", "", a)
      if (grepl("=", a2)) {
        kv <- strsplit(a2, "=", fixed = TRUE)[[1]]
        out[[kv[1]]] <- kv[2]
      } else {
        if (i < length(args) && !grepl("^--", args[i + 1])) {
          out[[a2]] <- args[i + 1]
          i <- i + 1L
        } else {
          out[[a2]] <- "true" 
        }
      }
    }
    i <- i + 1L
  }
  out
}

opts <- parse_args(args)

get_opt <- function(opts, key, default = NULL) {
  if (!is.null(opts[[key]])) opts[[key]] else default
}

# Parse common truthy strings to logical
to_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)) %in% c("1","true","t","yes","y")
}

raw_path      <- get_opt(opts, "Rawgenotypes", "Chickpea_ALL_RawGenotypes_filtered.namesTrimmed.tsv")
map_path      <- get_opt(opts, "Briochemappings", "AVR_Ta_Hv_40K_v1_2_with_cicar.CDCFrontier.gnm3.QT0PBrioche_all_markers1to1stagingforvcf.csv")
is_vcfraws    <- to_bool(get_opt(opts, "IsVCFraws", "false"))
is_snpchip    <- to_bool(get_opt(opts, "IsSNPchip", "true"))
isDArTfile    <- to_bool(get_opt(opts, "isDArTfile", "false"))
dobiallelic   <- to_bool(get_opt(opts, "dobiallelic", "false")) 
Outputfilename<- get_opt(opts, "Outputfilename", "outfile")
Accession     <- get_opt(opts, "Accession",  "")
genomename    <- get_opt(opts, "genomename", "")

message(sprintf("Rawgenotypes   : %s", raw_path))
message(sprintf("Briochemappings: %s", map_path))
message(sprintf("IsVCFraws      : %s", is_vcfraws))
message(sprintf("IsSNPchip      : %s", is_snpchip))
message(sprintf("isDArTfile     : %s", isDArTfile))
message(sprintf("dobiallelic   : %s", dobiallelic))
message(sprintf("Outputfilename : %s", Outputfilename))
message(sprintf("Accession      : %s", Accession))
message(sprintf("genomename     : %s", genomename))

# Functions

nzchar_safe <- function(x) !is.na(x) & x != "" & x != "."

# ALT CSV utilities
split_alts <- function(x) {
  x <- as.character(x)
  lapply(x, function(s) if (is.na(s) || s == "" || s == ".") character(0) else strsplit(s, ",", fixed = TRUE)[[1]])
}
join_alts <- function(v) if (!length(v)) "." else paste(v, collapse = ",")

# Build 0/1/2… index remap given old (REF + ALTs) and new (REF + ALTs)
# Returns an integer vector mapping old indices -> new indices (0-based).
build_symbol_maps <- function(ref_old, alts_old_chr, ref_new, alts_new_chr) {
  alts_old_chr <- if (is.null(alts_old_chr)) character(0) else alts_old_chr
  alts_new_chr <- if (is.null(alts_new_chr)) character(0) else alts_new_chr
  sym_old <- c(ref_old, alts_old_chr)
  sym_new <- c(ref_new, alts_new_chr)
  idx_map <- integer(length(sym_old))
  for (i in seq_along(sym_old)) {
    m <- match(toupper(sym_old[i]), toupper(sym_new))
    idx_map[i] <- ifelse(is.na(m), NA_integer_, m - 1L)  # 0-based
  }
  idx_map
}

upsert_info_key_vec <- function(x, key, val) {
  x <- as.character(x); x[is.na(x)] <- ""
  stopifnot(length(x) == length(val))
  out <- x
  pat <- paste0("(^|;)", key, "=[^;]*")
  for (i in seq_along(x)) {
    xi <- x[i]; vi <- val[i]
    if (!nzchar(vi)) vi <- ""  # allow empty if ever needed
    has <- grepl(paste0("(^|;)", key, "="), xi, perl = TRUE)
    if (has) {
      out[i] <- sub(pat, paste0("\\1", key, "=", vi), xi, perl = TRUE)
    } else {
      out[i] <- if (xi == "") paste0(key, "=", vi) else paste0(xi, ";", key, "=", vi)
    }
  }
  out
}

# Recode a vector of GT fields using an index map; preserves phasing and trailing subfields
# vec like "0/1:12" or "1|0" or "./." ; idx_map is 0-based integer vector from build_symbol_maps(...)
recode_gt_vec <- function(vec, idx_map) {
  v <- as.character(vec)
  v[is.na(v)] <- ""
  if (!length(idx_map)) return(v)

  gt   <- sub("^([^:]+).*", "\\1", v)
  rest <- sub("^[^:]*", "", v)     # includes leading ":" if present

  is_single <- gt %in% c("0","1","2","3","4")
  gt[is_single] <- paste0(gt[is_single], "/", gt[is_single])

  recode_one <- function(g) {
    if (g == "" || g == "." || g == "./.") return(if (g == "") "" else "./.")
    sep <- if (grepl("\\|", g)) "|" else "/"
    a  <- strsplit(g, "[/|]", perl = TRUE)[[1]]
    map_allele <- function(x) {
      if (x == ".") return(".")
      xi <- suppressWarnings(as.integer(x))
      if (is.na(xi) || xi + 1L > length(idx_map)) return(".")
      y <- idx_map[xi + 1L]
      if (is.na(y)) "." else as.character(y)
    }
    a2 <- vapply(a, map_allele, character(1))
    if (any(a2 == ".")) return("./.")
    paste(a2, collapse = sep)
  }

  out_gt <- vapply(gt, recode_one, character(1))
  out_gt[out_gt == "."] <- "./."
  paste0(out_gt, rest)
}

# INFO parsers/upserters (vectorised)
get_info_key <- function(x, key) {
  x <- as.character(x); x[is.na(x)] <- ""
  m <- regexec(paste0("(^|;)", key, "=([^;]+)"), x, perl = TRUE)
  out <- regmatches(x, m)
  vapply(out, function(h) if (length(h) >= 3) h[[3]] else NA_character_, character(1))
}
upsert_info_key <- function(x, key, val) {
  x <- as.character(x); x[is.na(x)] <- ""
  has <- grepl(paste0("(^|;)", key, "="), x, perl = TRUE)
  x[has]  <- sub(paste0("(^|;)", key, "=[^;]*"), paste0("\\1", key, "=", val), x[has], perl = TRUE)
  x[!has] <- ifelse(x[!has] == "", paste0(key, "=", val), paste0(x[!has], ";", key, "=", val))
  x
}

# Count allele freqs 
gt2ac <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "[|]", "/")
  x <- trimws(x)
  
  out <- rep(NA_integer_, length(x))
  out[x %in% c("./.", ".")] <- NA_integer_
  out[x == "0/0"] <- 0L
  out[x == "1/1"] <- 2L
  out[x %in% c("0/1","1/0")] <- 1L
  out
}

gt2ac_from_cell <- function(x){
  gt <- sub("^([^:]+).*", "\\1", as.character(x))
  gt <- stringr::str_replace_all(gt, "[|]", "/")
  gt <- trimws(gt)
  
  out <- rep(NA_integer_, length(gt))
  out[gt %in% c("./.", ".")] <- NA_integer_
  out[gt == "0/0"] <- 0L
  out[gt == "1/1"] <- 2L
  out[gt %in% c("0/1","1/0")] <- 1L
  out
}
# when there are multiple alt alleles need to convert brioche out to standard vcf multi allele
normalize_alt <- function(s) {
  if (is.null(s)) return(s)
  s <- as.character(s)
  
  # keep NAs as NAs
  is_na <- is.na(s)
  
  s <- gsub("[\\[\\]\\(\\)\\{\\}]", "", s, perl = TRUE)
  
  
  s <- gsub("\\s+", "", s, perl = TRUE)
  
  
  parts <- strsplit(s, "[/|,]+", perl = TRUE)
  
  s <- vapply(parts, function(v) {
    v <- toupper(v[nzchar(v)])
    if (length(v) == 0) return(NA_character_)
    v <- v[!duplicated(v)]
    paste(v, collapse = ",")
  }, character(1))
  
  s[is_na] <- NA_character_
  s
}
# Final genotyping function, in scenarios where the new reference is not the same as either the past ref or past alt we need to shift the values of all genotypes aling.
shift_gt_vec <- function(vec) {
  v    <- as.character(vec)
  gt   <- sub("^([^:]+).*", "\\1", v)
  rest <- sub("^[^:]*", "", v)      # keeps leading ':' and the remainder (or "" if none)
  gt   <- gsub("\\|", "/", gt)
  
  parts <- strsplit(gt, "/", fixed = TRUE)
  newgt <- vapply(parts, function(p) {
    p[p == "0"] <- "1"
    p[p == "1"] <- "2"
    paste(p, collapse = "/")
  }, character(1))
  
  paste0(newgt, rest)
}

normalize_strand <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(startsWith(x, "m") | x == "-", "-",
         ifelse(startsWith(x, "p") | x == "+", "+", "+"))
}

# IUPAC reverse-complement (works for single- or multi-base alleles)
revcomp_iupac <- function(x) {
  map <- c(A="T", C="G", G="C", T="A", M="K", R="Y", W="W", S="S", Y="R", K="M",
           V="B", H="D", D="H", B="V", N="N")
  rc1 <- function(s) {
    if (is.na(s) || s == "") return(NA_character_)
    s <- toupper(s)
    ch <- strsplit(s, "", fixed = TRUE)[[1]]
    rc <- rev(ifelse(ch %in% names(map), map[ch], ch))
    paste0(rc, collapse = "")
  }
  vapply(x, rc1, character(1))
}

# Token utilities
csv_tokens <- function(csv) if (is.na(csv) || csv == "") character(0) else strsplit(csv, ",", fixed=TRUE)[[1]]
has_any_token <- function(alt_csv, candidates) {
  if (is.na(alt_csv) || alt_csv == "") return(FALSE)
  altsU <- toupper(strsplit(alt_csv, ",", fixed=TRUE)[[1]])
  any(toupper(candidates) %in% altsU)
}
first_present_token <- function(alt_csv, candidates) {
  if (is.na(alt_csv) || alt_csv == "") return(NA_character_)
  alts  <- strsplit(alt_csv, ",", fixed=TRUE)[[1]]
  altsU <- toupper(alts)
  for (c in candidates) {
    j <- match(toupper(c), altsU)
    if (!is.na(j)) return(alts[j])
  }
  NA_character_
}

# Reorders the ALT CSV with 'want_first' first (if present)
reorder_alt_first <- function(alt_csv, want_first) {
  if (is.na(alt_csv) || is.na(want_first)) return(alt_csv)
  alts  <- strsplit(alt_csv, ",", fixed=TRUE)[[1]]
  altsU <- toupper(alts)
  m1 <- match(toupper(want_first), altsU)
  if (is.na(m1)) return(alt_csv)
  paste(c(alts[m1], alts[setdiff(seq_along(alts), m1)]), collapse=",")
}

# Put two targets first (if present), preserving the rest
reorder_alt_two <- function(alt_csv, firstA, secondA) {
  if (is.na(alt_csv)) return(alt_csv)
  alts  <- strsplit(alt_csv, ",", fixed=TRUE)[[1]]
  altsU <- toupper(alts)
  m1 <- match(toupper(firstA),  altsU)
  m2 <- match(toupper(secondA), altsU)
  ord <- c(m1, m2, setdiff(seq_along(alts), c(m1, m2)))
  paste(alts[!is.na(ord)], collapse=",")
}

# This will probably need to be expanded with more variation during testing. I think I got the bulk of pattern variation here already though
order_contigs <- function(x) {
  x <- unique(as.character(x))
  norm <- tolower(trimws(x))
  
  # unknown bucket
  unk <- grepl("^.*(chrunk|chrun|unk|un|chr0|chrNA|ChrNA)$", norm, perl = TRUE)
  
  # capture number and optional single letter right after the number
  # prefixes allowed: chromosome|chrom|chr|contig|scaffold|segment|seg
  # examples matched: "1", "chr01", "Chromosome12", "contig_7", "scaffold-0003" , "Seg1" Chrom1a, Chromosome2D
  m <- stringr::str_match(
    norm,
    "^.*(?:chromosome|chrom|chr|contig|scaffold|segment|seg)?[ _-]*0*([0-9]+)([a-z])?(?![a-z0-9]).*"
  )
  # columns: 1=full, 2=number, 3=letter (if present)
  digits     <- m[, 2]
  letter_sfx <- m[, 3]
  
  num <- suppressWarnings(as.integer(digits))
  # treat 0 as invalid
  num[!is.na(num) & num == 0L] <- NA_integer_
  
  # sorting key for the letter (empty string sorts before 'a')
  letter_key <- tolower(ifelse(is.na(letter_sfx), "", letter_sfx))
  
  idx_num   <- !is.na(num) & !unk
  idx_other <-  is.na(num) & !unk
  
  # sort numeric first by num, then by letter, then by original string for stability
  ord_num <- order(num[idx_num], letter_key[idx_num], x[idx_num], na.last = TRUE)
  
  c(
    x[idx_num][ord_num],  # numeric + optional letter, ordered
    x[idx_other],         # non-numeric, non-unknown (left as-is)
    x[unk]                # unknowns at the end
  )
}

# Return FIRST ALT (uppercased) from a CSV ALT string (NA if none)
first_alt_uc <- function(x) {
  L <- split_alts(x)
  vapply(L, function(v) if (length(v)) toupper(v[1]) else NA_character_, character(1))
}

# Swap biallelic genotypes (0<->1) while preserving ./., phasing/rest of subfields
swap_biallelic_gt_vec <- function(vec) {
  v  <- as.character(vec); v[is.na(v)] <- ""
  gt <- sub("^([^:]+).*", "\\1", v)
  rest <- sub("^[^:]*", "", v)
  gt <- gsub("\\|", "/", gt)
  gt[gt=="0"] <- "0/0"; gt[gt=="1"] <- "0/1"; gt[gt=="2"] <- "1/1"
  swap_one <- function(g) {
    if (g %in% c("", ".", "./.")) return(if (g=="") "" else "./.")
    if (g == "0/0") return("1/1")
    if (g == "1/1") return("0/0")
    if (g %in% c("0/1","1/0")) return("0/1")
    g
  }
  out <- vapply(gt, swap_one, character(1))
  paste0(out, rest)
}



# For genome info gathering function to search NCBI

fetch_ncbi_meta <- function(accession = "", genomename = "") {
  out <- list(ok = FALSE)
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  tryCatch({
    acc_in  <- trimws(accession)
    name_in <- trimws(genomename)
    
    pick_id <- NULL
    used    <- ""
    
    normalize_summaries <- function(sums) {
      # rentrez returns either a single esummary object or a list of them
      if (is.list(sums) && !is.null(sums$uid)) list(sums)
      else if (is.list(sums)) sums
      else list(sums)
    }
    
    if (nzchar(acc_in)) {
      ids <- character(0)
      
      q1 <- sprintf('%s[Assembly Accession]', acc_in)
      s1 <- rentrez::entrez_search(db = "assembly", term = q1, retmax = 20)
      ids <- unique(c(ids, s1$ids))
      
      if (length(ids) == 0) {
        q2 <- sprintf('"%s"[Assembly Accession]', acc_in)
        s2 <- rentrez::entrez_search(db = "assembly", term = q2, retmax = 20)
        ids <- unique(c(ids, s2$ids))
      }
      if (length(ids) == 0) {
        q3 <- sprintf('"%s"', acc_in)
        s3 <- rentrez::entrez_search(db = "assembly", term = q3, retmax = 20)
        ids <- unique(c(ids, s3$ids))
      }
      
      if (length(ids) > 0) {
        sums <- normalize_summaries(rentrez::entrez_summary(db = "assembly", id = ids))
        matches <- vapply(
          sums,
          function(x) {
            acc_sum <- x[["assemblyaccession"]] %||% x[["accn"]] %||% ""
            toupper(acc_sum) == toupper(acc_in)
          },
          logical(1)
        )
        if (sum(matches) == 1) {
          pick_id <- ids[which(matches)]
          used <- "acc"
        }
      }
    }
    
    
    if (is.null(pick_id) && nzchar(name_in)) {
      qn <- sprintf('"%s"[Assembly Name]', gsub('"', '\\"', name_in))
      sn <- rentrez::entrez_search(db = "assembly", term = qn, retmax = 20)
      if (length(sn$ids) > 0) {
        sums <- normalize_summaries(rentrez::entrez_summary(db = "assembly", id = sn$ids))
        matches <- vapply(
          sums,
          function(x) {
            aname <- x[["assemblyname"]] %||% ""
            tolower(aname) == tolower(name_in)
          },
          logical(1)
        )
        if (sum(matches) == 1) {
          pick_id <- sn$ids[which(matches)]
          used <- "name"
        }
      }
    }
    
    
    if (is.null(pick_id)) {
      message("Accession or genome name were nonspecific unable to populate genomic information")
      return(out)
    }
    
    
    sum <- rentrez::entrez_summary(db = "assembly", id = pick_id)
    
    if (used == "acc") {
      acc_sum <- sum[["assemblyaccession"]] %||% sum[["accn"]] %||% ""
      if (!nzchar(acc_sum) || toupper(acc_sum) != toupper(acc_in)) {
        message("Accession or genome name were nonspecific unable to populate genomic information")
        return(out)
      }
    } else if (used == "name") {
      aname <- sum[["assemblyname"]] %||% ""
      if (!nzchar(aname) || tolower(aname) != tolower(name_in)) {
        message("Accession or genome name were nonspecific unable to populate genomic information")
        return(out)
      }
    }
    
    acc_sum <- sum[["assemblyaccession"]] %||% sum[["accn"]] %||% ""
    aname   <- sum[["assemblyname"]]      %||% ""
    org     <- sum[["organism"]]          %||% sum[["speciesname"]] %||% ""
    taxid   <- sum[["taxid"]]             %||% NA
    ftp     <- sum[["ftppath_refseq"]]    %||% sum[["ftppath_genbank"]] %||% ""
    
    if (!nzchar(acc_sum) || !nzchar(ftp)) {
      message("Accession or genome name were nonspecific unable to populate genomic information")
      return(out)
    }
    
    ftp_https <- sub("^ftp://", "https://", ftp)
    base      <- basename(ftp_https)
    rpt_url   <- paste0(ftp_https, "/", base, "_assembly_report.txt")
    asm_url   <- paste0("https://www.ncbi.nlm.nih.gov/datasets/genome/", acc_sum)
    
    
    
    parse_assembly_report <- function(rpt_url) {
      lines <- readr::read_lines(rpt_url, progress = FALSE)
      hdr_i <- which(grepl("^#\\s*Sequence-Name\\t", lines))[1]
      if (is.na(hdr_i)) stop("assembly_report header line not found")
      header <- sub("^#\\s*", "", lines[hdr_i])
      data   <- lines[(hdr_i + 1L):length(lines)]
      data   <- data[!grepl("^#", data)]  # drop any remaining comment/meta lines
      txt    <- paste(c(header, data), collapse = "\n")
      readr::read_tsv(I(txt), show_col_types = FALSE, progress = FALSE,
                      col_types = readr::cols(.default = "c"))
    }
    
    
    
    ar <- parse_assembly_report(rpt_url)
    
    col_len      <- intersect(c("Sequence-Length","sequence-length","length"), names(ar))
    col_assigned <- intersect(c("Assigned-Molecule","assigned-molecule"), names(ar))
    col_md5      <- intersect(c("MD5","md5"), names(ar))
    col_seqname  <- intersect(c("Sequence-Name","sequence-name"), names(ar))
    
    key_col <- if (length(col_assigned)) col_assigned[1] else col_seqname[1]
    num_key <- suppressWarnings(as.integer(ar[[key_col]]))
    len_val <- suppressWarnings(as.numeric(ar[[col_len[1]]]))
    md5_val <- if (length(col_md5)) ar[[col_md5[1]]] else rep(NA_character_, nrow(ar))
    
    keep <- !is.na(num_key) & num_key >= 1L & !is.na(len_val)
    len_by_num <- tapply(len_val[keep], num_key[keep], function(v) v[1])
    md5_by_num <- tapply(md5_val[keep], num_key[keep], function(v) v[1])
    
    list(
      ok          = TRUE,
      accession   = acc_sum,
      asm_name    = aname,
      organism    = org,
      taxid       = as.integer(taxid),
      asm_url     = asm_url,
      len_by_num  = len_by_num,
      md5_by_num  = md5_by_num
    )
  }, error = function(e) {
    message("Accession or genome name were nonspecific unable to populate genomic information")
    out
  })
}

build_extra_meta <- function(ncbi_meta, genomename) {
  gn <- trimws(ifelse(is.null(genomename), "", genomename))

  esc <- function(x) gsub('"', '\\\\"', x)

  if (isTRUE(ncbi_meta$ok)) {
    # If organism missing but a genomename was provided, use it as a species fallback
    sp <- ncbi_meta$organism
    if (!nzchar(ifelse(is.null(sp), "", sp)) && nzchar(gn)) sp <- gn

    c(
      sprintf('##reference=%s', ncbi_meta$accession),
      sprintf(
        '##assembly=<accession=%s,name="%s"%s%s%s>',
        ncbi_meta$accession,
        esc(ncbi_meta$asm_name),
        if (nzchar(ifelse(is.null(sp), "", sp))) paste0(',species="', esc(sp), '"') else "",
        if (!is.null(ncbi_meta$taxid) && !is.na(ncbi_meta$taxid)) paste0(",taxonomy=", ncbi_meta$taxid) else "",
        if (!is.null(ncbi_meta$asm_url) && nzchar(ncbi_meta$asm_url)) paste0(',url="', esc(ncbi_meta$asm_url), '"') else ""
      )
    )
  } else if (nzchar(gn)) {
    # No NCBI metadata — still record the user-provided genome name
    c(
      sprintf('##reference=%s', gn),
      sprintf('##assembly=<name="%s",species="%s">', esc(gn), esc(gn))
    )
  } else {
    character(0)
  }
}

# Add contig info and build contig lines based on the contigs in the Brioche output
build_contig_lines <- function(ids, meta) {
  ids <- unique(as.character(ids))
  mk_attr <- function(id, len = NA, md5 = NA, meta) {
    parts <- c(paste0("ID=", id))
    if (!is.na(len)) parts <- c(parts, paste0("length=", as.integer(len)))
    if (!is.null(meta$organism) && nzchar(meta$organism)) parts <- c(parts, paste0('species="', meta$organism, '"'))
    if (!is.null(meta$taxid) && !is.na(meta$taxid)) parts <- c(parts, paste0("taxonomy=", meta$taxid))
    paste0("##contig=<", paste(parts, collapse = ","), ">")
  }
  
  # derive numeric key from id: accept chromosome|chrom|chr|contig|scaffold prefixes, case-insensitive
  num_from_id <- function(s) {
    ss <- tolower(trimws(s))
    # unknowns anywhere before, but ending with these tokens
    if (grepl("^.*(chrunk|chrun|unk|un|chr0)$", ss, perl = TRUE)) return(NA_integer_)
    # capture numeric after common prefixes; allow extra leading chars
    hit <- stringr::str_match(
      ss,
      "^.*(?:chromosome|chrom|chr|contig|scaffold|segment|seg)?[ _-]*0*([0-9]+)\\b"
    )[, 2]
    if (is.na(hit)) return(NA_integer_)
    k <- suppressWarnings(as.integer(hit))
    if (is.na(k) || k == 0L) return(NA_integer_) else k
  }
  
  out <- character(length(ids))
  for (i in seq_along(ids)) {
    id <- ids[i]
    if (isTRUE(meta$ok)) {
      k <- num_from_id(id)
      len <- if (!is.na(k) && !is.null(meta$len_by_num[[as.character(k)]]))
        as.numeric(meta$len_by_num[[as.character(k)]]) else NA_real_
      md5 <- if (!is.na(k) && !is.null(meta$md5_by_num[[as.character(k)]]))
        as.character(meta$md5_by_num[[as.character(k)]]) else NA_character_
      out[i] <- mk_attr(id, len = len, md5 = md5, meta = meta)
    } else {
      out[i] <- paste0("##contig=<ID=", id, ">")
    }
  }
  out
}

# If raws genotypes is not a vcf file  
if (!is_vcfraws & !isDArTfile) {
  
  genotypesRAWS  <- read.table(raw_path, header = TRUE, sep = ifelse(grepl("\\.tsv$", raw_path, ignore.case=TRUE), "\t", ""),
                               check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
  
}


Finalmappings  <- read.csv(map_path, check.names = FALSE, stringsAsFactors = FALSE)

#VCFchroms <- unique(Finalmappings$saccver)
VCFchroms <- order_contigs(Finalmappings$saccver)




ncbi_meta <- list(ok = FALSE)
if (nzchar(Accession) || nzchar(genomename)) {
  ncbi_meta <- fetch_ncbi_meta(Accession, genomename)
  if (isTRUE(ncbi_meta$ok)) {
    message(sprintf("NCBI assembly  : %s (%s) taxid=%s",
                    ncbi_meta$accession, ncbi_meta$asm_name, ncbi_meta$taxid))
    message(sprintf("NCBI URL       : %s", ncbi_meta$asm_url))
  } else {
    message("NCBI metadata   : none found (proceeding without enrichment)")
  }
}
if (nzchar(genomename)) {
  org_present <- !is.null(ncbi_meta$organism) && nzchar(ncbi_meta$organism)
  if (!org_present) ncbi_meta$organism <- genomename
}

# SNP-chip pathway (raw genotypes table to vcf)
if (!is_vcfraws & is_snpchip) {

  ## ---- 0) Delimiter + header (column names) ----
  delim <- if (grepl("\\.(tsv|txt)$", raw_path, ignore.case = TRUE)) "\t" else ","
  hdr0  <- utils::read.table(
    file = raw_path, header = TRUE, sep = delim, nrows = 0,
    check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE
  )
  cn_raw   <- names(hdr0)
  cn_fixed <- make.unique(cn_raw, sep = "...")  # ensure unique for VCF/sample IDs

  ## Required columns in SNP-chip table
  idx_Name <- match("Name", cn_raw)
  idx_REF  <- match("REF",  cn_raw)
  idx_ALT  <- match("ALT",  cn_raw)
  if (any(is.na(c(idx_Name, idx_REF, idx_ALT))))
    stop("SNP chip file must have columns: Name, REF, ALT (in that order before samples).")

  sample_idx <- setdiff(seq_along(cn_raw), c(idx_Name, idx_REF, idx_ALT))
  if (length(sample_idx) == 0L) stop("No sample genotype columns found after ALT.")
  sample_ids <- cn_fixed[sample_idx]  # VCF sample IDs

  ## Mapping frame (single hit per Name and Brioche REF/ALT targets)
  fm_base <- Finalmappings %>%
    dplyr::transmute(
      Name      = qaccver,
      CHROM_map = saccver,
      POS_map   = SNPpos,
      REF_map   = toupper(Ref),
      ALT_map   = normalize_alt(ALT)
    )

  fm_hits <- Finalmappings %>% dplyr::count(qaccver, name = "MAP_HITS")
  fm <- fm_base %>%
    dplyr::left_join(fm_hits, by = c("Name" = "qaccver")) %>%
    dplyr::mutate(MAP_HITS = ifelse(is.na(MAP_HITS), 0L, MAP_HITS))
  fm_Name <- fm$Name

  ## ---- 1) Write VCF header ----
  fileDate     <- format(Sys.Date(), "%Y%m%d")
  VCFchroms    <- unique(c(order_contigs(Finalmappings$saccver), "chrUnk"))
  contig_lines <- unique(build_contig_lines(VCFchroms, ncbi_meta))
  extra_meta   <- build_extra_meta(ncbi_meta, genomename)

  header_lines <- c(
    "##fileformat=VCFv4.3",
    sprintf("##fileDate=%s", fileDate),
    "##source=Brioche-VCF-build",
    extra_meta,
    "##FILTER=<ID=LowQual,Description=\"Low quality or ambiguous mapping\">",
    '##INFO=<ID=MAPSTATUS,Number=1,Type=String,Description="Mapping status from anchoring_script.R (Unique_mapping|Failed_to_map_uniquely)">',
    '##INFO=<ID=ORIENTATIONSTATUS,Number=1,Type=String,Description="Whether the marker has been oriented w.r.t. reference alleles (Yes|No)">',
    "##INFO=<ID=MAF,Number=A,Type=Float,Description=\"Minor Allele frequency per ALT\">",
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count per ALT\">",
    "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=NU,Number=1,Type=Integer,Description=\"Null allele flag from SNP chip (1=NC; 0=otherwise[N or genotype])\">",
    contig_lines
  )
  header_lines <- header_lines[!duplicated(header_lines)]
  writeLines(header_lines, con = Outputfilename)

  vcf_header <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", sample_ids)
  data.table::fwrite(as.list(vcf_header), file = Outputfilename, sep = "\t",
                     append = TRUE, col.names = FALSE, quote = FALSE)

  ## ---- 2) Streaming chunked processing (base R) ----
  chunk_size <- 1500L
  is_blank <- function(x) is.na(x) | x == "" | x == "."
  con_chunk <- if (grepl("\\.gz$", raw_path, ignore.case = TRUE)) gzfile(raw_path, "rt") else file(raw_path, "rt")
  on.exit(try(close(con_chunk), silent = TRUE), add = TRUE)
  invisible(readLines(con_chunk, n = 1L, warn = FALSE))  # consume header line

  repeat {
    raw_lines <- readLines(con_chunk, n = chunk_size, warn = FALSE)
    if (!length(raw_lines)) break

    df <- utils::read.table(text = raw_lines, sep = delim, header = FALSE,
                            quote = "", comment.char = "", stringsAsFactors = FALSE,
                            check.names = FALSE, col.names = cn_fixed, fill = TRUE)

    Name    <- df[[idx_Name]]
    REF_old <- toupper(df[[idx_REF]])
    ALT_old <- normalize_alt(df[[idx_ALT]])

    # Enforce biallelic on INPUT side
    ALT_old <- sub(",.*$", "", ALT_old)

    m     <- match(Name, fm_Name)
    hasM  <- !is.na(m)
    hits  <- ifelse(hasM, fm$MAP_HITS[m], 0L)
    u_map <- hasM & hits == 1L

    ## Coordinates (only when uniquely mapped)
    CHROM <- ifelse(u_map, fm$CHROM_map[m], "chrUnk")
    POS   <- ifelse(u_map, fm$POS_map[m],   0L)

    ## Candidate Brioche targets (only used if uniquely mapped)
    BriocheREF <- ifelse(u_map, toupper(fm$REF_map[m]), NA_character_)
    BriocheALT <- ifelse(u_map, toupper(sub(",.*$", "", fm$ALT_map[m])), NA_character_)  # single ALT
    # RC of Brioche alleles
    BriocheREFRC <- ifelse(!is.na(BriocheREF), revcomp_iupac(BriocheREF), NA_character_)
    BriocheALTRC <- ifelse(!is.na(BriocheALT), revcomp_iupac(BriocheALT), NA_character_)

    # Biallelic viability: need exactly one ALT on both sides
    old_alt_len <- ifelse(is.na(ALT_old) | ALT_old == "" | ALT_old == ".", 0L, 1L)
    new_alt_len <- ifelse(is.na(BriocheALT) | BriocheALT == "" | BriocheALT == ".", 0L, 1L)
    bial_ok     <- u_map & (old_alt_len == 1L) & (new_alt_len == 1L)

    # Identity / RC-identity / Swap / RC-swap tests (case-insensitive)
    genotypeREF <- REF_old
    genotypeALT <- ALT_old

    ident      <- bial_ok & !is.na(BriocheREF) & (BriocheREF  == genotypeREF)
    ident_rc   <- bial_ok & !is.na(BriocheREFRC) & (BriocheREFRC == genotypeREF)
    swap       <- bial_ok & !is.na(BriocheREF) & (BriocheREF  == genotypeALT)
    swap_rc    <- bial_ok & !is.na(BriocheREFRC) & (BriocheREFRC == genotypeALT)

    oriented_now <- ident | ident_rc | swap | swap_rc

    # Output REF/ALT: when oriented, force to Brioche REF/ALT; otherwise keep input
    REF_out <- ifelse(oriented_now, BriocheREF, REF_old)
    ALT_out <- ifelse(oriented_now, BriocheALT, ALT_old)

    # Always keep ALT single (strict biallelic)
    ALT_out <- sub(",.*$", "", ALT_out)
    ALT_out[is.na(ALT_out) | ALT_out == ""] <- "."

    # QUAL / FILTER / INFO (MAP only here; ORIENTATIONSTATUS appended later)
    pos_int     <- suppressWarnings(as.integer(POS))
    coord_bad   <- is.na(pos_int) | pos_int <= 0L | is_blank(CHROM)
    unknown_now <- !u_map | coord_bad

    QUAL   <- ifelse(unknown_now, 0L, 100L)
    FILTER <- ifelse(unknown_now, "LowQual", "PASS")
    INFO   <- ifelse(unknown_now, "MAPSTATUS=Failed_to_map_uniquely", "MAPSTATUS=Unique_mapping")

    ## ## Genotypes + NU flag (1=NC, 0 otherwise; “N” scores 0 as does true genotype data)
    geno_raw <- as.matrix(df[, sample_idx, drop = FALSE])
    geno_raw_upper <- toupper(trimws(geno_raw))
    NU <- (geno_raw_upper == "NC")
    storage.mode(NU) <- "integer"
    NU[is.na(NU)] <- 0L

    geno <- as.matrix(df[, sample_idx, drop = FALSE])
    geno[] <- toupper(gsub("\\|", "/", geno, perl = TRUE))
    miss <- is.na(geno) | geno == "" | geno == "." | geno == "-" | geno == "N" | geno == "NC"
    if (any(miss)) geno[miss] <- "./."
    dbl <- grepl("^[012]{2}$", geno, perl = TRUE)
    if (any(dbl)) geno[dbl] <- sub("^([012])([012])$", "\\1/\\2", geno[dbl], perl = TRUE)

    # SAME/RC-SAME: no genotype change; SWAP/RC-SWAP: flip 0<->2, keep 1 hetero (0/1 stays 0/1 as we don't actually know phasing so keep it consistent and use / for uncertainty vs |)
    rsame <- ifelse(is.na(ident | ident_rc), FALSE, (ident | ident_rc))
    rswap <- ifelse(is.na(swap | swap_rc),  FALSE, (swap  | swap_rc))

    if (any(rsame)) {
      g <- geno
      d0 <- (g == "0") & matrix(rep(rsame, times = ncol(geno)), nrow = nrow(geno))
      d1 <- (g == "1") & matrix(rep(rsame, times = ncol(geno)), nrow = nrow(geno))
      d2 <- (g == "2") & matrix(rep(rsame, times = ncol(geno)), nrow = nrow(geno))
      if (any(d0)) g[d0] <- "0/0"
      if (any(d1)) g[d1] <- "0/1"
      if (any(d2)) g[d2] <- "1/1"
      geno[matrix(rep(rsame, times = ncol(geno)), nrow = nrow(geno))] <-
        g[matrix(rep(rsame, times = ncol(geno)), nrow = nrow(geno))]
    }
    if (any(rswap)) {
      g <- geno
      d0 <- (g == "0") & matrix(rep(rswap, times = ncol(geno)), nrow = nrow(geno))
      d1 <- (g == "1") & matrix(rep(rswap, times = ncol(geno)), nrow = nrow(geno))
      d2 <- (g == "2") & matrix(rep(rswap, times = ncol(geno)), nrow = nrow(geno))
      if (any(d0)) g[d0] <- "1/1"
      if (any(d1)) g[d1] <- "0/1"
      if (any(d2)) g[d2] <- "0/0"
      geno[matrix(rep(rswap, times = ncol(geno)), nrow = nrow(geno))] <-
        g[matrix(rep(rswap, times = ncol(geno)), nrow = nrow(geno))]
    }

    # Normalize any lingering single digits
    i0 <- geno == "0"; if (any(i0)) geno[i0] <- "0/0"
    i1 <- geno == "1"; if (any(i1)) geno[i1] <- "0/1"
    i2 <- geno == "2"; if (any(i2)) geno[i2] <- "1/1"

    ## Append ORIENTATIONSTATUS
    INFO <- ifelse(oriented_now,
                   paste0(INFO, ";ORIENTATIONSTATUS=Yes"),
                   paste0(INFO, ";ORIENTATIONSTATUS=No"))

    ## AC/AN/MAF (strictly biallelic)
    nr <- nrow(geno)
    ac1 <- integer(nr); an <- integer(nr)
    for (jj in seq_len(ncol(geno))) {
      x <- geno[, jj]
      m <- (x == "./.")
      a1 <- substr(x, 1L, 1L); a1[m] <- NA_character_
      a2 <- substr(x, 3L, 3L); a2[m] <- NA_character_
      ac1 <- ac1 + as.integer((!is.na(a1)) & (a1 == "1")) + as.integer((!is.na(a2)) & (a2 == "1"))
      an  <- an  + as.integer(!m) * 2L
    }
    # ensure ALT_out is single; AC is for allele "1"
    AC_str  <- as.character(ac1)
    AF1     <- ifelse(an > 0L, ac1 / an, NA_real_)
    MAF_str <- sprintf("%.6g", AF1)
    INFO <- ifelse(
      is.na(INFO) | INFO == ".",
      paste0("MAF=", MAF_str, ";AC=", AC_str, ";AN=", an),
      paste0(INFO, ";MAF=", MAF_str, ";AC=", AC_str, ";AN=", an)
    )

    ## Join GT:NU
    NU_chr <- matrix(as.character(NU), nrow = nrow(NU), ncol = ncol(NU))
    NU_chr[is.na(NU_chr)] <- "0"
    geno_df <- as.data.frame(geno, stringsAsFactors = FALSE, check.names = FALSE)
    for (jj in seq_len(ncol(geno_df))) {
      geno_df[[jj]] <- paste0(geno_df[[jj]], ":", NU_chr[, jj])
    }

    body <- data.table::data.table(
      `#CHROM` = CHROM, POS = POS, ID = Name,
      REF = REF_out, ALT = ALT_out, QUAL = QUAL, FILTER = FILTER, INFO = INFO, FORMAT = "GT:NU"
    )
    body <- cbind(body, data.table::as.data.table(geno_df))

    data.table::fwrite(body, file = Outputfilename, sep = "\t",
                       append = TRUE, col.names = FALSE, quote = FALSE)

    rm(df, geno, geno_df, NU, NU_chr, body, ac1, an); gc(FALSE)
  }

  message(sprintf("Wrote: %s", Outputfilename))
}


 
# vcf pathway 
# vcf pathway 
if (is_vcfraws) {

  ## ---- Local helpers (define only if missing) ----
  if (!exists("nzchar_safe", mode = "function")) {
    nzchar_safe <- function(x) !is.na(x) & x != "" & x != "."
  }

  if (!exists("normalize_alt", mode = "function")) {
    # Normalize ALT csv (strip brackets/whitespace, uppercase, unique tokens, keep order)
    normalize_alt <- function(s) {
      if (is.null(s)) return(s)
      s0 <- as.character(s)
      is_na <- is.na(s0)
      s0 <- gsub("[\\[\\]\\(\\)\\{\\}]", "", s0, perl = TRUE)
      s0 <- gsub("\\s+", "", s0, perl = TRUE)
      parts <- strsplit(s0, "[/|,]+", perl = TRUE)
      out <- vapply(parts, function(v) {
        v <- toupper(v[nzchar(v)])
        if (!length(v)) return(NA_character_)
        v <- v[!duplicated(v)]
        paste(v, collapse = ",")
      }, character(1))
      out[is_na] <- NA_character_
      out
    }
  }

  if (!exists("first_alt_uc", mode = "function")) {
    first_alt_uc <- function(alt_csv) {
      x <- as.character(alt_csv)
      vapply(x, function(s) {
        if (is.na(s) || s == "" || s == ".") return(NA_character_)
        toupper(strsplit(s, ",", fixed = TRUE)[[1]][1])
      }, character(1))
    }
  }

  if (!exists("swap_biallelic_gt_vec", mode = "function")) {
    # Swap 0<->1 for diploid biallelic genotypes; preserve phasing and trailing subfields
    swap_biallelic_gt_vec <- function(vec) {
      v <- as.character(vec); v[is.na(v)] <- ""
      gt   <- sub("^([^:]+).*", "\\1", v)
      rest <- sub("^[^:]*", "", v)
      is_single <- gt %in% c("0","1")
      gt[is_single] <- paste0(gt[is_single], "/", gt[is_single])
      rec1 <- function(g) {
        if (g == "" || g == "." || g == "./.") return(if (g == "") "" else "./.")
        sep <- if (grepl("\\|", g)) "|" else "/"
        a <- strsplit(g, "[/|]", perl = TRUE)[[1]]
        a[a == "0"] <- "1"
        a[a == "1"] <- "0"
        paste(a, collapse = sep)
      }
      out_gt <- vapply(gt, rec1, character(1))
      paste0(out_gt, rest)
    }
  }

  if (!exists("info_get", mode = "function")) {
    info_get <- function(info_vec, key){
      x <- as.character(info_vec); x[is.na(x)] <- ""
      m <- regexec(paste0("(^|;)", key, "=([^;]+)"), x, perl = TRUE)
      regmatches(x, m) |>
        lapply(function(hit) if (length(hit) >= 3) hit[3] else NA_character_) |>
        unlist(use.names = FALSE)
    }
  }

  if (!exists("info_strip_keys", mode = "function")) {
    info_strip_keys <- function(info_vec, keys) {
      x <- as.character(info_vec); x[is.na(x)] <- ""
      for (k in keys) x <- gsub(paste0("(^|;)", k, "=[^;]*"), "\\1", x, perl = TRUE)
      x <- gsub(";{2,}", ";", x, perl = TRUE)
      gsub("^;|;$", "", x, perl = TRUE)
    }
  }

  if (!exists("upsert_info_key_vec", mode = "function")) {
    upsert_info_key_vec <- function(x, key, val) {
      y <- info_strip_keys(x, key)
      ifelse(y == "" | is.na(y), paste0(key, "=", val), paste0(y, ";", key, "=", val))
    }
  }

  if (!exists("revcomp_iupac", mode = "function")) {
    revcomp_iupac <- function(x) {
      map <- c(A="T", C="G", G="C", T="A", M="K", R="Y", W="W", S="S", Y="R", K="M",
               V="B", H="D", D="H", B="V", N="N")
      rc1 <- function(s) {
        if (is.na(s) || s == "") return(NA_character_)
        s <- toupper(s)
        ch <- strsplit(s, "", fixed = TRUE)[[1]]
        rc <- rev(ifelse(ch %in% names(map), map[ch], ch))
        paste0(rc, collapse = "")
      }
      vapply(x, rc1, character(1))
    }
  }

  if (!exists("recode_gt_vec", mode = "function")) {
    # Recode GT with 0-based allele index map; preserves phasing & subfields
    recode_gt_vec <- function(vec, idx_map) {
      v <- as.character(vec); v[is.na(v)] <- ""
      if (!length(idx_map)) return(v)
      gt   <- sub("^([^:]+).*", "\\1", v)
      rest <- sub("^[^:]*", "", v)
      is_single <- gt %in% c("0","1","2","3","4","5","6","7","8","9")
      gt[is_single] <- paste0(gt[is_single], "/", gt[is_single])
      out_gt <- vapply(gt, function(g){
        if (g == "" || g == "." || g == "./.") return(if (g=="") "" else "./.")
        sep <- if (grepl("\\|", g)) "|" else "/"
        a <- strsplit(g, "[/|]", perl = TRUE)[[1]]
        map_allele <- function(x) {
          if (x == ".") return(".")
          xi <- suppressWarnings(as.integer(x))
          if (is.na(xi) || xi + 1L > length(idx_map)) return(".")
          y <- idx_map[xi + 1L]
          if (is.na(y)) "." else as.character(y)
        }
        a2 <- vapply(a, map_allele, character(1))
        if (any(a2 == ".")) return("./.")
        paste(a2, collapse = sep)
      }, character(1))
      paste0(out_gt, rest)
    }
  }

  if (!exists("build_symbol_maps", mode = "function")) {
    build_symbol_maps <- function(ref_old, alts_old_chr, ref_new, alts_new_chr) {
      alts_old_chr <- if (length(alts_old_chr)) alts_old_chr else character(0)
      alts_new_chr <- if (length(alts_new_chr)) alts_new_chr else character(0)
      sym_old <- c(ref_old, alts_old_chr)
      sym_new <- c(ref_new, alts_new_chr)
      idx_map <- integer(length(sym_old))
      for (i in seq_along(sym_old)) {
        m <- match(toupper(sym_old[i]), toupper(sym_new))
        idx_map[i] <- ifelse(is.na(m), NA_integer_, m - 1L)
      }
      idx_map
    }
  }

  split_alts_csv <- function(s) {
    if (is.na(s) || s == "" || s == ".") character(0) else strsplit(s, ",", fixed = TRUE)[[1]]
  }

  if (!exists("order_contigs", mode = "function")) {
    order_contigs <- function(x) unique(as.character(x))
  }
  if (!exists("build_contig_lines", mode = "function")) {
    build_contig_lines <- function(ids, meta) {
      ids <- unique(as.character(ids))
      paste0("##contig=<ID=", ids, ">")
    }
  }
  if (!exists("build_extra_meta", mode = "function")) {
    build_extra_meta <- function(ncbi_meta, genomename) character(0)
  }

  ## -------- mapping frame (REF/ALT targets) ----------
  req <- c("qaccver","saccver","SNPpos","Ref","ALT")
  missing <- setdiff(req, names(Finalmappings))
  if (length(missing))
    stop(sprintf("Finalmappings is missing required columns: %s", paste(missing, collapse=", ")))

  suppressPackageStartupMessages(library(dplyr))
  fm_one <- Finalmappings %>%
    transmute(qaccver, saccver, SNPpos,
              REF_map = toupper(Ref),
              ALT_map = normalize_alt(ALT)) %>%
    group_by(qaccver) %>% slice(1L) %>% ungroup()
  fm_hits <- Finalmappings %>% count(qaccver, name = "MAP_HITS")

  fm <- fm_one %>%
    left_join(fm_hits, by = "qaccver") %>%
    transmute(
      ID        = qaccver,
      CHROM_map = saccver,
      POS_map   = SNPpos,
      REF_map   = REF_map,
      ALT_map   = ALT_map,
      MAP_HITS  = ifelse(is.na(MAP_HITS), 0L, MAP_HITS)
    )

  ## -------- read input header ----------
  read_vcf_header <- function(path) {
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    meta_lines <- character(0); header_line <- NULL; header_n <- 0L
    repeat {
      ln <- readLines(con, n = 1L, warn = FALSE); if (!length(ln)) break
      header_n <- header_n + 1L; ln <- sub("^\ufeff", "", ln, perl = TRUE)
      if (startsWith(ln, "##")) meta_lines <- c(meta_lines, ln)
      else if (startsWith(ln, "#")) { header_line <- ln; break }
      else stop("Malformed VCF: missing #CHROM header line.")
    }
    if (is.null(header_line)) stop("Malformed VCF: missing #CHROM header line (EOF).")
    list(meta = meta_lines, header = header_line, n = header_n)
  }

  raw_in <- raw_path
  H <- read_vcf_header(raw_in)
  meta_lines  <- H$meta
  header_line <- H$header
  header_n    <- H$n

  colnames_vcf <- strsplit(sub("^#", "", header_line), "\t", fixed = TRUE)[[1]]
  ix_CHROM <- match("CHROM",  colnames_vcf)
  ix_POS   <- match("POS",    colnames_vcf)
  ix_ID    <- match("ID",     colnames_vcf)
  ix_REF   <- match("REF",    colnames_vcf)
  ix_ALT   <- match("ALT",    colnames_vcf)
  ix_QUAL  <- match("QUAL",   colnames_vcf)
  ix_FIL   <- match("FILTER", colnames_vcf)
  ix_INFO  <- match("INFO",   colnames_vcf)
  ix_FMT   <- match("FORMAT", colnames_vcf)
  stopifnot(!any(is.na(c(ix_CHROM, ix_POS, ix_ID, ix_REF, ix_ALT))))
  sample_cols_idx <- if (is.na(ix_FMT)) integer(0) else (ix_FMT + 1L):length(colnames_vcf)

  ## -------- write new header (preserve sample header line) ----------
  meta_out   <- meta_lines[!grepl("^##(contig=|reference=|assembly=<)", meta_lines)]
  meta_out   <- meta_out[!duplicated(meta_out)]
  extra_meta <- build_extra_meta(ncbi_meta, genomename)

  contigs_vec  <- unique(c(order_contigs(fm$CHROM_map), "chrUnk"))
  contig_lines <- unique(build_contig_lines(contigs_vec, ncbi_meta))

  if (!any(grepl('^##INFO=<ID=MAPSTATUS,', meta_out))) {
    meta_out <- c(meta_out,
      '##INFO=<ID=MAPSTATUS,Number=1,Type=String,Description="Mapping status from anchoring_script.R (Unique_mapping|Failed_to_map_uniquely)">')
  }
  if (!any(grepl('^##INFO=<ID=ORIENTATIONSTATUS,', meta_out))) {
    meta_out <- c(meta_out,
      '##INFO=<ID=ORIENTATIONSTATUS,Number=1,Type=String,Description="Sticky orientation flag (Yes|No)">')
  }

  writeLines(c(meta_out, extra_meta, contig_lines, header_line), con = Outputfilename)

  ## -------- stream body (chunk = 5000) ----------
  ct_all_as_char <- rep("character", length(colnames_vcf))
  names(ct_all_as_char) <- colnames_vcf

  chunk_size <- 5000L
  is_gz      <- grepl("\\.gz$", raw_in, ignore.case = TRUE)
  zcat_cmd   <- if (nzchar(Sys.which("zcat"))) "zcat" else "gzip -dc"

  get_total_lines <- function(path, is_gz, zcmd){
    cmd <- if (is_gz) sprintf("%s %s | wc -l", zcmd, shQuote(path))
           else        sprintf("wc -l < %s", shQuote(path))
    out <- tryCatch(system(cmd, intern = TRUE), error = function(e) NA_character_)
    as.integer(gsub("[^0-9]", "", out[1]))
  }
  total_lines <- get_total_lines(raw_in, is_gz, zcat_cmd)
  if (is.na(total_lines) || total_lines < header_n)
    stop(sprintf("Unable to count lines or file shorter than header: total=%s header_n=%s",
                 as.character(total_lines), as.character(header_n)))

  first_data_line <- header_n + 1L
  data_lines      <- max(0L, total_lines - header_n)
  consumed        <- 0L

  repeat {
    if (consumed >= data_lines) break
    n_to_read <- min(chunk_size, data_lines - consumed)

    if (!is_gz) {
      df <- data.table::fread(
        file       = raw_in, sep = "\t", header = FALSE,
        skip       = (first_data_line - 1L) + consumed, nrows = n_to_read,
        col.names  = colnames_vcf, colClasses = ct_all_as_char,
        data.table = FALSE, integer64 = "character", showProgress = FALSE
      )
    } else {
      cmd <- sprintf("%s %s | tail -n +%d | head -n %d",
                     zcat_cmd, shQuote(raw_in), first_data_line + consumed, n_to_read)
      df <- data.table::fread(
        cmd        = cmd, sep = "\t", header = FALSE,
        col.names  = colnames_vcf, colClasses = ct_all_as_char,
        data.table = FALSE, integer64 = "character", showProgress = FALSE
      )
    }
    if (!nrow(df)) break
    df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)

    IDv      <- df[[ix_ID]]
    REF_in   <- toupper(df[[ix_REF]])
    ALT_in   <- as.character(df[[ix_ALT]])   # keep original case for output
    ALT_inU  <- toupper(ALT_in)

    # Current-map lookup
    j          <- match(IDv, fm$ID)
    in_map     <- !is.na(j)
    hits       <- ifelse(in_map, fm$MAP_HITS[j], 0L)
    chr_ok     <- in_map & nzchar_safe(fm$CHROM_map[j])
    pos_ok     <- in_map & !is.na(fm$POS_map[j]) & suppressWarnings(as.numeric(fm$POS_map[j]) > 0)
    mapped_now <- in_map & (hits == 1L) & chr_ok & pos_ok

    # Overwrite per-iteration coords + QUAL/FILTER (NOT STICKY)
    df[[ix_CHROM]] <- ifelse(mapped_now, fm$CHROM_map[j], "chrUnk")
    df[[ix_POS]]   <- ifelse(mapped_now, as.character(fm$POS_map[j]), "0")
    if (!is.na(ix_QUAL)) df[[ix_QUAL]] <- ifelse(mapped_now, "100", "0")
    if (!is.na(ix_FIL )) df[[ix_FIL ]] <- ifelse(mapped_now, "PASS", "LowQual")

    # Sticky ORIENTATIONSTATUS read (assume "No" if missing)
    orient_in_raw <- if (!is.na(ix_INFO)) df[[ix_INFO]] else rep("", nrow(df))
    orient_in     <- info_get(orient_in_raw, "ORIENTATIONSTATUS")
    orient_in     <- ifelse(is.na(orient_in) | orient_in == "" | orient_in == ".", "No", orient_in)
    orient_in_yes <- toupper(orient_in) == "YES"

    # Prepare Brioche alleles (targets) for all rows
    BriocheREF   <- ifelse(mapped_now, toupper(fm$REF_map[j]), NA_character_)
    BriocheALT   <- ifelse(mapped_now, first_alt_uc(fm$ALT_map[j]), NA_character_)
    BriocheREFRC <- revcomp_iupac(BriocheREF)
    BriocheALTRC <- revcomp_iupac(BriocheALT)

    # ---------- First-iteration orientation (only when ORIENTATIONSTATUS = No) ----------
    gen_REF <- REF_in
    gen_ALT <- first_alt_uc(ALT_inU)

    old_alt_single <- !grepl(",", ALT_inU, fixed = TRUE) & nzchar_safe(ALT_inU)
    bri_alt_single <- mapped_now & nzchar_safe(BriocheALT) &
                      !grepl(",", ifelse(mapped_now, fm$ALT_map[j], ""), fixed = TRUE)

    eligible <- mapped_now & !orient_in_yes & old_alt_single & bri_alt_single & nzchar_safe(BriocheREF)

    cond1 <- eligible & (BriocheREF   == gen_REF)  # keep REF/ALT, no GT change
    cond2 <- eligible & (BriocheREFRC == gen_REF)  # set REF/ALT to Brioche, no GT change (flip to Brioche alleles)
    cond3 <- eligible & (BriocheREF   == gen_ALT)  # set REF/ALT to Brioche, swap GT
    cond4 <- eligible & (BriocheREFRC == gen_ALT)  # set REF/ALT to Brioche, swap GT

    swap_rows    <- which(cond3 | cond4)
    oriented_now <- cond1 | cond2 | cond3 | cond4

    # Apply REF/ALT updates for first-iteration orientation
    REF_out <- df[[ix_REF]]
    ALT_out <- df[[ix_ALT]]
    idx_bri <- which(cond2 | cond3 | cond4)
    if (length(idx_bri)) {
      REF_out[idx_bri] <- BriocheREF[idx_bri]
      ALT_out[idx_bri] <- BriocheALT[idx_bri]
    }
    df[[ix_REF]] <- REF_out
    df[[ix_ALT]] <- ALT_out

    # Genotype swap for rows needing it (biallelic diploid only)
    if (length(swap_rows) && length(sample_cols_idx)) {
      for (jj in sample_cols_idx) {
        colv <- df[[jj]]
        colv[swap_rows] <- swap_biallelic_gt_vec(colv[swap_rows])
        df[[jj]] <- colv
      }
    }

    # ---------- Subsequent-iteration logic (ORIENTATIONSTATUS = Yes) : allow tri+ ----------
    tri_rows <- which(mapped_now & orient_in_yes & nzchar(BriocheREF))
    if (length(tri_rows)) {
      REF_cur <- df[[ix_REF]]
      ALT_cur <- df[[ix_ALT]]

      idx_map_by_row <- list()
      recode_needed  <- logical(length(tri_rows))

      for (k in seq_along(tri_rows)) {
        i <- tri_rows[k]

        roC <- REF_cur[i]
        roU <- toupper(roC)

        altC_vec <- split_alts_csv(as.character(ALT_cur[i]))
        altU_vec <- toupper(altC_vec)

        rn <- BriocheREF[i]
        if (!nzchar(rn)) next

        if (roU == rn) {
          # Identity: keep as-is
          next
        }

        m_in_old <- match(rn, altU_vec)

        if (!is.na(m_in_old)) {
          # swap-to-ref with multi-ALT preserved
          new_alt_vecC <- c(roC, altC_vec[-m_in_old])   # old REF becomes first ALT
          REF_cur[i]   <- rn
          ALT_cur[i]   <- if (length(new_alt_vecC)) paste(new_alt_vecC, collapse = ",") else "."
          idx_map_by_row[[as.character(i)]] <- build_symbol_maps(roU, altU_vec, rn, toupper(new_alt_vecC))
          recode_needed[k] <- TRUE
        } else {
          # tri+ expansion: add rn as new REF; push old REF and all old ALTs to ALT
          new_alt_vecC <- c(roC, altC_vec)
          REF_cur[i]   <- rn
          ALT_cur[i]   <- if (length(new_alt_vecC)) paste(new_alt_vecC, collapse = ",") else "."
          idx_map_by_row[[as.character(i)]] <- build_symbol_maps(roU, altU_vec, rn, toupper(new_alt_vecC))
          recode_needed[k] <- TRUE
        }
      }

      # Install updated REF/ALT
      df[[ix_REF]] <- REF_cur
      df[[ix_ALT]] <- ALT_cur

      # Recode genotypes for rows that changed allele indexing
      if (length(sample_cols_idx) && any(recode_needed)) {
        rows_to_recode <- tri_rows[recode_needed]
        for (jj in sample_cols_idx) {
          colv <- df[[jj]]
          for (rr in rows_to_recode) {
            idx_map <- idx_map_by_row[[as.character(rr)]]
            if (!is.null(idx_map)) colv[rr] <- recode_gt_vec(colv[rr], idx_map)
          }
          df[[jj]] <- colv
        }
      }
    }

    # ---------- Maintain sticky ORIENTATIONSTATUS + per-iteration MAPSTATUS ----------
    status_val <- ifelse(mapped_now, "Unique_mapping", "Failed_to_map_uniquely")
    orient_out <- ifelse(orient_in_yes | oriented_now, "Yes", "No")

    if (!is.na(ix_INFO)) {
      df[[ix_INFO]] <- info_strip_keys(df[[ix_INFO]], c("MAPSTATUS","ORIENTATIONSTATUS"))
      df[[ix_INFO]] <- upsert_info_key_vec(df[[ix_INFO]], "MAPSTATUS",         status_val)
      df[[ix_INFO]] <- upsert_info_key_vec(df[[ix_INFO]], "ORIENTATIONSTATUS", orient_out)
    }

    data.table::fwrite(df, file = Outputfilename, sep = "\t",
                       append = TRUE, col.names = FALSE, quote = FALSE)

    consumed <- consumed + nrow(df)
    rm(df); gc(FALSE)
  }

  message(sprintf("Wrote: %s", Outputfilename))
}



# dart file 
if (isDArTfile) {

  # ---- Helpers required by the DArT pathway ----

  detect_dart_layout_base <- function(path) {
    # Find the header line containing AlleleID & SNP; detect delimiter; find where genotypes start
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)

    topskip <- 0L
    header  <- NULL

    # scan for header
    for (i in seq_len(10000L)) {
      ln <- readLines(con, n = 1L, warn = FALSE)
      if (!length(ln)) break
      ln <- sub("^\ufeff", "", ln, perl = TRUE)  # strip BOM if present
      if (!nzchar(ln)) { topskip <- topskip + 1L; next }
      if (!(grepl("\t", ln, fixed = TRUE) || grepl(",", ln, fixed = TRUE))) { topskip <- topskip + 1L; next }

      delim <- if (grepl("\t", ln, fixed = TRUE)) "\t" else ","
      toks  <- strsplit(ln, delim, fixed = TRUE)[[1]]
      toks2 <- trimws(gsub('^"|"$', '', toks))
      if (any(tolower(toks2) == "alleleid") && any(tolower(toks2) == "snp")) {
        header <- ln
        break
      }
      topskip <- topskip + 1L
    }
    if (is.null(header)) stop("Could not locate DArT header line containing 'AlleleID' and 'SNP'.")

    delim <- if (grepl("\t", header, fixed = TRUE)) "\t" else ","
    cn    <- strsplit(header, delim, fixed = TRUE)[[1]]
    cn    <- trimws(gsub('^"|"$', '', cn))

    # read the first data line after header to detect where genotypes start
    first_data <- ""
    while (nzchar(first_data) == FALSE) {
      first_data <- readLines(con, n = 1L, warn = FALSE)
      if (!length(first_data)) break
    }

    # default: if no data, assume everything after header could be genotypes
    if (!length(first_data)) {
      lastmetric_idx  <- max(which(toupper(cn) %in% c("SNP","ALLELEID")), na.rm = TRUE)
      if (!is.finite(lastmetric_idx)) lastmetric_idx <- min(3L, length(cn))
      return(list(topskip = topskip,
                  lastmetric_idx = lastmetric_idx,
                  lastmetric_name = cn[lastmetric_idx],
                  cn = cn))
    }

    toks <- strsplit(first_data, delim, fixed = TRUE)[[1]]
    if (length(toks) < length(cn)) toks <- c(toks, rep("", length(cn) - length(toks)))
    toks <- toks[seq_len(length(cn))]

    is_gtlike <- function(x) {
      x <- toupper(trimws(x))
      if (!nzchar(x)) return(FALSE)
      grepl("^([012]|[012]/[012]|\\./\\.|N|NC|-|AA|TT|GG|CC|[ACGT](/[ACGT])?)$", x, perl = TRUE)
    }

    idx_snp <- {
      i <- match("SNP", cn); if (is.na(i)) match("SNP", toupper(cn)) else i
    }
    start_scan <- if (!is.na(idx_snp)) idx_snp + 1L else 3L

    first_gt <- NA_integer_
    for (j in seq.int(start_scan, length(cn))) {
      if (is_gtlike(toks[j])) { first_gt <- j; break }
    }

    if (is.na(first_gt)) {
      # fallback: last known metric-like column
      metric_candidates <- c("CallRate","Reproducibility","AvgCountRef","AvgCountSnp",
                             "OneRatioRef","OneRatioSnp","RepAvg","PIC","AvgPIC",
                             "Pavg","PICavg","AvgDP","AvgDepth","AvgMAPQ")
      mi <- which(toupper(cn) %in% toupper(metric_candidates))
      lastmetric_idx <- if (length(mi)) max(mi) else if (!is.na(idx_snp)) idx_snp else 2L
    } else {
      lastmetric_idx <- first_gt - 1L
    }

    lastmetric_idx  <- max(min(lastmetric_idx, length(cn) - 1L), 2L)
    lastmetric_name <- cn[lastmetric_idx]

    list(topskip = topskip,
         lastmetric_idx = lastmetric_idx,
         lastmetric_name = lastmetric_name,
         cn = cn)
  }

  parse_snp_ref_alt <- function(s) {
    # DArT SNP strings like "[A/G]" or "A/G" -> REF=A, ALT=G
    s <- toupper(trimws(gsub("[\\[\\](){} ]", "", s, perl = TRUE)))
    REF <- rep("N", length(s))
    ALT <- rep(".", length(s))

    ok <- grepl("^[ACGTN]+/[ACGTN]+$", s)
    if (any(ok)) {
      spl <- strsplit(s[ok], "/", fixed = TRUE)
      REF[ok] <- vapply(spl, `[[`, character(1), 1L)
      ALT[ok] <- vapply(spl, `[[`, character(1), 2L)
    }
    list(REF = REF, ALT = ALT)
  }

  nzchar_safe <- function(x) !is.na(x) & x != "" & x != "."
  is_blank    <- function(x) is.na(x) | x == "" | x == "."

  ## Layout detection
  lay <- detect_dart_layout_base(raw_path)
  message(sprintf("DArT layout: topskip=%s, lastmetric=%s (idx=%s)",
                  lay$topskip, lay$lastmetric_name, ifelse(is.na(lay$lastmetric_idx), "NA", lay$lastmetric_idx)))

  cn  <- lay$cn
  idx_allele <- match("AlleleID", cn)
  idx_snp    <- match("SNP",      cn)
  if (is.na(idx_allele) || is.na(idx_snp))
    stop("Could not find required columns 'AlleleID' and/or 'SNP' in DArT file header.")

  g_start  <- lay$lastmetric_idx + 1L
  if (g_start <= 0L || g_start > length(cn))
    stop("Computed genotype block start index is invalid.")
  geno_idx <- g_start:length(cn)

  ## Header
  fileDate     <- format(Sys.Date(), "%Y%m%d")
  VCFchroms    <- unique(c(order_contigs(Finalmappings$saccver), "chrUnk"))
  contig_lines <- unique(build_contig_lines(VCFchroms, ncbi_meta))
  extra_meta   <- build_extra_meta(ncbi_meta, genomename)

  header_lines <- c(
    "##fileformat=VCFv4.3",
    sprintf("##fileDate=%s", fileDate),
    "##source=Brioche-VCF-build",
    extra_meta,
    "##FILTER=<ID=LowQual,Description=\"Low quality or ambiguous mapping\">",
    '##INFO=<ID=MAPSTATUS,Number=1,Type=String,Description="Mapping status from anchoring_script.R (Unique_mapping|Failed_to_map_uniquely)">',
    '##INFO=<ID=ORIENTATIONSTATUS,Number=1,Type=String,Description="Whether the marker has been oriented w.r.t. reference alleles (Yes|No)">',
    "##INFO=<ID=MAF,Number=A,Type=Float,Description=\"Minor Allele frequency\">",
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">",
    "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    contig_lines
  )
  header_lines <- header_lines[!duplicated(header_lines)]
  writeLines(header_lines, con = Outputfilename)

  cn_fixed   <- make.unique(lay$cn, sep = "...")
  sample_ids <- cn_fixed[geno_idx]
  vcf_header <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", sample_ids)
  data.table::fwrite(as.list(vcf_header), file = Outputfilename, sep = "\t",
                     append = TRUE, col.names = FALSE, quote = FALSE)

  ## Mapping frame (unique-hit counts)
  fm <- Finalmappings %>%
    dplyr::transmute(
      Name      = qaccver,
      CHROM_map = saccver,
      POS_map   = SNPpos,
      REF_map   = toupper(Ref),
      ALT_map   = normalize_alt(ALT)
    )
  fm_hits <- Finalmappings %>% dplyr::count(qaccver, name = "MAP_HITS")
  fm <- fm %>%
    dplyr::left_join(fm_hits, by = c("Name" = "qaccver")) %>%
    dplyr::mutate(MAP_HITS = ifelse(is.na(MAP_HITS), 0L, MAP_HITS))
  fm_Name <- fm$Name

  ## Chunked reader (base R) — strictly biallelic
  chunk_size <- 1500L
  delim <- if (grepl("\\.(tsv|txt)$", raw_path, ignore.case = TRUE)) "\t" else ","
  con_dart <- if (grepl("\\.gz$", raw_path, ignore.case = TRUE)) gzfile(raw_path, "rt") else file(raw_path, "rt")
  on.exit(try(close(con_dart), silent = TRUE), add = TRUE)

  # Skip the preamble + header line that utils::read.* saw
  invisible(readLines(con_dart, n = lay$topskip + 1L, warn = FALSE))

  repeat {
    raw_lines <- readLines(con_dart, n = chunk_size, warn = FALSE)
    if (!length(raw_lines)) break

    df_all <- utils::read.table(text = raw_lines, sep = delim, header = FALSE,
                                quote = "", comment.char = "", stringsAsFactors = FALSE,
                                check.names = FALSE, col.names = cn_fixed, fill = TRUE)

    # keep only AlleleID, SNP and genotype columns
    df <- df_all[, c(idx_allele, idx_snp, geno_idx), drop = FALSE]

    Name <- df[[1L]]
    SNP  <- df[[2L]]
    geno <- as.matrix(df[, -(1:2), drop = FALSE])

    pa   <- parse_snp_ref_alt(SNP)
    REF_old <- toupper(pa$REF)
    ALT_old <- toupper(pa$ALT)         # already single ALT from DArT SNP field

    idx   <- match(Name, fm_Name)
    has_m <- !is.na(idx)

    hits   <- ifelse(has_m, fm$MAP_HITS[idx], 0L)
    u_map  <- has_m & hits == 1L

    CHROM <- ifelse(u_map, fm$CHROM_map[idx], "chrUnk")
    POS   <- ifelse(u_map, fm$POS_map[idx],   0L)

    ## Candidate Brioche targets (strictly one ALT)
    BriocheREF <- ifelse(u_map, toupper(fm$REF_map[idx]), NA_character_)
    BriocheALT <- ifelse(u_map, toupper(sub(",.*$", "", fm$ALT_map[idx])), NA_character_)
    BriocheREFRC <- ifelse(!is.na(BriocheREF), revcomp_iupac(BriocheREF), NA_character_)
    BriocheALTRC <- ifelse(!is.na(BriocheALT), revcomp_iupac(BriocheALT), NA_character_)

    # Biallelic viability: one ALT on both sides
    old_alt_len <- ifelse(is.na(ALT_old) | ALT_old == "" | ALT_old == ".", 0L, 1L)
    new_alt_len <- ifelse(is.na(BriocheALT) | BriocheALT == "" | BriocheALT == ".", 0L, 1L)
    bial_ok     <- u_map & (old_alt_len == 1L) & (new_alt_len == 1L)

    genotypeREF <- REF_old
    genotypeALT <- ALT_old

    ident    <- bial_ok & !is.na(BriocheREF)   & (BriocheREF   == genotypeREF)
    ident_rc <- bial_ok & !is.na(BriocheREFRC) & (BriocheREFRC == genotypeREF)
    swap     <- bial_ok & !is.na(BriocheREF)   & (BriocheREF   == genotypeALT)
    swap_rc  <- bial_ok & !is.na(BriocheREFRC) & (BriocheREFRC == genotypeALT)

    oriented_now <- ident | ident_rc | swap | swap_rc

    # Output REF/ALT: when oriented, force to Brioche REF/ALT; otherwise keep input
    REF_out <- ifelse(oriented_now, BriocheREF, REF_old)
    ALT_out <- ifelse(oriented_now, BriocheALT, ALT_old)
    ALT_out <- sub(",.*$", "", ALT_out)           # enforce single ALT
    ALT_out[is.na(ALT_out) | ALT_out == ""] <- "."

    ## FILTER/QUAL/INFO based on mapping
    pos_int     <- suppressWarnings(as.integer(POS))
    coord_bad   <- is.na(pos_int) | pos_int <= 0L | is_blank(CHROM)
    unknown_now <- !u_map | coord_bad

    QUAL   <- ifelse(unknown_now, 0L, 100L)
    FILTER <- ifelse(unknown_now, "LowQual", "PASS")
    INFO   <- ifelse(
      unknown_now,
      "MAPSTATUS=Failed_to_map_uniquely;ORIENTATIONSTATUS=No",
      ifelse(oriented_now, "MAPSTATUS=Unique_mapping;ORIENTATIONSTATUS=Yes",
                        "MAPSTATUS=Unique_mapping;ORIENTATIONSTATUS=No")
    )

    ## DArT genotypes recode (0,1,2, "./.", "N","NC","-")
    geno[] <- toupper(gsub("\\|", "/", geno, perl = TRUE))
    miss <- is.na(geno) | geno == "" | geno == "." | geno == "-" | geno == "N" | geno == "NC"
    if (any(miss)) geno[miss] <- "./."
    dbl <- grepl("^[012]{2}$", geno, perl = TRUE)
    if (any(dbl)) geno[dbl] <- sub("^([012])([012])$", "\\1/\\2", geno[dbl], perl = TRUE)

    # SAME/RC-SAME: no GT change; SWAP/RC-SWAP: flip 0<->2, keep hetero
    rsame <- ifelse(is.na(ident | ident_rc), FALSE, (ident | ident_rc))
    rswap <- ifelse(is.na(swap  | swap_rc ), FALSE, (swap  | swap_rc ))

    if (any(rsame)) {
      M_same <- matrix(rep(rsame, times = ncol(geno)), nrow = nrow(geno), ncol = ncol(geno))
      g <- geno
      d0 <- (g == "0") & M_same;  if (any(d0)) g[d0] <- "0/0"
      d1 <- (g == "1") & M_same;  if (any(d1)) g[d1] <- "0/1"
      d2 <- (g == "2") & M_same;  if (any(d2)) g[d2] <- "1/1"
      geno[M_same] <- g[M_same]
    }
    if (any(rswap)) {
      M_swap <- matrix(rep(rswap, times = ncol(geno)), nrow = nrow(geno), ncol = ncol(geno))
      g <- geno
      d0 <- (g == "0") & M_swap;  if (any(d0)) g[d0] <- "1/1"
      d1 <- (g == "1") & M_swap;  if (any(d1)) g[d1] <- "0/1"
      d2 <- (g == "2") & M_swap;  if (any(d2)) g[d2] <- "0/0"
      geno[M_swap] <- g[M_swap]
    }

    # Normalize any lingering single digits
    i0 <- geno == "0"; if (any(i0)) geno[i0] <- "0/0"
    i1 <- geno == "1"; if (any(i1)) geno[i1] <- "0/1"
    i2 <- geno == "2"; if (any(i2)) geno[i2] <- "1/1"

    ## AC/AN/MAF (strictly biallelic, count allele '1')
    nr <- nrow(geno); ac1 <- integer(nr); an <- integer(nr)
    for (j in seq_len(ncol(geno))) {
      x <- geno[, j]
      m <- (x == "./.")
      a1 <- substr(x, 1L, 1L); a1[m] <- NA_character_
      a2 <- substr(x, 3L, 3L); a2[m] <- NA_character_
      ac1 <- ac1 + as.integer((!is.na(a1)) & (a1 == "1")) + as.integer((!is.na(a2)) & (a2 == "1"))
      an  <- an  + as.integer(!m) * 2L
    }
    AC_str  <- as.character(ac1)
    AF1     <- ifelse(an > 0L, ac1 / an, NA_real_)
    MAF_str <- sprintf("%.6g", AF1)
    INFO <- ifelse(
      is.na(INFO) | INFO == ".",
      paste0("MAF=", MAF_str, ";AC=", AC_str, ";AN=", an),
      paste0(INFO, ";MAF=", MAF_str, ";AC=", AC_str, ";AN=", an)
    )

    body <- data.table::data.table(
      `#CHROM` = CHROM, POS = POS, ID = Name,
      REF = REF_out, ALT = ALT_out, QUAL = QUAL, FILTER = FILTER, INFO = INFO, FORMAT = "GT"
    )
    body <- cbind(body, data.table::as.data.table(geno))

    data.table::fwrite(body, file = Outputfilename, sep = "\t",
                       append = TRUE, col.names = FALSE, quote = FALSE)

    rm(df_all, df, geno, body, ac1, an); gc(FALSE)
  }

  message(sprintf("Wrote: %s", Outputfilename))
  quit(save = "no")
}



# Come back to this section to add in more diverse datatypes for non SNP chip output. 
if(!is_snpchip & !is_vcfraws & !isDArTfile) {
  
  message(sprintf("Wrote: %s", "Hello world"))
  
}

