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
#     --genomename ""\
#     --reference_species ""\
#     --reference_url ""\
#     --droplist ""\
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
#   - --reference_species "" :User can input a specific reference species to be added to the header metadata of the vcf file. (Useful for when not working with an NCBI reference)
#   - --reference_url "" :User can input a specific reference url to be added to the header metadata of the vcf file. (Useful for when the genome was downloaded from a repository that isn't NCBI)
#   - --droplist "" :user can provide drop list of markernames which should be changed to unknown Chromosome and position [allow for duplicates to be dealt with]. 


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
Outputfilename<- get_opt(opts, "Outputfilename", "outfile")
droplist_path <- get_opt(opts, "droplist",   "")
Accession     <- get_opt(opts, "Accession",  "")
genomename    <- get_opt(opts, "genomename", "")
user_species   <- get_opt(opts, "reference_species", "")
user_genome_url<- get_opt(opts, "reference_url",   "")

message(sprintf("Rawgenotypes   : %s", raw_path))
message(sprintf("Briochemappings: %s", map_path))
message(sprintf("IsVCFraws      : %s", is_vcfraws))
message(sprintf("IsSNPchip      : %s", is_snpchip))
message(sprintf("isDArTfile     : %s", isDArTfile))
message(sprintf("Outputfilename : %s", Outputfilename))
message(sprintf("Accession      : %s", Accession))
message(sprintf("genomename     : %s", genomename))
message(sprintf("reference_species: %s", user_species))
message(sprintf("reference_url    : %s", user_genome_url))
message(sprintf("droplist      : %s", droplist_path))

# Functions

nzchar_safe <- function(x) !is.na(x) & x != "" & x != "."

# ALT CSV utilities
split_alts <- function(x) {
  x <- as.character(x)
  lapply(x, function(s) if (is.na(s) || s == "" || s == ".") character(0) else strsplit(s, ",", fixed = TRUE)[[1]])
}
join_alts <- function(v) if (!length(v)) "." else paste(v, collapse = ",")

# Build 0/1/2 index remap given old (REF + ALTs) and new (REF and ALTs)
# Returns an integer vector mapping old indices new indices (0-based)
build_symbol_maps <- function(ref_old, alts_old_chr, ref_new, alts_new_chr) {
  alts_old_chr <- if (is.null(alts_old_chr)) character(0) else alts_old_chr
  alts_new_chr <- if (is.null(alts_new_chr)) character(0) else alts_new_chr
  sym_old <- c(ref_old, alts_old_chr)
  sym_new <- c(ref_new, alts_new_chr)
  idx_map <- integer(length(sym_old))
  for (i in seq_along(sym_old)) {
    m <- match(toupper(sym_old[i]), toupper(sym_new))
    idx_map[i] <- ifelse(is.na(m), NA_integer_, m - 1L) 
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

info_get <- function(info_vec, key){
  x <- as.character(info_vec); x[is.na(x)] <- ""
  m <- regexec(paste0("(^|;)", key, "=([^;]+)"), x, perl = TRUE)
  regmatches(x, m) |>
    lapply(function(hit) if (length(hit) >= 3) hit[3] else NA_character_) |>
    unlist(use.names = FALSE)
}

info_strip_keys <- function(info_vec, keys){
  x <- as.character(info_vec); x[is.na(x)] <- ""
  for (k in keys) x <- gsub(paste0("(^|;)", k, "=[^;]*"), "\\1", x, perl = TRUE)
  x <- gsub(";{2,}", ";", x, perl = TRUE)
  gsub("^;|;$", "", x, perl = TRUE)
}

norm_prior <- function(x){
  x <- tolower(trimws(as.character(x)))
  ifelse(x %in% c("+","plus"), "plus",
         ifelse(x %in% c("-","minus"), "minus", "none"))
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

load_drop_ids <- function(path) {
  if (!nzchar(path)) return(character(0))
  if (!file.exists(path)) {
    warning(sprintf("[droplist] File not found: %s (ignored)", path))
    return(character(0))
  }
  # Flexible reader: try readr first, fall back to base if needed
  read_lines_safely <- function(p) {
    tryCatch(readLines(p, warn = FALSE), error = function(e) character(0))
  }
  raw <- read_lines_safely(path)
  if (!length(raw)) return(character(0))
  # Strip BOM, trim, drop comments/empties
  raw <- sub("^\ufeff", "", raw, perl = TRUE)
  raw <- trimws(raw)
  raw <- raw[!grepl("^\\s*#", raw)]         # drop comment lines
  raw <- raw[nzchar(raw)]                   # drop empty lines
  # If file has tabs/commas and more than 1 col, keep only first field
  first_field <- sub("[\t,].*$", "", raw)
  unique(first_field[nzchar(first_field)])
}

DROP_IDS <- load_drop_ids(droplist_path)

# Fast membership test function
`%in_drop%` <- if (requireNamespace("fastmatch", quietly = TRUE)) {
  function(x, y) !is.na(fastmatch::fmatch(x, y))
} else {
  function(x, y) x %in% y
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
      if (is.list(sums) && !is.null(sums$uid)) list(sums)
      else if (is.list(sums)) sums
      else list(sums)
    }

    if (nzchar(acc_in)) {
      ids <- character(0)
      q1 <- sprintf('%s[Assembly Accession]', acc_in)
      s1 <- rentrez::entrez_search(db = "assembly", term = q1, retmax = 20)
      ids <- unique(c(ids, s1$ids))
      if (!length(ids)) {
        q2 <- sprintf('"%s"[Assembly Accession]', acc_in)
        s2 <- rentrez::entrez_search(db = "assembly", term = q2, retmax = 20)
        ids <- unique(c(ids, s2$ids))
      }
      if (!length(ids)) {
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
      data   <- data[!grepl("^#", data)]
      txt    <- paste(c(header, data), collapse = "\n")
      readr::read_tsv(I(txt), show_col_types = FALSE, progress = FALSE,
                      col_types = readr::cols(.default = "c"))
    }

    ar <- parse_assembly_report(rpt_url)

    # Flexible column lookups
    col_len      <- intersect(c("Sequence-Length","sequence-length","length"), names(ar))
    col_assigned <- intersect(c("Assigned-Molecule","assigned-molecule"), names(ar))
    col_md5      <- intersect(c("MD5","md5"), names(ar))
    col_seqname  <- intersect(c("Sequence-Name","sequence-name"), names(ar))
    col_role     <- intersect(c("Sequence-Role","sequence-role"), names(ar))
    col_refseq   <- intersect(c("RefSeq-Accn","refseq-accn","RefSeq Accn"), names(ar))
    col_genbank  <- intersect(c("GenBank-Accn","genbank-accn","GenBank Accn"), names(ar))

    # 1) Prefer Sequence-Name as key; fallback to Assigned-Molecule
    key_col <- if (length(col_seqname)) col_seqname[1] else if (length(col_assigned)) col_assigned[1] else NULL
    if (is.null(key_col)) stop("assembly_report lacks Sequence-Name / Assigned-Molecule columns")

    seq_name <- ar[[key_col]]
    assigned <- if (length(col_assigned)) ar[[col_assigned[1]]] else rep(NA_character_, nrow(ar))
    len_val  <- suppressWarnings(as.numeric(ar[[col_len[1]]]))
    md5_val  <- if (length(col_md5)) ar[[col_md5[1]]] else rep(NA_character_, nrow(ar))
    role_val <- if (length(col_role)) ar[[col_role[1]]] else rep(NA_character_, nrow(ar))
    ref_acc  <- if (length(col_refseq)) ar[[col_refseq[1]]] else rep(NA_character_, nrow(ar))
    gb_acc   <- if (length(col_genbank)) ar[[col_genbank[1]]] else rep(NA_character_, nrow(ar))

    # numeric chromosome index (use Sequence-Name first, else Assigned-Molecule)
    extract_num <- function(x) {
      s <- tolower(trimws(x))
      m <- stringr::str_match(s, "^.*(?:chromosome|chrom|chr|contig|scaffold|segment|seg)?[ _-]*0*([0-9]+)([a-z])?(?![a-z0-9]).*")
      suppressWarnings(as.integer(m[,2]))
    }

    chrom_num <- extract_num(seq_name)
    chrom_num <- ifelse(
      is.na(chrom_num),
      extract_num(assigned),
      chrom_num
    )
    preferred <- ifelse(nzchar(ref_acc), ref_acc, gb_acc)
    nuccore   <- ifelse(nzchar(preferred), paste0("https://www.ncbi.nlm.nih.gov/nuccore/", preferred), NA_character_)

    # keep rows that actually have a length and some label
    keep <- !is.na(len_val) & nzchar(seq_name)
    chr_table <- data.frame(
      chrom_num        = chrom_num,
      seq_name         = seq_name,
      assigned_molecule= assigned,
      sequence_role    = role_val,
      refseq_accn      = ref_acc,
      genbank_accn     = gb_acc,
      preferred_accn   = preferred,
      length           = len_val,
      md5              = md5_val,
      nuccore_url      = nuccore,
      stringsAsFactors = FALSE
    )[keep, , drop = FALSE]

    # name lengths / md5 by Sequence-Name (lowercased)
    key_seq <- tolower(chr_table$seq_name)
    len_by_seqname <- stats::setNames(chr_table$length, key_seq)
    md5_by_seqname <- stats::setNames(chr_table$md5,    key_seq)

    list(
      ok               = TRUE,
      accession        = acc_sum,
      asm_name         = aname,
      organism         = org,
      taxid            = as.integer(taxid),
      asm_url          = asm_url,
      chr_table        = chr_table,
      len_by_seqname   = len_by_seqname,
      md5_by_seqname   = md5_by_seqname
    )
  }, error = function(e) {
    message("Accession or genome name were nonspecific unable to populate genomic information")
    out
  })
}

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

build_extra_meta <- function(ncbi_meta, genomename, override_species = NULL, override_url = NULL) {
  esc <- function(x) gsub('"', '\\\\"', x, fixed = TRUE)

  gn <- trimws(genomename %||% "")
  # Species preference order:
  #   1) NCBI organism (when ncbi_meta OK)
  #   2) override_species (user-provided)
  #   3) genomename (pipeline input)
  #   4) "unknown"
  species_fallback <- function() {
    s <- override_species %||% gn %||% "unknown"
    trimws(s)
  }

  if (isTRUE(ncbi_meta$ok)) {
    sp <- (ncbi_meta$organism %||% species_fallback())
    acc <- ncbi_meta$accession %||% (gn %||% "unknown")
    asm_url <- ncbi_meta$asm_url %||% ""

    c(
      sprintf("##reference=%s", acc),
      sprintf(
        '##assembly=<accession=%s,name="%s"%s%s%s>',
        acc,
        esc(ncbi_meta$asm_name %||% "unknown"),
        if (nzchar(sp)) sprintf(',species="%s"', esc(sp)) else "",
        if (!is.null(ncbi_meta$taxid) && !is.na(ncbi_meta$taxid)) sprintf(",taxonomy=%s", ncbi_meta$taxid) else "",
        if (nzchar(asm_url)) sprintf(',url="%s"', esc(asm_url)) else ""
      ),
      # Convenience: also surface the URL as a simple key=value line
      if (nzchar(asm_url)) sprintf("##genome_url=%s", asm_url) else NULL
    )
  } else {
    # No NCBI metadata — fall back to user-provided species/url (or unknowns)
    sp <- species_fallback()
    url_out <- override_url %||% ""  # only write this when NCBI meta is absent

    ref_val <- gn %||% sp %||% "unknown"
    c(
      sprintf("##reference=%s", ref_val),
      sprintf(
        '##assembly=<name="%s",species="%s"%s>',
        esc(gn %||% "unknown"),
        esc(sp %||% "unknown"),
        if (nzchar(url_out)) sprintf(',url="%s"', esc(url_out)) else ""
      ),
      if (nzchar(url_out)) sprintf("##genome_url=%s", url_out) else NULL
    )
  }
}

build_contig_lines <- function(ids, meta, override_species = NULL) {
  ids <- unique(na.omit(as.character(ids)))
  esc <- function(x) gsub('"', '\\\\"', x, fixed = TRUE)

  # choose species string for contigs
  contig_species <- function() {
    if (isTRUE(meta$ok) && !is.null(meta$organism) && nzchar(meta$organism)) {
      meta$organism
    } else {
      override_species %||% "unknown"
    }
  }

  mk_attr <- function(id, len = NA_real_, acc = NA_character_, alias = NA_character_, tax = NA) {
    parts <- c(sprintf("ID=%s", id))
    if (!is.na(len)) parts <- c(parts, sprintf("length=%d", as.integer(len)))
    parts <- c(parts, sprintf('species="%s"', esc(contig_species())))
    # Taxonomy: integer when known, otherwise literal unknown (kept unquoted)
    parts <- c(parts, if (!is.null(tax) && !is.na(tax)) sprintf("taxonomy=%s", tax) else "taxonomy=unknown")
    parts <- c(parts, sprintf("Accession=%s", if (!is.na(acc) && nzchar(acc)) acc else "unknown"))
    parts <- c(parts, sprintf('Alias="%s"', if (!is.na(alias) && nzchar(alias)) esc(alias) else "unknown"))
    sprintf("##contig=<%s>", paste(parts, collapse = ","))
  }

  # If we have no chrom table, produce unknown-but-complete lines
  if (!isTRUE(meta$ok) || is.null(meta$chr_table) || !NROW(meta$chr_table)) {
    tax <- if (!is.null(meta$taxid) && !is.na(meta$taxid)) meta$taxid else NA
    return(vapply(ids, function(id) mk_attr(id, tax = tax), FUN.VALUE = character(1)))
  }

  ct <- meta$chr_table
  n  <- nrow(ct)

  # Lowercased keys for matching
  key_seq  <- tolower(trimws(if ("seq_name" %in% names(ct)) ct$seq_name else rep(NA_character_, n)))
  key_assn <- tolower(trimws(if ("assigned_molecule" %in% names(ct)) ct$assigned_molecule else rep(NA_character_, n)))

  # Named lookups (safe; missing -> NA)
  len_seq   <- setNames(suppressWarnings(as.numeric(if ("length" %in% names(ct)) ct$length else rep(NA_real_, n))), key_seq)
  len_ass   <- setNames(suppressWarnings(as.numeric(if ("length" %in% names(ct)) ct$length else rep(NA_real_, n))), key_assn)
  acc_seq   <- setNames(if ("preferred_accn" %in% names(ct)) ct$preferred_accn else rep(NA_character_, n), key_seq)
  acc_ass   <- setNames(if ("preferred_accn" %in% names(ct)) ct$preferred_accn else rep(NA_character_, n), key_assn)
  alias_seq <- setNames(if ("assigned_molecule" %in% names(ct)) ct$assigned_molecule else rep(NA_character_, n), key_seq)
  alias_ass <- setNames(if ("assigned_molecule" %in% names(ct)) ct$assigned_molecule else rep(NA_character_, n), key_assn)

  tax <- if (!is.null(meta$taxid) && !is.na(meta$taxid)) meta$taxid else NA

  vapply(ids, function(id) {
    lid <- tolower(trimws(id))

    # Prefer Sequence-Name ? fall back to Assigned-Molecule
    len   <- if (lid %in% names(len_seq))   len_seq[lid]   else if (lid %in% names(len_ass))   len_ass[lid]   else NA_real_
    acc   <- if (lid %in% names(acc_seq)   && nzchar(acc_seq[lid]))   acc_seq[lid]   else if (lid %in% names(acc_ass))   acc_ass[lid]   else NA_character_
    alias <- if (lid %in% names(alias_seq) && nzchar(alias_seq[lid])) alias_seq[lid] else if (lid %in% names(alias_ass)) alias_ass[lid] else NA_character_

    mk_attr(id, len = len, acc = acc, alias = alias, tax = tax)
  }, FUN.VALUE = character(1))
}


# If raws genotypes is not a vcf file  
if (!is_vcfraws & !isDArTfile) {
  
  genotypesRAWS  <- read.table(raw_path, header = TRUE, sep = ifelse(grepl("\\.tsv$", raw_path, ignore.case=TRUE), "\t", ""),
                               check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
  
}


Finalmappings  <- read.csv(map_path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "#")


# Newly added metadata in lines 1-10 of Brioche (might expand further)
# Ln 1= generic, Ln2 brioche commit version, Ln3 brioche repo, Ln4 brioche branch, ln5 pident ln6 coverage, ln7 blast options, ln8 target chromprior? ln9 sharedmarkersmap used?,ln10 ld mapping used?, ln11 genetic map used?, ln 12 rundate
Metadata <- readLines(map_path,n=12L)
# remove excess comment chars 
Metadata <- sub(pattern="## ", replacement="",x=Metadata, fixed =TRUE)
# Remove extra column at start
Metadata <- Metadata[-1]



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
    message("NCBI metadata   : none found (proceeding with user/unknown fallbacks)")
  }
}

# if NCBI failed and user gave species
if (!isTRUE(ncbi_meta$ok) && nzchar(user_species)) {
  ncbi_meta$organism <- user_species
}

# Header meta lines (uses NCBI URL when present; otherwise uses user_genome_url)
extra_meta <- build_extra_meta(
  ncbi_meta,
  genomename,
  override_species = if (nzchar(user_species)) user_species else NULL,
  override_url     = if (nzchar(user_genome_url)) user_genome_url else NULL
)

#Contig lines (always full-shape; uses overrides only when NCBI meta is missing)
VCFchroms <- order_contigs(Finalmappings$saccver)
contig_lines <- build_contig_lines(
  ids   = VCFchroms,
  meta  = ncbi_meta,
  override_species = if (nzchar(user_species)) user_species else NULL
)



# SNP-chip pathway (raw genotypes table to vcf)

if (!is_vcfraws & is_snpchip) {

  # Helpers
  strand_to_word <- function(x) {
    z <- tolower(trimws(as.character(x)))
    ifelse(is.na(z) | z == "" | z == ".", "none",
      ifelse(z %in% c("+","plus","pos","positive","p"), "plus",
        ifelse(z %in% c("-","minus","neg","negative","m"), "minus", "none")))
  }
  is_blank <- function(x) is.na(x) | x == "" | x == "."

  
  delim <- if (grepl("\\.(tsv|txt)$", raw_path, ignore.case = TRUE)) "\t" else ","
  hdr0  <- utils::read.table(
    file = raw_path, header = TRUE, sep = delim, nrows = 0,
    check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE
  )
  cn_raw   <- names(hdr0)
  cn_fixed <- make.unique(cn_raw, sep = "...")

  idx_Name <- match("Name", cn_raw)
  idx_REF  <- match("REF",  cn_raw)
  idx_ALT  <- match("ALT",  cn_raw)
  if (any(is.na(c(idx_Name, idx_REF, idx_ALT))))
    stop("SNP chip file must have columns: Name, REF, ALT.")

  sample_idx <- setdiff(seq_along(cn_raw), c(idx_Name, idx_REF, idx_ALT))
  if (!length(sample_idx)) stop("No sample genotype columns found after ALT.")
  sample_ids <- cn_fixed[sample_idx]

  ## Mapping frame (include sstrand)
  fm_base <- Finalmappings %>%
    dplyr::transmute(
      Name      = qaccver,
      CHROM_map = saccver,
      POS_map   = SNPpos,
      REF_map   = toupper(Ref),
      ALT_map   = normalize_alt(ALT),
      STRAND    = sstrand
    )
  fm_hits <- Finalmappings %>% dplyr::count(qaccver, name = "MAP_HITS")
  fm <- fm_base %>%
    dplyr::left_join(fm_hits, by = c("Name" = "qaccver")) %>%
    dplyr::mutate(MAP_HITS = ifelse(is.na(MAP_HITS), 0L, MAP_HITS))
  fm_Name <- fm$Name

  ## ---- 1) Write VCF header ----
  fileDate     <- format(Sys.Date(), "%Y%m%d")
  VCFchroms    <- unique(c(order_contigs(Finalmappings$saccver), "chrUnk"))
  contig_lines #<- unique(build_contig_lines(VCFchroms, ncbi_meta))
  extra_meta   #<- build_extra_meta(ncbi_meta, genomename)

header_lines <- c(
  "##fileformat=VCFv4.3",
  sprintf("##fileDate=%s", fileDate),
  "##source=Brioche-VCF-build",
  paste0("##",Metadata[1]),
  paste0("##",Metadata[2]),
  paste0("##",Metadata[3]),
  paste0("##Brioche_filtering_",Metadata[4]),
  paste0("##Brioche_filtering_",Metadata[5]),
  paste0("##Brioche_filtering_",Metadata[6]),
  paste0("##Brioche_priors_files:",Metadata[7]),
  paste0("##Brioche_priors_files:",Metadata[8]),
  paste0("##Brioche_priors_files:",Metadata[9]),
  paste0("##Brioche_priors_files:",Metadata[10]),
  paste0("##Brioche_",Metadata[11]),
  extra_meta,
  "##FILTER=<ID=LowQual,Description=\"Low quality or ambiguous mapping\">",
  '##INFO=<ID=MAPSTATUS,Number=1,Type=String,Description="Unique_mapping|Failed_to_map_uniquely">',
  '##INFO=<ID=PriorORIENT,Number=1,Type=String,Description="Accumulated strand orientation across runs (plus|minus|none)">',
  "##INFO=<ID=MAF,Number=A,Type=Float,Description=\"Minor Allele frequency per ALT\">",
  "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count per ALT\">",
  "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "##FORMAT=<ID=NU,Number=1,Type=Integer,Description=\"Null allele flag from SNP chip (1= Null allele; 0=otherwise)\">",
  contig_lines
)
  header_lines <- header_lines[!duplicated(header_lines)]
  writeLines(header_lines, con = Outputfilename)

  vcf_header <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", sample_ids)
  data.table::fwrite(as.list(vcf_header), file = Outputfilename, sep = "\t",
                     append = TRUE, col.names = FALSE, quote = FALSE)

  ## ---- 2) Stream rows ----
  chunk_size <- 1500L
  con_chunk <- if (grepl("\\.gz$", raw_path, ignore.case = TRUE)) gzfile(raw_path, "rt") else file(raw_path, "rt")
  on.exit(try(close(con_chunk), silent = TRUE), add = TRUE)
  invisible(readLines(con_chunk, n = 1L, warn = FALSE))  # consume header

  repeat {
    raw_lines <- readLines(con_chunk, n = chunk_size, warn = FALSE)
    if (!length(raw_lines)) break

    df <- utils::read.table(text = raw_lines, sep = delim, header = FALSE,
                            quote = "", comment.char = "", stringsAsFactors = FALSE,
                            check.names = FALSE, col.names = cn_fixed, fill = TRUE)

    Name    <- df[[idx_Name]]
    force_drop <- if (length(DROP_IDS)) Name %in_drop% DROP_IDS else rep(FALSE, length(Name))
    REF_old <- toupper(df[[idx_REF]])
    ALT_old <- normalize_alt(df[[idx_ALT]])
    ALT_old <- sub(",.*$", "", ALT_old)  # enforce single ALT

    m     <- match(Name, fm_Name)
    hasM  <- !is.na(m)
    hits  <- ifelse(hasM, fm$MAP_HITS[m], 0L)
    u_map <- hasM & hits == 1L & !force_drop

    CHROM <- ifelse(u_map, fm$CHROM_map[m], "chrUnk")
    POS   <- ifelse(u_map, fm$POS_map[m],   0L)

    BriocheREF <- ifelse(u_map, toupper(fm$REF_map[m]), NA_character_)
    BriocheALT <- ifelse(u_map, toupper(sub(",.*$", "", fm$ALT_map[m])), NA_character_)

    # Orientation decision (same as before for SNPchip), but compute PriorORIENT from sstrand
    BriocheREFRC <- ifelse(!is.na(BriocheREF), revcomp_iupac(BriocheREF), NA_character_)
    genotypeREF  <- REF_old
    genotypeALT  <- ALT_old

    ident    <- u_map & !is.na(BriocheREF)   & (BriocheREF   == genotypeREF)
    ident_rc <- u_map & !is.na(BriocheREFRC) & (BriocheREFRC == genotypeREF)
    swap     <- u_map & !is.na(BriocheREF)   & (BriocheREF   == genotypeALT)
    swap_rc  <- u_map & !is.na(BriocheREFRC) & (BriocheREFRC == genotypeALT)

    oriented_now <- ident | ident_rc | swap | swap_rc

    REF_out <- ifelse(oriented_now, BriocheREF, REF_old)
    ALT_out <- ifelse(oriented_now, BriocheALT, ALT_old)
    ALT_out <- sub(",.*$", "", ALT_out)
    ALT_out[is.na(ALT_out) | ALT_out == ""] <- "."

    # QUAL/FILTER/INFO
    pos_int     <- suppressWarnings(as.integer(POS))
    coord_bad   <- is.na(pos_int) | pos_int <= 0L | is_blank(CHROM)
    unknown_now <- !u_map | coord_bad

    QUAL   <- ifelse(unknown_now, 0L, 100L)
    FILTER <- ifelse(unknown_now, "LowQual", "PASS")
    INFO   <- ifelse(unknown_now, "MAPSTATUS=Failed_to_map_uniquely", "MAPSTATUS=Unique_mapping")

    # PriorORIENT from current sstrand when uniquely mapped; else 'none'
    cur_strand  <- ifelse(u_map, strand_to_word(fm$STRAND[m]), "none")
    INFO        <- paste0(INFO, ";PriorORIENT=", cur_strand)

    ## Genotypes + NU
    geno_raw <- as.matrix(df[, sample_idx, drop = FALSE])
    geno_raw_upper <- toupper(trimws(geno_raw))
    NU <- (geno_raw_upper == "N"); storage.mode(NU) <- "integer"; NU[is.na(NU)] <- 0L

    geno <- as.matrix(df[, sample_idx, drop = FALSE])
    geno[] <- toupper(gsub("\\|", "/", geno, perl = TRUE))
    miss <- is.na(geno) | geno == "" | geno == "." | geno == "-" | geno == "N" | geno == "NC"
    if (any(miss)) geno[miss] <- "./."
    dbl <- grepl("^[012]{2}$", geno, perl = TRUE)
    if (any(dbl)) geno[dbl] <- sub("^([012])([012])$", "\\1/\\2", geno[dbl], perl = TRUE)

    # SAME/RC-SAME: none; SWAP/RC-SWAP: 0<->2 flip
    rsame <- ifelse(is.na(ident | ident_rc), FALSE, (ident | ident_rc))
    rswap <- ifelse(is.na(swap  | swap_rc ), FALSE, (swap  | swap_rc ))
    if (any(rswap)) {
      g <- geno
      d0 <- (g == "0") & matrix(rep(rswap, times = ncol(geno)), nrow = nrow(geno))
      d1 <- (g == "1") & matrix(rep(rswap, times = ncol(geno)), nrow = nrow(geno))
      d2 <- (g == "2") & matrix(rep(rswap, times = ncol(geno)), nrow = nrow(geno))
      if (any(d0)) g[d0] <- "1/1"
      if (any(d1)) g[d1] <- "0/1"
      if (any(d2)) g[d2] <- "0/0"
      geno[matrix(rep(rswap, times = ncol(geno)), nrow = nrow(geno))] <- g[matrix(rep(rswap, times = ncol(geno)), nrow = nrow(geno))]
    }
    i0 <- geno == "0"; if (any(i0)) geno[i0] <- "0/0"
    i1 <- geno == "1"; if (any(i1)) geno[i1] <- "0/1"
    i2 <- geno == "2"; if (any(i2)) geno[i2] <- "1/1"

    ## AC/AN/MAF (biallelic)
    nr <- nrow(geno); ac1 <- integer(nr); an <- integer(nr)
    for (jj in seq_len(ncol(geno))) {
      x <- geno[, jj]
      msk <- x == "./."
      a1 <- substr(x, 1L, 1L); a1[msk] <- NA_character_
      a2 <- substr(x, 3L, 3L); a2[msk] <- NA_character_
      ac1 <- ac1 + as.integer((!is.na(a1)) & (a1 == "1")) + as.integer((!is.na(a2)) & (a2 == "1"))
      an  <- an  + as.integer(!msk) * 2L
    }
    AC_str  <- as.character(ac1)
    AF1     <- ifelse(an > 0L, ac1 / an, NA_real_)
    MAF_str <- sprintf("%.6g", AF1)
    INFO <- paste0(INFO, ";MAF=", MAF_str, ";AC=", AC_str, ";AN=", an)

    ## Write
    NU_chr <- matrix(as.character(NU), nrow = nrow(NU), ncol = ncol(NU)); NU_chr[is.na(NU_chr)] <- "0"
    geno_df <- as.data.frame(geno, stringsAsFactors = FALSE, check.names = FALSE)
    for (jj in seq_len(ncol(geno_df))) geno_df[[jj]] <- paste0(geno_df[[jj]], ":", NU_chr[, jj])

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
# VCF pathway (VCF -> re-anchored VCF using Brioche sstrand and PriorORIENT)
if (is_vcfraws) {

  # Helpers
  strand_to_word <- function(x) {
    z <- tolower(trimws(as.character(x)))
    ifelse(is.na(z) | z == "" | z == ".", NA_character_,
      ifelse(z %in% c("+","plus","pos","positive","p"), "plus",
        ifelse(z %in% c("-","minus","neg","negative","m"), "minus", NA_character_)))
  }
  pr_norm <- function(x) {
    z <- tolower(trimws(as.character(x)))
    ifelse(is.na(z) | z == "" | z == ".", "none",
      ifelse(z %in% c("plus","+","pos"), "plus",
        ifelse(z %in% c("minus","-","neg"), "minus", "none")))
  }
  is_blank <- function(x) is.na(x) | x == "" | x == "."

  split_alts_csv <- function(s) {
    if (is.na(s) || s == "" || s == ".") character(0) else strsplit(s, ",", fixed = TRUE)[[1]]
  }
  rc_csv <- function(csv) {
    v <- split_alts_csv(csv)
    if (!length(v)) return(csv)
    paste(revcomp_iupac(v), collapse = ",")
  }

  # Mapping frame with sstrand
  req <- c("qaccver","saccver","SNPpos","Ref","ALT","sstrand")
  missing <- setdiff(req, names(Finalmappings))
  if (length(missing))
    stop(sprintf("Finalmappings missing: %s", paste(missing, collapse=", ")))

  fm_one <- Finalmappings %>%
    transmute(qaccver, saccver, SNPpos,
              REF_map = toupper(Ref),
              ALT_map = normalize_alt(ALT),
              STRAND  = sstrand) %>%
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
      STRAND    = STRAND,
      MAP_HITS  = ifelse(is.na(MAP_HITS), 0L, MAP_HITS)
    )

  # Read input header
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

  fileDate     <- format(Sys.Date(), "%Y%m%d")

  # New header. Remove old brioche metadata, remove all assembly data, remove old file date, replace with newest. deduplicate etc
  meta_lines  <- meta_lines[!grepl("^##(fileDate=)", meta_lines)] 
  meta_lines  <- meta_lines[!grepl("^##(Brioche_)", meta_lines)]
  meta_lines  <- meta_lines[!grepl("^##(assembly=)", meta_lines)] 
  meta_lines  <- meta_lines[!grepl("^##(reference=)", meta_lines)] 
  meta_lines  <- meta_lines[!grepl("^##(genome_url=)", meta_lines)] 
  meta_lines  <- meta_lines[!grepl("^##(contig=)", meta_lines)]
  meta_lines  <- meta_lines[!grepl("^##(fileformat=)", meta_lines)]


    meta_lines <-c(
      "##fileformat=VCFv4.3",
      sprintf("##fileDate=%s", fileDate),
      meta_lines,
      paste0("##",Metadata[1]),
      paste0("##",Metadata[2]),
      paste0("##",Metadata[3]),
      paste0("##Brioche_filtering_",Metadata[4]),
      paste0("##Brioche_filtering_",Metadata[5]),
      paste0("##Brioche_filtering_",Metadata[6]),
      paste0("##Brioche_priors_files:",Metadata[7]),
      paste0("##Brioche_priors_files:",Metadata[8]),
      paste0("##Brioche_priors_files:",Metadata[9]),
      paste0("##Brioche_priors_files:",Metadata[10]),
      paste0("##Brioche_",Metadata[11]))


  meta_out   <- meta_out[!duplicated(meta_out)]
  meta_out <- meta_out[!grepl('^##FORMAT=<ID=ORIENTATIONSTATUS,', meta_out)]

  #extra_meta #<- build_extra_meta(ncbi_meta, genomename)

  contigs_vec  <- unique(c(order_contigs(fm$CHROM_map), "chrUnk"))
  contig_lines #<- unique(build_contig_lines(contigs_vec, ncbi_meta))

  add_if_missing <- function(vec, line_glob, line_exact) {
    if (!any(grepl(line_glob, vec))) c(vec, line_exact) else vec
  }
  meta_out <- add_if_missing(meta_out, '^##INFO=<ID=MAPSTATUS,',       '##INFO=<ID=MAPSTATUS,Number=1,Type=String,Description="Unique_mapping|Failed_to_map_uniquely">')
  meta_out <- add_if_missing(meta_out, '^##INFO=<ID=PriorORIENT,',     '##INFO=<ID=PriorORIENT,Number=1,Type=String,Description="Accumulated strand orientation across runs (plus|minus|none)">')


  writeLines(c(meta_out, extra_meta, contig_lines, header_line), con = Outputfilename)

  ct_all_as_char <- setNames(rep("character", length(colnames_vcf)), colnames_vcf)

  mf <- if (requireNamespace("fastmatch", quietly = TRUE)) fastmatch::fmatch else base::match

  # Use a single connection and stream forward
  con <- if (grepl("\\.gz$", raw_in, ignore.case = TRUE)) gzfile(raw_in, "rt") else file(raw_in, "rt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  # Skip header lines already parsed
  if (header_n > 0L) invisible(readLines(con, n = header_n, warn = FALSE))

  # tune chunk size (bigger = faster until memory says stop)
  chunk_size <- getOption("ANCHOR_VCF_CHUNK", 5000L)  # you can export ANCHOR_VCF_CHUNK to tweak on HPC

  # helper: read next N data lines as a data.frame via fread(text=...)
  read_chunk_df <- function(n) {
    lines <- readLines(con, n = n, warn = FALSE)
    if (!length(lines)) return(NULL)
    txt <- paste(lines, collapse = "\n")
    data.table::fread(
      text       = txt,
      sep        = "\t",
      header     = FALSE,
      col.names  = colnames_vcf,
      colClasses = ct_all_as_char,
      data.table = FALSE,
      integer64  = "character",
      showProgress = FALSE
    )
  }

  repeat {
    df <- read_chunk_df(chunk_size)
    if (is.null(df) || !nrow(df)) break
    df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)

    IDv    <- df[[ix_ID]]
    force_drop <- if (length(DROP_IDS)) IDv %in_drop% DROP_IDS else rep(FALSE, length(IDv))
    REF_in <- toupper(df[[ix_REF]])
    ALT_in <- as.character(df[[ix_ALT]])

    # fast (f)match into fm
    j          <- mf(IDv, fm$ID)
    in_map     <- !is.na(j)
    hits       <- ifelse(in_map, fm$MAP_HITS[j], 0L)
    chr_ok     <- in_map & nzchar_safe(fm$CHROM_map[j])
    pos_ok     <- in_map & !is.na(fm$POS_map[j]) & suppressWarnings(as.numeric(fm$POS_map[j]) > 0)
    mapped_now <- in_map & (hits == 1L) & chr_ok & pos_ok & !force_drop

    # Update CHROM/POS/QUAL/FILTER for this iteration
    df[[ix_CHROM]] <- ifelse(mapped_now, fm$CHROM_map[j], "chrUnk")
    df[[ix_POS]]   <- ifelse(mapped_now, as.character(fm$POS_map[j]), "0")
    if (!is.na(ix_QUAL)) df[[ix_QUAL]] <- ifelse(mapped_now, "100", "0")
    if (!is.na(ix_FIL )) df[[ix_FIL ]] <- ifelse(mapped_now, "PASS", "LowQual")

    # Pull prior & current strand
    info_in    <- if (!is.na(ix_INFO)) df[[ix_INFO]] else rep("", nrow(df))
    prior_raw  <- info_get(info_in, "PriorORIENT")
    prior      <- pr_norm(prior_raw)                                  # plus|minus|none
    s_now      <- strand_to_word(ifelse(mapped_now, fm$STRAND[j], NA_character_))  # plus|minus|NA

    # Newly mapped this round? set PriorORIENT = s_now; else keep prior
    prior_out  <- prior
    set_pr     <- !is.na(s_now)
    prior_out[set_pr] <- s_now[set_pr]

    # Brioche alleles
    BriocheREF   <- ifelse(mapped_now, toupper(fm$REF_map[j]), NA_character_)
    BriocheALT   <- ifelse(mapped_now, first_alt_uc(fm$ALT_map[j]), NA_character_)
    BriocheREFRC <- revcomp_iupac(BriocheREF)

    # Flip REF/ALT if prior != none and orientation changed
    need_flip <- mapped_now & prior != "none" & !is.na(s_now) & prior != s_now
    if (any(need_flip)) {
      REF_cur <- df[[ix_REF]]
      ALT_cur <- df[[ix_ALT]]
      REF_cur[need_flip] <- revcomp_iupac(REF_cur[need_flip])
      # reverse-complement CSV of ALTs
      if (any(need_flip)) {
        ALT_cur[need_flip] <- vapply(ALT_cur[need_flip], rc_csv, character(1))
      }
      df[[ix_REF]] <- REF_cur
      df[[ix_ALT]] <- ALT_cur
    }

    # FIRST mapping: prior == 'none' (treat as biallelic initial orientation)
    ALT_curU       <- toupper(as.character(df[[ix_ALT]]))
    old_alt_single <- !grepl(",", ALT_curU, fixed = TRUE) & nzchar_safe(ALT_curU)
    bri_alt_single <- mapped_now & nzchar_safe(BriocheALT) &
                      !grepl(",", ifelse(mapped_now, fm$ALT_map[j], ""), fixed = TRUE)

    eligible_first <- mapped_now & (prior == "none") & old_alt_single & bri_alt_single & nzchar_safe(BriocheREF)

    gen_REF <- toupper(df[[ix_REF]])
    gen_ALT <- first_alt_uc(ALT_curU)

    cond1 <- eligible_first & (BriocheREF   == gen_REF)
    cond2 <- eligible_first & (BriocheREFRC == gen_REF)
    cond3 <- eligible_first & (BriocheREF   == gen_ALT)
    cond4 <- eligible_first & (BriocheREFRC == gen_ALT)

    swap_rows    <- which(cond3 | cond4)
    oriented_now <- cond1 | cond2 | cond3 | cond4

    if (any(cond2 | cond3 | cond4)) {
      idx_bri <- which(cond2 | cond3 | cond4)
      tmpREF <- df[[ix_REF]]; tmpALT <- df[[ix_ALT]]
      tmpREF[idx_bri] <- BriocheREF[idx_bri]
      tmpALT[idx_bri] <- BriocheALT[idx_bri]
      df[[ix_REF]] <- tmpREF; df[[ix_ALT]] <- tmpALT
    }

    # Genotype swap where needed (apply column-wise; vectorized per column)
    if (length(swap_rows) && length(sample_cols_idx)) {
      for (jjj in sample_cols_idx) {
        colv <- df[[jjj]]
        colv[swap_rows] <- swap_biallelic_gt_vec(colv[swap_rows])
        df[[jjj]] <- colv
      }
    }

    # SUBSEQUENT mapping (same orientation) — allow tri+ remapping
    tri_rows <- which(mapped_now & prior != "none" & !is.na(s_now) & (prior == s_now) & nzchar(BriocheREF))
    if (length(tri_rows)) {
      REF_cur2 <- df[[ix_REF]]
      ALT_cur2 <- df[[ix_ALT]]

      # prepare once to avoid repeated object growth
      idx_map_by_row <- vector("list", length(tri_rows))
      names(idx_map_by_row) <- as.character(tri_rows)
      recode_needed  <- logical(length(tri_rows))

      for (k in seq_along(tri_rows)) {
        i <- tri_rows[k]
        roC <- REF_cur2[i]; roU <- toupper(roC)
        altC_vec <- split_alts_csv(as.character(ALT_cur2[i]))
        altU_vec <- toupper(altC_vec)
        rn <- BriocheREF[i]
        if (!nzchar(rn)) next

        if (roU == rn) next

        m_in_old <- match(rn, altU_vec)
        if (!is.na(m_in_old)) {
          new_alt_vecC <- c(roC, altC_vec[-m_in_old])     # swap-to-ref
        } else {
          new_alt_vecC <- c(roC, altC_vec)                # expand
        }
        REF_cur2[i]  <- rn
        ALT_cur2[i]  <- if (length(new_alt_vecC)) paste(new_alt_vecC, collapse = ",") else "."
        idx_map_by_row[[as.character(i)]] <- build_symbol_maps(roU, altU_vec, rn, toupper(new_alt_vecC))
        recode_needed[k] <- TRUE
      }

      df[[ix_REF]] <- REF_cur2
      df[[ix_ALT]] <- ALT_cur2

      if (length(sample_cols_idx) && any(recode_needed)) {
        rows_to_recode <- tri_rows[recode_needed]
        for (jjj in sample_cols_idx) {
          colv <- df[[jjj]]
          for (rr in rows_to_recode) {
            idx_map <- idx_map_by_row[[as.character(rr)]]
            if (!is.null(idx_map)) colv[rr] <- recode_gt_vec(colv[rr], idx_map)
          }
          df[[jjj]] <- colv
        }
      }
    }

    # INFO fields
    status_val <- ifelse(mapped_now, "Unique_mapping", "Failed_to_map_uniquely")

    if (!is.na(ix_INFO)) {
      df[[ix_INFO]] <- info_strip_keys(df[[ix_INFO]], c("MAPSTATUS","PriorORIENT"))
      df[[ix_INFO]] <- upsert_info_key_vec(df[[ix_INFO]], "MAPSTATUS",   status_val)
      df[[ix_INFO]] <- upsert_info_key_vec(df[[ix_INFO]], "PriorORIENT", prior_out)
    }

    data.table::fwrite(df, file = Outputfilename, sep = "\t",
                       append = TRUE, col.names = FALSE, quote = FALSE)

    rm(df); gc(FALSE)
  }
  message(sprintf("Wrote: %s", Outputfilename))
}


# dart file 
# DArT pathway (raw DArT table -> VCF; record PriorORIENT from sstrand)
if (isDArTfile) {

  strand_to_word <- function(x) {
    z <- tolower(trimws(as.character(x)))
    ifelse(is.na(z) | z == "" | z == ".", "none",
      ifelse(z %in% c("+","plus","pos","positive","p"), "plus",
        ifelse(z %in% c("-","minus","neg","negative","m"), "minus", "none")))
  }
  is_blank <- function(x) is.na(x) | x == "" | x == "."

  # ---- existing DArT helpers (detect_dart_layout_base, parse_snp_ref_alt, etc.) remain as in your script ----

  lay <- detect_dart_layout_base(raw_path)
  cn  <- lay$cn
  idx_allele <- match("AlleleID", cn)
  idx_snp    <- match("SNP",      cn)
  if (is.na(idx_allele) || is.na(idx_snp))
    stop("DArT header must include AlleleID and SNP.")

  g_start  <- lay$lastmetric_idx + 1L
  geno_idx <- g_start:length(cn)

  fileDate     <- format(Sys.Date(), "%Y%m%d")
  VCFchroms    <- unique(c(order_contigs(Finalmappings$saccver), "chrUnk"))
  contig_lines <- #unique(build_contig_lines(VCFchroms, ncbi_meta))
  extra_meta   <- #build_extra_meta(ncbi_meta, genomename)

if (nzchar(droplist_path)) {
  header_lines <- c(header_lines,
    sprintf('##Brioche_droplist="%s"', basename(droplist_path)))
}

  header_lines <- c(
    "##fileformat=VCFv4.3",
    sprintf("##fileDate=%s", fileDate),
    "##source=Brioche-VCF-build",
    paste0("##",Metadata[1]),
    paste0("##",Metadata[2]),
    paste0("##",Metadata[3]),
    paste0("##Brioche_filtering_",Metadata[4]),
    paste0("##Brioche_filtering_",Metadata[5]),
    paste0("##Brioche_filtering_",Metadata[6]),
    paste0("##Brioche_priors_files:",Metadata[7]),
    paste0("##Brioche_priors_files:",Metadata[8]),
    paste0("##Brioche_priors_files:",Metadata[9]),
    paste0("##Brioche_priors_files:",Metadata[10]),
    paste0("##Brioche_",Metadata[11]),
    extra_meta,
    "##FILTER=<ID=LowQual,Description=\"Low quality or ambiguous mapping\">",
    '##INFO=<ID=MAPSTATUS,Number=1,Type=String,Description="Unique_mapping|Failed_to_map_uniquely">',
    '##INFO=<ID=PriorORIENT,Number=1,Type=String,Description="Accumulated strand orientation across runs (plus|minus|none)">',
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

  # Mapping frame (include sstrand)
  fm <- Finalmappings %>%
    dplyr::transmute(
      Name      = qaccver,
      CHROM_map = saccver,
      POS_map   = SNPpos,
      REF_map   = toupper(Ref),
      ALT_map   = normalize_alt(ALT),
      STRAND    = sstrand
    )
  fm_hits <- Finalmappings %>% dplyr::count(qaccver, name = "MAP_HITS")
  fm <- fm %>%
    dplyr::left_join(fm_hits, by = c("Name" = "qaccver")) %>%
    dplyr::mutate(MAP_HITS = ifelse(is.na(MAP_HITS), 0L, MAP_HITS))
  fm_Name <- fm$Name

  # Stream
  chunk_size <- 1500L
  delim <- if (grepl("\\.(tsv|txt)$", raw_path, ignore.case = TRUE)) "\t" else ","
  con_dart <- if (grepl("\\.gz$", raw_path, ignore.case = TRUE)) gzfile(raw_path, "rt") else file(raw_path, "rt")
  on.exit(try(close(con_dart), silent = TRUE), add = TRUE)
  invisible(readLines(con_dart, n = lay$topskip + 1L, warn = FALSE))  # skip preamble+header

  repeat {
    raw_lines <- readLines(con_dart, n = chunk_size, warn = FALSE)
    if (!length(raw_lines)) break

    df_all <- utils::read.table(text = raw_lines, sep = delim, header = FALSE,
                                quote = "", comment.char = "", stringsAsFactors = FALSE,
                                check.names = FALSE, col.names = cn_fixed, fill = TRUE)
    df <- df_all[, c(idx_allele, idx_snp, geno_idx), drop = FALSE]

    Name <- df[[1L]]
    force_drop <- if (length(DROP_IDS)) Name %in_drop% DROP_IDS else rep(FALSE, length(Name))
    SNP  <- df[[2L]]
    geno <- as.matrix(df[, -(1:2), drop = FALSE])

    pa   <- parse_snp_ref_alt(SNP)
    REF_old <- toupper(pa$REF)
    ALT_old <- toupper(pa$ALT)  # single ALT

    idx   <- match(Name, fm_Name)
    has_m <- !is.na(idx)
    hits   <- ifelse(has_m, fm$MAP_HITS[idx], 0L)
    u_map  <- has_m & hits == 1L & !force_drop

    CHROM <- ifelse(u_map, fm$CHROM_map[idx], "chrUnk")
    POS   <- ifelse(u_map, fm$POS_map[idx],   0L)

    BriocheREF   <- ifelse(u_map, toupper(fm$REF_map[idx]), NA_character_)
    BriocheALT   <- ifelse(u_map, toupper(sub(",.*$", "", fm$ALT_map[idx])), NA_character_)
    BriocheREFRC <- ifelse(!is.na(BriocheREF), revcomp_iupac(BriocheREF), NA_character_)

    genotypeREF <- REF_old
    genotypeALT <- ALT_old

    ident    <- u_map & !is.na(BriocheREF)   & (BriocheREF   == genotypeREF)
    ident_rc <- u_map & !is.na(BriocheREFRC) & (BriocheREFRC == genotypeREF)
    swap     <- u_map & !is.na(BriocheREF)   & (BriocheREF   == genotypeALT)
    swap_rc  <- u_map & !is.na(BriocheREFRC) & (BriocheREFRC == genotypeALT)

    oriented_now <- ident | ident_rc | swap | swap_rc

    REF_out <- ifelse(oriented_now, BriocheREF, REF_old)
    ALT_out <- ifelse(oriented_now, BriocheALT, ALT_old)
    ALT_out <- sub(",.*$", "", ALT_out); ALT_out[is.na(ALT_out) | ALT_out == ""] <- "."

    pos_int     <- suppressWarnings(as.integer(POS))
    coord_bad   <- is.na(pos_int) | pos_int <= 0L | is_blank(CHROM)
    unknown_now <- !u_map | coord_bad

    QUAL   <- ifelse(unknown_now, 0L, 100L)
    FILTER <- ifelse(unknown_now, "LowQual", "PASS")

    # PriorORIENT from current sstrand when uniquely mapped; else 'none'
    cur_strand <- ifelse(u_map, strand_to_word(fm$STRAND[idx]), "none")
    INFO <- ifelse(
      unknown_now,
      paste0("MAPSTATUS=Failed_to_map_uniquely;PriorORIENT=", cur_strand),
      paste0("MAPSTATUS=Unique_mapping;PriorORIENT=", cur_strand)
    )

    # Genotype normalization & swap where needed
    geno[] <- toupper(gsub("\\|", "/", geno, perl = TRUE))
    miss <- is.na(geno) | geno == "" | geno == "." | geno == "-" | geno == "N" | geno == "NC"
    if (any(miss)) geno[miss] <- "./."
    dbl <- grepl("^[012]{2}$", geno, perl = TRUE)
    if (any(dbl)) geno[dbl] <- sub("^([012])([012])$", "\\1/\\2", geno[dbl], perl = TRUE)

    rsame <- ifelse(is.na(ident | ident_rc), FALSE, (ident | ident_rc))
    rswap <- ifelse(is.na(swap  | swap_rc ), FALSE, (swap  | swap_rc ))
    if (any(rswap)) {
      M_swap <- matrix(rep(rswap, times = ncol(geno)), nrow = nrow(geno), ncol = ncol(geno))
      g <- geno
      d0 <- (g == "0") & M_swap; if (any(d0)) g[d0] <- "1/1"
      d1 <- (g == "1") & M_swap; if (any(d1)) g[d1] <- "0/1"
      d2 <- (g == "2") & M_swap; if (any(d2)) g[d2] <- "0/0"
      geno[M_swap] <- g[M_swap]
    }
    i0 <- geno == "0"; if (any(i0)) geno[i0] <- "0/0"
    i1 <- geno == "1"; if (any(i1)) geno[i1] <- "0/1"
    i2 <- geno == "2"; if (any(i2)) geno[i2] <- "1/1"

    # AC/AN/MAF
    nr <- nrow(geno); ac1 <- integer(nr); an <- integer(nr)
    for (j in seq_len(ncol(geno))) {
      x <- geno[, j]
      msk <- x == "./."
      a1 <- substr(x, 1L, 1L); a1[msk] <- NA_character_
      a2 <- substr(x, 3L, 3L); a2[msk] <- NA_character_
      ac1 <- ac1 + as.integer((!is.na(a1)) & (a1 == "1")) + as.integer((!is.na(a2)) & (a2 == "1"))
      an  <- an  + as.integer(!msk) * 2L
    }
    AC_str  <- as.character(ac1)
    AF1     <- ifelse(an > 0L, ac1 / an, NA_real_)
    MAF_str <- sprintf("%.6g", AF1)
    INFO <- paste0(INFO, ";MAF=", MAF_str, ";AC=", AC_str, ";AN=", an)

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

