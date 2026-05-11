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
#     --Outputfilename my_output.vcf \
#     --Accession "" \
#     --genomename ""\
#     --reference_species "" \
#     --reference_url "" \
#     --droplist "" \
#     --mappriors "marker_unique_map_priors.tsv"
#
# Notes:
#   - --Rawgenotypes : path to the raw genotype table/vcf file. Must have minimum values of Name REF, ALT, and a genotype matrix of data 
#   - --Briochemappings : path to Brioche mappings CSV *Brioche_all_markers1to1stagingforvcf
#   - --IsVCFraws : if TRUE, sample raw genotype file is already in a vcf format
#   - --IsSNPchip : if TRUE, data being input is from an ilumina SNP chip
#   - --isDArTfile false : If TRUE, sample raw genotype file is in DArT SNP format. File must be 1 row format from Dart e.g., 1 allele with genotypes coded as 0,1,2 (as opposed to two row where there are two rows for the same allele with one representing the ref and one the alt and genotypes are coded as 1 0 binary
#   - --Outputfilename : path/filename to write VCF (default: output.vcf)
#   - --Accession : Accession of the reference genome used for mapping in Brioche e.g., GCF_000331145.2 (don't do GCA)  If provided, metadata from NCBI will be incorporated into vcf header info
#   - --genomename :genome name (title) of the reference genome used for mapping in Brioche e.g., Cicar.CDCFrontier_v2.0  If provided, metadata from NCBI will be incorporated into vcf header info. If Accession and genomename are provided, Accession is used first then genomename afterwards if no metadata was found using accession
#   - --reference_species "" :User can input a specific reference species to be added to the header metadata of the vcf file. (Useful for when not working with an NCBI reference)
#   - --reference_url "" :User can input a specific reference url to be added to the header metadata of the vcf file. (Useful for when the genome was downloaded from a repository that isn't NCBI)
#   - --droplist "" :user can provide drop list of markernames which should be changed to unknown Chromosome and position [allow for duplicates to be dealt with]. 
#   - --mappriors : path to a TSV with columns:qaccver, Trueunique, Chrommarkermap, proximatemarkermap,linkagemarkermap, geneticmapmarkermap, Duplicate_region. Used to differentiate priors derrived unique markers from non priors markers.
#       



conda_lib <- Sys.getenv("CONDA_PREFIX")
if (nzchar(conda_lib)) {
  conda_r_lib <- file.path(conda_lib, "lib", "R", "library")
  .libPaths(c(conda_r_lib, .libPaths()))
}



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

if (!requireNamespace("rentrez", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("rentrez", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
if (!requireNamespace("readr", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("readr", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
# data.table intentionally NOT loaded; optional via --use_data_table / BRIOCHE_USE_DT
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("tidyverse", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
if (!requireNamespace("stringi", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("stringi", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
if (!requireNamespace("tidyselect", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("tidyselect", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
if (!requireNamespace("magrittr", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("magrittr", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("tidyr", repos = "https://cloud.r-project.org", quiet = TRUE)))
}


suppressPackageStartupMessages({
  library(getopt)
  library(dplyr)
  library(stringr)
  library(rentrez)
  library(readr)
  # DO NOT: library(data.table)
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

raw_path      <- get_opt(opts, "Rawgenotypes", "")
map_path      <- get_opt(opts, "Briochemappings", "")
is_vcfraws    <- to_bool(get_opt(opts, "IsVCFraws", "false"))
is_snpchip    <- to_bool(get_opt(opts, "IsSNPchip", "true"))
isDArTfile    <- to_bool(get_opt(opts, "isDArTfile", "false"))
Outputfilename<- get_opt(opts, "Outputfilename", "outfile")
droplist_path <- get_opt(opts, "droplist",   "")
Accession     <- get_opt(opts, "Accession",  "")
genomename    <- get_opt(opts, "genomename", "")
user_species   <- get_opt(opts, "reference_species", "")
user_genome_url<- get_opt(opts, "reference_url",   "")
mappriors_path <- get_opt(opts, "mappriors", "")
# Optional toggle to allow data.table on capable nodes
use_dt <- tolower(get_opt(opts, "use_data_table", Sys.getenv("BRIOCHE_USE_DT", "false"))) %in% c("1","true","t","yes","y")

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
message(sprintf("mappriors      : %s", mappriors_path))
message(sprintf("use_data_table : %s", if (use_dt) "true" else "false"))

# ---------- Safe I/O shims (prefer base/readr; optionally data.table if explicitly enabled) ----------
dt_fread <- function(path, sep="\t", header=TRUE) {
  if (use_dt && requireNamespace("data.table", quietly=TRUE)) {
    data.table::fread(path, sep=sep, header=header, data.table=FALSE, showProgress=FALSE)
  } else {
    if (identical(sep, "\t")) readr::read_tsv(path, show_col_types=FALSE, progress=FALSE)
    else if (identical(sep, ",")) readr::read_csv(path, show_col_types=FALSE, progress=FALSE)
    else utils::read.table(path, sep=sep, header=header, quote="", comment.char="", stringsAsFactors=FALSE, check.names=FALSE)
  }
}
dt_fwrite <- function(df, file, append=TRUE, col.names=FALSE) {
  if (use_dt && requireNamespace("data.table", quietly=TRUE)) {
    data.table::fwrite(df, file=file, sep="\t", append=append, col.names=col.names, quote=FALSE)
  } else {
    utils::write.table(df, file=file, sep="\t", append=append, col.names=col.names, row.names=FALSE, quote=FALSE)
  }
}
write_header_row <- function(vec, file) {
  cat(paste(vec, collapse="\t"), "\n", file=file, append=TRUE)
}

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



`%||%` <- function(a, b) {
  if (!is.null(a) &&
      length(a) > 0 &&
      !all(is.na(a)) &&
      any(nzchar(as.character(a)))) {
    a
  } else {
    b
  }
}


# For genome info gathering function to search NCBI
# Testing Update which incorporates 1. GCA in assembly, BioProject codes, and searches assembly/genomes/bioproject deferring to assembly but returning others as partial data if all that is avaliable
fetch_ncbi_meta <- function(accession = "", genomename = "",
                            search_genome_db = TRUE,
                            search_bioproject_db = TRUE,
                            retmax = 10) {
  out <- list(ok = FALSE)

  empty_to_na <- function(x) {
    x <- as.character(x %||% NA_character_)
    x <- trimws(x)

    x_norm <- tolower(x)
    x_norm <- gsub("[[:space:]]+", " ", x_norm)
    x_norm <- gsub("^[[:punct:]]+|[[:punct:]]+$", "", x_norm)

    missing_tokens <- c(
      "",
      "na",
      "n/a",
      "n.a",
      "nan",
      "null",
      "none",
      "nil",
      "missing",
      "unknown",
      "unk",
      "not known",
      "not applicable",
      "not available",
      "unavailable",
      "unspecified",
      "undefined",
      "undetermined",
      "unassigned",
      "not assigned",
      "no data",
      "no information",
      "no value"
    )

    x[x_norm %in% missing_tokens] <- NA_character_
    x
  }

  first_nonempty <- function(...) {
    vals <- list(...)
    for (v in vals) {
      v <- empty_to_na(v)
      if (length(v) > 0 && any(!is.na(v) & nzchar(v))) {
        return(v[which(!is.na(v) & nzchar(v))[1]])
      }
    }
    NA_character_
  }

  nz <- function(x) {
    x <- empty_to_na(x)
    length(x) == 1 && !is.na(x) && nzchar(x)
  }

  any_nz <- function(x) {
    x <- empty_to_na(x)
    any(!is.na(x) & nzchar(x))
  }

  strip_asm_version <- function(x) {
    sub("\\.[0-9]+$", "", toupper(trimws(x)))
  }

  is_assembly_accession <- function(x) {
    grepl("^(GCA|GCF)_", trimws(x), ignore.case = TRUE)
  }

  is_bioproject_accession <- function(x) {
    grepl("^PRJ[A-Z]{2,}[0-9]+$", trimws(x), ignore.case = TRUE)
  }

  esc_query <- function(x) {
    gsub('"', '\\"', x)
  }

  normalize_summaries <- function(sums) {
    if (is.null(sums)) {
      list()
    } else if (is.list(sums) && !is.null(sums$uid)) {
      list(sums)
    } else if (is.list(sums)) {
      sums
    } else {
      list(sums)
    }
  }

  safe_search <- function(db, term, retmax = 10) {
    tryCatch(
      rentrez::entrez_search(db = db, term = term, retmax = retmax),
      error = function(e) list(ids = character(0), count = 0)
    )
  }

  safe_summary <- function(db, ids) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]

    if (!length(ids)) return(list())

    tryCatch(
      normalize_summaries(rentrez::entrez_summary(db = db, id = ids)),
      error = function(e) list()
    )
  }

  get_asm_acc <- function(x) {
    empty_to_na(x[["assemblyaccession"]] %||% x[["accn"]] %||% NA_character_)
  }

  get_asm_name <- function(x) {
    empty_to_na(x[["assemblyname"]] %||% NA_character_)
  }

  get_ftp <- function(x, prefer = c("auto", "refseq", "genbank")) {
    prefer <- match.arg(prefer)

    ftp_ref <- empty_to_na(x[["ftppath_refseq"]] %||% NA_character_)
    ftp_gb  <- empty_to_na(x[["ftppath_genbank"]] %||% NA_character_)

    if (prefer == "refseq") {
      if (!is.na(ftp_ref)) ftp_ref else ftp_gb
    } else if (prefer == "genbank") {
      if (!is.na(ftp_gb)) ftp_gb else ftp_ref
    } else {
      acc <- get_asm_acc(x)

      if (!is.na(acc) && grepl("^GCA_", acc, ignore.case = TRUE)) {
        if (!is.na(ftp_gb)) ftp_gb else ftp_ref
      } else {
        if (!is.na(ftp_ref)) ftp_ref else ftp_gb
      }
    }
  }

  unique_or_null <- function(ids) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]

    if (length(ids) == 1) ids else NULL
  }

  choose_assembly_by_accession <- function(ids, acc_in) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]

    if (!length(ids)) return(NULL)

    sums <- safe_summary("assembly", ids)
    if (!length(sums)) return(NULL)

    acc_in_full  <- toupper(trimws(acc_in))
    acc_in_short <- strip_asm_version(acc_in)

    asm_accs <- vapply(sums, get_asm_acc, character(1))
    asm_accs <- empty_to_na(asm_accs)

    full_match <- !is.na(asm_accs) & toupper(asm_accs) == acc_in_full
    if (sum(full_match) == 1) {
      return(ids[which(full_match)])
    }

    short_match <- !is.na(asm_accs) & strip_asm_version(asm_accs) == acc_in_short
    if (sum(short_match) == 1) {
      return(ids[which(short_match)])
    }

    NULL
  }

  choose_assembly_by_name <- function(ids, name_in) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]

    if (!length(ids)) return(NULL)

    sums <- safe_summary("assembly", ids)
    if (!length(sums)) return(NULL)

    asm_names <- vapply(sums, get_asm_name, character(1))
    asm_names <- empty_to_na(asm_names)

    exact <- !is.na(asm_names) & tolower(asm_names) == tolower(trimws(name_in))

    if (sum(exact) == 1) {
      return(ids[which(exact)])
    }

    NULL
  }

  link_database_to_assembly <- function(ids, dbfrom) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]

    if (!length(ids)) return(character(0))

    linked <- character(0)

    for (id_i in ids) {
      lk <- tryCatch(
        rentrez::entrez_link(dbfrom = dbfrom, db = "assembly", id = id_i),
        error = function(e) NULL
      )

      if (!is.null(lk) && !is.null(lk$links)) {
        linked <- unique(c(
          linked,
          unlist(lk$links, use.names = FALSE)
        ))
      }
    }

    linked <- unique(as.character(linked))
    linked[!is.na(linked) & nzchar(linked)]
  }

  link_genome_to_assembly <- function(genome_ids) {
    link_database_to_assembly(genome_ids, dbfrom = "genome")
  }

  link_bioproject_to_assembly <- function(bioproject_ids) {
    link_database_to_assembly(bioproject_ids, dbfrom = "bioproject")
  }

  search_genome_ids <- function(query_text, mode = c("accession", "name")) {
    mode <- match.arg(mode)

    query_text <- trimws(query_text)
    if (!nzchar(query_text)) return(character(0))

    if (mode == "accession") {
      g_queries <- unique(c(
        sprintf('"%s"[All Fields]', esc_query(query_text)),
        sprintf('%s[All Fields]', query_text),
        if (grepl("\\.[0-9]+$", query_text)) {
          sprintf('"%s"[All Fields]', strip_asm_version(query_text))
        } else {
          character(0)
        }
      ))
    } else {
      g_queries <- unique(c(
        sprintf('"%s"[All Fields]', esc_query(query_text)),
        sprintf('"%s"[Title]', esc_query(query_text))
      ))
    }

    genome_ids <- character(0)

    for (q in g_queries) {
      gs <- safe_search("genome", q, retmax = retmax)
      genome_ids <- unique(c(genome_ids, gs$ids))
    }

    genome_ids <- unique(as.character(genome_ids))
    genome_ids[!is.na(genome_ids) & nzchar(genome_ids)]
  }

  search_bioproject_ids <- function(query_text, mode = c("accession", "name")) {
    mode <- match.arg(mode)

    query_text <- trimws(query_text)
    if (!nzchar(query_text)) return(character(0))

    if (mode == "accession") {
      bp_queries <- unique(c(
        sprintf('"%s"[All Fields]', esc_query(query_text)),
        sprintf('%s[All Fields]', query_text),
        sprintf('"%s"[Project Accession]', esc_query(query_text)),
        sprintf('"%s"[BioProject Accession]', esc_query(query_text))
      ))
    } else {
      bp_queries <- unique(c(
        sprintf('"%s"[All Fields]', esc_query(query_text)),
        sprintf('"%s"[Title]', esc_query(query_text)),
        sprintf('"%s"[Project Name]', esc_query(query_text))
      ))
    }

    bioproject_ids <- character(0)

    for (q in bp_queries) {
      bp <- safe_search("bioproject", q, retmax = retmax)
      bioproject_ids <- unique(c(bioproject_ids, bp$ids))
    }

    bioproject_ids <- unique(as.character(bioproject_ids))
    bioproject_ids[!is.na(bioproject_ids) & nzchar(bioproject_ids)]
  }

  search_genome_for_assemblies <- function(query_text, mode = c("accession", "name")) {
    genome_ids <- search_genome_ids(query_text, mode = match.arg(mode))
    link_genome_to_assembly(genome_ids)
  }

  search_bioproject_for_assemblies <- function(query_text, mode = c("accession", "name")) {
    bioproject_ids <- search_bioproject_ids(query_text, mode = match.arg(mode))
    link_bioproject_to_assembly(bioproject_ids)
  }

  extract_genome_partial <- function(genome_id, query_text = "") {
    sums <- safe_summary("genome", genome_id)
    if (!length(sums)) return(NULL)

    s <- sums[[1]]

    title <- first_nonempty(
      s[["title"]],
      s[["name"]],
      s[["genome"]],
      s[["description"]],
      query_text
    )

    organism <- first_nonempty(
      s[["organism"]],
      s[["organismname"]],
      s[["speciesname"]],
      s[["scientificname"]],
      s[["taxname"]]
    )

    taxid <- first_nonempty(
      s[["taxid"]],
      s[["taxonomyid"]]
    )

    acc <- first_nonempty(
      s[["accession"]],
      s[["assemblyaccession"]],
      s[["accn"]],
      query_text
    )

    url <- paste0("https://www.ncbi.nlm.nih.gov/genome/", genome_id)

    has_any <- any_nz(c(title, organism, taxid, acc))
    if (!has_any) return(NULL)

    list(
      ok                = TRUE,
      partial           = TRUE,
      assembly_resolved = FALSE,
      source_db         = "genome",
      source_uid        = as.character(genome_id),
      source_accession  = acc,
      accession         = acc,
      asm_name          = title,
      organism          = organism,
      taxid             = suppressWarnings(as.integer(taxid)),
      asm_url           = url,
      genome_url        = url,
      ftp_path          = NA_character_,
      ncbi_uid          = as.character(genome_id),
      lookup_used       = "partial",
      lookup_via        = "genome",
      partial_reason    = "Genome record found but no unique linked Assembly record was resolved",
      chr_table         = data.frame(),
      len_by_seqname    = setNames(numeric(0), character(0)),
      md5_by_seqname    = setNames(character(0), character(0))
    )
  }

  extract_bioproject_partial <- function(bioproject_id, query_text = "") {
    sums <- safe_summary("bioproject", bioproject_id)
    if (!length(sums)) return(NULL)

    s <- sums[[1]]

    bp_acc <- first_nonempty(
      s[["project_acc"]],
      s[["projectacc"]],
      s[["projectaccession"]],
      s[["bioprojectaccn"]],
      s[["accession"]],
      query_text
    )

    title <- first_nonempty(
      s[["project_name"]],
      s[["projectname"]],
      s[["title"]],
      s[["name"]],
      s[["description"]],
      query_text
    )

    organism <- first_nonempty(
      s[["organism"]],
      s[["organismname"]],
      s[["speciesname"]],
      s[["scientificname"]],
      s[["taxname"]]
    )

    taxid <- first_nonempty(
      s[["taxid"]],
      s[["taxonomyid"]]
    )

    url <- if (!is.na(bp_acc) && nzchar(bp_acc)) {
      paste0("https://www.ncbi.nlm.nih.gov/bioproject/", bp_acc)
    } else {
      paste0("https://www.ncbi.nlm.nih.gov/bioproject/", bioproject_id)
    }

    has_any <- any_nz(c(bp_acc, title, organism, taxid))
    if (!has_any) return(NULL)

    list(
      ok                = TRUE,
      partial           = TRUE,
      assembly_resolved = FALSE,
      source_db         = "bioproject",
      source_uid        = as.character(bioproject_id),
      source_accession  = bp_acc,
      accession         = bp_acc,
      asm_name          = title,
      organism          = organism,
      taxid             = suppressWarnings(as.integer(taxid)),
      asm_url           = url,
      genome_url        = url,
      ftp_path          = NA_character_,
      ncbi_uid          = as.character(bioproject_id),
      lookup_used       = "partial",
      lookup_via        = "bioproject",
      partial_reason    = "BioProject record found but no unique linked Assembly record was resolved",
      chr_table         = data.frame(),
      len_by_seqname    = setNames(numeric(0), character(0)),
      md5_by_seqname    = setNames(character(0), character(0))
    )
  }

  extract_assembly_summary_partial <- function(asm_sum, asm_uid, reason = "") {
    acc <- first_nonempty(
      asm_sum[["assemblyaccession"]],
      asm_sum[["accn"]]
    )

    title <- first_nonempty(
      asm_sum[["assemblyname"]],
      asm_sum[["title"]],
      asm_sum[["name"]]
    )

    organism <- first_nonempty(
      asm_sum[["organism"]],
      asm_sum[["speciesname"]],
      asm_sum[["scientificname"]],
      asm_sum[["taxname"]]
    )

    taxid <- first_nonempty(
      asm_sum[["taxid"]],
      asm_sum[["taxonomyid"]]
    )

    url <- if (!is.na(acc) && nzchar(acc)) {
      paste0("https://www.ncbi.nlm.nih.gov/datasets/genome/", acc)
    } else {
      paste0("https://www.ncbi.nlm.nih.gov/assembly/", asm_uid)
    }

    has_any <- any_nz(c(acc, title, organism, taxid))
    if (!has_any) return(NULL)

    list(
      ok                = TRUE,
      partial           = TRUE,
      assembly_resolved = FALSE,
      source_db         = "assembly",
      source_uid        = as.character(asm_uid),
      source_accession  = acc,
      accession         = acc,
      asm_name          = title,
      organism          = organism,
      taxid             = suppressWarnings(as.integer(taxid)),
      asm_url           = url,
      genome_url        = url,
      ftp_path          = NA_character_,
      ncbi_uid          = as.character(asm_uid),
      lookup_used       = "partial",
      lookup_via        = "assembly",
      partial_reason    = reason,
      chr_table         = data.frame(),
      len_by_seqname    = setNames(numeric(0), character(0)),
      md5_by_seqname    = setNames(character(0), character(0))
    )
  }

  find_partial_meta <- function(acc_in, name_in) {
    if (nzchar(acc_in) && isTRUE(search_bioproject_db) && is_bioproject_accession(acc_in)) {
      bp_ids <- search_bioproject_ids(acc_in, mode = "accession")
      bp_id <- unique_or_null(bp_ids)

      if (!is.null(bp_id)) {
        partial <- extract_bioproject_partial(bp_id, query_text = acc_in)
        if (!is.null(partial)) return(partial)
      }
    }

    if (nzchar(acc_in) && isTRUE(search_genome_db)) {
      genome_ids <- search_genome_ids(acc_in, mode = "accession")
      genome_id <- unique_or_null(genome_ids)

      if (!is.null(genome_id)) {
        partial <- extract_genome_partial(genome_id, query_text = acc_in)
        if (!is.null(partial)) return(partial)
      }
    }

    if (nzchar(acc_in) && isTRUE(search_bioproject_db)) {
      bp_ids <- search_bioproject_ids(acc_in, mode = "accession")
      bp_id <- unique_or_null(bp_ids)

      if (!is.null(bp_id)) {
        partial <- extract_bioproject_partial(bp_id, query_text = acc_in)
        if (!is.null(partial)) return(partial)
      }
    }

    if (nzchar(name_in) && isTRUE(search_genome_db)) {
      genome_ids <- search_genome_ids(name_in, mode = "name")
      genome_id <- unique_or_null(genome_ids)

      if (!is.null(genome_id)) {
        partial <- extract_genome_partial(genome_id, query_text = name_in)
        if (!is.null(partial)) return(partial)
      }
    }

    if (nzchar(name_in) && isTRUE(search_bioproject_db)) {
      bp_ids <- search_bioproject_ids(name_in, mode = "name")
      bp_id <- unique_or_null(bp_ids)

      if (!is.null(bp_id)) {
        partial <- extract_bioproject_partial(bp_id, query_text = name_in)
        if (!is.null(partial)) return(partial)
      }
    }

    NULL
  }

  find_assembly_id <- function(acc_in, name_in) {
    if (nzchar(acc_in)) {
      acc_in <- trimws(acc_in)

      queries <- unique(c(
        sprintf('%s[Assembly Accession]', acc_in),
        sprintf('"%s"[Assembly Accession]', esc_query(acc_in)),
        sprintf('"%s"[All Fields]', esc_query(acc_in)),
        if (grepl("\\.[0-9]+$", acc_in)) {
          sprintf('"%s"[All Fields]', strip_asm_version(acc_in))
        } else {
          character(0)
        }
      ))

      ids <- character(0)

      for (q in queries) {
        s <- safe_search("assembly", q, retmax = retmax)
        ids <- unique(c(ids, s$ids))
      }

      if (is_assembly_accession(acc_in)) {
        pick_id <- choose_assembly_by_accession(ids, acc_in)
      } else {
        pick_id <- unique_or_null(ids)
      }

      if (!is.null(pick_id)) {
        return(list(id = pick_id, used = "acc", via = "assembly"))
      }

      if (isTRUE(search_genome_db)) {
        asm_ids <- search_genome_for_assemblies(acc_in, mode = "accession")

        if (is_assembly_accession(acc_in)) {
          pick_id <- choose_assembly_by_accession(asm_ids, acc_in)
        } else {
          pick_id <- unique_or_null(asm_ids)
        }

        if (!is.null(pick_id)) {
          return(list(id = pick_id, used = "acc", via = "genome_link"))
        }
      }

      if (isTRUE(search_bioproject_db)) {
        asm_ids <- search_bioproject_for_assemblies(acc_in, mode = "accession")

        if (is_assembly_accession(acc_in)) {
          pick_id <- choose_assembly_by_accession(asm_ids, acc_in)
        } else {
          pick_id <- unique_or_null(asm_ids)
        }

        if (!is.null(pick_id)) {
          return(list(id = pick_id, used = "acc", via = "bioproject_link"))
        }
      }
    }

    if (nzchar(name_in)) {
      name_in <- trimws(name_in)

      queries <- unique(c(
        sprintf('"%s"[Assembly Name]', esc_query(name_in)),
        sprintf('"%s"[All Fields]', esc_query(name_in))
      ))

      ids <- character(0)

      for (q in queries) {
        s <- safe_search("assembly", q, retmax = retmax)
        ids <- unique(c(ids, s$ids))
      }

      pick_id <- choose_assembly_by_name(ids, name_in)

      if (!is.null(pick_id)) {
        return(list(id = pick_id, used = "name", via = "assembly"))
      }

      if (isTRUE(search_genome_db)) {
        asm_ids <- search_genome_for_assemblies(name_in, mode = "name")
        pick_id <- choose_assembly_by_name(asm_ids, name_in)

        if (is.null(pick_id)) {
          pick_id <- unique_or_null(asm_ids)
        }

        if (!is.null(pick_id)) {
          return(list(id = pick_id, used = "name", via = "genome_link"))
        }
      }

      if (isTRUE(search_bioproject_db)) {
        asm_ids <- search_bioproject_for_assemblies(name_in, mode = "name")
        pick_id <- choose_assembly_by_name(asm_ids, name_in)

        if (is.null(pick_id)) {
          pick_id <- unique_or_null(asm_ids)
        }

        if (!is.null(pick_id)) {
          return(list(id = pick_id, used = "name", via = "bioproject_link"))
        }
      }
    }

    list(id = NULL, used = "", via = "")
  }

  tryCatch({
    acc_in  <- trimws(accession %||% "")
    name_in <- trimws(genomename %||% "")

    found <- find_assembly_id(acc_in, name_in)

    if (is.null(found$id)) {
      partial <- find_partial_meta(acc_in, name_in)

      if (!is.null(partial) && isTRUE(partial$ok)) {
        message("No unique Assembly record resolved; returning partial metadata from ",
                partial$source_db)
        return(partial)
      }

      message("Accession or genome name were nonspecific unable to populate genomic information")
      return(out)
    }

    sum <- rentrez::entrez_summary(db = "assembly", id = found$id)

    if (found$used == "acc") {
      acc_sum <- get_asm_acc(sum)

      if (is_assembly_accession(acc_in)) {
        if (is.na(acc_sum) ||
            strip_asm_version(acc_sum) != strip_asm_version(acc_in)) {
          message("Accession or genome name were nonspecific unable to populate genomic information")
          return(out)
        }
      }
    } else if (found$used == "name") {
      aname <- get_asm_name(sum)

      if (found$via == "assembly") {
        if (is.na(aname) || tolower(aname) != tolower(name_in)) {
          message("Accession or genome name were nonspecific unable to populate genomic information")
          return(out)
        }
      }
    }

    acc_sum <- get_asm_acc(sum)
    aname   <- get_asm_name(sum)
    org     <- empty_to_na(sum[["organism"]] %||% sum[["speciesname"]] %||% NA_character_)
    taxid   <- sum[["taxid"]] %||% NA

    ftp <- get_ftp(sum, prefer = "auto")

    if (is.na(acc_sum) || is.na(ftp)) {
      partial <- extract_assembly_summary_partial(
        sum,
        found$id,
        reason = "Assembly summary record found but FTP path or assembly report was unavailable"
      )

      if (!is.null(partial) && isTRUE(partial$ok)) {
        message("Assembly summary resolved but assembly report/FTP path was unavailable; returning partial metadata from assembly")
        return(partial)
      }

      partial <- find_partial_meta(acc_in, name_in)

      if (!is.null(partial) && isTRUE(partial$ok)) {
        message("Assembly record resolved but assembly report/FTP path was unavailable; returning partial metadata from ",
                partial$source_db)
        return(partial)
      }

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

      if (is.na(hdr_i)) {
        stop("assembly_report header line not found")
      }

      header <- sub("^#\\s*", "", lines[hdr_i])
      data   <- lines[(hdr_i + 1L):length(lines)]
      data   <- data[!grepl("^#", data)]

      txt <- paste(c(header, data), collapse = "\n")

      readr::read_tsv(
        I(txt),
        show_col_types = FALSE,
        progress = FALSE,
        col_types = readr::cols(.default = "c")
      )
    }

    ar <- parse_assembly_report(rpt_url)

    col_len      <- intersect(c("Sequence-Length", "sequence-length", "length"), names(ar))
    col_assigned <- intersect(c("Assigned-Molecule", "assigned-molecule"), names(ar))
    col_md5      <- intersect(c("MD5", "md5"), names(ar))
    col_seqname  <- intersect(c("Sequence-Name", "sequence-name"), names(ar))
    col_role     <- intersect(c("Sequence-Role", "sequence-role"), names(ar))
    col_refseq   <- intersect(c("RefSeq-Accn", "refseq-accn", "RefSeq Accn"), names(ar))
    col_genbank  <- intersect(c("GenBank-Accn", "genbank-accn", "GenBank Accn"), names(ar))

    if (!length(col_len)) {
      stop("assembly_report lacks Sequence-Length column")
    }

    key_col <- if (length(col_seqname)) {
      col_seqname[1]
    } else if (length(col_assigned)) {
      col_assigned[1]
    } else {
      NULL
    }

    if (is.null(key_col)) {
      stop("assembly_report lacks Sequence-Name / Assigned-Molecule columns")
    }

    seq_name <- empty_to_na(ar[[key_col]])

    assigned <- if (length(col_assigned)) {
      empty_to_na(ar[[col_assigned[1]]])
    } else {
      rep(NA_character_, nrow(ar))
    }

    len_val <- suppressWarnings(as.numeric(ar[[col_len[1]]]))

    md5_val <- if (length(col_md5)) {
      empty_to_na(ar[[col_md5[1]]])
    } else {
      rep(NA_character_, nrow(ar))
    }

    role_val <- if (length(col_role)) {
      empty_to_na(ar[[col_role[1]]])
    } else {
      rep(NA_character_, nrow(ar))
    }

    ref_acc <- if (length(col_refseq)) {
      empty_to_na(ar[[col_refseq[1]]])
    } else {
      rep(NA_character_, nrow(ar))
    }

    gb_acc <- if (length(col_genbank)) {
      empty_to_na(ar[[col_genbank[1]]])
    } else {
      rep(NA_character_, nrow(ar))
    }

    extract_num <- function(x) {
      s <- tolower(trimws(x))

      m <- stringr::str_match(
        s,
        "^.*(?:chromosome|chrom|chr|contig|scaffold|segment|seg)?[ _-]*0*([0-9]+)([a-z])?(?![a-z0-9]).*"
      )

      suppressWarnings(as.integer(m[, 2]))
    }

    chrom_num <- extract_num(seq_name)

    chrom_num <- ifelse(
      is.na(chrom_num),
      extract_num(assigned),
      chrom_num
    )

    preferred <- if (grepl("^GCA_", acc_sum, ignore.case = TRUE)) {
      ifelse(!is.na(gb_acc) & nzchar(gb_acc), gb_acc, ref_acc)
    } else {
      ifelse(!is.na(ref_acc) & nzchar(ref_acc), ref_acc, gb_acc)
    }

    nuccore <- ifelse(
      !is.na(preferred) & nzchar(preferred),
      paste0("https://www.ncbi.nlm.nih.gov/nuccore/", preferred),
      NA_character_
    )

    keep <- !is.na(len_val) & !is.na(seq_name) & nzchar(seq_name)

    chr_table <- data.frame(
      chrom_num         = chrom_num,
      seq_name          = seq_name,
      assigned_molecule = assigned,
      sequence_role     = role_val,
      refseq_accn       = ref_acc,
      genbank_accn      = gb_acc,
      preferred_accn    = preferred,
      length            = len_val,
      md5               = md5_val,
      nuccore_url       = nuccore,
      stringsAsFactors  = FALSE
    )[keep, , drop = FALSE]

    key_seq <- tolower(chr_table$seq_name)

    len_by_seqname <- stats::setNames(chr_table$length, key_seq)
    md5_by_seqname <- stats::setNames(chr_table$md5, key_seq)

    list(
      ok                = TRUE,
      partial           = FALSE,
      assembly_resolved = TRUE,
      accession         = acc_sum,
      asm_name          = aname %||% "",
      organism          = org %||% "",
      taxid             = suppressWarnings(as.integer(taxid)),
      asm_url           = asm_url,
      genome_url        = asm_url,
      ftp_path          = ftp_https,
      ncbi_uid          = found$id,
      lookup_used       = found$used,
      lookup_via        = found$via,
      source_db         = "assembly",
      source_uid        = found$id,
      source_accession  = acc_sum,
      partial_reason    = NA_character_,
      chr_table         = chr_table,
      len_by_seqname    = len_by_seqname,
      md5_by_seqname    = md5_by_seqname
    )
  }, error = function(e) {
    partial <- find_partial_meta(
      trimws(accession %||% ""),
      trimws(genomename %||% "")
    )

    if (!is.null(partial) && isTRUE(partial$ok)) {
      message("Assembly metadata extraction failed; returning partial metadata from ",
              partial$source_db)
      return(partial)
    }

    message("Accession or genome name were nonspecific unable to populate genomic information")
    out
  })
}


build_extra_meta <- function(ncbi_meta, genomename,
                             override_species = NULL,
                             override_url = NULL) {
  esc <- function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(trimws(x))] <- "unknown"
    gsub('"', '\\\\"', x, fixed = TRUE)
  }

  clean <- function(x, missing = "unknown") {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return(missing)
    }

    x <- as.character(x[1])
    x <- trimws(x)

    if (is.na(x) || !nzchar(x)) {
      missing
    } else {
      x
    }
  }

  clean_taxid <- function(x) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return("unknown")
    }

    tx <- suppressWarnings(as.integer(x[1]))

    if (is.na(tx)) {
      "unknown"
    } else {
      as.character(tx)
    }
  }

  clean_url <- function(...) {
    vals <- list(...)

    for (v in vals) {
      u <- clean(v, missing = "")
      if (nzchar(u)) {
        return(u)
      }
    }

    "unknown"
  }

  gn <- clean(genomename, missing = "unknown")

  is_ok <- isTRUE(ncbi_meta$ok)

  is_partial <- is_ok && (
    isTRUE(ncbi_meta$partial) ||
      identical(ncbi_meta$assembly_resolved, FALSE)
  )

  species_fallback <- function() {
    sp <- clean(override_species, missing = "")

    if (!nzchar(sp)) {
      sp <- clean(gn, missing = "")
    }

    if (!nzchar(sp)) {
      sp <- "unknown"
    }

    sp
  }

  if (is_ok) {
    sp <- clean(ncbi_meta$organism, missing = "")

    if (!nzchar(sp)) {
      sp <- species_fallback()
    }

    acc <- clean(ncbi_meta$accession, missing = "")

    if (!nzchar(acc)) {
      acc <- clean(ncbi_meta$source_accession, missing = "")
    }

    if (!nzchar(acc)) {
      acc <- gn
    }

    if (!nzchar(acc) || is.na(acc)) {
      acc <- "unknown"
    }

    asm_name <- clean(ncbi_meta$asm_name, missing = "")

    if (!nzchar(asm_name)) {
      asm_name <- gn
    }

    if (!nzchar(asm_name) || is.na(asm_name)) {
      asm_name <- "unknown"
    }

    taxid <- clean_taxid(ncbi_meta$taxid)

    meta_url <- clean_url(
      ncbi_meta$asm_url,
      ncbi_meta$genome_url,
      override_url
    )

    source_db <- clean(ncbi_meta$source_db, missing = "")

    if (!nzchar(source_db)) {
      source_db <- clean(ncbi_meta$lookup_via, missing = "")
    }

    if (!nzchar(source_db)) {
      source_db <- "unknown"
    }

    automated_lookup <- if (is_partial) {
      "partial metadata found only"
    } else {
      "resolved by automated process"
    }

    c(
      sprintf("##reference=%s", acc),
      sprintf(
        '##assembly=<accession=%s,name="%s",species="%s",taxonomy=%s,url="%s",source="%s",automated_lookup="%s">',
        acc,
        esc(asm_name),
        esc(sp),
        taxid,
        esc(meta_url),
        esc(source_db),
        esc(automated_lookup)
      ),
      sprintf("##genome_url=%s", meta_url)
    )
  } else {
    sp <- species_fallback()

    url_out <- clean_url(override_url)

    ref_val <- if (!is.na(gn) && nzchar(gn) && gn != "unknown") {
      gn
    } else if (!is.na(sp) && nzchar(sp) && sp != "unknown") {
      sp
    } else {
      "unknown"
    }

    c(
      sprintf("##reference=%s", ref_val),
      sprintf(
        '##assembly=<accession=%s,name="%s",species="%s",taxonomy=%s,url="%s",source="%s",automated_lookup="%s">',
        "unknown",
        esc(gn),
        esc(sp),
        "unknown",
        esc(url_out),
        "unknown",
        "not resolved by automated process"
      ),
      sprintf("##genome_url=%s", url_out)
    )
  }
}


build_contig_lines <- function(ids, meta, override_species = NULL) {
  clean <- function(x, missing = "unknown") {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return(missing)
    }

    x <- as.character(x[1])
    x <- trimws(x)

    if (is.na(x) || !nzchar(x)) {
      missing
    } else {
      x
    }
  }

  clean_num <- function(x) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return(NA_real_)
    }

    y <- suppressWarnings(as.numeric(x[1]))

    if (is.na(y)) {
      NA_real_
    } else {
      y
    }
  }

  clean_taxid <- function(x) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return("unknown")
    }

    tx <- suppressWarnings(as.integer(x[1]))

    if (is.na(tx)) {
      "unknown"
    } else {
      as.character(tx)
    }
  }

  has_value <- function(x) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return(FALSE)
    }

    x <- as.character(x[1])
    !is.na(x) && nzchar(trimws(x))
  }

  esc <- function(x) {
    x <- clean(x, missing = "unknown")
    gsub('"', '\\\\"', x, fixed = TRUE)
  }

  strip_seq_version <- function(x) {
    sub("\\.[0-9]+$", "", trimws(as.character(x)))
  }

  make_keys <- function(x) {
    x <- trimws(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]

    unique(c(
      tolower(x),
      tolower(strip_seq_version(x))
    ))
  }

  ids <- unique(trimws(as.character(ids)))
  ids <- ids[!is.na(ids) & nzchar(ids)]

  if (!length(ids)) {
    return(character(0))
  }

  contig_species <- function() {
    sp <- clean(meta$organism, missing = "")

    if (!nzchar(sp)) {
      sp <- clean(override_species, missing = "")
    }

    if (!nzchar(sp)) {
      sp <- "unknown"
    }

    sp
  }

  tax <- clean_taxid(meta$taxid)

  mk_attr <- function(id,
                      len = NA_real_,
                      acc = NA_character_,
                      alias = NA_character_,
                      tax = "unknown") {
    id    <- clean(id, missing = "unknown")
    acc   <- clean(acc, missing = "unknown")
    alias <- clean(alias, missing = "unknown")

    parts <- c(sprintf("ID=%s", id))

    # Important:
    # Do NOT write length=unknown.
    # VCF/BCF parsers expect contig length, if present, to be numeric.
    if (!is.na(len)) {
      parts <- c(parts, sprintf("length=%d", as.integer(len)))
    }

    parts <- c(parts, sprintf('species="%s"', esc(contig_species())))
    parts <- c(parts, sprintf("taxonomy=%s", clean(tax, missing = "unknown")))
    parts <- c(parts, sprintf("Accession=%s", acc))
    parts <- c(parts, sprintf('Alias="%s"', esc(alias)))

    sprintf("##contig=<%s>", paste(parts, collapse = ","))
  }

  has_resolved_assembly <- isTRUE(meta$ok) &&
    !isTRUE(meta$partial) &&
    !identical(meta$assembly_resolved, FALSE)

  has_chr_table <- has_resolved_assembly &&
    !is.null(meta$chr_table) &&
    NROW(meta$chr_table) > 0

  if (!has_chr_table) {
    return(vapply(
      ids,
      function(id) {
        mk_attr(
          id    = id,
          len   = NA_real_,
          acc   = NA_character_,
          alias = NA_character_,
          tax   = tax
        )
      },
      FUN.VALUE = character(1)
    ))
  }

  ct <- meta$chr_table

  get_col <- function(df, col, default = NA_character_) {
    if (col %in% names(df)) {
      df[[col]]
    } else {
      rep(default, nrow(df))
    }
  }

  seq_col       <- get_col(ct, "seq_name")
  assn_col      <- get_col(ct, "assigned_molecule")
  len_col       <- suppressWarnings(as.numeric(get_col(ct, "length", NA_real_)))
  preferred_col <- get_col(ct, "preferred_accn")
  genbank_col   <- get_col(ct, "genbank_accn")
  refseq_col    <- get_col(ct, "refseq_accn")

  alias_col <- ifelse(
    !is.na(assn_col) & nzchar(trimws(as.character(assn_col))),
    as.character(assn_col),
    as.character(seq_col)
  )

  accession_col <- ifelse(
    !is.na(preferred_col) & nzchar(trimws(as.character(preferred_col))),
    as.character(preferred_col),
    ifelse(
      !is.na(genbank_col) & nzchar(trimws(as.character(genbank_col))),
      as.character(genbank_col),
      as.character(refseq_col)
    )
  )

  lookup <- new.env(parent = emptyenv())

  add_lookup <- function(keys, row_i) {
    keys <- make_keys(keys)

    for (k in keys) {
      if (!exists(k, envir = lookup, inherits = FALSE)) {
        assign(k, row_i, envir = lookup)
      }
    }

    invisible(NULL)
  }

  for (i in seq_len(nrow(ct))) {
    add_lookup(seq_col[i], i)
    add_lookup(assn_col[i], i)
    add_lookup(preferred_col[i], i)
    add_lookup(genbank_col[i], i)
    add_lookup(refseq_col[i], i)
  }

  vapply(
    ids,
    function(id) {
      id_keys <- make_keys(id)

      matched_row <- NA_integer_

      for (k in id_keys) {
        if (exists(k, envir = lookup, inherits = FALSE)) {
          matched_row <- get(k, envir = lookup, inherits = FALSE)
          break
        }
      }

      if (!is.na(matched_row)) {
        len <- clean_num(len_col[matched_row])

        acc <- if (has_value(accession_col[matched_row])) {
          accession_col[matched_row]
        } else {
          NA_character_
        }

        alias <- if (has_value(alias_col[matched_row])) {
          alias_col[matched_row]
        } else {
          NA_character_
        }

        mk_attr(
          id    = id,
          len   = len,
          acc   = acc,
          alias = alias,
          tax   = tax
        )
      } else {
        mk_attr(
          id    = id,
          len   = NA_real_,
          acc   = NA_character_,
          alias = NA_character_,
          tax   = tax
        )
      }
    },
    FUN.VALUE = character(1)
  )
}



# If raws genotypes is not a vcf file  
.open_text_con <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    return(gzfile(path, open = "rt"))
  } else {
    return(file(path, open = "rt"))
  }
}

# Helper: infer separator from filename (and fallback to sniffing first line)
.infer_sep <- function(path) {
  # strip a trailing .gz (if present) before checking extension
  base <- sub("\\.gz$", "", basename(path), ignore.case = TRUE)
  ext  <- tolower(tools::file_ext(base))

  if (ext %in% c("tsv", "tab")) return("\t")
  if (ext %in% c("csv"))        return(",")
  if (ext %in% c("txt")) {
    # Try to guess for .txt (could be tab, comma, semicolon, or whitespace)
    con <- .open_text_con(path)
    on.exit(close(con), add = TRUE)
    first_line <- readLines(con, n = 1, warn = FALSE)

    # count candidate delimiters
    n_tab <- lengths(regmatches(first_line, gregexpr("\t", first_line, fixed = TRUE)))
    n_com <- lengths(regmatches(first_line, gregexpr(",",  first_line, fixed = TRUE)))
    n_sem <- lengths(regmatches(first_line, gregexpr(";",  first_line, fixed = TRUE)))

    best <- max(c(n_tab, n_com, n_sem), na.rm = TRUE)
    if (best == 0) return("")          # whitespace-delimited
    if (best == n_tab) return("\t")
    if (best == n_com) return(",")
    return(";")
  }

  # Unknown extension: sniff first line (same logic as .txt)
  con <- .open_text_con(path)
  on.exit(close(con), add = TRUE)
  first_line <- readLines(con, n = 1, warn = FALSE)

  n_tab <- lengths(regmatches(first_line, gregexpr("\t", first_line, fixed = TRUE)))
  n_com <- lengths(regmatches(first_line, gregexpr(",",  first_line, fixed = TRUE)))
  n_sem <- lengths(regmatches(first_line, gregexpr(";",  first_line, fixed = TRUE)))

  best <- max(c(n_tab, n_com, n_sem), na.rm = TRUE)
  if (best == 0) return("")
  if (best == n_tab) return("\t")
  if (best == n_com) return(",")
  return(";")
}

# Main reader: base R, works for plain/gz, tsv/csv/txt + delimiter sniffing
# Updated delimiter detection to not rely on any file extension including zipped status (note, it searches whitespace vs tab as delimiter. Anything else may stuff up from sample names e.g., inclusion of | 
## so we needs to be in simple formats. 


read_genotypes_auto <- function(path) {

  # Helper to open either gzipped or plain-text file
  open_connection <- function(path) {
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
      gzfile(path, open = "rt")
    } else {
      file(path, open = "rt")
    }
  }

  # -----------------------------
  # Step 1: Peek first non-empty line
  # -----------------------------
  con_peek <- open_connection(path)

  first_line <- ""
  repeat {
    line <- readLines(con_peek, n = 1, warn = FALSE)

    if (length(line) == 0) {
      close(con_peek)
      stop("Input file appears to be empty: ", path)
    }

    if (nzchar(trimws(line))) {
      first_line <- line
      break
    }
  }

  close(con_peek)

  # -----------------------------
  # Step 2: Infer delimiter
  # -----------------------------
  if (grepl("\t", first_line)) {
    sep <- "\t"
    delim_name <- "TAB"
  } else if (grepl("\\s+", first_line)) {
    sep <- ""
    delim_name <- "WHITESPACE"
  } else {
    stop(
      "Unable to infer delimiter (expected tab or whitespace) for file: ", path,
      "\nHeader line: ", first_line
    )
  }

  # -----------------------------
  # Step 3: Reopen and read data
  # -----------------------------
  con_read <- open_connection(path)
  on.exit(try(close(con_read), silent = TRUE), add = FALSE)

  dat <- tryCatch(
    read.table(
      con_read,
      header = TRUE,
      sep = sep,
      check.names = FALSE,
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE
    ),
    error = function(e) {
      stop(
        "Failed to parse genotype file: ", path,
        "\nInferred delimiter: ", delim_name,
        "\nOriginal error: ", e$message
      )
    }
  )

  return(dat)
}



# If raw genotypes input is not a vcf file
if (!is_vcfraws & !isDArTfile) {
  genotypesRAWS <- read_genotypes_auto(raw_path)
}

Finalmappings  <- read.csv(map_path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "#")

MAPPRIORS <- NULL
if (nzchar(mappriors_path)) {
  MAPPRIORS <- read.table(
    mappriors_path,
    sep = "\t",
    header = TRUE
  )
}


# Newly added metadata in lines 1-10 of Brioche (might expand further)
# Ln 1= generic, Ln2 brioche commit version, Ln3 brioche repo, Ln4 brioche branch, ln5 pident ln6 coverage, ln7 blast options, ln8 target chromprior? ln9 sharedmarkersmap used?,ln10 ld mapping used?, ln11 genetic map used?, ln 12 rundate
Metadata <- readLines(map_path,n=12L)
# remove excess comment chars 
Metadata <- sub(pattern="## ", replacement="",x=Metadata, fixed =TRUE)
# Remove extra column at start
Metadata <- Metadata[-1]

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

# Contig lines (always full-shape; uses overrides only when NCBI meta is missing)
VCFchroms <- order_contigs(Finalmappings$saccver)
contig_lines <- build_contig_lines(
  ids   = VCFchroms,
  meta  = ncbi_meta,
  override_species = if (nzchar(user_species)) user_species else NULL
)



# SNP-chip pathway (raw genotypes table to vcf)

if (!is_vcfraws & is_snpchip) {

  # ---------- helpers ----------
  strand_to_word <- function(x) {
    z <- tolower(trimws(as.character(x)))
    ifelse(is.na(z) | z == "" | z == ".", "none",
      ifelse(z %in% c("+","plus","pos","positive","p"), "plus",
        ifelse(z %in% c("-","minus","neg","negative","m"), "minus", "none")))
  }
  is_blank <- function(x) is.na(x) | x == "" | x == "."

  #  precompute MAPPRIORS lookups
  pri_flag_by_id <- logical(0)     # named logical: qaccver -> TRUE when qualifies for "Unique_mapping_with_priors_information"
  dup_label_by_id <- character(0)  # named char:   qaccver -> "LocallyDuplicatedRegion" | "SingleCopy"

  if (is.data.frame(MAPPRIORS)) {
    n  <- nrow(MAPPRIORS)
    tp <- MAPPRIORS
    key <- as.character(tp$qaccver)

    # normalize to vector length n, even if column is NULL/absent
    is_true_chr <- function(v) tolower(trimws(as.character(v))) == "true"
    norm_true   <- function(col) if (is.null(col)) rep(FALSE, n) else rep_len(is_true_chr(col), n)
    norm_is_na  <- function(col) if (is.null(col)) rep(TRUE,  n) else rep_len(is.na(col),        n)

    # prior-based uniqueness flag: Trueunique is NA AND any prior* == "true"
    tu_na <- norm_is_na(tp$Trueunique)
    any_prior_true <- norm_true(tp$Chrommarkermap) |
                      norm_true(tp$proximatemarkermap) |
                      norm_true(tp$linkagemarkermap) |
                      norm_true(tp$geneticmapmarkermap)

    # safe setNames helper (no mismatch crashes)
    safe_set_names <- function(x, nm) {
      lx <- length(x); ln <- length(nm)
      if (!lx || !ln) return(x)
      if (lx == ln)   return(stats::setNames(x, nm))
      k <- min(lx, ln)
      stats::setNames(x[seq_len(k)], nm[seq_len(k)])
    }

    pri_flag_by_id  <- safe_set_names(tu_na & any_prior_true, key)

    # duplication label
    dup_is_true     <- norm_true(tp$Duplicate_region)
    dup_label_by_id <- safe_set_names(
      ifelse(dup_is_true, "LocallyDuplicatedRegion", "SingleCopy"),
      key
    )
  }

  # ---------- header scaffold ----------
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
  fileDate  <- format(Sys.Date(), "%Y%m%d")
  VCFchroms <- unique(c(order_contigs(Finalmappings$saccver), "chrUnk"))

  header_lines <- c(
    "##fileformat=VCFv4.2",
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
    '##FILTER=<ID=LowQual,Description="Low quality or ambiguous mapping">',
    '##INFO=<ID=MAPSTATUS,Number=1,Type=String,Description="Unique_mapping|Unique_mapping_with_priors_information|Failed_to_map_uniquely">',
    '##INFO=<ID=PriorORIENT,Number=1,Type=String,Description="Accumulated strand orientation across runs (plus|minus|none)">',
    '##INFO=<ID=DUP,Number=1,Type=String,Description="SingleCopy|LocallyDuplicatedRegion">',
    '##INFO=<ID=MAF,Number=A,Type=Float,Description="Minor Allele frequency per ALT">',
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count per ALT">',
    '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=NU,Number=1,Type=Integer,Description="Null allele flag from SNP chip (1= Null allele; 0=otherwise)">',
    '##FORMAT=<ID=DU,Number=1,Type=Integer,Description="Number of local duplications of marker detected (0= Marker is single copy; 1+=copies, .=Unknown)">',
    contig_lines
  )
  header_lines <- header_lines[!duplicated(header_lines)]
  writeLines(header_lines, con = Outputfilename)

  vcf_header <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", sample_ids)
  write_header_row(vcf_header, Outputfilename)

  #2 Stream rows
  chunk_size <- 1500L
  con_chunk <- if (grepl("\\.gz$", raw_path, ignore.case = TRUE)) gzfile(raw_path, "rt") else file(raw_path, "rt")
  on.exit(try(close(con_chunk), silent = TRUE), add = TRUE)
  invisible(readLines(con_chunk, n = 1L, warn = FALSE))  # consume header

  repeat {
    raw_lines <- readLines(con_chunk, n = chunk_size, warn = FALSE)
    if (!length(raw_lines)) break

    # NEW: drop empty/whitespace-only lines (common at the tail of files)
    raw_lines <- raw_lines[ nzchar(raw_lines) & !grepl("^\\s*$", raw_lines) ]
    if (!length(raw_lines)) next

    df <- utils::read.table(
      text = raw_lines, sep = delim, header = FALSE,
      quote = "", comment.char = "", stringsAsFactors = FALSE,
      check.names = FALSE, col.names = cn_fixed, fill = TRUE
    )

    # NEW: if cleaning leaves no rows, skip quietly
    if (!NROW(df)) next

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

    # Orientation decision
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

    # QUAL/FILTER/INFO (with MAPPRIORS + DUP)
    pos_int     <- suppressWarnings(as.integer(POS))
    coord_bad   <- is.na(pos_int) | pos_int <= 0L | is_blank(CHROM)
    unknown_now <- !u_map | coord_bad

    QUAL   <- ifelse(unknown_now, 0L, 100L)
    FILTER <- ifelse(unknown_now, "LowQual", "PASS")

    # priors flag lookup per Name
    if (length(pri_flag_by_id)) {
      pri_flag <- unname(pri_flag_by_id[Name])
      pri_flag[is.na(pri_flag)] <- FALSE
    } else {
      pri_flag <- rep(FALSE, length(Name))
    }

    # status: only upgrade to "Unique_mapping_with_priors_information" when we *mapped* this round
    status_val <- ifelse(
      unknown_now, "Failed_to_map_uniquely",
      ifelse(pri_flag, "Unique_mapping_with_priors_information", "Unique_mapping")
    )

    # DUP label per Name
    if (length(dup_label_by_id)) {
      dup_lab <- unname(dup_label_by_id[Name])
      dup_lab[is.na(dup_lab)] <- "SingleCopy"
    } else {
      dup_lab <- rep("SingleCopy", length(Name))
    }

    # PriorORIENT from current sstrand when uniquely mapped; else 'none'
    cur_strand <- ifelse(u_map, strand_to_word(fm$STRAND[m]), "none")

    INFO <- paste0("MAPSTATUS=", status_val, ";PriorORIENT=", cur_strand, ";DUP=", dup_lab)

    ## Genotypes + NU + DU (DU presently unknown => '.')
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
      M <- matrix(rep(rswap, times = ncol(geno)), nrow = nrow(geno))
      d0 <- (g == "0") & M; if (any(d0)) g[d0] <- "1/1"
      d1 <- (g == "1") & M; if (any(d1)) g[d1] <- "0/1"
      d2 <- (g == "2") & M; if (any(d2)) g[d2] <- "0/0"
      geno[M] <- g[M]
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
    NU_chr <- matrix(as.character(NU), nrow = nrow(geno), ncol = ncol(geno)); NU_chr[is.na(NU_chr)] <- "0"
    DU_chr <- matrix(".", nrow = nrow(geno), ncol = ncol(geno))  # unknown duplication count per-sample for now

    geno_df <- as.data.frame(geno, stringsAsFactors = FALSE, check.names = FALSE)
    for (jj in seq_len(ncol(geno_df))) {
      geno_df[[jj]] <- paste0(geno_df[[jj]], ":", NU_chr[, jj], ":", DU_chr[, jj])
    }

    body <- data.frame(
      `#CHROM` = CHROM, POS = POS, ID = Name,
      REF = REF_out, ALT = ALT_out, QUAL = QUAL, FILTER = FILTER, INFO = INFO, FORMAT = "GT:NU:DU",
      stringsAsFactors = FALSE, check.names = FALSE
    )
    body <- cbind(body, geno_df, stringsAsFactors = FALSE)

    if (NROW(body)) {
      dt_fwrite(body, file = Outputfilename, append = TRUE, col.names = FALSE)
    }

    rm(df, geno, geno_df, NU, NU_chr, body, ac1, an); gc(FALSE)
  }

  message(sprintf("Wrote: %s", Outputfilename))
}



# VCF pathway (VCF -> re-anchored VCF using Brioche sstrand and PriorORIENT)
if (is_vcfraws) {

  # ---- small helpers ----
  strand_to_word <- function(x) {
    z <- tolower(trimws(as.character(x)))
    ifelse(
      is.na(z) | z == "" | z == ".",
      NA_character_,
      ifelse(
        z %in% c("+","plus","pos","positive","p"), "plus",
        ifelse(z %in% c("-","minus","neg","negative","m"), "minus", NA_character_)
      )
    )
  }
  pr_norm <- function(x) {
    z <- tolower(trimws(as.character(x)))
    ifelse(
      is.na(z) | z == "" | z == ".",
      "none",
      ifelse(
        z %in% c("plus","+","pos"), "plus",
        ifelse(z %in% c("minus","-","neg"), "minus", "none")
      )
    )
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
  add_if_missing <- function(vec, line_glob, line_exact) {
    if (!any(grepl(line_glob, vec))) c(vec, line_exact) else vec
  }

  # header reordering helper as it's getting a bit out of order in insilico
  reorder_vcf_header <- function(meta_lines, contig_lines) {
    meta <- meta_lines[!duplicated(meta_lines)]
    pick <- function(pat, ignore.case = TRUE) which(grepl(pat, meta, ignore.case = ignore.case))

    i_fileformat   <- pick("^##fileformat=")
    i_fileDate     <- pick("^##filedate=")
    i_phasing      <- pick("^##phasing=")

    i_bri_ver      <- pick("^##brioche_version")
    i_bri_repo     <- pick("^##brioche_repo=")
    i_bri_branch   <- pick("^##brioche_branch=")

    idx_taken      <- unique(c(i_bri_ver, i_bri_repo, i_bri_branch))
    i_bri_all      <- setdiff(pick("^##brioche"), idx_taken)

    i_filter       <- pick("^##filter")
    i_info         <- pick("^##info")
    i_format       <- pick("^##format")
    i_reference    <- pick("^##reference=")
    i_assembly     <- pick("^##assembly=")
    i_genomeurl    <- pick("^##genome_url=")
    i_bcftools     <- pick("^##bcftools")

    i_all   <- seq_along(meta)
    i_other <- setdiff(i_all, unique(c(
      i_fileformat, i_fileDate, i_phasing,
      i_bri_ver, i_bri_repo, i_bri_branch, i_bri_all,
      i_filter, i_info, i_format,
      i_reference, i_assembly, i_genomeurl, i_bcftools
    )))

    out <- c(
      meta[i_fileformat],
      meta[i_fileDate],
      meta[i_phasing],
      meta[i_bri_ver],
      meta[i_bri_repo],
      meta[i_bri_branch],
      meta[i_bri_all],
      meta[i_filter],
      meta[i_info],
      meta[i_format],
      meta[i_reference],
      meta[i_assembly],
      meta[i_genomeurl],
      contig_lines,
      meta[i_bcftools],
      meta[i_other]
    )

    out <- out[nzchar(out)]
    out[!duplicated(out)]
  }

  # ---- priors lookups (safe against length/name mismatch) ----
  pri_flag_by_id  <- logical(0)   # named logical: qaccver -> TRUE for "Unique_mapping_with_priors_information"
  dup_label_by_id <- character(0) # named char:   qaccver -> "LocallyDuplicatedRegion" | "SingleCopy"

  if (is.data.frame(MAPPRIORS)) {
    n  <- nrow(MAPPRIORS)
    tp <- MAPPRIORS
    key <- as.character(tp$qaccver)

    is_true_chr <- function(v) tolower(trimws(as.character(v))) == "true"
    norm_true   <- function(col) if (is.null(col)) rep(FALSE, n) else rep_len(is_true_chr(col), n)
    norm_is_na  <- function(col) if (is.null(col)) rep(TRUE,  n) else rep_len(is.na(col),        n)

    tu_na <- norm_is_na(tp$Trueunique)
    any_prior_true <- norm_true(tp$Chrommarkermap) |
                      norm_true(tp$proximatemarkermap) |
                      norm_true(tp$linkagemarkermap) |
                      norm_true(tp$geneticmapmarkermap)

    safe_set_names <- function(x, nm) {
      lx <- length(x); ln <- length(nm)
      if (!lx || !ln) return(x)
      if (lx == ln)   return(stats::setNames(x, nm))
      k <- min(lx, ln)
      stats::setNames(x[seq_len(k)], nm[seq_len(k)])
    }

    pri_flag_by_id  <- safe_set_names(tu_na & any_prior_true, key)
    dup_is_true     <- norm_true(tp$Duplicate_region)
    dup_label_by_id <- safe_set_names(ifelse(dup_is_true, "LocallyDuplicatedRegion", "SingleCopy"), key)
  }

  req <- c("qaccver","saccver","SNPpos","Ref","ALT","sstrand")
  missing <- setdiff(req, names(Finalmappings))
  if (length(missing))
    stop(sprintf("Finalmappings missing: %s", paste(missing, collapse=", ")))

  fm_one <- Finalmappings %>%
    dplyr::transmute(
      qaccver,
      saccver,
      SNPpos,
      REF_map = toupper(Ref),
      ALT_map = normalize_alt(ALT),
      STRAND  = sstrand
    ) %>%
    dplyr::group_by(qaccver) %>% dplyr::slice(1L) %>% dplyr::ungroup()
  fm_hits <- Finalmappings %>% dplyr::count(qaccver, name = "MAP_HITS")

  fm <- fm_one %>%
    dplyr::left_join(fm_hits, by = "qaccver") %>%
    dplyr::transmute(
      ID        = qaccver,
      CHROM_map = saccver,
      POS_map   = SNPpos,
      REF_map   = REF_map,
      ALT_map   = ALT_map,
      STRAND    = STRAND,
      MAP_HITS  = ifelse(is.na(MAP_HITS), 0L, MAP_HITS)
    )

  # read input heade 
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

  # ---- rebuild meta header: ensure GT/NU/DU present; standard Brioche/INFO/FILTER; reorder ----
  fileDate <- format(Sys.Date(), "%Y%m%d")

  # strip older versions of things we control
  meta_lines <- meta_lines[!grepl("^##(fileDate=)",            meta_lines)]
  meta_lines <- meta_lines[!grepl("^##(Brioche_)",             meta_lines)]
  meta_lines <- meta_lines[!grepl("^##(assembly=)",            meta_lines)]
  meta_lines <- meta_lines[!grepl("^##(reference=)",           meta_lines)]
  meta_lines <- meta_lines[!grepl("^##(genome_url=)",          meta_lines)]
  meta_lines <- meta_lines[!grepl("^##(contig=)",              meta_lines)]
  meta_lines <- meta_lines[!grepl("^##(fileformat=)",          meta_lines)]
  meta_lines <- meta_lines[!grepl("^##INFO=<ID=MAPSTATUS,",     meta_lines)]
  meta_lines <- meta_lines[!grepl("^##INFO=<ID=DUP,",           meta_lines)]
  meta_lines <- meta_lines[!grepl("^##INFO=<ID=PriorORIENT,",   meta_lines)]
  meta_lines <- meta_lines[!grepl("^##FILTER=<ID=LowQual,",     meta_lines)]
  meta_lines <- meta_lines[!grepl("^##FORMAT=<ID=DU,",          meta_lines)]
  meta_lines <- meta_lines[!grepl("^##FORMAT=<ID=NU,",          meta_lines)]
  # keep any existing GT line; we'll ensure it exists below if missing

  meta_out <- c(
    "##fileformat=VCFv4.2",
    sprintf("##fileDate=%s", fileDate),
    meta_lines,
    paste0("##", Metadata[1]),
    paste0("##", Metadata[2]),
    paste0("##", Metadata[3]),
    paste0("##Brioche_filtering_", Metadata[4]),
    paste0("##Brioche_filtering_", Metadata[5]),
    paste0("##Brioche_filtering_", Metadata[6]),
    paste0("##Brioche_priors_files:", Metadata[7]),
    paste0("##Brioche_priors_files:", Metadata[8]),
    paste0("##Brioche_priors_files:", Metadata[9]),
    paste0("##Brioche_priors_files:", Metadata[10]),
    paste0("##Brioche_", Metadata[11])
  )
  meta_out <- meta_out[!duplicated(meta_out)]

  # enforce presence of our INFO/FILTER/FORMAT lines
  meta_out <- add_if_missing(meta_out, '^##FILTER=<ID=LowQual,', '##FILTER=<ID=LowQual,Description="Low quality or ambiguous mapping">')
  meta_out <- add_if_missing(meta_out, '^##INFO=<ID=MAPSTATUS,', '##INFO=<ID=MAPSTATUS,Number=1,Type=String,Description="Unique_mapping|Unique_mapping_with_priors_information|Failed_to_map_uniquely">')
  meta_out <- add_if_missing(meta_out, '^##INFO=<ID=PriorORIENT,', '##INFO=<ID=PriorORIENT,Number=1,Type=String,Description="Accumulated strand orientation across runs (plus|minus|none)">')
  meta_out <- add_if_missing(meta_out, '^##INFO=<ID=DUP,', '##INFO=<ID=DUP,Number=1,Type=String,Description="SingleCopy|LocallyDuplicatedRegion">')
  meta_out <- add_if_missing(meta_out, '^##FORMAT=<ID=GT,', '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
  meta_out <- add_if_missing(meta_out, '^##FORMAT=<ID=NU,', '##FORMAT=<ID=NU,Number=1,Type=Integer,Description="Null allele flag from SNP chip (1= Null allele; 0=otherwise)">')
  meta_out <- add_if_missing(meta_out, '^##FORMAT=<ID=DU,', '##FORMAT=<ID=DU,Number=1,Type=Integer,Description="Number of local duplications of marker detected (0= Marker is single copy; 1+=copies, .=Unknown)">')

  # assemble + reorder
  meta_all <- unique(c(meta_out, extra_meta))
  header_ordered <- reorder_vcf_header(meta_all, contig_lines)
  writeLines(c(header_ordered, header_line), con = Outputfilename)

  # ---- stream & transform records ----
  ct_all_as_char <- setNames(rep("character", length(colnames_vcf)), colnames_vcf)
  mf <- if (requireNamespace("fastmatch", quietly = TRUE)) fastmatch::fmatch else base::match

  con <- if (grepl("\\.gz$", raw_in, ignore.case = TRUE)) gzfile(raw_in, "rt") else file(raw_in, "rt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  if (header_n > 0L) invisible(readLines(con, n = header_n, warn = FALSE))

  chunk_size <- getOption("ANCHOR_VCF_CHUNK", 5000L)

  read_chunk_df <- function(n) {
    lines <- readLines(con, n = n, warn = FALSE)
    if (!length(lines)) return(NULL)
    txt <- paste(lines, collapse = "\n")
    utils::read.table(
      text         = txt,
      sep          = "\t",
      header       = FALSE,
      col.names    = colnames_vcf,
      colClasses   = "character",
      quote        = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      check.names  = FALSE,
      fill         = TRUE
    )
  }

  repeat {
    df <- read_chunk_df(chunk_size)
    if (is.null(df) || !nrow(df)) break
    df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)

    IDv    <- df[[ix_ID]]
    # map status
    j          <- mf(IDv, fm$ID)
    in_map     <- !is.na(j)
    hits       <- ifelse(in_map, fm$MAP_HITS[j], 0L)
    chr_ok     <- in_map & nzchar_safe(fm$CHROM_map[j])
    pos_ok     <- in_map & !is.na(fm$POS_map[j]) & suppressWarnings(as.numeric(fm$POS_map[j]) > 0)
    mapped_now <- in_map & (hits == 1L) & chr_ok & pos_ok

    # CHROM/POS/QUAL/FILTER (overwrite FILTER always)
    df[[ix_CHROM]] <- ifelse(mapped_now, fm$CHROM_map[j], "chrUnk")
    df[[ix_POS]]   <- ifelse(mapped_now, as.character(fm$POS_map[j]), "0")
    if (!is.na(ix_QUAL)) df[[ix_QUAL]] <- ifelse(mapped_now, "100", "0")
    if (!is.na(ix_FIL )) df[[ix_FIL ]] <- ifelse(mapped_now, "PASS", "LowQual")

    # prior vs current strand
    info_in    <- if (!is.na(ix_INFO)) df[[ix_INFO]] else rep("", nrow(df))
    prior_raw  <- info_get(info_in, "PriorORIENT")
    prior      <- pr_norm(prior_raw)
    s_now      <- strand_to_word(ifelse(mapped_now, fm$STRAND[j], NA_character_))

    # newly mapped => set current orientation, else keep prior
    prior_out  <- prior
    set_pr     <- !is.na(s_now)
    prior_out[set_pr] <- s_now[set_pr]

    # allele updates from Brioche
    BriocheREF   <- ifelse(mapped_now, toupper(fm$REF_map[j]), NA_character_)
    BriocheALT   <- ifelse(mapped_now, first_alt_uc(fm$ALT_map[j]), NA_character_)
    BriocheREFRC <- revcomp_iupac(BriocheREF)

    # Flip REF/ALT if prior != none and orientation changed
    need_flip <- mapped_now & prior != "none" & !is.na(s_now) & prior != s_now
    if (any(need_flip)) {
      REF_cur <- df[[ix_REF]]
      ALT_cur <- df[[ix_ALT]]
      REF_cur[need_flip] <- revcomp_iupac(REF_cur[need_flip])
      if (any(need_flip)) {
        ALT_cur[need_flip] <- vapply(ALT_cur[need_flip], rc_csv, character(1))
      }
      df[[ix_REF]] <- REF_cur
      df[[ix_ALT]] <- ALT_cur
    }

    prior_eff <- prior
    prior_eff[need_flip] <- s_now[need_flip]

    # FIRST mapping: treat as biallelic initial orientation when prior == none
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
      tmpREF <- df[[ix_REF]]
      tmpALT <- df[[ix_ALT]]
      tmpREF[idx_bri] <- BriocheREF[idx_bri]
      tmpALT[idx_bri] <- BriocheALT[idx_bri]
      df[[ix_REF]] <- tmpREF
      df[[ix_ALT]] <- tmpALT
    }

    # SUBSEQUENT mapping (same orientation) — allow tri+ remapping
   tri_rows <- which(mapped_now & prior_eff != "none" & !is.na(s_now) & (prior_eff == s_now) & nzchar(BriocheREF))
    idx_map_by_row <- vector("list", length(tri_rows))
    names(idx_map_by_row) <- as.character(tri_rows)
    recode_needed  <- logical(length(tri_rows))

    if (length(tri_rows)) {
      REF_cur2 <- df[[ix_REF]]
      ALT_cur2 <- df[[ix_ALT]]

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
    }

    #  Build INFO fresh (overwrite any prior content)
    status_val <- ifelse(mapped_now, "Unique_mapping", "Failed_to_map_uniquely")
    if (length(pri_flag_by_id)) {
      pri_flag <- unname(pri_flag_by_id[IDv]); pri_flag[is.na(pri_flag)] <- FALSE
      status_val[mapped_now & pri_flag] <- "Unique_mapping_with_priors_information"
    }
    dup_lab <- if (length(dup_label_by_id)) {
      out <- unname(dup_label_by_id[IDv]); out[is.na(out)] <- "SingleCopy"; out
    } else {
      rep("SingleCopy", length(IDv))
    }
    # enforce overwrite of INFO (ignore any existing values)
    if (!is.na(ix_INFO)) {
      df[[ix_INFO]] <- paste0(
        "MAPSTATUS=", status_val,
        ";PriorORIENT=", ifelse(is.na(s_now), prior_eff, s_now),
        ";DUP=", dup_lab
      )
    }

    #- Normalize FORMAT to exactly "GT:NU:DU" and rebuild samples 
    if (!is.na(ix_FMT) && length(sample_cols_idx)) {

      # Cache incoming FORMAT before overwriting it
      fmt_in <- as.character(df[[ix_FMT]])
      fmt_in[is.na(fmt_in) | fmt_in == ""] <- "GT"

      # For speed: if FORMAT is constant across this chunk, compute positions once.
      # Otherwise compute per-row positions (still correct).
      fmt_unique <- unique(fmt_in)

      if (length(fmt_unique) == 1L) {
        fields <- strsplit(fmt_unique, ":", fixed = TRUE)[[1]]
        pos_NU_const <- match("NU", fields)
        pos_DU_const <- match("DU", fields)
        pos_GT_const <- match("GT", fields); if (is.na(pos_GT_const)) pos_GT_const <- 1L

        pos_NU_by_row <- rep(pos_NU_const, nrow(df))
        pos_DU_by_row <- rep(pos_DU_const, nrow(df))
        pos_GT_by_row <- rep(pos_GT_const, nrow(df))
      } else {
        split_fields <- strsplit(fmt_in, ":", fixed = TRUE)
        pos_NU_by_row <- vapply(split_fields, function(x) match("NU", x), integer(1))
        pos_DU_by_row <- vapply(split_fields, function(x) match("DU", x), integer(1))
        pos_GT_by_row <- vapply(split_fields, function(x) { p <- match("GT", x); if (is.na(p)) 1L else p }, integer(1))
      }

      # We will output a standard FORMAT for all rows
      df[[ix_FMT]] <- "GT:NU:DU"

      # rows to recode (same as your logic)
      if (length(swap_rows) || length(tri_rows)) {
        rows_to_recode <- tri_rows[recode_needed]
      } else {
        rows_to_recode <- integer(0)
      }

      for (jjj in sample_cols_idx) {

        colv <- df[[jjj]]
        colv[is.na(colv)] <- "."

        # Split existing sample field once (fast C implementation)
        # returns an N x K matrix (K = max fields present in this sample column)
        parts <- stringi::stri_split_fixed(colv, ":", simplify = TRUE)

        # Helper to pull a field by (possibly variable) position per row
        pull_by_pos <- function(parts_mat, pos_vec) {
          out <- rep(NA_character_, nrow(parts_mat))
          ok <- !is.na(pos_vec) & pos_vec >= 1L & pos_vec <= ncol(parts_mat)
          if (any(ok)) {
            ii <- which(ok)
            out[ii] <- parts_mat[cbind(ii, pos_vec[ii])]
          }
          out[out == ""] <- NA_character_
          out
        }

        # Extract existing NU/DU if present in the incoming FORMAT
        nu_old <- pull_by_pos(parts, pos_NU_by_row)
        du_old <- pull_by_pos(parts, pos_DU_by_row)

        # Default only when NU/DU truly absent from FORMAT
        # (if NU/DU present but missing in a cell, we treat as 0 / . respectively)
        if (all(is.na(pos_NU_by_row))) {
          nu_out <- rep("0", nrow(df))
        } else {
          nu_out <- ifelse(is.na(nu_old) | nu_old == ".", "0", nu_old)
        }

        if (all(is.na(pos_DU_by_row))) {
          du_out <- rep(".", nrow(df))
        } else {
          du_out <- ifelse(is.na(du_old) | du_old == "", ".", du_old)
        }

        # --- Build GT exactly like you already do ---
        gt <- sub("^([^:]+).*", "\\1", colv, perl = TRUE)
        gt[gt == "" | gt == "."] <- "./."
        gt <- gsub("\\|", "/", gt)
        gt[gt == "0"] <- "0/0"; gt[gt == "1"] <- "0/1"; gt[gt == "2"] <- "1/1"

        # apply biallelic swap (only to GT)
        if (length(swap_rows)) {
          tmp <- gt
          tmp[swap_rows] <- swap_biallelic_gt_vec(tmp[swap_rows])
          gt <- tmp
        }

        # apply tri+ recode where REF changed
        if (length(rows_to_recode)) {
          for (rr in rows_to_recode) {
            idx_map <- idx_map_by_row[[as.character(rr)]]
            if (!is.null(idx_map)) gt[rr] <- recode_gt_vec(gt[rr], idx_map)
          }
        }

        # compose final GT:NU:DU (preserving NU/DU)
        df[[jjj]] <- paste0(gt, ":", nu_out, ":", du_out)
      }
    }
    # write out
    dt_fwrite(df, file = Outputfilename, append = TRUE, col.names = FALSE)

    rm(df); gc(FALSE)
  }

  message(sprintf("Wrote: %s", Outputfilename))
}








# DArT pathway (raw DArT table -> VCF; record PriorORIENT from sstrand)
# DArT pathway (raw DArT table -> VCF; record PriorORIENT from sstrand)
if (isDArTfile) {

  # ----------------------------
  # DArT helpers
  # ----------------------------
  detect_dart_layout_base <- function(path, max_scan = 400L, n_sample = 25L, geno_prop = 0.80) {

    .open_text_con_local <- function(p) {
      if (grepl("\\.gz$", p, ignore.case = TRUE)) gzfile(p, open = "rt") else file(p, open = "rt")
    }

    con <- .open_text_con_local(path)
    on.exit(try(close(con), silent = TRUE), add = TRUE)

    lines <- readLines(con, n = max_scan, warn = FALSE)
    if (!length(lines)) stop(sprintf("DArT file appears empty: %s", path))
    lines <- sub("^\ufeff", "", lines, perl = TRUE)

    # Header line must contain an AlleleID token (SNP may or may not exist)
    has_alleleid_token <- function(ln) {
      grepl("(^|[,\t])\\s*allele[_ ]?id\\s*([,\t]|$)", ln, ignore.case = TRUE, perl = TRUE)
    }
    hdr_i <- which(vapply(lines, has_alleleid_token, logical(1)))[1]
    if (is.na(hdr_i)) {
      stop("Could not find a DArT header line containing column AlleleID (case-insensitive).")
    }

    header_line <- lines[hdr_i]

    # delimiter guess from header line
    n_tab <- lengths(regmatches(header_line, gregexpr("\t", header_line, fixed = TRUE)))
    n_com <- lengths(regmatches(header_line, gregexpr(",",  header_line, fixed = TRUE)))
    delim <- if (n_tab >= n_com && n_tab > 0) "\t" else if (n_com > 0) "," else "\t"

    # parse column names
    cn <- strsplit(header_line, delim, fixed = TRUE)[[1]]
    cn <- gsub('^"|"$', "", trimws(cn))

    # canonicalize only *exact* names (hard match; do NOT match SNPsomething)
    cn_trim_lc <- tolower(trimws(cn))
    cn[cn_trim_lc %in% c("alleleid", "allele_id", "allele id")] <- "AlleleID"
    cn[cn_trim_lc == "snp"] <- "SNP"

    # read a few data lines after header to infer where genotype columns start
    data_lines <- lines[(hdr_i + 1L):length(lines)]
    data_lines <- data_lines[nzchar(data_lines) & !grepl("^\\s*$", data_lines)]
    if (!length(data_lines)) stop("DArT file has a header but no data rows.")
    data_lines <- utils::head(data_lines, n_sample)

    split_row <- function(ln) {
      v <- strsplit(ln, delim, fixed = TRUE)[[1]]
      v <- gsub('^"|"$', "", trimws(v))
      if (length(v) < length(cn)) v <- c(v, rep("", length(cn) - length(v)))
      if (length(v) > length(cn)) v <- v[seq_len(length(cn))]
      v
    }

    mat <- do.call(rbind, lapply(data_lines, split_row))
    matU <- toupper(trimws(mat))

    is_geno_like <- function(x) {
      x <- toupper(trimws(x))
      x == "" | x %in% c(".", "./.", "NA", "N", "NC", "-") |
        grepl("^[0-3]$", x, perl = TRUE) |
        grepl("^[0-3][\\/|][0-3]$", x, perl = TRUE)
    }

    props <- vapply(seq_len(ncol(matU)), function(j) mean(is_geno_like(matU[, j])), numeric(1))

    j <- length(props)
    while (j >= 1L && is.finite(props[j]) && props[j] >= geno_prop) j <- j - 1L
    g_start <- j + 1L

    # fallback: if heuristic fails, start right after SNP if present, else after AlleleID
    if (g_start > length(cn)) {
      idx_snp    <- match("SNP", cn)        # exact only
      idx_allele <- match("AlleleID", cn)
      if (!is.na(idx_snp)) {
        g_start <- idx_snp + 1L
      } else if (!is.na(idx_allele)) {
        g_start <- idx_allele + 1L
      } else {
        stop("Header parsing failed to find AlleleID column.")
      }
    }

    lastmetric_idx <- g_start - 1L
    if (lastmetric_idx < 1L) lastmetric_idx <- 1L

    list(
      cn = cn,
      topskip = hdr_i - 1L,
      lastmetric_idx = lastmetric_idx,
      delim = delim
    )
  }

  # Extract REF/ALT from:
  #  1) SNP column format:  position:REF>ALT  (uses last ':' before REF)
  #  2) If SNP column absent OR row doesn't parse: AlleleID suffix :REF>ALT$
  dart_extract_ref_alt <- function(snp_vec, allele_vec) {
    snp_vec    <- as.character(snp_vec);    snp_vec[is.na(snp_vec)] <- ""
    allele_vec <- as.character(allele_vec); allele_vec[is.na(allele_vec)] <- ""

    # position:REF>ALT  (hard anchor at end)
    m_snp <- stringr::str_match(toupper(trimws(snp_vec)), ".*:([A-Z]+)>([A-Z]+)$")
    ref_s <- m_snp[, 2]
    alt_s <- m_snp[, 3]

    ok_snp <- !is.na(ref_s) & !is.na(alt_s) & nzchar(ref_s) & nzchar(alt_s)

    # AlleleID suffix :REF>ALT$
    m_id  <- stringr::str_match(toupper(trimws(allele_vec)), ".*:([A-Z]+)>([A-Z]+)$")
    ref_i <- m_id[, 2]
    alt_i <- m_id[, 3]
    ok_id <- !is.na(ref_i) & !is.na(alt_i) & nzchar(ref_i) & nzchar(alt_i)

    ref <- ifelse(ok_snp, ref_s, ifelse(ok_id, ref_i, NA_character_))
    alt <- ifelse(ok_snp, alt_s, ifelse(ok_id, alt_i, NA_character_))

    list(REF = ref, ALT = alt)
  }

  dart_is_minus <- function(x) {
    z <- tolower(trimws(as.character(x)))
    z[is.na(z)] <- ""
    z == "-" | startsWith(z, "m") | z %in% c("minus", "neg", "negative", "reverse", "rev")
  }

  strand_to_word <- function(x) {
    z <- tolower(trimws(as.character(x)))
    ifelse(is.na(z) | z == "" | z == ".", "none",
      ifelse(z %in% c("+","plus","pos","positive","p"), "plus",
        ifelse(z %in% c("-","minus","neg","negative","m"), "minus", "none")))
  }
  is_blank <- function(x) is.na(x) | x == "" | x == "."

  # ----------------------------
  # MAPPRIORS lookups (optional)
  # ----------------------------
  pri_flag_by_id  <- logical(0)
  dup_label_by_id <- character(0)

  if (is.data.frame(MAPPRIORS)) {
    n  <- nrow(MAPPRIORS)
    tp <- MAPPRIORS
    key <- as.character(tp$qaccver)

    is_true_chr <- function(v) tolower(trimws(as.character(v))) == "true"
    norm_true   <- function(col) if (is.null(col)) rep(FALSE, n) else rep_len(is_true_chr(col), n)
    norm_is_na  <- function(col) if (is.null(col)) rep(TRUE,  n) else rep_len(is.na(col),        n)

    tu_na <- norm_is_na(tp$Trueunique)
    any_prior_true <- norm_true(tp$Chrommarkermap) |
                      norm_true(tp$proximatemarkermap) |
                      norm_true(tp$linkagemarkermap) |
                      norm_true(tp$geneticmapmarkermap)

    safe_set_names <- function(x, nm) {
      lx <- length(x); ln <- length(nm)
      if (!lx || !ln) return(x)
      if (lx == ln)   return(stats::setNames(x, nm))
      k <- min(lx, ln)
      stats::setNames(x[seq_len(k)], nm[seq_len(k)])
    }

    pri_flag_by_id  <- safe_set_names(tu_na & any_prior_true, key)
    dup_is_true     <- norm_true(tp$Duplicate_region)
    dup_label_by_id <- safe_set_names(ifelse(dup_is_true, "LocallyDuplicatedRegion", "SingleCopy"), key)
  }

  # ----------------------------
  # detect DArT layout
  # ----------------------------
  lay <- detect_dart_layout_base(raw_path)
  cn  <- lay$cn

  # hard match AlleleID and SNP (exact), and Strand* (prefix match)
  idx_allele <- match("AlleleID", cn)
  idx_snp    <- match("SNP", cn)  # exact only; may be NA
  idx_strand <- {
    j <- which(grepl("^Strand", cn, ignore.case = TRUE))
    if (length(j)) j[1] else NA_integer_
  }

  if (is.na(idx_allele))
    stop("DArT header must include AlleleID (exact after trimming/canonicalization).")

  g_start  <- lay$lastmetric_idx + 1L
  geno_idx <- g_start:length(cn)

  fileDate  <- format(Sys.Date(), "%Y%m%d")
  VCFchroms <- unique(c(order_contigs(Finalmappings$saccver), "chrUnk"))

  # ----------------------------
  # header scaffold (STANDARDISED)
  # ----------------------------
  header_lines <- c(
    "##fileformat=VCFv4.2",
    sprintf("##fileDate=%s", fileDate),
    "##source=Brioche-VCF-build",
    paste0("##", Metadata[1]),
    paste0("##", Metadata[2]),
    paste0("##", Metadata[3]),
    paste0("##Brioche_filtering_", Metadata[4]),
    paste0("##Brioche_filtering_", Metadata[5]),
    paste0("##Brioche_filtering_", Metadata[6]),
    paste0("##Brioche_priors_files:", Metadata[7]),
    paste0("##Brioche_priors_files:", Metadata[8]),
    paste0("##Brioche_priors_files:", Metadata[9]),
    paste0("##Brioche_priors_files:", Metadata[10]),
    paste0("##Brioche_", Metadata[11]),
    extra_meta,
    '##FILTER=<ID=LowQual,Description="Low quality or ambiguous mapping">',
    '##INFO=<ID=MAPSTATUS,Number=1,Type=String,Description="Unique_mapping|Unique_mapping_with_priors_information|Failed_to_map_uniquely">',
    '##INFO=<ID=PriorORIENT,Number=1,Type=String,Description="Accumulated strand orientation across runs (plus|minus|none)">',
    '##INFO=<ID=DUP,Number=1,Type=String,Description="SingleCopy|LocallyDuplicatedRegion">',
    '##INFO=<ID=MAF,Number=A,Type=Float,Description="Minor Allele frequency per ALT">',
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count per ALT">',
    '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=NU,Number=1,Type=Integer,Description="Null allele flag from SNP chip (1= Null allele; 0=otherwise)">',
    '##FORMAT=<ID=DU,Number=1,Type=Integer,Description="Number of local duplications of marker detected (0= Marker is single copy; 1+=copies, .=Unknown)">',
    contig_lines
  )
  if (nzchar(droplist_path)) {
    header_lines <- c(header_lines, sprintf('##Brioche_droplist="%s"', basename(droplist_path)))
  }
  header_lines <- header_lines[!duplicated(header_lines)]
  writeLines(header_lines, con = Outputfilename)

  cn_fixed   <- make.unique(lay$cn, sep = "...")
  sample_ids <- cn_fixed[geno_idx]
  vcf_header <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", sample_ids)
  write_header_row(vcf_header, Outputfilename)

  # ----------------------------
  # mapping frame (include sstrand)
  # ----------------------------
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

  # ----------------------------
  # stream input
  # ----------------------------
  chunk_size <- 1500L
  delim <- lay$delim

  con_dart <- if (grepl("\\.gz$", raw_path, ignore.case = TRUE)) gzfile(raw_path, "rt") else file(raw_path, "rt")
  on.exit(try(close(con_dart), silent = TRUE), add = TRUE)
  invisible(readLines(con_dart, n = lay$topskip + 1L, warn = FALSE))  # skip preamble+header

  repeat {
    raw_lines <- readLines(con_dart, n = chunk_size, warn = FALSE)
    if (!length(raw_lines)) break

    raw_lines <- raw_lines[nzchar(raw_lines) & !grepl("^\\s*$", raw_lines)]
    if (!length(raw_lines)) next

    df_all <- utils::read.table(
      text = raw_lines, sep = delim, header = FALSE,
      quote = "", comment.char = "", stringsAsFactors = FALSE,
      check.names = FALSE, col.names = cn_fixed, fill = TRUE
    )
    if (!NROW(df_all)) next

    Name <- df_all[[idx_allele]]  # qaccver / AlleleID
    force_drop <- if (length(DROP_IDS)) Name %in_drop% DROP_IDS else rep(FALSE, length(Name))

    SNP <- if (!is.na(idx_snp)) df_all[[idx_snp]] else rep(NA_character_, nrow(df_all))

    # DArT "Strand*" column (optional)
    dart_strand <- if (!is.na(idx_strand)) df_all[[idx_strand]] else rep(NA_character_, nrow(df_all))

    geno <- as.matrix(df_all[, geno_idx, drop = FALSE])

    # ---- NEW: parse REF/ALT from SNP (position:REF>ALT) OR AlleleID suffix :REF>ALT$ ----
    pa      <- dart_extract_ref_alt(SNP, Name)
    REF_old <- toupper(pa$REF)
    ALT_old <- toupper(pa$ALT)

    # ---- NEW: if Strand is Minus, reverse complement REF/ALT before comparing/reanchoring ----
    if (!is.na(idx_strand)) {
      is_min <- dart_is_minus(dart_strand)
      if (any(is_min, na.rm = TRUE)) {
        REF_old[is_min] <- revcomp_iupac(REF_old[is_min])
        ALT_old[is_min] <- revcomp_iupac(ALT_old[is_min])
      }
    }

    # enforce single ALT
    ALT_old <- sub(",.*$", "", ALT_old)
    ALT_old[is.na(ALT_old) | ALT_old == ""] <- "."

    idx   <- match(Name, fm_Name)
    has_m <- !is.na(idx)
    hits  <- ifelse(has_m, fm$MAP_HITS[idx], 0L)
    u_map <- has_m & hits == 1L & !force_drop

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
    ALT_out <- sub(",.*$", "", ALT_out)
    ALT_out[is.na(ALT_out) | ALT_out == ""] <- "."

    pos_int     <- suppressWarnings(as.integer(POS))
    coord_bad   <- is.na(pos_int) | pos_int <= 0L | is_blank(CHROM)
    unknown_now <- !u_map | coord_bad

    QUAL   <- ifelse(unknown_now, 0L, 100L)
    FILTER <- ifelse(unknown_now, "LowQual", "PASS")

    # PriorORIENT from Brioche sstrand when uniquely mapped; else 'none'
    cur_strand <- ifelse(u_map, strand_to_word(fm$STRAND[idx]), "none")

    status_val <- ifelse(unknown_now, "Failed_to_map_uniquely", "Unique_mapping")
    if (length(pri_flag_by_id)) {
      pri_flag <- unname(pri_flag_by_id[Name])
      pri_flag[is.na(pri_flag)] <- FALSE
      status_val[u_map & pri_flag] <- "Unique_mapping_with_priors_information"
    }

    dup_lab <- if (length(dup_label_by_id)) {
      out <- unname(dup_label_by_id[Name]); out[is.na(out)] <- "SingleCopy"; out
    } else {
      rep("SingleCopy", length(Name))
    }

    INFO <- paste0(
      "MAPSTATUS=", status_val,
      ";PriorORIENT=", cur_strand,
      ";DUP=", dup_lab
    )

    # Genotype normalization & swap where needed
    geno[] <- toupper(gsub("\\|", "/", geno, perl = TRUE))
    miss <- is.na(geno) | geno == "" | geno == "." | geno == "-" | geno == "N" | geno == "NC"
    if (any(miss)) geno[miss] <- "./."
    dbl <- grepl("^[012]{2}$", geno, perl = TRUE)
    if (any(dbl)) geno[dbl] <- sub("^([012])([012])$", "\\1/\\2", geno[dbl], perl = TRUE)

    rswap <- ifelse(is.na(swap | swap_rc), FALSE, (swap | swap_rc))
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

    # AC/AN/MAF (biallelic)
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

    # STANDARDISED FORMAT: GT:NU:DU
    # DArT has no chip-derived null allele flag => NU=0 always; DU unknown => '.'
    FORMAT <- "GT:NU:DU"
    NU_chr <- matrix("0", nrow = nrow(geno), ncol = ncol(geno))
    DU_chr <- matrix(".", nrow = nrow(geno), ncol = ncol(geno))

    geno_df <- as.data.frame(geno, stringsAsFactors = FALSE, check.names = FALSE)
    for (jj in seq_len(ncol(geno_df))) {
      geno_df[[jj]] <- paste0(geno_df[[jj]], ":", NU_chr[, jj], ":", DU_chr[, jj])
    }

    body <- data.frame(
      `#CHROM` = CHROM, POS = POS, ID = Name,
      REF = REF_out, ALT = ALT_out, QUAL = QUAL, FILTER = FILTER, INFO = INFO, FORMAT = FORMAT,
      stringsAsFactors = FALSE, check.names = FALSE
    )
    body <- cbind(body, geno_df, stringsAsFactors = FALSE)

    dt_fwrite(body, file = Outputfilename, append = TRUE, col.names = FALSE)

    rm(df_all, geno, geno_df, body, ac1, an, NU_chr, DU_chr); gc(FALSE)
  }

  message(sprintf("Wrote: %s", Outputfilename))
  quit(save = "no")
}

# Come back to this section to add in more diverse datatypes for non SNP chip output. 
if(!is_snpchip & !is_vcfraws & !isDArTfile) {
  message(sprintf("Wrote: %s", "Hello world"))
}
