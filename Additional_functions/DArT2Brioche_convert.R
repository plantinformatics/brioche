#!/usr/bin/env Rscript

# Usage:
# Rscript DArT2Brioche_convert.R --dartfile "dart_testing_build_briochefiles.csv" --targetmarker "D"
# Converts 1-row-per-marker DArT reports to Brioche input (ID, Sequence, Target.bp, Target.base).
# Notes:
#   * Column names in DArT exports can vary. This script:
#       - Locates the header row by finding the first row that contains "AlleleID".
#       - Picks the sequence column using the ordered fallback:
#           AlleleSequence -> AlleleSequenceRef -> AlleleSequenceSnp ->
#           TrimmedSequence -> TrimmedSequenceRef -> TrimmedSequenceSnp
#   * If "CallRate" is present, columns are trimmed up to CallRate for speed; otherwise no trimming.

.t0 <- Sys.time()

# Quiet-install stringi if missing (used for fast substring ops and regex)
if (!requireNamespace("stringi", quietly = TRUE)) {
  suppressWarnings(suppressMessages(
    install.packages("stringi", repos = "https://cloud.r-project.org", quiet = TRUE)
  ))
}

suppressPackageStartupMessages({
  library(stringi)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "$"), args)
  if (length(hit) == 1 && hit < length(args)) return(args[hit + 1])
  # also allow --flag=value
  hit2 <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit2) == 1) return(sub(paste0("^", flag, "="), "", hit2))
  if (!is.null(default)) return(default)
  stop(sprintf("Missing required argument: %s", flag), call. = FALSE)
}

dartfile     <- get_arg("--dartfile")
targetmarker <- get_arg("--targetmarker")

stop_if_missing <- function(df, cols) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) {
    stop(sprintf("Required column(s) not found: %s", paste(miss, collapse = ", ")), call. = FALSE)
  }
}

# Parse [A/C] from "SNP" column (or similar two-base annotation)
parse_target_base <- function(snp) {
  m <- stri_match_first_regex(snp, "([ACGTacgt]).*?([ACGTacgt])")
  a <- toupper(m[, 2])
  b <- toupper(m[, 3])
  out <- ifelse(!is.na(a) & !is.na(b), paste0("[", a, "/", b, "]"), NA_character_)
  out
}

# Replace a single character at 1-based position 'pos'
replace_at <- function(strings, pos, repl) {
  out <- strings
  ok  <- !is.na(out) & !is.na(pos) & pos >= 1
  lens <- nchar(out, type = "chars", allowNA = TRUE, keepNA = TRUE)
  ok[lens < pos & !is.na(lens)] <- FALSE
  if (length(repl) == 1L) repl <- rep_len(repl, length(out))
  idx <- which(ok)
  if (length(idx)) {
    stri_sub(out[idx], from = pos[idx], to = pos[idx]) <- repl[idx]
  }
  out
}

dartdata <- tryCatch(
  read.csv(dartfile, header = FALSE, stringsAsFactors = FALSE, fill = TRUE, check.names = FALSE),
  error = function(e) stop(sprintf("Failed to read file '%s': %s", dartfile, e$message), call. = FALSE)
)

# Find first row containing literal "AlleleID" in any column; treat that as header row
hit <- Reduce(`|`, lapply(dartdata, function(col) as.character(col) == "AlleleID"))
i_header <- which(hit)[1]
if (is.na(i_header)) stop("'AlleleID' not found in any column.", call. = FALSE)

df <- dartdata[i_header:nrow(dartdata), , drop = FALSE]
rm(dartdata); gc()

names(df) <- df[1, ]
df <- df[-1, , drop = FALSE]

# Optional column trimming: only if CallRate exists (for speed). Otherwise keep all.
idx_callrate <- match("CallRate", names(df))
if (!is.na(idx_callrate)) {
  df <- df[, seq_len(idx_callrate), drop = FALSE]
}

seq_candidates <- c(
  "AlleleSequence",
  "AlleleSequenceRef",
  "AlleleSequenceSnp",
  "TrimmedSequence",
  "TrimmedSequenceRef",
  "TrimmedSequenceSnp"
)

seq_col <- NULL
for (nm in seq_candidates) {
  if (nm %in% names(df)) { seq_col <- nm; break }
}
if (is.null(seq_col)) {
  stop(sprintf(
    "None of the expected sequence columns were found. Tried: %s",
    paste(seq_candidates, collapse = ", ")
  ), call. = FALSE)
}
message(sprintf("[INFO] Using sequence column: %s", seq_col))

required_cols <- c("AlleleID", "SNP", "SnpPosition", seq_col)
stop_if_missing(df, required_cols)

# Keep only the needed columns (order them explicitly)
df <- df[, required_cols, drop = FALSE]

df$AlleleID    <- as.character(df$AlleleID)
df$SNP         <- as.character(df$SNP)
df[[seq_col]]  <- as.character(df[[seq_col]])
df$SnpPosition <- suppressWarnings(as.integer(df$SnpPosition))

# Filter out rows missing SNP info
df <- df[!is.na(df$SNP) & nzchar(df$SNP), , drop = FALSE]

Target.bp   <- df$SnpPosition + 1L
Target.base <- parse_target_base(df$SNP)
Sequence    <- replace_at(df[[seq_col]], Target.bp, targetmarker)

Brioche.out <- data.frame(
  ID          = df$AlleleID,
  Sequence    = Sequence,
  Target.bp   = Target.bp,
  Target.base = Target.base,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Free memory
rm(df); gc()

outfile <- sub("(\\.[^./\\\\]+)?$", "_brioche.tsv", dartfile)
write.table(Brioche.out, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %d rows to %s", nrow(Brioche.out), outfile))

.dt <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
message(sprintf("Total time: %.2f seconds", .dt))
