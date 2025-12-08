#!/usr/bin/env Rscript

# Usage:
# Rscript DArT2Brioche_convert.R --dartfile "dart_testing_build_briochefiles.csv" --targetmarker "D"

.t0 <- Sys.time()

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


parse_target_base <- function(snp) {
  m <- stri_match_first_regex(snp, "([ACGTacgt]).*?([ACGTacgt])")
  a <- toupper(m[, 2])
  b <- toupper(m[, 3])
  out <- ifelse(!is.na(a) & !is.na(b), paste0("[", a, "/", b, "]"), NA_character_)
  out
}


replace_at <- function(strings, pos, repl) {
  out <- strings
  ok  <- !is.na(out) & !is.na(pos) & pos >= 1

  lens <- nchar(out, type = "chars", allowNA = TRUE, keepNA = TRUE)
  ok[lens < pos & !is.na(lens)] <- FALSE
  if (any(!ok)) {
    # (silent skip for invalid positions; uncomment to warn)
    # warning("Some positions invalid or exceed sequence length; leaving those unchanged.")
  }
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


hit <- Reduce(`|`, lapply(dartdata, function(col) as.character(col) == "AlleleID"))
i_header <- which(hit)[1]
if (is.na(i_header)) stop("'AlleleID' not found in any column.", call. = FALSE)

df <- dartdata[i_header:nrow(dartdata), , drop = FALSE]
rm(dartdata); gc()


names(df) <- df[1, ]
df <- df[-1, , drop = FALSE]

# Remove excess columns for speed up
idx <- match("CallRate", names(df))
if (is.na(idx)) stop("'CallRate' column not found after header normalization.", call. = FALSE)
df <- df[, seq_len(idx), drop = FALSE]

required_cols <- c("AlleleID", "SNP", "SnpPosition", "AlleleSequence")
stop_if_missing(df, required_cols)
df <- df[, required_cols, drop = FALSE]
gc()


df <- df[!is.na(df$SNP) & nzchar(df$SNP), , drop = FALSE]
df$AlleleID        <- as.character(df$AlleleID)
df$AlleleSequence  <- as.character(df$AlleleSequence)
df$SnpPosition     <- suppressWarnings(as.integer(df$SnpPosition))

# build output
Target.bp   <- df$SnpPosition + 1L
Target.base <- parse_target_base(df$SNP)
Sequence    <- replace_at(df$AlleleSequence, Target.bp, targetmarker)

Brioche.out <- data.frame(
  ID          = df$AlleleID,
  Sequence    = Sequence,
  Target.bp   = Target.bp,
  Target.base = Target.base,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# df is no longer needed
rm(df); gc()
#write output
outfile <- sub("(\\.[^./\\\\]+)?$", "_brioche.tsv", dartfile)
write.table(Brioche.out, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("Wrote %d rows to %s", nrow(Brioche.out), outfile))




.dt <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
message(sprintf("Total time: %.2f seconds", .dt))




