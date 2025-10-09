#!/usr/bin/env Rscript

# summarise_results.R
# Renders a tabbed HTML summary for a Nextflow run.
# Usage:
#   Rscript summarise_results.R <meta_json> <results_dir> <out_html>
#       [--timeline FILE] [--report FILE] [--dag FILE]
#       [--fig_dir DIR] [--rmd_out FILE]
#
# Notes:
# - Requires packages: jsonlite, rmarkdown, (optional) htmltools

# ------------------------ arg parsing ------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste(
    "Usage:",
    "Rscript summarise_results.R <meta_json> <results_dir> <out_html>",
    "[--timeline FILE] [--report FILE] [--dag FILE]",
    "[--fig_dir DIR] [--rmd_out FILE]",
    sep = "\n  "
  ))
}

.meta_json <- args[1]
.results   <- args[2]
.out_html  <- args[3]

.timeline <- NA_character_
.report   <- NA_character_
.dag      <- NA_character_
.fig_dir  <- NA_character_
.rmd_out  <- NA_character_
.params_cfg <- NA_character_

if (length(args) > 3) {
  i <- 4
  while (i <= length(args)) {
    flag <- args[i]
    if (flag == "--timeline") { .timeline <- args[i+1]; i <- i + 2; next }
    if (flag == "--report")   { .report   <- args[i+1]; i <- i + 2; next }
    if (flag == "--dag")      { .dag      <- args[i+1]; i <- i + 2; next }
    if (flag == "--fig_dir")  { .fig_dir  <- args[i+1]; i <- i + 2; next }
    if (flag == "--rmd_out")  { .rmd_out  <- args[i+1]; i <- i + 2; next }
    if (flag == "--params_cfg") { .params_cfg <- args[i+1]; i <- i + 2; next }
    stop(sprintf("Unknown flag: %s", flag))
  }
}

# ------------------------ deps check ------------------------
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required. Install it (e.g., conda/mamba: r-jsonlite).")
}
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required. Install it (e.g., conda/mamba: r-rmarkdown).")
}


.rel_to_out <- function(path, out_html) {
  if (is.na(path) || !nzchar(path)) return("")
  out_dir <- normalizePath(dirname(out_html), winslash = "/", mustWork = FALSE)
  target  <- tryCatch(normalizePath(path, winslash = "/", mustWork = TRUE),
                      error = function(e) normalizePath(path, winslash = "/", mustWork = FALSE))
  out_parts <- strsplit(out_dir, "/")[[1]]
  tgt_parts <- strsplit(target,  "/")[[1]]
  maxlen <- min(length(out_parts), length(tgt_parts))
  k <- 0L
  while (k < maxlen && out_parts[k + 1L] == tgt_parts[k + 1L]) k <- k + 1L
  up   <- if ((length(out_parts) - k) > 0) rep("..", length(out_parts) - k) else character(0)
  down <- if (k < length(tgt_parts)) tgt_parts[(k + 1L):length(tgt_parts)] else character(0)
  rel  <- paste(c(up, down), collapse = "/")
  if (!nzchar(rel)) rel <- "."
  rel
}
`%||%` <- function(a, b) if (is.null(a) || is.na(a) || !nzchar(a)) b else a

  save.image("test000.Rdata") 
.describe_file <- function(p) {
  bn <- basename(p)
  if (grepl("summary_core\\.html$", bn, TRUE)) return("Core pipeline summary (HTML)")
  if (grepl("summary_report\\.html$", bn, TRUE)) return("Final pipeline summary (HTML)")
  if (grepl("_timeline\\.html$", bn, TRUE)) return("Nextflow execution timeline")
  if (grepl("_report\\.html$", bn, TRUE)) return("Nextflow resources/tasks report")
  if (grepl("_dag\\.", bn, TRUE)) return("Nextflow DAG")
  if (grepl("\\.(tsv|csv)(\\.gz)?$", bn, TRUE)) return("Tabular results")
  if (grepl("\\.(fa|fna|fasta|fastq|fq)(\\.gz)?$", bn, TRUE)) return("Sequence data")
  if (grepl("\\.(bam|cram|bai)$", bn, TRUE)) return("Alignment artefact")
  if (grepl("\\.(png|svg|jpg|jpeg)$", bn, TRUE)) return("Figure/image")
  if (grepl("\\.(rds|rdata)$", bn, TRUE)) return("R data object")
  "Output file"
}

.safe_exists <- function(x) {
  tryCatch(!is.na(x) && nzchar(x) && file.exists(x), error = function(e) FALSE)
}

.escape_md <- function(x) {
  x <- gsub("\\|", "\\\\|", x)
  x
}

.make_include_chunk <- function(path, out_html, height = 900) {
  rel <- .rel_to_out(path, out_html)
  c(
    "```{r, echo=FALSE, results='asis'}",
    "if (requireNamespace('htmltools', quietly = TRUE)) {",
    sprintf("  htmltools::includeHTML('%s')", rel),
    "} else {",
    sprintf("  cat('<iframe src=\"%s\" width=\"100%%\" height=\"%d\" style=\"border:1px solid #ddd;\"></iframe>')", rel, height),
    "}",
    "```"
  )
}

# A generic placeholder paragraph (easy to Ctrl+F and edit later)
.placeholder_para <- function(tag, name = "") {
  label <- if (nzchar(name)) sprintf(" (%s)", name) else ""
  sprintf("_**Description%s:** <!-- %s -->_", label, tag)
}



  save.image("test00.Rdata") 

# metadata
if (!file.exists(.meta_json)) {
  stop(sprintf("meta_json does not exist: %s", .meta_json))
}
meta <- jsonlite::fromJSON(.meta_json)
meta_df <- data.frame(
  Field = names(meta),
  Value = vapply(meta, function(x) paste0(x, collapse = ", "), character(1)),
  stringsAsFactors = FALSE
)

params_cfg_path <- .params_cfg
if (is.na(params_cfg_path) || !file.exists(params_cfg_path)) {
  cand <- c(
    file.path(meta$launchDir %||% "", "params.config"),
    file.path(meta$projectDir %||% "", "params.config"),
    file.path(.results,                    "params.config")
  )
  cand <- unique(cand[file.exists(cand)])
  if (length(cand)) params_cfg_path <- cand[1] else params_cfg_path <- NA_character_
}

params_cfg_text <- NULL
if (!is.na(params_cfg_path) && file.exists(params_cfg_path)) {
  params_cfg_text <- paste0(readLines(params_cfg_path, warn = FALSE), collapse = "\n")
}

# figs (ordered: marker_filtering_stages, BLASTnhits_filtering_stage, histogram_hits, then others)
fig_paths <- character(0)
if (!is.na(.fig_dir) && nzchar(.fig_dir) && dir.exists(.fig_dir)) {
  fig_paths <- list.files(
    .fig_dir,
    pattern     = "\\.(png|svg|jpg|jpeg)$",
    full.names  = TRUE,
    ignore.case = TRUE
  )
  
  if (length(fig_paths)) {
    bn <- basename(fig_paths)
    
    prio <- rep(99L, length(bn))  # default lowest priority (others)
    prio[grepl("marker_filtering_stages",    bn, ignore.case = TRUE)] <- 1
    prio[grepl("BLASTnhits_filtering_stage", bn, ignore.case = TRUE)] <- 2
    prio[grepl("histogram_hits",             bn, ignore.case = TRUE)] <- 3
    
    ord <- order(prio, bn)  # stable within group by filename
    fig_paths <- fig_paths[ord]
  }
}

all_files <- list.files(
  .results,
  pattern     = "\\.(tsv|csv|xlsx)$",
  all.files   = FALSE,
  full.names  = TRUE,
  recursive   = TRUE,
  ignore.case = TRUE,
  no..        = TRUE
)

if (length(all_files)) {
  bn   <- basename(all_files)
  prio <- rep(999L, length(bn)) 
  

  prio[grepl("unfiltered",                           bn, ignore.case = TRUE)] <- 1
  prio[grepl("all_mappings",                         bn, ignore.case = TRUE)] <- 2
  prio[grepl("filtered_mappings\\.csv$",            bn, ignore.case = TRUE)] <- 3
  prio[grepl("mappings\\.tsv$",                     bn, ignore.case = TRUE)] <- 4
  prio[grepl("intermediate_filtering_hits\\.csv$",  bn, ignore.case = TRUE)] <- 5
  prio[grepl("intermediate_filtering_mappings\\.tsv$", bn, ignore.case = TRUE)] <- 6
  prio[grepl("mappings-pretzel",                     bn, ignore.case = TRUE)] <- 7
  prio[grepl("strict_filtering_hits",                bn, ignore.case = TRUE)] <- 8  
  prio[grepl("strict_filtering_mappings\\.tsv$",    bn, ignore.case = TRUE)] <- 9
  prio[grepl("summary_filtering\\.csv$",            bn, ignore.case = TRUE)] <- 10
  
  ord <- order(prio, bn)
  all_files <- all_files[ord]


}
  save.image("test0.Rdata") 



this_out_abs <- normalizePath(.out_html, winslash = "/", mustWork = FALSE)
all_files <- all_files[normalizePath(all_files, winslash = "/", mustWork = FALSE) != this_out_abs]

sizes <- tryCatch(file.info(all_files)$size, error = function(e) rep(NA_real_, length(all_files)))
desc  <- vapply(all_files, .describe_file, character(1))
relp  <- vapply(all_files, function(p) .rel_to_out(p, .out_html), character(1))

outputs_df <- data.frame(
  File = relp,
  Size_Bytes = sizes,
  Description = desc,
  stringsAsFactors = FALSE
)
save.image("test1.Rdata")
# ------------------------ build Rmd ------------------------
if (is.na(.rmd_out) || !nzchar(.rmd_out)) {
  .rmd_out <- sub("\\.html?$", ".Rmd", .out_html, TRUE)
}
dir.create(dirname(.rmd_out), showWarnings = FALSE, recursive = TRUE)

rmd <- character()

rmd <- c(
  rmd,
  "---",
  "title: \"Pipeline Summary\"",
  "output:",
  "  html_document:",
  "    toc: true",
  "    toc_depth: 3",
  "    theme: readable",
  "    df_print: paged",
  "    self_contained: true",
  "---",
  "",
  "# Summary {.tabset}",
  "",
  "## Overview",
  "",
  "### Run metadata",
  "",
  "| Field | Value |",
  "|------:|:------|"
)

# metadata rows
if (nrow(meta_df)) {
  for (i in seq_len(nrow(meta_df))) {
    rmd <- c(rmd, sprintf("| %s | %s |",
                          .escape_md(meta_df$Field[i]),
                          .escape_md(meta_df$Value[i])))
  }
} else {
  rmd <- c(rmd, "_No metadata available._")
}

rmd <- c(rmd, "", "### Parameters config (source)")
if (!is.null(params_cfg_text)) {
  rmd <- c(
    rmd,
    sprintf("_Detected at_: `%s`", .rel_to_out(params_cfg_path, .out_html)),
    "",
    "```conf",
    params_cfg_text,
    "```"
  )
} else {
  rmd <- c(rmd, "_No params.config file was found at launchDir/projectDir/results._")
}

# ---------- Figures section ----------
rmd <- c(rmd, "", "### Generated figures")
if (length(fig_paths)) {
  for (fp in fig_paths) {
    link <- .rel_to_out(fp, .out_html)
    nm   <- basename(fp)
    rmd <- c(
      rmd,
      "",
      sprintf("**%s**", .escape_md(nm)),
      # >>> placeholder BEFORE the figure <<<
      .placeholder_para(tag = paste0("FIGDESC:", nm), name = nm),
      "",
      sprintf("![](%s)", link)
    )
  }
} else {
  rmd <- c(rmd, "", "_No figures were found in the provided directory._")
}

# ---------- Outputs table ----------
rmd <- c(rmd, "", "### Output files (click to open)", "",
         "| File | Size (bytes) | Description |",
         "|:-----|------------:|:------------|")
if (nrow(outputs_df)) {
  for (i in seq_len(nrow(outputs_df))) {
    link <- sprintf("[`%s`](%s)", .escape_md(outputs_df$File[i]), outputs_df$File[i])
    size <- ifelse(is.na(outputs_df$Size_Bytes[i]), "",
                   format(outputs_df$Size_Bytes[i], big.mark = ",", scientific = FALSE))
    rmd <- c(rmd, sprintf("| %s | %s | %s |", link, size, .escape_md(outputs_df$Description[i])))
  }
} else {
  rmd <- c(rmd, "| _No files found_ |  |  |")
}

# Optional: per-file notes (placeholders) right after the table
if (nrow(outputs_df)) {
  rmd <- c(rmd, "", "#### Notes for output files")
  for (i in seq_len(nrow(outputs_df))) {
    nm <- basename(outputs_df$File[i])
    rmd <- c(
      rmd,
      "",
      sprintf("**%s**", .escape_md(nm)),
      # >>> placeholder BEFORE the link description <<<
      .placeholder_para(tag = paste0("FILEDESC:", nm), name = nm)
    )
  }
}

# ---------- Optional artefacts (Timeline/Report/DAG) ----------
if (.safe_exists(.timeline) || .safe_exists(.report) || .safe_exists(.dag)) {
  rmd <- c(rmd, "", "## Additional run artefacts")
  
  if (.safe_exists(.timeline)) {
    rmd <- c(
      rmd, "",
      "### Timeline",
      # >>> placeholder BEFORE the embed <<<
      .placeholder_para(tag = "TIMELINE_DESC", name = basename(.timeline)),
      .make_include_chunk(.timeline, .out_html, height = 900),
      "",
      sprintf("> Open directly: [%s](%s)", basename(.timeline), .rel_to_out(.timeline, .out_html))
    )
  }
  
  if (.safe_exists(.report)) {
    rmd <- c(
      rmd, "",
      "### Report",
      # >>> placeholder BEFORE the embed <<<
      .placeholder_para(tag = "REPORT_DESC", name = basename(.report)),
      .make_include_chunk(.report, .out_html, height = 900),
      "",
      sprintf("> Open directly: [%s](%s)", basename(.report), .rel_to_out(.report, .out_html))
    )
  }
  
  if (.safe_exists(.dag)) {
    rmd <- c(rmd, "", "### DAG")
    # >>> placeholder BEFORE the DAG embed <<<
    rmd <- c(rmd, .placeholder_para(tag = "DAG_DESC", name = basename(.dag)))
    if (grepl("\\.html?$", .dag, TRUE)) {
      rmd <- c(
        rmd,
        .make_include_chunk(.dag, .out_html, height = 900),
        "",
        sprintf("> Open directly: [%s](%s)", basename(.dag), .rel_to_out(.dag, .out_html))
      )
    } else if (grepl("\\.(svg|png|jpg|jpeg)$", .dag, TRUE)) {
      rmd <- c(rmd, sprintf("![](%s)", .rel_to_out(.dag, .out_html)))
    } else if (grepl("\\.pdf$", .dag, TRUE)) {
      rmd <- c(rmd, sprintf("<embed src=\"%s\" type=\"application/pdf\" width=\"100%%\" height=\"900px\">",
                            .rel_to_out(.dag, .out_html)))
    } else {
      rmd <- c(rmd, "_DAG file found but format not supported for embedding._")
    }
  }
}


save.image("test2.Rdata")
fig_desc <- c(
  "marker_filtering_stages.png"    = "\n Shows % PASS/FAIL for markers  \n at each filtering stage.",
  "BLASTnhits_filtering_stages.png" = "\n Proportion of BLASTn hits retained and dropped at each filtering stage.",
  "histogram_hits.png"             = "\n Frequency distribution of BLAST hit counts \n per marker at each filtering stage."
)

file_desc <- c(
  "unfiltered"                            = "\n Blast hits file with identified SNP locations prior to any filtering to reduce blast hits. This contains every blast hit returned for each marker and will be very large",
  "all_mappings"                          = "\n Blast hits file with identified SNP locations filtered to remove hits which failed the user given pidentity and coverage filters (subset of unfiltered)",
  "filtered_mappings.csv"                 = "\n Blast hits file with identified SNP locations filtered to remove hits which which failed hybridisation filter (if applied) and filtered to remove all markers where there were still 10 or more hits to different locations on the reference genome (subset of all_mappings).",
  "mappings.tsv"                          = "\n Tab-delimited mappings file (in short format, SNP_ID,Unique_SNPID, SNP_position,REF_nt,ALT,Hybridised,Coverage,pident) for all hits which passed the hybridisation/secondary hits filter",
  "intermediate_filtering_hits.csv"       = "\n Blast hits file with identified SNP locations filtered under intermediate filtering stage where if no prior information was provided (chromosome locations/linked markers/ld maps) all hits which showed the highest bitscore per chromosome are kept if priors were given, all hits for the known target chromosome were kept and all hits of equal bitscore to the maximum bitscore of a hit for the known target chromosome are kept. (subset of filtered_mappings.csv)",
  "intermediate_filtering_mappings.tsv"   = "\n Tab-delimited mappings file (in short format, SNP_ID,Unique_SNPID, SNP_position,REF_nt,ALT,Hybridised,Coverage,pident) for all hits which passed the intermediate filtering stage.",
  "mappings-pretzel"                      = "\n Mappings formatted for Pretzel visualization. This file contains the markers/hits filtered to the intermediate_filtering stage.",
  "strict_filtering_hits"                 = "\n Blast hits file with identified SNP locations filtered under strict filtering stage.where if no prior information was provided (chromosome locations/linked markers/ld maps) any marker with 2 or more hits is removed. if priors were given the following process was followed \n A) target chromosome given: if only one hit is on the target chromosome, all other hits are removed and that hit is kept. If multiple/zero hits to target chromosome for the given marker, all hits for the marker are removed \n B) If linked markers/ld maps were provided the exact position of all linked markers is used to identify the closest hit to linked markers and all other hits are removed (subset of intermediate_filtering_hits.csv)",
  "strict_filtering_mappings.tsv"         = "\n Tab-delimited mappings file (in short format, SNP_ID,Unique_SNPID, SNP_position,REF_nt,ALT,Hybridised,Coverage,pident) for all hits which passed the strict filtering stage.",
  "summary_filtering.csv"                 = "\n Per-marker summary across all filters showing how many hits per marker were present at each stage and at what stage all hits for a marker were dropped."
)

art_desc <- c(
  "TIMELINE_DESC" = "Interactive runtime timeline of tasks and resources.",
  "REPORT_DESC"   = "Resource & performance report for each process.",
  "DAG_DESC"      = "Directed acyclic graph of process dependencies."
)

if (length(fig_paths)) {
  for (nm in names(fig_desc)) {
    tag <- paste0("FIGDESC:.*", nm)
    idx <- grep(tag,x=rmd)
    rmd[idx] <- gsub("Description.*\\)",paste0("Description: ",fig_desc[nm]),x=rmd[idx]) 
    } 
}
if (exists("outputs_df") && nrow(outputs_df)) {
  for (nm in names(file_desc)) {
    tag <- paste0("FILEDESC:.*", nm)
    idx <- grep(tag,x=rmd)
    rmd[idx] <- gsub("Description.*\\)",paste0("Description: ",file_desc[nm]),x=rmd[idx]) 

    } 
}

rmd <- gsub("TIMELINE_DESC:.*",
            art_desc[["TIMELINE_DESC"]], rmd)
rmd <- gsub("REPORT_DESC:.*",
            art_desc[["REPORT_DESC"]],   rmd)
rmd <- gsub("DAG_DESC:.*",
            art_desc[["DAG_DESC"]],      rmd)

save.image("test3.Rdata") 

# ------------------------ write Rmd ------------------------
dir.create(dirname(.rmd_out), showWarnings = FALSE, recursive = TRUE)
cat(paste0(rmd, collapse = "\n"), file = .rmd_out)

# ------------------------ render HTML ------------------------
res <- tryCatch({
  rmarkdown::render(
    input       = .rmd_out,
    output_file = basename(.out_html),
    output_dir  = dirname(.out_html),
    quiet       = TRUE
  )
  TRUE
}, error = function(e) {
  message("rmarkdown::render failed: ", conditionMessage(e))
  FALSE
})

if (!res) {
  stop("Failed to render summary HTML. Check that 'rmarkdown' is installed and writable paths exist.")
}

message("Summary written to: ", .out_html)
