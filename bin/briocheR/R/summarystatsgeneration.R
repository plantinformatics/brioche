# Script to collate final results and generate summary statistics for filtering processes 

#' @title Document summary stats 
#' @description This function is used to collate results across all filtering steps and generate summary graphs/csv files for easier interpretation
#' @author James ODwyer
#' @param probe.name Character. Name of the probe set (used in file names).
#' @param genome.name Character. Name of the reference genome (used in file names).
#' @param markers_priorsinformed.hits Table. Table tracking whether each marker was informed through priors and which priors did so
#' @param duplicationmap.hits Table. Table tracking the number of duplications per marker which were present in the reference genome mapped to
#' @param output.path Character. Directory where results are read/written (a `Figs/` subdir is created).
#' @param blast.hits Character. Path to strict-filtered BLAST hits CSV. \emph{Currently not used} (inputs are read from `output.path` by convention).
#' @param mappings.file Character. Path to strict-filtered mappings TSV. \emph{Currently not used}.
#' @param target.file Character. Path to the target-design file (CSV/TSV); the first column contains marker IDs.
#' @return No return value; called for side effects. Writes:
#' \itemize{
#'   \item \code{<probe>_with_<genome>_summary_filtering.csv} (in \code{output.path})
#'   \item \code{Figs/<probe>_with_<genome>_marker_filtering_stages.png}
#'   \item \code{Figs/<probe>_with_<genome>_BLASTnhits_filtering_stages.png}
#'   \item \code{Figs/<probe>_with_<genome>_histogram_hits.png}
#' }
#' @keywords Summary, Data visualisation, Prereport
#' @export
Dosumstats <-
  function(probe.name,
           genome.name,
           markers_priorsinformed.hits,
           duplicationmap.hits,
           output.path,
           blast.hits,
           mappings.file,
           target.file,
           metadata = NULL) {
    if (!dir.exists(paste0(output.path,"/Figs")))
      dir.create(paste0(output.path,"/Figs"), recursive = TRUE)
    
    file <- paste0(output.path,"/Figs/",probe.name, "_with_", genome.name)
    
    ext <- tolower(tools::file_ext(target.file))
    if (ext %in% c("tsv", "tab", "txt")) {
      target_df <- utils::read.delim(target.file, check.names = FALSE, stringsAsFactors = FALSE)
    } else {
      target_df <- utils::read.csv(target.file, check.names = FALSE, stringsAsFactors = FALSE)
    }
    
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  normalize_metadata <- function(md) {
    if (is.null(md)) md <- list()
    list(
      brioche_version   = md$brioche_version   %||% "unknown@unknown",
      brioche_repo_url   = md$brioche_repo_url   %||% "unknown@unknown",
      brioche_commit_url   = md$brioche_commit_url   %||% "unknown@unknown",
      coverage          = md$coverage          %||% NA,
      pident            = md$pident            %||% NA,
      otherblastoptions = md$otherblastoptions %||% "",
      usetargetchrom    = md$usetargetchrom    %||% "No",
      usesharedmarkersmap = md$usesharedmarkersmap %||% "No",
      useldedgemap      = md$useldedgemap      %||% "No",
      usegeneticmap     = md$usegeneticmap     %||% "No",
      run_date          = as.character(md$run_date %||% Sys.time())
    )
  }

  md <- normalize_metadata(metadata)

  md$brioche_version <- trimws(sub("^\\s*local@", "", as.character(md$brioche_version)))

    markersall <- unique(stringr::str_trim(target_df[[1]]))
    
    sumstatsfile <- as.data.frame(matrix(nrow = length(markersall), ncol = 10))
    colnames(sumstatsfile) <- c("Markername", "Blast_hits", "Filter_coverage-identity", "Cov-ident_BLASThits",
                                "Filter_secondaries_and_hybridisation", "Secondaries_hybridisation_BLASThits",
                                "Filter_intermediate_sensitivity", "Intermediate_BLASThits",
                                "Filter_strict_sensitivity", "Strict_BLASThits")
    
    sumstatsfile$Markername <- markersall
    
    # Unfiltered BLAST hit counts
    unf_file <- file.path(output.path,
                          paste0(probe.name, "_with_", genome.name, "_blastmappings_unfiltered.csv"))
    
    counts <- data.table::fread(
      cmd          = sprintf("grep -v '^#' %s", shQuote(unf_file)),
      select       = "qaccver",
      showProgress = TRUE
    ) |>
      dplyr::filter(!is.na(qaccver), qaccver != "") |>
      dplyr::count(qaccver, name = "Blast_hits")    

    # remember original position of Blast_hits
    .orig_names <- names(sumstatsfile)
    .pos_bh     <- match("Blast_hits", .orig_names)
    
    # join and fill zeros (zeros reflective of the marker being dropped completely)
    sumstatsfile <- sumstatsfile |>
      dplyr::select(-dplyr::any_of("Blast_hits")) |>
      dplyr::left_join(counts, by = c("Markername" = "qaccver")) |>
      dplyr::mutate(Blast_hits = dplyr::coalesce(.data$Blast_hits, 0L))
    
    # restore Blast_hits to its original position if it existed before
    if (!is.na(.pos_bh)) {
      .cols_now  <- names(sumstatsfile)
      .cols_wo   <- setdiff(.cols_now, "Blast_hits")
      .new_order <- append(.cols_wo, "Blast_hits", after = .pos_bh - 1)
      sumstatsfile <- sumstatsfile[, .new_order, drop = FALSE]
    }
    
    rm(counts)
    
    # Coverage/identity filter counts 
    cov_file <- file.path(output.path,
                      paste0(probe.name, "_with_", genome.name, "_all_mappings.csv"))

    cov_counts <- data.table::fread(
      cmd          = sprintf("grep -v '^#' %s", shQuote(cov_file)),
      select       = "qaccver",
      showProgress = TRUE
    ) |>
      dplyr::filter(!is.na(qaccver), qaccver != "") |>
      dplyr::count(qaccver, name = "Cov-ident_BLASThits")

    .orig_names <- names(sumstatsfile)
    .pos_filter <- match("Filter_coverage-identity", .orig_names)
    .pos_covhit <- match("Cov-ident_BLASThits", .orig_names)
    
    .present <- sumstatsfile$Markername %in% cov_counts$qaccver
    
    # apply updates
    sumstatsfile <- sumstatsfile |>
      dplyr::select(-dplyr::any_of(c("Filter_coverage-identity", "Cov-ident_BLASThits"))) |>
      dplyr::mutate(`Filter_coverage-identity` = ifelse(.present, "PASS", "FAIL")) |>
      dplyr::left_join(cov_counts, by = c("Markername" = "qaccver")) |>
    dplyr::mutate(`Cov-ident_BLASThits` = dplyr::coalesce(.data$`Cov-ident_BLASThits`, 0L))    
    
.cols_now <- names(sumstatsfile)
    
    # move Filter_coverage_identity back
    if (!is.na(.pos_filter)) {
      .wo <- setdiff(.cols_now, "Filter_coverage-identity")
      sumstatsfile <- sumstatsfile[, append(.wo, "Filter_coverage-identity", after = .pos_filter - 1), drop = FALSE]
      .cols_now <- names(sumstatsfile)
    }
    
    # move Cov-ident_BLASThits back
    if (!is.na(.pos_covhit)) {
      .wo <- setdiff(.cols_now, "Cov-ident_BLASThits")
      sumstatsfile <- sumstatsfile[, append(.wo, "Cov-ident_BLASThits", after = .pos_covhit - 1), drop = FALSE]
    }
    
    # no rm(cov_counts) here, we use it again later
    
    #Secondaries /hybridisation filter counts
sec_file <- file.path(output.path,
                      paste0(probe.name, "_with_", genome.name, "_filtered_mappings.csv"))

    sec_counts <- data.table::fread(
      cmd          = sprintf("grep -v '^#' %s", shQuote(sec_file)),
      select       = "qaccver",
      showProgress = TRUE
    ) |>
      dplyr::filter(!is.na(qaccver), qaccver != "") |>
      dplyr::count(qaccver, name = "Secondaries_hybridisation_BLASThits")
    
    .orig_names  <- names(sumstatsfile)
    .pos_secflag <- match("Filter_secondaries_and_hybridisation", .orig_names)
    .pos_sechits <- match("Secondaries_hybridisation_BLASThits", .orig_names)
    
    .present_sec <- sumstatsfile$Markername %in% sec_counts$qaccver
    
    .already_failed <- if ("Filter_coverage-identity" %in% names(sumstatsfile)) {
      sumstatsfile$`Filter_coverage-identity` != "PASS"
    } else {
      rep(FALSE, nrow(sumstatsfile))
    }
    
    sumstatsfile <- sumstatsfile |>
      dplyr::select(-dplyr::any_of(c("Filter_secondaries_and_hybridisation", "Secondaries_hybridisation_BLASThits"))) |>
      dplyr::mutate(
        Filter_secondaries_and_hybridisation = dplyr::case_when(
          .already_failed ~ "ALREADY REMOVED",
          .present_sec    ~ "PASS",
          TRUE            ~ "FAIL"
        )
      ) |>
      dplyr::left_join(sec_counts, by = c("Markername" = "qaccver")) |>
      dplyr::mutate(Secondaries_hybridisation_BLASThits =
                  dplyr::coalesce(.data$Secondaries_hybridisation_BLASThits, 0L))
    
    .cols_now <- names(sumstatsfile)
    
    if (!is.na(.pos_secflag)) {
      .wo <- setdiff(.cols_now, "Filter_secondaries_and_hybridisation")
      sumstatsfile <- sumstatsfile[, append(.wo, "Filter_secondaries_and_hybridisation", after = .pos_secflag - 1), drop = FALSE]
      .cols_now <- names(sumstatsfile)
    }
    if (!is.na(.pos_sechits)) {
      .wo <- setdiff(.cols_now, "Secondaries_hybridisation_BLASThits")
      sumstatsfile <- sumstatsfile[, append(.wo, "Secondaries_hybridisation_BLASThits", after = .pos_sechits - 1), drop = FALSE]
    }
    
    rm(sec_counts)
    

int_file <- file.path(output.path,
                      paste0(probe.name, "_with_", genome.name, "_intermediate_filtering_hits.csv"))

    int_counts <- data.table::fread(
      cmd          = sprintf("grep -v '^#' %s", shQuote(int_file)),
      select       = "qaccver",
      showProgress = TRUE
    ) |>
      dplyr::filter(!is.na(qaccver), qaccver != "") |>
      dplyr::count(qaccver, name = "Intermediate_BLASThits")


    #Intermediate sensitivity filter
    .orig_names <- names(sumstatsfile)
    .pos_iflag  <- match("Filter_intermediate_sensitivity", .orig_names)
    .pos_ihits  <- match("Intermediate_BLASThits", .orig_names)
    
    .present_int <- sumstatsfile$Markername %in% int_counts$qaccver
    
    .already_removed <- {
      failed_cov <- if ("Filter_coverage-identity" %in% names(sumstatsfile)) {
        sumstatsfile$`Filter_coverage-identity` != "PASS"
      } else rep(FALSE, nrow(sumstatsfile))
      failed_sec <- if ("Filter_secondaries_and_hybridisation" %in% names(sumstatsfile)) {
        sumstatsfile$Filter_secondaries_and_hybridisation != "PASS"
      } else rep(FALSE, nrow(sumstatsfile))
      failed_cov | failed_sec
    }
    
    sumstatsfile <- sumstatsfile |>
      dplyr::select(-dplyr::any_of(c("Filter_intermediate_sensitivity", "Intermediate_BLASThits"))) |>
      dplyr::mutate(
        `Filter_intermediate_sensitivity` = dplyr::case_when(
          .already_removed ~ "ALREADY REMOVED",
          .present_int     ~ "PASS",
          TRUE             ~ "FAIL"
        )
      ) |>
      dplyr::left_join(int_counts, by = c("Markername" = "qaccver")) |>
      dplyr::mutate(`Intermediate_BLASThits` =
                  dplyr::coalesce(.data$`Intermediate_BLASThits`, 0L))   
 
    .cols_now <- names(sumstatsfile)
    
    if (!is.na(.pos_iflag)) {
      .wo <- setdiff(.cols_now, "Filter_intermediate_sensitivity")
      sumstatsfile <- sumstatsfile[, append(.wo, "Filter_intermediate_sensitivity", after = .pos_iflag - 1), drop = FALSE]
      .cols_now <- names(sumstatsfile)
    }
    if (!is.na(.pos_ihits)) {
      .wo <- setdiff(.cols_now, "Intermediate_BLASThits")
      sumstatsfile <- sumstatsfile[, append(.wo, "Intermediate_BLASThits", after = .pos_ihits - 1), drop = FALSE]
    }
    
    rm(int_counts)
    
    filteredstrict <- read.csv(
      file         = paste0(output.path,"/", probe.name, "_with_", genome.name,"_strict_filtering_hits.csv"),
      comment.char = "#"
    )
    
    #Strict sensitivity filter
    .orig_names <- names(sumstatsfile)
    .pos_sflag  <- match("Filter_strict_sensitivity", .orig_names)
    .pos_shits  <- match("Strict_BLASThits", .orig_names)
    
    .present_str <- sumstatsfile$Markername %in% filteredstrict$qaccver
    str_counts <- filteredstrict |>
      dplyr::count(qaccver, name = "Strict_BLASThits")
    
    .already_removed <- {
      failed_cov <- if ("Filter_coverage-identity" %in% names(sumstatsfile)) {
        sumstatsfile$`Filter_coverage-identity` != "PASS"
      } else rep(FALSE, nrow(sumstatsfile))
      failed_sec <- if ("Filter_secondaries_and_hybridisation" %in% names(sumstatsfile)) {
        sumstatsfile$Filter_secondaries_and_hybridisation != "PASS"
      } else rep(FALSE, nrow(sumstatsfile))
      failed_int <- if ("Filter_intermediate_sensitivity" %in% names(sumstatsfile)) {
        sumstatsfile$`Filter_intermediate_sensitivity` != "PASS"
      } else rep(FALSE, nrow(sumstatsfile))
      failed_cov | failed_sec | failed_int
    }
    
    sumstatsfile <- sumstatsfile |>
      dplyr::select(-dplyr::any_of(c("Filter_strict_sensitivity", "Strict_BLASThits"))) |>
      dplyr::mutate(
        `Filter_strict_sensitivity` = dplyr::case_when(
          .already_removed ~ "ALREADY REMOVED",
          .present_str     ~ "PASS",
          TRUE             ~ "FAIL"
        )
      ) |>
      dplyr::left_join(str_counts, by = c("Markername" = "qaccver")) |>
      dplyr::mutate(`Strict_BLASThits` =
                      dplyr::coalesce(.data$`Strict_BLASThits`, 0L))
    
    .cols_now <- names(sumstatsfile)
    
    if (!is.na(.pos_sflag)) {
      .wo <- setdiff(.cols_now, "Filter_strict_sensitivity")
      sumstatsfile <- sumstatsfile[, append(.wo, "Filter_strict_sensitivity", after = .pos_sflag - 1), drop = FALSE]
      .cols_now <- names(sumstatsfile)
    }
    if (!is.na(.pos_shits)) {
      .wo <- setdiff(.cols_now, "Strict_BLASThits")
      sumstatsfile <- sumstatsfile[, append(.wo, "Strict_BLASThits", after = .pos_shits - 1), drop = FALSE]
    }
    
    ## Add BLAST fail flag and re-do coverage/identity summary using cov_counts
    sumstatsfile <- sumstatsfile |>
      dplyr::mutate(Filter_BLAST = ifelse(.data$Blast_hits > 0L, "PASS", "FAIL")) |>
      dplyr::relocate(Filter_BLAST, .after = "Markername")
    
    sumstatsfile <- sumstatsfile |>
      dplyr::select(-dplyr::any_of(c("Filter_coverage-identity", "Cov-ident_BLASThits"))) |>
      dplyr::left_join(cov_counts, by = c("Markername" = "qaccver")) |>
      dplyr::mutate(
        `Cov-ident_BLASThits` = dplyr::coalesce(.data$`Cov-ident_BLASThits`, 0L),
        `Filter_coverage-identity` = dplyr::case_when(
          .data$Filter_BLAST == "FAIL" ~ "ALREADY REMOVED",
          .present                     ~ "PASS",
          TRUE                         ~ "FAIL"
        )
      ) |>
      dplyr::relocate(`Filter_coverage-identity`, .after = "Filter_BLAST") |>
      dplyr::relocate(`Cov-ident_BLASThits`, .after = "Filter_coverage-identity") |>
      dplyr::relocate(`Blast_hits`, .after = "Filter_BLAST")
    
    # Plotting & histogram 
    stages <- c(
      "Failed_BLAST",
      "Filter_coverage-identity",
      "Filter_secondaries_and_hybridisation",
      "Filter_intermediate_sensitivity",
      "Filter_strict_sensitivity"
    )
    failed_blast <- sumstatsfile$Filter_BLAST == "FAIL"
    
    first_fail_stage <- sapply(seq_len(nrow(sumstatsfile)), function(i) {
      if (failed_blast[i]) return("Failed_BLAST")
      row <- as.character(sumstatsfile[i, c("Filter_coverage-identity",
                                            "Filter_secondaries_and_hybridisation",
                                            "Filter_intermediate_sensitivity",
                                            "Filter_strict_sensitivity")])
      idx <- which(row != "PASS")
      if (length(idx) == 0) return(NA_character_)
      c("Filter_coverage-identity",
        "Filter_secondaries_and_hybridisation",
        "Filter_intermediate_sensitivity",
        "Filter_strict_sensitivity")[min(idx)]
    })
    
    total_markers <- nrow(sumstatsfile)
    
    markers_dropped_df <- dplyr::tibble(Stage = stages) |>
      dplyr::left_join(
        dplyr::tibble(Stage = first_fail_stage) |>
          dplyr::count(Stage, name = "Dropped"),
        by = "Stage"
      ) |>
      dplyr::mutate(
        Dropped    = dplyr::coalesce(.data$Dropped, 0L),
        Pct_of_all = 100 * .data$Dropped / total_markers
      )
    
    markers_remaining     <- sum(is.na(first_fail_stage))
    markers_remaining_pct <- 100 * markers_remaining / total_markers
    
    markers_plot_df <- dplyr::bind_rows(
      markers_dropped_df |> dplyr::transmute(Category = Stage, Count = Dropped),
      dplyr::tibble(Category = "Remaining_after_strict", Count = markers_remaining)
    ) |>
      dplyr::mutate(Pct_of_all = 100 * .data$Count / total_markers)
    
    markers_plot_df$Category <- factor(
      markers_plot_df$Category,
      levels = c(stages, "Remaining_after_strict")
    )
    
    tot_unfiltered <- sum(sumstatsfile$Blast_hits)
    tot_cov        <- sum(sumstatsfile$`Cov-ident_BLASThits`)
    tot_sec        <- sum(sumstatsfile$`Secondaries_hybridisation_BLASThits`)
    tot_int        <- sum(sumstatsfile$Intermediate_BLASThits)
    tot_strict     <- sum(sumstatsfile$Strict_BLASThits)
    
    drop_cov <- max(0, tot_unfiltered - tot_cov)
    drop_sec <- max(0, tot_cov        - tot_sec)
    drop_int <- max(0, tot_sec        - tot_int)
    drop_str <- max(0, tot_int        - tot_strict)
    
    hits_plot_df <- dplyr::tibble(
      Category = factor(
        c("Dropped_at_coverage-identity",
          "Dropped_at_secondaries/hybridisation",
          "Dropped_at_intermediate",
          "Dropped_at_strict",
          "Remaining_after_strict"),
        levels = c("Dropped_at_coverage-identity",
                   "Dropped_at_secondaries/hybridisation",
                   "Dropped_at_intermediate",
                   "Dropped_at_strict",
                   "Remaining_after_strict")
      ),
      Count = c(drop_cov, drop_sec, drop_int, drop_str, tot_strict)
    ) |>
      dplyr::mutate(
        Pct_of_initial = if (tot_unfiltered > 0) 100 * .data$Count / tot_unfiltered else 0
      )
    
    palette5 <- c("#0072B2", "#E69F00", "#999999", "#D55E00", "#009E73")
    palette6 <- c("#CC79A7", "#0072B2", "#E69F00", "#999999", "#D55E00", "#009E73")
    
    p1 <- ggplot2::ggplot(markers_plot_df, ggplot2::aes(x = Category, y = Pct_of_all, fill = Category)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = palette6, name = "Stage") +
      ggplot2::labs(
        x = NULL,
        y = "Percentage of markers",
        title = "Marker outcomes across filtering stages",
        subtitle = sprintf("Remaining after strict: %.1f%%", markers_remaining_pct)
      ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")
    
    print(p1)
    
    ggplot2::ggsave(file = paste0(output.path,"/Figs/",probe.name, "_with_", genome.name,"_marker_filtering_stages.png"),
                    dpi = 400)
    
    p2 <- ggplot2::ggplot(hits_plot_df, ggplot2::aes(x = Category, y = Pct_of_initial, fill = Category)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = palette5, name = "Stage") +
      ggplot2::labs(
        x = NULL,
        y = "Percentage of initial BLAST hits",
        title = "BLAST hits dropped at each stage",
        subtitle = if (tot_unfiltered > 0)
          sprintf("Initial hits: %s | Remaining after strict: %.1f%%",
                  format(tot_unfiltered, big.mark = ","),
                  100 * tot_strict / tot_unfiltered) else
                    "No initial hits"
      ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")
    
    print(p2)
    
    ggplot2::ggsave(file = paste0(output.path,"/Figs/",probe.name, "_with_", genome.name,"_BLASTnhits_filtering_stages.png"),
                    dpi = 400)
    
    .safe_int <- function(x) as.integer(dplyr::coalesce(x, 0L))
    
    df_unf <- dplyr::tibble(
      Stage = "Unfiltered",
      Hits  = .safe_int(sumstatsfile$Blast_hits)
    )
    
    df_cov <- sumstatsfile |>
      dplyr::filter(`Filter_coverage-identity` == "PASS") |>
      dplyr::transmute(Stage = "Coverage-Identity",
                       Hits  = .safe_int(.data$`Cov-ident_BLASThits`))
    
    df_sec <- sumstatsfile |>
      dplyr::filter(`Filter_secondaries_and_hybridisation` == "PASS") |>
      dplyr::transmute(Stage = "Secondaries/Hybridisation",
                       Hits  = .safe_int(.data$`Secondaries_hybridisation_BLASThits`))
    
    df_int <- sumstatsfile |>
      dplyr::filter(`Filter_intermediate_sensitivity` == "PASS") |>
      dplyr::transmute(Stage = "Intermediate sensitivity",
                       Hits  = .safe_int(.data$`Intermediate_BLASThits`))
    
    df_str <- sumstatsfile |>
      dplyr::filter(`Filter_strict_sensitivity` == "PASS") |>
      dplyr::transmute(Stage = "Strict sensitivity",
                       Hits  = .safe_int(.data$`Strict_BLASThits`))
    
    hist_long_pass <- dplyr::bind_rows(df_unf, df_cov, df_sec, df_int, df_str)
    
    hist_long_pass$Stage <- factor(
      hist_long_pass$Stage,
      levels = c("Unfiltered",
                 "Coverage-Identity",
                 "Secondaries/Hybridisation",
                 "Intermediate sensitivity",
                 "Strict sensitivity")
    )
    
    bins <- c("1","2","3","4","5","6","7","8","9","10+")
    hist_long_pass_binned <- hist_long_pass |>
      dplyr::filter(Hits >= 1L) |>
      dplyr::mutate(
        Hits_bin = dplyr::if_else(Hits >= 10L, "10+", as.character(Hits)),
        Hits_bin = factor(Hits_bin, levels = bins, ordered = TRUE)
      )
    
    y_breaks <- c(1, 5, 20, 100, 500, 2000, 10000, 50000, 200000)
    
    p_hist_pass <- ggplot2::ggplot(hist_long_pass_binned, ggplot2::aes(x = Hits_bin, fill = Stage)) +
      ggplot2::geom_bar() +
      ggplot2::facet_wrap(~ Stage, ncol = 1, scales = "fixed") +
      ggplot2::scale_fill_manual(values = palette5, guide = "none") +
      ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::scale_y_continuous(trans = "log10", breaks = y_breaks) +
      ggplot2::labs(
        x = "BLAST hits per marker (1-9, 10+)",
        y = "Markers (frequency, log10)",
        title = "BLASTn hit count distributions for markers"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        strip.text  = ggplot2::element_text(face = "bold"),
        plot.title  = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank()
      )
    
    print(p_hist_pass)
    
    ggplot2::ggsave(file = paste0(output.path,"/Figs/",probe.name, "_with_", genome.name,"_histogram_hits.png"),
                    dpi = 400)
    

   # Add in information for duplication and NULL markers

 duphitsinfo <- read.table(duplicationmap.hits, sep='\t', header=TRUE)

# weird spelling because I really don't want to write null into a variable name in any way haha.
 duphitsinfo$nulcall <-NA 
 
 missing_markers <- setdiff(unique(sumstatsfile$Markername), unique(duphitsinfo$qaccver))

# add NA-filled rows for all missing markers. 
if (length(missing_markers) > 0) {
  na_rows <- duphitsinfo[0, , drop = FALSE]              # keep same columns/types
  na_rows <- na_rows[rep(1, length(missing_markers)), ]  # make N empty rows
  na_rows$qaccver <- missing_markers

  duphitsinfo <- dplyr::bind_rows(duphitsinfo, na_rows)
}

# Begin populating the required values for NULLS
duphitsinfo <- dplyr::select(
  dplyr::mutate(
    dplyr::left_join(
      duphitsinfo,
      dplyr::select(sumstatsfile, Markername, Blast_hits, `Cov-ident_BLASThits`),
      by = c("qaccver" = "Markername")
    ),
    nulcall = ifelse(
      (!is.na(Blast_hits) & Blast_hits == 0) |
        (!is.na(`Cov-ident_BLASThits`) & `Cov-ident_BLASThits` == 0),
      "Yes",
      NA_character_
    )
  ),
  -Blast_hits,
  -`Cov-ident_BLASThits`
)
    write.table(
      duphitsinfo,
      file=(paste0(probe.name, "_with_", genome.name, "_marker_localdups_NULLS_counts.tsv")),
      quote = FALSE,
      sep="\t",
      row.names = FALSE
    )


    write.csv(
      sumstatsfile,
      file.path(paste0(probe.name, "_with_", genome.name, "_summary_filtering.csv")),
      quote = FALSE,
      row.names = FALSE
    )
    
    # prepping VCF file creation
    target_ids <- markersall
    
    target_info <- dplyr::tibble(
      qaccver = target_ids,
      qseq    = as.character(target_df$Sequence),
      tbase   = as.character(target_df$Target.base)
    )
    
    ab <- stringr::str_match(target_info$tbase, "\\[([ACGT])\\/([ACGT])\\]")
    target_info$ClusterA_NT <- ab[, 2]
    target_info$ClusterB_NT <- ab[, 3]
    target_info$tbase <- NULL
    
    strict_unique <- filteredstrict |>
      dplyr::distinct(qaccver, .keep_all = TRUE)
    
    present_ids <- intersect(target_info$qaccver, strict_unique$qaccver)
    failed_ids  <- setdiff(target_info$qaccver,  strict_unique$qaccver)
    
    strict_present <- strict_unique |>
      dplyr::semi_join(dplyr::tibble(qaccver = present_ids), by = "qaccver") |>
      dplyr::arrange(match(qaccver, target_info$qaccver))
    
    template <- filteredstrict[0, , drop = FALSE]
    n_fail   <- length(failed_ids)
    
    na_like <- function(x, n) {
      if (is.integer(x)) rep(NA_integer_, n)
      else if (is.numeric(x)) rep(NA_real_, n)
      else if (inherits(x, "Date")) rep(as.Date(NA), n)
      else if (inherits(x, "POSIXt")) rep(as.POSIXct(NA), n)
      else rep(NA_character_, n)
    }
    
    failed_full <- as.data.frame(lapply(template, na_like, n = n_fail), stringsAsFactors = FALSE)
    
    saveRDS(failed_full,"failedfull.RDS")
    saveRDS(target_info,"targetinfo.RDS")
    saveRDS(template,"template.RDS")
    saveRDS(strict_present,"strict_present.RDS")
    saveRDS(failed_ids,"failed_ids.RDS")
    
    failed_full$qaccver <- failed_ids
    if ("saccver" %in% names(failed_full)) failed_full$saccver <- "chrUnk"
    if ("SNPpos" %in% names(failed_full)) {
      failed_full$SNPpos <- if (is.integer(template$SNPpos)) as.integer(0) else 0
    }
    
    failed_info <- target_info |>
      dplyr::semi_join(dplyr::tibble(qaccver = failed_ids), by = "qaccver")
    
    if ("qseq" %in% names(failed_full))        failed_full$qseq        <- failed_info$qseq
    if ("ClusterA_NT" %in% names(failed_full)) failed_full$ClusterA_NT <- failed_info$ClusterA_NT
    if ("ClusterB_NT" %in% names(failed_full)) failed_full$ClusterB_NT <- failed_info$ClusterB_NT
    
    strict_present$SNPpos <- as.numeric(strict_present$SNPpos)
    
    all_markers_1to1 <- dplyr::bind_rows(
      strict_present[, names(template), drop = FALSE],
      failed_full[,   names(template), drop = FALSE]
    )
    
all_markers_1to1_subset <- dplyr::select(all_markers_1to1,qaccver, saccver, SNPpos,Ref,ALT,sstrand)



out_1to1 <- file.path(
  output.path,
  paste0(probe.name, "_with_", genome.name, "Brioche_all_markers1to1stagingforvcf.csv")
)

# Build metadata header lines
meta_lines <- c(
  "## Brioche mapping metadata",
  paste0("## brioche_version (commit version)=",   md$brioche_version),
  paste0("## brioche_repo=",md$brioche_repo_url),
  paste0("## brioche_branch=",md$brioche_commit_url),
  paste0("## pident_threshold=",  md$pident),
  paste0("## coverage_threshold=",md$coverage),
  paste0("## additional_blast_options=", md$otherblastoptions),
  paste0("## used_target_chromosome=",   md$usetargetchrom),
  paste0("## used_shared_markers_map=",  md$usesharedmarkersmap),
  paste0("## used_ld_edge_map=",         md$useldedgemap),
  paste0("## used_genetic_map=",         md$usegeneticmap),
  paste0("## run_date=", md$run_date)
)

# Write comments + CSV
con <- file(out_1to1, open = "wt")
on.exit(try(close(con), silent = TRUE), add = TRUE)
writeLines(meta_lines, con = con)
utils::write.table(
  all_markers_1to1_subset,
  file      = con,
  sep       = ",",
  quote     = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

  }
