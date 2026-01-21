#' @title Read in file with comment lines and attach those lines to read in object
#' @description This function is used to read in csv/tsv files in the Brioche
#'   pipelines to preserve their appended headers of filtering information.
#'   It avoids loading the whole file into a single string, so it is safe for
#'   very large files.
#' @author James O'Dwyer
#' @param path Character. Path to the file to read.
#' @param sep Field separator for the data section (default `","`). Use `"\t"` for TSV.
#' @param header Logical; whether the first non-header line contains column names. Default `TRUE`.
#' @param parse_kv Logical; if `TRUE`, parse `## key=value` lines into `attr(df, "header_meta")`.
#' @param ... Additional arguments passed to `utils::read.csv()`.
#' @keywords genome, csv, comment lines
#'
#' @return A data.frame like from `read.csv`, with header lines stored in
#'   attributes: `attr(df, "raw_header_lines")` and (optionally) `attr(df, "header_meta")`.
#'
#' @export
read_mixed_csv_base <- function(path,
                                sep = ",",
                                header = TRUE,
                                parse_kv = TRUE,
                                ...) {

  # Read the leading "##"
  con <- file(path, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  header_lines <- character()
  n_header <- 0L

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break  # EOF
    if (startsWith(line, "##")) {
      header_lines <- c(header_lines, line)
      n_header <- n_header + 1L
    } else {
      # First non-header line reached; stop scanning.
      break
    }
  }

  # If there were header lines, skip them; otherwise read from top.
  if (n_header > 0L) {
    df <- utils::read.csv(
      file = path,
      sep = sep,
      header = header,
      stringsAsFactors = FALSE,
      skip = n_header,
      ...
    )
  } else {
    df <- utils::read.csv(
      file = path,
      sep = sep,
      header = header,
      stringsAsFactors = FALSE,
      ...
    )
  }

  # parse ##  into header_meta
  if (parse_kv && length(header_lines)) {
    kv_lines <- sub("^##\\s*", "", header_lines)
    kv_lines <- kv_lines[grepl("=", kv_lines, fixed = TRUE)]

    if (length(kv_lines)) {
      # strip inline comments that start with '#'
      kv_nocomment <- sub("\\s+#.*$", "", kv_lines)
      keys <- trimws(sub("^([^=]+)=.*$", "\\1", kv_nocomment))
      vals <- trimws(sub("^[^=]+=(.*)$", "\\1", kv_nocomment))

      # coerce to simple types when obvious (numeric/logical), else keep as character
      coerce1 <- function(x) {
        if (grepl("^[-+]?[0-9]*\\.?[0-9]+$", x)) return(as.numeric(x))
        if (tolower(x) %in% c("true","false","yes","no")) {
          return(tolower(x) %in% c("true","yes"))
        }
        x
      }
      vals <- lapply(vals, coerce1)

      meta <- as.list(vals)
      names(meta) <- make.names(keys, unique = TRUE)
      attr(df, "header_meta") <- meta
    }
  }

  attr(df, "raw_header_lines") <- header_lines
  df
}

# ---- Usage ----
# df <- read_mixed_csv_base("your_file.csv")
# attributes(df)$header_meta         # list of parsed parameters (if any)
# attributes(df)$raw_header_lines    # original '##' lines

# If you want to write the CSV back out and re-prepend the saved header:
# Generic writer with header stitching for CSV/TSV (no extra packages)

#' @title write out file with comment lines
#' @description This function is used to write out csv files in the Brioche pipelines to preserve their appended headers of filtering information
#' @author James O'Dwyer
#' @param df Data frame to write. If it has `attr(df, "raw_header_lines")`, those lines
#'   are written before the header row unless `header_override` is supplied.
#' @param path Character. Output file path.
#' @param sep Field separator. If `NULL` (default), inferred from `path`
#' @param row.names Logical; passed to `utils::write.table()`. Default `FALSE`.
#' @param na String used for missing values. Default `""`.
#' @param flatten_lists Logical; collapse list-columns to `;`-separated strings. Default `TRUE`.
#' @param header_override Optional character vector of header lines to write instead of
#'   `attr(df, "raw_header_lines")`.
#' @param add_blank_after_header Logical; write a blank line after the header lines. Default `FALSE`.
#' @param quote Logical; quote character columns. Passed to `utils::write.table()`. Default `TRUE`.
#' @param ... Additional arguments passed to `utils::write.table()`.
#'
#' @return Invisibly returns `normalizePath(path)`. Side effect: writes `path`.
#'
#' @keywords genome, csv, comment lines, write
#' @export
#'
write_with_header_base <- function(df, path,
                                   sep = NULL,              # NULL => infer from file extension
                                   row.names = FALSE,
                                   na = "",
                                   flatten_lists = TRUE,
                                   header_override = NULL,  # character vector of header lines
                                   add_blank_after_header = FALSE,
                                   quote = TRUE,            # quote character columns
                                   ...) {
  
  if (is.null(sep)) {
    ext <- tolower(sub(".*\\.(\\w+)$", "\\1", path))
    sep <- if (ext %in% c("tsv", "tab", "txt")) "\t" else ","
  }
  

  if (!inherits(df, "data.frame")) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
  } else {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
  }
  

  if (flatten_lists) {
    for (nm in names(df)) {
      if (is.list(df[[nm]]) && !is.data.frame(df[[nm]])) {
        df[[nm]] <- vapply(
          df[[nm]],
          function(x) paste(x, collapse = ";"),
          FUN.VALUE = character(1),
          USE.NAMES = FALSE
        )
      }
    }
  }
  

  header_lines <- header_override
  if (is.null(header_lines)) header_lines <- attr(df, "raw_header_lines")
  if (is.null(header_lines)) header_lines <- character(0)
  

  tmp <- tempfile(fileext = ".tbl")
  on.exit(unlink(tmp), add = TRUE)
  
  utils::write.table(
    df,
    file = tmp,
    sep = sep,
    row.names = row.names,
    col.names = TRUE,
    qmethod = "double",
    na = na,
    quote = quote,
    ...
  )
  
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  
  if (length(header_lines)) {
    writeLines(header_lines, con = con, sep = "\n")
    if (add_blank_after_header) writeChar("\n", con, eos = NULL)
  }
  writeLines(readLines(tmp, warn = FALSE), con = con, sep = "\n")
  
  invisible(normalizePath(path))
}



# Drop groups, flatten list/data.frame columns, ensure each column has length nrow(df)

#' @title write out file with comment lines
#' @description This function is used to sanitise brioche blast hit objects to ensure each column or row has the same length
#' @author James O'Dwyer
#' @param x Data-frame-like object to sanitise.
#' @return A `data.frame` with list/data.frame columns collapsed to one string per row.
#' @keywords genome, sanitize, rectangular
#' @export
sanitize_rectangular <- function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  n <- NROW(x)
  
  for (nm in names(x)) {
    col <- x[[nm]]
    
    # Flatten list-cols (including NULL elements)
    if (is.list(col) && !is.data.frame(col)) {
      x[[nm]] <- vapply(seq_len(n), function(i) {
        el <- col[[i]]
        if (is.null(el)) "" else paste(as.character(el), collapse = ";")
      }, FUN.VALUE = character(1), USE.NAMES = FALSE)
    }
    
    # Collapse embedded data.frame columns to one string per row
    if (is.data.frame(x[[nm]])) {
      subdf <- x[[nm]]
      x[[nm]] <- vapply(seq_len(n), function(i) {
        paste(unlist(subdf[i, , drop = FALSE]), collapse = ";")
      }, FUN.VALUE = character(1), USE.NAMES = FALSE)
    }    
    if (length(x[[nm]]) == 1L) x[[nm]] <- rep(x[[nm]], n)
    if (length(x[[nm]]) != n)  x[[nm]] <- as.character(x[[nm]])[seq_len(n)]
  }
  
  rownames(x) <- NULL
  x
}




#' @title DoIntermediatefiltering 
#' @description This function is used to undertake an intermediate level of filtering producing output files with only the best blast hits per marker chromosome
#' @author James O'Dwyer
#' @param probe.name name of the probe file being analysed
#' @param genome.name name of the reference genome being compared
#' @param output.path path were the results are saved
#' @param dotargetfilt whether to do target chromosome based filtering
#' @param chrom.comp Location of chromosome comparison pairs between old an new references
#' @param target.data 2 column table with first column showin allele ID and second expected target chrom to preference
#' @param domapmarkers whether to do key linked marker based filtering
#' @param similarity.maps X column table where the first column is the target marker and n columns after contain a marker known to be linked to the target marker
#' @param doldedge whether to do linkage disequilibrium based filtering
#' @param ldmapp 10 column table (standard .ld format) for pairwise ld scores and locations per marker
#' @param blast.hits Low filtered blast hits file as csv
#' @param mappings.file Low filtered mappings file as tsv
#' @param dogeneticmap whether to use genetic mapping prior
#' @param geneticmap.file CSV path to genetic mapping prior (MarkerName, Chromosomes_mapped, etc)
#' @param dup.dist distance for detection of local duplications
#' @param keeplocalduppos yes or no for whether to keep markers which contain local duplications if there are no other hits and all hits show a consensus variant site
#' @keywords genome, filter blast intermediate
#' @return Returns an object which has been filtered using given priors to remove/keep best hits per marker
#' @export
DoIntermediatefiltering <-
  function(probe.name,
           genome.name,
           output.path,
           dotargetfilt,
           chrom.comp,
           target.data,
           domapmarkers,
           similarity.maps,
           doldedge,
           ldmapp,
           blast.hits,
           mappings.file,
           dogeneticmap,
           geneticmap.file,
           dup.dist,
           keeplocalduppos)
  {
    
    hdr <- readLines(blast.hits, n = 100L, warn = FALSE)
    headerinfo <- hdr[grepl("^##", hdr)]


    blast.out <- read.csv(
      blast.hits,
      comment.char = "#",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )


    mappings.out <- read.table(mappings.file,sep="\t",header = TRUE,comment.char="#",fill=TRUE)
    mappings.out <- subset(mappings.out, mappings.out$Hybridised!="No" | is.na(mappings.out$Hybridised))
    
    blast.out$bitscore <- as.numeric(blast.out$bitscore)
    
    # The labels are wrong for Alternate_SNP_ID first thing is to fix them uo 
    
    blast.out <- blast.out |>
      dplyr::arrange(qaccver, saccver, Coverage, SNPpos) |>
      dplyr::group_by(qaccver) |>
      dplyr::mutate(
        .alt_n = dplyr::row_number(),
        Alternate_SNP_ID = paste0(qaccver, ".", .alt_n)
      ) |>
      dplyr::ungroup() |>
      dplyr::select(-.alt_n)
    
    mappings.out <- mappings.out |>
      dplyr::arrange(SNP_ID, Chr, Coverage, SNP_position) |>
      dplyr::group_by(SNP_ID) |>
      dplyr::mutate(
        .alt_n = dplyr::row_number(),
        Alternate_SNP_ID = paste0(SNP_ID, ".", .alt_n)
      ) |>
      dplyr::ungroup() |>
      dplyr::select(-.alt_n)

    saveRDS(blast.out, "blastout.RDS")
    saveRDS(mappings.out, "mappingsout.RDS")

    # Need to remove any < or > present SNP rows as they should not make it past intermediate filtering. 
    # just need to remove from blast.out as mappings.out is filtered through blast remaining regardless 

    blast.out <- blast.out[ !grepl("[<>]", as.character(blast.out$SNPpos)), , drop = FALSE ]        


    if (!dir.exists(output.path))
      dir.create(output.path, recursive = T)
    
    # Generate two separate filtering processes depending on the presence/absence of chrom info 
    
    if(dotargetfilt=="No" & domapmarkers=="No" & doldedge=="No" & dogeneticmap=="No") {

     # Filter to have only the top bitscore hit per chrom
      top_hits <- blast.out |>
        dplyr::group_by(qaccver, saccver) |>
        dplyr::filter(bitscore == max(bitscore, na.rm = TRUE)) |>
        dplyr::ungroup()
      
      # Filter to have only SNPs present in top_hits present in the mappings file
      filtered_mappings <- mappings.out |>
        dplyr::filter(Alternate_SNP_ID %in% top_hits$Alternate_SNP_ID)
    
      
      attributes(top_hits)$raw_header_lines <- headerinfo
      top_hits_clean <- sanitize_rectangular(top_hits)
  
  
      attributes(filtered_mappings)$raw_header_lines <- headerinfo
      filtered_mappings_clean <- sanitize_rectangular(filtered_mappings)
  
    }
    # There are 3 possible data streams to incorporate. First round of ifs will load in requisite components. Second will do analysis
    
    if(dotargetfilt=="Yes") { 
      
      # Read in Chromosome comparison file
      chromcomparisontable <- read.table(chrom.comp,header=TRUE,sep=",")
      # Remove excess empty rows
      chromcomparisontable[chromcomparisontable == ""] <- NA               # treat "" as missing (character cols)
      chromcomparisontable <- chromcomparisontable[rowSums(is.na(chromcomparisontable)) != ncol(chromcomparisontable), ]
      
      # Read in list of markers with a known chromosome to preference. 
      knowntargetchroms <- read.table(target.data,header=TRUE,sep=",")
      # Remove excess empty rows
      knowntargetchroms[knowntargetchroms == ""] <- NA               # treat "" as missing (character cols)
      knowntargetchroms <- knowntargetchroms[rowSums(is.na(knowntargetchroms)) != ncol(knowntargetchroms), ]
      #convert old chrom to new chrom data for knowntargetchroms
      lkp <- setNames(chromcomparisontable$New_reference_chromosome,
                      chromcomparisontable$Original_reference_chromsome)
      
      # map; unmatched stay NA
      knowntargetchroms$newChrom <- unname(lkp[knowntargetchroms$Chrom])

    saveRDS(chromcomparisontable, "chromecomparisontable.RDS")
    saveRDS(knowntargetchroms, "knowntargetchroms.RDS")

      
    }
    if(domapmarkers=="Yes") { 
      
      neighbourmarkers <- read.table(similarity.maps,header=TRUE,sep=",")
      # Remove excess empty rows
      neighbourmarkers[neighbourmarkers == ""] <- NA               # treat "" as missing (character cols)
      neighbourmarkers <- neighbourmarkers[rowSums(is.na(neighbourmarkers)) != ncol(neighbourmarkers), ]
      
    saveRDS(neighbourmarkers, "neighbourmarkers.RDS")
    }
    
    if(doldedge=="Yes") { 
      
      pwldmappings <- read.table(ldmapp, header = TRUE, sep="")
      # Remove excess empty rows
      pwldmappings[pwldmappings == ""] <- NA               # treat "" as missing (character cols)
      pwldmappings <- pwldmappings[rowSums(is.na(pwldmappings)) != ncol(pwldmappings), ]
      
      if ("R2" %in% names(pwldmappings)) {
        pwldmappings[["R2"]] <- suppressWarnings(as.numeric(pwldmappings[["R2"]]))
        pwldmappings <- pwldmappings[ is.na(pwldmappings[["R2"]]) | pwldmappings[["R2"]] >= 0.2, , drop = FALSE ]
      } else if ("R" %in% names(pwldmappings)) {
        pwldmappings[["R"]] <- suppressWarnings(as.numeric(pwldmappings[["R"]]))
        pwldmappings <- pwldmappings[ is.na(pwldmappings[["R"]]) | pwldmappings[["R"]] >= 0.45, , drop = FALSE ]
      }

      keep <- c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B",
                if ("R2" %in% names(pwldmappings)) "R2" else if ("R" %in% names(pwldmappings)) "R" else NULL)
      
      pwldmappings <- pwldmappings[ , intersect(keep, names(pwldmappings)), drop = FALSE ]
      
      ## choose score column & thresholds
      ## pick score column and thresholds
      score_col <- if ("R2" %in% names(pwldmappings)) "R2" else "R"
      thr <- if (score_col == "R2") c(0.8, 0.6, 0.4) else c(0.9, 0.78, 0.63)
      
      ## coerce score to numeric
      pwldmappings[[score_col]] <- suppressWarnings(as.numeric(pwldmappings[[score_col]]))
      
      ## split indices by SNP_A
      idx_by_snp <- split(seq_len(nrow(pwldmappings)), pwldmappings$SNP_A)
      
      keep <- logical(nrow(pwldmappings))
      for (idx in idx_by_snp) {
        v <- pwldmappings[[score_col]][idx]
        
        # choose the highest threshold with >5 qualifying rows
        chosen <- NA_real_
        for (t in thr) {
          if (sum(v >= t, na.rm = TRUE) > 5L) { chosen <- t; break }
        }
        
        if (is.na(chosen)) {
          # no threshold reached: keep ALL rows for this SNP_A (unchanged behavior)
          keep[idx] <- TRUE
        } else {
          # keep exactly the TOP 5 rows among those meeting the chosen threshold
          qual <- idx[which(v >= chosen)]
          # order those by descending score (break ties by original order)
          ord <- qual[order(-pwldmappings[[score_col]][qual], seq_along(qual))]
          keep[ord[seq_len(min(5L, length(ord)))]] <- TRUE
          # (rows under threshold are dropped; NAs are excluded from top-5 selection)
        }
      }
      
      pwldmappings <- pwldmappings[keep, , drop = FALSE]
          saveRDS(pwldmappings, "LDS.RDS")
      
      
       
    }    
        
    




    # first level, chrom mapping information is present but specific marker neighbours or LD aren't

    if(dotargetfilt=="Yes") {
      
      # Lookup: qaccver -> target chromosome
      lkp <- knowntargetchroms |>
        dplyr::transmute(qaccver = MarkerID, target_chr = newChrom)
      
      
      targetchrommarkers <- knowntargetchroms$MarkerID
      
      
      # Attach target chr to blast, and compute per-(qaccver, saccver) maxima
      tmp <- blast.out |>
        dplyr::left_join(lkp, by = "qaccver") |>
        dplyr::group_by(qaccver, saccver) |>
        dplyr::mutate(max_in_group = max(bitscore, na.rm = TRUE)) |>
        dplyr::ungroup()
      
      # Compute target_max per qaccver from rows ON the target chromosome only
      target_max_df <- tmp |>
        dplyr::filter(!is.na(target_chr), saccver == target_chr) |>
        dplyr::group_by(qaccver) |>
        dplyr::summarise(target_max = max(bitscore, na.rm = TRUE), .groups = "drop") |>
        # if all NA on target, max(...) is -Inf -> treat as NA (no usable target max)
        dplyr::mutate(target_max = dplyr::na_if(target_max, -Inf))
      
      # Join target_max back and flag usable target
      tmp <- tmp |>
        dplyr::left_join(target_max_df, by = "qaccver") |>
        dplyr::mutate(has_target = !is.na(target_chr) & !is.na(target_max))
      
      # Keep:
      # - all rows on the target chr (when usable),
      # - plus rows on other chrs only if they tie the target chr max,
      # - otherwise (no usable target) keep per-(qaccver, saccver) max rows
      # - saving output as linked to dataset type for combined analysis below in scenarios where multiple priors are given
      top_hits <- top_hits_knowntargetchroms <- tmp |>
        dplyr::filter(
          dplyr::case_when(
            has_target ~ (saccver == target_chr) | (bitscore == target_max),
            TRUE       ~ bitscore == max_in_group
          )
        ) |>
        dplyr::select(-max_in_group, -has_target)
      
      # Filter to have only SNPs present in top_hits present in the mappings file
      
    }
    
    if(domapmarkers=="Yes") {
     
      
      targetneighbourmarkers <- neighbourmarkers$TargetMarker

      
      
      nm_clean <- neighbourmarkers |>
        dplyr::mutate(
          dplyr::across(
            dplyr::everything(),
            ~ { x <- base::as.character(.); x <- base::trimws(x); dplyr::na_if(x, "") }
          )
        ) |>
        dplyr::filter(!base::is.na(TargetMarker))
      
      ## 1) Long pivot only the linkedmarker* columns (drops NA linked markers) (Have not decided on max cols for number of linked markers yet but 4-6 seems sufficient)
      neigh_long <- nm_clean |>
        tidyr::pivot_longer(
          cols = tidyselect::matches("^linkedmarker\\d+$"),
          names_to = "linked_col",
          values_to = "linked_marker",
          values_drop_na = TRUE
        ) |>
        dplyr::rename(target = TargetMarker)
      
      ## 2) Join linked markers to their BLASTn hits
      link_hits <- neigh_long |>
        dplyr::inner_join(blast.out, by = c("linked_marker" = "qaccver"))
      
      ## 3) Recurrent chromosome per target (tie-break: most hits, then higher max bitscore, then saccver)
      recurrent <- link_hits |>
        dplyr::group_by(target, saccver) |>
        dplyr::summarise(
          n_hits   = dplyr::n(),
          max_bits = base::max(bitscore, na.rm = TRUE),
          .groups  = "drop_last"
        ) |>
        dplyr::arrange(
          dplyr::desc(n_hits),
          dplyr::desc(dplyr::coalesce(max_bits, -Inf)),
          saccver,
          .by_group = TRUE
        ) |>
        dplyr::slice(1) |>
        dplyr::ungroup() |>
        dplyr::transmute(qaccver = target, recurrent_chr = saccver)
      
      ## 4) Attach recurrent chr; compute per-(qaccver, saccver) maxima to use as base for all other hits to other chroms
      tmp <- blast.out |>
        dplyr::left_join(recurrent, by = "qaccver") |>
        dplyr::group_by(qaccver, saccver) |>
        dplyr::mutate(max_in_group = base::max(bitscore, na.rm = TRUE)) |>
        dplyr::ungroup()
      
      ## 5) Target max on the recurrent chr 
      target_max_df <- tmp |>
        dplyr::filter(!base::is.na(recurrent_chr), saccver == recurrent_chr) |>
        dplyr::group_by(qaccver) |>
        dplyr::summarise(target_max = base::max(bitscore, na.rm = TRUE), .groups = "drop") |>
        dplyr::mutate(target_max = dplyr::na_if(target_max, -Inf))
      
      tmp2 <- tmp |>
        dplyr::left_join(target_max_df, by = "qaccver") |>
        dplyr::mutate(has_recurrent = !base::is.na(recurrent_chr) & !base::is.na(target_max))
      
      ## 6) Filter with prior; use per-chromosome max if marker not present in prior marker list
      top_hits <- top_hits_neighbourfilt <- tmp2 |>
        dplyr::filter(
          dplyr::case_when(
            has_recurrent ~ (saccver == recurrent_chr) | (bitscore == target_max),
            TRUE          ~ bitscore == max_in_group
          )
        ) |>
        dplyr::select(-max_in_group, -has_recurrent)
      
    }

    
    if(doldedge=="Yes") {
      
      ## 1) Build prior: for each SNP_A, pick the most frequent chromosome of its SNP_B partners (note less meningful relationships have been prefiltered already)
      ld_hits <- pwldmappings |>
        dplyr::select(SNP_A, SNP_B) |>
        dplyr::filter(!base::is.na(SNP_A), !base::is.na(SNP_B)) |>
        dplyr::inner_join(blast.out, by = c("SNP_B" = "qaccver"))   # brings in saccver, bitscore for SNP_Bs
      
      # Extract the number of unique markers that are present in the Brioche output from the other 
      targetLDmarkers <- unique(intersect(blast.out$qaccver,pwldmappings$SNP_A))
      
      if(nrow(ld_hits) >=1) {
      
      prior <- ld_hits |>
        dplyr::group_by(SNP_A, saccver) |>
        dplyr::summarise(
          n_hits   = dplyr::n(),
          max_bits = base::max(bitscore, na.rm = TRUE),
          .groups  = "drop_last"
        ) |>
        # most frequent; break ties by higher max bitscore, then by saccver
        dplyr::arrange(
          dplyr::desc(n_hits),
          dplyr::desc(dplyr::coalesce(max_bits, -Inf)),
          saccver,
          .by_group = TRUE
        ) |>
        dplyr::slice(1) |>
        dplyr::ungroup() |>
        dplyr::transmute(qaccver = SNP_A, prior_chr = saccver)
      
      ## 2) Attach prior to all BLAST rows; compute per-(qaccver, saccver) maxima for fallback
      tmp <- blast.out |>
        dplyr::left_join(prior, by = "qaccver") |>
        dplyr::group_by(qaccver, saccver) |>
        dplyr::mutate(max_in_group = base::max(bitscore, na.rm = TRUE)) |>
        dplyr::ungroup()
      
      ## 3) For targets with a prior chromosome, get that targets max bitscore ON that chr (no warnings)
      target_max_df <- tmp |>
        dplyr::filter(!base::is.na(prior_chr), saccver == prior_chr) |>
        dplyr::group_by(qaccver) |>
        dplyr::summarise(
          target_max = {
            bs <- bitscore
            if (base::all(base::is.na(bs))) NA_real_
            else base::max(bs, na.rm = TRUE)
          },
          .groups = "drop"
        )
      
      tmp2 <- tmp |>
        dplyr::left_join(target_max_df, by = "qaccver") |>
        dplyr::mutate(has_prior = !base::is.na(prior_chr) & !base::is.na(target_max))
      
      ## 4) Filter:
      ##   - with prior: keep ALL rows on prior chr + any rows on other chrs that TIE that chroms max
      ##   - without prior: keep per chrom-(qaccver, saccver) max rows 
      top_hits <- top_hits_ld_filt <- tmp2 |>
        dplyr::filter(
          dplyr::case_when(
            has_prior ~ (saccver == prior_chr) | (bitscore == target_max),
            TRUE      ~ bitscore == max_in_group
          )
        ) |>
        dplyr::select(-max_in_group, -has_prior)
      
      }
      else{
        
        top_hits <- top_hits_ld_filt <- blast.out |>
          dplyr::group_by(qaccver, saccver) |>
          dplyr::filter(bitscore == max(bitscore, na.rm = TRUE)) |>
          dplyr::ungroup()
        
      }
    }
    # combine results of targe chroms + neighbour markers
    if(dotargetfilt=="Yes" & domapmarkers=="Yes" & doldedge=="No") {
    
    # Value of similar markers is higher than target filt chrom so the base will be targetfilt
      
      # confirm there are no duplicates in the marker list
      replace_ids <- unique(as.character(targetneighbourmarkers))
      
      ## 1) Take the rows for the markers which were directly informed by the nearby marker prior information
      neigh_rows <- top_hits_neighbourfilt |>
        dplyr::semi_join(
          dplyr::tibble(qaccver = replace_ids),
          by = "qaccver"
        )
      

      known_remaining <- top_hits_knowntargetchroms |>
        dplyr::anti_join(
          dplyr::tibble(qaccver = replace_ids),
          by = "qaccver"
        )
      
      filtered_mappings <- dplyr::bind_rows(known_remaining, neigh_rows)
      
      
      
    }


if (dogeneticmap == "Yes") {

  # 1) Read priors (CSV with commas, header)
  gmap <- read.table(geneticmap.file, sep = ",", header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)

  # 2) For each row, pick the ChrmapNo* corresponding to the largest No_maps*
  chr_cols <- paste0("ChrmapNo", 1:5)
  num_cols <- paste0("No_maps",  1:5)


  no_mat  <- as.matrix(sapply(gmap[num_cols], as.numeric))
  best_ix <- max.col(no_mat, ties.method = "first")

  chr_mat <- as.matrix(gmap[chr_cols])
  chosen_chr_oldref <- chr_mat[cbind(seq_len(nrow(chr_mat)), best_ix)]

  # 3) Map old->new chromosomes via chrom.comp (CSV with commas, header)
  chromcomparisontable <- read.table(chrom.comp, sep = ",", header = TRUE,
                                     stringsAsFactors = FALSE, check.names = FALSE)
  lkp <- setNames(chromcomparisontable$New_reference_chromosome,
                  chromcomparisontable$Original_reference_chromsome)
  chosen_chr_newref <- unname(lkp[chosen_chr_oldref])

  # 4) Build prior table in new-ref space
  prior_gmap <- data.frame(
    qaccver = gmap$MarkerName,
    geneticmap_chr = chosen_chr_newref,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  targetgeneticmapmarkers <- prior_gmap$qaccver[!is.na(prior_gmap$geneticmap_chr)]

  # 5) Apply prior: keep rows on the prior chr + any ties to its max; else per-chr max fallback
  tmp <- blast.out |>
    dplyr::left_join(prior_gmap, by = "qaccver") |>
    dplyr::group_by(qaccver, saccver) |>
    dplyr::mutate(max_in_group = max(bitscore, na.rm = TRUE)) |>
    dplyr::ungroup()

  target_max_df <- tmp |>
    dplyr::filter(!is.na(geneticmap_chr), saccver == geneticmap_chr) |>
    dplyr::group_by(qaccver) |>
    dplyr::summarise(target_max = max(bitscore, na.rm = TRUE), .groups = "drop")

  tmp2 <- tmp |>
    dplyr::left_join(target_max_df, by = "qaccver") |>
    dplyr::mutate(has_gmap = !is.na(geneticmap_chr) & !is.na(target_max))

  top_hits <- top_hits_gmap_filt <- tmp2 |>
    dplyr::filter(
      dplyr::case_when(
        has_gmap ~ (saccver == geneticmap_chr) | (bitscore == target_max),
        TRUE     ~ bitscore == max_in_group
      )
    ) |>
    dplyr::select(-max_in_group, -has_gmap)

}


## New revised method which directly layers based on what priors are given cutting down on the if statements and code duplication. 
# test this method and confirm the results are identical (will test with genetic mapping for wheat but also need to test again for barley to make sure stuff comes out identical)
overlay_by_ids <- function(base_df, add_df, replace_ids) {
  if (is.null(base_df)) return(add_df)
  if (is.null(add_df) || length(replace_ids) == 0L) return(base_df)
  keep_base <- dplyr::anti_join(base_df, dplyr::tibble(qaccver = replace_ids), by = "qaccver")
  add_rows  <- add_df |>
    dplyr::semi_join(dplyr::tibble(qaccver = replace_ids), by = "qaccver")
  dplyr::bind_rows(keep_base, add_rows)
}

has_ld    <- identical(doldedge,      "Yes")
has_gmap  <- identical(dogeneticmap,  "Yes")
has_tchr  <- identical(dotargetfilt,  "Yes")
has_neigh <- identical(domapmarkers,  "Yes")

# Decide baseline (lowest available priority); if none, we keep whatever 'top_hits'
# was produced in the earlier "no priors" branch.
if (has_ld || has_gmap || has_tchr || has_neigh) {
  base_hits <- NULL
  if (has_ld)        base_hits <- top_hits_ld_filt
  else if (has_gmap) base_hits <- top_hits_gmap_filt
  else if (has_tchr) base_hits <- top_hits_knowntargetchroms
  else               base_hits <- top_hits_neighbourfilt

  # Overlay in ascending priority so later steps override earlier ones
  if (has_gmap)  base_hits <- overlay_by_ids(
                    base_hits, top_hits_gmap_filt,
                    unique(as.character(targetgeneticmapmarkers))
                  )
  if (has_tchr)  base_hits <- overlay_by_ids(
                    base_hits, top_hits_knowntargetchroms,
                    unique(as.character(targetchrommarkers))
                  )
  if (has_neigh) base_hits <- overlay_by_ids(
                    base_hits, top_hits_neighbourfilt,
                    unique(as.character(targetneighbourmarkers))
                  )

  top_hits <- base_hits
}
    
    filtered_mappings <- mappings.out |>
      dplyr::filter(Alternate_SNP_ID %in% top_hits$Alternate_SNP_ID)
    
    
    attributes(top_hits)$raw_header_lines <- headerinfo
    top_hits_clean <- sanitize_rectangular(top_hits)
    
    top_hits_clean <- top_hits_clean[,1:34]

    
    attributes(filtered_mappings)$raw_header_lines <- headerinfo
    filtered_mappings_clean <- sanitize_rectangular(filtered_mappings)
             

    hits.outpath <- file.path(output.path,
                                   paste0(probe.name, "_with_", genome.name, "_intermediate_filtering_hits.csv"))
    
    
    
    write_with_header_base(df = top_hits_clean,
              path = paste0(probe.name, "_with_", genome.name, "_intermediate_filtering_hits.csv"),sep=",",quote=FALSE)
    
    
    mappings.outpath <- file.path(output.path,
                             paste0(probe.name, "_with_", genome.name, "_intermediate_filtering_mappings.tsv"))
    
    write_with_header_base(df = filtered_mappings_clean,
                           path = paste0(probe.name, "_with_", genome.name, "_intermediate_filtering_mappings.tsv"),sep="\t",quote=FALSE)
    
    
    cat("\nNumber of rows in intermediate filtering results blast table:", nrow(top_hits_clean), "\n")
    cat("\nMapping results written to ", output.path)


   ###################
   # Adding new section to track the putative presence of local duplication 
   ###################

  
   duplications <- as.data.frame(matrix(nrow=(length(unique(top_hits$qaccver))),ncol=4))
   colnames(duplications) <- c("qaccver","copy_number","consensusbase", "keep")


top_hits2 <- dplyr::mutate(
  top_hits,
  .row_id = dplyr::row_number()
)

group_info <- dplyr::summarise(
  dplyr::group_by(top_hits2, qaccver),
  n_rows    = dplyr::n(),
  n_saccver = dplyr::n_distinct(saccver),
  .groups   = "drop"
)

# Choose highest bitscore; if tied, keep the earliest row in the original data
anchor <- dplyr::transmute(
  dplyr::ungroup(
    dplyr::slice(
      dplyr::arrange(
        dplyr::filter(
          dplyr::group_by(top_hits2, qaccver),
          dplyr::n() > 1,
          dplyr::n_distinct(saccver) == 1
        ),
        dplyr::desc(bitscore),
        .row_id,
        .by_group = TRUE
      ),
      1
    )
  ),
  qaccver,
  saccver,
  anchor_sstart = sstart
)

dup.dist <- as.numeric(dup.dist)

in_range <- dplyr::filter(
  dplyr::inner_join(top_hits2, anchor, by = c("qaccver", "saccver")),
  !is.na(sstart),
  !is.na(anchor_sstart),
  sstart >= pmax(anchor_sstart - dup.dist, 0),
  sstart <= (anchor_sstart + dup.dist)
)

counts_in_range <- dplyr::count(in_range, qaccver, name = "copy_number")

# Consensus base check in range (TRUE only if all Ref values identical AND none are NA)
consensus_in_range <- dplyr::summarise(
  dplyr::group_by(in_range, qaccver),
  consensusbase = dplyr::if_else(
    all(!is.na(Ref)) && dplyr::n_distinct(Ref) == 1,
    "true",
    "false"
  ),
  .groups = "drop"
)

# Build duplications table with required NA behaviour when multiple saccver exist
duplications <- dplyr::select(
  dplyr::mutate(
    dplyr::left_join(
      dplyr::left_join(group_info, counts_in_range, by = "qaccver"),
      consensus_in_range,
      by = "qaccver"
    ),
    copy_number = dplyr::case_when(
      n_rows == 1 ~ 1L,
      n_rows > 1 & n_saccver == 1 ~ as.integer(copy_number),
      n_rows > 1 & n_saccver > 1 ~ NA_integer_,
      TRUE ~ NA_integer_
    ),
    consensusbase = dplyr::case_when(
      n_rows == 1 ~ "true",
      n_rows > 1 & n_saccver == 1 ~ consensusbase,
      n_rows > 1 & n_saccver > 1 ~ NA_character_,
      TRUE ~ NA_character_
    )
  ),
  qaccver,
  copy_number,
  consensusbase
)

if (keeplocalduppos == "yes") {
  duplications$keep <- ifelse(
    !is.na(duplications$copy_number) &
      duplications$copy_number > 0 &
      duplications$consensusbase == "true",
    "yes",
    "no"
  )
} else {
  duplications$keep <- "NA"
}
       
   write.table(x = duplications, file = paste0(probe.name, "_with_", genome.name, "_marker_localduplications_counts.tsv"),sep="\t",quote=FALSE,row.names=FALSE)


  }
