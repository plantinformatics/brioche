#' @title Read in file with comment lines and attach those lines to read in object
#' @description This function is used to read in csv files in the Brioche pipelines to preserve their appended headers of filtering information
#' @author James O'Dwyer
#' @param path Character. Path to the file to read.
#' @param sep Field separator for the data section (default `","`). Use `"\t"` for TSV.
#' @param header Logical; whether the first non-header line contains column names. Default `TRUE`.
#' @param parse_kv Logical; if `TRUE`, parse `## key=value` lines into `attr(df, "header_meta")`.
#' @param ... Additional arguments passed to `utils::read.csv()`.
#' @keywords genome, csv, comment lines
#'
#' @return Returns an object like a read.table or read.csv function with comment.char applied but with the comment lines attached to the attributes of the object so they are not lost
#'
#' @export
read_mixed_csv_base <- function(path, sep = ",", header = TRUE, parse_kv = TRUE, ...) {
  lines <- readLines(path, warn = FALSE)
  if (!length(lines)) stop("Empty file: ", path)
  
  # 1) Split header vs CSV: header lines start with '##'
  is_header <- grepl("^##", lines)
  header_lines <- lines[is_header]
  data_lines   <- lines[!is_header]
  
  if (!length(data_lines)) stop("No CSV body detected after header lines.")
  
  # 2) Read CSV body
  csv_text <- paste(data_lines, collapse = "\n")
  df <- read.csv(text = csv_text, sep = sep, header = header, stringsAsFactors = FALSE, ...)
  
  # 3) Optional: parse '## key = value # comment' into a named list
  if (parse_kv && length(header_lines)) {
    # keep only lines that look like key=value after the leading ##
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
        if (tolower(x) %in% c("true","false","yes","no")) return(tolower(x) %in% c("true","yes"))
        x
      }
      vals <- lapply(vals, coerce1)
      
      meta <- as.list(vals)
      names(meta) <- make.names(keys, unique = TRUE)
      attr(df, "header_meta") <- meta
    }
  }
  
  # Always keep the raw header lines too
  attr(df, "raw_header_lines") <- header_lines
  df
}

# usage
# df <- read_mixed_csv_base("your_file.csv")
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
  
  # ---- 0) Decide separator ----
  if (is.null(sep)) {
    ext <- tolower(sub(".*\\.(\\w+)$", "\\1", path))
    sep <- if (ext %in% c("tsv", "tab", "txt")) "\t" else ","
  }
  
  # ---- 1) Coerce to plain data.frame ----
  if (!inherits(df, "data.frame")) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
  } else {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
  }
  
  # ---- 2) Flatten list-columns so write.table can always serialize ----
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
  
  # ---- 3) Gather header lines (any length allowed) ----
  header_lines <- header_override
  if (is.null(header_lines)) header_lines <- attr(df, "raw_header_lines")
  if (is.null(header_lines)) header_lines <- character(0)
  
  # ---- 4) Write the table body to a temp file ----
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
  
  # ---- 5) Stitch header + body to final path ----
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
    
    # After coercions, enforce column length == n
    if (length(x[[nm]]) == 1L) x[[nm]] <- rep(x[[nm]], n)
    if (length(x[[nm]]) != n)  x[[nm]] <- as.character(x[[nm]])[seq_len(n)]
  }
  
  rownames(x) <- NULL
  x
}








#' @title Do strict filtering 
#' @description This function is used to undertake an strict level of filtering producing output files with only one hit per marker allowed
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
#' @param blast.hits Intermediate filtered blast hits file as csv
#' @param mappings.file Low filtered mappings file as tsv
#' @param dogeneticmap whether to use genetic mapping prior
#' @param geneticmap.file CSV path to genetic mapping prior (MarkerName, Chromosomes_mapped, etc)
#' @param doorientation Whether to update marker orientation file ("Yes"/"No", case insensitive).
#' @param orientation.file TSV file with columns ID, Orientation giving current marker orientation.
#' @param project.dir Project directory where Updated_orientation_file.tsv will be written.
#' @keywords genome, filter blast strict
#' @export
#'
#' @return Returns an object which has been filtered using given priors to remove/keep a maximum of 1 hit (the best) per marker
#'
DoStrictfiltering <-
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
           doorientation     = "No",
           orientation.file  = NA,
           project.dir       = output.path)
  {    
    # Remember that the header Lines will be stored as an attribute for reattaching later
    blast.out <- read_mixed_csv_base(blast.hits)
    headerinfo <- attributes(blast.out)
    mappings.out <- read.table(mappings.file,sep="\t",header = TRUE)
 
    blast.out$bitscore <- as.numeric(blast.out$bitscore)
    
    saveRDS(blast.out, "blastout.RDS")
    saveRDS(mappings.out, "mappingsout.RDS")

    
    if (!dir.exists(output.path))
      dir.create(output.path, recursive = T)
    
    remove_namescombined <- vector()
    keep_namescombined <- vector()
    # Generate two separate filtering processes depending on the presence/absence of priors info 
    # Currently set to just remove any marker which has 2+ remaining hits. May revised and bring back but change the Chrom
    # to be Chrom or something
    if(dotargetfilt=="No" & domapmarkers=="No" & doldedge=="No" & dogeneticmap=="No") {
      
      
      top_hits <- blast.out |>
        dplyr::add_count(qaccver, name = "n") |>
        dplyr::filter(n == 1) |>
        dplyr::select(-n)
      
      
      remove_names <- setdiff(blast.out$Alternate_SNP_ID, top_hits$Alternate_SNP_ID)
      remove_namescombined <- append(remove_namescombined,remove_names)

      keep_names <- intersect(top_hits$Alternate_SNP_ID, blast.out$Alternate_SNP_ID)
      keep_namescombined <- append(keep_namescombined, keep_names)

      
      
    }
    
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
      
    }
    if(domapmarkers=="Yes") { 
      
      neighbourmarkers <- read.table(similarity.maps,header=TRUE,sep=",")
      # Remove excess empty rows
      neighbourmarkers[neighbourmarkers == ""] <- NA               # treat "" as missing (character cols)
      neighbourmarkers <- neighbourmarkers[rowSums(is.na(neighbourmarkers)) != ncol(neighbourmarkers), ]
      
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
      
    }
    
    if(dotargetfilt=="Yes") {
      # System here will be to 
      #1. Keep all markers where only 1 hit was present
      #2. If multiple hits are present for a marker but only one is on the target chromosome, keep that marker hit only
      #3. If multiple hits are present for a marker on the target chromosome, drop that marker
      #4. If multiple hits are present and no target chrom info is avaliable, drop that marker. 

      lkp <- knowntargetchroms |>
        dplyr::transmute(qaccver = MarkerID, target_chr = newChrom)
      
      targetchrommarkers <- knowntargetchroms$MarkerID
      
      tmp <- blast.out |>
        dplyr::left_join(lkp, by = "qaccver")
      
      # 2) Per-marker stats
      per_q <- tmp |>
        dplyr::group_by(qaccver) |>
        dplyr::summarise(
          total_hits  = dplyr::n(),
          target_chr  = dplyr::first(target_chr),
          n_on_target = base::sum(saccver == dplyr::first(target_chr), na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::mutate(has_target = !base::is.na(target_chr))
      
      # 3) Keep rules
      tmp2 <- tmp |>
        dplyr::left_join(per_q, by = "qaccver", suffix = c("", ".perq")) |>
        dplyr::mutate(
          # unify target_chr (prefer the one already on tmp)
          target_chr = dplyr::coalesce(target_chr, target_chr.perq),
          keep = dplyr::case_when(
            !base::is.na(target_chr) & total_hits == 1L ~ TRUE,                                   # single hit keep
            !base::is.na(target_chr) & total_hits > 1L & n_on_target == 1L ~ (saccver == target_chr), # one on target keep that row
            !base::is.na(target_chr) & n_on_target > 1L ~ FALSE,                                   # multiple on target drop marker
            base::is.na(target_chr) & total_hits > 1L ~ FALSE,                                     # no target & >1 hits drop marker
            base::is.na(target_chr) & total_hits == 1L ~ TRUE,                                     # no target & single hit keep
            TRUE ~ FALSE
          )
        ) |>
        # clean up any per_q duplicates
        dplyr::select(-dplyr::matches("\\.perq$"))
      
      # 4) Output
      top_hits <- top_hits_knowntargetchroms_strict <- tmp2 |>
        dplyr::filter(keep) |>
        dplyr::select(-total_hits, -n_on_target, -has_target, -keep)
      
      
      remove_nameschromchrom <- setdiff(blast.out$Alternate_SNP_ID, top_hits$Alternate_SNP_ID)
      remove_namescombined <- append(remove_namescombined,remove_nameschromchrom)

      keep_nameschromchrom <- intersect(top_hits$Alternate_SNP_ID, blast.out$Alternate_SNP_ID)
      keep_namescombined <- append(keep_namescombined, keep_nameschromchrom)

      
    }
    
 

if (dogeneticmap == "Yes") {

  # Read inputs
  gmap <- read.table(geneticmap.file, sep = ",", header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)
  chromcomparisontable <- read.table(chrom.comp, sep = ",", header = TRUE,
                                     stringsAsFactors = FALSE, check.names = FALSE)

  chr_cols <- paste0("ChrmapNo", 1:5)

 saveRDS(gmap , "gmap.RDS")
 saveRDS(blast.out , "blastout.RDS")
 saveRDS(chromcomparisontable, "chromchrom.RDS")

  # Count mapped chromosomes; prefer provided column
  nmapped <- if ("Chromosomes_mapped" %in% names(gmap)) {
    suppressWarnings(as.numeric(gmap$Chromosomes_mapped))
  } else {
    rowSums(!is.na(gmap[chr_cols]))
  }
  nmapped[is.na(nmapped)] <- 0

  # Remove priors with >=2 mapped chromosomes (ambiguous as priors)
  gmap_slim <- gmap[nmapped < 2, , drop = FALSE]

  # Partition IDs
  zero_ids <- gmap$MarkerName[nmapped == 0]  # no chromosome in gmap
  one_ids  <- gmap$MarkerName[nmapped == 1]  # exactly one chromosome in gmap

  # For one_ids: pick the single non-NA ChrmapNo* (old ref)
  chosen_chr_oldref <- rep(NA_character_, length(one_ids))
  if (length(one_ids) > 0) {
    chr_sub <- gmap[gmap$MarkerName %in% one_ids, chr_cols, drop = FALSE]
    chosen_chr_oldref <- apply(chr_sub, 1, function(r) {
      r <- as.character(r); r <- trimws(r)
      r[which(!is.na(r) & r != "")[1]]
    })
  }

  # Map old -> new reference (needed only for multi-hit markers in one_ids)
  lkp <- setNames(chromcomparisontable$New_reference_chromosome,
                  chromcomparisontable$Original_reference_chromsome)
  mapped_chr_new <- unname(lkp[chosen_chr_oldref])

  prior_gmap <- data.frame(
    qaccver = one_ids,
    geneticmap_chr = mapped_chr_new,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  prior_gmap <- prior_gmap[!is.na(prior_gmap$geneticmap_chr), , drop = FALSE]

  # Total BLAST hit counts per marker
  tgt_counts <- blast.out |>
    dplyr::count(qaccver, name = "total_hits")

  # 1) Single-hit markers: keep regardless of gmap (covers zero_ids, one_ids, or no gmap row)
  winners_single <- blast.out |>
    dplyr::left_join(tgt_counts, by = "qaccver") |>
    dplyr::filter(total_hits == 1L) |>
    dplyr::group_by(qaccver) |>
    dplyr::slice(1) |>
    dplyr::ungroup()

  # 2) Multi-hit markers with exactly one gmap chromosome:
  #    keep ONLY if exactly one BLAST hit lies on the mapped chromosome
  multi_with_prior <- blast.out |>
    dplyr::left_join(tgt_counts, by = "qaccver") |>
    dplyr::filter(total_hits > 1L) |>
    dplyr::inner_join(prior_gmap, by = "qaccver") |>
    dplyr::filter(saccver == geneticmap_chr)

  winners_multi <- multi_with_prior |>
    dplyr::group_by(qaccver) |>
    dplyr::filter(dplyr::n() == 1L) |>
    dplyr::ungroup()

  # Note: multi-hit markers in zero_ids (Chromosomes_mapped==0) have no prior_gmap,
  #       so they are NOT in winners_multi and are therefore dropped.

  # 3) Combine winners ? one row per kept marker
  top_hits_gmap_strict <- dplyr::bind_rows(winners_single, winners_multi) |>
    dplyr::distinct(qaccver, .keep_all = TRUE)

  # 4) Update keep/remove name collections
  keep_names_gmap   <- intersect(top_hits_gmap_strict$Alternate_SNP_ID, blast.out$Alternate_SNP_ID)
  remove_names_gmap <- setdiff(blast.out$Alternate_SNP_ID, keep_names_gmap)

  keep_namescombined   <- append(keep_namescombined,   keep_names_gmap)
  remove_namescombined <- append(remove_namescombined, remove_names_gmap)
}

    
    if(domapmarkers=="Yes") {
      
      # So two possibilities for how to do the strict filtering here
      # 1. I could do a chrom chrom match like I did in intermediate filtering and then apply stricter thresholds like above
      # 2. I could extract the chrom + bp mapping position of each of the linked markers and choose which blast hit is the closest. I could also 
      # pre-filter the linked here so that it only contains hits with 1 match. 
      
      targetneighbourmarkers <- neighbourmarkers$TargetMarker
      
      
      blast.out[["SNPpos"]] <- suppressWarnings(as.numeric(blast.out[["SNPpos"]]))
      
      neigh_long <- neighbourmarkers |>
        dplyr::mutate(dplyr::across(dplyr::everything(), ~{
          x <- base::as.character(.); x <- base::trimws(x); dplyr::na_if(x, "")
        })) |>
        dplyr::filter(!base::is.na(TargetMarker)) |>
        tidyr::pivot_longer(
          cols = tidyselect::matches("^linkedmarker\\d+$"),
          names_to = "linked_col", values_to = "linked_marker",
          values_drop_na = TRUE
        ) |>
        dplyr::rename(target = TargetMarker)
      
      # Keep only linked markers that have a SINGLE BLAST hit so that ambig in the linked markers is not an issue
      link_counts <- blast.out |>
        dplyr::count(qaccver, name = "n_hits_link")
      
      good_links <- neigh_long |>
        dplyr::inner_join(link_counts, by = c("linked_marker" = "qaccver")) |>
        dplyr::filter(n_hits_link == 1L) |>
        dplyr::select(target, linked_marker)
      
      ## Pull their single BLAST hit to get chr & bp. We now have a chromosome list and list of exact bps of nearest marker
      link_hits <- good_links |>
        dplyr::inner_join(blast.out, by = c("linked_marker" = "qaccver")) |>
        dplyr::transmute(
          target,
          neigh_chr = saccver,
          neigh_bp  = SNPpos
        )
      
      # Neighbour prior per target:
      # most frequent chromosome among remaining links
      # keep a reference bp as the median(neigh_bp) for tie-breaking & distance (very very unlikely to ever be needed)
      prior <- link_hits |>
        dplyr::group_by(target, neigh_chr) |>
        dplyr::summarise(
          n_links = dplyr::n(),
          med_bp  = stats::median(neigh_bp, na.rm = TRUE),
          .groups = "drop_last"
        ) |>
        dplyr::arrange(dplyr::desc(n_links), neigh_chr, .by_group = TRUE) |>
        dplyr::slice(1) |>
        dplyr::ungroup() |>
        dplyr::transmute(qaccver = target, prior_chr = neigh_chr, prior_bp = med_bp)
      
      ## Per-target total blast hit counts
      tgt_counts <- blast.out |>
        dplyr::count(qaccver, name = "total_hits")
      
      ## Attach prior & counts to all target BLAST rows
      tmp <- blast.out |>
        dplyr::left_join(prior,      by = "qaccver") |>
        dplyr::left_join(tgt_counts, by = "qaccver") |>
        dplyr::mutate(on_prior_chr = !base::is.na(prior_chr) & (saccver == prior_chr))
      
      # winners WITH prior:
      #  1. if zero hits on prior chr  no row will be selected (marker dropped)
      #  2.if one hit on prior chr   that one wins
      #  3. if >1 on prior chr  choose the one with SNPpos closest to prior_bp
      # winners with prior: choose the one on prior chr closest to prior_bp
      winners_with_prior <- tmp |>
        dplyr::filter(!base::is.na(prior_chr), saccver == prior_chr) |>
        dplyr::filter(!base::is.na(SNPpos), !base::is.na(prior_bp)) |>
        dplyr::group_by(qaccver) |>
        dplyr::arrange(base::abs(SNPpos - prior_bp),
                       dplyr::desc(dplyr::coalesce(bitscore, -Inf))) |>
        dplyr::slice(1) |>
        dplyr::ungroup() |>
        dplyr::select(qaccver, saccver, SNPpos, dplyr::everything())
      
      ## winners with no prior: keep only markers with exactly one total hit
      winners_no_prior_single <- tmp |>
        dplyr::filter(base::is.na(prior_chr), total_hits == 1L) |>
        dplyr::group_by(qaccver) |>
        dplyr::slice(1) |>
        dplyr::ungroup() |>
        dplyr::select(qaccver, saccver, SNPpos, dplyr::everything())
      
      ## Combine the two winner sets (exactly one row per kept qaccver)
      top_hits <- filteredtop_hits_neighbour_strict <- dplyr::bind_rows(
        winners_with_prior,
        winners_no_prior_single
      ) |>
        dplyr::distinct(qaccver, .keep_all = TRUE) |>
        dplyr::select(-total_hits, -prior_chr, -prior_bp, -on_prior_chr)
     
      remove_namessimmarkers <- setdiff(blast.out$Alternate_SNP_ID, top_hits$Alternate_SNP_ID)
      remove_namescombined <- append(remove_namescombined,remove_namessimmarkers)      
       

      keep_namessimmarkers <- intersect(top_hits$Alternate_SNP_ID, blast.out$Alternate_SNP_ID)
      keep_namescombined <- append(keep_namescombined, keep_namessimmarkers)


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
  for ( i in c(1:3)) {
	pwld_link <- pwldmappings |>
	  dplyr::transmute(
	    SNP_A,
	    SNP_B,
	    w = .data[[score_col]]  # use chosen score column (R2 or R)
	  ) |>
	  dplyr::filter(
	    !base::is.na(SNP_A),
	    !base::is.na(SNP_B),
	    !base::is.na(w)
	  )         
        ## Keep only neighbours (SNP_B) that map to a single BLAST site
        b_counts <- blast.out |>
          dplyr::count(qaccver, name = "n_hits")
        
        single_b <- b_counts |>
          dplyr::filter(n_hits == 1L) |>
          dplyr::select(qaccver)
        

         
        ## SNP_B (chr, bp) for single-hit neighbours
        b_hits <- blast.out |>
          dplyr::semi_join(single_b, by = "qaccver") |>
          dplyr::select(qaccver, saccver, SNPpos)
          ## insert the single b add here for iterative loop. Apply end added markers into single_b and then deduplicate. simple statement check if iter is >1 and if it is add markers identified as newly unique to SNP b allowed for next iter 
        if (i>1) {
          b_hits_add_prioritr <- as.data.frame(multi_with_prior[,1:3])
          b_hits <- rbind(b_hits,b_hits_add_prioritr)

        }

        ## Link SNP_A with its usable neighbours (SNP_B single-hit),
        ## keep top 10 by weight per SNP_A (in intermediate I used 5 but here it is 10 incase of threat of translocation issues)
        neigh10 <- pwld_link |>
          dplyr::inner_join(b_hits, by = c("SNP_B" = "qaccver")) |>
          dplyr::rename(neigh_chr = saccver, neigh_bp = SNPpos) |>
          dplyr::group_by(SNP_A) |>
          dplyr::arrange(dplyr::desc(w), .by_group = TRUE) |>
          dplyr::slice_head(n = 10) |>
          dplyr::ungroup()

        prior_chr <- neigh10 |>
          dplyr::group_by(SNP_A, neigh_chr) |>
          dplyr::summarise(
            w_sum = base::sum(w, na.rm = TRUE),
            n     = dplyr::n(),
            .groups = "drop_last"
          ) |>
          dplyr::arrange(dplyr::desc(w_sum), dplyr::desc(n), neigh_chr, .by_group = TRUE) |>
          dplyr::slice(1) |>
          dplyr::ungroup() |>
          dplyr::transmute(qaccver = SNP_A, prior_chr = neigh_chr)
        
        # For distance weighting, we only need neighbour (bp, weight) on the chosen prior chr
        neigh_on_prior <- neigh10 |>
          dplyr::inner_join(prior_chr, by = c("SNP_A" = "qaccver", "neigh_chr" = "prior_chr")) |>
          dplyr::select(qaccver = SNP_A, neigh_bp, w)

        ## Per-target BLAST counts
        t_counts <- blast.out |>
          dplyr::count(qaccver, name = "total_hits")
        
        ## Attach counts and prior chr to target hits
        tmp <- blast.out |>
          dplyr::left_join(t_counts, by = "qaccver") |>
          dplyr::left_join(prior_chr, by = "qaccver")
        
        ## Rule 1 & 6: markers with a single BLAST hit are kept outright
        winners_single <- tmp |>
          dplyr::filter(total_hits == 1L) |>
          dplyr::group_by(qaccver) |>
          dplyr::slice(1) |>
          dplyr::ungroup()
        
        # For markers with multiple hits AND a prior chr:
        # if zero hits on that chr  drop
        #  if one hit on that chr  keep it
        #  if >1 on that chr  compute weighted mean distance to neighbours on prior chr
        #                         dist = sum_i w_i * |SNPpos_hit - neigh_bp_i| / sum_i w_i
        #                         pick the hit with MIN dist; tie higher bitscore
        multi_with_prior <- tmp |>
          dplyr::filter(total_hits > 1L, !base::is.na(prior_chr)) |>
          dplyr::filter(saccver == prior_chr) |>
          dplyr::inner_join(neigh_on_prior, by = "qaccver") |>
          dplyr::mutate(d = w * base::abs(SNPpos - neigh_bp)) |>
          dplyr::group_by(qaccver, saccver, SNPpos, bitscore) |>
          dplyr::summarise(
            dist_wmean = base::sum(d, na.rm = TRUE) / base::sum(w, na.rm = TRUE),
            .groups = "drop_last"
          ) |>
          dplyr::ungroup() |>
          dplyr::group_by(qaccver) |>
          dplyr::arrange(dist_wmean, dplyr::desc(dplyr::coalesce(bitscore, -Inf))) |>
          dplyr::slice(1) |>
          dplyr::ungroup()
        
        ## Combine winners:
        ## - single-hit targets (regardless of neighbour info)
        ## - multi-hit targets with prior-chr winner as above
        ## (Markers with multiple hits and no prior chr are dropped.)
        top_hits <- linkage_strict_hits <- dplyr::bind_rows(
          winners_single |>
            dplyr::select(qaccver, saccver, SNPpos, bitscore, dplyr::everything()),
          multi_with_prior |>
            dplyr::left_join(
              tmp |> dplyr::select(qaccver, saccver, SNPpos, bitscore, dplyr::everything()),
              by = c("qaccver", "saccver", "SNPpos", "bitscore")
            )
        ) |>
          dplyr::distinct(qaccver, .keep_all = TRUE)
        
        
      }
      }
      else{
        
        top_hits <- blast_unique <- blast.out |>
          dplyr::add_count(qaccver, name = "n") |>
          dplyr::filter(n == 1) |>
          dplyr::select(-n)
      }
      
      remove_namesld <- setdiff(blast.out$Alternate_SNP_ID, top_hits$Alternate_SNP_ID)
      remove_namescombined <- append(remove_namescombined,remove_namesld)

      keep_namesld <- intersect(top_hits$Alternate_SNP_ID, blast.out$Alternate_SNP_ID)
      keep_namescombined <- append(keep_namescombined, keep_namesld)

    }
    
    # Here is where the intermitten deviates from the strict. Intermittent was a filter keep, here we want a remove all so we will tally all the names of failed markers
    # apply them, and take out of the dataframe completely. 
    
    keep_namescombined <- unique(keep_namescombined)

       
    #failed_reads <- blast.out |>
     # dplyr::filter((qaccver %in% remove_namescombined))
    
    #blast_unique <- blast.out |>
     # dplyr::filter(!(qaccver %in% remove_namescombined))
    

    failed_reads <- blast.out |>
      dplyr::filter(!(Alternate_SNP_ID %in% keep_namescombined ))
    
    blast_unique <- blast.out |>
      dplyr::filter((Alternate_SNP_ID %in% keep_namescombined))


    # Filter to have only SNPs present in blast_unique present in the mappings file
    filtered_mappings <- mappings.out |>
      dplyr::filter(Alternate_SNP_ID %in% blast_unique$Alternate_SNP_ID)
    
    
    attributes(blast_unique)$raw_header_lines <- headerinfo$raw_header_lines
    blast_unique_clean <- sanitize_rectangular(blast_unique)
    
    
    attributes(filtered_mappings)$raw_header_lines <- headerinfo$raw_header_lines
    filtered_mappings_clean <- sanitize_rectangular(filtered_mappings)
    

# orientation update file section

do_orient <- isTRUE(doorientation) ||
      tolower(as.character(doorientation)) == "yes"

    if (do_orient &&
        !is.null(orientation.file) &&
        !is.na(orientation.file) &&
        nzchar(orientation.file)) {

      if (!file.exists(orientation.file)) {
        warning("Orientation file not found: ", orientation.file,
                ". Skipping orientation update.")
      } else if (nrow(blast_unique) > 0) {

        orient_df <- utils::read.delim(
          orientation.file,
          header = TRUE,
          sep = "\t",
          stringsAsFactors = FALSE,
          check.names = FALSE
        )

        if (!all(c("ID", "Orientation") %in% colnames(orient_df))) {
          warning(
            "Orientation file is missing required columns 'ID' and 'Orientation'; ",
            "skipping orientation update."
          )
        } else {
          # Map from marker ID -> sstrand in final unique mappings
          # (using qaccver as the marker name, consistent with addSNPdetails)
          sstrand_map <- stats::setNames(
            tolower(as.character(blast_unique$sstrand)),
            blast_unique$qaccver
          )

          # Orientation current values
          ori_lower <- tolower(as.character(orient_df$Orientation))
          marker_sstrand <- sstrand_map[orient_df$ID]

          # Update only where:
          #  - Orientation is "unknown" (case-insensitive)
          #  - We have a unique mapping strand for that ID
          #  - Strand is plus/minus
          to_update <- !is.na(marker_sstrand) &
            ori_lower == "unknown" &
            marker_sstrand %in% c("plus", "minus")

          if (any(to_update)) {
            orient_df$Orientation[to_update] <- marker_sstrand[to_update]
          }

          # Decide where to write the updated file
          proj_dir <- if (!is.null(project.dir) && nzchar(project.dir)) {
            project.dir
          } else {
            output.path
          }

          if (!dir.exists(proj_dir)) {
            dir.create(proj_dir, recursive = TRUE, showWarnings = FALSE)
          }

          out_orient <- file.path(proj_dir, "Updated_orientation_file.tsv")

          # Overwrite allowed (default behaviour of write.table)
          utils::write.table(
            orient_df,
            file      = out_orient,
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE
          )
        }
      }
    }
    
    
    
    hits.outpath <- file.path(paste0(probe.name, "_with_", genome.name, "_strict_filtering_hits.csv"))
    
    
    
    write_with_header_base(df = blast_unique_clean,
                           path = hits.outpath,sep=",",quote=FALSE)
    
    
    mappings.outpath <- file.path(paste0(probe.name, "_with_", genome.name, "_strict_filtering_mappings.tsv"))
    
    write_with_header_base(df = filtered_mappings_clean,
                           path = mappings.outpath,sep="\t",quote=FALSE)
    
    
    cat("\nNumber of rows in final filtering results blast table:", nrow(blast_unique_clean), "\n")
    
  }