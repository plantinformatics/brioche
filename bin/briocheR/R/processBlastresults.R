# Define the State enumeration (using a named vector for simplicity)
State <- c(NONE = 0, QUERY_GAP = 1, TARGET_GAP = 2, MATCH = 3)
# Main parsing function

#' @title Parse BTOP from BLAST
#' @author David Chisanga
#' @description Function to separate variants within the btop string outpt by blast
#' @param btop String btop value from blast
#' @examples
#' parse_btop("4A-40-AGC25TA5")
#' @keywords btop,blast
#' @export
parse_btop <- function(btop) {
  target_coordinates <- 0  # Initialize with 0
  query_coordinates <- 0
  state <- State["NONE"]

  # Extract tokens using regular expressions
  tokens <- regmatches(btop, gregexpr("([A-Z-]*{2}|\\d+)", btop))[[1]]

  for (token in tokens) {
    if (startsWith(token, "-")) {
      if (state != State["QUERY_GAP"]) {
        target_coordinates <- c(target_coordinates, tail(target_coordinates, 1))
        query_coordinates <- c(query_coordinates, tail(query_coordinates, 1))
        state <- State["QUERY_GAP"]
      }
      target_coordinates[length(target_coordinates)] <- target_coordinates[length(target_coordinates)] + 1
    } else if (endsWith(token, "-")) {
      if (state != State["TARGET_GAP"]) {
        target_coordinates <- c(target_coordinates, tail(target_coordinates, 1))
        query_coordinates <- c(query_coordinates, tail(query_coordinates, 1))
        state <- State["TARGET_GAP"]
      }
      query_coordinates[length(query_coordinates)] <- query_coordinates[length(query_coordinates)] + 1
    } else {
      length <- suppressWarnings(as.numeric(token))  # Handle potential non-numeric tokens
      if (is.na(length)) {
        length <- 1
      }

      if (state == State["MATCH"]) {
        target_coordinates[length(target_coordinates)] <- target_coordinates[length(target_coordinates)] + length
        query_coordinates[length(query_coordinates)] <- query_coordinates[length(query_coordinates)] + length
      } else {
        target_coordinates <- c(target_coordinates, tail(target_coordinates, 1) + length)
        query_coordinates <- c(query_coordinates, tail(query_coordinates, 1) + length)
        state <- State["MATCH"]
      }
    }
  }

  coordinates <- rbind(target_coordinates, query_coordinates)
  return(coordinates)
}


#' @title Check input data
#' @description This function is used to process and check for input details
#' @author David Chisanga
#' @param blast.file path to the blast output file
#' @keywords blast, alignment
#'
#' @return Returns a named list
#' @noRd
getTarget <- function(targets) {
  # Check if targets file exists
  if (!file.exists(targets)) {
    stop(paste0("Targets file was not found at ", targets))
  }

  # Read data based on file type, prioritize RDS for efficiency
  if (grepl("RDS$", targets, ignore.case = TRUE)) {
    targets <- readRDS(targets)
  } else {
    targets <- read.table(
      targets,
      stringsAsFactors = FALSE, # Avoid factors for efficiency
      header = FALSE,
      strip.white = TRUE,
      fill = TRUE,
      comment.char = "" # Skip commented lines if present
    )

    # Handle header if present
    if ("ID" %in% targets[1,]) {
      colnames(targets) <- targets[1,]
      targets <- targets[-1,]
    } else {
      colnames(targets) <- c("ID", "Target.Chr", "Target.bp", "Targeted.SNP")
    }
  }

  # Check for data and remove duplicates
  if (nrow(targets) == 0) {
    stop("Targets file provided has null content, check file and try again")
  }

  targets <- unique(targets) # More efficient than !duplicated
  return(targets)
}


#' @title Function to check if a mapping is extendable
#' @author David Chisanga
#' @details This function is used to check if a mapping is extendable
#' @param istarget3primend Logic to check if a chip used is 3 prime-end
#' @param extandable.site.bps Integer with the number of allowed base pairs
#' @param min.matchbp Integer specifying the number of minimum
#' @param btop
#' @noRd
checkIfExtendable <- function(istarget3primeend,
                              extendable.site.bps,
                              min.matchbp,
                              target.bp,
                              marker.char,
                              btop) {
  # Extract matching bases from 3' end using regex
  matching_bases <- sub(".*_",
                        "",
                        gsub(
                          "([A-Z][A-Z])|([A-Z]-)|(-[A-Z])",
                          "_",
                          strsplit(btop, marker.char)[[1]][1]
                        ))

  return(matching_bases)
}

#'@title Check hybridisation status of each hit
#'@description Function used to check the hybridisation status of a probe sequence
#' @param marker.data Row list of marker data
#' @param variants list of variants and their positions
#' @param min.coverage minimum coverage as percentage
#' @param min.pident minimum percentage identity
#' @param istarget3primeend Boolean value showing whether the target marker is on the 3 prime end
#' @param min.matchbp Integer showing the number of expected correctly matched base pairs
#' @author David Chisanga

checkHybridisation <- function(coverage,
                               pident,
                               var_type,
                               variants,
                               var.position,
                               istarget3primeend = FALSE,
                               min.coverage,
                               min.pident,
                               min.matchbp,
                               marker.char) {
  # Initialize hybridization status
  hybridize <- 'No'

  if (!istarget3primeend) {
    # Check coverage and percent identity
    if (coverage >= min.coverage &&
        pident >= min.pident) {

      # Find marker index
      marker.index <- grep(marker.char, variants)
      if (length(marker.index) > 0) {
        # Check for gaps before and after marker
        b4var <- var_type[marker.index - 1]
        aftervar <- var_type[marker.index + 1]

        if (!is.na(b4var) && !is.na(aftervar)) {
          if (b4var == "M" && aftervar == "M") {
            hybridize <- "Yes"
          } else {
            hybridize <- "Maybe"
          }
        }
      }
    }
  } else {
    # Check for 3' end hybridization
    if (marker.data[["qend"]] == marker.data[["qlen"]] &&
        marker.data[["length"]] >= min.matchbp) {
      hybridize <- "Yes"
    }
  }

  return(hybridize)
}

#' @title Get Variants from BLAST BTOP
#' @author David Chisanga
#' @description Function to separate variants within the btop string outpt by blast
#' @param btop String btop value from blast
#' @examples
#' getVariantsFromBTOP("4A-40-AGC25TA5")
#' @keywords btop,variants,blast
#' @export
getVariantsFromBTOP <- function(btop) {
  prevchar <- NA
  i <- 1
  #replace numbers of matching bases with M
  btop <- gsub("[0-9]+", "M", btop)
  bpletters <- LETTERS[LETTERS != "M"]
  variants <- NULL
  for (char in strsplit(as.character(btop), "")[[1]]) {
    if (char == "M")
    {
      variants[i] <- char
      i <- i + 1
    }
    else if (is.na(prevchar)) {
      prevchar <- char
    }
    else{
      if (prevchar %in% bpletters & char %in% bpletters) {
        variants[i] <- paste0(char, ">", prevchar)
        i <- i + 1
        prevchar <- NA
      } else if ((prevchar %in% bpletters &
                  char == "-") || (char %in% bpletters & prevchar == "-")) {
        variants[i] <- paste0(prevchar, char)
        i <- i + 1
        prevchar <- NA
      }
      else {
        prevchar <- char
      }
    }
  }
  variant.type <- variants
  variant.type[grepl(">", variants)] <- "SNP"
  variant.type[grepl("-", variants)] <- "G"
  return(list(vars = variants, type = variant.type))
}


#' @title Function to get variant positions from blast btop
#' @author David Chisanga
#' @description The function is used to obtain the positions of variants
#' from the btop output from blast
#' @param btop String btop value
#' @keywords btop,blast,variant,position
#' @examples
#' getPositionsFromBTOP("4A-40-AGC25TA5")
#' @export
getPositionsFromBTOP <- function(btop,qstart)
{
  #get variant positions
  patt <- "_1_"
  positions2print <- ""
  tryCatch({
    prevchar <- NA
    positions <- NULL
    nums <- as.character(0:9)
    for (char in strsplit(as.character(btop), "")[[1]])
    {
      if (is.na(prevchar))
      {
        if (char %in% nums)
        {
          positions <- paste0(positions, char)
        }
        else{
          prevchar <- char
        }
      }
      else if ((prevchar %in% LETTERS &
                char %in% LETTERS) ||
               (prevchar == "-" &
                char %in% LETTERS) || (prevchar %in% LETTERS & char == "-"))
      {
        positions <- paste0(positions, patt)
        prevchar <- NA
      }
      else if (prevchar == "-" & char == "-")
      {
        positions <- paste0(positions, patt)
        prevchar <- NA
      }
      else{
        stop(paste0("No match found for '", prevchar, "' and '", char, "'"))
      }
    }

    positions <- strsplit(gsub("(_)+", "_", positions), "_")[[1]]

    if (length(positions) > 0)
      return(cumsum(as.numeric(positions))-as.numeric(qstart))
    else
    {
      return(".")
    }

  },
  error = function(e) {
    stop(e)
  },
  warning = function(w) {
    stop(w)
  },
  finally = function(f) {
    print(f)
  })
}



#' @title Process alignment results from blast
#' @description This function is used to process and format the blast output
#' @author David Chisanga
#' @param blast.file path to the blast output file
#' @param min.matchbp minimum length of HSP to be considered fully Hybridize
#' @param target table with target annotation
#' @param outfmt vector of format specifiers from the supported format specifiers for 6,7 and 10 in blastn's 'outfmt' parameter
#' @param path.2save.coords path where SNP coordinates to get base in reference genome using samtools are to be saved
#' @param extendable.site.bps number of matching base pairs from the 3 prime end for a probe to be considered as extendable
#' @param istarget3primeend logical value showing if target is on the 3' prime end
#' @param num.threads integer value of the number of threads passed to reading the data file
#' @param keep.duplicates logical value, whether to retain lines with multiple hits
#' @keywords blast, alignment
#' @export
#'
#' @return Returns a data frame of processed blast results
#'
#'
processBlastResults <- function(blast.file,
                                min.length = 40,
                                path.2save.coords = "SNPcoordinates",
                                extendable.site.bps = 3,
                                targets,
                                outfmt = c(
                                  "qaccver",
                                  "saccver",
                                  "pident",
                                  "qlen",
                                  "length",
                                  "mismatch",
                                  "gapopen",
                                  "qstart",
                                  "qend",
                                  "sstart",
                                  "send",
                                  "evalue",
                                  "bitscore",
                                  "btop",
                                  "qseq",
                                  "sseq",
                                  "sstrand"
                                ),
                                istarget3primeend,
                                marker.char = "N",
                                min.coverage = 70,
                                min.pident = 90) {
  tryCatch({
    message("Marker passed is ", marker.char) # Use message instead of print

    # Input Validation
    stopifnot(file.exists(blast.file), file.exists(targets))
    if (missing(istarget3primeend)) {
      stop("argument 'istarget3primeend' is missing!")
    }
    istarget3primeend <- toupper(as.character(istarget3primeend)) == "TRUE"

    # Load Data
    targets <- getTarget(targets)
    data <- data.table::fread(
      input = blast.file,
      strip.white = TRUE,
      fill = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      showProgress = TRUE
    )
    colnames(data) <- outfmt

    # Calculate Coverage
    data$Coverage <- 100 * data$length / (data$qlen-1) #Substract 1 to account for the extra base added to the query

    # Add Target Information
    data <- merge(data,
                  targets,
                  by.x = "qaccver",
                  by.y = "ID",
                  all.x = TRUE)
    # Extract SNP Clusters - using stringr for better readability
snp_col <- if ("Target.base" %in% names(data)) "Target.base" else
           if ("Targeted.SNP" %in% names(data)) "Targeted.SNP" else
           stop("Neither 'Target.base' nor 'Targeted.SNP' present after merge().")

# Use fully qualified calls ? no library() needed
data$ClusterA_NT <- stringr::str_remove_all(as.character(data[[snp_col]]), "\\[|\\]|/[ACTG]")
data$ClusterB_NT <- stringr::str_remove_all(as.character(data[[snp_col]]), "[ACTG]/|\\[|\\]")    # Extract SNP Clusters
    # data$ClusterA_NT <- gsub("\\[|\\]|/[ACTG]",
    #                          "",
    #                          data$Target.base,
    #                          perl = TRUE,
    #                          ignore.case = TRUE)
    # data$ClusterB_NT <- gsub("[ACTG]/|\\[|\\]",
    #                          "",
    #                          data$Target.base,
    #                          perl = TRUE,
    #                          ignore.case = TRUE)

    # Calculate SNP Position
# Map a single row's Target.bp (query coordinate) to subject coordinate
calculateSNPpos <- function(qstart, qend, sstart, send, sstrand, qseq, sseq, target_bp) {
  qstart   <- as.integer(qstart)
  qend     <- as.integer(qend)
  sstart   <- as.integer(sstart)
  send     <- as.integer(send)
  target_bp<- as.integer(target_bp)
  
  # normalise strand check (accepts "plus"/"minus" any case)
  is_plus  <- !is.na(sstrand) && grepl("^\\s*plus\\s*$", sstrand, ignore.case = TRUE)
  
  # subject "first aligned" genomic coordinate depends on strand/order
  s_first  <- if (is_plus) min(sstart, send) else max(sstart, send)
  
  # k = 1-based index of the query base *within this HSP* (counting only non-gaps in qseq)
  k <- target_bp - qstart + 1L
  
  # Precompute non-gap columns in qseq
  q_cols_m <- stringi::stri_locate_all_regex(qseq, "[^-]")[[1L]]
  q_cols   <- if (is.null(dim(q_cols_m))) integer(0) else q_cols_m[,1L]
  q_len    <- length(q_cols)  # number of aligned query bases in this HSP
  
  # Fast boundary cases relative to the alignment on the query:
  if (k == 0L) {
    # SNP is immediately BEFORE qstart
    return(if (is_plus) s_first - 1L else s_first + 1L)
  }
  if (k == q_len + 1L) {
    # SNP is immediately AFTER qend
    s_len <- stringi::stri_count_regex(sseq, "[^-]")  # aligned subject length (non-gaps)
    return(if (is_plus) s_first + s_len else s_first - s_len)
  }
  
  # If Target.bp lies outside this HSP (beyond the ±1 boundary), we can't map it here. This is strict, maybe its best to incorporate something else her
  # something that conveys some information about SNP position approximation. 
  if (k < 1L || k > q_len) return(NA_integer_)
  
  # Alignment column of the k-th non-gap in qseq
  col_idx <- q_cols[k]
  
  # Subject bases consumed up to that column (non-gaps in sseq up to col_idx)
  s_bases <- stringi::stri_count_regex(stringi::stri_sub(sseq, 1L, col_idx), "[^-]")
  
  # Convert to genomic subject coordinate
  if (is_plus) {
    pos <- s_first + s_bases - 1L
  } else {
    pos <- s_first - s_bases + 1L
  }
  
  # target bp is found somewhere not on the periphery of the marker. Need to count the number of gaps within a target region. for qseq, adjust sseq region by gaps
  # Count sseq region gaps and apply q and s gaps to absolute reference position. I can't test this with the avaliable data so will need to come back to it with new 
  # dataset
  if (k > 1L && k < q_len) {
    # alignment column of the k-th non-gap in qseq
    col_idx <- stringi::stri_locate_all_regex(qseq, "[^-]")[[1L]][,1L][k]
    
    # gaps up to that column
    q_gaps <- stringi::stri_count_fixed(stringi::stri_sub(qseq, 1L, col_idx), "-")
    s_gaps <- stringi::stri_count_fixed(stringi::stri_sub(sseq, 1L, col_idx), "-")
    
    # use offset inside HSP: (k - 1L)
    if (is_plus) {
      pos <- s_first + (k - 1L) + q_gaps - s_gaps
    } else {
      pos <- s_first - (k - 1L) - q_gaps + s_gaps
    }
  }
  
  pos
}

# Vectorized across rows; creates/overwrites data$SNPpos
data$SNPpos <- mapply(
  calculateSNPpos,
  qstart    = data$qstart,
  qend      = data$qend,
  sstart    = data$sstart,
  send      = data$send,
  sstrand   = data$sstrand,
  qseq      = data$qseq,
  sseq      = data$sseq,
  target_bp = data$Target.bp,
  USE.NAMES = FALSE
)

    cat("Finished calculating SNP position\n",
        file = "logger.log",
        append = TRUE)


    # Process Variants and Gaps - pre-allocate vectors
    SNPorGap <- SNPorGap_pos <- SNPorGAPpos_Ref <- Variant <- SNPpos_all <-
      Identical_bps_from_target <- vector("character", nrow(data))

    for (i in 1:nrow(data)) {
      variantinfo <- getVariantsFromBTOP(data$btop[i])
      var_position <- getPositionsFromBTOP(data$btop[i],data$qstart[i])
      var_type <- variantinfo$type
      variants <- variantinfo$vars

      # Filter out matches ("M")
      keep <- var_type != "M"
      var_position <- var_position[keep]
      var_type <- var_type[keep]
      variants <- variants[keep]

      var_type <- paste0(var_type, collapse = ";")
      SNPorGAPpos_Ref[i] <- "."
      SNPpos_all[i] <- "."
      # If target is not 3 prime end SNPpos is removed? Not sure why this is done and whether the back end non 3prime is as robust to redetect
      if (!istarget3primeend)
        data$SNPpos[i] <- "."
      # detect position of other variant sites identified using getVariantsFromBTOP. Based on their position in the the sequence simply add or
      # subtract ssttart to get the exact location of variant (plus vs minus strand)
      if (length(var_position) > 0) {
        if (class(var_position) != "character") {
          compute_posN_for_row <- function(i, var_position, variants, data) {
            # normalize inputs
            vpos <- as.integer(var_position)
            vtyp <- toupper(trimws(variants))
            strand <- tolower(trimws(data$sstrand[i]))
            sstart <- as.integer(data$sstart[i])
            send   <- as.integer(data$send[i])
            
            # anchor = first aligned subject coordinate for this HSP on the genome
            is_plus <- identical(strand, "plus")
            s_first <- if (is_plus) min(sstart, send) else max(sstart, send)
            
            # identify sseq gaps (only A/C/G/T followed by '-')
            is_sseq_gap <- grepl("^[ACGT]-$", vtyp)
            
            # alignment columns where sseq has gaps; duplicates allowed (consecutive gaps)
            gap_cols <- sort(vpos[is_sseq_gap])
            
            # count number of prior sseq gaps strictly before each variant's column
            prior_gap_count <- if (length(gap_cols)) {
              # use -0.5 to implement "strictly less than y" for integer columns
              findInterval(vpos - 0.5, gap_cols)
            } else {
              integer(length(vpos))
            }
            
            # corrected subject genomic positions (vector)
            if (is_plus) {
              s_first + vpos - prior_gap_count
            } else {
              s_first - vpos + prior_gap_count
            }
          }
          posN <- compute_posN_for_row(i, var_position, variants, data)
          # Do a check. If one of the variants is marked as a SNP with the marker character (in situations where the SNP is present part way
          # into the probe sequence not at the ends) then take the absolute position of that SNP and replace the earlier SNPpos
          # Note this is unnessary after I fixed the above command for finding the SNP position anyway. 
          targetpos <- posN[grepl(marker.char, variants)]
          if (length(targetpos) > 0)
            data$SNPpos[i] <- targetpos

          SNPorGAPpos_Ref[i] <- paste0(sapply(posN, function(xx) {
            if (is.na(xx))
              return(".")
            paste0(data$saccver[i], ":", xx, "-", xx)
          }), collapse = ";")

          if (!istarget3primeend)
            SNPpos_all[i] <- paste0(posN, collapse = ";")
        }
      }

      if (var_type == "")
        var_type <- "."
      SNPorGap[i] <- var_type
      SNPorGap_pos[i] <- paste0(var_position, collapse = ";")
      Variant[i] <- paste0(variants, collapse = ";")

      cat("BTOP for index",i,"is",data$btop[i],"\n",file="logger.log")
      # Identical_bps_from_target estimates the number of correct bases from the SNP of interest (target SNP).
      Identical_bps_from_target[i] <- gsub(".*_",
                                           "",
                                           gsub(
                                             "([A-Z][A-Z])|([A-Z]-)|(-[A-Z])",
                                             "_",
                                             strsplit(as.character(data$btop[i]), as.character(marker.char))[[1]][1]
                                           ))

    }

    # Add Processed Columns
    data$SNPorGap <- SNPorGap
    data$SNPorGap_pos <- SNPorGap_pos
    data$SNPorGAPpos_Ref <- SNPorGAPpos_Ref
    data$Variant <- Variant
    data$SNP.Refpos <- paste0(data$saccver, ":", data$SNPpos, "-", data$SNPpos)
    # Simple if else command it is asking, is the stretch of identical bases from target larger than the minimum set extension length
    # e.g., in this case larger than 3
    if (istarget3primeend) {
      qend_matches <- data$qend == (as.numeric(data$Target.bp) - 1)
      data$Identical_bps_from_target <- Identical_bps_from_target
      data$ExtendableSite <- ifelse(
        as.numeric(data$Identical_bps_from_target) >= extendable.site.bps &
          qend_matches,
        "Yes",
        "."
      )
      # This is making multiple or/ or+ and statements to generate the hybridized label. 
      data$Hybridized <- ifelse((abs(data$qend - data$qlen) == 1 |
                                   data$qend == data$qlen) &
                                  data$length >= min.length &
                                  data$gapopen == 0 &
                                  data$ExtendableSite == "Yes",
                                "Yes",
                                "No"
      )
    }

    # Filter and Arrange - using dplyr for better readability
    # May need to revise this bit here as We may want to report na values with the rows. Take it out temporarily # & !is.na(SNPpos)
    data <- data |>
      dplyr::filter(SNP.Refpos != ".") |>
      dplyr::arrange(dplyr::desc(Coverage), dplyr::desc(pident), qaccver)

    # Save Results
    if (!dir.exists(path.2save.coords))
      dir.create(path.2save.coords, recursive = TRUE)
    if ("Target.Chr" %in% colnames(data))
      data$Target.Chr <- NULL
    saveRDS(data,
            file.path(path.2save.coords, "processed_blast_results.RDS"))
 # Added in a filter for NA for markers for when the reference alignment wasn't complete.
    write.table(
      subset(data, SNP.Refpos != paste0(saccver, ":.-.") & SNP.Refpos != paste0(saccver,":NA-NA"))$SNP.Refpos,
      file = file.path(path.2save.coords, "samtools_marker_positions.txt"),
      quote = FALSE,
      sep = "\t",
      col.names = FALSE,
      row.names = FALSE
    )
    return(data)
  },
  error = function(e) {
    stop(e)
  },
  warning = function(w) {
    warning(w)
  })
}