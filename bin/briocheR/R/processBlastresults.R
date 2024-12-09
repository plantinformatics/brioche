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
    if(!require(stringr)) {
      install.packages("stringr")
      library(stringr)
    }

    data$ClusterA_NT <- str_remove_all(data$Target.base, "\\[|\\]|/[ACTG]")
    data$ClusterB_NT <- str_remove_all(data$Target.base, "[ACTG]/|\\[|\\]")
    # Extract SNP Clusters
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
    calculateSNPpos <- function(x) {
      targetbp <- as.numeric(x[["Target.bp"]]) - as.numeric(x[["qstart"]])
      pos <- as.numeric(x[["sstart"]]) + targetbp
      if (x[["sstrand"]] == "minus")
        pos <- as.numeric(x[["sstart"]]) - targetbp
      return(pos)
    }
    data$SNPpos <- sapply(1:nrow(data), function(i)
      calculateSNPpos(data[i, ]))
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

      if (!istarget3primeend)
        data$SNPpos[i] <- "."

      if (length(var_position) > 0) {
        if (class(var_position) != "character") {
          posN <- sapply(var_position, function(y) {
            ifelse(
              data$sstrand[i] == "plus",
              as.numeric(data$sstart[i]) + as.numeric(y),
              as.numeric(data$sstart[i]) - as.numeric(y)
            )
          })

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

    if (istarget3primeend) {
      qend_matches <- data$qend == (as.numeric(data$Target.bp) - 1)
      data$Identical_bps_from_target <- Identical_bps_from_target
      data$ExtendableSite <- ifelse(
        as.numeric(data$Identical_bps_from_target) >= extendable.site.bps &
          qend_matches,
        "Yes",
        "."
      )
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
    if(!require(dplyr)) {
      install.packages("dplyr")
      library(dplyr)
    }
    data <- data %>%
      filter(SNP.Refpos != "." & !is.na(SNPpos)) %>%
      arrange(desc(Coverage), desc(pident), qaccver)

    # Save Results
    if (!dir.exists(path.2save.coords))
      dir.create(path.2save.coords, recursive = TRUE)
    if ("Target.Chr" %in% colnames(data))
      data$Target.Chr <- NULL
    saveRDS(data,
            file.path(path.2save.coords, "processed_blast_results.RDS"))

    write.table(
      subset(data, SNP.Refpos != paste0(saccver, ":.-."))$SNP.Refpos,
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


