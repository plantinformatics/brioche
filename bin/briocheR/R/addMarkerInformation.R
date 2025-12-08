#' @title Get the alternative base for a probe ID
#' @description This function used is to get the alternative base pair from blast results
#' @author David Chisanga
#' @param NT_A Design cluster A nucleotide
#' @param NT_B Design cluster B nucleotide
#' @param Ref_NT Nucleotide at SNP position from reference genome
#' @keywords genome, alternative allele, reference allele
#' @export
#'
#' @return Returns a character string indicating the alternative base pair
#'
getAltBP <- function(NT_A, NT_B, Ref_NT)
{
  NT_A<-toupper(trimws(NT_A))
  NT_B<-toupper(trimws(NT_B))
  Ref_NT<-toupper(trimws(Ref_NT))

  map_table <-
    c(
      "TACG",
      "TAGC",
      "ATGC",
      "ATCG",
      "CGAT",
      "CGTA",
      "GCAT",
      "GCTA",
      "CATG",
      "GATC",
      "GTAC",
      "CTAG",
      "AGCT",
      "TGCA",
      "ACGT",
      "TCGA")

  map_tableshort <-
    c(
      "TAC",
      "TAG",
      "ATG",
      "ATC",
      "CGA",
      "CGT",
      "GCA",
      "GCT",
      "CAT",
      "GAT",
      "GTA",
      "CTA",
      "AGC",
      "TGC",
      "ACG",
      "TCG")
  sapply(seq_along(NT_A), function(i) {
    out <- ifelse(NT_A[i] == Ref_NT[i],
                  NT_B[i],
                  ifelse(Ref_NT[i] == NT_B[i], NT_A[i],
                         map_table[match(paste0(NT_A[i], Ref_NT[i], NT_B[i]),
                                         map_tableshort)]))
    out <- substr(out, nchar(out), nchar(out))
    if (is.na(out))
      out <- "."
    return(out)
  }, simplify = TRUE, USE.NAMES = FALSE)
}

#' @title Check if target
#' @description This function is used to add target details to the blast output returned by the `processBlastResults()` function
#' @author David Chisanga
#' @param target.chrom Target chromosome obtained from the probe design matrix
#' @param extendable String showing whether the target SNP called is extendable
#' @param SNP.position Position of SNP called by samtools
#' @param target.bp The target base pair obtained from the probe design matrix
#' @param ref.chrom The chromosome on the reference genome as reported by blast
#' @keywords target, genome,
#'
#' @return Returns a character vector showing the target results
#'
#' @noRd

checkTarget <-
  function(sstart,
           ssend,
           strand,
           variants,
           var.pos,
           extendable,
           SNP.position,
           minmismatch)
  {
    #=IF(extendable!="Yes","n/a",IF(LEN(AY4)<>2,"unknown",IF(AND(AY4=RIGHT(E4,2),AH4=AM4),"target","off-target")))
    sapply(1:length(target.chrom), function(x) {
      out <- "off-target"
      chr<-gsub("^chr","",target.chrom[x],ignore.case = T)
      ref.chr<-gsub("^chr","",ref.chrom[x],ignore.case = T)
      if (extendable[x] != "Yes")
        out <- "n/a"
      else if (chr==ref.chr & SNP.position[x] == target.bp[x])
        out <- "target"
      else
        out <- "unknown"
      return(out)
    }, simplify = T, USE.NAMES = F)
  }



#' @title Add SNP details to the blast output
#' @description This function is used to add more details to the blast output
#'   returned by the `processBlastResults()` function
#' @author David Chisanga
#' @param blast.path Dataframe or path to blast results returned by the
#'   `processBlastResults()` function (RDS file).
#' @param reference.bases Path to the reference bases obtained from samtools.
#' @param probe.name Name of the probe file being analysed.
#' @param genome.name Name of the reference genome being compared.
#' @param output.path Path where the results are saved.
#' @param pident.filter Double value to filter hits with small percentage identity.
#' @param coverage.filter Double value to filter hits with small coverage.
#' @param is3prime Logical. Does the SNP exist as 3' as part of SNP chip marker mapping.
#' @param doorientation Logical or character ("Yes"/"No"). If TRUE or "Yes"
#'   (case-insensitive), marker orientation harmonisation using `orientation.file`
#'   will be applied.
#' @param isbiallelic Logical or character ("Yes"/"No"). If TRUE or "Yes"
#'   (case-insensitive), REF/ALT are re-oriented to be consistent with the
#'   Cluster A/B nucleotides (using reverse complement if needed).
#' @param orientation.file Optional path to a TSV file describing the current state
#'   of marker orientation. Must contain columns `ID` and `Orientation`.
#' @keywords genome, alternative allele, reference allele
#' @export
#'
#' @return Returns a dataframe with the updated blast results. The returned object
#'   has attributes `is3prime`, `doorientation`, `isbiallelic`, and
#'   `orientation.file` attached for downstream use.
#'
addSNPdetails <-
  function(blast.path,
           reference.bases,
           probe.name,
           genome.name,
           output.path = "Mapping_results",
           pident.filter = 90,
           coverage.filter = 50,
           is3prime,
           doorientation = FALSE,
           isbiallelic = FALSE,
           orientation.file = NULL)
  {
    if (!dir.exists(output.path))
      dir.create(output.path, recursive = TRUE)

    ref <- utils::read.delim(reference.bases, header = FALSE)

    #
    if (inherits(blast.path, c("data.frame", "matrix", "data.table"))) {
      blast.out <- as.data.frame(blast.path)
    } else {
      if (!file.exists(blast.path)) {
        stop("Blast results should be a dataframe or path to an RDS file")
      }
      blast.out <- readRDS(blast.path)
    }

    #Attach orientation / biallelic settings as attributes for downstream use
    attr(blast.out, "is3prime")         <- is3prime
    attr(blast.out, "doorientation")    <- doorientation
    attr(blast.out, "isbiallelic")      <- isbiallelic
    attr(blast.out, "orientation.file") <- orientation.file

    #Add Ref coordinates to data frame
    blast.out$Ref <- ref$V2[match(blast.out$SNP.Refpos, ref$V1)]
    blast.out$Ref[is.na(blast.out$Ref)] <- "."

    if ("ClusterA_NT" %in% colnames(blast.out)) {

      #3prime design
      if (isTRUE(is3prime)) {
        blast.out$ALT <- with(
          blast.out,
          getAltBP(
            NT_A   = ClusterA_NT,
            NT_B   = ClusterB_NT,
            Ref_NT = Ref
          )
        )
      }

      # Non 3prime design
      if (!isTRUE(is3prime)) {

        # normalise inputs
        ss   <- tolower(as.character(blast.out$sstrand))   # "plus" or "minus"
        refb <- toupper(as.character(blast.out$Ref))
        a_nt <- toupper(as.character(blast.out$ClusterA_NT))
        b_nt <- toupper(as.character(blast.out$ClusterB_NT))

        # convert to plus-strand bases (complement if on minus)
        plusA <- ifelse(ss == "minus",
                        chartr("ACGTacgt", "TGCAtgca", a_nt),
                        a_nt)
        plusB <- ifelse(ss == "minus",
                        chartr("ACGTacgt", "TGCAtgca", b_nt),
                        b_nt)

        # set ALT: if one equals Ref, ALT is the other; else "[A/B]" using plus-strand bases
        blast.out$ALT <- ifelse(
          is.na(plusA) | is.na(plusB) | is.na(refb),
          NA_character_,
          ifelse(
            plusA == refb, plusB,
            ifelse(
              plusB == refb, plusA,
              paste0("[", plusA, "/", plusB, "]")
            )
          )
        )
      }

      ## Orientation harmonisation using orientation.file if set by user 
      do_orient <- isTRUE(doorientation) ||
        tolower(as.character(doorientation)) == "yes"

      if (do_orient &&
          !is.null(orientation.file) &&
          !is.na(orientation.file) &&
          nzchar(orientation.file)) {

        if (!file.exists(orientation.file)) {
          warning("orientation.file does not exist: ", orientation.file)
        } else {
          # Expect columns: ID, Orientation
          orient_df <- utils::read.delim(
            orientation.file,
            header = TRUE,
            stringsAsFactors = FALSE
          )

          if (!all(c("ID", "Orientation") %in% colnames(orient_df))) {
            warning(
              "orientation.file does not have required columns 'ID' and 'Orientation'; ",
              "skipping orientation harmonisation."
            )
          } else {
            # Named vector: marker ID -> orientation
            ori_map <- orient_df$Orientation
            names(ori_map) <- orient_df$ID

            marker_ids <- as.character(blast.out$qaccver)
            marker_ori <- tolower(as.character(ori_map[marker_ids]))
            sstrand    <- tolower(as.character(blast.out$sstrand))

            # we only care about plus/minus vs plus/minus; unknown = no change
            valid_strand <- sstrand %in% c("plus", "minus")
            known_ori    <- marker_ori %in% c("plus", "minus")

            # Different orientation and orientation is known (not 'unknown')
            mismatch_known <-
              !is.na(marker_ori) &
              !is.na(sstrand) &
              valid_strand &
              known_ori &
              (marker_ori != sstrand)

            if (any(mismatch_known)) {
              # Complement bases A<->T, C<->G in REF and ALT;
              # punctuation/brackets etc are preserved.
              comp_vec <- function(x) {
                x <- as.character(x)
                na_idx <- is.na(x) | x == "."
                x[!na_idx] <- chartr("ACGTacgt", "TGCAtgca", x[!na_idx])
                x
              }

              blast.out$Ref[mismatch_known] <-
                comp_vec(blast.out$Ref[mismatch_known])

              if ("ALT" %in% names(blast.out)) {
                blast.out$ALT[mismatch_known] <-
                  comp_vec(blast.out$ALT[mismatch_known])
              }
            }
          }
        }
      }


      ## biallelic REF/ALT re-orientation to match ClusterA_NT / ClusterB_NT

      do_bial <- isTRUE(isbiallelic) ||
        tolower(as.character(isbiallelic)) == "yes"

      if (do_bial &&
          all(c("ClusterA_NT", "ClusterB_NT") %in% colnames(blast.out))) {

        A    <- toupper(as.character(blast.out$ClusterA_NT))
        B    <- toupper(as.character(blast.out$ClusterB_NT))
        refb <- toupper(as.character(blast.out$Ref))
        altb <- toupper(as.character(blast.out$ALT))

        comp_base <- function(x) {
          x <- as.character(x)
          na_idx <- is.na(x) | x == "."
          x[!na_idx] <- chartr("ACGTacgt", "TGCAtgca", x[!na_idx])
          x
        }

        # check if current REF/ALT already match {A,B} (either order)
        same_pair <- (!is.na(refb) & !is.na(altb)) &
          (
            (refb == A & altb == B) |
            (refb == B & altb == A)
          )

        # candidate complemented REF/ALT
        refc <- comp_base(refb)
        altc <- comp_base(altb)

        comp_pair <- (!is.na(refc) & !is.na(altc)) &
          (
            (refc == A & altc == B) |
            (refc == B & altc == A)
          )

        # flip when original doesn't match A/B but complemented pair does
        flip <- !same_pair & comp_pair

        if (any(flip)) {
          blast.out$Ref[flip] <- refc[flip]
          blast.out$ALT[flip] <- altc[flip]
        }
      }
    }

    # Generate a copy of blast pre filtering (after any orientation changes)
    blast.outraws <- blast.out

    # Drop rows with missing ALT
    blast.out <- blast.out[!is.na(blast.out$ALT), ]

    # Apply identity & coverage filters
    blast.out <- blast.out[
      as.numeric(blast.out$pident) >= pident.filter &
        blast.out$Coverage >= coverage.filter,
    ]

    if (nrow(blast.out) > 1) {
      rst_new <- data.frame(
        SNP_ID              = blast.out$qaccver,
        Unique_SNP_Locus_ID = blast.out$qaccver,
        Chr                 = blast.out$saccver,
        SNP_position        = format(as.numeric(blast.out$SNPpos), scientific = FALSE),
        REF_nt              = blast.out$Ref,
        ALT                 = blast.out$ALT,
        stringsAsFactors    = FALSE
      )

      if (!"Target_or_off-target?" %in% colnames(blast.out)) {
        rst_new$"Target_or_off-target?" <- "."
      } else {
        rst_new$"Target_or_off-target?" <- blast.out$"Target_or_off-target?"
      }

      blast.out     <- as.data.frame(lapply(blast.out,     function(x) as.character(format(x, scientific = FALSE))))
      blast.outraws <- as.data.frame(lapply(blast.outraws, function(x) as.character(format(x, scientific = FALSE))))

      # create minimal data frame of all the data
      if ("Hybridized" %in% colnames(blast.out)) {
        rst_new$Hybridised <- blast.out$Hybridized
      } else {
        rst_new$Hybridised <- NA
      }

      rst_new$Coverage <- blast.out$Coverage
      rst_new$pident   <- blast.out$pident

      if (!"ALT" %in% colnames(blast.out)) {
        rst_new$ALT <- "."
      }
    } else {
      rst_new <- blast.out
      colnames(rst_new) <- c(
        "SNP_ID",
        "Unique_SNP_Locus_ID",
        "Chr",
        "SNP_position",
        "REF_nt",
        "ALT",
        "Target_or_off-target",
        "Hybridised",
        "pident",
        "Coverage"
      )
    }

    nofilters.outpath <- file.path(
      output.path,
      paste0(probe.name, "_with_", genome.name, "_complete_unfiltered_blast_results.csv")
    )

    utils::write.csv(
      blast.outraws,
      file      = nofilters.outpath,
      quote     = FALSE,
      row.names = FALSE
    )

    all.outpath <- file.path(
      output.path,
      paste0(probe.name, "_with_", genome.name, "_all_mappings.csv")
    )

    utils::write.csv(
      blast.out,
      file      = all.outpath,
      quote     = FALSE,
      row.names = FALSE
    )

    out.path <- file.path(
      output.path,
      paste0(probe.name, "_with_", genome.name, "_mapping.tsv")
    )

    cat("\nNumber of rows in blast results:", nrow(blast.out), "\n")
    cat("\nMapping results written to ", out.path)

    # Rare edge case: if no markers present, rst_new can get NA column names (happens if few markers align and all which do are filtered leaving zero rows)
    if (is.data.frame(rst_new)) {
      cn <- colnames(rst_new)
      if (is.null(cn)) {
        cn <- rep("", ncol(rst_new))
        colnames(rst_new) <- cn
      }

      bad_name <- is.na(cn) | cn == "" | tolower(trimws(cn)) == "na"
      drop_cols <- bad_name

      if (any(drop_cols)) {
        rst_new <- rst_new[, !drop_cols, drop = FALSE]
      }
    }

    utils::write.table(
      rst_new,
      file      = out.path,
      quote     = FALSE,
      row.names = FALSE
    )

    return(blast.out)
  }
