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
#' @description This function is used to add more details to the blast output returned by the `processBlastResults()` function
#' @author David Chisanga
#' @param blast.results Dataframe or path to blast results returned by the `processBlastResults()` function
#' @param reference.bases Path to the reference bases obtained from samtools
#' @param probe.name name of the probe file being analysed
#' @param genome.name name of the reference genome being compared
#' @param output.path path were the results are saved
#' @param coverage.filter double value to filter hits with small coverage
#' @param is3prime Does the SNP exist as 3' as part of SNP chip marker mapping
#' @keywords genome, alternative allele, reference allele
#' @export
#'
#' @return Returns a dataframe with the updated blast results
#'
addSNPdetails <-
  function(blast.path,
           reference.bases,
           probe.name,
           genome.name,
           output.path = "Mapping_results",
           pident.filter = 90,
           coverage.filter = 50,
           is3prime)
  {
    #Create the output directory used to store the output files
    if (!dir.exists(output.path))
      dir.create(output.path, recursive = T)

    #Read samtools output file with the coordinates of the markers on the reference genome
    ref <- read.delim(reference.bases, header = F)

    #Check the file type that the blast results are saved in
    if ((!class(blast.path) %in% c("data.frame", "matrix", "data.table")))
      if (!file.exists(blast.path))
        stop("Blast results should be a dataframe or path to an RDS file")
    else
      blast.out <- readRDS(blast.path)
    #Add Ref coordinates to data frame
    blast.out$Ref <- ref$V2[match(blast.out$SNP.Refpos, ref$V1)]
    blast.out$Ref[is.na(blast.out$Ref)] <- "."
    if ("ClusterA_NT" %in% colnames(blast.out))
    {
    
# Add here an if for if istarget3primeend and then below for !is3prime
# Will need to add in the params variable into the function + add it into the module.nf and finally the main.n its a param so val param.is3prime will suffice for the chanel  
if(is3prime) {
      blast.out$ALT <- with(blast.out,
                            getAltBP(
                              NT_A = ClusterA_NT,
                              NT_B = ClusterB_NT,
                              Ref_NT = Ref
                            ))
}

if (!is3prime) {

    # normalise inputs
    ss   <- tolower(as.character(blast.out$sstrand))   # "plus" or "minus"
    ref  <- toupper(as.character(blast.out$Ref))
    a_nt <- toupper(as.character(blast.out$ClusterA_NT))
    b_nt <- toupper(as.character(blast.out$ClusterB_NT))

    # convert to plus-strand bases (complement if on minus)
    plusA <- ifelse(ss == "minus", chartr("ACGT", "TGCA", a_nt), a_nt)
    plusB <- ifelse(ss == "minus", chartr("ACGT", "TGCA", b_nt), b_nt)

    # set ALT: if one equals Ref, ALT is the other; else "[A/B]" using plus-strand bases
    blast.out$ALT <- ifelse(
        is.na(plusA) | is.na(plusB) | is.na(ref),
        NA_character_,
        ifelse(
            plusA == ref, plusB,
            ifelse(
                plusB == ref, plusA,
                paste0("[", plusA, "/", plusB, "]")
            )
        )
    )
}
    }
    # Generate a copy of blast pre filtering
    blast.outraws <- blast.out
    blast.out<-blast.out[!is.na(blast.out$ALT),]
    blast.out <-
      blast.out[as.numeric(blast.out$pident) >= pident.filter &
                  blast.out$Coverage >= coverage.filter, ]
    
    #blast.out <- blast.out[blast.out$Hybridized!="No",]

    if (nrow(blast.out) > 1)
    {
      rst_new <- data.frame(
        SNP_ID = blast.out$qaccver,
        Unique_SNP_Locus_ID = blast.out$qaccver,
        Chr = blast.out$saccver,
        SNP_position = format(as.numeric(blast.out$SNPpos), scientific = F),
        REF_nt = blast.out$Ref,
        ALT = blast.out$ALT,
        stringsAsFactors = F
      )

      if (!"Target_or_off-target?" %in% colnames(blast.out))
         rst_new$"Target_or_off-target?" <- "."
      else
         rst_new$"Target_or_off-target?" <- blast.out$"Target_or_off-target?"
      blast.out <- as.data.frame(apply(blast.out, 2, function(x)
        as.character(format(x, scientific = F))))
      blast.outraws <- as.data.frame(apply(blast.outraws, 2, function(x)
        as.character(format(x, scientific = F))))
      #Create minimal data frame of all the data
      if("Hybridized" %in% colnames(blast.out))
      {
          rst_new$Hybridised<-blast.out$Hybridized
      }
      if(!"Hybridized" %in% colnames(blast.out))
      {
          rst_new$Hybridised<-NA
      }

      rst_new$Coverage = blast.out$Coverage
      rst_new$pident = blast.out$pident
      if (!"ALT" %in% colnames(blast.out))
        rst_new$ALT <- "."
    }
    else{
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

    nofilters.outpath <- file.path(output.path,
                                   paste0(probe.name, "_with_", genome.name, "_complete_unfiltered_blast_results.csv"))
    
    write.csv(blast.outraws,
              file = nofilters.outpath,
              quote = F,
              row.names = F)
    
    
    all.outpath <- file.path(output.path,
                             paste0(probe.name, "_with_", genome.name, "_all_mappings.csv"))
    write.csv(blast.out,
              file = all.outpath,
              quote = F,
              row.names = F)
    out.path <-
      file.path(output.path,
                paste0(probe.name, "_with_", genome.name, "_mapping.tsv"))

    cat("\nNumber of rows in blast results:", nrow(blast.out), "\n")
    cat("\nMapping results written to ", out.path)


# There is a rare edge case where if no markers are present in a subset the mappings file gets populated by NA
#Here is a quick remove NA cols to sort this out 

if (is.data.frame(rst_new)) {
  # Normalize/inspect column names
  cn <- colnames(rst_new)
  if (is.null(cn)) {
    cn <- rep("", ncol(rst_new))
    colnames(rst_new) <- cn
  }

  # columns with bad names: NA, "", or literal "NA" (case-insensitive)
  bad_name <- is.na(cn) | cn == "" | tolower(trimws(cn)) == "na"


  drop_cols <- bad_name

  if (any(drop_cols)) {
    rst_new <- rst_new[, !drop_cols, drop = FALSE]
  }

}

    write.table(rst_new,
                file = out.path,
                quote = F,
                row.names = F)
    return(blast.out)
  }
