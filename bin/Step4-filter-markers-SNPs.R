#Get input arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied")
}else {
  cat("\nParameters passed:\n")
  cat(args, sep = "\n")
}
args.list <- as.list(gsub(".*=", "", args))
names(args.list) <- gsub("=.*", "", args)
if (!"savedblast" %in% names(args.list))
{
  stop("Results from blast not found. Make sure that the blast step ran successfully before running this step!")
}
blast.results <- args.list$savedblast
reference.bases <- args.list$refSNPs
probe.name <- args.list$probename
genome.name <- args.list$genomename
output.path <- args.list$resultsdir #"Mapping_results"

#' @title Get the alternative base for a probe ID
#' @description This function used is to get the alternative base pair from blast results
#' @author David Chisanga
#' @param NT_A Design cluster A nucleotide
#' @param NT_B Design cluster B nucleotide
#' @param Ref_NT Nucleotide at SNP position from reference genome
#' @keywords genome, alternative allele, reference allele
#' @export
#' @examples
#'
#' @return Returns a character string indicating the alternative base pair
#'
getAltBP <- function(NT_A, NT_B, Ref_NT)
{
  #convert nucleotides to uppercase
  NT_A <- toupper(trimws(NT_A))
  NT_B <- toupper(trimws(NT_B))
  Ref_NT <- toupper(trimws(Ref_NT))
  sapply(1:length(NT_A), function(i) {

    if (NT_A[i] == Ref_NT[i])
    {
      out <- NT_B[i]
    }
    else if (Ref_NT[i] == NT_B[i])
    {
      out <- NT_A[i]
    }
    else if (NT_A[i] == "T" &
      Ref_NT[i] == "A" &
      NT_B[i] == "C")
    {
      out <- "G"
    }
    else if (NT_A[i] == "T" &
      Ref_NT[i] == "A" &
      NT_B[i] == "G")
    {
      out <- "C"
    }
    else if (NT_A[i] == "A" &
      Ref_NT[i] == "T" &
      NT_B[i] == "G")
    {
      out <- "C"
    }
    else if (NT_A[i] == "A" &
      Ref_NT[i] == "T" &
      NT_B[i] == "C")
    {
      out <- "G"
    }
    else if (NT_A[i] == "C" &
      Ref_NT[i] == "G" &
      NT_B[i] == "A")
    {
      out <- "T"
    }
    else if (NT_A[i] == "C" &
      Ref_NT[i] == "G" &
      NT_B[i] == "T")
    {
      out <- "A"
    }
    else if (NT_A[i] == "G" &
      Ref_NT[i] == "C" &
      NT_B[i] == "A")
    {
      out <- "T"
    }
    else if (NT_A[i] == "G" &
      Ref_NT[i] == "C" &
      NT_B[i] == "T")
    {
      out <- "A"
    }
    else if (Ref_NT[i] == NT_A[i])
    {
      out <- NT_B[i]
    }
    else if (NT_B[i] == "T" &
      Ref_NT[i] == "A" &
      NT_A[i] == "C")
    {
      out <- "G"
    }
    else if (NT_B[i] == "T" &
      Ref_NT[i] == "A" &
      NT_A[i] == "G")
    {
      out <- "C"
    }
    else if (NT_B[i] == "A" &
      Ref_NT[i] == "T" &
      NT_A[i] == "G")
    {
      out <- "C"
    }
    else if (NT_B[i] == "A" &
      Ref_NT[i] == "T" &
      NT_A[i] == "C")
    {
      out <- "G"
    }
    else if (NT_B[i] == "C" &
      Ref_NT[i] == "G" &
      NT_A[i] == "A")
    {
      out <- "T"
    }
    else if (NT_B[i] == "C" &
      Ref_NT[i] == "G" &
      NT_A[i] == "T")
    {
      out <- "A"
    }
    else if (NT_B[i] == "G" &
      Ref_NT[i] == "C" &
      NT_A[i] == "A")
    {
      out <- "T"
    }
    else if (NT_B[i] == "G" &
      Ref_NT[i] == "C" &
      NT_A[i] == "T")
    {
      out <- "A"
    }
    else
      out <- "."
    return(out)
  }, simplify = T, USE.NAMES = F)
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
#' @examples
#'
#' @return Returns a character vector showing the target results
#'
#'

checkTarget <-
  function(target.chrom,
           extendable,
           SNP.position,
           target.bp,
           ref.chrom)
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
#' @description This function is used to filter probes that had multiple hits `addSNPdetails()` function
#' @author David Chisanga
#' @param blast.results Dataframe or path to blast results returned by the `addSNPdetails()` function
#' @param filters
#' @keywords genome, alternative allele, reference allele
#' @examples
#'
#' @return Returns a dataframe with the updated blast results
#'
filterMarkers <- function(blast.results, filters = NULL)
{
  #blast.results$Keep_or_delete <- "keep"

  #Get a subset of markers with multiple hits and rank them
  dups <-
    unique(blast.results$qaccver[duplicated(blast.results$qaccver)])

  if ("Identical_nts_from_3prime_end" %in% colnames(blast.results)) {
    dups <- do.call(rbind, lapply(dups, function(x) {
      xx <- subset(blast.results, qaccver == x)
      xx <-
        xx[with(
          xx,
          order(
            pident,
            Identical_nts_from_3prime_end,
            qlen,
            length,
            mismatch,
            gapopen,
            evalue,
            bitscore,
            method = "radix",
            decreasing = c(T, T, T, F, F, F, T)
          )
        ),]
      xx <- rbind(xx[xx$`Target_or_off-target?` == "target",],
                  xx[xx$`Target_or_off-target?` != "target",])
      # xx[1, "Keep_or_delete"] <- "keep"
      # xx[-1, "Keep_or_delete"] <- "delete"
      return(xx)
    })) }else
  {
    dups <- do.call(rbind, lapply(dups, function(x) {
      xx <- subset(blast.results, qaccver == x)
      xx <-
        xx[with(
          xx,
          order(
            pident,
            qlen,
            length,
            mismatch,
            gapopen,
            evalue,
            bitscore,
            method = "radix",
            decreasing = c(T, T, F, F, F, T)
          )
        ),]
      return(xx)
    }))
  }

  blast.results <- blast.results[!blast.results$qaccver %in% dups$qaccver,]
  blast.results <- rbind(blast.results, dups)
  if ("Identical_nts_from_3prime_end" %in% colnames(blast.results))
    blast.results <- blast.results[with(
      blast.results,
      order(
        pident,
        Identical_nts_from_3prime_end,
        qlen,
        length,
        mismatch,
        gapopen,
        evalue,
        bitscore,
        method = "radix",
        decreasing = c(T, T, T, F, F, F, T)
      )
    ),]
  else
    blast.results <- blast.results[with(
      blast.results,
      order(
        pident,
        qlen,
        length,
        mismatch,
        gapopen,
        evalue,
        bitscore,
        method = "radix",
        decreasing = c(T, T, T, F, F, F, T)
      )
    ),]
  
  return(blast.results)
}


#' @title Add SNP details to the blast output
#' @description This function is used to add more details to the blast output returned by the `processBlastResults()` function
#' @author David Chisanga
#' @param blast.results Dataframe or path to blast results returned by the `processBlastResults()` function
#' @param reference.bases Path to the reference bases obtained from samtools
#' @param probe.name name of the probe file being analysed
#' @param genome.name name of the reference genome being compared
#' @param output.path path were the results are saved
#' @keywords genome, alternative allele, reference allele
#' @export
#' @examples
#'
#' @return Returns a dataframe with the updated blast results
#'
addSNPdetails <-
  function(blast.results,
           reference.bases,
           probe.name,
           genome.name,
           output.path = "Mapping_results")
  {
    #Check the file type that the blast results are saved in
    if ((!class(blast.results) %in% c("data.frame", "matrix")))
      if (!file.exists(blast.results))
        stop("Blast results should be a dataframe or path to an RDS file")
      else
        blast.results <- readRDS(blast.results)
    #Create the output directory used to store the 2 output files
    if (!dir.exists(output.path))
      dir.create(output.path, recursive = T)
    #Read samtools output file with the coordinates of the markers on the reference genome
    ref <- read.delim(reference.bases, header = F)

    #Add Ref coordinates to dataframe
    blast.results$Ref <-
      ref$V2[match(blast.results$SNP_position, ref$V1)]
    blast.results$Ref[is.na(blast.results$Ref)] <- "."
    if ("ClusterA_NT" %in% colnames(blast.results))
    {
      blast.results$ALT <- with(blast.results,
                                getAltBP(
                                  NT_A = ClusterA_NT,
                                  NT_B = ClusterB_NT,
                                  Ref_NT = Ref
                                ))
    }


    if ("Target.Chr" %in% colnames(blast.results)) {

      blast.results[, "Target_or_off-target?"]<- with(
        blast.results,
        checkTarget(
          target.chrom = Target.Chr,
          extendable = ExtendableSite,
          SNP.position = SNPpos,
          target.bp = Target.bp,
          ref.chrom = saccver
        )
      )
      blast.results$TargetHit <-
        as.character(apply(blast.results, 1, function(x) {
          #if(AL2=".",".",IF(AND(AL2=RIGHT(E2,2),AM2=AH2),1,0))
          out <- "."
          if (x[["Target.Chr"]] == ".")
            out <- "."
          else if (gsub("^chr","",x[["Target.Chr"]],ignore.case = T) == gsub("^chr","",x[["saccver"]],ignore.case = T) &
            x[["SNPpos"]] == x[["Target.bp"]])
            out <- 1
          else
            out <- 0
        }))
      print(blast.results[blast.results$qaccver=="uat135",])
      blast.results$Hybridised<-"Yes"
      blast.results$Hybridised[blast.results$ExtendableSite=="."]<-"No"
      # blast.results <-
      #   subset(blast.results, ExtendableSite != "." & SNP_position != ".")
      blast.results$Total_no_extendable_hits_for_probe <-
        as.numeric(with(blast.results, table(paste0(
          qaccver, ExtendableSite
        ))[paste0(qaccver, ExtendableSite)]))
      blast.results$Total_no_extendable_hits_for_probe[blast.results$ExtendableSite !=
                                                         "Yes"] <- 0
    }

    blast.results <- filterMarkers(blast.results)
    blast.results$qaccver.dup <-blast.results$qaccver
    blast.results <-do.call(rbind, lapply(unique(blast.results$qaccver.dup),function(x){
    x.data<-subset(blast.results,qaccver.dup==x)
    if(nrow(x.data)>1)
    {
       print(paste0("Duplicates for ",x)) 
       x.data$qaccver.dup<-paste0(x.data$qaccver.dup,".",seq(1, nrow(x.data), 1))     
    }
    return(x.data)
    }))
    rst_new <- blast.results
    if(!"ALT" %in% colnames(rst_new))
      rst_new$ALT<- "."
    if(!"Target_or_off-target?" %in% colnames(rst_new))
      rst_new$"Target_or_off-target?"<- "."
    rst_new <- data.frame(
      SNP_ID = rst_new$qaccver,
      Alternate_SNP_ID = rst_new$qaccver.dup,
      Unique_SNP_Locus_ID = rst_new$qaccver,
      Chr = rst_new$saccver,
      SNP_position =format(as.numeric(rst_new$SNPpos),scientific = F),
      REF_nt = rst_new$Ref,
      ALT = rst_new$ALT,
      "Target_or_off-target?" = rst_new$"Target_or_off-target?",
      Hybridize=rst_new$Hybridize
    )
    all.outpath <- file.path(output.path,
                             paste0(probe.name, "_with_", genome.name, "_all_mappings.csv"))
    blast.results <- apply(blast.results,2,function (x) format(x,scientific = F))
    blast.results <- apply(blast.results, 2, as.character)
    write.csv(blast.results, file = all.outpath, quote = F, row.names = F)
    out.path <-file.path(output.path, paste0(probe.name, "_with_", genome.name, "_mapping.tsv"))

    cat("\nNumber of rows in blast results:", nrow(blast.results), "\n")
    cat("\nMapping results written to ", out.path)

    write.table(rst_new,
                file = out.path,
                quote = F,
                row.names = F)
    return(blast.results)
  }

addSNPdetails(blast.results,
              reference.bases,
              probe.name,
              genome.name,
              output.path
)