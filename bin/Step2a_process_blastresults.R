args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied")
}else {
  cat("\nParameters passed:\n")
  cat(args, sep = "\n")
}
args.list <- as.list(gsub(".*=", "", args))
names(args.list) <- gsub("=.*", "", args)
args.list$blastformat <- strsplit(args.list$blastformat, " ")[[1]][-1]
if ("blastfile" %in% names(args.list))
{
  blast.file <- args.list$blastfile
  if (is.na(blast.file))
    stop("Results from blast not found. Make sure that the blast step ran successfully before running this step!")
}else {
  stop("Results from blast not found. Make sure that the blast step ran successfully before running this step!")
}
min.matchbp <- args.list$minmatches
extendable.site.bps <- as.numeric(args.list$extBPs)
path.2save.coords <- args.list$outsam #"SNPcoordinates"
istarget3primeend <- args.list$istarget3primeend
#Check if targets file has been provided
targets <- args.list$target
target.provided<-F
if (targets!= "none")
{
  #Table with a list of targets ordered as ID;Target.chromosome;TargetPosition;TargetSNP
  if (file.exists(targets))
  {
    
    if (grepl("RDS$", targets, ignore.case = T))
      targets <- readRDS(targets)
    else
    {
      targets <- read.table(targets, stringsAsFactors = F, header = F)
      if ("ID" %in% targets[1,])
      {
        colnames(targets)<-targets[1,]
        targets<-targets[-1,]
      }
      else
      {
        colnames(targets)<-c("ID","Target.Chr","Target.bp","Targetted.SNP")
      }
    }
    if(nrow(targets)>0)
    {
       target.provided=T 
    }
    else{
       stop("Targets file provided has null content, check file and try again") 
    }
  }else {
    stop(paste0("Targets file was not found at ", targets))
  }
}

path2BlastOut <- args.list$blastoutpath
blastOutFormat <- args.list$blastformat

#' @title Process alignment results from blast
#' @description This function is used to process and format the blast output
#' @author David Chisanga
#' @param blast.file path to the blast output file
#' @param min.matchbp minimum length of HSP to be considered fully Hybridize
#' @param target table with target annotation
#' @param outfmt vector of format specifiers from the supported format specifiers for 6,7 and 10 in blastn's 'outfmt' parameter
#' @param path.2save.coords path where SNP coordinates to get base in reference genome using samtools are to be saved
#' @param extendable.site.bps number of matching base pairs from the 3 prime end for a probe to be considered as extendable
#' @param save.as.RDS Boolean value showing whether to save results as an RDS object which can later be used in further downstream analysis
#' @param path2BlastOut If `save.as.RDS` is True, provide path where object will be stored. File is named after the blaSt file
#' @keywords blast, alignment
#' @export
#' @examples
#'
#' @return Returns a vector list of Discontinued and Current GeneIds from the list entered
#'
#'
processBlastResults <- function(blast.file,
                                min.matchbp = 40,
                                path.2save.coords = "SNPcoordinates",
                                extendable.site.bps = 3,
                                save.as.RDS = T,
                                path2BlastOut = "ProcessedRDS",
                                targets = NULL, #readRDS("v1.1_TargetSNP.RDs"),
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
                                ))
{

  #Read blast file to a dataframe
  blast_out <- read.table(
    file = blast.file,
    row.names = NULL,
    strip.white = T,
    fill = T,
    sep = "\t",
    nrows = 5000
  )

  #Rename the columns to use the blast output format used
  colnames(blast_out) <- outfmt
  blast_out$Hybridize <-
    ifelse(
      as.numeric(blast_out$qend) == as.numeric(blast_out$qlen) &
        as.numeric(blast_out$length) >= min.matchbp,
      "Yes",
      "No"
    )

  #Get the number of matching bases from the 3 prime end of the probe if target is provided
  if (target.provided)
  {

    blast_out$Identical_nts_from_3prime_end <-
      gsub(".*_",
           "",
           gsub("([A-Z][A-Z])|([A-Z]-)|(-[A-Z])", "_", blast_out$btop))


    #Check if there are no mismatches in at least 3bps from the 3 prime end of the probe
    blast_out$No_missmatch3bp <-unlist(
      sapply(as.numeric(as.numeric(blast_out$Identical_nts_from_3prime_end)),function(x)
             if(x>=extendable.site.bps)
               return("Yes")
             else
               return(".")
      ),use.names=F)

    #If the number of matching bps from the 3prime end is  more
    blast_out$ExtendableSite <-unlist(sapply(blast_out$Identical_nts_from_3prime_end,function(x){
      if(as.numeric(x)>=extendable.site.bps)
        return("Yes")
      else
        return(".")
    }),use.names = F)
  }
  #Get position of SNPs to be used in samtools to get basepairs
  if (istarget3primeend=="yes")
  {
    pos <-
      ifelse(
        blast_out$sstrand == "plus",
        as.numeric(blast_out$send) + 1,
        as.numeric(blast_out$send) - 1
      )
  }else{
    pos <- "."
  }

  #Check the “Blast trace-back operations” (BTOP) string which describes the alignment produced by BLAST.
  blast_out <- as.data.frame(t(apply(blast_out, 1, function(x) {
    #Split the btop string to check for SNPs or other variant types
    z <- strsplit(gsub("([A-Z][A-Z])|([A-Z]-)|(-[A-Z])", "_", x[["btop"]]), "_")[[1]]
    #Specify
    zz.type <- strsplit(gsub("[0-9]", "", gsub("([A-Z]-)|(-[A-Z])", "G",
                                               gsub("([A-Z][A-Z])", "m", x[["btop"]]))), split = "")[[1]]
    z[z == ""] <- 0
    z <- as.numeric(z) + 1
    zz.type <- paste0(zz.type, collapse = ";")
    z <- cumsum(z[-length(z)])
    if (length(z) == 0) {
      z <- "-"
      SNPorGAPpos_Ref <- "-"
    }else {
      SNPorGAPpos_Ref <- sapply(z, function(xx)
      {
        xx <- ifelse(x[["sstrand"]] == "plus",
                     as.numeric(x[["send"]]) + xx,
                     as.numeric(x[["send"]]) - xx)
        paste0(x[["saccver"]], ":", xx, "-", xx)
      }, simplify = T, USE.NAMES = F)
      SNPorGAPpos_Ref <- paste0(SNPorGAPpos_Ref, collapse = ";")
    }
    if (zz.type == "")
      zz.type <- "-"

    x[["SNPorGap"]] <- zz.type
    x[["SNPorGap_pos"]] <- paste0(z, collapse = ";")
    x[["SNPorGAPpos_Ref"]] <- SNPorGAPpos_Ref
    return(x)
  })))

  blast_out$SNP_position <-paste0(blast_out$saccver,
                  ":",
                  pos,
                  "-",
                  pos)
  blast_out$SNPpos <- pos
    
  if (target.provided)
  {
    #Include target information if available
    targets <- targets[!duplicated(targets),]
    chr.saccver<-blast_out$saccver[1]
    chr.saccver<-grepl("chr",chr.saccver,ignore.case = T)
    chr.qaccver <- with(targets,any(grepl("chr",Target.Chr,ignore.case = T)))
    # if(chr.saccver & !chr.qaccver)
    # {
    #   targets$Target.Chr<-paste0("chr",targets$Target.Chr)
    # }
    #Match targets using ID, target chrom, and position
    #targets$Target.Chr<-with(targets,ifelse(grepl("^chr",))
    target <-
      with(targets, split(Targetted.SNP, ID))
    blast_out$Target <- as.character(target[blast_out$qaccver])
    blast_out$ClusterA_NT <-
      gsub("\\[|\\]|/[ACTG]", "", blast_out$Target)
    blast_out$ClusterB_NT <-
      gsub("[ACTG]/|\\[|\\]", "", blast_out$Target)
    
    #Include target chromosome information from the probe design matrix
    target.chr <- with(targets, split(Target.Chr, ID))
    blast_out$Target.Chr <-
      as.character(target.chr[blast_out$qaccver])
    
    #Include target chromosome information from the probe design matrix
    target.bp <- with(targets, split(Target.bp, ID))
    blast_out$Target.bp <- as.character(target.bp[blast_out$qaccver])
  }

  #Reorder results
  blast_out[order(
    blast_out$Hybridize,
    blast_out$pident,
    blast_out$qaccver,
    decreasing = c(T, T, F),
    method = "radix"
  ),]

  #Save SNP coordinates to a tab delimited file
  #sam.coord <-file.path(path.2save.coords,"samtools_coordinates.tab")
  cat("\nSNP coordinates written to", path.2save.coords, "\n")
  write.table(
    subset(blast_out, SNP_position != ".")$SNP_position,
    file = path.2save.coords,
    quote = F,
    col.names = F,
    row.names = F
  )

  if (save.as.RDS)
  {
    # if (!dir.exists(path2BlastOut))
    #   dir.create(path2BlastOut, recursive = T)
    # path2BlastOut <- file.path(path2BlastOut,"blast_results.RDS")
    # cat("\nProcessed results are written to", path2BlastOut)
    saveRDS(blast_out,
            file = path2BlastOut)
  }
  else
    return(blast_out)
}


#### Process blast results #####
processBlastResults(
  blast.file = blast.file,
  min.matchbp = min.matchbp,
  path.2save.coords = path.2save.coords,
  extendable.site.bps = extendable.site.bps,
  save.as.RDS = T,
  path2BlastOut = path2BlastOut,
  targets = targets,
  outfmt = blastOutFormat
)

