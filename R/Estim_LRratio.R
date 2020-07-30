# SIRFindeR: Intron Retention Estimation and Detection
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
#########################################################


#_______________________________________________________#
#' Compute observed IR rates from a Long Read experiment
#'
#' @description This function computes the observed rates of intron
#' retention from a set of long read (genomic) alignments.
#'
#' @import magrittr
#' 
#' @importFrom dplyr filter group_by summarise n
#' @importFrom data.table fread fwrite
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomicAlignments readGappedReads
#'
#' @param gtf Path to a reference transcriptome is a gtf file.
#' @param bamFile Path to a bam file contraining long read alignments.
#' @param saveDir Path to the directory where results sould be saved.
#' 
#' @param flankWidth Width of the area on each side of the intron that a read should overlap
#' to be taken into account (Default: 20bp).
#' @param keepSecondaryAlignments Whether secondary alignments should be considered (Default: FALSE).
#'
#' @return Saves a \code{data.frame} containing observed intron retention rates
#' for each annotated intron in a file named "IRrates.txt" in \code{saveDir}.
#' 
#' @export
#'
#--------------------------------------------------------#

computeIRrates <- function(gtf, bamFile, saveDir,
                           flankWidth = 20, 
                           keepSecondaryAlignments = FALSE){
  
  #---------------------------------------------------------#
  IRevents.file <- paste(saveDir, "IRevents.txt", sep = "/")
  IRrates.file <- paste(saveDir, "IRrates.txt", sep = "/")
  
  param <- ScanBamParam(scanBamFlag(isUnmappedQuery = F, 
                                    isSecondaryAlignment = keepSecondaryAlignments))
  #---------------------------------------------------------#
  
  #Extract annotated introns
  cat("Creating reference intron annotation...")
  intronGR <- createIntronGR(gtf = gtf)
  GenomeInfoDb::seqlevelsStyle(intronGR) <- "UCSC"
  cat("OK.\n")
  
  #Count reads
  cat("Importing read alignments...")
  bam <- readGappedReads(file = bamFile, param = param, use.names = TRUE)
  GenomeInfoDb::seqlevelsStyle(bam) <-"UCSC"
  cat("OK.\n")
  
  cat("Counting intron spanning reads...")
  hits <- findOverlaps(query = intronGR, subject = GRanges(bam), type = "within")
  hits <- data.frame(hits) %>% group_by(queryHits) %>% summarise(readCount = n()) %>% data.frame()
  intronGR$readCount <- 0
  intronGR$readCount[ hits$queryHits ] <- hits$readCount
  cat("OK.\n")
  
  cat("Calling IR events...")
  IRevents <- callIRevents(bam, intronGR, flankWidth = flankWidth, keepSecondaryAlignments = keepSecondaryAlignments)
  cat("Saving...")
  fwrite(IRevents, file = IRevents.file, sep = "\t")
  cat("OK.\n")
  
  cat("Computing IR rates...")  
  IRevents <- IRevents %>% 
              filter(intron_fracoverlap==1) %>% 
              group_by(intron_chr, intron_start, intron_end, intron_strand) %>%
              summarise(intronicAbundance = n()) %>% data.frame()
  
  hits <- findOverlaps(query = GRanges(IRevents), subject = intronGR, type = "equal")
  
  intronGR$intronicAbundance <- 0
  intronGR$intronicAbundance[ subjectHits(hits) ] <- IRevents$intronicAbundance[ queryHits(hits) ]
  
  IRrates <- data.frame(intronGR)   
  IRrates <- IRrates %>% mutate(ratio = ifelse(readCount==0, 0, intronicAbundance/readCount)) %>% data.frame()
  
  # Save results
  cat("Saving...")
  fwrite(IRrates, file = IRrates.file, sep = "\t")
  
  cat("OK.\n")
}