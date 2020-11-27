# SIRFindeR: Intron Retention Estimation and Detection
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
#########################################################
#
# Calling IR events from RNA Nanopore data
#
#########################################################

#_______________________________________________________#
#' Compute the fraction of an intron contained in a long read
#'
#' @description This function computes the proportion
#' of an input genomic interval (eg: an intron or an exon) that is contained in a read.
#'
#' NB: it takes into account possible gaps and splicing in read alignments.
#'
#' @importFrom IRanges width ranges
#' @importFrom GenomicRanges findOverlaps pintersect
#' @importFrom GenomicAlignments extractAlignmentRangesOnReference cigar
#' @importFrom S4Vectors Hits queryHits subjectHits
#'
#' @param align.gr A \code{GRanges} object describing a read alignment (eg: extracted from a bam file)
#' @param intron.gr A \code{GRanges} object describing a genomic interval to compare to \code{align.gr}
#'
#'
#' @return The fraction of \code{intron.gr} contained in the sequence \code{align.gr}
#--------------------------------------------------------#
#

computeRegionOverlap <- function(align.gr, intron.gr){

  #cat(paste(bam_query_hit, ".", sep =""))
  granges <- extractAlignmentRangesOnReference(cigar( align.gr ),
                                               pos = start( align.gr ),
                                               drop.D.ranges=FALSE, f = NULL)

  z1 <- granges[[1]]
  z2 <- ranges( intron.gr )

  hits <- findOverlaps(query = z1, subject = z2, minoverlap = 3)

  x <- z1[ queryHits(hits) ]
  y <- z2[ subjectHits(hits) ]

  relative_overlap <- sum( width( pintersect(x, y)) ) / sum ( width(y) )

  return( relative_overlap )
}

#_______________________________________________________#
#' Call Intron Retention events from RNA long read alignments
#'
#' @description This function calls intron retention events from long read data
#' by crossing reference intron intervals and long read alignments.
#' It will:
#'
#' 1. detect reads which overlap at least one intron;
#'
#' 2. compute the proportion of each overlapped intron that is contained in a read;
#'
#' 3. assess whether the read extends beyond left and right intron borders
#'
#' NB: it takes into account possible gaps and splicing in read alignments.
#'
#' For more details and usage see vignette('IR-detection-from-Nanopore-data')
#'
#' @import magrittr
#'
#' @importFrom dplyr select
#' @importFrom S4Vectors Hits queryHits subjectHits
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom Rsamtools ScanBamParam
#' @importFrom GenomicAlignments readGappedReads
#' @importFrom GenomicRanges GRanges findOverlaps resize shift
#'
#' @param bam Either a path to the long read alignments (bam file) from which to call IR events,
#' or a GappedReads object.
#' @param intronGR A \code{GRanges} object listing intron genomic intervals to analyse (eg: generated using ['createIntronGR'])
#' @param flankWidth Minimal overlap width read should have over both exon-intron junctions.
#'                    Should be greater than 0.
#'                    (Default value: 20 bp - that is 10bp on each side of the junction position)
#' @param keepSecondaryAlignments Whether secondary alignments should be considered (Default: FALSE).
#'
#' @return A \code{data.frame} listing aligned long read overlapping at least one intron.
#'
#'
#' @examples
#'
#' ## First, create a GRanges with intron intervals from a reference transcriptome:
#' \dontrun{
#' gtf <- "Homo_sapiens.GRCh38.97.gtf"
#` intronGR <- createintronGR(gtf = gtf)
#' }
#' ## Then, call intron retention events using alignments in a bam file:
#' \dontrun{
#' bamFile <- "myLongReads.bam"
#' IR_events <- callIRevents(bamFile, intronGR)
#' }
#'
#' @export
#--------------------------------------------------------#
# TESTED

callIRevents <- function(bam, intronGR, flankWidth = 20, keepSecondaryAlignments = FALSE){
  
  fromFile <- FALSE
  
  ## Prior verifications
  if( flankWidth < 1 ) stop("flankWidth must be a positive integer.")
  if( !is(intronGR, 'GRanges') ) stop("intronGR must be a GRanges object.")
  
  if(is.character(bam)){
    if( !file.exists(bam) ) stop("the bam file does not exist.")
    fromFile <- TRUE
  }

  intron_minoverlap <- 5  # Minimum overlap between alignment and an intron (used for primary filtering)

  ## Importing splice alignments
  if( fromFile ){
    param <- ScanBamParam(scanBamFlag(isUnmappedQuery = F, isSecondaryAlignment = keepSecondaryAlignments))
    bam <- readGappedReads(file = bam, param = param, use.names = TRUE)
  }

  GenomeInfoDb::seqlevelsStyle(bam) <- GenomeInfoDb::seqlevelsStyle(intronGR) <- "UCSC"

  ## A quick findOverlap to extract potentially relevant alignments
  ## (allows to filter-out many uninformative short truncated reads)
  hits <- findOverlaps(query = bam, subject = intronGR, minoverlap = intron_minoverlap)

  ## Just reduce data in memory
  bam <- bam[ unique(queryHits(hits)) ]
  intronGR <- intronGR[ unique(subjectHits(hits)) ]

  hits <- findOverlaps(query = bam, subject = intronGR, minoverlap = intron_minoverlap)

  bam <- bam[ queryHits(hits) ]
  intronGR <- intronGR[ subjectHits(hits) ]

  start_time <- Sys.time()
  result <- tapply(1:length(bam), 1:length(bam), function(i) computeRegionOverlap(bam[i], intronGR[i]))
  end_time <- Sys.time()

  #ToDo: Possibly integrate TSS mark information
  overlap <- data.frame( data.frame(hits),
                         read_name = names(bam),
                         qwidth = qwidth(bam),
                         chr = seqnames(bam),
                         start = start(bam),
                         end = end(bam),
                         strand = strand(bam),
                         gene_id = mcols(intronGR)$gene_id,
                         intron_chr = seqnames(intronGR),
                         intron_start = start(intronGR),
                         intron_end = end(intronGR),
                         intron_strand = strand(intronGR),
                         intron_fracoverlap = result)

  ## Check 5' and 3' flanks

  # Check intron inner 5'SS
  intron_flank <- resize(intronGR, width = flankWidth, fix = "start")
  intron_flank <- shift(intron_flank, shift = -ceiling(flankWidth/2))

  flank_hits <- findOverlaps(query = bam, subject = intron_flank, minoverlap = flankWidth/2+5)
  flank_hits <- data.frame(flank_hits)
  flank_hits$intron_startOverlap <- TRUE

  overlap <- merge(overlap, flank_hits, by.x = c("queryHits", "subjectHits"), by.y=c("queryHits", "subjectHits"), all.x = T)
  overlap$intron_startOverlap[ which(is.na(overlap$intron_startOverlap)) ] <- FALSE

  # Check intron inner 3'SS
  intron_flank <- resize(intronGR, width = flankWidth, fix = "end")
  intron_flank <- shift(intron_flank, shift = ceiling(flankWidth/2))

  flank_hits <- findOverlaps(query = bam, subject = intron_flank, minoverlap = flankWidth/2+5)
  flank_hits <- data.frame(flank_hits)
  flank_hits$intron_endOverlap <- TRUE

  overlap <- merge(overlap, flank_hits, by = c("queryHits", "subjectHits"), all.x = T)
  overlap$intron_endOverlap[ which(is.na(overlap$intron_endOverlap)) ] <- FALSE

  overlap <- overlap %>% select(-queryHits, -subjectHits)

  return( overlap )

}



