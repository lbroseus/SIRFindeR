# IRFindeR2: Intron Retention Estimation and Detection
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
#########################################################

#_______________________________________________________#
#' Internal function for genome liftover
#'
#' @description This function performs genome lift over.
#' It is simply a wrapper for ['rtracklayer::liftover'].
#'
#' @importFrom rtracklayer import.chain liftOver
#' @importFrom GenomeInfoDb seqlevelsStyle
#'
#' @param chainPath A \code{Chain} object, usually imported with ['import.chain'], or something coercible to one.
#' @param oldIntervals A \code{GRanges} object with genomic intervals to lift over.
#'
#'
#' @return If out_file=NULL, a GRanges object containing intervals after lift over;
#'         otherwise the GRanges object will be saved in out_file.
#--------------------------------------------------------#
# TESTED

liftItOver <- function(chainPath, oldIntervals){

  chain <- import.chain(chainPath)

  seqlevelsStyle_old <- seqlevelsStyle(oldIntervals)
  GenomeInfoDb::seqlevelsStyle(oldIntervals) <- "UCSC"  # necessary

  newIntervals <- liftOver(oldIntervals, chain)
  newIntervals <- unlist( newIntervals )

  GenomeInfoDb::seqlevelsStyle(newIntervals) <- seqlevelsStyle_old

  return( newIntervals )

}

#_______________________________________________________#
#' Internal function for excluding low mappability regions from intronic intervals
#'
#' @description This function removes low mappability regions from intronic intervals.
#' Blacklisted regions are those having mappability score below a given threshold and spanning a minimal width (eg: read length).
#' In certain cases (depending on read length and your reference genome) pre-computed mappability scores can be downloaded
#' from ENCODE databases.
#'
#' @import stringr
#' @import GenomicFeatures
#'
#' @importFrom S4Vectors Hits queryHits subjectHits endoapply
#' @importFrom rlang .data
#' @importFrom data.table fread fwrite
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomicRanges GRanges makeGRangesListFromDataFrame setdiff findOverlaps 
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics strand start end 
#'
#' @param mappabilityFile Path to the (bed?) file containing mappability score to use.
#'                              If NULL, no mappability curation will be performed (Default: NULL).
#' @param chainPath (Default: NULL). A \code{Chain} object, usually imported with ['import.chain'], or something coercible to one.
#'                  Specify only if genome version liftover of mappability scores is needed.
#' @param minMapScore Masq regions with mappability score below this value (Default: 0.5)
#' @param readLength Consider only low mappability region longer than this value.
#'
#' @param saveDir Path to the directory where all sample results are saved.
#' @param verbose Whether information messages should be displayed (DEFAULT: TRUE).
#'
#' @return A \code{GRangesList} having same length as \code{indIntrons} listing sub-intervals
#'         with mappability score and width above the specified thresholds.
#'
#' @export
#--------------------------------------------------------#

curateIntrons <- function(mappabilityFile = NULL,
                          chainPath = NULL,
                          readLength,
                          minMapScore = 0.5,
                          saveDir,
                          verbose = TRUE){

  indIntronsFile <- paste(saveDir, "independent_introns.txt", sep = "/")

  indIntrons <- fread(file = indIntronsFile, data.table = F)
  indIntrons <- GRanges(indIntrons)
  indIntrons$blacklisted <- 0

  curatedIntronsFile <- paste(saveDir, "curated_introns.txt", sep = "/")

  curated_introns <- makeGRangesListFromDataFrame( data.frame(indIntrons) )

  if( !is.null(mappabilityFile) ){

    min.gapwidth <- 25

    if( verbose ) cat("Importing Mappability scores from", mappabilityFile, "\n")
    mappabilityScores <- import(mappabilityFile)

    ## Extract long enough blacklisted regions
    mappabilityScores <- mappabilityScores[ mappabilityScores$score<minMapScore ]
    mappabilityScores <- reduce(mappabilityScores, min.gapwidth = min.gapwidth)

    mappabilityScores <- mappabilityScores[ width(mappabilityScores)>=readLength ]
    mappabilityScores$score <- NULL

    ## lift over mappability scores
    if( !is.null(chainPath) ){

      if( verbose ) cat("Performing lift over...")

      mappabilityScores <- liftItOver(chainPath, mappabilityScores)

      if( verbose ) cat("OK.\n")

    }

    GenomeInfoDb::seqlevelsStyle(mappabilityScores) <- GenomeInfoDb::seqlevelsStyle(indIntrons) <- "Ensembl"

    within.hits <- suppressWarnings( findOverlaps(query = indIntrons, subject = mappabilityScores, type = "within") )

    indIntrons$blacklisted[ unique(queryHits(within.hits)) ] <- 1

    hits <- suppressWarnings( findOverlaps(mappabilityScores, indIntrons) )

    if( length(hits)> 0 ){

      hits <- setdiff(unique(subjectHits(hits)), unique(queryHits(within.hits)))

      print( length(hits) )

      #Avoid copying mappabilityScores
      curate <- function(granges){

        strand <- as.vector(granges@strand)[1]
        granges <- suppressWarnings( setdiff(granges, mappabilityScores, ignore.strand = TRUE) )
        BiocGenerics::strand(granges) <- strand

        return( granges )
      }

      curated_introns[ hits ] <- endoapply(X = curated_introns[ hits ],
                                           FUN = function(x) curate(x))

    }
  }else if( verbose ) cat("No curation of low mappability regions. \n")

    curated_introns <- unlist( curated_introns )
    curated_introns$intron_group <- names(curated_introns)

    fwrite(data.frame(curated_introns), file = curatedIntronsFile)


}
