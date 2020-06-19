# SIRFindeR: Stabilized estimation of Intron Retention levels
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
#########################################################

#_______________________________________________________#
#' Internal function for merging introns detected by STAR
#'
#' @description This junctions merges junction files
#'
#' @importFrom dplyr group_by summarise mutate
#' @importFrom data.table fread fwrite
#'
#' @param junctionFiles A character vector with STAR output files to be merged (.tab)
#' @param min_map_cross_junc Minimal number of reads supporting a junction to consider a novel intron
#' @param min_novel_intron_width Minimal width for an intronic interval to be kept
#' @param verbose (Default TRUE) Whether to print summary messages
#'
#' @seealso
#' * ['defineIntronicIntervals()'] to define the intronic intervals (...)
#'
#' @return A data.frame containing intronic intervals detected by STAR
#'         and junction counts summarised over all input junction files
#'
#--------------------------------------------------------#
# TESTED -
# POSSIBLE IMPROVEMENTS: include results from other junction calls

mergeJunctionFiles <- function(junctionFiles,
                               min_map_cross_junc = 5,
                               min_novel_intron_width = 30,
                               verbose = TRUE){

  mergedJunctions <- data.frame()

  mergedJunctions <- tapply(X = junctionFiles,
                            INDEX = 1:length(junctionFiles),
                            FUN = function(x){
                              mergedJunctions <- rbind.data.frame(mergedJunctions, data.frame(fread(x)))
                            })

  mergedJunctions <- do.call(mergedJunctions, what = rbind.data.frame)

  colnames(mergedJunctions) <- c("chr",
                                 "start",
                                 "end",
                                 "strand",
                                 "intron_motif",
                                 "annotated",
                                 "uniq_map_cross_junc",
                                 "mult_map_cross_jun",
                                 "max_splice_align_overhang")

  mergedJunctions <- mergedJunctions %>%
    filter(uniq_map_cross_junc>=min_map_cross_junc & end-start+1>=min_novel_intron_width) %>%
    group_by(chr, start, end, strand, intron_motif, annotated) %>%
    summarise(uniq_map_cross_junc = sum(uniq_map_cross_junc), nbSamples = n()) %>%
    data.frame()

  mergedJunctions <- mergedJunctions %>% mutate(strand = ifelse(strand == 1, "+", "-"))

  if( verbose ) cat("Number of potential introns detected in the experiment:", nrow(mergedJunctions), "\n")

  return( mergedJunctions )

}

#_______________________________________________________#
#' Internal function to define sample-specific intronic intervals for read counting
#'
#' @description This functions determines sample-specific intron borders used
#' for counting intronic reads.
#'
#' @import stringr
#' @import GenomicFeatures
#'
#' @importFrom GenomicRanges findOverlaps reduce
#' @importFrom S4Vectors Hits queryHits subjectHits mcols
#' @importFrom data.table fread fwrite
#'
#' @param gtf Path towards a reference transcriptome (gtf or gff format)
#' @param saveDir Path to where output files will be saved
#'
#' @param junctionFile Junction file output by STAR (.tab files)
#' @param keepAnnotatedRI Whether introns annotated as "retained" in the gtf should be considered.
#'                          (ie: transcripts annotated as "retained_intron" will be removed prior to defining intronic intervals)
#'                          (Default: TRUE)
#'
#' @param min_map_cross_junc Minimal number of reads supporting a junction to consider a novel intron
#' @param min_novel_intron_width Minimal width for an intronic interval to be kept
#'
#' @param verbose (Default: TRUE)
#'
#' @return Outputs curated and independent intronic intervals
#'
#--------------------------------------------------------#
# TESTED on ENSEMBL

defineIntronicIntervals <- function(gtf,
                                    saveDir,
                                    junctionFile,
                                    keepAnnotatedRI = FALSE,
                                    min_map_cross_junc = 5,
                                    min_novel_intron_width = 30,
                                    verbose = TRUE){

  if( ! dir.exists(saveDir) )  dir.create(saveDir)

  # Curated and Independent intronic intervals for read counting
  indIntronFile <- paste(saveDir, "independent_introns.txt", sep = "/")
  # Augmented annotation of all introns
  allIntronFile <-   paste(saveDir, "all_introns.txt", sep = "/")

  #-------------------------------------------------------------------#
  ## IMPORT ANNOTATED INTRONS

  if( verbose ) cat("Importing transcriptome annotation...")

  ref <- rtracklayer::import(gtf)
  colnames(S4Vectors::mcols(ref)) <- str_replace(string = colnames(S4Vectors::mcols(ref)),
                                                 pattern = "biotype", replacement = "type")
  ref <- ref[ !ref$transcript_type %in% c("retained_intron") ]

  txdb <- suppressWarnings( makeTxDbFromGRanges(ref) )

  introns <- intronsByTranscript(txdb, use.names = TRUE)
  introns <- unlist(introns)
  introns <- unique(introns)
  introns$Annotated <- 1

  if( verbose ) cat("OK. \n")

  #-------------------------------------------------------------------#
  ## INTEGRATE DETECTED INTRONS

  if( !is.null(junctionFiles) ){

    if( verbose ) cat("Importing introns detected in RNA-seq...")

    star_introns <- mergeJunctionFiles(junctionFiles = junctionFile, verbose = F)
    star_introns <- GRanges(star_introns)

    matches <- findOverlaps(introns, star_introns, type = "equal")
    introns$Detected <- 0
    introns$Detected[ queryHits(matches) ] <- 1

    introns$uniq_map_cross_junc <- 0
    introns$uniq_map_cross_junc[ queryHits(matches) ] <- star_introns$uniq_map_cross_junc[ subjectHits(matches) ]

    introns$intron_motif <- NA
    introns$intron_motif[ queryHits(matches) ] <- star_introns$intron_motif[ subjectHits(matches) ]

    # Unannotated BUT detected intron -> novel intron
    containing <- findOverlaps(introns, star_introns, type = "within")
    overlapping <- findOverlaps(introns, star_introns, type = "any")

    if( verbose ) cat("OK. \n")
  }
  # END OF INTEGRATION
  #-------------------------------------------------------------------#

  #-------------------------------------------------------------------#
  ## KEEP ONLY INTRONS THAT DO NOT CONTAIN KNOWN EXON
  ## GROUP OVERLAPPING INTRONS INTO THEIR "INDEPENDENT FORM"

  annotated_exons <- exonsBy(txdb, by = "tx")
  annotated_exons <- unlist(annotated_exons)
  annotated_exons <- unique(annotated_exons)

  hits <- findOverlaps(query = annotated_exons, subject = introns, type = "within")

  introns$AnnotExonInside <- 0
  introns$AnnotExonInside[ subjectHits(hits) ] <- 1
  introns <- introns[ which(introns$AnnotExonInside != 1) ]

  #Take the intersection of introns from a same group
  intron_grouping <- reduce( introns )

  hits <- findOverlaps(query = intron_grouping, subject = introns)
  #All introns must be part of exactly one group:
  stopifnot( length(subjectHits(hits)) == length(introns) )

  introns$intron_group <- 0
  introns$intron_group[ subjectHits(hits) ] <- queryHits(hits)

  intron_grouping$intron_group <- 0
  intron_grouping$intron_group[ queryHits(hits) ] <- queryHits(hits)

  genes <- rtracklayer::import( gtf )
  genes <- genes[ genes$type == "gene" & genes$gene_id != "" ]

  hits <- findOverlaps(query = intron_grouping, subject = genes, type = "within")

  intron_grouping$gene_id <- intron_grouping$gene_name <- NA

  intron_grouping$gene_id[ queryHits(hits) ] <- genes$gene_id[ subjectHits(hits) ]
  intron_grouping$gene_name[ queryHits(hits) ] <- genes$gene_name[ subjectHits(hits) ]

  #-------------------------------------------------------------------#
  ## OUTPUT

  if( verbose ) cat("Saving results into", saveDir, "\n")

  fwrite(data.frame(intron_grouping), file = indIntronFile, sep = "\t")
  fwrite(data.frame(introns), file = allIntronFile, sep = "\t")

}
