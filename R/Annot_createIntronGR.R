# SIRFindeR: Intron Retention Estimation and Detection
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
#########################################################
#
# Extract intron intervals from a reference transcriptome
# Possibly add TSS sites information
#
#########################################################

#_______________________________________________________#
#' Deduce independent intron intervals from a transcriptome annotation
#'
#' @description This function extracts introns from a transcriptome annotation
#'
#' @import GenomicRanges
#'
#' @importFrom GenomicFeatures makeTxDbFromGRanges genes exonicParts
#' @importFrom S4Vectors Hits queryHits subjectHits
#' @importFrom rtracklayer import
#' @importFrom stringr str_remove
#'
#' @param gtf Path to a reference transcriptome in gtff/gff format
#'
#' @return A \code{GRanges} object with independent intron intervals
#'
#' @export
#'
#--------------------------------------------------------#

createIntronGR <- function(gtf){

  gtf <- import(con = gtf)

  col <- GenomicRanges::intersect(grep(colnames(mcols(gtf)), pattern = "transcript"),
                   grep(colnames(mcols(gtf)), pattern = "type"))

  colnames(mcols(gtf))[col] <- "transcript_type"

  ir_transcripts <- gtf[ which(gtf$transcript_type == "retained_intron") ]
  gtf <- gtf[ -which(gtf$transcript_type == "retained_intron") ]

  txdb <- suppressWarnings( makeTxDbFromGRanges(gr = gtf) )
  rm(gtf)

  #---------------------------------------------------------#
  #cat("Derive independent introns from annotation. \n")
  #---------------------------------------------------------#

  genes <- genes(x = txdb)
  exons <- exonicParts(txdb = txdb)
  exons <- reduce(exons)
  intronic_parts <- GenomicRanges::setdiff(genes, exons)
  rm(exons)

  #---------------------------------------------------------#
  # cat("Annotate independent introns with gene_ids. \n")
  #---------------------------------------------------------#

  hits <- findOverlaps(query = intronic_parts, subject = genes, type = "within")

  intronic_parts$gene_id <- ""
  intronic_parts$gene_id[ queryHits(hits) ] <- paste(intronic_parts$gene_id[ queryHits(hits) ], genes$gene_id[ subjectHits(hits) ], sep = "+")
  intronic_parts$gene_id[ queryHits(hits) ] <- str_remove(intronic_parts$gene_id[ queryHits(hits) ], pattern = "\\+")

  hits <- findOverlaps(query = intronic_parts, subject = ir_transcripts, type = "within")
  intronic_parts$annotated_intron <- 0
  intronic_parts$annotated_intron[ queryHits(hits) ] <- 1
  intronic_parts$annotated_intron <- as.factor(intronic_parts$annotated_intron)


  return( intronic_parts )

}
