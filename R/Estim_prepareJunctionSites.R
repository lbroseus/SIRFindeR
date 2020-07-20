# IRFindeR2: Intron Retention Estimation and Detection
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
#########################################################

star_junc_names  <- c("chr", "start", "end", "strand",
                      "intron_motif","annotated",
                      "uniq_map_cross_junc",
                      "mult_map_cross_jun",
                      "max_splice_align_overhang")

#_______________________________________________________#
#' Internal function for extracting junction counts from a STAR junction file
#'
#' @description This function extracts and merges junctions counts information (for each intron_group) from a STAR junction file.
#' It is intended to be used after calling ['defineIntrons()'], with the same \code{saveDir}.
#'
#' @import magrittr
#'
#' @importFrom rlang .data
#' @importFrom dplyr mutate select filter group_by summarise distinct arrange
#' @importFrom data.table fread fwrite
#' @importFrom rtracklayer import export
#'
#' @param gtf Path to a reference transcriptome in a gtf/gff file.
#' @param saveDir Path to a Project directory.
#' @param junctionFile Junction file output by STAR.
#' @param min_map_cross_junc Minimal number of supporting reads an intron detected from splice-aware alignments should have (Default: 3).
#'
#' @param indIntronsFile Path to a file with independent intron intervals as output by ['defineIntrons()'] (Default: NULL).
#'                       Need not be specified when \code{saveDir} is given and contains a file ""independent_introns.txt".
#' @param annotIntronsFile Path to a file with annotated intron intervals as output by ['defineIntrons()'] (Default: NULL).
#'                       Need not be specified when \code{saveDir} is given and contains a file ""all_introns.txt".
#'
#' @param verbose (DEFAULT: TRUE)
#'
#' @return Creates and saves files "junction_clusters.txt" and "junction_sites.txt" in \code{saveFolder}.
#'
#--------------------------------------------------------#

prepareJunctionSites <- function(gtf,
                                 saveDir,
                                 junctionFile,
                                 min_map_cross_junc = 5,
                                 indIntronsFile = NULL,
                                 annotIntronsFile = NULL,
                                 verbose = TRUE){

  ## Checks
  if( is.null(indIntronsFile) ){
    indIntronsFile <- paste(saveDir, "independent_introns.txt", sep = "/")
  }
  if( !file.exists(indIntronsFile) ) stop(paste("File", indIntronsFile, "does not exist. \n"))


  if( is.null(annotIntronsFile) ){
    annotIntronsFile <-   paste(saveDir, "all_introns.txt", sep = "/")
  }
  if( !file.exists(annotIntronsFile) ) stop(paste("File", annotIntronsFile, "does not exist. \n"))

  #--------------------------------------------------------------#

  ## Import data: annotated and novel introns

  Introns <- fread(file = indIntronsFile) %>% data.frame()

  Introns <- data.frame( Introns ) %>%
    dplyr::select(chr = seqnames, start, end, strand, intron_group)

  junctionSites <- fread( junctionFile ) %>% data.frame()
  colnames(junctionSites) <- star_junc_names

  junctionSites <- junctionSites %>%
    filter(min_map_cross_junc >= 5) %>%
    mutate(strand = ifelse(strand == 1, "+", "-")) %>%
    dplyr::select(-mult_map_cross_jun, -max_splice_align_overhang)

  ### Aggregate junction counts per intron cluster (SpliceLeft, SpliceRight, SpliceMax values)

  if( verbose ) cat("Extract Exon to Exon counts from junction file", junctionFile, "\n")

  # LEFT
  junctions <- junctionSites %>%
    group_by(chr, start, strand) %>%
    mutate(SpliceLeft = sum(uniq_map_cross_junc)) %>%
    data.frame()

  # RIGHT
  junctions <- junctions %>%
    group_by(chr, end, strand) %>%
    mutate(SpliceRight = sum(uniq_map_cross_junc)) %>%
    data.frame()

  junctions <- merge(Introns, junctions,
                     by = c("chr", "end", "start", "strand"),
                     all.x = T, all.y = T)

  junctions$SpliceRight[ is.na(junctions$SpliceRight) ] <- 0
  junctions$SpliceLeft[ is.na(junctions$SpliceLeft) ] <- 0

  junctions <- junctions %>%
    filter( !is.na(intron_group) ) %>%
    group_by(intron_group) %>%
    summarise(SpliceLeft = sum(SpliceLeft),
              SpliceRight = sum(SpliceRight)) %>%
    data.frame()

  junctions <- junctions %>% mutate(SpliceMax = ifelse(SpliceLeft>SpliceRight, SpliceLeft, SpliceRight))

  fwrite(junctions, file = paste(saveDir, "junction_clusters.txt", sep = "/"))

  ### STEP: prepare junction counts and ExonToIntron positions (used later for bootstrap)

  junctionLeft <- junctionSites %>%
                  group_by(chr, start, strand) %>%
                  summarise(SpliceCount = sum(uniq_map_cross_junc)) %>%
                  mutate(end = start) %>%
                  data.frame()

  junctionRight <- junctionSites %>%
                    group_by(chr, end, strand) %>%
                    summarise(SpliceCount = sum(uniq_map_cross_junc)) %>%
                    mutate(start = end) %>%
                    data.frame()

  # Turn it to SAF format GeneID Chr Start End Strand for featureCount
  junctionSites <- rbind.data.frame(junctionLeft, junctionRight) %>%
                arrange(chr, start) %>%
                mutate(GeneID = paste(chr, ":", start, ":", strand, sep = "")) %>%
                dplyr::select(GeneID, Chr = chr, Start = start, End = end, Strand = strand, 
                              SpliceCount) %>%
                mutate(Start = Start-1, End = End + 1)

  genes <- import( gtf )
  genes <- genes[ genes$type == "gene" & genes$gene_id != "" ]

  hits <- findOverlaps(query = GRanges(junctionSites), subject = genes, type = "within")

  junctionSites$gene_id <- junctionSites$gene_name <- NA

  junctionSites$gene_id[ queryHits(hits) ] <- genes$gene_id[ subjectHits(hits) ]
  junctionSites$gene_name[ queryHits(hits) ] <- genes$gene_name[ subjectHits(hits) ]
  junctionSites <- junctionSites[ -which(is.na(junctionSites$gene_id)), ]
  junctionSites <- distinct(junctionSites, .keep_all = T)

  fwrite(junctionSites, file = paste(saveDir, "junction_sites.txt", sep = "/"))

}
