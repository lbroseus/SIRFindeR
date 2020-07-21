# SIRFindeR: Stabilized estimation of Intron Retention levels
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
#########################################################

#_______________________________________________________#
#' Internal function for merging introns detected by STAR
#'
#' @description This junctions merges junction files
#'
#' @importFrom rlang .data
#' @importFrom dplyr group_by summarise mutate n
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
#' @importFrom GenomicRanges findOverlaps reduce setdiff 
#' @importFrom S4Vectors Hits queryHits subjectHits mcols
#' @importFrom data.table fread fwrite
#' @importFrom dplyr summarize group_by
#' @importFrom BiocGenerics width start end 
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
  curatedIntronFile <- paste(saveDir, "independent_introns.txt", sep = "/")
  # Annotated Independent intronic intervals to merge samples in an experiment
  pseudointronFile <- paste(saveDir, "pseudointrons.txt", sep = "/")
  # Augmented annotation of all introns
  allIntronFile <-   paste(saveDir, "all_introns.txt", sep = "/")
  
  #-------------------------------------------------------------------#
  ## IMPORT ANNOTATED INTRONS
  
  if( verbose ) cat("Importing transcriptome annotation...\n")
  
  pseudowidth.thr <- 20
  
  ref <- rtracklayer::import(gtf)
  colnames(S4Vectors::mcols(ref)) <- str_replace(string = colnames(S4Vectors::mcols(ref)),
                                                 pattern = "biotype", replacement = "type")
  ref <- ref[ !ref$transcript_type %in% c("retained_intron") ]
  
  txdb <- suppressWarnings( makeTxDbFromGRanges(ref) )
  
  # INTRONS ANNOTES
  introns <- intronsByTranscript(txdb, use.names = TRUE)
  introns <- unlist(introns)
  introns <- unique(introns)
  introns$Annotated <- 1
  
  #--------------------------------------------------------#
  # NOUVEAUX INTRONS
  
  if( verbose ) cat("Integrating novel introns...\n")
  
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
  
  index <- union(subjectHits(matches), subjectHits(containing))
  index <- setdiff(subjectHits(overlapping), index) 
  novel_introns <- star_introns[ index ]
  
  # Filter new intervals
  novel_introns <- novel_introns[novel_introns$uniq_map_cross_junc>min_map_cross_junc 
                                 & width(novel_introns)>min_novel_intron_width]
  
  stopifnot( length( findOverlaps(introns, novel_introns, type = "equal") ) ==0 )
  
  novel_introns$Annotated <- novel_introns$annotated*2
  novel_introns$annotated <- NULL
  novel_introns$mult_map_cross_jun <- NULL
  novel_introns$max_splice_align_overhang <- NULL
  novel_introns$Detected <- 1
  
  # All anotated and/or detectable introns
  introns <- c(introns, novel_introns)
  
  # PSEUDOINTRONS
  # used as main reference (to identify features common to all samples in an study)
  
  if( verbose ) cat("Defining (independent) pseudointrons...\n")
  
  pseudointrons <- intronicParts(txdb = txdb, linked.to.single.gene.only = TRUE)
  pseudointrons <- pseudointrons[ width(pseudointrons) >= pseudowidth.thr ]
  
  ### OVERLAPPING ANNOTATED EXON
  annotated_exons <- exonsBy(txdb, by = "tx")
  annotated_exons <- unlist(annotated_exons)
  annotated_exons <- unique(annotated_exons)
  
  ## Tag introns with known exon
  hits <- findOverlaps(query = annotated_exons, subject = introns, type = "within")
  introns$AnnotExonInside <- 0
  introns$AnnotExonInside[ subjectHits(hits) ] <- 1
  introns <- introns[ which(introns$AnnotExonInside != 1) ]
  
  ## Remove intronicParts containing an exon
  hits1 <- findOverlaps(query = annotated_exons, subject = pseudointrons, type = "within")
  index <- subjectHits(hits1)
  pseudointrons <- pseudointrons[ -index ]
  #Overlap with exon but not a true intron
  hits1 <- findOverlaps(query = pseudointrons, subject = annotated_exons)
  hits2 <- findOverlaps(query =  pseudointrons, subject = introns, type = "equal")
  index <- setdiff(queryHits(hits1), queryHits(hits2))
  pseudointrons <- pseudointrons[ -index ]
  
  #INTRON AVEC ou SANS PARTIE IDENTIFIABLE:
  hits <- findOverlaps(query = pseudointrons, subject = introns, type = "within")
  introns$Identifiable <- 0
  introns$Identifiable[ subjectHits(hits) ] <- 1
  
  ## Record intron index in gene
  introns <- introns %>% arrange(seqnames, start, end) %>% 
    group_by(gene_id, gene_name) %>% mutate(intron_number = 1:n()) %>% data.frame()
  
  table(introns$Identifiable)
  
  #DETECTER LES PSEUDOINTRONS COMPLEXES
  hits <- findOverlaps(query = novel_introns, subject = pseudointrons, type = "within")
  hits <- data.frame(hits) %>% group_by(subjectHits) %>% summarize(nbHits = n()) %>% data.frame()
  
  pseudointrons$nbEvents <- 0
  pseudointrons$nbEvents[ hits$subjectHits ] <- hits$nbHits
  
  table(pseudointrons$nbEvents)
  
  #Take the intersection of introns from a same group
  #novel_intron already defined
  if( verbose ) cat("Refining pseudointrons...\n")
  
  novel_introns <- novel_introns[ novel_introns$Annotated == 0 ]
  
  curateIntron <- function(granges){
    
    hits <- findOverlaps(granges, novel_introns, type = "within")
    
    S4Vectors::mcols(granges) <- NULL
    
    if(length(hits)>0){
      
      unchanged <- granges[ -unique(queryHits(hits)) ]
    
      changed <- granges[ queryHits(hits) ]
      hits <- findOverlaps(changed, novel_introns, type = "within")
      
      pool <- unique( queryHits(hits) )
      for(p in pool){
        
        novel <- subjectHits(hits)[queryHits(hits) == p]
        novel <- novel_introns[ novel ]
        shrunknovel <- novel[1]
        BiocGenerics::start(shrunknovel) <- max(BiocGenerics::start(novel))
        BiocGenerics::end(shrunknovel) <- min(BiocGenerics::end(novel))
        
        S4Vectors::mcols(shrunknovel) <- NULL
        changed[ p ] <- shrunknovel
        
      }
      changed <- c(unchanged, changed)
      changed <- sort(changed)
      
      return( changed )
    
    }else return(granges)
    
  }
  
  curatedIntrons <- curateIntron(pseudointrons)
  
  hits <- findOverlaps(query = curatedIntrons, subject = pseudointrons, type = "within")
  
  pseudointrons$identifiable <- 0
  pseudointrons$identifiable[ subjectHits(hits) ] <- 1
  
  table(pseudointrons$identifiable)
  
  # Les introns non-identifiables sont exclus de l'analyse
  # Les introns identifiables peuvent avoir plusieurs parties identifiables?
  
  pseudointrons$intron_group <- 0
  pseudointrons$intron_group[ subjectHits(hits) ] <- queryHits(hits)
  S4Vectors::mcols(pseudointrons)$tx_id <- NULL
  
  curatedIntrons$intron_group <- 0
  curatedIntrons$intron_group[ queryHits(hits) ] <- queryHits(hits)
  
  genes <- rtracklayer::import( gtf )
  genes <- genes[ genes$type == "gene" & genes$gene_id != "" ]
  
  hits <- findOverlaps(query = curatedIntrons, subject = genes, type = "within")
  
  curatedIntrons$gene_id <- curatedIntrons$gene_name <- NA
  
  curatedIntrons$gene_id[ queryHits(hits) ] <- genes$gene_id[ subjectHits(hits) ]
  curatedIntrons$gene_name[ queryHits(hits) ] <- genes$gene_name[ subjectHits(hits) ]
  
  #-------------------------------------------------------------------#
  ## OUTPUT
  
  if( verbose ) cat("Saving results into", saveDir, "\n")
  
  fwrite(data.frame(curatedIntrons), file = curatedIntronFile, sep = "\t")
  fwrite(data.frame(pseudointrons), file = pseudointronFile, sep = "\t")
  fwrite(data.frame(introns), file = allIntronFile, sep = "\t")
  
}
