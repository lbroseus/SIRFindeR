# IRFindeR2: Intron Retention Estimation and Detection
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
#########################################################


#_______________________________________________________#
#' Internal function appendIntronicCounts()
#'
#' @description This function counts reads overlapping
#' independent introns intervals using ['Rsubread::featureCount()'].
#' These counts will be merged to intron intervals provided in file "curated_introns.txt".
#'
#' @import magrittr
#'
#' @importFrom dplyr mutate group_by summarise
#' @importFrom data.table fread fwrite
#' @importFrom Rsubread featureCounts
#'
#' @param bamFile Path to read alignments in a bam file.
#' @param saveDir Name for the Project Directory will be saved.
#'        The directory should contain a file "curated_introns.txt", as output by function ['curateIntrons()'].
#' @param libraryType Read type. Either "SE" (single-end) or "PE" (paired-end)
#' @param verbose logical indicating if verbose information for ['featureCount()'] debugging will be generated.
#'
#'
#' @export
#'
#' @return Intron read counts in a file names "intronic_counts.txt", saved in directory \code{saveDir}.
#--------------------------------------------------------#
#

appendIntronicCounts <- function(bamFile, saveDir, libraryType, verbose = FALSE){

  ### Checks

  if( !file.exists(bamFile) ) stop("Input bam file does not exist...\n")

  if( !libraryType %in% c("SE", "PE") ) stop("Invalid libraryType... \n")

  curatedIntronsFile <-  paste(saveDir, "curated_introns.txt", sep = "/")
  if( !file.exists(curatedIntronsFile) ) stop("Input curatedIntronsFile does not exist...\n")

  ###

  intronicCountsFile <-  paste(saveDir, "intronic_counts.txt", sep = "/")

  curatedIntrons <- fread(file = curatedIntronsFile, data.table = F)
  colnames(curatedIntrons) <- c("Chr", "Start", "End", "width", "Strand", "intron_group")

  curatedIntrons <- curatedIntrons %>% mutate(GeneID = paste(Chr, ":", Start, ":", Strand, sep = ""))

  counts <- featureCounts(files = bamFile,
                          # Intervals to quantify
                          isGTFAnnotationFile = FALSE,
                          annot.ext = curatedIntrons,
                          minOverlap = 3,
                          # Counting parameters
                          primaryOnly = TRUE,
                          strandSpecific = 1,
                          isPairedEnd = (libraryType=="PE"),
                          allowMultiOverlap = FALSE,
                          countChimericFragments = FALSE,
                          verbose = verbose)

  stopifnot( nrow(curatedIntrons) ==  length(as.vector( counts$counts )) )

  curatedIntrons$count <- as.vector( counts$counts )

  curatedIntrons <- curatedIntrons %>%
    group_by(intron_group) %>%
    summarise(intronicCount = sum(count)) %>%
    data.frame()

  fwrite(curatedIntrons, file = intronicCountsFile, sep = "\t")

}


#_______________________________________________________#
#' Internal function bootstrapSample()
#'
#' @description This function simply creates a bootstrap sample from a vector of counts,
#' and computes the sample mean and variance.
#'
#' @param counts Vector of counts to sample from.
#' @param sboot Size of bootstrap samples.
#'
#' @return A 2-vector with the mean and variance of a sample of size \code{sboot} drawn with replacement in \code{counts}.
#--------------------------------------------------------#

bootstrapSample <- function(counts, sboot){

  echantillon <- sample(x = counts, size = sboot, replace = T)

  return( c( mean(echantillon), var(echantillon) ) )

}

#_______________________________________________________#
#' Internal function bootstrapMean()
#'
#' @description This function does something useful
#'
#' @param SpliceCounts Vector with exon to exon counts.
#' @param ExonToIntronCounts Vector with exon to intron counts. Should have same size as \code{SpliceCounts}.
#' @param nboot Number of bootstrap sample to draw (Default: 100).
#' @param sboot Size of bootstrap samples (Default 50).
#'
#' @return
#--------------------------------------------------------#

bootstrapMean <- function(SpliceCounts, ExonToIntronCounts, nboot = 100, sboot = 50){

  if( length(SpliceCounts)<2 ) return( data.frame(X1 = SpliceCounts+ExonToIntronCounts, X2 = 0 ) )
  else{

    boot_samples <- replicate(n =  nboot,
                              simplify = TRUE,
                              expr = bootstrapSample(SpliceCounts+ExonToIntronCounts, sboot))

    return( data.frame( t(rowMeans(boot_samples)), nSites = length(SpliceCounts) ) )
  }
}

#_______________________________________________________#
#' Internal function bootstrapJunctionCounts()
#'
#' @description This function computes a bootstrap estimate of gene expression by using exon-exon junction counts.
#' It is not intended to be used alone, but as a part of the IRratio2 estimation procedure,
#' wrapped in function ['computeIRratio2()'].
#'
#' @import magrittr
#'
#' @importFrom dplyr mutate filter group_by summarise do select
#' @importFrom data.table fread fwrite
#' @importFrom Rsubread featureCounts
#'
#' @param bamFile Path to read alignments in a bam file.
#' @param saveDir Working project directory.
#' @param libraryType Type of the library of RNA-seq data. Either "SE" (single-end) or "PE" (paired-end).
#' @param nboot Number of bootstrap sample to draw (Default: 100).
#'
#' @param verbose logical indicating if verbose information for ['featureCount()'] debugging will be generated.
#'
#' @seealso ['featureCount()'], ['computeIRratio2()']
#'
#' @return
#'
#--------------------------------------------------------#

bootstrapJunctionCounts <- function(bamFile, saveDir, libraryType, nboot = 100, verbose = FALSE){

  ### Checks

  if( !libraryType %in% c("SE", "PE") ) stop("Invalid libraryType... \n")

  ## Needed input files

  if( !file.exists(bamFile) ) stop("Input bam file does not exist...\n")

  junctionSitesFile <- paste(saveDir, "junction_sites.txt", sep = "/")
  if( !file.exists(junctionSitesFile) ) stop("Input junctionSitesFile file does not exist...\n")

  ## Output files

  junctionCountsFile <- paste(saveDir, "junction_counts.txt", sep = "/")

  bootFile <- paste(saveDir, "bootstrapped_gene_coverage.txt", sep = "/")

  ### Count Exon To Intron Reads

  if( verbose ) cat("Count Exon to Intron junction reads from bam file", bamFile, "\n")

  junctionSites <- fread(junctionSitesFile, data.table = F)

  junctionSites <- junctionSites %>% mutate(GeneID = paste(Chr, ":", Start, ":", Strand, sep = ""))

  ExonToIntronCounts <-  featureCounts(files = bamFile,
                                       # Intervals to quantify (width = 3)
                                       isGTFAnnotationFile = FALSE,
                                       annot.ext = junctionSites,
                                       minOverlap = 3,
                                       # Counting parameters
                                       primaryOnly = TRUE,
                                       strandSpecific = 1,
                                       isPairedEnd = (libraryType=="PE"),
                                       allowMultiOverlap = TRUE,
                                       countChimericFragments = FALSE,
                                       verbose = verbose)

  ExonToIntronCounts <- data.frame(ExonToIntronCounts$annot,
                                   ExonToIntronCount = as.vector(ExonToIntronCounts$counts))

  junctionSites <- merge(junctionSites, ExonToIntronCounts, by = c("GeneID", "Chr", "Start", "End", "Strand"))
  junctionSites <- junctionSites %>% dplyr::select(Chr, Start, End, Strand, gene_name, gene_id, SpliceCount, ExonToIntronCount)

  ### Save junction sites with corresponding spliced and exon-to-intron counts

  fwrite(junctionSites, file = junctionCountsFile)

  ### Bootstrap Splice and ExonToIntron counts to estimate gene coverage

  if( verbose ) cat("Compute gene coverage estimate by bootstrap \n")

  system.time(

    boot <- junctionSites %>%
      group_by(gene_id) %>%
      mutate(zeroGene = (sum(SpliceCount)+sum(ExonToIntronCount)==0), nbSites = n()) %>%
      filter(!zeroGene) %>%
      do(bootstrapMean(.$SpliceCount, .$ExonToIntronCount, nboot)) %>%
      data.frame()

  )

  colnames(boot) <-  c("gene_id", "boot_mean", "boot_var", "nSites")

  boot <- boot %>% filter(gene_id != "")

  fwrite(boot, file = bootFile)

}
