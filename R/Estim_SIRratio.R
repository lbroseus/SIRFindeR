# IRFindeR2: Intron Retention Estimation and Detection
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
#########################################################


#_______________________________________________________#
#' Internal function to estimate the SIRratio
#'
#' @description This function merges intronic and exonic
#' abundances and compute the IRratio2.
#'
#' @importFrom dplyr mutate
#' @importFrom data.table fread fwrite
#'
#' @param saveDir Path to the directory where sample results are saved.
#' @param readLength Total read length.
#'
#' @return Saves a \code{data.frame} containing IRratio2 estimates
#' for each intron group in a file named "IRratio2.txt" in \code{saveDir}.
#'
#' @seealso ['computeIRratio2()']
#--------------------------------------------------------#
#

estimateSIRratio <- function(saveDir, readLength){

  # Details to be carefully checked
  Rcpp::cppFunction('double computeMeanCoveredArea(int& intronWidth, int& readLength) {

    int nbSP( intronWidth + readLength - 1);
    int trim_right(0);
    int trim_left(0);

    double meanArea(0);

    for(int sp = 0; sp < nbSP+1; ++sp){

      trim_right = sp+readLength-1 - (readLength-1+intronWidth);
      if( trim_right<0 ) trim_right = 0;

      trim_left = sp - readLength + 1;
      if( trim_left>0 ) trim_left = 0;

      meanArea += readLength + trim_left - trim_right ;
    }

    meanArea /= nbSP;

    return meanArea;

  }')

  indIntronsFile <- paste(saveDir, "independent_introns.txt", sep = "/")
  intronicCountsFile <- paste(saveDir, "intronic_counts.txt", sep = "/")
  junctionCountsFile <- paste(saveDir, "junction_clusters.txt", sep = "/")
  bootStatsFile <- paste(saveDir, "bootstrapped_gene_coverage.txt", sep = "/")

  resultsFile <- paste(saveDir, "SIRratio.txt", sep = "/")

  indIntrons <- fread(indIntronsFile, data.table = F)
  # Merge intronic counts, solice counts and bootstrap estimates
  counts <- fread(intronicCountsFile, data.table = F)
  indIntrons <- merge(indIntrons, counts, by = "intron_group")

  counts <- fread(junctionCountsFile, data.table = F)
  indIntrons <- merge(indIntrons, counts, by = "intron_group")

  counts <- fread(bootStatsFile, data.table = F)
  indIntrons <- merge(indIntrons, counts, by = "gene_id", all.x = T, all.y = F)

  indIntrons <- indIntrons %>%
                mutate(cv = ifelse(is.na(boot_var) | boot_var*boot_mean == 0, 0,
                       log(boot_var/boot_mean+ exp(1))),
                       nSites = ifelse(is.na(nSites),0, nSites),
                       boot_mean = ifelse(is.na(boot_mean),0, boot_mean))

  indIntrons$meanArea <- tapply(indIntrons$width, 1:nrow(indIntrons), function(x) computeMeanCoveredArea(x, readLength))
  indIntrons <- indIntrons %>% mutate(intronicDepth = intronicCount/(floor(width/meanArea)+1))

  #indIntrons <- indIntrons %>% mutate(boot_mean = boot_mean-intronicDepth)
  #indIntrons <- indIntrons %>% mutate(boot_mean = ifelse(boot_mean<0, 0, boot_mean))

  indIntrons <- indIntrons %>% mutate(lambda = nSites/(nSites+1)*abs(boot_mean-intronicDepth)/(boot_mean+intronicDepth+1))

  indIntrons <- indIntrons %>% mutate(exonicAbundance = lambda*boot_mean + (1-lambda)*SpliceMax)

  indIntrons <- indIntrons %>% mutate(width = end - start + 1,
                                      SIRratio = ifelse(intronicCount+exonicAbundance == 0,
                                                        0,  intronicCount/(intronicCount+exonicAbundance)))

  indIntrons <- indIntrons %>% mutate(SIRratio = SIRratio/(floor(width/readLength)+1 - floor(width/readLength)*SIRratio))

  fwrite(indIntrons, file = resultsFile, sep = "\t")

}
