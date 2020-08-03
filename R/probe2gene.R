#' Convert the EPIC methylation array probe-level DNA methylation to gene-level 
#'
#' This function extracts the methylation probes located in the promoter region and CpG islands for each gene and use their methylation levels to summerize the 
#' gene-level methylation by mean or median.
#' @param probeMat the probe-level DNA methylation matrix (probe by sample), entries are beta values.
#' @param nThread the number of threads to be used for the calcluation. If NULL, the number will be equal to the (# of CPU threads)-1. Default is NULL.
#' @param use specify how the multiple probe values be converted to the single gene-level value. Need to be either 'mean' or 'median'.
#' @param ... addtional parameter to for makeCluster. If you used a Linux or MacOS, it is proper to add type='FORK'.
#'
#' @return the gene-level DNA methylation matrix (gene by sample), entries are beta values.
#' @export
#'
#' @examples hnsccSegfile <- system.file("extdata", "hnscc_GDC_methylation_part.tsv", package = "genomicWidgets")
#' hnsccMethy <- read.table(file = hnsccSegfile, sep = '\t', row.names = 1, header = TRUE)
#' probe2gene(probeMat=hnsccMethy)
probe2gene = function(probeMat, nThread = NULL, use='median', ...){
  if(is.null(nThread))
    nThread = detectCores() - 1
  cl = makeCluster(nThread, ...)
  clusterExport(cl, "probeMat", envir = environment())
  if(use=='median')
    genelevel = parLapply(epicProbe2GeneMappingClean, function(x)apply(probeMat[x,],2,median,na.rm=T))
  else if(use=='mean')
    genelevel = parLapply(epicProbe2GeneMappingClean, function(x)apply(probeMat[x,],2,mean,na.rm=T))
  genelevel = do.call(rbind, genelevel)
  stop(cl)
  genelevel
}
