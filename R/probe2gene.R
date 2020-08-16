#' Convert the EPIC methylation array probe-level DNA methylation to gene-level 
#'
#' This function extracts the methylation probes located in the promoter region and CpG islands for each gene and use their methylation levels to summerize the 
#' gene-level methylation by mean or median.
#' @param probeMat the probe-level DNA methylation matrix (probe by sample), entries are beta values.
#' @param nThread the number of threads to be used for the calculation. If NULL, no parallel computation will be performed.
#' @param aggr to which level should the probe be aggregated, need to be either 'gene' or 'transcript'. Default is 'gene'.
#' @param use specify how the multiple probe values be converted to the single gene-level value. Need to be either 'mean' or 'median'. Default is 'median'.
#' @param ... addtional parameter to makeCluster. If you used a Linux or MacOS, it is proper to add type='FORK'.
#'
#' @return the gene-level or transcript-level DNA methylation matrix (gene by sample), entries are beta values.
#' @export
#'
#' @examples hnsccSegfile <- system.file("extdata", "hnscc_GDC_methylation_part.tsv", package = "genomicWidgets")
#' hnsccMethy <- read.table(file = hnsccSegfile, sep = '\t', row.names = 1, header = TRUE)
#' probe2gene(probeMat=hnsccMethy) #no parallel computation
#' probe2gene(probeMat=hnsccMethy, nThread = 4) #parallel computation using 4 cores
probe2gene = function(probeMat, nThread = NULL, aggr = 'gene', use='median', ...){
  if(is.null(nThread))
    myLapply = lapply
  else{
    cl = makeCluster(nThread, ...)
    clusterExport(cl, "probeMat", envir = environment())
    myLapply = function(X, fun) parLapply(cl = cl, X=X, fun=fun)
  }
  if(aggr == 'gene') 
    aggrList = epicProbe2GeneMappingClean
  else if(aggr == 'transcript')
    aggrList = epicProbe2TranscriptMappingClean
  else
    stop('aggr should be either "gene" or "transcript"')
  
  if(use=='median')
    aggrLevel = myLapply(aggrList, function(x)apply(probeMat[x,],2,median,na.rm=T))
  else if(use=='mean')
    aggrLevel = myLapply(aggrList, function(x)apply(probeMat[x,],2,mean,na.rm=T))
  
  aggrLevel = do.call(rbind, aggrLevel)
  if(!is.null(nThread)) stopCluster(cl)
  aggrLevel
}
