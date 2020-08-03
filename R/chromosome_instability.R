#' Calculate the chromosome instability from the segment-level SCNA data
#'
#' @param segDf segment-level SCNA data frame. This is usually the same input as GISTIC2. Make sure that there are these columns: "sample", "chromosome", "start", "end", "log2". 
#' @param genomeVersion the genome version by which SCNA data were generated, should be either 'hg19' or 'hg38'. Default is 'hg38'.  
#' @param option choose how to calculate chromosome instability, should be one of 'abs', 'amp', 'del' or 'origin'.
#' If 'abs' (default), both CN deletions and amplifications will contribute to the instability score.
#' If 'amp', only CN amplifications will be considered. If 'del', only CN deletions will be considered.
#' If 'origin', CN deletions and amplifications will be canceled out for each other, compared to 'abs' where the absolute CN deletions and CN amplifications will be added up.
#' @param chrDf a data frame specifying what genomic ranges should be considered for the instability score calculation.
#' Three columns should be specified: "chromosome", "start" and "end. If NULL (default), the whole chromosome will be used to calculated the instability score.
#' @para nThread the number of CPU cores used for the analysis, default is (# of cores - 1)
#' @param ... additional parameter to for makeCluster. If you used a Linux or MacOS, it is proper to add type='FORK'.
#' @return a data frame with columns representing each samples in segDf, and rows representing instability scores specified by   the input chrDf  
#' @export
#' @examples
#' hnsccSegfile <- system.file("extdata", "CPTAC3_HNSCC_SCNA_segment_level.tsv", package = "genomicWidgets")
#' hnsccSegDf <- read.table(hnsccSegfile, header = TRUE)
#' weightAveChr(hnsccSegDf, genomeVersion = 'hg38') 
#'  
weightAveChr = function(segDf, genomeVersion='hg38', option = "abs", chrDf=NULL, nThread = NULL, ...){
  if(is.null(nThread)) nThread = nThread = detectCores() - 1
  if(genomeVersion=='hg38') chrLength = hg38ChrLengthDf
  else if (genomeVersion=='hg19') chrLength = hg19ChrLengthDf
  else stop("Current only support genome version hg19 and hg38")
  
  #harmonize the chromosome to be the style of "Y"
  segDf$chromosome = as.character(segDf$chromosome)
  if(grepl(pattern = "(chr)|(Chr)",segDf$chromosome[1]))
    segDf$chromosome = substr(segDf$chromosome, 4, nchar(segDf$chromosome))
  segDf$chromosome[segDf$chromosome == "23"] = "X"
  segDf$chromosome[segDf$chromosome == "24"] = "Y"
  
  if(is.null(chrDf))
    chrDf=data.frame(chromosome=chrLength$chromosome, start=0, end=chrLength$length)
  segDfList = split(segDf, f = segDf$sample)
  changePerChr = list()
  cl = makeCluster(nThread, ...)
  clusterExport(cl, c("segDfList",".weightAveSeg","chrDf","option"))
  for(idx in 1:nrow(chrDf)){
    changePerChr[[chrDf$chromosome[idx]]] = parLapply(cl = cl, X = segDfList, function(x).weightAveSeg(segDf = x, 
                                                                                      chr = chrDf$chromosome[idx], 
                                                                                      start = chrDf$start[idx], 
                                                                                      end = chrDf$end[idx],
                                                                                      chrLength = chrLength,
                                                                                      option = option))
  }
  changePerChr = lapply(changePerChr, unlist)
  changePerChr = do.call(rbind, changePerChr)
  stop(cl)
  return(changePerChr)
}

.weightAveSeg = function(segDf, 
                        chr, 
                        start = NULL, 
                        end = NULL, 
                        option = "abs",
                        chrLength){

  if(is.null(start)){
    start = 0
  }
  if(is.null(end)){
    end = chrLength[chr, "length"]
  }
  #remove na rows
  segDf = segDf[complete.cases(segDf),]
  segDf$chromosome = gsub(segDf$chromosome, pattern = "chr", replacement = "")
  chr = gsub(chr, pattern = "chr", replacement = "")
  
  refRangeObj = GRanges(seqnames = chr, 
                        ranges = IRanges(start = start, end = end),
                        strand = Rle(strand("*"), 1))
  
  #take all the segments for a single sample and convert to a range Obj
  rangeObj = GRanges(seqnames = segDf$chromosome, 
                     ranges = IRanges(start = segDf$start, end = segDf$end),
                     strand = Rle(strand("*"), nrow(segDf)),
                     log2 = segDf$log2)
  #subset the specific chr
  rangeObj.sub = rangeObj[seqnames(rangeObj) == chr]
  
  #need to throw out warning: if refRangeObj is smaller than rangeObj
  overlapIndex = findOverlaps(query = rangeObj.sub,subject = refRangeObj, type = 'any')
  overlapIndex = as.data.frame(overlapIndex)
  overlapRegion = overlapsRanges(ranges(rangeObj.sub), ranges(refRangeObj))
  overlapRegion = as.data.frame(overlapRegion)
  overlapRegion$log2 = rangeObj.sub$log2[overlapIndex$queryHits]
  if(option == "abs"){
    overlapRegion$log2 = abs(overlapRegion$log2)
    totalWidth = sum(overlapRegion$width)
    weightAve = sum(overlapRegion$width * overlapRegion$log2)/totalWidth
  }
  else if(option == "amp"){
    #only summerize the amplifcation
    amp = overlapRegion$log2 >0
    if(sum(amp) == 0){
      weightAve = 0
    }else{
      totalWidth = sum(overlapRegion$width)
      weightAve = sum(overlapRegion$width[amp] * overlapRegion$log2[amp])/totalWidth
    }
  }else if(option == "del"){
    #only summerize the deletion
    del = overlapRegion$log2 <0
    if(sum(del) == 0){
      weightAve = 0
    }else{
      totalWidth = sum(overlapRegion$width)
      weightAve = sum(overlapRegion$width[del] * overlapRegion$log2[del])/totalWidth
    }
  }else if(option == "origin"){
    totalWidth = sum(overlapRegion$width)
    weightAve = sum(overlapRegion$width * overlapRegion$log2)/totalWidth
  }
  #cat(paste0("done with chromosome ",chr," for sample ",segDf$sample[1], "\n"))
  weightAve                  
}


