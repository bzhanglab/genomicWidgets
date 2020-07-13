#' plot SCNA landscape as heatmap
#'
#' plotCNV is used to plot the segment-level SCNA data into heatmap
#' @param segDf segment-level SCNA data frame. This is usually the same input as GISTIC2. Make sure there are these columns "sample", "chromosome", "start", "end", "log2".
#' @param geneVersion optionally specify the genome version of SCNA data were generated, should be either 'hg19' or 'hg38'. If this option is specified, the function will fill any gaps between the data and the whole genome. 
#' If geneVersion is NULL, the no gap will be filled. Default is NULL.
#' @param sampleOrder a optional character vector to specify the sample order for the display. If NULL, the samples will be ranked as their order in the segDf. Default is NULL. 
#' @param chrOrder a optional character vector to specify the chromosome order for the display. If NULL, the chromosomes will be ranked as their order in the segDf. Default is NULL. 
#' @return a ggplot2 plot object
#' @examples 
#' hnsccSegfile <- system.file("extdata", "CPTAC3_HNSCC_SCNA_segment_level.tsv", package = "genomicWidgets")
#' hnsccSegDf <- read.table(hnsccSegfile)
#' plotCNV(hnsccSegDf, geneVersion = 'hg38')
#' @export
plotCNV = function(segDf, geneVersion = NULL, sampleOrder = NULL, chrOrder = NULL){
  require(ggplot2, quietly = TRUE)
  require(scales, quietly = TRUE)
  if(!all(c("sample","chromosome","start","end","log2") %in% colnames(segFile))){
    stop("Error: the segFiles should contain 5 columns:\n
         \t\"sample\",\"chromosome\",\"start\",\"end\",\"log2\"")
  }
  
  #fill the gap as NAs, as required
  if(!is.null(geneVersion)){
    if(geneVersion=='hg38')
      chrLengthObj = .seg2Ranges(hg38ChrLengthDf)
    else if(geneVersion=='hg19')
      chrLengthObj = .seg2Ranges(hg19ChrLengthDf)
    else
      stop("Current only support genome version hg19 and hg38")
    segDf = .fillwithNA.2(segDf = segFile, totalChrRangeObj = chrLengthObj)
  }
  
  if(!is.factor(segDf$sample)){
    if(is.null(sampleOrder)){
      segDf$sample=factor(segDf$sample, levels = unique(segDf$sample))
    }else{
      segDf$sample=factor(segDf$sample, levels = sampleOrder)
    }
  }
  if(!is.factor(segDf$chromosome)){
    if(is.null(chrOrder)){
      segDf$chromosome=factor(segDf$chromosome, levels = unique(segDf$chromosome))
    }else{
      segDf$chromosome=factor(segDf$chromosome, levels = chrOrder)
    }
  }
  segDf$name = paste(segDf$sample, segDf$chromosome, segDf$start, sep = '-')
  segDf$width = segDf$end - segDf$start
  segDf$log2 = ifelse(segDf$log2 > 1, 1, segDf$log2)
  segDf$log2 = ifelse(segDf$log2 < -1, -1, segDf$log2)
  plotCNV = ggplot(segDf, aes(1, width, fill=log2)) + 
    geom_bar(stat="identity")+
    scale_fill_gradient2(high = 'red', low = 'blue')+
    geom_hline(yintercept = 0, size = 1)+
    scale_y_continuous(expand = c(0,0))+
    facet_grid(sample ~ chromosome, space = "free",scales = "free") +
    theme(axis.title =element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.spacing = unit(0, "lines"),
          strip.background = element_rect(colour="black", fill="white", 
                                          size=1, linetype="blank"),
          strip.text.y = element_text(angle = 0),
          panel.background = element_blank(),
          panel.grid = element_blank())+
    coord_flip()
  plotCNV
}


#' fill the gaps for segment-level SCNA data with NA
#'
#' fill the gaps for single sample's segment-level SCNA data with NA
#' @param segDf segment-level SCNA data frame. This is usually the same input as GISTIC2. Make sure there are these columns "sample", "chromosome", "start", "end", "log2".
#' @param totalChrRangeObj a GenomicRanges object, specifying the lengths of all chromosomes
#'
#' @return a data frame with gap-filled segment-level SCNA data
#'
.fillwithNA = function(segDf, totalChrRangeObj){

  segDf$chromosome = as.character(segDf$chromosome)
  #harmonize the chromosome to be the style of "Y"
  if(grepl(pattern = "chr",segDf$chromosome[1]))
    segDf$chromosome = substr(segDf$chromosome, 4, nchar(segDf$chromosome))
  segDf$chromosome[segDf$chromosome == "23"] = "X"
  segDf$chromosome[segDf$chromosome == "24"] = "Y"
  
  segObj = GRanges(seqnames = segDf$chromosome, 
                   ranges = IRanges(start = segDf$start, end = segDf$end),
                   strand = Rle(strand("*"), nrow(segDf)),
                   log2 = segDf$log2)
  ##the segmenets not covered by the segObj
  naRanges = setdiff(totalChrRangeObj, segObj)
  if(length(naRanges)!=0){
    naRanges$log2 = NA
    #combining naRanges to segObj
    segObj = c(segObj, naRanges)
  }
  segObj = sort(segObj)
  segObj = as.data.frame(segObj)
  colnames(segObj)[1] = "chromosome"
  segObj
}


#' fill the gaps for segment-level SCNA data with NA
#'
#' fill the gaps for multiple samples' segment-level SCNA data with NA
#' @param segDf segment-level SCNA data frame. This is usually the same input as GISTIC2. Make sure there are these columns "sample", "chromosome", "start", "end", "log2".
#' @param totalChrRangeObj a GenomicRanges object, specifying the lengths of all chromosomes
#'
#' @return a data frame with gap-filled segment-level SCNA data
#'
.fillwithNA.2 = function(segDf, totalChrRangeObj){
  #this version works for multiple samples
  if(length(unique(segDf$sample)) == 1)
    fillwithNA(segDf, totalChrRangeObj)
  else{
    segList = split(segDf, f = segDf$sample)
    segList = lapply(segList, function(x)fillwithNA(segDf = x, totalChrRangeObj = totalChrRangeObj))
    newSegDf = do.call(rbind, segList)
    newSegDf$sample = rep(names(segList), sapply(segList, nrow))
    newSegDf
  }
}


#' convert a data frame of chromosome length to GRanges Obj
#'
#' @param totalChrLengthDf 
#'
#' @return a GenomicRanges Obj
#'
#' @examples
.seg2Ranges = function(totalChrLengthDf){
  require(GenomicRanges, quietly = TRUE)
  rangeObj = GRanges(seqnames = totalChrLengthDf$chromosome, 
                     ranges = IRanges(start = 0, end = totalChrLengthDf$length),
                     strand = Rle(strand("*"), nrow(totalChrLengthDf)))
  rangeObj
}

