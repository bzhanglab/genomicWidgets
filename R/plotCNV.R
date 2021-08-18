#' Plot SCNA landscape as heatmap
#'
#' plotCNV is used to plot the segment-level SCNA data into heatmap.
#' @param segDf segment-level SCNA data frame. This is usually the same input as GISTIC2. Make sure that there are these columns "sample", "chromosome", "start", "end", "log2".
#' @param genomeVersion specify the genome version by which SCNA data were generated, should be either 'hg19' or 'hg38'.
#' If genomeVersion is NULL, the no gap will be filled. Default is NULL.
#' @param samples a optional character vector to specify the samples for the plot.
#' If NULL, all the samples will be ploted and ranked as their order in the segDf. Default is NULL. 
#' @param chr a optional character vector to specify the chromosomes for the plot.
#' If NULL, the chromosomes in the segDf will be all ploted and ranked as their order in the segDf. Default is NULL. 
#' @param start a integer vector specifying the start positions of each chromosome.
#' If NULL, the position will be the begining of the chromosome(s). Default is NULL.
#' @param end a integer vector specifying the end positions of each chromosome. 
#' If NULL, the position will be the end of the chromosome(s). Default is NULL.
#' @return a ggplot2 plot object
#' @examples 
#' hnsccSegfile <- system.file("extdata", "CPTAC3_HNSCC_SCNA_segment_level.tsv", package = "genomicWidgets")
#' hnsccSegDf <- read.table(hnsccSegfile, header = TRUE)
#' plotCNV(hnsccSegDf, genomeVersion = 'hg38') #plot the whole CNV landscape for the cohort
#' plotCNV(hnsccSegDf, genomeVersion = 'hg38', samples = c("C3L-00995","C3L-00977","C3L-00987","C3L-00994")) #plot SCNA heatmap for selected samples 
#' plotCNV(hnsccSegDf, genomeVersion = 'hg38', samples = c("C3L-00977","C3L-00987","C3L-00994","C3L-00995"), chr=c('1','2','X'), end=c(1000000,2000000,10000000)) #plot SCNA heatmap for selected samples and genomic ranges
#' @export
plotCNV = function(segDf, 
                   genomeVersion = 'hg38', 
                   samples = NULL, 
                   chr = NULL,
                   start = NULL,
                   end = NULL){
  require(ggplot2, quietly = TRUE)
  require(scales, quietly = TRUE)
  if(!all(c("sample","chromosome","start","end","log2") %in% colnames(segDf))){
    stop("Error: the input segDf should contain 5 columns:\n
         \t\"sample\",\"chromosome\",\"start\",\"end\",\"log2\"")
  }
  
  #harmonize the chromosome to be the style of "Y"
  segDf$chromosome = as.character(segDf$chromosome)
  if(grepl(pattern = "(chr)|(Chr)",segDf$chromosome[1]))
    segDf$chromosome = substr(segDf$chromosome, 4, nchar(segDf$chromosome))
  segDf$chromosome[segDf$chromosome == "23"] = "X"
  segDf$chromosome[segDf$chromosome == "24"] = "Y"
  
  segDf = segDf[segDf$chromosome %in% c(1:22, 'X', 'Y'), ]
  
  #subset the samples as required, also specify the sample order
  if(!is.null(samples)){
    if(!all(samples %in% segDf$sample))
      stop('Samples to be shown should all be included in the segment-level data input!')
    segDf = segDf[segDf$sample %in% samples,]
  }
  
  #fill the gap as NAs, as required
  if(genomeVersion=='hg38')
    chrLengthObj = .seg2Ranges(hg38ChrLengthDf)
  else if(genomeVersion=='hg19')
    chrLengthObj = .seg2Ranges(hg19ChrLengthDf)
  else
    stop("Current only support genome version hg19 and hg38")
  segDf = .fillwithNA.2(segDf = segDf, totalChrRangeObj = chrLengthObj)
  
  #select chromosomes and range if necessary
  if(!is.null(chr) | !is.null(start) | !is.null(end))
    segDf = .selectRange.2(segDf, 
                        chr = chr, 
                        start = start, 
                        end = end, 
                        genomeVersion = genomeVersion)
  
  #specify the sample order for plotting
  if(!is.null(samples))
    segDf$sample=factor(segDf$sample, levels = samples)
  
  segDf$chromosome = factor(segDf$chromosome, levels = c(1:22, 'X', 'Y'))
  segDf$name = paste(segDf$sample, segDf$chromosome, segDf$start, sep = '-')
  segDf$width = segDf$end - segDf$start
  segDf$log2 = ifelse(segDf$log2 > 1, 1, segDf$log2)
  segDf$log2 = ifelse(segDf$log2 < -1, -1, segDf$log2)
  plotCNV = ggplot(segDf, aes(1, width, fill=log2)) + 
    geom_bar(stat="identity")+
    scale_fill_gradientn(values=c(-1, 0, 1), 
                         colours=c("blue", "white", "red"))+
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
#' @param segDf segment-level SCNA data frame containing only one sample. This is usually the same input as GISTIC2. Make sure there are these columns "sample", "chromosome", "start", "end", "log2".
#' @param totalChrRangeObj a GenomicRanges object, specifying the lengths of all chromosomes
#'
#' @return a data frame with gap-filled segment-level SCNA data
#'
.fillwithNA = function(segDf, totalChrRangeObj){
  segObj = GRanges(seqnames = segDf$chromosome, 
                   ranges = IRanges(start = segDf$start, end = segDf$end),
                   strand = Rle(strand("*"), nrow(segDf)),
                   log2 = segDf$log2)
  ##the segmenets not covered by the segObj
  naRanges = GenomicRanges::setdiff(totalChrRangeObj, segObj)
  if(length(naRanges)!=0){
    naRanges$log2 = NA
    #combining naRanges to segObj
    segObj = c(segObj, naRanges)
  }
  segObj = sort(segObj)
  segObj = as.data.frame(segObj)
  segObj$sample = segDf$sample[1]
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
    .fillwithNA(segDf, totalChrRangeObj)
  else{
    segList = split(segDf, f = segDf$sample)
    segList = lapply(segList, function(x).fillwithNA(segDf = x, totalChrRangeObj = totalChrRangeObj))
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
  rangeObj = GRanges(seqnames = totalChrLengthDf$chromosome, 
                     ranges = IRanges(start = 0, end = totalChrLengthDf$length),
                     strand = Rle(strand("*"), nrow(totalChrLengthDf)))
  rangeObj
}


#' Subset genomic ranges from a single sample's segment-level SCNA data
#'
#' @param segDf segment-level SCNA data frame containing only one sample. This is usually the same input as GISTIC2. Make sure there are these columns: "sample", "chromosome", "start", "end", "log2".
#' @param chr a character vector specifying which chromosome(s) to select
#' @param start a integer vector specifying the start positions of each chromosome.
#' If NULL, the position will be the begining of the chromosome(s). Default is NULL.
#' @param end a integer vector specifying the end positions of each chromosome. 
#' If NULL, the position will be the end of the chromosome(s). Default is NULL.
#' @param genomeVerison which genomeVersion are you working on? Can be either 'hg38' and 'hg19'. Default is 'hg38'. 
#'
#' @return a data frame representing the subsetted segment-level SCNA data from the input. 
#'
.selectRange = function(segDf, 
                       chr = NULL, 
                       start = NULL, 
                       end = NULL, 
                       genomeVerison = 'hg38'){
  
  if(genomeVerison == 'hg38')
    chrLength = hg38ChrLengthDf
  else if(genomeVerison == 'hg19')
    chrLength = hg19ChrLengthDf
  if(is.null(chr))
    chr = chrLength$chromosome
  if(is.null(start)){
    start = rep(0, length(chr))
  }
  if(is.null(end)){
    end = chrLength[chr,"length"]
  }
  if(any(start > chrLength[chr,"length"]) | any(start < 0) | any(end > chrLength[chr,"length"]))
    stop("Some displaying ranges exceed the chromosome length, please corret it.")
  
  #convert segDf to Grange Obj
  if(grepl(pattern = "chr",segDf$chromosome[1]))
    segDf$chromosome = substr(segDf$chromosome, 4, nchar(segDf$chromosome))
  segDf$chromosome[segDf$chromosome == "23"] = "X"
  segDf$chromosome[segDf$chromosome == "24"] = "Y"
  segObj = GRanges(seqnames = segDf$chromosome, 
                   ranges = IRanges(start = segDf$start, end = segDf$end),
                   strand = Rle(strand("*"), nrow(segDf)),
                   log2 = segDf$log2)
  
  findOverlapRegions = function(singleChr, start, end, queryRegion){
    refRangeObj = GRanges(seqnames = singleChr, 
                          ranges = IRanges(start = start, end = end),
                          strand = Rle(strand("*"), length(singleChr)))
    overlapIndex = findOverlaps(query = queryRegion, subject = refRangeObj, type = 'any')
    overlapIndex = as.data.frame(overlapIndex)
    overlapDf = as.data.frame(queryRegion[overlapIndex$queryHit])
    overlapDf[1,"start"] = start
    overlapDf[1, "width"] = overlapDf[1,"end"] -overlapDf[1,"start"]+1
    overlapDf[nrow(overlapDf), "end"] = end
    overlapDf[nrow(overlapDf), "width"] = overlapDf[nrow(overlapDf),"end"] -overlapDf[nrow(overlapDf),"start"]+1
    overlapDf
  }
  if(length(chr) == 1) subSegDf = findOverlapRegions(chr, start, end, segObj)
  else{
    subSegDf = lapply(1:length(chr), function(x) findOverlapRegions(chr[x], start[x], end[x], segObj))
    subSegDf = do.call(rbind, subSegDf)
  }
  colnames(subSegDf)[1] = "chromosome"
  subSegDf
}

#' Subset genomic ranges from multiple samples' segment-level SCNA data
#'
#' @param segDf segment-level SCNA data frame containing one or multiple samples. This is usually the same input as GISTIC2. 
#' Make sure there are these columns: "sample", "chromosome", "start", "end", "log2".
#' @param chr a character vector specifying which chromosome(s) to select
#' @param start a integer vector specifying the start positions of each chromosome.
#' If NULL, the position will be the begining of the chromosome(s). Default is NULL.
#' @param end a integer vector specifying the end positions of each chromosome. 
#' If NULL, the position will be the end of the chromosome(s). Default is NULL.
#' @param genomeVersion which genomeVersion are you working on? Can be either 'hg38' and 'hg19'. Default is 'hg38'. 
#'
#' @return a data frame representing the subsetted segment-level SCNA data from the input. 
#'
.selectRange.2 = function(segDf, 
                          chr = NULL, 
                          start = NULL, 
                          end = NULL, 
                          genomeVersion = 'hg38'){
  if(length(unique(segDf$sample)) == 1)
    .selectRange(segDf, chr, start, end, genomeVersion)
  else{
    segList = split(segDf, f = segDf$sample)
    segList = lapply(segList, function(x).selectRange(x, chr, start, end, genomeVersion))
    newSegDf = do.call(rbind, segList)
    newSegDf$sample = rep(names(segList), sapply(segList, nrow))
    newSegDf
  }
}

