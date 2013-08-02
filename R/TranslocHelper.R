createGenomicRanges <- function (chrlen,start,end,binsize) {
  if (start == 0) start <- 1
  if (end == 0) end <- chrlen
  chrs <- rep(names(chrlen),c(ceiling((end-start+1)/binsize)))
  strands <- rep(c("+","-"),each=length(chrs))
  ends <- unlist(lapply(end,function(x){rev(seq(from=x,to=start,by=-binsize))}))
  starts <- unlist(lapply(ends,function(end){ max(start,end-binsize+1) } ))
  gr <- GRanges(seqnames=rep(chrs,2),ranges=IRanges(start=rep(starts,2),end=rep(ends,2)),strand=strands,seqlengths=chrlen)
  seqlevels(gr) <- names(chrlen)
  return(gr)
} 

getChromLens <- function (assembly) {
  chrlen <- read.delim(paste(Sys.getenv('GENOME_DB'),assembly,'annotation/ChromInfo.txt',sep="/"),header=F,as.is=T,col.names=c('Name','Length'))
  # Convert data frame to vector with element names set as chromosome names
  chrlen <- structure(chrlen$Length,names=as.character(chrlen$Name))
  
  # Organize chromosome names into the correct order
  chrnum <- as.numeric(sub("chr","",names(chrlen[grep("chr[0-9]+",names(chrlen),perl=T)])))
  chrlet <- sub("chr","",names(chrlen[grep("chr[A-Z]+",names(chrlen),perl=T)]))
  chrnum <- paste("chr",chrnum[order(chrnum)],sep="")
  chrlet <- paste("chr",chrlet[match(c("X","Y","M"),chrlet)],sep="")
  chrlevels <- names(chrlen)[c(match(chrnum,names(chrlen)),match(chrlet,names(chrlen)))]
  # Order the chrlen object by this order
  chrlen <- chrlen[chrlevels]
  return(chrlen)
}

getCytoBands <- function (assembly) {
  cyto <- read.delim(paste(Sys.getenv('GENOME_DB'),assembly,'annotation/cytoBand.txt',sep="/"),header=F,as.is=T,col.names=c('Chr','Start','End','ID','Type'))
  cytocolor <- getCytoColor()
  # Create new color column in cytoband data
  cyto$Color <- cytocolor[match(cyto$Type,names(cytocolor))]
  return(cyto)
}

getCytoColor <- function() {
  cytocolor <- c()
  
  cytocolor["gpos100"]  <- rgb(0,0,0,max=255)
  cytocolor["gpos"]     <- rgb(0,0,0,max=255)
  cytocolor["gpos75"]   <- rgb(130,130,130,max=255)
  cytocolor["gpos66"]   <- rgb(160,160,160,max=255)
  cytocolor["gpos50"]   <- rgb(200,200,200,max=255)
  cytocolor["gpos33"]   <- rgb(210,210,210,max=255)
  cytocolor["gpos25"]   <- rgb(200,200,200,max=255)
  cytocolor["gvar"]     <- rgb(220,220,220,max=255)
  cytocolor["gneg"]     <- rgb(255,255,255,max=255)
  cytocolor["acen"]     <- rgb(217,47,39,max=255)
  cytocolor["stalk"]    <- rgb(100,127,164,max=255)
  
  return(cytocolor)
}