plotJunctions <- function (gr,binsize,strand=1,plottype="dot",plotshape="arrow",pal=NULL,rotateVP=0) {

  if (sum(gr$hits) < 1) return()
  
  if (plottype == "dot") {
    dots <- data.frame(value=unlist(gr$hitvec),mids=rep(gr$mids,gr$hitveclen),stackpos=unlist(lapply(gr$hitveclen,function(x) { if (x>0) seq(1,x) })))
    if (rotateVP) {
      dotwidth <- convertHeight(-1*unit(binsize,"native"),"mm",valueOnly=T)
      vpheight <- convertWidth(unit(1,"npc"),"mm",valueOnly=T)
      dotheight <- convertWidth(strand*unit(min(dotwidth,vpheight/max(gr$hitveclen)),"mm"),"native",valueOnly=T)
      
    } else {
      dotwidth <- convertWidth(unit(binsize,"native"),"mm",valueOnly=T)
      vpheight <- convertHeight(unit(1,"npc"),"mm",valueOnly=T)
      dotheight <- convertHeight(strand*unit(min(dotwidth,vpheight/max(gr$hitveclen)),"mm"),"native",valueOnly=T)
      
    }
    
    if (plotshape == "arrow") {
      xpoints <- unit(c(dots$mids-strand*binsize/2,dots$mids-strand*binsize/5,dots$mids-strand*binsize/5,dots$mids+strand*binsize/2,dots$mids-strand*binsize/5,dots$mids-strand*binsize/5,dots$mids-strand*binsize/2),"native")
      ypoints <- unit(c(dots$stackpos-1/3,dots$stackpos-1/3,dots$stackpos,dots$stackpos-1/2,dots$stackpos-1,dots$stackpos-2/3,dots$stackpos-2/3)*dotheight,"native")
      vertices <- 7
    } else if (plotshape == "triangle") {
    xpoints <- unit(c(dots$mids-strand*binsize/2,dots$mids-strand*binsize/2,dots$mids+strand*binsize/2),"native")
    ypoints <- unit(c(dots$stackpos-1,dots$stackpos,dots$stackpos-0.5)*dotheight,"native")
    vertices <- 3
    } else if (plotshape == "diamond") {
      xpoints <- unit(c(dots$mids,dots$mids-strand*binsize/2,dots$mids,dots$mids+strand*binsize/2),"native")
      ypoints <- unit(c(dots$stackpos-1,dots$stackpos-0.5,dots$stackpos,dots$stackpos-0.5)*dotheight,"native")
      vertices <- 4
    }
    if (rotateVP) {
      tmp <- xpoints
      xpoints <- ypoints
      ypoints <- tmp
    }
    
    grid.polygon(x=xpoints,y=ypoints,id=rep(1:nrow(dots),vertices),gp=gpar(fill=pal[dots$value],lty=0))
    
    
  } else {
    xpoints <- gr$mids
    if (plottype == "linear") {
      ypoints <- gr$hits
    } else {
      ypoints <- log10(gr$hits)
      ypoints <- ifelse(ypoints == -Inf,-1,ypoints)
    }
    if (rotateVP) {
      tmp <- xpoints
      xpoints <- ypoints
      ypoints <- tmp
    }
    if (strand == 1) {
      linecolor <- pal[4]
    } else {
      linecolor <- pal[2]
    }
    grid.lines(x=unit(xpoints,"native"),y=unit(ypoints,"native"),gp=gpar(col=linecolor,lwd=2))
  }
}

createGenomicRanges <- function (chrlen,start=0,end=0,mid=0,window=0,binsize=0,binnum=0) {
  if (length(chrlen) > 1) {
    ends <- unlist(lapply(chrlen,function(x){rev(seq(from=x,to=1,by=-binsize,))}))
    starts <- unlist(lapply(ends,function(end){ max(1,end-binsize+1) } ))
    chrs <- rep(names(chrlen),c(ceiling((chrlen)/binsize)))
    strands <- rep(c("+","-"),each=length(chrs))
  } else if (mid != 0 && window != 0) {
    if (binsize == 0) binsize <- ceiling(2*window/binnum)
    starts <- c(rev(seq(from=mid-binsize,to=mid-window,by=-binsize)),seq(from=mid,to=mid+window-binsize,by=binsize))
    ends <- starts + binsize - 1
    chrs <- rep(names(chrlen),length(starts))
    strands <- rep(c("+","-"),each=length(chrs))
  } else if (start != 0 && end != 0 ) {
    if (binsize == 0) binsize <- ceiling((end-start+1)/binnum)
    starts <- seq(from=start,to=end-binsize+1,by=binsize)
    ends <- starts + binsize - 1
    chrs <- rep(names(chrlen),length(starts))
    strands <- rep(c("+","-"),each=length(chrs))
  } else {
    if (binsize == 0) binsize <- ceiling((end-start+1)/binnum)
    ends <- unlist(lapply(chrlen,function(x){rev(seq(from=x,to=1,by=-binsize,))}))
    starts <- unlist(lapply(ends,function(end){ max(1,end-binsize+1) } ))
    chrs <- rep(names(chrlen),length(starts))
    strands <- rep(c("+","-"),each=length(chrs))
  }
  
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

formatBP <- function(x,pts=1) {
  unlist(lapply(x, function (bp) {
    if (bp > 100000000) paste(round(bp/1000000000,pts),"Gb") else
    if (bp > 100000) paste(round(bp/1000000,pts),"Mb") else
    if (bp > 100) paste(round(bp/1000,pts),"Kb") else
    paste(bp,"bp")
  }))
}
