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

plotXScale <- function(chr,rstart,rend) {
  
  sizeArray <- c(50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000,200000000)
  grid.xaxis(at=c(rstart,rend),label=c(" "," "))
  grid.text(label=formatBP(rend-rstart))
  grid.text(label=prettyNum(rstart,big.mark=","),x=unit(0,"npc"),just="left")
  grid.text(label=prettyNum(rend,big.mark=","),x=unit(1,"npc"),just="right")
}


plotFeatures <- function(features,chr,rstart,rend) {
  grid.text(features$Name,x=unit((features$Start+features$End)/2,"native"),y=unit(1-(((0:(nrow(features)-1))%%3))*0.33,"npc"),just="top",gp=gpar(cex=0.75))
#   grid.text(features$Name,x=unit((features$Start+features$End)/2,"native"),y=unit(0,"npc"),just="bottom",gp=gpar(cex=0.75))
  
}

printHeader <- function(tlxfile,tlxdisp,tlxtot,assembly,chr,denom,pal,plottype,plotshape) {
  
  if (plottype == "dot") {
    pushViewport(viewport(x=unit(1,"npc"),width=unit(1.5,"inches"),xscale=c(0,2),just="right"))
    x_legend <- unit(floor((1:length(denom)-1)/3) + 0.5,"native")
    y_legend <- unit((1:length(denom)-1)%%3 + 0.5,"lines")
    grid.text(label=denom,x=unit(x_legend,"native"),y=unit(y_legend,"lines"),just="left")    
    
    if (plotshape == "arrow") {
      vertices <- 7
      xpoints <- rep(x_legend-unit(1,"mm"),each=vertices) - unit(c(4,8/3,8/3,0,8/3,8/3,4),"mm")
      ypoints <- rep(y_legend,each=vertices) + unit(c(4/5,4/5,2,0,-2,-4/5,-4/5),"mm")
    } else if (plotshape == "triangle") {
      vertices <- 3
      xpoints <- rep(x_legend-unit(1,"mm"),each=vertices) - unit(c(4,0,4),"mm")
      ypoints <- rep(y_legend,each=vertices) + unit(c(2,0,-2),"mm")
    } else if (plotshape == "diamond") {
      vertices <- 4      
      xpoints <- rep(x_legend-unit(1,"mm"),each=vertices) - unit(c(4,2,0,2),"mm")
      ypoints <- rep(y_legend,each=vertices) + unit(c(0,2,0,-2),"mm")
    }
    
    grid.polygon(x=xpoints,y=ypoints,id=rep(1:length(denom),each=vertices),gp=gpar(fill=pal,lty=0))
    popViewport()
    
    textwidth <- unit(1,"npc")-unit(1.5,"inches")
    
  } else {
    textwidth <- unit(1,"npc")    
  }
  
  spacer <- unit(3,"mm")
  
  pushViewport(viewport(x=unit(0.5,"npc"),width=unit(0.5,"npc"),just="right"))
  titletext <- sub(".tlx","",basename(tlxfile))
  grid.text(titletext,x=unit(1,"npc")-spacer,just="right",gp=gpar(cex=min(2.5,1/convertHeight(stringHeight(titletext),"npc",valueOnly=T),(1-convertWidth(spacer,"npc",valueOnly=T))/convertWidth(stringWidth(titletext),"npc",valueOnly=T))))
  popViewport()
  
  pushViewport(viewport(x=unit(0.5,"npc"),width=textwidth,y=unit(1,"npc"),height=unit(0.5,"npc"),just=c("left","top")))
  hittext <- paste("Displaying",prettyNum(tlxdisp,big.mark=","),"of",prettyNum(tlxtot,big.mark=","),"hits")
  grid.text(hittext,x=spacer,just="left",gp=gpar(cex=min(1.25,1/convertHeight(stringHeight(hittext),"npc",valueOnly=T),(1-convertWidth(spacer,"npc",valueOnly=T))/convertWidth(stringWidth(hittext),"npc",valueOnly=T))))
  popViewport()
  
  pushViewport(viewport(x=unit(0.5,"npc"),width=unit(0.5,"npc"),y=unit(0,"npc"),height=unit(0.5,"npc"),just=c("left","bottom")))
  
  if (chr != "") {
    displaytext <- paste(assembly,"-",chr,"-",formatBP(binsize,1),"bins")
  } else {
    displaytext <- paste(assembly,"-",formatBP(binsize,1),"bins")
  }
  grid.text(displaytext,x=spacer,just="left",gp=gpar(cex=min(1.25,1/convertHeight(stringHeight(displaytext),"npc",valueOnly=T),(1-convertWidth(spacer,"npc",valueOnly=T))/convertWidth(stringWidth(displaytext),"npc",valueOnly=T))))
  popViewport()
  
  
  
  
#   
#   
#   
#   bintext <- paste(formatBP(binsize,1),"bins")
#   
#   
#   
#   titletext <- paste(sub(".tlx","",basename(tlxfile))," - Displaying ",prettyNum(tlxdisp,big.mark=",")," of ",prettyNum(tlxtot,big.mark=",")," hits - ",formatBP(binsize,1)," bins",sep="")
#   grid.text(titletext,gp=gpar(cex=min(2,1/convertHeight(stringHeight(titletext),"npc",valueOnly=T),1/convertWidth(stringWidth(titletext),"npc",valueOnly=T))))
}

createGenomicRanges <- function (chrlen,rstart=0,rend=0,rmid=0,rwindow=0,binsize=0,binnum=0) {
  if (length(chrlen) > 1) {
    rends <- unlist(lapply(chrlen,function(x){rev(seq(from=x,to=1,by=-binsize,))}))
    rstarts <- unlist(lapply(rends,function(rend){ max(1,rend-binsize+1) } ))
    chrs <- rep(names(chrlen),c(ceiling((chrlen)/binsize)))
    strands <- rep(c("+","-"),each=length(chrs))
  } else if (rmid != 0 && rwindow != 0) {
    if (binnum != 0) binsize <- ceiling(2*rwindow/binnum)
    rstarts <- c(rev(seq(from=rmid-binsize,to=rmid-rwindow,by=-binsize)),seq(from=rmid,to=rmid+rwindow-binsize,by=binsize))
    rends <- rstarts + binsize - 1
    chrs <- rep(names(chrlen),length(rstarts))
    strands <- rep(c("+","-"),each=length(chrs))
  } else if (rstart != 0 && rend != 0 ) {
    if (binnum != 0) binsize <- ceiling((rend-rstart+1)/binnum)
    rstarts <- seq(from=rstart,to=rend-binsize+1,by=binsize)
    rends <- rstarts + binsize - 1
    chrs <- rep(names(chrlen),length(rstarts))
    strands <- rep(c("+","-"),each=length(chrs))
  } else {
    rstart <- 1
    rend <- chrlen
    if (binnum != 0) binsize <- ceiling((rend-rstart+1)/binnum)
    rends <- rev(seq(from=rend,to=rstart,by=-binsize))
    rstarts <- unlist(lapply(rends,function(rend){max(1,rend-binsize+1)}))
    chrs <- rep(names(chrlen),length(rstarts))
    strands <- rep(c("+","-"),each=length(chrs))
  }
    
  gr <- GRanges(seqnames=rep(chrs,2),ranges=IRanges(start=rep(rstarts,2),end=rep(rends,2)),strand=strands,seqlengths=chrlen)
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

getFeatures <- function (assembly,featurefile) {
  if (featurefile == "") featurefile <- paste(Sys.getenv('GENOME_DB'),assembly,'annotation/refGene.bed',sep="/")
  features <- read.delim(featurefile,header=F,as.is=T,col.names=c('Chr','Start','End','Name'))
  return(features)
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
    if (bp > 1000000000) paste(round(bp/1000000000,pts),"Gb") else
    if (bp > 1000000) paste(round(bp/1000000,pts),"Mb") else
    if (bp > 1000) paste(round(bp/1000,pts),"Kb") else
    paste(bp,"bp")
  }))
}
