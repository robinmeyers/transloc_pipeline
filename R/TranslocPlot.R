#!/usr/bin/Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfile", "character", "",
    "output","character", "file path to plot to"
  )
  
  OPTS <- c(
    "binsize","integer",2000000,"bps per bin",
    "assembly","character","mm9","genome assembly",
    "chr","character","","",
    "start","integer",0,"start",
    "end","integer",0,"end",
    "strand","integer",0,"1 for positive strand, -1 for minus strand, 0 for both strands, 2 for combined strands"
  )
  
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  source_local("TranslocHelper.R")
  parseArgs("TranslocPlot.R", ARGS, OPTS)
  
} else {
  source("~/Pipelines/R/Rsub.R")
  tlxfile <- "/Volumes/AltLab/Translocation/RawData/Alt024-20130429/NewPipelineTest/results-full/CC004_Alt024/CC004_Alt024.tlx"
  output <- "~/Working/TranslocTesting/TranslocPlot.pdf"
  binsize <- 2000000
  assembly <- "mm9"
  chr <- ""
  strand <- 0
  start <- 0
  end <- 0
  showM <- 0
}

suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(RColorBrewer))

denom <- c(1,5,20,100,500,2000,10000)
pal <- brewer.pal(9,"Set1")
pal <- pal[c(9,1,5,2,3,4,7)]

# Read in cytogenetic band data
cyto <- getCytoBands(assembly)
  
# Read in chomrosome length data
chrlen <- getChromLens(assembly)

if (chr != "") {
  if (! chr %in% names(chrlen)) stop("Error: chromosome ",chr," not found")
  chrlen <- chrlen[chr]
} else {
  if (!showM) chrlen <- chrlen[names(chrlen)!="chrM"]
}

gr <- createGenomicRanges(chrlen,start,end,binsize)

header <- readHeader(tlxfile)

columnsToRead <- c("Rname","Junction","Strand")

colClasses <- proColumnClasses(header,columnsToRead)

tlx <- read.delim(tlxfile,header=T,colClasses=colClasses)
tlxtot <- nrow(tlx)
  
tlx <- tlx[tlx$Rname %in% names(chrlen),]
if (strand == 1 || strand == -1) {
  tlx <- tlx[tlx$Strand == strand,]
} else if (strand == 2) {
  tlx$Strand <- 1
}
tlxgr <- GRanges(seqnames=tlx$Rname,ranges=IRanges(start=tlx$Junction,end=tlx$Junction),strand=ifelse(tlx$Strand==1,"+","-"),seqlengths=chrlen)

gr$hits <- countOverlaps(gr,tlxgr)

tlxdisp <- sum(gr$hits)





tmphits <- gr$hits
gr$hitvec <- list(c())

for (i in length(denom):1) {
  gr$hitvec <- mapply(function(x,y) { if (y > 0) c(x,rep(i,y)) else x },gr$hitvec,tmphits%/%denom[i])
  tmphits <- tmphits - tmphits%/%denom[i]*denom[i]
}

if (length(chrlen) > 1) {
  chrpos <- rep(0,length(chrlen))
  names(chrpos) <- names(chrlen)
  gr$ypos <- 0
  
  for (i in 1:length(chrlen)) {
    
    gr[seqnames(gr) == chrdisp[i]]$ypos <- match(start(gr[seqnames(gr) == chrdisp[i]]),unique(start(gr[seqnames(gr) == chrdisp[i]])))
    
    if (i == 1) {
      
      chrpos[i] <- max(unlist(lapply(gr[seqnames(gr) == chrdisp[i] & strand(gr) == "-"]$hitvec,length))) + 1
    } else {
      
      chrnegbins <- gr[seqnames(gr) == chrdisp[i] & strand(gr) == "-"]
      chrnegbins <- chrnegbins[order(chrnegbins$ypos)]
      
      chrposbins <- gr[seqnames(gr) == chrdisp[i-1] & strand(gr) == "+"]
      chrposbins <- chrposbins[order(chrposbins$ypos)]
      
      chrnegbins <- chrnegbins[chrnegbins$ypos <=  min(max(chrnegbins$ypos),max(chrposbins$ypos))]
      chrposbins <- chrposbins[chrposbins$ypos <=  min(max(chrnegbins$ypos),max(chrposbins$ypos))]
      
      chrpos[i] <- chrpos[i-1] + max(unlist(mapply(function(x,y,w,z) {if (w==z) {length(x) + length(y)} else 0},chrposbins$hitvec,chrnegbins$hitvec,chrposbins$ypos,chrnegbins$ypos))) + 2
      
    }
  }
  
  rightlim <- chrpos[length(chrpos)] + max(unlist(lapply(gr[seqnames(gr) == chrdisp[length(chrdisp)] & strand(gr) == "+"]$hitvec,length))) + 1
}


pdf(output,width=11,height=8.5)

pageVP <- viewport(name="page",width=unit(1,"npc")-unit(2*0.25,"inches"),height=unit(1,"npc")-unit(2*0.25,"inches"))
pushViewport(pageVP)
grid.rect()

titleVP <- viewport(y=unit(1,"npc"),name="title",width=unit(1,"npc"),height=unit(3,"lines"),just="top")
pushViewport(titleVP)

titletext <- paste(sub(".tlx","",basename(tlxfile))," - ",prettyNum(tlxtot,big.mark=",")," hits - ",formatBP(binsize,1)," bins",sep="")
grid.text(titletext,gp=gpar(cex=min(2,1/convertHeight(stringHeight(titletext),"npc",valueOnly=T),1/convertWidth(stringWidth(titletext),"npc",valueOnly=T))))

popViewport()
plotVP <- viewport(name="plot",y=0,height=unit(1,"npc")-unit(3,"lines"),just="bottom")
pushViewport(plotVP)
grid.rect()

genomeVP <- viewport(name="genome",y=unit(2,"lines"),height=unit(1,"npc")-unit(2,"lines"),just="bottom")
pushViewport(genomeVP)
grid.rect()

dev.off()

ymax <- max(chrlen[names(chrlen) %in% chrdisp])

pushViewport(viewport(name="PlotArea",y=unit(2,"lines"),height=unit(1,"npc")-unit(8,"lines"),width=unit(1,"npc")-unit(2,"lines"),just="bottom",xscale=c(0,rightlim),yscale=c(0,ymax)))


pushViewport(viewport(name="Title",y=unit(1,"npc"),height=unit(6,"lines"),just="bottom"))
grid.text(length(gr),y=unit(1,"lines"),gp=gpar(cex=1.5))
grid.text(sub(".tlx","",basename(tlxfile)),y=unit(2,"lines"),gp=gpar(cex=2))
popViewport()
pushViewport(viewport(name="Legend",x=unit(0.8,"npc"),y=unit(1,"npc")-unit(5,"lines"),width=unit(4,"lines"),height=unit(8,"lines"),just=c("right","bottom")))
grid.rect()
grid.text(label=denom,x=unit(0.9,"npc"),y=unit(1:length(denom),"lines"),just="right")
grid.polygon(x=unit(rep(c(0.5,0.5,1),each=length(denom)),"lines"),y=unit(c(1:length(denom)-0.25,1:length(denom)+0.25,1:length(denom)),"lines"),id=rep(1:length(denom),3),gp=gpar(fill=pal,lty=0))
popViewport()

dataheight <- convertY(unit(binsize,"native"),"mm")
xwidth <- min(dataheight*rightlim,convertX(unit(rightlim,"native"),"mm"))

pushViewport(viewport(name="Genome",width=xwidth,xscale=c(0,rightlim),yscale=c(0,ymax)))
grid.rect(x=unit(chrpos,"native"),y=unit(0,"npc"),height=unit(chrlen[names(chrlen) %in% chrdisp],"native"),width=unit(0.2,"lines"),just=c("center","bottom"))

cyto <- cyto[cyto$Chr %in% names(chrpos),]
cyto$Xpos <- chrpos[match(cyto$Chr,names(chrpos))]
grid.rect(x=unit(cyto$Xpos,"native"),y=unit(chrlen[cyto$Chr] - cyto$Start,"native"),width=unit(0.2,"lines"),height=unit(cyto$End-cyto$Start,"native"),just="top",gp=gpar(fill=cyto$Color,lty=0))
pushViewport(viewport(name="ChrLabels",y=unit(0,"npc"),height=unit(1,"lines"),just="top",xscale=c(0,rightlim)))
grid.text(sub("chr","",names(chrpos)),x=unit(chrpos,"native"),y=0.5,just="center")
popViewport()

for (i in 1:length(chrdisp)) {
  
  chrnegbins <- gr[seqnames(gr) == chrdisp[i] & strand(gr) == "-"]
  datapoints <- data.frame(type=unlist(chrnegbins$hitvec),start=rep(start(chrnegbins),unlist(lapply(chrnegbins$hitvec,length))),end=rep(end(chrnegbins),unlist(lapply(chrnegbins$hitvec,length))),stackpos=unlist(lapply(lapply(chrnegbins$hitvec,length),function(x) { if (x>0) seq(1,x) })))
  
  if (nrow(datapoints) > 0) {  
    
    xmax <- max(unlist(lapply(chrnegbins$hitvec,length)))
    
    dataheight <- convertY(unit(binsize,"native"),"mm")
    xwidth <- min(unit(xmax,"native"),xmax*dataheight)
    
    pushViewport(viewport(x=unit(chrpos[i],"native")-unit(0.1,"lines"),y=unit(0,"npc"),width=xwidth,height=unit(chrlen[names(chrlen) %in% chrdisp][i],"native"),just=c("right","bottom"),xscale=c(xmax+0.5,0.5),yscale=c(chrlen[names(chrlen) %in% chrdisp][i],0)))
    
    grid.polygon(x=unit(c(datapoints$stackpos-0.5,datapoints$stackpos+0.5,datapoints$stackpos),"native"),y=unit(c(datapoints$start+binsize/2,datapoints$start+binsize/2,datapoints$start-binsize/2),"native"),id=rep(1:nrow(datapoints),3),gp=gpar(fill=pal[datapoints$type],lty=0))
    
    popViewport()
  }
  
  chrposbins <- gr[seqnames(gr) == chrdisp[i] & strand(gr) == "+"]
  datapoints <- data.frame(type=unlist(chrposbins$hitvec),start=rep(start(chrposbins),unlist(lapply(chrposbins$hitvec,length))),end=rep(end(chrposbins),unlist(lapply(chrposbins$hitvec,length))),stackpos=unlist(lapply(lapply(chrposbins$hitvec,length),function(x) { if (x>0) seq(1,x) })))
  
  if (nrow(datapoints) > 0) {
    
    
    xmax <- max(unlist(lapply(chrposbins$hitvec,length)))
    
    dataheight <- convertY(unit(binsize,"native"),"mm")
    xwidth <- min(unit(xmax,"native"),xmax*dataheight)
    
    pushViewport(viewport(x=unit(chrpos[i],"native")+unit(0.1,"lines"),y=unit(0,"npc"),width=xwidth,height=unit(chrlen[names(chrlen) %in% chrdisp][i],"native"),just=c("left","bottom"),xscale=c(0.5,xmax+0.5),yscale=c(chrlen[names(chrlen) %in% chrdisp][i],0)))
    
    grid.polygon(x=unit(c(datapoints$stackpos-0.5,datapoints$stackpos+0.5,datapoints$stackpos),"native"),y=unit(c(datapoints$start+binsize/2,datapoints$start+binsize/2,datapoints$start+1.5*binsize),"native"),id=rep(1:nrow(datapoints),3),gp=gpar(fill=pal[datapoints$type],lty=0))
    
    popViewport()
  }
  
  
}

dev.off()
