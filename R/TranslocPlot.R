#!/usr/bin/Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfile", "character", "",
    "output","character", "file path to plot to"
  )
  
  OPTS <- c(
    "binsize","integer",2000000,"bps per bin",
    "assembly","character","mm9","genome assembly",
    "chrdisp","character","","comma-seperated list of included chromosomes",
    "start","integer",0,"start",
    "end","integer",0,"end",
    "stranddisp","integer",1,"1 for separate 0 for combined"
  )
  
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
} else {
  source("~/Pipelines/R/Rsub.R")
  tlxfile <- "/Volumes/AltLab/Translocation/RawData/Alt024-20130429/NewPipelineTest/results-full/CC004_Alt024/CC004_Alt024.tlx"
  output <- "~/Working/TranslocTesting/TranslocPlot.pdf"
  binsize <- 2000000
  assembly <- "mm9"
  chrdisp <- ""
  stranddisp <- 1
  start <- 0
  end <- 0
}

#suppressPackageStartupMessages()
library(Rsamtools)
library(GenomicRanges)
library(grid)
library(RColorBrewer)

denom <- c(1,5,20,100,500,2000,10000)
pal <- brewer.pal(9,"Set1")
pal <- pal[c(9,1,5,2,3,4,7)]

cyto <- read.delim(paste(Sys.getenv('GENOME_DB'),assembly,'annotation/cytoBand.txt',sep="/"),header=F,as.is=T,col.names=c('Chr','Start','End','ID','Type'))
cytocolor <- getCytoColor()
cyto$Color <- cytocolor[match(cyto$Type,names(cytocolor))]

chrlen <- read.delim(paste(Sys.getenv('GENOME_DB'),assembly,'annotation/ChromInfo.txt',sep="/"),header=F,as.is=T,col.names=c('Name','Length'))
chrdisp <- strsplit(chrdisp,",")[[1]]

chrlen <- structure(chrlen$Length,names=as.character(chrlen$Name))

chrnum <- as.numeric(sub("chr","",names(chrlen[grep("chr[0-9]+",names(chrlen),perl=T)])))
chrlet <- sub("chr","",names(chrlen[grep("chr[A-Z]+",names(chrlen),perl=T)]))
chrnum <- paste("chr",chrnum[order(chrnum)],sep="")
chrlet <- paste("chr",chrlet[match(c("X","Y","M"),chrlet)],sep="")
chrlevels <- names(chrlen)[c(match(chrnum,names(chrlen)),match(chrlet,names(chrlen)))]

chrlen <- chrlen[chrlevels]

if (length(chrdisp) == 0) {
  chrdisp <- seqlevels(gr)[seqlevels(gr) != "chrM"]
} else {
  chrdisp <- chrdisp[order(match(chrdisp,seqlevels(gr)))]
}

chrs <- rep(names(chrlen),c(ceiling(chrlen/binsize)))
strands <- rep(c("+","-"),each=length(chrs))
ends <- unlist(lapply(chrlen,function(len){seq(from=len,to=1,by=-binsize)}))
starts <- unlist(lapply(ends,function(end){ max(1,end-binsize+1) } ))

gr <- GRanges(seqnames=rep(chrs,2),ranges=IRanges(start=rep(starts,2),end=rep(ends,2)),strand=strands,seqlengths=chrlen)
seqlevels(gr) <- chrlevels

con  <- file(tlxfile, open = "r")
header <- unlist(strsplit(readLines(con, n = 1),"\t"))
close(con)

headersToRead <- c("Rname","Junction","Strand")

colClasses <- rep("NULL",length(header))
colClasses[match(headersToRead,header)] <- NA

tlx <- read.delim(tlxfile,header=T,colClasses=colClasses)
tlx <- GRanges(seqnames=tlx$Rname,ranges=IRanges(start=tlx$Junction,end=tlx$Junction),strand=ifelse(tlx$Strand==1,"+","-"))
if (stranddisp == 0) {
  strand(tlx) <- "+"
}

gr$hits <- countOverlaps(gr,tlx)

tmphits <- gr$hits
gr$hitvec <- list(c())

for (i in length(denom):1) {
  gr$hitvec <- mapply(function(x,y) { if (y > 0) c(x,rep(i,y)) else x },gr$hitvec,tmphits%/%denom[i])
  tmphits <- tmphits - tmphits%/%denom[i]*denom[i]
}

chrpos <- rep(0,length(chrdisp))
names(chrpos) <- chrdisp
gr$ypos <- 0

for (i in 1:length(chrdisp)) {
  
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


pdf(output,width=10,height=8)

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
    
    grid.polygon(x=unit(c(datapoints$stackpos-0.5,datapoints$stackpos+0.5,datapoints$stackpos),"native"),y=unit(c(datapoints$start+binsize,datapoints$start+binsize,datapoints$start),"native"),id=rep(1:nrow(datapoints),3),gp=gpar(fill=pal[datapoints$type],lty=0))
    
    popViewport()
  }
  
  chrposbins <- gr[seqnames(gr) == chrdisp[i] & strand(gr) == "+"]
  datapoints <- data.frame(type=unlist(chrposbins$hitvec),start=rep(start(chrposbins),unlist(lapply(chrposbins$hitvec,length))),end=rep(end(chrposbins),unlist(lapply(chrposbins$hitvec,length))),stackpos=unlist(lapply(lapply(chrposbins$hitvec,length),function(x) { if (x>0) seq(1,x) })))
  
  if (nrow(datapoints) > 0) {
    
    
    xmax <- max(unlist(lapply(chrposbins$hitvec,length)))
    
    dataheight <- convertY(unit(binsize,"native"),"mm")
    xwidth <- min(unit(xmax,"native"),xmax*dataheight)
    
    pushViewport(viewport(x=unit(chrpos[i],"native")+unit(0.1,"lines"),y=unit(0,"npc"),width=xwidth,height=unit(chrlen[names(chrlen) %in% chrdisp][i],"native"),just=c("left","bottom"),xscale=c(0.5,xmax+0.5),yscale=c(chrlen[names(chrlen) %in% chrdisp][i],0)))
    
    grid.polygon(x=unit(c(datapoints$stackpos-0.5,datapoints$stackpos+0.5,datapoints$stackpos),"native"),y=unit(c(datapoints$start,datapoints$start,datapoints$start+binsize),"native"),id=rep(1:nrow(datapoints),3),gp=gpar(fill=pal[datapoints$type],lty=0))
    
    popViewport()
  }
  
  
}

dev.off()
