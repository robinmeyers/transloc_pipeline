#!/usr/bin/Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfile", "character", "",
    "output","character", "file path to plot to"
  )
  
  OPTS <- c(
    "binsize","integer",2500000,"bps per bin",
    "assembly","character","mm9","genome assembly",
    "strand","integer",0,"1 for positive strand, -1 for minus strand, 0 for both strands, 2 for combined strands",    
    "brkchr","character","","",
    "brksite","integer",0,"",
    "brkstrand","integer",0,"",
    "featurefile","character","","e.g. RefGene",
    "chr","character","","",
    "rstart","integer",0,"start",
    "rend","integer",0,"end",
    "rmid","integer",0,"",
    "rwindow","integer",0,"",
    "binnum","integer",0,"",
    "plottype","character","dot","",
    "plotshape","character","arrow","",
    "showM","integer",0,"",
    "showY","integer",0,"",
    "ymax","integer",0,""
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
  source("~/TranslocPipeline/R/Rsub.R")
  source("~/TranslocPipeline/R/TranslocHelper.R")
  tlxfile <- "/Volumes/AltLab/Translocation/RawData/Alt024-20130429/NewPipelineTest/results-full//CC004_Alt024/CC004_Alt024.tlx"
  output <- "~/Working/TranslocTesting/TranslocPlot.pdf"
  binsize <- 2500000
  assembly <- "mm9"
  featurefile <- ""
  chr <- "chr15"
  strand <- 0
  rstart <- 0
  rend <- 0
  rmid <- 61818880
  rwindow <- 10000
  binnum <- 100
  showM <- 0
  showY <- 0
  plottype <- "linear"
  plotshape <- "diamond"
  ymax <- 0
  brkchr <- "chr15"
  brksite <- 61818880
  brkstrand <- 1
}

#argument checking
if (chr == "") plottype <- "dot"
if (strand == 2) plotshape <- "diamond"


suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(RColorBrewer))

# denom <- c(1,5,20,100,500,2000,10000)
# pal <- brewer.pal(9,"Set1")
# pal <- c("black",pal[c(1,5,2,3,4,7)])

denom <- c(1,5,20,100,500,2000)
pal <- brewer.pal(9,"Set1")
pal <- c("black",pal[c(1,5,2,3,4)])

# Read in cytogenetic band data
cyto <- getCytoBands(assembly)
features <- getFeatures(assembly,featurefile)

# Read in chomrosome length data
chrlen <- getChromLens(assembly)

if (chr != "") {
  if (! chr %in% names(chrlen)) stop("Error: chromosome ",chr," not found")
  chrlen <- chrlen[chr]
} else {
  if (!showM) chrlen <- chrlen[names(chrlen)!="chrM"]
  if (!showY) chrlen <- chrlen[names(chrlen)!="chrY"]
}

gr <- createGenomicRanges(chrlen,rstart=rstart,rend=rend,rmid=rmid,rwindow=rwindow,binsize=binsize,binnum=binnum)
binsize <- end(gr)[length(gr)] - start(gr)[length(gr)] + 1

header <- readHeader(tlxfile)

columnsToRead <- c("Qname","Rname","Junction","Strand","B_Rend","B_Qend","Qstart")

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

gr$mids <- end(gr) - binsize/2

tlxdisp <- sum(gr$hits)

tmphits <- gr$hits
gr$hitvec <- list(c())
for (i in length(denom):1) {
  gr$hitvec <- mapply(function(x,y) { if (y > 0) c(x,rep(i,y)) else x },gr$hitvec,tmphits%/%denom[i])
  tmphits <- tmphits - tmphits%/%denom[i]*denom[i]
}
gr$hitveclen <- unlist(lapply(gr$hitvec,length))


pdf(output,width=11,height=8.5)

marginsize <- 0.5
marginunits <- "inches"

pageVP <- viewport(name="page",width=unit(1,"npc")-unit(2*marginsize,marginunits),height=unit(1,"npc")-unit(2*marginsize,marginunits))
pushViewport(pageVP)

headerheight <- 3
headerunits <- "lines"

headerVP <- viewport(y=unit(1,"npc"),name="header",width=unit(1,"npc"),height=unit(headerheight,headerunits),just="top")
pushViewport(headerVP)

printHeader(tlxfile,tlxdisp,tlxtot,assembly,chr,denom,pal,plottype,plotshape)

popViewport()
plotVP <- viewport(name="plot",y=0,height=unit(1,"npc")-unit(headerheight,headerunits),just="bottom")
pushViewport(plotVP)


if (length(chrlen) > 1) {

  rotateVP <- 1
  
  chrwidth <- 2
  chrwidthunit <- "mm"
  
  chrpos <- rep(0,length(chrlen))
  names(chrpos) <- names(chrlen)
  gr$ypos <- 0
  for (i in 1:length(chrlen)) {
    gr[seqnames(gr)==names(chrlen)[i] & strand(gr)=="-"]$ypos <- rev(1:length(gr[seqnames(gr)==names(chrlen)[i] & strand(gr)=="-"]))
    gr[seqnames(gr)==names(chrlen)[i] & strand(gr)=="+"]$ypos <- rev(1:length(gr[seqnames(gr)==names(chrlen)[i] & strand(gr)=="+"]))
    
    
    if (i == 1) {
      chrpos[i] <- max(gr[seqnames(gr)==names(chrlen)[i] & strand(gr)=="-"]$hitveclen)
      
    } else {
      
      negbins <- rev(gr[seqnames(gr) == names(chrlen)[i] & strand(gr) == "-"])
      
      posbins <- rev(gr[seqnames(gr) == names(chrlen)[i-1] & strand(gr) == "+"])
      
      bothbins <- min(length(posbins),length(negbins))
      
#       negbins <- negbins[negbins$ypos <=  min(max(negbins$ypos),max(posbins$ypos))]
#       posbins <- posbins[posbins$ypos <=  min(max(negbins$ypos),max(posbins$ypos))]
      
      chrpos[i] <- chrpos[i-1] + max(negbins$hitveclen[1:bothbins] + posbins$hitveclen[1:bothbins],
                                     negbins$hitveclen[1:min(bothbins,length(posbins)-1)]+posbins$hitveclen[2:min(bothbins+1,length(posbins))],
                                     negbins$hitveclen[2:min(bothbins+1,length(negbins))]+posbins$hitveclen[1:min(bothbins,length(negbins)-1)]) + 1
    }
  }
  chrposlim <- chrpos[length(chrpos)] + max(gr[seqnames(gr) == names(chrlen)[length(chrlen)] & strand(gr) == "+"]$hitveclen)
  
  totalwidthmm <- convertX(unit(1,"npc"),"mm",valueOnly=T)
  chrwidthnative <- chrwidth*chrposlim/(totalwidthmm-length(chrlen)*chrwidth)
  chrposlim <- chrposlim+chrwidthnative*length(chrlen)
  
  genomeVP <- viewport(name="genome",y=unit(2,"lines"),height=unit(1,"npc")-unit(2,"lines"),just="bottom",xscale=c(0,chrposlim),yscale=c(1,max(chrlen)))
  pushViewport(genomeVP)
  
  chrpos <- 0:(length(chrlen)-1)*chrwidthnative + chrpos
  
  grid.rect(x=unit(chrpos,"native"),y=unit(0,"npc"),width=unit(chrwidth,chrwidthunit),height=unit(chrlen,"native"),just=c("left","bottom"))
  cyto$Xpos <-chrpos[match(cyto$Chr,names(chrpos))]
  grid.rect(x=unit(cyto$Xpos,"native"),y=unit(chrlen[cyto$Chr] - cyto$Start,"native"),width=unit(chrwidth,chrwidthunit),height=unit(cyto$End-cyto$Start,"native"),just=c("left","top"),gp=gpar(fill=cyto$Color,lty=0))
  
  
  if (brkchr %in% names(chrpos) && brksite > 0) {
    arrowlen <- 5
    arrowlenunit <- "mm"
    arrowlennative <- convertHeight(unit(arrowlen,arrowlenunit),"native",valueOnly=T)
    
    
    grid.rect(x=unit(chrpos[brkchr],"native"),y=unit(chrlen[brkchr]-brksite,"native"),width=unit(chrwidth,chrwidthunit),just="left",height=unit(1,"mm"),gp=gpar(fill="yellow",linejoin="mitre"))
    if (brkstrand == 1 || brkstrand == -1) {
      ypoints <- unit(chrlen[brkchr]-c(brksite-arrowlennative,brksite-arrowlennative/2,brksite-arrowlennative/2,brksite,brksite-arrowlennative/2,brksite-arrowlennative/2,brksite-arrowlennative),"native")
      xpoints <- unit(c(chrpos[brkchr]+chrwidthnative*3/4,chrpos[brkchr]+chrwidthnative*3/4,chrpos[brkchr]+chrwidthnative,chrpos[brkchr]+chrwidthnative/2,chrpos[brkchr],chrpos[brkchr]+chrwidthnative*1/4,chrpos[brkchr]+chrwidthnative*1/4),"native") 
      grid.polygon(x=xpoints,y=ypoints,gp=gpar(fill="yellow",linejoin="mitre"))
    }
  }
  
  chrVPs <- list()
  
  for (i in 1:length(chrlen)) {
    negVP <- viewport(x=unit(chrpos[i],"native"), y=unit(0,"npc"), 
                      width=unit(max(gr[seqnames(gr) == names(chrlen)[i] & strand(gr) == "-"]$hitveclen),"native"),
                      height=unit(chrlen[i],"native"), just=c("right","bottom"),
                      xscale=c(1,0),yscale=c(chrlen[i],1),clip="on")
    posVP <- viewport(x=unit(chrpos[i]+chrwidthnative,"native"), y=unit(0,"npc"), 
                      width=unit(max(gr[seqnames(gr) == names(chrlen)[i] & strand(gr) == "+"]$hitveclen),"native"),
                      height=unit(chrlen[i],"native"), just=c("left","bottom"),
                      xscale=c(0,1),yscale=c(chrlen[i],1),clip="on")
    
    
    chrVPs[[i]] <- list(posVP,negVP)
  }
  pushViewport(viewport(name="chrlabel",y=unit(0,"npc"),height=unit(2,"lines"),just="top",xscale=c(0,chrposlim)))
  grid.text(sub("chr","",names(chrlen)),x=unit(chrpos+0.75/2,"native"),gp=gpar(cex=1.25))
  popViewport()
  
} else {
  
  rotateVP <- 0
  
  rstart <- min(start(gr))
  rend <- max(end(gr))
  
  if (plottype == "dot") {
    yscalewidth = unit(0,"lines")
  } else {
    yscalewidth = unit(2,"lines")
  }
  
  xscaleVP <- viewport(name="xscale",x=yscalewidth,width=unit(1,"npc")-yscalewidth,y=unit(1,"npc"),height=unit(2,"lines"),just=c("left","top"),xscale=c(rstart,rend),clip="off")
  pushViewport(xscaleVP)
  plotXScale(chr,rstart,rend)
  popViewport()
  
 
  
  genomeVP <- viewport(name="genome",x=yscalewidth,width=unit(1,"npc")-yscalewidth,height=unit(1,"npc")-unit(4,"lines"),just=c("left"),xscale=c(rstart,rend),yscale=c(-1,1))
  pushViewport(genomeVP)
  
  chrwidth <- 3
  chrwidthunit <- "mm"
  
  grid.rect(y=unit(0,"native"),height=unit(chrwidth,chrwidthunit))  
    
  if (rend - rstart < 2000000) {
    features <- subset(features,Chr == chr & End >= rstart & Start <= rend)
    features <- features[!duplicated(features$Name),]
    features <- features[with(features,order(Start)),]
    features$Start <- ifelse(features$Start < rstart, rstart, features$Start)
    features$End <- ifelse(features$End > rend, rend, features$End)
    grid.rect(x=unit(features$Start,"native"),y=unit(0,"native"),width=unit(features$End-features$Start,"native"),height=unit(chrwidth,chrwidthunit),just="left",gp=gpar(fill=getCytoColor()["gpos25"]))
    featureVP <- viewport(name="feature",y=unit(0,"npc"),height=unit(2,"lines"),just="top",xscale=c(rstart,rend),clip="off")
    pushViewport(featureVP)
    plotFeatures(features,chr,rstart,rend)
    popViewport()
      
  } else {
    cyto <- subset(cyto, Chr == chr & End >= rstart & Start <= rend)
    cyto$Start <- ifelse(cyto$Start < rstart, rstart, cyto$Start)
    cyto$End <- ifelse(cyto$End > rend, rend, cyto$End)
    grid.rect(x=unit(cyto$Start,"native"),y=unit(0,"native"),width=unit(cyto$End-cyto$Start,"native"),height=unit(chrwidth,chrwidthunit),just="left",gp=gpar(fill=cyto$Color,lty=0))
  }
  
  
  
  if (names(chrlen)[1] == brkchr && brksite >= rstart && brksite <= rend) {
    grid.rect(x=unit(brksite,"native"),y=unit(0,"native"),width=unit(1,"mm"),height=unit(chrwidth,chrwidthunit),gp=gpar(fill="yellow",linejoin="mitre"))
    
    if (brkstrand == 1 || brkstrand == -1) {
      arrowlen <- 5
      arrowlenunit <- "mm"
      arrowlennative <- convertWidth(unit(arrowlen,arrowlenunit),"native",valueOnly=T)
      
      xpoints <- unit(c(brksite-arrowlennative,brksite-arrowlennative/2,brksite-arrowlennative/2,brksite,brksite-arrowlennative/2,brksite-arrowlennative/2,brksite-arrowlennative),"native")
      ypoints <- unit(0,"native")+unit(c(chrwidth/4,chrwidth/4,chrwidth/2,0,-chrwidth/2,-chrwidth/4,-chrwidth/4),chrwidthunit)
      grid.polygon(x=xpoints,y=ypoints,gp=gpar(fill="yellow",linejoin="mitre"))
    }
  }
  
  ymin <- 0
  if (plottype == "dot") {
    ymax <- 1
  } else {
    if (ymax == 0) ymax <- max(gr$hits)
    if (plottype == "log") {
      ymax <- log10(ymax)
      ymin <- -1
    }
  }
    
  posVP <- viewport(y=unit(0,"native")+unit(chrwidth/2,chrwidthunit),height=unit(0.45,"npc")-unit(chrwidth/2,chrwidthunit),just="bottom",xscale=c(rstart,rend),yscale=c(ymin,ymax),clip="on")
  negVP <- viewport(y=unit(0,"native")-unit(chrwidth/2,chrwidthunit),height=unit(0.45,"npc")-unit(chrwidth/2,chrwidthunit),just="top",xscale=c(rstart,rend),yscale=c(ymax,ymin),clip="on")
  chrVPs <- list(list(posVP,negVP))
  
  if (plottype != "dot") {
    pushViewport(viewport(x=unit(-2,"mm"),y=unit(0,"native")+unit(chrwidth/2,chrwidthunit),height=unit(0.45,"npc")-unit(chrwidth/2,chrwidthunit),just=c("left","bottom"),yscale=c(ymin,ymax),clip="off"))
    grid.yaxis(gp=gpar(cex=0.75))
    popViewport()
    
    pushViewport(viewport(x=unit(-2,"mm"),y=unit(0,"native")-unit(chrwidth/2,chrwidthunit),height=unit(0.45,"npc")-unit(chrwidth/2,chrwidthunit),just=c("left","top"),yscale=c(ymax,ymin),clip="off"))
    grid.yaxis(gp=gpar(cex=0.75))
    popViewport()
  }
  
}

for (i in 1:length(chrlen)) {
  posVP <- chrVPs[[i]][[1]]
  negVP <- chrVPs[[i]][[2]]
  
  
  
  pushViewport(posVP)
  plotJunctions(gr[seqnames(gr) == names(chrlen)[i] & strand(gr) == "+",],binsize,strand=1,plottype=plottype,plotshape=plotshape,pal=pal,rotateVP=rotateVP)
  popViewport()
  
  pushViewport(negVP)
  plotJunctions(gr[seqnames(gr) == names(chrlen)[i] & strand(gr) == "-",],binsize,strand=-1,plottype=plottype,plotshape=plotshape,pal=pal,rotateVP=rotateVP)
  popViewport()
  
}
  


dev.off()
# 
# ymax <- max(chrlen[names(chrlen) %in% chrdisp])
# 
# pushViewport(viewport(name="PlotArea",y=unit(2,"lines"),height=unit(1,"npc")-unit(8,"lines"),width=unit(1,"npc")-unit(2,"lines"),just="bottom",xscale=c(0,rightlim),yscale=c(0,ymax)))
# 
# 
# pushViewport(viewport(name="Title",y=unit(1,"npc"),height=unit(6,"lines"),just="bottom"))
# grid.text(length(gr),y=unit(1,"lines"),gp=gpar(cex=1.5))
# grid.text(sub(".tlx","",basename(tlxfile)),y=unit(2,"lines"),gp=gpar(cex=2))
# popViewport()
# pushViewport(viewport(name="Legend",x=unit(0.8,"npc"),y=unit(1,"npc")-unit(5,"lines"),width=unit(4,"lines"),height=unit(8,"lines"),just=c("right","bottom")))
# grid.rect()
# grid.text(label=denom,x=unit(0.9,"npc"),y=unit(1:length(denom),"lines"),just="right")
# grid.polygon(x=unit(rep(c(0.5,0.5,1),each=length(denom)),"lines"),y=unit(c(1:length(denom)-0.25,1:length(denom)+0.25,1:length(denom)),"lines"),id=rep(1:length(denom),3),gp=gpar(fill=pal,lty=0))
# popViewport()
# 
# dataheight <- convertY(unit(binsize,"native"),"mm")
# xwidth <- min(dataheight*rightlim,convertX(unit(rightlim,"native"),"mm"))
# 
# pushViewport(viewport(name="Genome",width=xwidth,xscale=c(0,rightlim),yscale=c(0,ymax)))
# grid.rect(x=unit(chrpos,"native"),y=unit(0,"npc"),height=unit(chrlen[names(chrlen) %in% chrdisp],"native"),width=unit(0.2,"lines"),just=c("center","bottom"))
# 
# cyto <- cyto[cyto$Chr %in% names(chrpos),]
# cyto$Xpos <- chrpos[match(cyto$Chr,names(chrpos))]
# grid.rect(x=unit(cyto$Xpos,"native"),y=unit(chrlen[cyto$Chr] - cyto$Start,"native"),width=unit(0.2,"lines"),height=unit(cyto$End-cyto$Start,"native"),just="top",gp=gpar(fill=cyto$Color,lty=0))
# pushViewport(viewport(name="ChrLabels",y=unit(0,"npc"),height=unit(1,"lines"),just="top",xscale=c(0,rightlim)))
# grid.text(sub("chr","",names(chrpos)),x=unit(chrpos,"native"),y=0.5,just="center")
# popViewport()
# 
# for (i in 1:length(chrdisp)) {
#   
#   chrnegbins <- gr[seqnames(gr) == chrdisp[i] & strand(gr) == "-"]
#   datapoints <- data.frame(type=unlist(chrnegbins$hitvec),start=rep(start(chrnegbins),unlist(lapply(chrnegbins$hitvec,length))),end=rep(end(chrnegbins),unlist(lapply(chrnegbins$hitvec,length))),stackpos=unlist(lapply(lapply(chrnegbins$hitvec,length),function(x) { if (x>0) seq(1,x) })))
#   
#   if (nrow(datapoints) > 0) {  
#     
#     xmax <- max(unlist(lapply(chrnegbins$hitvec,length)))
#     
#     dataheight <- convertY(unit(binsize,"native"),"mm")
#     xwidth <- min(unit(xmax,"native"),xmax*dataheight)
#     
#     pushViewport(viewport(x=unit(chrpos[i],"native")-unit(0.1,"lines"),y=unit(0,"npc"),width=xwidth,height=unit(chrlen[names(chrlen) %in% chrdisp][i],"native"),just=c("right","bottom"),xscale=c(xmax+0.5,0.5),yscale=c(chrlen[names(chrlen) %in% chrdisp][i],0)))
#     
#     grid.polygon(x=unit(c(datapoints$stackpos-0.5,datapoints$stackpos+0.5,datapoints$stackpos),"native"),y=unit(c(datapoints$start+binsize/2,datapoints$start+binsize/2,datapoints$start-binsize/2),"native"),id=rep(1:nrow(datapoints),3),gp=gpar(fill=pal[datapoints$type],lty=0))
#     
#     popViewport()
#   }
#   
#   chrposbins <- gr[seqnames(gr) == chrdisp[i] & strand(gr) == "+"]
#   datapoints <- data.frame(type=unlist(chrposbins$hitvec),start=rep(start(chrposbins),unlist(lapply(chrposbins$hitvec,length))),end=rep(end(chrposbins),unlist(lapply(chrposbins$hitvec,length))),stackpos=unlist(lapply(lapply(chrposbins$hitvec,length),function(x) { if (x>0) seq(1,x) })))
#   
#   if (nrow(datapoints) > 0) {
#     
#     
#     xmax <- max(unlist(lapply(chrposbins$hitvec,length)))
#     
#     dataheight <- convertY(unit(binsize,"native"),"mm")
#     xwidth <- min(unit(xmax,"native"),xmax*dataheight)
#     
#     pushViewport(viewport(x=unit(chrpos[i],"native")+unit(0.1,"lines"),y=unit(0,"npc"),width=xwidth,height=unit(chrlen[names(chrlen) %in% chrdisp][i],"native"),just=c("left","bottom"),xscale=c(0.5,xmax+0.5),yscale=c(chrlen[names(chrlen) %in% chrdisp][i],0)))
#     
#     grid.polygon(x=unit(c(datapoints$stackpos-0.5,datapoints$stackpos+0.5,datapoints$stackpos),"native"),y=unit(c(datapoints$start+binsize/2,datapoints$start+binsize/2,datapoints$start+1.5*binsize),"native"),id=rep(1:nrow(datapoints),3),gp=gpar(fill=pal[datapoints$type],lty=0))
#     
#     popViewport()
#   }
#   
#   
# }
# 
# dev.off()
