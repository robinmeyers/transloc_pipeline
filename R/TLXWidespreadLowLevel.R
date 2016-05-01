#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfiles", "character", "comma-separated list of files or dir and will grab all *.tlx",
    "outdir","character", "file path to plot to"
  )
  
  OPTS <- c(
    "meta","character","","",
    "by.percent","character","0.02,0.05,0.1","",
    "by.distance","character","","Kb",
    "tlxlabels","character","","",
    "ignore.strand","logical",T,""
    
  )
  
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  parseArgs("TLXWidespreadLowLevel.R", ARGS, OPTS)
  
} else {
  #   source("~/TranslocPipeline/R/Rsub.R")
  #   source("~/TranslocPipeline/R/TranslocHelper.R")
  tlxfiles <- "./"
  outdir <- "./WSLL"
  
  by.percent <- "0.02,0.05,0.1"
  by.distance <- ""
  tlxlabels <- ""
  ignore.strand <- T

  
}

suppressPackageStartupMessages(library(readr, quietly=TRUE))
suppressPackageStartupMessages(library(dplyr, quietly=TRUE))
suppressPackageStartupMessages(library(GenomicRanges, quietly=TRUE))
suppressPackageStartupMessages(library(ggplot2, quietly=TRUE))



if (file.info(tlxfiles)[["isdir"]]) {
  tlxfiles <- list.files(tlxfiles,pattern = "*.tlx",full.names = T)
} else {
  tlxfiles <- unlist(strsplit(tlxfiles,","))
}

names(tlxfiles) <- sub(".tlx","",basename(tlxfiles))

if (tlxlabels != "") {
  tlxlabels <- unlist(strsplit(tlxlabels,","))
  if (length(tlxlabels) != length(tlxfiles)) {
    stop("Error: label list has different length than tlxfile list")
  }
  names(tlxfiles) <- tlxlabels
}

if (!file.exists(outdir)) {
  dir.create(outdir)
}

by.distance <- sort(as.numeric(unlist(strsplit(by.distance,","))),decreasing=T)
by.percent <- sort(as.numeric(unlist(strsplit(by.percent,","))))

for (tlxfile in names(tlxfiles)) {
  
  if (!file.exists(file.path(outdir,tlxfile))) {
    dir.create(file.path(outdir,tlxfile))
  }
  
  tlx <- read_tsv(tlxfiles[tlxfile])
  
  gr <-  with(tlx,GRanges(seqnames=Rname,ranges=IRanges(start=Junction,width=1,names=Qname),strand=Strand))
  hits <- distanceToNearest(gr,ignore.strand=ignore.strand)
  tlx$Distance2Nearest <- as.data.frame(hits)$distance
  
  if (length(by.distance) > 0) {
    distance.quant <- by.distance*1000
    distance.labels <- c()
    for (i in 1:length(distance.quant)) {
      if (i == 1) {
        distance.labels <- c(distance.labels,paste("gt",by.distance[i],"kb",sep=""))
      } else {
        distance.labels <- c(distance.labels,paste(by.distance[i],"-",by.distance[i-1],"kb",sep=""))
      }
    }
    remainder.label <- paste("0-",by.distance[i],"kb",sep="")
  } else {
    distance.quant <- quantile(tlx$Distance2Nearest,1-by.percent,na.rm=T)
    distance.labels <- c()
    for (i in 1:length(by.percent)) {
      if (i == 1) {
        distance.labels <- c(distance.labels,paste("0-",100*by.percent[i],"pc",sep=""))
      } else {
        distance.labels <- c(distance.labels,paste(100*by.percent[i-1],"-",100*by.percent[i],"pc",sep=""))
      }
    }
    remainder.label <- paste(100*by.percent[i],"-100pc",sep="")
  }
  
  distances <- data.frame(quant = distance.quant, name = distance.labels)
  
  gg <- ggplot(mutate(tlx,Distance2Nearest = ifelse(Distance2Nearest==0,1,Distance2Nearest)),aes(x=Distance2Nearest,y=..density..)) + geom_histogram() + scale_x_log10()
  gg <- gg + geom_vline(data=distances,mapping=aes(xintercept=quant),color="blue")
  gg <- gg + geom_text(data=distances,mapping=aes(label=name,x=quant,y=0),angle=90,vjust=1,hjust=0,color="red")
  gg <- gg + ggtitle(tlxfile)
  
  pdf(file.path(outdir,paste(tlxfile,"_dis2nearest.pdf",sep="")))
  print(gg)
  dev.off()
  
  for (i in 1:nrow(distances)) {
    tlx.wsll <- filter(tlx, Distance2Nearest >= distances$quant[i])
    tlx <- filter(tlx, Distance2Nearest < distances$quant[i])
    outfile <- file.path(outdir,tlxfile,paste(tlxfile,"_WSLL_",distances$name[i],".tlx",sep=""))
    write.table(tlx.wsll,outfile,sep="\t",row.names=F,quote=F,na="")
  }

outfile <- file.path(outdir,tlxfile,paste(tlxfile,"_WSLL_",remainder.label,".tlx",sep=""))
  write.table(tlx,outfile,sep="\t",row.names=F,quote=F,na="")
                             
}
