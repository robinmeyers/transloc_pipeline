#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfiles", "character", "comma-separated list of files or dir and will grab all *.tlx",
    "outdir","character", "file path to plot to"
  )
  
  OPTS <- c(
    "metafile","character","","explicitly set breaksite, otherwise guess",
    "by.percent","character","","",
    "by.distance","character","1,10,200,1000","Kb",
    "tlxlabels","character","","",
    "resection.only","logical",F,""
    
  )
  
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  parseArgs("TLXBrkChr.R", ARGS, OPTS)
  
} else {
  #   source("~/TranslocPipeline/R/Rsub.R")
  #   source("~/TranslocPipeline/R/TranslocHelper.R")
  metafile <- "./metadata.txt"
  tlxfiles <- "./"
  outdir <- "./BrkChr"
  
  by.percent <- ""
  by.distance <- "1,10,200,1000"
  tlxlabels <- ""
  resection.only <- F

  
}

suppressPackageStartupMessages(library(data.table, quietly=TRUE))
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

by.distance <- sort(as.numeric(unlist(strsplit(by.distance,","))))
by.percent <- sort(as.numeric(unlist(strsplit(by.percent,","))))

if (metafile != "") {
  meta <- fread(metafile,sep="\t",header=T)
  meta <- mutate(meta,Name = paste(Library,"_",Sequencing,sep=""))
}

for (tlxfile in names(tlxfiles)) {

  if (metafile != "") {
    metarow <- grep(tlxfile,meta$Name)
  }

  if (metafile != "" && length(metarow) > 0) {
    brk.chr <- meta$Chr[metarow[1]]
    brk.strand <- ifelse(meta$Strand[metarow[1]] == "+", 1, -1)
    brk.site <- ifelse(brk.strand == 1,meta$End[metarow[1]],meta$Start[metarow[1]])
  } else {
    # Code in guessing at breaksite if no metadata
    stop("Error: need metadata file for now")
  }

  if (!file.exists(file.path(outdir,tlxfile))) {
    dir.create(file.path(outdir,tlxfile))
  }
  
  tlx <- fread(tlxfiles[tlxfile],sep="\t",header=T)
  
  tlx.off <- filter(tlx,Rname != brk.chr)
  tlx <- filter(tlx,Rname == brk.chr)

  outfile <- file.path(outdir,tlxfile,paste(tlxfile,"_InterChr.tlx",sep=""))
  write.table(tlx.off,outfile,sep="\t",row.names=F,quote=F,na="")

  tlx <- mutate(tlx, Distance2Break = Junction - brk.site)

  if (resection.only) {

    tlx.nonresect <- tlx %>%
#       filter(!(Strand == brk.strand & abs(Distance2Break)/Distance2Break == brk.strand)) %>%
      filter(!(abs(Distance2Break)/Distance2Break == brk.strand)) %>%
      mutate(Distance2Break = abs(Distance2Break))
    tlx <- tlx %>%
      filter(abs(Distance2Break)/Distance2Break == brk.strand) %>%
      mutate(Distance2Break = abs(Distance2Break))

    outfile <- file.path(outdir,paste(tlxfile,"_BrkChrNonResection.tlx",sep=""))
    write.table(tlx.nonresect,outfile,sep="\t",row.names=F,quote=F,na="")

  } else {
    tlx <- mutate(tlx,Distance2Break = abs(Distance2Break))
  }

  
  if (length(by.percent) > 0) {
    distance.quant <- quantile(tlx$Distance2Break,by.percent,na.rm=T)
    distance.labels <- c()
    for (i in 1:length(by.percent)) {
      if (i == 1) {
        distance.labels <- c(distance.labels,paste("0-",100*by.percent[i],"pc",sep=""))
      } else {
        distance.labels <- c(distance.labels,paste(100*by.percent[i-1],"-",100*by.percent[i],"pc",sep=""))
      }
    }
    remainder.label <- paste(100*by.percent[i],"-100pc",sep="")
  } else {
    distance.quant <- by.distance*1000
    distance.labels <- c()
    for (i in 1:length(distance.quant)) {
      if (i == 1) {
        distance.labels <- c(distance.labels,paste("0-",by.distance[i],"kb",sep=""))
      } else {
        distance.labels <- c(distance.labels,paste(by.distance[i-1],"-",by.distance[i],"kb",sep=""))
      }
    }
    remainder.label <- paste("gt",by.distance[i],"kb",sep="")
  }
  
  distances <- data.frame(quant = distance.quant, name = distance.labels)
  
  gg <- ggplot(mutate(tlx,Distance2Break = ifelse(Distance2Break==0,1,Distance2Break)),aes(x=Distance2Break,y=..density..)) + geom_histogram() + scale_x_log10()
  gg <- gg + geom_vline(data=distances,mapping=aes(xintercept=quant),color="blue")
  gg <- gg + geom_text(data=distances,mapping=aes(label=name,x=quant,y=0),angle=90,vjust=0,hjust=0,color="red")
  gg <- gg + ggtitle(tlxfile)
  
  pdf(file.path(outdir,paste(tlxfile,"_dis2brk.pdf",sep="")))
  print(gg)
  dev.off()
  
  for (i in 1:nrow(distances)) {
    tlx.brk <- filter(tlx, Distance2Break <= distances$quant[i])
    tlx <- filter(tlx, Distance2Break > distances$quant[i])
    outfile <- file.path(outdir,tlxfile,paste(tlxfile,"_BrkChr_",distances$name[i],".tlx",sep=""))
    write.table(tlx.brk,outfile,sep="\t",row.names=F,quote=F,na="")
  }
  
  outfile <- file.path(outdir,tlxfile,paste(tlxfile,"_BrkChr_",remainder.label,".tlx",sep=""))
  
  write.table(tlx,outfile,sep="\t",row.names=F,quote=F,na="")
                             
}
