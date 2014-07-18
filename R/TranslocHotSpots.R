#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfile", "character", " ",
    "output","character", " "
  )
  
  OPTS <- c(
    "assembly","character","mm9"," ",
    "bin.width","numeric",500," ",
    "bg.width","numeric",500000," ",
    "alpha","numeric",0.01," ",
    "hits.min","numeric",3," ",
    "strand.min","numeric",1," ",
    "cores","numeric",4," "
  )
  
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  source_local("TranslocHelper.R")
  parseArgs("TranslocHotSpots.R", ARGS, OPTS)
  
} else {
  source("~/TranslocPipeline//R/Rsub.R")
  source("~/TranslocPipeline/R/TranslocHelper.R")
  
  tlxfile <-  "~/Working/hotspots//ATM_4in1_G30.tlx"
  output <- "~/Working/hotspots//ATM_4in1_G30"
  
  assembly <- "mm9"
  bin.width <- 500
  bg.width <- 500000
  alpha <- 0.01
  hits.min <- 3
  strand.min <- 1
  cores <- 4
  
}

suppressPackageStartupMessages(library(rtracklayer, quietly=TRUE))
suppressPackageStartupMessages(library(GenomicRanges, quietly=TRUE))

tlx <- readTLX(tlxfile,columnsToRead=c("Qname","Rname","Junction","Strand"))
tlx <- tlx[(with(tlx,order(Rname,Junction))),]

chrlen <- getChromLens(assembly)
chr.stats <- data.frame(len=chrlen)

tlx.gr <- tlxToGR(tlx,chrlen)

# Create Initial GR Bins

bins1 <- suppressWarnings(resize(tlx.gr[distance(tlx.gr,tlx.gr[nearest(tlx.gr)]) <= bin.width/2],width=bin.width,fix="center"))
# 
# bins1 <- tlx.gr[2:length(tlx.gr)]
# start(bins1) <- start(tlx.gr[1:(length(tlx.gr)-1)])
# start(bins1) <- pmax(1,ifelse(strand(tlx.gr) == "+",start(tlx.gr)-bin.width/2,start(tlx.gr)-bin.width/2+1))
# end(bins1) <- pmin(seqlengths(tlx.gr)[as.character(seqnames(tlx.gr))],ifelse(strand(tlx.gr) == "+",start(tlx.gr)+bin.width/2-1,start(tlx.gr)+bin.width/2))
strand(bins1) <- "*"

bins1 <- unique(bins1)

bins1.by.chr <- split.by.chr(bins1)
tlx.gr.by.chr <- split.by.chr(tlx.gr)

tlx.cov <- coverage(tlx.gr)

tlx.cov.cumsum <- lapply(tlx.cov,cumsum)

if (any(names(bins1.by.chr) != names(tlx.cov))) stop("Error: chromosome names out of order")

bins1.by.chr <- GRangesList(mcmapply(calculate.local.significance,bins=bins1.by.chr,tlx.cumsum=tlx.cov.cumsum,bin.width=bin.width,bg.width=bg.width,mc.cores=cores))
bins1.by.chr <- GRangesList(mcmapply(calculate.chr.significance,bins=bins1.by.chr,tlx.cumsum=tlx.cov.cumsum,bin.width=bin.width,mc.cores=cores))

# bins1 <- unlist(bins1.by.chr)


# Trim peaks to last translocation on each edge

peaks.local <- GRangesList(mcmapply(function(b,c) {
  trimPeaks(reduce(b[b$p.local < alpha]),c)
}, bins1.by.chr, tlx.gr.by.chr,mc.cores=cores))

peaks.chr <- GRangesList(mcmapply(function(b,c) {
  trimPeaks(reduce(b[b$p.chr < alpha]),c)
}, bins1.by.chr, tlx.gr.by.chr,mc.cores=cores))

peaks.local <- GRangesList(mcmapply(calculate.local.significance,bins=peaks.local,tlx.cumsum=tlx.cov.cumsum,bin.width=bin.width,bg.width=bg.width,mc.cores=cores))
peaks.chr <- GRangesList(mcmapply(calculate.chr.significance,bins=peaks.chr,tlx.cumsum=tlx.cov.cumsum,bin.width=bin.width,mc.cores=cores))

peaks.local <- unlist(peaks.local)
peaks.chr <- unlist(peaks.chr)

peaks.local$q.value <- p.adjust(peaks.local$p.local,method="BY")
peaks.chr$q.value <- p.adjust(peaks.chr$p.chr,method="BY")

peaks.local <- peaks.local[peaks.local$q.value <= alpha]
peaks.chr <- peaks.chr[peaks.chr$q.value <= alpha]

peaks.local$hits <- countOverlaps(peaks.local,tlx.gr)
peaks.local$hits.plus <- countOverlaps(peaks.local,tlx.gr[strand(tlx.gr) == "+"])
peaks.local$hits.minus <- countOverlaps(peaks.local,tlx.gr[strand(tlx.gr) == "-"])

peaks.chr$hits <- countOverlaps(peaks.chr,tlx.gr)
peaks.chr$hits.plus <- countOverlaps(peaks.chr,tlx.gr[strand(tlx.gr) == "+"])
peaks.chr$hits.minus <- countOverlaps(peaks.chr,tlx.gr[strand(tlx.gr) == "-"])


peaks.local <- peaks.local[ peaks.local$hits >= hits.min &
                              peaks.local$hits.plus >= strand.min &
                              peaks.local$hits.minus >= strand.min ]
peaks.chr <- peaks.chr[ peaks.chr$hits >= hits.min &
                          peaks.chr$hits.plus >= strand.min &
                          peaks.chr$hits.minus >= strand.min ]

peaks.local$p.local <- log10(peaks.local$p.local)
peaks.chr$p.chr <- log10(peaks.chr$p.chr)

peaks.local$q.value <- log10(peaks.local$q.value)
peaks.chr$q.value <- log10(peaks.chr$q.value)

export.bed(peaks.local,paste(output,"_localpeaks.bed",sep=""),ignore.strand=T)
export.bed(peaks.chr,paste(output,"_chrpeaks.bed",sep=""),ignore.strand=T)

write.table(as.data.frame(peaks.local,row.names=NULL),paste(output,"_localpeaks.txt",sep=""),quote=F,sep="\t",row.names=F)
write.table(as.data.frame(peaks.chr,row.names=NULL),paste(output,"_chrpeaks.txt",sep=""),quote=F,sep="\t",row.names=F)

