#!/usr/bin/Rscript

ARGS <- c(
  "tlxfile", "character", "file path of bowtie ref tags",
  "output","character", "file path to plot to"
)

OPTS <- c(
  "binsize","integer",2000000,"bps per bin"
  "assembly","character","mm9","genome assembly"
)


source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local("Rsub.R")

#suppressPackageStartupMessages()
library(Rsamtools)
library(GenomicRanges)


chrlen <- read.delim(paste(Sys.getenv('GENOME_DB'),assembly,'annotation/ChromInfo.txt',sep="/"),header=F,col.names=c('Name','Length'))
chrlen <- structure(chrlen$Length,names=as.character(chrlen$Name))



chrs <- rep(names(chrlen),c(ceiling(chrlen/bin)))
strands <- rep(c("+","-"),c(length(chrs),length(chrs)))
starts <- unlist(lapply(chrlen,function(len){seq(from=1,to=len,by=bin)}))
ends <- unlist(lapply(chrlen,function(len){if (len > bin) c(seq(from=bin,to=len,by=bin),len) else len } ))




gr <- GRanges(seqnames=rep(chrs,2),ranges=IRanges(start=rep(starts,2),end=rep(ends,2)),strand=strands,seqlengths=chrlen)

con  <- file(tlxfile, open = "r")
header <- unlist(strsplit(readLines(con, n = 1),"\t"))
close(con)

headersToRead <- c("Rname","Junction","Strand")

colClasses <- rep("NULL",length(header))
colClasses[match(headersToRead,header)] <- NA

tlx <- read.delim(tlxfile,header=T,colClasses=colClasses)
tlx <- GRanges(seqnames=tlx$Rname,ranges=IRanges(start=tlx$Junction,end=tlx$Junction),strand=ifelse(tlx$Strand==1,"+","-"))

gr$hits <- countOverlaps(gr,tlx)