#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfile", "character", " ",
    "output","character", " "
  )
  
  OPTS <- c(
    
    "cores","numeric",1,"Number of compute nodes to run on"
  )
  
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  
} else {
  source("~/Pipelines/R/Rsub.R")
  tlxfile <- "/Volumes/AltLab/Translocation/RawData/Alt024-20130429/NewPipelineTest/test/CC004_Alt024/CC004_Alt024.tlx"
  output <- "/Volumes/AltLab/Translocation/RawData/Alt024-20130429/NewPipelineTest/test/CC004_Alt024/CC004_Alt024_dedup.txt"
  rdif <- 10
  qdif <- 2
  cores <- 4
}

library(plyr)
library(parallel)

con  <- file(tlxfile, open = "r")
header <- unlist(strsplit(readLines(con, n = 1),"\t"))
close(con)

headersToSkip <- c("Seq")

colClasses <- rep(NA,length(header))
colClasses[match(headersToSkip,header)] <- "NULL"

tlxs <- read.delim(tlxfile,header=T,colClasses=colClasses,as.is=T)
tlxs <- tlxs[with(tlxs,order(Rname,Strand,Junction)),]
tlxs_by_chr <- split(tlxs,tlxs$Rname)

findDuplicates <- function(tlx,tlxs,tlxs_by_chr) {
  tlx <- tlxs[tlx,]
  tlxs_by_chr <- tlxs_by_chr[[tlx$Rname]]
  matches <- subset(tlxs_by_chr, Qname > tlx$Qname & Strand == tlx$Strand & Junction == tlx$Junction)
  if (nrow(matches) > 0) {
    return(paste(matches$Qname,sep=","))
  } else {
    return("")
  }
}

dups <- mclapply(1:nrow(tlxs),findDuplicates,tlxs,tlxs_by_chr,mc.cores=cores)
tlxs$Dups <- unlist(dups)
write.table(tlxs[tlxs$Dups != "",c("Qname","Dups")],output,sep="\t",quote=F,na="",row.names=F,col.names=F)
