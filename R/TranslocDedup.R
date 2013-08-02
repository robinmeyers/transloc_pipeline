#!/usr/bin/Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfile", "character", " ",
    "output","character", " "
  )
  
  OPTS <- c(
    "qdist","numeric",2," ",
    "rdist","numeric",10," ",
    "cores","numeric",4,"Number of compute nodes to run on"
  )
  
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  parseArgs("TranslocDedup.R", ARGS, OPTS)
  
} else {
  source("~/Pipelines/R/Rsub.R")
  tlxfile <- "/Volumes/AltLab/Translocation/RawData/Alt024-20130429/NewPipelineTest/test/CC004_Alt024/CC004_Alt024.tlx"
  output <- "/Volumes/AltLab/Translocation/RawData/Alt024-20130429/NewPipelineTest/test/CC004_Alt024/CC004_Alt024_dedup.txt"
  rdist <- 10
  qdist <- 2
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

tlxs$Offset <- with(tlxs,ifelse(Strand==1,Junction-Qstart,Junction+Qstart))
tlxs$B_Junction <- with(tlxs,ifelse(B_Strand==1,B_Rend,B_Rstart))
tlxs$B_Offset <- with(tlxs,ifelse(B_Strand==1,B_Junction-B_Qend,B_Junction+B_Qend))

tlxs_by_chr <- split(tlxs,tlxs$Rname)

findDuplicates <- function(tlx,tlxs,tlxs_by_chr) {
  tlx <- tlxs[tlx,]
  tlxs_by_chr <- tlxs_by_chr[[tlx$Rname]]
  matches <- subset(tlxs_by_chr, Qname>tlx$Qname & Strand==tlx$Strand &
                      abs(Offset-tlx$Offset)<=qdist & abs(Junction-tlx$Junction)<=rdist &
                      B_Rname == tlx$B_Rname & B_Strand == tlx$B_Strand &
                      abs(B_Offset-tlx$B_Offset)<=qdist & abs(B_Junction-tlx$B_Junction)<=rdist)
  if (nrow(matches) > 0) {
    return(paste(matches$Qname,collapse=","))
  } else {
    return("")
  }
}

dups <- mclapply(1:nrow(tlxs),findDuplicates,tlxs,tlxs_by_chr,mc.cores=cores)
tlxs$Dups <- unlist(dups)
write.table(tlxs[tlxs$Dups != "",c("Qname","Dups")],output,sep="\t",quote=F,na="",row.names=F,col.names=F)
