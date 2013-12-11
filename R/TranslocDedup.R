#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfile", "character", " ",
    "output","character", " "
  )
  
  OPTS <- c(
    "qdist","numeric",2," ",
    "rdist","numeric",10," ",
    "cores","numeric",0,"Number of compute nodes to run on"
  )
  
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  parseArgs("TranslocDedup.R", ARGS, OPTS)
  
} else {
  source("~/TranslocPipeline//R/Rsub.R")
  tlxfile <- "~/Working/NewPipelineValidations/DedupTesting/YZ101_r.tlx"
  output <- "~/Working/NewPipelineValidations/DedupTesting/YZ101_r_dedup.txt"
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

tlxs_by_chr_and_strand <- split(tlxs,list(tlxs$Rname,tlxs$Strand))

findDuplicates <- function(n,tlxs,split_tlxs) {
  tlx <- tlxs[n,]
  tlxs_by_chr_and_strand <- split_tlxs[[paste(tlx$Rname,tlx$Strand,sep='.')]]
  matches <- subset(tlxs_by_chr_and_strand, Qname>tlx$Qname &
                      abs(Offset-tlx$Offset)<=qdist & abs(Junction-tlx$Junction)<=rdist &
                      abs(B_Offset-tlx$B_Offset)<=qdist & abs(B_Junction-tlx$B_Junction)<=rdist)
  if (nrow(matches) > 0) {
    return(paste(matches$Qname,"(",matches$B_Junction-tlx$B_Junction,",",matches$Junction-tlx$Junction,")",sep="",collapse=","))
  } else {
    return("")
  }
}

if (cores == 0) {
  cores <- detectCores()
}

cat("Deduplicating junctions on",cores,"cores\n")
dups <- mclapply(1:nrow(tlxs),findDuplicates,tlxs,tlxs_by_chr_and_strand,mc.cores=cores)
tlxs$Dups <- unlist(dups)
write.table(tlxs[tlxs$Dups != "",c("Qname","Dups")],output,sep="\t",quote=F,na="",row.names=F,col.names=F)
