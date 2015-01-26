#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfile", "character", "",
    "bedfile","character","",
    "output","character", ""
  )
  
  OPTS <- c(
  
  )

  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  parseArgs("TranslocRepeatSeq.R", ARGS, OPTS)
  
} else {
  #   source("~/TranslocPipeline/R/Rsub.R")
  #   source("~/TranslocPipeline/R/TranslocHelper.R")
  
}

suppressPackageStartupMessages(library(data.table, quietly=TRUE))
suppressPackageStartupMessages(library(dplyr, quietly=TRUE))
suppressPackageStartupMessages(library(GenomicRanges, quietly=TRUE))
suppressPackageStartupMessages(library(rtracklayer, quietly=TRUE))

bed <- import.bed(bedfile)

tlx <- fread(tlxfiles[tlxfile],sep="\t",header=T,select=c("Rname","Junction","Strand","junction"))

tlx <- tlx %>% group_by(Qname) %>% mutate(JuncID = 1:n()) %>% filter(tlx, junction)

gr <-  with(tlx,GRanges(seqnames=Rname,ranges=IRanges(start=Junction,width=1,names=Qname),strand=Strand))

tlx <- mutate(tlx, Overlap = overlapsAny(gr,bed,type="any"))

tlx <- filter(tlx, Overlap) %>% select(Qname,JuncID)

write.table(tlx,output,sep="\t",row.names=F,quote=F,na="")

