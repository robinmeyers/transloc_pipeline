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
  source("~/TranslocPipeline/R/Rsub.R")
  source("~/TranslocPipeline/R/TranslocHelper.R")
  tlxfile <- "./VK020_Alt133.tlx"
  output <- "./VK020_Alt133_repeatseq.txt"
  bedfile <- "/Volumes/AltLab/Genomes/mm9/annotation/repeatSeq.bed"
}

suppressPackageStartupMessages(library(data.table, quietly=TRUE))
suppressPackageStartupMessages(library(dplyr, quietly=TRUE))
suppressPackageStartupMessages(library(GenomicRanges, quietly=TRUE))
suppressPackageStartupMessages(library(rtracklayer, quietly=TRUE))

bed <- import.bed(bedfile)

tlx <- fread(tlxfile,sep="\t",header=T,select=c("Qname","JuncID","Rname","Junction"))

tlx <- filter(tlx, ! is.na(Rname) & ! is.na(Junction))

gr <-  with(tlx,GRanges(seqnames=Rname,ranges=IRanges(start=Junction,width=1,names=Qname)))

ol <- suppressWarnings(overlapsAny(gr,bed,type="any"))

tlx <- mutate(tlx, Overlap = ol) %>% filter(Overlap == T) %>% select(Qname,JuncID)

write.table(tlx,output,sep="\t",col.names=F,row.names=F,quote=F,na="")

