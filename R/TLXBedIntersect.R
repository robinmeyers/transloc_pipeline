#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfiles", "character", "comma-separated list of files or dir and will grab all *.tlx",
    "bedfile","character","",
    "outdir","character", "file path to plot to"
  )
  
  OPTS <- c(

    "tlxlabels","character","",""

    
  )
  
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  parseArgs("TLXBedIntersect.R", ARGS, OPTS)
  
} else {
  #   source("~/TranslocPipeline/R/Rsub.R")
  #   source("~/TranslocPipeline/R/TranslocHelper.R")
  tlxfiles <- "../richard/brkchr/Cas9A_combined/Cas9A_combined_InterChr.tlx"
  outdir <- "./"
  bedfile <- "./Cas9A_OT.bed"
  
  tlxlabels <- ""
  
}

suppressPackageStartupMessages(library(data.table, quietly=TRUE))
suppressPackageStartupMessages(library(dplyr, quietly=TRUE))
suppressPackageStartupMessages(library(GenomicRanges, quietly=TRUE))
suppressPackageStartupMessages(library(rtracklayer, quietly=TRUE))
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

bedname <- sub(".bed","",basename(bedfile))
bed <- import.bed(bedfile)


for (tlxfile in names(tlxfiles)) {

  
  tlx <- fread(tlxfiles[tlxfile],sep="\t",header=T)
  
  gr <-  with(tlx,GRanges(seqnames=Rname,ranges=IRanges(start=Junction,width=1,names=Qname),strand=Strand))

  tlx <- mutate(tlx, Overlap = overlapsAny(gr,bed,type="any"))

  tlx.in <- filter(tlx, Overlap == T)
  tlx.out <- filter(tlx, Overlap == F)

  in.outfile <- file.path(outdir,paste(tlxfile,"_",bedname,".tlx",sep=""))
  write.table(tlx.in,in.outfile,sep="\t",row.names=F,quote=F,na="")

  out.outfile <- file.path(outdir,paste(tlxfile,"_complement.tlx",sep=""))
  write.table(tlx.out,out.outfile,sep="\t",row.names=F,quote=F,na="")
                             
}
