#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfile", "character", "",
    "output","character",""
  )
  
  OPTS <- c(
    "filter","logical",TRUE,"set to TRUE to start with all and filter out; set to FALSE to start with none and filter in",
    "f.unaligned","integer",0,"set to 1 to keep unaligned",
    "f.baitonly","integer",0,"set to 1 to keep baitonly",
    "f.uncut","integer",0,"set to 1 to keep uncut",
    "f.misprimed","integer",10,"min bases aligned after primer; set to 0 to keep all",
    "f.freqcut","integer",0,"set to 1 to keep freqcut",
    "f.largegap","integer",30,"max bases between bait and prey; set to 0 to keep all",
    "f.mapqual","integer",10,"min mapqual score; set to 0 to keep all",
    "f.breaksite","integer",0,"set to 1 to keep breaksite",
    "f.sequential","integer",0,"set to 1 to keep sequential",
    "f.repeatseq","integer",0,"set to 1 to keep repeatseq",
    "f.duplicate","integer",0,"set to 1 to keep duplicate"
  )

  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  parseArgs("TranslocFilter.R", ARGS, OPTS)
  
} else {
  source("~/TranslocPipeline/R/Rsub.R")
  source("~/TranslocPipeline/R/TranslocHelper.R")
  

  tlxfile <- ""
  output <- ""

  filter <- TRUE
  f.unaligned <- 0
  f.baitonly <- 0
  f.uncut <- 0
  f.misprimed <- 10
  f.freqcut <- 0
  f.largegap <- 30
  f.mapqual <- 10
  f.breaksite <- 0
  f.sequential <- 0
  f.repeatseq <- 0
  f.duplicate <- 0
}

suppressPackageStartupMessages(library(data.table, quietly=TRUE))
suppressPackageStartupMessages(library(dplyr, quietly=TRUE))


tlx <- fread(tlxfile,sep="\t",header=T)



write.table(tlx,output,sep="\t",col.names=F,row.names=F,quote=F,na="")

