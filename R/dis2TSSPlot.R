#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tss.file", "character",""
  )
  
  OPTS <- c(
    "pdf.file","character","","",
    "bin.width","integer",1000,"",
    "x.min","integer",-50000,"",
    "x.max","integer",50000,""
  )
  
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  parseArgs("dis2TSS.R", ARGS, OPTS)
  
} else {
  #   source("~/TranslocPipeline/R/Rsub.R")
  #   source("~/TranslocPipeline/R/TranslocHelper.R")
  tss.file <- "~/AltLab/Simulations/mapqual/results/dis2TSS/results_dis2TSS.txt"
  pdf.file <- ""
  bin.width <- 2000
  x.min <- -100000
  x.max <- 100000
}

suppressPackageStartupMessages(library(data.table, quietly=TRUE))
suppressPackageStartupMessages(library(ggplot2, quietly=TRUE))

if (pdf.file == "") {
  pdf.file <- sub(paste(file_ext(tss.file),"$",sep=""),"pdf",tss.file)
}

dis2tss <- fread(tss.file,header=T,sep="\t")

gg <- ggplot(dis2tss,aes(x=Dis2TSS,y=..density..,color=Library))
gg <- gg + geom_freqpoly(binwidth=bin.width)
gg <- gg + scale_x_continuous(limits = c(x.min,x.max))

pdf(pdf.file,width=10,height=7)
print(gg)
dev.off()
