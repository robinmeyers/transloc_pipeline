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
  tlxfile <- "/Volumes//AltLab/Translocation/NewPipeline/Alt046/results/PW002_Alt046/PW002_Alt046.tlx"
  output <- "~/Working/NewPipelineValidations/DedupTesting/PW002_dedup.txt"
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


findDuplicates <- function(n,tlxs) {
  tlx <- tlxs[n,]
  matches <- subset(tlxs, Qname>tlx$Qname &
                      abs(Offset-tlx$Offset)<=qdist & abs(Junction-tlx$Junction)<=rdist &
                      abs(B_Offset-tlx$B_Offset)<=qdist & abs(B_Junction-tlx$B_Junction)<=rdist)
  if (nrow(matches) > 0) {
    return(paste(paste(matches$Qname,"(",matches$B_Junction-tlx$B_Junction,",",matches$Junction-tlx$Junction,")",sep="")[1:3],collapse=","))
  } else {
    return("")
  }
}

if (cores == 0) {
  cores <- detectCores()
}

cat("Deduplicating junctions on",cores,"cores\n")
print(object.size(tlxs_by_chr_and_strand),units="Mb")
tlxs <- ldply(lapply(1:length(tlxs_by_chr_and_strand),function (n,tlxs) {
  if (nrow(tlxs[[n]]) > 0) {
    
    cat(nrow(tlxs[[n]])," - ")
    print(object.size(tlxs[[n]]),units="Mb")
    
    dups <- mclapply(1:nrow(tlxs[[n]]),findDuplicates,tlxs[[n]][,c("Qname","Offset","Junction","B_Offset","B_Junction")],mc.cores=cores)
    
    print(object.size(dups),units="Mb")
    
    tlxs[[n]]$Dups <- unlist(dups)
    return(tlxs[[n]])
  }
},tlxs_by_chr_and_strand))

write.table(tlxs[tlxs$Dups != "",c("Qname","Dups")],output,sep="\t",quote=F,na="",row.names=F,col.names=F)
