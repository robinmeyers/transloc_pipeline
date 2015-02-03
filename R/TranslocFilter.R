#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfile", "character", "",
    "output","character",""
  )
  
  OPTS <- c(
    "as.filter","logical",TRUE,"set to FALSE to use as.selector",
    "f.unaligned","character","","",
    "f.baitonly","character","","",
    "f.uncut","character","","",
    "f.misprimed","character","","",
    "f.freqcut","character","","",
    "f.largegap","character","","",
    "f.mapqual","character","","",
    "f.breaksite","character","","",
    "f.sequential","character","","",
    "f.repeatseq","character","","",
    "f.duplicate","character","",""
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
  

  tlxfile <- "~/Working/TranslocTesting/results_mid/RF204_Alt055//RF204_Alt055.tlx"
  output <- "~/Working/TranslocTesting/results_mid/RF204_Alt055/RF204_Alt055_filtered.txt"

  as.filter <- TRUE
  f.unaligned <- "G0"
  f.baitonly <- "G0"
  f.uncut <- "G0"
  f.misprimed <- "L10"
  f.freqcut <- "G0"
  f.largegap <- "G30"
  f.mapqual <- "L10"
  f.breaksite <- "G0"
  f.sequential <- "G0"
  f.repeatseq <- "G0"
  f.duplicate <- "G0"
}

suppressPackageStartupMessages(library(data.table, quietly=TRUE))
suppressPackageStartupMessages(library(dplyr, quietly=TRUE))

stats.file <- sub(paste(".",file_ext(output),sep=""),"_stats.txt",output)


filter.names <- c("unaligned","baitonly","uncut","misprimed","freqcut","largegap","mapqual","breaksite","sequential","repeatseq","duplicate")
filter.values <- c()

for (filter.name in filter.names) {
  tmp.value <- get(paste("f.",filter.name,sep=""))
  if (grepl("^[0-9]+$",tmp.value)) {
    filter.values[filter.name] <- paste("==",tmp.value,sep="")
    
  } else if (grepl("^[GL]?[E]?[0-9]+$",tmp.value)) {
    tmp.value <- sub("[Gg][Ee]",">=",tmp.value)
    tmp.value <- sub("[Ll][Ee]","<=",tmp.value)
    tmp.value <- sub("[Gg]",">",tmp.value)
    tmp.value <- sub("[Ll]","<",tmp.value)
    tmp.value <- sub("[Ee]","==",tmp.value)
    
    filter.values[filter.name] <- tmp.value
  } else if (grepl("^[0-9]+$",tmp.value)) {
    filter.values[filter.name] <- paste("==",tmp.value,sep="")
  } else if (tmp.value != "") {
    stop(paste("Error:",filter.name,"filter entered in wrong format"))
  }
}


tlx <- fread(tlxfile,sep="\t",header=T,select=c("Qname","JuncID","Rname",filter.names))
tlx <- filter(tlx,!(JuncID > 1 & Rname == "Adapter"))

stats.names <- c("total",filter.names,"result")
filter.stats <- data.frame(reads=rep(0,length(stats.names)),junctions=rep(0,length(stats.names)),row.names=stats.names)
                           
reads.total <- as.integer(summarize(tlx,n_distinct(Qname)))
junctions.total <- nrow(filter(tlx,Rname != "" & Rname != "Adapter"))

filter.stats["total",] <- c(reads.total,junctions.total)

tlx.filt.list <- list()

for (filter.name in filter.names) {
  filter.value <- filter.values[filter.name]
  if (filter.value == "") next

  filter.text <- paste(filter.name,filter.value,sep="")

  tlx.filt.list[[filter.name]] <- filter(tlx,eval(parse(text=filter.text))) %>% select(Qname,JuncID,Rname)
  
  junctions.count <- nrow(filter(tlx.filt.list[[filter.name]],Rname != "" & Rname != "Adapter"))
  reads.count <- as.integer(summarize(tlx.filt.list[[filter.name]],n_distinct(Qname)))
  
  filter.stats[filter.name,] <- c(reads.count,junctions.count)

}

tlx.filt <- as.data.table(rbind_all(tlx.filt.list))

if (as.filter) {
  tlx <- anti_join(tlx,tlx.filt,by=c("Qname","JuncID"))

} else {
  tlx <- semi_join(tlx,tlx.filt,by=c("Qname","JuncID"))
}

reads.count <- as.integer(summarize(tlx,n_distinct(Qname)))
junctions.count <- nrow(filter(tlx,Rname != "" & Rname != "Adapter"))
filter.stats["result",] <- c(reads.count,junctions.count)

tlx <- select(tlx,Qname,JuncID)

write.table(tlx,output,sep="\t",col.names=F,row.names=F,quote=F,na="")
write.table(filter.stats,stats.file,sep="\t",col.names=NA,row.names=T,quote=F,na="")

