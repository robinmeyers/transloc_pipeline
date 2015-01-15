#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfiles", "character", "comma-separated list of files or dir and will grab all *.tlx",
    "outputstub","character", "file path to plot to"
  )
  
  OPTS <- c(
    "mh.ymax","integer",NA,"",
    "mh.xmin","integer",0,"",
    "mh.xmax","integer",8,"",
    "bait.ymax","integer",NA,"",
    "bait.xmin","integer",NA,"",
    "bait.xmax","integer",NA,"",
    "bait.binsize","integer",1,"",
    "tlxlabels","character","","",
    "normalize","logical",T,"",
    "normalize.by","integer",0,"",
    "rm.ins","logical",T,""
  )
  
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  parseArgs("TLXBaitMH.R", ARGS, OPTS)
  
} else {
#   source("~/TranslocPipeline/R/Rsub.R")
#   source("~/TranslocPipeline/R/TranslocHelper.R")
  tlxfiles <- "~/Working/MH_BaitLen/rm_G30/"
  outputstub <- "~/Working/MH_BaitLen/VK_test"
  
  mh.ymax <- NA
  mh.xmin <- 0
  mh.xmax <- 8
  bait.ymax <- NA
  bait.xmin <- NA
  bait.xmax <- NA
  bait.binsize <- 1
  tlxlabels <- ""
  normalize <- T
  normalize.by <- 0 
  rm.ins <- T
  
}

suppressPackageStartupMessages(library(data.table, quietly=TRUE))
suppressPackageStartupMessages(library(dplyr, quietly=TRUE))
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

tlx <- list()

for (tlxfile in names(tlxfiles)) {
  tlx[[tlxfile]] <- fread(tlxfiles[tlxfile],sep="\t",header=T,select=c("Qname","B_Rstart","B_Rend","B_Strand","B_Qend","Qstart")) %>% mutate(tlxlabel=tlxfile)
}

tlx <- rbind_all(tlx)

bait.anchor <- tlx %>% mutate(BaitAnchor = ifelse(B_Strand == 1, B_Rstart, B_Rend)) %>%
  group_by(tlxlabel,BaitAnchor) %>%
  summarize(count=n()) %>%
  group_by(tlxlabel) %>%
  filter(count == max(count))

tlx <- inner_join(tlx,bait.anchor,by="tlxlabel")

tlx <- mutate(tlx,BaitLen = ifelse(B_Strand == 1, B_Rend-BaitAnchor+1, BaitAnchor-B_Rstart+1 ))
tlx <- mutate(tlx,MicroHom = B_Qend - Qstart + 1)

if (normalize) {
  freqpoly.y <- aes(y=..density..)
} else {
  freqpoly.y <- aes()
}

bait.bins = seq(min(tlx$BaitLen)-0.5,max(tlx$BaitLen)+0.5,by=bait.binsize)


bait.gg <- ggplot(tlx,aes(x=BaitLen,color=tlxlabel))
bait.gg <- bait.gg + geom_freqpoly(freqpoly.y,breaks=bait.bins,lwd=1)
bait.gg <- bait.gg + ylim(c(0,bait.ymax)) + xlim(c(bait.xmin,bait.xmax))
bait.gg <- bait.gg + xlab("Bait Length") + theme(legend.title = element_blank())


if (rm.ins) {
  tlx <- filter(tlx,MicroHom >= 0)
}

mh.bins = seq(min(tlx$MicroHom)-0.5,max(tlx$MicroHom)+0.5,by=1)

mh.gg <- ggplot(tlx,aes(x=MicroHom,color=tlxlabel))
mh.gg <- mh.gg + geom_freqpoly(freqpoly.y,breaks=mh.bins,lwd=1)
mh.gg <- mh.gg + ylim(c(0,mh.ymax)) + xlim(c(mh.xmin,mh.xmax))
mh.gg <- mh.gg + xlab("Micro-Homology") + theme(legend.title = element_blank())


pdf(paste(outputstub,"_baitlen.pdf",sep=""))
print(bait.gg)
dev.off()
pdf(paste(outputstub,"_mh.pdf",sep=""))
print(mh.gg)
dev.off()