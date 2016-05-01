#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfiles", "character", "comma-separated list of files or dir and will grab all *.tlx",
    "outputstub","character", "file path to plot to"
  )
  
  OPTS <- c(
    "mh.ymax","integer",NA_integer_,"",
    "mh.xmin","integer",0,"",
    "mh.xmax","integer",8,"",
    "bait.ymax","integer",NA_integer_,"",
    "bait.xmin","integer",NA_integer_,"",
    "bait.xmax","integer",NA_integer_,"",
    "bait.binsize","integer",1,"",
    "tlxlabels","character","","",
    "normalize","logical",T,"",
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
  
  mh.ymax <- NA_integer_
  mh.xmin <- 0
  mh.xmax <- 8
  bait.ymax <- NA_integer_
  bait.xmin <- NA_integer_
  bait.xmax <- NA_integer_
  bait.binsize <- 1
  tlxlabels <- ""
  normalize <- T
  rm.ins <- T
  
}

suppressPackageStartupMessages(library(readr, quietly=TRUE))
suppressPackageStartupMessages(library(dplyr, quietly=TRUE))
suppressPackageStartupMessages(library(grid, quietly=TRUE))
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
  tlx[[tlxfile]] <- read_tsv(tlxfiles[tlxfile],
                           col_types=cols_only("Qname" = "c",
                                               "B_Rstart" = "i",
                                               "B_Rend" = "i",
                                               "B_Strand" = "i",
                                               "B_Qend"="i",
                                               "Qstart"="i")) %>% 
    mutate(tlxlabel=tlxfile)
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

tlx <- tlx %>% group_by(tlxlabel) %>% summarize(bait.total = prettyNum(n(),big.mark=",")) %>% inner_join(tlx,by="tlxlabel")
bait.bins = seq(min(tlx$BaitLen)-0.5,max(tlx$BaitLen)+0.5,by=bait.binsize)

bait.gg <- ggplot(tlx,aes(x=BaitLen,color=paste(tlxlabel,"\nn=",bait.total,sep="")))
bait.gg <- bait.gg + geom_freqpoly(freqpoly.y,breaks=bait.bins,lwd=1)
bait.gg <- bait.gg + ylim(c(0,bait.ymax)) + xlim(c(bait.xmin,bait.xmax))
bait.gg <- bait.gg + xlab("Bait Length") + theme(legend.title = element_blank(),legend.key.height=unit(2,"line"))


if (rm.ins) {
  tlx <- filter(tlx,MicroHom >= 0)
}

tlx <- tlx %>% group_by(tlxlabel) %>% summarize(mh.total = prettyNum(n(),big.mark=",")) %>% inner_join(tlx,by="tlxlabel")
mh.bins = seq(min(tlx$MicroHom)-0.5,max(tlx$MicroHom)+0.5,by=1)

mh.gg <- ggplot(tlx,aes(x=MicroHom,color=paste(tlxlabel,"\nn=",mh.total,sep="")))
mh.gg <- mh.gg + geom_freqpoly(freqpoly.y,breaks=mh.bins,lwd=1)
mh.gg <- mh.gg + ylim(c(0,mh.ymax)) + xlim(c(mh.xmin,mh.xmax))
mh.gg <- mh.gg + xlab("Micro-Homology") + theme(legend.title = element_blank(),legend.key.height=unit(2,"line"))


pdf(paste(outputstub,"_baitlen.pdf",sep=""),width=10)
print(bait.gg)
dev.off()
pdf(paste(outputstub,"_mh.pdf",sep=""),width=10)
print(mh.gg)
dev.off()
