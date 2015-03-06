#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlx.files", "character", "comma-separated list of files or dir and will grab all *.tlx",
    "ref.file","character",""
  )
  
  OPTS <- c(

    "outstub","character","./","file path to plot to",
    "ref.format","character","auto","format of reference file - 'refseq' or 'bed'",
    "TSS","logical",T,"use only start of feature",
    "tlx.labels","character","",""
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
  tlx.files <- "~/AltLab/Simulations/mapqual/results/tlx"
  outstub <- "~/AltLab/Simulations/mapqual/results/dis2TSS/results"
  ref.file <- "/Volumes/AltLab/Genomes/mm9/annotation/refGene.txt"
  ref.format <- "auto"
  tlx.labels <- ""
  TSS <- T
}

suppressPackageStartupMessages(library(data.table, quietly=TRUE))
suppressPackageStartupMessages(library(dplyr, quietly=TRUE))
suppressPackageStartupMessages(library(GenomicRanges, quietly=TRUE))

tss.outfile <- paste(outstub,"_dis2TSS.txt",sep="")
tss.pdffile <- paste(outstub,"_dis2TSS.pdf",sep="")


ref.data <- fread(ref.file,header=F,sep="\t")

if (ref.format == "auto") {
  if (is.character(ref.data$V1)) {
    ref.format <- "bed"  
  } else {
    ref.format <- "refseq"
  }
}

if (ref.format == "refseq") {
  ref.gr <- GRanges(seqnames = ref.data$V3,
                    ranges = IRanges(start=ref.data$V5,
                                     end=ref.data$V6,
                                     names = ref.data$V13),
                    strand = ref.data$V4)
} else if (ref.format == "bed") {
  ref.gr <- GRanges(seqnames = ref.data$V1,
                    ranges = IRanges(start=ref.data$V2,
                                     end=ref.data$V3,
                                     names = ref.data$V4),
                    strand = ref.data$V6)
} else {
  stop("Error: do not recognize ref.format")
}

if (TSS == T) {
  start(ref.gr) <- ifelse(strand(ref.gr) == "+",start(ref.gr),end(ref.gr))
  end(ref.gr) <- ifelse(strand(ref.gr) == "+",start(ref.gr),end(ref.gr))
}


if (file.info(tlx.files)[["isdir"]]) {
  tlx.files <- list.files(tlx.files,pattern = "*.tlx",full.names = T)
} else {
  tlx.files <- unlist(strsplit(tlx.files,","))
}

names(tlx.files) <- sub(".tlx","",basename(tlx.files))

if (tlx.labels != "") {
  tlx.labels <- unlist(strsplit(tlx.labels,","))
  if (length(tlx.labels) != length(tlx.files)) {
    stop("Error: label list has different length than tlxfile list")
  }
  names(tlx.files) <- tlx.labels
}

if (!file.exists(dirname(outstub))) {
  dir.create(dirname(outstub),recursive=T)
}

tlx.dis2TSS <- list()

for (tlx.file in names(tlx.files)) {
  
  tlx <- fread(tlx.files[tlx.file],sep="\t",header=T)
  
  tlx.gr <- GRanges(seqnames = tlx$Rname,
                  ranges = IRanges(start=tlx$Junction,
                                   end=tlx$Junction,
                                   names = tlx$Qname),
                  strand = tlx$Strand)

#   suppressWarnings(
#     hits <- as.data.frame(distanceToNearest(tlx.gr,ref.gr))
#   )
  suppressWarnings(
    tlx <- tlx %>% mutate(Nearest.IDX = nearest(tlx.gr,ref.gr)) %>%
      filter(!is.na(Nearest.IDX))
  )
  

  tlx.dis2TSS[[tlx.file]] <- tlx %>%
    mutate(NearestTSS = names(ref.gr)[Nearest.IDX],
           Dis2TSS = ifelse(strand(ref.gr)[Nearest.IDX] == "+",
                            pmin(Junction - start(ref.gr)[Nearest.IDX],
                                 Junction - end(ref.gr)[Nearest.IDX]),
                            pmin(start(ref.gr)[Nearest.IDX] - Junction,
                                 end(ref.gr)[Nearest.IDX] - Junction)),
           Library = tlx.file) %>%
    select(Library,Qname,Rname,Junction,Strand,NearestTSS,Dis2TSS)
           
}

tlx.dis2TSS <- rbind_all(tlx.dis2TSS)

write.table(tlx.dis2TSS, tss.outfile, quote = F, sep = "\t", row.names = F)


system(paste("dis2TSSPlot.R",tss.outfile))

