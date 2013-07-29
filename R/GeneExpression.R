#!/usr/bin/Rscript

ARGS <- c(
  "tagfile", "character", "file path of bowtie ref tags",
  "output","character", "file path to write summary output to",
  "refgenefile","character","file path of refseq gene list",
  "nonmapfile","character", "file path of bed file designating unmappable regions in genome"
)

OPTS <- c(
  "tagsize","numeric",35,"Length of tags in basepairs",
  "lambda","numeric",0.1,"Mean background tags per kb (both strands) from negative binomial distribution",
  "alpha","numeric",2,"Shape/size variable from negative binomial distribution",
  "profile_file","character",NA,"File stub to write profile to",
  "upstream","numeric",5,"Kilobases upstream from TSS to begin profile",
  "downstream","numeric",10,"Kilobases downstream from TSS to begin profile",
  "binsize","numeric",50,"Bin size in bp for profile",
  "nodes","numeric",8,"Number of compute nodes to run on"
)


source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local("Rsub.R")
source_local("calcGeneExpr.R")
require(parallel,quietly=T,warn.conflicts=F)
require(plyr,quietly=T,warn.conflicts=F)

parseArgs("GeneExpression.R", ARGS, OPTS)
options(scipen=5)

tags <- read.delim(tagfile, sep="\t", header=FALSE, stringsAsFactors=F)

colnames(tags) <- c("Tname","Strand","Tstart")
# bowtie ref tag coordinates were 0 based; change to 1
tags$Tstart <- tags$Tstart + 1
#tags <- tags[with(tags,order(Strand,Tstart)),]


refgene <- read.delim.save(refgenefile, header=T, as.is=T, sep="\t")
colnames(refgene) <- c("Gene","Tname","Tstart","Tend","Length","RefSeq","Strand")
refgene$Tname <- paste("chr",refgene$Tname,sep="")
refgene <- refgene[refgene$Tname == tags$Tname[1] & refgene$Length >= 3000,]
refgene <- refgene[with(refgene,order(Tstart,Tend)),]

nonmap <- read.delim.save(nonmapfile, header=F, as.is=T, sep="\t")
colnames(nonmap) <- c("Tname","Tstart","Tend")

#convert lambda into tags/bp
lambda <- lambda/1000

cat("\nCalculating activity for",nrow(refgene),"genes on",tags$Tname[1],"with",nodes,"compute nodes\n")
# cl <- makeSOCKcluster(rep("localhost",nodes))
# clusterExport(cl,c("refgene","tags","nonmap","tagsize","lambda","alpha","bedCount","bedIntersect"))

ngenes <- nrow(refgene)

summary <- ldply(mclapply(1:ngenes, calcGeneExpr, mc.cores=nodes))

write.table(summary,output,sep="\t",quote=F,row.names=F)


if (profile_file == "NA") {
  profile_sense_file <- NA
}

if (! is.na(profile_file) ) {
  


  x1 <- seq(-upstream*1000,downstream*1000,binsize)
  x2 <- x1 + binsize - 1
  bins <- data.frame(x1=x1,x2=x2)

#   clusterExport(cl, c("summary","bins","binsize"))


  profile <- ldply(mclapply(1:ngenes, profileGene,mc.cores=nodes))

  write.table(profile,profile_file,sep="\t",quote=F,row.names=F)

}

# stopCluster(cl)



