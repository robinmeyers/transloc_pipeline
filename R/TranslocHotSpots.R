suppressPackageStartupMessages(library(GenomicRanges, quietly=TRUE))

source("~/TranslocPipeline/R/Rsub.R")
source("~/TranslocPipeline/R/TranslocHelper.R")

calculate.bin.significance <- function(bins,tlx.gr,bg.scale) {
  bg.bins <- resize(bins,width=bg.scale*width(bins),fix="center")
  bins$hits <- countOverlaps(bins,tlx.gr)
  bg.bins$hits <- countOverlaps(bg.bins,tlx.gr)
  bg.bins$prob <- width(bg.bins)/(width(bg.bins)+bg.bins$hits)
  bins$p <- pnbinom(bins$hits-1,size=width(bins),prob=bg.bins$prob,lower.tail=F)
  return(bins$p)
}

bin.width <- 500
bg.scale <- 10

tlxfile <- "/Volumes/AltLab/Translocation/NewPipeline/Alt055/results/RF204_Alt055/RF204_Alt055.tlx"
tlx <- readTLX(tlxfile,columnsToRead=c("Qname","Rname","Junction","Strand"))

chrlen <- getChromLens("hg19")
chr.stats <- data.frame(len=chrlen)

tlx.gr <- tlxToGR(tlx,chrlen)
chr.hit.table <- table(seqnames(tlx.gr))
chr.stats[names(chr.hit.table),"hits"] <- chr.hit.table

# Create Initial GR Bins
bins1 <- resize(tlx.gr,width=bin.width,fix="center")
# start(bins1) <- pmax(1,ifelse(strand(tlx.gr) == "+",start(tlx.gr)-bin.width/2,start(tlx.gr)-bin.width/2+1))
# end(bins1) <- pmin(seqlengths(tlx.gr)[as.character(seqnames(tlx.gr))],ifelse(strand(tlx.gr) == "+",start(tlx.gr)+bin.width/2-1,start(tlx.gr)+bin.width/2))
strand(bins1) <- "*"

# Merge overlapping bins

bins2 <- reduce(bins1)

# For each chromosome, run a merging algorithm

bins3 <- bins2[0,0]

for (chr in seqlevels(bins2)) {
  
  chr.bins <- bins2[seqnames(bins2)==chr]
  
  repeat {
    nbins <- length(chr.bins)
        
    # Find left merge
#     chr.bins.left <- chr.bins
#     start(chr.bins.left)[2:nbins] <- start(chr.bins)[1:(nbins-1)]
#     chr.bins.left$hits <- countOverlaps(chr.bins.left,tlx.gr)
    
    # Perform right merge
    chr.bins.right <- chr.bins
    end(chr.bins.right)[1:(nbins-1)] <- end(chr.bins)[2:nbins]
    chr.bins.right$hits <- countOverlaps(chr.bins.right,tlx.gr)
    
    chr.bins$p <- calculate.bin.significance(chr.bins,tlx.gr,bg.scale)
    chr.bins$p.right <- calculate.bin.significance(chr.bins.right,tlx.gr,bg.scale)
    chr.bins$p.left[1] <- 1
    chr.bins$p.left[2:nbins] <- chr.bins$p.right[1:(nbins-1)]
    
    chr.bins$merge <- 0
    chr.bins$merge <- ifelse(chr.bins$p.left <= chr.bins$p & chr.bins$p.left < chr.bins$p.right, -1, chr.bins$best)
    chr.bins$merge <- ifelsechr.bins$p.right <= chr.bins$p & chr.bins$p.right <= chr.bins$p.left, 1, chr.bins$best)

    chr.bins.new <- chr.bins[0,0]
    for (i in 1:nbins) {
      if (i < nbins &)
      
      c(chr.bins.new,chr.bins[i])
    }
    
    
    
    
    
    if (length(chr.bins.new) == length(chr.bins) ) break
    chr.bins <- chr.bins.new
  }
  
  while (length(chr.bins.new) < chr.bins)
  # merge with left bin
  
  
}





