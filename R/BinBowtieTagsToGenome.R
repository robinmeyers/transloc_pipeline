ARGS <- c(
  "tagfile", "character", "file path of bowtie ref tags",
  "output","character", "file path to write output to",
  "nonmapfile","character", "file path of bed file designating unmappable regions in genome",
  "chrlen","integer","Length of chromosome to be binned"
)

OPTS <- c(
  "binsize","integer",500,"Binsize to divy chromosome into (in kb)"
)

TOOLSDIR <- Sys.getenv('TOOLS')
GROSEQ <- Sys.getenv('GROSEQ')
source(paste(TOOLSDIR,"Rsub.R",sep="/"))

parseArgs("BinBowtieTags.R", ARGS, OPTS)
options(scipen=10)


tags <- read.delim(tagfile, sep="\t", header=FALSE, stringsAsFactors=F)
colnames(tags) <- c("Tname","Strand","Tstart")



tags$Tstart <- tags$Tstart + 1

nonmap <- read.delim.save(nonmapfile, header=F, as.is=T, sep="\t")
colnames(nonmap) <- c("Tname","Tstart","Tend")


binsize <- binsize * 1000
x1 <- seq(1,chrlen,binsize)
x2 <- x1 + binsize - 1
x2[length(x2)] <- chrlen


map <- x2-x1+1 - apply(as.array(cbind(x1,x2)),1,bedIntersect,nonmap)

binnedtags <- apply(as.array(cbind(x1,x2)),1,bedCount,tags)

result <- cbind(Tname=nonmap$Tname[1],Tstart=x1,Tend=x2,Map=map,Tags=binnedtags)

write.table(result,output,sep="\t",quote=F,row.names=F)

