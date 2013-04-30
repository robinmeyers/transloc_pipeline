ARGS <- c(
  "infile", "character", "file path of unfiltered tlx file",
  "outfile","character","file path to output filtered tlx file",
  "rest1","character", "restriction site expected in forward primer",
  "rest2","character", "restriction site expected in reverse primer"
)

OPTS <- c(
  "maxoverlap","numeric",30,"Maximum overlap between locus sequence and alignment",
  "cleanMarg","numeric",30,"Number of bp from beginning of alignment where rest2 site is not allowed"
)

source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local("Rsub.R")

parseArgs("4CFilter.R", ARGS, OPTS)

con  <- file(infile, open = "r")
header <- unlist(strsplit(readLines(con, n = 1),"\t"))
close(con)

headersToSkip <- c("RedpQstart","RedpQend","BlupQstart","BlupQend","Margin","blockSizes","qStarts","tStarts")

colClasses <- rep(NA,length(header))
colClasses[match(headersToSkip,header)] <- "NULL"

tlx <- read.delim(infile,header=T,as.is=T,colClasses=colClasses)

tlx <- tlx[which(!is.na(tlx$BluQstart)),]
tlx <- tlx[which(!is.na(tlx$RedQend)),]

# Overlap and gap filters
tlx <- tlx[with(tlx, RedQend - Qstart < maxoverlap) ,]
tlx <- tlx[with(tlx, RedQend - Qstart >= nchar(rest1)) ,]
tlx <- tlx[with(tlx, Qend - BluQstart < maxoverlap) ,]


# Check for correct restriction site at end of locus alignments
tlx <- tlx[with(tlx,adist(rest1,substr(Raw,RedQend-5,RedQend))) < 2,]
tlx <- tlx[with(tlx,adist(rest2,substr(Raw,BluQstart,BluQstart+3))) < 2,]

tlx <- tlx[with(tlx,regexpr(rest2,Raw,ignore.case=T) >= Qstart + cleanMarg),]


tlx_filtered <- tlx[0,]

for (qname in unique(tlx$Qname)) {
  reads <- tlx[tlx$Qname == qname,]
  reads <- reads[order(reads$Mismatch-reads$Match),]
  if (nrow(reads) > 1) {
   if (reads$Match[1]-reads$Mismatch[1] - (reads$Match[2]-reads$Mismatch[2]) > 3) {
     tlx_filtered <- rbind(tlx_filtered,reads[1,])
    }
  } else {
    tlx_filtered <- rbind(tlx_filtered,reads[1,])
  }
}

write.table(tlx_filtered,outfile,sep="\t",row.names=F,quote=F)


