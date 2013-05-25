# ARGS <- c(
#   "targetsfile","character","file path of target file",
#   "dir","character","file path to directory where data is stored",
#   "outstub","character","file path to write output files to"
# )
# 
# OPTS <- c(
#   "perc_similar","numeric",75,"percent similarity threshold for clones to be repeats"
# )
# 
# source_local <- function(fname){
#   argv <- commandArgs(trailingOnly = FALSE)
#   base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
#   source(paste(base_dir, fname, sep="/"))
# }
# 
# source_local("Rsub.R")
# 
# parseArgs("DGE-GROseq.R", ARGS, OPTS)

library(edgeR)

targets <- readTargets(targetsfile)
setwd(dir)


d <- readDGE(targets)
colnames(d) <- targets$files
d <- d[rowSums(cpm(d) > 5 ) >= 2,]
d <- calcNormFactors(d)
d <- estimateCommonDisp(d,verbose=T)
d <- estimateTagwiseDisp(d)
plotBCV(d)
et <- exactTest(d,pair=c("MEF","ESC"))
de <- decideTestsDGE(et, p=0.000001, adjust.method="BH")
summary(de)
detags <- rownames(d)[as.logical(de)]
pdf(paste(outstub,".pdf",sep=""))
plotSmear(et,de.tags=detags)
abline(h=c(-2,2),col="blue")
dev.off()

write.table(et$table[detags,],paste(outstub,".txt",sep=""),quote=F,sep="\t")


