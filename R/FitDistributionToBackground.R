#!/usr/bin/Rscript --vanilla

ARGS <- c(
  "bkgdfile", "character", "file path of groseq background bins",
  "outdir","character", "file path of reference data (TSSs) serving as the reference for the distances",
  "win","integer", "Bin size used for binning alignments across genome (in kb)"
)

OPTS <- c(
  "background","numeric",1,"Calculate poisson using this assumed percentage of reads from background",
  "alpha","numeric",2,"Size/shape parameter in negative binomial distribution",
  "xmax","numeric",1,"X-axis limit for histogram",
  "ymax","numeric",0,"Y-axis limit for histogram; 0 - autoscale",
  "binnum","integer",200,"Number of bins for histogram",
  "expt","character","","Experiment title for labeling the histogram"
)

TOOLSDIR <- "~/Tools"
TLPDIR <- "~/TLPipeline"
source(paste(TOOLSDIR,"Rsub.R",sep="/"))

parseArgs("FitPoissonToBKGD.R", ARGS, OPTS)
options(scipen=0)


data <- read.delim(bkgdfile, sep="\t",header=TRUE, stringsAsFactors=F)

data$Density <- data$Tags/(data$Map)*1000
totalReads <- sum(as.numeric(data$Tags))
mappableBP <- sum(as.numeric(data$Map))
genomebins <- nrow(data[data$Map > 0,])
data <- data[data$Map > 0 & data$Density <= xmax & data$Density > 0,]
#genomebins <- nrow(data)


bins <- seq(0,xmax,xmax/binnum)
#data <- data[data$Density > bins[2],]



background <- background/100
lambda <- background*totalReads/(mappableBP)*1000
#poisson <- exp(bins*(window)*log(lambda) - lambda - lfactorial(bins*(window)))
x <- 0:(max(bins)*win)
#poisson <- dpois(x,lambda)
nbinom <- dnbinom(x,size=alpha,prob=alpha/(alpha+lambda*win))
#nbinom <- nbinom*genomebins/length(bins)*length(x)

#smalldata <- data[data$Density < .004,]

cat("Expt:",expt,"\n")
cat("Total Mappable Bases (both strands):",mappableBP,"\n")
cat("Genome Bins:",genomebins,"\n")
cat("Total Hits:",totalReads,"\n")
cat("Number of Reads from Background:",background*totalReads,"\n")
cat("Background Tag Density (per kb):",lambda,"\n")
cat("Length bins:",length(bins),"Length x:", length(x), "\n")
#cat("Mean:",mean(smalldata$Density*window),"Var:",var(smalldata$Density*window),"\n")


pdf(paste(outdir,"/background-poisson.pdf",sep=""))
  if(ymax==0) {ylimopt <- NULL} else { ylimopt <- c(0,ymax)}
  denshist <- hist(data$Density,breaks=bins,ylim=ylimopt,col="grey",xlab="GRO-seq Density (reads/kb)",ylab=paste("Number of",win,"kb Windows",sep=" "),main=paste(expt,"GRO-seq Background Calculation",sep=" "))

  #lines(x/window,nbinom,col="red",lwd=3)
  par(new=T)
  plot(x,nbinom,axes=F,type='l',col="red",xlab="",ylab="",lwd=3)
  

  legend("topright",legend=paste("background = ",background*100,"%\nlambda = ",signif(lambda,6)," tags/kb\nalpha = ",alpha,sep=""),lwd=3,col="red",bty="n")
dev.off()



