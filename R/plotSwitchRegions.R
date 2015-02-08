#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {
  
  ARGS <- c(
    "tlxfile","character","file path of tlx file",
    "output","character","file path ofoutput plot"
    
  )
  
  OPTS <- c(
    "facetscales","character","free","allow free or free_x",
    "numbins","numeric",100,"",
    "flanks","numeric",0,"",
    "force.regions","character","",""
    
  )
  
  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }
  
  source_local("Rsub.R")
  
  parseArgs("removeDupClones.R", ARGS, OPTS)
  
} else {
  
  tlxfile <- "~/Working/NewPipelineValidations/Junchao/IsceI-Switch_joins/JD003_dedup1.tlx"
  output <- sub(".tlx","_switch.pdf",tlxfile)
  numbins <- 100
  flanks <- 5000
  facetscales <- "free"
  source("~/TranslocPipeline//R/Rsub.R")
  force.regions <- ""

}


suppressPackageStartupMessages(library(ggplot2,quietly=TRUE))
suppressPackageStartupMessages(library(gridExtra,quietly=TRUE))


chr <- "chr12"
switch_regions <- data.frame(Name=c("Smu","Sg3","Sg1","Sg2b","Sg2a","Se","Sa"),
                             Start=c(114661032,114599706,114568854,114546284,114527282,114512435,114498531),
                             End=c(114664779,114603915,114577055,114552458,114533996,114515374,114502837))

  

cores <- 8


make_bins <- function(start,end,n) {
  starts <- seq(start,end,length.out=n)
  ends <- starts + diff(starts)[1]
  strands <- rep(c(1,-1),each=length(starts))
  
  bins <- data.frame(Start=rep(starts,2),End=rep(ends,2),Strand=strands)
  bins <- within(bins, Mid <- (Start + End) / 2)
  return(bins)
}

data <- read.delim(tlxfile,header=T,as.is=T)

xmin <- min(switch_regions$Start) - 10000
xmax <- max(switch_regions$End) + 10000
bins <- make_bins(xmin,xmax,numbins)

bins$TLX <- unlist(lapply(1:nrow(bins),function(i,bins,data) {
  return(nrow(subset(data,Rname == chr & Strand == bins$Strand[i] &
                       Junction >= bins$Start[i] & Junction < bins$End[i])))
},bins,data))

bins$TLX <- bins$TLX * bins$Strand

switch_regions_data <- rbind(switch_regions,switch_regions)
switch_regions_data <- within(switch_regions_data, {
  Strand <- rep(c(1,-1),each=nrow(switch_regions))
  Mid <- (Start + End) / 2
})

# switch_regions_data$y <- unlist(mclapply(1:nrow(switch_regions_data),function(i,swtch,bins) {
#   if(swtch$Strand[i]==1){
#     return(max(bins$TLX[bins$End>swtch$Start[i] & bins$Start<swtch$End[i]]))
#   } else {
#     return(min(bins$TLX[bins$End>swtch$Start[i] & bins$Start<swtch$End[i]]))
#   }
# },switch_regions_data,bins,mc.cores=cores))
# 

  

switch_regions_data$TLX <- unlist(lapply(1:nrow(switch_regions_data),function(i,bins,data) {
  return(nrow(subset(data,Rname == chr & Strand == bins$Strand[i] &
                       Junction >= bins$Start[i] & Junction < bins$End[i])))
},switch_regions_data,data))

# switch_regions_data$y <- switch_regions_data$TLX * switch_regions_data$Strand
switch_regions_data$y <- ifelse(switch_regions_data$Strand==1,max(bins$TLX),min(bins$TLX))
switch_regions_data$Strand <- factor(switch_regions_data$Strand)
bins$Strand <- factor(bins$Strand)

# totals <- data.frame(variable="TLX",Strand=as.factor(c(1,-1)),x=min(newdata$Mid),y=c(max(newdata$TLX),min(newdata$TLX)),lab=c(sum(newdata$TLX[newdata$Strand == 1]),-sum(newdata$TLX[newdata$Strand == -1])))
# 
# 
# master <- melt(newdata,id.vars=c("Start","End","Mid","Strand"))
# master$Strand <- as.factor(master$Strand)

p_main <- ggplot() +
  geom_rect(aes(xmin=Start,xmax=End,ymin=-Inf,ymax=Inf),data=switch_regions,fill="black",alpha=0.15) +
  geom_ribbon(aes(x=Mid,ymin=0,ymax=TLX,fill=Strand),data=bins,color="black",alpha=0.5,size=0.2,show_guide = FALSE) +
  scale_x_continuous(trans="reverse") +
  scale_y_continuous(trans="reverse") +
  ylab("") +
  xlab("Position") +
  ggtitle(sub(".tlx","",basename(tlxfile))) +
#   geom_text(aes(x=Start,y=y,label=TLX),data=switch_regions_data,colour="black",vjust=0,hjust=0,size=2,show_guide = FALSE) +
  geom_text(aes(x=End,label=Name),data=switch_regions,y=-min(bins$TLX),colour="black",hjust=1,vjust=1,size=3,show_guide = FALSE) +
  theme(axis.text.x = element_text(size=8),axis.title.x = element_blank())



switch_agg <- aggregate(TLX ~ Name, data=switch_regions_data,sum)
switch_agg <- with(switch_agg,switch_agg[order(-TLX),])

force.regions <- unlist(strsplit(force.regions,","))
which.regions <- switch_agg[switch_agg$Name %in% force.regions,]

for (i in 1:nrow(switch_agg)) {
  if (nrow(which.regions) < 4 && !switch_agg$Name[i] %in% which.regions$Name) {
    which.regions <- rbind(which.regions,switch_agg[i,])
  }
}


master <- data.frame()
totals <- data.frame()

for (srname in rev(switch_regions$Name[switch_regions$Name %in% which.regions$Name])) {
  sr <- switch_regions[switch_regions$Name == srname,]
  subbins <- make_bins(sr$Start-flanks,sr$End+flanks,numbins/2)
  subbins$Name <- srname
  subbins$TLX <- unlist(lapply(1:nrow(subbins),function(i,bins,data) {
    return(nrow(subset(data,Rname == chr & Strand == bins$Strand[i] &
                         Junction >= bins$Start[i] & Junction < bins$End[i])))
  },subbins,data))
  
  subbins$TLX <- subbins$TLX * subbins$Strand
  subbins$Strand <- factor(subbins$Strand)
  
  tally <- data.frame(Name=srname,variable="TLX",Strand=as.factor(c(1,-1)),x=min(subbins$Mid),y=c(max(subbins$TLX),min(subbins$TLX)),lab=c(sum(subbins$TLX[subbins$Strand == 1]),-sum(subbins$TLX[subbins$Strand == -1])))
  
  master <- rbind(master,subbins)
  totals <- rbind(totals,tally)
}

master$Name <- factor(master$Name,ordered=T,levels=switch_regions$Name)

p_sub <- ggplot() +
  facet_wrap(~ Name, ncol=2, scales=facetscales) +
  geom_ribbon(aes(x=Mid,ymin=0,ymax=TLX,fill=Strand),data=master,color="black",alpha=0.5,size=0.2,show_guide = FALSE) +
  geom_rect(aes(xmin=Start,xmax=End,ymin=-Inf,ymax=Inf),data=switch_regions[switch_regions$Name %in% totals$Name,],fill="black",alpha=0.15) +
    
  scale_x_reverse() +
  scale_y_reverse() +
  ylab("") +
  xlab("Position") +
  geom_text(aes(x=x,y=y,label=lab,ymax=NULL,color=Strand),data=totals,size=3,hjust=1,show_guide = FALSE) +
  theme(axis.text.x = element_text(size=8),axis.title.x = element_blank())


pdf(output)

grid.arrange(p_main,p_sub,nrow=2,heights=c(2,3))


dev.off()
