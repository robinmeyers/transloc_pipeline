#!/usr/bin/env Rscript

library(argparser)
library(magrittr)

parser <- arg_parser("mark duplicate juntions",
                     name="TranslocDedup.R") %>%
    add_argument("tlxfile",
                 "",
                 type="character") %>%
    add_argument("output",
                 "",
                 type="character") %>%
    add_argument("--offset.dist",
                 "",
                 default=0,
                 type="integer") %>%
    add_argument("--break.dist",
                 "",
                 default=0,
                 type="integer") %>%
    add_argument("--barcode",
                 "",
                 default=NA,
                 type="integer") %>%
    add_argument("--cores",
                 "",
                 default=1,
                 type="integer")

if (interactive()) {
    argv <- parse_args(parser,argv=c("~/AltLab/Data/Alt055/preresults-sample/RF204_Alt055/RF204_Alt055.tlx",
                                     "~/AltLab/Data/Alt055/preresults-sample/RF204_Alt055/RF204_Alt055_duplicate1.txt",
                                     "--cores", 4))
} else {
    argv <- parse_args(parser)
}

library(readr)
library(plyr)
library(dplyr)

l_ply(names(argv), function(i) {
  if (i != "help" && i != "opts") {
      assign(i, argv[[i]], envir=.GlobalEnv)
  }
})


if (break.dist > 0 || offset.dist > 0) {
    stop("Error: still working on implementation for non-zero break.dist and offset.dist")
}


if (cores > 1) {
    library(doMC)
    registerDoMC(cores=cores)
    doParallel <- TRUE
} else {
    doParallel <- FALSE
}

con  <- file(tlxfile, open = "r")
tlxheader <- unlist(strsplit(readLines(con, n = 1),"\t"))
close(con)

col_collectors <- rep(list(col_skip()), length(tlxheader)) %>%
    set_names(tlxheader)

col_collectors[["Qname"]] <-
    col_collectors[["Rname"]] <-
    col_collectors[["Barcode"]] <- col_character()
col_collectors[["Junction"]] <-
    col_collectors[["Strand"]] <-
    col_collectors[["Qstart"]] <-
    col_collectors[["Qend"]] <-
    col_collectors[["B_Strand"]] <-
    col_collectors[["B_Rstart"]] <-
    col_collectors[["B_Rend"]]  <- col_integer()

tlxs <- read_tsv(tlxfile, col_types=col_collectors, progress=F)

if (nrow(tlxs) == 0) {
    cat("Qname\tDups",file=output)
    quit()
}

tlx.juncs <- filter(tlxs, !is.na(Junction) & Rname != "Adapter") %>%
    arrange(Qname) %>%
    mutate(Offset = ifelse(Strand==1, Junction - Qstart, Junction + Qstart),
           B_Junction = ifelse(B_Strand==1, B_Rend, B_Rstart))

dup.rows <- tlx.juncs %>% group_indices(Rname, Strand, Offset, B_Junction)

out.df <- split(tlx.juncs, dup.rows) %>% ldply(function(d) {
    if (nrow(d) > 1) {
        return(as.data.frame(cbind(d$Qname[2:nrow(d)] ,d$Qname[1])))
    } else {
        return(NULL)
    }
}, .id=NULL, .parallel=doParallel)

write.table(out.df, output, sep="\t", quote=F, row.names=F, col.names=F)

quit()

# if (commandArgs()[1] != "RStudio") {

#   ARGS <- c(
#     "tlxfile", "character", " ",
#     "output","character", " "
#   )

#   OPTS <- c(
#     "offset.dist","numeric",0," ",
#     "break.dist","numeric",0," ",
#     "random.barcode","numeric",0," ",
#     "cores","numeric",1,"Number of compute nodes to run on"
#   )


#   source_local <- function(fname){
#     argv <- commandArgs(trailingOnly = FALSE)
#     base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
#     source(paste(base_dir, fname, sep="/"))
#   }

#   source_local("Rsub.R")
#   parseArgs("TranslocDedup.R", ARGS, OPTS)

# } else {
#   source("~/TranslocPipeline//R/Rsub.R")
#   tlxfile <- "/Volumes//AltLab/Translocation/NewPipeline/Alt046/results/PW002_Alt046/PW002_Alt046.tlx"
#   output <- "~/Working/NewPipelineValidations/DedupTesting/PW002_dedup.txt"

#   cores <- 4
# }

suppressPackageStartupMessages(library(plyr, quietly=TRUE))
suppressPackageStartupMessages(library(parallel, quietly=TRUE))


con  <- file(tlxfile, open = "r")
header <- unlist(strsplit(readLines(con, n = 1),"\t"))
close(con)

headersToSkip <- c("Seq")

colClasses <- rep(NA,length(header))
colClasses[match(headersToSkip,header)] <- "NULL"

tlxs <- read.delim(tlxfile,header=T,colClasses=colClasses,as.is=T)


tlxs <- tlxs[with(tlxs,order(Rname,Strand,Junction)),]

tlxs$Offset <- with(tlxs,ifelse(Strand==1,Junction-Qstart,Junction+Qstart))
tlxs$B_Junction <- with(tlxs,ifelse(B_Strand==1,B_Rend,B_Rstart))
# tlxs$B_Offset <- with(tlxs,ifelse(B_Strand==1,B_Junction-B_Qend,B_Junction+B_Qend))

if (nrow(tlxs) > 0) {
  tlxs_by_chr_and_strand <- split(tlxs,list(tlxs$Rname,tlxs$Strand))


  findDuplicates <- function(n,tlxs) {
    tlx <- tlxs[n,]
    matches <- subset(tlxs,
                        Qname>tlx$Qname &
                        abs(Offset-tlx$Offset) <= offset.dist &
                        # abs(Junction-tlx$Junction)<=junc_dist &
                        abs(B_Junction-tlx$B_Junction) <= break.dist &
                        abs(B_Qend-tlx$B_Qend) <= break.dist)

    if (random.barcode > 0 && nchar(tlx$Barcode) > 0) {
      matches <- subset(matches, adist(Barcode,tlx$Barcode) <= 2)
    }

    if (nrow(matches) > 0) {
      return(paste(paste(matches$Qname,"(",matches$B_Junction-tlx$B_Junction,",",matches$Junction-tlx$Junction,")",sep="")[1:min(nrow(matches),3)],collapse=","))
    } else {
      return("")
    }
  }

  if (cores == 0) {
    cores <- detectCores()
  }



  cat("Deduplicating junctions on",cores,"cores\n")
  print(object.size(tlxs_by_chr_and_strand),units="Mb")
  tlxs <- ldply(lapply(1:length(tlxs_by_chr_and_strand),function (n,tlxs) {
    if (nrow(tlxs[[n]]) > 0) {

      cat(nrow(tlxs[[n]])," - ")
      print(object.size(tlxs[[n]]),units="Mb")

      dups <- mclapply(1:nrow(tlxs[[n]]),findDuplicates,tlxs[[n]][,c("Qname","Offset","Junction","B_Junction","B_Qend","Barcode")],mc.cores=cores)

      print(object.size(dups),units="Mb")

      tlxs[[n]]$Dups <- unlist(dups)
      return(tlxs[[n]])
    }
  },tlxs_by_chr_and_strand))

write.table(tlxs[tlxs$Dups != "",c("Qname","Dups")],output,sep="\t",quote=F,na="",row.names=F,col.names=F)

} else {
  cat("Qname\tDups",file=output)
}
