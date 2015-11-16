#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(magrittr))


parser <- arg_parser("detect translocation hotspots using scan statistics",
                     name="TranslocHotSpots.R") %>%
    add_argument("tlxfile",
                 "input tlx file of translocation junctions",
                 type="character") %>%
    add_argument("output",
                 "output prefix - several files will be written using this",
                 type="character") %>%
    add_argument("--assembly",
                 "genome assembly to use",
                 default="hg19",
                 type="character") %>%
    add_argument("--bin.width",
                 "initial bin width in basepairs",
                 default=500,
                 type="integer") %>%
    add_argument("--bg.width",
                 "background width in basepairs for caculating enrichment",
                 default=500000,
                 type="integer") %>%
    add_argument("--alpha",
                 "significance threshold",
                 default=0.01,
                 type="numeric") %>%
    add_argument("--hits.min",
                 "minimum number of hits to be called a hotspot",
                 default=3,
                 type="integer") %>%
    add_argument("--strand.min",
                 "minimum number of hits per strand to be called a hotspot",
                 default=1,
                 type="integer") %>%
    add_argument("--cores",
                 "number of compute cores to run on",
                 default=1,
                 type="integer")

if (interactive()) {
    argv <- parse_args(parser,
                       argv=c("~/AltLab/Data/Alt055/results-new/RF204_Alt055/RF204_Alt055_result.tlx",
                              "~/AltLab/Data/Alt055/results-new/RF204_Alt055/RF204_Alt055_hotspots",
                              "--cores", 4,
                              "--assembly", "hg19"))

    source("./TranslocHelper.R")



} else {
    argv <- parse_args(parser)

    source_local <- function(fname){
        argv <- commandArgs(trailingOnly = FALSE)
        base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
        source(paste(base_dir, fname, sep="/"))
    }

    source_local("TranslocHelper.R")

}

if (argv$cores > 1) {
    do_parallel <- TRUE
    suppressPackageStartupMessages(library(doMC, quietly=TRUE))
    registerDoMC(cores=argv$cores)
} else {
    do_parallel <- FALSE
}

suppressPackageStartupMessages(library(readr, quietly=TRUE))
suppressPackageStartupMessages(library(stringr, quietly=TRUE))
suppressPackageStartupMessages(library(plyr, quietly=TRUE))
suppressPackageStartupMessages(library(dplyr, quietly=TRUE))
suppressPackageStartupMessages(library(BSgenome, quietly=TRUE))

tlx <- read_tsv(argv$tlxfile,
                col_types=cols_only(Qname = 'c',
                                    Rname = 'c',
                                    Junction = 'd',
                                    Strand = 'd')) %>%
    arrange(Rname, Junction)


if (str_detect(argv$assembly, "^(mm9|mm10|hg19|hg38)$")) {
    seq_info <- Seqinfo(genome=argv$assembly)
    seq_info <- seq_info[str_subset(names(seq_info),"^chr([0-9]+|[XYM])$")]
    chr.stats <- data_frame(Chr = names(seq_info), Len=seqlengths(seq_info))
} else {
    if (is.na(argv$chrlens)) {
        argv$chrlens <- file.path(Sys.getenv("GENOME_DB"),
                                  argv$assembly,
                                  "annotation",
                                  "ChromInfo.txt")
    }
    if (file.exists(argv$chrlens)) {
        chr.stats <- read_tsv(chr.lens.file, col_names=c("Chr", "Len"))
    } else {
    stop("Error: cannot get chromosome length info")
    }
}

tlx.gr <- tlx %>% mutate(Chr = Rname,
                         Start = Junction,
                         End = Junction,
                         StrandInt = Strand,
                         Strand = ifelse(StrandInt==1, "+", "-")) %>%
    makeGRangesFromDataFrame(keep.extra.columns=T,
                             seqinfo=chr.stats$Len %>% set_names(chr.stats$Chr))



# Create Initial GR Bins

# bins1 <- suppressWarnings(
#     resize(tlx.gr[distance(tlx.gr,tlx.gr[nearest(tlx.gr)]) <= argv$bin.width/2],
#            width=argv$bin.width, fix="center"))

tlx.dis2nearest <- mcols(distanceToNearest(tlx.gr, ignore.strand=T))$distance

bins1 <- suppressWarnings(tlx.gr[tlx.dis2nearest <= argv$bin.width/2] %>%
                              resize(width=argv$bin.width, fix="center") %>%
                              trim())


strand(bins1) <- "*"
bins1 <- unique(bins1)

bins1.by.chr <- split(bins1, seqnames(bins1)) %>% as.list()
tlx.gr.by.chr <- split(tlx.gr, seqnames(tlx.gr)) %>% as.list()

# bins1.by.chr <- split.by.chr(bins1)
# tlx.gr.by.chr <- split.by.chr(tlx.gr)


tlx.cov <- coverage(tlx.gr)
tlx.cov.cumsum <- llply(tlx.cov, cumsum)

if (any(names(bins1.by.chr) != names(tlx.cov))) stop("Error: chromosome names out of order")


bins1.p.local <- mlply(cbind(bins=bins1.by.chr, tlx.cumsum=tlx.cov.cumsum),
                 calculate_local_significance,
                 bin.width=argv$bin.width, bg.width=argv$bg.width,
                 .parallel=do_parallel, .expand=F) %>%
    set_names(names(bins1.by.chr))

bins1.p.chr <- mlply(cbind(bins=bins1.by.chr, tlx.cumsum=tlx.cov.cumsum),
                 calculate_chr_significance,
                 bin.width=argv$bin.width,
                 .parallel=do_parallel, .expand=F) %>%
    set_names(names(bins1.by.chr))





# bins1.by.chr <- GRangesList(mcmapply(calculate.local.significance,bins=bins1.by.chr,tlx.cumsum=tlx.cov.cumsum,bin.width=bin.width,bg.width=bg.width,mc.cores=cores))
# bins1.by.chr <- GRangesList(mcmapply(calculate.chr.significance,bins=bins1.by.chr,tlx.cumsum=tlx.cov.cumsum,bin.width=bin.width,mc.cores=cores))





# Trim peaks to last translocation on each edge

peaks_local_list <- mlply(cbind(b=bins1.by.chr, c=tlx.gr.by.chr, p=bins1.p.local),
                     function(b, c, p) {
                         trim_peaks(reduce(b[p < argv$alpha]), c)
                     }, .parallel=do_parallel, .expand=F) %>%
    set_names(names(bins1.by.chr))

peaks_chr_list <- mlply(cbind(b=bins1.by.chr, c=tlx.gr.by.chr, p=bins1.p.chr),
                     function(b, c, p) {
                         trim_peaks(reduce(b[p < argv$alpha]), c)
                     }, .parallel=do_parallel, .expand=F) %>%
    set_names(names(bins1.by.chr))
#
# peaks.local <- GRangesList(mcmapply(function(b,c) {
#   trimPeaks(reduce(b[b$p.local < alpha]),c)
# }, bins1.by.chr, tlx.gr.by.chr,mc.cores=cores))
#
# peaks.chr <- GRangesList(mcmapply(function(b,c) {
#   trimPeaks(reduce(b[b$p.chr < alpha]),c)
# }, bins1.by.chr, tlx.gr.by.chr,mc.cores=cores))




p.local <- mlply(cbind(bins=peaks_local_list, tlx.cumsum=tlx.cov.cumsum),
                       calculate_local_significance,
                       bin.width=argv$bin.width, bg.width=argv$bg.width,
                       .parallel=do_parallel, .expand=F) %>%
    unlist()

p.chr <- mlply(cbind(bins=peaks_chr_list, tlx.cumsum=tlx.cov.cumsum),
                     calculate_chr_significance,
                     bin.width=argv$bin.width,
                     .parallel=do_parallel, .expand=F) %>%
    unlist()

#
# p.local <- GRangesList(mcmapply(calculate.local.significance,bins=peaks.local,tlx.cumsum=tlx.cov.cumsum,bin.width=bin.width,bg.width=bg.width,mc.cores=cores))
# p.chr <- GRangesList(mcmapply(calculate.chr.significance,bins=peaks.chr,tlx.cumsum=tlx.cov.cumsum,bin.width=bin.width,mc.cores=cores))

peaks_local <- unlist(GRangesList(peaks_local_list))
peaks_chr <- unlist(GRangesList(peaks_chr_list))


peaks_local$q.value <- p.adjust(p.local, method="BY")
peaks_chr$q.value <- p.adjust(p.chr, method="BY")

peaks_local <- peaks_local[peaks_local$q.value <= argv$alpha]
peaks_chr <- peaks_chr[peaks_chr$q.value <= argv$alpha]

peaks_local$hits <- countOverlaps(peaks_local, tlx.gr)
peaks_local$hits.plus <- countOverlaps(peaks_local, tlx.gr[strand(tlx.gr) == "+"])
peaks_local$hits.minus <- countOverlaps(peaks_local, tlx.gr[strand(tlx.gr) == "-"])

peaks_chr$hits <- countOverlaps(peaks_chr,tlx.gr)
peaks_chr$hits.plus <- countOverlaps(peaks_chr, tlx.gr[strand(tlx.gr) == "+"])
peaks_chr$hits.minus <- countOverlaps(peaks_chr, tlx.gr[strand(tlx.gr) == "-"])


peaks_local <- peaks_local[ peaks_local$hits >= argv$hits.min &
                              peaks_local$hits.plus >= argv$strand.min &
                              peaks_local$hits.minus >= argv$strand.min ]
peaks_chr <- peaks_chr[ peaks_chr$hits >= argv$hits.min &
                          peaks_chr$hits.plus >= argv$strand.min &
                          peaks_chr$hits.minus >= argv$strand.min ]

# peaks.local$p.local <- log10(peaks.local$p.local)
# peaks.chr$p.chr <- log10(peaks.chr$p.chr)

peaks_local$q.value <- log10(peaks_local$q.value)
peaks_chr$q.value <- log10(peaks_chr$q.value)

export.bed(peaks_local, str_c(argv$output, "_localpeaks.bed"),
           ignore.strand=T)
export.bed(peaks_chr, str_c(argv$output, "_chrpeaks.bed"),
           ignore.strand=T)

as.data.frame(peaks_local, row.names=NULL) %>%
    write_tsv(str_c(argv$output, "_localpeaks.txt"))

as.data.frame(peaks_chr, row.names=NULL) %>%
    write_tsv(str_c(argv$output, "_chrpeaks.txt"))

