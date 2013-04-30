
profileGene <- function (idx) {


  gene <- summary[idx,]
  strand <- gene$Strand

  if (gene$Prom.P < 0.001) {
    if (strand == "+") {
      TSS <- gene$Prom.Tstart
    } else {
      TSS <- gene$Prom.Tend
    }
  } else {
    if (strand == "+") {
      TSS <- gene$Ref.Tstart
    } else {
      TSS <- gene$Ref.Tend
    }
  }

  # Determine crude boundaries
  if (strand == "+") {
    T1 <- TSS + min(bins$x1) - 1000
    T2 <- TSS + max(bins$x2) + 1000
  } else {
    T1 <- TSS - max(bins$x2) - 1000
    T2 <- TSS - min(bins$x1) + 1000
  }

  # Limit tags and map to these boundaries
  gtags <- tags[tags$Tstart >= T1 & tags$Tstart <= T2,]
  gnonmap <- nonmap[nonmap$Tend >= T1 & nonmap$Tstart <= T2,]

  # Shift tags on - strand
  gtags$Tstart <- ifelse(gtags$Strand == "+", gtags$Tstart, gtags$Tstart+tagsize-1)

  # Make coordinates relative to TSS
  if (strand == "+") {
    gtags$Tstart <- gtags$Tstart - TSS
    gnonmap$Tstart <- gnonmap$Tstart - TSS
    gnonmap$Tend <- gnonmap$Tend - TSS
  } else {
    gtags$Tstart <- TSS - gtags$Tstart
    gnonmap$Tstart <- TSS - gnonmap$Tstart
    gnonmap$Tend <- TSS - gnonmap$Tend
  }

  tags_sense <- gtags[gtags$Strand == strand,]
  tags_anti <- gtags[gtags$Strand != strand,]

  gnonmap_p <- gnonmap
  gnonmap_m <- gnonmap
  gnonmap_m$Tstart <- gnonmap_m$Tstart + tagsize - 1
  gnonmap_m$Tend <- gnonmap_m$Tend + tagsize - 1

  bintags_sense <- apply(bins,1,bedCount,tags_sense)
  bintags_anti <- apply(bins,1,bedCount,tags_anti)
  if (gene$Strand == "+") {
    map_sense <- binsize - apply(bins,1,bedIntersect,gnonmap_p)
    map_anti <- binsize - apply(bins,1,bedIntersect,gnonmap_m)
  } else {
    map_sense <- binsize - apply(bins,1,bedIntersect,gnonmap_m)
    map_anti <- binsize - apply(bins,1,bedIntersect,gnonmap_p)
  }

  profile <- data.frame(Gene=gene$Gene,Tname=gene$Tname,Strand=strand,TSS=TSS,
    BinStart=bins$x1,BinEnd=bins$x2,SenseTags=bintags_sense,SenseMap=map_sense,
    AntiTags=bintags_anti,AntiMap=map_anti)



  



}

#profileGene <- function (idx, str, ) {
#
#  if ( ! str %in% c("sense","anti","both") ) {
#    cat("Error: strand must be specified as \"sense\", \"anti\", or \"both\"")
#    return(0)
#  }
#
#  gene <- summary[idx,]
#  strands <- c("+","-")
#  
#  if (str == "sense") {
#    strand <- gene$Strand
#  } else if (str == "anti") {
#    strand <- strands[strands != gene$Strand]
#  } else {
#    strand <- strands
#  }
#
#  TSS <- gene$Ref.Tstart
#
#  gtags <- tags[tags$Strand %in% strand,]
#  gtags$Tstart <- ifelse(gtags$Strand == "-", gtags$Tstart+tagsize-1, gtags$Tstart)
#  if (gene$Strand == "+") {
#    gtags$Tstart <- gtags$Tstart - TSS
#  } else {
#    gtags$Tstart <- TSS - gtags$Tstart
#  }
#  gtags <- gtags[gtags$Tstart >= min(bins$x1) & gtags$Tstart <= max(bins$x2),]
#
#  profile <- matrix(unlist(apply(bins,1,bedCount,gtags)),byrow=T,nrow=1)
#  colnames(profile) <- (bins$x1 + bins$x2) / 2
#  rownames(profile) <- gene$Gene
#
#  return(profile)
#
#}
#


# Calculate gene expression information according to Core, Waterfall, Lis 2008
# Items that must be in scope:
# refgene, tags, nonmap, lambda, alpha, tagsize, bedCount, bedIntersect

calcGeneExpr <- function (idx) {
  # idx is refgene index to be used in calculation

  if (length(idx) != 1) {
    cat("Error: idx in calcGeneExpr must be an integer of length 1")
    return(0)
  }

  gene <- refgene[idx,]
  T1 <- gene$Tstart - 2500
  T2 <- gene$Tend + 2500


  tags_sense <- tags[tags$Strand == gene$Strand & tags$Tstart > T1 & tags$Tstart < T2,]
  tags_anti <- tags[tags$Strand != gene$Strand & tags$Tstart > T1 & tags$Tstart < T2,]
  gnonmap_p <- nonmap[nonmap$Tend > T1 & nonmap$Tstart < T2, ]
  gnonmap_m <- gnonmap_p
  gnonmap_m$Tstart <- gnonmap_m$Tstart + tagsize - 1
  gnonmap_m$Tend <- gnonmap_m$Tend + tagsize - 1

  # Shift tag on (-) strand by the tag length
  if (gene$Strand == "+") {
    tags_anti$Tstart <- tags_anti$Tstart + tagsize - 1
  } else {
    tags_sense$Tstart <- tags_sense$Tstart + tagsize - 1
  }



  ### Step 1: Identification of Promoter Proximal Peak ###

  # Construct promoter proximal windows
  if (gene$Strand == "+") {
    ppstart <- seq(gene$Tstart-1000,gene$Tstart+1000,5)
  } else {
    ppstart <- seq(gene$Tend-1000,gene$Tend+1000,5)
  }
  ppend <- ppstart + 50 - 1
  ppbins <- data.frame(Tstart=ppstart,Tend=ppend)
  
  # Count tags per window and number of mappable bases
  pptags <- apply(ppbins,1,bedCount,tags_sense)

  if (gene$Strand == "+") {
    ppmap <- 50 - apply(ppbins,1,bedIntersect,gnonmap_p)
  } else {
    ppmap <- 50 - apply(ppbins,1,bedIntersect,gnonmap_m)
  }


  # Calculate density of windows
  ppdens <- ifelse(ppmap > 0, pptags/ppmap, 0)
  # Calculate p-value of window
  ppp <- ifelse(ppmap > 0, apply(data.frame(pptags,ppmap),1,function(pp){
      pnbinom(pp[1]-1,alpha,alpha/(alpha+lambda*pp[2]),lower.tail=F)
    }), 1)

  ppbins <- cbind(ppbins,Tags=pptags,Map=ppmap,Density=ppdens,Pvalue=ppp)

  if (gene$Strand == "+") {
    ppbins <- ppbins[with(ppbins,order(Pvalue,-Density,-Tags,Tstart)),]
  } else {
    ppbins <- ppbins[with(ppbins,order(Pvalue,-Density,-Tags,-Tstart)),]
  }
  # The top window according to this ordering is defined as the peak
  ppeak <- ppbins[1,]


  ### Step 2: Identification of Divergent Peak ###

  # Construct divergent windows
  if (ppeak$Pvalue < 0.001) {
    if (gene$Strand == "+") {
      dpstart <- seq(ppeak$Tstart-1000,ppeak$Tstart+1000,5)
    } else {
      dpstart <- seq(ppeak$Tend-1000,ppeak$Tend+1000,5)
    }
  } else {
    if (gene$Strand == "+") {
      dpstart <- seq(gene$Tstart-1000,gene$Tstart+1000,5)
    } else {
      dpstart <- seq(gene$Tend-1000,gene$Tend+1000,5)
    }
  }
  dpend <- dpstart + 50 - 1
  dpbins <- data.frame(Tstart=dpstart,Tend=dpend)

  # Count tags per window and number of mappable bases
  dptags <- apply(dpbins,1,bedCount,tags_anti)

  if (gene$Strand == "+") {
    dpmap <- 50 - apply(dpbins,1,bedIntersect,gnonmap_m)
  } else {
    dpmap <- 50 - apply(dpbins,1,bedIntersect,gnonmap_p)
  }
  # Calculate density of windows
  dpdens <- ifelse(dpmap > 0, dptags/dpmap, 0)
  # Calculate p-value of window
  dpp <- ifelse(dpmap > 0, apply(data.frame(dptags,dpmap),1,function(dp){
      pnbinom(dp[1]-1,alpha,alpha/(alpha+lambda*dp[2]),lower.tail=F)
    }), 1)

  dpbins <- cbind(dpbins,Tags=dptags,Map=dpmap,Density=dpdens,Pvalue=dpp)

  if (gene$Strand == "+") {
    dpbins <- dpbins[with(dpbins,order(Pvalue,-Density,-Tags,-Tstart)),]
  } else {
    dpbins <- dpbins[with(dpbins,order(Pvalue,-Density,-Tags,Tstart)),]
  }
  # The top window according to this ordering is defined as the peak
  dpeak <- dpbins[1,]


  ### Step 3: Calculation of Gene Activity ###
  # Define Gene head as from promoter to 1kb downstream
  # Define Gene body as 1kb downstream from promoter
  # Promoter defined as promoter peak if significant
  # else it is the refgene coordinates

  if (ppeak$Pvalue < 0.001) {
    head_start <- ifelse(gene$Strand == "+", ppeak$Tstart, ppeak$Tend - 1000)
    head_end <- ifelse(gene$Strand == "+", ppeak$Tstart + 1000, ppeak$Tend)
    body_start <- ifelse(gene$Strand == "+", ppeak$Tstart + 1000, gene$Tstart)
    body_end <- ifelse(gene$Strand == "+", gene$Tend, ppeak$Tend - 1000)
  } else {
    head_start <- ifelse(gene$Strand == "+", gene$Tstart, gene$Tend - 1000)
    head_end <- ifelse(gene$Strand == "+", gene$Tstart + 1000, gene$Tend)
    body_start <- ifelse(gene$Strand == "+", gene$Tstart + 1000, gene$Tstart)
    body_end <- ifelse(gene$Strand == "+", gene$Tend, gene$Tend - 1000)
  }

  bodytags <- bedCount(c(body_start,body_end),tags_sense)
  headtags_sense <- bedCount(c(head_start,head_end),tags_sense)
  headtags_anti <- bedCount(c(head_start,head_end),tags_anti)

  if (gene$Strand == "+") {
    bodymap <- body_end - body_start + 1 - bedIntersect(c(body_start,body_end),gnonmap_p)
    headmap_sense <- head_end - head_start + 1 - bedIntersect(c(head_start,head_end),gnonmap_p)
    headmap_anti <- head_end - head_start + 1 - bedIntersect(c(head_start,head_end),gnonmap_m)
  } else {
    bodymap <- body_end - body_start + 1 - bedIntersect(c(body_start,body_end),gnonmap_m)
    headmap_sense <- head_end - head_start + 1 - bedIntersect(c(head_start,head_end),gnonmap_m)
    headmap_anti <- head_end - head_start + 1 - bedIntersect(c(head_start,head_end),gnonmap_p)
  }

  bodydens <- ifelse(bodymap > 0, bodytags/bodymap, 0)
  bodyp <- ifelse(bodymap > 0, pnbinom(bodytags-1,alpha,alpha/(alpha+lambda*bodymap),lower.tail=F), 1)

  headdens_sense <- ifelse(headmap_sense > 0, headtags_sense/headmap_sense, 0)
  headdens_anti <- ifelse(headmap_anti > 0, headtags_anti/headmap_anti, 0)
  headp_sense <- ifelse(headmap_sense > 0, pnbinom(headtags_sense-1,alpha,alpha/(alpha+lambda*headmap_sense),lower.tail=F), 1)
  headp_anti <- ifelse(headmap_anti > 0, pnbinom(headtags_anti-1,alpha,alpha/(alpha+lambda*headmap_anti),lower.tail=F), 1)


  ### Step 4: Identification of Paused Genes ###

  expectedpeak <- floor((ppeak$Tags + bodytags) * ppeak$Map / (ppeak$Map + bodymap))
  expectedbody <- ceiling((ppeak$Tags + bodytags) * bodymap / (ppeak$Map + bodymap))

  # Construct 2x2 matrix for Fisher Exact Test
  F <- matrix( c(ppeak$Tags, bodytags, expectedpeak, expectedbody),
    ncol=2, nrow=2, byrow=T)
# P-value indicating significant pausing
  ppaused <- ifelse(bodymap > 0 || ppeak$Map > 0, fisher.test(F, alternative='g')$p.value, 1)
  
  pindex <- ppeak$Density/bodydens

  type <- ifelse(bodyp < 0.001 && ppaused >=0.001, "I",
          ifelse(bodyp < 0.001 && ppaused < 0.001, "II",
          ifelse(bodyp >=0.001 && ppaused < 0.001, "III",
          ifelse(bodyp >=0.001 && ppaused >=0.001, "IV",
          "NA"))))

  generesult <- data.frame(Gene=gene$Gene,Tname=gene$Tname,Strand=gene$Strand, Ref.Tstart=gene$Tstart,Ref.Tend=gene$Tend,
    Prom.Tstart=ppeak$Tstart,Prom.Tend=ppeak$Tend,Prom.Tags=ppeak$Tags,Prom.Map=ppeak$Map,Prom.Dens=ppeak$Density,Prom.P=ppeak$Pvalue,
    Div.Tstart=dpeak$Tstart,Div.Tend=dpeak$Tend,Div.Tags=dpeak$Tags,Div.Map=dpeak$Map,Div.Dens=dpeak$Density,Div.P=dpeak$Pvalue,
    Head.Tstart=head_start,Head.Tend=head_end,Head.Sense.Tags=headtags_sense,Head.Sense.Map=headmap_sense,Head.Sense.Dens=headdens_sense,Head.Sense.P=headp_sense,
    Head.Anti.Tags=headtags_anti,Head.Anti.Map=headmap_anti,Head.Anti.Dens=headdens_anti,Head.Anti.P=headp_anti,
    Body.Tstart=body_start,Body.Tend=body_end,Body.Tags=bodytags,Body.Map=bodymap,Body.Dens=bodydens,Body.P=bodyp,
    Pausing.Idx=pindex,Pausing.P=ppaused,Type=type)


  return(generesult);
}