genefile1 <- "/Volumes//AltLab/GRO-seq//Results//Neuron//mm9-trim16-v3//gene_expression-new//totalGeneExpression.txt"
genefile2 <- "/Volumes/AltLab/GRO-seq/Results/Tviral1/mm9-trim15-v3/gene-expression-new//totalGeneExpression.txt"

gene1 <- read.delim(genefile1,header=T,as.is=T)
gene2 <- read.delim(genefile2,header=T,as.is=T)

enrich1 <- gene1$Type %in% c("I","II") & gene2$Type %in% c("III","IV")
enrich2 <- gene2$Type %in% c("I","II") & gene1$Type %in% c("III","IV")

name1 <- "Neuron"
name2 <- "B-cell"

part1 <- gene1[,c("Type","Prom.Dens","Prom.P","Head.Sense.Dens","Head.Sense.P","Body.Dens","Body.P")]
part2 <- gene2[,c("Type","Prom.Dens","Prom.P","Head.Sense.Dens","Head.Sense.P","Body.Dens","Body.P")]
colnames(part1) <- paste(name1,".",colnames(part1),sep="")
colnames(part2) <- paste(name2,".",colnames(part2),sep="")

output1 <- cbind(gene1[,1:5],part1,part2)[enrich1,]
output2 <- cbind(gene1[,1:5],part1,part2)[enrich2,]

outfile1 <- "/Volumes//AltLab/GRO-seq/Results/Comparison/Neuron-Bcell.txt"
outfile2 <- "/Volumes//AltLab/GRO-seq/Results/Comparison/Bcell-Neuron.txt"

write.table(output1,outfile1,col.names=T,row.names=F,quote=F,sep="\t")
write.table(output2,outfile2,col.names=T,row.names=F,quote=F,sep="\t")

