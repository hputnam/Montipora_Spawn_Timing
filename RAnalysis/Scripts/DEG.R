#last modified 20170224
#Copyright 2017 HM Putnam
#Use of this code must be accompanied by a citation to XXXX
#Data should not be used without permission from HM Putnam
#See Readme

rm(list=ls()) # removes all prior objects

#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("genefilter")
library("ggplot2")
library("gplots")

#############################################################
setwd("/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Data") #set working directory


counts <- read.csv(file="isoforms_counts_matrix.counts.matrix", header=T, sep="") #Load expressin matrix from trinity

thresh <- 10 #threshold for expression average across all samples for a transcript to be retained
#filt <- filterfun(pOverA(0.90,10))
#tfil <- genefilter(counts, filt)
#counts.10x <- counts[tfil,]
counts.10x <- as.matrix(counts[rowMeans(counts[,-1])> thresh,]) #data filtered on a expression average threshold across all samples
storage.mode(counts.10x) = "integer" #store counts data as integer
sample.info <- read.csv(file="sample_description.csv", header=T, sep=",", row.names=1) #load sample info
sample.info$group <- factor(paste0(sample.info$Spawn, sample.info$Time)) #merge condition and time into group
data <- DESeqDataSetFromMatrix(countData = counts.10x, colData = sample.info, design = ~ group) #create a DESeqDataSet object

# Expression Visualization
rld <- rlog(data, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalizes wrt library size
head(assay(rld), 3) #view data
sampleDists <- dist(t(assay(rld))) #calculate distance matix
sampleDistMatrix <- as.matrix(sampleDists) #distance matrix
rownames(sampleDistMatrix) <- paste(rld$Time, rld$Spawn, sep="-") #assign row names
colnames(sampleDistMatrix) <- NULL #assign col names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
pheatmap(sampleDistMatrix, #plot matrix
         clustering_distance_rows=sampleDists, #cluster rows
         clustering_distance_cols=sampleDists, #cluster columns
         col=colors) #set colors

plotPCA(rld, intgroup = c("Time", "Spawn"))

# Differential Gene Expression Analysis
#Interaction Test
DEG.int <- DESeq(data)
DEG.int.res <- results(DEG.int)
resultsNames(DEG.int)
sig.num <- sum(DEG.int.res$padj <0.1, na.rm=T) #the number of significant p values with 10%FDR
sig <- subset(DEG.int.res, padj<0.1,)
sig.list <- data[which(rownames(data) %in% rownames(sig)),]
rsig <- rlog(sig.list, blind=FALSE)

##### View DEG Data Heatmap and PCA
plotPCA(rsig, intgroup = c("Time", "Spawn")) #Plot PCA of all samples for DEG only
topVarGenes <- head(order(rowVars(assay(rsig)),decreasing=TRUE),sig.num)
mat <- assay(rsig)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rsig)[,c("Spawn","Time")])
pheatmap(mat, annotation_col=df, show_rownames =F, show_colnames =F)

#Mains and Interaction
DEG.18 <- results(DEG.int, contrast=c("group", "YES18:00", "NO18:00"))
DEG.18
sig.num.18 <- sum(DEG.18$padj <0.1, na.rm=T) #the number of significant p values with 10%FDR
sig.18 <- subset(DEG.18, padj<0.1,)
sig.list.18 <- data[which(rownames(data) %in% rownames(sig.18)),]

DEG.20 <- results(DEG.int, contrast=c("group", "YES20:00", "NO20:00"))
DEG.20
sum(DEG.20$padj <0.1, na.rm=T) #the number of significant p values with 10%FDR
sig.20 <- subset(DEG.20, padj<0.1,)
sig.list.20 <- data[which(rownames(data) %in% rownames(sig.20)),]

DEG.00 <- results(DEG.int, contrast=c("group", "YES0:00", "NO0:00"))
DEG.00
sum(DEG.00$padj <0.1, na.rm=T) #the number of significant p values with 10%FDR
sig.00 <- subset(DEG.00, padj<0.1,)
sig.list.00 <- data[which(rownames(data) %in% rownames(sig.00)),]

#Look for match between lists
shared <- Reduce(intersect, list(rownames(sig.list.18),rownames(sig.list.20),rownames(sig.list.00)))
shared
#the majority of DEG between spawn and not were shared between time points
#the expression patterns of spawning and not spawning corals was very similar
#is this due to egg transcripts? We have 2 libraries we could use to address this...

L18 <- data.frame(rownames(sig.list.18))
L18$L18.Count <- 1
colnames(L18) <- c("Transcript.ID", "L18.Count")
L20 <- data.frame(rownames(sig.list.20))
L20$L20.Count <- 1
colnames(L20) <- c("Transcript.ID", "L20.Count")
L00 <- data.frame(rownames(sig.list.00))
L00$L00.Count <- 1
colnames(L00) <- c("Transcript.ID", "L00.Count")

L18.20 <- merge(L18, L20, by="Transcript.ID", all=T )
L18.20.00 <- merge(L18.20, L00, by="Transcript.ID", all=T ) 

L18.20.00$L18.Count[is.na(L18.20.00$L18.Count)] <- 0
L18.20.00$L20.Count[is.na(L18.20.00$L20.Count)] <- 0
L18.20.00$L00.Count[is.na(L18.20.00$L00.Count)] <- 0

library("limma")
vennData<-L18.20.00[,2:4]
a <- vennCounts(vennData)
colnames(a) <- c('18:00', '20:00', '00:00', "Counts")
vennDiagram(a, main='DEG betwen spawning and not spawning')
title(main='for each time point', cex.main=0.9)

## 
GroupA <- sample(geneNames, 400, replace=FALSE)
GroupB <- sample(geneNames, 750, replace=FALSE)
GroupC <- sample(geneNames, 250, replace=FALSE)
GroupD <- sample(geneNames, 300, replace=FALSE)

venn(list(GrpA=GroupA,GrpB=GroupB,GrpC=GroupC,GrpD=GroupD))

venn()

#Expression clusters
#Plot Expression by Cluster as a function of time



plotMA(DEG.int, main="DESeq2", ylim=c(-2,2))
plotCounts(DEG, gene=which.min(DEG.res$padj), intgroup="var")


