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
library("limma")
library("spdep") 
library("adegenet") 

#############################################################
setwd("/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Data") #set working directory

counts <- read.csv(file="isoforms_counts_matrix.counts.matrix", header=T, sep="") #Load expressin matrix from trinity
annot <- read.csv(file="trinotate_annotation_report.csv", header=T)
#annotation <- annotation[match(rownames(counts), annot$Pfam),]

#want to pull up until the 3rd occurrance of ^ in Pfam column

filt <- filterfun(pOverA(0.25,5)) #set filter values for PoverA, P percent of the samples have counts over A
tfil <- genefilter(counts, filt) #create filter for the counts data
keep <- counts[tfil,] #identify transcripts to keep by count filter
gn.keep <- rownames(keep) #identify transcript list
counts.10x <- as.matrix(counts[which(rownames(counts) %in% gn.keep),]) #data filtered in PoverA, P percent of the samples have counts over A
storage.mode(counts.10x) = "integer" #store counts data as integer
sample.info <- read.csv(file="sample_description.csv", header=T, sep=",", row.names=1) #load sample info
sample.info$group <- factor(paste0(sample.info$Spawn, sample.info$Time)) #merge condition and time into group
data <- DESeqDataSetFromMatrix(countData = counts.10x, colData = sample.info, design = ~ group) #create a DESeqDataSet object

# Expression Visualization
rld <- rlog(data, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
head(assay(rld), 3) #view data
sampleDists <- dist(t(assay(rld))) #calculate distance matix
sampleDistMatrix <- as.matrix(sampleDists) #distance matrix
rownames(sampleDistMatrix) <- colnames(rld) #assign row names
colnames(sampleDistMatrix) <- NULL #assign col names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
pheatmap(sampleDistMatrix, #plot matrix
         clustering_distance_rows=sampleDists, #cluster rows
         clustering_distance_cols=sampleDists, #cluster columns
         col=colors) #set colors

plotPCA(rld, intgroup = c("Time", "Spawn")) #plot PCA of samples with all data

# Differential Gene Expression Analysis
#Interaction Test: test of the factor of "group" with all combinations of the original factors as groups
DEG.int <- DESeq(data) #run differential expression test by group using the wald test
DEG.int.res <- results(DEG.int) #save DE results
resultsNames(DEG.int) #view DE results
sig.num <- sum(DEG.int.res$padj <0.1, na.rm=T) #identify the number of significant p values with 10%FDR (padj<0.1)
sig <- subset(DEG.int.res, padj<0.1,) #identify signficant pvalues with 10%FDR
sig.list <- data[which(rownames(data) %in% rownames(sig)),] #subset list of sig transcripts from original count data
rsig <- rlog(sig.list, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size

##### Run PCA
#princomp(rsig$)

##### View DEG Data Heatmap and PCA
PCA.plot <- plotPCA(rsig, intgroup = c("Time", "Spawn")) #Plot PCA of all samples for DEG only
PCA.plot #view plot
PC.info <-PCA.plot$data #extract plotting data
jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/PCA.DEG.jpg")
plot(PC.info$PC1, PC.info$PC2, xlim=c(-50,50), ylim=c(-20, 20), xlab="PC1 82%", ylab="PC2 5%", pch = c(15, 16, 17)[as.numeric(sample.info$Time)], col=c("black", "gray")[sample.info$Spawn], cex=1.3)
legend(x="topleft", 
       bty="n",
       legend = c("00:00", "18:00", "20:00"),
       pch = c(15, 16, 17),
       col=c("black", "black", "black"),
       cex=1)
dev.off()

topVarGenes <- head(order(rowVars(assay(rsig)),decreasing=TRUE),sig.num) #can choose a subset of transcripts for viewing
mat <- assay(rsig)[ topVarGenes, ] #make an expression object
mat <- mat - rowMeans(mat) #difference in expression compared to average across all samples
df <- as.data.frame(colData(rsig)[,c("Spawn","Time")]) #make dataframe
jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Heatmap.DEG.jpg")
pheatmap(mat, annotation_col=df, 
         show_rownames =F, 
         show_colnames =F) #plot heatmap of all DEG by group
dev.off()

# Time point contrasts
DEG.18 <- results(DEG.int, contrast=c("group", "YES18:00", "NO18:00")) #contrasts between SP NSP at each time
DEG.18 #view results
sig.num.18 <- sum(DEG.18$padj <0.1, na.rm=T) #the number of significant p values with 10%FDR
sig.18 <- subset(DEG.18, padj<0.1,) #identify signficant pvalues with 10%FDR
sig.list.18 <- data[which(rownames(data) %in% rownames(sig.18)),] #subset list of sig transcripts from original count data

DEG.20 <- results(DEG.int, contrast=c("group", "YES20:00", "NO20:00")) #contrasts between SP NSP at each time
DEG.20 #view results
sig.num.20 <- sum(DEG.20$padj <0.1, na.rm=T) #the number of significant p values with 10%FDR
sig.20 <- subset(DEG.20, padj<0.1,) #identify signficant pvalues with 10%FDR
sig.list.20 <- data[which(rownames(data) %in% rownames(sig.20)),] #subset list of sig transcripts from original count data

DEG.00 <- results(DEG.int, contrast=c("group", "YES0:00", "NO0:00")) #contrasts between SP NSP at each time
DEG.00 #view results
sig.num.00 <- sum(DEG.00$padj <0.1, na.rm=T) #the number of significant p values with 10%FDR
sig.00 <- subset(DEG.00, padj<0.1,) #identify signficant pvalues with 10%FDR
sig.list.00 <- data[which(rownames(data) %in% rownames(sig.00)),] #subset list of sig transcripts from original count data


#Prepare Data for Venn Diagram
L18 <- data.frame(rownames(sig.list.18)) #identify DE transcript list from 18:00
L18$L18.Count <- 1 #assign all transcripts count of 1
colnames(L18) <- c("Transcript.ID", "L18.Count") #name columns
L20 <- data.frame(rownames(sig.list.20)) #identify DE transcript list from 20:00
L20$L20.Count <- 1 #assign all transcripts count of 1
colnames(L20) <- c("Transcript.ID", "L20.Count") #name columns
L00 <- data.frame(rownames(sig.list.00)) #identify DE transcript list from 00:00
L00$L00.Count <- 1 #assign all transcripts count of 1
colnames(L00) <- c("Transcript.ID", "L00.Count") #name columns
L18.20 <- merge(L18, L20, by="Transcript.ID", all=T ) #merge lists sequentially with no removal
L18.20.00 <- merge(L18.20, L00, by="Transcript.ID", all=T ) #merge lists sequentially with no removal
L18.20.00$L18.Count[is.na(L18.20.00$L18.Count)] <- 0 #assign NA=0 count i.e., not DE
L18.20.00$L20.Count[is.na(L18.20.00$L20.Count)] <- 0 #assign NA=0 count i.e., not DE
L18.20.00$L00.Count[is.na(L18.20.00$L00.Count)] <- 0 #assign NA=0 count i.e., not DE

# Plot Venn Diagram of shared DE between time points
vennData<-L18.20.00[,2:4] #set venn data to DE counts only
a <- vennCounts(vennData) # compute classification counts
colnames(a) <- c('18:00', '20:00', '00:00', "Counts") #set catagory names
jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Venn.DEG.jpg")
vennDiagram(a, main='DEG between groups') #draw venn diagram
title(main='for each time point', cex.main=0.9) #add a title
dev.off()

#Look for match between all lists
shared <- Reduce(intersect, list(rownames(sig.list.18),rownames(sig.list.20),rownames(sig.list.00))) #idenitfy shared DE transcripts between each time point 
shared #view results
write.csv(shared, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/DEG_shared_by_time.csv")

#Found in all time points
intersctions <- subset(L18.20.00, L18.Count==1 & L20.Count==1 & L00.Count==1) #list genes found in all time points

#Unique to each time point
unique.18 <- subset(L18.20.00, L18.Count==1 & L20.Count==0 & L00.Count==0) #list genes unique to 18:00
unique.20 <- subset(L18.20.00, L18.Count==0 & L20.Count==1 & L00.Count==0) #list genes unique to 20:00
unique.00 <- subset(L18.20.00, L18.Count==0 & L20.Count==0 & L00.Count==1) #list genes unique to 00:00

#Shared unique between pairs
unique.18.20 <- subset(L18.20.00, L18.Count==1 & L20.Count==1 & L00.Count==0) #list genes common to 18:00 and 20:00
unique.18.00 <- subset(L18.20.00, L18.Count==1 & L20.Count==0 & L00.Count==1) #list genes common to 18:00 and 00:00
unique.20.00 <- subset(L18.20.00, L18.Count==0 & L20.Count==1 & L00.Count==1) #list genes common to 20:00 and 00:00

#combine ID and rlog normalized counts of genes found in all time points
norm.sig.counts <- as.data.frame(mat)
norm.sig.counts$Transcript.ID <- rownames(norm.sig.counts)
inters <- merge(intersctions, norm.sig.counts, by="Transcript.ID") 
#### Why only 162 and not full 165???

#Need to Decide on plot type for 165 DEG shared between all times
#what do we want to portray?
#heatmap?


#GOSEQ Analysis


#shared transcripts with egg transcriptomes?







