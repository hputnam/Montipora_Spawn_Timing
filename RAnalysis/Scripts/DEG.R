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
#biocLite("GSEABase")
library(GSEABase)
#biocLite("pathview")
library("pathview")
#biocLite("goseq")
library("goseq")
library("GO.db")
library("UpSetR")
library("reshape2")
library(lattice)
library(latticeExtra)
#library("grid")
#library("gridExtra")
#biocLite("KEGGgraph")
#library("KEGGgraph")


#############################################################
setwd("/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Data") #set working directory
#############################################################

##### Differential Expression analysis #####
counts <- read.csv(file="adult_bundle_isoforms_counts_matrix.counts.matrix", header=T, sep="") #Load expressin matrix from trinity
counts <- counts[,1:18] #Adult data only
filt <- filterfun(pOverA(0.50,5)) #set filter values for PoverA, P percent of the samples have counts over A
tfil <- genefilter(counts, filt) #create filter for the counts data
keep <- counts[tfil,] #identify transcripts to keep by count filter
gn.keep <- rownames(keep) #identify transcript list
counts.5x <- as.matrix(counts[which(rownames(counts) %in% gn.keep),]) #data filtered in PoverA, P percent of the samples have counts over A
write.csv(counts.5x, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/filtered_counts.csv")
#ft.counts <- rlog(counts.5x, blind=FALSE)
#write.csv(ft.counts, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/filtered_transformed_counts.csv")
storage.mode(counts.5x) = "integer" #store counts data as integer
sample.info <- read.csv(file="sample_description.csv", header=T, sep=",", row.names=1) #load sample info
sample.info$group <- factor(paste0(sample.info$Spawn, sample.info$Time)) #merge condition and time into group
annot <- read.csv("trinotate_annotation_report_wo_seq.csv", header=T, sep=",") #biological annotation information
colnames(annot)[2] <- "Transcript.ID"
coral <- read.csv("Taxa/coral.outfmt6", header=F, sep="") #taxonomic annotation information
colnames(coral) <- c("Transcript.ID", "coral_subject_id", "coral_pct_identity", "coral_aln_length", "coral_n_of_mismatches", "coral_gap_openings", "coral_q_start", "coral_q_end","coral_s_start", "coral_s_end", "coral_e_value", "coral_bit_score")
sym <- read.csv("Taxa/sym.outfmt6", header=F, sep="") #taxonomic annotation information
colnames(sym) <- c("Transcript.ID", "sym_subject_id", "sym_pct_identity", "sym_aln_length", "sym_n_of_mismatches", "sym_gap_openings", "sym_q_start", "sym_q_end","sym_s_start", "sym_s_end", "sym_e_value", "sym_bit_score")
bact <- read.csv("Taxa/bacteria.outfmt6", header=F, sep="") #taxonomic annotation information
colnames(bact) <- c("Transcript.ID", "bact_subject_id", "bact_pct_identity", "bact_aln_length", "bact_n_of_mismatches", "bact_gap_openings", "bact_q_start", "bact_q_end","bact_s_start", "bact_s_end", "bact_e_value", "bact_bit_score")
vir <-  read.csv("Taxa/viruses.outfmt6", header=F, sep="") #taxonomic annotation information
colnames(vir) <- c("Transcript.ID", "vir_subject_id", "vir_pct_identity", "vir_aln_length", "vir_n_of_mismatches", "vir_gap_openings", "vir_q_start", "vir_q_end","vir_s_start", "vir_s_end", "vir_e_value", "vir_bit_score")

data <- DESeqDataSetFromMatrix(countData = counts.5x, colData = sample.info, design = ~ group) #create a DESeqDataSet object

# Expression Visualization
rld <- rlog(data, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
head(assay(rld), 3) #view data
sampleDists <- dist(t(assay(rld))) #calculate distance matix
sampleDistMatrix <- as.matrix(sampleDists) #distance matrix
rownames(sampleDistMatrix) <- colnames(rld) #assign row names
colnames(sampleDistMatrix) <- NULL #assign col names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
pheatmap(sampleDistMatrix, #plot matrix of expression similarity
         clustering_distance_rows=sampleDists, #cluster rows
         clustering_distance_cols=sampleDists, #cluster columns
         col=colors) #set colors

plotPCA(rld, intgroup = c("Spawn", "Time")) #plot PCA of samples with all data

# Differential Gene Expression Analysis
#Interaction Test: test of the factor of "group" with all combinations of the original factors as groups
DEG.int <- DESeq(data) #run differential expression test by group using the wald test 
DEG.int.res <- results(DEG.int) #save DE results
resultsNames(DEG.int) #view DE results
sig.num <- sum(DEG.int.res$padj <0.05, na.rm=T) #identify the number of significant p values with 5%FDR (padj<0.05)
sig <- subset(DEG.int.res, padj<0.05,) #identify signficant pvalues with 5%FDR
sig.list <- data[which(rownames(data) %in% rownames(sig)),] #subset list of sig transcripts from original count data
rsig <- rlog(sig.list, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
write.csv(counts(sig.list), file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/DEG_5FDR.all.csv")


##### View DEG Data Heatmap and PCA #####
PCA.plot <- plotPCA(rsig, intgroup = c("Spawn", "Time")) #Plot PCA of all samples for DEG only
PCA.plot #view plot
PC.info <-PCA.plot$data #extract plotting data
pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/PCA.DEG.pdf")
plot(PC.info$PC1, PC.info$PC2, xlim=c(-50,50), ylim=c(-16, 10), xlab="PC1 86%", ylab="PC2 5%", col = c("lightpink2", "steelblue1","yellow3")[as.numeric(PC.info$Time)], pch=c(16, 17)[as.numeric(PC.info$Spawn)], cex=1.3)
legend(x="top", 
       bty="n",
       legend = c("00:00", "18:00", "20:00", "Not Spawning", "Spawning"),
       text.col = c("lightpink2","steelblue1","yellow3", "black", "black"),
       pch = c(15, 15, 15, 16, 17),
       col = c("white","white","white", "black", "black"),
       cex=1)
dev.off()

topVarGenes <- head(order(rowVars(assay(rsig)),decreasing=TRUE),sig.num) #sort by decreasing sig
mat <- assay(rsig)[ topVarGenes, ] #make an expression object
mat <- mat - rowMeans(mat) #difference in expression compared to average across all samples
col.order <- c("T4.17.includes.adapter_S1","T4.1.includes.adapter_S1","T4.6.includes.adapter_S1",
               "T5.17.includes.adapter_S2","T5.1.includes.adapter_S2","T5.6.includes.adapter_S2",
               "T7.17.includes.adapter_S3","T7.1.includes.adapter_S3","T7.6.includes.adapter_S3",
               "T4.10.includes.adapter_S1","T4.16.includes.adapter_S1","T4.8.includes.adapter_S1",
               "T5.10.includes.adapter_S2","T5.16.includes.adapter_S2","T5.8.includes.adapter_S2",
               "T7.10.includes.adapter_S3","T7.16.includes.adapter_S3","T7.8.includes.adapter_S3")
mat <- mat[,col.order]
df <- as.data.frame(colData(rsig)[,c("Spawn","Time")]) #make dataframe
df <- df[order(df$Spawn),]
pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/DEG_Heatmap.pdf")
pheatmap(mat, annotation_col=df, clustering_method = "average", 
         clustering_distance_rows="euclidean", show_rownames =FALSE, cluster_cols=F,
         show_colnames =F) #plot heatmap of all DEG by group
dev.off()


#check for optimal number of clusters for reporting expression patterns
set.seed(123) 
k.max <- 10 # Compute and plot wss for k = 2 to k = 10
wss <- sapply(1:k.max, function(k){kmeans(mat, k, nstart=50,iter.max = 15 )$tot.withinss}) #within sum of squares
#plot percentage of variance explained as a function of the number of clusters
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


#set annotation colors on heatmap factors
ann_colors <- list(Spawn = c(YES="darkorange", NO="burlywood4"), Time = c("18:00"="steelblue1", "20:00"="yellow3", "0:00"="lightpink2"))
#"lightpink2", "yellow3", "steelblue1"

#save heatmap results
DEG.Heat.res <- pheatmap(mat, color = colorRampPalette(c("darkblue", "cyan", "yellow"))(100), annotation_col=df, annotation_colors =ann_colors, clustering_method = "average", 
                clustering_distance_rows="euclidean", show_rownames =FALSE, 
                show_colnames =F)
#set optimal cluster number based on results above
knum <-6
DEG.clust <- cbind(mat, cluster = cutree(DEG.Heat.res$tree_row, k = knum))[DEG.Heat.res$tree_row[["order"]]]
DEG.clust <- as.data.frame(DEG.clust)
cluster.order <- cbind(mat[c(DEG.Heat.res$tree_row[["order"]]),DEG.Heat.res$tree_col[["order"]]],cluster=cutree(DEG.Heat.res$tree_row, k=knum)[DEG.Heat.res$tree_row[["order"]]])

cluster.order <- as.data.frame(cluster.order)
DEG.clust1 <- cluster.order[cluster.order$cluster == '1',]
DEG.clust1 <- DEG.clust1[,-19] 
DEG.clust2 <- cluster.order[cluster.order$cluster == '2',]
DEG.clust2 <- DEG.clust2[,-19] 
DEG.clust3 <- cluster.order[cluster.order$cluster == '3',]
DEG.clust3 <- DEG.clust3[,-19] 
DEG.clust4 <- cluster.order[cluster.order$cluster == '4',]
DEG.clust4 <- DEG.clust4[,-19] 
DEG.clust5 <- cluster.order[cluster.order$cluster == '5',]
DEG.clust5 <- DEG.clust5[,-19] 
DEG.clust6 <- cluster.order[cluster.order$cluster == '6',]
DEG.clust6 <- DEG.clust6[,-19] 

CL1 <- t(DEG.clust1)
CL1 <- CL1[ order(row.names(CL1)), ]
sample.info <- sample.info[order(row.names(sample.info)), ]
CL1 <- cbind(CL1, sample.info)

CL2 <- t(DEG.clust2)
CL2 <- CL2[ order(row.names(CL2)), ]
sample.info <- sample.info[order(row.names(sample.info)), ]
CL2 <- cbind(CL2, sample.info)

CL3 <- t(DEG.clust3)
CL3 <- CL3[ order(row.names(CL3)), ]
sample.info <- sample.info[order(row.names(sample.info)), ]
CL3 <- cbind(CL3, sample.info)

CL4 <- t(DEG.clust4)
CL4 <- CL4[ order(row.names(CL4)), ]
sample.info <- sample.info[order(row.names(sample.info)), ]
CL4 <- cbind(CL4, sample.info)

CL5 <- t(DEG.clust5)
CL5 <- CL5[ order(row.names(CL5)), ]
sample.info <- sample.info[order(row.names(sample.info)), ]
CL5 <- cbind(CL5, sample.info)

CL6 <- t(DEG.clust6)
CL6 <- CL6[ order(row.names(CL6)), ]
sample.info <- sample.info[order(row.names(sample.info)), ]
CL6 <- cbind(CL6, sample.info)

cluster1 <- melt(CL1, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp") 
cluster2 <- melt(CL2, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp") 
cluster3 <- melt(CL3, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp") 
cluster4 <- melt(CL4, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster5 <- melt(CL5, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster6 <- melt(CL6, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")

C1.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster1, mean, na.rm = TRUE)
C2.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster2, mean, na.rm = TRUE)
C3.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster3, mean, na.rm = TRUE)
C4.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster4, mean, na.rm = TRUE)
C5.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster5, mean, na.rm = TRUE)
C6.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster6, mean, na.rm = TRUE)

#Want to plot expression plots for each cluster with all data points and lines of mean of cluster
#need to change order of x axis to 18, 20, 00
#need to create set colors and legends for SP vs NSP
#YES="burlywood4", NO="darkorange"
pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Clusters.pdf", width=2, height=11, useDingbats = FALSE)
par(mfrow=c(6,1))
par(mai=c(0.3,0.1,0,0.3))
#sets the bottom, left, top and right margins respectively of the plot region in number of lines of text.

boxplot(Rel.Exp ~ Spawn * TimePoint, data=cluster1, col=c("burlywood4", "darkorange") ,  xaxt="n", yaxt="n", ylim=c(-6,6))
axis(4)
#stripchart(cluster1$Rel.Exp ~ cluster1$Spawn * cluster1$TimePoint, add=T, method = "jitter", jitter = 0.1, vertical = T, pch=16, cex=0.2, col="black")
#points(C1.means$Rel.Exp,   col="black", pch=15, cex=1)

boxplot(Rel.Exp ~ Spawn * TimePoint, data=cluster3, col=c("burlywood4", "darkorange") ,  xaxt="n", yaxt="n", ylim=c(-6,6))
axis(4)
#stripchart(cluster3$Rel.Exp ~ cluster3$Spawn * cluster3$TimePoint, add=T, method = "jitter", jitter = 0.1, vertical = T, pch=16, cex=0.2, col="black")
#points(C3.means$Rel.Exp,   col="black",pch=15, cex=1)

boxplot(Rel.Exp ~ Spawn * TimePoint, data=cluster2, col=c("burlywood4", "darkorange") ,  xaxt="n", yaxt="n", ylim=c(-6,6))
axis(4)
#stripchart(cluster2$Rel.Exp ~ cluster2$Spawn * cluster2$TimePoint, add=T, method = "jitter", jitter = 0.1, vertical = T, pch=16, cex=0.2, col="black")
#points(C2.means$Rel.Exp,   col="black",pch=15, cex=1)

boxplot(Rel.Exp ~ Spawn * TimePoint, data=cluster5, col=c("burlywood4", "darkorange") ,  xaxt="n", yaxt="n", ylim=c(-6,6))
axis(4)
#stripchart(cluster5$Rel.Exp ~ cluster5$Spawn * cluster5$TimePoint, add=T, method = "jitter", jitter = 0.1, vertical = T, pch=16, cex=0.2, col="black")
#points(C5.means$Rel.Exp,   col="black",pch=15, cex=1)

boxplot(Rel.Exp ~ Spawn * TimePoint, data=cluster6, col=c("burlywood4", "darkorange") ,  xaxt="n", yaxt="n", ylim=c(-6,6))
axis(4)
#stripchart(cluster6$Rel.Exp ~ cluster6$Spawn * cluster6$TimePoint, add=T, method = "jitter", jitter = 0.1, vertical = T, pch=16, cex=0.2, col="black")
#points(C6.means$Rel.Exp,   col="black",pch=15, cex=1)

boxplot(Rel.Exp ~ Spawn * TimePoint, data=cluster4, col=c("burlywood4", "darkorange") ,  xaxt="n", yaxt="n", ylim=c(-6,6))
axis(4)
axis(1, at = c(1.5 , 3.5 , 5.5), labels = c("18:00","20:00","00:00"), tick=FALSE , cex=2)
#stripchart(cluster4$Rel.Exp ~ cluster4$Spawn * cluster4$TimePoint, add=T, method = "jitter", jitter = 0.1, vertical = T, pch=16, cex=0.2, col="black")
#points(C4.means$Rel.Exp,   col="black",pch=15, cex=1)
dev.off()

row.annots <- as.data.frame(cluster.order$cluster)
rownames(row.annots) <- rownames(cluster.order)

pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/DEG_Heatmap.pdf", width=4, height=11)
pheatmap(mat, color = colorRampPalette(c("darkblue", "cyan", "yellow"))(100), annotation_col=df, annotation_colors =ann_colors, clustering_method = "average", 
         clustering_distance_rows="euclidean", show_rownames =FALSE,  cutree_rows=6, cluster_cols=F,
         show_colnames =F) #plot heatmap of all DEG by group
dev.off()

#Cluster Counts and Annotations
DEG.clust1
C1.counts <- as.data.frame(counts.5x[which(rownames(counts.5x) %in% rownames(DEG.clust1)),])
C1.counts$Transcript.ID <- rownames(C1.counts)
C1.counts.annot <- merge(C1.counts, annot[!duplicated(annot$Transcript.ID),], by="Transcript.ID")
C1.counts.annot <- merge(C1.counts.annot, coral[!duplicated(coral$Transcript.ID),], by="Transcript.ID")
#C1.counts.annot <- merge(C1.counts.annot, sym[!duplicated(sym$Transcript.ID),], by="Transcript.ID")
#C1.counts.annot <- merge(C1.counts.annot, bact[!duplicated(bact$Transcript.ID),], by="Transcript.ID")
#C1.counts.annot <- merge(C1.counts.annot, vir[!duplicated(vir$Transcript.ID),], by="Transcript.ID")
write.csv(C1.counts.annot, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Cluster1_rawcounts_annotated.csv")

DEG.clust2
C2.counts <- as.data.frame(counts.5x[which(rownames(counts.5x) %in% rownames(DEG.clust2)),])
C2.counts$Transcript.ID <- rownames(C2.counts)
C2.counts.annot <- merge(C2.counts, annot[!duplicated(annot$Transcript.ID),], by="Transcript.ID")
C2.counts.annot <- merge(C2.counts.annot, coral[!duplicated(coral$Transcript.ID),], by="Transcript.ID")
#Lost one transcript ID during merge, need to retain
#There is no taxonomic annotation for TRINITY_DN57574_c0_g1_i3 
#C2.counts.annot <- merge(C2.counts.annot, sym[!duplicated(sym$Transcript.ID),], by="Transcript.ID")
#C2.counts.annot <- merge(C2.counts.annot, bact[!duplicated(bact$Transcript.ID),], by="Transcript.ID")
#C2.counts.annot <- merge(C2.counts.annot, vir[!duplicated(vir$Transcript.ID),], by="Transcript.ID")
write.csv(C2.counts.annot, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Cluster2_rawcounts_annotated_X.csv")

DEG.clust3
C3.counts <- as.data.frame(counts.5x[which(rownames(counts.5x) %in% rownames(DEG.clust3)),])
C3.counts$Transcript.ID <- rownames(C3.counts)
C3.counts.annot <- merge(C3.counts, annot[!duplicated(annot$Transcript.ID),], by="Transcript.ID")
C3.counts.annot <- merge(C3.counts.annot, coral[!duplicated(coral$Transcript.ID),], by="Transcript.ID")
#C3.counts.annot <- merge(C3.counts.annot, sym[!duplicated(sym$Transcript.ID),], by="Transcript.ID")
#C3.counts.annot <- merge(C3.counts.annot, bact[!duplicated(bact$Transcript.ID),], by="Transcript.ID")
#C3.counts.annot <- merge(C3.counts.annot, vir[!duplicated(vir$Transcript.ID),], by="Transcript.ID")
write.csv(C3.counts.annot, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Cluster3_rawcounts_annotated.csv")

DEG.clust4
C4.counts <- as.data.frame(counts.5x[which(rownames(counts.5x) %in% rownames(DEG.clust4)),])
C4.counts$Transcript.ID <- rownames(C4.counts)
C4.counts.annot <- merge(C4.counts, annot[!duplicated(annot$Transcript.ID),], by="Transcript.ID")
C4.counts.annot <- merge(C4.counts.annot, coral[!duplicated(coral$Transcript.ID),], by="Transcript.ID")
#C4.counts.annot <- merge(C4.counts.annot, sym[!duplicated(sym$Transcript.ID),], by="Transcript.ID")
#C4.counts.annot <- merge(C4.counts.annot, bact[!duplicated(bact$Transcript.ID),], by="Transcript.ID")
#C4.counts.annot <- merge(C4.counts.annot, vir[!duplicated(vir$Transcript.ID),], by="Transcript.ID")
write.csv(C4.counts.annot, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Cluster4_rawcounts_annotated.csv")

DEG.clust5
C5.counts <- as.data.frame(t(counts.5x[which(rownames(counts.5x) %in% rownames(DEG.clust5)),]))
rownames(C5.counts) <- rownames(DEG.clust5)
C5.counts$Transcript.ID <- rownames(DEG.clust5)
C5.counts.annot <- merge(C5.counts, annot[!duplicated(annot$Transcript.ID),], by="Transcript.ID")
C5.counts.annot <- merge(C5.counts.annot, coral[!duplicated(coral$Transcript.ID),], by="Transcript.ID")
#C5.counts.annot <- merge(C5.counts.annot, sym[!duplicated(sym$Transcript.ID),], by="Transcript.ID")
#C5.counts.annot <- merge(C5.counts.annot, bact[!duplicated(bact$Transcript.ID),], by="Transcript.ID")
#C5.counts.annot <- merge(C5.counts.annot, vir[!duplicated(vir$Transcript.ID),], by="Transcript.ID")
write.csv(C5.counts.annot, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Cluster5_rawcounts_annotated.csv")

DEG.clust6
C6.counts <- as.data.frame(counts.5x[which(rownames(counts.5x) %in% rownames(DEG.clust6)),])
C6.counts$Transcript.ID <- rownames(C6.counts)
C6.counts.annot <- merge(C6.counts, annot[!duplicated(annot$Transcript.ID),], by="Transcript.ID")
C6.counts.annot <- merge(C6.counts.annot, coral[!duplicated(coral$Transcript.ID),], by="Transcript.ID")
#C6.counts.annot <- merge(C6.counts.annot, sym[!duplicated(sym$Transcript.ID),], by="Transcript.ID")
#C6.counts.annot <- merge(C6.counts.annot, bact[!duplicated(bact$Transcript.ID),], by="Transcript.ID")
#C6.counts.annot <- merge(C6.counts.annot, vir[!duplicated(vir$Transcript.ID),], by="Transcript.ID")
write.csv(C6.counts.annot, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Cluster6_rawcounts_annotated.csv")



##### CONTRASTING SPAWNING AND NOT SPAWNING CORALS AT EACH TIME POINT #####
#pairwise comparisons of DEGs for SP and NSP corals at each collection point
DEG.18 <- results(DEG.int, contrast=c("group", "YES18:00", "NO18:00")) #contrasts between SP NSP at each time
DEG.18 #view results
sig.num.18 <- sum(DEG.18$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
sig.18 <- subset(DEG.18, padj<0.05,) #identify signficant pvalues with 1%FDR
sig.list.18 <- data[which(rownames(data) %in% rownames(sig.18)),] #subset list of sig transcripts from original count data

DEG.20 <- results(DEG.int, contrast=c("group", "YES20:00", "NO20:00")) #contrasts between SP NSP at each time
DEG.20 #view results
sig.num.20 <- sum(DEG.20$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
sig.20 <- subset(DEG.20, padj<0.05,) #identify signficant pvalues with 1%FDR
sig.list.20 <- data[which(rownames(data) %in% rownames(sig.20)),] #subset list of sig transcripts from original count data

DEG.00 <- results(DEG.int, contrast=c("group", "YES0:00", "NO0:00")) #contrasts between SP NSP at each time
DEG.00 #view results
sig.num.00 <- sum(DEG.00$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
sig.00 <- subset(DEG.00, padj<0.05,) #identify signficant pvalues with 1%FDR
sig.list.00 <- data[which(rownames(data) %in% rownames(sig.00)),] #subset list of sig transcripts from original count data


#Prepare DEG Data for Venn Diagram or Intersection Plot
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
#pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Venn.DEG.pdf")
#vennDiagram(a, main='DEG between Spawning and NonSpawning Corals', circle.col=c("yellow", "red", "green")) #draw venn diagram
#dev.off()

#Group Intersections
colnames(vennData) <- c("18:00","20:00", "00:00") #set catagory names
pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/SP_NSP_DEG_Intersections.pdf", useDingbats = TRUE)
upset(vennData, sets = c('18:00','20:00', '00:00'), sets.bar.color = c("lightpink2", "yellow3", "steelblue1"), main.bar.color = c("darkgrey", "black", "black", "black", "steelblue1", "yellow3", "lightpink2"), order.by = "degree")
dev.off()

#Look for match between all lists
shared.all <- Reduce(intersect, list(rownames(sig.list.18),rownames(sig.list.20),rownames(sig.list.00))) #idenitfy shared DE transcripts between each time point 
shared.all #view results
shared.all <- as.data.frame(shared.all)
colnames(shared.all)[1] <- "Transcript.ID"
write.csv(shared.all, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/DEG_shared_18_20_00.csv")

#DEG between spawning and not spawning corals Found in all time points
intersections <- subset(L18.20.00, L18.Count==1 & L20.Count==1 & L00.Count==1) #list DEG found in all time points

#DEG between spawning and not spawning corals Unique to each time point
unique.18 <- subset(L18.20.00, L18.Count==1 & L20.Count==0 & L00.Count==0) #list genes unique to 18:00
unique.20 <- subset(L18.20.00, L18.Count==0 & L20.Count==1 & L00.Count==0) #list genes unique to 20:00
unique.00 <- subset(L18.20.00, L18.Count==0 & L20.Count==0 & L00.Count==1) #list genes unique to 00:00

#DEG between spawning and not spawning corals Shared unique between pairs
unique.18.20 <- subset(L18.20.00, L18.Count==1 & L20.Count==1 & L00.Count==0) #list genes common to 18:00 and 20:00
unique.18.00 <- subset(L18.20.00, L18.Count==1 & L20.Count==0 & L00.Count==1) #list genes common to 18:00 and 00:00
unique.20.00 <- subset(L18.20.00, L18.Count==0 & L20.Count==1 & L00.Count==1) #list genes common to 20:00 and 00:00

#combine ID and rlog normalized counts of genes found in all time points
norm.sig.counts <- as.data.frame(mat)
norm.sig.counts$Transcript.ID <- rownames(norm.sig.counts)
inters <- merge(intersections, norm.sig.counts, by="Transcript.ID") 

###### Combine with DE list with annotation data #####
anot.inters <- merge(inters,  annot, by="Transcript.ID") #merge sig intersection list and annotations
# anot.inters$KG <- sapply(strsplit(as.character(anot.inters$Kegg), split="\\`"), "[", 1)
# anot.inters$KG <- gsub(".*:","",anot.inters$KG)
# anot.inters$KO <- sapply(strsplit(as.character(anot.inters$Kegg), split="\\`"), "[", 2)
# anot.inters$KO <- gsub(".*:","",anot.inters$KO)

annot.intersect_18_20_00 <- merge(shared.all, annot, by="Transcript.ID")
write.csv(annot.intersect_18_20_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_shared_DEG_18_20_00.csv")

annot.intersect_18_20 <- merge(unique.18.20, annot, by="Transcript.ID")
unique(annot.intersect_18_20$Transcript.ID)
write.csv(annot.intersect_18_20, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_shared_DEG_18_20.csv")

annot.intersect_18_00 <- merge(unique.18.00, annot, by="Transcript.ID")
unique(annot.intersect_18_00$Transcript.ID)
write.csv(annot.intersect_18_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_shared_DEG_18_00.csv")

annot.intersect_20_00 <- merge(unique.20.00, annot, by="Transcript.ID")
unique(annot.intersect_20_00$Transcript.ID)
write.csv(annot.intersect_20_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_shared_DEG_20_00.csv")

annot.unique.18 <- merge(unique.18, annot, by="Transcript.ID")
write.csv(annot.unique.18, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_unique_18.csv")

annot.unique.20 <- merge(unique.20, annot, by="Transcript.ID")
write.csv(annot.unique.20, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_unique_20.csv")

annot.unique.00 <- merge(unique.00, annot, by="Transcript.ID")
write.csv(annot.unique.00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_unique_00.csv")


##### GO ENRICHMENT ANALYSIS #####

#Analysis for All DEGs (both Host and Sym together)
#read in data files

##### INTERSECTIONS #####
ALL<-row.names(data) #set the all transcripts list
DEG <- as.character(intersections$Transcript.ID) #set the enrichment test list
LENGTH <-read.table("Trinity.fasta.seq_lens", sep = "", header=F) #reads in a list of gene lengths
gn.keep <- as.data.frame(gn.keep) #filter length to subset same count filter as above
colnames(gn.keep) <- c("V1") #name columns
lens <- as.data.frame(LENGTH) #set lengths
LENGTH <- merge(as.data.frame(LENGTH), gn.keep, by="V1") #merge lengths and transcripts
GO <- read.table("go_annotations.txt", sep = "", header=F, stringsAsFactors=FALSE) #id GO annotations
GO <- as.data.frame(GO) #filter GO to subset same count filter as above
GOs <- merge(as.data.frame(GO), gn.keep, by="V1", all=T) #merge GO, ID
GOs[is.na(GOs)] <- "unknown" #list NA as unknown
GOs <- merge(GOs, LENGTH, by="V1") #merge
GOs <- GOs[,1:2] #remove length 
splitted <- strsplit(as.character(GOs$V2.x), ",") #slit into multiple GO ids
GO.terms <- data.frame(v1 = rep.int(GOs$V1, sapply(splitted, length)), v2 = unlist(splitted)) #list all GOs with each assigned transcript
IDs <- row.names(data) #set ID names

#change contig lists to vectors
ALL.vector <-c(t(ALL))
DEG.vector <-c(t(DEG))
ID.vector <- LENGTH$V1
LENGTH.vector <-LENGTH$V2

gene.vector=as.integer(ALL.vector%in%DEG.vector) #Construct new vector with 1 for DEG and 0 for others
names(gene.vector)=ALL.vector #set names
DEG.pwf<-nullp(gene.vector, ID.vector, bias.data=LENGTH.vector) #weight vector by length of gene

#Find enriched GO terms, "selection-unbiased testing for category enrichment amongst differentially expressed (DE) genes for RNA-seq data"
GO.wall<-goseq(DEG.pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)

#How many enriched GO terms do we have
class(GO.wall)
head(GO.wall)
tail(GO.wall)
nrow(GO.wall)

#Find only enriched GO terms that are statistically significant at cutoff 
enriched.GO.05.a<-GO.wall$category[GO.wall$over_represented_pvalue<.05]
enriched.GO.05<-data.frame(enriched.GO.05.a)
colnames(enriched.GO.05) <- c("category")
enriched.GO.05 <- merge(enriched.GO.05, GO.wall, by="category")

MF.INT <- subset(enriched.GO.05, ontology=="MF")
MF.INT <- MF.INT[order(-MF.INT$numDEInCat),]
CC.INT <- subset(enriched.GO.05, ontology=="CC")
CC.INT <- CC.INT[order(-CC.INT$numDEInCat),]
BP.INT <- subset(enriched.GO.05, ontology=="BP")
BP.INT <- BP.INT[order(-BP.INT$numDEInCat),]
write.csv(MF.INT , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/MF_Sig_Enriched_GO.05_INT.csv")
write.csv(CC.INT , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/CC_Sig_Enriched_GO.05_INT.csv")
write.csv(BP.INT , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/BP_Sig_Enriched_GO.05_INT.csv")


##### UNIQUE 18 #####
ALL<-row.names(data)
DEG <- as.character(unique.18$Transcript.ID) #set the enrichment test list
LENGTH <-read.table("Trinity.fasta.seq_lens", sep = "", header=F) #reads in a list of gene lengths
gn.keep <- as.data.frame(gn.keep) #filter length to subset same count filter as above
colnames(gn.keep) <- c("V1") #name columns
lens <- as.data.frame(LENGTH) #set lengths
LENGTH <- merge(as.data.frame(LENGTH), gn.keep, by="V1") #merge lengths and transcripts
GO <- read.table("go_annotations.txt", sep = "", header=F, stringsAsFactors=FALSE) #id GO annotations
GO <- as.data.frame(GO) #filter GO to subset same count filter as above
GOs <- merge(as.data.frame(GO), gn.keep, by="V1", all=T) #merge GO, ID
GOs[is.na(GOs)] <- "unknown" #list NA as unknown
GOs <- merge(GOs, LENGTH, by="V1") #merge
GOs <- GOs[,1:2] #remove length 
splitted <- strsplit(as.character(GOs$V2.x), ",") #slit into multiple GO ids
GO.terms <- data.frame(v1 = rep.int(GOs$V1, sapply(splitted, length)), v2 = unlist(splitted)) #list all GOs with each assigned transcript
IDs <- row.names(data) #set ID names

#change contig lists to vectors
ALL.vector <-c(t(ALL))
DEG.vector <-c(t(DEG))
ID.vector <- LENGTH$V1
LENGTH.vector <-LENGTH$V2

gene.vector=as.integer(ALL.vector%in%DEG.vector) #Construct new vector with 1 for DEG and 0 for others
names(gene.vector)=ALL.vector #set names
DEG.pwf<-nullp(gene.vector, ID.vector, bias.data=LENGTH.vector) #weight vector by length of gene

#Find enriched GO terms, "selection-unbiased testing for category enrichment amongst differentially expressed (DE) genes for RNA-seq data"
GO.wall<-goseq(DEG.pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)

#How many enriched GO terms do we have
class(GO.wall)
head(GO.wall)
tail(GO.wall)
nrow(GO.wall)

#Find only enriched GO terms that are statistically significant at cutoff 
enriched.GO.05.a<-GO.wall$category[GO.wall$over_represented_pvalue<.05]
enriched.GO.05<-data.frame(enriched.GO.05.a)
colnames(enriched.GO.05) <- c("category")
enriched.GO.05 <- merge(enriched.GO.05, GO.wall, by="category")

MF.18 <- subset(enriched.GO.05, ontology=="MF")
MF.18 <- MF.18[order(-MF.18$numDEInCat),]
CC.18 <- subset(enriched.GO.05, ontology=="CC")
CC.18 <- CC.18[order(-CC.18$numDEInCat),]
BP.18 <- subset(enriched.GO.05, ontology=="BP")
BP.18 <- BP.18[order(-BP.18$numDEInCat),]
write.csv(MF.18 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/MF_Sig_Enriched_GO.05_18.csv")
write.csv(CC.18 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/CC_Sig_Enriched_GO.05_18.csv")
write.csv(BP.18 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/BP_Sig_Enriched_GO.05_18.csv")



##### UNIQUE 20 #####
ALL<-row.names(data)
DEG <- as.character(unique.20$Transcript.ID) #row.names(sig) 
LENGTH <-read.table("Trinity.fasta.seq_lens", sep = "", header=F) #reads in a list of gene lengths
gn.keep <- as.data.frame(gn.keep) #filter length to subset same count filter as above
colnames(gn.keep) <- c("V1") #name columns
lens <- as.data.frame(LENGTH) #set lengths
LENGTH <- merge(as.data.frame(LENGTH), gn.keep, by="V1") #merge lengths and transcripts
GO <- read.table("go_annotations.txt", sep = "", header=F, stringsAsFactors=FALSE) #id GO annotations
GO <- as.data.frame(GO) #filter GO to subset same count filter as above
GOs <- merge(as.data.frame(GO), gn.keep, by="V1", all=T) #merge GO, ID
GOs[is.na(GOs)] <- "unknown" #list NA as unknown
GOs <- merge(GOs, LENGTH, by="V1") #merge
GOs <- GOs[,1:2] #remove length 
splitted <- strsplit(as.character(GOs$V2.x), ",") #slit into multiple GO ids
GO.terms <- data.frame(v1 = rep.int(GOs$V1, sapply(splitted, length)), v2 = unlist(splitted)) #list all GOs with each assigned transcript
IDs <- row.names(data) #set ID names

#change contig lists to vectors
ALL.vector <-c(t(ALL))
DEG.vector <-c(t(DEG))
ID.vector <- LENGTH$V1
LENGTH.vector <-LENGTH$V2

gene.vector=as.integer(ALL.vector%in%DEG.vector) #Construct new vector with 1 for DEG and 0 for others
names(gene.vector)=ALL.vector #set names
DEG.pwf<-nullp(gene.vector, ID.vector, bias.data=LENGTH.vector) #weight vector by length of gene

#Find enriched GO terms, "selection-unbiased testing for category enrichment amongst differentially expressed (DE) genes for RNA-seq data"
GO.wall<-goseq(DEG.pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)

#How many enriched GO terms do we have
class(GO.wall)
head(GO.wall)
tail(GO.wall)
nrow(GO.wall)

#Find only enriched GO terms that are statistically significant at cutoff 
enriched.GO.05.a<-GO.wall$category[GO.wall$over_represented_pvalue<.05]
enriched.GO.05<-data.frame(enriched.GO.05.a)
colnames(enriched.GO.05) <- c("category")
enriched.GO.05 <- merge(enriched.GO.05, GO.wall, by="category")

MF.20 <- subset(enriched.GO.05, ontology=="MF")
MF.20 <- MF.20[order(-MF.20$numDEInCat),]
CC.20 <- subset(enriched.GO.05, ontology=="CC")
CC.20 <- CC.20[order(-CC.20$numDEInCat),]
BP.20 <- subset(enriched.GO.05, ontology=="BP")
BP.20 <- BP.20[order(-BP.20$numDEInCat),]
write.csv(MF.20 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/MF_Sig_Enriched_GO.05_20.csv")
write.csv(CC.20 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/CC_Sig_Enriched_GO.05_20.csv")
write.csv(BP.20 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/BP_Sig_Enriched_GO.05_20.csv")

##### UNIQUE 00 #####
ALL<-row.names(data)
DEG <- as.character(unique.00$Transcript.ID) #row.names(sig) 
LENGTH <-read.table("Trinity.fasta.seq_lens", sep = "", header=F) #reads in a list of gene lengths
gn.keep <- as.data.frame(gn.keep) #filter length to subset same count filter as above
colnames(gn.keep) <- c("V1") #name columns
lens <- as.data.frame(LENGTH) #set lengths
LENGTH <- merge(as.data.frame(LENGTH), gn.keep, by="V1") #merge lengths and transcripts
GO <- read.table("go_annotations.txt", sep = "", header=F, stringsAsFactors=FALSE) #id GO annotations
GO <- as.data.frame(GO) #filter GO to subset same count filter as above
GOs <- merge(as.data.frame(GO), gn.keep, by="V1", all=T) #merge GO, ID
GOs[is.na(GOs)] <- "unknown" #list NA as unknown
GOs <- merge(GOs, LENGTH, by="V1") #merge
GOs <- GOs[,1:2] #remove length 
splitted <- strsplit(as.character(GOs$V2.x), ",") #slit into multiple GO ids
GO.terms <- data.frame(v1 = rep.int(GOs$V1, sapply(splitted, length)), v2 = unlist(splitted)) #list all GOs with each assigned transcript
IDs <- row.names(data) #set ID names

#change contig lists to vectors
ALL.vector <-c(t(ALL))
DEG.vector <-c(t(DEG))
ID.vector <- LENGTH$V1
LENGTH.vector <-LENGTH$V2

gene.vector=as.integer(ALL.vector%in%DEG.vector) #Construct new vector with 1 for DEG and 0 for others
names(gene.vector)=ALL.vector #set names
DEG.pwf<-nullp(gene.vector, ID.vector, bias.data=LENGTH.vector) #weight vector by length of gene

#Find enriched GO terms, "selection-unbiased testing for category enrichment amongst differentially expressed (DE) genes for RNA-seq data"
GO.wall<-goseq(DEG.pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)

#How many enriched GO terms do we have
class(GO.wall)
head(GO.wall)
tail(GO.wall)
nrow(GO.wall)

#Find only enriched GO terms that are statistically significant at cutoff 
enriched.GO.05.a<-GO.wall$category[GO.wall$over_represented_pvalue<.05]
enriched.GO.05<-data.frame(enriched.GO.05.a)
colnames(enriched.GO.05) <- c("category")
enriched.GO.05 <- merge(enriched.GO.05, GO.wall, by="category")

MF.00 <- subset(enriched.GO.05, ontology=="MF")
MF.00 <- MF.00[order(-MF.00$numDEInCat),]
CC.00 <- subset(enriched.GO.05, ontology=="CC")
CC.00 <- CC.00[order(-CC.00$numDEInCat),]
BP.00 <- subset(enriched.GO.05, ontology=="BP")
BP.00 <- BP.00[order(-BP.00$numDEInCat),]
write.csv(MF.00 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/MF_Sig_Enriched_GO.05_00.csv")
write.csv(CC.00 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/CC_Sig_Enriched_GO.05_00.csv")
write.csv(BP.00 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/BP_Sig_Enriched_GO.05_00.csv")


##### GO SLIM #####

##### BP #####
BP.Ids <- as.character(BP.INT$category)
myCollection <- GOCollection(BP.Ids)
fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl, evidenceCode="TAS")
BP.slim.INT <- goSlim(myCollection, slim, "BP")
BP.slim.INT <-BP.slim.INT[,c(1,3)]

BP.Ids <- as.character(BP.18$category)
myCollection <- GOCollection(BP.Ids)
fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl, evidenceCode="TAS")
BP.slim.18 <- goSlim(myCollection, slim, "BP")
BP.slim.18 <-BP.slim.18[,c(1,3)]

BP.Ids <- as.character(BP.20$category)
myCollection <- GOCollection(BP.Ids)
fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl, evidenceCode="TAS")
BP.slim.20 <- goSlim(myCollection, slim, "BP")
BP.slim.20 <-BP.slim.20[,c(1,3)]

BP.Ids <- as.character(BP.00$category)
myCollection <- GOCollection(BP.Ids)
fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl, evidenceCode="TAS")
BP.slim.00 <- goSlim(myCollection, slim, "BP")
BP.slim.00 <-BP.slim.00[,c(1,3)]

BP.slim.all <- merge(BP.slim.INT, BP.slim.18, by="Term", all = T )
BP.slim.all <- merge(BP.slim.all, BP.slim.20, by="Term", all = T )
colnames(BP.slim.all) <- c("Term", "Intersection", "18:00", "20:00")
BP.slim.all <- merge(BP.slim.all, BP.slim.00, by="Term", all = T )
colnames(BP.slim.all) <- c("Term", "Intersection", "18:00", "20:00", "00:00")
BP.slim.all 
BP.slim.all  <- BP.slim.all[order(BP.slim.all$Intersection),]
BP.slim.all$Term <- gsub("biological_process", "unknown", BP.slim.all$Term)
BP.slim.all
BP.slim.all <- BP.slim.all[rowSums(BP.slim.all[, -1])>3, ]
BP.slim.all <-head(BP.slim.all,-1)


pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/BP_enrichment_SPvsNO.pdf")
par(mfrow=c(1,4))
par(mar=c(3,1,1,0), oma = c(0.1, 10, 0.1, 0.5))
barplot(BP.slim.all$Intersection, horiz = T, col="blue", names.arg=BP.slim.all$Term, las=1, cex.names = 0.6, xlim=c(0,105), main = "ALL")
barplot(BP.slim.all$`18:00`, horiz = T, col="yellow", xlim=c(0,105), main = "18:00")
barplot(BP.slim.all$`20:00`, horiz = T, col="red", xlim=c(0,105), main = "20:00")
barplot(BP.slim.all$`00:00`, horiz = T, col="green", xlim=c(0,105), main = "00:00")
dev.off()

##### MF #####
MF.Ids <- as.character(MF.INT$category)
myCollection <- GOCollection(MF.Ids)
fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl, evidenceCode="TAS")
MF.slim.INT <- goSlim(myCollection, slim, "MF")
MF.slim.INT <-MF.slim.INT[,c(1,3)]

MF.Ids <- as.character(MF.18$category)
myCollection <- GOCollection(MF.Ids)
fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl, evidenceCode="TAS")
MF.slim.18 <- goSlim(myCollection, slim, "MF")
MF.slim.18 <-MF.slim.18[,c(1,3)]

MF.Ids <- as.character(MF.20$category)
myCollection <- GOCollection(MF.Ids)
fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl, evidenceCode="TAS")
MF.slim.20 <- goSlim(myCollection, slim, "MF")
MF.slim.20 <-MF.slim.20[,c(1,3)]

MF.Ids <- as.character(MF.00$category)
myCollection <- GOCollection(MF.Ids)
fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl, evidenceCode="TAS")
MF.slim.00 <- goSlim(myCollection, slim, "MF")
MF.slim.00 <-MF.slim.00[,c(1,3)]

MF.slim.all <- merge(MF.slim.INT, MF.slim.18, by="Term", all = T )
MF.slim.all <- merge(MF.slim.all, MF.slim.20, by="Term", all = T )
colnames(MF.slim.all) <- c("Term", "Intersection", "18:00", "20:00")
MF.slim.all <- merge(MF.slim.all, MF.slim.00, by="Term", all = T )
colnames(MF.slim.all) <- c("Term", "Intersection", "18:00", "20:00", "00:00")
MF.slim.all 
MF.slim.all  <- MF.slim.all [order(MF.slim.all$Intersection),]
MF.slim.all$Term <- gsub("molecular_function", "unknown", MF.slim.all$Term)
MF.slim.all
MF.slim.all <- MF.slim.all[rowSums(MF.slim.all[, -1])>1, ]
MF.slim.all <-head(MF.slim.all,-1)

pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/MF_enrichment_SPvsNO.pdf")
par(mfrow=c(1,4))
par(mar=c(3,1,1,0), oma = c(0.1, 10, 0.1, 0.5))
barplot(MF.slim.all$Intersection, horiz = T, col="blue", names.arg=MF.slim.all$Term, las=1, cex.names = 0.6, xlim=c(0,65), main = "ALL")
barplot(MF.slim.all$`18:00`, horiz = T, col="yellow", xlim=c(0,65), main = "18:00")
barplot(MF.slim.all$`20:00`, horiz = T, col="red", xlim=c(0,65), main = "20:00")
barplot(MF.slim.all$`00:00`, horiz = T, col="green", xlim=c(0,65), main = "00:00")
dev.off()

##### CC #####
CC.Ids <- as.character(CC.INT$category)
myCollection <- GOCollection(CC.Ids)
fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl, evidenceCode="TAS")
CC.slim.INT <- goSlim(myCollection, slim, "CC")
CC.slim.INT <-CC.slim.INT[,c(1,3)]

CC.Ids <- as.character(CC.18$category)
myCollection <- GOCollection(CC.Ids)
fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl, evidenceCode="TAS")
CC.slim.18 <- goSlim(myCollection, slim, "CC")
CC.slim.18 <-CC.slim.18[,c(1,3)]

CC.Ids <- as.character(CC.20$category)
myCollection <- GOCollection(CC.Ids)
fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl, evidenceCode="TAS")
CC.slim.20 <- goSlim(myCollection, slim, "CC")
CC.slim.20 <-CC.slim.20[,c(1,3)]

CC.Ids <- as.character(CC.00$category)
myCollection <- GOCollection(CC.Ids)
fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl, evidenceCode="TAS")
CC.slim.00 <- goSlim(myCollection, slim, "CC")
CC.slim.00 <-CC.slim.00[,c(1,3)]

CC.slim.all <- merge(CC.slim.INT, CC.slim.18, by="Term", all = T )
CC.slim.all <- merge(CC.slim.all, CC.slim.20, by="Term", all = T )
colnames(CC.slim.all) <- c("Term", "Intersection", "18:00", "20:00")
CC.slim.all <- merge(CC.slim.all, CC.slim.00, by="Term", all = T )
colnames(CC.slim.all) <- c("Term", "Intersection", "18:00", "20:00", "00:00")
CC.slim.all 
CC.slim.all  <- CC.slim.all[order(CC.slim.all$Intersection),]
CC.slim.all$Term <- gsub("cellular_component", "unknown", CC.slim.all$Term)
CC.slim.all
CC.slim.all <- CC.slim.all[rowSums(CC.slim.all[, -1])>1, ]
CC.slim.all <-head(CC.slim.all,-1)

pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/CC_enrichment_SPvsNO.pdf")
par(mfrow=c(1,4))
par(mar=c(3,1,1,0), oma = c(0.1, 10, 0.1, 0.5))
barplot(CC.slim.all$Intersection, horiz = T, col="blue", names.arg=CC.slim.all$Term, las=1, cex.names = 0.6, xlim=c(0,20), main = "ALL")
barplot(CC.slim.all$`18:00`, horiz = T, col="yellow", xlim=c(0,30), main = "18:00")
barplot(CC.slim.all$`20:00`, horiz = T, col="red", xlim=c(0,30), main = "20:00")
barplot(CC.slim.all$`00:00`, horiz = T, col="green", xlim=c(0,30), main = "00:00")
dev.off()

##### All Enriched GO Plots #####
#Plot all Enriched GOs
pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/GO_enrichment_SPvsNO.pdf")
par(mfrow=c(2,4))
par(mar=c(3,1,1,0), oma = c(0.1, 15, 0.1, 0.5))
barplot(BP.slim.all$Intersection, horiz = T, col="black", names.arg=BP.slim.all$Term, las=1, cex.names = 0.4, xlim=c(0,25), main = "ALL")
text(labels = 'Biological Process', xpd = NA, srt = 90, x=-26, y=15, cex=1)
barplot(BP.slim.all$`18:00`, horiz = T, col="yellow3", xlim=c(0,25), main = "18:00")
barplot(BP.slim.all$`20:00`, horiz = T, col="steelblue1", xlim=c(0,25), main = "20:00")
barplot(BP.slim.all$`00:00`, horiz = T, col="lightpink2", xlim=c(0,25), main = "00:00")

barplot(MF.slim.all$Intersection, horiz = T, col="black", names.arg=MF.slim.all$Term, las=1, cex.names = 0.6, xlim=c(0,25), main = "ALL")
text(labels = 'Molecular Function', xpd = NA, srt = 90, x=-26, y=8, cex=1)
barplot(MF.slim.all$`18:00`, horiz = T, col="yellow3", xlim=c(0,25), main = "18:00")
barplot(MF.slim.all$`20:00`, horiz = T, col="steelblue1", xlim=c(0,25), main = "20:00")
barplot(MF.slim.all$`00:00`, horiz = T, col="lightpink2", xlim=c(0,25), main = "00:00")

# barplot(CC.slim.all$Intersection, horiz = T, col="black", names.arg=CC.slim.all$Term, las=1, cex.names = 0.6, xlim=c(0,25), main = "ALL")
# text(labels = 'Cellular Component', xpd = NA, srt = 90, x=-26, y=8, cex=1)
# barplot(CC.slim.all$`18:00`, horiz = T, col="yellow3", xlim=c(0,25), main = "18:00")
# barplot(CC.slim.all$`20:00`, horiz = T, col="steelblue1", xlim=c(0,25), main = "20:00")
# barplot(CC.slim.all$`00:00`, horiz = T, col="lightpink2", xlim=c(0,25), main = "00:00")
dev.off()


########## CONTRASTING CORALS BETWEEN EACH TIME POINT #####
resultsNames(DEG.int)
#"Intercept"     "groupNO0.00"   "groupNO18.00"  "groupNO20.00"  "groupYES0.00"  "groupYES18.00" "groupYES20.00"

# Time point Yes No contrasts
###### NO = DEG between different time points for NONspawning corals #####
DEG.NO_18_20 <- results(DEG.int, contrast=c("group", "NO18:00", "NO20:00")) #contrasts between SP NSP at each time
DEG.NO_18_20 #view results

sig.num.NO_18_20 <- sum(DEG.NO_18_20$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
sig.NO_18_20 <- subset(DEG.NO_18_20, padj<0.05,) #identify signficant pvalues with 1%FDR
sig.list.NO_18_20 <- data[which(rownames(data) %in% rownames(sig.NO_18_20)),] #subset list of sig transcripts from original count data
sig.list.NO_18_20

DEG.NO_18_00 <- results(DEG.int, contrast=c("group", "NO18:00", "NO0:00")) #contrasts between SP NSP at each time
DEG.NO_18_00 #view results

sig.num.NO_18_00 <- sum(DEG.NO_18_00$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
sig.NO_18_00 <- subset(DEG.NO_18_00, padj<0.05,) #identify signficant pvalues with 1%FDR
sig.list.NO_18_00 <- data[which(rownames(data) %in% rownames(sig.NO_18_00)),] #subset list of sig transcripts from original count data
sig.list.NO_18_00

DEG.NO_20_00 <- results(DEG.int, contrast=c("group", "NO20:00", "NO0:00")) #contrasts between SP NSP at each time
DEG.NO_20_00 #view results

sig.num.NO_20_00 <- sum(DEG.NO_20_00$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
sig.NO_20_00 <- subset(DEG.NO_20_00, padj<0.05,) #identify signficant pvalues with 1%FDR
sig.list.NO_20_00 <- data[which(rownames(data) %in% rownames(sig.NO_20_00)),] #subset list of sig transcripts from original count data
sig.list.NO_20_00

times_18_20 <- merge(intersections, norm.sig.counts, by="Transcript.ID")

###### Prepare Data for Venn Diagram #####
NO_18_20 <- data.frame(rownames(sig.list.NO_18_20)) #identify DE transcript list from NO betwee 18:00 and 20
NO_18_20$NO_18_20.Count <- 1 #assign all transcripts count of 1
colnames(NO_18_20) <- c("Transcript.ID", "NO_18_20.Count") #name columns

NO_18_00 <- data.frame(rownames(sig.list.NO_18_00)) #identify DE transcript list from NO betwee 18:00 and 00
NO_18_00$NO_18_00.Count <- 1 #assign all transcripts count of 1
colnames(NO_18_00) <- c("Transcript.ID", "NO_18_00.Count") #name columns

NO_20_00 <- data.frame(rownames(sig.list.NO_20_00)) #identify DE transcript list from NO betwee 20:00 and 00
NO_20_00$NO_20_00.Count <- 1 #assign all transcripts count of 1
colnames(NO_20_00) <- c("Transcript.ID", "NO_20_00.Count") #name columns

L_18_20 <- merge(NO_18_20, NO_18_00, by="Transcript.ID", all=T ) #merge lists sequentially with no removal
L_18_20_00 <- merge(L_18_20, NO_20_00, by="Transcript.ID", all=T ) #merge lists sequentially with no removal
L_18_20_00$NO_18_20.Count[is.na(L_18_20_00$NO_18_20.Count)] <- 0 #assign NA=0 count i.e., not DE
L_18_20_00$NO_18_00.Count[is.na(L_18_20_00$NO_18_00.Count)] <- 0 #assign NA=0 count i.e., not DE
L_18_20_00$NO_20_00.Count[is.na(L_18_20_00$NO_20_00.Count)] <- 0 #assign NA=0 count i.e., not DE

# Plot Venn Diagram of shared DE between time points
vennData.NO<-L_18_20_00[,2:4] #set venn data to DE counts only
a <- vennCounts(vennData.NO) # compute classification counts
colnames(a) <- c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00', "Counts") #set catagory names
jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Venn.DEG.No.Spawn.jpg")
vennDiagram(a, main='DEG between Time Points for Nonspawning Corals') #draw venn diagram
dev.off()

colnames(vennData.NO) <- c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00') #set catagory names
jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/NOSP_Timepoints_DEG_Intersections.jpg")
upset(vennData.NO, sets = c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00'), order.by = "degree")
dev.off()

#####

#Look for match between all lists
NO.shared.all <- Reduce(intersect, list(rownames(sig.list.NO_18_20),rownames(sig.list.NO_18_00),rownames(sig.list.NO_20_00))) #idenitfy shared DE transcripts between each time point 
NO.shared.all #view results
NO.shared.all <- as.data.frame(NO.shared.all)
colnames(NO.shared.all)[1] <- "Transcript.ID"
write.csv(NO.shared.all, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/DEG_shared_NOSpawn_18_20_00.csv")

#DEG between spawning and not spawning corals Found in all time points
intersections.NO <- subset(L_18_20_00, NO_18_20.Count==1 & NO_18_00.Count==1 & NO_20_00.Count==1) #list DEG found in all time points

#DEG between not spawning corals Unique to each time point
unique.NO_18_20 <- subset(L_18_20_00, NO_18_20.Count==1 & NO_18_00.Count==0 & NO_20_00.Count==0) #list genes unique to 18:00
unique.NO_18_00 <- subset(L_18_20_00, NO_18_20.Count==0 & NO_18_00.Count==1 & NO_20_00.Count==0) #list genes unique to 20:00
unique.NO_20_00 <- subset(L_18_20_00, NO_18_20.Count==0 & NO_18_00.Count==0 & NO_20_00.Count==1) #list genes unique to 00:00

#DEG between spawning and not spawning corals Shared unique between pairs
shared.NO_18_20 <- subset(L_18_20_00, NO_18_20.Count==1 & NO_18_00.Count==1 & NO_20_00.Count==0) #list genes common to 18:00 and 20:00
shared.NO_18_00 <- subset(L_18_20_00, NO_18_20.Count==1 & NO_18_00.Count==0 & NO_20_00.Count==1) #list genes common to 18:00 and 00:00
shared.NO_20_00 <- subset(L_18_20_00, NO_18_20.Count==0 & NO_18_00.Count==1 & NO_20_00.Count==1) #list genes common to 20:00 and 00:00

##### Combine with DE list with annotation data #####

annot.intersect_18_20_00 <- merge(NO.shared.all, annot, by="Transcript.ID")
write.csv(annot.intersect_18_20_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_shared_DEG_18_20_00_NO.csv")

annot.intersect_18_20 <- merge(unique.NO_18_20, annot, by="Transcript.ID")
unique(annot.intersect_18_20$Transcript.ID)
write.csv(annot.intersect_18_20, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_unique_DEG_18_20_NO.csv")

annot.intersect_18_00 <- merge(unique.NO_18_00, annot, by="Transcript.ID")
unique(annot.intersect_18_00$Transcript.ID)
write.csv(annot.intersect_18_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_unique_DEG_18_00_NO.csv")

annot.intersect_20_00 <- merge(unique.NO_20_00, annot, by="Transcript.ID")
unique(annot.intersect_20_00$Transcript.ID)
write.csv(annot.intersect_20_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_unique_DEG_20_00_NO.csv")

##### Timepoint Intersections #####

# YES = DEG between different time points for spawning corals
DEG.YES_18_20 <- results(DEG.int, contrast=c("group", "YES18:00", "YES20:00")) #contrasts between SP NSP at each time
DEG.YES_18_20 #view results

sig.num.YES_18_20 <- sum(DEG.YES_18_20$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
sig.YES_18_20 <- subset(DEG.YES_18_20, padj<0.05,) #identify signficant pvalues with 1%FDR
sig.list.YES_18_20 <- data[which(rownames(data) %in% rownames(sig.YES_18_20)),] #subset list of sig transcripts from original count data
sig.list.YES_18_20

DEG.YES_18_00 <- results(DEG.int, contrast=c("group", "YES18:00", "YES0:00")) #contrasts between SP NSP at each time
DEG.YES_18_00 #view results

sig.num.YES_18_00 <- sum(DEG.YES_18_00$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
sig.YES_18_00 <- subset(DEG.YES_18_00, padj<0.05,) #identify signficant pvalues with 1%FDR
sig.list.YES_18_00 <- data[which(rownames(data) %in% rownames(sig.YES_18_00)),] #subset list of sig transcripts from original count data
sig.list.YES_18_00

DEG.YES_20_00 <- results(DEG.int, contrast=c("group", "YES20:00", "YES0:00")) #contrasts between SP NSP at each time
DEG.YES_20_00 #view results

sig.num.YES_20_00 <- sum(DEG.YES_20_00$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
sig.YES_20_00 <- subset(DEG.YES_20_00, padj<0.05,) #identify signficant pvalues with 1%FDR
sig.list.YES_20_00 <- data[which(rownames(data) %in% rownames(sig.YES_20_00)),] #subset list of sig transcripts from original count data
sig.list.YES_20_00


##### Prepare Data for Venn Diagram #####
YES_18_20 <- data.frame(rownames(sig.list.YES_18_20)) #identify DE transcript list from YES betwee 18:00 and 20
YES_18_20$YES_18_20.Count <- 1 #assign all transcripts count of 1
colnames(YES_18_20) <- c("Transcript.ID", "YES_18_20.Count") #name columns

YES_18_00 <- data.frame(rownames(sig.list.YES_18_00)) #identify DE transcript list from YES betwee 18:00 and 00
YES_18_00$YES_18_00.Count <- 1 #assign all transcripts count of 1
colnames(YES_18_00) <- c("Transcript.ID", "YES_18_00.Count") #name columns

YES_20_00 <- data.frame(rownames(sig.list.YES_20_00)) #identify DE transcript list from YES betwee 20:00 and 00
YES_20_00$YES_20_00.Count <- 1 #assign all transcripts count of 1
colnames(YES_20_00) <- c("Transcript.ID", "YES_20_00.Count") #name columns

Y_18_20 <- merge(YES_18_20, YES_18_00, by="Transcript.ID", all=T ) #merge lists sequentially with YES removal
Y_18_20_00 <- merge(Y_18_20, YES_20_00, by="Transcript.ID", all=T ) #merge lists sequentially with YES removal
Y_18_20_00$YES_18_20.Count[is.na(Y_18_20_00$YES_18_20.Count)] <- 0 #assign NA=0 count i.e., YESt DE
Y_18_20_00$YES_18_00.Count[is.na(Y_18_20_00$YES_18_00.Count)] <- 0 #assign NA=0 count i.e., YESt DE
Y_18_20_00$YES_20_00.Count[is.na(Y_18_20_00$YES_20_00.Count)] <- 0 #assign NA=0 count i.e., YESt DE

# Plot Venn Diagram of shared DE between time points
vennData.YES<-Y_18_20_00[,2:4] #set venn data to DE counts only
a <- vennCounts(vennData.YES) # compute classification counts
#a<- a[rowSums(a[, -1])>0, ]
colnames(a) <- c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00', "Counts") #set catagory names
#a<- a[,-4]
jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Venn.DEG.YES.Spawn.jpg")
vennDiagram(a, main='DEG between Time Points for YES spawning Corals') #draw venn diagram
dev.off()

colnames(vennData.YES) <- c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00') #set catagory names
jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/YESSP_Timepoints_DEG_Intersections.jpg")
upset(vennData.YES, sets = c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00'), order.by = "degree")
dev.off()


# jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Timepoints_DEG_Intersections.jpg")
# layout(mat=matrix(c(1,2),nrow=1,ncol=2,byrow=FALSE))
# upset(vennData.NO, sets = c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00'), order.by = "degree")
# upset(vennData.YES, sets = c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00'), order.by = "degree")
# dev.off()


#####
#Look for match between all lists
YES.shared.all <- Reduce(intersect, list(rownames(sig.list.YES_18_20),rownames(sig.list.YES_18_00),rownames(sig.list.YES_20_00))) #idenitfy shared DE transcripts between each time point 
YES.shared.all #view results
YES.shared.all <- as.data.frame(YES.shared.all)
colnames(YES.shared.all)[1] <- "Transcript.ID"
write.csv(YES.shared.all, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/DEG_shared_YESSpawn_18_20_00.csv")

#DEG between spawning and YESt spawning corals Found in all time points
intersections.YES <- subset(Y_18_20_00, YES_18_20.Count==1 & YES_18_00.Count==1 & YES_20_00.Count==1) #list DEG found in all time points

#DEG between YESt spawning corals Unique to each time point
unique.YES_18_20 <- subset(Y_18_20_00, YES_18_20.Count==1 & YES_18_00.Count==0 & YES_20_00.Count==0) #list genes unique to 18:00
unique.YES_18_00 <- subset(Y_18_20_00, YES_18_20.Count==0 & YES_18_00.Count==1 & YES_20_00.Count==0) #list genes unique to 20:00
unique.YES_20_00 <- subset(Y_18_20_00, YES_18_20.Count==0 & YES_18_00.Count==0 & YES_20_00.Count==1) #list genes unique to 00:00

#DEG between spawning and YESt spawning corals Shared unique between pairs
shared.YES_18_20 <- subset(Y_18_20_00, YES_18_20.Count==1 & YES_18_00.Count==1 & YES_20_00.Count==0) #list genes common to 18:00 and 20:00
shared.YES_18_00 <- subset(Y_18_20_00, YES_18_20.Count==1 & YES_18_00.Count==0 & YES_20_00.Count==1) #list genes common to 18:00 and 00:00
shared.YES_20_00 <- subset(Y_18_20_00, YES_18_20.Count==0 & YES_18_00.Count==1 & YES_20_00.Count==1) #list genes common to 20:00 and 00:00

##### Combine with DE list with annotation data #####

annot.intersect_18_20_00 <- merge(YES.shared.all, annot, by="Transcript.ID")
write.csv(annot.intersect_18_20_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_shared_DEG_18_20_00_YES.csv")

annot.intersect_18_20 <- merge(unique.YES_18_20, annot, by="Transcript.ID")
unique(annot.intersect_18_20$Transcript.ID)
write.csv(annot.intersect_18_20, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_shared_DEG_18_20_YES.csv")

annot.intersect_18_00 <- merge(unique.YES_18_00, annot, by="Transcript.ID")
unique(annot.intersect_18_00$Transcript.ID)
write.csv(annot.intersect_18_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_shared_DEG_18_00_YES.csv")

annot.intersect_20_00 <- merge(unique.YES_20_00, annot, by="Transcript.ID")
unique(annot.intersect_20_00$Transcript.ID)
write.csv(annot.intersect_20_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_shared_DEG_20_00_NO.csv")





##### KEGG #####
# 
# #KEGG
# annot$KO <- sapply(strsplit(as.character(annot$Kegg), split="\\`"), "[", 2)
# annot$KO <- gsub(".*:","",annot$KO)
# annots <- annot$KO[!is.na(annot$KO),]
# KOs <- annot$KO
# KOs <-KOs[!is.na(KOs)]
# write.csv(KOs, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/All.KO.List.csv")
# 
# write.csv(anot.inters$KO, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/KO.List.csv")
# 
# anot.inters <- anot.inters[!is.na(anot.inters$KO),]
# #rownames(anot.inters) <- anot.inters$KO
# #anot.inters <- anot.inters[!duplicated(anot.inters$KG),]
# 
# keggcode <- "hsa"
# all.DEG.annot <- pathview(gene.data = anot.inters[, 38], gene.idtype = "entrez", pathway.id = "04713", species = keggcode, out.suffix = "intersection")
# 
# #example
# #pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",species = "hsa", out.suffix = "gse16873")
