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
storage.mode(counts.5x) = "integer" #store counts data as integer
sample.info <- read.csv(file="sample_description.csv", header=T, sep=",", row.names=1) #load sample info
sample.info$group <- factor(paste0(sample.info$Spawn, sample.info$Time)) #merge condition and time into group
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
DEG.int <- DESeq(data) #run differential expression test by group using the wald test (Does this filter out the low read counts, or is this DESeq2?)
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
plot(PC.info$PC1, PC.info$PC2, xlim=c(-50,50), ylim=c(-16, 10), xlab="PC1 86%", ylab="PC2 5%", pch = c(15, 16, 17, 18)[as.numeric(sample.info$Time)], col=c("black", "gray")[sample.info$Spawn], cex=1.3)
legend(x="topleft", 
       bty="n",
       legend = c("00:00", "18:00", "20:00"),
       pch = c(15, 16, 17, 18),
       col=c("black", "black", "black"),
       cex=1)
dev.off()

topVarGenes <- head(order(rowVars(assay(rsig)),decreasing=TRUE),sig.num) #sort by decreasing sig
mat <- assay(rsig)[ topVarGenes, ] #make an expression object
mat <- mat - rowMeans(mat) #difference in expression compared to average across all samples
df <- as.data.frame(colData(rsig)[,c("Spawn","Time")]) #make dataframe
jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/DEG_Heatmap.jpg")
pheatmap(mat, annotation_col=df, clustering_method = "average", 
         clustering_distance_rows="euclidean", show_rownames =F, 
         show_colnames =F) #plot heatmap of all DEG by group
dev.off()

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
colnames(vennData) <- c("18:00",'20:00', '00:00') #set catagory names
jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/SP_NSP_DEG_Intersections.jpg")
upset(vennData, sets = c('18:00','20:00', '00:00'), sets.bar.color = c("lightpink2", "yellow3", "steelblue1"), main.bar.color = c("blue", "black", "black", "black", "yellow", "red", "green"), order.by = "degree")
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
annot <- read.csv("trinotate_annotation_report_wo_seq.csv", header=T, sep=",")
colnames(annot)[2] <- "Transcript.ID"

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

#Analysis for All DE (both Host and Sym together)

#Would like to make this a for loop to run through 4 files intersection, 18, 20, 00

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
BP.slim.all <- BP.slim.all[rowSums(BP.slim.all[, -1])>1, ]

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

pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/CC_enrichment_SPvsNO.pdf")
par(mfrow=c(1,4))
par(mar=c(3,1,1,0), oma = c(0.1, 10, 0.1, 0.5))
barplot(CC.slim.all$Intersection, horiz = T, col="blue", names.arg=CC.slim.all$Term, las=1, cex.names = 0.6, xlim=c(0,20), main = "ALL")
barplot(CC.slim.all$`18:00`, horiz = T, col="yellow", xlim=c(0,30), main = "18:00")
barplot(CC.slim.all$`20:00`, horiz = T, col="red", xlim=c(0,30), main = "20:00")
barplot(CC.slim.all$`00:00`, horiz = T, col="green", xlim=c(0,30), main = "00:00")
dev.off()

#Plot all Enriched GOs
pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/GO_enrichment_SPvsNO.pdf")
par(mfrow=c(3,4))
par(mar=c(3,1,1,0), oma = c(0.1, 10, 0.1, 0.5))
barplot(BP.slim.all$Intersection, horiz = T, col="blue", names.arg=BP.slim.all$Term, las=1, cex.names = 0.4, xlim=c(0,105), main = "ALL")
barplot(BP.slim.all$`18:00`, horiz = T, col="yellow", xlim=c(0,105), main = "18:00")
barplot(BP.slim.all$`20:00`, horiz = T, col="red", xlim=c(0,105), main = "20:00")
barplot(BP.slim.all$`00:00`, horiz = T, col="green", xlim=c(0,105), main = "00:00")

barplot(MF.slim.all$Intersection, horiz = T, col="blue", names.arg=MF.slim.all$Term, las=1, cex.names = 0.6, xlim=c(0,105), main = "ALL")
barplot(MF.slim.all$`18:00`, horiz = T, col="yellow", xlim=c(0,105), main = "18:00")
barplot(MF.slim.all$`20:00`, horiz = T, col="red", xlim=c(0,105), main = "20:00")
barplot(MF.slim.all$`00:00`, horiz = T, col="green", xlim=c(0,105), main = "00:00")

barplot(CC.slim.all$Intersection, horiz = T, col="blue", names.arg=CC.slim.all$Term, las=1, cex.names = 0.6, xlim=c(0,105), main = "ALL")
barplot(CC.slim.all$`18:00`, horiz = T, col="yellow", xlim=c(0,105), main = "18:00")
barplot(CC.slim.all$`20:00`, horiz = T, col="red", xlim=c(0,105), main = "20:00")
barplot(CC.slim.all$`00:00`, horiz = T, col="green", xlim=c(0,105), main = "00:00")
dev.off()

########## Timing contrasts #####
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

L_18_20 <- merge(YES_18_20, YES_18_00, by="Transcript.ID", all=T ) #merge lists sequentially with YES removal
L_18_20_00 <- merge(L_18_20, YES_20_00, by="Transcript.ID", all=T ) #merge lists sequentially with YES removal
L_18_20_00$YES_18_20.Count[is.na(L_18_20_00$YES_18_20.Count)] <- 0 #assign NA=0 count i.e., YESt DE
L_18_20_00$YES_18_00.Count[is.na(L_18_20_00$YES_18_00.Count)] <- 0 #assign NA=0 count i.e., YESt DE
L_18_20_00$YES_20_00.Count[is.na(L_18_20_00$YES_20_00.Count)] <- 0 #assign NA=0 count i.e., YESt DE

# Plot Venn Diagram of shared DE between time points
vennData.YES<-L_18_20_00[,2:4] #set venn data to DE counts only
a <- vennCounts(vennData.YES) # compute classification counts
a<- a[rowSums(a[, -1])>0, ]
colnames(a) <- c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00', "Counts") #set catagory names
a<- a[,-4]
jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Venn.DEG.YES.Spawn.jpg")
vennDiagram(a, main='DEG between Time Points for YES spawning Corals') #draw venn diagram
dev.off()

#####
#Look for match between all lists
YES.shared.all <- Reduce(intersect, list(rownames(sig.list.YES_18_20),rownames(sig.list.YES_18_00),rownames(sig.list.YES_20_00))) #idenitfy shared DE transcripts between each time point 
YES.shared.all #view results
YES.shared.all <- as.data.frame(YES.shared.all)
colnames(YES.shared.all)[1] <- "Transcript.ID"
write.csv(YES.shared.all, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/DEG_shared_YESSpawn_18_20_00.csv")

#DEG between spawning and YESt spawning corals Found in all time points
intersections.YES <- subset(L_18_20_00, YES_18_20.Count==1 & YES_18_00.Count==1 & YES_20_00.Count==1) #list DEG found in all time points

#DEG between YESt spawning corals Unique to each time point
unique.YES_18_20 <- subset(L_18_20_00, YES_18_20.Count==1 & YES_18_00.Count==0 & YES_20_00.Count==0) #list genes unique to 18:00
unique.YES_18_00 <- subset(L_18_20_00, YES_18_20.Count==0 & YES_18_00.Count==1 & YES_20_00.Count==0) #list genes unique to 20:00
unique.YES_20_00 <- subset(L_18_20_00, YES_18_20.Count==0 & YES_18_00.Count==0 & YES_20_00.Count==1) #list genes unique to 00:00

#DEG between spawning and YESt spawning corals Shared unique between pairs
shared.YES_18_20 <- subset(L_18_20_00, YES_18_20.Count==1 & YES_18_00.Count==1 & YES_20_00.Count==0) #list genes common to 18:00 and 20:00
shared.YES_18_00 <- subset(L_18_20_00, YES_18_20.Count==1 & YES_18_00.Count==0 & YES_20_00.Count==1) #list genes common to 18:00 and 00:00
shared.YES_20_00 <- subset(L_18_20_00, YES_18_20.Count==0 & YES_18_00.Count==1 & YES_20_00.Count==1) #list genes common to 20:00 and 00:00

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

#KEGG

anot.inters <- anot.inters[!is.na(anot.inters$KO),]
#rownames(anot.inters) <- anot.inters$KO
#anot.inters <- anot.inters[!duplicated(anot.inters$KG),]

keggcode <- "hsa"
all.DEG.annot <- pathview(gene.data = anot.inters[, 38], gene.idtype = "entrez", pathway.id = "04713", species = keggcode, out.suffix = "intersection")

#example
#pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",species = "hsa", out.suffix = "gse16873")
