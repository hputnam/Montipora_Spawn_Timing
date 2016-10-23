library(cluster)
library(Biobase)
library(qvalue)

# try to reuse earlier-loaded data if possible
if (file.exists("diffExpr.P0.05_C0.matrix.RData")) {
    print('RESTORING DATA FROM EARLIER ANALYSIS')
    load("diffExpr.P0.05_C0.matrix.RData")
} else {
    print('Reading matrix file.')
    primary_data = read.table("diffExpr.P0.05_C0.matrix", header=T, com='', sep="\t", row.names=1, check.names=F)
    primary_data = as.matrix(primary_data)
}
source("/usr/local/opt/trinityrnaseq_r20140717/Analysis/DifferentialExpression/R/heatmap.3.R")
source("/usr/local/opt/trinityrnaseq_r20140717/Analysis/DifferentialExpression/R/misc_rnaseq_funcs.R")
source("/usr/local/opt/trinityrnaseq_r20140717/Analysis/DifferentialExpression/R/pairs3.R")
data = primary_data
samples_data = read.table("/home/hputnam/Mcap_Spawn/Refs/sample_description.txt", header=F, check.names=F)
sample_types = as.character(unique(samples_data[,1]))
data = data[, colnames(data) %in% samples_data[,2], drop=F ]
nsamples = length(sample_types)
sample_colors = rainbow(nsamples)
names(sample_colors) = sample_types
sample_type_list = list()
for (i in 1:nsamples) {
    samples_want = samples_data[samples_data[,1]==sample_types[i], 2]
    sample_type_list[[sample_types[i]]] = as.vector(samples_want)
}
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
# reorder according to sample type.
tmp_sample_reordering = order(sample_factoring)
data = data[,tmp_sample_reordering,drop=F]
sample_factoring = sample_factoring[tmp_sample_reordering]
data = log2(data+1)
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
sampleAnnotations = matrix(ncol=ncol(data),nrow=nsamples)
for (i in 1:nsamples) {
  sampleAnnotations[i,] = colnames(data) %in% sample_type_list[[sample_types[i]]]
}
sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
rownames(sampleAnnotations) = as.vector(sample_types)
colnames(sampleAnnotations) = colnames(data)
data = as.matrix(data) # convert to matrix
write.table(data, file="diffExpr.P0.05_C0.matrix.log2.dat", quote=F, sep='	');
sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
sample_dist = dist(t(data), method='euclidean')
hc_samples = hclust(sample_dist, method='complete')
pdf("diffExpr.P0.05_C0.matrix.log2.sample_cor_matrix.pdf")
heatmap.3(sample_cor, dendrogram='both', Rowv=as.dendrogram(hc_samples), Colv=as.dendrogram(hc_samples), col = greenred(75), scale='none', symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, margins=c(10,10), cexCol=1, cexRow=1, cex.main=0.75, main=paste("sample correlation matrix
", "diffExpr.P0.05_C0.matrix.log2") , ColSideColors=sampleAnnotations, RowSideColors=t(sampleAnnotations))
dev.off()
gene_cor = NULL
gene_dist = dist(data, method='euclidean')
if (nrow(data) <= 1) { message('Too few genes to generate heatmap'); quit(status=0); }
hc_genes = hclust(gene_dist, method='complete')
myheatcol = greenred(75)
data = t(scale(t(data), scale=F)) # center rows, mean substracted
write.table(data, file="diffExpr.P0.05_C0.matrix.log2.centered.dat", quote=F, sep='	');
# redo the sample clustering according to the centered distance values.
sample_dist = dist(t(data), method='euclidean')
# redo the sample clustering
hc_samples = hclust(sample_dist, method='complete')
# redo the gene clustering according to centered distance values.
gene_dist = dist(data, method='euclidean')
# redo the gene clustering
hc_genes = hclust(gene_dist, method='complete')
heatmap_data = data
pdf("diffExpr.P0.05_C0.matrix.log2.centered.genes_vs_samples_heatmap.pdf")
heatmap.3(heatmap_data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=1, margins=c(10,10), cex.main=0.75, main=paste("samples vs. features
", "diffExpr.P0.05_C0.matrix.log2.centered" ) , ColSideColors=sampleAnnotations)
dev.off()
save(list=ls(all=TRUE), file="diffExpr.P0.05_C0.matrix.RData")
