

## Load libraries
library(DESeq2)
library(dplyr)
library(pheatmap)
library(tidyverse)
library(ggplot2)
library(calibrate)
library(cluster)

# READ COUNT DATA FROM FEATURECOUNTS

countData<-read.csv("LambAllSamples.featureCounts.Rmatrix.txt", sep="\t", header=T,check.names=F) # This file has gene ids and sample files only
head(countData)
countData$Chr<-NULL # run this only if chr column is present
countData<-countData[ , !names(countData) %in% c("Low7085.bam","Medium7073.bam")]# Remove 7085 and 7073 as they had poor mapping rates
dim(countData)
colnames(countData)<-gsub(".bam","",colnames(countData))# Remove the .bam from the column names
head(countData)
orig_names <- names(countData) # keep a back-up copy of the original names
geneID <- countData$Geneid# Convert count data to a matrix of appropriate form that DEseq2 can read
countData <- as.matrix(countData[ , -1]) # removing first column geneID from the table
# make sure the rownames are gene ids and first column onwards should be samples. any other columns then remove.otherwise deseq2 returns assys negative error in next steps
sampleIndex <- colnames(countData)
countData <- as.matrix(countData[,sampleIndex])
rownames(countData) <- geneID
head(countData)

# check library size (total mapped reads)
colSums(countData[,2:ncol(countData)])

# check the number of detected genes
apply(countData[,2:ncol(countData)], 2, function(c)sum(c!=0))

# filter genes
filter_genes <- apply(
    countData[,2:ncol(countData)],
    1,
    function(x) length(x[x > 2]) >= 2
)

mx_filterGenes <- countData[filter_genes,]
head(mx_filterGenes[,1:5])
dim(mx_filterGenes)

# READ THE METADATA TABLE AND TRANSFORM IT IN A WAY DESEQ2 ACCEPTS

# Convert sample variable mappings to an appropriate form that DESeq2 can read
metaData <-read.csv("metadata.csv",sep=",",header=T)
dim(metaData)
head(metaData)
rownames(metaData) <- metaData$ID
metaData$ID <- factor(metaData$ID)
head(metaData)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
countData <- countData[,unique(rownames(metaData))]
all(colnames(countData) == rownames(metaData)) # should return true if both tables have same info

# CREATE THE DESEQ2DATASET OBJECT

# design formula is used to estimate the dispersions and to estimate the log2 fold changes of the model.
# we are interested in the sampletype column which has four factor levels, which tells DESeq2 that for each gene we want to evaluate gene expression change with respect to these different levels.
deseq2Data <- DESeqDataSetFromMatrix(countData=countData, colData=metaData, design= ~sampletype)
dim(deseq2Data) #35117 genes

keep <- rowSums(counts(deseq2Data)) >= 4
deseq2Data <- deseq2Data[keep,]

# remove rows with only 0
dds <- dds[rowSums(assay(dds)) > 0, ]
dim(dds)

# normlize using rlog mathod
norm <- rlog(deseq2Data,blind=FALSE)
norm_matrix <- assay(norm)
norm_df <- data.frame(Gene=rownames(norm_matrix), norm_matrix)
write.table(norm_df, "DESeq2.rlog.data", quote=F, row.names = FALSE,sep="\t")

# It is also advisable to investigate any systematic bias in the sequencing data, such as whether one sample has been sequenced more deeply than others. 
# One can extract size factors using the sizeFactors() function. Usually these size factors should vary around 1, indicating comparable sequencing depth.
#sizeFactors(dds)

# compare sizefactors to seq depth
pdf("sizeFactors_deseqdata.pdf", width=15)
tibble(sample = colnames(dds), SF = sizeFactors(dds)) %>%
  ggplot(aes(x = sample, y = SF)) +
  geom_col()
dev.off()

pdf("seqdepth_deseqdata.pdf", width=15)
tibble(sample = colnames(dds), SeqDepth = colSums(assay(dds))) %>%
  ggplot(aes(x = sample, y = SeqDepth)) +
  geom_col()
dev.off()

#You can use DESeq-specific functions to access the different slots and retrieve information, if you wish.
#counts(deseq2Data)# to view original count matrix

# GENERTAE NORMALIZED COUNTS
# DESeq2 has a function estimateSizefactors which generate the size factors to perform median of ratios normalization.
# For an RNA seq analysis, this is automatically done by the DESeq() function.
#deseq2Data <- estimateSizeFactors(deseq2Data)
#sizeFactors(deseq2Data)# take a look at the normalization factor applied to each sample
#normalized_counts <- counts(deseq2Data, normalized=TRUE)# to retrieve the normalized counts matrix
#write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)# save this normalized data matrix to file for later use

##### NOTE ####
#NOTE: DESeq2 doesn't actually use normalized counts, rather it uses the raw counts and models the normalization inside the Generalized Linear Model (GLM). 
#These normalized counts will be useful for downstream visualization of results, but cannot be used as input to DESeq2 or any other tools that peform differential 
#expression analysis which use the negative binomial model.

# DIFFERENTIAL EXPRESSION ANALYSIS WITH DESEQ2
# Null hypothesis: no diff expn between groups, ie log2FC=0. A statistical test, the Wald test, will determine whether the data provides sufficient evidence to conclude that this value is really different from zero.
# Pvalue generated from wald test result in many false positives, hence needs to be adjusted which DESEq2 adjust using BH method and generate adj pvalues.
deseq2Data <- DESeq(deseq2Data, betaPrior = FALSE)#The betaPrior = T option tells DESeq2 to squeeze the log Fold Changes of lowly expressed genes toward 0.
# betaPrior=F is the default. this F option should be used if you want to perfrom LFC cutoff downstram
# The above code does
# 1. sequencing depth normalization between the samples (estimateSizeFactors)
# 2. gene-wise dispersion estimates across all samples (estimateDispersions)
# 3. fits a negative binomial GLM and applies Wald statistics to each gene (nbinomWaldTest)

resultsNames(deseq2Data) # tells you which types of values can be extracted with results()

# extract results for different conditions. 
# The default approach of DESeq uses a null hypothesis of 0 LFC. In this case, each gene is compared across treatment for this null hypothesis, 
# and considered significantly different if smaller than the FDR. When you use a fold change cut off (of lets say 1.5) or LFC 0.584, your null 
# hypothesis is that genes between the experiment differ more than 1.5 fold change, and you will be testing that hypothesis with FDR cut off.
res <- results(deseq2Data,
               name = "sampletype_Low_vs_Control")

#The filtering threshold that has been used to filter low count genes can be extracted from the results metadata.
metadata(res)$filterThreshold
res_tbl <- as_tibble(res, rownames = "GENE")
# histogram of p-values
pdf("p-value_histogram.pdf")
hist(res_tbl$pvalue)
dev.off()
# We observe that most of the genes that have a pvalue around 0.8 (corresponding to genes belonging to the second peak) have a pvalue, but no padjusted value. 
# These are hence genes that were filtered out by the independent filtering process. This means that DESeq2 calculates pvalue for every gene,
# but the pvalue adjustement is computed later only, on genes that pass the independent filtering threshold. 
#If we plot again the pvalues histogram using filtered genes only, then we observe that the pvalues behaviour is now much nicer!

pfilt<- res_tbl %>%
    filter(baseMean > metadata(res)$filterThreshold) %>%
    pull(pvalue)
pdf("p-value_filtered_histogram.pdf")
hist(pfilt)
dev.off()

## MA plot
# a scatter plot of log2FC versus the mean of normalised counts. Genes with a padjusted value lower than 0.05 are colored. 
# The plot highlights the fact that genes with low read counts show substantially higher variability than highly expressed genes, resulting in a strong log2FC even though are likely not really differentially expressed. 
# In the MA-plot, we hope to observe some genes located in the upper/lower right quadrant of the plots (the most interesting candidates).
pdf("MAplot.pdf")
plotMA(res,alpha=0.1,colNonSig = "black",
       colSig = "red",
       colLine = "blue",)
dev.off()

# Extracting and Filetring results
LC_res <- results(deseq2Data, independentFiltering = TRUE,contrast = c("sampletype","Low","Control"), alpha=0.1)#build results table
summary(LC_res)
LC_shrunken <- lfcShrink(deseq2Data, coef="sampletype_Low_vs_Control", type="apeglm", res=LC_res)
summary(LC_shrunken)# will summarize the results using the alpha threshold: FDR < 0.1 (padj/FDR is used even though the output says p-value < 0.05).
table(LC_shrunken$padj < 0.1)
table(LC_shrunken$padj < 0.1 & LC_shrunken$log2FoldChange > 0.584) #upreg
table(LC_shrunken$padj < 0.1 & LC_shrunken$log2FoldChange < -0.584) #downreg
# writing sig genes to file
LC_resSigUP = LC_shrunken[ which(LC_shrunken$padj < 0.1 & LC_shrunken$log2FoldChange > 0.584), ]
LC_resSigDown = LC_shrunken[ which(LC_shrunken$padj < 0.1 & LC_shrunken$log2FoldChange < -0.584), ]
LC_resSig = rbind(LC_resSigUP, LC_resSigDown)
write.csv(LC_resSig,"LC.DEGs_0.1pval_0.58log2fc.csv", row.names = FALSE)

MC_res <- results(deseq2Data, independentFiltering = TRUE,contrast = c("sampletype","Medium","Control"), alpha=0.1)#build results table
summary(MC_res)
MC_shrunken <- lfcShrink(deseq2Data, coef="sampletype_Medium_vs_Control", type="apeglm", res=MC_res)
summary(MC_shrunken)# will summarize the results using the alpha threshold: FDR < 0.1 (padj/FDR is used even though the output says p-value < 0.05).
table(MC_shrunken$padj < 0.1)
table(MC_shrunken$padj < 0.1 & MC_shrunken$log2FoldChange > 0.584) #upreg
table(MC_shrunken$padj < 0.1 & MC_shrunken$log2FoldChange < -0.584) #downreg
# writing sig genes to file
MC_resSigUP = MC_shrunken[ which(MC_shrunken$padj < 0.1 & MC_shrunken$log2FoldChange > 0.584), ]
MC_resSigDown = MC_shrunken[ which(MC_shrunken$padj < 0.1 & MC_shrunken$log2FoldChange < -0.584), ]
MC_resSig = rbind(MC_resSigUP, MC_resSigDown)
write.csv(LC_resSig,"MC.DEGs_0.1pval_0.58log2fc.csv", row.names = FALSE)

HC_res <- results(deseq2Data, independentFiltering = TRUE,contrast = c("sampletype","High","Control"), alpha=0.1)#build results table
summary(HC_res)
HC_shrunken <- lfcShrink(deseq2Data, coef="sampletype_High_vs_Control", type="apeglm", res=HC_res)
summary(HC_shrunken)# will summarize the results using the alpha threshold: FDR < 0.1 (padj/FDR is used even though the output says p-value < 0.05).
table(HC_shrunken$padj < 0.1)
table(HC_shrunken$padj < 0.1 & HC_shrunken$log2FoldChange > 0.584) #upreg
table(HC_shrunken$padj < 0.1 & HC_shrunken$log2FoldChange < -0.584) #downreg
# writing sig genes to file
HC_resSigUP = HC_shrunken[ which(HC_shrunken$padj < 0.1 & HC_shrunken$log2FoldChange > 0.584), ]
HC_resSigDown = HC_shrunken[ which(HC_shrunken$padj < 0.1 & HC_shrunken$log2FoldChange < -0.584), ]
HC_resSig = rbind(HC_resSigUP, HC_resSigDown)
write.csv(HC_resSig,"HC.DEGs_0.1pval_0.58log2fc.csv", row.names = FALSE)

# Another way to extract sig DEGs
### Set thresholds
padj.cutoff = 0.1
lfc.cutoff = 0.584 # this translates to an actual fold change of 1.5 which is pretty reasonable.

library(tidyverse)
LC_shrunken_tb <- LC_shrunken %>%
  data.frame() %>%
  rownames_to_column(var="Gene") %>% 
  as_tibble() # convert result table into a tibble

MC_shrunken_tb <- MC_shrunken %>%
  data.frame() %>%
  rownames_to_column(var="Gene") %>% 
  as_tibble()
  
HC_shrunken_tb <- HC_shrunken %>%
  data.frame() %>%
  rownames_to_column(var="Gene") %>% 
  as_tibble()

# subset that table to only keep the significant genes using our pre-defined thresholds:
sig_up <- LC_shrunken_tb %>%
        filter(padj < padj.cutoff & log2FoldChange > lfc.cutoff)#8 for LC, 20 for HC, 120 for MC
sig_down <- LC_shrunken_tb %>%
        filter(padj < padj.cutoff & log2FoldChange < -lfc.cutoff)#0 for LC, HC, MC
sig_up
sig_down

## HEATMAP
# 1. LC
LC_res <- LC_res[order(LC_res$padj),]
#Create file for count visualization by vst normalization
rld <- rlog(deseq2Data)

mat <- assay(rld)[head( match(row.names(LC_res), row.names(rld)) , 30), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[, c("sampletype")])
df <- data.frame(sample=df$`colData(rld)[, c("sampletype")]`)
rownames(df) <- colnames(mat)

#Reorder sample column for heatmap
#                   WT                 KO   
#mat <- mat[,c( 5, 6, 7, 8,     1, 2, 3, 4)]

#Produce heatmap
png("LC_TopGenesHeatmap.png")
pheatmap(mat, annotation_col = df, cluster_rows = T , cluster_cols = T ,show_colnames = F)
dev.off()

# 2. MC
MC_res <- MC_res[order(MC_res$padj),]
#Create file for count visualization by vst normalization
rld <- rlog(deseq2Data)

mat <- assay(rld)[head( match(row.names(MC_res), row.names(rld)) , 30), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[, c("sampletype")])
df <- data.frame(sample=df$`colData(rld)[, c("sampletype")]`)
rownames(df) <- colnames(mat)

#Reorder sample column for heatmap
#                   WT                 KO   
#mat <- mat[,c( 5, 6, 7, 8,     1, 2, 3, 4)]

#Produce heatmap
png("MC_TopGenesHeatmap.png")
pheatmap(mat, annotation_col = df, cluster_rows = T , cluster_cols = T ,show_colnames = F)
dev.off()

# 3. HC
HC_res <- HC_res[order(HC_res$padj),]
#Create file for count visualization by vst normalization
rld <- rlog(deseq2Data)

mat <- assay(rld)[head( match(row.names(HC_res), row.names(rld)) , 30), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[, c("sampletype")])
df <- data.frame(sample=df$`colData(rld)[, c("sampletype")]`)
rownames(df) <- colnames(mat)

#Reorder sample column for heatmap
#                   WT                 KO   
#mat <- mat[,c( 5, 6, 7, 8,     1, 2, 3, 4)]

#Produce heatmap
png("HC_TopGenesHeatmap.png")
pheatmap(mat, annotation_col = df, cluster_rows = T , cluster_cols = T ,show_colnames = F)
dev.off()


# filtering threshold that has been used to filter low count genes can be extracted from the results metadata. (if independentFiltering = TRUE mentioned)
# independent filtering : some genes may have low counts and they may have high dispersion which results in less significance of these genes. Such genes can be removed
# prior to FDR procedure which improves the power of the test. Independent filtering removes these weakly expressed genes from the input to the FDR procedure, 
# more significant genes can be found among those that are kept, and this improves the power of the test. This approach is known as independent filtering

# The mean expression threshold used by DESeq2 for independentfiltering is defined automatically by the software.

metadata(LC_res)$filterThreshold # This means that genes whith basemean < 4.215383 have been filtered. This represents 48.46939% of all tested genes!
#NOTE: on p-values set to NA
metadata(MC_res)$filterThreshold
metadata(HC_res)$filterThreshold

#If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p-value and adjusted p-value will all be set to NA.
#If a row contains a sample with an extreme count outlier then the p-value and adjusted p-value will be set to NA. These outlier counts are detected by Cook’s distance.
#If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p-value will be set to NA.
###

# QC FOR DE ANALYSIS USING DESEQ2
# TRANSFORM NORMALIZED COUNTS USING THE RLOG FUNCTION

# 1. PCA
# To improve the distances/clustering for the PCA and heirarchical clustering visualization methods, we need to moderate the variance across the mean by applying 
# the rlog transformation to the normalized counts.
rld <- rlog(dds, blind=TRUE)# The blind=TRUE argument results in a transformation unbiased to sample condition information. When performing quality assessment, it is important to include this option.
pdf("PCA.pdf")
pca_data <- plotPCA(rld, intgroup = "sampletype", returnData = TRUE, ntop=50)
ggplot(pca_data, aes(x = PC1, y = PC2, color = sampletype, label = name)) +
  geom_point() +
  geom_text()+
   stat_ellipse()
dev.off()


fviz_eig(rld, addlabels = T, ylim = c(0, 100))

 
# INVESTIGATING DE RESULTS
# The MA-plot provides a global view of the differential genes, with the log2 fold change on the y-axis over the mean of normalized counts.
# Genes that pass the significance threshold (adjusted p.value 0.05) are colored in red.

png("LC_MAplot_0.1.png")
plotMA(LC_shrunken, alpha = 0.1, 
       main = "Low vs Control: p.adj.value < 0.1", colSig="red", ylim = c(-4,4))
dev.off()

png("MC_MAplot_0.1.png")
plotMA(MC_shrunken, alpha = 0.01, 
       main = "Medium vs Control: p.adj.value < 0.1", colSig="red", ylim = c(-4,4))
dev.off()

png("HC_MAplot_0.1.png")
plotMA(HC_shrunken, alpha = 0.1, 
       main = "High vs Control: p.adj.value < 0.1", colSig="red", ylim = c(-4,4))
dev.off()


# p-value histograms
png("Low_vs_Control_PvalueHistogram_0.1.png")
hist(LC_res$padj, 
     col="grey", border="white", xlab="", ylab="", main="frequencies of adj. p-values\n(all genes)")
dev.off()
png("Medium_vs_Control_PvalueHistogram_0.1.png")
hist(MC_res$padj, 
     col="grey", border="white", xlab="", ylab="", main="frequencies of adj. p-values\n(all genes)")
dev.off()

png("High_vs_Control_PvalueHistogram_0.1.png")
hist(HC_res$padj, 
     col="grey", border="white", xlab="", ylab="", main="frequencies of adj. p-values\n(all genes)")
dev.off()

# HEATMAP OF DE genes.
# This requires normalized gene COUNTS
library(ComplexHeatmap)
ntd <- normTransform(deseq2Data, f = log2, pc = 1)
ntd
#we can extract normalized and log2-transformed expression values via:
assay(ntd)[1:5,]

top.genes <- order(rowVars(assay(ntd)), decreasing = TRUE)[1:50]
top.ntd <- ntd[top.genes,]
assay(top.ntd) <- assay(top.ntd) - rowMeans(assay(top.ntd))

pdf("hhhh.pdf")
Heatmap(assay(top.ntd), 
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        name = "Expression")
dev.off()

pdf("genecounts.pdf")
plotCounts(deseq2Data, gene=which.min(res$padj), intgroup="sampletype")
dev.off()


# 1. Low vs control 
#heatmap of the genes that show differential expression with adjusted p-value 0.1 :
# identify genes with the desired adjusted p-value cut-off
LC.DGEgenes <- rownames(subset(LC_DGEs, padj < 0.1))
rlog.dge <- DESeq.rlog[LC.DGEgenes,] %>% assay # extract rlog-transformed values into a matrix
#rlog.dge <- assay(rlog(LC.DGE.results.sorted, blind=FALSE))
pdf("Low_vs_Control_Heatmap_0.1.pdf")
pheatmap(rld_mat, scale="none", show_rownames = TRUE, main = "DGE (no scaling)")# heatmap of DEG sorted by p.adjust
dev.off()

LC_res<-LC_res[order(LC_res$padj),]
rld <- rlogTransformation(deseq2Data, blind=TRUE)
vd <- varianceStabilizingTransformation(deseq2Data, blind=TRUE)
library("genefilter")
topVarGenes <- head(order(rowVars(assay(vd)),decreasing=TRUE),200)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("ID","sampletype", "sizeFactor")])
pdf("Low_vs_Control_Heatmap_0.1.pdf")
pheatmap(mat, show_rownames=FALSE)
dev.off()



# 2. Medium vs Control
MC.DGEgenes <- rownames(subset(MC.DGE.results.sorted, padj < 0.1))
rlog.dge <- DESeq.rlog[MC.DGEgenes,] %>% assay# extract rlog-transformed values into a matrix
pdf("Medium_vs_Control_Heatmap_0.1.pdf")
pheatmap(rlog.dge, scale="none", show_rownames = FALSE, main = "DGE (no scaling)")# heatmap of DEG sorted by p.adjust
dev.off()

# 3. High vs Control
HC.DGEgenes <- rownames(subset(HC.DGE.results.sorted, padj < 0.1))
rlog.dge <- DESeq.rlog[HC.DGEgenes,] %>% assay# extract rlog-transformed values into a matrix
pdf("High_vs_Control_Heatmap_0.1.pdf")
pheatmap(rlog.dge, scale="none", show_rownames = FALSE, main = "DGE (no scaling)")# heatmap of DEG sorted by p.adjust
dev.off()

## Volcano plot with "significant" genes labeled
# 1. LC

## Merge LC_res with normalized count data
LC_resdata <- merge(as.data.frame(LC_res), as.data.frame(counts(deseq2Data, normalized=TRUE)), by="row.names", sort=FALSE)
names(LC_resdata)[1] <- "Gene"
head(LC_resdata)
## Write results
#write.csv(resdata, file="LC_res_diffexpr-results_default.csv")

volcanoplot <- function (LC_res, lfcthresh=0.58, sigthresh=0.1, main="Low vs Control Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(LC_res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(LC_res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(LC_res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(LC_res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(LC_res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("LC_diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(LC_resdata, lfcthresh=0.58, sigthresh=0.1, textcx=.6, xlim=c(-2.3, 2))
dev.off()

# 2. MC

## Merge MC_res with normalized count data
MC_resdata <- merge(as.data.frame(MC_res), as.data.frame(counts(deseq2Data, normalized=TRUE)), by="row.names", sort=FALSE)
names(MC_resdata)[1] <- "Gene"
head(MC_resdata)
## Write results
#write.csv(resdata, file="LC_res_diffexpr-results_default.csv")

volcanoplot <- function (MC_res, lfcthresh=0.58, sigthresh=0.1, main="Low vs Control Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(MC_res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(MC_res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(MC_res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(MC_res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(MC_res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("MC_diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(MC_resdata, lfcthresh=0.58, sigthresh=0.1, textcx=.6, xlim=c(-2.3, 2))
dev.off()

# 3. HC

## Merge HC_res with normalized count data
HC_resdata <- merge(as.data.frame(HC_res), as.data.frame(counts(deseq2Data, normalized=TRUE)), by="row.names", sort=FALSE)
names(HC_resdata)[1] <- "Gene"
head(HC_resdata)
## Write results
#write.csv(resdata, file="LC_res_diffexpr-results_default.csv")

volcanoplot <- function (HC_res, lfcthresh=0.58, sigthresh=0.1, main="Low vs Control Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(HC_res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(HC_res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(HC_res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(HC_res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(HC_res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("HC_diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(HC_resdata, lfcthresh=0.58, sigthresh=0.1, textcx=.6, xlim=c(-2.3, 2))
dev.off()


##GSEA
# load the GTF file using GenomicFeatures and group the exons based on the gene annotation.
library(GenomicFeatures)
txdb = makeTxDbFromGFF("/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/genomic.gtf", format="gtf")
(ebg = exonsBy(txdb, by="gene"))

#### GSEA
search_kegg_organism('oas', by='kegg_code')
   kegg_code scientific_name common_name
66       oas      Ovis aries       sheep
