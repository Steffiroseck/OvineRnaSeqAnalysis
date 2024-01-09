library(DESeq2)
library(dplyr)
library(pheatmap)
library(tidyverse)
library(ggplot2)
library(calibrate)
library(cluster)
library(RColorBrewer)
library(gplots)
library(genefilter)
library(grDevices)

# Set the working directory
#setwd("Z:/SHEEP_TRANSCRIPTOMICS/RNA-Seq/RNA-SEQ-11-12-2023/DESEQ2/Deseq2_featurecounts_Hisat2/")

# prefix names for writing to files
prefix1="Low_vs_Control"
prefix2="Medium_vs_Control"
prefix3="High_vs_Control"
deseqprefix="Deseq2"

# READ COUNT DATA FROM FEATURECOUNTS
# In this analysis, 2 samples were removed. These two samples info need to be removed from the count matrix file and metadata as well. 
countData<-read.csv("LMH_Control_Lamb_rawCounts.csv",sep=",", header=T, check.names=F)
head(countData)
#countData$Chr<-NULL # run this only if chr column is present
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

# READ THE METADATA TABLE AND TRANSFORM IT IN A WAY DESEQ2 ACCEPTS

# Convert sample variable mappings to an appropriate form that DESeq2 can read
metaData <-read.csv("LMH_C_metadata.csv",sep=",",header=T)
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

# remove rows with only 0
deseq2Data <- deseq2Data[rowSums(assay(deseq2Data)) > 0, ]
dim(deseq2Data)#24405 remaining for LC

keep <- rowSums(counts(deseq2Data)) >= 10
deseq2Data <- deseq2Data[keep,]
dim(deseq2Data)#22721 remaining for LMH_C HISAT, 22733 for STAR

# set control condition as reference
deseq2Data$sampletype <- relevel(deseq2Data$sampletype, ref = "Control")

# It is also advisable to investigate any systematic bias in the sequencing data, such as whether one sample has been sequenced more deeply than others. 
# One can extract size factors using the sizeFactors() function. Usually these size factors should vary around 1, indicating comparable sequencing depth.
deseq2Data<-estimateSizeFactors(deseq2Data)

# compare sizefactors to seq depth
pdf(paste0(deseqprefix,"sizeFactors_deseqdata.pdf"), width=15)
tibble(sample = colnames(deseq2Data), SF = sizeFactors(deseq2Data)) %>%
  ggplot(aes(x = sample, y = SF)) +
  geom_col()
dev.off()

pdf(paste0(deseqprefix,"seqdepth_deseqdata.pdf"), width=15)
tibble(sample = colnames(deseq2Data), SeqDepth = colSums(assay(deseq2Data))) %>%
  ggplot(aes(x = sample, y = SeqDepth)) +
  geom_col()
dev.off()

# DIFFERENTIAL EXPRESSION ANALYSIS WITH DESEQ2
# Null hypothesis: no diff expn between groups, ie log2FC=0. A statistical test, the Wald test, will determine whether the data provides sufficient evidence to conclude that this value is really different from zero.
# Pvalue generated from wald test result in many false positives, hence needs to be adjusted which DESEq2 adjust using BH method and generate adj pvalues.

deseq2Data <- DESeq(deseq2Data, betaPrior = FALSE)

#The betaPrior = T option tells DESeq2 to squeeze the log Fold Changes of lowly expressed genes toward 0.
# betaPrior=F is the default. this F option should be used if you want to perfrom LFC cutoff downstream
# The above code does
# 1. sequencing depth normalization between the samples (estimateSizeFactors)
# 2. gene-wise dispersion estimates across all samples (estimateDispersions)
# 3. fits a negative binomial GLM and applies Wald statistics to each gene (nbinomWaldTest)

resultsNames(deseq2Data) # tells you which types of values can be extracted with results()

# send normalized counts to tab delimited file for GSEA, etc.
write.csv(as.data.frame(counts(deseq2Data),normalized=T), "LMH_vs_Control_deseq2Data_normalized_counts.csv", sep = '\t', row.names=F)

# extract results for different conditions. 
# get gene expression table
# at this step independent filtering is applied by default to remove low count genes
# independent filtering can be turned off by passing independentFiltering=FALSE to results

# The default approach of DESeq uses a null hypothesis of 0 LFC. In this case, each gene is compared across treatment for this null hypothesis, 
# and considered significantly different if smaller than the FDR. When you use a fold change cut off (of lets say 1.5) or LFC 0.584, your null 
# hypothesis is that genes between the experiment differ more than 1.5 fold change, and you will be testing that hypothesis with FDR cut off.

# 1. Low vs Control results
LC_res <- results(deseq2Data, contrast=c("sampletype","Low", "Control"))#build results table
summary(LC_res)
MC_res <- results(deseq2Data, contrast=c("sampletype","Medium", "Control"))#build results table
summary(MC_res)
HC_res <- results(deseq2Data, contrast=c("sampletype","High", "Control"))#build results table
summary(HC_res)
# alpha=0.1 and lfc =0 : 10 up and 34 DOWN for LC

# Get summary of differential gene expression with adjusted p value cut-off at 0.05
summary(results(deseq2Data, alpha=0.05))

#Order gene expression table by adjusted p value (Benjamini-Hochberg FDR method)
LC_resOrdered <- LC_res[order(LC_res$padj),]
LC_resOrdered
LC_res_df <- as.data.frame(LC_resOrdered)

MC_resOrdered <- MC_res[order(MC_res$padj),]
MC_resOrdered
MC_res_df <- as.data.frame(MC_resOrdered)

HC_resOrdered <- HC_res[order(HC_res$padj),]
HC_resOrdered
HC_res_df <- as.data.frame(HC_resOrdered)

#Note: You may get some genes with p value set to NA. 
#This is due to all samples have zero counts for a gene or there is extreme outlier count for a gene or that gene is subjected to independent filtering by DESeq2.

sum(LC_res$padj < 0.1, na.rm=TRUE)# 31 
sum(LC_res$padj < 0.05, na.rm=TRUE)#9
sum(LC_res$padj < 0.05 & LC_res$log2FoldChange < -0.58, na.rm=TRUE) # 9
LCresSig <- subset(LC_resOrdered, padj < 0.1 & log2FoldChange > 0 | padj < 0.1 & log2FoldChange <0)# These are the significant genes with p <0.1 and lfc < or > 0
dim(LCresSig)#31
write.csv(LCresSig,file = paste0(prefix1, ".DEGs_0.1pval_log2fc_zero.csv"), row.names = T)#44 sig genes written to a file
LCresSigUP <- subset(LCresSig, log2FoldChange >0.584)#3
LCresSigDOWN <- subset(LCresSig, log2FoldChange < -0.584)#25
dim(LCresSigUP)
dim(LCresSigDOWN)
LCresSigUP_DOWN<-rbind(LCresSigUP,LCresSigDOWN) 
dim(LCresSigUP_DOWN)
write.csv(LCresSigUP_DOWN,file = paste0(prefix1, ".DEGs_0.1pval_log2fc_0.58.csv"), row.names = T)# 26 sig genes written to a file

sum(MC_res$padj < 0.1, na.rm=TRUE)# 850
sum(MC_res$padj < 0.05, na.rm=TRUE)#333
sum(MC_res$padj < 0.05 & MC_res$log2FoldChange < -0.58, na.rm=TRUE) # 7
MCresSig <- subset(MC_resOrdered, padj < 0.1 & log2FoldChange > 0 | padj < 0.1 & log2FoldChange <0)# These are the significant genes with p <0.1 and lfc < or > 0
dim(MCresSig)#850
write.csv(MCresSig,file = paste0(prefix2, ".DEGs_0.1pval_log2fc_zero.csv"), row.names = T)#44 sig genes written to a file
MCresSigUP <- subset(MCresSig, log2FoldChange >0.584)#48
MCresSigDOWN <- subset(MCresSig, log2FoldChange < -0.584)#271
dim(MCresSigUP)
dim(MCresSigDOWN)
MCresSigUP_DOWN<-rbind(MCresSigUP,MCresSigDOWN) 
dim(MCresSigUP_DOWN)#319
write.csv(MCresSigUP_DOWN,file = paste0(prefix2, ".DEGs_0.1pval_log2fc_0.58.csv"), row.names = T)# 26 sig genes written to a file


sum(HC_res$padj < 0.1, na.rm=TRUE)# 44 for HC
sum(HC_res$padj < 0.05, na.rm=TRUE)#28 for HC
sum(HC_res$padj < 0.05 & HC_res$log2FoldChange < -0.58, na.rm=TRUE) # 20 for HC
HCresSig <- subset(HC_resOrdered, padj < 0.1 & log2FoldChange > 0 | padj < 0.1 & log2FoldChange <0)# These are the significant genes with p <0.1 and lfc < or > 0
dim(HCresSig)#44
write.csv(HCresSig,file = paste0(prefix3, ".DEGs_0.1pval_log2fc_zero.csv"), row.names = T)#44 sig genes written to a file
HCresSigUP <- subset(HCresSig, log2FoldChange >0.584)#1
HCresSigDOWN <- subset(HCresSig, log2FoldChange < -0.584)#25
dim(HCresSigUP)
dim(HCresSigDOWN)
HCresSigUP_DOWN<-rbind(HCresSigUP,HCresSigDOWN) 
dim(HCresSigUP_DOWN)
write.csv(HCresSigUP_DOWN,file = paste0(prefix3, ".DEGs_0.1pval_log2fc_0.58.csv"), row.names = T)# 26 sig genes written to a file


# Shrinkage estimations
# The shrinkage of effect size (LFC) helps to remove the low count genes (by shrinking towards zero). 
# The low or highly variable read count genes can give large estimates of LFCs which may not represent 
# true difference in changes in gene expression between two conditions.
HC_resLFC <- lfcShrink(deseq2Data, coef="sampletype_Low_vs_Control", type="apeglm")
head(HC_resLFC)
dim(HC_resLFC)#21785

#To find genes which change in either condition, you would use a likelihood ratio test. You just need to use ~1 as the reduced formula. so
dds <- DESeq(deseq2Data, test="LRT", reduced=~1)
res <- results(dds)
# plot heatmap of genes that shows differences in conditions
gene_vsd <- vst(dds)
pdf(paste0(deseqprefix, "-genes-HEATMAP.pdf"))
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:20]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=50)
heatmap.2(assay(gene_vsd)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=T,
          main="Top 20 Genes Heatmap")
dev.off()


#Then define the set of genes passing an FDR cutoff:
select <- which(res$padj < 0.05)

###################################################################################
# Exploratory data analysis of RNAseq data with DESeq2
#
# these next R scripts are for a variety of visualization, QC and other plots to
# get a sense of what the RNAseq data looks like based on DESEq2 analysis
#
# 1) MA plot
# 2) rlog stabilization and variance stabiliazation
# 3) variance stabilization plot
# 4) heatmap of clustering analysis
# 5) PCA plot
#
####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red

pdf(paste0(prefix1, "-MAplot.pdf"))
plotMA(LC_res, ylim=c(-8,8),main = "RNAseq experiment", alpha=0.1,colNonSig = "grey", colSig = "red")
dev.off()

pdf(paste0(prefix2, "-MAplot.pdf"))
plotMA(MC_res, ylim=c(-8,8),main = "RNAseq experiment", alpha=0.1,colNonSig = "grey", colSig = "red")
dev.off()

pdf(paste0(prefix3, "-MAplot.pdf"))
plotMA(HC_res, ylim=c(-8,8),main = "RNAseq experiment", alpha=0.1,colNonSig = "grey", colSig = "red")
dev.off()

## Plot dispersion estimates
pdf(paste0(deseqprefix, "-DispersionEstimates.pdf"))
plotDispEsts(deseq2Data, main ="Dispersion plot")
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(deseq2Data, blind=FALSE)
head(assay(rld))
hist(assay(rld))
vsd <- varianceStabilizingTransformation(deseq2Data, blind=FALSE)
head(assay(vsd))
hist(assay(vsd))

pdf(paste0(deseqprefix, "-PCA.pdf"))
PCA<-plotPCA(vsd,intgroup='sampletype')
PCA + geom_text(aes(label=name),size=2, position = position_dodge(width=0.9))+ggtitle('PCA plot')
dev.off()

# plot heatmap of genes that shows differences in conditions (top 50)
top50 <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat <- as.matrix(assay(vsd)[top50,])
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd))
#anno <- as.data.frame(colData(gene_vsd)[, c("genotype","condition")])
pdf(paste0(deseqprefix, "Top 50-genes-HEATMAP.pdf"))
pheatmap(mat, annotation_col = anno,  fontsize_row =8)
dev.off()

# save normalized values
write.csv(assay(rld),file = paste0(deseqprefix, "-rlog-transformed-counts.txt"), sep = '\t', row.names=F)
write.csv(assay(vsd),file = paste0(deseqprefix, "-vst-transformed-counts.txt"), sep = '\t',  row.names=F)

# heatmap of clustering analysis
# excerpts from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/

pdf(paste0(deseqprefix, "-clustering using heatmap.pdf"))
distsRL <- dist(t(assay(vsd)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(deseq2Data), paste(sampletype, ID, sep=" : "))
#Or if you want conditions use:
#rownames(mat) <- colnames(mat) <- with(colData(deseq2Data),condition)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13))
dev.off()

## Enhanced Volcano Plot
library("EnhancedVolcano")

pdf(paste0(prefix1, "-EnhancedVolcano.pdf"), width=10, height=9)
EnhancedVolcano(LC_res,
	lab=rownames(LC_res),
	x='log2FoldChange',
	y='padj',
	title='Low vs Control',
	pCutoff=0.1,
	FCcutoff=0.584,
	pointSize=3.0,
	labSize=2,
	legendPosition = 'right',
    legendLabSize = 12,
    drawConnectors = TRUE,
    widthConnectors = 0.3)
dev.off()

pdf(paste0(prefix2, "-EnhancedVolcano.pdf"), width=10, height=9)
EnhancedVolcano(MC_res,
	lab=rownames(MC_res),
	x='log2FoldChange',
	y='padj',
	title='Medium vs Control',
	pCutoff=0.1,
	FCcutoff=0.584,
	pointSize=3.0,
	labSize=2,
	legendPosition = 'right',
    legendLabSize = 12,
    drawConnectors = TRUE,
    widthConnectors = 0.3)
dev.off()

pdf(paste0(prefix3, "-EnhancedVolcano.pdf"), width=10, height=9)
EnhancedVolcano(HC_res,
	lab=rownames(HC_res),
	x='log2FoldChange',
	y='padj',
	title='High vs Control',
	pCutoff=0.1,
	FCcutoff=0.584,
	pointSize=3.0,
	labSize=2,
	legendPosition = 'right',
    legendLabSize = 12,
    drawConnectors = TRUE,
    widthConnectors = 0.3)
dev.off()

## Plot gene counts for 3 common genes in all comparisons

deseq2Data$time <- sample(1:22, ncol(deseq2Data), replace=TRUE)
gene_set <- c("SCD", "CCR3", "SMOC1")
names(gene_set) <- gene_set
df <- lapply(gene_set, \(x) {
    y <- plotCounts(deseq2Data, x, c("sampletype", "time"), returnData=TRUE)
    y$feature <- x
    return(y)
})

df <- do.call(rbind, df)
pdf(paste0(deseqprefix, "-Genecounts_of_common_3_genes.pdf"))
ggplot(df, aes(x=sampletype, y=count, color=sampletype)) +
  geom_jitter(
    position=position_jitterdodge(dodge.width=0.75),
    size=0.75) +
  geom_boxplot(
    position=position_dodge(width=0.5),
    fill=NA, width=0.25, outlier.shape=NA) +
  facet_grid(feature~.)
  dev.off()
  
pdf(paste0(deseqprefix, "ALTERNATE-Genecounts_of_common_3_genes.pdf"))
ggplot(df, aes(x=sampletype, y=count, color=sampletype)) +
  geom_jitter(
    position=position_jitterdodge(dodge.width=0.75),
    size=0.75) +
  geom_boxplot(
    fill=NA, width=0.25, outlier.shape=NA) +
  facet_grid(feature~.)
dev.off()
#################################################################
#################################################################
#
#
#          edgeR
#
#
#################################################################
library(edgeR)
# We will now create a DGEList object to hold our read counts. This object is a container for the counts themsleves, and
#also for all the associated metadata - these include, for example, sample names, gene names and normalisation factors
#once these are computed. The DGEList is an example of the custom task-specific structures that are frequently used in
#Bioconductor to make analyses easier.

dgList <- DGEList(counts=countData, genes=rownames(countData))
dgList$samples
head(dgList$counts) #Many rows!
head(dgList$genes)

# There are approximately 26000 genes in this dataset. However, many of them will not be expressed, or will not be
#represented by enough reads to contribute to the analysis. Removing these genes means that we have ultimately have
#fewer tests to perform, thereby reducing the problems associated with multiple testing. Here, we retain only those genes
#that are represented at least 1cpm reads in at least two samples (cpm=counts per million).

countsPerMillion <- cpm(dgList)
summary(countsPerMillion)
#'summary' is a useful function for exploring numeric data; eg. summary(1:100)
countCheck <- countsPerMillion > 1
head(countCheck)
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]
summary(cpm(dgList)) #compare this to the original summary

#  it is important to normalise RNA-seq both within and between samples. edgeR implements the
# trimmed mean of M-values (TMM) method.
?calcNormFactors
dgList <- calcNormFactors(dgList, method="TMM")

#We can examine inter-sample relationships by producing a plot based on mutlidimensional scaling.
pdf(paste0(prefix3,"-edgeR-MDS.pdf"))
plotMDS(dgList)
dev.off()

# We are now ready to set up the model! We first need to specify our design matrix, which describes the setup of the experiment.
sampleType<- rep("Control", ncol(dgList)) #N=normal; T=tumour
sampleType[grep("Low", colnames(dgList))] <- "Low"
#'grep' is a string matching function.
designMat <- model.matrix(~sampleType)
designMat

# Estimate dispersion
dgList <- estimateGLMCommonDisp(dgList, design=designMat)
dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)

#We can plot the estimates and see how they differ. The biological coefficient of variation (BCV) is the square root of
#the dispersion parameter in the negative binomial model.
pdf(paste0(prefix3,"-edgeR-BCV.pdf"))
plotBCV(dgList)
dev.off()

# We can now find our differentially expressed genes. After fitting the model, we can use the topTags() function to explore
#the results, and set theresholds to identify subsets of differentially expressed genes.

fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit, coef=4)
edgeR_result <- topTags(lrt)
result_edgeR <- as.data.frame(topTags(lrt, n=nrow(countData)))
table(result_edgeR$FDR < 0.05)
table(result_edgeR$FDR < 0.05 & result_edgeR$logFC > 0.584) #2
table(result_edgeR$FDR < 0.05 & result_edgeR$logFC < -0.584) # 21

subset(result_edgeR,FDR < 0.05 & logFC > 0.584)

library("EnhancedVolcano")
pdf(paste0(prefix3, "-edgeR-EnhancedVolcano.pdf"), width=10, height=9)
EnhancedVolcano(result_edgeR,
	lab=rownames(result_edgeR),
	x='logFC',
	y='FDR',
	title='Low vs Control',
	pCutoff=0.1,
	FCcutoff=0.584,
	pointSize=3.0,
	labSize=2,
	legendPosition = 'right',
    legendLabSize = 12,
    drawConnectors = TRUE,
    widthConnectors = 0.3)
dev.off()
