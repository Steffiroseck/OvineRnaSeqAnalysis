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
setwd("Z:/SHEEP_TRANSCRIPTOMICS/RNA-Seq/RNA-SEQ-11-12-2023/DESEQ2/Deseq2_featurecounts_Hisat2/")

# prefix names for writing to files
prefix1="Low_vs_Control"
prefix2="Medium_vs_Control"
prefix3="High_vs_Control"
deseqprefix="Deseq2"

# READ COUNT DATA FROM FEATURECOUNTS
# In this analysis, 2 samples were removed. These two samples info need to be removed from the count matrix file and metadata as well. 
countData<-read.csv("LambAllSamples.featureCounts.RmatrixNew.csv", sep=",", header=T,check.names=F) # This file has gene ids and sample files only
countData<-read.csv("LambAllSamples.featureCounts.Rmatrix.txt", sep="\t", header=T,check.names=F)# for star out
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

# remove rows with only 0
deseq2Data <- deseq2Data[rowSums(assay(deseq2Data)) > 0, ]
dim(deseq2Data)# 24950 for hisat and 25043 for star

keep <- rowSums(counts(deseq2Data)) >= 3
deseq2Data <- deseq2Data[keep,] #24053 for hisat and 24119 for star

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
# betaPrior=F is the default. this F option should be used if you want to perfrom LFC cutoff downstram
# The above code does
# 1. sequencing depth normalization between the samples (estimateSizeFactors)
# 2. gene-wise dispersion estimates across all samples (estimateDispersions)
# 3. fits a negative binomial GLM and applies Wald statistics to each gene (nbinomWaldTest)

resultsNames(deseq2Data) # tells you which types of values can be extracted with results()

# send normalized counts to tab delimited file for GSEA, etc.
write.csv(as.data.frame(counts(deseq2Data),normalized=T), "deseq2Data_normalized_counts.txt", sep = '\t', row.names=F)

# extract results for different conditions. 
# The default approach of DESeq uses a null hypothesis of 0 LFC. In this case, each gene is compared across treatment for this null hypothesis, 
# and considered significantly different if smaller than the FDR. When you use a fold change cut off (of lets say 1.5) or LFC 0.584, your null 
# hypothesis is that genes between the experiment differ more than 1.5 fold change, and you will be testing that hypothesis with FDR cut off.

# 1. Low vs Control results

LC_res <- results(deseq2Data, independentFiltering = TRUE,contrast = c("sampletype","Low","Control"), alpha=0.1)#build results table
summary(LC_res)
dim(LC_res)
LC_shrunken <- lfcShrink(deseq2Data, coef="sampletype_Low_vs_Control", type="apeglm", res=LC_res)
summary(LC_shrunken)# will summarize the results using the alpha threshold: FDR < 0.1 (padj/FDR is used even though the output says p-value < 0.05).
dim(LC_shrunken)
# save data results and normalized reads to csv
LC_resdata <- merge(as.data.frame(LC_shrunken), as.data.frame(counts(deseq2Data,normalized =TRUE)), by = 'row.names', sort = FALSE)
# write this to a file as it is a file of normalized counts
write.csv(LC_resdata, file = paste0(deseqprefix, "-results-with-normalized.csv"),row.names=F)

# order results by padj value (most significant to least)
LC_shrunken= subset(LC_shrunken, padj<0.1)#only 30 genes remaining
LC_shrunken <- LC_shrunken[order(LC_shrunken$padj),]# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj
LC_shrunken_resdata <- merge(as.data.frame(LC_shrunken), as.data.frame(counts(deseq2Data,normalized =TRUE)), by = 'row.names', sort = FALSE)
# remove columns with Medium and High info since we need this information for only Low vs Control 
LC_res_df <- LC_shrunken_resdata %>% dplyr::select(-contains(c("Medium", "High")))
names(LC_res_df)[1] <- 'gene'
write.csv(LC_res_df, file = paste0(prefix1, "-results-with-normalized.csv"),row.names=F)

# produce DataFrame of results of statistical tests
mcols(LC_shrunken, use.names = T)
write.csv(as.data.frame(mcols(LC_shrunken, use.name = T)),file = paste0(prefix1, "-test-conditions.csv"))

table(LC_shrunken$padj < 0.1)
table(LC_shrunken$padj < 0.1 & LC_shrunken$log2FoldChange > 0.584) #upreg
table(LC_shrunken$padj < 0.1 & LC_shrunken$log2FoldChange < -0.584) #downreg
# writing sig genes to file
LC_resSigUP = LC_shrunken[ which(LC_shrunken$padj < 0.1 & LC_shrunken$log2FoldChange > 0.584), ]
LC_resSigDown = LC_shrunken[ which(LC_shrunken$padj < 0.1 & LC_shrunken$log2FoldChange < -0.584), ]
LC_resSig = rbind(LC_resSigUP, LC_resSigDown)
write.csv(LC_resSig,file = paste0(prefix1, ".DEGs_0.1pval_0.58log2fc.csv"), row.names = T)

# 2. Medium vs Control results

MC_res <- results(deseq2Data, independentFiltering = TRUE,contrast = c("sampletype","Medium","Control"), alpha=0.1)#build results table
summary(MC_res)
MC_shrunken <- lfcShrink(deseq2Data, coef="sampletype_Medium_vs_Control", type="apeglm", res=MC_res)
summary(MC_shrunken)# will summarize the results using the alpha threshold: FDR < 0.1 (padj/FDR is used even though the output says p-value < 0.05).
dim(MC_shrunken)

# order results by padj value (most significant to least)
MC_shrunken= subset(MC_shrunken, padj<0.1)
MC_shrunken <- MC_shrunken[order(MC_shrunken$padj),]# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj
# save data results and normalized reads to csv
MC_resdata <- merge(as.data.frame(MC_shrunken), as.data.frame(counts(deseq2Data,normalized =TRUE)), by = 'row.names', sort = FALSE)
# remove columns with Medium and High info 
MC_res_df <- MC_resdata %>% dplyr::select(-contains(c("Low", "High")))
names(MC_res_df)[1] <- 'gene'
write.csv(MC_res_df, file = paste0(prefix2, "-results-with-normalized.csv"))

# produce DataFrame of results of statistical tests
mcols(MC_shrunken, use.names = T)
write.csv(as.data.frame(mcols(MC_shrunken, use.name = T)),file = paste0(prefix2, "-test-conditions.csv"))

table(MC_shrunken$padj < 0.1)
table(MC_shrunken$padj < 0.1 & MC_shrunken$log2FoldChange > 0.584) #upreg
table(MC_shrunken$padj < 0.1 & MC_shrunken$log2FoldChange < -0.584) #downreg
# writing sig genes to file
MC_resSigUP = MC_shrunken[ which(MC_shrunken$padj < 0.1 & MC_shrunken$log2FoldChange > 0.584), ]
MC_resSigDown = MC_shrunken[ which(MC_shrunken$padj < 0.1 & MC_shrunken$log2FoldChange < -0.584), ]
MC_resSig = rbind(MC_resSigUP, MC_resSigDown)
write.csv(MC_resSig,file = paste0(prefix2, ".DEGs_0.1pval_0.58log2fc.csv"), row.names = T)

# 3. High vs Control results

HC_res <- results(deseq2Data, independentFiltering = TRUE,contrast = c("sampletype","High","Control"), alpha=0.1)#build results table
summary(HC_res)
HC_shrunken <- lfcShrink(deseq2Data, coef="sampletype_High_vs_Control", type="apeglm", res=HC_res)
summary(HC_shrunken)# will summarize the results using the alpha threshold: FDR < 0.1 (padj/FDR is used even though the output says p-value < 0.05).

# order results by padj value (most significant to least)
HC_shrunken= subset(HC_shrunken, padj<0.1)
HC_shrunken <- HC_shrunken[order(HC_shrunken$padj),]# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj
# save data results and normalized reads to csv
HC_resdata <- merge(as.data.frame(HC_shrunken), as.data.frame(counts(deseq2Data,normalized =TRUE)), by = 'row.names', sort = FALSE)
# remove columns with Medium and Low info 
HC_res_df <- HC_resdata %>% dplyr::select(-contains(c("Low", "Medium")))
names(HC_res_df)[1] <- 'gene'
write.csv(HC_res_df, file = paste0(prefix3, "-results-with-normalized.csv"), row.names=F)

# produce DataFrame of results of statistical tests
mcols(HC_shrunken, use.names = T)
write.csv(as.data.frame(mcols(HC_shrunken, use.name = T)),file = paste0(prefix3, "-test-conditions.csv"))

table(HC_shrunken$padj < 0.1)
table(HC_shrunken$padj < 0.1 & HC_shrunken$log2FoldChange > 0.584) #upreg
table(HC_shrunken$padj < 0.1 & HC_shrunken$log2FoldChange < -0.584) #downreg
# writing sig genes to file
HC_resSigUP = HC_shrunken[ which(HC_shrunken$padj < 0.1 & HC_shrunken$log2FoldChange > 0.584), ]
HC_resSigDown = HC_shrunken[ which(HC_shrunken$padj < 0.1 & HC_shrunken$log2FoldChange < -0.584), ]
HC_resSig = rbind(HC_resSigUP, HC_resSigDown)
write.csv(HC_resSig,file = paste0(prefix3, ".DEGs_0.1pval_0.58log2fc.csv"), row.names = T)

####################################################################################
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

pdf(paste0(prefix1, "-MAplot_initial_analysis.pdf"))
plotMA(LC_res, ylim=c(-8,8),main = "RNAseq experiment", alpha=0.1,colNonSig = "black", colSig = "red")
dev.off()

pdf(paste0(prefix2, "-MAplot_initial_analysis.pdf"))
plotMA(MC_res, ylim=c(-8,8),main = "RNAseq experiment", alpha=0.1,colNonSig = "black", colSig = "red")
dev.off()

pdf(paste0(prefix3, "-MAplot_initial_analysis.pdf"))
plotMA(HC_res, ylim=c(-8,8),main = "RNAseq experiment", alpha=0.1,colNonSig = "black", colSig = "red")
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(deseq2Data, blind=T)
vsd <- varianceStabilizingTransformation(deseq2Data, blind=T)

# save normalized values
write.csv(assay(rld),file = paste0(deseqprefix, "-rlog-transformed-counts.txt"), sep = '\t', row.names=F)
write.csv(assay(vsd),file = paste0(deseqprefix, "-vst-transformed-counts.txt"), sep = '\t',  row.names=F)

# plot to show effect of transformation
# axis is square root of variance over the mean for all samples
pdf(paste0(deseqprefix, "-variance_stabilizing.pdf"))
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(deseq2Data)[,1] / sizeFactors(deseq2Data)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.off()

# heatmap of clustering analysis
# excerpts from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/

pdf(paste0(deseqprefix, "-clustering.pdf"))
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(deseq2Data), paste(sampletype, ID, sep=" : "))
#Or if you want conditions use:
#rownames(mat) <- colnames(mat) <- with(colData(deseq2Data),condition)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13))
dev.off()

#Principal components plot shows additional but rough clustering of samples
rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))

# set condition
treatments = c("Control","Control","Control","Control","Control","Control","Low", "Low","Low","Low","Low","Medium", "Medium","Medium","Medium","Medium","High","High","High","High","High","High")
condition <- treatments
scores <- data.frame(pc$x, condition)

(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 5)
  + ggtitle("Principal Components")
  + scale_colour_brewer(name = " ", palette = "Set1")
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.position = c(.9,.2),
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA)
  ))

ggsave(pcaplot,file=paste0(deseqprefix, "PCA-ggplot2.pdf"))

# heatmap of data

# 1000 top expressed genes with heatmap.2
pdf(paste0(deseqprefix, "-HEATMAP.pdf"))
select <- order(rowMeans(counts(deseq2Data,normalized=T)),decreasing=T)[1:1000]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
heatmap.2(assay(vsd)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="1000 Top Expressed Genes Heatmap")
dev.off()

## Volcano plot with "significant" genes labeled
# 1. LC

## Merge LC_res with normalized count data
LC_resdata <- merge(as.data.frame(LC_res), as.data.frame(counts(deseq2Data, normalized=TRUE)), by="row.names", sort=FALSE)
names(LC_resdata)[1] <- "Gene"
head(LC_resdata)
dim(LC_resdata)

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
pdf(paste0(prefix1, "-Volcano.pdf"),pointsize=10, width=12, height=11)
volcanoplot(LC_resdata, lfcthresh=0.58, sigthresh=0.1, textcx=.6, xlim=c(-2.3, 2))
dev.off()

# 2. MC

## Merge MC_res with normalized count data
MC_resdata <- merge(as.data.frame(MC_res), as.data.frame(counts(deseq2Data, normalized=TRUE)), by="row.names", sort=FALSE)
names(MC_resdata)[1] <- "Gene"
head(MC_resdata)

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
pdf(paste0(prefix2, "-Volcano.pdf"), pointsize=10, width=12, height=11)
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
pdf(paste0(prefix3, "-Volcano.pdf"),pointsize=10, width=12, height=11)
volcanoplot(HC_resdata, lfcthresh=0.58, sigthresh=0.1, textcx=.6, xlim=c(-2.3, 2))
dev.off()

# Plot gene counts
pdf("genecounts.pdf")
plotCounts(deseq2Data, gene=which.min(res$padj), intgroup="sampletype")
dev.off()

pdf("SCD_genecounts.pdf")
plotCounts(deseq2Data, gene="SCD", intgroup="sampletype") 
dev.off()

pdf("CCR3_genecounts.pdf")
plotCounts(deseq2Data, gene="CCR3", intgroup="sampletype") 
dev.off()

# p-value histograms
pdf(paste0(prefix1,"_PvalueHistogram_0.1.pdf"))
hist(LC_res$padj, 
     col="grey", border="white", xlab="", ylab="", main="frequencies of adj. p-values\n(all genes)")
dev.off()

pdf(paste0(prefix2,"_PvalueHistogram_0.1.pdf"))
hist(MC_res$padj, 
     col="grey", border="white", xlab="", ylab="", main="frequencies of adj. p-values\n(all genes)")
dev.off()

pdf(paste0(prefix3,"_PvalueHistogram_0.1.pdf"))
hist(HC_res$padj, 
     col="grey", border="white", xlab="", ylab="", main="frequencies of adj. p-values\n(all genes)")
dev.off()

## Another code for heatmap (Low vs Control)
LC_resSig <- LC_resSig[order(LC_resSig$padj),]
#Create file for count visualization by vst normalization
mat <- assay(vsd)[head( match(row.names(LC_res), row.names(vsd))), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[, c("sampletype")])
df <- data.frame(sample=df$`colData(vsd)[, c("sampletype")]`)
rownames(df) <- colnames(mat)

#Reorder sample column for heatmap
#                   WT                 KO   
#mat <- mat[,c( 5, 6, 7, 8,     1, 2, 3, 4)]

#Produce heatmap
pdf("LC_SigGenesHeatmap.pdf")
pheatmap(mat, annotation_col = df, cluster_rows = T , cluster_cols = T ,show_colnames = F)
dev.off()
