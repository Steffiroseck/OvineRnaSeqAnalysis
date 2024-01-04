##GSEA

# 1- obtaining stats for ranking genes in your experiment,
# 2- creating a named vector out of the DESeq2 result
# 3- Obtaining a gene set from mysigbd
# 4- doing analysis


gffpath="/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0"
gfffile="/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/genomic.gff"

library(clusterProfiler)
library(pathview)
library(stringr)
library(AnnotationHub)
library(ggridges)
library(enrichplot)

# extract annotation DB of sheep
ah <- AnnotationHub()
AnnotationHub::query(ah, c("Ovis", "aries"))
Oaries <- ah[["AH111978"]]

# 1. Low vs Control
# we want the log2 fold change 
prefix1="Low_vs_Control"

# create geneList, which is a list of all genes from deseq2 with p<0.1
entrez<-read.csv("/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/entrezids_new",sep="\t")#35105
LC_res_GO <- LC_res
LC_res_GO <- LC_res_GO[!is.na(LC_res_GO$padj),]
LC_res_GO<-LC_res_GO[order(LC_res_GO$padj),]
#convert rownames to a column  and name it as GeneName so as to merge it with entrez
LC_res_GO$GeneName <- rownames(LC_res_GO)
merge_res_with_entrez<-merge(as.data.frame(LC_res_GO), entrez, by="GeneName")
# again change GeneName to Locus so as to merge with geneTable. Then remove all occurences of "LOC" from GeneName column
colnames(merge_res_with_entrez)[colnames(merge_res_with_entrez) == 'GeneID'] <- 'entrezgene'
merge_res_with_entrez$GeneName <- gsub("LOC", "", merge_res_with_entrez$GeneName)
rownames(merge_res_with_entrez) <- merge_res_with_entrez$entrezgene
original_gene_list <- merge_res_with_entrez$log2FoldChange

# name the vector
names(original_gene_list) <- merge_res_with_entrez$entrezgene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# gene set enrichment
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = Oaries, 
             pAdjustMethod = "BH")
write.csv(gse, file = paste0(prefix1, "-Gene-set-enrichment.csv"),row.names=F)

# when no FDR correction perfrmed, some enriched terms are found. Even after specifying p of 0.1 no enriched terms

pdf(paste0(prefix1,"-Dotplot.pdf"), width=10, height=12)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

# enrichment map
pdf(paste0(prefix1,"-Emapplot.pdf"), width=10, height=12)
x2 <- pairwise_termsim(gse)
emapplot(x2, showCategory = 10)
dev.off()

#category net plot
pdf(paste0(prefix1,"-Cnetplot.pdf"), width=10, height=12)
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
dev.off()

# Ridgeplot
pdf(paste0(prefix1,"-Ridgeplot.pdf"), width=10, height=12)
ridgeplot(gse) + labs(x = "enrichment distribution")
dev.off()


# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
pdf(paste0(prefix1,"-GSEAplot.pdf"))
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
dev.off()

# KEGG ENRICHMENT
# Create a vector of the gene unuiverse
kegg_gene_list <- merge_res_with_entrez$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- merge_res_with_entrez$entrezgene

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "oas"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")
write.csv(kk2, file = paste0(prefix1, "-Kegg-enrichment.csv"),row.names=F)

pdf(paste0(prefix1,"Kegg-Dotplot.pdf"), width=10, height=12)
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
dev.off()

pdf(paste0(prefix1,"Kegg-Emapplot.pdf"), width=10, height=12)
x2 <- pairwise_termsim(kk2)
emapplot(x2, showCategory = 10)
dev.off()

# categorySize can be either 'pvalue' or 'geneNum'
pdf(paste0(prefix1,"Kegg-Cnetplot.pdf"), width=10, height=12)
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
dev.off()

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
pdf(paste0(prefix1,"Kegg-Gseaplot.pdf"), width=10, height=12)
gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)
dev.off()


# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="oas04610", species = kegg_organism)

# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="oas04610", species = kegg_organism, kegg.native = F)
knitr::include_graphics("oas04610.pathview.png")


# 2. Medium vs Control
# we want the log2 fold change 
prefix2="Medium_vs_Control"

# create geneList, which is a list of all genes from deseq2 with p<0.1
entrez<-read.csv("/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/entrezids_new",sep="\t")#35105
MC_res_GO <- MC_resSig#24119
MC_res_GO <- MC_res_GO[!is.na(MC_res_GO$padj),]
MC_res_GO<-MC_res_GO[order(MC_res_GO$padj),]
#convert rownames to a column  and name it as GeneName so as to merge it with entrez
MC_res_GO$GeneName <- rownames(MC_res_GO)
merge_res_with_entrez<-merge(as.data.frame(MC_res_GO), entrez, by="GeneName")
# again change GeneName to Locus so as to merge with geneTable. Then remove all occurences of "LOC" from GeneName column
colnames(merge_res_with_entrez)[colnames(merge_res_with_entrez) == 'GeneID'] <- 'entrezgene'
merge_res_with_entrez$GeneName <- gsub("LOC", "", merge_res_with_entrez$GeneName)

rownames(merge_res_with_entrez) <- merge_res_with_entrez$entrezgene
original_gene_list <- merge_res_with_entrez$log2FoldChange

# name the vector
names(original_gene_list) <- merge_res_with_entrez$entrezgene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# gene set enrichment
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = Oaries, 
             pAdjustMethod = "BH")
write.csv(gse, file = paste0(prefix2, "-Gene-set-enrichment.csv"),row.names=F)

# when no FDR correction perfrmed, some enriched terms are found. Even after specifying p of 0.1 no enriched terms

pdf(paste0(prefix2,"-Dotplot.pdf"), width=10, height=12)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

# enrichment map
pdf(paste0(prefix2,"-Emapplot.pdf"), width=10, height=12)
x2 <- pairwise_termsim(gse)
emapplot(x2, showCategory = 10)
dev.off()

#category net plot
pdf(paste0(prefix2,"-Cnetplot.pdf"), width=10, height=12)
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
dev.off()

# Ridgeplot
pdf(paste0(prefix2,"-Ridgeplot.pdf"), width=10, height=12)
ridgeplot(gse) + labs(x = "enrichment distribution")
dev.off()


# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
pdf(paste0(prefix2,"-GSEAplot.pdf"))
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
dev.off()

# KEGG ENRICHMENT
# Create a vector of the gene unuiverse
kegg_gene_list <- merge_res_with_entrez$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- merge_res_with_entrez$entrezgene

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "oas"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")
write.csv(kk2, file = paste0(prefix2, "-Kegg-enrichment.csv"),row.names=F)

pdf(paste0(prefix2,"Kegg-Dotplot.pdf"), width=10, height=12)
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
dev.off()

pdf(paste0(prefix2,"Kegg-Emapplot.pdf"), width=10, height=12)
x2 <- pairwise_termsim(kk2)
emapplot(x2, showCategory = 10)
dev.off()

# categorySize can be either 'pvalue' or 'geneNum'
pdf(paste0(prefix2,"Kegg-Cnetplot.pdf"), width=10, height=12)
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
dev.off()

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
pdf(paste0(prefix2,"Kegg-Gseaplot.pdf"), width=10, height=12)
gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)
dev.off()


# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="oas04060", species = kegg_organism)

# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="oas04060", species = kegg_organism, kegg.native = F)
knitr::include_graphics("oas04060.pathview.png")

# 3. High vs Control

# we want the log2 fold change 
prefix3="High_vs_Control"

# create geneList, which is a list of all genes from deseq2 with p<0.1
entrez<-read.csv("/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/entrezids_new",sep="\t")#35105
HC_res_GO <- HC_res#24119
HC_res_GO <- HC_res_GO[!is.na(HC_res_GO$padj),]
HC_res_GO<-HC_res_GO[order(HC_res_GO$padj),]
#convert rownames to a column  and name it as GeneName so as to merge it with entrez
HC_res_GO$GeneName <- rownames(HC_res_GO)
merge_res_with_entrez<-merge(as.data.frame(HC_res_GO), entrez, by="GeneName")
# again change GeneName to Locus so as to merge with geneTable. Then remove all occurences of "LOC" from GeneName column
colnames(merge_res_with_entrez)[colnames(merge_res_with_entrez) == 'GeneID'] <- 'entrezgene'
merge_res_with_entrez$GeneName <- gsub("LOC", "", merge_res_with_entrez$GeneName)

rownames(merge_res_with_entrez) <- merge_res_with_entrez$entrezgene
original_gene_list <- merge_res_with_entrez$log2FoldChange

# name the vector
names(original_gene_list) <- merge_res_with_entrez$entrezgene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# gene set enrichment
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = Oaries, 
             pAdjustMethod = "BH")
write.csv(gse, file = paste0(prefix3, "-Gene-set-enrichment.csv"),row.names=F)

# when no FDR correction perfrmed, some enriched terms are found. Even after specifying p of 0.1 no enriched terms

pdf(paste0(prefix3,"-Dotplot.pdf"), width=10, height=12)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

# enrichment map
pdf(paste0(prefix3,"-Emapplot.pdf"), width=10, height=12)
x2 <- pairwise_termsim(gse)
emapplot(x2, showCategory = 10)
dev.off()

#category net plot
pdf(paste0(prefix3,"-Cnetplot.pdf"), width=10, height=12)
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
dev.off()

# Ridgeplot
pdf(paste0(prefix3,"-Ridgeplot.pdf"), width=10, height=12)
ridgeplot(gse) + labs(x = "enrichment distribution")
dev.off()


# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
pdf(paste0(prefix3,"-GSEAplot.pdf"))
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
dev.off()

# KEGG ENRICHMENT
# Create a vector of the gene unuiverse
kegg_gene_list <- merge_res_with_entrez$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- merge_res_with_entrez$entrezgene

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "oas"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")
write.csv(kk2, file = paste0(prefix3, "-Kegg-enrichment.csv"),row.names=F)

pdf(paste0(prefix3,"Kegg-Dotplot.pdf"), width=10, height=12)
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
dev.off()

pdf(paste0(prefix3,"Kegg-Emapplot.pdf"), width=10, height=12)
x2 <- pairwise_termsim(kk2)
emapplot(x2, showCategory = 10)
dev.off()

# categorySize can be either 'pvalue' or 'geneNum'
pdf(paste0(prefix3,"Kegg-Cnetplot.pdf"), width=10, height=12)
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
dev.off()

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
pdf(paste0(prefix3,"Kegg-Gseaplot.pdf"), width=10, height=12)
gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)
dev.off()


# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="oas04613", species = kegg_organism)

# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="oas04613", species = kegg_organism, kegg.native = F)
knitr::include_graphics("oas04613.pathview.png")
