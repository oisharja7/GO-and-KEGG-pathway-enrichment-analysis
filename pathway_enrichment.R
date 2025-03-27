# install and load packages 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))

library(clusterProfiler)
library(org.Hs.eg.db)   # Human gene annotation
library(enrichplot)     # Visualization for enrichment analysis
library(DOSE)           # For additional enrichment analysis functions

# prim vs norm 

# working directory 
setwd("D:/DEG analysis_oisharja/Significant csvs/miRNA Target Predictions")

# load 
common_genes_prim_norm <- read.csv("Common_Genes_prim_norm.csv", stringsAsFactors = FALSE)

# 
gene_vector_1 <- common_genes_prim_norm$Gene

# trim extra space 
gene_vector_1 <- trimws(gene_vector_1)

# convert gene symbol to entrez ID
gene_ids_prim_norm <- bitr(gene_vector_1,
                           fromType = "SYMBOL",
                           toType = "ENTREZID", 
                           OrgDb = "org.Hs.eg.db")

head(common_genes_prim_norm)

# extract entrex ID vector for enrichment 
entrez_list_1 <- gene_ids_prim_norm$ENTREZID

# GO enrichment analysis for biological process 
ego_1 <- enrichGO(gene = entrez_list_1,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)

head(ego_1)


# KEGG pathway enrichment analysis 
ekegg_1 <- enrichKEGG(gene = entrez_list_1,
                      organism = "hsa",
                      pvalueCutoff = 0.05)

head(ekegg_1)


# visualise GO enrichment result 
dotplot(ego_1, showCategory = 10) + ggtitle("GO Biological Process Enrichment of Primary Ovarian Cancer vs Normal Ovary")

# visualise KEGG enrichment result 
barplot(ekegg_1, showCategory = 10) + ggtitle("KEGG Pathway Enrichment of Primary Ovarian Cancer vs Normal Ovary")


# rec vs norm 


# load 
common_genes_rec_norm <- read.csv("common_genes_rec_norm_df.csv", stringsAsFactors = FALSE)

# extract vector 
gene_vector_2 <- common_genes_rec_norm$Gene

# trim extra space 
gene_vector_2 <- trimws(gene_vector_2)

# convert gene symbol to entrez ID
gene_ids_rec_norm <- bitr(gene_vector_2,
                           fromType = "SYMBOL",
                           toType = "ENTREZID", 
                           OrgDb = "org.Hs.eg.db")

head(common_genes_prim_norm)

# extract entrez ID vector for enrichment 
entrez_list_2 <- gene_ids_rec_norm$ENTREZID

# GO enrichment analysis for biological process 
ego_2 <- enrichGO(gene = entrez_list_2,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)

head(ego_2)


# KEGG pathway enrichment analysis 
ekegg_2 <- enrichKEGG(gene = entrez_list_2,
                      organism = "hsa",
                      pvalueCutoff = 0.05)

head(ekegg_2)


# visualise GO enrichment result 
dotplot(ego_2, showCategory = 10) + ggtitle("GO Biological Process Enrichment of Recurrent Ovarian Cancer vs Normal Ovary")

# visualise KEGG enrichment result 
barplot(ekegg_2, showCategory = 10) + ggtitle("KEGG Pathway Enrichment of Recurrent Ovarian Cancer vs Normal Ovary")

