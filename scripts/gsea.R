##GENE AND PATHWAY ENRICHMENT ANALYSIS
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)

## BASAL ENRICHMENT ANALYSIS
Basal= read.csv("NormalvsBasal.csv", header = TRUE)
head(Basal)

##Identify significant genes
sig_basal <- Basal[
  !is.na(Basal$padj) &
    Basal$padj < 0.05 &
    abs(Basal$log2FoldChange) > 1,
]

sig_genes_basal <- sub("\\..*$", "", sig_basal$X)
universe_genes <- sub("\\..*$", "", Basal$X)

## Overrepresentaion analysis
enrichedBasal <- enrichGO(gene = sig_genes_basal,
                          universe= universe_genes,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENSEMBL",
                          ont= "ALL",
                          pAdjustMethod = "BH")
head(enrichedBasal)
dotplot(enrichedBasal,font.size=6, showCategory = 10, split= "ONTOLOGY") +facet_grid(~ONTOLOGY)
write.csv(enrichedBasal, file= "Basal_ORA.csv")

enrichedBasal_BP <- enrichGO(gene = sig_genes_basal,
                             universe = universe_genes,
                             OrgDb = org.Hs.eg.db,
                             keyType = "ENSEMBL",
                             ont = "BP",
                             pAdjustMethod = "BH")
dotplot(enrichedBasal_BP, showCategory= 15, font.size= 6)
ORAplotBasal <- pairwise_termsim(enrichedBasal)
emapplot(ORAplotBasal)
dotplot(enrichedBasal_BP, font.size= 6)

##Geneset enrichment analysis
Basal_genes = Basal$log2FoldChange

#name the vector
names(Basal_genes) <- Basal$X
Basal_genes <- na.omit(Basal_genes)
Basal_genes <- sort(Basal_genes, decreasing = TRUE)
head(Basal_genes)
genes_ensembl <- sub("\\..*$", "", names(Basal_genes))
names(Basal_genes) <- genes_ensembl

gseaBasal <-gseGO(geneList = Basal_genes,
                  ont = "ALL",
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENSEMBL")
dotplot(gseaBasal, showCategory=10, font.size= 8, split= "ONTOLOGY") +facet_grid(~ONTOLOGY)
gseaBasal_result <- gseaBasal@result[, c("ID", "Description", "NES", "p.adjust")]
write.csv(gseaBasal_result, file = "GenesetAnalysis_GOBasal.csv")

gseaBasal_BP <- gseGO(
  geneList = Basal_genes,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL"
)
gseaBasal_BPresult <- gseaBasal_BP@result[, c("ID", "Description", "NES", "p.adjust")]
write.csv(gseaBasal_BPresult, file = "GSEABasal_BP.csv")
dotplot(gseaBasal_BP, showCategory = 15, font.size= 8, split=".sign") + facet_grid(.~.sign)

##NormalvsHER2+
her2 <- read.csv("HER2.csv")

##Identify significant genes
sig_her2 <- her2[
  !is.na(her2$padj) &
    her2$padj < 0.05 &
    abs(her2$log2FoldChange) > 1,
]

sig_genes_her2 <- sub("\\..*$", "", sig_her2$X)
universe_genes <- sub("\\..*$", "", her2$X)

##Over-representation analysis
ora_her2 <- enrichGO(gene = sig_genes_her2,
                     universe = universe_genes,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENSEMBL",
                     ont = "BP")
head(ora_her2)
dotplot(ora_her2, showCategory = 15, font.size= 6)
enrichedHER <- enrichGO(gene =sig_genes_her2,
                        universe = universe_genes,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENSEMBL",
                        ont = "ALL")
dotplot(enrichedHER, font.size= 6, split= "ONTOLOGY") + facet_grid(~ONTOLOGY)

##Gene Set Enrichment analysis
her2_genes <-her2$log2FoldChange 
names(her2_genes) <- her2$X 
her2_genes <- na.omit(her2_genes) 
her2_genes <- sort(her2_genes, decreasing = TRUE) 
names(her2_genes) <- sub("\\..*$", "", names(her2_genes))


gseaher_BP <- gseGO(
  geneList = her2_genes,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL"
)
gseaHER_BPresult <- gseaher_BP@result[, c("ID", "Description", "NES", "p.adjust")]
dotplot(gseaher_BP, showCategory = 10, font.size= 7, split=".sign") + facet_grid(.~.sign)

##HER2+ vs Basal
herbasal <- read.csv("HER2vsBasal.csv")
head(herbasal)
##Identify significant genes
sig_hvb <- herbasal[
  !is.na(herbasal$padj) &
    herbasal$padj < 0.05 &
    abs(herbasal$log2FoldChange) > 1,
]

sig_genes_hvb <- sub("\\..*$", "", sig_hvb$X)
universe_genes <- sub("\\..*$", "", herbasal$X)

## Overrepresentaion analysis
enrichedhvb <- enrichGO(gene = sig_genes_hvb,
                        universe= universe_genes,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENSEMBL",
                        ont= "ALL",
                        pAdjustMethod = "BH")
head(enrichedhvb)
dotplot(enrichedhvb,font.size=6, showCategory = 10, split= "ONTOLOGY") +facet_grid(~ONTOLOGY)


enrichedhvb_BP <- enrichGO(gene = sig_genes_hvb,
                           universe = universe_genes,
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENSEMBL",
                           ont = "BP",
                           pAdjustMethod = "BH")
dotplot(enrichedhvb_BP, showCategory= 15, font.size= 6)

##Geneset enrichment analysis
HervsBasal = herbasal$log2FoldChange

#name the vector
names(HervsBasal) <- herbasal$X
HervsBasal <- na.omit(HervsBasal)
HervsBasal <- sort(HervsBasal, decreasing = TRUE)
head(HervsBasal)
genes_ensembl <- sub("\\..*$", "", names(HervsBasal))
names(HervsBasal) <- genes_ensembl

gseahvb <-gseGO(geneList = HervsBasal,
                ont = "ALL",
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL")
dotplot(gseahvb, showCategory=10, font.size= 8, split= "ONTOLOGY") +facet_grid(~ONTOLOGY)
gseahvb_result <- gseahvb@result[, c("ID", "Description", "NES", "p.adjust")]
write.csv(gseahvb_result, file = "GenesetAnalysis_HervsBasal.csv")

gseahvb_BP <-gseGO(geneList = HervsBasal,
                   ont = "BP",
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENSEMBL")
dotplot(gseahvb_BP, showCategory=10, font.size= 8, split= ".sign") +facet_grid(.~.sign)
