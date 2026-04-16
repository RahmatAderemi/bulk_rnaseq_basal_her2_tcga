##GENE AND PATHWAY ENRICHMENT ANALYSIS
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)

#A. NORMAL VS BASAL-LIKE BREAST CANCER ENRICHMENT ANALYSIS
##1. Read the DGE data
Basal= read.csv("results/NormalvsBasal.csv", header = TRUE)
head(Basal)

##2. Identify significant genes
sig_basal <- Basal[
  !is.na(Basal$padj) &
    Basal$padj < 0.05 &
    abs(Basal$log2FoldChange) > 1,
]

sig_genes_basal <- sub("\\..*$", "", sig_basal$X)
universe_genes <- sub("\\..*$", "", Basal$X)

##3.Overrepresentation analysis

Upgenes <- read.csv("results/upregulated_Basal.csv", header= TRUE)
Downgenes <- read.csv("results/downregulated_Basal.csv", header = TRUE)

up_vec <- Upgenes$X
down_vec <- Downgenes$X

up_vec <- sub("\\..*$", "", up_vec)
down_vec <- sub("\\..*$", "", down_vec)

up_vec <- unique(na.omit(up_vec))
down_vec <- unique(na.omit(down_vec))

gene_list <- list(
  Upregulated = up_vec,
  Downregulated = down_vec
)

### Comparing Biology Process
cc <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP"
)

dotplot(cc, showCategory = 10, font.size= 6)

### Comparing Cellular Component
cc1 <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "CC"
)

dotplot(cc1, showCategory = 10, font.size= 6)

### Comparing Molecular Function
cc2 <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "MF"
)

dotplot(cc2, showCategory = 10, font.size= 6)

##4. Geneset enrichment analysis
Basal_genes = Basal$log2FoldChange

#name the vector
names(Basal_genes) <- Basal$X
Basal_genes <- na.omit(Basal_genes)
Basal_genes <- sort(Basal_genes, decreasing = TRUE)
head(Basal_genes)
genes_ensembl <- sub("\\..*$", "", names(Basal_genes))
names(Basal_genes) <- genes_ensembl

# Comparing BP
gseaBasal<-gseGO(geneList = Basal_genes,
                  ont = "BP",
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENSEMBL")
dotplot(gseaBasal, showCategory=10, font.size= 7, split=".sign") + facet_grid(.~.sign)

# Comparing CC
gseaBasal_cc<-gseGO(geneList = Basal_genes,
                 ont = "CC",
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENSEMBL")
dotplot(gseaBasal_cc, showCategory=10, font.size= 7, split=".sign") + facet_grid(.~.sign)


# Comparing MF
gseaBasal_mf<-gseGO(geneList = Basal_genes,
                    ont = "MF",
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENSEMBL")
dotplot(gseaBasal_mf, showCategory=10, font.size= 7, split=".sign") + facet_grid(.~.sign)

#B. NORMAL VS HER2 ENRICHED BREAST CANCER
##1. Read the DGE data
Her2= read.csv("results/NormalvsHER2.csv", header = TRUE)
head(Her2)

##2. Identify significant genes
sig_her2 <- Her2[
  !is.na(Her2$padj) &
    Her2$padj < 0.05 &
    abs(Her2$log2FoldChange) > 1,
]

sig_genes_Her2 <- sub("\\..*$", "", sig_her2$X)
universe_genes <- sub("\\..*$", "", Her2$X)

##3.Overrepresentation analysis

Upgenes <- read.csv("results/upregulated_HER2.csv", header= TRUE)
Downgenes <- read.csv("results/downregulated_HER2.csv", header = TRUE)

up_vec <- Upgenes$X
down_vec <- Downgenes$X

up_vec <- sub("\\..*$", "", up_vec)
down_vec <- sub("\\..*$", "", down_vec)

up_vec <- unique(na.omit(up_vec))
down_vec <- unique(na.omit(down_vec))

gene_list <- list(
  Upregulated = up_vec,
  Downregulated = down_vec
)

## Comparing ORA for Biological Process
cc <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP"
)

dotplot(cc, showCategory = 10, font.size= 6)

### Comparing Cellular Component
cc1 <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "CC"
)

dotplot(cc1, showCategory = 10, font.size= 6)

### Comparing Molecular Function
cc2 <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "MF"
)

dotplot(cc2, showCategory = 10, font.size= 6)

##4. Geneset enrichment analysis
Her2_genes = Her2$log2FoldChange

#name the vector
names(Her2_genes) <- Her2$X
Her2_genes <- na.omit(Her2_genes)
Her2_genes <- sort(Her2_genes, decreasing = TRUE)
head(Her2_genes)
genes_ensembl <- sub("\\..*$", "", names(Her2_genes))
names(Her2_genes) <- genes_ensembl

# Comparing BP
gseaHer2<-gseGO(geneList = Her2_genes,
                 ont = "BP",
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENSEMBL")

dotplot(gseaHer2, showCategory=10, font.size= 7, split=".sign") + facet_grid(.~.sign)

# Comparing CC
gseaHer2_cc<-gseGO(geneList = Her2_genes,
                ont = "CC",
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL")

dotplot(gseaHer2_cc, showCategory=10, font.size= 7, split=".sign") + facet_grid(.~.sign)

# Comparing MF
gseaHer2_mf<-gseGO(geneList = Her2_genes,
                ont = "MF",
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL")

dotplot(gseaHer2_mf, showCategory=10, font.size= 7, split=".sign") + facet_grid(.~.sign)

#C. HER2 ENRICHED VS BASAL-LIKE BREAST CANCER
##1. Read the DGE data
HvB= read.csv("results/HER2vsBasal.csv", header = TRUE)
head(HvB)

##2. Identify significant genes
sig_hvb <- HvB[
  !is.na(HvB$padj) &
    HvB$padj < 0.05 &
    abs(HvB$log2FoldChange) > 1,
]

sig_genes_HvB <- sub("\\..*$", "", sig_hvb$X)
universe_genes <- sub("\\..*$", "", HvB$X)

##3.Overrepresentation analysis

Upgenes <- read.csv("results/upregulated_HVB.csv", header= TRUE)
Downgenes <- read.csv("results/downregulated_HVB.csv", header = TRUE)

up_vec <- Upgenes$X
down_vec <- Downgenes$X

up_vec <- sub("\\..*$", "", up_vec)
down_vec <- sub("\\..*$", "", down_vec)

up_vec <- unique(na.omit(up_vec))
down_vec <- unique(na.omit(down_vec))

gene_list <- list(
  Upregulated = up_vec,
  Downregulated = down_vec
)

# Comparing BP
cc <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP"
)

dotplot(cc, showCategory = 10, font.size= 6)

# Comparing CC
cc1 <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "CC"
)

dotplot(cc1, showCategory = 10, font.size= 6)

# Comparing MF
cc2 <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "MF"
)

dotplot(cc2, showCategory = 10, font.size= 6)

##4. Geneset enrichment analysis
HvB_genes = HvB$log2FoldChange

#name the vector
names(HvB_genes) <- HvB$X
HvB_genes <- na.omit(HvB_genes)
HvB_genes <- sort(HvB_genes, decreasing = TRUE)
head(HvB_genes)
genes_ensembl <- sub("\\..*$", "", names(HvB_genes))
names(HvB_genes) <- genes_ensembl

# Comparing BP
gseaHvB<-gseGO(geneList = HvB_genes,
                ont = "BP",
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL")

dotplot(gseaHvB, showCategory=10, font.size= 7, split=".sign") + facet_grid(.~.sign)

# Comparing CC
gseaHv_CC<-gseGO(geneList = HvB_genes,
               ont = "CC",
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL")

dotplot(gseaHv_CC, showCategory=10, font.size= 7, split=".sign") + facet_grid(.~.sign)

# Comparing MF
gseaHvB_MF<-gseGO(geneList = HvB_genes,
               ont = "MF",
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL")

dotplot(gseaHvB_MF, showCategory=10, font.size= 7, split=".sign") + facet_grid(.~.sign)

sessionInfo()
# R version 4.4.3 (2025-02-28 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26200)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# time zone: Europe/London
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] org.Hs.eg.db_3.20.0    AnnotationDbi_1.68.0   IRanges_2.40.1         S4Vectors_0.44.0      
# [5] Biobase_2.66.0         BiocGenerics_0.52.0    ggplot2_4.0.2          enrichplot_1.26.6     
# [9] clusterProfiler_4.14.6
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.3.0               gson_0.1.0              rlang_1.2.0             magrittr_2.0.5         
# [5] DOSE_4.0.1              otel_0.2.0              compiler_4.4.3          RSQLite_2.4.6          
# [9] png_0.1-8               vctrs_0.6.5             reshape2_1.4.5          stringr_1.6.0          
# [13] pkgconfig_2.0.3         crayon_1.5.3            fastmap_1.2.0           XVector_0.46.0         
# [17] labeling_0.4.3          rmarkdown_2.31          UCSC.utils_1.2.0        purrr_1.2.0            
# [21] bit_4.6.0               xfun_0.57               zlibbioc_1.52.0         cachem_1.1.0           
# [25] aplot_0.2.9             GenomeInfoDb_1.42.3     jsonlite_2.0.0          blob_1.3.0             
# [29] BiocParallel_1.40.2     parallel_4.4.3          R6_2.6.1                stringi_1.8.7          
# [33] RColorBrewer_1.1-3      GOSemSim_2.32.0         Rcpp_1.1.1              knitr_1.51             
# [37] snow_0.4-4              ggtangle_0.1.1          R.utils_2.13.0          Matrix_1.7-2           
# [41] splines_4.4.3           igraph_2.2.3            tidyselect_1.2.1        qvalue_2.38.0          
# [45] rstudioapi_0.18.0       yaml_2.3.12             codetools_0.2-20        lattice_0.22-6         
# [49] tibble_3.3.0            plyr_1.8.9              treeio_1.30.0           withr_3.0.2            
# [53] KEGGREST_1.46.0         S7_0.2.1                evaluate_1.0.5          gridGraphics_0.5-1     
# [57] Biostrings_2.74.1       pillar_1.11.1           ggtree_3.14.0           ggfun_0.2.0            
# [61] generics_0.1.4          scales_1.4.0            tidytree_0.4.7          glue_1.8.0             
# [65] lazyeval_0.2.3          tools_4.4.3             data.table_1.17.8       fgsea_1.32.4           
# [69] fs_2.0.1                fastmatch_1.1-8         cowplot_1.2.0           grid_4.4.3             
# [73] tidyr_1.3.2             ape_5.8-1               colorspace_2.1-2        nlme_3.1-167           
# [77] GenomeInfoDbData_1.2.13 patchwork_1.3.2         cli_3.6.5               rappdirs_0.3.3         
# [81] dplyr_1.1.4             gtable_0.3.6            R.methodsS3_1.8.2       yulab.utils_0.2.4      
# [85] digest_0.6.39           ggrepel_0.9.8           ggplotify_0.1.3         farver_2.1.2           
# [89] memoise_2.0.1           htmltools_0.5.9         R.oo_1.27.1             lifecycle_1.0.5        
# [93] httr_1.4.8              GO.db_3.20.0            bit64_4.6.0-1