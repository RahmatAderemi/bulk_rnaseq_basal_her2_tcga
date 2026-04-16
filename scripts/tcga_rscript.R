# DATA DOWNLOAD

library(TCGAbiolinks)
library(maftools)
library(pheatmap)
library(tidyverse)
library(SummarizedExperiment)

## Find the projects available on GDC
projects <- getGDCprojects()
# View(projects)

## Look at project summary for the TCGA-BRCA project
summary <- getProjectSummary("TCGA-BRCA")
summary

## Query GDC for the Transcriptomic data from TCGA-BRCA
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Transcriptome Profiling",
                  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "STAR - Counts",
                  access = "open")

output_query <- getResults(query)

## Download data
GDCdownload(query)

## Prepare data
tcga_obj <- GDCprepare(query = query, summarizedExperiment = TRUE)

# Moving the TCGA data into another variable so any changes I make does not affect the initial GDC object.
data <- tcga_obj

# SELECTING DATA FOR NORMAL BREAST TISSUE, BASAL AND HER2+
colData(tcga_obj)

### substr(sampleIDs, 1,12 because the expression data has 12 unique identifiers while the subtype data only the patients IDs)
sample_ids <- colnames(data)
patient_ids <- substr(sample_ids, 1, 12)

## First get normal and tumour data type separated in the sample type variable. 
## Then create a tibble of patient IDs and their PAM50 subtypes
sample_type <- colData(data)$sample_type
subtype <- TCGAquery_subtype("BRCA")
subtype <- subtype[, c("patient", "BRCA_Subtype_PAM50")]

##Create a variable condition to that segregate the samples into normal, basal and her2
condition <- ifelse(
  sample_type == "Solid Tissue Normal", 
  "Normal",
  ifelse(
    patient_ids %in% subtype$patient[subtype$BRCA_Subtype_PAM50 == "Basal"],
    "Basal",
    ifelse(
      patient_ids %in% subtype$patient[subtype$BRCA_Subtype_PAM50 == "Her2"],
      "HER2",
      NA
    )
  )
)

colData(data)$condition <- condition

# EXPLORATORY DATA ANALYSIS
table(colData(data)$condition)

##1. Remove NAs
filtered_data <- data[, !is.na(colData(data)$condition)]

#2. Convert the conditions to factors and relevel them so Normal becomes the reference.
filtered_data$condition <- as.factor(filtered_data$condition)
levels(filtered_data$condition)

filtered_data$condition <- relevel(filtered_data$condition, ref = "Normal")
filtered_data$condition


#3. Convert data to a DESeq object and remove low expression counts.
library(DESeq2)
ddSE <- DESeqDataSet(filtered_data, design = ~condition)

keep <- rowSums(counts(ddSE))>=10
ddSE <- ddSE[keep,]

ddSE <- readRDS("ddSE.rds")
#4. Normalization
dds<- estimateSizeFactors(ddSE)
saveRDS(dds, "dds.rds")
plot(sizeFactors(dds), colSums(counts(dds)), 
     xlab = "Size Factors", 
     ylab= "Total Counts", 
     main = "QC: Size Factors Vs Total Counts")
abline(lm(colSums(counts(dds))~ sizeFactors(dds)))

#5. Variance stabilization
vst <- vst(dds, blind = FALSE)

library(BiocManager)
BiocManager::install("vsn")
library(vsn)

norm_count <- counts(dds, normalized = TRUE)
meanSdPlot(norm_count)

meanSdPlot(assay(vst))

##6. PCA
plotPCA(vst, intgroup= "condition")

##7. Hierachical clustering
plot(hclust(dist(t(assay(vst)))), labels= colData(vst)$condition)


# DIFFERENTIAL GENE EXPRESSION ANALYSIS
ddsSE <- DESeq(ddSE)

resultsNames(ddsSE)
#[1] "Intercept"                 "condition_Basal_vs_Normal" "condition_HER2_vs_Normal" 

## 1. Normal Vs Basal
NvB <- results(ddsSE, name = "condition_Basal_vs_Normal")
NvB

all(rownames(NvB) == rownames(ddsSE))
NvB$gene_name <- rowData(ddsSE)$gene_name
summary(NvB)
write.csv(NvB, file = "NormalvsBasal.csv")

library(EnhancedVolcano)
EnhancedVolcano(NvB, lab= NvB$gene_name, x= "log2FoldChange", y= "padj", title = "Normal vs Basal")

upregulated_Basal <- NvB[!is.na(NvB$padj) & NvB$log2FoldChange >1 & NvB$padj < 0.01,]
write.csv(upregulated_Basal, file = "upregulated_Basal.csv")

downregulated_Basal <- NvB[!is.na(NvB$padj) & NvB$log2FoldChange < 0 & NvB$padj < 0.01,]
write.csv(downregulated_Basal, file = "downregulated_Basal.csv")

## HEATMAP
###1. Significant genes
sig_basal <- rownames(NvB)[
  !is.na(NvB$padj) &
    NvB$padj < 0.01 &
    abs(NvB$log2FoldChange) > 1
]

###Draw heatmap
library(pheatmap)
table(colData(vst)$condition)
keep <- colData(vst)$condition %in% c("Normal", "Basal")
vst_nb <- vst[, keep]
mat_nb <- assay(vst_nb)
mat_nb_basal <- mat_nb[sig_basal, ]
rownames(mat_nb_basal) <- NvB[sig_basal, "gene_name"]

sample_subtype <- colData(vst_nb)$condition

ann_col <- data.frame(
  Subtype = factor(sample_subtype, levels = c("Normal", "Basal"))
)

rownames(ann_col) <- colnames(mat_nb_basal)
pheatmap(
  mat_nb_basal[1:50, ],
  fontsize = 5,
  annotation_col = ann_col,
  annotation_colors = list(
    Subtype = c(Normal = "#4DAF4A", Basal = "#E41A1C")
  ),
  show_colnames = FALSE
)

#Normal vs HER2
NvH <- results(ddsSE, name = "condition_HER2_vs_Normal")
NvH
write.csv(NvH, file = "NormalvsHER2.csv")
summary(NvH)

all(rownames(NvH) == rownames(ddsSE))
NvH$gene_name <- rowData(ddsSE)$gene_name

upregulated_HER2 <- NvH[!is.na(NvH$padj) & NvH$log2FoldChange > 1 & NvH$padj < 0.01,]
write.csv(upregulated_HER2, file = "upregulated_HER2.csv")

downregulated_HER2 <- NvH[!is.na(NvH$padj) & NvH$log2FoldChange < 0 & NvH$padj < 0.01,]
write.csv(downregulated_HER2, file = "downregulated_HER2.csv")

##Volcano plot
EnhancedVolcano(NvH, lab = NvH$gene_name, x= "log2FoldChange", y= "padj", title = "Normal vs HER2+")

## HEATMAP
###1. Significant genes
sig_her2 <- rownames(NvH)[
  !is.na(NvH$padj) &
    NvH$padj < 0.01 &
    abs(NvH$log2FoldChange) > 1
]

###Draw heatmap of top 50 significant genes.
keep2 <- colData(vst)$condition %in% c("Normal", "HER2")
vst_nh <- vst[, keep2]
mat_nh <- assay(vst_nh)
mat_nh_her2 <- mat_nh[sig_her2, ]
rownames(mat_nh_her2) <- NvH[sig_her2, "gene_name"]

sample_subtype <- colData(vst_nh)$condition

ann_col2 <- data.frame(
  Subtype = factor(sample_subtype, levels = c("Normal", "HER2"))
)

rownames(ann_col2) <- colnames(mat_nh_her2)
pheatmap(
  mat_nh_her2[1:50, ],
  fontsize = 5,
  annotation_col = ann_col2,
  annotation_colors = list(
    Subtype = c(Normal = "#4DAF4A", HER2 = "#E41A1C")
  ),
  show_colnames = FALSE
)

## HER2vsBasal

HvB <- results(ddsSE, contrast = c('condition', "HER2", "Basal"))
summary(HvB)
all(rownames(HvB) == rownames(ddsSE))
HvB$gene_name <- rowData(ddsSE)$gene_name
write.csv(HvB, file = "HER2vsBasal.csv")

upregulated_HVB <- HvB[!is.na(HvB$padj) & HvB$log2FoldChange > 1 & HvB$padj < 0.01,]
write.csv(upregulated_HVB, file = "upregulated_HVB.csv")

downregulated_HVB <- HvB[!is.na(HvB$padj) & HvB$log2FoldChange < 0 & HvB$padj < 0.01,]
write.csv(downregulated_HVB, file = "downregulated_HVB.csv")

##Volcano plot comparing HER2+ and Basal Subtypes.
EnhancedVolcano(HvB, lab = HvB$gene_name, x= "log2FoldChange", y="padj", title= "HER2 vs Basal")

## HEATMAP
###1. Significant genes
sig_hvb <- rownames(HvB)[
  !is.na(HvB$padj) &
    HvB$padj < 0.01 &
    abs(HvB$log2FoldChange) > 1
]

###Draw heatmap of top 50 significant genes.
keep3 <- colData(vst)$condition %in% c("HER2", "Basal")
vst_hb <- vst[, keep3]
mat_hb <- assay(vst_hb)
mat_hb_her2vsbasal <- mat_hb[sig_hvb, ]
rownames(mat_hb_her2vsbasal) <- HvB[sig_hvb, "gene_name"]

sample_subtype <- colData(vst_hb)$condition

ann_col3 <- data.frame(
  Subtype = factor(sample_subtype, levels = c("HER2", "Basal"))
)

rownames(ann_col3) <- colnames(mat_hb_her2vsbasal)
pheatmap(
  mat_hb_her2vsbasal[1:50, ],
  fontsize = 5,
  annotation_col = ann_col3,
  annotation_colors = list(
    Subtype = c(HER2 = "#4DAF4A", Basal = "#E41A1C")
  ),
  show_colnames = FALSE
)

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
#   [1] EnhancedVolcano_1.24.0      ggrepel_0.9.8               vsn_3.74.0                 
# [4] lubridate_1.9.4             forcats_1.0.1               stringr_1.6.0              
# [7] dplyr_1.1.4                 purrr_1.2.0                 readr_2.1.6                
# [10] tidyr_1.3.1                 tibble_3.3.0                ggplot2_4.0.1              
# [13] tidyverse_2.0.0             pheatmap_1.0.13             maftools_2.22.0            
# [16] DESeq2_1.46.0               SummarizedExperiment_1.36.0 Biobase_2.66.0             
# [19] GenomicRanges_1.58.0        GenomeInfoDb_1.42.3         IRanges_2.40.1             
# [22] S4Vectors_0.44.0            BiocGenerics_0.52.0         MatrixGenerics_1.18.1      
# [25] matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.6            lattice_0.22-6          tzdb_0.5.0              vctrs_0.6.5            
# [5] tools_4.4.3             generics_0.1.4          parallel_4.4.3          pkgconfig_2.0.3        
# [9] Matrix_1.7-2            data.table_1.17.8       RColorBrewer_1.1-3      S7_0.2.1               
# [13] lifecycle_1.0.5         GenomeInfoDbData_1.2.13 compiler_4.4.3          farver_2.1.2           
# [17] statmod_1.5.1           codetools_0.2-20        preprocessCore_1.68.0   hexbin_1.28.5          
# [21] pillar_1.11.1           crayon_1.5.3            BiocParallel_1.40.0     limma_3.62.2           
# [25] affy_1.84.0             DelayedArray_0.32.0     abind_1.4-8             tidyselect_1.2.1       
# [29] locfit_1.5-9.12         stringi_1.8.7           labeling_0.4.3          splines_4.4.3          
# [33] grid_4.4.3              colorspace_2.1-2        cli_3.6.5               SparseArray_1.6.2      
# [37] magrittr_2.0.4          S4Arrays_1.6.0          survival_3.8-3          withr_3.0.2            
# [41] scales_1.4.0            UCSC.utils_1.2.0        timechange_0.3.0        XVector_0.46.0         
# [45] httr_1.4.8              affyio_1.76.0           hms_1.1.4               DNAcopy_1.80.0         
# [49] rlang_1.1.6             Rcpp_1.1.0              glue_1.8.0              BiocManager_1.30.27    
# [53] rstudioapi_0.18.0       jsonlite_2.0.0          R6_2.6.1                zlibbioc_1.52.0


