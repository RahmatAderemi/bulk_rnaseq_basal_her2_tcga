I collected STAR counts data from TCGA for Normal breast tissue, Basal-like breast cancer and HER2 enriched breast cancer.

I want to compare the differential gene expression levels and gene set enrichments in this order: 
* Normal breast tissue Vs HER2 enriched breast cancer
* Normal breast tissue Vs Basal-like breast cancer
* HER2 enriched breast cancer Vs Basal-like breast cancer

I picked these two BC PAM50 subtypes because of their immune evasion strategies and I wanted to see the different strategies they employ in immune evasion.

# Data Collection-
Using TCGA Biolinks, I queried TCGA for STAR counts from Transcriptomic profiling data for Primary Tumour and Solid Normal Tissue. These samples were collected from the TCGA-BRCA project.
I downloaded the samples and used GDC prepare to get all the metadata concerning the samples as well.

# Data Segregation- 
The data from TCGA-BRCA project included Primary Tumour PAM50 subtypes- LumA, LumB, Normal-like, Basal and Her2.
I wanted to exclude the other subtypes and only focus on those I want to analyse. To do this, I created patient's ID variable that extracts the first 12 unique identifier from the sample IDs (column names of the TCGA summarised experiment data).
Then I created a subtype variable using the TCGAquery_subtype function and asked for the BRCA subtypes. This returned an output of all the BRCA_PAM50 subtypes. Then I created a tibble with the patient's IDs and their PAM50 subtypes.
Using this information, I was able to separate the Basal and Her2 samples. To get the normal breast tissue, I just used subset samples whose sample type was "Solid Normal Tissue" from "Primary Tumor". 
Using ifelse statement, I was able to create a variable condition which grouped all patient IDs with sample type "Solid Tissue Normal" as Normal, patient IDs with BRCA_PAM50_subtype Basal and HER2 as Basal and Her2. The rest were classified NA.
The condition table was then added to my samples' column data.

# Exploratory data analysis-
First step was to remove all the NAs. After removing the NAs from my data, I converted the condition variable to factors. Then made Normal the reference level.
I converted the data into a DESEqDataSet then began to filter out low count rows. I only kept genes with row sums >= 10 as I believe they are the informative genes.
## Normalisation- I used DESEq2 to estimate size factors and use this information to normalise my data. This helps remove sequencing depth bias and ensures counts are on the same scale.
After estimating size factors, I plotted my size factors against the total gene counts for each samples. They mostly followed the expected trend but there were a few outliers with small size factors and high counts. The plot for comparing size factors against the Total counts was saved as Size_Factors_Vs_Total_Counts.

## Variance stabilisation- I used the vst function to stabilise the variance around the mean. To confirm my VST did what it was supposed to do, I plotted a mean vs SD plot of the normalised raw counts first. The SD was flat for most genes and then there was a spike in genes with high gene expression. This showed that the variance exploded for highly expressed genes. 
The data was strongly heteroscedastic with highly expressed genes dominating everything and the low expressed genes were compressed to near zero.
I also plotted a meanSdplot after vst. After vst, I observed that the SD spread is a bit more even, the low expression genes are more spread out and the extreme spike is gone. 

## Principal Component Analysis (PCA)- I carried out PCA on the variance stabilised counts and used the condition as the group of interest. PC1 accounted for 28% of the total variance observed in the data.
PC1 was able to clearly separate the normal breast tissue data from the diseased state. PC2 accounted for 12% of the total variance and it separared Basal and HER2 data but there are some data points from the Basal data that were seen clustering with HER2 data points.
The figure for the PCA was saved at Principal_Component_Analysis.

## Hierarchical clustering- I observed some overlap between the HER2 and Basal data as observed in the PCA. Plot was saved a  Cluster_Dendogram.

# Differential Gene Expression Analysis-
Using DESeq, I ran DGE on the samples. First result I checked was Normal vs Basal. The result was:
adjusted p-value < 0.1
LFC > 0 (up)       : 22788, 41%
LFC < 0 (down)     : 10760, 20%
outliers [1]       : 0, 0%
low counts [2]     : 8547, 16%

Normal vs HER2
out of 55095 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 14669, 27%
LFC < 0 (down)     : 13040, 24%
outliers [1]       : 0, 0%
low counts [2]     : 9615, 17%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

HER2 vs Basal
out of 55095 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 8002, 15%
LFC < 0 (down)     : 20237, 37%
outliers [1]       : 0, 0%
low counts [2]     : 11751, 21%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

# GENE SET ENRICHMENT ANALYSIS

Using Enhanced Volcano, I constructed a volcano plot and I used pheatmap package to construct the heatmaps.