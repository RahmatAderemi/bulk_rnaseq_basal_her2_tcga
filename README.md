# TRANSCRIPTOMIC PROFILING AND GENE SET ENRICHMENT ANALYSIS OF BASAL-LIKE BC AND HER2-ENRICHED BC.
## INTRODUCTION
Breast cancer (BC) is one of the most prevalent cancers affecting women accounting for about 30% of new cancer cases yearly. Next to lung cancer, BC is the highest cause of cancer-related deaths in women with a mortality rate of 2.3% (American Cancer Society, 2025). The molecular characteristics and clinical manifestations of BC are varied, and these inform the prognosis and treatment options.
The PAM50 assay is a 50-gene signature test that classifies breast cancer into five intrinsic molecular subtypes— Normal-like, Luminal A, Luminal B, HER2-enriched, and Basal-like—based on gene expression patterns.
This study will focus on analysing Basal-like and HER2 enriched BC data collected from the TCGA-BRCA project because these subtypes exhibit immune cell infiltration and I want to see how they do this.
I will be comparing the differential gene expression levels and gene set enrichments in this order: 
* Normal breast tissue Vs Basal-like breast cancer
* Normal breast tissue Vs HER2 enriched breast cancer
* HER2 enriched breast cancer Vs Basal-like breast cancer

## METHODOLOGY
### Data Collection-
Using TCGA Biolinks, I queried TCGA for STAR counts from Transcriptomic profiling data for Primary Tumour and Solid Normal Tissue. These samples were collected from the TCGA-BRCA project.
I downloaded the samples and used GDC prepare to get all the metadata concerning the samples as well.

### Data Segregation- 
The data from TCGA-BRCA project included Primary Tumour PAM50 subtypes- LumA, LumB, Normal-like, Basal and Her2.
I wanted to exclude the other subtypes and only focus on those I want to analyse. To do this, I created patient's ID variable that extracts the first 12 unique identifier from the sample IDs (column names of the TCGA summarised experiment data).
Then I created a subtype variable using the TCGAquery_subtype function and asked for the BRCA subtypes. This returned an output of all the BRCA_PAM50 subtypes. Then I created a tibble with the patient's IDs and their PAM50 subtypes.
Using this information, I was able to separate the Basal and Her2 samples. To get the normal breast tissue, I just used subset samples whose sample type was "Solid Normal Tissue" from "Primary Tumor". 
Using ifelse statement, I was able to create a variable condition which grouped all patient IDs with sample type "Solid Tissue Normal" as Normal, patient IDs with BRCA_PAM50_subtype Basal and HER2 as Basal and Her2. The rest were classified NA.
The condition table was then added to my samples' column data.

### Exploratory data analysis-
First step was to remove all the NAs. After removing the NAs from my data, I converted the condition variable to factors. Then made Normal the reference level.
I converted the data into a DESEqDataSet then began to filter out low count rows. I only kept genes with row sums >= 10 as I believe they are the informative genes.

#### Normalisation- 
I used DESEq2 to estimate size factors and use this information to normalise my data. This helps remove sequencing depth bias and ensures counts are on the same scale.
After estimating size factors, I plotted my size factors against the total gene counts for each samples. They mostly followed the expected trend but there were a few outliers with small size factors and high counts. The plot for comparing size factors against the Total counts was saved as Size_Factors_Vs_Total_Counts.

#### Variance stabilisation- 
I used the vst function to stabilise the variance around the mean. To confirm my VST did what it was supposed to do, I plotted a mean vs SD plot of the normalised raw counts first. The SD was flat for most genes and then there was a spike in genes with high gene expression. This showed that the variance exploded for highly expressed genes. 
The data was strongly heteroscedastic with highly expressed genes dominating everything and the low expressed genes were compressed to near zero.
I also plotted a meanSdplot after vst. After vst, I observed that the SD spread is a bit more even, the low expression genes are more spread out and the extreme spike is gone. 

#### Principal Component Analysis (PCA)- 
I carried out PCA on the variance stabilised counts and used the condition as the group of interest. PC1 accounted for 28% of the total variance observed in the data.
PC1 was able to clearly separate the normal breast tissue data from the diseased state. PC2 accounted for 12% of the total variance and it separared Basal and HER2 data but there are some data points from the Basal data that were seen clustering with HER2 data points.
The figure for the PCA was saved at Principal_Component_Analysis.

#### Hierarchical clustering- 
I observed some overlap between the HER2 and Basal data as observed in the PCA. Plot was saved a  Cluster_Dendogram.

### Differential Gene Expression Analysis-
Using DESeq, I ran DGE on the samples. I classified genes with logfoldchange >1 and padjust < 0.01 as upregulated and genes with logfoldchange < 0 and padjust < 0.01as downregulated.
Using Enhanced Volcano, I constructed a volcano plot and I used pheatmap package to construct the heatmaps.

### Gene Set Enrichment Analysis-
Using clusterProfile, I carried out an overrepresentation analysis and gene set enrichment analysis of the toptables from NormalvsBasal, NormalvsHER2 and HER2vsBasal.
I visualised the upregulated and downregulated biological processes from the overrepresentation analysis using dotplots and I also visualised the activated and suppressed gene sets using dotplots.

### Immune deconvolution-
I looked at the immune landscape of the tumour samples using the immunedeconv package that takes bulk RNA-seq data and try to deconvolute the immune landscape. I used EPIC (Experimental Predictive Immunoinformatics Cell-type deconvolution) to try to get the characteristics of the immune cells in the tumour microenvironment.
I chose EPIC as the method because it uses a pre-defined "signature matrix" (profiles of known immune and stromal cells) and a mathematical process called constrained least-squares regression to estimate the fraction of each cell type in the sample.
This was my first time using immunedeconv so this code was mostly written by Google gemini. 

## RESULTS
### Normal Vs Basal
Out of the 55095 genes compared between normal breast tissue and basal-line breast cancer, 41% were upregulated and 20% were downregulated.
Some genes that were downregulated in BLBC include PDK4, MYOC,HSPB6,PLIN1 and LEP. Some upregulated genes in BLBC include FTHL17, DBNDD1, TLX1, MMP1 and DNMT3B.

#### Over-representation Analysis (ORA)-
Using the over-representation analysis, I compared Molecular functions, Cellular components and biological processes of the toptable gotten from the differential gene expression analysis.
The results paint a clear picture of a highly proliferative, immune-infiltrated, and metabolically reprogrammed environment.
##### Cellular Components
This plot identifies where the gene products are active. There is a massive enrichment for nuclear and chromosomal components (nucleosomes, condensed chromosomes, centromeric regions). This indicates high levels of DNA organization and replication. There are also immune-related complexes (T-cell receptor and immunoglobulin complexes), suggesting a significant presence of immune cells within the tumor microenvironment.
The genes associated with muscle and contractile structures (sarcomeres, myofibrils, Z-discs) are downregulated. This likely reflects the loss of normal myoepithelial cell function or the displacement of normal breast tissue architecture by the tumor.
##### Molecular functions 
This identifies the biochemical activities of the genes. The dominant upregulated functions involve DNA and protein handling, such as DNA helicase activity (unzipping DNA for replication) and serine-type endopeptidase activity (enzymes that break down proteins). The presence of "structural constituent of skin epidermis" is a classic hallmark of the "basal-like" subtype, which expresses genes typically found in basal epithelial layers.
There is a significant loss of hormone and nuclear receptor activity. Since basal-like cancers are frequently "triple-negative," they lack the estrogen and progesterone receptors found in normal breast tissue. There is also a decrease in metabolic enzymes (oxidoreductases and aldehyde dehydrogenases).
##### Biological processes
This describes the larger biological "programs" being turned on or off. The upregulated results are dominated by the Cell Cycle. Processes like DNA replication, chromosome segregation, and nuclear division are highly active, which is expected in aggressive, fast-growing basal-like tumors. "Keratinization" is also upregulated, reinforcing the basal/squamous identity of these cells.
The tumor has suppressed normal metabolic pathways, specifically fatty acid and lipid oxidation. Additionally, processes related to muscle contraction and vascular maintenance (regulating blood vessel diameter) are turned down, indicating a shift away from the homeostatic functions of normal breast tissue.

#### Gene set enrichment analysis (GSEA)-
##### Cellular components
Strong enrichment is seen in immune-related complexes (T cell receptor and immunoglobulin complexes) and nuclear components (nucleosomes, kinetochores, and condensed chromosomes). This indicates a tumor microenvironment with high immune infiltration and intense nuclear activity for replication.
The suppressed components are largely related to muscle and structural integrity, including the sarcolemma, Z disc, and I band. Peroxisomes and microbodies, involved in metabolic breakdown, are also notably suppressed.
##### Molecular functions 
The activated functions are centered on structural identity (structural constituent of skin epidermis) and DNA handling (DNA helicase activity). There is also significant antigen binding activity, further supporting the presence of active immune responses within the tumor.
There is a significant loss of hormone and nuclear receptor activity, as well as various oxidoreductase activities. The suppression of lipid transporter and transmembrane transporter binding points toward a breakdown in normal cellular transport and communication.
##### Biological processes
The dominant theme is cell division and proliferation. Processes such as chromosome separation, nucleosome assembly, and mitotic sister chromatid segregation are highly active. Keratinization is also a key activated feature, which is a hallmark of the basal-like identity.
There is a massive downregulation of metabolic processes, specifically fatty acid beta-oxidation and lipid oxidation. Additionally, pathways related to muscle contraction and renal system processes are suppressed, reflecting the loss of normal tissue homeostatic functions.

The results characterize BLBC as a highly aggressive state where the cell has abandoned normal metabolic "housekeeping" (like fatty acid oxidation) and structural identity (muscle-related components) in favor of rapid, uncontrolled growth (cell cycle and DNA replication) and basal epithelial identity (keratinization). The data also highlights a significant immune response within these tumors that is absent in normal tissue.

### Normal vs HER2
Out of the 55,095 genes analysed between normal breast tissue and HER2 enriched BC, a total of 27,709 genes (approximately 51% of the analyzed transcriptome) are significantly differentially expressed at an adjusted p-value of < 0.1. This indicates that HER2-enriched tumors represent a profound biological departure from normal breast tissue.
There is a slightly higher proportion of upregulated genes (27%) compared to downregulated genes (24%).
Some of the upregulated genes include TPX2, MMP11, COX7B2 etc. Some downregulated genes include MYOC, VEGFD,EZH1,TNMD, MEOX1 etc.
#### Over-representation analysis-
The overrepresentation analysis (ORA) comparing Normal breast tissue and HER2-enriched breast cancer reveals a biological landscape marked by hyperactive cell division, a strong immune presence, and a significant loss of normal tissue structural markers.
##### Cellular Component
High enrichment is seen in nuclear and immune structures. This includes nucleosomes, condensed chromosomes, and centromeric regions (related to mitosis), alongside immunoglobulin complexes and T cell receptor complexes. Notably, collagen-containing extracellular matrix is also upregulated, suggesting significant tissue remodeling around the tumor.
The suppressed components are structural elements of muscle and contractile fibers, such as the sarcomere, myofibril, I band, and Z disc. This highlights the loss of the normal myoepithelial or stromal architecture found in healthy breast tissue.
##### Molecular Function
The upregulated functions are centered on binding and replication. Key activities include antigen binding (immune response), structural constituent of chromatin, and DNA helicase activity (DNA unzipping for replication). There is also activity related to ligand-gated ion channels and hormone activity.
The loss of function is tied to growth factor signaling and structural binding. Significant terms include actin binding, glycosaminoglycan binding, and growth factor activity. The downregulation of transmembrane receptor protein tyrosine kinase activity and Wnt-protein binding suggests that the normal signaling pathways that maintain healthy tissue have been shut down or bypassed.
##### Biological Processes
The dominant up-regulated processes are cell cycle progression and DNA replication. Terms like chromosome segregation, DNA-templated DNA replication, and nuclear chromosome segregation indicate that the HER2-enriched cells are rapidly dividing. Interestingly, there is also a clear immune signature involving B cell mediated immunity and immunoglobulin mediated immune response.
The suppressed processes are overwhelmingly related to muscle and organ development, such as muscle system process, muscle contraction, and heart morphogenesis. This reflects the displacement of normal, structured breast tissue by the invading tumor.

#### Gene set enrichment analysis-
##### Cellular Component
The tumor is highly enriched for nuclear components (nucleosomes, kinetochores, and condensed chromosomes) and immune complexes. The strong presence of the T cell receptor complex and immunoglobulin complex suggests significant immune infiltration within the HER2-enriched tumor.
he suppressed components are almost exclusively related to muscle architecture (sarcomere, Z disc, myofibril, and sarcolemma) and the extracellular matrix (collagen-containing matrix). This reinforces the idea that the structural integrity of normal breast tissue is lost in the cancerous state.
##### Molecular Function
The primary functions involve DNA/Chromatin handling (structural constituent of chromatin, helicase activity) and immune recognition (antigen binding). These functions provide the fuel for rapid growth and facilitate the interaction between the tumor and the immune system.
here is a notable suppression of signaling and transport functions, specifically fibroblast growth factor (FGF) receptor binding, peptide hormone binding, and G protein-coupled receptor (GPCR) activity. This suggests the tumor has become less responsive to normal physiological signaling pathways.
##### Biological Processes
There is a massive enrichment for mitosis and the cell cycle. Key activated pathways include protein localization to the centromeric region, regulation of mitotic sister chromatid separation, and nucleosome assembly. This indicates a high rate of DNA replication and cell division, characteristic of aggressive cancer.
Normal biological functions are "turned off." These include muscle system processes (like contraction), fat cell differentiation, and regulation of endothelial cell proliferation (vessel maintenance). This suggests that the tumor has displaced the normal fat and myoepithelial cells of the breast.

The GSEA results show that HER2-enriched breast cancer functions like a high-speed engine that has discarded its "normal" structural parts (muscle and fat markers) to focus entirely on replicating DNA and managing a complex immune microenvironment. The suppression of growth factor and hormone binding suggests that HER2-enriched tumors may rely more on their internal oncogenic signaling (like the HER2/PI3K/AKT pathway) rather than normal external growth signals.

### HER2 vs Basal
Out of the 55,095 genes analysed, 52% were differentially expressed with 15% upregulated in HER2-enriched breast cancer and 37% upregulated in BLBC.This suggests that the Basal-like subtype has a much more "active" or distinct transcriptomic program than HER2-enriched.
Some of the genes upregulated in HER2-enriched BC include PNMT, TBX10, TRIM3, XBP1 etc. Some of the genes upregulated in BLBC include SOX10, RGMA, TNMD etc.
#### Over-representation analysis-
HER2-enriched tumors are defined by active metabolic and signaling pathways, while BLBC is significantly more aggressive in its cell division and maintains a strong basal epithelial identity.
##### Cellular Component
HER2 cells are enriched for specialized structural and metabolic organelles. This includes motile cilia, axonemes, and peroxisomes/microbodies. The apical plasma membrane enrichment suggests a higher degree of cellular polarity compared to Basal cells.
BLBC is defined by its mitotic machinery. The enrichments are almost entirely nuclear: chromosomal regions, condensed chromosomes, and the MCM complex (essential for DNA replication). This reinforces the finding that BLBC is the more proliferative of the two aggressive subtypes.
##### Molecular Function
HER2 tumors are characterized by transmembrane transport and enzymatic activity. Specifically, metal ion transmembrane transporter activity and various oxidoreductase activities are prominent. This points to a cell that is actively moving ions and performing complex redox chemistry.
BLBC functions are centered on structural integrity and DNA replication. Key terms include DNA helicase activity (unzipping DNA), microtubule/tubulin binding (for spindle formation during mitosis), and structural constituent of skin epidermis.
##### Biological Processes
HER2-enriched tumors show a major focus on metabolic processing. High-level terms include hormone, steroid, and fatty acid metabolic processes, as well as xenobiotic metabolism. This suggests HER2 tumors are more biochemically active in processing lipids and hormones than their Basal counterparts.
The terms associated with BLBC are heavily focused on rapid cell division. There are mitotic nuclear division, regulation of nuclear division, and positive regulation of cell cycle. Additionally, there is a strong signal for skin/epidermis development, which is the classic transcriptomic signature of the Basal subtype.

#### Gene set enrichment analysis-
##### Cellular Component
There is a very specific activation of peroxisomes and microbodies. These organelles are critical for the breakdown of very-long-chain fatty acids, which aligns perfectly with the lipid metabolism.
BLBC shows higher expression of the perineuronal net, GABA receptor complexes, and keratin filaments. The presence of "synapse-associated" terms in BLBC is an interesting finding that suggests some neuro-epithelial signaling mimicking might be occurring in this subtype.
##### Molecular Function
High activity in FAD/FAD binding and oxidoreductase activity indicates a high redox metabolic rate. There is also a strong signal for antigen binding and transmembrane transporter activity (ATP-dependent), showing that HER2 cells are actively moving molecules across their membranes.
The most prominent function is the structural constituent of skin epidermis. This is the molecular "fingerprint" of the Basal-like subtype, reflecting the high concentration of keratins and intermediate filaments that give these cells their distinct structural properties.
##### Biological Processes
The metabolic machinery is significantly higher in HER2 tumors, particularly for monocarboxylic acid catabolism, cholesterol biosynthesis, and estrogen/steroid metabolism. This suggests HER2 tumors are more reliant on complex lipid and hormone processing for energy and signaling.
The BLBC side dominates terms related to keratinization and epidermis development.

#### Immune deconvolution-
Cancer Associated Fibroblasts (CAFs) are markedly and significantly higher in HER2 compared to Basal. This is consistent with the "mesenchymal" and stromal-heavy nature of HER2-enriched tumors. CAFs are responsible for synthesizing the extracellular matrix (ECM) and driving desmoplasia (tissue hardening), which is a prominent feature in HER2+ cancers that can sometimes act as a physical barrier to immune cell entry.
CD8+ T cells and Macrophages are both significantly higher in the Basal subtype. Literature consistently identifies Basal-like/Triple-Negative Breast Cancer (TNBC) as the most immunogenic and "inflamed" subtype. Higher CD8+ T cells infiltration in Basal-like tumors is a well-documented prognostic marker associated with better response to chemotherapy.
B Cells & CD4+ T Cells show no significant difference (ns) which is consistent with both being "immune-hot" compared to Luminal A/B; they share a similar baseline of general immune infiltration.
Uncharacterized Cells are higher in BLBC than HER2 enriched. This massive difference reflects the distinct cellularity of Basal-like cancers.
From this result, I observed that HER2 tumors have more fibroblasts and stroma, whereas Basal tumors exist in a state of high inflammatory signaling with a more direct infiltration of cytotoxic T cells and macrophages.
