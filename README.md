# RNA-seq Differential Expression Analysis of SARS-CoV-2-Infected Calu-3 Cells

A reproducible RNA-seq analysis pipeline using publicly available NCBI SRA data to investigate transcriptional differences between SARS-CoV-2-infected and control human Calu-3 cells at 6 hours post-infection.

## Project overview

This project analyzes paired-end RNA-seq data from the publicly available GEO dataset **GSE151513**.

The analysis focuses on six samples collected at **6 hours post-infection (6 hpi)**:

* 3 SARS-CoV-2-infected samples
* 3 control samples

The objective was to build an end-to-end RNA-seq workflow starting from raw sequencing data and progressing through quality assessment, transcript quantification, gene-level expression analysis, dimensionality reduction, differential expression, and visualization.

The project was developed as a practical bioinformatics portfolio project using Python and R.

## Biological question

**How does SARS-CoV-2 infection affect gene expression in human Calu-3 cells at 6 hours post-infection?**

The primary comparison is:

**Infected vs control**

with three biological replicates per condition.

## Dataset

Source: NCBI Gene Expression Omnibus / Sequence Read Archive

* GEO accession: **GSE151513**
* Organism: *Homo sapiens*
* Cell model: Calu-3
* Experimental condition: SARS-CoV-2 infection
* Time point analyzed: 6 hours post-infection
* Sequencing: paired-end RNA-seq

### Samples analyzed

| Sample  | SRA accession | GEO accession | Condition | Time point |
| ------- | ------------- | ------------- | --------- | ---------- |
| sample1 | SRR11884716   | GSM4579971    | Infected  | 6 h        |
| sample2 | SRR11884717   | GSM4579972    | Infected  | 6 h        |
| sample3 | SRR11884718   | GSM4579973    | Infected  | 6 h        |
| sample4 | SRR11884719   | GSM4579974    | Control   | 6 h        |
| sample5 | SRR11884720   | GSM4579975    | Control   | 6 h        |
| sample6 | SRR11884721   | GSM4579976    | Control   | 6 h        |

## Analysis workflow

```text
NCBI SRA
   ↓
SRA Toolkit
   ↓
Paired-end FASTQ files
   ↓
Custom Python quality control
   ↓
GENCODE v49 transcript reference
   ↓
Salmon transcript quantification
   ↓
tximport
   ↓
Gene-level expression matrix
   ↓
DESeq2
   ↓
 ┌───────────────┬───────────────┬────────────────┐
 │ PCA           │ Differential  │ Visualization  │
 │               │ expression    │                │
 │               │               │                │
 │               │ Volcano plot  │ Heatmap        │
 └───────────────┴───────────────┴────────────────┘
```

## Software and tools

### Data acquisition

* NCBI SRA Toolkit 3.4.1
* `prefetch`
* `fasterq-dump`

### Python

* Python 3
* pandas
* NumPy
* SciPy
* Matplotlib

Python was used for:

* FASTQ quality assessment
* Per-position Phred quality calculation
* GC-content calculation
* QC report generation
* QC visualization
* Gene-symbol mapping

### R

* R 4.4.2
* DESeq2
* tximport
* ggplot2
* pheatmap

R was used for:

* Transcript-to-gene expression import
* Gene-level expression analysis
* Variance-stabilizing transformation
* PCA
* Differential expression analysis
* Heatmap generation
* Volcano plot generation

### Transcript quantification

* Salmon 2.8.0

Salmon was used for transcript-level quantification from paired-end FASTQ files.

## Quality control

A custom Python FASTQ parser was developed to calculate:

* Number of reads
* Total sequenced bases
* Average read length
* Mean Phred quality
* Percentage of bases below Q20
* GC content
* Average Phred quality at each read position

The six samples produced 12 paired-end FASTQ files in total.

The quality-control workflow was designed to assess sequencing quality before downstream transcript quantification.

### Example QC output

Quality-by-position plots are available in:

`results/qc/`

The corresponding numerical summaries are stored as CSV files.

### QC limitation

The custom QC analysis provides basic sequencing-quality metrics but does not replace comprehensive tools such as FastQC/MultiQC.

In particular, these metrics alone cannot establish the absence of adapter contamination, sequence duplication, rRNA contamination, or other library-specific technical biases.

## Transcript quantification

Transcript abundance was estimated using **Salmon 2.8.0** with the GENCODE v49 human transcript reference.

The resulting `quant.sf` files were imported into R using **tximport** and summarized to gene-level expression.

The GENCODE reference used in this project is based on **GRCh38** and differs from the reference annotation used in the original study. This reference choice is documented explicitly to improve reproducibility.

Large reference files and Salmon output directories are excluded from Git because they are too large for practical repository storage.

## Differential expression analysis

Gene-level expression was analyzed using **DESeq2**.

Low-count genes were filtered using:

* Minimum count of 10
* Present in at least 3 samples

After filtering:

* **18,619 genes** remained for the DESeq2 dataset
* **18,609 genes** were included in the final differential-expression results

The primary statistical comparison was:

**Infected vs control**

### Differential expression results

Using an adjusted p-value threshold of:

**padj < 0.05**

24 genes showed statistically significant differential expression.

Using the additional effect-size criterion:

**|log2 fold change| > 1**

9 genes met both criteria:

* 6 upregulated
* 3 downregulated

The complete DESeq2 results are available in:

`results/expression/DESeq2_infected_vs_control.csv`

A table containing gene symbols for the significant genes is available in:

`results/expression/DESeq2_significant_genes_with_symbols.csv`

## PCA

Principal component analysis was performed using a variance-stabilizing transformation of the gene-count data.

The first two principal components explained:

* **PC1: 28.1%**
* **PC2: 22.6%**
* **Combined: 50.7%**

The PCA provides an assessment of overall transcriptional variation among the six samples.

The samples do not show complete separation of infected and control groups across the whole transcriptome. This indicates that biological or technical variation remains substantial in this small dataset.

PCA results:

`results/pca/PCA_coordinates.csv`

PCA visualization:

`results/pca/PCA_plot.png`

## Differential-expression visualization

### Volcano plot

The volcano plot displays the relationship between log2 fold change and adjusted statistical significance.

![Volcano plot](results/volcano/volcano_plot.png)

The plot highlights genes meeting the predefined thresholds of:

* adjusted p-value < 0.05
* log2 fold change > 1 or < -1

### Heatmap

A heatmap was generated for the 9 genes satisfying both the adjusted p-value and effect-size thresholds.

![Significant gene heatmap](results/heatmap/significant_genes_heatmap.png)

The selected genes show a clearer condition-associated expression pattern than the whole-transcriptome PCA.

This difference is expected because the heatmap specifically focuses on genes identified as strongly associated with the infected-versus-control comparison.

## Significant genes

The nine genes meeting both statistical and effect-size thresholds were:

| Gene symbol | Ensembl gene ID    | log2 fold change | Direction |
| ----------- | ------------------ | ---------------: | --------- |
| ARHGEF35    | ENSG00000288422.1  |             8.45 | Up        |
| AATF        | ENSG00000276072.2  |             8.02 | Up        |
| APOM        | ENSG00000235754.7  |            -5.13 | Down      |
| MUC20-OT1   | ENSG00000281915.2  |            -7.59 | Down      |
| HLA-DMB     | ENSG00000241296.8  |             7.33 | Up        |
| HSPA1L      | ENSG00000226704.6  |             7.13 | Up        |
| PLEKHG5     | ENSG00000171680.24 |             1.96 | Up        |
| APOM        | ENSG00000204444.11 |             3.22 | Up        |
| H2AC19      | ENSG00000288859.2  |            -1.89 | Down      |

Gene symbols are provided for biological readability, while Ensembl gene IDs are retained as the primary identifiers.

Two different Ensembl gene identifiers map to **APOM** in the selected annotation. They were intentionally retained as separate identifiers rather than manually merged.

## Repository structure

```text
rna-seq-analysis/
│
├── data/
│   └── metadata.csv
│
├── reference/
│   └── Salmon reference files
│
├── results/
│   ├── expression/
│   │   ├── DESeq2_infected_vs_control.csv
│   │   ├── DESeq2_significant_genes_with_symbols.csv
│   │   └── sample_metadata.csv
│   │
│   ├── heatmap/
│   │   ├── significant_genes_heatmap.png
│   │   └── significant_genes_zscore.csv
│   │
│   ├── pca/
│   │   ├── PCA_coordinates.csv
│   │   └── PCA_plot.png
│   │
│   ├── qc/
│   │   ├── QC summary files
│   │   └── quality-by-position plots
│   │
│   └── volcano/
│       └── volcano_plot.png
│
└── scripts/
    ├── run_qc.py
    ├── plot_qc_from_csv.py
    ├── run_salmon.py
    ├── tximport_analysis.R
    ├── pca_analysis.R
    ├── volcano_plot.R
    ├── heatmap_analysis.R
    └── map_gene_symbols.py
```

## Reproducibility

The project separates large raw/reference data from analysis code and lightweight results.

Raw sequencing files, reference files, Salmon indexes, and Salmon quantification directories are excluded through `.gitignore`.

The analysis scripts are retained in the repository so that the workflow can be inspected and reproduced.

The six SRA accessions required to reconstruct the dataset are documented in:

`data/metadata.csv`

## Limitations

Several limitations should be considered when interpreting these results:

1. Only six samples were analyzed, with three biological replicates per condition.
2. The dataset represents a single time point: 6 hours post-infection.
3. The custom QC workflow provides basic sequencing metrics rather than a complete MultiQC/FastQC assessment.
4. No adapter trimming was performed in the current workflow.
5. The PCA shows substantial variation among samples and does not produce perfect separation of the two conditions.
6. Only a small number of genes meet both the statistical-significance and effect-size thresholds.
7. The analysis uses a modern GENCODE v49 reference rather than exactly reproducing the original study's reference annotation.
8. Differential expression should therefore be interpreted as an analysis of this selected subset of samples rather than a definitive characterization of the complete SARS-CoV-2 response.

## Future improvements

Potential extensions include:

* FastQC/MultiQC-based sequencing QC
* Adapter detection and trimming
* Alignment-based quantification for comparison with Salmon
* Gene ontology and pathway enrichment
* Ranked gene-set enrichment analysis
* Integration of additional time points
* Inclusion of additional biological replicates
* Comparison against the original study's reference and analysis pipeline
* More extensive batch and technical-variable assessment

## Project objective

The primary purpose of this project was to develop practical experience with an end-to-end RNA-seq workflow, including:

**raw sequencing data → quality control → transcript quantification → gene-level expression → statistical testing → visualization → biological interpretation**

The project emphasizes reproducibility, explicit documentation of analytical decisions, and cautious interpretation of results rather than treating statistical significance alone as evidence of biological importance.
## Ranked Pathway Enrichment Analysis

### Objective

To investigate biological processes associated with SARS-CoV-2 infection, I performed preranked Gene Set Enrichment Analysis (GSEA) using the complete set of gene-level differential-expression statistics rather than restricting the analysis to individually significant genes.

### Method

* **Input:** DESeq2 differential-expression results for 18,609 genes.
* **Ranking metric:** DESeq2 Wald test statistic, ranked from highest to lowest.
* **Gene sets:** MSigDB Hallmark collection for *Homo sapiens*, accessed using `msigdbr`.
* **Enrichment analysis:** `fgsea`, with a minimum gene-set size of 15 and a maximum of 500.
* **Multiple testing:** Benjamini–Hochberg-adjusted p-values.
* **Significance threshold:** Adjusted p-value < 0.05.
* **Visualization:** Normalized enrichment scores (NES) plotted for significant pathways.

The analysis used the Wald statistic to retain information about both the direction and strength of differential expression across the ranked gene list. Positive NES values indicate enrichment toward the infected condition; negative NES values indicate enrichment toward the control condition.

### Results

Eighteen Hallmark pathways met the adjusted p-value threshold of 0.05.

The strongest enrichment patterns included:

| Pathway                         |    NES | Adjusted p-value | Direction |
| ------------------------------- | -----: | ---------------: | --------- |
| Oxidative phosphorylation       | -2.666 |     7.52 × 10⁻¹⁸ | Control   |
| MYC targets V1                  | -2.522 |     8.87 × 10⁻¹⁵ | Control   |
| TNFα signaling via NF-κB        |  1.541 |      8.66 × 10⁻³ | Infected  |
| Hedgehog signaling              |  1.727 |      2.24 × 10⁻² | Infected  |
| Glycolysis                      | -1.553 |      5.81 × 10⁻³ | Control   |
| Fatty-acid metabolism           | -1.620 |      7.59 × 10⁻³ | Control   |
| Unfolded protein response       | -1.450 |      3.58 × 10⁻² | Control   |
| Reactive oxygen species pathway | -1.604 |      4.01 × 10⁻² | Control   |

*Negative NES indicates enrichment toward the control side of the ranked gene list, not necessarily that every gene in the pathway is downregulated.*

### Interpretation

The results suggest several patterns worth investigating further:

1. **Inflammatory signaling:** TNFα signaling via NF-κB was enriched toward the infected condition, consistent with a potential difference in inflammatory-response-associated transcription.
2. **Cellular energy metabolism:** Oxidative phosphorylation, glycolysis, and fatty-acid metabolism were enriched toward the control condition, suggesting differences in metabolic gene-expression programs.
3. **Cellular growth and stress responses:** MYC targets, the unfolded protein response, and the reactive oxygen species pathway were also enriched toward the control condition.
4. **Additional signaling differences:** Hedgehog signaling was enriched toward the infected condition.

These findings are exploratory associations from this dataset. They do not establish that infection directly caused the observed pathway patterns or identify the underlying mechanisms.

### Reproducibility and outputs

The analysis scripts and selected outputs are available in this repository:

* `scripts/ranked_pathway_analysis.R` — performs preranked Hallmark GSEA.
* `scripts/plot_gsea.R` — generates the pathway enrichment plot.
* `scripts/create_gsea_interpretation.py` — generates a readable pathway interpretation table.
* `results/pathway_analysis/hallmark_gsea_plot.png` — visualization of significant pathways.
* `results/pathway_analysis/hallmark_gsea_significant.csv` — significant pathway results.
* `results/pathway_analysis/hallmark_gsea_interpretation.csv` — pathway directions, statistics, and descriptive interpretations.

### Limitations

The comparison includes three infected and three control samples at six hours post-infection. The small sample size limits statistical power and generalizability. Principal component analysis did not show clear separation of all samples by condition, so the pathway findings should be interpreted cautiously. Independent validation and additional biological replicates would strengthen the conclusions.

The Hallmark gene sets are broad biological signatures. Enrichment indicates that genes from a set are disproportionately represented toward one end of the ranked list; it does not mean every gene in that set changes in the same direction or that the entire pathway is activated or inhibited.
