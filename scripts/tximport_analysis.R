library(tximport)
library(DESeq2)

# ============================================================
# RNA-seq expression analysis
# Dataset: GSE151513
# Samples: six 6-hour Calu-3 samples
# Quantification: Salmon 2.8.0
# Reference: GENCODE v49
# ============================================================

# ------------------------------------------------------------
# 1. Project directories
# ------------------------------------------------------------

project_dir <- getwd()

salmon_dir <- file.path(
    project_dir,
    "results",
    "salmon"
)

output_dir <- file.path(
    project_dir,
    "results",
    "expression"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Sample metadata
# ------------------------------------------------------------

metadata <- data.frame(
    sample_id = c(
        "sample1",
        "sample2",
        "sample3",
        "sample4",
        "sample5",
        "sample6"
    ),

    sra_accession = c(
        "SRR11884716",
        "SRR11884717",
        "SRR11884718",
        "SRR11884719",
        "SRR11884720",
        "SRR11884721"
    ),

    condition = factor(
        c(
            "infected",
            "infected",
            "infected",
            "control",
            "control",
            "control"
        ),
        levels = c(
            "control",
            "infected"
        )
    ),

    timepoint = c(
        "6h",
        "6h",
        "6h",
        "6h",
        "6h",
        "6h"
    )
)

rownames(metadata) <- metadata$sample_id

# ------------------------------------------------------------
# 3. Locate Salmon quantification files
# ------------------------------------------------------------

quant_files <- file.path(
    salmon_dir,
    metadata$sra_accession,
    "quant.sf"
)

names(quant_files) <- metadata$sample_id

cat("\nChecking Salmon files...\n\n")

for (file in quant_files) {

    if (!file.exists(file)) {
        stop(
            "Missing Salmon quantification file: ",
            file
        )
    }

    cat("FOUND:", file, "\n")
}

# ------------------------------------------------------------
# 4. Read transcript names from Salmon
# ------------------------------------------------------------

cat("\nReading transcript identifiers...\n")

first_quant <- read.delim(
    quant_files[1],
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

transcript_names <- first_quant$Name

# ------------------------------------------------------------
# 5. Build transcript-to-gene mapping
# ------------------------------------------------------------

cat("\nBuilding transcript-to-gene mapping...\n")

tx_fields <- strsplit(
    transcript_names,
    "\\|"
)

gene_ids <- vapply(
    tx_fields,
    function(x) {
        if (length(x) >= 2) {
            return(x[2])
        } else {
            return(NA_character_)
        }
    },
    character(1)
)

tx2gene <- data.frame(
    transcript = transcript_names,
    gene = gene_ids,
    stringsAsFactors = FALSE
)

# Remove any invalid mappings

tx2gene <- tx2gene[
    !is.na(tx2gene$gene) &
    tx2gene$gene != "",
]

# ------------------------------------------------------------
# 6. Validate transcript-to-gene mapping
# ------------------------------------------------------------

if (nrow(tx2gene) == 0) {
    stop(
        "No valid transcript-to-gene mappings were created"
    )
}

cat(
    "\nTranscript records:",
    nrow(tx2gene),
    "\n"
)

cat(
    "Unique genes:",
    length(unique(tx2gene$gene)),
    "\n"
)

cat(
    "Example mappings:\n"
)

print(
    head(
        tx2gene,
        5
    )
)

# Save mapping

write.csv(
    tx2gene,
    file.path(
        output_dir,
        "tx2gene.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 7. Import Salmon quantifications
# ------------------------------------------------------------

cat(
    "\nImporting Salmon quantifications with tximport...\n"
)

txi <- tximport(
    quant_files,
    type = "salmon",
    tx2gene = tx2gene,
    countsFromAbundance = "no"
)

# ------------------------------------------------------------
# 8. Rename samples
# ------------------------------------------------------------

colnames(txi$counts) <- metadata$sample_id
colnames(txi$abundance) <- metadata$sample_id
colnames(txi$length) <- metadata$sample_id

# ------------------------------------------------------------
# 9. Report dimensions
# ------------------------------------------------------------

cat(
    "\nGene-level count matrix dimensions:\n"
)

cat(
    "Genes:",
    nrow(txi$counts),
    "\n"
)

cat(
    "Samples:",
    ncol(txi$counts),
    "\n"
)

# ------------------------------------------------------------
# 10. Save expression matrices
# ------------------------------------------------------------

write.csv(
    txi$counts,
    file.path(
        output_dir,
        "gene_counts.csv"
    )
)

write.csv(
    txi$abundance,
    file.path(
        output_dir,
        "gene_abundance.csv"
    )
)

write.csv(
    txi$length,
    file.path(
        output_dir,
        "gene_effective_lengths.csv"
    )
)

write.csv(
    metadata,
    file.path(
        output_dir,
        "sample_metadata.csv"
    ),
    row.names = TRUE
)

cat(
    "\nExpression matrices saved to:\n",
    output_dir,
    "\n"
)

# ------------------------------------------------------------
# 11. Create DESeq2 dataset
# ------------------------------------------------------------

cat(
    "\nCreating DESeq2 dataset...\n"
)

dds <- DESeqDataSetFromTximport(
    txi = txi,
    colData = metadata,
    design = ~ condition
)

# ------------------------------------------------------------
# 12. Filter low-count genes
# ------------------------------------------------------------

keep <- rowSums(
    counts(dds) >= 10
) >= 3

dds <- dds[
    keep,
]

cat(
    "\nGenes retained after filtering:",
    nrow(dds),
    "\n"
)

# ------------------------------------------------------------
# 13. Run DESeq2
# ------------------------------------------------------------

cat(
    "\nRunning DESeq2...\n"
)

dds <- DESeq(
    dds
)

# ------------------------------------------------------------
# 14. Infected vs control
# ------------------------------------------------------------

cat(
    "\nCalculating infected vs control results...\n"
)

results_de <- results(
    dds,
    contrast = c(
        "condition",
        "infected",
        "control"
    )
)

# Order by adjusted p-value

results_de <- results_de[
    order(
        results_de$padj,
        na.last = NA
    ),
]

results_table <- as.data.frame(
    results_de
)

results_table$gene_id <- rownames(
    results_table
)

results_table <- results_table[
    ,
    c(
        "gene_id",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj"
    )
]

# ------------------------------------------------------------
# 15. Save differential-expression results
# ------------------------------------------------------------

write.csv(
    results_table,
    file.path(
        output_dir,
        "DESeq2_infected_vs_control.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 16. Summarize significant genes
# ------------------------------------------------------------

significant <- (
    !is.na(results_table$padj) &
    results_table$padj < 0.05
)

upregulated <- (
    significant &
    results_table$log2FoldChange > 1
)

downregulated <- (
    significant &
    results_table$log2FoldChange < -1
)

cat(
    "\n========================================\n"
)

cat(
    "Differential-expression summary\n"
)

cat(
    "========================================\n"
)

cat(
    "Genes tested:",
    nrow(results_table),
    "\n"
)

cat(
    "Significant genes (padj < 0.05):",
    sum(significant),
    "\n"
)

cat(
    "Upregulated (padj < 0.05, log2FC > 1):",
    sum(upregulated),
    "\n"
)

cat(
    "Downregulated (padj < 0.05, log2FC < -1):",
    sum(downregulated),
    "\n"
)

cat(
    "\nAnalysis completed successfully\n"
)