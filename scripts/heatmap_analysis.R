library(DESeq2)
library(pheatmap)

# ============================================================
# Heatmap of significant genes
# Dataset: GSE151513
# Six 6-hour Calu-3 samples
# ============================================================

project_dir <- getwd()

expression_dir <- file.path(
    project_dir,
    "results",
    "expression"
)

output_dir <- file.path(
    project_dir,
    "results",
    "heatmap"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ------------------------------------------------------------
# Load gene-level counts
# ------------------------------------------------------------

counts_file <- file.path(
    expression_dir,
    "gene_counts.csv"
)

gene_counts <- read.csv(
    counts_file,
    row.names = 1,
    check.names = FALSE
)

# ------------------------------------------------------------
# Load sample metadata
# ------------------------------------------------------------

metadata_file <- file.path(
    expression_dir,
    "sample_metadata.csv"
)

metadata <- read.csv(
    metadata_file,
    row.names = 1,
    check.names = FALSE
)

metadata$condition <- factor(
    metadata$condition,
    levels = c("control", "infected")
)

# Make sure sample order matches
gene_counts <- gene_counts[, rownames(metadata)]

# ------------------------------------------------------------
# Load significant genes
# ------------------------------------------------------------

significant_file <- file.path(
    expression_dir,
    "DESeq2_significant_genes_with_symbols.csv"
)

significant_genes <- read.csv(
    significant_file,
    stringsAsFactors = FALSE
)

cat("Significant genes:", nrow(significant_genes), "\n")

# ------------------------------------------------------------
# Create DESeq2 object
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
    countData = round(gene_counts),
    colData = metadata,
    design = ~ condition
)

# ------------------------------------------------------------
# Filter low-count genes
# ------------------------------------------------------------

keep <- rowSums(counts(dds) >= 10) >= 3

dds <- dds[keep, ]

cat(
    "Genes retained after filtering:",
    nrow(dds),
    "\n"
)

# ------------------------------------------------------------
# Run DESeq2
# ------------------------------------------------------------

dds <- DESeq(dds)

# ------------------------------------------------------------
# Variance-stabilizing transformation
# ------------------------------------------------------------

vsd <- vst(
    dds,
    blind = FALSE
)

vsd_matrix <- assay(vsd)

# ------------------------------------------------------------
# Select significant genes
# ------------------------------------------------------------

gene_ids <- significant_genes$gene_id

gene_ids <- gene_ids[
    gene_ids %in% rownames(vsd_matrix)
]

cat(
    "Significant genes available for heatmap:",
    length(gene_ids),
    "\n"
)

heatmap_matrix <- vsd_matrix[
    gene_ids,
    ,
    drop = FALSE
]

# ------------------------------------------------------------
# Replace gene IDs with gene symbols for display
# ------------------------------------------------------------

gene_symbols <- significant_genes$gene_symbol[
    match(
        gene_ids,
        significant_genes$gene_id
    )
]

display_names <- paste(
    gene_symbols,
    gene_ids,
    sep = " | "
)

rownames(heatmap_matrix) <- display_names

# ------------------------------------------------------------
# Z-score each gene
# ------------------------------------------------------------

heatmap_matrix_scaled <- t(
    scale(
        t(heatmap_matrix)
    )
)

# ------------------------------------------------------------
# Sample annotation
# ------------------------------------------------------------

annotation_col <- data.frame(
    Condition = metadata$condition
)

rownames(annotation_col) <- rownames(metadata)

# ------------------------------------------------------------
# Save heatmap matrix
# ------------------------------------------------------------

write.csv(
    heatmap_matrix_scaled,
    file.path(
        output_dir,
        "significant_genes_zscore.csv"
    )
)

# ------------------------------------------------------------
# Create heatmap
# ------------------------------------------------------------

heatmap_file <- file.path(
    output_dir,
    "significant_genes_heatmap.png"
)

pheatmap(
    heatmap_matrix_scaled,
    annotation_col = annotation_col,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "none",
    fontsize_row = 9,
    fontsize_col = 10,
    border_color = NA,
    main = "Significant genes: infected vs control",
    filename = heatmap_file,
    width = 9,
    height = 7
)

cat(
    "\nHeatmap saved to:",
    heatmap_file,
    "\n"
)

cat(
    "\nHeatmap analysis completed successfully\n"
)