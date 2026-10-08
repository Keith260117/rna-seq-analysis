library(DESeq2)
library(ggplot2)

# ============================================================
# PCA analysis of RNA-seq samples
# Dataset: GSE151513
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
    "pca"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ------------------------------------------------------------
# Recreate DESeq2 dataset from saved expression data
# ------------------------------------------------------------

gene_counts <- read.csv(
    file.path(
        expression_dir,
        "gene_counts.csv"
    ),
    row.names = 1,
    check.names = FALSE
)

metadata <- read.csv(
    file.path(
        expression_dir,
        "sample_metadata.csv"
    ),
    row.names = 1,
    check.names = FALSE
)

metadata$condition <- factor(
    metadata$condition,
    levels = c(
        "control",
        "infected"
    )
)

# Make sure count columns match metadata

gene_counts <- gene_counts[
    ,
    rownames(metadata)
]

# ------------------------------------------------------------
# Build DESeq2 object
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
    countData = round(gene_counts),
    colData = metadata,
    design = ~ condition
)

# ------------------------------------------------------------
# Filter low-count genes
# ------------------------------------------------------------

keep <- rowSums(
    counts(dds) >= 10
) >= 3

dds <- dds[
    keep,
]

cat(
    "Genes retained:",
    nrow(dds),
    "\n"
)

# ------------------------------------------------------------
# Run DESeq2 normalization
# ------------------------------------------------------------

dds <- DESeq(
    dds
)

# ------------------------------------------------------------
# Variance-stabilizing transformation
# ------------------------------------------------------------

vsd <- vst(
    dds,
    blind = FALSE
)

# ------------------------------------------------------------
# Calculate PCA
# ------------------------------------------------------------

pca_data <- plotPCA(
    vsd,
    intgroup = "condition",
    returnData = TRUE
)

percent_variance <- round(
    100 * attr(
        pca_data,
        "percentVar"
    ),
    1
)

# Add sample identifiers

pca_data$sample_id <- rownames(
    pca_data
)

# ------------------------------------------------------------
# Print PCA information
# ------------------------------------------------------------

cat(
    "\nPCA variance explained:\n"
)

cat(
    "PC1:",
    percent_variance[1],
    "%\n"
)

cat(
    "PC2:",
    percent_variance[2],
    "%\n"
)

cat(
    "\nPCA sample coordinates:\n"
)

print(
    pca_data[
        ,
        c(
            "sample_id",
            "condition",
            "PC1",
            "PC2"
        )
    ]
)

# ------------------------------------------------------------
# Save PCA coordinates
# ------------------------------------------------------------

write.csv(
    pca_data,
    file.path(
        output_dir,
        "PCA_coordinates.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Create PCA plot
# ------------------------------------------------------------

pca_plot <- ggplot(
    pca_data,
    aes(
        x = PC1,
        y = PC2,
        shape = condition,
        label = sample_id
    )
) +
    geom_point(
        size = 4
    ) +
    geom_text(
        vjust = -1,
        size = 3
    ) +
    xlab(
        paste0(
            "PC1 (",
            percent_variance[1],
            "%)"
        )
    ) +
    ylab(
        paste0(
            "PC2 (",
            percent_variance[2],
            "%)"
        )
    ) +
    ggtitle(
        "PCA of RNA-seq samples"
    ) +
    theme_minimal(
        base_size = 12
    )

# ------------------------------------------------------------
# Save plot
# ------------------------------------------------------------

ggsave(
    filename = file.path(
        output_dir,
        "PCA_plot.png"
    ),
    plot = pca_plot,
    width = 8,
    height = 6,
    dpi = 300
)

cat(
    "\nPCA plot saved to:",
    file.path(
        output_dir,
        "PCA_plot.png"
    ),
    "\n"
)

cat(
    "\nPCA analysis completed successfully\n"
)