library(ggplot2)

project_dir <- getwd()

expression_dir <- file.path(
    project_dir,
    "results",
    "expression"
)

output_dir <- file.path(
    project_dir,
    "results",
    "volcano"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

results_file <- file.path(
    expression_dir,
    "DESeq2_infected_vs_control.csv"
)

results_table <- read.csv(
    results_file,
    stringsAsFactors = FALSE
)

cat("Loaded DESeq2 results\n")
cat("Genes in results table:", nrow(results_table), "\n")

# ------------------------------------------------------------
# Classify genes
# ------------------------------------------------------------

results_table$category <- "Not significant"

results_table$category[
    !is.na(results_table$padj) &
    results_table$padj < 0.05 &
    results_table$log2FoldChange > 1
] <- "Upregulated"

results_table$category[
    !is.na(results_table$padj) &
    results_table$padj < 0.05 &
    results_table$log2FoldChange < -1
] <- "Downregulated"

# ------------------------------------------------------------
# Calculate -log10 adjusted p-value
# ------------------------------------------------------------

results_table$neg_log10_padj <- NA_real_

valid_padj <- (
    !is.na(results_table$padj) &
    results_table$padj > 0
)

results_table$neg_log10_padj[valid_padj] <- (
    -log10(results_table$padj[valid_padj])
)

# ------------------------------------------------------------
# Print summary
# ------------------------------------------------------------

cat("\nVolcano plot categories\n")
cat("----------------------------------------\n")

cat(
    "Not significant:",
    sum(results_table$category == "Not significant"),
    "\n"
)

cat(
    "Upregulated:",
    sum(results_table$category == "Upregulated"),
    "\n"
)

cat(
    "Downregulated:",
    sum(results_table$category == "Downregulated"),
    "\n"
)

# ------------------------------------------------------------
# Save classified results
# ------------------------------------------------------------

write.csv(
    results_table,
    file.path(
        output_dir,
        "DESeq2_volcano_data.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Create volcano plot
# ------------------------------------------------------------

volcano_plot <- ggplot(
    results_table,
    aes(
        x = log2FoldChange,
        y = neg_log10_padj,
        color = category
    )
) +

    geom_point(
        alpha = 0.7,
        size = 2
    ) +

    geom_vline(
        xintercept = c(-1, 1),
        linetype = "dashed"
    ) +

    geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed"
    ) +

    xlab("Log2 fold change") +

    ylab("-Log10 adjusted p-value") +

    ggtitle(
        "Differential expression: infected vs control"
    ) +

    theme_minimal(
        base_size = 12
    )

ggsave(
    filename = file.path(
        output_dir,
        "volcano_plot.png"
    ),
    plot = volcano_plot,
    width = 8,
    height = 6,
    dpi = 300
)

cat(
    "\nVolcano plot saved to:",
    file.path(
        output_dir,
        "volcano_plot.png"
    ),
    "\n"
)

cat("\nVolcano analysis completed successfully\n")