library(ggplot2)

project_dir <- getwd()

input_file <- file.path(
    project_dir,
    "results",
    "pathway_analysis",
    "hallmark_gsea_significant.csv"
)

output_dir <- file.path(
    project_dir,
    "results",
    "pathway_analysis"
)

if (!file.exists(input_file)) {
    stop(
        "GSEA significant results not found: ",
        input_file
    )
}

gsea <- read.csv(
    input_file,
    stringsAsFactors = FALSE
)

if (nrow(gsea) == 0) {
    stop(
        "No significant pathways were found."
    )
}

gsea <- gsea[
    order(gsea$NES),
]

gsea$pathway <- gsub(
    "^HALLMARK_",
    "",
    gsea$pathway
)

gsea$pathway <- gsub(
    "_",
    " ",
    gsea$pathway
)

gsea$pathway <- tools::toTitleCase(
    tolower(gsea$pathway)
)

gsea$pathway <- factor(
    gsea$pathway,
    levels = gsea$pathway
)

gsea$direction <- ifelse(
    gsea$NES > 0,
    "Enriched in infected",
    "Enriched in control"
)

plot <- ggplot(
    gsea,
    aes(
        x = NES,
        y = pathway
    )
) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed"
    ) +
    geom_point(
        aes(
            size = -log10(padj)
        )
    ) +
    labs(
        title = "Hallmark pathway enrichment",
        subtitle = "Preranked GSEA of DESeq2 Wald statistics",
        x = "Normalized enrichment score (NES)",
        y = NULL,
        size = "-log10 adjusted p-value"
    ) +
    theme_minimal(
        base_size = 12
    ) +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(
            size = 9
        )
    )

ggsave(
    filename = file.path(
        output_dir,
        "hallmark_gsea_plot.png"
    ),
    plot = plot,
    width = 10,
    height = 8,
    dpi = 300
)

cat(
    "GSEA plot saved to:",
    file.path(
        output_dir,
        "hallmark_gsea_plot.png"
    ),
    "\n"
)