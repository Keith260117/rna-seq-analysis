library(fgsea)
library(msigdbr)

project_dir <- getwd()

input_file <- file.path(
    project_dir,
    "results",
    "expression",
    "DESeq2_infected_vs_control.csv"
)

output_dir <- file.path(
    project_dir,
    "results",
    "pathway_analysis"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!file.exists(input_file)) {
    stop(
        "DESeq2 results file not found: ",
        input_file
    )
}

cat("Reading DESeq2 results...\n")

de_results <- read.csv(
    input_file,
    stringsAsFactors = FALSE
)

required_columns <- c(
    "gene_id",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj"
)

missing_columns <- setdiff(
    required_columns,
    colnames(de_results)
)

if (length(missing_columns) > 0) {
    stop(
        "Missing required columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

cat(
    "Total DESeq2 rows:",
    nrow(de_results),
    "\n"
)

de_results <- de_results[
    !is.na(de_results$stat) &
    !is.na(de_results$gene_id),
]

cat(
    "Rows with valid Wald statistics:",
    nrow(de_results),
    "\n"
)

de_results$gene_id <- sub(
    "\\..*$",
    "",
    de_results$gene_id
)

de_results <- de_results[
    !duplicated(de_results$gene_id),
]

ranking <- de_results$stat

names(ranking) <- de_results$gene_id

ranking <- sort(
    ranking,
    decreasing = TRUE
)

cat(
    "Genes in ranked list:",
    length(ranking),
    "\n"
)

cat(
    "Top positively ranked genes:\n"
)

print(
    head(ranking, 10)
)

cat(
    "\nTop negatively ranked genes:\n"
)

print(
    head(
        sort(ranking),
        10
    )
)

cat("\nLoading MSigDB Hallmark gene sets...\n")

hallmark <- msigdbr(
    species = "Homo sapiens",
    collection = "H"
)

cat(
    "Hallmark gene-set records:",
    nrow(hallmark),
    "\n"
)

hallmark_sets <- split(
    hallmark$ensembl_gene,
    hallmark$gs_name
)

hallmark_sets <- lapply(
    hallmark_sets,
    function(x) {
        unique(
            x[
                !is.na(x) &
                x != ""
            ]
        )
    }
)

hallmark_sets <- hallmark_sets[
    lengths(hallmark_sets) > 0
]

cat(
    "Hallmark pathways:",
    length(hallmark_sets),
    "\n"
)

cat("\nRunning preranked GSEA...\n")

set.seed(123)

gsea_results <- fgsea(
    pathways = hallmark_sets,
    stats = ranking,
    minSize = 15,
    maxSize = 500,
    eps = 0
)

gsea_results <- as.data.frame(
    gsea_results
)

gsea_results <- gsea_results[
    order(
        gsea_results$padj,
        -abs(gsea_results$NES)
    ),
]

significant_results <- gsea_results[
    !is.na(gsea_results$padj) &
    gsea_results$padj < 0.05,
]

cat(
    "\nSignificant Hallmark pathways (padj < 0.05):",
    nrow(significant_results),
    "\n"
)

if (nrow(significant_results) > 0) {

    cat("\nSignificant pathways:\n")

    display_results <- data.frame(
        pathway = significant_results$pathway,
        NES = significant_results$NES,
        pval = significant_results$pval,
        padj = significant_results$padj,
        stringsAsFactors = FALSE
    )

    print(
        display_results,
        row.names = FALSE
    )

} else {

    cat(
        "\nNo Hallmark pathways reached padj < 0.05.\n"
    )
}

cat(
    "\nSaving pathway results...\n"
)

gsea_output <- data.frame(
    pathway = gsea_results$pathway,
    pval = gsea_results$pval,
    padj = gsea_results$padj,
    log2err = gsea_results$log2err,
    ES = gsea_results$ES,
    NES = gsea_results$NES,
    size = gsea_results$size,
    leadingEdge = sapply(
        gsea_results$leadingEdge,
        function(x) {
            paste(x, collapse = ";")
        }
    ),
    stringsAsFactors = FALSE
)

write.csv(
    gsea_output,
    file.path(
        output_dir,
        "hallmark_gsea_results.csv"
    ),
    row.names = FALSE
)

significant_output <- data.frame(
    pathway = significant_results$pathway,
    pval = significant_results$pval,
    padj = significant_results$padj,
    log2err = significant_results$log2err,
    ES = significant_results$ES,
    NES = significant_results$NES,
    size = significant_results$size,
    leadingEdge = sapply(
        significant_results$leadingEdge,
        function(x) {
            paste(x, collapse = ";")
        }
    ),
    stringsAsFactors = FALSE
)

write.csv(
    significant_output,
    file.path(
        output_dir,
        "hallmark_gsea_significant.csv"
    ),
    row.names = FALSE
)

cat(
    "\nSaving ranked gene list...\n"
)

ranked_gene_table <- data.frame(
    gene_id = names(ranking),
    wald_statistic = as.numeric(ranking)
)

write.csv(
    ranked_gene_table,
    file.path(
        output_dir,
        "ranked_gene_list.csv"
    ),
    row.names = FALSE
)

cat(
    "\nPathway analysis completed successfully.\n"
)

cat(
    "Results saved to:",
    output_dir,
    "\n"
)