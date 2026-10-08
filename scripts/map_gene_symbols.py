from pathlib import Path
import pandas as pd

# ------------------------------------------------------------
# File locations
# ------------------------------------------------------------

project_dir = Path.cwd()

reference_file = (
    project_dir
    / "reference"
    / "gencode.v49.transcripts.fa"
)

results_file = (
    project_dir
    / "results"
    / "volcano"
    / "DESeq2_volcano_data.csv"
)

output_dir = (
    project_dir
    / "results"
    / "expression"
)

output_file = (
    output_dir
    / "DESeq2_significant_genes_with_symbols.csv"
)

# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

if not reference_file.is_file():
    raise FileNotFoundError(
        f"Reference FASTA not found: {reference_file}"
    )

if not results_file.is_file():
    raise FileNotFoundError(
        f"DESeq2 results not found: {results_file}"
    )

output_dir.mkdir(
    parents=True,
    exist_ok=True
)

# ------------------------------------------------------------
# Load DESeq2 results
# ------------------------------------------------------------

results = pd.read_csv(results_file)

significant = results[
    results["padj"].notna()
    & (results["padj"] < 0.05)
    & (results["log2FoldChange"].abs() > 1)
].copy()

print("Significant genes found:", len(significant))

# ------------------------------------------------------------
# Remove Ensembl version suffix
#
# Example:
# ENSG00000288422.1
# becomes
# ENSG00000288422
# ------------------------------------------------------------

significant["gene_id_base"] = (
    significant["gene_id"]
    .str.split(".")
    .str[0]
)

# ------------------------------------------------------------
# Parse GENCODE transcript FASTA
# ------------------------------------------------------------

gene_mapping = {}

with open(reference_file, "r") as file:

    for line in file:

        if not line.startswith(">"):
            continue

        header = line[1:].strip()

        fields = header.split("|")

        if len(fields) < 6:
            continue

        gene_id = fields[1].split(".")[0]
        gene_symbol = fields[5]

        if gene_id not in gene_mapping:
            gene_mapping[gene_id] = gene_symbol

# ------------------------------------------------------------
# Map gene symbols
# ------------------------------------------------------------

significant["gene_symbol"] = (
    significant["gene_id_base"]
    .map(gene_mapping)
)

# ------------------------------------------------------------
# Report mapping
# ------------------------------------------------------------

print("\nGene ID mapping")
print("----------------------------------------")

print(
    significant[
        [
            "gene_id",
            "gene_symbol",
            "log2FoldChange",
            "pvalue",
            "padj"
        ]
    ].to_string(index=False)
)

mapped = significant["gene_symbol"].notna()

print("\nMapped genes:", mapped.sum())
print("Unmapped genes:", (~mapped).sum())

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

significant.to_csv(
    output_file,
    index=False
)

print(
    "\nMapped results saved to:",
    output_file
)