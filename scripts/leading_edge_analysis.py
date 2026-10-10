from pathlib import Path
import pandas as pd

# ------------------------------------------------------------
# File locations
# ------------------------------------------------------------

project_dir = Path.cwd()

gsea_file = (
    project_dir
    / "results"
    / "pathway_analysis"
    / "hallmark_gsea_results.csv"
)

deseq_file = (
    project_dir
    / "results"
    / "expression"
    / "DESeq2_infected_vs_control.csv"
)

reference_file = (
    project_dir
    / "reference"
    / "gencode.v49.transcripts.fa"
)

output_dir = project_dir / "results" / "pathway_analysis"

leading_edge_file = output_dir / "leading_edge_genes.csv"
summary_file = output_dir / "leading_edge_gene_summary.csv"

# ------------------------------------------------------------
# Check inputs
# ------------------------------------------------------------

for file_path in [gsea_file, deseq_file, reference_file]:
    if not file_path.is_file():
        raise FileNotFoundError(f"Required input not found: {file_path}")

output_dir.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Load GSEA and DESeq2 results
# ------------------------------------------------------------

gsea = pd.read_csv(gsea_file)
deseq = pd.read_csv(deseq_file)

required_gsea = {"pathway", "NES", "padj", "leadingEdge"}
required_deseq = {"gene_id", "log2FoldChange", "stat", "padj"}

if not required_gsea.issubset(gsea.columns):
    raise ValueError(f"GSEA file must contain: {sorted(required_gsea)}")

if not required_deseq.issubset(deseq.columns):
    raise ValueError(f"DESeq2 file must contain: {sorted(required_deseq)}")

# Analyze pathways meeting the existing adjusted p-value threshold.
gsea = gsea[gsea["padj"].notna() & (gsea["padj"] < 0.05)].copy()

# ------------------------------------------------------------
# Build Ensembl gene ID -> gene symbol mapping from GENCODE
# ------------------------------------------------------------

gene_mapping = {}

with open(reference_file, "r") as fasta:
    for line in fasta:
        if not line.startswith(">"):
            continue

        fields = line[1:].strip().split("|")

        if len(fields) < 6:
            continue

        gene_id = fields[1].split(".")[0]
        gene_symbol = fields[5]

        if gene_id and gene_symbol and gene_symbol != ".":
            gene_mapping.setdefault(gene_id, gene_symbol)

print("GENCODE gene IDs available for mapping:", len(gene_mapping))

# ------------------------------------------------------------
# Prepare DESeq2 results
# ------------------------------------------------------------

deseq = deseq.copy()
deseq["gene_id_base"] = deseq["gene_id"].astype(str).str.split(".").str[0]

# If duplicate base IDs occur, retain the row with the smallest
# non-missing adjusted p-value.
deseq = deseq.sort_values("padj", na_position="last")
deseq = deseq.drop_duplicates("gene_id_base", keep="first")

deseq_lookup = deseq[
    ["gene_id_base", "log2FoldChange", "stat", "padj"]
].copy()

# ------------------------------------------------------------
# Expand each pathway's leading-edge list into individual rows
# ------------------------------------------------------------

records = []

for _, row in gsea.iterrows():
    gene_ids = str(row["leadingEdge"]).split(";")

    for gene_id in gene_ids:
        gene_id = gene_id.strip()

        if not gene_id or gene_id.lower() == "nan":
            continue

        gene_id_base = gene_id.split(".")[0]

        records.append({
            "pathway": row["pathway"],
            "NES": row["NES"],
            "pathway_padj": row["padj"],
            "gene_id": gene_id,
            "gene_id_base": gene_id_base,
        })

leading = pd.DataFrame(records)

if leading.empty:
    raise ValueError("No leading-edge genes were found in significant pathways.")

leading["gene_symbol"] = leading["gene_id_base"].map(gene_mapping)

leading = leading.merge(
    deseq_lookup,
    on="gene_id_base",
    how="left",
    validate="many_to_one",
)

leading = leading[
    [
        "pathway",
        "NES",
        "pathway_padj",
        "gene_id",
        "gene_symbol",
        "log2FoldChange",
        "stat",
        "padj",
    ]
].rename(columns={"padj": "gene_padj"})

leading.to_csv(leading_edge_file, index=False)

# ------------------------------------------------------------
# Summarize genes recurring across pathways
# ------------------------------------------------------------

mapped = leading[leading["gene_symbol"].notna()].copy()

summary = (
    mapped.groupby(["gene_id", "gene_symbol"], as_index=False)
    .agg(
        pathway_count=("pathway", "nunique"),
        pathways=("pathway", lambda values: ";".join(sorted(set(values)))),
        log2FoldChange=("log2FoldChange", "first"),
        stat=("stat", "first"),
        gene_padj=("gene_padj", "first"),
    )
    .sort_values(
        ["pathway_count", "gene_padj"],
        ascending=[False, True],
        na_position="last",
    )
)

summary.to_csv(summary_file, index=False)

# ------------------------------------------------------------
# Report
# ------------------------------------------------------------

print("\nLeading-edge analysis complete")
print("----------------------------------------")
print("Significant pathways analyzed:", gsea["pathway"].nunique())
print("Pathway-gene rows:", len(leading))
print("Unique leading-edge genes:", leading["gene_id"].nunique())
print("Rows with mapped gene symbols:", leading["gene_symbol"].notna().sum())
print("Unique mapped genes:", mapped["gene_id"].nunique())
print("Leading-edge table:", leading_edge_file)
print("Recurring-gene summary:", summary_file)

print("\nTop genes recurring across pathways:")
print(summary.head(15).to_string(index=False))
