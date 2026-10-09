
import pandas as pd
from pathlib import Path

# 1. Define input and output files
project_dir = Path.cwd()

input_file = (
    project_dir
    / "results"
    / "pathway_analysis"
    / "hallmark_gsea_significant.csv"
)

output_file = (
    project_dir
    / "results"
    / "pathway_analysis"
    / "hallmark_gsea_interpretation.csv"
)

# 2. Check that the input file exists
if not input_file.exists():
    raise FileNotFoundError(
        f"Could not find the GSEA results file: {input_file}\n"
        "Run this script from the rna-seq-analysis project folder."
    )

# 3. Load the significant pathway results
gsea = pd.read_csv(input_file)

required_columns = ["pathway", "NES", "pval", "padj"]

missing_columns = [
    column for column in required_columns
    if column not in gsea.columns
]

if missing_columns:
    raise ValueError(
        f"Missing columns in the input file: {missing_columns}"
    )

if gsea.empty:
    raise ValueError("The GSEA results file contains no pathways.")

# 4. Convert pathway names into readable labels
gsea["pathway_name"] = (
    gsea["pathway"]
    .str.replace("^HALLMARK_", "", regex=True)
    .str.replace("_", " ", regex=False)
    .str.title()
)

# 5. Determine the direction of enrichment
# Positive NES = enriched toward infected samples.
# Negative NES = enriched toward control samples.
gsea["direction"] = gsea["NES"].apply(
    lambda value: (
        "Enriched in infected" if value > 0
        else "Enriched in control" if value < 0
        else "No directional enrichment"
    )
)

# 6. Add cautious biological descriptions
interpretations = {
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB":
        "Inflammatory signaling and NF-kB-associated transcription.",

    "HALLMARK_HEDGEHOG_SIGNALING":
        "Cell signaling involved in development and differentiation.",

    "HALLMARK_OXIDATIVE_PHOSPHORYLATION":
        "Mitochondrial energy production and respiratory metabolism.",

    "HALLMARK_MYC_TARGETS_V1":
        "MYC-associated transcription, growth, and proliferation.",

    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION":
        "Cellular programs associated with epithelial state and cell adhesion.",

    "HALLMARK_ADIPOGENESIS":
        "Lipid storage and metabolic regulation-associated genes.",

    "HALLMARK_E2F_TARGETS":
        "Cell-cycle progression and DNA replication-associated genes.",

    "HALLMARK_PROTEIN_SECRETION":
        "Protein processing, trafficking, and secretion.",

    "HALLMARK_GLYCOLYSIS":
        "Glucose utilization and glycolytic energy metabolism.",

    "HALLMARK_FATTY_ACID_METABOLISM":
        "Fatty-acid utilization and lipid metabolism.",

    "HALLMARK_COMPLEMENT":
        "Innate immune and complement-associated processes.",

    "HALLMARK_MTORC1_SIGNALING":
        "Nutrient sensing, protein synthesis, and metabolic regulation.",

    "HALLMARK_COAGULATION":
        "Coagulation and hemostasis-associated genes.",

    "HALLMARK_XENOBIOTIC_METABOLISM":
        "Metabolism of foreign compounds and detoxification.",

    "HALLMARK_KRAS_SIGNALING_UP":
        "Genes associated with KRAS-related signaling states.",

    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE":
        "Cellular responses to protein-folding and ER stress.",

    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY":
        "Cellular responses associated with oxidative stress.",

    "HALLMARK_UV_RESPONSE_DN":
        "Genes associated with a cellular stress-response signature."
}

gsea["biological_interpretation"] = (
    gsea["pathway"].map(interpretations)
    .fillna("Pathway signature requiring further investigation.")
)

# 7. Label statistical significance
gsea["significance"] = gsea["padj"].apply(
    lambda value: (
        "Highly significant" if value < 0.001
        else "Significant"
    )
)

# 8. Select and sort the final table
output = gsea[
    [
        "pathway",
        "pathway_name",
        "direction",
        "NES",
        "pval",
        "padj",
        "significance",
        "biological_interpretation"
    ]
].copy()

output = output.sort_values("padj", ascending=True)

# 9. Save the interpretation table
output.to_csv(output_file, index=False)

# 10. Print a summary
print()
print("GSEA interpretation table created successfully.")
print(f"Significant pathways: {len(output)}")
print(f"Output file: {output_file}")
print()
print("Pathway summary:")
print("-" * 80)

for _, row in output.iterrows():
    print(
        f"{row['pathway_name']}: "
        f"{row['direction']}; "
        f"NES={row['NES']:.3f}; "
        f"adjusted p-value={row['padj']:.3e}"
    )

print()
print("Important: enrichment indicates an association in this dataset,")
print("not proof that infection caused the pathway change.")