from pathlib import Path
import csv
import json
import sys


salmon_directory = Path("results/salmon")
metadata_file = Path("data/metadata.csv")
output_directory = Path("results/expression")

if not salmon_directory.is_dir():
    print(
        "ERROR: Salmon results directory not found:",
        salmon_directory
    )
    sys.exit(1)

if not metadata_file.is_file():
    print(
        "ERROR: Metadata file not found:",
        metadata_file
    )
    sys.exit(1)

output_directory.mkdir(
    parents=True,
    exist_ok=True
)


# ---------------------------------------------------------
# Read sample metadata
# ---------------------------------------------------------

metadata = {}

with open(
    metadata_file,
    "r",
    newline=""
) as file:

    reader = csv.DictReader(file)

    if reader.fieldnames is None:
        print(
            "ERROR: Metadata file has no header"
        )
        sys.exit(1)

    required_columns = {
        "sample_id",
        "sra_accession",
        "geo_accession",
        "condition",
        "timepoint"
    }

    missing_columns = (
        required_columns
        - set(reader.fieldnames)
    )

    if missing_columns:
        print(
            "ERROR: Missing metadata columns:",
            ", ".join(
                sorted(missing_columns)
            )
        )
        sys.exit(1)

    for row in reader:

        sra_accession = (
            row["sra_accession"].strip()
        )

        metadata[sra_accession] = {
            "sample_id":
                row["sample_id"].strip(),

            "geo_accession":
                row["geo_accession"].strip(),

            "condition":
                row["condition"].strip(),

            "timepoint":
                row["timepoint"].strip()
        }


# ---------------------------------------------------------
# Find Salmon sample directories
# ---------------------------------------------------------

sample_directories = sorted(
    path
    for path in salmon_directory.iterdir()
    if path.is_dir()
)

if not sample_directories:

    print(
        "ERROR: No Salmon sample directories found in",
        salmon_directory
    )

    sys.exit(1)


results = []


# ---------------------------------------------------------
# Analyze each Salmon result
# ---------------------------------------------------------

for sample_directory in sample_directories:

    sample = sample_directory.name

    quant_file = (
        sample_directory
        / "quant.sf"
    )

    meta_info_file = (
        sample_directory
        / "aux_info"
        / "meta_info.json"
    )

    if not quant_file.is_file():

        print(
            f"WARNING: quant.sf missing for {sample}"
        )

        continue

    if not meta_info_file.is_file():

        print(
            f"WARNING: meta_info.json missing for {sample}"
        )

        continue


    # -----------------------------------------------------
    # Count transcript records in quant.sf
    # -----------------------------------------------------

    transcript_count = 0

    with open(
        quant_file,
        "r",
        newline=""
    ) as file:

        reader = csv.DictReader(
            file,
            delimiter="\t"
        )

        if reader.fieldnames is None:

            print(
                f"WARNING: No header found in {quant_file}"
            )

            continue

        required_quant_columns = {
            "Name",
            "Length",
            "EffectiveLength",
            "TPM",
            "NumReads"
        }

        missing_quant_columns = (
            required_quant_columns
            - set(reader.fieldnames)
        )

        if missing_quant_columns:

            print(
                f"WARNING: Missing columns in {quant_file}:",
                ", ".join(
                    sorted(missing_quant_columns)
                )
            )

            continue

        for row in reader:
            transcript_count += 1


    if transcript_count == 0:

        print(
            f"WARNING: No transcript records found for {sample}"
        )

        continue


    # -----------------------------------------------------
    # Read Salmon meta_info.json
    # -----------------------------------------------------

    with open(
        meta_info_file,
        "r"
    ) as file:

        salmon_metadata = json.load(file)


    num_processed = salmon_metadata.get(
        "num_processed",
        0
    )

    num_mapped = salmon_metadata.get(
        "num_mapped",
        0
    )

    percent_mapped = salmon_metadata.get(
        "percent_mapped",
        0.0
    )

    detected_library_type = salmon_metadata.get(
        "detected_library_type",
        "unknown"
    )

    salmon_version = salmon_metadata.get(
        "salmon_version",
        "unknown"
    )

    quant_errors = salmon_metadata.get(
        "quant_errors",
        []
    )


    # -----------------------------------------------------
    # Match sample to metadata
    # -----------------------------------------------------

    if sample not in metadata:

        print(
            f"WARNING: No metadata found for {sample}"
        )

        continue

    metadata_row = metadata[sample]


    # -----------------------------------------------------
    # Store result
    # -----------------------------------------------------

    results.append({

        "sample":
            sample,

        "sample_id":
            metadata_row["sample_id"],

        "geo_accession":
            metadata_row["geo_accession"],

        "condition":
            metadata_row["condition"],

        "timepoint":
            metadata_row["timepoint"],

        "transcripts":
            transcript_count,

        "processed_fragments":
            num_processed,

        "mapped_fragments":
            num_mapped,

        "mapping_rate_percent":
            round(
                percent_mapped,
                2
            ),

        "library_type":
            detected_library_type,

        "salmon_version":
            salmon_version,

        "quant_errors":
            len(quant_errors)
    })


# ---------------------------------------------------------
# Verify expected samples
# ---------------------------------------------------------

expected_samples = {
    "SRR11884716",
    "SRR11884717",
    "SRR11884718",
    "SRR11884719",
    "SRR11884720",
    "SRR11884721"
}

found_samples = {
    row["sample"]
    for row in results
}

missing_samples = (
    expected_samples
    - found_samples
)

if missing_samples:

    print(
        "\nWARNING: Expected samples missing:",
        ", ".join(
            sorted(missing_samples)
        )
    )


# ---------------------------------------------------------
# Save summary
# ---------------------------------------------------------

output_file = (
    output_directory
    / "salmon_quantification_summary.csv"
)

fieldnames = [

    "sample",

    "sample_id",

    "geo_accession",

    "condition",

    "timepoint",

    "transcripts",

    "processed_fragments",

    "mapped_fragments",

    "mapping_rate_percent",

    "library_type",

    "salmon_version",

    "quant_errors"
]


with open(
    output_file,
    "w",
    newline=""
) as file:

    writer = csv.DictWriter(
        file,
        fieldnames=fieldnames
    )

    writer.writeheader()

    writer.writerows(
        results
    )


# ---------------------------------------------------------
# Print summary
# ---------------------------------------------------------

print(
    "\nSalmon Quantification Summary"
)

print(
    "=" * 90
)

for row in results:

    print(

        f"{row['sample']} | "

        f"{row['condition']} | "

        f"processed: "
        f"{row['processed_fragments']:,} | "

        f"mapped: "
        f"{row['mapped_fragments']:,} | "

        f"mapping: "
        f"{row['mapping_rate_percent']:.2f}%"
    )


print(
    "\nSamples analyzed:",
    len(results)
)

print(
    "Summary saved to:",
    output_file
)