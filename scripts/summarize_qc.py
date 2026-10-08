from pathlib import Path
import csv


qc_directory = Path("results/qc")
output_file = qc_directory / "all_samples_qc_summary.csv"


if not qc_directory.is_dir():
    raise FileNotFoundError(
        f"QC directory not found: {qc_directory}"
    )


qc_files = sorted(
    qc_directory.glob("*_qc.csv")
)


if not qc_files:
    raise FileNotFoundError(
        "No QC summary CSV files found"
    )


rows = []


for qc_file in qc_files:

    metrics = {}

    with open(qc_file, "r", newline="") as file:

        reader = csv.DictReader(file)

        for row in reader:
            metrics[row["metric"]] = row["value"]


    rows.append({
        "sample": qc_file.stem.replace("_qc", ""),
        "number_of_reads": metrics["number_of_reads"],
        "total_bases": metrics["total_bases"],
        "average_read_length": metrics["average_read_length"],
        "average_phred_quality": metrics["average_phred_quality"],
        "bases_below_q20_percent": metrics[
            "bases_below_q20_percent"
        ],
        "gc_content_percent": metrics["gc_content_percent"]
    })


with open(output_file, "w", newline="") as file:

    fieldnames = [
        "sample",
        "number_of_reads",
        "total_bases",
        "average_read_length",
        "average_phred_quality",
        "bases_below_q20_percent",
        "gc_content_percent"
    ]

    writer = csv.DictWriter(
        file,
        fieldnames=fieldnames
    )

    writer.writeheader()
    writer.writerows(rows)


print(
    "QC files found:",
    len(qc_files)
)

print(
    "Summary saved to:",
    output_file
)

print("\nQC summary")
print("-" * 90)

for row in rows:

    print(
        row["sample"],
        "| reads:",
        row["number_of_reads"],
        "| avg length:",
        row["average_read_length"],
        "| avg Phred:",
        row["average_phred_quality"],
        "| Q<20:",
        row["bases_below_q20_percent"] + "%",
        "| GC:",
        row["gc_content_percent"] + "%"
    )