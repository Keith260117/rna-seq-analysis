from pathlib import Path
import csv


qc_directory = Path("results/qc")
output_file = qc_directory / "sample_qc_summary.csv"


samples = [
    "SRR11884716",
    "SRR11884717",
    "SRR11884718",
    "SRR11884719",
    "SRR11884720",
    "SRR11884721",
]


rows = []


for sample in samples:

    row = {
        "sample": sample
    }

    for read in ["1", "2"]:

        qc_file = (
            qc_directory
            / f"{sample}_{read}_qc.csv"
        )

        if not qc_file.is_file():
            raise FileNotFoundError(
                f"Missing QC file: {qc_file}"
            )

        metrics = {}

        with open(qc_file, "r", newline="") as file:

            reader = csv.DictReader(file)

            for item in reader:
                metrics[item["metric"]] = item["value"]

        row[f"R{read}_reads"] = metrics[
            "number_of_reads"
        ]

        row[f"R{read}_avg_phred"] = metrics[
            "average_phred_quality"
        ]

        row[f"R{read}_Q20_below_percent"] = metrics[
            "bases_below_q20_percent"
        ]

        row[f"R{read}_GC_percent"] = metrics[
            "gc_content_percent"
        ]

    rows.append(row)


fieldnames = [
    "sample",
    "R1_reads",
    "R2_reads",
    "R1_avg_phred",
    "R2_avg_phred",
    "R1_Q20_below_percent",
    "R2_Q20_below_percent",
    "R1_GC_percent",
    "R2_GC_percent"
]


with open(output_file, "w", newline="") as file:

    writer = csv.DictWriter(
        file,
        fieldnames=fieldnames
    )

    writer.writeheader()
    writer.writerows(rows)


print(
    "Created six-sample QC table:"
)

print(output_file)

print("\n")

print(
    "Sample | R1 Phred | R2 Phred | "
    "R1 Q<20 | R2 Q<20 | R1 GC | R2 GC"
)

print("-" * 85)


for row in rows:

    print(
        row["sample"],
        "|",
        row["R1_avg_phred"],
        "|",
        row["R2_avg_phred"],
        "|",
        row["R1_Q20_below_percent"] + "%",
        "|",
        row["R2_Q20_below_percent"] + "%",
        "|",
        row["R1_GC_percent"] + "%",
        "|",
        row["R2_GC_percent"] + "%"
    )