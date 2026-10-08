from pathlib import Path
import csv
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Directories
# ------------------------------------------------------------

qc_directory = Path("results/qc")

if not qc_directory.is_dir():
    raise FileNotFoundError(
        f"QC directory not found: {qc_directory}"
    )

# ------------------------------------------------------------
# Find position-quality CSV files
# ------------------------------------------------------------

position_files = sorted(
    qc_directory.glob("*_position_quality.csv")
)

if not position_files:
    raise FileNotFoundError(
        "No position-quality CSV files found"
    )

print(
    "Found",
    len(position_files),
    "position-quality files"
)

# ------------------------------------------------------------
# Create one plot for each FASTQ file
# ------------------------------------------------------------

for csv_file in position_files:

    positions = []
    average_quality = []

    with open(csv_file, "r", newline="") as file:

        reader = csv.DictReader(file)

        for row in reader:
            positions.append(
                int(row["position"])
            )

            average_quality.append(
                float(row["average_phred_quality"])
            )

    if not positions:
        print("Skipping empty file:", csv_file)
        continue

    output_file = (
        qc_directory
        / f"{csv_file.stem}.png"
    )

    plt.figure(figsize=(10, 6))

    plt.plot(
        positions,
        average_quality,
        marker="o",
        markersize=3
    )

    plt.axhline(
        y=20,
        linestyle="--",
        label="Q20 threshold"
    )

    plt.xlabel("Read position")
    plt.ylabel("Average Phred quality")

    plt.title(
        f"Per-position sequencing quality: {csv_file.stem}"
    )

    plt.legend()
    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=300
    )

    plt.close()

    print(
        "Plot saved to:",
        output_file
    )

# ------------------------------------------------------------
# Finished
# ------------------------------------------------------------

print("\nQC plotting completed successfully")
