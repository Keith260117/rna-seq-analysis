from pathlib import Path
import subprocess
import sys

# ------------------------------------------------------------
# Directory locations
# ------------------------------------------------------------

raw_directory = Path("data/raw")
index_directory = Path("reference/salmon_index")
output_directory = Path("results/salmon")

# ------------------------------------------------------------
# Samples
# ------------------------------------------------------------

samples = [
    "SRR11884716",
    "SRR11884717",
    "SRR11884718",
    "SRR11884719",
    "SRR11884720",
    "SRR11884721",
]

# ------------------------------------------------------------
# Check directories
# ------------------------------------------------------------

if not raw_directory.is_dir():
    print("ERROR: Raw data directory not found:", raw_directory)
    sys.exit(1)

if not index_directory.is_dir():
    print("ERROR: Salmon index not found:", index_directory)
    sys.exit(1)

output_directory.mkdir(
    parents=True,
    exist_ok=True
)

# ------------------------------------------------------------
# Run Salmon for each sample
# ------------------------------------------------------------

for sample in samples:

    read1 = raw_directory / f"{sample}_1.fastq"
    read2 = raw_directory / f"{sample}_2.fastq"
    sample_output = output_directory / sample

    if not read1.is_file():
        print(f"ERROR: Missing Read 1 for {sample}")
        sys.exit(1)

    if not read2.is_file():
        print(f"ERROR: Missing Read 2 for {sample}")
        sys.exit(1)

    if (sample_output / "quant.sf").is_file():
        print(
            f"\nSkipping {sample}: "
            "quant.sf already exists"
        )
        continue

    sample_output.mkdir(
        parents=True,
        exist_ok=True
    )

    print("\n" + "=" * 70)
    print(f"Quantifying {sample}")
    print("=" * 70)

    command = [
        "salmon",
        "quant",
        "-i",
        str(index_directory),
        "-l",
        "A",
        "-1",
        str(read1),
        "-2",
        str(read2),
        "-p",
        "4",
        "-o",
        str(sample_output),
    ]

    result = subprocess.run(command)

    if result.returncode != 0:
        print(
            f"\nERROR: Salmon failed for {sample}"
        )
        sys.exit(result.returncode)

    print(f"\nCompleted: {sample}")

# ------------------------------------------------------------
# Finished
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("Salmon quantification completed for all samples")
print("=" * 70)