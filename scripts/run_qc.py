from pathlib import Path
import subprocess
import sys

# Location of the raw FASTQ files
raw_directory = Path("data/raw")

# Location of the QC script
qc_script = Path("scripts/fastq_quality.py")

# Find every FASTQ file in data/raw
fastq_files = sorted(raw_directory.glob("*.fastq"))

# Stop if no FASTQ files are found
if not fastq_files:
    print("No FASTQ files found in", raw_directory)
    sys.exit(1)

print("Found", len(fastq_files), "FASTQ files")

# Run the QC script on every FASTQ file
for fastq_file in fastq_files:

    print("\n" + "=" * 60)
    print("Processing:", fastq_file)
    print("=" * 60)

    result = subprocess.run([
        sys.executable,
        str(qc_script),
        str(fastq_file)
    ])

    # Stop immediately if QC fails
    if result.returncode != 0:
        print("QC failed for:", fastq_file)
        sys.exit(result.returncode)

print("\nQC completed for all FASTQ files.")