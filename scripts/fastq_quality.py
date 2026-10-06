
from pathlib import Path
import sys
import csv

# Check that a FASTQ file was provided
if len(sys.argv) != 2:
    print("Usage: python scripts/fastq_quality.py <fastq_file>")
    sys.exit(1)

fastq_file = Path(sys.argv[1])

# Check that the input file exists
if not fastq_file.is_file():
    print("Error: FASTQ file not found:", fastq_file)
    sys.exit(1)

print("Analyzing:", fastq_file)

# Initialize summary statistics
read_count = 0
total_bases = 0
total_quality = 0
low_quality_bases = 0
gc_bases = 0

# Store running quality totals for each read position
position_quality_sum = []
position_quality_count = []

# Read the FASTQ file one record at a time
with open(fastq_file, "r") as file:
    while True:
        header = file.readline()

        if not header:
            break

        sequence = file.readline().strip()
        separator = file.readline()
        quality = file.readline().strip()

        # Check that the FASTQ record is complete
        if not separator or not quality:
            raise ValueError("Incomplete FASTQ record encountered")

        if len(sequence) != len(quality):
            raise ValueError("Sequence and quality lengths do not match")

        read_count += 1
        total_bases += len(sequence)

        # Calculate GC content
        sequence_upper = sequence.upper()
        gc_bases += sequence_upper.count("G")
        gc_bases += sequence_upper.count("C")

        # Calculate base quality statistics
        for position, score in enumerate(quality):
            phred_score = ord(score) - 33

            total_quality += phred_score

            if phred_score < 20:
                low_quality_bases += 1

            if position >= len(position_quality_sum):
                position_quality_sum.append(0)
                position_quality_count.append(0)

            position_quality_sum[position] += phred_score
            position_quality_count[position] += 1

# Avoid division by zero if the file contains no reads
if read_count == 0 or total_bases == 0:
    raise ValueError("No reads found in the FASTQ file")

# Calculate summary metrics
average_read_length = total_bases / read_count
average_phred = total_quality / total_bases
low_quality_percentage = (low_quality_bases / total_bases) * 100
gc_percentage = (gc_bases / total_bases) * 100

# Print the results
print("\nQC Summary")
print("-" * 35)
print("Number of reads:", read_count)
print("Total bases:", total_bases)
print("Average read length:", round(average_read_length, 2))
print("Average Phred quality:", round(average_phred, 2))
print("Bases below Phred 20:", round(low_quality_percentage, 2), "%")
print("GC content:", round(gc_percentage, 2), "%")

print("\nPer-position average quality:")

for position, total in enumerate(position_quality_sum, start=1):
    average_position_quality = (
        total / position_quality_count[position - 1]
    )

    print(
        "Position",
        position,
        ":",
        round(average_position_quality, 2)
    )

# Create the QC output directory
output_directory = Path("results/qc")
output_directory.mkdir(parents=True, exist_ok=True)

# Write summary statistics to CSV
output_file = output_directory / f"{fastq_file.stem}_qc.csv"

with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)

    writer.writerow(["metric", "value"])
    writer.writerow(["input_file", str(fastq_file)])
    writer.writerow(["number_of_reads", read_count])
    writer.writerow(["total_bases", total_bases])
    writer.writerow(["average_read_length", round(average_read_length, 2)])
    writer.writerow(["average_phred_quality", round(average_phred, 2)])
    writer.writerow([
        "bases_below_q20_percent",
        round(low_quality_percentage, 2)
    ])
    writer.writerow(["gc_content_percent", round(gc_percentage, 2)])

# Write per-position quality to a separate CSV
position_output_file = (
    output_directory / f"{fastq_file.stem}_position_quality.csv"
)

with open(position_output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)

    writer.writerow(["position", "average_phred_quality"])

    for position, total in enumerate(position_quality_sum, start=1):
        average_position_quality = (
            total / position_quality_count[position - 1]
        )

        writer.writerow([
            position,
            round(average_position_quality, 2)
        ])

print("\nQC summary saved to:", output_file)
print("Position quality saved to:", position_output_file)