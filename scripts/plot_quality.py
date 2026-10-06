import sys
import matplotlib.pyplot as plt


if len(sys.argv) != 2:
    print("Usage: python3 scripts/plot_quality.py <fastq_file>")
    sys.exit(1)


fastq_file = sys.argv[1]

position_quality_sum = []
position_quality_count = []


with open(fastq_file, "r") as file:

    while True:

        header = file.readline()
        sequence = file.readline()
        separator = file.readline()
        quality = file.readline()

        if not quality:
            break

        quality = quality.strip()

        for position, score in enumerate(quality):

            phred_score = ord(score) - 33

            if position >= len(position_quality_sum):
                position_quality_sum.append(0)
                position_quality_count.append(0)

            position_quality_sum[position] += phred_score
            position_quality_count[position] += 1


positions = []
average_quality = []


for position, total in enumerate(position_quality_sum, start=1):

    positions.append(position)

    average_quality.append(
        total / position_quality_count[position - 1]
    )


plt.figure(figsize=(10, 6))

plt.plot(
    positions,
    average_quality,
    marker="o",
    markersize=3
)

plt.xlabel("Read position")
plt.ylabel("Average Phred quality")
plt.title(f"Per-position sequencing quality: {fastq_file}")

plt.axhline(
    y=20,
    linestyle="--",
    label="Q20 threshold"
)

plt.legend()

plt.tight_layout()

output_file = (
    f"results/qc/{fastq_file.split('/')[-1]}_quality.png"
)

plt.savefig(output_file, dpi=300)

print("Plot saved to:", output_file)