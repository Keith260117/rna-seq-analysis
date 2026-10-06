
# RNA-seq Quality Control

## Sample

- SRA accession: SRR11884716
- Data type: paired-end RNA sequencing
- Read length: 51 bp

## Metrics calculated

The custom Python script calculates:

- Total read count
- Total number of bases
- Average read length
- Mean Phred quality score
- Percentage of bases below Q20
- GC content
- Average Phred quality at each read position

## Quality score calculation

The script uses Phred+33 encoding:

Phred score = ASCII character value - 33

A Phred Q20 score corresponds to an estimated base-call error probability of 1%.

## Output files

- `*_qc.csv`: summary statistics
- `*_position_quality.csv`: average quality at each read position
- `*.fastq_quality.png`: quality-by-position plots

## Limitations

Mean base quality and GC content alone cannot establish that a library is free from contamination, adapter sequences, or other technical biases. Further assessment and alignment-based quality checks are required.