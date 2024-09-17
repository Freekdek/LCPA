import os
import re


def calculate_coverage_from_contigs(contigs_file):
    total_length = 0
    total_coverage = 0

    with open(contigs_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                # Extract contig length and coverage from the header
                match = re.search(r'_length_(\d+)_cov_([\d\.]+)', line)
                if match:
                    length = int(match.group(1))
                    coverage = float(match.group(2))

                    # Update total length and weighted coverage
                    total_length += length
                    total_coverage += length * coverage

    if total_length == 0:
        return 0

    # Calculate average coverage
    avg_coverage = total_coverage / total_length
    return avg_coverage


def process_assemblies(input_dir):
    # Traverse all subdirectories
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file == "contigs.fasta":
                contigs_path = os.path.join(root, file)
                avg_coverage = calculate_coverage_from_contigs(contigs_path)
                print(f"Average coverage for {contigs_path}: {avg_coverage:.2f}x")


input_directory = "/mnt/d/Fastq_files/assembly/duplicated/"
process_assemblies(input_directory)

