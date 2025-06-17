import os
import subprocess
from pathlib import Path
import pandas as pd
from Bio import SeqIO


def count_reads(file_path):
    result = subprocess.run(f"pigz -cd {file_path} | wc -l", shell=True, capture_output=True, text=True)
    return int(result.stdout.strip()) // 4  # Each read has 4 lines


def trim_reads(input_dir, output_dir):
    """Trim raw reads using fastp and rename output files."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    for file in os.listdir(input_dir):
        if file.endswith(".fastq.gz"):
            # Extract sample name and read identifier
            parts = file.split("_")
            sample_name = parts[0]
            read_id = parts[4]  # R1 or R2
            if read_id == "R2": continue  # fastp for PE, so only run once

            input_path1 = os.path.join(input_dir, file)
            input_path2 = os.path.join(input_dir, file.replace("R1", "R2"))
            output_name1 = f"{sample_name}_R1.fastq.gz"
            output_path1 = os.path.join(output_dir, output_name1)
            output_name2 = f"{sample_name}_R2.fastq.gz"
            output_path2 = os.path.join(output_dir, output_name2)
            
            if os.path.exists(output_path1) and os.path.exists(output_path2):
                print(f"{sample_name} already exists")
                initial_reads = count_reads(input_path1)
                trimmed_reads = count_reads(output_path1)
                df_count.loc[sample_name, "initial_reads"] = initial_reads * 2  # multiply by two because 2 reads, R1 and R2
                df_count.loc[sample_name, "trimmed_reads"] = trimmed_reads * 2
                continue

            # Run fastp
            subprocess.run(f"mamba run -n LCPA fastp -i {input_path1} -I {input_path2} -o {output_path1} -O {output_path2} --compression 6 "
                           f"-w 6 -q 30 -l 100 -D"
            , shell=True)
            initial_reads = count_reads(input_path1)
            trimmed_reads = count_reads(output_path1)
            df_count.loc[sample_name, "initial_reads"] = initial_reads * 2  # multiply by two because 2 reads, R1 and R2
            df_count.loc[sample_name, "trimmed_reads"] = trimmed_reads * 2
            print(f"Trimmed {file} -> {output_name1}")


def filter_human(input_dir, bwa_index, output_dir_filtered):
    Path(output_dir_filtered).mkdir(parents=True, exist_ok=True)
    for file in os.listdir(input_dir):
        if file.endswith(".fastq.gz"):
            parts = file.split("_")
            sample_name_bam = parts[0] + ".bam"
            sample_name_fastq = parts[0] + ".fastq.gz"
            sample_name1 = parts[0] + "_R1.fastq.gz"
            sample_name2 = parts[0] + "_R2.fastq.gz"

            input_path1 = os.path.join(input_dir, sample_name1)
            input_path2 = os.path.join(input_dir, sample_name2)
            output_path_filtered_bam = os.path.join(output_dir_filtered, sample_name_bam)
            output_path_filtered = os.path.join(output_dir_filtered, sample_name_fastq)

            if os.path.exists(output_path_filtered):
                print(f"{file} already exists")
                filtered_reads = count_reads(output_path_filtered)
                df_count.loc[parts[0], "filtered_reads"] = filtered_reads
                continue

            # Run bwa to map human reads
            if not os.path.exists(output_path_filtered_bam):
                print(f"Filtering human reads for {parts[0]}")
                subprocess.run(f"bwa mem -t 8 {bwa_index} {input_path1} {input_path2}| samtools view -o {output_path_filtered_bam}", shell=True)

            # Get unmatched reads
            print(f"Getting unmapped reads (non-human) for {parts[0]}")
            subprocess.run(f"samtools fastq -f 4 {output_path_filtered_bam} -s unpaired.fastq | pigz > {output_path_filtered}", shell=True)
            filtered_reads = count_reads(output_path_filtered)
            df_count.loc[parts[0], "filtered_reads"] = filtered_reads


def assemble_reads(input_dir, output_dir):
    """Assemble reads using metaSPAdes."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    for file in os.listdir(input_dir):
        if file.endswith(".fastq.gz"):
            sample_name = file.split(".fastq")[0]
            output_path = os.path.join(output_dir, sample_name)
            if os.path.exists(output_path):
                print(f"{file} already exists")
                continue
            Path(output_path).mkdir(parents=True, exist_ok=True)

            subprocess.run(f"mamba run -n metagenome spades.py --meta -t 6 -m 10 --12 {input_dir}{file} -o {output_path}", shell=True)
            print(f"Assembly completed for {sample_name} -> {output_path}")


def run_phabox2(input_dir, output_dir, phabox_db):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    for sample_name in os.listdir(input_dir):
        input_path = os.path.join(input_dir, sample_name, "contigs.fasta")
        sample_output_dir = os.path.join(output_dir, sample_name)
        if os.path.exists(sample_output_dir):
            print(f"{sample_name} already exists")
            try:
                with open(f"{sample_output_dir}/filtered_contigs.fa", "r") as f:
                    read_count = 0
                    for record in SeqIO.parse(f, "fasta"):
                        header = record.id.split("_")
                        length = int(header[3])
                        coverage = float(header[5])
                        read_count += length * coverage
                df_count.loc[
                    sample_name, "virus_reads"] = read_count / 150  # coverage * length = bp, so divide by 150 for read count estimate
            except FileNotFoundError:
                df_count.loc[sample_name, "virus_reads"] = 0
            continue

        subprocess.run(f"mamba run -n phabox2 phabox2 --task end_to_end --dbdir {phabox_db} \
        --outpth  {sample_output_dir} --contigs {input_path} --threads 8 --len 500", shell=True)
        try:
            with open(f"{sample_output_dir}/filtered_contigs.fa", "r") as f:
                read_count = 0
                for record in SeqIO.parse(f, "fasta"):
                    header = record.id.split("_")
                    length = int(header[3])
                    coverage = float(header[5])
                    read_count += length * coverage
            df_count.loc[sample_name, "virus_reads"] = read_count / 150  # coverage * length = bp, so divide by 150 for read count estimate
        except FileNotFoundError:
            df_count.loc[sample_name, "virus_reads"] = 0
        print(f"Phabox2 completed for {sample_name} -> {sample_output_dir}")


if __name__ == "__main__":
    # Directories
    raw_reads_dir = "data/metagenome/reads/"
    trimmed_reads_dir = "data/metagenome/trimmed_reads/"
    contigs_dir = "data/metagenome/contigs/"
    phabox2_output_dir = "data/metagenome/phabox2_output/"
    phabox2_db_dir = "data/metagenome/phabox_db_v2/"
    bwa_index_path = "data/metagenome/human_genome/genome.fa"
    filtered_reads_dir = "data/metagenome/filtered_reads/"
    read_counts_out = "data/metagenome/read_counts.xlsx"

    # Count reads dataframe
    df_count = pd.DataFrame()

    # Pipeline
    print("Trimming reads...")
    trim_reads(raw_reads_dir, trimmed_reads_dir)

    print("Filtering reads...")
    filter_human(trimmed_reads_dir, bwa_index_path, filtered_reads_dir)

    print("Assembling reads...")
    assemble_reads(filtered_reads_dir, contigs_dir)

    print("Running Phabox2...")
    run_phabox2(contigs_dir, phabox2_output_dir, phabox2_db_dir)

    print(f"Saving read counts to {read_counts_out}...")
    df_count.to_excel(read_counts_out)

    print("Pipeline complete!")
