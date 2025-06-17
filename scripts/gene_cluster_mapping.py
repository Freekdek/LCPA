import os
import subprocess
import glob
from collections import defaultdict
import pysam
# currently edited so it runs the RL28 fragment as reference for both datasets


force_output = False  # True to force rerunning of all samples

# Directory paths
reads_dir_RL = "/mnt/d/Fastq_files/RL_stammen/trimmed/"
assemblies_dir_RL = "results/assembly_results/RL_stammen/"
reads_dir_vrouwen = "/mnt/d/Fastq_files/trimmed_fastq/stringent/"
assemblies_dir_vrouwen = "results/assembly_results/stringent/filtered_genomes/"
output_dir = "results/gene_cluster_genes/mapping_output_RL28_expanded/"
reference_seq = "results/gene_cluster_genes/CTV_05_fragment.fasta"
# reference_seq = "results/gene_cluster_genes/RL28_fragment_expanded.fasta"
# output_dir = "results/plasmid/output_fastas/"
# reference_seq = "results/plasmid/plasmid_bella.fasta"

# Create output directories if they don’t exist
os.makedirs(f"{output_dir}", exist_ok=True)
os.makedirs(f"{output_dir}/bams", exist_ok=True)
os.makedirs(f"{output_dir}/assembly", exist_ok=True)
os.makedirs(f"{output_dir}/fastas", exist_ok=True)


# Function to run shell commands
def run_command(command):
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e.cmd}\n{e.output}")


def test_glob(reads_dir, assemblies_dir):
    assemblies = os.listdir(assemblies_dir)
    assemblies_list = [asmbly for asmbly in assemblies if asmbly.endswith(".fasta")]

    for reads_R1, assembly in zip(sorted(glob.glob(f"{reads_dir}*R1*.fastq")), assemblies_list):
        sample = os.path.splitext(assembly)[0]
        reads_R2 = reads_R1.replace("_R1", "_R2")

        print(f"Mapping sample: {sample}")
        print("Reads R1:", reads_R1)
        print("Reads R2", reads_R2)
        print("Assembly:", assembly, "\n")


def reference_mapping(reads_dir, assemblies_dir):
    # Process each sample based on read files
    assemblies = os.listdir(assemblies_dir)
    assemblies_list = [asmbly for asmbly in assemblies if asmbly.endswith(".fasta")]
    # Loop each sample and map against reference sequence, output is a bam file which is used for downstream analysis.
    for reads_R1, assembly in zip(sorted(glob.glob(f"{reads_dir}*R1*.fastq")), assemblies_list):
        sample = os.path.splitext(assembly)[0]
        reads_R2 = reads_R1.replace("_R1", "_R2")

        print(f"Mapping sample: {sample}")

        # Check if read files exist
        if not os.path.isfile(reads_R2):
            print(f"Read pair for {sample} not found. Skipping...")
            continue
        if not force_output and os.path.exists(f"{output_dir}fastas/{sample}.fasta"):
            print(f"{sample} has already been processed! Skipping...")
            continue

        # Step 2: Map reads to the reference gene cluster
        sam_path = f"{output_dir}bams/{sample}.sam"
        run_command(f"bwa-mem2 mem -t 8 {reference_seq} {reads_R1} {reads_R2} > {sam_path}")

        # Step 3: Convert SAM to BAM, sort, and index
        bam_path = f"{output_dir}bams/{sample}_sorted.bam"
        run_command(f"samtools view -Sb {sam_path} | samtools sort -o {bam_path}")
        run_command(f"samtools index {bam_path}")

        # Step 4: Extract consensus FASTA from BAM file
        run_command(
            f"samtools consensus {bam_path} --threads 4 -o {output_dir}fastas/{sample}.fasta")

        # Optional: remove SAM file to save space
        os.remove(sam_path)
    return print("All samples mapped.")


def split_read_positions(min_split_read_support=10, min_mapq=30):
    positions_summary_file = f"{output_dir}/insertion_positions_summary.tsv"
    with open(positions_summary_file, "w") as summary:
        summary.write("Sample\tInsertion_Position\tSplit_Read_Count\n")  # Header for summary file
        # Process each sample
        samples_dir = f"{output_dir}bams/"
        for sample in os.listdir(samples_dir):
            if sample.endswith(".bai"):
                continue
            sample_name = sample.split("_sorted.bam")[0]
            bam_file = os.path.join(samples_dir, sample)
            print(f"Processing sample: {sample_name}")

            # Identify split reads
            split_read_positions = defaultdict(int)
            with pysam.AlignmentFile(bam_file, "rb") as bam:
                for read in bam.fetch():
                    if read.has_tag("SA"):  # Split read tag
                        if read.mapping_quality >= min_mapq:
                            position = f"{read.reference_name}:{read.reference_start}"
                            split_read_positions[position] += 1

            # Filter valid insertions and log them
            valid_insertions = [(pos, count) for pos, count in split_read_positions.items() if
                                count >= min_split_read_support]
            for pos, count in valid_insertions:
                summary.write(f"{sample_name}\t{pos}\t{count}\n")

            # Convert bam to FASTA (get the consensus)
            run_command(
                f"samtools-1.21/samtools consensus {bam_file} -a --threads 6 -o {output_dir}fastas/{sample_name}.fasta")
        # Rename header of each FASTA
        for fasta_file in os.listdir(f"{output_dir}fastas/"):
            sample_name = fasta_file.split(".fasta")[0]  # Get filename without extension

            # Read and modify the FASTA file content
            modified_lines = []
            with open(fasta_file, "r") as file:
                for line in file:
                    if line.startswith(">"):  # Update the header
                        modified_lines.append(f">{sample_name}\n")
                    else:  # Keep the sequence as is
                        modified_lines.append(line)

            # Overwrite the original FASTA file
            with open(fasta_file, "w") as file:
                file.writelines(modified_lines)
    return print(f"All samples processed. Summary saved to {positions_summary_file}")


def main():
    # Step 1: Index the reference gene cluster
    run_command(f"bwa-mem2 index {reference_seq}")

    # Step 2, 3 and 4: Map reads to the reference, convert SAM to BAM, extract consensus from BAM
    # test_glob(reads_dir_RL, assemblies_dir_RL)
    reference_mapping(reads_dir_vrouwen, assemblies_dir_vrouwen)
    reference_mapping(reads_dir_RL, assemblies_dir_RL)

    # Step 5: Use split reads to determine possible gene insertions
    # split_read_positions()

    # Step 6: Annotate fasta
    # Step 7: Modify gff files by running BLAST and adding the hits to gff file
    # Step 8: Modify gff files by adding the possible gene insertions based on the split reads to gff file
    # Step 9: Run lovis4u
    # Step 10: Edit lovis4u feature annotation table: changing order of samples and changing colors of annotations


main()
