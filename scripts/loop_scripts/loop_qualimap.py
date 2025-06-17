import os
import subprocess
import glob

genome_mapping = False
force_output = False
qualimap_only = True

# Define directory paths (adjust these based on your data structure)
if genome_mapping:
    assemblies_dir = "results/assembly_results/stringent/filtered_genomes/"         # Folder with assembled genomes
    reads_dir = "/mnt/d/Fastq_files/trimmed_fastq/stringent/"                   # Folder with Illumina reads (assumes R1 and R2 fastq pairs)
    output_dir = "results/assembly_results/stringent/QC/filtered/"                 # Folder to store BAM files and Qualimap reports
else:
    assemblies_dir = "results/assembly_results/stringent/filtered_genomes/"   # "results/assembly_results/RL_stammen/"       # Folder with assembled genomes
    reads_dir = "/mnt/d/Fastq_files/trimmed_fastq/stringent/"   # "/mnt/d/Fastq_files/RL_stammen/trimmed/"                  # Folder with Illumina reads (assumes R1 and R2 fastq pairs)
    output_dir = "results/gene_cluster_genes/mapping_output_RL28/"                 # Folder to store BAM files and Qualimap reports
    reference_seq = "results/plasmid/plasmid_bella.fasta"

# Create output directories if they don’t exist
os.makedirs(f"{output_dir}/bams", exist_ok=True)
os.makedirs(f"{output_dir}/qualimap_reports", exist_ok=True)

# Function to run shell commands
def run_command(command):
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e.cmd}\n{e.output}")


# Process each assembly file in the assemblies_dir
for assembly in os.listdir(assemblies_dir):
    if assembly.endswith(".fasta"):
        sample = os.path.splitext(assembly)[0]
        print(f"Processing sample: {sample}")

        reads_R1 = glob.glob(f"{reads_dir}{sample}*_R1*.fastq")[0]
        reads_R2 = glob.glob(f"{reads_dir}{sample}*_R2*.fastq")[0]

        # Check if read files exist or if has already been run (existing qualimap file)
        if not (os.path.isfile(reads_R1) and os.path.isfile(reads_R2)):
            print(f"Read files not found for {sample}, skipping...")
            continue
        if os.path.exists(f"{output_dir}/qualimap_reports/{sample}_report") and not force_output:
            print("Qualimap has already been processed")
            continue

        # Index the assembly
        if genome_mapping:
            run_command(f"bwa index {os.path.join(assemblies_dir, assembly)}")

        # Map reads to the assembly
        sam_path = f"{output_dir}bams/{sample}.sam"
        if genome_mapping:
            run_command(f"bwa-mem2 mem -t 8 {os.path.join(assemblies_dir, assembly)} {reads_R1} {reads_R2} > {sam_path}")
        elif qualimap_only:
            print(f"Executing qualimap only on {sample}")
        else:
            run_command(f"bwa-mem2 mem -t 8 {reference_seq} {reads_R1} {reads_R2} > {sam_path}")

        # Convert SAM to BAM, sort, and index
        bam_path = f"{output_dir}bams/{sample}_sorted.bam"
        if not qualimap_only:
            run_command(f"samtools view -Sb {sam_path} | samtools sort -o {bam_path}")
            run_command(f"samtools index {bam_path}")

        # Run Qualimap on the sorted BAM file
        qualimap_output_dir = f"{output_dir}/qualimap_reports/{sample}_report"
        run_command(f"qualimap bamqc -bam {bam_path} -outdir {qualimap_output_dir}")

        # Optional: remove SAM file to save space
        if not qualimap_only:
            os.remove(sam_path)
        print(f"Finished processing {sample}. Output available in {qualimap_output_dir}")

print("All samples processed.")
