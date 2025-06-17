import os
import subprocess

# Set input and output dirs
# input_dir = "/mnt/d/Fastq_files/BWS_whole_genome_sequence_raw_reads/"
# output_dir = "/mnt/d/Fastq_files/trimmed_fastq/stringent/"
input_dir = "/mnt/d/Fastq_files/RL_stammen/raw_reads/todo/"
output_dir = "/mnt/d/Fastq_files/RL_stammen/trimmed/"

os.makedirs(output_dir, exist_ok=True)

# Loop over all input files
for filename in os.listdir(input_dir):
    if "_R1_" in filename and filename.endswith(".fastq"):
        in1 = os.path.join(input_dir, filename)
        out1 = os.path.join(output_dir, filename)
        in2 = in1.replace("_R1_", "_R2_")
        out2 = out1.replace("_R1_", "_R2_")

        # Now run fastp command
        if os.path.exists(in2):
            # default: fastp --in1 {in1} --out1 {out1} --in2 {in2} --out2 {out2} -w 10
            command = f"fastp --in1 {in1} --out1 {out1} --in2 {in2} --out2 {out2} -w 6 -q 30 -l 100 -D"
            print(f"Running command: {command}")
            subprocess.run(command, shell=True)
        else:
            print(f"Matching _R2_ file not found for {in1}")
    else:
        print(f"{filename} does not include '_R1_' or is not a fastq file")
