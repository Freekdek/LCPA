import os
import subprocess

# Set input and output dirs
input_dir = "/mnt/d/Fastq_files/trimmed_fastq/stringent/"
output_dir = "/mnt/d/Fastq_files/assembly/stringent/"

os.makedirs(output_dir, exist_ok=True)

# Loop over all input files
for filename in os.listdir(input_dir):
    if "_R1_" in filename and filename.endswith(".fastq"):
        in1 = os.path.join(input_dir, filename)
        in2 = in1.replace("_R1_", "_R2_")
        out = os.path.join(output_dir, filename.replace("_R1_", ""))

        # Now run fastp command
        if os.path.exists(in2):
            if os.path.exists(out):
                print("Assembly already exists")
                continue
            command = f"python3 SPAdes-4.0.0-Linux/bin/spades.py -1 {in1} -2 {in2} -o {out} -t 12 --isolate -m 8"
            print(f"Running command: {command}")
            subprocess.run(command, shell=True)
        else:
            print(f"Matching _R2_ file not found for {in1}")
    else:
        print(f"{filename} does not include '_R1_' or is not a fastq file")
