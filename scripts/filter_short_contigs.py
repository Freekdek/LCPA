import os
import subprocess

# Set input and output dirs
input_dir = "/mnt/c/Users/freek/Documents/BioSb/Internship/PangenomeProject/results/assembly_results/stringent/scaffolds/"
output_dir = "/mnt/c/Users/freek/Documents/BioSb/Internship/PangenomeProject/results/assembly_results/stringent/filtered_genomes/"

os.makedirs(output_dir, exist_ok=True)

# Loop over all input files
for filename in os.listdir(input_dir):
    if filename.endswith(".fasta"):
        in1 = os.path.join(input_dir, filename)
        out = os.path.join(output_dir, filename)

        # Now run fastp command
        if os.path.exists(out):
            print("Annotation already exists")
            continue
        command = f"seqtk seq -L 1000 {in1} > {out}"
        print(f"Running command: {command}")
        subprocess.run(command, shell=True)
    else:
        print(f"{filename} is not a fasta file")
