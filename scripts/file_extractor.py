import os
import subprocess
import re

# Set input and output dirs
input_dir = "/mnt/d/Fastq_files/assembly/stringent/"
output_dir = "/mnt/c/Users/freek/Documents/BioSb/Internship/PangenomeProject/results/assembly_results/stringent/scaffolds/"
output_dir_graphs = "/mnt/c/Users/freek/Documents/BioSb/Internship/PangenomeProject/results/assembly_results/stringent/graphs/"
# assembly_graph_with_scaffolds.gfa
os.makedirs(output_dir, exist_ok=True)
os.makedirs(output_dir_graphs, exist_ok=True)

# Loop over all input files
for file_dir in os.listdir(input_dir):
    in1 = os.path.join(input_dir, file_dir)
    filename = re.match(r'^[^_]+_[^_]+', file_dir).group(0)
    filename_graphs = filename + "_assembly_graph.gfa"
    filename = filename + ".fasta"
    out = os.path.join(output_dir, filename)
    out_graphs = os.path.join(output_dir_graphs, filename_graphs)

    # Now run copy scaffolds.fasta and scaffolds.gfa files
    if os.path.exists(out):
        print("Assembly already exists")
        continue
    command = f"cp {in1}/scaffolds.fasta {out}"
    subprocess.run(command, shell=True)

    command_graphs = f"cp {in1}/assembly_graph_with_scaffolds.gfa {out_graphs}"
    subprocess.run(command_graphs, shell=True)
