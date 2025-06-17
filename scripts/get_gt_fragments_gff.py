import os
import subprocess
import pandas as pd
import re

# Inputs
fasta_folder = "results/gene_cluster_genes/mapping_output_RL28_expanded/fastas/"
custom_genes_file = "data/gene_cluster_genes/genes.fasta"
gff_folder = "results/gene_cluster_genes/mapping_output_RL28_expanded/annotations/"
intermediate_results_dir = "results/gene_cluster_genes/mapping_output_RL28_expanded/lovis4u/temp/"
output_gff_folder = "results/gene_cluster_genes/mapping_output_RL28_expanded/annotations/edited/"
insertion_positions_file = "results/gene_cluster_genes/mapping_output_RL28_expanded/insertion_positions_summary.tsv"
gene_lengths = {
    "GT1": 513,
    "GT2": 228,
    "GT3": 328,
    "GTA": 891,
    "GTB": 1119,
    "Flippase": 1421,
    "UDP-galactopyranosemutase": 1116,
}
gene_dbxref = {
    "GT1": "Dbxref=RefSeq:WP_005727395.1,SO:0001217,UniParc:UPI0001EC2972,UniRef:UniRef100_K1MUX3,UniRef:UniRef50_K1MUX3,UniRef:UniRef90_K1MUX3",
    "GT2": "",
    "GT3": "Dbxref=SO:0001217,UniParc:UPI000280BB18,UniRef:UniRef100_C7XKG0,UniRef:UniRef50_C7XKG0,UniRef:UniRef90_C7XKG0",
    "GTA": "Dbxref=RefSeq:WP_005723850.1,SO:0001217,UniParc:UPI0001B29BC6,UniRef:UniRef100_K1NSM6,UniRef:UniRef50_K1NSM6,UniRef:UniRef90_K1NSM6",
    "GTB": "Dbxref=RefSeq:WP_005723852.1,SO:0001217,UniParc:UPI0001B29BC7,UniRef:UniRef100_K1MPT4,UniRef:UniRef50_K1MPT4,UniRef:UniRef90_K1MPT4",
    "Flippase": "Dbxref=RefSeq:WP_005727389.1,SO:0001217,UniParc:UPI0001EC296D,UniRef:UniRef100_K1M3B3,UniRef:UniRef50_Q74JS1,UniRef:UniRef90_K1M3B3",
    "UDP-galactopyranosemutase": "Dbxref=EC:5.4.99.9,GO:0008767,GO:0009273,RefSeq:WP_005723855.1,SO:0001217,UniParc:UPI0001B29BCA,UniRef:UniRef100_E3R249,UniRef:UniRef50_A0A6I4Q4C7,UniRef:UniRef90_A8YX78",
}

# Ensure directories exist
os.makedirs(intermediate_results_dir, exist_ok=True)
os.makedirs(output_gff_folder, exist_ok=True)


# Create BLAST database for each FASTA file
def make_blast_db(fasta_file, intermediate_results_dir):
    db_name = os.path.join(intermediate_results_dir, os.path.basename(fasta_file) + "_db")
    command = f"makeblastdb -in {fasta_file} -dbtype nucl -out {db_name}"
    subprocess.run(command, shell=True, check=True)
    return db_name


# Run BLAST for a custom gene against a database
def run_blast(query_file, db_name, output_file):
    command = (
        f"blastn -query {query_file} -db {db_name} -out {output_file} "
        f"-outfmt '6 qseqid sseqid sstart send length pident' -max_hsps 1 -max_target_seqs 1"
    )
    subprocess.run(command, shell=True, check=True)


# Parse BLAST results
def parse_blast_results(blast_output_file, gene_lengths):
    results = []
    with open(blast_output_file, "r") as file:
        for line in file:
            qseqid, sseqid, sstart, send, length, pident = line.strip().split("\t")
            sstart, send, length = int(sstart), int(send), int(length)
            expected_length = gene_lengths[qseqid]

            # Adjust start and end positions
            if sstart > send:  # end is smaller, so we switch the coordinates
                sstart, send = send, sstart



            results.append({
                "qseqid": qseqid,
                "sseqid": sseqid,
                "sstart": sstart,
                "send": send,
                "length": length,
                "strand": "+" if sstart < send else "-",
                "pident": float(pident),
            })
    return results


# Modify GFF file with new BLAST hit features
def modify_gff_file(gff_file, blast_results, output_gff_folder):
    with open(gff_file, "r") as file:
        gff_data = file.readlines()

    # Identify the position where the ##FASTA section begins
    fasta_index = next((i for i, line in enumerate(gff_data) if line.startswith("##FASTA")), len(gff_data))

    # Extract feature lines before the ##FASTA section
    file_header = [line for line in gff_data[:fasta_index] if line.startswith("#")]
    feature_lines = [line for line in gff_data[:fasta_index] if not line.startswith("#")]

    # Extract the last feature line before ##FASTA for ID calculation
    last_feature_line = feature_lines[-1] if feature_lines else None
    if last_feature_line:
        last_id_match = re.search(r"ID=([\w_]+)", last_feature_line)
        if last_id_match:
            base_id = last_id_match.group(1)
            prefix, num = re.match(r"([A-Za-z_]+)(\d+)", base_id).groups()
            num = int(num)
        else:
            raise ValueError("No valid ID found in the last feature line.")
    else:
        raise ValueError("No feature lines found in the GFF file.")

    # Filter out overlapping features and prepare new features
    new_features = []
    non_overlapping_features = []

    for feature_line in feature_lines:
        parts = feature_line.split("\t")
        if len(parts) < 9:  # Skip malformed lines
            continue
        feature_start = int(parts[3])
        feature_end = int(parts[4])

        # Check for overlaps with BLAST results
        overlaps = any(
            max(result["sstart"], feature_start) <= min(result["send"], feature_end)
            for result in blast_results
        )

        if not overlaps:
            non_overlapping_features.append(feature_line)

    # Add new features based on BLAST results
    for i, result in enumerate(blast_results):
        new_id = f"{prefix}{num + (i + 1) * 5}"  # Increment ID by 5 for each feature
        attributes = f"ID={new_id};Name={result['qseqid']};locus_tag={new_id};product={result['qseqid']};{gene_dbxref.get(result['qseqid'], '')}"
        new_feature = "\t".join([
            "contig_1",  # seqid (chromosome/contig ID)
            "BLAST",  # source
            "CDS",  # type
            str(result["sstart"]),  # start
            str(result["send"]),  # end
            ".",  # score
            result["strand"],  # strand
            "0",  # phase
            attributes  # attributes
        ]) + "\n"
        new_features.append(new_feature)

    # Merge non-overlapping features and new features
    updated_gff_data = file_header + non_overlapping_features + new_features + gff_data[fasta_index:]

    # Save modified GFF file
    output_gff_file = os.path.join(output_gff_folder, os.path.basename(gff_file))
    with open(output_gff_file, "w") as file:
        file.writelines(updated_gff_data)
    print(f"Modified GFF file saved to {output_gff_file}")


# Add features based on insertion positions
def add_insertion_positions_to_gff(insertion_file, gff_folder, output_gff_folder):
    # Read insertion_positions_summary.tsv
    insertion_data = pd.read_csv(insertion_file, sep="\t")

    # Process each row in the insertion positions summary
    for _, row in insertion_data.iterrows():
        sample_name = row["Sample"]
        insertion_pos = row["Insertion_Position"]
        split_read_count = row["Split_Read_Count"]

        # Parse insertion position
        contig, start = insertion_pos.split(":")
        start = int(start)
        end = start + 5

        # Locate the corresponding GFF file for the sample
        gff_file = os.path.join(gff_folder, f"{sample_name}.gff3")
        if not os.path.exists(gff_file):
            print(f"GFF file for sample {sample_name} not found, skipping...")
            continue

        # Read the GFF file
        with open(gff_file, "r") as file:
            gff_data = file.readlines()

        # Identify the position where the ##FASTA section begins
        fasta_index = next((i for i, line in enumerate(gff_data) if line.startswith("##FASTA")), len(gff_data))

        # Extract feature lines before the ##FASTA section
        feature_lines = [line for line in gff_data[:fasta_index] if not line.startswith("#")]

        # Extract the last feature line before ##FASTA for ID calculation
        last_feature_line = feature_lines[-1] if feature_lines else None
        if last_feature_line:
            last_id_match = re.search(r"ID=([\w_]+)", last_feature_line)
            if last_id_match:
                base_id = last_id_match.group(1)
                prefix, num = re.match(r"([A-Za-z_]+)(\d+)", base_id).groups()
                num = int(num)
            else:
                raise ValueError("No valid ID found in the last feature line.")
        else:
            raise ValueError("No feature lines found in the GFF file.")

        # Create a new feature line
        new_id = f"{prefix}{num + 5}"  # Increment ID by 5 for each feature
        new_feature = "\t".join([
            "contig_1",  # Contig ID from insertion position
            "Insertion",  # Source
            "CDS",  # Type
            str(start),  # Start position
            str(end),  # End position
            ".",  # Score
            "+",  # Strand (default to "+" for simplicity)
            "0",  # Phase
            f"ID=ID={new_id};Name=Possible_gene_insertion;product={split_read_count}"  # Attributes
        ]) + "\n"

        # Insert the new feature before the ##FASTA section
        updated_gff_data = gff_data[:fasta_index] + [new_feature] + gff_data[fasta_index:]

        # Save the updated GFF file
        output_gff_file = os.path.join(output_gff_folder, os.path.basename(gff_file))
        with open(output_gff_file, "w") as file:
            file.writelines(updated_gff_data)
        print(f"Added insertion to GFF for sample {sample_name}, saved to {output_gff_file}")


# Main function
def gt_pipeline(fasta_folder, custom_genes_file, gff_folder, intermediate_results_dir, output_gff_folder, gene_lengths):
    # Iterate over FASTA files and corresponding GFF files
    for fasta_file in os.listdir(fasta_folder):
        if fasta_file.endswith(".fasta"):
            fasta_path = os.path.join(fasta_folder, fasta_file)
            sample_name = os.path.splitext(fasta_file)[0]

            gff_file = os.path.join(gff_folder, f"{sample_name}.gff3")
            if not os.path.exists(gff_file):
                print(f"No GFF file found for {sample_name}, skipping...")
                continue

            # Create BLAST database and run for each gt cluster gene
            db_name = make_blast_db(fasta_path, intermediate_results_dir)
            blast_output_file = os.path.join(intermediate_results_dir, f"{sample_name}_blast_results.tsv")
            run_blast(custom_genes_file, db_name, blast_output_file)

            # Parse BLAST results and modify GFF file with the BLAST hits
            blast_results = parse_blast_results(blast_output_file, gene_lengths)
            modify_gff_file(gff_file, blast_results, output_gff_folder)
            # add_insertion_positions_to_gff(insertion_positions_file, gff_folder, output_gff_folder)


# Run the pipeline
gt_pipeline(fasta_folder, custom_genes_file, gff_folder, intermediate_results_dir, output_gff_folder, gene_lengths)
