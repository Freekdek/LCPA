import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def combine_fasta(directory, output_file):
    """
    Combines all FASTA files in a directory into a single multi-FASTA file.
    Each FASTA file's contigs are merged into a single sequence.

    Parameters:
        directory (str): Path to the directory containing FASTA files.
        output_file (str): Path to the output multi-FASTA file.
    """
    combined_records = []

    # Iterate over all files in the directory
    for file_name in os.listdir(directory):
        if file_name.endswith(".fasta") or file_name.endswith(".fa"):
            file_path = os.path.join(directory, file_name)
            sample_name = os.path.splitext(file_name)[0]  # Extract sample name (file name without extension)

            # Concatenate all sequences from the FASTA file
            concatenated_sequence = ""
            for record in SeqIO.parse(file_path, "fasta"):
                concatenated_sequence += str(record.seq)

            # Create a single SeqRecord for the merged sequence
            combined_record = SeqRecord(
                Seq(concatenated_sequence),
                id=sample_name,  # Use the sample name as the FASTA header
                description=""  # Clear any description
            )
            combined_records.append(combined_record)

    # Write all combined records to the output multi-FASTA file
    with open(output_file, "w") as output_handle:
        SeqIO.write(combined_records, output_handle, "fasta")


# Example usage
input_directory = "/mnt/c/Users/freek/Documents/BioSb/Internship/PangenomeProject/results/assembly_results/stringent/filtered_genomes/"
output_fasta = "/mnt/c/Users/freek/Documents/BioSb/Internship/PangenomeProject/results/probiotic_properties/combined_filtered_genomes.fasta"
combine_fasta(input_directory, output_fasta)
print(f"Combined FASTA file saved to: {output_fasta}")
