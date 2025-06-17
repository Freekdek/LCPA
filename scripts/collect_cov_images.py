import os
import shutil
import glob


def copy_genome_coverage_images(input_folder, output_folder):
    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Find all matching files in the specified structure
    image_paths = glob.glob(
        os.path.join(input_folder, '*_report/images_qualimapReport/genome_coverage_across_reference.png'))

    for image_path in image_paths:
        # Extract the sample name from the parent folder of the report
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(image_path)))

        # Define a new name for each file based on its sample name
        new_image_name = f"{sample_name}_genome_coverage_across_reference.png"
        destination_path = os.path.join(output_folder, new_image_name)

        # Copy the file to the output folder
        shutil.copy(image_path, destination_path)
        print(f"Copied {image_path} to {destination_path}")


# Example usage
input_folder = 'results/gene_cluster_genes/mapping_output_RL28/qualimap_reports'
output_folder = 'results/gene_cluster_genes/mapping_output_RL28/coverage_images/'
copy_genome_coverage_images(input_folder, output_folder)
