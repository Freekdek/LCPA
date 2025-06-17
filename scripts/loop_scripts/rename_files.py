import os
import shutil

# directories
input_dir = 'results/genome_comparison/bella_ordered/to_be_renamed/'
output_dir = 'results/genome_comparison/bella_ordered/'
sample_list = 'results/sample_names_vvv_incomplete.txt'

# variables
file_list = sorted(os.listdir(input_dir), key=len)

# loop over the files and rename them
with open(sample_list) as f:
    for file in file_list:
        file_name = f.readline().split()[0]
        if file != file_name:
            print(file, file_name)
        shutil.copyfile(f'{input_dir}{file}', f'{output_dir}{file_name}.gbff')