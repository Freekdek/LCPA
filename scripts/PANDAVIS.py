"""
Made by Freek de Kreek, 2024
PANDAVIS --> PAngenome DAta VISualizer
This script takes a PIRATE.gene_families.ordered.tsv as input and an EggNOG database e5
(version 5.0.0, it will probably? work with other versions) as input.

It creates a general overview of the pangenome by producing multiple plots:
- Average gene length and average gene length by gene type (pangenome part: core, soft core, shell and cloud)*
- Pie chart of parts of the pangenome
- Gene counts per genomes, stacked bar chart visualizing the number of genes by gene type*
- Stacked bar chart of functional categories (COG) by gene type*
- Heatmap of functional group counts within individuals and all genomes
- Clustermap of gene differences within individuals and all genomes - requires gene-dists or providing a distance matrix

The figures are meant to visualize the overall genomes in the pangenome analysis, but also the variation withing genomes
from the same host (within host comparison).


* Gene type or pangenome part:
Core = genes present in > 95% of all genomes
Soft core =  genes present in 90 - 95% of all genomes
Shell =  genes present in 10 - 90% of all genomes
Cloud = genes present in < 10% of all genomes
"""

import os.path
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import scipy
import pandas as pd
from rapidfuzz import process, fuzz
from joblib import Parallel, delayed
import openpyxl
import argparse
import Bio from SeqIO

sns.set_theme(style="whitegrid")

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", help="Provide input file, a PIRATE.gene_families.ordered.tsv, with -i or --input")
parser.add_argument("--output", "-o", help="Provide output directory with -o or --output")
parser.add_argument("--database", "-db", help="Provide EggNOG database file -db or --database")
parser.add_argument("--func_only", "-f", help="Functional only when provided -f or --func_only")
parser.add_argument("--dist_matrix", "-d", help="Distance matrix of pangenome absence presence multiple sequence alignment")
parser.add_argument("--snp_dist", "-snp", help="Pangenome absence presence multiple sequence alignment, require snp-dists to be installed")
args = parser.parse_args()

gene_data = args.input
output_dir = args.output
output_dir_fig = os.path.join(output_dir, "figures")
eggNOG_db = args.database
func_only = args.func_only

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    os.makedirs(output_dir_fig)

gene_data = pd.read_csv("../results/pangenome_results/stringent/PIRATE.gene_families.ordered.tsv", sep='\t')

# Set gene_type, based on total dataset size of 53.
gene_data.loc[gene_data.number_genomes > 50, "gene_type"] = "core" # edit
gene_data.loc[(gene_data.number_genomes > 47) & (gene_data.number_genomes < 51), "gene_type"] = "soft_core" # edit
gene_data.loc[(gene_data.number_genomes > 5) & (gene_data.number_genomes < 48), "gene_type"] = "shell" # edit
gene_data.loc[gene_data.number_genomes < 6, "gene_type"] = "cloud" # edit

########### Histograms of average lengths #############
sns.set_theme(style="whitegrid")
gene_data.hist("average_length(bp)", color="skyblue")
gene_data.hist("average_length(bp)", by="gene_type", color="skyblue")

############# Pie chart pangenome components ###########
# Data for the gene families
labels = ['Core genes', 'Soft-core genes', 'Shell genes', 'Cloud genes']
sizes = [1597, 57, 821, 981] # edit
colors = sns.color_palette("pastel")[0:4]  # Seaborn pastel color palette

# Explode parameter to add spacing between the slices
explode = (0.05, 0.12, 0.05, 0.05)  # Fraction to "explode" each slice


# Function to display both percentage and actual values
def func(pct, allvals):
    absolute = int(pct / 100. * sum(allvals))
    return f"{pct:.1f}%\n({absolute:d})"


# Plotting the pie chart with explode to separate slices
plt.figure(figsize=(7, 7))
plt.pie(sizes, labels=labels, colors=colors, explode=explode, autopct=lambda pct: func(pct, sizes), startangle=90,
        wedgeprops={"ls": ""}, counterclock=False)
plt.title('Parts of the pangenome')
plt.axis('equal')
plt.savefig(f"{output_dir_fig}Pangenome_pie_chart")

################### Stacked bar chart gene counts by pangenome part ####################
# Make a subselection based on the pangenome part/type, called gene_type
core_genes = gene_data[gene_data.gene_type == "core"]
soft_core_genes = gene_data[gene_data.gene_type == "soft_core"]
shell_genes = gene_data[gene_data.gene_type == "shell"]
cloud_genes = gene_data[gene_data.gene_type == "cloud"]


# Function to count non-empty samples in the last 53 columns for a given DataFrame
def count_non_empty_samples(df):
    sample_columns = df.iloc[:, -54:-1] # edit
    return sample_columns.notna().sum()


# Count non-empty samples for each subselection
core_counts = count_non_empty_samples(core_genes)
soft_core_counts = count_non_empty_samples(soft_core_genes)
shell_counts = count_non_empty_samples(shell_genes)
cloud_counts = count_non_empty_samples(cloud_genes)

# Combine results into a DataFrame for plotting
results = pd.DataFrame({
    'Samples': core_counts.index,
    'Core': core_counts.values,
    'Soft core': soft_core_counts.values,
    'Shell': shell_counts.values,
    'Cloud': cloud_counts.values
})

results.set_index('Samples', inplace=True)

# Plotting the stacked bar chart
plt.figure(figsize=(14, 8))
results.plot(kind='bar', stacked=True, figsize=(14, 8), color=sns.color_palette("husl", 4))
plt.title('Gene counts per genome (by pangenome part)')
plt.xlabel('Samples')
plt.ylabel('Gene count')
plt.xticks(rotation=90, ha='right')  # Rotate the sample names for better readability
plt.legend(title='Gene type', loc='lower left')
plt.tight_layout()
plt.savefig(f"{output_dir_fig}Gene_counts.png")


############### COG category classification ##################
# Matching to the EggNOG e5 database with exact matching, fuzzy matching and based on key terms to assign functional groups
threshold_n = 80

# Key terms for each functional category
functional_keywords = {
    # INFORMATION STORAGE AND PROCESSING
    'J': ['ribosomal', 'ribosome', 'biogenesis', 'translation', 'tRNA', 'rRNA', 'protein synthesis',
          'peptidyl transferase', 'initiation', 'elongation', 'termination'],
    'A': ['RNA processing', 'splicing', 'modification', 'editing', 'maturation', 'ribonuclease', 'snoRNA',
          'RNA binding', 'polyadenylation', '5’ capping'],
    'K': ['transcription', 'RNA polymerase', 'transcription factor', 'promoter binding', 'enhancer', 'silencer',
          'gene regulation', 'RNA synthesis'],
    'L': ['DNA replication', 'recombination', 'repair', 'helicase', 'ligase', 'topoisomerase', 'polymerase',
          'DNA damage response', 'replication fork', 'checkpoint', 'tranposase'],
    'B': ['chromatin', 'histone', 'nucleosome', 'DNA packaging', 'epigenetic', 'chromatin remodeling',
          'DNA methylation', 'transcriptional regulation', 'nuclear structure', 'helix-turn-helix', 'DNA-binding'],

    # CELLULAR PROCESSES AND SIGNALING
    'D': ['cell cycle', 'division', 'mitosis', 'meiosis', 'partitioning', 'septum formation', 'checkpoints',
          'cytokinesis', 'cyclin', 'CDK'],
    'Y': ['nuclear envelope', 'nucleolus', 'nuclear pore', 'nuclear transport', 'nuclear import', 'nuclear export'],
    'V': ['defense', 'toxin', 'immunity', 'antimicrobial', 'restriction-modification', 'anti-toxin', 'resistance',
          'virulence factors', 'immune response', 'bacteriocin'],
    'T': ['signal transduction', 'kinase', 'phosphatase', 'receptor', 'second messenger', 'cAMP', 'GTPase', 'sensor',
          'signal cascade', 'feedback regulation', 'signal peptide'],
    'M': ['cell wall', 'membrane', 'envelope', 'peptidoglycan', 'lipopolysaccharide', 'outer membrane',
          'membrane potential', 'transport proteins', 'membrane fusion', 'surface'],
    'N': ['flagella', 'pili', 'motility', 'chemotaxis', 'movement', 'locomotion', 'twitching', 'swimming', 'adhesion'],
    'Z': ['cytoskeleton', 'microtubule', 'actin', 'filament', 'tubulin', 'structural support', 'cell shape',
          'cell motility'],
    'W': ['extracellular', 'biofilm', 'matrix', 'capsule', 'adhesin', 'quorum sensing', 'extracellular matrix'],
    'U': ['secretion', 'vesicle', 'trafficking', 'exocytosis', 'endocytosis', 'transport', 'membrane transport',
          'vesicular transport', 'cellular export'],
    'O': ['chaperone', 'protein folding', 'ubiquitin', 'degradation', 'modification', 'proteasome',
          'protein quality control', 'post-translational modification'],

    # METABOLISM
    'C': ['ATP', 'oxidation', 'respiration', 'glycolysis', 'fermentation', 'electron transport', 'photosynthesis',
          'metabolic pathway', 'energy metabolism'],
    'G': ['sugar', 'carbohydrate', 'glycolysis', 'polysaccharide', 'monosaccharide', 'glucose', 'fructose', 'transport',
          'glycogen', 'starch'],
    'E': ['amino acid', 'protein degradation', 'ammonia', 'transaminase', 'deaminase', 'glutamine', 'lysine',
          'proteolysis', 'nitrogen metabolism'],
    'F': ['nucleotide', 'nucleoside', 'purine', 'pyrimidine', 'DNA synthesis', 'RNA synthesis', 'nucleotide metabolism',
          'deoxyribonucleotide'],
    'H': ['coenzyme', 'vitamin', 'cofactor', 'biotin', 'NAD', 'FAD', 'vitamin B', 'folate', 'enzyme activation',
          'cofactor binding', 'phosphate oxidase'],
    'I': ['lipid', 'fatty acid', 'phospholipid', 'steroid', 'cholesterol', 'beta-oxidation', 'membrane lipid',
          'lipid metabolism', 'triglyceride', 'lipoprotein', 'hydrolase'],
    'P': ['ion', 'transport', 'metal ion', 'sodium', 'potassium', 'calcium', 'iron', 'magnesium', 'zinc', 'ion channel',
          'ion transport', 'reductase'],
    'Q': ['secondary metabolite', 'antibiotic', 'toxin', 'pigment', 'biosynthesis', 'natural products',
          'plant metabolite', 'terpenes'],

    # POORLY CHARACTERIZED
    'R': ['general function'],
    'S': ['unknown', 'hypothetical protein', 'DUF', 'uncharacterized', 'novel', 'molecular function',
          'experimental validation', 'hypothetical', 'predicted protein', 'conserved domain', 'uncharacterized protein',
          'unknown mechanism'],

    # PHAGE PROTEINS
    'X': ['phage', 'prophage', 'head', 'capsid', 'tail', 'integrase', 'viral', 'bacteriophage', 'lytic', 'lysogenic',
          'viral replication', 'viral proteins'],

    # MOBILE GENETIC ELEMENTS -- Transposases
    'MGE': ['transposase', 'transposon', 'mobile genetic element', 'plasmid']
}


# Function to assign functional category based on key terms
def assign_by_key_terms(description):
    for category, keywords in functional_keywords.items():
        for keyword in keywords:
            if keyword.lower() in description.lower():
                return category
    return None


# Function to match gene descriptions with fuzzy matching if no exact match
def fuzzy_match(row, right_data, right_column, threshold=threshold_n):
    match_result = process.extractOne(row, right_data[right_column].values, scorer=fuzz.partial_ratio)
    if match_result is None:
        return None
    best_match, score, _ = match_result
    if score >= threshold:
        return best_match
    else:
        return None


# Function to prioritize exact match, then fuzzy match, and finally key terms
def match_functional_category(gene_data, func_annot_data, threshold=threshold_n):
    # Step 1: Exact match and hypothetical, domain of unknown function and transposases classification
    merged_df = gene_data.merge(func_annot_data, how='left', left_on='consensus_product', right_on='gene_description')
    merged_df.loc[merged_df['consensus_product'] == 'hypothetical protein', 'functional_category'] = 'S'
    merged_df.loc[merged_df['consensus_product'] == 'hypothetical protein', 'gene_description'] = 'hypothetical protein'
    merged_df.loc[merged_df['consensus_product'].str.contains('DUF'), 'functional_category'] = 'S'
    merged_df.loc[merged_df['consensus_product'].str.contains('DUF'), 'gene_description'] = 'hypothetical protein'
    merged_df.loc[merged_df['consensus_product'].str.contains('transposase', case=False), 'functional_category'] = 'MGE'
    merged_df.loc[
        merged_df['consensus_product'].str.contains('transposase', case=False), 'gene_description'] = 'transposase'
    print("After exact matching:")
    print("NA values:", merged_df["functional_category"].isna().sum())
    print("Number of unknown functions:", merged_df["functional_category"].value_counts()["S"], "\n")

    # Step 2: Apply fuzzy matching where no exact match was found
    no_match_indices = merged_df[merged_df['functional_category'].isna()].index
    fuzzy_matches = parallel_fuzzy_merge(gene_data.loc[no_match_indices], func_annot_data,
                                         left_on='consensus_product', right_on='gene_description', threshold=threshold)

    # Update the dataframe with fuzzy matches
    merged_df.loc[no_match_indices, 'fuzzy_match'] = fuzzy_matches
    fuzzy_non_empty = merged_df[merged_df['fuzzy_match'].notna()]
    fuzzy_merged_df = fuzzy_non_empty.merge(func_annot_data, how='left', left_on='fuzzy_match',
                                            right_on='gene_description', suffixes=('', '_fuzzy'))

    # Fill in missing functional categories from fuzzy matching
    merged_df['functional_category'] = merged_df['functional_category'].fillna(
        fuzzy_merged_df['functional_category_fuzzy'])
    merged_df['gene_description'] = merged_df['gene_description'].fillna(fuzzy_merged_df['gene_description_fuzzy'])

    print("After fuzzy matching:")
    print("NA values:", merged_df["functional_category"].isna().sum())
    print("Number of unknown functions:", merged_df["functional_category"].value_counts()["S"], "\n")

    # Step 3: Apply key term matching for remaining unmatched rows
    unknown_func_indices = merged_df[merged_df['functional_category'] == 'S'].index
    for idx in unknown_func_indices:
        key_term_category = assign_by_key_terms(gene_data.loc[idx, 'consensus_product'])
        if key_term_category is not None:
            merged_df.loc[idx, 'functional_category'] = key_term_category  # Only update if a valid term match is found
    still_no_match_indices = merged_df[merged_df['functional_category'].isna()].index
    merged_df.loc[still_no_match_indices, 'functional_category'] = gene_data.loc[
        still_no_match_indices, 'consensus_product'].apply(assign_by_key_terms)
    merged_df['functional_category'].fillna('S')
    merged_df.loc[merged_df['functional_category'].isna(), 'functional_category'] = 'S'

    # Remove auxiliary columns
    merged_df.drop(columns=['fuzzy_match', 'gene_description_fuzzy', 'functional_category_fuzzy'], errors='ignore',
                   inplace=True)

    return merged_df


# Parallel fuzzy matching to improve performance
def parallel_fuzzy_merge(left_data, right_data, left_on, right_on, threshold=threshold_n, n_jobs=-1):
    results = Parallel(n_jobs=n_jobs)(
        delayed(fuzzy_match)(row, right_data, right_on, threshold) for row in left_data[left_on]
    )
    return results


# Load data
func_annot_data = pd.read_csv(eggNOG_db), sep="\t",
                              names=["OMA_group", "eggNOG_OG", "functional_category", "gene_description"])
func_annot_data['functional_category'] = func_annot_data['functional_category'].str[0]

# Some duplicates are present, and for simplification only one functional group is assigned (based on majority)
category_counts = func_annot_data.groupby(['gene_description', 'functional_category']).size().reset_index(name='count')
majority_category = category_counts.loc[category_counts.groupby('gene_description')['count'].idxmax()]

# Start matching the functional groups
merged_gene_data = match_functional_category(gene_data, majority_category)

# Check the number of NaN values after merging
print("After keyword matching:")
print("NA values:", merged_gene_data["functional_category"].isna().sum())
print("Number of unknown functions:", merged_gene_data["functional_category"].value_counts()["S"])

################### Stacked bar chart of functional categories by gene type ####################
# Map one-letter functional category codes to orthologous cluster names for better readability
category_map = {
    'J': 'Translation, ribosomal structure and biogenesis',
    'A': 'RNA processing and modification',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'B': 'Chromatin structure and dynamics',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'Y': 'Nuclear structure',
    'V': 'Defense mechanisms',
    'T': 'Signal transduction mechanisms',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'Z': 'Cytoskeleton',
    'W': 'Extracellular structures',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'O': 'Posttranslational modification, protein turnover, chaperones',
    'C': 'Energy production and conversion',
    'G': 'Carbohydrate transport and metabolism',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
    'R': 'General function prediction only',
    'S': 'Function unknown',
    'X': 'Phage genes',
    'MGE': 'Transposases',
    'Hypothetical protein': 'Hypothetical protein'
}

# Replace one-letter functional category codes with orthologous cluster names
merged_gene_data.loc[merged_gene_data['gene_description'] == 'hypothetical protein', 'functional_category'] = 'Hypothetical protein'
merged_gene_data['functional_category'] = merged_gene_data['functional_category'].map(category_map)

# Gene type in order
gene_type_order = ['core', 'soft_core', 'shell', 'cloud']
merged_gene_data['gene_type'] = pd.Categorical(merged_gene_data['gene_type'], categories=gene_type_order, ordered=True)

# Pivot table
pivot_table = pd.pivot_table(
    merged_gene_data,
    index='functional_category',
    columns='gene_type',
    aggfunc='size',
    fill_value=0,
    observed=True
)

# Sort the pivot table by total gene count across all gene types (sum along rows)
pivot_table['Total'] = pivot_table.sum(axis=1)
pivot_table = pivot_table.sort_values(by='Total', ascending=True)
pivot_table = pivot_table.drop(columns='Total')  # Remove the helper 'Total' column

# Calculate the total gene count per functional category and percentage
grand_total_genes = pivot_table.sum(axis=1).sum()
if func_only:
    pivot_table = pivot_table.drop(index=['Hypothetical protein', 'Function unknown'])
total_genes_per_category = pivot_table.sum(axis=1)
percentages = (total_genes_per_category / grand_total_genes) * 100

# Plot the stacked bar chart with the sorted functional categories and gene_type order
husl_colors = sns.color_palette("husl", 4)
husl_cmap = ListedColormap(husl_colors)
ax = pivot_table.plot(kind='barh', stacked=True, figsize=(12, 8), colormap=husl_cmap)

# Add percentage labels next to each functional category
for i, (total, percentage) in enumerate(zip(total_genes_per_category, percentages)):
    ax.text(total + 2, i, f'{percentage:.1f}%', va='center', fontsize=10, color='black')

plt.title("Functional categories (COG) by gene type", fontsize=16)
plt.xlabel("Number of genes", fontsize=12)
plt.ylabel("Functional category (COG)", fontsize=12)
plt.xticks(rotation=0, ha='right')  # Rotate x-axis labels for better readability
plt.legend(title="Gene type", title_fontsize='13', fontsize='11')
plt.tight_layout()

# Save files, fig and table, with or without functional groups only
plt.savefig(f"{output_dir_fig}functional_categories_by_gene_type.png")
merged_gene_data.to_excel(f"{output_dir}COG_matched_gene_families.xlsx", index=False)
if func_only:
    merged_gene_data = merged_gene_data[(merged_gene_data['functional_category'] != 'Hypothetical protein') & (merged_gene_data['functional_category'] != 'Function unknown')]
    merged_gene_data.to_excel(f"{output_dir}COG_matched_gene_families_func_only.xlsx", index=False)

################ Heatmap functional groups ######################
# Group by functional category and sum the number of genes across samples
sample_columns = merged_gene_data.columns[merged_gene_data.columns.get_loc('Bella1_1'):merged_gene_data.columns.get_loc('Tonna1_2') + 1] # edit
functional_group_totals = merged_gene_data.groupby('functional_category')[sample_columns].count()

intra_individuals = ["Bella", "Campanarius", "Nefesh", "Rhea", "Rivka"] # edit

for indv in intra_individuals:
    plt.figure(figsize=(12, 8))
    sample_df = functional_group_totals.filter(regex=f"{indv}*")
    sample_df_sorted = sample_df.sort_values(by=sample_df.columns[0], ascending=False)
    sns.heatmap(sample_df_sorted, cmap="coolwarm", annot=True, fmt=".0f", linewidths=.5)
    plt.title("Functional group counts", fontsize=16)
    plt.xlabel("Samples", fontsize=12)
    plt.ylabel("Functional category", fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0, ha='right')
    plt.tight_layout()
    plt.savefig(f"{output_dir_fig}/functional_group_heatmap_{indv}.png")

# Heatmap for all samples
plt.figure(figsize=(22, 8))
sample_df = functional_group_totals
sample_df_sorted = sample_df.sort_values(by=sample_df.columns[0], ascending=False)
sns.heatmap(sample_df_sorted, cmap="coolwarm", annot=False, fmt=".0f", linewidths=.5)
plt.title("Functional group counts", fontsize=16)
plt.xlabel("Samples", fontsize=12)
plt.ylabel("Functional category", fontsize=12)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0, ha='right')
plt.tight_layout()
plt.savefig(f"{output_dir_fig}/functional_group_heatmap_all.png")

# Load the merged gene data
merged_gene_data = pd.read_excel("merged_gene_data_func_filtered_mge.xlsx")

# Define functional groups to exclude
exclude_groups = ['Hypothetical protein', 'Function unknown']

# Create a mask to filter out these functional categories
mask = ~merged_gene_data['functional_category'].isin(exclude_groups)

# Filter the merged_gene_data to only include the rows you want
filtered_gene_data = merged_gene_data[mask]

# Get the indices of the genes to keep (since merged_gene_data is ordered the same as the alignment)
indices_to_keep = filtered_gene_data.index.tolist()
print(len(indices_to_keep))

# Load the binary gene presence/absence FASTA file
input_fasta = args.pangenome_alignment
output_fasta = f"{output_dir}/filtered_binary_presence_absence.fasta"

# Use SeqIO to read the FASTA and write only the sequences at the indices to keep
with open(output_fasta, "w") as output_handle:
    for idx, record in enumerate(SeqIO.parse(input_fasta, "fasta")):
        record.seq = "".join([record.seq[i] for i in indices_to_keep])
        # record.seq = record.seq[indices_to_keep]
        SeqIO.write(record, output_handle, "fasta")

print(f"Filtered alignment written to {output_fasta}")

################### Plot clustermaps #############################
if args.dist_matrix not None:
    total_dist_matrix = pd.read_csv(args.dist_matrix)
elif args.snp_diff not None:
    subprocess.run(f"snp-dists {args.snp_diff} > {output_dir}pangenome_dist_matix.csv")
    total_dist_matrix - pd.read_csv(f"{output_dir}pangenome_dist_matrix.csv")
total_dist_matrix.index.name = ""

sns.set_theme(style="whitegrid")

for indv in intra_individuals:
    dist_matrix = total_dist_matrix.filter(regex=f"{indv}*", axis='index').filter(regex=f"{indv}*", axis='columns')

clustermap = sns.clustermap(dist_matrix, annot=True, cmap="mako", fmt="3")
clustermap.savefig(f"{output_dir_fig}/{indv}_clustermap.png")

# Total clustermap
sample_names = total_dist_matrix.index
unique_groups = set([name.split('_')[0] for name in sample_names])
palette = sns.color_palette("husl", len(unique_groups))
group_colors = {group: palette[i] for i, group in enumerate(unique_groups)}
row_colors = sample_names.map(lambda x: group_colors[x.split('_')[0]])
col_colors = row_colors.copy()

clustermap = sns.clustermap(total_dist_matrix, cmap="mako", figsize=(30,30),row_colors=row_colors, col_colors=col_colors, linewidths=.5)
colorbar = clustermap.cax
colorbar.tick_params(labelsize=30)
clustermap.savefig(f"{output_dir_fig}/total_clustermap.png")
