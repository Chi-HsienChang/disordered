import os
import re
import matplotlib.pyplot as plt
from collections import Counter
import pandas as pd

# Function to read the FASTA file and group sequences by clade
def read_fasta_grouped_by_clade(filename, name_to_clade):
    sequences_by_clade = {}
    current_clade = None
    current_sequence = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # If there is a current sequence, append it to the correct clade
                if current_sequence and current_clade:
                    sequences_by_clade.setdefault(current_clade, []).append(''.join(current_sequence))
                    current_sequence = []

                # Parse the header to extract the name and clade
                identifier = line[1:]  # Remove '>' from the start of the header
                parts = identifier.split('|')
                name = parts[0]  # Assume the first part of the header is the name
                current_clade = name_to_clade.get(name, 'Unknown')

            else:
                # Append sequence lines to the current sequence
                current_sequence.append(line)

        # Append the last sequence to the clade group
        if current_sequence and current_clade:
            sequences_by_clade.setdefault(current_clade, []).append(''.join(current_sequence))

    # Concatenate all sequences for each clade into a single string
    for clade in sequences_by_clade:
        sequences_by_clade[clade] = ''.join(sequences_by_clade[clade])

    return sequences_by_clade

# Function to count amino acid pairs in a sequence
def count_amino_acid_pairs(sequence):
    pairs = [sequence[i:i+2] for i in range(len(sequence) - 1)]
    pair_counts = Counter(pairs)
    return pair_counts

# Function to count amino acid pairs for each clade
def count_amino_acid_pairs_by_clade(sequences_by_clade):
    clade_pair_counts = {}

    for clade, sequence in sequences_by_clade.items():
        pair_counts = count_amino_acid_pairs(sequence)
        clade_pair_counts[clade] = pair_counts

    return clade_pair_counts

# Function to save amino acid pair proportions for each clade
def save_proportions_by_clade(clade_pair_counts, output_dir='output'):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for clade, pair_counts in clade_pair_counts.items():
        total_pairs = sum(pair_counts.values())
        output_filename = os.path.join(output_dir, f"{clade}_amino_acid_pair_proportions.txt")
        with open(output_filename, 'w') as file:
            for pair, count in sorted(pair_counts.items()):
                proportion = count * 100 / total_pairs
                file.write(f"{pair}: {proportion:.6f}\n")

    print(f"Proportions saved for each clade in the '{output_dir}' directory.")

# Function to plot and save figures for amino acid pairs of each clade
def plot_and_save_multiple_figures_by_clade(clade_pair_counts, output_dir='clade_plots', top_n=10):
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # For each clade, generate a plot and save it as a separate figure
    for clade, pair_counts in clade_pair_counts.items():
        total_pairs = sum(pair_counts.values())
        # Sort pairs by count and select the top N pairs
        top_pairs = pair_counts.most_common(top_n)
        pairs, counts = zip(*top_pairs)
        proportions = [count * 100 / total_pairs for count in counts]

        # Create a new figure for the current clade
        plt.figure(figsize=(10, 5))
        plt.bar(pairs, proportions, color='blue')
        plt.xlabel('Amino Acid Pairs')
        plt.ylabel('Proportion (%)')
        plt.title(f'Top {top_n} Amino Acid Pairs for Clade: {clade}')
        plt.xticks(rotation=90)
        plt.tight_layout()

        # Save the figure to a file named after the clade
        output_filename = os.path.join(output_dir, f"{clade}_amino_acid_pair_proportions.png")
        plt.savefig(output_filename)
        plt.close()  # Close the figure to avoid memory issues

    print(f"All figures have been saved in the '{output_dir}' directory.")

# Load the TSV file and create a mapping from TaxID to Named Lineage
def load_taxid_to_lineage_mapping(tsv_file_path):
    tsv_df = pd.read_csv(tsv_file_path, sep='\t')
    taxid_to_lineage = {str(row['Taxid']): row['Named Lineage'] for index, row in tsv_df.iterrows()}
    return taxid_to_lineage

# Map FASTA names to clades using TaxID
def map_fasta_names_to_clades(fasta_file_path, taxid_to_lineage):
    name_to_clade = {}
    with open(fasta_file_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                identifier = line[1:]  # Remove '>' from the start of the header
                parts = identifier.split('|')
                name = parts[0]  # Assume the first part of the header is the name
                taxID = parts[-1].split(':')[-1]  # Extract taxID from the last part
                named_lineage = taxid_to_lineage.get(taxID, 'NA,NA,NA,NA,NA')
                name_to_clade[name] = named_lineage.split(',')[3]  # Map to clade
    return name_to_clade




# File paths
fasta_file_path = 'aligned_real.fasta'
tsv_file_path = 'LUC7p_species_phylogeny.tsv'

# Read the TSV file containing taxID and Named Lineage
tsv_df = pd.read_csv(tsv_file_path, sep='\t')
# Create a dictionary to map taxID to Named Lineage
taxid_to_lineage = {str(row['Taxid']): row['Named Lineage'] for index, row in tsv_df.iterrows()}

# Initialize a dictionary to hold the mapping from FASTA name to Named Lineage
name_to_clade = {}
name_to_taxID = {}

# Process the FASTA file
with open(fasta_file_path, 'r') as fasta_file:
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            # Extract the identifier and remove the leading '>'
            identifier = line[1:]  # Remove '>' from the start of the header
            parts = identifier.split('|')
            name = parts[0]  # Assume the first part of the header is the name, now cleaned of '>'
            taxID = parts[-1].split(':')[-1]  # Extract taxID from the last part
            # Retrieve the Named Lineage using taxID
            named_lineage = taxid_to_lineage.get(taxID, 'NA,NA,NA,NA,NA')
            # Map the name to its corresponding Named Lineage
            name_to_clade[name] = named_lineage.split(',')[3]
            # print(named_lineage.split(',')[3])
            name_to_taxID[name] = taxID


# Read the sequences grouped by clade
sequences_by_clade = read_fasta_grouped_by_clade(fasta_file_path, name_to_clade)

# Count amino acid pairs for each clade
clade_pair_counts = count_amino_acid_pairs_by_clade(sequences_by_clade)

# Save the proportions to text files
save_proportions_by_clade(clade_pair_counts, output_dir='amino_acid_proportions')

# Generate and save the figures for each clade
plot_and_save_multiple_figures_by_clade(clade_pair_counts, output_dir='clade_plots', top_n=10)


