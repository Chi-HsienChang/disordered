from collections import Counter
import re


# Function to read the FASTA file and return the sequences as a single string
def read_fasta(filename):
    with open(filename, 'r') as file:
        content = file.read()
    # Extract sequences and concatenate them
    sequences = re.sub(r'^>.*$', '', content, flags=re.M).replace('\n', '')
    return sequences

# Function to count all possible consecutive amino acid pairs
def count_amino_acid_pairs(sequence):
    # Generate consecutive amino acid pairs
    pairs = [sequence[i:i+2] for i in range(len(sequence) - 1)]
    # Count the pairs
    pair_counts = Counter(pairs)
    return pair_counts

# Function to calculate proportions and save to a file
def save_proportions_to_file(pair_counts, output_filename):
    total_pairs = sum(pair_counts.values())
    with open(output_filename, 'w') as file:
        for pair, count in sorted(pair_counts.items()):
            proportion = count*100 / total_pairs
            file.write(f"{pair}: {proportion:.6f}\n")



import pandas as pd

# File paths
fasta_file_path = 'aligned_real.fasta'
tsv_file_path = 'LUC7p_species_phylogeny.tsv'

# Read the TSV file containing taxID and Named Lineage
tsv_df = pd.read_csv(tsv_file_path, sep='\t')
# Create a dictionary to map taxID to Named Lineage
taxid_to_lineage = {str(row['Taxid']): row['Named Lineage'] for index, row in tsv_df.iterrows()}

# Initialize a dictionary to hold the mapping from FASTA name to Named Lineage
name_to_clade = {}
# name_to_taxID = {}

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
            # name_to_taxID[name] = taxID

# Read the FASTA file
sequence = read_fasta('aligned_real.fasta')
# Count the pairs
pair_counts = count_amino_acid_pairs(sequence)

# Save the proportions to a file
save_proportions_to_file(pair_counts, 'amino_acid_pair_proportions.txt')

print("Proportions have been saved to 'amino_acid_pair_proportions.txt'.")
