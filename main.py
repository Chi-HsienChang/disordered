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

# Read the FASTA file
sequence = read_fasta('aligned_real.fasta')

# Count the pairs
pair_counts = count_amino_acid_pairs(sequence)

# Save the proportions to a file
save_proportions_to_file(pair_counts, 'amino_acid_pair_proportions.txt')

print("Proportions have been saved to 'amino_acid_pair_proportions.txt'.")
