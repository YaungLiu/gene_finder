import sys

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        sequence = ''
        for line in file:
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence

def reverse_complement(sequence):
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
        'H': 'D', 'V': 'B', 'N': 'N'
    }
    return ''.join(complement.get(base, 'N') for base in reversed(sequence))

def has_rbs(sequence, start_index, rbs_sequence, rbs_distance):
    """Check if there's an RBS sequence within the specified distance upstream of the start codon."""
    upstream_start = max(0, start_index - rbs_distance)
    upstream_region = sequence[upstream_start:start_index]
    return rbs_sequence in upstream_region

def find_genes(sequence, min_length, rbs_sequence, rbs_distance):
    stop_codons = ['TAA', 'TAG', 'TGA']
    genes = []

    # Function to find genes in a given reading frame
    def search_in_frame(seq, frame):
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon == 'ATG':  # Start codon found
                for j in range(i+3, len(seq) - 2, 3):
                    stop_codon = seq[j:j+3]
                    if stop_codon in stop_codons:
                        gene = seq[i:j+3]
                        if len(gene) / 3 >= min_length:  # Only include if length is >= min_length codons
                            # Check for RBS within the upstream region
                            if has_rbs(seq, i, rbs_sequence, rbs_distance):
                                genes.append(gene)
                        break

    # Search in the forward sequence (3 frames)
    for frame in range(3):
        search_in_frame(sequence, frame)

    # Get reverse complement and search in reverse sequence (3 frames)
    rev_comp_sequence = reverse_complement(sequence)
    for frame in range(3):
        search_in_frame(rev_comp_sequence, frame)

    return genes

def main():
    if len(sys.argv) != 5:
        print("Usage: python gene_finder.py <fasta_file> <min_length> <rbs_sequence> <rbs_distance>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    min_length = int(sys.argv[2])  # Convert minimum length to an integer
    rbs_sequence = sys.argv[3]  # Ribosome binding site sequence (e.g., AGGAGG)
    rbs_distance = int(sys.argv[4])  # Maximum distance upstream of start codon for the RBS

    sequence = read_fasta(fasta_file)
    genes = find_genes(sequence, min_length, rbs_sequence, rbs_distance)

    if genes:
        print(f"Found genes longer than {min_length} codons with RBS ({rbs_sequence}) within {rbs_distance}bp upstream:")
        for gene in genes:
            print(gene)
    else:
        print(f"No genes found matching the criteria.")

if __name__ == "__main__":
    main()
