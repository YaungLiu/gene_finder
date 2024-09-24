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

def find_genes(sequence, min_length):
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
    if len(sys.argv) != 3:
        print("Usage: python gene_finder.py <fasta_file> <min_length>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    min_length = int(sys.argv[2])  # Convert minimum length to an integer
    sequence = read_fasta(fasta_file)
    genes = find_genes(sequence, min_length)

    if genes:
        print(f"Found genes longer than {min_length} codons:")
        for gene in genes:
            print(gene)
    else:
        print(f"No genes longer than {min_length} codons found.")

if __name__ == "__main__":
    main()
