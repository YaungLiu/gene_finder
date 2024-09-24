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

def find_genes(sequence):
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
    if len(sys.argv) != 2:
        print("Usage: python gene_finder.py <fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    sequence = read_fasta(fasta_file)
    genes = find_genes(sequence)

    if genes:
        print("Found genes:")
        for gene in genes:
            print(gene)
    else:
        print("No genes found.")

if __name__ == "__main__":
    main()

