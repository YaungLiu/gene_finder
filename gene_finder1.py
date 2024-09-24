import sys

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        sequence = ''
        for line in file:
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence

def find_genes(sequence):
    stop_codons = ['TAA', 'TAG', 'TGA']
    genes = []

    # Check each of the three reading frames
    for frame in range(3):
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon == 'ATG':  # Start codon found
                for j in range(i+3, len(sequence) - 2, 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:
                        gene = sequence[i:j+3]
                        genes.append(gene)
                        break
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
