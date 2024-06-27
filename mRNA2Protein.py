import sys


def translate_sequence(cdna_sequence):
    # Define the codon to amino acid mapping
    codon_table = {
        'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'ACU': 'T',
        'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    # Convert DNA sequence to RNA sequence
    rna_sequence = cdna_sequence.replace('T', 'U')

    # Initialize the protein sequence
    protein_sequence = ''

    # Iterate over the RNA sequence in codon chunks
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        if len(codon) == 3:
            amino_acid = codon_table.get(codon, '?')
            if amino_acid == '*':  # Stop codon
                break
            protein_sequence += amino_acid

    return protein_sequence

def translate_sequence_from_file(sequence_file):

    f = open(sequence_file, "r")
    f.readline()
    mRNA_sequence = f.readline().strip()
    protein_sequence = translate_sequence(mRNA_sequence)
    print(protein_sequence)

translate_sequence_from_file(sys.argv[1])