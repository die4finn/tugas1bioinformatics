rna_to_amino_acid = {
    'AUG': 'M', 'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'UAU': 'Y',
    'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop', 'UGU': 'C', 'UGC': 'C',
    'UGA': 'Stop', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L',
    'CUG': 'L', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R',
    'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AUU': 'I', 'AUC': 'I',
    'AUA': 'I', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGU': 'S',
    'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V', 'GUC': 'V',
    'GUA': 'V', 'GUG': 'V', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A',
    'GCG': 'A', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def dna_to_rna(dna):
    return dna.replace('T', 'U')

def complement_dna_to_rna(dna):
    complement = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
    rna1 = ''.join(complement[base] for base in dna)
    return rna1

def rna_to_amino_acids(rna):
    amino_acids = []
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        if len(codon) == 3:  # Ensure codon is complete
            amino_acid = rna_to_amino_acid.get(codon, '')
            if amino_acid == 'Stop':
                break
            amino_acids.append(amino_acid)
    return ''.join(amino_acids)

dna_input = input("Enter DNA sequence: ").upper()
rna1 = complement_dna_to_rna(dna_input)
rna2 = dna_to_rna(dna_input)
amino_acids = rna_to_amino_acids(rna1)

print("DNA: ", dna_input)
print("RNA1: ", rna1)
print("RNA2: ", rna2)
print("AMINO ACIDS: ", amino_acids)
