################################################################################
####################### EXERCISES Laura Ciaran 192802 ##########################
################################################################################

class Sequence (object):

    molecular_weight = {}
    alphabet = ""

    def __init__(self, identifier, sequence):
        self.__identifier = identifier
        self.__sequence = sequence

    def check_alphabet(self):
        wrong_letters = []
        [wrong_letters.append(x) for x in self.__sequence if x not in self.alphabet]
        if len(wrong_letters) != 0:
            for x in wrong_letters:
                raise ValueError("Impossible to create instance: "+x+" not possible")

    def get_identifier(self):
        """Function that returns the id"""
        return self.__identifier

    def get_sequence(self):
        """Function that returns the sequence"""
        return self.__sequence

    def get_mw(self):
        """Function that returns the molecular weight of a sequence"""
        mw = sum(self.molecular_weight[aa] for aa in self.__sequence)
        return (mw)

    def has_subsequence(self, sequence_obj):
        """Function that returns a boolean depending on if a given subsequence is
        inside a sequence"""
        return (sequence_obj.get_sequence() in self.__sequence)

class ProteinSequence(Sequence):

    def __init__(self, identifier, sequence):
        Sequence.__init__(self, identifier, sequence)
        self.molecular_weight = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
        self.alphabet = 'ACDEFGHIKLMNPQRSTVWY'

class NucleotideSequence(Sequence):

    table = {}
    stop_codons = ""
    start_codons = ""

    def __init__(self, identifier, sequence):
        Sequence.__init__(self, identifier, sequence)

    def translate(self):
        """This function returns a protein sequence given a nucleotide sequence"""
        n=0
        prot_seq=""
        sequence = self.get_sequence()
        identifier = self.get_identifier()
        if len(sequence) % 3 == 0:
            for i in range(0, len(sequence), 3):
                codon = sequence[i:i+3]
                if codon in self.start_codons:
                    n = 1
                if n == 1:
                    if codon not in self.stop_codons:
                        prot_seq += self.table[codon]
                else:
                    break
        return ProteinSequence(identifier, prot_seq)

class RNASequence (NucleotideSequence):

    def __init__(self, identifier, sequence):
        NucleotideSequence.__init__(self, identifier, sequence)
        self.molecular_weight = {'A': 363.0, 'C': 339.0, 'U': 340.0, 'G': 379.0}
        self.alphabet = 'GAUC'
        self.start_codons=['UUG', 'CUG', 'AUG']
        self.stop_codons=['UAA', 'UAG', 'UGA']
        self.table={'GUC': 'V', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GUU': 'V', 'AAC': 'N', 'AGG': 'R', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'AGA': 'R', 'AAU': 'N', 'ACU': 'T', 'CAC': 'H', 'GUG': 'V', 'CCG': 'P', 'CCA': 'P', 'AGU': 'S', 'CCC': 'P', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'CGA': 'R', 'CAG': 'Q', 'CGC': 'R', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'CCU': 'P', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GAG': 'E', 'UCC': 'S', 'UAC': 'Y', 'CGU': 'R', 'GAA': 'E', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'UCA': 'S', 'AUG': 'M', 'CUG': 'L', 'AUU': 'I', 'CAU': 'H', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'GAC': 'D', 'GUA': 'V', 'UGC': 'C', 'GCU': 'A', 'UGU': 'C', 'CUC': 'L', 'UUG': 'L', 'UUA': 'L', 'GAU': 'D', 'UUC': 'F'}

    def reverse_transcribe(self):
        """Function that returns the DNA sequence given a RNA sequence"""
        identifier =self.get_identifier()
        sequence = self.get_sequence()
        return DNASequence (identifier, sequence.replace("U", "T"))

class DNASequence (NucleotideSequence):

    def __init__(self, identifier, sequence):
        NucleotideSequence.__init__(self, identifier, sequence)
        self.molecular_weight = {'A': 347.0, 'C': 323.0, 'T': 322.0, 'G': 363.0}
        self.alphabet = 'GATC'
        self.start_codons = ['TTG', 'CTG', 'ATG']
        self.stop_codons = ['TAA', 'TAG', 'TGA']
        self.table = {'CTT': 'L', 'ATG': 'M', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'AGC': 'S', 'AGA': 'R', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'ACT': 'T', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'TAC': 'Y', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GAC': 'D', 'GAA': 'E', 'AAG': 'K', 'AAA': 'K', 'AAC': 'N', 'CTC': 'L', 'CAT': 'H', 'AAT': 'N', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'TGT': 'C', 'TCT': 'S', 'GAT': 'D', 'TTT': 'F', 'TGC': 'C', 'TGG': 'W', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TCA': 'S', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A'}

    def transcribe(self):
        """Function that returns the RNA sequence given a DNA sequence"""
        identifier =self.get_identifier()
        sequence = self.get_sequence()
        return RNASequence (identifier, sequence.replace("T", "U"))


# Vicen=RNASequence("Laura", "ATG")
# Vicen.check_alphabet()
# id=Vicen.get_identifier()
# seq=Vicen.get_sequence()
# mw=Vicen.get_mw()
# Prot_seq = Vicen.translate()
# RNA_seq = Vicen.transcribe()
#
# print(id)
# print(seq)
# print(mw)
# print (RNA_seq)
# print (Prot_seq)
