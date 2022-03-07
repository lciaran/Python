################################################################################
####################### EXERCISES Laura Ciaran 192802 ##########################
################################################################################
import sys
import os

class Sequence(object):

    molecular_weight = {}
    alphabet = ""

    def __init__(self, identifier, sequence):
        self.__identifier = identifier
        self.__sequence = sequence

        for letter in self.__sequence:
            if letter not in self.alphabet:
                raise IncorrectSequenceLetter(letter, self.__class__.__name__)

    def get_identifier(self):
        """Function that returns the id"""
        return self.__identifier

    def get_sequence(self):
        """Function that returns the sequence"""
        return self.__sequence

    def get_mw(self):
        """Function that returns the molecular weight of a sequence"""
        mw = sum(self.molecular_weight[aa] for aa in self.__sequence)
        return mw

    def has_subsequence(self, sequence_obj):
        """Function that returns a boolean depending on if a given subsequence is
        inside a sequence"""
        return (sequence_obj.get_sequence() in self.__sequence)

    def __len__(self):
        """Function that returns the length of the sequence"""
        return len(self.__sequence)

    def __eq__(self, other_sequence):
        """Function that compares two sequences and returns a True if equal"""
        if self.__class__.__name__ == other_sequence.__class__.__name__:
            return self.__sequence == other_sequence.__sequence
        else:
            raise TypeError ("Comparison can only be performed in sequences of the same class")

    def __ne__(self, other_sequence):
        """Function that compares two sequences and returns True if different"""
        if self.__class__.__name__ == other_sequence.__class__.__name__:
            return self.__sequence != other_sequence.__sequence
        else:
            raise TypeError ("Comparison can only be performed in sequences of the same class")

    def __add__(self, other_sequence):
        if self.__class__ == other_sequence.__class__:
            new_seq = self.__sequence + other_sequence.__sequence
            new_id = self.__identifier+"+"+other_sequence.__identifier
        else:
            raise TypeError ("Concatenation can only be performed in sequences of the same class")
        return type(self)(new_id, new_seq)

    def __getitem__(self, key):
        return self.__sequence[key]

    def __contains__(self, substring):
        if self.__class__.__name__ == substring.__class__.__name__:
            return substring in self.__sequence
        else:
            raise TypeError ("Function contains can only be performed in sequences of the same class")

    def __lt__ (self, other_sequence):
        if self.__class__.__name__ ==other_sequence.__class__.__name__:
            return self.get_mw() < other_sequence.get_mw()
        else:
            raise TypeError ("Function less than can only be performed in sequences of the same class")

    def __hash__(self):
        key=self.__identifier+self.__sequence
        return key.__hash__()


class ProteinSequence(Sequence):

    molecular_weight = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'

class NucleotideSequence(Sequence):

    table = {}
    stop_codons = ""
    start_codons = ""

    def translate(self):
        """This function returns a protein sequence given a nucleotide sequence"""
        started = False
        translated_sequence = ""
        sequence = self.get_sequence()
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3]
            if started is False:
                if codon in self.start_codons:
                    started = True
                    translated_sequence = self.table[codon]
            elif codon in self.stop_codons:
                break
            else:
                translated_sequence += self.table[codon]
        return ProteinSequence (identifier = self.get_identifier()+"_translated", sequence = translated_sequence)

class RNASequence(NucleotideSequence):

    molecular_weight = {'A': 363.0, 'C': 339.0, 'U': 340.0, 'G': 379.0}
    alphabet = 'GAUC'
    start_codons=['UUG', 'CUG', 'AUG']
    stop_codons=['UAA', 'UAG', 'UGA']
    table={'GUC': 'V', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GUU': 'V', 'AAC': 'N', 'AGG': 'R', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'AGA': 'R', 'AAU': 'N', 'ACU': 'T', 'CAC': 'H', 'GUG': 'V', 'CCG': 'P', 'CCA': 'P', 'AGU': 'S', 'CCC': 'P', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'CGA': 'R', 'CAG': 'Q', 'CGC': 'R', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'CCU': 'P', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GAG': 'E', 'UCC': 'S', 'UAC': 'Y', 'CGU': 'R', 'GAA': 'E', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'UCA': 'S', 'AUG': 'M', 'CUG': 'L', 'AUU': 'I', 'CAU': 'H', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'GAC': 'D', 'GUA': 'V', 'UGC': 'C', 'GCU': 'A', 'UGU': 'C', 'CUC': 'L', 'UUG': 'L', 'UUA': 'L', 'GAU': 'D', 'UUC': 'F'}

    def reverse_transcribe(self):
        """Function that returns the DNA sequence given a RNA sequence"""
        identifier =self.get_identifier()
        sequence = self.get_sequence()
        return DNASequence (identifier, sequence.replace("U", "T"))

class DNASequence(NucleotideSequence):

    molecular_weight = {'A': 347.0, 'C': 323.0, 'T': 322.0, 'G': 363.0}
    alphabet = 'GATC'
    start_codons = ['TTG', 'CTG', 'ATG']
    stop_codons = ['TAA', 'TAG', 'TGA']
    table = {'CTT': 'L', 'ATG': 'M', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'AGC': 'S', 'AGA': 'R', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'ACT': 'T', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'TAC': 'Y', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GAC': 'D', 'GAA': 'E', 'AAG': 'K', 'AAA': 'K', 'AAC': 'N', 'CTC': 'L', 'CAT': 'H', 'AAT': 'N', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'TGT': 'C', 'TCT': 'S', 'GAT': 'D', 'TTT': 'F', 'TGC': 'C', 'TGG': 'W', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TCA': 'S', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A'}

    def transcribe(self):
        """Function that returns the RNA sequence given a DNA sequence"""
        identifier =self.get_identifier()
        sequence = self.get_sequence()
        return RNASequence (identifier, sequence.replace("T", "U"))

class IncorrectSequenceLetter(Exception):
    def __init__(self, letter, class_name):
        self.letter = letter
        self.class_name = class_name

    def __str__(self):
        return "The sequence item %s is not found in the alphabet of %s\n" %(self.letter, self.class_name)

################################################################################

def FASTA_iterator(fasta_filename, name_class):
	fd = open(fasta_filename,"r")
	sequence = ""
	for line in fd:
		if line[0]==">":
			if len(sequence)>0:
				try:
					yield name_class(identifier, sequence)
				except IncorrectSequenceLetter as e:
					print(e, file=sys.stderr)
			identifier = line[1:].strip()
			sequence = ""
		else:
			sequence+=line.strip()
	fd.close()

	if len(sequence)>0:
		try:
			yield name_class(identifier, sequence)
		except IncorrectSequenceLetter as e:
			print(e, file=sys.stderr)

################################################################################
# #the following commands are to create a list of files depending on the terminal input
if __name__=="__main__":
    list_fa = []
    if len(sys.argv) == 1:
        list_all = os.listdir()
        for file in list_all:
            if file.endswith(".fasta") or file.endswith(".fa"):
                list_fa.append(file)
    elif len(sys.argv) > 1:
        if os.path.isdir(sys.argv[1]) == True:
            list_all = os.listdir(sys.argv[1])
            for file in list_all:
                if file.endswith(".fasta") or file.endswith(".fa"):
                    list_fa.append(sys.argv[1]+os.sep+file)
        elif os.path.isfile(sys.argv[1]) == True:
            list_fa.append(sys.argv[1])
    if len(sys.argv) == 3:
        output_file_name = sys.argv[2]
        output_file=open(output_file_name, "w")

list_unsorted = []
list_sorted = []

# now the FASTA_iterator is performed in all the files of the list and the objects
# are added in a list of tupples that is not sorted by molecular weight
sys.stderr.write("%s FASTA files found\n" %str(len(list_fa)))
n=0
for file_fasta in list_fa:
    for object in FASTA_iterator(file_fasta, DNASequence):
        prot_seq=object.translate()
        list_unsorted.append(prot_seq)
        n=n+1
    sys.stderr.write("%s finished\n" %file_fasta)

sys.stderr.write("%s sequences found\n" %str(n))
sys.stderr.write("Sorting the sequences...\n")
# now the list of tupple is sorted by molecular weight
list_sorted=sorted(list_unsorted)
sys.stderr.write("Sort process finished.\n")

#once all the objects are introduced in the list and sorted depending on if there
# is an output file or not, the results are printed on the file or on the terminal
for object in list_sorted:
    if len(sys.argv) == 3:
        output_file.write(object.get_identifier()+"\t"+str(object.__len__())+"\t"+"{:.2f}".format(object.get_mw())+"\n")
    else:
        sys.stdout.write(object.get_identifier()+"\t"+str(object.__len__())+"\t"+"{:.2f}".format(object.get_mw())+"\n")

if len(sys.argv) == 3:
    output_file.close()

sys.stderr.write("Program finished correctly.\n")

# Vicen=DNASequence("Laura", "ATG")
# Vicen2=RNASequence("Clara", "AUG")
#
# ne = Vicen.__contains__(Vicen2)
# print(ne)

# for element in FASTA_iterator("1.fa", DNASequence):
#     print (element)


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
