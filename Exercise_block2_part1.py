################################################################################
####################### EXERCISES Laura Ciaran 192802 ##########################
################################################################################

### Exercise 1 ###
class Protein (object):

    def __init__(self, identifier, sequence):
        self.identifier = identifier
        self.sequence = sequence

    def get_identifier(self):
        """Function that returns the id"""
        return self.identifier

    def get_sequence(self):
        """Function that returns the sequence"""
        return self.sequence

    def get_mw(self):
        """Function that returns the molecular weight of a sequence"""
        aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
        mw = sum(aminoacid_mw[aa] for aa in self.sequence)
        return (mw)

    def has_subsequence(self,subsequence):
        """Function that returns a boolean depending on if a given subsequence is
        inside a sequence"""
        return (subsequence.upper() in self.sequence)

    def get_length(self):
        """Function that returns the length of a sequence"""
        return (len(self.sequence))

# Vicen=Protein("B6I411", "MSVPLILTILAGAATFIGAFLGVLGQKPSNRLLAFSLGFAAGIMLLISLMEMLPAALAAE")
# id=Vicen.get_identifier()
# seq=Vicen.get_sequence()
# mw=Vicen.get_mw()
# subs=Vicen.has_subsequence("LAGaa")
# length=Vicen.get_length()
#
# print(id)
# print(seq)
# print(mw)
# print(subs)
# print(length)


### Exercise 2###
def FASTA_iterator(fasta_filename):
    """Generator Function that reads a Fasta file. In each iteration, returns a
    tuple with the following format: (identifier, sequence)."""
    with open (fasta_filename, "r") as file1:
        sequence = ""
        identifier = ""
        for line in file1:
            if line[0] != ">":
                sequence += line.strip("\n")
            else:
                if len(sequence) > 0:
                    yield (Protein(identifier, sequence))
                sequence = ""
                identifier = line.strip("\n>")
        yield (Protein(identifier, sequence))


# for element in FASTA_iterator("1.fa"):
#     print (element)

# if __name__ == "__main__":
# 	for element in FASTA_iterator( "uniprot_sprot_sample.fasta" ) :
# 		print(element)
# 		print(type(element))
# 		print (element.get_identifier())
# 		print (element.get_sequence())
# 		print (element.get_length())
# 		print (element.get_mw())
# 		print (element.has_subsequence("LVE"))
# 		input("")
