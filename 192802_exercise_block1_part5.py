################################################################################
####################### EXERCISES Laura Ciaran 192802 ##########################
################################################################################

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
                    yield (identifier, sequence)
                sequence = ""
                identifier = line.strip("\n>")
        yield (identifier, sequence)

#for element in FASTA_iterator("1.fa"):
#    print (element)

################################################################################

### Exercise 1 ###
def get_proteins_ratio_by_residue_threshold (filename, residue,
relative_threshold=0.1, absolute_threshold=10):
    """This function returns a float corresponding to the ratio of proteins
    in the fasta file having a relative frequency higher or equal than a given
    threshold and an absolute frequency of the same residue higher or equal
    than a given threshold"""
    prot = {}
    correct_proteins = 0
    total_proteins = 0
    for id, seq in FASTA_iterator(filename):
        total_proteins += 1
        abs_freq = seq.count(residue)
        rel_freq = abs_freq / len (seq)
        if rel_freq > relative_threshold and abs_freq > absolute_threshold:
            correct_proteins += 1
    return (correct_proteins / total_proteins)

#results = get_proteins_ratio_by_residue_threshold("example_fasta_file.fa", "A")
#print (results)

def print_sequence_summary(filename, output_filename, first_n=3, last_m=5):
    """This function saves on a output file the protein identifier, the first
    N-aminoacids, the last M-aminoacids and the absolute frequency in the
    protein of each of the first N-aminoacids and the last M-aminoacids, given
    a protein FASTA file (filename)."""
    with open(output_filename, 'w') as file2:
        for id, seq in FASTA_iterator(filename):
            file2.write(id+"\t")
            start = seq[0:first_n]
            end = seq[-last_m:]
            file2.write(start+"\t")
            file2.write(end+"\t")
            for a in start:
                abs_freq = seq.count(a)
                file2.write(a+":"+str(abs_freq)+",")
            for a in end:
                abs_freq = seq.count(a)
                file2.write(a+":"+str(abs_freq)+",")
            file2.write("\n")

#print_sequence_summary("example_fasta_file.fa", "results.txt")

################################################################################

### Exercise 2 ###
def get_max_sequence_length_from_FASTA_file (fasta_filename):
    """Given a multiline FASTA file, the function returns the length of the
    sequence with the maximum length"""
    lengths = [len(seq) for id, seq in FASTA_iterator(fasta_filename)]
    return(max(lengths))

#results = get_max_sequence_length_from_FASTA_file ("example.fa")
#print (results)

################################################################################

### Exercise 3 ###
def get_min_sequence_length_from_FASTA_file (fasta_filename):
    """Given a multiline FASTA file, the function returns the length of the
    sequence with the minimum length"""
    lengths = [len(seq) for id, seq in FASTA_iterator(fasta_filename)]
    return(min(lengths))

#results = get_min_sequence_length_from_FASTA_file ("example.fa")
#print (results)

################################################################################

### Exercise 4 ###
def get_longest_sequences_from_FASTA_file(fasta_filename):
    """Given a FASTA file, the function returns a list of tuples corresponding
    to the sequence(s) with maximum length."""
    max_length = get_max_sequence_length_from_FASTA_file (fasta_filename)
    longest_seq = [(id, seq) for id, seq in FASTA_iterator(fasta_filename) if len(seq)==max_length]
    return sorted(longest_seq, key=lambda x:x[0].upper())

#results = get_longest_sequences_from_FASTA_file ("example.fa")
#print (results)

################################################################################

### Exercise 5 ###
def get_shortest_sequences_from_FASTA_file(fasta_filename):
    """Given a FASTA file, the function returns a list of tuples corresponding
    to the sequence(s) with minimum length."""
    min_length = get_min_sequence_length_from_FASTA_file (fasta_filename)
    shortest_seq = [(id, seq) for id, seq in FASTA_iterator(fasta_filename) if len(seq)==min_length]
    return sorted(shortest_seq, key=lambda x:x[0].upper())

#results = get_shortest_sequences_from_FASTA_file ("example.fa")
#print (results)

################################################################################

### Exercise 6 ###
def get_molecular_weights(fasta_filename):
    """Given a protein FASTA file, the function returns a dictionary with the
    molecular weights of all the proteins in the file."""
    aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
    dic = {}
    for id, seq in FASTA_iterator(fasta_filename):
        dic[id] = 0
        for aa in seq:
            if aa in aminoacid_mw:
                dic[id] += aminoacid_mw[aa]
    return (dic)

#results = get_molecular_weights("uniprot_sprot_sample.fasta")
#print (results)

################################################################################

### Exercise 7 ###
def get_sequence_with_min_molecular_weight(fasta_filename):
    """Given a protein FASTA file, the function returns a tuple with (identifier,
    sequence) of the protein with the lowest molecular weight."""
    weights = get_molecular_weights(fasta_filename)
    min_weight_seq = [(id, seq) for id,seq in FASTA_iterator(fasta_filename) if weights[id]==min(weights.values())]
    return sorted(min_weight_seq, key=lambda x:x[0].upper())

#results = get_sequence_with_min_molecular_weight("uniprot_sprot_sample.fasta")
#print (results)

################################################################################

### Exercise 8 ###
def get_mean_molecular_weight(fasta_filename):
    weights = get_molecular_weights(fasta_filename)
    return (sum(weights.values()) / len(weights))

#results = get_mean_molecular_weight("uniprot_sprot_sample.fasta")
#print (results)
