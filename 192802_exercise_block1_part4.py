################################################################################
####################### EXERCISES Laura Ciaran 192802 ##########################
################################################################################

### Exercise 1 ###
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

### Exercise 2 ###
def compare_fasta_file_identifiers(fasta_filenames_list):
    """With a list of FASTA files, the function returns a dictionary that
    contains the 4 following keys with the associated values:
    “intersection”: a set with the common identifiers found in all the files
    “union”: a set with all the identifiers (unique) found in all the files
    “frequency”: a dictionary with all the identifiers as keys and the number
    of files in which it appears as values (int)
    “specific”: a dictionary with the name of the input files as keys and a set
    with the specific identifiers as values"""

    final_dic = {"Intersection": set(),
                "Union": set(),
                "Frequency": dict(),
                "Specific": dict()}

    ids = {}

    for file in fasta_filenames_list:
        ids[file] = set()
        for id, seq in FASTA_iterator(file):
            ids[file].add(id.upper())

    n = 0

    for file, prot in ids.items():
        n += 1
        if n == 1:
            intersection = set(prot)
            union = set(prot)
        else:
            intersection.intersection_update(prot)
            union.update(prot)
        for element in prot:
            if element not in final_dic["Frequency"]:
                final_dic["Frequency"][element] = 1
            else:
                final_dic["Frequency"][element] += 1

        specific = set(prot)
        for file2, prot2 in ids.items():
            if file2 == file:
                continue
            else:
                 specific.difference_update(prot2)

        final_dic["Specific"][file] = specific

    final_dic["Intersection"] = intersection

    final_dic["Union"] = union

    return (final_dic)



#results = compare_fasta_file_identifiers(["1.fa", "2.fa"])
#print (results)
