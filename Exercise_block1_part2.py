### Exercise 1 ###
def get_proteins_ratio_by_residue_threshold (filename, residue,
relative_threshold=0.1, absolute_threshold=10):
    """This function returns a float corresponding to the ratio of proteins
    in the fasta file having a relative frequency higher or equal than a given
    threshold and an absolute frequency of the same residue higher or equal
    than a given threshold"""
    with open (filename, "r") as file1:
        name = ""
        dic = {}
        correct_proteins = 0
        total_proteins = 0
        for line in file1:
            if ">" in line:
                name = line.strip("\n>")
                seq = ""
                total_proteins += 1
            else:
                seq = seq + line.strip("\n")
                dic[name] = seq
        for prot, aa in dic.items():
            abs_freq = aa.count(residue)
            rel_freq = abs_freq / len (aa)
            if rel_freq > relative_threshold and abs_freq > absolute_threshold:
                correct_proteins += 1
        return (correct_proteins / total_proteins)

results = get_proteins_ratio_by_residue_threshold ("example_fasta_file.fa", "A")
print (results)

### Exercise 2 ###
def print_sequence_summary(filename, output_filename, first_n=3, last_m=5):
    """This function saves on a output file the protein identifier, the first
    N-aminoacids, the last M-aminoacids and the absolute frequency in the
    protein of each of the first N-aminoacids and the last M-aminoacids, given
    a protein FASTA file (filename)."""
    with open (filename, "r") as file1:
        name = ""
        dic = {}
        for line in file1:
            if ">" in line:
                name = line.strip("\n>")
                seq = ""
            else:
                seq = seq + line.strip("\n")
                dic[name] = seq
    with open(output_filename, 'w') as file2:
        for prot, aa in dic.items():
            file2.write(prot+"\t")
            start = aa[0:first_n]
            end = aa[-last_m:]
            file2.write(start+"\t")
            file2.write(end+"\t")
            for a in start:
                abs_freq = aa.count(a)
                file2.write(a+":"+str(abs_freq)+",")
            for a in end:
                abs_freq = aa.count(a)
                file2.write(a+":"+str(abs_freq)+",")
            file2.write("\n")

print_sequence_summary("test.txt", "results.txt")
