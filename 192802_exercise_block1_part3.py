################################################################################
####################### EXERCISES Laura Ciaran 192802 ##########################
################################################################################

### Exercise 1 ###
def calculate_aminoacid_frequencies (fasta_filename, subsequences_filename, number_of_repetitions, output_filename):
    """This function given a multi-line protein FASTA file and a “sub-sequences”
    file (one sequence in each line), calculates the proportion of proteins in
    the FASTA file containing at least N-times each of the sub-sequences."""
    subseqs = []
    repetitions = []
    with open (subsequences_filename, "r") as fd:
        for line in fd:
            subseqs.append(line.strip())
    with open (fasta_filename, "r") as fs:
        name = ""
        dic = {}
        num_prot = 0
        for line in fs:
            if ">" in line:
                num_prot += 1
                name = line.strip("\n>")
                seq = ""
            else:
                seq = seq + line.strip("\n")
                dic[name] = seq
    for seq in subseqs:
        counter = 0
        for prot, aa in dic.items():
            if seq in aa:
                if aa.count(seq) >= number_of_repetitions:
                    counter += 1
        repetitions.append(counter)
    num_subseq = len(repetitions)
    with open(output_filename, "w") as results_file:
        right_aligned_num_prot = f"{num_prot:>20}"
        right_aligned_num_subseq = f"{num_subseq:>16}"
        results_file.write("#Number of proteins:"+right_aligned_num_prot+"\n")
        results_file.write("#Number of subsequences:"+right_aligned_num_subseq+"\n")
        results_file.write("#Subsequence proportions:\n")
        merged = {}
        for key in subseqs:
            for value in repetitions:
                merged[key] = value
                repetitions.remove(value)
                break
        sort_merged = sorted(merged.items(), key=lambda x: x[1], reverse=True)
        for prot in sort_merged:
            proportions = "{:.4f}".format(prot[1]/num_prot)
            results_file.write(str(prot[0])+"\t"+f"{prot[1]:>10}"+"\t"+f"{proportions:>17}"+"\n")

calculate_aminoacid_frequencies ("example_fasta_file.fa", "sequence_fragments.txt", 5, "results.txt")
