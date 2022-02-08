################################################################################
####################### EXERCISES Laura Ciaran 192802 ##########################
################################################################################

from math import sqrt
from sys import argv
from sys import stdout
from sys import stdin

def calculate_pdb_chain_mean_minimum_distances(pdb_file_path):
    """Python function that calculates the mean of the minimum distance between
    any two residues pairs found in the same chain of a PDB. It uses a single
    argument corresponding to the PDB file path to use. But is optional. If the
    PDB file path is not defined, read the PDB file from standard input."""
    with open (pdb_file_path, "r") as file1:
        dic = dict()
        output = dict()
        for line in file1:
            if line.startswith("ATOM"):
                chain = line[21]
                residue_id = line[23:27]
                atom_id = line[7:11]
                X_coord = float(line[31:38])
                Y_coord = float(line[39:46])
                Z_coord = float(line[47:54])
                coords = (X_coord, Y_coord, Z_coord)
                if chain not in dic:
                    dic[chain] = {}
                if residue_id not in dic[chain]:
                    dic[chain][residue_id] = {}
                if atom_id not in dic[chain][residue_id]:
                    dic[chain][residue_id][atom_id] = coords
    for chain in dic:
        min_dist_tot = []
        for residue in dic[chain]:
            for residue2 in dic[chain]:
                if residue != residue2:
                    distances = []
                    for atom, coords in dic[chain][residue].items():
                        dist = [sqrt((coords[0]-coords2[0])**2+(coords[1]-coords2[1])**2+(coords[2]-coords2[2])**2) for atom2, coords2 in dic[chain][residue2].items()]
                        distances.extend(dist)
                    min_dist = min(distances)
                    min_dist_tot.append(min_dist)
        mean_min_dist = sum(min_dist_tot) / len(min_dist_tot)
        output[chain] = round(mean_min_dist,4)
    return (output)


if __name__=="__main__":
    if len(argv) == 2:
        pdb_file_path = argv[1]
    else:
        with open ("input.pdb", "w") as file2:
            file_lines = stdin.readlines()
            for line in file_lines:
                file2.write(line)
            pdb_file_path = "input.pdb"

# final_dict = calculate_pdb_chain_mean_minimum_distances(pdb_file_path)
# for chain, mean in final_dict.items():
#     stdout.write(chain+":"+str(mean)+"\n")
