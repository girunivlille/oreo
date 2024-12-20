from Bio import SeqIO
from statistics import mean 
from statistics import median
import sys

def make_dict_from_ids(reads_file):
    dict_ids = {}
    i = 0
    for elem in SeqIO.parse(reads_file, "fasta"):
        dict_ids[elem.id] = i
        i += 1
    return dict_ids

def calculate_distance_by_element(base_file, sorted_file):
    dict_ids = make_dict_from_ids(base_file)
    distances = []
    i = 0
    for read in SeqIO.parse(sorted_file, "fasta"):
        distances.append(abs(i-dict_ids[read.id]))
        i += 1
    
    return distances

def calculate_distance_by_bloc(base_file, sorted_file):
    dict_ids = make_dict_from_ids(base_file)
    distances = []
    i = -1
    begin_sort = 0
    begin_base = 0
    last = -2
    bloc_size = 1
    for read in SeqIO.parse(sorted_file, "fasta"):
        i += 1
        if dict_ids[read.id] == last+1:
            last = dict_ids[read.id]
            bloc_size += 1
        else:
            if last != -2:
                distances.append((bloc_size, abs(begin_base-begin_sort)))
            last = dict_ids[read.id]
            bloc_size = 1
            begin_sort = i
            begin_base = dict_ids[read.id]
    distances.append((bloc_size, abs(begin_base-begin_sort)))
    return distances   

dist_by_elem = calculate_distance_by_element(sys.argv[1], sys.argv[2])
dist_by_bloc = calculate_distance_by_bloc(sys.argv[1], sys.argv[2])

print(dist_by_elem)
print(dist_by_bloc)

print("moyenne:")
print(mean(dist_by_elem))
print("mediane:")
print(median(dist_by_elem))

print("\n\nnombre de blocs:")
print(len(dist_by_bloc))
print("taille moyenne des blocs:")
print(mean([elem[0] for elem in dist_by_bloc]))
print("taille mediane des blocs:")
print(median([elem[0] for elem in dist_by_bloc]))
print("distance moyenne des blocs:")
print(mean([elem[1] for elem in dist_by_bloc]))
print("distance mediane des blocs:")
print(median([elem[1] for elem in dist_by_bloc]))
