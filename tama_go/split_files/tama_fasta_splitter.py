import re
import sys
import time

#
# This script is used to split fasta files for running blastp in parallel
# 
#


fasta_file = sys.argv[1]
fasta_file_contents = open(fasta_file).read().rstrip("\n").split("\n")

outfile_name = sys.argv[2]

split_number = int(sys.argv[3])

seq_list = [] #list of fasta headers
seq_dict = {} #key fasta header value fasta object
seq_count = 0
current_header = ""

for line in fasta_file_contents:
    if line.startswith(">"):
        seq_list.append(line)
        seq_dict[line] = [line]
        current_header = line
        seq_count += 1
    else:
        seq_dict[current_header].append(line)



split_size_float = float(seq_count) / float(split_number)
split_size_int = int(float(seq_count) / float(split_number))

if split_size_float > split_size_int:
    split_size_int = split_size_int + 1

split_size = split_size_int

new_split_index = 0

seq_count = 0
for seq in seq_list:
    if seq_count % split_size == 0:
        new_split_index += 1
        new_outfile_name = outfile_name + "_" + str(new_split_index) + ".fa"
        new_outfile = open(new_outfile_name,"w")
    
    for item in seq_dict[seq]:
        new_outfile.write(item)
        new_outfile.write("\n")
        
    seq_count += 1
    

        