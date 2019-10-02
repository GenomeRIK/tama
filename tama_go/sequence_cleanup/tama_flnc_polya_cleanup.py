import re
import sys
import time
from Bio import SeqIO
import argparse


#
# This script removes tailing poly A from FLNC fasta files
# This is necessary for FLNC generated from Teloprime libraries
#


ap = argparse.ArgumentParser(description='This script removes tailing poly A from FLNC fasta files')


ap.add_argument('-f', type=str, nargs=1, help='FLNC Fasta file')
ap.add_argument('-p', type=str, nargs=1, help='Output prefix')


opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0


if not opts.f:
    print("FLNC Fasta file missing")
    missing_arg_flag = 1
if not opts.p:
    print("Output prefix missing")
    missing_arg_flag = 1

if missing_arg_flag == 1:
    print("Please try again with complete arguments")


fasta_file = opts.f[0]

outfile_prefix = opts.p[0]


outfile_name = outfile_prefix + ".fa"
outfile = open(outfile_name,"w")

outfile_report_name = outfile_prefix + "_polya_flnc_report.txt"
outfile_report = open(outfile_report_name,"w")

polya_dict = {} # polya_dict[polya number] = count

print("going through fasta")
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    ##  seq_name = str(seq_record.id)
    seq_name = str(seq_record.description)
    seq_string = str(seq_record.seq)
    seq_string = seq_string.upper()
    
    
    seq_string_list = list(seq_string)
    
    polya_count = 0
    polya_flag = 0
    last_index = -1
    while polya_flag == 0:

        if seq_string_list[last_index] == "A":
            polya_count += 1
        else:
            polya_flag = 1
        
        last_index = last_index - 1
        
    end_index = len(seq_string_list) - polya_count
    
    trim_seq_list = seq_string_list[0:end_index]
    
    trim_seq_string = "".join(trim_seq_list)
    
    if polya_count not in polya_dict:
        polya_dict[polya_count] = 0
    
    polya_dict[polya_count] += 1
    

    outline = ">" + seq_name

    outfile.write(outline)
    outfile.write("\n")
    
    outline = trim_seq_string

    outfile.write(outline)
    outfile.write("\n")

#report summary of poly-A found

polya_num_list = list(polya_dict.keys())
polya_num_list.sort()

outline = "polya_num" + "\t" + "polya_num_count"
outfile_report.write(outline)
outfile_report.write("\n")

for polya_num in polya_num_list:
    
    polya_num_count = polya_dict[polya_num]
    
    # print(str(polya_num) +"\t" + str(polya_num_count))
    
    outline = str(polya_num) +"\t" + str(polya_num_count)
    outfile_report.write(outline)
    outfile_report.write("\n")














