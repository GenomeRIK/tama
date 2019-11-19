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

non_a_buffer_threshold = 2
a_percent_threshold = 0.7


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

    tail_index = -1

    non_a_buffer_count = 0

    while polya_flag == 0:


        if seq_string_list[last_index] == "A":
            polya_count += 1
            tail_index = last_index
            # non_a_buffer_count = 0 # reset non A buffer count

        elif non_a_buffer_count < non_a_buffer_threshold: # allow for non A buffer of up to 2 nucleotides of non A
            non_a_buffer_count += 1

        else:
            polya_flag = 1
        
        last_index = last_index - 1
        
    # end_index = len(seq_string_list) - polya_count

    # check that we did not go into the real sequence
    pre_tail_index = tail_index - 1
    # for i in xrange(non_a_buffer_threshold):

    check_boundary_count = 0
    end_poly_a_flag = 0
    stop_checking_flag = 0
    while stop_checking_flag == 0:
    # while end_poly_a_flag == 0 and check_boundary_count < non_a_buffer_threshold:

        pre_tail_index = pre_tail_index + 1
        check_boundary_count += 1
        print(str(pre_tail_index) + ": " + str(seq_string_list[pre_tail_index]))


        if pre_tail_index > -1: # dont let pre tail index wrap around
            pre_tail_index = -1
            end_poly_a_flag = 1
            check_boundary_count = non_a_buffer_threshold
            stop_checking_flag = 1

            if seq_string_list[-1] == "A":
                tail_index = -1
            else:
                tail_index = 0

        elif seq_string_list[pre_tail_index] != "A":
            tail_index = pre_tail_index + 1
            end_poly_a_flag = 0
        else:
            end_poly_a_flag = 1

        if end_poly_a_flag == 1 and check_boundary_count >= non_a_buffer_threshold:
            stop_checking_flag = 1



    # check that there is enough A in the tail string to remove
    tail_string = seq_string_list[tail_index:]
    a_count = tail_string.count("A")
    a_percent = float(a_count) / float(len(tail_string))

    if a_percent < a_percent_threshold:
        tail_index = 0


    if tail_index > -1:
        trim_seq_list = seq_string_list[0:]
    else:
        trim_seq_list = seq_string_list[0:tail_index]
    
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














