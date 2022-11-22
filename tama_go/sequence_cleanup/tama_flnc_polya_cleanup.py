import re
import sys
import time
from Bio import SeqIO
import argparse


#
# This script removes tailing poly A from FLNC fasta files
# This is necessary for FLNC generated from Teloprime libraries
# Updated 2022/06/03
#


ap = argparse.ArgumentParser(description='This script removes tailing poly A from FLNC fasta files')


ap.add_argument('-f', type=str, nargs=1, help='FLNC Fasta file')
ap.add_argument('-p', type=str, nargs=1, help='Output prefix')
ap.add_argument('-m', type=str, nargs=1, help='Minimum read length to keep (default is 200)')


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


if not opts.m:
    print("Minimum read length to keep (default is 200)")
    min_length = 200
else:
    min_length = int(opts.m[0])

fasta_file = opts.f[0]

outfile_prefix = opts.p[0]


outfile_name = outfile_prefix + ".fa"
outfile = open(outfile_name,"w")

outfile_name = outfile_prefix + "_tails.fa"
outfile_tails = open(outfile_name,"w")

outfile_report_name = outfile_prefix + "_polya_flnc_report.txt"
outfile_report = open(outfile_report_name,"w")

outfile_filtered_name = outfile_prefix + "_discarded_reads.txt"
outfile_filtered = open(outfile_filtered_name,"w")

outfile_summary_name = outfile_prefix + "_summary.txt"
outfile_summary = open(outfile_summary_name,"w")

polya_dict = {} # polya_dict[polya number] = count

non_a_buffer_threshold = 2
a_percent_threshold = 0.7



class Block:
    def __init__(self, block_index,block_type,block_count,block_start,block_end):
        self.block_index = block_index
        self.block_type = block_type
        self.block_count = block_count
        self.block_start = block_start
        self.block_end = block_end


total_read_count = 0
pass_read_count = 0
discarded_read_count = 0
polya_zero_read_count = 0
polya_poly_five_read_count = 0
polya_poly_sixmore_read_count = 0
longest_trimmed_read_length = 0
longest_trimmed_read_id = ""
peak_read_length = 0
longest_polya_length = 0

read_length_dict = {} # read_length_dict[length bin] = counts

length_bin_size = 200

print("going through fasta")
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    ##  seq_name = str(seq_record.id)
    seq_name = str(seq_record.description)
    seq_string = str(seq_record.seq)
    seq_string = seq_string.upper()
    
    
    seq_string_list = list(seq_string)

    total_read_count += 1
    
    #########################################
    a_block_dict = {} # a_block_dict[block index] = block_obj
    block_index = 0
    this_block_flag = "NA"
    block_start = -1
    block_end = -1
    this_block_count = 0
    
    #print(seq_name)
    
    # use A and non A block method
    for i in range(len(seq_string_list)):
        
        block_end += 1
        this_block_count += 1
        
        last_index = len(seq_string_list) - 1 - i
        seq_nuc = seq_string_list[last_index]
        
        if this_block_flag == "NA" :
            if seq_nuc == "A":
                this_block_flag = "a_block"
            else:
                this_block_flag = "not_a_block"
            
            block_start = 0
            this_block_count = 0
            continue
    
        if seq_nuc == "A" and this_block_flag == "a_block": # same a block
            continue
        
        if seq_nuc != "A" and this_block_flag == "not_a_block": # same not a block
            continue
        
        if seq_nuc == "A" and this_block_flag == "not_a_block": # block change
            
            block_type = this_block_flag
            
            block_obj = Block(block_index,block_type,this_block_count,block_start,block_end)
            
            a_block_dict[block_index] = block_obj
            
            block_index += 1
            this_block_count = 0
            block_start = block_end
            this_block_flag = "a_block"
            continue
        
        if seq_nuc != "A" and this_block_flag == "a_block": # block change
            block_type = this_block_flag
            
            block_obj = Block(block_index,block_type,this_block_count,block_start,block_end)
            a_block_dict[block_index] = block_obj
            
            block_index += 1
            this_block_count = 0
            block_start = block_end
            this_block_flag = "not_a_block"
            continue
    
    tail_index = 0
    this_block_index = -1
    stop_checking_flag = 0
    polya_count = 0
    tail_count = 0
    this_a_percent = 1.0
    
    simple_a_tail_index = 0

    block_overrun_flag = 0
    
    while stop_checking_flag == 0:
        
        this_block_index += 1

        # check if we have run out blocks
        # this happens when the reads are too short
        if this_block_index not in a_block_dict:
            stop_checking_flag = 1
            block_overrun_flag = 1
            continue
        
        this_block_obj = a_block_dict[this_block_index]
        
        this_block_type = this_block_obj.block_type
        this_block_count = this_block_obj.block_count
        this_block_start = this_block_obj.block_start
        this_block_end = this_block_obj.block_end
        
        if this_block_index > 0:
            this_a_percent = float(polya_count) / float(tail_count)
        tail_count = tail_count + this_block_count
        
        #print(this_block_type + ": " + str(this_block_count))
        
        if this_block_type == "a_block":
            
            if this_block_index == 0:
                simple_a_tail_index = -1 * this_block_end
            
            if this_block_index+1 not in a_block_dict: # check for running out of index
                tail_index = -1 * this_block_end 
                stop_checking_flag = 1
                continue
            
            next_block_obj = a_block_dict[this_block_index+1]
            next_block_type = next_block_obj.block_type
            next_block_count = next_block_obj.block_count
            
            polya_count = polya_count + this_block_count
            
            if next_block_count > 2: # Stop tail region 
                tail_index = -1 * this_block_end 
                stop_checking_flag = 1
        
        if this_block_type == "not_a_block":
            
            if this_block_index+1 not in a_block_dict: # check for running out of index
                tail_index = -1 * this_block_start 
                stop_checking_flag = 1
                continue
                
            next_block_obj = a_block_dict[this_block_index+1]
            next_block_type = next_block_obj.block_type
            next_block_count = next_block_obj.block_count
            
            if this_block_count > 2: # Stop tail region 
                tail_index = -1 * this_block_start 
                stop_checking_flag = 1
                continue
            elif next_block_count < 2:
                
                if this_block_index+3 not in a_block_dict: # check for running out of index
                    tail_index = -1 * this_block_start 
                    stop_checking_flag = 1
                    continue
            
                b_next_block_obj = a_block_dict[this_block_index+2]
                b_next_block_type = b_next_block_obj.block_type
                b_next_block_count = b_next_block_obj.block_count
                
                c_next_block_obj = a_block_dict[this_block_index+3]
                c_next_block_type = c_next_block_obj.block_type
                c_next_block_count = c_next_block_obj.block_count
                
                if b_next_block_count > 1 or c_next_block_count < 3:    
                    tail_index =  -1 * this_block_start 
                    stop_checking_flag = 1
        
        if this_block_index > 2:
            if this_a_percent < a_percent_threshold:
                if this_block_type == "a_block":
                    tail_index =  -1 * this_block_end 
                if this_block_type == "not_a_block":
                    tail_index =  -1 * this_block_start 
                stop_checking_flag = 1
        
        if this_block_index >= block_index:
            stop_checking_flag = 1
                
  
                
    #########################################
        


    # check that there is enough A in the tail string to remove
    tail_string = seq_string_list[tail_index:]
    a_count = tail_string.count("A")
    a_percent = float(a_count) / float(len(tail_string))
    #print(seq_string_list)
    #print(tail_string)


    if a_percent < a_percent_threshold:
        if simple_a_tail_index != 0:
            tail_index = simple_a_tail_index
        else:
            tail_index = 0


    if tail_index > -1:
        trim_seq_list = seq_string_list[0:]
        tail_seq_list = []
    else:
        trim_seq_list = seq_string_list[0:tail_index]
        tail_seq_list = seq_string_list[tail_index:]
    
    trim_seq_string = "".join(trim_seq_list)
    tail_seq_string = "".join(tail_seq_list)
    
    if polya_count not in polya_dict:
        polya_dict[polya_count] = 0
    
    polya_dict[polya_count] += 1



    if len(trim_seq_list) <  min_length:

        outline = ">" + seq_name + "\t" + "short"  + "\t" + "trimlen:"+ str(len(trim_seq_list)) + "\t" + "prelen:"+ str(len(seq_string_list))
        outfile_filtered.write(outline)
        outfile_filtered.write("\n")

        outline = seq_string
        outfile_filtered.write(outline)
        outfile_filtered.write("\n")

        discarded_read_count += 1

    elif block_overrun_flag == 1:

        outline = ">" + seq_name + "\t" + "overrun"  + "\t" + "trimlen:"+ str(len(trim_seq_list)) + "\t" + "prelen:"+ str(len(seq_string_list))
        outfile_filtered.write(outline)
        outfile_filtered.write("\n")

        outline = seq_string
        outfile_filtered.write(outline)
        outfile_filtered.write("\n")

        discarded_read_count += 1


    else:

        outline = ">" + seq_name


        outfile.write(outline)
        outfile.write("\n")

        outline = trim_seq_string

        outfile.write(outline)
        outfile.write("\n")

        pass_read_count += 1

        this_trimmed_read_length = len(trim_seq_string)

        if this_trimmed_read_length > longest_trimmed_read_length:
            longest_trimmed_read_length = this_trimmed_read_length
            longest_trimmed_read_id = seq_name

        bin_num = this_trimmed_read_length // length_bin_size

        bin_length = bin_num * length_bin_size

        if bin_length not in  read_length_dict:
            read_length_dict[bin_length] = 0

        read_length_dict[bin_length] += 1





    outline = ">tail_" + seq_name
    
    outfile_tails.write(outline)
    outfile_tails.write("\n")
    
    outline = tail_seq_string
    
    outfile_tails.write(outline)
    outfile_tails.write("\n")

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

    if polya_num == 0 :

        polya_zero_read_count = polya_num_count

    elif polya_num > 0 and polya_num < 6 :
        polya_poly_five_read_count = polya_poly_five_read_count + polya_num_count

    elif polya_num > 6 :
        polya_poly_sixmore_read_count = polya_poly_sixmore_read_count + polya_num_count


    if polya_num > longest_polya_length:
        longest_polya_length = polya_num



read_length_bin_list = list(read_length_dict.keys())


max_bin_count = 0
max_bin_length = 0

for bin_length in read_length_bin_list:

    bin_count = read_length_dict[bin_length]

    if bin_count > max_bin_count:
        max_bin_count = bin_count
        max_bin_length = bin_length



peak_read_length = max_bin_length
peak_read_count = max_bin_count


# write out summary report


outline = "total_read_count:" + "\t" + str(total_read_count)

outfile_summary.write(outline)
outfile_summary.write("\n")

outline = "pass_read_count:" + "\t" + str(pass_read_count)

outfile_summary.write(outline)
outfile_summary.write("\n")

outline = "discarded_read_count:" + "\t" + str(discarded_read_count)

outfile_summary.write(outline)
outfile_summary.write("\n")

outline = "polya_zero_read_count:" + "\t" + str(polya_zero_read_count)

outfile_summary.write(outline)
outfile_summary.write("\n")

outline = "polya_upto_five_read_count:" + "\t" + str(polya_poly_five_read_count)

outfile_summary.write(outline)
outfile_summary.write("\n")

outline = "polya_sixormore_read_count:" + "\t" + str(polya_poly_sixmore_read_count)

outfile_summary.write(outline)
outfile_summary.write("\n")


outline = "longest_trimmed_read_length:" + "\t" + str(longest_trimmed_read_length) + "\t" + "longest_trimmed_read_id:" + "\t" + str(longest_trimmed_read_id)

outfile_summary.write(outline)
outfile_summary.write("\n")

outline = "peak_read_length:" + "\t" + str(peak_read_length) + "\t" + "peak_read_count:" + "\t" + str(peak_read_count)

outfile_summary.write(outline)
outfile_summary.write("\n")

outline = "longest_polya_length:" + "\t" + str(longest_polya_length)

outfile_summary.write(outline)
outfile_summary.write("\n")

