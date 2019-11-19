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

outfile_name = outfile_prefix + "_tails.fa"
outfile_tails = open(outfile_name,"w")

outfile_report_name = outfile_prefix + "_polya_flnc_report.txt"
outfile_report = open(outfile_report_name,"w")

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


print("going through fasta")
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    ##  seq_name = str(seq_record.id)
    seq_name = str(seq_record.description)
    seq_string = str(seq_record.seq)
    seq_string = seq_string.upper()
    
    
    seq_string_list = list(seq_string)
    
    
    #########################################
    a_block_dict = {} # a_block_dict[block index] = block_obj
    block_index = 0
    this_block_flag = "NA"
    block_start = -1
    block_end = -1
    this_block_count = 0
    
    #print(seq_name)
    
    # use A and non A block method
    for i in xrange(len(seq_string_list)):
        
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
    
    while stop_checking_flag == 0:
        
        this_block_index += 1
        
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
    

    outline = ">" + seq_name


    outfile.write(outline)
    outfile.write("\n")
    
    outline = trim_seq_string

    outfile.write(outline)
    outfile.write("\n")
    
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














