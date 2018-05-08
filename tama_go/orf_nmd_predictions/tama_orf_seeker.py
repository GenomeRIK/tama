import re
import sys
import time
from Bio import SeqIO
import argparse


#
# This script finds open reading frames from transcript sequences 
# It outputs the location and protein sequence for each open reading frame
#

#
# First counts from beginnning in all frames
# The counts starting from M and ending in stop codon
# If no stop codon is found the reading frame is not reported (due to poly A tail selection)
#

#####################################################################
#####################################################################

ap = argparse.ArgumentParser(description='This script finds open reading frames from transcript sequences')

ap.add_argument('-f', type=str, nargs=1, help='Fasta file (required)')
ap.add_argument('-o', type=str, nargs=1, help='Output file name (required)')


opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0

if not opts.f:
    print("Fasta file missing")
    missing_arg_flag = 1
if not opts.o:
    print("Output name missing")
    missing_arg_flag = 1

if missing_arg_flag == 1:
    print("Please try again with complete arguments")

fasta_file = opts.f[0]
outfile_name = opts.o[0]

outfile = open(outfile_name,"w")


amino_dict = {} #amino_dict[GCA] = A letter amino acid code

amino_dict["GTT"] = "V"
amino_dict["GTC"] = "V"
amino_dict["GTA"] = "V"
amino_dict["GTG"] = "V"

amino_dict["GCT"] = "A"
amino_dict["GCC"] = "A"
amino_dict["GCA"] = "A"
amino_dict["GCG"] = "A"

amino_dict["GAT"] = "D"
amino_dict["GAC"] = "D"

amino_dict["GAA"] = "E"
amino_dict["GAG"] = "E"

amino_dict["GGT"] = "G"
amino_dict["GGC"] = "G"
amino_dict["GGG"] = "G"
amino_dict["GGA"] = "G"

amino_dict["TTT"] = "F"
amino_dict["TTC"] = "F"

amino_dict["TTA"] = "L"
amino_dict["TTG"] = "L"

amino_dict["TCT"] = "S"
amino_dict["TCC"] = "S"
amino_dict["TCA"] = "S"
amino_dict["TCG"] = "S"

amino_dict["TAT"] = "Y"
amino_dict["TAC"] = "Y"
amino_dict["TAA"] = "$"  #stop
amino_dict["TAG"] = "$"  #stop

amino_dict["TGT"] = "C"
amino_dict["TGC"] = "C"
amino_dict["TGA"] = "$"  #stop
amino_dict["TGG"] = "W"

amino_dict["CTT"] = "L"
amino_dict["CTC"] = "L"
amino_dict["CTA"] = "L"
amino_dict["CTG"] = "L"

amino_dict["CCT"] = "P"
amino_dict["CCC"] = "P"
amino_dict["CCA"] = "P"
amino_dict["CCG"] = "P"

amino_dict["CAT"] = "H"
amino_dict["CAC"] = "H"

amino_dict["CAA"] = "Q"
amino_dict["CAG"] = "Q"

amino_dict["CGT"] = "R"
amino_dict["CGC"] = "R"
amino_dict["CGA"] = "R"
amino_dict["CGG"] = "R"

amino_dict["ATT"] = "I"
amino_dict["ATC"] = "I"
amino_dict["ATA"] = "I"

amino_dict["ATG"] = "M"  #start

amino_dict["ACT"] = "T"
amino_dict["ACC"] = "T"
amino_dict["ACA"] = "T"
amino_dict["ACG"] = "T"

amino_dict["AAT"] = "N"
amino_dict["AAC"] = "N"

amino_dict["AAA"] = "K"
amino_dict["AAG"] = "K"

amino_dict["AGT"] = "S"
amino_dict["AGC"] = "S"
amino_dict["AGA"] = "R"
amino_dict["AGG"] = "R"



class Orf:
    def __init__(self,orf_seq,a_start,a_end,n_start,n_end,frame):
        self.orf_seq = orf_seq
        self.a_start = a_start
        self.a_end = a_end
        self.frame = frame
        self.start_codon = orf_seq[0]
        self.orf_length = len(orf_seq)

        self.n_start = n_start
        self.n_end = n_end
        
        #self.n_start = (a_start - 1) * 3 + frame
        #self.n_end = (a_end - 1) * 3 + frame  + 2

        self.orf_id = "_".join([str(frame), str(a_start), str(a_end), str(self.n_start), str(self.n_end), str(self.orf_length),self.start_codon])
        



#make amino acid conversion for one frame
def frame_iterate(seq,frame):
    seq_length = len(seq)
    pos = frame - 1
    prot_seq = []
    n_flag = 0
    n_pos_list = [] #keep track of nuc positions

    while pos < seq_length - 2:
        n_pos_list.append(pos)
        codon = seq[pos:pos+3]
        if codon in amino_dict:
            prot_seq.append(amino_dict[codon])
        else:
            print("sequence has N's")
            n_flag = 1
            break
            
        pos = pos + 3
    return [prot_seq,n_flag,n_pos_list]

#use prot seq from frame_iterate to find orf's
def orf_seeker(prot_seq,frame,n_pos_list):
    
    #orf_list = [] # list of orf objects

    orf_dict = {} # orf_dict[orf id] = orf obj
    orf_seq = []
    orf_start = 0
    
    #find first orf, as in if there is 5' degradation
    for i,amino_acid in enumerate(prot_seq):
        if amino_acid == "$":
            if i > 0:
                n_start = n_pos_list[0]
                n_end = n_pos_list[i] + 2

                orf_start = orf_start + 1

                orf_end = i
                orf_obj = Orf(orf_seq,orf_start,orf_end,n_start,n_end,frame)
                orf_id = orf_obj.orf_id
                orf_dict[orf_id] = orf_obj
                #orf_list.append(Orf(orf_seq,orf_start,orf_end,frame))
                break
            else:
                break
        
        orf_seq.append(amino_acid)
    
    
    orf_seq = []
    orf_start = 0
    #find orf's starting with Met
    for i,amino_acid in enumerate(prot_seq):
        #start orf seq
        if len(orf_seq) ==  0:
            if amino_acid == "M":
                orf_start = i
                orf_seq.append(amino_acid)
            continue

        #end orf seq
        if amino_acid == "$":
            if len(orf_seq) >  0:

                n_start = n_pos_list[orf_start]
                n_end = n_pos_list[i] + 2

                orf_start = orf_start + 1

                orf_end = i
                orf_obj = Orf(orf_seq, orf_start, orf_end,n_start,n_end, frame)
                orf_id = orf_obj.orf_id
                orf_dict[orf_id] = orf_obj
                #orf_list.append(Orf(orf_seq,orf_start,orf_end,frame))
                #refresh orf variables
                orf_seq = []
                orf_start = i+1
            continue    
        #continue adding rof seq
        orf_seq.append(amino_acid)

    orf_id_list = list(orf_dict.keys())
    orf_id_list.sort()

    orf_list = []  # list of orf objects
    for orf_id in orf_id_list:

        this_orf_obj = orf_dict[orf_id]
        orf_list.append(this_orf_obj)

    return orf_list


def sort_orf_list(orf_list1,orf_list2,orf_list3):
    orf_size_dict = {} # orf_size_dict[size] = list of orf objects
    
    all_orf_list = orf_list1 + orf_list2 + orf_list3
    
    
    for orf_obj in all_orf_list:
        orf_length = orf_obj.orf_length
        
        if orf_length not in orf_size_dict:
            orf_size_dict[orf_length] = []
            
        orf_size_dict[orf_length].append(orf_obj)
    
    orf_size_list = orf_size_dict.keys()
    orf_size_list.sort(reverse=True)
    
    top_orf_list = []
    
    #get longest 4 orf's, may be less if there aren't that many due to the stop codon assumption
    for i in xrange(4):
        if i < len(orf_size_list):
            for orf_obj in orf_size_dict[orf_size_list[i]]:
                
                top_orf_list.append(orf_obj)
    
    return top_orf_list


count = 0
print("Going through sequences")
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    if count % 1000 == 0:
        print(count)
    count += 1
    
    trans_id = str(seq_record.id)
    trans_seq = str(seq_record.seq)
    trans_seq = trans_seq.upper()
    
    [f1_prot_seq,f1_n_flag,f1_n_pos_list] = frame_iterate(trans_seq,1)
    [f2_prot_seq,f2_n_flag,f2_n_pos_list] = frame_iterate(trans_seq,2)
    [f3_prot_seq,f3_n_flag,f3_n_pos_list] = frame_iterate(trans_seq,3)
    
    all_n_flag = f1_n_flag + f2_n_flag + f3_n_flag
    
    if all_n_flag > 0:
        outline = ">" + trans_id + ":" + "missing_nucleotides"
        outfile.write(outline)
        outfile.write("\n")
        continue
    
    f1_orf_list = orf_seeker(f1_prot_seq,1,f1_n_pos_list)
    f2_orf_list = orf_seeker(f2_prot_seq,2,f2_n_pos_list)
    f3_orf_list = orf_seeker(f3_prot_seq,3,f3_n_pos_list)
    
    top_orf_list = sort_orf_list(f1_orf_list,f2_orf_list,f3_orf_list)
    
    for orf_obj in top_orf_list:
        frame_flag = "F" + str(orf_obj.frame)
        outline = ":".join([trans_id,frame_flag,str(orf_obj.n_start),str(orf_obj.n_end),str(orf_obj.a_start),str(orf_obj.a_end),str(orf_obj.orf_length),orf_obj.start_codon])
        #outline = ":".join([trans_id,frame_flag,str(orf_obj.a_start),str(orf_obj.a_end),str(orf_obj.orf_length),orf_obj.start_codon])
        outline = ">" + outline
        outfile.write(outline)
        outfile.write("\n")
        
        outline = "".join(orf_obj.orf_seq)
        outfile.write(outline)
        outfile.write("\n")

    
    
    





