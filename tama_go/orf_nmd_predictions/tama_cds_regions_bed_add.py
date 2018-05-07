import re
import sys
import time
from Bio import SeqIO



#
# This script uses data from the blastp parse file and the pacbio annotation to assign the locations of the UTR/CDS regions to the bed file
# Use sequence file to double check the sequence lengths
#


print("opening blastp parse file")
parse_file = sys.argv[1]
parse_file_contents = open(parse_file).read().rstrip("\n").split("\n")

print("opening pbri anno file")
pbri_file = sys.argv[2]
pbri_file_contents = open(pbri_file).read().rstrip("\n").split("\n")


print("opening pbri fasta file")
fasta_file = sys.argv[3]
#fasta_file_contents = open(fasta_file).read().rstrip("\n").split("\n")

outfile_name = sys.argv[4]
outfile = open(outfile_name,"w")

trans_dict = {} # trans_dict[trans_id] = trans line list

for line in parse_file_contents:
    line_split = line.split("\t")
    trans_id = line_split[0]
    trans_dict[trans_id] = line_split
    
trans_len_dict = {} # trans_len_dict[trans id] = length
print("Going through sequences")
for seq_record in SeqIO.parse(fasta_file, "fasta"):

    trans_id = str(seq_record.id).split(":")[0]
    trans_seq = str(seq_record.seq)
    trans_seq = trans_seq.upper()
    trans_length = len(trans_seq)
    trans_len_dict[trans_id] = trans_length


    
    
for line in pbri_file_contents:
    line_split = line.split("\t")
    trans_id = line_split[3]
    block_list = line_split[10].split(",")


    if trans_id not in trans_dict:
        line_split[6] = "0"
        line_split[7] = "0"

        degrade_flag = "missing"
        nmd_flag = "missing"

        prot_id = "none"
        match_flag = "no_orf"
        
        frame = "na"

        line_split[3] = line_split[3] + ";" + prot_id + ";" + degrade_flag + ";" + match_flag + ";" + nmd_flag + ";" + frame

        outline = "\t".join(line_split)
        outfile.write(outline)
        outfile.write("\n")


        continue


    frame = trans_dict[trans_id][1]
    prot_id = trans_dict[trans_id][6]
    match_flag = trans_dict[trans_id][7]
    trans_start = int(line_split[1])
    trans_end = int(line_split[2])
    block_start_list = line_split[11].split(",")
    strand = line_split[5]
    
    prot_start = int(trans_dict[trans_id][4]) - 1
    prot_end = int(trans_dict[trans_id][5]) - 1

    nuc_start = int(trans_dict[trans_id][2])
    nuc_end = int(trans_dict[trans_id][3])

    #if this transcript has missing nucleotides frame will be "no_frame"
    if match_flag == "missing_nucleotides":
        
        line_split[6] = "0"
        line_split[7] = "0"
        
        degrade_flag = "missing"
        nmd_flag = "missing"
        
        frame = "missing"
        
        
        line_split[3] = line_split[3] + ";" + prot_id + ";" + degrade_flag + ";" + match_flag + ";" + nmd_flag + ";" + frame

        outline = "\t".join(line_split)
        outfile.write(outline)
        outfile.write("\n")
        
        continue
    
    pos_adjust = 0
    #if frame == "F2":
    #    pos_adjust = 1
    #elif frame == "F3":
    #    pos_adjust = 2
    
    if strand == "+":
        #cds_rel_start =  prot_start * 3 + pos_adjust # minus 1 to adjust to bed coords
        #cds_rel_end =  prot_end * 3 + pos_adjust + 3  # minus 1 to adjust to bed coords

        cds_rel_start = nuc_start + pos_adjust
        cds_rel_end = nuc_end + pos_adjust - 2

    elif strand == "-":
        #cds_rel_start = trans_len_dict[trans_id] - (prot_end * 3 + pos_adjust + 3)
        #cds_rel_end = trans_len_dict[trans_id] - (prot_start * 3 + pos_adjust + 1)

        cds_rel_start = trans_len_dict[trans_id] - (nuc_end + pos_adjust - 2)
        cds_rel_end = trans_len_dict[trans_id] - (nuc_start + pos_adjust )
        
        
    cds_coord_start = 0
    cds_coord_end = 0
        
    #get sequence length by summing up blocks
    block_sum = 0
    exon_cds_start = 0
    exon_cds_end = 0
    
    for i, block_size in enumerate(block_list):
        prev_block_sum = block_sum
        block_sum = block_sum + int(block_size)
        if cds_rel_start >= prev_block_sum and cds_rel_start < block_sum:
            exon_cds_start = i
            cds_coord_start = trans_start + int(block_start_list[i]) + cds_rel_start - prev_block_sum
        if cds_rel_end >= prev_block_sum and cds_rel_end <= block_sum:
            exon_cds_end = i
            cds_coord_end = trans_start + int(block_start_list[i]) + cds_rel_end - prev_block_sum
    
    # look into why sometimes cds_coord_end is 0 !!!!!!
    if cds_coord_end == 0:
        print("")
        
    #Check to see length is the same as sequence file
    if block_sum != trans_len_dict[trans_id]:
        print("Mismatching lengths!")
        print(trans_id + "\t" + str(block_sum) + "\t" + str(trans_len_dict[trans_id]))
        sys.exit()
        
    exon_nums = len(block_list)
    
    nmd_flag = "prot_ok"
    # NMD candidates, stop codon not on last exon
    if strand == "+":
        if exon_nums > exon_cds_end + 1:
            nmd_flag = str( exon_nums - (exon_cds_end + 1))
            nmd_flag = "NMD" + nmd_flag
    elif strand == "-":
        if exon_cds_start > 0:
            nmd_flag = str(exon_cds_start)
            nmd_flag = "NMD" + nmd_flag
    
    line_split[6] = str(cds_coord_start)
    line_split[7] = str(cds_coord_end)
    
    degrade_flag = "full_length"
    if prot_start == 0:
        degrade_flag = "5prime_degrade"
    

    line_split[3] = line_split[3] + ";" + prot_id + ";" + degrade_flag + ";" + match_flag + ";" + nmd_flag + ";" + frame 
    
    
    outline = "\t".join(line_split)
    outfile.write(outline)
    outfile.write("\n")
    
    
    
    
    
    
    
    
    
    
    
    
    