import re
import sys
import time
from Bio import SeqIO
import argparse


#
# This script uses data from the blastp parse file and the pacbio annotation to assign the locations of the UTR/CDS regions to the bed file
# Use sequence file to double check the sequence lengths
# Last updated 2020/12/11
#


ap = argparse.ArgumentParser(description='This script uses data from the blastp parse file and the original annotation to assign the locations of the UTR/CDS regions to the bed file')

ap.add_argument('-p', type=str, nargs=1, help='Blastp parse file (required)')
ap.add_argument('-a', type=str, nargs=1, help='Annotation bed file (required)')
ap.add_argument('-f', type=str, nargs=1, help='Fasta for annotation file (required)')
ap.add_argument('-o', type=str, nargs=1, help='Output file name (required)')
ap.add_argument('-s', type=str, nargs=1, help='Include stop codon in CDS region (include_stop), default is to remove stop codon from CDS region')
ap.add_argument('-d', type=str, nargs=1, help='Distance from last splice junction to call NMD (default 50bp)')

opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0

if not opts.p:
    print("Blastp parse file missing")
    missing_arg_flag = 1
if not opts.a:
    print("Annotation bed file missing")
    missing_arg_flag = 1
if not opts.f:
    print("Fasta file missing")
    missing_arg_flag = 1
if not opts.o:
    print("output name missing")
    missing_arg_flag = 1

if missing_arg_flag == 1:
    print("Please try again with complete arguments")


if not opts.d:
    print("Using default SJ distance of 50bp")
    sj_dist_threshold = 50
else:
    sj_dist_threshold = int(opts.d[0])

parse_file = opts.p[0]
pbri_file = opts.a[0]
fasta_file = opts.f[0]
outfile_name = opts.o[0]

if not opts.s:
    include_stop_flag = "no_stop_codon"
else:
    include_stop_flag = opts.s[0]


print("opening blastp parse file")
parse_file_contents = open(parse_file).read().rstrip("\n").split("\n")

print("opening anno file")
pbri_file_contents = open(pbri_file).read().rstrip("\n").split("\n")

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
    trans_id = trans_id.split("(")[0]
    trans_seq = str(seq_record.seq)
    trans_seq = trans_seq.upper()
    trans_length = len(trans_seq)
    trans_len_dict[trans_id] = trans_length


#check for issues with matching transcript ID's
matching_trans_id_count = 0
for line in pbri_file_contents:
    line_split = line.split("\t")
    trans_id = line_split[3]
    
    if trans_id not in trans_dict:
        continue

    matching_trans_id_count += 1

if matching_trans_id_count == 0:
    print("Error with transcript ID's not matching.")
    print("Check the ID's in all the input files to see if they have matching formats.")
    sys.exit()

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

        if include_stop_flag == "include_stop":
            cds_rel_end = nuc_end + pos_adjust + 1
        else:
            cds_rel_end = nuc_end + pos_adjust - 2

    elif strand == "-":
        #cds_rel_start = trans_len_dict[trans_id] - (prot_end * 3 + pos_adjust + 3)
        #cds_rel_end = trans_len_dict[trans_id] - (prot_start * 3 + pos_adjust + 1)
        
        if include_stop_flag == "include_stop":
            cds_rel_start = trans_len_dict[trans_id] - (nuc_end + pos_adjust + 1)
        else:
            cds_rel_start = trans_len_dict[trans_id] - (nuc_end + pos_adjust - 2)

        cds_rel_end = trans_len_dict[trans_id] - (nuc_start + pos_adjust )
        
        
    cds_coord_start = 0
    cds_coord_end = 0
        
    #get sequence length by summing up blocks
    block_sum = 0
    exon_cds_start = 0
    exon_cds_end = 0

    exon_start_list = []
    exon_end_list = []
    
    for i, block_size in enumerate(block_list):
        prev_block_sum = block_sum
        block_sum = block_sum + int(block_size)


        exon_start = trans_start + int(block_start_list[i])
        exon_start_list.append(exon_start)

        exon_end = exon_start + int(block_size)
        exon_end_list.append(exon_end)


        if cds_rel_start >= prev_block_sum and cds_rel_start < block_sum:
            exon_cds_start = i
            cds_coord_start = trans_start + int(block_start_list[i]) + cds_rel_start - prev_block_sum
        if cds_rel_end >= prev_block_sum and cds_rel_end <= block_sum:
            exon_cds_end = i
            cds_coord_end = trans_start + int(block_start_list[i]) + cds_rel_end - prev_block_sum
    
    # look into why sometimes cds_coord_end is 0 !!!!!!
    if cds_coord_end == 0:
        print("Error with cds_coord_end == 0" +trans_id)
        sys.exit()
        
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

            stop_dist_sj = exon_start_list[-1] - cds_coord_end

            if stop_dist_sj > sj_dist_threshold:
                nmd_flag = str( exon_nums - (exon_cds_end + 1))
                nmd_flag = "NMD" + nmd_flag

    elif strand == "-":
        if exon_cds_start > 0:

            stop_dist_sj = cds_coord_start - exon_end_list[0]

            if stop_dist_sj > sj_dist_threshold:
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
    
    
    
    
    
    
    
    
    
    
    
    
    
