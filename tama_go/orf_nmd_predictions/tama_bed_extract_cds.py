import re
import sys
import time
from Bio import SeqIO

import os
import argparse


#This script takes a bed file with cds information and creates a bed file with only cds regions
# Use this to get CDS sequence using bedtools getfasta

ap = argparse.ArgumentParser(description='This script takes a bed file with cds information and creates a bed file with only cds regions')

ap.add_argument('-b', type=str, nargs=1, help='Bed file (required)')
ap.add_argument('-s', type=str, nargs=1, help='Stop codon include flag (required)')
ap.add_argument('-o', type=str, nargs=1, help='Output file name (required)')


opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0


if not opts.b:
    print("Annotation bed file missing")
    missing_arg_flag = 1
if not opts.s:
    print("Stop codon include flag missing")
    missing_arg_flag = 1
if not opts.o:
    print("output name missing")
    missing_arg_flag = 1

if missing_arg_flag == 1:
    print("Please try again with complete arguments")

bed_file = opts.b[0]
stop_codon_flag = opts.s[0]
outfile_name = opts.o[0]

print("opening bed file")
#bed_file = sys.argv[1]
bed_file_contents = open(bed_file).read().rstrip("\n").split("\n")


#stop_codon_flag = sys.argv[2] # either include_stop or exclude_stop

#outfile_name = sys.argv[3]
outfile = open(outfile_name,"w")



gene_source = "tama"
trans_source = "tama"

def calc_end(start,block_size):
    end = int(start) + int(block_size) - 1 # adjust for bed 0 base and gtf 1 base coords
    #end = int(start) + int(block_size)  #
    return str(end)

def calc_exon_start(t_start,e_start):
    coordinate_start = int(e_start) + int(t_start)
    return str(coordinate_start)

def convert_str_list_to_int(str_list):
    int_list = []
    for string in str_list:
        int_list.append(int(string))
    return int_list

def convert_int_list_to_str(int_list):
    str_list = []
    for integer in int_list:
        str_list.append(str(integer))
    return str_list


class Transcript:
    def __init__(self, bed_line):
        line_split = bed_line.split("\t")
        chrom = line_split[0]
        t_start = int(line_split[1])
        t_end = int(line_split[2])
        id_line = line_split[3]
        strand = line_split[5]
        num_exons = int(line_split[9])
        blocks = line_split[10]
        starts = line_split[11]

        cds_start = int(line_split[6])
        cds_end = int(line_split[7])

        id_split = id_line.split(";")
        gene_id = id_split[0]
        trans_id = id_split[1]
        prot_id = id_split[2]
        degrade_flag = id_split[3]
        match_flag = id_split[4]
        nmd_flag = id_split[5]


        self.id_list = id_split
        self.prot_id = prot_id
        self.degrade_flag = degrade_flag
        self.match_flag = match_flag
        self.nmd_flag = nmd_flag

        self.trans_id = trans_id
        self.gene_id = gene_id
        self.chrom = chrom
        self.t_start = str(t_start)
        self.t_end = str(t_end)
        self.strand = strand

        t_start_list = [int(t_start)] * int(num_exons)
        start_list = starts.split(",")
        start_list = filter(None, start_list)
        block_list = blocks.split(",")
        block_list = filter(None, block_list)

        self.bed_starts = starts
        self.bed_blocks = blocks

        # coordinate starts and ends
        self.start_list = map(calc_exon_start, t_start_list, start_list)
        self.end_list = map(calc_end, self.start_list, block_list)
        self.num_exons = num_exons

        self.cds_start = cds_start  #################################################
        self.cds_end = cds_end  #################################################
        
        # make position list

        self.trans_coord_list = []
        self.trans_coord_dict = {} # trans_coord_dict[coord] = index
        for i in xrange(int(self.num_exons)):

            e_index = i
            e_num = e_index + 1

            e_start = int(self.start_list[e_index])
            e_end = int(self.end_list[e_index])

            e_length = e_end - e_start

            for j in xrange(e_length + 1):
                pos_coord = e_start + j
                self.trans_coord_list.append(pos_coord)
                self.trans_coord_dict[pos_coord] = len(self.trans_coord_list) - 1

        
        # adjust cds end because of bed format 0 1 number method
        # correct for neg strand when cds goes to end
        if self.strand == "-" and self.cds_end == self.trans_coord_list[-1] + 1:
                self.cds_end_adj = self.trans_coord_list[-1]
                
        else: # this is the normal condition for CDS
            try:
                cds_end_index = self.trans_coord_dict[self.cds_end]
            except:
                print(self.trans_coord_list)
                print(len(self.trans_coord_list))
                print(self.cds_end)
                print(self.trans_id)
                print("error with self.cds_end")
                sys.exit()
                
            self.cds_end_adj = self.trans_coord_list[cds_end_index - 1]

        # this is deprecated because I do this check before calling this function
        # just keep it here for sanity
        no_cds_flag = "cds_exists"
        if self.cds_start == 0 and self.cds_end == 0:
            no_cds_flag = "no_cds"

        if no_cds_flag == "cds_exists":
            if self.strand == "+":
                if stop_codon_flag == "include_stop":

                    try:
                        cds_end_index = self.trans_coord_dict[self.cds_end_adj]
                    except:
                        print(self.trans_coord_list)
                        print(self.cds_end_adj)
                        print(self.trans_id)
                        print("error with self.cds_end_adj")
                        sys.exit()

                    stop_codon_index = cds_end_index + 3

                    #if stop_codon_index >= len(self.trans_coord_list):
                    #    stop_codon_index = len(self.trans_coord_list) - 1
                        
                        #print(self.trans_coord_list)
                        #print(stop_codon_index)
                        #print(len(self.trans_coord_list))
                        #print(self.trans_id)
                        #sys.exit()

                    self.cds_end_adj = self.trans_coord_list[stop_codon_index]


            if self.strand == "-":
                if stop_codon_flag == "include_stop":
                    cds_start_index = self.trans_coord_dict[self.cds_start]
                    stop_codon_index = cds_start_index - 3
                    self.cds_start = self.trans_coord_list[stop_codon_index]
                    
#        if self.trans_id == "G1594.5":
#            print(self.trans_coord_list)
#            print(stop_codon_index)
#            print(len(self.trans_coord_list))
#            print(self.trans_id)
#            print(self.cds_end_adj)
#            print(self.cds_end)
#            sys.exit()


    def cds_extract(self):

        if int(self.num_exons) != len(self.start_list):
            print("Error with number of exons")
            sys.exit()

        cds_e_start_list = []
        cds_e_end_list = []

        for i in xrange(int(self.num_exons)):

            e_index = i
            e_num = e_index + 1

            e_start = int(self.start_list[e_index])
            e_end = int(self.end_list[e_index])

            cds_e_start = e_start
            cds_e_end = e_end

            if self.cds_start >= e_start and self.cds_start <= e_end:
                cds_e_start = self.cds_start
                #######################


            if self.cds_end_adj >= e_start and self.cds_end_adj <= e_end:
                cds_e_end = self.cds_end_adj
                ######################

            #if self.cds_end_adj == e_start:
            #    continue

            if self.cds_start > e_end: # cds start is after this exon
                continue

            if self.cds_end_adj < e_start: # cds end is before this exon
                continue


            cds_e_start_list.append(cds_e_start)
            cds_e_end_list.append(cds_e_end)

        self.cds_start = cds_e_start_list[0]
        self.cds_end_adj = cds_e_end_list[-1]

        # convert to new bed line
        new_block_list = []
        new_starts_list = []
        for i in xrange(len(cds_e_start_list)):
            cds_e_start = cds_e_start_list[i]
            cds_e_end = cds_e_end_list[i]

            new_start = str(cds_e_start - self.cds_start)
            new_block = str(1 + cds_e_end - cds_e_start)

            new_starts_list.append(new_start)
            new_block_list.append(new_block)

        new_starts_line = ",".join(new_starts_list)
        new_blocks_line = ",".join(new_block_list)

        new_bed_list = line_split
        new_bed_list[10] = new_blocks_line
        new_bed_list[11] = new_starts_line

        new_bed_list[1] = str(self.cds_start)
        new_bed_list[2] = str(self.cds_end_adj + 1)

        new_bed_list[6] = str(self.cds_start)
        new_bed_list[7] = str(self.cds_end_adj + 1)

        new_num_exons = str(len(new_starts_list))

        new_bed_list[9] = new_num_exons

        new_id_list = ["cds"]
        new_id_list.extend(self.id_list)

        new_id_line = ";".join(new_id_list)

        new_bed_list[3] = new_id_line

        new_bed_line = "\t".join(new_bed_list)

        return new_bed_line


count = 0

print("Going through bed file")
for line in bed_file_contents:

    count += 1
    if count % 10000 == 0 :
        print(count)


    line_split = line.split("\t")

    chrom = line_split[0]
    t_start = line_split[1]
    t_end = line_split[2]
    id_line = line_split[3]
    strand = line_split[5]
    num_exons = line_split[9]
    block_sizes = line_split[10]
    block_starts = line_split[11]

    cds_start = line_split[6]
    cds_end = line_split[7]

    id_split = id_line.split(";")
    gene_id = id_split[0]
    trans_id = id_split[1]



    if cds_start != "0" and cds_end != "0":
        trans_obj = Transcript(line)

        new_bed_line = trans_obj.cds_extract()

        outfile.write(new_bed_line)
        outfile.write("\n")











