import re
import sys
import time


import os
import argparse


#This script takes the tama collapse files from a nocap run and a capped run to calculate degradation signature

ap = argparse.ArgumentParser(description='This script takes the tama collapse trans_read.bed files from a nocap run and a capped run to calculate degradation signature')

ap.add_argument('-c', type=str, nargs=1, help='Bed file from capped TAMA Collapse run (required)')
ap.add_argument('-nc', type=str, nargs=1, help='Bed file from no_cap TAMA Collapse run (required)')
ap.add_argument('-o', type=str, nargs=1, help='Output file name (required)')


opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0


if not opts.c:
    print("Capped bed file missing")
    missing_arg_flag = 1
if not opts.nc:
    print("No cap bed file missing")
    missing_arg_flag = 1
if not opts.o:
    print("output name missing")
    missing_arg_flag = 1

if missing_arg_flag == 1:
    print("Please try again with complete arguments")

capped_file = opts.c[0]
nocap_file = opts.nc[0]

outfile_name = opts.o[0]

print("opening capped file")
#bed_file = sys.argv[1]
capped_file_contents = open(capped_file).read().rstrip("\n").split("\n")

print("opening nocap file")
#bed_file = sys.argv[1]
nocap_file_contents = open(nocap_file).read().rstrip("\n").split("\n")

#outfile_name = sys.argv[3]
outfile = open(outfile_name,"w")

capped_gene_exon_dict = {} # capped_gene_exon_dict[gene_id] = max exons
capped_gene_read_dict = {} # capped_gene_read_dict[gene_id][read_id] = 1

capped_gene_trans_dict = {} # capped_gene_trans_dict[gene_id][trans_id] = 1

capped_gene_list = []

capped_trans_exon_dict = {} # capped_trans_exon_dict[trans_id] = num exons

for line in capped_file_contents:
    line_split = line.split("\t")

    id_line = line_split[3]

    id_split = id_line.split(";")

    trans_id = id_split[0]
    read_id = id_split[1]

    gene_id = trans_id.split(".")[0]

    num_exons = int(line_split[9])

    if gene_id not in capped_gene_exon_dict:
        capped_gene_exon_dict[gene_id] = 0
        capped_gene_read_dict[gene_id] = {}
        capped_gene_trans_dict[gene_id] = {}

        capped_gene_list.append(gene_id)


    if capped_gene_exon_dict[gene_id] < num_exons:
        capped_gene_exon_dict[gene_id] = num_exons

    capped_gene_read_dict[gene_id][read_id] = 1

    capped_gene_trans_dict[gene_id][trans_id] = 1

    if trans_id not in capped_trans_exon_dict:
        capped_trans_exon_dict[trans_id] = 0

    if capped_trans_exon_dict[trans_id] < num_exons:
        capped_trans_exon_dict[trans_id] = num_exons



nocap_gene_exon_dict = {} # nocap_gene_exon_dict[gene_id] = max exons
nocap_gene_read_dict = {} # nocap_gene_read_dict[gene_id][read_id] = 1

nocap_gene_trans_dict = {} # nocap_gene_trans_dict[gene_id][trans_id] = 1

nocap_gene_list = []

nocap_trans_exon_dict = {} # nocap_trans_exon_dict[trans_id] = num exons

for line in nocap_file_contents:
    line_split = line.split("\t")

    id_line = line_split[3]

    id_split = id_line.split(";")

    trans_id = id_split[0]
    read_id = id_split[1]

    gene_id = trans_id.split(".")[0]

    num_exons = int(line_split[9])

    if gene_id not in nocap_gene_exon_dict:
        nocap_gene_exon_dict[gene_id] = 0
        nocap_gene_read_dict[gene_id] = {}
        nocap_gene_trans_dict[gene_id] = {}

        nocap_gene_list.append(gene_id)


    if nocap_gene_exon_dict[gene_id] < num_exons:
        nocap_gene_exon_dict[gene_id] = num_exons

    nocap_gene_read_dict[gene_id][read_id] = 1

    nocap_gene_trans_dict[gene_id][trans_id] = 1

    if trans_id not in nocap_trans_exon_dict:
        nocap_trans_exon_dict[trans_id] = 0

    if nocap_trans_exon_dict[trans_id] < num_exons:
        nocap_trans_exon_dict[trans_id] = num_exons


capped_trans_count = 0

capped_single_exon_gene_count = 0
capped_multi_exon_gene_count = 0

capped_single_exon_single_read_gene_count = 0
capped_multiexon_single_read_gene_count = 0

for gene_id in capped_gene_list:

    if capped_gene_exon_dict[gene_id] < 2:
        capped_single_exon_gene_count += 1
        num_reads = len(list(capped_gene_read_dict[gene_id].keys()))
        if num_reads == 1:
            capped_single_exon_single_read_gene_count += 1
        continue

    capped_multi_exon_gene_count += 1

    num_reads = len(list(capped_gene_read_dict[gene_id].keys()))
    if num_reads < 2:
        capped_multiexon_single_read_gene_count += 1
        continue

    num_trans = len(list(capped_gene_trans_dict[gene_id].keys()))

    capped_trans_count = capped_trans_count + num_trans


nocap_trans_count = 0
nocap_single_exon_gene_count = 0

nocap_multiexon_single_read_gene_count = 0
nocap_single_exon_single_read_gene_count = 0

nocap_multi_exon_gene_count = 0

for gene_id in nocap_gene_list:

    if nocap_gene_exon_dict[gene_id] < 2:
        nocap_single_exon_gene_count += 1

        num_reads = len(list(nocap_gene_read_dict[gene_id].keys()))
        if num_reads == 1:
            nocap_single_exon_single_read_gene_count += 1
        continue

    nocap_multi_exon_gene_count += 1

    num_reads = len(list(nocap_gene_read_dict[gene_id].keys()))
    if num_reads < 2:
        nocap_multiexon_single_read_gene_count += 1
        continue

    num_trans = len(list(nocap_gene_trans_dict[gene_id].keys()))

    nocap_trans_count = nocap_trans_count + num_trans



deg_sig = (float(capped_trans_count) - float(nocap_trans_count))/float(capped_trans_count)

# print out summary

outline = "Degradation Signature = " + str(deg_sig)
outfile.write(outline)
outfile.write("\n")

## trans level

outline = "Capped multi-exon, multi-read, transcript count = " + str(capped_trans_count)
outfile.write(outline)
outfile.write("\n")

outline = "No-cap multi-exon, multi-read, transcript count = " + str(nocap_trans_count)
outfile.write(outline)
outfile.write("\n")


total_capped_num_trans = len(list(capped_trans_exon_dict.keys()))
total_nocap_num_trans = len(list(nocap_trans_exon_dict.keys()))

outline = "Capped total transcript count = " + str(total_capped_num_trans)
outfile.write(outline)
outfile.write("\n")

outline = "No-cap total transcript count = " + str(total_nocap_num_trans)
outfile.write(outline)
outfile.write("\n")


capped_single_exon_trans_count = 0
capped_multi_exon_trans_count = 0

for trans_id in capped_trans_exon_dict:
    if capped_trans_exon_dict[trans_id] == 1:
        capped_single_exon_trans_count += 1
    else:
        capped_multi_exon_trans_count += 1

nocap_single_exon_trans_count = 0
nocap_multi_exon_trans_count = 0

for trans_id in nocap_trans_exon_dict:
    if nocap_trans_exon_dict[trans_id] == 1:
        nocap_single_exon_trans_count += 1
    else:
        nocap_multi_exon_trans_count += 1


outline = "Capped single exon trans count = " + str(capped_single_exon_trans_count)
outfile.write(outline)
outfile.write("\n")

outline = "No-cap single exon trans count = " + str(nocap_single_exon_trans_count)
outfile.write(outline)
outfile.write("\n")

outline = "Capped multi exon trans count = " + str(capped_multi_exon_trans_count)
outfile.write(outline)
outfile.write("\n")

outline = "No-cap multi exon trans count = " + str(nocap_multi_exon_trans_count)
outfile.write(outline)
outfile.write("\n")



# gene level

capped_gene_count = len(capped_gene_list)
nocap_gene_count = len(nocap_gene_list)

outline = "Capped total gene count = " + str(capped_gene_count)
outfile.write(outline)
outfile.write("\n")

outline = "No-cap total gene count = " + str(nocap_gene_count)
outfile.write(outline)
outfile.write("\n")


outline = "Capped single exon gene count = " + str(capped_single_exon_gene_count)
outfile.write(outline)
outfile.write("\n")

outline = "No-cap single exon gene count = " + str(nocap_single_exon_gene_count)
outfile.write(outline)
outfile.write("\n")

outline = "Capped multi exon gene count = " + str(capped_multi_exon_gene_count)
outfile.write(outline)
outfile.write("\n")

outline = "No-cap multi exon gene count = " + str(nocap_multi_exon_gene_count)
outfile.write(outline)
outfile.write("\n")



outline = "Capped single exon single read gene count = " + str(capped_single_exon_single_read_gene_count)
outfile.write(outline)
outfile.write("\n")

outline = "No-cap single exon single read gene count = " + str(nocap_single_exon_single_read_gene_count)
outfile.write(outline)
outfile.write("\n")


outline = "Capped multi exon single read gene count = " + str(capped_multiexon_single_read_gene_count)
outfile.write(outline)
outfile.write("\n")

outline = "No-cap multi exon single read gene count = " + str(nocap_multiexon_single_read_gene_count)
outfile.write(outline)
outfile.write("\n")











