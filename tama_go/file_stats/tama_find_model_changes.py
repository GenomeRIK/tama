import re
import sys
import time


import os
import argparse


# This script looks at read support files and finds reads that have mapped to different genes

ap = argparse.ArgumentParser(description='This script looks at read support files and finds reads that have mapped to different genes')

ap.add_argument('-b', type=str, nargs=1, help='Annotation bed file')
ap.add_argument('-r', type=str, nargs=1, help='Read support file')
ap.add_argument('-o', type=str, nargs=1, help='Output file prefix')

ap.add_argument('-ref', type=str, nargs=1, help='Reference source')
ap.add_argument('-alt', type=str, nargs=1, help='Alternative source')

opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0

if not opts.b:
    print("Annotation bed file missing")
    missing_arg_flag = 1
if not opts.r:
    print("Read Support file missing")
    missing_arg_flag = 1
if not opts.o:
    print("output file prefix missing")
    missing_arg_flag = 1

if missing_arg_flag == 1:
    print("Please try again with complete arguments")


ref_source = "NA"
alt_source = "NA"

if not opts.ref:
    print("Default is no ref source")
    ref_source = "NA"
else:
    ref_source = opts.ref[0]

if not opts.alt:
    print("Default is no alt source")
    alt_source = "NA"
else:
    alt_source = opts.alt[0]


bed_file = opts.b[0]
read_file = opts.r[0]

outfile_prefix = opts.o[0]
outfile_name_gene = outfile_prefix + "_diff_genes.txt"
outfile_gene = open(outfile_name_gene,"w")

outfile_name_trans = outfile_prefix + "_diff_trans.txt"
outfile_trans = open(outfile_name_trans,"w")

outfile_name_report = outfile_prefix + "_diff_report.txt"
outfile_report = open(outfile_name_report,"w")

outfile_prefix = opts.o[0]
outfile_name_onegene = outfile_prefix + "_diff_one_source_genes.txt"
outfile_onegene = open(outfile_name_onegene,"w")

outfile_name_onetrans = outfile_prefix + "_diff_one_source_trans.txt"
outfile_onetrans = open(outfile_name_onetrans,"w")


print("opening bed file")
#bed_file = sys.argv[1]
bed_file_contents = open(bed_file).read().rstrip("\n").split("\n")

print("opening read support file")
#bed_file = sys.argv[1]
read_file_contents = open(read_file).read().rstrip("\n").split("\n")


gene_trans_dict = {} # gene_trans_dict[gene_id][trans_id] = trans_split

gene_pos_dict = {} # gene_pos_dict[gene_id] = chrom start end

for line in bed_file_contents:
    line_split = line.split("\t")

    id_line = line_split[3]

    id_split = id_line.split(";")

    gene_id = id_split[0]
    trans_id = id_split[1]

    if gene_id not in gene_trans_dict:
        gene_trans_dict[gene_id] = {}

    gene_trans_dict[gene_id][trans_id] = line_split

    num_exons = int(line_split[9])
    chrom = line_split[0]
    t_start = int(line_split[1])
    t_end = int(line_split[2])
    
    if gene_id not in gene_pos_dict:
        gene_pos_dict[gene_id] = []
        gene_pos_dict[gene_id].append(chrom)
        gene_pos_dict[gene_id].append(t_start)
        gene_pos_dict[gene_id].append(t_end)
    else:
        if gene_pos_dict[gene_id][0] != chrom:
            print("Error with mismatch in chrom")
            sys.exit()
        
        if t_start < gene_pos_dict[gene_id][1] :
            gene_pos_dict[gene_id][1]  = t_start
        if t_end > gene_pos_dict[gene_id][2] :
            gene_pos_dict[gene_id][2] = t_end



read_source_gene_dict = {} # read_source_gene_dict[read_id][source][gene_id] = 1
read_gene_dict = {} # read_gene_dict[read_id][gene_id] = 1

read_source_trans_dict = {} # read_source_trans_dict[read_id][source][trans_id] = 1
read_trans_dict = {} # read_trans_dict[read_id][trans_id] = 1

gene_source_dict = {} # gene_source_dict[gene_id][source] = 1
trans_source_dict = {} # trans_source_dict[trans_id][source] = 1

source_gene_merge_dict = {} # source_gene_merge_dict[source][gene_id] = merge_gene_id
source_trans_merge_dict = {} # source_trans_merge_dict[source][trans_id] = merge_trans_id


source_dict = {}  # source_dict[source] = 1

all_read_list = []

for line in read_file_contents:
    #"merge_gene_id","merge_trans_id","gene_read_count","trans_read_count","source_line","support_line"

    if line.startswith("merge_gene_id"):
        continue
        
    line_split = line.split("\t")

    merge_gene_id = line_split[0]
    merge_trans_id = line_split[1]
    gene_read_count = line_split[2]
    num_read_support = line_split[3]
    this_source_line = line_split[4]
    this_read_line = line_split[5]

    if merge_gene_id not in gene_source_dict:
        gene_source_dict[merge_gene_id] = {}

    if merge_trans_id not in trans_source_dict:
        trans_source_dict[merge_trans_id] = {}


    this_source_split = this_source_line.split(",")
    for this_source in this_source_split:
        gene_source_dict[merge_gene_id][this_source] = 1
        trans_source_dict[merge_trans_id][this_source] = 1

        source_dict[this_source] = 1

    this_read_line_split = this_read_line.split(";")

    for this_source_read_line in this_read_line_split:
        this_source_name = this_source_read_line.split(":")[0]

        this_read_list = this_source_read_line.split(":")[1].split(",")

        for read_id in this_read_list:
            if read_id not in read_source_gene_dict:
                read_source_gene_dict[read_id] = {}
                all_read_list.append(read_id)
            
            if read_id not in read_gene_dict:
                read_gene_dict[read_id] = {}
            
            if read_id not in read_source_trans_dict: 
                read_source_trans_dict[read_id] = {}
                read_trans_dict[read_id] = {}
            
            if this_source_name not in read_source_gene_dict[read_id]:
                read_source_gene_dict[read_id][this_source_name] = {}
            
            if this_source_name not in read_source_trans_dict[read_id]:
                read_source_trans_dict[read_id][this_source_name] = {}
            
            read_source_gene_dict[read_id][this_source_name][merge_gene_id] = 1
            read_source_trans_dict[read_id][this_source_name][merge_trans_id] = 1
            
            read_gene_dict[read_id][merge_gene_id] = 1
            read_trans_dict[read_id][merge_trans_id] = 1


#write out headers
outline =  "\t".join(['read_id','num_genes','all_gene_line','all_pos_line','all_trans_line'])
outfile_gene.write(outline)
outfile_gene.write("\n")

outline =  "\t".join(['read_id','alt_trans_diff_count','alt_diff_trans_id_list_line','alt_trans_id_list_line','ref_trans_id_list_line'])
outfile_trans.write(outline)
outfile_trans.write("\n")

read_diff_gene_dict = {} # read_diff_gene_dict[read_id] = 1
read_diff_trans_dict = {} # read_diff_trans_dict[read_id] = 1

# use to find diff genes and trans by source
diff_gene_source_dict = {} # diff_gene_source_dict[source][gene_id] = 1
diff_trans_source_dict = {} # diff_trans_source_dict[source][trans_id] = 1

# use to find total number of diff genes and trans
merge_diff_gene_dict = {} # merge_diff_gene_dict[gene_id] = 1
merge_diff_trans_dict = {} # merge_diff_trans_dict[trans_id] = 1


all_source_list = list(source_dict.keys())
all_source_list.sort()

for this_source in all_source_list:
    diff_gene_source_dict[this_source] = {}
    diff_trans_source_dict[this_source] = {}


for read_id in all_read_list:
    
    # look for gene differences
    num_genes = len(list(read_gene_dict[read_id].keys()))
    
    if num_genes > 1:
        gene_line_list = []
        pos_line_list = []
        
        trans_line_list = []
        
        for source_name in read_source_gene_dict[read_id]:
            for merge_gene_id in read_source_gene_dict[read_id][source_name]:

                chrom = gene_pos_dict[merge_gene_id][0]
                g_start = gene_pos_dict[merge_gene_id][1]
                g_end = gene_pos_dict[merge_gene_id][2]
                
                pos_line = chrom + ":" + str(g_start) + "-" + str(g_end)
                
                pos_line_list.append(pos_line)
                
                gene_line = source_name + ":" + merge_gene_id
                
                gene_line_list.append(gene_line)

                diff_gene_source_dict[source_name][merge_gene_id] = 1
                merge_diff_gene_dict[merge_gene_id] = 1
        
        for source_name in read_source_trans_dict[read_id]:
            for trans_id in read_source_trans_dict[read_id][source_name]:
                trans_line = source_name + ":" + trans_id
                trans_line_list.append(trans_line)
                
        all_gene_line = ",".join(gene_line_list)
        all_pos_line = ",".join(pos_line_list)
        
        all_trans_line = ",".join(trans_line_list)
        
        read_diff_gene_dict[read_id] = 1
        
        outline =  "\t".join([read_id,str(num_genes),all_gene_line,all_pos_line,all_trans_line])
        outfile_gene.write(outline)
        outfile_gene.write("\n")

    # look for transcript differences ###############################################
    
    #skip trans assessment if it is already a diff gene read
    if read_id in read_diff_gene_dict:
        continue
    
    ref_trans_dict = {} # ref_trans_dict[trans_id] = 1
    
    for source_name in read_source_trans_dict[read_id]:
        
        if source_name == ref_source:
            for trans_id in read_source_trans_dict[read_id][source_name]:
                ref_trans_dict[trans_id] = 1

    diff_trans_flag = 0
    
    alt_trans_id_list = []
    alt_diff_trans_id_list = []
    alt_trans_diff_count = 0
    
    for source_name in read_source_trans_dict[read_id]:
        
        if source_name == alt_source:
            for trans_id in read_source_trans_dict[read_id][source_name]:
                
                alt_trans_id_list.append(trans_id)

                # read neads to have different modles but be in both sources
                if trans_id not in ref_trans_dict:
                    if ref_source in read_source_trans_dict[read_id]:
                        diff_trans_flag = 1
                        alt_trans_diff_count += 1
                        alt_diff_trans_id_list.append(trans_id)


                        # Add merge trans id for alt transcript
                        diff_trans_source_dict[source_name][trans_id] = 1
                        merge_diff_trans_dict[trans_id] = 1

    if diff_trans_flag == 0:
        continue

    read_diff_trans_dict[read_id] = 1

    alt_trans_id_list_line = ",".join(alt_trans_id_list)
    
    alt_diff_trans_id_list_line = ",".join(alt_diff_trans_id_list)
    
    ref_trans_id_list = []
    if ref_source in read_source_trans_dict[read_id]:
        for ref_trans_id in read_source_trans_dict[read_id][ref_source]:
            ref_trans_id_list.append(ref_trans_id)

            # Add merge trans id for ref transcript
            diff_trans_source_dict[ref_source][ref_trans_id] = 1
            merge_diff_trans_dict[ref_trans_id] = 1
            
            
    else:
        ref_trans_id_list.append("NA")
        #skip this read because it was discarded in ref isoseq dataset
        # polish may have rescued the read but that does not mean it was good
        continue
    
    ref_trans_id_list_line = ",".join(ref_trans_id_list)
    
    outline =  "\t".join([read_id,str(alt_trans_diff_count),alt_diff_trans_id_list_line,alt_trans_id_list_line,ref_trans_id_list_line])
    outfile_trans.write(outline)
    outfile_trans.write("\n")


################################################
# write out report file


num_diff_gene_reads = len(list(read_diff_gene_dict.keys()))
num_diff_trans_reads = len(list(read_diff_trans_dict.keys()))

outline = "num_diff_gene_reads: " + str(num_diff_gene_reads)
outfile_report.write(outline)
outfile_report.write("\n")

outline = "num_diff_trans_reads: " + str(num_diff_trans_reads)
outfile_report.write(outline)
outfile_report.write("\n")

############

num_merge_diff_gene =  len(list(merge_diff_gene_dict.keys()))
num_merge_diff_trans =  len(list(merge_diff_trans_dict.keys()))


outline = "num_merge_diff_gene: " + str(num_merge_diff_gene)
outfile_report.write(outline)
outfile_report.write("\n")

outline = "num_merge_diff_trans: " + str(num_merge_diff_trans)
outfile_report.write(outline)
outfile_report.write("\n")

############


for this_source in all_source_list:
    this_source_diff_genes = len(list(diff_gene_source_dict[this_source].keys()))
    this_source_diff_trans = len(list(diff_trans_source_dict[this_source].keys()))

    outline = "this_source_diff_genes "+ str(this_source) + ": " + str(this_source_diff_genes)
    outfile_report.write(outline)
    outfile_report.write("\n")

    outline = "this_source_diff_trans " + str(this_source) + ": " + str(this_source_diff_trans)
    outfile_report.write(outline)
    outfile_report.write("\n")


############

# find genes and transcripts only in one source

only_source_gene_dict = {} # only_source_gene_dict[source][gene_id] = 1
only_source_trans_dict = {} # only_source_trans_dict[source][trans_id] = 1


# write geadrer for one source gene and trans files

outline = "\t".join(['merge_source','merge_gene_id'])
outfile_onegene.write(outline)
outfile_onegene.write("\n")

outline = "\t".join(['merge_source','merge_trans_id'])
outfile_onetrans.write(outline)
outfile_onetrans.write("\n")

for this_source in all_source_list:
    only_source_gene_dict[this_source] = {}
    only_source_trans_dict[this_source] = {}

for merge_gene_id in merge_diff_gene_dict:
    merge_source_list = list(gene_source_dict[merge_gene_id].keys())

    num_sources = len(merge_source_list)

    if num_sources == 2:
        continue

    elif num_sources == 1:
        only_source_gene_dict[merge_source_list[0]][merge_gene_id] = 1

        outline = "\t".join([merge_source_list[0],merge_gene_id])
        outfile_onegene.write(outline)
        outfile_onegene.write("\n")

    else:
        print("Error with number of sources")
        sys.exit()


for merge_trans_id in merge_diff_trans_dict:
    merge_source_list = list(trans_source_dict[merge_trans_id].keys())

    num_sources = len(merge_source_list)

    if num_sources == 2:
        continue

    elif num_sources == 1:

        merge_gene_id = merge_trans_id.split(".")[0]
        # only count one source trans if it is not involved in a one source gene! Very important to rememeber this!!!!!!
        if merge_gene_id not in only_source_gene_dict[merge_source_list[0]]:
            only_source_trans_dict[merge_source_list[0]][merge_trans_id] = 1

            ########################
            # only count one source trans if it is not involved in a one source gene! Very important to rememeber this!!!!!!
            outline = "\t".join([merge_source_list[0], merge_trans_id])
            outfile_onetrans.write(outline)
            outfile_onetrans.write("\n")

    else:
        print("Error with number of sources")
        sys.exit()


total_one_source_genes_count = 0
total_one_source_trans_count = 0


for this_source in all_source_list:
    only_source_num_genes = len(list(only_source_gene_dict[this_source].keys()))
    only_source_num_trans = len(list(only_source_trans_dict[this_source].keys()))

    total_one_source_genes_count = total_one_source_genes_count + only_source_num_genes
    total_one_source_trans_count = total_one_source_trans_count + only_source_num_trans

    outline = "only_source_num_genes "+ str(this_source) + ": " + str(only_source_num_genes)
    outfile_report.write(outline)
    outfile_report.write("\n")

    outline = "only_source_num_trans " + str(this_source) + ": " + str(only_source_num_trans)
    outfile_report.write(outline)
    outfile_report.write("\n")


outline = "total_one_source_genes_count: " + str(total_one_source_genes_count)
outfile_report.write(outline)
outfile_report.write("\n")

outline = "total_one_source_trans_count: " + str(total_one_source_trans_count)
outfile_report.write(outline)
outfile_report.write("\n")





















