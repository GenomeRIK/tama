
import re
import sys
import time
# from Bio import SeqIO


#
# This script finds all read support for transcripts in tama merge output
#

print("opening merge")
merge_file = sys.argv[1]
merge_file_contents = open(merge_file).read().rstrip("\n").split("\n")

print("opening filelist")
filelist_file = sys.argv[2]
filelist_file_contents = open(filelist_file).read().rstrip("\n").split("\n")


outfile_name = sys.argv[3]
outfile = open(outfile_name,"w")

trans_read_dict = {} #  trans_read_dict[source][trans_id] =  read count

gene_trans_read_dict = {} # gene_trans_read_dict[source][gene id][trans id][read id] = 1
prefix_list = []


for fline in filelist_file_contents:
    fline_split = fline.split("\t")

    #read_support_tama_collapse_oldbrain_nocap.txt   oldbrain        /exports/cmvm/eddie/eb/groups/burt_chicken_isoseq/dev_rk/4_tama_collapse_all/1_oldbrain/
    #read_support_tama_collapse_embryo.txt   embryo  /exports/cmvm/eddie/eb/groups/burt_chicken_isoseq/dev_rk/4_tama_collapse_all/2_embryo/


    filename = fline_split[0]
    prefix = fline_split[1]

    fpath = fline_split[2]

    prefix_list.append(prefix)

    support_file = fpath + filename

    if prefix in trans_read_dict:
        print("Error with duplicate prefix")
        sys.exit()

    trans_read_dict[prefix] = {}
    gene_trans_read_dict[prefix] = {}

    # filename check
    support_file_contents = open(support_file).read().rstrip("\n").split("\n")

    for line in support_file_contents:
        # gene_id trans_id        gene_num_reads  trans_num_reads cluster_line
        # G1      G1.1    3       3       rec_c6109:m54041_161210_051053/74842291/1157_29_CCS,m54041_161209_145041/37749615/30_1156_CCS,m54041_161210_051053/56034011/1128_31_CCS

        if line.startswith("gene_id"):
            continue

        line_split = line.split("\t")

        gene_id = line_split[0]
        trans_id = line_split[1]
        gene_num_reads = line_split[2]
        trans_num_reads = line_split[3]
        cluster_line = line_split[4]
        
        #collect all reads
        cluster_split = cluster_line.split(";")
        total_read_list = []
        for cluster_group in cluster_split:
            read_list = cluster_group.split(":")[1].split(",")
            total_read_list.extend(read_list)

        if trans_id in trans_read_dict[prefix]:
            print("Error with duplicate trans id")
            sys.exit()

        trans_read_dict[prefix][trans_id] = int(trans_num_reads)
        
        if gene_id not in gene_trans_read_dict[prefix]:
            gene_trans_read_dict[prefix][gene_id] = {}
        if trans_id not in gene_trans_read_dict[prefix][gene_id]:
            gene_trans_read_dict[prefix][gene_id][trans_id] = {}
        
        for this_read_id in total_read_list:
            gene_trans_read_dict[prefix][gene_id][trans_id][this_read_id] = 1
        

merge_read_dict = {}  # merge_read_dict[merge trans][prefix] = read count
merge_gene_trans_read_dict = {} # merge_gene_trans_read_dict[merge gene][merge trans][prefix][source trans][read id] = 1

merge_gene_trans_dict = {} # merge_gene_trans_dict[merge gene] = list of merge trans

merge_trans_dict = {}  # merge_trans_dict[merge trans][source trans] = 1

merge_trans_list = []
merge_gene_list = []

for line in merge_file_contents:
    line_split = line.split("\t")

    # 1       219     3261    G1.1;spleen_G1.1        40      +       219     3261    255,0,0 5       98,93,181,107,714       0,1457,1757,2132,2328

    id_line = line_split[3]
    id_split = id_line.split(";")
    trans_id = id_split[0]
    support_id = id_split[1]
    
    gene_id = trans_id.split(".")[0]

    num_exons = int(line_split[9])

    prefix = support_id.split("_")[0]

    collapse_id = support_id.split("_")[1]

    if trans_id not in merge_read_dict:
        merge_read_dict[trans_id] = {}

        for all_prefix in prefix_list:
            merge_read_dict[trans_id][all_prefix] = 0

        merge_trans_list.append(trans_id)

        merge_trans_dict[trans_id] = {}

    merge_trans_dict[trans_id][support_id] = 1

    read_count = trans_read_dict[prefix][collapse_id]

    merge_read_dict[trans_id][prefix] = merge_read_dict[trans_id][prefix] + int(read_count)
    
    # get full dict of read support
    
    if gene_id not in merge_gene_trans_read_dict:
        merge_gene_list.append(gene_id)
        merge_gene_trans_read_dict[gene_id] = {}
        merge_gene_trans_dict[gene_id] = []
    if trans_id not in merge_gene_trans_read_dict[gene_id]:
        merge_gene_trans_read_dict[gene_id][trans_id] = {}
        merge_gene_trans_dict[gene_id].append(trans_id)
        
    if prefix not in merge_gene_trans_read_dict[gene_id][trans_id]:
        merge_gene_trans_read_dict[gene_id][trans_id][prefix] = {}
    if collapse_id not in merge_gene_trans_read_dict[gene_id][trans_id][prefix]:
        merge_gene_trans_read_dict[gene_id][trans_id][prefix][collapse_id] = {}
    
    source_gene_id = collapse_id.split(".")[0]
    for support_read_id in gene_trans_read_dict[prefix][source_gene_id][collapse_id]:
        merge_gene_trans_read_dict[gene_id][trans_id][prefix][collapse_id][support_read_id] = 1



outline = "\t".join(["merge_gene_id","merge_trans_id", "gene_read_support", "trans_read_support","source_prefix", "source_trans_line","source_read_line"])

outfile.write(outline)
outfile.write("\n")

for merge_gene_id in merge_gene_list:
    
    total_read_count_gene = 0
    
    this_gene_read_dict = {} # this_gene_read_dict[read id] = 1
    
    #get gene read support
    for merge_trans_id in merge_gene_trans_read_dict[merge_gene_id]:
        for this_prefix in merge_gene_trans_read_dict[merge_gene_id][merge_trans_id]:
            for support_trans_id in merge_gene_trans_read_dict[merge_gene_id][merge_trans_id][this_prefix]:
                for support_read_id in merge_gene_trans_read_dict[merge_gene_id][merge_trans_id][this_prefix][support_trans_id]:
                    this_gene_read_dict[support_read_id] = 1
    
    total_read_count_gene = len(list(this_gene_read_dict.keys()))
    
    #get transcript read support
    for merge_trans_id in merge_gene_trans_dict[merge_gene_id]:
        total_read_count_trans = 0
        this_trans_read_dict = {} # this_trans_read_dict[read id] = 1
        this_prefix_list = []
        this_support_trans_list = []
        
        this_all_support_read_list = []
        
        for this_prefix in merge_gene_trans_read_dict[merge_gene_id][merge_trans_id]:
            this_prefix_list.append(this_prefix)
            for support_trans_id in merge_gene_trans_read_dict[merge_gene_id][merge_trans_id][this_prefix]:
                
                prefix_support_trans_id = this_prefix + "_" + support_trans_id
                this_support_trans_list.append(prefix_support_trans_id)
                
                this_trans_support_read_list = list(merge_gene_trans_read_dict[merge_gene_id][merge_trans_id][this_prefix][support_trans_id].keys())
                this_trans_support_read_line = ",".join(this_trans_support_read_list)
                this_all_support_read_list.append(this_trans_support_read_line)
                
                for support_read_id in merge_gene_trans_read_dict[merge_gene_id][merge_trans_id][this_prefix][support_trans_id]:
                    this_trans_read_dict[support_read_id] = 1
        
        total_read_count_trans = len(list(this_trans_read_dict.keys()))
        this_all_support_read_line = ";".join(this_all_support_read_list)
        this_prefix_line = ",".join(this_prefix_list)
        this_support_trans_line = ",".join(this_support_trans_list)
        
        outline = "\t".join([merge_gene_id,merge_trans_id,str(total_read_count_gene),str(total_read_count_trans),this_prefix_line,this_support_trans_line,this_all_support_read_line])
        outfile.write(outline)
        outfile.write("\n")

