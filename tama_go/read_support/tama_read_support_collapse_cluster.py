
import re
import sys
import time

#
#This script gets all read support for TAMA collapse transcripts from clustering output
#


print("USAGE: python tama_collapse_trans_support_cluster.py prefix_trans_read.bed cluster_file output_file")

print("opening collapse")
collapse_file = sys.argv[1]
collapse_file_contents = open(collapse_file).read().rstrip("\n").split("\n")

print("opening cluster")
cluster_file = sys.argv[2]
cluster_file_contents = open(cluster_file).read().rstrip("\n").split("\n")

outfile_name = sys.argv[3]
outfile = open(outfile_name,"w")


cluster_read_dict = {} #  cluster_read_dict[cluster][read_id] = 1

cluster_type_flag = "na"

line = cluster_file_contents[0]
bed_split = line.split("\t")
if line.split(",") >= 2 and len(bed_split) != 12:
    cluster_type_flag = "v1"
    print("Interpretting as Iso-Seq1 cluster")
elif line.split() == 2:
    cluster_type_flag = "v3"
    print("Interpretting as Iso-Seq3 cluster")
elif len(bed_split) == 12:
    if bed_split[3].startswith("G") and len(bed_split[3].split(".")) == 2:
        cluster_type_flag = "no_cluster"
        print("Interpretting as no cluster")

if cluster_type_flag == "na":
    print("Error: Cannot understand file type.")
    sys.exit()


print("Going through cluster")
for line in cluster_file_contents:
    # iso-seq3 cluster file: m54041_180706_222522.unpolished.primer_5p--primer_3p.cluster
    #from to
    #m54041_180706_222522/32637046/ccs m54041_180706_222522/32637046/ccs
    #m54041_180706_222522/56099655/ccs m54041_180706_222522/32637046/ccs
    
    #iso-seq1 cluster file: cluster_report.csv
    #cluster_id,read_id,read_type
    #c0,m54041_170613_200235/46924722/6127_30_CCS,FL

    # 1       5332    5907    G1.1;15_25_c122441/1/591        40      +       5332    5907    255,0,0 1       575     0


    if line.startswith("cluster_id"):
        continue
    elif line.startswith("from"):
        continue

    if cluster_type_flag == "v1":

        line_split = line.split(",")
    
        cluster_id = line_split[0]
        read_id = line_split[1]
        
    elif cluster_type_flag == "v3":
        line_split = line.split()
    
        cluster_id = line_split[1]
        read_id = line_split[0]
    elif cluster_type_flag == "no_cluster":
        line_split = line.split("\t")

        cluster_id = line_split[3].split(";")[1]
        read_id = line_split[3].split(";")[1]
    else:
        print("Error: Cluster version is not recognized. Please check cluster file.")
        sys.exit()
    
    if cluster_id not in cluster_read_dict:
        cluster_read_dict[cluster_id] = {}
    
    cluster_read_dict[cluster_id][read_id] = 1


trans_cluster_dict = {} # trans_cluster_dict[trans id][cluster id] = 1
trans_list = []
gene_list = []

gene_trans_dict = {} # gene_trans_dict[gene id] = list of trans id


outline = "\t".join(["gene_id","trans_id","gene_num_reads","trans_num_reads","cluster_line"])

outfile.write(outline)
outfile.write("\n")

for line in collapse_file_contents:
    line_split = line.split("\t")

    #1       5332    5907    G1.1;15_25_c122441/1/591        40      +       5332    5907    255,0,0 1       575     0
    #1       8599    9314    G2.1;1_14_c58121/1/714  40      +       8599    9314    255,0,0 1       715     0
    #1       8713    9314    G2.1;15_25_c122146/3/601        40      +       8713    9314    255,0,0 1       601     0

    id_line = line_split[3]
    id_split = id_line.split(";")
    trans_id = id_split[0]

    if cluster_type_flag == "v1":
        cluster_id = id_split[1].split("/")[0]
    elif cluster_type_flag == "v3":
        cluster_id = id_split[1]
    elif cluster_type_flag == "no_cluster":
        cluster_id = id_split[1]
    
    gene_id = trans_id.split(".")[0]
    
    if gene_id not in gene_trans_dict:
        gene_trans_dict[gene_id] = []
        gene_list.append(gene_id)
    
    if trans_id not in trans_cluster_dict:
        trans_cluster_dict[trans_id] = {}
        trans_list.append(trans_id)
        gene_trans_dict[gene_id].append(trans_id)
    
    trans_cluster_dict[trans_id][cluster_id] = 1
    
    
for gene_id in gene_list:
    
    gene_trans_list = gene_trans_dict[gene_id]
    
    gene_num_reads = 0
    
    # get num reads for each gene
    for trans_id in gene_trans_list:
        trans_cluster_list = list(trans_cluster_dict[trans_id].keys())
        
        for cluster_id in trans_cluster_list:
            cluster_read_list = list(cluster_read_dict[cluster_id].keys())
            
            num_reads = len(cluster_read_list)
            
            gene_num_reads = gene_num_reads + num_reads

    # now do each transcript
    for trans_id in gene_trans_list:
        trans_cluster_list = list(trans_cluster_dict[trans_id].keys())
        
        trans_num_reads = 0
        
        cluster_list = []
        
        for cluster_id in trans_cluster_list:
            cluster_read_list = list(cluster_read_dict[cluster_id].keys())
            
            num_reads = len(cluster_read_list)
            
            trans_num_reads = trans_num_reads + num_reads
            
            read_line = ",".join(cluster_read_list)
            cluster_read_line = ":".join([cluster_id,read_line])
            cluster_list.append(cluster_read_line)
        
        cluster_line = ";".join(cluster_list)
        
        outline = "\t".join([gene_id,trans_id,str(gene_num_reads),str(trans_num_reads),cluster_line])

        outfile.write(outline)
        outfile.write("\n")
        
        
        
        









