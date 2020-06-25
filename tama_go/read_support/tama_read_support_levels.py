import re
import sys
import time


import os
import argparse


# This script produces a read support file for identifying what models are supported by which reads

ap = argparse.ArgumentParser(description='This script produces a read support file for identifying what models are supported by which reads')


ap.add_argument('-f', type=str, nargs=1, help='Filelist file with pre-merge trans_read.bed file names')
ap.add_argument('-m', type=str, nargs=1, help='Merge.txt file from after merging. Use "no_merge" if there is no merge file.')
ap.add_argument('-o', type=str, nargs=1, help='Output file prefix')

ap.add_argument('-d', type=str, nargs=1, help='Ignore duplicate read name warning with -d dup_ok, default is to flag duplicates and terminate early.')

ap.add_argument('-mt', type=str, nargs=1, help='Merge type flag indicates the type of merge file used. Use -mt cupcake for cupcake file. Default is TAMA output.')


opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0



if not opts.m:
    print("Merge file missing")
    missing_arg_flag = 1
if not opts.f:
    print("Filelist file missing")
    missing_arg_flag = 1
if not opts.o:
    print("output file prefix missing")
    missing_arg_flag = 1

if missing_arg_flag == 1:
    print("Please try again with complete arguments")

duplicate_flag = "no_dup"

if not opts.d:
    print("Default is to flag duplicates and terminate early.")
    duplicate_flag = "no_dup"
else:
    duplicate_flag = opts.d[0]
    if duplicate_flag != "no_dup" and duplicate_flag != "dup_ok":
        print("Error with -d input. Should be either no_dup or dup_ok.")
        print("Terminating early.")
        sys.exit()

merge_type_flag = "tama"

if not opts.mt:
    print("Default is tama flag for TAMA merge file.")
    duplicate_flag = "tama"
else:
    merge_type_flag = opts.mt[0]
    if merge_type_flag != "tama" and merge_type_flag != "cupcake" and  merge_type_flag != "filter":
        print("Error with -mt input. Should be tama or cupcake or filter.")
        print("Terminating early.")
        sys.exit()

filelist_file = opts.f[0]
merge_file = opts.m[0]

outfile_prefix = opts.o[0]
outfile_name = outfile_prefix + "_read_support.txt"
outfile = open(outfile_name,"w")


print("opening filelist file")
#bed_file = sys.argv[1]
filelist_file_contents = open(filelist_file).read().rstrip("\n").split("\n")


if merge_file != "no_merge":
    print("opening merge file")
    #bed_file = sys.argv[1]
    merge_file_contents = open(merge_file).read().rstrip("\n").split("\n")


source_trans_read_dict = {}  # source_trans_read_dict[source][trans_id][read_id] = 1
source_list = []

source_trans_list = []

for line in filelist_file_contents:
    # source_name	trans_read	file_type
    line_split = line.split("\t")

    source_name = line_split[0]
    transread_file = line_split[1]

    file_type = line_split[2]


    transread_file_contents = open(transread_file).read().rstrip("\n").split("\n")

    source_trans_read_dict[source_name] = {}
    source_list.append(source_name)

    if file_type == "trans_read":

        for transread_line in transread_file_contents:
            transread_line_split = transread_line.split("\t")

            id_line = transread_line_split[3]

            id_split = id_line.split(";")

            trans_id = id_split[0]
            read_id = id_split[1]

            gene_id = trans_id.split(".")[0]

            chrom = transread_line_split[0]
            t_start = transread_line_split[1]
            t_end = transread_line_split[2]
            strand = transread_line_split[5]
            num_exons = int(transread_line_split[9])

            if trans_id not in source_trans_read_dict[source_name]:
                source_trans_read_dict[source_name][trans_id] = {}
                source_trans_list.append(trans_id)

            source_trans_read_dict[source_name][trans_id][read_id] = 1

    elif file_type == "ref_anno":

        for transread_line in transread_file_contents:
            transread_line_split = transread_line.split("\t")

            id_line = transread_line_split[3]

            id_split = id_line.split(";")

            gene_id = id_split[0]
            trans_id = id_split[1]
            read_id = id_split[1]


            chrom = transread_line_split[0]
            t_start = transread_line_split[1]
            t_end = transread_line_split[2]
            strand = transread_line_split[5]
            num_exons = int(transread_line_split[9])

            if trans_id not in source_trans_read_dict[source_name]:
                source_trans_read_dict[source_name][trans_id] = {}
                source_trans_list.append(trans_id)

            source_trans_read_dict[source_name][trans_id][read_id] = 1

    elif file_type == "read_support":
        
        #"merge_gene_id","merge_trans_id","gene_read_count","trans_read_count","source_line","support_line"

        for transread_line in transread_file_contents:
            if transread_line.startswith("merge_gene_id"):
                continue
                
            transread_line_split = transread_line.split("\t")

            gene_id = transread_line_split[0]
            trans_id = transread_line_split[1]
            gene_read_count = transread_line_split[2]
            num_read_support = transread_line_split[3]
            this_source_line = transread_line_split[4]
            this_read_line = transread_line_split[5]
            

            if trans_id not in source_trans_read_dict[source_name]:
                source_trans_read_dict[source_name][trans_id] = {}
                source_trans_list.append(trans_id)

            this_read_line_split = this_read_line.split(";")

            for this_source_read_line in this_read_line_split:
                this_source_name = this_source_read_line.split(":")[0]


                #############
                # use this to deal with nanopore read ID's with the ":" character in them
                this_source_read_line_split = this_source_read_line.split(":")

                if len(this_source_read_line_split) > 2:

                    new_this_source_read_line_split = []
                    for i in xrange(len(this_source_read_line_split) - 1):
                        new_this_source_read_line_split.append(this_source_read_line_split[i+1])

                    new_this_source_read_line = ":".join(new_this_source_read_line_split)
                    read_list = new_this_source_read_line.split(",")

                else:

                    read_list = this_source_read_line.split(":")[1].split(",")
                #############


                for read_id in read_list:
                    source_trans_read_dict[source_name][trans_id][read_id] = 1
    
    elif file_type == "cluster":
        #cluster_id,read_id,read_type
        #transcript/0,m64012_181221_231243/12191118/ccs,FL
        #transcript/0,m64012_181221_231243/118686392/ccs,FL


        for transread_line in transread_file_contents:
            if transread_line.startswith("cluster_id"):
                continue
                
            transread_line_split = transread_line.split(",")

            trans_id = transread_line_split[0]
            read_id = transread_line_split[1]
            read_type = transread_line_split[2]
            
            if trans_id not in source_trans_read_dict[source_name]:
                source_trans_read_dict[source_name][trans_id] = {}
                source_trans_list.append(trans_id)

            source_trans_read_dict[source_name][trans_id][read_id] = 1


    else:
        print("Error with file type")
        print(file_type)
        sys.exit()




#HLphyDis2_00000001      3567978 3572092 G7.1;clontech_G7.1      40      -       3567978 3572092 255,0,0 1       4114    0
#HLphyDis2_00000001      3614048 3615498 G8.1;clontech_G8.1      40      -       3614048 3615498 255,0,0 1       1450    0
#HLphyDis2_00000001      3639112 3641626 G9.1;clontech_G9.1      40      -       3639112 3641626 255,0,0 1       2514    0
#HLphyDis2_00000001      3755241 3759319 G10.1;clontech_G10.1    40      -       3755241 3759319 255,0,0 1       4078    0
#HLphyDis2_00000001      3776405 3780736 G11.1;clontech_G11.1    40      -       3776405 3780736 255,0,0 1       4331    0

#print header
outline = "\t".join(["merge_gene_id","merge_trans_id","gene_read_count","trans_read_count","source_line","support_line"])
outfile.write(outline)
outfile.write("\n")

trans_gene_dict = {} # trans_gene_dict[trans_id] = gene id

gene_read_dict = {} # gene_read_dict[gene_id][source][read_id] = 1

# this is for when a merge file is given
if merge_file != "no_merge":
    merge_trans_dict = {} # merge_trans_dict[merge_trans_id][source name][source_trans_id] = read list
    merge_trans_list = []
    
    for line in merge_file_contents:

        if merge_type_flag == "tama" :
            line_split = line.split("\t")
    
            id_line = line_split[3]

            id_split = id_line.split(";")

            merge_trans_id = id_split[0]
            source_trans_id_line = id_split[1]


            if "trans_read.bed" in merge_file:
                source_name = source_list[0]
                source_trans_id = source_trans_id_line
            else:

                num_underscore_fields = len(source_trans_id_line.split("_"))

                if num_underscore_fields == 2:
                    source_name = source_trans_id_line.split("_")[0]
                    source_trans_id = source_trans_id_line.split("_")[1]
                else:
                    source_trans_id = source_trans_id_line.split("_")[num_underscore_fields - 1]
                    source_name_underscore_list = []
                    for i in xrange(num_underscore_fields-1):
                        source_name_underscore_list.append(source_trans_id_line.split("_")[i])

                    source_name = "_".join(source_name_underscore_list)

            merge_gene_id = merge_trans_id.split(".")[0]

            # collect gene id and gene read num info
            if merge_trans_id not in trans_gene_dict:
                trans_gene_dict[merge_trans_id] = merge_gene_id
            else:
                if merge_gene_id != trans_gene_dict[merge_trans_id]:
                    print("Error with trans_id having multiple gene_id.")
                    print("Terminating early.")
                    sys.exit()

            chrom = line_split[0]
            t_start = line_split[1]
            t_end = line_split[2]
            strand = line_split[5]

            if merge_trans_id not in merge_trans_dict:
                merge_trans_dict[merge_trans_id] = {}
                merge_trans_list.append(merge_trans_id)
            if source_name not in merge_trans_dict[merge_trans_id] :
                merge_trans_dict[merge_trans_id][source_name] = {}

            #print(source_name)
            #print(source_trans_id)
            read_list = list(source_trans_read_dict[source_name][source_trans_id].keys())

            merge_trans_dict[merge_trans_id][source_name][source_trans_id] = read_list

            if merge_gene_id not in gene_read_dict:
                gene_read_dict[merge_gene_id] = {}
            if source_name not in gene_read_dict[merge_gene_id]:
                gene_read_dict[merge_gene_id][source_name] = {}

            for read_id in read_list:
                gene_read_dict[merge_gene_id][source_name][read_id] = 1

        elif  merge_type_flag == "cupcake" :
            line_split = line.split("\t")
            #PB.2.1  transcript/111304,transcript/110485,transcript/89421,transcript/75502
            #PB.2.15 transcript/95975,transcript/96463
            merge_trans_id = line_split[0]
            source_trans_id_line = line_split[1]
            source_trans_id_list = source_trans_id_line.split(",")
            
            source_name = source_list[0] # should only be one source

            for source_trans_id in source_trans_id_list:
                merge_gene_id = merge_trans_id.split(".")[0] + "." + merge_trans_id.split(".")[1]

                # collect gene id and gene read num info
                if merge_trans_id not in trans_gene_dict:
                    trans_gene_dict[merge_trans_id] = merge_gene_id
                else:
                    if merge_gene_id != trans_gene_dict[merge_trans_id]:
                        print("Error with trans_id having multiple gene_id.")
                        print("Terminating early.")
                        sys.exit()

                if merge_trans_id not in merge_trans_dict:
                    merge_trans_dict[merge_trans_id] = {}
                    merge_trans_list.append(merge_trans_id)
                if source_name not in merge_trans_dict[merge_trans_id] :
                    merge_trans_dict[merge_trans_id][source_name] = {}

                read_list = list(source_trans_read_dict[source_name][source_trans_id].keys())

                merge_trans_dict[merge_trans_id][source_name][source_trans_id] = read_list

                if merge_gene_id not in gene_read_dict:
                    gene_read_dict[merge_gene_id] = {}
                if source_name not in gene_read_dict[merge_gene_id]:
                    gene_read_dict[merge_gene_id][source_name] = {}

                for read_id in read_list:
                    gene_read_dict[merge_gene_id][source_name][read_id] = 1
        

        elif  merge_type_flag == "filter" : #########################################
            line_split = line.split("\t")
            
            if line.startswith("old_gene_id"):
                continue
            
            #old_gene_id     old_trans_id    source_line     num_reads       new_gene_id     new_trans_id    num_exons
            #G1      G1.1    a,b     6       G1      G1.1    1
            #G1      G1.2    b       6       G1      removed_transcript      1

            merge_trans_id = line_split[5]
            source_trans_id = line_split[1]
            merge_gene_id = line_split[4]
            
            if merge_trans_id == "removed_transcript": # skip the transcripts that were removed in filtering
                continue

            source_name = source_list[0] # should only be one source
            

            # collect gene id and gene read num info
            if merge_trans_id not in trans_gene_dict:
                trans_gene_dict[merge_trans_id] = merge_gene_id
            else:
                if merge_gene_id != trans_gene_dict[merge_trans_id]:
                    print("Error with trans_id having multiple gene_id.")
                    print("Terminating early.")
                    print(merge_gene_id)
                    print(trans_gene_dict[merge_trans_id])
                    sys.exit()

            if merge_trans_id not in merge_trans_dict:
                merge_trans_dict[merge_trans_id] = {}
                merge_trans_list.append(merge_trans_id)
            if source_name not in merge_trans_dict[merge_trans_id] :
                merge_trans_dict[merge_trans_id][source_name] = {}

            read_list = list(source_trans_read_dict[source_name][source_trans_id].keys())

            merge_trans_dict[merge_trans_id][source_name][source_trans_id] = read_list

            if merge_gene_id not in gene_read_dict:
                gene_read_dict[merge_gene_id] = {}
            if source_name not in gene_read_dict[merge_gene_id]:
                gene_read_dict[merge_gene_id][source_name] = {}

            for read_id in read_list:
                gene_read_dict[merge_gene_id][source_name][read_id] = 1
    
    
    for merge_trans_id in merge_trans_list:
    
        this_read_dict = {} # this_read_dict[read id] = 1
        this_source_read_dict = {} # this_source_read_dict[source][read id] = 1
        this_read_source_dict = {} # this_read_source_dict[read id][source] = 1
        this_source_list = []
    
        for source_name in source_list:
    
            if source_name in merge_trans_dict[merge_trans_id]:
                if source_name not in this_source_read_dict:
                    this_source_read_dict[source_name] = {}
    
                this_source_list.append(source_name)
                for source_trans_id in merge_trans_dict[merge_trans_id][source_name]:
                    this_read_list = merge_trans_dict[merge_trans_id][source_name][source_trans_id]
    
                    for this_read_id in this_read_list:
                        this_read_dict[this_read_id] = 1
                        this_source_read_dict[source_name][this_read_id] = 1
                        
                        if this_read_id not in this_read_source_dict:
                            this_read_source_dict[this_read_id] = {}
                        
                        this_read_source_dict[this_read_id][source_name] = 1
                        
                        if len(list(this_read_source_dict[this_read_id].keys())) > 1:
                            
                            if duplicate_flag == "no_dup":
                                source_list = list(this_read_source_dict[this_read_id].keys())
                                print("Error with duplicate read ID's from different sources")
                                print(merge_trans_id)
                                print(source_list)
                                print(source_trans_id)
                                print(this_read_id)
                            
                                sys.exit()
                        
        num_read_support = len(list(this_read_dict.keys()))
    
        # get gene id and gene read support
        merge_gene_id = trans_gene_dict[merge_trans_id]
        gene_read_num = 0
        
        support_list = []
        for source_name in source_list:
            if source_name in this_source_read_dict:
                this_source_read_list = list(this_source_read_dict[source_name].keys())
                this_source_read_line = ",".join(this_source_read_list)
    
                this_source_read_line =  source_name + ":" + this_source_read_line
    
                support_list.append(this_source_read_line)
            
            if source_name in gene_read_dict[merge_gene_id]:
                # get gene id and gene read support
                this_gene_read_num = len(list(gene_read_dict[merge_gene_id][source_name].keys()))
                gene_read_num = gene_read_num + this_gene_read_num
    
        support_line = ";".join(support_list)
    
        this_read_line = ",".join(list(this_read_dict.keys()))
    
        this_source_line = ",".join(this_source_list)

        outline = "\t".join([merge_gene_id,merge_trans_id,str(gene_read_num),str(num_read_support),this_source_line,support_line])
        outfile.write(outline)
        outfile.write("\n")

else: # this is for no merge read support 
    num_sources = len(source_list)
    source_name = source_list[0]
    
    if num_sources > 1:
        print("Error: with no_merge option you can only supply one source.")
        print("Terminating early.")
        sys.exit()
    
    if source_name != "cluster":
        # collect gene id and read num
        for merge_trans_id in source_trans_list:
            merge_gene_id = merge_trans_id.split(".")[0]
            
            if merge_trans_id not in trans_gene_dict:
                trans_gene_dict[merge_trans_id] = merge_gene_id
            else:
                if merge_gene_id != trans_gene_dict[merge_trans_id]:
                    print("Error with trans_id having multiple gene_id.")
                    print("Terminating early.")
                    sys.exit()
        
            if merge_gene_id not in gene_read_dict:
                gene_read_dict[merge_gene_id] = {}
            if source_name not in gene_read_dict[merge_gene_id]:
                gene_read_dict[merge_gene_id][source_name] = {}
            
            read_list = list(source_trans_read_dict[source_name][merge_trans_id].keys())
            
            for read_id in read_list:
                gene_read_dict[merge_gene_id][source_name][read_id] = 1
            


    for merge_trans_id in source_trans_list:
            
        num_read_support = len(list(source_trans_read_dict[source_name][merge_trans_id].keys()))
    
        this_read_line = ",".join(list(source_trans_read_dict[source_name][merge_trans_id].keys()))
        support_line = source_name + ":" + this_read_line
    
        this_source_line = source_name

        # get gene id and gene read support
        if source_name != "cluster":
            merge_gene_id = trans_gene_dict[merge_trans_id]
            gene_read_num = len(list(gene_read_dict[merge_gene_id][source_name].keys()))
        else:
            merge_gene_id = "NA"
            gene_read_num = "NA"
    
        outline = "\t".join([merge_gene_id,merge_trans_id,str(gene_read_num),str(num_read_support),this_source_line,support_line])
        outfile.write(outline)
        outfile.write("\n")










