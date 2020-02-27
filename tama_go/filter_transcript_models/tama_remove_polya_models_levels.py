import re
import sys
import time


import os
import argparse


# This script uses the TAMA collapse and TAMA merge outputs to remove Poly-A models
# Poly-A models are genomic contaminaton or internal truncation 

ap = argparse.ArgumentParser(description='This script uses the TAMA collapse and TAMA merge outputs to remove Poly-A models')


ap.add_argument('-b', type=str, nargs=1, help='Annotation bed file')
ap.add_argument('-f', type=str, nargs=1, help='Filelist file with Poly-A file names')
ap.add_argument('-r', type=str, nargs=1, help='Read support file')
ap.add_argument('-o', type=str, nargs=1, help='Output prefix (required)')

ap.add_argument('-p', type=str, nargs=1, help='Percent poly-A threshold (default of 75.0)')

ap.add_argument('-l', type=str, nargs=1, help='Level of removal (gene or transcript level)')
ap.add_argument('-a', type=str, nargs=1, help='Remove all models with Poly-A (all_polya or singleton_polya). Default is singleton_polya.')
ap.add_argument('-k', type=str, nargs=1, help='Keep all multi-exon models (keep_multi or remove_multi)')

opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0


if not opts.b:
    print("Annotation bed file missing")
    missing_arg_flag = 1
if not opts.r:
    print("Read Support file missing")
    missing_arg_flag = 1
if not opts.f:
    print("Filelist file missing")
    missing_arg_flag = 1
if not opts.o:
    print("output prefix missing")
    missing_arg_flag = 1

if missing_arg_flag == 1:
    print("Please try again with complete arguments")

filter_flag = "NA"

if not opts.l:
    print("Default to gene level filtering.")
    filter_flag = "gene"
else:
    filter_flag = opts.l[0]

multi_flag = "NA"

if not opts.k:
    print("Default to remove multi-exon for transcript level.")
    multi_flag = "remove_multi"
else:
    multi_flag = opts.k[0]


support_flag = "NA"

if not opts.a:
    print("Default is to remove only singleton Poly-A models.")
    support_flag = "singleton_polya"
else:
    support_flag = opts.a[0]
    
polya_threshold = 75.0

if not opts.p:
    print("Default to 75.0 Poly-A threshold.")
else:
    polya_threshold = float(opts.p[0])

bed_file = opts.b[0]
filelist_file = opts.f[0]
read_file = opts.r[0]

outfile_prefix = opts.o[0]

print("opening bed file")
#bed_file = sys.argv[1]
bed_file_contents = open(bed_file).read().rstrip("\n").split("\n")

print("opening filelist file")
#bed_file = sys.argv[1]
filelist_file_contents = open(filelist_file).read().rstrip("\n").split("\n")

print("opening read support file")
#bed_file = sys.argv[1]
read_file_contents = open(read_file).read().rstrip("\n").split("\n")


outbed_name = outfile_prefix + ".bed"
report_name = outfile_prefix + "_polya_report.txt"
polya_name = outfile_prefix + "_trash_polya.bed"

support_name = outfile_prefix + "_polya_support.txt"

outfile_bed = open(outbed_name,"w")
outfile_report = open(report_name,"w")
outfile_polya = open(polya_name,"w")

outfile_support = open(support_name,"w")


gene_trans_dict = {} # gene_trans_dict[gene_id][trans_id] = trans_split

gene_trans_list_dict = {} # gene_trans_list_dict[gene_id] = list of trans

trans_exon_dict = {} # trans_exon_dict[trans_id] =  num exons

gene_list = []

for line in bed_file_contents:
    line_split = line.split("\t")

    id_line = line_split[3]

    id_split = id_line.split(";")

    gene_id = id_split[0]
    trans_id = id_split[1]

    if gene_id not in gene_trans_dict:
        gene_trans_dict[gene_id] = {}
        gene_trans_list_dict[gene_id] = []
        gene_list.append(gene_id)
    
    gene_trans_dict[gene_id][trans_id] = line_split
    gene_trans_list_dict[gene_id].append(trans_id)

    num_exons = int(line_split[9])

    trans_exon_dict[trans_id] = num_exons



#HLphyDis2_00000001      3567978 3572092 G7.1;clontech_G7.1      40      -       3567978 3572092 255,0,0 1       4114    0
#HLphyDis2_00000001      3614048 3615498 G8.1;clontech_G8.1      40      -       3614048 3615498 255,0,0 1       1450    0
#HLphyDis2_00000001      3639112 3641626 G9.1;clontech_G9.1      40      -       3639112 3641626 255,0,0 1       2514    0
#HLphyDis2_00000001      3755241 3759319 G10.1;clontech_G10.1    40      -       3755241 3759319 255,0,0 1       4078    0
#HLphyDis2_00000001      3776405 3780736 G11.1;clontech_G11.1    40      -       3776405 3780736 255,0,0 1       4331    0


merge_trans_read_dict = {} # merge_trans_read_dict[merge_trans_id][read id] = 1
merge_source_read_dict = {} # merge_source_read_dict[merge trans id][source][read_id] = 1

for line in read_file_contents:
    # merge_trans_id  support_line    num_read_support        source_line

    #merge_gene_id   merge_trans_id  gene_read_count trans_read_count        source_line     support_line
    #G1      G1.1    1039    8       a,b     a:m64012_181221_231243/168495354/ccs,m64012_181221_231243/6686664/ccs

    if line.startswith("merge_gene_id"):
        continue

    line_split = line.split("\t")

    merge_gene_id = line_split[0]
    merge_trans_id = line_split[1]
    this_read_line = line_split[5]
    # merge_gene_id = merge_trans_id.split(".")[0]

    if merge_trans_id not in merge_trans_read_dict:
        merge_trans_read_dict[merge_trans_id] = {}
        merge_source_read_dict[merge_trans_id] = {}

    this_read_line_split = this_read_line.split(";")

    for this_source_read_line in this_read_line_split:
        this_source_name = this_source_read_line.split(":")[0]

        merge_source_read_dict[merge_trans_id][this_source_name] = {}

        this_source_read_line_split = this_source_read_line.split(":")

        if len(this_source_read_line_split) > 2:

            new_this_source_read_line_split = []
            for i in xrange(len(this_source_read_line_split) - 1):
                new_this_source_read_line_split.append(this_source_read_line_split[i+1])

            new_this_source_read_line = ":".join(new_this_source_read_line_split)
            read_list = new_this_source_read_line.split(",")

        else:

            read_list = this_source_read_line.split(":")[1].split(",")

        for read_id in read_list:
            merge_trans_read_dict[merge_trans_id][read_id] = 1
            merge_source_read_dict[merge_trans_id][this_source_name][read_id] = 1


source_polya_read_dict = {} # source_polya_read_dict[source][read_id] = polya line split
#  source_polya_trans_dict = {} # source_polya_trans_dict[source][trans_id] = polya line

for line in filelist_file_contents:
    # source_name	polya
    line_split = line.split("\t")
    
    source_name = line_split[0]
    polya_file = line_split[1]
    

    polya_file_contents = open(polya_file).read().rstrip("\n").split("\n")
    
    #cluster_id      trans_id        strand  a_percent       a_count sequence
    #m54041_170531_132524/67174556/30_403_CCS        G5.1    +       90.0    18      TAAAAAAAAAAAAGAAAAAA
    #m54041_170609_082817/20119745/30_408_CCS        G5.1    +       90.0    18      TAAAAAAAAAAAAGAAAAAA
    

    if source_name not in source_polya_read_dict:
        source_polya_read_dict[source_name] = {}
    
    for polya_line in polya_file_contents:
        if polya_line.startswith("cluster_id"):
            continue

        polya_line_split = polya_line.split("\t")

        read_id = polya_line_split[0]
        trans_id = polya_line_split[1]
        
        source_polya_read_dict[source_name][read_id] = polya_line_split



def write_report(gene_id, trans_id,source_line,total_trans_num_reads, new_gene_id, new_trans_id,num_exons):

    outline = "\t".join([gene_id,trans_id,source_line,str(total_trans_num_reads),new_gene_id,new_trans_id,str(num_exons)])
    outfile_report.write(outline)
    outfile_report.write("\n")

def write_polya(trans_split):
    # write out removed models
    trans_line = "\t".join(trans_split)
    outfile_polya.write(trans_line)
    outfile_polya.write("\n")


def write_support(this_merge_trans_id,this_support_source,polya_support_list):
    outlist = []
    outlist.append(this_merge_trans_id)
    outlist.append(this_support_source)
    outlist.extend(polya_support_list)

    outline = "\t".join(outlist)

    outfile_support.write(outline)
    outfile_support.write("\n")


def write_bed(trans_split):
    trans_line = "\t".join(trans_split)
    outfile_bed.write(trans_line)
    outfile_bed.write("\n")


removed_trans_count = 0

removed_gene_count = 0

new_gene_num = 0

new_trans_dict = {} # new_trans_dict[new trans] = old trans
new_gene_dict = {}  # new_gene_dict[new gene] = old gene

new_gene_list = []

# write out header
outline = "\t".join(["old_gene_id","old_trans_id","source_line","num_reads","new_gene_id","new_trans_id","num_exons"])
outfile_report.write(outline)
outfile_report.write("\n")

# write out header
#m54041_170609_082817/20119745/30_408_CCS        G5.1    +       90.0    18      TAAAAAAAAAAAAGAAAAAA
outline = "\t".join(["trans_id","source","read_id","source_trans_id","strand","percent_polya","a_count","polya_seq"])
outfile_support.write(outline)
outfile_support.write("\n")

for gene_id in gene_list:
    output_flag = "NA"
    
    trans_list = gene_trans_list_dict[gene_id]


    #########################################################
    #######################################################

    this_removed_trans_count = 0
    this_total_trans_count = 0

    new_trans_num = 0

    ############################
    gene_max_exons = 0

    single_exon_trans_dict = {} # single_exon_trans_dict[trans_id] = num reads
    multi_exon_trans_dict = {} # multi_exon_trans_dict[trans_id] = num reads

    this_trans_num_reads_dict = {} # this_trans_num_reads_dict[trans_id] = num reads

    single_read_trans_dict = {} # single_read_trans_dict[trans_id] = num exons
    multi_read_trans_dict = {}  # multi_read_trans_dict[trans_id] = num exons

    this_polya_trans_dict = {} # this_polya_trans_dict[trans_id] = 1

    this_trans_read_dict = {} # this_trans_read_dict[trans_id][read_id] = 1


    this_read_dict = {} # this_read_dict[read_id] = 1
    this_polya_read_dict = {} # this_polya_read_dict[read_id] = 1

    ##################################
    this_trans_read_polya_support_dict = {}  # this_trans_read_polya_support_dict[trans_id][read_id] = polya support list
    this_read_source_dict = {}  # this_read_source_dict[read_id] = source


    for trans_id in trans_list:
        num_exons = trans_exon_dict[trans_id]

        if num_exons > gene_max_exons:
            gene_max_exons = num_exons

        this_total_trans_count += 1
        trans_num_sources = len(list(merge_source_read_dict[trans_id].keys()))
        this_source_list = list(merge_source_read_dict[trans_id].keys())


        total_trans_num_reads = len(list(merge_trans_read_dict[trans_id].keys()))

        polya_trans_num_reads = 0

        # print("Trans: " + trans_id) #######################

        for this_source_name in this_source_list:

            # print(this_source_name) #######################

            this_read_id_list = list(merge_source_read_dict[trans_id][this_source_name].keys())

            for this_read_id in this_read_id_list:

                # print(this_read_id) #######################

                if trans_id not in this_trans_read_dict:
                    this_trans_read_dict[trans_id] = {}
                if this_read_id not in this_trans_read_dict[trans_id]:
                    this_trans_read_dict[trans_id][this_read_id] = 1

                ##################################
                if this_read_id not in this_read_source_dict:
                    this_read_source_dict[this_read_id] = this_source_name
                elif this_read_source_dict[this_read_id] != this_source_name:
                    print("Error with duplicate read ID in different sources")
                    sys.exit()


                this_read_dict[this_read_id] = 1

                if this_read_id in source_polya_read_dict[this_source_name]:

                    polya_percent = float(source_polya_read_dict[this_source_name][this_read_id][3])

                    # print(polya_percent) ##################################

                    if polya_percent >= polya_threshold:
                        polya_trans_num_reads += 1

                        this_polya_read_dict[this_read_id] = 1

                        ##################################
                        if trans_id not in this_trans_read_polya_support_dict:
                            this_trans_read_polya_support_dict[trans_id] = {}

                        this_trans_read_polya_support_dict[trans_id][this_read_id] = source_polya_read_dict[this_source_name][this_read_id]

        #if trans_id == "G2.1":
        #    print("Pre checks")
        #    print("Trans: " + trans_id)
        #    print(polya_trans_num_reads)
        #    print(total_trans_num_reads)
        #    sys.exit()

        if polya_trans_num_reads == total_trans_num_reads:
            this_polya_trans_dict[trans_id] = 1

        if num_exons == 1:
            single_exon_trans_dict[trans_id] = total_trans_num_reads
        elif num_exons > 1:
            multi_exon_trans_dict[trans_id] = total_trans_num_reads

        this_trans_num_reads_dict[trans_id] = total_trans_num_reads

        if total_trans_num_reads == 1:
            single_read_trans_dict[trans_id] = num_exons
        elif total_trans_num_reads > 1:
            multi_read_trans_dict[trans_id] = num_exons


    this_remove_trans_dict = {} # this_remove_trans_dict[trans_id] = 1
    this_keep_trans_dict = {} # this_keep_trans_dict[trans_id] = 1

    for this_trans_id in trans_list:
        trans_split = gene_trans_dict[gene_id][this_trans_id]
        this_num_exons = trans_exon_dict[this_trans_id]
        source_line = ",".join(list(merge_source_read_dict[trans_id].keys()))

        if this_trans_id in this_polya_trans_dict:

            this_trans_num_reads = this_trans_num_reads_dict[this_trans_id]


            # resolve gene level filtration
            if filter_flag == "gene":
                if len(trans_list) > 1: # keep if gene has multiple transcripts
                    this_keep_trans_dict[this_trans_id] = 1
                    continue
                elif len(trans_list) == 1:
                    if multi_flag == "keep_multi" and this_num_exons > 1:# keep if multiple exons
                        this_keep_trans_dict[this_trans_id] = 1
                        continue

                    if support_flag == "singleton_polya" and this_trans_num_reads > 1:  # keep trans if multiple reads
                        this_keep_trans_dict[this_trans_id] = 1
                        continue

                    if support_flag == "all_polya":
                        this_polya_read_count = 0
                        this_all_read_count = 0

                        for this_read_id in this_trans_read_dict[this_trans_id]:
                            this_all_read_count += 1
                            if this_read_id in this_polya_read_dict:
                                this_polya_read_count += 1

                        if this_polya_read_count < this_all_read_count:
                            this_keep_trans_dict[this_trans_id] = 1
                            continue

                # what ever makes it through the above conditions will be removed
                this_remove_trans_dict[this_trans_id] = 1

            # resolve transcript level filtration
            if filter_flag == "transcript":

                if multi_flag == "keep_multi" and this_num_exons > 1:  # multi flag is for keeping multiple exon trans
                    this_keep_trans_dict[this_trans_id] = 1
                    continue

                if support_flag == "singleton_polya" and this_trans_num_reads > 1: # keep trans if multiple reads
                    this_keep_trans_dict[this_trans_id] = 1
                    continue

                if support_flag == "all_polya":
                    this_polya_read_count = 0
                    this_all_read_count = 0

                    for this_read_id in this_trans_read_dict[this_trans_id]:
                        this_all_read_count += 1
                        if this_read_id in this_polya_read_dict:
                            this_polya_read_count += 1

                    if this_polya_read_count < this_all_read_count:

                        this_keep_trans_dict[this_trans_id] = 1
                        continue



                # what ever makes it through the above conditions will be removed
                this_remove_trans_dict[this_trans_id] = 1

        else:

            this_keep_trans_dict[this_trans_id] = 1

    # make new gene id if there are keepers
    if len(list(this_keep_trans_dict.keys())) > 0:
        new_gene_num += 1
        new_gene_id = "G" + str(new_gene_num)
        new_gene_list.append(new_gene_id)

    new_trans_num = 0

    for this_trans_id in trans_list:
        trans_split = gene_trans_dict[gene_id][this_trans_id]
        this_num_exons = trans_exon_dict[this_trans_id]
        source_line = ",".join(list(merge_source_read_dict[trans_id].keys()))

        if this_trans_id in this_remove_trans_dict:
            write_polya(trans_split)

            if filter_flag == "gene":
                write_report(gene_id, this_trans_id,source_line,total_trans_num_reads, "removed_gene","removed_transcript",this_num_exons)

            elif filter_flag == "transcript":
                if len(this_keep_trans_dict.keys()) == 0:
                    write_report(gene_id, this_trans_id,source_line,total_trans_num_reads, "removed_gene","removed_transcript",this_num_exons)

                else:
                    write_report(gene_id, this_trans_id,source_line,total_trans_num_reads, new_gene_id,"removed_transcript",this_num_exons)


            if this_trans_id in this_trans_read_polya_support_dict:
                for this_read_id in this_trans_read_polya_support_dict[this_trans_id]:
                    polya_support_list = this_trans_read_polya_support_dict[this_trans_id][this_read_id]

                    this_support_source = this_read_source_dict[this_read_id]

                    write_support(this_trans_id, this_support_source, polya_support_list)

        elif this_trans_id in this_keep_trans_dict:

            new_trans_num += 1

            new_trans_id = new_gene_id + "." + str(new_trans_num)

            trans_split[3] = new_gene_id + ";" + new_trans_id

            write_bed(trans_split)

            write_report(gene_id, this_trans_id,source_line,total_trans_num_reads, new_gene_id,new_trans_id,this_num_exons)






    













