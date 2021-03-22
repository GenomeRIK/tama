import re
import sys
import time


import os
import argparse


# This script uses the TAMA collapse and TAMA merge outputs to remove single read models


ap = argparse.ArgumentParser(description='This script uses the TAMA collapse and TAMA merge outputs to remove single read models')


ap.add_argument('-b', type=str, nargs=1, help='Annotation bed file')
ap.add_argument('-r', type=str, nargs=1, help='Read support file')
ap.add_argument('-o', type=str, nargs=1, help='Output prefix (required)')

ap.add_argument('-l', type=str, nargs=1, help='Level of removal (gene or transcript level). Gene level will only remove genes with a single read, transcript level will remove all singleton transcripts.')

ap.add_argument('-k', type=str, nargs=1, help='Default to keep all multi-exon models (keep_multi or remove_multi)')

ap.add_argument('-s', type=str, nargs=1, help='Requires models to have support from at least this number of sources. Default is 1')

ap.add_argument('-n', type=str, nargs=1, help='Requires models to have support from at least this number of reads. Default is 2')

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
    print("Default to keep multi-exon models.")
    multi_flag = "keep_multi"
else:
    multi_flag = opts.k[0]


source_support_flag = "NA"

if not opts.s:
    print("Default only requires single source support")
    source_support_flag = 1
else:
    source_support_flag = int(opts.s[0])

if not opts.n:
    print("Default only requires read support of 2")
    read_support_threshold = 2
else:
    read_support_threshold = int(opts.n[0])

bed_file = opts.b[0]

read_file = opts.r[0]

outfile_prefix = opts.o[0]

print("opening bed file")
#bed_file = sys.argv[1]
bed_file_contents = open(bed_file).read().rstrip("\n").split("\n")
##############################################################################continue here! 2019/06/22

print("opening read support file")
#bed_file = sys.argv[1]
read_file_contents = open(read_file).read().rstrip("\n").split("\n")


outbed_name = outfile_prefix + ".bed"
report_name = outfile_prefix + "_singleton_report.txt"
single_name = outfile_prefix + "_singleton.bed"

outfile_bed = open(outbed_name,"w")
outfile_report = open(report_name,"w")

outfile_single = open(single_name,"w")



#HLphyDis2_00000001      3567978 3572092 G7.1;clontech_G7.1      40      -       3567978 3572092 255,0,0 1       4114    0
#HLphyDis2_00000001      3614048 3615498 G8.1;clontech_G8.1      40      -       3614048 3615498 255,0,0 1       1450    0
#HLphyDis2_00000001      3639112 3641626 G9.1;clontech_G9.1      40      -       3639112 3641626 255,0,0 1       2514    0
#HLphyDis2_00000001      3755241 3759319 G10.1;clontech_G10.1    40      -       3755241 3759319 255,0,0 1       4078    0
#HLphyDis2_00000001      3776405 3780736 G11.1;clontech_G11.1    40      -       3776405 3780736 255,0,0 1       4331    0


merge_trans_read_dict = {} # merge_trans_read_dict[merge_trans_id][read id] = 1
merge_source_read_dict = {} # merge_source_read_dict[mnerge trans id][source][read_id] = 1

for line in read_file_contents:
    # merge_trans_id  support_line    num_read_support        source_line

    # merge_gene_id   merge_trans_id  gene_read_count trans_read_count        source_line     support_line
    # G1      G1.1    1039    8       a,b     a:m64012_181221_231243/168495354/ccs,m64012_181221_231243/6686664/ccs

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

        this_source_read_line_split = this_source_name.split(":")

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

    num_exons = int(line_split[9])

    total_trans_num_reads = 0

    for source_name in merge_source_read_dict[trans_id]:
        this_source_num_reads = len(list(merge_source_read_dict[trans_id][source_name].keys()))
        total_trans_num_reads = total_trans_num_reads + this_source_num_reads

    if gene_id not in gene_trans_dict:
        gene_trans_dict[gene_id] = {}
        gene_trans_list_dict[gene_id] = []
        gene_list.append(gene_id)
    
    gene_trans_dict[gene_id][trans_id] = line_split
    gene_trans_list_dict[gene_id].append(trans_id)



    trans_exon_dict[trans_id] = num_exons




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

for gene_id in gene_list:
    output_flag = "NA"
    
    trans_list = gene_trans_list_dict[gene_id]

    #print(trans_list)
    
    #remove on gene level for single read gene support
    if len(trans_list) == 1:
        trans_id = trans_list[0]

        num_exons = trans_exon_dict[trans_id]
        num_sources = len(list(merge_source_read_dict[trans_id].keys()))
        if num_sources == 1:
            source_name = list(merge_source_read_dict[trans_id].keys())[0]

            num_reads = len(list(merge_source_read_dict[trans_id][source_name].keys()))

            # num_reads = len(list(source_trans_read_dict[source_name][sub_trans_id].keys()))
            
            if num_reads == 1:
                if multi_flag == "keep_multi" and num_exons > 1:  # multi flag is for keeping multiple exon trans
                    output_flag = "keep"

                else: # remove gene if it is single exon or if multi flag is set to remove

                    output_flag = "remove"
                    removed_trans_count += 1
                    removed_gene_count += 1

                    source_line = ",".join(list(merge_source_read_dict[trans_id].keys()))

                    outline = "\t".join([gene_id,trans_id,source_line,str(num_reads),"removed_gene","removed_transcript",str(num_exons)])
                    outfile_report.write(outline)
                    outfile_report.write("\n")

                    # write out removed models
                    trans_split = gene_trans_dict[gene_id][trans_id]
                    trans_line = "\t".join(trans_split)
                    outfile_single.write(trans_line)
                    outfile_single.write("\n")

                    continue

    #new_gene_num += 1
    #new_gene_id = "G" + str(new_gene_num)
    new_trans_num = 0
    # remove on transcript level

    keep_trans_dict = {} # keep_trans_dict[trans_id] = 1
    remove_trans_dict = {} # remove_trans_dict[trans_id] = 1
    if filter_flag == "transcript":

        new_trans_list = []
        new_report_list = []

        this_removed_trans_count = 0
        this_total_trans_count = 0
    
        for trans_id in trans_list:
            num_exons = trans_exon_dict[trans_id]
            this_total_trans_count += 1
            num_sources = len(list(merge_source_read_dict[trans_id].keys()))
            this_source_list = list(merge_source_read_dict[trans_id].keys())


            total_trans_num_reads = len(list(merge_trans_read_dict[trans_id].keys()))

            #print(total_trans_num_reads)
            #print(trans_id)
            #print(list(merge_trans_read_dict[trans_id].keys()))
            #sys.exit()

            if source_support_flag > 1:
                if num_sources >= source_support_flag:
                    output_flag = "keep"

                else:
                    # note that keep_multi flag trunmps num sources flag!!! ###################################
                    if multi_flag == "keep_multi" and num_exons > 1:  # multi flag is for keeping multiple exon trans
                        output_flag = "keep"
                    else:  # remove if single exon or if multi flag is set to remove multi exon transcripts
                        output_flag = "remove"

            elif source_support_flag == 1:

                if num_sources > 1 and total_trans_num_reads >= read_support_threshold:
                    output_flag = "keep"

                else:

                    if total_trans_num_reads >= read_support_threshold:
                        output_flag = "keep"
                    else:

                        if multi_flag == "keep_multi" and num_exons > 1:  # multi flag is for keeping multiple exon trans
                            output_flag = "keep"
                        else: # remove if single exon or if multi flag is set to remove multi exon transcripts
                            output_flag = "remove"


            if output_flag ==  "keep":
                keep_trans_dict[trans_id] = 1
                print("keep "+ trans_id + " " + str(total_trans_num_reads) +" "+ str(read_support_threshold) + " numexons "+ str(num_exons))
                #sys.exit()

            elif output_flag ==  "remove":
                remove_trans_dict[trans_id] =  1


            else:
                print("Error with output flag")
                print(output_flag)
                print(gene_id)
                sys.exit()

        if len(list(keep_trans_dict.keys())) > 0:
            new_gene_num += 1
            new_gene_id = "G" + str(new_gene_num)
            new_gene_list.append(new_gene_id)
        else:
            new_gene_id = "removed_trans_level"

        for trans_id in trans_list:
            if trans_id in keep_trans_dict:
                trans_split = gene_trans_dict[gene_id][trans_id]

                new_trans_num += 1

                new_trans_id = new_gene_id + "." + str(new_trans_num)

                trans_split[3] = new_gene_id + ";" + new_trans_id

                trans_line = "\t".join(trans_split)

                new_trans_list.append(trans_line)
                #outfile_bed.write(trans_line)
                #outfile_bed.write("\n")

                source_line = ",".join(list(merge_source_read_dict[trans_id].keys()))

                # source_trans_line = ",".join(list(merge_trans_dict[trans_id].keys()))

                total_trans_num_reads = len(list(merge_trans_read_dict[trans_id].keys()))
                num_exons = trans_exon_dict[trans_id]

                outline = "\t".join([gene_id, trans_id,source_line, str(total_trans_num_reads),new_gene_id, new_trans_id, str(num_exons)])
                #print(trans_id)
                #print(total_trans_num_reads)

                new_report_list.append(outline)
                #outfile_report.write(outline)
                #outfile_report.write("\n")
            elif trans_id in remove_trans_dict:

                removed_trans_count += 1
                this_removed_trans_count += 1

                source_line = ",".join(list(merge_source_read_dict[trans_id].keys()))

                # source_trans_line = ",".join(list(merge_trans_dict[trans_id].keys()))

                total_trans_num_reads = len(list(merge_trans_read_dict[trans_id].keys()))
                num_exons = trans_exon_dict[trans_id]

                outline = "\t".join([gene_id, trans_id,source_line, str(total_trans_num_reads),new_gene_id, "removed_transcript",str(num_exons)])
                #outfile_report.write(outline)
                #outfile_report.write("\n")
                new_report_list.append(outline)

                # write out removed models
                trans_split = gene_trans_dict[gene_id][trans_id]
                trans_line = "\t".join(trans_split)
                outfile_single.write(trans_line)
                outfile_single.write("\n")


        # deal with writing out
        if this_removed_trans_count == this_total_trans_count:
            for report_line in new_report_list:
                report_split = report_line.split("\t")
                report_split[4] = "all_trans_removed"

                new_report_line = "\t".join(report_split)
                outfile_report.write(new_report_line)
                outfile_report.write("\n")
        else:
            for report_line in new_report_list:

                outfile_report.write(report_line)
                outfile_report.write("\n")

        for new_trans_line in new_trans_list:
            outfile_bed.write(new_trans_line)
            outfile_bed.write("\n")
            
    elif filter_flag == "gene":

        if source_support_flag > 1:
            gene_num_sources = 0
            gene_num_source_dict = {}  # gene_num_source_dict[source name] = 1

            for trans_id in trans_list:
                trans_split = gene_trans_dict[gene_id][trans_id]
                this_source_list = list(merge_source_read_dict[trans_id].keys())

                for this_source in this_source_list:
                    gene_num_source_dict[this_source] = 1

            gene_num_sources = len(list(gene_num_source_dict.keys()))

            if gene_num_sources < source_support_flag: # skip this gene if it does not meet number of sources requirement
                removed_gene_count += 1
                for trans_id in trans_list:
                    removed_trans_count += 1

                    source_line = ",".join(list(merge_source_read_dict[trans_id].keys()))

                    total_trans_num_reads = len(list(merge_trans_read_dict[trans_id].keys()))


                    outline = "\t".join([gene_id,trans_id,source_line,str(total_trans_num_reads),"removed_gene","removed_transcript",str(num_exons)])
                    outfile_report.write(outline)
                    outfile_report.write("\n")

                    # write out removed models
                    trans_split = gene_trans_dict[gene_id][trans_id]
                    trans_line = "\t".join(trans_split)
                    outfile_single.write(trans_line)
                    outfile_single.write("\n")
                continue

        for trans_id in trans_list:
            trans_split = gene_trans_dict[gene_id][trans_id]
            num_exons = trans_exon_dict[trans_id]

            this_source_list = list(merge_source_read_dict[trans_id].keys())

            total_trans_num_reads = len(list(merge_trans_read_dict[trans_id].keys()))

            #print(total_trans_num_reads)
            #print(list(merge_trans_read_dict[trans_id].keys()))


            new_trans_num += 1

            if new_trans_num == 1:
                new_gene_num += 1
                new_gene_id = "G" + str(new_gene_num)

                new_gene_list.append(new_gene_id)

            new_trans_id = new_gene_id + "." + str(new_trans_num)

            trans_split[3] = new_gene_id + ";" + new_trans_id
            
            trans_line = "\t".join(trans_split)
            outfile_bed.write(trans_line)
            outfile_bed.write("\n")

            source_line = ",".join(list(merge_source_read_dict[trans_id].keys()))

            # source_trans_line = ",".join(list(merge_trans_dict[trans_id].keys()))

            outline = "\t".join([gene_id, trans_id,source_line,str(total_trans_num_reads), new_gene_id, new_trans_id,str(num_exons)])
            outfile_report.write(outline)
            outfile_report.write("\n")
    else:
        print("Error with filter flag: " + filter_flag)
        print("Exiting prematurely.")
        sys.exit()
        
    













