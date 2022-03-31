import re
import sys
import time


import os
import argparse


# This script uses the TAMA read levels support file for sampling curves


ap = argparse.ArgumentParser(description='This script uses the TAMA read support levels file to create a saturation curve')


ap.add_argument('-r', type=str, nargs=1, help='Read support file')

ap.add_argument('-b', type=str, nargs=1, help='Read bin size')

ap.add_argument('-o', type=str, nargs=1, help='Output file name')


opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0


if not opts.r:
    print("Read support file missing")
    missing_arg_flag = 1

if not opts.b:
    print("Read bin size missing")
    missing_arg_flag = 1

if not opts.o:
    print("Output file name missing")
    missing_arg_flag = 1


if missing_arg_flag == 1:
    print("Please try again with complete arguments")




report_file = opts.r[0]

read_bin = int(opts.b[0])

outfile_name = opts.o[0]

print("opening report file")

report_file_contents = open(report_file).read().rstrip("\n").split("\n")
##############################################################################continue here! 2019/06/22


outfile = open(outfile_name,"w")


read_gene_dict = {} # read_gene_dict[read_id] = gene_id

all_gene_dict = {} # all_gene_dict[gene_id] = 1



for line in report_file_contents:


#merge_gene_id   merge_trans_id  gene_read_count trans_read_count        source_line     support_line
#G1      G1.1    50      27      1,2,3,5,7,8,9,10,12,13,14,15,16,20      1:85:1177|9fd7a7bb-6855-4f7f-ab6b-6dbc05749667,
# 84:1370|1766e447-0938-4b76-9597-cd91cb920157,88:1207|e1815fab-159e-4cf6-9c82-b5138ff5707c,81:1130|63452154-836a-4229-80ae-de246dc52309;2:
# 86:1168|a1a63a71-47e6-44d4-b68f-dca0928ae38b,78:1384|9ad1788f-8f7e-4691-b425-ea42368b0e6a;3:79:1168|786f5f4c-2fb7-4905-a259-3316b303eccc,
# 83:1248|43fa11c3-a11a-46d6-843c-18ec6c191fb9,85:1328|2177475f-a9bf-4c9b-b891-3786e7be9aae;5:75:1187|2400629e-b440-446d-a058-0d32b4530b90;7:
# 82:1394|53a9db5a-2645-457c-b144-79323261bca3,80:1155|483860c7-d11b-4b65-bc36-9e313ca15f16;8:86:1098|874034f7-c5d2-473e-a215-4bdf20aa189c,
# 80:1255|5b287434-9cb1-4cef-91f9-00ccb23c78b5;9:88:1380|94e73afb-2ecf-4c70-87ac-4b9aed5fa436,90:1270|aaf8848b-0b9b-4fd6-85e2-18f6dc94ee6b;
# 10:78:1155|84d1e045-ab4c-4fb4-af37-8f71cdfef937;12:77:1167|6ddf4f6b-dd54-46df-88f8-d83982e092d8;13:89:1422|49ef2889-1b0b-4877-b74b-6dc480a955c8,
# 74:1238|4ba7116f-c238-4a2b-9a73-fc5508f32133,83:1239|296f397a-0783-46b8-a22b-8f6ffc8f3af1;14:88:1257|81cc9e2d-3333-494e-b1b6-9b49cbfc1e6a,
# 85:1166|5566136d-bb13-4419-b639-705aa1319349;15:89:1387|3198bd0b-2a7e-4002-bcaa-1ab3f95d731f;16:71:1134|8088bd4f-4748-4da0-8fa8-6d45c5251e25,
# 73:1229|870bd315-8ea4-4c87-86d5-83533e32c4cd;20:79:1137|1543aa02-ab4d-4bca-92c6-d434ed094461
#G1      G1.2    50      2       12,13   12:77:1598|e1bd00e1-90ca-4aff-b23c-c118bc8cd28b;13:88:1654|1c3c1984-88ca-4d53-a660-ee03952bede2
#G1      G1.3    50      1       16      16:91:1031|5778176c-374e-483b-a18b-8881be593047


    if line.startswith("gene_id"):
        continue

    line_split = line.split("\t")

    merge_gene_id = line_split[0]
    num_reads = line_split[2]
    support_line = line_split[5]

    support_list = support_line.split(";")

    if merge_gene_id not in all_gene_dict:
        all_gene_dict[merge_gene_id] = 1

    for source_line in support_list:
        source_split = source_line.split(":")

        source_name = source_split[0]

        source_split.pop(0)

        read_line = ":".join(source_split)

        read_list = read_line.split(",")

        for read_id in read_list:

            if read_id not in read_gene_dict:
                read_gene_dict[read_id] = merge_gene_id
            elif merge_gene_id != read_gene_dict[read_id] :
                print("Error with read support multiple genes!")
                print(read_id +"\t" + read_gene_dict[read_id]+"\t" +merge_gene_id)
                sys.exit()

all_read_list = list(read_gene_dict.keys())

read_count = 0
gene_count = 0

check_gene_dict = {} # check_gene_dict[gene_id] = 1

outline = "read_count" + "\t" + "gene_count"

outfile.write(outline)
outfile.write("\n")

for read_id in all_read_list:
    read_count += 1

    if read_count % read_bin == 0:
        outline = "\t".join([str(read_count), str(gene_count)])

        outfile.write(outline)
        outfile.write("\n")

    this_gene_id = read_gene_dict[read_id]

    if this_gene_id not in check_gene_dict:
        gene_count += 1
        check_gene_dict[this_gene_id] = 1





