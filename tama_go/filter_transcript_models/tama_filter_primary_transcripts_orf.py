import re
import sys
import time

import argparse


#
# This script uses the ORF/NMD output bed file and filters to have only 1 transcript per gene
#
#



ap = argparse.ArgumentParser(description='This script uses the ORF/NMD output bed file and filters to have only 1 transcript per gene')

ap.add_argument('-b', type=str, nargs=1, help='bed file (required)')
ap.add_argument('-o', type=str, nargs=1, help='Output file name (required)')


opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0

if not opts.b:
    print("bed file missing")
    missing_arg_flag = 1
if not opts.o:
    print("output name missing")
    missing_arg_flag = 1


if missing_arg_flag == 1:
    print("Please try again with complete arguments")

bed_file = opts.b[0]
outfile_name = opts.o[0]


print("opening bed file")
#blastp_file = sys.argv[1]
bed_file_contents = open(bed_file).read().rstrip("\n").split("\n")

#outfile_name = sys.argv[2]
outfile = open(outfile_name,"w")



match_quality_dict = {} # match_quality_dict[match qual]  = score
match_quality_dict["full_match"] = 900000
match_quality_dict["90_match"] = 800000
match_quality_dict["50_match"] = 700000
match_quality_dict["bad_match"] = 600000
#match_quality_dict["no_hit"] = 500000
#match_quality_dict["no_orf"] = 400000


class Transcript:
    def __init__(self, trans_id,bed_line):
        self.trans_id = trans_id
        self.gene_id = ""

        self.scaffold = "none"
        self.trans_start = -1
        self.trans_end = 0
        self.exon_start_list = []
        self.exon_end_list = []

        self.id_line = ""

        self.num_exons = 0
        self.strand = ""

        self.prot_id = ""
        self.match_length = ""
        self.match_quality = ""
        self.nmd_flag = ""
        self.frame_flag = ""

        self.bed_line = bed_line

        self.trans_score = 0

        self.trans_length = 0

    def add_bed_info(self,bed_line):
        line_split = bed_line.split("\t")
        self.scaffold = line_split[0]
        self.trans_start = int(line_split[1])
        self.trans_end = int(line_split[2])
        self.id_line = line_split[3]
        id_line_split = self.id_line.split(";")
        self.gene_id = id_line_split[0]

        self.strand = line_split[5]
        self.num_exons = int(line_split[9])

        self.trans_length = self.trans_end - self.trans_start


        block_size_list = line_split[10].split(",")
        block_start_list = line_split[11].split(",")

        #add this for commas at the end of the exon block and start strings
        if block_size_list[-1] == "":
            block_size_list.pop(-1)

        if block_start_list[-1] == "":
            block_start_list.pop(-1)

        for i in range(len(block_size_list)):
            rel_exon_start = int(block_start_list[i])
            rel_exon_end = rel_exon_start + int(block_size_list[i])

            exon_start = self.trans_start + rel_exon_start
            exon_end = self.trans_start + rel_exon_end

            self.exon_end_list.append(exon_end)
            self.exon_start_list.append(exon_start)


    def id_parser(self,id_line):

        # for ORF NMD output
        # 1       14361   29359   G1;G1.1;UniRef50_B4DXR4;full_length;full_match;NMD3;F3  40      -
        # 16747   17310   200,0,255       10      468,69,152,159,198,510,147,112,154,39   0,608,1434,2245,2496,2871,3553,3906,10376,14959


        id_split = id_line.split(";")

        tama_gene_id = id_split[0]
        tama_trans_id = id_split[1]

        self.prot_id = id_split[2]
        self.match_length = id_split[3]
        self.match_quality = id_split[4]
        self.nmd_flag = id_split[5]
        self.frame_flag = id_split[6]

        if self.match_length == "full_length":
            self.trans_score =  self.trans_score + 1000000

        if self.match_quality in match_quality_dict:
            self.trans_score =  self.trans_score + match_quality_dict[self.match_quality]

        if self.nmd_flag == "prot_ok":
            self.trans_score = self.trans_score + 10000

def select_best(trans_dict):

    select_trans_dict = {} # select_trans_dict[trans_id] = 1
    trans_list = []
    score_trans_dict = {} # score_trans_dict[score][trans id] = length

    for this_trans_id in trans_dict:
        select_trans_dict[this_trans_id] = 1
        trans_list.append(this_trans_id)

        this_trans_obj =  trans_dict[this_trans_id]

        if this_trans_obj.trans_score not in score_trans_dict:
            score_trans_dict[this_trans_obj.trans_score] = {}

        score_trans_dict[this_trans_obj.trans_score][this_trans_id] = this_trans_obj.trans_length


    score_list = list(score_trans_dict.keys())
    score_list.sort()

    high_score = score_list[-1]

    num_high_score_trans = len(list(score_trans_dict[high_score].keys()))

    length_trans_dict = {} # length_trans_dict[trans_length][trans id] = 1

    for this_trans_id in score_trans_dict[high_score]:

        this_trans_length = score_trans_dict[high_score][this_trans_id]

        if this_trans_length not in length_trans_dict:
            length_trans_dict[this_trans_length] = {}

        length_trans_dict[this_trans_length][this_trans_id] = 1


    length_list = list(length_trans_dict.keys())
    length_list.sort()

    longest_length = length_list[-1]

    longest_trans_list = list(length_trans_dict[longest_length].keys())
    num_longest = len(longest_trans_list)

    longest_trans_list.sort()

    best_trans_id = longest_trans_list[0]

    return best_trans_id




gene_trans_dict = {} # gene_trans_dict[gene_id][trans_id] = trans_obj
gene_list = []

for line in bed_file_contents:
    
    line_split = line.split("\t")
    chrom = line_split[0]
    t_start = line_split[1]
    t_end = line_split[2]
    id_line = line_split[3]

    id_split = id_line.split(";")

    gene_id = id_split[0]
    trans_id = id_split[1]

    trans_obj = Transcript(trans_id,line)

    if gene_id not in gene_trans_dict:
        gene_trans_dict[gene_id] = {}
        gene_list.append(gene_id)

    trans_obj.add_bed_info(line)
    trans_obj.id_parser(id_line)


    gene_trans_dict[gene_id][trans_id] = trans_obj


for gene_id in gene_list:

    trans_list = list(gene_trans_dict[gene_id].keys())

    trans_dict = gene_trans_dict[gene_id]

    best_trans_id = select_best(trans_dict)

    new_bed_line = gene_trans_dict[gene_id][best_trans_id].bed_line

    outfile.write(new_bed_line)
    outfile.write("\n")




