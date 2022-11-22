import re
import sys
import time

import argparse


#
# This script makes the Ensembl ID's the primary ID's and allows for filtering
#
#



ap = argparse.ArgumentParser(description='This script makes the Ensembl IDs the primary IDs and allows for filtering')

ap.add_argument('-b', type=str, nargs=1, help='bed file (required)')
ap.add_argument('-o', type=str, nargs=1, help='Output file name (required)')
ap.add_argument('-f', type=str, nargs=1, help='Filter level (default "none", use "only_match" to only include models with a match)')
ap.add_argument('-s', type=str, nargs=1, help='Sub-field management method (default "ensembl_merge" for restructuring sub-fields from Ensembl ID, use "custom" to define sub-field shuffling) ')
ap.add_argument('-r', type=str, nargs=1, help='Sub-field reshuffle parameter (default "none")')
ap.add_argument('-d', type=str, nargs=1, help='Sub-field reshuffle delimiters (default ";")')


opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0

if not opts.b:
    print("bed file missing")
    missing_arg_flag = 1
if not opts.o:
    print("output name missing")
    missing_arg_flag = 1

if not opts.f:
    filter_level = "none"
else:
    filter_level = opts.f[0]


if not opts.s:
    subfield_method = "ensembl_merge"
else:
    subfield_method = opts.s[0]


if not opts.r:
    reshuffle_params = "none"
else:
    reshuffle_params = opts.r[0]

if not opts.d:
    delim_line = ";"
else:
    delim_line = opts.d[0]


if missing_arg_flag == 1:
    print("Please try again with complete arguments")

bed_file = opts.b[0]
outfile_name = opts.o[0]


print("opening bed file")
#blastp_file = sys.argv[1]
bed_file_contents = open(bed_file).read().rstrip("\n").split("\n")

#outfile_name = sys.argv[2]
outfile = open(outfile_name,"w")


def id_parser(id_line):

    if subfield_method == "ensembl_merge":

        output_flag = "pass"

        delim_string = delim_line

        id_split = id_line.split(delim_string)

        tama_gene_id = id_split[0]
        tama_trans_id = id_split[1]

        ens_gene_id = ""
        ens_trans_id = ""
        other_id_list = []

        if len(id_split) > 3:
            ens_gene_id = id_split[2]
            ens_trans_id = id_split[3]

            # check that the right parameters are being used
            if "ENS" not in ens_gene_id:
                print("Ensembl gene ID is not right. You may be using the wrong parameters. Please see wiki for manual.")
                print(ens_gene_id)

                new_line = "Ensembl gene ID is not right. You may be using the wrong parameters. Please see wiki for manual."
                outfile.write(new_line)
                outfile.write("\n")

                outfile.write(ens_gene_id)
                outfile.write("\n")



                sys.exit()

            # check that the right parameters are being used
            if "ENS" not in ens_trans_id:
                print("Ensembl transcript ID is not right. You may be using the wrong parameters. Please see wiki for manual.")
                print(ens_trans_id)

                new_line = "Ensembl transcript ID is not right. You may be using the wrong parameters. Please see wiki for manual."
                outfile.write(new_line)
                outfile.write("\n")

                outfile.write(ens_trans_id)
                outfile.write("\n")

                sys.exit()

            other_id_list = id_split[4:]

        if filter_level == "only_match":
            if "ENS" not in ens_gene_id:
                output_flag = "fail"


        new_output_list = []
        new_output_list.append(ens_gene_id)
        new_output_list.append(ens_trans_id)
        new_output_list.append(tama_gene_id)
        new_output_list.append(tama_trans_id)

        new_output_list.extend(other_id_list)

        new_output_line = ";".join(new_output_list)

    elif subfield_method == "ensembl_orf":

        output_flag = "pass"

        delim_string = delim_line

        id_split = id_line.split(delim_string)

        tama_gene_id = id_split[0]
        tama_trans_id = id_split[1]

        ens_gene_id = ""
        ens_trans_id = ""
        other_id_list = []

        if len(id_split) > 3:
            ens_id_line = id_split[2]

            ens_id_split = ens_id_line.split(",")

            if len(ens_id_split) > 1:

                ens_gene_id = ens_id_split[0]
                ens_trans_id = ens_id_split[1]

                # check that the right parameters are being used
                if "ENS" not in ens_gene_id:
                    print("Ensembl gene ID is not right. You may be using the wrong parameters. Please see wiki for manual.")
                    print(ens_gene_id)

                    new_line = "Ensembl gene ID is not right. You may be using the wrong parameters. Please see wiki for manual."
                    outfile.write(new_line)
                    outfile.write("\n")

                    outfile.write(ens_gene_id)
                    outfile.write("\n")

                    sys.exit()

                # check that the right parameters are being used
                if "ENS" not in ens_trans_id:
                    print("Ensembl transcript ID is not right. You may be using the wrong parameters. Please see wiki for manual.")
                    print(ens_trans_id)

                    new_line = "Ensembl transcript ID is not right. You may be using the wrong parameters. Please see wiki for manual."
                    outfile.write(new_line)
                    outfile.write("\n")

                    outfile.write(ens_trans_id)
                    outfile.write("\n")

                    sys.exit()


            else:
                ens_gene_id = tama_gene_id
                ens_trans_id = tama_trans_id

            other_id_list = id_split[3:]

        if filter_level == "only_match":
            if "ENS" not in  ens_gene_id:
                output_flag = "fail"


        new_output_list = []
        new_output_list.append(ens_gene_id)
        new_output_list.append(ens_trans_id)
        new_output_list.append(tama_gene_id)
        new_output_list.append(tama_trans_id)

        new_output_list.extend(other_id_list)

        new_output_line = ";".join(new_output_list)


    elif subfield_method == "custom":

        output_flag = "pass"

        delim_list = []

        for i in range(len(delim_line)):
            delim_list.append(delim_line[i])

        delim_string = "|".join(delim_list)

        #print(delim_string)

        # id_split = id_line.split(delim_string)

        id_split = re.split(delim_string,id_line)

        # print(id_split)

        # example reshuffle_params 3,4,1,2 this means the order of the new subfields
        reshuffle_list = reshuffle_params.split(",")

        new_output_list = []

        for reshufle_index in reshuffle_list:
            reshufle_index = int(reshufle_index) - 1

            new_output_list.append(id_split[reshufle_index])

        new_output_line = ";".join(new_output_list)





    return new_output_line,output_flag



for line in bed_file_contents:
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

    new_output_line,output_flag = id_parser(id_line)
    new_line_split = line_split

    if output_flag == "pass":
        new_line_split[3] = new_output_line

        new_line = "\t".join(new_line_split)

        outfile.write(new_line)
        outfile.write("\n")



