#!/usr/bin/env python

import re
import sys
import time

import os
import argparse

"""
Transcriptome Annotation by Modular Algorithms (TAMA)

Author: Richard I. Kuo

Last changed: 2020/11/27

"""


#####################################################################
#####################################################################

ap = argparse.ArgumentParser(description='This script absorbs transcriptomes.')

ap.add_argument('-f', type=str, nargs=1, help='Bed file')
ap.add_argument('-o', type=str, nargs=1, help='Output file prefix')
ap.add_argument('-m', type=str, nargs=1, help='Exon ends threshold/ splice junction threshold (Default is 10)')
ap.add_argument('-e', type=str, nargs=1, help='Trans ends wobble threshold (Default is 500)')
ap.add_argument('-s', type=str, nargs=1, help='Single exon overlap percent threshold (Default is 20 percent)')
ap.add_argument('-id', type=str, nargs=1, help='Use original ID line original_id (Default is tama_id line based on gene_id;transcript_id structure')
ap.add_argument('-cds', type=str, nargs=1, help='Pull CDS option. Default is tama_cds where CDS regions matching TSS and TTS are ignored if another CDS is found. Use longest_cds to pick the longest CDS')


opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0

if not opts.f:
    print("Bed file missing")
    missing_arg_flag = 1
if not opts.o:
    print("Output file prefix missing")
    missing_arg_flag = 1


if not opts.m:
    print("Default exon end/splice junction threshold: 10")
    wobble_threshold = 10
else:
    wobble_threshold = int(opts.m[0])

if not opts.e:
    print("Trans ends wobble threshold (Default is 100)")
    ends_wobble_threshold = 500
else:
    ends_wobble_threshold = int(opts.e[0])

if not opts.s:
    print("Single exon overlap percent threshold (Default is 20 percent)")
    overlap_percent_threshold = 20
else:
    overlap_percent_threshold = int(opts.s[0])

if not opts.id:
    print("Single exon overlap percent threshold (Default is 20 percent)")
    id_use_flag = "tama_id"
else:
    id_use_flag = opts.id[0]

    if id_use_flag != "tama_id" and id_use_flag != "original_id":
        print("Error with ID input. Please use either tama_id or original_id")
        sys.exit()

if not opts.cds:
    print("Default is tama cds (ignore CDS predictions from TAMA Collapse)")
    cds_flag = "tama_cds"
else:
    cds_flag = opts.cds[0]

    if cds_flag != "tama_cds" and cds_flag != "longest_cds":
        print("Error with CDS input. Please use either tama_cds or longest_cds")
        sys.exit()

if missing_arg_flag == 1:
    print("Please try again with complete arguments")


bed_file = opts.f[0]
outfile_prefix = opts.o[0]





#####################################################################
#####################################################################


print("opening bed file")
bed_file_contents = open(bed_file).read().rstrip("\n").split("\n")

bed_outfile_name = outfile_prefix + ".bed"
outfile_bed = open(bed_outfile_name,"w")

discard_outfile_name = outfile_prefix + "_discarded.txt"
outfile_discard = open(discard_outfile_name,"w")


####################################################################################################



class Transcript:
    def __init__(self, trans_id):
        self.trans_id = trans_id
        self.uniq_trans_id = ""
        self.gene_id = ""
        self.merge_gene_id = ""
        self.scaffold = "none"
        self.trans_start = -1
        self.start_pos = self.trans_start
        self.trans_end = 0
        self.end_pos = self.trans_end
        self.exon_start_list = []
        self.exon_end_list = []

        self.id_line = ""

        self.cds_start = 0
        self.cds_end = 0
        self.num_exons = 0
        self.strand = ""
        self.source_id = ""



    def add_bed_info(self,bed_line):
        line_split = bed_line.split("\t")
        self.scaffold = line_split[0]
        self.trans_start = int(line_split[1])
        self.start_pos = self.trans_start
        self.trans_end = int(line_split[2])
        self.end_pos = self.trans_end
        self.id_line = line_split[3]
        id_line_split = self.id_line.split(";")
        self.gene_id = id_line_split[0]

        self.strand = line_split[5]
        self.num_exons = int(line_split[9])

        self.cds_start = int(line_split[6])
        self.cds_end = int(line_split[7])

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


    def make_exon_start_end_lines(self):

        exon_start_string_list = []
        exon_end_string_list = []

        for i in range(len(self.exon_start_list)):

            exon_start_string_list.append(str(self.exon_start_list[i]))
            exon_end_string_list.append(str(self.exon_end_list[i]))

        exon_start_string_line = ",".join(exon_start_string_list)
        exon_end_string_line = ",".join(exon_end_string_list)

        return exon_start_string_line,exon_end_string_line

    def format_bed_line(self):
        bed_list = []
        bed_list.append(str(self.scaffold))

        # correct for shifts in trans start due to 5' longer fragment model
        if self.exon_start_list[0] < self.trans_start:
            self.trans_start = self.exon_start_list[0]

        # correct for shifts in trans end due to 3' longer fragment model
        if self.exon_end_list[-1] > self.trans_end:
            self.trans_end = self.exon_end_list[-1]

        bed_list.append(str(self.trans_start))
        bed_list.append(str(self.trans_end))



        #gene_id = self.trans_id.split(".")[0]

        if id_use_flag == "tama_id":
            id_line = ";".join([self.gene_id,self.trans_id])
        elif id_use_flag == "original_id":
            id_line = self.id_line

        bed_list.append(str(id_line))
        bed_list.append("40")
        bed_list.append(self.strand)

        if cds_flag == "tama_cds":
            bed_list.append(str(self.trans_start))
            bed_list.append(str(self.trans_end))
        elif cds_flag == "longest_cds":
            bed_list.append(str(self.cds_start))
            bed_list.append(str(self.cds_end))
        else:
            print("Error with cds flag")
            print(cds_flag)
            sys.exit()

        bed_list.append("255,0,0")

        bed_list.append(str(self.num_exons))

        relative_exon_start_list = []
        exon_length_list = []
        for i in range(self.num_exons):
            exon_start = self.exon_start_list[i]
            exon_end = self.exon_end_list[i]
            exon_length = exon_end - exon_start
            relative_exon_start = exon_start - self.trans_start

            relative_exon_start_list.append(str(relative_exon_start))
            exon_length_list.append(str(exon_length))

        relative_exon_start_line = ",".join(relative_exon_start_list)
        exon_length_line = ",".join(exon_length_list)

        bed_list.append(exon_length_line)
        bed_list.append(relative_exon_start_line)

        bed_line = "\t".join(bed_list)

        return bed_line


####################################################################################################


#wobble_threshold = 10
#ends_wobble_threshold = 100

def compare_absorb_transcripts(a_trans_obj, b_trans_obj):

    a_exon_start_list = a_trans_obj.exon_start_list
    b_exon_start_list = b_trans_obj.exon_start_list

    a_exon_end_list = a_trans_obj.exon_end_list
    b_exon_end_list = b_trans_obj.exon_end_list

    a_num_exons = len(a_exon_start_list)
    b_num_exons = len(b_exon_start_list)

    #####
    a_trans_start = a_exon_start_list[0]
    b_trans_start = b_exon_start_list[0]

    a_trans_end = a_exon_end_list[-1]
    b_trans_end = b_exon_end_list[-1]

    a_trans_length = a_trans_end - a_trans_start
    b_trans_length = b_trans_end - b_trans_start

    a_cds_start = a_trans_obj.cds_start   #################################################################################AAAAAAAAAAAAA
    b_cds_start = b_trans_obj.cds_start

    a_cds_end = a_trans_obj.cds_end  #################################################################################AAAAAAAAAAAAA
    b_cds_end = b_trans_obj.cds_end

    final_cds_start = 0
    final_cds_end = 0

    if cds_flag == "tama_cds":
        # pick best CDS
        # Choose longest if both have CDS info
        if a_cds_start == a_trans_start and a_cds_end == a_trans_end:
            if b_cds_start == b_trans_start and b_cds_end == b_trans_end: #both seem to be just the TSS and TTS
                final_cds_start = 0
                final_cds_end = 0
            elif b_cds_start != b_trans_start or b_cds_end != b_trans_end:
                final_cds_start = b_cds_start
                final_cds_end = b_cds_end

        elif b_cds_start == b_trans_start and b_cds_end == b_trans_end:
            final_cds_start = a_cds_start
            final_cds_end = a_cds_end

        elif b_cds_start != b_trans_start or b_cds_end != b_trans_end:
            a_cds_length = a_cds_end - a_cds_start
            b_cds_length = b_cds_end - b_cds_start

            if a_cds_length >= b_cds_length:
                final_cds_start = a_cds_start
                final_cds_end = a_cds_end
            else:
                final_cds_start = b_cds_start
                final_cds_end = b_cds_end

    elif cds_flag == "longest_cds":
        a_cds_length = a_cds_end - a_cds_start
        b_cds_length = b_cds_end - b_cds_start

        if a_cds_length >= b_cds_length:
            final_cds_start = a_cds_start
            final_cds_end = a_cds_end
        else:
            final_cds_start = b_cds_start
            final_cds_end = b_cds_end




    #####

    # find index of matching SJ

    no_match_flag = "na"
    absorb_match_flag = "na"

    long_trans_id = "na"
    short_trans_id = "na"
    long_trans_obj = "na"


    # find long trans and short trans
    if a_num_exons == b_num_exons: # works for both single and both multi
        ##################
        # choose longest single exon transacript model
        if a_trans_length > b_trans_length:
            long_trans_obj = a_trans_obj
            short_trans_obj = b_trans_obj
        elif a_trans_length < b_trans_length:
            long_trans_obj = b_trans_obj
            short_trans_obj = a_trans_obj
        else:
            if a_trans_start <= b_trans_start:
                long_trans_obj = a_trans_obj
                short_trans_obj = b_trans_obj
            else:
                long_trans_obj = b_trans_obj
                short_trans_obj = a_trans_obj

        ####################
    else:
        if a_num_exons > b_num_exons:
            long_trans_obj = a_trans_obj
            short_trans_obj = b_trans_obj
        elif a_num_exons < b_num_exons:
            long_trans_obj = b_trans_obj
            short_trans_obj = a_trans_obj
        else:
            print("Error with number of exons!")
            sys.exit()

    ###### setup long short variables!! #################################

    long_num_exons = len(long_trans_obj.exon_start_list)
    short_num_exons = len(short_trans_obj.exon_start_list)

    long_trans_start = long_trans_obj.exon_start_list[0]
    short_trans_start = short_trans_obj.exon_start_list[0]

    long_trans_end = long_trans_obj.exon_end_list[-1]
    short_trans_end = short_trans_obj.exon_end_list[-1]

    long_trans_length = long_trans_end - long_trans_start
    short_trans_length = short_trans_end - short_trans_start

    long_index_start_match_dict = {} # long_index_start_match_dict[index] = 1
    short_index_start_match_dict = {}  # short_index_start_match_dict[index] = 1

    long_index_end_match_dict = {}  # long_index_end_match_dict[index] = 1
    short_index_end_match_dict = {}  # short_index_end_match_dict[index] = 1

    long_trans_id = long_trans_obj.trans_id
    short_trans_id = short_trans_obj.trans_id

    long_exon_start_list = long_trans_obj.exon_start_list
    long_exon_end_list = long_trans_obj.exon_end_list

    short_exon_start_list = short_trans_obj.exon_start_list
    short_exon_end_list = short_trans_obj.exon_end_list


    # in case both are single exon models #############################################################
    if long_num_exons == 1 and short_num_exons == 1:

        # find overlap percentage

        if long_trans_start <= short_trans_start:
            overlap_length = long_trans_end - short_trans_start
        elif long_trans_start > short_trans_start:
            overlap_length = short_trans_end - long_trans_start

        if overlap_length < 0:
            no_match_flag = "no_match"

        else:
            long_overlap_percent = (overlap_length * 100) // long_trans_length
            short_overlap_percent = (overlap_length * 100) // short_trans_length

            # combine if overlap meets threshold
            if long_overlap_percent > overlap_percent_threshold or short_overlap_percent > overlap_percent_threshold:

                no_match_flag = "trans_match"

                if long_trans_start < short_trans_start:
                    long_trans_start = long_trans_start
                else:
                    long_trans_start = short_trans_start

                long_trans_obj.exon_start_list[0] = long_trans_start

                if long_trans_end > short_trans_end:
                    long_trans_end = long_trans_end
                else:
                    long_trans_end = short_trans_end

                long_trans_obj.exon_end_list[-1] = long_trans_end

            else: # if overlap percent is not high enough this is not a match
                no_match_flag = "no_match"


    # in case one is a single exon model ####################################################################
    # find if shorter single exon model fits in larger multi exon model
    elif short_num_exons == 1:


        for exon_index in range(long_num_exons):
            this_exon_start = long_exon_start_list[exon_index]
            this_exon_end = long_exon_end_list[exon_index]

            # check for exon overlap
            if short_trans_start < this_exon_end and short_trans_end > this_exon_start:
                start_wobble = this_exon_start - short_trans_start
                end_wobble = short_trans_end - this_exon_end

                # check if this is first exon
                if exon_index == 0:
                    start_wobble_threshold = ends_wobble_threshold
                    end_wobble_theshold = wobble_threshold

                    if start_wobble < start_wobble_threshold:
                        if end_wobble < end_wobble_theshold:

                            no_match_flag = "trans_match"

                            if start_wobble > 0:
                                long_trans_obj.exon_start_list[0] = short_trans_start #use the earliest start
                                long_trans_obj.trans_start = short_trans_start

                elif exon_index == long_num_exons-1:
                    start_wobble_threshold = wobble_threshold
                    end_wobble_theshold = ends_wobble_threshold

                    if start_wobble < start_wobble_threshold:
                        if end_wobble < end_wobble_theshold:

                            no_match_flag = "trans_match"

                            if end_wobble > 0:
                                long_trans_obj.exon_end_list[-1] = short_trans_end # use the latest end
                                long_trans_obj.trans_end = short_trans_end ###############################################################################################
                else:
                    start_wobble_threshold = wobble_threshold
                    end_wobble_theshold = wobble_threshold

                    if start_wobble < start_wobble_threshold:
                        if end_wobble < end_wobble_theshold:

                            no_match_flag = "trans_match"
        # for single exon one i reverse the logic on the flags so i need this hear to set all to no match
        if no_match_flag != "trans_match":
            no_match_flag = "no_match"

    elif long_num_exons > 1 and short_num_exons > 1: # if both are multi exon models! ########################################
        for i in range(long_num_exons):
            for j in range(short_num_exons):


                long_this_exon_start = long_exon_start_list[i]
                short_this_exon_start = short_exon_start_list[j]

                long_this_exon_end = long_exon_end_list[i]
                short_this_exon_end = short_exon_end_list[j]


                # at trans start for both
                if i == 0 and j == 0:


                    # check that exons overlap
                    if long_this_exon_start <= short_this_exon_end and long_this_exon_end >= short_this_exon_start:

                        # no abs to account for TSS variability
                        this_start_wobble = long_this_exon_start - short_this_exon_start

                        if this_start_wobble <= ends_wobble_threshold:
                            long_index_start_match_dict[i] = 1
                            short_index_start_match_dict[j] = 1

                        ################

                        # abs because this is a splice junction
                        this_end_wobble = abs(short_this_exon_end - long_this_exon_end)

                        if this_end_wobble <= wobble_threshold:
                            long_index_end_match_dict[i] = 1
                            short_index_end_match_dict[j] = 1

                elif i > 0 and j == 0: # start only for j

                    # check that exons overlap
                    if long_this_exon_start <= short_this_exon_end and long_this_exon_end >= short_this_exon_start:

                        # no abs to account for fragment end variability
                        this_start_wobble = long_this_exon_start - short_this_exon_start

                        if this_start_wobble <= wobble_threshold:
                            long_index_start_match_dict[i] = 1
                            short_index_start_match_dict[j] = 1

                        ################
                        # abs because this is a splice junction
                        this_end_wobble = abs(short_this_exon_end - long_this_exon_end)

                        if this_end_wobble <= wobble_threshold:
                            long_index_end_match_dict[i] = 1
                            short_index_end_match_dict[j] = 1

                elif i == long_num_exons-1 and j == short_num_exons-1: # at the end of both

                    # check that exons overlap
                    if long_this_exon_start <= short_this_exon_end and long_this_exon_end >= short_this_exon_start:

                        # abs because this is a splice junction
                        this_start_wobble = abs(long_this_exon_start - short_this_exon_start)

                        if this_start_wobble <= wobble_threshold:
                            long_index_start_match_dict[i] = 1
                            short_index_start_match_dict[j] = 1

                        ################
                        # no abs to account for fragment end variability
                        this_end_wobble = short_this_exon_end - long_this_exon_end

                        if this_end_wobble <= ends_wobble_threshold:
                            long_index_end_match_dict[i] = 1
                            short_index_end_match_dict[j] = 1

                elif i < long_num_exons-1 and j == short_num_exons-1: # at the end of short

                    # check that exons overlap
                    if long_this_exon_start <= short_this_exon_end and long_this_exon_end >= short_this_exon_start:

                        # abs because this is a splice junction
                        this_start_wobble = abs(long_this_exon_start - short_this_exon_start)

                        if this_start_wobble <= wobble_threshold:
                            long_index_start_match_dict[i] = 1
                            short_index_start_match_dict[j] = 1

                        ################
                        # no abs to account for fragment end variability
                        this_end_wobble = short_this_exon_end - long_this_exon_end

                        if this_end_wobble <= wobble_threshold: # use sj wobble because this is sj for long
                            long_index_end_match_dict[i] = 1
                            short_index_end_match_dict[j] = 1
                else:

                    # check that exons overlap
                    if long_this_exon_start <= short_this_exon_end and long_this_exon_end >= short_this_exon_start:

                        # abs because this is a splice junction
                        this_start_wobble = abs(long_this_exon_start - short_this_exon_start)

                        if this_start_wobble <= wobble_threshold:
                            long_index_start_match_dict[i] = 1
                            short_index_start_match_dict[j] = 1

                        ################
                            # abs because this is a splice junction
                        this_end_wobble = abs(short_this_exon_end - long_this_exon_end)

                        if this_end_wobble <= wobble_threshold: # use sj wobble
                            long_index_end_match_dict[i] = 1
                            short_index_end_match_dict[j] = 1




        long_index_start_list = list(long_index_start_match_dict.keys())
        long_index_start_list.sort()

        short_index_start_list = list(short_index_start_match_dict.keys())
        short_index_start_list.sort()

        long_index_end_list = list(long_index_end_match_dict.keys())
        long_index_end_list.sort()

        short_index_end_list = list(short_index_end_match_dict.keys())
        short_index_end_list.sort()

        ##########################################
        #print(long_index_start_list)
        #print(long_index_end_list)
        #print(short_index_start_list)
        #print(short_index_end_list)

        # if the number of start and end matches are not equal, this is not a match
        if len(long_index_start_list) != len(long_index_end_list):
            no_match_flag = "no_match"

        # if the number of start and end matches are not equal, this is not a match
        if len(short_index_start_list) != len(short_index_end_list):
            no_match_flag = "no_match"

        # no matches
        if len(long_index_start_list) == 0:
            no_match_flag = "no_match"

        # no matches
        if len(short_index_start_list) == 0:
            no_match_flag = "no_match"


        # check that all short SJ and outer ends match
        if len(short_index_start_list) != short_num_exons:
            no_match_flag = "no_match"
        if len(short_index_end_list) != short_num_exons:
            no_match_flag = "no_match"


        if no_match_flag != "no_match":

            # check that indices are matching for SJ####################################
            # check for consecutive indices ######################################
            # for long
            prev_index_start = -1
            prev_index_end = -1
            for i in range(len(long_index_start_list)):
                long_index_start = long_index_start_list[i]
                long_index_end = long_index_end_list[i]

                if long_index_start - long_index_end != 0:
                    no_match_flag = "no_match"

                if i > 0 :
                    if long_index_start-prev_index_start != 1:
                        no_match_flag = "no_match"
                    if long_index_end-prev_index_end != 1:
                        no_match_flag = "no_match"

                prev_index_start = long_index_start
                prev_index_end = long_index_end


            # check that indices are matching for SJ
            # for short
            prev_index_start = -1
            prev_index_end = -1
            for i in range(len(short_index_start_list)):
                short_index_start = short_index_start_list[i]
                short_index_end = short_index_end_list[i]

                if short_index_start - short_index_end != 0:
                    no_match_flag = "no_match"

                if i > 0 :
                    if short_index_start-prev_index_start != 1:
                        no_match_flag = "no_match"
                    if short_index_end-prev_index_end != 1:
                        no_match_flag = "no_match"

                prev_index_start = short_index_start
                prev_index_end = short_index_end

            # if matching find longest model
            if no_match_flag != "no_match":

                # check short TSS and TES to make sure it fits in the longer model

                # get index of a exon start right before sj match
                long_sj_start_match_index = long_index_end_list[0]
                short_trans_start = short_exon_start_list[0]

                # dont need this anymore as I take care of in the above section
                #if long_exon_start_list[long_sj_start_match_index] - short_trans_start > ends_wobble_threshold:
                #    no_match_flag = "no_match"

                ####################
                # match with one short and one long
                if long_num_exons > short_num_exons:
                    # check that all SJ in shorter model are matching ##################
                    num_short_match_start = len(short_index_start_list)
                    if num_short_match_start < short_num_exons:
                        no_match_flag = "no_match"


                else: # when both have same number of exons

                    ##################
                    # check that all SJ in both models are matching ##################
                    num_long_match_starts = len(long_index_start_list)

                    if num_long_match_starts < long_num_exons:
                        no_match_flag = "no_match"

                    # check that all SJ in both models are matching ##################
                    num_short_match_starts = len(short_index_start_list)

                    if num_short_match_starts < short_num_exons:
                        no_match_flag = "no_match"

                    #################

                    if no_match_flag != "no_match":

                        if long_trans_start < short_trans_start:
                            long_trans_start = long_trans_start
                        else:
                            long_trans_start = short_trans_start

                        long_trans_obj.exon_start_list[0] = long_trans_start

                        if long_trans_end > short_trans_end:
                            long_trans_end = long_trans_end
                        else:
                            long_trans_end = short_trans_end

                        long_trans_obj.exon_end_list[-1] = long_trans_end


    else:

        print("Error with comparing transcripts!!")
        sys.exit()

        #absorb_match_flag = no_match_flag
        #long_trans_id = "na"
        #short_trans_id = "na"
        #long_trans_obj = "na"


    if no_match_flag == "no_match":
        absorb_match_flag = no_match_flag
    else:
        absorb_match_flag = "trans_match"

        if final_cds_start == 0 and final_cds_end == 0: # if both seem to have TSS and TTS for CDS then just use new TSS and TTS
            long_trans_obj.cds_start = long_trans_obj.trans_start
            long_trans_obj.cds_end = long_trans_obj.trans_end
        else:
            long_trans_obj.cds_start = final_cds_start
            long_trans_obj.cds_end = final_cds_end


    return absorb_match_flag,long_trans_id, short_trans_id, long_trans_obj

############################################################################################

bed_dict = {} # bed_dict[scaffold][start][end][strand] = list of bed lines with added source id

gene_trans_obj_dict = {} # gene_trans_obj_dict[gene id][trans id] = Trans Obj
gene_trans_order_dict = {} # gene_trans_order_dict[gene_id] = list of trans id
gene_list = []

for line in bed_file_contents:
    line_split = line.split("\t")

    scaffold = line_split[0]
    trans_start = int(line_split[1])
    trans_end = int(line_split[2])
    id_line = line_split[3]

    id_line_split = id_line.split(";")

    if len(id_line_split) < 2:
        print("Error with bed file ID field")
        print(filename)
        print(line)
        print("bed12 files must have the gene ID's and transcript ID's formatted as such \"gene_id;transcript_id\" in the 4th column.")
        print("The gene ID must be the first subfield and the subfields must be delimited with a semicolon (;).")
        sys.exit()

    strand = line_split[5]
    num_exon = int(line_split[9])
    block_size_list = line_split[10].split(",")
    block_start_list = line_split[11].split(",")

    gene_id = id_line_split[0]
    trans_id = id_line_split[1]

    trans_obj = Transcript(trans_id)
    trans_obj.add_bed_info(line)

    if gene_id not in gene_trans_obj_dict:
        gene_trans_obj_dict[gene_id] = {}
        gene_list.append(gene_id)

        gene_trans_order_dict[gene_id] = []

    gene_trans_obj_dict[gene_id][trans_id] = trans_obj

    gene_trans_order_dict[gene_id].append(trans_id)


for gene_id in gene_list:


    #trans_id_list = list(gene_trans_obj_dict[gene_id].keys())

    trans_id_list = gene_trans_order_dict[gene_id]

    remove_trans_dict = {} # remove_trans_dict[trans_id] = 1

    # check each transcript against each other

    for a_trans_id in trans_id_list:
        for b_trans_id in trans_id_list:

            if a_trans_id == b_trans_id:
                continue

            if a_trans_id in remove_trans_dict:
                continue

            if b_trans_id in remove_trans_dict:
                continue


            a_trans_obj = gene_trans_obj_dict[gene_id][a_trans_id]
            b_trans_obj = gene_trans_obj_dict[gene_id][b_trans_id]

            absorb_match_flag, long_trans_id, short_trans_id, long_trans_obj = compare_absorb_transcripts(a_trans_obj,b_trans_obj)

            ########################################################################################
            #print(absorb_match_flag + "\t"+ long_trans_id  + "\t"+short_trans_id )
            #if gene_id == "G17":
            #    sys.exit()
            ##########################################################################################


            if absorb_match_flag == "trans_match":
                remove_trans_dict[short_trans_id] = 1

                if long_trans_id == a_trans_id:
                    gene_trans_obj_dict[gene_id][a_trans_id] = long_trans_obj
                elif long_trans_id == b_trans_id:
                    gene_trans_obj_dict[gene_id][b_trans_id] = long_trans_obj


    for this_trans_id in trans_id_list:
        if this_trans_id in remove_trans_dict:
            dis_trans_obj = gene_trans_obj_dict[gene_id][this_trans_id]
            dis_bed_line = dis_trans_obj.format_bed_line()

            outfile_discard.write(dis_bed_line)
            outfile_discard.write("\n")

            continue

        this_trans_obj = gene_trans_obj_dict[gene_id][this_trans_id]
        this_bed_line = this_trans_obj.format_bed_line()

        outfile_bed.write(this_bed_line)
        outfile_bed.write("\n")






    ####################################################################################################






####################################################################################################