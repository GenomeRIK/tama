#!/usr/bin/env python

import re
import sys
import time
from Bio import SeqIO
from StringIO import StringIO
from Bio import AlignIO
import os
import argparse

from __init__ import __version__

"""
Transcriptome Annotation by Modular Algorithms (TAMA)
TAMA Merge

Author: Richard I. Kuo

This script merges transcriptome/genome annotations.


"""

tm_version = 'tm0.0'


start_time = time.time()
prev_time = start_time


#####################################################################
#####################################################################

ap = argparse.ArgumentParser(description='This script merges transcriptomes.')

ap.add_argument('-f', type=str, nargs=1, help='File list')
ap.add_argument('-p', type=str, nargs=1, help='Output prefix')

ap.add_argument('-e', type=str, nargs=1, help='Collapse exon ends flag: common_ends or longest_ends  (Default is common_ends)')

ap.add_argument('-a', type=str, nargs=1, help='5 prime threshold (Default is 10)')
ap.add_argument('-m', type=str, nargs=1, help='Exon ends threshold/ splice junction threshold (Default is 10)')
ap.add_argument('-z', type=str, nargs=1, help='3 prime threshold (Default is 10)')

opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0

if not opts.f:
    print("Fasta file missing")
    missing_arg_flag = 1
if not opts.p:
    print("Output prefix name missing")
    missing_arg_flag = 1
  

if not opts.e:
    print("Default collapse exon ends flag will be used: common_ends")
    collapse_flag = "common_ends"
else:
    collapse_flag = opts.e[0]
    
if not opts.a:
    print("Default 5 prime threshold: 20")
    fiveprime_threshold = 20
else:
    fiveprime_threshold = int(opts.a[0])
    
if not opts.m:
    print("Default exon end/splice junction threshold: 10")
    exon_diff_threshold = 10
else:
    exon_diff_threshold = int(opts.m[0])

if not opts.z:
    print("Default 3 prime threshold: 20")
    threeprime_threshold = 20
else:
    threeprime_threshold = int(opts.z[0])

if missing_arg_flag == 1:
    print("Please try again with complete arguments")


filelist_file = opts.f[0]
outfile_prefix = opts.p[0]

#####################################################################
#####################################################################


print("opening file list")
#filelist_file = sys.argv[1]
filelist_file_contents = open(filelist_file).read().rstrip("\n").split("\n")

#outfile_prefix = sys.argv[2]


bed_outfile_name = outfile_prefix + ".bed"
outfile_bed = open(bed_outfile_name,"w")

trans_report_outfile_name = outfile_prefix + "_trans_report.txt"
outfile_trans_report = open(trans_report_outfile_name,"w")
trans_report_line = "\t".join(["transcript_id","num_clusters","sources","start_wobble_list","end_wobble_list","exon_start_support","exon_end_support"])
outfile_trans_report.write(trans_report_line)
outfile_trans_report.write("\n")

gene_report_outfile_name = outfile_prefix + "_gene_report.txt"
outfile_gene_report = open(gene_report_outfile_name,"w")
gene_report_line = "\t".join(["gene_id","num_clusters","num_final_trans","sources","chrom", "start","end"])
outfile_gene_report.write(gene_report_line)
outfile_gene_report.write("\n")

merge_outfile_name = outfile_prefix + "_merge.txt"
outfile_merge = open(merge_outfile_name,"w")
#merge_line = "\t".join(["transcript_id","cluster_id","scaffold","strand","start","end","exon_starts","exon_ends"])
#outfile_merge.write(merge_line)
#outfile_merge.write("\n")



bed_dict = {} # bed_dict[scaffold][start][end][strand] = list of bed lines with added source id
source_dict = {} # source_dict[source id][filename, seq type, priority rank] = value

   

def track_time(start_time,prev_time):
    end_time = time.time()
    time_taken = int(end_time - prev_time)
    tt_hours = time_taken / 60
    tt_hours = tt_hours /60
    leftover_time = time_taken - (tt_hours * 60 * 60)
    tt_minutes = leftover_time / 60
    leftover_time = time_taken - (tt_minutes * 60)
    tt_seconds = leftover_time
    print("time taken since last check:\t" +  str(tt_hours) + ":" + str(tt_minutes) + ":" + str(tt_seconds) )
    
    time_total = int(end_time - start_time)
    tt_hours = time_total / 60
    tt_hours = tt_hours /60
    leftover_time = time_total - (tt_hours * 60 * 60)
    tt_minutes = leftover_time / 60
    leftover_time = time_total - (tt_minutes * 60)
    tt_seconds = leftover_time
    print("time taken since beginning:\t" +  str(tt_hours) + ":" + str(tt_minutes) + ":" + str(tt_seconds) )
    
    this_time = end_time
    
    return this_time


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
        
        self.start_priority = 0
        self.end_priority = 0
        self.junction_priority = 0
    
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
        
        for i in xrange(len(block_size_list)):
            rel_exon_start = int(block_start_list[i])
            rel_exon_end = rel_exon_start + int(block_size_list[i])
            
            exon_start = self.trans_start + rel_exon_start
            exon_end = self.trans_start + rel_exon_end
            
            self.exon_end_list.append(exon_end)
            self.exon_start_list.append(exon_start)
            
            
    
    def add_source_id(self,source_id):
        self.source_id = source_id
        self.uniq_trans_id = "_".join([source_id,self.trans_id])
    
    def add_priority(self,start_priority,junction_priority,end_priority):
        self.start_priority = int(start_priority)
        self.end_priority = int(end_priority)
        self.junction_priority = int(junction_priority)
    
    
    def make_exon_start_end_lines(self):
        
        exon_start_string_list = []
        exon_end_string_list = []
        
        for i in xrange(len(self.exon_start_list)):
            
            exon_start_string_list.append(str(self.exon_start_list[i]))
            exon_end_string_list.append(str(self.exon_end_list[i]))
        
        exon_start_string_line = ",".join(exon_start_string_list)
        exon_end_string_line = ",".join(exon_end_string_list)
        
        return exon_start_string_line,exon_end_string_line
        
    def format_bed_line(self,final_trans_id):
        bed_list = []
        bed_list.append(str(self.scaffold))
        bed_list.append(str(self.trans_start))
        bed_list.append(str(self.trans_end))
        
        #gene_id = self.trans_id.split(".")[0]
        id_line = ";".join([final_trans_id,self.uniq_trans_id])
        
        bed_list.append(str(id_line))
        bed_list.append("40")
        bed_list.append(self.strand)
        
        bed_list.append(str(self.trans_start))
        bed_list.append(str(self.trans_end))
        
        bed_list.append("255,0,0")
        
        bed_list.append(str(self.num_exons))
        
        relative_exon_start_list = []
        exon_length_list = []
        for i in xrange(self.num_exons):
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

rgb_colour_dict = {} # rgb_colour_dict[number] = rgb values ie (255,0,0)
rgb_colour_dict[1] = "255,0,0" #red
rgb_colour_dict[2] = "255,100,0" #orange
rgb_colour_dict[3] = "255,200,0" #yellow
rgb_colour_dict[4] = "200,255,0" #lime
rgb_colour_dict[5] = "0,255,200" # light turquoise 
rgb_colour_dict[6] = "0,200,255" # light blue
rgb_colour_dict[7] = "0,100,255" # royal blue
rgb_colour_dict[8] = "0,0,255" #dark blue
rgb_colour_dict[9] = "100,0,255" # dark purple
rgb_colour_dict[10] = "200,0,255" # magenta

class Merged:
    def __init__(self, trans_id):
        self.trans_id = trans_id
        
        #same as collapse start and end list but used for flexibility in calling
        self.exon_start_list = [] 
        self.exon_end_list = []
        
        self.scaffold = "none"
        
        self.merged_trans_dict = {} # merged_trans_dict[source_trans id] = trans obj
        self.trans_list = []
        self.trans_obj_list = []
        
        self.strand = "none"
        self.num_trans = 0
        
        self.collapse_start_list = []
        self.collapse_end_list = []
        self.start_wobble_list = [] 
        self.end_wobble_list = []
        
        self.start_pos = 0
        self.end_pos = 0
        
        self.num_exons = 0
        
        #use these to see which transcripts support which coordinates
        self.e_start_trans_dict = {}
        self.e_end_trans_dict = {}
        
        self.e_start_trans_list = []
        self.e_end_trans_list = []

    
    def add_merged_trans(self,trans_obj):
        merged_trans_id = trans_obj.uniq_trans_id
        
        self.trans_list.append(merged_trans_id)
        self.trans_obj_list.append(trans_obj)
        
        self.merged_trans_dict[merged_trans_id] = trans_obj
        
        merged_trans_id_list = self.merged_trans_dict.keys()
        self.num_trans = len(merged_trans_id_list)
        
        if self.num_exons < int(trans_obj.num_exons):
            self.num_exons = int(trans_obj.num_exons)
        
        self.scaffold = trans_obj.scaffold
        
        if self.strand == "none":
            self.strand = trans_obj.strand
        elif self.strand != trans_obj.strand:
            print("Error with merged trans not on the same strand")
            sys.exit()

    def add_merge_info(self,collapse_start_list,collapse_end_list,start_wobble_list,end_wobble_list,e_start_trans_dict,e_end_trans_dict ):
        self.collapse_start_list = collapse_start_list
        self.collapse_end_list = collapse_end_list
        self.start_wobble_list = start_wobble_list
        self.end_wobble_list = end_wobble_list
        
        self.exon_start_list = collapse_start_list
        self.exon_end_list = collapse_end_list

        self.start_pos = collapse_start_list[0]
        self.end_pos = collapse_end_list[-1]
        
        
        self.e_start_trans_dict = e_start_trans_dict
        self.e_end_trans_dict = e_end_trans_dict
        
        for i in xrange(len(collapse_start_list)):
            e_start = collapse_start_list[i]
            e_end = collapse_end_list[i]
            
            e_start_trans = list(e_start_trans_dict[e_start].keys())
            e_end_trans = list(e_end_trans_dict[e_end].keys())
            
            e_start_trans.sort()
            e_end_trans.sort()
            
            e_start_trans_line = ",".join(e_start_trans)
            e_end_trans_line = ",".join(e_end_trans)
            
            self.e_start_trans_list.append(e_start_trans_line)
            self.e_end_trans_list.append(e_end_trans_line)
            
            
    
    def format_bed_line(self):
        bed_list = []
        bed_list.append(str(self.scaffold))
        bed_list.append(str(self.start_pos))
        bed_list.append(str(self.end_pos))
        
        gene_id = self.trans_id.split(".")[0]
        id_line = ";".join([gene_id,self.trans_id])
        
        bed_list.append(str(id_line))
        bed_list.append("40")
        bed_list.append(self.strand)
        
        bed_list.append(str(self.start_pos))
        bed_list.append(str(self.end_pos))
        
        
        this_source_dict = {} # this_source_dict[source] = 1
        for source_trans_id in self.merged_trans_dict:
            source_trans_obj = self.merged_trans_dict[source_trans_id]
            this_source = source_trans_obj.source_id
            this_source_dict[this_source] = 1
        
        this_source_list = list(this_source_dict.keys())
        if len(this_source_list) > 1:
            rgb_colour_code = rgb_colour_dict[10] # magenta for mixed sources
        elif len(this_source_list) == 1:
            colour_num = source_colour_dict[this_source_list[0]]
            rgb_colour_code = rgb_colour_dict[colour_num]

        #bed_list.append("255,0,0")
        bed_list.append(rgb_colour_code)
        
        bed_list.append(str(self.num_exons))
        
        relative_exon_start_list = []
        exon_length_list = []
        for i in xrange(self.num_exons):
            exon_start = self.collapse_start_list[i]
            exon_end = self.collapse_end_list[i]
            exon_length = exon_end - exon_start
            relative_exon_start = exon_start - self.start_pos
            
            relative_exon_start_list.append(str(relative_exon_start))
            exon_length_list.append(str(exon_length))
            
        relative_exon_start_line = ",".join(relative_exon_start_list)
        exon_length_line = ",".join(exon_length_list)
                
        bed_list.append(exon_length_line)
        bed_list.append(relative_exon_start_line)
        
        bed_line = "\t".join(bed_list)
        
        return bed_line
    
    def format_trans_report_line(self):
        trans_report_list = []
        trans_report_list.append(self.trans_id)
        
        trans_report_list.append(str(self.num_trans))
        
        #get all sources
        merged_source_dict = {} #merged_source_dict[source id] = 1
        for merged_trans_id in self.merged_trans_dict:
            a_trans_obj = self.merged_trans_dict[merged_trans_id]
            source_id = a_trans_obj.source_id
            merged_source_dict[source_id] = 1
        
        merged_source_list = list(merged_source_dict.keys())
        merged_source_list.sort()
        merged_source_line = ",".join(merged_source_list)
        trans_report_list.append(merged_source_line)
        
        start_wobble_string_list = []
        end_wobble_string_list = []
        for i in xrange(len(self.start_wobble_list)):
            start_wobble_string = str(self.start_wobble_list[i])
            start_wobble_string_list.append(start_wobble_string)
            end_wobble_string = str(self.end_wobble_list[i])
            end_wobble_string_list.append(end_wobble_string)
            
        start_wobble_line = ",".join(start_wobble_string_list)
        end_wobble_line = ",".join(end_wobble_string_list)
        
        trans_report_list.append(start_wobble_line)
        trans_report_list.append(end_wobble_line)
        
        e_start_trans_support = ";".join(self.e_start_trans_list)
        e_end_trans_support = ";".join(self.e_end_trans_list)
        
        trans_report_list.append(e_start_trans_support)
        trans_report_list.append(e_end_trans_support)
       
        trans_report_line = "\t".join(trans_report_list)
        
        return trans_report_line

#above this line are def's used for looping through sam file
####################################################################################################
####################################################################################################
####################################################################################################
#below this line are def's use for post sam pocessing



####################################################################################################
####################################################################################################

def fuzzy_match(coord1, coord2, diff_threshold):  # use this to allow for fuzzy matches of splice junctions
    diff_num = 0
    match_flag = "none"
    # diff_threshold = 10

    if coord1 == coord2:
        match_flag = "perfect_match"
    else:
        diff_num = coord1 - coord2

    if match_flag == "none":
        if abs(diff_num) <= diff_threshold:
            match_flag = "wobbly_match"
        else:
            match_flag = "no_match"

    if match_flag == "none":
        print("Error with match flag")
        sys.exit()

    return match_flag, diff_num
####################################################################################################

def compare_transcripts_both_capped(trans_obj,o_trans_obj,fivecap_flag,o_fivecap_flag ,strand): #use this to compare two transcripts
    
    #check cap flag, both transcripts should be capped
    if fivecap_flag != "capped" or o_fivecap_flag != "capped":
        print("Error, both transcripts need to be from capped libraries!")
        print(fivecap_flag + " " + o_fivecap_flag)
        sys.exit()
    
    diff_num_exon_flag = 0
    max_exon_num = 0
    min_exon_num = 0
    
    e_start_list = trans_obj.exon_start_list
    o_e_start_list = o_trans_obj.exon_start_list
    
    e_end_list = trans_obj.exon_end_list
    o_e_end_list = o_trans_obj.exon_end_list
    
    if len(e_start_list) != len(o_e_start_list):
        diff_num_exon_flag = 1
    
    if diff_num_exon_flag == 1: # if 5prime capped then should have same number of exons
        trans_comp_flag = "diff_transcripts"
        start_match_list = []
        start_diff_list = []
        end_match_list = []
        end_diff_list = []
        short_trans = "none"
        min_exon_num = 0
    else:
        short_trans = "same"        
        min_exon_num = len(e_start_list)
        
        start_match_list = []
        start_diff_list = []
        end_match_list = []
        end_diff_list = []
        
        all_match_flag = 1 # 1 if all matching and 0 if at least one not matching
        
        for i in xrange(min_exon_num): # iterate from 3' end of transcript, strand corrected
            
            if strand == "+":
                j = -1 * (i + 1) #iterate from last exon to account for possible 5' degradation for forward strand
            elif strand == "-":
                j = i # iterate from first exon for reverse strand
            
            
            e_start = e_start_list[j]
            o_e_start = o_e_start_list[j]
            e_end = e_end_list[j]
            o_e_end = o_e_end_list[j]
            
            # set as default before deciding which thresholds to use depending on situation
            # Need to know which exon this if for threshold setting
            start_threshold = exon_diff_threshold
            end_threshold = exon_diff_threshold
            
            if strand == "+":
                if i == 0: #use three prime threshold if this is last exon
                    end_threshold = threeprime_threshold
                if diff_num_exon_flag == 0 and i == min_exon_num-1: #use 5 prime threshold if this is first exon
                    start_threshold = fiveprime_threshold

            elif strand == "-":
                if i == 0: #use three prime threshold if this is last exon
                    start_threshold = threeprime_threshold
                if diff_num_exon_flag == 0 and i == min_exon_num-1: #use 5 prime threshold if this is first exon
                    end_threshold = fiveprime_threshold
            
            start_match_flag,start_diff_num = fuzzy_match(e_start,o_e_start,start_threshold)
            end_match_flag,end_diff_num = fuzzy_match(e_end,o_e_end,end_threshold)
            
            #remember this is for both 5cap
            if start_match_flag == "no_match" or end_match_flag == "no_match":
                all_match_flag = 0
            
            
            start_match_list.append(start_match_flag)
            start_diff_list.append(start_diff_num)
            end_match_list.append(end_match_flag)
            end_diff_list.append(end_diff_num)
            
        
        trans_comp_flag = "none"
        
        if all_match_flag == 1:        
            trans_comp_flag = "same_transcript"
        else:
            trans_comp_flag = "diff_transcripts"
        
        #Keep in mind that the lists are ordered from 3' end to 5' end
        
    return trans_comp_flag,start_match_list,start_diff_list,end_match_list,end_diff_list,short_trans,min_exon_num

####################################################################################################

def compare_transcripts_capped_nocap(trans_obj,o_trans_obj,fivecap_flag,o_fivecap_flag ,strand): #use this to compare two transcripts
    
    #check cap situation
    if fivecap_flag == "capped" and o_fivecap_flag == "no_cap":
        c_trans_obj = trans_obj
        n_trans_obj = o_trans_obj
        #should be implied but just for checking
        c_fivecap_flag = fivecap_flag
        n_fivecap_flag = o_fivecap_flag
    elif fivecap_flag == "no_cap" and o_fivecap_flag == "capped":
        c_trans_obj = o_trans_obj
        n_trans_obj = trans_obj
        #should be implied but just for checking
        c_fivecap_flag = o_fivecap_flag
        n_fivecap_flag = fivecap_flag
    else:
        print("Error with cap flags, one needs to be capped and the other no cap")
        print(fivecap_flag +" "+o_fivecap_flag)
        sys.exit()
    
    diff_num_exon_flag = 0
    max_exon_num = 0
    min_exon_num = 0
    
    c_e_start_list = c_trans_obj.exon_start_list
    n_e_start_list = n_trans_obj.exon_start_list
    
    c_e_end_list = c_trans_obj.exon_end_list
    n_e_end_list = n_trans_obj.exon_end_list
    
    if len(c_e_start_list) != len(n_e_start_list):
        diff_num_exon_flag = 1

    trans_comp_flag = "unassigned"
    c_num_exons = len(c_e_start_list)
    n_num_exons = len(n_e_start_list)

    if diff_num_exon_flag == 1:  # nocap should not have more exons than capped

        if n_num_exons > c_num_exons:
            trans_comp_flag = "diff_transcripts"
            start_match_list = []
            start_diff_list = []
            end_match_list = []
            end_diff_list = []
            short_trans = "none"
            min_exon_num = 0

    if trans_comp_flag != "diff_transcripts":

        max_exon_num = c_num_exons
        min_exon_num = n_num_exons

        short_trans = n_trans_obj.trans_id


        start_match_list = []
        start_diff_list = []
        end_match_list = []
        end_diff_list = []

        all_match_flag = 1  # 1 if all matching and 0 if at least one not matching

        for i in xrange(min_exon_num):  # iterate from 3' end of transcript, strand corrected

            if strand == "+":
                j = -1 * (i + 1)  # iterate from last exon to account for possible 5' degradation for forward strand
            elif strand == "-":
                j = i  # iterate from first exon for reverse strand

            c_e_start = c_e_start_list[j]
            n_e_start = n_e_start_list[j]
            c_e_end = c_e_end_list[j]
            n_e_end = n_e_end_list[j]

            start_threshold = exon_diff_threshold
            end_threshold = exon_diff_threshold

            if strand == "+":
                if i == 0:  # use three prime threshold if this is last exon
                    end_threshold = threeprime_threshold
                if i == max_exon_num - 1:  # use 5 prime threshold if this is first exon, never reaches max_exon_num if diff exon num
                    start_threshold = fiveprime_threshold

            elif strand == "-":
                if i == 0:  # use three prime threshold if this is last exon
                    start_threshold = threeprime_threshold
                if i == max_exon_num - 1:  # use 5 prime threshold if this is first exon, never reaches max_exon_num if diff exon num
                    end_threshold = fiveprime_threshold

            start_match_flag, start_diff_num = fuzzy_match(c_e_start, n_e_start, start_threshold)
            end_match_flag, end_diff_num = fuzzy_match(c_e_end, n_e_end, end_threshold)

            if i < min_exon_num - 1 : # not final exon
                if end_match_flag == "no_match" or start_match_flag == "no_match":
                    all_match_flag = 0
                    trans_comp_flag = "diff_transcripts"
                    break

            # Note that these occur at last loop, so no need ot break after match decision
            elif max_exon_num > min_exon_num and i == min_exon_num - 1: # diff num exons and at first min exon
                if strand == "+":
                    if end_match_flag == "no_match":
                        all_match_flag = 0
                        trans_comp_flag = "diff_transcripts"
                    elif start_match_flag == "no_match":
                        if c_e_start > n_e_start: #if capped start is later than nocap start then no match
                            all_match_flag = 0
                            trans_comp_flag = "diff_transcripts"
                        else: # same three prime exons, with first no cap exon not as long as matching capped exon
                            trans_comp_flag = "same_three_prime_diff_exons"
                    else: #  matches all exons for nocap but diff num exons
                        trans_comp_flag = "same_three_prime_diff_exons"
                elif strand == "-":
                    if start_match_flag == "no_match":
                        all_match_flag = 0
                        trans_comp_flag = "diff_transcripts"
                    elif end_match_flag == "no_match":
                        if c_e_end < n_e_end: # for reverse strand, if capped start is not as long as nocap start then no match
                            all_match_flag = 0
                            trans_comp_flag = "diff_transcripts"
                        else: # same three prime exons, with first no cap exon not as long as matching capped exon
                            trans_comp_flag = "same_three_prime_diff_exons"
                    else:  # matches all exons for nocap but diff num exons
                        trans_comp_flag = "same_three_prime_diff_exons"

            elif max_exon_num == min_exon_num and i == min_exon_num - 1: # same number of exons and first exon
                if strand == "+":
                    if end_match_flag == "no_match":
                        all_match_flag = 0
                        trans_comp_flag = "diff_transcripts"
                    elif start_match_flag == "no_match":
                        if c_e_start > n_e_start: #if capped start is later than nocap start then no match
                            all_match_flag = 0
                            trans_comp_flag = "diff_transcripts"
                        else: # start does not match but start exon matches
                            trans_comp_flag = "same_three_prime_same_exons"
                    else: # complete match
                        trans_comp_flag = "same_transcript"
                elif strand == "-":
                    if start_match_flag == "no_match":
                        all_match_flag = 0
                        trans_comp_flag = "diff_transcripts"
                    elif end_match_flag == "no_match":
                        if c_e_end < n_e_end: # for reverse strand, if capped start is not as long as nocap start then no match
                            all_match_flag = 0
                            trans_comp_flag = "diff_transcripts"
                        else: # start does not match but start exon matches
                            trans_comp_flag = "same_three_prime_same_exons"
                    else:  # complete match
                        trans_comp_flag = "same_transcript"


            start_match_list.append(start_match_flag)
            start_diff_list.append(start_diff_num)
            end_match_list.append(end_match_flag)
            end_diff_list.append(end_diff_num)

    # note that comparison stops at first no match signal
    # so the match lists may not be as long as number of exons if there is an internal no match
    # Keep in mind that the lists are ordered from 3' end to 5' end
    return trans_comp_flag, start_match_list, start_diff_list, end_match_list, end_diff_list, short_trans, min_exon_num

####################################################################################################
####################################################################################################

# use this for compare_transcripts_both_nocap
def assign_short_long_trans(s_trans_obj,l_trans_obj):
    # short trans obj is first in arguments
    l_e_start_list = l_trans_obj.exon_start_list # long
    s_e_start_list = s_trans_obj.exon_start_list # short
            
    l_e_end_list = l_trans_obj.exon_end_list #long
    s_e_end_list = s_trans_obj.exon_end_list #short
    
    max_exon_num = len(l_e_start_list)
    min_exon_num = len(s_e_start_list)
    
    return(s_e_start_list,s_e_end_list,l_e_start_list,l_e_end_list,max_exon_num,min_exon_num)

####################################################################################################

def compare_transcripts_both_nocap(trans_obj,o_trans_obj,fivecap_flag,o_fivecap_flag ,strand): #use this to compare two transcripts
    diff_num_exon_flag = 0
    max_exon_num = 0
    min_exon_num = 0

    #check cap situation
    if fivecap_flag != "no_cap" or o_fivecap_flag != "no_cap":
        print("At least one transcript is capped, both need to be nocap")
        sys.exit()
    
    e_start_list = trans_obj.exon_start_list
    o_e_start_list = o_trans_obj.exon_start_list
    
    e_end_list = trans_obj.exon_end_list
    o_e_end_list = o_trans_obj.exon_end_list
    
    ############################
    #assign short and long trans
    if len(e_start_list) != len(o_e_start_list):
        diff_num_exon_flag = 1
    
        if len(e_start_list) > len(o_e_start_list):
            short_trans = o_trans_obj.trans_id
            
            # short trans is first in argument
            [s_e_start_list,s_e_end_list,l_e_start_list,l_e_end_list,max_exon_num,min_exon_num] = assign_short_long_trans(o_trans_obj,trans_obj)
            
        elif len(e_start_list) < len(o_e_start_list):
            short_trans = trans_obj.trans_id
            # short trans is first in argument
            [s_e_start_list,s_e_end_list,l_e_start_list,l_e_end_list,max_exon_num,min_exon_num] = assign_short_long_trans(trans_obj,o_trans_obj)
            
    elif len(e_start_list) == len(o_e_start_list):
        if strand == "+":
            # assigan first as long trans if longer or equal
            if e_start_list[0] <= o_e_start_list[0]:
                short_trans = o_trans_obj.trans_id
                # short trans is first in argument
                [s_e_start_list,s_e_end_list,l_e_start_list,l_e_end_list,max_exon_num,min_exon_num] = assign_short_long_trans(o_trans_obj,trans_obj)
            elif e_start_list[0] > o_e_start_list[0]:
                short_trans = trans_obj.trans_id
                # short trans is first in argument
                [s_e_start_list,s_e_end_list,l_e_start_list,l_e_end_list,max_exon_num,min_exon_num] = assign_short_long_trans(trans_obj,o_trans_obj)
            else:
                print("Error with short long trans assigning")
                sys.exit()
        
        elif strand == "-":
            # assigan first as long trans if longer or equal
            if e_end_list[-1] >= o_e_end_list[-1]:
                short_trans = o_trans_obj.trans_id
                # short trans is first in argument
                [s_e_start_list,s_e_end_list,l_e_start_list,l_e_end_list,max_exon_num,min_exon_num] = assign_short_long_trans(o_trans_obj,trans_obj)
            elif e_end_list[-1] < o_e_end_list[-1]:
                short_trans = trans_obj.trans_id
                # short trans is first in argument
                [s_e_start_list,s_e_end_list,l_e_start_list,l_e_end_list,max_exon_num,min_exon_num] = assign_short_long_trans(trans_obj,o_trans_obj)
            else:
                print("Error with short long trans assigning")
                sys.exit()
    else:
        print("Error with short long trans assigning")
        sys.exit()
    #assign short and long trans above
    ############################
    ############################
    
    start_match_list = []
    start_diff_list = []
    end_match_list = []
    end_diff_list = []

    all_match_flag = 1  # 1 if all matching and 0 if at least one not matching

    for i in xrange(min_exon_num):  # iterate from 3' end of transcript, strand corrected

        if strand == "+":
            j = -1 * (i + 1)  # iterate from last exon to account for possible 5' degradation for forward strand
        elif strand == "-":
            j = i  # iterate from first exon for reverse strand

        l_e_start = l_e_start_list[j]
        s_e_start = s_e_start_list[j]
        l_e_end = l_e_end_list[j]
        s_e_end = s_e_end_list[j]

        start_threshold = exon_diff_threshold
        end_threshold = exon_diff_threshold

        if strand == "+":
            if i == 0:  # use three prime threshold if this is last exon
                end_threshold = threeprime_threshold
            if i == max_exon_num - 1:  # use 5 prime threshold if this is first exon, never reaches max_exon_num if diff exon num
                start_threshold = fiveprime_threshold

        elif strand == "-":
            if i == 0:  # use three prime threshold if this is last exon
                start_threshold = threeprime_threshold
            if i == max_exon_num - 1:  # use 5 prime threshold if this is first exon, never reaches max_exon_num if diff exon num
                end_threshold = fiveprime_threshold

        start_match_flag, start_diff_num = fuzzy_match(l_e_start, s_e_start, start_threshold)
        end_match_flag, end_diff_num = fuzzy_match(l_e_end, s_e_end, end_threshold)
    #################################################################

        if i < min_exon_num - 1 : # not final exon
            if end_match_flag == "no_match" or start_match_flag == "no_match":
                all_match_flag = 0
                trans_comp_flag = "diff_transcripts"
                break

        # Note that these occur at last loop, so no need ot break after match decision
        elif max_exon_num > min_exon_num and i == min_exon_num - 1: # diff num exons and at first min exon, last loop
            if strand == "+":
                if end_match_flag == "no_match":
                    all_match_flag = 0
                    trans_comp_flag = "diff_transcripts"
                elif start_match_flag == "no_match":
                    if l_e_start > s_e_start: #if capped start is later than nocap start then no match
                        all_match_flag = 0
                        trans_comp_flag = "diff_transcripts"
                    else: # same three prime exons, with first no cap exon not as long as matching capped exon
                        trans_comp_flag = "same_three_prime_diff_exons"
                else: #  matches all exons for nocap but diff num exons
                    trans_comp_flag = "same_three_prime_diff_exons"
            elif strand == "-":
                if start_match_flag == "no_match":
                    all_match_flag = 0
                    trans_comp_flag = "diff_transcripts"
                elif end_match_flag == "no_match":
                    if l_e_end < s_e_end: # for reverse strand, if capped start is not as long as nocap start then no match
                        all_match_flag = 0
                        trans_comp_flag = "diff_transcripts"
                    else: # same three prime exons, with first no cap exon not as long as matching capped exon
                        trans_comp_flag = "same_three_prime_diff_exons"
                else:  # matches all exons for nocap but diff num exons
                    trans_comp_flag = "same_three_prime_diff_exons"

        elif max_exon_num == min_exon_num and i == min_exon_num - 1: # same number of exons and first exon, last loop
            if strand == "+":
                if end_match_flag == "no_match":
                    all_match_flag = 0
                    trans_comp_flag = "diff_transcripts"
                elif start_match_flag == "no_match":
                    if l_e_start > s_e_start: #if capped start is later than nocap start then no match
                        all_match_flag = 0
                        trans_comp_flag = "diff_transcripts"
                    else: # start does not match but start exon matches
                        trans_comp_flag = "same_three_prime_same_exons"
                else: # complete match
                    trans_comp_flag = "same_transcript"
            elif strand == "-":
                if start_match_flag == "no_match":
                    all_match_flag = 0
                    trans_comp_flag = "diff_transcripts"
                elif end_match_flag == "no_match":
                    if l_e_end < s_e_end: # for reverse strand, if capped start is not as long as nocap start then no match
                        all_match_flag = 0
                        trans_comp_flag = "diff_transcripts"
                    else: # start does not match but start exon matches
                        trans_comp_flag = "same_three_prime_same_exons"
                else:  # complete match
                    trans_comp_flag = "same_transcript"


        start_match_list.append(start_match_flag)
        start_diff_list.append(start_diff_num)
        end_match_list.append(end_match_flag)
        end_diff_list.append(end_diff_num)

    # note that comparison stops at first no match signal
    # so the match lists may not be as long as number of exons if there is an internal no match
    # Keep in mind that the lists are ordered from 3' end to 5' end
    return trans_comp_flag, start_match_list, start_diff_list, end_match_list, end_diff_list, short_trans, min_exon_num

    ####################################################################################################
    #################################################################################################### Work here!! 2017/06/30 !!!
    ####################################################################################################
####################################################################################################

####################################################################################################

def collapse_transcripts(trans_obj_list,collapse_flag): #use this to collapse transcripts
    # all supplied transcripts will be merged
    # create supplied transcripts by comparing with compare transcripts
    
    try:
        collapse_flag
    except NameError:
        print "collapse_flag not defined, using default of most commond ends"
        collapse_flag == "common_ends"
    
    max_exon_num = 0
    strand = "none"
    num_trans = len(trans_obj_list)
    
    # assume they would never use a priority ranking up to 1000
    best_start_priority = 1000
    best_end_priority = 1000
    best_junction_priority = 1000
    
    #check strand and get max exon num
    for trans_obj in trans_obj_list:
        e_start_list = trans_obj.exon_start_list
        if strand == "none":
            strand = trans_obj.strand
        elif trans_obj.strand != strand:
            print("mismatch in strand from trans_obj_list for collapsing def")
            sys.exit()
        exon_num = len(e_start_list)
        if exon_num > max_exon_num:
            max_exon_num = exon_num
        
        # figure out best priority for different features
        start_priority = trans_obj.start_priority
        end_priority = trans_obj.end_priority 
        junction_priority = trans_obj.junction_priority
        
        if best_start_priority > start_priority:
            best_start_priority = start_priority
        if best_end_priority > end_priority:
            best_end_priority = end_priority
        if best_junction_priority > junction_priority:
            best_junction_priority = junction_priority
        
        
        
    
    collapse_start_list = []
    collapse_end_list = []
    e_start_trans_dict = {} # start_trans_dict[e start][uniq_trans_id] = 1
    e_end_trans_dict = {} # e_end_trans_dict[e end][uniq_trans_id] = 1
    
    #track how much wobble for the starts and end in the collapse
    start_wobble_list = []
    end_wobble_list = []
    for i in xrange(max_exon_num): #go from 3 prime end
        if strand == "+":
            j = -1 * (i + 1) #iterate from last exon to account for possible 5' degradation for forward strand
        elif strand == "-":
            j = i # iterate from first exon for reverse strand
        
        e_start_dict = {} # e_start_dict[start] = number of occurrences
        e_end_dict = {} # e_end_dict[end] = number of occurrences
        for trans_obj in trans_obj_list:
            e_start_list = trans_obj.exon_start_list
            e_end_list = trans_obj.exon_end_list
            
            # figure out priority for different features
            start_priority = trans_obj.start_priority
            end_priority = trans_obj.end_priority 
            junction_priority = trans_obj.junction_priority 
            
            if i >= len(e_start_list):# use for no cap when exon numbers may not match
                continue
            
            e_start = int(e_start_list[j])
            e_end = int(e_end_list[j])
            
            #do not use 5' end if this is not priority level one for start
            if i == len(e_start_list)-1 and start_priority != best_start_priority:
                if strand == "+":
                    e_start = -1
                elif strand == "-":
                    e_end = -1

            #do not use splice junctions if not priority
            if junction_priority != best_junction_priority and max_exon_num > 1:
                if i > 0 and i < len(e_start_list)-1 :
                    continue 
                if i == 0:
                    if strand == "+":
                        e_end = -1
                    elif strand == "-":
                        e_start = -1
                    
                if i == len(e_start_list)-1:
                    if strand == "+":
                        e_start = -1
                    elif strand == "-":
                        e_end = -1
            
            # do not use 3' end if not priority
            if i == 0 and end_priority != best_end_priority:
                if strand == "+":
                    e_end = -1
                elif strand == "-":
                    e_start = -1
                
                
            
            if e_start != -1:
                if e_start not in e_start_dict:
                    e_start_dict[e_start] = 0
                    e_start_trans_dict[e_start] = {}
                e_start_dict[e_start] += 1
                e_start_trans_dict[e_start][trans_obj.uniq_trans_id] = 1
            
            if e_end != -1:
                if e_end not in e_end_dict:
                    e_end_dict[e_end] = 0
                    e_end_trans_dict[e_end] = {}
                e_end_dict[e_end] += 1
                e_end_trans_dict[e_end][trans_obj.uniq_trans_id] = 1
        
        
        ##########################################
        best_e_start = -1
        long_e_start = -1
        short_e_start = -1
        num_starts = 0
        for e_start in e_start_dict:
            if e_start_dict[e_start] > num_starts:
                best_e_start = e_start
                num_starts = e_start_dict[e_start]
            
            if long_e_start == -1:
                long_e_start = e_start
            if e_start < long_e_start:
                long_e_start = e_start
            
            if short_e_start == -1:
                short_e_start = e_start
            if e_start > short_e_start:
                short_e_start = e_start
        
        # if there are multiple most num e starts then choose the longest one
        most_starts = num_starts
        num_most_starts = 0
        most_long_e_start = -1
        for e_start in e_start_dict:
            if e_start_dict[e_start] == most_starts:
                num_most_starts += 1
                if most_long_e_start == -1:
                    most_long_e_start = e_start
                elif most_long_e_start > e_start:
                    most_long_e_start = e_start
        if num_most_starts > 1:
            best_e_start = most_long_e_start
            if num_trans > 2:
                print("more than one best e start! " + str(best_e_start) + " num_trans: " + str(num_trans))
        ##########################################
        
        e_start_wobble = short_e_start - long_e_start
        start_wobble_list.append(e_start_wobble)
        
        best_e_end = 0
        long_e_end = -1
        short_e_end = -1
        num_ends = 0
        for e_end in e_end_dict:
            if e_end_dict[e_end] > num_ends:
                best_e_end = e_end
                num_ends = e_end_dict[e_end]
            
            if long_e_end == -1:
                long_e_end = e_end
            if e_end > long_e_end:
                long_e_end = e_end
            
            if short_e_end == -1:
                short_e_end = e_end
            if e_end < short_e_end:
                short_e_end = e_end
        
        # if there are multiple most num e ends then choose the longest one
        most_ends = num_ends
        num_most_ends = 0
        most_long_e_end = -1
        for e_end in e_end_dict:
            if e_end_dict[e_end] == most_ends:
                num_most_ends += 1
                if most_long_e_end == -1:
                    most_long_e_end = e_end
                elif most_long_e_end < e_end:
                    most_long_e_end = e_end
        if num_most_ends > 1:
            best_e_end = most_long_e_end
            if num_trans > 2:
                print("more than one best e end! " + str(best_e_end) + " num_trans: " + str(num_trans))
        ##########################################
        
        e_end_wobble = long_e_end - short_e_end
        end_wobble_list.append(e_end_wobble)
        
        if collapse_flag == "longest_ends": #allow user to choose whether to use the most common ends or the longest ends, default is most common ends
            if i+1 == max_exon_num:
                if strand == "+":
                    best_e_start = long_e_start
                elif strand == "-":
                    best_e_end = long_e_end
                    
            if i == 0:
                if strand == "+":
                    best_e_end = long_e_end
                elif strand == "-":
                    best_e_start = long_e_start
        
        collapse_start_list.append(best_e_start)
        collapse_end_list.append(best_e_end)

    #put the coords in the right order maintaining order with wobble lists
    collapse_start_list, start_wobble_list = zip(*sorted(zip(collapse_start_list, start_wobble_list)))
    collapse_end_list, end_wobble_list = zip(*sorted(zip(collapse_end_list, end_wobble_list)))
    
    collapse_start_list = list(collapse_start_list)
    start_wobble_list = list(start_wobble_list)
    collapse_end_list = list(collapse_end_list)
    end_wobble_list = list(end_wobble_list)
    
    return collapse_start_list,collapse_end_list,start_wobble_list,end_wobble_list,e_start_trans_dict,e_end_trans_dict

####################################################################################################

def gene_group(trans_obj_list): #groups trans into genes, does not take into account strand
    

    gene_trans_dict = {} # gene_trans_dict[gene id][trans id] = 1
    trans_gene_dict = {} # trans_gene_dict[trans id] = gene group
    gene_start_dict = {} # gene_start_dict[gene num] = gene start
    start_gene_dict = {} # start_gene_dict[start] = gene num

    gene_count = 0
    
    if len(trans_obj_list) == 1:
        gene_count += 1
        uniq_trans_id = trans_obj_list[0].uniq_trans_id
        single_gene_start = trans_obj_list[0].exon_start_list[0]
        gene_start_dict[gene_count] = single_gene_start
        gene_trans_dict[gene_count] = {}
        gene_trans_dict[gene_count][uniq_trans_id] = 1
        trans_gene_dict[uniq_trans_id] = gene_count
        
    
    for i in xrange(len(trans_obj_list)):
        trans_obj = trans_obj_list[i]
        for j in xrange(i+1,len(trans_obj_list)):
            o_trans_obj = trans_obj_list[j]

            uniq_trans_id  = trans_obj.uniq_trans_id
            o_uniq_trans_id = o_trans_obj.uniq_trans_id
            
            #if uniq_trans_id == o_uniq_trans_id:#skip if same
            #    continue
            
            if uniq_trans_id in trans_gene_dict and o_uniq_trans_id in trans_gene_dict: # skip if already in same group
                if trans_gene_dict[uniq_trans_id] == trans_gene_dict[o_uniq_trans_id]:
                    continue
            
            exon_start_list = trans_obj.exon_start_list
            o_exon_start_list = o_trans_obj.exon_start_list
            
            exon_end_list = trans_obj.exon_end_list
            o_exon_end_list = o_trans_obj.exon_end_list
            
            num_exons = len(exon_start_list)
            o_num_exons = len(o_exon_start_list)
            
            overlap_flag = 0
            
            for i in xrange(num_exons): #search for overlapping exons
                for j in xrange(o_num_exons):
                    exon_start = exon_start_list[i]
                    exon_end = exon_end_list[i]
                    o_exon_start = o_exon_start_list[j]
                    o_exon_end = o_exon_end_list[j]
                    
                    if exon_start <= o_exon_end and exon_end >= o_exon_start:
                        overlap_flag = 1
            
            if overlap_flag == 0: # no overlap make new gene groups

                
                if uniq_trans_id not in trans_gene_dict: #if no gene groups make new one
                    gene_count += 1
                    trans_gene_dict[uniq_trans_id] = gene_count

                    gene_trans_dict[gene_count] = {}
                    gene_trans_dict[gene_count][uniq_trans_id] = 1

                    #add gene start
                    gene_start_dict[gene_count] = exon_start_list[0]
                
                if o_uniq_trans_id not in trans_gene_dict: #if no gene groups make new one
                    gene_count += 1
                    trans_gene_dict[o_uniq_trans_id] = gene_count

                    gene_trans_dict[gene_count] = {}
                    gene_trans_dict[gene_count][o_uniq_trans_id] = 1

                    #add gene start
                    gene_start_dict[gene_count] = o_exon_start_list[0]
                
                    
                
            if overlap_flag == 1:
                if uniq_trans_id not in trans_gene_dict and o_uniq_trans_id not in trans_gene_dict: #if no gene groups make new one
                    gene_count += 1
                    trans_gene_dict[uniq_trans_id] = gene_count
                    trans_gene_dict[o_uniq_trans_id] = gene_count
                    gene_trans_dict[gene_count] = {}
                    gene_trans_dict[gene_count][uniq_trans_id] = 1
                    gene_trans_dict[gene_count][o_uniq_trans_id] = 1
                    #add gene start
                    min_gene_start = exon_start_list[0]
                    if min_gene_start > o_exon_start_list[0]:
                        min_gene_start = o_exon_start_list[0]
                    gene_start_dict[gene_count] = min_gene_start
                    
                    
                elif uniq_trans_id not in trans_gene_dict: # add to other gene group

                    gene_num = trans_gene_dict[o_uniq_trans_id]
                    trans_gene_dict[uniq_trans_id] = gene_num
                    gene_trans_dict[gene_num][uniq_trans_id] = 1
                    
                    min_gene_start = exon_start_list[0]
                    if min_gene_start > o_exon_start_list[0]:
                        min_gene_start = o_exon_start_list[0]
                    gene_start_dict[gene_num] = min_gene_start
                elif o_uniq_trans_id not in trans_gene_dict:# add to other gene group
                    
                    gene_num = trans_gene_dict[uniq_trans_id]
                    trans_gene_dict[o_uniq_trans_id] = gene_num
                    gene_trans_dict[gene_num][o_uniq_trans_id] = 1
                    
                    min_gene_start = exon_start_list[0]
                    if min_gene_start > o_exon_start_list[0]:
                        min_gene_start = o_exon_start_list[0]
                    gene_start_dict[gene_num] = min_gene_start
                elif uniq_trans_id in trans_gene_dict and o_uniq_trans_id in trans_gene_dict:

                    gene_num = trans_gene_dict[uniq_trans_id]
                    o_gene_num = trans_gene_dict[o_uniq_trans_id]
                    
                    if gene_num != o_gene_num: #merge gene groups
                        m_uniq_trans_id_list = list(gene_trans_dict[o_gene_num].keys())
                        for m_uniq_trans_id in m_uniq_trans_id_list:
                            
                            trans_gene_dict[m_uniq_trans_id] = gene_num
                            gene_trans_dict[gene_num][m_uniq_trans_id] = 1
                        #delete old gene num
                        gene_trans_dict.pop(o_gene_num, None)
                        
                        min_gene_start = gene_start_dict[gene_num]
                        if min_gene_start > gene_start_dict[o_gene_num]:
                            min_gene_start = gene_start_dict[o_gene_num]
                        gene_start_dict[gene_num] = min_gene_start
                        #delete old gene num
                        gene_start_dict.pop(o_gene_num, None)
                        
                    if gene_num == o_gene_num: #same gene groups
                        continue
                else:
                    print("Unknown condition in gene grouping")
                    sys.exit()
                    
                    
    for gene_num in gene_start_dict: #make dict for coordinate to gene num
        gene_start = gene_start_dict[gene_num]
        if gene_start in start_gene_dict:
            print("multiple gene starts!")
            sys.exit()
        start_gene_dict[gene_start] = gene_num
    
    start_gene_list = start_gene_dict.keys()
    start_gene_list.sort()
    
    gene_start_trans_dict = {} # gene_start_trans_dict[gene start][trans id] = 1
    for gene_start in start_gene_list: # make duct for gene starts to trans
        gene_num = start_gene_dict[gene_start]
        gene_start_trans_dict[gene_start] = {}
        for uniq_trans_id in gene_trans_dict[gene_num]:
            gene_start_trans_dict[gene_start][uniq_trans_id] = 1
    
        
    return gene_start_trans_dict,start_gene_list

####################################################################################################

def sort_transcripts(trans_obj_list,trans_obj_dict):
    #sort transcripts by start-end-exon starts
    #pad with 0 for numerical sort (prevents 23 from being after 203)
    
    pos_trans_dict = {} # pos_trans_dict[pos] = trans obj
    pos_trans_list = []
    
    for trans_obj in trans_obj_list:
        trans_exon_start_list = trans_obj.exon_start_list
        trans_exon_end_list = trans_obj.exon_end_list
        trans_exon_start_list.sort()
        trans_exon_end_list.sort()
        trans_start = trans_exon_start_list[0]
        trans_end = trans_exon_end_list[-1]

        #max_digit_length = len(str(trans_end))
        max_digit_length = 10
        max_pos_line_length = 2000
        
        trans_pos_list = []
        
        trans_start_pad = str(trans_start).rjust(max_digit_length,'0') #pad with 0 on left side
        trans_end_pad = str(trans_end).rjust(max_digit_length,'0') #pad with 0 on left side
        
        
        trans_pos_list.append("1")
        trans_pos_list.append(str(trans_start_pad))
        trans_pos_list.append(str(trans_end_pad))
        
        for exon_start in trans_exon_start_list:
            exon_start_pad = str(exon_start).rjust(max_digit_length,'0') #pad with 0 on left side
            trans_pos_list.append(exon_start_pad)
        
        for exon_end in trans_exon_end_list:
            exon_end_pad = str(exon_end).rjust(max_digit_length,'0') #pad with 0 on left side
            trans_pos_list.append(exon_end_pad)
        
        trans_pos_line = "".join(trans_pos_list)
        
        if len(trans_pos_line) > max_pos_line_length:
            print("padding for trans pos line is insufficient. Use larger max_pos_line_length")
            sys.exit()
            
        trans_pos_line = trans_pos_line.ljust(max_pos_line_length,'0')
        
        trans_pos_line = int(trans_pos_line)

        
        if trans_pos_line in pos_trans_dict:
            old_merge_id = pos_trans_dict[trans_pos_line].trans_id
            new_merge_id = trans_obj.trans_id
            
            old_trans_list = pos_trans_dict[trans_pos_line].trans_list
            new_trans_list = trans_obj.trans_list            

            print("Duplicate transcript positions in transcript sorting!")
            print(trans_obj.merged_trans_dict.keys())
            print(str(trans_start)+" "+str(trans_end))
            print(pos_trans_dict[trans_pos_line].merged_trans_dict.keys())
            this_bed_line = trans_obj.format_bed_line()
            other_bed_line = pos_trans_dict[trans_pos_line].format_bed_line()
            print(this_bed_line)
            print(other_bed_line)
            print("a###########################################")
            for a_uniq_trans_id in trans_obj.merged_trans_dict:
                a_bed_line = trans_obj_dict[a_uniq_trans_id].format_bed_line(a_uniq_trans_id)
                print(a_bed_line)
            print("b###########################################")
            for b_uniq_trans_id in pos_trans_dict[trans_pos_line].merged_trans_dict:
                b_bed_line = trans_obj_dict[b_uniq_trans_id].format_bed_line(b_uniq_trans_id)
                print(b_bed_line)
            
            

            sys.exit()
        
        pos_trans_dict[trans_pos_line] = trans_obj
        pos_trans_list.append(trans_pos_line)
    
    pos_trans_list.sort()
    
    
    sorted_trans_obj_list = []
    for pos_trans in pos_trans_list:
        
        trans_obj = pos_trans_dict[pos_trans]
        sorted_trans_obj_list.append(trans_obj)
        tmp_id = trans_obj.trans_id

    
    # sorted_trans_obj_list list of trana obj that have been sorted by position
    return sorted_trans_obj_list

##############################################################################


def longest_transcript(trans_id,trans_id_list,trans_obj_dict):
    #determines if transcript is one of the longest in group on the 5 prime end
    # if longest_trans == "longest" than this is one of the longest trans
    # if equal to "short" it is not
    
    trans_obj = trans_obj_dict[trans_id]
    strand = trans_obj.strand

    trans_start = trans_obj.start_pos
    trans_end = trans_obj.end_pos
    num_exons = trans_obj.num_exons
    
    longest_trans = "none"
    
    #exon_diff_threshold = 10
    #fiveprime_threshold = 10
    #threeprime_threshold = 10
    longest_comp_list = []
    
    for o_trans_id in trans_id_list:
        o_trans_obj = trans_obj_dict[o_trans_id]

        o_trans_start = o_trans_obj.start_pos
        o_trans_end = o_trans_obj.end_pos
        o_num_exons = o_trans_obj.num_exons
        
        if trans_id == o_trans_id:
            longest_trans = "longest"
            continue
        
        if num_exons > o_num_exons:
            longest_trans = "longest"
            longest_comp_list.append(longest_trans)
        elif num_exons < o_num_exons: # this is not the longest trans in group
            longest_trans = "short"
            break
        elif num_exons == o_num_exons:
            if strand == "+":
                # only allow exon diff threshold in case of nocap super long 5 prime threshold
                start_match_flag,start_diff_num = fuzzy_match(trans_start,o_trans_start,exon_diff_threshold)
                if start_diff_num <= 0 :
                    longest_trans = "longest"
                elif start_match_flag != "no_match":
                    longest_trans = "longest"
                elif start_match_flag == "no_match":
                    longest_trans = "short"
                    break
                else:
                    print("Error with finding longest trans")
            elif strand == "-":
                # only allow exon diff threshold in case of nocap super long 5 prime threshold
                end_match_flag,end_diff_num = fuzzy_match(trans_end,o_trans_end,exon_diff_threshold)
                if end_diff_num >= 0 :
                    longest_trans = "longest"
                elif end_match_flag != "no_match":
                    longest_trans = "longest"
                elif end_match_flag == "no_match":
                    longest_trans = "short"
                    break
                else:
                    print("Error with finding longest trans")
    
    return longest_trans    
    
    
##############################################################################

####################################################################################################

def no_cap_check(this_trans_dict,trans_obj_dict): # this checks to see if all transcripts in a group are from no_cap
    
    no_cap_flag = "all_no_cap"
    for uniq_trans_id in this_trans_dict:
        trans_obj = trans_obj_dict[uniq_trans_id]
        fivecap_flag = source_dict[trans_obj.source_id]['seq_type']
        if fivecap_flag != "no_cap":
            no_cap_flag = "cap_found"
    
    return no_cap_flag


    
##############################################################################
# Use this class to manage group merging
class TransGroup:
    def __init__(self, group_name):
        self.group_name = group_name 
        self.trans_group_dict = {}  #trans_group_dict[trans id][trans group] = "longest" or "short" for no cap
        self.group_trans_dict = {}  # group_trans_dict[group num][trans id] = "longest" or "short" for no cap
        self.group_longest_dict = {} #group_longest_dict[group num][longest/long][trans] = 1
        
        self.group_count = 0
    
    def check_trans_status(self,trans_a):
        # may add on for future features
        group_check = 0

        if trans_a in self.trans_group_dict:
            group_check = 1
        
        return group_check
    
    def check_same_group(self,trans_a,trans_b):
        group_match_flag = 0
                    
        for a_trans_group in self.trans_group_dict[trans_a]:
            if a_trans_group in self.trans_group_dict[trans_b]:
                group_match_flag = 1
        
        return group_match_flag
    
    def check_nocap_group(self,trans_a):
        assoc_trans_dict = {} # assoc_trans_dict[trans id] = 1
        
        for a_trans_group in self.trans_group_dict[trans_a]:
            for trans_b in self.group_trans_dict[a_trans_group]:
                if trans_b != trans_a:
                    assoc_trans_dict[trans_b] = 1
        
        num_assoc_trans = len(list(assoc_trans_dict.keys()))
        
        return num_assoc_trans

#####################################################################
    
    def delete_trans(self,trans_a):
        
        #delete trans from all its associated groups
        for a_trans_group in self.trans_group_dict[trans_a]:
            self.group_trans_dict[a_trans_group].pop(trans_a,None)
            
            if len(list(self.group_trans_dict[a_trans_group])) == 0:# if ther group is now empty, delete it
                self.group_trans_dict.pop(a_trans_group,None)
        
        #delete trans from trans group dict
        self.trans_group_dict.pop(trans_a,None)

#####################################################################

    def new_group_a(self,trans_a):
        self.group_count += 1
        if trans_a in self.trans_group_dict:
            print("Error in new group, trans_a already in group")
            sys.exit()
        
        if self.group_count in self.group_trans_dict:
            print("group num already used")
            sys.exit()
                
        self.trans_group_dict[trans_a] = {}
        self.trans_group_dict[trans_a][self.group_count] = "longest"
        
        self.group_trans_dict[self.group_count] = {}
        self.group_trans_dict[self.group_count][trans_a] = "longest"
        

#####################################################################

   
    def add_a_to_b_group(self,trans_a,trans_b,trans_obj_dict):
        # only cap libs should be used for b group and they only have one group
        # does not take all of a_group just uses a_trans
        if len(list(self.trans_group_dict[trans_b].keys())) > 1:
            print("multiple groups")
            sys.exit()
            
        b_group_num = list(self.trans_group_dict[trans_b].keys())[0]
        
        #remove initial nocap group that is a self identity group
        if len(list(self.trans_group_dict[trans_a].keys())) == 1:
            a_trans_group = list(self.trans_group_dict[trans_a].keys())[0]
            if len(list(self.group_trans_dict[a_trans_group].keys())) == 1:
                if list(self.group_trans_dict[a_trans_group].keys())[0] == trans_a:
                    self.group_trans_dict.pop(a_trans_group,None)
                    self.trans_group_dict.pop(trans_a,None)
                else:
                    print("Error with group trans in add a to b group")
                    sys.exit()

        if trans_a not in self.trans_group_dict:
            self.trans_group_dict[trans_a] = {}
        else:
            # this happens if trans_a is a nocap trans in which case it can be in multiple groups
            print("trans_a already in group, should be nocap trans: " + trans_a)
            
        self.trans_group_dict[trans_a][b_group_num] = 1
        self.group_trans_dict[b_group_num][trans_a] = 1
        
        b_trans_id_list = list(self.group_trans_dict[b_group_num].keys())
        
        #redo longest trans flags
        for trans_c in self.group_trans_dict[b_group_num]:
            longest_trans_flag = longest_transcript(trans_c,b_trans_id_list,trans_obj_dict)
            self.trans_group_dict[trans_c][b_group_num] = longest_trans_flag
            self.group_trans_dict[b_group_num][trans_c] = longest_trans_flag

##############################################################################

    def add_a_to_b_group_both_nocap(self,trans_a,trans_b,trans_obj_dict):
        print("invoke add_a_to_b_group_both_nocap " + trans_a + " " + trans_b)
        # does not take all of a_group just uses a_trans
        if len(list(self.trans_group_dict[trans_b].keys())) > 1:
            print("multiple groups")
            #sys.exit()

        #check that trans b has more exons than trans a
        trans_obj_a = trans_obj_dict[trans_a]
        trans_obj_b = trans_obj_dict[trans_b]
        
        if trans_obj_a.num_exons >= trans_obj_b.num_exons:
            print("Error trans_a does not have fewer exons than trans_b")
            print(trans_a + " " + trans_b)
            print(str(trans_obj_a.num_exons) + " " + str(trans_obj_b.num_exons))
            sys.exit()
        
        #remove initial nocap group that is a self identity group
        if len(list(self.trans_group_dict[trans_a].keys())) == 1:# if only one group
            a_trans_group = list(self.trans_group_dict[trans_a].keys())[0]
            print(trans_a)
            if len(list(self.group_trans_dict[a_trans_group].keys())) == 1: #if only this in one group
                if list(self.group_trans_dict[a_trans_group].keys())[0] == trans_a:
                    self.group_trans_dict.pop(a_trans_group,None)
                    self.trans_group_dict.pop(trans_a,None)

        if trans_a not in self.trans_group_dict:
            self.trans_group_dict[trans_a] = {}
            
        else:
            # this happens if trans_a is a nocap trans in which case it can be in multiple groups
            print("trans_a already in group, should be nocap trans: " + trans_a)
        
        for b_group_num in list(self.trans_group_dict[trans_b].keys()):
            # add a trans to b group
            self.trans_group_dict[trans_a][b_group_num] = 1
            self.group_trans_dict[b_group_num][trans_a] = 1
   
            # search through a groups for group mergings
            for a_group_num in list(self.trans_group_dict[trans_a].keys()):
                if a_group_num == b_group_num:
                    continue
                
                print(str(a_group_num) + "-a and b group num-" + str(b_group_num))
                a_trans_id_list = list(self.group_trans_dict[a_group_num].keys())
                a_longest_trans_flag = longest_transcript(trans_a,a_trans_id_list,trans_obj_dict)
                
                if a_longest_trans_flag == "longest":
                    for a_group_trans in a_trans_id_list:
                        self.trans_group_dict[a_group_trans][b_group_num] = 1
                        self.group_trans_dict[b_group_num][a_group_trans] = 1
                        self.trans_group_dict[a_group_trans].pop(a_group_num,None)
                        self.group_trans_dict.pop(a_group_num,None)
        
            b_trans_id_list = list(self.group_trans_dict[b_group_num].keys())
            
            #redo longest trans flags
            for trans_c in self.group_trans_dict[b_group_num]:
                longest_trans_flag = longest_transcript(trans_c,b_trans_id_list,trans_obj_dict)
                self.trans_group_dict[trans_c][b_group_num] = longest_trans_flag
                self.group_trans_dict[b_group_num][trans_c] = longest_trans_flag #####################continue here!!
            
##############################################################################

####################################################################################
        
    def merge_a_b_groups(self,trans_a,trans_b):
        self.group_count += 1
        #only cap lib trans should be used for merging groups
        if len(list(self.trans_group_dict[trans_a].keys())) > 1:
            print("multiple groups a")
            sys.exit()
        if len(list(self.trans_group_dict[trans_b].keys())) > 1:
            print("multiple groups b")
            sys.exit()
            
        a_group_num = list(self.trans_group_dict[trans_a].keys())[0]
        b_group_num = list(self.trans_group_dict[trans_b].keys())[0]
        
        if self.group_count in self.group_trans_dict:
            print("group num already used")
            sys.exit()
        
        #make new group
        self.group_trans_dict[self.group_count] = {}
        
        for group_trans in self.group_trans_dict[a_group_num]:
            self.group_trans_dict[self.group_count][group_trans] = 1
            self.trans_group_dict[group_trans].pop(a_group_num, None)
            self.trans_group_dict[group_trans][self.group_count] = 1
        
        for group_trans in self.group_trans_dict[b_group_num]:
            self.group_trans_dict[self.group_count][group_trans] = 1
            self.trans_group_dict[group_trans].pop(b_group_num, None)
            self.trans_group_dict[group_trans][self.group_count] = 1
        
        #remove old groups
        self.group_trans_dict.pop(a_group_num, None)
        self.group_trans_dict.pop(b_group_num, None)
        
        
#####################################################################

    def merge_a_b_groups_nocap(self,trans_a,trans_b,trans_obj_dict):
        
        self.group_count += 1
        print("invoke merge_a_b_groups_nocap " + str(self.group_count )+ " " + trans_a + " " +trans_b )
        #only cap lib trans should be used for merging groups
        if len(list(self.trans_group_dict[trans_a].keys())) > 1:
            print("multiple groups a nocap")

        if len(list(self.trans_group_dict[trans_b].keys())) > 1:
            print("multiple groups b nocap")

        #check that trans b has more exons than trans a
        trans_obj_a = trans_obj_dict[trans_a]
        trans_obj_b = trans_obj_dict[trans_b]
        
        if trans_obj_a.num_exons != trans_obj_b.num_exons:
            print("Error trans_a does not same num exons as trans_b")
            print(trans_a + " " + trans_b)
            print(str(trans_obj_a.num_exons) + " " + str(trans_obj_b.num_exons))
            sys.exit()
        
            
        a_group_num_list = list(self.trans_group_dict[trans_a].keys())
        b_group_num_list = list(self.trans_group_dict[trans_b].keys())
        
        if self.group_count in self.group_trans_dict:
            print("group num already used")
            sys.exit()
        
        #make new group
        self.group_trans_dict[self.group_count] = {}
        
        for a_group_num in a_group_num_list:
            if self.group_trans_dict[a_group_num][trans_a] == "longest":
                for group_trans in self.group_trans_dict[a_group_num]:
                    self.group_trans_dict[self.group_count][group_trans] = 1
                    self.trans_group_dict[group_trans].pop(a_group_num, None)
                    self.trans_group_dict[group_trans][self.group_count] = 1
                #remove old group
                self.group_trans_dict.pop(a_group_num, None)
        
        for b_group_num in b_group_num_list:
            if self.group_trans_dict[b_group_num][trans_b] == "longest":
                for group_trans in self.group_trans_dict[b_group_num]:
                    self.group_trans_dict[self.group_count][group_trans] = 1
                    self.trans_group_dict[group_trans].pop(b_group_num, None)
                    self.trans_group_dict[group_trans][self.group_count] = 1
                #remove old group
                self.group_trans_dict.pop(b_group_num, None)
        
        #add longest trans information
        trans_id_list = list(self.group_trans_dict[self.group_count].keys())
        
        
        if len(trans_id_list) > 0:# new group was made
            #redo longest trans flags
            for trans_c in self.group_trans_dict[self.group_count]:
                longest_trans_flag = longest_transcript(trans_c,trans_id_list,trans_obj_dict)
                self.trans_group_dict[trans_c][self.group_count] = longest_trans_flag
                self.group_trans_dict[self.group_count][trans_c] = longest_trans_flag
        else: # now new group was made, merge did not happen, 2 short transcripts
            self.group_trans_dict.pop(self.group_count, None)
        
#####################################################################

####################################################################################################
def simplify_gene(trans_obj_list,trans_obj_dict): # goes through transcripts in gene and groups transcripts for collapsing

    transgroup = TransGroup("5cap_mix")
    
    #keep track of all no cap trans
    nocap_trans_dict = {} # nocap_trans_dict[trans_id] = 1
    
    group_num = 0
    
    for i in xrange(len(trans_obj_list)):
        trans_obj = trans_obj_list[i]
        uniq_trans_id = trans_obj.uniq_trans_id
        strand = trans_obj.strand
        
        #assume source_dict[trans_obj.source_id]['seq_type'] == no_cap or capped
        fivecap_flag = source_dict[trans_obj.source_id]['seq_type']
        
        #collect all nocap transcripts
        if fivecap_flag == "no_cap":
            nocap_trans_dict[uniq_trans_id] = 1
        
        a_group_check = transgroup.check_trans_status(uniq_trans_id)
        #make groups for each transcript if no group
        if a_group_check != 1:
            transgroup.new_group_a(uniq_trans_id)
        
        for j in xrange(i+1,len(trans_obj_list)):
            o_trans_obj = trans_obj_list[j]

            o_uniq_trans_id = o_trans_obj.uniq_trans_id
            o_strand = o_trans_obj.strand
            #assume source_dict[trans_obj.source_id]['seq_type'] == no_cap or capped
            o_fivecap_flag = source_dict[o_trans_obj.source_id]['seq_type']
            
            if uniq_trans_id == o_uniq_trans_id:
                continue
            
            #check strand match
            if strand != o_strand:
                print("Strand of transcripts within gene do not match")
                sys.exit()
            
            cap_match_flag = "none"
            
            if fivecap_flag == "capped" and o_fivecap_flag == "capped":
                cap_match_flag = "bothcap"
                trans_comp_flag,start_match_list,start_diff_list,end_match_list,end_diff_list,short_trans,min_exon_num = compare_transcripts_both_capped(trans_obj,o_trans_obj,fivecap_flag,o_fivecap_flag,strand)
            elif fivecap_flag == "capped" and o_fivecap_flag == "no_cap":
                cap_match_flag = "capnocap"
                trans_comp_flag,start_match_list,start_diff_list,end_match_list,end_diff_list,short_trans,min_exon_num = compare_transcripts_capped_nocap(trans_obj,o_trans_obj,fivecap_flag,o_fivecap_flag,strand)
            elif fivecap_flag == "no_cap" and o_fivecap_flag == "capped":
                cap_match_flag = "capnocap"
                trans_comp_flag,start_match_list,start_diff_list,end_match_list,end_diff_list,short_trans,min_exon_num = compare_transcripts_capped_nocap(trans_obj,o_trans_obj,fivecap_flag,o_fivecap_flag,strand)
            elif fivecap_flag == "no_cap" and o_fivecap_flag == "no_cap":
                #ignore nocap nocap pairings for now, will deal with single nocaps later
                continue
            else:
                print("Error with cap flags!")
                print(fivecap_flag + "\t" + o_fivecap_flag)
                sys.exit()
            
            #if uniq_trans_id == "ovary_G19720.2" and o_uniq_trans_id == "ovary_G19720.1" or o_uniq_trans_id == "ovary_G19720.2" and uniq_trans_id == "ovary_G19720.1":
            #    print(uniq_trans_id)
            #    print(o_uniq_trans_id)
            #    print(trans_comp_flag)
            #    print(start_match_list)
            #    print(end_match_list)
            #    print(end_diff_list)
            #    sys.exit()
            
            trans_match_flag = 0
            if trans_comp_flag == "same_transcript":
                trans_match_flag = 1
            elif trans_comp_flag == "same_three_prime_same_exons":
                trans_match_flag = 1
            elif trans_comp_flag == "same_three_prime_diff_exons":
                trans_match_flag = 1
            
 
            b_group_check = transgroup.check_trans_status(o_uniq_trans_id)           
            #make groups for each transcript if no group
            if b_group_check != 1:
                transgroup.new_group_a(o_uniq_trans_id)
                
            a_group_check = transgroup.check_trans_status(uniq_trans_id)
            b_group_check = transgroup.check_trans_status(o_uniq_trans_id)
            
            if trans_match_flag == 1:
                if a_group_check == 1 and b_group_check == 1: #if both are in groups already
                    
                    group_match_flag = transgroup.check_same_group(uniq_trans_id,o_uniq_trans_id)
                    
                    #if trans_group_dict[uniq_trans_id] == trans_group_dict[o_uniq_trans_id]: # if they are both in the same group
                    if group_match_flag == 1:# if they are both in the same group
                        continue
                    elif fivecap_flag == "no_cap" and o_fivecap_flag == "no_cap":# dont merge no cap trans, do this later for no cap only trans
                        continue
                    elif fivecap_flag == "no_cap" and o_fivecap_flag == "capped":# only add no cap transcript to capped group
                        transgroup.add_a_to_b_group(uniq_trans_id,o_uniq_trans_id,trans_obj_dict)
                        
                    elif fivecap_flag == "capped" and o_fivecap_flag == "no_cap":# only add no cap transcript to capped group
                        transgroup.add_a_to_b_group(o_uniq_trans_id,uniq_trans_id,trans_obj_dict)
                        
                    else: # if they are in different groups, merge groups, applies to both caps
                        transgroup.merge_a_b_groups(uniq_trans_id,o_uniq_trans_id)

                else:
                    print("Error with transcript with no group!")
                    print(uniq_trans_id)
                    print(o_uniq_trans_id)
                    sys.exit()
    
    ##########################################################################################
    ## Now deal with the single nocap transcripts
    ## only group them with themselves
    single_nocap_trans_dict = {} # single_nocap_trans_dict[nocap_trans_id] = 1
    single_nocap_trans_list = []
    
    for nocap_trans_id in nocap_trans_dict: # check for nocap trans that did not group, they can only group to cap trans
        num_assoc_trans = transgroup.check_nocap_group(nocap_trans_id)
        
        if num_assoc_trans == 0:
            #if nocap_trans_id == "oldbrain_G10.1":
            #    print(nocap_trans_id)
            #    print("no cap solo")
            #    sys.exit()
            
            single_nocap_trans_dict[nocap_trans_id] = 1
            single_nocap_trans_list.append(nocap_trans_id)
            transgroup.delete_trans(nocap_trans_id)
    
    ##########################################################################################
    ##########################################################################################
    #need to then check for nocap to nocap grouping using the single nocap trans vs all nocap trans
    for i in xrange(len(single_nocap_trans_list)):
        uniq_trans_id = single_nocap_trans_list[i]
        trans_obj = trans_obj_dict[uniq_trans_id]

        strand = trans_obj.strand
        num_exons = trans_obj.num_exons
        
        a_group_check = transgroup.check_trans_status(uniq_trans_id)
        #make groups for each transcript if no group
        # make sure to delete groups from previous grouping#################
        if a_group_check != 1:
            transgroup.new_group_a(uniq_trans_id)
                
        for j in xrange(i+1,len(single_nocap_trans_list)):
            o_uniq_trans_id = single_nocap_trans_list[j]
            o_trans_obj = trans_obj_dict[o_uniq_trans_id]
            
            o_strand = o_trans_obj.strand
            o_num_exons = o_trans_obj.num_exons

            if uniq_trans_id == o_uniq_trans_id:
                continue
            
            #check strand match
            if strand != o_strand:
                print("Strand of transcripts within gene do not match")
                sys.exit()
            
            #these are only no cap trans
            fivecap_flag = "no_cap"
            o_fivecap_flag = "no_cap"
            
            trans_comp_flag,start_match_list,start_diff_list,end_match_list,end_diff_list,short_trans,min_exon_num = compare_transcripts_both_nocap(trans_obj,o_trans_obj,fivecap_flag,o_fivecap_flag,strand)

            trans_match_flag = 0
            if trans_comp_flag == "same_transcript":
                trans_match_flag = 1
            elif trans_comp_flag == "same_three_prime_same_exons":
                trans_match_flag = 1
            elif trans_comp_flag == "same_three_prime_diff_exons":
                trans_match_flag = 1
            
            b_group_check = transgroup.check_trans_status(o_uniq_trans_id)
            #make groups for each transcript if no group
            # make sure to delete groups from previous grouping#################
            if b_group_check != 1:
                transgroup.new_group_a(o_uniq_trans_id)
                
            a_group_check = transgroup.check_trans_status(uniq_trans_id)
            b_group_check = transgroup.check_trans_status(o_uniq_trans_id)
            
            if trans_match_flag == 1:
                if a_group_check == 1 and b_group_check == 1: #if both are in groups already
                    
                    group_match_flag = transgroup.check_same_group(uniq_trans_id,o_uniq_trans_id)
                    
                    #if trans_group_dict[uniq_trans_id] == trans_group_dict[o_uniq_trans_id]: # if they are both in the same group
                    if group_match_flag == 1:# if they are both in the same group
                        continue
                    else: # if they are in different groups
                        if num_exons == o_num_exons:   # same number of exons
                            transgroup.merge_a_b_groups_nocap(uniq_trans_id,o_uniq_trans_id,trans_obj_dict)
                            
                        elif num_exons > o_num_exons: #add shorter to longer
                            transgroup.add_a_to_b_group_both_nocap(o_uniq_trans_id,uniq_trans_id,trans_obj_dict)
                        elif o_num_exons > num_exons: # add shorter to longer
                            transgroup.add_a_to_b_group_both_nocap(uniq_trans_id,o_uniq_trans_id,trans_obj_dict)
                        else:
                            print("Error with nocap comparison, num exons issue")
                            sys.exit()
                            
                else:
                    print("Error with transcript with no group! - Nocap")
                    print(uniq_trans_id)
                    print(o_uniq_trans_id)
                    sys.exit()
    ###############################################################################
    ############################################################################### Edge of work!!!

    trans_group_dict = transgroup.trans_group_dict
    group_trans_dict = transgroup.group_trans_dict
        
    return trans_group_dict,group_trans_dict

####################################################################################################
def reverse_complement(seq_list):
    comp_dict = {}
    comp_dict["A"] = "T"
    comp_dict["T"] = "A"
    comp_dict["G"] = "C"
    comp_dict["C"] = "G"
    comp_dict["N"] = "N"
    
    reverse_seq = seq_list[::-1] # that's a neat way of reversing the string
    
    rev_comp_list = []
    for base in reverse_seq:
        rev_comp_list.append(comp_dict[base])
    
    
    #rev_comp_string = "".join(rev_comp_list)
    
    return rev_comp_list

def format_gene_report_line(trans_obj_list,total_gene_count,total_final_trans_count):
    
    gene_id = "G" + str(total_gene_count)
    
    
    gene_start = -1
    gene_end = 0
    gene_chrom = "none"
    
    gene_source_dict = {} # gene_source_dict[source] = 1
    
    num_source_trans = 0
    
    for trans_obj in trans_obj_list:
        num_source_trans += 1
        
        if gene_start == -1:
            gene_start = int(trans_obj.start_pos)
        if int(trans_obj.start_pos) < gene_start:
            gene_start = int(trans_obj.start_pos)
        
        if gene_end == 0 :
            gene_end = int(trans_obj.end_pos)
        if int(trans_obj.end_pos) > gene_end:
            gene_end = int(trans_obj.end_pos)
        
        if gene_chrom == "none":
            gene_chrom = trans_obj.scaffold
        
        if gene_chrom != trans_obj.scaffold:
            print("Error with different chroms for gene group")
            sys.exit()
        
        gene_source_dict[trans_obj.source_id] = 1
    
    gene_source_list = list(gene_source_dict.keys())
    gene_source_list.sort()
    
    gene_source_line = ",".join(gene_source_list)
    
    gene_report_list = []
    gene_report_list.append(gene_id)
    gene_report_list.append(str(num_source_trans))
    gene_report_list.append(str(total_final_trans_count))
    gene_report_list.append(gene_source_line)
    gene_report_list.append(gene_chrom)
    gene_report_list.append(str(gene_start))
    gene_report_list.append(str(gene_end))
    
    gene_report_line = "\t".join(gene_report_list)
    
    
    return gene_report_line


def process_trans_group(trans_line_list, total_gene_count):
    # this takes grouped transcripts and processes them to produce bed file output
    # deals with gene group, strand separation, collapsing, and sorting
    
    forward_trans_list = []
    reverse_trans_list = []
    
    trans_obj_dict = {} # trans_obj_dict[trans id] = trans obj
    #merged_obj_dict = {}
    
    
    for trans_line_split in trans_line_list:
        trans_line = "\t".join(trans_line_split)
        id_line = trans_line_split[3]
        id_line_split = id_line.split(";")
        gene_id = id_line_split[0]
        trans_id = id_line_split[1]
         
        source_id = trans_line_split[12]
        uniq_trans_id = "_".join([source_id,trans_id])
        
        [start_priority,junction_priority,end_priority] = source_dict[source_id]['priority_rank'].split(",")
        
        trans_obj = Transcript(trans_id)
        trans_obj.add_bed_info(trans_line)
        trans_obj.add_source_id(source_id)
        trans_obj.add_priority(start_priority,junction_priority,end_priority)
        
        if uniq_trans_id in trans_obj_dict:
            print("Error with duplicate trans id")
            sys.exit()
            
        trans_obj_dict[uniq_trans_id] = trans_obj
        
        if trans_obj.strand == "+":
            forward_trans_list.append(trans_obj)
        elif trans_obj.strand == "-":
            reverse_trans_list.append(trans_obj)
        else:
            print("Error with strand information")
            sys.exit()
    
    forward_gene_start_trans_dict,forward_start_gene_list = gene_group(forward_trans_list)
    reverse_gene_start_trans_dict,reverse_start_gene_list = gene_group(reverse_trans_list)

        
    all_start_gene_dict = {} #     all_start_gene_dict[start] = 1
    
    #collect all starts forward and reverse strands
    for gene_start in forward_start_gene_list:
        all_start_gene_dict[gene_start] = 1
    for gene_start in reverse_gene_start_trans_dict:
        all_start_gene_dict[gene_start] = 1
    
    all_start_list = all_start_gene_dict.keys()
    
    all_start_list.sort()
    
    for gene_start in all_start_list:
        gene_trans_obj_list = [] #list of trans obj lists

        #if a forward and reverse gene start at the same place use this to make the the forward strand gene is represented first
        if gene_start in forward_gene_start_trans_dict:
            uniq_trans_id_list = forward_gene_start_trans_dict[gene_start].keys()
            trans_obj_list = []
            for uniq_trans_id in uniq_trans_id_list:
                trans_obj_list.append(trans_obj_dict[uniq_trans_id])
            gene_trans_obj_list.append(trans_obj_list)
        
        if gene_start in reverse_gene_start_trans_dict:
            uniq_trans_id_list = reverse_gene_start_trans_dict[gene_start].keys()
            trans_obj_list = []
            for uniq_trans_id in uniq_trans_id_list:
                trans_obj_list.append(trans_obj_dict[uniq_trans_id])
            gene_trans_obj_list.append(trans_obj_list)
            
        #loop through list of trans obj lists, usually only one list since unlikely for forward and reverse strand genes to coincide
        for trans_obj_list in gene_trans_obj_list:      
            total_gene_count += 1
            
            
            
            #group transcripts by collapsability
            match_trans_group_dict,match_group_trans_dict = simplify_gene(trans_obj_list,trans_obj_dict)

            #print(list(trans_obj_dict.keys()))
            #print(list(match_trans_group_dict.keys()))
            #print(match_group_trans_dict)
            
            merge_obj_list = []
            tmp_count = 0
            #go through matched groups that should be collapase
            for match_group_num in match_group_trans_dict:
                tmp_count += 1
                tmp_trans_id = "G" + str(total_gene_count) + ".tmp." + str(tmp_count)
                merged_obj = Merged(tmp_trans_id)
                
                match_trans_id_list = match_group_trans_dict[match_group_num].keys()
                match_trans_obj_list = []
                for match_trans_id in match_trans_id_list:
                    match_trans_obj = trans_obj_dict[match_trans_id]
                    match_trans_obj_list.append(match_trans_obj)
                    
                    merged_obj.add_merged_trans(match_trans_obj)
                    
                    print(match_trans_id)
                
                redundant_trans_flag = 0
                if len(match_trans_obj_list) > 1: #if there are redundant transcripts, collapse
                    redundant_trans_flag = 1
                    collapse_flag = "common_ends"
                    collapse_start_list,collapse_end_list,start_wobble_list,end_wobble_list,e_start_trans_dict,e_end_trans_dict = collapse_transcripts(match_trans_obj_list,collapse_flag)
                    
                    merged_obj.add_merge_info(collapse_start_list,collapse_end_list,start_wobble_list,end_wobble_list,e_start_trans_dict,e_end_trans_dict )
                else: # if only one transcript
                    #print(tmp_trans_id)
                    #print(len(match_trans_obj_list))
                    exon_start_list = match_trans_obj_list[0].exon_start_list
                    exon_end_list = match_trans_obj_list[0].exon_end_list
                    start_wobble_list = [0] * len(exon_start_list)
                    end_wobble_list = [0] * len(exon_start_list)
                    
                    exon_start_list.sort()
                    exon_end_list.sort()
                    
                    match_trans_obj_list[0].uniq_trans_id
                    e_start_trans_dict = {}
                    e_end_trans_dict = {}
                    for z in xrange(len(exon_start_list)):
                        e_start_trans_dict[exon_start_list[z]] = {}
                        e_start_trans_dict[exon_start_list[z]][match_trans_obj_list[0].uniq_trans_id] = 1
                        
                        e_end_trans_dict[exon_end_list[z]] = {}
                        e_end_trans_dict[exon_end_list[z]][match_trans_obj_list[0].uniq_trans_id] = 1
                        
                        

                    merged_obj.add_merge_info(exon_start_list,exon_end_list,start_wobble_list,end_wobble_list,e_start_trans_dict,e_end_trans_dict)
                
                merge_obj_list.append(merged_obj)
                
            
            sorted_merge_obj_list = sort_transcripts(merge_obj_list,trans_obj_dict)
            
            trans_count = 0
            for merged_obj in sorted_merge_obj_list:
                trans_count += 1
                final_trans_id = "G" + str(total_gene_count) + "." + str(trans_count)
                merged_obj.trans_id = final_trans_id
                
                #merged_obj_dict[final_trans_id] = merged_obj
                
                #write out to bed file
                bed_line = merged_obj.format_bed_line()
                outfile_bed.write(bed_line)
                outfile_bed.write("\n")
                
                #write out to transcript report file
                trans_report_line = merged_obj.format_trans_report_line()
                outfile_trans_report.write(trans_report_line)
                outfile_trans_report.write("\n")

                #write out to transcript merge report file
                for merged_trans_id in merged_obj.merged_trans_dict:
                    merged_trans_obj = merged_obj.merged_trans_dict[merged_trans_id]
                    trans_id = merged_trans_obj.trans_id
                    scaffold = merged_trans_obj.scaffold
                    strand = merged_trans_obj.strand
                    trans_start = merged_trans_obj.trans_start
                    trans_end = merged_trans_obj.trans_end
                    
                    exon_start_line,exon_end_line = merged_trans_obj.make_exon_start_end_lines()
                    merge_list = []
                    merge_list.append(final_trans_id)
                    merge_list.append(merged_trans_id)
                    merge_list.append(scaffold)
                    merge_list.append(strand)
                    merge_list.append(str(trans_start))
                    merge_list.append(str(trans_end))
                    merge_list.append(exon_start_line)
                    merge_list.append(exon_end_line)
                
                    merge_line = "\t".join(merge_list)
                    
                    merge_trans_bed_line = merged_trans_obj.format_bed_line(final_trans_id)
                    
                    outfile_merge.write(merge_trans_bed_line)
                    outfile_merge.write("\n")
                    
            gene_report_line = format_gene_report_line(trans_obj_list,total_gene_count, trans_count) #######################################
            outfile_gene_report.write(gene_report_line)
            outfile_gene_report.write("\n")
                           
    return total_gene_count
                

        
        
        ##############################################convert into trans object !!!!!!!!!!!!!!!!
    
#end of def
###################################################################################################
###################################################################################################
###################################################################################################


for file_line in filelist_file_contents:
    file_line_split = file_line.split("\t")
    
    filename = file_line_split[0]
    seq_type = file_line_split[1]
    priority_rank = file_line_split[2]
    source_id = file_line_split[3]
    
    if source_id not in source_dict:
        source_dict[source_id] = {}
        source_dict[source_id]['filename'] = filename
        source_dict[source_id]['seq_type'] = seq_type
        source_dict[source_id]['priority_rank'] = priority_rank # by type of merge event 1,1,1 ->  start,junction,end
        
        if seq_type != "capped" and seq_type != "no_cap":
            print("Incorrect seq types given. Please use capped or  no_cap")
    
    print("opening bed list")
    bed_file_contents = open(filename).read().rstrip("\n").split("\n")
    
    for line in bed_file_contents:
        line_split = line.split("\t")
        
        scaffold = line_split[0]
        trans_start = int(line_split[1])
        trans_end = int(line_split[2])
        id_line = line_split[3]
        strand = line_split[5]
        num_exon = int(line_split[9])
        block_size_list = line_split[10].split(",")
        block_start_list = line_split[11].split(",")
        
        line_split.append(source_id)
        
        if scaffold not in bed_dict:
            bed_dict[scaffold] = {}
        if trans_start not in bed_dict[scaffold]:
            bed_dict[scaffold][trans_start] = {}
        if trans_end not in bed_dict[scaffold][trans_start]:
            bed_dict[scaffold][trans_start][trans_end] = []
        
        bed_dict[scaffold][trans_start][trans_end].append(line_split)
        
scaff_list = list(bed_dict.keys())
scaff_list.sort()


#####################################################
colour_source_list = list(source_dict.keys())
colour_source_list.sort()
colour_source_count = 0
source_colour_dict = {} # source_colour_dict[source] = colour number
for colour_source in colour_source_list:
    colour_source_count += 1
    source_colour_dict[colour_source] = colour_source_count
#####################################################



total_gene_count = 0

for scaffold in scaff_list:
    
    trans_start_list = list(bed_dict[scaffold].keys())
    trans_start_list.sort()
    
    num_trans_start = len(trans_start_list)
    
    trans_start_count  = 0
    group_trans_list = []
    
    last_trans_end = 0
    for trans_start in trans_start_list:
        trans_start_count += 1
        trans_end_list = list(bed_dict[scaffold][trans_start].keys())
        trans_end_list.sort()
        if last_trans_end == 0:
            last_trans_end = trans_end_list[0]
        
        if trans_start <= last_trans_end: # matches previous group

            for trans_end in trans_end_list:
                # add transcripts to group
                group_trans_list.extend(bed_dict[scaffold][trans_start][trans_end])
                if last_trans_end < trans_end: # update last trans end
                    last_trans_end = trans_end
        
        elif trans_start > last_trans_end: # if this is a new group, process past group
            
            total_gene_count = process_trans_group(group_trans_list,total_gene_count)
            
            #refresh group list
            group_trans_list = []
            #update last trans end
            for trans_end in trans_end_list:
                # add transcripts to group
                group_trans_list.extend(bed_dict[scaffold][trans_start][trans_end])
                if last_trans_end < trans_end: # update last trans end
                    last_trans_end = trans_end
        
        
        # after dealing with groups check if this is the last group in scaffold
        if trans_start_count == num_trans_start:  # last trans
            #update last trans end            
            total_gene_count = process_trans_group(group_trans_list,total_gene_count)
            
            
        
            
    


        

    
    
