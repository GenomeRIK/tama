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

Last Updated: 2020/12/17

Added more informative error message for issues with strand information from input bed files. 

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

ap.add_argument('-d', type=str, nargs=1, help='Flag for merging duplicate transcript groups (default no_merge quits when duplicates are found, merge_dup will merge duplicates)')

ap.add_argument('-s', type=str, nargs=1, help='Use gene and transcript ID from a merge source. Specify source name from filelist file here.')

ap.add_argument('-cds', type=str, nargs=1, help='Use CDS from a merge source. Specify source name from filelist file here.')

opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0

if not opts.f:
    print("Filelist file missing")
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

if not opts.d:
    print("Default duplicate merge flag: no_merge")
    duplicate_flag = "no_merge"
else:
    duplicate_flag = str(opts.d[0])

if not opts.s:
    print("Default source ID merge flag: no_source_id")
    source_id_flag = "no_source_id"
else:
    source_id_flag = str(opts.s[0])

    source_id_flag_list = source_id_flag.split(",")

if not opts.cds:
    print("Default CDS merge flag: no_cds")
    source_cds_flag = "no_cds"
else:
    source_cds_flag = str(opts.cds[0])

    source_cds_flag_list = source_cds_flag.split(",")


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
trans_report_line = "\t".join(["transcript_id","num_clusters","sources","start_wobble_list","end_wobble_list","exon_start_support","exon_end_support","all_source_trans"])
outfile_trans_report.write(trans_report_line)
outfile_trans_report.write("\n")

gene_report_outfile_name = outfile_prefix + "_gene_report.txt"
outfile_gene_report = open(gene_report_outfile_name,"w")
gene_report_line = "\t".join(["gene_id","num_clusters","num_final_trans","sources","chrom", "start","end","source_genes","source_summary"])
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
        
        #add this for commas at the end of the exon block and start strings
        if block_size_list[-1] == "":
            block_size_list.pop(-1)
        
        if block_start_list[-1] == "":
            block_start_list.pop(-1)
        
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

        self.extra_trans_id_list = "NA"
        self.extra_gene_id_list = "NA"
        
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

        self.start_cds = 0
        self.end_cds = 0
        
        self.num_exons = 0
        
        #use these to see which transcripts support which coordinates
        self.e_start_trans_dict = {}
        self.e_end_trans_dict = {}
        
        self.e_start_trans_list = []
        self.e_end_trans_list = []

    
    def add_merged_trans(self,trans_obj):

        merged_trans_id = trans_obj.uniq_trans_id

        if merged_trans_id in self.merged_trans_dict:
            print("Transcript is already in this group")
            return


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

        #print(self.trans_list)
        #print(collapse_start_list)
        #print(collapse_end_list)
        
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


        if self.extra_trans_id_list == "NA":
            id_line = ";".join([gene_id,self.trans_id])
        else:
            this_id_line_list = []
            this_id_line_list.append(gene_id)
            this_id_line_list.append(self.trans_id)

            for k in xrange(len(self.extra_gene_id_list)):
                this_id_line_list.append(self.extra_gene_id_list[k])
                this_id_line_list.append(self.extra_trans_id_list[k])

            id_line = ";".join(this_id_line_list)


        bed_list.append(str(id_line))
        bed_list.append("40")
        bed_list.append(self.strand)

        if self.start_cds != 0:
            bed_list.append(str(self.start_cds))
            bed_list.append(str(self.end_cds))
        
        else:
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

        ####################
        # 2020/05/31
        all_source_trans_line = ",".join(self.trans_list)

        trans_report_list.append(all_source_trans_line)

        # 2020/05/31
        ###################
       
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
        
        e_start_dict = {} # e_start_dict[priority number][start] = number of occurrences
        e_end_dict = {} # e_end_dict[priority number][end] = number of occurrences

        e_start_priority_dict = {} #e_start_priority_dict[e start priority] = 1
        e_end_priority_dict = {} #e_end_priority_dict[e end priority] = 1

        for trans_obj in trans_obj_list:
            e_start_list = trans_obj.exon_start_list
            e_end_list = trans_obj.exon_end_list
            
            # figure out priority for different features
            start_priority = trans_obj.start_priority
            end_priority = trans_obj.end_priority 
            junction_priority = trans_obj.junction_priority

            this_cap_flag = source_dict[trans_obj.source_id]['seq_type']
            
            if i >= len(e_start_list):# use for no cap when exon numbers may not match
                continue
            
            e_start = int(e_start_list[j])
            e_end = int(e_end_list[j])

            
            # if this is the first 5' exon of group, use 5' start priority, and junction end priority
            if i == max_exon_num - 1:
                if strand == "+":
                    e_start_priority = start_priority
                    e_end_priority = junction_priority

                    e_start_priority_dict[e_start_priority] = 1
                    e_end_priority_dict[e_end_priority] = 1

                if strand == "-": # because e end is trans start
                    e_start_priority = junction_priority
                    e_end_priority = start_priority

                    e_start_priority_dict[e_start_priority] = 1
                    e_end_priority_dict[e_end_priority] = 1

            # if this is the first exon for this read and not the first for the transcript group then set priority to last for 5' end of exon
            elif i <  max_exon_num - 1 and i == len(e_start_list)-1 :
                if this_cap_flag != "no_cap": # this situation should only happen with nocap reads
                    print("Error with capped transcript treated like no_cap")
                    print(trans_obj.trans_id)
                    sys.exit()

                if strand == "+":
                    e_start_priority = 999 # force e start priority to last
                    e_end_priority = junction_priority

                    e_start_priority_dict[e_start_priority] = 1
                    e_end_priority_dict[e_end_priority] = 1

                if strand == "-":

                    e_start_priority = junction_priority
                    e_end_priority = 999 # force e end priority to last

                    e_start_priority_dict[e_start_priority] = 1
                    e_end_priority_dict[e_end_priority] = 1


            # if this is a junction exon then use junction priorities
            elif i < max_exon_num - 1 and i < len(e_start_list)-1 and i > 0:
                e_start_priority = junction_priority
                e_end_priority = junction_priority

                e_start_priority_dict[e_start_priority] = 1
                e_end_priority_dict[e_end_priority] = 1

            
            # if this is the 3' end of the transcript use end priority and junction priority for start
            elif i == 0 :
                if strand == "+":
                    e_start_priority = junction_priority
                    e_end_priority = end_priority

                    e_start_priority_dict[e_start_priority] = 1
                    e_end_priority_dict[e_end_priority] = 1

                if strand == "-":

                    e_start_priority = end_priority
                    e_end_priority = junction_priority

                    e_start_priority_dict[e_start_priority] = 1
                    e_end_priority_dict[e_end_priority] = 1
            else:
                print("no assignment of exon number within transcript group")
                print(trans_obj.trans_id)
                sys.exit()

            if e_start != -1:
                if e_start_priority not in e_start_dict:
                    e_start_dict[e_start_priority] = {}
                if e_start not in e_start_dict[e_start_priority]:
                    e_start_dict[e_start_priority][e_start] = 0
                    e_start_trans_dict[e_start] = {}
                e_start_dict[e_start_priority][e_start] += 1
                e_start_trans_dict[e_start][trans_obj.uniq_trans_id] = 1
            
            if e_end != -1:
                if e_end_priority not in e_end_dict:
                    e_end_dict[e_end_priority] = {}
                if e_end not in e_end_dict[e_end_priority]:
                    e_end_dict[e_end_priority][e_end] = 0
                    e_end_trans_dict[e_end] = {}
                e_end_dict[e_end_priority][e_end] += 1
                e_end_trans_dict[e_end][trans_obj.uniq_trans_id] = 1
        
        
        ##########################################
        best_e_start = -1
        long_e_start = -1
        short_e_start = -1
        num_starts = 0

        start_priority_list = list(e_start_priority_dict.keys())
        start_priority_list.sort()
        best_start_priority = start_priority_list[0]

        end_priority_list = list(e_end_priority_dict.keys())
        end_priority_list.sort()
        best_end_priority = end_priority_list[0]

        for e_start in e_start_dict[best_start_priority]:
            if e_start_dict[best_start_priority][e_start] > num_starts:
                best_e_start = e_start
                num_starts = e_start_dict[best_start_priority][e_start]

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


        for e_start in e_start_dict[best_start_priority]:
            if e_start_dict[best_start_priority][e_start] == most_starts:
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
        for e_end in e_end_dict[best_end_priority]:
            if e_end_dict[best_end_priority][e_end] > num_ends:
                best_e_end = e_end
                num_ends = e_end_dict[best_end_priority][e_end]
            
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
        for e_end in e_end_dict[best_end_priority]:
            if e_end_dict[best_end_priority][e_end] == most_ends:
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

####################################################################################################

def iterate_sort_list(list_trans_pos_list,pos_index):
    # sort the list by each element
    sort_flag = 0

    blah_flag = 0

    while sort_flag == 0:

        pre_pos = -1

        same_order_index_dict = {}  # same_order_index_dict[pos][index] = 1

        # collect positions and index where the sort was equal
        for j in xrange(len(list_trans_pos_list)):

            trans_pos_line_split = list_trans_pos_list[j]
            pos_element = trans_pos_line_split[pos_index]

            if pos_element == "971115":
                blah_flag = 1


            if pos_element == pre_pos:
                if pre_pos not in same_order_index_dict:
                    same_order_index_dict[pre_pos] = {}
                    same_order_index_dict[pre_pos][j-1] = 1

                if pos_element not in same_order_index_dict:
                    same_order_index_dict[pos_element] = {}
                same_order_index_dict[pos_element][j] = 1

            pre_pos = pos_element

        same_order_index_list = list(same_order_index_dict.keys())
        same_order_index_list.sort()

        if len(same_order_index_list) == 0:
            sort_flag = 1
        else:
            for pos_element in same_order_index_list:
                slice_pos_index_list = list(same_order_index_dict[pos_element].keys())
                slice_pos_index_list.sort()

                # check for iterative index
                past_pos = "na"
                for this_pos in slice_pos_index_list:
                    if past_pos == "na":
                        past_pos = this_pos

                    elif this_pos != past_pos + 1:
                        print("Error with non-consequtive indices for same pos in trans sorting.")
                        print(slice_pos_index_list)
                        sys.exit()
                    else:
                        past_pos = this_pos

                min_index = slice_pos_index_list[0]
                max_index = slice_pos_index_list[-1] + 1

                # replace slice with ordered slice
                list_slice = list_trans_pos_list[min_index:max_index]
                list_slice.sort(key=lambda x: int(x[pos_index + 1]))

                if blah_flag == 1:
                    print("blah blah 1")
                    print(list_slice)
                    print(min_index)
                    print(max_index)

                [list_slice, sort_flag] = iterate_sort_list(list_slice, pos_index + 1)

                if blah_flag == 1:
                    print("blah blah 2")
                    print(list_slice)
                    print(pos_index)

                list_trans_pos_list[min_index:max_index] = list_slice

    if blah_flag == 1:
        print("blah blah 3")
        print(list_trans_pos_list)
#        sys.exit()

    return [list_trans_pos_list,sort_flag]


####################################################################################################

def sort_pos_trans_list(pos_trans_list,pos_trans_dict):
    # sort the pos_trans_list according to position of transcripts on chromosome

    sorted_trans_list = []

    new_pos_trans_dict = {}

    list_trans_pos_list = [] # this is a list of lists

    max_pos_num = 0

    #get max number of elements for pos lists
    for trans_pos_line in pos_trans_list:
        trans_pos_line_split = trans_pos_line.split(",")

        if max_pos_num < len(trans_pos_line_split):
            max_pos_num = len(trans_pos_line_split)

    for trans_pos_line in pos_trans_list:
        trans_pos_line_split = trans_pos_line.split(",")

        diff_pos = max_pos_num - len(trans_pos_line_split)

        # pad out list so all pos lists have same number of elements
        for i in xrange(diff_pos):
            trans_pos_line_split.append(0)

        trans_pos_line_split_str = []
        for pos_element in trans_pos_line_split:
            trans_pos_line_split_str.append(str(pos_element))

        new_trans_pos_line = ",".join(trans_pos_line_split_str)

        new_pos_trans_dict[new_trans_pos_line] = pos_trans_dict[trans_pos_line]

        list_trans_pos_list.append(trans_pos_line_split)

    #sort the list by each element
    i = 0

    list_trans_pos_list.sort(key=lambda x: int(x[i]))

    pos_index = i

    [list_trans_pos_list, sort_flag] = iterate_sort_list(list_trans_pos_list, pos_index)

    if sort_flag != 1:
        print("Error with sort flag!")
        print(sort_flag)
        sys.exit()

    for trans_pos_line_split in list_trans_pos_list:
        trans_pos_line_split_str = []
        for pos_element in trans_pos_line_split:
            trans_pos_line_split_str.append(str(pos_element))

        new_trans_pos_line = ",".join(trans_pos_line_split_str)
        sorted_trans_list.append(new_trans_pos_line)


    return [sorted_trans_list,new_pos_trans_dict]


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

        trans_pos_list = []

        #trans_pos_list.append(str(trans_start))
        #trans_pos_list.append(",")
        #trans_pos_list.append(str(trans_end))
        #trans_pos_list.append(",")

        num_exons = len(trans_exon_start_list)

        for i in xrange(num_exons):
            exon_start = trans_exon_start_list[i]
            trans_pos_list.append(str(exon_start))
            trans_pos_list.append(",")

            exon_end = trans_exon_end_list[i]
            trans_pos_list.append(str(exon_end))
            trans_pos_list.append(",")

        # remove last element because it is a comma
        trans_pos_list.pop(-1)

        trans_pos_line = "".join(trans_pos_list)
        
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

            if duplicate_flag == "no_merge":
                print("By default TAMA merge does not allow merging of duplicate transcript groups.")
                print("Duplicate transcript groups occur when different groupings of transcripts results in the same collapsed model.")
                print("If you would like to merge duplicate transcript groups please add -d merge_dup to the arguments.")
                sys.exit()
            elif duplicate_flag == "merge_dup":

                # add trans obj to merge group
                for a_uniq_trans_id in trans_obj.merged_trans_dict:

                    pos_trans_dict[trans_pos_line].add_merged_trans(trans_obj.merged_trans_dict[a_uniq_trans_id])

                match_trans_obj_list = []

                # collect trans obj in list
                for b_uniq_trans_id in pos_trans_dict[trans_pos_line].merged_trans_dict:
                    match_trans_obj_list.append(pos_trans_dict[trans_pos_line].merged_trans_dict[b_uniq_trans_id])

                collapse_start_list, collapse_end_list, start_wobble_list, end_wobble_list, e_start_trans_dict, e_end_trans_dict = collapse_transcripts(match_trans_obj_list, collapse_flag)

                # update merge info
                pos_trans_dict[trans_pos_line].add_merge_info(collapse_start_list, collapse_end_list, start_wobble_list,end_wobble_list, e_start_trans_dict, e_end_trans_dict)

                trans_obj = pos_trans_dict[trans_pos_line]

            else:
                print("Error with duplicate transcript group flag")
                sys.exit()
        
        pos_trans_dict[trans_pos_line] = trans_obj
        #pos_trans_list.append(trans_pos_line)

    pos_trans_list = list(pos_trans_dict.keys())

    [new_pos_trans_list, new_pos_trans_dict] = sort_pos_trans_list(pos_trans_list, pos_trans_dict)
    
    
    sorted_trans_obj_list = []
    for pos_trans in new_pos_trans_list:
        
        trans_obj = new_pos_trans_dict[pos_trans]
        sorted_trans_obj_list.append(trans_obj)
        tmp_id = trans_obj.trans_id

    
    # sorted_trans_obj_list list of trana obj that have been sorted by position
    return sorted_trans_obj_list

##############################################################################


def longest_transcript_old(trans_id,trans_id_list,trans_obj_dict):
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


def longest_transcript(trans_id,trans_id_list,trans_obj_dict):
    #determines if transcript is one of the longest in group on the 5 prime end
    # if longest_trans == "longest" than this is one of the longest trans
    # if "long" then it is within range3 of the longest but it is not the longest
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


    ##############################################
    #new algorithm
    longest_num_exons = 0
    longest_tss_pos = -1
    longest_trans_id_dict = {}  # longest_trans_id_dict[trans id] = longest_tss_pos
    longest_trans_id = "none"

    #find the longest in the list first
    for o_trans_id in trans_id_list:
        o_trans_obj = trans_obj_dict[o_trans_id]

        o_trans_start = o_trans_obj.start_pos
        o_trans_end = o_trans_obj.end_pos
        o_num_exons = o_trans_obj.num_exons

        if o_num_exons > longest_num_exons:
            longest_num_exons = o_num_exons

        if o_num_exons == longest_num_exons:
            if strand == "+":
                if longest_tss_pos == -1:
                    longest_tss_pos = o_trans_start
                    longest_trans_id = o_trans_id
                elif o_trans_start < longest_tss_pos:
                    longest_tss_pos = o_trans_start
                    longest_trans_id = o_trans_id
            elif strand == "-":
                if longest_tss_pos == -1:
                    longest_tss_pos = o_trans_end
                    longest_trans_id = o_trans_id
                elif o_trans_end > longest_tss_pos:
                    longest_tss_pos = o_trans_end
                    longest_trans_id = o_trans_id
            else:
                print("Error with strands in longest_transcript")
                sys.exit()

    #compare the longest in the list to the trans_id
    if strand == "+":
        if num_exons > longest_num_exons:
            longest_trans = "longest"
        elif num_exons == longest_num_exons:
            if trans_start < longest_tss_pos:
                longest_trans = "longest"
            elif trans_start == longest_tss_pos:
                longest_trans = "longest"
            elif trans_start > longest_tss_pos:
                start_match_flag,start_diff_num = fuzzy_match(trans_start,longest_tss_pos,exon_diff_threshold)
                if start_match_flag == "wobbly_match":
                    longest_trans = "long"
                elif start_match_flag == "no_match":
                    longest_trans = "short"
                else:
                    print("Error with fuzzy_match in longest_transcript")
                    print(start_match_flag)
                    sys.exit()
            else:
                print("Error with trans_start/longest_tss_pos comparison in longest_transcript ")
                sys.exit()
        elif num_exons < longest_num_exons:
            longest_trans = "short"
        else:
            print("Error with num exons comparison in longest_transcript ")
            sys.exit()

    elif strand == "-":
        if num_exons > longest_num_exons:
            longest_trans = "longest"
        elif num_exons == longest_num_exons:
            if trans_end > longest_tss_pos:
                longest_trans = "longest"
            elif trans_end == longest_tss_pos:
                longest_trans = "longest"
            elif trans_end < longest_tss_pos:
                start_match_flag,start_diff_num = fuzzy_match(trans_end,longest_tss_pos,exon_diff_threshold)
                if start_match_flag == "wobbly_match":
                    longest_trans = "long"
                elif start_match_flag == "no_match":
                    longest_trans = "short"
                else:
                    print("Error with fuzzy_match in longest_transcript")
                    print(start_match_flag)
                    sys.exit()
            else:
                print("Error with trans_start/longest_tss_pos comparison in longest_transcript ")
                sys.exit()
        elif num_exons < longest_num_exons:
            longest_trans = "short"
        else:
            print("Error with num exons comparison in longest_transcript ")
            sys.exit()
    else:
        print("Error with strand in longest_transcript")
        sys.exit()

    #new algorithm
    ##############################################

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

        self.group_max_exon_dict = {} # group_max_exon_dict[group num] = max exon num


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
    
    def check_trans_group(self,trans_a):
        # may add on for future features
        single_group_check = 0

        for trans_a in self.trans_group_dict:
            
            if len(list(self.trans_group_dict[trans_a].keys())) == 1:
                group_a = list(self.trans_group_dict[trans_a].keys())[0]
                
                if len(list(self.group_trans_dict[group_a].keys())) == 1:
                    if trans_a in self.group_trans_dict[group_a]:
                        single_group_check = 1
                
        return single_group_check

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

    def new_group_a(self,trans_a,trans_obj_dict):

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

        self.group_longest_dict[self.group_count] = {}
        self.group_longest_dict[self.group_count]["longest"] = {}
        self.group_longest_dict[self.group_count]["longest"][trans_a] = 1

        self.group_max_exon_dict[self.group_count] =  trans_obj_dict[trans_a].num_exons ######################

        #print("new_group_a")
        #print(self.group_max_exon_dict)

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
        # add trans_a to b group and if a trans_a is longest nocap trans in it's group add other nocap trans
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
            #print(trans_a)
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

    def add_a_to_b_group_nocap(self,trans_a,trans_b):
        #if log_flag == "log_on":
        #    print("invoke add_a_to_b_group " + trans_a + " " + trans_b)
        # only cap libs should be used for b group and they only have one group
        # does not take all of a_group just uses a_trans
        #if len(list(self.trans_group_dict[trans_b].keys())) > 1:
        #    if log_flag == "log_on":
        #        print("multiple groups")
                #sys.exit()

        #check that trans b has more exons than trans a
        #trans_obj_a = trans_obj_dict[trans_a]
        #trans_obj_b = trans_obj_dict[trans_b]
        
        
        #remove initial nocap group that is a self identity group
        if len(list(self.trans_group_dict[trans_a].keys())) == 1:# if only one group for trans_a
            a_trans_group = list(self.trans_group_dict[trans_a].keys())[0]
            #if log_flag == "log_on":
            #    print(trans_a)
            if len(list(self.group_trans_dict[a_trans_group].keys())) == 1: #if only one transcript in one group
                if list(self.group_trans_dict[a_trans_group].keys())[0] == trans_a: # if trans_a is just a group by itself
                    self.group_trans_dict.pop(a_trans_group,None)
                    self.trans_group_dict.pop(trans_a,None)

        if trans_a not in self.trans_group_dict:
            self.trans_group_dict[trans_a] = {}
            
        #else:
            # this happens if trans_a is a nocap trans in which case it can be in multiple groups
            #if log_flag == "log_on":
            #    print("trans_a already in group, should be nocap trans: " + trans_a)
        

        for b_group_num in list(self.trans_group_dict[trans_b].keys()):
            # add a trans to b group
            self.trans_group_dict[trans_a][b_group_num] = "short"
            self.group_trans_dict[b_group_num][trans_a] = "short"
            # a_trans has fewer exons than b_trans so it must be short
                
   

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

        # add group max exons
        a_group_max_exons = self.group_max_exon_dict[a_group_num]
        b_group_max_exons = self.group_max_exon_dict[b_group_num]

        if a_group_max_exons < b_group_max_exons:
            self.group_max_exon_dict[self.group_count] = b_group_max_exons
        else:
            self.group_max_exon_dict[self.group_count] = a_group_max_exons
        
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

        print("merge_a_b_groups")
        print(self.group_count)
        print(trans_a)
        print(trans_b)
        
        
#####################################################################

    def merge_a_b_groups_nocap_old(self,trans_a,trans_b,trans_obj_dict):
        
        self.group_count += 1
        print("invoke merge_a_b_groups_nocap " + str(self.group_count)+ " " + trans_a + " " +trans_b )
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

##############################################################################

    def merge_a_b_groups_nocap(self,trans_a,trans_b,trans_obj_dict):
        # need to create new group and cant use time saving thing in capped merge
        # this is because nocaps can have multiple groups
        self.group_count += 1


        print("invoke merge_a_b_groups_nocap " + str(self.group_count )+ " " + trans_a + " " +trans_b )
        #only cap lib trans should be used for merging groups
        if len(list(self.trans_group_dict[trans_a].keys())) > 1:
        #    if log_flag == "log_on":
            print("multiple groups a nocap")

        if len(list(self.trans_group_dict[trans_b].keys())) > 1:

            print("multiple groups b nocap")

        #check that trans b has same num exons as trans a
        trans_obj_a = trans_obj_dict[trans_a]
        trans_obj_b = trans_obj_dict[trans_b]

        if trans_obj_a.num_exons != trans_obj_b.num_exons:
            print("Error trans_a does not same num exons as trans_b")
            print(trans_a + " " + trans_b)
            print(str(trans_obj_a.num_exons) + " " + str(trans_obj_b.num_exons))
            sys.exit()

        a_group_num_list = list(self.trans_group_dict[trans_a].keys())

        if self.group_count in self.group_trans_dict:
            print("group num already used")
            sys.exit()

        #make new group
        self.group_trans_dict[self.group_count] = {}

        self.group_max_exon_dict[self.group_count] = 0


        for a_group_num in a_group_num_list:
            #print("group_max_exon_dict")
            #print(trans_a)
            #print(self.group_max_exon_dict)
            #print(a_group_num_list)
            #print(self.trans_group_dict)
            #print(self.group_trans_dict[a_group_num])
            if self.group_max_exon_dict[a_group_num] == trans_obj_dict[trans_a].num_exons: # if this trans has as many exons as longest trans in group
            #if self.group_trans_dict[a_group_num][trans_a].startswith("long"): # could be longest or long
                for group_trans in self.group_trans_dict[a_group_num]:
                    self.group_trans_dict[self.group_count][group_trans] = 1 #make new group
                    self.trans_group_dict[group_trans].pop(a_group_num, None) # remove old group num
                    self.trans_group_dict[group_trans][self.group_count] = 1 #add new group num

                self.group_max_exon_dict[self.group_count] = self.group_max_exon_dict[a_group_num]
                #remove old group
                self.group_trans_dict.pop(a_group_num, None)
                self.group_longest_dict.pop(a_group_num, None)

        # get b group list after a group processing because
        # some of b groups could have been same as a groups
        # thus they would have been deleted
        b_group_num_list = list(self.trans_group_dict[trans_b].keys())

        for b_group_num in b_group_num_list:
            if self.group_max_exon_dict[b_group_num] == trans_obj_dict[trans_b].num_exons: # if this trans has as many exons as longest trans in group
            #if self.group_trans_dict[b_group_num][trans_b].startswith("long"): # could be longest or long
                for group_trans in self.group_trans_dict[b_group_num]:
                    self.group_trans_dict[self.group_count][group_trans] = 1 # make new group
                    self.trans_group_dict[group_trans].pop(b_group_num, None) # remove old group num
                    self.trans_group_dict[group_trans][self.group_count] = 1 # add new group num

                if self.group_max_exon_dict[self.group_count] < self.group_max_exon_dict[b_group_num]:
                    self.group_max_exon_dict[self.group_count] = self.group_max_exon_dict[b_group_num]
                #remove old group
                self.group_trans_dict.pop(b_group_num, None)
                self.group_longest_dict.pop(b_group_num, None)

        #add longest trans information
        trans_id_list = list(self.group_trans_dict[self.group_count].keys())

        if len(trans_id_list) > 0:# new group was made
            #redo longest trans flags

            #get longest trans in new group
            longest_trans_id = trans_id_list[0]
            for check_trans_id in trans_id_list:
                longest_trans_flag = longest_transcript(check_trans_id,[longest_trans_id],trans_obj_dict)
                if longest_trans_flag == "longest":
                    longest_trans_id = check_trans_id

            #refresh trans flags in new group
            self.group_longest_dict[self.group_count] = {}
            self.group_longest_dict[self.group_count]["longest"] = {}

            for trans_c in trans_id_list:
                longest_trans_flag = longest_transcript(trans_c,[longest_trans_id],trans_obj_dict)
                self.trans_group_dict[trans_c][self.group_count] = longest_trans_flag
                self.group_trans_dict[self.group_count][trans_c] = longest_trans_flag

                if longest_trans_flag == "longest":
                    self.group_longest_dict[self.group_count]["longest"][trans_c] = 1
        else: # no new group was made, merge did not happen, 2 short transcripts
            self.group_trans_dict.pop(self.group_count, None)

        ###################################### <<< edge of work
#####################################################################
####################################################################################################

def hunter_prey_capped(capped_trans_obj_list,trans_obj_dict,transgroup):
    #new cluster grouping algorithm
    ############################################################################################################
    ###########################
    # convert capped_trans_obj_list to trans_obj_dict for while looping

    # clusters that have yet been grouped
    ungrouped_trans_obj_dict = {}  # ungrouped_trans_obj_dict[cluster_id] =  trans_obj
    # clusters that have been grouped
    grouped_trans_obj_dict = {}  # grouped_trans_obj_dict[cluster_id] =  trans_obj
    # grouped clusters that have been not been searched/been used as a hunter
    unsearched_trans_obj_dict = {}  # unsearched_trans_obj_dict[cluster_id] =  trans_obj

    for trans_obj in capped_trans_obj_list:
        ungrouped_trans_obj_dict[trans_obj.uniq_trans_id] = trans_obj


    ###########################
    
    all_trans_id_dict = {} # all_trans_id_dict[trans id] = 1

    ungrouped_count = len(capped_trans_obj_list)
    unsearched_count = 0

    while ungrouped_count > 0:

        if unsearched_count == 0:
            ungrouped_cluster_list = list(ungrouped_trans_obj_dict.keys())
            ungrouped_cluster_list.sort()

            # hunter trans and prey trans, hunter used to look for prey
            hunter_cluster_id = ungrouped_cluster_list[0]
            hunter_trans_obj = ungrouped_trans_obj_dict[hunter_cluster_id]

            hunter_fivecap_flag = source_dict[hunter_trans_obj.source_id]['seq_type']

            # remove hunter from ungrouped
            ungrouped_trans_obj_dict.pop(hunter_cluster_id)

            unsearched_count = 1

        while unsearched_count > 0 and ungrouped_count > 0:

            if hunter_cluster_id == "new_hunter":
                unsearched_cluster_list = list(unsearched_trans_obj_dict.keys())
                unsearched_cluster_list.sort()
                # hunter trans and prey trans, hunter used to look for prey
                hunter_cluster_id = unsearched_cluster_list[0]
                hunter_trans_obj = unsearched_trans_obj_dict[hunter_cluster_id]

                unsearched_trans_obj_dict.pop(hunter_cluster_id)

            all_trans_id_dict[hunter_cluster_id] = 1
            hunter_strand = hunter_trans_obj.strand
            hunter_num_exons = hunter_trans_obj.num_exons

            a_group_check = transgroup.check_trans_status(hunter_cluster_id)
            # make groups for each transcript if no group
            if a_group_check != 1:
                transgroup.new_group_a(hunter_cluster_id,trans_obj_dict)

            for prey_cluster_id in ungrouped_trans_obj_dict:
                prey_trans_obj = ungrouped_trans_obj_dict[prey_cluster_id]

                prey_fivecap_flag = source_dict[prey_trans_obj.source_id]['seq_type']

                prey_strand = prey_trans_obj.strand
                prey_num_exons = prey_trans_obj.num_exons

                # this condition should not happen anymore because i delete hunters after I use them
                if hunter_cluster_id == prey_cluster_id:
                    continue

                # check strand match
                if hunter_strand != prey_strand:
                    print("Strand of transcripts within gene do not match")
                    sys.exit()

                b_group_check = transgroup.check_trans_status(prey_cluster_id)
                # make groups for each transcript if no group
                if b_group_check != 1:
                    transgroup.new_group_a(prey_cluster_id,trans_obj_dict)

                group_match_flag = transgroup.check_same_group(hunter_cluster_id, prey_cluster_id)

                # this shoudn't be needed anymore due to new dict based group search
                if group_match_flag == 1:  # if they are both in the same group
                    continue

                # trans_comp_flag, start_match_list, start_diff_list, end_match_list, end_diff_list, short_trans, long_trans, min_exon_num, diff_num_exon_flag
                trans_comp_flag,start_match_list,start_diff_list,end_match_list,end_diff_list,short_trans,min_exon_num = compare_transcripts_both_capped(hunter_trans_obj, prey_trans_obj, hunter_fivecap_flag,prey_fivecap_flag, hunter_strand)

                trans_match_flag = 0
                # same_transcript means clusters should be grouped for collapsing!
                if trans_comp_flag == "same_transcript":
                    trans_match_flag = 1

                a_group_check = transgroup.check_trans_status(hunter_cluster_id)
                b_group_check = transgroup.check_trans_status(prey_cluster_id)

                ##########################################Affects all downstream code!
                if trans_match_flag != 1:  # skip if there is no match
                    continue

                else:  # if they are in different groups, merge groups, applies to both caps
                    transgroup.merge_a_b_groups(hunter_cluster_id, prey_cluster_id)

                    # remove the prey cluster from dict to avoid redundant searching
                    # ungrouped_trans_obj_dict.pop(prey_cluster_id) # remove this outside of for loop
                    # add prey to unsearched dict
                    unsearched_trans_obj_dict[prey_cluster_id] = prey_trans_obj

            # remove grouped prey from ungrouped dict
            for unsearched_cluster_id in unsearched_trans_obj_dict:
                if unsearched_cluster_id in ungrouped_trans_obj_dict:
                    ungrouped_trans_obj_dict.pop(unsearched_cluster_id)

            unsearched_count = len(unsearched_trans_obj_dict)
            # reset hunter id
            hunter_cluster_id = "new_hunter"

            ungrouped_count = len(ungrouped_trans_obj_dict)
    
    #trans_group_dict = transgroup.trans_group_dict
    #group_trans_dict = transgroup.group_trans_dict

  
    return transgroup
            
####################################################################################################

####################################################################################################

# compare capped groups with no_cap groups
def hunter_prey_nocap(nocap_trans_obj_list,trans_obj_dict,transgroup):
    # clusters that have yet been grouped
    ungrouped_trans_obj_dict = {}  # ungrouped_trans_obj_dict[cluster_id] =  trans_obj
    # clusters that have been grouped
    grouped_trans_obj_dict = {}  # grouped_trans_obj_dict[cluster_id] =  trans_obj
    # grouped clusters that have been not been searched/been used as a hunter
    unsearched_trans_obj_dict = {}  # unsearched_trans_obj_dict[cluster_id] =  trans_obj

    used_hunter_dict = {} # used_hunter_dict[hunter_id] = 1  # use this to skip hunters that have already been used
    
    # use this dict to organize clusters by num exons
    exon_trans_obj_dict = {} # exon_trans_obj_dict[num exons][cluster id] = trans_obj
    # use this to refresh sub exon dicts
    sub_exon_cluster_dict = {} # sub_exon_cluster_dict[num exons][cluster id] = 1

    print("running hunter_prey_nocap") ###################################################  troubleshooting
    print(str(len(nocap_trans_obj_list))) ###################################################  troubleshooting

    for trans_obj in nocap_trans_obj_list:        
        this_strand = trans_obj.strand
        this_num_exons = trans_obj.num_exons
        
        if this_num_exons not in exon_trans_obj_dict:
            exon_trans_obj_dict[this_num_exons] = {}
            sub_exon_cluster_dict[this_num_exons] = {}

        exon_trans_obj_dict[this_num_exons][trans_obj.uniq_trans_id] = trans_obj
        sub_exon_cluster_dict[this_num_exons][trans_obj.uniq_trans_id] = 1

    ###########################
    
    all_trans_id_dict = {} # all_trans_id_dict[trans id] = 1

    ungrouped_count = len(nocap_trans_obj_list)
    unsearched_count = 0
    
    ###########################

    exon_num_list = list(exon_trans_obj_dict.keys())
    exon_num_list.sort(reverse=True)
    
    #exon_num_level = exon_num_list[0]
    exon_num_index = 0
    
    sub_length_cluster_dict = {} # sub_length_cluster_dict[cluster id] = 1  use this to mark degraded clusters
    
    while exon_num_index < len(exon_num_list):
        
        exon_num_level = exon_num_list[exon_num_index]
        #initialize ungrouped for this exon num level
        ungrouped_trans_obj_dict = {}
        for ungroup_cluster_id in exon_trans_obj_dict[exon_num_level]:
            ungrouped_trans_obj_dict[ungroup_cluster_id] = exon_trans_obj_dict[exon_num_level][ungroup_cluster_id]

        ungrouped_count = len(ungrouped_trans_obj_dict)
    
        while ungrouped_count > 0:
    
            if unsearched_count == 0:
                ungrouped_cluster_list = list(ungrouped_trans_obj_dict.keys())
                ungrouped_cluster_list_sorted = []

                ungrouped_strand = ungrouped_trans_obj_dict[ungrouped_cluster_list[0]].strand
                # sort by 5' end length to deal with issue of long based nocap merging
                ungrouped_cluster_sort_dict = {} # ungrouped_cluster_sort_dict[start][cluster id] = 1
                for ungrouped_cluster_id in ungrouped_cluster_list:
                    if ungrouped_strand == "+":
                        ungrouped_five_prime = ungrouped_trans_obj_dict[ungrouped_cluster_id].start_pos  ####################################################
                    elif ungrouped_strand == "-":
                        ungrouped_five_prime = ungrouped_trans_obj_dict[ungrouped_cluster_id].end_pos
                    else:
                        print("Error with ungrouped_strand")
                        sys.exit()

                    if ungrouped_five_prime not in ungrouped_cluster_sort_dict:
                        ungrouped_cluster_sort_dict[ungrouped_five_prime] = {}

                    ungrouped_cluster_sort_dict[ungrouped_five_prime][ungrouped_cluster_id] = 1

                ungrouped_cluster_five_list = list(ungrouped_cluster_sort_dict.keys())
                if ungrouped_strand == "+":
                    ungrouped_cluster_five_list.sort()
                elif ungrouped_strand == "-":
                    ungrouped_cluster_five_list.sort(reverse=True)
                else:
                    print("Error with ungrouped_strand")
                    sys.exit()

                for ungrouped_five_coord in ungrouped_cluster_five_list:
                    for ungrouped_cluster_id in list(ungrouped_cluster_sort_dict[ungrouped_five_coord].keys()):
                        ungrouped_cluster_list_sorted.append(ungrouped_cluster_id)

                ungrouped_cluster_list = ungrouped_cluster_list_sorted
    
                # hunter trans and prey trans, hunter used to look for prey
                hunter_cluster_id = ungrouped_cluster_list[0]
                hunter_trans_obj = ungrouped_trans_obj_dict[hunter_cluster_id]
                # remove hunter from ungrouped
                ungrouped_trans_obj_dict.pop(hunter_cluster_id)

                used_hunter_dict[hunter_cluster_id] = 1

                print("hunter id unsearched_count: " + str(hunter_cluster_id)) ###################################################  troubleshooting
                print(ungrouped_cluster_list) ###################################################  troubleshooting
    
                unsearched_count = 1
                  
            #use this to keep track of sub clusters that have not already been grouped with this group
            this_sub_exon_cluster_dict = sub_exon_cluster_dict 
    
            while unsearched_count > 0:
                
                if hunter_cluster_id == "new_hunter":
                    unsearched_cluster_list = list(unsearched_trans_obj_dict.keys())
                    unsearched_cluster_list.sort()
                    
                    if len(unsearched_cluster_list) == 0:
                        print("empty unsearched_cluster_list")
                        print(hunter_cluster_id)
                        print(ungrouped_trans_obj_dict)
                        print(unsearched_count)
                        print(ungrouped_count)

                    # hunter trans and prey trans, hunter used to look for prey
                    hunter_cluster_id = unsearched_cluster_list[0]
                    hunter_trans_obj = unsearched_trans_obj_dict[hunter_cluster_id]

                    used_hunter_dict[hunter_cluster_id] = 1

                    # ungrouped_trans_obj_dict.pop(hunter_cluster_id)
                    unsearched_trans_obj_dict.pop(hunter_cluster_id)

                    this_ungrouped_trans_obj_list = list(ungrouped_trans_obj_dict.keys())

                    print("hunter id new_hunter: " + str(hunter_cluster_id)) ###################################################  troubleshooting
                    print(unsearched_cluster_list) ###################################################  troubleshooting
                    print(this_ungrouped_trans_obj_list) ###################################################  troubleshooting
                    print(ungrouped_count) ###################################################  troubleshooting



                # print("hunter: " + hunter_cluster_id )
                #print(exon_num_list)
    
                all_trans_id_dict[hunter_cluster_id] = 1
                hunter_strand = hunter_trans_obj.strand
                hunter_num_exons = hunter_trans_obj.num_exons
                hunter_fivecap_flag = source_dict[hunter_trans_obj.source_id]['seq_type']
    
                a_group_check = transgroup.check_trans_status(hunter_cluster_id)
                # make groups for each transcript if no group
                if a_group_check != 1:
                    transgroup.new_group_a(hunter_cluster_id,trans_obj_dict)
    
                ungrouped_trans_obj_list = list(ungrouped_trans_obj_dict.keys())
                unsearched_cluster_list = list(unsearched_trans_obj_dict.keys())
    
                # search at same exon num level
                ###############################################################################
                if hunter_cluster_id not in sub_length_cluster_dict: # only search at same level if this has not been grouped with longer transcript
                    for prey_cluster_id in ungrouped_trans_obj_dict:
        
                        prey_trans_obj = ungrouped_trans_obj_dict[prey_cluster_id]

                        prey_fivecap_flag = source_dict[prey_trans_obj.source_id]['seq_type']
    
                        prey_strand = prey_trans_obj.strand
                        prey_num_exons = prey_trans_obj.num_exons
        
                        # this condition should not happen anymore because i delete hunters after I use them
                        if hunter_cluster_id == prey_cluster_id:
                            continue
        
                        # check strand match
                        if hunter_strand != prey_strand:
                            print("Strand of transcripts within gene do not match")
                            sys.exit()
        
                        b_group_check = transgroup.check_trans_status(prey_cluster_id)
                        # make groups for each transcript if no group
                        if b_group_check != 1:
                            transgroup.new_group_a(prey_cluster_id,trans_obj_dict)
        
                        group_match_flag = transgroup.check_same_group(hunter_cluster_id, prey_cluster_id)
        
                        # this shoudn't be needed anymore due to new dict based group search
                        if group_match_flag == 1:  # if they are both in the same group
                            continue
        
                        # trans_comp_flag, start_match_list, start_diff_list, end_match_list, end_diff_list, short_trans, long_trans, min_exon_num, diff_num_exon_flag
                        trans_comp_flag, start_match_list, start_diff_list, end_match_list, end_diff_list, short_trans, min_exon_num = compare_transcripts_both_nocap(hunter_trans_obj, prey_trans_obj, hunter_fivecap_flag,prey_fivecap_flag, hunter_strand)
                        
                        # print("compare_transcripts: "+hunter_trans_obj.cluster_id +"\t"+prey_trans_obj.cluster_id+"\t" +trans_comp_flag + "\t" + str(hunter_trans_obj.num_exons)+ "\t" + str(prey_trans_obj.num_exons) + "\t" + str(diff_num_exon_flag))

                        #For nocap only!!!!
                        trans_match_flag = 0
                        if trans_comp_flag == "same_transcript":
                            trans_match_flag = 1
                        elif trans_comp_flag == "same_three_prime_same_exons" :
                            trans_match_flag = 1
                        #elif trans_comp_flag == "same_three_prime_diff_exons":
                        #    trans_match_flag = 1
        
                        a_group_check = transgroup.check_trans_status(hunter_cluster_id)
                        b_group_check = transgroup.check_trans_status(prey_cluster_id)
        
                        ##########################################Affects all downstream code!
                        ###For no cap!!!
                        # only merge groups if they have the same number of exons
                        # if diff num exons then only add shorter one to longer one
                        if trans_match_flag != 1: # skip if there is no match
                            continue
        
                        else: # if they are in different groups, but match
                            if hunter_num_exons == prey_num_exons:   # same number of exons

                                if trans_comp_flag == "same_transcript": ###################### 2019/06/07
                                    #print(transgroup.group_trans_dict)
                                    transgroup.merge_a_b_groups_nocap(hunter_cluster_id,prey_cluster_id,trans_obj_dict)
                                elif trans_comp_flag == "same_three_prime_same_exons" : ###################### 2019/06/07
                                    transgroup.add_a_to_b_group_nocap(prey_cluster_id,hunter_cluster_id)
                                else:
                                    print("Error with match flag")
                                    print(trans_comp_flag)
                                    sys.exit()

                                # make sure to only add unsearched prey if it has not already been used as a hunter
                                if prey_cluster_id not in used_hunter_dict:
                                    unsearched_trans_obj_dict[prey_cluster_id] = prey_trans_obj

                # search at same exon num level end
                ###############################################################################

                ################################################################################
                # search at lower exon num level
                prey_exon_index = exon_num_index
                
                #this_sub_exon_cluster_dict.pop(exon_num_list[prey_exon_index])

                while prey_exon_index < len(exon_num_list):
                    
                    prey_exon_index += 1
                    
                    if prey_exon_index >= len(exon_num_list):
                        continue
                    
                    prey_exon_level = exon_num_list[prey_exon_index]
                    
                    #initialize ungrouped for this exon num level
                    subgrouped_trans_obj_dict = {}
                    #for subgroup_cluster_id in exon_trans_obj_dict[prey_exon_level]:
                    for subgroup_cluster_id in this_sub_exon_cluster_dict[prey_exon_level]:
  
                        subgrouped_trans_obj_dict[subgroup_cluster_id] = exon_trans_obj_dict[prey_exon_level][subgroup_cluster_id]
                    
                    for prey_cluster_id in subgrouped_trans_obj_dict:
        
                        prey_trans_obj = subgrouped_trans_obj_dict[prey_cluster_id]

                        prey_fivecap_flag = source_dict[prey_trans_obj.source_id]['seq_type']

                        prey_strand = prey_trans_obj.strand
                        prey_num_exons = prey_trans_obj.num_exons
        
                        # this condition should not happen anymore because i delete hunters after I use them
                        if hunter_cluster_id == prey_cluster_id:
                            continue
        
                        # check strand match
                        if hunter_strand != prey_strand:
                            print("Strand of transcripts within gene do not match")
                            sys.exit()
        
                        b_group_check = transgroup.check_trans_status(prey_cluster_id)
                        # make groups for each transcript if no group
                        if b_group_check != 1:
                            transgroup.new_group_a(prey_cluster_id,trans_obj_dict)
        
                        group_match_flag = transgroup.check_same_group(hunter_cluster_id, prey_cluster_id)
        
                        # this shoudn't be needed anymore due to new dict based group search
                        if group_match_flag == 1:  # if they are both in the same group
                            continue
        
                        # trans_comp_flag, start_match_list, start_diff_list, end_match_list, end_diff_list, short_trans, long_trans, min_exon_num, diff_num_exon_flag \
                        trans_comp_flag, start_match_list, start_diff_list, end_match_list, end_diff_list, short_trans, min_exon_num = compare_transcripts_both_nocap(hunter_trans_obj, prey_trans_obj, hunter_fivecap_flag,prey_fivecap_flag , hunter_strand)
        
                        #For nocap only!!!!
                        trans_match_flag = 0
                        if trans_comp_flag == "same_transcript":
                            #trans_match_flag = 1
                            print("Error with subgroup search same_transcript")
                            sys.exit()
                        elif trans_comp_flag == "same_three_prime_same_exons" :
                            #trans_match_flag = 1
                            print("Error with subgroup search same_three_prime_same_exons")
                            sys.exit()
                        elif trans_comp_flag == "same_three_prime_diff_exons":
                            trans_match_flag = 1
        
                        a_group_check = transgroup.check_trans_status(hunter_cluster_id)
                        b_group_check = transgroup.check_trans_status(prey_cluster_id)
        
                        ##########################################Affects all downstream code!
                        ###For no cap!!!
                        # only merge groups if they have the same number of exons
                        # if diff num exons then only add shorter one to longer one
                        if trans_match_flag != 1: # skip if there is no match
                            continue
        
                        else: # if they are in different groups, but match
                            if hunter_num_exons == prey_num_exons:   # same number of exons
                                #transgroup.merge_a_b_groups_nocap(hunter_cluster_id,prey_cluster_id)
                                print("Error with subgroup exon match: equal")
                                sys.exit()
                            elif hunter_num_exons > prey_num_exons: #add shorter to longer
                                transgroup.add_a_to_b_group_nocap(prey_cluster_id,hunter_cluster_id)
                                sub_length_cluster_dict[prey_cluster_id] = 1
                                print("sub_length_cluster_dict : hunter - " +  hunter_cluster_id + "; prey - " + prey_cluster_id)
                                
                                #this_sub_exon_cluster_dict[prey_exon_level].pop(prey_cluster_id)  ######### 2019/06/06
                                
                            elif prey_num_exons > hunter_num_exons: # add shorter to longer
                                #transgroup.add_a_to_b_group(hunter_cluster_id,prey_cluster_id)
                                print("Error with subgroup exon match: hunter smaller")
                                sys.exit()

                # reset hunter id
                # sub_length_cluster_dict[hunter_cluster_id] = 1
                hunter_cluster_id = "new_hunter"
                unsearched_count = len(unsearched_trans_obj_dict)
                ungrouped_count = len(ungrouped_trans_obj_dict)


    
        exon_num_index += 1
                ################################################################################



    ############################################################################################################
    # End of new cluster grouping algorithm

    #trans_group_dict = transgroup.trans_group_dict
    #group_trans_dict = transgroup.group_trans_dict
    
    return transgroup


####################################################################################################

# compare capped groups with no_cap groups
# should run hunter_prey_capped before running this to get capped groups ready
def hunter_prey_mixed(capped_trans_obj_list,nocap_trans_obj_list,trans_obj_dict,transgroup):
    #new cluster grouping algorithm
    ############################################################################################################
    ###########################
    # convert capped_trans_obj_list to trans_obj_dict for while looping

    # clusters that have yet been grouped
    ungrouped_trans_obj_dict = {}  # ungrouped_trans_obj_dict[uniq_trans_id] =  trans_obj
    # clusters that have been grouped
    grouped_trans_obj_dict = {}  # grouped_trans_obj_dict[uniq_trans_id] =  trans_obj
    # grouped clusters that have been not been searched/been used as a hunter
    unsearched_trans_obj_dict = {}  # unsearched_trans_obj_dict[uniq_trans_id] =  trans_obj
    
    # use this dict to organize clusters by num exons
    exon_trans_obj_dict = {} # exon_trans_obj_dict[num exons][cluster id] = trans_obj
    # use this to refresh sub exon dicts
    sub_exon_cluster_dict = {} # sub_exon_cluster_dict[num exons][cluster id] = 1

    for trans_obj in nocap_trans_obj_list:        
        this_strand = trans_obj.strand
        this_num_exons = trans_obj.num_exons
        
        if this_num_exons not in exon_trans_obj_dict:
            exon_trans_obj_dict[this_num_exons] = {}
            sub_exon_cluster_dict[this_num_exons] = {}

        exon_trans_obj_dict[this_num_exons][trans_obj.uniq_trans_id] = trans_obj
        sub_exon_cluster_dict[this_num_exons][trans_obj.uniq_trans_id] = 1


    for hunter_trans_obj in capped_trans_obj_list:
        hunter_id = hunter_trans_obj.uniq_trans_id
        hunter_fivecap_flag = source_dict[hunter_trans_obj.source_id]['seq_type']
        hunter_strand = hunter_trans_obj.strand
        hunter_num_exons = hunter_trans_obj.num_exons

        
        for prey_trans_obj in nocap_trans_obj_list:
            prey_id = prey_trans_obj.uniq_trans_id
            prey_fivecap_flag = source_dict[prey_trans_obj.source_id]['seq_type']
            prey_num_exons = prey_trans_obj.num_exons
            
            #check if prey in hunter group already
            group_match_flag = transgroup.check_same_group(hunter_id,prey_id) ################################################
            
            if group_match_flag == 1: # skip this if hunter and prey are already associated
                continue
        
            trans_comp_flag,start_match_list,start_diff_list,end_match_list,end_diff_list,short_trans,min_exon_num = compare_transcripts_capped_nocap(hunter_trans_obj,prey_trans_obj,hunter_fivecap_flag,prey_fivecap_flag,hunter_strand)
        
            
            #For nocap only!!!!
            trans_match_flag = 0
            if trans_comp_flag == "same_transcript":
                trans_match_flag = 1
            elif trans_comp_flag == "same_three_prime_same_exons" :
                trans_match_flag = 1
            elif trans_comp_flag == "same_three_prime_diff_exons":
                trans_match_flag = 1


            ##########################################Affects all downstream code!
            ###For no cap!!!
            # only merge groups if they have the same number of exons
            # if diff num exons then only add shorter one to longer one
            if trans_match_flag != 1: # skip if there is no match
                continue

            else: # if they are in different groups, but match
                if hunter_num_exons >= prey_num_exons:   # same number of exons or fewer for prey

                    transgroup.add_a_to_b_group_nocap(prey_id,hunter_id)

    ############################################################################################################
    # End of new cluster grouping algorithm
    #trans_group_dict = transgroup.trans_group_dict
    #group_trans_dict = transgroup.group_trans_dict
    

    return transgroup


####################################################################################################
def simplify_gene(trans_obj_list,trans_obj_dict): # goes through transcripts in gene and groups transcripts for collapsing

    transgroup = TransGroup("5cap_mix")
    
    #keep track of all no cap trans
    nocap_trans_dict = {} # nocap_trans_dict[trans_id] = 1
    capped_trans_dict = {} # capped_trans_dict[trans_id] = 1
    
    nocap_trans_obj_list = []
    capped_trans_obj_list = []
    
    #### Collect transcripts into groups of capped and no cap
    ############################################
    
    for i in xrange(len(trans_obj_list)):
        trans_obj = trans_obj_list[i]
        uniq_trans_id = trans_obj.uniq_trans_id
        strand = trans_obj.strand
        
        #assume source_dict[trans_obj.source_id]['seq_type'] == no_cap or capped
        fivecap_flag = source_dict[trans_obj.source_id]['seq_type']
        
        #collect all nocap transcripts
        if fivecap_flag == "no_cap":
            nocap_trans_dict[uniq_trans_id] = 1
            nocap_trans_obj_list.append(trans_obj)
        elif fivecap_flag == "capped":
            capped_trans_dict[uniq_trans_id] = 1
            capped_trans_obj_list.append(trans_obj)
        
        a_group_check = transgroup.check_trans_status(uniq_trans_id)
        #make groups for each transcript if no group
        if a_group_check != 1:
            transgroup.new_group_a(uniq_trans_id,trans_obj_dict)
    
    ############################################
    # go through capped transcripts first
    
    if len(capped_trans_obj_list) > 0 :
        transgroup = hunter_prey_capped(capped_trans_obj_list,trans_obj_dict,transgroup)
    ############################################
    # compare capped groups with no-cap transcripts

    if len(capped_trans_obj_list) > 0 and len(nocap_trans_obj_list) > 0 :
        transgroup = hunter_prey_mixed(capped_trans_obj_list,nocap_trans_obj_list,trans_obj_dict,transgroup)
    
    #############################################
    # go through nocap transcripts that were not grouped
    
    ungrouped_nocap_trans_obj_list = []
    
    for nocap_trans_obj in nocap_trans_obj_list:
        nocap_trans_id = nocap_trans_obj.uniq_trans_id
        
        single_group_check = transgroup.check_trans_group(nocap_trans_id)
        
        if single_group_check == 1:
            ungrouped_nocap_trans_obj_list.append(nocap_trans_obj)

    transgroup = hunter_prey_nocap(ungrouped_nocap_trans_obj_list,trans_obj_dict,transgroup)


    #############################################
 
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

def format_gene_report_line(trans_obj_list,total_gene_count,total_final_trans_count,track_gene_source):
    
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


    ###########################
    # add gene source info 2020/05/31

    all_track_source_dict = {} # all_track_source_dict[source] = count
    track_gene_source_split = track_gene_source.split(",")

    source_order_list = []

    for track_source_gene in track_gene_source_split:
        track_source = track_source_gene.split("_")[0]
        if track_source not in all_track_source_dict:
            all_track_source_dict[track_source] = 0
            source_order_list.append(track_source)

        all_track_source_dict[track_source] += 1

    track_source_list = []
    for track_source in source_order_list:
        source_gene_count = all_track_source_dict[track_source]

        this_track_source_line = track_source + ":" + str(source_gene_count)

        track_source_list.append(this_track_source_line)

    all_track_source_line = ",".join(track_source_list)

    gene_report_list.append(str(track_gene_source))
    gene_report_list.append(str(all_track_source_line))

    # add gene source info 2020/05/31
    ###########################
    
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
        
        if len(id_line_split) < 2:
            print("Error with bed file ID field")
            print(trans_line_list)
            print("bed12 files must have the gene ID's and transcript ID's formatted as such \"gene_id;transcript_id\" in the 4th column.")
            print("The gene ID must be the first subfield and the subfields must be delimited with a semicolon (;).")
            sys.exit()
            
        
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
            issue_trans_id = trans_obj.trans_id
            issue_strand = trans_obj.strand
            issue_source = trans_obj.source_id
            this_error_line = "Source: " + issue_source + " Trans: " + issue_trans_id + " Strand: " + issue_strand
            print(this_error_line)
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

            # keep track of gene sources
            track_gene_source_dict = {} # track_gene_source_dict[source_gene_id] = 1

            trans_count = 0
            for merged_obj in sorted_merge_obj_list:
                trans_count += 1
                final_trans_id = "G" + str(total_gene_count) + "." + str(trans_count)
                merged_obj.trans_id = final_trans_id
                
                #merged_obj_dict[final_trans_id] = merged_obj

                #write out to transcript report file
                trans_report_line = merged_obj.format_trans_report_line()
                outfile_trans_report.write(trans_report_line)
                outfile_trans_report.write("\n")

                #write out to transcript merge report file
                # check for merge source id in case of using source ID's for merge ID's
                for merged_trans_id in merged_obj.merged_trans_dict:

                    #########################
                    # use this for keeping track of source genes for gene report
                    track_merge_gene_id_split = merged_trans_id.split(".")
                    if len(track_merge_gene_id_split) > 2:
                        track_merge_gene_id_list = []

                        for t_index in xrange(len(track_merge_gene_id_split)-1):
                            track_merge_gene_id_list.append(track_merge_gene_id_split[t_index])

                        track_merge_gene_id = ".".join(track_merge_gene_id_list)
                    else:
                        track_merge_gene_id = track_merge_gene_id_split[0]

                    track_gene_source_dict[track_merge_gene_id] = 1
                    # use this for keeping track of source genes for gene report
                    #########################

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

                    ########
                    # get source id in case of source id flag
                    if source_id_flag != "no_source_id":
                        this_source_name_split = merged_trans_id.split("_")
                        this_source_name = this_source_name_split[0]

                        # this accounts for underscores used in the original trans id's
                        this_source_name_split.pop(0)
                        this_source_trans_id = "_".join(this_source_name_split)

                        for this_source_id_flag in source_id_flag_list:
                            if this_source_id_flag == this_source_name:
                                this_source_gene_id = source_trans_gene_dict[this_source_name][this_source_trans_id]##################################################

                                if merged_obj.extra_trans_id_list == "NA":
                                    merged_obj.extra_trans_id_list = []
                                    merged_obj.extra_gene_id_list = []

                                merged_obj.extra_trans_id_list.append(this_source_trans_id)
                                merged_obj.extra_gene_id_list.append(this_source_gene_id)

                    ### add CDS information from a source file
                    if source_cds_flag !=  "no_cds": ##################################################AXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                        this_source_name_split = merged_trans_id.split("_")
                        this_source_name = this_source_name_split[0]

                        # this accounts for underscores used in the original trans id's
                        this_source_name_split.pop(0)
                        this_source_trans_id = "_".join(this_source_name_split)

                        for this_source_cds_flag in source_cds_flag_list:
                            if this_source_cds_flag == this_source_name:
                                merged_obj.start_cds = merged_trans_obj.cds_start
                                merged_obj.end_cds = merged_trans_obj.cds_end


                # write out to bed file
                bed_line = merged_obj.format_bed_line()
                outfile_bed.write(bed_line)
                outfile_bed.write("\n")

            track_gene_source_list = list(track_gene_source_dict.keys())
            track_gene_source_list.sort()
            track_gene_source = ",".join(track_gene_source_list)

            gene_report_line = format_gene_report_line(trans_obj_list,total_gene_count, trans_count,track_gene_source) #######################################
            outfile_gene_report.write(gene_report_line)
            outfile_gene_report.write("\n")
                           
    return total_gene_count
                

        
        
        ##############################################convert into trans object !!!!!!!!!!!!!!!!
    
#end of def
###################################################################################################
###################################################################################################
###################################################################################################
source_trans_gene_dict = {} # source_trans_gene_dict[source name][trans_id] = gene_id

for file_line in filelist_file_contents:
    file_line_split = file_line.split("\t")


    # check for dos ^M new line character
    if "\r" in file_line:
        print("Error with " + filelist_file)
        print(file_line)
        print("Please make sure the filelist file does not have DOS new line characters")
        sys.exit()

    # check for issues with the use of spaces in filelist file
    space_line_split = file_line.split(" ")
    if len(space_line_split) > 1:
        print("Error with " + filelist_file)
        print(file_line)
        print("Please make sure it is tab separated with no empty lines and no spaces")
        sys.exit()

    if len(file_line_split) != 4:
        print("Error with " + filelist_file )
        print(file_line)
        print("Please make sure it is tab separated with no empty lines")
        sys.exit()
    
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
            print("Incorrect seq types given in filelist file. Please use capped or  no_cap")
            sys.exit()

        source_trans_gene_dict[source_id] = {}
    
    print("opening bed list")
    bed_file_contents = open(filename).read().rstrip("\n").split("\n")

    ###################
    #2020/05/31
    # check for sorted bed file

    check_scaffold = "none"
    check_prev_start = 0

    #for line in bed_file_contents:
    #    line_split = line.split("\t")

    #    scaffold = line_split[0]
    #    trans_start = int(line_split[1])
    #    trans_end = int(line_split[2])
    #    id_line = line_split[3]

    #    if check_scaffold ==  "none":
    #        check_scaffold = scaffold

    #    if trans_start < check_prev_start:
    #        if check_scaffold == scaffold:
    #            print("bed lines out of order")
    #            print(id_line)
    #            sys.exit()

    #    check_prev_start = trans_start

    #    if scaffold != check_scaffold:
    #        check_scaffold = scaffold
    #        check_prev_start = 0


    # 2020/05/31
    ###################

    
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
        
        line_split.append(source_id)
        
        if scaffold not in bed_dict:
            bed_dict[scaffold] = {}
        if trans_start not in bed_dict[scaffold]:
            bed_dict[scaffold][trans_start] = {}
        if trans_end not in bed_dict[scaffold][trans_start]:
            bed_dict[scaffold][trans_start][trans_end] = []
        
        bed_dict[scaffold][trans_start][trans_end].append(line_split)

        this_gene_id = id_line_split[0]
        this_trans_id = id_line_split[1]

        source_trans_gene_dict[source_id][this_trans_id] = this_gene_id
        
scaff_list = list(bed_dict.keys())
scaff_list.sort()

# if source id flag used, check that source given is in filelist
if source_id_flag != "no_source_id":
    for this_source_id_flag in source_id_flag_list:
        if this_source_id_flag not in source_dict:
            print("Error: Source name given in -n parameter is not found in filelist file.")
            print("Terminating early.")
            sys.exit()

#####################################################
colour_source_list = list(source_dict.keys())
colour_source_list.sort()
colour_source_count = 0
source_colour_dict = {} # source_colour_dict[source] = colour number
for colour_source in colour_source_list:
    colour_source_count += 1
    if colour_source_count > 10:
        colour_source_count = 10
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
        
        if trans_start < last_trans_end: # matches previous group

            for trans_end in trans_end_list:
                # add transcripts to group
                group_trans_list.extend(bed_dict[scaffold][trans_start][trans_end])
                if last_trans_end < trans_end: # update last trans end
                    last_trans_end = trans_end
        
        elif trans_start >= last_trans_end: # if this is a new group, process past group
            
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
            
            
        
            
    

print("TAMA Merge has completed successfully!")
    
