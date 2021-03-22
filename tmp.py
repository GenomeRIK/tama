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
TAMA Collapse

Author: Richard I. Kuo

This script collapses transcripts and groups transcripts into genes for long reads mapped onto a genome assembly.

"""

tc_version = 'tc0.0'
tc_date = 'tc_version_date_2020_12_14'

### Notes on changes
# Fixed issue with coordinates of soft clipped variants in the variant output file. 

#####################################################################
#####################################################################

ap = argparse.ArgumentParser(description='This script collapses mapped transcript models')

ap.add_argument('-s', type=str, nargs=1, help='Sorted sam file (required)')
ap.add_argument('-f', type=str, nargs=1, help='Genome fasta file (required)')
ap.add_argument('-p', type=str, nargs=1, help='Output prefix (required)')
ap.add_argument('-x', type=str, nargs=1, help='Capped flag: capped or no_cap')
ap.add_argument('-e', type=str, nargs=1, help='Collapse exon ends flag: common_ends or longest_ends (default common_ends)')
ap.add_argument('-c', type=str, nargs=1, help='Coverage (default 99)')
ap.add_argument('-i', type=str, nargs=1, help='Identity (default 85)')
ap.add_argument('-icm', type=str, nargs=1, help='Identity calculation method (default ident_cov for including coverage) (alternate is ident_map for excluding hard and soft clipping)')
ap.add_argument('-a', type=str, nargs=1, help='5 prime threshold (default 10)')
ap.add_argument('-m', type=str, nargs=1, help='Exon/Splice junction threshold (default 10)')
ap.add_argument('-z', type=str, nargs=1, help='3 prime threshold (default 10)')
ap.add_argument('-d', type=str, nargs=1, help='Flag for merging duplicate transcript groups (default is merge_dup will merge duplicates ,no_merge quits when duplicates are found)')
ap.add_argument('-sj', type=str, nargs=1, help='Use error threshold to prioritize the use of splice junction information from collapsing transcripts(default no_priority, activate with sj_priority)')
ap.add_argument('-sjt', type=str, nargs=1, help='Threshold for detecting errors near splice junctions (default is 10bp)')
ap.add_argument('-lde', type=str, nargs=1, help='Threshold for amount of local density error near splice junctions that is allowed (default is 1000 errors which practically means no threshold is applied)')
ap.add_argument('-ses', type=str, nargs=1, help='Simple error symbol. Use this to pick the symbol used to represent matches in the simple error string for LDE output.')

ap.add_argument('-b', type=str, nargs=1, help='Use BAM instead of SAM')
ap.add_argument('-log', type=str, nargs=1, help='Turns off log output to screen of collapsing process. (default on, use log_off to turn off)')
ap.add_argument('-v', type=str, nargs=1, help='Prints out version date and exits.')

ap.add_argument('-rm', type=str, nargs=1, help='Run mode allows you to use original or low_mem mode, default is original')
ap.add_argument('-vc', type=str, nargs=1, help='Variation covwerage threshold: Default 5 reads')


opts = ap.parse_args()

#check for version request
if not opts.v:
    print(tc_date)
else:
    print(tc_date)
    print("Program did not run")
    sys.exit()

#check for missing args
missing_arg_flag = 0
if not opts.s:
    print("Sam file is missing")
    missing_arg_flag = 1
if not opts.f:
    print("Fasta file missing")
    missing_arg_flag = 1
if not opts.p:
    print("Output prefix name missing")
    missing_arg_flag = 1
    
if not opts.x:
    print("Default capped flag will be used: Capped")
    fiveprime_cap_flag = "capped"
else:
    fiveprime_cap_flag = opts.x[0]
    
    if fiveprime_cap_flag != "no_cap" and fiveprime_cap_flag != "capped":
        print("Error with cap flag. Should be capped or no_cap.")
        sys.exit()

if not opts.e:
    print("Default collapse exon ends flag will be used: common_ends")
    collapse_flag = "common_ends"
else:
    collapse_flag = opts.e[0]

if not opts.c:
    print("Default coverage: 99")
    coverage_threshold = 99.0
else:
    coverage_threshold = float(opts.c[0])
    
if not opts.i:
    print("Default identity: 85")
    identity_threshold = 85.0
else:
    identity_threshold = float(opts.i[0])
    
if not opts.icm:
    print("Default identity calculation method: ident_cov")
    ident_calc_method = 'ident_cov'
else:
    ident_calc_method = str(opts.icm[0])
    if ident_calc_method != "ident_cov" and ident_calc_method != "ident_map":
        print("Error with -icm input. Should be ident_cov or ident_map. Run terminated.")
        print(ident_calc_method)
        sys.exit()
    
if not opts.a:
    print("Default 5 prime threshold: 10")
    fiveprime_threshold  = 10
else:
    fiveprime_threshold  = int(opts.a[0])
    
if not opts.m:
    print("Default exon/splice junction threshold: 10")
    exon_diff_threshold = 10
else:
    exon_diff_threshold = int(opts.m[0])

if not opts.z:
    print("Default 3 prime threshold: 10")
    threeprime_threshold = 10
else:
    threeprime_threshold = int(opts.z[0])
    
if not opts.d:
    print("Default duplicate merge flag: merge_dup")
    duplicate_flag = "merge_dup"
else:
    duplicate_flag = str(opts.d[0])

if not opts.sj:
    print("Default splice junction priority: no_priority")
    sj_priority_flag = "no_priority"
else:
    sj_priority_flag = str(opts.sj[0])

if not opts.sjt:
    print("Default splice junction error threshold: 10")
    sj_err_threshold = 10
else:
    sj_err_threshold = int(opts.sjt[0])
    
if not opts.lde:
    print("Default splice junction local density error threshold: 1000")
    lde_threshold = 1000
else:
    lde_threshold = int(opts.lde[0])

if not opts.ses:
    print("Default simple error symbol for matches is the underscore \"_\" .")
    ses_match_char = "_"
else:
    ses_match_char = int(opts.ses[0])

if not opts.b:
    print("Using SAM format for reading in.")
    bam_flag = "SAM"
else:
    print("Using BAM format for reading in.")
    import pysam
    bam_flag = str(opts.b[0])

if not opts.log:
    print("Default log output on")
    log_flag = "log_on"
else:
    log_flag = str(opts.log[0])
    if log_flag != "log_off":
        print("Please use log_off to turn off log prints to screen")
        sys.exit()

if not opts.rm:
    print("Default run mode original")
    run_mode_flag = "original"
else:
    run_mode_flag = str(opts.rm[0])
    if run_mode_flag != "original" and run_mode_flag != "low_mem" :
        print("Please use original or low_mem for -rm setting")
        sys.exit()

if not opts.vc:
    print("Default 5 read threshold")
    var_support_threshold = 5
else:
    var_support_threshold = int(opts.vc[0])




if missing_arg_flag == 1:
    print("Please try again with complete arguments")

sam_file = opts.s[0]
fasta_file_name = opts.f[0]
outfile_prefix = opts.p[0]

input_sambam_flag = "na"

if sam_file.endswith("bam"):
    input_sambam_flag = "BAM"
    if bam_flag == "SAM":
        print("You designated SAM input but are supplying BAM. Please use the same format for input as specified.")
        sys.exit()

if sam_file.endswith("sam"):
    input_sambam_flag = "SAM"
    if bam_flag == "BAM":
        print("You designated BAM input but are supplying SAM. Please use the same format for input as specified.")
        sys.exit()

if input_sambam_flag == "na":
    print("Input SAM/BAM file not recoginized from extension format designation.")


#####################################################################
#####################################################################

start_time = time.time()
prev_time = start_time

#print("opening sam file")
#sam_file = sys.argv[1]
#sam_file_contents = open(sam_file).read().rstrip("\n").split("\n")

#print("opening fasta file")
#fasta_file_name = sys.argv[2]

#outfile_prefix = sys.argv[3]

#fiveprime_cap_flag = "capped"
#collapse_flag = "common_ends"

#default threshold, will add option to change via arguments
#coverage_threshold = 99.0
#identity_threshold = 85.0

#default poly A threshold
a_window = 20
a_perc_thresh = 70.0

no_mismatch_flag = "0" # use this for showing no mismatch near splice junction
# see calc_error_rate and  sj_error_priority_start and sj_error_priority_end

bed_outfile_name = outfile_prefix + ".bed"
outfile_bed = open(bed_outfile_name,"w")

cluster_outfile_name = outfile_prefix + "_read.txt"
outfile_cluster = open(cluster_outfile_name,"w")
cluster_line = "\t".join(["read_id","mapped_flag","accept_flag","percent_coverage","percent_identity","error_line<h;s;i;d;m>", "length", "cigar"])
outfile_cluster.write(cluster_line)
outfile_cluster.write("\n")

trans_report_outfile_name = outfile_prefix + "_trans_report.txt"
outfile_trans_report = open(trans_report_outfile_name,"w")
trans_report_line = "\t".join(["transcript_id","num_clusters","high_coverage","low_coverage","high_quality_percent","low_quality_percent","start_wobble_list","end_wobble_list","collapse_sj_start_err","collapse_sj_end_err","collapse_error_nuc"])
#trans_report_line = "\t".join(["transcript_id","num_clusters","high_coverage","low_coverage","high_quality_percent","low_quality_percent","start_wobble_list","end_wobble_list","collapse_sj_start_err","collapse_sj_end_err","collapse_error_nuc","sj_error_simple"])
outfile_trans_report.write(trans_report_line)
outfile_trans_report.write("\n")

trans_clust_outfile_name = outfile_prefix + "_trans_read.bed"
outfile_trans_clust_report = open(trans_clust_outfile_name,"w")
#trans_clust_line = "\t".join(["transcript_id","cluster_id","scaffold","strand","start","end","exon_starts","exon_ends"])
#outfile_trans_clust_report.write(trans_clust_line)
#outfile_trans_clust_report.write("\n")

if run_mode_flag == "original":
    variant_outfile_name = outfile_prefix + "_variants.txt"
    outfile_variant = open(variant_outfile_name,"w")
    variant_file_line = "\t".join(["scaffold","position","type","ref_allele","alt_allele","count","cov_count","cluster_list"])
    outfile_variant.write(variant_file_line)
    outfile_variant.write("\n")

    varcov_outfile_name = outfile_prefix + "_varcov.txt"
    outfile_varcov = open(varcov_outfile_name,"w")
    varcov_file_line = "\t".join(["positions","overlap_clusters"])
    outfile_varcov.write(varcov_file_line)
    outfile_varcov.write("\n")

polya_outfile_name = outfile_prefix + "_polya.txt"
outfile_polya = open(polya_outfile_name,"w")
polya_file_line = "\t".join(["cluster_id","trans_id","strand","a_percent","a_count","sequence"])
outfile_polya.write(polya_file_line)
outfile_polya.write("\n")

#rtswitch_outfile_name = outfile_prefix + "_rtswitch.txt"
#outfile_rtswitch = open(rtswitch_outfile_name,"w")
#rtswitch_file_line = "\t".join(["trans_id","junct_num","first_seq","second_seq","rev_comp_first_seq"])
#outfile_rtswitch.write(rtswitch_file_line)
#outfile_rtswitch.write("\n")


strand_outfile_name = outfile_prefix + "_strand_check.txt"
outfile_strand = open(strand_outfile_name,"w")
strand_file_line = "\t".join(["read_id","scaff_name","start_pos","cigar","strands"])
outfile_strand.write(strand_file_line)
outfile_strand.write("\n")


lde_outfile_name = outfile_prefix + "_local_density_error.txt"
outfile_lde = open(lde_outfile_name,"w")
lde_file_line = "\t".join(["cluster_id","lde_flag","scaff_name","start_pos","end_pos","strand","num_exons","bad_sj_num_line","bad_sj_error_count_line","sj_error_profile_idmsh","sj_error_nuc","sj_error_simple","cigar"])
outfile_lde.write(lde_file_line)
outfile_lde.write("\n")

variation_dict = {} # variation_dict[scaffold][position][variant type][alt allele][cluster id] = 1
var_coverage_dict = {} # var_coverage_dict[scaffold][position][trans_id] = 1

## sj hash 2020/07/27

sj_hash_read_threshold = 20

check_trans_id = '11_c110717/1/696' #########################################################################debugging

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


#convert cigar into list of digits and list of characters
def cigar_list(cigar):
    # fix issue with X and = used in SAM format from mininmap2
    cigar = re.sub('=', 'M', cigar)
    cigar = re.sub('X', 'M', cigar)

    cig_char = re.sub('\d', ' ', cigar)
    cig_char_list = cig_char.split()

    cig_digit = re.sub("[a-zA-Z]+", ' ', cigar)
    cig_dig_list = cig_digit.split()
    
    return(cig_dig_list,cig_char_list)

####################################################################################################

#get mapped sequence length
def mapped_seq_length(cigar):
    [cig_dig_list,cig_char_list] = cigar_list(cigar)
    map_seq_length = 0
    
    for i in xrange(len(cig_dig_list)):
        cig_flag = cig_char_list[i]
        
        if cig_flag == "H":
            continue
        elif cig_flag == "S":
            continue
        elif cig_flag == "M":
            map_seq_length = map_seq_length + int(cig_dig_list[i])
            continue
        elif cig_flag == "I":
            continue
        elif cig_flag == "D":
            map_seq_length = map_seq_length + int(cig_dig_list[i])
            continue
        elif cig_flag == "N":
            continue
        
    return map_seq_length


####################################################################################################

#get the coordinate for the end of the transcript
def trans_coordinates(start_pos,cigar):
    [cig_dig_list,cig_char_list] = cigar_list(cigar)
    end_pos = int(start_pos)
    exon_start_list = []
    exon_end_list = []

    exon_start_list.append(int(start_pos))

    for i in xrange(len(cig_dig_list)):
        cig_flag = cig_char_list[i]
        
        if cig_flag == "H":
            continue
        elif cig_flag == "S":
            continue
        elif cig_flag == "M":
            end_pos = end_pos + int(cig_dig_list[i])
            continue
        elif cig_flag == "I":
            continue
        elif cig_flag == "D":
            end_pos = end_pos + int(cig_dig_list[i])
            continue
        elif cig_flag == "N":
            #add exon end position to list (must be before updating end pos info)
            exon_end_list.append(end_pos)
            end_pos = end_pos + int(cig_dig_list[i])
            #add exon start
            exon_start_list.append(end_pos)
            continue
    #add last exon end position
    exon_end_list.append(end_pos)
    
    return end_pos,exon_start_list,exon_end_list

####################################################################################################

#deal with wildcards

nuc_char_dict = {} # nuc_char_dict[nuc char][nuc single] = 1

nuc_char_dict["A"] = {}
nuc_char_dict["A"]["A"] = 1

nuc_char_dict["T"] = {}
nuc_char_dict["T"]["T"] = 1

nuc_char_dict["C"] = {}
nuc_char_dict["C"]["C"] = 1

nuc_char_dict["G"] = {}
nuc_char_dict["G"]["G"] = 1

nuc_char_dict["S"] = {}
nuc_char_dict["S"]["C"] = 1
nuc_char_dict["S"]["G"] = 1

nuc_char_dict["W"] = {}
nuc_char_dict["W"]["A"] = 1
nuc_char_dict["W"]["T"] = 1

nuc_char_dict["K"] = {}
nuc_char_dict["K"]["A"] = 1
nuc_char_dict["K"]["C"] = 1

nuc_char_dict["M"] = {}
nuc_char_dict["M"]["G"] = 1
nuc_char_dict["M"]["T"] = 1

nuc_char_dict["Y"] = {}
nuc_char_dict["Y"]["A"] = 1
nuc_char_dict["Y"]["G"] = 1

nuc_char_dict["R"] = {}
nuc_char_dict["R"]["C"] = 1
nuc_char_dict["R"]["T"] = 1

nuc_char_dict["V"] = {}
nuc_char_dict["V"]["C"] = 1
nuc_char_dict["V"]["G"] = 1
nuc_char_dict["V"]["T"] = 1

nuc_char_dict["H"] = {}
nuc_char_dict["H"]["A"] = 1
nuc_char_dict["H"]["G"] = 1
nuc_char_dict["H"]["T"] = 1

nuc_char_dict["D"] = {}
nuc_char_dict["D"]["A"] = 1
nuc_char_dict["D"]["C"] = 1
nuc_char_dict["D"]["T"] = 1

nuc_char_dict["B"] = {}
nuc_char_dict["B"]["A"] = 1
nuc_char_dict["B"]["C"] = 1
nuc_char_dict["B"]["G"] = 1

nuc_char_dict["N"] = {}
nuc_char_dict["N"]["A"] = 1
nuc_char_dict["N"]["T"] = 1
nuc_char_dict["N"]["C"] = 1
nuc_char_dict["N"]["G"] = 1
#nuc_char_dict["Z"] = {}



#use this to find mismatch between two aligned sequences
def mismatch_seq(genome_seq,query_seq,genome_pos,seq_pos):
    
    if len(genome_seq) != len(query_seq):
        print("Genome seq is not the same length as query seq")
        print(genome_seq)
        print(query_seq)
        print(genome_pos)
        print(seq_pos)
        
        sys.exit()
    
    genome_mismatch_list = []
    seq_mismatch_list = []
    nuc_mismatch_list = []
    
    for i in xrange(len(genome_seq)):

        genome_nuc = genome_seq[i]
        read_nuc = query_seq[i]

        read_nuc = read_nuc.upper()

        nuc_match_flag = 0
        for g_nuc_char in nuc_char_dict[genome_nuc]:
            for t_nuc_char in nuc_char_dict[read_nuc]:
                if g_nuc_char == t_nuc_char:
                    nuc_match_flag = 1


        #if genome_seq[i] != query_seq[i]: #########################################need to fix this for wildcard situations
        if nuc_match_flag == 0:
            genome_mismatch_coord = genome_pos + i
            seq_mismatch_coord = seq_pos + i
            genome_mismatch_list.append(genome_mismatch_coord)
            seq_mismatch_list.append(seq_mismatch_coord)
            
            nuc_mismatch_str = query_seq[i] + "." + genome_seq[i]
            nuc_mismatch_list.append(nuc_mismatch_str)
    
    return genome_mismatch_list,seq_mismatch_list,nuc_mismatch_list


####################################################################################################

def update_variation_dict(scaffold,var_pos,var_type,var_seq,read_id):
    
    if type(var_seq) is list:
        var_seq = "".join(var_seq)

    if scaffold not in variation_dict:
        variation_dict[scaffold] = {}        
    if var_pos not in variation_dict[scaffold]:
        variation_dict[scaffold][var_pos] = {}        
    if var_type not in variation_dict[scaffold][var_pos]:
        variation_dict[scaffold][var_pos][var_type] = {}
    if var_seq not in variation_dict[scaffold][var_pos][var_type]:
        variation_dict[scaffold][var_pos][var_type][var_seq] = {}
    if read_id not in variation_dict[scaffold][var_pos][var_type][var_seq]:
        variation_dict[scaffold][var_pos][var_type][var_seq][read_id] = 1
    
    #For variation coverage
    if scaffold not in var_coverage_dict:
        var_coverage_dict[scaffold] = {}###############for var cov
    if var_pos not in var_coverage_dict[scaffold]:
        var_coverage_dict[scaffold][var_pos] = {}###############for var cov
    var_coverage_dict[scaffold][var_pos][read_id] = 1

#used to calculate error rate of mapping
def calc_error_rate(start_pos,cigar,seq_list,scaffold,read_id):
    [cig_dig_list,cig_char_list] = cigar_list(cigar)
    
    h_count = 0
    s_count = 0
    i_count = 0
    d_count = 0
    mis_count = 0
    
    all_genome_mismatch_list = []
    all_nuc_mismatch_list = [] # this shows the nuc mismatch for the all_genome_mismatch_list

    insertion_list = []
    insertion_length_list = []
    deletion_list = []
    deletion_length_list = []

    sj_pre_error_list = []
    sj_post_error_list = []

    
    nomatch_dict = {} # nomatch_dict[non match id] = list
    nomatch_dict['mismatch_list'] = all_genome_mismatch_list
    nomatch_dict['insertion_list'] = insertion_list
    nomatch_dict['insertion_length_list'] = insertion_length_list
    nomatch_dict['deletion_list'] = deletion_list
    nomatch_dict['deletion_length_list'] = deletion_length_list
    
    #walk through both the query seq and the mapped seq
    genome_pos = start_pos - 1 #adjust for 1 base to 0 base coordinates
    seq_pos = 0
    
    for i in xrange(len(cig_dig_list)):
        cig_flag = cig_char_list[i]
        
        #print(seq_pos)
        
        if cig_flag == "H":
            h_count = h_count + int(cig_dig_list[i])
            
            ### for variation detection start
            var_pos = genome_pos - h_count 
            var_seq = "N" * h_count
            var_type = cig_flag
            update_variation_dict(scaffold,var_pos,var_type,var_seq,read_id)
            ### for variation detection end
            
            continue
        elif cig_flag == "S":
            s_count = s_count + int(cig_dig_list[i])
            
            ### for variation detection start
            seq_start = seq_pos
            seq_end = seq_pos + int(cig_dig_list[i])
            var_pos = genome_pos - s_count
            var_seq = seq_list[seq_start:seq_end]
            var_type = cig_flag
            
            ###### change made 2019/11/19
            #correct for soft clips at end of mapping
            num_cig_entities = len(cig_dig_list)
            if i == num_cig_entities - 1:
                var_pos = genome_pos + 1 - s_count
            ###### change made 2019/11/19
            
            update_variation_dict(scaffold,var_pos,var_type,var_seq,read_id)
            ### for variation detection end
            
            #adjust seq position to account for missing soft mask in genome coord
            seq_pos = seq_pos + int(cig_dig_list[i])

            continue
        elif cig_flag == "M":
            match_length = int(cig_dig_list[i])
            seq_start = seq_pos
            seq_end = seq_pos + match_length
            
            genome_start = genome_pos
            genome_end = genome_pos + match_length
            
            #seq_slice = seq_list[seq_start:seq_end+1]
            #genome_slice = fasta_dict[scaffold][genome_start:genome_end+1]
            
            seq_slice = seq_list[seq_start:seq_end]
            genome_slice = fasta_dict[scaffold][genome_start:genome_end]
            
            seq_diff = seq_end - seq_start
            gen_diff = genome_end - genome_start

            ######################################################
            #print(str(match_length)+cig_flag)
            ######################################################

            #get number and location of mis matches
            genome_mismatch_list,seq_mismatch_list,nuc_mismatch_list = mismatch_seq(genome_slice,seq_slice,genome_start,seq_start)
            
            mis_count = mis_count + len(genome_mismatch_list)
            all_genome_mismatch_list.extend(genome_mismatch_list)
            all_nuc_mismatch_list.extend(nuc_mismatch_list)
            
            ### for variation detection start
            for mismatch_index in xrange(len(seq_mismatch_list)):
                var_pos = genome_mismatch_list[mismatch_index]
                seq_var_pos = seq_mismatch_list[mismatch_index]
                var_seq = seq_list[seq_var_pos]
                var_type = cig_flag
                update_variation_dict(scaffold,var_pos,var_type,var_seq,read_id)
            ### for variation detection end
            
            seq_pos = seq_pos + match_length
            genome_pos = genome_pos + match_length
            
            continue
        elif cig_flag == "I":

            ### for variation detection start
            seq_start = seq_pos
            seq_end = seq_pos + int(cig_dig_list[i])
            var_pos = genome_pos
            var_seq = seq_list[seq_start:seq_end]
            var_type = cig_flag
            update_variation_dict(scaffold,var_pos,var_type,var_seq,read_id)
            ### for variation detection end
            
            insertion_length = int(cig_dig_list[i])
            seq_pos = seq_pos + insertion_length
            i_count = i_count + insertion_length
            
            insertion_list.append(genome_pos)
            insertion_length_list.append(insertion_length)
            
            continue
        elif cig_flag == "D":
            
            deletion_list.append(genome_pos)
            deletion_length = int(cig_dig_list[i])
            deletion_length_list.append(deletion_length)
            
            ### for variation detection start
            var_pos = genome_pos
            var_seq = str(deletion_length) #var seq for deletion is length of deletion
            var_type = cig_flag
            update_variation_dict(scaffold,var_pos,var_type,var_seq,read_id)
            ### for variation detection end
            
            genome_pos = genome_pos + deletion_length
            d_count = d_count + deletion_length
            
            continue
        elif cig_flag == "N":

            sj_pre_error_builder_list = [] # list of cigar string before splice junctions
            sj_post_error_builder_list = [] # list of cigar string after splice junctions


            # check for errors before splice junction pre
            #########################################
            prev_cig_flag = cig_char_list[i - 1]
            prev_cig_length = int(cig_dig_list[i - 1])

            #############################################################################################################################
            # Add sj error info
            
            prev_total_mismatch_length = 0
            this_mismatch_length = 0 # keep track of where in  the error list of cig

            all_genome_mismatch_index = len(all_genome_mismatch_list) - 1

            prev_total_cig_length = prev_cig_length
            prev_cig_index = i - 1
            
            this_cig_genome_pos = genome_pos

            prev_sj_flag = 0

            while prev_total_mismatch_length <= sj_err_threshold and prev_sj_flag == 0:
                this_mismatch_length = this_mismatch_length + prev_cig_length
                prev_total_mismatch_length = this_mismatch_length

                if prev_cig_flag == "M":
                    # if no mismatch
                    if len(all_genome_mismatch_list) < 1: # no mismatch from M cigar regions up to this point 
                        #last_genome_mismatch = -2 * sj_err_threshold
                        
                        if prev_total_mismatch_length <= sj_err_threshold: # no mis matches for this gene
                            sj_pre_error_builder_list.append(str(prev_cig_length) + prev_cig_flag)

                    elif all_genome_mismatch_index >= 0:
                        last_genome_mismatch = all_genome_mismatch_list[all_genome_mismatch_index]
                        last_nuc_mismatch = all_nuc_mismatch_list[all_genome_mismatch_index]

                        dist_prev_mismatch = genome_pos - last_genome_mismatch

                        m_builder_add_count = 0 

                        #while dist_prev_mismatch <= this_mismatch_length and dist_prev_mismatch <= sj_err_threshold:
                        while last_genome_mismatch >= this_cig_genome_pos - prev_cig_length and last_genome_mismatch <= this_cig_genome_pos and dist_prev_mismatch <= sj_err_threshold: # dist_prev_mismatch is still in range of sj threshold

                            sj_pre_error_builder_list.append(str(dist_prev_mismatch) + "." + last_nuc_mismatch)

                            m_builder_add_count += 1

                            all_genome_mismatch_index -= 1

                            last_genome_mismatch = all_genome_mismatch_list[all_genome_mismatch_index]
                            last_nuc_mismatch = all_nuc_mismatch_list[all_genome_mismatch_index]
                            dist_prev_mismatch = genome_pos - last_genome_mismatch



                            if all_genome_mismatch_index < 0:
                                break
                        
                        if m_builder_add_count == 0:
                            if prev_total_mismatch_length <= sj_err_threshold: #  still within in sj threshold so output
                                sj_pre_error_builder_list.append(str(prev_cig_length) + prev_cig_flag)
                        
                    elif all_genome_mismatch_index < 0: #  length of mismatch list is not 0 but index has gone past that
                        #last_genome_mismatch = -2 * sj_err_threshold
                        if prev_total_mismatch_length <= sj_err_threshold: #  still within in sj threshold so output
                            sj_pre_error_builder_list.append(str(prev_cig_length) + prev_cig_flag)
                    else:
                        print("Error with all_genome_mismatch_index")
                        print(all_genome_mismatch_list)
                        print(all_genome_mismatch_index)
                        sys.exit()
                elif prev_cig_flag != "N":

                    prev_cig_flag = cig_char_list[prev_cig_index]
                    prev_cig_length = int(cig_dig_list[prev_cig_index])

                    sj_pre_error_builder_list.append(str(prev_cig_length) + prev_cig_flag)

                elif prev_cig_flag == "N":
                    prev_sj_flag = 1

                prev_cig_index -= 1

                if prev_cig_index < 0 :
                    break
                
                this_cig_genome_pos = this_cig_genome_pos - prev_cig_length
                
                prev_cig_flag = cig_char_list[prev_cig_index]
                prev_cig_length = int(cig_dig_list[prev_cig_index])

            if len(sj_pre_error_builder_list) == 0: # in case no errors are found
                sj_pre_error_builder_list.append(no_mismatch_flag)

            # end of adding sj error info prev
            #############################################################################################################################

            #########################################
            # not related to sj error should probably move this up or down code of this block
            #########################################
            intron_length = int(cig_dig_list[i])
            genome_pos = genome_pos + intron_length
            #########################################

            #check for errors after splice junction post
            #########################################
            next_cig_flag = cig_char_list[i + 1]
            next_cig_length = int(cig_dig_list[i + 1])

            #############################################################################################################################

            next_total_mismatch_length = 0
            this_mismatch_length = 0 # keep track of where in  the error list of cig

            all_genome_mismatch_index = len(all_genome_mismatch_list) - 1

            next_total_cig_length = next_cig_length
            next_cig_index = i + 1

            this_next_seq_pos = seq_pos
            this_genome_pos = genome_pos

            #genome_mismatch_index = 0
            next_sj_flag = 0

            while next_total_mismatch_length <= sj_err_threshold and next_sj_flag == 0:
                this_mismatch_length = this_mismatch_length + next_cig_length
                next_total_mismatch_length = this_mismatch_length

                if next_cig_flag == "M":
                    match_length = next_cig_length
                    seq_start = this_next_seq_pos
                    seq_end = this_next_seq_pos + match_length

                    genome_start = this_genome_pos
                    genome_end = this_genome_pos + match_length

                    seq_slice = seq_list[seq_start:seq_end]
                    genome_slice = fasta_dict[scaffold][genome_start:genome_end]

                    #######################################################
                    #print(str(match_length)+next_cig_flag)
                    #print(seq_start)
                    #print(seq_end)
                    ######################################################

                    genome_mismatch_list, seq_mismatch_list, nuc_mismatch_list = mismatch_seq(genome_slice, seq_slice, genome_start, seq_start)

                    genome_mismatch_index = 0

                    # if no mismatch
                    if len(genome_mismatch_list) < 1: # no mismatch in this cig entry
                        #dist_next_mismatch = 2 * sj_err_threshold

                        if next_total_mismatch_length <= sj_err_threshold:

                            sj_post_error_builder_list.append(str(next_cig_length) + next_cig_flag)

                    elif genome_mismatch_index <= len(genome_mismatch_list):
                        nuc_next_mismatch = nuc_mismatch_list[0]
                        #dist_next_mismatch = genome_mismatch_list[0] - this_genome_pos
                        dist_next_mismatch = genome_mismatch_list[0] - genome_pos

                        m_builder_add_count = 0 

                        while dist_next_mismatch <= this_mismatch_length and dist_next_mismatch < sj_err_threshold:
                            next_genome_mismatch = genome_mismatch_list[genome_mismatch_index]
                            next_nuc_mismatch = nuc_mismatch_list[genome_mismatch_index]

                            #dist_next_mismatch = next_genome_mismatch - this_genome_pos
                            dist_next_mismatch = next_genome_mismatch - genome_pos

                            if dist_next_mismatch <= sj_err_threshold: # this mismatch goes beyond this M length

                                sj_post_error_builder_list.append(str(dist_next_mismatch) + "." + next_nuc_mismatch)

                                m_builder_add_count += 1


                            genome_mismatch_index += 1
                            if genome_mismatch_index >= len(genome_mismatch_list):
                                break
                        
                        if m_builder_add_count == 0:
                            if next_total_mismatch_length <= sj_err_threshold:
                                sj_post_error_builder_list.append(str(next_cig_length) + next_cig_flag)
                        

                    else:
                        print("Error with genome_mismatch_index")
                        sys.exit()
                elif next_cig_flag != "N":

                    next_cig_flag = cig_char_list[next_cig_index]
                    next_cig_length = int(cig_dig_list[next_cig_index])

                    sj_post_error_builder_list.append(str(next_cig_length) + next_cig_flag)
                elif next_cig_flag == "N":
                    next_sj_flag = 1



                if next_cig_flag != "D": # increase this seq pos only if not a deletion
                    this_next_seq_pos = this_next_seq_pos + next_cig_length

                if next_cig_flag != "I":
                        this_genome_pos = this_genome_pos + next_cig_length

                next_cig_index += 1

                if next_cig_index >= len(cig_char_list):
                    break

                next_cig_flag = cig_char_list[next_cig_index]
                next_cig_length = int(cig_dig_list[next_cig_index])



            if len(sj_post_error_builder_list) == 0: # in case no errors are found
                sj_post_error_builder_list.append(no_mismatch_flag)

            # end of adding sj error info post
            #############################################################################################################################
            #############################################################################################################################


            sj_post_error_builder_line = "_".join(sj_post_error_builder_list)

            sj_pre_error_builder_list_reverse = list(reversed(sj_pre_error_builder_list))
            sj_pre_error_builder_line = "_".join(sj_pre_error_builder_list_reverse)

            sj_pre_error_list.append(sj_pre_error_builder_line)
            sj_post_error_list.append(sj_post_error_builder_line)


            #########################################



            continue

        else:
            match_length = int(cig_dig_list[i])

            print("Error with cigar flag")
            print(str(match_length) + cig_flag)
            print(cig_dig_list)
            print(cig_char_list)
            sys.exit()
    
    
    #print("calc_error_rate")
    #print(cigar)
    #blah = ",".join([str(h_count),str(s_count),str(i_count),str(d_count),str(mis_count)])
    #print(blah)
    #print(cig_dig_list)
    #print(cig_char_list)
        
    #sys.exit()
    
    return h_count,s_count,i_count,d_count,mis_count,nomatch_dict,sj_pre_error_list,sj_post_error_list


####################################################################################################

#get variation information without calc_error_rate

#used to calculate error rate of mapping
def calc_variation(start_pos,cigar,seq_list,scaffold,read_id):
    [cig_dig_list,cig_char_list] = cigar_list(cigar)

    h_count = 0
    s_count = 0
    i_count = 0
    d_count = 0
    mis_count = 0



    #walk through both the query seq and the mapped seq
    genome_pos = start_pos - 1 #adjust for 1 base to 0 base coordinates
    seq_pos = 0

    for i in xrange(len(cig_dig_list)):
        cig_flag = cig_char_list[i]

        #print(seq_pos)

        if cig_flag == "H":
            h_count = h_count + int(cig_dig_list[i])

            ### for variation detection start
            var_pos = genome_pos - h_count
            var_seq = "N" * h_count
            var_type = cig_flag
            update_variation_dict(scaffold,var_pos,var_type,var_seq,read_id)
            ### for variation detection end

            continue
        elif cig_flag == "S":
            s_count = s_count + int(cig_dig_list[i])

            ### for variation detection start
            seq_start = seq_pos
            seq_end = seq_pos + int(cig_dig_list[i])
            var_pos = genome_pos - s_count
            var_seq = seq_list[seq_start:seq_end]
            var_type = cig_flag

            ###### change made 2019/11/19
            #correct for soft clips at end of mapping
            num_cig_entities = len(cig_dig_list)
            if i == num_cig_entities - 1:
                var_pos = genome_pos + 1 - s_count
            ###### change made 2019/11/19

            update_variation_dict(scaffold,var_pos,var_type,var_seq,read_id)
            ### for variation detection end

            #adjust seq position to account for missing soft mask in genome coord
            seq_pos = seq_pos + int(cig_dig_list[i])

            continue
        elif cig_flag == "M":
            match_length = int(cig_dig_list[i])
            seq_start = seq_pos
            seq_end = seq_pos + match_length

            genome_start = genome_pos
            genome_end = genome_pos + match_length

            #seq_slice = seq_list[seq_start:seq_end+1]
            #genome_slice = fasta_dict[scaffold][genome_start:genome_end+1]

            seq_slice = seq_list[seq_start:seq_end]
            genome_slice = fasta_dict[scaffold][genome_start:genome_end]

            seq_diff = seq_end - seq_start
            gen_diff = genome_end - genome_start

            ######################################################
            #print(str(match_length)+cig_flag)
            ######################################################

            #get number and location of mis matches
            genome_mismatch_list,seq_mismatch_list,nuc_mismatch_list = mismatch_seq(genome_slice,seq_slice,genome_start,seq_start)


            ### for variation detection start
            for mismatch_index in xrange(len(seq_mismatch_list)):
                var_pos = genome_mismatch_list[mismatch_index]
                seq_var_pos = seq_mismatch_list[mismatch_index]
                var_seq = seq_list[seq_var_pos]
                var_type = cig_flag
                update_variation_dict(scaffold,var_pos,var_type,var_seq,read_id)
            ### for variation detection end

            seq_pos = seq_pos + match_length
            genome_pos = genome_pos + match_length

            continue
        elif cig_flag == "I":

            ### for variation detection start
            seq_start = seq_pos
            seq_end = seq_pos + int(cig_dig_list[i])
            var_pos = genome_pos
            var_seq = seq_list[seq_start:seq_end]
            var_type = cig_flag
            update_variation_dict(scaffold,var_pos,var_type,var_seq,read_id)
            ### for variation detection end

            insertion_length = int(cig_dig_list[i])
            seq_pos = seq_pos + insertion_length



            continue
        elif cig_flag == "D":

            deletion_length = int(cig_dig_list[i])

            ### for variation detection start
            var_pos = genome_pos
            var_seq = str(deletion_length) #var seq for deletion is length of deletion
            var_type = cig_flag
            update_variation_dict(scaffold,var_pos,var_type,var_seq,read_id)
            ### for variation detection end

            genome_pos = genome_pos + deletion_length
            d_count = d_count + deletion_length

            continue
        elif cig_flag == "N":


            #####################################################################################################################

            #########################################
            # not related to sj error should probably move this up or down code of this block
            #########################################
            intron_length = int(cig_dig_list[i])
            genome_pos = genome_pos + intron_length
            #########################################

            #############################################################################################################################

            continue

        else:
            match_length = int(cig_dig_list[i])

            print("Error with cigar flag")
            print(str(match_length) + cig_flag)
            print(cig_dig_list)
            print(cig_char_list)
            sys.exit()




#####################################################################################################

## New calc error for low mem mode
#used to calculate error rate of mapping
def calc_error_rate_lowmem(start_pos,cigar,seq_list,scaffold,read_id):
    [cig_dig_list,cig_char_list] = cigar_list(cigar)

    h_count = 0
    s_count = 0
    i_count = 0
    d_count = 0
    mis_count = 0

    all_genome_mismatch_list = []
    all_nuc_mismatch_list = [] # this shows the nuc mismatch for the all_genome_mismatch_list

    insertion_list = []
    insertion_length_list = []
    deletion_list = []
    deletion_length_list = []

    sj_pre_error_list = []
    sj_post_error_list = []


    nomatch_dict = {} # nomatch_dict[non match id] = list
    nomatch_dict['mismatch_list'] = all_genome_mismatch_list
    nomatch_dict['insertion_list'] = insertion_list
    nomatch_dict['insertion_length_list'] = insertion_length_list
    nomatch_dict['deletion_list'] = deletion_list
    nomatch_dict['deletion_length_list'] = deletion_length_list

    #walk through both the query seq and the mapped seq
    genome_pos = start_pos - 1 #adjust for 1 base to 0 base coordinates
    seq_pos = 0

    for i in xrange(len(cig_dig_list)):
        cig_flag = cig_char_list[i]

        #print(seq_pos)

        if cig_flag == "H":
            h_count = h_count + int(cig_dig_list[i])

            continue
        elif cig_flag == "S":
            s_count = s_count + int(cig_dig_list[i])

            #adjust seq position to account for missing soft mask in genome coord
            seq_pos = seq_pos + int(cig_dig_list[i])

            continue
        elif cig_flag == "M":
            match_length = int(cig_dig_list[i])
            seq_start = seq_pos
            seq_end = seq_pos + match_length

            genome_start = genome_pos
            genome_end = genome_pos + match_length

            seq_slice = seq_list[seq_start:seq_end]
            genome_slice = fasta_dict[scaffold][genome_start:genome_end]

            seq_diff = seq_end - seq_start
            gen_diff = genome_end - genome_start


            #get number and location of mis matches
            genome_mismatch_list,seq_mismatch_list,nuc_mismatch_list = mismatch_seq(genome_slice,seq_slice,genome_start,seq_start)

            mis_count = mis_count + len(genome_mismatch_list)
            all_genome_mismatch_list.extend(genome_mismatch_list)
            all_nuc_mismatch_list.extend(nuc_mismatch_list)

            seq_pos = seq_pos + match_length
            genome_pos = genome_pos + match_length

            continue
        elif cig_flag == "I":

            insertion_length = int(cig_dig_list[i])
            seq_pos = seq_pos + insertion_length
            i_count = i_count + insertion_length

            insertion_list.append(genome_pos)
            insertion_length_list.append(insertion_length)

            continue
        elif cig_flag == "D":

            deletion_list.append(genome_pos)
            deletion_length = int(cig_dig_list[i])
            deletion_length_list.append(deletion_length)

            genome_pos = genome_pos + deletion_length
            d_count = d_count + deletion_length

            continue
        elif cig_flag == "N":

            sj_pre_error_builder_list = [] # list of cigar string before splice junctions
            sj_post_error_builder_list = [] # list of cigar string after splice junctions


            # check for errors before splice junction pre
            #########################################
            prev_cig_flag = cig_char_list[i - 1]
            prev_cig_length = int(cig_dig_list[i - 1])

            #############################################################################################################################
            # Add sj error info

            prev_total_mismatch_length = 0
            this_mismatch_length = 0 # keep track of where in  the error list of cig

            all_genome_mismatch_index = len(all_genome_mismatch_list) - 1

            prev_total_cig_length = prev_cig_length
            prev_cig_index = i - 1

            this_cig_genome_pos = genome_pos

            prev_sj_flag = 0

            while prev_total_mismatch_length <= sj_err_threshold and prev_sj_flag == 0:
                this_mismatch_length = this_mismatch_length + prev_cig_length
                prev_total_mismatch_length = this_mismatch_length

                if prev_cig_flag == "M":
                    # if no mismatch
                    if len(all_genome_mismatch_list) < 1: # no mismatch from M cigar regions up to this point
                        #last_genome_mismatch = -2 * sj_err_threshold

                        if prev_total_mismatch_length <= sj_err_threshold: # no mis matches for this gene
                            sj_pre_error_builder_list.append(str(prev_cig_length) + prev_cig_flag)

                    elif all_genome_mismatch_index >= 0:
                        last_genome_mismatch = all_genome_mismatch_list[all_genome_mismatch_index]
                        last_nuc_mismatch = all_nuc_mismatch_list[all_genome_mismatch_index]

                        dist_prev_mismatch = genome_pos - last_genome_mismatch

                        m_builder_add_count = 0

                        #while dist_prev_mismatch <= this_mismatch_length and dist_prev_mismatch <= sj_err_threshold:
                        while last_genome_mismatch >= this_cig_genome_pos - prev_cig_length and last_genome_mismatch <= this_cig_genome_pos and dist_prev_mismatch <= sj_err_threshold: # dist_prev_mismatch is still in range of sj threshold

                            sj_pre_error_builder_list.append(str(dist_prev_mismatch) + "." + last_nuc_mismatch)

                            m_builder_add_count += 1

                            all_genome_mismatch_index -= 1

                            last_genome_mismatch = all_genome_mismatch_list[all_genome_mismatch_index]
                            last_nuc_mismatch = all_nuc_mismatch_list[all_genome_mismatch_index]
                            dist_prev_mismatch = genome_pos - last_genome_mismatch



                            if all_genome_mismatch_index < 0:
                                break

                        if m_builder_add_count == 0:
                            if prev_total_mismatch_length <= sj_err_threshold: #  still within in sj threshold so output
                                sj_pre_error_builder_list.append(str(prev_cig_length) + prev_cig_flag)

                    elif all_genome_mismatch_index < 0: #  length of mismatch list is not 0 but index has gone past that
                        #last_genome_mismatch = -2 * sj_err_threshold
                        if prev_total_mismatch_length <= sj_err_threshold: #  still within in sj threshold so output
                            sj_pre_error_builder_list.append(str(prev_cig_length) + prev_cig_flag)
                    else:
                        print("Error with all_genome_mismatch_index")
                        print(all_genome_mismatch_list)
                        print(all_genome_mismatch_index)
                        sys.exit()
                elif prev_cig_flag != "N":

                    prev_cig_flag = cig_char_list[prev_cig_index]
                    prev_cig_length = int(cig_dig_list[prev_cig_index])

                    sj_pre_error_builder_list.append(str(prev_cig_length) + prev_cig_flag)

                elif prev_cig_flag == "N":
                    prev_sj_flag = 1

                prev_cig_index -= 1

                if prev_cig_index < 0 :
                    break

                this_cig_genome_pos = this_cig_genome_pos - prev_cig_length

                prev_cig_flag = cig_char_list[prev_cig_index]
                prev_cig_length = int(cig_dig_list[prev_cig_index])

            if len(sj_pre_error_builder_list) == 0: # in case no errors are found
                sj_pre_error_builder_list.append(no_mismatch_flag)

            # end of adding sj error info prev
            #############################################################################################################################

            #########################################
            # not related to sj error should probably move this up or down code of this block
            #########################################
            intron_length = int(cig_dig_list[i])
            genome_pos = genome_pos + intron_length
            #########################################

            #check for errors after splice junction post
            #########################################
            next_cig_flag = cig_char_list[i + 1]
            next_cig_length = int(cig_dig_list[i + 1])

            #############################################################################################################################

            next_total_mismatch_length = 0
            this_mismatch_length = 0 # keep track of where in  the error list of cig

            all_genome_mismatch_index = len(all_genome_mismatch_list) - 1

            next_total_cig_length = next_cig_length
            next_cig_index = i + 1

            this_next_seq_pos = seq_pos
            this_genome_pos = genome_pos

            #genome_mismatch_index = 0
            next_sj_flag = 0

            while next_total_mismatch_length <= sj_err_threshold and next_sj_flag == 0:
                this_mismatch_length = this_mismatch_length + next_cig_length
                next_total_mismatch_length = this_mismatch_length

                if next_cig_flag == "M":
                    match_length = next_cig_length
                    seq_start = this_next_seq_pos
                    seq_end = this_next_seq_pos + match_length

                    genome_start = this_genome_pos
                    genome_end = this_genome_pos + match_length

                    seq_slice = seq_list[seq_start:seq_end]
                    genome_slice = fasta_dict[scaffold][genome_start:genome_end]

                    #######################################################
                    #print(str(match_length)+next_cig_flag)
                    #print(seq_start)
                    #print(seq_end)
                    ######################################################

                    genome_mismatch_list, seq_mismatch_list, nuc_mismatch_list = mismatch_seq(genome_slice, seq_slice, genome_start, seq_start)

                    genome_mismatch_index = 0

                    # if no mismatch
                    if len(genome_mismatch_list) < 1: # no mismatch in this cig entry
                        #dist_next_mismatch = 2 * sj_err_threshold

                        if next_total_mismatch_length <= sj_err_threshold:

                            sj_post_error_builder_list.append(str(next_cig_length) + next_cig_flag)

                    elif genome_mismatch_index <= len(genome_mismatch_list):
                        nuc_next_mismatch = nuc_mismatch_list[0]
                        #dist_next_mismatch = genome_mismatch_list[0] - this_genome_pos
                        dist_next_mismatch = genome_mismatch_list[0] - genome_pos

                        m_builder_add_count = 0

                        while dist_next_mismatch <= this_mismatch_length and dist_next_mismatch < sj_err_threshold:
                            next_genome_mismatch = genome_mismatch_list[genome_mismatch_index]
                            next_nuc_mismatch = nuc_mismatch_list[genome_mismatch_index]

                            #dist_next_mismatch = next_genome_mismatch - this_genome_pos
                            dist_next_mismatch = next_genome_mismatch - genome_pos

                            if dist_next_mismatch <= sj_err_threshold: # this mismatch goes beyond this M length

                                sj_post_error_builder_list.append(str(dist_next_mismatch) + "." + next_nuc_mismatch)

                                m_builder_add_count += 1


                            genome_mismatch_index += 1
                            if genome_mismatch_index >= len(genome_mismatch_list):
                                break

                        if m_builder_add_count == 0:
                            if next_total_mismatch_length <= sj_err_threshold:
                                sj_post_error_builder_list.append(str(next_cig_length) + next_cig_flag)


                    else:
                        print("Error with genome_mismatch_index")
                        sys.exit()
                elif next_cig_flag != "N":

                    next_cig_flag = cig_char_list[next_cig_index]
                    next_cig_length = int(cig_dig_list[next_cig_index])

                    sj_post_error_builder_list.append(str(next_cig_length) + next_cig_flag)
                elif next_cig_flag == "N":
                    next_sj_flag = 1



                if next_cig_flag != "D": # increase this seq pos only if not a deletion
                    this_next_seq_pos = this_next_seq_pos + next_cig_length

                if next_cig_flag != "I":
                        this_genome_pos = this_genome_pos + next_cig_length

                next_cig_index += 1

                if next_cig_index >= len(cig_char_list):
                    break

                next_cig_flag = cig_char_list[next_cig_index]
                next_cig_length = int(cig_dig_list[next_cig_index])



            if len(sj_post_error_builder_list) == 0: # in case no errors are found
                sj_post_error_builder_list.append(no_mismatch_flag)

            # end of adding sj error info post
            #############################################################################################################################
            #############################################################################################################################


            sj_post_error_builder_line = "_".join(sj_post_error_builder_list)

            sj_pre_error_builder_list_reverse = list(reversed(sj_pre_error_builder_list))
            sj_pre_error_builder_line = "_".join(sj_pre_error_builder_list_reverse)

            sj_pre_error_list.append(sj_pre_error_builder_line)
            sj_post_error_list.append(sj_post_error_builder_line)


            #########################################



            continue

        else:
            match_length = int(cig_dig_list[i])

            print("Error with cigar flag")
            print(str(match_length) + cig_flag)
            print(cig_dig_list)
            print(cig_char_list)
            sys.exit()

    return h_count,s_count,i_count,d_count,mis_count,nomatch_dict,sj_pre_error_list,sj_post_error_list

####################################################################################################

class Transcript:
    def __init__(self, cluster_id):
        self.cluster_id = cluster_id
        self.trans_id = cluster_id #for flexible calling
        self.sam_flag = "none"
        self.scaff_name = "none"
        self.start_pos = -1
        self.cigar = "none"
        self.read_seq = "none"
        self.seq_length = 0
        self.seq_list = []
        
        self.map_seq_length = 0
        self.h_count = "none"
        self.s_count = "none"
        self.i_count = "none"
        self.d_count = "none"
        self.mis_count = "none"
        self.nomatch_dict = {}
        
        self.end_pos = 0
        self.exon_start_list = [] 
        self.exon_end_list = []
        
        self.num_exons = 0
        
        self.strand = ""
        
        self.downstream_seq = []
        self.dseq_length = 0
        self.a_count = 0
        self.n_count = 0
        self.a_percent = 0.0
        self.n_percent = 0.0
        
        self.percent_cov = 0.0
        self.percent_identity = 0.0

        self.sj_hash = "none"  # use this to speed up matching
        
    
    def add_sam_info(self,sam_flag,scaff_name,start_pos,cigar,read_seq,seq_list):
        self.sam_flag = sam_flag
        self.scaff_name = scaff_name
        self.start_pos = int(start_pos)
        self.cigar = cigar
        self.read_seq = read_seq
        self.seq_list = seq_list
        
        #need to add hard clipped seq later if there is any
        self.seq_length = len(read_seq)
        
        if sam_flag_dict[sam_flag] == "forward_strand":
            self.strand = "+"
        elif sam_flag_dict[sam_flag] == "reverse_strand":
            self.strand = "-"
        else:
            print("Error with interpretting SAM flag")
            sys.exit()
    
    def add_map_seq_length(self,map_seq_length):
        self.map_seq_length = map_seq_length
    
    def add_exon_coords(self,end_pos,exon_start_list,exon_end_list):
        self.exon_start_list = exon_start_list
        self.exon_end_list = exon_end_list
        self.end_pos = end_pos
        
        self.num_exons = len(exon_start_list)


    def make_sj_hash_string(self):
        ################### 2020/07/26
        # make sj hash
        exon_start_string_list = []
        exon_end_string_list = []

        for exon_index in xrange(len(self.exon_start_list)):
            exon_start_string_list.append(str(self.exon_start_list[exon_index]))
            exon_end_string_list.append(str(self.exon_end_list[exon_index]))

        sj_right_list = exon_start_string_list
        sj_right_list.pop(0)

        sj_left_list = exon_end_string_list
        sj_left_list.pop(-1)

        sj_left_line = ",".join(sj_left_list)
        sj_right_line = ",".join(sj_right_list)
        self.sj_hash = ";".join([sj_left_line,sj_right_line])
        ###################

    def make_sj_hash_int(self):
        ################### 2020/07/26
        # make sj hash
        self.sj_hash = 0


        for exon_index in xrange(len(self.exon_start_list)-1):
            self.sj_hash = self.sj_hash + (self.exon_start_list[exon_index+1] * self.exon_end_list[exon_index])


        ###################

    
    def add_mismatch(self,h_count,s_count,i_count,d_count,mis_count,nomatch_dict,sj_pre_error_list,sj_post_error_list):
        self.h_count = int(h_count)
        self.s_count = int(s_count)
        self.i_count = int(i_count)
        self.d_count = int(d_count)
        self.mis_count = int(mis_count)
        self.nomatch_dict = nomatch_dict

        self.sj_pre_error_list = sj_pre_error_list

        self.sj_post_error_list = sj_post_error_list
        
        #hard clipped seq is missing in given read seq so need to add
        self.seq_length = self.seq_length + h_count
        
    
    def calc_coverage(self):
        
        percent_cov = float(self.seq_length - self.h_count - self.s_count) / float(self.seq_length)
        percent_cov = percent_cov * 100.0
        
        self.percent_cov = percent_cov
        
        return percent_cov
    
    def calc_identity(self):
        
        map_seq_length = self.seq_length - self.h_count - self.s_count
        
        if ident_calc_method == "ident_cov":
        
            nonmatch_count = self.h_count + self.s_count + self.i_count + self.d_count + self.mis_count
            percent_identity = float(self.seq_length - nonmatch_count) / float(self.seq_length)
            percent_identity = percent_identity * 100.0
        
        elif ident_calc_method == "ident_map":
            nonmatch_count = self.i_count + self.d_count + self.mis_count
            percent_identity = float(map_seq_length - nonmatch_count) / float(map_seq_length)
            percent_identity = percent_identity * 100.0
            
        self.percent_identity = percent_identity
        
        return percent_identity
    
    def make_error_line(self):
        
        error_list = []
        error_list.append(str(self.h_count))
        error_list.append(str(self.s_count))
        error_list.append(str(self.i_count))
        error_list.append(str(self.d_count))
        error_list.append(str(self.mis_count))
        
        error_line = ";".join(error_list)
        
        return error_line
    
    def make_exon_start_end_lines(self):
        
        exon_start_string_list = []
        exon_end_string_list = []
        
        for i in xrange(len(self.exon_start_list)):
            
            exon_start_string_list.append(str(self.exon_start_list[i]))
            exon_end_string_list.append(str(self.exon_end_list[i]))
        
        exon_start_string_line = ",".join(exon_start_string_list)
        exon_end_string_line = ",".join(exon_end_string_list)
        
        return exon_start_string_line,exon_end_string_line
    
    def add_polya_info(self,downstream_seq,dseq_length,a_count,n_count,a_percent,n_percent):
        self.downstream_seq = downstream_seq
        self.dseq_length = dseq_length
        self.a_count = a_count
        self.n_count = n_count
        self.a_percent = a_percent
        self.n_percent = n_percent
    
    
    def format_bed_line(self,final_trans_id):
        bed_list = []
        bed_list.append(str(self.scaff_name))
        bed_list.append(str(self.start_pos-1))
        bed_list.append(str(self.end_pos-1))
        
        gene_id = self.trans_id.split(".")[0]
        id_line = ";".join([final_trans_id,self.trans_id])
        
        bed_list.append(str(id_line))
        bed_list.append("40")
        bed_list.append(self.strand)
        
        bed_list.append(str(self.start_pos-1))
        bed_list.append(str(self.end_pos-1))
        
        bed_list.append("255,0,0")
        
        bed_list.append(str(self.num_exons))
        
        relative_exon_start_list = []
        exon_length_list = []
        for i in xrange(self.num_exons):
            exon_start = self.exon_start_list[i]
            exon_end = self.exon_end_list[i]
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


####################################################################################################

class Merged:
    def __init__(self, trans_id):
        self.trans_id = trans_id
        
        #same as collapse start and end list but used for flexibility in calling
        self.exon_start_list = [] 
        self.exon_end_list = []
        
        self.scaff_name = "none"
        
        self.merged_trans_dict = {} # merged_trans_dict[trans id] = trans obj
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

    
    def add_merged_trans(self,trans_obj):
        merged_trans_id = trans_obj.cluster_id
        
        #self.trans_list.append(merged_trans_id)
        self.trans_obj_list.append(trans_obj)
        
        self.merged_trans_dict[merged_trans_id] = trans_obj
        
        merged_trans_id_list = self.merged_trans_dict.keys()
        self.trans_list = self.merged_trans_dict.keys()

        self.num_trans = len(merged_trans_id_list)
        
        if self.num_exons < trans_obj.num_exons:
            self.num_exons = trans_obj.num_exons
        
        self.scaff_name = trans_obj.scaff_name
        
        if self.strand == "none":
            self.strand = trans_obj.strand
        elif self.strand != trans_obj.strand:
            print("Error with merged trans not on the same strand")
            sys.exit()

    def add_merge_info(self,collapse_start_list,collapse_end_list,start_wobble_list,end_wobble_list,collapse_sj_start_err_list,collapse_sj_end_err_list,collapse_start_error_nuc_list,collapse_end_error_nuc_list ):
        self.collapse_start_list = collapse_start_list
        self.collapse_end_list = collapse_end_list
        self.start_wobble_list = start_wobble_list
        self.end_wobble_list = end_wobble_list
        
        self.exon_start_list = collapse_start_list
        self.exon_end_list = collapse_end_list

        self.start_pos = collapse_start_list[0]
        self.end_pos = collapse_end_list[-1]

        self.collapse_sj_start_err_list = collapse_sj_start_err_list
        self.collapse_sj_end_err_list = collapse_sj_end_err_list

        self.collapse_start_error_nuc_list = collapse_start_error_nuc_list
        self.collapse_end_error_nuc_list = collapse_end_error_nuc_list

    def format_bed_line(self):
        bed_list = []
        bed_list.append(str(self.scaff_name))
        bed_list.append(str(self.start_pos-1))
        bed_list.append(str(self.end_pos-1))
        
        gene_id = self.trans_id.split(".")[0]
        id_line = ";".join([gene_id,self.trans_id])
        
        bed_list.append(str(id_line))
        bed_list.append("40")
        bed_list.append(self.strand)
        
        bed_list.append(str(self.start_pos-1))
        bed_list.append(str(self.end_pos-1))
        
        bed_list.append("255,0,0")
        
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
            
            if(relative_exon_start < 0):
                print("negative exon_length")
                print(exon_length)
                print(id_line)
                sys.exit()
            
            if(relative_exon_start < 0):
                print("negative relative_exon_start")
                print(relative_exon_start)
                print(id_line)
                sys.exit()
            
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
        
        high_quality_percent = -1.0
        low_quality_percent = -1.0
        high_cov_percent = -1.0
        low_cov_percent = -1.0
        
        for merged_trans_id in self.merged_trans_dict:
            merged_trans_obj = self.merged_trans_dict[merged_trans_id]
            
            quality_percent = merged_trans_obj.calc_identity()
            coverage_percent = merged_trans_obj.calc_coverage()
            
            if high_quality_percent == -1.0:
                high_quality_percent = quality_percent
                low_quality_percent = quality_percent
                high_cov_percent = coverage_percent
                low_cov_percent = coverage_percent
            else:
                if high_quality_percent < quality_percent:
                    high_quality_percent = quality_percent
                if low_quality_percent > quality_percent:
                    low_quality_percent = quality_percent
                if high_cov_percent < coverage_percent:
                    high_cov_percent = coverage_percent
                if low_cov_percent > coverage_percent:
                    low_cov_percent = coverage_percent
                    
        high_quality_percent_str = str(round(high_quality_percent,2))
        low_quality_percent_str = str(round(low_quality_percent,2))
        high_cov_percent_str = str(round(high_cov_percent,2))
        low_cov_percent_str = str(round(low_cov_percent,2))
        
        
        trans_report_list.append(high_cov_percent_str)
        trans_report_list.append(low_cov_percent_str)
        trans_report_list.append(high_quality_percent_str)
        trans_report_list.append(low_quality_percent_str)
        
        
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

        collapse_sj_start_err_list_str = []
        collapse_sj_end_err_list_str = []
        for i in xrange(len(self.collapse_sj_start_err_list)):
            collapse_sj_start_err_list_str.append(str(self.collapse_sj_start_err_list[i]))
            collapse_sj_end_err_list_str.append(str(self.collapse_sj_end_err_list[i]))

        if self.strand == "+":
            collapse_sj_start_err_list_str = list(reversed(collapse_sj_start_err_list_str))
            collapse_sj_end_err_list_str = list(reversed(collapse_sj_end_err_list_str))

        collapse_sj_start_err_list_line = ",".join(collapse_sj_start_err_list_str)
        collapse_sj_end_err_list_line = ",".join(collapse_sj_end_err_list_str)

        trans_report_list.append(collapse_sj_start_err_list_line)
        trans_report_list.append(collapse_sj_end_err_list_line)


        collapse_error_nuc_list = self.collapse_start_error_nuc_list
        
        if len(collapse_error_nuc_list) > 1:

            if self.strand == "+":
                if collapse_error_nuc_list[-1] ==  "na":
                    collapse_error_nuc_list.pop(-1)
                    collapse_error_nuc_list = list(reversed(collapse_error_nuc_list))
                else:
                    print("Error with collapse_error_nuc_list")
                    sys.exit()
            elif self.strand == "-":
                if collapse_error_nuc_list[0] == "na":
                    collapse_error_nuc_list.pop(0)
                else:
                    print("Error with collapse_error_nuc_list")
                    sys.exit()
            else:
                print("Error with strand in format_trans_report_line")
                sys.exit

        ####################################################
        # create simple error line

#        all_sj_both_error_simple_list = []
#        if collapse_error_nuc_list[0] != "na":
#            for sj_error_nuc in collapse_error_nuc_list:
#                sj_pre_error_line = sj_error_nuc.split(">")[0]
#                sj_post_error_line = sj_error_nuc.split(">")[1]
#
#                sj_pre_error_split = sj_pre_error_line.split("_")
#                sj_post_error_split = sj_post_error_line.split("_")
#
#                sj_pre_error_simple_string,sj_post_error_simple_string = simple_sj_error(sj_pre_error_split, sj_post_error_split)
#
#                sj_both_error_simple_string = sj_pre_error_simple_string + ">" + sj_post_error_simple_string
#
#                all_sj_both_error_simple_list.append(sj_both_error_simple_string)
#
#            all_sj_both_error_simple_list_line = ";".join(all_sj_both_error_simple_list)
#        else:
#            all_sj_both_error_simple_list_line = "na"

        ####################################################

        collapse_error_nuc_list_line = ";".join(collapse_error_nuc_list)

        trans_report_list.append(collapse_error_nuc_list_line)

#        trans_report_list.append(all_sj_both_error_simple_list_line)

        trans_report_line = "\t".join(trans_report_list)
        
        return trans_report_line

#above this line are def's used for looping through sam file
####################################################################################################
####################################################################################################
####################################################################################################
#below this line are def's use for post sam pocessing



####################################################################################################

def fuzzy_match(coord1,coord2,diff_threshold): # use this to allow for fuzzy matches of splice junctions
    diff_num = 0
    match_flag = "none"
    #diff_threshold = 10 
    
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
    
    return match_flag,diff_num
    
####################################################################################################

## added 2020/07/26
# use this for speeding up transcript matching
def exact_match_capped(trans_obj,o_trans_obj,strand):

    exact_match_flag = "not_exact"

    a_trans_start = trans_obj.start_pos
    a_trans_end = trans_obj.end_pos

    b_trans_start = o_trans_obj.start_pos
    b_trans_end = o_trans_obj.end_pos

    this_start_diff = abs(a_trans_start - b_trans_start)
    this_end_diff = abs(a_trans_end - b_trans_end)

    if strand == "+":
        if this_start_diff <= fiveprime_threshold and this_end_diff <= threeprime_threshold:

            sj_hash_a = trans_obj.sj_hash
            sj_hash_b = o_trans_obj.sj_hash

            if sj_hash_a == sj_hash_b:
                exact_match_flag = "exact_match"
    elif strand == "-":
        if this_start_diff <= threeprime_threshold and this_end_diff <= fiveprime_threshold:

            sj_hash_a = trans_obj.sj_hash
            sj_hash_b = o_trans_obj.sj_hash

            if sj_hash_a == sj_hash_b:
                exact_match_flag = "exact_match"


    return exact_match_flag

####################################################################################################

## added 2020/07/26
# use this for speeding up transacript matching , NO CAP
def exact_match_nocap(trans_obj,o_trans_obj,strand):

    exact_match_flag = "not_exact"

    if strand == "+":
        a_trans_end = trans_obj.end_pos
        b_trans_end = o_trans_obj.end_pos
        this_end_diff = abs(a_trans_end - b_trans_end)


        if this_end_diff <= threeprime_threshold:
            sj_hash_a = trans_obj.sj_hash
            sj_hash_b = o_trans_obj.sj_hash

            if sj_hash_a == sj_hash_b:
                exact_match_flag = "exact_match"

    elif strand == "-":
        a_trans_start = trans_obj.start_pos
        b_trans_start = o_trans_obj.start_pos
        this_start_diff = abs(a_trans_start - b_trans_start)

        if this_start_diff <= threeprime_threshold:

            sj_hash_a = trans_obj.sj_hash
            sj_hash_b = o_trans_obj.sj_hash

            if sj_hash_a == sj_hash_b:
                exact_match_flag = "exact_match"

    return exact_match_flag


####################################################################################################

def compare_transcripts(trans_obj,o_trans_obj,fiveprime_cap_flag,strand): #use this to compare two transcripts
    diff_num_exon_flag = 0
    max_exon_num = 0
    min_exon_num = 0
    
    #exon_diff_threshold = 10
    #fiveprime_threshold = 10
    #threeprime_threshold = 10
    
    a_trans_id = trans_obj.cluster_id
    b_trans_id = o_trans_obj.cluster_id
    
    e_start_list = trans_obj.exon_start_list
    o_e_start_list = o_trans_obj.exon_start_list
    
    e_end_list = trans_obj.exon_end_list
    o_e_end_list = o_trans_obj.exon_end_list
    
    if len(e_start_list) != len(o_e_start_list):
        diff_num_exon_flag = len(e_start_list) - len(o_e_start_list)
    
    trans_comp_flag = "none"
    
    short_trans = "none"
    long_trans = "none"
    
    if fiveprime_cap_flag == "capped" and diff_num_exon_flag != 0: # if 5prime capped then should have same number of exons
        trans_comp_flag = "diff_transcripts"
        start_match_list = []
        start_diff_list = []
        end_match_list = []
        end_diff_list = []
        short_trans = "none"
        min_exon_num = 0
    else:

        #get max and min number of exons
        if len(e_start_list) > len(o_e_start_list):
            max_exon_num = len(e_start_list)
            min_exon_num = len(o_e_start_list)
            short_trans = o_trans_obj.cluster_id
        elif len(e_start_list) < len(o_e_start_list):
            max_exon_num = len(o_e_start_list)
            min_exon_num = len(e_start_list)
            short_trans = trans_obj.cluster_id
        elif len(e_start_list) == len(o_e_start_list):
            max_exon_num = len(o_e_start_list)
            min_exon_num = len(e_start_list)
            
            #do this for nocap libraries
            if strand == "+":
                if e_start_list[0] < o_e_start_list[0]:
                    short_trans = o_trans_obj.cluster_id
                    long_trans = trans_obj.cluster_id
                elif e_start_list[0] > o_e_start_list[0]:
                    short_trans = trans_obj.cluster_id
                    long_trans = o_trans_obj.cluster_id
                elif e_start_list[0] == o_e_start_list[0]:
                    short_trans = "same"
                    long_trans = "same"
                else:
                    print("Error with short and long trans identification")
                    print(trans_obj.cluster_id + " " + o_trans_obj.cluster_id)
                    sys.exit()
            elif strand == "-":
                if e_end_list[0] > o_e_end_list[0]:
                    short_trans = o_trans_obj.cluster_id
                    long_trans = trans_obj.cluster_id
                elif e_end_list[0] < o_e_end_list[0]:
                    short_trans = trans_obj.cluster_id
                    long_trans = o_trans_obj.cluster_id
                elif e_end_list[0] == o_e_end_list[0]:
                    short_trans = "same"
                    long_trans = "same"
                else:
                    print("Error with short and long trans identification")
                    print(trans_obj.cluster_id + " - " + o_trans_obj.cluster_id + " strand: " + strand)
                    sys.exit()
                    

        
        start_match_list = []
        start_diff_list = []
        end_match_list = []
        end_diff_list = []
        
        all_match_flag = 1 # 1 if all matching and 0 if at least one not matching
        
        for i in xrange(min_exon_num):
            
            if strand == "+":
                j = -1 * (i + 1) #iterate from last exon to account for possible 5' degradation for forward strand
            elif strand == "-":
                j = i # iterate from first exon for reverse strand
            
            
            e_start = e_start_list[j]
            o_e_start = o_e_start_list[j]
            e_end = e_end_list[j]
            o_e_end = o_e_end_list[j]

            # check for micro exons which do not overlap but fit in wobble range
            if e_start >= o_e_end:
                trans_comp_flag = "diff_transcripts"
                start_match_list = []
                start_diff_list = []
                end_match_list = []
                end_diff_list = []
                short_trans = "none"
                min_exon_num = 0
                continue

            if o_e_start >= e_end:
                trans_comp_flag = "diff_transcripts"
                start_match_list = []
                start_diff_list = []
                end_match_list = []
                end_diff_list = []
                short_trans = "none"
                min_exon_num = 0
                continue

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
            
            #use this for the condition that no 5 prime cap and shorter 5 prime end
            if fiveprime_cap_flag == "no_cap" and i == min_exon_num-1: # if no 5 cap and this is last 5 prime exon
                if strand == "+":
                    if end_match_flag == "no_match":
                        all_match_flag = 0
                    elif start_match_flag == "no_match" and all_match_flag == 1:
                        # check that shorter transcript has shorter 5' end
                        if b_trans_id == short_trans and start_diff_num < 0 and diff_num_exon_flag != 0:
                            trans_comp_flag = "same_three_prime_diff_exons"
                        elif  a_trans_id == short_trans and start_diff_num > 0 and diff_num_exon_flag != 0:
                            trans_comp_flag = "same_three_prime_diff_exons"
                        #elif short_trans == "same" and start_diff_num < 0 and diff_num_exon_flag == 0: # if short trans is same then same number of exons
                        elif start_diff_num < 0 and diff_num_exon_flag == 0:  # if short trans is same then same number of exons
                            short_trans = b_trans_id
                            trans_comp_flag = "same_three_prime_same_exons"
                        #elif short_trans == "same" and start_diff_num > 0 and diff_num_exon_flag == 0: # if short trans is same then same number of exons
                        elif start_diff_num > 0 and diff_num_exon_flag == 0:  # if short trans is same then same number of exons
                            short_trans = a_trans_id
                            trans_comp_flag = "same_three_prime_same_exons"
                        else:
                            all_match_flag = 0
                        #################################################################################################################
                        ########################## Continue here! 2017-06-13 Figure out how to differentiate shorter transcript!!!!! Need to define shorter/longer roles for trans grouping!! 
                        #################################################################################################################
                elif strand == "-":
                    if start_match_flag == "no_match":
                        all_match_flag = 0
                    elif end_match_flag == "no_match" and all_match_flag == 1:
                        # check that shorter transcript has shorter 5' end
                        if b_trans_id == short_trans and end_diff_num > 0 and diff_num_exon_flag != 0:
                            trans_comp_flag = "same_three_prime_diff_exons"
                        elif  a_trans_id == short_trans and end_diff_num < 0 and diff_num_exon_flag != 0:
                            trans_comp_flag = "same_three_prime_diff_exons"
                        #elif short_trans == "same" and end_diff_num > 0 and diff_num_exon_flag == 0: # if short trans is same then same number of exons
                        elif end_diff_num > 0 and diff_num_exon_flag == 0: # if short trans is same then same number of exons
                            short_trans = b_trans_id
                            trans_comp_flag = "same_three_prime_same_exons"
                        #elif short_trans == "same" and end_diff_num < 0 and diff_num_exon_flag == 0: # if short trans is same then same number of exons
                        elif end_diff_num < 0 and diff_num_exon_flag == 0:  # if short trans is same then same number of exons
                            short_trans = a_trans_id
                            trans_comp_flag = "same_three_prime_same_exons"
                        else:
                            all_match_flag = 0
                            
            # if capped or if nocap without being the last exon use normal matching procedure
            else:      
                if start_match_flag == "no_match":
                    all_match_flag = 0
                if end_match_flag == "no_match":
                    all_match_flag = 0
                    
            
            start_match_list.append(start_match_flag)
            start_diff_list.append(start_diff_num)
            end_match_list.append(end_match_flag)
            end_diff_list.append(end_diff_num)
            
            
            
        
        #cleanup trans_comp_flag
        if trans_comp_flag == "none":
            if all_match_flag == 1:        
                if diff_num_exon_flag == 0 and fiveprime_cap_flag == "capped":
                    trans_comp_flag = "same_transcript"
                elif diff_num_exon_flag == 0 and fiveprime_cap_flag == "no_cap":
                    trans_comp_flag = "same_three_prime_same_exons"
                elif diff_num_exon_flag != 0 and fiveprime_cap_flag == "no_cap":
                    trans_comp_flag = "same_three_prime_diff_exons"
                else:
                    print("Error with missing trans_comp_flag ")
                    print(trans_comp_flag)
                    print(diff_num_exon_flag)
                    print(max_exon_num)
                    print(min_exon_num)
                    print(start_match_list)
                    print(end_match_list)
                    sys.exit()
            else:
                trans_comp_flag = "diff_transcripts"
        
        
        
        
        #Keep in mind that the lists are ordered from 3' end to 5' end
        
#        if trans_obj.trans_id == '1_14_c81301/2/1370' and o_trans_obj.trans_id == '1_14_c147778/1/1455':
#            print("blah blah")
#            print(trans_comp_flag)
#            print(start_match_list)
#            print(end_match_list)
#            print(start_diff_list)
#            print(end_diff_list)
#            print(min_exon_num)
#            sys.exit()
#        if trans_obj.trans_id == '1_14_c147778/1/1455' and o_trans_obj.trans_id == '1_14_c81301/2/1370':
#            print("blah blah")
#            print(trans_comp_flag)
#           print(start_match_list)
#            print(end_match_list)
#            print(start_diff_list)
#            print(end_diff_list)
#            print(min_exon_num)
#            sys.exit()
        
    return trans_comp_flag,start_match_list,start_diff_list,end_match_list,end_diff_list,short_trans,long_trans,min_exon_num,diff_num_exon_flag

####################################################################################################


def sj_error_priority_start(sj_pre_error,sj_post_error):

    if sj_post_error == no_mismatch_flag and sj_pre_error == no_mismatch_flag:
        e_start_priority = 0
    elif sj_post_error == no_mismatch_flag and sj_pre_error != no_mismatch_flag:
        e_start_priority = 1
    elif sj_post_error != no_mismatch_flag and sj_pre_error == no_mismatch_flag:
        e_start_priority = 2
    elif sj_post_error != no_mismatch_flag and sj_pre_error != no_mismatch_flag:
        e_start_priority = 3
    else:
        print("Error with splice junction priority start")
        sys.exit()

    return e_start_priority


def sj_error_priority_end(sj_pre_error, sj_post_error):
    if sj_post_error == no_mismatch_flag and sj_pre_error == no_mismatch_flag:
        e_end_priority = 0
    elif sj_post_error != no_mismatch_flag and sj_pre_error == no_mismatch_flag:
        e_end_priority = 1
    elif sj_post_error == no_mismatch_flag and sj_pre_error != no_mismatch_flag:
        e_end_priority = 2
    elif sj_post_error != no_mismatch_flag and sj_pre_error != no_mismatch_flag:
        e_end_priority = 3
    else:
        print("Error with splice junction priority end")
        sys.exit()

    return e_end_priority


####################################################################################################
def sj_error_priority_finder(trans_obj,i,max_exon_num):
    #figure out errors near splice junctions
    ######################################################

    e_start_list = trans_obj.exon_start_list
    e_end_list = trans_obj.exon_end_list

    strand = trans_obj.strand
    exon_num = len(e_start_list)
    if max_exon_num == "none":
        max_exon_num = exon_num

    sj_pre_error_list = trans_obj.sj_pre_error_list
    sj_post_error_list = trans_obj.sj_post_error_list

    
    #use these to record the error type
    e_start_priority_error = "na"
    e_end_priority_error = "na"

    priority_error_delimit = ">"

    if len(e_start_list) > 1: # if there are splice junctions to cover

        if strand == "+":
            if i == 0:  # last 3' exon accoridng to genome
                e_end_priority = 0

                sj_pre_error = sj_pre_error_list[-1]
                sj_post_error = sj_post_error_list[-1]

                e_start_priority = sj_error_priority_start(sj_pre_error, sj_post_error)
                
                e_start_priority_error = sj_pre_error + priority_error_delimit + sj_post_error

            elif i > 0 and i < max_exon_num - 1: # middle exons

                sj_pre_error_start = sj_pre_error_list[-1 - i]
                sj_post_error_start = sj_post_error_list[-1 - i]

                sj_pre_error_end = sj_pre_error_list[-i]
                sj_post_error_end = sj_post_error_list[-i]

                e_start_priority = sj_error_priority_start(sj_pre_error_start, sj_post_error_start)
                e_end_priority = sj_error_priority_end(sj_pre_error_end, sj_post_error_end)
                
                e_start_priority_error = sj_pre_error_start + priority_error_delimit + sj_post_error_start
                e_end_priority_error = sj_pre_error_end + priority_error_delimit + sj_post_error_end
                


            elif i == max_exon_num - 1: # first 5' exon
                e_start_priority = 0
                sj_pre_error_end = sj_pre_error_list[-i]
                sj_post_error_end = sj_post_error_list[-i]

                e_end_priority = sj_error_priority_end(sj_pre_error_end, sj_post_error_end)
                
                e_end_priority_error = sj_pre_error_end + priority_error_delimit + sj_post_error_end
                

            else:
                print("Error with plus strand start and end priority")
                print(max_exon_num)
                print(i)
                sys.exit()


        elif strand == "-":
            if i == 0:  # 5' first exon according to genome
                e_start_priority = 0

                sj_post_error = sj_post_error_list[0]
                sj_pre_error = sj_pre_error_list[0]

                e_end_priority = sj_error_priority_end(sj_pre_error, sj_post_error)
                
                e_end_priority_error = sj_pre_error + priority_error_delimit + sj_post_error

            elif i > 0 and i < max_exon_num - 1: # middle exons

                sj_pre_error_start = sj_pre_error_list[i-1]
                sj_post_error_start = sj_post_error_list[i-1]

                sj_pre_error_end = sj_pre_error_list[i]
                sj_post_error_end = sj_post_error_list[i]

                e_start_priority = sj_error_priority_start(sj_pre_error_start, sj_post_error_start)
                e_end_priority = sj_error_priority_end(sj_pre_error_end, sj_post_error_end)
                
                e_start_priority_error = sj_pre_error_start + priority_error_delimit + sj_post_error_start
                e_end_priority_error = sj_pre_error_end + priority_error_delimit + sj_post_error_end


            elif i == max_exon_num - 1: # last 3' exon according to genome
                e_end_priority = 0
                sj_pre_error_start = sj_pre_error_list[i - 1]
                sj_post_error_start = sj_post_error_list[i - 1]

                e_start_priority = sj_error_priority_start(sj_pre_error_start, sj_post_error_start)
        
                e_start_priority_error = sj_pre_error_start + priority_error_delimit + sj_post_error_start

            else:
                print("Error with minus strand start and end priority")
                print(max_exon_num)
                print(i)
                sys.exit()

    else: # if 1 exon read
        e_start_priority = 0
        e_end_priority = 0

    # you can turn priority on or off with arguments
    if sj_priority_flag == "no_priority":
        e_start_priority = 0
        e_end_priority = 0
    elif sj_priority_flag != "sj_priority":
        print("Error with splice junction priority flag. Please use no_priority or sj_priority.")
        sys.exit()

    return e_start_priority,e_end_priority,e_start_priority_error,e_end_priority_error


####################################################################################################
def length_error_type(error_string,error_report):
    # 0>0.C.A;0>0;1I_5M>0;1I_7M>0;8.A.T_7.T.A>8M_1D;0>0

    if len(error_string.split(".")) == 3:
        error_type = "mismatch"
    elif error_string == "0":
        error_type = "0"
    elif "I" in error_string:
        error_type = "I"
    elif "D" in error_string:
        error_type = "D"
    elif "M" in error_string:
        error_type = "M"
    elif "S" in error_string:
        error_type = "S"
    elif "H" in error_string:
        error_type = "H"
    else:
        print("Error with error_string: error type not recognized")
        print(error_string)
        print(error_report)
        sys.exit()


    if error_type == "0":
        error_length = 0
    if error_type == "I":
    
        if len(error_string.split("I")) != 2:
            print("Error with error_string: Insertion line issue")
            print(error_string)
            print(error_report)
            sys.exit()
        
        error_length = int(error_string.split("I")[0])
    
    elif error_type == "D":
    
        if len(error_string.split("D")) != 2:
            print("Error with error_string: Deletion line issue")
            print(error_string)
            print(error_report)
            sys.exit()
        
        error_length = int(error_string.split("D")[0])
    
    elif error_type == "S":
    
        if len(error_string.split("S")) != 2:
            print("Error with error_string: Soft clipping line issue")
            print(error_string)
            print(error_report)
            sys.exit()
        
        error_length = int(error_string.split("S")[0])
    
    elif error_type == "H":
    
        if len(error_string.split("H")) != 2:
            print("Error with error_string: Hard clipping line issue")
            print(error_string)
            print(error_report)
            sys.exit()
        
        error_length = int(error_string.split("H")[0])
    
    elif error_type == "mismatch": # this is a mismatch error representing only one nt position
        error_length = 1
        
    elif error_type == "M": # This is a match and not an error
        if len(error_string.split("M")) != 2:
            print("Error with error_string: M line issue")
            print(error_string)
            print(error_report)
            sys.exit()
        
        error_length = int(error_string.split("M")[0])
    
    return error_length,error_type

def convert_int_list_to_string(int_list):
    
    str_list = []
    
    for int_val in int_list:
        str_list.append(str(int_val))
    
    return str_list


def simple_sj_error(sj_pre_error_split,sj_post_error_split):

    #ses_match_char = "_"

    ################################################################################
    ################################################################################

    sj_pre_error_simple_list = []
    sj_pre_error_count = 1

    sj_pre_error_split_reverse = sj_pre_error_split

    sj_pre_error_split_reverse.reverse()

    for sj_pre_error_string in sj_pre_error_split_reverse:
        # sj_pre_error_count += 1

        if len(sj_pre_error_string.split(".")) == 3:
            mismatch_position = int(sj_pre_error_string.split(".")[0])
            mismatch_position += 1
            if sj_pre_error_count == 1:
                for j in xrange(mismatch_position - 1):
                    sj_pre_error_simple_list.append(ses_match_char)
                    sj_pre_error_count += 1
                sj_pre_error_simple_list.append("X")
                sj_pre_error_count += 1
            elif sj_pre_error_count > 1:
                m_pos_diff = mismatch_position - sj_pre_error_count
                # for j in xrange(m_pos_diff-1):
                for j in xrange(m_pos_diff):
                    sj_pre_error_simple_list.append(ses_match_char)
                    sj_pre_error_count += 1
                sj_pre_error_simple_list.append("X")
                sj_pre_error_count += 1
        elif len(sj_pre_error_string.split(".")) == 1:

            if sj_pre_error_string == "0":
                for j in xrange(sj_err_threshold):
                    sj_pre_error_simple_list.append(ses_match_char)

            else:
                [cig_dig_list, cig_char_list] = cigar_list(sj_pre_error_string)

                cig_dig = int(cig_dig_list[0])
                cig_char = cig_char_list[0]

                if cig_char == "M":
                    for j in xrange(cig_dig):
                        sj_pre_error_simple_list.append(ses_match_char)
                        sj_pre_error_count += 1
                elif cig_char == "I":
                    for j in xrange(cig_dig):
                        sj_pre_error_simple_list.append("I")
                        sj_pre_error_count += 1
                elif cig_char == "D":
                    for j in xrange(cig_dig):
                        sj_pre_error_simple_list.append("D")
                        sj_pre_error_count += 1
                elif cig_char == "S":
                    for j in xrange(cig_dig):
                        sj_pre_error_simple_list.append("S")
                        sj_pre_error_count += 1
                elif cig_char == "H":
                    for j in xrange(cig_dig):
                        sj_pre_error_simple_list.append("H")
                        sj_pre_error_count += 1
        else:
            print("Error with LDE error char")
            sys.exit()

    if len(sj_pre_error_simple_list) < sj_err_threshold:
        for j in xrange(sj_err_threshold - len(sj_pre_error_simple_list)):
            sj_pre_error_simple_list.append(ses_match_char)

    sj_pre_error_simple_list_reverse = sj_pre_error_simple_list
    sj_pre_error_simple_list_reverse.reverse()

    sj_pre_error_simple_string = "".join(sj_pre_error_simple_list_reverse)

    ################################################################################
    ################################################################################

    sj_post_error_simple_list = []
    sj_post_error_count = 1
    for sj_post_error_string in sj_post_error_split:
        # sj_post_error_count += 1

        if len(sj_post_error_string.split(".")) == 3:
            mismatch_position = int(sj_post_error_string.split(".")[0])
            mismatch_position += 1
            if sj_post_error_count == 1:
                for j in xrange(mismatch_position - 1):
                    sj_post_error_simple_list.append(ses_match_char)
                    sj_post_error_count += 1
                sj_post_error_simple_list.append("X")
                sj_post_error_count += 1
            elif sj_post_error_count > 1:
                m_pos_diff = mismatch_position - sj_post_error_count
                # for j in xrange(m_pos_diff-1):
                for j in xrange(m_pos_diff):
                    sj_post_error_simple_list.append(ses_match_char)
                    sj_post_error_count += 1
                sj_post_error_simple_list.append("X")
                sj_post_error_count += 1
                #################################################################################################Continue here
        elif len(sj_post_error_string.split(".")) == 1:

            if sj_post_error_string == "0":
                for j in xrange(sj_err_threshold):
                    sj_post_error_simple_list.append(ses_match_char)

            else:
                [cig_dig_list, cig_char_list] = cigar_list(sj_post_error_string)

                cig_dig = int(cig_dig_list[0])
                cig_char = cig_char_list[0]

                if cig_char == "M":
                    for j in xrange(cig_dig):
                        sj_post_error_simple_list.append(ses_match_char)
                        sj_post_error_count += 1
                elif cig_char == "I":
                    for j in xrange(cig_dig):
                        sj_post_error_simple_list.append("I")
                        sj_post_error_count += 1
                elif cig_char == "D":
                    for j in xrange(cig_dig):
                        sj_post_error_simple_list.append("D")
                        sj_post_error_count += 1
                elif cig_char == "S":
                    for j in xrange(cig_dig):
                        sj_post_error_simple_list.append("S")
                        sj_post_error_count += 1
                elif cig_char == "H":
                    for j in xrange(cig_dig):
                        sj_post_error_simple_list.append("H")
                        sj_post_error_count += 1
        else:
            print("Error with LDE error char")
            sys.exit()

    if len(sj_post_error_simple_list) < sj_err_threshold:
        for j in xrange(sj_err_threshold - len(sj_post_error_simple_list)):
            sj_post_error_simple_list.append(ses_match_char)

    sj_post_error_simple_string = "".join(sj_post_error_simple_list)


    return sj_pre_error_simple_string,sj_post_error_simple_string


def sj_error_local_density(trans_obj):
    # figure out errors near splice junctions
    ######################################################

    e_start_list = trans_obj.exon_start_list
    e_end_list = trans_obj.exon_end_list

    strand = trans_obj.strand
    exon_num = len(e_start_list)
    
    trans_cigar = trans_obj.cigar

    max_exon_num = exon_num

    sj_pre_error_list = trans_obj.sj_pre_error_list
    sj_post_error_list = trans_obj.sj_post_error_list

    # use these to record the error type
    e_start_priority_error = "na"
    e_end_priority_error = "na"

    priority_error_delimit = ">"
    
    bad_sj_num_list = [] # list of the splice junctions that failed to pass quality threshold
    
    bad_sj_error_count_list = [] #  number of errors at each SJ
    
    bad_sj_num_pre_list = []
    
    bad_sj_num_post_list = []
    
    bad_sj_flag = 0 # if 0 then no bad SJ, if higher than 0 then there are bad SJ
    
    sj_error_list = []

    sj_error_nuc_list = []

    all_sj_post_error_simple_list = []
    all_sj_pre_error_simple_list = []

    all_sj_both_error_simple_list = []

    for i in xrange(max_exon_num-1):
        
        this_bad_sj_flag = 0

        sj_pre_error_i_count = 0
        sj_post_error_i_count = 0

        sj_pre_error_d_count = 0
        sj_post_error_d_count = 0

        sj_pre_error_m_count = 0
        sj_post_error_m_count = 0
        
        sj_pre_error_s_count = 0
        sj_post_error_s_count = 0
        
        sj_pre_error_h_count = 0
        sj_post_error_h_count = 0
        
        sj_pre_error_all_line = ""
        sj_post_error_all_line = ""

        sj_pre_error = sj_pre_error_list[i]
        sj_post_error = sj_post_error_list[i]

        sj_error_nuc_string = sj_pre_error + ">" + sj_post_error
        sj_error_nuc_list.append(sj_error_nuc_string)

        # 0>0.C.A;0>0;1I_5M>0;1I_7M>0;8.A.T_7.T.A>8M_1D;0>0
        
        sj_pre_error_split = sj_pre_error.split("_")
        sj_post_error_split = sj_post_error.split("_")

        sj_pre_error_simple_string,sj_post_error_simple_string = simple_sj_error(sj_pre_error_split, sj_post_error_split)

        all_sj_post_error_simple_list.append(sj_post_error_simple_string)
        all_sj_pre_error_simple_list.append(sj_pre_error_simple_string)

        all_sj_both_error_simple_list.append(sj_pre_error_simple_string + ">" +sj_post_error_simple_string)


        #print(e_start_list)
        #print(e_end_list)
        #print(sj_pre_error_list)
        #print(sj_post_error_list)

        for sj_pre_error_char in sj_pre_error_split:


            error_length,error_type = length_error_type(sj_pre_error_char,trans_cigar)
            #########################################################################################
            #########################################################################################
            #########################################################################################Continue here RK 2018/10/09

            if error_type == "mismatch":
                sj_pre_error_m_count = sj_pre_error_m_count + error_length
            elif error_type == "I":
                sj_pre_error_i_count = sj_pre_error_i_count + error_length
            elif error_type == "D":
                sj_pre_error_d_count = sj_pre_error_d_count + error_length
            elif error_type == "S":
                sj_pre_error_s_count = sj_pre_error_s_count + error_length
            elif error_type == "H":
                sj_pre_error_h_count = sj_pre_error_h_count + error_length

        
        for sj_post_error_char in sj_post_error_split:
            
            error_length,error_type = length_error_type(sj_post_error_char,trans_cigar)
            
            #####################################################Continue here RK 2018/10/09

            if error_type == "mismatch":
                sj_post_error_m_count = sj_post_error_m_count + error_length
            elif error_type == "I":
                sj_post_error_i_count = sj_post_error_i_count + error_length
            elif error_type == "D":
                sj_post_error_d_count = sj_post_error_d_count + error_length
            elif error_type == "S":
                sj_post_error_s_count = sj_post_error_s_count + error_length
            elif error_type == "H":
                sj_post_error_h_count = sj_post_error_h_count + error_length

        
        sj_pre_error_all_count = sj_pre_error_i_count + sj_pre_error_d_count + sj_pre_error_m_count + sj_pre_error_s_count + sj_pre_error_h_count
        sj_post_error_all_count = sj_post_error_i_count + sj_post_error_d_count + sj_post_error_m_count + sj_post_error_s_count + sj_post_error_h_count
        
        sj_pre_error_all_line = ",".join([str(sj_pre_error_i_count),str(sj_pre_error_d_count),str(sj_pre_error_m_count),str(sj_pre_error_s_count),str(sj_pre_error_h_count)])
        sj_post_error_all_line = ",".join([str(sj_post_error_i_count),str(sj_post_error_d_count),str(sj_post_error_m_count),str(sj_post_error_s_count),str(sj_post_error_h_count)])
        
        sj_all_error_all_line = ">".join([sj_pre_error_all_line,sj_post_error_all_line])
        
        sj_error_list.append(sj_all_error_all_line)
        
        
        if strand == "+":
            sj_num = i + 1
        elif strand == "-":
            sj_num = max_exon_num - i - 1
        
        
        if sj_pre_error_all_count > lde_threshold:
            bad_sj_flag += 1
            this_bad_sj_flag = 1
            
            bad_sj_num_pre_list.append(sj_num)
        
        if sj_post_error_all_count > lde_threshold:
            bad_sj_flag += 1
            this_bad_sj_flag = 1
            
            bad_sj_num_post_list.append(sj_num)
        
        if this_bad_sj_flag > 0 :
            bad_sj_num_list.append(sj_num)
        
        sj_error_string =  str(sj_pre_error_all_count) + ">" + str(sj_post_error_all_count)
        
        bad_sj_error_count_list.append(sj_error_string)




    sj_lde_flag = "na"
    if bad_sj_flag ==  0:
        sj_lde_flag = "lde_pass"
    elif bad_sj_flag >  0:
        sj_lde_flag = "lde_fail"
    else:
        print("Error with sj_lde_flag")
        sys.exit()

    #prepare lde outline
    if len(sj_error_list) > 0:
        sj_error_line = ";".join(sj_error_list)
    elif len(sj_error_list) == 0:
        sj_error_line = "na"
    else:
        print("Error with sj_error_line")
        sys.exit()

    bad_sj_num_str_list = convert_int_list_to_string(bad_sj_num_list)

    if len(bad_sj_num_str_list) > 0:
        bad_sj_num_line = ",".join(bad_sj_num_str_list)
    elif len(bad_sj_num_str_list) == 0:
        bad_sj_num_line = "0"
    else:
        print("Error with bad_sj_num_line")
        sys.exit()

    if len(bad_sj_error_count_list) > 0 :
        bad_sj_error_count_line = ",".join(bad_sj_error_count_list)
    elif len(bad_sj_error_count_list) == 0 :
        bad_sj_error_count_line = "na"
    else:
        print("Error with bad_sj_error_count_line")
        sys.exit()

    num_exon_str = str(max_exon_num)
    
    #####################################################
    #generate detailed error profile around splice junctions

    if len(sj_error_nuc_list) > 0:
        collapse_error_nuc_list_line = ";".join(sj_error_nuc_list)
    else:
        collapse_error_nuc_list_line = "na"

    if len(all_sj_both_error_simple_list) > 0:
        all_sj_both_error_simple_line = ";".join(all_sj_both_error_simple_list)
    else:
        all_sj_both_error_simple_line = "na"

    #####################################################
    
    lde_outlist = []
    lde_outlist.append(trans_obj.cluster_id)
    lde_outlist.append(sj_lde_flag)
    lde_outlist.append(trans_obj.scaff_name)
    lde_outlist.append(str(trans_obj.start_pos))
    lde_outlist.append(str(trans_obj.end_pos))
    lde_outlist.append(trans_obj.strand)
    lde_outlist.append(num_exon_str)
    lde_outlist.append(bad_sj_num_line)
    lde_outlist.append(bad_sj_error_count_line)
    lde_outlist.append(sj_error_line)
    lde_outlist.append(collapse_error_nuc_list_line)
    lde_outlist.append(all_sj_both_error_simple_line)
    lde_outlist.append(trans_obj.cigar)
    
    
    
    lde_outline = "\t".join(lde_outlist)


    return bad_sj_flag,bad_sj_num_list,bad_sj_num_pre_list,bad_sj_num_post_list,bad_sj_error_count_list,lde_outline

##############################################################


def collapse_transcripts(trans_obj_list,fiveprime_cap_flag,collapse_flag): #use this to collapse transcripts
    # all supplied transcripts will be merged
    # create supplied transcripts by comparing with compare transcripts
    
    try:
        collapse_flag
    except NameError:
        print("collapse_flag not defined, using default of most commond ends")
        collapse_flag == "common_ends"
    
    max_exon_num = 0
    strand = "none"
    num_trans = len(trans_obj_list)
    collapse_trans_id_list = []
    
    #get strand and max exon num
    for trans_obj in trans_obj_list:
        collapse_trans_id_list.append(trans_obj.trans_id)
        
        e_start_list = trans_obj.exon_start_list
        if strand == "none":
            strand = trans_obj.strand
        elif trans_obj.strand != strand:
            print("mismatch in strand from trans_obj_list for collapsing def")
            sys.exit()
        exon_num = len(e_start_list)
        if exon_num > max_exon_num:
            max_exon_num = exon_num
    
    collapse_start_list = []
    collapse_end_list = []
    
    #use these to record actualy nuc mismatch info for collapsed models
    collapse_start_error_nuc_list = []
    collapse_end_error_nuc_list = []

    #use these to keep track of errors near splice junctions
    collapse_sj_start_err_list = []
    collapse_sj_end_err_list = []

    
    
    #track how much wobble for the starts and end in the collapse
    start_wobble_list = []
    end_wobble_list = []
    for i in xrange(max_exon_num): #go from 3 prime end
        if strand == "+":
            j = -1 * (i + 1) #iterate from last exon to account for possible 5' degradation for forward strand
        elif strand == "-":
            j = i # iterate from first exon for reverse strand

        e_start_dict = {}  # e_start_dict[priority number][start] = number of occurrences
        e_end_dict = {}  # e_end_dict[priority number][end] = number of occurrences
        
        priority_error_start_dict = {} #priority_error_start_dict[priority number][start][error string] = 1
        priority_error_end_dict = {} #priority_error_end_dict[priority number][end][error string] = 1


        #e_start_dict = {} # e_start_dict[start] = number of occurrences
        #e_end_dict = {} # e_end_dict[end] = number of occurrences

        #use these to get wobble despite sj priority algorithm
        e_start_range_list = []
        e_end_range_list = []

        for trans_obj in trans_obj_list:
            e_start_list = trans_obj.exon_start_list
            e_end_list = trans_obj.exon_end_list


            # figure out errors near splice junctions
            ######################################################
            sj_pre_error_list = trans_obj.sj_pre_error_list
            sj_post_error_list = trans_obj.sj_post_error_list

            this_max_exon_num = len(e_start_list)

            if i >= len(e_start_list):# use for no cap when exon numbers may not match
                continue

            e_start_priority, e_end_priority, e_start_priority_error,e_end_priority_error = sj_error_priority_finder(trans_obj, i, this_max_exon_num) ####################################


            e_start = int(e_start_list[j])
            e_end = int(e_end_list[j])

            #Use this to prevent using truncated 5' end of trans without max exon num
            if fiveprime_cap_flag == "no_cap":####################################################################
                #do not use 5' end if this is not priority level one for start
                if i == len(e_start_list)-1 and i < max_exon_num-1:
                    if strand == "+":
                        e_start = -1
                    elif strand == "-":
                        e_end = -1
            
            if e_start != -1: # check that it isnt truncated 5' end
                if e_start_priority not in e_start_dict:
                    e_start_dict[e_start_priority] = {}
                    priority_error_start_dict[e_start_priority] = {} #####################################################
                if e_start not in e_start_dict[e_start_priority]:
                    e_start_dict[e_start_priority][e_start] = 0
                    priority_error_start_dict[e_start_priority][e_start] = {} #####################################################

                e_start_range_list.append(e_start) #use these to get wobble despite sj priority algorithm

                e_start_dict[e_start_priority][e_start] += 1
                
                priority_error_start_dict[e_start_priority][e_start][e_start_priority_error] = 1
            else:
                if log_flag == "log_on":
                    print("truncated transcript")
                #e_start_dict[e_start_priority] = {}
                
                #priority_error_start_dict[e_start_priority] = {}


            
            if e_end != -1: # check that it isnt truncated 5' end
                if e_end_priority not in e_end_dict:
                    e_end_dict[e_end_priority] = {}
                    priority_error_end_dict[e_end_priority] = {}#####################################################
                if e_end not in e_end_dict[e_end_priority]:
                    e_end_dict[e_end_priority][e_end] = 0
                    priority_error_end_dict[e_end_priority][e_end] = {} #####################################################

                e_end_range_list.append(e_end) #use these to get wobble despite sj priority algorithm

                e_end_dict[e_end_priority][e_end] += 1
                
                priority_error_end_dict[e_end_priority][e_end][e_end_priority_error] = 1
            else:
                if log_flag == "log_on":
                    print("truncated transcript")
                #e_end_dict[e_end_priority] = {}
                
                #priority_error_end_dict[e_end_priority] = {}

        ##########################################

        best_e_start = -1
        long_e_start = -1
        short_e_start = -1
        num_starts = 0

        priority_start_list = list(e_start_dict.keys())
        priority_end_list = list(e_end_dict.keys())
        
        priority_start_list.sort()
        priority_end_list.sort()

        best_start_priority = priority_start_list[0]
        best_end_priority = priority_end_list[0]

        collapse_sj_start_err_list.append(best_start_priority)
        collapse_sj_end_err_list.append(best_end_priority)

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
                if log_flag == "log_on":
                    print("more than one best e start! " + str(best_e_start) + " num_trans: " + str(num_trans))
        ##########################################

        e_start_range_list.sort()

        e_start_wobble = e_start_range_list[-1] - e_start_range_list[0]
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
                if log_flag == "log_on":
                    print("more than one best e end! " + str(best_e_end) + " num_trans: " + str(num_trans))
        ##########################################

        e_end_range_list.sort()

        e_end_wobble = e_end_range_list[-1] - e_end_range_list[0]
        end_wobble_list.append(e_end_wobble)
        
        if fiveprime_cap_flag == "no_cap" and i+1 == max_exon_num: # use earliest 5' end for no cap libraries
            if strand == "+":
                best_e_start = long_e_start
            elif strand == "-":
                best_e_end = long_e_end
        
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
        
        if best_e_start > best_e_end:
            print("Error with collapsing, e_start bigger than e_end")
            print(best_e_start)
            print(best_e_end)
            print(trans_obj.trans_id)
            print(trans_obj.strand)
            print(e_start_dict)
            print(e_end_dict)
            for trans_obj in trans_obj_list:
                print(str(trans_obj.trans_id) + " ##########################")
                print(trans_obj.exon_start_list)
                print(trans_obj.exon_end_list)

            if long_e_start < long_e_end:
                best_e_start = long_e_start
                best_e_end = long_e_end
            else:
                print("Error long_e_start is greater than long_e_end")
                print(long_e_start)
                print(long_e_end)

                sys.exit()
        
        collapse_start_list.append(best_e_start)
        collapse_end_list.append(best_e_end)
        
        ##########################################################
        ###debugging
#        print("debugging ==================")
#        print(best_e_start)
#        print(best_e_end)
#        print(best_start_priority)
#        print(long_e_start)
#        print(most_long_e_start)
#        print(i)
#        print(max_exon_num)
#        print(priority_error_start_dict)
#        print(strand)
#        print(trans_obj_list)
#        print("debugging ==================")
        
        ##########################################################
        priority_error_start_list = list(priority_error_start_dict[best_start_priority][best_e_start].keys()) ###########################################################################
        priority_error_end_list = list(priority_error_end_dict[best_end_priority][best_e_end].keys())
        
        priority_error_start_line = "-".join(priority_error_start_list)
        priority_error_end_line = "-".join(priority_error_end_list)
        
        collapse_start_error_nuc_list.append(priority_error_start_line)
        collapse_end_error_nuc_list.append(priority_error_end_line)

    #put the coords in the right order maintaining order with wobble lists
    collapse_start_list, start_wobble_list = zip(*sorted(zip(collapse_start_list, start_wobble_list)))
    collapse_end_list, end_wobble_list = zip(*sorted(zip(collapse_end_list, end_wobble_list)))
    
    collapse_start_list = list(collapse_start_list)
    start_wobble_list = list(start_wobble_list)
    collapse_end_list = list(collapse_end_list)
    end_wobble_list = list(end_wobble_list)
    
    #######################################################
    #Below: check start and end list to make sure there are no overlapping coordinates
    prev_start = -1
    prev_end = -1
    for i in xrange(len(collapse_start_list)):
        check_start = collapse_start_list[i]
        check_end = collapse_end_list[i]
        
        if check_end <= check_start: # exon end must always be greater than exon start
            print("Error with exon end earlier than exon start")
            print(str(check_start) + "\t" + str(check_end))
            print(collapse_trans_id_list)
            sys.exit()
        
        if check_start <= prev_start: # next start should always be later than prev start
            print("Error with this exon start not later than previous start")
            print(str(prev_start) + "\t" + str(check_start))
            print(collapse_trans_id_list)
            sys.exit()
        
        if check_end <= prev_end: # next start should always be later than prev start
            print("Error with this exon end not later than previous end")
            print(str(prev_end) + "\t" + str(check_end))
            print(collapse_trans_id_list)
            sys.exit()
        
        prev_start = check_start
        prev_end = check_end
    # Above: check start and end list to make sure there are no overlapping coordinates
    #######################################################      

    #flip order of sj err list if the strand is positive since we started from the 3' end.
#    if strand == "+":
#        collapse_sj_start_err_list = list(reversed(collapse_sj_start_err_list))
#        collapse_sj_end_err_list = list(reversed(collapse_sj_end_err_list))

        #collapse_start_error_nuc_list = list(reversed(collapse_start_error_nuc_list))
        #collapse_end_error_nuc_list = list(reversed(collapse_end_error_nuc_list))

#    elif strand == "-":
#        collapse_sj_start_err_list = collapse_sj_start_err_list
#        collapse_sj_end_err_list = collapse_sj_end_err_list
#    else:
#        print("Issue with strand value in collapse transcripts")
#        print(strand)
#        sys.exit()


    
    return collapse_start_list,collapse_end_list,start_wobble_list,end_wobble_list,collapse_sj_start_err_list,collapse_sj_end_err_list,collapse_start_error_nuc_list,collapse_end_error_nuc_list

####################################################################################################

def gene_group(trans_list): #groups trans into genes, does not take into account strand
    
    trans_obj_list = []
    
    for trans_id in trans_list:
        trans_obj_list.append(trans_obj_dict[trans_id])

    gene_trans_dict = {} # gene_trans_dict[gene id][trans id] = 1
    trans_gene_dict = {} # trans_gene_dict[trans id] = gene group
    gene_start_dict = {} # gene_start_dict[gene num] = gene start
    start_gene_dict = {} # start_gene_dict[start] = gene num

    gene_count = 0
    
    if len(trans_list) == 1:
        gene_count += 1
        trans_id = trans_list[0]
        single_gene_start = trans_obj_list[0].exon_start_list[0]
        gene_start_dict[gene_count] = single_gene_start
        gene_trans_dict[gene_count] = {}
        gene_trans_dict[gene_count][trans_id] = 1
        trans_gene_dict[trans_id] = gene_count
        
    
    for i in xrange(len(trans_obj_list)):
        trans_obj = trans_obj_list[i]
        for j in xrange(i+1,len(trans_obj_list)):
            o_trans_obj = trans_obj_list[j]

            trans_id  = trans_obj.cluster_id
            o_trans_id = o_trans_obj.cluster_id
            
            #if trans_id == o_trans_id:#skip if same
            #    continue
            
            if trans_id in trans_gene_dict and o_trans_id in trans_gene_dict: # skip if already in same group
                if trans_gene_dict[trans_id] == trans_gene_dict[o_trans_id]:
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

                
                if trans_id not in trans_gene_dict: #if no gene groups make new one
                    gene_count += 1
                    trans_gene_dict[trans_id] = gene_count

                    gene_trans_dict[gene_count] = {}
                    gene_trans_dict[gene_count][trans_id] = 1

                    #add gene start
                    gene_start_dict[gene_count] = exon_start_list[0]
                
                if o_trans_id not in trans_gene_dict: #if no gene groups make new one
                    gene_count += 1
                    trans_gene_dict[o_trans_id] = gene_count

                    gene_trans_dict[gene_count] = {}
                    gene_trans_dict[gene_count][o_trans_id] = 1

                    #add gene start
                    gene_start_dict[gene_count] = o_exon_start_list[0]
                
                    
                
            if overlap_flag == 1:
                if trans_id not in trans_gene_dict and o_trans_id not in trans_gene_dict: #if no gene groups make new one
                    gene_count += 1
                    trans_gene_dict[trans_id] = gene_count
                    trans_gene_dict[o_trans_id] = gene_count
                    gene_trans_dict[gene_count] = {}
                    gene_trans_dict[gene_count][trans_id] = 1
                    gene_trans_dict[gene_count][o_trans_id] = 1
                    #add gene start
                    min_gene_start = exon_start_list[0]
                    if min_gene_start > o_exon_start_list[0]:
                        min_gene_start = o_exon_start_list[0]
                    gene_start_dict[gene_count] = min_gene_start
                    
                    
                elif trans_id not in trans_gene_dict: # add to other gene group

                    gene_num = trans_gene_dict[o_trans_id]
                    trans_gene_dict[trans_id] = gene_num
                    gene_trans_dict[gene_num][trans_id] = 1
                    
                    min_gene_start = exon_start_list[0]
                    if min_gene_start > o_exon_start_list[0]:
                        min_gene_start = o_exon_start_list[0]
                    gene_start_dict[gene_num] = min_gene_start
                elif o_trans_id not in trans_gene_dict:# add to other gene group
                    
                    gene_num = trans_gene_dict[trans_id]
                    trans_gene_dict[o_trans_id] = gene_num
                    gene_trans_dict[gene_num][o_trans_id] = 1
                    
                    min_gene_start = exon_start_list[0]
                    if min_gene_start > o_exon_start_list[0]:
                        min_gene_start = o_exon_start_list[0]
                    gene_start_dict[gene_num] = min_gene_start
                elif trans_id in trans_gene_dict and o_trans_id in trans_gene_dict:

                    gene_num = trans_gene_dict[trans_id]
                    o_gene_num = trans_gene_dict[o_trans_id]
                    
                    if gene_num != o_gene_num: #merge gene groups
                        m_trans_id_list = list(gene_trans_dict[o_gene_num].keys())
                        for m_trans_id in m_trans_id_list:
                            
                            trans_gene_dict[m_trans_id] = gene_num
                            gene_trans_dict[gene_num][m_trans_id] = 1
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
        for trans_id in gene_trans_dict[gene_num]:
            gene_start_trans_dict[gene_start][trans_id] = 1
    
        
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

                [list_slice, sort_flag] = iterate_sort_list(list_slice, pos_index + 1)

                list_trans_pos_list[min_index:max_index] = list_slice

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


####################################################################################################


def sort_transcripts(trans_obj_list):
    #sort transcripts by start-end-exon starts

    pos_trans_dict = {} # pos_trans_dict[pos] = trans obj
    pos_trans_list = []
    
    for trans_obj in trans_obj_list:
        trans_scaff = trans_obj.scaff_name
        trans_exon_start_list = trans_obj.exon_start_list
        trans_exon_end_list = trans_obj.exon_end_list
        trans_exon_start_list.sort()
        trans_exon_end_list.sort()
        trans_start = trans_exon_start_list[0]
        trans_end = trans_exon_end_list[-1]
        
        trans_pos_list = []

        trans_pos_list.append(str(trans_start))
        trans_pos_list.append(",")
        trans_pos_list.append(str(trans_end))
        trans_pos_list.append(",")

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

        dup_detect_flag = 0 # use for signaling that a duplicate has been detected

        #if this trans model is already present in the dict then this is a duplicate
        if trans_pos_line in pos_trans_dict:

            dup_detect_flag = 1

            old_merge_id = pos_trans_dict[trans_pos_line].trans_id
            new_merge_id = trans_obj.trans_id
            
            #old_trans_list = pos_trans_dict[trans_pos_line].trans_list
            #new_trans_list = trans_obj.trans_list

            if log_flag == "log_on":

                print("Duplicate transcript positions in transcript sorting!")
                print(trans_obj.merged_trans_dict.keys())
                print(str(trans_start)+" "+str(trans_end))
                print(pos_trans_dict[trans_pos_line].merged_trans_dict.keys())
            this_bed_line = trans_obj.format_bed_line()
            other_bed_line = pos_trans_dict[trans_pos_line].format_bed_line()

            if log_flag == "log_on":
                print(this_bed_line)
                print(other_bed_line)
                print("a###########################################")
            for a_uniq_trans_id in trans_obj.merged_trans_dict:
                a_bed_line = trans_obj_dict[a_uniq_trans_id].format_bed_line(a_uniq_trans_id)
                if log_flag == "log_on":
                    print(a_bed_line)
            if log_flag == "log_on":
                print("b###########################################")
            for b_uniq_trans_id in pos_trans_dict[trans_pos_line].merged_trans_dict:
                b_bed_line = trans_obj_dict[b_uniq_trans_id].format_bed_line(b_uniq_trans_id)
                if log_flag == "log_on":
                    print(b_bed_line)
            if log_flag == "log_on":
                print("end duplicate###########################################")

            #sys.exit()
            
            #################################################################################
            #################################################################################
            if duplicate_flag == "no_merge":
                print("By default TAMA collapse does not allow merging of duplicate transcript groups.")
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

                #collapse_start_list, collapse_end_list, start_wobble_list, end_wobble_list, e_start_trans_dict, e_end_trans_dict = collapse_transcripts(match_trans_obj_list, collapse_flag)
                collapse_start_list,collapse_end_list,start_wobble_list,end_wobble_list,collapse_sj_start_err_list,collapse_sj_end_err_list,collapse_start_error_nuc_list,collapse_end_error_nuc_list = collapse_transcripts(match_trans_obj_list,fiveprime_cap_flag,collapse_flag)

                # update merge info
                #pos_trans_dict[trans_pos_line].add_merge_info(collapse_start_list, collapse_end_list, start_wobble_list,end_wobble_list, e_start_trans_dict, e_end_trans_dict)
                pos_trans_dict[trans_pos_line].add_merge_info(collapse_start_list,collapse_end_list,start_wobble_list,end_wobble_list,collapse_sj_start_err_list,collapse_sj_end_err_list,collapse_start_error_nuc_list,collapse_end_error_nuc_list )

                trans_obj = pos_trans_dict[trans_pos_line]

            else:
                print("Error with duplicate transcript group flag")
                sys.exit()
            #################################################################################
            #################################################################################

        if dup_detect_flag == 0:
            pos_trans_dict[trans_pos_line] = trans_obj
            pos_trans_list.append(trans_pos_line)
        elif dup_detect_flag == 1:
            pos_trans_dict[trans_pos_line] = trans_obj
        else:
            print("Error with dup_detect_flag.")
            sys.exit()

    [new_pos_trans_list,new_pos_trans_dict] = sort_pos_trans_list(pos_trans_list, pos_trans_dict)
    
    sorted_trans_obj_list = []
    for pos_trans in new_pos_trans_list:
        
        trans_obj = new_pos_trans_dict[pos_trans]
        sorted_trans_obj_list.append(trans_obj)
        tmp_id = trans_obj.trans_id

    
    # sorted_trans_obj_list list of trana obj that have been sorted by position
    return sorted_trans_obj_list

##############################################################################
def longest_transcript(trans_id,trans_id_list):
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
# Use this class to manage group merging
class TransGroup:
    def __init__(self, group_name):
        self.group_name = group_name 
        self.trans_group_dict = {}  #trans_group_dict[trans id][trans group] = "longest" or "short" for no cap
        self.group_trans_dict = {}  # group_trans_dict[group num][trans id] = "longest" or "short" for no cap
        self.group_count = 0
        self.group_longest_dict = {} #group_longest_dict[group num][longest/long][trans] = 1

        self.group_max_exon_dict = {} # group_max_exon_dict[group num] = max exon num
    
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
        #checks the number of associated trans to a nocap group
        assoc_trans_dict = {} # assoc_trans_dict[trans id] = 1
        
        for a_trans_group in self.trans_group_dict[trans_a]:
            for trans_b in self.group_trans_dict[a_trans_group]:
                if trans_b != trans_a:
                    assoc_trans_dict[trans_b] = 1
        
        num_assoc_trans = len(list(assoc_trans_dict.keys()))
        
        return num_assoc_trans
    
    def delete_trans(self,trans_a):
        
        #delete trans from all its associated groups
        for a_trans_group in self.trans_group_dict[trans_a]:
            self.group_trans_dict[a_trans_group].pop(trans_a,None)
            
            if len(list(self.group_trans_dict[a_trans_group])) == 0:# if ther group is now empty, delete it
                self.group_trans_dict.pop(a_trans_group,None)
        
        #delete trans from trans group dict
        self.trans_group_dict.pop(trans_a,None)
    
    def new_group_a(self,trans_a):

        if log_flag == "log_on":
            print("invoke new_group_a " + str(self.group_count) + " newgroup " + trans_a)
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
        
##############################################################################

    def add_a_to_b_group(self,trans_a,trans_b):
        if log_flag == "log_on":
            print("invoke add_a_to_b_group " + trans_a + " " + trans_b)
        # only cap libs should be used for b group and they only have one group
        # does not take all of a_group just uses a_trans
        if len(list(self.trans_group_dict[trans_b].keys())) > 1:
            if log_flag == "log_on":
                print("multiple groups")
                #sys.exit()

        #check that trans b has more exons than trans a
        trans_obj_a = trans_obj_dict[trans_a]
        trans_obj_b = trans_obj_dict[trans_b]
        
        #if trans_obj_a.num_exons >= trans_obj_b.num_exons: ######################################## 2019/06/07
        #    print("Error trans_a does not have fewer exons than trans_b")
        #    print(trans_a + " " + trans_b)
        #    print(str(trans_obj_a.num_exons) + " " + str(trans_obj_b.num_exons))
        #    sys.exit()
        
        #remove initial nocap group that is a self identity group
        if len(list(self.trans_group_dict[trans_a].keys())) == 1:# if only one group
            a_trans_group = list(self.trans_group_dict[trans_a].keys())[0]
            if log_flag == "log_on":
                print(trans_a)
            if len(list(self.group_trans_dict[a_trans_group].keys())) == 1: #if only this in one group
                if list(self.group_trans_dict[a_trans_group].keys())[0] == trans_a:
                    self.group_trans_dict.pop(a_trans_group,None)
                    self.trans_group_dict.pop(trans_a,None)

        if trans_a not in self.trans_group_dict:
            self.trans_group_dict[trans_a] = {}
            
        else:
            # this happens if trans_a is a nocap trans in which case it can be in multiple groups
            if log_flag == "log_on":
                print("trans_a already in group, should be nocap trans: " + trans_a)
        

        for b_group_num in list(self.trans_group_dict[trans_b].keys()):
            # add a trans to b group
            self.trans_group_dict[trans_a][b_group_num] = "short"
            self.group_trans_dict[b_group_num][trans_a] = "short"
            # a_trans has fewer exons than b_trans so it must be short
                
   
            ## search through a groups for group mergings
            ## if a is longest in any of it's groups then you can add other trans from those groups to b
            #for a_group_num in list(self.trans_group_dict[trans_a].keys()):
            #    if a_group_num == b_group_num:  ######################################## 2019/06/07
            #        continue
            #
            #    if log_flag == "log_on":
            #        print(str(a_group_num) + "-a and b group num-" + str(b_group_num))
            #    # if trans_a is the longest in group
            #    #add all shorter transcripts from a group to b group too
            #    if trans_a in self.group_longest_dict[a_group_num]["longest"]:
            #        a_trans_id_list = list(self.group_trans_dict[a_group_num].keys())
            #        for a_group_trans in a_trans_id_list:
            #            self.trans_group_dict[a_group_trans][b_group_num] = "short"
            #            self.group_trans_dict[b_group_num][a_group_trans] = "short"
            #            self.trans_group_dict[a_group_trans].pop(a_group_num,None)
            #            self.group_trans_dict.pop(a_group_num,None)
            #            self.group_longest_dict.pop(a_group_num,None)
        
            # dont need to redo longest and short because added a trans is all short compared to b group


##############################################################################

    def merge_a_b_groups(self,trans_a,trans_b):
        if log_flag == "log_on":
            print("invoke merge_a_b_groups "+ str(self.group_count )+ " " + trans_a + " " +trans_b )
        
        #self.group_count += 1
        #only cap lib trans should be used for merging groups
        if len(list(self.trans_group_dict[trans_a].keys())) > 1:
            print("multiple groups a")
            sys.exit()
        if len(list(self.trans_group_dict[trans_b].keys())) > 1:
            print("multiple groups b")
            sys.exit()
            
        a_group_num = list(self.trans_group_dict[trans_a].keys())[0]
        b_group_num = list(self.trans_group_dict[trans_b].keys())[0]
        
        if a_group_num == b_group_num:
            print("Error, groups are the same, no need to merge!")
            sys.exit()
        
        #if self.group_count in self.group_trans_dict:
        #    print("group num already used")
        #    sys.exit()
        
        #find bigger group
        num_trans_group_a = len(list(self.group_trans_dict[a_group_num].keys()))
        num_trans_group_b = len(list(self.group_trans_dict[b_group_num].keys()))
        
        merge_group_num = -1
        if num_trans_group_a > num_trans_group_b:
            merge_group_num = a_group_num
            self.group_max_exon_dict[merge_group_num] = self.group_max_exon_dict[a_group_num]
            
            for group_trans in self.group_trans_dict[b_group_num]:
                self.group_trans_dict[merge_group_num][group_trans] = 1
                self.trans_group_dict[group_trans].pop(b_group_num, None)
                self.trans_group_dict[group_trans][merge_group_num] = 1
            
            #remove old group
            self.group_trans_dict.pop(b_group_num, None)

        elif num_trans_group_b >= num_trans_group_a:
            merge_group_num = b_group_num
            self.group_max_exon_dict[merge_group_num] = self.group_max_exon_dict[b_group_num]
            
            for group_trans in self.group_trans_dict[a_group_num]:
                self.group_trans_dict[merge_group_num][group_trans] = 1
                self.trans_group_dict[group_trans].pop(a_group_num, None)
                self.trans_group_dict[group_trans][merge_group_num] = 1
            #remove old group
            self.group_trans_dict.pop(a_group_num, None)
            
        else:
            print("Error with comparing group trans counts")
            sys.exit()
        
        #add longest trans information
        #trans_id_list = list(self.group_trans_dict[merge_group_num].keys())
        
        #redo longest trans flags
        #for trans_c in self.group_trans_dict[merge_group_num]:
        #    longest_trans_flag = longest_transcript(trans_c,trans_id_list)
        #    self.trans_group_dict[trans_c][merge_group_num] = longest_trans_flag
        #    self.group_trans_dict[merge_group_num][trans_c] = longest_trans_flag
        
        if merge_group_num == -1:
            print("Error with merge groups, merge_group_num == -1")
            sys.exit()



##############################################################################
 
    def merge_a_b_groups_nocap(self,trans_a,trans_b):
        # need to create new group and cant use time saving thing in capped merge
        # this is because nocaps can have multiple groups
        self.group_count += 1

        if log_flag == "log_on":
            print("invoke merge_a_b_groups_nocap " + str(self.group_count )+ " " + trans_a + " " +trans_b )
        #only cap lib trans should be used for merging groups
        if len(list(self.trans_group_dict[trans_a].keys())) > 1:
            if log_flag == "log_on":
                print("multiple groups a nocap")

        if len(list(self.trans_group_dict[trans_b].keys())) > 1:
            if log_flag == "log_on":
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
                longest_trans_flag = longest_transcript(check_trans_id,[longest_trans_id])
                if longest_trans_flag == "longest":
                    longest_trans_id = check_trans_id
            
            #refresh trans flags in new group
            self.group_longest_dict[self.group_count] = {}
            self.group_longest_dict[self.group_count]["longest"] = {}
            
            for trans_c in trans_id_list:
                longest_trans_flag = longest_transcript(trans_c,[longest_trans_id])
                self.trans_group_dict[trans_c][self.group_count] = longest_trans_flag
                self.group_trans_dict[self.group_count][trans_c] = longest_trans_flag
                
                if longest_trans_flag == "longest":
                    self.group_longest_dict[self.group_count]["longest"][trans_c] = 1
        else: # no new group was made, merge did not happen, 2 short transcripts
            self.group_trans_dict.pop(self.group_count, None)
        
        ###################################### <<< edge of work


        

####################################################################################################

def simplify_gene_capped(trans_obj_list,fiveprime_cap_flag): # goes through transcripts in gene and groups transcripts for collapsing
    #for capped only!!
    if log_flag == "log_on":
        print("invoking simplify_gene_capped")
    
    transgroup = TransGroup("transgroup")
    
    #trans_group_dict = {} # trans_group_dict[trans id] = trans group
    #group_trans_dict = {} # group_trans_dict[group num][trans id] = 1
    
    #group_num = 0

    #new cluster grouping algorithm
    ############################################################################################################
    ###########################
    # convert trans_obj_list to trans_obj_dict for while looping

    # clusters that have yet been grouped
    ungrouped_trans_obj_dict = {}  # ungrouped_trans_obj_dict[cluster_id] =  trans_obj
    # clusters that have been grouped
    grouped_trans_obj_dict = {}  # grouped_trans_obj_dict[cluster_id] =  trans_obj
    # grouped clusters that have been not been searched/been used as a hunter
    unsearched_trans_obj_dict = {}  # unsearched_trans_obj_dict[cluster_id] =  trans_obj

    #############
    # create SJ Hash for transcript models if there are many reads 2020/07/27
    num_group_reads = len(trans_obj_list)
    #if num_group_reads > sj_hash_read_threshold - 1:
    #    for trans_obj in trans_obj_list:
    #        #trans_obj.make_sj_hash_int()
    #        trans_obj.make_sj_hash_string()
    #############

    for trans_obj in trans_obj_list:
        ungrouped_trans_obj_dict[trans_obj.cluster_id] = trans_obj


    ###########################
    
    #all_trans_id_dict = {} # all_trans_id_dict[trans id] = 1

    ungrouped_count = len(trans_obj_list)
    unsearched_count = 0

    while ungrouped_count > 0:

        if unsearched_count == 0:
            ungrouped_cluster_list = list(ungrouped_trans_obj_dict.keys())
            ungrouped_cluster_list.sort()

            # hunter trans and prey trans, hunter used to look for prey
            hunter_cluster_id = ungrouped_cluster_list[0]
            hunter_trans_obj = ungrouped_trans_obj_dict[hunter_cluster_id]
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

            #all_trans_id_dict[hunter_cluster_id] = i
            hunter_strand = hunter_trans_obj.strand
            hunter_num_exons = hunter_trans_obj.num_exons

            a_group_check = transgroup.check_trans_status(hunter_cluster_id)
            # make groups for each transcript if no group
            if a_group_check != 1:
                transgroup.new_group_a(hunter_cluster_id)

            for prey_cluster_id in ungrouped_trans_obj_dict:
                prey_trans_obj = ungrouped_trans_obj_dict[prey_cluster_id]

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
                    transgroup.new_group_a(prey_cluster_id)

                group_match_flag = transgroup.check_same_group(hunter_cluster_id, prey_cluster_id)

                # this shoudn't be needed anymore due to new dict based group search
                if group_match_flag == 1:  # if they are both in the same group
                    continue


                trans_comp_flag, start_match_list, start_diff_list, end_match_list, end_diff_list, short_trans, long_trans, min_exon_num, diff_num_exon_flag = compare_transcripts(hunter_trans_obj, prey_trans_obj, fiveprime_cap_flag, hunter_strand)

                ######added to speed things up


                #if num_group_reads < sj_hash_read_threshold or hunter_num_exons < 3:
                #    trans_comp_flag, start_match_list, start_diff_list, end_match_list, end_diff_list, short_trans, long_trans, min_exon_num, diff_num_exon_flag = compare_transcripts(hunter_trans_obj, prey_trans_obj, fiveprime_cap_flag, hunter_strand)
                #
                #else:
                #
                #    exact_match_flag = exact_match_capped(hunter_trans_obj, prey_trans_obj,hunter_strand)
                #
                #    if exact_match_flag == "exact_match":
                #        trans_comp_flag = "same_transcript"
                #
                #    elif hunter_num_exons == 1: # exact match is the same as compare for single exon transcripts
                #        trans_comp_flag = 'diff_transcripts'
                #
                #    elif exact_match_flag == "not_exact":
                #        trans_comp_flag, start_match_list, start_diff_list, end_match_list, end_diff_list, short_trans, long_trans, min_exon_num, diff_num_exon_flag = compare_transcripts(hunter_trans_obj, prey_trans_obj, fiveprime_cap_flag, hunter_strand)
                #
                #    else:
                #        print("Error with exact match output")
                #        sys.exit()



                # old system used to generalize for nocap mode but this slows down capped mode
                #trans_match_flag = 0
                # same_transcript means clusters should be grouped for collapsing!
                #if trans_comp_flag == "same_transcript":
                #    trans_match_flag = 1

                #a_group_check = transgroup.check_trans_status(hunter_cluster_id)
                #b_group_check = transgroup.check_trans_status(prey_cluster_id)

                ##########################################Affects all downstream code!
                if trans_comp_flag != "same_transcript":  # skip if there is no match
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



    ############################################################################################################
    # End of new cluster grouping algorithm


    trans_group_dict = transgroup.trans_group_dict
    group_trans_dict = transgroup.group_trans_dict

  
    return trans_group_dict,group_trans_dict


#####################################################################################################

def simplify_gene_nocap(trans_obj_list,fiveprime_cap_flag): # goes through transcripts in gene and groups transcripts for collapsing
    # For nocap only!
    if log_flag == "log_on":
        print("invoking simplify_gene_nocap")
    
    transgroup = TransGroup("transgroup")

    #############
    # create SJ Hash for transcript models if there are many reads 2020/07/27
    num_group_reads = len(trans_obj_list)
    #if num_group_reads > sj_hash_read_threshold - 1:
    #    for trans_obj in trans_obj_list:
    #        #trans_obj.make_sj_hash_int()
    #        trans_obj.make_sj_hash_string()
    #############

    #new cluster grouping algorithm
    ############################################################################################################
    ###########################
    # convert trans_obj_list to trans_obj_dict for while looping

    # clusters that have yet been grouped
    ungrouped_trans_obj_dict = {}  # ungrouped_trans_obj_dict[cluster_id] =  trans_obj
    # clusters that have been grouped
    grouped_trans_obj_dict = {}  # grouped_trans_obj_dict[cluster_id] =  trans_obj
    # grouped clusters that have been not been searched/been used as a hunter
    unsearched_trans_obj_dict = {}  # unsearched_trans_obj_dict[cluster_id] =  trans_obj
    # use this dict to organize clusters by num exons
    exon_trans_obj_dict = {} # exon_trans_obj_dict[num exons][cluster id] = trans_obj
    # use this to refresh sub exon dicts
    sub_exon_cluster_dict = {} # sub_exon_cluster_dict[num exons][cluster id] = 1

    for trans_obj in trans_obj_list:
        #ungrouped_trans_obj_dict[trans_obj.cluster_id] = trans_obj
        
        this_strand = trans_obj.strand
        this_num_exons = trans_obj.num_exons
        
        if this_num_exons not in exon_trans_obj_dict:
            exon_trans_obj_dict[this_num_exons] = {}
            sub_exon_cluster_dict[this_num_exons] = {}
            
        
        exon_trans_obj_dict[this_num_exons][trans_obj.cluster_id] = trans_obj
        sub_exon_cluster_dict[this_num_exons][trans_obj.cluster_id] = 1


    ###########################

    #all_trans_id_dict = {} # all_trans_id_dict[trans id] = 1

    ungrouped_count = len(trans_obj_list)
    unsearched_count = 0
    
    exon_num_list = list(exon_trans_obj_dict.keys())
    exon_num_list.sort(reverse=True)

    max_num_exons = exon_num_list[0]
    
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
    
                    unsearched_trans_obj_dict.pop(hunter_cluster_id)

                # print("hunter: " + hunter_cluster_id )
                #print(exon_num_list)
    
                #all_trans_id_dict[hunter_cluster_id] = i
                hunter_strand = hunter_trans_obj.strand
                hunter_num_exons = hunter_trans_obj.num_exons
    
                a_group_check = transgroup.check_trans_status(hunter_cluster_id)
                # make groups for each transcript if no group
                if a_group_check != 1:
                    transgroup.new_group_a(hunter_cluster_id)
    
                ungrouped_trans_obj_list = list(ungrouped_trans_obj_dict.keys())
                unsearched_cluster_list = list(unsearched_trans_obj_dict.keys())
    
                # search at same exon num level
                ###############################################################################
                if hunter_cluster_id not in sub_length_cluster_dict: # only search at same level if this has not been grouped with longer transcript
                    for prey_cluster_id in ungrouped_trans_obj_dict:
        
                        prey_trans_obj = ungrouped_trans_obj_dict[prey_cluster_id]
    
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
                            transgroup.new_group_a(prey_cluster_id)
        
                        group_match_flag = transgroup.check_same_group(hunter_cluster_id, prey_cluster_id)
        
                        # this shoudn't be needed anymore due to new dict based group search
                        if group_match_flag == 1:  # if they are both in the same group
                            continue

                        trans_comp_flag, start_match_list, start_diff_list, end_match_list, end_diff_list, short_trans, long_trans, min_exon_num, diff_num_exon_flag = compare_transcripts(hunter_trans_obj, prey_trans_obj, fiveprime_cap_flag, hunter_strand)

                        # this is for same number of exon matches
                        #if num_group_reads < sj_hash_read_threshold or hunter_num_exons < 3:
                        #    trans_comp_flag, start_match_list, start_diff_list, end_match_list, end_diff_list, short_trans, long_trans, min_exon_num, diff_num_exon_flag = compare_transcripts(hunter_trans_obj, prey_trans_obj, fiveprime_cap_flag, hunter_strand)
                        #else:
                        #    exact_match_flag = exact_match_nocap(hunter_trans_obj, prey_trans_obj,hunter_strand)
                        #
                        #    if exact_match_flag == "exact_match":
                        #        trans_comp_flag = "same_three_prime_same_exons"  # Use this for long-long nocap comparison
                        #    elif hunter_num_exons == 1:
                        #        trans_comp_flag = 'diff_transcripts'
                        #    elif exact_match_flag == "not_exact":
                        #        trans_comp_flag, start_match_list, start_diff_list, end_match_list, end_diff_list, short_trans, long_trans, min_exon_num, diff_num_exon_flag = compare_transcripts(hunter_trans_obj, prey_trans_obj, fiveprime_cap_flag, hunter_strand)
                        #
                        #    else:
                        #        print("Error with exact match output")
                        #        sys.exit()


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
                                    transgroup.merge_a_b_groups_nocap(hunter_cluster_id,prey_cluster_id)
                                elif trans_comp_flag == "same_three_prime_same_exons" : ###################### 2019/06/07
                                    transgroup.add_a_to_b_group(prey_cluster_id,hunter_cluster_id)
                                else:
                                    print("Error with match flag")
                                    print(trans_comp_flag)
                                    sys.exit()

                                #print(transgroup.group_trans_dict)
        
                            #    elif hunter_num_exons > prey_num_exons: #add shorter to longer
                            #        transgroup.add_a_to_b_group(prey_cluster_id,hunter_cluster_id)
                            #    elif prey_num_exons > hunter_num_exons: # add shorter to longer
                            #        transgroup.add_a_to_b_group(hunter_cluster_id,prey_cluster_id)

                                # remove the prey cluster from dict to avoid redundant searching
                                # ungrouped_trans_obj_dict.pop(prey_cluster_id) # remove this outside of for loop
                                # add prey to unsearched dict
                                unsearched_trans_obj_dict[prey_cluster_id] = prey_trans_obj
        
                    # remove grouped prey from ungrouped dict
                    #for unsearched_cluster_id in unsearched_trans_obj_dict: #2019_06_07
                    #    if unsearched_cluster_id in ungrouped_trans_obj_dict: #2019_06_07
                    #        ungrouped_trans_obj_dict.pop(unsearched_cluster_id) #2019_06_07 
        
                    
                
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
                            transgroup.new_group_a(prey_cluster_id)
        
                        group_match_flag = transgroup.check_same_group(hunter_cluster_id, prey_cluster_id)
        
                        # this shoudn't be needed anymore due to new dict based group search
                        if group_match_flag == 1:  # if they are both in the same group
                            continue
        
                        trans_comp_flag, start_match_list, start_diff_list, end_match_list, end_diff_list, short_trans, long_trans, min_exon_num, diff_num_exon_flag = compare_transcripts(hunter_trans_obj, prey_trans_obj, fiveprime_cap_flag, hunter_strand)
        
                        #For nocap only!!!!
                        trans_match_flag = 0
                        if trans_comp_flag == "same_transcript":
                            #trans_match_flag = 1
                            print("Error with subgroup seach same_transcript")
                            sys.exit()
                        elif trans_comp_flag == "same_three_prime_same_exons" :
                            #trans_match_flag = 1
                            print("Error with subgroup seach same_three_prime_same_exons")
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
                                transgroup.add_a_to_b_group(prey_cluster_id,hunter_cluster_id)
                                sub_length_cluster_dict[prey_cluster_id] = 1
                                
                                #this_sub_exon_cluster_dict[prey_exon_level].pop(prey_cluster_id)  ######### 2019/06/06
                                
                            elif prey_num_exons > hunter_num_exons: # add shorter to longer
                                #transgroup.add_a_to_b_group(hunter_cluster_id,prey_cluster_id)
                                print("Error with subgroup exon match: hunter smaller")
                                sys.exit()

                # reset hunter id
                hunter_cluster_id = "new_hunter"
                unsearched_count = len(unsearched_trans_obj_dict)
                ungrouped_count = len(ungrouped_trans_obj_dict)
    
        exon_num_index += 1
                ################################################################################



    ############################################################################################################
    # End of new cluster grouping algorithm

    trans_group_dict = transgroup.trans_group_dict
    group_trans_dict = transgroup.group_trans_dict
    

    return trans_group_dict,group_trans_dict

####################################################################################################
####################################################################################################
####################################################################################################
def reverse_complement(seq_list):
    comp_dict = {}
    comp_dict["A"] = "T"
    comp_dict["T"] = "A"
    comp_dict["G"] = "C"
    comp_dict["C"] = "G"
    comp_dict["N"] = "N"

    comp_dict["R"] = "N"
    comp_dict["Y"] = "N"
    comp_dict["K"] = "N"

    comp_dict["M"] = "N"
    comp_dict["S"] = "N"
    comp_dict["W"] = "N"

    comp_dict["B"] = "N"
    comp_dict["D"] = "N"
    comp_dict["H"] = "N"
    comp_dict["V"] = "N"

    
    reverse_seq = seq_list[::-1] # that's a neat way of reversing the string
    
    rev_comp_list = []
    for base in reverse_seq:
        rev_comp_list.append(comp_dict[base])
    
    
    #rev_comp_string = "".join(rev_comp_list)
    
    return rev_comp_list

def detect_polya(trans_obj,a_window): # looks for a stretch of poly A in the genome region 3' of the end of the transcript
    strand = trans_obj.strand
    scaffold = trans_obj.scaff_name 

    if strand == "+":
        trans_end = trans_obj.end_pos
        
        downstream_seq = fasta_dict[scaffold][trans_end:trans_end+a_window]
        dseq_length = len(downstream_seq)
        if dseq_length == 0:
            if log_flag == "log_on":
                print("dseq_length == 0")
                print(trans_obj.trans_id)
                print(scaffold + " " + str(trans_obj.start_pos) + " " +  str(trans_obj.end_pos) + " " + strand)
            dseq_length = 1
        a_count = downstream_seq.count("A")
        n_count = downstream_seq.count("N")
        a_percent = float(a_count) / float(dseq_length)
        n_percent = float(n_count) / float(dseq_length)
    elif strand == "-":
        trans_end = trans_obj.start_pos
        a_window_start = trans_end-a_window
        if a_window_start < 0:
            if log_flag == "log_on":
                print("Window start less than 0")
                print(trans_obj.trans_id)
                print(scaffold + " " + str(trans_obj.start_pos) + " " +  str(trans_obj.end_pos)+ " " + strand)
            a_window_start = 0
        downstream_seq = fasta_dict[scaffold][a_window_start:trans_end]
        rev_comp_seq = reverse_complement(downstream_seq)
        downstream_seq = rev_comp_seq
        
        dseq_length = len(downstream_seq)
        if dseq_length == 0:

            if log_flag == "log_on":
                print("dseq_length == 0")
                print(trans_obj.trans_id)
                print(scaffold + " " + str(trans_obj.start_pos) + " " +  str(trans_obj.end_pos)+ " " + strand)
            dseq_length = 1
        
        a_count = downstream_seq.count("A")
        n_count = downstream_seq.count("N")
        a_percent = float(a_count) / float(dseq_length)
        n_percent = float(n_count) / float(dseq_length)
        
        #a_count = downstream_seq.count("T")
        #n_count = downstream_seq.count("N")
        #a_percent = float(a_count) / float(dseq_length)
        #n_percent = float(n_count) / float(dseq_length)
    else:
        print("Error with strand information for poly a detection")
        sys.exit()
    
    return downstream_seq,dseq_length,a_count,n_count,a_percent,n_percent
    
        


def detect_rt_switch(trans_obj): # looks for complementary structure in intronic region that might cause rt switch
    rt_window = 20 #window to search
    bind_length = 8 #kmer used to represent binding
    
    strand = trans_obj.strand
    scaffold = trans_obj.scaff_name
    
    this_exon_starts = trans_obj.exon_start_list
    this_exon_ends = trans_obj.exon_end_list
    
    num_junct_bind = 0
    
    bind_seq_dict = {} # bind_seq_dict[splice junction]['end seq'/'start seq'] = seq
    
    for i in xrange(len(this_exon_starts)-1):
        
        bind_flag = 0
        start_index = i + 1
        end_index = i
        end_seq = fasta_dict[scaffold][end_index:end_index+rt_window]
        start_seq = fasta_dict[scaffold][start_index-rt_window:start_index]
        
        rev_comp_end_seq = reverse_complement(end_seq)
        
        binding_dict = {} # binding_dict[bind seq] = 1
        for j in xrange(rt_window-bind_length):
            bind_seq = start_seq[j:j+bind_length]
            bind_seq_string = "".join(bind_seq)
            binding_dict[bind_seq_string] = 1
        
        for j in xrange(rt_window-bind_length):
            bind_seq = rev_comp_end_seq[j:j+bind_length]
            bind_seq_string = "".join(bind_seq)
            if bind_seq_string in binding_dict:
                bind_flag += 1
                
        
        if bind_flag > 0:
            num_junct_bind += 1
            if i not in bind_seq_dict:
                bind_seq_dict[i] = {}
            bind_seq_dict[i][1] = end_seq
            bind_seq_dict[i][2] = start_seq
            bind_seq_dict[i][3] = rev_comp_end_seq
    
    return num_junct_bind,bind_seq_dict

####################################################################################################

def compare_multimaps(trans_obj_a,trans_obj_b): ### Added this 2019/03/04
    
    best_map_id = "na"
    a_percent_cov = trans_obj_a.percent_cov
    b_percent_cov = trans_obj_b.percent_cov
    
    a_percent_identity = trans_obj_a.percent_identity
    b_percent_identity = trans_obj_b.percent_identity
    
    a_pass_flag = 0
    b_pass_flag = 0
    
    if a_percent_cov > coverage_threshold and a_percent_identity > identity_threshold:
        a_pass_flag = 1
    
    if b_percent_cov > coverage_threshold  and b_percent_identity > identity_threshold:
        b_pass_flag = 1
    
    if a_pass_flag == 0 and b_pass_flag == 0:
        best_trans_obj = trans_obj_a # both will not pass thresholds so just choose A
        best_map_id = "A"
    elif a_pass_flag > 0 and b_pass_flag == 0:
        best_trans_obj = trans_obj_a # B does not pass thresholds
        best_map_id = "A"
    elif a_pass_flag == 0 and b_pass_flag > 0:
        best_trans_obj = trans_obj_b # A does not pass thresholds
        best_map_id = "B"
    elif a_pass_flag > 0 and b_pass_flag > 0:
    
        # just use coverage to choose best mapping
        if a_percent_cov >= b_percent_cov:
            best_trans_obj = trans_obj_a
            best_map_id = "A"
        elif a_percent_cov < b_percent_cov:
            best_trans_obj = trans_obj_b
            best_map_id = "B"
    else:
        print("Error with compare_multimaps")
        print("Early termination of TAMA collapse run!!")
        sys.exit()
    
    return best_trans_obj,best_map_id
        
    
    
    

####################################################################################################
####################################################################################################
####################################################################################################

####################################################################################### Loop through fasta file

fasta_dict = {} # fasta_dict[scaffold name] = array for seq
fasta_header_dict = {} # fasta_header_dict[scaffold name] = fasta header
fasta_scaffold_list = [] # list of fatsa seq names to be compared to SAM file header

prev_time = track_time(start_time,prev_time)
#Create fasta lookup dict
print("going through fasta")
for seq_record in SeqIO.parse(fasta_file_name, "fasta"):
    seq_name = str(seq_record.id)
    seq_desc = str(seq_record.description)
    
    seq_string = str(seq_record.seq)
    seq_string = seq_string.upper()
    seq_length = len(seq_string)
    
    fasta_dict[seq_name] = list(seq_string)
    
    fasta_header_dict[seq_name] = seq_desc
    fasta_scaffold_list.append(seq_name)


sam_flag_dict = {} #sam_flag_dict[flag number] = meaning
sam_flag_dict[0] = "forward_strand"
sam_flag_dict[4] = "unmapped"
sam_flag_dict[16] = "reverse_strand"
sam_flag_dict[2048] = "chimeric"
sam_flag_dict[2064] = "chimeric"
sam_flag_dict[256] = "not_primary"
sam_flag_dict[272] = "not_primary"

unmapped_dict = {} # unmapped_dict[cluster id] = 1

sam_scaffold_list = []
sam_scaffold_dict = {} # sam_scaffold_dict[seq name] = seq length

####################################################################################################


########################################################################### loop through sam file
trans_obj_dict = {} # trans_obj_dict[cluster id] = trans obj

group_trans_list_dict = {} # group_trans_list_dict[group id] = list of trans
trans_group_dict = {} # trans_group_dict[trans id] = group id

this_scaffold = "none"
group_start_pos = 0
group_end_pos = 0
group_count = 0

scaffold_list = []

sam_count = 0

prev_time = track_time(start_time,prev_time)


#################################################################################################### SAM or BAM


print("going through sam file")

if bam_flag == "BAM":
    
    from subprocess import Popen, PIPE
    sam_file_contents = []

    samtools_path = "samtools"
    pline = [samtools_path, 'view', sam_file]
    try:
        p = Popen(pline, bufsize=-1, stdout=PIPE, stderr=PIPE)
    except OSError:
        raise OSError('Samtools not found!\n')

    sam_file_list = p.communicate()
    sam_file_contents = sam_file_list[0].split("\n")
    
    print(len(sam_file_contents))

elif bam_flag == "SAM":
    sam_file_obj = open(sam_file)
    sam_file_contents = sam_file_obj.read().rstrip("\n").split("\n")

##########################


############################################################################################################
############################################################################################################
############################################################################################################
# original mode start ##############################################

if run_mode_flag == "original":
    for line in sam_file_contents:

        #if sam_count == 0:
        #    print(line)

        line_split = line.split("\t")

    #    if line_split[0] == "@SQ":
    #        seq_name = line_split[1].split(":")[1]
    #        seq_length = line_split[2].split(":")[1]
    #
    #        sam_scaffold_dict[seq_name] = seq_length
    #        sam_scaffold_list.append(seq_name)
    #
        if line.startswith("@"):
            continue

        if line == "":
            continue

        sam_count += 1
        if sam_count % 5000 == 0:
            print("sam count " + str(sam_count))

        read_id = line_split[0]
        sam_flag = int(line_split[1])
        scaff_name = line_split[2]
        start_pos = int(line_split[3])

        cigar = line_split[5]
        read_seq = line_split[9]
        seq_list = list(read_seq)
        mapped_flag = sam_flag_dict[sam_flag]

        ####################################
        #Check sam and gmap strand info!!!
        #get strand information from gmap flag

        xs_flag = "na"
        for field in line_split:
            if "XS:A:" in field:
                xs_flag = field.split(":")[-1]

        if mapped_flag == "forward_strand"  and xs_flag == "-":
            outline_strand = "\t".join([read_id,scaff_name,str(start_pos),cigar,"+-"])
            outfile_strand.write(outline_strand)
            outfile_strand.write("\n")
        elif mapped_flag == "reverse_strand" and xs_flag == "+":
            outline_strand = "\t".join([read_id,scaff_name,str(start_pos),cigar,"-+"])
            outfile_strand.write(outline_strand)
            outfile_strand.write("\n")

        #
        # Above: Check sam and gmap strand info!!!
        ####################################

        if mapped_flag == "unmapped" or mapped_flag == "not_primary" or mapped_flag == "chimeric" :
            unmapped_dict[read_id] = 1
            accept_flag = mapped_flag # added this 2019/03/04
            percent_coverage = "NA"
            percent_identity = "NA"
            error_line = "NA"
            quality_percent = "NA"
            length = "NA"
            strand = "NA"
            cigar = "NA"

            cluster_line = "\t".join([read_id,mapped_flag,accept_flag,percent_coverage,percent_identity,error_line, length, cigar])
            outfile_cluster.write(cluster_line)
            outfile_cluster.write("\n")
            continue


        map_seq_length = mapped_seq_length(cigar)

        [end_pos,exon_start_list,exon_end_list] = trans_coordinates(start_pos,cigar)

        [h_count,s_count,i_count,d_count,mis_count,nomatch_dict,sj_pre_error_list,sj_post_error_list] =  calc_error_rate(start_pos,cigar,seq_list,scaff_name,read_id)

        trans_obj = Transcript(read_id)
        trans_obj.add_sam_info(sam_flag,scaff_name,start_pos,cigar,read_seq,seq_list)
        trans_obj.add_map_seq_length(map_seq_length)
        trans_obj.add_exon_coords(end_pos,exon_start_list,exon_end_list)
        trans_obj.add_mismatch(h_count,s_count,i_count,d_count,mis_count,nomatch_dict,sj_pre_error_list,sj_post_error_list)

        ##### 2020/07/27 sj hash
        #trans_obj.make_sj_hash_string()
        ##### 2020/07/27 sj hash


        percent_coverage = trans_obj.calc_coverage()
        percent_identity = trans_obj.calc_identity()

        percent_coverage_str = str(round(percent_coverage,2))
        percent_identity_str = str(round(percent_identity,2))

        error_line = trans_obj.make_error_line()
        seq_length = trans_obj.seq_length
        strand = trans_obj.strand

        multimap_flag = 0

        if percent_coverage < coverage_threshold or percent_identity < identity_threshold:
            accept_flag = "discarded"
            cluster_line = "\t".join([read_id,mapped_flag,accept_flag,percent_coverage_str,percent_identity_str,error_line, str(seq_length), cigar])
            outfile_cluster.write(cluster_line)
            outfile_cluster.write("\n")
            #skip the transcript because the mapping is poor
            continue

        else:

            bad_sj_flag,bad_sj_num_list,bad_sj_num_pre_list,bad_sj_num_post_list,bad_sj_error_count_list,lde_outline = sj_error_local_density(trans_obj)

            outfile_lde.write(lde_outline)
            outfile_lde.write("\n")

            if bad_sj_flag > 0:
                #sadfsdfs

                #bad_sj_num_str_list =convert_int_list_to_string(bad_sj_num_list)

                #bad_sj_num_line = ",".join(bad_sj_num_str_list)
                #bad_sj_error_count_line = ",".join(bad_sj_error_count_list)

                #lde_file_line = "\t".join([trans_obj.cluster_id,trans_obj.scaff_name,str(trans_obj.start_pos),str(trans_obj.end_pos),trans_obj.strand,bad_sj_num_line,bad_sj_error_count_line,trans_obj.cigar])
                #outfile_lde.write(lde_outline)
                #outfile_lde.write("\n")

                #######################################

                #add to cluster file

                accept_flag = "local_density_error"
                cluster_line = "\t".join([read_id,mapped_flag,accept_flag,percent_coverage_str,percent_identity_str,error_line, str(seq_length), cigar])
                outfile_cluster.write(cluster_line)
                outfile_cluster.write("\n")

                continue

            accept_flag = "accepted"
            cluster_line = "\t".join([read_id,mapped_flag,accept_flag,percent_coverage_str,percent_identity_str,error_line, str(seq_length), cigar])
            outfile_cluster.write(cluster_line)
            outfile_cluster.write("\n")

            #only run poly detetcion on accepted transcripts
            downstream_seq,dseq_length,a_count,n_count,a_percent,n_percent = detect_polya(trans_obj,a_window)
            trans_obj.add_polya_info(downstream_seq,dseq_length,a_count,n_count,a_percent,n_percent)

            #check for multi maps
            if read_id in trans_obj_dict:
                print("Read has multi map")
                print(line)
                print(percent_coverage)
                print(percent_identity)

                trans_obj_a = trans_obj_dict[read_id]
                trans_obj_b = trans_obj

                best_trans_obj,best_map_id = compare_multimaps(trans_obj_a,trans_obj_b)

                #only re-assign if the new map is better, otherwise old map is aready processed
                if best_map_id == "B":
                    trans_obj_dict[read_id] = best_trans_obj
                    multimap_flag = 1
                else:
                    # if this new map is not going to be used we can skip the rest of the loop
                    continue

            else:
                trans_obj_dict[read_id] = trans_obj


        #check if a read has multi mapped!
        # remove old read info if new map is better
        if multimap_flag == 1:
            old_group_count = trans_group_dict[read_id]

            #check that the old group is not only made up of this read mapping
            if len(group_trans_list_dict[old_group_count]) > 1:
                group_trans_list_dict[old_group_count].remove(read_id) # remove read from old group
                trans_group_dict.pop(read_id, None) # remove read from trans_group_dict, will be re-assigned later

            elif len(group_trans_list_dict[old_group_count]) == 1: # group is only made of this mapping
                group_trans_list_dict.pop(old_group_count, None) # remove group
                trans_group_dict.pop(read_id, None) # remove read from trans_group_dict, will be re-assigned later
            else:
                print("Error with dealing with multimap management")
                print("Warning: temrinated early!")
                sys.exit()

        #if read_id in trans_group_dict:
        #    print("cluster multi mapped!")
        #    print(line)
        #    print(percent_coverage)
        #    print(percent_identity)
        #    print(trans_obj.h_count)
        #    sys.exit()

        #group trans by start and end coords
        if this_scaffold == "none":
            this_scaffold = scaff_name
            group_start_pos = start_pos
            group_end_pos = end_pos

            group_trans_list_dict[group_count] = []
            group_trans_list_dict[group_count].append(read_id)
            trans_group_dict[read_id] = group_count

            scaffold_list.append(this_scaffold)

            continue

        if scaff_name == this_scaffold:

            if start_pos >= group_start_pos and start_pos <= group_end_pos: #add to group
                group_trans_list_dict[group_count].append(read_id)
                trans_group_dict[read_id] = group_count
                #update group end position
                if end_pos > group_end_pos:
                    group_end_pos = end_pos

            elif start_pos > group_end_pos: #start new group
                group_count += 1

                group_start_pos = start_pos
                group_end_pos = end_pos

                group_trans_list_dict[group_count] = []
                group_trans_list_dict[group_count].append(read_id)
                trans_group_dict[read_id] = group_count

            elif start_pos < group_start_pos: #check if sam sorted
                print("Sam file not sorted!")
                print(read_id)
                sys.exit()


        else: #start new group
            this_scaffold = scaff_name
            group_start_pos = start_pos
            group_end_pos = end_pos
            group_count += 1

            group_trans_list_dict[group_count] = []
            group_trans_list_dict[group_count].append(read_id)
            trans_group_dict[read_id] = group_count

            scaffold_list.append(this_scaffold)

        #check read id add
        if group_trans_list_dict[group_count][-1] != read_id:
            print("cluster not added to group_trans_list_dict")
            print(str(group_count) + " " + read_id)
            sys.exit()

    if bam_flag == "SAM":
        sam_file_obj.close()

    total_group_count = group_count
    ####################################################################################################

    ########################################################################### loop through groups

    merged_obj_dict = {} # merged_obj_dict[final trans id] = merged obj
    gene_count = 0
    trans_check_count = 0 ##########################################################################debugging

    prev_time = track_time(start_time,prev_time)
    print("going through groups: " + str(total_group_count))

    if len(list(group_trans_list_dict.keys())) == 0:
    #  if total_group_count == 0:
        print("Error, no groups found!")
        sys.exit()

    multimap_missing_group_flag = 0

    for i in xrange(total_group_count+1):

        if i not in group_trans_list_dict:
            print("Missing group num, check for multi-maps in SAM file")
            print("This should only occur if you have a multi-map site that no reads are preferring.")
            multimap_missing_group_flag = 1
            continue

        trans_list = group_trans_list_dict[i]

        forward_trans_list = []
        reverse_trans_list = []

        first_trans_id = trans_list[0]

        #separate into forward and reverse
        for trans_id in trans_list:

            if trans_check_count % 1000 == 0:
                print(trans_check_count)
            trans_check_count += 1

            trans_obj = trans_obj_dict[trans_id]
            if trans_obj.strand == "+":
                forward_trans_list.append(trans_id)

            elif trans_obj.strand == "-":
                reverse_trans_list.append(trans_id)

            ################################################# For variation coverage
            scaffold = trans_obj.scaff_name
            if scaffold not in var_coverage_dict:
                print("scaffold not in var_coverage_dict")
                print(scaffold)
                #sys.exit()
                continue

            this_exon_start_list = trans_obj.exon_start_list
            this_exon_end_list = trans_obj.exon_end_list
            for exon_index in xrange(len(this_exon_start_list)):
                this_exon_start = this_exon_start_list[exon_index]
                this_exon_end = this_exon_end_list[exon_index]
                for this_coord in range(this_exon_start,this_exon_end):
                    if this_coord in var_coverage_dict[scaffold]:
                        var_coverage_dict[scaffold][this_coord][trans_id] = 1
                        if trans_id == "":
                            print("Issue with trans id in var cov dict")
                            print(trans_list)
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
                trans_id_list = forward_gene_start_trans_dict[gene_start].keys()
                trans_obj_list = []
                for trans_id in trans_id_list:
                    trans_obj_list.append(trans_obj_dict[trans_id])
                gene_trans_obj_list.append(trans_obj_list)

            if gene_start in reverse_gene_start_trans_dict:
                trans_id_list = reverse_gene_start_trans_dict[gene_start].keys()
                trans_obj_list = []
                for trans_id in trans_id_list:
                    trans_obj_list.append(trans_obj_dict[trans_id])
                gene_trans_obj_list.append(trans_obj_list)

            #loop through list of trans obj lists, usually only one list since unlikely for forward and reverse strand genes to coincide
            for trans_obj_list in gene_trans_obj_list:
                gene_count += 1

                #group transcripts by collapsability
                if fiveprime_cap_flag == "capped":
                    match_trans_group_dict,match_group_trans_dict = simplify_gene_capped(trans_obj_list,fiveprime_cap_flag)
                elif fiveprime_cap_flag == "no_cap":
                    match_trans_group_dict,match_group_trans_dict = simplify_gene_nocap(trans_obj_list,fiveprime_cap_flag)
                else:
                    print("Error with cap flag " + fiveprime_cap_flag)
                    sys.exit()

                merge_obj_list = []
                tmp_count = 0
                for match_group_num in match_group_trans_dict:
                    tmp_count += 1
                    tmp_trans_id = "G" + str(gene_count) + ".tmp." + str(tmp_count)
                    merged_obj = Merged(tmp_trans_id)

                    match_trans_id_list = match_group_trans_dict[match_group_num].keys()
                    match_trans_obj_list = []
                    for match_trans_id in match_trans_id_list:
                        match_trans_obj = trans_obj_dict[match_trans_id]
                        match_trans_obj_list.append(match_trans_obj)

                        merged_obj.add_merged_trans(match_trans_obj)

                    redundant_trans_flag = 0
                    if len(match_trans_obj_list) > 1: #if there are redundant transcripts, collapse
                        redundant_trans_flag = 1
                        collapse_start_list,collapse_end_list,start_wobble_list,end_wobble_list,collapse_sj_start_err_list,collapse_sj_end_err_list,collapse_start_error_nuc_list,collapse_end_error_nuc_list = collapse_transcripts(match_trans_obj_list,fiveprime_cap_flag,collapse_flag)

                        merged_obj.add_merge_info(collapse_start_list,collapse_end_list,start_wobble_list,end_wobble_list,collapse_sj_start_err_list,collapse_sj_end_err_list,collapse_start_error_nuc_list,collapse_end_error_nuc_list )
                    else: # if only one transcript
                        exon_start_list = match_trans_obj_list[0].exon_start_list
                        exon_end_list = match_trans_obj_list[0].exon_end_list
                        start_wobble_list = [0] * len(exon_start_list)
                        end_wobble_list = [0] * len(exon_start_list)

                        exon_start_list.sort()
                        exon_end_list.sort()

                        collapse_sj_start_err_list = []
                        collapse_sj_end_err_list = []
                        solo_trans_obj = match_trans_obj_list[0]
                        max_exon_num = len(exon_start_list)

                        collapse_start_error_nuc_list = []
                        collapse_end_error_nuc_list = []

                        for exon_index in xrange(len(exon_start_list)):  # go from 3 prime end
                            e_start_priority, e_end_priority, e_start_priority_error,e_end_priority_error = sj_error_priority_finder(solo_trans_obj, exon_index, max_exon_num)  ####################################

                            collapse_sj_start_err_list.append(e_start_priority)
                            collapse_sj_end_err_list.append(e_end_priority)

                            collapse_start_error_nuc_list.append(e_start_priority_error)
                            collapse_end_error_nuc_list.append(e_end_priority_error)

                        #collapse_sj_start_err_list = trans_obj.sj_pre_error_list
                        #collapse_sj_end_err_list = trans_obj.sj_post_error_list

                        merged_obj.add_merge_info(exon_start_list,exon_end_list,start_wobble_list,end_wobble_list,collapse_sj_start_err_list,collapse_sj_end_err_list,collapse_start_error_nuc_list,collapse_end_error_nuc_list )

                    merge_obj_list.append(merged_obj)


                sorted_merge_obj_list = sort_transcripts(merge_obj_list)

                trans_count = 0
                for merged_obj in sorted_merge_obj_list:
                    trans_count += 1
                    final_trans_id = "G" + str(gene_count) + "." + str(trans_count)
                    merged_obj.trans_id = final_trans_id
                    print(final_trans_id)

                    merged_obj_dict[final_trans_id] = merged_obj

                    #write out to bed file
                    bed_line = merged_obj.format_bed_line()
                    outfile_bed.write(bed_line)
                    outfile_bed.write("\n")

                    #write out to transcript report file
                    trans_report_line = merged_obj.format_trans_report_line()
                    outfile_trans_report.write(trans_report_line)
                    outfile_trans_report.write("\n")

                    #write out to rt switch file
                    num_junct_bind,bind_seq_dict = detect_rt_switch(merged_obj)
                    if num_junct_bind > 0 :
                        for junct_num in bind_seq_dict:
                            first_seq = "".join(bind_seq_dict[junct_num][1])
                            second_seq = "".join(bind_seq_dict[junct_num][2])
                            rec_comp_seq = "".join(bind_seq_dict[junct_num][3])

                            rtswitch_line = "\t".join([final_trans_id,str(junct_num),first_seq,second_seq,rec_comp_seq])

                            #outfile_rtswitch.write(rtswitch_line)
                            #outfile_rtswitch.write("\n")

                    #write out to transcript cluster report file
                    for merged_trans_id in merged_obj.merged_trans_dict:
                        merged_trans_obj = merged_obj.merged_trans_dict[merged_trans_id]
                        cluster_id = merged_trans_obj.cluster_id
                        scaff_name = merged_trans_obj.scaff_name
                        strand = merged_trans_obj.strand
                        start_pos = merged_trans_obj.start_pos
                        end_pos = merged_trans_obj.end_pos

                        exon_start_line,exon_end_line = merged_trans_obj.make_exon_start_end_lines()

                        merge_trans_bed_line = merged_trans_obj.format_bed_line(final_trans_id)
                        #trans_clust_list = []
                        #trans_clust_list.append(final_trans_id)
                        #trans_clust_list.append(cluster_id)
                        #trans_clust_list.append(scaff_name)
                        #trans_clust_list.append(strand)
                        #trans_clust_list.append(str(start_pos))
                        #trans_clust_list.append(str(end_pos))
                        #trans_clust_list.append(exon_start_line)
                        #trans_clust_list.append(exon_end_line)

                        #trans_clust_line = "\t".join(trans_clust_list)
                        outfile_trans_clust_report.write(merge_trans_bed_line)
                        outfile_trans_clust_report.write("\n")

                        #write out to polya file
                        downstream_seq = "".join(merged_trans_obj.downstream_seq)
                        dseq_length = merged_trans_obj.dseq_length
                        a_count = merged_trans_obj.a_count
                        a_percent = merged_trans_obj.a_percent * 100


                        if a_percent > a_perc_thresh:
                            a_percent_string =str(round(a_percent,2))
                            polya_file_line = "\t".join([cluster_id,final_trans_id,strand,a_percent_string,str(a_count),downstream_seq])
                            outfile_polya.write(polya_file_line)
                            outfile_polya.write("\n")


    # original mode end ##############################################
    ############################################################################################################
    ############################################################################################################
    ############################################################################################################

############################################################################################################
############################################################################################################
# no multimap mode start (low mem) ##############################################

##################################################################################################
##################################################################################################
# variation functions for low mem mode
def process_variation(scaffold,variation_dict,var_coverage_dict,var_support_threshold):

    cov_group_var_dict = {} # cov_group_var_dict[cov group line][position line] = 1
    cov_group_var_list = []

    var_type_list = []
    var_type_list.append("H")
    var_type_list.append("S")
    var_type_list.append("M")
    var_type_list.append("I")
    var_type_list.append("D")

    if scaffold not in variation_dict:
        print("error with scaffold match in loci_variation")
        print(variation_dict.keys())
        sys.exit()


    position_list = []
    position_list = list(variation_dict[scaffold].keys())
    position_list.sort()

    for var_pos in position_list:

        var_cov_trans_id_list = list(var_coverage_dict[scaffold][var_pos].keys())
        var_cov_trans_id_list.sort()
        var_coverage = len(var_cov_trans_id_list)

        var_pos_accept_flag = 0 # Use this to signal if a variation ahs passed threshold for this position
        for var_type in var_type_list:

            if var_type not in variation_dict[scaffold][var_pos]:
                continue

            for alt_seq in variation_dict[scaffold][var_pos][var_type]:

                read_list = list(variation_dict[scaffold][var_pos][var_type][alt_seq].keys())

                var_support_count = len(read_list)

                if var_support_count >= var_support_threshold:
                    var_pos_accept_flag = 1

                    #scaffold        position        type    alt_allele      count  cov_count   cluster_list
                    var_outlist = []
                    var_outlist.append(scaffold)
                    var_outlist.append(str(var_pos))
                    var_outlist.append(var_type)
                    var_outlist.append(alt_seq)
                    var_outlist.append(str(var_support_count))
                    var_outlist.append(str(var_coverage))

                    read_line = ",".join(read_list)
                    var_outlist.append(read_line)

                    var_outline = "\t".join(var_outlist)
                    outfile_variant.write(var_outline)
                    outfile_variant.write("\n")

        #Update variant coverage

        if var_pos_accept_flag == 1:
            var_cov_trans_line = ",".join(var_cov_trans_id_list)
            position_line = "_".join([scaffold,str(var_pos)])
            if var_cov_trans_line not in cov_group_var_dict:
                cov_group_var_dict[var_cov_trans_line] = {}
                cov_group_var_list.append(var_cov_trans_line)
            cov_group_var_dict[var_cov_trans_line][position_line] = 1



    ################################################################################# write to var coverage file


    for cov_line in cov_group_var_list:
        position_list = list(cov_group_var_dict[cov_line].keys())
        position_list.sort()

        all_pos_line = ",".join(position_list)

        varcov_file_line = "\t".join([all_pos_line,cov_line])

        outfile_varcov.write(varcov_file_line)
        outfile_varcov.write("\n")

    del variation_dict
    del var_coverage_dict



##################################################################################################
##################################################################################################


def process_loci(this_trans_obj_dict,trans_list,this_gene_count):

    merged_obj_dict = {}  # merged_obj_dict[final trans id] = merged obj


    forward_trans_list = []
    reverse_trans_list = []

    first_trans_id = trans_list[0]

    #separate into forward and reverse
    for trans_id in trans_list:

        #if trans_check_count % 1000 == 0:
        #    print(trans_check_count)
        #trans_check_count += 1

        trans_obj = this_trans_obj_dict[trans_id]
        if trans_obj.strand == "+":
            forward_trans_list.append(trans_id)

        elif trans_obj.strand == "-":
            reverse_trans_list.append(trans_id)


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
            trans_id_list = forward_gene_start_trans_dict[gene_start].keys()
            trans_obj_list = []
            for trans_id in trans_id_list:
                trans_obj_list.append(this_trans_obj_dict[trans_id])
            gene_trans_obj_list.append(trans_obj_list)

        if gene_start in reverse_gene_start_trans_dict:
            trans_id_list = reverse_gene_start_trans_dict[gene_start].keys()
            trans_obj_list = []
            for trans_id in trans_id_list:
                trans_obj_list.append(this_trans_obj_dict[trans_id])
            gene_trans_obj_list.append(trans_obj_list)

        #loop through list of trans obj lists, usually only one list since unlikely for forward and reverse strand genes to coincide
        for trans_obj_list in gene_trans_obj_list:
            this_gene_count += 1

            #group transcripts by collapsability
            if fiveprime_cap_flag == "capped":
                match_trans_group_dict,match_group_trans_dict = simplify_gene_capped(trans_obj_list,fiveprime_cap_flag)
            elif fiveprime_cap_flag == "no_cap":
                match_trans_group_dict,match_group_trans_dict = simplify_gene_nocap(trans_obj_list,fiveprime_cap_flag)
            else:
                print("Error with cap flag " + fiveprime_cap_flag)
                sys.exit()

            merge_obj_list = []
            tmp_count = 0
            for match_group_num in match_group_trans_dict:
                tmp_count += 1
                tmp_trans_id = "G" + str(this_gene_count) + ".tmp." + str(tmp_count)
                merged_obj = Merged(tmp_trans_id)

                match_trans_id_list = match_group_trans_dict[match_group_num].keys()
                match_trans_obj_list = []
                for match_trans_id in match_trans_id_list:
                    match_trans_obj = this_trans_obj_dict[match_trans_id]
                    match_trans_obj_list.append(match_trans_obj)

                    merged_obj.add_merged_trans(match_trans_obj)

                redundant_trans_flag = 0
                if len(match_trans_obj_list) > 1: #if there are redundant transcripts, collapse
                    redundant_trans_flag = 1
                    collapse_start_list,collapse_end_list,start_wobble_list,end_wobble_list,collapse_sj_start_err_list,collapse_sj_end_err_list,collapse_start_error_nuc_list,collapse_end_error_nuc_list = collapse_transcripts(match_trans_obj_list,fiveprime_cap_flag,collapse_flag)

                    merged_obj.add_merge_info(collapse_start_list,collapse_end_list,start_wobble_list,end_wobble_list,collapse_sj_start_err_list,collapse_sj_end_err_list,collapse_start_error_nuc_list,collapse_end_error_nuc_list )
                else: # if only one transcript
                    exon_start_list = match_trans_obj_list[0].exon_start_list
                    exon_end_list = match_trans_obj_list[0].exon_end_list
                    start_wobble_list = [0] * len(exon_start_list)
                    end_wobble_list = [0] * len(exon_start_list)

                    exon_start_list.sort()
                    exon_end_list.sort()

                    collapse_sj_start_err_list = []
                    collapse_sj_end_err_list = []
                    solo_trans_obj = match_trans_obj_list[0]
                    max_exon_num = len(exon_start_list)

                    collapse_start_error_nuc_list = []
                    collapse_end_error_nuc_list = []

                    for exon_index in xrange(len(exon_start_list)):  # go from 3 prime end
                        e_start_priority, e_end_priority, e_start_priority_error,e_end_priority_error = sj_error_priority_finder(solo_trans_obj, exon_index, max_exon_num)  ####################################

                        collapse_sj_start_err_list.append(e_start_priority)
                        collapse_sj_end_err_list.append(e_end_priority)

                        collapse_start_error_nuc_list.append(e_start_priority_error)
                        collapse_end_error_nuc_list.append(e_end_priority_error)

                    #collapse_sj_start_err_list = trans_obj.sj_pre_error_list
                    #collapse_sj_end_err_list = trans_obj.sj_post_error_list

                    merged_obj.add_merge_info(exon_start_list,exon_end_list,start_wobble_list,end_wobble_list,collapse_sj_start_err_list,collapse_sj_end_err_list,collapse_start_error_nuc_list,collapse_end_error_nuc_list )

                merge_obj_list.append(merged_obj)


            sorted_merge_obj_list = sort_transcripts(merge_obj_list)

            trans_count = 0
            for merged_obj in sorted_merge_obj_list:
                trans_count += 1
                final_trans_id = "G" + str(this_gene_count) + "." + str(trans_count)
                merged_obj.trans_id = final_trans_id
                print(final_trans_id)

                merged_obj_dict[final_trans_id] = merged_obj

                #write out to bed file
                bed_line = merged_obj.format_bed_line()
                outfile_bed.write(bed_line)
                outfile_bed.write("\n")

                #write out to transcript report file
                trans_report_line = merged_obj.format_trans_report_line()
                outfile_trans_report.write(trans_report_line)
                outfile_trans_report.write("\n")

                #write out to rt switch file
                num_junct_bind,bind_seq_dict = detect_rt_switch(merged_obj)
                if num_junct_bind > 0 :
                    for junct_num in bind_seq_dict:
                        first_seq = "".join(bind_seq_dict[junct_num][1])
                        second_seq = "".join(bind_seq_dict[junct_num][2])
                        rec_comp_seq = "".join(bind_seq_dict[junct_num][3])

                        rtswitch_line = "\t".join([final_trans_id,str(junct_num),first_seq,second_seq,rec_comp_seq])

                        #outfile_rtswitch.write(rtswitch_line)
                        #outfile_rtswitch.write("\n")

                #write out to transcript cluster report file
                for merged_trans_id in merged_obj.merged_trans_dict:
                    merged_trans_obj = merged_obj.merged_trans_dict[merged_trans_id]
                    cluster_id = merged_trans_obj.cluster_id
                    scaff_name = merged_trans_obj.scaff_name
                    strand = merged_trans_obj.strand
                    start_pos = merged_trans_obj.start_pos
                    end_pos = merged_trans_obj.end_pos

                    exon_start_line,exon_end_line = merged_trans_obj.make_exon_start_end_lines()

                    merge_trans_bed_line = merged_trans_obj.format_bed_line(final_trans_id)

                    #trans_clust_line = "\t".join(trans_clust_list)
                    outfile_trans_clust_report.write(merge_trans_bed_line)
                    outfile_trans_clust_report.write("\n")

                    #write out to polya file
                    downstream_seq = "".join(merged_trans_obj.downstream_seq)
                    dseq_length = merged_trans_obj.dseq_length
                    a_count = merged_trans_obj.a_count
                    a_percent = merged_trans_obj.a_percent * 100


                    if a_percent > a_perc_thresh:
                        a_percent_string =str(round(a_percent,2))
                        polya_file_line = "\t".join([cluster_id,final_trans_id,strand,a_percent_string,str(a_count),downstream_seq])
                        outfile_polya.write(polya_file_line)
                        outfile_polya.write("\n")


    #del this_trans_obj_dict
    #del trans_list

    return this_gene_count


############################################################################################################
############################################################################################################

if run_mode_flag == "low_mem":
    gene_count = 0

    trans_obj_dict = {} # trans_obj_dict[cluster id] = trans obj
    group_trans_list = []


    for line in sam_file_contents:

        #if sam_count == 0:
        #    print(line)

        line_split = line.split("\t")

        if line.startswith("@"):
            continue

        if line == "":
            continue

        sam_count += 1
        if sam_count % 5000 == 0:
            print("sam count " + str(sam_count))

        read_id = line_split[0]
        sam_flag = int(line_split[1])
        scaff_name = line_split[2]
        start_pos = int(line_split[3])

        cigar = line_split[5]
        read_seq = line_split[9]
        seq_list = list(read_seq)
        mapped_flag = sam_flag_dict[sam_flag]

        ####################################
        #Check sam and gmap strand info!!!
        #get strand information from gmap flag

        xs_flag = "na"
        for field in line_split:
            if "XS:A:" in field:
                xs_flag = field.split(":")[-1]

        if mapped_flag == "forward_strand"  and xs_flag == "-":
            outline_strand = "\t".join([read_id,scaff_name,str(start_pos),cigar,"+-"])
            outfile_strand.write(outline_strand)
            outfile_strand.write("\n")
        elif mapped_flag == "reverse_strand" and xs_flag == "+":
            outline_strand = "\t".join([read_id,scaff_name,str(start_pos),cigar,"-+"])
            outfile_strand.write(outline_strand)
            outfile_strand.write("\n")

        #
        # Above: Check sam and gmap strand info!!!
        ####################################

        if mapped_flag == "unmapped" or mapped_flag == "not_primary" or mapped_flag == "chimeric" :
            unmapped_dict[read_id] = 1
            accept_flag = mapped_flag # added this 2019/03/04
            percent_coverage = "NA"
            percent_identity = "NA"
            error_line = "NA"
            quality_percent = "NA"
            length = "NA"
            strand = "NA"
            cigar = "NA"

            cluster_line = "\t".join([read_id,mapped_flag,accept_flag,percent_coverage,percent_identity,error_line, length, cigar])
            outfile_cluster.write(cluster_line)
            outfile_cluster.write("\n")
            continue


        map_seq_length = mapped_seq_length(cigar)

        [end_pos,exon_start_list,exon_end_list] = trans_coordinates(start_pos,cigar)

        [h_count,s_count,i_count,d_count,mis_count,nomatch_dict,sj_pre_error_list,sj_post_error_list] =  calc_error_rate_lowmem(start_pos,cigar,seq_list,scaff_name,read_id)

        trans_obj = Transcript(read_id)
        trans_obj.add_sam_info(sam_flag,scaff_name,start_pos,cigar,read_seq,seq_list)
        trans_obj.add_map_seq_length(map_seq_length)
        trans_obj.add_exon_coords(end_pos,exon_start_list,exon_end_list)
        trans_obj.add_mismatch(h_count,s_count,i_count,d_count,mis_count,nomatch_dict,sj_pre_error_list,sj_post_error_list)

        ##### 2020/07/27 sj hash
        #trans_obj.make_sj_hash_string()
        ##### 2020/07/27 sj hash


        percent_coverage = trans_obj.calc_coverage()
        percent_identity = trans_obj.calc_identity()

        percent_coverage_str = str(round(percent_coverage,2))
        percent_identity_str = str(round(percent_identity,2))

        error_line = trans_obj.make_error_line()
        seq_length = trans_obj.seq_length
        strand = trans_obj.strand

        multimap_flag = 0

        if percent_coverage < coverage_threshold or percent_identity < identity_threshold:
            accept_flag = "discarded"
            cluster_line = "\t".join([read_id,mapped_flag,accept_flag,percent_coverage_str,percent_identity_str,error_line, str(seq_length), cigar])
            outfile_cluster.write(cluster_line)
            outfile_cluster.write("\n")
            #skip the transcript because the mapping is poor
            continue

        else:

            bad_sj_flag,bad_sj_num_list,bad_sj_num_pre_list,bad_sj_num_post_list,bad_sj_error_count_list,lde_outline = sj_error_local_density(trans_obj)

            outfile_lde.write(lde_outline)
            outfile_lde.write("\n")

            if bad_sj_flag > 0:


                accept_flag = "local_density_error"
                cluster_line = "\t".join([read_id,mapped_flag,accept_flag,percent_coverage_str,percent_identity_str,error_line, str(seq_length), cigar])
                outfile_cluster.write(cluster_line)
                outfile_cluster.write("\n")

                continue

            accept_flag = "accepted"
            cluster_line = "\t".join([read_id,mapped_flag,accept_flag,percent_coverage_str,percent_identity_str,error_line, str(seq_length), cigar])
            outfile_cluster.write(cluster_line)
            outfile_cluster.write("\n")

            #only run poly detection on accepted transcripts
            downstream_seq,dseq_length,a_count,n_count,a_percent,n_percent = detect_polya(trans_obj,a_window)
            trans_obj.add_polya_info(downstream_seq,dseq_length,a_count,n_count,a_percent,n_percent)

            trans_obj_dict[read_id] = trans_obj

        #group trans by start and end coords
        if this_scaffold == "none":
            this_scaffold = scaff_name
            group_start_pos = start_pos
            group_end_pos = end_pos

            group_trans_list = []
            group_trans_list.append(read_id)

            scaffold_list.append(this_scaffold)

            continue

        if scaff_name == this_scaffold:

            if start_pos >= group_start_pos and start_pos <= group_end_pos: #add to group

                # this is so we can get variation information after confirming the read belongs to this group
                # 2020/07/31
                calc_variation(start_pos, cigar, seq_list, scaff_name, read_id)

                group_trans_list.append(read_id)
                #update group end position
                if end_pos > group_end_pos:
                    group_end_pos = end_pos

            elif start_pos > group_end_pos: #start new group #################################

                # clear new read id from new group
                trans_obj_dict.pop(read_id, None)

                gene_count = process_loci(trans_obj_dict,group_trans_list, gene_count) #################################### 2020/07/30
                #if variation_dict:
                #    if len(group_trans_list) >= var_support_threshold:
                #        process_variation(this_scaffold, variation_dict, var_coverage_dict,var_support_threshold)

                #        #refresh variation dict
                #        #del variation_dict
                #        #del var_coverage_dict
                #        variation_dict = {}
                #        var_coverage_dict = {}

                # add variation information from first read of new group
                # 2020/07/31
                #calc_variation(start_pos, cigar, seq_list, scaff_name, read_id)
                [h_count,s_count,i_count,d_count,mis_count,nomatch_dict,sj_pre_error_list,sj_post_error_list] =  calc_error_rate_lowmem(start_pos,cigar,seq_list,scaff_name,read_id)

                group_count += 1

                group_start_pos = start_pos
                group_end_pos = end_pos

                #del trans_obj_dict
                trans_obj_dict = {}  # trans_obj_dict[cluster id] = trans obj # refresh this dict to save memory
                trans_obj_dict[read_id] = trans_obj

                #del group_trans_list
                group_trans_list = []  # refresh this dict to save memory
                group_trans_list.append(read_id)



            elif start_pos < group_start_pos: #check if sam sorted
                print("Sam file not sorted!")
                print(read_id)
                sys.exit()


        else: #start new group #################################

            # clear new read id from new group
            trans_obj_dict.pop(read_id, None)

            gene_count = process_loci(trans_obj_dict, group_trans_list, gene_count)  #################################### 2020/07/30

            #if variation_dict:
            #    if len(group_trans_list) >= var_support_threshold:
            #        process_variation(this_scaffold, variation_dict, var_coverage_dict,var_support_threshold)
            #
            #        #refresh variation dict
            #        #del variation_dict
            #        #del var_coverage_dict
            #        variation_dict = {}
            #        var_coverage_dict = {}

            # add variation information from first read of new group
            # 2020/07/31
            #calc_variation(start_pos, cigar, seq_list, scaff_name, read_id)
            [h_count,s_count,i_count,d_count,mis_count,nomatch_dict,sj_pre_error_list,sj_post_error_list] =  calc_error_rate_lowmem(start_pos,cigar,seq_list,scaff_name,read_id)

            this_scaffold = scaff_name
            group_start_pos = start_pos
            group_end_pos = end_pos
            group_count += 1

            scaffold_list.append(this_scaffold)

            #del trans_obj_dict
            trans_obj_dict = {}  # trans_obj_dict[cluster id] = trans obj # refresh this dict to save memory
            trans_obj_dict[read_id] = trans_obj

            #del group_trans_list
            group_trans_list = []  # refresh this dict to save memory
            group_trans_list.append(read_id)

    # make sure to handle the last group for the low mem mode
    # this is not an issue with the original mode because you collect groups as you go along.
    # but in low mem mode you only process a group when a new group is found.
    gene_count = process_loci(trans_obj_dict, group_trans_list, gene_count)  #################################### 2020/07/30
    #if variation_dict:
    #    if len(group_trans_list) >= var_support_threshold:
    #        process_variation(this_scaffold, variation_dict, var_coverage_dict,var_support_threshold)

    #        #refresh variation dict
    #        #del variation_dict
    #        #del var_coverage_dict
    #        variation_dict = {}
    #        var_coverage_dict = {}

    if bam_flag == "SAM":
        sam_file_obj.close()

    total_group_count = group_count
    ####################################################################################################

    ########################################################################### loop through groups

    prev_time = track_time(start_time,prev_time)




# no multimap mode end (low mem) ##############################################
############################################################################################################
############################################################################################################

################################################################################# write to variation file


if run_mode_flag == "original":

    cov_group_var_dict = {} # cov_group_var_dict[cov group line][position line] = 1
    cov_group_var_list = []

    var_type_list = []
    var_type_list.append("H")
    var_type_list.append("S")
    var_type_list.append("M")
    var_type_list.append("I")
    var_type_list.append("D")

    var_support_threshold = 5

    prev_time = track_time(start_time,prev_time)
    print("Writing variant file")
    for scaffold in scaffold_list:

        if scaffold not in variation_dict:
            continue

        position_list = []
        position_list = list(variation_dict[scaffold].keys())
        position_list.sort()

        for var_pos in position_list:

            var_cov_trans_id_list = list(var_coverage_dict[scaffold][var_pos].keys())
            var_cov_trans_id_list.sort()
            var_coverage = len(var_cov_trans_id_list)


            ########################################################################################2020/12/14
            if var_pos > len(fasta_dict[scaffold]) or var_pos < 0:
                print("Read mapping off scaffold")
                print(scaffold +" : "+ str(var_pos))
                print(var_cov_trans_id_list)
                continue
            #    print(variation_dict[scaffold][var_pos])

            ########################################################################################


            ref_allele = fasta_dict[scaffold][var_pos]

            var_pos_accept_flag = 0 # Use this to signal if a variation ahs passed threshold for this position
            for var_type in var_type_list:

                if var_type not in variation_dict[scaffold][var_pos]:
                    continue

                if var_type != "M":
                    ref_allele = "NA"


                for alt_seq in variation_dict[scaffold][var_pos][var_type]:

                    read_list = list(variation_dict[scaffold][var_pos][var_type][alt_seq].keys())

                    var_support_count = len(read_list)

                    if var_support_count >= var_support_threshold:
                        var_pos_accept_flag = 1

                        #scaffold        position        type    alt_allele      count  cov_count   cluster_list
                        var_outlist = []
                        var_outlist.append(scaffold)
                        var_outlist.append(str(var_pos))
                        var_outlist.append(var_type)
                        var_outlist.append(ref_allele)
                        var_outlist.append(alt_seq)
                        var_outlist.append(str(var_support_count))
                        var_outlist.append(str(var_coverage))

                        read_line = ",".join(read_list)
                        var_outlist.append(read_line)

                        var_outline = "\t".join(var_outlist)
                        outfile_variant.write(var_outline)
                        outfile_variant.write("\n")

            #Update variant coverage

            if var_pos_accept_flag == 1:
                var_cov_trans_line = ",".join(var_cov_trans_id_list)
                position_line = "_".join([scaffold,str(var_pos)])
                if var_cov_trans_line not in cov_group_var_dict:
                    cov_group_var_dict[var_cov_trans_line] = {}
                    cov_group_var_list.append(var_cov_trans_line)
                cov_group_var_dict[var_cov_trans_line][position_line] = 1



    prev_time = track_time(start_time,prev_time)

    ################################################################################# write to var coverage file


    for cov_line in cov_group_var_list:
        position_list = list(cov_group_var_dict[cov_line].keys())
        position_list.sort()

        all_pos_line = ",".join(position_list)

        varcov_file_line = "\t".join([all_pos_line,cov_line])

        outfile_varcov.write(varcov_file_line)
        outfile_varcov.write("\n")





prev_time = track_time(start_time,prev_time)

if run_mode_flag == "original":
    if multimap_missing_group_flag == 1:
        print("Missing group num, check for multi-maps in SAM file")
        print("This should only occur if you have a multi-map site that no reads are preferring.")

print("TAMA Collapse has successfully finished running!")