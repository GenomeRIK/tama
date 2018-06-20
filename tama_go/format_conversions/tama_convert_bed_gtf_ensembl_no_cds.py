import re
import sys
import time
#from Bio import SeqIO


#
# This script is used to convert the pacbio bed format file into a gtf file that mimics Ensembl's format
# for use on bed files wihtout CDS information
# For use with bed files that are an output of tama merge or tama collapse
#
# Note bed is 0 based and gtf is 1 based
#


print("opening pacbio bed file")
bed_file = sys.argv[1]
bed_file_contents = open(bed_file).read().rstrip("\n").split("\n")

outfile_name = sys.argv[2]
outfile = open(outfile_name,"w")


#outline = "#!genome-build galGal4\n"
#outfile.write(outline)
#outline = "#!genome-version galGal4\n"
#outfile.write(outline)
#outline = "#!genome-date 2011-11\n"
#outfile.write(outline)
#outline = "#!genome-build-accession NCBI:GCA_000002315.2\n"
#outfile.write(outline)
#outline = "#!genebuild-last-updated 2015-03\n"
#outfile.write(outline)
#outline = "#!genebuild-version 1.0\n"
#outfile.write(outline)



def calc_end(start,block_size):
    end = int(start) + int(block_size)
    return str(end)

def calc_exon_start(t_start,e_start):
    coordinate_start = int(e_start) + int(t_start)
    return str(coordinate_start)

def convert_str_list_to_int(str_list):
    int_list = []
    for string in str_list:
        int_list.append(int(string))
    return int_list

def convert_int_list_to_str(int_list):
    str_list = []
    for integer in int_list:
        str_list.append(str(integer))
    return str_list

class Transcript:
    def __init__(self,trans_id,gene_id, uniq_trans_id, chrom,strand,t_start,t_end,starts,blocks,num_exons):
        self.trans_id = trans_id
        self.gene_id = gene_id
        self.uniq_trans_id = uniq_trans_id
        self.chrom = chrom
        self.t_start = str(t_start)
        self.t_end = str(t_end)
        self.strand = strand
        t_start_list = [int(t_start)] * int(num_exons)
        start_list = starts.split(",")
        start_list = filter(None, start_list)
        block_list = blocks.split(",")
        block_list = filter(None, block_list)
        
        self.bed_starts = starts
        self.bed_blocks = blocks
        
        #coordinate starts and ends
        self.start_list = map(calc_exon_start,t_start_list,start_list)
        self.end_list =  map(calc_end,self.start_list,block_list)
        self.num_exons = num_exons



gene_list = [] # used to keep order of genes
gene_trans_dict = {} #gene_trans_dict[gene_id] = trans list
trans_dict = {} # trans_dict[trans_id] = trans_obj

print("Going through bed file")
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
    
    id_split = id_line.split(";")
    gene_id = id_split[0]
    trans_id = id_split[1]
#    uniq_trans_id = id_split[2]
    uniq_trans_id = id_split[1]
    

    if gene_id not in gene_trans_dict:
        gene_trans_dict[gene_id] = []
        gene_list.append(gene_id)
    
    gene_trans_dict[gene_id].append(trans_id)
    
    trans_dict[trans_id] = Transcript(trans_id,gene_id, uniq_trans_id, chrom,strand,t_start,t_end,block_starts,block_sizes,num_exons)
    
source = "PBRI"

def format_gene_line(chrom,gene_start,gene_end,strand,gene_id,source):
    #anno_line = "gene_id \"" + gene_id + "\"; gene_source \"" + source + "\";"
    anno_line = "gene_id \"" + gene_id + "\";"
    # convert 0 coord bed to 1 coord gtf
    gene_start = gene_start + 1
    
    outline = "\t".join([chrom,source,"gene",str(gene_start),str(gene_end),".",strand,".",anno_line])
    return outline

def format_trans_line(chrom,t_start,t_end,strand,gene_id,trans_id,uniq_trans_id,source):
    # convert 0 coord bed to 1 coord gtf
    t_start = int(t_start) + 1
    
    anno_line = "gene_id \"" + gene_id + "\"; transcript_id \"" + trans_id + "\"; uniq_trans_id \"" + uniq_trans_id + "\";"
    outline = "\t".join([chrom,source,"transcript",str(t_start),str(t_end),".",strand,".",anno_line])
    return outline

def format_exon_line(chrom,e_start,e_end,strand,gene_id,trans_id,uniq_trans_id,source,e_num):
    # convert 0 coord bed to 1 coord gtf
    e_start = int(e_start) + 1
    
    anno_line = "gene_id \"" + gene_id + "\"; transcript_id \"" + trans_id + "\"; exon_number \"" + str(e_num) +"\"; uniq_trans_id \"" + uniq_trans_id + "\";"
    outline = "\t".join([chrom,source,"exon",str(e_start),str(e_end),".",strand,".",anno_line])
    return outline

print("Writing out gtf file")
for gene_id in gene_list:
    min_start = -1
    max_end = -1
    chrom = -1
    strand = ""
    for trans_id in gene_trans_dict[gene_id]:
        #beginning condition
        t_start = int(trans_dict[trans_id].t_start)
        t_end = int(trans_dict[trans_id].t_end)
        if min_start == -1:
            min_start = t_start
            max_end = t_end
            chrom = trans_dict[trans_id].chrom
            strand = trans_dict[trans_id].strand
            continue
        
        if t_start < min_start:
            min_start = t_start
        if t_end > max_end:
            max_end = t_end
            
        
    outline = format_gene_line(chrom,min_start,max_end,strand,gene_id,source)
    outfile.write(outline)
    outfile.write("\n")
    
    for trans_id in gene_trans_dict[gene_id]:
        t_start = trans_dict[trans_id].t_start
        t_end = trans_dict[trans_id].t_end
        uniq_trans_id = trans_dict[trans_id].uniq_trans_id
        
        #write transcript line
        outline = format_trans_line(chrom,t_start,t_end,strand,gene_id,trans_id,uniq_trans_id,source)
        outfile.write(outline)
        outfile.write("\n")
        
        #write exon lines
        for i,start in enumerate(trans_dict[trans_id].start_list):
            e_start = start
            e_end = trans_dict[trans_id].end_list[i]
            num_exons = int(trans_dict[trans_id].num_exons)
            if strand == "+":
                e_num = i + 1
            elif strand == "-":
                e_num = num_exons - i
            else:
                print("error with strand " + trans_id)
                sys.exit()
            
            outline = format_exon_line(chrom,e_start,e_end,strand,gene_id,trans_id,uniq_trans_id,source,e_num)
            outfile.write(outline)
            outfile.write("\n")
            

            
            
        
        

