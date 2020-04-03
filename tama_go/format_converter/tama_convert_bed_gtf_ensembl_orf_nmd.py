
import re
import sys
import time
#from Bio import SeqIO


#
# This script is used to convert the pacbio bed format file into a gtf file that mimics Ensembl's format
# For use with bed files that are an output of tama orf/nmd prediction pipeline
# includes CDS information
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

gene_source = "tama"
trans_source = "tama"

def calc_end(start,block_size):
    end = int(start) + int(block_size) - 1 # adjust for bed 0 base and gtf 1 base coords
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
    def __init__(self,bed_line):

        line_split = bed_line.split("\t")
        chrom = line_split[0]
        t_start = int(line_split[1]) + 1
        t_end = int(line_split[2])
        id_line = line_split[3]
        strand = line_split[5]
        num_exons = int(line_split[9])
        blocks = line_split[10]
        starts = line_split[11]

        cds_start = int(line_split[6]) + 1
        cds_end = int(line_split[7])

        id_split = id_line.split(";")
        gene_id = id_split[0]
        trans_id = id_split[1]
        prot_id = id_split[2]
        degrade_flag = id_split[3]
        match_flag = id_split[4]
        nmd_flag = id_split[5]

        self.prot_id = prot_id
        self.degrade_flag = degrade_flag
        self.match_flag = match_flag
        self.nmd_flag = nmd_flag


        self.trans_id = trans_id
        self.gene_id = gene_id
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


        self.cds_start = cds_start#################################################
        self.cds_end = cds_end #################################################
        

    def write_trans(self):
        
        if int(self.num_exons) != len(self.start_list):
            print("Error with number of exons")
            sys.exit()


        outline = format_trans_line(self.chrom,self.t_start,self.t_end,self.strand,self.gene_id,self.trans_id,source)
        outfile.write(outline)
        outfile.write("\n")

        five_utr_start_list = []
        five_utr_end_list = []
        five_utr_num_list = []

        three_utr_start_list = []
        three_utr_end_list = []
        three_utr_num_list = []

        
        for i in xrange(int(self.num_exons)):
            if self.strand == "+":
                e_index = i
                e_num = e_index + 1
                
                e_start = int(self.start_list[e_index])
                e_end = int(self.end_list[e_index])
                
                region_type = "exon"
                frame = "."

                outline = format_exon_line(region_type, frame, self.chrom, str(e_start), str(e_end), self.strand,self.gene_id, self.trans_id, source, e_num, self.prot_id, self.degrade_flag,self.match_flag, self.nmd_flag)
                outfile.write(outline)
                outfile.write("\n")

                ####################################start here!!!!!

                if self.cds_start == 1 and self.cds_end == 0: # no CDS region
                    continue

                if e_start > int(self.cds_end):  # exon is completely 3' utr for pos strand
                    three_utr_start = e_start
                    three_utr_start_list.append(three_utr_start)
                    three_utr_end = e_end
                    three_utr_end_list.append(three_utr_end)
                    three_utr_num_list.append(e_num)
                    continue

                if int(self.cds_start) > e_end:    # exon is completely 5' utr on pos strand
                    five_utr_start_list.append(e_start)
                    five_utr_end_list.append(e_end)
                    five_utr_num_list.append(e_num)
                    continue


                if e_start > int(self.cds_start) and e_end < int(self.cds_end):  # exon is within coding region but does not contain start or end
                    region_type = "CDS"
                    frame = "."
                    outline = format_exon_line(region_type, frame, self.chrom, str(e_start), str(e_end), self.strand,self.gene_id, self.trans_id, source, e_num, self.prot_id, self.degrade_flag, self.match_flag, self.nmd_flag)
                    outfile.write(outline)
                    outfile.write("\n")
                    continue

                ## all following conditions have cds start/stop regions in exon

                stop_codon_flag = 0
                start_codon_flag = 0

                cds_e_start = e_start
                cds_e_end = e_end

                if e_start < int(self.cds_end) and e_end >= int(self.cds_end):  # cds end is within this exon

                    stop_codon_start = int(self.cds_end) - 2
                    stop_codon_end = int(self.cds_end)

                    cds_e_end = stop_codon_start - 1  # pos cds stop codon, cds does not includes it

                    three_utr_start = stop_codon_end + 1
                    three_utr_end = e_end
                    three_utr_start_list.append(three_utr_start)
                    three_utr_end_list.append(three_utr_end)
                    three_utr_num_list.append(e_num)

                    ####################################################################################


                    stop_codon_flag = 1


                if int(self.cds_start) > e_start and int(self.cds_start) <= e_end:  # cds start (start codon) is within this exon
                    start_codon_start = self.cds_start
                    start_codon_end = int(self.cds_start) + 2

                    cds_e_start = start_codon_start # pos cds start codon, cds includes it

                    five_utr_start = e_start
                    five_utr_end = start_codon_start - 1
                    five_utr_start_list.append(five_utr_start)
                    five_utr_end_list.append(five_utr_end)
                    five_utr_num_list.append(e_num)

                    start_codon_flag = 1


                if int(self.cds_start) == int(self.t_start) :  # no start codon because cds begins at beginning of transcript
                    start_codon_flag = 0

                # CDS region starts or ends in this exon
                region_type = "CDS"
                frame = "."

                # check that the cds start is not after the cds end #############
                if cds_e_end >= cds_e_start:
                    outline = format_exon_line(region_type, frame, self.chrom, str(cds_e_start), str(cds_e_end), self.strand,self.gene_id, self.trans_id, source, e_num, self.prot_id,self.degrade_flag, self.match_flag, self.nmd_flag)
                    outfile.write(outline)
                    outfile.write("\n")

                if start_codon_flag == 1:
                    region_type = "start_codon"
                    frame = "."
                    outline = format_exon_line(region_type, frame, self.chrom, str(start_codon_start),str(start_codon_end), self.strand, self.gene_id, self.trans_id, source,e_num, self.prot_id, self.degrade_flag, self.match_flag, self.nmd_flag)
                    outfile.write(outline)
                    outfile.write("\n")

                if stop_codon_flag == 1:
                    region_type = "stop_codon"
                    frame = "."
                    outline = format_exon_line(region_type, frame, self.chrom, str(stop_codon_start),str(stop_codon_end), self.strand, self.gene_id, self.trans_id, source,e_num, self.prot_id, self.degrade_flag, self.match_flag, self.nmd_flag)
                    outfile.write(outline)
                    outfile.write("\n")

                
                
            elif self.strand == "-":
                e_index = int(self.num_exons) - i - 1

                ###########################################################################################

                e_num = i + 1

                e_start = int(self.start_list[e_index])
                e_end = int(self.end_list[e_index])

                region_type = "exon"
                frame = "."

                outline = format_exon_line(region_type, frame, self.chrom, str(e_start), str(e_end), self.strand,self.gene_id, self.trans_id, source, e_num, self.prot_id, self.degrade_flag,self.match_flag, self.nmd_flag)
                outfile.write(outline)
                outfile.write("\n")

                ####################################start here!!!!!

                if self.cds_start == 1 and self.cds_end == 0: # no CDS region
                    continue

                if e_start > int(self.cds_end):  # exon is completely 5' utr on neg strand
                    five_utr_start_list.append(e_start)
                    five_utr_end_list.append(e_end)
                    five_utr_num_list.append(e_num)
                    continue

                if int(self.cds_start) > e_end:  # exon is completely 3' utr for neg strand
                    three_utr_start = e_start
                    three_utr_start_list.append(three_utr_start)
                    three_utr_end = e_end
                    three_utr_end_list.append(three_utr_end)
                    three_utr_num_list.append(e_num)
                    continue


                if e_start > int(self.cds_start) and e_end < int(self.cds_end):  # exon is within coding region but does not contain start or end
                    region_type = "CDS"
                    frame = "."
                    outline = format_exon_line(region_type, frame, self.chrom, str(e_start), str(e_end), self.strand,self.gene_id, self.trans_id, source, e_num, self.prot_id, self.degrade_flag, self.match_flag, self.nmd_flag)
                    outfile.write(outline)
                    outfile.write("\n")
                    continue

                ## all following conditions have cds start/stop regions in exon

                stop_codon_flag = 0
                start_codon_flag = 0

                cds_e_start = e_start
                cds_e_end = e_end

                if e_start < int(self.cds_end) and e_end >= int(self.cds_end):  # cds end is within this exon neg strand
                    cds_e_end = int(self.cds_end)  # neg cds start codon, cds includes it


                    start_codon_start = int(self.cds_end) - 2
                    start_codon_end = int(self.cds_end)

                    five_utr_start = start_codon_end + 1
                    five_utr_start_list.append(five_utr_start)
                    five_utr_end = e_end
                    five_utr_end_list.append(five_utr_end)
                    five_utr_num_list.append(e_num)

                    start_codon_flag = 1


                if int(self.cds_start) > e_start and int(self.cds_start) <= e_end:  # cds start (stop codon) is within this exon
                    stop_codon_start = self.cds_start
                    stop_codon_end = int(self.cds_start) + 2

                    cds_e_start = stop_codon_end + 1

                    three_utr_start = e_start
                    three_utr_end = stop_codon_start
                    three_utr_start_list.append(three_utr_start)
                    three_utr_end_list.append(three_utr_end)
                    three_utr_num_list.append(e_num)

                    stop_codon_flag = 1


                if int(self.cds_end) == int(self.t_end) :  # no start codon because cds begins at beginning of transcript
                    start_codon_flag = 0

                # CDS region starts or ends in this exon
                region_type = "CDS"
                frame = "."
                # check that the cds start is not after the cds end #############
                if cds_e_end >= cds_e_start:
                    outline = format_exon_line(region_type, frame, self.chrom, str(cds_e_start), str(cds_e_end), self.strand,self.gene_id, self.trans_id, source, e_num, self.prot_id,self.degrade_flag, self.match_flag, self.nmd_flag)
                    outfile.write(outline)
                    outfile.write("\n")

                if start_codon_flag == 1:
                    region_type = "start_codon"
                    frame = "."
                    outline = format_exon_line(region_type, frame, self.chrom, str(start_codon_start),str(start_codon_end), self.strand, self.gene_id, self.trans_id, source,e_num, self.prot_id, self.degrade_flag, self.match_flag, self.nmd_flag)
                    outfile.write(outline)
                    outfile.write("\n")

                if stop_codon_flag == 1:
                    region_type = "stop_codon"
                    frame = "."
                    outline = format_exon_line(region_type, frame, self.chrom, str(stop_codon_start),str(stop_codon_end), self.strand, self.gene_id, self.trans_id, source,e_num, self.prot_id, self.degrade_flag, self.match_flag, self.nmd_flag)
                    outfile.write(outline)
                    outfile.write("\n")


                ###########################################################################################



            else:
                print("Error with strand")
                sys.exit()


        for i in xrange(len(five_utr_start_list)):
            five_utr_start = five_utr_start_list[i]
            five_utr_end = five_utr_end_list[i]
            five_utr_e_num = five_utr_num_list[i]

            region_type = "five_prime_utr"
            frame = "."

            if five_utr_end >= five_utr_start:
                outline = format_exon_line(region_type, frame, self.chrom, str(five_utr_start), str(five_utr_end), self.strand,self.gene_id, self.trans_id, source, five_utr_e_num, self.prot_id, self.degrade_flag,self.match_flag, self.nmd_flag)
                outfile.write(outline)
                outfile.write("\n")

        for i in xrange(len(three_utr_start_list)):
            three_utr_start = three_utr_start_list[i]
            three_utr_end = three_utr_end_list[i]
            three_utr_e_num = three_utr_num_list[i]

            region_type = "three_prime_utr"
            frame = "."
            if three_utr_end >= three_utr_start:
                outline = format_exon_line(region_type, frame, self.chrom, str(three_utr_start), str(three_utr_end), self.strand,self.gene_id, self.trans_id, source, three_utr_e_num, self.prot_id, self.degrade_flag,self.match_flag, self.nmd_flag)
                outfile.write(outline)
                outfile.write("\n")
            
            
            
        


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
    
    cds_start = line_split[6]
    cds_end = line_split[7]
    
    id_split = id_line.split(";")
    gene_id = id_split[0]
    trans_id = id_split[1]

    if gene_id not in gene_trans_dict:
        gene_trans_dict[gene_id] = []
        gene_list.append(gene_id)
    
    gene_trans_dict[gene_id].append(trans_id)
    
    trans_dict[trans_id] = Transcript(line)
    
source = "PBRI"

def format_gene_line(chrom,gene_start,gene_end,strand,gene_id,source):
    #anno_line = "gene_id \"" + gene_id + "\"; gene_source \"" + source + "\";"
    #anno_line = "gene_id \"" + gene_id + "\";"
    # convert 0 coord bed to 1 coord gtf
    gene_start = gene_start
    
    anno_list = []
    gene_id_line = "gene_id \"" + gene_id + "\""
    anno_list.append(gene_id_line)
    
    gene_source_line = "gene_source \"" + gene_source + "\""
    anno_list.append(gene_source_line)
    
    #gene_biotype = "misc"
    #gene_biotype_line = "gene_biotype \"" + gene_biotype + "\""
    #anno_list.append(gene_biotype_line)
    
    anno_list.append("")
    
    anno_line = "; ".join(anno_list)
    
    outline = "\t".join([chrom,source,"gene",str(gene_start),str(gene_end),".",strand,".",anno_line])
    return outline

def format_trans_line(chrom,t_start,t_end,strand,gene_id,trans_id,source):
    # convert 0 coord bed to 1 coord gtf
    t_start = int(t_start)
    
    ################################################
    
    anno_list = []
    gene_id_line = "gene_id \"" + gene_id + "\""
    anno_list.append(gene_id_line)
    
    trans_id_line = "transcript_id \"" + trans_id + "\""
    anno_list.append(trans_id_line)
    
    gene_source_line = "gene_source \"" + gene_source + "\""
    anno_list.append(gene_source_line)
    
    trans_source_line = "transcript_source \"" + trans_source + "\""
    anno_list.append(trans_source_line)
    
    #gene_biotype = "misc"
    #gene_biotype_line = "gene_biotype \"" + gene_biotype + "\""
    #anno_list.append(gene_biotype_line)
    
    anno_list.append("")
    
    anno_line = "; ".join(anno_list)
    
    ################################################
    
    #anno_line = "gene_id \"" + gene_id + "\"; transcript_id \"" + trans_id + "\"; uniq_trans_id \"" + uniq_trans_id + "\";"
    outline = "\t".join([chrom,source,"transcript",str(t_start),str(t_end),".",strand,".",anno_line])
    return outline

def format_exon_line(region_type, frame,chrom,e_start,e_end,strand,gene_id,trans_id,source,e_num,prot_id,degrade_flag,match_flag,nmd_flag):
    # convert 0 coord bed to 1 coord gtf
    e_start = int(e_start)
    
    ################################################
    
    anno_list = []
    gene_id_line = "gene_id \"" + gene_id + "\""
    anno_list.append(gene_id_line)
    
    trans_id_line = "transcript_id \"" + trans_id + "\""
    anno_list.append(trans_id_line)
    
    exon_num = str(e_num)
    exon_num_line = "exon_number \"" + exon_num + "\""
    anno_list.append(exon_num_line)
    
    gene_source_line = "gene_source \"" + gene_source + "\""
    anno_list.append(gene_source_line)
    
    trans_source_line = "transcript_source \"" + trans_source + "\""
    anno_list.append(trans_source_line)
    
    prot_id_line = "prot_id \"" + prot_id + "\""
    anno_list.append(prot_id_line)
    
    degrade_flag_line = "degrade_flag \"" + degrade_flag + "\""
    anno_list.append(degrade_flag_line)
    
    match_flag_line = "match_flag \"" + match_flag + "\""
    anno_list.append(match_flag_line)
    
    nmd_flag_line = "nmd_flag \"" + nmd_flag + "\""
    anno_list.append(nmd_flag_line)

    
    #gene_biotype = "misc"
    #gene_biotype_line = "gene_biotype \"" + gene_biotype + "\""
    #anno_list.append(gene_biotype_line)
    
    anno_list.append("")
    
    anno_line = "; ".join(anno_list)
    
    ################################################
    
    #anno_line = "gene_id \"" + gene_id + "\"; transcript_id \"" + trans_id + "\"; exon_number \"" + str(e_num) +"\"; uniq_trans_id \"" + uniq_trans_id + "\";"
    outline = "\t".join([chrom,source,region_type,str(e_start),str(e_end),".",strand,frame,anno_line])
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

        trans_dict[trans_id].write_trans()


            
            
        
        

