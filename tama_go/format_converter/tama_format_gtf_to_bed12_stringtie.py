import re
import sys
import time



#
# This script converts cufflinks/stringtie gtf into bed format file
#


print("opening gtf file")
gtf_file = sys.argv[1]
gtf_file_contents = open(gtf_file).read().rstrip("\n").split("\n")


outfile_name = sys.argv[2]
outfile = open(outfile_name,"w")


class Transcript:
    def __init__(self,trans_id,gene_id,chrom,strand):
        self.trans_id = trans_id
        self.gene_id = gene_id
        self.chrom = chrom
        self.t_start = 0
        self.t_end = 0
        self.strand = strand
    
        
        self.bed_starts = []
        self.bed_blocks = []
        
        self.bed_starts_string = ""
        self.bed_blocks_string = ""

        
        #coordinate starts and ends
        self.start_list = []
        self.end_list =  []
        self.num_exons = 0
        self.exon_dict = {} # exon_dict[exon num][start/end] = start/end 
        
        
    def add_exon(self,e_num,e_start,e_end):
        if e_num not in self.exon_dict:
            self.exon_dict[e_num] = {}
            self.exon_dict[e_num]["start"] = e_start
            self.exon_dict[e_num]["end"] = e_end
        else:
            print("Error with duplicate exon number")
            print(e_num + "\t" + e_start + "\t" +e_end)
            sys.exit()
            
        
    
    def finish_bed(self):
        
        
        exon_num_list = list(self.exon_dict.keys())
        self.num_exons = len(exon_num_list)
        exon_num_list.sort()
        
        last_start = 0 #use to check that exons are in order
        
        for e_num in exon_num_list:
            e_start = self.exon_dict[e_num]["start"]
            e_end = self.exon_dict[e_num]["end"]
            self.start_list.append(e_start)
            self.end_list.append(e_end)
            
            if e_start < last_start:
                print("Error with exons not in order by starts")
                sys.exit()
            last_start = e_start
            
        self.t_start = self.start_list[0] - 1 # adjust for bed indexing
        self.t_end = self.end_list[-1]
        
        for i in range(self.num_exons):
            e_start = self.start_list[i]
            e_rel_start = e_start - self.t_start - 1 # adjust for bed indexing
            
            e_end = self.end_list[i]
            
            e_block_size = e_end + 1 - e_start # add one to get number of bases
            
            self.bed_starts.append(str(e_rel_start))
            self.bed_blocks.append(str(e_block_size))
        
        self.bed_starts_string = ",".join(self.bed_starts)
        self.bed_blocks_string = ",".join(self.bed_blocks)
        


exon_trans_dict = {} # exon_trans_dict[exon id] = trans id

gene_trans_dict = {} # gene_trans_dict[gene_id][trans id] = trans obj

gene_list = []
gene_dict = {} # use this to add non-redundant gene id's to gene list

print("Going through gtf file")
for line in gtf_file_contents:
    if line.startswith("#"):
        continue
    line_split = line.split("\t")
    
    chrom = line_split[0]
    e_start = int(line_split[3])
    e_end = int(line_split[4])
    strand = line_split[6]
    id_line = line_split[8]
    
    id_split = id_line.split(";")
    
    for id_field in id_split:
        if id_field == "":
            continue
        if "\"" not in id_field:
            continue
        
        id_code = id_field.split("\"")[1]
        
        if "gene_id" in id_field:
            gene_id = id_code
            if gene_id not in gene_dict:
                
                gene_dict[gene_id] = []
                
                gene_list.append(gene_id)
                gene_trans_dict[gene_id] = {}
                
                
        elif "transcript_id" in id_field:
            trans_id = id_code
            if trans_id not in gene_dict[gene_id]:
                gene_dict[gene_id].append(trans_id)
                
                gene_trans_dict[gene_id][trans_id] = Transcript(trans_id,gene_id,chrom,strand)
        
        elif "exon_number" in id_field:
            exon_num = int(id_code)
            
            gene_trans_dict[gene_id][trans_id].add_exon(exon_num,e_start,e_end)
            
        

print("writing out file")
for gene_id in gene_list:
    
    trans_list = gene_dict[gene_id]
    
    for trans_id in trans_list:
        
        gene_trans_dict[gene_id][trans_id].finish_bed()
        
        trans_obj = gene_trans_dict[gene_id][trans_id]
        
        chrom = trans_obj.chrom
        t_start = str(trans_obj.t_start)
        t_end = str(trans_obj.t_end)
        
        id_line = ";".join([gene_id,trans_id])
        strand = trans_obj.strand
        num_exons = str(trans_obj.num_exons)
        blocks = trans_obj.bed_blocks_string
        starts = trans_obj.bed_starts_string
        
        bed_line = "\t".join([chrom,t_start,t_end,id_line,"40",strand,t_start,t_end,"255,0,0",num_exons,blocks,starts])
        
        outfile.write(bed_line)
        outfile.write("\n")
    
    
    

    
    
    
