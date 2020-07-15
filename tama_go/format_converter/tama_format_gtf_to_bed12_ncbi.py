import re
import sys
import time



#
# This script converts NCBI gtf to bed format
# It is designed to work with all current versions.
# It includes CDS boundaries for the 7th and 8th columns of the bed file.
#

print("opening gtf file")
gtf_file = sys.argv[1]
gtf_file_contents = open(gtf_file).read().rstrip("\n").split("\n")

outfile_name = sys.argv[2]
outfile = open(outfile_name,"w")



class Gene:
    def __init__(self,gene_line):
        line_split = gene_line.split("\t")
        anno_line = line_split[8]
        anno_split = anno_line.split(";")
        
        for subfield in anno_split:
            if "gene_id" in subfield:
                self.gene_id = subfield.split("\"")[1]
        
        self.chrom = line_split[0]
        #self.g_start = int(line_split[3]) - 1
        #self.g_end = int(line_split[4])
        self.strand = line_split[6]
        

# all coords are non relative. need to convert to relative coords for exons
class Transcript:
    def __init__(self,trans_line):
        line_split = trans_line.split("\t")
        anno_line = line_split[8]
        anno_split = anno_line.split(";")
        
        for subfield in anno_split:
            if "gene_id" in subfield:
                self.gene_id = subfield.split("\"")[1]
            if "transcript_id" in subfield:
                self.trans_id = subfield.split("\"")[1]
            if "gene_biotype" in subfield:
                self.gene_class = subfield.split("\"")[1]
            if "transcript_biotype" in subfield:
                self.trans_class = subfield.split("\"")[1]
        
        
        self.chrom = line_split[0]
        self.t_start = int(line_split[3]) - 1
        self.t_end = int(line_split[4])
        self.strand = line_split[6]
        
        #self.gene_class = line_split[1]
        
        self.e_start_list = []
        self.e_end_list = []
        self.e_obj_list = []
        #self.e_block_size_list = []
        
        #cds region
        self.c_start_list = []
        self.c_end_list = []
        self.c_obj_list = []
        
        self.start_codon = "none"
        self.stop_codon = "none"
        
        self.utr_list = []
        
        
    def add_exon(self,exon_line):
        exon_obj = Exon(exon_line)
        
        #make sure exon belongs to this transcript
        if exon_obj.trans_id != self.trans_id:
            print("Exon does not match transcript!!!")
            sys.exit()
        
        
        self.e_start_list.append(exon_obj.e_start)
        self.e_start_list.sort()
        self.e_end_list.append(exon_obj.e_end)
        self.e_end_list.sort()
        self.e_obj_list.append(exon_obj)

        if self.t_start > exon_obj.e_start:
            self.t_start = exon_obj.e_start

        if self.t_end < exon_obj.e_end:
            self.t_end = exon_obj.e_end
    
    def add_cds(self,cds_line):
        cds_obj = Cds(cds_line)
        
        #make sure exon belongs to this transcript
        if cds_obj.trans_id != self.trans_id:
            print("CDS does not match transcript!!!")
            sys.exit()
        
        
        self.c_start_list.append(cds_obj.c_start)
        self.c_end_list.append(cds_obj.c_end)
        self.c_obj_list.append(cds_obj)

    def add_start(self,start_line):
        start_obj = StartCodon(start_line)
        
        #make sure exon belongs to this transcript
        if start_obj.trans_id != self.trans_id:
            print("CDS does not match transcript!!!")
            sys.exit()
        
        if self.start_codon == "none":
            self.start_codon = start_obj
        else:
            print("Multiple start codons")
            print(self.trans_id)
            self.start_codon = "multiple"

        
    def add_stop(self,stop_line):
        stop_obj = StopCodon(stop_line)
        
        #make sure exon belongs to this transcript
        if stop_obj.trans_id != self.trans_id:
            print("CDS does not match transcript!!!")
            sys.exit()
        
        if self.stop_codon == "none":
            self.stop_codon = stop_obj
        else:
            print("Multiple stop codons")
            print(self.trans_id)
            self.stop_codon = "multiple"
            
        
    def add_utr(self,utr_line):
        utr_obj = Utr(utr_line)
        
        #make sure exon belongs to this transcript
        if utr_obj.trans_id != self.trans_id:
            print("CDS does not match transcript!!!")
            sys.exit()
        
        self.utr_list.append(utr_obj)



      
class Exon:
    def __init__(self,exon_line):
        line_split = exon_line.split("\t")
        anno_line = line_split[8]
        anno_split = anno_line.split(";")
        
        region_type = line_split[2]
        if region_type != "exon":
            print("Error with region type exon")
            print(exon_line)
            sys.exit()
        
        for subfield in anno_split:
            if "gene_id" in subfield:
                self.gene_id = subfield.split("\"")[1]
            if "transcript_id" in subfield:
                self.trans_id = subfield.split("\"")[1]
            if "exon_number" in subfield:
                self.exon_num = subfield.split("\"")[1]
        
        self.e_start = int(line_split[3]) - 1
        self.e_end = int(line_split[4])
        self.strand = line_split[6]
        
class Cds:
    def __init__(self,cds_line):
        line_split = cds_line.split("\t")
        anno_line = line_split[8]
        anno_split = anno_line.split(";")
        
        region_type = line_split[2]
        if region_type != "CDS":
            print("Error with region type CDS")
            print(cds_line)
            sys.exit()
        
        for subfield in anno_split:
            if "gene_id" in subfield:
                self.gene_id = subfield.split("\"")[1]
            if "transcript_id" in subfield:
                self.trans_id = subfield.split("\"")[1]
        
        #self.c_start = int(line_split[3]) - 1
        self.c_start = int(line_split[3])
        self.c_end = int(line_split[4])
        self.strand = line_split[6]
        
class StartCodon:
    def __init__(self,start_line):
        line_split = start_line.split("\t")
        anno_line = line_split[8]
        anno_split = anno_line.split(";")
        
        region_type = line_split[2]
        if region_type != "start_codon":
            print("Error with region type start_codon")
            print(start_line)
            sys.exit()
        
        for subfield in anno_split:
            if "gene_id" in subfield:
                self.gene_id = subfield.split("\"")[1]
            if "transcript_id" in subfield:
                self.trans_id = subfield.split("\"")[1]
        
        #self.bc_start = int(line_split[3]) - 1
        self.bc_start = int(line_split[3])
        self.bc_end = int(line_split[4])
        self.strand = line_split[6]

class StopCodon:
    def __init__(self,stop_line):
        line_split = stop_line.split("\t")
        anno_line = line_split[8]
        anno_split = anno_line.split(";")
        
        region_type = line_split[2]
        if region_type != "stop_codon":
            print("Error with region type stop_codon")
            print(stop_line)
            sys.exit()
        
        for subfield in anno_split:
            if "gene_id" in subfield:
                self.gene_id = subfield.split("\"")[1]
            if "transcript_id" in subfield:
                self.trans_id = subfield.split("\"")[1]
        
        #self.ec_start = int(line_split[3]) - 1
        self.ec_start = int(line_split[3])
        self.ec_end = int(line_split[4])
        self.strand = line_split[6]
        
class Utr:
    def __init__(self,utr_line):
        line_split = utr_line.split("\t")
        anno_line = line_split[8]
        anno_split = anno_line.split(";")
        
        region_type = line_split[2]
        if region_type != "UTR" and "utr" not in region_type:
            print("Error with region type UTR")
            print(utr_line)
            sys.exit()
        
        for subfield in anno_split:
            if "gene_id" in subfield:
                self.gene_id = subfield.split("\"")[1]
            if "transcript_id" in subfield:
                self.trans_id = subfield.split("\"")[1]
        
        #self.u_start = int(line_split[3]) - 1
        self.u_start = int(line_split[3])
        self.u_end = int(line_split[4])
        self.strand = line_split[6]


gene_trans_dict = {} # gene_trans_exon_dict[gene id][trans id] = 1
gene_dict = {} # gene_dict[gene id] = gene object
trans_dict = {} # trans_dict[trans id] = trans object
exon_dict = {} # exon_dict[exon id] = exon object

trans_list = [] # keep ordering of transcripts


print("Going through GTF genes")
for line in gtf_file_contents:
    if line.startswith("#"):
        continue
    
    line_split = line.split("\t")
    region_type = line_split[2]
    
    if region_type != "exon":
        continue
    
    anno_line = line_split[8]
    anno_split = anno_line.split(";")
        
    for subfield in anno_split:
        if "gene_id" in subfield:
            gene_id = subfield.split("\"")[1]

    # continue if we covered this gene already
    if gene_id in gene_dict:
        continue

    gene_dict[gene_id] = Gene(line)

ncbi_trans_dict = {} # ncbi_trans_dict[trans id] = 1   used to check if transcript has already been logged in


print("Going through GTF transcripts")
for line in gtf_file_contents:
    if line.startswith("#"):
        continue
    
    line_split = line.split("\t")
    region_type = line_split[2]
    
    if region_type != "exon":
        continue
    
    anno_line = line_split[8]
    anno_split = anno_line.split(";")
        
    for subfield in anno_split:
        if "gene_id" in subfield:
            gene_id = subfield.split("\"")[1]
        if "transcript_id" in subfield:
            trans_id = subfield.split("\"")[1]

    # skip if this transcript has already been added
    if trans_id in trans_dict:
        continue

    trans_dict[trans_id] = Transcript(line)
    if gene_id not in gene_trans_dict:
        gene_trans_dict[gene_id] = {}
    
    if trans_id not in gene_trans_dict[gene_id]:
        gene_trans_dict[gene_id][trans_id] = {}
        trans_list.append(trans_id)

print("Going through GTF exons/cds/start codon/stop codon/utr")
for line in gtf_file_contents:
    if line.startswith("#"):
        continue
    
    line_split = line.split("\t")
    region_type = line_split[2]
    
    if region_type == "gene" or region_type == "transcript":
        continue
    
    anno_line = line_split[8]
    anno_split = anno_line.split(";")
        
    for subfield in anno_split:
        if "gene_id" in subfield:
            gene_id = subfield.split("\"")[1]
        if "transcript_id" in subfield:
            trans_id = subfield.split("\"")[1]

    # skip if this is an unknown transcript
    if trans_id == "unknown_transcript_1":
        continue

    if region_type == "exon":
        trans_dict[trans_id].add_exon(line)
    elif region_type == "CDS":
        trans_dict[trans_id].add_cds(line)
    elif region_type == "start_codon":
        trans_dict[trans_id].add_start(line)
    elif region_type == "stop_codon":
        trans_dict[trans_id].add_stop(line)
    elif region_type == "UTR" or "utr" in region_type:
        trans_dict[trans_id].add_utr(line)
    elif region_type == "Selenocysteine":
        continue
    else:
        print("Region type unknown")
        print(line)
        sys.exit()

####################################################<<<<< RK left off here. Start here!!! Add prcoessing into bed format!!!!

####################################################
    
print("write bed12 file")        
for trans_id in trans_list:
    e_starts = sorted(trans_dict[trans_id].e_start_list)
    e_ends = sorted(trans_dict[trans_id].e_end_list)
    e_block_sizes = []
    rel_starts = []
    
    
    
    for i,start in enumerate(e_starts):
        block_size = e_ends[i] - e_starts[i]
        e_block_sizes.append(str(block_size))
        
        rel_start = start - trans_dict[trans_id].t_start
        rel_starts.append(str(rel_start))
        
    e_num = str(len(e_block_sizes))
    
    chrom = trans_dict[trans_id].chrom
    t_start = str(trans_dict[trans_id].t_start)
    t_end = str(trans_dict[trans_id].t_end)
    gene_id = trans_dict[trans_id].gene_id
    trans_id = trans_dict[trans_id].trans_id
    # gene_class = trans_dict[trans_id].gene_class
    # trans_class = trans_dict[trans_id].trans_class
    # id_line = ";".join([gene_id,trans_id,gene_class,trans_class])
    id_line = ";".join([gene_id, trans_id])
    strand = trans_dict[trans_id].strand
    blocks_line = ",".join(e_block_sizes)
    starts_line = ",".join(rel_starts)
    
    #get CDS/UTR information for columns 7 and 8
    c_starts = sorted(trans_dict[trans_id].c_start_list)
    c_ends = sorted(trans_dict[trans_id].c_end_list)

    start_codon_obj = trans_dict[trans_id].start_codon
    stop_codon_obj = trans_dict[trans_id].stop_codon
    
    
    start_codon_flag = 1
    stop_codon_flag = 1
    if start_codon_obj == "none" :
        start_codon_flag = 0
        print("Error missing start codon")
        print(trans_id)
        #sys.exit()
    if start_codon_obj == "multiple":
        start_codon_flag = 0
        print("Error multiple start codon")
        print(trans_id)
    if stop_codon_obj == "none":
        stop_codon_flag = 0
        print("Error missing stop codon")
        print(trans_id)
        #sys.exit()
    if stop_codon_obj == "multiple":
        stop_codon_flag = 0
        print("Error multiple stop codon")
        print(trans_id)
    
    #check UTR against start/stop codons
    utr_list = trans_dict[trans_id].utr_list
    u_start_list = []
    u_end_list = []
    for utr_obj in utr_list:
        u_start_list.append(utr_obj.u_start)
        u_end_list.append(utr_obj.u_end)
        
    u_start_list.sort()
    u_end_list.sort()
    
    
    # clean up utr's which are on outskirts and match exons
#    for i,u_start in enumerate(u_start_list):
#        u_end = u_end_list[i]
#        for j,e_start in enumerate(e_starts):
#            e_start = e_start + 1 #correct for bed format conversion done in exon class
#            e_end = e_ends[j]                
#            if int(u_start) == int(e_start) and int(u_end) == int(e_end):
#                u_start_list.pop(i)
#                u_end_list.pop(i)
    
    # clean up utr's which are on outskirts and match exons
    u_index = 0
    while u_index < len(u_start_list) :
        u_start = u_start_list[u_index]
        u_end = u_end_list[u_index]
        
        for j,e_start in enumerate(e_starts):
            e_start = e_start + 1 #correct for bed format conversion done in exon class
            e_end = e_ends[j]                
            if int(u_start) == int(e_start) and int(u_end) == int(e_end):
                u_start_list.pop(u_index)
                u_end_list.pop(u_index)
                if u_index < len(u_start_list):
                    u_start = u_start_list[u_index]
                    u_end = u_end_list[u_index]
                    continue
                else:
                    break
            
                
                
            
        u_index += 1

    
    utr_flag = 0 #if utr flag == 1 then utr check will be done elsewise no
    # utr number check
    if len(u_start_list) == 2:
        utr_flag = 1
    elif len(u_start_list) == 0:
        utr_flag = 0
    else:
        print("Error with UTR numbers")
        print(trans_id)
        print(len(u_start_list))
        #print(u_start_list)
        #print(u_end_list)
        #print(e_starts)
        #print(e_ends)
        #sys.exit()
    
    bed_cds_start = 0
    bed_cds_end = 0
    
    cds_flag = 1
    if len(c_starts) == 0:
        cds_flag = 0
        print("no cds")
        print(trans_id)
        #sys.exit()

    if cds_flag == 1:
        #cds versus start/stop codons check and utr versus cds 
        if strand == "+":
            cds_start = c_starts[0]
            cds_end = c_ends[-1]
            
            if start_codon_flag == 1:
                start_codon_coord = start_codon_obj.bc_start
                if cds_start != start_codon_coord:
                    print("Error with cds versus start codon")
                    print(trans_id)
                    #sys.exit()
            
            if stop_codon_flag == 1:
                stop_codon_coord = stop_codon_obj.ec_start
                if cds_end + 1 != stop_codon_coord:
                    print("Error with cds versus stop codon")
                    print(trans_id)
                    print(cds_end)
                    print(stop_codon_coord)
                    #sys.exit()
                
            #utr check
            # stop codon adjustment value
            if stop_codon_flag == 0:
                stop_codon_adjust = 0
            elif stop_codon_flag == 1:
                stop_codon_adjust = 3
                
            if utr_flag == 1:
                first_utr_end = u_end_list[0]
                last_utr_start =  u_start_list[1]
                
                if first_utr_end + 1 != cds_start:
                    print("Error with utr versus cds_start")
                    print(trans_id)
                    print(first_utr_end)
                    print(cds_start)
                    print(u_start_list)
                    print(u_end_list)
                    print(e_starts)
                    print(e_ends)
                    #sys.exit()
                if last_utr_start - stop_codon_adjust != cds_end + 1:
                    print("Error with utr versus cds_end")
                    print(trans_id)
                    print(str(last_utr_start) + "\t" + str(cds_end))
                    #sys.exit()
            
            bed_cds_start = cds_start
            # Note that we are not counting stop codons here!!!
            bed_cds_end = cds_end  
            #bed_cds_end = cds_end + 3 + 1 # note that indexing here is not corrected for by bed format so we add a one for the later coord
            #bed_cds_end = stop_codon_obj.ec_end + 1 # note that indexing here is not corrected for by bed format so we add a one for the later coord
            
        elif strand == "-":
            cds_start = c_ends[-1]
            cds_end = c_starts[0]
            
            if start_codon_flag == 1:
                start_codon_coord = start_codon_obj.bc_end
                
                if cds_start != start_codon_coord:
                    print("Error with cds versus start codon")
                    print(trans_id)
                    #sys.exit()
            
            if stop_codon_flag == 1:
                stop_codon_coord = stop_codon_obj.ec_end
                if cds_end - 1 != stop_codon_coord:
                    print("Error with cds versus stop codon")
                    print(trans_id)
                    print(cds_end)
                    print(stop_codon_coord)
                    #sys.exit()
                
            #utr check
            # stop codon adjustment value
            if stop_codon_flag == 0:
                stop_codon_adjust = 0
            elif stop_codon_flag == 1:
                stop_codon_adjust = 3

            if utr_flag == 1:
                first_utr_end = u_start_list[1]
                last_utr_start = u_end_list[0] 
                
                if first_utr_end - 1 != cds_start:
                    print("Error with utr versus cds_start")
                    print(trans_id)
                    print(first_utr_end)
                    print(cds_start)
                    #sys.exit()
            
                if last_utr_start + stop_codon_adjust != cds_end - 1:
                    print("Error with utr versus cds_end")
                    print(trans_id)
                    print(last_utr_start)
                    print(cds_end)
                    #sys.exit()
            
            #bed_cds_start = c_starts[0] - 3
            bed_cds_start = c_starts[0] # Note that we are not counting stop codons here!!!
            bed_cds_end = c_ends[-1]  
     
            #bed_cds_start = stop_codon_obj.ec_start
            #bed_cds_end = start_codon_obj.bc_end + 1 # note that indexing here is not corrected for by bed format so we add a one for the later coord
  
        #rel_bed_cds_start = bed_cds_start - int(t_start)
        #rel_bed_cds_end = bed_cds_end - int(t_start) + 1 # note that indexing here is not corrected for by bed format so we add a one for the later coord
        bed_cds_start = bed_cds_start - 1 #convert to bed format indexing rules
    
    elif cds_flag == 0: # not a coding gene!
        bed_cds_start = t_start
        bed_cds_end = t_start
    else:
        print("Error with cds flag!")
        print(trans_id)
        sys.exit()

    outline = "\t".join([chrom,t_start,t_end,id_line,"40",strand,str(bed_cds_start),str(bed_cds_end),"255,0,0",e_num,blocks_line,starts_line])
    outfile.write(outline)
    outfile.write("\n")
    
    
        

    