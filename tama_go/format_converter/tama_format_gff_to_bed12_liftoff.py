import re
import sys
import time


## from Bio import SeqIO



#
# This script converts  gff to bed format for a liftoff gff3 file
# 2022_04_20 Richard Kuo


#

print("opening gtf file")
gtf_file = sys.argv[1]
gtf_file_contents = open(gtf_file).read().rstrip("\n").split("\n")

outfile_name = sys.argv[2]
outfile = open(outfile_name,"w")



class Gene:
    def __init__(self,gene_line):

        #ChrC    TAIR10  gene    107599  107701  .       +       .       ID=ATCG00960;Name=ATCG00960;Name=ATCG00960

        line_split = gene_line.split("\t")
        anno_line = line_split[8]
        anno_split = anno_line.split(";")
        
        for subfield in anno_split:
            subfield_split = subfield.split("=")
            if subfield_split[0] == "ID":
                self.gene_id = subfield.split("=")[1]
        
        self.chrom = line_split[0]
        self.g_start = int(line_split[3]) - 1 
        self.g_end = int(line_split[4])
        self.strand = line_split[6]
        
        self.gene_type = line_split[2]
        

# all coords are non relative. need to convert to relative coords for exons
class Transcript:
    def __init__(self,trans_line):

        #ChrC    TAIR10  transcript      107599  107701  .       +       .       ID=ATCG00960.1;Parent=ATCG00960

        line_split = trans_line.split("\t")
        anno_line = line_split[8]
        anno_split = anno_line.split(";")
        
        region_type = line_split[2]
        
        self.gene_id = "NA"
        self.trans_id = "NA"
        self.trans_name = "NA"
        self.trans_type = "NA"

        self.trans_type = region_type
        
        for subfield in anno_split:
            subfield_split = subfield.split("=")
            if subfield_split[0] == "Parent":
                self.gene_id = subfield.split("=")[1]
            if subfield_split[0] == "ID":
                self.trans_id = subfield.split("=")[1]
            if subfield_split[0] == "Name":
                self.trans_name = subfield.split("=")[1]

        # correct for pseudogenes which do not have transcript lines
        if region_type == "exon":
            if gene_dict[self.gene_id].gene_type == "pseudogene":
                self.trans_id = self.gene_id
            #if "gene_biotype" in subfield:
            #    self.gene_class = subfield.split("\"")[1]
            #if "transcript_biotype" in subfield:
            #    self.trans_class = subfield.split("\"")[1]
        
        if self.gene_id == "NA":
            print("Error with gene id")
            print(trans_line)
            
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
       #
        #self.start_codon = "none"
        #self.stop_codon = "none"
        
        self.utr_list = []
        
        
    def add_exon(self,exon_line):
        exon_obj = Exon(exon_line)
        
        #make sure exon belongs to this transcript
        #if exon_obj.trans_id != self.trans_id:
        #if self.trans_id not in exon_obj.trans_id:
        #    print("Exon does not match transcript!!!")
        #    print(exon_obj.trans_id)
        #    print(self.trans_id)
        #    sys.exit()
        
        
        self.e_start_list.append(exon_obj.e_start)
        self.e_end_list.append(exon_obj.e_end)
        self.e_obj_list.append(exon_obj)
    
    def add_cds(self,cds_line):
        cds_obj = Cds(cds_line)
        
        #make sure exon belongs to this transcript
        if cds_obj.trans_id != self.trans_id:
            print("CDS does not match transcript!!!")
            sys.exit()
        
        
        self.c_start_list.append(cds_obj.c_start)
        self.c_end_list.append(cds_obj.c_end)
        self.c_obj_list.append(cds_obj)

        self.c_start_list.sort()
        self.c_end_list.sort()

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
        #ChrC    TAIR10  exon    107599  107701  .       +       .       ID=exon:ATCG00960.1:1;Parent=ATCG00960.1

        line_split = exon_line.split("\t")
        anno_line = line_split[8]
        anno_split = anno_line.split(";")
        
        region_type = line_split[2]
        if region_type != "exon":
            print("Error with region type exon")
            print(exon_line)
            sys.exit()

        for subfield in anno_split:
            subfield_split = subfield.split("=")
            
            if subfield_split[0] == "Parent":
                trans_id = subfield.split("=")[1]
            if subfield_split[0] == "ID":
                e_id = subfield.split("=")[1]

        self.trans_id = trans_id
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
            subfield_split = subfield.split("=")

            if subfield_split[0] == "Parent":
                self.trans_id = subfield.split("=")[1].split(",")[0]


        
        #self.c_start = int(line_split[3]) - 1
        self.c_start = int(line_split[3]) - 1 # adjust from gff 1 base to bed 0 base
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
        if region_type != "five_prime_UTR" and region_type != "three_prime_UTR":
            print("Error with region type UTR")
            print(utr_line)
            sys.exit()

        for subfield in anno_split:
            subfield_split = subfield.split("=")

            if subfield_split[0] == "Parent":
                self.trans_id = subfield.split("=")[1].split(",")[0]


        
        #self.u_start = int(line_split[3]) - 1
        self.u_start = int(line_split[3]) - 1 # adjust from gff 1 base to bed 0 base
        self.u_end = int(line_split[4])
        self.strand = line_split[6]


gene_trans_dict = {} # gene_trans_exon_dict[gene id][trans id] = 1
gene_dict = {} # gene_dict[gene id] = gene object
trans_dict = {} # trans_dict[trans id] = trans object
exon_dict = {} # exon_dict[exon id] = exon object

trans_list = [] # keep ordering of transcripts


print("Going through GFF genes")
for line in gtf_file_contents:
    if line.startswith("#"):
        continue
    
    line_split = line.split("\t")
    region_type = line_split[2]
    
    if region_type == "gene" or  region_type == "pseudogene":
        
        anno_line = line_split[8]
        anno_split = anno_line.split(";")
            
        for subfield in anno_split:
            subfield_split = subfield.split("=")
            if subfield_split[0] == "ID":
                gene_id = subfield.split("=")[1]
    
        gene_dict[gene_id] = Gene(line)
        


#use this to match transcript types from NCBI
trans_type_dict = {} # trans_type_dict[trans type] = 1

trans_type_dict["transcript"] = 1
trans_type_dict["mRNA"] = 1
trans_type_dict["ncRNA"] = 1
trans_type_dict["tRNA"] = 1
trans_type_dict["rRNA"] = 1
trans_type_dict["snRNA"] = 1
trans_type_dict["miRNA"] = 1
trans_type_dict["primary_transcript"] = 1
trans_type_dict["snoRNA"] = 1
trans_type_dict["V_gene_segment"] = 1
trans_type_dict["guide_RNA"] = 1
trans_type_dict["C_gene_segment"] = 1
trans_type_dict["telomerase_RNA"] = 1
trans_type_dict["SRP_RNA"] = 1
trans_type_dict["lnc_RNA"] = 1

region_type_dict = {} # region_type_dict[region type] = 1
#region_type_dict["gene"] = 1
#region_type_dict["mRNA"] = 1
region_type_dict["exon"] = 1
region_type_dict["CDS"] = 1
region_type_dict["five_prime_UTR"] = 1
region_type_dict["three_prime_UTR"] = 1



print("Going through GFF transcripts")
for line in gtf_file_contents:
    if line.startswith("#"):
        continue
    
    line_split = line.split("\t")
    region_type = line_split[2]
    
    if region_type in trans_type_dict:
        trans_flag = 1
    else:
        continue
    
    anno_line = line_split[8]
    anno_split = anno_line.split(";")
        
    for subfield in anno_split:
        subfield_split = subfield.split("=")
        if subfield_split[0] == "Parent":
            gene_id = subfield.split("=")[1]
        if subfield_split[0] == "ID":
            trans_id = subfield.split("=")[1]
            

    trans_dict[trans_id] = Transcript(line)
    if gene_id not in gene_trans_dict:
        gene_trans_dict[gene_id] = {}
    
    if trans_id not in gene_trans_dict[gene_id]:
        gene_trans_dict[gene_id][trans_id] = {}
        trans_list.append(trans_id)

print("Going through GFF exons/cds/start codon/stop codon/utr")
for line in gtf_file_contents:
    if line.startswith("#"):
        continue
    
    line_split = line.split("\t")
    region_type = line_split[2]
    
    #if region_type == "gene" or region_type == "transcript" or region_type == "mRNA":
    #    continue

    if region_type not in region_type_dict:
        #print("error with unrecognized region type")
        #print(region_type)
        #sys.exit()
        continue

    anno_line = line_split[8]
    anno_split = anno_line.split(";")


#Chr1_RagTag_polished    Liftoff gene    6983    9251    .       +       .       ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010;coverage=1.0;sequence_ID=1.0;valid_ORFs=1;extra_copy_number=0;copy_num_ID=AT1G01010_0
#Chr1_RagTag_polished    Liftoff mRNA    6983    9251    .       +       .       ID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1;Index=1;matches_ref_protein=True;valid_ORF=True;extra_copy_number=0
#Chr1_RagTag_polished    Liftoff exon    6983    7265    .       +       .       ID=exon_1;Parent=AT1G01010.1;extra_copy_number=0
#Chr1_RagTag_polished    Liftoff exon    7348    7628    .       +       .       ID=exon_2;Parent=AT1G01010.1;extra_copy_number=0
#Chr1_RagTag_polished    Liftoff exon    7838    7957    .       +       .       ID=exon_3;Parent=AT1G01010.1;extra_copy_number=0
#Chr1_RagTag_polished    Liftoff exon    8058    8447    .       +       .       ID=exon_4;Parent=AT1G01010.1;extra_copy_number=0
#Chr1_RagTag_polished    Liftoff exon    8526    8678    .       +       .       ID=exon_5;Parent=AT1G01010.1;extra_copy_number=0
#Chr1_RagTag_polished    Liftoff exon    8791    9251    .       +       .       ID=exon_6;Parent=AT1G01010.1;extra_copy_number=0
#Chr1_RagTag_polished    Liftoff CDS     7112    7265    .       +       .       ID=CDS_1;Parent=AT1G01010.1,AT1G01010.1-Protein;extra_copy_number=0
#Chr1_RagTag_polished    Liftoff CDS     7348    7628    .       +       .       ID=CDS_2;Parent=AT1G01010.1,AT1G01010.1-Protein;extra_copy_number=0
#Chr1_RagTag_polished    Liftoff CDS     7838    7957    .       +       .       ID=CDS_3;Parent=AT1G01010.1,AT1G01010.1-Protein;=;extra_copy_number=0
#Chr1_RagTag_polished    Liftoff CDS     8058    8447    .       +       .       ID=CDS_4;Parent=AT1G01010.1,AT1G01010.1-Protein;=;extra_copy_number=0
#Chr1_RagTag_polished    Liftoff CDS     8526    8678    .       +       .       ID=CDS_5;Parent=AT1G01010.1,AT1G01010.1-Protein;=;extra_copy_number=0
#Chr1_RagTag_polished    Liftoff CDS     8791    8982    .       +       .       ID=CDS_6;Parent=AT1G01010.1,AT1G01010.1-Protein;=;extra_copy_number=0
#Chr1_RagTag_polished    Liftoff five_prime_UTR  6983    7111    .       +       .       ID=five_prime_UTR_1;Parent=AT1G01010.1;extra_copy_number=0
#Chr1_RagTag_polished    Liftoff three_prime_UTR 8983    9251    .       +       .       ID=three_prime_UTR_1;Parent=AT1G01010.1;extra_copy_number=0



    for subfield in anno_split:
        subfield_split = subfield.split("=")
        
        if subfield_split[0] == "Parent":
            trans_id = subfield.split("=")[1].split(",")[0]


    if region_type == "exon":

        #print(line)
        if trans_id in trans_dict:
            trans_dict[trans_id].add_exon(line)
        else:
            trans_dict[trans_id] = Transcript(line)
            trans_dict[trans_id].add_exon(line)

    elif region_type == "CDS":
        trans_dict[trans_id].add_cds(line)
    #elif region_type == "start_codon":
    #    trans_dict[trans_id].add_start(line)
    #elif region_type == "stop_codon":
    #    trans_dict[trans_id].add_stop(line)
    elif region_type == "five_prime_UTR":
        trans_dict[trans_id].add_utr(line)
    elif region_type == "three_prime_UTR":
        trans_dict[trans_id].add_utr(line)
    #elif region_type == "Selenocysteine":
    #    continue
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
    trans_name = trans_dict[trans_id].trans_name
    trans_type = trans_dict[trans_id].trans_type
    #gene_class = trans_dict[trans_id].gene_class
    #trans_class = trans_dict[trans_id].trans_class
    #id_line = ";".join([gene_id,trans_id,gene_class,trans_class])
    id_line = ";".join([gene_id, trans_id,trans_name,trans_type])
    strand = trans_dict[trans_id].strand
    blocks_line = ",".join(e_block_sizes)
    starts_line = ",".join(rel_starts)

    if len(trans_dict[trans_id].c_start_list) > 0 :
        cds_start = str(trans_dict[trans_id].c_start_list[0])
        cds_end = str(trans_dict[trans_id].c_end_list[-1])
    else:
        cds_start = "0"
        cds_end = "0"

    ##################
    # blocked out lines used to go here but I moved them to the end for ease of reading
    # may want to use them again for CDS stuff
    ##################



    outline = "\t".join([chrom,t_start,t_end,id_line,"40",strand,cds_start,cds_end,"255,0,0",e_num,blocks_line,starts_line])
    outfile.write(outline)
    outfile.write("\n")
