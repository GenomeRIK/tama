import re
import sys
import time
from Bio import SeqIO
import argparse


#
# This script parses information from the default output of blastp
# 
#



ap = argparse.ArgumentParser(description='This script parses information from the default output of blastp')

ap.add_argument('-b', type=str, nargs=1, help='blastp file (required)')
ap.add_argument('-o', type=str, nargs=1, help='Output file name (required)')
ap.add_argument('-f', type=str, nargs=1, help='Format of input DB ID (default is UniRef, use "ensembl" for Ensembl generated DB)')

opts = ap.parse_args()

#check for missing args
missing_arg_flag = 0

if not opts.b:
    print("blastp file missing")
    missing_arg_flag = 1
if not opts.o:
    print("output name missing")
    missing_arg_flag = 1
if not opts.f:
    format_type = "uniref"
else:
    format_type = opts.f[0]


if missing_arg_flag == 1:
    print("Please try again with complete arguments")

blastp_file = opts.b[0]
outfile_name = opts.o[0]


print("opening blastp file")
#blastp_file = sys.argv[1]
#blastp_file_contents = open(blastp_file).read().rstrip("\n").split("\n")

#outfile_name = sys.argv[2]
outfile = open(outfile_name,"w")




class Hit:
    def __init__(self,s_name,q_name, s_len, q_len ,align_len,ident_len, ident_percent,q_start, q_end, s_start, s_end, e_val, pos_numerator, pos_denominator ):
        self.s_name = s_name
        self.q_name = q_name
        self.s_len = int(s_len)
        self.q_len = int(q_len)
        self.align_len = int(align_len)
        self.ident_len = int(ident_len)
        self.q_start = int(q_start)
        self.q_end = int(q_end)
        self.s_start = int(s_start)
        self.s_end = int(s_end)
        self.s_align_len = self.s_end + 1 - self.s_start
        self.q_align_len = self.q_end + 1 - self.q_start
        
        #self.ident_numerator = ident_numerator #same as ident_len
        #self.ident_denominator = ident_denominator # same as align_len
        self.pos_numerator = int(pos_numerator)
        self.pos_denominator = int(pos_denominator)
        
        #if self.s_align_len != self.q_align_len:
        #    print("Mismatching alignment lengths")
        #    print(q_name + "\t" + str(self.s_align_len))
        #    print(s_name + "\t" + str(self.q_align_len))
        #    sys.exit()
        
        self.ident_percent = ident_percent
        self.e_val = e_val
        

q_name = ""
q_len = 0
s_len = 0
align_start_flag = 0

q_start = 0
q_end = 0
s_start = 0
s_end = 0

ident_len = 0
align_len = 0
ident_percent = 0
#ident_numerator = 0 # same as ident_len
#ident_denominator = 0 # same as align_len
pos_numerator = 0
pos_denominator = 0

hit_start_flag = 0
empty_line_count = 0

s_name_flag = 0
s_name_list = []

query_list = [] # list of dicts for query hits
query_dict = {}

passed_headers_flag = 0
split_header_flag = 0

q_name = ""
s_name = ""

prev_line_length = 0

trans_id_list = []
trans_id_dict = {}

query_id_line_flag = 0

print("opening blastp file")
with open(blastp_file) as blastp_file_contents:
    for line in blastp_file_contents:

        line = line.rstrip("\n")

        # check that we are not in the middle of getting the subject name
        if s_name_flag == 1:

            # >E1C721 Uncharacterized protein OS=Gallus gallus GN=LOC425783 PE=4 SV=2

            # > sp|Q9HC56|PCDH9_HUMAN Protocadherin-9 OS=Homo sapiens OX=9606

            # >ENSGALP00000050219.1 pep chromosome:GRCg6a:MT:13071:14888:1 gene:ENSGALG00000029500.1
            # transcript:ENSGALT00000047001.1 gene_biotype:protein_coding
            # transcript_biotype:protein_coding gene_symbol:ND5 description:NADH-ubiquinone
            # oxidoreductase chain 5  [Source:UniProtKB/Swiss-Prot;Acc:P18940]
            # Length=605


            # length of query is taken care of after this
            if line.startswith("Length="):

                s_name = "".join(s_name_list)

                if format_type == "uniref":
                    s_name = s_name.split()[0]

                elif format_type == "ensembl":
                    this_gene_id = ""
                    this_trans_id = ""
                    this_desc_id = ""
                    this_transcript_biotype = ""

                    s_name_split = s_name.split()

                    for s_name_field in s_name_split:
                        if s_name_field.startswith("gene:"):
                            this_gene_id = s_name_field.split(":")[1]
                        if s_name_field.startswith("transcript:"):
                            this_trans_id = s_name_field.split(":")[1]

                        if s_name_field.startswith("description:"):
                            this_desc_id = s_name_field.split(":")[1]

                        if s_name_field.startswith("transcript_biotype:"):
                            this_transcript_biotype = s_name_field.split(":")[1]

                    s_name = ",".join([this_gene_id,this_trans_id])





                s_name_flag = 0
                s_name_list = []

                if s_len == 0:
                    s_len = line.split("=")[1]

                continue

            else:
                s_name_list.append(line)



        #because the query name wraps around to new lines this is needed to collect full names
        if query_id_line_flag > 0:
            if line.startswith("Length="):
                if q_len == 0:
                    q_len = line.split("=")[1]

                query_list.append(q_name)
                query_dict[q_name] = []

                # Bedtools 2.26 or earlier
                #>G1;G1.1::1:219-3261(+)
                #>G1;G1.2::1:1713-3246(+)

                # bedtools 2.27 or later
                #G1;G1.1(-)

                trans_id = q_name.split(":")[0]
                trans_id = trans_id.split("(")[0]
                if trans_id not in trans_id_dict:
                    trans_id_dict[trans_id] = {}
                    trans_id_list.append(trans_id)

                trans_id_dict[trans_id][q_name] = []

                query_id_line_flag = 0
            elif len(line) > 0:
                q_name = q_name + line

            continue


        #indicates that headers are over and can start parsing
        if passed_headers_flag > 0:

            #print(trans_id + "\t" + "passed_headers_flag" + "\t" + str(empty_line_count) + "\t" + str(hit_start_flag)  + "\t" + line)
            if empty_line_count > 1 and hit_start_flag == 1:

                hit_start_flag = 0
                empty_line_count = 0


                trans_id_dict[trans_id][q_name].append(Hit(s_name,q_name,s_len,q_len,align_len,ident_len,ident_percent,q_start,q_end,s_start,s_end,e_val,pos_numerator,pos_denominator))

                query_dict[q_name].append(Hit(s_name,q_name,s_len,q_len,align_len,ident_len,ident_percent,q_start,q_end,s_start,s_end,e_val,pos_numerator,pos_denominator))
                check_line = ",".join([s_name,q_name,str(s_len),str(q_len),str(align_len),str(ident_len),str(ident_percent),str(q_start),str(q_end),str(s_start),str(s_end),str(e_val),pos_numerator,pos_denominator])
                #print(check_line)
                s_name = ""
                s_len = 0
                align_len = 0
                ident_len = 0
                ident_percent = 0
                q_start = 0
                q_end = 0
                s_start = 0
                s_end = 0
                e_val = 0
                pos_numerator = 0
                pos_denominator = 0

            if line == "" and hit_start_flag == 1:
                #print("prev_line_length: " + str(prev_line_length))
                if prev_line_length == 0:
                    empty_line_count += 1
                else:
                    empty_line_count = 1

        prev_line_length = len(line)
        if prev_line_length > 0:
            empty_line_count = 0


        if split_header_flag > 0:
            q_name = query_line + line
            query_id_line_flag = 1
            passed_headers_flag = 1
            split_header_flag = 0

            #print(q_name)

            if "missing_nucleotides" in q_name:
                q_name = query_line + line
                query_id_line_flag = 0
                passed_headers_flag = 0
                split_header_flag = 0

                continue

        if not line.startswith((" S"," I","Query","Sbjct","Length",">")):
            continue


        if line.startswith("Query="):

            if "missing_nucleotides" in line:
                continue
            q_name = line.split()[1]
            #print(q_name)
            #sys.exit()
            q_len = 0
            hit_dict = {} #hit_dict
            passed_headers_flag = 1 # used to flag when we have passed all the header info
            query_id_line_flag = 1

            if line.endswith("-"):
                query_id_line_flag = 0
                passed_headers_flag = 0
                split_header_flag = 1
                query_line = q_name

                #print(query_line)

            continue



        if line.startswith(">"):
            # >E1C721 Uncharacterized protein OS=Gallus gallus GN=LOC425783 PE=4 SV=2

            # > sp|Q9HC56|PCDH9_HUMAN Protocadherin-9 OS=Homo sapiens OX=9606

            # >ENSGALP00000050219.1 pep chromosome:GRCg6a:MT:13071:14888:1 gene:ENSGALG00000029500.1
            # transcript:ENSGALT00000047001.1 gene_biotype:protein_coding
            # transcript_biotype:protein_coding gene_symbol:ND5 description:NADH-ubiquinone
            # oxidoreductase chain 5  [Source:UniProtKB/Swiss-Prot;Acc:P18940]
            # Length=605

            # use this to add onto to subject name in case of text wrap around in file
            s_name_list = []


            s_name_line = line[1:] # remove > character from beginning
            s_name_line = s_name_line.lstrip()
            # s_name_line = s_name_line.split()[0]
            s_len = 0
            hit_start_flag = 1
            empty_line_count = 0

            s_name_list.append(s_name_line)

            s_name_flag = 1

            continue

        # moved this statement to beginning of loop
        # length of query is take care of in first statement of loop
    #    if line.startswith("Length="):
    #        if s_len == 0:
    #            s_len = line.split("=")[1]
    #
    #        continue

        if line.startswith(" Score"):
            e_val = line.split(",")[1].split("= ")[1]
            continue

        if line.startswith(" Identities"):
            ident_line_split = line.split(",")
            ident_line = line.split(",")[0].split("= ")[1]
            ident_len = ident_line.split("/")[0]
            align_len = ident_line.split("/")[1].split()[0]
            ident_percent = ident_line.split("/")[1].split()[1].split("(")[1].split("%")[0]

            #ident_numerator = ident_line_split[0].split()[2].split("/")[0]
            #ident_denominator = ident_line_split[0].split()[2].split("/")[1]

            pos_numerator = ident_line_split[1].split()[2].split("/")[0]
            pos_denominator = ident_line_split[1].split()[2].split("/")[1]

            align_start_flag = 1
            continue


        if line.startswith("Query "):

            if len(line.split()) <  4:
                continue
            elif align_start_flag == 1:
                q_start = line.split()[1]
                q_end = line.split()[3]
            elif align_start_flag == 0:
                q_end = line.split()[3]
            continue

        if line.startswith("Sbjct "):
            if align_start_flag == 1:
                s_start = line.split()[1]
                s_end = line.split()[3]
                align_start_flag = 0
            elif align_start_flag == 0:
                s_end = line.split()[3]
            continue

        if line == "***** No hits found *****":
            s_name = ""
            s_len = 0
            align_len = 0
            ident_len = 0
            ident_percent = 0
            q_start = 0
            q_end = 0
            s_start = 0
            s_end = 0
            e_val = 0
            pos_numerator = 0
            pos_denominator = 0
            continue


# check protein alignments for each query and define CDS based on that
for trans_id in trans_id_list:
    
    best_match_line = ""
    best_align_percent = 0.0
    
    for q_name in trans_id_dict[trans_id]:
        q_name_split = q_name.split(":")
        trans_id = q_name_split[0]
        trans_id = trans_id.split("(")[0] #account for newer bedtools 2.27 formatting
        
        q_frame = q_name_split[-7]
        
        #for missing nucleotides the protein will not be complete so cannot do ORF prediction
        if q_frame == "missing_nucleotides":
            best_match_line = "\t".join([trans_id,"no_frame","-1","-1","-1","-1","none","missing_nucleotides","-1","-1"])
            break

        #>G1;G1.4::1:2350-3116(+):F2:1:74:74:E


        #print(q_name)
        q_rel_start = int(q_name_split[-4])
        q_rel_end = int(q_name_split[-3])

        q_nuc_start = int(q_name_split[-6])
        q_nuc_end = int(q_name_split[-5])
        
        #print(q_name)

        for s_obj in query_dict[q_name]:
            #print(s_obj.s_name +  "\t" + str(s_obj.s_align_len) + "\t" + str(s_obj.s_len))
            s_align_percent = float(s_obj.s_align_len) / float(s_obj.s_len)
            
            ident_percent = str(s_obj.ident_percent)
            
            s_align_percent_format = str(int(s_align_percent*100))
            # if there is a complete protein match use protein sequence for cds
            if int(s_align_percent) == 1:
                #new_q_rel_start = s_obj.q_start - 1 + q_rel_start # minus 1 to adjust for coordinates
                #new_q_rel_end = s_obj.q_end - 1 + q_rel_start

                new_q_rel_start = q_rel_start
                new_q_rel_end = q_rel_end

                if s_align_percent > best_align_percent:

                    match_flag = "full_match"
                    #best_match_line = "\t".join([trans_id,q_frame,str(q_rel_start),str(q_rel_end),str(new_q_rel_start),str(new_q_rel_end),s_obj.s_name,match_flag,s_align_percent_format,ident_percent])
                    best_match_line = "\t".join([trans_id, q_frame, str(q_nuc_start), str(q_nuc_end), str(new_q_rel_start), str(new_q_rel_end),s_obj.s_name, match_flag, s_align_percent_format, ident_percent])

                    best_align_percent = s_align_percent

            
            ########################################## Start up here. Need to find best match then use that as the output!!!!!
    
            #If there is 90% match then use the start site of the protein and the end position of the ORF        
            if s_align_percent >= 0.9:
                #new_q_rel_start = s_obj.q_start - 1 + q_rel_start # minus 1 to adjust for coordinates
                #new_q_rel_end = q_rel_end

                new_q_rel_start = q_rel_start
                new_q_rel_end = q_rel_end

                
                if s_align_percent > best_align_percent:
                    match_flag = "90_match"
                    #best_match_line = "\t".join([trans_id,q_frame,str(q_rel_start),str(q_rel_end),str(new_q_rel_start),str(new_q_rel_end),s_obj.s_name,match_flag,s_align_percent_format,ident_percent])
                    best_match_line = "\t".join([trans_id, q_frame, str(q_nuc_start), str(q_nuc_end), str(new_q_rel_start), str(new_q_rel_end),s_obj.s_name, match_flag, s_align_percent_format, ident_percent])
                    best_align_percent = s_align_percent
     
            
            #If there is 50% match then use the start site of the protein match on the query and the end position of the ORF
            #If there is 50% match then use ORF start and end 
            if s_align_percent >= 0.5:
                #new_q_rel_start = s_obj.q_start - 1 + q_rel_start # minus 1 to adjust for coordinates
                new_q_rel_start = q_rel_start 
                new_q_rel_end = q_rel_end
                
                if s_align_percent > best_align_percent:
                    match_flag = "50_match"
                    #best_match_line = "\t".join([trans_id,q_frame,str(q_rel_start),str(q_rel_end),str(new_q_rel_start),str(new_q_rel_end),s_obj.s_name,match_flag,s_align_percent_format,ident_percent])
                    best_match_line = "\t".join([trans_id, q_frame, str(q_nuc_start), str(q_nuc_end), str(new_q_rel_start), str(new_q_rel_end),s_obj.s_name, match_flag, s_align_percent_format, ident_percent])
                    best_align_percent = s_align_percent
    
            
            #If there is less than 50% match use ORF start and end    
            if s_align_percent < 0.5:
                new_q_rel_start = q_rel_start 
                new_q_rel_end = q_rel_end
                
                if s_align_percent > best_align_percent:
                    match_flag = "bad_match"
                    #best_match_line = "\t".join([trans_id,q_frame,str(q_rel_start),str(q_rel_end),str(new_q_rel_start),str(new_q_rel_end),s_obj.s_name,match_flag,s_align_percent_format,ident_percent])
                    best_match_line = "\t".join([trans_id, q_frame, str(q_nuc_start), str(q_nuc_end), str(new_q_rel_start), str(new_q_rel_end),s_obj.s_name, match_flag, s_align_percent_format, ident_percent])
                    best_align_percent = s_align_percent
                    
            #break here so that we only use the top hit
            #break
    
    # when no hits are found 
    if len(best_match_line) == 0:
        max_len = 0
        max_query = ""
        for q_name in trans_id_dict[trans_id]:
            q_name_split = q_name.split(":")
            trans_id = q_name_split[0]
            trans_id = trans_id.split("(")[0] #account for newer bedtools 2.27 formatting

            q_frame = q_name_split[-7]

            q_rel_start = int(q_name_split[-4])
            q_rel_end = int(q_name_split[-3])

            q_nuc_start = int(q_name_split[-6])
            q_nuc_end = int(q_name_split[-5])

            q_len = int(q_name_split[-2])
            if q_len > max_len:
                max_len = q_len
                max_query = q_name
                max_frame = q_frame
                max_q_rel_start = q_rel_start
                max_q_rel_end = q_rel_end

                max_q_nuc_start = q_nuc_start
                max_q_nuc_end = q_nuc_end
            
        best_match_line = "\t".join([trans_id,max_frame,str(max_q_nuc_start),str(max_q_nuc_end),str(max_q_rel_start),str(max_q_rel_end),"none","no_hit","0","0"])
        
        
    outfile.write(best_match_line)
    outfile.write("\n")
        




   