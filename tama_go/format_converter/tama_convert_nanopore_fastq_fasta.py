
import re
import sys
import time

#
# This script comverts nanopore fastq to fasta
#

print("opening fastq")
fastq_file = sys.argv[1]
fastq_file_contents = open(fastq_file).read().rstrip("\n").split("\n")


outfile_name = sys.argv[2]
outfile = open(outfile_name,"w")



#@8d5574f8-f8bd-42a8-a7cd-48c2d861d28f runid=4bc3cf4ba4a58c18a40a60b88bd63beebc1d69bf read=1009 ch=111 start_time=2017-11-15T14:17:18Z
#TTGGTATTAGCGTTCAGATTTGGCAGCCGCGCCTGTCGCTCATCTTCTTTTTAGGTTCTTAAATCTGAAACTACGTGAGCCAGACCAAGGTCAATTCCTGAACCCTTATCCGGCCCGGCCTTGCTGTGTACGTCTCGGGTCTAGAACCAAAAAGGAAAAACCTGGCTGAAAAACCGAAATCACAAAACAAAAAGCATAAAAAAGAAAAATTCTTTCAATTACTTACTTTATTAACTTATTTATTATATAAATATATATTCACCTATCACCGTCTCCTGCTCATAATATAGCTGGAGCTATGCTGATGTGCGTCGGGCGTGTAAACCTCCAGAGGTGGGCTGATGGCCTGCTGAGGCGGTTTTGTTTCCCCTTTTTTTGAGTTTAATTCTGATTGATTTTCCTCTGGGGACGGATAATAAAACATGTAATATTTTTATAAGTGAAACGCA
#+
#(',-*%%%%')'+,(&()(*())&%#$$%'(),*/30((*,+)'))%0242/+*..7+3+01))*+*)'$$&%)')(&&*21214&'$$#&')*)(,,%&052/%(''-50'-/63-51./0')+)+''#%%&&(*$-).''('+-%$0662+)),/32),'%%##$$')+-+($"-/2,*+--5:,/,8984$*(+.8940+&,2/*(*)$&)*))'*&%'+(),282,0+/+'))(-6505<997189;35-+*(*&&%'%'&(&#&(((()()*'))+)+.'+((%$&$$&(%'%$&%(,--*(%%&-.6/3.-((&'''%*)(-()+%*(*/-)*),,/*')&%%(()''%./.774)#()&*/0*)::642*)'/093+'/.223425-'+(&&'&%)),3:/,,)3222811:82./.-)*()&/576((*(&%%''*('$$%

count_lines = 0

for line in fastq_file_contents:

    count_lines += 1

    # grab fasta header
    if line.startswith("@") and count_lines == 1:
        fasta_header = ">" + line[1:]

    elif count_lines == 2:
        fasta_seq = line

    if count_lines == 4:
        count_lines = 0

        outfile.write(fasta_header)
        outfile.write("\n")
        outfile.write(fasta_seq)
        outfile.write("\n")

        fasta_header = "na"
        fasta_seq = "na"



