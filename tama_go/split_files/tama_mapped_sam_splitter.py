import re
import sys
import time


#
# This script splits mapped sam files by chromosome
#



print("opening sam file")
sam_file = sys.argv[1]
sam_file_contents = open(sam_file).read().rstrip("\n").split("\n")

num_files = int(sys.argv[2])

outfile_prefix = sys.argv[3]



if sam_file.endswith("sam"):
    bam_flag = "SAM"
elif sam_file.endswith("bam"):
    bam_flag = "BAM"
else:
    print("Please use Sam or Bam input.")
    sys.exit()


print("going through sam/bam file")

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

chrom_read_dict = {} # chrom_read_dict[chrom] = list of read line
chrom_numread_dict = {} # chrom_numread_dict[chrom] = number of reads

chrom_list = []

sam_header_list = []

sam_count = 0

for line in sam_file_contents:


    line_split = line.split("\t")

#    if line_split[0] == "@SQ":
#        seq_name = line_split[1].split(":")[1]
#        seq_length = line_split[2].split(":")[1]
#
#        sam_scaffold_dict[seq_name] = seq_length
#        sam_scaffold_list.append(seq_name)
#
    if line.startswith("@"):
        sam_header_list.append(line)
        continue

    if line == "":
        continue

    sam_count += 1
    if sam_count % 10000 == 0:
        print("sam count " + str(sam_count))

    read_id = line_split[0]
    sam_flag = int(line_split[1])
    scaff_name = line_split[2]
    start_pos = int(line_split[3])

    cigar = line_split[5]
    read_seq = line_split[9]
    seq_list = list(read_seq)
    # mapped_flag = sam_flag_dict[sam_flag]

    if scaff_name not in chrom_read_dict:
        chrom_read_dict[scaff_name] = []
        chrom_numread_dict[scaff_name] = 0
        chrom_list.append(scaff_name)

    chrom_read_dict[scaff_name].append(line)
    chrom_numread_dict[scaff_name] += 1

# calculate split of files

# num_files = 10

num_file_reads = sam_count / num_files

prev_file_read_count = 0
file_read_count = 0
total_read_count = 0
file_count = 1
chrom_count = 0
prev_chrom_count = 0

last_file_chrom = 0

# outfile_name = outfile_prefix + "_" + str(prev_chrom_count) + "_" + str(chrom_count) + ".sam"
outfile_name = outfile_prefix + "_" + str(file_count)  + ".sam"
outfile = open(outfile_name, "w")

for header_line in sam_header_list:
    outfile.write(header_line)
    outfile.write("\n")

while total_read_count <  sam_count:
    chrom = chrom_list[chrom_count]
    chrom_numreads = chrom_numread_dict[chrom]

    prev_file_read_count = file_read_count
    file_read_count = file_read_count + chrom_numreads

    chrom_count += 1

    if file_read_count > num_file_reads and prev_file_read_count > 0:


        file_read_count = chrom_numreads
        file_count += 1
        outfile.close()
        outfile_name = outfile_prefix + "_" + str(file_count) + ".sam"
        outfile = open(outfile_name, "w")

        last_file_chrom = chrom_count
        prev_chrom_count = chrom_count

        for header_line in sam_header_list:
            outfile.write(header_line)
            outfile.write("\n")

    for read_line in chrom_read_dict[chrom]:
        outfile.write(read_line)
        outfile.write("\n")

        total_read_count += 1


outfile.close()










