import pysam
import sys
import argparse
from Bio.Seq import reverse_complement

def recovery_read(output,read):
	if read.is_read1 == True:
		if read.query_name in fwd_reads:
			seq = fwd_reads.get(read.query_name)
			# All reads were forward
			qual = read.query_qualities
			mpq = read.mapping_quality
			segment_output = read
			segment_output.query_sequence = seq
			segment_output.mapping_quality = mpq
			segment_output.query_qualities = qual
			output.write(segment_output)
		else:
			output.write(read)
	elif read.is_read2 == True: 
		if read.query_name in rev_reads:
			seq = rev_reads.get(read.query_name)
			# here read2 is 100% reverse
			seq = reverse_complement(seq)
			qual = read.query_qualities
			mpq = read.mapping_quality
			segment_output = read
			segment_output.query_sequence = seq
			segment_output.mapping_quality = mpq
			segment_output.query_qualities = qual
			output.write(segment_output)
		else:
			output.write(read)


parser = argparse.ArgumentParser(fromfile_prefix_chars='@',description=__doc__,formatter_class=argparse.RawTextHelpFormatter)
#Require
group_required = parser.add_argument_group("Required")
group_required.add_argument("-i","--input",dest="input",required=True,help="input sam")
group_required.add_argument("-f","--forward",dest="forward",required=True,help="input forward fastq")
group_required.add_argument("-r","--reverse",dest="reverse",required=True,help="input reverse fastq")
group_required.add_argument("-o","--output",dest="output",required=True,help="output bam")
options = parser.parse_args()

fwd_reads = {}
with open(options.forward,"r") as input:
	n = 0
	line = input.readline()
	while (line):
		if n == 0:
			name = line.strip()[1:].split(" ")[0]
			n += 1
			line = input.readline()
		elif n == 1:
			seq = line.strip()
			if "C" in seq:
				fwd_reads[name] = seq
			n += 1
			line = input.readline()
		elif n == 2:
			n += 1
			line = input.readline()
		elif n == 3:
			n = 0
			line = input.readline()
			
rev_reads = {}
with open(options.reverse,"r") as input:
	n = 0
	line = input.readline()
	while (line):
		if n == 0:
			name = line.strip()[1:].split(" ")[0]
			n += 1
			line = input.readline()
		elif n == 1:
			seq = line.strip()
			if "G" in seq:
				rev_reads[name] = seq
			n += 1
			line = input.readline()
		elif n == 2:
			n += 1
			line = input.readline()
		elif n == 3:
			n = 0
			line = input.readline()

if options.input.endswith("sam"):
	with pysam.AlignmentFile(options.input,"r") as input, pysam.AlignmentFile(options.output,"wb",template=input) as output:
		for read in input:
			recovery_read(output,read)

else:
	with pysam.AlignmentFile(options.input,"rb") as input, pysam.AlignmentFile(options.output,"wb",template=input) as output:
		for read in input:
			recovery_read(output,read)
