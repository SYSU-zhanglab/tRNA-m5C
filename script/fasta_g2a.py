#!/usr/bin/env python
import sys
import os
import getopt
from Bio import SeqIO

opts,args=getopt.getopt(sys.argv[1:],"i:h")
for op,value in opts:
	if op == "-i":
		file = value
	elif op == "-h":
		sys.exit()

for seq_record in SeqIO.parse(file,"fasta"):
	print ">"+seq_record.description
	print str(seq_record.seq).replace("G","A")