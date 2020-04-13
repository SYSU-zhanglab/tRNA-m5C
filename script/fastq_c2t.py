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

with open(file,'r') as input:
	line = input.readline()
	i = 0
	while (line):
		if i == 0:
			print line.strip("\n")
			i += 1
		elif i == 1:
			print line.strip("\n").replace("C","T")
			i += 1
		elif i == 2:
			print line.strip("\n")
			i +=1
		elif i == 3:
			print line.strip("\n")
			i = 0
		line = input.readline()