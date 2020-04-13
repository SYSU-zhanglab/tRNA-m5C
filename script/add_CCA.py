from Bio import SeqIO
import sys

for seq in SeqIO.parse(sys.argv[1], "fasta"):
	print ">"+seq.id
	print str(seq.seq) + "CCA"
