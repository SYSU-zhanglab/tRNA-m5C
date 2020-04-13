from Bio import SeqIO
import sys
data  = {}
for seq in SeqIO.parse(sys.argv[1],"fasta"):
	# tRNA-Ala-AGC-1-1
	ID = seq.id.split("-")
	new_id = "-".join(ID[:-1])
	if new_id not in data:
		data[new_id] = str(seq.seq).upper()
	else:
		if data[new_id] != str(seq.seq).upper():
			data[seq.id] = str(seq.seq).upper()

for key, values in data.items():
	print ">"+key
	print values
