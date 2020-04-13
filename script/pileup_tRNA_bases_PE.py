from Bio import SeqIO
import pandas as pd
import argparse
import pysam

def report(row):
	if row["ref_base"] == "C" and row["C"] + row["T"] >= 10:
		return round(row["C"]/(row["C"]+row["T"]+0.0),5)
	else:
		return "NaN"

if __name__ == "__main__":
	description = """pileup_tRNA_bases"""
	parser = argparse.ArgumentParser(prog="pileup_tRNA_bases",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	#Require
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("-r","--ref",dest="reference",required=True,help="reference fasta")
	group_required.add_argument("-b","--bam",dest="bam",required=True,help="input bam, sorted")
	group_required.add_argument("-o","--output",dest="output",required=True,help="output")
	# group_optional = parser.add_argument_group("Optional")
	options = parser.parse_args()
	
	reference = options.reference
	BAM = options.bam
	
	# start = 0
	# chr,start=start,stop=end,
	
	reference_genome = {}
	output = {}
	
	for seq in SeqIO.parse(reference,"fasta"):
		reference_genome[seq.id] = str(seq.seq).upper()

	with pysam.AlignmentFile(BAM, 'rb') as samfile:
		for pileupcolumn in samfile.pileup(max_depth=999999,truncate=True, ignore_orphans=False, ignore_overlaps=False, flag_filter=1536): 
			ref_base = reference_genome[pileupcolumn.reference_name][pileupcolumn.pos]
			KEY = (pileupcolumn.reference_name,pileupcolumn.pos+1) # 1-based
			if KEY not in output:
				output[KEY] = {"ref_base":ref_base,"A":0,"T":0,"C":0,"G":0,"del":0}
			TEMP = {}
			for pileupread in pileupcolumn.pileups:
				query_name = pileupread.alignment.query_name
				
				if pileupread.is_refskip == 1: # insertion
					continue

				if pileupread.is_del == 1:
					if query_name not in TEMP:
						TEMP[query_name] = "del"
					else:
						if TEMP[query_name] == "del":
							continue
						else:
							TEMP[query_name] = "conflict"
				else:
					quality = pileupread.alignment.query_qualities[pileupread.query_position]
					if quality < 30:
						continue
					base = pileupread.alignment.query_sequence[pileupread.query_position]
					if query_name not in TEMP:
						TEMP[query_name] = base
					else:
						if TEMP[query_name] == base:
							continue
						elif TEMP[query_name] == "N":
							TEMP[query_name] = base
						else:
							TEMP[query_name] = "conflict"
			for value in TEMP.values():
				if (value != "conflict") and (value != "N"):
					output[KEY][value] += 1
					

	df = pd.DataFrame(output).T[["ref_base","A","T","C","G","del"]]
	# df["query_pos"] = df.index
	#	A	T	C	G	del
	#0	xx	xx	xx	xx	xx
	# print df
	df.index.names = ["tRNA", "pos_1"]
	df["level"] = df.apply(report,axis=1)
	df.to_csv(options.output,sep="\t")
	# df.apply(lambda x: cons_sequence(x), axis=1)

