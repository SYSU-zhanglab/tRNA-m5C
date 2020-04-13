import pysam
import argparse

def report(output, primary_alignments, all_alignments):
	temp = []
	best_alignments = []
	best_alignments_ref = []
	primary_ref = None
	for i in all_alignments:
		if i.reference_name.startswith("MT"):
			temp.append(i.reference_name)
		elif i.reference_name.startswith("ENST"):
			temp.append(i.reference_name)
		elif i.reference_name.startswith("tRNA"):
			temp.append("-".join(i.reference_name.split("-")[0:3])) # to isodecoder level
		else:
			temp.append(i.reference_name)
	if len(all_alignments) == len(primary_alignments):
		for i in primary_alignments:
			if i.reference_name.startswith("tRNA") or i.reference_name.startswith("MT"):
				output.write(i)
	else:
		primary_AS = 0
		for i in primary_alignments:
			primary_AS += i.get_tag("AS")
			if i.reference_name.startswith("MT"):
				primary_ref = i.reference_name
			elif i.reference_name.startswith("ENST"):
				primary_ref = i.reference_name
			elif i.reference_name.startswith("tRNA"):
				primary_ref = "-".join(i.reference_name.split("-")[0:3])
			else:
				primary_ref = i.reference_name
			#primary_ref = "-".join((i.reference_name.split("-")[0:3]))
			best_alignments.append(i)
			best_alignments_ref.append(primary_ref)			

		for i in all_alignments:
			if i.is_secondary is True:
				score = i.get_tag("AS") + i.get_tag("YS")
				if score == primary_AS:
					best_alignments.append(i)
					if i.reference_name.startswith("MT"):
						best_alignments_ref.append(i.reference_name)
					elif i.reference_name.startswith("ENST"):
						best_alignments_ref.append(i.reference_name)
					elif i.reference_name.startswith("tRNA"):
						best_alignments_ref.append("-".join(i.reference_name.split("-")[0:3]))
					else:
						best_alignments_ref.append(i.reference_name)
		names_zipped = list(set(best_alignments_ref)) 
		#print names_zipped
		if len(names_zipped) == 1:
			if names_zipped[0].startswith("tRNA") or names_zipped[0].startswith("MT"):
				for i in best_alignments:
					#print names_zipped
					output.write(i)

if __name__ == "__main__":
	description = """filter_tRNA_alignments"""
	parser = argparse.ArgumentParser(prog="filter_tRNA_alignments",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	#Require
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("-i","--input",dest="bam",required=True,help="input bam")
	group_required.add_argument("-o","--output",dest="output",required=True,help="output")
	options = parser.parse_args()
	
	with pysam.AlignmentFile(options.bam, 'rb') as input, pysam.AlignmentFile(options.output, 'wb', template=input) as output:
		lastName = None
		primary_alignments = []
		all_alignments = []
		for read in input:
			if lastName is None:
				lastName = read.query_name
				if read.is_secondary is False:
					primary_alignments.append(read)
				all_alignments.append(read)
			else:
				if lastName == read.query_name:
					if read.is_secondary is False:
						primary_alignments.append(read)
					all_alignments.append(read)
				else:
					report(output, primary_alignments, all_alignments)
					lastName = read.query_name
					primary_alignments = []
					all_alignments = []
					if read.is_secondary is False:
						primary_alignments.append(read)
					all_alignments.append(read)
		report(output, primary_alignments, all_alignments)
					
