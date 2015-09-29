#! /usr/bin/env python

def process_ref_line(line):
	# extract last column where the attributes are
	
	attribute = line[:-1].split("\t")[-1]
	attribute = attribute.split(";")
	
	return attribute

def update_gene_set(attribute, gene_set):
	#determine the association between gene_id  and (gene_type, gene_name)
	#update in gene_set: {gene_id:[gene_type, gene_name]}
	
	gene_id, gene_type, gene_name = [item.split("=")[-1]
									 for item in attribute
									 if item.split("=")[0] in
									 ("gene_id", "gene_type",  "gene_name")]
	if gene_id not in gene_set:
		gene_set[gene_id] = (gene_type, gene_name)

def process_cor_line(line):
	#process each line,
	#return coordinate info & feature ids 
	
	cor_info, feature_ids = line[:-1].split("|")

	return cor_info, feature_ids
	
def determine_feature(feature_ids, gene_set):
	# determine appropriate features based on the reference ids  
	
	# genecode use following rule to rgister different features
	# for gene and transcript: only use  ID
	# for other features such as exon, UTR, intron, etc. use feature:transcript_ID:number
	
	# if there are no overlap with known gene coordinartes
	# it is annotated as intergenic
	
	# if there are overlap with known gene coordinates
	# make step-wise determination follwoing order: UTR3, UTR5, exon
	# if the above features are not associated, then it is annoated as intron
	
	if feature_ids == "":
		gene_type, gene_name, feature, gene_id  = "unknown", "unknown",  "intergenic", "unknown"
		
	else:
		
		if "UTR5" in feature_ids:
			feature = "UTR5"
			
		elif "UTR3" in feature_ids:
			feature = "UTR3"
			
		elif "CDS" in feature_ids:
			feature = "CDS"

		elif "exon" in feature_ids:	
			feature = "exon"
		
		else:
			feature = "intron"
			
		for feature_id in feature_ids.split(";"):
			
			if gene_set.has_key(feature_id):
				gene_id = feature_id
				gene_type, gene_name =  gene_set[gene_id]
				break

	return "\t".join([gene_type, gene_name, feature, gene_id])	

def check_argv(argv):
	
	usage = """
	parse_anot.py
	
	version: 0.0.1
	authors: Dilmurat Yusuf
	
	Usage:
	
	$ parse_anot.py -r track < bedmapBED > outputBed
  
	$ bedmap_upstream_process ... | parse_anot.py -r track > outputBed
	
	The output consists of coordinate info, gene_type, gene_name, feature,  gene_id.
	The format is tab delimited bed format.

	Parse feature ids from bedmap output and then map to genomic featurs
	according to given track which is specified by '-r' or '--track'.
	
	The script assums each cooridinate is only assocaited with a single gene_id.
	
	The reported features include exon, UTR, CDS, intron, intergenic,
	gene_type, gene_name, gene_id.
	
	For a coordinate which is assocaited with a gene_id
	but NOT with exons, it is annotated as intron.

	For a coordiate which is NOT assocaited with a genes_ids,
	it is annotated as intergenic.
	
	"""

	try:
		opts, args = getopt.getopt(argv,"hr:",["track="])
		
	except getopt.GetoptError:
		print usage
		sys.exit(1)
		
	if opts == []:
		print usage
		sys.exit(1)	
	
	for opt, arg in opts:
		if opt == '-h' :
			print  usage
			sys.exit()
		elif opt in ("-r", "--track"):
			track = arg
			
			return track
		
		
if __name__ == '__main__':
	import sys, getopt
	
	"""
	the script consists of two components
	
	component 1: retrieve reference info
	
	update gene_set = {gene_id:[gene_type, gene_name]}
	"""
	
	track = check_argv(sys.argv[1:])
	
	# extract attribute columns
	with open(track) as read_track:
		attribute_set = [process_ref_line(line) for line in read_track]
	
	# gene_set: {gene_id:[gene_type, gene_name]}
	gene_set = {}
	
	# update gene set 
	[update_gene_set(attribute, gene_set) for attribute in attribute_set]
	
	
	#----------------------------------------------------------
	
	"""
	component 2: map ids to features
	"""
	
	read_cor_ids = sys.stdin	
	
	#extract coordinate info & feature ids
	cor_ID_set = [process_cor_line(line) for line in read_cor_ids]

	#determine appropriate features based on the reference ids
	# cor_id[0] :  coordinate info
	# cor_id[1] : feature ids
	cor_feature_set = ["\t".join([
					  cor_id[0], determine_feature(cor_id[1], gene_set)
					  ])
					  for cor_id in cor_ID_set]
	
	#transform list into strings
	print "\n".join(cor_feature_set)
	
	
		
	
