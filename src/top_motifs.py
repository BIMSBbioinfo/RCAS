#! /usr/bin/env python

def check_argv(argv):

	usage = """
	top_motifs.py

	version: 0.0.1
	author: Dilmurat Yusuf

	Usage:

	$ top_motifs.py -m centrimo.html -c cor_flank -a annotation_table -n number_of_top_motifs > annotation_motifs.tsv

	Default number of top motifs: 10

	"""

	try:
		opts, args = getopt.getopt(argv,"m:c:a:h:n:",[""])

	except getopt.GetoptError:
		print usage
		sys.exit(1)

	if set([opt for opt, arg in opts]) not in [set(["-m", "-a", "-c"]), set(["-m", "-a", "-c", '-n'])]:
		print usage
		sys.exit(1)

	number_of_top_motifs = 10

	for opt, arg in opts:
		if opt == '-h' :
			print  usage
			sys.exit()

		if opt == '-m':
			centrimo_html = arg

		if opt == "-c":
			coordinates_original_flank = arg

		if opt == '-a':
			annotation = arg

		if opt == '-n':
			number_of_top_motifs = int(arg)

	return 	centrimo_html, coordinates_original_flank, annotation, number_of_top_motifs

def extract_motif_data(line, on, motif_data):
	if "//@JSON_VAR data" in line:
		on = True

	if  line.endswith("};\n"):
		on = False

	if on:
		motif_data.append(line)

	return on, motif_data

def show_nucleotides(pwm):
		# order of PWM: A C G T

		def assess(site_pwm):
			nucleotide_composition = []
			if site_pwm[0] >= 0.1:
				nucleotide_composition.append("A")
			if site_pwm[1] >= 0.1:
				nucleotide_composition.append("C")
			if site_pwm[2] >= 0.1:
				nucleotide_composition.append("G")
			if site_pwm[3] >= 0.1:
				nucleotide_composition.append("T")
			return nucleotide_composition

		return [tuple(sorted(assess(site_pwm))) for site_pwm in pwm]

def pwm2IUPAC(pwm):

	IUPAC_consensus = {
		('A',): 'A',
		('C',): 'C',
		('G',): 'G',
		('T',): 'T',
		('A', 'G') : 'R',
		('C', 'T') : 'Y',
		('C', 'G') : 'S',
		('A', 'T') : 'W',
		('G', 'T') : 'K',
		('A', 'C') : 'M',
		('C', 'G', 'T') : 'B',
		('A', 'G', 'T') : 'D',
		('A', 'C', 'T') : 'H',
		('A', 'C', 'G') : 'V',
		('A', 'C', 'G', 'T') : 'N'
	}

	compositions = show_nucleotides(pwm)

	consensus = [IUPAC_consensus[site_composition] for site_composition in compositions]

	return "".join(consensus)

def update_mapped_coordinates(line, mapped_coordinates):
		line = line[:-1]
		line = line.split()

		seq_id = line[0]
		start_flank = line[1]
		end_flank = line[2]
		start_peak = line[-2]
		end_peak = line[-1]
		strand = line[5]

		coordinates_flank = "%s:%s-%s(%s)" % (seq_id, start_flank, end_flank, strand)
		coordinates_peak = "%s:%s-%s(%s)" % (seq_id, start_peak, end_peak, strand)

		mapped_coordinates[coordinates_flank] = coordinates_peak

		return mapped_coordinates

def update_anot_info(line, anot_info):
		line = line[:-1]
		line = line.split()

		seq_id = line[0]
		start_peak = line[1]
		end_peak = line[2]
		strand = line[5]

		coordinates_peak = "%s:%s-%s(%s)" % (seq_id, start_peak, end_peak, strand)
		annotation = line[8:]

		try:
			anot_info[coordinates_peak].append(annotation)
		except:
			anot_info[coordinates_peak] = [annotation]

		return anot_info

def update_table(motif, sequneces, mapped_coordinates, anot_info, table):
	matched_sites = "%.f" % motif[0]
	motif_id = motif[1]
	IUPAC_consensus = motif[2]
	matched_seq_ids = motif[3]

	for seq_id in matched_seq_ids:
		coordinates_flank = sequneces[seq_id]
		coordinates_peak = mapped_coordinates[coordinates_flank]

		anot_features = anot_info[coordinates_peak]

		chromosome_id, rest = coordinates_peak.split(":")
		coordinate, strand = rest[:-1].split("(")
		start, end = coordinate.split("-")

		for feature in anot_features:
			feature = "\t".join(feature)

			table.append("\t".join([chromosome_id, start, end, strand, feature, IUPAC_consensus, motif_id, matched_sites]))

	return table

if __name__ == '__main__':
	import sys, getopt, json

	argv = sys.argv[1:]

	#check commandline options
	centrimo_html, coordinates_peak_flank, annotation, number_of_top_motifs = check_argv(argv)

	###########process for motif info, start
	#extract motif data
	with open(centrimo_html) as handle:
		on = False
		motif_data = []

		for line in handle:
			on, motif_data = extract_motif_data(line, on, motif_data)

	motif_data = motif_data[2:]
	motif_data[0] = "{\n"
	motif_data.append("}\n")
	motif_data = "".join(motif_data)
	motif_data = json.loads(motif_data)
	#data stucture of motif_data is following:
	#[u'motif_dbs', u'cmd', u'seqlen', u'tested', u'program', u'motifs', u'sequences', u'release', u'sequence_db', u'options', u'revision']
	#u'motifs': [u'peaks', u'score_threshold', u'total_sites', u'db', u'sites', u'len', u'motif_nsites', u'n_tested', u'seqs', u'alt', u'pwm', u'motif_evalue', u'id']

	#there can be duplicated motifs
	#use set to remove duplications
	motif_info = set([(motif['peaks'][0]['sites'], motif['id'],  pwm2IUPAC(motif['pwm']), tuple(motif['seqs'])) for motif in motif_data['motifs']])
	motif_info = sorted(motif_info, reverse=True)


	sequneces = motif_data['sequences']
	###########process for motif info, end

	###########coordinate mapping, start
	#map original coordinates and 100-nt coordinates
	mapped_coordinates = {}

	with open(coordinates_peak_flank) as handle:
		for line in handle:
			mapped_coordinates = update_mapped_coordinates(line, mapped_coordinates)
	###########coordinate mapping, end

	###########extract annotation information, start
	#extract anot info: start, end, strand, feature, feature_id
	anot_info = {}

	with open(annotation) as handle:
		handle.next()

		for line in handle:
			anot_info = update_anot_info(line, anot_info)
	###########extract annotation information, end

	#generate table
	header = 'chromosome_id\tstart_position\tend_position\tstrand\tfeature\tfeature_id\tgene_id\ttranscript_id\tgene_type\tgene_name\tgene_strand\tmotif\tmotif_id\tmatched_sites'
	table = [header]

	for motif in motif_info[:number_of_top_motifs]:
		table = update_table(motif, sequneces, mapped_coordinates, anot_info, table)

	print "\n".join(table)

