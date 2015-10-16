import glob, os

TRACK_gff = config["ref"]

infile = config["infile"]
infile = glob.glob(infile)
outfile = [out.split(".")[0] + ".anot.tsv" for out in infile]

rule target:
	 input: outfile

rule intersect:
	 #obtain intersect between infile and TRACK_gff
	 input: infile
	 output: "{sample}.intersect.bed"
	 shell: "bedtools intersect -b {TRACK_gff} -a {input[0]}  -wao > {output}"
	 
rule anot_cor:
	#annotate cors with features
	input: "{sample}.intersect.bed", TRACK_gff
	output: "{sample}.anot.tsv"
	shell: "parse_anot.py < {input[0]}  > {output}"
