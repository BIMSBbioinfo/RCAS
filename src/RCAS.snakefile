import glob, os

TRACK_gff = config["ref"]
TRACK_bed = ".".join(TRACK_gff.split(".")[:-1]) + ".bed"

infile = config["infile"]
outfile = [out.split(".")[0] + ".anot.bed" for out in glob.glob(infile)]


rule target:
	 input: outfile

rule clean:
     shell: "rm -f *sorted.bed *anot.bed"
	 
rule gff2bed:
	#convert gff to bed
	input: TRACK_gff
	output: TRACK_bed
	shell: 'gff2bed < {input} > {output}'
	 
rule sortbed:
	#sort bed file
	input: "{sample}.bed"
	output: "{sample}.sorted.bed"
	shell: 'sort-bed  {input} > {output}'

rule anot_cor:
	#annotate cors with features
	input: "{sample}.sorted.bed", TRACK_bed
	output: "{sample}.anot.bed"
	shell: "bedmap --echo --echo-map-id --fraction-ref .95  {input[0]} {TRACK_bed} | parse_anot.py -r {TRACK_bed}  > {output}"
