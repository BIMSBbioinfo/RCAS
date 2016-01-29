import glob, os

RCAS_path = config["RCAS_path"]
anot = RCAS_path  + "/src/RCAS.anot"
motif = RCAS_path  + "/src/RCAS.motif"
GOrich = RCAS_path  + "/src/RCAS.GOrich"
PATHrich = RCAS_path  + "/src/RCAS.PATHrich"

TRACK_gff = config["gff3"]
genome_reference = config["genome"]

infile = config["infile"]
infile = glob.glob(infile)
outfile = [out.split(".")[0] + ".rcas.html" for out in infile]

rule target:
	 input:
		  outfile

#default step
include: anot

#optionl step
include: motif

#optional step
include: PATHrich

#optional step
include: GOrich

#default step
rule html_report:
	 input:
		  "{sample}.anot.tsv",
		  infile,
		  {TRACK_gff},
		  "{sample}-GO-term",
		  "{sample}.msigdb.results.tsv",
		  "{sample}_memechip_output",
		  "{sample}.anot-motif.tsv"
	 output:
		  "{sample}.rcas.html"
	 shell:
		  "bash {RCAS_path}/src/generate_report.sh"
		  " --output_filename={output} --annot={input[0]} --peaks={input[1]} --gff3={input[2]}"
		  " --go_bp={input[3]}/BP.GO.results.tsv --go_mf={input[3]}/MF.GO.results.tsv"
		  " --go_cc={input[3]}/CC.GO.results.tsv --msigdb={input[4]}"
		  " --meme_out={input[5]}/meme_out/ --motif_annot={input[6]}"
