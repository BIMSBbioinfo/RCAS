import glob, os

RCAS_path = config["RCAS_path"]

TRACK_gff = config["gff3"]
genome_reference = config["genome"]

infile = config["infile"]
infile = glob.glob(infile)
outfile = [out.split(".")[0] + ".rcas.html" for out in infile]

rule target:
	 input:
		  outfile

rule intersect:
	 #obtain intersect between infile and TRACK_gff
	 input:
		  infile
	 output:
		  "{sample}.intersect.bed"
	 shell:
		  "bedtools intersect -b {TRACK_gff} -a {input[0]}  -wao > {output}"

rule anot_cor:
	 #annotate coordinates with features
	 input:
		  "{sample}.intersect.bed"
	 output:
		  "{sample}.anot.tsv"
	 shell:
		  "python2 {RCAS_path}/src/parse_anot.py < {input}  > {output}"

rule get_flanking_coordinates:
	 #prepare flanking coordinates centering on binding site
	 input:
		  "{sample}.bed"
	 output:
		  "{sample}-summit-100bp.bed"
	 shell:
		  "awk 'BEGIN{{ OFS=\"\t\";}}"
		  "{{ midPos=int(($2+$3)/2); start= midPos-50; end = midPos+50;"
		  "if (start <0) {{start = 1 ; end=100}};"
		  "print $1, start, end, $4, $5, $6, $2, $3;}}'  {input} > {output}"

rule get_fasta:
	 input:
		  "{sample}-summit-100bp.bed",
		  genome_reference
	 output:
		  "{sample}-summit-100bp.fa"
	 shell:
		  "fastaFromBed -s -fi {input[1]} -bed  {input[0]} -fo {output}"

rule run_meme_chip:
	 input:
		  "{sample}-summit-100bp.fa"
	 output:
		  "{sample}_memechip_output"
	 shell:
		  "meme-chip -meme-maxw 8 -norc -oc {output}"
		  " -db {RCAS_path}/src/base/Homo_sapiens-U2T.meme"
		  " -db {RCAS_path}/src/base/Mus_musculus-U2T.meme"
		  " -db {RCAS_path}/src/base/Drosophila_melanogaster-U2T.meme"
		  " -db {RCAS_path}/src/base/Caenorhabditis_elegans-U2T.meme {input}"

rule profile_top_motifs:
	 input:
		  "{sample}_memechip_output",
		  "{sample}-summit-100bp.bed",
		  "{sample}.anot.tsv"
	 output:
		  "{sample}.anot-motif.tsv"
	 shell:
		  "python2 {RCAS_path}/src/top_motifs.py"
		  " -m {input[0]}/centrimo_out/centrimo.html"
		  " -e {input[0]}/centrimo_out/centrimo.txt"
		  " -c {input[1]} -a {input[2]} > {output}"
	 
rule report_msigd:
	 input:
		  TRACK_gff,
		  "{sample}.anot.tsv"
	 output:
		  "{sample}.msigdb.results.tsv"
	 shell:
		  "Rscript {RCAS_path}/src/rcas.msigdb.R --gmt={RCAS_path}/src/base/c2.cp.v5.0.entrez.gmt"
		  " --gff3={input[0]} --anot={input[1]} --out={output}"

rule report_GO:
	 input:
		  TRACK_gff,
		  "{sample}.anot.tsv"
	 output:
		  "{sample}-GO-term"
	 shell:
		  "Rscript {RCAS_path}/src/rcas.GO.R --gff3={input[0]} --anot={input[1]} --out={output}"

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
		  " -o={output} -a={input[0]} -p={input[1]} -g={input[2]}"
		  " -b={input[3]}/BP.GO.results.tsv -m={input[3]}/MF.GO.results.tsv"
		  " -c={input[3]}/CC.GO.results.tsv -s={input[4]}"
		  " -e={input[5]}/meme_out/ -t={input[6]}"
