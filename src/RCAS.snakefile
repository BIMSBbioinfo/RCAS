import glob, os

def report_arguments(run_motif, run_PATHrich, run_GOrich, run_coverage):
	 if run_motif:
		  meme_out = "{input[5]}/meme_out/"
		  motif_annot = "{input[6]}"
		  memechip_out = "{sample}_memechip_output"
		  anot_motif = "{sample}.anot-motif.tsv"
	 else:
		  meme_out = "NOT_RUN"
		  motif_annot = "NOT_RUN"
		  memechip_out = None
		  anot_motif = None

	 if run_PATHrich:
		  msigdb = "{input[4]}"
		  msigdb_results = "{sample}.msigdb.results.tsv"
	 else:
		  msigdb = "NOT_RUN"
		  msigdb_results = None

	 if run_GOrich:
		  go_bp = "{input[3]}/BP.GO.results.tsv"
		  go_mf = "{input[3]}/MF.GO.results.tsv"
		  go_cc = "{input[3]}/CC.GO.results.tsv"
		  GO_term = "{sample}-GO-term"
	 else:
		  go_bp = "NOT_RUN"
		  go_mf = "NOT_RUN"
		  go_cc = "NOT_RUN"
		  GO_term = None

	 if run_coverage:
		  coverage_profile = "RUN"
	 else:
		  coverage_profile = "NOT_RUN"

	 cmd = ["bash {RCAS_path}/src/generate_report.sh",
		  "--output_filename={output} --annot={input[0]} --peaks={input[1]} --gff3={input[2]}",
		  "--go_bp=%s" % go_bp,
		  "--go_mf=%s" % go_mf,
		  "--go_cc=%s" % go_cc,
		  "--msigdb=%s" % msigdb,
		  "--meme_out=%s --motif_annot=%s" % (meme_out, motif_annot),
		  "--coverage_profile_option=%s" % coverage_profile]

	 cmd = " ".join(cmd)
	 imput_args = ["{sample}.anot.tsv",
					infile,
					{TRACK_gff},
					GO_term,
					msigdb_results,
					memechip_out,
					anot_motif]

	 imput_args = [arg for arg in imput_args if arg != None]

	 return imput_args, cmd


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

run_motif = False
run_PATHrich = False
run_GOrich = False
run_coverage = True

imput_args, cmd = report_arguments(run_motif, run_PATHrich, run_GOrich, run_coverage)

rule target:
	 input:
		  outfile

#default step
include: anot

#optionl step
if run_motif:
	 include: motif

#optional step
if run_PATHrich:
	 include: PATHrich

#optional step
if run_GOrich:
	 include: GOrich

#default step
rule html_report:
	 input:
		  imput_args
	 output:
		  "{sample}.rcas.html"
	 shell:
		  cmd
