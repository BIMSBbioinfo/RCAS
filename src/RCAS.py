#!/usr/bin/python2

def get_argument_parser():
    parser = argparse.ArgumentParser(
        description="RCAS provides intuitive reports and publication ready graphics"
        " from input peak intervals in BED foramt,"
        " which are detected in clip-seq data."
        )

    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s 0.1')

    parser.add_argument("BED",
                        metavar="BED",
                        nargs="*",
                        help="Target intervals in BED format.")

    parser.add_argument("--RCAS_path", "-r",
                        metavar="path/to/RCAS",
                        required=True,
                        help="Path to RCAS.")

    parser.add_argument("--genome", "-g",
                        metavar="FILE",
                        required=True,
                        help="Reference genome whose version should conform with"
                        " generation of input peak invervals.")

    parser.add_argument("--gff3", "-f",
                        metavar="FILE",
                        required=True,
                        help="Annotation reference in gff3 format"
                        " whose version should conform with"
                        " genome version.")

    parser.add_argument("--run_motif", "-m",
                        default="False",
                        choices=["True", "False"],
                        help="True: run motif search."
                        " False (default): not run.")

    parser.add_argument("--run_PATHrich", "-p",
                        default="False",
                        choices=["True", "False"],
                        help="True: run pathway enrichment."
                        " False (default): not run.")

    parser.add_argument("--run_GOrich", "-t",
                        default="False",
                        choices=["True", "False"],
                        help="True: run GO-term enrichment."
                        " False (default): not run.")

    parser.add_argument("--run_coverage", "-c",
                        default="False",
                        choices=["True", "False"],
                        help="True: run coverage profile."
                        " False (default): not run.")

    parser.add_argument("--species", "-s",
                        default="human",
                        choices=["human", "fly", "worm", "mouse"],
                        help="Required for running GO-term or pathway enrichment."
                        " Default: human")

    parser.add_argument("--forcerun", "-F",
                        action='store_true',
                        help=("Force the re-execution of snakemake."
                              " Use this option if you want to have all "
                              "output in your workflow updated."))

    return parser

def extract_key(filename):
    base = os.path.basename(filename)
    key = base.split(".")[:-1]
    key = ".".join(key)
    return key

def generate_config(args):
    BED_files = args.BED
    infiles = {}

    for BED in BED_files:
        infiles[extract_key(BED)] = BED

    config = {
      "RCAS_path": os.path.abspath(args.RCAS_path),

      "gff3": args.gff3,

      "genome": args.genome,

      "species": args.species,

      "infile": infiles,

      "switch": {
        "run_motif": args.run_motif,
        "run_PATHrich": args.run_PATHrich,
        "run_GOrich": args.run_GOrich,
        "run_coverage": args.run_coverage
      }
    }

    # Writing JSON data
    with open('config.json', 'w') as f:
         json.dump(config, f, sort_keys=True, indent=4)

    print "wrote config.json.\n"

def call_snakemake(RCAS_path):
    print "start snakemake:\n"

    if args.forcerun:
        forcerun = "F"
    else:
        forcerun = ""

    subprocess.call("snakemake -%sp -s %s/src/RCAS.snakefile" % (forcerun, RCAS_path), shell=True)

if __name__ == '__main__':
    import argparse
    import json
    import os
    import subprocess

    #process commandline Arguments
    parser = get_argument_parser()
    args = parser.parse_args()

    #dump argument to config.json
    generate_config(args)

    #run snakemake
    call_snakemake(args.RCAS_path)
