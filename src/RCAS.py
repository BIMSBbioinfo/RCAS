#!/usr/bin/python2

def get_argument_parser():
    parser = argparse.ArgumentParser(
        description="RCAS provides intuitive reports and publication ready graphics"
        " from input peak intervals in BED foramt, generated from clip-seq.")

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
                        help="Annotation reference in gff3 format.")

    parser.add_argument("--run_motif", "-m",
                        default="False",
                        choices=["True", "False"],
                        help="True: run motif search."
                        " False (default): not run motif search.")

    parser.add_argument("--run_PATHrich", "-p",
                        default="False",
                        choices=["True", "False"],
                        help="True: run pathway enrichment."
                        " False (default): not run pathway enrichment.")

    parser.add_argument("--run_GOrich", "-t",
                        default="False",
                        choices=["True", "False"],
                        help="True: run GO-term enrichment."
                        " False (default): not run GO-term enrichment.")

    parser.add_argument("--run_coverage", "-c",
                        default="False",
                        choices=["True", "False"],
                        help="True: run coverage profile."
                        " False (default): not run coverage_profile profile.")

    return parser

def generate_config():
    pass

def call_snakemake():
    pass

if __name__ == '__main__':

    import argparse

    #process commandline Arguments
    parser = get_argument_parser()
    args = parser.parse_args()
