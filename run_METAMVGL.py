#!/usr/bin/env python
import argparse
import datetime
import os
import subprocess
import sys


def logging(message):
    print(f"[{datetime.datetime.now()}] {message}", flush=True)


def run_cmd(command, name, verbose=True):
    if verbose:
        logging(f'{name} starts: {" ".join(command)}')
    ret = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    if ret.returncode:
        logging(f'{name} fails: {" ".join(command)}')
        sys.exit(1)
    elif verbose:
        logging(f'{name} ends: {" ".join(command)}')


parser = argparse.ArgumentParser()
parser.add_argument(
    "-a",
    "--assembler",
    required=True,
    help="the assembler that was used to generate contigs (choose from metaSPAdes and MEGAHIT)",
)
parser.add_argument(
    "-b",
    "--initial_binning",
    required=False,
    default="MetaBat2",
    help="the initial binning tool to create initial binning results (choose from MetaBat2 and MaxBin2)",
)
parser.add_argument(
    "-r1",
    "--reads1",
    required=True,
    help="the path to reads1",
)
parser.add_argument(
    "-r2",
    "--reads2",
    required=True,
    help="the path to reads2",
)
parser.add_argument(
    "-c",
    "--contigs",
    required=True,
    help="the path to the contigs",
)
parser.add_argument(
    "-p",
    "--paths",
    required=False,
    help="the path to the contigs.paths file under metaSPAdes output folder, needed for metaSPAdes",
)
parser.add_argument(
    "-g",
    "--assembly_graph",
    required=True,
    help="the path to the assembly_graph.fastg file under metaSPAdes output folder (metaSPAdes), or fastg format file ceretad by megahit_toolkit (MEGAHIT)",
)
parser.add_argument(
    "-m",
    "--mapping_quality",
    required=False,
    default=10,
    help="the threshold of mapping quality for reads alignment",
)
parser.add_argument(
    "-s",
    "--identity",
    required=False,
    default=0.95,
    help="the threshold of alignment identity for reads alignment",
)
parser.add_argument(
    "-i",
    "--insert_size",
    required=False,
    default=270,
    help="the insert size of paired-end reads",
)
parser.add_argument(
    "-n",
    "--pe",
    required=False,
    default=3,
    help="the minimum number of paired-end reads to support a PE link",
)
parser.add_argument(
    "-u",
    "--max_iter",
    required=False,
    default=100,
    help="the maximun number of iteration of label propagation",
)
parser.add_argument(
    "-d",
    "--threshold",
    required=False,
    default=0.00000001,
    help="the threshold to stop iteration of label propagation",
)
parser.add_argument(
    "-t",
    "--threads",
    required=False,
    default=16,
    help="the number of threads for initial binning tools",
)
parser.add_argument(
    "-o",
    "--output",
    required=True,
    default=3,
    help="output dir",
)
args = parser.parse_args()

# check assembler
if args.assembler.lower() not in ("metaspades", "megahit"):
    logging("Please check the assembler name, which should be metaSPAdes or MEGAHIT.")
    sys.exit(1)
if args.assembler.lower() == "metaspades" and args.paths is None:
    logging("Please provide --paths or -p for metaSPAdes assembler.")
    sys.exit(1)

# check initial binning tool
if args.initial_binning.lower() not in ("maxbin2", "metabat2"):
    logging(
        "Please check the --initial_binning or -i, which should be MetaBat2 or MaxBin2."
    )
    sys.exit(1)

# check paths
if not os.path.isfile(args.reads1):
    logging("Please input a valid path for --reads1 or -r1.")
    sys.exit(1)
if not os.path.isfile(args.reads2):
    logging("Please input a valid path for --reads2 or -r2.")
    sys.exit(1)
if not os.path.isfile(args.contigs):
    logging("Please input a valid path for --contigs or -c.")
    sys.exit(1)
if not os.path.isfile(args.assembly_graph):
    logging("Please input a valid path for --assembly_graph or -g.")
    sys.exit(1)
if args.assembler.lower() == "metaspades" and not os.path.isfile(args.paths):
    logging("Please input a valid path for --paths or -p.")
    sys.exit(1)

# create output dir
if not os.path.isdir(args.output):
    run_cmd(["mkdir", args.output], "mkdir")
else:
    logging("The output dir already exists. Please choose another dir.")
    sys.exit(1)

# initial binning
bindir = ""
bamfile = ""
if args.initial_binning.lower() == "maxbin2":
    run_cmd(["mkdir", f"{args.output}/MaxBin2"], "mkdir")
    run_cmd(
        [
            "run_MaxBin.pl",
            "-contig",
            f"{args.contigs}",
            "-reads",
            args.reads1,
            "-reads2",
            args.reads2,
            "-out",
            f"{args.output}/MaxBin2/bin",
            "-thread",
            str(args.threads),
        ],
        "MaxBin2",
    )
    bindir = f"{args.output}/MaxBin2"
else:
    run_cmd(["mkdir", f"{args.output}/MetaBat2"], "mkdir")
    logging("Create bam file starts.")
    if not os.path.isfile(f"{args.contigs}.bwt"):
        subprocess.run(["bwa", "index", args.contigs])
    pipe = subprocess.Popen(
        ["bwa", "mem", "-t", str(args.threads), args.contigs, args.reads1, args.reads2],
        stdout=subprocess.PIPE,
    )
    subprocess.check_output(
        [
            "samtools",
            "sort",
            "-@",
            str(args.threads),
            "-",
            "-o",
            f"{args.output}/MetaBat2/align.megahit.bam",
        ],
        stdin=pipe.stdout,
    )
    pipe.communicate()
    logging("Create bam file ends.")
    run_cmd(
        [
            "jgi_summarize_bam_contig_depths",
            "--outputDepth",
            f"{args.output}/MetaBat2/align.megahit.depth",
            f"{args.output}/MetaBat2/align.megahit.bam",
        ],
        "Metabat2_bam_depth",
    )
    run_cmd(
        [
            "metabat2",
            "-i",
            args.contigs,
            "-a",
            f"{args.output}/MetaBat2/align.megahit.depth",
            "-m",
            str(1500),
            "-o",
            f"{args.output}/MetaBat2/bin",
            "-t",
            str(args.threads),
        ],
        "Metabat2_main",
    )
    bindir = f"{args.output}/MetaBat2"
    bamfile = f"{args.output}/MetaBat2/align.megahit.bam"
if args.assembler.lower() == "metaspades":
    run_cmd(
        [
            "prepResult.py",
            "--binned",
            bindir,
            "--output",
            bindir,
            "--assembler",
            "spades",
        ],
        "prepResult",
    )
else:
    run_cmd(
        [
            "prepResult.py",
            "--binned",
            bindir,
            "--output",
            bindir,
            "--assembler",
            "megahit",
        ],
        "prepResult",
    )
initial_results = f"{bindir}/initial_contig_bins.csv"

# prepare bam file
if bamfile == "":
    logging("Create bam file starts.")
    if not os.path.isfile(f"{args.contigs}.bwt"):
        subprocess.run(["bwa", "index", args.contigs])
    pipe = subprocess.Popen(
        ["bwa", "mem", "-t", str(args.threads), args.contigs, args.reads1, args.reads2],
        stdout=subprocess.PIPE,
    )
    subprocess.check_output(
        [
            "samtools",
            "sort",
            "-@",
            str(args.threads),
            "-o",
            f"{args.output}/align.contigs.bam",
        ],
        stdin=pipe.stdout,
    )
    pipe.communicate()
    logging("Create bam file ends.")
    bamfile = f"{args.output}/align.contigs.bam"

# prep_graph and METAMVGL
if args.assembler.lower() == "metaspades":
    run_cmd(
        [
            "prep_graph",
            "-a",
            "metaSPAdes",
            "-p",
            args.paths,
            "-g",
            args.assembly_graph,
            "-b",
            bamfile,
            "-o",
            f"{args.output}/graph",
        ],
        "prep_graph",
    )
    run_cmd(
        [
            "METAMVGL.py",
            "--assembler",
            "metaSPAdes",
            "--contigs",
            args.contigs,
            "--assembly_graph",
            f"{args.output}/graph.ag",
            "--PE_graph",
            f"{args.output}/graph.pe",
            "--binned",
            initial_results,
            "--output",
            f"{args.output}/METAMVGL",
        ],
        "METAMVGL",
    )
else:
    run_cmd(
        [
            "prep_graph",
            "-a",
            "MEGAHIT",
            "-c",
            args.contigs,
            "-g",
            args.assembly_graph,
            "-b",
            bamfile,
            "-o",
            f"{args.output}/graph",
        ],
        "prep_graph",
    )
    run_cmd(
        [
            "METAMVGL.py",
            "--assembler",
            "MEGAHIT",
            "--contigs",
            args.contigs,
            "--assembly_graph",
            f"{args.output}/graph.ag",
            "--PE_graph",
            f"{args.output}/graph.pe",
            "--binned",
            initial_results,
            "--output",
            f"{args.output}/METAMVGL",
        ],
        "METAMVGL",
    )
