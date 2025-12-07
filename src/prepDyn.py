#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# auxiliary.py
# prepDyn - Copyright (C) 2025
# Daniel Y. M. Nakamura
# GNU General Public Licence version 3.0
# Contact: dani_ymn@outlook.com

COPYRIGHT="""
prepDyn - Data preprocessing for dynamic homology
Copyright (C) 2025 - Daniel Y. M. Nakamura
    prepDyn comes with ABSOLUTELY NO WARRANTY!
    This is a free software, and you are welcome to
    redistribute them under certain conditions
    For additional information, please visit:
    http://opensource.org/licenses/GPL-3.0
"""
print(COPYRIGHT)

###########
# MODULES #
###########

# The programs pip and mafft should be installed beforehand.
# The Python modules can be automatically installed using the
# following script. If already installed, it will load them.

import subprocess
import sys
import importlib

def install_and_import(package_name, import_path, from_list=None, alias=None):
    try:
        if from_list:
            module = importlib.import_module(import_path)
            for name in from_list:
                globals()[name] = getattr(module, name)
        else:
            globals()[alias or import_path] = importlib.import_module(import_path)
    except ModuleNotFoundError:
        print(f"Installing: {package_name}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package_name])
        # Try importing again after installation
        if from_list:
            module = importlib.import_module(import_path)
            for name in from_list:
                globals()[name] = getattr(module, name)
        else:
            globals()[alias or import_path] = importlib.import_module(import_path)

# Biopython modules
install_and_import("biopython", "Bio.AlignIO", from_list=["read", "write"])
install_and_import("biopython", "Bio.Entrez", from_list=["efetch", "email"])
install_and_import("biopython", "Bio.SeqIO", from_list=["parse", "write"])
install_and_import("biopython", "Bio.Align", from_list=["MultipleSeqAlignment"])
install_and_import("biopython", "Bio.Seq", from_list=["Seq"])
install_and_import("biopython", "Bio.SeqRecord", from_list=["SeqRecord"])

# Other useful packages
install_and_import("matplotlib", "matplotlib.pyplot", alias="plt")
install_and_import("numpy", "numpy", alias="np")
install_and_import("termcolor", "termcolor", from_list=["colored"])

# Standard libraries (should not need installation)
import argparse
import ast
import csv
from io import StringIO
import os
import pathlib
import re
import subprocess
import tempfile
import time

# prepDyn libraries
from prepDyn_auxiliary import prepDyn
from Bio import SeqIO
from Bio.Seq import Seq

###########
# PARSERS #
###########

def parse_internal_leaves(value):
    if value == "all":
        return "all"
    try:
        return value.split(",")
    except Exception:
        raise argparse.ArgumentTypeError("internal_leaves must be 'all' or a comma-separated list (e.g., seq1,seq2)")

def parse_internal_column_ranges(value):
    if value == "all":
        return "all"
    try:
        return ast.literal_eval(value)
    except Exception:
        raise argparse.ArgumentTypeError("internal_column_ranges must be 'all' or a valid Python list (e.g., [(10, 20), (30, 40)])")

def parse_partitioning_round(value):
    if value == "max":
        return value
    try:
        return int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("partitioning_round must be an integer or 'max'")

def parse_n2question_leaves(value):
    if value.lower() == "none":
        return None
    elif value.lower() == "all":
        return "all"
    else:
        return [leaf.strip() for leaf in value.split(",")]

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
    
########
# MAIN #
########

def main():
    parser = argparse.ArgumentParser(
        description="Preprocess sequences for dynamic homology in PhyG. The four steps are (1) data collection from GB, (2) trimming, (3) identification of missing data, and (4) successive partitioning.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""\
Examples:
  # Given a CSV file with GenBank accession numbers, download and align sequences, and identify terminal and internal missing data
  python prepDyn.py -gb accessions.csv -o output

  # Given a FASTA alignment, identify terminal missing data and delete orphan nucleotides of length < 10.
  python prepDyn.py -i aln.fasta -o output -om integer -ot 10

  # Given a FASTA alignment with hDNA sequences in sp1 and sp4, replace IUPAC N with ?
  python prepDyn.py -i aln.fasta -o output -n2q sp1,sp4

  # Given a CSV file with GenBank accession numbers, download and align sequences, identify terminal and internal missing data, trim invariants and orphan nucleotides  
  python prepDyn.py -gb accessions.csv -o output -di True -g2q semi 
""")
    # Parsing
    parser.add_argument("-i", "--input_file", help="Path to input alignment file or directory containing multiple files. Ignored if GB_input is provided.", default=None)
    parser.add_argument("-gb", "--GB_input", help="Path to a dataframe containing GenBank accession numbers (CSV/TSV). If provided, sequences will be downloaded from GenBank and aligned using MAFFT before preprocessing. Ignored if input_file is provided.", default=None)
    parser.add_argument("-if", "--input_format", help="Input file format. Options: 'fasta' (default), 'clustal', 'phylip', or any format accepted by Biopython.", default="fasta")
    parser.add_argument("-o", "--output_file", help="Path (including prefix) for output file(s)", default=None)
    parser.add_argument("-of", "--output_format", help="Output format [default: fasta]", default="fasta")
    parser.add_argument("-l", "--log", default=True, type=str2bool, help="Write time log")
    parser.add_argument("-msa", "--MSA", default=False, type=str2bool, help="Perform a MSA. Only use it if the sequences specified in input_file are unaligned. Ignore if GB_input is used.")
    parser.add_argument("-s", "--sequence_names", default=True, type=str2bool, help="Write sequence names. Useful to manage taxon sampling in POY/PhyG.")

    # Trimming
    parser.add_argument("-om", "--orphan_method", help="Method to trim orphan nucleotides. Options: 'none' (default), 'percentile' (define a threshold using percentile), 'integer' (define a threshold using an integer)", choices=["integer", "percentile"], default=None)
    parser.add_argument("-ot", "--orphan_threshold", type=int, help="Threshold integer if orphan_method='integer' (default: 10)", default=10)
    parser.add_argument("-op", "--percentile", type=float, help="Percentile of gap lengths to define the orphan threshold if orphan_method='percentile' (default: 25).", default=25.0)
    parser.add_argument("-di", "--del_inv", default=True, type=str2bool, help="Trim invariant terminal columns (default: True)")

    # Missing data
    parser.add_argument("-g2q", "--internal_method", help="Method to handle internal missing data: 'manual', 'semi', or 'none' (default)", default=None)
    parser.add_argument("-g2q_c", "--internal_column_ranges", help="Column ranges (Python list format) if internal_method=manual", type=parse_internal_column_ranges, default="all")
    parser.add_argument("-g2q_l", "--internal_leaves", help="Sequence names to apply the parameters if internal_method='semi' or 'manual'. Use 'all' (default) or comma-separated list e.g. seq1,seq2", type=parse_internal_leaves, default="all")
    parser.add_argument("-g2q_t", "--internal_threshold", type=int, help="Threshold for internal gaps if method='semi'", default=None)
    parser.add_argument("-n2q", "--n2question", type=parse_n2question_leaves, default=None, help="Replace IUPAC N with ?. Options: 'none' (default), 'all' (apply to all leaves), single leaf name, or list of leaf names ['sp1', 'sp2'].")

    # Partitioning
    parser.add_argument("-pm", "--partitioning_method", type=str, default="balanced", choices=["balanced", "conservative", "equal", "max", "None"], help="Method of partitioning. Options: (1) conservative (given blocks of contiguous invariants sorted by length, partition the n-largest block(s); define n using partitioning_round), (2) equal (insert # to divide the alignment into equal-length partitions; define the size of partitions using partitioning_size or the round of partitioning using partitioning_round); (3) max (insert '#' columns around blocks of contiguous missing data i.e. before and after every instance of '?' opening/closure); (4) balanced (insert '#' around the n largest block of missing data; define n using partitioning_round).")
    parser.add_argument("-pr", "--partitioning_round", type=parse_partitioning_round, help="Round of successive partitioning. Use it if partitioning_method is 'conservative' or 'equal'.", default=0)
    parser.add_argument("-ps", "--partitioning_size", type=int, default=None, help="Size of equal-length partitions if partitioning_method = 'equal'.")

    # Default
    parser.set_defaults(log=True)
    parser.set_defaults(del_inv=True)
    parser.set_defaults(msa=False)
    parser.set_defaults(sequence_names=True)

    args = parser.parse_args()
    # Error messages
    if not args.input_file and not args.GB_input:
        parser.error("You must provide either --input_file or --GB_input.")

    # prepDyn
    prepDyn(input_file=args.input_file,
            GB_input=args.GB_input,
            input_format=args.input_format,
            MSA=args.MSA,
            output_file=args.output_file,
            output_format=args.output_format,
            log=args.log,
            sequence_names=args.sequence_names,
            # Trimming parameters
            orphan_method=args.orphan_method,
            orphan_threshold=args.orphan_threshold,
            percentile=args.percentile,
            del_inv=args.del_inv,
            # Missing data parameters
            n2question=args.n2question,
            internal_method=args.internal_method,
            internal_column_ranges=args.internal_column_ranges,
            internal_leaves=args.internal_leaves,
            internal_threshold=args.internal_threshold,
            # Partitioning parameters
            partitioning_method=args.partitioning_method,
            partitioning_round=args.partitioning_round,
            partitioning_size=args.partitioning_size)

if __name__ == "__main__":
    main()
