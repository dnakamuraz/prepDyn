#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# auxiliary.py
# prepDyn - Copyright (C) 2025
# Daniel Y. M. Nakamura
# GNU General Public Licence version 3.0
# Contact: dani_ymn@outlook.com

COPYRIGHT=""""
addSeq.py is part of prepDyn.

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

# The program mafft should be installed in $PATH beforehand.
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
from prepDyn_auxiliary import addSeq
from Bio import SeqIO
from Bio.Seq import Seq

##################
# CUSTOM PARSERS #
##################

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
Examples: python addSeq.py --alignment aln.fas --new_seqs 12s_sp_new.fas --output aln_updated.fas --write_names
""")
    
    parser.add_argument("-a", "--alignment", type=str, required=True, help="Path to FASTA input alignment. If question marks and pound signs are present, they will be maintained.", default=None)
    parser.add_argument("-n", "--new_seqs", type=str, required=True, help="Path to FASTA input sequence(s) to be added to the alignment.", default=None)     
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to the output file with the new sequences aligned to the core alignment.", default=None)
    parser.add_argument("-w", "--write_names", type=str2bool, default=False, required=False, help="Write sequence names in a separate file, which can be used as input data in POY/PhyG to select taxon sample (default: False).")
    parser.add_argument("-ot", "--orphan_threshold", type=int, required=False, default=0, help="Threshold (int) to detect and remove orphan DNA blocks. Default = 0.")
    parser.add_argument("-n2q", "--n2question", type=parse_n2question_leaves, required=False, default=None, help="Replace IUPAC N with ?. Options: 'none' (default), 'all' (apply to all added leaves), single added leaf, or list of added leaves ['sp1', 'sp2'].")
    parser.add_argument("-g2q", "--gaps2question", default=None, required=False, type=int, help="gaps2question (int or None): Replace contiguous gap blocks larger than this threshold with '?'. Only applied to added sequences.")
    parser.add_argument("-l", "--log", type=str2bool, default=True, required=False, help="Write log tracking all operations and reporting runtime (default: true)")

    args = parser.parse_args()

    addSeq(
        alignment=args.alignment, 
        new_seqs=args.new_seqs, 
        output=args.output,
        write_names=args.write_names,
        log=args.log,
        # Trimming
        orphan_threshold=args.orphan_threshold,
        # Missing data identification
        gaps2question=args.gaps2question,
        n2question=args.n2question)

if __name__ == "__main__":
    main()