#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# prepDyn_auxiliary.py
# prepDyn - Copyright (C) 2025
# Daniel Y. M. Nakamura
# GNU General Public Licence version 3.0
# Contact: dani_ymn@outlook.com

COPYRIGHT=""""
prepDyn_auxialiary.py is part of prepDyn.

prepDyn - Data preprocessing for dynamic homology
Copyright (C) 2025 - Daniel Y. M. Nakamura
    prepDyn comes with ABSOLUTELY NO WARRANTY!
    This is a free software, and you are welcome to
    redistribute them under certain conditions
    For additional information, please visit:
    http://opensource.org/licenses/GPL-3.0
"""

###########
# MODULES #
###########

# The programs pip and Mafft should be installed beforehand.
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

# Load packages
from Bio import AlignIO, Entrez, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
from io import StringIO
import matplotlib.pyplot as plt
import numpy as np
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
from termcolor import colored
import time
from collections import Counter, defaultdict

#######################
# AUXILIARY FUNCTIONS #
#######################

def normalize_aligner(aligner):
    aligner = (aligner or "mafft").lower().replace("-", "").replace("_", "")
    if aligner == "mafft":
        return "mafft"
    if aligner == "clustalw":
        return "clustalw"
    raise ValueError("aligner must be 'mafft' or 'clustalw'.")

def run_multiple_sequence_alignment(input_fasta, output_fasta, aligner="mafft"):
    aligner = normalize_aligner(aligner)
    if aligner == "mafft":
        if shutil.which("mafft") is None:
            raise EnvironmentError("MAFFT not found. Please install MAFFT and ensure it is in your PATH.")
        with open(output_fasta, "w") as out_handle:
            process = subprocess.run(
                ["mafft", "--auto", "--quiet", input_fasta],
                stdout=out_handle,
                stderr=subprocess.PIPE,
                text=True,
                check=False
            )
        if process.returncode != 0:
            raise RuntimeError(f"MAFFT failed for {input_fasta}:\n{process.stderr}")
        return

    clustalw_exe = shutil.which("clustalw") or shutil.which("clustalw2")
    if clustalw_exe is None:
        raise EnvironmentError(
            "ClustalW not found. Install it with 'conda install bioconda::clustalw' "
            "and ensure 'clustalw' is in your PATH."
        )

    original_records = list(SeqIO.parse(input_fasta, "fasta"))
    if not original_records:
        raise ValueError(f"No sequences found in {input_fasta}.")
    nonempty_records = [
        record for record in original_records
        if any(base not in "-?" for base in str(record.seq).strip())
    ]
    empty_ids = {record.id for record in original_records if record.id not in {rec.id for rec in nonempty_records}}
    if not nonempty_records:
        raise ValueError(f"All sequences in {input_fasta} are empty or contain only gaps/missing data.")

    clustal_input = input_fasta
    clustal_output = output_fasta
    temp_paths = []
    if empty_ids:
        temp_dir = os.path.dirname(output_fasta) or tempfile.gettempdir()
        os.makedirs(temp_dir, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode="w", delete=False, dir=temp_dir, suffix=".fasta") as tmp_in:
            SeqIO.write(nonempty_records, tmp_in, "fasta")
            clustal_input = tmp_in.name
        clustal_output = os.path.join(temp_dir, f"{os.path.basename(clustal_input)}_clustalw.fasta")
        temp_paths.extend([clustal_input, clustal_output])

    try:
        if len(nonempty_records) == 1:
            SeqIO.write(nonempty_records, clustal_output, "fasta")
        else:
            process = subprocess.run(
                [
                    clustalw_exe,
                    f"-INFILE={clustal_input}",
                    f"-OUTFILE={clustal_output}",
                    "-OUTPUT=FASTA",
                    "-QUIET",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False
            )
            if process.returncode != 0:
                raise RuntimeError(
                    f"ClustalW failed for {input_fasta}:\n{process.stderr or process.stdout}"
                )

        if empty_ids:
            aligned_records = {record.id: record for record in SeqIO.parse(clustal_output, "fasta")}
            aligned_length = len(next(iter(aligned_records.values())).seq)
            restored_records = []
            for record in original_records:
                if record.id in empty_ids:
                    restored_records.append(SeqRecord(Seq("-" * aligned_length), id=record.id, description=""))
                else:
                    restored_records.append(aligned_records[record.id])
            SeqIO.write(restored_records, output_fasta, "fasta")
    finally:
        for tree_path in {
            os.path.splitext(input_fasta)[0] + ".dnd",
            os.path.splitext(output_fasta)[0] + ".dnd",
            os.path.splitext(clustal_input)[0] + ".dnd",
            os.path.splitext(clustal_output)[0] + ".dnd",
        }:
            if os.path.exists(tree_path):
                os.remove(tree_path)
        for temp_path in temp_paths:
            if os.path.exists(temp_path):
                os.remove(temp_path)

def read_delimited_rows(input_file, delimiter=",", encodings=None):
    """
    Read a CSV/TSV-style table using a small set of encoding fallbacks.

    UTF-8 remains the default, but some spreadsheets exported on older systems
    may use legacy single-byte encodings.
    """
    encodings = encodings or ("utf-8", "utf-8-sig", "latin-1", "mac_roman")
    last_error = None

    for encoding in encodings:
        try:
            with open(input_file, newline="", encoding=encoding) as file:
                reader = csv.reader(file, delimiter=delimiter)
                rows = list(reader)
            if encoding not in ("utf-8", "utf-8-sig"):
                print(f"Warning: {input_file} was decoded using fallback encoding '{encoding}'.")
            return rows
        except UnicodeDecodeError as exc:
            last_error = exc

    raise UnicodeDecodeError(
        last_error.encoding,
        last_error.object,
        last_error.start,
        last_error.end,
        f"Could not decode {input_file} using encodings: {', '.join(encodings)}"
    )

def report_duplicate_genbank_accessions(rows, csv_dir):
    """
    Inspect CSV/TSV cells before any downloads and report duplicated GenBank accessions.
    """
    if not rows or len(rows) < 2:
        return

    header = rows[0]
    accession_counts = Counter()
    accession_locations = defaultdict(list)

    for row_idx, row in enumerate(rows[1:], start=2):
        if not row:
            continue
        sample_name = row[0].strip() if row[0].strip() else f"row_{row_idx}"
        for col_idx, cell in enumerate(row[1:], start=1):
            cell = cell.strip()
            if not cell or cell.upper() == "NA" or cell == "-":
                continue

            gene_name = header[col_idx].strip() if col_idx < len(header) and header[col_idx].strip() else f"column_{col_idx + 1}"
            try:
                sources, _ = _parse_cell_sources(cell, csv_dir)
            except FileNotFoundError:
                # Defer local path validation errors to the normal workflow.
                continue

            for source in sources:
                if source["type"] != "genbank":
                    continue
                accession = source["label"]
                accession_counts[accession] += 1
                accession_locations[accession].append(f"{sample_name}/{gene_name}")

    duplicated = sorted(acc for acc, count in accession_counts.items() if count > 1)
    if not duplicated:
        return

    print("Warning: duplicated GenBank accession numbers were found in the CSV input:")
    for accession in duplicated:
        locations = ", ".join(accession_locations[accession])
        print(f"  - {accession} ({accession_counts[accession]} occurrences): {locations}")

# Function 1. Visualize colored alignments
# Define colors for each nucleotide (case insensitive)
def color_nucleotide(nucleotide):
    color_map = {
        "A": "red",
        "T": "blue",
        "G": "green",
        "C": "yellow",
        "-": "white",
        "#": "black",
        "?": "magenta"
    }
    return colored(nucleotide, color_map.get(nucleotide.upper(), "white"))
# Print the alignment
def print_colored_alignment(alignment):
    """
    Prints a DNA alignment with each nucleotide in a different color.

    Parameters:
        alignment (dict): A dictionary with sequence names as keys and DNA sequences as values.
    """
    for name, sequence in alignment.items():
        colored_seq = ''.join(color_nucleotide(n) for n in sequence)
        print(f"{name}: {colored_seq}")

# Function 2. Convert dict to MultipleSeqAlignment
def dict_to_multiple_seq_alignment(seq_dict):
    """
    Converts a dictionary of sequences into a MultipleSeqAlignment object.
    
    Args:
        seq_dict (dict): A dictionary where keys are sequence identifiers and values are sequences (strings).
        
    Returns:
        MultipleSeqAlignment: A Biopython MultipleSeqAlignment object.
    """
    # Create a list of SeqRecord objects from the dictionary
    seq_records = []
    
    for seq_id, seq_str in seq_dict.items():
        # Create a SeqRecord for each sequence
        seq = Seq(seq_str)
        seq_record = SeqRecord(seq, id=seq_id, description="")
        seq_records.append(seq_record)
    
    # Create and return the MultipleSeqAlignment object
    alignment = MultipleSeqAlignment(seq_records)
    return alignment

# Function 3. List lengths of blocks of contiguous gaps in internal and terminal positions
def list_gap_blocks_by_type(alignment, plot_distribution=False):
    """
    Identify all blocks of contiguous gaps in the DNA alignment.
    Classify gap blocks into terminal and internal blocks.
    Optionally, plot the distributions of gap block lengths for terminal and internal blocks.
    
    Args:
        alignment (dict): A dictionary with sequence IDs as keys and sequences as values.
        plot_distribution (bool): If True, plot the distributions of terminal and internal gap block lengths.
        
    Returns:
        tuple: Two lists:
            - terminal_blocks: A list of gap block lengths at the start or end of sequences.
            - internal_blocks: A list of gap block lengths in the middle of sequences.
    """
    terminal_blocks = []
    internal_blocks = []
    
    # Iterate over each sequence in the alignment
    for seq_id, sequence in alignment.items():
        sequence_length = len(sequence)
        gap_count = 0
        in_gap = False

        # Identify gap blocks and separate terminal vs internal
        for i, nucleotide in enumerate(sequence):
            if nucleotide == '-':  # We are in a gap
                if not in_gap:
                    in_gap = True
                    gap_count = 1  # Start a new gap block
                else:
                    gap_count += 1  # Continue counting the current gap block
            else:  # We encountered a non-gap nucleotide
                if in_gap:
                    # End of a gap block
                    # Check if this gap block is terminal (at start or end of the sequence)
                    if i == sequence_length or (i - gap_count == 0):
                        terminal_blocks.append(gap_count)
                    else:
                        internal_blocks.append(gap_count)
                    in_gap = False
                    gap_count = 0  # Reset gap count for the next block

        # If the sequence ends with a gap block, add it to the appropriate list
        if in_gap:
            if sequence[-1] == '-':  # Last character is a gap
                terminal_blocks.append(gap_count)
            else:
                internal_blocks.append(gap_count)

    # Optionally plot the distributions of terminal and internal blocks
    if plot_distribution:
        plt.figure(figsize=(8, 6))
        # Plotting terminal blocks
        plt.hist(terminal_blocks, bins=20, alpha=0.5, label='Terminal Blocks', color='blue')
        # Plotting internal blocks
        plt.hist(internal_blocks, bins=20, alpha=0.5, label='Internal Blocks', color='red')
        
        plt.xlabel('Gap Block Length')
        plt.ylabel('Frequency')
        plt.title('Distribution of Gap Block Lengths')
        plt.legend(loc='upper right')
        plt.show()

    return terminal_blocks, internal_blocks

# Function 4. Remove underscores in the beggining of file names
def remove_leading_underscores(file_path):
    """
    Remove leading contiguous underscores from the beginning of the file name.
    If a directory is provided, it renames all files in that directory.
    
    Args:
        file_path (str): The full path of the file or directory.
    
    Returns:
        str or None: The new file path with leading underscores removed if a file, 
                     or None if a directory (modifies in place).
    """
    if os.path.isdir(file_path):
        # If it's a directory, rename all files inside it
        for root, dirs, files in os.walk(file_path):
            for file in files:
                old_file_path = os.path.join(root, file)
                new_file_name = file.lstrip('_')
                new_file_path = os.path.join(root, new_file_name)
                
                if old_file_path != new_file_path:
                    os.rename(old_file_path, new_file_path)
        return None  # No return for directories, as the renaming is in place
    else:
        # If it's a single file, rename it
        dir_name, file_name = os.path.split(file_path)
        new_file_name = file_name.lstrip('_')
        new_file_path = os.path.join(dir_name, new_file_name)
        if file_path != new_file_path:
            os.rename(file_path, new_file_path)
        return new_file_path
    
# Function 5. Merge overlapping sequences
IUPAC_TO_BASES = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    #"U": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

BASES_TO_IUPAC = {
    frozenset(bases): code for code, bases in IUPAC_TO_BASES.items() if code != "U"
}

LOCAL_SEQUENCE_SUFFIXES = {
    ".fa", ".fas", ".fasta", ".fna", ".ffn", ".faa", ".frn",
    ".gb", ".gbk", ".genbank", ".embl", ".fastq", ".fq", ".seq",
    ".seqxml", ".nex", ".nexus", ".phy", ".phylip"
}

LOCAL_SEQUENCE_FORMAT_CANDIDATES = [
    "fasta",
    "genbank",
    "embl",
    "fastq",
    "nexus",
    "phylip",
    "seqxml",
    "ig",
    "tab"
]

EXTENSION_TO_SEQIO_FORMAT = {
    ".fa": "fasta",
    ".fas": "fasta",
    ".fasta": "fasta",
    ".fna": "fasta",
    ".ffn": "fasta",
    ".faa": "fasta",
    ".frn": "fasta",
    ".gb": "genbank",
    ".gbk": "genbank",
    ".genbank": "genbank",
    ".embl": "embl",
    ".fastq": "fastq",
    ".fq": "fastq",
    ".seqxml": "seqxml",
    ".nex": "nexus",
    ".nexus": "nexus",
    ".phy": "phylip",
    ".phylip": "phylip"
}

def _looks_like_local_sequence_path(source_text):
    expanded = os.path.expanduser(source_text.strip())
    suffix = pathlib.Path(expanded).suffix.lower()
    return (
        os.path.isabs(expanded)
        or expanded.startswith(".")
        or expanded.startswith("~")
        or "\\" in source_text
        or suffix in LOCAL_SEQUENCE_SUFFIXES
    )

def _resolve_local_sequence_path(source_text, csv_dir):
    expanded = os.path.expanduser(source_text.strip())
    candidates = []
    if os.path.isabs(expanded):
        candidates.append(expanded)
    else:
        candidates.append(os.path.join(csv_dir, expanded))
        candidates.append(expanded)

    for candidate in candidates:
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)
    return None

def _read_single_sequence_file(file_path):
    suffix = pathlib.Path(file_path).suffix.lower()
    candidate_formats = []
    preferred_format = EXTENSION_TO_SEQIO_FORMAT.get(suffix)
    if preferred_format:
        candidate_formats.append(preferred_format)
    candidate_formats.extend(
        file_format for file_format in LOCAL_SEQUENCE_FORMAT_CANDIDATES
        if file_format not in candidate_formats
    )

    parse_errors = []
    for file_format in candidate_formats:
        try:
            with open(file_path, "r") as handle:
                records = list(SeqIO.parse(handle, file_format))
        except Exception as exc:
            parse_errors.append(f"{file_format}: {exc}")
            continue

        if len(records) == 1:
            return str(records[0].seq)
        if len(records) > 1:
            raise ValueError(
                f"Local sequence file '{file_path}' contains {len(records)} records in {file_format} format. "
                "Each CSV cell must resolve to exactly one sequence per amplicon."
            )

    raise ValueError(
        f"Could not parse local sequence file '{file_path}'. Tried formats: {', '.join(candidate_formats)}. "
        f"Last errors: {'; '.join(parse_errors[:3])}"
    )

def _extract_manual_amplicon_classification(cell):
    """
    Extract manual classification marker (O) or (N) from a cell string.
    - (O) means force overlapping treatment
    - (N) means force non-overlapping treatment
    - No marker means automatic classification
    
    The markers are replaced with '|' delimiter so that multi-amplicons are properly separated.
    
    Returns:
        tuple: (cleaned_cell, classification) where classification is None, "overlapping", or "non-overlapping"
    """
    classification = None
    
    # Check for (O) marker (force overlapping) - replace with pipe delimiter
    if "(O)" in cell:
        cell = cell.replace("(O)", "|")
        classification = "overlapping"
    # Check for (N) marker (force non-overlapping) - replace with pipe delimiter
    elif "(N)" in cell:
        cell = cell.replace("(N)", "|")
        classification = "non-overlapping"
    
    return cell.strip(), classification

def _parse_cell_sources(cell, csv_dir):
    """
    Parse a CSV cell to extract sequence sources (local files or GenBank accessions).
    Supports manual classification markers: (O) for overlapping, (N) for non-overlapping.
    
    Returns:
        tuple: (sources_list, classification_override)
        - sources_list: list of dicts with "type", "label", and "path" (if local)
        - classification_override: None (automatic), "overlapping", or "non-overlapping"
    """
    # Extract manual classification marker if present
    cell, classification = _extract_manual_amplicon_classification(cell)
    
    sources = []
    for part in cell.split("|"):
        token = part.strip()
        if not token or token.upper() == "NA" or token == "-":
            continue

        local_path = _resolve_local_sequence_path(token, csv_dir)
        if local_path:
            sources.append({
                "type": "local",
                "label": token,
                "path": local_path
            })
            continue

        if _looks_like_local_sequence_path(token):
            raise FileNotFoundError(
                f"Local sequence file '{token}' was not found relative to '{csv_dir}'."
            )

        for accession in token.split("/"):
            accession = accession.strip()
            if accession and accession.upper() != "NA" and accession != "-":
                sources.append({
                    "type": "genbank",
                    "label": accession
                })

    return sources, classification

def _iupac_consensus_char(base1, base2):
    bases1 = IUPAC_TO_BASES.get(base1.upper(), {base1.upper()})
    bases2 = IUPAC_TO_BASES.get(base2.upper(), {base2.upper()})
    consensus_bases = frozenset(bases1 | bases2)
    return BASES_TO_IUPAC.get(consensus_bases, "N")

def _build_overlap_consensus(seq_left, seq_right):
    return "".join(_iupac_consensus_char(left, right) for left, right in zip(seq_left, seq_right))

def _merge_multi_via_alignment(seq_list, labels, multi_amplicon_action, aligner="mafft"):
    """
    Merge multiple sequences using MAFFT alignment to identify and handle overlaps.
    
    Parameters:
    -----------
    seq_list : list
        List of sequences to merge
    labels : list
        Labels for each sequence
    multi_amplicon_action : str
        "trim" to remove overlapping region from one sequence, "consensus" to use IUPAC consensus
    aligner : str
        Alignment tool to use (default "mafft")
    
    Returns:
    --------
    tuple: (merged_seq, merge_events, fragments_labels)
        - merged_seq: the merged sequence result
        - merge_events: list of merge event dictionaries
        - fragments_labels: list of label lists showing how sequences were combined
    """
    import tempfile
    import os
    
    merge_events = []
    
    # Create temporary directory and files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Write input FASTA
        input_fasta = os.path.join(temp_dir, "input.fasta")
        with open(input_fasta, "w") as f:
            for i, seq in enumerate(seq_list):
                f.write(f">{i}\n{seq}\n")
        
        # Run alignment
        output_fasta = os.path.join(temp_dir, "aligned.fasta")
        run_multiple_sequence_alignment(input_fasta, output_fasta, aligner=aligner)
        
        # Parse aligned sequences
        aligned_records = list(SeqIO.parse(output_fasta, "fasta"))
        if len(aligned_records) != len(seq_list):
            raise RuntimeError(f"Alignment failed: expected {len(seq_list)} sequences, got {len(aligned_records)}")
        
        # Extract alignment as list of strings (uppercase)
        alignment = [str(record.seq).upper() for record in aligned_records]
        
        # Find the overlapping region
        # Overlapping region = columns where all sequences have non-gap characters
        aln_len = len(alignment[0])
        
        # For each sequence, find start and end positions (positions with actual bases, not gaps)
        seq_positions = []
        for aln_seq in alignment:
            start_pos = None
            end_pos = None
            for pos in range(aln_len):
                if aln_seq[pos] != '-':
                    if start_pos is None:
                        start_pos = pos
                    end_pos = pos
            seq_positions.append((start_pos, end_pos))
        
        # Find overlapping column range: where multiple sequences have non-gaps
        overlap_start = None
        overlap_end = None
        sequences_with_bases_per_col = []
        
        for col in range(aln_len):
            seqs_with_bases = sum(1 for aln_seq in alignment if aln_seq[col] != '-')
            sequences_with_bases_per_col.append(seqs_with_bases)
            
            # Overlapping region is where >= 2 sequences have bases
            if seqs_with_bases >= 2:
                if overlap_start is None:
                    overlap_start = col
                overlap_end = col
        
        # If no overlap found, concatenate with Ns
        if overlap_start is None or overlap_end is None:
            merge_events.append({
                "type": "classification_override",
                "classification": "overlapping",
                "labels": labels[:],
                "reason": "Forced via (O) marker; no overlap detected in alignment, concatenating"
            })
            merged = "".join(seq_list)
            return merged, merge_events, [labels]
        
        overlap_len = overlap_end - overlap_start + 1
        
        # Extract overlap region sequences (without gaps)
        overlap_seqs = []
        for i, aln_seq in enumerate(alignment):
            overlap_region = aln_seq[overlap_start:overlap_end+1].replace('-', '')
            overlap_seqs.append(overlap_region)
        
        # Build merged sequence based on action
        # Start with first sequence up to (and including) overlap
        merged_parts = []
        
        # Add sequences before overlap
        for i in range(len(alignment)):
            start_pos, end_pos = seq_positions[i]
            # Add non-overlapping part before overlap column range
            before_overlap = alignment[i][:overlap_start].replace('-', '')
            if before_overlap:
                merged_parts.append(before_overlap)
                break  # Only take the first sequence's non-overlapping prefix
        
        # Handle the overlap region
        if multi_amplicon_action == "consensus" and len(overlap_seqs) > 1:
            # Build consensus from all overlap sequences
            max_overlap_len = max(len(os) for os in overlap_seqs)
            consensus = []
            for pos in range(max_overlap_len):
                chars = []
                for os in overlap_seqs:
                    if pos < len(os):
                        chars.append(os[pos])
                if chars:
                    if len(chars) == 1:
                        consensus.append(chars[0])
                    else:
                        consensus.append(_iupac_consensus_char(chars[0], chars[1] if len(chars) > 1 else chars[0]))
            merged_parts.append("".join(consensus))
        else:
            # Trim action: just use the first sequence's overlap
            merged_parts.append(overlap_seqs[0] if overlap_seqs else "")
        
        # Add sequences after overlap
        for i in range(len(alignment)):
            start_pos, end_pos = seq_positions[i]
            # Add non-overlapping part after overlap column range
            after_overlap = alignment[i][overlap_end+1:].replace('-', '')
            if after_overlap:
                merged_parts.append(after_overlap)
                break  # Only take the last sequence's non-overlapping suffix
        
        merged = "".join(merged_parts)
        
        # Log the merge event
        merge_events.append({
            "type": "classification_override",
            "classification": "overlapping",
            "labels": labels[:],
            "reason": f"Forced via (O) marker; merged via MAFFT alignment with {multi_amplicon_action} action (overlap length: {overlap_len}bp)"
        })
        
        return merged, merge_events, [labels]

def merge_overlapping_sequences(seq_list, multi_amplicon_min_overlap=10, multi_amplicon_mismatch_rate=0.05, multi_amplicon_action="trim", labels=None, forced_classification=None, aligner="mafft"):
    """
    Tries to merge a list of sequences by identifying overlaps.
    If some sequences don't overlap, they remain separate.
    Returns merged sequences, structured merge events, and the final amplicon groups.
    
    Parameters:
    -----------
    seq_list : list
        List of sequences to merge
    multi_amplicon_min_overlap : int
        Minimum overlap length to consider as overlapping
    multi_amplicon_mismatch_rate : float
        Maximum mismatch rate allowed in overlaps
    multi_amplicon_action : str
        Action for overlaps: "trim" or "consensus"
    labels : list
        Labels for each sequence
    forced_classification : str or None
        Manual override for classification:
        - None: automatic classification based on overlap detection
        - "overlapping": force attempt to merge sequences (uses relaxed overlap detection)
        - "non-overlapping": keep sequences separate without attempting merge
    """
    merge_log = []
    if not seq_list:
        return [], merge_log, []
    if multi_amplicon_action not in {"trim", "consensus"}:
        raise ValueError("multi_amplicon_action must be either 'trim' or 'consensus'.")

    if labels is None:
        labels = [f"amplicon_{idx + 1}" for idx in range(len(seq_list))]
    if len(labels) != len(seq_list):
        raise ValueError("labels and seq_list must have the same length.")

    # If forced to non-overlapping, return sequences as-is
    if forced_classification == "non-overlapping":
        for idx, label in enumerate(labels):
            merge_log.append({
                "type": "classification_override",
                "classification": "non-overlapping",
                "labels": [label],
                "reason": "Forced to non-overlapping via (N) marker in CSV"
            })
        return seq_list, merge_log, [[label] for label in labels]

    # When forced to overlapping, use MAFFT alignment to robustly detect and merge overlaps
    if forced_classification == "overlapping":
        try:
            merged_seq, align_events, align_groups = _merge_multi_via_alignment(
                seq_list, labels, multi_amplicon_action, aligner=aligner
            )
            merge_log.extend(align_events)
            return [merged_seq], merge_log, align_groups
        except Exception as e:
            # If alignment fails, fall back to raw sequence detection with relaxed parameters
            print(f"Warning: Alignment-based merge failed ({e}), falling back to raw sequence detection", file=sys.stderr)

    # When NOT forced to overlapping, use automatic classification based on raw sequence overlap detection
    effective_min_overlap = multi_amplicon_min_overlap
    effective_mismatch_rate = multi_amplicon_mismatch_rate

    # Operate on upper case to ensure case-insensitive matching
    fragments = [
        {"seq": seq.upper(), "labels": [label]}
        for seq, label in zip(seq_list, labels)
    ]
        
    merged_any = True
    while merged_any:
        merged_any = False
        for i in range(len(fragments)):
            for j in range(i + 1, len(fragments)):
                fragment1 = fragments[i]
                fragment2 = fragments[j]
                seq1 = fragment1["seq"]
                seq2 = fragment2["seq"]
                
                # 1. Check containment first
                if seq1 in seq2:
                    merge_log.append({
                        "type": "containment",
                        "classification": "overlapping",
                        "removed_labels": fragment1["labels"][:],
                        "kept_labels": fragment2["labels"][:],
                        "deleted_string": seq1,
                        "deleted_length": len(seq1)
                    })
                    fragments[j]["labels"].extend(fragment1["labels"])
                    fragments.pop(i)
                    merged_any = True
                    break
                if seq2 in seq1:
                    merge_log.append({
                        "type": "containment",
                        "classification": "overlapping",
                        "removed_labels": fragment2["labels"][:],
                        "kept_labels": fragment1["labels"][:],
                        "deleted_string": seq2,
                        "deleted_length": len(seq2)
                    })
                    fragments[i]["labels"].extend(fragment2["labels"])
                    fragments.pop(j)
                    merged_any = True
                    break
                    
                # 2. Check overlap: seq1 suffix with seq2 prefix
                overlap_len = 0
                best_merged = None
                deleted_string = ""
                mismatch_count = 0
                merge_direction = None
                
                len1, len2 = len(seq1), len(seq2)
                max_k = min(len1, len2)
                
                # Check seq1 suffix matching seq2 prefix
                for k in range(max_k, effective_min_overlap - 1, -1):
                    s1_suffix = seq1[-k:]
                    s2_prefix = seq2[:k]
                    
                    if s1_suffix == s2_prefix:
                        overlap_len = k
                        best_merged = seq1 + seq2[k:]
                        deleted_string = seq2[:k]
                        mismatch_count = 0
                        merge_direction = "seq1_suffix_seq2_prefix"
                        break
                        
                    max_mismatches = max(2, int(effective_mismatch_rate * k))
                    mismatches = 0
                    for a, b in zip(s1_suffix, s2_prefix):
                        if a != b:
                            mismatches += 1
                            if mismatches > max_mismatches:
                                break
                    if mismatches <= max_mismatches:
                        overlap_len = k
                        best_merged = seq1 + seq2[k:]
                        deleted_string = seq2[:k]
                        mismatch_count = mismatches
                        merge_direction = "seq1_suffix_seq2_prefix"
                        break
                        
                # Check seq2 suffix matching seq1 prefix
                if not best_merged:
                    for k in range(max_k, effective_min_overlap - 1, -1):
                        s2_suffix = seq2[-k:]
                        s1_prefix = seq1[:k]
                        
                        if s2_suffix == s1_prefix:
                            overlap_len = k
                            best_merged = seq2 + seq1[k:]
                            deleted_string = seq1[:k]
                            mismatch_count = 0
                            merge_direction = "seq2_suffix_seq1_prefix"
                            break
                            
                        max_mismatches = max(2, int(effective_mismatch_rate * k))
                        mismatches = 0
                        for a, b in zip(s2_suffix, s1_prefix):
                            if a != b:
                                mismatches += 1
                                if mismatches > max_mismatches:
                                    break
                        if mismatches <= max_mismatches:
                            overlap_len = k
                            best_merged = seq2 + seq1[k:]
                            deleted_string = seq1[:k]
                            mismatch_count = mismatches
                            merge_direction = "seq2_suffix_seq1_prefix"
                            break
                        
                if best_merged:
                    overlap_seq_left = seq1[-overlap_len:] if merge_direction == "seq1_suffix_seq2_prefix" else seq2[-overlap_len:]
                    overlap_seq_right = seq2[:overlap_len] if merge_direction == "seq1_suffix_seq2_prefix" else seq1[:overlap_len]
                    consensus_string = _build_overlap_consensus(overlap_seq_left, overlap_seq_right)
                    if multi_amplicon_action == "consensus":
                        if merge_direction == "seq1_suffix_seq2_prefix":
                            best_merged = seq1[:-overlap_len] + consensus_string + seq2[overlap_len:]
                        else:
                            best_merged = seq2[:-overlap_len] + consensus_string + seq1[overlap_len:]

                    merge_log.append({
                        "type": "overlap",
                        "classification": "overlapping",
                        "left_labels": fragment1["labels"][:],
                        "right_labels": fragment2["labels"][:],
                        "overlap_length": overlap_len,
                        "action": multi_amplicon_action,
                        "deleted_string": deleted_string,
                        "deleted_length": len(deleted_string),
                        "consensus_string": consensus_string,
                        "mismatch_count": mismatch_count,
                        "direction": merge_direction
                    })
                    fragments.pop(j)
                    if merge_direction == "seq1_suffix_seq2_prefix":
                        fragments[i] = {
                            "seq": best_merged,
                            "labels": fragment1["labels"] + fragment2["labels"]
                        }
                    else:
                        fragments[i] = {
                            "seq": best_merged,
                            "labels": fragment2["labels"] + fragment1["labels"]
                        }
                    merged_any = True
                    break
            if merged_any:
                break
                
    return [fragment["seq"] for fragment in fragments], merge_log, [fragment["labels"][:] for fragment in fragments]

def _format_deleted_string_for_log(sequence, max_len=120):
    if len(sequence) <= max_len:
        return sequence
    flank_len = max((max_len - 3) // 2, 1)
    return f"{sequence[:flank_len]}...{sequence[-flank_len:]}"

def _format_amplicon_groups_for_log(groups):
    return "; ".join("/".join(group) for group in groups)

def _build_source_origin_log(local_sources, genbank_sources):
    lines = []
    if local_sources:
        lines.append(f"Local file sources: {' | '.join(local_sources)}")
    if genbank_sources:
        lines.append(f"GenBank sources: {'/'.join(genbank_sources)}")
    return lines

def _build_multi_amplicon_log(requested_sources, loaded_sources, failed_sources, merge_events, final_groups):
    lines = [f"Multi-amplicon input sources: {' | '.join(requested_sources)}"]

    if failed_sources:
        lines.append(f"Sources not loaded: {' | '.join(failed_sources)}")

    if len(loaded_sources) <= 1:
        lines.append("Overlap classification: not enough loaded amplicons to evaluate overlap.")
        return lines

    # Track if we've seen forced non-overlapping classification
    has_forced_non_overlapping = False

    if merge_events:
        for event in merge_events:
            if event["type"] == "classification_override":
                # Handle forced classification events with simplified format
                sources_str = "/".join(event["labels"])
                if event["classification"] == "non-overlapping":
                    lines.append(f"Non-overlapping (forced): {sources_str}")
                    has_forced_non_overlapping = True
                elif event["classification"] == "overlapping":
                    lines.append(f"Overlapping (forced): {sources_str}")
            elif event["type"] == "containment":
                removed = "/".join(event["removed_labels"])
                kept = "/".join(event["kept_labels"])
                deleted = _format_deleted_string_for_log(event["deleted_string"])
                lines.append(
                    f"Overlapping: {removed} was fully contained in {kept}; deleted string ({event['deleted_length']} bp): {deleted}"
                )
            else:
                left = "/".join(event["left_labels"])
                right = "/".join(event["right_labels"])
                mismatch_note = ""
                if event["mismatch_count"] > 0:
                    mismatch_note = f", mismatches allowed/used: {event['mismatch_count']}"
                if event.get("action") == "consensus":
                    consensus = _format_deleted_string_for_log(event["consensus_string"])
                    lines.append(
                        f"Overlapping: {left} with {right}; overlap {event['overlap_length']} bp{mismatch_note}; consensus string ({event['overlap_length']} bp): {consensus}"
                    )
                else:
                    deleted = _format_deleted_string_for_log(event["deleted_string"])
                    lines.append(
                        f"Overlapping: {left} with {right}; overlap {event['overlap_length']} bp{mismatch_note}; deleted string ({event['deleted_length']} bp): {deleted}"
                    )

    # Only report non-overlapping groups if there are actually multiple groups remaining
    # (meaning some amplicons did not merge). Skip if forced non-overlapping was already reported.
    if len(final_groups) > 1 and not has_forced_non_overlapping:
        lines.append(f"Non-overlapping: {_format_amplicon_groups_for_log(final_groups)}")

    return lines

############################################   
# MAIN FUNCTIONS: STEP 1. DATA COLLECTION  #
############################################

def GB2MSA_1(input_file, output_prefix, delimiter=',', write_names=True, multi_amplicon_min_overlap=10, multi_amplicon_mismatch_rate=0.05, multi_amplicon_action="trim", aligner="mafft"):
    """
    Loads locus sequences from GenBank accession numbers and/or local sequence files listed in a CSV/TSV
    file and aligns them by gene using the selected aligner. If two fragments of the same locus are concatenated with no overlap between them, the space
    between them will be treaed as missing data (15 Ws will flank these blocks of missing data, which will
    be used to track these regions and be replaced with question marks in GB2MSA_2).

    Parameters:
    -----------
    input_file : str
        Path to the CSV or TSV input file. The first column should contain sequence names (sample identifiers).
        The first row should contain gene names starting from the second column. Cells may contain GenBank
        accession numbers (one or more separated by slashes and/or pipes), local sequence-file paths (absolute
        or relative), or both across the same table. Multi-amplicon local files should be separated with
        pipes ("|").
        "NA", empty cells, or dashes are ignored.

    output_prefix : str
        Prefix used for naming intermediate FASTA files and final aligned output files.

    delimiter : str, optional (default=',')
        Delimiter used in the input file (e.g., ',' for CSV or '\t' for TSV).

    write_names : bool, optional (default=True)
        If True, writes a TXT file listing all sequence names (from the first column).

    multi_amplicon_min_overlap : int, optional (default=10)
        Minimum overlap length to merge multi-amplicons downlaoded from GenBank.

    multi_amplicon_mismatch_rate : float, optional (default=0.05)
        Maximum mismatch rate allowed when merging overlapping multi-amplicons.

    multi_amplicon_action : str, optional (default="trim")
        If "trim", remove one copy of overlapping regions when merging multi-amplicons.
        If "consensus", keep the overlap and replace it with an IUPAC consensus string.

    aligner : str, optional (default="mafft")
        Multiple sequence aligner to use: "mafft" or "clustalw".

    Returns:
    --------
    aligned_files : list of str
        List of file paths to the aligned FASTA files for each gene.
    gb1_logs : dict
        A dictionary containing merge logs for each sequence and gene.
    gene_names : list
        List of gene names parsed from the input file.
    """
    csv_dir = os.path.dirname(os.path.abspath(input_file)) or os.getcwd()

    # Open the input CSV/TSV file with common spreadsheet encoding fallbacks.
    rows = read_delimited_rows(input_file, delimiter=delimiter)
    report_duplicate_genbank_accessions(rows, csv_dir)

    # Replace spaces in sequence names with underscores
    sequence_names = [row[0].replace(" ", "_") for row in rows[1:]]
    
    # Extract gene names from the header row (excluding first column) and strip whitespace
    gene_names = [name.strip() for name in rows[0][1:]]
    
    # Extract gene accession data for each sequence (excluding first column)
    gene_columns = [row[1:] for row in rows[1:]]

    # If requested, write the sequence names to a text file
    if write_names:
        names_file = f"{output_prefix}_sequence_names.txt"
        with open(names_file, 'w') as nf:
            for name in sequence_names:
                nf.write(f"{name}\n")

    aligned_files = []  # List to store paths of aligned output files
    gb1_logs = {gene: {} for gene in gene_names}

    # Iterate over each gene (i.e., each column after the first)
    for gene_idx, gene_name in enumerate(gene_names):
        fasta_file = f"{output_prefix}_{gene_name}.fasta"  # Name of temporary FASTA file
        
        with open(fasta_file, 'w') as fasta_out:
            # Iterate through each row (sample/sequence)
            for i, seq_name in enumerate(sequence_names):
                cell = gene_columns[i][gene_idx].strip()
                
                # Skip cells with missing data ("NA", empty, or dash)
                if cell.upper() == "NA" or not cell or cell == "-":
                    continue

                sources, amplicon_classification = _parse_cell_sources(cell, csv_dir)
                if not sources:
                    continue
                sequences = []
                loaded_sources = []
                failed_sources = []
                local_sources = [source["label"] for source in sources if source["type"] == "local"]
                genbank_sources = [source["label"] for source in sources if source["type"] == "genbank"]
                requested_sources = [source["label"] for source in sources]

                for source in sources:
                    if source["type"] == "local":
                        try:
                            sequences.append(_read_single_sequence_file(source["path"]))
                            loaded_sources.append(source["label"])
                        except Exception as e:
                            failed_sources.append(source["label"])
                            print(f"Error loading local sequence file {source['label']}: {e}")
                    else:
                        acc = source["label"]
                        try:
                            handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
                            seq_record = SeqIO.read(handle, "fasta")
                            handle.close()
                            sequences.append(str(seq_record.seq))
                            loaded_sources.append(acc)
                        except Exception as e:
                            failed_sources.append(acc)
                            print(f"Error fetching {acc}: {e}")

                # Merge overlapping sequences if possible
                sequences, m_log, final_groups = merge_overlapping_sequences(
                    sequences,
                    multi_amplicon_min_overlap=multi_amplicon_min_overlap,
                    multi_amplicon_mismatch_rate=multi_amplicon_mismatch_rate,
                    multi_amplicon_action=multi_amplicon_action,
                    labels=loaded_sources,
                    forced_classification=amplicon_classification,
                    aligner=aligner
                )

                gb1_logs[gene_name][seq_name] = {
                    "requested_sources": requested_sources,
                    "loaded_sources": loaded_sources,
                    "failed_sources": failed_sources,
                    "local_sources": local_sources,
                    "genbank_sources": genbank_sources,
                    "merge_events": m_log,
                    "final_groups": final_groups
                }

                # Combine multiple sequences with 'W' delimiters to mark junctions (to be handled later)
                combined_seq = "WWWWWWWWWWWWWWW".join(sequences)
                
                # Write the sequence to the FASTA file with its name as header
                fasta_out.write(f">{seq_name}\n{combined_seq}\n")

        # Define the name for the output alignment file
        aligned_file = f"{output_prefix}_{gene_name}_aligned.fasta"

        # Run the selected aligner on the generated FASTA file and save the alignment
        run_multiple_sequence_alignment(fasta_file, aligned_file, aligner=aligner)

        # Append the aligned file path to the result list
        aligned_files.append(aligned_file)

    return aligned_files, gb1_logs, gene_names  

def GB2MSA_2(alignment_file):
    """
    If internal missing data were identified by GB2MSA_1, 15 Ws flank the blocks of missing data. 
    For each sequence in the alignment, GB2MSA_2:
    - Replaces internal blocks of 15 'w' or spaced 'w' (e.g., w-w-w) with question marks.
    - Removes columns with only '?' or '-' in all rows.
    - Replaces dash blocks flanked by a nucleotide and a question mark (e.g., A??--AAC) with 
    question marks. These dashes can be artifacts from spaced 'w' and are missing data.

    Parameters:
    -----------
    alignment_file : str
        Path to the aligned FASTA file to be processed.
    
    Returns:
    --------
    str
        Path to the cleaned alignment file.
    """
    alignment = list(SeqIO.parse(alignment_file, "fasta"))
    updated_records = []

    # Step 1: Replace exact block of 15 'w' or 'W' with 15 '?'
    for record in alignment:
        seq = str(record.seq)
        seq_cleaned = seq.replace("w" * 15, "?" * 15).replace("W" * 15, "?" * 15)
        record.seq = Seq(seq_cleaned)
        updated_records.append(record)

    # Step 2: Replace non-contiguous 'w' blocks (e.g., w-w-w) with '?'
    for record in updated_records:
        chars = list(str(record.seq))
        i = 0
        while i < len(chars):
            if chars[i].lower() == 'w':
                count = 1
                indices = [i]
                j = i + 1
                while j < len(chars) and count < 15:
                    if chars[j] == '-':
                        indices.append(j)
                    elif chars[j].lower() == 'w':
                        indices.append(j)
                        count += 1
                    else:
                        break
                    j += 1
                if count == 15:
                    for idx in indices:
                        chars[idx] = '?'
                i = j
            else:
                i += 1
        record.seq = Seq(''.join(chars))

    # Step 3: Remove columns with only '?' or '-' in all rows
    sequences = [list(str(record.seq)) for record in updated_records]
    
    if not sequences:
        cleaned_file = alignment_file.replace(".fasta", "_GB2MSA.fasta")
        with open(cleaned_file, "w") as out_handle:
            pass
        return cleaned_file

    if len(set(len(seq) for seq in sequences)) > 1:
        raise ValueError("Sequences are not of the same length!")

    valid_columns = []
    for i in range(len(sequences[0])):
        column = [seq[i] for seq in sequences]
        if any(base not in ['?', '-'] for base in column):
            valid_columns.append(i)

    # Rebuild records with cleaned sequences
    cleaned_records = []
    for record in updated_records:
        cleaned_seq = ''.join(str(record.seq)[i] for i in valid_columns)
        record.seq = Seq(cleaned_seq)
        cleaned_records.append(record)

    # Step 4: Replace dash blocks flanked by nucleotide and '?' with '?'
    for record in cleaned_records:
        chars = list(str(record.seq))
        i = 0
        while i < len(chars):
            if chars[i] == '-':
                start = i
                while i < len(chars) and chars[i] == '-':
                    i += 1
                end = i - 1

                # Check flanking characters safely
                left = chars[start - 1] if start > 0 else ''
                right = chars[end + 1] if end + 1 < len(chars) else ''

                if ((left in 'ACGTacgt' and right in '?Nn') or 
                    (left in '?Nn' and right in 'ACGTacgt')):
                    for j in range(start, end + 1):
                        chars[j] = '?'
            else:
                i += 1
        record.seq = Seq(''.join(chars))

    # Step 5: Write to cleaned output
    cleaned_file = alignment_file.replace(".fasta", "_GB2MSA.fasta")
    with open(cleaned_file, "w") as out_handle:
        SeqIO.write(cleaned_records, out_handle, "fasta")

    return cleaned_file

def GB2MSA_3(alignment_dict, orphan_threshold=10, log=False):
    """
    Replaces specific blocks in a DNA alignment with '?':
    - Contiguous gap blocks (length >= orphan_threshold) that are adjacent
      to contiguous nucleotide blocks (length < orphan_threshold), where the
      other side of that nucleotide block touches a '?'

    Parameters:
    alignment_dict (dict): {sequence_name: aligned_sequence}
    orphan_threshold (int): Minimum size to define a valid gap block
    log (bool): If True, print the start and end positions of replaced blocks

    Returns:
    dict: Cleaned alignment with selective replacements
    """
    cleaned_alignment = {}

    for seq_name, original_seq in alignment_dict.items():
        sequence = list(original_seq)
        seq_len = len(sequence)
        updated_seq = ''.join(sequence)

        gap_matches = list(re.finditer(r'-+', updated_seq))
        valid_gap_blocks = [(m.start(), m.end() - 1) for m in gap_matches
                            if (m.end() - m.start()) >= orphan_threshold]

        to_replace = set()

        if log:
            print(f"Sequence: {seq_name}")

        for gap_start, gap_end in valid_gap_blocks:
            replaced = False

            # Check left nucleotide orphan
            left_end = gap_start - 1
            i = left_end
            while i >= 0 and updated_seq[i] not in '-?':
                i -= 1
            left_start = i + 1
            left_len = left_end - left_start + 1

            if left_len > 0 and left_len < orphan_threshold and i >= 0 and updated_seq[i] == '?':
                to_replace.update(range(left_start, left_end + 1))
                to_replace.update(range(gap_start, gap_end + 1))
                replaced = True
                if log:
                    print(f"  Gap block: {gap_start}-{gap_end}")
                    print(f"  Left orphan nucleotide block: {left_start}-{left_end}")

            if not replaced:
                # Check right nucleotide orphan
                right_start = gap_end + 1
                i = right_start
                while i < seq_len and updated_seq[i] not in '-?':
                    i += 1
                right_end = i - 1
                right_len = right_end - right_start + 1

                if right_len > 0 and right_len < orphan_threshold and i < seq_len and updated_seq[i] == '?':
                    to_replace.update(range(right_start, right_end + 1))
                    to_replace.update(range(gap_start, gap_end + 1))
                    if log:
                        print(f"  Gap block: {gap_start}-{gap_end}")
                        print(f"  Right orphan nucleotide block: {right_start}-{right_end}")

        # Replace in sequence
        for i in to_replace:
            sequence[i] = '?'

        cleaned_alignment[seq_name] = ''.join(sequence)

    return cleaned_alignment

def GB2MSA_4(alignment_dict):
    cleaned_alignment = {}

    for seq_name, seq in alignment_dict.items():
        sequence = list(seq)

        # Use regex to find blocks with at least 15 characters of w/W/- combined
        # but must contain at least 15 w or W
        pattern = re.finditer(r'([wW\-]{15,})', ''.join(sequence))

        for match in pattern:
            block = match.group()
            start, end = match.start(), match.end()

            # Count number of w/W in the block
            w_count = sum(1 for c in block if c in 'wW')
            if w_count >= 15:
                for i in range(start, end):
                    sequence[i] = '?'

        cleaned_alignment[seq_name] = ''.join(sequence)

    return cleaned_alignment

####################################
# MAIN FUNCTIONS: STEP 2. TRIMMING #
####################################

def remove_all_gap_columns(alignment, empty=True):
    """
    Remove columns from the alignment where all sequences have only gaps ('-' or '?').
    Optionally remove sequences composed entirely of gaps.

    Parameters:
        alignment (dict): A dictionary where keys are sequence names and values are sequence strings.
        empty (bool): If True, remove sequences composed only of gap characters.

    Returns:
        dict: Updated alignment with gap-only columns and optionally gap-only sequences removed.
    """
    if not alignment:
        return {}

    gap_chars = {'-', '?'}
    sequences = list(alignment.values())
    seq_length = len(sequences[0])  # Assumes all sequences are same length

    # Identify columns where all sequences have a gap character
    columns_to_remove = [
        col for col in range(seq_length)
        if all(seq[col] in gap_chars for seq in sequences)
    ]

    # Remove the identified columns from each sequence
    alignment = {
        name: ''.join(seq[col] for col in range(seq_length) if col not in columns_to_remove)
        for name, seq in alignment.items()
    }

    # Optionally remove sequences that are entirely gap characters
    if empty:
        alignment = {
            name: seq for name, seq in alignment.items()
            if any(base not in gap_chars for base in seq)
        }

    return alignment

def calculate_orphan_threshold_from_percentile(alignment, percentile=25, log=False, terminal_only=True):
    """
    Calculate orphan threshold based on a specific percentile of gap lengths in the alignment.
    
    Args:
        alignment (dict): A dictionary with sequence IDs as keys and sequences as values.
        percentile (float): The percentile to use for setting the orphan threshold (e.g., 75 for the 75th percentile).
        log (bool): If True, print the list of gap lengths and the computed orphan threshold.
        terminal_only (bool): If True, only consider gap blocks at the terminal positions (start and end of sequences).
        
    Returns:
        int: Calculated orphan threshold based on the specified percentile.
    """
    gap_lengths = []

    # Loop through each sequence to calculate gap lengths
    for sequence in alignment.values():
        sequence_length = len(sequence)
        gap_count = 0
        in_gap = False

        # If terminal_only is True, we'll check only the first and last blocks of gaps
        if terminal_only:
            # Check for gaps starting at the beginning of the sequence
            if sequence[0] == '-':
                gap_count = 1
                in_gap = True
                for nucleotide in sequence[1:]:
                    if nucleotide == '-':
                        gap_count += 1
                    else:
                        break  # End of the first terminal gap block
                gap_lengths.append(gap_count)  # Add the terminal gap block at the start
                
            # Check for gaps starting at the end of the sequence
            gap_count = 0
            in_gap = False
            if sequence[-1] == '-':
                gap_count = 1
                in_gap = True
                for nucleotide in reversed(sequence[:-1]):
                    if nucleotide == '-':
                        gap_count += 1
                    else:
                        break  # End of the last terminal gap block
                gap_lengths.append(gap_count)  # Add the terminal gap block at the end
        else:
            # Loop through the sequence to find contiguous blocks of gaps (not restricted to terminals)
            for nucleotide in sequence:
                if nucleotide == '-':
                    if not in_gap:
                        in_gap = True  # Start of a new gap block
                        gap_count = 1  # Start counting the length of the new gap block
                    else:
                        gap_count += 1  # Continue counting the gap block length
                else:
                    if in_gap:
                        gap_lengths.append(gap_count)  # End of a gap block, save its length
                        in_gap = False
                        gap_count = 0  # Reset the gap count for the next block

            # If the sequence ends with a gap block, make sure to add the last block's length
            if in_gap:
                gap_lengths.append(gap_count)
        
    # Calculate the specified percentile of gap lengths
    orphan_threshold = int(np.percentile(gap_lengths, percentile))  # Use the given percentile (e.g., 75th percentile)
    
    # Print the list of gap lengths
    if log:
        print("List of gap lengths:", gap_lengths)
        print("Orphan threshold:", orphan_threshold)
        
    return orphan_threshold

def delete_orphan_nucleotides2(alignment, orphan_threshold, log_changes=False, orphan_action="trim", orphan_method="integer", orphan_limit=0.05):
    """
    Iteratively eliminates or pushes orphan nucleotide blocks from the start and end of each sequence.
    An orphan block is a short contiguous run of nucleotides near the terminal ends separated by many gaps.

    Args:
        alignment (dict): {sequence_id: sequence_string}
        orphan_threshold (int): Max block length and max gap tolerance (for static methods). 
                                Also acts as the maximum cap limit for adaptive methods.
        log_changes (bool): Whether to return a log of changes made.
        orphan_action (str): Action to perform ('trim' to delete or 'push' to move adjacent to next block).
        orphan_method (str): "integer", "percentile", "adaptive_1", "adaptive_2", or "adaptive_3".
                             - adaptive_1: no gap check, no shared homology check
                             - adaptive_2: gap check, no shared homology check
                             - adaptive_3: gap check + shared homology check
        orphan_limit (float): Modification limit as a fraction of sequence length for adaptive methods (default: 0.05 = 5%).

    Returns:
        dict: Cleaned alignment.
        str (optional): Log of changes made.
    """
    change_log = []

    def find_blocks(seq):
        """Find contiguous non-gap blocks as (start, end) tuples."""
        blocks = []
        i = 0
        while i < len(seq):
            if seq[i] != '-':
                start = i
                while i < len(seq) and seq[i] != '-':
                    i += 1
                end = i
                blocks.append((start, end))
            else:
                i += 1
        return blocks

    def find_terminal_blocks(seq):
        """Find only the outer two non-gap blocks on each side."""
        seq_len = len(seq)
        left_blocks = []
        i = 0
        while i < seq_len and len(left_blocks) < 2:
            if seq[i] != '-':
                start = i
                while i < seq_len and seq[i] != '-':
                    i += 1
                left_blocks.append((start, i))
            else:
                i += 1

        if len(left_blocks) < 2:
            return None

        right_blocks = []
        i = seq_len - 1
        while i >= 0 and len(right_blocks) < 2:
            if seq[i] != '-':
                end = i + 1
                while i >= 0 and seq[i] != '-':
                    i -= 1
                right_blocks.append((i + 1, end))
            else:
                i -= 1

        return left_blocks[0], left_blocks[1], right_blocks[1], right_blocks[0]

    # --- REFINEMENT STEP: Prevent deletion of shared isolated blocks ---
    def has_shared_homology(start, end, block_str, current_seq_id):
        """
        Checks if at least one other sequence shares this exact block at this exact position,
        AND ensures that the block is similarly isolated (flanked by gaps/missing data) in the other sequence.
        """
        for other_id, other_seq in alignment.items():
            if other_id != current_seq_id:
                if other_seq[start:end].lower() == block_str.lower():
                    # Ensure it is actually an isolated block in the other sequence too
                    left_is_gap = (start == 0) or (other_seq[start - 1] in ['-', '?', '#'])
                    right_is_gap = (end == len(other_seq)) or (other_seq[end] in ['-', '?', '#'])
                    
                    if left_is_gap and right_is_gap:
                        return True
        return False
    # -------------------------------------------------------------------

    adaptive_methods = ["adaptive_1", "adaptive_2", "adaptive_3"]

    # For tracking maximum allowed modifications in adaptive methods
    seq_mods = {seq_id: 0 for seq_id in alignment}
    # orphan_limit multiplied by sequence length, with a generous baseline of 3 allowed nucleotides for short sequences
    seq_limits = {
        seq_id: max(3, int(orphan_limit * len(seq.replace('-', '').replace('?', '').replace('#', '')))) 
        for seq_id, seq in alignment.items()
    }

    if orphan_method in adaptive_methods:
        current_threshold = 1
    else:
        current_threshold = orphan_threshold

    alignment_changed = True
    
    # Adaptive loop stopping condition for adaptive methods
    while alignment_changed or (orphan_method in adaptive_methods and current_threshold <= orphan_threshold):
        alignment_changed = False

        for seq_id, sequence in alignment.items():
            seq_list = list(sequence)
            changed_this_seq = False
            
            while True:
                terminal_blocks = find_terminal_blocks(seq_list)
                if terminal_blocks is None:
                    break
                
                dynamic_threshold = current_threshold

                # Left side (strictly targets the outermost left block)
                first_block, next_block, prev_block, last_block = terminal_blocks
                first_start, first_end = first_block
                next_start = next_block[0]
                gap_count_left = seq_list[first_end:next_start].count('-')
                left_len = first_end - first_start

                # Right side (strictly targets the outermost right block)
                last_start, last_end = last_block
                prev_end = prev_block[1]
                gap_count_right = seq_list[prev_end:last_start].count('-')
                right_len = last_end - last_start

                made_change_left = False
                made_change_right = False

                # Evaluate Left side condition
                if orphan_method == "adaptive_1":
                    cond_left = (left_len <= dynamic_threshold) and (seq_mods[seq_id] + left_len <= seq_limits[seq_id])
                elif orphan_method in ["adaptive_2", "adaptive_3"]:
                    cond_left = (left_len <= dynamic_threshold) and (gap_count_left > dynamic_threshold) and (seq_mods[seq_id] + left_len <= seq_limits[seq_id])
                else:
                    cond_left = (left_len < dynamic_threshold) and (gap_count_left > dynamic_threshold)

                # APPLY REFINEMENT VETO (only for adaptive_3)
                if cond_left and orphan_method == "adaptive_3":
                    block_str_left = ''.join(seq_list[first_start:first_end])
                    if has_shared_homology(first_start, first_end, block_str_left, seq_id):
                        cond_left = False # Veto! Another terminal has this exact isolated block.

                if cond_left:
                    if orphan_action == "trim":
                        deleted = ''.join(seq_list[first_start:first_end])
                        seq_list[first_start:first_end] = ['-'] * left_len
                        if log_changes:
                            change_log.append(f"{seq_id}: Left block {first_start}-{first_end} deleted ('{deleted}')")
                    elif orphan_action == "push":
                        block_str = ''.join(seq_list[first_start:first_end])
                        seq_list[first_start:next_start] = ['-'] * gap_count_left + list(block_str)
                        if log_changes:
                            change_log.append(f"{seq_id}: Left block {first_start}-{first_end} pushed to {next_start - left_len}-{next_start} ('{block_str}')")
                    
                    seq_mods[seq_id] += left_len
                    made_change_left = True
                    changed_this_seq = True

                if made_change_left:
                    continue

                # Evaluate Right side condition
                if orphan_method == "adaptive_1":
                    cond_right = (right_len <= dynamic_threshold) and (seq_mods[seq_id] + right_len <= seq_limits[seq_id])
                elif orphan_method in ["adaptive_2", "adaptive_3"]:
                    cond_right = (right_len <= dynamic_threshold) and (gap_count_right > dynamic_threshold) and (seq_mods[seq_id] + right_len <= seq_limits[seq_id])
                else:
                    cond_right = (right_len < dynamic_threshold) and (gap_count_right > dynamic_threshold)

                # APPLY REFINEMENT VETO (only for adaptive_3)
                if cond_right and orphan_method == "adaptive_3":
                    block_str_right = ''.join(seq_list[last_start:last_end])
                    if has_shared_homology(last_start, last_end, block_str_right, seq_id):
                        cond_right = False # Veto! Another terminal has this exact isolated block.

                if cond_right:
                    if orphan_action == "trim":
                        deleted = ''.join(seq_list[last_start:last_end])
                        seq_list[last_start:last_end] = ['-'] * right_len
                        if log_changes:
                            change_log.append(f"{seq_id}: Right block {last_start}-{last_end} deleted ('{deleted}')")
                    elif orphan_action == "push":
                        block_str = ''.join(seq_list[last_start:last_end])
                        seq_list[prev_end:last_end] = list(block_str) + ['-'] * gap_count_right
                        if log_changes:
                            change_log.append(f"{seq_id}: Right block {last_start}-{last_end} pushed to {prev_end}-{prev_end + right_len} ('{block_str}')")
                    
                    seq_mods[seq_id] += right_len
                    made_change_right = True
                    changed_this_seq = True

                if not made_change_left and not made_change_right:
                    break

            if changed_this_seq:
                alignment_changed = True
                alignment[seq_id] = ''.join(seq_list)

        if orphan_method in adaptive_methods:
            if alignment_changed:
                current_threshold = 1 # reset because merging might have created slightly larger new orphans
            else:
                current_threshold += 1
                if current_threshold <= orphan_threshold:
                    # check if all sequences hit their max modification limit
                    if all(seq_mods[sid] >= seq_limits[sid] for sid in alignment):
                        break

    if log_changes:
        return alignment, "\n".join(change_log)
    else:
        return alignment

def is_parsimony_non_informative(column):
    """
    Check if a column is non-informative.
    A column is non-informative if:
    - All characters are the same.
    - The characters include '?' and only one other character (e.g., 'A' and '?').

    Parameters:
        column (list): A list of characters in a column (e.g., ['A', 'A', 'A', '?']).
    
    Returns:
        bool: True if the column is non-informative, False if informative.
    """
    unique_characters = set(column)
        
    # Check if the column has only one unique character (parsimony non-informative)
    if len(unique_characters) == 1:
        return True
    
    # Check if the column contains '?' and exactly one other unique character
    if '?' in unique_characters and len(unique_characters) == 2:
        return True
    
    # Otherwise, the column is informative (multiple unique characters without '?')
    return False

def remove_non_informative_positions(alignment, removed_indices=None):
    """
    Remove parsimony non-informative positions from both the start and the end of the alignment.
    
    Parameters:
        alignment (dict): Dictionary where keys are sequence names and values are sequences (strings).
        removed_indices (list, optional): If provided, the function will append the indices of removed columns.
        
    Returns:
        dict: Updated alignment with non-informative positions removed.
    """
    # Convert the alignment to a list of sequences
    sequences = list(alignment.values())
    seq_length = len(sequences[0])

    # Remove non-informative positions from the start
    first_position = 0
    while first_position < seq_length and is_parsimony_non_informative([seq[first_position] for seq in sequences]):
        first_position += 1

    # Remove non-informative positions from the end
    last_position = seq_length - 1
    while last_position >= first_position and is_parsimony_non_informative([seq[last_position] for seq in sequences]):
        last_position -= 1

    # If logging, record the indices of removed columns
    if removed_indices is not None:
        removed_start = list(range(0, first_position))
        removed_end = list(range(last_position + 1, seq_length))
        removed_indices.extend(removed_start + removed_end)

    # Build the updated alignment
    for seq_name, seq in alignment.items():
        alignment[seq_name] = seq[first_position:last_position + 1]

    return alignment

##########################################################
# MAIN FUNCTIONS: STEP 3. IDENTIFICATION OF MISSING DATA #
##########################################################

def replace_terminal_gaps_dict(alignment):
    """
    Replaces terminal gaps ('-') with '?' in a sequence alignment stored as a dictionary,
    while keeping internal gaps as '-'.

    Parameters:
        alignment (dict): Dictionary with sequence names as keys and sequences as values.

    Returns:
        dict: Modified alignment with terminal gaps replaced by '?'.
    """
    # Rename 'modified_alignment' to 'alignment'
    for name, sequence in alignment.items():
        # Find the first and last non-gap characters
        first_non_gap = next((i for i, char in enumerate(sequence) if char != "-"), None)
        last_non_gap = next((i for i, char in enumerate(reversed(sequence), 1) if char != "-"), None)

        if first_non_gap is not None and last_non_gap is not None:
            last_non_gap = len(sequence) - last_non_gap  # Adjust reversed index
            
            # Replace terminal gaps with '?' and keep internal gaps as '-'
            new_sequence = (
                "?" * first_non_gap +
                sequence[first_non_gap:last_non_gap + 1] +
                "?" * (len(sequence) - last_non_gap - 1)
            )
            alignment[name] = new_sequence
        else:
            # Handle sequences with only gaps
            alignment[name] = "?" * len(sequence)
    
    return alignment

def replace_dashes_with_question_marks(alignment, internal_column_ranges=None, internal_leaves="all", internal_method="manual", internal_threshold=None):
    """
    Replace dashes with question marks in the specified column ranges for each sequence in the alignment.
    
    Args:
        alignment (dict): A dictionary with sequence IDs as keys and sequences as values.
        internal_column_ranges (list of tuples, optional): A list of tuples where each tuple defines a range of columns 
                                                          (inclusive) to check for dashes. E.g., [(5, 10), (50, 60)].
        internal_leaves (str or list, optional): If "all", replace dashes in all sequences. If a list of sequence IDs 
                                                   is provided, replace dashes in those sequences only.
        internal_method (str, optional): Defines how to replace dashes.
            - "manual": Specify column ranges and terminal sequences to replace dashes.
            - "semi": Replace internal blocks of contiguous dashes larger than the threshold with question marks.
        internal_threshold (int, optional): The threshold for "semi" method. Only internal blocks of contiguous gaps 
                                             larger than this threshold are replaced with question marks.
        
    Returns:
        dict: Updated alignment with dashes replaced by question marks in the specified columns.
    """
    
    # If internal_leaves is a list, only consider those sequences
    if internal_leaves != "all":
        sequences_to_process = set(internal_leaves)
    else:
        sequences_to_process = set(alignment.keys())
    
    # Convert the alignment into a list of sequences for easier indexing
    alignment = {seq_id: list(seq) for seq_id, seq in alignment.items()}  # Convert sequences to lists for mutability

    if internal_method == "manual":
        # Replace dashes with question marks in the specified column ranges
        for seq_id, seq in alignment.items():
            if seq_id not in sequences_to_process:
                continue  # Skip sequences not in the internal_leavesinternal_terminals list
            
            for start, end in internal_column_ranges:
                # Ensure the range is within the bounds of the sequence length
                start = max(0, start)
                end = min(len(seq), end)

                # Replace dashes with question marks within the specified range
                for i in range(start, end + 1):  # +1 because the end is inclusive
                    if seq[i] == '-':
                        seq[i] = '?'
        
    elif internal_method == "semi" and internal_threshold is not None:
        # Replace internal blocks of contiguous dashes larger than the internal_threshold with question marks
        for seq_id, seq in alignment.items():
            if seq_id not in sequences_to_process:
                continue  # Skip sequences not in the internal_leavesinternal_terminals list

            # Identify contiguous blocks of gaps (internal and terminal)
            gap_block_start = None
            for i in range(len(seq)):
                if seq[i] == '-':
                    if gap_block_start is None:
                        gap_block_start = i  # Start of a new gap block
                else:
                    if gap_block_start is not None:
                        # End of a gap block
                        gap_length = i - gap_block_start
                        if gap_length > internal_threshold and gap_block_start != 0 and gap_block_start != len(seq) - gap_length:
                            # It's an internal block larger than threshold, replace with '?'
                            for j in range(gap_block_start, i):
                                seq[j] = '?'
                        gap_block_start = None  # Reset for the next block
            # Check for a gap block at the end of the sequence
            if gap_block_start is not None:
                gap_length = len(seq) - gap_block_start
                if gap_length > internal_threshold and gap_block_start != 0:
                    for j in range(gap_block_start, len(seq)):
                        seq[j] = '?'
        
    # Convert the list back to a string
    alignment = {seq_id: ''.join(seq) for seq_id, seq in alignment.items()}

    return alignment

def n2question_func(alignment: dict, leaves='all', log=False):
    """
    Replace all ambiguous nucleotides 'N' or 'n' with '?' in selected sequences.

    Parameters:
    - alignment (dict): Dictionary where keys are sequence names and values are DNA sequences (str).
    - leaves (str or list): Sequence name(s) to apply the replacement. Use 'all' to apply to all sequences.
    - log (bool): If True, also return a log of replaced blocks with their positions.

    Returns:
    - dict: Modified alignment with 'N'/'n' replaced with '?'.
    - list (optional): List of tuples (seq_name, start, end) for each replaced block.
    """
    if leaves == 'all':
        leaves_to_process = alignment.keys()
    elif isinstance(leaves, str):
        leaves_to_process = [leaves]
    else:
        leaves_to_process = leaves

    updated_alignment = {}
    replacement_log = []

    for name, seq in alignment.items():
        if name in leaves_to_process:
            new_seq = []
            i = 0
            while i < len(seq):
                if seq[i] in ('N', 'n'):
                    start = i
                    while i < len(seq) and seq[i] in ('N', 'n'):
                        i += 1
                    end = i - 1
                    replacement_log.append((name, start, end))
                    new_seq.extend(['?'] * (end - start + 1))
                else:
                    new_seq.append(seq[i])
                    i += 1
            updated_alignment[name] = ''.join(new_seq)
        else:
            updated_alignment[name] = seq

    if log:
        return updated_alignment, replacement_log
    return updated_alignment

###################################################
# MAIN FUNCTIONS: STEP 4. SUCCESSIVE PARTITIONING #
###################################################

def add_breaks_terminal(alignment):
    """
    Add # in all instances of terminal gap opening/closure (indicated by ?).
       
    Parameters:
        alignment (dict): Dictionary where keys are sequence names and values are sequences.
    
    Returns:
        dict: Updated alignment with '#' before terminal gap opening and after terminal gap closure.
    """
    # Determine the length of the sequences
    seq_length = len(next(iter(alignment.values())))

    # Initialize a list to keep track of positions that need '#' in all sequences
    hash_positions = [False] * seq_length

    # Iterate through each sequence to find gap regions
    for seq in alignment.values():
        i = 0
        while i < seq_length:
            if seq[i] == '?':
                # Found the start of a gap region
                start = i
                while i < seq_length and seq[i] == '?':
                    i += 1
                end = i
                # Mark the positions for this gap region (avoid marking the first column)
                if start > 0:
                    hash_positions[start] = True
                if end < seq_length:
                    hash_positions[end] = True
            else:
                i += 1

    # Update each sequence with '#' at the identified positions
    for key in alignment:
        new_seq = []
        for i in range(seq_length):
            if hash_positions[i]:
                new_seq.append('#')
            new_seq.append(alignment[key][i])
        # Handle the case where the last position is a gap
        if hash_positions[-1]:
            new_seq.append('#')
        alignment[key] = ''.join(new_seq)

def classify_and_insert_hashtags(alignment, 
                                 partitioning_round=1, 
                                 partitioning_conservative="midpoint",
                                 log_csv_output=False, 
                                 csv_file_path="contiguous_invariant_blocks.csv"):
    # Step 1: Classify columns as invariant or variant
    num_sequences = len(alignment)
    num_columns = len(next(iter(alignment.values())))  # Get the number of columns from one sequence

    column_types = []  # To store the type of each column (invariant/variant)
    contiguous_invariant_blocks = []  # To store lengths and positions of invariant blocks

    for col_idx in range(num_columns):
        column = [seq[col_idx] for seq in alignment.values()]
        unique_values = set(column)
        if len(unique_values - {'?'}) == 1:
            column_types.append('invariant')
        else:
            column_types.append('variant')

    # Step 2: Identify contiguous invariant columns and their lengths
    current_invariant_block = None
    for col_idx in range(num_columns):
        if column_types[col_idx] == 'invariant':
            if current_invariant_block is None:
                current_invariant_block = {'start': col_idx, 'length': 1}
            else:
                current_invariant_block['length'] += 1
        else:
            if current_invariant_block:
                contiguous_invariant_blocks.append(current_invariant_block)
                current_invariant_block = None
    if current_invariant_block:
        contiguous_invariant_blocks.append(current_invariant_block)

    if partitioning_conservative not in {"midpoint", "flank"}:
        raise ValueError("partitioning_conservative must be 'midpoint' or 'flank'.")

    # Step 3: Log and optionally process blocks
    contiguous_invariant_blocks.sort(key=lambda x: x['length'], reverse=True)

    if partitioning_round == "max":
        add_breaks_terminal(alignment)  # Call the other function instead of inserting hashtags
    else:
        # Step 4: Track the positions for inserting hashtags
        hashtag_positions = []
        selected_blocks = contiguous_invariant_blocks[:partitioning_round]
        for block in selected_blocks:
            start_idx = block['start']
            end_idx = start_idx + block['length'] - 1

            if partitioning_conservative == "midpoint":
                hashtag_positions.append((start_idx + end_idx) // 2)
            else:
                hashtag_positions.extend([start_idx, end_idx + 1])

        # Step 5: Insert hashtags
        for seq_id, seq in alignment.items():
            sorted_positions = sorted(set(hashtag_positions))
            shift = 0
            for insert_idx in sorted_positions:
                adjusted_idx = insert_idx + shift
                seq = seq[:adjusted_idx] + '#' + seq[adjusted_idx:]
                shift += 1
            alignment[seq_id] = seq

    # Step 6: Optionally log to CSV
    if log_csv_output:
        file_dir = os.path.dirname(csv_file_path)
        if not os.path.exists(file_dir) and file_dir:
            os.makedirs(file_dir)

        with open(csv_file_path, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Start', 'End', 'Length'])
            for block in contiguous_invariant_blocks:
                start_idx = block['start']
                end_idx = start_idx + block['length'] - 1
                writer.writerow([start_idx, end_idx, block['length']])
        print(f"Log of contiguous invariant blocks written to {csv_file_path}")

    return alignment, contiguous_invariant_blocks

def refinement_question2hyphen(alignment):
    """
    Replaces contiguous '?' characters flanked by '#' with '-' in the sequence alignment.
    This includes blocks surrounded by '#' as well as the first and last blocks that are flanked only on one side.

    Parameters:
        alignment (dict): Dictionary with sequence names as keys and sequences as values.

    Returns:
        dict: Modified alignment with '?' replaced by '-' where flanked by '#'.
    """
    # Iterate through each sequence in the alignment
    for name, sequence in alignment.items():
        # Replace '?' surrounded by '#' on both sides (internal blocks)
        modified_sequence = re.sub(r'(?<=#)\?+(?=#)', lambda m: '-' * len(m.group(0)), sequence)
        
        # Replace '?' in the first block (flanked by # on the right)
        modified_sequence = re.sub(r'^(\?+)(?=#)', lambda m: '-' * len(m.group(0)), modified_sequence)
        
        # Replace '?' in the last block (flanked by # on the left)
        modified_sequence = re.sub(r'(?<=#)(\?+)$', lambda m: '-' * len(m.group(0)), modified_sequence)

        # Update the alignment dictionary with the modified sequence
        alignment[name] = modified_sequence

    return alignment

def remove_columns_with_W(alignment: dict) -> dict:
    """
    Remove problematic columns in a DNA alignment:
    
    1. For each sequence, identify 15 contiguous 'W' or 'w' characters. 
       Mark those columns for deletion.
    2. Remove columns where all values are:
       - only 'W'
       - only 'W' and '?'
       - only 'W', '?', and '-'
    
    Args:
        alignment (dict): Dictionary of {sequence_id: sequence_string}

    Returns:
        dict: Cleaned alignment with columns removed
    """
    import numpy as np

    # Convert alignment to matrix
    seq_ids = list(alignment.keys())
    seqs = list(alignment.values())
    alignment_array = np.array([list(seq) for seq in seqs])
    n_rows, n_cols = alignment_array.shape

    columns_to_remove = set()

    # Step 1: Find 15 contiguous W/w in any sequence
    for row in alignment_array:
        upper_row = [char.upper() for char in row]
        for i in range(n_cols - 14):
            if all(base == 'W' for base in upper_row[i:i+15]):
                columns_to_remove.update(range(i, i + 15))

    # Step 2: Remove columns with only W / W+? / W+?+-
    for col_idx in range(n_cols):
        col_bases = set(alignment_array[:, col_idx].astype(str).flatten().tolist())
        col_bases_upper = {b.upper() for b in col_bases}
        if col_bases_upper.issubset({'W'}) or \
           col_bases_upper.issubset({'W', '?'}) or \
           col_bases_upper.issubset({'W', '?', '-'}):
            columns_to_remove.add(col_idx)

    # Remove columns
    columns_to_keep = sorted(set(range(n_cols)) - columns_to_remove)
    cleaned_array = alignment_array[:, columns_to_keep]

    # Convert back to dictionary
    cleaned_alignment = {
        seq_id: ''.join(cleaned_array[i]) for i, seq_id in enumerate(seq_ids)
    }

    return cleaned_alignment

def remove_adjacent_pound_columns(alignment):
    """
    Remove adjacent columns of '#' in an alignment dictionary.
    Keeps only one '#' per contiguous block of # columns.

    Args:
        alignment (dict): DNA alignment {sequence_id: sequence_string}

    Returns:
        dict: Cleaned alignment with single '#' per contiguous block
    """
    # Transpose alignment to column-wise format
    columns = list(zip(*alignment.values()))

    # Identify positions to keep (i.e., not redundant '#')
    keep_indices = []
    prev_pound = False
    for i, col in enumerate(columns):
        if all(c == '#' for c in col):
            if not prev_pound:
                keep_indices.append(i)
                prev_pound = True
            # else skip this redundant pound column
        else:
            keep_indices.append(i)
            prev_pound = False

    # Rebuild alignment from kept columns
    cleaned_alignment = {}
    for seq_id in alignment:
        new_seq = ''.join([alignment[seq_id][i] for i in keep_indices])
        cleaned_alignment[seq_id] = new_seq

    return cleaned_alignment

def normalize_pound_columns(alignment):
    """
    Ensure that any alignment column containing at least one '#' becomes a full
    pound-sign column across all sequences.

    Args:
        alignment (dict): DNA alignment {sequence_id: sequence_string}

    Returns:
        dict: Alignment with normalized '#' columns.
    """
    if not alignment:
        return alignment

    seq_ids = list(alignment.keys())
    matrix = [list(alignment[seq_id]) for seq_id in seq_ids]
    num_cols = len(matrix[0])

    if any(len(row) != num_cols for row in matrix):
        raise ValueError("All sequences must have the same length.")

    for col_idx in range(num_cols):
        if any(row[col_idx] == '#' for row in matrix):
            for row in matrix:
                row[col_idx] = '#'

    return {seq_ids[i]: ''.join(matrix[i]) for i in range(len(seq_ids))}

def get_pound_column_indices(alignment):
    """
    Return the indices of true partition columns, defined as columns made only of '#'.

    Args:
        alignment (dict): DNA alignment {sequence_id: sequence_string}

    Returns:
        list[int]: Zero-based column indices in the current alignment.
    """
    if not alignment:
        return []

    columns = list(zip(*alignment.values()))
    return [i for i, col in enumerate(columns) if all(c == '#' for c in col)]

def get_pound_column_indices_from_output_file(output_path, output_format):
    """
    Read the saved alignment file and return the indices of '#' in the final output.

    Args:
        output_path (str): Path to the written alignment file.
        output_format (str): Biopython format used to read the alignment.

    Returns:
        list[int]: One-based column indices in the saved output alignment.
    """
    records = list(SeqIO.parse(output_path, output_format))
    if not records:
        return []

    first_sequence = str(records[0].seq)
    return [i + 1 for i, char in enumerate(first_sequence) if char == '#']

def get_partition_sizes_from_output_file(output_path, output_format):
    """
    Read the saved alignment file and return the sizes of partitions delimited by '#'.

    Args:
        output_path (str): Path to the written alignment file.
        output_format (str): Biopython format used to read the alignment.

    Returns:
        list[int]: Length of each partition in the saved output alignment.
    """
    records = list(SeqIO.parse(output_path, output_format))
    if not records:
        return []

    first_sequence = str(records[0].seq)
    return [len(part) for part in first_sequence.split('#')]

def equal_length_partitioning(alignment, partitioning_size=None, partitioning_round=None, log=False):
    """
    Insert pound signs '#' into a DNA alignment at regular intervals, either:
      - every `partitioning_size` bp (e.g. every 100 bp), or
      - using `partitioning_round` to divide into equal-length partitions.

    Args:
        alignment (dict): Dictionary of sequences {seq_id: sequence_str}
        partitioning_size (int): Length of partitions (in base pairs)
        partitioning_round (int): Number of pound signs to insert (creates N+1 partitions)
        log (bool): Whether to print insertion positions and runtime

    Returns:
        dict: Modified alignment with pound signs inserted
    """
    start_time = time.time()
    
    if not alignment:
        raise ValueError("Alignment is empty")

    aln_length = len(next(iter(alignment.values())))
    positions = []

    # Mode 1: divide by fixed size
    if partitioning_size:
        positions = list(range(partitioning_size, aln_length, partitioning_size))

    # Mode 2: divide by number of equal partitions
    elif partitioning_round:
        if partitioning_round >= aln_length:
            raise ValueError("Too many partitions requested for alignment length.")
        interval = aln_length // (partitioning_round + 1)
        positions = [(i + 1) * interval for i in range(partitioning_round)]

    else:
        raise ValueError("You must provide either 'partitioning_size' or 'partitioning_round'")

    # Insert '#' at the selected positions
    updated_alignment = {}
    for seq_id, seq in alignment.items():
        seq_list = list(seq)
        for offset, pos in enumerate(positions):
            seq_list.insert(pos + offset, "#")  # account for shifting due to insertions
        updated_alignment[seq_id] = ''.join(seq_list)

    runtime = time.time() - start_time

    if log:
        print("--- equal_length_partitioning ---")
        print(f"Total time: {runtime:.4f} seconds")
        print(f"Pound signs inserted at columns: {positions}\n")

    return updated_alignment

def split_alignment_into_equal_chunks(alignment, chunk_size):
    """
    Split an alignment dictionary into consecutive equal-length chunks.

    Args:
        alignment (dict): Dictionary of sequences {seq_id: sequence_str}
        chunk_size (int): Maximum size of each chunk

    Returns:
        list[dict]: List of alignment chunks
    """
    if not alignment:
        return []
    if chunk_size is None or chunk_size <= 0:
        raise ValueError("chunk_size must be a positive integer.")

    aln_length = len(next(iter(alignment.values())))
    return [
        {seq_id: seq[start:start + chunk_size] for seq_id, seq in alignment.items()}
        for start in range(0, aln_length, chunk_size)
    ]

def join_alignment_chunks_with_pounds(chunks, return_metadata=False):
    """
    Join alignment chunks, separating them with full '#' columns.

    Args:
        chunks (list[dict]): List of alignment chunks
        return_metadata (bool): If True, also return the origin of inserted '#'
                                columns in final one-based coordinates.

    Returns:
        dict or tuple: Joined alignment, optionally with metadata
    """
    if not chunks:
        return ({}, {"max_size_positions": [], "method_positions": []}) if return_metadata else {}

    seq_ids = list(chunks[0].keys())
    joined = {seq_id: "" for seq_id in seq_ids}
    max_size_positions = []
    method_positions = []
    offset = 0

    for chunk_idx, chunk in enumerate(chunks):
        method_positions.extend([offset + idx + 1 for idx in get_pound_column_indices(chunk)])
        for seq_id in seq_ids:
            joined[seq_id] += chunk[seq_id]

        chunk_length = len(next(iter(chunk.values()))) if chunk else 0
        if chunk_idx < len(chunks) - 1:
            max_size_positions.append(offset + chunk_length + 1)
            for seq_id in seq_ids:
                joined[seq_id] += "#"
            offset += chunk_length + 1
        else:
            offset += chunk_length

    metadata = {
        "max_size_positions": max_size_positions,
        "method_positions": method_positions,
    }
    return (joined, metadata) if return_metadata else joined

def apply_partitioning_strategy(alignment,
                                partitioning_method,
                                partitioning_round=0,
                                partitioning_conservative="midpoint",
                                partitioning_size=None):
    """
    Apply a single partitioning strategy to one alignment dictionary.

    Returns:
        tuple: (updated_alignment, conservative_blocks, partitioning_log_entry, balanced_metadata)
    """
    partitioning_log_entry = None
    conservative_blocks = []
    balanced_metadata = None

    if partitioning_method == "conservative" and partitioning_round > 0:
        alignment, conservative_blocks = classify_and_insert_hashtags(
            alignment,
            partitioning_round=partitioning_round,
            partitioning_conservative=partitioning_conservative
        )

    elif partitioning_method == "equal":
        if partitioning_round > 0:
            alignment = equal_length_partitioning(
                alignment=alignment,
                partitioning_round=partitioning_round,
                partitioning_size=None,
                log=False
            )
        elif partitioning_size:
            alignment = equal_length_partitioning(
                alignment=alignment,
                partitioning_round=None,
                partitioning_size=partitioning_size,
                log=False
            )

    elif partitioning_method == "max":
        alignment = insert_pound_around_questions(alignment)

    elif partitioning_method == "balanced":
        if partitioning_round > 0:
            try:
                alignment, balanced_metadata = balanced_partitioning(
                    alignment,
                    log=False,
                    partitioning_round=partitioning_round,
                    return_metadata=True
                )
            except (IndexError, ValueError):
                partitioning_log_entry = (
                    f"WARNING: 'balanced' partitioning with partitioning_round={partitioning_round} was skipped for this alignment. "
                    "This typically happens when the alignment has fewer blocks of missing data ('?') than required. "
                    "The process will continue without partitioning this file."
                )

    return alignment, conservative_blocks, partitioning_log_entry, balanced_metadata

def insert_pound_around_questions(alignment):
    """
    Insert '#' columns around the opening and closure of '?' blocks.
    
    - If a '?' block starts in the first column, insert '#' in the second column.
    - If a '?' block ends in the last column, insert '#' in the penultimate column.
    - Elsewhere, insert '#' before or after transition points.
    
    Parameters:
        alignment (dict): {sequence_name: aligned_sequence}
    
    Returns:
        dict: Updated alignment with '#' inserted as new columns.
    """
    if not alignment:
        return alignment

    names = list(alignment.keys())
    matrix = [list(alignment[name]) for name in names]
    num_rows = len(matrix)
    num_cols = len(matrix[0])

    # Check consistency
    if any(len(row) != num_cols for row in matrix):
        raise ValueError("All sequences must have the same length.")

    insert_indices = set()

    for col in range(num_cols):
        column = [matrix[r][col] for r in range(num_rows)]

        if '?' not in column:
            continue

        # Check for start of '?' block
        if col == 0:
            next_col = [matrix[r][col + 1] for r in range(num_rows)]
            if any(c == '?' and nc != '?' for c, nc in zip(column, next_col)) or \
               any(c != '?' and nc == '?' for c, nc in zip(column, next_col)):
                insert_indices.add(1)
        elif col == num_cols - 1:
            prev_col = [matrix[r][col - 1] for r in range(num_rows)]
            if any(pc != '?' and c == '?' for pc, c in zip(prev_col, column)) or \
               any(pc == '?' and c != '?' for pc, c in zip(prev_col, column)):
                insert_indices.add(num_cols - 1)
        else:
            prev_col = [matrix[r][col - 1] for r in range(num_rows)]
            next_col = [matrix[r][col + 1] for r in range(num_rows)]

            # Transition into '?' block → insert '#' before
            if any(pc != '?' and c == '?' for pc, c in zip(prev_col, column)):
                insert_indices.add(col)

            # Transition out of '?' block → insert '#' after
            if any(c == '?' and nc != '?' for c, nc in zip(column, next_col)):
                insert_indices.add(col + 1)

    # Insert pound signs (from right to left to maintain positions)
    for idx in sorted(insert_indices, reverse=True):
        for row in matrix:
            row.insert(idx, '#')

    # Return as dict
    return {names[i]: ''.join(row) for i, row in enumerate(matrix)}

def balanced_partitioning(alignment, log=False, partitioning_round=1, return_metadata=False):
    """
    Partition alignment by merging blocks flanked by '#' based on a length threshold.

    Steps:
    1. Call insert_pound_around_questions to place '#' around blocks of '?'.
    2. Identify blocks flanked by '#'.
    3. Use the N-th longest block (based on partitioning_round) as the merging threshold.
    4. Merge adjacent blocks (delete '#' between them) if their combined size < N.
    5. First pass left-to-right, then right-to-left.

    Parameters:
        alignment (dict): {sequence_name: aligned_sequence}
        log (bool): If True, print logs of changes.
        partitioning_round (int): Use the N-th longest block as the threshold.

    Returns:
        dict or tuple: Updated alignment with balanced '#' partitioning, optionally
                       with metadata about the reference missing-data block.
    """
    from copy import deepcopy

    alignment = insert_pound_around_questions(alignment)
    names = list(alignment.keys())
    matrix = [list(alignment[name]) for name in names]
    num_cols = len(matrix[0])

    # Identify block boundaries (flanked by '#')
    pound_indices = [i for i in range(num_cols) if all(row[i] == '#' for row in matrix)]
    pound_indices = [-1] + pound_indices + [num_cols]  # Add virtual boundaries
    blocks = [(pound_indices[i] + 1, pound_indices[i + 1]) for i in range(len(pound_indices) - 1)]

    # Compute lengths and positions of all blocks
    block_lengths = [(i, end - start, start, end - 1) for i, (start, end) in enumerate(blocks)]  # (index, size, start, end)
    sorted_blocks = sorted(block_lengths, key=lambda x: x[1], reverse=True)

    if partitioning_round > len(sorted_blocks):
        raise ValueError(f"partitioning_round={partitioning_round} is greater than total number of blocks.")

    block_idx, N, start_col, end_col = sorted_blocks[partitioning_round - 1]

    merged_log = []
    balanced_metadata = {
        "reference_block_rank": partitioning_round,
        "reference_block_size": N,
        "reference_block_start": start_col,
        "reference_block_end": end_col,
    }

    if log:
        suffix = {1: "st", 2: "nd", 3: "rd"}.get(partitioning_round, "th")
        merged_log.append(
            f"Threshold = size of the {partitioning_round}{suffix} largest block: {N} columns "
            f"(columns {start_col} to {end_col})"
        )

    def merge_blocks(direction):
        nonlocal matrix
        while True:
            # Recompute pound indices and blocks after each merge
            num_cols = len(matrix[0])
            pound_indices = [i for i in range(num_cols) if all(row[i] == '#' for row in matrix)]
            pound_indices = [-1] + pound_indices + [num_cols]
            blocks = [(pound_indices[i] + 1, pound_indices[i + 1]) for i in range(len(pound_indices) - 1)]
            changed = False

            if direction == "left":
                i_range = range(0, len(blocks) - 1)
            else:
                i_range = range(len(blocks) - 2, -1, -1)

            for i in i_range:
                left_start, left_end = blocks[i]
                right_start, right_end = blocks[i + 1]
                combined_size = (left_end - left_start) + (right_end - right_start)

                if combined_size < N:
                    # Remove pound between blocks
                    pound_col_index = left_end
                    for row in matrix:
                        del row[pound_col_index]
                    merged_log.append(f"Merged blocks at columns {left_start}-{left_end - 1} and {right_start}-{right_end - 1} "
                                      f"(deleted '#' at col {pound_col_index}) → size={combined_size}")
                    changed = True
                    break  # Restart loop

            if not changed:
                break

    merge_blocks("left")
    merge_blocks("right")

    if log:
        print("\n".join(merged_log))

    result_alignment = {names[i]: ''.join(row) for i, row in enumerate(matrix)}
    if return_metadata:
        return result_alignment, balanced_metadata
    return result_alignment

##############################
# AUXILIARY FUNCTIONS TO LOG #
##############################

def compute_summary_after(alignment):
    num_seqs = len(alignment)

    # Transpose to columns
    columns = list(zip(*alignment.values()))
    
    # Count columns that contain pound signs
    total_pound = len(get_pound_column_indices(alignment))
    
    # Alignment length excluding columns of only pound signs
    aln_length = len(columns)


    # Count nucleotide and gap characters
    total_nt = sum(c in "ACGTacgt" for seq in alignment.values() for c in seq)
    total_gaps = sum(seq.count("-") for seq in alignment.values())
    total_ns = sum(c in "Nn" for seq in alignment.values() for c in seq)
    total_qm = sum(seq.count("?") for seq in alignment.values())  # <- includes everything now

    # Count additional missing data as gap-only blocks between #
    missing_by_partition = 0
    for seq in alignment.values():
        parts = seq.split("#")
        for part in parts:
            if all(c == '-' for c in part):
                missing_by_partition += len(part)

    total_missing = total_qm + missing_by_partition

    return {
        "num_seqs": num_seqs,
        "aln_length": aln_length,
        "total_nt": total_nt,
        "total_gaps": total_gaps,
        "total_ns": total_ns,
        "total_qm": total_qm,
        "total_pound": total_pound,
        "missing_by_partition": missing_by_partition,
        "total_missing": total_missing
    }

def detect_fully_missing_partitions(alignment):
    """
    Logs all `?` and `-` from fully missing partitions.
    A fully missing partition is a region (between #) in which a sequence has only dashes.
    """
    log_entries = []
    total_question_marks = 0
    total_dash_from_missing_partitions = 0

    for seq_id, seq in alignment.items():
        total_question_marks += seq.count("?")

        parts = seq.split("#")
        col_index = 0  # absolute position tracker
        for i, part in enumerate(parts):
            if all(c == '-' for c in part):
                dash_count = len(part)
                total_dash_from_missing_partitions += dash_count
                start = col_index
                end = col_index + dash_count - 1
                log_entries.append(f"{seq_id}: partition {i} ({start}–{end}, length {dash_count}) fully missing (all '-')")
            col_index += len(part) + 1  # +1 for the '#' removed in split

    summary_lines = [
        f"Total '?' characters: {total_question_marks}",
        f"Total '-' characters in fully missing partitions: {total_dash_from_missing_partitions}",
        f"Combined total: {total_question_marks + total_dash_from_missing_partitions}"
    ]

    if log_entries:
        return "\n".join(summary_lines + ["", "List of zero-length partitions:"] + log_entries)
    return "\n".join(summary_lines)

#######################
# INTEGRATED FUNCTIONS #
#######################

import subprocess
import tempfile
import re
import os
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def UP2AP(input_fasta, output_fasta, keep_unusual=False):
    """
    1. Replaces '#' with 15 'w's.
    2. Aligns using MAFFT.
    3. Replaces aligned 'w' blocks (contiguous or gapped) back to '#'.
    4. Rectifies alignment: pads partitions to ensure '#' are vertically aligned.
    5. Cleans alignment: removes columns that contain ONLY gaps (-).
    6. If keep_unusual=False, removes '#' columns and replaces any remaining '?' with '-'.
    """
    
    # Configuration
    PLACEHOLDER_CHAR = 'w'
    REPEAT_COUNT = 15
    PLACEHOLDER_SEQ = PLACEHOLDER_CHAR * REPEAT_COUNT
    
    # Check if MAFFT is installed
    if shutil.which("mafft") is None:
        raise EnvironmentError("MAFFT not found. Please install MAFFT and ensure it is in your PATH.")

    # Temp files
    temp_input = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta")
    temp_input_name = temp_input.name
    temp_output = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta")
    temp_output_name = temp_output.name
    temp_output.close()

    try:
        # --- STEP 1: Pre-processing (Replace # with www...) ---
        original_records = list(SeqIO.parse(input_fasta, "fasta"))
        if not original_records:
            raise ValueError("Input file is empty or not a valid FASTA.")

        modified_records = []
        for record in original_records:
            seq_str = str(record.seq)
            if "#" in seq_str:
                seq_str = seq_str.replace("#", PLACEHOLDER_SEQ)
            
            new_record = SeqRecord(
                Seq(seq_str),
                id=record.id,
                description=record.description
            )
            modified_records.append(new_record)
        
        SeqIO.write(modified_records, temp_input, "fasta")
        temp_input.close()

        # --- STEP 2: Alignment (Run MAFFT) ---
        cmd = ["mafft", "--auto", "--quiet"]
        if keep_unusual:
            cmd.append("--anysymbol")
        cmd.append(temp_input_name)

        with open(temp_output_name, "w") as out_handle:
            process = subprocess.run(cmd, stdout=out_handle, stderr=subprocess.PIPE, text=True)
            
        if process.returncode != 0:
            raise RuntimeError(f"MAFFT failed: {process.stderr}")

        # --- STEP 3: Restore '#' (Regex substitution) ---
        # Matches 'w' followed by exactly 14 instances of (optional hyphens + 'w')
        pattern = re.compile(f"{PLACEHOLDER_CHAR}(?:-*{PLACEHOLDER_CHAR}){{{REPEAT_COUNT - 1}}}", re.IGNORECASE)

        aligned_records = list(SeqIO.parse(temp_output_name, "fasta"))
        restored_seqs = []
        
        for record in aligned_records:
            seq_str = str(record.seq)
            new_seq_str = pattern.sub("#", seq_str)
            restored_seqs.append(new_seq_str)

        # --- STEP 4: Rectify Alignment (Force Perfect Columns for #) ---
        split_seqs = [s.split('#') for s in restored_seqs]
        
        # Check consistency
        num_partitions = len(split_seqs[0])
        is_consistent = all(len(parts) == num_partitions for parts in split_seqs)

        if not is_consistent:
            print("Warning: Mismatch in '#' count. Skipping rectification step.")
            rectified_seqs = restored_seqs
        else:
            # Transpose to iterate by partition group
            columns = zip(*split_seqs)
            normalized_columns = []
            
            for col_group in columns:
                # Pad to max length in this partition
                max_len = max(len(segment) for segment in col_group)
                padded_group = [segment.ljust(max_len, '-') for segment in col_group]
                normalized_columns.append(padded_group)
            
            # Transpose back and join
            rectified_seqs = ["#".join(parts) for parts in zip(*normalized_columns)]

        # --- STEP 5: Remove Invariant Gap Columns ---
        # rectified_seqs is a list of strings of equal length.
        # zip(*rectified_seqs) creates a tuple for every column in the alignment.
        
        if not rectified_seqs:
            final_seq_strs = []
        else:
            # Keep the column ONLY IF it is not composed entirely of '-'
            valid_columns = [col for col in zip(*rectified_seqs) if not all(c == '-' for c in col)]
            
            if not keep_unusual:
                # '#' forms full columns due to Step 4, so dropping columns with '#' perfectly removes them.
                # '?' are normally ignored and removed by MAFFT without --anysymbol. 
                # If any survive, replace them safely with '-' to protect alignment length.
                valid_columns = [col for col in valid_columns if '#' not in col]
                valid_columns = [tuple('-' if c == '?' else c for c in col) for col in valid_columns]
            
            # zip(*valid_columns) transposes back to rows (sequences)
            if valid_columns:
                final_seq_strs = ["".join(row) for row in zip(*valid_columns)]
            else:
                # In the rare case that the entire alignment is gaps
                final_seq_strs = [""] * len(rectified_seqs)

        # --- STEP 6: Write Output ---
        final_records = []
        for i, record in enumerate(aligned_records):
            final_record = SeqRecord(
                Seq(final_seq_strs[i]),
                id=record.id,
                description=record.description
            )
            final_records.append(final_record)

        SeqIO.write(final_records, output_fasta, "fasta")
        print(f"Process complete. {len(final_records)} sequences processed.")
        print(f"Output saved to: {output_fasta}")

    finally:
        if os.path.exists(temp_input_name):
            os.remove(temp_input_name)
        if os.path.exists(temp_output_name):
            os.remove(temp_output_name)

def GB2MSA(input_file, 
           output_prefix, 
           delimiter=',', 
           write_names=False, 
           log=False, 
           orphan_threshold=10,
           multi_amplicon_min_overlap=10,
           multi_amplicon_mismatch_rate=0.05,
           multi_amplicon_action="trim",
           aligner="mafft",
           return_gene_logs=False,
           write_run_log=True):
    """
    Complete GenBank-to-MSA pipeline:
    1. Downloads sequences from GenBank and aligns them by gene using the selected aligner.
    2. Cleans the alignments by replacing internal missing data and removing empty columns.
    3. Applies GB2MSA_3 to replace selected gap and orphan nucleotide blocks with '?'.
    4. Applies GB2MSA_4 to replace blocks of 15 or more w/W (with or without interspersed gaps) with '?'.
    5. Replaces terminal '?' in sequences with '-'.
    6. Deletes intermediate files ending with '_aligned.fasta'.
    7. Optionally logs wall clock and CPU time to a log file named '<output_prefix>_log.txt'.
    """
    start_wall = time.time()
    start_cpu = time.process_time()

    # Step 1: Generate aligned FASTA files
    aligned_files, gb1_logs, gene_names = GB2MSA_1(
        input_file,
        output_prefix,
        delimiter=delimiter,
        write_names=write_names,
        multi_amplicon_min_overlap=multi_amplicon_min_overlap,
        multi_amplicon_mismatch_rate=multi_amplicon_mismatch_rate,
        multi_amplicon_action=multi_amplicon_action,
        aligner=aligner
    )

    # Step 2: Clean each aligned FASTA file
    cleaned_files = []
    for aligned_file in aligned_files:
        cleaned_file = GB2MSA_2(aligned_file)
        cleaned_files.append(cleaned_file)

    # Step 3: Apply GB2MSA_3 to clean orphan gap/nucleotide blocks
    for cleaned_file in cleaned_files:
        records = list(SeqIO.parse(cleaned_file, "fasta"))
        alignment_dict = {record.id: str(record.seq) for record in records}
        updated_dict = GB2MSA_3(alignment_dict, orphan_threshold=orphan_threshold, log=log)
        
        # Step 4: Apply GB2MSA_4 to handle w/W blocks
        updated_dict = GB2MSA_4(updated_dict)

        updated_records = []
        for record in records:
            record.seq = Seq(updated_dict[record.id])
            updated_records.append(record)
        with open(cleaned_file, "w") as out_handle:
            SeqIO.write(updated_records, out_handle, "fasta")

    # Step 5: Replace terminal '?' with '-' in each sequence
    for cleaned_file in cleaned_files:
        records = list(SeqIO.parse(cleaned_file, "fasta"))
        updated_records = []
        for record in records:
            seq = str(record.seq)
            left = len(seq) - len(seq.lstrip('?'))
            right = len(seq) - len(seq.rstrip('?'))
            new_seq = '-' * left + seq[left:len(seq)-right] + '-' * right if right > 0 else '-' * left + seq[left:]
            record.seq = Seq(new_seq)
            updated_records.append(record)
        with open(cleaned_file, "w") as out_handle:
            SeqIO.write(updated_records, out_handle, "fasta")

    # Step 6: Delete intermediate *_aligned.fasta files
    for aligned_file in aligned_files:
        if aligned_file.endswith("_aligned.fasta") and os.path.exists(aligned_file):
            os.remove(aligned_file)

    # Step 7: Log timing and amplicons if requested
    if log and write_run_log:
        end_wall = time.time()
        end_cpu = time.process_time()
        wall_time = end_wall - start_wall
        cpu_time = end_cpu - start_cpu

        log_file = f"{output_prefix}_log.txt"
        with open(log_file, "a") as lf:
            lf.write(f"--- GB2MSA run for '{output_prefix}' ---\n")
            
            # Log amplicon merges and ? insertions per gene
            for gene_name, cleaned_file in zip(gene_names, cleaned_files):
                lf.write(f"\n--- Locus: {gene_name} ---\n")
                gene_logs = gb1_logs.get(gene_name, {})
                
                # Scan final files to determine where '?' ended up natively in the alignment
                records = list(SeqIO.parse(cleaned_file, "fasta"))
                for record in records:
                    seq = str(record.seq)
                    q_blocks = []
                    for m in re.finditer(r'\?+', seq):
                        q_blocks.append((m.start() + 1, m.end()))  # 1-based index logging
                    
                    metadata = gene_logs.get(record.id)
                    events = []
                    if metadata:
                        events.extend(
                            _build_multi_amplicon_log(
                                metadata["requested_sources"],
                                metadata["loaded_sources"],
                                metadata["failed_sources"],
                                metadata["merge_events"],
                                metadata["final_groups"]
                            )
                        )
                    # Append alignment coordinate logs
                    for start, end in q_blocks:
                        events.append(f"Internal missing data (?) at alignment columns {start}-{end} (length: {end-start+1})")
                        
                    if events:
                        lf.write(f"{record.id}:\n")
                        for ev in events:
                            lf.write(f"  - {ev}\n")
                            
            lf.write(f"\nWall clock time: {wall_time:.2f} seconds\n")
            lf.write(f"CPU time: {cpu_time:.2f} seconds\n\n")

    if return_gene_logs:
        return cleaned_files, gb1_logs, gene_names
    return cleaned_files

def addSeq(
    alignment,
    new_seqs,
    output,
    write_names=True,
    orphan_threshold=0,
    log=False,
    n2question=None,
    gaps2question=None
):
    """
    Add new sequences to an existing alignment using MAFFT, clean and standardize the result.

    Steps performed:
    1. Remove '#' columns from original alignment.
    2. Align new sequences with MAFFT using --add.
    3. Trim orphan nucleotide blocks from new sequences.
    4. Remove short DNA blocks near # in first and last partitions.
    5. Replace terminal '-' with '?'.
    6. Optionally replace all N/n with '?' in selected sequences.
    7. Optionally replace long gap blocks with '?'.
    8. Reinsert '#' columns.
    9. Replace all-'?' blocks between '#' with '-'.
    10. Write the result and an optional log file.

    Parameters:
        alignment (str or dict): Existing alignment file path (FASTA) or dictionary {id: sequence}.
        new_seqs (str or dict): New sequences to add (FASTA path or dict {id: sequence}).
        output (str): Path to write the updated alignment in FASTA format.
        write_names (bool): Whether to write a _terminal_names.txt file listing sequence IDs.
        orphan_threshold (int): Threshold to detect and remove orphan DNA blocks.
        log (bool): If True, writes a log file with trimming and runtime information.
        n2question (str, list or None): Replace 'N/n' with '?' in specific sequences:
            - 'all': apply to all sequences
            - str: apply to a single sequence ID
            - list: apply to listed sequence IDs
        gaps2question (int or None): Replace contiguous gap blocks larger than this threshold with '?'. Only applied to added sequences.

    Returns:
        None

    Example usage:
        addSeq("alignment.fasta", "new.fasta", "updated.fasta", n2question="seq123", log=True)
        addSeq(alignment_dict, new_dict, "out.fas", n2question='all')
    """

    # Start timing the execution
    start_time = time.time()
    temp_files_to_remove = []  # Temporary files to be removed after execution
    log_lines = [] 

    # Log the function call and parameters for reproducibility
    if log:
        alignment_repr = alignment if isinstance(alignment, str) else "<dict>"
        new_seqs_repr = new_seqs if isinstance(new_seqs, str) else "<dict>"
        cmd_used = (
            f"addSeq(alignment={alignment_repr}, new_seqs={new_seqs_repr}, "
            f"output='{output}', write_names={write_names}, "
            f"orphan_threshold={orphan_threshold}, log={log}, "
            f"n2question={n2question}, gaps2question={gaps2question})"
        )
        log_lines.append(f"Command used: {cmd_used}")
        log_lines.append("")

    # === Step 1: Load and clean the alignment ===
    def write_dict_to_temp_fasta(seq_dict):
        records = [SeqRecord(Seq(seq), id=str(seq_id), description="") for seq_id, seq in seq_dict.items()]
        tmp = tempfile.NamedTemporaryFile("w+", delete=False)
        SeqIO.write(records, tmp, "fasta")
        tmp.close()
        return tmp.name

    if isinstance(alignment, dict):
        alignment_path = write_dict_to_temp_fasta(alignment)
        temp_files_to_remove.append(alignment_path)
    elif isinstance(alignment, str):
        alignment_path = alignment
    else:
        raise ValueError("alignment must be a FASTA file path or a dictionary")

    records = list(SeqIO.parse(alignment_path, "fasta"))
    if not records:
        raise ValueError("Input alignment is empty or not found")

    aln_len = len(records[0].seq)
    # Identify columns that are '#' characters to temporarily remove them for alignment
    pound_cols = [i for i in range(aln_len) if any(rec.seq[i] == '#' for rec in records)]

    # Log input alignment info
    if log:
        log_lines.append(f"Input alignment: {len(records)} sequences")
        log_lines.append(f"Input alignment: {len(pound_cols)} # columns")

    def remove_cols(seq, cols):
        return ''.join(seq[i] for i in range(len(seq)) if i not in cols)

    # Remove '#' columns
    aln_no_pound = [
        SeqRecord(Seq(remove_cols(str(rec.seq), pound_cols)), id=rec.id, description="")
        for rec in records
    ]
    # Write cleaned alignment to a temporary file
    with tempfile.NamedTemporaryFile("w+", delete=False) as aln_tmp:
        SeqIO.write(aln_no_pound, aln_tmp, "fasta")
        aln_path = aln_tmp.name
        temp_files_to_remove.append(aln_path)

    # === Step 2: Load new sequences ===
    if isinstance(new_seqs, dict):
        new_seqs_path = write_dict_to_temp_fasta(new_seqs)
        new_seq_count = len(new_seqs)
        temp_files_to_remove.append(new_seqs_path)
    elif isinstance(new_seqs, str):
        new_seq_count = sum(1 for _ in SeqIO.parse(new_seqs, "fasta"))
        new_seqs_path = new_seqs
    else:
        raise ValueError("new_seqs must be a FASTA file path or a dictionary")
    # Log new sequence info
    if log:
        log_lines.append(f"Input new sequences: {new_seq_count} sequences")

    # === Step 3: Align new sequences with MAFFT ===
    with tempfile.NamedTemporaryFile("w+", delete=False) as out_tmp:
        out_path = out_tmp.name
        temp_files_to_remove.append(out_path)

    try:
        subprocess.run(
            ['mafft', '--add', new_seqs_path, '--keeplength', '--preservecase', aln_path],
            check=True,
            stdout=open(out_path, 'w'),
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        for file in temp_files_to_remove:
            try:
                os.remove(file)
            except Exception:
                pass
        raise RuntimeError(f"MAFFT failed!\nCommand: {e.cmd}\nExit status: {e.returncode}\nMAFFT error output:\n{e.stderr}")

    mafft_aligned_records = list(SeqIO.parse(out_path, "fasta"))
    original_ids = {rec.id for rec in records}
    new_records = [rec for rec in mafft_aligned_records if rec.id not in original_ids]


    def replace_gap_blocks(seq, threshold, seq_id=None):
        seq_list = list(seq)
        replaced_log = []
        i = 0
        while i < len(seq_list):
            if seq_list[i] == '-':
                start = i
                while i < len(seq_list) and seq_list[i] == '-':
                    i += 1
                if (i - start) > threshold:
                    for j in range(start, i):
                        seq_list[j] = '?'
                    if seq_id:
                        replaced_log.append(f"{seq_id}: {i - start} contiguous '-' replaced with '?' at {start}–{i}")
            else:
                i += 1
        return ''.join(seq_list), replaced_log

    def find_dna_blocks(seq_list, start, end):
        blocks = []
        i = start
        while i < end:
            if seq_list[i] not in "-?#":
                s = i
                while i < end and seq_list[i] not in "-?#":
                    i += 1
                e = i
                blocks.append((s, e))
            else:
                i += 1
        return blocks

    def trim_orphan_blocks(seq, threshold, seq_id=None):
        seq_list = list(seq)
        trimmed_log = []

        def find_blocks(seq_list):
            blocks = []
            i = 0
            while i < len(seq_list):
                if seq_list[i] not in "-?#":
                    start = i
                    while i < len(seq_list) and seq_list[i] not in "-?#":
                        i += 1
                    end = i
                    blocks.append((start, end))
                else:
                    i += 1
            return blocks

        changed = True
        while changed:
            changed = False

            while True:
                blocks = find_blocks(seq_list)
                if len(blocks) < 2:
                    break
                first_start, first_end = blocks[0]
                next_start = blocks[1][0]
                gap_count = seq_list[first_end:next_start].count('-') + seq_list[first_end:next_start].count('?')
                size = first_end - first_start
                if size < threshold and gap_count > threshold:
                    deleted = ''.join(seq_list[first_start:first_end])
                    seq_list[first_start:first_end] = ['-'] * size
                    if seq_id:
                        trimmed_log.append(f"{seq_id}: Left {first_start}–{first_end} (size={size}, '{deleted}')")
                    changed = True
                    continue
                break

            while True:
                blocks = find_blocks(seq_list)
                if len(blocks) < 2:
                    break
                last_start, last_end = blocks[-1]
                prev_end = blocks[-2][1]
                gap_count = seq_list[prev_end:last_start].count('-') + seq_list[prev_end:last_start].count('?')
                size = last_end - last_start
                if size < threshold and gap_count > threshold:
                    deleted = ''.join(seq_list[last_start:last_end])
                    seq_list[last_start:last_end] = ['-'] * size
                    if seq_id:
                        trimmed_log.append(f"{seq_id}: Right {last_start}–{last_end} (size={size}, '{deleted}')")
                    changed = True
                    continue
                break

        if pound_cols:
            first_hash = min(pound_cols)
            last_hash = max(pound_cols)
            blocks = find_dna_blocks(seq_list, 0, first_hash)
            if len(blocks) == 1:
                s, e = blocks[0]
                if e == first_hash and (e - s) < threshold:
                    deleted = ''.join(seq_list[s:e])
                    seq_list[s:e] = ['-'] * (e - s)
                    if seq_id:
                        trimmed_log.append(f"{seq_id}: Left {s}–{e} (size={e - s}, '{deleted}')")
            blocks = find_dna_blocks(seq_list, last_hash + 1, len(seq_list))
            if len(blocks) == 1:
                s, e = blocks[0]
                if s == last_hash + 1 and (e - s) < threshold:
                    deleted = ''.join(seq_list[s:e])
                    seq_list[s:e] = ['-'] * (e - s)
                    if seq_id:
                        trimmed_log.append(f"{seq_id}: Right {s}–{e} (size={e - s}, '{deleted}')")

        return ''.join(seq_list), trimmed_log

    trimmed_new_records = []
    all_trim_logs = []
    gaps2q_log = []  # This will no longer be used for logging replaced gaps

    for rec in new_records:
        trimmed_seq, seq_log = trim_orphan_blocks(str(rec.seq), orphan_threshold, seq_id=rec.id)
        # Remove gaps2question here to not log replaced gaps multiple times
        trimmed_new_records.append(SeqRecord(Seq(trimmed_seq), id=rec.id, description=""))
        all_trim_logs.extend(seq_log)

    processed_records = [rec for rec in mafft_aligned_records if rec.id in original_ids] + trimmed_new_records

    updated_records = []
    for rec in processed_records:
        seq_chars = list(str(rec.seq))
        for i in range(len(seq_chars)):
            if seq_chars[i] == '-':
                seq_chars[i] = '?'
            else:
                break
        for i in range(len(seq_chars) - 1, -1, -1):
            if seq_chars[i] == '-':
                seq_chars[i] = '?'
            else:
                break
        updated_records.append(SeqRecord(Seq(''.join(seq_chars)), id=rec.id, description=""))

    n2q_log = []
    if n2question:
        if isinstance(n2question, str) and n2question not in ('all', 'flanking'):
            target_ids = {n2question}
        elif isinstance(n2question, list):
            target_ids = set(n2question)
        elif n2question in ('all', 'flanking'):
            target_ids = {rec.id for rec in updated_records}
        else:
            target_ids = set()

        for rec in updated_records:
            if rec.id in target_ids:
                if n2question == 'flanking':
                    seq_chars = list(str(rec.seq))
                    count_n = 0
                    for i in range(len(seq_chars)):
                        if seq_chars[i] in ('N', 'n', '?', '-'):
                            if seq_chars[i] in ('N', 'n'):
                                seq_chars[i] = '?'
                                count_n += 1
                        else:
                            break
                    for i in range(len(seq_chars) - 1, -1, -1):
                        if seq_chars[i] in ('N', 'n', '?', '-'):
                            if seq_chars[i] in ('N', 'n'):
                                seq_chars[i] = '?'
                                count_n += 1
                        else:
                            break
                    if count_n > 0:
                        rec.seq = Seq(''.join(seq_chars))
                        n2q_log.append(f"{rec.id}: {count_n} flanking N/n replaced with ?")
                else:
                    seq_str = str(rec.seq)
                    count_n = seq_str.count('N') + seq_str.count('n')
                    if count_n > 0:
                        rec.seq = Seq(seq_str.replace('N', '?').replace('n', '?'))
                        n2q_log.append(f"{rec.id}: {count_n} N/n replaced with ?")

    # === Now apply gaps2question as the very last step on added sequences ===
    def replace_gap_blocks(seq, threshold, seq_id=None):
        seq_list = list(seq)
        replaced_log = []
        i = 0
        while i < len(seq_list):
            if seq_list[i] == '-':
                start = i
                while i < len(seq_list) and seq_list[i] == '-':
                    i += 1
                if (i - start) > threshold:
                    for j in range(start, i):
                        seq_list[j] = '?'
                    if seq_id:
                        replaced_log.append(f"{seq_id}: {i - start} contiguous '-' replaced with '?' at {start}–{i}")
            else:
                i += 1
        return ''.join(seq_list), replaced_log

    # Apply gaps2question only if specified
    if gaps2question:
        updated_records_dict = {rec.id: rec for rec in updated_records}
        for rec in trimmed_new_records:
            seq = str(updated_records_dict[rec.id].seq)
            new_seq, _ = replace_gap_blocks(seq, gaps2question, seq_id=rec.id)
            updated_records_dict[rec.id].seq = Seq(new_seq)
        updated_records = list(updated_records_dict.values())

    final_records = []
    for rec in updated_records:
        seq_list = list(str(rec.seq))
        for col in sorted(pound_cols):
            seq_list.insert(col, '#')
        final_records.append(SeqRecord(Seq(''.join(seq_list)), id=rec.id, description=""))

    def process_blocks(seq):
        seq_chars = list(seq)
        blocks = []
        start = 0
        for i, c in enumerate(seq_chars):
            if c == '#':
                blocks.append((start, i))
                start = i + 1
        blocks.append((start, len(seq_chars)))
        for (start, end) in blocks:
            if all(seq_chars[i] == '?' for i in range(start, end)):
                for i in range(start, end):
                    seq_chars[i] = '-'
        return ''.join(seq_chars)

    final_output = [
        SeqRecord(Seq(process_blocks(str(rec.seq))), id=rec.id, description="")
        for rec in final_records
    ]

    output_dir = os.path.dirname(output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # --- New function to find contiguous '?' blocks ---
    def find_question_blocks(seq):
        blocks = []
        seq_len = len(seq)
        i = 0
        while i < seq_len:
            if seq[i] == '?':
                start = i
                while i < seq_len and seq[i] == '?':
                    i += 1
                end = i
                blocks.append((start, end))
            else:
                i += 1
        return blocks

    # Identify gap ('?') blocks in added sequences for final logging
    gap_question_blocks_log = []
    added_ids = {rec.id for rec in trimmed_new_records}  # IDs of added sequences

    for rec in final_output:
        if rec.id in added_ids:
            q_blocks = find_question_blocks(str(rec.seq))
            for (start, end) in q_blocks:
                length = end - start
                gap_question_blocks_log.append(f"{rec.id}: ? block at positions {start}-{end} (length={length})")

    SeqIO.write(final_output, output, "fasta")
    if write_names:
        with open(output + "_terminal_names.txt", "w") as f:
            for rec in final_output:
                f.write(rec.id + "\n")

    if log:
        elapsed = time.time() - start_time
        log_lines.append(f"Final alignment: {len(final_output)} sequences")
        log_lines.append(f"Final alignment: {len(final_output[0].seq)} columns")
        log_lines.append(f"Final alignment: {sum(1 for i in range(len(final_output[0].seq)) if any(rec.seq[i] == '#' for rec in final_output))} # columns")
        log_lines.append("")
        log_lines.append("Trimmed orphan blocks from new sequences:")
        log_lines.extend(all_trim_logs or ["None"])
        if gap_question_blocks_log:
            log_lines.append("")
            log_lines.append("Gap block replacements:")
            log_lines.extend(gap_question_blocks_log)
        if n2q_log:
            log_lines.append("")
            log_lines.append("N/n to ? replacements:")
            log_lines.extend(n2q_log)
        log_lines.append("")
        log_lines.append(f"Runtime: {elapsed:.2f} seconds")
        with open(output + ".log", "w") as log_file:
            log_file.write("\n".join(log_lines))

    for file in temp_files_to_remove:
        try:
            os.remove(file)
        except Exception:
            pass

def _expand_partitioning_rounds(partitioning_round):
    if isinstance(partitioning_round, range):
        return list(partitioning_round)
    if isinstance(partitioning_round, (list, tuple)):
        return list(partitioning_round)
    if isinstance(partitioning_round, str):
        if partitioning_round == "max":
            return ["max"]
        if re.fullmatch(r"\d+-\d+", partitioning_round):
            start, end = partitioning_round.split("-", 1)
            return list(range(int(start), int(end) + 1))
    return [partitioning_round]

def _expand_partitioning_methods(partitioning_method):
    if isinstance(partitioning_method, str) and partitioning_method.lower() == "all":
        return ["conservative", "balanced", "max", "equal"]
    return [partitioning_method]

def _build_batch_output_prefix(output_file, method_name, round_value, method_batch_active, round_batch_active):
    root_output = output_file or "output"
    normalized_root = os.path.normpath(root_output)
    prefix_base = os.path.basename(normalized_root)
    if not prefix_base or prefix_base == os.path.sep:
        prefix_base = "output"

    run_directory = root_output
    if method_batch_active:
        run_directory = os.path.join(run_directory, str(method_name))
        if method_name != "max":
            run_directory = os.path.join(run_directory, f"round_{round_value}")
    elif round_batch_active:
        if method_name != "max":
            run_directory = os.path.join(run_directory, f"round_{round_value}")

    return os.path.join(run_directory, prefix_base)

def prepDyn(input_file=None,
            CSV_input=None,
            input_format="fasta",
            MSA=False,
            aligner="mafft",
            output_file=None,
            output_format="fasta",
            log=False,
            sequence_names=True,
            # Trimming parameters
            orphan_method=None,
            orphan_threshold=10,
            orphan_limit=0.05,
            orphan_action="trim",
            percentile=25,
            del_inv=False,
            multi_amplicon_min_overlap=10,
            multi_amplicon_mismatch_rate=0.05,
            multi_amplicon_action="trim",
            # Missing data parameters
            internal_method=None,
            internal_column_ranges=None,
            internal_leaves="all",
            internal_threshold=None,
            n2question=None,
            # Partitioning parameters
            partitioning_round=0,
            partitioning_method="conservative",
            partitioning_conservative="midpoint",
            partitioning_max_size=None,
            partitioning_size=None
            ):
    """
    Preprocess missing data for dynamic homology in PhyG. First, columns containing
    only gaps, orphan nucleotides, and invariant columns can be trimmed. Second,
    missing data is coded with question marks. Third, partitions are delimited in
    highly conserved regions.

    Args:
        input_file (str): Path to the input alignment file or directory. Ignored if CSV_input is provided.
        CSV_input (str): Path to a CSV/TSV file containing GenBank accession numbers. If provided,
                        sequences will be downloaded from GenBank and aligned before preprocessing.
        input_format (str): Format of the input alignment. Options: 'fasta' (default),
                            'clustal', 'phylip', or any format accepted by Biopython.
        output_file (str): Custom prefix for output files. If None, base_name from input_file is used.
        output_format (str): Output format. Default is 'fasta'.
        log (bool): Whether to write a log with wall-clock time. Default is False.
        sequence_names (bool): If True, writes a TXT file listing all sorted unique sequence names. Default is True.
        MSA (bool): Whether to perform MSA if input sequences specified in input_file are unaligned
        aligner (str): Multiple sequence aligner used when alignment is performed. Options are
                       "mafft" (default) and "clustalw".
        orphan_method (str): The trimming method. By default, trimming orphan nucleotides
                             is not performed. Options:
                            - 'percentile': trim using the 25th percentile;
                            - 'integer': trim with a manual threshold.
                            - 'adaptive_1': adaptive threshold with modification cap only.
                            - 'adaptive_2': adaptive threshold with modification cap and gap check.
                            - 'adaptive_3': adaptive threshold with modification cap, gap check, and shared homology check.
        orphan_threshold (int): Threshold used to trim orphan nucleotides if orphan_method = 'integer'.
        orphan_action (str): Action for orphan nucleotides. 'trim' (default) removes them. 
                             'push' moves them adjacent to the next block iteratively.
        percentile (float): Used with orphan_method = 'percentile' to define trimming threshold.
        del_inv (bool): Whether to trim invariant terminal columns. Default is False.
        multi_amplicon_min_overlap (int): Minimum overlap length to merge multi-amplicons when downloading from GenBank. Default is 10.
        multi_amplicon_mismatch_rate (float): Maximum mismatch rate allowed when merging overlapping multi-amplicons from GenBank. Default is 0.05.
        multi_amplicon_action (str): How to handle overlapping multi-amplicons. "trim" (default)
                                     removes one overlapping copy; "consensus" replaces the overlap
                                     with IUPAC consensus nucleotides.
        internal_method (str): Defines how to identify internal missing data. Automatic identificaton
                               of missing data is made if CSV_input is provided. Otherwise, naive
                               options to identify internal missing data are:
                               - "manual": Use column ranges;
                               - "semi": Use a threshold for gaps.
        internal_column_ranges (list): Column ranges (inclusive) if internal_method = 'manual'.
        internal_leaves (str or list): Sequences to apply internal missing data replacement
                                       if internal_method is not "None".
        internal_threshold (int): Used with internal_method = 'semi' to define gap threshold.
                                  Contiguous '-' larger than the threshold are replaced with '?'.
        n2question (str or list): If specified, replaces ambiguous nucleotide 'N' or 'n' with '?'. If None (default), n2question is not performed. If 'all', apply to all sequences. If you want to apply to only one sequence, write the name of this sequence. If you want to apply to multiple sequences but no all, wrie the list of sequences.    
        partitioning_method (str): Method of partitioning:
                                   - 'balanced': Based on initial '#' inserted with 'max', iteratively 
                                   merges adjacent blocks flanked by # if their combined length is below
                                   a threshold (the n-largest block of missing data).
                                   - 'conservative': Blocks containing only invariant columns are sorted
                                   by length and '#' column(s) inserted either at the midpoint or around
                                   the flanks of the n-largest block(s). Must define n using
                                   partitioning_round.
                                   - 'equal: Equal-length partitions are created by specifying their size
                                   or their round. If partitioning_round = 1, only 1 '#' column is
                                   inserted; if partitioning_round = 2, then 2 '#' columns are inserted.
                                   - 'max': '#' columns are inserted around blocks of missing data (every
                                   instance of '?' opening/closure but not '?' extension).
                                   - 'all': Run conservative, balanced, max, and equal in separate directories.
        partitioning_conservative (str): Conservative partition placement mode. Use 'midpoint'
                                         (default) to insert one '#' in the middle of each selected
                                         invariant block, or 'flank' to insert '#' columns on both sides
                                         of each selected invariant block.
        partitioning_max_size (int): Initial maximum partition size. If specified, the alignment is first
                                     split into equal-length partitions of this size, and the selected
                                     partitioning method is then applied independently within each partition.
        partitioning_round (int or str): Number of partitioning round. Invariant regions are sorted by length
                                         in descendant order and the n-largest block(s) partitioned using '#'.
                                         Ranges such as "0-10" generate one run per round in separate directories.
                                         If "max" is specified, pound signs are inserted arund all blocks of
                                         missing data.
        partitioning_size (int): Size of equal-length partitions if partitioning_method = 'equal'.

    Returns:
        dict: The preprocessed unaligned sequences.
    """

    # --- Start overall timer for the top-level call ---
    aligner = normalize_aligner(aligner)
    overall_start_wall_time = time.time()
    overall_start_cpu_time = time.process_time()

    # --- Initialize the shared sequence ID set for the top-level call ---
    _all_sequence_ids_shared = set()

    # --- Store original output_file for final sequence_names.txt path ---
    original_output_file_arg = output_file # Store the exact argument passed to prepDyn
    
    # --- Generate the original command line string here ---
    cmd_parts = ["prepDyn("]
    params = {
        "input_file": input_file,
        "CSV_input": CSV_input,
        "input_format": input_format,
        "MSA": MSA,
        "aligner": aligner,
        "output_file": output_file,
        "output_format": output_format,
        "log": log,
        "sequence_names": sequence_names,
        "orphan_method": orphan_method,
        "orphan_threshold": orphan_threshold,
        "orphan_action": orphan_action,
        "percentile": percentile,
        "del_inv": del_inv,
        "multi_amplicon_min_overlap": multi_amplicon_min_overlap,
        "multi_amplicon_mismatch_rate": multi_amplicon_mismatch_rate,
        "multi_amplicon_action": multi_amplicon_action,
        "internal_method": internal_method,
        "internal_column_ranges": internal_column_ranges,
        "internal_leaves": internal_leaves,
        "internal_threshold": internal_threshold,
        "n2question": n2question,
        "partitioning_method": partitioning_method,
        "partitioning_round": partitioning_round,
        "partitioning_conservative": partitioning_conservative,
        "partitioning_max_size": partitioning_max_size,
        "partitioning_size": partitioning_size,
    }

    param_strs = []
    for k, v in params.items():
        # Represent strings with quotes, lists/tuples as is, others repr()
        if isinstance(v, str):
            param_strs.append(f"{k}='{v}'")
        elif isinstance(v, (list, tuple)):
            param_strs.append(f"{k}={v}")
        elif v is not None: 
            param_strs.append(f"{k}={repr(v)}")
    
    cmd_parts.append(", ".join(param_strs))
    cmd_parts.append(")")
    original_cmd_line = "".join(cmd_parts)

    expanded_methods = _expand_partitioning_methods(partitioning_method)
    expanded_rounds = _expand_partitioning_rounds(partitioning_round)
    method_batch_active = len(expanded_methods) > 1
    round_batch_active = len(expanded_rounds) > 1

    if method_batch_active or round_batch_active:
        batch_start_wall_time = time.time()
        batch_start_cpu_time = time.process_time()
        batch_root_output = output_file or "output"
        batch_root_dir = os.path.normpath(batch_root_output)
        batch_root_name = os.path.basename(batch_root_dir) or "output"
        batch_runtime_log_path = os.path.join(batch_root_output, f"{batch_root_name}_overall_log.txt")
        batch_input_file = input_file
        batch_csv_input = CSV_input
        batch_input_format = input_format
        batch_msa = MSA
        cache_root = None

        try:
            if CSV_input is not None:
                cache_root = os.path.join(batch_root_output, "_gb2msa_cache")
                cache_output_prefix = os.path.join(cache_root, "cache")
                cache_alignment_dir = os.path.join(cache_root, "alignments")
                os.makedirs(cache_alignment_dir, exist_ok=True)

                print("Running GB2MSA on GenBank input once for batch partitioning...")
                cleaned_files, _, gene_names = GB2MSA(
                    CSV_input,
                    output_prefix=cache_output_prefix,
                    write_names=False,
                    log=False,
                    multi_amplicon_min_overlap=multi_amplicon_min_overlap,
                    multi_amplicon_mismatch_rate=multi_amplicon_mismatch_rate,
                    multi_amplicon_action=multi_amplicon_action,
                    aligner=aligner,
                    return_gene_logs=True,
                    write_run_log=False
                )

                for gene_name, cleaned_file in zip(gene_names, cleaned_files):
                    cached_alignment_path = os.path.join(cache_alignment_dir, f"{gene_name}.fasta")
                    shutil.copy2(cleaned_file, cached_alignment_path)

                batch_input_file = cache_alignment_dir
                batch_csv_input = None
                batch_input_format = "fasta"
                batch_msa = False

            for method_name in expanded_methods:
                current_rounds = [expanded_rounds[0]] if method_name == "max" else expanded_rounds
                for round_value in current_rounds:
                    batch_output_prefix = _build_batch_output_prefix(
                        output_file,
                        method_name,
                        round_value,
                        method_batch_active=method_batch_active,
                        round_batch_active=round_batch_active
                    )
                    prepDyn(
                        input_file=batch_input_file,
                        CSV_input=batch_csv_input,
                        input_format=batch_input_format,
                        MSA=batch_msa,
                        aligner=aligner,
                        output_file=batch_output_prefix,
                        output_format=output_format,
                        log=log,
                        sequence_names=sequence_names,
                        orphan_method=orphan_method,
                        orphan_threshold=orphan_threshold,
                        orphan_action=orphan_action,
                        percentile=percentile,
                        del_inv=del_inv,
                        multi_amplicon_min_overlap=multi_amplicon_min_overlap,
                        multi_amplicon_mismatch_rate=multi_amplicon_mismatch_rate,
                        multi_amplicon_action=multi_amplicon_action,
                        internal_method=internal_method,
                        internal_column_ranges=internal_column_ranges,
                        internal_leaves=internal_leaves,
                        internal_threshold=internal_threshold,
                        n2question=n2question,
                        partitioning_round=round_value,
                        partitioning_method=method_name,
                        partitioning_conservative=partitioning_conservative,
                        partitioning_max_size=partitioning_max_size,
                        partitioning_size=partitioning_size
                    )
        finally:
            if cache_root and os.path.isdir(cache_root):
                shutil.rmtree(cache_root, ignore_errors=True)

        if log:
            batch_end_wall_time = time.time()
            batch_end_cpu_time = time.process_time()
            total_batch_wall_time = batch_end_wall_time - batch_start_wall_time
            total_batch_cpu_time = batch_end_cpu_time - batch_start_cpu_time
            os.makedirs(os.path.dirname(batch_runtime_log_path) or ".", exist_ok=True)
            with open(batch_runtime_log_path, "w") as batch_log_file:
                batch_log_file.write("--- Overall prepDyn Batch Execution Summary ---\n")
                batch_log_file.write(f"Command used:\n{original_cmd_line}\n\n")
                batch_log_file.write(f"Partitioning methods: {expanded_methods}\n")
                batch_log_file.write(f"Partitioning rounds: {expanded_rounds}\n\n")
                batch_log_file.write(f"Total Wall-clock time: {total_batch_wall_time:.8f} seconds\n")
                batch_log_file.write(f"Total CPU time: {total_batch_cpu_time:.8f} seconds\n")
                batch_log_file.write("\nNote: Individual run logs provide detailed information.\n")
        return None


    def _prepDyn_recursive(input_val, # Renamed to input_val to be more generic for file path or dict
                           CSV_input,
                           input_format,
                           MSA,
                           aligner,
                           output_file, # This is crucial: the output prefix for THIS specific call
                           output_format,
                           log,
                           sequence_names, # This will be False for recursive calls
                           _all_sequence_ids, # This is the shared set
                           _is_top_level_call, # Flag for controlling final write
                           _original_cmd_line, # Parameter to pass the original command
                           _gb_gene_log_metadata,
                           orphan_method,
                           orphan_threshold,
                           orphan_action,
                           percentile,
                           del_inv,
                           multi_amplicon_min_overlap,
                           multi_amplicon_mismatch_rate,
                           multi_amplicon_action,
                           internal_method,
                           internal_column_ranges,
                           internal_leaves,
                           internal_threshold,
                           n2question,
                           partitioning_round,
                           partitioning_method,
                           partitioning_conservative,
                           partitioning_max_size,
                           partitioning_size
                           ):

        # Start timers if logging is enabled (only for current processing, not recursive overall)
        if log:
            start_wall_time = time.time()
            start_cpu_time = time.process_time()

        # --- Pre-processing for output_file path handling for *this* recursive call ---
        current_output_dir = os.path.dirname(output_file)
        if not current_output_dir:
            current_output_dir = "." # Default to current directory if no path in output_file
        os.makedirs(current_output_dir, exist_ok=True)


        # Step 1: Run GB2MSA if GenBank input is provided
        if CSV_input is not None:
            # (This part of the logic remains unchanged)
            print("Running GB2MSA on GenBank input...")
            
            gb_output_prefix_for_gb2msa = output_file 
            cleaned_files, gb1_logs, gene_names = GB2MSA(
                CSV_input,
                output_prefix=gb_output_prefix_for_gb2msa,
                write_names=False,
                log=log,
                multi_amplicon_min_overlap=multi_amplicon_min_overlap,
                multi_amplicon_mismatch_rate=multi_amplicon_mismatch_rate,
                multi_amplicon_action=multi_amplicon_action,
                aligner=aligner,
                return_gene_logs=True,
                write_run_log=False
            )
            
            for gene_name, file_path_from_gb2msa in zip(gene_names, cleaned_files):
                file_basename_no_ext = os.path.splitext(os.path.basename(file_path_from_gb2msa))[0]
                gene_name_part = file_basename_no_ext.replace("_aligned", "").replace("_GB2MSA", "")
                base_from_original_output = os.path.basename(original_output_file_arg) if original_output_file_arg else ""
                if base_from_original_output and base_from_original_output != os.path.basename(os.path.normpath(original_output_file_arg)):
                    gene_specific_prefix_base = f"{base_from_original_output}_{gene_name_part}"
                elif os.path.isdir(original_output_file_arg):
                    gene_specific_prefix_base = gene_name_part
                else:
                    gene_specific_prefix_base = gene_name_part

                specific_output_prefix_for_recursion = os.path.join(current_output_dir, gene_specific_prefix_base)
                try:
                    alignment = AlignIO.read(file_path_from_gb2msa, "fasta")
                    alignment_dict = {record.id: str(record.seq) for record in alignment}
                except ValueError as e:
                    if "Sequences must all be the same length" in str(e):
                        # Read sequences and pad shorter ones with gaps
                        sequences = list(SeqIO.parse(file_path_from_gb2msa, "fasta"))
                        max_length = max(len(record.seq) for record in sequences)
                        alignment_dict = {}
                        for record in sequences:
                            seq_str = str(record.seq)
                            padded_seq = seq_str + "-" * (max_length - len(seq_str))
                            alignment_dict[record.id] = padded_seq
                    else:
                        raise
                _all_sequence_ids.update(alignment_dict.keys())
                _prepDyn_recursive(input_val=alignment_dict,
                                CSV_input=None, input_format="dict", MSA=MSA, aligner=aligner,
                                orphan_method=orphan_method, orphan_threshold=orphan_threshold, orphan_action=orphan_action, percentile=percentile, del_inv=del_inv, multi_amplicon_min_overlap=multi_amplicon_min_overlap, multi_amplicon_mismatch_rate=multi_amplicon_mismatch_rate, multi_amplicon_action=multi_amplicon_action,
                                internal_method=internal_method, internal_column_ranges=internal_column_ranges, internal_leaves=internal_leaves, internal_threshold=internal_threshold,
                                n2question=n2question, partitioning_method=partitioning_method, partitioning_round=partitioning_round, partitioning_conservative=partitioning_conservative, partitioning_max_size=partitioning_max_size, partitioning_size=partitioning_size,
                                output_format=output_format, log=log, sequence_names=False,
                                _all_sequence_ids=_all_sequence_ids, _is_top_level_call=False, _original_cmd_line=_original_cmd_line,
                                _gb_gene_log_metadata=gb1_logs.get(gene_name, {}),
                                output_file=specific_output_prefix_for_recursion)
            return

        # Step 2: If a folder is provided, process each alignment inside
        if isinstance(input_val, str) and os.path.isdir(input_val):
            # (This part of the logic remains unchanged)
            processed_any_file = False
            base_name_for_output_prefix = ""
            if original_output_file_arg and os.path.isdir(original_output_file_arg):
                base_name_for_output_prefix = os.path.basename(os.path.normpath(original_output_file_arg))
            elif original_output_file_arg:
                 base_name_for_output_prefix = os.path.basename(original_output_file_arg)

            for file_name in os.listdir(input_val):
                file_extension = os.path.splitext(file_name)[1].lstrip('.')
                if file_extension == input_format:
                    processed_any_file = True
                    file_path = os.path.join(input_val, file_name)
                    base_name_of_file = os.path.splitext(file_name)[0]
                    
                    if base_name_for_output_prefix:
                        specific_output_prefix = os.path.join(current_output_dir, f"{base_name_for_output_prefix}_{base_name_of_file}")
                    else: 
                        specific_output_prefix = os.path.join(current_output_dir, base_name_of_file)

                    current_file_alignment = None
                    if MSA:
                        # (MSA logic for folder input unchanged)
                        print(f"Processing unaligned file with {aligner}: {file_path}")
                        temp_dir = current_output_dir if current_output_dir != "." else tempfile.gettempdir()
                        os.makedirs(temp_dir, exist_ok=True)
                        with tempfile.NamedTemporaryFile(mode="w", delete=False, dir=temp_dir, suffix=f".{input_format}") as tmp_in:
                            sequences = list(SeqIO.parse(file_path, input_format))
                            SeqIO.write(sequences, tmp_in, "fasta")
                            tmp_in_path = tmp_in.name
                        tmp_out_path = os.path.join(temp_dir, f"{os.path.basename(tmp_in_path)}_aligned.fasta")
                        try:
                            run_multiple_sequence_alignment(tmp_in_path, tmp_out_path, aligner=aligner)
                            try:
                                current_file_alignment = {record.id: str(record.seq) for record in AlignIO.read(tmp_out_path, "fasta")}
                            except ValueError as e:
                                if "Sequences must all be the same length" in str(e):
                                    # Read sequences and pad shorter ones with gaps
                                    sequences = list(SeqIO.parse(tmp_out_path, "fasta"))
                                    max_length = max(len(record.seq) for record in sequences)
                                    current_file_alignment = {}
                                    for record in sequences:
                                        seq_str = str(record.seq)
                                        padded_seq = seq_str + "-" * (max_length - len(seq_str))
                                        current_file_alignment[record.id] = padded_seq
                            except Exception as e:
                                if os.path.exists(tmp_in_path): os.remove(tmp_in_path)
                                if os.path.exists(tmp_out_path): os.remove(tmp_out_path)
                        finally:
                            if os.path.exists(tmp_in_path): os.remove(tmp_in_path)
                            if os.path.exists(tmp_out_path): os.remove(tmp_out_path)
                    else:
                        try:
                            alignment_temp = AlignIO.read(file_path, input_format)
                            current_file_alignment = {record.id: str(record.seq) for record in alignment_temp}
                        except ValueError as e:
                            if "Sequences must all be the same length" in str(e):
                                # Read sequences and pad shorter ones with gaps
                                sequences = list(SeqIO.parse(file_path, input_format))
                                max_length = max(len(record.seq) for record in sequences)
                                current_file_alignment = {}
                                for record in sequences:
                                    seq_str = str(record.seq)
                                    padded_seq = seq_str + "-" * (max_length - len(seq_str))
                                    current_file_alignment[record.id] = padded_seq
                            else:
                                raise

                    if current_file_alignment:
                        _all_sequence_ids.update(current_file_alignment.keys())
                        _prepDyn_recursive(input_val=current_file_alignment,
                                CSV_input=None, input_format="dict", MSA=False, aligner=aligner,
                                orphan_method=orphan_method, orphan_threshold=orphan_threshold, orphan_action=orphan_action, percentile=percentile, del_inv=del_inv, multi_amplicon_min_overlap=multi_amplicon_min_overlap, multi_amplicon_mismatch_rate=multi_amplicon_mismatch_rate, multi_amplicon_action=multi_amplicon_action,
                                internal_method=internal_method, internal_column_ranges=internal_column_ranges, internal_leaves=internal_leaves, internal_threshold=internal_threshold,
                                n2question=n2question, partitioning_method=partitioning_method, partitioning_round=partitioning_round, partitioning_conservative=partitioning_conservative, partitioning_max_size=partitioning_max_size, partitioning_size=partitioning_size,
                                output_format=output_format, log=log, sequence_names=False,
                                _all_sequence_ids=_all_sequence_ids, _is_top_level_call=False, _original_cmd_line=_original_cmd_line,
                                _gb_gene_log_metadata=None,
                                output_file=specific_output_prefix)
            if not processed_any_file:
                print(f"WARNING: No files with extension '.{input_format}' found in directory: {input_val}")
            return


        # Step 3: Read and process alignment
        alignment = None
        if isinstance(input_val, dict):
            alignment = input_val
        elif isinstance(input_val, str) and os.path.isfile(input_val):
            # (MSA logic for single file input unchanged)
            if MSA:
                temp_dir = current_output_dir if current_output_dir != "." else tempfile.gettempdir()
                os.makedirs(temp_dir, exist_ok=True)
                with tempfile.NamedTemporaryFile(mode="w", delete=False, dir=temp_dir, suffix=f".{input_format}") as tmp_in:
                    sequences = list(SeqIO.parse(input_val, input_format))
                    SeqIO.write(sequences, tmp_in, "fasta")
                    tmp_in_path = tmp_in.name
                tmp_out_path = os.path.join(temp_dir, f"{os.path.basename(tmp_in_path)}_aligned.fasta")
                try:
                    run_multiple_sequence_alignment(tmp_in_path, tmp_out_path, aligner=aligner)
                    try:
                        alignment_obj = AlignIO.read(tmp_out_path, "fasta")
                    except ValueError as e:
                        if "Sequences must all be the same length" in str(e):
                            # Read sequences and pad shorter ones with gaps
                            sequences = list(SeqIO.parse(tmp_out_path, "fasta"))
                            max_length = max(len(record.seq) for record in sequences)
                            alignment_dict_padded = {}
                            for record in sequences:
                                seq_str = str(record.seq)
                                padded_seq = seq_str + "-" * (max_length - len(seq_str))
                                alignment_dict_padded[record.id] = padded_seq
                            alignment = alignment_dict_padded
                            _all_sequence_ids.update(alignment.keys())
                            alignment_obj = None
                        else:
                            raise
                finally:
                    if os.path.exists(tmp_in_path): os.remove(tmp_in_path)
                    if os.path.exists(tmp_out_path): os.remove(tmp_out_path)
                
                if alignment_obj is not None:
                    alignment = {record.id: str(record.seq) for record in alignment_obj}
                    _all_sequence_ids.update(alignment.keys())
            else:
                try:
                    alignment_obj = AlignIO.read(input_val, input_format)
                except ValueError as e:
                    if "Sequences must all be the same length" in str(e):
                        # Read sequences and pad shorter ones with gaps
                        sequences = list(SeqIO.parse(input_val, input_format))
                        max_length = max(len(record.seq) for record in sequences)
                        alignment_dict_padded = {}
                        for record in sequences:
                            seq_str = str(record.seq)
                            padded_seq = seq_str + "-" * (max_length - len(seq_str))
                            alignment_dict_padded[record.id] = padded_seq
                        alignment = alignment_dict_padded
                        _all_sequence_ids.update(alignment.keys())
                    else:
                        raise
                else:
                    alignment = {record.id: str(record.seq) for record in alignment_obj}
                    _all_sequence_ids.update(alignment.keys())
        else:
            if isinstance(input_val, str):
                raise FileNotFoundError(
                    f"Input path not found: {input_val} "
                    f"(current working directory: {os.getcwd()})"
                )
            raise ValueError(f"Invalid input_val type: {type(input_val).__name__}.")

        alignment = normalize_pound_columns(alignment)


        ### FIX: PART 1 - Capture "before" statistics right after loading ###
        # Initialize variables to store the "before" summary
        num_seqs_before, aln_length_before, total_nt_before, total_gaps_before, total_ns_before = 0, 0, 0, 0, 0
        
        # If logging is enabled, calculate and store the initial state of the alignment
        if log:
            num_seqs_before = len(alignment)
            aln_length_before = len(next(iter(alignment.values())) if alignment else 0)
            total_nt_before = sum(c.upper() in "ACGT" for seq in alignment.values() for c in seq)
            total_gaps_before = sum(seq.count("-") for seq in alignment.values())
            total_ns_before = sum(c in "Nn" for seq in alignment.values() for c in seq)
        ### END FIX PART 1 ###


        # 3.1 Remove columns with gaps in all leaves
        alignment = remove_all_gap_columns(alignment)

        # 3.2 Trim orphan nucleotides
        orphan_log = None
        if orphan_method == "percentile":
            orphan_threshold_val = calculate_orphan_threshold_from_percentile(alignment, percentile, terminal_only=True)
            if log:
                alignment, orphan_log = delete_orphan_nucleotides2(alignment, orphan_threshold_val, log_changes=True, orphan_action=orphan_action, orphan_method="integer", orphan_limit=orphan_limit)
            else:
                alignment = delete_orphan_nucleotides2(alignment, orphan_threshold_val, orphan_action=orphan_action, orphan_method="integer", orphan_limit=orphan_limit)
        elif orphan_method in ["integer", "adaptive_1", "adaptive_2", "adaptive_3"]:
            if log:
                alignment, orphan_log = delete_orphan_nucleotides2(alignment, orphan_threshold, log_changes=True, orphan_action=orphan_action, orphan_method=orphan_method, orphan_limit=orphan_limit)
            else:
                alignment = delete_orphan_nucleotides2(alignment, orphan_threshold, orphan_action=orphan_action, orphan_method=orphan_method, orphan_limit=orphan_limit)


        # 3.3 Replace terminal gaps with ?
        alignment = replace_terminal_gaps_dict(alignment)

        # 3.4 Trim invariant columns
        removed_cols = []
        if del_inv:
            alignment = remove_non_informative_positions(alignment, removed_indices=removed_cols)

        alignment = replace_terminal_gaps_dict(alignment)

        # 3.5 Replace internal gaps with ?
        if internal_method == "manual":
            alignment = replace_dashes_with_question_marks(alignment=alignment,
                                                           internal_column_ranges=internal_column_ranges,
                                                           internal_leaves=internal_leaves,
                                                           internal_method="manual")
        elif internal_method == "semi":
            alignment = replace_dashes_with_question_marks(alignment=alignment,
                                                           internal_leaves=internal_leaves,
                                                           internal_method="semi",
                                                           internal_threshold=internal_threshold)

        # 3.6 Replace ambiguous nucleotides N/n with ?
        n_blocks = []
        if n2question is not None:
            alignment, n_blocks = n2question_func(alignment, leaves=n2question, log=True)


        # 3.7 Partitioning
        partitioning_log_entry = None
        conservative_blocks = []
        conservative_blocks_by_chunk = []
        balanced_metadata = None
        balanced_metadata_by_chunk = []
        partitioning_max_size_positions = []
        partitioning_method_positions = []

        if partitioning_max_size:
            chunked_alignments = split_alignment_into_equal_chunks(alignment, partitioning_max_size)
            processed_chunks = []
            chunk_warnings = []

            for chunk_alignment in chunked_alignments:
                processed_chunk, chunk_conservative_blocks, chunk_log_entry, chunk_balanced_metadata = apply_partitioning_strategy(
                    chunk_alignment,
                    partitioning_method=partitioning_method,
                    partitioning_round=partitioning_round,
                    partitioning_conservative=partitioning_conservative,
                    partitioning_size=partitioning_size
                )
                processed_chunks.append(processed_chunk)
                conservative_blocks_by_chunk.append(chunk_conservative_blocks)
                balanced_metadata_by_chunk.append(chunk_balanced_metadata)
                if chunk_log_entry:
                    chunk_warnings.append(chunk_log_entry)

            alignment, partitioning_chunk_metadata = join_alignment_chunks_with_pounds(
                processed_chunks,
                return_metadata=True
            )
            partitioning_max_size_positions = partitioning_chunk_metadata["max_size_positions"]
            partitioning_method_positions = partitioning_chunk_metadata["method_positions"]
            if chunk_warnings:
                partitioning_log_entry = "\n".join(chunk_warnings)
                print(f"\n{partitioning_log_entry}\n")
        else:
            alignment, conservative_blocks, partitioning_log_entry, balanced_metadata = apply_partitioning_strategy(
                alignment,
                partitioning_method=partitioning_method,
                partitioning_round=partitioning_round,
                partitioning_conservative=partitioning_conservative,
                partitioning_size=partitioning_size
            )
            if partitioning_log_entry:
                print(f"\n{partitioning_log_entry}\n")

        refinement_question2hyphen(alignment)
        alignment = remove_columns_with_W(alignment)
        alignment = normalize_pound_columns(alignment)
        alignment = remove_adjacent_pound_columns(alignment)


        # Step 4: Write output file
        records = [SeqRecord(Seq(seq), id=key, description="") for key, seq in alignment.items()]
        final_output_path_prefix = output_file
        output_directory_for_final_write = os.path.dirname(final_output_path_prefix)
        if output_directory_for_final_write and not os.path.exists(output_directory_for_final_write):
            os.makedirs(output_directory_for_final_write, exist_ok=True)

        output_path = f"{final_output_path_prefix}_preprocessed.{output_format}"
        with open(output_path, "w") as output_handle:
            SeqIO.write(records, output_handle, output_format)

        # Step 5: Write log (local to this gene/file)
        if log:
            end_wall_time = time.time()
            end_cpu_time = time.process_time()
            wall_time = end_wall_time - start_wall_time
            cpu_time = end_cpu_time - start_cpu_time

            log_path = f"{final_output_path_prefix}_log.txt"
            with open(log_path, "w") as log_file:
                log_file.write("--- Command used ---\n")
                log_file.write(f"{_original_cmd_line}\n\n")

                ### FIX: PART 2 - Use the stored "before" values for logging ###
                log_file.write("--- Step 1: Summary before preprocessing ---\n")
                log_file.write(f"No. sequences: {num_seqs_before}\n")
                log_file.write(f"No. columns: {aln_length_before}\n")
                log_file.write(f"Total no. nucleotides (A/C/G/T only): {total_nt_before} bp\n")
                log_file.write(f"Total no. gaps (-): {total_gaps_before}\n")
                log_file.write(f"Total no. IUPAC N: {total_ns_before}\n\n")
                ### END FIX PART 2 ###

                if _gb_gene_log_metadata:
                    source_summary_lines = []
                    multi_amplicon_lines = []
                    for seq_name, metadata in _gb_gene_log_metadata.items():
                        source_lines = _build_source_origin_log(
                            metadata["local_sources"],
                            metadata["genbank_sources"]
                        )
                        if source_lines:
                            source_summary_lines.append(f"{seq_name}:")
                            source_summary_lines.extend([f"  - {line}" for line in source_lines])

                        if len(metadata["requested_sources"]) > 1:
                            events = _build_multi_amplicon_log(
                                metadata["requested_sources"],
                                metadata["loaded_sources"],
                                metadata["failed_sources"],
                                metadata["merge_events"],
                                metadata["final_groups"]
                            )
                            multi_amplicon_lines.append(f"{seq_name}:")
                            multi_amplicon_lines.extend([f"  - {event}" for event in events])
                    if source_summary_lines:
                        log_file.write("--- Step 1: Input sequence sources ---\n")
                        log_file.write("\n".join(source_summary_lines))
                        log_file.write("\n\n")
                    if multi_amplicon_lines:
                        log_file.write("--- Step 1: Multi-amplicon processing ---\n")
                        log_file.write("\n".join(multi_amplicon_lines))
                        log_file.write("\n\n")

                if del_inv:
                    log_file.write("--- Step 2: Trimming (invariant columns) ---\n")
                    log_file.write(f"{removed_cols}\n\n")
                if orphan_log:
                    log_file.write("--- Step 2: Trimming (orphan nucleotides) ---\n")
                    log_file.write(f"{orphan_log}\n\n")

                if n_blocks:
                    log_file.write("--- Step 3: Missing data identification (Ns replaced with '?') ---\n")
                    for seq_name, start, end in n_blocks:
                        log_file.write(f"{seq_name}: {start}–{end}\n")
                    log_file.write("\n")

                missing_partition_log = detect_fully_missing_partitions(alignment)
                if missing_partition_log:
                    log_file.write("--- Step 3: Missing data identification (gaps replaced with '?') ---\n")
                    log_file.write(f"{missing_partition_log}\n\n")

                log_file.write("--- Step 4: Partitioning ---\n")
                if partitioning_log_entry:
                    # If partitioning failed, write the stored warning message.
                    log_file.write(f"{partitioning_log_entry}\n\n")
                else:
                    # Otherwise, report success as usual.
                    pound_indices = get_pound_column_indices_from_output_file(output_path, output_format)
                    if pound_indices:
                        if partitioning_max_size:
                            log_file.write(f"Method used: initial partitioning (partitioning_max_size={partitioning_max_size})")
                            secondary_applied = (
                                partitioning_method == "max" or
                                (partitioning_method in ["conservative", "balanced"] and partitioning_round > 0) or
                                (partitioning_method == "equal" and (partitioning_round > 0 or partitioning_size))
                            )
                            if secondary_applied:
                                log_file.write(f" + {partitioning_method}")
                                if partitioning_method == "equal" and partitioning_size:
                                    log_file.write(f" (partitioning_size={partitioning_size})")
                                elif partitioning_method in ["conservative", "balanced"]:
                                    log_file.write(f" (partitioning_round={partitioning_round})")
                                    if partitioning_method == "conservative":
                                        log_file.write(f" (partitioning_conservative={partitioning_conservative})")
                                elif partitioning_method == "max":
                                    log_file.write(" (inserted at '?' block boundaries)")
                        else:
                            log_file.write(f"Method used: {partitioning_method}")
                            if partitioning_method == "equal" and partitioning_size:
                                log_file.write(f" (partitioning_size={partitioning_size})")
                            elif partitioning_method in ["conservative", "balanced"]:
                                log_file.write(f" (partitioning_round={partitioning_round})")
                                if partitioning_method == "conservative":
                                    log_file.write(f" (partitioning_conservative={partitioning_conservative})")
                            elif partitioning_method == "max":
                                log_file.write(" (inserted at '?' block boundaries)")
                        log_file.write("\n")
                        if partitioning_max_size:
                            log_file.write(
                                "Position of pound sign columns '#' from initial partitioning "
                                f"(partitioning_max_size): {partitioning_max_size_positions}\n"
                            )
                            if partitioning_method_positions:
                                log_file.write(
                                    "Position of pound sign columns '#' from subsequent "
                                    f"partitioning_method='{partitioning_method}': {partitioning_method_positions}\n"
                                )
                            else:
                                log_file.write(
                                    "Position of pound sign columns '#' from subsequent "
                                    f"partitioning_method='{partitioning_method}': []\n"
                                )
                        if partitioning_method == "conservative":
                            if partitioning_max_size and conservative_blocks_by_chunk:
                                for chunk_idx, chunk_blocks in enumerate(conservative_blocks_by_chunk, start=1):
                                    top_blocks = chunk_blocks[:partitioning_round]
                                    block_sizes = [block['length'] for block in top_blocks]
                                    log_file.write(
                                        f"Size of the {partitioning_round} largest invariant block(s) "
                                        f"in partition {chunk_idx}: {block_sizes}\n"
                                    )
                            elif conservative_blocks:
                                top_blocks = conservative_blocks[:partitioning_round]
                                block_sizes = [block['length'] for block in top_blocks]
                                log_file.write(
                                    f"Size of the {partitioning_round} largest invariant block(s): "
                                    f"{block_sizes}\n"
                                )
                        if partitioning_method == "balanced":
                            if partitioning_max_size and balanced_metadata_by_chunk:
                                for chunk_idx, chunk_metadata in enumerate(balanced_metadata_by_chunk, start=1):
                                    if chunk_metadata is not None:
                                        log_file.write(
                                            f"Size of the {chunk_metadata['reference_block_rank']} largest block of missing data "
                                            f"used as reference in partition {chunk_idx}: {chunk_metadata['reference_block_size']}\n"
                                        )
                            elif balanced_metadata is not None:
                                log_file.write(
                                    f"Size of the {balanced_metadata['reference_block_rank']} largest block of missing data "
                                    f"used as reference: {balanced_metadata['reference_block_size']}\n"
                                )
                        if partitioning_method == "equal" and partitioning_round > 0 and not partitioning_size:
                            partition_sizes = get_partition_sizes_from_output_file(output_path, output_format)
                            log_file.write(
                                f"Size of each equal-length partition: {partition_sizes}\n"
                            )
                        log_file.write(
                            "Position of pound sign columns '#' in final output alignment: "
                            f"{pound_indices}\n\n"
                        )
                    elif partitioning_method and partitioning_method != 'none':
                        log_file.write(f"Method used: {partitioning_method}\nNo partitions were inserted based on the criteria.\n\n")
                    else:
                        log_file.write("No partitioning method was specified.\n\n")

                summary_post = compute_summary_after(alignment)
                log_file.write("--- Summary after preprocessing ---\n")
                log_file.write(f"No. sequences: {summary_post['num_seqs']}\n")
                log_file.write(f"No. columns: {summary_post['aln_length']}\n")
                log_file.write(f"No. pound sign columns (#): {summary_post['total_pound']}\n")
                log_file.write(f"Total no. nucleotides (A/C/G/T): {summary_post['total_nt']} bp\n")
                log_file.write(f"Total no. gaps (-): {summary_post['total_gaps']}\n")
                log_file.write(f"Total no. IUPAC N: {summary_post['total_ns']}\n")
                log_file.write(f"Total no. missing values (?): {summary_post['total_missing']}\n\n")

                log_file.write("--- Run time ---\n")
                log_file.write(f"Wall-clock time: {wall_time:.8f} seconds\n")
                log_file.write(f"CPU time: {cpu_time:.8f} seconds\n")
        
        return alignment

    # --- Initial call to the recursive helper function (Unchanged) ---
    final_processed_alignment = _prepDyn_recursive(input_val=input_file,
                                                   CSV_input=CSV_input,
                                                   input_format=input_format,
                                                   MSA=MSA,
                                                   aligner=aligner,
                                                   output_file=output_file,
                                                   output_format=output_format,
                                                   log=log,
                                                   sequence_names=sequence_names,
                                                   _all_sequence_ids=_all_sequence_ids_shared,
                                                   _is_top_level_call=True,
                                                   _original_cmd_line=original_cmd_line,
                                                   _gb_gene_log_metadata=None,
                                                   orphan_method=orphan_method,
                                                   orphan_threshold=orphan_threshold,
                                                   orphan_action=orphan_action,
                                                   percentile=percentile,
                                                   del_inv=del_inv,
                                                   multi_amplicon_min_overlap=multi_amplicon_min_overlap,
                                                   multi_amplicon_mismatch_rate=multi_amplicon_mismatch_rate,
                                                   multi_amplicon_action=multi_amplicon_action,
                                                   internal_method=internal_method,
                                                   internal_column_ranges=internal_column_ranges,
                                                   internal_leaves=internal_leaves,
                                                   internal_threshold=internal_threshold,
                                                   n2question=n2question,
                                                   partitioning_round=partitioning_round,
                                                   partitioning_method=partitioning_method,
                                                   partitioning_conservative=partitioning_conservative,
                                                   partitioning_max_size=partitioning_max_size,
                                                   partitioning_size=partitioning_size)

    # --- Final sequence_names.txt and overall log writing (Unchanged) ---
    if sequence_names:
        sorted_unique_names = sorted(list(_all_sequence_ids_shared))
        names_file_path = None
        if original_output_file_arg:
            output_dir_for_final_names = os.path.dirname(original_output_file_arg)
            if not output_dir_for_final_names: output_dir_for_final_names = "."
            output_base_for_final_names = os.path.basename(original_output_file_arg)
            if os.path.isdir(original_output_file_arg) and not output_base_for_final_names:
                output_base_for_final_names = os.path.basename(os.path.normpath(original_output_file_arg))
            names_file_path = os.path.join(output_dir_for_final_names, f"{output_base_for_final_names}_sequence_names.txt")
        elif CSV_input: names_file_path = "output_sequence_names.txt"
        elif isinstance(input_file, str) and os.path.isdir(input_file):
            base_for_names = os.path.basename(os.path.normpath(input_file))
            names_file_path = f"{base_for_names}_sequence_names.txt"
        elif isinstance(input_file, str) and not os.path.isdir(input_file):
            base_for_names = os.path.splitext(os.path.basename(input_file))[0]
            names_file_path = f"{base_for_names}_sequence_names.txt"
        else: names_file_path = "alignment_sequence_names.txt"
        
        if names_file_path:
            os.makedirs(os.path.dirname(names_file_path), exist_ok=True)
            with open(names_file_path, 'w') as nf:
                for name in sorted_unique_names:
                    nf.write(f"{name}\n")

    if log and (CSV_input is not None or (isinstance(input_file, str) and os.path.isdir(input_file))):
        overall_end_wall_time = time.time()
        overall_end_cpu_time = time.process_time()
        total_wall_time = overall_end_wall_time - overall_start_wall_time
        total_cpu_time = overall_end_cpu_time - overall_start_cpu_time

        overall_log_path = None
        if original_output_file_arg:
            output_dir_for_overall_log = os.path.dirname(original_output_file_arg)
            if not output_dir_for_overall_log: output_dir_for_overall_log = "."
            output_base_for_overall_log = os.path.basename(original_output_file_arg)
            if os.path.isdir(original_output_file_arg) and not output_base_for_overall_log:
                output_base_for_overall_log = os.path.basename(os.path.normpath(original_output_file_arg))
            overall_log_path = os.path.join(output_dir_for_overall_log, f"{output_base_for_overall_log}_overall_log.txt")
        elif CSV_input: overall_log_path = "overall_prepDyn_log.txt"
        elif isinstance(input_file, str) and os.path.isdir(input_file):
            base_for_overall_log = os.path.basename(os.path.normpath(input_file))
            overall_log_path = f"{base_for_overall_log}_overall_log.txt"
        
        if overall_log_path:
            os.makedirs(os.path.dirname(overall_log_path), exist_ok=True)
            with open(overall_log_path, 'w') as of:
                of.write("--- Overall prepDyn Execution Summary ---\n")
                of.write(f"Command used:\n{original_cmd_line}\n\n")
                of.write(f"Total Wall-clock time: {total_wall_time:.8f} seconds\n")
                of.write(f"Total CPU time: {total_cpu_time:.8f} seconds\n")
                of.write("\nNote: Individual gene/alignment logs provide detailed information.\n")

    if CSV_input or (isinstance(input_file, str) and os.path.isdir(input_file)):
        return None 
    else:
        return final_processed_alignment
