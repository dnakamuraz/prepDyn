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

#######################
# AUXILIARY FUNCTIONS #
#######################

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
    
############################################   
# MAIN FUNCTIONS: STEP 1. DATA COLLECTION  #
############################################

def GB2MSA_1(input_file, output_prefix, delimiter=',', write_names=True):
    """
    Downloads GenBank sequences based on accession numbers in a CSV/TSV file and aligns them by gene using 
    MAFFT. If two fragments of the same locus are concatenated with no overlap between them, the space
    between them will be treaed as missing data (15 Ws will flank these blocks of missing data, which will
    be used to track these regions and be replaced with question marks in GB2MSA_2).

    Parameters:
    -----------
    input_file : str
        Path to the CSV or TSV input file. The first column should contain sequence names (sample identifiers).
        The first row should contain gene names starting from the second column. Cells contain GenBank 
        accession numbers (one or more separated by slashes). "NA", empty cells, or dashes are ignored.

    output_prefix : str
        Prefix used for naming intermediate FASTA files and final aligned output files.

    delimiter : str, optional (default=',')
        Delimiter used in the input file (e.g., ',' for CSV or '\t' for TSV).

    write_names : bool, optional (default=True)
        If True, writes a TXT file listing all sequence names (from the first column).

    Returns:
    --------
    aligned_files : list of str
        List of file paths to the MAFFT-aligned FASTA files for each gene.
    """
    # Open the input CSV/TSV file
    with open(input_file, newline='') as file:
        reader = csv.reader(file, delimiter=delimiter)
        rows = list(reader)

    # Replace spaces in sequence names with underscores
    sequence_names = [row[0].replace(" ", "_") for row in rows[1:]]
    
    # Extract gene names from the header row (excluding first column)
    gene_names = rows[0][1:]
    
    # Extract gene accession data for each sequence (excluding first column)
    gene_columns = [row[1:] for row in rows[1:]]

    # If requested, write the sequence names to a text file
    if write_names:
        names_file = f"{output_prefix}_sequence_names.txt"
        with open(names_file, 'w') as nf:
            for name in sequence_names:
                nf.write(f"{name}\n")

    aligned_files = []  # List to store paths of aligned output files

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

                # Split accession numbers by '/' and filter out invalid entries
                accessions = [acc for acc in cell.split('/') if acc.upper() != "NA" and acc != "" and acc != "-"]
                sequences = []  # To hold the sequences retrieved from GenBank

                # Fetch each sequence from GenBank
                for acc in accessions:
                    try:
                        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
                        seq_record = SeqIO.read(handle, "fasta")
                        handle.close()
                        sequences.append(str(seq_record.seq))  # Store the sequence string
                    except Exception as e:
                        print(f"Error fetching {acc}: {e}")

                # Combine multiple sequences with 'W' delimiters to mark junctions (to be handled later)
                combined_seq = "WWWWWWWWWWWWWWW".join(sequences)
                
                # Write the sequence to the FASTA file with its name as header
                fasta_out.write(f">{seq_name}\n{combined_seq}\n")

        # Define the name for the output alignment file
        aligned_file = f"{output_prefix}_{gene_name}_aligned.fasta"

        # Run MAFFT on the generated FASTA file and save the alignment
        with open(aligned_file, 'w') as aligned_out:
            subprocess.run(["mafft", "--auto", fasta_file], stdout=aligned_out)

        # Append the aligned file path to the result list
        aligned_files.append(aligned_file)

    return aligned_files  # Return list of aligned output file paths

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
        Path to the MAFFT-aligned FASTA file to be processed.
    
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

def GB2MSA_3(alignment_dict, orphan_threshold=6, log=False):
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
    
def delete_orphan_nucleotides2(alignment, orphan_threshold, log_changes=False):
    """
    Iteratively eliminates orphan nucleotide blocks from the start and end of each sequence.
    An orphan block is a short contiguous run of nucleotides near the terminal ends separated by many gaps.

    Args:
        alignment (dict): {sequence_id: sequence_string}
        orphan_threshold (int): Max block length and max gap tolerance.
        log_changes (bool): Whether to return a log of changes made.

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

    alignment_changed = True
    while alignment_changed:
        alignment_changed = False

        for seq_id, sequence in alignment.items():
            seq_list = list(sequence)
            changed = False

            # Iteratively check from left side
            while True:
                blocks = find_blocks(seq_list)
                if len(blocks) < 2:
                    break
                first_start, first_end = blocks[0]
                next_start = blocks[1][0]
                gap_count = seq_list[first_end:next_start].count('-')

                if (first_end - first_start < orphan_threshold) and (gap_count > orphan_threshold):
                    deleted = ''.join(seq_list[first_start:first_end])
                    seq_list[first_start:first_end] = ['-'] * (first_end - first_start)
                    changed = True
                    if log_changes:
                        change_log.append(
                            f"{seq_id}: Left block {first_start}-{first_end} deleted ('{deleted}')"
                        )
                else:
                    break

            # Iteratively check from right side
            while True:
                blocks = find_blocks(seq_list)
                if len(blocks) < 2:
                    break
                last_start, last_end = blocks[-1]
                prev_end = blocks[-2][1]
                gap_count = seq_list[prev_end:last_start].count('-')

                if (last_end - last_start < orphan_threshold) and (gap_count > orphan_threshold):
                    deleted = ''.join(seq_list[last_start:last_end])
                    seq_list[last_start:last_end] = ['-'] * (last_end - last_start)
                    changed = True
                    if log_changes:
                        change_log.append(
                            f"{seq_id}: Right block {last_start}-{last_end} deleted ('{deleted}')"
                        )
                else:
                    break

            if changed:
                alignment_changed = True
                alignment[seq_id] = ''.join(seq_list)

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

    # Step 3: Log and optionally process blocks
    contiguous_invariant_blocks.sort(key=lambda x: x['length'], reverse=True)
    block_lengths = {}
    for block in contiguous_invariant_blocks:
        block_lengths.setdefault(block['length'], []).append(block)

    if partitioning_round == "max":
        add_breaks_terminal(alignment)  # Call the other function instead of inserting hashtags
    else:
        # Step 4: Track the positions for inserting hashtags
        hashtag_positions = []
        block_lengths_sorted = sorted(block_lengths.keys(), reverse=True)
        for block_length in block_lengths_sorted[:partitioning_round]:
            blocks = block_lengths[block_length]
            for block in blocks:
                start_idx = block['start']
                end_idx = start_idx + block['length'] - 1
                middle_idx = (start_idx + end_idx) // 2
                hashtag_positions.append(middle_idx)

        # Step 5: Insert hashtags
        for seq_id, seq in alignment.items():
            sorted_positions = sorted(hashtag_positions)
            shift = 0
            for middle_idx in sorted_positions:
                adjusted_idx = middle_idx + shift
                seq = seq[:adjusted_idx + 1] + '#' + seq[adjusted_idx + 1:]
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

def balanced_partitioning(alignment, log=False, partitioning_round=1):
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
        dict: Updated alignment with balanced '#' partitioning.
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

    return {names[i]: ''.join(row) for i, row in enumerate(matrix)}

##############################
# AUXILIARY FUNCTIONS TO LOG #
##############################

def compute_summary_after(alignment):
    num_seqs = len(alignment)

    # Transpose to columns
    columns = list(zip(*alignment.values()))
    
    # Count columns that contain pound signs
    total_pound = sum('#' in col for col in columns)
    
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

    summary = (
        f"Total '?' characters: {total_question_marks}\n"
        f"Total '-' characters in fully missing partitions: {total_dash_from_missing_partitions}\n"
        f"Combined total: {total_question_marks + total_dash_from_missing_partitions}\n"
    )

    return summary + "\n" + "\n".join(log_entries)

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

def UP2AP(input_fasta, output_fasta):
    """
    1. Replaces '#' with 15 'w's.
    2. Aligns using MAFFT.
    3. Replaces aligned 'w' blocks (contiguous or gapped) back to '#'.
    4. Rectifies alignment: pads partitions to ensure '#' are vertically aligned.
    5. Cleans alignment: removes columns that contain ONLY gaps (-).
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
        cmd = ["mafft", "--auto", "--quiet", temp_input_name]
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
           orphan_threshold=6):
    """
    Complete GenBank-to-MSA pipeline:
    1. Downloads sequences from GenBank and aligns them by gene using MAFFT.
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
    aligned_files = GB2MSA_1(input_file, output_prefix, delimiter=delimiter, write_names=write_names)

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

    # Step 7: Log timing if requested
    if log:
        end_wall = time.time()
        end_cpu = time.process_time()
        wall_time = end_wall - start_wall
        cpu_time = end_cpu - start_cpu

        log_file = f"{output_prefix}_log.txt"
        with open(log_file, "a") as lf:
            lf.write(f"--- GB2MSA run for '{output_prefix}' ---\n")
            lf.write(f"Wall clock time: {wall_time:.2f} seconds\n")
            lf.write(f"CPU time: {cpu_time:.2f} seconds\n\n")

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
        cmd_used = f"addSeq(alignment=..., new_seqs=..., output='{output}', write_names={write_names}, orphan_threshold={orphan_threshold}, log={log}, n2question={n2question}, gaps2question={gaps2question})"
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
        if isinstance(n2question, str) and n2question != 'all':
            target_ids = {n2question}
        elif isinstance(n2question, list):
            target_ids = set(n2question)
        elif n2question == 'all':
            target_ids = {rec.id for rec in updated_records}
        else:
            target_ids = set()

        for rec in updated_records:
            if rec.id in target_ids:
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

def prepDyn(input_file=None,
            GB_input=None,
            input_format="fasta",
            MSA=False,
            output_file=None,
            output_format="fasta",
            log=False,
            sequence_names=True,
            # Trimming parameters
            orphan_method=None,
            orphan_threshold=10,
            percentile=25,
            del_inv=True,
            # Missing data parameters
            internal_method=None,
            internal_column_ranges=None,
            internal_leaves="all",
            internal_threshold=None,
            n2question=None,
            # Partitioning parameters
            partitioning_round=0,
            partitioning_method="balanced",
            partitioning_size=None
            ):
    """
    Preprocess missing data for dynamic homology in PhyG. First, columns containing
    only gaps, orphan nucleotides, and invariant columns can be trimmed. Second,
    missing data is coded with question marks. Third, partitions are delimited in
    highly conserved regions.

    Args:
        input_file (str): Path to the input alignment file or directory. Ignored if GB_input is provided.
        GB_input (str): Path to a CSV/TSV file containing GenBank accession numbers. If provided,
                        sequences will be downloaded from GenBank and aligned before preprocessing.
        input_format (str): Format of the input alignment. Options: 'fasta' (default),
                            'clustal', 'phylip', or any format accepted by Biopython.
        output_file (str): Custom prefix for output files. If None, base_name from input_file is used.
        output_format (str): Output format. Default is 'fasta'.
        log (bool): Whether to write a log with wall-clock time. Default is False.
        sequence_names (bool): If True, writes a TXT file listing all sorted unique sequence names. Default is True.
        MSA (bool): Whether to perform MSA if input sequences specified in input_file are unaligned
        orphan_method (str): The trimming method. By default, trimming orphan nucleotides
                             is not performed. Options:
                            - 'percentile': trim using the 25th percentile;
                            - 'integer': trim with a manual threshold.
        orphan_threshold (int): Threshold used to trim orphan nucleotides if orphan_method = 'integer'.
        percentile (float): Used with orphan_method = 'percentile' to define trimming threshold.
        del_inv (bool): Whether to trim invariant terminal columns. Default is True.
        internal_method (str): Defines how to identify internal missing data. Automatic identificaton
                               of missing data is made if GB_input is provided. Otherwise, naive
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
                                   by length and '#' column(s) inserted at the midpoint of the n-largest
                                   block(s). Must define n using partitioning_round.
                                   - 'equal: Equal-length partitions are created by specifying their size
                                   or their round. If partitioning_round = 1, only 1 '#' column is
                                   inserted; if partitioning_round = 2, then 2 '#' columns are inserted.
                                   - 'max': '#' columns are inserted around blocks of missing data (every
                                   instance of '?' opening/closure but not '?' extension).
        partitioning_round (int): Number of partitioning round. Invariant regions are sorted by length
                                  in descendant order and the n-largest block(s) partitioned using '#'.
                                  If "max" is specified, pound signs are inserted arund all blocks of
                                  missing data.
        partitioning_size (int): Size of equal-length partitions if partitioning_method = 'equal'.

    Returns:
        dict: The preprocessed unaligned sequences.
    """

    # --- Start overall timer for the top-level call ---
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
        "GB_input": GB_input,
        "input_format": input_format,
        "MSA": MSA,
        "output_file": output_file,
        "output_format": output_format,
        "log": log,
        "sequence_names": sequence_names,
        "orphan_method": orphan_method,
        "orphan_threshold": orphan_threshold,
        "percentile": percentile,
        "del_inv": del_inv,
        "internal_method": internal_method,
        "internal_column_ranges": internal_column_ranges,
        "internal_leaves": internal_leaves,
        "internal_threshold": internal_threshold,
        "n2question": n2question,
        "partitioning_method": partitioning_method,
        "partitioning_round": partitioning_round,
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


    def _prepDyn_recursive(input_val, # Renamed to input_val to be more generic for file path or dict
                           GB_input,
                           input_format,
                           MSA,
                           output_file, # This is crucial: the output prefix for THIS specific call
                           output_format,
                           log,
                           sequence_names, # This will be False for recursive calls
                           _all_sequence_ids, # This is the shared set
                           _is_top_level_call, # Flag for controlling final write
                           _original_cmd_line, # Parameter to pass the original command
                           orphan_method,
                           orphan_threshold,
                           percentile,
                           del_inv,
                           internal_method,
                           internal_column_ranges,
                           internal_leaves,
                           internal_threshold,
                           n2question,
                           partitioning_round,
                           partitioning_method,
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
        if GB_input is not None:
            # (This part of the logic remains unchanged)
            print("Running GB2MSA on GenBank input...")
            
            gb_output_prefix_for_gb2msa = output_file 
            cleaned_files = GB2MSA(GB_input, output_prefix=gb_output_prefix_for_gb2msa, write_names=False, log=False)
            
            for file_path_from_gb2msa in cleaned_files:
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
                alignment = AlignIO.read(file_path_from_gb2msa, "fasta")
                alignment_dict = {record.id: str(record.seq) for record in alignment}
                _all_sequence_ids.update(alignment_dict.keys())
                _prepDyn_recursive(input_val=alignment_dict,
                                GB_input=None, input_format="dict", MSA=MSA,
                                orphan_method=orphan_method, orphan_threshold=orphan_threshold, percentile=percentile, del_inv=del_inv,
                                internal_method=internal_method, internal_column_ranges=internal_column_ranges, internal_leaves=internal_leaves, internal_threshold=internal_threshold,
                                n2question=n2question, partitioning_method=partitioning_method, partitioning_round=partitioning_round, partitioning_size=partitioning_size,
                                output_format=output_format, log=log, sequence_names=False,
                                _all_sequence_ids=_all_sequence_ids, _is_top_level_call=False, _original_cmd_line=_original_cmd_line,
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
                        print(f"Processing unaligned file for MAFFT: {file_path}")
                        temp_dir = current_output_dir if current_output_dir != "." else tempfile.gettempdir()
                        os.makedirs(temp_dir, exist_ok=True)
                        with tempfile.NamedTemporaryFile(mode="w", delete=False, dir=temp_dir, suffix=f".{input_format}") as tmp_in:
                            sequences = list(SeqIO.parse(file_path, input_format))
                            SeqIO.write(sequences, tmp_in, "fasta")
                            tmp_in_path = tmp_in.name
                        tmp_out_path = os.path.join(temp_dir, f"{os.path.basename(tmp_in_path)}_aligned.fasta")
                        try:
                            mafft_result = subprocess.run(["mafft", "--auto", tmp_in_path], capture_output=True, text=True, check=False)
                            if mafft_result.returncode != 0: raise RuntimeError(f"MAFFT alignment failed for {file_name}.")
                            with open(tmp_out_path, "w") as f_out: f_out.write(mafft_result.stdout)
                            current_file_alignment = {record.id: str(record.seq) for record in AlignIO.read(tmp_out_path, "fasta")}
                        finally:
                            if os.path.exists(tmp_in_path): os.remove(tmp_in_path)
                            if os.path.exists(tmp_out_path): os.remove(tmp_out_path)
                    else:
                        alignment_temp = AlignIO.read(file_path, input_format)
                        current_file_alignment = {record.id: str(record.seq) for record in alignment_temp}

                    if current_file_alignment:
                        _all_sequence_ids.update(current_file_alignment.keys())
                        _prepDyn_recursive(input_val=current_file_alignment,
                                GB_input=None, input_format="dict", MSA=False,
                                orphan_method=orphan_method, orphan_threshold=orphan_threshold, percentile=percentile, del_inv=del_inv,
                                internal_method=internal_method, internal_column_ranges=internal_column_ranges, internal_leaves=internal_leaves, internal_threshold=internal_threshold,
                                n2question=n2question, partitioning_method=partitioning_method, partitioning_round=partitioning_round, partitioning_size=partitioning_size,
                                output_format=output_format, log=log, sequence_names=False,
                                _all_sequence_ids=_all_sequence_ids, _is_top_level_call=False, _original_cmd_line=_original_cmd_line,
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
                    mafft_result = subprocess.run(["mafft", "--auto", tmp_in_path], capture_output=True, text=True, check=False)
                    if mafft_result.returncode != 0: raise RuntimeError(f"MAFFT alignment failed for {os.path.basename(input_val)}.")
                    with open(tmp_out_path, "w") as f_out: f_out.write(mafft_result.stdout)
                    alignment_obj = AlignIO.read(tmp_out_path, "fasta")
                finally:
                    if os.path.exists(tmp_in_path): os.remove(tmp_in_path)
                    if os.path.exists(tmp_out_path): os.remove(tmp_out_path)
            else:
                alignment_obj = AlignIO.read(input_val, input_format)
            
            alignment = {record.id: str(record.seq) for record in alignment_obj}
            _all_sequence_ids.update(alignment.keys())
        else:
            raise ValueError(f"Invalid input_val type or path: {input_val}.")


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
            orphan_threshold = calculate_orphan_threshold_from_percentile(alignment, percentile, terminal_only=True)
            if log:
                alignment, orphan_log = delete_orphan_nucleotides2(alignment, orphan_threshold, log_changes=True)
            else:
                alignment = delete_orphan_nucleotides2(alignment, orphan_threshold)
        elif orphan_method == "integer":
            if log:
                alignment, orphan_log = delete_orphan_nucleotides2(alignment, orphan_threshold, log_changes=True)
            else:
                alignment = delete_orphan_nucleotides2(alignment, orphan_threshold)


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

        if partitioning_method == "conservative" and partitioning_round > 0:
            classify_and_insert_hashtags(alignment, partitioning_round=partitioning_round)

        elif partitioning_method == "equal":
            if partitioning_round > 0:
                alignment = equal_length_partitioning(alignment=alignment, partitioning_round=partitioning_round, partitioning_size=None, log=False)
            elif partitioning_size:
                alignment = equal_length_partitioning(alignment=alignment, partitioning_round=None, partitioning_size=partitioning_size, log=False)

        elif partitioning_method == "max":
            alignment = insert_pound_around_questions(alignment)
        
        elif partitioning_method == "balanced":
            try:
                # Attempt to run the balanced partitioning. This may fail if there are not enough '?' blocks.
                processed_alignment = balanced_partitioning(
                    alignment,
                    log=False,
                    partitioning_round=partitioning_round
                )
                alignment = processed_alignment
            except (IndexError, ValueError) as e:
                # If it fails, catch the error, prepare a log message, and continue.
                warning_message = (
                    f"WARNING: 'balanced' partitioning with partitioning_round={partitioning_round} was skipped for this alignment. "
                    "This typically happens when the alignment has fewer blocks of missing data ('?') than required. "
                    "The process will continue without partitioning this file."
                )
                print(f"\n{warning_message}\n") # Print to console for immediate user feedback.
                partitioning_log_entry = warning_message # Save the message for the log file.

        refinement_question2hyphen(alignment)
        alignment = remove_columns_with_W(alignment)
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
                    columns = list(zip(*alignment.values()))
                    pound_indices = [i for i, col in enumerate(columns) if '#' in col]
                    if pound_indices:
                        log_file.write(f"Method used: {partitioning_method}")
                        if partitioning_method == "equal" and partitioning_size:
                            log_file.write(f" (partitioning_size={partitioning_size})")
                        elif partitioning_method in ["conservative", "balanced"]:
                            log_file.write(f" (partitioning_round={partitioning_round})")
                        elif partitioning_method == "max":
                            log_file.write(" (inserted at '?' block boundaries)")
                        log_file.write("\n")
                        log_file.write(f"Columns with '#' inserted: {pound_indices}\n\n")
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
                                                   GB_input=GB_input,
                                                   input_format=input_format,
                                                   MSA=MSA,
                                                   output_file=output_file,
                                                   output_format=output_format,
                                                   log=log,
                                                   sequence_names=sequence_names,
                                                   _all_sequence_ids=_all_sequence_ids_shared,
                                                   _is_top_level_call=True,
                                                   _original_cmd_line=original_cmd_line,
                                                   orphan_method=orphan_method,
                                                   orphan_threshold=orphan_threshold,
                                                   percentile=percentile,
                                                   del_inv=del_inv,
                                                   internal_method=internal_method,
                                                   internal_column_ranges=internal_column_ranges,
                                                   internal_leaves=internal_leaves,
                                                   internal_threshold=internal_threshold,
                                                   n2question=n2question,
                                                   partitioning_round=partitioning_round,
                                                   partitioning_method=partitioning_method,
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
        elif GB_input: names_file_path = "output_sequence_names.txt"
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

    if log and (GB_input is not None or (isinstance(input_file, str) and os.path.isdir(input_file))):
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
        elif GB_input: overall_log_path = "overall_prepDyn_log.txt"
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

    if GB_input or (isinstance(input_file, str) and os.path.isdir(input_file)):
        return None 
    else:
        return final_processed_alignment