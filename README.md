# prepDyn: Preprocessing sequences for dynamic homology

[![language](https://img.shields.io/badge/language-python-green?style=flat&logo=python&logoColor=green)](https://www.python.org)
[![author](https://img.shields.io/badge/author-DYM_Nakamura-green?logo=googlescholar&logoColor=green)](https://scholar.google.com/citations?user=c0W8Cm8AAAAJ&hl=en)
[![license](https://img.shields.io/badge/license-GPL_v3-green?logo=gnu&logoColor=green)](https://www.gnu.org/licenses/gpl-3.0.html)

A collection of Python scripts to facilitate the preprocessing of input sequences for dynamic homology. 

In dynamic homology, data should be preprocessed to distinguish differences in sequence length resulting from missing data or insertion-deletion events to avoid grouping from artifacts. However, previous empirical studies using POY/PhyG manually preprocessed data with varying approaches. Here we present **prepDyn**, a collection of Python scripts to facilitate the preprocessing of input sequences to POY/PhyG. **prepDyn** comprises four steps: (1) data collection from GenBank, (2) trimming, (3) identification of missing data, and (4) partitioning.

Copyright (C) Daniel Y. M. Nakamura 2025

## Installation

The two dependencies that should be installed beforehand by the user are:
- Python v. 3.10.9 (or newer), including *argparse*, *ast*, *csv*, *importlib*, *re*, *StringIO*, *subprocess*, *sys*, *tempfile*, and *time*, which are usually part of recent versions of Python.
- MAFFT v. 7.5.2 (or newer), installed in $PATH as 'mafft'.

```
conda create -n new_env python=3.10 --yes
conda install bioconda::mafft
```

Other dependencies are Python modules that will be automatically installed by **prepDyn** when you run it for the first time:
- Bio v. 1.73 (or newer), including *AlignIO*, *Entrez*, *SeqIO*, *Align*, *Seq*, and *SeqRecord*.
- matplotlib v. 3.7.0 (or newer)
- numpy v. 1.23.5 (or newer)
- termcolor

If the  modules are not installed automatically, try:

```
conda install conda-forge::biopython
conda install conda-forge::matplotlib
conda install anaconda::numpy
conda install conda-forge::termcolor
```

Finally, clone the **prepDyn** repository using the command:

```
git clone https://github.com/danimelsz/PrepDyn.git
```

## Usage

**prepDyn** is organized in three Python files in the directory src:
- prepDyn.py: main script integrating the pipeline.
- GB2MSA.py: script to download sequences from GenBank and identify internal missing data.
- addSeq.py: script to align one or a few sequence(s) to a previously preprocessed (profile) alignment.

A summary of parameters used in prepDyn.py:

| **Parameter**            | **Type**              | **Default**  | **Description**                                                                   |
| ------------------------ | --------------------- | ------------ | --------------------------------------------------------------------------------- |
| `input_file`             | `str`                 | –            | Path to an alignment file or directory of alignments.                             |
| `GB_input`               | `str`                 | –            | Path to a CSV/TSV of GenBank accessions. Overrides `input_file`.                  |
| `input_format`           | `str`                 | `"fasta"`    | Format of input file(s) (e.g., `"phylip"`, `"clustal"`).                          |
| `MSA`                    | `bool`                | `False`      | If `True`, perform multiple sequence alignment on unaligned input sequences.      |
| `output_file`            | `str`                 | –            | Custom prefix for output files.                                                   |
| `output_format`          | `str`                 | `"fasta"`    | Format for the output alignment.                                                  |
| `log`                    | `bool`                | `False`      | Write a detailed log file.                                                        |
| `sequence_names`         | `bool`                | `True`       | Write a file with all unique sequence names.                                      |
| `orphan_method`          | `str`                 | `None`       | Method to trim orphan nucleotides: `'percentile'`, `'integer'`, or `None`.                 |
| `orphan_threshold`       | `int`                 | `10`         | Manual length threshold for `orphan_method='integer'`.                               |
| `percentile`             | `float`               | `25`         | Percentile for `orphan_method='auto'`.                                            |
| `del_inv`                | `bool`                | `True`       | Trim invariant columns from alignment ends.                                       |
| `internal_method`        | `str`                 | `None`       | Method to replace internal gaps with `?`: `'manual'`, `'semi'`, or `None`.        |
| `internal_column_ranges` | `list`                | –            | Column ranges (e.g., `[[10, 20]]`) for `'manual'` method.                         |
| `internal_leaves`        | `str` or `list`       | `"all"`      | Sequences to apply internal gap replacement to.                                   |
| `internal_threshold`     | `int`                 | –            | Gap length threshold for `'semi'` method.                                         |
| `n2question`             | `str`, `list`, `None` | `None`       | Replace ambiguous `'N'` with `?`. Options: `'all'`, list of names, or `None`.     |
| `partitioning_method`    | `str`                 | `"balanced"` | Method to insert `#` markers: `'balanced'`, `'conservative'`, `'equal'`, `'max'`. |
| `partitioning_round`     | `int`                 | `0`          | Number of partitions/rounds for relevant partitioning methods.                    |
| `partitioning_size`      | `int`                 | –            | Partition size for `partitioning_method='equal'`.                                 |

Parameters can be either specified with long or short options. For more information:

```
python src/prepDyn.py -h
python src/GB2MSA.py -h
python src/addSeq.py -h
```

The following examples are designed for users with little experience on Unix. If you have questions, send a message using **GitHub issues**. Do not move the scripts from the directory *src*, otherwise the modular structure will break.

### Example 1: Basic

The basic use of **prepDyn** is running all four steps using a single command. Given an input CSV, whose first column is called *Terminals* and the other columns are the names of genes (each cell containing the correspondent GenBank accession number), the following command will download sequences, trim invariants and orphan nucleotides <10 bp in terminal positions, and identify missing data as *?* (all differences in sequence length in terminal positions are missing data). In the CSV file, if more than one GenBank accession number is specified in the same cell refering to non-overlapping fragments of the same gene (e.g. MT893619/MT895696), the space between them is automatically identified as internal missing data (?).

```
python src/prepDyn.py \
    --GB_input test_data/tutorial/ex1.1/ex1.1_input.csv \
    --output_file test_data/tutorial/ex1.1/ex1.1 \
    --del_inv T \
    --orphan_method integer \
    --orphan_threshold 10 \
    --partitioning_method None \
    --log T 
```

We specified *--paritioning_round 0*, which means that partitioning was not performed. As a heuristic, we recommend testing the impact of adding pound signs to the tree optimality scores using a successive partitioning strategy. For instance, if you specify *partitioning_method conservative* and *--partitioning_round 1*, the largest block(s) of contiguous invariants will be partitioned.

```
python src/prepDyn.py \
    --input_file test_data/tutorial/ex1.2/ex1.2_input.fasta \
    --output_file test_data/tutorial/ex1.2/ex1.2 \
    --partitioning_method balanced \
    --partitioning_round 1 \
    --log T
```

This process can continue until tree costs reported by POY/PhyG remain stationary (e.g. *--partitioning_round 2* inserts pound signs in the 1- and 2-largest block(s) of contiguous invariants). Other methods of partitioning are also available and the user should explore whether they can reduce tree costs.

### Example 2: Multiple alignments

Suppose you have a phylogenomic dataset with hundreds of gene alignmens in the directory *./data/*. Phylogenomic datasets are usually unavailable in GenBank, but are available in repositories like Dryad and Zenodo. You can preprocess all unaligned gene alignments in FASTA format using a single command:

```
python src/prepDyn.py \
    --input_file test_data/tutorial/ex3.1/ \
    --input_format fasta \
    --output_file test_data/tutorial/ex3.1/out \
    --MSA T \
    --del_inv T \
    --orphan_method integer --orphan_threshold 10 \
    --internal_method semi --internal_threshold 15 \
    --partitioning_method max
```

If the input files are already aligned, just change the boolean parameter MSA to False:

```
python src/prepDyn.py \
    --input_file test_data/tutorial/ex3.2/ \
    --input_format fasta \
    --output_file test_data/tutorial/ex3.2/ \
    --MSA F \
    --del_inv T \
    --orphan_method integer --orphan_threshold 10 \
    --internal_method semi --internal_threshold 15 \
    --log T
```

### Example 3: Appending new sequences

 MAFFT is unable to align sequences if pound signs or question marks are present. This is a problem when we try to align new sequences to a prevously preprocessed profile alignment. To avoid manual alignment by eye, addSeq.py allows aligning new sequences to profile alignments. Gaps, missing data, and pound signs are not modified for the sequences present in the profile alignment. Gaps, missing data, and pound signs are only inserted in the new sequences.

A simple example:

```
python src/addSeq.py \
    --alignment test_data/tutorial/ex4.1/ex4.1_aln.fas \
    --new_seqs test_data/tutorial/ex4.1/ex4.1_new_seqs.fas \
    --output test_data/tutorial/ex4.1/ex4.1_out.fas \
    --log True
```

A more complex example, where new sequences were preprocessed using trimming of blocks of orphan nucleotides of length lesser than 45 bp, replacement of internal blocks of gaps longer than 20 with question marks, and replacement of all IUPAC N with question marks in the sequence *Thoropa_miliaris_CFBH10125*:

```
python src/addSeq.py \
    --alignment test_data/tutorial/ex4.2/ex4.2_aln.fas \
    --new_seqs test_data/tutorial/ex4.2/ex4.2_new_seqs.fas \
    --output test_data/tutorial/ex4.2/ex4.2_out.fas \
    --orphan_threshold 45 \
    --gaps2question 20 \
    --n2question Thoropa_miliaris_CFBH10125 \
    --write_names True \
    --log True
```

Warning: The input *--new_seqs* cannot be longer than the input profile *--alignment*.

## Cite

If you use **prepDyn** in your research, cite this repository.