# Preprocessing dynamic characters <img src="figures/logo.jpg" align="right" width="120">


[![language](https://img.shields.io/badge/language-python-green?style=flat&logo=python&logoColor=green)](https://www.python.org)
[![author](https://img.shields.io/badge/author-DYM_Nakamura-green?logo=googlescholar&logoColor=green)](https://scholar.google.com/citations?user=c0W8Cm8AAAAJ&hl=en)
[![license](https://img.shields.io/badge/license-GPL_v3-green?logo=gnu&logoColor=green)](https://www.gnu.org/licenses/gpl-3.0.html)

**prepDyn** is a collection of Python scripts to facilitate the preprocessing of input DNA sequences for dynamic homology. 

In dynamic homology, data should be preprocessed to distinguish differences in sequence length resulting from missing data or insertion-deletion events to avoid grouping from artifacts. However, previous empirical studies using POY and PhyG manually preprocessed data with varying approaches. Here we present **prepDyn**, a collection of Python scripts to facilitate the preprocessing of input sequences to POY/PhyG. The main script `prepDyn.py` comprises four steps: (1) data collection from GenBank, (2) trimming, (3) identification of missing data, and (4) partitioning.

Copyright (C) Daniel Y. M. Nakamura 2026

## Installation

### Conda

The easiest way to install **prepDyn** is creating a new environment in Conda:

```bash
# Create environment
conda create -n prepdyn -c dnakamuraz -c bioconda -c conda-forge prepdyn

# Activate
conda activate prepdyn

# Test
python -m prepDyn
```

Keep it updated:

```bash
conda deactivate prepdyn
conda update -c dnakamuraz prepdyn
conda activate prepdyn
```

### Docker

You can also run **prepDyn** inside a container with Docker. First, open the Docker software. Then, pull the published image from GitHub Container Registry (GHCR):

```bash
docker pull ghcr.io/dnakamuraz/prepdyn:latest
```

If you want to keep it updated, just pull it again. 

### Manual installation

If you prefer to install all dependencies manually, the external dependencies that should be installed are:
- Python v. 3.10.9 (or newer), including *argparse*, *ast*, *csv*, *importlib*, *re*, *StringIO*, *subprocess*, *sys*, *tempfile*, and *time*, which are usually part of recent versions of Python.
- MAFFT v. 7.5.2 (or newer), installed in $PATH as 'mafft'. This is the default aligner.
- Optionally, ClustalW installed in $PATH as 'clustalw', if you want to run `--aligner clustalw`.

```
# Create a conda environmnt called 'prepdyn'
conda create -n prepdyn python=3.10 --yes

# Inside the newly created environment, install 'mafft'
conda activate prepdyn
conda install bioconda::mafft

# Optional: install ClustalW
conda install bioconda::clustalw
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
git clone https://github.com/dnakamuraz/prepDyn.git
```

## Usage

<img src="figures/fig4_workflow.jpg" align="center" width="1200">

**prepDyn** is organized in four stand-alone Python scripts in the directory src:

| Script        | Description                                                                 |
|----------------|-----------------------------------------------------------------------------|
| `prepDyn.py`   | The main script integrating the pipeline.                                   |
| `GB2MSA.py`    | Downloads sequences from GenBank and identifies internal missing data.     |
| `addSeq.py`    | Aligns one or a few sequence(s) to a previously preprocessed alignment.    |
| `UP2AP.py`     | Aligns sequences containing pound signs.                                    |

The main script is `prepDyn.py`, which comprises four steps: 

- (1) Data collection: Based on a CSV dataframe containing GenBank accession numbers or FASTA sequences in a local directory

- (2) Trimming: Deletion of flanking invariants and orphan nucleotides

- (3) Identification of missing data: Internal missing data identified in the first step or specified by the user; flanking gaps are automatically corrected to missing characters

- (4) Successive partitioning: Pound signs are inserted sucessively until tree costs stabilize. Position of pound signs defined by partitioning strategies (balanced, conservative, equal-length, and maximum), which are competitive via tree costs. Recommended for large datasets.

A summary of parameters available in `prepDyn.py` is [here](tables/parameters.md).

## Citation

If you use **prepDyn** in your research, cite this repository.

## Tutorial

Check the [**Wiki**](https://github.com/dnakamuraz/prepDyn/wiki) page for tutorials using the Python scripts or the Docker image. If you have questions, send a message using **GitHub issues**.

## FAQ

<details>
<summary>What is prepDyn used for?</summary>

**prepDyn** is used to preprocess DNA sequences for dynamic homology in **POY/PhyG**,
including trimming orphan nucleotides, handling missing data,
and generating partitions.

</details>

<details>
<summary>What is dynamic homology?</summary>

Two main strategies are used to align DNA sequences and determine homology. In static homology, sequences are first aligned using a similarity-based (= phenetic) multiple sequence alignment, and that fixed alignment is then evaluated across trees to find the lowest-cost tree. In dynamic homology, unaligned sequences are optimized directly on each tree during the search, and homologous nucleotides are inferred afterward from the best tree via implied alignments. Thus, dynamic homology accounts for alignment uncertainty, whereas static homology relies on a single possible alignment. Furthermore, dynamic homology uses the same optimality criterion throughout all steps, whereas static homology does not.

Finding the lowest-cost tree from unaligned sequences is NP-hard, so empirical analyses rely on heuristic methods (e.g. direct optimization and iterative-pass). Dynamic homology frequently finds more optimal hypotheses than static homology at the cost of computational resources (runtime and memory).

</details>

<details>
<summary>What do the different output files mean?</summary>

**prepDyn** generates multiple output files for each gene or dataset. Understanding these files is essential for downstream analyses:

- **output_*.fasta**: Unaligned sequences. These are the raw sequences extracted from GenBank or provided as input, before any alignment or preprocessing. Use these if you plan to perform your own alignment.

- **output_*_aligned_GB2MSA.fasta**: Aligned but not preprocessed sequences. These sequences have been aligned using the selected aligner during the GB2MSA step, but no preprocessing has been applied yet. They may contain flanking gaps, orphan nucleotides, or internal missing data that need to be addressed.

- **output_*_aligned_preprocessed.fasta**: Fully preprocessed sequences. These are the final output files after all preprocessing steps have been completed, including trimming of flanking invariants, removal/realignment of orphan nucleotides, correction of internal missing data, and identification of gaps. These files are ready for downstream phylogenetic analyses in POY/PhyG or other programs.

- **output_*_log.txt**: Log file containing a summary of preprocessing steps applied to the alignment, including the number of sequences, alignment length before and after each step, and any warnings or issues encountered.

</details>

<details>
<summary>In addition to POY/PhyG, can prepDyn be used to preprocess input data for other phylogenetic programs?</summary>

Yes. If partitioning is skipped, **prepDyn** can still preprocess DNA sequences for other phylogenetic programs. The Step 3 is particularly useful to identify missing data and avoid downstream problems in software that treat gaps as a fifth character-state (e.g., TNT).

</details>

<details>
<summary>My sequences are too long and POY/PhyG is unable to start phylogenetic analyses. What should I do?</summary>

Dynamic homology implemented in **POY/PhyG** is NP-hard. If sequences are too long, use the parameter `partitioning_max_size`, so that sequences are initially split into equal-length partitions of size X before applying the partitioning methods (balanced, conservative, equal, maximum) in each resulting chunk.

</details>

<details>
<summary>What are multi-amplicons?</summary>

In many cases, sequences are available for a few genes but not others. In this case, missing data should be indicated with the string "NA" in the CSV cell (or empty cells). CSV cells may contain either GenBank accession numbers or local sequence-file paths (absolute or relative), and the same table can mix both sources. For multi-amplicons, GenBank accession numbers may be delimited with `/` or `|`, while local file paths should be delimited with `|`. These multi-amplicons can be classified as non-overlapping or overlapping. If overlapping, the overlapping nucleotides from one of the amplicons are deleted or consensus is called. If non-overlapping, the intersequence dashes are replaced with question marks (internal missing data).

It is also possible to force the classification of multi-amplicons as **overlapping** by replacing the standard delimiter with `(O)` (e.g., `MF624199(O)MF624174` for GenBank accessions or `file1.fas(O)file2.fas` for local files) or as **non-overlapping** by using `(N)` (e.g., `MF624199(N)MF624174` or `file1.fas(N)file2.fas`). Without these markers, classification is automatic based on sequence overlap detection.

When using the `(O)` marker:
- The sequences are **automatically aligned using MAFFT**
- The alignment is used to identify the overlapping region
- Overlapping nucleotides are deleted (if `-maa trim`) or consensus is called (if `-maa consensus`)
- This approach is robust to sequence divergence and handles overlaps that automatic detection might miss

Without using markers:
- Overlap detection is automatic using the parameters `-mr` (mismatch rate) and `-mo` (minimum overlap length)
- Classification depends on detecting sequence overlaps in the raw (unaligned) sequences

</details>

<details>
<summary>Which partitioning strategies are available?</summary>

Balanced: Inserts pound signs (`#`) around all blocks of missing data (similar to maximum partitioning) and then merge two adjacent blocks when smaller than the n-largest block. The value of *n* is defined using the `-pr` (`partitioning_round`) parameter.

Conservative: Given blocks of contiguous invariants sorted by length, partitions the *n*-largest block(s). The value of *n* is defined using the `-pr` (`partitioning_round`) parameter. Pound signs can be inserted in the middle or flanking the blocks.

Equal-length: Inserts pound signs (`#`) to divide the alignment into equal-length partitions. The size is defined using `-ps` (`partitioning_size`) or the round using `-pr` (`partitioning_round`).

Maximum: Inserts pound signs (`#`) columns around all blocks of contiguous missing data (i.e., before and after every instance of `?` opening/closure).

</details>

<details>
<summary>What is the best partitioning strategy?</summary>

The best partitioning strategy is dataset-dependent and the user must test it empirically. Empirical analyses indicate that conservative, equal-length, and maximum partitioning perform better.

</details>

<details>
<summary>Can I run several partitioning strategies or rounds at once?</summary>

Yes. `partitioning_method` accepts `all`, which runs `conservative`, `balanced`, `max`, and `equal` in separate directories. `partitioning_round` also accepts ranges such as `0-10`, which generate one run per round. For example, `prepDyn.py -pm all -pr 0-2 -o output` creates one directory per method under `output`, and within each method directory creates subdirectories for `round_0`, `round_1`, and `round_2`. If `CSV_input` is used in one of these batch runs, GenBank sequences are downloaded and aligned only once, and all subsequent batch runs reuse the cached aligned FASTA files. The temporary `_gb2msa_cache` directory is deleted automatically after the batch finishes, and a top-level overall runtime log is written in the root output directory.

</details>

<details>
<summary>What are orphan nucleotides and which strategy is better to handle them?</summary>

Orphan nucleotides, defined as short stretches of nucleotides that appear separated from the main sequence block by long runs of gaps. These nucleotides can be artifacts from sequencing errors or alignment errors. As such ,orphan nucleotides should be either trimmed (`orphan_action trim`) or realigned (`orphan_action push`). 

Operationally, orphan nucleotides can be identified as contiguous nucleotide segments shorter than a user-defined threshold x, located at the flanks of a sequence and separated from the nearest substantial nucleotide block by gap regions longer than x. Because the optimal value of x depends on the characteristics of the dataset, it should be specified by the user (orphan_threshold) via visual inspection of alignment. Based on our experience, values between 10 and 30 generally perform well across many datasets. See also automatic methods available in `orphan_method`.

</details>

<details>
<summary>Can I detect orphan nucleotides automatically?</summary>

When a single orphan threshold is not feasible or visual inspection is too laborious in large datasets, adaptive orphan threshold methods are available using `orphan_method adaptive_1`, `orphan_method adaptive_2`, or `orphan_method adaptive_3`. These methods use an iterative approach where the threshold is updated dynamically. The key difference is:

- **`adaptive_1`**: Uses block size and modification-budget criteria only.
- **`adaptive_2`**: Uses block size, modification-budget, and flanking-gap criteria.
- **`adaptive_3`**: Uses block size, modification-budget, flanking-gap, and shared-string criteria.

**How adaptive methods work:**

a. **Budgeting**: Before doing anything, it calculates a fraction of the length of each sequence (defined by `orphan_limit`, default 5%), which is the maximum percentage of sequence allowed to be trimmed or realigned.

b. **Starting threshold**: Instead of immediately using the user-provided `orphan_threshold`, the adaptive methods start their dynamic threshold at 1. The user-defined `orphan_threshold` acts as the maximum value of this dynamic threshold.

c. **Iterative growth**: It enters a loop where it looks strictly at the outermost left and outermost right contiguous blocks of nucleotides for every sequence. For a block to be considered an orphan, it must meet two to four conditions depending on the adaptive method: (1) the length of the block must be less than or equal to the current dynamic threshold, and (2) modifying this block must not cause the sequence to exceed its modification budget (defined by `orphan_limit`).

   For **`adaptive_2`** and **`adaptive_3`**, there is an additional condition: (3) the gap regions flanking the block must be longer than the dynamic threshold.

   For **`adaptive_3`** only, there is an additional condition: (4) other orphan blocks sharing exactly the same sequence string cannot occur at the same position in other sequences (shared-string check).

   If these conditions are met, the `orphan_action` is conducted (either trimming or realignment). If a change was made anywhere in the alignment, the script resets the dynamic threshold back to 1. This is because trimming or pushing an outer block exposes a new outer block, which might be a tiny 1-nucleotide orphan. If no changes were made, it increments the dynamic_threshold by 1 and scans the alignment again.

</details>

<details>
<summary>Can I specify more than one gene in the same column in the CSV file?</summary>

We recommend only specifying a single gene in each column in the CSV input file. However, some systematists are used to treat 12S + tRNAVal + 16S as a single fragment called H1. If analyzing H1 amplified with more than one set of primers (multi-amplicons), you should always list 12S, tRNAVal, and 16S in the same order (e.g., "12S_accession/16S_accession") for all sequences in the column refering to H1 in CSV. If some sequences have them as "12S/16S" while others have "16S/12S", the sequences will be misaligned, leading to incorrect orthology assignments and biased phylogenetic results.

</details>

<details>
<summary>Why should I annotate mitogenomes?</summary>

Mitochondrial genomes, especially in organisms with rearranged gene orders, can cause serious orthology assignment problems if not properly annotated. When sequences from multiple species are downloaded from GenBank and aligned, gene order rearrangements can lead to misalignment of non-orthologous sequences. Therefore, before running **prepDyn**, annotate mitogenomes and split genes that will be used as input (i.e. one input FASTA file for each gene specified in CSV as local files instead of specifying a GenBank accession number of a whole mitogenome). Only include sequences from the same gene in the same column in the input CSV file.

Proper annotation ensures that orthologous sequences are correctly identified and aligned, resulting in more reliable phylogenetic inferences.

</details>
