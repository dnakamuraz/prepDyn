List of parameters in prepDyn scripts. Parameters can be either specified with long or short options.

## `prepDyn.py`

| **Parameter**                       | **Type**              | **Default**      | **Description** |
| ----------------------------------- | --------------------- | ---------------- | --------------- |
| `input_file` or `i`                 | `str`                 | –                | Path to an input alignment file or a directory containing multiple files. Required unless `CSV_input` is provided. |
| `CSV_input` or `csv`                | `str`                 | –                | Path to a CSV/TSV dataframe with GenBank accession numbers and/or local sequence-file paths. Required unless `input_file` is provided. Cells may contain GenBank multi-amplicons delimited by `/` or `\|`. Optionally, use `(O)` to force overlapping treatment (e.g., `MF624199(O)MF624174`) or `(N)` for non-overlapping treatment (e.g., `MF624199(N)MF624174`). |
| `input_format` or `if`              | `str`                 | `"fasta"`        | Input format, such as `fasta`, `clustal`, `phylip`, or any format accepted by Biopython. |
| `output_file` or `o`                | `str`                 | –                | Path or prefix for output file(s). |
| `output_format` or `of`             | `str`                 | `"fasta"`        | Output alignment format. |
| `log` or `l`                        | `bool`                | `True`           | Write a time log. |
| `MSA` or `msa`                      | `bool`                | `False`          | Perform a multiple sequence alignment when `input_file` contains unaligned sequences. Ignored if `CSV_input` is used. |
| `aligner` or `a`                    | `str`                 | `"mafft"`        | Multiple sequence aligner to use when alignment is performed: `mafft` or `clustalw`. |
| `sequence_names` or `s`             | `bool`                | `True`           | Write a file with sequence names. Useful for taxon sampling in POY/PhyG. |
| `multi_amplicon_min_overlap` or `mo`| `int`                 | `10`             | Minimum overlap length required to merge overlapping multi-amplicons downloaded from GenBank. Only used when overlap detection is automatic (no `(O)` marker). When `(O)` is specified in the CSV, MAFFT alignment is used instead. |
| `multi_amplicon_mismatch_rate` or `mr` | `float`            | `0.05`           | Maximum mismatch rate allowed when merging overlapping multi-amplicons from GenBank. Only used when overlap detection is automatic (no `(O)` marker). When `(O)` is specified in the CSV, MAFFT alignment is used instead. |
| `multi_amplicon_action` or `maa`    | `str`                 | `"trim"`         | How to handle overlapping multi-amplicons from GenBank: `trim` removes one overlapping copy, and `consensus` replaces the overlap with IUPAC consensus nucleotides. This applies both to automatic overlap detection and when `(O)` is used in the CSV. |
| `orphan_method` or `om`             | `str`                 | `None`           | Method to handle orphan nucleotides: `integer`, `percentile`, `adaptive_1`, `adaptive_2`, `adaptive_3`, or `None`. |
| `orphan_threshold` or `ot`          | `int`                 | `10`             | Threshold used when `orphan_method='integer'`, and also the maximum dynamic threshold when `orphan_method` is `adaptive_1`, `adaptive_2`, or `adaptive_3`. |
| `orphan_limit` or `ol`              | `float`               | `0.05`           | Modification limit as a fraction of sequence length for adaptive methods (default: 0.05 = 5%). Specifies the maximum proportion of a sequence that can be trimmed or realigned. |
| `orphan_action` or `oa`             | `str`                 | `"trim"`         | Action for orphan nucleotides: `trim` removes them, and `push` moves them adjacent to the next block iteratively. |
| `percentile` or `op`                | `float`               | `25.0`           | Percentile used to define the orphan threshold when `orphan_method='percentile'`. |
| `del_inv` or `di`                   | `bool`                | `False`          | Trim invariant terminal columns. |
| `internal_method` or `g2q`          | `str`                 | `None`           | Method to handle internal missing data: `manual`, `semi`, or `None`. |
| `internal_column_ranges` or `g2q_c` | `list` or `str`       | `"all"`          | Column ranges in Python list format for `internal_method='manual'`. |
| `internal_leaves` or `g2q_l`        | `str` or `list`       | `"all"`          | Sequence names to which internal missing-data handling is applied. Use `all` or a comma-separated list. |
| `internal_threshold` or `g2q_t`     | `int`                 | `None`           | Gap-length threshold used when `internal_method='semi'`. |
| `n2question` or `n2q`               | `str` or `list`       | `None`           | Replace IUPAC `N` with `?`. Use `all`, a single leaf name, a list of leaf names, or `None`. |
| `partitioning_method` or `pm`       | `str`                 | `"conservative"` | Partitioning method: `balanced`, `conservative`, `equal`, `max`, `None`, or `all`. `all` runs the four partitioning methods in separate directories. When `CSV_input` is used, GenBank download/alignment is performed once and reused across the batch runs, the temporary cache is deleted afterward, and a top-level overall runtime log is written. |
| `partitioning_round` or `pr`        | `int` or `str`        | `0`              | Round or number of selected blocks used by `balanced`, `conservative`, or `equal` partitioning. Ranges such as `0-10` generate one run per round in separate directories. When `CSV_input` is used in batch mode, these runs reuse the same cached aligned FASTA files and still produce a root-level batch runtime log. |
| `partitioning_conservative` or `pc` | `str`                 | `"midpoint"`     | Placement mode for conservative partitioning: `midpoint` inserts one `#` in the middle of each selected invariant block, and `flank` inserts `#` columns around each selected invariant block. |
| `partitioning_max_size` or `pms`    | `int`                 | `None`           | Initial maximum partition size. If set, the alignment is first split into equal-length partitions of this size, and the selected partitioning method is then applied independently within each partition. |
| `partitioning_size` or `ps`         | `int`                 | `None`           | Size of equal-length partitions when `partitioning_method='equal'`. |
| `keep_partitions` or `pk`           | `bool`                | `False`          | If `True`, keep existing `#` partition columns found in the input alignment and proceed with preprocessing. If `False` (default), delete any previous `#` columns before preprocessing. |

## `GB2MSA.py`

| **Parameter**                            | **Type** | **Default** | **Description** |
| ---------------------------------------- | -------- | ----------- | --------------- |
| `CSV_input` or `csv`                     | `str`    | –           | Path to a CSV/TSV file with GenBank accession numbers and/or local sequence-file paths. Each cell may contain one accession, multiple slash- or pipe-delimited GenBank accessions, one local file path, or multiple pipe-delimited local file paths for the same locus. |
| `output_prefix` or `o`                   | `str`    | –           | Path or prefix for output FASTA files. |
| `delimiter` or `d`                       | `str`    | `","`       | Delimiter used in the input file. |
| `write_names` or `w`                     | `bool`   | `True`      | Write sequence names in a separate file for taxon sampling in POY/PhyG. |
| `log` or `l`                             | `bool`   | `True`      | Write wall and CPU time to a log file. |
| `orphan_threshold` or `ot`               | `int`    | `10`        | Threshold used to clean orphan nucleotides. |
| `aligner` or `a`                         | `str`    | `"mafft"`   | Multiple sequence aligner to use: `mafft` or `clustalw`. Install ClustalW with `conda install bioconda::clustalw`. |
| `multi_amplicon_min_overlap` or `mo`     | `int`    | `10`        | Minimum overlap length required to merge overlapping multi-amplicons. Only used when overlap detection is automatic (no `(O)` marker). When `(O)` is specified in the CSV, MAFFT alignment is used instead. |
| `multi_amplicon_mismatch_rate` or `mr`   | `float`  | `0.05`      | Maximum mismatch rate allowed when merging overlapping multi-amplicons. Only used when overlap detection is automatic (no `(O)` marker). When `(O)` is specified in the CSV, MAFFT alignment is used instead. |
| `multi_amplicon_action` or `maa`         | `str`    | `"trim"`    | How to handle overlapping multi-amplicons: `trim` removes one overlapping copy, and `consensus` replaces the overlap with IUPAC consensus nucleotides. This applies both to automatic overlap detection and when `(O)` is used in the CSV. |

## `UP2AP.py`

| **Parameter**            | **Type** | **Default** | **Description** |
| ------------------------ | -------- | ----------- | --------------- |
| `input_fasta` or `i`     | `str`    | –           | Path to FASTA input containing unaligned sequences with pound signs. |
| `output_fasta` or `o`    | `str`    | –           | Path to FASTA output of aligned sequences containing pound signs. |
| `keep_unusual` or `k`    | `bool`   | `False`     | Keep pound signs (`#`) and question marks (`?`) in the output. |

## `addSeq.py`

| **Parameter**         | **Type**        | **Default** | **Description** |
| --------------------- | --------------- | ----------- | --------------- |
| `alignment` or `a`    | `str`           | –           | Path to FASTA input alignment. If question marks and pound signs are present, they are maintained. |
| `new_seqs` or `n`     | `str`           | –           | Path to FASTA input sequence(s) to be added to the alignment. |
| `output` or `o`       | `str`           | –           | Path to the output file with the new sequences aligned to the core alignment. |
| `write_names` or `w`  | `bool`          | `False`     | Write sequence names in a separate file for taxon sampling in POY/PhyG. |
| `orphan_threshold` or `ot` | `int`      | `0`         | Threshold used to detect and remove orphan DNA blocks. |
| `n2question` or `n2q` | `str` or `list` | `None`      | Replace IUPAC `N` with `?`. Use `all`, a single added leaf name, a list of added leaf names, or `None`. |
| `gaps2question` or `g2q` | `int` or `None` | `None`   | Replace contiguous gap blocks larger than this threshold with `?`. Only applied to added sequences. |
| `log` or `l`          | `bool`          | `True`      | Write a log tracking operations and runtime. |
