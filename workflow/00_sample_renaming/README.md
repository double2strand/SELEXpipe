# Sample renaming

First stage of the pipeline, before [preprocessing](../01_preprocessing/README.md). Sequencing runs arrive under Sanger sample IDs such as `8030STDY17011914`. This stage renames those IDs to the sample and assay names used for the rest of the analysis.

Example:

```text
8030STDY17011914/8030STDY17011914_L1_R1.fastq.gz
    -> HTS/SOX2R1/SOX2R1c1/SOX2R1c1_L1_R1.fastq.gz
```

The script is [`selex_rename.py`](selex_rename.py). An example map is [`examples/example_map.tsv`](examples/example_map.tsv).

## Step 1. Build the map from the sequencing manifest

The renamer does not read the manifest. It reads a two-column map you build from it.

Take the manifest column that holds the Sanger sample ID and the column that holds the sample name, and write one sample per line:

```text
<SangerSampleID><TAB><SampleName>
```

`#` comments and blank lines are ignored. Columns are split on whitespace, so a tab or spaces both work. A line that is not two columns is skipped with a warning. A repeated Sanger ID keeps the first row and warns.

Sample names in this project look like:

| Example | Pieces |
|---------|--------|
| `C1HTS200SOX2R1` | cycle `C1`, assay `HTS`, ligand length `200`, TF `SOX2`, replicate `R1` |
| `C1NCAP200DUX4_100R2` | cycle `C1`, assay `NCAP`, ligand length `200`, TF `DUX4_100`, replicate `R2` |

The pattern is `C<cycle><ASSAY><ligand length><TF>R<replicate>`. The ligand length is the number inside the assay token (`HTS200`, `NCAP200`), the same `200` that marks that design. Derive each name from the manifest's sample-name field for that Sanger ID. The map does not invent the name from the FASTQ filename.

`examples/example_map.tsv` is a reference map in this shape: 768 rows, 384 `HTS200` and 384 `NCAP200`, including `8030STDY17011914	C1HTS200SOX2R1`.

## Step 2. Run the renamer

From the repository root, preview first:

```bash
python3 workflow/00_sample_renaming/selex_rename.py \
  --base-dir <run_dir> \
  --map workflow/00_sample_renaming/examples/example_map.tsv \
  --dry-run
```

`--dry-run` prints the planned renames and does not create or move anything. Run the same command without `--dry-run` to apply them.

`<run_dir>` is the folder whose children are the Sanger ID directories (`8030STDY17011914/`, and so on). Each of those holds the lane FASTQs.

## Final layout

With `--assay-subdir` (the default):

```text
<base>/HTS|NCAP/<TF+replicate>/<TF+replicate>c<cycle>/<same>_L<lane>_R<read>.fastq.gz
```

`8030STDY17011914` with name `C1HTS200SOX2R1` becomes:

```text
<base>/HTS/SOX2R1/SOX2R1c1/SOX2R1c1_L1_R1.fastq.gz
```

The sample replicate (`R1` in `SOX2R1`) and the sequencing read (`_R1` / `_R2`) both appear in the filename. Read 1 and read 2 are the suffixes `_R1` and `_R2`. The replicate is part of the sample name that comes from the map, not that read suffix.

`--no-assay-subdir` drops the `HTS` or `NCAP` directory. An HT-SELEX sample and an NCAP-SELEX sample with the same TF, replicate, and cycle then land on the same folder and the run stops as a clash. Leave the assay directory on unless those pairs cannot occur.

## Options

`python3 workflow/00_sample_renaming/selex_rename.py --help`

| Option | Role |
|--------|------|
| `--map` | Required. Two columns, TAG and NAME (for example `8030STDY17011914 C1HTS200SOX2R1`) |
| `--base-dir` | Required. Folder containing the Sanger-ID folders. Cycle-pattern folders under it are searched too |
| `--assay-pattern` | Middle pattern between the cycle and the sample name. Repeatable. Default `HTS200`, `NCAP200` |
| `--assay-subdir` / `--no-assay-subdir` | Put outputs under `HTS/` and `NCAP/` so the two assays never collide. Default on |
| `--no-recursive-file-rename` | Rename files only directly inside each folder, not in sub-folders |
| `--ids-only` | Only Sanger-ID folders and folders named by the map, not other folders that merely match the cycle pattern |
| `--skip-existing` | Skip items whose target already exists or clashes instead of aborting |
| `--dry-run` | Show what would be done without changing anything |

A trailing lowercase `hts` is stripped from the sample name (`MAFG_2147c4hts` becomes `MAFG_2147c4`). Other characters outside `[A-Za-z0-9._-]` become `_`, runs of `_` collapse, and `_` is stripped from the ends. The name is then parsed as `C<cycle>` + assay pattern + the rest. The folder name is `<TF+replicate>c<cycle>` with the cycle in lower case (`SOX2R1c1`). A name that does not match that pattern is used as-is. Two map rows with the same name get `_2`, `_3`, and so on, in map order.

Only files named `<TAG>_*.fastq.gz` or `<TAG>_*.fq.gz` are renamed from the Sanger ID, and the tag match is case-insensitive. Other files in the folder move with it and keep their names, unless the name contains the old sample or folder token as its own word, in which case that token is replaced.

Every move is checked before anything is touched. A clash (target already there, two items aiming at one path, a file blocking a folder) changes nothing and exits 2. The script never overwrites a file. A second run skips work that is already in the final layout, and it will finish FASTQs left inside a folder that was only partly renamed.

| Exit | Meaning |
|------|---------|
| 0 | Finished, including a dry run |
| 1 | `--base-dir` is not a directory, or `--map` is not a file |
| 2 | A clash was found before any change, or a target appeared during the run |

Preprocessing reads the renamed FASTQs. For `SOX2R1c1`, `--input-dir` is the cycle folder and `--sample` is `SOX2R1c1`. See [preprocessing](../01_preprocessing/README.md).
