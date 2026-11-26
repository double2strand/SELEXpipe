# Config

Reference files for preprocessing. The scripts load the adapter FASTA when you leave `--adapter-fasta` unset. They do not load the ligand table.

| File | Role |
|------|------|
| `adapters.fa` | Example fastp adapter sequences. Default for `run_preprocess.sh` and `run_fastp.sh` |
| `ligands.tsv` | Merged construct lengths for lig147, lig200, ligN40, and ligN70 |

`run_preprocess.sh` and `run_fastp.sh` set the repository root to two levels above the script and open `config/adapters.fa` from there. With the scripts in `workflow/01_preprocessing/`, that path is this directory. Pass `--adapter-fasta` (or the fourth argument of `run_fastp.sh`) to use another file. The sequences here are examples. For custom adapters, edit `adapters.fa` or pass your own file.

`ligands.tsv` is a reference copy of the same lengths. `ligand_length` in `workflow/01_preprocessing/run_preprocess.sh` is what the driver uses. `filter_by_length.py` has no ligand names. The comment in `ligands.tsv` points at that driver.

The Sanger-ID to sample-name map is not in this directory. The example is [`workflow/00_sample_renaming/examples/example_map.tsv`](../workflow/00_sample_renaming/examples/example_map.tsv). Preprocessing does not load it.

ligN40 and ligN70 use different constant flanks, so the gap between 187 bp and 247 bp is 60 bp.
