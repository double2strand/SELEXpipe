# SEQ extract

Fourth process in [preprocessing](../README.md), and the last one. It runs after the [length filter](../03_length_filter/README.md).

[`../trim_to_kmer.py`](../trim_to_kmer.py) writes one fixed-width window per FASTA record into a SEQ file: one sequence per line, no header. The default width is 40 bp. The driver calls `python3 "${SCRIPT_DIR}/trim_to_kmer.py"`, so the script stays in the parent directory.

## What it does

One window per molecule. The default window is `--start 0` and `--kmer-size 40`. `--center` centers the window instead of using `--start`. `--length` is an alias of `--kmer-size`.

A window that does not fit is dropped. Random bases are written only when `--allow-pad` and `--seed` are both set. A read the length filter has already made 101 bp or 154 bp is long enough for the default 40-mer, so this step writes the window from the read. Filter padding sits on the right, so a 5′ 40-mer (`--start 0`) is original sequence when the retained read was at least 40 bp.

## Inputs

| Input | Role |
|-------|------|
| `{SAMPLE}_{TARGET_LEN}.fa` | Length-filtered FASTA |
| `--kmer-size` | Window length. Default 40 |
| `--start` or `--center` | Where the window sits. The driver passes `--start` unless `--center-kmer` was set |
| `--allow-pad`, `--seed` | Optional. Both are required to pad a short window |

## Outputs

| File | Role |
|------|------|
| `{SAMPLE}_{TARGET_LEN}.seq` | One line per kept read. Motif discovery reads this file |

From the repository root, continuing the lig147 synthetic example:

```bash
python3 workflow/01_preprocessing/trim_to_kmer.py \
  /tmp/SAMPLE_101.fa \
  /tmp/SAMPLE_101.seq \
  --kmer-size 40
```

The 101 bp synthetic record is `A`×40 + `C`×40 + `G`×21, so the SEQ line is forty `A`s. The 8 bp and 105 bp records are outside a 101–101 window and are dropped by the length filter before this step.

Details and the centered and `--allow-pad` forms: [preprocessing.md](../preprocessing.md).

Next stage: [motif discovery](../../02_motif_discovery/README.md).
