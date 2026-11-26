# Length filter

Third process in [preprocessing](../README.md). It runs after [FASTQ to FASTA](../02_fastq_to_fasta/README.md) and before the [SEQ extract](../04_seq_extract/README.md).

[`../filter_by_length.py`](../filter_by_length.py) keeps reads inside a length window, then makes every kept read exactly the target length. The driver calls it as `python3 "${SCRIPT_DIR}/filter_by_length.py"`, so the script stays in the parent directory.

## What it does

Reads shorter than `--min_input_len` or longer than `--max_input_len` are dropped. A kept read shorter than `--target_length` is padded on the right with random A/T/C/G. A kept read longer than the target is trimmed. `--trim_from` is `right` (default), `left`, or `center`. `center` drops `extra // 2` bases from the left and keeps the next `target_length` bases.

The script's own defaults are target 154, window 145–163, trim from the right. That is the lig200 case. The script has no `--ligand` flag. [`../run_preprocess.sh`](../run_preprocess.sh) maps a ligand design to a target length and, when the bounds are omitted, uses a window of target±9.

| Ligand | Target | Window the driver uses when bounds are omitted |
|--------|--------|-------------------------------------------------|
| lig147 | 101 | 92–110 |
| lig200 | 154 | 145–163 |
| ligN40 | 187 | 178–196 |
| ligN70 | 247 | 238–256 |

Those lengths are merged construct lengths, hard-coded in `ligand_length` in the driver. ligN40 and ligN70 use different constant flanks, so the gap is 60 bp. [`config/ligands.tsv`](../../../config/ligands.tsv) repeats the table. The scripts do not read that file.

`--seed` is optional. Set it and the padding repeats. Omit it and the padding draw is unseeded. This padding is separate from `--allow-pad` on the SEQ extract.

## Inputs

| Input | Role |
|-------|------|
| `{SAMPLE}.fa` | FASTA from the previous process |
| `--target_length` | Length of every kept read. Script default 154 |
| `--min_input_len` | Shortest read to keep. Script default 145 |
| `--max_input_len` | Longest read to keep. Script default 163 |
| `--trim_from` | `right`, `left`, or `center` |
| `--seed` | Optional padding seed |

The driver passes the FASTA paths as positionals, then `--target_length`, `--min_input_len`, `--max_input_len`, `--trim_from`, and `--seed` when a seed was set.

## Outputs

| File | Role |
|------|------|
| `{SAMPLE}_{TARGET_LEN}.fa` | Reads inside the window, each exactly the target length |

The script prints:

```text
[OK] Total reads:     N
[OK] Reads kept:     N
[OK] Reads padded:   N
[OK] Reads trimmed:  N
[OK] Reads dropped:  N
```

## Run on its own

The synthetic example is 8 bp, 101 bp, and 105 bp. The script default window 145–163 drops all three records. Pin the window to the lig147 length to keep the 101 bp read. From the repository root:

```bash
python3 workflow/01_preprocessing/filter_by_length.py \
  data/examples/lig147_synthetic.fa \
  /tmp/SAMPLE_101.fa \
  --target_length 101 \
  --min_input_len 101 \
  --max_input_len 101
```

Needs Python ≥ 3.9 and Biopython. Argument table: [preprocessing.md](../preprocessing.md).

Next process: [SEQ extract](../04_seq_extract/README.md).
