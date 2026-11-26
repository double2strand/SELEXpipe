# Motif discovery

Builds motif models from a background SEQ library and a signal (enriched) SEQ library with `totalautomodelIC` 0.22 and the `spacek40` binary shipped in this directory.

[`run_motif.sh`](run_motif.sh) is the driver. [`motif.md`](motif.md) has the argument table and the commands written out.

## Processes

| Order | Process | What runs |
|-------|---------|-----------|
| 1 | [SEQ width check](01_seq_width_check/README.md) | `validate_seq_width` inside `run_motif.sh` |
| 2 | [Model building](02_model_building/README.md) | `./totalautomodelIC`, which calls `./spacek40` |

## Scripts in this directory

| File | Role |
|------|------|
| `run_motif.sh` | Width check, then `totalautomodelIC` |
| `totalautomodelIC` | Upstream recipe, version 0.22 (03 May 2018). Unmodified |
| `spacek40` | ELF binary. `totalautomodelIC` runs it as `./spacek40` |

`run_motif.sh` defaults `--bin-dir` to its own directory and expects `totalautomodelIC` and `spacek40` beside it. The recipe writes relative to the current directory and calls `./spacek40`, so the driver symlinks both files into `--work-dir` and runs there. The numbered folders describe the two processes. The three files stay in this directory so that lookup still works.

## Inputs

SEQ files from [preprocessing](../01_preprocessing/README.md): one sequence per line, every line the same width. The default width is 40 (`--seq-length`).

| Input | Meaning |
|-------|---------|
| `--background` | Naïve or input library, for example `processed/BACKGROUND/BACKGROUND_101.seq` |
| `--signal` | An enriched round, for example `processed/SIGNAL_r4/SIGNAL_r4_101.seq` |

Run preprocessing once per library to produce those files.

## Outputs

Written in `--work-dir`:

| File | Meaning |
|------|---------|
| `{prefix}_logos.svg` | Accepted, unrefined, and rejected models |
| `{prefix}_refined_seeds.tmp.txt` | Refined seeds |
| `{prefix}_refined_model*.tmp.pfm` and matching `.svg` | Per-model files from `totalautomodelIC` |

If `--output-prefix` is omitted, the prefix is the signal filename without `.seq`, plus `v1`.

## Run

From the repository root:

```bash
workflow/02_motif_discovery/run_motif.sh \
  --background ./processed/BACKGROUND/BACKGROUND_101.seq \
  --signal ./processed/SIGNAL_r4/SIGNAL_r4_101.seq \
  --output-prefix SIGNAL_r4v1 \
  --work-dir ./motif_out/SIGNAL_r4
```

The `--help` text still prints `scripts/02_motif/run_motif.sh` and describes `--bin-dir` as `scripts/02_motif/`. Those strings are inside the script. The command above is the path in this tree.

Defaults: multinomial `1`, seed length `6` to `6`, length cutoff `0.25`, IUPAC cutoff `-`, local-max cutoff `20`, logo number `50`, similarity `0.1`. `--all-gaps` passes the length flag as the single argument `-${SEQ_LENGTH}N -allgaps`.

Needs `bash`, `python3` for the width check, and `bc`, `awk`, `sed`, and `grep` for `totalautomodelIC`. `spacek40` is the binary in this directory. It is not a conda package.
