# Model building

Second process in [motif discovery](../README.md). It runs after the [SEQ width check](../01_seq_width_check/README.md) passes.

[`../totalautomodelIC`](../totalautomodelIC) is the upstream recipe, version 0.22 (03 May 2018). The file in this repository is unmodified. It calls `./spacek40` and writes into the current directory. [`../run_motif.sh`](../run_motif.sh) turns the SEQ paths and the bin directory into absolute paths, links `totalautomodelIC` and [`../spacek40`](../spacek40) into `--work-dir`, changes to that directory, and runs `./totalautomodelIC`.

## What the recipe does

1. Find local-maximum k-mer seeds (`spacek40 -local`).
2. Build hit logos and position frequency matrices for each seed (`spacek40 --f`).
3. Filter by complexity and by pairwise logo similarity.
4. Refine accepted seeds. The information content of a seed is limited to between 10 and 20 bits.
5. Write accepted, unrefined, and rejected models into one SVG logo file.

## Inputs

The driver makes both SEQ paths absolute before it changes directory.

| Argument | Default in `run_motif.sh` | Role |
|----------|---------------------------|------|
| length flag | `-40N` | `-${SEQ_LENGTH}N`, passed as one argument. With `--all-gaps` that argument is `-${SEQ_LENGTH}N -allgaps` |
| background SEQ | required | Background library |
| signal SEQ | required | Enriched library |
| multinomial | `1` | Multinomial setting |
| min / max seed | `6` / `6` | Seed length range |
| length cutoff | `0.25` | Length cutoff |
| IUPAC cutoff | `-` | `-` leaves IUPAC seeds off |
| local-max cutoff | `20` | Minimum count for a local-maximum k-mer |
| logo number | `50` | Maximum number of logos |
| logo similarity | `0.1` | Logo similarity cutoff |
| output prefix | `<signal_basename>v1` | Prefix for the output files |

## Outputs

In `--work-dir`:

| File | Role |
|------|------|
| `{prefix}_logos.svg` | Accepted, unrefined, and rejected models |
| `{prefix}_refined_seeds.tmp.txt` | Refined seeds |
| `{prefix}_*.tmp.pfm`, `{prefix}_*.tmp.svg`, `{prefix}_*.tmp.txt` | Working files written by `totalautomodelIC` |

`spacek40` stays next to `totalautomodelIC` in the parent directory. `--bin-dir` points the driver at a different copy. The help text still names `scripts/02_motif/`. In this tree the directory is `workflow/02_motif_discovery/`.

The full command and a shorter hand-written `totalautomodelIC` example (maximum seed `12`, logo number `100`) are in [motif.md](../motif.md).
