# SEQ width check

First process in [motif discovery](../README.md). It runs before [model building](../02_model_building/README.md), and before the driver looks for `totalautomodelIC` or `spacek40`.

The check is `validate_seq_width` inside [`../run_motif.sh`](../run_motif.sh). It is a Python heredoc in that script, so there is no separate file for this process.

## What it does

The driver calls the check on the background file, then on the signal file.

- `--seq-length` has to be a positive integer. The default is 40.
- The file has to contain at least one sequence line.
- After the newline is stripped, every line has to be non-empty and exactly that length.

The process exits on any of those failures. An empty library or a mixed-width library does not reach `totalautomodelIC`.

## Inputs

| Input | Role |
|-------|------|
| `--background` | Background SEQ path |
| `--signal` | Signal SEQ path |
| `--seq-length` | Expected line width. Default 40 |

`python3` has to be on `PATH`.

## Outputs

No output file. A passing check continues to model building. A failure prints a `[motif] error:` line on stderr and exits non-zero.

The same check, written out, is in [motif.md](../motif.md).
