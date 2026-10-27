# Example data

`lig147_synthetic.fa` is a three-record FASTA made for the length filter and
the one-window extract. It is not a sequencing run.

| Record | Length | With `--target_length 101` and window 101–101 |
|--------|--------|-----------------------------------------------------|
| `synthetic_exact_lig147` | 101 bp | kept, unchanged |
| `synthetic_short` | 8 bp | dropped |
| `synthetic_long_lig147` | 105 bp | dropped |

The filter script's own defaults (target 154, window 145–163) drop all three
records. Those lengths are outside that window.

The 101 bp record is `A`×40 + `C`×40 + `G`×21, so a default 40-mer extract is
forty `A`s. No adapter FASTA and no lane FASTQs are included; those stay with
the experiment. `spacek40` is not included either.
