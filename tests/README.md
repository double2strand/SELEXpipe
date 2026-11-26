# Tests

`test_safe_extract.py` checks the SEQ extract, the length filter, and the motif driver's width guard.

From the repository root, after `python3 -m pip install -r envs/requirements.txt`:

```bash
python3 -m unittest discover -s tests -v
```

The file covers:

- an exact-length read yields one SEQ line of the requested width
- a short read is dropped, and the default output has no random padding
- the length filter applies the window and the `right`, `left`, and `center` trim modes
- `--seed` repeats padding
- a 154 bp filtered read still yields one 40-mer
- `run_motif.sh` rejects a SEQ line of the wrong width, and an empty SEQ file, before it runs `spacek40`

## Paths inside the test file

`test_safe_extract.py` was not edited. It still loads and invokes the old locations, and one assertion looks for the string `scripts/02_motif` in the driver log (that log prints the absolute `--bin-dir`). The scripts now live under `workflow/`. Until those references are updated, `unittest` fails while importing the module.

| Line | Reference that still uses the old path |
|------|----------------------------------------|
| 36 | `scripts/01_preprocess/trim_to_kmer.py` |
| 37 | `scripts/01_preprocess/filter_by_length.py` |
| 235 | `scripts/01_preprocess/filter_by_length.py` |
| 309 | `scripts/01_preprocess/run_preprocess.sh` |
| 324 | `scripts/02_motif/run_motif.sh` |
| 371 | `scripts/02_motif/run_motif.sh` |
| 390 | stdout must contain `scripts/02_motif` |

`ROOT` is the parent of this directory (`Path(__file__).resolve().parents[1]`), which is still the repository root while this file stays in `tests/`.
