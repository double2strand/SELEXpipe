# Tests

`test_sample_renaming.py` checks Sanger-ID renaming: a dry run leaves the tree alone, a real run builds the assay / TF / cycle paths from a few rows of `workflow/00_sample_renaming/examples/example_map.tsv`, and a flat HT-SELEX plus NCAP-SELEX clash exits 2 with nothing moved.

`test_safe_extract.py` checks the SEQ extract, the length filter, and the motif driver's width guard.

From the repository root, after `python3 -m pip install -r envs/requirements.txt`:

```bash
python3 -m unittest discover -s tests -v
```

`test_safe_extract.py` covers:

- an exact-length read yields one SEQ line of the requested width
- a short read is dropped, and the default output has no random padding
- the length filter applies the window and the `right`, `left`, and `center` trim modes
- `--seed` repeats padding
- a 154 bp filtered read still yields one 40-mer
- `run_motif.sh` rejects a SEQ line of the wrong width, and an empty SEQ file, before it runs `spacek40`

Both test modules load scripts under `workflow/`. `ROOT` is the parent of this directory (`Path(__file__).resolve().parents[1]`), the repository root.
