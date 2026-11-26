# Environment

Installs for a workstation. This directory has a conda environment and a pip requirements file. A Softpack module file for the Sanger FARM cluster is not in this repository.

| File | What it installs |
|------|------------------|
| `environment.yml` | Conda env `selexpipe`: Python ≥ 3.9, Biopython ≥ 1.81, fastp, seqtk, bc |
| `requirements.txt` | `biopython>=1.81`, enough for the length filter, the SEQ extract, and the tests |

[fastp](https://github.com/opengene/fastp) and seqtk were not version-pinned in the scripts, so `environment.yml` does not pin them either. `spacek40` is the binary at `workflow/02_motif_discovery/spacek40`. It is not a conda package. The comment in `environment.yml` names that path. `requirements.txt` is this directory. `selex_rename.py` needs only the Python standard library.

From the repository root:

```bash
conda env create -f envs/environment.yml
conda activate selexpipe
```

Python only:

```bash
python3 -m pip install -r envs/requirements.txt
```

`totalautomodelIC` also uses `bash`, `awk`, `sed`, and `grep`, which come from the system. `bc` is in the conda environment and is also expected on `PATH` for the motif recipe.
