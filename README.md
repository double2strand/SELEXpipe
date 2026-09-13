# SELEXpipe

**SELEXpipe** is a computational pipeline for collating, managing, and updating methods used to analyse Systematic Evolution of Ligands by Exponential Enrichment (**SELEX**) sequencing data. It provides a structured, extensible framework that spans the full analysis lifecycle: raw-read ingestion and quality control, sequence filtering and reformatting, motif discovery, motif position matching against enriched libraries, and mutual information (**MI**) analyses of sequence–binding relationships.

This repository is intended as a living methods hub. As preprocessing strategies, motif models, matching criteria, and information-theoretic analyses evolve, SELEXpipe will be updated so that workflows remain reproducible, comparable across experiments, and easy to extend.

---

## Scientific context

SELEX (and related protocols such as HT-SELEX) iteratively enriches nucleic acid libraries for sequences that bind a target of interest (typically a protein). High-throughput sequencing of successive rounds yields large oligonucleotide datasets in which binding-competent motifs are progressively over-represented relative to the naïve library.

Downstream analysis therefore requires:

1. **Reliable preprocessing** of multiplexed or adapter-containing reads into clean, analysis-ready sequence sets.
2. **Principled filtering** to remove low-quality, off-target, or artifactual sequences that would bias enrichment statistics.
3. **Consistent reformatting** so that heterogeneous inputs (FASTQ, FASTA, count tables, round-wise libraries) can feed a common analysis interface.
4. **Motif discovery** to recover position weight matrices (PWMs), consensus patterns, or related models that describe enriched binding sites.
5. **Motif position matching** to locate putative sites within individual sequences and quantify occupancy, position preferences, and enrichment across rounds.
6. **Mutual information (MI) analyses** to measure statistical dependencies between sequence features (e.g. base identity at paired positions, motif presence vs. enrichment rank) beyond what simple frequency counting captures.

SELEXpipe organises these stages into a coherent pipeline so that methods can be versioned, swapped, and compared without ad hoc scripting for each dataset.

---

## Pipeline overview

At a high level, the intended workflow is:

```text
Raw SELEX reads / libraries
        │
        ▼
┌───────────────────┐
│  Pre-processing   │  fastp adapter trim + PE merge (per lane),
│                   │  lane concat, FASTQ→FASTA, exact-length filter,
│                   │  one k-mer extract per molecule (default 40)
└─────────┬─────────┘
          │
          ▼
┌───────────────────┐
│  Filtering        │  enrichment-oriented filters, duplicate handling,
│                   │  contaminant removal, round-wise retention rules
└─────────┬─────────┘
          │
          ▼
┌───────────────────┐
│  Re-formatting    │  FASTA/FASTQ ↔ count tables, alphabet checks,
│                   │  standardised sequence IDs and metadata
└─────────┬─────────┘
          │
          ▼
┌───────────────────┐
│  Motif discovery  │  totalautomodelIC + spacek40 on SEQ files
│                   │  (background vs signal); seed refine; logos
└─────────┬─────────┘
          │
          ▼
┌───────────────────┐
│  Motif matching   │  scan sequences for motif instances, score
│                   │  thresholds, position and strand annotation
└─────────┬─────────┘
          │
          ▼
┌───────────────────┐
│  MI analyses      │  pairwise / higher-order mutual information,
│                   │  feature–enrichment associations, summaries
└───────────────────┘
```

Each stage is designed to be modular: inputs and outputs should be explicit, intermediate artefacts should be inspectable, and alternative methods for a given stage should be pluggable as the repository grows.

---

## Functional scope (current and planned)

### 1. Pre-processing

**Status: implemented** (see [`preprocessing/`](./preprocessing/)). Default path is exact-length + one extract per read; it does not inject synthetic sequence.

Pre-processing converts raw paired-end SELEX FASTQ reads (typically two lanes per sample) into analysis-ready SEQ files. The current workflow:

1. **`fastp`** — quality filtering (`-q 10`), adapter trimming (user-supplied adapter FASTA, e.g. TruSeq3 PE), overlap-based merging of R1/R2 (`-m`), polyG/complexity filters, and per-lane HTML/JSON QC reports.
2. **Lane concatenation** — `cat` the two lane-merged `.fq.gz` files into a single sample-level FASTQ.gz.
3. **Decompress** — `gzip -d -f -c` writes an uncompressed FASTQ **copy**; the `.fq.gz` is kept (do not use in-place `gzip -d`).
4. **FASTQ → FASTA** — `seqtk seq -a`.
5. **Length filter (`filter_by_length.py`)** — keep reads whose length equals `--target_length`, or the length implied by `--ligand` (see ligand table below; no hard-coded default). Default window is exact (`min = max = target`). Shorts are dropped. Random ATCG padding is **not** the default; it requires `--allow-pad` and `--seed`.
6. **One-window extract** — `trim_to_kmer.py` writes **one** SEQ line per retained molecule (default **40 bp**, `--start` or `--center`). Short windows are dropped. This is not a tiled 40-mer library and does not invent trailing k-mers.

#### Ligand design expected lengths

These are merged construct lengths. They are correct; do not change them.

| Ligand design | Expected merged length |
|---------------|------------------------|
| lig147        | 101 bp                 |
| lig200        | 154 bp                 |
| ligN40        | 187 bp                 |
| ligN70        | 247 bp                 |

**Flanks differ:** ligN40 and ligN70 use different constant flanks. The 60 bp merged-length delta (187 vs 247) is expected — it is not a 30 bp N-region typo. Do not "correct" 187 / 247.

Inputs use a generic sample prefix (`SAMPLE_L1_R1.fastq.gz`, …). Full usage is documented in [`preprocessing/README.md`](./preprocessing/README.md).

### 2. Filtering

Filtering removes sequences that would distort motif and MI estimates. Planned and evolving criteria include:

- Length constraints matching the randomised region design.
- Minimum quality or maximum ambiguous-base (`N`) thresholds.
- Optional deduplication or UMI-aware collapse.
- Exclusion of known contaminants, primer dimers, or off-design constructs.
- Round-aware retention policies (e.g. minimum count thresholds in later rounds).

### 3. Re-formatting

Re-formatting standardises heterogeneous SELEX artefacts into pipeline-native representations, for example:

- Sequence–count tables keyed by round.
- FASTA exports for motif discovery tools.
- Metadata tables linking sample, target, round, and experimental conditions.
- Validation of alphabet (DNA/RNA), fixed-flank consistency, and variable-region extraction.

### 4. Motif discovery

**Status: implemented** (see [`motif/`](./motif/)). The driver refuses SEQ files whose line width ≠ `--seq-length` (default 40).

Motif discovery consumes SEQ files from preprocessing (`trim_to_kmer.py`, one extract per molecule) and runs **`totalautomodelIC`** (with **`spacek40`**) to:

- Discover local-maximum k-mer seeds in signal vs. background libraries.
- Build and trim PFMs / logos for each seed.
- Filter by complexity and logo similarity.
- Refine seeds to an information-content window (10–20 bits).
- Emit accepted / unrefined / rejected models as SVG logos and PFM artefacts.

Generic naming replaces experiment-specific IDs, e.g.:

```bash
./totalautomodelIC -40N BACKGROUND_101.seq SIGNAL_r4_101.seq \
  1 6 6 0.25 - 20 50 0.1 SIGNAL_r4v1
```

See [`motif/README.md`](./motif/README.md) for the full argument table, driver usage, and the requirement to supply the `spacek40` binary locally.

### 5. Motif position matching

Once motifs are available, sequences are scanned to annotate putative binding sites:

- Score each window against one or more motif models.
- Apply significance or score thresholds.
- Record match position, strand, and score per sequence.
- Aggregate match statistics across rounds (occupancy, positional bias, enrichment of match-containing sequences).

### 6. Mutual information (MI) analyses

MI analyses quantify dependencies that simple motif logos may miss:

- Pairwise MI between positions within the variable region (covariation / structural or specificity coupling).
- MI between motif features and enrichment-related variables (e.g. round, count rank, binder vs. non-binder labels when available).
- Summaries and visualisations that support interpretation of specificity and higher-order sequence constraints.

Exact estimators, regularisation choices, and finite-sample corrections will be documented as implementations land in the repository.

---

## Repository status

This repository is in an early stage. The documentation above describes the **intended architecture and scientific scope**. Implementation modules, command-line interfaces, configuration schemas, and example datasets will be added incrementally. The README will be updated in lockstep with those changes so that it remains an accurate technical description of available functionality.

| Area                    | Status                                      |
|-------------------------|---------------------------------------------|
| Pre-processing          | Implemented (`preprocessing/`; default path drops shorts, no pad) |
| Filtering               | Planned                                     |
| Re-formatting           | Planned                                     |
| Motif discovery         | Implemented (`motif/`; refuses width-mismatched SEQ) |
| Motif position matching | Planned                                     |
| MI analyses             | Planned                                     |
| End-to-end CLI / config | Partial (preprocess + motif drivers)        |
| Example datasets / tests| Partial (`tests/` covers extract / drop / no-pad) |

---

## Design principles

- **Methods collation**: Prefer a curated set of well-defined analysis steps over one-off scripts, so alternative algorithms can be compared under the same I/O contracts.
- **Reproducibility**: Pin inputs, parameters, and motif model versions; retain intermediate artefacts needed to regenerate figures and tables.
- **Modularity**: Keep preprocessing, filtering, motif, and MI stages separable so individual methods can be updated without rewriting the full workflow.
- **Extensibility**: New filters, motif engines, or MI estimators should plug into existing stage boundaries rather than forking the pipeline.
- **Transparency**: Document assumptions (library design, alphabet, round structure, scoring thresholds) alongside code.

---

## Getting started

### Preprocessing

Requires `fastp`, `seqtk`, `gzip`, `bash`, `python3`, and Biopython (`pip install -r requirements.txt`). Place lane FASTQs using the generic naming scheme, then run:

```bash
./preprocessing/run_preprocess.sh \
  --sample SAMPLE \
  --adapter-fasta ./Adapters_TruSeq3PE.fa \
  --ligand lig147 \
  --input-dir ./raw \
  --output-dir ./processed
```

See [`preprocessing/README.md`](./preprocessing/README.md) for file naming, `fastp` defaults, standalone helper usage, and the `--allow-pad` opt-in.

Helper contract (no FASTQ required):

```bash
python3 -m unittest discover -s tests -v
```

### Motif processing

Requires `bash`, `bc`, and the `spacek40` binary placed alongside `motif/totalautomodelIC`. SEQ inputs come from preprocessing. Example:

```bash
./motif/run_motif.sh \
  --background ./processed/BACKGROUND/BACKGROUND_101.seq \
  --signal ./processed/SIGNAL_r4/SIGNAL_r4_101.seq \
  --output-prefix SIGNAL_r4v1 \
  --work-dir ./motif_out/SIGNAL_r4
```

See [`motif/README.md`](./motif/README.md). Further stages (motif position matching, MI) will gain install/usage docs as they are implemented.

---

## Contributing / evolution of this document

SELEXpipe is expected to grow as SELEX analysis methods are collated and refined. When functionality is added or behaviour changes:

1. Update the relevant section under **Functional scope**.
2. Adjust the **Repository status** table.
3. Document new interfaces, parameters, and file formats in dedicated usage sections.
4. Keep scientific terminology precise (rounds, enrichment, PWM/PSSM, MI) so the README remains useful as a technical reference.

---

## License

License information will be added when a license is chosen for the project.

---

## Citation / acknowledgement

Citation guidance will be provided once the pipeline reaches a citable release. Until then, please refer to this repository (`https://github.com/tw599/SELEXpipe`) when using or adapting its methods.
