# coverage-plots-nf

Nextflow pipeline for post-capture QC plots on HCS v2 sequencing runs. Runs on Seqera Platform (AWS Batch + Wave + Fusion FS).

## What it does

For each sample BAM, the pipeline runs three plots:

1. **Chromosome coverage** — per-chromosome mean coverage across HCS v2 targets
2. **GC bias** — coverage as a function of target GC content
3. **Coverage heterogeneity** — fraction of targets with normalized coverage < 0.2 or > 5.0, plotted across all samples in the run, color-coded by group

---

## Inputs

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `bam_dir` | Yes | — | S3 folder containing `*.bam` and `*.bam.bai` files |
| `outdir` | Yes | — | S3 path for output plots |
| `run_name` | No | `HCS_QC` | Label used in plot titles and output folder name |
| `groups_tsv` | No | — | S3 path to a groups TSV (see below). If provided, acts as both a sample allowlist and group assignment |
| `bed` | No | HCS v2 BED on S3 | S3 path to target BED file |
| `gc_cache` | No | HCS v2 GC cache on S3 | S3 path to per-target GC content CSV |

### BAM file naming

BAMs must follow the naming convention used by the Volta Labs sequencing pipeline. The pipeline auto-strips the following suffixes to derive sample names:

```
SGULP12_1_aligned_sorted_S3.bam  →  SGULP12_1
HG002_Manual_131_aligned_sorted.bam  →  HG002_Manual_131
```

BAI files must be co-located in the same S3 folder with the `.bam.bai` extension (e.g. `SGULP12_1_aligned_sorted.bam.bai`).

### Groups TSV

This is the same TSV generated during sequencing run setup (the "Group ID" file). It must have two tab-separated columns with a header row:

```
sample	group
SGULP12_1	32-1
SGULP12_2	32-1
SGULP15_1	32-2
HG002_Manual_131	control
```

When `--groups_tsv` is provided:
- Only samples listed in the TSV are processed (acts as an allowlist — unlisted BAMs are skipped)
- Each sample is assigned its group ID for the heterogeneity plot
- Groups are sorted naturally in the plot (32-1, 32-2, ..., 33-1, 33-2, ...)

When omitted, all BAMs in `bam_dir` are processed and assigned to a single group named after `run_name`.

---

## Outputs

```
<outdir>/<run_name>/
├── coverage_heterogeneity_<run_name>.png   # run-level heterogeneity plot
└── per_sample/
    ├── <sample>_chromosome_coverage.png
    ├── <sample>_gc_bias.png
    └── ...
```

---

## Running on Seqera Platform

1. Add this repo as a pipeline in your Seqera workspace (Pipeline → Add Pipeline → GitHub URL).
2. Set the branch to `master`.
3. Select the `vl-nextflow-ssd_copy` compute environment.
4. Fill in the launch parameters:

```json
{
  "bam_dir":    "s3://volta-labs-sequencing-data/illumina/ILMN159/bam/",
  "groups_tsv": "s3://volta-labs-sequencing-data/illumina/ILMN159/ilmn159_groups.tsv",
  "run_name":   "ILMN159",
  "outdir":     "s3://volta-labs-sequencing-data/illumina/ILMN159/coverage_plots"
}
```

`bed` and `gc_cache` default to the HCS v2 reference files and do not need to be set unless you are using a different target panel.

---

## Running locally (for testing)

```bash
nextflow run . \
  -profile local \
  --bam_dir   /path/to/bam/folder \
  --groups_tsv /path/to/groups.tsv \
  --run_name  TEST \
  --outdir    /path/to/output \
  --bed       /path/to/targets.bed \
  --gc_cache  /path/to/gc_cache.csv
```

The `local` profile disables Wave and Fusion; conda environments are built locally.

---

## Reference files (HCS v2 defaults)

| File | S3 path |
|------|---------|
| Target BED | `s3://volta-labs-sequencing-data/reference-genome-index/hybrid-capture-target-files/hcs_v2_inferred_final_v1.bed` |
| GC content cache | `s3://volta-labs-sequencing-data/reference-genome-index/hybrid-capture-target-files/hcs_v2_inferred_target_gc_v1.csv` |

If you use a different BED file, you must generate a matching GC cache CSV using `fetch_all_gc.py` (once per BED) and supply both via `--bed` and `--gc_cache`.

---

## Resource requirements

| Process | CPUs | Memory | Time |
|---------|------|--------|------|
| BEDTOOLS_COVERAGE | 4 | 16 GB | 2 h |
| PLOT_QC | 2 | 8 GB | 30 min |
| PLOT_HETEROGENEITY | 2 | 8 GB | 30 min |

---

## Repository structure

```
coverage-plots-nf/
├── main.nf                   # Pipeline definition
├── nextflow.config           # AWS Batch, Wave, Fusion, resource config
├── bin/
│   ├── plot_hcs_qc.py        # Chromosome coverage + GC bias plots
│   └── plot_coverage_heterogeneity.py  # Run-level heterogeneity plot
└── assets/
    └── samplesheet_template.csv
```
