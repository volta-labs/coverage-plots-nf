#!/usr/bin/env python3
"""
build_gc_cache.py — Compute GC content for each interval in a BED file.

Outputs a CSV with columns: chrom, start, end, gc
where gc = fraction of G+C bases in that interval (NaN if region has >50% Ns).

Usage (local, with samtools on PATH and reference accessible):
    python3 build_gc_cache.py \
        --bed SG_MSKACCESSFLEX_v1_target_hg38.bed \
        --ref /path/to/hg38.fa \
        --out SG_MSKACCESSFLEX_v1_target_hg38_gc.csv

Usage (S3 reference via samtools with htslib S3 support):
    export AWS_DEFAULT_REGION=us-east-1
    python3 build_gc_cache.py \
        --bed SG_MSKACCESSFLEX_v1_target_hg38.bed \
        --ref s3://volta-labs-sequencing-data/reference-genome-index/bwa-mem2/hg38.fasta \
        --out SG_MSKACCESSFLEX_v1_target_hg38_gc.csv

Requirements:
    samtools >= 1.9 (with htslib; S3 access requires htslib compiled with libcurl)
    pandas
"""

import argparse
import subprocess
import sys
import pandas as pd


def compute_gc(ref: str, bed_path: str) -> pd.DataFrame:
    """
    Use samtools faidx to fetch each interval's sequence, then compute GC fraction.
    Handles both local FASTA and s3:// paths (requires htslib S3 support).
    """
    # Read BED
    records = []
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            records.append((chrom, start, end))

    print(f"Processing {len(records)} intervals...", flush=True)

    rows = []
    for i, (chrom, start, end) in enumerate(records):
        if i % 100 == 0:
            print(f"  {i}/{len(records)}", end='\r', flush=True)

        # samtools faidx uses 1-based closed coordinates
        region = f"{chrom}:{start + 1}-{end}"
        try:
            result = subprocess.run(
                ["samtools", "faidx", ref, region],
                capture_output=True, text=True, check=True
            )
            # FASTA output: first line is header, rest is sequence
            seq_lines = result.stdout.strip().split('\n')[1:]
            seq = ''.join(seq_lines).upper()

            if not seq:
                rows.append({'chrom': chrom, 'start': start, 'end': end, 'gc': float('nan')})
                continue

            n_count = seq.count('N')
            if n_count / len(seq) > 0.5:
                rows.append({'chrom': chrom, 'start': start, 'end': end, 'gc': float('nan')})
                continue

            gc = (seq.count('G') + seq.count('C')) / (len(seq) - n_count)
            rows.append({'chrom': chrom, 'start': start, 'end': end, 'gc': round(gc, 4)})

        except subprocess.CalledProcessError as e:
            print(f"\nWARNING: samtools faidx failed for {region}: {e.stderr.strip()}")
            rows.append({'chrom': chrom, 'start': start, 'end': end, 'gc': float('nan')})

    print(f"\nDone. {sum(1 for r in rows if pd.notna(r['gc']))} intervals with valid GC.")
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="Compute per-interval GC content for coverage-plots-nf gc_cache.")
    parser.add_argument('--bed', required=True, help='Input BED file (hg38, chr-prefixed)')
    parser.add_argument('--ref', required=True, help='Reference FASTA (local path or s3:// URL; must be indexed)')
    parser.add_argument('--out', required=True, help='Output CSV path')
    args = parser.parse_args()

    # Verify samtools available
    try:
        subprocess.run(['samtools', '--version'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        sys.exit("ERROR: samtools not found on PATH. Install via: conda install -c bioconda samtools")

    df = compute_gc(args.ref, args.bed)

    df.to_csv(args.out, index=False)
    print(f"Saved: {args.out} ({len(df)} rows)")

    # Summary stats
    valid = df['gc'].dropna()
    print(f"GC content range: {valid.min():.3f} – {valid.max():.3f}")
    print(f"GC content mean:  {valid.mean():.3f}")
    print(f"NaN intervals:    {df['gc'].isna().sum()}")


if __name__ == '__main__':
    main()
