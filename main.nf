#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ── Parameters ────────────────────────────────────────────────────────────────
params.bam_dir    = null   // required: S3 folder containing *.bam and *.bam.bai files
params.groups_tsv = null   // optional: S3 path to groups TSV (sample<TAB>group)
params.bed        = 's3://volta-labs-sequencing-data/reference-genome-index/hybrid-capture-target-files/hcs_v2_inferred_final_v1.bed'
params.gc_cache   = 's3://volta-labs-sequencing-data/reference-genome-index/hybrid-capture-target-files/hcs_v2_inferred_target_gc_v1.csv'
params.run_name   = 'HCS_QC'
params.outdir     = null   // required: S3 path for outputs

if (!params.bam_dir) { error "Please provide --bam_dir <s3://bucket/path/to/bam/folder/>" }
if (!params.outdir)  { error "Please provide --outdir <s3://bucket/path>" }

// ── Process: per-target coverage ─────────────────────────────────────────────
process BEDTOOLS_COVERAGE {
    tag "$sample"

    conda 'bioconda::bedtools=2.31.1 bioconda::samtools=1.19'

    input:
    tuple val(sample), val(group), path(bam), path(bai)
    path bed

    publishDir "${params.outdir}/${params.run_name}/markdup_stats", mode: 'copy', pattern: "*.markdup.txt"

    output:
    tuple val(sample), val(group), path("${sample}.coverage.bed")
    path "${sample}.markdup.txt"
    path "${sample}.mapq20.coverage.bed"
    path "${sample}.mapq30.coverage.bed"

    script:
    """
    # Normalise BED chromosome names to match BAM (add 'chr' prefix if absent).
    # Handles both chr-prefixed BEDs and non-prefixed BEDs (e.g. SOPHiA HCS v2.0).
    # Idempotent: lines already starting with 'chr' pass through unchanged.
    awk 'BEGIN{OFS="\\t"} /^[^#]/ && \$1 !~ /^chr/ {\$1="chr"\$1} {print}' ${bed} > normalised_design.bed

    # Remove duplicate reads before coverage calculation to match SOPHiA's molecule-based counting.
    # Input BAMs are pre-deduplication (_aligned_sorted.bam).
    # markdup requires mate scores: collate (name-sort) -> fixmate -m -> coord-sort -> markdup -r.
    samtools collate -O -u -@ ${task.cpus} ${bam} tmp_collate | \
        samtools fixmate -m -u - - | \
        samtools sort -u -@ ${task.cpus} - | \
        samtools markdup -r -s -f ${sample}.markdup.txt -@ ${task.cpus} - dedup.bam

    # No MAPQ filter
    bedtools coverage \\
        -a normalised_design.bed \\
        -b dedup.bam \\
        -mean \\
        > ${sample}.coverage.bed

    # MAPQ >= 20 filtered coverage
    samtools view -bq 20 -@ ${task.cpus} dedup.bam | \\
        bedtools coverage -a normalised_design.bed -b - -mean \\
        > ${sample}.mapq20.coverage.bed

    # MAPQ >= 30 filtered coverage
    samtools view -bq 30 -@ ${task.cpus} dedup.bam | \\
        bedtools coverage -a normalised_design.bed -b - -mean \\
        > ${sample}.mapq30.coverage.bed
    """
}

// ── Process: chromosome coverage + GC bias plots (per sample) ────────────────
process PLOT_QC {
    tag "$sample"

    conda 'conda-forge::python=3.11 conda-forge::numpy conda-forge::pandas conda-forge::matplotlib-base conda-forge::scipy'

    publishDir "${params.outdir}/${params.run_name}/per_sample", mode: 'copy'

    input:
    tuple val(sample), val(group), path(coverage_bed)
    path gc_cache
    path plot_script

    output:
    path "*.png"

    script:
    """
    mkdir -p out
    python3 ${plot_script} \\
        --targets ${coverage_bed} \\
        --gc      ${gc_cache} \\
        --outdir  out \\
        --sample  ${sample}
    mv out/*.png .
    """
}

// ── Import PLOT_HETEROGENEITY three times with aliases (DSL2 multi-call pattern) ─
include { PLOT_HETEROGENEITY as PLOT_HETEROGENEITY_NOFILTER } from './modules/plot_heterogeneity'
include { PLOT_HETEROGENEITY as PLOT_HETEROGENEITY_MAPQ20   } from './modules/plot_heterogeneity'
include { PLOT_HETEROGENEITY as PLOT_HETEROGENEITY_MAPQ30   } from './modules/plot_heterogeneity'

// ── Workflow ──────────────────────────────────────────────────────────────────
workflow {

    // Load groups TSV into a map: { sample -> group }
    def group_map = [:]
    if (params.groups_tsv) {
        file(params.groups_tsv).readLines().drop(1).each { line ->
            def cols = line.split('\t')
            if (cols.size() >= 2) group_map[cols[0].trim()] = cols[1].trim()
        }
    }

    // Auto-discover BAMs; if groups_tsv provided, only process samples listed in it
    Channel
        .fromPath("${params.bam_dir}/*.bam")
        .map { bam ->
            def name = bam.name
                .replaceAll(/_aligned_sorted\.bam$/, '')
                .replaceAll(/_S\d+$/, '')
            def bai  = bam.parent.resolve("${bam.name}.bai")
            def group = group_map.containsKey(name) ? group_map[name] : params.run_name
            tuple(name, group, bam, bai)
        }
        .filter { name, group, bam, bai ->
            params.groups_tsv ? group_map.containsKey(name) : true
        }
        .set { ch_input }

    bed        = file(params.bed)
    gc_cache   = file(params.gc_cache)
    qc_script  = file("${projectDir}/bin/plot_hcs_qc.py")
    het_script = file("${projectDir}/bin/plot_coverage_heterogeneity.py")

    BEDTOOLS_COVERAGE(ch_input, bed)

    PLOT_QC(BEDTOOLS_COVERAGE.out[0], gc_cache, qc_script)

    // Build sample sheet TSV from coverage output — avoids Groovy list flattening issues
    BEDTOOLS_COVERAGE.out[0]
        .map { sample, group, bed -> "${sample}\t${group}" }
        .collectFile(name: 'sample_sheet.tsv', newLine: true, seed: 'sample\tgroup')
        .set { ch_sample_sheet }

    // Collect all coverage beds into one process (three sets: no filter, mapq20, mapq30)
    BEDTOOLS_COVERAGE.out[0]
        .map { sample, group, bed -> bed }
        .collect()
        .set { ch_all_beds }

    BEDTOOLS_COVERAGE.out[2]
        .collect()
        .set { ch_mapq20_beds }

    BEDTOOLS_COVERAGE.out[3]
        .collect()
        .set { ch_mapq30_beds }

    PLOT_HETEROGENEITY_NOFILTER(ch_all_beds,    ch_sample_sheet, het_script, params.run_name)
    PLOT_HETEROGENEITY_MAPQ20(ch_mapq20_beds,  ch_sample_sheet, het_script, "${params.run_name}_mapq20")
    PLOT_HETEROGENEITY_MAPQ30(ch_mapq30_beds,  ch_sample_sheet, het_script, "${params.run_name}_mapq30")
}
