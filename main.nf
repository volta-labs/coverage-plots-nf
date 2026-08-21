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

    output:
    tuple val(sample), val(group), path("${sample}.coverage.bed")

    script:
    """
    bedtools coverage \\
        -a ${bed} \\
        -b ${bam} \\
        -mean \\
        > ${sample}.coverage.bed
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

    output:
    path "*.png"

    script:
    """
    mkdir -p out
    plot_hcs_qc.py \\
        --targets ${coverage_bed} \\
        --gc      ${gc_cache} \\
        --outdir  out \\
        --sample  ${sample}
    mv out/*.png .
    """
}

// ── Process: coverage heterogeneity plot (across all samples in run) ──────────
process PLOT_HETEROGENEITY {
    conda 'conda-forge::python=3.11 conda-forge::numpy conda-forge::pandas conda-forge::matplotlib-base'

    publishDir "${params.outdir}/${params.run_name}", mode: 'copy'

    input:
    val(meta_list)   // collected list of [sample, group] pairs
    path(beds)       // collected coverage beds, named ${sample}.coverage.bed

    output:
    path "*.png"

    script:
    def setup_dirs = meta_list.collect { sample, group ->
        "mkdir -p qc_${sample} && ln -sf \"\${PWD}/${sample}.coverage.bed\" qc_${sample}/hcs_coverage_raw.bed"
    }.join("\n")

    def tsv_lines = meta_list.collect { sample, group -> "${sample}\t${group}" }.join("\\n")

    """
    ${setup_dirs}

    printf "sample\\tgroup\\n${tsv_lines}\\n" > sample_sheet.tsv

    plot_coverage_heterogeneity.py \\
        --dir          . \\
        --run-name     ${params.run_name} \\
        --sample-sheet sample_sheet.tsv
    """
}

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

    // Auto-discover BAM files from bam_dir; extract sample name from filename
    // Strips _aligned_sorted and _S## suffixes (e.g. SGULP15_1_S1_aligned_sorted.bam -> SGULP15_1)
    // Auto-discover BAMs; if groups_tsv provided, only process samples listed in it (allowlist)
    Channel
        .fromPath("${params.bam_dir}/*.bam")
        .map { bam ->
            def name = bam.name
                .replaceAll(/_aligned_sorted\.bam$/, '')
                .replaceAll(/_S\d+$/, '')
            def bai = bam.sibling("${bam.name}.bai")
            def group = group_map.containsKey(name) ? group_map[name] : params.run_name
            tuple(name, group, bam, bai)
        }
        .filter { name, group, bam, bai ->
            // If groups_tsv supplied, skip samples not listed in it
            params.groups_tsv ? group_map.containsKey(name) : true
        }
        .set { ch_input }

    bed      = file(params.bed)
    gc_cache = file(params.gc_cache)

    BEDTOOLS_COVERAGE(ch_input, bed)

    PLOT_QC(BEDTOOLS_COVERAGE.out, gc_cache)

    BEDTOOLS_COVERAGE.out
        .multiMap { sample, group, bed ->
            meta: [sample, group]
            bed:  bed
        }
        .set { cov_ch }

    cov_ch.meta.collect().set { ch_meta }
    cov_ch.bed.collect().set  { ch_beds }

    PLOT_HETEROGENEITY(ch_meta, ch_beds)
}
