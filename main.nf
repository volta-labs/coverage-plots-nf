#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ── Parameters ────────────────────────────────────────────────────────────────
params.input      = null   // required: S3 path to samplesheet CSV (sample,bam,bai)
params.groups_tsv = null   // optional: S3 path to groups TSV (sample<TAB>group) — same file used for sequencing pipeline
params.bed        = 's3://volta-labs-sequencing-data/reference-genome-index/hybrid-capture-target-files/hcs_v2_inferred_final_v1.bed'
params.gc_cache   = 's3://volta-labs-sequencing-data/reference-genome-index/hybrid-capture-target-files/hcs_v2_inferred_target_gc_v1.csv'
params.run_name   = 'HCS_QC'
params.outdir     = null   // required: S3 path for outputs, e.g. s3://volta-labs-sequencing-data/qc_results

if (!params.input)  { error "Please provide --input <samplesheet.csv>" }
if (!params.outdir) { error "Please provide --outdir <s3://bucket/path>" }

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
    // Build qc_{sample}/ directory structure expected by the discovery script
    def setup_dirs = meta_list.collect { sample, group ->
        "mkdir -p qc_${sample} && ln -sf \"\${PWD}/${sample}.coverage.bed\" qc_${sample}/hcs_coverage_raw.bed"
    }.join("\n")

    // Build sample sheet TSV (sample<TAB>group)
    def tsv_lines = meta_list.collect { sample, group -> "${sample}\t${group}" }.join("\\n")

    """
    # Reconstruct directory structure expected by plot_coverage_heterogeneity.py
    ${setup_dirs}

    # Write sample sheet
    printf "sample\\tgroup\\n${tsv_lines}\\n" > sample_sheet.tsv

    plot_coverage_heterogeneity.py \\
        --dir          . \\
        --run-name     ${params.run_name} \\
        --sample-sheet sample_sheet.tsv
    """
}

// ── Workflow ──────────────────────────────────────────────────────────────────
workflow {

    // Load groups TSV into a map if provided: { sample -> group }
    def group_map = [:]
    if (params.groups_tsv) {
        file(params.groups_tsv).readLines().drop(1).each { line ->
            def cols = line.split('\t')
            if (cols.size() >= 2) group_map[cols[0].trim()] = cols[1].trim()
        }
    }

    // Parse samplesheet — group comes from groups_tsv if supplied, else samplesheet group col, else run_name
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def group = group_map.containsKey(row.sample)
                ? group_map[row.sample]
                : (row.containsKey('group') && row.group ? row.group : params.run_name)
            tuple(row.sample, group, file(row.bam), file(row.bai))
        }
        .set { ch_input }

    bed      = file(params.bed)
    gc_cache = file(params.gc_cache)

    // Step 1: bedtools coverage per sample
    BEDTOOLS_COVERAGE(ch_input, bed)

    // Step 2: per-sample plots
    PLOT_QC(BEDTOOLS_COVERAGE.out, gc_cache)

    // Step 3: run-level heterogeneity plot
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
