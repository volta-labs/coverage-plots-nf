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
        samtools markdup -r -@ ${task.cpus} - dedup.bam

    bedtools coverage \\
        -a normalised_design.bed \\
        -b dedup.bam \\
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

// ── Process: coverage heterogeneity plot (across all samples in run) ──────────
process PLOT_HETEROGENEITY {
    conda 'conda-forge::python=3.11 conda-forge::numpy conda-forge::pandas conda-forge::matplotlib-base'

    publishDir "${params.outdir}/${params.run_name}", mode: 'copy'

    input:
    path(beds)         // collected coverage beds, named ${sample}.coverage.bed
    path(sample_sheet) // TSV with sample<TAB>group, built by collectFile in workflow
    path het_script

    output:
    path "*.png"

    script:
    """
    # Reconstruct qc_{sample}/ directory structure from bed filenames
    for bed in *.coverage.bed; do
        sample=\${bed%.coverage.bed}
        mkdir -p "qc_\${sample}"
        ln -sf "\${PWD}/\${bed}" "qc_\${sample}/hcs_coverage_raw.bed"
    done

    python3 ${het_script} \\
        --dir          . \\
        --run-name     ${params.run_name} \\
        --sample-sheet ${sample_sheet}
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

    PLOT_QC(BEDTOOLS_COVERAGE.out, gc_cache, qc_script)

    // Build sample sheet TSV from coverage output — avoids Groovy list flattening issues
    BEDTOOLS_COVERAGE.out
        .map { sample, group, bed -> "${sample}\t${group}" }
        .collectFile(name: 'sample_sheet.tsv', newLine: true, seed: 'sample\tgroup')
        .set { ch_sample_sheet }

    // Collect all coverage beds into one process
    BEDTOOLS_COVERAGE.out
        .map { sample, group, bed -> bed }
        .collect()
        .set { ch_all_beds }

    PLOT_HETEROGENEITY(ch_all_beds, ch_sample_sheet, het_script)
}
