process PLOT_HETEROGENEITY {
    conda 'conda-forge::python=3.11 conda-forge::numpy conda-forge::pandas conda-forge::matplotlib-base'

    publishDir { "${params.outdir}/${run_name}" }, mode: 'copy'

    input:
    path(beds)          // collected coverage beds (any suffix: .coverage.bed, .mapq20.coverage.bed, etc.)
    path(sample_sheet)  // TSV with sample<TAB>group, built by collectFile in workflow
    path het_script
    val  run_name       // label used for plot title, output filename, and publishDir subfolder

    output:
    path "*.png"

    script:
    """
    # Reconstruct qc_{sample}/ directory structure from bed filenames.
    # Strip .mapqNN suffix from sample names so groups_tsv lookup still works.
    for bed in *.coverage.bed; do
        sample=\${bed%.coverage.bed}
        sample=\$(echo "\${sample}" | sed 's/\\.mapq[0-9]*\$//')
        mkdir -p "qc_\${sample}"
        ln -sf "\${PWD}/\${bed}" "qc_\${sample}/hcs_coverage_raw.bed"
    done

    python3 ${het_script} \\
        --dir          . \\
        --run-name     ${run_name} \\
        --sample-sheet ${sample_sheet}
    """
}
