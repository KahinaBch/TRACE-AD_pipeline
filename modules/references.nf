// Index the phased reference panel and keep it in ref_dir for later runs.
// The panel itself was downloaded by Nextflow if it wasn't already in ref_dir.
process INDEX_PANEL {
    label 'bcftools'
    publishDir params.ref_dir, mode: 'copy', overwrite: false

    input:
    path panel

    output:
    tuple path(panel, includeInputs: true), path("${panel}.tbi")

    script:
    """
    bcftools index -ft --threads ${task.cpus} ${panel}
    """
}
