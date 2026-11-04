// Index a bgzipped VCF (the cohort input, or a helper's output).
process INDEX_VCF {
    label 'bcftools'

    input:
    path vcf

    output:
    tuple path(vcf, includeInputs: true), path("${vcf}.tbi")

    script:
    """
    bcftools index -ft --threads ${task.cpus} ${vcf} \\
      || { echo "ERROR: ${vcf} must be a bgzip-compressed VCF (use bgzip, not gzip)" >&2; exit 1; }
    """
}

// Pre-flight checks; never fails itself, the workflow reads STATUS from the report.
process SANITY_CHECK {
    label 'python'

    input:
    tuple path(vcf), path(tbi)
    tuple path(panel), path(panel_tbi)
    path targets
    val chrom

    output:
    path 'sanity_report.txt'

    script:
    """
    sanity_check.py --vcf ${vcf} --panel ${panel} --targets ${targets} \\
      --protein ${params.protein} --chrom ${chrom} ${params.ignore_flips ? '--ignore-flips' : ''} \\
      --out sanity_report.txt
    """
}
