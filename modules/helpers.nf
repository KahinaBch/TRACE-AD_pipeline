// Helper: plink bed/bim/fam -> bgzipped VCF (same plink2 export as the
// A4/NATA sample-preparation jobs). Indexed afterwards by INDEX_EXPORT.
process PLINK_TO_VCF {
    label 'plink2'

    input:
    tuple path(bed), path(bim), path(fam)
    val out_prefix

    output:
    path "${out_prefix}.vcf.gz"

    script:
    def chr = params.chrom ? "--chr ${params.chrom}" : ''
    """
    plink2 --bfile ${bed.baseName} ${chr} \\
           --export vcf-4.2 bgz --out ${out_prefix} \\
           --threads ${task.cpus} --memory ${(task.memory.toMega() * 0.9) as long}
    """
}

// Helper: correct strand flips / REF-ALT swaps on the protein's chromosome.
process FIX_STRANDS {
    label 'python'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    tuple path(panel), path(panel_tbi)
    val chrom
    val out_name

    output:
    tuple path(out_name), path("${out_name}.tbi")

    script:
    """
    fix_strands.py --vcf ${vcf} --panel ${panel} --chrom ${chrom} --out ${out_name}
    """
}
