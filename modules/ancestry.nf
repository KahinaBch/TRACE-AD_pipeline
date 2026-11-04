// 1. Rename contigs (chrN -> N), extract the chromosome, keep only markers
//    shared with the reference panel: bin/pre_phasing_flare.sh (unchanged from sample_preparation/).
process PREPARE {
    label 'bcftools'

    input:
    val sanity_passed
    tuple path(vcf), path(tbi)
    tuple path(panel), path(panel_tbi)
    path contigs
    val chrom
    val prefix

    output:
    tuple path("${prefix}_query_com.vcf.gz"), path("${prefix}_ref_com.vcf.gz")

    script:
    """
    pre_phasing_flare.sh -sample ${vcf} -chr ${chrom} -ref ${panel} \\
      -contigs ${contigs} -out ${prefix} -threads ${task.cpus}
    """
}

// 2. Phase the cohort against the reference subset with Beagle 5.2
//    (same java command as the former sample_preparation/phasing_flare.sh).
process PHASE {
    label 'java'

    input:
    tuple path(query), path(ref)
    path genetic_map
    path beagle_jar
    val prefix

    output:
    tuple path("${prefix}_phased.vcf.gz"), path(ref, includeInputs: true), emit: vcf
    path "${prefix}_phased.log", emit: log

    script:
    """
    java -Xmx${params.beagle_mem_gb}g -jar ${beagle_jar} \\
         ref=${ref} \\
         gt=${query} \\
         map=${genetic_map} \\
         out=${prefix}_phased \\
         nthreads=${task.cpus}
    """
}

// 3. Index the phased output and verify the pair is FLARE-ready (fully phased,
//    no missing calls, same markers): bin/prepare_for_flare.sh (unchanged from sample_preparation/).
process CHECK_PHASED {
    label 'bcftools'

    input:
    tuple path(phased), path(ref)
    val chrom

    output:
    tuple path(phased, includeInputs: true), path("${phased}.tbi"), path(ref, includeInputs: true)

    script:
    """
    bcftools index -ft ${phased}
    bcftools index -ft ${ref}
    prepare_for_flare.sh -gt ${phased} -ref ${ref} -chr ${chrom}
    """
}

// 4. Local + global ancestry with FLARE 0.6.0 (same commands as the former
//    local_ancestry_flare/run_flare_*_chr*.job + run_flare.sh). The reference panel map is sorted by
//    ancestry label so FLARE's ancestry codes come out in a deterministic order.
process FLARE {
    label 'java'

    input:
    tuple path(phased), path(phased_tbi), path(ref)
    path sample_map
    path genetic_map
    path flare_jar
    val prefix

    output:
    path "${prefix}_flare.anc.vcf.gz",    emit: anc
    path "${prefix}_flare.global.anc.gz", emit: global
    path "${prefix}_flare.model",         emit: model
    path "${prefix}_flare.log",           emit: log

    script:
    """
    sort -k2,2 -k1,1 ${sample_map} > ${prefix}_refpanel.txt
    java -Xmx${params.flare_mem_gb}g -jar ${flare_jar} \\
         ref=${ref} \\
         ref-panel=${prefix}_refpanel.txt \\
         gt=${phased} \\
         map=${genetic_map} \\
         out=${prefix}_flare \\
         nthreads=${task.cpus} \\
         probs=false
    """
}
