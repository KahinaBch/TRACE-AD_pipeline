// 5. Wide local-ancestry table at the protein's target variants.
process EXTRACT_LOCAL {
    label 'python'

    input:
    tuple path(phased), path(phased_tbi), path(ref)
    path anc
    path targets
    val chrom

    output:
    path "${params.cohort}_${params.protein}.local_ancestry.tsv", emit: table
    path 'targets_notes.tsv',                                      emit: notes

    script:
    """
    extract_local_ancestry.py --phased ${phased} --anc ${anc} --targets ${targets} \\
      --protein ${params.protein} --chrom ${chrom} --max-distance ${params.max_proxy_distance} \\
      --out ${params.cohort}_${params.protein}.local_ancestry.tsv --notes targets_notes.tsv
    """
}

// 6. run_report.txt: what was run, with which versions, and what came out.
process REPORT {
    label 'python'

    input:
    path sanity
    path notes
    path beagle_log
    path flare_log
    val chrom
    val panel_name

    output:
    path 'run_report.txt'

    script:
    """
    build_report.py --cohort ${params.cohort} --protein ${params.protein} --chrom ${chrom} \\
      --input '${params.vcf}' --panel '${panel_name}' --sanity ${sanity} --notes ${notes} \\
      --beagle-log ${beagle_log} --flare-log ${flare_log} \\
      --pipeline-version '${workflow.commitId ?: 'local copy'}' --out run_report.txt
    """
}
