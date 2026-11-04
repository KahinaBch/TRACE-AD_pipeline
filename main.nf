#!/usr/bin/env nextflow
/*
 * TRACE-AD local ancestry pipeline (FLARE)
 *
 *   nextflow run . --vcf cohort.vcf.gz --cohort ADNI --protein TREM2 -profile snellius
 *     -> results_ADNI_TREM2/  global ancestry, local ancestry at the protein's
 *                             variants, FLARE model, run report
 * Helpers:
 *   --mode plink_to_vcf --bfile path/to/prefix [--chrom 6] [--out name]
 *   --mode fix_strands  --vcf cohort.vcf.gz --protein TREM2 [--out name.vcf.gz]
 */

include { INDEX_VCF; INDEX_VCF as INDEX_EXPORT; SANITY_CHECK } from './modules/checks'
include { INDEX_PANEL }                                        from './modules/references'
include { PREPARE; PHASE; CHECK_PHASED; FLARE }                from './modules/ancestry'
include { EXTRACT_LOCAL; REPORT }                              from './modules/results'
include { PLINK_TO_VCF; FIX_STRANDS }                          from './modules/helpers'

def helpMessage() {
    def proteins = file(params.targets).splitCsv(header: true)
        .groupBy { row -> row.protein }
        .collect { protein, rows -> "              ${protein.padRight(8)} chr${rows[0].chrom}: ${rows*.label.join(', ')}" }
        .join('\n')
    return """
    TRACE-AD local ancestry pipeline (FLARE)
    Global ancestry, and local ancestry at a protein's target variants, for one cohort.
    Quick start: README.md   Details: https://github.com/KahinaBch/TRACE-AD_pipeline/wiki

    Check the installation first (synthetic data, no input files needed):
      nextflow run <pipeline> -profile test,docker        (or test,apptainer / test,singularity)

    Main run:
      nextflow run <pipeline> --vcf cohort.vcf.gz --cohort NAME --protein PROTEIN -profile <profiles>

      --vcf       cohort genotypes: bgzipped VCF (.vcf.gz), GRCh38, with a GT field.
                  Unphased is fine; one chromosome or the whole genome; chromosomes
                  named '19' or 'chr19'. Plink files? Convert them first (see helpers).
      --cohort    name used for the output folder and files, e.g. ADNI
      --protein   one of (assets/targets.csv):
${proteins}

    Output: <outdir>/results_<cohort>_<protein>/
      <cohort>_<chr>_flare.global.anc.gz      global ancestry per participant (FLARE)
      <cohort>_<protein>.local_ancestry.tsv   one row per participant: allele and local
                                              ancestry of each haplotype at each target variant
      <cohort>_<chr>_flare.model              FLARE model parameters
      sanity_report.txt, run_report.txt       input checks, versions, how each variant was resolved

    Optional:
      --outdir DIR              where results_<cohort>_<protein>/ is written (default: .)
      --ref_dir DIR             cache for the phased reference panel (default: ./trace_ad_references)
      --ignore_flips            continue even if the sanity check finds strand flips
      --threads N               CPUs for phasing and FLARE (default: ${params.threads})
      --beagle_mem_gb N / --flare_mem_gb N   Java memory (default: ${params.beagle_mem_gb} / ${params.flare_mem_gb})
      --flare_jar / --beagle_jar PATH        use a local jar instead of downloading it

    Helpers:
      --mode plink_to_vcf --bfile path/prefix [--chrom N] [--out NAME]
                        plink .bed/.bim/.fam -> bgzipped, indexed VCF in --outdir
      --mode fix_strands  --vcf cohort.vcf.gz --protein PROTEIN [--out NAME.vcf.gz]
                        corrects strand flips / REF-ALT swaps found by the sanity check

    Profiles (-profile, comma-separated):
      where steps run    slurm | pbspro | pbs | sge | lsf | condor | local
      container engine   apptainer | singularity | docker
      ready-made         snellius (SURF Snellius) | test (installation test)
    """.stripIndent()
}

def require(names) {
    names.each { name -> if (!params[name]) error "Missing required parameter --${name} (see --help)" }
}

// protein -> its chromosome, checked against assets/targets.csv before anything runs
def proteinChrom() {
    require(['protein'])
    def rows = file(params.targets).splitCsv(header: true)
    def mine = rows.findAll { row -> row.protein.equalsIgnoreCase(params.protein) }
    if (!mine)
        error "Unknown protein '${params.protein}'. Available in ${params.targets}: ${rows*.protein.unique().join(', ')}"
    def chroms = mine*.chrom.unique()
    if (chroms.size() != 1)
        error "${params.protein}: target variants span several chromosomes (${chroms}); one per run is supported"
    if (!(chroms[0] in params.supported_chroms.toString().tokenize(',')))
        error "${params.protein} is on chr${chroms[0]}, which has no published phased reference panel (supported: ${params.supported_chroms})"
    return chroms[0]
}

// phased reference panel: reused from ref_dir, otherwise downloaded once and kept there
def panelSource(chrom) {
    def local = file("${params.ref_dir}/ref_${chrom}_phased.vcf.gz")
    if (local.exists()) return local
    def url = params.panel_url.replace('{chr}', chrom)
    log.info "chr${chrom} reference panel not in ${params.ref_dir}: downloading it once from ${url}\n" +
             "  (if this fails -- no internet on this node, or no access -- put ref_${chrom}_phased.vcf.gz in ${params.ref_dir} yourself)"
    return file(url)
}

// stop the run with a clear message (and what to do next) if the sanity check failed
def checkSanity(report) {
    def text = report.text
    if (text.contains('STATUS: PASS')) return true
    def msg = "Sanity check failed -- full report: ${params.outdir}/results_${params.cohort}_${params.protein}/sanity_report.txt\n  " +
              text.readLines().findAll { line -> line.startsWith('FAIL:') }.join('\n  ')
    if (!text.contains('BUILD: GRCh38')) {
        msg += "\n\nThe pipeline needs GRCh38 input. Lift the VCF over to GRCh38 first (best option if available: " +
               "re-impute on a GRCh38 panel such as TOPMed), then re-run -- see https://github.com/KahinaBch/TRACE-AD_pipeline/wiki/Input-on-GRCh37"
    } else if (text =~ /RECOVERABLE_FLIPS: [1-9]/) {
        def pipeline = workflow.repository ? "${workflow.repository} -r ${workflow.revision}" : workflow.projectDir
        msg += "\n\nFix them with:\n  nextflow run ${pipeline} --mode fix_strands --vcf ${params.vcf} --protein ${params.protein} -profile ${workflow.profile}\n" +
               "then re-run with --vcf <the .fixed.vcf.gz it writes> (or add --ignore_flips to continue anyway)."
    }
    error msg
}

workflow RUN {
    require(['vcf', 'cohort'])
    def chrom   = proteinChrom()
    def prefix  = "${params.cohort}_${chrom}"
    def targets = file(params.targets)
    def gmap    = file(params.genetic_map.replace('{chr}', chrom), checkIfExists: true)
    def source  = panelSource(chrom)

    panel  = INDEX_PANEL(source)
    input  = INDEX_VCF(file(params.vcf, checkIfExists: true))
    sanity = SANITY_CHECK(input, panel, targets, chrom)

    PREPARE(sanity.map { report -> checkSanity(report) }, input, panel, file(params.contigs), chrom, prefix)
    PHASE(PREPARE.out, gmap, file(params.beagle_jar), prefix)
    CHECK_PHASED(PHASE.out.vcf, chrom)
    FLARE(CHECK_PHASED.out, file(params.sample_map), gmap, file(params.flare_jar), prefix)
    EXTRACT_LOCAL(CHECK_PHASED.out, FLARE.out.anc, targets, chrom)
    REPORT(sanity, EXTRACT_LOCAL.out.notes, PHASE.out.log, FLARE.out.log, chrom, source.name)
}

workflow PLINK2VCF {
    require(['bfile'])
    def out = (params.out ?: file(params.bfile).name).replaceAll(/\.vcf\.gz$/, '')
    def plink = ['bed', 'bim', 'fam'].collect { ext -> file("${params.bfile}.${ext}", checkIfExists: true) }
    PLINK_TO_VCF(plink, out)
    INDEX_EXPORT(PLINK_TO_VCF.out)
}

workflow FIXSTRANDS {
    require(['vcf'])
    def chrom = proteinChrom()
    def vcf   = file(params.vcf, checkIfExists: true)
    def out   = params.out ?: "${vcf.simpleName}.chr${chrom}.fixed.vcf.gz"
    FIX_STRANDS(INDEX_VCF(vcf), INDEX_PANEL(panelSource(chrom)), chrom, out)
}

workflow {
    if (params.help)                        log.info(helpMessage())
    else if (params.mode == 'run')          RUN()
    else if (params.mode == 'plink_to_vcf') PLINK2VCF()
    else if (params.mode == 'fix_strands')  FIXSTRANDS()
    else error "Unknown --mode '${params.mode}' (expected run, plink_to_vcf or fix_strands)"
}
