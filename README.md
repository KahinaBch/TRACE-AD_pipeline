# TRACE-AD local ancestry pipeline

Gives each participant of a cohort their **global ancestry** and their **local
ancestry at a protein's risk variants** (APOE, TREM2), with FLARE and the gnomAD
HGDP + 1000 Genomes reference panel.

All parameters, cluster settings, adding a protein, GRCh37 input: **[the wiki](https://github.com/KahinaBch/TRACE-AD_pipeline/wiki)**. Options: `--help`.

## 1. Install (once)

You need Java 17+ and a container engine (Apptainer/Singularity or Docker; most
HPCs have Apptainer). Then:

```
curl -s https://get.nextflow.io | bash
```

## 2. Test that it works (no data needed)

```
./nextflow run KahinaBch/TRACE-AD_pipeline -r main -profile test,apptainer
```

Use `test,singularity` or `test,docker` depending on you container engine. 
It takes a few minutes on built-in synthetic data. **It works if** it finishes without error and
`pipeline_test/results_TEST_APOE/` contains 5 files (see results Section).

To find your container engine, check which of these commands works on the machine where
you run the pipeline, and use the first one that does:

| Command | Profile |
|---|---|
| `apptainer --version` | `apptainer` |
| `singularity --version` | `singularity` |
| `docker --version` | `docker` (laptops and servers; HPCs rarely allow it) |

None works? Ask your HPC support to make Apptainer available (often `module load apptainer`).

## 3. Check your input file

Your cohort must be **one VCF file** that is:

- [x] compressed with `bgzip` (ends in `.vcf.gz`)
- [x] on genome build **GRCh38** (hg38)
- [x] with genotypes (a `GT` field), one column per participant
- [x] containing the protein's chromosome: **chr19 for APOE, chr6 for TREM2**

Phased or not, imputed or not, one chromosome or the whole genome: all fine.
The pipeline checks all of this before starting and tells you what to fix.

## 4. Run

```
./nextflow run KahinaBch/TRACE-AD_pipeline -r main -profile slurm,apptainer \
    --vcf cohort.vcf.gz --cohort ADNI --protein TREM2
```

Replace the values after:
- `--vcf`: the path to your VCF file
- `--cohort`: a short name for your cohort, e.g. `ADNI` (your choice; it is only used to name the results)
- `--protein`: one of

  | `--protein` | Chromosome | Target variants |
  |---|---|---|
  | `APOE` | chr19 | rs429358, rs7412 |
  | `TREM2` | chr6 | R47H (rs75932628), R62H (rs143332484) |

  To add a protein, see the [wiki](https://github.com/KahinaBch/TRACE-AD_pipeline/wiki/Adding-a-protein).

`-profile` takes two names: **your cluster's scheduler**, then **your container engine**
(`apptainer`, `singularity` or `docker`). To find your scheduler, check which of these
commands exists on your cluster (e.g. `which sbatch`):

| Command available | Scheduler | Profile |
|---|---|---|
| `sbatch` | SLURM | `slurm` |
| `qsub` | PBS Pro, PBS/Torque or Grid Engine | `pbspro`, `pbs` or `sge` (ask your HPC support which one) |
| `bsub` | LSF | `lsf` |
| `condor_submit` | HTCondor | `condor` |
| none (laptop, single server) | | `local` |

On SURF Snellius, `-profile snellius` sets both. Queue names, accounts, and other
cluster-specific settings: see the [wiki](https://github.com/KahinaBch/TRACE-AD_pipeline/wiki/Running-the-pipeline).

## 5. Results

In `results_ADNI_TREM2/`:

| File | Content |
|---|---|
| `ADNI_6_flare.global.anc.gz` | global ancestry of each participant (proportion per region) |
| `ADNI_TREM2.local_ancestry.tsv` | for each participant and each variant: the allele and local ancestry of both haplotypes |
| `ADNI_6_flare.model` | FLARE model parameters |
| `sanity_report.txt`, `run_report.txt` | input checks, versions, notes on each variant |

## Additionally: If the pipeline stops

It prints the reason and what to do. The three common cases:

- **Wrong input format**: the file is not a bgzipped VCF, has no genotypes (`GT`), or lacks
  the protein's chromosome.
  - VCF compressed with plain `gzip`, or not compressed: recompress it with `bgzip`:
    ```
    gunzip -c cohort.vcf.gz | bgzip -c > cohort.bgz.vcf.gz     # was gzip-compressed
    bgzip -c cohort.vcf > cohort.vcf.gz                         # was not compressed
    ```
  - Plink files (`.bed/.bim/.fam`): convert them first (same `-profile` as in step 4);
    this writes `cohort.vcf.gz` in the current folder:
    ```
    ./nextflow run KahinaBch/TRACE-AD_pipeline -r main -profile slurm,apptainer \
        --mode plink_to_vcf --bfile path/to/cohort
    ```
- **Strand flips found**: run the fix command it prints (`--mode fix_strands ...`), then
  re-run step 4 with the fixed file it creates.
- **Input is GRCh37**: convert it to GRCh38 first ([wiki: Input on GRCh37](https://github.com/KahinaBch/TRACE-AD_pipeline/wiki/Input-on-GRCh37)).

An interrupted run continues where it stopped if you add `-resume` to the same command.
Example: 
```
./nextflow run KahinaBch/TRACE-AD_pipeline -r main -profile slurm,apptainer \
    --vcf cohort.vcf.gz --cohort ADNI --protein TREM2 -resume
```
