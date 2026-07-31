# MTB WGS Pipeline

Paired-end Illumina WGS pipeline for *Mycobacterium tuberculosis*: QC → trim → align → dedup → variant call (ploidy=1) → hard-filter → annotate → lineage/drug-resistance profiling.

## Pipeline steps

| # | Step | Tool |
|---|------|------|
| 1 | Raw read QC | FastQC |
| 2 | Adapter/quality trimming | fastp |
| 3 | Alignment to H37Rv | BWA-MEM → sorted BAM |
| 4 | Mark duplicates + coverage stats | GATK MarkDuplicates, mosdepth |
| 5 | Variant calling (haploid) | GATK HaplotypeCaller (`-ploidy 1`) → GenotypeGVCFs |
| 6 | Hard-filtering | GATK VariantFiltration (SNP/indel separately, no VQSR — MTB has no truth-set VCF) |
| 7 | Functional annotation | SnpEff (optional, if installed) |
| 8 | Lineage + drug resistance | TB-Profiler |

## Tuned for your laptop (simranthakur-X556UR)

Intel i5-6200U (2C/4T @2.30GHz), 11GB RAM (~6.5GB available), 221GB free disk — MTB's ~4.4Mb genome is ~700x smaller than human WGS, so this laptop handles it fine, just slower than a server:

| Step | Rough time (1 sample, 30-60x coverage) |
|------|------------------------------------------|
| BWA-MEM alignment | ~10-20 min |
| HaplotypeCaller | ~5-10 min |
| Full pipeline, 1 sample | ~30-45 min |

Defaults already set in `config.sh`: `THREADS=3` (leaves 1 thread free), `GATK_MEM=4g` (JVM heap capped so GATK doesn't starve the OS — your available RAM is only ~6.5GB).

**For batch runs of many samples**, don't run them all on the laptop sequentially overnight — offload to your Arjun server (dual Xeon Gold 6248R) instead, same as your RH50/60/70 production runs. This pipeline works unchanged there — just bump `THREADS` and `GATK_MEM` in `config.sh` (e.g. `THREADS=16 GATK_MEM=16g`).

## Requirements

```bash
conda create -n mtb-ngs -c bioconda -c conda-forge \
    fastqc fastp bwa samtools gatk4 mosdepth tb-profiler snpeff -y
conda activate mtb-ngs
```

## Reference genome

Download H37Rv (RefSeq: `GCF_000195955.2`) and place as `reference/H37Rv.fasta`.
The pipeline auto-indexes it (BWA index, `.fai`, sequence dictionary) on first run.

```bash
datasets download genome accession GCF_000195955.2 --include genome
# unzip and move the .fna to reference/H37Rv.fasta
```

## Directory layout expected

```
raw_reads/
  sample1_R1.fastq.gz
  sample1_R2.fastq.gz
  sample2_R1.fastq.gz
  sample2_R2.fastq.gz
reference/
  H37Rv.fasta
```

## Downloading reads (ENA)

`wget` kept aborting/retrying on this network (EBI FTP drops connections repeatedly). **aria2c** worked reliably instead — multi-connection + auto-resume.

```bash
sudo apt install aria2 -y

# Current test sample: SRR786503
aria2c -x 8 -s 8 -c https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR786/SRR786503/SRR786503_1.fastq.gz
aria2c -x 8 -s 8 -c https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR786/SRR786503/SRR786503_2.fastq.gz

mkdir -p raw_reads
mv SRR786503_1.fastq.gz raw_reads/SRR786503_R1.fastq.gz
mv SRR786503_2.fastq.gz raw_reads/SRR786503_R2.fastq.gz
```

(Same commands are kept as comments at the top of `run_mtb_pipeline.sh` for reference — copy/uncomment for new samples.)

## Usage

Single sample (current downloaded sample):
```bash
./run_mtb_pipeline.sh -s SRR786503 -1 raw_reads/SRR786503_R1.fastq.gz -2 raw_reads/SRR786503_R2.fastq.gz
```

Batch mode — all `*_R1*/*_R2*` pairs in `raw_reads/`:
```bash
./run_mtb_pipeline.sh -a
```

Override defaults via env vars before running (or edit `config.sh`):
```bash
THREADS=16 RAW_DIR=/data/mtb_reads OUT_DIR=/data/mtb_results ./run_mtb_pipeline.sh -a
```

## Output structure

```
results/
  01_qc/                  FastQC + fastp reports
  02_trimmed/             trimmed FASTQs
  03_aligned/             sorted BAMs
  04_dedup/                dedup BAMs + duplicate metrics + mosdepth coverage
  05_variants/             raw GVCF/VCF
  06_filtered/             hard-filtered VCFs (SNP+indel merged), *.PASS.vcf.gz = final callset
  07_annotation/           SnpEff-annotated VCF + summary (if snpEff installed)
  08_resistance_lineage/   TB-Profiler per-sample lineage + drug-resistance report (txt/csv/json)
  logs/                    per-sample logs + pipeline.log
```

## Notes / rationale

- **Ploidy = 1**: MTB is haploid, set explicitly in HaplotypeCaller — this is the single most important non-default GATK flag for bacterial variant calling.
- **Hard filters, not VQSR/BQSR**: there's no curated truth-set VCF for MTB, so GATK's standard best-practice hard-filter thresholds (QD, FS, MQ, SOR, DP) are used instead, same approach as your dissertation pipeline.
- **TB-Profiler** replaces a custom resistance-lookup step — it calls lineage (via SNP barcode) and predicts resistance to first/second-line drugs directly from FASTQs, cross-referenced against curated resistance mutation databases (rpoB, katG, gyrA, etc.).
- Pipeline is idempotent per-sample and safe to re-run — reference indexing is skipped if already done.

## Extending

- Swap in `Snakemake`/`Nextflow` wrapper if you want DAG-based parallel scheduling across many samples (useful for your NGS SaaS idea).
- Add `qualimap` or `multiqc` in step 1/4 for a consolidated cross-sample QC report.
- Add a `PE_SAMPLE_SHEET.csv` driver if sample naming isn't uniform `_R1/_R2`.
