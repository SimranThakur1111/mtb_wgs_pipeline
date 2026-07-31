#!/usr/bin/env bash
# ============================================================
# config.sh — MTB WGS Pipeline Configuration
# Edit paths below before running run_mtb_pipeline.sh
# ============================================================

# ---- Reference genome (H37Rv, NC_000962.3) ----
REF_DIR="${REF_DIR:-./reference}"
REF_FASTA="${REF_DIR}/H37Rv.fasta"          # download from NCBI: GCF_000195955.2
REF_GFF="${REF_DIR}/H37Rv.gff"              # for SnpEff / annotation (optional)

# ---- Directories ----
RAW_DIR="${RAW_DIR:-./raw_reads}"           # input FASTQs: {sample}_R1.fastq.gz, {sample}_R2.fastq.gz
                                              # current sample: SRR786503
                                              # (SRR786503_R1.fastq.gz / SRR786503_R2.fastq.gz)
OUT_DIR="${OUT_DIR:-./results}"
QC_DIR="${OUT_DIR}/01_qc"
TRIM_DIR="${OUT_DIR}/02_trimmed"
ALIGN_DIR="${OUT_DIR}/03_aligned"
DEDUP_DIR="${OUT_DIR}/04_dedup"
VCF_DIR="${OUT_DIR}/05_variants"
FILT_DIR="${OUT_DIR}/06_filtered"
ANNOT_DIR="${OUT_DIR}/07_annotation"
RESIST_DIR="${OUT_DIR}/08_resistance_lineage"
LOG_DIR="${OUT_DIR}/logs"

# ---- Run parameters ----
# Tuned for: Intel i5-6200U (2C/4T, 2.30GHz), 11GB RAM (~6.5GB available), 221GB free disk
THREADS="${THREADS:-3}"                      # leave 1 thread free for OS/background
GATK_MEM="${GATK_MEM:-4g}"                   # cap JVM heap so it doesn't starve the OS
PLOIDY=1                                     # MTB is haploid
MIN_QUAL=20                                  # fastp quality trimming threshold
MIN_LEN=50                                   # min read length post-trim

# NOTE: MTB genome (~4.4 Mb) is ~700x smaller than human WGS, so this laptop
# handles it fine, just slower than a server. Rough expectations:
#   BWA-MEM alignment:   ~10-20 min / sample (30-60x coverage, paired 150bp)
#   HaplotypeCaller:     ~5-10 min / sample
#   Full single sample:  ~30-45 min end-to-end
# For batch runs of many samples, consider offloading to your Arjun server instead.

# ---- GATK hard-filter thresholds (no truth VCF exists for MTB -> hard filtering, not VQSR/BQSR) ----
QD_FILTER=2.0
FS_FILTER=60.0
MQ_FILTER=40.0
SOR_FILTER=3.0
MIN_DP=10                                    # min depth to keep a variant

# ---- Tool paths (assumes on $PATH; override if needed) ----
FASTP=fastp
BWA=bwa
SAMTOOLS=samtools
GATK=gatk
TBPROFILER=tb-profiler
MOSDEPTH=mosdepth
