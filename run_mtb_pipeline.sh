#!/usr/bin/env bash
# ============================================================
# run_mtb_pipeline.sh
# Mycobacterium tuberculosis WGS pipeline
# QC -> Trim -> Align (BWA-MEM) -> Dedup -> Call (GATK, ploidy=1)
#   -> Hard-filter -> Annotate (SnpEff) -> Lineage/Resistance (TB-Profiler)
#
# Usage:
#   ./run_mtb_pipeline.sh -s <sample_name> [-1 R1.fastq.gz -2 R2.fastq.gz]
#   ./run_mtb_pipeline.sh -a   (run all samples found in RAW_DIR, paired *_R1/_R2)
#
# Requires: fastp, bwa, samtools, gatk4, tb-profiler, mosdepth, snpEff (optional)
# ============================================================
#
# ------------------------------------------------------------
# SAMPLE DOWNLOAD (ENA) — run manually before the pipeline, NOT auto-run here.
# Uncomment and run these lines in your terminal to fetch reads.
# aria2c is used because wget kept aborting/retrying on this network (EBI FTP
# drops connections) — aria2c's multi-connection resume handled it reliably.
# ------------------------------------------------------------
#
# --- Install aria2c (one-time) ---
# sudo apt install aria2 -y
#
# --- SRR786503 (current test sample) ---
# aria2c -x 8 -s 8 -c https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR786/SRR786503/SRR786503_1.fastq.gz
# aria2c -x 8 -s 8 -c https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR786/SRR786503/SRR786503_2.fastq.gz
#
# --- Move + rename into RAW_DIR with the naming this script expects ---
# mkdir -p raw_reads
# mv SRR786503_1.fastq.gz raw_reads/SRR786503_R1.fastq.gz
# mv SRR786503_2.fastq.gz raw_reads/SRR786503_R2.fastq.gz
#
# --- Add more samples below the same way, e.g.: ---
# aria2c -x 8 -s 8 -c https://ftp.sra.ebi.ac.uk/vol1/fastq/SRRxxx/SRRxxxxxx/SRRxxxxxx_1.fastq.gz
# aria2c -x 8 -s 8 -c https://ftp.sra.ebi.ac.uk/vol1/fastq/SRRxxx/SRRxxxxxx/SRRxxxxxx_2.fastq.gz
# mv SRRxxxxxx_1.fastq.gz raw_reads/SRRxxxxxx_R1.fastq.gz
# mv SRRxxxxxx_2.fastq.gz raw_reads/SRRxxxxxx_R2.fastq.gz
# ------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

usage() {
    echo "Usage: $0 -s SAMPLE -1 R1.fastq.gz -2 R2.fastq.gz"
    echo "       $0 -a                 # process all pairs in ${RAW_DIR}"
    exit 1
}

SAMPLE=""; R1=""; R2=""; RUN_ALL=false
while getopts "s:1:2:ah" opt; do
    case $opt in
        s) SAMPLE="$OPTARG" ;;
        1) R1="$OPTARG" ;;
        2) R2="$OPTARG" ;;
        a) RUN_ALL=true ;;
        h|*) usage ;;
    esac
done

mkdir -p "$QC_DIR" "$TRIM_DIR" "$ALIGN_DIR" "$DEDUP_DIR" "$VCF_DIR" "$FILT_DIR" "$ANNOT_DIR" "$RESIST_DIR" "$LOG_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_DIR}/pipeline.log" >&2; }

# ------------------------------------------------------------
check_reference() {
    if [[ ! -f "$REF_FASTA" ]]; then
        log "ERROR: Reference FASTA not found at $REF_FASTA"
        log "Download H37Rv: datasets download genome accession GCF_000195955.2"
        exit 1
    fi
    if [[ ! -f "${REF_FASTA}.bwt" ]]; then
        log "Indexing reference (BWA)..."
        $BWA index "$REF_FASTA"
    fi
    if [[ ! -f "${REF_FASTA}.fai" ]]; then
        $SAMTOOLS faidx "$REF_FASTA"
    fi
    local dict="${REF_FASTA%.fasta}.dict"
    if [[ ! -f "$dict" ]]; then
        log "Creating sequence dictionary..."
        $GATK CreateSequenceDictionary -R "$REF_FASTA" -O "$dict"
    fi
}

# ------------------------------------------------------------
step_qc() {
    local sample=$1 r1=$2 r2=$3
    log "[$sample] Step 1: Raw QC (FastQC)"
    fastqc -t "$THREADS" -o "$QC_DIR" "$r1" "$r2" &>> "${LOG_DIR}/${sample}.log" || true
}

# ------------------------------------------------------------
step_trim() {
    local sample=$1 r1=$2 r2=$3
    local out1="${TRIM_DIR}/${sample}_R1.trim.fastq.gz"
    local out2="${TRIM_DIR}/${sample}_R2.trim.fastq.gz"
    log "[$sample] Step 2: Adapter/quality trimming (fastp)"
    $FASTP -i "$r1" -I "$r2" -o "$out1" -O "$out2" \
        --qualified_quality_phred "$MIN_QUAL" \
        --length_required "$MIN_LEN" \
        --thread "$THREADS" \
        --json "${QC_DIR}/${sample}_fastp.json" \
        --html "${QC_DIR}/${sample}_fastp.html" \
        &>> "${LOG_DIR}/${sample}.log"
    echo "$out1 $out2"
}

# ------------------------------------------------------------
step_align() {
    local sample=$1 t1=$2 t2=$3
    local bam="${ALIGN_DIR}/${sample}.sorted.bam"
    log "[$sample] Step 3: Alignment (BWA-MEM) + sort"
    $BWA mem -t "$THREADS" \
        -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}" \
        "$REF_FASTA" "$t1" "$t2" 2>> "${LOG_DIR}/${sample}.log" \
        | $SAMTOOLS sort -@ "$THREADS" -o "$bam" -
    $SAMTOOLS index "$bam"
    echo "$bam"
}

# ------------------------------------------------------------
step_dedup() {
    local sample=$1 bam=$2
    local dedup_bam="${DEDUP_DIR}/${sample}.dedup.bam"
    local metrics="${DEDUP_DIR}/${sample}.dup_metrics.txt"
    log "[$sample] Step 4: Mark duplicates (GATK)"
    $GATK --java-options "-Xmx${GATK_MEM}" MarkDuplicates \
        -I "$bam" -O "$dedup_bam" -M "$metrics" \
        &>> "${LOG_DIR}/${sample}.log"
    $SAMTOOLS index "$dedup_bam"

    log "[$sample] Step 4b: Coverage stats (mosdepth)"
    $MOSDEPTH -t "$THREADS" -x "${DEDUP_DIR}/${sample}" "$dedup_bam" &>> "${LOG_DIR}/${sample}.log" || true

    echo "$dedup_bam"
}

# ------------------------------------------------------------
step_call_variants() {
    local sample=$1 bam=$2
    local gvcf="${VCF_DIR}/${sample}.g.vcf.gz"
    local vcf="${VCF_DIR}/${sample}.vcf.gz"
    log "[$sample] Step 5: Variant calling (GATK HaplotypeCaller, ploidy=$PLOIDY)"
    $GATK --java-options "-Xmx${GATK_MEM}" HaplotypeCaller \
        -R "$REF_FASTA" -I "$bam" \
        -ploidy "$PLOIDY" \
        -ERC GVCF \
        -O "$gvcf" \
        &>> "${LOG_DIR}/${sample}.log"

    $GATK --java-options "-Xmx${GATK_MEM}" GenotypeGVCFs \
        -R "$REF_FASTA" -V "$gvcf" -O "$vcf" \
        &>> "${LOG_DIR}/${sample}.log"
    echo "$vcf"
}

# ------------------------------------------------------------
step_filter_variants() {
    local sample=$1 vcf=$2
    local snp_vcf="${FILT_DIR}/${sample}.snp.vcf.gz"
    local indel_vcf="${FILT_DIR}/${sample}.indel.vcf.gz"
    local final_vcf="${FILT_DIR}/${sample}.filtered.vcf.gz"

    log "[$sample] Step 6: Hard-filtering (no truth set exists for MTB -> hard filters, not VQSR)"

    $GATK SelectVariants -R "$REF_FASTA" -V "$vcf" --select-type-to-include SNP -O "$snp_vcf" &>> "${LOG_DIR}/${sample}.log"
    $GATK VariantFiltration -R "$REF_FASTA" -V "$snp_vcf" \
        --filter-expression "QD < ${QD_FILTER}" --filter-name "QD_fail" \
        --filter-expression "FS > ${FS_FILTER}" --filter-name "FS_fail" \
        --filter-expression "MQ < ${MQ_FILTER}" --filter-name "MQ_fail" \
        --filter-expression "SOR > ${SOR_FILTER}" --filter-name "SOR_fail" \
        --filter-expression "DP < ${MIN_DP}" --filter-name "LowDepth" \
        -O "${FILT_DIR}/${sample}.snp.filtered.vcf.gz" &>> "${LOG_DIR}/${sample}.log"

    $GATK SelectVariants -R "$REF_FASTA" -V "$vcf" --select-type-to-include INDEL -O "$indel_vcf" &>> "${LOG_DIR}/${sample}.log"
    $GATK VariantFiltration -R "$REF_FASTA" -V "$indel_vcf" \
        --filter-expression "QD < ${QD_FILTER}" --filter-name "QD_fail" \
        --filter-expression "FS > 200.0" --filter-name "FS_fail" \
        --filter-expression "DP < ${MIN_DP}" --filter-name "LowDepth" \
        -O "${FILT_DIR}/${sample}.indel.filtered.vcf.gz" &>> "${LOG_DIR}/${sample}.log"

    $GATK MergeVcfs \
        -I "${FILT_DIR}/${sample}.snp.filtered.vcf.gz" \
        -I "${FILT_DIR}/${sample}.indel.filtered.vcf.gz" \
        -O "$final_vcf" &>> "${LOG_DIR}/${sample}.log"

    # PASS-only VCF for downstream annotation
    $GATK SelectVariants -R "$REF_FASTA" -V "$final_vcf" --exclude-filtered \
        -O "${FILT_DIR}/${sample}.PASS.vcf.gz" &>> "${LOG_DIR}/${sample}.log"

    echo "${FILT_DIR}/${sample}.PASS.vcf.gz"
}

# ------------------------------------------------------------
step_annotate() {
    local sample=$1 vcf=$2
    log "[$sample] Step 7: Functional annotation (SnpEff, if configured)"
    if command -v snpEff &> /dev/null; then
        snpEff -v Mycobacterium_tuberculosis_h37rv "$vcf" \
            > "${ANNOT_DIR}/${sample}.snpeff.vcf" 2>> "${LOG_DIR}/${sample}.log"
        mv snpEff_summary.html "${ANNOT_DIR}/${sample}_snpEff_summary.html" 2>/dev/null || true
        mv snpEff_genes.txt "${ANNOT_DIR}/${sample}_snpEff_genes.txt" 2>/dev/null || true
    else
        log "[$sample] snpEff not found on PATH — skipping functional annotation step"
    fi
}

# ------------------------------------------------------------
step_resistance_lineage() {
    local sample=$1 t1=$2 t2=$3
    log "[$sample] Step 8: Lineage + drug-resistance profiling (TB-Profiler)"
    if command -v $TBPROFILER &> /dev/null; then
        local t1_abs t2_abs logfile
        t1_abs="$(realpath "$t1")"
        t2_abs="$(realpath "$t2")"
        logfile="$(realpath "$LOG_DIR")/${sample}.log"
        ( cd "$RESIST_DIR" && \
          $TBPROFILER profile -1 "$t1_abs" -2 "$t2_abs" -p "$sample" \
            -t "$THREADS" --txt --csv &>> "$logfile" )
    else
        log "[$sample] tb-profiler not found on PATH — skipping resistance/lineage step"
        log "Install: pip install tb-profiler  OR  conda install -c bioconda tb-profiler"
    fi
}

# ------------------------------------------------------------
run_sample() {
    local sample=$1 r1=$2 r2=$3
    log "===== Starting pipeline for sample: $sample ====="

    step_qc "$sample" "$r1" "$r2"
    read -r t1 t2 <<< "$(step_trim "$sample" "$r1" "$r2")"
    bam=$(step_align "$sample" "$t1" "$t2")
    dedup_bam=$(step_dedup "$sample" "$bam")
    raw_vcf=$(step_call_variants "$sample" "$dedup_bam")
    pass_vcf=$(step_filter_variants "$sample" "$raw_vcf")
    step_annotate "$sample" "$pass_vcf"
    step_resistance_lineage "$sample" "$t1" "$t2"

    log "===== Finished sample: $sample | Final PASS VCF: $pass_vcf ====="
}

# ------------------------------------------------------------
main() {
    check_reference

    if $RUN_ALL; then
        shopt -s nullglob
        for r1file in "${RAW_DIR}"/*_R1*.fastq.gz; do
            base=$(basename "$r1file")
            sample=${base%%_R1*}
            r2file="${RAW_DIR}/${sample}_R2${base#*_R1}"
            if [[ -f "$r2file" ]]; then
                run_sample "$sample" "$r1file" "$r2file"
            else
                log "WARNING: no matching R2 for $r1file — skipping"
            fi
        done
    else
        [[ -z "$SAMPLE" || -z "$R1" || -z "$R2" ]] && usage
        run_sample "$SAMPLE" "$R1" "$R2"
    fi

    log "Pipeline complete. Results in: $OUT_DIR"
}

main "$@"
