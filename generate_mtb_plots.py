#!/usr/bin/env python3
"""
generate_mtb_plots.py
Post-analysis plots for the MTB WGS pipeline (run_mtb_pipeline.sh output).

Produces:
  1. genome_coverage.png     - depth across the H37Rv genome (from mosdepth per-base bed)
  2. depth_distribution.png  - histogram of per-base depth
  3. variant_summary.png     - SNP/Indel counts, PASS vs filtered (from filtered VCFs via bcftools stats)
  4. resistance_profile.png  - per-drug resistance status (from TB-Profiler JSON)
  5. fastp_qc.png            - read quality before/after trimming (from fastp JSON)

Usage:
    conda activate mtb-ngs
    python generate_mtb_plots.py -s SRR786503 -r ./results -o ./results/09_plots

Requires: matplotlib, pandas (pip/conda install if missing:
    conda install -n mtb-ngs -c conda-forge matplotlib pandas -y)
Uses `bcftools stats` (already installed via bioconda in mtb-ngs env) as a subprocess call.
"""
import argparse
import json
import gzip
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


# ------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Generate post-analysis plots for MTB pipeline output")
    p.add_argument("-s", "--sample", required=True, help="Sample name, e.g. SRR786503")
    p.add_argument("-r", "--results_dir", default="./results", help="Pipeline results/ directory")
    p.add_argument("-o", "--out_dir", default=None, help="Where to save plots (default: results/09_plots)")
    return p.parse_args()


# ------------------------------------------------------------------
def plot_genome_coverage(sample, results_dir, out_dir):
    """Uses mosdepth's per-base bed.gz (results/04_dedup/{sample}.per-base.bed.gz)."""
    bed_path = results_dir / "04_dedup" / f"{sample}.per-base.bed.gz"
    if not bed_path.exists():
        print(f"[skip] {bed_path} not found — skipping genome coverage plot")
        return

    positions, depths = [], []
    with gzip.open(bed_path, "rt") as f:
        for line in f:
            chrom, start, end, depth = line.strip().split("\t")
            start, end, depth = int(start), int(end), float(depth)
            # downsample: take midpoint of each interval, this keeps plotting fast
            positions.append((start + end) / 2)
            depths.append(depth)

    df = pd.DataFrame({"pos": positions, "depth": depths})
    # bin into ~2000 windows across genome for a smooth, fast-to-render line
    n_bins = 2000
    df["bin"] = pd.cut(df["pos"], bins=n_bins)
    binned = df.groupby("bin", observed=True).agg(pos=("pos", "mean"), depth=("depth", "mean")).dropna()

    mean_depth = df["depth"].mean()

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.fill_between(binned["pos"] / 1e6, binned["depth"], color="#3B6E8F", alpha=0.4)
    ax.plot(binned["pos"] / 1e6, binned["depth"], color="#1F4E5F", linewidth=0.8)
    ax.axhline(mean_depth, color="crimson", linestyle="--", linewidth=1, label=f"Mean depth = {mean_depth:.1f}x")
    ax.set_xlabel("Genome position (Mb)")
    ax.set_ylabel("Read depth (x)")
    ax.set_title(f"{sample} — Sequencing depth across H37Rv genome")
    ax.legend(loc="upper right")
    ax.set_xlim(0, binned["pos"].max() / 1e6)
    fig.tight_layout()
    out_path = out_dir / "genome_coverage.png"
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"[saved] {out_path}")


# ------------------------------------------------------------------
def plot_depth_distribution(sample, results_dir, out_dir):
    bed_path = results_dir / "04_dedup" / f"{sample}.per-base.bed.gz"
    if not bed_path.exists():
        print(f"[skip] {bed_path} not found — skipping depth distribution plot")
        return

    depths = []
    with gzip.open(bed_path, "rt") as f:
        for line in f:
            _, start, end, depth = line.strip().split("\t")
            depths.extend([float(depth)] * (int(end) - int(start)))

    depths_series = pd.Series(depths)
    # cap x-axis at 99th percentile for readability (avoid extreme outlier compression)
    cap = depths_series.quantile(0.99)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(depths_series[depths_series <= cap], bins=80, color="#3B6E8F", edgecolor="white", linewidth=0.3)
    ax.axvline(depths_series.mean(), color="crimson", linestyle="--", linewidth=1.2,
               label=f"Mean = {depths_series.mean():.1f}x")
    ax.axvline(depths_series.median(), color="orange", linestyle="--", linewidth=1.2,
               label=f"Median = {depths_series.median():.1f}x")
    ax.set_xlabel("Depth (x)")
    ax.set_ylabel("Number of bases")
    ax.set_title(f"{sample} — Depth distribution across genome")
    ax.legend()
    fig.tight_layout()
    out_path = out_dir / "depth_distribution.png"
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"[saved] {out_path}")


# ------------------------------------------------------------------
def plot_variant_summary(sample, results_dir, out_dir):
    """Runs bcftools stats on SNP and indel filtered VCFs to count PASS/filtered variants."""
    filt_dir = results_dir / "06_filtered"
    snp_vcf = filt_dir / f"{sample}.snp.filtered.vcf.gz"
    indel_vcf = filt_dir / f"{sample}.indel.filtered.vcf.gz"

    counts = {"SNP": {"PASS": 0, "Filtered": 0}, "Indel": {"PASS": 0, "Filtered": 0}}

    for label, vcf in [("SNP", snp_vcf), ("Indel", indel_vcf)]:
        if not vcf.exists():
            print(f"[skip] {vcf} not found — skipping {label} counts")
            continue
        try:
            out = subprocess.run(
                ["bcftools", "view", "-H", str(vcf)],
                capture_output=True, text=True, check=True
            ).stdout
        except FileNotFoundError:
            print("[error] bcftools not found on PATH — skipping variant_summary plot")
            return
        for line in out.strip().split("\n"):
            if not line:
                continue
            fields = line.split("\t")
            filter_col = fields[6]
            if filter_col == "PASS":
                counts[label]["PASS"] += 1
            else:
                counts[label]["Filtered"] += 1

    df = pd.DataFrame(counts).T  # rows: SNP, Indel; cols: PASS, Filtered
    if df.sum().sum() == 0:
        print("[skip] no variant counts found — skipping variant_summary plot")
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    df.plot(kind="bar", stacked=True, ax=ax, color=["#2E7D32", "#C62828"], edgecolor="white")
    ax.set_ylabel("Number of variants")
    ax.set_title(f"{sample} — Variant calls (PASS vs filtered)")
    ax.set_xticklabels(df.index, rotation=0)
    ax.legend(title="")
    for i, (idx, row) in enumerate(df.iterrows()):
        total = row.sum()
        ax.text(i, total + max(total * 0.02, 1), f"n={int(total)}", ha="center", fontsize=9)
    fig.tight_layout()
    out_path = out_dir / "variant_summary.png"
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"[saved] {out_path}")


# ------------------------------------------------------------------
def plot_resistance_profile(sample, results_dir, out_dir):
    json_path = results_dir / "08_resistance_lineage" / "results" / f"{sample}.results.json"
    if not json_path.exists():
        print(f"[skip] {json_path} not found — skipping resistance profile plot")
        return

    with open(json_path) as f:
        data = json.load(f)

    dr_variants = data.get("dr_variants", [])
    all_drugs = ["rifampicin", "isoniazid", "ethambutol", "pyrazinamide", "streptomycin",
                 "fluoroquinolones", "amikacin", "capreomycin", "kanamycin", "cycloserine",
                 "ethionamide", "clofazimine", "para-aminosalicylic_acid", "delamanid",
                 "bedaquiline", "linezolid"]

    resistant_drugs = set()
    for v in dr_variants:
        for d in v.get("drugs", []):
            drug_name = d.get("drug", "").lower().replace(" ", "_")
            resistant_drugs.add(drug_name)

    status = []
    for drug in all_drugs:
        label = drug.replace("_", " ").title()
        status.append((label, "Resistant" if drug in resistant_drugs else "Susceptible"))

    df = pd.DataFrame(status, columns=["Drug", "Status"])
    colors = df["Status"].map({"Resistant": "#C62828", "Susceptible": "#2E7D32"})

    fig, ax = plt.subplots(figsize=(8, 7))
    y_pos = range(len(df))
    ax.barh(y_pos, [1] * len(df), color=colors, edgecolor="white")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df["Drug"])
    ax.invert_yaxis()
    ax.set_xticks([])
    strain = data.get("sublin", data.get("main_lin", "N/A"))
    dr_type = data.get("drtype", "N/A")
    ax.set_title(f"{sample} — Drug resistance profile\nLineage: {strain} | Classification: {dr_type}")

    from matplotlib.patches import Patch
    legend_elems = [Patch(facecolor="#C62828", label="Resistant"),
                    Patch(facecolor="#2E7D32", label="Susceptible")]
    ax.legend(handles=legend_elems, loc="lower right")
    fig.tight_layout()
    out_path = out_dir / "resistance_profile.png"
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"[saved] {out_path}")


# ------------------------------------------------------------------
def plot_fastp_qc(sample, results_dir, out_dir):
    json_path = results_dir / "01_qc" / f"{sample}_fastp.json"
    if not json_path.exists():
        print(f"[skip] {json_path} not found — skipping fastp QC plot")
        return

    with open(json_path) as f:
        data = json.load(f)

    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})

    metrics = ["total_reads", "q20_rate", "q30_rate", "gc_content"]
    labels = ["Total reads (M)", "Q20 rate", "Q30 rate", "GC content"]

    before_vals = [
        before.get("total_reads", 0) / 1e6,
        before.get("q20_rate", 0),
        before.get("q30_rate", 0),
        before.get("gc_content", 0),
    ]
    after_vals = [
        after.get("total_reads", 0) / 1e6,
        after.get("q20_rate", 0),
        after.get("q30_rate", 0),
        after.get("gc_content", 0),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(14, 4))
    for i, (ax, label) in enumerate(zip(axes, labels)):
        ax.bar(["Before", "After"], [before_vals[i], after_vals[i]],
               color=["#90A4AE", "#3B6E8F"], edgecolor="white")
        ax.set_title(label, fontsize=10)
        for j, v in enumerate([before_vals[i], after_vals[i]]):
            ax.text(j, v, f"{v:.2f}" if v < 10 else f"{v:.1f}", ha="center", va="bottom", fontsize=8)
    fig.suptitle(f"{sample} — Read QC before vs after trimming (fastp)")
    fig.tight_layout()
    out_path = out_dir / "fastp_qc.png"
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"[saved] {out_path}")


# ------------------------------------------------------------------
def main():
    args = parse_args()
    results_dir = Path(args.results_dir).resolve()
    out_dir = Path(args.out_dir).resolve() if args.out_dir else results_dir / "09_plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Sample:      {args.sample}")
    print(f"Results dir: {results_dir}")
    print(f"Output dir:  {out_dir}\n")

    plot_genome_coverage(args.sample, results_dir, out_dir)
    plot_depth_distribution(args.sample, results_dir, out_dir)
    plot_variant_summary(args.sample, results_dir, out_dir)
    plot_resistance_profile(args.sample, results_dir, out_dir)
    plot_fastp_qc(args.sample, results_dir, out_dir)

    print(f"\nDone. Plots saved in: {out_dir}")


if __name__ == "__main__":
    main()
