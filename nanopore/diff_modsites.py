#!/usr/bin/env python3
import argparse
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import numpy as np
import re
import os

topN=20

mod_code_to_name = {
    "m": "5-Methylcytosine",
    "h": "5-Hydroxymethylcytosine",
    "f": "5-Formylcytosine",
    "c": "5-Carboxylcytosine",
    "C": "Ambiguitycode;anyCmod",
    "g": "5-Hydroxymethyluracil",
    "e": "5-Formyluracil",
    "b": "5-Carboxyluracil",
    "T": "Ambiguitycode;anyTmod",
    "U": "Ambiguitycode;anyUmod",
    "a": "6-Methyladenine",
    "A": "Ambiguitycode;anyAmod",
    "o": "8-Oxoguanine",
    "G": "Ambiguitycode;anyGmod",
    "n": "Xanthosine",
    "N": "Ambiguitycode;anymod"
}

def load_pileup(path,depth_filter, mod_fraction):
    """
    读取 modkit pileup.bed
    假设列顺序为：
    chr start end base depth strand ... mod_fraction mod_count unmod_count ...
    """
    df = pd.read_csv(path, sep="\t", header=None)
    df = df.iloc[:, :15]  # 截取前15列，避免版本差异
    df.columns = ["chr","start","end","base","depth","strand",
                  "s2","s3","color","depth2","mod_fraction",
                  "mod_count","unmod_count","x1","x2"]
    #add curated
    ref = pd.read_table('../evaluation/curated.gene.list',header=None, names=['chr','start','end','gene','score','strand'])
    result = []
    for _, row in ref.iterrows():
        sub = df[
            (df["chr"] == row["chr"]) &
            (df["start"] >= row["start"]) &
            (df["end"] <= row["end"])
        ].copy()
        result.append(sub)
    df = pd.concat(result, ignore_index=True)
    #
    if depth_filter>0 or mod_fraction > 0:
        df = df.loc[(df['depth'] >= depth_filter) & (df['mod_fraction'] >= mod_fraction), ]
    df = df[["chr","start","end","base","strand",
             "mod_count","unmod_count","mod_fraction"]]
    return df

def merge_sites(ctrl, treat,summary_base_h):
    merged = pd.merge(ctrl, treat,
                      on=["chr","start","end","base","strand"],
                      how="inner",
                      suffixes=("_ctrl","_treat"))
    summary_base_h.write('%s sites in control\n' % ctrl.shape[0])
    summary_base_h.write('%s sites in treat\n' % treat.shape[0])
    summary_base_h.write('%s common sites between control and treat used for analysis\n' % merged.shape[0])
    return merged

def fisher_test(row):
    table = [[row["mod_count_ctrl"], row["unmod_count_ctrl"]],
             [row["mod_count_treat"], row["unmod_count_treat"]]]
    try:
        _, p = fisher_exact(table)
    except Exception:
        p = 1.0
    return p

def chrom_order(chrom):
    # 主染色体
    if chrom.startswith("chr"):
        core = chrom[3:]  # 去掉 'chr'
        if core.isdigit():
            return int(core)
        elif core == "X":
            return 23
        elif core == "Y":
            return 24
        elif core == "M":
            return 25
    # 非主染色体放后面
    return 1000 + hash(chrom) % 1000  # 保持唯一顺序

def manhattan_plot(df, outpng, base, fdr_thresh=0.05):
    # 为每个染色体分配累积位置
    # chroms = sorted(df["chr"].unique(), key=lambda x: int(re.sub(r"\D", "", x)))
    chroms = sorted(df["chr"].unique(), key=chrom_order)
    chrom_offsets = {}
    current_offset = 0
    xticks, xticklabels = [], []

    for chrom in chroms:
        chr_len = df.loc[df["chr"] == chrom, "end"].max()
        chrom_offsets[chrom] = current_offset
        xticks.append(current_offset + chr_len / 2)
        xticklabels.append(chrom)
        current_offset += chr_len

    df["pos_cum"] = df.apply(lambda r: r["start"] + chrom_offsets[r["chr"]], axis=1)

    # 画图
    plt.figure(figsize=(12, 5))
    colors = ["#4daf4a", "#377eb8"]
    for i, chrom in enumerate(chroms):
        subset = df[df["chr"] == chrom]
        plt.scatter(subset["pos_cum"], -np.log10(subset["fdr"]+1e-300),
                    c=colors[i % 2], s=8, alpha=0.6)

    # 显著位点标红
    sig = df["fdr"] < fdr_thresh
    plt.scatter(df.loc[sig, "pos_cum"], -np.log10(df.loc[sig, "fdr"]+1e-300),
                c="red", s=10, alpha=0.7)
    top_sites = df.nsmallest(topN, 'fdr')
    for _, row in top_sites.iterrows():
        plt.text(row['pos_cum'], -np.log10(row['fdr']),f"{row['chr']}:{row['start']}", fontsize=6, rotation=45)
    plt.axhline(-np.log10(fdr_thresh), color="grey", linestyle="--", lw=1)
    plt.xticks(xticks, xticklabels, rotation=90)
    plt.ylabel("-log10(FDR)")
    plt.title("Manhattan plot of %s modification (n = %s)" % (mod_code_to_name.get(base,'Unknown'), df.shape[0]))
    plt.tight_layout()
    plt.savefig(outpng, dpi=300)
    plt.close()

def volcano_plot(df, outpng, base, fdr_thresh=0.05):
    plt.figure(figsize=(6,6))
    plt.scatter(df["log2FC"], -np.log10(df["fdr"]+1e-300),
                c="grey", alpha=0.6, s=10)
    sig = df["fdr"] < fdr_thresh
    plt.scatter(df.loc[sig,"log2FC"], -np.log10(df.loc[sig,"fdr"]+1e-300),
                c="red", alpha=0.7, s=12)
    top_sites = df.nsmallest(topN, 'fdr')
    for _, row in top_sites.iterrows():
        plt.text(row['log2FC'], -np.log10(row['fdr']),f"{row['chr']}:{row['start']}", fontsize=6, rotation=45)
    plt.axhline(-np.log10(fdr_thresh), color="blue", linestyle="--", lw=1)
    plt.xlabel("log2FC (treat/control)")
    plt.ylabel("-log10(FDR)")
    plt.title("Volcano plot of %s modification (n = %s)"  % (mod_code_to_name.get(base,'Unknown'), df.shape[0]))
    plt.tight_layout()
    plt.savefig(outpng, dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Compare differential modification sites between control and treatment.")
    parser.add_argument("control", help="control pileup.bed")
    parser.add_argument("treatment", help="treatment pileup.bed")
    parser.add_argument("-o", "--output", default="diff_modsites.tsv", help="output file")
    parser.add_argument("--volcano", default="volcano.png", help="volcano plot PNG")
    parser.add_argument("--manhattan", default="manhattan.png", help="manhattan plot PNG")
    parser.add_argument("--summary", default="summary.txt", help="summary of sites")
    parser.add_argument("--fdr", type=float, default=0.05, help="FDR threshold")
    parser.add_argument("--depth_filter", type=int, default=0, help="depth threshold")
    parser.add_argument("--mod_fraction", type=float, default=0.0, help="modification fraction threshold")
    parser.add_argument("--folder", help="output folder")
    args = parser.parse_args()

    print("Arguments:")
    for arg, value in vars(args).items():
        print(f"  {arg}: {value}")

    os.makedirs(args.folder, exist_ok=True)

    # 1. load
    ctrl = load_pileup(args.control, args.depth_filter, args.mod_fraction)
    treat = load_pileup(args.treatment, args.depth_filter, args.mod_fraction)

    # per modification
    for base in treat['base'].unique():
        ctrl_base = ctrl.loc[ctrl['base']==base,]
        treat_base = treat.loc[treat['base']==base,]
        # 2. merge
        summary_base_h = open(args.folder + '/' + base + '.' + args.summary,'wt')
        merged = merge_sites(ctrl_base, treat_base,summary_base_h)
        

        # 3. fisher test
        merged["pval"] = merged.apply(fisher_test, axis=1)
        merged["fdr"] = multipletests(merged["pval"], method="fdr_bh")[1]
        # merged["fdr"] = multipletests(merged["pval"], method="bonferroni")[1]

        # 4. calculate effect size
        merged["frac_ctrl"] = merged["mod_count_ctrl"] / (merged["mod_count_ctrl"] + merged["unmod_count_ctrl"]).replace(0, np.nan)
        merged["frac_treat"] = merged["mod_count_treat"] / (merged["mod_count_treat"] + merged["unmod_count_treat"]).replace(0, np.nan)
        merged["delta_frac"] = merged["frac_treat"] - merged["frac_ctrl"]
        merged["log2FC"] = np.log2((merged["frac_treat"]+1e-6) / (merged["frac_ctrl"]+1e-6))

        # 5. save results
        merged.to_csv(args.folder + '/' + base + '.' + args.output, sep="\t", index=False)
        summary_base_h.write('\n')
        summary_base_h.write('treat: %s  vs control: %s ' % (args.treatment, args.control ))
        summary_base_h.write('FDR p-value < 0.05\n')
        up = sum((merged["fdr"]<=0.05) & (merged["log2FC"] > 0 ))
        down = sum((merged["fdr"]<=0.05) & (merged["log2FC"] < 0 ))
        summary_base_h.write('LFC > 0 (up)       : %s, %s\n' % ( up, up/merged.shape[0] ) )
        summary_base_h.write('LFC < 0 (down)     : %s, %s\n' % ( down, down/merged.shape[0] ))
        summary_base_h.close()

        # 6. plots
        volcano_plot(merged, args.folder + '/' + base + '.' + args.volcano, base, fdr_thresh=args.fdr)
        manhattan_plot(merged, args.folder + '/' + base + '.' + args.manhattan, base, fdr_thresh=args.fdr)

if __name__ == "__main__":
    main()
