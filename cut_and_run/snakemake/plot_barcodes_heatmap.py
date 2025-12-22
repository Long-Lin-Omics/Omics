import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

table_file = snakemake.input[0]
out_png = snakemake.output[0]

df = pd.read_csv(table_file, sep="\t")

# 选出修饰列：排除 meta 列
meta_cols = ["sample", "marker",
             "percent_total_barcode_reads",
             "sample_sequencing_depth"]
mod_cols = [c for c in df.columns if c not in meta_cols]

# 行索引可以用 sample，也可以拼上 marker
df_indexed = df.set_index("sample")[mod_cols]

plt.figure(figsize=(10, max(4, 0.4 * df_indexed.shape[0])))  # 高度随样本数略调

sns.heatmap(
    df_indexed,
    cmap="viridis",       # 颜色方案可自己换
    annot=True,           # 显示数字
    fmt=".3g",            # 有效数字（3 位）："0.00374" 会变 "0.00374"
    annot_kws={"size": 6}, 
    linewidths=0.3,
    linecolor="white",
    cbar_kws={"label": "ratio"}
)

plt.xlabel("Modification with barcodes")
plt.ylabel("Sample")
plt.title("Histone modification ratios per sample")

plt.tight_layout()
plt.savefig(out_png, dpi=300)
plt.close()
