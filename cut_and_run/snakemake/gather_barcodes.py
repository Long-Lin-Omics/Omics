import pandas as pd

input_files = snakemake.input
output_file = snakemake.output[0]

rows = []

for f in sorted(input_files):
    with open(f) as fh:
        lines = [l.strip() for l in fh if l.strip()]

    # 第一行就是 sample name
    sample = lines[0]
    marker = ""
    modifications = {}
    percent_total_str = None  # "0.374%"
    percent_total = None      # 0.00374

    # 从第二行开始解析
    for line in lines[1:]:
        # Marker 行
        if line.startswith("Marker:"):
            marker = line.split(":", 1)[1].strip()
            continue

        # % total barcode reads 行
        if line.startswith("% total barcode reads"):
            val = line.split(":", 1)[1].strip()
            percent_total_str = val
            continue

        # 其他三列格式: 修饰名 counts ratio
        parts = line.split()
        
        if line.startswith("Uniq align reads"):
            sample_sequencing_depth = line.split(":", 1)[1].strip()
            continue

        if line.startswith("Total barcode reads"):
            continue

        if len(parts) == 3:
            mod = parts[0]          # H3K27me1 ...
            ratio = f"{float(parts[2]):.3g}"  # 0.15548...
            modifications[mod] = ratio

    row = {
        "sample": sample,
        "marker": marker,
        **modifications,
        "percent_total_barcode_reads": percent_total_str,
        "sample_sequencing_depth": sample_sequencing_depth
        
    }
    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(output_file, sep="\t", index=False)
