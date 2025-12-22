import pysam
import csv
from collections import defaultdict
from scipy.stats import fisher_exact

def parse_modbam(bam_path, threshold=0.5):
    """
    从 modBAM 统计每个位点的 修饰/未修饰 read 数
    threshold: ML 概率阈值 (0~1)，默认 >=0.5 当作修饰
    返回: dict[(chrom, pos, mod_type)] -> [modified, unmodified]
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    site_counts = defaultdict(lambda: [0, 0])  # (modified, unmodified)
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        try:
            mm_tags = read.get_tag("MM").split(";")  # 可能有多个修饰
            ml = read.get_tag("ML")
        except KeyError:
            continue
        seq = read.query_sequence
        ref_pos = read.get_reference_positions(full_length=True)
        probs = [p / 255.0 for p in ml]
        # 遍历每个 MM 子标签
        for mm in mm_tags:
            if not mm:
                continue
            # e.g. "C+m,0,0,5" 或 "A+a,10,0,0"
            mod_type = mm.split(",")[0]  # "C+m" / "A+a" / "T+17802"
            # 找到该修饰对应的碱基
            base = mod_type.split("+")[0]  # "C" / "A" / "T"
            base_positions = [i for i, b in enumerate(seq) if b == base]
    # wrong here
            for ci, prob in zip(base_positions, probs):
                ref_coord = ref_pos[ci]
                if ref_coord is None:
                    continue
                chrom = bam.get_reference_name(read.reference_id)
                key = (chrom, ref_coord, mod_type)
                if prob >= threshold:
                    site_counts[key][0] += 1
                else:
                    site_counts[key][1] += 1
    bam.close()
    return site_counts

