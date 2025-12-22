#!/usr/bin/env python3
import sys
import pysam
import re
from multiprocessing import Pool
from collections import defaultdict

def normalize_chrom_name(chrom, bam_refs):
    # 1. direct match
    if chrom in bam_refs:
        return chrom
    # 2. add 'chr'
    if "chr" + chrom in bam_refs:
        return "chr" + chrom
    # 3. remove 'chr'
    if chrom.startswith("chr") and chrom[3:] in bam_refs:
        return chrom[3:]
    # 4. strip .数字 后匹配
    chrom_base = re.sub(r"\.\d+$", "", chrom)
    for bchr in bam_refs:
        if chrom_base in bchr:
            return bchr
    if "chr" + chrom_base in bam_refs:
        return "chr" + chrom_base
    # no match
    return None

def parse_gtf(gtf_file, bam_references):
    """
    Parse GTF and return gene_exons dict.
    Will normalize chromosome names to match BAM references.
    Also report BAM chromosomes not found in GTF.
    """
    bam_refs = set(bam_references)
    matched_bam_chroms = set()
    gene_exons = defaultdict(list)

    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] != "exon":
                continue
            chrom, start, end, info = fields[0], int(fields[3]) - 1, int(fields[4]), fields[8]
            gene_id = None
            for kv in info.split(";"):
                if "gene_id" in kv:
                    gene_id = kv.strip().split(" ")[1].replace('"', "")
                    break
            if gene_id is None:
                continue
            chrom_norm = normalize_chrom_name(chrom,bam_refs)
            if chrom_norm:
                gene_exons.setdefault(gene_id, {}).setdefault(chrom_norm, []).append((start, end))
                matched_bam_chroms.add(chrom_norm)
    unmapped_bam_chr = bam_refs - matched_bam_chroms
    if unmapped_bam_chr:
        print("\n[BAM chroms without GTF match]:")
        for bchr in sorted(unmapped_bam_chr):
            print(" ", bchr)

    return gene_exons


def classify_read(aln, exons):
    """
    简单分类逻辑：
    - 如果 CIGAR 里有 N -> MATURE
    - 如果无 N，但 read 覆盖了基因内一个 intron -> IR
    - 其他情况 -> UNSURE
    """
    if aln.is_unmapped:
        return "UNMAPPED_but_overlapped"

    # CIGAR op: 3 = N (skipped region, splice junction)
    # has_splice = any(op == 3 for op, _ in aln.cigartuples or [])

    read_start = aln.reference_start
    read_end = aln.reference_end

    # assigned_gene = None
    # for gene_id, exons in gene_exons.items():
    #     for (s, e) in exons:
    #         if not (read_end < s or read_start > e):  # overlap
    #             assigned_gene = gene_id
    #             break
    #     if assigned_gene:
    #         break

    # if assigned_gene is None:
    #     return "UNSURE"

    # if has_splice:
        # return "MATURE"

    # 检查 IR: read 覆盖 intron 区间
    exons = sorted(exons)
    for i in range(len(exons) - 1):
        intron_start = exons[i][1] + 1
        intron_end = exons[i+1][0] - 1
        if read_start <= intron_start and read_end >= intron_end:
            return "IR"

    return "UNSURE"


def process_region(chrom, start, end, bam_file, gene_exons):
    """处理一个染色体区间"""
    bam = pysam.AlignmentFile(bam_file, "rb")
    results = []
    for read in bam.fetch(chrom, start, end):
        no_exon_anchor = True
        if read.is_unmapped:
            label = 'UNMAPPED'
            no_exon_anchor = False
            results.append((read.query_name, label))
            continue
        elif any(op == 3 for op, _ in read.cigartuples or []):
            label = "MATURE"
            no_exon_anchor = False
            results.append((read.query_name, label))
            continue
        else:
            for gene, chrom_dict in gene_exons.items():
                if chrom not in chrom_dict:
                    continue
                if any(s <= read.reference_start < e or s < read.reference_end <= e for s, e in chrom_dict[chrom]):
                    label = classify_read(read, chrom_dict[chrom])
                    no_exon_anchor = False
                    results.append((read.query_name, label))
                    break
        if no_exon_anchor:
            label = 'No_exon_anchor'
            results.append((read.query_name, label))
    bam.close()
    return results


def main():
    if len(sys.argv) != 6:
        sys.stderr.write("Usage: python3 IR_from_miminap.py <bam> <gtf> <output> <nproc> <block_size>\n")
        sys.exit(1)

    bam_file, gtf_file, output_file = sys.argv[1:4]
    nproc = int(sys.argv[4])
    block_size = int(sys.argv[5])

    sys.stderr.write("[INFO] Parsing GTF...\n")

    bam = pysam.AlignmentFile(bam_file, "rb")
    gene_exons = parse_gtf(gtf_file, bam.references)
    tasks = []
    for chrom, length in zip(bam.references, bam.lengths):
        for start in range(0, length, block_size):
            end = min(start + block_size, length)
            tasks.append((chrom, start, end, bam_file, gene_exons))
    bam.close()

    sys.stderr.write(f"[INFO] Created {len(tasks)} tasks, running with {nproc} processes...\n")

    with Pool(processes=nproc) as pool:
        results = pool.starmap(process_region, tasks)

    # flatten results
    with open(output_file, "w") as out:
        for res in results:
            for rname, label in res:
                out.write(f"{rname}\t{label}\n")

    sys.stderr.write(f"[INFO] Done! Results written to {output_file}\n")


if __name__ == "__main__":
    main()
