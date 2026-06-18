#!/usr/bin/env python3
"""Resolve GEO sample names to SRA run accessions and download FASTQ.gz files.

Inputs:
  1) A TSV manifest with columns:
       gse\tpattern\tlabel
     - gse: GEO series accession, e.g. GSE70355
     - pattern: a regular expression matched against sample/title fields
     - label: optional friendly label used in the output folder name

  2) NCBI public GEO/SRA pages.

The script:
  - opens the GEO series page
  - extracts the linked SRA study accession (PRJNA... or SRP...)
  - opens the SRA Run Selector page
  - parses the run table and matches rows by your sample patterns
  - downloads each SRR with fasterq-dump
  - gzips the resulting FASTQ files

Requirements:
  - Python 3.9+
  - requests
  - pandas
  - lxml (or html5lib) for pandas.read_html
  - SRA Toolkit in PATH: fasterq-dump
  - gzip or pigz in PATH (gzip is standard)

Examples:
  python geo_sra_fastq_download.py --manifest targets.tsv --outdir downloads --threads 8

References:
  - GEO series pages list sample titles and an SRA Run Selector.
  - NCBI SRA docs describe using fasterq-dump / prefetch to obtain FASTQ from SRR accessions.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import pandas as pd
import requests
from pysradb import SRAweb
db = SRAweb()


UA = "Mozilla/5.0 (compatible; GEO-SRA-downloader/1.0; +https://www.ncbi.nlm.nih.gov/)"

GEO_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"
# In GEO pages, the SRA Run Selector link usually resolves to the BioProject / Study accession.
SRA_STUDY_URL = "https://www.ncbi.nlm.nih.gov/Traces/study/?acc={acc}"
SRR_URL= "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR24458366"

@dataclass
class Target:
    gse: str
    pattern: str
    label: str
    comparison: str


@dataclass
class Match:
    gse: str
    label: str
    pattern: str
    run: str
    matched_text: str


class DownloadError(RuntimeError):
    pass


def fetch_html(url: str, timeout: int = 60) -> str:
    resp = requests.get(url, headers={"User-Agent": UA}, timeout=timeout)
    resp.raise_for_status()
    return resp.text


def extract_study_accession(geo_html: str) -> Optional[str]:
    # Prefer PRJNA if present; otherwise accept SRP.
    m = re.search(r"\b(PRJNA\d+)\b", geo_html)
    if m:
        return m.group(1)
    m = re.search(r"\b(SRP\d+)\b", geo_html)
    if m:
        return m.group(1)
    return None


def load_manifest(path: Path) -> List[Target]:
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str).fillna("")
    required = {"gse", "pattern", "comparison"}
    if not required.issubset(df.columns):
        raise ValueError(f"Manifest must contain columns: {sorted(required)}")
    if "label" not in df.columns:
        df["label"] = ""
    targets: List[Target] = []
    for _, row in df.iterrows():
        gse = row["gse"].strip()
        pat = row["pattern"].strip()
        label = row["label"].strip() or pat
        comparison = row["comparison"].strip()
        if not gse or not pat:
            continue
        targets.append(Target(gse=gse, pattern=pat, label=label, comparison=comparison))
    return targets


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    return df


def find_run_table(study_html: str) -> pd.DataFrame:
    """Parse the SRA Run Selector page and return the most likely run table.

    NCBI pages occasionally render multiple tables. We select the table with a
    Run/Accession-like column and at least one metadata column.
    """
    tables = pd.read_html(study_html)
    if not tables:
        raise DownloadError("No HTML tables found on SRA study page")

    best = None
    best_score = -1
    for tbl in tables:
        tbl = normalize_columns(tbl)
        cols = [c.lower() for c in tbl.columns]
        has_run = any(c in {"run", "run accession", "run_accession", "srr"} or "run" in c for c in cols)
        has_meta = any(any(k in c for k in ["sample", "title", "bioproject", "biosample", "experiment", "library"]) for c in cols)
        score = int(has_run) * 10 + int(has_meta) * 3 + len(tbl)
        if score > best_score:
            best = tbl
            best_score = score
    if best is None:
        raise DownloadError("Could not identify a run table on the SRA study page")
    return best


def pick_run_column(df: pd.DataFrame) -> str:
    candidates = [c for c in df.columns if re.search(r"(?i)\b(run|srr|accession)\b", str(c))]
    if not candidates:
        raise DownloadError(f"Could not find a run accession column in {list(df.columns)}")
    # Prefer exact-looking columns.
    for preferred in ["Run", "run", "Run accession", "run accession", "Run_Accession", "SRR"]:
        if preferred in df.columns:
            return preferred
    return candidates[0]


def row_text(row: pd.Series) -> str:
    parts = [str(v) for v in row.values if pd.notna(v)]
    return " | ".join(parts)


def match_rows(df: pd.DataFrame, pattern: str, gse: str, label: str) -> List[Match]:
    regex = re.compile(pattern, flags=re.I)
    run_col = pick_run_column(df)
    matches: List[Match] = []
    for _, row in df.iterrows():
        text = row_text(row)
        if regex.search(text):
            run = str(row[run_col]).strip()
            if re.match(r"^SRR\d+$", run):
                matches.append(Match(gse=gse, label=label, pattern=pattern, run=run, matched_text=text))
    return matches


def gzip_fastq(path: Path) -> Path:
    gz = path.with_suffix(path.suffix + ".gz")
    with open(path, "rb") as src, gzip.open(gz, "wb", compresslevel=6) as dst:
        shutil.copyfileobj(src, dst)
    path.unlink()
    return gz


def run_fasterq_dump(run: str, outdir: Path, threads: int) -> List[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    before = set(outdir.glob(f"{run}*.fastq")) | set(outdir.glob(f"{run}*.fq"))

    cmd = [
        "fasterq-dump",
        "--split-files",
        "--threads",
        str(threads),
        "--outdir",
        str(outdir),
        run,
    ]
    subprocess.run(cmd, check=True)

    after = set(outdir.glob(f"{run}*.fastq")) | set(outdir.glob(f"{run}*.fq"))
    new_files = sorted(after - before)
    if not new_files:
        # Some versions write exactly run_1.fastq/run_2.fastq, so we also scan the folder.
        new_files = sorted(outdir.glob(f"{run}*.fastq")) + sorted(outdir.glob(f"{run}*.fq"))

    gz_files: List[Path] = []
    for fq in new_files:
        if fq.suffix in {".fastq", ".fq"}:
            gz_files.append(gzip_fastq(fq))
    return gz_files

def find_gsm_by_regex(html, pattern):
    regex = re.compile(pattern)

    matches = re.findall(
        r'acc=(GSM\d+)[^>]*>.*?</a></td>\s*<td[^>]*>([^<]+)</td>',
        html,
        flags=re.S
    )

    return [
        (gsm, sample_name)
        for gsm, sample_name in matches
        if regex.search(sample_name)
    ]

def priority(name, key):
    if re.search(key, name):
        return 0
    return 1

def ena_to_ftp(url_or_path: str) -> str:
    """
    Convert ENA fasp-style path:
      era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/...
    to:
      ftp://ftp.sra.ebi.ac.uk/vol1/fastq/...
    Also passes through already-ftp/http paths unchanged.
    """

    if not url_or_path or str(url_or_path).lower() == "nan":
        return ""

    s = str(url_or_path).strip()

    if s.startswith("ftp://") or s.startswith("https://") or s.startswith("http://"):
        return s

    if s.startswith("era-fasp@fasp.sra.ebi.ac.uk:"):
        return "ftp://ftp.sra.ebi.ac.uk/" + s.split(":", 1)[1].lstrip("/")

    return s

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, type=Path, help="TSV with columns gse, pattern, label")
    ap.add_argument("--outdir", required=True, type=Path, help="Output directory")
    args = ap.parse_args()

    args.outdir = os.path.abspath(args.outdir)

    targets = load_manifest(args.manifest)
    if not targets:
        raise SystemExit("Manifest is empty")

    all_matches: List[Match] = []

    # Cache by GSE to avoid repeated downloads.
    geo_cache = {}

    for t in targets:
        print('# gse: ' + t.gse + '\tpattern: ' + t.pattern + '\tlabel: ' + t.label + '\tcomparison: ' + t.comparison)
        if t.gse not in geo_cache:
            geo_url = GEO_URL.format(gse=t.gse)
            geo_html = fetch_html(geo_url)
            geo_cache[t.gse] = geo_html
        else:
            geo_html = geo_cache[t.gse]
        gsm_list = find_gsm_by_regex(geo_html, t.pattern)
        if not gsm_list:
            raise DownloadError(f"Could not find GSMs accession on GEO page for {t.gse}")
        prior = t.comparison.split('_vs_')[0]
        gsm_list = sorted(gsm_list, key=lambda x: priority(x[1], prior))
        gsm_list = [
            [item[0], item[1].replace(" ", "_").replace('(','').replace(')','')]
            for item in gsm_list
        ]
        final_outdir = str(args.outdir) + '/' + t.label
        # print(final_outdir + '\t' + str(len(gsm_list)))
        if not os.path.exists(final_outdir + '/data'):
            os.makedirs(final_outdir + '/data')
        seq_config = open(final_outdir + '/rnaseq.config.yaml', 'wt')
        seq_config.write("""
identifier: "{identifier}"
output_dir: "{outdir}"

spike_in: false

cases: """.format(identifier=t.label, outdir=final_outdir))
        for gsm, sample_name in gsm_list:
            # sample_name = sample_name.replace(" ", "_")
            srr = db.gsm_to_srr(gsm).run_accession[0]
            layout = db.sra_metadata(srr, detailed=True)['library_layout'].iloc[0]
            print("# gsm: " + gsm + '\tsample_name: ' + sample_name + '\tsrr: ' + srr + '\tlibray_layout: ' + layout)
            if layout == 'SINGLE':
                fq1 = final_outdir + '/data/' + sample_name + '.fastq.gz'
                fastq_link = ena_to_ftp(db.sra_metadata(srr,detailed=True)['ena_fastq_ftp_1'].iloc[0])
                print("wget -c -nv -O {fq} {fq_link}; echo $?".format(fq=fq1,fq_link=fastq_link))
            else:
                fq1 = final_outdir + '/data/' + sample_name + '.1.fastq.gz'
                fq2 = final_outdir + '/data/' + sample_name + '.2.fastq.gz'
                fastq_link = ena_to_ftp(db.sra_metadata(srr,detailed=True)['ena_fastq_ftp_1'].iloc[0])
                fastq_link2 = ena_to_ftp(db.sra_metadata(srr,detailed=True)['ena_fastq_ftp_2'].iloc[0])
                print("wget -c -nv -O {fq} {fq_link}; echo $?".format(fq=fq1,fq_link=fastq_link))
                print("wget -c -nv -O {fq} {fq_link}; echo $?".format(fq=fq2,fq_link=fastq_link2))
            
## for rnaseq.config.yaml
            seq_config.write("""
    {sample}:
        fastq1: "{fq}" """.format(sample=sample_name,fq=fq1))
            if not layout == 'SINGLE':
                seq_config.write("""
        fastq2: "{fq}"
            """.format(fq=fq2))
        
        seq_config.write("""
comparisons:
    {comp}: {sample_list}

genome: "/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/ERCC92_Spike_In/STAR_mm10_75bp_ucsc/"  # Reference genome
gtf: "/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/ERCC92_Spike_In/mm10_plus_Ercc.gtf"
bit: "/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/ERCC92_Spike_In/mm10_plus_Ercc.2bit"
transcriptome: "/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/ERCC92_Spike_In/mm10_plus_Ercc.transcripts.fa"
tx2gene: "/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/ERCC92_Spike_In/tx2gene.tsv"
fragment_length: 300
effectGenomeSize: 2652866255
scripts_folder: '/ddn/gs1/home/linl7/bin/scripts'
""".format(comp=t.comparison, sample_list = '["' + '", "'.join([item[1] for item in gsm_list]) + '"]'))


    return 0


if __name__ == "__main__":
    raise SystemExit(main())
