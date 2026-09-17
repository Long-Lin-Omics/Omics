#!/usr/bin/env python3

import argparse

import scanpy as sc

import pandas as pd


def normalize_barcodes(barcodes, source):
    """
    Normalize Cell Ranger / AnnData / velocyto loom cell barcodes.

    Parameters
    ----------
    barcodes : pd.Index, pd.Series, list-like
        Input cell barcodes.
    source : {"loom", "adata"}
        Format of the input barcodes.

    Returns
    -------
    clean_barcodes : pd.Index
        Normalized barcodes used for matching.
    raw_barcodes : pd.Index
        Original barcodes, unchanged.
    meta : dict
        Information about the transformations performed.

    Notes
    -----
    loom example:
        possorted_genome_bam_CM0X7:AAATGCTTCCTGCGGAx
        -> AAATGCTTCCTGCGGA

    adata example:
        AAACCAAAGGAGTCGA-1_WT-EB-12h
        -> AAACCAAAGGAGTCGA-1
        -> AAACCAAAGGAGTCGA   # if all well suffixes are -1
    """

    if source not in {"loom", "adata"}:
        raise ValueError("source must be 'loom' or 'adata'")

    # Always preserve the original input
    raw_barcodes = pd.Index(barcodes).astype(str)

    if len(raw_barcodes) == 0:
        return raw_barcodes.copy(), raw_barcodes.copy(), {
            "source": source,
            "removed_sample_suffix": False,
            "well_suffixes": [],
            "removed_well_suffix": False,
            "removed_loom_x": False,
        }

    if source == "loom":
        # Example:
        # possorted_genome_bam_CM0X7:AAATGCTTCCTGCGGAx
        # -> AAATGCTTCCTGCGGAx
        # -> AAATGCTTCCTGCGGA

        clean = raw_barcodes.str.rsplit(":", n=1).str[-1]

        # Remove exactly one trailing x
        has_trailing_x = clean.str.endswith("x")

        if has_trailing_x.all():
            clean = clean.str[:-1]
            removed_loom_x = True
        else:
            # Do not silently modify a mixed-format input
            clean = pd.Index(clean)
            removed_loom_x = False

        meta = {
            "source": source,
            "removed_sample_suffix": False,
            "well_suffixes": [],
            "removed_well_suffix": False,
            "removed_loom_x": removed_loom_x,
        }

    else:  # adata
        # Example:
        # AAACCAAAGGAGTCGA-1_WT-EB-12h
        # -> AAACCAAAGGAGTCGA-1

        clean = raw_barcodes.str.split("_", n=1).str[0]
        removed_sample_suffix = (clean != raw_barcodes).any()

        # Look for the 10x GEM well suffix
        # e.g. AAAC...-1
        well_match = clean.str.extract(r"-(\d+)$", expand=False)

        well_suffixes = sorted(well_match.dropna().unique().tolist())

        # Only remove the suffix when:
        # 1. every barcode has a well suffix
        # 2. all suffixes are exactly "1"
        all_have_well = well_match.notna().all()
        all_well_is_one = all_have_well and (well_match == "1").all()

        if all_well_is_one:
            clean = clean.str.replace(r"-1$", "", regex=True)
            removed_well_suffix = True
        else:
            removed_well_suffix = False

        meta = {
            "source": source,
            "removed_sample_suffix": removed_sample_suffix,
            "well_suffixes": well_suffixes,
            "removed_well_suffix": removed_well_suffix,
            "removed_loom_x": False,
        }

    clean_barcodes = pd.Index(clean.astype(str))

    # Sanity check: normalization should not create duplicates
    n_unique_raw = raw_barcodes.nunique()
    n_unique_clean = clean_barcodes.nunique()

    meta["n_barcodes"] = len(raw_barcodes)
    meta["n_unique_raw"] = n_unique_raw
    meta["n_unique_clean"] = n_unique_clean
    meta["introduced_duplicates"] = n_unique_clean < n_unique_raw

    print('input barcode example: ' + raw_barcodes[0])
    print(meta)
    print('clean barcode example: ' + clean_barcodes[0])

    return clean_barcodes


def main(
    loom_file,
    adata_file,
):

    print(
        f"Reading loom:\n{loom_file}"
    )

    loom = sc.read_loom(
        loom_file,
    )

    print("\n=== LOOM ===")
    print(loom)

    print(
        "\nloom layers:"
    )

    for key in loom.layers.keys():
        print(
            f"  {key}: {loom.layers[key].shape}"
        )

    print(
        "\nFirst loom barcodes:"
    )

    print(
        loom.obs_names[:10].tolist()
    )

    print(
        f"\nReading AnnData:\n{adata_file}"
    )

    adata = sc.read_h5ad(
        adata_file
    )

    print("\n=== RNA AnnData ===")
    print(adata)

    print(
        "\nFirst AnnData barcodes:"
    )

    print(
        adata.obs_names[:10].tolist()
    )

    # ----------------------------------------------------------
    # Remove common velocyto barcode suffix if present.
    # ----------------------------------------------------------

    loom_barcodes=normalize_barcodes( loom.obs_names.astype(str),source='loom')

    adata_barcodes=normalize_barcodes(adata.obs_names.astype(str),source='adata')

    # ----------------------------------------------------------
    # Matching diagnostics
    # ----------------------------------------------------------

    print(
        "\nBarcode overlap using exact names:"
    )

    exact = (
        set(loom_barcodes)
        & set(adata_barcodes)
    )

    print(
        f"  {len(exact):,}"
    )

    if "original_barcode" in adata.obs.columns:

        original = set(
            adata.obs[
                "original_barcode"
            ]
            .astype(str)
        )

        overlap_original = (
            set(loom_barcodes)
            & original
        )

        print(
            "\nBarcode overlap using "
            "`adata.obs['original_barcode']`:"
        )

        print(
            f"  {len(overlap_original):,}"
        )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--loom",
        required=True,
    )

    parser.add_argument(
        "--adata",
        required=True,
    )

    args = parser.parse_args()

    main(
        loom_file=args.loom,
        adata_file=args.adata,
    )