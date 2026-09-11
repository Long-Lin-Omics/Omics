#!/usr/bin/env python3

import argparse

import scanpy as sc


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

    loom_barcodes = (
        loom.obs_names
        .astype(str)
    )

    adata_barcodes = (
        adata.obs_names
        .astype(str)
    )

    loom_simple = (
        loom_barcodes
        .str.replace(
            ":",
            "-",
            regex=False,
        )
    )

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