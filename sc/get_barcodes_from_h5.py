#!/usr/bin/env python3

import os
import argparse

import h5py


def read_h5_barcodes(h5_file):
    """
    Extract 10x cell barcodes from filtered_feature_bc_matrix.h5.

    Supports standard 10x HDF5 layout:
        matrix/barcodes
    """

    with h5py.File(h5_file, "r") as f:

        if "matrix" not in f:
            raise RuntimeError(
                f"`matrix` group not found in:\n{h5_file}"
            )

        if "barcodes" not in f["matrix"]:
            raise RuntimeError(
                f"`matrix/barcodes` not found in:\n{h5_file}"
            )

        barcodes = f["matrix"]["barcodes"][:]

    decoded = []

    for x in barcodes:
        if isinstance(x, bytes):
            x = x.decode("utf-8")
        else:
            x = str(x)

        decoded.append(x)

    return decoded


def main(
    h5_file,
    outfile,
):

    if not os.path.exists(h5_file):
        raise FileNotFoundError(
            h5_file
        )

    barcodes = read_h5_barcodes(
        h5_file
    )

    os.makedirs(
        os.path.dirname(
            os.path.abspath(outfile)
        ),
        exist_ok=True,
    )

    with open(outfile, "w") as f:

        for barcode in barcodes:
            f.write(
                barcode + "\n"
            )

    print(
        f"Input H5 : {h5_file}"
    )

    print(
        f"Barcodes : {len(barcodes):,}"
    )

    print(
        f"Output   : {outfile}"
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--h5",
        required=True,
    )

    parser.add_argument(
        "--out",
        required=True,
    )

    args = parser.parse_args()

    main(
        h5_file=args.h5,
        outfile=args.out,
    )