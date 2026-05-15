#!/usr/bin/env python3
"""
Add stage labels to SQANTI3 classification output.

The script links SQANTI3 isoforms/PBIDs to barcodes in the abundance file,
then maps barcodes to developmental stage labels. If an isoform is supported
by barcodes from multiple stages, the most frequent stage is assigned.
"""

import argparse
import os
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--classification", required=True, help="SQANTI3 classification.txt")
    parser.add_argument("--abundance", required=True, help="Mapped abundance file containing pbid and cell_barcodes")
    parser.add_argument("--barcode-stage", required=True, help="Two-column barcode to stage table")
    parser.add_argument("--output", required=True, help="Output classification file with stage column")
    return parser.parse_args()


def main():
    args = parse_args()

    for f in [args.classification, args.abundance, args.barcode_stage]:
        if not os.path.isfile(f):
            raise FileNotFoundError(f)

    cls = pd.read_csv(args.classification, sep="\t", dtype=str)
    iso_col = cls.columns[0]

    ab = pd.read_csv(args.abundance, sep="\t", comment="#", dtype=str)
    if "pbid" not in ab.columns or "cell_barcodes" not in ab.columns:
        raise ValueError("Abundance file must contain columns: pbid, cell_barcodes")

    ab["cell_barcodes"] = ab["cell_barcodes"].fillna("")

    exploded = (
        ab[["pbid", "cell_barcodes"]]
        .assign(cell_barcodes=lambda df: df["cell_barcodes"].str.split(","))
        .explode("cell_barcodes")
        .rename(columns={"cell_barcodes": "barcode", "pbid": iso_col})
    )
    exploded["barcode"] = exploded["barcode"].str.strip()

    bc2stage = pd.read_csv(args.barcode_stage, sep="\t", dtype=str)

    if "barcode_rc" in bc2stage.columns and "barcode" not in bc2stage.columns:
        bc2stage = bc2stage.rename(columns={"barcode_rc": "barcode"})

    required_cols = {"barcode", "stage"}
    missing = required_cols - set(bc2stage.columns)
    if missing:
        raise ValueError(f"Barcode-stage file missing columns: {', '.join(missing)}")

    merged = exploded.merge(bc2stage[["barcode", "stage"]], on="barcode", how="left")

    def pick_stage(series):
        s = series.dropna()
        if len(s) == 0:
            return None
        return s.value_counts().idxmax()

    stage_per_iso = (
        merged.groupby(iso_col)["stage"]
        .agg(pick_stage)
        .reset_index()
    )

    out = cls.merge(stage_per_iso, on=iso_col, how="left")

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
