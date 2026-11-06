"""
Annotate regions with ENCODE SCREEN promoter / enhancer cCREs.

For each region in REGIONS_BED, the region center ((start + end) // 2) is
tested for overlap with promoter-like (PLS), proximal enhancer-like (pELS),
and distal enhancer-like (dELS) cCREs. Output is a TSV with the original
region columns plus boolean columns:
    promoter, proximal_enhancer, distal_enhancer, both, neither
where both = promoter & (proximal_enhancer | distal_enhancer) and
neither = not (promoter | proximal_enhancer | distal_enhancer).

Note: SCREEN cCREs are non-overlapping and each has exactly one class, so a
single center base can fall in at most one cCRE. With these files `both` is
therefore always False; it is kept for schema completeness.

Usage:
    python annotate_regions_screen.py GRCh38-PLS.bed GRCh38-ELS.bed regions.bed out_dir
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def read_bed(path, names=None):
    df = pd.read_csv(path, sep="\t", header=None, comment="#", low_memory=False)
    if names is None:
        names = ["chrom", "start", "end"] + [f"col{i}" for i in range(4, df.shape[1] + 1)]
    df.columns = names[: df.shape[1]]
    return df


def read_ccres(path):
    return read_bed(path, names=["chrom", "start", "end", "dhs_id", "ccre_id", "ccre_class"])


def center_in(intervals, chroms, centers):
    """Boolean array: is each center inside any interval on its chromosome?

    Assumes intervals do not overlap one another (true for SCREEN cCREs), so the
    only candidate is the interval with the largest start <= center.
    """
    hits = np.zeros(len(centers), dtype=bool)
    for chrom, grp in intervals.groupby("chrom"):
        mask = chroms == chrom
        if not mask.any():
            continue
        grp = grp.sort_values("start")
        starts = grp["start"].to_numpy()
        ends = grp["end"].to_numpy()
        c = centers[mask]
        idx = np.searchsorted(starts, c, side="right") - 1
        valid = idx >= 0
        hit = np.zeros(len(c), dtype=bool)
        hit[valid] = c[valid] < ends[idx[valid]]
        hits[mask] = hit
    return hits


def output_path(regions_path, out_dir):
    name = Path(regions_path).name
    gz = name.endswith(".gz")
    if gz:
        name = name[:-3]
    stem = name[:-4] if name.endswith(".bed") else Path(name).stem
    return Path(out_dir) / f"{stem}.screen_annotated.tsv{'.gz' if gz else ''}"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("promoters", help="SCREEN promoter BED (GRCh38-PLS.bed)")
    parser.add_argument("enhancers", help="SCREEN enhancer BED (GRCh38-ELS.bed)")
    parser.add_argument("regions", help="Regions BED (.bed or .bed.gz)")
    parser.add_argument("out_dir", help="Output directory (created if missing)")
    args = parser.parse_args()

    regions = read_bed(args.regions)
    pls = read_ccres(args.promoters)
    els = read_ccres(args.enhancers)
    pels = els[els["ccre_class"].str.startswith("pELS")]
    dels = els[els["ccre_class"].str.startswith("dELS")]
    print(f"Loaded {len(regions)} regions, {len(pls)} PLS, {len(pels)} pELS, {len(dels)} dELS",
          file=sys.stderr)

    chroms = regions["chrom"].astype(str).to_numpy()
    centers = ((regions["start"] + regions["end"]) // 2).to_numpy()

    out = regions.copy()
    out["promoter"] = center_in(pls, chroms, centers)
    out["proximal_enhancer"] = center_in(pels, chroms, centers)
    out["distal_enhancer"] = center_in(dels, chroms, centers)
    any_enh = out["proximal_enhancer"] | out["distal_enhancer"]
    out["both"] = out["promoter"] & any_enh
    out["neither"] = ~(out["promoter"] | any_enh)

    out_path = output_path(args.regions, args.out_dir)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)

    for col in ["promoter", "proximal_enhancer", "distal_enhancer", "both", "neither"]:
        print(f"{col}: {int(out[col].sum())}", file=sys.stderr)
    print(f"Wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
