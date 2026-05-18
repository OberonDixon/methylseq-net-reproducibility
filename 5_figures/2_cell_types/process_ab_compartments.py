### Load A/B compartment data from Rao et al. 2014, liftover coordinates, and build query index
# Requirements: pip install pyliftover openpyxl pandas

import pandas as pd
import numpy as np
from pyliftover import LiftOver
from tqdm.auto import tqdm

# ============================================================================
# CONFIGURE THESE THREE THINGS
# ============================================================================

# 1. Path to the xlsx file (Rao et al. 2014 Supplementary Table S2)
XLSX_PATH = "NIHMS828671-supplement-2.xlsx"

# 2. Path to chain file for liftover (e.g., hg19ToHg38.over.chain.gz)
#    Download from: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
#    Set to None to skip liftover (keep hg19 coordinates)
CHAIN_PATH = "/clusterfs/nilah/oberon/genomes/hg19ToHg38.over.chain.gz"  # or None

# 3. Column-to-cell-type mapping. Every column in the spreadsheet must have an entry.
#    Values are your desired display names; keys are the original column headers.
CELL_TYPE_MAP = {
    "GM12878": "GM12878",   # lymphoblastoid
    "H1":     "H1-hESC",   # human embryonic stem cell
    "IMR90":  "IMR90",      # fetal lung fibroblast
    "MES":    "MES",        # mesendoderm
    "MSC":    "MSC",        # mesenchymal stem cell
    "NPC":    "NPC",        # neural progenitor cell
    "TRO":    "TRO",        # trophoblast
    "AD":     "adrenal",
    "BL":     "bladder",
    "CO":     "Neuron",     # dorsolateral prefrontal cortex
    "HC":     "hippocampus",
    "LG":     "Lung-Alveolar-Epithel",
    "OV":     "ovary",
    "PA":     "pancreas",
    "PO":     "Skeletal-Muscle",      # psoas muscle
    "SB":     "Small-int-Epithel",
    "AO":     "Coronary-Artery-Smooth-Muscle",
    "RV":     "right_ventricle",
    "LV":     "Heart-Fibroblasts",
    "LI":     "Liver-Hepatocytes",
    "SX":     "Pancreas-Acinar",
}

# ============================================================================
# LOAD & PROCESS
# ============================================================================

# Load the A/B labels sheet
df = pd.read_excel(XLSX_PATH, sheet_name="A_B LABELS")

# Fix chromosome names: integers -> "chrN", 23 -> "chrX"
chr_map = {i: f"chr{i}" for i in range(1, 23)}
chr_map[23] = "chrX"
df["chr"] = df["chr"].map(chr_map)

# Rename cell-type columns
missing = set(CELL_TYPE_MAP.keys()) - set(df.columns)
extra = set(df.columns) - set(CELL_TYPE_MAP.keys()) - {"chr", "start", "end"}
if missing:
    raise ValueError(f"CELL_TYPE_MAP keys not found in spreadsheet: {missing}")
if extra:
    print(f"⚠️  Columns in spreadsheet without mapping (will be dropped): {extra}")
df = df.rename(columns=CELL_TYPE_MAP)
cell_types = list(CELL_TYPE_MAP.values())

# Replace string 'NA' with actual NaN
for ct in cell_types:
    df[ct] = df[ct].replace("NA", np.nan)

# ============================================================================
# LIFTOVER (hg19 -> target assembly)
# ============================================================================
if CHAIN_PATH is not None:
    print(f"Lifting over coordinates using {CHAIN_PATH} ...")
    lo = LiftOver(CHAIN_PATH)

    new_chroms, new_starts, new_ends, kept = [], [], [], []
    for idx, row in df.iterrows():
        chrom = row["chr"]
        # Lift the start (0-based for pyliftover)
        res_start = lo.convert_coordinate(chrom, int(row["start"] - 1))
        # Lift the end
        res_end = lo.convert_coordinate(chrom, int(row["end"]))

        if res_start and res_end:
            new_chr_s, new_pos_s = res_start[0][0], res_start[0][1]
            new_chr_e, new_pos_e = res_end[0][0], res_end[0][1]
            # Only keep if both endpoints map to same chromosome
            if new_chr_s == new_chr_e:
                new_chroms.append(new_chr_s)
                new_starts.append(new_pos_s + 1)  # back to 1-based
                new_ends.append(new_pos_e)
                kept.append(idx)
                continue
        # If liftover failed, drop this bin
        pass

    df = df.loc[kept].copy()
    df["chr"] = new_chroms
    df["start"] = new_starts
    df["end"] = new_ends
    print(f"  Liftover: kept {len(kept)}/{len(kept) + (len(chr_map) * 0 + len(df) - len(df))} bins "
          f"({len(df)} rows)")
else:
    print("Skipping liftover (keeping hg19 coordinates)")

# ============================================================================
# BUILD QUERY INDEX: (chrom, position) -> compartment label per cell type
# ============================================================================

# Store as an IntervalIndex-based structure for fast overlap queries
# We build one sorted dataframe indexed by (chr, start) for bisect-based lookup

df = df.sort_values(["chr", "start"]).reset_index(drop=True)

print(df)

# Precompute per-chromosome arrays for fast numpy-based lookup
_chrom_data = {}
for chrom, grp in df.groupby("chr"):
    _chrom_data[chrom] = {
        "starts": grp["start"].values,
        "ends":   grp["end"].values,
        "labels": grp[cell_types].values,  # shape: (n_bins, n_cell_types)
    }
_cell_type_idx = {ct: i for i, ct in enumerate(cell_types)}


def query_compartment(chrom: str, position: int, cell_type: str = None) -> dict | str:
    """
    Query A/B compartment for a genomic position.

    Parameters
    ----------
    chrom : str
        Chromosome name, e.g. "chr1"
    position : int
        1-based genomic coordinate
    cell_type : str, optional
        If provided, return just the compartment for that cell type ("A", "B", or "NA").
        If None, return a dict of {cell_type: compartment} for all cell types.

    Returns
    -------
    dict or str
        Compartment assignments. Returns "NA" / dict of "NA"s if position not in any bin.
    """
    if chrom not in _chrom_data:
        if cell_type:
            return "NA"
        return {ct: "NA" for ct in cell_types}

    cd = _chrom_data[chrom]
    # Binary search: find bins where start <= position
    idx = np.searchsorted(cd["starts"], position, side="right") - 1
    if idx >= 0 and cd["starts"][idx] <= position <= cd["ends"][idx]:
        if cell_type:
            ci = _cell_type_idx.get(cell_type)
            if ci is None:
                return "NA"
            val = cd["labels"][idx, ci]
            return val if pd.notna(val) else "NA"
        else:
            row = cd["labels"][idx]
            return {ct: (row[i] if pd.notna(row[i]) else "NA") for i, ct in enumerate(cell_types)}
    else:
        if cell_type:
            return "NA"
        return {ct: "NA" for ct in cell_types}


# ============================================================================
# USAGE EXAMPLES
# ============================================================================

# Query a specific cell type:
label = query_compartment("chr1", 2_500_000, "Neuron")