"""
Background (Vollger et al. 2025, Nature Genetics)
---------------------------------------------------
UDN318336 carries a de novo balanced t(X;13) translocation. NHEJ produced two
derivative chromosomes with slightly different breakpoint coordinates at each
junction:
 
    der(X)  junction: chr13:35,480,299 | chrX:24,535,842
    der(13) junction: chr13:35,480,297 | chrX:24,535,860
 
The small coordinate offsets (2 bp on chr13, 18 bp on chrX) reflect the NHEJ
repair, with no microhomology at the junctions.

── Derivative chromosome construction ─────────────────────────────────────────

Reference chromosomes (plus strand, 5'→3'):

  chr13:  |=====A=====|=====B=====|
                   bp13_der13^ ^bp13_derX

  chrX:   |=====C=====|=====D=====|
                   bpX_der13^  ^bpX_derX

Balanced t(X;13) produces two derivative chromosomes. On the forward
(plus) strand of each derivative contig:

  der(13)_fwd:  |---A--->|---rc(C)--->|
                chr13 native   chrX revcomp

  der(X)_fwd:   |---rc(B)--->|---D--->|
                chr13 revcomp   chrX native

Methylation (CpG bigWig):
  Values are stored at the C of each CpG. Reverse-complemented segments
  require flip + left-shift-by-1 (flip_and_shift_meth) to place values
  at the new plus-strand CpG C position.

Fiber-seq accessibility (non-CpG-specific bigWig):
  Reverse-complemented segments require only a simple array reversal
  (no shift), since values are not tied to a strand-specific dinucleotide.

Breakpoint coordinates (0-based half-open, from Vollger et al. 2025):
  der(13) junction: chr13:35,480,297 | chrX:24,535,860
  der(X)  junction: chr13:35,480,299 | chrX:24,535,842
───────────────────────────────────────────────────────────────────────────────
"""

import numpy as np
import pyBigWig
from pysam import FastaFile
from Bio.Seq import Seq

from methylseqnet_repro.paths import haplotypes, pacbio_5mC_tracks, genomes, fiberseq_tracks

# ── Configuration ──────────────────────────────────────────────────────────────
 
INPUT_FASTA  = str(haplotypes / "UDN318336.hg38.HP1.fa")
OUTPUT_FASTA = str(haplotypes / "UDN318336.hg38.derXder13.fa")

INPUT_CPG_BIGWIG = str(pacbio_5mC_tracks / "UDN318336_retinal.aligned_bam_to_cpg_scores.hap1.bw")
OUTPUT_CPG_BIGWIG = str(pacbio_5mC_tracks / "UDN318336_retinal.hg38.derXder13.cpg_scores.hap1.bw")

INPUT_FIBER_BIGWIG = str(fiberseq_tracks / "UDN318336_retinal_trackHub" / "bw" / "hap1.acc.bw")
OUTPUT_FIBER_BIGWIG = str(fiberseq_tracks / "UDN318336_retinal_trackHub" / "bw" / "hap1.derXder13.acc.bw")
 
first_chrom  = "chr13"
second_chrom = "chrX"
first_der    = "der13"
second_der   = "derX"
 
# Per-junction breakpoints from Vollger et al. Fig. 2c.
# 0-based half-open: the value is the first base EXCLUDED from the
# pre-breakpoint segment (equivalently, first base INCLUDED in the
# post-breakpoint segment).
#
# der(13) junction: chr13:35,480,297 | chrX:24,535,860
bp13_der13 = 35_480_297   # end of chr13 segment A on der(13)
bpX_der13  = 24_535_860   # end of chrX segment C on der(13)
#
# der(X) junction: chrX:24,535,842 | chr13:35,480,299
bpX_derX   = 24_535_842   # start of chrX segment D on der(X)
bp13_derX  = 35_480_299   # start of chr13 segment B on der(X)
 
 
# == Helpers ===================================================================
 
def revcomp(seq_str: str) -> str:
    """Reverse complement a DNA string."""
    return str(Seq(seq_str).reverse_complement())
 
 
def write_fasta(path: str, contigs: dict, line_width: int = 80):
    """Write contigs to a FASTA file."""
    with open(path, "w") as f:
        for name, seq in contigs.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + "\n")
 
def load_bigwig_array(bw, chrom: str, start: int, end: int) -> np.ndarray:
    """Load a region from a bigWig as a numpy array. Missing values become NaN."""
    vals = bw.values(chrom, start, end)
    return np.array(vals, dtype=np.float64)
 
 
def flip_and_shift_meth(arr: np.ndarray) -> np.ndarray:
    """
    Reverse a methylation array and shift left by 1 to account for the CpG
    coordinate convention (value stored at the C, which moves by -1 relative
    to a naive array flip when going from plus-strand CG to revcomp plus-strand CG).
 
    After reversal, position len-i-1 holds the value from position i.
    The correct position is len-i-2 (the C of the new CpG).
    Shifting left by 1: new[j] = flipped[j+1], with the last position becoming NaN.
    """
    flipped = arr[::-1].copy()
    shifted = np.empty_like(flipped)
    shifted[:-1] = flipped[1:]
    shifted[-1] = np.nan
    return shifted
 
 
def get_cpg_positions(seq: str) -> set:
    """Return set of 0-based positions of C in CpG dinucleotides (plus strand)."""
    seq_upper = seq.upper()
    positions = set()
    for i in range(len(seq_upper) - 1):
        if seq_upper[i] == 'C' and seq_upper[i + 1] == 'G':
            positions.add(i)
    return positions
 
 
def write_bigwig(path: str, contig_data: dict):
    """
    Write a bigWig file from contig methylation data.
 
    contig_data: dict of {contig_name: (length, meth_array)}
    Only writes intervals where the value is not NaN.
    """
    # Collect contig names and lengths in order
    header = [(name, length) for name, (length, _) in contig_data.items()]
 
    bw_out = pyBigWig.open(path, "w")
    bw_out.addHeader(header)
 
    for name, (length, arr) in contig_data.items():
        # Find non-NaN positions
        valid = ~np.isnan(arr)
        indices = np.where(valid)[0]
 
        if len(indices) == 0:
            continue
 
        # Build bedGraph-style entries (individual 1-bp intervals)
        chroms = [name] * len(indices)
        starts = indices.astype(np.int64).tolist()
        ends = (indices + 1).astype(np.int64).tolist()
        values = arr[indices].astype(np.float64).tolist()
 
        bw_out.addEntries(chroms, starts, ends=ends, values=values)
 
    bw_out.close()
 
 
# == Main ======================================================================
 
def main():
    fasta = FastaFile(INPUT_FASTA)
 
    chr13_len = fasta.get_reference_length(first_chrom)
    chrX_len  = fasta.get_reference_length(second_chrom)
 
    # -- Extract sequence segments ---------------------------------------------
 
    A_seq = fasta.fetch(first_chrom,  0,          bp13_der13)
    C_seq = fasta.fetch(second_chrom, 0,          bpX_der13)
    D_seq = fasta.fetch(second_chrom, bpX_derX,   chrX_len)
    B_seq = fasta.fetch(first_chrom,  bp13_derX,  chr13_len)
 
    print("Segment lengths:")
    print(f"  A  chr13[0:{bp13_der13}]            = {len(A_seq):>12,} bp")
    print(f"  C  chrX[0:{bpX_der13}]              = {len(C_seq):>12,} bp")
    print(f"  D  chrX[{bpX_derX}:{chrX_len}]      = {len(D_seq):>12,} bp")
    print(f"  B  chr13[{bp13_derX}:{chr13_len}]   = {len(B_seq):>12,} bp")
 
    # -- Build FASTA contigs ---------------------------------------------------
 
    der13_fwd_seq = A_seq + revcomp(C_seq)
    derX_fwd_seq  = revcomp(B_seq) + D_seq
    der13_rev_seq = revcomp(der13_fwd_seq)
    derX_rev_seq  = revcomp(derX_fwd_seq)
 
    print(f"\nDerived contig lengths:")
    print(f"  {first_der}_fwd  = {len(der13_fwd_seq):>12,} bp")
    print(f"  {second_der}_fwd = {len(derX_fwd_seq):>12,} bp")
 
    fasta_contigs = {
        f"{first_der}_fwd":  der13_fwd_seq,
        f"{first_der}_rev":  der13_rev_seq,
        f"{second_der}_fwd": derX_fwd_seq,
        f"{second_der}_rev": derX_rev_seq,
    }
 
    write_fasta(OUTPUT_FASTA, fasta_contigs)
    print(f"\nWrote FASTA to {OUTPUT_FASTA}")
 
    # -- Load methylation arrays -----------------------------------------------
 
    bw_in = pyBigWig.open(INPUT_CPG_BIGWIG)
    
    print("\nLoading methylation from bigWig...")
 
    # Segments for der(13) = chr13[0:bp13_der13] + revcomp(chrX[0:bpX_der13])
    meth_chr13_pre  = load_bigwig_array(bw_in, first_chrom,  0,         bp13_der13)
    meth_chrX_pre   = load_bigwig_array(bw_in, second_chrom, 0,         bpX_der13)
 
    # Segments for der(X) = revcomp(chr13[bp13_derX:end]) + chrX[bpX_derX:end]
    meth_chrX_post  = load_bigwig_array(bw_in, second_chrom, bpX_derX,  chrX_len)
    meth_chr13_post = load_bigwig_array(bw_in, first_chrom,  bp13_derX, chr13_len)
 
    bw_in.close()
 
    print(f"  chr13 pre  (der13 native):   {np.sum(~np.isnan(meth_chr13_pre)):>10,} values in {len(meth_chr13_pre):,} bp")
    print(f"  chrX  pre  (der13 revcomp):  {np.sum(~np.isnan(meth_chrX_pre)):>10,} values in {len(meth_chrX_pre):,} bp")
    print(f"  chrX  post (derX native):    {np.sum(~np.isnan(meth_chrX_post)):>10,} values in {len(meth_chrX_post):,} bp")
    print(f"  chr13 post (derX revcomp):   {np.sum(~np.isnan(meth_chr13_post)):>10,} values in {len(meth_chr13_post):,} bp")
 
    # -- Assemble derivative methylation arrays --------------------------------
 
    # der13_fwd = chr13_pre (as-is) + revcomp(chrX_pre) (flip+shift)
    der13_fwd_meth = np.concatenate([meth_chr13_pre, flip_and_shift_meth(meth_chrX_pre)])
 
    # derX_fwd = revcomp(chr13_post) (flip+shift) + chrX_post (as-is)
    derX_fwd_meth = np.concatenate([flip_and_shift_meth(meth_chr13_post), meth_chrX_post])
 
    # _rev contigs: flip+shift the entire _fwd array
    der13_rev_meth = flip_and_shift_meth(der13_fwd_meth)
    derX_rev_meth  = flip_and_shift_meth(derX_fwd_meth)
 
    # -- Validate: non-NaN values should only occur at CpG positions -----------
 
    print("\nValidating methylation positions against CpG sites in FASTA...")
 
    validation_pairs = [
        (f"{first_der}_fwd",  der13_fwd_seq,  der13_fwd_meth),
        (f"{first_der}_rev",  der13_rev_seq,  der13_rev_meth),
        (f"{second_der}_fwd", derX_fwd_seq,   derX_fwd_meth),
        (f"{second_der}_rev", derX_rev_seq,   derX_rev_meth),
    ]
 
    all_valid = True
    for name, seq, meth in validation_pairs:
        cpg_pos = get_cpg_positions(seq)
        meth_pos = set(np.where(~np.isnan(meth))[0])
 
        off_target = meth_pos - cpg_pos
        on_target  = meth_pos & cpg_pos
        cpg_no_val = cpg_pos - meth_pos
 
        print(f"\n  {name}:")
        print(f"    CpG sites in sequence:          {len(cpg_pos):>10,}")
        print(f"    Non-NaN methylation positions:   {len(meth_pos):>10,}")
        print(f"    Methylation at CpG sites:        {len(on_target):>10,}")
        print(f"    CpG sites without methylation:   {len(cpg_no_val):>10,}")
        print(f"    Methylation at non-CpG sites:    {len(off_target):>10,}  {'*** ERROR ***' if off_target else '(OK)'}")
 
        if off_target:
            all_valid = False
            # Show a few examples
            examples = sorted(off_target)[:5]
            for pos in examples:
                local_seq = seq[max(0, pos-2):pos+3].upper()
                print(f"      pos {pos}: context ...{local_seq}... val={meth[pos]:.3f}")
    
    if all_valid:
        print("\n  All contigs pass validation!")
    else:
        print("\n  *** VALIDATION FAILURES DETECTED -- check CpG shift logic ***")
 
    # -- Write bigWig ----------------------------------------------------------
 
    contig_data = {
        f"{first_der}_fwd":  (len(der13_fwd_seq),  der13_fwd_meth),
        f"{first_der}_rev":  (len(der13_rev_seq),  der13_rev_meth),
        f"{second_der}_fwd": (len(derX_fwd_seq),   derX_fwd_meth),
        f"{second_der}_rev": (len(derX_rev_seq),   derX_rev_meth),
    }
 
    write_bigwig(OUTPUT_CPG_BIGWIG, contig_data)
    print(f"\nWrote bigWig to {OUTPUT_CPG_BIGWIG}")
 
    # -- Fiber-seq accessibility bigWig ----------------------------------------

    print("\nLoading Fiber-seq accessibility from bigWig...")
    bw_fiber_in = pyBigWig.open(INPUT_FIBER_BIGWIG)

    fiber_chr13_pre  = load_bigwig_array(bw_fiber_in, first_chrom,  0,         bp13_der13)
    fiber_chrX_pre   = load_bigwig_array(bw_fiber_in, second_chrom, 0,         bpX_der13)
    fiber_chrX_post  = load_bigwig_array(bw_fiber_in, second_chrom, bpX_derX,  chrX_len)
    fiber_chr13_post = load_bigwig_array(bw_fiber_in, first_chrom,  bp13_derX, chr13_len)

    bw_fiber_in.close()

    print(f"  chr13 pre  (der13 native):   {np.sum(~np.isnan(fiber_chr13_pre)):>10,} values in {len(fiber_chr13_pre):,} bp")
    print(f"  chrX  pre  (der13 revcomp):  {np.sum(~np.isnan(fiber_chrX_pre)):>10,} values in {len(fiber_chrX_pre):,} bp")
    print(f"  chrX  post (derX native):    {np.sum(~np.isnan(fiber_chrX_post)):>10,} values in {len(fiber_chrX_post):,} bp")
    print(f"  chr13 post (derX revcomp):   {np.sum(~np.isnan(fiber_chr13_post)):>10,} values in {len(fiber_chr13_post):,} bp")

    # Simple flip (no shift) for non-CpG-specific tracks
    def flip_array(arr: np.ndarray) -> np.ndarray:
        return arr[::-1].copy()

    der13_fwd_fiber = np.concatenate([fiber_chr13_pre, flip_array(fiber_chrX_pre)])
    derX_fwd_fiber  = np.concatenate([flip_array(fiber_chr13_post), fiber_chrX_post])
    der13_rev_fiber = flip_array(der13_fwd_fiber)
    derX_rev_fiber  = flip_array(derX_fwd_fiber)

    fiber_contig_data = {
        f"{first_der}_fwd":  (len(der13_fwd_seq),  der13_fwd_fiber),
        f"{first_der}_rev":  (len(der13_rev_seq),  der13_rev_fiber),
        f"{second_der}_fwd": (len(derX_fwd_seq),   derX_fwd_fiber),
        f"{second_der}_rev": (len(derX_rev_seq),   derX_rev_fiber),
    }

    write_bigwig(OUTPUT_FIBER_BIGWIG, fiber_contig_data)
    print(f"Wrote Fiber-seq bigWig to {OUTPUT_FIBER_BIGWIG}")

if __name__ == "__main__":
    main()