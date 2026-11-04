"""
Orient the reads of a Kinnex Iso-Seq BAM by their segment-adapter tags.

Why: the GM12878 Kinnex BAM never went through `isoseq refine`, so about 37% of its reads carry a strand label opposite to
their true transcript strand (checked with the splice motifs, on chr22). Isoseq collapse takes the transcript strand, and so
the TSS, from that label. The segment-adapter IDs in the dl and dr tags give the direction of the array: dr - dl = +1 means
the read is correct (0.1% disagree with the splice motif), dr - dl = -1 means it is flipped (85.5% disagree). Orienting by
this tag alone still leaves about 6% of spliced reads wrong (the 14.5% of descending reads that are in fact fine), and it does
not replicate the primer and polyA filtering of refine.

What: toggle the strand flag (0x10) of every mapped read with dr < dl and tag it fo:i:1. Position, CIGAR and SEQ are untouched,
and all other tags (HP, PS, ...) are kept. The MM/ML methylation tags stay keyed to the original read orientation, so do not use
the output for Fiber-seq analysis. The output BAM is indexed. Reads without dl and dr tags are written unchanged.

Usage: python orient_bam_dldr.py input.bam output.bam
"""
import sys
import time

import pysam


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    inp, outp = sys.argv[1:3]
    t0 = time.time()
    n = flipped = notag = 0
    with pysam.AlignmentFile(inp, "rb", threads=3) as fin, pysam.AlignmentFile(outp, "wb", template=fin, threads=3) as fout:
        for r in fin:
            n += 1
            if not r.is_unmapped:
                if r.has_tag("dl") and r.has_tag("dr"):
                    if r.get_tag("dr") < r.get_tag("dl"):
                        r.flag = r.flag ^ 0x10
                        r.set_tag("fo", 1, "i")
                        flipped += 1
                else:
                    notag += 1
            fout.write(r)
            if n % 500000 == 0:
                print(f"{n} records, {flipped} flipped [{time.time()-t0:.0f}s]", flush=True)
    pysam.index(outp)
    print(f"DONE {inp.split('/')[-1]}: {n} records, {flipped} flipped ({100*flipped/max(n, 1):.1f}%), "
          f"{notag} mapped without dl/dr [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    main()
