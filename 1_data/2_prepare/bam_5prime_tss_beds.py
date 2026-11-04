"""
TSS count beds from the 5' ends of the reads of an Iso-Seq BAM.

The TSS of a read is the first aligned base at its 5' end: the start of the alignment for a read on the + strand, the end of
the alignment for a read on the - strand (soft-clipped bases are not counted). A read's TSS depends on nothing but the read,
so the unphased and the haplotype beds share their TSS positions by construction. The reads are not clustered into
isoforms: `isoseq collapse` merges the reads that share their 3' end and splice junctions but differ at the 5' end into one
isoform with one TSS, which moves the TSS of most reads in deep libraries away from where they start (on HepG2 chr22,
only 19% of the TSS of the collapsed reads fall within 100 bp of an annotated TSS, against 88% of the reads' own 5' ends).

Every primary alignment of a mapped read counts once; unmapped, secondary and supplementary alignments are skipped. There is no
filtering by mapping quality or by how much of the read aligns. The strand is the strand flag of the alignment, so the BAM
has to hold reads in the orientation of their transcripts (see orient_bam_dldr.py for Kinnex BAMs that never went through
`isoseq refine`). The BAM does not need an index.

Outputs are bgzipped, tabix-indexed beds with one TSS per row: chrom, tss, tss+1, name ('.'), reads, strand, sorted by position.
  --unphased-output   all reads
  --hp1-output, --hp2-output   the reads with that HP tag (as written by whatshap haplotag); a read without an HP tag is in
                      the unphased bed only, so the HP1 and HP2 reads of a TSS add up to at most its unphased count
"""
import argparse
import collections
import sys

import pysam


def count_five_prime_ends(bam_path, threads=1):
    """
    Returns ({(chrom, tss, strand): [all reads, HP1 reads, HP2 reads]}, statistics) for the primary alignments of the BAM.
    """
    counts = collections.defaultdict(lambda: [0, 0, 0])
    stats = collections.Counter()
    with pysam.AlignmentFile(bam_path, 'rb', check_sq=False, threads=threads) as bam:
        for read in bam.fetch(until_eof=True):
            stats['alignments'] += 1
            if read.is_unmapped:
                stats['unmapped'] += 1
                continue
            if read.is_secondary or read.is_supplementary:
                stats['secondary or supplementary'] += 1
                continue
            if read.is_reverse:
                tss, strand = read.reference_end - 1, '-'
            else:
                tss, strand = read.reference_start, '+'
            entry = counts[(read.reference_name, tss, strand)]
            entry[0] += 1
            stats['reads'] += 1
            if read.has_tag('HP'):
                haplotype = read.get_tag('HP')
                if haplotype == 1 or haplotype == 2:
                    entry[haplotype] += 1
                    stats[f'HP{haplotype}'] += 1
    return counts, stats


def write_bedgz(counts, column, output_path):
    """
    Write the TSS with reads in `column` of the counts (0 all, 1 HP1, 2 HP2) as a bgzipped, tabix-indexed bed file
    (chrom, tss, tss+1, '.', reads, strand), sorted by position. Returns the number of TSS written.
    """
    if not output_path.endswith('.gz'):
        raise ValueError(f"the output path has to end in .gz: {output_path}")
    plain_path = output_path[:-len('.gz')]
    rows = sorted((chrom, tss, strand, c[column]) for (chrom, tss, strand), c in counts.items() if c[column] > 0)
    with open(plain_path, 'w') as f:
        for chrom, tss, strand, reads in rows:
            f.write(f"{chrom}\t{tss}\t{tss + 1}\t.\t{reads}\t{strand}\n")
    pysam.tabix_index(plain_path, preset='bed', force=True)   # bgzips the file to plain_path + '.gz' and indexes it
    return len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bam', help='BAM of aligned reads in the orientation of their transcripts')
    parser.add_argument('--unphased-output', required=True, help='Output bed.gz for all reads')
    parser.add_argument('--hp1-output', help='Output bed.gz for the reads with HP tag 1 (needs --hp2-output)')
    parser.add_argument('--hp2-output', help='Output bed.gz for the reads with HP tag 2 (needs --hp1-output)')
    parser.add_argument('--threads', type=int, default=1, help='BAM decompression threads')
    args = parser.parse_args()
    if (args.hp1_output is None) != (args.hp2_output is None):
        parser.error("--hp1-output and --hp2-output go together")

    counts, stats = count_five_prime_ends(args.bam, args.threads)
    print(f"{stats['alignments']} alignments: {stats['reads']} counted, {stats['unmapped']} unmapped, "
          f"{stats['secondary or supplementary']} secondary or supplementary; reads with an HP tag: "
          f"{stats['HP1']} HP1, {stats['HP2']} HP2, {stats['reads'] - stats['HP1'] - stats['HP2']} without", file=sys.stderr)
    if args.hp1_output and stats['HP1'] + stats['HP2'] == 0:
        raise SystemExit("The BAM has no reads with an HP tag, so there are no haplotype counts (was it haplotagged?).")

    outputs = [(0, args.unphased_output, 'all reads')]
    if args.hp1_output:
        outputs += [(1, args.hp1_output, 'HP1 reads'), (2, args.hp2_output, 'HP2 reads')]
    for column, output, description in outputs:
        num_written = write_bedgz(counts, column, output)
        print(f"Wrote {num_written} TSS with {description} to {output}", file=sys.stderr)


if __name__ == '__main__':
    main()
