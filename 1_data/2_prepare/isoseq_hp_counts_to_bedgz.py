"""
Split the TSS counts of an IsoSeq collapse of ALL reads by haplotype.

Collapsing the HP1 and HP2 reads separately (and the unphased reads separately again) does not give comparable TSS
positions: the position of an isoform's TSS depends on which reads were collapsed together, so the same reads can land
in different 128 bp bins in the unphased and haplotype beds. Instead this script keeps the isoforms of the collapse of
all reads and attributes each of their reads to a haplotype, so the unphased and haplotype beds share their TSS positions
by construction:

  read -> isoform   from the collapse's read_stat.txt (one row per read: id, pbid)
  read -> haplotype from the HP tag that whatshap haplotag wrote to the BAM, or to several BAMs (reads without an HP tag count
                    for neither)
  isoform -> TSS    from the collapse's GFF, as in isoseq_collapsed_to_bedgz.py

The outputs are bgzipped, tabix-indexed beds in the same format as isoseq_collapsed_to_bedgz.py writes
(chrom, tss, tss+1, isoform, reads, strand), one for each haplotype, with a row for every isoform that has reads in that
haplotype. The reads of an isoform in the two haplotype beds add up to at most its count in the unphased bed; the rest are
reads without an HP tag.
"""
import argparse
import collections
import sys

import pysam

from isoseq_collapsed_to_bedgz import parse_gff, parse_abundance, write_bedgz


def read_haplotypes(bam_paths, threads=1):
    """Return ({read name: HP tag} for the reads that have one, number of alignments) over the BAMs."""
    haplotype_of_read = {}
    num_alignments = 0
    for bam_path in bam_paths:
        with pysam.AlignmentFile(bam_path, 'rb', check_sq=False, threads=threads) as bam:
            for read in bam.fetch(until_eof=True):
                num_alignments += 1
                if read.has_tag('HP'):
                    haplotype_of_read[read.query_name] = read.get_tag('HP')
    return haplotype_of_read, num_alignments


def count_reads_by_haplotype(read_stat_path, haplotype_of_read):
    """Return {isoform: [reads in HP1, reads in HP2, all reads]} from the collapse's read_stat.txt."""
    counts = collections.defaultdict(lambda: [0, 0, 0])
    with open(read_stat_path) as f:
        header = next(f).rstrip('\n').split('\t')
        if header[:2] != ['id', 'pbid']:
            raise ValueError(f"{read_stat_path}: expected a header starting with id, pbid, got {header}")
        for line in f:
            read, pbid = line.rstrip('\n').split('\t')[:2]
            isoform_counts = counts[pbid]
            isoform_counts[2] += 1
            haplotype = haplotype_of_read.get(read)
            if haplotype == 1 or haplotype == 2:
                isoform_counts[haplotype - 1] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gff', required=True, help='GFF from IsoSeq collapse of all reads')
    parser.add_argument('--abundance', required=True, help='abundance.txt from the same collapse')
    parser.add_argument('--read-stat', required=True, help='read_stat.txt from the same collapse')
    parser.add_argument('--bam', required=True, nargs='+',
                        help='BAM of the same reads with whatshap HP tags. Several can be given, e.g. BAMs that were split by '
                             'haplotype: it is enough that together they contain every read that has an HP tag')
    parser.add_argument('--hp1-output', required=True, help='Output bed.gz for the HP1 reads')
    parser.add_argument('--hp2-output', required=True, help='Output bed.gz for the HP2 reads')
    parser.add_argument('--threads', type=int, default=1, help='BAM decompression threads')
    args = parser.parse_args()

    tss_dict = parse_gff(args.gff)
    count_dict = parse_abundance(args.abundance)
    print(f"{len(tss_dict)} transcripts in the GFF, {len(count_dict)} in the abundance file", file=sys.stderr)

    haplotype_of_read, num_alignments = read_haplotypes(args.bam, args.threads)
    num_hp1 = sum(1 for hp in haplotype_of_read.values() if hp == 1)
    num_hp2 = sum(1 for hp in haplotype_of_read.values() if hp == 2)
    print(f"{num_alignments} alignments in the BAM; reads with an HP tag: {num_hp1} HP1, {num_hp2} HP2", file=sys.stderr)
    counts = count_reads_by_haplotype(args.read_stat, haplotype_of_read)

    # Consistency checks between the collapse outputs and the BAM
    num_reads = sum(c[2] for c in counts.values())
    attributed = [sum(c[0] for c in counts.values()), sum(c[1] for c in counts.values())]
    if num_hp1 + num_hp2 > 0 and sum(attributed) == 0:
        raise SystemExit("None of the BAM's HP-tagged reads are in read_stat.txt: the read names do not match "
                         "(is the BAM the input of this collapse?).")
    mismatched = sum(1 for pbid, c in counts.items() if count_dict.get(pbid) != c[2])
    print(f"{num_reads} reads in read_stat.txt over {len(counts)} isoforms: {attributed[0]} attributed to HP1 "
          f"({attributed[0] / max(num_hp1, 1):.1%} of the HP1-tagged reads in the BAM), {attributed[1]} to HP2 "
          f"({attributed[1] / max(num_hp2, 1):.1%} of the HP2-tagged reads), {num_reads - sum(attributed)} without an HP tag; "
          f"{mismatched} isoforms whose read_stat count differs from the abundance file", file=sys.stderr)
    if mismatched:
        print("WARNING: read_stat.txt and abundance.txt disagree for some isoforms; are they from the same collapse?", file=sys.stderr)
    missing_tss = [pbid for pbid in counts if pbid not in tss_dict]
    if missing_tss:
        print(f"WARNING: {len(missing_tss)} isoforms in read_stat.txt are not in the GFF and are skipped", file=sys.stderr)

    for haplotype, output in ((1, args.hp1_output), (2, args.hp2_output)):
        entries = sorted(
            (tss_dict[pbid][0], tss_dict[pbid][1], pbid, c[haplotype - 1], tss_dict[pbid][2])
            for pbid, c in counts.items() if c[haplotype - 1] > 0 and pbid in tss_dict
        )
        num_written = write_bedgz(entries, output)
        print(f"Wrote {num_written} TSS sites with HP{haplotype} reads to {output}", file=sys.stderr)


if __name__ == '__main__':
    main()
