"""
Extract TSS positions with counts from IsoSeq collapse output.

The TSS of an isoform is the start of its transcript in the collapse GFF (the end for the - strand), and its count is
the number of full-length reads associated with it in the abundance file. Also used as a library by
isoseq_hp_counts_to_bedgz.py (parse_gff, parse_abundance, write_bedgz).
"""
import argparse
import subprocess
import sys

def parse_gff(gff_path):
    tss_dict = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'transcript':
                continue
            chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            for attr in fields[8].split(';'):
                if 'transcript_id' in attr:
                    transcript_id = attr.split('"')[1]
                    break
            else:
                continue
            tss_pos = start - 1 if strand == '+' else end - 1
            tss_dict[transcript_id] = (chrom, tss_pos, strand)
    return tss_dict

def parse_abundance(abundance_path):
    count_dict = {}
    with open(abundance_path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('pbid'):
                continue
            fields = line.strip().split('\t')
            count_dict[fields[0]] = int(fields[1])
    return count_dict

def write_bedgz(entries, output_path):
    """
    Write (chrom, tss_pos, id, count, strand) entries as a bgzipped, tabix-indexed bed file with one TSS per row
    (chrom, tss_pos, tss_pos+1, id, count, strand), sorted by position.
    """
    entries = sorted(entries, key=lambda x: (x[0], x[1]))
    with open(output_path, 'wb') as output_file:
        bgzip_proc = subprocess.Popen(['bgzip', '-c'], stdin=subprocess.PIPE, stdout=output_file, text=True)
        for chrom, tss_pos, tid, count, strand in entries:
            bgzip_proc.stdin.write(f"{chrom}\t{tss_pos}\t{tss_pos+1}\t{tid}\t{count}\t{strand}\n")
        bgzip_proc.stdin.close()
        if bgzip_proc.wait() != 0:
            raise RuntimeError(f"bgzip failed writing {output_path}")
    subprocess.run(['tabix', '-f', '-p', 'bed', output_path], check=True)
    return len(entries)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('gff', help='GFF file from IsoSeq collapse')
    parser.add_argument('abundance', help='Abundance file from IsoSeq collapse')
    parser.add_argument('output', help='Output bed.gz path')
    args = parser.parse_args()

    tss_dict = parse_gff(args.gff)
    print(f"Found {len(tss_dict)} transcripts in GFF", file=sys.stderr)
    count_dict = parse_abundance(args.abundance)
    print(f"Found {len(count_dict)} transcripts in abundance file", file=sys.stderr)

    entries = [
        (chrom, tss_pos, tid, count_dict[tid], strand)
        for tid, (chrom, tss_pos, strand) in tss_dict.items()
        if tid in count_dict
    ]
    num_written = write_bedgz(entries, args.output)
    print(f"Wrote {num_written} TSS sites to {args.output}", file=sys.stderr)

if __name__ == '__main__':
    main()