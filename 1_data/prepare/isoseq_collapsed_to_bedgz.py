"""
Extract TSS positions with counts from IsoSeq collapse output
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
    entries.sort(key=lambda x: (x[0], x[1]))

    bgzip_proc = subprocess.Popen(
        ['bgzip', '-c'], stdin=subprocess.PIPE,
        stdout=open(args.output, 'wb'), text=True
    )
    for chrom, tss_pos, tid, count, strand in entries:
        bgzip_proc.stdin.write(f"{chrom}\t{tss_pos}\t{tss_pos+1}\t{tid}\t{count}\t{strand}\n")
    bgzip_proc.stdin.close()
    bgzip_proc.wait()

    print(f"Wrote {len(entries)} TSS sites to {args.output}", file=sys.stderr)
    subprocess.run(['tabix', '-p', 'bed', args.output], check=True)

if __name__ == '__main__':
    main()