"""
Compute, once and in parallel over the cell lines, the depth-corrected library size of the Iso-seq RNA track that
preprocessing normalizes (see MultiBedGzLabelHandler in methylseqnet.builders), so that the preprocessing jobs, which would
all otherwise compute it themselves, find it cached next to the bed file.

For each cell line of a cell lines TSV it takes the unphased bed and the depth bigwigs of all reads that the preprocessing
corrects it with: unphased_rna_file with unphased_rna_depth_files for a phased TSV, or rna_files with rna_depth_files for an
unphased TSV (the two give the same result for a cell line, whose unphased bed and all-reads depth files are the same).

The clamp and exponent must be those of the preprocessing config (depth_ratio_min, depth_ratio_max, depth_exponent in its 'rna'
kwargs; the defaults here are the ones of the _rnadepth configs), otherwise the preprocessing finds a different cache key and
computes the library size again. The same goes for the files themselves: the cache is keyed by their paths, sizes and times.

Usage: python precompute_rna_cn_totals.py configs/longread_cell_lines_phased.tsv [--processes 4] [--cache-dir DIR] [--lines HepG2 K562]
"""
import argparse
import json
import multiprocessing
import time

from methylseqnet.builders import MultiBedGzLabelHandler, PhasedFiberRNA


def bed_and_depth_files(cell_line):
    """The unphased bed and the all-reads depth files that preprocessing normalizes the RNA with, or None."""
    if cell_line['unphased_rna_file'] is not None and cell_line['unphased_rna_depth_files'] is not None:
        return cell_line['unphased_rna_file'], cell_line['unphased_rna_depth_files']
    if cell_line['rna_files'] is not None and cell_line['rna_depth_files'] is not None and len(cell_line['rna_files']) == 1:
        return cell_line['rna_files'][0], cell_line['rna_depth_files'][0]
    return None


def compute(job):
    name, bed, depth_files, kwargs = job
    start = time.time()
    handler = MultiBedGzLabelHandler([bed], normalize_counts_per=1e9, depth_bigwig_files=depth_files, **kwargs)
    with open(handler._depth_corrected_total_cache_file(bed)) as f:
        totals = json.load(f)
    return name, totals['raw_total'], totals['corrected_total'], totals['uncorrected_counts'], time.time() - start


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cell_lines_tsv')
    parser.add_argument('--processes', type=int, default=4)
    parser.add_argument('--lines', nargs='+', help='only these state_names')
    parser.add_argument('--cache-dir', help='cache here instead of next to the beds')
    parser.add_argument('--label-bin-size', type=int, default=128)
    parser.add_argument('--depth-ratio-min', type=float, default=0.25)
    parser.add_argument('--depth-ratio-max', type=float, default=32.0)
    parser.add_argument('--depth-exponent', type=float, default=1.0)
    args = parser.parse_args()

    with open(args.cell_lines_tsv) as f:
        rows = [line for line in f if line.strip() and not line.startswith('#')]
    num_phases = rows[1].split('\t')[1].count(';') + 1   # the methylation_files column has one file per phase
    cell_lines = PhasedFiberRNA._parse_cell_lines_tsv(args.cell_lines_tsv, num_phases)
    kwargs = dict(label_bin_size=args.label_bin_size, depth_ratio_min=args.depth_ratio_min, depth_ratio_max=args.depth_ratio_max,
                  depth_exponent=args.depth_exponent, depth_total_cache_dir=args.cache_dir)
    jobs = []
    for cell_line in cell_lines:
        if args.lines and cell_line['state_name'] not in args.lines:
            continue
        pair = bed_and_depth_files(cell_line)
        if pair is None:
            print(f"{cell_line['state_name']}: no unphased bed with depth files, skipping")
            continue
        jobs.append((cell_line['state_name'], pair[0], pair[1], kwargs))
    with multiprocessing.Pool(min(args.processes, max(len(jobs), 1))) as pool:
        results = pool.map(compute, jobs, chunksize=1)
    print(f"\n{'line':10} {'raw total':>12} {'corrected total':>16} {'corrected/raw':>14} {'uncorrected':>12} {'time':>7}")
    for name, raw, corrected, uncorrected, seconds in results:
        print(f"{name:10} {raw:12.0f} {corrected:16.0f} {corrected / raw:14.3f} {uncorrected / raw:11.2%} {seconds:6.0f}s")


if __name__ == '__main__':
    main()
