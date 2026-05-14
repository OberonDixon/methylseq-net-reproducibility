import pickle
import pysam
import pyBigWig
from tqdm.auto import tqdm
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from pathlib import Path
import pandas as pd

from methylseqnet.predict import Predictor
from methylseqnet.callbacks import HaplotypedPredLogger
from methylseqnet.readers import load_track
from methylseqnet_repro.paths import model_checkpoints, pacbio_5mC_tracks, fiberseq_tracks, rna_tracks, haplotypes, genomes, analysis_pickles, annotations, atac_atlas, cage_atlas, configs

atac_file = str(atac_atlas / "bigwigs/Memory%20B.bw")
cage_file = str(cage_atlas / "hg38_v9/CD19+ B Cells, donor1.CNhs12343.11544-120B5.hg38.nobarcode.ctss.sorted.bed.gz")
hp1_genome_fasta = str(haplotypes / "NA12878.hg38.hp1.fa")
hp2_genome_fasta = str(haplotypes / "NA12878.hg38.hp2.fa")
hp1_cpg = str(pacbio_5mC_tracks / "GM12878_WGS-pb-5mC.hap1.bw")
hp2_cpg = str(pacbio_5mC_tracks / "GM12878_WGS-pb-5mC.hap2.bw")
hp1_fiber = str(fiberseq_tracks / "GM12878_trackHub/bw/hap1.acc.bw")
hp2_fiber = str(fiberseq_tracks / "GM12878_trackHub/bw/hap2.acc.bw")
rna_file = str(rna_tracks / "GM12878.kinnex.no5exon.tss.counts.bed.gz")
hp1_rna = str(rna_tracks / "GM12878.kinnex.HP2.no5exon.tss.counts.bed.gz")
hp2_rna = str(rna_tracks / "GM12878.kinnex.HP1.no5exon.tss.counts.bed.gz")
rna_file_hp1 = hp1_rna
rna_file_hp2 = hp2_rna

def load_bed_intervals(bed_file, descriptions=None):
    """Load a BED file into a dict of {chrom: [(start, end), ...]}."""
    intervals = {}
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if descriptions is not None and parts[3].upper() not in [desc.upper() for desc in descriptions]:
                continue
            intervals.setdefault(parts[0], []).append((int(parts[1]), int(parts[2])))
    return intervals

def point_in_intervals(intervals, chrom, pos):
    """Check if a position falls within any interval for the given chrom."""
    return any(s <= pos < e for s, e in intervals.get(chrom, []))

def fetch_tss_for_chrom(
    tabix_path,
    fetch_params,
    subset,
    intervals=None,
    ):
    gtf = pysam.TabixFile(tabix_path)
    tss_list = []
    for record in gtf.fetch(*fetch_params):
        fields = record.split("\t")
        if fields[2] == "transcript":
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]

            # Parse attributes
            attr_dict = {}
            for attr in attributes.strip().split(';'):
                attr = attr.strip()
                if attr:
                    parts = attr.split(' ', 1)
                    if len(parts) == 2:
                        key = parts[0]
                        val = parts[1].strip('"')
                        attr_dict[key] = val
            
            # TSS is start for + strand, end for - strand
            tss = start if strand == '+' else end

            tss_dict = {
                'chrom': chrom,
                'tss': tss,
                'strand': strand,
                'gene_id': attr_dict.get('gene_id', ''),
                'gene_name': attr_dict.get('gene_name', '').upper(),
                'transcript_id': attr_dict.get('transcript_id', ''),
                'transcript_type': attr_dict.get('transcript_type', ''),
                'canonical': 'Ensembl_canonical' in attributes
            }

            if (
                    (len(subset)==0 or (tss_dict['gene_name'] in subset))
                    and tss_dict['canonical']
                    # and tss_dict['transcript_type'] == 'protein_coding'
                    and (intervals is None or point_in_intervals(intervals, chrom, tss))
                ):
                tss_list.append(tss_dict)
    return tss_list

imprinted_df = pd.read_csv("imprinted_genes.csv", encoding='mac_roman')
imprinted_df['All_Names'] = imprinted_df.apply(
    lambda row: [row['Gene'].strip().upper()] + 
                [alias.strip().upper() for alias in row['Aliases'].split(', ')] 
                if pd.notna(row['Aliases']) and row['Aliases'] else [row['Gene'].strip().upper()], 
    axis=1
)
imprinted_subset = set([name for names in imprinted_df['All_Names'] for name in names])
normalize_counts_per = None

for folds in "test",: #"allfolds":
    training_splits_bed = configs / "training-splits" / "borzoi-human" / "sequences_human.bed"
    if folds == "heldoutfolds":
        intervals = load_bed_intervals(training_splits_bed, descriptions=["fold3","fold4"])
    elif folds == "test":
        intervals = load_bed_intervals(training_splits_bed, descriptions=["fold3"])
    else:
        intervals = None
    for channels_selection in ['rna-fiber']: #,'all','blood-b'
        for checkpoint_id in [
            "slurm32260895task0",
            "slurm32260895task1",
        ]: 
            best_ckpt = max(
                (model_checkpoints / f"{checkpoint_id}/checkpoints/").glob('best*.ckpt'),
                key=lambda p: p.stat().st_mtime
            )    
            predictor = Predictor(best_ckpt, supplemental_outputs={"true_conditioning_state_rep","conditional_seq_rep","unconditional_seq_rep"})

            for label, fetch_params, subset in zip(
                ['chrX','chr13','imprinted'],
                [("chrX",),("chr13",),tuple()],
                [(),(),imprinted_subset],
            ):
                tss_list = fetch_tss_for_chrom(
                    str(annotations / "gencode.v49.basic.annotation.sorted.gtf.gz"),
                    # "/clusterfs/nilah/oberon/datasets/gencode/gencode.v49.basic.annotation.sorted.gtf.gz",
                    fetch_params=fetch_params,
                    subset=subset,
                    intervals=intervals,
                )
                io_mappings_df = predictor.model.get_io_mappings_df()
                if channels_selection == 'all':
                    atac_channels = io_mappings_df[io_mappings_df["data_type"] == "ATAC-seq"]["absolute_channel"].values
                    cage_channels = io_mappings_df[io_mappings_df["data_type"] == "CAGE-seq"]["absolute_channel"].values
                elif channels_selection == 'blood-b':
                    cage_channels = io_mappings_df[io_mappings_df["methylation_files"].str.contains("Blood-B") & (io_mappings_df["data_type"] == "CAGE-seq")]["absolute_channel"].values
                    atac_channels = io_mappings_df[io_mappings_df["methylation_files"].str.contains("Blood-B") & (io_mappings_df["data_type"] == "ATAC-seq")]["absolute_channel"].values
                elif channels_selection == 'rna-fiber':
                    cage_channels = io_mappings_df[io_mappings_df["data_type"] == "RNA-seq"]["absolute_channel"].values
                    atac_channels = io_mappings_df[io_mappings_df["data_type"] == "Fiber-seq"]["absolute_channel"].values
                else:
                    raise ValueError(f"Unknown channels_selection: {channels_selection}")
                
                receptive_window_size = 524288

                atac_bw = pyBigWig.open(atac_file)
                # rna_bam = pysam.AlignmentFile(rna_file, "rb")
                # rna_bam_hp1 = pysam.AlignmentFile(rna_file_hp1, "rb")
                # rna_bam_hp2 = pysam.AlignmentFile(rna_file_hp2, "rb")
                cage_bed = pysam.TabixFile(cage_file)

                chrom_sizes = atac_bw.chroms()
                acc_capture_window = 2
                rna_capture_window = 8
                checkpoint_file = f"pickles/checkpoint_{checkpoint_id}_{folds}_{label}_capture{acc_capture_window},{rna_capture_window}_{channels_selection}_snp_sequence.pkl"

                save_every = 20

                atac_cts_list = []
                rna_cts_list = []
                rna_cts_list_hp1 = []
                rna_cts_list_hp2 = []
                cage_cts_list = []

                hp1_rna_preds = []
                hp2_rna_preds = []
                hp1_acc_preds = []
                hp2_acc_preds = []
                hp1_methylations = []
                hp2_methylations = []
                hp1_rna_methylations = []
                hp2_rna_methylations = []
                hp1_rna_targets = []
                hp2_rna_targets = []
                hp1_acc_targets = []
                hp2_acc_targets = []

                hp1_acc_conditional_seq_reps = []
                hp2_acc_conditional_seq_reps = []
                hp1_rna_conditional_seq_reps = []
                hp2_rna_conditional_seq_reps = []
                hp1_acc_unconditional_seq_reps = []
                hp2_acc_unconditional_seq_reps = []
                hp1_rna_unconditional_seq_reps = []
                hp2_rna_unconditional_seq_reps = []
                

                for i, tss_dict in enumerate(tqdm(tss_list,desc=f"Running TSS on {label} through {checkpoint_id}")):
                    start = tss_dict['tss'] - receptive_window_size//2
                    end = tss_dict['tss'] + receptive_window_size//2
                    chrom = tss_dict['chrom']
                    gene_name = tss_dict['gene_name']
                    chrom_length = chrom_sizes.get(chrom, 0)
                    if start<0 or end>chrom_length:
                        continue
                    hp1_predictions_dict = predictor.predict_locus(
                        chromosome=chrom,
                        start=start,
                        end=end,
                        sequence_path=hp1_genome_fasta,
                        methylation_paths=hp1_cpg,
                        methylation_load_kwargs={"extend_cpg_sites": True, "cpg_values_rescale": 0.01},
                    )
                    hp1_pred_atac = hp1_predictions_dict["predictions"][...,atac_channels,:].mean(dim=-2).squeeze()
                    hp1_pred_cage = hp1_predictions_dict["predictions"][...,cage_channels,:].mean(dim=-2).squeeze()
                    hp1_methylation = hp1_predictions_dict["true_conditioning_state_rep"].squeeze()
                    hp1_conditional_seq_rep = hp1_predictions_dict["conditional_seq_rep"].squeeze()
                    hp1_unconditional_seq_rep = hp1_predictions_dict["unconditional_seq_rep"].squeeze()
                    hp2_predictions_dict = predictor.predict_locus(
                        chromosome=chrom,
                        start=start,
                        end=end,
                        sequence_path=hp2_genome_fasta,
                        methylation_paths=hp2_cpg,
                        methylation_load_kwargs={"extend_cpg_sites": True, "cpg_values_rescale": 0.01},
                    )
                    hp2_pred_atac = hp2_predictions_dict["predictions"][...,atac_channels,:].mean(dim=-2).squeeze()
                    hp2_pred_cage = hp2_predictions_dict["predictions"][...,cage_channels,:].mean(dim=-2).squeeze()
                    hp2_methylation = hp2_predictions_dict["true_conditioning_state_rep"].squeeze()
                    hp2_conditional_seq_rep = hp2_predictions_dict["conditional_seq_rep"].squeeze()
                    hp2_unconditional_seq_rep = hp2_predictions_dict["unconditional_seq_rep"].squeeze()
                    hp1_target_atac = predictor.load_targets(
                        chromosome=chrom,
                        start=start,
                        end=end,
                        target_paths=[hp1_fiber],
                        normalize_counts_per=normalize_counts_per,
                    ).squeeze()
                    hp2_target_atac = predictor.load_targets(
                        chromosome=chrom,
                        start=start,
                        end=end,
                        target_paths=[hp2_fiber],
                        normalize_counts_per=normalize_counts_per,
                    ).squeeze()
                    hp1_target_rna = predictor.load_targets(
                        chromosome=chrom,
                        start=start,
                        end=end,
                        target_paths=[hp1_rna],
                        normalize_counts_per=normalize_counts_per,
                    ).squeeze()
                    hp2_target_rna = predictor.load_targets(
                        chromosome=chrom,
                        start=start,
                        end=end,
                        target_paths=[hp2_rna],
                        normalize_counts_per=normalize_counts_per,
                    ).squeeze()
                    atac_cts = predictor.load_targets(
                        chromosome=chrom,
                        start=start,
                        end=end,
                        target_paths=[atac_file],
                        normalize_counts_per=normalize_counts_per,
                    ).squeeze()
                    rna_cts = predictor.load_targets(
                        chromosome=chrom,
                        start=start,
                        end=end,
                        target_paths=[rna_file],
                        normalize_counts_per=normalize_counts_per,
                    ).squeeze()
                    rna_cts_hp1 = predictor.load_targets(
                        chromosome=chrom,
                        start=start,
                        end=end,
                        target_paths=[rna_file_hp1],
                        normalize_counts_per=normalize_counts_per,
                    ).squeeze()
                    rna_cts_hp2 = predictor.load_targets(
                        chromosome=chrom,
                        start=start,
                        end=end,
                        target_paths=[rna_file_hp2],
                        normalize_counts_per=normalize_counts_per,
                    ).squeeze()
                    cage_cts = predictor.load_targets(
                        chromosome=chrom,
                        start=start,
                        end=end,
                        target_paths=[cage_file],
                        normalize_counts_per=normalize_counts_per,
                    ).squeeze()
                    pred_len = (hp1_pred_cage.shape[0])
                    acc_trimmed_slice = slice(pred_len//2 - acc_capture_window, pred_len//2 + acc_capture_window)
                    rna_trimmed_slice = slice(pred_len//2 - rna_capture_window, pred_len//2 + rna_capture_window)
                    hp1_rna_preds.append(hp1_pred_cage[rna_trimmed_slice].mean())
                    hp2_rna_preds.append(hp2_pred_cage[rna_trimmed_slice].mean())
                    hp1_acc_preds.append(hp1_pred_atac[acc_trimmed_slice].mean())
                    hp2_acc_preds.append(hp2_pred_atac[acc_trimmed_slice].mean())
                    hp1_rna_targets.append(hp1_target_rna[rna_trimmed_slice].mean())
                    hp2_rna_targets.append(hp2_target_rna[rna_trimmed_slice].mean())
                    hp1_acc_targets.append(hp1_target_atac[acc_trimmed_slice].mean())
                    hp2_acc_targets.append(hp2_target_atac[acc_trimmed_slice].mean())
                    hp1_methylations.append(hp1_methylation[acc_trimmed_slice].mean())
                    hp2_methylations.append(hp2_methylation[acc_trimmed_slice].mean())
                    hp1_rna_methylations.append(hp1_methylation[rna_trimmed_slice].mean())
                    hp2_rna_methylations.append(hp2_methylation[rna_trimmed_slice].mean())
                    atac_cts_list.append(atac_cts[acc_trimmed_slice].mean())
                    rna_cts_list.append(rna_cts[rna_trimmed_slice].mean())
                    rna_cts_list_hp1.append(rna_cts_hp1[rna_trimmed_slice].mean())
                    rna_cts_list_hp2.append(rna_cts_hp2[rna_trimmed_slice].mean())
                    cage_cts_list.append(cage_cts[rna_trimmed_slice].mean())
                    hp1_acc_conditional_seq_reps.append(hp1_conditional_seq_rep[...,acc_trimmed_slice].mean())
                    hp2_acc_conditional_seq_reps.append(hp2_conditional_seq_rep[...,acc_trimmed_slice].mean())
                    hp1_rna_conditional_seq_reps.append(hp1_conditional_seq_rep[...,rna_trimmed_slice].mean())
                    hp2_rna_conditional_seq_reps.append(hp2_conditional_seq_rep[...,rna_trimmed_slice].mean())
                    hp1_acc_unconditional_seq_reps.append(hp1_unconditional_seq_rep[...,acc_trimmed_slice].mean())
                    hp2_acc_unconditional_seq_reps.append(hp2_unconditional_seq_rep[...,acc_trimmed_slice].mean())
                    hp1_rna_unconditional_seq_reps.append(hp1_unconditional_seq_rep[...,rna_trimmed_slice].mean())
                    hp2_rna_unconditional_seq_reps.append(hp2_unconditional_seq_rep[...,rna_trimmed_slice].mean())
                    # Save checkpoint every N iterations
                    if (i + 1) % save_every == 0:
                        checkpoint = {
                            'atac_cts_list': atac_cts_list,
                            'rna_cts_list': rna_cts_list,
                            'rna_cts_list_hp1': rna_cts_list_hp1,
                            'rna_cts_list_hp2': rna_cts_list_hp2,
                            'cage_cts_list': cage_cts_list,
                            'hp1_preds': hp1_rna_preds,
                            'hp2_preds': hp2_rna_preds,
                            'hp1_acc_preds': hp1_acc_preds,
                            'hp2_acc_preds': hp2_acc_preds,
                            'hp1_methylations': hp1_methylations,
                            'hp2_methylations': hp2_methylations,
                            'hp1_rna_methylations': hp1_rna_methylations,
                            'hp2_rna_methylations': hp2_rna_methylations,
                            'hp1_rna_targets': hp1_rna_targets,
                            'hp2_rna_targets': hp2_rna_targets,
                            'hp1_acc_targets': hp1_acc_targets,
                            'hp2_acc_targets': hp2_acc_targets,
                            'hp1_acc_conditional_seq_reps': hp1_acc_conditional_seq_reps,
                            'hp2_acc_conditional_seq_reps': hp2_acc_conditional_seq_reps,
                            'hp1_rna_conditional_seq_reps': hp1_rna_conditional_seq_reps,
                            'hp2_rna_conditional_seq_reps': hp2_rna_conditional_seq_reps,
                            'hp1_acc_unconditional_seq_reps': hp1_acc_unconditional_seq_reps,
                            'hp2_acc_unconditional_seq_reps': hp2_acc_unconditional_seq_reps,
                            'hp1_rna_unconditional_seq_reps': hp1_rna_unconditional_seq_reps,
                            'hp2_rna_unconditional_seq_reps': hp2_rna_unconditional_seq_reps,
                            'gene_name': gene_name,
                            'iteration': i + 1
                        }
                        # for key, value in checkpoint.items():
                        #     print(f"{key}: {len(value) if isinstance(value, list) else value}, dtype: {type(value[0]) if isinstance(value, list) and value else type(value)}")
                        with open(checkpoint_file, 'wb') as f:
                            pickle.dump(checkpoint, f)

                # Save final results
                final_results = {
                    'atac_cts_list': atac_cts_list,
                    'rna_cts_list': rna_cts_list,
                    'rna_cts_list_hp1': rna_cts_list_hp1,
                    'rna_cts_list_hp2': rna_cts_list_hp2,
                    'cage_cts_list': cage_cts_list,
                    'hp1_preds': hp1_rna_preds,
                    'hp2_preds': hp2_rna_preds,
                    'hp1_acc_preds': hp1_acc_preds,
                    'hp2_acc_preds': hp2_acc_preds,
                    'hp1_methylations': hp1_methylations,
                    'hp2_methylations': hp2_methylations,
                    'hp1_rna_methylations': hp1_rna_methylations,
                    'hp2_rna_methylations': hp2_rna_methylations,
                    'hp1_rna_targets': hp1_rna_targets,
                    'hp2_rna_targets': hp2_rna_targets,
                    'hp1_acc_targets': hp1_acc_targets,
                    'hp2_acc_targets': hp2_acc_targets,
                    'hp1_acc_conditional_seq_reps': hp1_acc_conditional_seq_reps,
                    'hp2_acc_conditional_seq_reps': hp2_acc_conditional_seq_reps,
                    'hp1_rna_conditional_seq_reps': hp1_rna_conditional_seq_reps,
                    'hp2_rna_conditional_seq_reps': hp2_rna_conditional_seq_reps,
                    'hp1_acc_unconditional_seq_reps': hp1_acc_unconditional_seq_reps,
                    'hp2_acc_unconditional_seq_reps': hp2_acc_unconditional_seq_reps,
                    'hp1_rna_unconditional_seq_reps': hp1_rna_unconditional_seq_reps,
                    'hp2_rna_unconditional_seq_reps': hp2_rna_unconditional_seq_reps,
                    'gene_name': gene_name,
                }
                with open(checkpoint_file, 'wb') as f:
                    pickle.dump(final_results, f)
