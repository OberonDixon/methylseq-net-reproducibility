from matplotlib import pyplot as plt
import re
from torch.utils.data import DataLoader
from tqdm.auto import tqdm
from collections import defaultdict
import numpy as np
import seaborn as sns
from scipy import stats
import pickle
import os
import copy

from methylseqnet.dataset import BaseHDF5Dataset
from methylseqnet_repro.paths import preprocessed_datasets

load_peaks_from_fixed_cell_type = 'random_celltype_cpg.05-.15_10cts'

model_passes = [
    "slurm32260895task0/motif_insertions_synthetic_0.03",
    "slurm32260895task0/motif_insertions_synthetic_0.95",
    "slurm32260895task0/motif_insertions_truecondweight_0.0",
    "slurm32260895task1/motif_insertions",
    "slurm32260895task0/motif_insertions_synthetic_0.5",
    ]
task_id = 0
center_slice = slice(63,65)

files_written = 0

dataloaders_dict = {}
methyl_predictions_paths = {
    model_pass: preprocessed_datasets / f"motif_insertion/preprocessed/{load_peaks_from_fixed_cell_type}_2048/{model_pass}/predictions.h5"
    for model_pass in model_passes
}
for model_pass, methyl_predictions_path in methyl_predictions_paths.items():
    try:
        methyl_dataset = BaseHDF5Dataset(methyl_predictions_path, datasets=set(["predictions","cpg_density","imputed_conditioning_state_rep","conditional_seq_rep","unconditional_seq_rep"]), return_specifiers=True, batch_size=100)
        io_mappings_df = methyl_dataset.get_io_mappings_df()
        dataloaders_dict[model_pass] = DataLoader(methyl_dataset, batch_size=None, shuffle=False, num_workers=24)
    except Exception as e:
        print(f"task{task_id} -- could not load methyl dataset for {load_peaks_from_fixed_cell_type} at methylation level {model_pass}. {'File exists' if os.path.exists(methyl_predictions_path) else 'File does not exist'}.")
    task_id += 1

task_channels = io_mappings_df['absolute_channel'].tolist()
conditioning_states = io_mappings_df['absolute_cell_type'].tolist()
cell_type_names = [methylation_file.translate(str.maketrans('', '', '"\'[]')) for methylation_file in io_mappings_df['methylation_files'].tolist()]
data_types = io_mappings_df['data_type'].tolist()

SAMPLE_ATTRIBUTES_DICT = {
    # Per-sample, stored once (indexed by sample order)
    "cpg_density": None,
    "conditional_rep": None,
    "unconditional_rep": None,
    "seq_name": None,
    "file_name": None,
    # Per-channel (indexed by channel_idx, then by sample order)
    "prediction_profile_dict": dict(),
    "imputed_methyl_profile_dict": dict(),
    "cell_type_dict": dict(),
    "data_type_dict": dict(),
}

for model_pass, dataloader in tqdm(dataloaders_dict.items(), desc="Processing model passes", unit="model pass"):
    outputs_dict = defaultdict(dict)

    for batch in tqdm(dataloader, desc=f"Extracting all channels at methylation level {model_pass}", unit="batch"):
        # Vectorize: pull entire tensors to numpy once
        predictions = batch['predictions'][:, 0, task_channels, center_slice].mean(dim=-1).cpu().numpy()       # (B, n_channels, L)
        cpg_density = batch['cpg_density'].cpu().numpy()                                # (B, L)
        imputed = batch['imputed_conditioning_state_rep'][:, conditioning_states, 0, center_slice].mean(dim=-1).cpu().numpy()  # (B, n_channels, L)
        cond_rep = batch['conditional_seq_rep'][:, :, center_slice].mean(dim=2).cpu().numpy()          # (B, D)
        uncond_rep = batch['unconditional_seq_rep'][:, :, center_slice].mean(dim=2).cpu().numpy()      # (B, D)
        specifiers = batch['specifier']  # list of strings, len B

        for sample_index in range(predictions.shape[0]):
            seq_name = batch['specifier'][sample_index].split(">")[-1]
            file_name = batch['specifier'][sample_index].split(">")[0].split("/")[-1]
            peak_id = seq_name.split("_")[0].split(":")[0]
            if "endogenous" in file_name:
                if "shuffled" in file_name:
                    sample_type = "endogenous_shuffled"
                else:
                    sample_type = "endogenous"
            else:
                sample_type = "motif_inserted"
                insertion_replicate = seq_name.split("_")[1]
                inserted_motif = seq_name.split("_")[2].split(":")[0]

            if sample_type in ["endogenous", "endogenous_shuffled"]:
                if sample_type not in outputs_dict[peak_id]:
                    outputs_dict[peak_id][sample_type] = copy.deepcopy(SAMPLE_ATTRIBUTES_DICT)
                sample_subdictionary = outputs_dict[peak_id][sample_type]
            else:  # sample_type == "motif_inserted"
                if sample_type not in outputs_dict[peak_id]:
                    outputs_dict[peak_id][sample_type] = dict()
                if inserted_motif not in outputs_dict[peak_id][sample_type]:
                    outputs_dict[peak_id][sample_type][inserted_motif] = dict()
                if insertion_replicate not in outputs_dict[peak_id][sample_type][inserted_motif]:
                    outputs_dict[peak_id][sample_type][inserted_motif][insertion_replicate] = copy.deepcopy(SAMPLE_ATTRIBUTES_DICT)
                sample_subdictionary = outputs_dict[peak_id][sample_type][inserted_motif][insertion_replicate]
            
            sample_subdictionary["cpg_density"] = cpg_density[sample_index]
            sample_subdictionary["conditional_rep"] = cond_rep[sample_index]
            sample_subdictionary["unconditional_rep"] = uncond_rep[sample_index]
            sample_subdictionary["seq_name"] = seq_name
            sample_subdictionary["file_name"] = file_name

            for ch_pos, (state_idx, channel_idx, cell_type_name, data_type) in enumerate(zip(conditioning_states, task_channels, cell_type_names, data_types)):
                sample_subdictionary["prediction_profile_dict"][channel_idx] = (predictions[sample_index, ch_pos])
                sample_subdictionary["imputed_methyl_profile_dict"][channel_idx] = (imputed[sample_index, ch_pos])
                sample_subdictionary["cell_type_dict"][channel_idx] = (cell_type_name)
                sample_subdictionary["data_type_dict"][channel_idx] = (data_type)

    safe_model_pass = model_pass.replace("/", "__")
    with open(preprocessed_datasets / f'pickles/{safe_model_pass}__outputs_dict_cpg,imputed,reps___{load_peaks_from_fixed_cell_type}.pkl', 'wb') as f:
        pickle.dump(outputs_dict, f)