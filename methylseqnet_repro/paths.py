import tomllib
from pathlib import Path
from importlib.resources import files

_CONFIG_DIR = files("methylseqnet_repro")
_CONFIG_PATH = _CONFIG_DIR / "paths.toml"
_EXAMPLE_PATH = _CONFIG_DIR / "paths.toml.example"

if not _CONFIG_PATH.exists():
    raise FileNotFoundError(
        f"Local paths config not found at {_CONFIG_PATH}\n"
        f"Copy the example and edit it for your machine:\n"
        f"  cp {_EXAMPLE_PATH} {_CONFIG_PATH}"
    )

with open(_CONFIG_PATH, "rb") as f:
    _config = tomllib.load(f)

configs = Path(_config["configs"])

pacbio_5mC_tracks = Path(_config["pacbio_5mC_tracks"])
fiberseq_tracks = Path(_config["fiberseq_tracks"])
rna_tracks = Path(_config["rna_tracks"])

methylation_atlas = Path(_config["methylation_atlas"])
atac_atlas = Path(_config["atac_atlas"])
cage_atlas = Path(_config["cage_atlas"])

genomes = Path(_config["genomes"])
haplotypes = Path(_config["haplotypes"])

preprocessed_datasets = Path(_config["preprocessed_datasets"])
model_checkpoints = Path(_config["model_checkpoints"])

analysis_pickles = Path(_config["analysis_pickles"])
annotations = Path(_config["annotations"])

def print_path():
    import sys
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="Path name to retrieve")
    args = parser.parse_args()
    paths = {
        "pacbio_5mC_tracks": pacbio_5mC_tracks,
        "fiberseq_tracks": fiberseq_tracks,
        "rna_tracks": rna_tracks,
        "methylation_atlas": methylation_atlas,
        "atac_atlas": atac_atlas,
        "cage_atlas": cage_atlas,
        "genomes": genomes,
        "haplotypes": haplotypes,
        "preprocessed_datasets": preprocessed_datasets,
        "model_checkpoints": model_checkpoints,
        "analysis_pickles": analysis_pickles,
        "annotations": annotations,
        "configs": configs,
    }
    print(paths[args.name])
