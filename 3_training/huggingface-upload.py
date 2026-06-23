from huggingface_hub import HfApi, create_repo

from methylseqnet_repro.paths import model_checkpoints

checkpoint_to_upload = model_checkpoints / "slurm32260895task0" / "checkpoints" / "best-checkpoint-train-factorized-pretrained-embeddings-v1.ckpt"

repo_id = "OberonDixon/methylseqnet"
create_repo(repo_id, repo_type="model", private=True, exist_ok=True)

HfApi().upload_file(
    path_or_fileobj=checkpoint_to_upload,
    path_in_repo="factorized-borzoi-rep0.ckpt",
    repo_id=repo_id,
)

# HfApi().create_tag(
#     repo_id="OberonDixon/methylseqnet",
#     repo_type="model",
#     tag="v1.0",
#     revision="ffc428c44553333c74bfcf73aa6a6c0ebdde7f42",
#     tag_message="Initial preprint posting",  # optional, makes it an annotated tag
# )