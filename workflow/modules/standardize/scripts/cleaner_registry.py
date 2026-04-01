from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CleanerSpec:
    env: str  # "py" | "r"
    input_type: str
    output_type: str
    path_keys: tuple[str, ...]


CLEANER_REGISTRY: dict[str, CleanerSpec] = {
    "count_matrix_10x": CleanerSpec(
        env="py",
        input_type="raw_triplet",
        output_type="h5ad",
        path_keys=("input.matrix", "input.genes", "input.barcodes", "input.metadata"),
    ),
    "tenx_h5_multi": CleanerSpec(
        env="py",
        input_type="samples_h5_or_h5ad",
        output_type="h5ad",
        path_keys=("input.metadata", "input.samples[].h5", "input.samples[].h5ad"),
    ),
    "metadata_join": CleanerSpec(
        env="py",
        input_type="h5ad_plus_metadata",
        output_type="h5ad",
        path_keys=("input.metadata",),
    ),
    "seurat_rds": CleanerSpec(
        env="r",
        input_type="seurat_rds",
        output_type="h5ad",
        path_keys=("input.rds",),
    ),
    "csv_handler": CleanerSpec(
        env="py",
        input_type="raw_counts_csv_dir",
        output_type="h5ad",
        path_keys=("input.raw_dir",),
    ),
    "hra007511_count_matrix_runs": CleanerSpec(
        env="py",
        input_type="runs_count_matrix",
        output_type="h5ad",
        path_keys=("input.counts_root", "input.metadata", "input.runs_csv"),
    ),
}


def normalize_cleaner_name(value: str) -> str:
    cleaner = str(value).strip()
    if cleaner.endswith(".py"):
        cleaner = cleaner[:-3]
    return cleaner
