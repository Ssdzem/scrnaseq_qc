import os
from pathlib import Path


QC_CONFIG_PATH = str(Path("config/qc.yaml"))


def _studies_to_process():
    studies = config.get("inputs", {}).get("studies", [])
    if studies:
        return list(studies)

    root = Path("data_backup")
    if not root.exists():
        return []
    return sorted([p.name for p in root.iterdir() if p.is_dir()])


def _detect_inputs_for_study(study: str) -> dict:
    """Auto-detect study input type(s) from data_backup/{study}/*.

    Priority:
      1) If a .rds exists -> treat as RDS study (single sample).
      2) If a .mtx/.mtx.gz exists AND genes+barcodes csv exist -> treat as MTX study.
      3) Any .h5 files -> treat as one-sample-per-h5.
    """

    study_dir = Path("data_backup") / study
    detected: dict = {"rds": None, "mtx": None, "h5": []}
    if not study_dir.exists():
        return detected

    rds_files = sorted(study_dir.glob("*.rds"))
    if rds_files:
        detected["rds"] = str(rds_files[0])

    mtx_candidates = sorted(list(study_dir.rglob("*.mtx")) + list(study_dir.rglob("*.mtx.gz")))
    if mtx_candidates:
        genes_candidates = sorted(list(study_dir.rglob("*genes*.csv")) + list(study_dir.rglob("*features*.tsv")) + list(study_dir.rglob("*features*.csv")))
        barcodes_candidates = sorted(list(study_dir.rglob("*cell*.csv")) + list(study_dir.rglob("*barcode*.tsv")) + list(study_dir.rglob("*barcode*.csv")))
        if genes_candidates and barcodes_candidates:
            detected["mtx"] = {
                "matrix": str(mtx_candidates[0]),
                "genes": str(genes_candidates[0]),
                "barcodes": str(barcodes_candidates[0]),
            }

    detected["h5"] = [str(p) for p in sorted(study_dir.glob("*.h5"))]
    return detected


def _samples_for_study(study: str):
    d = _detect_inputs_for_study(study)
    samples: list[str] = []
    if d["rds"]:
        samples.append(f"{study}_rds")
    if d["mtx"]:
        samples.append(f"{study}_mtx")
    for h5 in d["h5"]:
        samples.append(Path(h5).stem)
    return sorted(set(samples))


STUDIES = _studies_to_process()


rule qc_all:
    input:
        expand(
            "data_backup/{study}/{sample}_qc.h5",
            zip,
            study=[s for st in STUDIES for s in [st] * len(_samples_for_study(st))],
            sample=[sa for st in STUDIES for sa in _samples_for_study(st)],
        )


rule rds_to_h5ad:
    input:
        rds=lambda wc: _detect_inputs_for_study(wc.study)["rds"]
    output:
        h5ad=temp("data_backup/{study}/{study}_rds.h5ad")
    params:
        study=lambda wc: wc.study
    log:
        "out/logs/qc_preproc/rds_to_h5ad/{study}.log"
    conda:
        "../env/qc_preproc_r.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        Rscript workflow/modules/qc_preproc/scripts/rds_to_h5ad.R \
          --input-rds {input.rds} \
          --output-h5ad {output.h5ad} \
          --log {log}
        """


rule qc_scanpy:
    input:
        raw=lambda wc: (
            _detect_inputs_for_study(wc.study)["rds"] if wc.sample == f"{wc.study}_rds" else
            _detect_inputs_for_study(wc.study)["mtx"]["matrix"] if wc.sample == f"{wc.study}_mtx" else
            str(Path("data_backup") / wc.study / f"{wc.sample}.h5")
        ),
        qc_config=QC_CONFIG_PATH,
        genes=lambda wc: [_detect_inputs_for_study(wc.study)["mtx"]["genes"]] if wc.sample == f"{wc.study}_mtx" else [],
        barcodes=lambda wc: [_detect_inputs_for_study(wc.study)["mtx"]["barcodes"]] if wc.sample == f"{wc.study}_mtx" else [],
        h5ad=lambda wc: [f"data_backup/{wc.study}/{wc.study}_rds.h5ad"] if wc.sample == f"{wc.study}_rds" else [],
    output:
        qc_h5="data_backup/{study}/{sample}_qc.h5",
    params:
        study=lambda wc: wc.study,
        sample=lambda wc: wc.sample,
    log:
        "out/logs/qc_preproc/qc_scanpy/{study}_{sample}.log"
    conda:
        "../env/qc_preproc_py.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})
        python workflow/modules/qc_preproc/scripts/qc_scanpy.py \
          --study {params.study} \
          --sample {params.sample} \
          --input {input.raw} \
          --genes '{input.genes}' \
          --barcodes '{input.barcodes}' \
          --h5ad '{input.h5ad}' \
          --qc-config {input.qc_config} \
          --output {output.qc_h5} \
          --log {log}
        """
