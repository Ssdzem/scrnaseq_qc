import os


METADATA_CSV = os.path.expanduser(
    "~/phd_coding/a1_p1_qc/config/liver_scRNAseq_metadata.merged.csv"
)
OUT_ROOT = os.path.expanduser("~/phd_coding/a1_p1_qc/data")


def list_studies(metadata_csv: str) -> list[str]:
    import json
    import subprocess

    cmd = [
        "python",
        "scripts/fetch_data.py",
        "--metadata",
        metadata_csv,
        "--out-root",
        OUT_ROOT,
        "list-studies-json",
    ]
    out = subprocess.check_output(cmd, text=True)
    return json.loads(out)


STUDIES = list_studies(METADATA_CSV)


rule all:
    input:
        expand(f"{OUT_ROOT}/{{study}}/.fetch_done", study=STUDIES)


rule fetch_study:
    input:
        metadata=METADATA_CSV,
        script="scripts/fetch_data.py",
    output:
        touch(f"{OUT_ROOT}/{{study}}/.fetch_done"),
    params:
        out_root=OUT_ROOT,
    conda:
        "env/fetch.yaml"
    shell:
        "python {input.script} --metadata {input.metadata} --out-root {params.out_root} download-study {wildcards.study}"
