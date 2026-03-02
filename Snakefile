# Project entrypoint Snakefile.

from pathlib import Path

configfile: str(Path(workflow.basedir) / "config" / "qc.yaml")

include: "workflow/modules/fetch/Snakefile"
include: "workflow/modules/qc_preproc/Snakefile"


# Choose what runs by default.
# - `snakemake`            -> runs both fetch + qc
# - `snakemake fetch_all`  -> runs only fetch module
# - `snakemake qc_all`     -> runs only qc module
rule all:
    input:
        rules.fetch_all.input,
        rules.qc_all.input,
