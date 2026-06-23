rule figures_all:
    input:
        h5ad=_figures_h5ad(),
        metadata_csv=lambda wc: str(_figures_cfg().get("metadata_csv", "config/liver_scRNAseq_metadata_everything.csv")),
    output:
        figures=P1_FIGURE_OUTPUTS,
        done=touch("out/figures/.figures_done"),
    log:
        "out/logs/figures/figures.log"
    conda:
        "../../exploration/env/a1p2def.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        {{
          echo "[$(date -Is)] [figures_all]"
          echo "[$(date -Is)] INPUT h5ad={input.h5ad:q}"
          echo "[$(date -Is)] INPUT metadata_csv={input.metadata_csv:q}"
          echo "[$(date -Is)] OUTPUT figures={output.figures:q}"
          $CONDA_PREFIX/bin/python workflow/modules/figures/scripts/make_figures.py \
            --h5ad {input.h5ad:q} \
            --metadata-csv {input.metadata_csv:q} \
            --out-dir "out/figures" \
            --manifest "out/figures/manifest.json"
        }} 2>&1 | tee -a {log:q}
        """
