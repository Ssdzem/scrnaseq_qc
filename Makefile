.DEFAULT_GOAL := help

SHELL := /usr/bin/bash

PROFILE ?= profiles/local
SNAKEFILE_FETCH := workflow/modules/fetch/Snakefile

CORES ?= 6
TARGET ?=

.PHONY: help
help:
	@echo "Targets:"
	@echo "  make fetch-study STUDY=<study_accession> [CORES=6]"
	@echo "  make fetch-run STUDY=<study_accession> DONOR=<donor_id> SAMPLE=<sample_id> RUN=<run_id> [CORES=1]"
	@echo "  make fetch-all [CORES=6]"
	@echo "  make dry TARGET=<snakemake-target> [CORES=6]"
	@echo ""
	@echo "Notes:"
	@echo "  - Uses --profile $(PROFILE) by default."
	@echo "  - For Aspera downloads, set config/fetch.yaml: aspera.enabled=true and paths."


.PHONY: fetch-study
fetch-study:
	@if [[ -z "$(STUDY)" ]]; then echo "Missing STUDY=<study_accession>"; exit 2; fi
	snakemake --profile $(PROFILE) --cores $(CORES) -s $(SNAKEFILE_FETCH) out/$(STUDY)/.fetch_done


.PHONY: fetch-run
fetch-run:
	@if [[ -z "$(STUDY)" || -z "$(DONOR)" || -z "$(SAMPLE)" || -z "$(RUN)" ]]; then \
		echo "Usage: make fetch-run STUDY=<study_accession> DONOR=<donor_id> SAMPLE=<sample_id> RUN=<run_id> [CORES=1]"; \
		exit 2; \
	fi
	snakemake --profile $(PROFILE) --cores $(CORES) -s $(SNAKEFILE_FETCH) out/$(STUDY)/fastqs/$(DONOR)/$(SAMPLE)/$(RUN)/.fetch_done


.PHONY: fetch-all
fetch-all:
	snakemake --profile $(PROFILE) --cores $(CORES) -s $(SNAKEFILE_FETCH) fetch_all


.PHONY: dry
dry:
	@if [[ -z "$(TARGET)" ]]; then echo "Missing TARGET=<snakemake-target>"; exit 2; fi
	snakemake --profile $(PROFILE) --cores $(CORES) -n -p -s $(SNAKEFILE_FETCH) $(TARGET)
