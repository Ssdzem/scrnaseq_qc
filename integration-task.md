# Remote Orchestrator Task Template

Fill what you know. Leave unknown items as `UNKNOWN`. The remote-orchestrator agent should question missing/unclear fields before coding or remote execution.

## 1. Goal / Intent

What do you want the code to do, and why?

```text
I want to integrate the pre integration anndata object to correct batch effect by study accession and assay columns of the .obs of the input anndata. I want to incorporate an option to use GPU with scvi. I finally want to get the integrated output that corrects the batch effect of the study accession and assay values.  
```

## 2. Remote Test Ground

SSH host alias from your local WSL `~/.ssh/config`:

```text
ssh -p 43370 root@38.247.78.7 -L 8080:localhost:8080
```

Absolute remote sandbox path where the agent may read/write:

```text
/workspace/intanno/
```

Important: the agent must not write outside this sandbox.

## 3. Code Location

GitHub repo URL or local repo path:

```text
/home/sosad/intanno/ We are working in a worktree named intanno with branch name intanno and remote branch named origin/intanno. The remote sandbox is also a worktree called intanno with branch intanno and remote branch origin/intanno
```

Branch or commit to work from:

```text
local is named intanno, the remote is origin/intanno
```

Entry point, script, notebook, module, or command if known:

```text
a1_p1_qc/workflow/modules/integration_annotation, we will work mainly in this lcoation of the repo
```

## 4. Primary Input(s)

List each input file/object that already exists or should exist on the server.

| Name | Remote path inside sandbox | Format/type | What it must contain | Known concerns |
|---|---|---|---|---|
| pre_integrate_data | /workspace/intanno/a1_p1_qc/out/all_studies/all_pre_integration.h5ad | AnnData `.h5ad` | pre integrated scRNA-seq object | UNKNOWN |

## 5. Desired Code Behavior

Describe what the code should do with the inputs.

```text
From the pre integrated object, run an integration with scvi GPU that contemplates batches from study accession and assay obs columns values. The code should us the scvi option to integrate and incorporate gpu usage. The output should be an integrated anndata obj that corrects batch effect for said columns
```

## 6. Desired Output(s)

| Output name | Remote path inside sandbox | Format/type | Required contents |
|---|---|---|---|
| annotated data | /workspace/intanno/a1_p1_qc/out/integration_annotation/all_studies/all_integrated.h5ad | AnnData `.h5ad` | latent representation and cell type labels in obs metadata |

## 7. Validation Rules

How should the agent know the output is good enough?

Structural checks:

- [X] Output file exists.
- [X] Output loads with the expected reader/tool.
- [X] Required fields/metadata/representations exist.
- [X] Expected dimensions/counts are plausible.

Sanity / degeneracy checks:

- [X] Required labels are not all `NaN`, empty, or null.
- [X] Numeric representations are not all zero or all `NaN`.
- [X] The output does not obviously contradict the task intent.

Domain-specific checks supplied by user, if any:

```text
UNKNOWN
```

## 8. Known Method Requirements

List requirements you already know. Leave `UNKNOWN` if the agent should discover/question them.

```text
You need to tinker the config yaml of this snakefile process (/workspace/intanno/a1_p1_qc/config/integration_annotation.yam) so that it uses the scvi integration option and incorporate a gpu option for the user. You need to use snakemake --use-conda for conda generation, tinker it as necesary for including all the necessary packages. Optimize the gpu usage of the run. Run a dry run to verify the inputs and outputs, then do a real run. 
UNKNOWN
```

## 9. Constraints / Non-goals

Keep these unless you intentionally override them:

- Do not commit or push to GitHub.
- Do not edit SSH credentials or `~/.ssh/config`.
- Do not use `sudo` or admin access.
- Do not change services, daemons, ports, system packages, nginx, systemd, etc.
- Do not claim full scientific/domain correctness unless explicit criteria are provided.

Additional constraints:

```text
UNKNOWN
```

## 10. What To Do If Blocked

If required inputs are missing or invalid, the agent should:

```text
Stop, explain what is missing/invalid, explain why it matters, and ask me for the corrected input or decision.
```
