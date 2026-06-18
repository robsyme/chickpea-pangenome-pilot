# Platform & compute environment (reference)

Durable reference for where this project runs and stores data on Seqera Platform.
Not an attempt log — update in place when the infrastructure changes. Verified
2026-06-18 via `tw`.

## Seqera Platform

- Endpoint: `https://api.cloud.seqera.io` (Seqera Cloud). Authenticated user: `robsyme`.
- Interact via the `tw` CLI. `TOWER_ACCESS_TOKEN` is set in the environment.
- Set `export TOWER_WORKSPACE_ID=66005590816034` before workspace `tw` commands.

| Item | Value |
|------|-------|
| Organization | `robsyme-research` (ID `151253561315404`) |
| Workspace | `chickpea-pangenome-pilot` (ID `66005590816034`), PRIVATE |
| AWS credentials | `aws-scidev-playground` (ID `1PTV0DyBVrpmMMuOptx5Ha`), AWS acct `307946633589` |

## Compute environment

`aws-batch-default` (ID `SD36aYdcBFm9gvLmCqu2I`) — **primary**. AWS Batch (Forge):

| Setting | Value |
|---|---|
| Region | `us-east-1` |
| Provisioning | SPOT, `SPOT_CAPACITY_OPTIMIZED`, maxCpus 1000 |
| Instance types | `c6id`, `m6id`, `r6id` (all local NVMe; r6id for memory-heavy assembly) |
| Fusion v2 / Wave | enabled |
| Fusion Snapshots | enabled (spot-interruption recovery) |
| EBS auto-scale | off (NVMe fast instance storage instead) |
| Head-job env | `NXF_SYNTAX_PARSER=v2` |

The `tw` CLI has **no compute-env update**: to change settings, delete (async — wait
for `DELETING` to clear) and re-add. Recreate recipe:

```bash
export TOWER_WORKSPACE_ID=66005590816034
tw compute-envs add aws-batch forge \
  --name aws-batch-default \
  --credentials aws-scidev-playground \
  --region us-east-1 \
  --work-dir s3://scidev-playground-us-east-1/robsyme/chickpea-pilot/work \
  --allow-buckets s3://scidev-playground-us-east-1 \
  --max-cpus 1000 --provisioning-model SPOT --alloc-strategy SPOT_CAPACITY_OPTIMIZED \
  --instance-types c6id,m6id,r6id \
  --fusion-v2 --wave --fast-storage --snapshots --no-ebs-auto-scale \
  --env NXF_SYNTAX_PARSER=v2 \
  --wait AVAILABLE
tw compute-envs primary set --name aws-batch-default
```

## Storage (S3, us-east-1)

| Purpose | Path |
|---|---|
| Nextflow work dir | `s3://scidev-playground-us-east-1/robsyme/chickpea-pilot/work` |
| Published results (`outputDir`) | `s3://scidev-playground-us-east-1/robsyme/chickpea-pilot/results` |

The CE work bucket is `s3://scidev-playground-us-east-1` (allow-listed on the CE).
Reading S3 from a laptop needs a valid AWS SSO token (`aws sso login` if expired).

## Launchpad pipeline

`chickpea-gate-zero` — repo `https://github.com/robsyme/chickpea-pangenome-pilot`,
CE `aws-batch-default`. Config: `outputDir` (S3 results) + `dumpHashes`. The
**"Enable Nextflow syntax parser v2" toggle is ON** (UI-only; the `tw` CLI does not
expose it, and a CLI `tw pipelines update` will silently turn it off — edit this
pipeline's config in the UI). Launch by name:

```bash
export TOWER_WORKSPACE_ID=66005590816034
tw launch chickpea-gate-zero --name <run-name>
```

## Nextflow version

The pipeline requires **≥ 26.04.4** (fix for nextflow-io/nextflow#7226). The CE
otherwise provisions 26.04.3, so pin the version. Per-launch pin via a pre-run
script:

```bash
printf 'export NXF_VER=26.04.4\nexport NXF_SYNTAX_PARSER=v2\n' > prerun.sh
tw launch <pipeline-or-url> --compute-env aws-batch-default --pre-run prerun.sh ...
```

Durable options (TODO — pick one): add `NXF_VER=26.04.4` to the CE env (rebuild via
the recipe above), or add it to the Launchpad pre-run in the UI.

## Useful operations

- Run status: `tw runs view -i <runId>`
- Per-task metadata: `tw runs dump -i <runId> -o run.tar.gz` then read
  `workflow-tasks.json`. (The dump currently 500s on the `nextflow.log` section.)
- **Head-job Nextflow log** (incl. `dumpHashes` output): it is written to the work
  dir root as `s3://.../chickpea-pilot/work/nf-<runId>.log` — pull with `aws s3 cp`.
- **Default practice: resume, don't relaunch fresh** (see AGENTS.md "Conventions") —
  reuse cached tasks (esp. `FETCH_READS`) and build iteratively. Resume onto newer
  code: `tw runs relaunch -i <runId> --revision main --commit-id $(git rev-parse HEAD)`
  (`--pull-latest` alone keeps the original commit).
