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

**Separation of concerns:** the workflow in Git stays infrastructure-agnostic; all
infra-specific settings live on the Launchpad pipeline / compute environment so the
repo stays portable. That means `outputDir` (the S3 results prefix), the work dir,
and allowed buckets are **not** in `nextflow.config` — they are set on the Launchpad
pipeline (`outputDir`) and the CE (`workDir`, `allowBuckets`). Consequence: results
publish to S3 only when launched **through the Launchpad pipeline**; a bare
`tw launch <git-url>` would default `outputDir` to `./results` in the ephemeral
head-job container and lose them — so pass infra config explicitly for ad-hoc runs,
or just use the Launchpad. For a local run, pass `-output-dir <path>`.

The CE work bucket is `s3://scidev-playground-us-east-1` (allow-listed on the CE).
Reading S3 from a laptop needs a valid AWS SSO token (`aws sso login` if expired).

## Launchpad pipeline

`chickpea-gate-zero` — repo `https://github.com/robsyme/chickpea-pangenome-pilot`,
CE `aws-batch-default`. Config: `outputDir` (S3 results prefix — infra-specific, so
it lives here, not in the repo). Its **pre-run script** exports `NXF_VER=26.04.4`
and `NXF_SYNTAX_PARSER=v2`. The
**"Enable Nextflow syntax parser v2" toggle is ON** (UI-only; the `tw` CLI does not
expose it, and a CLI `tw pipelines update` will silently turn it off — edit this
pipeline's config in the UI). Launch by name:

```bash
export TOWER_WORKSPACE_ID=66005590816034
tw launch chickpea-gate-zero --name <run-name>
```

## Nextflow version

The pipeline requires **≥ 26.04.4** (fix for nextflow-io/nextflow#7226). The CE
otherwise provisions 26.04.3, so the version is pinned.

**Durable pin (in place):** the `chickpea-gate-zero` Launchpad pipeline's **pre-run
script** exports both `NXF_VER=26.04.4` and `NXF_SYNTAX_PARSER=v2`, so
`tw launch chickpea-gate-zero` (and resumes of its runs) use 26.04.4 automatically.

For **ad-hoc** `tw launch <git-url>` runs that bypass the Launchpad pipeline, pass
the same pre-run explicitly:

```bash
printf 'export NXF_VER=26.04.4\nexport NXF_SYNTAX_PARSER=v2\n' > prerun.sh
tw launch <git-url> --compute-env aws-batch-default --pre-run prerun.sh ...
```

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
