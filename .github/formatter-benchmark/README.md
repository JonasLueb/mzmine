# Formatter CI benchmark

This is an experiment for choosing the Java formatter. It does not establish a
runtime winner until **both engines complete** on GitHub Actions. The candidate
Spotless integration remains separate from the decision to adopt its output.

Push these files to a feature branch to run **Formatter benchmark**. No merge to
`master` is needed. Subsequent manual runs are available once GitHub knows the
workflow on the default branch. Re-run **all jobs**, not only failed jobs: each
run attempt deliberately uses a new cache key.

## What is measured

Both engines check the same deterministic selection of 50 real Java files from
this repository. An unmeasured preparation job formats separate copies with each
engine, then verifies them. This models a successful check on code already using
the chosen formatter, rather than timing one engine's expected formatting errors.
The source hashes, fixed Git baseline, canonical hashes and exact target count
are checked. Temporary copies live in a separate `src/formatterBenchmark/java`
directory and are removed after the run. Production Java files are not rewritten.

| Scenario | Measurements | What the runner starts with |
|---|---:|---|
| Cold | 1 per engine | Empty Gradle and formatter-tool caches |
| Warm | 3 per engine | Exact cache saved by that engine's cold job |

Every observation runs on a fresh `ubuntu-latest` GitHub-hosted VM with Java
26.0.1, a fresh Gradle JVM, no Gradle build cache and no parallel Gradle tasks.
The two engines run in separate jobs; their caches never mix. Warm cache misses
fail rather than silently being reported as warm timings. IntelliJ 2026.2.3 is
downloaded from JetBrains and SHA-256 verified; no interactive IDE installation
or paid IDE license is supplied. Eclipse uses the proposed pinned Spotless/JDT
configuration. IntelliJ uses one batch command from a Gradle task and the existing
IntelliJ XML export. The batch task is an experimental integration, not a finished
replacement for the Spotless integration. Its file list is preselected; a production
IntelliJ integration would also need Git change selection. Eclipse includes its
existing ratchet scan, so these compare candidate workflows, not engine speed alone.

Read **Compare workflow durations** in the run summary or download the
`formatter-comparison` artifact. It reports total job duration, Gradle check time,
IntelliJ download/extraction, cache restoration and cache saving separately.
Warm results show median and range. The total includes checkout, Java setup,
cache transfer, benchmark fixture/result transfer and cleanup. Queue time,
canonical fixture preparation and report generation are excluded. Therefore it
is the measured duration of a dedicated formatting job, not the whole existing
build/test pipeline. The artifact retains individual results and GitHub's job
and step timestamps. Failed or missing runs make the comparison incomplete.

A single cold observation is descriptive, not a precise distribution of network
performance. Repeat the whole workflow if the decision depends on a small gap.
Runner CPU allocation, hosted image, package-server latency and cache throughput
can vary. Do not substitute these numbers for native interactive IDE latency.

## Scope and cleanup

The fixed manifest is intentionally tied to this baseline; source changes fail
validation instead of silently changing the workload. Update the manifest and
baseline together for a new experiment. The full-source formatting comparison
and the two known Eclipse failures must be considered alongside these timings;
passing the 50-file sample does not resolve those failures.

Each run creates two dependency/tool caches and temporary result artifacts.
The IntelliJ cache is large. After downloading the comparison, remove caches
whose keys begin `formatter-benchmark-v1-<run-id>-` from the repository's Actions
cache page. Results expire after 14 days; canonical fixture artifacts after 7.
Private repository jobs consume the account's normal Actions allowance. This
benchmark does not provision a persistent runner or change the billing plan.
