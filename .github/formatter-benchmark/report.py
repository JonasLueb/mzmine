#!/usr/bin/env python3
"""Summarize observed Actions job durations; never substitute estimates for missing runs."""
import argparse
from datetime import datetime
import json
import math
import os
from pathlib import Path
import re
import statistics
import urllib.request


def seconds(start, end):
    if not start or not end:
        return None
    return (datetime.fromisoformat(end.replace('Z', '+00:00')) -
            datetime.fromisoformat(start.replace('Z', '+00:00'))).total_seconds()


def fetch_jobs():
    base = os.environ.get('GITHUB_API_URL', 'https://api.github.com')
    repo, run, attempt = (os.environ[k] for k in
                          ('GITHUB_REPOSITORY', 'GITHUB_RUN_ID', 'GITHUB_RUN_ATTEMPT'))
    jobs = []
    for page in range(1, 20):
        url = f'{base}/repos/{repo}/actions/runs/{run}/attempts/{attempt}/jobs?per_page=100&page={page}'
        request = urllib.request.Request(url, headers={
            'Authorization': 'Bearer ' + os.environ['GH_TOKEN'],
            'Accept': 'application/vnd.github+json',
            'X-GitHub-Api-Version': '2022-11-28',
        })
        with urllib.request.urlopen(request, timeout=60) as response:
            batch = json.load(response)['jobs']
        jobs.extend(batch)
        if len(batch) < 100:
            return jobs
    raise RuntimeError('Unexpected job count; pagination incomplete')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--jobs-json', type=Path, help='Offline validation with saved API response')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    jobs = json.loads(args.jobs_json.read_text()) if args.jobs_json else fetch_jobs()
    if isinstance(jobs, dict):
        jobs = jobs['jobs']
    (args.output/'jobs.json').write_text(json.dumps(jobs, indent=2)+'\n')
    records, problems = [], []
    for job in jobs:
        match = re.search(r'bench-(cold|warm)-(eclipse|intellij)-(\d+)$', job['name'])
        if not match:
            continue
        phase, engine, repeat = match.groups()
        artifact = args.results / f'formatter-result-{phase}-{engine}-{repeat}'
        result_path = artifact / 'results.json'
        result = json.loads(result_path.read_text()) if result_path.exists() else {}
        duration = seconds(job.get('started_at'), job.get('completed_at'))
        valid = (job.get('conclusion') == 'success' and result.get('success') is True
                 and result.get('target_count') == 50 and duration is not None
                 and result.get('engine') == engine and result.get('formatter_phase') == phase
                 and result.get('source_sha') == job.get('head_sha'))
        if not valid:
            problems.append(f'{phase}/{engine}/{repeat}: job or 50-file measurement failed/missing')
        steps = {step['name']: seconds(step.get('started_at'), step.get('completed_at'))
                 for step in job.get('steps', []) if step.get('conclusion') != 'skipped'}
        records.append(dict(phase=phase, engine=engine, repeat=int(repeat), valid=valid,
                            job_seconds=duration, job_url=job.get('html_url'),
                            billed_minutes=math.ceil(duration/60) if duration else None,
                            steps_seconds=steps, measurement=result))
    expected = {(p,e,r) for p, rs in [('cold',[1]),('warm',[1,2,3])]
                for e in ['eclipse','intellij'] for r in rs}
    seen = {(r['phase'],r['engine'],r['repeat']) for r in records}
    for missing in sorted(expected-seen):
        problems.append('Missing job: ' + '/'.join(map(str,missing)))
    if len(seen) != len(records):
        problems.append('Duplicate benchmark jobs')
    lines = ['# Formatter CI comparison', '',
             '50 real Java files per engine, each already formatted by that engine. '
             'Cold: one observation with empty Gradle/tool caches. Warm: three fresh runners '
             'restoring the exact cache produced by the cold job.', '',
             '| Engine | Cache | Valid runs | Whole job | Gradle check | Tool provisioning | Cache restore | Cache save |',
             '|---|---|---:|---:|---:|---:|---:|---:|']
    for engine in ['eclipse','intellij']:
        for phase, count in [('cold',1),('warm',3)]:
            group = [r for r in records if r['engine']==engine and r['phase']==phase and r['valid']]
            def format_values(getter):
                values = [getter(r) for r in group]
                if len(values) != count or any(v is None for v in values):
                    return 'incomplete'
                if count == 1:
                    return f'{values[0]:.1f}s'
                return f'{statistics.median(values):.1f}s ({min(values):.1f}–{max(values):.1f})'
            job_time = format_values(lambda r:r['job_seconds'])
            check = format_values(lambda r:r['measurement'].get('check_seconds'))
            provision = format_values(lambda r:r['measurement'].get('provision_seconds'))
            restore = format_values(lambda r:r['steps_seconds'].get('Restore formatter cache')) if phase=='warm' else '—'
            save = format_values(lambda r:r['steps_seconds'].get('Save formatter cache')) if phase=='cold' else '—'
            lines.append(f'| {engine} | {phase} | {len(group)}/{count} | {job_time} | {check} | {provision} | {restore} | {save} |')
    lines += ['', 'Warm values are median (minimum–maximum). Whole job is GitHub’s observed '
              'start-to-completion time, including checkout, Java setup, artifact transfer, '
              'cache transfer and cleanup; queue time is excluded. Individual columns overlap: '
              'provisioning and Gradle check are inside the measured script step. '
              'Cold Gradle check includes Gradle/dependency downloads; tool provisioning measures '
              'the separate IntelliJ distribution download/extraction. These are dedicated formatting jobs, '
              'not the duration of the entire existing build/test pipeline.', '',
              'Preparing canonical fixtures and generating this report are benchmark setup, excluded '
              'from the comparison. Both engines use the same 50 source selections; their canonical '
              'input text differs. The IntelliJ path uses one batch process per Gradle task. '
              'The Eclipse path uses the proposed Spotless task and its ratchet scan. '
              'IntelliJ receives a preselected list; production Git change selection is not included. '
              'Full-source formatting compatibility is a separate measurement.', '']
    if problems:
        lines += ['**Incomplete: do not use this run to choose an engine.**', '']
        lines += [f'- {problem}' for problem in problems]
    lines += ['', '## Individual runs', '', '| Engine | Cache | Repeat | Whole job | Status |',
              '|---|---|---:|---:|---|']
    for record in sorted(records,key=lambda r:(r['engine'],r['phase'],r['repeat'])):
        duration = record['job_seconds']
        label = f'{duration:.1f}s' if duration is not None else 'missing'
        url = record['job_url']
        link = f'[{label}]({url})' if url else label
        lines.append(f"| {record['engine']} | {record['phase']} | {record['repeat']} | {link} | {'valid' if record['valid'] else 'invalid'} |")
    text = '\n'.join(lines)+'\n'
    (args.output/'comparison.md').write_text(text)
    (args.output/'comparison.json').write_text(json.dumps(dict(records=records,problems=problems),indent=2)+'\n')
    if os.environ.get('GITHUB_STEP_SUMMARY'):
        with open(os.environ['GITHUB_STEP_SUMMARY'],'a') as summary:
            summary.write(text)
    print(text)
    return 1 if problems else 0


if __name__ == '__main__':
    raise SystemExit(main())
