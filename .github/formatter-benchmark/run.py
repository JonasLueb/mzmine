#!/usr/bin/env python3
"""Prepare and check isolated, hash-verified Java formatter benchmark fixtures."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import traceback
import urllib.request
from pathlib import Path


ROOT = Path.cwd().resolve()
BENCHMARK = ROOT / ".github" / "formatter-benchmark"
CONFIG_PATH = BENCHMARK / "config.json"
INIT = BENCHMARK / "formatter.init.gradle"
COUNT = 50


def fail(message: str) -> None:
    raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inside(path: Path, root: Path = ROOT) -> Path:
    resolved = path.resolve()
    try:
        resolved.relative_to(root)
    except ValueError:
        fail(f"path escapes repository: {path}")
    return resolved


def below(path: Path, root: Path) -> Path:
    resolved = path.resolve()
    try:
        resolved.relative_to(root.resolve())
    except ValueError:
        fail(f"path escapes its fixture root: {path}")
    return resolved


def load_config() -> dict:
    if not CONFIG_PATH.is_file():
        fail(f"missing benchmark config: {CONFIG_PATH}")
    data = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    required = {"baseline", "idea_version", "idea_linux_url", "idea_linux_sha256", "idea_profile", "files"}
    missing = required - data.keys()
    if missing or not isinstance(data["files"], list) or len(data["files"]) != COUNT:
        fail(f"invalid benchmark config; missing={sorted(missing)}, files must contain {COUNT} entries")
    for entry in data["files"]:
        if set(("source", "sha256", "sample")) - entry.keys():
            fail("each config file requires source, sha256, and sample")
        source = inside(ROOT / entry["source"])
        if not source.is_file() or source.suffix != ".java":
            fail(f"invalid Java source: {entry['source']}")
        sample = Path(entry["sample"])
        if sample.is_absolute() or ".." in sample.parts or sample.suffix != ".java":
            fail(f"unsafe sample path: {entry['sample']}")
    return data


def java_worktree_clean(baseline: str) -> None:
    command = ["git", "diff", "--name-only", baseline, "--", "*.java"]
    changed = subprocess.check_output(command, cwd=ROOT, text=True).splitlines()
    if changed:
        fail("tracked Java differs from baseline: " + ", ".join(changed))


def verify_sources(config: dict) -> None:
    for entry in config["files"]:
        source = inside(ROOT / entry["source"])
        if sha256(source) != entry["sha256"].lower():
            fail(f"source SHA-256 changed: {entry['source']}")


def owner(source: Path) -> Path:
    relative = source.relative_to(ROOT)
    if len(relative.parts) < 3 or relative.parts[1] != "src":
        fail(f"source is not in a module src tree: {relative}")
    module = ROOT / relative.parts[0]
    if not module.is_dir():
        fail(f"source module does not exist: {module}")
    return module


def scratch_paths(config: dict) -> list[tuple[dict, Path]]:
    paths = []
    seen = set()
    for entry in config["files"]:
        source = inside(ROOT / entry["source"])
        scratch = owner(source) / "src" / "formatterBenchmark" / "java" / Path(entry["sample"])
        scratch = inside(scratch)
        if scratch in seen:
            fail(f"duplicate fixture destination: {scratch}")
        seen.add(scratch)
        paths.append((entry, scratch))
    return paths


def ensure_absent(paths: list[tuple[dict, Path]]) -> None:
    roots = scratch_roots(paths)
    existing = [path for path in roots if path.exists()]
    if existing:
        fail("refusing pre-existing formatter scratch directories: " + ", ".join(map(str, existing)))


def scratch_roots(paths: list[tuple[dict, Path]]) -> set[Path]:
    return {path.parents[len(Path(entry["sample"]).parts)] for entry, path in paths}


def cleanup(roots: set[Path]) -> None:
    for root in roots:
        if root.exists():
            shutil.rmtree(root)
        parent = root.parent
        if parent.exists() and not any(parent.iterdir()):
            parent.rmdir()


def copy_sources(paths: list[tuple[dict, Path]], source_root: Path | None = None) -> None:
    for entry, destination in paths:
        source = inside(ROOT / entry["source"]) if source_root is None else below(source_root / entry["fixture"], source_root)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)


def runner_temp() -> Path:
    return Path(os.environ.get("RUNNER_TEMP", tempfile.gettempdir())).resolve()


def download(url: str, destination: Path) -> None:
    with urllib.request.urlopen(url) as response, destination.open("wb") as stream:
        shutil.copyfileobj(response, stream)


def product_version(idea_home: Path) -> str:
    info = next((path for path in (idea_home / "product-info.json", idea_home / "Contents" / "Resources" / "product-info.json") if path.is_file()), None)
    if info is None:
        fail(f"missing IntelliJ product-info.json below: {idea_home}")
    return str(json.loads(info.read_text(encoding="utf-8"))["version"])


def provision_idea(config: dict) -> tuple[Path, float]:
    override = os.environ.get("FORMATTER_IDEA_BINARY")
    started = time.monotonic()
    if override:
        binary = Path(override).resolve()
        if not binary.is_file():
            fail(f"FORMATTER_IDEA_BINARY is not a file: {binary}")
        version = product_version(binary.parents[2])
        if version != config["idea_version"]:
            fail(f"IntelliJ version {version} does not equal {config['idea_version']}")
        return binary, time.monotonic() - started
    tools = Path(os.environ.get("FORMATTER_TOOLS_DIR", runner_temp() / "formatter-tools")).resolve()
    home = tools / f"idea-{config['idea_version']}"
    binary = home / "bin" / "idea.sh"
    if binary.is_file() and product_version(home) == config["idea_version"]:
        return binary, time.monotonic() - started
    tools.mkdir(parents=True, exist_ok=True)
    archive = tools / f"idea-{config['idea_version']}.tar.gz"
    download(config["idea_linux_url"], archive)
    expected = config["idea_linux_sha256"]
    if sha256(archive) != expected:
        fail("IntelliJ archive SHA-256 mismatch")
    staging = tools / f"idea-{config['idea_version']}.staging"
    shutil.rmtree(staging, ignore_errors=True)
    staging.mkdir()
    with tarfile.open(archive, "r:gz") as tar:
        tar.extractall(staging, filter="data")
    roots = [path for path in staging.iterdir() if (path / "bin" / "idea.sh").is_file()]
    if len(roots) != 1:
        fail("IntelliJ archive did not contain one IDEA home")
    shutil.rmtree(home, ignore_errors=True)
    roots[0].replace(home)
    shutil.rmtree(staging)
    if product_version(home) != config["idea_version"]:
        fail("downloaded IntelliJ product version differs from config")
    archive.unlink(missing_ok=True)
    return binary, time.monotonic() - started


def idea_properties() -> Path:
    root = runner_temp() / "formatter-idea-properties"
    root.mkdir(parents=True, exist_ok=True)
    for name in ("config", "system", "log", "plugins"):
        (root / name).mkdir(exist_ok=True)
    properties = root / "idea.properties"
    properties.write_text("\n".join(f"idea.{name}.path={root / name}" for name in ("config", "system", "log", "plugins")) + "\n", encoding="utf-8")
    return properties


def gradle_config(engine: str, config: dict, files: list[Path], idea: Path | None) -> Path:
    payload = {"engine": engine, "files": [str(path) for path in files], "ideaProfile": config["idea_profile"], "ideaBinary": str(idea) if idea else None}
    handle = tempfile.NamedTemporaryFile(mode="w", suffix=".json", prefix="formatter-benchmark-", delete=False, encoding="utf-8")
    with handle:
        json.dump(payload, handle)
    return Path(handle.name)


def run_gradle(engine: str, phase: str, config: dict, files: list[Path], idea: Path | None, results: Path) -> tuple[int, float, str]:
    runtime_config = gradle_config(engine, config, files, idea)
    task = "spotlessApply" if engine == "eclipse" and phase == "apply" else "spotlessCheck" if engine == "eclipse" else "benchmarkIntellijApply" if phase == "apply" else "benchmarkIntellijCheck"
    command = ["./gradlew", "--no-daemon", "--no-parallel", "--no-build-cache", "-I", str(INIT), f"-Dformatter.benchmark.config={runtime_config}", task]
    started = time.monotonic()
    environment = os.environ.copy()
    if engine == "intellij":
        environment["IDEA_PROPERTIES"] = str(idea_properties())
    completed = subprocess.run(command, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=environment)
    elapsed = time.monotonic() - started
    results.mkdir(parents=True, exist_ok=True)
    log = results / f"{engine}-{phase}.log"
    log.write_text(completed.stdout, encoding="utf-8")
    runtime_config.unlink(missing_ok=True)
    return completed.returncode, elapsed, completed.stdout


def load_manifest(fixtures: Path, engine: str, config: dict) -> list[dict]:
    manifest = json.loads((fixtures / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("engine") != engine or manifest.get("idea_version") != config["idea_version"] or len(manifest.get("files", [])) != COUNT:
        fail("fixture manifest engine, version, or count mismatch")
    return manifest["files"]


def write_result(destination: Path, result: dict) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", required=True, choices=("prepare", "check"))
    parser.add_argument("--engine", required=True, choices=("eclipse", "intellij"))
    parser.add_argument("--fixtures", required=True)
    parser.add_argument("--results", required=True)
    args = parser.parse_args()
    fixtures, results = Path(args.fixtures).resolve(), Path(args.results).resolve()
    config = load_config()
    source_sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    result = {
        "engine": args.engine,
        "phase": args.phase,
        "formatter_phase": os.environ.get("FORMATTER_PHASE", "unspecified"),
        "repeat": os.environ.get("FORMATTER_REPEAT", "unspecified"),
        "source_sha": source_sha,
        "exit": 1,
        "success": False,
        "valid": False,
        "provision_seconds": 0.0,
        "check_seconds": 0.0,
        "apply_seconds": 0.0,
        "target_count": 0,
        "idea_version": config["idea_version"],
        "files": [],
    }
    paths: list[tuple[dict, Path]] = []
    owned_roots: set[Path] = set()
    try:
        java_worktree_clean(config["baseline"])
        verify_sources(config)
        paths = scratch_paths(config)
        ensure_absent(paths)
        owned_roots = scratch_roots(paths)
        if args.phase == "check":
            manifest_files = load_manifest(fixtures, args.engine, config)
            by_source = {entry["sourcepath"]: entry for entry in manifest_files}
            if set(by_source) != {entry["source"] for entry, _ in paths}:
                fail("fixture manifest source set mismatch")
            for entry, _ in paths:
                manifest = by_source[entry["source"]]
                if manifest.get("source_sha256") != entry["sha256"] or manifest.get("sourcepath") != entry["source"] or manifest.get("fixture") != entry["sample"]:
                    fail(f"fixture manifest source mismatch: {entry['source']}")
                entry["fixture"] = manifest["fixture"]
                if sha256(below(fixtures / entry["fixture"], fixtures)) != manifest["expected_sha256"]:
                    fail(f"fixture hash mismatch: {entry['fixture']}")
            copy_sources(paths, fixtures)
        else:
            copy_sources(paths)
        files = [path for _, path in paths]
        idea = None
        if args.engine == "intellij":
            idea, duration = provision_idea(config)
            result["provision_seconds"] = duration
        before = {path: sha256(path) for path in files}
        if args.phase == "prepare":
            code, duration, output = run_gradle(args.engine, "apply", config, files, idea, results)
            result["apply_seconds"] = duration
            if code:
                fail(f"Gradle apply failed ({code}); see {results}")
        code, duration, output = run_gradle(args.engine, "check", config, files, idea, results)
        result["check_seconds"] = duration
        if code:
            fail(f"Gradle check failed ({code}); see {results}")
        after = {path: sha256(path) for path in files}
        if args.phase == "check" and before != after:
            fail("check changed benchmark fixtures")
        if args.engine == "eclipse":
            audited = {Path(line.split("=", 1)[1]).resolve() for line in output.splitlines() if line.startswith("FORMATTER_BENCHMARK_ECLIPSE_TARGET=")}
            if audited != set(files):
                fail("Eclipse target audit did not match the 50 fixture files")
            result["target_count"] = len(audited)
        if args.engine == "intellij" and len(files) != COUNT:
            fail("IntelliJ batch did not receive 50 fixture files")
        if args.engine == "intellij":
            scan = re.search(r"\b(\d+) file\(s\) scanned\.", output)
            if scan is None or int(scan.group(1)) != COUNT:
                fail("IntelliJ output did not report exactly 50 scanned files")
            result["target_count"] = int(scan.group(1))
        if args.phase == "prepare":
            if fixtures.exists():
                fail(f"fixtures already exist: {fixtures}")
            fixtures.mkdir(parents=True)
            manifest_files = []
            for entry, path in paths:
                fixture = Path(entry["sample"])
                output_path = fixtures / fixture
                output_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(path, output_path)
                manifest_files.append({"sourcepath": entry["source"], "source_sha256": entry["sha256"], "fixture": fixture.as_posix(), "expected_sha256": sha256(path)})
            (fixtures / "manifest.json").write_text(json.dumps({"engine": args.engine, "idea_version": config["idea_version"], "files": manifest_files}, indent=2) + "\n", encoding="utf-8")
        result.update({"exit": 0, "success": True, "valid": True, "files": [str(path.relative_to(ROOT)) for path in files], "fixture_hashes": {str(path.relative_to(ROOT)): sha256(path) for path in files}})
    except Exception as error:
        result["error"] = str(error)
        result["traceback"] = traceback.format_exc()
    finally:
        if owned_roots:
            cleanup(owned_roots)
        write_result(results, result)
    if result["exit"]:
        print(result.get("error", "benchmark failed"), file=sys.stderr)
    return result["exit"]


if __name__ == "__main__":
    raise SystemExit(main())
