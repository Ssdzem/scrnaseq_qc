"""Fetch helper for Snakemake.

Reads the project metadata CSV and downloads assets referenced by columns:

- `fastq_r1`, `fastq_r2` -> data_nature `fastqs`
- `raw_counts` -> data_nature `raw`
- `processed_counts` -> data_nature `qconly`

Files are written under:

`out/{study_accession}/{data_nature}/{donor_id}/{sample_id}/{run_id}/<files>`

This module intentionally does not attempt to unpack/reshape study-level
archives (e.g. GEO `*_RAW.tar`). It stores what it downloads deterministically;
downstream workflow steps can unpack and standardize as needed.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import time
import shlex
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set


def _dir_total_bytes(path: Path) -> int:
    total = 0
    for root, _, files in os.walk(path):
        for name in files:
            try:
                total += (Path(root) / name).stat().st_size
            except OSError:
                # Files can disappear/change while we traverse.
                continue
    return total


def _normalize_header(name: str) -> str:
    return re.sub(r"\s+", " ", name).strip()


def _open_csv(path: Path):
    return path.open(newline="", encoding="utf-8", errors="replace")


DATA_COLUMNS = {
    "fastqs": ("fastq_r1", "fastq_r2"),
    "raw": ("raw_counts",),
    "qconly": ("processed_counts",),
}


def _is_missing(value: str) -> bool:
    return (value or "").strip() == "" or (value or "").strip().lower() == "na"


def read_rows(metadata_csv: Path) -> List[Dict[str, str]]:
    with _open_csv(metadata_csv) as handle:
        reader = csv.reader(handle)
        try:
            raw_headers = next(reader)
        except StopIteration:
            raise ValueError(f"Empty CSV: {metadata_csv}")

        headers = [_normalize_header(h) for h in raw_headers]
        index = {h: i for i, h in enumerate(headers)}

        required = {"study_accession", "donor_id", "sample_id", "run_id"}
        missing = required - set(index)
        if missing:
            raise ValueError(
                f"Missing columns {sorted(missing)} in {metadata_csv}. Found: {headers}"
            )

        rows: List[Dict[str, str]] = []
        for raw_row in reader:
            if not raw_row:
                continue
            row: Dict[str, str] = {}
            for h, i in index.items():
                row[h] = raw_row[i].strip() if i < len(raw_row) else ""
            rows.append(row)
    return rows


def list_studies(metadata_csv: Path) -> List[str]:
    studies: Set[str] = set()
    for row in read_rows(metadata_csv):
        study = (row.get("study_accession") or "").strip()
        if study and not _is_missing(study):
            studies.add(study)
    return sorted(studies)


def _truthy(val: Optional[str]) -> bool:
    return (val or "").strip().lower() in {"1", "true", "yes", "y", "on"}


def _infer_gsa_run_id_from_url(url: str) -> Optional[str]:
    m = re.search(r"/(HRR\d{6,})/", url)
    if m:
        return m.group(1)
    return None


def _aspera_download(
    ascp_path: str,
    key_path: str,
    host: str,
    remote_root: str,
    study: str,
    run_id: str,
    dest_dir: Path,
    port: int = 33001,
    rate_limit: str = "100m",
    retries: int = 8,
    retry_sleep_seconds: int = 5,
    slow_retry_sleep_seconds: int = 45,
    slow_retries: int = 3,
    retry_sleep_schedule_seconds: Optional[List[int]] = None,
    stall_window_seconds: int = 15 * 60,
    stall_min_bytes: int = 50 * 1024 * 1024,
    stall_check_interval_seconds: int = 30,
    term_grace_seconds: int = 10,
) -> None:
    dest_dir.mkdir(parents=True, exist_ok=True)

    # Trailing slash on source is important: it downloads directory contents.
    src = f"{host}:{remote_root}/{study}/{run_id}/"
    cmd = [
        ascp_path,
        "-P",
        str(port),
        "-i",
        key_path,
        "-QT",
        "-l",
        rate_limit,
        "-k1",
        "-d",
        src,
        str(dest_dir),
    ]

    print(f"[fetch] ascp command: {shlex.join(cmd)}", flush=True)

    # Large transfers can time out (e.g. Error 36). In practice, re-running ascp
    # typically resumes and completes (-k1). Add retries so snakemake jobs are
    # robust to transient network hiccups.
    attempt = 0

    # Custom multi-tier retry schedule: values are the sleep time *before* the
    # next retry. Total attempts = 1 initial try + len(schedule) retries.
    if retry_sleep_schedule_seconds is None:
        # Backwards-compatible default: (retries - slow_retries - 1) quick
        # retries at 5s, then slow_retries at 45s.
        max_attempts = max(1, retries)
        quick_retries = max(0, max_attempts - 1 - max(0, slow_retries))
        schedule = ([retry_sleep_seconds] * quick_retries) + (
            [slow_retry_sleep_seconds] * max(0, slow_retries)
        )
    else:
        schedule = list(retry_sleep_schedule_seconds)
        max_attempts = 1 + len(schedule)

    while True:
        try:
            proc = subprocess.Popen(cmd)

            last_bytes = _dir_total_bytes(dest_dir)
            last_progress_at = time.monotonic()

            while True:
                try:
                    returncode = proc.wait(timeout=stall_check_interval_seconds)
                    if returncode != 0:
                        raise subprocess.CalledProcessError(returncode, cmd)
                    return
                except subprocess.TimeoutExpired:
                    current_bytes = _dir_total_bytes(dest_dir)
                    if current_bytes >= last_bytes + stall_min_bytes:
                        last_bytes = current_bytes
                        last_progress_at = time.monotonic()
                        continue

                    idle_for = time.monotonic() - last_progress_at
                    if idle_for >= stall_window_seconds:
                        print(
                            f"[fetch] ascp appears stalled: < {stall_min_bytes} bytes in "
                            f"{stall_window_seconds}s; terminating pid={proc.pid}",
                            flush=True,
                        )
                        proc.terminate()
                        try:
                            proc.wait(timeout=term_grace_seconds)
                        except subprocess.TimeoutExpired:
                            proc.kill()
                            proc.wait()
                        raise subprocess.CalledProcessError(124, cmd)
        except subprocess.CalledProcessError as e:
            attempt += 1
            if attempt >= max_attempts:
                raise

            sleep_for = schedule[min(attempt - 1, len(schedule) - 1)]
            print(
                f"[fetch] ascp failed (attempt {attempt}/{max_attempts}) with exit code {e.returncode}; "
                f"retrying in {sleep_for}s...",
                flush=True,
            )
            time.sleep(sleep_for)

def _dest_dir(out_root: Path, row: Dict[str, str], data_nature: str) -> Path:
    return (
        out_root
        / row["study_accession"]
        / data_nature
        / row["donor_id"]
        / row["sample_id"]
        / row["run_id"]
    )


def _ensure_run_dir(out_root: Path, row: Dict[str, str], data_nature: str = "fastqs") -> Path:
    dest_dir = _dest_dir(out_root, row, data_nature)
    dest_dir.mkdir(parents=True, exist_ok=True)
    return dest_dir


def _dest_from_url(url: str, outdir: Path) -> Path:
    name = url.rstrip("/").split("/")[-1]
    if not name:
        raise ValueError(f"Could not infer filename from url: {url}")
    return outdir / name


def _best_level_for_url(url: str, row: Dict[str, str]) -> str:
    """Return the most specific hierarchy level implied by `url`.

    Levels from most->least specific:
    run_id < sample_id < donor_id < study_accession

    If none of the tokens appear in the URL, returns "unknown".
    """

    url_lc = url.lower()

    def _token_present(token: str) -> bool:
        token_lc = token.lower()
        # Match token as a standalone chunk bounded by non-alphanumeric chars.
        # This avoids false positives like donor_id="3" matching inside
        # ".../<study>/...".
        pattern = rf"(?<![a-z0-9]){re.escape(token_lc)}(?![a-z0-9])"
        return re.search(pattern, url_lc) is not None

    for level in ("run_id", "sample_id", "donor_id", "study_accession"):
        token = (row.get(level) or "").strip()
        if token and not _is_missing(token) and _token_present(token):
            return level
    return "unknown"


def _dest_dir_with_level(out_root: Path, row: Dict[str, str], data_nature: str, level: str) -> Path:
    base = out_root / row["study_accession"] / data_nature
    if level == "study_accession":
        return base
    if level == "donor_id":
        return base / row["donor_id"]
    if level == "sample_id":
        return base / row["donor_id"] / row["sample_id"]
    if level == "run_id":
        return base / row["donor_id"] / row["sample_id"] / row["run_id"]
    raise ValueError(f"Unknown level: {level}")


def _curl_download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)

    # If a previous attempt created an empty file, remove it so curl can resume
    # cleanly.
    try:
        if dest.exists() and dest.stat().st_size == 0:
            dest.unlink()
    except OSError:
        pass

    cmd = [
        "curl",
        "-L",
        "--fail",
        "--retry",
        # Large GEO downloads occasionally fail mid-transfer; be more patient.
        "10",
        "--retry-delay",
        "2",
        "--connect-timeout",
        "30",
        "--max-time",
        "0",
        "-C",
        "-",
        "-o",
        str(dest),
        url,
    ]
    # `--retry-all-errors` only exists on newer curl; pass it conditionally so we
    # don't hard-fail on older binaries.
    if _curl_supports_retry_all_errors():
        cmd.insert(7, "--retry-all-errors")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        # curl exit 22 == HTTP >= 400. With `-C -` (resume), a fully downloaded
        # file can still trigger a spurious 416/22 on the resume probe; treat a
        # size-complete file as success.
        if exc.returncode == 22 and dest.exists():
            remote_size = _http_content_length(url)
            try:
                local_size = dest.stat().st_size
            except OSError:
                local_size = -1
            if remote_size is not None and local_size == remote_size:
                return
        raise


@lru_cache(maxsize=1)
def _curl_supports_retry_all_errors() -> bool:
    try:
        result = subprocess.run(
            ["curl", "--help", "all"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return False
    return "--retry-all-errors" in result.stdout


def _http_content_length(url: str) -> Optional[int]:
    try:
        result = subprocess.run(
            ["curl", "-sSIL", url],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None

    content_length: Optional[int] = None
    for line in result.stdout.splitlines():
        if line.lower().startswith("content-length:"):
            _, _, value = line.partition(":")
            value = value.strip()
            try:
                content_length = int(value)
            except ValueError:
                continue
    return content_length


def _resolve_aspera_param(val: str, env_key: str) -> str:
    if (val or "").strip():
        return os.path.expanduser(val)
    return os.path.expanduser(os.environ.get(env_key, ""))


def cmd_download_study(
    metadata_csv: Path,
    out_root: Path,
    study: str,
    use_aspera: bool = False,
    ascp_path: str = "",
    key_path: str = "",
    aspera_host: str = "",
    aspera_remote_root: str = "",
    aspera_port: int = 33001,
    aspera_rate_limit: str = "100m",
) -> None:
    rows = [r for r in read_rows(metadata_csv) if r.get("study_accession") == study]
    if not rows:
        raise SystemExit(
            f"Study '{study}' not found in {metadata_csv}. Available: {list_studies(metadata_csv)}"
        )

    if use_aspera:
        if not ascp_path or not key_path:
            raise SystemExit(
                "Aspera requested but ascp/key paths missing. "
                "Set them in config/fetch.yaml or env vars ASCP_PATH/ASCP_KEY."
            )
        if not Path(ascp_path).exists():
            raise SystemExit(f"ascp not found: {ascp_path}")
        if not Path(key_path).exists():
            raise SystemExit(f"Aspera key not found: {key_path}")
        if not aspera_host:
            raise SystemExit("Aspera host missing (config.aspera.host)")
        if not aspera_remote_root:
            raise SystemExit("Aspera remote_root missing (config.aspera.remote_root)")

        run_ids: Set[str] = set()
        for row in rows:
            rid = (row.get("run_id") or "").strip()
            if not rid or _is_missing(rid):
                continue
            run_ids.add(rid)
            for col in ("fastq_r1", "fastq_r2"):
                url = (row.get(col) or "").strip()
                if not _is_missing(url):
                    inferred = _infer_gsa_run_id_from_url(url)
                    if inferred and inferred != rid:
                        raise SystemExit(
                            f"Run id mismatch for study {study}: metadata run_id={rid} but url implies {inferred}: {url}"
                        )

        if not run_ids:
            raise SystemExit(f"No run_id values found for study {study}; cannot use Aspera")

        # Download per-run directories; snakemake parallelism can schedule
        # multiple runs at once.
        retry_schedule = (
            [5] * 3
            + [15] * 3
            + [60] * 3
            + [360] * 2
            + [720]
        )

        for rid in sorted(run_ids):
            dest_dir = out_root / study / "fastqs" / rid
            _aspera_download(
                ascp_path=ascp_path,
                key_path=key_path,
                host=aspera_host,
                remote_root=aspera_remote_root,
                study=study,
                run_id=rid,
                dest_dir=dest_dir,
                port=aspera_port,
                rate_limit=aspera_rate_limit,
                retry_sleep_schedule_seconds=retry_schedule,
            )
        return

    # Curl mode: download all runs belonging to the study using the
    # per-run downloader (which applies URL-level hierarchy heuristics and
    # deduplication across data_nature columns).
    run_ids: Set[str] = set()
    for row in rows:
        rid = (row.get("run_id") or "").strip()
        if rid and not _is_missing(rid):
            run_ids.add(rid)

    if not run_ids:
        raise SystemExit(f"No run_id values found for study {study}")

    for rid in sorted(run_ids):
        cmd_download_run(
            metadata_csv=metadata_csv,
            out_root=out_root,
            study=study,
            run_id=rid,
            use_aspera=False,
            ascp_path="",
            key_path="",
            aspera_host="",
            aspera_remote_root="",
            aspera_port=aspera_port,
            aspera_rate_limit=aspera_rate_limit,
        )


def cmd_download_run(
    metadata_csv: Path,
    out_root: Path,
    study: str,
    run_id: str,
    use_aspera: bool = False,
    ascp_path: str = "",
    key_path: str = "",
    aspera_host: str = "",
    aspera_remote_root: str = "",
    aspera_port: int = 33001,
    aspera_rate_limit: str = "100m",
) -> None:
    rows = [
        r
        for r in read_rows(metadata_csv)
        if r.get("study_accession") == study and (r.get("run_id") or "").strip() == run_id
    ]
    if not rows:
        raise SystemExit(
            f"Run '{run_id}' for study '{study}' not found in {metadata_csv}."
        )

    if use_aspera:
        if not ascp_path or not key_path:
            raise SystemExit(
                "Aspera requested but ascp/key paths missing. "
                "Set them in config/fetch.yaml or env vars ASCP_PATH/ASCP_KEY."
            )
        if not Path(ascp_path).exists():
            raise SystemExit(f"ascp not found: {ascp_path}")
        if not Path(key_path).exists():
            raise SystemExit(f"Aspera key not found: {key_path}")
        if not aspera_host:
            raise SystemExit("Aspera host missing (config.aspera.host)")
        if not aspera_remote_root:
            raise SystemExit("Aspera remote_root missing (config.aspera.remote_root)")

        for row in rows:
            for col in ("fastq_r1", "fastq_r2"):
                url = (row.get(col) or "").strip()
                if not _is_missing(url):
                    inferred = _infer_gsa_run_id_from_url(url)
                    if inferred and inferred != run_id:
                        raise SystemExit(
                            f"Run id mismatch for study {study}: run_id={run_id} but url implies {inferred}: {url}"
                        )

        # Match the project output hierarchy.
        dest_dir = _ensure_run_dir(out_root, rows[0], data_nature="fastqs")
        retry_schedule = (
            [5] * 3
            + [15] * 3
            + [60] * 3
            + [360] * 2
            + [720]
        )
        _aspera_download(
            ascp_path=ascp_path,
            key_path=key_path,
            host=aspera_host,
            remote_root=aspera_remote_root,
            study=study,
            run_id=run_id,
            dest_dir=dest_dir,
            port=aspera_port,
            rate_limit=aspera_rate_limit,
            retry_sleep_schedule_seconds=retry_schedule,
        )
        return

    # Precompute a fallback "level" per URL based on how unique it is across
    # different hierarchy keys. This handles study-level links that do not
    # include tokens like sample_id/run_id in the URL.
    url_to_rows: Dict[str, List[Dict[str, str]]] = {}
    for row in rows:
        for _, cols in DATA_COLUMNS.items():
            for col in cols:
                url = (row.get(col) or "").strip()
                if _is_missing(url):
                    continue
                url_to_rows.setdefault(url, []).append(row)

    url_fallback_level: Dict[str, str] = {}
    for url, urows in url_to_rows.items():
        for level in ("run_id", "sample_id", "donor_id", "study_accession"):
            values = {(r.get(level) or "").strip() for r in urows}
            values = {v for v in values if v and not _is_missing(v)}
            if len(values) == 1:
                url_fallback_level[url] = level
                break
        else:
            url_fallback_level[url] = "study_accession"

    seen: Set[tuple[str, Path]] = set()
    for row in rows:
        for data_nature, cols in DATA_COLUMNS.items():
            for col in cols:
                url = (row.get(col) or "").strip()
                if _is_missing(url):
                    continue

                level = _best_level_for_url(url, row)
                if level == "unknown":
                    level = url_fallback_level.get(url, "study_accession")

                dest_dir = _dest_dir_with_level(out_root, row, data_nature, level)
                dest = _dest_from_url(url, dest_dir)
                key = (url, dest)
                if key in seen:
                    continue
                seen.add(key)
                _curl_download(url, dest)


def cmd_list_studies(metadata_csv: Path) -> None:
    for s in list_studies(metadata_csv):
        print(s)


def cmd_list_studies_json(metadata_csv: Path) -> None:
    import json

    print(json.dumps(list_studies(metadata_csv)))


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--metadata",
        default="config/samples.csv",
        help="Metadata CSV with study_accession/raw_uri/processed_uri",
    )
    parser.add_argument(
        "--out-root",
        default="data",
        help="Output root directory",
    )

    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list-studies")
    sub.add_parser("list-studies-json")
    dl = sub.add_parser("download-study")
    dl.add_argument("study")
    dl.add_argument(
        "--aspera",
        action="store_true",
        help="Use Aspera Connect (ascp) to download per-run directories",
    )
    dl.add_argument("--ascp-path", default="", help="Path to ascp (or env ASCP_PATH)")
    dl.add_argument("--ascp-key", default="", help="Path to Aspera SSH key (or env ASCP_KEY)")
    dl.add_argument("--aspera-host", default="", help="e.g. aspera01@download.cncb.ac.cn")
    dl.add_argument("--aspera-remote-root", default="", help="e.g. gsa-human")
    dl.add_argument("--aspera-port", type=int, default=33001)
    dl.add_argument("--aspera-rate", default="100m")

    dr = sub.add_parser("download-run")
    dr.add_argument("study")
    dr.add_argument("run_id")
    dr.add_argument(
        "--aspera",
        action="store_true",
        help="Use Aspera Connect (ascp) to download a single run directory",
    )
    dr.add_argument("--ascp-path", default="", help="Path to ascp (or env ASCP_PATH)")
    dr.add_argument("--ascp-key", default="", help="Path to Aspera SSH key (or env ASCP_KEY)")
    dr.add_argument("--aspera-host", default="", help="e.g. aspera01@download.cncb.ac.cn")
    dr.add_argument("--aspera-remote-root", default="", help="e.g. gsa-human")
    dr.add_argument("--aspera-port", type=int, default=33001)
    dr.add_argument("--aspera-rate", default="100m")

    args = parser.parse_args(argv)
    metadata_csv = Path(os.path.expanduser(args.metadata))
    out_root = Path(os.path.expanduser(args.out_root))

    if args.cmd == "list-studies":
        cmd_list_studies(metadata_csv)
        return 0
    if args.cmd == "list-studies-json":
        cmd_list_studies_json(metadata_csv)
        return 0
    if args.cmd == "download-study":
        cmd_download_study(
            metadata_csv,
            out_root,
            args.study,
            use_aspera=args.aspera,
            ascp_path=_resolve_aspera_param(args.ascp_path, "ASCP_PATH"),
            key_path=_resolve_aspera_param(args.ascp_key, "ASCP_KEY"),
            aspera_host=(args.aspera_host or "").strip(),
            aspera_remote_root=(args.aspera_remote_root or "").strip(),
            aspera_port=args.aspera_port,
            aspera_rate_limit=args.aspera_rate,
        )
        return 0

    if args.cmd == "download-run":
        cmd_download_run(
            metadata_csv,
            out_root,
            args.study,
            args.run_id,
            use_aspera=args.aspera,
            ascp_path=_resolve_aspera_param(args.ascp_path, "ASCP_PATH"),
            key_path=_resolve_aspera_param(args.ascp_key, "ASCP_KEY"),
            aspera_host=(args.aspera_host or "").strip(),
            aspera_remote_root=(args.aspera_remote_root or "").strip(),
            aspera_port=args.aspera_port,
            aspera_rate_limit=args.aspera_rate,
        )
        return 0

    raise AssertionError("unreachable")


if __name__ == "__main__":
    raise SystemExit(main())
