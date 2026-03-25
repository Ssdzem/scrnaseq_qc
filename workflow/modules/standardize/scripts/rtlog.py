from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable, TextIO


def _ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log_line(log: TextIO, line: str, *, prefix: str | None = None, with_ts: bool = True) -> None:
    msg = line.rstrip("\n")
    if prefix:
        msg = f"{prefix}{msg}"
    if with_ts:
        msg = f"[{_ts()}] {msg}"
    log.write(msg + "\n")
    log.flush()
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


def log_kv(log: TextIO, event: str, /, **fields: object) -> None:
    if fields:
        ordered = " ".join(f"{k}={fields[k]}" for k in sorted(fields))
        log_line(log, f"{event} {ordered}")
        return
    log_line(log, event)


def _write_raw(log: TextIO, line: str) -> None:
    log.write(line)
    if not line.endswith("\n"):
        log.write("\n")
    log.flush()


def run_stream(
    cmd: list[str] | tuple[str, ...],
    *,
    log: TextIO,
    cwd: str | Path | None = None,
    env: dict[str, str] | None = None,
    prefix: str | None = None,
) -> int:
    """Run a command and stream stdout/stderr to log in real time.

    This avoids `subprocess.run(..., capture_output=True)` buffering, which can
    look like a hang on long-running steps.
    """

    cmd_list = list(cmd)
    merged_env = os.environ.copy()
    if env:
        merged_env.update({k: str(v) for k, v in env.items()})

    log_line(
        log,
        f"[{_ts()}] run: {' '.join(cmd_list)}" + (f" (cwd={cwd})" if cwd else ""),
        with_ts=False,
    )

    proc = subprocess.Popen(
        cmd_list,
        cwd=str(cwd) if cwd is not None else None,
        env=merged_env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
    )
    assert proc.stdout is not None
    for line in iter(proc.stdout.readline, ""):
        if prefix:
            log_line(log, line.rstrip("\n"), prefix=prefix)
        else:
            log_line(log, line.rstrip("\n"))
    proc.stdout.close()
    return_code = proc.wait()
    log_line(log, f"[{_ts()}] exit={return_code}: {' '.join(cmd_list)}", with_ts=False)
    return return_code


def tee_lines(lines: Iterable[str], *, log: TextIO) -> None:
    for line in lines:
        _write_raw(log, line)
