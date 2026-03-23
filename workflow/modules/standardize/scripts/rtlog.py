from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path
from typing import Iterable, TextIO


def _ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def _write_line(log: TextIO, line: str) -> None:
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

    _write_line(
        log,
        f"[{_ts()}] run: {' '.join(cmd_list)}" + (f" (cwd={cwd})" if cwd else ""),
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
            _write_line(log, f"{prefix}{line.rstrip()}\n")
        else:
            _write_line(log, line.rstrip("\n"))
    proc.stdout.close()
    return_code = proc.wait()
    _write_line(log, f"[{_ts()}] exit={return_code}: {' '.join(cmd_list)}")
    return return_code


def tee_lines(lines: Iterable[str], *, log: TextIO) -> None:
    for line in lines:
        _write_line(log, line)
