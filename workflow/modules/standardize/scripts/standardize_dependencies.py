from __future__ import annotations

from pathlib import PurePosixPath
from typing import Iterable


def unpack_requires_fetch_done(*, study: str, required_paths: Iterable[str]) -> bool:
    """Return True when standardize unpack should wait on fetch markers.

    Only studies whose required inputs live under out/{study}/raw or out/{study}/qconly
    need fetch/unpack coordination. Studies that read directly from other ready paths such
    as out/{study}/counts/... should not be blocked on missing fetch marker files.
    """

    for raw_path in required_paths:
        parts = PurePosixPath(str(raw_path).replace('\\', '/')).parts
        for idx in range(len(parts) - 2):
            if parts[idx] == 'out' and parts[idx + 1] == study and parts[idx + 2] in {'raw', 'qconly'}:
                return True
    return False
