from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from tests.support import build_validation_project

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "0check_data.sh"


def _git_bash() -> Path | None:
    candidates = [
        Path(os.environ.get("GIT_BASH", "")),
        Path("D:/Git/bin/bash.exe"),
        Path("C:/Program Files/Git/bin/bash.exe"),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    bash = shutil.which("bash")
    return Path(bash) if bash and "System32" not in bash else None


def _run_script(project: Path, *extra: str) -> subprocess.CompletedProcess[str]:
    bash = _git_bash()
    if bash is None:
        pytest.skip("Git Bash is required for 0check_data.sh tests")
    return subprocess.run(
        [str(bash), str(SCRIPT), "--project-root", str(project), *extra],
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )


def test_shell_wrapper_passes_valid_manifest_without_compatibility_links(
    tmp_path: Path,
) -> None:
    paths = build_validation_project(tmp_path, compatible_links=False)
    result = _run_script(paths["workspace"], "--no-compat-links")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "1 PASS / 0 FAIL" in result.stdout


def test_shell_wrapper_fails_empty_manifest(tmp_path: Path) -> None:
    paths = build_validation_project(tmp_path, species_count=0)
    result = _run_script(paths["workspace"])
    assert result.returncode == 1
    assert "empty_manifest" in result.stdout + result.stderr

