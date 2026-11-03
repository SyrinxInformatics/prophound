# SPDX-License-Identifier: Apache-2.0
#
# Copyright 2025 Jacob Alford <jalford0000@gmail.com>
#                and Benjamin Alford <ssolar5.dev@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""The commit guard (scripts/check_commit.py) refuses what it should and nothing else."""

import importlib.util
from pathlib import Path

import pytest

_spec = importlib.util.spec_from_file_location(
    "check_commit", Path(__file__).resolve().parents[1] / "scripts" / "check_commit.py")
check_commit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(check_commit)

# assembled at run time so this file itself never holds a token-shaped string
GITHUB = "ghp_" + "a1" * 18
AWS = "AKIA" + "ABCDEFGHIJKLMNOP"


@pytest.mark.parametrize("text", [
    f"TOKEN={GITHUB}",
    f"aws = {AWS}",
    "-----BEGIN RSA " + "PRIVATE KEY-----",
    'password = "' + "hunter22secret" + '"',
    "api_key: '" + "abcdef123456" + "'",
])
def test_secrets_are_refused(text):
    assert any("possible secret" in p for p in check_commit.content_problems(text.encode()))


@pytest.mark.parametrize("text", [
    "x = 1",
    "password = os.environ['PASSWORD']",
    'token = "' + "notarealtoken1" + '"  # check-commit: allow',
    "tokens = tokenize(text)",
])
def test_ordinary_text_passes(text):
    assert check_commit.content_problems(text.encode()) == []


def test_size_limit():
    assert check_commit.content_problems(b"\0" * (check_commit.MAX_BYTES + 1))
    assert not check_commit.content_problems(b"\0" * 1000)


@pytest.mark.parametrize("path", [
    ".env", "config/.env.local", "deploy/server.pem", "id_ed25519", ".netrc", "aws/credentials",
    "CLAUDE.md", "SESSIONS.md", "docs/release-checklist.md",
    ".DS_Store", "prophound/__pycache__/cli.cpython-314.pyc", ".Rhistory",
    "arc_data/raw/genome.fna", ".pixi/envs/x", "prophound/models/prophage/runs/r1/best.pt",
    "genome.fna", "data/reads.fastq.gz", "prophound/models/prophage/analysis/x/best.pt",
])
def test_paths_are_refused(path):
    assert check_commit.path_problems(path, local_repo=False)


@pytest.mark.parametrize("path", [
    "README.md", "prophound/cli.py", "docs/benchmarks.md",
    "examples/hin_rd_kw20_region.fna",
    "prophound/models/prophage/model/prophage_v2_pool_r1.pt",
    "src/prophound/models/prophage/model/prophage_v2_pool_r1.pt",   # the pre-2026-11-02 layout, in history
    "src/mge_detect/models/prophage/scripts/spikein_small/candidates.fna",
    "prophound/models/prophage/analysis/lockbox2/results.json",
])
def test_shipped_paths_pass(path):
    assert check_commit.path_problems(path, local_repo=False) == []


def test_local_repo_holds_markdown_and_text_only():
    assert check_commit.path_problems("CLAUDE.md", local_repo=True) == []
    assert check_commit.path_problems("docs/release-checklist.md", local_repo=True) == []
    assert check_commit.path_problems("prophound/models/prophage/notes/RESEARCH_NOTES.txt", local_repo=True) == []
    assert check_commit.path_problems("scripts/x.py", local_repo=True)
    assert check_commit.path_problems(".env", local_repo=True)


def test_files_the_local_repo_tracks_are_refused_under_any_layout():
    local_only = frozenset(check_commit.LOCAL_ONLY + ["docs/private.md", "prophound/models/prophage/notes/N.md"])
    for path in ["docs/private.md", "prophound/models/prophage/notes/N.md",
                 "src/prophound/models/prophage/notes/N.md", "src/mge_detect/models/prophage/notes/N.md"]:
        assert check_commit.path_problems(path, local_repo=False, local_only=local_only), path
    assert check_commit.path_problems("docs/benchmarks.md", local_repo=False, local_only=local_only) == []


def test_local_repo_listing_ignores_the_hook_environment(tmp_path, monkeypatch):
    """Inside `git commit`, GIT_INDEX_FILE names the main index; the .git-local listing must not use it."""
    import subprocess

    def run(*args, **kw):
        subprocess.run(["git", *args], cwd=tmp_path, check=True, capture_output=True, **kw)

    run("init", "-q")
    (tmp_path / "public.md").write_text("x\n")
    (tmp_path / "private.md").write_text("y\n")
    run("add", "public.md")
    run("--git-dir=.git-local", f"--work-tree={tmp_path}", "init", "-q")
    run("--git-dir=.git-local", f"--work-tree={tmp_path}", "add", "private.md")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("GIT_INDEX_FILE", str(tmp_path / ".git" / "index"))
    monkeypatch.setenv("GIT_DIR", str(tmp_path / ".git"))
    tracked = check_commit.local_only_paths()
    assert "private.md" in tracked
    assert "public.md" not in tracked
