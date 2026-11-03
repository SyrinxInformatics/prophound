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
"""Refuse commits and pushes that carry secrets or files that do not belong in the repository.

Standard library only, so it runs from a git hook, from CI or by hand:

  python3 scripts/check_commit.py --staged          # the files a commit would add (pre-commit hook)
  python3 scripts/check_commit.py --range A..B      # every file version in commits A..B (pre-push hook)
  python3 scripts/check_commit.py --all             # every tracked file (CI)

Checks, in the main repository:
  - secrets: private keys and the token formats of GitHub, AWS, Slack, OpenAI, Anthropic, Google,
    Hugging Face and GitLab, plus password/token/key assignments with a quoted literal;
  - credential files by name (.env, *.pem, *.key, id_rsa, .netrc, .pypirc, *.p12, ...);
  - files that stay local: CLAUDE.md, SESSIONS.md, docs/release-checklist.md and .git-local/;
  - files over 10 MB;
  - sequence data and model checkpoints outside the places that ship them;
  - OS, editor and cache files, and local data directories (arc_data/, .pixi/, runs/).

In the local-only repository (.git-local), only Markdown files are allowed, and the secret check applies.

A line containing `check-commit: allow` is exempt from the secret patterns (for documented false
positives). Install the hooks once per clone with `pixi run install-hooks`.
"""

from __future__ import annotations

import argparse
import fnmatch
import os
import re
import subprocess
import sys

MAX_BYTES = 10_000_000

SECRET = re.compile(
    rb"BEGIN [A-Z ]*PRIVATE KEY"
    rb"|ghp_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}|gh[ousr]_[A-Za-z0-9]{30,}"
    rb"|AKIA[0-9A-Z]{16}|ASIA[0-9A-Z]{16}"
    rb"|xox[baprs]-[A-Za-z0-9-]{10,}"
    rb"|sk-ant-[A-Za-z0-9_-]{20,}|sk-[A-Za-z0-9]{32,}"
    rb"|AIza[0-9A-Za-z_-]{35}"
    rb"|hf_[A-Za-z0-9]{30,}"
    rb"|glpat-[A-Za-z0-9_-]{20}"
    rb"|(?i:(?:password|passwd|passphrase|secret|token|api_?key|access_?key)\s*[:=]\s*[\"'][^\"'\s]{6,}[\"'])"
)
ALLOW_MARK = b"check-commit: allow"

CREDENTIAL_NAMES = [
    ".env", ".env.*", "*.pem", "*.key", "*.p12", "*.pfx", "*.keystore", "*.jks",
    "id_rsa*", "id_dsa*", "id_ecdsa*", "id_ed25519*", ".netrc", "_netrc", ".pypirc",
    ".git-credentials", "credentials", "credentials.*", "*credentials.json", "secrets.*", ".htpasswd",
]
JUNK_NAMES = [".DS_Store", "Thumbs.db", "*.pyc", "*.pyo", "*.swp", "*.swo", "*~", ".Rhistory", ".RData"]
JUNK_DIRS = ["__pycache__", ".ipynb_checkpoints", ".pytest_cache", ".ruff_cache", ".mypy_cache", ".idea", ".vscode"]
LOCAL_ONLY = ["CLAUDE.md", "SESSIONS.md", "docs/release-checklist.md"]
LOCAL_DATA_DIRS = ["arc_data", ".pixi", ".git-local"]

SEQUENCE = re.compile(
    r"\.(fna|fa|fasta|fas|ffn|faa|gb|gbk|gbff|embl|fastq|fq|sam|bam|cram|2bit)(\.gz|\.bz2|\.xz)?$", re.I)
# the package moved from src/mge_detect/ (renamed 2026-09-29) and src/prophound/ (2026-11-02); history keeps them
PACKAGE_PREFIX = r"(?:src/(?:prophound|mge_detect)/|prophound/)"
SEQUENCE_OK = re.compile(rf"^(?:examples/|{PACKAGE_PREFIX}models/prophage/scripts/spikein_small/)")
CHECKPOINT = re.compile(r"\.(pt|pth|ckpt|safetensors|h5|onnx)$", re.I)
CHECKPOINT_OK = re.compile(rf"^{PACKAGE_PREFIX}models/prophage/model/")


def git(*args: str) -> bytes:
    return subprocess.run(["git", *args], capture_output=True, check=True).stdout


def is_local_repo() -> bool:
    return os.path.basename(os.path.normpath(git("rev-parse", "--git-dir").decode().strip())) == ".git-local"


def path_problems(path: str, local_repo: bool) -> list[str]:
    name = os.path.basename(path)
    parts = path.split("/")
    out = []
    if any(fnmatch.fnmatch(name, pat) for pat in CREDENTIAL_NAMES):
        out.append("credential file name")
    if local_repo:
        if not path.endswith(".md"):
            out.append("the local-only repository holds Markdown files only")
        return out
    if path in LOCAL_ONLY:
        out.append("local-only file (tracked in .git-local)")
    if any(fnmatch.fnmatch(name, pat) for pat in JUNK_NAMES) or any(p in JUNK_DIRS for p in parts[:-1]):
        out.append("OS, editor or cache file")
    if parts[0] in LOCAL_DATA_DIRS or "runs" in parts[:-1]:
        out.append("local data directory")
    if SEQUENCE.search(path) and not SEQUENCE_OK.match(path):
        out.append("sequence data outside examples/ and the spike-in set")
    if CHECKPOINT.search(path) and not CHECKPOINT_OK.match(path):
        out.append("model checkpoint outside the bundled model/ directory")
    return out


def content_problems(data: bytes) -> list[str]:
    out = []
    if len(data) > MAX_BYTES:
        out.append(f"{len(data) / 1e6:.1f} MB, over the {MAX_BYTES / 1e6:.0f} MB limit")
    for line_no, line in enumerate(data.splitlines(), 1):
        if ALLOW_MARK in line:
            continue
        m = SECRET.search(line)
        if m:
            shown = m.group(0)[:12].decode(errors="replace")
            out.append(f"possible secret on line {line_no} ({shown}...)")
    return out


def blobs_staged() -> list[tuple[str, str]]:
    """(path, blob id) for every file a commit would add or change."""
    names = git("diff", "--cached", "--name-only", "--diff-filter=ACMR", "-z").split(b"\0")
    out = []
    for raw in filter(None, names):
        path = raw.decode()
        sha = git("rev-parse", f":{path}").decode().strip()
        out.append((path, sha))
    return out


def blobs_range(rev_range: str) -> list[tuple[str, str]]:
    """(path, blob id) for every file version added or changed by a commit in the range."""
    out, seen = [], set()
    for commit in git("rev-list", *rev_range.split()).decode().split():
        listing = git("diff-tree", "-r", "-z", "--no-commit-id", "--root", "--diff-filter=ACMR", commit)
        fields = listing.split(b"\0")
        for meta, raw in zip(fields[0::2], fields[1::2]):
            if not meta:
                continue
            sha = meta.split()[3].decode()
            key = (raw.decode(), sha)
            if key not in seen:
                seen.add(key)
                out.append(key)
    return out


def blobs_tracked() -> list[tuple[str, str]]:
    out = []
    for entry in filter(None, git("ls-files", "-s", "-z").split(b"\0")):
        meta, path = entry.split(b"\t", 1)
        out.append((path.decode(), meta.split()[1].decode()))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--staged", action="store_true", help="check the files staged for commit")
    mode.add_argument("--range", metavar="A..B", help="check every file version in a commit range")
    mode.add_argument("--all", action="store_true", help="check every tracked file")
    args = ap.parse_args()

    local_repo = is_local_repo()
    if args.staged:
        blobs = blobs_staged()
    elif args.range:
        blobs = blobs_range(args.range)
    else:
        blobs = blobs_tracked()

    failures = []
    for path, sha in blobs:
        problems = path_problems(path, local_repo)
        if git("cat-file", "-t", sha).strip() == b"blob":
            problems += content_problems(git("cat-file", "blob", sha))
        failures += [f"  {path}: {p}" for p in problems]

    where = ".git-local" if local_repo else "main repository"
    if failures:
        print(f"[check-commit] refused ({where}, {len(blobs)} files checked):", file=sys.stderr)
        print("\n".join(failures), file=sys.stderr)
        print("[check-commit] fix or unstage these files. For a documented false positive, add "
              "'check-commit: allow' on that line.", file=sys.stderr)
        return 1
    print(f"[check-commit] ok ({where}, {len(blobs)} files checked)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
