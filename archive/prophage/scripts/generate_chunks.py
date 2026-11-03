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

from pathlib import Path
import math
import sys

ROOT = Path("/work/pke835/mge_explorer/plasmids/selected_fna")  # root of all .fna
ROOT = ROOT.resolve()
CHUNK_SIZE = 1000
CHUNK_DIR = ROOT / "chunks"
CHUNK_DIR.mkdir(parents=True, exist_ok=True)

# Collect all files recursively as absolute paths.
files = []
for p in ROOT.rglob("*"):
    if p.is_file():
        files.append(p.resolve())

if not files:
    print(f"[ERROR] No files found under {ROOT}", file=sys.stderr)
    sys.exit(1)

paths = [str(p) for p in sorted(files)]

# Sanity-check that every path stays under ROOT.
bad = [p for p in paths if not p.startswith(str(ROOT))]
if bad:
    print("[ERROR] Found paths that are not under ROOT. Example:", bad[:5], file=sys.stderr)
    sys.exit(1)

# Write the master file list.
master_file = CHUNK_DIR / "files.txt"
with master_file.open("w", encoding="utf-8") as f:
    for path in paths:
        f.write(path + "\n")

print(f"[OK] Wrote {len(paths)} absolute paths to {master_file}")

print("[SAMPLE] first 5 paths:")
for p in paths[:5]:
    print("   ", p)

# Split the master list into chunks.
num_chunks = math.ceil(len(paths) / CHUNK_SIZE)
for i in range(num_chunks):
    chunk_paths = paths[i * CHUNK_SIZE : (i + 1) * CHUNK_SIZE]
    chunk_file = CHUNK_DIR / f"files_chunk_{i+1:04d}.list"
    # Guard against writing bare filenames instead of full paths.
    for cp in chunk_paths:
        if "/" not in cp:
            print(f"[ERROR] would write non-path entry: {cp}", file=sys.stderr)
            sys.exit(1)
    with chunk_file.open("w", encoding="utf-8") as f:
        for path in chunk_paths:
            f.write(path + "\n")
    print(f"  -> {chunk_file} ({len(chunk_paths)} lines)")

print(f"[DONE] Created {num_chunks} chunk files in {CHUNK_DIR}")
