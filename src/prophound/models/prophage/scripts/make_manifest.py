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
import csv
from pathlib import Path

allc = "windows/bed/all_win.contig.sorted.bed"
lab = "windows/labels/all_512k.labels.tsv"
tr = "splits/train.contig.bed"
va = "splits/val.contig.bed"
te = "splits/test.contig.bed"
outp = "windows/manifest/all_512k_manifest.csv"


def idset(path):
    s = set()
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            s.add(line.rstrip("\n").split("\t")[3])
    return s


split_map = {}
for i, (p, tag) in enumerate([(tr, "train"), (va, "val"), (te, "test")]):
    for k in idset(p):
        split_map[k] = tag

labels = {}
with open(lab) as f:
    for line in f:
        if not line.strip():
            continue
        _id, y = line.rstrip("\n").split("\t")
        labels[_id] = int(y)

rows = []
with open(allc) as f:
    for line in f:
        if not line.strip():
            continue
        contig, start, end, _id = line.rstrip("\n").split("\t")
        rows.append(
            {
                "id": _id,
                "contig": contig,
                "start": int(start),
                "end": int(end),
                "label": labels.get(_id, ""),
                "split": split_map.get(_id, ""),
            }
        )

Path("windows/manifest").mkdir(parents=True, exist_ok=True)
with open(outp, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["id", "contig", "start", "end", "label", "split"])
    w.writeheader()
    w.writerows(rows)
print(outp, len(rows))

