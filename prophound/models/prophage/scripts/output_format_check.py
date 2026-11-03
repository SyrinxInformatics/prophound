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
"""Check a declared output-format change (UI pass step 2): same calls, new presentation.

Compares two `reference_outputs.py run` directories, the old format (`a`) and the new
(`b`). The calls must be identical in coordinates and every value, probabilities within
4 dp rounding. windows.tsv must be identical apart from the dropped fasta_path. BED, GFF3
and FASTA must carry the same calls under call_id. The metadata may differ only in the
declared keys. regions.filtered.tsv must hold exactly the calls the metadata counts as set
aside (docs/output-files.md).

  python output_format_check.py <old run dir> <new run dir> <inputs dir>
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

old_root, new_root, inputs = map(Path, sys.argv[1:4])
SHORT = re.compile(r"_short(\d+)")
RT = {"sep": "\t", "float_precision": "round_trip"}
CHANGED_META = {"checkpoint", "checkpoint_sha256", "input_fasta", "outputs", "num_tiles"}
NEW_META = {"config", "input_sha256", "started_utc", "elapsed_s", "command", "platform"}


def contig_order(fna: Path) -> dict:
    ids = [ln[1:].split()[0] for ln in fna.read_text().splitlines() if ln.startswith(">")]
    return {c: i for i, c in enumerate(ids)}


def fasta_records(p: Path) -> list:
    recs, cur = [], None
    for ln in p.read_text().splitlines():
        if ln.startswith(">"):
            cur = [ln[1:].split()[1], []]
            recs.append(cur)
        else:
            cur[1].append(ln)
    return sorted((loc, "".join(seq)) for loc, seq in recs)


manifest = json.loads((new_root / "manifest.json").read_text())["cases"]
old_manifest = json.loads((old_root / "manifest.json").read_text())["cases"]
bad = 0
for case, info in sorted(manifest.items()):
    diffs = []
    if info["exit_code"] != old_manifest[case]["exit_code"]:
        diffs.append("exit code")
    if info["exit_code"] != 0:
        print(f"{case}: exit {info['exit_code']} (both)" if not diffs else f"{case}: {diffs}")
        bad += bool(diffs)
        continue
    o, n = old_root / case, new_root / case
    order = contig_order(inputs / f"{info['input']}.fna")
    meta_o = json.loads((o / "run_metadata.json").read_text())
    meta_n = json.loads((n / "run_metadata.json").read_text())
    wl = meta_o["window_len"]

    ro, rn = pd.read_csv(o / "regions.tsv", **RT), pd.read_csv(n / "regions.tsv", **RT)
    if len(ro):
        ro["window_bp"] = [int(m.group(1)) if (m := SHORT.search(k)) else wl for k in ro.boundary_kind.astype(str)]
        ro = ro.drop(columns="boundary_kind")
        ro = ro.assign(_o=ro.contig_id.map(order)).sort_values(["_o", "start", "end"], kind="stable")
        ro = ro.drop(columns="_o").reset_index(drop=True)
    if list(rn.call_id) != [f"prophound_{i}" for i in range(1, len(rn) + 1)]:
        diffs.append("call_id not prophound_1..N")
    rn = rn.drop(columns="call_id")
    if sorted(ro.columns) != sorted(rn.columns) or len(ro) != len(rn):
        diffs.append(f"regions shape/columns {list(ro.columns)} vs {list(rn.columns)}")
    else:
        for c in ro.columns:
            if c in ("mean_prob", "max_prob"):
                if len(ro) and (ro[c] - rn[c]).abs().max() > 5.0001e-5:
                    diffs.append(f"regions {c} beyond rounding")
            elif not ro[c].equals(rn[c]):
                diffs.append(f"regions {c}")
    rn_full = pd.read_csv(n / "regions.tsv", **RT)

    wo, wn = pd.read_csv(o / "windows.tsv", **RT).drop(columns="fasta_path"), pd.read_csv(n / "windows.tsv", **RT)
    if not wo.equals(wn):
        diffs.append("windows.tsv")
    if (o / "fine_windows.tsv").exists():
        fo = pd.read_csv(o / "fine_windows.tsv", **RT).drop(columns="fasta_path", errors="ignore")
        if not fo.equals(pd.read_csv(n / "fine_windows.tsv", **RT)):
            diffs.append("fine_windows.tsv")

    bed = [ln.split("\t") for ln in (n / "regions.bed").read_text().splitlines()]
    exp = [[r.contig_id, str(r.start), str(r.end), r.call_id, str(int(round(1000 * r.mean_prob))), "."]
           for r in rn_full.itertuples()]
    if bed != exp:
        diffs.append("regions.bed")
    gff = [ln.split("\t") for ln in (n / "regions.gff3").read_text().splitlines() if not ln.startswith("#")]
    if [(g[0], int(g[3]) - 1, int(g[4]), g[8].split(";")[0]) for g in gff] != \
            [(r.contig_id, r.start, r.end, f"ID={r.call_id}") for r in rn_full.itertuples()]:
        diffs.append("regions.gff3")
    if (o / "regions.fna").exists() and fasta_records(o / "regions.fna") != fasta_records(n / "regions.fna"):
        diffs.append("regions.fna sequences")

    for k in sorted((set(meta_o) | set(meta_n)) - CHANGED_META - NEW_META):
        if meta_o.get(k) != meta_n.get(k):
            diffs.append(f"metadata[{k}]")
    if not set(meta_n) >= NEW_META:
        diffs.append(f"metadata lacks {NEW_META - set(meta_n)}")
    if sorted(meta_n["checkpoint_sha256"].values()) != sorted(meta_o["checkpoint_sha256"].values()):
        diffs.append("checkpoint hashes")
    if "/Users/" in json.dumps(meta_n):
        diffs.append("home path in metadata")
    if set(meta_n["outputs"]) != {p.name for p in n.iterdir() if p.name not in ("stderr.txt",)}:
        diffs.append(f"outputs list {meta_n['outputs']}")

    filt = pd.read_csv(n / "regions.filtered.tsv", **RT)
    expected = ((meta_n.get("seed_support") or {}).get("calls_dropped_by_vote", 0)
                + meta_n["n_gap_filter"]["calls_dropped"])
    columns_ok = list(filt.columns) == [*rn_full.columns, "filter"]
    if len(filt) != expected or not columns_ok:
        diffs.append(f"filtered: {len(filt)} rows vs {expected} set aside, columns ok={columns_ok}")
    bad += bool(diffs)
    print(f"{case}: {'same calls' if not diffs else 'DIFF ' + '; '.join(diffs)} "
          f"({len(rn)} calls, {len(filt)} filtered)")
print("all same" if not bad else f"{bad} cases differ")
