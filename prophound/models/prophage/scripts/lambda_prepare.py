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
"""Build the LAMBDA clean-subset benchmark inputs (docs/claims-and-evidence.md, gate G1-G4).

Reads LAMBDA's per-genome dashboard JSONs (github.com/leannmlindsey/lambda-benchmark,
public/data/<assembly>_{2k,4k,8k}.json) and its reference table
(Zenodo LAMBDA_v1, metadata/prophage_reference_locations.csv). Writes, for the
clean subset (analysis/external_benchmarks/lambda_genome_overlap.tsv: genomes not
in our corpus, plus test-split genomes):

  lambda_genomes.tsv  assembly, the one sequence carrying its curated elements
                      (LAMBDA scores that sequence only: genome_length equals its
                      length rounded up to 1 kb, checked for all 39), corpus status
  lambda_truth.tsv    curated elements, 0-based half-open on that sequence
  lambda_calls.tsv    every other tool's published calls, same frame. Genomic
                      language models vary by segment size, so all three sizes
                      are kept (tool label "<model> <size>"); annotation tools
                      are identical across sizes and kept once.
  lambda_provenance.json  source versions and checksums

LAMBDA coordinates are 1-based inclusive (segments start at 1).

Example:
  python scripts/lambda_prepare.py --json_dir /tmp/lambda --reference LAMBDA_v1/metadata/prophage_reference_locations.csv \\
      --overlap analysis/external_benchmarks/lambda_genome_overlap.tsv --out_dir analysis/external_benchmarks/lambda
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

SIZES = ("2k", "4k", "8k")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json_dir", required=True, help="holds <assembly>_{2k,4k,8k}.json")
    ap.add_argument("--reference", required=True, help="LAMBDA metadata/prophage_reference_locations.csv")
    ap.add_argument("--overlap", required=True, help="lambda_genome_overlap.tsv (corpus status per assembly)")
    ap.add_argument("--source_note", default="", help="free text: repo commit / Zenodo record used")
    ap.add_argument("--out_dir", required=True)
    a = ap.parse_args()

    ov = pd.read_csv(a.overlap, sep="\t")
    clean = ov[ov.corpus_status.isin(["not_in_corpus", "test"])]
    ref = pd.read_csv(a.reference)
    genomes, truth, calls, sums, dup_truth, ran, cands = [], [], [], {}, {}, [], []
    for r in clean.itertuples():
        js = {}
        for s in SIZES:
            p = Path(a.json_dir) / f"{r.assembly}_{s}.json"
            js[s] = json.loads(p.read_text())
            sums[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
        j = js["8k"]
        seqs = ref[ref.Assembly == r.assembly]["NCBI Id"].unique()
        if len(seqs) != 1:
            raise SystemExit(f"{r.assembly}: {len(seqs)} truth sequences in the reference table")
        seq = seqs[0]
        genomes.append({"assembly": r.assembly, "sequence": seq, "organism": j["organism"],
                        "genome_length": j["genome_length"], "corpus_status": r.corpus_status,
                        "genus": j.get("taxonomy", {}).get("genus", "")})
        # the dashboard JSON repeats some elements (Photorhabdus TTO1: 25 entries, 20 unique);
        # the reference table is the deduplicated source, and the two must agree as sets
        gt = {(int(x["start"]), int(x["end"])) for x in j["ground_truth"]}
        rt = ref[ref.Assembly == r.assembly]
        rset = set(zip(rt.start.astype(int), rt.end.astype(int)))
        if gt != rset:
            raise SystemExit(f"{r.assembly}: JSON ground_truth disagrees with the reference table")
        dup_truth[r.assembly] = len(j["ground_truth"]) - len(gt)
        for s0, e0 in sorted(rset):
            truth.append({"assembly": r.assembly, "sequence": seq, "start": s0 - 1, "end": e0})
        for x in j.get("candidate_prophages", []):
            cands.append({"assembly": r.assembly, "sequence": seq, "start": int(x["start"]) - 1, "end": int(x["end"]),
                          "classification": x.get("classification", ""), "num_models": x.get("num_models", ""),
                          "novel_flag": x.get("novel_flag", "")})
        seen_annot = set()
        for s in SIZES:
            for tool in js[s].get("metrics", {}):   # a tool listed in metrics was run on this genome
                if tool.endswith(f" {s}") or (tool, "any") not in {(t["tool"], "any") for t in ran if t["assembly"] == r.assembly}:
                    ran.append({"assembly": r.assembly, "tool": tool})
            for tool, lst in js[s]["clustered_predictions"].items():
                is_lm = tool.endswith(f" {s}")
                if not is_lm and tool in seen_annot:
                    continue
                seen_annot.add(tool)
                for x in lst:
                    calls.append({"assembly": r.assembly, "sequence": seq, "tool": tool,
                                  "start": int(x["start"]) - 1, "end": int(x["end"])})
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    g, t, c = pd.DataFrame(genomes), pd.DataFrame(truth), pd.DataFrame(calls)
    g.to_csv(out / "lambda_genomes.tsv", sep="\t", index=False)
    pd.DataFrame(ran).drop_duplicates().to_csv(out / "lambda_tools_run.tsv", sep="\t", index=False)
    pd.DataFrame(cands).to_csv(out / "lambda_candidates.tsv", sep="\t", index=False)
    # 1-based inclusive copy in the eval-manifest format, for eval_containment.py
    t.assign(accession=t.assembly, contig=t.sequence, start=t.start + 1).assign(len_bp=lambda d: d.end - d.start + 1)[
        ["accession", "contig", "start", "end", "len_bp"]].to_csv(out / "lambda_truth_manifest.csv", index=False)
    t.to_csv(out / "lambda_truth.tsv", sep="\t", index=False)
    c.to_csv(out / "lambda_calls.tsv", sep="\t", index=False)
    (out / "lambda_provenance.json").write_text(json.dumps({
        "source_note": a.source_note,
        "duplicate_truth_entries_dropped": {k: v for k, v in dup_truth.items() if v}, "reference_sha256": hashlib.sha256(Path(a.reference).read_bytes()).hexdigest(),
        "json_sha256": sums, "coordinates": "0-based half-open on `sequence` (LAMBDA is 1-based inclusive)"}, indent=2) + "\n")
    print(f"[lambda] {len(g)} genomes ({g.corpus_status.value_counts().to_dict()}), {len(t)} curated elements, "
          f"{len(c)} tool calls from {c.tool.nunique()} tool labels", flush=True)


if __name__ == "__main__":
    main()
