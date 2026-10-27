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
"""Are att sites real? An independent check: do the att cores sit on tRNA/tmRNA genes?

docs/att-boundaries.md "Are the att sites real?". Many temperate phages and other
integrative elements integrate into a tRNA or tmRNA gene: the core is the gene's 3' end, so
one copy of it overlaps the gene. The att search never sees genes, so this signature is
independent of how att pairs are chosen (analysis only; nothing here enters the pipeline).

For each sample, three groups of accepted att pairs, all under the shipped rule (src/att.py):
  - `call`: att_boundaries.py's call items, selected offline with the shipped thresholds;
  - `decoy`: the same search and rule on host zones of the same geometry (the false-att
    control att_boundaries.py already uses);
  - `shipped`: the shipped pipeline's att calls (`<sample>_pipeline.calls.tsv`,
    boundary = att; attL = [start, start + att_len), attR = [end - att_len, end)).
A pair counts as tRNA-borne when either core copy overlaps an ARAGORN tRNA or tmRNA gene.
Each group is also split by whether its call overlaps a truth element, and compared with the
rate expected if the same cores sat at random positions on the same contig.

ARAGORN 1.2.41 runs through `pixi exec` (bioconda), so the project environment is unchanged;
its output is cached in --out/aragorn/.

  python scripts/att_trna_check.py \\
      --sample tune42 analysis/att/tune42 ../../../../arc_data/raw/genomes scripts/eval_sizestrat_elements.csv \\
      --out analysis/att/trna
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import att as att_mod                      # noqa: E402
from att_boundaries import select          # noqa: E402

GENE = re.compile(r"^\s*\d+\s+(tRNA-\S+|tmRNA\S*)\s+c?\[(\d+),(\d+)\]")


def aragorn(fna: Path, cache: Path) -> pd.DataFrame:
    if not cache.exists():
        out = subprocess.run(["pixi", "exec", "-c", "bioconda", "-c", "conda-forge", "-s", "aragorn",
                              "aragorn", "-t", "-m", "-gcbact", "-l", "-w", str(fna)],
                             check=True, capture_output=True, text=True, cwd=cache.parent).stdout
        cache.write_text(out)
    rows, contig = [], None
    for line in cache.read_text().splitlines():
        if line.startswith(">"):
            contig = line[1:].split()[0]
            continue
        m = GENE.match(line)
        if m and contig:
            a, b = int(m.group(2)), int(m.group(3))
            if b >= a:                                     # skip genes wrapped across the origin
                rows.append((contig, a - 1, b, m.group(1)))
    return pd.DataFrame(rows, columns=["contig_id", "start", "end", "gene"])


def contig_lengths(fna: Path) -> dict[str, int]:
    fai = Path(str(fna) + ".fai")
    if fai.exists():
        return {l.split("\t")[0]: int(l.split("\t")[1]) for l in fai.read_text().splitlines()}
    out, name, n = {}, None, 0
    for line in open(fna):
        if line.startswith(">"):
            if name:
                out[name] = n
            name, n = line[1:].split()[0], 0
        else:
            n += len(line.strip())
    if name:
        out[name] = n
    return out


def hits(contig_genes: pd.DataFrame, lo: int, hi: int) -> str:
    g = contig_genes[(contig_genes.end > lo) & (contig_genes.start < hi)]
    return ";".join(g.gene) if len(g) else ""


def null_rate(contig_genes: pd.DataFrame, core: int, L: int) -> float:
    """P(at least one of two random copies of a `core`-bp core overlaps a gene on this contig)."""
    p = min(1.0, float((contig_genes.end - contig_genes.start + core - 1).sum()) / max(L, 1))
    return 1 - (1 - p) ** 2


def on_element(truth: pd.DataFrame, genome: str, contig: str, s: int, e: int) -> bool:
    t = truth[(truth.accession == genome) & (truth.contig.astype(str).str.split(".").str[0]
                                             == str(contig).split(".")[0])]
    return bool(((t.end > s) & (t.start < e)).any())


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", nargs=4, action="append", required=True,
                    metavar=("NAME", "PREFIX", "GENOMES_DIR", "ELEMENTS"),
                    help="att_boundaries.py run prefix (<prefix>.items.tsv, .cands.tsv and "
                         "<prefix>_pipeline.calls.tsv; items optional), its genomes and truth CSV")
    ap.add_argument("--modes", nargs="+", default=["default", "exploratory", "conservative"])
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    (a.out / "aragorn").mkdir(parents=True, exist_ok=True)

    pairs = []
    for name, prefix, gdir, elements in a.sample:
        truth = pd.read_csv(elements)
        frames = []
        if Path(f"{prefix}.items.tsv").exists():
            items = pd.read_csv(f"{prefix}.items.tsv", sep="\t")
            cands = pd.read_csv(f"{prefix}.cands.tsv", sep="\t")
            sel = select(items, cands, att_mod.ACCEPT_IN_BP, att_mod.ACCEPT_OUT_BP, att_mod.MAX_CORE_LEN,
                         att_mod.MAX_COPIES, att_mod.MIN_ENTROPY_BITS, None, att_mod.MAX_EVALUE)
            sel = sel.merge(items[["item_id", "mode", "genome_id", "contig_id", "kind", "call_start", "call_end"]],
                            on="item_id")
            frames.append(sel.rename(columns={"kind": "group"}).assign(
                l0=lambda d: d.att_start, l1=lambda d: d.att_start + d.att_len,
                r0=lambda d: d.att_end - d.att_len, r1=lambda d: d.att_end))
        pc = pd.read_csv(f"{prefix}_pipeline.calls.tsv", sep="\t")
        pc = pc[pc.boundary == "att"].copy()
        pc["att_len"] = pc.att_len.astype(int)
        frames.append(pc.assign(group="shipped", call_start=pc.start, call_end=pc.end,
                                l0=pc.start, l1=pc.start + pc.att_len, r0=pc.end - pc.att_len, r1=pc.end)
                      [["mode", "genome_id", "contig_id", "group", "call_start", "call_end", "att_len",
                        "l0", "l1", "r0", "r1"]])
        df = pd.concat(frames, ignore_index=True)
        df = df[df["mode"].isin(a.modes)]
        for genome, g in df.groupby("genome_id"):
            fna = next(Path(gdir).glob(f"{genome}*.fna"), None) or Path(gdir) / f"{genome}.fna"
            genes = aragorn(fna.resolve(), (a.out / "aragorn" / f"{genome}.txt").resolve())
            lengths = contig_lengths(fna)
            for _, r in g.iterrows():
                cg = genes[genes.contig_id == r.contig_id]
                hl, hr = hits(cg, int(r.l0), int(r.l1)), hits(cg, int(r.r0), int(r.r1))
                pairs.append(dict(sample=name, mode=r["mode"], group=r.group, genome_id=genome,
                                  contig_id=r.contig_id, call_start=int(r.call_start), call_end=int(r.call_end),
                                  att_len=int(r.att_len), attL=f"{int(r.l0)}-{int(r.l1)}",
                                  attR=f"{int(r.r0)}-{int(r.r1)}", trna_left=hl, trna_right=hr,
                                  trna=bool(hl or hr),
                                  null=null_rate(cg, int(r.att_len), lengths.get(r.contig_id, 0)),
                                  on_element=(r.group != "decoy") and on_element(
                                      truth, genome, r.contig_id, int(r.call_start), int(r.call_end))))
    p = pd.DataFrame(pairs)
    p.to_csv(a.out / "att_trna.pairs.tsv", sep="\t", index=False)

    rows = []
    p["where"] = np.where(p.group == "decoy", "host", np.where(p.on_element, "on element", "off element"))
    for keys, g in list(p.groupby(["mode", "group", "where"])) + list(p.groupby(["mode", "group"])):
        k, n = int(g.trna.sum()), len(g)
        lo, hi = wilson(k, n)
        rows.append(dict(mode=keys[0], group=keys[1], where=keys[2] if len(keys) > 2 else "all",
                         n=n, trna=k, rate=k / n, lo=lo, hi=hi, null=g.null.mean(),
                         samples=",".join(sorted(g["sample"].unique()))))
    s = pd.DataFrame(rows).sort_values(["mode", "group", "where"])
    s.to_csv(a.out / "att_trna.summary.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 200, "display.float_format", "{:.3f}".format):
        print(s.to_string(index=False))


if __name__ == "__main__":
    main()
