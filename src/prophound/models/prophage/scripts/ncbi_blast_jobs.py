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
"""Submit remote BLAST searches to NCBI and fetch them later, from any session.

Uses NCBI's BLAST URL API (https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html) instead
of `blastn -remote`, which keeps the request ID only in memory. The jobs file is the only
state: one row per search, with its request ID (RID), the query file's SHA-256, the search's
settings and the submission time. `submit` writes a row as soon as NCBI accepts a search.
`fetch` reads only the jobs file, polls each RID at most once per minute (NCBI's guideline),
and saves the tabular result next to the jobs file. NCBI keeps results for about 36 hours.
After that, `fetch` reports the RID as expired, and `submit --resubmit_expired` sends that
search again with the same settings.

  python scripts/ncbi_blast_jobs.py submit --jobs JOBS.tsv --name NAME --query Q.fna \\
      [--entrez "txid10239[Organism:exp]"] [--db core_nt] [--evalue 1e-10] [--hits 20]
  python scripts/ncbi_blast_jobs.py fetch --jobs JOBS.tsv
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
COLS = ["name", "rid", "query", "query_sha256", "db", "entrez", "evalue", "hits", "submitted_utc", "status"]


def read_jobs(path: Path) -> list[dict]:
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    return [dict(zip(lines[0].split("\t"), l.split("\t"))) for l in lines[1:] if l.strip()]


def write_jobs(path: Path, rows: list[dict]) -> None:
    path.write_text("\t".join(COLS) + "\n" + "".join("\t".join(r.get(c, "") for c in COLS) + "\n" for r in rows))


def http(params: dict, data: bool = False) -> str:
    body = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(URL, data=body if data else None, method="POST" if data else "GET") if data \
        else urllib.request.Request(URL + "?" + body.decode())
    with urllib.request.urlopen(req, timeout=300) as r:
        return r.read().decode(errors="replace")


def put(row: dict) -> str:
    q = Path(row["query"]).read_text()
    params = {"CMD": "Put", "PROGRAM": "blastn", "MEGABLAST": "on", "DATABASE": row["db"], "QUERY": q,
              "EXPECT": row["evalue"], "HITLIST_SIZE": row["hits"], "TOOL": "prophound"}
    if row["entrez"]:
        params["ENTREZ_QUERY"] = row["entrez"]
    page = http(params, data=True)
    m = re.search(r"RID = (\S+)", page)
    if not m:
        raise RuntimeError(f"no RID in NCBI's response: {page[:200]!r}")
    return m.group(1)


def cmd_submit(a) -> None:
    jobs = Path(a.jobs)
    rows = read_jobs(jobs)
    if a.resubmit_expired:
        todo = [r for r in rows if r["status"] in ("expired", "failed")]
    else:
        if any(r["name"] == a.name for r in rows):
            print(f"[submit] {a.name}: already in {jobs}"); return
        row = {"name": a.name, "query": str(Path(a.query)), "db": a.db, "entrez": a.entrez or "",
               "evalue": str(a.evalue), "hits": str(a.hits), "status": "submitted",
               "query_sha256": hashlib.sha256(Path(a.query).read_bytes()).hexdigest()}
        rows.append(row); todo = [row]
    for r in todo:
        r["rid"] = put(r)
        r["submitted_utc"] = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        r["status"] = "submitted"
        write_jobs(jobs, rows)
        print(f"[submit] {r['name']}: RID {r['rid']}", flush=True)
        time.sleep(11)


def cmd_fetch(a) -> None:
    jobs = Path(a.jobs)
    rows = read_jobs(jobs)
    for r in rows:
        out = jobs.parent / f"blast_{r['name']}.tabular.txt"
        if r["status"] == "done" and out.exists():
            continue
        info = http({"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": r["rid"]})
        st = re.search(r"Status=(\w+)", info)
        st = st.group(1) if st else "UNKNOWN"
        if st == "READY":
            res = http({"CMD": "Get", "FORMAT_TYPE": "Tabular", "RID": r["rid"], "ALIGNMENT_VIEW": "Tabular",
                        "DESCRIPTIONS": r["hits"], "ALIGNMENTS": r["hits"]})
            out.write_text(res)
            r["status"] = "done"
        elif st == "WAITING":
            r["status"] = "waiting"
        elif st == "FAILED":
            r["status"] = "failed"
        else:
            r["status"] = "expired"
        print(f"[fetch] {r['name']}: {r['status']}", flush=True)
        write_jobs(jobs, rows)
        time.sleep(61)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("submit")
    s.add_argument("--jobs", required=True); s.add_argument("--name"); s.add_argument("--query")
    s.add_argument("--db", default="core_nt"); s.add_argument("--entrez", default="")
    s.add_argument("--evalue", type=float, default=1e-10); s.add_argument("--hits", type=int, default=20)
    s.add_argument("--resubmit_expired", action="store_true")
    f = sub.add_parser("fetch"); f.add_argument("--jobs", required=True)
    a = ap.parse_args()
    if a.cmd == "submit" and not a.resubmit_expired and not (a.name and a.query):
        ap.error("submit needs --name and --query (or --resubmit_expired)")
    {"submit": cmd_submit, "fetch": cmd_fetch}[a.cmd](a)


if __name__ == "__main__":
    main()
