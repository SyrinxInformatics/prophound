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
"""Submit genomes to the PHASTEST API and fetch the results later, from any session.

The jobs file is the only state: one row per submitted sequence (a whole genome, or one contig
with --targets) with its PHASTEST job ID, the input's SHA-256 and the submission time. `submit`
appends a row as soon as each job is accepted, skips sequences already in the file, and logs a
rejected submission (HTTP error) without stopping, so a re-run retries only those. `fetch` reads only the jobs file, saves each finished
job's JSON (and its zip of region details) under --out, and can be re-run until every job is
complete. Results also stay viewable at https://phastest.ca/submissions/<job_id>.

API (docs/data-reliability.md): raw POST of the FASTA body to https://phastest.ca/phastest_api
(`?contigs=1` for a multi-FASTA), then GET https://phastest.ca/phastest_api?acc=<job_id>.

  python scripts/phastest_jobs.py submit --genomes LIST --fasta_dir DIR --jobs JOBS.tsv
  python scripts/phastest_jobs.py submit --targets GENOME_CONTIG.tsv --fasta_dir DIR --jobs JOBS.tsv
  python scripts/phastest_jobs.py fetch --jobs JOBS.tsv --out DIR
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import time
import urllib.request
from pathlib import Path

API = "https://phastest.ca/phastest_api"
COLS = ["genome", "contig", "job_id", "input_sha256", "n_contigs", "submitted_utc"]


def read_jobs(path: Path) -> list[dict]:
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    return [dict(zip(lines[0].split("\t"), l.split("\t"))) for l in lines[1:] if l.strip()]


def read_fasta(path: Path) -> dict[str, bytes]:
    out, name, buf = {}, None, []
    for line in path.read_bytes().splitlines(keepends=True):
        if line.startswith(b">"):
            if name:
                out[name] = b"".join(buf)
            name, buf = line[1:].split()[0].decode(), [line]
        else:
            buf.append(line)
    if name:
        out[name] = b"".join(buf)
    return out


def cmd_submit(a) -> None:
    jobs = Path(a.jobs)
    done = {(r["genome"], r.get("contig", "")) for r in read_jobs(jobs)}
    if not jobs.exists():
        jobs.parent.mkdir(parents=True, exist_ok=True)
        jobs.write_text("\t".join(COLS) + "\n")
    if a.targets:
        targets = [tuple(l.rstrip("\n").split("\t")[:2]) for l in open(a.targets) if l.strip() and not l.startswith("genome\t")]
    else:
        targets = [(l.strip(), "") for l in open(a.genomes) if l.strip()]
    cache: dict[str, dict[str, bytes]] = {}
    for g, c in targets:
        if (g, c) in done:
            continue
        path = Path(a.fasta_dir) / f"{g}.fna"
        if c:
            cache.setdefault(g, read_fasta(path))
            body = cache[g][c]
        else:
            body = path.read_bytes()
        n = body.count(b">")
        url = API + ("?contigs=1" if n > 1 else "")
        req = urllib.request.Request(url, data=body, method="POST", headers={"Content-Type": "text/plain"})
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                resp = json.loads(r.read().decode())
        except (OSError, ValueError) as e:
            print(f"[submit] {g} {c}: rejected ({e}); re-run to retry", flush=True)
            time.sleep(a.pause)
            continue
        if "job_id" not in resp:
            print(f"[submit] {g} {c}: no job_id in response: {str(resp)[:200]}", flush=True)
            continue
        row = [g, c, resp["job_id"], hashlib.sha256(body).hexdigest(), str(n), dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")]
        with open(jobs, "a") as fh:
            fh.write("\t".join(row) + "\n")
        print(f"[submit] {g} {c}: {resp['job_id']}", flush=True)
        time.sleep(a.pause)


def cmd_fetch(a) -> None:
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    states = {}
    for r in read_jobs(Path(a.jobs)):
        key = r["genome"] + (f"__{r['contig']}" if r.get("contig") else "")
        f = out / f"{key}.json"
        if f.exists() and json.loads(f.read_text()).get("status") == "Complete":
            states[key] = "Complete"
            continue
        try:
            with urllib.request.urlopen(f"{API}?acc={r['job_id']}", timeout=120) as resp:
                d = json.loads(resp.read().decode())
        except (OSError, ValueError) as e:
            states[key] = f"fetch error: {e}"
            continue
        st = str(d.get("status", "?"))
        states[key] = st
        if st == "Complete":
            f.write_text(json.dumps(d, indent=1))
            z = d.get("zip")
            if z:
                try:
                    urllib.request.urlretrieve(("https://" + z) if not z.startswith("http") else z, out / f"{key}.zip")
                except OSError as e:
                    print(f"[fetch] {key}: zip not saved ({e})")
        time.sleep(a.pause)
    counts = {}
    for s in states.values():
        counts[s.split(":")[0][:40]] = counts.get(s.split(":")[0][:40], 0) + 1
    print(f"[fetch] {len(states)} jobs: {counts}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("submit")
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--genomes", help="one genome per line: whole genomes")
    g.add_argument("--targets", help="genome<TAB>contig per line: one submission per contig")
    s.add_argument("--fasta_dir", required=True)
    s.add_argument("--jobs", required=True); s.add_argument("--pause", type=float, default=10.0)
    f = sub.add_parser("fetch")
    f.add_argument("--jobs", required=True); f.add_argument("--out", required=True)
    f.add_argument("--pause", type=float, default=2.0)
    a = ap.parse_args()
    {"submit": cmd_submit, "fetch": cmd_fetch}[a.cmd](a)


if __name__ == "__main__":
    main()
