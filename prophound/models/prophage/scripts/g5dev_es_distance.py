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
"""Edge-snap distance vs added sequence on corpus MAG shreds (freeze decision, 2026-10-05).

For each mode and snap distance d, MAG covered50 (three shred seeds pooled), total bp
called, and bp called outside geNomad elements, on the 42 (tune) and 18 (check).
Covered50 and false/Mb cannot see how much longer calls get; this does.
Run from prophound/models/prophage:  python scripts/g5dev_es_distance.py
"""
from pathlib import Path

import evalkit as ek
import pandas as pd
from evalkit import cov_frac
from g5dev_compare import off_element_bp
from g5dev_eval import edge_snap, load, to_genome

Lt=pd.read_csv('analysis/g5dev/contig_lengths.tsv',sep='\t',dtype={'contig_id':str})
for sname,path in [('tune','scripts/eval_sizestrat_elements.csv'),('check','scripts/eval_confirm18_elements.csv')]:
    T=pd.read_csv(path).rename(columns={'accession':'assembly','contig':'sequence'})[['assembly','sequence','start','end']].reset_index(drop=True)
    G=sorted(T.assembly.unique()); acc={}
    for cond in ['mag_s0','mag_s1','mag_s2']:
        ens,seeds=load(Path('analysis/g5dev/prophound'),cond,G)
        lc=Lt[(Lt.condition==cond)&Lt.genome.isin(G)]; L={(g,c):int(n) for g,c,n in zip(lc.genome,lc.contig_id,lc.length)}
        for mode in ['default','exploratory','conservative']:
            base=ek.shipped_calls(ens,seeds,mode)
            for d in [0,1000,2000,5000,10000]:
                c=to_genome(edge_snap(base,L,d)); h=cov_frac(T,c)>=0.5
                a=acc.setdefault((mode,d),[0,0,0,0]); a[0]+=h.sum(); a[1]+=len(h); a[2]+=(c.end-c.start).sum(); a[3]+=off_element_bp(c,T).sum()
    print('==',sname)
    for (mode,d),a in acc.items():
        b=acc[(mode,0)]
        print(f'{mode:12s} ES {d:5d}: MAG {a[0]/a[1]:.3f}  called {a[2]/1e6:6.2f} Mb (+{(a[2]-b[2])/b[2]*100:4.1f}%)  off-element +{(a[3]-b[3])/1e6:.2f} Mb')
