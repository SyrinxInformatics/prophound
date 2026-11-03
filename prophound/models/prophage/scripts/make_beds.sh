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

cd /work/pke835/mge_explorer

find genomad_out -type f -name '*_provirus.tsv' -print0 \
| while IFS= read -r -d '' t; do
    base=$(basename "$(dirname "$(dirname "$t")")")
    awk -F'\t' -v OFS='\t' -v G="$base" '
      FNR==1{
        seq=1; b=2; e=3;
        for(i=1;i<=NF;i++){
          if($i ~ /^(seq_name|sequence|contig)$/) seq=i
          if($i ~ /^(start|begin)$/)               b=i
          if($i ~ /^(end|stop)$/)                  e=i
        }
        next
      }
      ($b~/^[0-9]+$/ && $e~/^[0-9]+$/){
        print G "|" $seq, $b-1, $e, "geNomad"
      }' "$t"
  done \
| sort -k1,1 -k2,2n > positives/geNomad.bed

echo "Lines in positives/geNomad.bed: $(wc -l < positives/geNomad.bed)"
head positives/geNomad.bed

