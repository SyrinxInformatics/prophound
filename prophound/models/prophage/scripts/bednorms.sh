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

awk -F'\t' -v OFS='\t' '
{
  # $1 = genome|seq_name     where seq_name may be "contig" or "contig|provirus_..."
  split($1,a,"|");
  genome=a[1];
  contig=a[2];
  prov=(a[3] ~ /^provirus_/) ? a[3] : "NA";
  name = (prov=="NA" ? "geNomad" : "geNomad|" prov);
  print genome "|" contig, $2, $3, name
}' positives/geNomad.bed \
| sort -k1,1 -k2,2n > positives/geNomad.norm.bed

echo "Normalized: $(wc -l < positives/geNomad.norm.bed)"
head positives/geNomad.norm.bed

