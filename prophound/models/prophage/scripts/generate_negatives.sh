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

WIN=256000

mkdir -p windows/bed windows/meta

awk -F'\t' -v OFS='\t' '{mid=int(($2+$3)/2); print $1, mid, mid+1, "pos|"NR}' positives/geNomad.norm.bed \
| bedtools slop -g labels/all.sizes -b $WIN \
| sort -k1,1 -k2,2n > windows/bed/pos_win.bed
NPOS=$(wc -l < windows/bed/pos_win.bed)
echo "positive windows: $NPOS"

cat labels/negatives/*.neg.bed > labels/negatives/all.neg.bed
bedtools shuffle -i windows/bed/pos_win.bed -g labels/all.sizes \
                 -incl labels/negatives/all.neg.bed -noOverlapping -chrom \
| awk -F'\t' -v OFS='\t' '{split($4,a,"|"); print $1,$2,$3,"neg|"a[2]}' \
> windows/bed/neg_win.bed
NNEG=$(wc -l < windows/bed/neg_win.bed)
echo "negative windows: $NNEG"

cat windows/bed/pos_win.bed windows/bed/neg_win.bed \
| sort -k1,1 -k2,2n > windows/bed/all_win.bed

