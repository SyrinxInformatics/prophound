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

cd /work/$USER/mge_explorer/plasmids/chunks
mkdir -p reruns
for CH in files_chunk_0001.list files_chunk_0002.list files_chunk_0003.list files_chunk_0004.list; do
  CHFILE="$CH"
  DONE="reruns/${CH}.done.txt"
  REM="reruns/${CH}.remain.txt"
  : > "$DONE"; : > "$REM"
  while read -r g; do
    b=$(basename "$g" .fna)
    t="/work/$USER/mge_explorer/plasmids/genome_chunks/genomad_out/$b/${b}_find_proviruses/${b}_provirus.tsv"
    if [ -s "$t" ]; then echo "$g" >> "$DONE"; else echo "$g" >> "$REM"; fi
  done < "$CHFILE"
  echo "$CH  done: $(wc -l < "$DONE")   remain: $(wc -l < "$REM")"
done
