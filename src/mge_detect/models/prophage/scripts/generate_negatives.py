import random
from bisect import bisect_left
from pathlib import Path

POS_WIN_BED = "windows/bed/pos_win.bed"
STRICT_NEG_BED = "labels/negatives/all.neg.strict.bed"
POS_BED = "positives/geNomad.norm.bed"
OUT_SOFT = "windows/bed/neg_win.soft.bed"
OUT_HARD = "windows/bed/neg_win.hard.bed"
OUT_ALL = "windows/bed/neg_win.rebuilt.bed"

SOFT_MIN = 50000
HARD_MAX = 5000
HARD_FRAC = 0.40
RNG_SEED = 42
MAX_TRIES_SAME_CONTIG = 200
MAX_TRIES_GLOBAL = 2000

random_generator = random.Random(RNG_SEED)

positives_by_chrom = {}
with open(POS_BED) as file_handle:
    for line in file_handle:
        if not line.strip():
            continue
        chrom, start_str, end_str, name = line.rstrip().split("\t")
        start, end = int(start_str), int(end_str)
        positives_by_chrom.setdefault(chrom, []).append((start, end))
for chrom in positives_by_chrom:
    positives_by_chrom[chrom].sort()

positive_start_positions = {}
for chrom, intervals in positives_by_chrom.items():
    starts_for_chrom = []
    for interval_start, interval_end in intervals:
        starts_for_chrom.append(interval_start)
    positive_start_positions[chrom] = starts_for_chrom


def overlaps(chrom, query_start, query_end):
    intervals = positives_by_chrom.get(chrom, [])
    if not intervals:
        return False
    insertion_index = bisect_left(positive_start_positions[chrom], query_end)
    for neighbor_index in (insertion_index - 1, insertion_index, insertion_index + 1):
        if 0 <= neighbor_index < len(intervals):
            pos_start, pos_end = intervals[neighbor_index]
            if not (pos_end <= query_start or pos_start >= query_end):
                return True
    return False


def min_distance_to_positive(chrom, query_start, query_end):
    intervals = positives_by_chrom.get(chrom, [])
    if not intervals:
        return float("inf")
    if overlaps(chrom, query_start, query_end):
        return 0
    insertion_index = bisect_left(positive_start_positions[chrom], query_start)
    distance_bp = float("inf")
    if insertion_index - 1 >= 0:
        prev_start, prev_end = intervals[insertion_index - 1]
        distance_bp = min(distance_bp, query_start - prev_end)
    if insertion_index < len(intervals):
        next_start, next_end = intervals[insertion_index]
        distance_bp = min(distance_bp, next_start - query_end)
    return max(0, distance_bp)


negative_regions_by_chrom = {}
with open(STRICT_NEG_BED) as file_handle:
    for line in file_handle:
        if not line.strip():
            continue
        chrom, start_str, end_str = line.rstrip().split("\t")[:3]
        start, end = int(start_str), int(end_str)
        negative_regions_by_chrom.setdefault(chrom, []).append((start, end))
for chrom in negative_regions_by_chrom:
    negative_regions_by_chrom[chrom].sort()


def sample_within_chrom(chrom, window_length_bp):
    regions = negative_regions_by_chrom.get(chrom, [])
    if not regions:
        return None
    candidate_regions = [
        (region_start, region_end)
        for (region_start, region_end) in regions
        if region_end - region_start >= window_length_bp
    ]
    if not candidate_regions:
        return None
    region_start, region_end = random_generator.choice(candidate_regions)
    neg_start = random_generator.randint(region_start, region_end - window_length_bp)
    return neg_start, neg_start + window_length_bp


all_chromosomes = list(negative_regions_by_chrom.keys())


def sample_globally(window_length_bp):
    for attempt_index in range(2000):
        selected_chromosome = random_generator.choice(all_chromosomes)
        sampled_interval = sample_within_chrom(selected_chromosome, window_length_bp)
        if sampled_interval:
            neg_start, neg_end = sampled_interval
            return selected_chromosome, neg_start, neg_end
    return None


positive_windows = []
with open(POS_WIN_BED) as file_handle:
    for window_id, line in enumerate(file_handle, start=1):
        chrom, start_str, end_str, window_name = line.rstrip().split("\t")
        start, end = int(start_str), int(end_str)
        positive_windows.append((window_id, chrom, start, end, end - start))
num_positive_windows = len(positive_windows)
target_num_hard_negatives = int(round(HARD_FRAC * num_positive_windows))

hard_negatives = []
for window_id, chrom, start, end, window_length_bp in positive_windows:
    placed = False
    for attempt_index in range(MAX_TRIES_SAME_CONTIG):
        sampled_interval = sample_within_chrom(chrom, window_length_bp)
        if not sampled_interval:
            break
        neg_start, neg_end = sampled_interval
        if overlaps(chrom, neg_start, neg_end):
            continue
        distance_bp = min_distance_to_positive(chrom, neg_start, neg_end)
        if distance_bp <= HARD_MAX:
            hard_negatives.append((chrom, neg_start, neg_end, f"negH|{window_id}"))
            placed = True
            break
    if placed:
        continue
    for attempt_index in range(MAX_TRIES_GLOBAL):
        global_sample = sample_globally(window_length_bp)
        if not global_sample:
            break
        neg_chrom, neg_start, neg_end = global_sample
        if overlaps(neg_chrom, neg_start, neg_end):
            continue
        distance_bp = min_distance_to_positive(neg_chrom, neg_start, neg_end)
        if distance_bp <= HARD_MAX:
            hard_negatives.append((neg_chrom, neg_start, neg_end, f"negH|{window_id}"))
            placed = True
            break
if len(hard_negatives) > target_num_hard_negatives:
    hard_negatives = hard_negatives[:target_num_hard_negatives]

hard_ids = {int(record[3].split("|")[1]) for record in hard_negatives}


def want_soft_negative(chrom, neg_start, neg_end):
    return (not overlaps(chrom, neg_start, neg_end)) and (
        min_distance_to_positive(chrom, neg_start, neg_end) >= SOFT_MIN
    )


soft_negatives = []
for window_id, chrom, start, end, window_length_bp in positive_windows:
    if window_id in hard_ids:
        continue
    placed = False
    for attempt_index in range(MAX_TRIES_SAME_CONTIG):
        sampled_interval = sample_within_chrom(chrom, window_length_bp)
        if not sampled_interval:
            break
        neg_start, neg_end = sampled_interval
        if want_soft_negative(chrom, neg_start, neg_end):
            soft_negatives.append((chrom, neg_start, neg_end, f"negS|{window_id}"))
            placed = True
            break
    if placed:
        continue
    for attempt_index in range(MAX_TRIES_GLOBAL):
        global_sample = sample_globally(window_length_bp)
        if not global_sample:
            break
        neg_chrom, neg_start, neg_end = global_sample
        if want_soft_negative(neg_chrom, neg_start, neg_end):
            soft_negatives.append((neg_chrom, neg_start, neg_end, f"negS|{window_id}"))
            placed = True
            break

relaxation_thresholds = [40000, 30000, 20000, 10000, 5000, 0]
relaxation_index = 0
while len(hard_negatives) + len(
    soft_negatives
) < num_positive_windows and relaxation_index < len(relaxation_thresholds):
    threshold = relaxation_thresholds[relaxation_index]
    relaxation_index += 1
    existing_ids = {int(x[3].split("|")[1]) for x in soft_negatives} | hard_ids
    for window_id, chrom, start, end, window_length_bp in positive_windows:
        if window_id in existing_ids:
            continue
        tries = 0
        while tries < 1000:
            global_sample = sample_globally(window_length_bp)
            tries += 1
            if not global_sample:
                break
            neg_chrom, neg_start, neg_end = global_sample
            if overlaps(neg_chrom, neg_start, neg_end):
                continue
            if min_distance_to_positive(neg_chrom, neg_start, neg_end) >= threshold:
                soft_negatives.append(
                    (neg_chrom, neg_start, neg_end, f"negS|{window_id}")
                )
                break

Path(OUT_SOFT).parent.mkdir(parents=True, exist_ok=True)
with open(OUT_SOFT, "w") as file_handle:
    for record in soft_negatives:
        file_handle.write("\t".join(map(str, record)) + "\n")
with open(OUT_HARD, "w") as file_handle:
    for record in hard_negatives:
        file_handle.write("\t".join(map(str, record)) + "\n")


def sort_key_by_window_id(record):
    return int(record[3].split("|")[1])


all_negatives = sorted(hard_negatives + soft_negatives, key=sort_key_by_window_id)
with open(OUT_ALL, "w") as file_handle:
    for record in all_negatives:
        file_handle.write("\t".join(map(str, record)) + "\n")

print(
    f"hard: {len(hard_negatives)}   soft: {len(soft_negatives)}   total: {len(all_negatives)}   target: {num_positive_windows}"
)

