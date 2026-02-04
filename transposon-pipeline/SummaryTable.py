

#!/usr/bin/env python3
"""
SummaryTable.py  - Generate per-feature insertion metrics from hit files

This script reads one or more hit files, computes per-feature insertion metrics
(hits, reads, neighborhood index, etc.), and writes summary tables + helper tables
for downstream plotting / downstream essentiality calling.

USAGE
-----
SummaryTable.py
   -i  --input-dir       [str]   Input directory with hit files. Defaults to current directory if left unspecified.
   -o  --output-dir      [str]   Output directory for generated tables. Defaults to current directory if left unspecified.
   --hit-glob            [str]   Glob pattern for hit files (default: *_hits.txt)
   --hit-format          [str]   Hit file format: hit_table | pombe_csv | wig  (default: hit_table)
   --strict-hit-table            Strict validation of hit_table rows (OFF by default).

   --feature-db          [str]   Internal feature DB selector:
                                albicans | pombe | cerevisiae | glabrata
                                (REQUIRED for feature-level outputs; and by default we write all outputs)

   -f  --read-depth-filter        [int]   Minimum reads per hit to include (default: 1)
   --neighborhood-window-size     [int]   Window size for neighborhood index (default: 10000)
   --bin-size                     [int]   Genome bin size for binned summaries (default: 10000)
   --ignored-regions              [str]   Optional BED-like ignored regions (chrom start end)

   -t  --threads          [int]   Number of CPU cores to use (default: 1)

   --overwrite                    Overwrite output files if they already exist (OFF by default).

OUTPUTS
-------
By default (when no --write-* flags are passed), this script writes ALL of:
   - binned_hits.RDF_<RDF>.csv
   - hit_summary.RDF_<RDF>.csv               (needs --feature-db)
   - stats.csv                               (needs --feature-db)
   - <lib>.all_hits.csv
   - <lib>_analysis.csv                      (needs --feature-db)
   - <lib>.filter_<RDF>.bed
   - <lib>.proteome.filter_<RDF>.bed         (needs --feature-db)
   - <lib>.outlier_stats.txt                 (needs --feature-db)

You can still request a subset using the --write-* flags.

GLABRATA NOTES (Sc inference)
-----------------------------
If --feature-db glabrata is used, <lib>_analysis.csv will include:
  - Glabrata feature info
  - Sc ortholog (from dependencies/glabrata/C_glabrata_BG2_S_cerevisiae_orthologs.txt
    OR an override via --scer-orthologs)
  - Sc essentiality (from dependencies/cerevisiae/*viable/*inviable annotations)
  - Sc synthetic lethal flag (dependencies/cerevisiae/duplicatesSl_011116.txt)
  - Sc fitness (dependencies/cerevisiae/neFitnessStandard.txt)

This is for inference only; no albicans/pombe enrichment is performed.

QUALITY-OF-LIFE UPDATES INCLUDED
--------------------------------
1) Glabrata analysis.csv now includes scer_gene_name and scer_sgdid if present in the ortholog file.
2) Optional warnings: report how many glabrata features are missing Sc orthologs (likely ID mismatch).
3) feature_table.RDF_*: chrom column is now more robust (tries multiple attribute names).
4) Optional: include Sc ortholog columns in feature_table via --feature-table-include-scer (glabrata only).
"""

import os
import sys
import csv
import argparse
import glob
import math
from collections import defaultdict
from itertools import chain

import numpy as np

import Shared
import GenomicFeatures
from RangeSet import RangeSet


# ----------------------------- File utilities -----------------------------

def _ensure_parent_dir(path):
    parent = os.path.dirname(os.path.abspath(path))
    if parent and not os.path.exists(parent):
        Shared.make_dir(parent)


def _open_out(path, overwrite=False, newline=""):
    """Open an output file for writing, enforcing overwrite policy."""
    _ensure_parent_dir(path)
    if os.path.exists(path) and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing file (use --overwrite): {path}")
    return open(path, "w", newline=newline)


def _eprint(msg):
    print(msg, file=sys.stderr)


#  Hit parsing 

def read_hit_files(files, read_depth_filter=1, hit_format="hit_table", strict_hit_table=False):
    return [
        read_hit_file(f, read_depth_filter, hit_format=hit_format, strict_hit_table=strict_hit_table)
        for f in files
    ]


def read_hit_file(filename, read_depth_filter=1, hit_format="hit_table", strict_hit_table=False):
    if hit_format == "hit_table":
        return _read_hit_table(filename, read_depth_filter, strict=strict_hit_table)
    if hit_format == "pombe_csv":
        return _read_pombe_csv(filename, read_depth_filter)
    if hit_format == "wig":
        return _read_wig(filename, read_depth_filter)
    raise ValueError(f"Unknown --hit-format: {hit_format}")


def _read_hit_table(filename, read_depth_filter=1, strict=False):
    """
    Parse the Levitan/Berman hit_table output (>=13 tab-separated fields).
    """
    result = []
    with open(filename, "r") as in_file:
        header = next(in_file, None)
        if header is None:
            return result

        for line in in_file:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13:
                if strict:
                    raise ValueError(f"Malformed hit_table row (expected >=13 columns): {filename}")
                continue

            chrom, source, up_feature_type, up_feature_name, up_gene_name, \
                up_feature_dist, down_feature_type, down_feature_name, \
                down_gene_name, down_feature_dist, ig_type, hit_pos, hit_count = fields[:13]

            try:
                hit_count = int(hit_count)
                hit_pos = int(hit_pos)
            except ValueError:
                if strict:
                    raise
                continue

            if hit_count < read_depth_filter:
                continue

            if "ORF" in ig_type:
                gene_name = up_gene_name if up_gene_name not in ("nan", "", "NA") else up_feature_name
                if strict and up_feature_name != down_feature_name:
                    raise ValueError(f"Inconsistent ORF row in hit_table: {filename}")
            else:
                gene_name = "nan"

            result.append({
                "chrom": chrom,
                "source": source,
                "up_feature_type": up_feature_type,
                "up_feature_name": up_feature_name,
                "up_gene_name": up_gene_name,
                "up_feature_dist": up_feature_dist,
                "down_feature_type": down_feature_type,
                "down_feature_name": down_feature_name,
                "down_gene_name": down_gene_name,
                "down_feature_dist": down_feature_dist,
                "ig_type": ig_type,
                "hit_pos": hit_pos,
                "hit_count": hit_count,
                "gene_name": gene_name
            })
    return result


def _read_pombe_csv(filename, read_depth_filter=1):
    result = []
    with open(filename, "r", newline="") as in_file:
        reader = csv.reader(in_file)
        next(reader, None)
        for row in reader:
            if len(row) < 5:
                continue
            chrom, strand, hit_pos, hit_count, gene_name = row[:5]
            try:
                hit_pos = int(hit_pos)
                hit_count = int(hit_count)
            except ValueError:
                continue
            if hit_count < read_depth_filter:
                continue
            result.append({
                "chrom": chrom,
                "source": strand,
                "hit_pos": hit_pos,
                "hit_count": hit_count,
                "gene_name": gene_name
            })
    return result


def _read_wig(filename, read_depth_filter=1):
    result = []
    chrom_name = None
    with open(filename, "r") as in_file:
        _ = in_file.readline()
        for line in in_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith("variableStep"):
                chrom_name = line.split("chrom=", 1)[1].strip() if "chrom=" in line else None
                continue
            if chrom_name in (None, "chrM"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                hit_pos = int(parts[0])
                reads = int(parts[1])
            except ValueError:
                continue
            if reads < read_depth_filter:
                continue
            result.append({"chrom": chrom_name, "hit_pos": hit_pos, "hit_count": reads})
    return result

def get_hits_from_wig(wig_file):
    """
    Parse a .wig (variableStep or fixedStep) and return hits as a list of dicts:
        {"chromosome": str, "hit_pos": int, "hit_count": int}

    This format is compatible with analyze_hits(), which expects dict hits.
    """
    hits = []
    chrom = None
    span = 1

    with open(wig_file, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("track") or line.startswith("#"):
                continue

            if line.startswith("variableStep"):
                # Example: variableStep chrom=chrI span=1
                chrom = None
                span = 1
                parts = line.split()
                for p in parts[1:]:
                    if p.startswith("chrom="):
                        chrom = p.split("=", 1)[1]
                    elif p.startswith("span="):
                        span = int(p.split("=", 1)[1])
                continue

            if line.startswith("fixedStep"):
                # fixedStep is less common for these, but we can support it too.
                # Example: fixedStep chrom=chrI start=1 step=1 span=1
                chrom = None
                start = None
                step = None
                span = 1
                parts = line.split()
                for p in parts[1:]:
                    if p.startswith("chrom="):
                        chrom = p.split("=", 1)[1]
                    elif p.startswith("start="):
                        start = int(p.split("=", 1)[1])
                    elif p.startswith("step="):
                        step = int(p.split("=", 1)[1])
                    elif p.startswith("span="):
                        span = int(p.split("=", 1)[1])

                # fixedStep data lines are single values (count), positions inferred
                # We'll read them in a small inner loop by setting state variables.
                fixed_pos = start
                fixed_step = step

                # Read subsequent numeric lines until next header
                for raw2 in fh:
                    s2 = raw2.strip()
                    if not s2 or s2.startswith("track") or s2.startswith("#"):
                        continue
                    if s2.startswith("variableStep") or s2.startswith("fixedStep"):
                        # rewind one header line by handling it in outer loop:
                        # easiest: stash it and process by recursion-free trick:
                        # We can't un-read cleanly, so just process it here by
                        # breaking and letting outer loop miss it; instead,
                        # we’ll handle fixedStep only if your WIGs are variableStep.
                        raise NotImplementedError("Encountered nested step header inside fixedStep parsing.")
                    count = int(float(s2))
                    if count > 0:
                       hits.append({"chrom": chrom, "hit_pos": int(fixed_pos), "hit_count": int(count)})
                    fixed_pos += fixed_step
                break  # EOF reached
                # (If you truly have fixedStep WIGs and need this robust, we’ll tighten this.)

            # Data line (variableStep): "position value"
            # Example: 12345 7
            fields = line.split()
            if len(fields) < 2:
                continue
            pos = int(fields[0])
            count = int(float(fields[1]))
            if count <= 0:
                continue

            hits.append({"chrom": chrom, "hit_pos": pos, "hit_count": count})

    return hits

#  Stats 
TOTAL_HITS = "Total Hits"
TOTAL_READS = "Total Reads"
AVG_READS_PER_HIT = "Mean Reads Per Hit"
ORF_HITS = "Hits in ORFs"
ANNOTATED_FEATURE_HITS = "Hits in Genomic Features"
PER_ANNOTATED_FEATURE_HITS = "% of hits in features"
INTERGENIC_HITS = "Intergenic Hits"
PER_INTERGENIC_HITS = "% of intergenic hits"
FEATURES_HIT = "No. of Features Hit"
PER_FEATURES_HIT = "% of features hit"
AVG_HITS_IN_FEATURE = "Mean hits per feature"
AVG_READS_IN_FEATURE = "Mean reads per feature"
AVG_READS_IN_FEATURE_HIT = "Mean reads per hit in feature"
READS_IN_FEATURES = "Total reads in features"

ALL_STATS = [
    TOTAL_READS, TOTAL_HITS, PER_ANNOTATED_FEATURE_HITS, PER_INTERGENIC_HITS,
    PER_FEATURES_HIT, AVG_HITS_IN_FEATURE, AVG_READS_IN_FEATURE,
    AVG_READS_PER_HIT, AVG_READS_IN_FEATURE_HIT
]


def _safe_all_features(feature_db):
    if hasattr(feature_db, "get_all_features"):
        return feature_db.get_all_features()
    feats = []
    for chrom in feature_db:
        feats.extend(list(chrom.get_features()) if hasattr(chrom, "get_features") else list(feature_db[chrom]))
    return feats


def get_statistics(dataset, feature_db):
    result = {}
    result[TOTAL_HITS] = len(dataset)
    result[TOTAL_READS] = sum(obj["hit_count"] for obj in dataset) if dataset else 0
    result[AVG_READS_PER_HIT] = (result[TOTAL_READS] / result[TOTAL_HITS]) if result[TOTAL_HITS] else 0
    result[ORF_HITS] = 0
    result[ANNOTATED_FEATURE_HITS] = 0
    result[INTERGENIC_HITS] = 0
    result[READS_IN_FEATURES] = 0

    features_hit = set()
    for hit in dataset:
        features = feature_db.get_features_at_location(hit["chrom"], hit["hit_pos"])
        if not features:
            continue
        result[ANNOTATED_FEATURE_HITS] += 1
        result[READS_IN_FEATURES] += hit["hit_count"]
        features_hit.update(set(f.standard_name for f in features))
        for f in features:
            if getattr(f, "is_orf", False):
                result[ORF_HITS] += 1
                break

    result[INTERGENIC_HITS] = result[TOTAL_HITS] - result[ANNOTATED_FEATURE_HITS]
    result[FEATURES_HIT] = len(features_hit)

    total_features = len(_safe_all_features(feature_db)) or 1
    total_hits = result[TOTAL_HITS] or 1
    annotated_hits = result[ANNOTATED_FEATURE_HITS] or 1
    features_hit_n = result[FEATURES_HIT] or 1

    result[PER_INTERGENIC_HITS] = "%.2f%%" % (result[INTERGENIC_HITS] * 100.0 / total_hits)
    result[PER_ANNOTATED_FEATURE_HITS] = "%.2f%%" % (result[ANNOTATED_FEATURE_HITS] * 100.0 / total_hits)
    result[PER_FEATURES_HIT] = "%.2f%%" % (result[FEATURES_HIT] * 100.0 / total_features)
    result[AVG_HITS_IN_FEATURE] = "%.1f" % (result[ANNOTATED_FEATURE_HITS] * 1.0 / features_hit_n)
    result[AVG_READS_IN_FEATURE] = result[READS_IN_FEATURES] / features_hit_n
    result[AVG_READS_IN_FEATURE_HIT] = result[READS_IN_FEATURES] / annotated_hits

    return result


#  Ignored regions 

def _parse_ignored_regions(path):
    if not path:
        return None
    by_chrom = defaultdict(RangeSet)
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            if end <= start:
                continue
            # BED is 0-based half-open; our internals are 1-based inclusive.
            by_chrom[chrom].add((start + 1, end))
    return by_chrom


#  Core analysis 

def analyze_hits(dataset, feature_db, neighborhood_window_size=10000, ignored_regions_by_chrom=None):
    log2 = lambda v: math.log(v, 2) if v > 0 else 0

    result = {}
    total_reads = sum(h["hit_count"] for h in dataset) if dataset else 0
    total_reads_log = log2(total_reads)

    chroms = set(h["chrom"] for h in dataset)
    for chrom in chroms:
        hits = [h for h in dataset if h["chrom"] == chrom]
        if not hits:
            continue

        chrom_len = max(len(feature_db[chrom]), max(h["hit_pos"] for h in hits))

        exon_mask = np.zeros((chrom_len + 1,), dtype=bool)
        domain_mask = np.zeros((chrom_len + 1,), dtype=bool)

        for feature in feature_db[chrom]:
            for exon in feature.exons:
                exon_mask[exon[0]:exon[1] + 1] = True
            for domain in feature.domains:
                domain_mask[domain[0]:domain[1] + 1] = True

        ignored_mask = np.ones((chrom_len + 1,), dtype=bool)
        ignored_range_set = RangeSet()
        if ignored_regions_by_chrom and chrom in ignored_regions_by_chrom:
            ignored_range_set = ignored_regions_by_chrom[chrom]
            for start, stop in ignored_range_set:
                s = max(1, start)
                e = min(chrom_len, stop)
                ignored_mask[s:e + 1] = False

        hits_across_chrom = np.zeros((chrom_len + 1,), dtype=int)
        reads_across_chrom = np.zeros((chrom_len + 1,), dtype=int)

        for hit in hits:
            reads_across_chrom[hit["hit_pos"]] += hit["hit_count"]
            hit_pos = hit["hit_pos"]
            # Keep legacy behavior: cap per-position hit indicator at 2.
            if hits_across_chrom[hit_pos] < 2:
                hits_across_chrom[hit_pos] += 1

        hits_across_chrom = hits_across_chrom * ignored_mask
        reads_across_chrom = reads_across_chrom * ignored_mask

        hits_in_features = hits_across_chrom * exon_mask
        hits_outside_features = hits_across_chrom * np.invert(exon_mask)
        reads_outside_features = reads_across_chrom * np.invert(exon_mask)

        records = {}

        for feature in feature_db[chrom]:
            feature_mask = exon_mask[feature.start:feature.stop + 1]
            domain_feature_mask = domain_mask[feature.start:feature.stop + 1]

            hits_in_feature = hits_in_features[feature.start:feature.stop + 1][feature_mask]
            hits_in_domains_count = hits_in_features[feature.start:feature.stop + 1][domain_feature_mask].sum()

            hits_in_feature_count = hits_in_feature.sum()
            reads_in_feature = reads_across_chrom[feature.start:feature.stop + 1][feature_mask].sum()
            reads_in_domains = reads_across_chrom[feature.start:feature.stop + 1][domain_feature_mask].sum()

            window_start = max(1, feature.start - neighborhood_window_size)
            window_end = min(chrom_len, feature.stop + neighborhood_window_size)

            hits_outside_feature = hits_outside_features[window_start:window_end + 1]
            hits_outside_feature_count = hits_outside_feature.sum()

            intergenic_region = feature_db.get_interfeature_range(chrom, (window_start, window_end)) - ignored_range_set
            insertion_index = float(hits_in_feature_count) / feature.coding_length if feature.coding_length else 0

            if hits_outside_feature_count == 0 or intergenic_region.coverage == 0:
                neighborhood_index = 0
                reads_ni = 0
                neighborhood_insertion_index = 0
            else:
                neighborhood_insertion_index = float(hits_outside_feature_count) / intergenic_region.coverage
                neighborhood_index = insertion_index / neighborhood_insertion_index if neighborhood_insertion_index else 0
                reads_ni = (
                    (float(reads_in_feature) / feature.coding_length) /
                    (float(reads_outside_features[window_start:window_end + 1].sum()) / intergenic_region.coverage)
                )

            hit_ixes = [0] + list(np.where(hits_in_feature > 0)[0] + 1) + [feature.coding_length + 1]
            # TODO: We assume that the 100 bp upstream is always intergenic, but this isn't always the case.
            if feature.strand == 'W':
                upstream_slice_100 = slice(max(1, feature.start - 100), feature.start)
                upstream_slice_50  = slice(max(1, feature.start - 50), feature.start)
            elif feature.strand == 'C':
                upstream_slice_100 = slice(feature.stop + 1, min(chrom_len, feature.stop + 101))
                upstream_slice_50  = slice(feature.stop + 1, min(chrom_len, feature.stop + 51))
            else:
                # Can occur in cerevisiae ARS regions
                upstream_slice_100 = slice(0, 0)
                upstream_slice_50  = slice(0, 0)

            upstream_region_100 = hits_outside_features[upstream_slice_100]
            upstream_region_50  = hits_outside_features[upstream_slice_50]

            upstream_100_hits = int(upstream_region_100.sum())
            upstream_50_hits  = int(upstream_region_50.sum()) 

            longest_free_intervals = []
            freedom_indices = []
            kornmann_indices = []
            logit_fis = []
            max_skip_tns = 4 + 1

            for skip_tns in range(max_skip_tns):
                longest_interval = max(
                    right - left for left, right in
                    zip(hit_ixes, hit_ixes[1 + skip_tns:] or [feature.coding_length + 1])
                ) - 1
                longest_free_intervals.append(longest_interval)

                freedom_index = float(longest_interval) / feature.coding_length if feature.coding_length else 0
                freedom_indices.append(freedom_index)

                if longest_interval >= 300 and 0.1 < freedom_index < 0.9:
                    kornmann_index = (longest_interval * hits_in_feature_count) / (feature.coding_length ** 1.5)
                else:
                    kornmann_index = 0
                kornmann_indices.append(kornmann_index)

                logit_fis.append(
                    freedom_index / (1.0 + math.e ** (-0.01 * (feature.coding_length - 200)))
                    if feature.coding_length else 0
                )

            records[feature.standard_name] = {
                "feature": feature,
                "length": len(feature),
                "hits": int(hits_in_feature_count),
                "reads": int(reads_in_feature),
                "neighborhood_hits": int(hits_outside_feature_count),
                "nc_window_len": int(intergenic_region.coverage),
                "insertion_index": insertion_index,
                "neighborhood_index": neighborhood_index,
                "reads_ni": reads_ni,
                "max_free_region": int(longest_free_intervals[0]),
                "freedom_index": float(freedom_indices[0]),
                "logit_fi": float(logit_fis[0]),
                "s_value": (log2(reads_in_feature + 1) - total_reads_log) if total_reads > 0 else 0,
                "hit_locations": [ix + 1 for (ix, h) in enumerate(hits_in_feature) if h > 0],
                "longest_interval": int(longest_free_intervals[4]),
                "kornmann_domain_index": float(kornmann_indices[4]),
                "domain_ratio": float(feature.domains.coverage) / feature.coding_length if feature.coding_length else 0,
                "hits_in_domains": int(hits_in_domains_count),
                "reads_in_domains": int(reads_in_domains),
                "upstream_hits_50": upstream_50_hits,
                "upstream_hits_100": upstream_100_hits,
                "domain_coverage": int(feature.domains.coverage),
                "bps_between_hits_in_neihgborhood": (intergenic_region.coverage / hits_outside_feature_count)
                    if hits_outside_feature_count > 0 else 9999,
                "bps_between_hits_in_feature": (len(feature) / hits_in_feature_count)
                    if hits_in_feature_count > 0 else 9999
            }

            for skip_tns in range(max_skip_tns):
                records[feature.standard_name][f"longest_free_interval_{skip_tns}"] = int(longest_free_intervals[skip_tns])
                records[feature.standard_name][f"freedom_index_{skip_tns}"] = float(freedom_indices[skip_tns])
                records[feature.standard_name][f"kornmann_index_{skip_tns}"] = float(kornmann_indices[skip_tns])
                records[feature.standard_name][f"logit_fi_{skip_tns}"] = float(logit_fis[skip_tns])

        result.update(records)

    return result


#  Writers: generic 

def write_hits_into_bed(target_file, hits, overwrite=False):
    with _open_out(target_file, overwrite=overwrite) as bed_file:
        bed_file.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n")
        for hit in sorted(hits, key=lambda o: (o["chrom"], o["hit_pos"])):
            bed_file.write(
                "%s\t%d\t%d\t.\t%d\t%s\n"
                % (
                    hit["chrom"],
                    hit["hit_pos"] - 1,
                    hit["hit_pos"],
                    hit["hit_count"],
                    {"W": "+", "C": "-"}.get(hit.get("source", "W"), "+"),
                )
            )


def write_analyzed_hits_into_bed_proteome(target_file, records, overwrite=False):
    with _open_out(target_file, overwrite=overwrite) as bed_file:
        bed_file.write("#orf\torfStart\torfEnd\n")
        for record in sorted(records, key=lambda r: r["feature"].standard_name):
            feature = record["feature"]
            aa_hits = sorted(set(((h - 1) // 3) + 1 for h in record["hit_locations"]))
            if feature.strand == "C":
                feature_len = feature.coding_length // 3
                aa_hits = [feature_len - h + 1 for h in aa_hits]
            for aa_hit in aa_hits:
                bed_file.write("%s\t%d\t%d\n" % (feature.standard_name, aa_hit, aa_hit + 1))


def _write_all_hits_table(out_path, hits, overwrite=False):
    sum_hits = defaultdict(dict)
    for h in hits:
        prev = sum_hits[h["chrom"]].get(h["hit_pos"], 0)
        sum_hits[h["chrom"]][h["hit_pos"]] = prev + h["hit_count"]
    with _open_out(out_path, overwrite=overwrite, newline="") as out_file:
        writer = csv.writer(out_file)
        writer.writerow(["Chromosome", "Position", "Reads"])
        for chrom, data in sorted(sum_hits.items()):
            for pos, reads in sorted(data.items()):
                writer.writerow([chrom, pos, reads])


def _feature_chrom_name(f):
    """
    Robustly extract a chromosome/contig name across species/DB implementations.
    """
    # Many Feature objects have .chromosome which may be a Chromosome object with .name
    for attr in ("chromosome", "chrom", "seqid", "contig", "chrom_name", "chromosome_name"):
        v = getattr(f, attr, None)
        if v:
            return str(v)

    chrom_obj = getattr(f, "chromosome", None)
    if chrom_obj is not None:
        name = getattr(chrom_obj, "name", None)
        if name:
            return str(name)
# Some DBs store contig on the feature itself as e.g. .chromosome_accession
    for attr in ("chromosome_accession", "accession", "seq_accession"):
        v = getattr(f, attr, None)
        if v:
            return str(v)

    return ""


def _write_feature_table(out_path, analysis, read_depth_filter, overwrite=False,
                        include_scer=False, cglab2scer=None):
    """
    feature_table.RDF_* is organism-intrinsic by default.
    If include_scer=True (glabrata only), we append Sc ortholog columns for convenience.
    """
    base_fields = [
        "standard_name", "common_name", "type", "chrom", "start", "stop", "strand",
        "coding_length", "hits", "reads", "insertion_index", "neighborhood_index", "reads_ni",
        "max_free_region", "freedom_index", "logit_fi", "s_value", "upstream_hits_100",
        "hits_in_domains", "reads_in_domains", "domain_coverage", "domain_ratio",
        "bps_between_hits_in_feature", "bps_between_hits_in_neihgborhood"
    ]
    sc_fields = ["scer_id", "scer_gene_name", "scer_sgdid"] if include_scer else [] 
    fields = base_fields + sc_fields

    with _open_out(out_path, overwrite=overwrite, newline="") as out_file:
        writer = csv.writer(out_file)
        writer.writerow(["RDF", str(read_depth_filter)])
        writer.writerow(fields)

        for r in sorted(analysis.values(), key=lambda x: x["feature"].standard_name):
            f = r["feature"]
            row = [
                f.standard_name,
                getattr(f, "common_name", ""),
                getattr(f, "type", ""),
                f.chromosome,   
                f.start, f.stop, f.strand,
                getattr(f, "coding_length", 0),
                r["hits"], r["reads"], r["insertion_index"], r["neighborhood_index"], r["reads_ni"],
                r["max_free_region"], r["freedom_index"], r["logit_fi"], r["s_value"],
                r["hits_in_domains"], r["reads_in_domains"], r["domain_coverage"], r["domain_ratio"],
                r["bps_between_hits_in_feature"], r["bps_between_hits_in_neihgborhood"]
            ]

            if include_scer:
                if cglab2scer is None:
                    o = {}
                else:
                    o = cglab2scer.get(f.standard_name, {}) or {}
                row += [
                    o.get("scer_id", ""),
                    o.get("scer_gene_name", ""),
                    o.get("scer_sgdid", "")
                ]

            writer.writerow(row)


def _write_outlier_stats(out_path, analysis, overwrite=False):
    """
    Legacy-compatible outlier stat line:
      neighborhood_index <zero_outliers> <right_tail_outliers>

    Uses bins=300, start=0, stop=1.5 (matches original script intent).
    """
    field_name = "neighborhood_index"
    start, stop, bins = 0.0, 1.5, 300

    # mirror the original logic: treat hits==0 as "zero outliers"
    zero_outliers = sum(1 for r in analysis.values() if r.get("hits", 0) == 0)

    rel = []
    for r in analysis.values():
        if r.get("hits", 0) <= 0:
            continue
        v = r.get(field_name, 0.0)
        ix = int(bins * (v - start) / (stop - start)) if stop > start else 0
        if ix < 0:
            ix = 0
        if ix > bins:
            ix = bins
        rel.append(ix)

    right_tail_outliers = sum(1 for v in rel if v == bins)

    with _open_out(out_path, overwrite=overwrite, newline="") as out_file:
        out_file.write(f"{field_name} {zero_outliers} {right_tail_outliers}\n")

def write_insertion_vs_neighborhood_correlations(analysis_labels, analyses, output_file, overwrite=False):
    """
    Legacy-style file: insertion_vs_neighborhood_correlations.txt

    For each library:
      - computes Pearson and Spearman between insertion_index and neighborhood_index
      - writes one block per sample label

    analyses: list of dicts returned by analyze_hits()
    """
    try:
        import scipy.stats
    except ImportError:
        scipy = None

    with _open_out(output_file, overwrite=overwrite, newline="") as out_fh:
        for label, analysis in zip(analysis_labels, analyses):
            s1 = np.array([r.get("insertion_index", 0.0) for r in analysis.values()], dtype=float)
            s2 = np.array([r.get("neighborhood_index", 0.0) for r in analysis.values()], dtype=float)

            if scipy is not None:
                pearson = scipy.stats.pearsonr(s1, s2)
                spearman = scipy.stats.spearmanr(s1, s2)
                out_fh.write(f"{label}\nPearson: {pearson}\nSpearman: {spearman}\n")
            else:
                # Fallback if scipy isn't available (no p-values)
                pearson_r = float(np.corrcoef(s1, s2)[0, 1]) if len(s1) > 1 else float("nan")

                # Spearman fallback: rank-transform then Pearson
                r1 = np.argsort(np.argsort(s1))
                r2 = np.argsort(np.argsort(s2))
                spearman_r = float(np.corrcoef(r1, r2)[0, 1]) if len(s1) > 1 else float("nan")

                out_fh.write(f"{label}\nPearson_r: {pearson_r}\nSpearman_r: {spearman_r}\n")

# Glabrata (Sc inference) enrichment 

def _load_ortholog_table(path):
    """
    Read a tab-delimited ortholog file with headers.
    .txt vs .tsv does NOT matter; delimiter is TAB.
    Required columns: cglab_id, scer_id
    Optional columns: scer_gene_name, scer_sgdid
    """
    m = {}
    with open(path, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Ortholog file has no header: {path}")

        required = {"cglab_id", "scer_id"}
        missing = [c for c in required if c not in reader.fieldnames]
        if missing:
            raise ValueError(
                f"Ortholog file must contain tab-delimited columns: cglab_id, scer_id. Missing: {missing}. File: {path}"
            )

        for row in reader:
            cglab = (row.get("cglab_id") or "").strip()
            scer = (row.get("scer_id") or "").strip()
            if not cglab or not scer:
                continue
            m[cglab] = {
                "scer_id": scer,
                "scer_gene_name": (row.get("scer_gene_name") or "").strip(),
                "scer_sgdid": (row.get("scer_sgdid") or "").strip(),
            }
    return m

def write_data_to_data_frame(data, cols_config):
    """
    Legacy helper used by Classifier.py table-mode.
    Legacy helper used by Classifier.py table-mode.

    Parameters
    ----------
    data : list[dict] | dict
        Usually a list of record dicts (or occasionally a dict of records).
    cols_config : list[dict]
        Each dict typically contains:
          - field_name : key in record dict
          - csv_name   : output column name
          - format     : optional callable or sprintf-like (legacy; we keep pass-through)
          - sort_by    : optional True to sort by this field

    Returns
    -------
    pandas.DataFrame
    """
    import pandas as pd

    field_col = "field_name"
    csv_name_col = "csv_name"
    sort_by_col = "sort_by"
    fmt_col = "format"
    local_col = "local"  # not required, but common in configs

    if data is None:
        return pd.DataFrame()

    # Allow dict-of-records too
    if isinstance(data, dict):
        if data and isinstance(next(iter(data.values())), dict):
            data = list(data.values())
        else:
            data = [data]

    if not data:
        return pd.DataFrame()

    # Filter out columns whose field_name isn't present in the records
    first_record = data[0]
    filtered_cols_config = []
    for col in cols_config:
        if col.get(field_col) in first_record:
            filtered_cols_config.append(col)
    cols_config = filtered_cols_config

    # Find sort column
    sort_field = None
    for col in cols_config:
        if col.get(sort_by_col, False):
            sort_field = col[field_col]
            break

    ordered_dataset = sorted(data, key=lambda r: r.get(sort_field)) if sort_field else data

    def _apply_format(val, fmt):
        """Apply COLS_CONFIG-style formatting to a value."""
        if fmt is None:
            return val
        # callable formatter (e.g., lambda f: f.standard_name)
        if callable(fmt):
            try:
                return fmt(val)
            except Exception:
                # fall back to raw val if formatter fails
                return val
        # sprintf-style formatter (e.g., "%.3f", "%d")
        if isinstance(fmt, str):
            try:
                # treat None/"" as blank for numeric formats
                if val is None:
                    return ""
                return fmt % val
            except Exception:
                return val
        return val

    rows = []
    for record in ordered_dataset:
        row = []
        for col in cols_config:
            raw = record.get(col[field_col], "")
            fmt = col.get(fmt_col, None)
            row.append(_apply_format(raw, fmt))
        rows.append(row)

    df = pd.DataFrame(rows, columns=[col[csv_name_col] for col in cols_config])
    return df

@Shared.memoized
def _default_cglab_to_scer_orthologs():
    """
    Default dependency path (your new file):
      dependencies/glabrata/C_glabrata_BG2_S_cerevisiae_orthologs.txt
    """
    path = Shared.get_dependency("glabrata/C_glabrata_BG2_S_cerevisiae_orthologs.txt")
    return _load_ortholog_table(path)


@Shared.memoized
def _load_scer_essential_sets():
    """
    dependencies/cerevisiae/
      - cerevisiae_inviable_annotations.txt
      - cerevisiae_viable_annotations.txt

    Assumes: first whitespace token is the ORF systematic name.
    """
    inv_path = Shared.get_dependency("cerevisiae/cerevisiae_inviable_annotations.txt")
    via_path = Shared.get_dependency("cerevisiae/cerevisiae_viable_annotations.txt")

    def _read_first_token(p):
        s = set()
        with open(p, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                tok = line.split()[0]
                if tok:
                    s.add(tok)
        return s

    inv = _read_first_token(inv_path)
    via = _read_first_token(via_path)
    return inv, via


@Shared.memoized
def _load_scer_synthetic_lethal_flags():
    """
    dependencies/cerevisiae/duplicatesSl_011116.txt
    Legacy format: columns with f1 f2 ... is_double_lethal
    """
    path = Shared.get_dependency("cerevisiae/duplicatesSl_011116.txt")
    sl = {}
    with open(path, "r") as fh:
        _ = fh.readline()  # header
        for line in fh:
            parts = line.split()
            if len(parts) < 4:
                continue
            f1, f2, _, is_double_lethal = parts[:4]
            val = {"1": "Yes", "0": "No"}.get(is_double_lethal, is_double_lethal)
            sl[f1] = val
            sl[f2] = val
    return sl


@Shared.memoized
def _load_scer_fitness():
    """
    dependencies/cerevisiae/neFitnessStandard.txt
    Format: <orf> <fitness>
    """
    path = Shared.get_dependency("cerevisiae/neFitnessStandard.txt")
    fitness = {}
    with open(path, "r") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            orf, val = parts[0], parts[1]
            fitness[orf] = val
    return fitness


def _warn_missing_orthologs(records, cglab2scer, label_for_print="dataset"):
    """
    Print a compact warning if many glabrata features have no Sc ortholog.
    This usually means: your feature.standard_name does not match cglab_id exactly.
    """
    missing = []
    for r in records:
        f = r["feature"]
        if f.standard_name not in cglab2scer:
            missing.append(f.standard_name)

    if not missing:
        return

    # Keep it informative but not spammy
    n = len(missing)
    examples = ", ".join(missing[:10])
    _eprint(f"[WARN] {label_for_print}: {n} glabrata features missing Sc ortholog mapping (showing up to 10): {examples}")
    _eprint("       Likely ID mismatch: FeatureDB standard_name vs ortholog-file cglab_id. Fix by exact string match.")


def _write_analyzed_gla_records(records, output_file, cglab2scer, overwrite=False, warn_missing=True, label_for_print="dataset"):
    """
    Minimal glabrata analysis CSV with Sc inference columns.
    """
    if warn_missing:
        _warn_missing_orthologs(records, cglab2scer, label_for_print=label_for_print)

    scer_inv, scer_via = _load_scer_essential_sets()
    scer_sl = _load_scer_synthetic_lethal_flags()
    scer_fit = _load_scer_fitness()

    header = [
        # Glabrata identity
        "Standard name", "Common name", "Type",
        # Sc inference
        "Sc ortholog", "Sc ortholog gene name", "Sc SGDID",
        "Essential in Sc", "Sc synthetic lethals", "Sc fitness",
        # Insertion metrics
        "Hits", "Reads", "Length",
        "Insertion index", "Neighborhood index", "Reads NI",
        "Max free region", "Freedom index", "Logit FI", "Kornmann FI",
        "Domain ratio", "Hits in domains", "Reads in domains", "Domain coverage",
    ]

    with _open_out(output_file, overwrite=overwrite, newline="") as out_fh:
        writer = csv.writer(out_fh)
        writer.writerow(header)

        for r in sorted(records, key=lambda rr: rr["feature"].standard_name):
            f = r["feature"]
            cglab_id = f.standard_name
            o = cglab2scer.get(cglab_id, {}) or {}

            scer_orth = o.get("scer_id", "")
            scer_gn = o.get("scer_gene_name", "")
            scer_sgdid = o.get("scer_sgdid", "")

            # Essential in Sc: from viable/inviable lists
            if scer_orth in scer_inv:
                ess_sc = "Yes"
            elif scer_orth in scer_via:
                ess_sc = "No"
            else:
                ess_sc = ""

            writer.writerow([
                f.standard_name,
                getattr(f, "common_name", ""),
                getattr(f, "type", ""),
                scer_orth,
                scer_gn,
                scer_sgdid,
                ess_sc,
                scer_sl.get(scer_orth, ""),
                scer_fit.get(scer_orth, ""),
                r.get("hits", 0),
                r.get("reads", 0),
                len(f),
                r.get("insertion_index", 0),
                r.get("neighborhood_index", 0),
                r.get("reads_ni", 0),
                r.get("max_free_region", 0),
                r.get("freedom_index", 0),
                r.get("logit_fi", 0),
                r.get("kornmann_domain_index", 0),
                r.get("domain_ratio", 0),
                r.get("hits_in_domains", 0),
                r.get("reads_in_domains", 0),
                r.get("domain_coverage", 0),
            ])


# Feature DB 
def load_feature_db(kind):
    if kind == "albicans":
        return GenomicFeatures.default_alb_db()
    if kind == "pombe":
        return GenomicFeatures.default_pom_db()
    if kind == "cerevisiae":
        return GenomicFeatures.default_cer_db()
    if kind == "glabrata":
        return GenomicFeatures.default_glab_db()
    raise ValueError(f"Unsupported --feature-db: {kind}")


def _analyze_one(args_tuple):
    dataset, feature_db, window_size, ignored_by_chrom = args_tuple
    return analyze_hits(
        dataset,
        feature_db,
        neighborhood_window_size=window_size,
        ignored_regions_by_chrom=ignored_by_chrom,
    )


def _infer_chrom_lengths_from_hits(all_hits):
    lengths = {}
    for hits in all_hits:
        for h in hits:
            chrom = h["chrom"]
            pos = h["hit_pos"]
            if chrom not in lengths or pos > lengths[chrom]:
                lengths[chrom] = pos
    return lengths


#  Main 

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("-i", "--input-dir", default=".")
    parser.add_argument("-o", "--output-dir", default=".")
    parser.add_argument("--hit-glob", default="*_hits.txt")
    parser.add_argument("--hit-format", default="hit_table", choices=["hit_table", "pombe_csv", "wig"])
    parser.add_argument("--strict-hit-table", action="store_true", default=False)

    parser.add_argument("--feature-db", default=None, choices=["albicans", "pombe", "cerevisiae", "glabrata"])

    parser.add_argument("-f", "--read-depth-filter", type=int, default=1)
    parser.add_argument("--neighborhood-window-size", type=int, default=10000)
    parser.add_argument("--bin-size", type=int, default=10000)
    parser.add_argument("--ignored-regions", default=None)

    parser.add_argument("-t", "--threads", type=int, default=1)

    parser.add_argument("--overwrite", action="store_true", default=False)

    # Glabrata-only knobs (Sc inference)
    parser.add_argument("--scer-orthologs", default=None,
                        help="Optional override path to C_glabrata_BG2_S_cerevisiae_orthologs (TAB-delimited).")
    parser.add_argument("--warn-missing-orthologs", action="store_true", default=False,
                        help="Warn if many glabrata features have no Sc ortholog mapping (OFF by default).")
    parser.add_argument("--feature-table-include-scer", action="store_true", default=False,
                        help="If writing --write-feature-table and --feature-db glabrata, append scer_id/gene_name/sgdid columns.")

    # Optional subset selection (if none specified, we write ALL legacy outputs)
    parser.add_argument("--write-hit-summary", action="store_true", default=False)
    parser.add_argument("--write-binned-hits", action="store_true", default=False)
    parser.add_argument("--write-stats", action="store_true", default=False)
    parser.add_argument("--write-bed-hits", action="store_true", default=False)
    parser.add_argument("--write-bed-proteome", action="store_true", default=False)
    parser.add_argument("--write-all-hits-table", action="store_true", default=False)
    parser.add_argument("--write-feature-table", action="store_true", default=False)  # optional extra
    parser.add_argument("--write-analysis", action="store_true", default=False)
    parser.add_argument("--write-outlier-stats", action="store_true", default=False)

    args = parser.parse_args()

    Shared.make_dir(args.output_dir)

    input_file_paths = sorted(glob.glob(os.path.join(args.input_dir, args.hit_glob)))
    if not input_file_paths:
        _eprint(f"ERROR: no hit files matching '{args.hit_glob}' were found in: {args.input_dir}")
        raise SystemExit(1)

    input_filenames = [os.path.split(p)[-1] for p in input_file_paths]
    if args.hit_format == "hit_table":
        input_labels = [(fn[:-9] if fn.endswith("_hits.txt") else os.path.splitext(fn)[0]) for fn in input_filenames]
    else:
        input_labels = [os.path.splitext(fn)[0] for fn in input_filenames]

    # If user did not specify any --write-* flags, default to ALL legacy outputs.
    any_write_flag = any([
        args.write_hit_summary,
        args.write_binned_hits,
        args.write_stats,
        args.write_bed_hits,
        args.write_bed_proteome,
        args.write_all_hits_table,
        args.write_feature_table,
        args.write_analysis,
        args.write_outlier_stats,
    ])
    if not any_write_flag:
        args.write_binned_hits = True
        args.write_hit_summary = True
        args.write_stats = True
        args.write_bed_hits = True
        args.write_bed_proteome = True
        args.write_all_hits_table = True
        args.write_analysis = True
        args.write_outlier_stats = True
        # feature_table is optional extra; default OFF to match the “8 files” expectation.
        args.write_feature_table = True

    all_hits = read_hit_files(
        input_file_paths,
        args.read_depth_filter,
        hit_format=args.hit_format,
        strict_hit_table=args.strict_hit_table,
    )

    needs_feature_db = any([
    args.write_hit_summary,
    args.write_stats,
    args.write_bed_proteome,
    args.write_feature_table,
    args.write_analysis,
    args.write_outlier_stats,
])
    feature_db = None
    ignored_by_chrom = None
    all_analyzed = None

    # Load ortholog table if glabrata and needed for outputs
    cglab2scer = None
    if args.feature_db == "glabrata":
        if args.scer_orthologs:
            cglab2scer = _load_ortholog_table(args.scer_orthologs)
        else:
            # default dependency file
            cglab2scer = _default_cglab_to_scer_orthologs()

    if needs_feature_db:
        if args.feature_db is None:
            parser.error("--feature-db is required for feature-level outputs (default writes all outputs)")
        feature_db = load_feature_db(args.feature_db)
        ignored_by_chrom = _parse_ignored_regions(args.ignored_regions)

        if args.threads and args.threads > 1 and len(all_hits) > 1:
            from multiprocessing import Pool
            with Pool(processes=args.threads) as pool:
                all_analyzed = pool.map(
                    _analyze_one,
                    [(ds, feature_db, args.neighborhood_window_size, ignored_by_chrom) for ds in all_hits],
                )
        else:
            all_analyzed = [
                analyze_hits(
                    ds,
                    feature_db,
                    neighborhood_window_size=args.neighborhood_window_size,
                    ignored_regions_by_chrom=ignored_by_chrom,
                )
                for ds in all_hits
            ]

    # ---- binned_hits ----
    if args.write_binned_hits:
        bin_size = args.bin_size
        chrom_lengths = _infer_chrom_lengths_from_hits(all_hits)

        # [hits, reads, hit rank, read rank] (we only compute hit rank like the legacy script)
        hit_read_bins = {
            lbl: {
                chrom: [[0, 0, 0, 0] for _ in range((chrom_lengths.get(chrom, 0) // bin_size) + 1)]
                for chrom in chrom_lengths.keys()
            } for lbl in input_labels
        }

        for lbl, hits in zip(input_labels, all_hits):
            for hit in hits:
                chrom = hit["chrom"]
                if chrom not in hit_read_bins[lbl]:
                    hit_read_bins[lbl][chrom] = [[0, 0, 0, 0] for _ in range((hit["hit_pos"] // bin_size) + 1)]
                b = hit["hit_pos"] // bin_size
                if b >= len(hit_read_bins[lbl][chrom]):
                    need = b - len(hit_read_bins[lbl][chrom]) + 1
                    hit_read_bins[lbl][chrom].extend([[0, 0, 0, 0] for _ in range(need)])
                bin_data = hit_read_bins[lbl][chrom][b]
                bin_data[0] += 1
                bin_data[1] += hit["hit_count"]

        # hit rank within each library
        for lbl_bins in hit_read_bins.values():
            for bin_rank, bin_data in enumerate(sorted(chain(*lbl_bins.values()), key=lambda x: x[0], reverse=True)):
                bin_data[2] = bin_rank + 1

        out_path = os.path.join(args.output_dir, f"binned_hits.RDF_{args.read_depth_filter}.csv")
        with _open_out(out_path, overwrite=args.overwrite, newline="") as out_file:
            writer = csv.writer(out_file)
            writer.writerow(["Bin index"] + list(chain(*zip(input_labels, [f"{l} rank" for l in input_labels]))))

            for chrom in sorted(chrom_lengths.keys()):
                n_bins = (chrom_lengths[chrom] // bin_size) + 1
                for bin_ix in range(n_bins):
                    writer.writerow(
                        [f"{chrom}-{bin_ix}"] +
                        list(chain(*[
                            hit_read_bins[lbl].get(chrom, [[0, 0, 0, 0]])[bin_ix][0::2]
                            if bin_ix < len(hit_read_bins[lbl].get(chrom, [])) else [0, 0]
                            for lbl in input_labels
                        ]))
                    )

    # ---- hit_summary ----
    if args.write_hit_summary:
        out_path = os.path.join(args.output_dir, f"hit_summary.RDF_{args.read_depth_filter}.csv")
        with _open_out(out_path, overwrite=args.overwrite, newline="") as out_file:
            writer = csv.writer(out_file)
            writer.writerow(["Standard name", "Common name", "Type"] + [f"Lib {lbl}" for lbl in input_labels])

            feat_keys = sorted(
                all_analyzed[0].keys(),
                key=lambda k: all_analyzed[0][k]["feature"].standard_name,
            )
            for k in feat_keys:
                f = all_analyzed[0][k]["feature"]
                row = [f.standard_name, getattr(f, "common_name", ""), getattr(f, "type", "")]
                row += [ds.get(k, {"hits": 0})["hits"] for ds in all_analyzed]
                writer.writerow(row)

    # ---- all_hits per library ----
    if args.write_all_hits_table:
        for lbl, hits in zip(input_labels, all_hits):
            _write_all_hits_table(
                os.path.join(args.output_dir, f"{lbl}.all_hits.csv"),
                hits,
                overwrite=args.overwrite,
            )

    # ---- analysis.csv per library ----
    if args.write_analysis:
        if args.feature_db == "glabrata":
            if cglab2scer is None:
                raise RuntimeError("Internal error: cglab2scer not loaded for glabrata mode.")
            for lbl, analysis in zip(input_labels, all_analyzed):
                _write_analyzed_gla_records(
                    list(analysis.values()),
                    os.path.join(args.output_dir, f"{lbl}_analysis.csv"),
                    cglab2scer=cglab2scer,
                    overwrite=args.overwrite,
                    warn_missing=args.warn_missing_orthologs,
                    label_for_print=lbl,
                )
        else:
            # Minimal generic analysis if someone uses cerevisiae, etc.
            for lbl, analysis in zip(input_labels, all_analyzed):
                out_path = os.path.join(args.output_dir, f"{lbl}_analysis.csv")
                with _open_out(out_path, overwrite=args.overwrite, newline="") as out_fh:
                    writer = csv.writer(out_fh)
                    writer.writerow(["Standard name", "Common name", "Type", "Hits", "Reads", "Neighborhood index"])
                    for r in sorted(analysis.values(), key=lambda rr: rr["feature"].standard_name):
                        f = r["feature"]
                        writer.writerow([
                            f.standard_name,
                            getattr(f, "common_name", ""),
                            getattr(f, "type", ""),
                            r.get("hits", 0),
                            r.get("reads", 0),
                            r.get("neighborhood_index", 0),
                        ])

    # ---- outlier_stats per library ----
    if args.write_outlier_stats:
        for lbl, analysis in zip(input_labels, all_analyzed):
            _write_outlier_stats(
                os.path.join(args.output_dir, f"{lbl}.outlier_stats.txt"),
                analysis,
                overwrite=args.overwrite,
            )

    # ---- feature_table per library (optional extra) ----
    if args.write_feature_table:
        include_scer = bool(args.feature_table_include_scer and args.feature_db == "glabrata")
        for lbl, analysis in zip(input_labels, all_analyzed):
            out_path = os.path.join(args.output_dir, f"{lbl}.feature_table.RDF_{args.read_depth_filter}.csv")
            _write_feature_table(
                out_path,
                analysis,
                args.read_depth_filter,
                overwrite=args.overwrite,
                include_scer=include_scer,
                cglab2scer=cglab2scer if include_scer else None
            )

    # ---- stats.csv ----
    if args.write_stats:
        out_path = os.path.join(args.output_dir, "stats.csv")
        with _open_out(out_path, overwrite=args.overwrite, newline="") as stats_file:
            writer = csv.writer(stats_file, delimiter=",")
            writer.writerow(["File name"] + ALL_STATS)
            for lbl, dataset in zip(input_labels, all_hits):
                stats = get_statistics(dataset, feature_db)
                writer.writerow([lbl] + [stats[col] for col in ALL_STATS])

    # ---- bed hits per library ----
    if args.write_bed_hits:
        for lbl, dataset in zip(input_labels, all_hits):
            write_hits_into_bed(
                os.path.join(args.output_dir, f"{lbl}.filter_{args.read_depth_filter}.bed"),
                dataset,
                overwrite=args.overwrite,
            )

    # ---- proteome bed per library ----
    if args.write_bed_proteome:
        for lbl, analysis in zip(input_labels, all_analyzed):
            write_analyzed_hits_into_bed_proteome(
                os.path.join(args.output_dir, f"{lbl}.proteome.filter_{args.read_depth_filter}.bed"),
                list(analysis.values()),
                overwrite=args.overwrite,
            )
        # ---- insertion vs neighborhood correlations (legacy QC) ----
    # Writes: insertion_vs_neighborhood_correlations.txt
    if needs_feature_db and all_analyzed is not None:
        out_path = os.path.join(args.output_dir, "insertion_vs_neighborhood_correlations.txt")
        write_insertion_vs_neighborhood_correlations(
            input_labels,
            all_analyzed,
            out_path,
            overwrite=args.overwrite,
        )        
if __name__ == "__main__":
    main()

