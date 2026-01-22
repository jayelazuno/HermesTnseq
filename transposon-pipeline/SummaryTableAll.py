#!/usr/bin/env python3
"""
SummaryTable_generic.py (01/21/2026)
Generic (non-albicans) adaptation of Berman-lab SummaryTable.py for Tn-seq.

Inputs:
- *_hits.txt (from CreateHitFile; tab-delimited with header)
- genome FASTA (for chromosome lengths; uses .fai if present)
- annotation GFF3 (for CDS/exon coordinates; gene_id or Name used as feature name)
Optional:
- ignored/unmappable BED (0-based half-open) to mask problematic regions

Outputs (in --output-dir):
- <sample>_analysis.csv : per-gene metrics table
- stats.csv            : per-sample global stats
- <sample>.filter_N.bed : hit sites as BED6
- hit_summary.RDF_N.csv (if multiple samples): nominal hits per gene
- binned_hits.RDF_N.csv (if multiple samples): bin hit ranks per chromosome

Notes:
- Hit position from hit file is assumed 1-based.
- GFF3 coordinates are 1-based inclusive.
- Internally we use arrays indexed 0..chrom_len; position p uses index p.
"""

import os, csv, glob, math, argparse, re
from collections import defaultdict
import numpy as np
import pandas as pd
import scipy.stats

# If the repo has RangeSet.py, keep using it (it’s handy for coverage).
from RangeSet import RangeSet

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------- IO helpers --------------------------------

def read_fasta_lengths(fasta_path: str) -> dict:
    """Return dict: chrom -> length. Uses .fai if present."""
    lengths = {}
    fai = fasta_path + ".fai"
    if os.path.exists(fai):
        with open(fai) as f:
            for line in f:
                chrom, L = line.split("\t")[0], int(line.split("\t")[1])
                lengths[chrom] = L
        return lengths

    # fallback: scan fasta
    name = None
    seq_len = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = seq_len
                name = line[1:].strip().split()[0]
                seq_len = 0
            else:
                seq_len += len(line.strip())
    if name is not None:
        lengths[name] = seq_len
    return lengths

def read_bed_ignored(bed_path: str, chrom_lengths: dict) -> dict:
    """
    Read 0-based half-open BED and return dict chrom -> RangeSet (1-based inclusive)
    so it matches the original masking logic (start:stop inclusive).
    """
    ignored = {c: RangeSet() for c in chrom_lengths.keys()}
    if not bed_path:
        return ignored
    with open(bed_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start0, end0 = line.split("\t")[:3]
            if chrom not in chrom_lengths:
                continue
            s0 = int(start0)
            e0 = int(end0)
            # BED [s0, e0) -> 1-based inclusive [s0+1, e0]
            s1 = s0 + 1
            e1 = e0
            if e1 < s1:
                continue
            ignored[chrom].add((max(1, s1), min(chrom_lengths[chrom], e1)))
    return ignored

# ----------------------------- Feature DB --------------------------------

class Feature:
    """
    Minimal feature object compatible with the analysis needs.
    Exons are stored as (start, stop) in 1-based inclusive coordinates.
    """
    def __init__(self, standard_name, chromosome, start, stop, strand, exons):
        self.standard_name = standard_name
        self.common_name = standard_name
        self.name = standard_name
        self.chromosome = chromosome
        self.start = int(start)
        self.stop = int(stop)
        # Use W/C like the original pipeline
        self.strand = "W" if strand == "+" else "C" if strand == "-" else strand
        self.exons = exons[:] if exons else [(self.start, self.stop)]
        self.type = "ORF"
        self.description = ""
        self.is_orf = True
        self.cerevisiae_orthologs = set()
        self.domains = RangeSet()  # empty by default

    def __len__(self):
        return abs(self.stop - self.start) + 1

    @property
    def coding_length(self):
        return sum((e2 - e1 + 1) for e1, e2 in self.exons)

def parse_gff3_features(gff_path: str, chrom_lengths: dict,
                        feature_types=("CDS",),
                        name_attr_candidates=("gene", "gene_id", "Name", "ID", "locus_tag")):
    """
    Parse GFF3 and return dict chrom -> list[Feature].
    Strategy:
    - Group CDS rows by Parent or by gene-like attribute
    - Build per-gene exon list from CDS intervals
    """
    gene_exons = defaultdict(lambda: defaultdict(list))  # chrom -> gene_key -> [(s,e),...]
    gene_strand = defaultdict(dict)  # chrom -> gene_key -> '+/-'
    gene_bounds = defaultdict(dict)  # chrom -> gene_key -> (min_start, max_stop)

    def parse_attrs(attr_str: str) -> dict:
        d = {}
        for part in attr_str.strip().split(";"):
            if not part:
                continue
            if "=" in part:
                k, v = part.split("=", 1)
                d[k] = v
        return d

    with open(gff_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _src, ftype, start, end, _score, strand, _phase, attrs = cols
            if chrom not in chrom_lengths:
                continue
            if ftype not in feature_types:
                continue

            s = int(start); e = int(end)
            if e < s:
                s, e = e, s

            ad = parse_attrs(attrs)
            gene_key = ad.get("Parent") or ad.get("gene_id") or ad.get("ID")
            if not gene_key:
                for cand in name_attr_candidates:
                    if cand in ad:
                        gene_key = ad[cand]
                        break
            if not gene_key:
                # last resort: make something stable-ish
                gene_key = f"{chrom}:{s}-{e}:{strand}"

            gene_exons[chrom][gene_key].append((s, e))
            gene_strand[chrom][gene_key] = strand
            if gene_key not in gene_bounds[chrom]:
                gene_bounds[chrom][gene_key] = (s, e)
            else:
                mn, mx = gene_bounds[chrom][gene_key]
                gene_bounds[chrom][gene_key] = (min(mn, s), max(mx, e))

    features_by_chrom = {}
    for chrom in gene_exons:
        feats = []
        for gene_key, exons in gene_exons[chrom].items():
            exons = sorted(exons)
            mn, mx = gene_bounds[chrom][gene_key]
            strand = gene_strand[chrom].get(gene_key, "+")
            # If Parent=transcript, clean a bit (optional)
            std_name = gene_key.split(",")[0]
            feats.append(Feature(std_name, chrom, mn, mx, strand, exons))
        # sort by coordinate
        feats.sort(key=lambda f: (f.start, f.stop, f.standard_name))
        features_by_chrom[chrom] = feats

    return features_by_chrom

class FeatureDB:
    def __init__(self, chrom_lengths: dict, features_by_chrom: dict, ignored_regions: dict):
        self.chrom_lengths = chrom_lengths
        self.features_by_chrom = features_by_chrom
        self.ignored_regions = ignored_regions or {c: RangeSet() for c in chrom_lengths}

        # Index by name
        self._by_name = {}
        for chrom, feats in self.features_by_chrom.items():
            for feat in feats:
                self._by_name[feat.standard_name] = feat

    def __iter__(self):
        # yield "chrom objects" minimal: (name, length, feature list)
        for chrom in sorted(self.chrom_lengths.keys()):
            yield chrom

    def __getitem__(self, chrom):
        return self.features_by_chrom.get(chrom, [])

    def get_all_features(self):
        return list(self._by_name.values())

    def get_feature_by_name(self, name):
        return self._by_name.get(name)

    def get_features_at_location(self, chrom, pos_1based):
        feats = self.features_by_chrom.get(chrom, [])
        out = []
        for f in feats:
            if f.start <= pos_1based <= f.stop:
                # check exon membership
                for e1, e2 in f.exons:
                    if e1 <= pos_1based <= e2:
                        out.append(f)
                        break
        return out

    def get_interfeature_range(self, chrom, window):
        """
        Return a RangeSet of 'interfeature' positions in [window_start, window_end] (1-based inclusive),
        defined here as positions NOT covered by any exon of any feature.
        """
        ws, we = window
        ws = max(1, ws)
        we = min(self.chrom_lengths.get(chrom, we), we)
        if we < ws:
            return RangeSet()

        covered = RangeSet()
        for f in self.features_by_chrom.get(chrom, []):
            # If feature doesn't overlap window, skip
            if f.stop < ws or f.start > we:
                continue
            for e1, e2 in f.exons:
                s = max(ws, e1)
                t = min(we, e2)
                if t >= s:
                    covered.add((s, t))

        window_rs = RangeSet()
        window_rs.add((ws, we))
        return window_rs - covered

# ----------------------------- Hit parsing --------------------------------

def read_hit_file(filename, read_depth_filter=1):
    result = []
    with open(filename, "r") as in_file:
        next(in_file)  # header
        for line in in_file:
            if not line.strip():
                continue
            chrom, source, _up_feature_type, up_feature_name, up_gene_name, \
                   up_feature_dist, _down_feature_type, down_feature_name, \
                   down_gene_name, down_feature_dist, ig_type, \
                   hit_pos, hit_count = line.rstrip("\n").split("\t")

            hit_count = int(hit_count)
            if hit_count < read_depth_filter:
                continue

            # ORF hits: take the ORF gene name if present; else feature name
            if "ORF" in ig_type:
                gene_name = up_gene_name if up_gene_name != "nan" else up_feature_name
            else:
                gene_name = "nan"

            result.append({
                "chrom": chrom,
                "source": source,  # W/C
                "hit_pos": int(hit_pos),  # 1-based
                "hit_count": hit_count,
                "gene_name": gene_name
            })
    return result

def read_hit_files(files, read_depth_filter=1):
    return [read_hit_file(f, read_depth_filter) for f in files]

# ----------------------------- Statistics ---------------------------------

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

def get_statistics(dataset, feature_db):
    result = {}
    result[TOTAL_HITS] = len(dataset)
    result[TOTAL_READS] = sum(obj["hit_count"] for obj in dataset)
    result[AVG_READS_PER_HIT] = (result[TOTAL_READS] / result[TOTAL_HITS]) if result[TOTAL_HITS] else 0
    result[ORF_HITS] = 0
    result[ANNOTATED_FEATURE_HITS] = 0
    result[INTERGENIC_HITS] = 0
    result[READS_IN_FEATURES] = 0

    features_hit = set()
    for hit in dataset:
        feats = feature_db.get_features_at_location(hit["chrom"], hit["hit_pos"])
        if not feats:
            continue
        result[ANNOTATED_FEATURE_HITS] += 1
        result[READS_IN_FEATURES] += hit["hit_count"]
        features_hit.update(set(f.standard_name for f in feats))
        for f in feats:
            if f.is_orf:
                result[ORF_HITS] += 1
                break

    result[INTERGENIC_HITS] = result[TOTAL_HITS] - result[ANNOTATED_FEATURE_HITS]
    result[FEATURES_HIT] = len(features_hit)

    if result[TOTAL_HITS]:
        result[PER_INTERGENIC_HITS] = "%.2f%%" % (result[INTERGENIC_HITS] * 100.0 / result[TOTAL_HITS])
        result[PER_ANNOTATED_FEATURE_HITS] = "%.2f%%" % (result[ANNOTATED_FEATURE_HITS] * 100.0 / result[TOTAL_HITS])
    else:
        result[PER_INTERGENIC_HITS] = "0.00%"
        result[PER_ANNOTATED_FEATURE_HITS] = "0.00%"

    all_feats = feature_db.get_all_features()
    denom_feats = len(all_feats) if all_feats else 1
    result[PER_FEATURES_HIT] = "%.2f%%" % (result[FEATURES_HIT] * 100.0 / denom_feats) if denom_feats else "0.00%"

    result[AVG_HITS_IN_FEATURE] = "%.1f" % (result[ANNOTATED_FEATURE_HITS] * 1.0 / result[FEATURES_HIT]) if result[FEATURES_HIT] else "0.0"
    result[AVG_READS_IN_FEATURE] = (result[READS_IN_FEATURES] / result[FEATURES_HIT]) if result[FEATURES_HIT] else 0
    result[AVG_READS_IN_FEATURE_HIT] = (result[READS_IN_FEATURES] / result[ANNOTATED_FEATURE_HITS]) if result[ANNOTATED_FEATURE_HITS] else 0

    return result

# ----------------------------- Core analysis ------------------------------

def analyze_hits(dataset, feature_db, neighborhood_window_size=10000):
    log2 = lambda v: math.log(v, 2)

    result = {}
    total_reads = sum(h["hit_count"] for h in dataset)
    total_reads_log = log2(total_reads + 1)

    chroms = set(h["chrom"] for h in dataset)
    for chrom in chroms:
        hits = [h for h in dataset if h["chrom"] == chrom]
        chrom_len = feature_db.chrom_lengths.get(chrom, 0)
        if chrom_len <= 0:
            continue

        # Exon mask and domain mask
        exon_mask = np.zeros((chrom_len + 1,), dtype=bool)
        domain_mask = np.zeros((chrom_len + 1,), dtype=bool)

        for feature in feature_db[chrom]:
            for e1, e2 in feature.exons:
                exon_mask[e1:e2+1] = True
            for d1, d2 in feature.domains:
                domain_mask[d1:d2+1] = True

        # ignored mask (True = keep)
        ignored_mask = np.ones((chrom_len + 1,), dtype=bool)
        ignored_range_set = feature_db.ignored_regions.get(chrom, RangeSet())
        for start, stop in ignored_range_set:
            ignored_mask[start:stop+1] = False

        hits_across_chrom = np.zeros((chrom_len + 1,), dtype=np.int32)
        reads_across_chrom = np.zeros((chrom_len + 1,), dtype=np.int32)

        for hit in hits:
            pos = hit["hit_pos"]
            if pos < 1 or pos > chrom_len:
                continue
            reads_across_chrom[pos] += hit["hit_count"]
            # cap hit-count-per-site at 2 (kept from original)
            if hits_across_chrom[pos] < 2:
                hits_across_chrom[pos] += 1

        hits_across_chrom = hits_across_chrom * ignored_mask
        reads_across_chrom = reads_across_chrom * ignored_mask

        hits_in_features = hits_across_chrom * exon_mask
        hits_outside_features = hits_across_chrom * (~exon_mask)
        reads_outside_features = reads_across_chrom * (~exon_mask)

        records = {}

        for feature in feature_db[chrom]:
            # Feature-local masks
            f_start, f_stop = feature.start, feature.stop
            if f_start < 1 or f_stop > chrom_len or f_stop < f_start:
                continue

            local_exon_mask = exon_mask[f_start:f_stop+1]
            local_domain_mask = domain_mask[f_start:f_stop+1]

            local_hits_all = hits_in_features[f_start:f_stop+1]
            local_hits_exonic = local_hits_all[local_exon_mask]
            hits_in_feature_count = int(local_hits_exonic.sum())
            reads_in_feature = int(reads_across_chrom[f_start:f_stop+1][local_exon_mask].sum())

            hits_in_domains_count = int(local_hits_all[local_domain_mask].sum())
            reads_in_domains = int(reads_across_chrom[f_start:f_stop+1][local_domain_mask].sum())

            window_start = max(1, f_start - neighborhood_window_size)
            window_end = min(chrom_len, f_stop + neighborhood_window_size)

            hits_outside_feature = hits_outside_features[window_start:window_end+1]
            hits_outside_feature_count = int(hits_outside_feature.sum())

            intergenic_region = feature_db.get_interfeature_range(chrom, (window_start, window_end)) - ignored_range_set
            insertion_index = float(hits_in_feature_count) / feature.coding_length if feature.coding_length else 0.0

            if hits_outside_feature_count == 0 or intergenic_region.coverage == 0:
                neighborhood_index = 0.0
                reads_ni = 0.0
                neighborhood_insertion_index = 0.0
            else:
                neighborhood_insertion_index = float(hits_outside_feature_count) / intergenic_region.coverage
                neighborhood_index = insertion_index / neighborhood_insertion_index if neighborhood_insertion_index else 0.0
                reads_ni = (
                    (float(reads_in_feature) / feature.coding_length) /
                    (float(reads_outside_features[window_start:window_end+1].sum()) / intergenic_region.coverage)
                ) if feature.coding_length else 0.0

            # upstream hits
            if feature.strand == "W":
                upstream_slice_100 = slice(max(1, f_start - 100), f_start)
                upstream_slice_50 = slice(max(1, f_start - 50), f_start)
            else:
                upstream_slice_100 = slice(f_stop + 1, min(chrom_len, f_stop + 101))
                upstream_slice_50 = slice(f_stop + 1, min(chrom_len, f_stop + 51))

            upstream_100_hits = int(hits_outside_features[upstream_slice_100].sum())
            upstream_50_hits = int(hits_outside_features[upstream_slice_50].sum())

            # longest hit-free interval in the exonic-collapsed vector
            hit_ixes = [0] + list(np.where(local_hits_exonic > 0)[0] + 1) + [feature.coding_length + 1]

            longest_free_intervals = []
            freedom_indices = []
            kornmann_indices = []
            logit_fis = []
            max_skip_tns = 5  # 0..4 inclusive

            for skip_tns in range(max_skip_tns):
                pairs = list(zip(hit_ixes, hit_ixes[1+skip_tns:] or [feature.coding_length + 1]))
                longest_interval = max((r - l) for l, r in pairs) - 1
                longest_free_intervals.append(int(longest_interval))

                fi = float(longest_interval) / feature.coding_length if feature.coding_length else 0.0
                freedom_indices.append(fi)

                if longest_interval >= 300 and 0.1 < fi < 0.9:
                    ki = (longest_interval * hits_in_feature_count) / (feature.coding_length ** 1.5) if feature.coding_length else 0.0
                else:
                    ki = 0.0
                kornmann_indices.append(ki)

                logit_fis.append(fi / (1.0 + math.e ** (-0.01 * (feature.coding_length - 200))) if feature.coding_length else 1.0)

            records[feature.standard_name] = {
                "feature": feature,
                "length": int(len(feature)),
                "hits": hits_in_feature_count,
                "reads": reads_in_feature,
                "neighborhood_hits": hits_outside_feature_count,
                "nc_window_len": int(intergenic_region.coverage),
                "insertion_index": insertion_index,
                "neighborhood_index": neighborhood_index,
                "reads_ni": reads_ni,
                "upstream_hits_50": upstream_50_hits,
                "upstream_hits_100": upstream_100_hits,
                "max_free_region": longest_free_intervals[0],
                "freedom_index": freedom_indices[0],
                "logit_fi": logit_fis[0],
                "kornmann_domain_index": kornmann_indices[4],
                "s_value": log2(reads_in_feature + 1) - total_reads_log,
                "hit_locations": [ix + 1 for (ix, hval) in enumerate(local_hits_exonic) if hval > 0],
                "longest_interval": longest_free_intervals[4],
                "domain_ratio": float(feature.domains.coverage) / feature.coding_length if feature.coding_length else 0.0,
                "hits_in_domains": hits_in_domains_count,
                "reads_in_domains": reads_in_domains,
                "domain_coverage": int(feature.domains.coverage),
                "bps_between_hits_in_neighborhood": (intergenic_region.coverage / hits_outside_feature_count) if hits_outside_feature_count > 0 else 9999,
                "bps_between_hits_in_feature": (len(feature) / hits_in_feature_count) if hits_in_feature_count > 0 else 9999,
            }

        result.update(records)

    return result

# ----------------------------- Output writers -----------------------------

def write_analyzed_records_generic(records, output_file):
    """
    Write a clean per-gene CSV without albicans-only enrichment.
    """
    rows = []
    for r in records:
        f = r["feature"]
        rows.append({
            "Standard name": f.standard_name,
            "Common name": getattr(f, "common_name", f.standard_name),
            "Type": getattr(f, "type", "ORF"),
            "Hits": r.get("hits", 0),
            "Reads": r.get("reads", 0),
            "Length": r.get("length", len(f)),
            "Neighborhood index": r.get("neighborhood_index", 0.0),
            "100 bp upstream hits": r.get("upstream_hits_100", 0),
            "Max free region": r.get("max_free_region", 0),
            "Freedom index": r.get("freedom_index", 0.0),
            "Logit FI": r.get("logit_fi", 0.0),
            "Insertion index": r.get("insertion_index", 0.0),
            "Reads NI": r.get("reads_ni", 0.0),
            "Neighborhood hits": r.get("neighborhood_hits", 0),
            "NC window len": r.get("nc_window_len", 0),
            "S value": r.get("s_value", 0.0),
        })
    df = pd.DataFrame(rows)
    df.to_csv(output_file, index=False)

def write_hits_into_bed(target_file, hits):
    with open(target_file, "w") as bed_file:
        bed_file.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n")
        for hit in sorted(hits, key=lambda o: (o["chrom"], o["hit_pos"])):
            bed_file.write(
                "%s\t%d\t%d\t.\t%d\t%s\n" %
                (hit["chrom"], hit["hit_pos"] - 1, hit["hit_pos"], hit["hit_count"],
                 {"W": "+", "C": "-"}.get(hit["source"], "+"))
            )

def compare_insertion_and_neighborhood(payload, output_file):
    analysis_labels, analyses = payload
    with open(output_file, "w") as out:
        for label, analysis in zip(analysis_labels, analyses):
            s1 = np.array([r["insertion_index"] for r in analysis.values()])
            s2 = np.array([r["neighborhood_index"] for r in analysis.values()])
            if len(s1) < 3:
                out.write(f"{label}\nPearson: NA\nSpearman: NA\n")
                continue
            pearson = scipy.stats.pearsonr(s1, s2)
            spearman = scipy.stats.spearmanr(s1, s2)
            out.write("%s\nPearson: %s\nSpearman: %s\n" % (label, pearson, spearman))

# ----------------------------- Main ---------------------------------------

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input-dir", default=".", help="Input dir containing *_hits.txt")
    parser.add_argument("-o", "--output-dir", default=".", help="Output dir")
    parser.add_argument("-f", "--read-depth-filter", type=int, default=1, help="Ignore hits < this read depth")
    parser.add_argument("-c", "--pairwise-correlations", default=False, action="store_true",
                        help="Perform pairwise correlations (requires multiple hit files).")

    parser.add_argument("--chr-fasta", required=True, help="Reference genome FASTA (chrom lengths)")
    parser.add_argument("--gff3", required=True, help="Annotation GFF3 (CDS-based gene models)")
    parser.add_argument("--ignored-bed", default=None, help="Optional BED of unmappable/ignored regions")

    parser.add_argument("--window", type=int, default=10000, help="Neighborhood window size (default 10000)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # build feature DB
    chrom_lengths = read_fasta_lengths(args.chr_fasta)
    ignored = read_bed_ignored(args.ignored_bed, chrom_lengths) if args.ignored_bed else {c: RangeSet() for c in chrom_lengths}
    feats_by_chrom = parse_gff3_features(args.gff3, chrom_lengths, feature_types=("CDS",))
    feature_db = FeatureDB(chrom_lengths, feats_by_chrom, ignored)

    hit_file_template = "*_hits.txt"
    input_file_paths = sorted(glob.glob(os.path.join(args.input_dir, hit_file_template)))
    input_filenames = [os.path.split(p)[-1].replace("_hits.txt", "") for p in input_file_paths]

    if not input_filenames:
        raise SystemExit(f"ERROR: no hit files matching '{hit_file_template}' found in {args.input_dir}")

    all_hits = read_hit_files(input_file_paths, args.read_depth_filter)
    all_analyzed = [analyze_hits(ds, feature_db, args.window) for ds in all_hits]

    # per-sample analysis table
    for fname, analysis in zip(input_filenames, all_analyzed):
        out_csv = os.path.join(args.output_dir, f"{fname}_analysis.csv")
        write_analyzed_records_generic(list(analysis.values()), out_csv)

    # correlations file
    compare_insertion_and_neighborhood((input_filenames, all_analyzed),
                                       os.path.join(args.output_dir, "insertion_vs_neighborhood_correlations.txt"))

    # stats
    with open(os.path.join(args.output_dir, "stats.csv"), "w", newline="") as stats_file:
        writer = csv.writer(stats_file)
        writer.writerow(["File name"] + ALL_STATS)
        for fname, dataset in zip(input_filenames, all_hits):
            stats = get_statistics(dataset, feature_db)
            writer.writerow([fname] + [stats[col] for col in ALL_STATS])

    # BED per sample (hit sites)
    for fname, dataset in zip(input_filenames, all_hits):
        write_hits_into_bed(os.path.join(args.output_dir, f"{fname}.filter_{args.read_depth_filter}.bed"), dataset)

    # Simple hit_summary across libraries
    if len(all_analyzed) > 1:
        # align on gene names intersection
        all_gene_sets = [set(a.keys()) for a in all_analyzed]
        common_genes = sorted(set.intersection(*all_gene_sets)) if all_gene_sets else []
        with open(os.path.join(args.output_dir, f"hit_summary.RDF_{args.read_depth_filter}.csv"), "w", newline="") as out:
            w = csv.writer(out)
            w.writerow(["Standard name"] + [f"Lib {n}" for n in input_filenames])
            for g in common_genes:
                w.writerow([g] + [a[g]["hits"] for a in all_analyzed])

if __name__ == "__main__":
    main()
