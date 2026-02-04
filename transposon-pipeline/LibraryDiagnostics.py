
#!/usr/bin/env python3
"""
LibraryDiagnostics.py  - Species-agnostic Tn-seq library QC diagnostics

WHAT IT DOES
------------
Given one or more insertion-site count tables (per library / replicate), compute:

1) MidLC complexity curve (Gale et al. 2020 style subsampling)
   - Outputs: <sample>.midlc.csv
   - Also estimates a MidLC scalar (half-saturation proxy) from the curve for reporting/normalization.

2) OPTIONAL: MidLC-based library depth normalization (Gale et al. "depth correction" idea)
   - If --normalize-to-midlc TARGET is provided:
       For each library, compute MidLC_est = depth where mean unique sites reaches ~50% of max
       Compute reads_target = TARGET * MidLC_est
       If total_reads > reads_target: downsample counts by binomial thinning p = reads_target/total_reads
       If total_reads <= reads_target: do nothing (never upsample)
   - Outputs:
       * If input mode is --hits-txt:
           <sample>_normalized_hits.txt   (mirrors original hits format; ready for SummaryTable.py)
       * Otherwise:
           <sample>.<suffix>.sites.bedgraph4
           <sample>.<suffix>.SiteCount.csv

3) Insertion-site sequence bias (+2 and +7 bases relative to insertion coordinate)
   - Outputs: <sample>.seqbias_2_7.tsv (+ background frequencies if --fasta given)

4) Centromere bias (1 kb bins from centromere, averaged across chromosome arms)
   - Outputs: <sample>.centromere_bins.tsv and <sample>.centromere_bias.png

5) OPTIONAL: Per-gene top-site removal (your older "remove top insertion per gene" behavior)
   - If --remove-top-site-per-gene is set AND --gene-gff is provided:
       For each gene, identify the insertion site within that gene with maximum read count.
       Remove (set to 0) that ONE site per gene (ties: remove one arbitrary max site).
   - This occurs BEFORE MidLC QC and BEFORE normalization.
   - Intended to mirror:
        newlistTopRemoved = list(newlist); sort desc; del top1

OUTPUT ORGANIZATION
-------------------
Outputs are written to per-sample subdirectories under --outdir:

  <outdir>/<sample>/<sample>.midlc.csv
  <outdir>/<sample>/<sample>.seqbias_2_7.tsv
  <outdir>/<sample>/<sample>.centromere_bins.tsv
  <outdir>/<sample>/<sample>.centromere_bias.png
  <outdir>/<sample>/<sample>_normalized_hits.txt  (if normalization in hits mode)

Per-sample summary CSV:
  <outdir>/<sample>/<sample>.summary.csv

Combined summary CSV for all samples (one level above sample dirs):
  <outdir>/library_diagnostics.summary.csv

NOTES
-----
- Per-gene top-site removal requires mapping sites to genes, therefore needs a GFF/GTF.
- GFF is assumed 1-based inclusive. We convert to 1-based intervals for membership tests.
- For speed, we keep per-chrom gene interval lists and do a simple scan (OK for yeast).

"""

from __future__ import annotations

import argparse
import os
import sys
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Tuple, List, Optional, Iterable

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ------------------------------ I/O helpers ------------------------------

def die(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)

def mkdir_p(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def infer_sample_name(path: str) -> str:
    base = os.path.basename(path)
    for suf in [
        ".sites.bedgraph4",
        ".sites.bedgraph",
        ".bedgraph4",
        ".bedgraph",
        "SiteCount.csv",
        ".csv",
        "_hits.txt",
        ".hits.txt",
        ".txt",
    ]:
        if base.endswith(suf):
            base = base[: -len(suf)]
            break
    return base


# ------------------------------ data model ------------------------------

SiteCounts = Dict[Tuple[str, int], int]  # (chrom, site_1based) -> total_reads

@dataclass
class LibraryQCResult:
    sample: str
    input_path: str
    total_reads: int
    unique_sites: int
    midlc_est: int
    depth_ratio: float
    per_gene_top_site_removed: int
    per_gene_sites_removed: int
    normalize_target: Optional[float]
    reads_target: Optional[int]
    thinning_p: Optional[float]
    reads_after_norm: Optional[int]
    normalized_written: int
    normalized_output: str
    midlc_points: int
    midlc_max_sampled: int
    seqbias_done: int
    centromere_done: int


# ------------------------------ readers ----------------------------------

def read_sites_bedgraph4(path: str) -> SiteCounts:
    sc: SiteCounts = defaultdict(int)
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom, start0, _end0, count = parts[0], parts[1], parts[2], parts[3]
            try:
                site_1based = int(start0) + 1
                c = int(float(count))
            except Exception:
                continue
            if c > 0:
                sc[(chrom, site_1based)] += c
    return dict(sc)

def read_sitecount_csv(path: str) -> SiteCounts:
    sc: SiteCounts = defaultdict(int)
    with open(path, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            row = [x.strip() for x in row]
            if len(row) < 4:
                continue
            chrom = row[0]
            try:
                site_1based = int(row[1])
                total = int(float(row[3]))
            except Exception:
                continue
            if total > 0:
                sc[(chrom, site_1based)] += total
    return dict(sc)

def read_hits_txt(path: str) -> SiteCounts:
    sc: SiteCounts = defaultdict(int)
    with open(path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        colmap = {name.strip(): i for i, name in enumerate(header)}
        need = ["Chromosome", "Hit position", "Hit count"]
        if not all(k in colmap for k in need):
            die(f"{path}: missing required columns {need}. Found: {header}")

        i_chr = colmap["Chromosome"]
        i_pos = colmap["Hit position"]
        i_cnt = colmap["Hit count"]

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) <= max(i_chr, i_pos, i_cnt):
                continue
            chrom = parts[i_chr]
            try:
                site_1based = int(parts[i_pos])
                total = int(float(parts[i_cnt]))
            except Exception:
                continue
            if total > 0:
                sc[(chrom, site_1based)] += total
    return dict(sc)


# ------------------------------ FASTA ------------------------------------

def read_fasta_as_dict(path: str) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    name = None
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in seqs:
                    die(f"FASTA has duplicate contig name: {name}")
                seqs[name] = []
            else:
                if name is None:
                    die("FASTA parse error: sequence before header")
                seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}

def genome_pair_background(seqs: Dict[str, str], offsets: Tuple[int, int]) -> Counter:
    off1, off2 = offsets
    bg = Counter()
    for _chrom, s in seqs.items():
        n = len(s)
        min_site = 1
        max_site = n
        min_site = max(min_site, 1 - (min(off1, off2) - 1))
        max_site = min(max_site, n - (max(off1, off2) - 1))
        if max_site < min_site:
            continue
        for site in range(min_site, max_site + 1):
            i1 = site + off1 - 1 - 1
            i2 = site + off2 - 1 - 1
            b1 = s[i1]
            b2 = s[i2]
            if b1 in "ACGT" and b2 in "ACGT":
                bg[b1 + b2] += 1
    return bg


# ------------------------------ GFF gene intervals ------------------------

# gene_intervals_by_chrom[chrom] = list of (start1, end1, gene_id)
GeneIntervalsByChrom = Dict[str, List[Tuple[int, int, str]]]

def _parse_attr(attr: str) -> Dict[str, str]:
    # GFF3: key=value;key2=value2
    d: Dict[str, str] = {}
    for chunk in attr.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" in chunk:
            k, v = chunk.split("=", 1)
            d[k.strip()] = v.strip()
    return d

def load_gene_intervals_from_gff(
    gff_path: str,
    feature_types: Iterable[str],
    gene_id_attr: str = "ID",
) -> GeneIntervalsByChrom:
    ft_set = {ft.lower() for ft in feature_types}
    out: GeneIntervalsByChrom = defaultdict(list)

    with open(gff_path, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _source, ftype, start1, end1, _score, _strand, _phase, attr = parts
            if ftype.strip().lower() not in ft_set:
                continue
            try:
                s1 = int(start1)
                e1 = int(end1)
            except Exception:
                continue
            if e1 < s1:
                s1, e1 = e1, s1

            attrs = _parse_attr(attr)
            gene_id = attrs.get(gene_id_attr)
            if gene_id is None and gene_id_attr != "gene_id":
                gene_id = attrs.get("gene_id")
            if gene_id is None:
                gene_id = attrs.get("Name") or attrs.get("locus_tag") or "UNKNOWN"

            out[chrom].append((s1, e1, gene_id))

    for chrom in list(out.keys()):
        out[chrom].sort(key=lambda x: (x[0], x[1], x[2]))
    return dict(out)

def compute_per_gene_max_sites(
    sc: SiteCounts,
    genes_by_chrom: GeneIntervalsByChrom,
) -> Tuple[set, int]:
    """
    Determine ONE max-count site per gene (ties: keep first encountered max).
    Returns:
      drop_sites: set of (chrom, site_1based) to remove (one per gene)
      n_genes_with_site: number of genes for which we removed a site (=len(drop_sites))
    """
    best: Dict[str, Tuple[int, Tuple[str, int]]] = {}  # gene -> (count, (chrom, site))

    for (chrom, site), c in sc.items():
        if c <= 0:
            continue
        intervals = genes_by_chrom.get(chrom)
        if not intervals:
            continue

        # Scan intervals; yeast-scale OK.
        for s1, e1, gid in intervals:
            if s1 <= site <= e1:
                prev = best.get(gid)
                if prev is None or c > prev[0]:
                    best[gid] = (c, (chrom, site))
                # If tie, keep earlier (do nothing).
                break
    drop = set(v[1] for v in best.values())
    return drop, len(drop)

def apply_drop_sites(sc: SiteCounts, drop_sites: set) -> SiteCounts:
    if not drop_sites:
        return dict(sc)
    out: SiteCounts = {}
    for k, c in sc.items():
        if k in drop_sites:
            continue
        if c > 0:
            out[k] = int(c)
    return out


# ------------------------------ MidLC ------------------------------------

def _expand_reads_exact(sc: SiteCounts) -> List[int]:
    sites = list(sc.keys())
    site_index = {k: i for i, k in enumerate(sites)}
    out: List[int] = []
    for k, c in sc.items():
        out.extend([site_index[k]] * int(c))
    return out

def midlc_curve(
    sc: SiteCounts,
    start_reads: int = 100,
    factor: float = 4.0,
    trials: int = 3,
    exact_max_reads: int = 2_000_000,
    seed: int = 1,
) -> Tuple[List[int], List[List[int]]]:
    rng = np.random.default_rng(seed)

    total_reads = int(sum(sc.values()))
    if total_reads <= 0:
        return [0], [[0] for _ in range(trials)]

    sampled_depths: List[int] = []
    uniques_by_trial: List[List[int]] = [[] for _ in range(trials)]

    n = int(start_reads)
    while n < total_reads:
        sampled_depths.append(n)
        n = int(round(n * factor))
        if n == sampled_depths[-1]:
            n += 1
    sampled_depths.append(total_reads)

    if total_reads <= exact_max_reads:
        expanded = _expand_reads_exact(sc)
        expanded = np.array(expanded, dtype=np.int32)
        for depth in sampled_depths:
            for t in range(trials):
                idx = rng.choice(expanded.shape[0], size=depth, replace=False)
                uniques_by_trial[t].append(int(np.unique(expanded[idx]).size))
        return sampled_depths, uniques_by_trial

    sites = list(sc.keys())
    weights = np.array([sc[k] for k in sites], dtype=np.float64)
    p = weights / weights.sum()
    m = len(sites)

    for depth in sampled_depths:
        for t in range(trials):
            draws = rng.choice(m, size=depth, replace=True, p=p)
            uniques_by_trial[t].append(int(np.unique(draws).size))

    return sampled_depths, uniques_by_trial

def write_midlc_csv(out_csv: str, depths: List[int], uniques: List[List[int]]) -> None:
    mkdir_p(os.path.dirname(out_csv) or ".")
    with open(out_csv, "w", newline="") as out:
        w = csv.writer(out)
        header = ["Reads Sampled"] + [f"Unique Sites Trial{i+1}" for i in range(len(uniques))]
        w.writerow(header)
        for i, n in enumerate(depths):
            row = [n] + [uniques[t][i] for t in range(len(uniques))]
            w.writerow(row)

def estimate_midlc_from_curve(depths: List[int], uniques: List[List[int]]) -> int:
    if not depths or not uniques or not uniques[0]:
        return 0
    mean_u = np.array(
        [np.mean([uniques[t][i] for t in range(len(uniques))]) for i in range(len(depths))],
        dtype=float,
    )
    umax = float(np.max(mean_u))
    if umax <= 0:
        return 0
    target = umax / 2.0
    for d, u in zip(depths, mean_u):
        if u >= target:
            return int(d)
    return int(depths[-1])


# ------------------------------ Normalization (MidLC) ---------------------

def downsample_sitecounts_binomial(sc: SiteCounts, p: float, seed: int = 1) -> SiteCounts:
    if p >= 1.0:
        return dict(sc)
    if p <= 0.0:
        return {}

    rng = np.random.default_rng(seed)
    out: SiteCounts = {}
    for k, c in sc.items():
        c = int(c)
        if c <= 0:
            continue
        out_c = int(rng.binomial(n=c, p=p))
        if out_c > 0:
            out[k] = out_c
    return out

def write_sites_bedgraph4(out_path: str, sc: SiteCounts) -> None:
    mkdir_p(os.path.dirname(out_path) or ".")
    items = sorted(sc.items(), key=lambda kv: (kv[0][0], kv[0][1]))
    with open(out_path, "w") as out:
        for (chrom, site_1based), c in items:
            start0 = site_1based - 1
            end0 = site_1based
            out.write(f"{chrom}\t{start0}\t{end0}\t{int(c)}\n")

def write_sitecount_csv_simple(out_path: str, sc: SiteCounts, orient: str = "NA") -> None:
    mkdir_p(os.path.dirname(out_path) or ".")
    items = sorted(sc.items(), key=lambda kv: (kv[0][0], kv[0][1]))
    with open(out_path, "w", newline="") as out:
        w = csv.writer(out)
        for (chrom, site_1based), c in items:
            w.writerow([chrom, site_1based, orient, int(c)])

def downsample_hits_txt_binomial(
    in_hits: str,
    out_hits: str,
    p: float,
    seed: int = 1,
    drop_sites: Optional[set] = None,
) -> Tuple[int, int]:
    """
    Write a normalized copy of a hits table where 'Hit count' is thinned:
      new_count ~ Binomial(old_count, p)

    If drop_sites is provided, any row with (Chromosome, Hit position) in drop_sites
    is forced to Hit count = 0 BEFORE thinning.
    """
    rng = np.random.default_rng(seed)
    drop_sites = drop_sites or set()

    total_before = 0
    total_after = 0

    with open(in_hits, "r") as fin:
        header = fin.readline().rstrip("\n").split("\t")
        colmap = {name.strip(): i for i, name in enumerate(header)}
        need = ["Chromosome", "Hit position", "Hit count"]
        if not all(k in colmap for k in need):
            die(f"{in_hits}: missing required columns {need}; cannot write normalized hits.")
        i_chr = colmap["Chromosome"]
        i_pos = colmap["Hit position"]
        i_cnt = colmap["Hit count"]

        mkdir_p(os.path.dirname(out_hits) or ".")
        with open(out_hits, "w") as fout:
            fout.write("\t".join(header) + "\n")

            for line in fin:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) <= max(i_chr, i_pos, i_cnt):
                    continue

                chrom = parts[i_chr]
                try:
                    site = int(parts[i_pos])
                    c = int(float(parts[i_cnt]))
                except Exception:
                    continue

                if c < 0:
                    c = 0

                # Per-gene top-site removal (force to 0) before thinning
                if (chrom, site) in drop_sites:
                    c = 0

                total_before += c

                if c == 0:
                    new_c = 0
                else:
                    new_c = int(rng.binomial(n=c, p=p))

                total_after += new_c
                parts[i_cnt] = str(new_c)
                fout.write("\t".join(parts) + "\n")

    return total_before, total_after


# ------------------------------ Seq bias (+2,+7) --------------------------

def seqbias_2_7(sc: SiteCounts, seqs: Dict[str, str]) -> Tuple[Counter, Counter]:
    obs = Counter()
    for (chrom, site_1based), reads in sc.items():
        if chrom not in seqs:
            continue
        s = seqs[chrom]
        i2 = site_1based + 2 - 1 - 1
        i7 = site_1based + 7 - 1 - 1
        if i2 < 0 or i7 < 0 or i2 >= len(s) or i7 >= len(s):
            continue
        b2 = s[i2]
        b7 = s[i7]
        if b2 in "ACGT" and b7 in "ACGT":
            obs[b2 + b7] += int(reads)

    bg = genome_pair_background(seqs, offsets=(2, 7))
    return obs, bg

def write_seqbias_tsv(out_tsv: str, obs: Counter, bg: Counter) -> None:
    mkdir_p(os.path.dirname(out_tsv) or ".")
    dinucs = [a + b for a in "ACGT" for b in "ACGT"]
    obs_total = sum(obs.values()) or 1
    bg_total = sum(bg.values()) or 1

    with open(out_tsv, "w") as out:
        out.write("pair\tobs_reads\tobs_frac\tbg_sites\tbg_frac\tenrichment\n")
        for d in dinucs:
            o = obs.get(d, 0)
            b = bg.get(d, 0)
            of = o / obs_total
            bf = b / bg_total
            enr = (of / bf) if bf > 0 else float("inf")
            out.write(f"{d}\t{o}\t{of:.6g}\t{b}\t{bf:.6g}\t{enr:.6g}\n")


# ------------------------------ Centromere bias ---------------------------

def read_centromeres_bed(path: str) -> Dict[str, int]:
    cent = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start0 = int(parts[1]); end0 = int(parts[2])
            except Exception:
                continue
            mid0 = (start0 + end0) // 2
            cent[chrom] = mid0
    return cent

def centromere_bins(
    sc: SiteCounts,
    centromeres_mid0: Dict[str, int],
    bin_size: int = 1000,
    max_bins: int = 500,
) -> Tuple[List[int], List[float], List[int]]:
    per_arm = defaultdict(Counter)

    for (chrom, site_1based), reads in sc.items():
        if chrom not in centromeres_mid0:
            continue
        mid0 = centromeres_mid0[chrom]
        pos0 = site_1based - 1
        dist = abs(pos0 - mid0)
        b = dist // bin_size
        if b > max_bins:
            continue
        side = "L" if pos0 < mid0 else "R"
        per_arm[(chrom, side)][b] += int(reads)

    arms = list(per_arm.values())
    if not arms:
        return [0], [0.0], [0]

    bin_starts = [i * bin_size for i in range(max_bins + 1)]
    means: List[float] = []
    n_arms: List[int] = []
    for b in range(max_bins + 1):
        vals = [c.get(b, 0) for c in arms]
        means.append(float(np.mean(vals)))
        n_arms.append(len(vals))
    return bin_starts, means, n_arms

def write_centromere_tsv(out_tsv: str, bin_starts: List[int], means: List[float], n_arms: List[int]) -> None:
    mkdir_p(os.path.dirname(out_tsv) or ".")
    with open(out_tsv, "w") as out:
        out.write("bin_start_bp\tmean_reads_across_arms\tn_arms\n")
        for b, m, n in zip(bin_starts, means, n_arms):
            out.write(f"{b}\t{m:.6g}\t{n}\n")

def plot_centromere_bias(out_png: str, bin_starts: List[int], means: List[float], sample: str) -> None:
    mkdir_p(os.path.dirname(out_png) or ".")
    xs_kb = [b / 1000.0 for b in bin_starts]
    plt.figure()
    plt.plot(xs_kb, means)
    plt.xlabel("Distance from centromere (kb)")
    plt.ylabel("Mean reads per 1kb bin (avg across arms)")
    plt.title(f"Centromere bias: {sample}")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


# ------------------------------ main driver ------------------------------

def load_sitecounts(mode: str, path: str) -> SiteCounts:
    if mode == "sites-bedgraph4":
        return read_sites_bedgraph4(path)
    if mode == "sitecount-csv":
        return read_sitecount_csv(path)
    if mode == "hits-txt":
        return read_hits_txt(path)
    die(f"Unknown mode: {mode}")

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Species-agnostic Tn-seq library QC: MidLC, sequence bias (+2/+7), centromere bias, and optional MidLC-based normalization."
    )

    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--sites-bedgraph4", nargs="+", help="One or more bedGraph4 files: chrom start end count")
    g.add_argument("--sitecount-csv", nargs="+", help="One or more SiteCount.csv files: chrom, site_1based, orient, total")
    g.add_argument("--hits-txt", nargs="+", help="One or more *_hits.txt files with Chromosome/Hit position/Hit count")

    p.add_argument("--outdir", required=True, help="Output directory (will create per-sample subdirs inside).")
    p.add_argument("--sample-name", default=None,
                   help="Optional single sample name (only valid if exactly one input file).")
    p.add_argument("--seed", type=int, default=1, help="Random seed")

    # MidLC params
    p.add_argument("--midlc-start-reads", type=int, default=100, help="Starting depth for MidLC subsampling")
    p.add_argument("--midlc-factor", type=float, default=4.0, help="Multiplicative factor for MidLC depth steps")
    p.add_argument("--midlc-trials", type=int, default=3, help="Number of random trials per depth")
    p.add_argument("--exact-max-reads", type=int, default=2_000_000,
                   help="Use exact per-read expansion up to this many reads; above uses weighted approximation")

    # MidLC-based normalization (optional)
    p.add_argument("--normalize-to-midlc", type=float, default=None,
                   help="Downsample libraries to TARGET × MidLC_est (MidLC units). Downsample only; never upsample.")
    p.add_argument("--normalized-suffix", default="normalized",
                   help="Suffix used for normalized site outputs (non-hits modes). Default: normalized")
    p.add_argument("--min-sites-for-normalization", type=int, default=1000,
                   help="Skip normalization if unique sites < this threshold (default: 1000).")
    p.add_argument("--normalized-hits-suffix", default="_normalized_hits.txt",
                   help="Suffix for normalized hits output in hits mode (default: _normalized_hits.txt).")

    # per-gene top-site removal (before QC/normalization)
    p.add_argument("--remove-top-site-per-gene", action="store_true",
                   help="Remove ONE highest-count insertion site per gene before QC/normalization (requires --gene-gff).")
    p.add_argument("--gene-gff", default=None,
                   help="GFF/GTF with gene (or CDS) intervals for per-gene top-site removal.")
    p.add_argument("--gene-feature-types", default="gene",
                   help="Comma-separated feature types to treat as gene intervals (default: gene). Examples: gene,CDS")
    p.add_argument("--gene-id-attr", default="ID",
                   help="GFF attribute key for gene identifier (default: ID; fallback tries gene_id/Name/locus_tag).")

    # sequence bias
    p.add_argument("--fasta", default=None, help="Reference FASTA for insertion-site sequence bias (+2/+7)")

    # centromere bias
    p.add_argument("--centromeres-bed", default=None, help="Centromere BED (chrom start0 end0 ...)")
    p.add_argument("--centro-bin-size", type=int, default=1000, help="Bin size (bp) for centromere bias")
    p.add_argument("--centro-max-bins", type=int, default=500, help="Max number of bins (limits max distance)")

    return p.parse_args()

def write_combined_summary_csv(path: str, results: List[LibraryQCResult]) -> None:
    mkdir_p(os.path.dirname(path) or ".")
    header = [
        "sample",
        "input_path",
        "total_reads",
        "unique_sites",
        "midlc_est",
        "depth_ratio_R_over_midlc",
        "per_gene_top_site_removed",
        "per_gene_sites_removed",
        "normalize_target",
        "reads_target",
        "thinning_p",
        "reads_after_norm",
        "normalized_written",
        "normalized_output",
        "midlc_points",
        "midlc_max_sampled",
        "seqbias_done",
        "centromere_done",
    ]
    with open(path, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(header)
        for r in results:
            w.writerow([
                r.sample,
                r.input_path,
                r.total_reads,
                r.unique_sites,
                r.midlc_est,
                f"{r.depth_ratio:.6g}",
                r.per_gene_top_site_removed,
                r.per_gene_sites_removed,
                "NA" if r.normalize_target is None else r.normalize_target,
                "NA" if r.reads_target is None else r.reads_target,
                "NA" if r.thinning_p is None else f"{r.thinning_p:.6g}",
                "NA" if r.reads_after_norm is None else r.reads_after_norm,
                r.normalized_written,
                r.normalized_output,
                r.midlc_points,
                r.midlc_max_sampled,
                r.seqbias_done,
                r.centromere_done,
            ])

def write_per_sample_summary_csv(path: str, r: LibraryQCResult) -> None:
    # Single-row CSV (same columns as combined)
    write_combined_summary_csv(path, [r])

def main() -> None:
    args = parse_args()
    mkdir_p(args.outdir)

    if args.sites_bedgraph4:
        inputs = args.sites_bedgraph4
        mode = "sites-bedgraph4"
    elif args.sitecount_csv:
        inputs = args.sitecount_csv
        mode = "sitecount-csv"
    else:
        inputs = args.hits_txt
        mode = "hits-txt"

    if args.sample_name and len(inputs) != 1:
        die("--sample-name can only be used with exactly one input file.")

    seqs = read_fasta_as_dict(args.fasta) if args.fasta else None
    cent = read_centromeres_bed(args.centromeres_bed) if args.centromeres_bed else None

    genes_by_chrom: Optional[GeneIntervalsByChrom] = None
    if args.remove_top_site_per_gene:
        if not args.gene_gff:
            die("--remove-top-site-per-gene requires --gene-gff")
        fts = [x.strip() for x in args.gene_feature_types.split(",") if x.strip()]
        genes_by_chrom = load_gene_intervals_from_gff(
            args.gene_gff,
            feature_types=fts,
            gene_id_attr=args.gene_id_attr,
        )
        if not genes_by_chrom:
            die(f"No gene intervals loaded from {args.gene_gff} with feature types: {fts}")

    results: List[LibraryQCResult] = []

    for path in inputs:
        sample = args.sample_name or infer_sample_name(path)
        sample_outdir = os.path.join(args.outdir, sample)
        mkdir_p(sample_outdir)

        sc_raw = load_sitecounts(mode, path)

        # ---- per-gene top-site removal (BEFORE QC/normalization) ----
        drop_sites = set()
        per_gene_top_site_removed = 0
        per_gene_sites_removed = 0
        sc = sc_raw

        if args.remove_top_site_per_gene and genes_by_chrom is not None:
            drop_sites, per_gene_sites_removed = compute_per_gene_max_sites(sc_raw, genes_by_chrom)
            sc = apply_drop_sites(sc_raw, drop_sites)
            per_gene_top_site_removed = 1

        total_reads = int(sum(sc.values()))
        uniq_sites = int(len(sc))

        # ---- MidLC curve ----
        depths, uniques = midlc_curve(
            sc,
            start_reads=args.midlc_start_reads,
            factor=args.midlc_factor,
            trials=args.midlc_trials,
            exact_max_reads=args.exact_max_reads,
            seed=args.seed,
        )
        midlc_csv = os.path.join(sample_outdir, f"{sample}.midlc.csv")
        write_midlc_csv(midlc_csv, depths, uniques)

        # ---- MidLC scalar estimate ----
        midlc_est = estimate_midlc_from_curve(depths, uniques)
        depth_ratio = (total_reads / float(midlc_est)) if midlc_est else 0.0

        # ---- OPTIONAL MidLC normalization (downsample only) ----
        norm_target = args.normalize_to_midlc
        reads_target: Optional[int] = None
        thinning_p: Optional[float] = None
        reads_after_norm: Optional[int] = None
        normalized_written = 0
        normalized_output = "NA"

        if norm_target is not None:
            if uniq_sites >= args.min_sites_for_normalization and midlc_est > 0 and total_reads > 0:
                reads_target = int(round(float(norm_target) * float(midlc_est)))
                if reads_target < total_reads:
                    thinning_p = reads_target / float(total_reads)

                    if mode == "hits-txt":
                        normalized_output = os.path.join(sample_outdir, f"{sample}{args.normalized_hits_suffix}")
                        _before, _after = downsample_hits_txt_binomial(
                            in_hits=path,
                            out_hits=normalized_output,
                            p=thinning_p,
                            seed=args.seed,
                            drop_sites=drop_sites if args.remove_top_site_per_gene else None,
                        )
                        reads_after_norm = _after
                        normalized_written = 1
                    else:
                        sc_norm = downsample_sitecounts_binomial(sc, thinning_p, seed=args.seed)
                        reads_after_norm = int(sum(sc_norm.values()))

                        suf = args.normalized_suffix
                        out_bg = os.path.join(sample_outdir, f"{sample}.{suf}.sites.bedgraph4")
                        out_csv = os.path.join(sample_outdir, f"{sample}.{suf}.SiteCount.csv")
                        write_sites_bedgraph4(out_bg, sc_norm)
                        write_sitecount_csv_simple(out_csv, sc_norm, orient="NA")

                        normalized_output = out_bg
                        normalized_written = 1

        # ---- Seq bias (+2/+7) ----
        seqbias_done = 0
        if seqs is not None:
            obs, bg = seqbias_2_7(sc, seqs)
            seq_tsv = os.path.join(sample_outdir, f"{sample}.seqbias_2_7.tsv")
            write_seqbias_tsv(seq_tsv, obs, bg)
            seqbias_done = 1

        # ---- Centromere bias ----
        centromere_done = 0
        if cent is not None:
            bin_starts, means, n_arms = centromere_bins(
                sc,
                centromeres_mid0=cent,
                bin_size=args.centro_bin_size,
                max_bins=args.centro_max_bins,
            )
            centro_tsv = os.path.join(sample_outdir, f"{sample}.centromere_bins.tsv")
            write_centromere_tsv(centro_tsv, bin_starts, means, n_arms)

            centro_png = os.path.join(sample_outdir, f"{sample}.centromere_bias.png")
            plot_centromere_bias(centro_png, bin_starts, means, sample)
            centromere_done = 1

        r = LibraryQCResult(
            sample=sample,
            input_path=os.path.abspath(path),
            total_reads=total_reads,
            unique_sites=uniq_sites,
            midlc_est=midlc_est,
            depth_ratio=depth_ratio,
            per_gene_top_site_removed=per_gene_top_site_removed,
            per_gene_sites_removed=per_gene_sites_removed,
            normalize_target=norm_target,
            reads_target=reads_target,
            thinning_p=thinning_p,
            reads_after_norm=reads_after_norm,
            normalized_written=normalized_written,
            normalized_output=normalized_output,
            midlc_points=len(depths),
            midlc_max_sampled=(max(depths) if depths else 0),
            seqbias_done=seqbias_done,
            centromere_done=centromere_done,
        )
        results.append(r)

        # per-sample summary CSV
        per_sample_summary = os.path.join(sample_outdir, f"{sample}.summary.csv")
        write_per_sample_summary_csv(per_sample_summary, r)

    # combined summary CSV (one level above sample dirs = args.outdir)
    combined_summary = os.path.join(args.outdir, "library_diagnostics.summary.csv")
    write_combined_summary_csv(combined_summary, results)

    print(f"Wrote diagnostics for {len(results)} sample(s) to: {args.outdir}")
    print(f"Combined summary: {combined_summary}")

if __name__ == "__main__":
    main()
