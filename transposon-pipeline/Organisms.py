

#!/usr/bin/env python3
"""
Organisms.py : should update to include all yeast species for the community when developing the app


Provides organism-specific resources used by downstream Tn-seq analysis:
- Feature DB
- Homologous/deleted regions (for masking / ignored-feature computation)
- Literature essentials / non-essentials
- Paralog lists

NOTE (Glabrata):
Following Gale et al.2020-style logic, we do NOT mask homologous/deleted regions by default.
So glabrata.homologous_regions and glabrata.deleted_regions are empty dicts unless you
later decide to add masking files explicitly.

Glabrata ↔ S. cerevisiae orthology:
We support direct mapping using:
  dependencies/glabrata/C_glabrata_BG2_S_cerevisiae_orthologs.txt

Columns we use
  - cglab_id  (glabrata feature id; matches FeatureDB standard_name)
  - scer_id   (S. cerevisiae ORF id; e.g., YNL318C)
"""

import os
import csv
import pandas as pd

import GenomicFeatures
from RangeSet import RangeSet
import Shared


class Organism(object):
    _KEYS = (
        "feature_db",
        "homologous_regions",
        "deleted_regions",
        "ignored_regions",
        "ignored_features",
        "literature_essentials",
        "literature_non_essentials",
        "genes_with_paralogs",
    )

    def __init__(self):
        self._attr_cache = {}

        self._ignore_region_threshold = 50
        self._min_feature_coverage = 0.95

        # Autocomplete placeholders
        self.feature_db = None
        self.homologous_regions = None
        self.deleted_regions = None
        self.ignored_regions = None
        self.ignored_features = None
        self.literature_essentials = None
        self.literature_non_essentials = None
        self.genes_with_paralogs = None

    @Shared.memoized
    def _read_range_data(self, file_name):
        """
        Reads a CSV with columns: chrom,start,stop (header row present).
        Returns: dict chrom -> RangeSet
        """
        result = {}
        with open(file_name, "r", newline="") as in_file:
            reader = csv.reader(in_file)
            next(reader, None)  # header
            for chrom, start, stop in reader:
                chrom = self.feature_db.get_std_chrom_name(chrom)
                result.setdefault(chrom, []).append((int(start), int(stop)))

        return {chrom: RangeSet(ranges) for chrom, ranges in result.items()}

    def __getattribute__(self, attr):
        if attr in Organism._KEYS:
            cache = super(Organism, self).__getattribute__("_attr_cache")
            if attr not in cache:
                getter = super(Organism, self).__getattribute__(f"_get_{attr}")
                cache[attr] = getter()
            return cache[attr]
        return super(Organism, self).__getattribute__(attr)

    def _get_ignored_regions(self):
        deleted_regions = self.deleted_regions or {}
        homologous_regions = self.homologous_regions or {}
        result = {}

        all_chroms = set(list(deleted_regions.keys()) + list(homologous_regions.keys()))
        for chrom in all_chroms:
            result[chrom] = deleted_regions.get(chrom, RangeSet()) + homologous_regions.get(chrom, RangeSet())

        return result

    def _get_ignored_features(self):
        """
        Ignore features if:
        - too much of the feature overlaps ignored regions (deleted/homologous masking), OR
        - not an ORF
        """
        result = set()
        for f in self.feature_db.get_all_features():
            chrom_ignored = self.ignored_regions.get(f.chromosome, RangeSet())
            frac_ignored = (chrom_ignored & RangeSet([(f.start, f.stop)])).coverage / float(len(f))
            if frac_ignored > (1.0 - self._min_feature_coverage) or (not f.is_orf):
                result.add(f.standard_name)
        return result

    def _get_genes_with_paralogs(self, paralog_filename):
        with open(paralog_filename, "r") as in_file:
            features = (self.feature_db.get_feature_by_name(line.strip()) for line in in_file)
            return {f.standard_name for f in features if f is not None}


class Calbicans(Organism):
    def _get_feature_db(self):
        return GenomicFeatures.default_alb_db()

    def _get_homologous_regions(self):
        ranges = self._read_range_data(Shared.get_dependency(os.path.join("albicans", "homologous_regions.csv")))
        return {
            chrom: RangeSet(r for r in range_set if (r[1] - r[0] + 1) >= self._ignore_region_threshold)
            for chrom, range_set in ranges.items()
        }

    def _get_deleted_regions(self):
        ranges = self._read_range_data(Shared.get_dependency(os.path.join("albicans", "deleted_regions.csv")))
        return {
            chrom: RangeSet(r for r in range_set if (r[1] - r[0] + 1) >= self._ignore_region_threshold)
            for chrom, range_set in ranges.items()
        }

    def _get_literature_essentials(self):
        return None

    def _get_literature_non_essentials(self):
        return None

    def _get_genes_with_paralogs(self):
        return Organism._get_genes_with_paralogs(
            self, Shared.get_dependency(os.path.join("albicans", "hasParalogs_ca.txt"))
        )


class Scerevisiae(Organism):
    def _get_feature_db(self):
        return GenomicFeatures.default_cer_db()

    def _get_homologous_regions(self):
        ranges = self._read_range_data(Shared.get_dependency(os.path.join("cerevisiae", "homologous_regions.csv")))
        return {
            chrom: RangeSet(r for r in range_set if (r[1] - r[0] + 1) >= self._ignore_region_threshold)
            for chrom, range_set in ranges.items()
        }

    def _get_deleted_regions(self):
        # data unavailable in original pipeline
        return {}

    def _get_literature_essentials(self):
        annotated_as_inviable, annotated_as_viable = self.get_nominal_annotations()
        return annotated_as_inviable - annotated_as_viable

    def _get_literature_non_essentials(self):
        annotated_as_inviable, annotated_as_viable = self.get_nominal_annotations()
        return annotated_as_viable - annotated_as_inviable
    
    @Shared.memoized
    def get_nominal_annotations(self):
        """
        Returns (essential_set, nonessential_set) as standard names.

        NOTE: Keeps original logic: uses SGD viability tables where "Mutant Information" == "null".
        """
        viable_filepath = Shared.get_dependency("cerevisiae", "cerevisiae_viable_annotations.txt")
        inviable_filepath = Shared.get_dependency("cerevisiae", "cerevisiae_inviable_annotations.txt")

        viable_table = pd.read_csv(viable_filepath, skiprows=8, delimiter="\t", na_filter=False)
        inviable_table = pd.read_csv(inviable_filepath, skiprows=8, delimiter="\t", na_filter=False)

        viable_genes = set(viable_table[viable_table["Mutant Information"] == "null"]["Gene"])
        inviable_genes = set(inviable_table[inviable_table["Mutant Information"] == "null"]["Gene"])

        annotated_as_viable = {self.feature_db.get_feature_by_name(g) for g in viable_genes} - {None}
        annotated_as_inviable = {self.feature_db.get_feature_by_name(g) for g in inviable_genes} - {None}

        consensus_viable_orfs = [f for f in annotated_as_viable if f.is_orf and f.feature_qualifier != "Dubious"]
        consensus_inviable_orfs = [f for f in annotated_as_inviable if f.is_orf and f.feature_qualifier != "Dubious"]

        # original returns (inviable, viable)
        return (
            {f.standard_name for f in consensus_inviable_orfs},
            {f.standard_name for f in consensus_viable_orfs},
        )

    @property
    @Shared.memoized
    def conflicting_essentials(self):
        return (self.literature_essentials or set()) & (self.literature_non_essentials or set())

    def _get_genes_with_paralogs(self):
        return Organism._get_genes_with_paralogs(
            self, Shared.get_dependency(os.path.join("cerevisiae", "hasParalogs_sc.txt"))
        )


class Spombe(Organism):
    def _get_feature_db(self):
        return GenomicFeatures.default_pom_db()

    def _get_homologous_regions(self):
        ranges = self._read_range_data(Shared.get_dependency(os.path.join("pombe", "homologous_regions.csv")))
        return {
            chrom: RangeSet(r for r in range_set if (r[1] - r[0] + 1) >= self._ignore_region_threshold)
            for chrom, range_set in ranges.items()
        }

    def _get_deleted_regions(self):
        return {}

    def _get_literature_essentials(self):
        return self._get_spom_essentials()[0]

    def _get_literature_non_essentials(self):
        return self._get_spom_essentials()[1]

    @Shared.memoized
    def _get_spom_essentials(self):
        viability_table = pd.read_csv(
            Shared.get_dependency("pombe/FYPOviability.tsv"),
            header=None,
            delimiter="\t",
            names=["pombe standard name", "essentiality"],
        )
        ess = {r[0] for _ix, r in viability_table.iterrows() if r[1] == "inviable"}
        non = {r[0] for _ix, r in viability_table.iterrows() if r[1] == "viable"}
        return ess, non

    def _get_genes_with_paralogs(self):
        return Organism._get_genes_with_paralogs(
            self, Shared.get_dependency(os.path.join("pombe", "hasParalogs_sp.txt"))
        )


class Cglabrata(Organism):
    """
    C. glabrata (BG2-style FeatureDB)

    By default, no masking files are applied (deleted/homologous regions empty),
    matching Gale et al.-style assumptions (no region masking).
    """

    def _get_feature_db(self):
        return GenomicFeatures.default_glab_db()

    def _get_homologous_regions(self):
        # default: no masking
        return {}

    def _get_deleted_regions(self):
        # default: no masking
        return {}

    def _get_literature_essentials(self):
        # Not wired yet. Add later if you compile Gale et al. essentials list.
        return None

    def _get_literature_non_essentials(self):
        return None

    def _get_genes_with_paralogs(self):
        # Optional: if you create a dependency file, wire it here.
        # If it doesn't exist, return empty set (safe default).
        try:
            path = Shared.get_dependency(os.path.join("glabrata", "hasParalogs_cgla.txt"))
        except Exception:
            print(
                "[WARN] glabrata paralog file not found in dependencies/glabrata/hasParalogs_cgla.txt; using empty set."
            )
            return set()
        return Organism._get_genes_with_paralogs(self, path)


# Lazy singletons:
alb = Calbicans()
cer = Scerevisiae()
pom = Spombe()
gla = Cglabrata()


@Shared.memoized
def get_genes_by_name(name, include_glabrata=False):
    """
    Backward compatible default: returns (alb, cer, pom).
    If include_glabrata=True, returns (alb, cer, pom, gla).
    """
    alb_f = alb.feature_db.get_feature_by_name(name)
    cer_f = cer.feature_db.get_feature_by_name(name)
    pom_f = pom.feature_db.get_feature_by_name(name)
    if include_glabrata:
        gla_f = gla.feature_db.get_feature_by_name(name)
        return (alb_f, cer_f, pom_f, gla_f)
    return (alb_f, cer_f, pom_f)


# ---- Existing Calb/Scer/Spom orthology helpers kept as-is (Python 3 syntax) ----

@Shared.memoized
def get_orths_by_name(name):
    alb_feature, cer_feature, pom_feature = get_genes_by_name(name)

    # If the name was not a Calb name, attempt to recover it:
    if alb_feature is None:
        alb_name = None
        if cer_feature:
            alb_name = cer_to_alb_map().get(cer_feature.standard_name)
        if (not alb_name) and pom_feature:
            alb_name = pom_to_alb_map().get(pom_feature.standard_name)

        if not alb_name:
            return (alb_feature, cer_feature, pom_feature)
        alb_feature = alb.feature_db.get_feature_by_name(alb_name)

    cer_ortholog = None
    if getattr(alb_feature, "cerevisiae_orthologs", None):
        cer_ortholog = cer.feature_db.get_feature_by_name(list(alb_feature.cerevisiae_orthologs)[0])

    pom_ortholog = pom.feature_db.get_feature_by_name(get_calb_orths_in_sp().get(alb_feature.standard_name, ""))

    return (alb_feature, cer_ortholog, pom_ortholog)


@Shared.memoized
def cer_to_alb_map():
    result = {}
    for f in alb.feature_db.get_all_features():
        if not getattr(f, "cerevisiae_orthologs", None):
            continue

        cer_orth = list(f.cerevisiae_orthologs)[0]
        cer_feature = cer.feature_db.get_feature_by_name(cer_orth)
        if cer_feature is None:
            print(f"WARNING: albicans ortholog {f.standard_name} doesn't exist in cerevisiae database!")
            continue

        result[cer_feature.standard_name] = f.standard_name

    return result


@Shared.memoized
def pom_to_alb_map():
    return {
        pom.feature_db.get_feature_by_name(pom_name).standard_name: alb_name
        for alb_name, pom_name in get_calb_orths_in_sp().items()
        if pom.feature_db.get_feature_by_name(pom_name) is not None
    }


@Shared.memoized
def get_calb_orths_in_sp():
    pom_db = GenomicFeatures.default_pom_db()

    # NOTE: leaving file parsing logic intact; assumes dependency files match prior expectations.
    ortholog_table = pd.read_csv(
        Shared.get_dependency("albicans", "C_albicans_SC5314_S_pombe_orthologs.txt"),
        skiprows=8,
        delimiter="\t",
        header=None,
    )

    best_hit_table = pd.read_csv(
        Shared.get_dependency("albicans", "C_albicans_SC5314_S_pombe_best_hits.txt"),
        skiprows=8,
        delimiter="\t",
        header=None,
    )

    joined_table = pd.concat([ortholog_table, best_hit_table], ignore_index=True)

    # Expect columns: [alb_std, alb_common, alb_id, pom_std, pom_common, pom_id] (legacy)
    # Build map alb_std -> pom_feature.name (as in original)
    result = {}
    for alb_feature in GenomicFeatures.default_alb_db().get_all_features():
        # alb std assumed in col 0, pom std assumed in col 3
        ortholog_rows = joined_table[joined_table.iloc[:, 0] == alb_feature.standard_name]
        if ortholog_rows.empty:
            continue

        pom_std = ortholog_rows.iloc[0, 3]
        pom_feature = pom_db.get_feature_by_name(pom_std)
        if pom_feature:
            result[alb_feature.standard_name] = pom_feature.name

    return result


@Shared.memoized
def get_all_orthologs():
    result = []
    for f in alb.feature_db.get_all_features():
        orths = get_orths_by_name(f.standard_name)
        if None not in orths:
            result.append(orths)
    return result


# ---- Glabrata ↔ S. cerevisiae orthology helpers ----

@Shared.memoized
def _load_gla_scer_ortholog_table():
    """
    Load glabrata↔scer ortholog table from dependencies.

    File: dependencies/glabrata/C_glabrata_BG2_S_cerevisiae_orthologs.txt

    We use:
      - cglab_id (glabrata feature id; matches glabrata FeatureDB standard_name)
      - scer_id  (S. cerevisiae ORF id; matches cerevisiae FeatureDB standard_name)

    Notes:
      - This mapping is often many(glabrata) -> one(scer).
      - Reverse mapping scer -> set(glabrata) is therefore 1-to-many.
    """
    fp = Shared.get_dependency("glabrata", "C_glabrata_BG2_S_cerevisiae_orthologs.txt")
    df = pd.read_csv(fp, sep="\t", dtype=str).fillna("")

    required = {"cglab_id", "scer_id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Ortholog table missing columns {sorted(missing)}: {fp}")

    df = df[["cglab_id", "scer_id"]].copy()
    df.columns = ["gla", "scer"]

    df["gla"] = df["gla"].astype(str).str.strip()
    df["scer"] = df["scer"].astype(str).str.strip()
    df = df[(df["gla"] != "") & (df["scer"] != "")]

    return df


@Shared.memoized
def gla_to_scer_map():
    """
    Map: glabrata_standard_name -> scer_standard_name

    If duplicates exist for the same glabrata id, we keep the first occurrence.
    (Your file is primarily many(glabrata)->one(scer), so this is usually safe.)
    """
    df = _load_gla_scer_ortholog_table()
    m = {}
    for gla_name, scer_name in df.itertuples(index=False):
        if gla_name not in m:
            m[gla_name] = scer_name
    return m


@Shared.memoized
def scer_to_gla_map():
    """
    Map: scer_standard_name -> set(glabrata_standard_name)

    Safe for 1-to-many.
    """
    df = _load_gla_scer_ortholog_table()
    m = {}
    for gla_name, scer_name in df.itertuples(index=False):
        m.setdefault(scer_name, set()).add(gla_name)
    return m


@Shared.memoized
def get_gla_scer_orths_by_name(name):
    """
    Given a gene name (either glabrata or scer), return (gla_feature, scer_feature).

    - If 'name' matches glabrata FeatureDB: return (gla_feature, mapped scer_feature).
    - Else if 'name' matches scer FeatureDB: return (one mapped gla_feature, scer_feature).
      (If scer has multiple glabrata orthologs, we return an arbitrary one. If you want
       ALL, use scer_to_gla_map() directly.)
    """
    # Try interpret as glabrata first
    gla_feature = gla.feature_db.get_feature_by_name(name)
    scer_feature = cer.feature_db.get_feature_by_name(name)

    if gla_feature is not None:
        scer_name = gla_to_scer_map().get(gla_feature.standard_name)
        scer_feature = cer.feature_db.get_feature_by_name(scer_name) if scer_name else None
        return (gla_feature, scer_feature)

    if scer_feature is not None:
        gla_names = scer_to_gla_map().get(scer_feature.standard_name, set())
        gla_name = next(iter(gla_names), None)
        gla_feature = gla.feature_db.get_feature_by_name(gla_name) if gla_name else None
        return (gla_feature, scer_feature)

    return (None, None)


def write_ignored_genes_table(organism, out_file):
    columns = [
        "Standard name",
        "Common name",
        "Type",
        "Deleted fraction",
        "Duplicated fraction",
        "Reason for exclusion",
    ]
    ignored_features = [organism.feature_db.get_feature_by_name(n) for n in sorted(organism.ignored_features)]
    data = [
        [
            f.standard_name,
            getattr(f, "common_name", ""),
            getattr(f, "type", ""),
            (f.exons & organism.deleted_regions.get(f.chromosome, RangeSet())).coverage / float(len(f)),
            (f.exons & organism.homologous_regions.get(f.chromosome, RangeSet())).coverage / float(len(f)),
        ]
        for f in ignored_features
        if f is not None
    ]

    for row in data:
        ignore_reasons = []
        if "ORF" not in str(row[2]).upper():
            ignore_reasons.append("Not ORF")
        if row[-2] > 0.05:
            ignore_reasons.append("More than 5% deleted")
        if row[-1] > 0.05:
            ignore_reasons.append("More than 5% duplicated")
        row.append("; ".join(ignore_reasons))

    with pd.ExcelWriter(out_file) as excel_writer:
        df = pd.DataFrame(columns=columns, data=data)
        df.to_excel(excel_writer, sheet_name="Ignored genes", index=False)


if __name__ == "__main__":
    write_ignored_genes_table(alb, os.path.join(Shared.get_script_dir(), "output", "ignored_alb_genes.xlsx"))

