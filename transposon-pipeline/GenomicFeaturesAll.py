"""
GenomicFeaturesAll.py

A generic genome-feature wrapper built from GFF3 (optionally FASTA lengths),
adapted to replace organism-specific GenomicFeatures modules.

Notes
-----
- All indices are 1-based (inclusive).
- Chromosome slicing and location queries match the original API.
- Exons are taken from GFF3 child features if present (exon/CDS); otherwise
  default to the full feature span.
"""

import os
import gffutils
from Bio import SeqIO

from RangeSet import RangeSet
from SortedCollection import SortedCollection
from operator import attrgetter


# ----------------------------- Core objects ------------------------------

class _Feature(object):
    """
    Minimal feature class compatible with SummaryTable-style analyses.

    Attributes
    ----------
    strand : {'W','C'}
    start, stop : int (1-based inclusive)
    exons : RangeSet of (start, stop) genomic coords
    domains : RangeSet
    coding_length : int
    """
    def __init__(self, standard_name, common_name, ftype, chrom, start, stop, strand, description="",
                 aliases=None):
        self.primary_name = self.standard_name = standard_name
        self.common_name = common_name or ""
        self.name = self.common_name or self.standard_name
        self.all_names = set([self.standard_name, self.name])

        if aliases:
            self.all_names |= set(a for a in aliases if a)

        self.type = ftype or "ORF"
        self.chromosome = chrom

        # normalize strand into W/C like original pipeline
        if strand in ("+", "W"):
            self.strand = "W"
        elif strand in ("-", "C"):
            self.strand = "C"
        else:
            self.strand = strand  # keep as-is if weird

        self.description = description or ""

        self.start = int(min(start, stop))
        self.stop = int(max(start, stop))
        self._len = self.stop - self.start + 1

        # heuristic: treat protein-coding as ORF-like
        self.is_orf = ("ORF" in self.type) or (self.type in ("gene", "mRNA", "CDS"))

        # by default, one exon covering whole span
        self._set_exons([(self.start, self.stop)])

        # default: no domains
        self.domains = RangeSet()

    def __len__(self):
        return self._len

    def _set_exons(self, exons=None):
        if not exons:
            exons = [(self.start, self.stop)]
        self.exons = RangeSet(exons)
        self.coding_length = self.exons.coverage

    def _set_domains(self, domains=None):
        if domains is None:
            domains = RangeSet()
        self.domains = domains & self.exons

    def aa_to_genomic_nt(self, aa_ix):
        """
        Convert 1-based amino-acid index to genomic nucleotide coordinate.

        Preserves original behavior; assumes CDS-like structure (exons/introns).
        """
        apparent_nt = (aa_ix - 1) * 3 + 1

        exon_lengths = [e - b + 1 for b, e in self.exons]
        intron_lengths = [e - b + 1 for b, e in self.exons.complement(self.start, self.stop)]
        if self.strand == "C":
            exon_lengths.reverse()
            intron_lengths.reverse()

        offset = 0
        apparent_length = 0
        for exon_length, intron_length in zip(exon_lengths, intron_lengths):
            apparent_length += exon_length
            if apparent_nt > apparent_length:
                offset += intron_length
            else:
                break

        if self.strand == "C":
            return self.stop - apparent_nt - offset + 1
        else:
            return self.start + apparent_nt + offset - 1


class Chromosome(object):
    """Chromosome wrapper to preserve original API behavior."""
    def __init__(self, name, length, db):
        self._name = self.name = name
        self._db = db
        self._len = db.get_last_effective_chrom_index(name) if length is None else length

    def __len__(self):
        return self._len

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self._db.get_features_at_range(self._name, (key.start, key.stop))
        return self._db.get_features_at_location(self._name, key)

    def __iter__(self):
        return iter(self.get_features())

    def get_features(self):
        # Keep stable ordering
        return sorted((f for f in self._db.get_all_features() if f.chromosome == self._name),
                      key=lambda f: f.start)


class _FeatureDB(object):
    """
    Base feature DB. Keeps original methods used throughout the pipeline.
    """
    def __init__(self, features, fasta_filename=None):
        self._features = tuple(features)
        self._name_map = {}
        self._chrom_features = {}

        duplicate_names = set()
        for feature in features:
            for name in feature.all_names:
                if name not in duplicate_names:
                    if name in self._name_map:
                        duplicate_names.add(name)
                        del self._name_map[name]
                    else:
                        self._name_map[name] = feature

            chrom_name = feature.chromosome
            if chrom_name not in self._chrom_features:
                self._chrom_features[chrom_name] = SortedCollection(key=attrgetter("start"))
            self._chrom_features[chrom_name].insert(feature)

        self._max_feature_lengths = {}
        for chrom, feats in self._chrom_features.items():
            self._max_feature_lengths[chrom] = max(len(f) for f in feats) if len(feats) else 1

        self._chrom_names = {}
        for chrom in self._chrom_features:
            self._chrom_names[chrom] = chrom

        self._create_chrom_cache(fasta_filename)

    def get_feature_by_name(self, feature_name):
        return self._name_map.get(feature_name)

    def get_all_features(self):
        return self._features

    def get_std_chrom_name(self, chrom_alias):
        return self._chrom_names[chrom_alias]

    def get_last_effective_chrom_index(self, chrom):
        chrom = self._chrom_names[chrom]
        return max(f.stop for f in self._chrom_features[chrom]) if self._chrom_features.get(chrom) else 0

    def get_features_at_location(self, chrom, location):
        return self.get_features_at_range(chrom, (location, location))

    def get_features_at_range(self, chrom, start_stop):
        start, stop = start_stop
        chrom = self._chrom_names[chrom]
        feats = self._chrom_features[chrom]

        left_border = start - self._max_feature_lengths[chrom]
        try:
            first_ix = feats.index(feats.find_lt(left_border))
        except ValueError:
            first_ix = 0

        out = []
        i = first_ix
        while i < len(feats):
            f = feats[i]
            if f.start > stop:
                break
            if f.start <= stop and f.stop >= start:
                out.append(f)
            i += 1
        return out

    def get_interfeature_range(self, chrom_name, start_stop):
        start, stop = start_stop
        chrom = self[chrom_name]
        start = max(start, 1)
        stop = min(len(chrom), stop)

        feats = chrom[start:stop]
        return RangeSet([(f.start, f.stop) for f in feats]).complement(start, stop)

    def _create_chrom_cache(self, fasta_filename=None):
        self._cached_chroms = {}

        chrom_lens = {}
        if fasta_filename is not None:
            with open(fasta_filename, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    chrom_lens[record.id] = len(record)

        for name in self._chrom_features:
            self._cached_chroms[name] = Chromosome(name, chrom_lens.get(name), self)

    def __getitem__(self, key):
        return self._cached_chroms.get(self._chrom_names.get(key))

    def __iter__(self):
        for chrom in self._cached_chroms.values():
            yield chrom


# ----------------------------- Generic GFF3 DB ----------------------------

class GenericGFFFeatureDB(_FeatureDB):
    """
    Build a feature DB from a GFF3 file, using gffutils.

    Parameters
    ----------
    gff_filename : str
        GFF3 path
    fasta_filename : str or None
        Optional FASTA path to set chromosome lengths
    featuretype : str
        Which feature type to treat as the primary "gene feature".
        Recommended: 'gene' (default). If your GFF lacks gene entries,
        set to 'mRNA' or similar.
    exon_types : tuple[str]
        Child featuretypes to use as exons (in priority order).
        Recommended: ('exon', 'CDS') or ('CDS',).
    name_priority : tuple[str]
        GFF attribute keys used to name features.
        We'll try in order, then fall back to the gffutils feature id.
    """
    def __init__(self,
                 gff_filename,
                 fasta_filename=None,
                 featuretype="gene",
                 exon_types=("exon", "CDS"),
                 name_priority=("gene", "gene_id", "Name", "locus_tag", "ID")):

        db_path = gff_filename + ".gffutils_db.sqlite"
        if not os.path.exists(db_path):
            gff_db = gffutils.create_db(
                gff_filename,
                dbfn=db_path,
                merge_strategy="create_unique",
                keep_order=True,
                force=True
            )
        else:
            gff_db = gffutils.FeatureDB(db_path, keep_order=True)

        features = []

        for gff_feat in gff_db.all_features(featuretype=featuretype):
            # choose name
            std_name = None
            for k in name_priority:
                v = gff_feat.attributes.get(k)
                if v:
                    std_name = v[0]
                    break
            if not std_name:
                std_name = gff_feat.id

            common_name = gff_feat.attributes.get("Name", [std_name])[0]
            desc = gff_feat.attributes.get("description", [""])[0]

            # gather exons/CDS children
            exons = []
            for et in exon_types:
                kids = list(gff_db.children(gff_feat, featuretype=et, order_by="start"))
                if kids:
                    exons = [(k.start, k.end) for k in kids]
                    break

            # if gene has no exon/CDS children, fall back to the feature span
            start, stop = gff_feat.start, gff_feat.end
            if exons:
                start = min(s for s, e in exons)
                stop = max(e for s, e in exons)

            f = _Feature(
                standard_name=std_name,
                common_name=common_name,
                ftype=("ORF" if featuretype in ("gene", "mRNA") else featuretype),
                chrom=gff_feat.seqid,
                start=start,
                stop=stop,
                strand=gff_feat.strand,
                description=desc,
                aliases=[gff_feat.id]  # keep gff id as alias
            )
            f._set_exons(exons if exons else [(f.start, f.stop)])
            features.append(f)

        super(GenericGFFFeatureDB, self).__init__(features, fasta_filename=fasta_filename)


# ----------------------------- Convenience ctor ---------------------------

def default_all_db(gff_filename, fasta_filename=None,
                   featuretype="gene",
                   exon_types=("exon", "CDS"),
                   name_priority=("gene", "gene_id", "Name", "locus_tag", "ID")):
    """
    Convenience wrapper so downstream scripts can do:
        import GenomicFeaturesAll as GenomicFeatures
        db = GenomicFeatures.default_all_db(gff, fasta)
    """
    return GenericGFFFeatureDB(
        gff_filename=gff_filename,
        fasta_filename=fasta_filename,
        featuretype=featuretype,
        exon_types=exon_types,
        name_priority=name_priority
    )
