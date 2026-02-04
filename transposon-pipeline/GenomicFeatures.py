
"""
A module for wrapping genome feature tables.

Can work with multiple gene aliases, find features by location (and range),
compute the inter-feature region at a given range, and generally allow easier
access to the features.

Attributes
----------
cer_db : CerevisiaeFeatureDB
    The default cerevisiae feature database.

alb_db : AlbicansFeatureDB
    The default albicans feature database.

pom_db : PombeFeatureDB
    The default pombe feature database.

glab_db : GlabrataFeatureDB
    The default glabrata (BG2) feature database.

Notes
-----
All indices are 1-based.
"""

import os
from operator import attrgetter

import gffutils
import pandas as pd
from Bio import SeqIO

import Shared
from RangeSet import RangeSet
from SortedCollection import SortedCollection


class _FeatureDB(object):
    """The base class for genome feature databases.

    Wraps the tab-delimited tables, such as those used by both CGD and SGD.

    Attributes
    ----------
    _features : sequence of _Feature objects
        A sequence of all of the genome features.
    _name_map : dict of str to _Feature
        A dictionary mapping all known feature names (including aliases) to the
        feature instance itself. Names that map to more than one feature are
        ignored.
    _chrom_features : dict of str to list of _Feature objects
        A map of chromosome names to the features that reside on them. The
        features are sorted by their .start attribute.
    _max_feature_lengths : dict of str to int
        A map of the longest features in each chromosome. This gives a maximum
        search window when searching for features in `get_features_at_range`.
    _chrom_names : dict of obj to str
        A map of chromosome name aliases to their "canonical" names.
    """

    def __init__(self, features, fasta_filename=None):
        self._features = tuple(features)  # Read-only copy
        self._name_map = {}
        self._chrom_features = {}

        duplicate_names = set()
        for feature in features:
            for name in feature.all_names:
                # We will only keep aliases that are unique to the feature:
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

        # Longest feature lengths are required for getting features in ranges:
        self._max_feature_lengths = {}
        for chrom, features_on_chrom in self._chrom_features.items():
            self._max_feature_lengths[chrom] = max(len(f) for f in features_on_chrom)

        # Create basis for chromosome name maps (subclasses can add their own aliases):
        self._chrom_names = {}
        for chrom in self._chrom_features:
            self._chrom_names[chrom] = chrom

        self._create_chrom_cache(fasta_filename)

    def get_features_at_location(self, chrom, location):
        """Get all the features at a particular location in the genome."""
        return self.get_features_at_range(chrom, (location, location))

    def get_features_at_range(self, chrom, start_stop):
        """
        Get features that intersect a given range.

        Parameters
        ----------
        chrom : str or int
            The target chromosome.
        start_stop : tuple(int, int)
            (start, stop), inclusive, 1-based.
        """
        start, stop = start_stop

        chrom = self._chrom_names[chrom]
        features = self._chrom_features[chrom]

        left_border = start - self._max_feature_lengths[chrom]

        try:
            first_feature_ix = features.index(features.find_lt(left_border))
        except ValueError:
            first_feature_ix = 0

        result = []

        # Walk right:
        walk_ix = first_feature_ix
        while walk_ix < len(features):
            f = features[walk_ix]
            if f.start > stop:
                break
            if f.start <= stop and f.stop >= start:
                result.append(f)
            walk_ix += 1

        return result

    def get_last_effective_chrom_index(self, chrom):
        """Return the last chromosome position that is accessible."""
        chrom = self._chrom_names[chrom]
        return max(f.stop for f in self._chrom_features[chrom])

    def get_feature_by_name(self, feature_name):
        return self._name_map.get(feature_name)

    def get_all_features(self):
        return self._features

    def get_interfeature_range(self, chrom_name, start_stop):
        """Return maximal ranges within a given region that contain no features."""
        start, stop = start_stop

        chrom = self[chrom_name]
        start = max(start, 1)
        stop = min(len(chrom), stop)
        features = chrom[start:stop]

        return RangeSet([(f.start, f.stop) for f in features]).complement(start, stop)

    def get_std_chrom_name(self, chrom_alias):
        return self._chrom_names[chrom_alias]

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
        # Use caching because creating the chromosome may trigger expensive computation.
        return self._cached_chroms.get(self._chrom_names.get(key))

    def __iter__(self):
        yield from self._cached_chroms.values()


class AlbicansFeatureDB(_FeatureDB):
    A21, A22 = (21, 22)

    def __init__(self, feature_filename, fasta_filename=None, gff_filename=None, domain_filename=None):
        if gff_filename is not None:
            gff_db_filename = gff_filename + ".gffutils_db.sqlite"
            if not os.path.exists(gff_db_filename):
                gff_db = gffutils.create_db(
                    gff_filename,
                    gff_db_filename,
                    merge_strategy="create_unique",
                    keep_order=True,
                )
            else:
                gff_db = gffutils.FeatureDB(gff_db_filename, keep_order=True)
        else:
            gff_db = None

        # Collect all protein domains. Later, when we create the features, match every feature to its domains.
        domains = {}
        if domain_filename is not None:
            domain_table = pd.read_table(domain_filename, header=None)
            for _, row in domain_table.iterrows():
                e_value = row[9 - 1]
                if e_value != "NA" and float(e_value) >= 0.0001:
                    continue

                feature_name = row[1 - 1]
                domain_start_aa = int(row[7 - 1])
                domain_stop_aa = int(row[8 - 1])

                domains.setdefault(feature_name, []).append((domain_start_aa, domain_stop_aa))

        features = []
        with open(feature_filename, "r") as in_file:
            for line in in_file:
                if not line:
                    continue
                if line.startswith("!"):
                    if "Genome version" in line:
                        if "A21" in line:
                            self.assembly = AlbicansFeatureDB.A21
                        elif "A22" in line:
                            self.assembly = AlbicansFeatureDB.A22
                    continue

                feature = AlbicansFeature(line)
                if (
                    "chrM" in feature.chromosome  # Ignore mitochondria
                    or feature.primary_name.endswith("_B")  # If A22, only keep the A haplotype
                    or feature.start == -1  # Unmapped features aren't useful
                ):
                    continue

                if gff_db is not None:
                    exons = []
                    try:
                        gff_feature = gff_db[feature.standard_name]
                        assert feature.start == gff_feature.start and feature.stop == gff_feature.end

                        for exon in gff_db.children(gff_feature, featuretype="exon"):
                            exons.append((exon.start, exon.end))
                    except gffutils.exceptions.FeatureNotFoundError:
                        pass

                    feature._set_exons(exons)

                feature_domains = RangeSet()
                for feature_name in feature.all_names:
                    if feature_name not in domains:
                        continue

                    for domain_start_aa, domain_stop_aa in domains[feature_name]:
                        start_nt = feature.aa_to_genomic_nt(domain_start_aa)
                        stop_nt = feature.aa_to_genomic_nt(domain_stop_aa)
                        feature_domains = feature_domains + RangeSet([sorted([start_nt, stop_nt])])

                    feature._set_domains(feature_domains)
                    break

                features.append(feature)

        super().__init__(features, fasta_filename)

        # Enrich the chromosome name maps:
        if self.assembly == AlbicansFeatureDB.A21:
            chrom_template = "Ca21chr%s_C_albicans_SC5314"
        else:
            chrom_template = "Ca22chr%sA_C_albicans_SC5314"

        for chrom_shorthand in list(range(1, 8)) + ["R"]:
            chrom_name = chrom_template % chrom_shorthand
            self._chrom_names[chrom_shorthand] = chrom_name
            self._chrom_names[str(chrom_shorthand)] = chrom_name

        # Construct a reverse map of cerevisiae_orthologs:
        self._cer_orths = {}
        for feature in features:
            for cer_orth in feature.cerevisiae_orthologs:
                self._cer_orths[cer_orth] = feature

    def get_feature_by_cerevisiae_ortholog(self, cer_ortholog):
        return self._cer_orths.get(cer_ortholog)


class CerevisiaeFeatureDB(_FeatureDB):
    def __init__(self, feature_filename, fasta_file=None):
        features = []
        introns = {}
        with open(feature_filename, "r") as in_file:
            for line in in_file:
                feature = CerevisiaeFeature(line)

                if feature.type == "intron":
                    parent = feature.parent_feature
                    introns.setdefault(parent, []).append(feature)
                    continue

                if (
                    feature.type in ("telomere", "ARS")  # Ignore unwanted feature types
                    or feature.chromosome == "17"  # Ignore mitochondria
                    or feature.start == -1  # Unmapped features begone!
                    or "chromosome" not in feature.parent_feature  # We don't care about sub-annotations
                ):
                    continue

                features.append(feature)

        # Go over all features and look up their introns:
        for feature in features:
            introns_for_feature = []
            for alias in feature.all_names:
                if alias in introns:
                    introns_for_feature = introns[alias]
                    break

            exons = RangeSet([(intron.start, intron.stop) for intron in introns_for_feature]).complement(
                feature.start, feature.stop
            )
            feature._set_exons(exons)

        super().__init__(features, fasta_file)

        self._chrom_names.update(
            {
                "chrI": "1",
                "chrII": "2",
                "chrIII": "3",
                "chrIV": "4",
                "chrV": "5",
                "chrVI": "6",
                "chrVII": "7",
                "chrVIII": "8",
                "chrIX": "9",
                "chrX": "10",
                "chrXI": "11",
                "chrXII": "12",
                "chrXIII": "13",
                "chrXIV": "14",
                "chrXV": "15",
                "chrXVI": "16",
            }
        )
        self._chrom_names.update(
            {
                "ref|NC_001133|": "1",
                "ref|NC_001134|": "2",
                "ref|NC_001135|": "3",
                "ref|NC_001136|": "4",
                "ref|NC_001137|": "5",
                "ref|NC_001138|": "6",
                "ref|NC_001139|": "7",
                "ref|NC_001140|": "8",
                "ref|NC_001141|": "9",
                "ref|NC_001142|": "10",
                "ref|NC_001143|": "11",
                "ref|NC_001144|": "12",
                "ref|NC_001145|": "13",
                "ref|NC_001146|": "14",
                "ref|NC_001147|": "15",
                "ref|NC_001148|": "16",
            }
        )


class PombeFeatureDB(_FeatureDB):
    def __init__(self, gff_filename, fasta_filename=None):
        gff_db_filename = gff_filename + ".gffutils_db.sqlite"
        if not os.path.exists(gff_db_filename):
            gff_db = gffutils.create_db(
                gff_filename,
                gff_db_filename,
                merge_strategy="create_unique",
                keep_order=True,
            )
        else:
            gff_db = gffutils.FeatureDB(gff_db_filename, keep_order=True)

        features = []
        for gff_feature in gff_db.all_features(featuretype="gene"):
            # Use CDS blocks as exons.
            cds_borders = [(cds.start, cds.end) for cds in gff_db.children(gff_feature, featuretype="CDS")]
            start = min(left for left, _ in cds_borders)
            stop = max(right for _, right in cds_borders)

            components = [
                gff_feature.attributes["gene_id"][0],
                gff_feature.attributes["Name"][0],
                "ORF" if gff_feature.featuretype == "gene" else gff_feature.featuretype,
                gff_feature.seqid,
                start,
                stop,
                {"+": "W", "-": "C"}[gff_feature.strand],
                gff_feature.attributes.get("description", [""])[0],
            ]

            feature = _Feature(components, 0, 1, 2, 3, 4, 5, 6, 7)
            feature._set_exons(cds_borders)
            features.append(feature)

        super().__init__(features, fasta_filename)


class GlabrataFeatureDB(_FeatureDB):
    """
    NCBI-style feature DB for Nakaseomyces glabratus (Candida glabrata) BG2.

    Builds one feature per coding transcript using locus_tag as the stable ID.
    Uses CDS blocks as exons (coding regions), for Tn-seq essentiality metrics.

    Optional: attach InterProScan domains by mapping protein_id -> locus_tag.
    """

    def __init__(
        self,
        gff_filename,
        fasta_filename=None,
        domain_tsv_filename=None,
        interpro_evalue_max=1e-4,
        keep_all_interpro_methods=True,
    ):
        gff_db_filename = gff_filename + ".gffutils_db.sqlite"
        if not os.path.exists(gff_db_filename):
            gff_db = gffutils.create_db(
                gff_filename,
                gff_db_filename,
                merge_strategy="create_unique",
                keep_order=True,
            )
        else:
            gff_db = gffutils.FeatureDB(gff_db_filename, keep_order=True)

        domains_by_protein = {}
        if domain_tsv_filename is not None and os.path.exists(domain_tsv_filename):
            with open(domain_tsv_filename, "r") as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if not line or line.startswith("#"):
                        continue

                    parts = line.split("\t")
                    if len(parts) < 8:
                        continue

                    prot = parts[0].strip()
                    if not prot:
                        continue

                    try:
                        aa_start = int(parts[6])
                        aa_stop = int(parts[7])
                    except ValueError:
                        continue

                    if len(parts) >= 9:
                        ev = parts[8].strip()
                        if ev not in ("", "-", "NA"):
                            try:
                                if float(ev) > float(interpro_evalue_max):
                                    continue
                            except ValueError:
                                pass

                    if keep_all_interpro_methods:
                        pass

                    domains_by_protein.setdefault(prot, []).append((aa_start, aa_stop))

        features = []
        transcript_types = ("mRNA", "transcript")

        for tx in gff_db.all_features(featuretype=transcript_types):
            attrs = tx.attributes

            locus_tag = attrs.get("locus_tag", [None])[0]
            if not locus_tag:
                continue

            gene_name = attrs.get("gene", [None])[0] or attrs.get("Name", [None])[0] or locus_tag
            product = attrs.get("product", [""])[0] if "product" in attrs else ""

            cds_children = list(gff_db.children(tx, featuretype="CDS", order_by="start"))
            if not cds_children:
                continue

            cds_borders = [(cds.start, cds.end) for cds in cds_children]
            start = min(b for b, _ in cds_borders)
            stop = max(e for _, e in cds_borders)

            protein_id = None
            for cds in cds_children:
                protein_id = cds.attributes.get("protein_id", [None])[0]
                if protein_id:
                    break

            strand = {"+": "W", "-": "C"}.get(tx.strand, "W")
            chrom = tx.seqid

            components = [
                locus_tag,
                gene_name,
                "ORF",
                chrom,
                start,
                stop,
                strand,
                product,
            ]
            feature = _Feature(components, 0, 1, 2, 3, 4, 5, 6, 7)

            feature._set_exons(cds_borders)

            feature.all_names.add(locus_tag)
            if gene_name:
                feature.all_names.add(gene_name)

            tx_id = attrs.get("ID", [None])[0]
            if tx_id:
                feature.all_names.add(tx_id)

            name_attr = attrs.get("Name", [None])[0]
            if name_attr:
                feature.all_names.add(name_attr)

            if protein_id:
                feature.protein_id = protein_id
                feature.all_names.add(protein_id)

            if protein_id and protein_id in domains_by_protein:
                feature_domains = RangeSet()
                for aa_start, aa_stop in domains_by_protein[protein_id]:
                    try:
                        start_nt = feature.aa_to_genomic_nt(aa_start)
                        stop_nt = feature.aa_to_genomic_nt(aa_stop)
                        feature_domains = feature_domains + RangeSet([sorted([start_nt, stop_nt])])
                    except Exception:
                        continue
                feature._set_domains(feature_domains)

            features.append(feature)

        super().__init__(features, fasta_filename)

        for chrom in self._chrom_features:
            self._chrom_names[chrom] = chrom


class Chromosome(object):
    """A chromosome in the genome (wrapper around _FeatureDB methods)."""

    def __init__(self, name, length, db):
        self._name = self.name = name
        self._db = db
        self._len = db.get_last_effective_chrom_index(name) if length is None else length

    def __len__(self):
        return self._len

    def __getitem__(self, key):
        if isinstance(key, slice):
            start = 1 if key.start is None else key.start
            stop = len(self) if key.stop is None else key.stop
            return self._db.get_features_at_range(self._name, (start, stop))
        return self._db.get_features_at_location(self._name, key)

    def __iter__(self):
        return iter(self.get_features())

    def get_features(self):
        return sorted(
            (f for f in self._db.get_all_features() if f.chromosome == self._name),
            key=lambda f: f.start,
        )


class _Feature(object):
    """A base class for a genomic feature."""

    def __init__(
        self,
        components,
        primary_name_col,
        common_name_col,
        type_col,
        chrom_col,
        start_col,
        stop_col,
        strand_col,
        description_col,
    ):
        self.primary_name = self.standard_name = components[primary_name_col]
        self.common_name = components[common_name_col]
        self.name = self.common_name or self.standard_name
        self.all_names = set((self.name, self.primary_name))

        self.type = components[type_col]
        self.chromosome = components[chrom_col]
        self.strand = components[strand_col]
        self.description = components[description_col]

        if components[start_col] == "":
            self.start = self.stop = -1
        else:
            start = int(components[start_col])
            stop = int(components[stop_col])
            self.start, self.stop = min(start, stop), max(start, stop)

        self._len = self.stop - self.start + 1
        self.is_orf = "ORF" in self.type

        self._set_exons()
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
        return self.start + apparent_nt + offset - 1


class AlbicansFeature(_Feature):
    def __init__(self, feature_line):
        components = feature_line.split("\t")
        super().__init__(components, 0, 1, 3, 4, 5, 6, 7, 10)

        self.primary_cgdid = components[9 - 1]  # Example: CAL0000184345. Used in GO annotations.

        aliases = [x for x in components[2].split("|") if x]
        aliases.append(self.primary_cgdid)
        self.all_names |= set(aliases)

        self.cerevisiae_orthologs = frozenset(c.strip() for c in components[17].split("|") if c.strip())


class CerevisiaeFeature(_Feature):
    def __init__(self, feature_line):
        components = feature_line.split("\t")
        super().__init__(components, 0, 4, 1, 8, 9, 10, 11, 15)

        self.parent_feature = components[7 - 1]
        self.feature_name = components[4 - 1]
        self.all_names.add(self.feature_name)
        self.name = self.common_name or self.feature_name or self.standard_name
        self.feature_qualifier = components[3 - 1]


@Shared.memoized
def default_cer_db():
    return CerevisiaeFeatureDB(Shared.get_dependency("cerevisiae/SGD_features.tab"))


@Shared.memoized
def default_alb_db():
    # NOTE: underscore path: reference_genome (NOT "reference genome")
    return AlbicansFeatureDB(
        Shared.get_dependency(
            "albicans/reference_genome/C_albicans_SC5314_version_A22-s07-m01-r08_chromosomal_feature.tab"
        ),
        Shared.get_dependency(
            "albicans/reference_genome/C_albicans_SC5314_version_A22-s07-m01-r08_chromosomes.fasta"
        ),
        Shared.get_dependency(
            "albicans/reference_genome/C_albicans_SC5314_version_A22-s07-m01-r08_features.gff"
        ),
        Shared.get_dependency(
            "albicans/reference_genome/C_albicans_SC5314_iprscan.out"
        ),
    )


@Shared.memoized
def default_pom_db():
    return PombeFeatureDB(
        Shared.get_dependency("pombe/schizosaccharomyces_pombe.chr.gff3"),
        Shared.get_dependency("pombe/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa"),
    )


@Shared.memoized
def default_glab_db():
    return GlabrataFeatureDB(
        Shared.get_dependency("glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.gff"),
        Shared.get_dependency("glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.fna"),
        Shared.get_dependency("glabrata/reference_genome/interproscan6/GCA_014217725.1_ASM1421772v1_protein.faa.tsv"),
    )


if __name__ == "__main__":
    pom_db = default_pom_db()
    alb_db = default_alb_db()
    cer_db = default_cer_db()
    glab_db = default_glab_db()

    assert alb_db.get_features_at_range(1, (2000, 3000)) == []
    assert len(alb_db.get_features_at_range(1, (4000, 4400))) == 1
    assert alb_db.get_features_at_range(1, (4000, 4400))[0].standard_name == "C1_00010W_A"
    assert list(alb_db.get_interfeature_range(1, (4100, 4700))) == [(4398, 4408)]
    assert len(alb_db.get_features_at_range(1, (1, 10000))) == 3
    assert len(alb_db.get_features_at_range(5, (1150249, 1152390))) == 3
    assert [f.standard_name for f in alb_db.get_features_at_range(5, (1150249, 1152390))] == [
        "C5_05280C_A",
        "C5_05290C_A",
        "C5_05300W_A",
    ]

    assert list(cer_db.get_feature_by_name("YAL001C").exons) == [(147594, 151006), (151097, 151166)]
    assert cer_db.get_feature_by_name("YAL001C").coding_length == 3483

    tub2 = alb_db.get_feature_by_name("C1_00710C_A")
    assert tub2.aa_to_genomic_nt(1) == 125597 + 1 - 1
    assert tub2.aa_to_genomic_nt(2) == 125597 + 1 - 4
    assert tub2.aa_to_genomic_nt(4) == 125597 + 1 - 10
    assert tub2.aa_to_genomic_nt(5) == 125597 + 1 - 228
    assert tub2.aa_to_genomic_nt(16) == 125597 + 1 - 261
    assert tub2.aa_to_genomic_nt(17) == 125597 + 1 - 428
    assert tub2.aa_to_genomic_nt(400) == 125597 + 1 - 428 - ((400 - 17) * 3)

    rad6 = alb_db.get_feature_by_name("C7_03870W_A")
    assert rad6.aa_to_genomic_nt(1) == 855127 - 1 + 1
    assert rad6.aa_to_genomic_nt(2) == 855127 - 1 + 4
    assert rad6.aa_to_genomic_nt(14) == 855127 - 1 + 40
    assert rad6.aa_to_genomic_nt(15) == 855127 - 1 + 152
    assert rad6.aa_to_genomic_nt(36) == 855127 - 1 + 215
    assert rad6.aa_to_genomic_nt(37) == 855127 - 1 + 324
    assert rad6.aa_to_genomic_nt(137) == 855127 - 1 + 324 + 300

    ade6 = pom_db.get_feature_by_name("ade6")
    assert ade6.standard_name == "SPCC1322.13"
    assert ade6.start == 1316337
    assert ade6.stop == 1317995
    assert ade6.strand == "W"

    assert glab_db.get_std_chrom_name("CP048230.1") == "CP048230.1"
    assert glab_db.get_feature_by_name("GWK60_A00033") is not None
    # ============================================================

    print("\n--- Feature DB summary checks ---")

    def summarize_db(label, db, n_show_chroms=5):
        chroms = sorted(db._chrom_features.keys())
        n_features = len(db.get_all_features())

        print(f"\n[{label}]")
        print("  Total features:", n_features)
        print("  Chromosomes / contigs:", len(chroms))
        print("  Example contigs:", chroms[:n_show_chroms])

        if chroms:
            print("  Mean features per contig:",
                  round(n_features / float(len(chroms)), 2))

        # name-resolution sanity check
        test_feat = db.get_all_features()[0]
        test_name = next(iter(test_feat.all_names))
        assert db.get_feature_by_name(test_name) is test_feat, \
            f"Name lookup failed for {label}"

    summarize_db("S. cerevisiae", cer_db)
    summarize_db("C. albicans", alb_db)
    summarize_db("S. pombe", pom_db)
    summarize_db("C. glabrata (BG2)", glab_db)

    # Print some statistics:
    print("Genes with introns in Sc:", len([f for f in cer_db.get_all_features() if f.coding_length != len(f)]))
    print("Genes with introns in Calb:", len([f for f in alb_db.get_all_features() if f.coding_length != len(f)]))
    print("Genes with introns in Cgla:", len([f for f in glab_db.get_all_features() if f.coding_length != len(f)]))

    print("Genes with domains in Calb:", len([f for f in alb_db.get_all_features() if f.domains]))
    print("Genes with domains in Cgla:", len([f for f in glab_db.get_all_features() if f.domains]))

    print("Genes with more than one exon in Sp:", len([f for f in pom_db.get_all_features() if len(f.exons) > 1]))
    print("Genes with more than one exon in Cgla:", len([f for f in glab_db.get_all_features() if len(f.exons) > 1]))

    print(
        "Genes which are wider than their exons in Sp:",
        len(
            [
                f
                for f in pom_db.get_all_features()
                if f.start > min(e[0] for e in f.exons) or f.stop < max(e[1] for e in f.exons)
            ]
        ),
    )

