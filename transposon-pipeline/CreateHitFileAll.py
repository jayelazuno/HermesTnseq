#!/usr/bin/env python3
import glob, os
import pandas as pd
import numpy as np
import pysam
import itertools

usage = """CreateHitsAll.py
  -i  --in-dir       [str]   Input directory with .bam files to parse. Defaults to current directory.
  -o  --out-dir      [str]   Output directory to which the hit file will be written. Defaults to current directory.
  -q  --min-mapq     [int]   Map Quality - hits to parse from the bam file (default is 20)
  -m  --merge-dist   [int]   Hits to merge with at most x nt distance between two hits. Default is 2
                                Example: Hits in positions 1 and 3 (3-1=2) will be united into a single hit

  --chr-fasta        [str]   REQUIRED. Reference FASTA (for chromosome lengths)
  --feature-tab      [str]   REQUIRED. CGD chromosomal_feature.tab (for ORF locations)
  --feature-skiprows [int]   Skip rows for feature table (default 8; matches original)
"""

# Same columns as original script
ChrFeatCols = [
    'FeatureName', 'GeneName', 'Aliases', 'FeatureType', 'Chromosome',
    'StartCoord', 'StopCoord', 'Strand', 'PrimaryCGDID', 'SecondaryCGDID',
    'Description', 'DateCreated', 'SeqCoordVerDate', 'Blank1', 'Blank2',
    'GeneNameReserDate', 'ReservedIsstandardName', 'SC_ortholog'
]

def read_chr_len_from_fasta_headers(fasta_path):
    """
    Match the original script behavior:
    expects headers like:
      >Ca22chr1A_C_albicans_SC5314 (3188363 nucleotides)
    Falls back to .fai if present, else sequence scan.
    """
    lengths = {}

    # Prefer .fai if it exists
    fai = fasta_path + ".fai"
    if os.path.exists(fai):
        with open(fai, "r") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    lengths[parts[0]] = int(parts[1])
        if lengths:
            return lengths

    # Try header-parsing like the original
    try:
        with open(fasta_path, "r") as f:
            for l in f:
                if not l:
                    continue
                if l[0] == '>':
                    # Example: >Ca22chr1A_C_albicans_SC5314 (3188363 nucleotides)
                    chr_name = l[1:l.find(' ')] if ' ' in l else l[1:].strip()
                    if '(' in l and 'nuc' in l:
                        try:
                            L = int(l[l.find('(') + 1:l.find('nuc')].strip())
                            lengths[chr_name] = L
                        except Exception:
                            pass
        if lengths:
            return lengths
    except Exception:
        pass

    # Fallback: scan sequence lengths
    name = None
    seq_len = 0
    with open(fasta_path, "r") as f:
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

def FindHitsPerSample(SamAlign, ChrFeatMap, Sep0N=2, MapQ=10):
    """
    Ported 1:1 from original, except Python3 syntax.
    Uses global ChrLen (like original).
    """
    Chromosomes = SamAlign.references
    hit_map = {}

    for line in SamAlign:
        # original had: assert line.tid == line.reference_id
        if line.is_unmapped:
            continue
        if line.mapq < MapQ:
            continue

        # reverse fragment "start" is at reference_end
        if line.is_reverse:
            pos = line.reference_end
            source = 'C'
        else:
            pos = line.reference_start + 1
            source = 'W'

        chrom = SamAlign.get_reference_name(line.reference_id)

        if chrom not in hit_map:
            hit_map[chrom] = {
                'W': np.zeros(ChrLen[chrom] + 1, dtype=int),
                'C': np.zeros(ChrLen[chrom] + 1, dtype=int)
            }
        hit_map[chrom][source][pos] += 1

    ChrHitList = {}
    TotalHits = 0

    for (Chr, source) in itertools.product(Chromosomes, ('W', 'C')):
        if Chr not in ChrFeatMap:
            print(f"{Chr} not found in Chromosome feature file")
            continue

        PosCount = hit_map[Chr][source]
        HitsPos = np.where(PosCount > 0)[0]  # positions with >=1 read
        TotalHits += len(HitsPos)

        if Chr not in ChrHitList:
            ChrHitList[Chr] = []
        i = 0

        while i < len(HitsPos):
            StartI = i
            while (i < len(HitsPos) - 1) and (HitsPos[i + 1] - HitsPos[i] <= Sep0N):
                i += 1
            ChrHitList[Chr].append((
                int(HitsPos[StartI]),
                int(HitsPos[i]),
                int(HitsPos[i] - HitsPos[StartI]),
                int(PosCount[HitsPos[StartI]:HitsPos[i] + 1].sum()),
                Chr,
                source
            ))
            i += 1

    TotalUniqueHits = sum(len(chrom_hits) for chrom_hits in ChrHitList.values())
    total_reads = sum(int(hits.sum()) for source in hit_map.values() for hits in source.values())
    return ChrHitList, TotalHits, TotalUniqueHits, total_reads

class NearestORF:
    def __init__(self, Ind, Dist, Strand, UpDown):
        self.Ind = Ind
        self.Dist = Dist
        self.Strand = Strand
        self.UpDown = UpDown

    def __gt__(self, ORF2):
        return self.Dist > ORF2.Dist

    def __lt__(self, ORF2):
        return self.Dist < ORF2.Dist

    def GetFeatureName(self):
        if self.Ind >= 0:
            return ChrFeature.loc[self.Ind].FeatureName
        else:
            return 'None'

    def GetGeneName(self):
        if self.Ind >= 0:
            return ChrFeature.loc[self.Ind].GeneName
        else:
            return 'None'

    def GetFeatureType(self):
        if self.Ind >= 0:
            return ChrFeature.loc[self.Ind].FeatureType
        else:
            return 'None'

    def GetStrand(self):
        if self.Dist == 0:
            return "ORF(" + self.Strand + ')'
        elif self.Ind < 0:
            return " - "
        else:
            return self.Strand

def FindNearestORFInStrand(StartI, StopI, ChrFeatMap, Strand):
    # ORF exists at hit position, first merged hit
    if not isinstance(ChrFeatMap[StartI], tuple):
        return NearestORF(int(ChrFeatMap[StartI]), 0, Strand, 'Up'), NearestORF(int(ChrFeatMap[StartI]), 0, Strand, 'Down')

    # ORF exists at hit position, last merged hit
    elif not isinstance(ChrFeatMap[StopI], tuple):
        return NearestORF(int(ChrFeatMap[StopI]), 0, Strand, 'Up'), NearestORF(int(ChrFeatMap[StopI]), 0, Strand, 'Down')

    else:
        # upstream is down the index
        upi = StartI - 1
        upF = -1
        if not isinstance(ChrFeatMap[upi], tuple):
            upF = int(ChrFeatMap[upi])
        else:
            new_upi = ChrFeatMap[upi][0]
            if new_upi != -1:
                upi = int(new_upi)
                upF = int(ChrFeatMap[upi])
            else:
                upi = 1

        # downstream
        downi = StopI + 1
        DownF = -1
        if not isinstance(ChrFeatMap[downi], tuple):
            DownF = int(ChrFeatMap[downi])
        else:
            new_downi = ChrFeatMap[downi][1]
            if new_downi != -1:
                downi = int(new_downi)
                DownF = int(ChrFeatMap[downi])
            else:
                downi = len(ChrFeatMap)

        return NearestORF(upF, StartI - upi, Strand, 'Up'), NearestORF(DownF, downi - StopI, Strand, 'Down')

def ListHitProp(ChrHitList, FileName, ChrFeatC, ChrFeatW):
    with open(FileName, 'w') as Featuref:
        Featuref.write(
            'Chromosome\tSource\tUp feature type\tUp feature name\tUp gene name\tUp feature dist\t'
            'Down feature type\tDown feature name\tDown gene name\tDown feature dist\t'
            'IntergenicType\tHit position\tHit count\n'
        )
        for Chr in ChrHitList.keys():
            for StartI, StopI, Len, Count, chr_, source in ChrHitList[Chr]:
                CORFUp, CORFDown = FindNearestORFInStrand(StartI, StopI, ChrFeatC[Chr], 'C')
                WORFUp, WORFDown = FindNearestORFInStrand(StartI, StopI, ChrFeatW[Chr], 'W')
                UpORF = CORFUp if CORFUp < WORFUp else WORFUp
                DownORF = CORFDown if CORFDown < WORFDown else WORFDown

                Featuref.write('\t'.join([
                    Chr, source,
                    UpORF.GetFeatureType(), UpORF.GetFeatureName(), str(UpORF.GetGeneName()), str(UpORF.Dist),
                    DownORF.GetFeatureType(), DownORF.GetFeatureName(), str(DownORF.GetGeneName()), str(DownORF.Dist),
                    (UpORF.GetStrand() + '-' + DownORF.GetStrand()),
                    str(StartI), str(Count)
                ]) + '\n')

def _strip_sorted_suffix(basename):
    """
    Original script expected *.bam.sorted and did:
        BaseName = BaseName[:-7]  (remove ".sorted")
        PrefixName = BaseName[:-4] (remove ".bam")
    We support both:
        foo.bam.sorted  -> PrefixName=foo
        foo.sorted.bam  -> PrefixName=foo
    """
    if basename.endswith(".bam.sorted"):
        base = basename[:-7]          # remove ".sorted"
        prefix = base[:-4]            # remove ".bam"
        return prefix
    if basename.endswith(".sorted.bam"):
        prefix = basename[:-11]       # remove ".sorted.bam"
        return prefix
    # fallback: mimic original behavior if ".sorted" appears anywhere
    if ".sorted" in basename:
        base = basename.replace(".sorted", "")
        if base.endswith(".bam"):
            return base[:-4]
        return os.path.splitext(base)[0]
    return None

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("-o", "--out-dir", default='.')
    parser.add_argument("-i", "--in-dir", default='.')
    parser.add_argument("-q", "--min-mapq", default=20, type=int)
    parser.add_argument("-m", "--merge-dist", default=2, type=int)

    parser.add_argument("--chr-fasta", required=True)
    parser.add_argument("--feature-tab", required=True)
    parser.add_argument("--feature-skiprows", default=8, type=int)

    args = parser.parse_args()

    SamFileDir = args.in_dir
    GeneListFileDir = args.out_dir
    MapQ = args.min_mapq
    MergeDist = args.merge_dist

    # ---- read chromosome lengths (ChrLen global like original) ----
    ChrLen = read_chr_len_from_fasta_headers(args.chr_fasta)

    # ---- read all features (ChrFeature global like original) ----
    ChrFeature = pd.read_table(args.feature_tab, skiprows=args.feature_skiprows, names=ChrFeatCols)
    Chromosomes = ChrFeature.Chromosome.unique()

    # ---- build ChrFeat maps (W and C), same structure as original ----
    ChrFeatW = {}
    for Chr in Chromosomes:
        if Chr not in ChrLen:
            continue
        ChrFeatW[Chr] = -1 * np.ones(ChrLen[Chr] + 1, dtype=int)
        Feat = ChrFeature[(ChrFeature.Chromosome == Chr) & (ChrFeature.Strand == 'W')]
        for row in Feat.iterrows():
            ChrFeatW[Chr][int(row[1].StartCoord): int(row[1].StopCoord) + 1] = row[0]

    ChrFeatC = {}
    for Chr in Chromosomes:
        if Chr not in ChrLen:
            continue
        ChrFeatC[Chr] = -1 * np.ones(ChrLen[Chr] + 1, dtype=int)
        Feat = ChrFeature[(ChrFeature.Chromosome == Chr) & (ChrFeature.Strand == 'C')]
        for row in Feat.iterrows():
            ChrFeatC[Chr][int(row[1].StopCoord): int(row[1].StartCoord) + 1] = row[0]

    # ---- preprocess maps into (left,right) tuples, same logic as original ----
    for feat_map in (ChrFeatC, ChrFeatW):
        for chrom_name in feat_map:
            chrom = list(map(int, feat_map[chrom_name]))  # Python3: map -> list
            prev_feat_ix = next_feat_ix = -1
            i = 1
            while i < len(chrom):
                if chrom[i] != -1:
                    prev_feat_ix = i
                    i += 1
                else:
                    next_feat_ix = -1
                    j = i + 1
                    while j < len(chrom):
                        if chrom[j] != -1:
                            next_feat_ix = j
                            break
                        j += 1
                    for k in range(i, j):
                        chrom[k] = (prev_feat_ix, next_feat_ix)
                    i = j
            feat_map[chrom_name] = chrom

    if len(GeneListFileDir) > 0 and not os.path.isdir(GeneListFileDir):
        os.makedirs(GeneListFileDir)

    fNames = glob.glob(os.path.join(SamFileDir, '*.bam'))

    for Name in fNames:
        BaseName = os.path.basename(Name)
        if ".sorted" not in BaseName:
            continue

        PrefixName = _strip_sorted_suffix(BaseName)
        if PrefixName is None:
            continue

        OutFileName = os.path.join(GeneListFileDir, PrefixName + '_hits.txt')
        if os.path.isfile(OutFileName):
            print(f'ERROR: Hit file for {BaseName} already exists.')
            continue

        Sami = pysam.AlignmentFile(Name, "rb")
        print(f"Parsing file {BaseName}...")

        MyList, TotalHits, TotalUniqueHits, TotalReads = FindHitsPerSample(Sami, ChrFeatW, MergeDist, MapQ)
        ListHitProp(MyList, OutFileName, ChrFeatC, ChrFeatW)

        UniqueHitPercent = round(float(TotalUniqueHits) / float(TotalHits) * 100, 2) if TotalHits else 0.0

        Log = (
            f"\r\n=== Finding hits ===\r\n"
            f"{TotalReads} reads found of map quality >= {MapQ}; in these:\r\n"
            f"  {TotalHits} hits were found; of these:\r\n"
            f"    {TotalUniqueHits} ({UniqueHitPercent}%) hit positions were found to be unique (Minimal distance = {MergeDist})\r\n"
        )
        print(Log)

        LogFile = open(os.path.join(GeneListFileDir, f'{PrefixName}_log.txt'), 'a')
        LogFile.write(Log)
        LogFile.close()
