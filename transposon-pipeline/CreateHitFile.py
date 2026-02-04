
#!/usr/bin/env python3

from collections import defaultdict
import glob, os
import pandas as pd
import numpy as np
import pysam
import itertools
import Shared


usage = '''CreateHitFile.py  
   -i  --in-dir       [str]   Input directory with .bam files to parse. Defaults to current directory if left unspecified.
   -b  --bam          [str]   One sorted BAM to parse (repeatable). If provided, overrides directory scanning.
   -o  --out-dir      [str]   Output directory to which per-sample outputs will be written. Defaults to current directory if left unspecified.
   --out-prefix       [str]   Prefix for output files (per BAM). Defaults to BAM basename (without .sorted.bam/.bam).
   -q  --min-mapq     [int]   Map Quality - hits to parse from the bam file (default is 20)
   -m  --merge-dist   [int]   Hits to merge with at most x nt distance between two hits. Default is 2 
                                Example: Hits in positions 1 and 3  (3-1=2) will be united into a single hit
   -f  --fasta        [str]   REQUIRED. Reference genome FASTA used to derive chromosome lengths.
   -g  --features     [str]   REQUIRED. Feature table used to map hits to nearest ORFs/features.
   --feature-format   [str]   Feature file format: ncbi_feature_table | cgd_tab | gff | gtf  (default: ncbi_feature_table)
   --ncbi-chrom-field [str]   For ncbi_feature_table: genomic_accession | chromosome (default: genomic_accession)
   --ncbi-keep        [str]   For ncbi_feature_table: gene | cds | all (default: gene)
   -t  --threads      [int]   Number of CPU cores to use (default: 1)

   --write-gene-counts        Write GeneCount.csv and removedGeneCount.csv (Gale-style).
   --drop-top-sites   [int]   For removedGeneCount: remove top N sites per gene (default 1).
   -h  --help                 Show this help message and exit 
'''


ChrFeatCols = ['FeatureName', 'GeneName','Aliases','FeatureType','Chromosome','StartCoord','StopCoord','Strand','PrimaryCGDID','SecondaryCGDID',\
        'Description','DateCreated','SeqCoordVerDate','Blank1','Blank2','GeneNameReserDate','ReservedIsstandardName','SC_ortholog']


def FindHitsPerSample(SamAlign, ChrFeatMap, Sep0N = 2,MapQ=10):
    """Goes through Bam file, checks for high confidence alignment, unites unique positions if they can be aligned with adjunct positions."""
    Chromosomes = SamAlign.references
    hit_map = {}

    for line in SamAlign:
        assert line.tid == line.reference_id
        if line.mapq < MapQ:
            continue
        # Since start < end always, in alignments which are reversed (along the Crick strand) the start of the fragment is actually at the 'end' point. 
        if line.is_reverse:
            pos = line.reference_end
            source = 'C'
        else:
            # BAM files use 0-based indexing, and we work in 1-based indexing, so we have to add one.
            pos = line.reference_start + 1
            source = 'W'
        chrom = SamAlign.get_reference_name(line.reference_id)
        if chrom not in hit_map:
            hit_map[chrom] = {'W': np.zeros(ChrLen[chrom]+1, dtype=int),
                              'C': np.zeros(ChrLen[chrom]+1, dtype=int)}
        hit_map[chrom][source][pos] += 1
    
    ChrHitList = {}
    TotalHits = 0
    for (Chr, source) in itertools.product(Chromosomes, ('W', 'C')):
        if Chr not in ChrFeatMap:
            print('{} not found in Chromosome feature file'.format(Chr))
            continue
            
        PosCount = hit_map[Chr][source]
        HitsPos = np.where(PosCount>0)[0]
        TotalHits += len(HitsPos)

        if Chr not in ChrHitList:
            ChrHitList[Chr] = []
        i=0
        
        while i < len(HitsPos):
            StartI=i
            while (i < len(HitsPos)-1) and (HitsPos[i+1] - HitsPos[i] <= Sep0N):
                i+=1
            ChrHitList[Chr].append((HitsPos[StartI],
                                    HitsPos[i],
                                    HitsPos[i] - HitsPos[StartI],
                                    sum(PosCount[HitsPos[StartI]:HitsPos[i]+1]),
                                    Chr,
                                    source))
            i+=1
    TotalUniqueHits = sum(len(chrom_hits) for chrom_hits in ChrHitList.values())
    total_reads = sum(hits.sum() for source in hit_map.values() for hits in source.values())
    return ChrHitList, TotalHits, TotalUniqueHits, total_reads


class NearestORF():
    """Class gets the position of a hit and finds its nearest ORF in the specified strand (W/C) up or down stream of that index"""
    def __init__(self,Ind, Dist, Strand, UpDown):
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
            return "ORF(" + self.Strand +')'
        elif self.Ind < 0:
            return " - "
        else :
            return self.Strand


def FindNearestORFInStrand(StartI, StopI, ChrFeatMap, Strand):
    if not isinstance(ChrFeatMap[StartI], tuple):
        return NearestORF(ChrFeatMap[StartI],0, Strand,'Up'), NearestORF(ChrFeatMap[StartI],0, Strand,'Down')
    elif not isinstance(ChrFeatMap[StopI], tuple):
        return NearestORF(ChrFeatMap[StopI],0, Strand,'Up'), NearestORF(ChrFeatMap[StopI],0, Strand,'Down')
    else:
        upi=StartI-1
        upF =-1
        if not isinstance(ChrFeatMap[upi], tuple):
            upF = ChrFeatMap[upi]
        else:
            new_upi = ChrFeatMap[upi][0]
            if new_upi != -1:
                upi = new_upi
                upF = ChrFeatMap[upi]
            else:
                upi = 1

        downi=StopI+1
        DownF = -1
        if not isinstance(ChrFeatMap[downi], tuple):
            DownF = ChrFeatMap[downi]
        else:
            new_downi = ChrFeatMap[downi][1]
            if new_downi != -1:
                downi = new_downi 
                DownF = ChrFeatMap[downi]
            else:
                downi = len(ChrFeatMap)-1

        return NearestORF(upF,StartI - upi, Strand,'Up'), NearestORF(DownF, downi-StopI, Strand,'Down')


def _open_maybe_gz(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        import gzip
        return gzip.open(path, mode)
    return open(path, mode)


def read_fasta_lengths(fasta_path: str):
    ChrLen = {}
    name = None
    seqlen = 0
    with _open_maybe_gz(fasta_path, "rt") as f:
        for l in f:
            if not l:
                continue
            if l[0] == '>':
                if name is not None:
                    ChrLen[name] = seqlen
                header = l[1:].strip()
                name = header.split()[0]
                seqlen = 0
            else:
                seqlen += len(l.strip())
        if name is not None:
            ChrLen[name] = seqlen
    return ChrLen


def _ncbi_pick_name(row):
    for k in ("symbol", "name", "locus_tag", "GeneID"):
        v = row.get(k, "")
        if isinstance(v, str) and v.strip():
            return v.strip()
    return "NA"


def load_features(feature_path: str, feature_format: str,
                  ncbi_chrom_field: str = "genomic_accession",
                  ncbi_keep: str = "gene"):
    if feature_format == "ncbi_feature_table":
        df = pd.read_table(feature_path, dtype=str)
        df.columns = [c.strip().lstrip("#").strip() for c in df.columns]
        for col in ("feature", "start", "end", "strand"):
            if col not in df.columns:
                raise ValueError(f"NCBI feature table missing required column: {col}")

        if ncbi_keep == "gene":
            df = df[df["feature"] == "gene"].copy()
        elif ncbi_keep == "cds":
            df = df[df["feature"] == "CDS"].copy()
        else:
            df = df.copy()

        if ncbi_chrom_field not in df.columns:
            raise ValueError(f"--ncbi-chrom-field {ncbi_chrom_field} not found in table columns")
        chrom_col = ncbi_chrom_field

        df["start"] = pd.to_numeric(df["start"], errors="coerce")
        df["end"] = pd.to_numeric(df["end"], errors="coerce")
        df = df.dropna(subset=["start", "end", "strand", chrom_col])
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)

        def _to_wc(s):
            s = str(s).strip()
            if s == "+":
                return "W"
            if s == "-":
                return "C"
            return ""

        df["StrandWC"] = df["strand"].map(_to_wc)
        df = df[df["StrandWC"].isin(["W", "C"])].copy()

        rows = []
        for _, r in df.iterrows():
            feat_name = (r.get("locus_tag", "") or "").strip()
            if not feat_name:
                feat_name = _ncbi_pick_name(r)

            gene_name = _ncbi_pick_name(r)

            feature_type = (r.get("class", "") or r.get("feature", "") or "").strip()
            if not feature_type:
                feature_type = str(r.get("feature", "gene")).strip()

            chrom = str(r[chrom_col]).strip()

            s = int(min(r["start"], r["end"]))
            e = int(max(r["start"], r["end"]))

            rows.append({
                "FeatureName": feat_name,
                "GeneName": gene_name,
                "Aliases": "",
                "FeatureType": feature_type,
                "Chromosome": chrom,
                "StartCoord": s,
                "StopCoord": e,
                "Strand": r["StrandWC"]
            })

        ChrFeature = pd.DataFrame(rows)
        for c in ChrFeatCols:
            if c not in ChrFeature.columns:
                ChrFeature[c] = ""
        ChrFeature["StartCoord"] = pd.to_numeric(ChrFeature["StartCoord"], errors="coerce").astype(int)
        ChrFeature["StopCoord"] = pd.to_numeric(ChrFeature["StopCoord"], errors="coerce").astype(int)
        return ChrFeature

    raise ValueError("Only ncbi_feature_table is enabled in this universalized version for now.")


def derive_prefix_from_bam(bam_path: str) -> str:
    base = os.path.basename(bam_path)
    if base.endswith(".sorted.bam"):
        base = base[:-11]
    elif base.endswith(".bam"):
        base = base[:-4]
    return base


def ListHitProp_and_gene_counts(ChrHitList, HitFileName, ChrFeatC, ChrFeatW,
                               write_gene_counts=False,
                               gene_count_path=None,
                               removed_gene_count_path=None,
                               drop_top_sites=1):
    Featuref = open(HitFileName,'w')
    Featuref.write('Chromosome\tSource\tUp feature type\tUp feature name\tUp gene name\tUp feature dist\tDown feature type\tDown feature name\tDown gene name\tDown feature dist\tIntergenicType\tHit position\tHit count\n')

    # For Gale-style gene summaries (sum counts per gene; optionally drop top N site counts)
    gene_sites = defaultdict(list)  # feature_name -> list of site counts (each merged hit count)

    for Chr in ChrHitList.keys():
        for StartI,StopI,Len,Count,chr_,source in ChrHitList[Chr]:
            CORFUp,CORFDown = FindNearestORFInStrand(StartI,StopI, ChrFeatC[Chr], 'C')
            WORFUp,WORFDown = FindNearestORFInStrand(StartI,StopI, ChrFeatW[Chr], 'W')
            UpORF = CORFUp if CORFUp < WORFUp else WORFUp
            DownORF = CORFDown if CORFDown< WORFDown else WORFDown

            Featuref.write('\t'.join([Chr,source,UpORF.GetFeatureType(), UpORF.GetFeatureName(), str(UpORF.GetGeneName()), str(UpORF.Dist),\
                                      DownORF.GetFeatureType(), DownORF.GetFeatureName(), str(DownORF.GetGeneName()), str(DownORF.Dist),\
                                      (UpORF.GetStrand()+'-'+DownORF.GetStrand()),str(StartI), str(Count)])+'\n')

            if write_gene_counts:
                # Assign to a gene only if the hit overlaps an ORF in either strand-map
                # (legacy: Dist==0 indicates "in ORF")
                in_feat = None
                if UpORF.Dist == 0 and UpORF.GetFeatureName() != 'None':
                    in_feat = UpORF.GetFeatureName()
                elif DownORF.Dist == 0 and DownORF.GetFeatureName() != 'None':
                    in_feat = DownORF.GetFeatureName()

                if in_feat is not None:
                    gene_sites[in_feat].append(int(Count))

    Featuref.close()

    if write_gene_counts:
        # GeneCount: sum of all site counts per gene
        # removedGeneCount: sum after removing top N site counts per gene
        with open(gene_count_path, "w") as gc, open(removed_gene_count_path, "w") as rgc:
            for feat_name, counts in gene_sites.items():
                counts_sorted = sorted(counts, reverse=True)
                total = sum(counts_sorted)

                counts_removed = counts_sorted[:]
                n = int(drop_top_sites)
                if n < 0:
                    n = 0
                if n > 0 and len(counts_removed) > 0:
                    del counts_removed[:min(n, len(counts_removed))]
                total_removed = sum(counts_removed)

                gc.write(f"{feat_name} , {total}\n")
                rgc.write(f"{feat_name} , {total_removed}\n")


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(usage=usage)
    
    parser.add_argument("-o", "--out-dir", default='.')
    parser.add_argument("-i", "--in-dir", default='.')
    parser.add_argument("-b", "--bam", action="append", default=[])
    parser.add_argument("--out-prefix", default=None)

    parser.add_argument("-q", "--min-mapq", default=20, type=int)
    parser.add_argument("-m", "--merge-dist", default=2, type=int)

    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-g", "--features", required=True)

    parser.add_argument("--feature-format", default="ncbi_feature_table",
                        choices=["ncbi_feature_table"])
    parser.add_argument("--ncbi-chrom-field", default="genomic_accession",
                        choices=["genomic_accession", "chromosome"])
    parser.add_argument("--ncbi-keep", default="gene", choices=["gene", "cds", "all"])

    parser.add_argument("-t", "--threads", type=int, default=4)

    parser.add_argument("--write-gene-counts", action="store_true", default=False)
    parser.add_argument("--drop-top-sites", type=int, default=1)

    args = parser.parse_args()
    
    SamFileDir = args.in_dir
    OutDir = args.out_dir
    MapQ = args.min_mapq
    MergeDist = args.merge_dist
    Threads = args.threads

    if len(OutDir) > 0 and not os.path.isdir(OutDir): 
        os.makedirs(OutDir)

    # read chromosome lengths from fasta
    ChrLen = read_fasta_lengths(args.fasta)

    # read features
    ChrFeature = load_features(args.features, args.feature_format,
                               ncbi_chrom_field=args.ncbi_chrom_field,
                               ncbi_keep=args.ncbi_keep)
    Chromosomes = ChrFeature.Chromosome.unique()

    # Create chromosome feature maps (1-based arrays)
    ChrFeatW={} 
    for Chr in Chromosomes:
        if Chr not in ChrLen.keys():
            continue
        ChrFeatW[Chr] = -1*np.ones(ChrLen[Chr]+1)
        Feat = ChrFeature[(ChrFeature.Chromosome==Chr) & (ChrFeature.Strand == 'W') ]
        for row in Feat.iterrows():
            s = int(min(row[1].StartCoord, row[1].StopCoord))
            e = int(max(row[1].StartCoord, row[1].StopCoord))
            if s < 1:
                s = 1
            if e > ChrLen[Chr]:
                e = ChrLen[Chr]
            ChrFeatW[Chr][s: e+1] = row[0]

    ChrFeatC={}
    for Chr in Chromosomes:
        if Chr not in ChrLen.keys():
            continue
        ChrFeatC[Chr] = -1*np.ones(ChrLen[Chr]+1)
        Feat = ChrFeature[(ChrFeature.Chromosome==Chr) & (ChrFeature.Strand == 'C')]
        for row in Feat.iterrows():
            s = int(min(row[1].StartCoord, row[1].StopCoord))
            e = int(max(row[1].StartCoord, row[1].StopCoord))
            if s < 1:
                s = 1
            if e > ChrLen[Chr]:
                e = ChrLen[Chr]
            ChrFeatC[Chr][s: e+1] = row[0]

    # Preprocess the ChrFeat maps to include (left, right) tuples of the nearest features for every index:
    for feat_map in (ChrFeatC, ChrFeatW):
        for chrom_name in feat_map:
            chrom = list(map(int, feat_map[chrom_name]))
            prev_feat_ix = next_feat_ix = -1
            i = 1
            while i < len(chrom):
                if chrom[i] != -1:
                    prev_feat_ix = i
                    i += 1
                else:
                    next_feat_ix = -1
                    j = i+1
                    while j < len(chrom):
                        if chrom[j] != -1:
                            next_feat_ix = j
                            break
                        j += 1
                    for k in range(i, j):
                        chrom[k] = (prev_feat_ix, next_feat_ix)
                    i = j
            feat_map[chrom_name] = chrom

    # Determine BAMs to process
    bam_list = []
    if args.bam and len(args.bam) > 0:
        bam_list = args.bam[:]
    else:
        fNames = glob.glob(os.path.join(SamFileDir, '*.bam'))
        for Name in fNames:
            BaseName = os.path.basename(Name)
            if BaseName.find(".sorted") == -1:
                continue
            bam_list.append(Name)

    for bam_path in bam_list:
        if not os.path.isfile(bam_path):
            print(f"ERROR: BAM not found: {bam_path}")
            continue

        prefix = args.out_prefix if args.out_prefix is not None else derive_prefix_from_bam(bam_path)

        # per-replicate subdir (so later aggregation is trivial)
        rep_dir = os.path.join(OutDir, prefix)
        if not os.path.isdir(rep_dir):
            os.makedirs(rep_dir)

        OutHitFile = os.path.join(rep_dir, prefix + '_hits.txt')
        if os.path.isfile(OutHitFile):
            print('ERROR: Hit file for %s already exists.' % (prefix))
            continue

        Sami = pysam.AlignmentFile(bam_path, "rb")
        print("Parsing file %s..." % (prefix))
        
        MyList,TotalHits, TotalUniqueHits, TotalReads = FindHitsPerSample(Sami,ChrFeatW, MergeDist,MapQ)

        gene_count_path = os.path.join(rep_dir, prefix + "_GeneCount.csv")
        removed_gene_count_path = os.path.join(rep_dir, "removed" + prefix + "_GeneCount.csv")

        ListHitProp_and_gene_counts(
            MyList,
            OutHitFile,
            ChrFeatC, ChrFeatW,
            write_gene_counts=args.write_gene_counts,
            gene_count_path=gene_count_path,
            removed_gene_count_path=removed_gene_count_path,
            drop_top_sites=args.drop_top_sites
        )

        UniqueHitPercent = round(float(TotalUniqueHits)/float(TotalHits)*100, 2)
        
        Log = '\r\n=== Finding hits ===\r\n%s reads found of map quality >= %s; in these:\r\n  %s hits were found; of these:\r\n    %s (%s%%) hit positions were found to be unique (Minimal distance = %s)\r\n' % (TotalReads, MapQ, TotalHits, TotalUniqueHits, UniqueHitPercent, MergeDist)
        print(Log)
        
        LogFile = open('%s_log.txt' % (prefix), 'a')
        LogFile.write(Log)
        LogFile.close()

