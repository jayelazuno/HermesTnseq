#!/usr/bin/env python3

import os
from pathlib import Path
from collections import defaultdict

# Snakemake variables
forward_fasta = snakemake.input.forward
reverse_fasta = snakemake.input.reverse
genome_sizes_file = snakemake.input.sizes
output_bed_all = snakemake.output[0]

species = snakemake.params.species
qscore_cutoff = int(snakemake.params.qscore)
kmer_size = int(snakemake.params.kmer)
output_dir = Path(snakemake.params.outdir) / f"{species}_GenomeFiles"

# Output files
bed_forward = output_dir / f"{species}{qscore_cutoff}UnmappableForward.bed"
bed_reverse = output_dir / f"{species}{qscore_cutoff}UnmappableReverse.bed"
sam_forward = output_dir / f"Mapped_{kmer_size}Kmers_{species}.sam"
sam_reverse = output_dir / f"Mapped_Reverse_{kmer_size}Kmers_{species}.sam"
kmer_forward_fa = output_dir / f"Kmers_{forward_fasta.name}"
kmer_reverse_fa = output_dir / f"Kmers_{reverse_fasta.name}"

os.makedirs(output_dir, exist_ok=True)

# -------------------------
# Step 1: Generate K-mers
# -------------------------

def generate_kmers(fasta_path, output_path, k):
    with open(fasta_path, "r") as f_in, open(output_path, "w") as f_out:
        for line in f_in:
            if line.startswith(">"):
                chrom = line.strip()[1:]
                continue
            seq = line.strip()
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                name = f">{chrom}_(pos){i+1}"
                f_out.write(f"{name}\n{kmer}\n")

print("🧬generating k-mers...")
generate_kmers(forward_fasta, kmer_forward_fa, kmer_size)
generate_kmers(reverse_fasta, kmer_reverse_fa, kmer_size)

# -------------------------
# Step 2: Bowtie2 Mapping
# -------------------------

def run_bowtie2(input_fa, output_sam, index_prefix):
    cmd = f"bowtie2 -p 4 -x {index_prefix} -f {input_fa} -S {output_sam}"
    print(f"🧬 Running: {cmd}")
    os.system(cmd)

forward_index = output_dir / "Forward" / "Forward"
reverse_index = output_dir / "Reverse" / "Reverse"

print(" Mapping forward k-mers...")
run_bowtie2(kmer_forward_fa, sam_forward, forward_index)

print("📌Mapping reverse k-mers...")
run_bowtie2(kmer_reverse_fa, sam_reverse, reverse_index)

# -------------------------
# Step 3: Parse SAM files to BED
# -------------------------

# Load chromosome sizes
chrom_sizes = {}
with open(genome_sizes_file, "r") as f:
    for line in f:
        chrom, size = line.strip().split("\t")
        chrom_sizes[chrom] = int(size)

def parse_sam(sam_file, strand, bed_file, reverse=False):
    with open(sam_file, "r") as f_in, open(bed_file, "w") as f_out:
        for line in f_in:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = fields[3]
            mapq = int(fields[4])
            mdtag = [f for f in fields if f.startswith("MD:Z:")]
            mdtag = mdtag[0] if mdtag else "MD:Z:0"

            # Extract intended chrom and position from qname
            if "_(pos)" not in qname:
                continue
            chrom_q = qname.split("_(pos)")[0]
            loc_q = int(qname.split("_(pos)")[1])

            # Conditions for unmappability
            is_mismatch = (mdtag == "MD:Z:0")
            is_lowq = (mapq <= qscore_cutoff)
            is_wrong_chrom = (rname != "*" and rname != chrom_q)
            is_wrong_pos = (rname != "*" and int(pos) != loc_q)
            is_unmapped = (rname == "*")
            is_reversed = (flag == 16)

            if any([is_mismatch, is_lowq, is_wrong_chrom, is_wrong_pos, is_unmapped, is_reversed]):
                chrom_out = chrom_q
                pos_out = loc_q
                if reverse:
                    genome_len = chrom_sizes.get(chrom_out, None)
                    if genome_len:
                        pos_out = abs(genome_len - loc_q) + 1
                f_out.write(f"{chrom_out}\t{pos_out}\t{pos_out}\t.\t500\t{strand}\n")

print("🧬 Parsing SAM files to BED...")
parse_sam(sam_forward, "+", bed_forward, reverse=False)
parse_sam(sam_reverse, "-", bed_reverse, reverse=True)

# Merge into output_bed_all
with open(output_bed_all, "w") as out_all, open(bed_forward, "r") as f1, open(bed_reverse, "r") as f2:
    out_all.writelines(f1.readlines())
    out_all.writelines(f2.readlines())

print("All unmappables written to:", output_bed_all)

