#!/usr/bin/env python3
"""
MapFastq.py

(1) trim 5' transposon/primer (+tail) using cutadapt (optionally discard-untrimmed)
(2) optionally trim 3' adapter
(3) bowtie2 align -> SAM -> BAM -> sorted BAM + index
(4) logs the key stats
"""

import os
import sys
import gzip
import shutil
import argparse
from pathlib import Path
from subprocess import PIPE, run

import Shared  # keep dependency + behavior (no organism-specific calls here)

usage = """MapFastq.py
   -o  --out-dir            [str]   Output directory. Defaults to current directory if left unspecified.
   -i  --input-file-name    [str]   REQUIRED (SE only). Input fastq file (Need to include path and .fastq.gz at end of filename)
   -a  --clean-adapters     [str]   Trim 3' adapter sequence provided by user (optional).
   -r  --reverse-strand             Search with reverse complement sequence for R2 files.
                                     (Adaptor cleaning works the same as with R1.)
   -d  --delete-originals           Delete input FASTQ files.
   -k  --keep-clean-fqs             Keep the cleaned FASTQ files.
   -p  --primer-check               Check primer specificity if percent transposon in reads is low.
   -t  --threads             [int]  Threads for bowtie2/samtools (default 32)

   -x  --bt2-index          [str]   REQUIRED. Bowtie2 index prefix (the value you would pass to bowtie2 -x).
                                     Can be a full path to the index prefix.
                                     NOTE: bowtie2 requires spaces to be escaped with a backslash for -x.

   --r1                    [str]   (PE only) Read 1 FASTQ(.gz)
   --r2                    [str]   (PE only) Read 2 FASTQ(.gz)
   --junction-mate         [str]   (PE only) {R1,R2,both}; which mate defines junction/insertion (default R1)

   --trim-5p-seq           [str]   5' transposon/primer sequence to trim (optional).
   --trim-5p-seq-r2        [str]   (PE only) 5' sequence for R2 (only used when trimming both mates)
   --trim-5p-overlap       [int]   cutadapt --overlap for 5' trimming (default 24 if -p else 37)
   --discard-untrimmed             If set, discard reads lacking the 5' sequence match.

   --mapq                  [int]   MAPQ cutoff for BAM QC filtering (default 20)
   --drop-flags            [str]   samtools -F flags to drop (default 0xD04)
"""

# Example sequences (Mariner/Tn5):
# TnPrimerAndTail = "GTATTTTACCGACCGTTACCGACCGTTTTCATCCCTA"
# TnRev           = "TAGGGATGAAAACGGTCGGTAACGGTCGGTAAAATAC"
# PrimerOnly      = "GTATTTTACCGACCGTTACCGACC"
# PrimerRev       = "GGTCGGTAACGGTCGGTAAAATAC"
# AdapterSeq      = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"


def which_or_die(prog: str) -> str:
    p = shutil.which(prog)
    if not p:
        raise RuntimeError(
            f"ERROR: '{prog}' not found on PATH.\n"
            f"Activate the tnseq env (or install '{prog}') and try again."
        )
    return p


def run_cmd(command: str):
    p = run(command, shell=True, stdout=PIPE, stderr=PIPE, text=True)
    return p.returncode, p.stdout, p.stderr


def GetReads(val: str, text: str) -> int:
    idx = text.find(val)
    if idx < 0:
        return -1
    sub = text[idx + len(val) + 1 :].lstrip(" ")
    sub = sub[: sub.find("\n")]
    vals = sub.split(" ")
    try:
        return int(vals[0].replace(",", ""))
    except Exception:
        return -1


def CutAdaptOutput(output: str):
    if "=== Summary ===" in output:
        summary = output[output.find("=== Summary ===") + len("=== Summary ===") :]
    else:
        summary = output

    total = GetReads("Total reads processed", summary)
    if total < 0:
        total = GetReads("Total read pairs processed", summary)

    written = GetReads("Reads written (passing filters)", summary)
    if written < 0:
        written = GetReads("Reads written", summary)

    if total and total > 0 and written >= 0:
        percent = round((float(written) / float(total)) * 100.0, 2)
    else:
        percent = 0.0
    return total, written, percent


def RemoveTn_SE(InputFName: str, Trim5Seq: str, overlap: int, CleanfName: str,
                CutAdaptPath: str, discard_untrimmed: bool):
    discard_flag = " --discard-untrimmed" if discard_untrimmed else ""
    cmd = (
        f'{CutAdaptPath} -m 2 -g {Trim5Seq} '
        f'-o "{CleanfName}" "{InputFName}"{discard_flag} --overlap {overlap}'
    )
    rc, out, err = run_cmd(cmd)
    return CutAdaptOutput(out + "\n" + err), out, err, rc


def RemoveTn_PE(R1: str, R2: str,
                Trim5Seq_R1: str, Trim5Seq_R2: str,
                overlap: int,
                CleanR1: str, CleanR2: str,
                CutAdaptPath: str,
                discard_untrimmed: bool,
                trim_mate_mode: str):
    discard_flag = " --discard-untrimmed" if discard_untrimmed else ""
    if trim_mate_mode == "R1":
        cmd = (
            f'{CutAdaptPath} -m 2 -g {Trim5Seq_R1} '
            f'-o "{CleanR1}" -p "{CleanR2}" "{R1}" "{R2}"'
            f'{discard_flag} --overlap {overlap}'
        )
    elif trim_mate_mode == "R2":
        cmd = (
            f'{CutAdaptPath} -m 2 -G {Trim5Seq_R2} '
            f'-o "{CleanR1}" -p "{CleanR2}" "{R1}" "{R2}"'
            f'{discard_flag} --overlap {overlap}'
        )
    else:
        cmd = (
            f'{CutAdaptPath} -m 2 -g {Trim5Seq_R1} -G {Trim5Seq_R2} '
            f'-o "{CleanR1}" -p "{CleanR2}" "{R1}" "{R2}"'
            f'{discard_flag} --overlap {overlap}'
        )
    rc, out, err = run_cmd(cmd)
    return CutAdaptOutput(out + "\n" + err), out, err, rc


def TrimAdapter3_SE(InputFName: str, Adapter: str, CleanfName: str, CutAdaptPath: str):
    cmd = f'{CutAdaptPath} -m 2 -a {Adapter} -o "{CleanfName}" "{InputFName}"'
    rc, out, err = run_cmd(cmd)
    return CutAdaptOutput(out + "\n" + err), out, err, rc


def TrimAdapter3_PE(R1: str, R2: str, Adapter: str, CleanR1: str, CleanR2: str, CutAdaptPath: str):
    cmd = f'{CutAdaptPath} -m 2 -a {Adapter} -A {Adapter} -o "{CleanR1}" -p "{CleanR2}" "{R1}" "{R2}"'
    rc, out, err = run_cmd(cmd)
    return CutAdaptOutput(out + "\n" + err), out, err, rc


def AlignFastq_SE(BowtiePath: str, bt2_index: str, CleanfName: str, SamfName: str, threads: int, LogFile):
    cmd = (
        f'{BowtiePath} -p {threads} -x "{bt2_index}" '
        f'-U "{CleanfName}" -S "{SamfName}"'
    )
    rc, out, err = run_cmd(cmd)
    log = "\r\n\r\n=== Sequence alignment ===\r\n" + err.replace("\n", "\r\n")
    LogFile.write(log)
    if rc != 0:
        raise RuntimeError(f"bowtie2 failed (exit {rc}). See log for details.")


def AlignFastq_PE(BowtiePath: str, bt2_index: str, R1: str, R2: str, SamfName: str, threads: int, LogFile):
    cmd = (
        f'{BowtiePath} -p {threads} -x "{bt2_index}" '
        f'-1 "{R1}" -2 "{R2}" -S "{SamfName}"'
    )
    rc, out, err = run_cmd(cmd)
    log = "\r\n\r\n=== Sequence alignment ===\r\n" + err.replace("\n", "\r\n")
    LogFile.write(log)
    if rc != 0:
        raise RuntimeError(f"bowtie2 failed (exit {rc}). See log for details.")


def derive_prefix(fq_path: str) -> str:
    p = Path(fq_path)
    name = p.name
    if name.endswith(".gz"):
        name = name[:-3]
    for suf in (".fastq", ".fq"):
        if name.endswith(suf):
            name = name[: -len(suf)]
    return name


def copy_to_clean_se(src_path: str, dst_path: str):
    if src_path.endswith(".gz"):
        with gzip.open(src_path, "rb") as src, open(dst_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copyfile(src_path, dst_path)


def main():
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("-o", "--out-dir", default=".")
    parser.add_argument("-i", "--input-file-name", default=None)
    parser.add_argument("-x", "--bt2-index", required=True)

    parser.add_argument("--r1", default=None)
    parser.add_argument("--r2", default=None)
    parser.add_argument("--junction-mate", default="R1", choices=["R1", "R2", "both"])

    parser.add_argument("-a", "--clean-adapters", nargs="?", const="", default=None)
    parser.add_argument("-r", "--reverse-strand", action="store_true")
    parser.add_argument("-d", "--delete-originals", action="store_true")
    parser.add_argument("-k", "--keep-clean-fastqs", action="store_true")
    parser.add_argument("-p", "--primer_check", action="store_true")
    parser.add_argument("-t", "--threads", type=int, default=32)

    parser.add_argument("--trim-5p-seq", default=None)
    parser.add_argument("--trim-5p-seq-r2", default=None)
    parser.add_argument("--trim-5p-overlap", type=int, default=None)
    parser.add_argument("--discard-untrimmed", action="store_true")

    parser.add_argument("--mapq", type=int, default=20)
    parser.add_argument("--drop-flags", default="0xD04")

    args = parser.parse_args()

    pe_mode = (args.r1 is not None) or (args.r2 is not None)
    if pe_mode:
        if not (args.r1 and args.r2):
            raise RuntimeError("ERROR: PE mode requires both --r1 and --r2.")
        if args.input_file_name is not None:
            raise RuntimeError("ERROR: In PE mode, do not use -i/--input-file-name (use --r1/--r2).")
        if args.reverse_strand:
            raise RuntimeError("ERROR: -r/--reverse-strand is not used in PE mode (use --junction-mate).")
    else:
        if args.input_file_name is None:
            raise RuntimeError("ERROR: SE mode requires -i/--input-file-name.")
        if args.junction_mate != "R1":
            raise RuntimeError("ERROR: --junction-mate is PE-only; use SE input + -r if needed.")

    os.makedirs(args.out_dir, exist_ok=True)

    CutAdaptPath = which_or_die("cutadapt")
    BowtiePath = which_or_die("bowtie2")
    which_or_die("samtools")

    if pe_mode:
        prefix = derive_prefix(args.r1)
    else:
        prefix = derive_prefix(args.input_file_name)

    log_path = os.path.join(args.out_dir, f"{prefix}_log.txt")
    LogFile = open(log_path, "w+", encoding="utf-8")

    # trimming parameters with cutadapt

    overlap_default = 24 if args.primer_check else 37
    overlap = args.trim_5p_overlap if args.trim_5p_overlap is not None else overlap_default

    if pe_mode:
        CleanR1 = os.path.join(args.out_dir, f"{prefix}.R1.clean.fq")
        CleanR2 = os.path.join(args.out_dir, f"{prefix}.R2.clean.fq")

        if args.trim_5p_seq is None:
            LogFile.write("=== Trimming step skipped ===\r\n")
            copy_to_clean_se(args.r1, CleanR1)
            copy_to_clean_se(args.r2, CleanR2)
        else:
            trim_mate_mode = args.junction_mate
            if trim_mate_mode == "both":
                if args.trim_5p_seq_r2 is None:
                    raise RuntimeError("ERROR: --trim-5p-seq-r2 is required when --junction-mate both.")
                Trim5_R1 = args.trim_5p_seq
                Trim5_R2 = args.trim_5p_seq_r2
            elif trim_mate_mode == "R1":
                Trim5_R1 = args.trim_5p_seq
                Trim5_R2 = args.trim_5p_seq_r2 if args.trim_5p_seq_r2 is not None else args.trim_5p_seq
            else:
                Trim5_R2 = args.trim_5p_seq
                Trim5_R1 = args.trim_5p_seq_r2 if args.trim_5p_seq_r2 is not None else args.trim_5p_seq

            header = "=== Trimming (PE) ==="
            label = "transposon" if not args.primer_check else "primer"

            CurrRes, out, err, rc = RemoveTn_PE(
                args.r1, args.r2,
                Trim5_R1, Trim5_R2,
                overlap,
                CleanR1, CleanR2,
                CutAdaptPath,
                args.discard_untrimmed,
                trim_mate_mode=trim_mate_mode
            )
            LogFile.write(
                f"{header}\r\n{CurrRes[0]} reads; of these:\r\n"
                f"  {CurrRes[1]} ({CurrRes[2]}%) contained the {label}"
            )
            if rc != 0:
                LogFile.write("\r\nWARNING: cutadapt returned non-zero exit code.\r\n")
                LogFile.write(out.replace("\n", "\r\n") + "\r\n")
                LogFile.write(err.replace("\n", "\r\n") + "\r\n")

        if args.clean_adapters is not None:
            if args.clean_adapters == "":
                LogFile.write("\r\n\r\n=== Adapter removal skipped ===\r\n")
                LogFile.write("No adapter sequence provided with -a/--clean-adapters.\r\n")
            else:
                tmpR1 = CleanR1 + ".no_adapters.fq"
                tmpR2 = CleanR2 + ".no_adapters.fq"
                CurrRes2, out, err, rc = TrimAdapter3_PE(
                    CleanR1, CleanR2, args.clean_adapters, tmpR1, tmpR2, CutAdaptPath
                )
                LogFile.write(
                    "\r\n\r\n=== Adapter removal ===\r\n"
                    f"{CurrRes2[0]} reads; of these:\r\n"
                    f"  {CurrRes2[1]} ({CurrRes2[2]}%) contained the adapter"
                )
                os.replace(tmpR1, CleanR1)
                os.replace(tmpR2, CleanR2)

        # mapping PE with bowtie2

        SamfName = os.path.join(args.out_dir, f"{prefix}.sam")
        AlignFastq_PE(BowtiePath, args.bt2_index, CleanR1, CleanR2, SamfName, args.threads, LogFile)

        bam = os.path.splitext(SamfName)[0] + ".bam"
        sorted_bam = os.path.splitext(SamfName)[0] + ".sorted.bam"

        include_flag = ""
        if args.junction_mate == "R1":
            include_flag = " -f 0x40"
        elif args.junction_mate == "R2":
            include_flag = " -f 0x80"

        cmd_view = (
            f'samtools view -@ {args.threads} -b -q {args.mapq} -F {args.drop_flags}'
            f'{include_flag} "{SamfName}" > "{bam}"'
        )
        rc, out, err = run_cmd(cmd_view)
        if rc != 0:
            raise RuntimeError(f"samtools view failed (exit {rc}):\n{err}")

        rc, out, err = run_cmd(f'samtools sort -@ {args.threads} "{bam}" -o "{sorted_bam}"')
        if rc != 0:
            raise RuntimeError(f"samtools sort failed (exit {rc}):\n{err}")

        rc, out, err = run_cmd(f'samtools index -@ {args.threads} "{sorted_bam}"')
        if rc != 0:
            raise RuntimeError(f"samtools index failed (exit {rc}):\n{err}")

        if args.delete_originals:
            try:
                os.remove(args.r1)
            except OSError:
                pass
            try:
                os.remove(args.r2)
            except OSError:
                pass

        if not args.keep_clean_fastqs:
            try:
                os.remove(CleanR1)
            except OSError:
                pass
            try:
                os.remove(CleanR2)
            except OSError:
                pass

    else:
        CleanfName = os.path.join(args.out_dir, f"{prefix}.clean.fq")

        if args.trim_5p_seq is None:
            LogFile.write("=== Trimming step skipped ===\r\n")
            copy_to_clean_se(args.input_file_name, CleanfName)
        else:
            header = (
                "Reserve strand primer removal (searched w/o Tn tail)"
                if args.reverse_strand and args.primer_check else
                "=== Reverse strand transposon removal (searched with primer+Tn tail) ==="
                if args.reverse_strand else
                "=== Primer removal (searched w/o Tn tail) ==="
                if args.primer_check else
                "=== Transposon removal (searched with primer+Tn tail) ==="
            )
            label = "primer" if args.primer_check else "transposon"

            CurrRes, out, err, rc = RemoveTn_SE(
                args.input_file_name,
                args.trim_5p_seq,
                overlap,
                CleanfName,
                CutAdaptPath,
                args.discard_untrimmed
            )

            LogFile.write(
                f"{header}\r\n{CurrRes[0]} reads; of these:\r\n"
                f"  {CurrRes[1]} ({CurrRes[2]}%) contained the {label}"
            )
            if rc != 0:
                LogFile.write("\r\nWARNING: cutadapt returned non-zero exit code.\r\n")
                LogFile.write(out.replace("\n", "\r\n") + "\r\n")
                LogFile.write(err.replace("\n", "\r\n") + "\r\n")

        if args.clean_adapters is not None:
            if args.clean_adapters == "":
                LogFile.write("\r\n\r\n=== Adapter removal skipped ===\r\n")
                LogFile.write("No adapter sequence provided with -a/--clean-adapters.\r\n")
            else:
                temp_fq = CleanfName + ".no_adapters.fq"
                CurrRes2, out, err, rc = TrimAdapter3_SE(
                    CleanfName,
                    args.clean_adapters,
                    temp_fq,
                    CutAdaptPath
                )
                LogFile.write(
                    "\r\n\r\n=== Adapter removal ===\r\n"
                    f"{CurrRes2[0]} reads; of these:\r\n"
                    f"  {CurrRes2[1]} ({CurrRes2[2]}%) contained the adapter"
                )
                os.replace(temp_fq, CleanfName)

        # mapping SE with bowtie2

        SamfName = os.path.join(args.out_dir, f"{prefix}.sam")
        AlignFastq_SE(BowtiePath, args.bt2_index, CleanfName, SamfName, args.threads, LogFile)

        bam = os.path.splitext(SamfName)[0] + ".bam"
        sorted_bam = os.path.splitext(SamfName)[0] + ".sorted.bam"

        cmd_view = (
            f'samtools view -@ {args.threads} -b -q {args.mapq} -F {args.drop_flags} "{SamfName}" > "{bam}"'
        )
        rc, out, err = run_cmd(cmd_view)
        if rc != 0:
            raise RuntimeError(f"samtools view failed (exit {rc}):\n{err}")

        rc, out, err = run_cmd(f'samtools sort -@ {args.threads} "{bam}" -o "{sorted_bam}"')
        if rc != 0:
            raise RuntimeError(f"samtools sort failed (exit {rc}):\n{err}")

        rc, out, err = run_cmd(f'samtools index -@ {args.threads} "{sorted_bam}"')
        if rc != 0:
            raise RuntimeError(f"samtools index failed (exit {rc}):\n{err}")

        if args.delete_originals:
            try:
                os.remove(args.input_file_name)
            except OSError:
                pass

        if not args.keep_clean_fastqs:
            try:
                os.remove(CleanfName)
            except OSError:
                pass

    # final log output

    LogFile.seek(0)
    sys.stdout.write(LogFile.read())
    LogFile.close()

    try:
        os.remove(SamfName)
    except OSError:
        pass
    try:
        os.remove(bam)
    except OSError:
        pass


if __name__ == "__main__":
    main()
