#!/usr/bin/env python3

import os
import sys
import argparse
from subprocess import PIPE, Popen

# -----------------------------------------------------------------------------
# MapFastqAll.py
# -----------------------------------------------------------------------------
# Universal FASTQ -> sorted/indexed BAM mapping for transposon insertion
# sequencing.
#
# Adapted from Berman-lab transposon-pipeline MapFastq.py.  When adapting, we
# carried over the original inline comments where they explain why a step exists.
#
# In our 2026 C. glabrata runs, reads start directly in genomic DNA (the Tn end
# is not present in the sequencing read), and primers/adapters were spiked so that
# many reads do not contain them. Therefore, trimming/filtering is OFF by default.
# -----------------------------------------------------------------------------

# ----------------------------- defaults from original -------------------------
# These are kept as defaults for optional trimming modes.
TnPrimerAndTail = 'GTATTTTACCGACCGTTACCGACCGTTTTCATCCCTA'
TnRev = 'TAGGGATGAAAACGGTCGGTAACGGTCGGTAAAATAC'
PrimerOnly = 'GTATTTTACCGACCGTTACCGACC'
PrimerRev = 'GGTCGGTAACGGTCGGTAAAATAC'
AdapterSeq = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'


def die(msg: str, code: int = 1):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def GetCmdPath(program: str) -> str:
    """Gets the path to a desired program file on given computer.

    Parameters
    ----------
    program : str
        Name of program wish to have path for.

    Returns
    -------
    str
        Path for calling program as a command in the POSIX terminal.

    Notes
    -----
    In the original pipeline, the script relied on `command -v`.
    """
    shell_path = os.popen('command -v ' + program).read()
    cmd_path = shell_path.strip('\n')
    return cmd_path


def cmdline(command: str):
    """Takes text of command and has shell process it.

    Parameters
    ----------
    command : str
        Text of command that needs to be processed by the POSIX shell.

    Returns
    -------
    (stdout, stderr) bytes

    Notes
    -----
    The original MapFastq.py used Popen(shell=True) and then process.communicate().
    We keep that behavior so command strings look familiar and log output is similar.
    """
    process = Popen(args=command, stdout=PIPE, stderr=PIPE, shell=True)
    return process.communicate()


def GetReads(val: str, text: str) -> int:
    """Gets number of reads found in given subcategory of cutadapt output."""
    sub = text[text.find(val) + len(val) + 1:].lstrip(' ')
    sub = sub[:sub.find('\n')]
    vals = sub.split(' ')
    return int(vals[0].replace(',', ''))


def CutAdaptOutput(output: str, paired: bool) -> tuple[int, int, float]:
    """Gets desired numerical values from CutAdapt program output.

    Returns
    -------
    TotalRead : int
        Total number of reads (or read pairs) processed.
    TotalWrite : int
        Total number of reads (or read pairs) written.
    Percent : float
        Percentage written.

    Notes
    -----
    Original MapFastq.py assumed PE and parsed "Total read pairs processed".
    Here we support SE and PE.
    """
    if '=== Summary ===' not in output:
        return 0, 0, 0.0

    summary = output[output.find('=== Summary ===') + 17:]
    if paired:
        total_read = GetReads('Total read pairs processed', summary)
    else:
        total_read = GetReads('Total reads processed', summary)

    total_write = GetReads('Reads written (passing filters)', summary)
    percent = round((float(total_write) / float(total_read) * 100) if total_read else 0.0, 2)
    return total_read, total_write, percent


def RemoveTn_SE(cutadapt_path: str, input_fq: str, tn_seq: str, overlap: int, clean_fq: str) -> tuple[int, int, float]:
    """Locates reads with transposon/primer and removes from 5' end using cutadapt.

    Parameters - RemoveTn
    ---------------------
    input_fq : str
        Input FASTQ(.gz) filename.
    tn_seq : str
        Sequence to find and remove.
    overlap : int
        Minimal number of bases from tn_seq required.
    clean_fq : str
        Output filename for temporary fq file created by program.

    Parameters - CutAdapt
    ---------------------
    -g : from the beginning
    -o : output
    --discard-untrimmed : keep only reads that contain the sequence
    --overlap : minimal overlap with tn_seq

    Returns
    -------
    (total_processed, total_written, percent_written)
    """
    cutadapt_cmd = (
        f'{cutadapt_path} -m 2 -g {tn_seq} -o "{clean_fq}" "{input_fq}" '
        f'--discard-untrimmed --overlap {overlap}'
    )
    out, err = cmdline(cutadapt_cmd)
    out_s = out.decode('utf-8', errors='replace')
    err_s = err.decode('utf-8', errors='replace')
    # cutadapt writes summary to stdout; if not, fallback to stderr
    summary_src = out_s if '=== Summary ===' in out_s else err_s
    return CutAdaptOutput(summary_src, paired=False)


def RemoveTn_PE(cutadapt_path: str, fq1: str, fq2: str, tn_seq: str, overlap: int, clean1: str, clean2: str) -> tuple[int, int, float]:
    """PE version of RemoveTn using cutadapt.

    Notes
    -----
    We trim/filter both mates in sync and keep only read pairs passing filters.
    """
    cutadapt_cmd = (
        f'{cutadapt_path} -m 2 -g {tn_seq} -o "{clean1}" -p "{clean2}" "{fq1}" "{fq2}" '
        f'--discard-untrimmed --overlap {overlap}'
    )
    out, err = cmdline(cutadapt_cmd)
    out_s = out.decode('utf-8', errors='replace')
    err_s = err.decode('utf-8', errors='replace')
    summary_src = out_s if '=== Summary ===' in out_s else err_s
    return CutAdaptOutput(summary_src, paired=True)


def RemoveAdapters_SE(cutadapt_path: str, input_fq: str, adapter_seq: str, output_fq: str) -> tuple[int, int, float]:
    """Remove 3' adapters from SE reads using cutadapt."""
    out, err = cmdline(f'{cutadapt_path} -m 2 -a {adapter_seq} -o "{output_fq}" "{input_fq}"')
    out_s = out.decode('utf-8', errors='replace')
    err_s = err.decode('utf-8', errors='replace')
    summary_src = out_s if '=== Summary ===' in out_s else err_s
    return CutAdaptOutput(summary_src, paired=False)


def RemoveAdapters_PE(cutadapt_path: str, fq1: str, fq2: str, adapter_seq: str, out1: str, out2: str) -> tuple[int, int, float]:
    """Remove 3' adapters from PE reads using cutadapt."""
    out, err = cmdline(
        f'{cutadapt_path} -m 2 -a {adapter_seq} -A {adapter_seq} -o "{out1}" -p "{out2}" "{fq1}" "{fq2}"'
    )
    out_s = out.decode('utf-8', errors='replace')
    err_s = err.decode('utf-8', errors='replace')
    summary_src = out_s if '=== Summary ===' in out_s else err_s
    return CutAdaptOutput(summary_src, paired=True)


def AlignFastq_SE(bowtie2_path: str, index_prefix: str, fq: str, sam: str, threads: int, extra: str, log_handle):
    """Run Bowtie2 alignment of SE reads to a reference genome index.

    Parameters - Bowtie
    -------------------
    -x : reference genome index prefix
    -U : fq file of unpaired reads
    -S : output SAM alignment file

    Writes
    ------
    Sequence alignment summary into logfile.
    """
    # Original script wrote Bowtie2 stderr to a temp file then appended to the log.
    bowtie_cmd = f'{bowtie2_path} -p {threads} -x "{index_prefix}" -U "{fq}" -S "{sam}" {extra} 2> temp.bowtie2'
    cmdline(bowtie_cmd)

    with open('temp.bowtie2', 'r') as tmp:
        temp_txt = tmp.read()

    log_handle.write('\r\n\r\n=== Sequence alignment ===\r\n' + temp_txt.replace('\n', '\r\n'))
    os.remove('temp.bowtie2')


def AlignFastq_PE(bowtie2_path: str, index_prefix: str, fq1: str, fq2: str, sam: str, threads: int, extra: str, log_handle):
    """Run Bowtie2 alignment of PE reads to a reference genome index."""
    bowtie_cmd = (
        f'{bowtie2_path} -p {threads} -x "{index_prefix}" -1 "{fq1}" -2 "{fq2}" -S "{sam}" {extra} 2> temp.bowtie2'
    )
    cmdline(bowtie_cmd)

    with open('temp.bowtie2', 'r') as tmp:
        temp_txt = tmp.read()

    log_handle.write('\r\n\r\n=== Sequence alignment ===\r\n' + temp_txt.replace('\n', '\r\n'))
    os.remove('temp.bowtie2')


def sam_to_sorted_bam(samtools_path: str, sam: str, sorted_bam: str, threads: int, mapq: int, drop_flags: str, log_handle):
    """Convert SAM -> BAM, filter, sort, index.

    Notes
    -----
    - We keep only alignments with MAPQ >= mapq.
    - We drop unmapped/secondary/duplicates/supplementary by default.

    The default DROP_FLAGS matches our current mapping script:
      0x4 (unmapped) + 0x100 (secondary) + 0x400 (duplicate) + 0x800 (supplementary)
      = 0xD04
    """
    bam = os.path.splitext(sorted_bam)[0] + '.unsorted.bam'

    # SAM -> QC BAM
    cmdline(
        f'{samtools_path} view -@ {threads} -b -q {mapq} -F {drop_flags} "{sam}" > "{bam}"'
    )

    # Sorting and indexing
    cmdline(f'{samtools_path} sort -@ {threads} "{bam}" -o "{sorted_bam}"')
    cmdline(f'{samtools_path} index -@ {threads} "{sorted_bam}"')

    log_handle.write(
        f'\r\n\r\n=== BAM processing ===\r\n'
        f'  MAPQ cutoff: {mapq}\r\n'
        f'  Dropped flags (-F): {drop_flags}\r\n'
        f'  Output: {sorted_bam}\r\n'
    )

    # remove intermediate BAM
    try:
        os.remove(bam)
    except OSError:
        pass


def main():
    usage = (
        "MapFastqAll.py\n"
        "  --index PREFIX            REQUIRED. Bowtie2 index prefix (the -x argument).\n"
        "  --fastq1 FILE             REQUIRED. Read1 FASTQ(.gz) (SE uses fastq1 only).\n"
        "  --fastq2 FILE             Optional. Read2 FASTQ(.gz) -> enables PE mode.\n"
        "  --out-dir DIR             Output directory (default: .).\n"
        "  --prefix NAME             Output prefix (default derived from fastq1).\n"
        "  --threads INT             Threads for bowtie2/samtools (default: 4).\n"
        "  --mapq INT                MAPQ cutoff for BAM (default: 20).\n"
        "  --drop-flags HEX          samtools -F flags to drop (default: 0xD04).\n"
        "\n"
        "Optional trimming/filtering (default OFF):\n"
        "  --trim-mode {none,tn,primer,tn_rev,primer_rev}  Sequence to filter/trim from 5'.\n"
        "  --overlap INT             Overlap for trim-mode (default: 37 for tn, 24 for primer).\n"
        "  --clean-adapters          Remove Illumina universal adapter from 3'.\n"
        "\n"
        "Housekeeping:\n"
        "  --delete-originals        Delete input FASTQ files.\n"
        "  --keep-clean-fastqs       Keep intermediate cleaned FASTQ(s).\n"
        "  --keep-sam                Keep SAM file (default removes).\n"
        "\n"
        "Advanced:\n"
        "  --bowtie2-extra \"...\"       Extra args appended to bowtie2 command (e.g. --very-sensitive).\n"
    )

    parser = argparse.ArgumentParser(add_help=False, usage=usage)
    parser.add_argument('--index', required=True)
    parser.add_argument('--fastq1', required=True)
    parser.add_argument('--fastq2', default=None)
    parser.add_argument('-o', '--out-dir', default='.')
    parser.add_argument('--prefix', default=None)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--mapq', type=int, default=20)
    parser.add_argument('--drop-flags', default='0xD04')

    parser.add_argument('--trim-mode', choices=['none', 'tn', 'primer', 'tn_rev', 'primer_rev'], default='none')
    parser.add_argument('--overlap', type=int, default=None)
    parser.add_argument('--clean-adapters', action='store_true', default=False)

    parser.add_argument('--delete-originals', action='store_true', default=False)
    parser.add_argument('--keep-clean-fastqs', action='store_true', default=False)
    parser.add_argument('--keep-sam', action='store_true', default=False)

    parser.add_argument('--bowtie2-extra', default='')
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS)

    args = parser.parse_args()

    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    fq1 = args.fastq1
    fq2 = args.fastq2
    paired = fq2 is not None

    if not os.path.exists(fq1):
        die(f"FASTQ not found: {fq1}")
    if paired and (not os.path.exists(fq2)):
        die(f"FASTQ2 not found: {fq2}")

    # Derive prefix if not provided
    if args.prefix:
        prefix = args.prefix
    else:
        base = os.path.basename(fq1)
        for suf in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
            if base.endswith(suf):
                base = base[:-len(suf)]
                break
        prefix = base

    log_path = os.path.join(out_dir, f"{prefix}_log.txt")
    log_handle = open(log_path, 'w+', newline='')

    cutadapt_path = GetCmdPath('cutadapt')
    bowtie2_path = GetCmdPath('bowtie2')
    samtools_path = GetCmdPath('samtools')

    if not bowtie2_path:
        die("bowtie2 not found in PATH")
    if not samtools_path:
        die("samtools not found in PATH")
    if (args.trim_mode != 'none' or args.clean_adapters) and (not cutadapt_path):
        die("cutadapt requested but not found in PATH")

    # ----------------- optional trimming/filtering -----------------
    # First removing the transposon head from the beginning and then removing the adapter.
    # This is because the sequencing tech, for whatever reason, removes the adapater from
    # the 5' Tn end, but keeps the tail adapters (sometimes), presumably if the reads are
    # too short. After we have all of the reads that have a transposon, we then trim the
    # adapter from the tail.
    #
    # In our pipeline, we default to doing NONE of this (trim_mode=none, clean_adapters=off).

    clean_fq1 = os.path.join(out_dir, f"{prefix}.clean.R1.fq")
    clean_fq2 = os.path.join(out_dir, f"{prefix}.clean.R2.fq") if paired else None

    # Start with inputs; if we trim, update curr_fq1/curr_fq2.
    curr_fq1 = fq1
    curr_fq2 = fq2

    if args.trim_mode != 'none':
        # Pick sequence + default overlap
        if args.trim_mode == 'tn':
            seq = TnPrimerAndTail
            default_ov = 37
            label = 'Transposon removal (searched with primer+Tn tail)'
        elif args.trim_mode == 'primer':
            seq = PrimerOnly
            default_ov = 24
            label = 'Primer removal (searched w/o Tn tail)'
        elif args.trim_mode == 'tn_rev':
            seq = TnRev
            default_ov = 37
            label = 'Reverse strand transposon removal (searched with primer+Tn tail)'
        else:  # primer_rev
            seq = PrimerRev
            default_ov = 24
            label = 'Reverse strand primer removal (searched w/o Tn tail)'

        overlap = args.overlap if args.overlap is not None else default_ov

        if paired:
            res = RemoveTn_PE(cutadapt_path, curr_fq1, curr_fq2, seq, overlap, clean_fq1, clean_fq2)
        else:
            res = RemoveTn_SE(cutadapt_path, curr_fq1, seq, overlap, clean_fq1)

        log_handle.write(
            f"=== {label} ===\r\n"
            f"{res[0]} reads{' pairs' if paired else ''}; of these:\r\n"
            f"  {res[1]} ({res[2]}%) contained the sequence\r\n"
        )

        curr_fq1 = clean_fq1
        curr_fq2 = clean_fq2

    if args.clean_adapters:
        # Remove adapters from the 3' end. For PE, do both mates.
        temp1 = os.path.join(out_dir, f"{prefix}.no_adapters.R1.fq")
        temp2 = os.path.join(out_dir, f"{prefix}.no_adapters.R2.fq") if paired else None

        if paired:
            res = RemoveAdapters_PE(cutadapt_path, curr_fq1, curr_fq2, AdapterSeq, temp1, temp2)
        else:
            res = RemoveAdapters_SE(cutadapt_path, curr_fq1, AdapterSeq, temp1)

        log_handle.write(
            f"\r\n\r\n=== Adapter removal ===\r\n"
            f"{res[0]} reads{' pairs' if paired else ''}; of these:\r\n"
            f"  {res[1]} ({res[2]}%) contained the adapter\r\n"
        )

        # Swap temp -> current
        if curr_fq1 != fq1:
            try:
                os.remove(curr_fq1)
            except OSError:
                pass
        if paired and curr_fq2 and curr_fq2 != fq2:
            try:
                os.remove(curr_fq2)
            except OSError:
                pass

        os.rename(temp1, clean_fq1)
        curr_fq1 = clean_fq1

        if paired:
            os.rename(temp2, clean_fq2)
            curr_fq2 = clean_fq2

    # ----------------- alignment -----------------
    sam_path = os.path.join(out_dir, f"{prefix}.sam")

    if paired:
        AlignFastq_PE(bowtie2_path, args.index, curr_fq1, curr_fq2, sam_path, args.threads, args.bowtie2_extra, log_handle)
    else:
        AlignFastq_SE(bowtie2_path, args.index, curr_fq1, sam_path, args.threads, args.bowtie2_extra, log_handle)

    # ----------------- SAM -> sorted BAM -----------------
    sorted_bam = os.path.join(out_dir, f"{prefix}.sorted.bam")
    sam_to_sorted_bam(
        samtools_path=samtools_path,
        sam=sam_path,
        sorted_bam=sorted_bam,
        threads=args.threads,
        mapq=args.mapq,
        drop_flags=args.drop_flags,
        log_handle=log_handle,
    )

    # Print log to stdout (like original script), then close.
    log_handle.flush()
    log_handle.seek(0)
    sys.stdout.write(log_handle.read())
    log_handle.close()

    # ----------------- cleanup -----------------
    # remove intermediates:
    if args.delete_originals:
        try:
            os.remove(fq1)
        except OSError:
            pass
        if paired:
            try:
                os.remove(fq2)
            except OSError:
                pass

    # Keep or delete cleaned FASTQs
    if not args.keep_clean_fastqs:
        # Only delete if they were created
        if curr_fq1 != fq1 and os.path.exists(curr_fq1):
            try:
                os.remove(curr_fq1)
            except OSError:
                pass
        if paired and curr_fq2 and curr_fq2 != fq2 and os.path.exists(curr_fq2):
            try:
                os.remove(curr_fq2)
            except OSError:
                pass

    if not args.keep_sam:
        try:
            os.remove(sam_path)
        except OSError:
            pass


if __name__ == '__main__':
    main()