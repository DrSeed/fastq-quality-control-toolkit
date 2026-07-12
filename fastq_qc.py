#!/usr/bin/env python3
# Pure-Python reimplementations of the classic FASTQ QC one-liners.
# Each function mirrors an awk/sed one-liner but is importable and testable.
import numpy as np


def parse_fastq(lines):
    # Yield (header, sequence, quality) from an iterable of FASTQ lines.
    # A FASTQ record is 4 lines: @header / sequence / + / quality.
    it = iter(lines)
    for header in it:
        seq = next(it).rstrip("\n")
        next(it)                       # the '+' separator line
        qual = next(it).rstrip("\n")
        yield header.rstrip("\n"), seq, qual


def read_length_distribution(records):
    # Equivalent to: awk 'NR%4==2 {lengths[length($0)]++}'
    counts = {}
    for _, seq, _ in records:
        counts[len(seq)] = counts.get(len(seq), 0) + 1
    return counts


def gc_content(seq):
    # Fraction of G/C bases in one read.
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GCgc")
    return gc / len(seq)


def phred_to_scores(qual):
    # Decode a Phred+33 quality string to integer quality scores.
    return np.array([ord(c) - 33 for c in qual], dtype=int)


def per_base_quality(records, max_len):
    # Mean quality at each read position (the FastQC per-base quality plot).
    total = np.zeros(max_len)
    count = np.zeros(max_len)
    for _, _, qual in records:
        q = phred_to_scores(qual)
        total[: len(q)] += q
        count[: len(q)] += 1
    with np.errstate(invalid="ignore"):
        return np.where(count > 0, total / count, np.nan)


def barcode_frequency(barcodes):
    # Equivalent to: sed -n '2~4p' | sort | uniq -c | sort -k1,1 -nr
    counts = {}
    for bc in barcodes:
        counts[bc] = counts.get(bc, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: kv[1], reverse=True))


if __name__ == "__main__":
    demo = ["@r1", "GGGGCCCC", "+", "IIIIFFFF", "@r2", "ATAT", "+", "IIII"]
    recs = list(parse_fastq(demo))
    print("lengths:", read_length_distribution(recs))
    print("gc r1:", gc_content(recs[0][1]))
