#!/usr/bin/env python3
# Self-contained demo: generate a synthetic FASTQ file, then reproduce the
# classic FASTQ QC plots (read-length distribution, per-base quality, GC
# content, index-barcode frequency) using pure-Python versions of the
# well-known awk/sed one-liners.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fastq_qc import (parse_fastq, read_length_distribution, gc_content,
                      per_base_quality, barcode_frequency, phred_to_scores)

RNG = np.random.default_rng(42)

N_READS = 5000
LENGTH_CHOICES = np.array([36, 50, 75, 100])
LENGTH_PROBS = np.array([0.1, 0.2, 0.3, 0.4])
BASES = np.array(list("ACGT"))
# A small index whitelist with deliberately skewed usage, plus a couple of
# off-whitelist "error" barcodes so the frequency table looks real.
WHITELIST = ["ACGTACGT", "TTGGCCAA", "GATCGATC", "CCCCGGGG", "AATTAATT", "TGCATGCA"]
WHITELIST_PROBS = np.array([0.34, 0.24, 0.17, 0.12, 0.08, 0.05])


def random_seq(length, gc_bias):
    # gc_bias in [0,1] tilts base composition toward G/C.
    p_gc = 0.5 * gc_bias
    p_at = 0.5 * (1 - gc_bias)
    probs = np.array([p_at, p_gc, p_gc, p_at])  # A C G T
    probs = probs / probs.sum()
    return "".join(RNG.choice(BASES, size=length, p=probs))


def random_quality(length):
    # Quality starts high and degrades toward the 3' end (typical Illumina),
    # then encode as Phred+33 ASCII.
    pos = np.arange(length)
    mean_q = 38 - 12 * (pos / max(length - 1, 1)) ** 2
    q = np.clip(RNG.normal(mean_q, 2.5), 2, 40).astype(int)
    return "".join(chr(int(v) + 33) for v in q)


def write_synthetic_fastq(path):
    lengths = RNG.choice(LENGTH_CHOICES, size=N_READS, p=LENGTH_PROBS)
    gc_bias = RNG.beta(5, 5, size=N_READS)          # centred near 0.5
    bc_idx = RNG.choice(len(WHITELIST), size=N_READS, p=WHITELIST_PROBS)
    with open(path, "w") as fh:
        for i in range(N_READS):
            L = int(lengths[i])
            seq = random_seq(L, gc_bias[i])
            qual = random_quality(L)
            bc = WHITELIST[bc_idx[i]]
            if RNG.random() < 0.03:                 # 3% sequencing error in index
                j = RNG.integers(len(bc))
                bc = bc[:j] + RNG.choice(BASES) + bc[j + 1:]
            fh.write(f"@read{i}:index={bc}\n{seq}\n+\n{qual}\n")


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    os.makedirs("data", exist_ok=True)
    fq = "data/synthetic.fastq"
    write_synthetic_fastq(fq)

    # Read the records once, materialise, and reuse for each metric.
    with open(fq) as fh:
        records = list(parse_fastq(fh))

    len_dist = read_length_distribution(records)
    gc = np.array([gc_content(seq) for _, seq, _ in records])
    max_len = max(len(seq) for _, seq, _ in records)
    pbq = per_base_quality(records, max_len)
    barcodes = [h.split("index=")[1] for h, _, _ in records]
    bc_freq = barcode_frequency(barcodes)

    fig, ax = plt.subplots(2, 2, figsize=(12, 9))

    # Panel 1: read-length distribution (Day 1 one-liner).
    a = ax[0, 0]
    xs = sorted(len_dist)
    a.bar([str(x) for x in xs], [len_dist[x] for x in xs], color="#4C72B0")
    a.set_xlabel("read length (bp)"); a.set_ylabel("number of reads")
    a.set_title("Read-length distribution")

    # Panel 2: per-base mean quality across read position.
    a = ax[0, 1]
    a.plot(np.arange(1, max_len + 1), pbq, color="#55A868")
    a.axhline(30, color="#C44E52", ls="--", lw=1, label="Q30")
    a.axhline(20, color="#DD8452", ls="--", lw=1, label="Q20")
    a.set_xlabel("position in read (bp)"); a.set_ylabel("mean quality (Phred)")
    a.set_title("Per-base quality decays toward 3' end"); a.legend(fontsize=8)

    # Panel 3: GC-content distribution per read.
    a = ax[1, 0]
    a.hist(gc * 100, bins=30, color="#8172B3", edgecolor="white")
    a.axvline(gc.mean() * 100, color="#C44E52", lw=1.5,
              label=f"mean {gc.mean()*100:.1f}%")
    a.set_xlabel("GC content per read (%)"); a.set_ylabel("number of reads")
    a.set_title("GC-content distribution"); a.legend(fontsize=8)

    # Panel 4: top index-barcode frequency (Day 20 one-liner).
    a = ax[1, 1]
    top = list(bc_freq.items())[:8]
    a.barh([b for b, _ in top][::-1], [c for _, c in top][::-1], color="#4C72B0")
    a.set_xlabel("read count"); a.set_ylabel("index barcode")
    a.set_title("Top index barcodes (whitelist stands out)")

    fig.suptitle("FASTQ quality-control toolkit (synthetic reads)", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig("figures/demo.png", dpi=120)

    all_q = np.concatenate([phred_to_scores(q) for _, _, q in records])
    summary = pd.DataFrame([{
        "n_reads": len(records),
        "mean_read_length": np.mean([len(s) for _, s, _ in records]),
        "mean_gc_percent": gc.mean() * 100,
        "mean_quality": all_q.mean(),
        "percent_bases_q30": np.mean(all_q >= 30) * 100,
        "distinct_barcodes": len(bc_freq),
        "top_barcode": next(iter(bc_freq)),
        "top_barcode_reads": next(iter(bc_freq.values())),
    }])
    summary.to_csv("results/summary.csv", index=False)

    print(f"Wrote {N_READS} reads to {fq}")
    print(f"Length distribution: {dict(sorted(len_dist.items()))}")
    print(f"Mean GC: {gc.mean()*100:.1f}%  Mean Q: {all_q.mean():.1f}  %Q30: {np.mean(all_q>=30)*100:.1f}")
    print(f"Distinct barcodes: {len(bc_freq)}, top = {next(iter(bc_freq))}")
    print("Wrote figures/demo.png and results/summary.csv")


if __name__ == "__main__":
    main()
