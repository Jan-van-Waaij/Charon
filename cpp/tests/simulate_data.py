#!/usr/bin/env python3
"""Simulate CHARON input data from the model, for tests.

Each locus gets an anchor frequency y and a coverage R. The number of derived alleles k
among the 2n chromosomes is drawn from P(k | n, y, tauC, tauA) = binom(2n, k)
(e^{tauA Q} e^{tauC Qdown} h(y))_k (Schraiber 2018), and the number of derived reads
from Binomial(R, k(1-eps)/(2n) + (2n-k)eps/(2n)). Matrix exponentials use scipy's
expm, independently of the eigendecomposition in the C++ code.

Writes a DICE-2 file (aggregated) and/or a base count + frequency file pair.

    python3 simulate_data.py --n 3 --tauC 0.1 --tauA 0.2 --eps 0.01 --loci 5000 \
        --seed 1 --dice out.dice.csv --basecount out_basecount.txt.gz --freq out_freq.txt.gz
"""

import argparse
import collections
import gzip

import numpy as np
from scipy.linalg import expm
from scipy.special import comb


def tridiagonal(lower, diagonal, upper):
    return np.diag(diagonal) + np.diag(lower, -1) + np.diag(upper, 1)


def make_q(n):
    m = 2 * n
    return tridiagonal(0.5 * np.arange(1.0, m + 1) * np.arange(0.0, m),
                       -np.arange(0.0, m + 1) * np.arange(m, -1, -1.0),
                       0.5 * np.arange(m, 0, -1.0) * np.arange(m - 1, -1, -1.0))


def make_qdown(n):
    m = 2 * n
    return tridiagonal(0.5 * np.arange(1.0, m + 1) * np.arange(0.0, m),
                       -np.arange(0.0, m + 1) * np.arange(m + 1, 0, -1.0),
                       0.5 * np.arange(m + 1, 1, -1.0) * np.arange(m, 0, -1.0))


def open_text(path):
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--tauC", type=float, required=True)
    parser.add_argument("--tauA", type=float, required=True)
    parser.add_argument("--eps", type=float, required=True)
    parser.add_argument("--loci", type=int, required=True)
    parser.add_argument("--coverage", type=float, default=10.0, help="mean (Poisson) coverage")
    parser.add_argument("--freq-denominator", type=int, default=200, help="anchor frequencies are multiples of 1/D")
    parser.add_argument("--fixation-fraction", type=float, default=0.05,
                        help="fraction of loci with anchor frequency 0 or 1 (filtered by CHARON)")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--dice")
    parser.add_argument("--basecount")
    parser.add_argument("--freq")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    n = args.n
    M = expm(args.tauA * make_q(n)) @ expm(args.tauC * make_qdown(n))
    ks = np.arange(2 * n + 1)
    binom2n = comb(2 * n, ks)
    p_read = ks * (1 - args.eps) / (2 * n) + (2 * n - ks) * args.eps / (2 * n)

    D = args.freq_denominator
    rows = []
    for _ in range(args.loci):
        if rng.random() < args.fixation_fraction:
            y = float(rng.integers(0, 2))
        else:
            # U-shaped site frequency spectrum, rounded to multiples of 1/D.
            y = min(max(round(rng.beta(0.4, 0.6) * D), 1), D - 1) / D
        R = int(rng.poisson(args.coverage))
        h = y ** ks * (1 - y) ** (2 * n - ks)
        pk = np.clip(binom2n * (M @ h), 0, None)
        k = rng.choice(ks, p=pk / pk.sum())
        d = int(rng.binomial(R, p_read[k])) if R > 0 else 0
        rows.append((R, d, y))

    if args.dice:
        counts = collections.Counter(rows)
        with open_text(args.dice) as f:
            f.write("ancestral,derived,freq,counts\n")
            for (R, d, y), c in sorted(counts.items()):
                f.write(f"{R - d},{d},{y!r},{c}\n")
    if args.basecount:
        with open_text(args.basecount) as f:
            for R, d, _ in rows:
                f.write(f"{d} {R}\n")
    if args.freq:
        with open_text(args.freq) as f:
            for _, _, y in rows:
                f.write(f"{y!r}\n")


if __name__ == "__main__":
    main()
