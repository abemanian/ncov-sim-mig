#!/usr/bin/env python3
"""Subsample sequences from a FASTA file, optionally stratified by time.

pairsnp is O(n^2) in the number of sequences, so a full simulated outbreak of
tens of thousands of tips dominates the runtime of a batch.

Uniform sampling keeps RR ratios roughly unbiased but allocates tips in
proportion to incidence, which leaves the early epidemic almost unsampled: in a
run of 46k hosts where only 500 infections occur before t=40, a uniform 8k
sample retains fewer than 90 of them. Every early-window RR estimate is then
built from a handful of sequences and comes out noisy or exactly zero.

Passing --meta switches to stratified sampling: hosts are binned by infection
time and the budget is shared evenly across bins, with any quota a sparse bin
cannot fill redistributed to the others. Sampling stays uniform *within* each
bin, so the deme composition of pairs inside any one time window is still
unbiased -- which is what RR is computed on. What changes is that early windows
get far more tips than their share of incidence would give them.
"""

import argparse
import csv
import random
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="Input FASTA")
    p.add_argument("--output", required=True, help="Output FASTA")
    p.add_argument("--n", type=int, required=True, help="Number of sequences to keep")
    p.add_argument("--seed", type=int, default=17, help="RNG seed")
    p.add_argument(
        "--meta",
        help="Optional meta.tsv with hosts.ID and inf.time; enables time stratification",
    )
    p.add_argument(
        "--bin-width",
        type=float,
        default=10.0,
        help="Width of each time stratum (default: 10)",
    )
    return p.parse_args()


def records(handle):
    """Yield (header, [lines]) pairs from a FASTA stream."""
    header = None
    body = []
    for line in handle:
        if line.startswith(">"):
            if header is not None:
                yield header, body
            header = line
            body = []
        else:
            body.append(line)
    if header is not None:
        yield header, body


def host_id(header):
    """Leading token of a FASTA header, minus '>'."""
    return header[1:].strip().split()[0]


def load_bins(meta_path, bin_width):
    """Map host id -> time-bin index."""
    bins = {}
    with open(meta_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                bins[row["hosts.ID"]] = int(float(row["inf.time"]) // bin_width)
            except (KeyError, ValueError):
                continue
    if not bins:
        raise SystemExit(f"subsample_fasta: no usable rows in {meta_path}")
    return bins


def allocate(counts, budget):
    """Split budget across bins as evenly as their populations allow.

    Bins too small to use a full share release the remainder to the others, and
    the pass repeats until nothing more can be redistributed.
    """
    quota = {}
    remaining = dict(counts)
    left = budget
    while remaining and left > 0:
        share = left // len(remaining)
        if share == 0:
            # Fewer slots left than bins: hand them out one at a time, largest
            # bins first, so the budget is spent exactly.
            for b in sorted(remaining, key=lambda k: -remaining[k])[:left]:
                quota[b] = quota.get(b, 0) + 1
            left = 0
            break
        filled = [b for b, c in remaining.items() if c <= share]
        if not filled:
            for b in remaining:
                quota[b] = quota.get(b, 0) + share
                left -= share
            # Any rounding remainder goes to the largest bins
            for b in sorted(remaining, key=lambda k: -remaining[k])[:left]:
                quota[b] = quota.get(b, 0) + 1
            left = 0
            break
        for b in filled:
            quota[b] = quota.get(b, 0) + remaining[b]
            left -= remaining[b]
            del remaining[b]
    return quota


def main():
    args = parse_args()
    rng = random.Random(args.seed)

    bins = load_bins(args.meta, args.bin_width) if args.meta else None

    if bins is None:
        # Uniform reservoir sampling over the whole file
        reservoir = []
        seen = 0
        with open(args.input) as fh:
            for rec in records(fh):
                seen += 1
                if len(reservoir) < args.n:
                    reservoir.append(rec)
                else:
                    j = rng.randrange(seen)
                    if j < args.n:
                        reservoir[j] = rec
        kept = reservoir
        note = "uniform"
    else:
        # Two passes: count per bin so the budget can be allocated, then
        # reservoir-sample each bin to its quota.
        counts = defaultdict(int)
        with open(args.input) as fh:
            for header, _ in records(fh):
                b = bins.get(host_id(header))
                if b is not None:
                    counts[b] += 1
        if not counts:
            raise SystemExit(
                "subsample_fasta: no FASTA headers matched hosts.ID in --meta"
            )

        quota = allocate(dict(counts), args.n)
        reservoirs = defaultdict(list)
        seen = defaultdict(int)
        with open(args.input) as fh:
            for header, body in records(fh):
                b = bins.get(host_id(header))
                if b is None or quota.get(b, 0) == 0:
                    continue
                seen[b] += 1
                res = reservoirs[b]
                if len(res) < quota[b]:
                    res.append((header, body))
                else:
                    j = rng.randrange(seen[b])
                    if j < quota[b]:
                        res[j] = (header, body)
        kept = [rec for b in sorted(reservoirs) for rec in reservoirs[b]]
        note = f"stratified over {len(counts)} bins of width {args.bin_width}"

    with open(args.output, "w") as out:
        for header, body in kept:
            out.write(header)
            out.writelines(body)

    print(
        f"subsample_fasta: kept {len(kept)} sequences ({note}, seed {args.seed})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
