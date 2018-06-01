#!/usr/bin/env python
import sys
import numpy as np


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--npairs", type=int, default=20,
                   help=("Number of pairs within bucket to store (to "
                         "account for potential inaccuracies"
                         " in our estimation."))
    p.add_argument("--linspace", "-l",
                   type=int, default=100,
                   help="Number of buckets between 0 and 1 to cover.")
    p.add_argument("-o", "--outpath", default="/dev/stdout")
    p.add_argument("distfile", help="Path to distance file")
    a = p.parse_args()
    p = a.npairs
    labels = None
    lsp = a.linspace
    idxpairs = [[] for i in range(lsp)]
    print("Len of idx pairs at start: %i" % len(idxpairs), file=sys.stderr)
    assert len(idxpairs) == lsp
    linenum = 0

    def unfilled_bucket_count():
        return sum(map(lambda x: len(x) < p, idxpairs))
    names = next(open(a.distfile)).strip().split('\t')[1:]
    for line in open(a.distfile):
        if line[0] == '#':
            if not labels:
                labels = line[1:].strip().split()
            continue
        gen = (i for i in line.strip().split('\t')[1:] if i != '-')
        for ind, dist in enumerate(map(float, gen), 1):
            bucket = int(lsp * dist)
            if bucket < lsp:
                dest = idxpairs[bucket]
                if len(dest) < p:
                    dest.append((linenum, ind + linenum))
        linenum += 1
        if linenum % 100 == 0:
            ubc = unfilled_bucket_count()
            if ubc == 0:
                break
            print(f"Unfilled bucket count: {unfilled_bucket_count()}. "
                  f"Total buckets to fill: {len(idxpairs)}",
                  file=sys.stderr)
    if unfilled_bucket_count():
        print(f"Unfilled bucket count: {unfilled_bucket_count()}. "
              "Emitting unfinished results.", file=sys.stderr)
    else:
        print("No unfilled bucket count any more, we got 'em all.",
              file=sys.stderr)
    with open(a.outpath, "w") as f:
        for pairlist in idxpairs:
            for x, y in pairlist:
                f.write("%s\t%s\n" % (names[x], names[y]))
        f.write("=" * 80)
        f.write("#Fraction\tIndices of pairs [0-based]...\n")
        for values, frac in zip(idxpairs, np.linspace(0., 1., lsp + 1)[:-1]):
            f.write("%f\t%s\n" % (frac, '\t'.join(map(str, values))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
