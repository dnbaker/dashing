#!/usr/bin/env python
import sys
import numpy as np

def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--npairs", type=int, default=20, help="Number of pairs within bucket to store (to account for potential inaccuracies in our estimation.")
    p.add_argument("--linspace", "-l", type=int, default=20, help="Number of buckets between 0 and 1 to cover.")
    p.add_argument("distfile", help="Path to distance file")
    a = p.parse_args()
    p = a.npairs
    labels = None
    spacing = np.linspace(0., 1., p.linspace + 1)
    idxpairs = {s:[] for s in np.linspace(0., 1., p.linspace + 1)}
    for line in open(a.distfile):
        if line[0] == '#':
            if not labels:
                labels = line[1:].strip().split()
            continue
        dists = list(map(float, line.strip().split('\t')[1:]))
    
