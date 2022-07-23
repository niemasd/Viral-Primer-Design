#! /usr/bin/env python3
'''
Compute Shannon Entropy at each position of a Multiple Sequence Alignment (one-line FASTA). Input from STDIN and output to STDOUT.
'''
from math import log2
from sys import stdin
MAX_ENTROPY = 2

if __name__ == "__main__":
    # count bases at each position
    counts = None # counts[pos][nuc] = count
    for line in stdin:
        if line.startswith('>'):
            continue
        l = line.strip().upper()
        if counts is None:
            counts = [dict() for _ in range(len(l))]
        for i, c in enumerate(l):
            if c not in counts[i]:
                counts[i][c] = 0
            counts[i][c] += 1

    # compute shannon entropy and print
    print("Position (1-based)\tShannon Entropy")
    for i, curr in enumerate(counts):
        ent = MAX_ENTROPY
        if len(curr) != 0:
            ent = 0; tot = sum(curr.values())
            for c in curr:
                p = curr[c]/tot
                ent -= (p*log2(p))
        print("%d\t%s" % (i+1,ent))
