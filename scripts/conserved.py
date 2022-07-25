#! /usr/bin/env python3
'''
Compute max base frequency and Shannon Entropy at each position of a Multiple Sequence Alignment (one-line FASTA). Input from STDIN and output to STDOUT.
'''
from math import log2
from sys import stdin, stderr
MAX_ENTROPY = 2
NUM_SEQS_PROGRESS = 1000
NUCS_SET = {'A','C','G','T'}

if __name__ == "__main__":
    # count bases at each position
    counts = None # counts[pos][nuc] = count
    num_seqs = 0
    for line in stdin:
        if line.startswith('>'):
            continue
        num_seqs += 1; l = line.strip().upper()
        if counts is None:
            counts = [dict() for _ in range(len(l))]
        for i, c in enumerate(l):
            if c not in NUCS_SET:
                continue
            if c not in counts[i]:
                counts[i][c] = 0
            counts[i][c] += 1
        if num_seqs % NUM_SEQS_PROGRESS == 0:
            stderr.write("Parsed %d sequences...\n" % num_seqs); stderr.flush()

    # compute shannon entropy and print
    print("Position (1-based)\tMax Base Frequency\tShannon Entropy")
    for i, curr in enumerate(counts):
        ent = MAX_ENTROPY
        if len(curr) != 0:
            ent = 0; tot = sum(curr.values()); m = max(curr.values())/tot
            for c in curr:
                p = curr[c]/tot
                ent -= (p*log2(p))
        print("%d\t%s\t%s" % (i+1,m,ent))
