#! /usr/bin/env python3
'''
Pipeline for designing primers from a collection of viral sequences
'''
VERSION = '0.0.1'

# imports
from datetime import datetime
from gzip import open as gopen
from math import log2
from os import mkdir
from os.path import abspath, expanduser, getsize, isdir, isfile
from seaborn import relplot
from subprocess import call, check_output, STDOUT
from sys import argv, stdin, stderr
import argparse
import matplotlib.pyplot as plt

# constants
MIN_MAFFT_OUTPUT_SIZE = 100
MIN_MAFFT_VERSION = 7.467
MIN_PRIMER3_OUTPUT_SIZE = 100
NUC_COLORS = {'A':'red', 'C':'blue', 'G':'purple', 'T':'yellow', '-':'black'}
NUCS_SORT = ['A', 'C', 'G', 'T', '-']; NUCS = set(NUCS_SORT)
NUM_SEQS_PROGRESS = 500

# defaults
DEFAULT_BLAST_WORD_SIZE = 11
DEFAULT_BUFSIZE = 1048576 # 1 MB
DEFAULT_PRIMER3_PRIMER_MAX_SIZE = 36
DEFAULT_PRIMER3_PRIMER_MIN_SIZE = 18
DEFAULT_PRIMER3_PRIMER_OPT_SIZE = 20
DEFAULT_PRIMER3_PRIMER_PRODUCT_MIN_SIZE = 50
DEFAULT_PRIMER3_PRIMER_PRODUCT_MAX_SIZE = 170
DEFAULT_PRIMER3_STEP_SIZE = 500
DEFAULT_PRIMER3_WINDOW_SIZE = 1000

# return the current time as a string
def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# helper class for logging
class Log:
    def __init__(self, loggern, quiet=False, bufsize=DEFAULT_BUFSIZE):
        self.log_f = open(loggern, 'w', buffering=bufsize)
        if quiet:
            self.ostream = None
        else:
            self.ostream = stderr
    def __del__(self):
        self.log_f.close()
    def write(self, s, include_time=True):
        if include_time:
            s = '[%s] %s' % (get_time(), s)
        self.log_f.write(s)
        if self.ostream is not None:
            self.ostream.write(s)
        self.flush()
    def flush(self):
        self.log_f.flush()
        if self.ostream is not None:
            self.ostream.flush()

# parse user args
def parse_args():
    # check for -v/--version
    if '-v' in argv or '--version' in argv:
        print("ViralPrimerDesign v%s" % VERSION); exit(0)

    # use argparse to parse user arguments
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--sequences', required=False, type=str, default='stdin', help="Input Sequences (FASTA format)")
    parser.add_argument('-r', '--reference', required=True, type=str, help="Reference Genome (FASTA format)")
    parser.add_argument('-o', '--outdir', required=True, type=str, help="Output Directory")
    parser.add_argument('--skip_alignment', action="store_true", help="Skip Alignment (if input FASTA is already aligned)")
    parser.add_argument('-t', '--threads', required=False, type=int, default=1, help="Number of Threads (for MSA)")
    parser.add_argument('--primer3_window_size', required=False, type=int, default=DEFAULT_PRIMER3_WINDOW_SIZE, help="Primer3 Sliding Window Size")
    parser.add_argument('--primer3_step_size', required=False, type=int, default=DEFAULT_PRIMER3_STEP_SIZE, help="Primer3 Sliding Window Step Size")
    parser.add_argument('--primer3_primer_opt_size', required=False, type=int, default=DEFAULT_PRIMER3_PRIMER_OPT_SIZE, help="Primer3 Optimal Primer Length (PRIMER_OPT_SIZE)")
    parser.add_argument('--primer3_primer_min_size', required=False, type=int, default=DEFAULT_PRIMER3_PRIMER_MIN_SIZE, help="Primer3 Minimum Primer Length (PRIMER_MIN_SIZE)")
    parser.add_argument('--primer3_primer_max_size', required=False, type=int, default=DEFAULT_PRIMER3_PRIMER_MAX_SIZE, help="Primer3 Maximum Primer Length (PRIMER_MAX_SIZE)")
    parser.add_argument('--primer3_primer_product_min_size', required=False, type=int, default=DEFAULT_PRIMER3_PRIMER_PRODUCT_MIN_SIZE, help="Primer3 Minimum Primer Product Length (left part of PRIMER_PRODUCT_SIZE_RANGE)")
    parser.add_argument('--primer3_primer_product_max_size', required=False, type=int, default=DEFAULT_PRIMER3_PRIMER_PRODUCT_MAX_SIZE, help="Primer3 Maximum Primer Product Length (right part of PRIMER_PRODUCT_SIZE_RANGE)")
    parser.add_argument('--blast_word_size', required=False, type=int, default=DEFAULT_BLAST_WORD_SIZE, help="BLAST Word Size")
    parser.add_argument('-q', '--quiet', action="store_true", help="Quiet (hide verbose messages)")
    parser.add_argument('-v', '--version', action="store_true", help="Show Version Number")
    args = parser.parse_args()

    # check user arguments for validity
    if args.sequences != 'stdin' and not args.sequences.startswith('/dev/fd/') and not isfile(args.sequences):
        raise ValueError("File not found: %s" % args.sequences)
    if not args.reference.startswith('/dev/fd/') and not isfile(args.reference):
        raise ValueError("File not found: %s" % args.reference)
    args.outdir = abspath(expanduser(args.outdir))
    if isdir(args.outdir) or isfile(args.outdir):
        raise ValueError("Output exists: %s" % args.outdir)
    if args.threads < 1:
        raise ValueError("Number of threads must be positive integer: %s" % args.threads)
    if args.primer3_window_size < 1:
        raise ValueError("Window size must be positive integer: %s" % args.primer3_window_size)
    if args.primer3_step_size > args.primer3_window_size:
        raise ValueError("Window step size must be <= window size: %s" % args.primer3_step_size)
    if args.primer3_primer_opt_size < 1:
        raise ValueError("Optimal primer length must be positive integer: %s" % args.primer3_primer_opt_size)
    if args.primer3_primer_min_size > args.primer3_primer_opt_size:
        raise ValueError("Minimum primer length must be at most optimal primer length: %s" % args.primer3_primer_min_size)
    if args.primer3_primer_max_size < args.primer3_primer_opt_size:
        raise ValueError("Maximum primer length must be at least optimal primer length: %s" % args.primer3_primer_max_size)
    if args.primer3_primer_product_min_size <= args.primer3_primer_max_size:
        raise ValueError("Minimum primer product length must be greater than maximum primer length: %s" % args.primer3_primer_product_min_size)
    if args.primer3_primer_product_max_size < args.primer3_primer_product_min_size:
        raise ValueError("Maximum primer product length must be at least minimum primer product length: %s" % args.primer3_primer_product_max_size)
    if args.blast_word_size < 4:
        raise ValueError("BLAST word size must be at least 4: %s" % args.blast_word_size)
    return args

# align using MAFFT
def align_mafft(in_fn, ref_fn, out_fn, threads=1, logger=None, bufsize=DEFAULT_BUFSIZE):
    if in_fn.lower().endswith('.gz'):
        raise ValueError("MAFFT doesn't support gzipped input sequences: %s" % in_fn)
    try:
        mafft_version = check_output(['mafft', '--version'], stderr=STDOUT).decode()
    except:
        raise RuntimeError("Unable to execute 'mafft'. Are you sure it's in your PATH?")
    mafft_version = float(mafft_version.split()[0].lstrip('v'))
    if mafft_version < MIN_MAFFT_VERSION:
        raise RuntimeError("Must use MAFFT v%s or higher, but detected v%s" % (MIN_MAFFT_VERSION, mafft_version))
    for fn in [in_fn, ref_fn]:
        if not fn.startswith('/dev/fd/') and not isfile(fn):
            raise ValueError("File not found: %s" % fn)
    if isfile(out_fn) or isdir(out_fn):
        raise ValueError("Output exists: %s" % fn)
    command = ['mafft', '--thread', str(threads), '--6merpair', '--addfragments']
    if logger is not None:
        logger.write("Aligning sequences from: %s\n" % in_fn)
        logger.write("Reference genome: %s\n" % ref_fn)
    e_fn = '%s/mafft.log' % '/'.join(out_fn.split('/')[:-1])
    if out_fn.lower().endswith('.gz'):
        command = '%s "%s" "%s" 2> "%s" | gzip -9 > "%s"' % (' '.join(command), in_fn, ref_fn, e_fn, out_fn)
        if logger is not None:
            logger.write("MAFFT command: %s\n" % command)
        call(command, shell=True)
    else:
        command += [in_fn, ref_fn]
        if logger is not None:
            logger.write("MAFFT command: %s\n" % ' '.join(command))
        o = open(out_fn, 'w', buffering=bufsize); e = open(e_fn, 'w', buffering=bufsize); call(command, stdout=o, stderr=e); o.close(); e.close()
    if getsize(out_fn) < MIN_MAFFT_OUTPUT_SIZE:
        raise ValueError("MAFFT crashed. See log: %s" % e_fn)
    if logger is not None:
        logger.write("Multiple sequence alignment written to: %s\n" % out_fn)

# iterate over sequences of FASTA
def iter_fasta(in_fn, bufsize=DEFAULT_BUFSIZE):
    if in_fn == 'stdin':
        in_f = stdin
    elif not in_fn.startswith('/dev/fd/') and not isfile(in_fn):
        raise ValueError("File not found: %s" % in_fn)
    elif in_fn.lower().endswith('.gz'):
        in_f = gopen(in_fn, 'rt')
    else:
        in_f = open(in_fn, 'r', buffering=bufsize)
    seq = None
    for line in in_f:
        l = line.strip()
        if len(l) == 0:
            continue
        if l.startswith('>'):
            if seq is not None:
                yield seq
            seq = ''
        else:
            seq += l.upper()
    in_f.close()

# count bases at each position of MSA
def count_bases(in_fn, logger=None):
    if logger is not None:
        logger.write("Counting bases from: %s\n" % in_fn)
    counts = None # counts[pos][nuc] = count
    for seq_ind, seq in enumerate(iter_fasta(in_fn)):
        if counts is None:
            counts = [{c:0 for c in NUCS_SORT} for _ in range(len(seq))]
        elif len(seq) != len(counts):
            raise ValueError("MSA sequences have differing lengths: %s" % in_fn)
        for i, c in enumerate(seq):
            if c not in NUCS:
                continue
            counts[i][c] += 1
        if logger is not None and (seq_ind+1) % NUM_SEQS_PROGRESS == 0:
            logger.write("Parsed %d sequences...\n" % (seq_ind+1))
    return counts

# write base counts
def write_counts(counts, out_fn, delim='\t', logger=None, bufsize=DEFAULT_BUFSIZE):
    if isfile(out_fn) or isdir(out_fn):
        raise ValueError("Output exists: %s" % out_fn)
    if logger is not None:
        logger.write("Writing base counts to file...\n")
    out_f = open(out_fn, 'w', buffering=bufsize); out_f.write('Position (1-indexed)%s%s\n' % (delim, delim.join(NUCS_SORT)))
    for i, curr in enumerate(counts):
        out_f.write('%d%s%s\n' % (i+1, delim, delim.join(str(curr[c]) if c in curr else '0' for c in NUCS_SORT)))
    out_f.close()
    if logger is not None:
        logger.write("Base counts written to: %s\n" % out_fn)

# compute consensus sequence
def compute_consensus(counts, out_fn, logger=None, bufsize=DEFAULT_BUFSIZE):
    if logger is not None:
        logger.write("Computing consensus sequence...\n")
    seq = ''.join(c for count in counts for c in 'ACGT' if count[c] == max(count.values()))
    out_f = open(out_fn, 'w', buffering=bufsize); out_f.write('>Consensus - ViralPrimerDesign v%s\n%s\n' % (VERSION, seq)); out_f.close()
    if logger is not None:
        logger.write("Consensus sequence written to: %s\n" % out_fn)
    return seq

# compute Shannon entorpy from MSA base counts
def compute_entropies(counts, logger=None):
    if logger is not None:
        logger.write("Computing Shannon entropies...\n")
    ents = [None for _ in range(len(counts))] # ents[pos] = (max base freq, entropy)
    for i, curr in enumerate(counts):
        ent = 0; tot = sum(curr.values()); m = max(curr.values())/tot
        for c in curr:
            if curr[c] != 0:
                p = curr[c]/tot; ent -= (p*log2(p))
        ents[i] = ent
    return ents

# write Shannon entropies
def write_entropies(ents, out_fn, delim='\t', logger=None, bufsize=DEFAULT_BUFSIZE):
    if isfile(out_fn) or isdir(out_fn):
        raise ValueError("Output exists: %s" % out_fn)
    if logger is not None:
        logger.write("Writing entropies to file...\n")
    out_f = open(out_fn, 'w', buffering=bufsize); out_f.write('Position (1-indexed)%sEntropy\n' % delim)
    for i, curr in enumerate(ents):
        out_f.write('%d%s%s\n' % (i+1, delim, curr))
    out_f.close()
    if logger is not None:
        logger.write("Entropies written to: %s\n" % out_fn)

# plot entropies
def plot_entropies(ents, out_fn, logger=None):
    if isfile(out_fn) or isdir(out_fn):
        raise ValueError("Output exists: %s" % out_fn)
    if logger is not None:
        logger.write("Plotting entropies...\n")
    x = [i+1 for i, e in enumerate(ents) if e is not None]
    y = [e for i, e in enumerate(ents) if e is not None]
    fig = relplot(x=x, y=y, kind='line')
    plt.xlabel("Position (1-indexed)")
    plt.ylabel("Shannon Entropy")
    fig.savefig(out_fn, format='pdf', bbox_inches='tight')
    if logger is not None:
        logger.write("Entropy plot saved to: %s\n" % out_fn)

# design primers using Primer3: https://primer3.org/manual.html
def design_primers(
    consensus, primer3_input_fn, primer3_output_fn, logger=None, bufsize=DEFAULT_BUFSIZE,
    window_size=DEFAULT_PRIMER3_WINDOW_SIZE, step_size=DEFAULT_PRIMER3_STEP_SIZE,
    primer_opt_size=DEFAULT_PRIMER3_PRIMER_OPT_SIZE, primer_min_size=DEFAULT_PRIMER3_PRIMER_MIN_SIZE, primer_max_size=DEFAULT_PRIMER3_PRIMER_MAX_SIZE,
    primer_product_min_size=DEFAULT_PRIMER3_PRIMER_PRODUCT_MIN_SIZE, primer_product_max_size=DEFAULT_PRIMER3_PRIMER_PRODUCT_MAX_SIZE,
):
    # run Primer3
    if isfile(primer3_input_fn) or isdir(primer3_input_fn):
        raise ValueError("Output exists: %s" % primer3_input_fn)
    if logger is not None:
        logger.write("Creating Primer3 input file: %s\n" % primer3_input_fn)
    f = open(primer3_input_fn, 'w', buffering=bufsize); num_windows = 0
    for start in range(0, len(consensus)-step_size, step_size):
        f.write("SEQUENCE_ID=%d_%d\n" % (start+1, min(start+window_size,len(consensus))))
        f.write("SEQUENCE_TEMPLATE=%s\n" % consensus[start:start+window_size])
        f.write("PRIMER_TASK=generic\n")
        f.write("PRIMER_PICK_LEFT_PRIMER=1\n")
        f.write("PRIMER_PICK_INTERNAL_OLIGO=1\n")
        f.write("PRIMER_PICK_RIGHT_PRIMER=1\n")
        f.write("PRIMER_OPT_SIZE=%d\n" % primer_opt_size)
        f.write("PRIMER_MIN_SIZE=%d\n" % primer_min_size)
        f.write("PRIMER_MAX_SIZE=%d\n" % primer_max_size)
        f.write("PRIMER_PRODUCT_SIZE_RANGE=%d-%d\n" % (primer_product_min_size, primer_product_max_size))
        f.write("=\n")
        num_windows += 1
    f.close()
    if logger is not None:
        logger.write("Running Primer3...\n")
    o = open(primer3_output_fn, 'w', buffering=bufsize); call(['primer3_core', primer3_input_fn], stdout=o); o.close()
    if getsize(primer3_output_fn) < MIN_PRIMER3_OUTPUT_SIZE:
        raise ValueError("Primer3 crashed")
    if logger is not None:
        logger.write("Primer3 output written to: %s\n" % primer3_output_fn)

    # parse Primer3 output and return
    primers = dict() # primers[(left_start, left_len, right_start, right_len)] = dict storing info about that primer pair
    for entry in open(primer3_output_fn, buffering=bufsize).read().rstrip().rstrip('=').split('\n=\n'):
        data = {k.strip():v.strip() for k,v in [l.split('=') for l in entry.splitlines()]}
        num_primers = int(data['PRIMER_PAIR_NUM_RETURNED'])
        for primer_ind in range(num_primers):
            primer_left_ind, primer_left_len = [int(v) for v in data['PRIMER_LEFT_%d' % primer_ind].split(',')]
            primer_right_ind, primer_right_len = [int(v) for v in data['PRIMER_RIGHT_%d' % primer_ind].split(',')]
            primers_key = (primer_left_ind, primer_left_len, primer_right_ind, primer_right_len)
            tmp = '_%d_' % primer_ind
            if primers_key not in primers:
                primers[primers_key] = {k.replace(tmp, '_'): data[k] for k in data if tmp in k}
    return primers

# write unique primers to FASTA file
def write_unique_primers(primers, unique_primer_fn, logger=None, bufsize=DEFAULT_BUFSIZE):
    if isfile(unique_primer_fn) or isdir(unique_primer_fn):
        raise ValueError("Output exists: %s" % unique_primer_fn)
    if logger is not None:
        logger.write("Generating unique primers...\n")
    primer_seqs = {curr[k] for curr in primers.values() for k in curr if k.endswith('_SEQUENCE')}
    if logger is not None:
        logger.write("Writing unique primers to FASTA file: %s\n" % unique_primer_fn)
    f = open(unique_primer_fn, 'w', buffering=bufsize)
    for seq in primer_seqs:
        f.write(">%s\n%s\n" % (seq, seq))
    f.close()

# blast all unique primer sequences
def blast_primer_seqs(unique_primer_fn, out_fn, word_size=DEFAULT_BLAST_WORD_SIZE, logger=None, bufsize=DEFAULT_BUFSIZE):
    command = ['blastn', '-query', unique_primer_fn, '-db', 'nt', '-remote', '-word_size', str(word_size), '-out', out_fn, '-outfmt', '6']# sallseqid sallacc evalue bitscore score pident']
    if logger is not None:
        logger.write("BLAST Command: %s\n" % ' '.join(command))
    e_fn = '%s/blastn.log' % '/'.join(out_fn.split('/')[:-1]); e = open(e_fn, 'w'); call(command, stderr=e); e.close()

# main program
def main():
    args = parse_args(); mkdir(args.outdir); logger = Log('%s/log.txt' % args.outdir, quiet=args.quiet)
    logger.write("=== ViralPrimerDesign v%s ===\n" % VERSION)
    logger.write("Command: %s\n" % ' '.join(argv))
    if not args.skip_alignment:
        orig_sequences = args.sequences; args.sequences = '%s/sequences.aln.gz' % args.outdir
        align_mafft(orig_sequences, args.reference, args.sequences, threads=args.threads, logger=logger)
    counts = count_bases(args.sequences, logger=logger)
    write_counts(counts, '%s/counts.tsv' % args.outdir, logger=logger)
    consensus = compute_consensus(counts, '%s/consensus.fas' % args.outdir, logger=logger)
    ents = compute_entropies(counts, logger=logger)
    write_entropies(ents, '%s/entropies.tsv' % args.outdir, logger=logger)
    plot_entropies(ents, '%s/entropies.pdf' % args.outdir, logger=logger)
    primers = design_primers(
        consensus, '%s/primer3_input.txt' % args.outdir, '%s/primer3_output.txt' % args.outdir, logger=logger,
        window_size=args.primer3_window_size, step_size=args.primer3_step_size,
        primer_opt_size=args.primer3_primer_opt_size, primer_min_size=args.primer3_primer_min_size, primer_max_size=args.primer3_primer_max_size,
        primer_product_min_size=args.primer3_primer_product_min_size, primer_product_max_size=args.primer3_primer_product_max_size,
    )
    unique_primer_fn = '%s/unique_primers.fas' % args.outdir
    write_unique_primers(primers, unique_primer_fn, logger=logger)
    blast_results = blast_primer_seqs(unique_primer_fn, '%s/blastn.tsv' % args.outdir, word_size=args.blast_word_size, logger=logger)
    #print(len(primer_seqs))
    #print(list(primers.keys())[0])
    #print(primers[list(primers.keys())[0]])
    #print(list(primer_seqs)[0])
    # blastn example: https://bioinformatics.stackexchange.com/a/19796/1115

# run from CLI
if __name__ == "__main__":
    main()
