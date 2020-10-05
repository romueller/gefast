#!/usr/bin/python3

# Conversion of maximising scoring function with affine gap costs
# to a minimising scoring function with a fixed match reward of 0.
# The conversion is based on the ideas presented by Smith et al.
# (https://www.ncbi.nlm.nih.gov/pubmed/7334527) in section "The Measures" (p. 39)
# and guarantees that alignments obtained by maximisation with the old and
# minimising with the new scoring function are identical.
#
# This script is provided to obtain the converted scoring function beforehand
# in order to choose sensible thresholds when using modes employing
# alignment-score thresholds.
#
# Requires Python 3.5 or higher.

import math


# Convert the scoring functions and show both next to each other.
def convert(match, mismatch, gap_open, gap_extend):

    penalty_mismatch = 2 * (match - mismatch)
    penalty_open = -2 * gap_open
    penalty_extend = match - 2 * gap_extend

    penalty_factor = math.gcd(math.gcd(penalty_mismatch, penalty_open), penalty_extend)

    penalty_mismatch = penalty_mismatch // penalty_factor
    penalty_open = penalty_open // penalty_factor
    penalty_extend = penalty_extend // penalty_factor

    print('=== Conversion of scoring function ===')
    print('\t\tGiven\t\tConverted')
    print('Match:\t\t%i\t->\t%i' % (match, 0))
    print('Mismatch:\t%i\t->\t%i' % (mismatch, penalty_mismatch))
    print('Open gap:\t%i\t->\t%i' % (gap_open, penalty_open))
    print('Extend gap:\t%i\t->\t%i' % (gap_extend, penalty_extend))


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('match', type = int, help = 'reward for matching nucleotides (positive integer)')
    argparser.add_argument('mismatch', type = int, help = 'penalty for mismatching nucleotides (negative integer)')
    argparser.add_argument('gap_open', type = int, help = 'penalty for opening a gap (negative integer)')
    argparser.add_argument('gap_extend', type = int, help = 'penalty for extending a gap (negative integer)')
    args = argparser.parse_args()

    convert(args.match, args.mismatch, args.gap_open, args.gap_extend)