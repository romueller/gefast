/*
 * GeFaST
 *
 * Copyright (C) 2016 - 2017 Robert Mueller
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact: Robert Mueller <romueller@techfak.uni-bielefeld.de>
 * Faculty of Technology, Bielefeld University,
 * PO box 100131, DE-33501 Bielefeld, Germany
 */

#ifndef GEFAST_VERIFICATIONGOTOH_HPP
#define GEFAST_VERIFICATIONGOTOH_HPP

#include "Base.hpp"
#include "Buffer.hpp"
#include "Relation.hpp"
#include "Utility.hpp"

namespace GeFaST {
namespace Verification {



/*
 * Representation of scoring function with affine gap costs.
 * Converted towards minimisation with a fixed match reward of 0.
 * The conversion is based on the ideas presented by Smith et al.
 * (https://www.ncbi.nlm.nih.gov/pubmed/7334527) in section "The Measures" (p. 39).
 */
struct Scoring {

//    val_t rewardMatch; // match reward
    val_t penMismatch; // mismatch penalty
    val_t penOpen; // gap-opening penalty
    val_t penExtend; // gap-extension penalty

    Scoring () {
        penMismatch = penOpen = penExtend = 0;
    }

    Scoring (long long m, long long p, long long g, long long e) {

        penMismatch = 2 * (m - p);
        penOpen = -2 * g;
        penExtend = m - 2 * e;

        auto penFactor = gcd(gcd(penMismatch, penOpen), penExtend);

        penMismatch /= penFactor;
        penOpen /= penFactor;
        penExtend /= penFactor;

    }

};


// "backtracking flags" for Gotoh's algorithm using three matrices
const char DIAGONAL_IN_D = 1;
const char JUMP_TO_Q = 2;
const char JUMP_TO_P = 4;
const char UP_TO_D = 8;
const char UP_IN_P = 16;
const char LEFT_TO_D = 32;
const char LEFT_IN_Q = 64;

const char IN_D = 0;
const char IN_P = 1;
const char IN_Q = 2;


//===========================================================
//                  Exact computation
//===========================================================
// Based on:
// Gotoh (1982), An improved algorithm for matching biological sequences, and
// Rose (2012), Bioinformatics I (lecture), http://www.bioinf.uni-freiburg.de/Lehre/Courses/2013_SS/V_Bioinformatik_1/lecture4.pdf (accessed on 22 Sep 2016)

/*
 * Classic dynamic-programming scheme for global alignment (score only) with affine gap costs.
 *  - full matrices
 *  - computes score of best alignment
 */
val_t computeGotohScoreFull(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring);

/*
 * Classic dynamic-programming scheme for global alignment with affine gap costs.
 *  - full matrices
 *  - computes best alignment and its score
 */
std::pair<val_t, std::string> computeGotohAlignmentFull(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring);


/*
 * Classic dynamic-programming scheme for global alignment (score only) with affine gap costs.
 *  - full matrices
 *  - computes number of differences (mismatches, insertions, deletions) in the best alignment
 */
lenSeqs_t computeGotohFull(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring);

/*
 * Space-reduced version of above classic dynamic-programming scheme.
 *  - keeps always only a single row of every matrix
 *  - computes number of differences (mismatches, insertions, deletions) in the best alignment
 */
lenSeqs_t computeGotohRow(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring,
                          val_t* D, val_t* P, val_t* Q, char* BT, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP, lenSeqs_t* cntDiffsQ);

/*
 * The same as computeGotohRow(...), but with early termination if detected that all paths imply too many differences (returns bound + 1 in this case).
 */
lenSeqs_t computeGotohEarlyRow(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const lenSeqs_t bound, const Scoring& scoring,
                               val_t* D, val_t* P, val_t* Q, char* BT, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP, lenSeqs_t* cntDiffsQ);



//===========================================================
//                  Bounded computation
//===========================================================
/*
 * Space-reduced version of above classic dynamic-programming scheme with early termination and restricted to some diagonals (number depends on bound).
 * Uses ideas from Li et al. (2013), A partition-based method for string similarity joins with edit-distance constraints.
 *
 * Contains several additional optimisations:
 *  - Handles one-diagonal case separately.
 *  - Renders the BT array redundant.
 *  - Removes case distinction between left-most, right-most and inner diagonals through a few extra assignments.
 *  - Uses casts and labs(...) instead of ternary operator in computation of early-flag.
 *  - Combines computation of maximum values for next entry in arrays D, P and Q with case distinction for differences-counting arrays.
 *  - Substitutes the arrays Q and cntDiffsQ with two integer variables.
 *  - Avoids some recomputations by using a few additional integer variables.
 */
lenSeqs_t computeGotohLengthAwareEarlyRow(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const lenSeqs_t bound,
                                          const Scoring& scoring, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP);


/*
 * Computes the number of differences in best alignments of all incoming candidates.
 * The verification lasts until everything in the buffer is worked off and it signals that no new candidates will be inserted.
 */
void verifyGotoh(const AmpliconCollection& ac, Matches& mat, Buffer<Candidate>& buf, lenSeqs_t width, lenSeqs_t t, const Scoring& scoring);



//===========================================================
//                  CIGAR computation
//===========================================================

struct AlignmentInformation {

    std::string cigar;
    lenSeqs_t length;
    lenSeqs_t numDiffs;
//    val_t score;

    AlignmentInformation(const std::string& c, const lenSeqs_t l, const lenSeqs_t n/*, const val_t s*/) {

        cigar = c;
        length = l;
        numDiffs = n;
//        score = s;

    }

};

/*
 * Computes one optimal global alignment with affine gap costs and returns that alignment in the CIGAR format,
 * its length and the number of differences (mismatches, insertions, deletions) in it.
 */
AlignmentInformation computeGotohCigarFull(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring);
AlignmentInformation computeGotohCigarFull1(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring,
                                            val_t* D, val_t* P, val_t* Q, char* BT);

/*
 * Computes one optimal global alignment with affine gap costs and returns that alignment in the CIGAR format,
 * its length and the number of differences (mismatches, insertions, deletions) in it.
 *
 * Only one row of each matrix (except the one for backtracking) is kept in memory.
 */
AlignmentInformation computeGotohCigarRow(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring);
AlignmentInformation computeGotohCigarRow1(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring,
                                           val_t* D, val_t* P, char* BT);

}
}

#endif //GEFAST_VERIFICATIONGOTOH_HPP
