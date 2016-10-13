/*
 * SCT-PJ
 *
 * Copyright (C) 2016 Robert Mueller
 *
 * TODO add licence text (e.g. GNU (Affero) General Public License)
 *
 * Contact: Robert Mueller <romueller@techfak.uni-bielefeld.de>
 * Faculty of Technology, Bielefeld University,
 * PO box 100131, DE-33501 Bielefeld, Germany
 */

#ifndef SCT_PJ_VERIFICATIONGOTOH_HPP
#define SCT_PJ_VERIFICATIONGOTOH_HPP

#include "Base.hpp"
#include "Buffer.hpp"
#include "Relation.hpp"

namespace SCT_PJ {
namespace Verification {


typedef long val_t;
const val_t NEG_INF = INT16_MIN;

// two-dimensional array types for below computations
typedef std::vector<std::vector<val_t>> ValArray; // stores alignment scores
typedef std::vector<std::vector<lenSeqs_t>> CountArray; // stores counts of differences
typedef std::vector<std::vector<char>> BacktrackingArray; // stores information for alignment reconstruction

/**
 * Representation of scoring function with affine gap costs.
 */
struct Scoring {

    val_t rewardMatch; // match reward
    val_t penMismatch; // mismatch penalty
    val_t penOpen; // gap-opening penalty
    val_t penExtend; // gap-extension penalty

    Scoring () {
        rewardMatch = penMismatch = penOpen = penExtend = 0;
    }

    Scoring (val_t m, val_t p, val_t g, val_t e) {

        rewardMatch = m;
        penMismatch = p;
        penOpen = g;
        penExtend = e;

    }

};

/*
// reduces number of arguments to the computeGotoh...() methods, but is slower than passing several arrays
struct OneRowArrays {

    std::vector<val_t> D;
    std::vector<val_t> P;
    std::vector<val_t> Q;

    std::vector<char> BT;

    std::vector<lenSeqs_t> cntDiffs;
    std::vector<lenSeqs_t> cntDiffsP;
    std::vector<lenSeqs_t> cntDiffsQ;

    OneRowArrays() {

        D = std::vector<val_t>(0);
        P = std::vector<val_t>(0);
        Q = std::vector<val_t>(0);

        BT = std::vector<char>(0);

        cntDiffs = std::vector<lenSeqs_t>(0);
        cntDiffsP = std::vector<lenSeqs_t>(0);
        cntDiffsQ = std::vector<lenSeqs_t>(0);

    }

    OneRowArrays(const lenSeqs_t size) {

        D = std::vector<val_t>(size);
        P = std::vector<val_t>(size);
        Q = std::vector<val_t>(size);

        BT = std::vector<char>(size);

        cntDiffs = std::vector<lenSeqs_t>(size);
        cntDiffsP = std::vector<lenSeqs_t>(size);
        cntDiffsQ = std::vector<lenSeqs_t>(size);

    }

};
*/

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

/**
 * Classic dynamic-programming scheme for global alignment (score only) with affine gap costs.
 *  - full matrices
 *  - computes score of best alignment
 */
val_t computeGotohScoreFull(const std::string& s, const std::string& t, const Scoring& scoring);

/**
 * Classic dynamic-programming scheme for global alignment with affine gap costs.
 *  - full matrices
 *  - computes best alignment and its score
 */
std::pair<val_t, std::string> computeGotohAlignmentFull(const std::string& s, const std::string& t, const Scoring& scoring);


/**
 * Classic dynamic-programming scheme for global alignment (score only) with affine gap costs.
 *  - full matrices
 *  - computes number of differences (mismatches, insertions, deletions) in the best alignment
 */
lenSeqs_t computeGotohFull(const std::string& s, const std::string& t, const Scoring& scoring);

/**
 * Space-reduced version of above classic dynamic-programming scheme.
 *  - keeps always only a single row of every matrix
 *  - computes number of differences (mismatches, insertions, deletions) in the best alignment
 */
lenSeqs_t computeGotohRow(const std::string& s, const std::string& t, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], char BT[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]);

/**
 * The same as computeGotohRow(...), but with early termination if detected that all paths imply too many differences (returns bound + 1 in this case).
 */
lenSeqs_t computeGotohEarlyRow(const std::string& s, const std::string& t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], char BT[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]);



//===========================================================
//                  Bounded computation
//===========================================================

/**
 * Space-reduced version of above classic dynamic-programming scheme with early termination and restricted to some diagonals (number depends on bound).
 */
lenSeqs_t computeGotohBoundedEarlyRow(const std::string& s, const std::string& t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], char BT[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]);

/**
 * Space-reduced version of above classic dynamic-programming scheme with early termination and restricted to some diagonals (number depends on bound and length difference).
 * More restricted w.r.t. the number of diagonals than computeGotohBoundedEarlyRow(...).
 *
 * Uses ideas from Li et al. (2013), A partition-based method for string similarity joins with edit-distance constraints.
 */
lenSeqs_t computeGotohLengthAwareEarlyRow(const std::string& s, const std::string& t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], char BT[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]);

/**
 * Improved version of computeGotohLengthAwareEarlyRow(...). Handles one-diagonal case separately.
 */
lenSeqs_t computeGotohLengthAwareEarlyRow2(const std::string& s, const std::string& t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], char BT[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]);

/**
 * Improved version of computeGotohLengthAwareEarlyRow2(...). Renders the BT array redundant.
 */
lenSeqs_t computeGotohLengthAwareEarlyRow3(const std::string& s, const std::string& t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]);

/**
 * Improved version of computeGotohLengthAwareEarlyRow3(...). Removes case distinction between left-most, right-most and inner diagonals through a few extra assignments.
 * Also uses casts and labs(...) instead of ternary operator in computation of early-flag.
 */
lenSeqs_t computeGotohLengthAwareEarlyRow4(const std::string& s, const std::string& t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]);

/**
 * Improved version of computeGotohLengthAwareEarlyRow4(...). Avoids recomputing (bound -/+ delta) / 2 often by storing it in local variables.
 */
lenSeqs_t computeGotohLengthAwareEarlyRow5(const std::string& s, const std::string& t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]);

/**
 * Improved version of computeGotohLengthAwareEarlyRow4(...). Combines computation of maximum values for next entry in arrays D, P and Q with case distinction for differences-counting arrays.
 *
 * First tests have indicated that this version is the fastest one.
 */
lenSeqs_t computeGotohLengthAwareEarlyRow6(const std::string& s, const std::string& t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]);


// compute differences in best alignments of all incoming candidates
// verification lasts until everything in the buffer is worked off and it signals that no new candidates will be inserted
void verifyGotoh(const AmpliconCollection& ac, Matches& mat, Buffer<Candidate>& buf, lenSeqs_t t, const Scoring& scoring);
}
}

#endif //SCT_PJ_VERIFICATIONGOTOH_HPP
