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

#ifndef SCT_PJ_VERIFICATION_HPP
#define SCT_PJ_VERIFICATION_HPP

#include <iostream>

#include "Base.hpp"
#include "Buffer.hpp"
#include "Relation.hpp"


namespace SCT_PJ {
namespace Verification {

//===========================================================
//                  Exact computation
//===========================================================


/*
 * Classic dynamic-programming scheme for edit distance.
 *  - full matrix
 *  - computes exact edit distance
 */
// classic_f / compute_edit_distance_full
lenSeqs_t computeClassicFull(const std::string& s, const std::string& t);


/*
 * Space-reduced version of classic dynamic-programming scheme.
 *  - keeps always only a single row of the matrix
 *  - computes exact edit distance
 */
// classic_1 / compute_edit_distance_one_row
lenSeqs_t computeClassicRow(const std::string& s, const std::string& t, lenSeqs_t M[]);




//===========================================================
//                  Bounded computation
//===========================================================

/*
 * Helper function for mapping column indices when using slim dynamic programming schemes.
 */
lenSeqs_t mapColIndex(const lenSeqs_t j, const lenSeqs_t i, const lenSeqs_t bound);


/*
 * Dynamic-programming scheme restricted to some diagonals.
 *  - full matrix but only some diagonals (number depends on bound) filled
 *  - computes exact edit distance if d_e(s,t) <= bound and returns bound + 1 otherwise
 *  - early termination if all values in row computed last already larger than bound
 */
// bounded_f / compute_bounded_edit_distance_with_early_termination
lenSeqs_t computeBoundedFull(const std::string& s, const std::string& t, const lenSeqs_t bound);


/*
 * Space-reduced dynamic-programming scheme restricted to some diagonals.
 *  - overall only some diagonals (number depends on bound) are computed,
 *      but keeps always only a single row of the matrix
 *  - computes exact edit distance if d_e(s,t) <= bound and returns bound + 1 otherwise
 *  - early termination if all values in row computed last already larger than bound
 */
// bounded_1 / compute_bounded_edit_distance_with_early_termination_one_row
lenSeqs_t computeBoundedRow(const std::string& s, const std::string& t, const lenSeqs_t bound, lenSeqs_t M[]);


/*
 * Dynamic-programming scheme restricted to some diagonals.
 *  - only some diagonals (number depends on bound) filled / computed
 *  - computes exact edit distance if d_e(s,t) <= bound and returns bound + 1 otherwise
 *  - early termination if all values in row computed last already larger than bound
 */
// bounded_fs / compute_bounded_edit_distance_with_early_termination_slim
lenSeqs_t computeBoundedFullSlim(const std::string& s, const std::string& t, const lenSeqs_t bound);


/*
 * Space-reduced dynamic-programming scheme restricted to some diagonals.
 *  - overall only some diagonals (number depends on bound) are computed,
 *      but keeps always only a section of a single row of the matrix
 *  - computes exact edit distance if d_e(s,t) <= bound and returns bound + 1 otherwise
 *  - early termination if all values in row computed last already larger than bound
 */
// bounded_1s / compute_bounded_edit_distance_with_early_termination_one_row_slim
lenSeqs_t computeBoundedRowSlim(const std::string& s, const std::string& t, const lenSeqs_t bound, lenSeqs_t M[]);




//===========================================================
//            Length-aware bounded computation
//===========================================================
// Based on:
// Li et al. (2013), A partition-based method for string similarity joins with edit-distance constraints

/*
 * Dynamic-programming scheme restricted to some diagonals.
 *  - full matrix but only some diagonals (number depends on bound) filled
 *  - uses tighter bounds by considering length difference of s and t
 *  - computes exact edit distance if d_e(s,t) <= bound and returns bound + 1 otherwise
 *  - improved early termination if all expected edit distance values for the row computed last already larger than bound
 */
// l-aware_f / compute_length_aware_bounded_edit_distance_with_early_termination
lenSeqs_t computeLengthAwareFull(const std::string& s, const std::string& t, const lenSeqs_t bound);


/*
 * Space-reduced dynamic-programming scheme restricted to some diagonals.
 *  - overall only some diagonals (number depends on bound) are computed,
 *      but keeps always only a single row of the matrix
 *  - uses tighter bounds by considering length difference of s and t
 *  - computes exact edit distance if d_e(s,t) <= bound and returns bound + 1 otherwise
 *  - improved early termination if all expected edit distance values for the row computed last already larger than bound
 */
// l-aware_1 / compute_length_aware_bounded_edit_distance_with_early_termination_one_row
lenSeqs_t computeLengthAwareRow(const std::string& s, const std::string& t, const lenSeqs_t bound, lenSeqs_t M[]);


/*
 * Dynamic-programming scheme restricted to some diagonals.
 *  - only some diagonals (number depends on bound) filled / computed
 *  - uses tighter bounds by considering length difference of s and t
 *  - computes exact edit distance if d_e(s,t) <= bound and returns bound + 1 otherwise
 *  - improved early termination if all expected edit distance values for the row computed last already larger than bound
 */
// l-aware_fs / compute_length_aware_bounded_edit_distance_with_early_termination_slim
lenSeqs_t computeLengthAwareFullSlim(const std::string& s, const std::string& t, const lenSeqs_t bound);


/*
 * Space-reduced dynamic-programming scheme restricted to some diagonals.
 *  - overall only some diagonals (number depends on bound) are computed,
 *      but keeps always only a section of a single row of the matrix
 *  - uses tighter bounds by considering length difference of s and t
 *  - computes exact edit distance if d_e(s,t) <= bound and returns bound + 1 otherwise
 *  - improved early termination if all expected edit distance values for the row computed last already larger than bound
 */
// l-aware_1s / compute_length_aware_bounded_edit_distance_with_early_termination_one_row_slim
lenSeqs_t computeLengthAwareRowSlim(const std::string& s, const std::string& t, const lenSeqs_t bound, lenSeqs_t M[]);


// compute edit distance of all incoming candidates
// verification lasts until everything in the buffer is worked off and it signals that no new candidates will be inserted
void verify(const AmpliconCollection& ac, Matches& mat, Buffer<Candidate>& buf, lenSeqs_t t);

}
}


#endif //SCT_PJ_VERIFICATION_HPP