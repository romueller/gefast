/*
 * SCT-PJ
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
lenSeqs_t computeClassicFull(const std::string& s, const std::string& t);


/*
 * Space-reduced version of classic dynamic-programming scheme.
 *  - keeps always only a single row of the matrix
 *  - computes exact edit distance
 */
lenSeqs_t computeClassicRow(const std::string& s, const std::string& t, lenSeqs_t* M);




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
lenSeqs_t computeBoundedFull(const std::string& s, const std::string& t, const lenSeqs_t bound);


/*
 * Space-reduced dynamic-programming scheme restricted to some diagonals.
 *  - overall only some diagonals (number depends on bound) are computed,
 *      but keeps always only a single row of the matrix
 *  - computes exact edit distance if d_e(s,t) <= bound and returns bound + 1 otherwise
 *  - early termination if all values in row computed last already larger than bound
 */
lenSeqs_t computeBoundedRow(const std::string& s, const std::string& t, const lenSeqs_t bound, lenSeqs_t* M);


/*
 * Dynamic-programming scheme restricted to some diagonals.
 *  - only some diagonals (number depends on bound) filled / computed
 *  - computes exact edit distance if d_e(s,t) <= bound and returns bound + 1 otherwise
 *  - early termination if all values in row computed last already larger than bound
 */
lenSeqs_t computeBoundedFullSlim(const std::string& s, const std::string& t, const lenSeqs_t bound);


/*
 * Space-reduced dynamic-programming scheme restricted to some diagonals.
 *  - overall only some diagonals (number depends on bound) are computed,
 *      but keeps always only a section of a single row of the matrix
 *  - computes exact edit distance if d_e(s,t) <= bound and returns bound + 1 otherwise
 *  - early termination if all values in row computed last already larger than bound
 */
lenSeqs_t computeBoundedRowSlim(const std::string& s, const std::string& t, const lenSeqs_t bound, lenSeqs_t* M);




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
lenSeqs_t computeLengthAwareFull(const std::string& s, const std::string& t, const lenSeqs_t bound);


/*
 * Space-reduced dynamic-programming scheme restricted to some diagonals.
 *  - overall only some diagonals (number depends on bound) are computed,
 *      but keeps always only a single row of the matrix
 *  - uses tighter bounds by considering length difference of s and t
 *  - computes exact edit distance if d_e(s,t) <= bound and returns bound + 1 otherwise
 *  - improved early termination if all expected edit distance values for the row computed last already larger than bound
 */
lenSeqs_t computeLengthAwareRow(const std::string& s, const std::string& t, const lenSeqs_t bound, lenSeqs_t* M);


/*
 * Dynamic-programming scheme restricted to some diagonals.
 *  - only some diagonals (number depends on bound) filled / computed
 *  - uses tighter bounds by considering length difference of s and t
 *  - computes exact edit distance if d_e(s,t) <= bound and returns bound + 1 otherwise
 *  - improved early termination if all expected edit distance values for the row computed last already larger than bound
 */
lenSeqs_t computeLengthAwareFullSlim(const std::string& s, const std::string& t, const lenSeqs_t bound);


/*
 * Space-reduced dynamic-programming scheme restricted to some diagonals.
 *  - overall only some diagonals (number depends on bound) are computed,
 *      but keeps always only a section of a single row of the matrix
 *  - uses tighter bounds by considering length difference of s and t
 *  - computes exact edit distance if d_e(s,t) <= bound and returns bound + 1 otherwise
 *  - improved early termination if all expected edit distance values for the row computed last already larger than bound
 */
lenSeqs_t computeLengthAwareRowSlim(const std::string& s, const std::string& t, const lenSeqs_t bound, lenSeqs_t* M);


// compute edit distance of all incoming candidates
// verification lasts until everything in the buffer is worked off and it signals that no new candidates will be inserted
void verify(const AmpliconCollection& ac, Matches& mat, Buffer<Candidate>& buf, lenSeqs_t width, lenSeqs_t t);

}
}


#endif //SCT_PJ_VERIFICATION_HPP