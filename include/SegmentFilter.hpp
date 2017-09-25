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

#ifndef GEFAST_SEGMENTFILTER_HPP
#define GEFAST_SEGMENTFILTER_HPP

#include "Base.hpp"
#include "Buffer.hpp"
#include "VerificationGotoh.hpp"


namespace GeFaST {
namespace SegmentFilter {

    /*
     * Concerning all four filter variants:
     * Currently, a segment contributes with each matched substring (instead of at most one), which can make an indexed sequence a candidate even though not enough different segments are matched.
     * This would increase the workload of the verification stage.
     * However, some runtime tests with means preventing this multiple counting increased the overall runtime, as the costs for avoiding multiple counting seem to exceed the additional verification costs.
     * Therefore, multiple counting is not prevented (as it also has no influence on the final result).
     *
     * The methods with the suffix 'Directly' verify the candidates themselves directly when they occur and do not hand them over to verifier threads through a buffer.
     *
     * All filter methods assume that the amplicons are sorted by increasing sequence length.
     */

    // (forward) segment filter for the general pigeonhole principle (t + k segments, k segments must be matched)
    void filterForward(const AmpliconCollection& ac, const Subpool& sp, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k);
    void filterForwardDirectly(const AmpliconCollection& ac, const Subpool& sp, Matches& matches, const lenSeqs_t t, const lenSeqs_t k, const bool useScore, const Verification::Scoring& scoring);

    // (backward) segment filter for the general pigeonhole principle (t + k segments, k segments must be matched)
    void filterBackward(const AmpliconCollection& ac, const Subpool& sp, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k);
    void filterBackwardDirectly(const AmpliconCollection& ac, const Subpool& sp, Matches& matches, const lenSeqs_t t, const lenSeqs_t k, const bool useScore, const Verification::Scoring& scoring);

    // (forward) segment filter for the general pigeonhole principle (t + k segments, k segments must be matched) + pipelined backward filtering
    void filterForwardBackward(const AmpliconCollection& ac, const Subpool& sp, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k);
    void filterForwardBackwardDirectly(const AmpliconCollection& ac, const Subpool& sp, Matches& matches, const lenSeqs_t t, const lenSeqs_t k, const bool useScore, const Verification::Scoring& scoring);

    // (backward) segment filter for the general pigeonhole principle (t + k segments, k segments must be matched) + pipelined forward filtering
    void filterBackwardForward(const AmpliconCollection& ac, const Subpool& sp, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k);
    void filterBackwardForwardDirectly(const AmpliconCollection& ac, const Subpool& sp, Matches& matches, const lenSeqs_t t, const lenSeqs_t k, const bool useScore, const Verification::Scoring& scoring);

    void filter(const AmpliconCollection& ac, const Subpool& sp, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k, const int mode);
    void filterDirectly(const AmpliconCollection& ac, const Subpool& sp, Matches& matches, const lenSeqs_t t, const lenSeqs_t k, const int mode, const bool useScore, const Verification::Scoring& scoring);

}
}

#endif //GEFAST_SEGMENTFILTER_HPP