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

#ifndef SCT_PJ_SEGMENTFILTER_HPP
#define SCT_PJ_SEGMENTFILTER_HPP

#include "Base.hpp"
#include "Buffer.hpp"


namespace SCT_PJ {
namespace SegmentFilter {

    // describes a set of substrings chosen for comparison with a segment
    struct Substrings {

        lenSeqs_t first; //start position of first substring to be checked
        lenSeqs_t last; //start position of last substring to be checked
        lenSeqs_t len; //common length of all substrings to be checked

        Substrings() {

            first = 0;
            last = 0;
            len = 0;

        }

        Substrings(lenSeqs_t fp, lenSeqs_t lp, lenSeqs_t l) {

            first = fp;
            last = lp;
            len = l;

        }

    };

    // describes a set of segments to be indexed
    typedef std::vector<std::pair<lenSeqs_t, lenSeqs_t>> Segments;


    // select 'substrings' (MMASS) for segment filter
    Substrings selectSubstrs(const lenSeqs_t selfLen, const lenSeqs_t partnerLen, const lenSeqs_t segIndex, const lenSeqs_t t, const lenSeqs_t k);
    Substrings selectSubstrsBackward(const lenSeqs_t selfLen, const lenSeqs_t partnerLen, const lenSeqs_t segIndex, const lenSeqs_t t, const lenSeqs_t k);

    // select 'segments' (to be stored in a parameter) for indexing step
    void selectSegments(Segments& segments, const lenSeqs_t seqLen, const lenSeqs_t t, const lenSeqs_t k);


    //TODO avoid counting multiple times per segment (during the full filter step), e.g. by collecting the found partners in a set and +1 per set element at the end of processing the segment

    // (forward) segment filter for the general pigeonhole principle (t + k segments, k segments must be matched)
    void filterForward(const AmpliconCollection& ac, const Subpool& sp, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k);

    // (backward) segment filter for the general pigeonhole principle (t + k segments, k segments must be matched)
    void filterBackward(const AmpliconCollection& ac, const Subpool& sp, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k);

    // (forward) segment filter for the general pigeonhole principle (t + k segments, k segments must be matched) + pipelined backward filtering
    void filterForwardBackward(const AmpliconCollection& ac, const Subpool& sp, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k);

    // (backward) segment filter for the general pigeonhole principle (t + k segments, k segments must be matched) + pipelined forward filtering
    void filterBackwardForward(const AmpliconCollection& ac, const Subpool& sp, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k);

    void filter(const AmpliconCollection& ac, const Subpool& sp, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k, const int mode);


    void writeMatches(std::string oFile, AmpliconPools& pools, std::vector<Matches*>& allMatches);

    void writeMatchesOneWay(std::string oFile, AmpliconPools& pools, std::vector<Matches*>& allMatches);

//    //old methods not capable of dealing with subpools
//    void filterForward(AmpliconCollection& ac, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k);
//    void filterBackward(AmpliconCollection& ac, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k);
//    void filterForwardBackward(AmpliconCollection& ac, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k);
//    void filterBackwardForward(AmpliconCollection& ac, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k);
}
}

#endif //SCT_PJ_SEGMENTFILTER_HPP