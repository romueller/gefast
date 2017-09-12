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

#include <fstream>
#include <unordered_set>

#include "../include/Relation.hpp"
#include "../include/SegmentFilter.hpp"
#include "../include/Verification.hpp"


namespace GeFaST {

void SegmentFilter::filterForward(const AmpliconCollection& ac, const Subpool& sp, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k) {

    RollingIndices<InvertedIndex> indices(t + 1, t + k, true);

    Substrings substrs[t + 1][t + k];
    Segments segments(t + k);
    StringIteratorPair sip;

    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    std::vector<Candidate> candColl;

    lenSeqs_t seqLen = 0;
    numSeqs_t curIntId = sp.beginIndex;

    // process index-only amplicons
    for (; curIntId < sp.beginMatch; curIntId++) {

        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            selectSegments(segments, seqLen, t, k);

        }

        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

    }



    // process remaining amplicons by matching and indexing
    seqLen = 0; // ensures that the substring information are computed
    for (; curIntId < sp.end; curIntId++) { // for every amplicon in the collection...

        // on reaching new length group, open new inverted indices and discard those no longer needed ...
        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            // ... and determine position information shared by all amplicons of this length
            for (lenSeqs_t lenDiff = 0; lenDiff < t + 1; lenDiff++) {
                for (lenSeqs_t segmentIndex = 0; segmentIndex < t + k; segmentIndex++) {
                    substrs[lenDiff][segmentIndex] = selectSubstrs(seqLen, seqLen - lenDiff, segmentIndex, t, k);
                }
            }
            selectSegments(segments, seqLen, t, k);

        }


        for (lenSeqs_t len = (seqLen > t) * (seqLen - t); len <= seqLen; len++) { // ... consider already indexed seqs with a feasible length...

            candCnts.clear();

            for (lenSeqs_t i = 0; i < t + k; i++) { // ... and apply segment filter for each segment

                Substrings& subs = substrs[seqLen - len][i];
                InvertedIndex& inv = indices.getIndex(len, i);
                sip.first = ac[curIntId].seq + subs.first;
                sip.second = sip.first + subs.len;

                for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++, sip.first++, sip.second++) {
                    inv.addLabelCountsOf(sip, candCnts);
                }

            }

            // general pigeonhole principle: for being a candidate, at least k segments have to be matched
            for (auto candIter = candCnts.begin(); candIter != candCnts.end(); candIter++) {

                if (candIter->second >= k) {
                    candColl.push_back(Candidate(curIntId, candIter->first));
                }

            }

        }

        // index sequence
        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

        cands.push(candColl);
        candColl = std::vector<Candidate>();

    }

}


void SegmentFilter::filterForwardDirectly(const AmpliconCollection& ac, const Subpool& sp, Matches& matches, const lenSeqs_t t, const lenSeqs_t k, const bool useScore, const Verification::Scoring& scoring) {

    RollingIndices<InvertedIndex> indices(t + 1, t + k, true);

    Substrings substrs[t + 1][t + k];
    Segments segments(t + k);
    StringIteratorPair sip;

    lenSeqs_t M[useScore ? 1 : ac.back().len + 1];
    val_t D[useScore? ac.back().len + 1 : 1];
    val_t P[useScore? ac.back().len + 1 : 1];
    lenSeqs_t cntDiffs[useScore? ac.back().len + 1 : 1];
    lenSeqs_t cntDiffsP[useScore? ac.back().len + 1 : 1];

    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;

    lenSeqs_t seqLen = 0;
    numSeqs_t curIntId = sp.beginIndex;
    lenSeqs_t dist;

    // process index-only amplicons
    for (; curIntId < sp.beginMatch; curIntId++) {

        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            selectSegments(segments, seqLen, t, k);

        }

        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

    }



    // process remaining amplicons by matching and indexing
    seqLen = 0; // ensures that the substring information are computed
    for (; curIntId < sp.end; curIntId++) { // for every amplicon in the collection...

        // on reaching new length group, open new inverted indices and discard those no longer needed ...
        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            // ... and determine position information shared by all amplicons of this length
            for (lenSeqs_t lenDiff = 0; lenDiff < t + 1; lenDiff++) {
                for (lenSeqs_t segmentIndex = 0; segmentIndex < t + k; segmentIndex++) {
                    substrs[lenDiff][segmentIndex] = selectSubstrs(seqLen, seqLen - lenDiff, segmentIndex, t, k);
                }
            }
            selectSegments(segments, seqLen, t, k);

        }


        for (lenSeqs_t len = (seqLen > t) * (seqLen - t); len <= seqLen; len++) { // ... consider already indexed seqs with a feasible length...

            candCnts.clear();

            for (lenSeqs_t i = 0; i < t + k; i++) { // ... and apply segment filter for each segment

                Substrings& subs = substrs[seqLen - len][i];
                InvertedIndex& inv = indices.getIndex(len, i);
                sip.first = ac[curIntId].seq + subs.first;
                sip.second = sip.first + subs.len;

                for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++, sip.first++, sip.second++) {
                    inv.addLabelCountsOf(sip, candCnts);
                }

            }

            // general pigeonhole principle: for being a candidate, at least k segments have to be matched
            for (auto candIter = candCnts.begin(); candIter != candCnts.end(); candIter++) {

                if (candIter->second >= k) {

                    dist = useScore ?
                             Verification::computeGotohLengthAwareEarlyRow8(ac[curIntId].seq, ac[curIntId].len, ac[candIter->first].seq, ac[candIter->first].len, t, scoring, D, P, cntDiffs, cntDiffsP)
                           : Verification::computeLengthAwareRow(ac[curIntId].seq, ac[curIntId].len, ac[candIter->first].seq, ac[candIter->first].len, t, M);

                    if (dist <= t){
                        matches.add(curIntId, candIter->first, dist);
                    }

                }

            }

        }

        // index sequence
        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

    }

}



void SegmentFilter::filterBackward(const AmpliconCollection& ac, const Subpool& sp, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k) {

    RollingIndices<InvertedIndex> indices(t + 1, t + k, false);

    Substrings substrs[t + 1][t + k];
    Segments segments(t + k);
    StringIteratorPair sip;

    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    std::vector<Candidate> candColl;

    lenSeqs_t seqLen = 0;
    numSeqs_t curIntId = sp.end - 1;

    // process index-only amplicons
    for (; curIntId >= sp.beginIndex; curIntId--) {

        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            selectSegments(segments, seqLen, t, k);

        }

        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

    }



    // process remaining amplicons by matching and indexing
    seqLen = 0; // ensures that the substring information are computed
    curIntId++; // courtesy of do-while loop and unsigned integer type
    do { // for every remaining amplicon in the collection...

        curIntId--;

        // on reaching new length group, open new inverted indices and discard those no longer needed ...
        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            // ... and determine position information shared by all amplicons of this length
            for (lenSeqs_t lenDiff = 0; lenDiff < t + 1; lenDiff++) {
                for (lenSeqs_t segmentIndex = 0; segmentIndex < t + k; segmentIndex++) {
                    substrs[lenDiff][segmentIndex] = selectSubstrsBackward(seqLen, seqLen + lenDiff, segmentIndex, t, k);
                }
            }
            selectSegments(segments, seqLen, t, k);

        }


        for (lenSeqs_t len = seqLen + t; len >= seqLen; len--) { // ... consider already indexed seqs with a feasible length...

            candCnts.clear();

            for (lenSeqs_t i = 0; i < t + k; i++) { // ... and apply segment filter for each segment

                Substrings& subs = substrs[len - seqLen][i];
                InvertedIndex& inv = indices.getIndex(len, i);
                sip.first = ac[curIntId].seq + subs.first;
                sip.second = sip.first + subs.len;

                for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++, sip.first++, sip.second++) {
                    inv.addLabelCountsOf(sip, candCnts);
                }

            }

            // general pigeonhole principle: for being a candidate, at least k segments have to be matched
            for (auto candIter = candCnts.begin(); candIter != candCnts.end(); candIter++) {

                if (candIter->second >= k) {
                    candColl.push_back(Candidate(curIntId, candIter->first));
                }

            }

        }

        // index sequence
        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

        cands.push(candColl);
        candColl = std::vector<Candidate>();

    } while (curIntId != sp.beginMatch);

}


void SegmentFilter::filterBackwardDirectly(const AmpliconCollection& ac, const Subpool& sp, Matches& matches, const lenSeqs_t t, const lenSeqs_t k, const bool useScore, const Verification::Scoring& scoring) {

    RollingIndices<InvertedIndex> indices(t + 1, t + k, false);

    Substrings substrs[t + 1][t + k];
    Segments segments(t + k);
    StringIteratorPair sip;

    lenSeqs_t M[useScore ? 1 : ac.back().len + 1];
    val_t D[useScore? ac.back().len + 1 : 1];
    val_t P[useScore? ac.back().len + 1 : 1];
    lenSeqs_t cntDiffs[useScore? ac.back().len + 1 : 1];
    lenSeqs_t cntDiffsP[useScore? ac.back().len + 1 : 1];

    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;

    lenSeqs_t seqLen = 0;
    numSeqs_t curIntId = sp.end - 1;
    lenSeqs_t dist;

    // process index-only amplicons
    for (; curIntId >= sp.beginIndex; curIntId--) {

        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            selectSegments(segments, seqLen, t, k);

        }

        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

    }



    // process remaining amplicons by matching and indexing
    seqLen = 0; // ensures that the substring information are computed
    curIntId++; // courtesy of do-while loop and unsigned integer type
    do { // for every remaining amplicon in the collection...

        curIntId--;

        // on reaching new length group, open new inverted indices and discard those no longer needed ...
        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            // ... and determine position information shared by all amplicons of this length
            for (lenSeqs_t lenDiff = 0; lenDiff < t + 1; lenDiff++) {
                for (lenSeqs_t segmentIndex = 0; segmentIndex < t + k; segmentIndex++) {
                    substrs[lenDiff][segmentIndex] = selectSubstrsBackward(seqLen, seqLen + lenDiff, segmentIndex, t, k);
                }
            }
            selectSegments(segments, seqLen, t, k);

        }


        for (lenSeqs_t len = seqLen + t; len >= seqLen; len--) { // ... consider already indexed seqs with a feasible length...

            candCnts.clear();

            for (lenSeqs_t i = 0; i < t + k; i++) { // ... and apply segment filter for each segment

                Substrings& subs = substrs[len - seqLen][i];
                InvertedIndex& inv = indices.getIndex(len, i);
                sip.first = ac[curIntId].seq + subs.first;
                sip.second = sip.first + subs.len;

                for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++, sip.first++, sip.second++) {
                    inv.addLabelCountsOf(sip, candCnts);
                }

            }

            // general pigeonhole principle: for being a candidate, at least k segments have to be matched
            for (auto candIter = candCnts.begin(); candIter != candCnts.end(); candIter++) {

                if (candIter->second >= k) {

                    dist = useScore ?
                             Verification::computeGotohLengthAwareEarlyRow8(ac[curIntId].seq, ac[curIntId].len, ac[candIter->first].seq, ac[candIter->first].len, t, scoring, D, P, cntDiffs, cntDiffsP)
                           : Verification::computeLengthAwareRow(ac[curIntId].seq, ac[curIntId].len, ac[candIter->first].seq, ac[candIter->first].len, t, M);

                    if (dist <= t) {
                        matches.add(curIntId, candIter->first, dist);
                    }
                }

            }

        }

        // index sequence
        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

    } while (curIntId != sp.beginMatch);

}



void SegmentFilter::filterForwardBackward(const AmpliconCollection& ac, const Subpool& sp, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k) {

    RollingIndices<InvertedIndex> indices(t + 1, t + k, true);

    Substrings substrs[t + 1][t + k];
    Segments segments(t + k);
    std::vector<std::string> segmentStrs(t + k);
    StringIteratorPair sip;

    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    std::vector<Candidate> candColl;

    lenSeqs_t seqLen = 0;
    numSeqs_t curIntId = sp.beginIndex;

    // process index-only amplicons
    for (; curIntId < sp.beginMatch; curIntId++) {

        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            selectSegments(segments, seqLen, t, k);

        }

        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

    }



    // process remaining amplicons by matching and indexing
    seqLen = 0; // ensures that the substring information are computed
    for (; curIntId < sp.end; curIntId++) { // for every amplicon in the collection...

        // on reaching new length group, open new inverted indices and discard those no longer needed ...
        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            // ... and determine position information shared by all amplicons of this length
            for (lenSeqs_t lenDiff = 0; lenDiff < t + 1; lenDiff++) {
                for (lenSeqs_t segmentIndex = 0; segmentIndex < t + k; segmentIndex++) {
                    substrs[lenDiff][segmentIndex] = selectSubstrs(seqLen, seqLen - lenDiff, segmentIndex, t, k);
                }
            }
            selectSegments(segments, seqLen, t, k);

        }

        // compute actual segments of the current amplicon here once (used in pipelined filtering and indexing steps)
        for (auto i = 0; i < t + k; i++) {
            segmentStrs[i] = std::string(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second);
        }

        for (lenSeqs_t len = (seqLen > t) * (seqLen - t); len <= seqLen; len++) { // ... consider already indexed seqs with a feasible length...

            candCnts.clear();

            for (lenSeqs_t i = 0; i < t + k; i++) { // ... and apply segment filter for each segment

                Substrings& subs = substrs[seqLen - len][i];
                InvertedIndex& inv = indices.getIndex(len, i);
                sip.first = ac[curIntId].seq + subs.first;
                sip.second = sip.first + subs.len;

                for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++, sip.first++, sip.second++) {
                    inv.addLabelCountsOf(sip, candCnts);
                }

            }


            // general pigeonhole principle: for being a candidate, at least k segments have to be matched
            //  + pipelined backward filtering
            //      advantages of performing it directly here:
            //      (a) initial candidate set stays small (more precisely, there is no actual intermediate candidate set between the filtering steps)
            //      (b) length of candidates is known and equal for all of them -> selectSubstrs() needed only once
            //      (c) segment information of current amplicon known -> reuse
            std::string candStr;
            lenSeqs_t cnt = 0; // number of substring-segment matches for the current candidate in the second filter step

            Substrings candSubs[t + k]; // common substring information of all candidates
            for (lenSeqs_t i = 0; i < t + k; i++) {
                candSubs[i] = selectSubstrsBackward(len, seqLen, i, t, k);
            }

            for (auto candIter = candCnts.begin(); candIter != candCnts.end(); candIter++) { // 1st component = integer id of amplicon, 2nd component = number of substring-segment matches during first filter step

                if (candIter->second >= k) { // initial candidate found -> immediate backward filtering

                    candStr = ac[candIter->first].seq;

                    for (lenSeqs_t i = 0; i < t + k && cnt < k; i++) {
                        cnt += (candStr.substr(candSubs[i].first, (candSubs[i].last - candSubs[i].first) + candSubs[i].len).find(segmentStrs[i]) < std::string::npos);
                    }

                    if (cnt == k) {
                        candColl.push_back(Candidate(curIntId, candIter->first));
                    }

                    cnt = 0;


                }

            }

        }

        // index sequence
        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

        cands.push(candColl);
        candColl = std::vector<Candidate>();

    }

}


void SegmentFilter::filterForwardBackwardDirectly(const AmpliconCollection& ac, const Subpool& sp, Matches& matches, const lenSeqs_t t, const lenSeqs_t k, const bool useScore, const Verification::Scoring& scoring) {

    RollingIndices<InvertedIndex> indices(t + 1, t + k, true);

    Substrings substrs[t + 1][t + k];
    Segments segments(t + k);
    std::vector<std::string> segmentStrs(t + k);
    StringIteratorPair sip;

    lenSeqs_t M[useScore ? 1 : ac.back().len + 1];
    val_t D[useScore? ac.back().len + 1 : 1];
    val_t P[useScore? ac.back().len + 1 : 1];
    lenSeqs_t cntDiffs[useScore? ac.back().len + 1 : 1];
    lenSeqs_t cntDiffsP[useScore? ac.back().len + 1 : 1];

    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;

    lenSeqs_t seqLen = 0;
    numSeqs_t curIntId = sp.beginIndex;
    lenSeqs_t dist;

    // process index-only amplicons
    for (; curIntId < sp.beginMatch; curIntId++) {

        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            selectSegments(segments, seqLen, t, k);

        }

        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

    }



    // process remaining amplicons by matching and indexing
    seqLen = 0; // ensures that the substring information are computed
    for (; curIntId < sp.end; curIntId++) { // for every amplicon in the collection...

        // on reaching new length group, open new inverted indices and discard those no longer needed ...
        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            // ... and determine position information shared by all amplicons of this length
            for (lenSeqs_t lenDiff = 0; lenDiff < t + 1; lenDiff++) {
                for (lenSeqs_t segmentIndex = 0; segmentIndex < t + k; segmentIndex++) {
                    substrs[lenDiff][segmentIndex] = selectSubstrs(seqLen, seqLen - lenDiff, segmentIndex, t, k);
                }
            }
            selectSegments(segments, seqLen, t, k);

        }

        // compute actual segments of the current amplicon here once (used in pipelined filtering and indexing steps)
        for (auto i = 0; i < t + k; i++) {
            segmentStrs[i] = std::string(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second);
        }

        for (lenSeqs_t len = (seqLen > t) * (seqLen - t); len <= seqLen; len++) { // ... consider already indexed seqs with a feasible length...

            candCnts.clear();

            for (lenSeqs_t i = 0; i < t + k; i++) { // ... and apply segment filter for each segment

                Substrings& subs = substrs[seqLen - len][i];
                InvertedIndex& inv = indices.getIndex(len, i);
                sip.first = ac[curIntId].seq + subs.first;
                sip.second = sip.first + subs.len;

                for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++, sip.first++, sip.second++) {
                    inv.addLabelCountsOf(sip, candCnts);
                }

            }


            // general pigeonhole principle: for being a candidate, at least k segments have to be matched
            //  + pipelined backward filtering
            //      advantages of performing it directly here:
            //      (a) initial candidate set stays small (more precisely, there is no actual intermediate candidate set between the filtering steps)
            //      (b) length of candidates is known and equal for all of them -> selectSubstrs() needed only once
            //      (c) segment information of current amplicon known -> reuse
            std::string candStr;
            lenSeqs_t cnt = 0; // number of substring-segment matches for the current candidate in the second filter step

            Substrings candSubs[t + k]; // common substring information of all candidates
            for (lenSeqs_t i = 0; i < t + k; i++) {
                candSubs[i] = selectSubstrsBackward(len, seqLen, i, t, k);
            }

            for (auto candIter = candCnts.begin(); candIter != candCnts.end(); candIter++) { // 1st component = integer id of amplicon, 2nd component = number of substring-segment matches during first filter step

                if (candIter->second >= k) { // initial candidate found -> immediate backward filtering

                    candStr = ac[candIter->first].seq;

                    for (lenSeqs_t i = 0; i < t + k && cnt < k; i++) {
                        cnt += (candStr.substr(candSubs[i].first, (candSubs[i].last - candSubs[i].first) + candSubs[i].len).find(segmentStrs[i]) < std::string::npos);
                    }

                    if (cnt == k) {

                        dist = useScore ?
                                 Verification::computeGotohLengthAwareEarlyRow8(ac[curIntId].seq, ac[curIntId].len, ac[candIter->first].seq, ac[candIter->first].len, t, scoring, D, P, cntDiffs, cntDiffsP)
                               : Verification::computeLengthAwareRow(ac[curIntId].seq, ac[curIntId].len, ac[candIter->first].seq, ac[candIter->first].len, t, M);

                        if (dist <= t){
                            matches.add(curIntId, candIter->first, dist);
                        }

                    }

                    cnt = 0;


                }

            }

        }

        // index sequence
        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

    }

}



void SegmentFilter::filterBackwardForward(const AmpliconCollection& ac, const Subpool& sp, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k) {

    RollingIndices<InvertedIndex> indices(t + 1, t + k, false);

    Substrings substrs[t + 1][t + k];
    Segments segments(t + k);
    std::vector<std::string> segmentStrs(t + k);
    StringIteratorPair sip;

    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    std::vector<Candidate> candColl;

    lenSeqs_t seqLen = 0;
    numSeqs_t curIntId = sp.end - 1;

    // process index-only amplicons
    for (; curIntId >= sp.beginIndex; curIntId--) {

        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            selectSegments(segments, seqLen, t, k);

        }

        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

    }



    // process remaining amplicons by matching and indexing
    seqLen = 0; // ensures that the substring information are computed
    curIntId++; // courtesy of do-while loop and unsigned integer type
    do { // for every remaining amplicon in the collection...

        curIntId--;

        // on reaching new length group, open new inverted indices and discard those no longer needed ...
        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            // ... and determine position information shared by all amplicons of this length
            for (lenSeqs_t lenDiff = 0; lenDiff < t + 1; lenDiff++) {
                for (lenSeqs_t segmentIndex = 0; segmentIndex < t + k; segmentIndex++) {
                    substrs[lenDiff][segmentIndex] = selectSubstrsBackward(seqLen, seqLen + lenDiff, segmentIndex, t, k);
                }
            }
            selectSegments(segments, seqLen, t, k);

        }

        // compute actual segments of the current amplicon here once (used in pipelined filtering and indexing steps)
        for (auto i = 0; i < t + k; i++) {
            segmentStrs[i] = std::string(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second);
        }

        for (lenSeqs_t len = seqLen + t; len >= seqLen; len--) { // ... consider already indexed seqs with a feasible length...

            candCnts.clear();

            for (lenSeqs_t i = 0; i < t + k; i++) { // ... and apply segment filter for each segment

                Substrings& subs = substrs[len - seqLen][i];
                InvertedIndex& inv = indices.getIndex(len, i);
                sip.first = ac[curIntId].seq + subs.first;
                sip.second = sip.first + subs.len;

                for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++, sip.first++, sip.second++) {
                    inv.addLabelCountsOf(sip, candCnts);
                }

            }


            // general pigeonhole principle: for being a candidate, at least k segments have to be matched
            //  + pipelined forward filtering
            //      advantages of performing it directly here:
            //      (a) initial candidate set stays small (more precisely, there is no actual intermediate candidate set between the filtering steps)
            //      (b) length of candidates is known and equal for all of them -> selectSubstrs() needed only once
            //      (c) segment information of current amplicon known -> reuse
            std::string candStr;
            lenSeqs_t cnt = 0; // number of substring-segment matches for the current candidate in the second filter step

            Substrings candSubs[t + k]; // common substring information of all candidates
            for (lenSeqs_t i = 0; i < t + k; i++) {
                candSubs[i] = selectSubstrs(len, seqLen, i, t, k);
            }

            for (auto candIter = candCnts.begin(); candIter != candCnts.end(); candIter++) { // 1st component = integer id of amplicon, 2nd component = number of substring-segment matches during first filter step

                if (candIter->second >= k) { // initial candidate found -> immediate forward filtering

                    candStr = ac[candIter->first].seq;

                    for (lenSeqs_t i = 0; i < t + k && cnt < k; i++) {
                        cnt += (candStr.substr(candSubs[i].first, (candSubs[i].last - candSubs[i].first) + candSubs[i].len).find(segmentStrs[i]) < std::string::npos);
                    }

                    if (cnt == k) {
                        candColl.push_back(Candidate(curIntId, candIter->first));
                    }

                    cnt = 0;


                }

            }

        }

        // index sequence
        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

        cands.push(candColl);
        candColl = std::vector<Candidate>();

    } while (curIntId != sp.beginMatch);

}


void SegmentFilter::filterBackwardForwardDirectly(const AmpliconCollection& ac, const Subpool& sp, Matches& matches, const lenSeqs_t t, const lenSeqs_t k, const bool useScore, const Verification::Scoring& scoring) {

    RollingIndices<InvertedIndex> indices(t + 1, t + k, false);

    Substrings substrs[t + 1][t + k];
    Segments segments(t + k);
    std::vector<std::string> segmentStrs(t + k);
    StringIteratorPair sip;

    lenSeqs_t M[useScore ? 1 : ac.back().len + 1];
    val_t D[useScore? ac.back().len + 1 : 1];
    val_t P[useScore? ac.back().len + 1 : 1];
    lenSeqs_t cntDiffs[useScore? ac.back().len + 1 : 1];
    lenSeqs_t cntDiffsP[useScore? ac.back().len + 1 : 1];

    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;

    lenSeqs_t seqLen = 0;
    numSeqs_t curIntId = sp.end - 1;
    lenSeqs_t dist;

    // process index-only amplicons
    for (; curIntId >= sp.beginIndex; curIntId--) {

        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            selectSegments(segments, seqLen, t, k);

        }

        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

    }



    // process remaining amplicons by matching and indexing
    seqLen = 0; // ensures that the substring information are computed
    curIntId++; // courtesy of do-while loop and unsigned integer type
    do { // for every remaining amplicon in the collection...

        curIntId--;

        // on reaching new length group, open new inverted indices and discard those no longer needed ...
        if (ac[curIntId].len != seqLen) {

            seqLen = ac[curIntId].len;
            indices.roll(seqLen);

            // ... and determine position information shared by all amplicons of this length
            for (lenSeqs_t lenDiff = 0; lenDiff < t + 1; lenDiff++) {
                for (lenSeqs_t segmentIndex = 0; segmentIndex < t + k; segmentIndex++) {
                    substrs[lenDiff][segmentIndex] = selectSubstrsBackward(seqLen, seqLen + lenDiff, segmentIndex, t, k);
                }
            }
            selectSegments(segments, seqLen, t, k);

        }

        // compute actual segments of the current amplicon here once (used in pipelined filtering and indexing steps)
        for (auto i = 0; i < t + k; i++) {
            segmentStrs[i] = std::string(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second);
        }

        for (lenSeqs_t len = seqLen + t; len >= seqLen; len--) { // ... consider already indexed seqs with a feasible length...

            candCnts.clear();

            for (lenSeqs_t i = 0; i < t + k; i++) { // ... and apply segment filter for each segment

                Substrings& subs = substrs[len - seqLen][i];
                InvertedIndex& inv = indices.getIndex(len, i);
                sip.first = ac[curIntId].seq + subs.first;
                sip.second = sip.first + subs.len;

                for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++, sip.first++, sip.second++) {
                    inv.addLabelCountsOf(sip, candCnts);
                }

            }


            // general pigeonhole principle: for being a candidate, at least k segments have to be matched
            //  + pipelined forward filtering
            //      advantages of performing it directly here:
            //      (a) initial candidate set stays small (more precisely, there is no actual intermediate candidate set between the filtering steps)
            //      (b) length of candidates is known and equal for all of them -> selectSubstrs() needed only once
            //      (c) segment information of current amplicon known -> reuse
            std::string candStr;
            lenSeqs_t cnt = 0; // number of substring-segment matches for the current candidate in the second filter step

            Substrings candSubs[t + k]; // common substring information of all candidates
            for (lenSeqs_t i = 0; i < t + k; i++) {
                candSubs[i] = selectSubstrs(len, seqLen, i, t, k);
            }

            for (auto candIter = candCnts.begin(); candIter != candCnts.end(); candIter++) { // 1st component = integer id of amplicon, 2nd component = number of substring-segment matches during first filter step

                if (candIter->second >= k) { // initial candidate found -> immediate forward filtering

                    candStr = ac[candIter->first].seq;

                    for (lenSeqs_t i = 0; i < t + k && cnt < k; i++) {
                        cnt += (candStr.substr(candSubs[i].first, (candSubs[i].last - candSubs[i].first) + candSubs[i].len).find(segmentStrs[i]) < std::string::npos);
                    }

                    if (cnt == k) {

                        dist = useScore ?
                                 Verification::computeGotohLengthAwareEarlyRow8(ac[curIntId].seq, ac[curIntId].len, ac[candIter->first].seq, ac[candIter->first].len, t, scoring, D, P, cntDiffs, cntDiffsP)
                               : Verification::computeLengthAwareRow(ac[curIntId].seq, ac[curIntId].len, ac[candIter->first].seq, ac[candIter->first].len, t, M);

                        if (dist <= t){
                            matches.add(curIntId, candIter->first, dist);
                        }

                    }

                    cnt = 0;


                }

            }

        }

        // index sequence
        for (lenSeqs_t i = 0; i < t + k; i++) {
            indices.getIndex(seqLen,i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

    } while (curIntId != sp.beginMatch);

}


void SegmentFilter::filter(const AmpliconCollection& ac, const Subpool& sp, RotatingBuffers<Candidate>& cands, const lenSeqs_t t, const lenSeqs_t k, const int mode) { //TODO set "best" default (see also config creation)

    switch (mode) {

        case 0:
            filterForward(ac, sp, cands, t, k);
            break;

        case 1:
            filterBackward(ac, sp, cands, t, k);
            break;

        case 2:
            filterForwardBackward(ac, sp, cands, t, k);
            break;

        case 3:
            filterBackwardForward(ac, sp, cands, t, k);
            break;

        default:
            filterForward(ac, sp, cands, t, k);

    }

}

void SegmentFilter::filterDirectly(const AmpliconCollection& ac, const Subpool& sp, Matches& matches, const lenSeqs_t t, const lenSeqs_t k, const int mode, const bool useScore, const Verification::Scoring& scoring) { //TODO set "best" default (see also config creation)

    switch (mode) {

        case 0:
            filterForwardDirectly(ac, sp, matches, t, k, useScore, scoring);
            break;

        case 1:
            filterBackwardDirectly(ac, sp, matches, t, k, useScore, scoring);
            break;

        case 2:
            filterForwardBackwardDirectly(ac, sp, matches, t, k, useScore, scoring);
            break;

        case 3:
            filterBackwardForwardDirectly(ac, sp, matches, t, k, useScore, scoring);
            break;

        default:
            filterForwardDirectly(ac, sp, matches, t, k, useScore, scoring);

    }

}



void SegmentFilter::writeMatches(std::string oFile, AmpliconPools& pools, std::vector<Matches*>& allMatches) {

    std::ofstream oStream(oFile);

    for (lenSeqs_t i = 0; i < pools.numPools(); i++) {

        AmpliconCollection* ac = pools.get(i);
        Matches* matches = allMatches[i];

        for (numSeqs_t j = 0; j < ac->size(); j++) {

            auto m = matches->getMatchesOfAugmented(j);

            for (auto matchIter = m.begin(); matchIter != m.end(); matchIter++) {
                oStream << (*ac)[j].id << ";" << (*ac)[matchIter->first].id << ";" << (matchIter->second) << std::endl;
            }

        }

    }

    oStream.close();

}


void SegmentFilter::writeMatchesOneWay(std::string oFile, AmpliconPools& pools, std::vector<Matches*>& allMatches) {

    std::ofstream oStream(oFile);

    for (lenSeqs_t i = 0; i < pools.numPools(); i++) {

        AmpliconCollection* ac = pools.get(i);
        Matches* matches = allMatches[i];

        for (numSeqs_t j = 0; j < ac->size(); j++) {

            auto m = matches->getMatchesOfAugmented(j);

            for (auto matchIter = m.begin(); matchIter != m.end(); matchIter++) {

                if (matchIter->first <= j) continue; // to avoid writing a match in both directions; a match is written with the shorter (or, if equally long, lexicographically 'smaller') amplicon first
                oStream << (*ac)[j].id << ";" << (*ac)[matchIter->first].id << ";" << (matchIter->second) << std::endl;
            }

        }

    }

    oStream.close();

}

}