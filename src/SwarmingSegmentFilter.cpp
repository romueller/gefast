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

#include "../include/SIMD.hpp"
#include "../include/SwarmingSegmentFilter.hpp"


namespace GeFaST {

void SegmentFilter::addCandCnts(const Amplicon& amplicon, lenSeqs_t childLen, lenSeqs_t numSegments, std::unordered_map<numSeqs_t, lenSeqs_t>& candCnts, SwarmingIndices& indices, std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive) {

    StringIteratorPair sip;
    for (lenSeqs_t i = 0; i < numSegments; i++) {

        const Substrings& subs = substrsArchive[amplicon.len][childLen][i];
        SwarmingInvertedIndex& inv = indices.getIndex(childLen, i);
        sip.first = amplicon.seq + subs.first;
        sip.second = sip.first + subs.len;

        for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++, sip.first++, sip.second++) {
            inv.addLabelCountsOf(sip, candCnts);
        }

    }

}

void SegmentFilter::verifyCands(const Amplicon& amplicon, const AmpliconCollection& ac, std::unordered_map<numSeqs_t, lenSeqs_t>& candCnts, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& matches, const SwarmClustering::SwarmConfig& sc, lenSeqs_t* M, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP) {

    // general pigeonhole principle: for being a candidate, at least k segments have to be matched
    for (auto candIter = candCnts.begin(); candIter != candCnts.end(); candIter++) {

#if QGRAM_FILTER
        if ((candIter->second >= sc.extraSegs) && (sc.noOtuBreaking || amplicon.abundance >= ac[candIter->first].abundance) && (qgram_diff(amplicon, ac[candIter->first]) <= sc.threshold)) {
#else
            if ((candIter->second >= sc.extraSegs) && (sc.noOtuBreaking || amplicon.abundance >= ac[candIter->first].abundance)) {
#endif

            lenSeqs_t dist = sc.useScore ?
                                  Verification::computeGotohLengthAwareEarlyRow(amplicon.seq, amplicon.len, ac[candIter->first].seq, ac[candIter->first].len, sc.threshold, sc.scoring, D, P, cntDiffs, cntDiffsP)
                                : Verification::computeLengthAwareRow(amplicon.seq, amplicon.len, ac[candIter->first].seq, ac[candIter->first].len, sc.threshold, M);

            if (dist <= sc.threshold) {
                matches.push_back(std::make_pair(candIter->first, dist));
            }

        }

    }

}

void SegmentFilter::verifyCandsTwoWay(const Amplicon& amplicon, std::vector<std::string>& segmentStrs, const AmpliconCollection& ac, std::unordered_map<numSeqs_t, lenSeqs_t>& candCnts, std::vector<Substrings>& candSubstrs, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& matches, const SwarmClustering::SwarmConfig& sc, lenSeqs_t* M, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP) {

    std::string candStr;
    lenSeqs_t cnt = 0; // number of substring-segment matches for the current candidate in the second filter step

    // general pigeonhole principle: for being a candidate, at least k segments have to be matched
    //  + pipelined backward filtering
    for (auto candIter = candCnts.begin(); candIter != candCnts.end(); candIter++) {

        if ((candIter->second >= sc.extraSegs) && (sc.noOtuBreaking || amplicon.abundance >= ac[candIter->first].abundance)) {

            candStr = ac[candIter->first].seq;

            for (lenSeqs_t i = 0; i < sc.threshold + sc.extraSegs && cnt < sc.extraSegs; i++) {
                cnt += (candStr.substr(candSubstrs[i].first, (candSubstrs[i].last - candSubstrs[i].first) + candSubstrs[i].len).find(segmentStrs[i]) < std::string::npos);
            }

#if QGRAM_FILTER
            if ((cnt == sc.extraSegs) && (qgram_diff(amplicon, ac[candIter->first]) <= sc.threshold)) {
#else
                if (cnt == sc.extraSegs) {
#endif

                lenSeqs_t dist = sc.useScore ?
                                      Verification::computeGotohLengthAwareEarlyRow(amplicon.seq, amplicon.len, ac[candIter->first].seq, ac[candIter->first].len, sc.threshold, sc.scoring, D, P, cntDiffs, cntDiffsP)
                                    : Verification::computeLengthAwareRow(amplicon.seq, amplicon.len, ac[candIter->first].seq, ac[candIter->first].len, sc.threshold, M);

                if (dist <= sc.threshold) {
                    matches.push_back(std::make_pair(candIter->first, dist));
                }

            }

            cnt = 0;

        }

    }

}


SegmentFilter::ChildrenFinder::ChildrenFinder(const AmpliconCollection& ac, SwarmingIndices& indices, std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive, const SwarmClustering::SwarmConfig& sc, lenSeqs_t* M, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP) : ac_(ac), indices_(indices), substrsArchive_(substrsArchive), sc_(sc) {

    M_ = M;
    D_ = D;
    P_ = P;
    cntDiffs_ = cntDiffs;
    cntDiffsP_ = cntDiffsP;

}

std::vector<std::pair<numSeqs_t, lenSeqs_t>> SegmentFilter::ChildrenFinder::getChildren(const numSeqs_t id) {

    std::vector<std::pair<numSeqs_t, lenSeqs_t>> matches;
    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    auto& amplicon = ac_[id];

    for (lenSeqs_t childLen = (amplicon.len > sc_.threshold) * (amplicon.len - sc_.threshold); childLen <= amplicon.len + sc_.threshold; childLen++) {

        candCnts.clear();

        SegmentFilter::addCandCnts(amplicon, childLen, sc_.threshold + sc_.extraSegs, candCnts, indices_, substrsArchive_);
        SegmentFilter::verifyCands(amplicon, ac_, candCnts, matches, sc_, M_, D_, P_, cntDiffs_, cntDiffsP_);

    }

    return matches;

}

void SegmentFilter::ChildrenFinder::getChildren(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children) {

    children.clear();

    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    auto& amplicon = ac_[id];

    for (lenSeqs_t childLen = (amplicon.len > sc_.threshold) * (amplicon.len - sc_.threshold); childLen <= amplicon.len + sc_.threshold; childLen++) {

        candCnts.clear();

        SegmentFilter::addCandCnts(amplicon, childLen, sc_.threshold + sc_.extraSegs, candCnts, indices_, substrsArchive_);
        SegmentFilter::verifyCands(amplicon, ac_, candCnts, children,sc_, M_, D_, P_, cntDiffs_, cntDiffsP_);

    }

}

std::vector<std::pair<numSeqs_t, lenSeqs_t>> SegmentFilter::ChildrenFinder::getChildrenTwoWay(const numSeqs_t id) {

    std::vector<std::pair<numSeqs_t, lenSeqs_t>> matches;
    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    auto& amplicon = ac_[id];

    Segments segments(sc_.threshold + sc_.extraSegs);
    selectSegments(segments, amplicon.len, sc_.threshold, sc_.extraSegs);

    std::vector<std::string> segmentStrs(sc_.threshold + sc_.extraSegs);
    for (auto i = 0; i < sc_.threshold + sc_.extraSegs; i++) {
        segmentStrs[i] = std::string(amplicon.seq + segments[i].first, amplicon.seq + segments[i].first + segments[i].second);
    }

    for (lenSeqs_t childLen = (amplicon.len > sc_.threshold) * (amplicon.len - sc_.threshold); childLen <= amplicon.len + sc_.threshold; childLen++) {

        candCnts.clear();

        SegmentFilter::addCandCnts(amplicon, childLen, sc_.threshold + sc_.extraSegs, candCnts, indices_, substrsArchive_);
        SegmentFilter::verifyCandsTwoWay(amplicon, segmentStrs, ac_, candCnts, substrsArchive_[childLen][amplicon.len], matches, sc_, M_, D_, P_, cntDiffs_, cntDiffsP_);

    }

    return matches;

}

void SegmentFilter::ChildrenFinder::getChildrenTwoWay(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children) {

    children.clear();

    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    auto& amplicon = ac_[id];

    Segments segments(sc_.threshold + sc_.extraSegs);
    selectSegments(segments, amplicon.len, sc_.threshold, sc_.extraSegs);

    std::vector<std::string> segmentStrs(sc_.threshold + sc_.extraSegs);
    for (auto i = 0; i < sc_.threshold + sc_.extraSegs; i++) {
        segmentStrs[i] = std::string(amplicon.seq + segments[i].first, amplicon.seq + segments[i].first + segments[i].second);
    }

    for (lenSeqs_t childLen = (amplicon.len > sc_.threshold) * (amplicon.len - sc_.threshold); childLen <= amplicon.len + sc_.threshold; childLen++) {

        candCnts.clear();

        SegmentFilter::addCandCnts(amplicon, childLen, sc_.threshold + sc_.extraSegs, candCnts, indices_, substrsArchive_);
        SegmentFilter::verifyCandsTwoWay(amplicon, segmentStrs, ac_, candCnts, substrsArchive_[childLen][amplicon.len], children, sc_, M_, D_, P_, cntDiffs_, cntDiffsP_);

    }

}


SegmentFilter::ParallelChildrenFinder::ParallelChildrenFinder(const AmpliconCollection& ac, SwarmingIndices& indices, std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive, const lenSeqs_t width, const SwarmClustering::SwarmConfig& sc) : ac_(ac), indices_(indices), substrsArchive_(substrsArchive), sc_(sc) {

    cbs_ = RotatingBuffers<Candidate>(sc.numThreadsPerCheck - 1);

    numDiscarded_ = 0;
    verifierThreads_ = std::vector<std::thread>(sc.numThreadsPerCheck - 1);
    auto fun = sc.useScore ? &ParallelChildrenFinder::verifyGotoh : &ParallelChildrenFinder::verify;

    for (unsigned long v = 0; v < sc.numThreadsPerCheck - 1; v++) {
        verifierThreads_[v] = std::thread(fun, this, std::ref(matches_), std::ref(cbs_.getBuffer(v)), width);
    }

}

SegmentFilter::ParallelChildrenFinder::~ParallelChildrenFinder() {

    cbs_.close();

    for (unsigned long v = 0; v < sc_.numThreadsPerCheck - 1; v++) {
        verifierThreads_[v].join();
    }

}

numSeqs_t SegmentFilter::ParallelChildrenFinder::sendCandsToVerification(const numSeqs_t id, const Amplicon& amplicon, std::unordered_map<numSeqs_t, lenSeqs_t>& candCnts) {

    numSeqs_t numCands = 0;

    // general pigeonhole principle: for being a candidate, at least k segments have to be matched
    for (auto candIter = candCnts.begin(); candIter != candCnts.end(); candIter++) {

#if QGRAM_FILTER
        if ((candIter->second >= sc_.extraSegs) && (sc_.noOtuBreaking || amplicon.abundance >= ac_[candIter->first].abundance) && (qgram_diff(amplicon, ac_[candIter->first]) <= sc_.threshold)) {
#else
            if ((candIter->second >= sc_.extraSegs) && (sc_.noOtuBreaking || amplicon.abundance >= ac_[candIter->first].abundance)) {
#endif

            Candidate cand(id, candIter->first);
            cbs_.push(cand);
            numCands++;

        }

    }

    return numCands;

}

numSeqs_t SegmentFilter::ParallelChildrenFinder::sendCandsToVerificationTwoWay(const numSeqs_t id, const Amplicon& amplicon, std::vector<std::string>& segmentStrs, std::unordered_map<numSeqs_t, lenSeqs_t>& candCnts, std::vector<Substrings>& candSubstrs) {

    std::string candStr;
    lenSeqs_t cnt = 0; // number of substring-segment matches for the current candidate in the second filter step
    numSeqs_t numCands = 0;

    // general pigeonhole principle: for being a candidate, at least k segments have to be matched
    //  + pipelined backward filtering
    for (auto candIter = candCnts.begin(); candIter != candCnts.end(); candIter++) {

        if ((candIter->second >= sc_.extraSegs) && (sc_.noOtuBreaking || amplicon.abundance >= ac_[candIter->first].abundance)) {

            candStr = ac_[candIter->first].seq;

            for (lenSeqs_t i = 0; i < sc_.threshold + sc_.extraSegs && cnt < sc_.extraSegs; i++) {
                cnt += (candStr.substr(candSubstrs[i].first, (candSubstrs[i].last - candSubstrs[i].first) + candSubstrs[i].len).find(segmentStrs[i]) < std::string::npos);
            }

#if QGRAM_FILTER
            if ((cnt == sc_.extraSegs) && (qgram_diff(amplicon, ac_[candIter->first]) <= sc_.threshold)) {
#else
                if (cnt == sc_.extraSegs) {
#endif

                Candidate cand(id, candIter->first);
                cbs_.push(cand);
                numCands++;

            }

            cnt = 0;

        }

    }

    return numCands;

}

void SegmentFilter::ParallelChildrenFinder::verify(std::vector<std::pair<numSeqs_t, lenSeqs_t>>& matches, Buffer<Candidate>& buf, lenSeqs_t width) {

    Candidate c;
    Buffer<Candidate> localBuffer;
    lenSeqs_t M[width]; // reusable DP-matrix (wide enough for all possible calculations for this AmpliconCollection)

    while (!buf.isClosed() || buf.syncSize() > 0) {

        buf.syncSwapContents(localBuffer);

        while (localBuffer.size() > 0) {

            c = localBuffer.pop();

            lenSeqs_t d = Verification::computeLengthAwareRow(ac_[c.first].seq, ac_[c.first].len, ac_[c.second].seq, ac_[c.second].len, sc_.threshold, M);

            if (d <= sc_.threshold) {

                std::lock_guard<std::mutex> lock(mtxMatches_);
                matches.push_back(std::make_pair(c.second, d));

            } else {

                std::lock_guard<std::mutex> lock(mtxDiscard_);
                numDiscarded_++;

            }

        }

    }

}

void SegmentFilter::ParallelChildrenFinder::verifyGotoh(std::vector<std::pair<numSeqs_t, lenSeqs_t>>& matches, Buffer<Candidate>& buf, lenSeqs_t width) {

    Candidate c;
    Buffer<Candidate> localBuffer;
    val_t D[width]; // reusable DP-matrix (wide enough for all possible calculations for this AmpliconCollection)
    val_t P[width];
    lenSeqs_t cntDiffs[width];
    lenSeqs_t cntDiffsP[width];

    while (!buf.isClosed() || buf.syncSize() > 0) {

        buf.syncSwapContents(localBuffer);

        while (localBuffer.size() > 0) {

            c = localBuffer.pop();

            lenSeqs_t d = Verification::computeGotohLengthAwareEarlyRow(ac_[c.first].seq, ac_[c.first].len, ac_[c.second].seq, ac_[c.second].len, sc_.threshold, sc_.scoring, D, P, cntDiffs, cntDiffsP);

            if (d <= sc_.threshold) {

                std::lock_guard<std::mutex> lock(mtxMatches_);
                matches.push_back(std::make_pair(c.second, d));

            } else {

                std::lock_guard<std::mutex> lock(mtxDiscard_);
                numDiscarded_++;

            }

        }

    }

}

std::vector<std::pair<numSeqs_t, lenSeqs_t>> SegmentFilter::ParallelChildrenFinder::getChildren(const numSeqs_t id) {

    {
        std::lock_guard<std::mutex> lock(mtxMatches_);
        matches_.clear();
    }
    {
        std::lock_guard<std::mutex> lock(mtxDiscard_);
        numDiscarded_ = 0;
    }

    numSeqs_t numCands = 0;
    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    auto& amplicon = ac_[id];

    for (lenSeqs_t childLen = (amplicon.len > sc_.threshold) * (amplicon.len - sc_.threshold); childLen <= amplicon.len + sc_.threshold; childLen++) {

        candCnts.clear();

        SegmentFilter::addCandCnts(amplicon, childLen, sc_.threshold + sc_.extraSegs, candCnts, indices_, substrsArchive_);
        numCands += sendCandsToVerification(id, amplicon, candCnts);

    }

    numSeqs_t numDone = 0;
    while (numDone < numCands) { // wait for all verifier threads until no work is pending

        {
            std::lock_guard<std::mutex> lock(mtxMatches_);
            numDone = matches_.size();
        }
        {
            std::lock_guard<std::mutex> lock(mtxDiscard_);
            numDone += numDiscarded_;
        }

    }

    return matches_;

}

void SegmentFilter::ParallelChildrenFinder::getChildren(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children) {

    children.clear();
    {
        std::lock_guard<std::mutex> lock(mtxMatches_);
        matches_.clear();
    }
    {
        std::lock_guard<std::mutex> lock(mtxDiscard_);
        numDiscarded_ = 0;
    }

    numSeqs_t numCands = 0;
    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    auto& amplicon = ac_[id];

    for (lenSeqs_t childLen = (amplicon.len > sc_.threshold) * (amplicon.len - sc_.threshold); childLen <= amplicon.len + sc_.threshold; childLen++) {

        candCnts.clear();

        SegmentFilter::addCandCnts(amplicon, childLen, sc_.threshold + sc_.extraSegs, candCnts, indices_, substrsArchive_);
        numCands += sendCandsToVerification(id, amplicon, candCnts);

    }

    numSeqs_t numDone = 0;
    while (numDone < numCands) { // wait for all verifier threads until no work is pending

        {
            std::lock_guard<std::mutex> lock(mtxMatches_);
            numDone = matches_.size();
        }
        {
            std::lock_guard<std::mutex> lock(mtxDiscard_);
            numDone += numDiscarded_;
        }

    }

    {
        std::lock_guard<std::mutex> lock(mtxMatches_);
        children.swap(matches_);
    }

}

std::vector<std::pair<numSeqs_t, lenSeqs_t>> SegmentFilter::ParallelChildrenFinder::getChildrenTwoWay(const numSeqs_t id) {

    {
        std::lock_guard<std::mutex> lock(mtxMatches_);
        matches_.clear();
    }
    {
        std::lock_guard<std::mutex> lock(mtxDiscard_);
        numDiscarded_ = 0;
    }

    numSeqs_t numCands = 0;
    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    auto& amplicon = ac_[id];

    Segments segments(sc_.threshold + sc_.extraSegs);
    selectSegments(segments, amplicon.len, sc_.threshold, sc_.extraSegs);

    std::vector<std::string> segmentStrs(sc_.threshold + sc_.extraSegs);
    for (auto i = 0; i < sc_.threshold + sc_.extraSegs; i++) {
        segmentStrs[i] = std::string(amplicon.seq + segments[i].first, amplicon.seq + segments[i].first + segments[i].second);
    }

    for (lenSeqs_t childLen = (amplicon.len > sc_.threshold) * (amplicon.len - sc_.threshold); childLen <= amplicon.len + sc_.threshold; childLen++) {

        candCnts.clear();

        SegmentFilter::addCandCnts(amplicon, childLen, sc_.threshold + sc_.extraSegs, candCnts, indices_, substrsArchive_);
        sendCandsToVerificationTwoWay(id, amplicon, segmentStrs, candCnts, substrsArchive_[childLen][amplicon.len]);

    }

    numSeqs_t numDone = 0;
    while (numDone < numCands) { // wait for all verifier threads until no work is pending

        {
            std::lock_guard<std::mutex> lock(mtxMatches_);
            numDone = matches_.size();
        }
        {
            std::lock_guard<std::mutex> lock(mtxDiscard_);
            numDone += numDiscarded_;
        }

    }

    return matches_;

}

void SegmentFilter::ParallelChildrenFinder::getChildrenTwoWay(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children) {

    children.clear();
    {
        std::lock_guard<std::mutex> lock(mtxMatches_);
        matches_.clear();
    }
    {
        std::lock_guard<std::mutex> lock(mtxDiscard_);
        numDiscarded_ = 0;
    }

    numSeqs_t numCands = 0;
    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    auto& amplicon = ac_[id];

    Segments segments(sc_.threshold + sc_.extraSegs);
    selectSegments(segments, amplicon.len, sc_.threshold, sc_.extraSegs);

    std::vector<std::string> segmentStrs(sc_.threshold + sc_.extraSegs);
    for (auto i = 0; i < sc_.threshold + sc_.extraSegs; i++) {
        segmentStrs[i] = std::string(amplicon.seq + segments[i].first, amplicon.seq + segments[i].first + segments[i].second);
    }

    for (lenSeqs_t childLen = (amplicon.len > sc_.threshold) * (amplicon.len - sc_.threshold); childLen <= amplicon.len + sc_.threshold; childLen++) {

        candCnts.clear();

        SegmentFilter::addCandCnts(amplicon, childLen, sc_.threshold + sc_.extraSegs, candCnts, indices_, substrsArchive_);
        sendCandsToVerificationTwoWay(id, amplicon, segmentStrs, candCnts, substrsArchive_[childLen][amplicon.len]);

    }

    numSeqs_t numDone = 0;
    while (numDone < numCands) { // wait for all verifier threads until no work is pending

        {
            std::lock_guard<std::mutex> lock(mtxMatches_);
            numDone = matches_.size();
        }
        {
            std::lock_guard<std::mutex> lock(mtxDiscard_);
            numDone += numDiscarded_;
        }

    }

    {
        std::lock_guard<std::mutex> lock(mtxMatches_);
        children.swap(matches_);
    }

}


std::vector<std::pair<numSeqs_t, lenSeqs_t>> SegmentFilter::getChildren(const numSeqs_t id, const AmpliconCollection& ac, SwarmingIndices& indices, std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive, lenSeqs_t* M, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP, const SwarmClustering::SwarmConfig& sc){

    std::vector<std::pair<numSeqs_t, lenSeqs_t>> matches;
    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    auto& amplicon = ac[id];

    for (lenSeqs_t childLen = (amplicon.len > sc.threshold) * (amplicon.len - sc.threshold); childLen <= amplicon.len + sc.threshold; childLen++) {

        candCnts.clear();

        SegmentFilter::addCandCnts(amplicon, childLen, sc.threshold + sc.extraSegs, candCnts, indices, substrsArchive);
        SegmentFilter::verifyCands(amplicon, ac, candCnts, matches, sc, M, D, P, cntDiffs, cntDiffsP);

    }

    return matches;

}

void SegmentFilter::getChildren(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children, const AmpliconCollection& ac, SwarmingIndices& indices, std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive, lenSeqs_t* M, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP, const SwarmClustering::SwarmConfig& sc){

    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    auto& amplicon = ac[id];

    for (lenSeqs_t childLen = (amplicon.len > sc.threshold) * (amplicon.len - sc.threshold); childLen <= amplicon.len + sc.threshold; childLen++) {

        candCnts.clear();

        SegmentFilter::addCandCnts(amplicon, childLen, sc.threshold + sc.extraSegs, candCnts, indices, substrsArchive);
        SegmentFilter::verifyCands(amplicon, ac, candCnts, children, sc, M, D, P, cntDiffs, cntDiffsP);

    }

}

std::vector<std::pair<numSeqs_t, lenSeqs_t>> SegmentFilter::getChildrenTwoWay(const numSeqs_t id, const AmpliconCollection& ac, SwarmingIndices& indices, std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive, lenSeqs_t* M, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP, const SwarmClustering::SwarmConfig& sc){

    std::vector<std::pair<numSeqs_t, lenSeqs_t>> matches;
    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    auto& amplicon = ac[id];

    Segments segments(sc.threshold + sc.extraSegs);
    selectSegments(segments, amplicon.len, sc.threshold, sc.extraSegs);
    std::vector<std::string> segmentStrs(sc.threshold + sc.extraSegs);
    for (auto i = 0; i < sc.threshold + sc.extraSegs; i++) {
        segmentStrs[i] = std::string(amplicon.seq + segments[i].first, amplicon.seq + segments[i].first + segments[i].second);
    }

    for (lenSeqs_t childLen = (amplicon.len > sc.threshold) * (amplicon.len - sc.threshold); childLen <= amplicon.len + sc.threshold; childLen++) {

        candCnts.clear();

        SegmentFilter::addCandCnts(amplicon, childLen, sc.threshold + sc.extraSegs, candCnts, indices, substrsArchive);
        SegmentFilter::verifyCandsTwoWay(amplicon, segmentStrs, ac, candCnts, substrsArchive[childLen][amplicon.len], matches, sc, M, D, P, cntDiffs, cntDiffsP);

    }

    return matches;

}

void SegmentFilter::getChildrenTwoWay(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children, const AmpliconCollection& ac, SwarmingIndices& indices, std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive, lenSeqs_t* M, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP, const SwarmClustering::SwarmConfig& sc){

    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    auto& amplicon = ac[id];

    Segments segments(sc.threshold + sc.extraSegs);
    selectSegments(segments, amplicon.len, sc.threshold, sc.extraSegs);
    std::vector<std::string> segmentStrs(sc.threshold + sc.extraSegs);
    for (auto i = 0; i < sc.threshold + sc.extraSegs; i++) {
        segmentStrs[i] = std::string(amplicon.seq + segments[i].first, amplicon.seq + segments[i].first + segments[i].second);
    }

    for (lenSeqs_t childLen = (amplicon.len > sc.threshold) * (amplicon.len - sc.threshold); childLen <= amplicon.len + sc.threshold; childLen++) {

        candCnts.clear();

        SegmentFilter::addCandCnts(amplicon, childLen, sc.threshold + sc.extraSegs, candCnts, indices, substrsArchive);
        SegmentFilter::verifyCandsTwoWay(amplicon, segmentStrs, ac, candCnts, substrsArchive[childLen][amplicon.len], children, sc, M, D, P, cntDiffs, cntDiffsP);

    }

}


#if SUCCINCT
void SegmentFilter::prepareIndices(const AmpliconCollection& ac, SwarmingIndices& indices, std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive, const SwarmClustering::SwarmConfig& sc) { // fill indices and substrsArchive

    SharingRollingIndices<RankedAscendingLabels, RelationPrecursor> tmpIndices(sc.threshold + 1, sc.threshold + sc.extraSegs, true, false);
    std::unordered_map<lenSeqs_t, Segments> segmentsArchive;
    Segments segments(sc.threshold + sc.extraSegs);
    lenSeqs_t seqLen;

    // index all amplicons
    for (numSeqs_t curIntId = 0; curIntId < ac.size(); curIntId++) {

        seqLen = ac[curIntId].len;

        if (!tmpIndices.contains(seqLen)) {

            // inverted index
            tmpIndices.roll(seqLen, ac.numSeqsOfLen(seqLen));

            // segments information
            selectSegments(segments, seqLen, sc.threshold, sc.extraSegs);
            segmentsArchive[seqLen] = segments;

            // substrings information
            std::unordered_map<lenSeqs_t, std::vector<Substrings>>& substrs = substrsArchive[seqLen];
            for (lenSeqs_t partnerLen = (seqLen > sc.threshold) * (seqLen - sc.threshold); partnerLen <= seqLen + sc.threshold; partnerLen++) {

                std::vector<Substrings>& vec = substrs[partnerLen];
                for (lenSeqs_t segmentIndex = 0; segmentIndex < sc.threshold + sc.extraSegs; segmentIndex++) {
                    if (partnerLen <= seqLen) {
                        vec.push_back(selectSubstrs(seqLen, partnerLen, segmentIndex, sc.threshold, sc.extraSegs));
                    } else {
                        vec.push_back(selectSubstrsBackward(seqLen, partnerLen, segmentIndex, sc.threshold, sc.extraSegs));
                    }
                }

            }

        } else {
            segments = segmentsArchive[seqLen];
        }

        auto& row = tmpIndices.getIndicesRow(seqLen);
        numSeqs_t rankedId = row.shared.add(curIntId);
        for (lenSeqs_t i = 0; i < sc.threshold + sc.extraSegs; i++) {
            row.indices[i].add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), rankedId);
        }

    }

    for (lenSeqs_t len = tmpIndices.minLength(); len <= tmpIndices.maxLength(); len++) {

        if (tmpIndices.contains(len)) {

            indices.roll(len, ac.numSeqsOfLen(len));
            auto& succRow = indices.getIndicesRow(len);
            auto& tmpRow = tmpIndices.getIndicesRow(len);
            succRow.shared.swap(tmpRow.shared);

            for (lenSeqs_t i = 0; i < sc.threshold + sc.extraSegs; i++) {

                succRow.indices[i] = SwarmingInvertedIndex(tmpRow.indices[i], succRow.shared, SuccinctConfig(0,2,5,0,0,0,0));
                tmpRow.indices[i].clear();

            }

        }

    }

}
#else
void SegmentFilter::prepareIndices(const AmpliconCollection& ac, SwarmingIndices& indices, std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive, const SwarmClustering::SwarmConfig& sc) {

    std::unordered_map<lenSeqs_t, Segments> segmentsArchive;
    Segments segments(sc.threshold + sc.extraSegs);
    lenSeqs_t seqLen;

    // index all amplicons
    for (numSeqs_t curIntId = 0; curIntId < ac.size(); curIntId++) {

        seqLen = ac[curIntId].len;

        if (!indices.contains(seqLen)) {

            // inverted index
            indices.roll(seqLen);

            // segments information
            selectSegments(segments, seqLen, sc.threshold, sc.extraSegs);
            segmentsArchive[seqLen] = segments;

            // substrings information
            std::unordered_map<lenSeqs_t, std::vector<Substrings>>& substrs = substrsArchive[seqLen];
            for (lenSeqs_t partnerLen = (seqLen > sc.threshold) * (seqLen - sc.threshold); partnerLen <= seqLen + sc.threshold; partnerLen++) {

                std::vector<Substrings>& vec = substrs[partnerLen];
                for (lenSeqs_t segmentIndex = 0; segmentIndex < sc.threshold + sc.extraSegs; segmentIndex++) {
                    if (partnerLen <= seqLen) {
                        vec.push_back(selectSubstrs(seqLen, partnerLen, segmentIndex, sc.threshold, sc.extraSegs));
                    } else {
                        vec.push_back(selectSubstrsBackward(seqLen, partnerLen, segmentIndex, sc.threshold, sc.extraSegs));
                    }
                }

            }

        } else {
            segments = segmentsArchive[seqLen];
        }

        for (lenSeqs_t i = 0; i < sc.threshold + sc.extraSegs; i++) {
            indices.getIndex(seqLen, i).add(StringIteratorPair(ac[curIntId].seq + segments[i].first, ac[curIntId].seq + segments[i].first + segments[i].second), curIntId);
        }

    }

}
#endif

void SegmentFilter::swarmFilter(const AmpliconCollection& ac, std::vector<SwarmClustering::Otu*>& otus, const SwarmClustering::SwarmConfig& sc) {

    SwarmingIndices indices(sc.threshold + 1, sc.threshold + sc.extraSegs, true, false);
    std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>> substrsArchive;

    prepareIndices(ac, indices, substrsArchive, sc);

    // determine order of amplicons based on abundance (descending) without invalidating the integer (position) ids of the amplicons
    SwarmClustering::Otu* curOtu = 0;
    std::vector<SwarmClustering::OtuEntryPrecursor> tmpMembers;
    std::vector<bool> visited(ac.size(), false); // visited amplicons are already included in an OTU

    SwarmClustering::OtuEntryPrecursor curSeed, newSeed;
    bool unique;
    std::unordered_set<StringIteratorPair, hashStringIteratorPair, equalStringIteratorPair> nonUniques;
    std::vector<std::pair<numSeqs_t, lenSeqs_t>> next;
    lenSeqs_t lastGen;
    numSeqs_t pos;

    lenSeqs_t width =  indices.maxLength() + 1;
    ParallelChildrenFinder cf(ac, indices, substrsArchive, width, sc);

    // open new OTU for the amplicon with the highest abundance that is not yet included in an OTU
    const Amplicon* begin = ac.begin();
    const Amplicon* seed = begin;
    for (numSeqs_t seedIter = 0; seedIter < ac.size(); seedIter++, seed++) {

        if (!visited[seedIter]) {

            /* (a) Initialise new OTU with seed */
            curOtu = new SwarmClustering::Otu();

            newSeed.member = seed;
            newSeed.parent = newSeed.member;
            newSeed.parentDist = 0;
            newSeed.gen = 0;
            newSeed.rad = 0;
            tmpMembers.push_back(newSeed);

            visited[seedIter] = true;
#if SUCCINCT
            indices.getIndicesRow(ac[seedIter].len).shared.remove(seedIter);
#else
            {
                auto& invs = indices.getIndicesRow(ac[seedIter].len);
                for (auto i = 0; i < sc.threshold + sc.extraSegs; i++) {
                    invs[i].removeLabel(seedIter);
                }
            }
#endif
            nonUniques.clear();

            lastGen = 0;


            /* (b) BFS through 'match space' */
            pos = 0;
            while (pos < tmpMembers.size()) { // expand current OTU until no further similar amplicons can be added

                if (lastGen != tmpMembers[pos].gen) { // work through generation by decreasing abundance
                    std::sort(tmpMembers.begin() + pos, tmpMembers.end(), SwarmClustering::CompareOtuEntryPrecursorsAbund());
                }

                // get next OTU (sub)seed
                curSeed = tmpMembers[pos];

                unique = true;

                // update OTU information
                curOtu->mass += curSeed.member->abundance;

                // Consider yet unseen (unvisited) amplicons to continue the exploration.
                // An amplicon is marked as 'visited' as soon as it occurs the first time as matching partner
                // in order to prevent the algorithm from queueing it more than once coming from different amplicons.
                next = sc.filterTwoWay ? cf.getChildrenTwoWay(curSeed.member - begin) : cf.getChildren(curSeed.member - begin);
//                sc.filterTwoWay ? cf.getChildrenTwoWay(curSeed.member - begin, next) : cf.getChildren(curSeed.member - begin, next);

                for (auto matchIter = next.begin(); matchIter != next.end(); matchIter++) {

                    unique &= (matchIter->second != 0);

                    if (!visited[matchIter->first]) {

                        newSeed.member = begin + matchIter->first;
                        newSeed.parent = curSeed.member;
                        newSeed.parentDist = matchIter->second;
                        newSeed.gen = curSeed.gen + 1;
                        newSeed.rad = curSeed.rad + matchIter->second;
                        tmpMembers.push_back(newSeed);
                        visited[matchIter->first] = true;

#if SUCCINCT
                        indices.getIndicesRow(ac[matchIter->first].len).shared.remove(matchIter->first);
#else
                        auto& invs = indices.getIndicesRow(ac[matchIter->first].len);
                        for (auto i = 0; i < sc.threshold + sc.extraSegs; i++) {
                            invs[i].removeLabel(matchIter->first);
                        }
#endif

                    }

                }

                // unique sequences contribute when they occur, non-unique sequences only at their first occurrence
                unique = unique || nonUniques.insert(StringIteratorPair(curSeed.member->seq, curSeed.member->seq + curSeed.member->len)).second;
                curOtu->numUniqueSequences += unique;

                lastGen = curSeed.gen;
                pos++;

            }

            /* (c) Close the no longer extendable OTU */
            curOtu->setMembers(tmpMembers);
            std::vector<SwarmClustering::OtuEntryPrecursor>().swap(tmpMembers);
            otus.push_back(curOtu);

        }

    }

}

void SegmentFilter::swarmFilterDirectly(const AmpliconCollection& ac, std::vector<SwarmClustering::Otu*>& otus, const SwarmClustering::SwarmConfig& sc) {

    SwarmingIndices indices(sc.threshold + 1, sc.threshold + sc.extraSegs, true, false);
    std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>> substrsArchive;

    prepareIndices(ac, indices, substrsArchive, sc);

    // determine order of amplicons based on abundance (descending) without invalidating the integer (position) ids of the amplicons
    SwarmClustering::Otu* curOtu = 0;
    std::vector<SwarmClustering::OtuEntryPrecursor> tmpMembers;
    std::vector<bool> visited(ac.size(), false); // visited amplicons are already included in an OTU

    SwarmClustering::OtuEntryPrecursor curSeed, newSeed;
    bool unique;
    std::unordered_set<StringIteratorPair, hashStringIteratorPair, equalStringIteratorPair> nonUniques;
    std::vector<std::pair<numSeqs_t, lenSeqs_t>> next;
    lenSeqs_t lastGen;
    numSeqs_t pos;

    lenSeqs_t M[sc.useScore ? 1 : indices.maxLength() + 1];
    val_t D[sc.useScore ? indices.maxLength() + 1 : 1];
    val_t P[sc.useScore ? indices.maxLength() + 1 : 1];
    lenSeqs_t cntDiffs[sc.useScore ? indices.maxLength() + 1 : 1];
    lenSeqs_t cntDiffsP[sc.useScore ? indices.maxLength() + 1 : 1];
#if CHILDREN_FINDER
    ChildrenFinder cf(ac, indices, substrsArchive, sc, M, D, P, cntDiffs, cntDiffsP);
#endif

    // open new OTU for the amplicon with the highest abundance that is not yet included in an OTU
    const Amplicon* begin = ac.begin();
    const Amplicon* seed = begin;
    for (numSeqs_t seedIter = 0; seedIter < ac.size(); seedIter++, seed++) {

        if (!visited[seedIter]) {

            /* (a) Initialise new OTU with seed */
            curOtu = new SwarmClustering::Otu();

            newSeed.member = seed;
            newSeed.parent = newSeed.member;
            newSeed.parentDist = 0;
            newSeed.gen = 0;
            newSeed.rad = 0;
            tmpMembers.push_back(newSeed);

            visited[seedIter] = true;
#if SUCCINCT
            indices.getIndicesRow(ac[seedIter].len).shared.remove(seedIter);
#else
            {
                auto& invs = indices.getIndicesRow(ac[seedIter].len);
                for (auto i = 0; i < sc.threshold + sc.extraSegs; i++) {
                    invs[i].removeLabel(seedIter);
                }
            }

#endif
            nonUniques.clear();

            lastGen = 0;


            /* (b) BFS through 'match space' */
            pos = 0;
            while (pos < tmpMembers.size()) { // expand current OTU until no further similar amplicons can be added

                if (lastGen != tmpMembers[pos].gen) { // work through generation by decreasing abundance
                    std::sort(tmpMembers.begin() + pos, tmpMembers.end(), SwarmClustering::CompareOtuEntryPrecursorsAbund());
                }

                // get next OTU (sub)seed
                curSeed = tmpMembers[pos];

                unique = true;

                // update OTU information
                curOtu->mass += curSeed.member->abundance;

                // Consider yet unseen (unvisited) amplicons to continue the exploration.
                // An amplicon is marked as 'visited' as soon as it occurs the first time as matching partner
                // in order to prevent the algorithm from queueing it more than once coming from different amplicons.
#if CHILDREN_FINDER
                next = sc.filterTwoWay ? cf.getChildrenTwoWay(curSeed.member - begin) : cf.getChildren(curSeed.member - begin);
//                sc.filterTwoWay ? cf.getChildrenTwoWay(curSeed.member - begin, next) : cf.getChildren(curSeed.member - begin, next);
#else
                next = sc.filterTwoWay ? getChildrenTwoWay(curSeed.member - begin, next, ac, indices, substrsArchive, M, D, P, cntDiffs, cntDiffsP, sc) : getChildren(curSeed.member - begin, ac, indices, substrsArchive, M, D, P, cntDiffs, cntDiffsP, sc);
//                sc.filterTwoWay ? getChildrenTwoWay(curSeed.member - begin, next, ac, indices, substrsArchive, M, D, P, cntDiffs, cntDiffsP, sc) : getChildren(curSeed.member - begin, ac, indices, substrsArchive, M, D, P, cntDiffs, cntDiffsP, sc);
#endif
                for (auto matchIter = next.begin(); matchIter != next.end(); matchIter++) {

                    unique &= (matchIter->second != 0);

                    if (!visited[matchIter->first]) {

                        newSeed.member = begin + matchIter->first;
                        newSeed.parent = curSeed.member;
                        newSeed.parentDist = matchIter->second;
                        newSeed.gen = curSeed.gen + 1;
                        newSeed.rad = curSeed.rad + matchIter->second;
                        tmpMembers.push_back(newSeed);
                        visited[matchIter->first] = true;

#if SUCCINCT
                        indices.getIndicesRow(ac[matchIter->first].len).shared.remove(matchIter->first);
#else
                        auto& invs = indices.getIndicesRow(ac[matchIter->first].len);
                        for (auto i = 0; i < sc.threshold + sc.extraSegs; i++) {
                            invs[i].removeLabel(matchIter->first);
                        }
#endif

                    }

                }

                // unique sequences contribute when they occur, non-unique sequences only at their first occurrence
                unique = unique || nonUniques.insert(StringIteratorPair(curSeed.member->seq, curSeed.member->seq + curSeed.member->len)).second;
                curOtu->numUniqueSequences += unique;

                lastGen = curSeed.gen;
                pos++;

            }

            /* (c) Close the no longer extendable OTU */
            curOtu->setMembers(tmpMembers);
            std::vector<SwarmClustering::OtuEntryPrecursor>().swap(tmpMembers);
            otus.push_back(curOtu);

        }

    }

}

}
