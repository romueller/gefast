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
#include "../include/SwarmClustering.hpp"
#include "../include/SwarmingSegmentFilter.hpp"

#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <thread>
#include <unordered_set>

namespace GeFaST {

void SwarmClustering::explorePool(const AmpliconCollection& ac, Matches& matches, std::vector<Otu*>& otus, const SwarmConfig& sc) {

    // determine order of amplicons based on abundance (descending) without invalidating the integer (position) ids of the amplicons
    std::vector<numSeqs_t> index(ac.size());
    std::iota(std::begin(index), std::end(index), 0);
    std::sort(index.begin(), index.end(), CompareIndicesAbund(ac));

    Otu* curOtu = 0;
    std::vector<SwarmClustering::OtuEntryPrecursor> tmpMembers;
    std::vector<bool> visited(ac.size(), false); // visited amplicons are already included in an OTU

    OtuEntryPrecursor curSeed, newSeed;
    bool unique;
    std::unordered_set<StringIteratorPair, hashStringIteratorPair, equalStringIteratorPair> uniqueSeqs;
    std::vector<std::pair<numSeqs_t, lenSeqs_t>> next;
    lenSeqs_t lastGen;
    numSeqs_t pos;

    // open new OTU for the amplicon with the highest abundance that is not yet included in an OTU
    const Amplicon* begin = ac.begin();
    for (auto seedIter = index.begin(); seedIter != index.end(); seedIter++) {

        if (!visited[*seedIter]) {

            /* (a) Initialise new OTU with seed */
            curOtu = new Otu();

            newSeed.member = &ac[*seedIter];
            newSeed.parent = newSeed.member;
            newSeed.parentDist = 0;
            newSeed.gen = 0;
            newSeed.rad = 0;
            tmpMembers.push_back(newSeed);

            visited[*seedIter] = true;
            uniqueSeqs.insert(StringIteratorPair(newSeed.member->seq, newSeed.member->seq + newSeed.member->len));

            lastGen = 0;


            /* (b) BFS through 'match space' */
            pos = 0;
            while (pos < tmpMembers.size()) { // expand current OTU until no further similar amplicons can be added

                if (lastGen != tmpMembers[pos].gen) { // work through generation by decreasing abundance

                    uniqueSeqs.clear();
                    std::sort(tmpMembers.begin() + pos, tmpMembers.end(), CompareOtuEntryPrecursorsAbund());

                }

                // get next OTU (sub)seed
                curSeed = tmpMembers[pos];

                // unique sequences contribute when they occur, non-unique sequences only at their first occurrence
                unique = (curSeed.parentDist != 0) && uniqueSeqs.insert(StringIteratorPair(curSeed.member->seq, curSeed.member->seq + curSeed.member->len)).second;

                // update OTU information
                curOtu->mass += curSeed.member->abundance;

                // Consider yet unseen (unvisited) amplicons to continue the exploration.
                // An amplicon is marked as 'visited' as soon as it occurs the first time as matching partner
                // in order to prevent the algorithm from queueing it more than once coming from different amplicons.
                next = matches.getMatchesOfAugmented(curSeed.member - begin);
                for (auto matchIter = next.begin(); matchIter != next.end(); matchIter++) {

                    unique &= (matchIter->second != 0);

                    if (!visited[matchIter->first] && (sc.noOtuBreaking || ac[matchIter->first].abundance <= curSeed.member->abundance)) {

                        newSeed.member = begin + matchIter->first;
                        newSeed.parent = curSeed.member;
                        newSeed.parentDist = matchIter->second;
                        newSeed.gen = curSeed.gen + 1;
                        newSeed.rad = curSeed.rad + matchIter->second;
                        tmpMembers.push_back(newSeed);
                        visited[matchIter->first] = true;

                    }

                }

                curOtu->numUniqueSequences += unique || (curSeed.gen == 0);

                lastGen = curSeed.gen;
                pos++;

            }

            /* (c) Close the no longer extendable OTU */
            uniqueSeqs.clear();
            curOtu->setMembers(tmpMembers);
            std::vector<SwarmClustering::OtuEntryPrecursor>().swap(tmpMembers);
            otus.push_back(curOtu);

        }

    }

}


void SwarmClustering::fastidiousIndexOtu(PrecursorIndices& indices, std::unordered_map<lenSeqs_t, Segments>& segmentsArchive, const AmpliconCollection& ac, Otu& otu, std::vector<GraftCandidate>& graftCands, const SwarmConfig& sc) {

    auto begin = ac.begin();
    for (numSeqs_t m = 0; m < otu.numMembers; m++) {

        auto ampl = otu.members[m].member;
        Segments& segments = segmentsArchive[ampl->len];

        if (segments.size() == 0) {

            indices.roll(ampl->len);
            segments = Segments(sc.fastidiousThreshold + sc.extraSegs);
            selectSegments(segments, ampl->len, sc.fastidiousThreshold, sc.extraSegs);

        }

        auto& row = indices.getIndicesRow(ampl->len);
        for (lenSeqs_t i = 0; i < sc.fastidiousThreshold + sc.extraSegs; i++) {
            row[i].add(StringIteratorPair(ampl->seq + segments[i].first, ampl->seq + segments[i].first + segments[i].second), ampl - begin);
        }

        graftCands[ampl - begin].childOtu = &otu;
        graftCands[ampl - begin].childMember = ampl;

    }

}

inline bool compareCandidates(const Amplicon& newCand, const Amplicon& oldCand) {
    return (newCand.abundance > oldCand.abundance) || ((newCand.abundance == oldCand.abundance) && (strcmp(newCand.id, oldCand.id) < 0));
}

void SwarmClustering::verifyFastidious(const AmpliconPools& pools, const AmpliconCollection& acOtus, const AmpliconCollection& acIndices, std::vector<GraftCandidate>& graftCands, Buffer<CandidateFastidious>& buf, const lenSeqs_t width, const lenSeqs_t t, std::mutex& mtx) {

    CandidateFastidious c;
    Buffer<CandidateFastidious> localBuffer;
    lenSeqs_t M[width]; // reusable DP-matrix (wide enough for all possible calculations for this AmpliconCollection)

    while (!buf.isClosed() || buf.syncSize() > 0) {

        buf.syncSwapContents(localBuffer);

        while (localBuffer.size() > 0) {

            c = localBuffer.pop();

            for (auto childIter = c.children.begin(); childIter != c.children.end(); childIter++) {

                std::unique_lock<std::mutex> lock(mtx);
                if ((graftCands[*childIter].parentOtu == 0) || compareCandidates(*c.parentMember->member, *graftCands[*childIter].parentMember->member)) {

                    lock.unlock();
                    if (Verification::computeLengthAwareRow(c.parentMember->member->seq, c.parentMember->member->len, acIndices[*childIter].seq, acIndices[*childIter].len, t, M) <= t) {

                        lock.lock();
                        if ((graftCands[*childIter].parentOtu == 0) || compareCandidates(*c.parentMember->member, *graftCands[*childIter].parentMember->member)) {

                            graftCands[*childIter].parentOtu = c.parentOtu;
                            graftCands[*childIter].parentMember = c.parentMember;

                        }
                        lock.unlock();

                    }

                }

            }

        }

    }

}

void SwarmClustering::verifyGotohFastidious(const AmpliconPools& pools, const AmpliconCollection& acOtus, const AmpliconCollection& acIndices, std::vector<GraftCandidate>& graftCands, Buffer<CandidateFastidious>& buf, const lenSeqs_t width, const lenSeqs_t t, const Verification::Scoring& scoring, std::mutex& mtx) {

    CandidateFastidious c;
    Buffer<CandidateFastidious> localBuffer;

    // reusable DP-matrices (wide enough for all possible calculations for this AmpliconCollection)
    val_t D[width];
    val_t P[width];
    lenSeqs_t cntDiffs[width];
    lenSeqs_t cntDiffsP[width];

    while (!buf.isClosed() || buf.syncSize() > 0) {

        buf.syncSwapContents(localBuffer);

        while (localBuffer.size() > 0) {

            c = localBuffer.pop();

            for (auto childIter = c.children.begin(); childIter != c.children.end(); childIter++) {

                std::unique_lock<std::mutex> lock(mtx);
                if ((graftCands[*childIter].parentOtu == 0) || compareCandidates(*c.parentMember->member, *graftCands[*childIter].parentMember->member)) {

                    lock.unlock();
                    if (Verification::computeGotohLengthAwareEarlyRow(c.parentMember->member->seq, c.parentMember->member->len, acIndices[*childIter].seq, acIndices[*childIter].len, t, scoring, D, P, cntDiffs, cntDiffsP) <= t) {

                        lock.lock();
                        if ((graftCands[*childIter].parentOtu == 0) || compareCandidates(*c.parentMember->member, *graftCands[*childIter].parentMember->member)) {

                            graftCands[*childIter].parentOtu = c.parentOtu;
                            graftCands[*childIter].parentMember = c.parentMember;

                        }
                        lock.unlock();

                    }

                }

            }

        }

    }

}

void SwarmClustering::fastidiousCheckOtus(RotatingBuffers<CandidateFastidious>& cbs, const std::vector<Otu*>& otus, const AmpliconCollection& acOtus, IndicesFastidious& indices, const AmpliconCollection& acIndices, std::vector<GraftCandidate>& graftCands, const SwarmConfig& sc) {

    std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>> substrsArchive;
    std::vector<numSeqs_t> candCnts;
    lenSeqs_t seqLen;
    StringIteratorPair sip;

    std::vector<CandidateFastidious> localCands;

    for (auto otuIter = otus.begin(); otuIter != otus.end(); otuIter++) {

        if ((*otuIter)->mass >= sc.boundary) { // for each heavy OTU of the pool ...

            for (numSeqs_t m = 0; m < (*otuIter)->numMembers; m++) { // ... consider every amplicon in the OTU and ...

                auto ampl = (*otuIter)->members[m].member;
                seqLen = ampl->len;

                std::unordered_map<lenSeqs_t, std::vector<Substrings>>& substrs = substrsArchive[seqLen];

                // on reaching new length group, open new inverted indices
                if (substrs.empty()) {

                    // ... and determine position information shared by all amplicons of this length
                    for (lenSeqs_t partnerLen = (seqLen > sc.fastidiousThreshold) * (seqLen - sc.fastidiousThreshold); partnerLen <= seqLen + sc.fastidiousThreshold; partnerLen++) {

                        std::vector<Substrings>& vec = substrs[partnerLen];
                        for (lenSeqs_t segmentIndex = 0; segmentIndex < sc.fastidiousThreshold + sc.extraSegs; segmentIndex++) {
                            if (partnerLen <= seqLen) {
                                vec.push_back(selectSubstrs(seqLen, partnerLen, segmentIndex, sc.fastidiousThreshold, sc.extraSegs));
                            } else {
                                vec.push_back(selectSubstrsBackward(seqLen, partnerLen, segmentIndex, sc.fastidiousThreshold, sc.extraSegs));
                            }
                        }

                    }

                }

                localCands.push_back(CandidateFastidious(*otuIter, (*otuIter)->members + m));

                for (lenSeqs_t len = (seqLen > sc.fastidiousThreshold) * (seqLen - sc.fastidiousThreshold); len <= seqLen + sc.fastidiousThreshold; len++) { // ... search for graft candidates among the amplicons in light OTUs

                    for (lenSeqs_t i = 0; i < sc.fastidiousThreshold + sc.extraSegs; i++) { // ... and apply segment filter for each segment

                        Substrings& subs = substrs[len][i];
                        InvertedIndexFastidious& inv = indices.getIndex(len, i);
                        sip.first = ampl->seq + subs.first;
                        sip.second = sip.first + subs.len;

                        for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++, sip.first++, sip.second++) {
                            inv.addLabelCountsOf(sip, candCnts);
                        }

                    }

                    std::sort(candCnts.begin(), candCnts.end());

                    // general pigeonhole principle: for being a candidate, at least sc.extraSegs segments have to be matched
                    lenSeqs_t cnt = 0;
                    numSeqs_t prevCand = (candCnts.size() > 0) ? candCnts.front() : 0;
                    for (auto candId : candCnts) {

                        if (prevCand != candId) {

#if QGRAM_FILTER
                            if ((cnt >= sc.extraSegs) && (qgram_diff(*ampl, acIndices[prevCand]) <= sc.fastidiousThreshold)) {
#else
                            if (cnt >= sc.extraSegs) {
#endif
                                localCands.back().children.push_back(prevCand);
                            }

                            cnt = 1;
                            prevCand = candId;

                        } else {
                            cnt++;
                        }

                    }

#if QGRAM_FILTER
                    if ((cnt >= sc.extraSegs) && (qgram_diff(*ampl, acIndices[prevCand]) <= sc.fastidiousThreshold)) {
#else
                    if (cnt >= sc.extraSegs) {
#endif
                        localCands.back().children.push_back(prevCand);
                    }

                    candCnts.clear();

                }

                cbs.push(localCands);
                localCands = std::vector<CandidateFastidious>();

            }

        }

    }

}

void SwarmClustering::fastidiousCheckOtusDirectly(const AmpliconPools& pools, const std::vector<Otu*>& otus, const AmpliconCollection& acOtus, IndicesFastidious& indices, const AmpliconCollection& acIndices, std::vector<GraftCandidate>& graftCands, const lenSeqs_t width, std::mutex& graftCandsMtx, const SwarmConfig& sc) {

    std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>> substrsArchive;
    std::vector<numSeqs_t> candCnts;
    lenSeqs_t seqLen;
    StringIteratorPair sip;

    lenSeqs_t M[sc.useScore ? 1 : width];
    val_t D[sc.useScore? width : 1];
    val_t P[sc.useScore? width : 1];
    lenSeqs_t cntDiffs[sc.useScore? width : 1];
    lenSeqs_t cntDiffsP[sc.useScore? width : 1];

    for (auto otuIter = otus.begin(); otuIter != otus.end(); otuIter++) {

        if ((*otuIter)->mass >= sc.boundary) { // for each heavy OTU of the pool ...

            for (numSeqs_t m = 0; m < (*otuIter)->numMembers; m++) { // ... consider every amplicon in the OTU and ...

                auto ampl = (*otuIter)->members[m].member;
                seqLen = ampl->len;

                std::unordered_map<lenSeqs_t, std::vector<Substrings>>& substrs = substrsArchive[seqLen];

                // on reaching new length group, open new inverted indices
                if (substrs.empty()) {

                    // ... and determine position information shared by all amplicons of this length
                    for (lenSeqs_t partnerLen = (seqLen > sc.fastidiousThreshold) * (seqLen - sc.fastidiousThreshold); partnerLen <= seqLen + sc.fastidiousThreshold; partnerLen++) {

                        std::vector<Substrings>& vec = substrs[partnerLen];
                        for (lenSeqs_t segmentIndex = 0; segmentIndex < sc.fastidiousThreshold + sc.extraSegs; segmentIndex++) {
                            if (partnerLen <= seqLen) {
                                vec.push_back(selectSubstrs(seqLen, partnerLen, segmentIndex, sc.fastidiousThreshold, sc.extraSegs));
                            } else {
                                vec.push_back(selectSubstrsBackward(seqLen, partnerLen, segmentIndex, sc.fastidiousThreshold, sc.extraSegs));
                            }
                        }

                    }

                }


                for (lenSeqs_t len = (seqLen > sc.fastidiousThreshold) * (seqLen - sc.fastidiousThreshold); len <= seqLen + sc.fastidiousThreshold; len++) { // ... search for graft candidates among the amplicons in light OTUs

                    for (lenSeqs_t i = 0; i < sc.fastidiousThreshold + sc.extraSegs; i++) { // ... and apply segment filter for each segment

                        Substrings& subs = substrs[len][i];
                        InvertedIndexFastidious& inv = indices.getIndex(len, i);
                        sip.first = ampl->seq + subs.first;
                        sip.second = sip.first + subs.len;

                        for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++, sip.first++, sip.second++) {
                            inv.addLabelCountsOf(sip, candCnts);
                        }

                    }

                    std::sort(candCnts.begin(), candCnts.end());

                    // general pigeonhole principle: for being a candidate, at least sc.extraSegs segments have to be matched
                    lenSeqs_t cnt = 0;
                    numSeqs_t prevCand = (candCnts.size() > 0) ? candCnts.front() : 0;
                    for (auto candId : candCnts) {

                        if (prevCand != candId) {

//                            if ((cnt >= sc.extraSegs)
//                                    &&((graftCands[prevCand].parentOtu == 0) || compareCandidates(*ampl, *graftCands[prevCand].parentMember))
//                                    && ((useScore ?
//                                            Verification::computeGotohLengthAwareEarlyRow(ampl->seq, ampl->len, acIndices[prevCand].seq, acIndices[prevCand].len, sc.fastidiousThreshold, sc.scoring, D, P, cntDiffs, cntDiffsP)
//                                          : Verification::computeLengthAwareRow(ampl->seq, ampl->len, acIndices[prevCand].seq, acIndices[prevCand].len, sc.fastidiousThreshold, M)) <= sc.fastidiousThreshold)) {
//
//                                        graftCands[prevCand].parentOtu = *otuIter;
//                                        graftCands[prevCand].parentMember = ampl;
//
//                            }

                            std::unique_lock<std::mutex> lock(graftCandsMtx);
#if QGRAM_FILTER
                            if ((cnt >= sc.extraSegs) && ((graftCands[prevCand].parentOtu == 0) || compareCandidates(*ampl, *graftCands[prevCand].parentMember->member)) && (qgram_diff(*ampl, acIndices[prevCand]) <= sc.fastidiousThreshold)) {
#else
                            if ((cnt >= sc.extraSegs) && ((graftCands[prevCand].parentOtu == 0) || compareCandidates(*ampl, *graftCands[prevCand].parentMember->member))) {
#endif

                                lock.unlock();
                                if ((sc.useScore ?
                                      Verification::computeGotohLengthAwareEarlyRow(ampl->seq, ampl->len, acIndices[prevCand].seq, acIndices[prevCand].len, sc.fastidiousThreshold, sc.scoring, D, P, cntDiffs, cntDiffsP)
                                    : Verification::computeLengthAwareRow(ampl->seq, ampl->len, acIndices[prevCand].seq, acIndices[prevCand].len, sc.fastidiousThreshold, M)) <= sc.fastidiousThreshold) {

                                    lock.lock();
                                    if (((graftCands[prevCand].parentOtu == 0) || compareCandidates(*ampl, *graftCands[prevCand].parentMember->member))) {

                                        graftCands[prevCand].parentOtu = *otuIter;
                                        graftCands[prevCand].parentMember = (*otuIter)->members + m;

                                    }
                                    lock.unlock();

                                }

                            }

                            cnt = 1;
                            prevCand = candId;

                        } else {
                            cnt++;
                        }

                    }

                    std::unique_lock<std::mutex> lock(graftCandsMtx);
#if QGRAM_FILTER
                    if ((cnt >= sc.extraSegs) && ((graftCands[prevCand].parentOtu == 0) || compareCandidates(*ampl, *graftCands[prevCand].parentMember->member)) && (qgram_diff(*ampl, acIndices[prevCand]) <= sc.fastidiousThreshold)) {
#else
                    if ((cnt >= sc.extraSegs) && ((graftCands[prevCand].parentOtu == 0) || compareCandidates(*ampl, *graftCands[prevCand].parentMember->member))) {
#endif

                        lock.unlock();
                        if ((sc.useScore ?
                               Verification::computeGotohLengthAwareEarlyRow(ampl->seq, ampl->len, acIndices[prevCand].seq, acIndices[prevCand].len, sc.fastidiousThreshold, sc.scoring, D, P, cntDiffs, cntDiffsP)
                             : Verification::computeLengthAwareRow(ampl->seq, ampl->len, acIndices[prevCand].seq, acIndices[prevCand].len, sc.fastidiousThreshold, M)) <= sc.fastidiousThreshold) {

                            lock.lock();
                            if (((graftCands[prevCand].parentOtu == 0) || compareCandidates(*ampl, *graftCands[prevCand].parentMember->member))) {

                                graftCands[prevCand].parentOtu = *otuIter;
                                graftCands[prevCand].parentMember = (*otuIter)->members + m;

                            }
                            lock.unlock();

                        }

                    }

                    candCnts.clear();

                }

            }

        }

    }

}

void SwarmClustering::checkAndVerify(const AmpliconPools& pools, const std::vector<Otu*>& otus, const AmpliconCollection& acOtus, IndicesFastidious& indices, const AmpliconCollection& acIndices, std::vector<GraftCandidate>& graftCands, const lenSeqs_t width, std::mutex& graftCandsMtx, const SwarmConfig& sc) {

    if (sc.numThreadsPerCheck == 1) {
        fastidiousCheckOtusDirectly(pools, otus, acOtus, indices, acIndices, graftCands, width, graftCandsMtx, sc);
    } else {

        RotatingBuffers<CandidateFastidious> cbs = RotatingBuffers<CandidateFastidious>(sc.numThreadsPerCheck);
        std::thread verifierThreads[sc.numThreadsPerCheck];

        for (unsigned long v = 0; v < sc.numThreadsPerCheck; v++) {
            verifierThreads[v] = sc.useScore ?
                                   std::thread(&SwarmClustering::verifyGotohFastidious, std::ref(pools), std::ref(acOtus), std::ref(acIndices), std::ref(graftCands), std::ref(cbs.getBuffer(v)), width, sc.fastidiousThreshold, std::ref(sc.scoring), std::ref(graftCandsMtx))
                                 : std::thread(&SwarmClustering::verifyFastidious, std::ref(pools), std::ref(acOtus), std::ref(acIndices), std::ref(graftCands), std::ref(cbs.getBuffer(v)), width, sc.fastidiousThreshold, std::ref(graftCandsMtx));
        }

        fastidiousCheckOtus(cbs, otus, acOtus, indices, acIndices, graftCands, sc);
        cbs.close();

        for (unsigned long v = 0; v < sc.numThreadsPerCheck; v++) {
            verifierThreads[v].join();
        }

    }

}

void SwarmClustering::determineGrafts(const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, std::vector<GraftCandidate>& allGraftCands, const numSeqs_t p, std::mutex& allGraftCandsMtx, const SwarmConfig& sc) {

    AmpliconCollection* ac = pools.get(p);
    IndicesFastidious indices(2 * sc.fastidiousThreshold + 1, sc.fastidiousThreshold + sc.extraSegs, true, false);
    std::vector<GraftCandidate> graftCands(ac->size()); // initially, graft candidates for all amplicons of the pool are "empty"

    // a) Index amplicons of all light OTUs of the current pool
#if SUCCINCT_FASTIDIOUS
    {
        SuccinctConfig succinctConfig(0, 2, 5, 0, 0, 0, 5);
        PrecursorIndices tmpIndices(2 * sc.fastidiousThreshold + 1, sc.fastidiousThreshold + sc.extraSegs, true, false);
        std::unordered_map<lenSeqs_t, Segments> segmentsArchive;
        for (auto otuIter = otus[p].begin(); otuIter != otus[p].end(); otuIter++) {

            if ((*otuIter)->mass < sc.boundary) {
                fastidiousIndexOtu(tmpIndices, segmentsArchive, *ac, *(*otuIter), graftCands, sc);
            }

        }

        for (lenSeqs_t len = tmpIndices.minLength(); len <= tmpIndices.maxLength(); len++) {

            if (tmpIndices.contains(len)) {

                indices.roll(len, ac->numSeqsOfLen(len));
                auto& succRow = indices.getIndicesRow(len);
                auto& tmpRow = tmpIndices.getIndicesRow(len);
                succRow.shared = RankedAscendingLabels(tmpRow[0].pairs);

                for (lenSeqs_t i = 0; i < sc.fastidiousThreshold + sc.extraSegs; i++) {

                    std::sort(tmpRow[i].pairs.begin(), tmpRow[i].pairs.end(),
                              [](const std::pair<numSeqs_t, numSeqs_t>& lhs, const std::pair<numSeqs_t, numSeqs_t>& rhs) {
                                  return lhs.second < rhs.second || (lhs.second == rhs.second && lhs.first < rhs.second);
                              }
                    );
                    for (numSeqs_t j = 0; j < tmpRow[i].pairs.size(); j++) {
                        tmpRow[i].pairs[j].second = j;
                    }

                    succRow.indices[i] = InvertedIndexFastidious(tmpRow[i], succRow.shared, succinctConfig);
                    tmpRow[i].clear();

                }

            }

        }
    }
#else
    {
        std::unordered_map<lenSeqs_t, Segments> segmentsArchive;
        for (auto otuIter = otus[p].begin(); otuIter != otus[p].end(); otuIter++) {

            if ((*otuIter)->mass < sc.boundary) {
                fastidiousIndexOtu(indices, segmentsArchive, *ac, *(*otuIter), graftCands, sc);
            }

        }
    }
#endif

    // determine maximum sequence length to adjust data structures in subsequently called methods
    lenSeqs_t maxLen = 0;
    for (auto iter = ac->begin(); iter != ac->end(); iter++) {
        maxLen = std::max(maxLen, iter->len);
    }

    // b) Search with amplicons of all heavy OTUs of current and neighbouring pools
    std::mutex graftCandsMtx;
    lenSeqs_t halfRange = sc.fastidiousThreshold / (sc.threshold + 1);
    lenSeqs_t minP = (p > halfRange) ? (p - halfRange) : 0;
    lenSeqs_t maxP = std::min(p + halfRange, pools.numPools() - 1);
#if FASTIDIOUS_PARALLEL_CHECK

    switch (sc.fastidiousCheckingMode) {

        case 0: {

            for (lenSeqs_t q = minP; q < p; q++) {
                checkAndVerify(pools, otus[q], *(pools.get(q)), indices, *ac, graftCands, maxLen + 1, graftCandsMtx, sc);
            }

            checkAndVerify(pools, otus[p], *ac, indices, *ac, graftCands, maxLen + 1, graftCandsMtx, sc);

            for (lenSeqs_t q = p + 1; q <= maxP; q++) {

                AmpliconCollection* succAc = pools.get(q);

                // adjust maxLen as successor amplicon collection contains longer sequences
                for (auto iter = succAc->begin(); iter != succAc->end(); iter++) {
                    maxLen = std::max(maxLen, iter->len);
                }

                checkAndVerify(pools, otus[q], *succAc, indices, *ac, graftCands, maxLen + 1, graftCandsMtx, sc);

            }

            break;

        }

        case 1: {

            std::thread t(&SwarmClustering::checkAndVerify, std::ref(pools), std::ref(otus[p]), std::ref(*ac), std::ref(indices), std::ref(*ac), std::ref(graftCands), maxLen + 1, std::ref(graftCandsMtx), std::ref(sc));

            for (lenSeqs_t q = minP; q < p; q++) {
                checkAndVerify(pools, otus[q], *(pools.get(q)), indices, *ac, graftCands, maxLen + 1, graftCandsMtx, sc);
            }

            for (lenSeqs_t q = p + 1; q <= maxP; q++) {

                AmpliconCollection* succAc = pools.get(q);

                // adjust maxLen as successor amplicon collection contains longer sequences
                for (auto iter = succAc->begin(); iter != succAc->end(); iter++) {
                    maxLen = std::max(maxLen, iter->len);
                }

                checkAndVerify(pools, otus[q], *succAc, indices, *ac, graftCands, maxLen + 1, graftCandsMtx, sc);

            }

            t.join();

            break;

        }

        default: {

            std::thread pred, succ;

            std::thread self(&SwarmClustering::checkAndVerify, std::ref(pools), std::ref(otus[p]), std::ref(*ac), std::ref(indices), std::ref(*ac), std::ref(graftCands), maxLen + 1, std::ref(graftCandsMtx), std::ref(sc));

            for (lenSeqs_t d = 1; d <= halfRange; d++) {

                if (d <= p - minP) {
                    pred = std::thread(&SwarmClustering::checkAndVerify, std::ref(pools), std::ref(otus[p - d]), std::ref(*(pools.get(p - d))), std::ref(indices), std::ref(*ac), std::ref(graftCands), maxLen + 1, std::ref(graftCandsMtx), std::ref(sc));
                }

                if (d <= maxP - p) {

                    AmpliconCollection* succAc = pools.get(p + d);

                    // adjust maxLen as successor amplicon collection contains longer sequences
                    for (auto iter = succAc->begin(); iter != succAc->end(); iter++) {
                        maxLen = std::max(maxLen, iter->len);
                    }

                    succ = std::thread(&SwarmClustering::checkAndVerify, std::ref(pools), std::ref(otus[p + d]), std::ref(*succAc), std::ref(indices), std::ref(*ac), std::ref(graftCands), maxLen + 1, std::ref(graftCandsMtx), std::ref(sc));

                }

                if (d <= p - minP) {
                    pred.join();
                }

                if (d <= maxP - p) {
                    succ.join();
                }

            }

            self.join();

            break;

        }

    }

#else

    if (p > 0) {
        checkAndVerify(pools, otus[p - 1], *(pools.get(p - 1)), indices, *ac, graftCands, maxLen + 1, graftCandsMtx, sc);
    }

    checkAndVerify(pools, otus[p], *ac, indices, *ac, graftCands, maxLen + 1, graftCandsMtx, sc);

    if (p < pools.numPools() - 1) {

        AmpliconCollection* succAc = pools.get(p + 1);

        // adjust maxLen as successor amplicon collection contains longer sequences
        for (auto iter = succAc->begin(); iter != succAc->end(); iter++) {
            maxLen = std::max(maxLen, iter->len);
        }

        checkAndVerify(pools, otus[p + 1], *succAc, indices, *ac, graftCands, maxLen + 1, graftCandsMtx, sc);

    }

#endif

    // c) Collect the (actual = non-empty) graft candidates for the current pool
    auto newEnd = std::remove_if(
            graftCands.begin(),
            graftCands.end(),
            [](GraftCandidate& gc) {
                return gc.parentOtu == 0;
            });

    std::lock_guard<std::mutex> lock(allGraftCandsMtx);
    allGraftCands.reserve(allGraftCands.size() + std::distance(graftCands.begin(), newEnd));
    std::move(graftCands.begin(), newEnd, std::back_inserter(allGraftCands));

}

void SwarmClustering::graftOtus(numSeqs_t& maxSize, numSeqs_t& numOtus, const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, const SwarmConfig& sc) {

    std::vector<GraftCandidate> allGraftCands;
    std::mutex allGraftCandsMtx;

#if FASTIDIOUS_PARALLEL_POOL

    std::thread grafters[sc.numGrafters];
    unsigned long r = 0;
    for (; r + sc.numGrafters <= pools.numPools(); r += sc.numGrafters) {

        for (unsigned long g = 0; g < sc.numGrafters; g++) {
            grafters[g] = std::thread(&SwarmClustering::determineGrafts, std::ref(pools), std::ref(otus), std::ref(allGraftCands), r + g, std::ref(allGraftCandsMtx), std::ref(sc));
        }
        for (unsigned long g = 0; g < sc.numGrafters; g++) {
            grafters[g].join();
        }

    }

    for (unsigned long g = 0; g < pools.numPools() % sc.numGrafters; g++) {
        grafters[g] = std::thread(&SwarmClustering::determineGrafts, std::ref(pools), std::ref(otus), std::ref(allGraftCands), r + g, std::ref(allGraftCandsMtx), std::ref(sc));
    }
    for (unsigned long g = 0; g < pools.numPools() % sc.numGrafters; g++) {
        grafters[g].join();
    }

#else

    for (numSeqs_t p = 0; p < pools.numPools(); p++) {
        determineGrafts(pools, otus, allGraftCands, p, allGraftCandsMtx, sc);
    }

#endif

    // Sort all graft candidates and perform actual grafting
    std::cout << "Got " << allGraftCands.size() << " graft candidates." << std::endl;
    std::sort(allGraftCands.begin(), allGraftCands.end(), CompareGraftCandidatesAbund());
    Otu* parentOtu = 0;
    Otu* childOtu = 0;
    numSeqs_t numGrafts = 0;
    for (auto graftIter = allGraftCands.begin(); graftIter != allGraftCands.end(); graftIter++) {

        if (!(graftIter->childOtu->attached())) {

            parentOtu = graftIter->parentOtu;
            childOtu = graftIter->childOtu;

            // "attach the seed of the light swarm to the tail of the heavy swarm"
            // OTU entries are moved unchanged (entry of seed of attached swarm stays 'incomplete', grafting 'link' is not recorded)
            parentOtu->attach(childOtu, graftIter->parentMember, graftIter->childMember);

            // update stats
            maxSize = std::max(maxSize, parentOtu->numTotalMembers());

            numGrafts++;
            numOtus--;

        }

    }

    std::cout << "Made " << numGrafts << " grafts." << std::endl;

}


void SwarmClustering::processOtus(const AmpliconPools& pools, std::vector<std::vector<Otu*>>& otus, const SwarmConfig& sc) {

    // make OTU IDs unique over all pools (so far IDs start at 1 in each pool) (currently commented out)
    // add pool IDs and determine some overall statistics
    numSeqs_t numOtus = 0;
    numSeqs_t numOtusAdjusted = 0;
    numSeqs_t numAmplicons = 0;
    numSeqs_t maxSize = 0;
    numSeqs_t maxGen = 0;

    for (numSeqs_t p = 0; p < pools.numPools(); p++) {

        for (auto otuIter = otus[p].begin(); otuIter != otus[p].end(); otuIter++) {

            maxSize = std::max(maxSize, (*otuIter)->numMembers); // sufficient prior to fastidious grafting
            maxGen = std::max(maxGen, (*otuIter)->maxGen());

        }

        numOtus += otus[p].size();
        numAmplicons += pools.get(p)->size();

    }

    numOtusAdjusted = numOtus;

    /* (b) Optional (second) clustering phase of swarm */
    if (sc.fastidious) {

        std::cout << "Results before fastidious processing: " << std::endl;
        std::cout << "Number of swarms: " << numOtus << std::endl;
        std::cout << "Largest swarms: " << maxSize << std::endl;

        std::cout << "Counting amplicons in heavy and light swarms..." << std::endl;
        numSeqs_t numLightOtus = 0;
        numSeqs_t numAmplLightOtus = 0;
        for (numSeqs_t p = 0; p < pools.numPools(); p++) {
            for (auto otuIter = otus[p].begin(); otuIter != otus[p].end(); otuIter++) {

                numLightOtus += ((*otuIter)->mass < sc.boundary);
                numAmplLightOtus += ((*otuIter)->mass < sc.boundary) * (*otuIter)->numMembers; // sufficient prior to fastidious grafting

            }
        }
        std::cout << "Heavy swarms: " << (numOtus - numLightOtus) << ", with " << (numAmplicons - numAmplLightOtus) << " amplicons" << std::endl;
        std::cout << "Light swarms: " << numLightOtus << ", with " << numAmplLightOtus << " amplicons" << std::endl;

        if ((numLightOtus == 0) || (numLightOtus == numOtus)) {
            std::cout << "Fastidious: Only light or only heavy OTUs. No further action." << std::endl;
        } else {
            graftOtus(maxSize, numOtusAdjusted, pools, otus, sc);
        }

    }

    /* (c) Generating results */
    std::vector<Otu*> flattened(numOtus);
    auto iter = flattened.begin();
    for (auto p = 0; p < pools.numPools(); p++) {

        iter = std::move(otus[p].begin(), otus[p].end(), iter);
        otus[p] = std::vector<Otu*>();

    }

    std::sort(flattened.begin(), flattened.end(), CompareOtusSeedAbund());

    if (sc.outInternals) outputInternalStructures(sc.oFileInternals, pools, flattened, sc);
    if (sc.outOtus) {
        (sc.outMothur) ?
          outputOtusMothur(sc.oFileOtus, pools, flattened, sc.threshold, numOtusAdjusted, sc.sepMothur, sc.sepMothurOtu, sc.sepAbundance)
        : outputOtus(sc.oFileOtus, pools, flattened, sc.sepOtus, sc.sepAbundance);
    }
    if (sc.outStatistics) outputStatistics(sc.oFileStatistics, pools, flattened, sc.sepStatistics);
    if (sc.outSeeds) outputSeeds(sc.oFileSeeds, pools, flattened, sc.sepAbundance);
    if (sc.outUclust) outputUclust(sc.oFileUclust, pools, flattened, sc);

    std::cout << "Number of swarms: " << numOtusAdjusted << std::endl;
    std::cout << "Largest swarm: " << maxSize << std::endl;
    std::cout << "Max generations: " << maxGen << std::endl;


    /* (d) Cleaning up */
    for (auto otuIter = flattened.begin(); otuIter != flattened.end(); otuIter++) {
        delete *otuIter;
    }

}


void SwarmClustering::cluster(const AmpliconPools& pools, const SwarmConfig& sc) {

    /* (a) Mandatory (first) clustering phase of swarm */
    // determine OTUs by exploring all pools
    std::vector<std::vector<Otu*>> otus(pools.numPools());
    std::thread explorers[sc.numExplorers];
    auto fun = (sc.numThreadsPerExplorer == 1) ? &SegmentFilter::swarmFilterDirectly : &SegmentFilter::swarmFilter;
    unsigned long r = 0;
    for (; r + sc.numExplorers <= pools.numPools(); r += sc.numExplorers) {

        for (unsigned long e = 0; e < sc.numExplorers; e++) {
            explorers[e] = std::thread(fun, std::ref(*(pools.get(r + e))), std::ref(otus[r + e]), std::ref(sc));
        }
        for (unsigned long e = 0; e < sc.numExplorers; e++) {
            explorers[e].join();
        }

    }

    for (unsigned long e = 0; e < pools.numPools() % sc.numExplorers; e++) {
        explorers[e] = std::thread(fun, std::ref(*(pools.get(r + e))), std::ref(otus[r + e]), std::ref(sc));
    }
    for (unsigned long e = 0; e < pools.numPools() % sc.numExplorers; e++) {
        explorers[e].join();
    }

    processOtus(pools, otus, sc);

}


void SwarmClustering::dereplicate(const AmpliconPools& pools, const SwarmConfig& sc) {

    struct lessCharArray {
        bool operator()(const char* lhs, const char* rhs) const {
            return strcmp(lhs, rhs) < 0;
        }
    };

    std::vector<Otu*> otus;

    for (numSeqs_t p = 0; p < pools.numPools(); p++) {

        AmpliconCollection* ac = pools.get(p);
        std::map<const char*, std::vector<Amplicon*>, lessCharArray> groups;

        for (auto iter = ac->begin(); iter != ac->end(); iter++) {
            groups[iter->seq].push_back(iter);
        }

        for (auto& g : groups) {

            Otu* otu = new Otu();

            otu->numMembers = g.second.size();
            otu->numUniqueSequences = otu->numMembers;
            otu->members = new OtuEntry[otu->numMembers];

            for (auto i = 0; i < otu->numMembers; i++) {

                otu->members[i] = OtuEntry(g.second[i], g.second[0], 0, 0);
                otu->mass += g.second[i]->abundance;

            }

            otus.push_back(otu);

        }

    }

    std::sort(otus.begin(), otus.end(), CompareOtusMass());
    outputDereplicate(pools, otus, sc);

    numSeqs_t maxSize = 0;
    for (auto& o : otus) {
        maxSize = std::max(maxSize, o->numMembers);
    }

    std::cout << "Number of swarms: " << otus.size() << std::endl;
    std::cout << "Largest swarm: " << maxSize << std::endl;
    std::cout << "Max generations: 0" << std::endl;

}


void SwarmClustering::outputInternalStructures(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const SwarmConfig& sc) {

    std::ofstream oStream(oFile);
    std::stringstream sStream;

    lenSeqs_t width = 0;
    for (lenSeqs_t p = 0; p < pools.numPools(); p++) {
        width = std::max(width, pools.get(p)->maxLen());
    }

    lenSeqs_t M[sc.useScore ? 1 : width];
    val_t D[sc.useScore? width : 1];
    val_t P[sc.useScore? width : 1];
    lenSeqs_t cntDiffs[sc.useScore? width : 1];
    lenSeqs_t cntDiffsP[sc.useScore? width : 1];

    Otu* otu = 0;
    numSeqs_t otuId = 0;

    for (auto i = 0; i < otus.size(); i++) {

        otu = otus[i];

        if (!otu->attached()) {

            otuId++;

            for (auto otuIter = otu; otuIter != 0; otuIter = otuIter->nextGraftedOtu) {

                for (auto memberIter = otuIter->members; memberIter != otuIter->members + otuIter->numMembers; memberIter++) {

                    if (otuIter->graftChild == memberIter->member) {

                        lenSeqs_t dist = (sc.useScore) ? Verification::computeGotohLengthAwareEarlyRow(otuIter->graftChild->seq, otuIter->graftChild->len, otuIter->graftParent->member->seq, otuIter->graftParent->member->len, sc.fastidiousThreshold, sc.scoring, D, P, cntDiffs, cntDiffsP)
                                                       : Verification::computeLengthAwareRow(otuIter->graftChild->seq, otuIter->graftChild->len, otuIter->graftParent->member->seq, otuIter->graftParent->member->len, sc.fastidiousThreshold, M);
                        sStream << otuIter->graftParent->member->id << sc.sepInternals << otuIter->graftChild->id << sc.sepInternals << dist << sc.sepInternals << otuId << sc.sepInternals << (otuIter->graftParent->gen + 1) << std::endl;

                    }

                    if (memberIter->gen != 0) {

                        sStream << memberIter->parent->id << sc.sepInternals << memberIter->member->id << sc.sepInternals << memberIter->parentDist << sc.sepInternals << otuId << sc.sepInternals << memberIter->gen << std::endl;
                        oStream << sStream.rdbuf();
                        sStream.str(std::string());

                    }

                }

            }

        }

    }

    oStream.close();

}

void SwarmClustering::outputOtus(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const char sep, const std::string sepAbundance) {

    std::ofstream oStream(oFile);
    std::stringstream sStream;

    Otu* otu = 0;

    for (auto i = 0; i < otus.size(); i++) {

        otu = otus[i];

        if (!otu->attached()) {

            for (auto otuIter = otu; otuIter != 0; otuIter = otuIter->nextGraftedOtu) {

                if (otuIter != otu) {
                    sStream << sep;
                }

                sStream << otuIter->seed()->id << sepAbundance << otuIter->seedAbundance();

                for (auto memberIter = otuIter->members + 1; memberIter != otuIter->members + otuIter->numMembers; memberIter++) {
                    sStream << sep << memberIter->member->id << sepAbundance << memberIter->member->abundance;
                }

                sStream << std::flush;
                oStream << sStream.rdbuf();
                sStream.str(std::string());

            }

            oStream << std::endl;

        }

    }

    oStream.close();

}

void SwarmClustering::outputOtusMothur(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const lenSeqs_t threshold, const numSeqs_t numOtusAdjusted, const char sep, const std::string sepOtu, const std::string sepAbundance) {

    std::ofstream oStream(oFile);
    std::stringstream sStream;

    Otu* otu = 0;

    oStream << "swarm_" << threshold << "\t" << numOtusAdjusted;

    for (auto i = 0; i < otus.size(); i++) {

        otu = otus[i];

        if (!otu->attached()) {

            for (auto otuIter = otu; otuIter != 0; otuIter = otuIter->nextGraftedOtu) {

                sStream << ((otuIter == otu) ? sepOtu : std::string(1, sep)) << otuIter->seed()->id << sepAbundance << otuIter->seed()->abundance;

                for (auto memberIter = otuIter->members + 1; memberIter != otuIter->members + otuIter->numMembers; memberIter++) {
                    sStream << sep << memberIter->member->id << sepAbundance << memberIter->member->abundance;
                }

                oStream << sStream.rdbuf();
                sStream.str(std::string());

            }

        }

    }

    oStream << std::endl;

    oStream.close();

}

void SwarmClustering::outputStatistics(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const char sep) {

    std::ofstream oStream(oFile);
    std::stringstream sStream;

    Otu* otu = 0;

    for (auto i = 0; i < otus.size(); i++) {

        otu = otus[i];

        if (!otu->attached()) {

            auto maxGenRad = otu->maxGenRad();
            sStream << otu->numUniqueSequences << sep << otu->mass << sep << otu->seed()->id << sep << otu->seedAbundance() << sep << otu->numTotalSingletons() << sep << maxGenRad.first << sep << maxGenRad.second << std::endl;
            oStream << sStream.rdbuf();
            sStream.str(std::string());

        }

    }

    oStream.close();

}

void SwarmClustering::outputSeeds(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const std::string sepAbundance) {

    std::ofstream oStream(oFile);
    std::stringstream sStream;

    Otu* otu = 0;

    for (auto i = 0; i < otus.size(); i++) {

        otu = otus[i];

        if (!otu->attached()) {

            sStream << ">" << otu->seed()->id << sepAbundance << otu->mass << std::endl << otu->seed()->seq << std::endl;
            oStream << sStream.rdbuf();
            sStream.str(std::string());

        }

    }

    oStream.close();

}

void SwarmClustering::outputUclust(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const SwarmConfig& sc) {

    std::ofstream oStream(oFile);
    std::stringstream sStream;
    sStream << std::fixed << std::setprecision(1);

    AmpliconCollection* ac = 0;
    Otu* otu = 0;
    numSeqs_t otuId = 0;

    ac = pools.get(pools.numPools() - 1);
    lenSeqs_t maxLen = 0;
    for (auto iter = ac->begin(); iter != ac->end(); iter++) {
        maxLen = std::max(maxLen, iter->len);
    }
    val_t D[maxLen + 1];
    val_t P[maxLen + 1];
    char BT[(maxLen + 1) * (maxLen + 1)];

    for (auto i = 0; i < otus.size(); i++) {

        otu = otus[i];

        if (!otu->attached()) {

            auto& seed = *otu->seed();

            sStream << 'C' << sc.sepUclust << otuId << sc.sepUclust << otu->numTotalMembers() << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust
                    << seed.id << sc.sepAbundance << seed.abundance << sc.sepUclust << '*' << '\n';
            sStream << 'S' << sc.sepUclust << otuId << sc.sepUclust << seed.len << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust
                    << seed.id << sc.sepAbundance << seed.abundance << sc.sepUclust << '*' << '\n';
            oStream << sStream.rdbuf();
            sStream.str(std::string());

            for (auto otuIter = otu; otuIter != 0; otuIter = otuIter->nextGraftedOtu) {

                for (auto memberIter = otuIter->members + 1; memberIter != otuIter->members + otuIter->numMembers; memberIter++) {

                    auto& member = *memberIter->member;
                    auto ai = Verification::computeGotohCigarRow1(seed.seq, seed.len, member.seq, member.len, sc.scoring, D, P, BT);

                    sStream << 'H' << sc.sepUclust << otuId << sc.sepUclust << member.len << sc.sepUclust << (100.0 * (ai.length - ai.numDiffs) / ai.length) << sc.sepUclust << '+' << sc.sepUclust << '0' << sc.sepUclust << '0' << sc.sepUclust << ((ai.numDiffs == 0) ? "=" : ai.cigar) << sc.sepUclust
                            << member.id << sc.sepAbundance << member.abundance << sc.sepUclust
                            << seed.id << sc.sepAbundance << seed.abundance << '\n';
                    oStream << sStream.rdbuf();
                    sStream.str(std::string());

                }

            }

            otuId++;

        }

    }

    oStream.close();

}

void SwarmClustering::outputDereplicate(const AmpliconPools& pools, const std::vector<Otu*>& otus, const SwarmConfig& sc) {

    std::ofstream oStreamInternals, oStreamOtus, oStreamStatistics, oStreamSeeds, oStreamUclust;
    std::stringstream sStreamInternals, sStreamOtus, sStreamStatistics, sStreamSeeds, sStreamUclust;
    sStreamUclust << std::fixed << std::setprecision(1);

    if (sc.outInternals) oStreamInternals.open(sc.oFileInternals);
    if (sc.outOtus) oStreamOtus.open(sc.oFileOtus);
    if (sc.outStatistics) oStreamStatistics.open(sc.oFileStatistics);
    if (sc.outSeeds) oStreamSeeds.open(sc.oFileSeeds);
    if (sc.outUclust) oStreamUclust.open(sc.oFileUclust);

    if (sc.outOtus && sc.outMothur) oStreamOtus << "swarm_" << sc.threshold << "\t" << otus.size();

    for (auto i = 0; i < otus.size(); i++) {

        Otu& otu = *(otus[i]);

        if (sc.outInternals) {

            for (auto memberIter = otu.members + 1; memberIter != otu.members + otu.numMembers; memberIter++) {

                sStreamInternals << otu.seed()->id << sc.sepInternals << memberIter->member->id << sc.sepInternals << 0 << sc.sepInternals << (i + 1) << sc.sepInternals << 0 << std::endl;
                oStreamInternals << sStreamInternals.rdbuf();
                sStreamInternals.str(std::string());

            }

        }

        if (sc.outOtus) {

            if (sc.outMothur) {

                sStreamOtus << sc.sepMothurOtu << otu.seed()->id << sc.sepAbundance << otu.seedAbundance();

                for (auto memberIter = otu.members + 1; memberIter != otu.members + otu.numMembers; memberIter++) {
                    sStreamOtus << sc.sepMothur << memberIter->member->id << sc.sepAbundance << memberIter->member->abundance;
                }

                oStreamOtus << sStreamOtus.rdbuf();
                sStreamOtus.str(std::string());

            } else {

                sStreamOtus << otu.seed()->id << sc.sepAbundance << otu.seedAbundance();

                for (auto memberIter = otu.members + 1; memberIter != otu.members + otu.numMembers; memberIter++) {
                    sStreamOtus << sc.sepOtus << memberIter->member->id << sc.sepAbundance << memberIter->member->abundance;
                }

                sStreamOtus << std::endl;
                oStreamOtus << sStreamOtus.rdbuf();
                sStreamOtus.str(std::string());

            }

        }

        if (sc.outStatistics) {

            sStreamStatistics << otu.numUniqueSequences << sc.sepStatistics << otu.mass << sc.sepStatistics << otu.seed()->id << sc.sepStatistics << otu.seedAbundance() << sc.sepStatistics << otu.numSingletons() << sc.sepStatistics << 0 << sc.sepStatistics << 0 << std::endl;
            oStreamStatistics << sStreamStatistics.rdbuf();
            sStreamStatistics.str(std::string());

        }

        if (sc.outSeeds) {

            sStreamSeeds << ">" << otu.seed()->id << sc.sepAbundance << otu.mass << std::endl << otu.seed()->seq << std::endl;
            oStreamSeeds << sStreamSeeds.rdbuf();
            sStreamSeeds.str(std::string());

        }

        if (sc.outUclust) {

            auto& seed = *otu.seed();

            sStreamUclust << 'C' << sc.sepUclust << i << sc.sepUclust << otu.numMembers << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust
                          << seed.id << sc.sepAbundance << seed.abundance << sc.sepUclust << '*' << '\n';
            sStreamUclust << 'S' << sc.sepUclust << i << sc.sepUclust << seed.len << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust << '*' << sc.sepUclust
                          << seed.id << sc.sepAbundance << seed.abundance << sc.sepUclust << '*' << '\n';
            oStreamUclust << sStreamUclust.rdbuf() << std::flush;
            sStreamUclust.str(std::string());

            for (auto memberIter = otu.members + 1; memberIter != otu.members + otu.numMembers; memberIter++) {

                sStreamUclust << 'H' << sc.sepUclust << i << sc.sepUclust << memberIter->member->len << sc.sepUclust << "100.0" << sc.sepUclust << '+' << sc.sepUclust << '0' << sc.sepUclust << '0' << sc.sepUclust << '=' << sc.sepUclust
                            << memberIter->member->id << sc.sepAbundance << memberIter->member->abundance << sc.sepUclust
                            << seed.id << sc.sepAbundance << seed.abundance << '\n';

                oStreamUclust << sStreamUclust.rdbuf() << std::flush;
                sStreamUclust.str(std::string());

            }


        }

    }

    if (sc.outOtus && sc.outMothur) oStreamOtus << std::endl;

    if (sc.outInternals) oStreamInternals.close();
    if (sc.outOtus) oStreamOtus.close();
    if (sc.outStatistics) oStreamStatistics.close();
    if (sc.outSeeds) oStreamSeeds.close();
    if (sc.outUclust) oStreamUclust.close();

}

}