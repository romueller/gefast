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

#include "../include/SwarmClustering.hpp"

#include <fstream>
#include <set>
#include <sstream>
#include <thread>
#include <unordered_set>

namespace SCT_PJ {

void SwarmClustering::explorePool(const AmpliconCollection& ac, Matches& matches, std::vector<Otu*>& otus, const SwarmConfig& sc) {

    // determine order of amplicons based on abundance (descending) without invalidating the integer (position) ids of the amplicons
    std::vector<numSeqs_t> index(ac.size());
    std::iota(std::begin(index), std::end(index), 0);
    std::sort(index.begin(), index.end(), CompareIndicesAbund(ac));

    Otu* curOtu = 0;
//    numSeqs_t numOtu = 0;
    std::vector<bool> visited(ac.size(), false); // visited amplicons are already included in an OTU

    OtuEntry curSeed, newSeed;
    bool unique;
    std::unordered_set<std::string> nonUniques;
    std::vector<std::pair<numSeqs_t, lenSeqs_t>> next;
    lenSeqs_t lastGen;
    numSeqs_t pos;

    // open new OTU for the amplicon with the highest abundance that is not yet included in an OTU
    for (auto seedIter = index.begin(); seedIter != index.end(); seedIter++) {

        if (!visited[*seedIter]) {

            /* (a) Initialise new OTU with seed */
            curOtu = new Otu(/*++numOtu, */*seedIter, ac[*seedIter].abundance);

            newSeed.id = *seedIter;
            newSeed.parentId = newSeed.id;
            newSeed.parentDist = 0;
            newSeed.gen = 0;
            newSeed.rad = 0;
            curOtu->members.push_back(newSeed);

            visited[*seedIter] = true;
            nonUniques.clear();

            lastGen = 0;


            /* (b) BFS through 'match space' */
            pos = 0;
            while (pos < curOtu->members.size()) { // expand current OTU until no further similar amplicons can be added

                if (lastGen != curOtu->members[pos].gen) { // work through generation by decreasing abundance
                    std::sort(curOtu->members.begin() + pos, curOtu->members.end(), CompareOtuEntriesAbund(ac));
                }

                // get next OTU (sub)seed
                curSeed = curOtu->members[pos];

                unique = true;

                // update OTU information
                curOtu->mass += ac[curSeed.id].abundance;
                curOtu->numSingletons += (ac[curSeed.id].abundance == 1);

                if (curSeed.gen > curOtu->maxGen) curOtu->maxGen = curSeed.gen;
                if (curSeed.rad > curOtu->maxRad) curOtu->maxRad = curSeed.rad;

                // Consider yet unseen (unvisited) amplicons to continue the exploration.
                // An amplicon is marked as 'visited' as soon as it occurs the first time as matching partner
                // in order to prevent the algorithm from queueing it more than once coming from different amplicons.
                next = matches.getMatchesOfAugmented(curSeed.id);
                for (auto matchIter = next.begin(); matchIter != next.end(); matchIter++) {

                    unique &= (matchIter->second != 0);

                    if (!visited[matchIter->first] && (sc.noOtuBreaking || ac[matchIter->first].abundance <= ac[curSeed.id].abundance)) {

                        newSeed.id = matchIter->first;
                        newSeed.parentId = curSeed.id;
                        newSeed.parentDist = matchIter->second;
                        newSeed.gen = curSeed.gen + 1;
                        newSeed.rad = curSeed.rad + matchIter->second;
                        curOtu->members.push_back(newSeed);
                        visited[matchIter->first] = true;

                    }
                }

                // unique sequences contribute when they occur, non-unique sequences only at their first occurrence
                // and when dereplicating each contributes (numUniqueSequences used to count the multiplicity of the sequence)
                unique = unique || sc.dereplicate || nonUniques.insert(ac[curSeed.id].seq).second;
                curOtu->numUniqueSequences += unique;

                lastGen = curSeed.gen;
                pos++;

            }

            /* (c) Close the no longer extendable OTU */
            otus.push_back(curOtu);

        }

    }

}


void SwarmClustering::fastidiousIndexOtu(RollingIndices<InvertedIndexFastidious>& indices, std::unordered_map<lenSeqs_t, SegmentFilter::Segments>& segmentsArchive, const AmpliconCollection& ac, Otu& otu, std::vector<GraftCandidate>& graftCands, const SwarmConfig& sc) {

    lenSeqs_t seqLen;

    for (auto memberIter = otu.members.begin(); memberIter != otu.members.end(); memberIter++) {

        seqLen = ac[memberIter->id].seq.length();
        SegmentFilter::Segments& segments = segmentsArchive[seqLen];

        if (segments.size() == 0) {

            indices.roll(seqLen);
            segments = SegmentFilter::Segments(2 * sc.threshold + sc.extraSegs);
            SegmentFilter::selectSegments(segments, seqLen, 2 * sc.threshold, sc.extraSegs);

        }

        for (lenSeqs_t i = 0; i < 2 * sc.threshold + sc.extraSegs; i++) {
            indices.getIndex(seqLen, i).add(ac[memberIter->id].seq.substr(segments[i].first, segments[i].second), &(*memberIter));
        }

        graftCands[memberIter->id].childOtu = &otu;
        graftCands[memberIter->id].childMember = &(*memberIter);

    }

}

inline bool compareCandidates(const Amplicon& newCand, const Amplicon& oldCand) {
    return (oldCand.abundance < newCand.abundance)
#if INPUT_RANK
           || ((oldCand.abundance == newCand.abundance) && (oldCand.rank > newCand.rank))
#endif
            ;
}

void SwarmClustering::verifyFastidious(const AmpliconCollection& acOtus, const AmpliconCollection& acIndices, std::vector<GraftCandidate>& graftCands, Buffer<CandidateFastidious>& buf, lenSeqs_t t, std::mutex& mtx) {

    CandidateFastidious c;
    Buffer<CandidateFastidious> localBuffer;
    lenSeqs_t M[std::max(acOtus.back().seq.length(), acIndices.back().seq.length()) + 1]; // reusable DP-matrix (wide enough for all possible calculations for this AmpliconCollection)

    while (!buf.isClosed() || buf.syncSize() > 0) {

        buf.syncSwapContents(localBuffer);

        while (localBuffer.size() > 0) {

            c = localBuffer.pop();

            for (auto childIter = c.children.begin(); childIter != c.children.end(); childIter++) {

                std::unique_lock<std::mutex> lock(mtx);
                if ((graftCands[*childIter].parentOtu == 0) || compareCandidates(acOtus[c.parent], acOtus[graftCands[*childIter].parentMember->id])) {

                    lock.unlock();
                    if (Verification::computeLengthAwareRow(acOtus[c.parent].seq, acIndices[*childIter].seq, t, M) <= t) {

                        lock.lock();
                        if (((graftCands[*childIter].parentOtu == 0) || compareCandidates(acOtus[c.parent], acOtus[graftCands[*childIter].parentMember->id]))) {

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

void SwarmClustering::fastidiousCheckOtus(RotatingBuffers<CandidateFastidious>& cbs, const std::vector<Otu*>& otus, const AmpliconCollection& acOtus, RollingIndices<InvertedIndexFastidious>& indices, const AmpliconCollection& acIndices, std::vector<GraftCandidate>& graftCands, std::mutex& mtx, const SwarmConfig& sc) {

    std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<SegmentFilter::Substrings>>> substrsArchive;
    std::vector<OtuEntry*> candMembers;
    std::unordered_map<numSeqs_t, lenSeqs_t> candCnts;
    lenSeqs_t M[std::max(acOtus.back().seq.length(), acIndices.back().seq.length()) + 1];
    lenSeqs_t seqLen;

    std::vector<CandidateFastidious> localCands;

    for (auto otuIter = otus.begin(); otuIter != otus.end(); otuIter++) {

        if ((*otuIter)->mass >= sc.boundary) { // for each heavy OTU of the pool ...

            for (auto memberIter = (*otuIter)->members.begin(); memberIter != (*otuIter)->members.end(); memberIter++) { // ... consider every amplicon in the OTU and ...

                seqLen = acOtus[memberIter->id].seq.length();

                std::unordered_map<lenSeqs_t, std::vector<SegmentFilter::Substrings>>& substrs = substrsArchive[seqLen];

                // on reaching new length group, open new inverted indices
                if (substrs.empty()) {

                    // ... and determine position information shared by all amplicons of this length
                    for (lenSeqs_t partnerLen = (seqLen > 2 * sc.threshold) * (seqLen - 2 * sc.threshold); partnerLen <= seqLen + 2 * sc.threshold; partnerLen++) {

                        std::vector<SegmentFilter::Substrings>& vec = substrs[partnerLen];
                        for (lenSeqs_t segmentIndex = 0; segmentIndex < 2 * sc.threshold + sc.extraSegs; segmentIndex++) {
                            if (partnerLen <= seqLen) {
                                vec.push_back(SegmentFilter::selectSubstrs(seqLen, partnerLen, segmentIndex, 2 * sc.threshold, sc.extraSegs));
                            } else {
                                vec.push_back(SegmentFilter::selectSubstrsBackward(seqLen, partnerLen, segmentIndex, 2 * sc.threshold, sc.extraSegs));
                            }
                        }

                    }

                }

                localCands.push_back(CandidateFastidious(memberIter->id, *otuIter, &(*memberIter)));

                for (lenSeqs_t len = (seqLen > 2 * sc.threshold) * (seqLen - 2 * sc.threshold); len <= seqLen + 2 * sc.threshold; len++) { // ... search for graft candidates among the amplicons in light OTUs

                    for (lenSeqs_t i = 0; i < 2 * sc.threshold + sc.extraSegs; i++) {//... and apply segment filter for each segment

                        SegmentFilter::Substrings& subs = substrs[len][i];
                        InvertedIndexFastidious& inv = indices.getIndex(len, i);

                        for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++) {

                            candMembers = inv.getLabelsOf(std::string(acOtus[memberIter->id].seq, substrPos, subs.len));

                            for (auto candIter = candMembers.begin(); candIter != candMembers.end(); candIter++) {
                                candCnts[(*candIter)->id]++;
                            }

                        }

                    }

                    // general pigeonhole principle: for being a candidate, at least sc.extraSegs segments have to be matched
                    for (auto candIter = candCnts.begin(); candIter != candCnts.end(); candIter++) {

                        if (candIter->second >= sc.extraSegs) {
                            localCands.back().children.push_back(candIter->first);
                        }

                    }

                    candCnts = std::unordered_map<numSeqs_t, lenSeqs_t>();

                }

                cbs.push(localCands);
                localCands = std::vector<CandidateFastidious>();

            }

        }

    }

}

void SwarmClustering::checkAndVerify(const std::vector<Otu*>& otus, const AmpliconCollection& acOtus, RollingIndices<InvertedIndexFastidious>& indices, const AmpliconCollection& acIndices, std::vector<GraftCandidate>& graftCands, std::mutex& graftCandsMtx, const SwarmConfig& sc) {

    RotatingBuffers<CandidateFastidious> cbs = RotatingBuffers<CandidateFastidious>(sc.numVerifiersPerChecker);
    std::thread verifierThreads[sc.numVerifiersPerChecker];

    for (unsigned long v = 0; v < sc.numVerifiersPerChecker; v++) {
        verifierThreads[v] = std::thread(&SwarmClustering::verifyFastidious, std::ref(acOtus), std::ref(acIndices), std::ref(graftCands), std::ref(cbs.getBuffer(v)), 2 * sc.threshold, std::ref(graftCandsMtx));
    }

    fastidiousCheckOtus(cbs, otus, acOtus, indices, acIndices, graftCands, graftCandsMtx, sc);
    cbs.close();

    for (unsigned long v = 0; v < sc.numVerifiersPerChecker; v++) {
        verifierThreads[v].join();
    }

}

void SwarmClustering::determineGrafts(const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, std::vector<GraftCandidate>& allGraftCands, const numSeqs_t p, std::mutex& allGraftCandsMtx, const SwarmConfig& sc) {

    AmpliconCollection* ac = pools.get(p);
    RollingIndices<InvertedIndexFastidious> indices = RollingIndices<InvertedIndexFastidious>(4 * sc.threshold + 1, 2 * sc.threshold + sc.extraSegs, true, false);
    std::vector<GraftCandidate> graftCands = std::vector<GraftCandidate>(ac->size()); // initially, graft candidates for all amplicons of the pool are "empty"
    std::unordered_map<lenSeqs_t, SegmentFilter::Segments> segmentsArchive = std::unordered_map<lenSeqs_t, SegmentFilter::Segments>();
    std::mutex graftCandsMtx;

    // a) Index amplicons of all light OTUs of the current pool
    for (auto otuIter = otus[p].begin(); otuIter != otus[p].end(); otuIter++) {

        if ((*otuIter)->mass < sc.boundary) {
            fastidiousIndexOtu(indices, segmentsArchive, *ac, *(*otuIter), graftCands, sc);
        }

    }

    // b) Search with amplicons of all heavy OTUs of current and directly neighbouring pools
#if FASTIDIOUS_PARALLEL_CHECK

    switch (sc.fastidiousCheckingMode) {

        case 0: {

            if (p > 0) {
                checkAndVerify(otus[p - 1], *(pools.get(p - 1)), indices, *ac, graftCands, graftCandsMtx, sc);
            }

            checkAndVerify(otus[p], *ac, indices, *ac, graftCands, graftCandsMtx, sc);

            if (p < pools.numPools() - 1) {
                checkAndVerify(otus[p + 1], *(pools.get(p + 1)), indices, *ac, graftCands, graftCandsMtx, sc);
            }

            break;

        }

        case 1: {

            std::thread t(&SwarmClustering::checkAndVerify, std::ref(otus[p]), std::ref(*ac), std::ref(indices), std::ref(*ac), std::ref(graftCands), std::ref(graftCandsMtx), std::ref(sc));

            if (p > 0) {
                checkAndVerify(otus[p - 1], *(pools.get(p - 1)), indices, *ac, graftCands, graftCandsMtx, sc);
            }

            if (p < pools.numPools() - 1) {
                checkAndVerify(otus[p + 1], *(pools.get(p + 1)), indices, *ac, graftCands, graftCandsMtx, sc);
            }

            t.join();

            break;

        }

        default: {

                std::thread pred, succ;

                if (p > 0) {
                    pred = std::thread(&SwarmClustering::checkAndVerify, std::ref(otus[p - 1]), std::ref(*(pools.get(p - 1))), std::ref(indices), std::ref(*ac), std::ref(graftCands), std::ref(graftCandsMtx), std::ref(sc));
                }

                if (p < pools.numPools() - 1) {
                    succ = std::thread(&SwarmClustering::checkAndVerify, std::ref(otus[p + 1]), std::ref(*(pools.get(p + 1))), std::ref(indices), std::ref(*ac), std::ref(graftCands), std::ref(graftCandsMtx), std::ref(sc));
                }

                checkAndVerify(otus[p], *ac, indices, *ac, graftCands, graftCandsMtx, sc);

                if (p > 0) {
                    pred.join();
                }
                if (p < pools.numPools() - 1) {
                    succ.join();
                }

            break;

        }

    }

#else

    if (p > 0) {
        checkAndVerify(otus[p - 1], *(pools.get(p - 1)), indices, *ac, graftCands, graftCandsMtx, sc);
    }

    checkAndVerify(otus[p], *ac, indices, *ac, graftCands, graftCandsMtx, sc);

    if (p < pools.numPools() - 1) {
        checkAndVerify(otus[p + 1], *(pools.get(p + 1)), indices, *ac, graftCands, graftCandsMtx, sc);
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
    std::sort(allGraftCands.begin(), allGraftCands.end(), CompareGraftCandidatesAbund(pools));
    Otu* parentOtu = 0;
    Otu* childOtu = 0;
    numSeqs_t numGrafts = 0;
    for (auto graftIter = allGraftCands.begin(); graftIter != allGraftCands.end(); graftIter++) {

        if (!(graftIter->childOtu->attached)) {

            parentOtu = graftIter->parentOtu;
            childOtu = graftIter->childOtu;

            // "attach the see of the light swarm to the tail of the heavy swarm"
            // OTU entries are moved unchanged (entry of seed of attached swarm stays 'incomplete', grafting 'link' is not recorded)
            parentOtu->members.reserve(parentOtu->members.size() + childOtu->members.size());
            std::move(std::begin(childOtu->members), std::end(childOtu->members), std::back_inserter(parentOtu->members));
            childOtu->members.clear();

            // update stats
            if (parentOtu->members.size() > maxSize) {
                maxSize = parentOtu->members.size();
            }

            parentOtu->numUniqueSequences += childOtu->numUniqueSequences;
            parentOtu->numSingletons += childOtu->numSingletons;
            parentOtu->mass += childOtu->mass;
            // maximum generation / radius are untouched

            childOtu->attached = true;
            numGrafts++;
            numOtus--;

        }

    }

    std::cout << "Made " << numGrafts << " grafts." << std::endl;

}


void SwarmClustering::prepareGraftInfos(const numSeqs_t poolSize, const std::vector<Otu*>& otus, std::vector<GraftCandidate>& curGraftCands, std::vector<GraftCandidate>& nextGraftCands, const SwarmConfig& sc) {

    nextGraftCands = std::vector<GraftCandidate>(0);
    curGraftCands = std::vector<GraftCandidate>(poolSize);

    for (auto otuIter = otus.begin(); otuIter != otus.end(); otuIter++) {

        for (auto memberIter = (*otuIter)->members.begin(); memberIter != (*otuIter)->members.end(); memberIter++) {

            curGraftCands[memberIter->id].childOtu = *otuIter;
            curGraftCands[memberIter->id].childMember = &(*memberIter);

        }

    }

}

AmpliconCollection::iterator SwarmClustering::shiftIndexWindow(RollingIndices<InvertedIndexFastidious2>& indices, const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, const numSeqs_t poolIndex, const AmpliconCollection::iterator indexIter,
                                                               std::vector<GraftCandidate>& curGraftCands, std::vector<GraftCandidate>& nextGraftCands, const lenSeqs_t len, const bool forerunning, const SwarmConfig& sc) {

    // Advance indexIter in the same pool and update the collection of inverted indices
    AmpliconCollection* ac = pools.get(poolIndex + forerunning);
    auto newIter = indexIter;
    auto a = std::distance(ac->begin(), newIter);
    lenSeqs_t seqLen = 0;
    SegmentFilter::Segments segments(2 * sc.threshold + sc.extraSegs);

    for (; newIter != ac->end() && newIter->seq.length() <= len + 2 * sc.threshold; newIter++, a++) {

        if (curGraftCands[a].childOtu->mass >= sc.boundary) {

            if (newIter->seq.length() != seqLen) {

                seqLen = newIter->seq.length();
                indices.roll(seqLen);

                SegmentFilter::selectSegments(segments, seqLen, 2 * sc.threshold, sc.extraSegs);

            }

            for (lenSeqs_t i = 0; i < 2 * sc.threshold + sc.extraSegs; i++) {
                indices.getIndex(seqLen, i).add(newIter->seq.substr(segments[i].first, segments[i].second), std::make_pair(curGraftCands[a].childOtu, curGraftCands[a].childMember));
            }

        }

    }


    // Move indexIter to the next pool, if ...
    // ... the end position of indexIter (= first amplicon which is too long) in the current shifting round is in the next pool
    // ... and such a next pool exists
    // ... and indexIter would not advance two pools ahead of the amplicons currently filtered.
    // If then possible, advance indexIter and update collection of inverted indices as above.
    if (newIter == ac->end() && (poolIndex + 1 < pools.numPools()) && !forerunning) {

        ac = pools.get(poolIndex + 1);
        prepareGraftInfos(ac->size(), otus[poolIndex + 1], nextGraftCands, nextGraftCands, sc);

        for (newIter = ac->begin(), a = 0; newIter != ac->end() && newIter->seq.length() <= len + 2 * sc.threshold; newIter++, a++) {

            if (nextGraftCands[a].childOtu->mass >= sc.boundary) {

                if (newIter->seq.length() != seqLen) {

                    seqLen = newIter->seq.length();
                    indices.roll(seqLen);

                    SegmentFilter::selectSegments(segments, seqLen, 2 * sc.threshold, sc.extraSegs);

                }

                for (lenSeqs_t i = 0; i < 2 * sc.threshold + sc.extraSegs; i++) {
                    indices.getIndex(seqLen, i).add(newIter->seq.substr(segments[i].first, segments[i].second), std::make_pair(nextGraftCands[a].childOtu, nextGraftCands[a].childMember));
                }

            }

        }

    }

    return newIter;

}

void SwarmClustering::graftOtus2(numSeqs_t& maxSize, numSeqs_t& numOtus, const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, const SwarmConfig& sc) {

    RollingIndices<InvertedIndexFastidious2> indices(4 * sc.threshold + 1, 2 * sc.threshold + sc.extraSegs, true, false);
    std::unordered_map<lenSeqs_t, std::vector<SegmentFilter::Substrings>> substrs;//TODO? RollingIndices?
    std::vector<std::pair<Otu*, OtuEntry*>> candMembers;
    std::unordered_map<Otu*, std::unordered_map<OtuEntry*, lenSeqs_t>> candCnts;
    std::vector<GraftCandidate> curGraftCands, nextGraftCands, allGraftCands;
    lenSeqs_t M[pools.get(pools.numPools() - 1)->back().seq.length() + 1];

    // iterate over all amplicons in order of increasing length (as if they were stored in one big pool)
    numSeqs_t a = 0;
    numSeqs_t p = 0;
    AmpliconCollection* ac = pools.get(0);
    AmpliconCollection* last = pools.get(pools.numPools() - 1);
    lenSeqs_t seqLen = 0;
    std::vector<SegmentFilter::Substrings> tmp(2 * sc.threshold + sc.extraSegs);

    prepareGraftInfos(ac->size(), otus[0], curGraftCands, nextGraftCands, sc);
    auto indexIter = ac->begin();
    nextGraftCands.resize(0);//TODO? unnecessary?

    for (auto amplIter = ac->begin(); amplIter != last->end();) {

        if (curGraftCands[a].childOtu->mass < sc.boundary) {

            // arriving at a new length, continue indexing
            if (amplIter->seq.length() != seqLen) {

                seqLen = amplIter->seq.length();

                if (nextGraftCands.size() == 0) { // indexIter still in the same pool
                    indexIter = shiftIndexWindow(indices, pools, otus, p, indexIter, curGraftCands, nextGraftCands, seqLen, false, sc);
                } else { // indexIter already in the next pool
                    indexIter = shiftIndexWindow(indices, pools, otus, p, indexIter, nextGraftCands, nextGraftCands, seqLen, true, sc);
                }

                // determine substring information for this length
                for (lenSeqs_t partnerLen = (seqLen > 2 * sc.threshold) * (seqLen - 2 * sc.threshold); partnerLen <= seqLen + 2 * sc.threshold; partnerLen++) {

                    std::vector<SegmentFilter::Substrings>& vec = substrs[partnerLen];

                    for (lenSeqs_t segmentIndex = 0; segmentIndex < 2 * sc.threshold + sc.extraSegs; segmentIndex++) {
                        if (partnerLen <= seqLen) {
                            tmp[segmentIndex] = SegmentFilter::selectSubstrs(seqLen, partnerLen, segmentIndex, 2 * sc.threshold, sc.extraSegs);
                        } else {
                            tmp[segmentIndex] = SegmentFilter::selectSubstrsBackward(seqLen, partnerLen, segmentIndex, 2 * sc.threshold, sc.extraSegs);
                        }
                    }

                    vec.swap(tmp);
                    tmp.reserve(2 * sc.threshold + sc.extraSegs);

                }

            }

            // search for partners and verify them directly
            GraftCandidate& gc = curGraftCands[a];
            for (lenSeqs_t len = (seqLen > 2 * sc.threshold) * (seqLen - 2 * sc.threshold); len <= seqLen + 2 * sc.threshold; len++) { // ... search for grafting candidates among the indexed amplicons

                for (lenSeqs_t i = 0; i < 2 * sc.threshold + sc.extraSegs; i++) {//... and apply segment filter for each segment

                    SegmentFilter::Substrings& subs = substrs[len][i];
                    InvertedIndexFastidious2& inv = indices.getIndex(len, i);

                    for (auto substrPos = subs.first; substrPos <= subs.last; substrPos++) {

                        candMembers = inv.getLabelsOf(std::string(amplIter->seq, substrPos, subs.len));

                        for (auto candIter = candMembers.begin(); candIter != candMembers.end(); candIter++) {
                            candCnts[candIter->first][candIter->second]++;
                        }

                    }

                }

                // general pigeonhole principle: for being a candidate, at least sc.extraSegs segments have to be matched
                for (auto otuIter = candCnts.begin(); otuIter != candCnts.end(); otuIter++) {

                    for (auto candIter = otuIter->second.begin(); candIter != otuIter->second.end(); candIter++) {

                        if (candIter->second >= sc.extraSegs) {

                            Amplicon& cand = (*pools.get(otuIter->first->poolId))[candIter->first->id];

                            if ((gc.parentOtu == 0 || compareCandidates(cand, (*pools.get(gc.parentOtu->poolId))[gc.parentMember->id]))
                                    && (Verification::computeBoundedRow(amplIter->seq, cand.seq, 2 * sc.threshold, M) <= 2 * sc.threshold)) {

                                gc.parentOtu = otuIter->first;
                                gc.parentMember = candIter->first;

                            }

                        }

                    }

                }

                candCnts.clear();

            }

        }

        // iterator handling (includes jumping to next pool)
        amplIter++;
        a++;

        if (amplIter == ac->end() && ac != last) {

            // collect and sort all (actual = non-empty) graft candidates for the current pool
            auto candEnd = std::remove_if(
                    curGraftCands.begin(),
                    curGraftCands.end(),
                    [](GraftCandidate& gc) {
                        return gc.parentOtu == 0;
                    });

            allGraftCands.reserve(allGraftCands.size() + std::distance(curGraftCands.begin(), candEnd));
            std::move(curGraftCands.begin(), candEnd, std::back_inserter(allGraftCands));

            ac = pools.get(++p);
            amplIter = ac->begin();
            a = 0;

            if (nextGraftCands.size() != 0) {

                curGraftCands.swap(nextGraftCands);
                nextGraftCands = std::vector<GraftCandidate>(0);

            } else {
                prepareGraftInfos(ac->size(), otus[p], curGraftCands, nextGraftCands, sc);
            }

        }
    }


    // collect and sort all (actual = non-empty) graft candidates for the last pool
    auto candEnd = std::remove_if(
            curGraftCands.begin(),
            curGraftCands.end(),
            [](GraftCandidate& gc) {
                return gc.parentOtu == 0;
            });

    allGraftCands.reserve(allGraftCands.size() + std::distance(curGraftCands.begin(), candEnd));
    std::move(curGraftCands.begin(), candEnd, std::back_inserter(allGraftCands));


    // perform actual grafting
    std::cout << "Got " << allGraftCands.size() << " graft candidates" << std::endl;
    std::sort(allGraftCands.begin(), allGraftCands.end(), CompareGraftCandidatesAbund(pools));

    numSeqs_t numGrafts = 0;
    Otu* parentOtu = 0;
    Otu* childOtu = 0;
    for (auto graftIter = allGraftCands.begin(); graftIter != allGraftCands.end(); graftIter++) {

        if (!graftIter->childOtu->attached) {

            parentOtu = graftIter->parentOtu;
            childOtu = graftIter->childOtu;

            // "attach the see of the light swarm to the tail of the heavy swarm"
            // OTU entries are moved unchanged (entry of seed of attached swarm stays 'incomplete', grafting 'link' is not recorded)
            parentOtu->members.reserve(parentOtu->members.size() + childOtu->members.size());
            std::move(std::begin(childOtu->members), std::end(childOtu->members), std::back_inserter(parentOtu->members));
            childOtu->members.clear();

            // update stats
            if (parentOtu->members.size() > maxSize) {
                maxSize = parentOtu->members.size();
            }

            parentOtu->numUniqueSequences += childOtu->numUniqueSequences;
            parentOtu->numSingletons += childOtu->numSingletons;
            parentOtu->mass += childOtu->mass;
            // maximum generation / radius are untouched

            childOtu->attached = true;
            numGrafts++;
            numOtus--;

        }

    }

    std::cout << "Made " << numGrafts << " grafts" << std::endl;

}

void SwarmClustering::cluster(const AmpliconPools& pools, std::vector<Matches*>& allMatches, const SwarmConfig& sc) {

    /* (a) Mandatory (first) clustering phase of swarm */
    // determine OTUs by exploring all pools
    std::vector<std::vector<Otu*>> otus(pools.numPools());
    std::thread explorers[sc.numExplorers];
    unsigned long r = 0;
    for (; r + sc.numExplorers <= pools.numPools(); r += sc.numExplorers) {

        for (unsigned long e = 0; e < sc.numExplorers; e++) {
            explorers[e] = std::thread(&SwarmClustering::explorePool, std::ref(*(pools.get(r + e))), std::ref(*(allMatches[r + e])), std::ref(otus[r + e]), std::ref(sc));
        }
        for (unsigned long e = 0; e < sc.numExplorers; e++) {
            explorers[e].join();
        }

    }

    for (unsigned long e = 0; e < pools.numPools() % sc.numExplorers; e++) {
        explorers[e] = std::thread(&SwarmClustering::explorePool, std::ref(*(pools.get(r + e))), std::ref(*(allMatches[r + e])), std::ref(otus[r + e]), std::ref(sc));
    }
    for (unsigned long e = 0; e < pools.numPools() % sc.numExplorers; e++) {
        explorers[e].join();
    }

    // make OTU IDs unique over all pools (so far IDs start at 1 in each pool) (currently commented out)
    // add pool IDs and determine some overall statistics
    numSeqs_t numOtus = otus[0].size();
    numSeqs_t numOtusAdjusted = 0;
    numSeqs_t numAmplicons = pools.get(0)->size();
    numSeqs_t maxSize = 0;
    numSeqs_t maxGen = 0;

    for (numSeqs_t p = 1; p < pools.numPools(); p++) {

        for (auto otuIter = otus[p].begin(); otuIter != otus[p].end(); otuIter++) {

//            (*otuIter)->id += numOtus;
            (*otuIter)->poolId = p;

            if ((*otuIter)->members.size() > maxSize){
                maxSize = (*otuIter)->members.size();
            }

            if ((*otuIter)->maxGen > maxGen){
                maxGen = (*otuIter)->maxGen;
            }

        }

        numOtus += otus[p].size();
        numAmplicons += pools.get(p)->size();

    }

    numOtusAdjusted = numOtus;

#if !PRINT_INTERNAL_MODIFIED
    if (!sc.dereplicate && sc.outInternals) {

        std::vector<Otu*> flattened;
        flattened.reserve(numOtus);

        for (auto p = 0; p < pools.numPools(); p++) {
            flattened.insert(flattened.end(), otus[p].begin(), otus[p].end());
        }

        std::sort(flattened.begin(), flattened.end(), CompareOtusSeedAbund(pools));
        outputInternalStructures(sc.oFileInternals, pools, flattened, sc.sepInternals);

    }
#endif

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
                numAmplLightOtus += ((*otuIter)->mass < sc.boundary) * (*otuIter)->members.size();

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

    if (sc.dereplicate) {

        std::sort(flattened.begin(), flattened.end(), CompareOtusMass(pools));

        outputDereplicate(pools, flattened, sc);

    } else {

        std::sort(flattened.begin(), flattened.end(), CompareOtusSeedAbund(pools));

#if PRINT_INTERNAL_MODIFIED
        if (sc.outInternals) outputInternalStructures(sc.oFileInternals, pools, flattened, sc.sepInternals);
#endif
        if (sc.outOtus) outputOtus(sc.oFileOtus, pools, flattened, sc.sepOtus, sc.sepAbundance);
        if (sc.outStatistics) outputStatistics(sc.oFileStatistics, pools, flattened, sc.sepStatistics);
        if (sc.outSeeds) outputSeeds(sc.oFileSeeds, pools, flattened, sc.sepAbundance);

    }

    std::cout << "Number of swarms: " << numOtusAdjusted << std::endl;
    std::cout << "Largest swarm: " << maxSize << std::endl;
    std::cout << "Max generations: " << maxGen << std::endl;


    /* (d) Cleaning up */
    for (auto otuIter = flattened.begin(); otuIter != flattened.end(); otuIter++) {
        delete *otuIter;
    }

}


void SwarmClustering::outputInternalStructures(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const char sep) {

    std::ofstream oStream(oFile);
    std::stringstream sStream;

    AmpliconCollection* ac = 0;
    Otu* otu = 0;

    for (auto i = 0; i < otus.size(); i++) {

        otu = otus[i];
        ac = pools.get(otu->poolId);

        if (!otu->attached) {

            for (auto memberIter = otu->members.begin() + 1; memberIter != otu->members.end(); memberIter++) {

#if PRINT_INTERNAL_MODIFIED
                if (memberIter->id == memberIter->parentId) continue;
#endif

                sStream << (*ac)[memberIter->parentId].id << sep << (*ac)[memberIter->id].id << sep << memberIter->parentDist << sep << (i + 1) << sep << memberIter->gen << std::endl;
                oStream << sStream.rdbuf();
                sStream.str(std::string());

            }

        }

    }

    oStream.close();

}

void SwarmClustering::outputOtus(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const char sep, const std::string sepAbundance) {

    std::ofstream oStream(oFile);
    std::stringstream sStream;

    AmpliconCollection* ac = 0;
    Otu* otu = 0;

    for (auto i = 0; i < otus.size(); i++) {

        otu = otus[i];
        ac = pools.get(otu->poolId);

        if (!otu->attached) {

            sStream << (*ac)[otu->seedId].id << sepAbundance << (*ac)[otu->seedId].abundance;

            for (auto memberIter = otu->members.begin() + 1; memberIter != otu->members.end(); memberIter++) {
                sStream << sep << (*ac)[memberIter->id].id << sepAbundance << (*ac)[memberIter->id].abundance;
            }

            sStream << std::endl;
            oStream << sStream.rdbuf();
            sStream.str(std::string());

        }

    }

    oStream.close();

}

void SwarmClustering::outputStatistics(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const char sep) {

    std::ofstream oStream(oFile);
    std::stringstream sStream;

    AmpliconCollection* ac = 0;
    Otu* otu = 0;

    for (auto i = 0; i < otus.size(); i++) {

        otu = otus[i];
        ac = pools.get(otu->poolId);

        if (!otu->attached) {

            sStream << otu->numUniqueSequences << sep << otu->mass << sep << (*ac)[otu->seedId].id << sep << otu->seedAbundance << sep << otu->numSingletons << sep << otu->maxGen << sep << otu->maxRad << std::endl;
            oStream << sStream.rdbuf();
            sStream.str(std::string());

        }

    }

    oStream.close();

}

void SwarmClustering::outputSeeds(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const std::string sepAbundance) {

    std::ofstream oStream(oFile);
    std::stringstream sStream;

    AmpliconCollection* ac = 0;
    Otu* otu = 0;

    for (auto i = 0; i < otus.size(); i++) {

        otu = otus[i];
        ac = pools.get(otu->poolId);

        if (!otu->attached) {
            sStream << ">" << (*ac)[otu->seedId].id << sepAbundance << otu->mass << std::endl << (*ac)[otu->seedId].seq << std::endl;
            oStream << sStream.rdbuf();
            sStream.str(std::string());
        }

    }

    oStream.close();

}

void SwarmClustering::outputDereplicate(const AmpliconPools& pools, const std::vector<Otu*>& otus, const SwarmConfig& sc) {

    std::ofstream oStreamInternals, oStreamOtus, oStreamStatistics, oStreamSeeds;
    std::stringstream sStreamInternals, sStreamOtus, sStreamStatistics, sStreamSeeds;

    if (sc.outInternals) oStreamInternals.open(sc.oFileInternals);
    if (sc.outOtus) oStreamOtus.open(sc.oFileOtus);
    if (sc.outStatistics) oStreamStatistics.open(sc.oFileStatistics);
    if (sc.outSeeds) oStreamSeeds.open(sc.oFileSeeds);

    for (auto i = 0; i < otus.size(); i++) {

        Otu& otu = *(otus[i]);
        AmpliconCollection& ac = *(pools.get(otu.poolId));

        if (sc.outInternals) {

            for (auto memberIter = otu.members.begin() + 1; memberIter != otu.members.end(); memberIter++) {

                sStreamInternals << ac[otu.seedId].id << sc.sepInternals << ac[memberIter->id].id << sc.sepInternals << 0 << sc.sepInternals << (i + 1) << sc.sepInternals << 0 << std::endl;
                oStreamInternals << sStreamInternals.rdbuf();
                sStreamInternals.str(std::string());

            }

        }

        if (sc.outOtus) {

            sStreamOtus << ac[otu.seedId].id << sc.sepAbundance << ac[otu.seedId].abundance;

            for (auto memberIter = otu.members.begin() + 1; memberIter != otu.members.end(); memberIter++) {
                sStreamOtus << sc.sepOtus << ac[memberIter->id].id << sc.sepAbundance << ac[memberIter->id].abundance;
            }

            sStreamOtus << std::endl;
            oStreamOtus << sStreamOtus.rdbuf();
            sStreamOtus.str(std::string());

        }

        if (sc.outStatistics) {

            sStreamStatistics << otu.numUniqueSequences << sc.sepStatistics << otu.mass << sc.sepStatistics << ac[otu.seedId].id << sc.sepStatistics << otu.seedAbundance << sc.sepStatistics << otu.numSingletons << sc.sepStatistics << 0 << sc.sepStatistics << 0 << std::endl;
            oStreamStatistics << sStreamStatistics.rdbuf();
            sStreamStatistics.str(std::string());

        }

        if (sc.outSeeds) {

            sStreamSeeds << ">" << ac[otu.seedId].id << sc.sepAbundance << otu.mass << std::endl << ac[otu.seedId].seq << std::endl;
            oStreamSeeds << sStreamSeeds.rdbuf();
            sStreamSeeds.str(std::string());

        }

    }

    if (sc.outInternals) oStreamInternals.close();
    if (sc.outOtus) oStreamOtus.close();
    if (sc.outStatistics) oStreamStatistics.close();
    if (sc.outSeeds) oStreamSeeds.close();

}

}