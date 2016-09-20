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

#include <algorithm>
#include <iostream>

#include "../include/Base.hpp"

namespace SCT_PJ {

// ===== amplicon comparer structures =====

struct AmpliconCompareAlph { // lexicographic, ascending
    bool operator()(Amplicon amplA, Amplicon amplB) {

        return (amplA.seq < amplB.seq);
    }
} ampliconCompareAlph;

struct AmpliconCompareLen { // by length, ascending
    bool operator()(Amplicon amplA, Amplicon amplB) {

        return (amplA.seq.length() < amplB.seq.length());
    }
} ampliconCompLen;

struct AmpliconCompareAbund { // by abundance, ascending
    bool operator()(Amplicon amplA, Amplicon amplB) {

        return (amplA.abundance < amplB.abundance);
    }
} ampliconCompAbund;

struct AmpliconSeqEqual {
    bool operator()(const Amplicon& amplA, const Amplicon& amplB) {
        return (amplA.seq == amplB.seq);
    }
} ampliconSeqEqual;


// ===== AmpliconPools =====

AmpliconPools::AmpliconPools() {
    // nothing to do
}

AmpliconPools::~AmpliconPools() {

    for (auto iter = pools_.begin(); iter != pools_.end(); iter++) {
        delete *iter;
    }

}

void AmpliconPools::add(AmpliconCollection* ac) {
    pools_.push_back(ac);
}


AmpliconCollection* AmpliconPools::get(const lenSeqs_t i) const {
    return (i < pools_.size()) ? pools_[i] : 0;
}

lenSeqs_t AmpliconPools::numPools() const {
    return pools_.size();
}

numSeqs_t AmpliconPools::numAmplicons() const {

    numSeqs_t sum = 0;
    for (auto iter = pools_.begin(); iter != pools_.end(); iter++) {
        sum += (*iter)->size();
    }

    return sum;

}


std::vector<Subpool> getSubpoolBoundaries(const AmpliconCollection& ac, const numSeqs_t num, const lenSeqs_t threshold) { // with even-partitioning scheme

    std::vector<Subpool> subpools;

    // use only one subpool if there are not enough sequences or every index-only part would start at the same position
    if ((ac.size() < num) || ((ac.back().seq.size() - ac.front().seq.size()) <= threshold)) {

        Subpool sp(0, 0, ac.size());
        subpools.push_back(sp);
        return subpools;

    }

    numSeqs_t l = ac.size() - (ac.size() / num) * num;

    // beginMatch and end for the first num - l subpools
    subpools.push_back(Subpool());
    subpools.back().end = ac.size() / num;
    for (numSeqs_t i = 2; i <= num - l; i++) {
        subpools.push_back(Subpool(0, subpools.back().end, subpools.back().end + ac.size() / num));
    }

    // beginMatch and end for the last l subpools
    for (numSeqs_t i = num - l + 1; i <= num; i++) {
        subpools.push_back(Subpool(0, subpools.back().end, subpools.back().end + ac.size() / num + 1));
    }

    // beginIndex for all subpools except the first one (beginIndex always 0 as initialised)
    // by moving backwards from beginMatch until sequences in the pool get too short w.r.t. the given threshold
    for (numSeqs_t i = 1; i < num; i++) {

        lenSeqs_t lowerBound = (ac[subpools[i].beginMatch].seq.length() > threshold) * (ac[subpools[i].beginMatch].seq.length() - threshold);

        numSeqs_t j = subpools[i].beginMatch;
        for (; j > 0 && ac[j - 1].seq.length() >= lowerBound; j--) { }
        subpools[i].beginIndex = j;

    }

    // combine subpools if their index-only part starts at the same position
    std::vector<Subpool> reducedSubpools;

    if (ac.size() > 1) {

        for (auto copyIter = subpools.begin(); copyIter != subpools.end();) {

            auto compIter = copyIter + 1;
            for (; (compIter != subpools.end()) && (compIter->beginIndex == copyIter->beginIndex); compIter++) {
                    copyIter->end = compIter->end;
            }

            reducedSubpools.push_back(*copyIter);
            copyIter = compIter;

        }

    }

    return reducedSubpools;

}


std::vector<Subpool> getSubpoolBoundariesBackward(const AmpliconCollection& ac, const numSeqs_t num, const lenSeqs_t threshold) { // with even-partitioning scheme

    std::vector<Subpool> subpools;

    // use only one subpool if there are not enough sequences or every index-only part would end at the same position
    if (ac.size() < num || ((ac.back().seq.size() - ac.front().seq.size()) <= threshold)) {

        Subpool sp(ac.size(), 0, ac.size());
        subpools.push_back(sp);
        return subpools;

    }


    numSeqs_t l = ac.size() - (ac.size() / num) * num;

    // beginMatch and beginIndex for the first num - l subpools
    subpools.push_back(Subpool(ac.size() / num, 0, ac.size()));
    for (numSeqs_t i = 2; i <= num - l; i++) {
        subpools.push_back(Subpool(subpools.back().beginIndex + ac.size() / num, subpools.back().beginIndex, ac.size()));
    }

    // beginMatch and beginIndex for the last l subpools
    for (numSeqs_t i = num - l + 1; i <= num; i++) {
        subpools.push_back(Subpool(subpools.back().beginIndex + ac.size() / num + 1, subpools.back().beginIndex, ac.size()));
    }

    // end for all subpools except the last one (end always ac.size() as initialised)
    // by moving backwards from beginIndex until sequences in the pool get too long w.r.t. the given threshold
    for (numSeqs_t i = 0; i < num - 1; i++) {

        lenSeqs_t upperBound = ac[subpools[i].beginIndex - 1].seq.length() + threshold;

        numSeqs_t j = subpools[i].beginIndex;
        for (; j < ac.size() && ac[j].seq.length() <= upperBound; j++) { }
        subpools[i].end = j;

    }

    // combine subpools if their index-only part ends at the same position
    std::vector<Subpool> reducedSubpools;

    if (ac.size() > 1) {

        for (auto copyIter = subpools.begin(); copyIter != subpools.end();) {

            auto compIter = copyIter + 1;
            for (; (compIter != subpools.end()) && (compIter->end == copyIter->end); compIter++) {
                    copyIter->beginIndex = compIter->beginIndex;
            }

            reducedSubpools.push_back(*copyIter);
            copyIter = compIter;

        }

    }

    return reducedSubpools;

}


// ===== LengthGroups =====

LengthGroups::LengthGroups() {
    // nothing to do
}

LengthGroups::LengthGroups(const AmpliconCollection& amplicons) {

    for (auto iter = amplicons.begin(); iter != amplicons.end(); iter++) {
        add(*iter);
    }

}

void LengthGroups::add(const Amplicon& ampl) {
    groups_[ampl.seq.length()].push_back(ampl);
}

std::map<lenSeqs_t, AmpliconCollection>& LengthGroups::getGroups() {
    return groups_;
}

numSeqs_t LengthGroups::size() {

    numSeqs_t res = 0;
    for (auto iter = groups_.begin(); iter != groups_.end(); iter++) {
        res += iter->second.size();
    }

    return res;

}

AmpliconPools* LengthGroups::pool(const lenSeqs_t threshold) {

    AmpliconPools* pools = new AmpliconPools();
    AmpliconCollection* curPool = 0;

    if (groups_.size() == 0) return pools;

    // sort alphabetically within the groups
    sortInGroupsByLength();


    // start the pool by initialising it with the length group comprising the shortest amplicons
    lenSeqs_t lastLen = groups_.begin()->first;
    curPool = new AmpliconCollection();
    curPool->insert(curPool->end(), groups_.begin()->second.begin(), groups_.begin()->second.end());
    groups_.begin()->second.erase(groups_.begin()->second.begin(), groups_.begin()->second.end());


    // iterate over the remaining length groups, appending to the current pool or starting a new one
    // if a break is detected
    for (auto iter = ++groups_.begin(); iter != groups_.end(); iter++) {

        if ((lastLen + threshold) < iter->first) {//new pool

            pools->add(curPool);
            curPool = new AmpliconCollection();

        }

        curPool->insert(curPool->end(), iter->second.begin(), iter->second.end());
        iter->second.erase(iter->second.begin(), iter->second.end());

        lastLen = iter->first;

    }

    pools->add(curPool);

    // prepare object for inserting new amplicons
    groups_ = std::map<lenSeqs_t, AmpliconCollection>();

    return pools;

}

void LengthGroups::sortInGroupsByLength() {

    for (auto iter = groups_.begin(); iter != groups_.end(); iter++) {
        sort(iter->second.begin(), iter->second.end(), ampliconCompareAlph);
    }

}

}