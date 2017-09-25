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

#include <algorithm>
#include <iostream>

#include "../include/Base.hpp"

namespace GeFaST {

#if QGRAM_FILTER
char acgtuMap[256] =
        {
                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, 4, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, 4, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
        };
#endif

// ===== amplicon comparer structures =====

bool AmpliconCompareAlph::operator()(const Amplicon& amplA, const Amplicon& amplB) {
    return strcmp(amplA.seq, amplB.seq) < 0;
}

bool AmpliconCompareLen::operator()(const Amplicon& amplA, const Amplicon& amplB) {
    return (amplA.len < amplB.len);
}

bool AmpliconCompareAbund::operator()(const Amplicon& amplA, const Amplicon& amplB) {
    return (amplA.abundance < amplB.abundance);
}

bool AmpliconSeqEqual::operator()(const Amplicon& amplA, const Amplicon& amplB) {
    return strcmp(amplA.seq, amplB.seq) == 0;
}

// ===== substrings & segments =====

Substrings selectSubstrs(const lenSeqs_t selfLen, const lenSeqs_t partnerLen, const lenSeqs_t segIndex, const lenSeqs_t t, const lenSeqs_t k) {

    lenSeqs_t p = 0;
    lenSeqs_t l = 0;

    lenSeqs_t i = segIndex + 1; // computations from PJ-paper use segment indices from 1 to t+1

    // by even-partitioning scheme
    // * first t + k - d segments: length floor(|s| / (t+l))
    // * last d segments: length floor(|s| / (t+k)) + 1
    lenSeqs_t d = partnerLen - (partnerLen / (t + k)) * (t + k);
    l = (partnerLen / (t + k)) + (i > t + k - d);
    p = 1 + (i - 1) * l - (i > t + k - d) * (t + k - d);


    // multimatch-aware substring selection boundaries
//    // calculations from paper (assume that positions p start at 1)
//    lenSeqs_t lowerL = std::max((lenSeqs_t)1, p - (i - 1));
//    lenSeqs_t upperL = std::min(s.length() - l + 1, p + (i - 1));
//
//    lenSeqs_t lowerR = std::max((lenSeqs_t)1, p + (s.length() - len) - (t + k - i));
//    lenSeqs_t upperR = std::min(s.length() - l + 1, p + (s.length() - len) + (t + k - i));
//
//    lenSeqs_t lower = std::max(lowerL, lowerR);
//    lenSeqs_t upper = std::min(upperL, upperR);


//    // problem: underflow of non-negative type lenSeqs_t
//    lenSeqs_t lower = std::max({(lenSeqs_t)1, p - (i - 1), p + (s.length() - len) - (t + k - i)});
    lenSeqs_t lower = 1;
    if (p > (i - 1)) {
        lower = p - (i - 1);
    }
    if (((p + (selfLen - partnerLen)) > (t + k - i)) && ((p + (selfLen - partnerLen) - (t + k - i)) > lower)) {
        lower = p + (selfLen - partnerLen) - (t + k - i);
    }

    // std::min okay
    lenSeqs_t upper = std::min({
                                       selfLen - l + 1,
                                       p + (i - 1),
                                       p + (selfLen - partnerLen) + (t + k - i)
                               });

    return Substrings(lower - 1, upper - 1, l); // -1 for translation to index starting at 0

}

Substrings selectSubstrsBackward(const lenSeqs_t selfLen, const lenSeqs_t partnerLen, const lenSeqs_t segIndex, const lenSeqs_t t, const lenSeqs_t k) {

    lenSeqs_t p = 0;
    lenSeqs_t l = 0;

    lenSeqs_t i = segIndex + 1; // computations from PJ-paper use segment indices from 1 to t+1

    // by even-partitioning scheme
    // * first t + k - d segments: length floor(|s| / (t+l))
    // * last d segments: length floor(|s| / (t+k)) + 1
    lenSeqs_t d = partnerLen - (partnerLen / (t + k)) * (t + k);
    l = (partnerLen / (t + k)) + (i > t + k - d);
    p = 1 + (i - 1) * l - (i > t + k - d) * (t + k - d);


    // multimatch-aware substring selection boundaries
//    // calculations from paper (assume that positions p start at 1)
//    lenSeqs_t lowerL = std::max((lenSeqs_t)1, p - (i - 1));
//    lenSeqs_t upperL = std::min(s.length() - l + 1, p + (i - 1));
//
//    lenSeqs_t lowerR = std::max((lenSeqs_t)1, p - (len - s.length()) - (t + k - i));
//    lenSeqs_t upperR = std::min(s.length() - l + 1, p + (len - s.length()) + (t + k - i));
//
//    lenSeqs_t lower = std::max(lowerL, lowerR);
//    lenSeqs_t upper = std::min(upperL, upperR);


//    // problem: underflow of non-negative type lenSeqs_t
//    lenSeqs_t lower = std::max({(lenSeqs_t)1, p - (i - 1), p + (s.length() - len) - (t + k - i)});

    lenSeqs_t lower = 1;
    if (p > (i - 1)) {
        lower = p - (i - 1);
    }
    if ((p > ((partnerLen - selfLen) + (t + k - i))) && ((p - (partnerLen - selfLen) - (t + k - i)) > lower)) {
        lower = p - (partnerLen - selfLen) - (t + k - i);
    }

    // std::min okay
    lenSeqs_t upper = std::min({
                                       selfLen - l + 1,
                                       p + (i - 1),
                                       p - (partnerLen - selfLen) + (t + k - i)
                               });

    return Substrings(lower - 1, upper - 1, l); // -1 for translation to index starting at 0

}

// The first t + k - d segments have length floor(seqLen / (t + k)), while the last d segments have length ceil(seqLen / (t + k)).
// Since d is calculated by seqLen - floor(seqLen / (t + k)) * (t + k), longer segments exist only if seqLen is not divisible by (t + k)
// and their lengths then higher by exactly one.
void selectSegments(Segments& segments, const lenSeqs_t seqLen, const lenSeqs_t t, const lenSeqs_t k) {

    lenSeqs_t d = seqLen - (seqLen / (t + k)) * (t + k);
    lenSeqs_t p = 0;
    lenSeqs_t i = 0;


    lenSeqs_t l = (seqLen / (t + k));
    for (; i < t + k - d; i++) {

        segments[i] = std::make_pair(p, l);
        p += l;

    }


    l++;
    for (; i < t + k; i++) {

        segments[i] = std::make_pair(p, l);
        p += l;

    }

}


// ===== AmpliconCollection =====

AmpliconCollection::AmpliconCollection(const numSeqs_t capacity, const std::vector<std::pair<lenSeqs_t, numSeqs_t>>& counts) {

    amplicons_ = new Amplicon[capacity];
    size_ = 0;
    capacity_ = capacity;
    numLengths_ = counts.size();
    counts_ = new std::pair<lenSeqs_t, numSeqs_t>[numLengths_];
    for (lenSeqs_t i = 0; i < numLengths_; i++) {
        counts_[i] = counts[i];
    }

}

AmpliconCollection::~AmpliconCollection() {

    delete[] counts_;
    delete[] amplicons_;

}

void AmpliconCollection::push_back(const Amplicon& ampl) {
    amplicons_[size_++] = ampl;
}

Amplicon& AmpliconCollection::operator[](const numSeqs_t i) {
    return amplicons_[i];
}

const Amplicon& AmpliconCollection::operator[](const numSeqs_t i) const {
    return amplicons_[i];
}

Amplicon& AmpliconCollection::front() const {
    return amplicons_[0];
}

Amplicon& AmpliconCollection::back() const {
    return amplicons_[size_ - 1];
}

Amplicon* AmpliconCollection::begin() const {
    return amplicons_;
}

Amplicon* AmpliconCollection::end() const {
    return amplicons_ + size_;
}

numSeqs_t AmpliconCollection::size() const {
    return size_;
}

numSeqs_t AmpliconCollection::numSeqsOfLen(const lenSeqs_t len) const {

    auto iter = std::lower_bound(counts_, counts_ + numLengths_, len,
                                 [](const std::pair<lenSeqs_t, numSeqs_t>& pair, lenSeqs_t l) {
                                     return pair.first < l;
                                 });

    return (iter != (counts_ + numLengths_) && iter->first == len) ? iter->second : 0;

}

void AmpliconCollection::reserve(const numSeqs_t newCapacity) {

    if (capacity_ >= newCapacity || size_ >= newCapacity) return;

    Amplicon* tmp = new Amplicon[newCapacity];
    for (numSeqs_t i = 0; i < size_; i++) {
        tmp[i] = amplicons_[i];
    }

    delete[] amplicons_;
    amplicons_ = tmp;

}


// ===== AmpliconPools =====

AmpliconPools::AmpliconPools(std::map<lenSeqs_t, numSeqs_t>& counts, const unsigned long long capacity, const lenSeqs_t threshold) {

    strings_ = new char[capacity];
    nextPos_ = strings_;
    capacity_ = capacity;

    if (counts.size() != 0) {

        numSeqs_t curPoolSize = 0;
        std::vector<std::pair<lenSeqs_t, numSeqs_t>> localCounts;

        // start the pool that will comprise the shortest amplicons
        lenSeqs_t lastLen = counts.begin()->first;
        curPoolSize += counts.begin()->second;
        localCounts.push_back(std::make_pair(lastLen, curPoolSize));
        counts.begin()->second = pools_.size();

        // iterate over the remaining length groups, starting a new one if a break is detected
        for (auto iter = ++counts.begin(); iter != counts.end(); iter++) {

            if ((lastLen + threshold) < iter->first) { // new pool

                pools_.push_back(new AmpliconCollection(curPoolSize, localCounts));

                curPoolSize = 0;
                localCounts.clear();

            }

            localCounts.push_back(std::make_pair(iter->first, iter->second));
            curPoolSize += iter->second;
            iter->second = pools_.size();
            lastLen = iter->first;

        }

        pools_.push_back(new AmpliconCollection(curPoolSize, localCounts));

    }

}

AmpliconPools::~AmpliconPools() {

    for (auto iter = pools_.begin(); iter != pools_.end(); iter++) {
        delete *iter;
    }

    delete[] strings_;

}

void AmpliconPools::add(const lenSeqs_t i, const std::string& header, const std::string& sequence, const numSeqs_t abundance) {

    char* h = nextPos_;
    strcpy(nextPos_, header.c_str());
    nextPos_ += header.length() + 1;

    char* s = nextPos_;
    nextPos_ = strcpy(nextPos_, sequence.c_str());
    nextPos_ += sequence.length() + 1;

    pools_[i]->push_back(Amplicon(h, s, sequence.length(), abundance));

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


//TODO? custom hash function (e.g. FNV) to avoid construction of temporary string
size_t hashStringIteratorPair::operator()(const StringIteratorPair& p) const {
    return std::hash<std::string>{}(std::string(p.first, p.second));
}

bool equalStringIteratorPair::operator()(const StringIteratorPair& lhs, const StringIteratorPair& rhs) const {
    return ((lhs.second - lhs.first) == (rhs.second - rhs.first)) && std::equal(lhs.first, lhs.second, rhs.first);
}

bool lessStringIteratorPair::operator()(const StringIteratorPair& a, const StringIteratorPair& b) const {
    return std::lexicographical_compare(a.first, a.second, b.first, b.second);
}


std::vector<Subpool> getSubpoolBoundaries(const AmpliconCollection& ac, const numSeqs_t num, const lenSeqs_t threshold) { // with even-partitioning scheme

    std::vector<Subpool> subpools;

    // use only one subpool if there are not enough sequences or every index-only part would start at the same position
    if ((ac.size() < num) || ((ac.back().len - ac.front().len) <= threshold)) {

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

        lenSeqs_t lowerBound = (ac[subpools[i].beginMatch].len > threshold) * (ac[subpools[i].beginMatch].len - threshold);

        numSeqs_t j = subpools[i].beginMatch;
        for (; j > 0 && ac[j - 1].len >= lowerBound; j--) { }
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
    if (ac.size() < num || ((ac.back().len - ac.front().len) <= threshold)) {

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

        lenSeqs_t upperBound = ac[subpools[i].beginIndex - 1].len + threshold;

        numSeqs_t j = subpools[i].beginIndex;
        for (; j < ac.size() && ac[j].len <= upperBound; j++) { }
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

}