#include "../include/RelationTest.hpp"

namespace GeFaST {


void RelationPrecursor::add(const StringIteratorPair seg, const numSeqs_t id) {

    auto iter = mapping.find(seg);
    numSeqs_t segId;

    if (iter != mapping.end()) {
        segId = iter->second;
    } else {

        segId = mapping.size();
        mapping[seg] = segId;

    }

    pairs.push_back(std::make_pair(segId, id));

}

void RelationPrecursor::clear() {

    std::vector<std::pair<numSeqs_t, numSeqs_t>>().swap(pairs);
    std::map<StringIteratorPair, numSeqs_t, lessStringIteratorPair>().swap(mapping);

}





struct CompareAbsoluteValues {

    bool operator() (const long long lhs, const long long rhs) {
        return llabs(lhs) < llabs(rhs);
    }

} cmpAbsVals;


RankedLabels::RankedLabels() {

    labels_ = 0;
    size_ = 0;
    capa_= 0;

}


RankedLabels::RankedLabels(const numSeqs_t capacity) {

    labels_ = new long long[capacity];
    size_ = 0;
    capa_= capacity;

}

RankedLabels& RankedLabels::operator=(const RankedLabels& other) {

    // check for self-assignment
    if (&other == this) {
        return *this;
    }

    delete[] labels_;

    size_ = other.size_;
    capa_ = other.capa_;

    labels_ = new long long[capa_];
    for (numSeqs_t i = 0; i < size_; i++) {
        labels_[i] = other.labels_[i];
    }

    return *this;

}

RankedLabels::~RankedLabels() {
    delete[] labels_;
}

numSeqs_t RankedLabels::add(const numSeqs_t lab) {

    labels_[size_] = lab + 1;

    return size_++;

}

numSeqs_t RankedLabels::unrank(const numSeqs_t r) {
    return llabs(labels_[r]) - 1;
}

bool RankedLabels::contains(const numSeqs_t lab) {

    auto iter = std::lower_bound(labels_, labels_ + size_, lab + 1, cmpAbsVals);

    return (iter != (labels_ + size_)) && (lab + 1 == *iter);

}

bool RankedLabels::containsRank(const numSeqs_t r) {
    return (r < size_) && (labels_[r] > 0);
}

void RankedLabels::remove(const numSeqs_t lab) {

    auto iter = std::lower_bound(labels_, labels_ + size_, lab + 1, cmpAbsVals);

    if (iter != (labels_ + size_)) {
        *iter = -llabs(*iter);
    }

}

void RankedLabels::swap(RankedLabels& other) {

    std::swap(labels_, other.labels_);
    std::swap(size_, other.size_);
    std::swap(capa_, other.capa_);

}

numSeqs_t RankedLabels::size() {
    return size_;
}

}