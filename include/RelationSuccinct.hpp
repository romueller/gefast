#ifndef GEFAST_RELATIONSUCCINCT_HPP
#define GEFAST_RELATIONSUCCINCT_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <unordered_map>
#include <StaticHybridTree.hpp>
#include <StaticBasicRectangularTree.hpp>
#include <StaticUnevenRectangularTree.hpp>

#include "Base.hpp"
#include "StaticBasicTree.hpp"
#include "StaticRowTree.hpp"
#include "StaticHybridRowTree.hpp"
#include "StaticMiniRowTree.hpp"
#include "StaticMiniK2Tree.hpp"
#include "StaticUnevenRectangularOrMiniTree.hpp"

namespace GeFaST {

struct RelationPrecursor {//TODO mapping: map or unordered_map?

    std::vector<std::pair<numSeqs_t, numSeqs_t>> pairs;
    std::map<StringIteratorPair, numSeqs_t, lessStringIteratorPair> mapping;

    void add(const StringIteratorPair seg, const numSeqs_t id);

    void clear();

};

struct SuccinctConfig {

    unsigned long succK, succKR, succKC, succUK, succUH, succLK, succMB;

    SuccinctConfig(unsigned long k, unsigned long kr, unsigned long kc, unsigned long uk, unsigned long uh, unsigned long lk, unsigned long mb) {
        succK = k;
        succKR = kr;
        succKC = kc;
        succUK = uk;
        succUH = uh;
        succLK = lk;
        succMB = mb;
    }

};

// maps arbitary ascending sequence of n unique (positive) integers onto [0:n-1]
class RankedAscendingLabels {

public:
    RankedAscendingLabels();

    RankedAscendingLabels(const numSeqs_t capacity);

    RankedAscendingLabels(std::vector<std::pair<numSeqs_t, numSeqs_t>>& pairs);

    RankedAscendingLabels& operator=(const RankedAscendingLabels& other);

    ~RankedAscendingLabels();

    numSeqs_t add(const numSeqs_t lab);

    numSeqs_t unrank(const numSeqs_t r);

    bool contains(const numSeqs_t lab);

    bool containsRank(const numSeqs_t r);

    void remove(const numSeqs_t lab);

    void swap(RankedAscendingLabels& other);

    numSeqs_t size();

private:
    long long* labels_;
    numSeqs_t size_;
    numSeqs_t capa_;

};


template<typename S, typename T>
class SharingRollingIndices {

public:
    struct Row {

        S shared;
        std::vector<T> indices;

        Row() {
            // nothing to do
        }

        Row(lenSeqs_t w, numSeqs_t sharedCapacity) : shared(S(sharedCapacity)) {
            indices = std::vector<T>(w);
        }

        Row& operator=(const Row& other) {

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            shared = other.shared;
            indices = other.indices;

            return *this;

        }

    };

    SharingRollingIndices(lenSeqs_t t, lenSeqs_t w, bool f, bool s = true) {

        threshold_ = t;
        width_ = w;
        forward_ = f;
        shrink_ = s;

        empty_ = T();
        emptyRow_ = Row(0, 0);

    }


    // return the indices for the specified length
    Row& getIndicesRow(const lenSeqs_t len) {

        auto iter = rows_.find(len);

        return (iter != rows_.end()) ? iter->second : emptyRow_;

    }


    // return the index corresponding to the specified length and segment
    T& getIndex(const lenSeqs_t len, const lenSeqs_t i) {

        if (i >= width_) return empty_;

        auto iter = rows_.find(len);

        return (iter != rows_.end()) ? (iter->second).indices[i] : empty_;

    }


    // add new row (and remove then outdated rows)
    void roll(const lenSeqs_t len, const numSeqs_t sharedCapacity) {

        if (rows_.find(len) == rows_.end()) {

            rows_[len] = Row(width_, sharedCapacity);
            if (shrink_) shrink(len);

        }


    }


    // remove outdated rows
    void shrink(const lenSeqs_t cur) {

        if (forward_) {

            auto bound = rows_.lower_bound(cur - threshold_);
            rows_.erase(rows_.begin(), bound);

        } else {

            auto bound = rows_.upper_bound(cur + threshold_);
            rows_.erase(bound, rows_.end());

        }

    }

    bool contains(const lenSeqs_t len) {
        return (rows_.find(len) != rows_.end());
    }

    lenSeqs_t minLength() const {
        return rows_.begin()->first;
    }
    lenSeqs_t maxLength() const {
        return rows_.rbegin()->first;
    }


private:

    lenSeqs_t threshold_; // limits number of rows when applying shrink()
    lenSeqs_t width_; // number of columns / segments per row

    std::map<lenSeqs_t, Row> rows_; // indices grid

    T empty_; // empty (dummy) index returned for out-of-bounds queries
    Row emptyRow_; // empty (dummy) row returned for out-of-bounds queries

    bool forward_; // flag indicating whether rolling forwards (increasingly larger lengths are 'inserted') or backwards (shorter lengths are 'inserted')
    bool shrink_; // flag indicating whether roll() automatically shrinks the index

};


template<typename O, typename S>
class AddingTree : public KrKcTree<bool> {

public:
    using KrKcTree<bool>::KrKcTree; // adds ctors

    void addSuccessorCountsOf(const O& p, std::unordered_map<numSeqs_t, lenSeqs_t>& candCnts, S& labels) {

        if (L_.empty()) return;

        std::queue<SubrowInfo> queue, nextLevelQueue;
        size_type lenT = T_.size();

        if (lenT == 0) {

            size_type offset = p * numCols_;
            for (size_type i = 0; i < numCols_; i++) {
                if (L_[offset + i] && labels.containsRank(i)) {
                    candCnts[labels.unrank(i)]++;
                }
            }

        } else {

            // successorsPosInit
            size_type nr = numRows_/ kr_;
            size_type nc = numCols_/ kc_;
            size_type relP = p;
            for (size_type j = 0, dq = 0, z = kc_ * (relP / nr); j < kc_; j++, dq += nc, z++) {
                queue.push(SubrowInfo(dq, z));
            }

            // successorsPos
            relP %= nr;
            nr /= kr_;
            nc /= kc_;
            for (; nr > 1; relP %= nr, nr /= kr_, nc /= kc_) {

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (T_[cur.z]) {

                        auto y = R_.rank(cur.z + 1) * kr_ * kc_ + kc_ * (relP / nr);

                        for (size_type j = 0, newDq = cur.dq; j < kc_; j++, newDq += nc, y++) {
                            nextLevelQueue.push(SubrowInfo(newDq, y));
                        }

                    }

                    queue.pop();

                }

                queue.swap(nextLevelQueue);

            }


            while (!queue.empty()) {

                auto& cur = queue.front();

                if (T_[cur.z]) {

                    auto y = R_.rank(cur.z + 1) * kr_ * kc_ + kc_ * (relP / nr) - lenT;

                    for (size_type j = 0, newDq = cur.dq; j < kc_; j++, newDq += nc, y++) {
                        if (L_[y] && labels.containsRank(newDq)) {
                            candCnts[labels.unrank(newDq)]++;
                        }
                    }

                }

                queue.pop();

            }

        }

    }

};

template<typename S>
class K2TreeBinaryRelation {

    typedef StringIteratorPair O;
    typedef numSeqs_t L;

public:
    K2TreeBinaryRelation() {
        // nothing to do
    }

    K2TreeBinaryRelation(RelationPrecursor& ir, S& labels, const SuccinctConfig& sc) {

//        binRel_ = BasicK2Tree<bool>(ir.pairs, sc.succK);
//        binRel_ = HybridK2Tree<bool>(ir.pairs, sc.succUK, sc.succUH, sc.succLK);
        binRel_ = KrKcTree<bool>(ir.pairs, sc.succKR, sc.succKC);
//        binRel_ = UnevenKrKcTree<bool>(ir.pairs, sc.succKR, sc.succKC);
//        binRel_ = AddingTree<numSeqs_t, S>(ir.pairs, sc.succKR, sc.succKC);
//        binRel_ = UnevenKrKcOrMiniTree<bool>(ir.pairs, sc.succKR, sc.succKC, sc.succMB);
        segIdMap_.swap(ir.mapping);
        labels_ = &labels;

    }

    K2TreeBinaryRelation& operator=(const K2TreeBinaryRelation& other) {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        binRel_ = other.binRel_;
        segIdMap_ = other.segIdMap_;
        labels_ = other.labels_;

        return *this;

    }

    ~K2TreeBinaryRelation() {
        // nothing to do
    }

    bool areRelated(const O& obj, const L& lab) {
        return containsLabel(lab) && containsObject(obj) && binRel_.areRelated(segIdMap_[obj], lab);
    }

    bool containsLabel(const L& lab) {
        return (labels_ != 0) && labels_->contains(lab);
    }

    bool containsObject(const O& obj) {
        return segIdMap_.find(obj) != segIdMap_.end();
    }

    std::vector<L> getLabelsOf(const O& obj) {

        std::vector<L> res;
        if (containsObject(obj) && labels_ != 0) {

            auto tmp = binRel_.getSuccessors(segIdMap_[obj]);
            res.reserve(tmp.size());

            for (auto& l : tmp) {
                if (labels_->containsRank(l)) {
                    res.push_back(labels_->unrank(l));
                }
            }

        }

        return res;

    }

    void addLabelCountsOf(const O& obj, std::unordered_map<numSeqs_t, lenSeqs_t>& candCnts) {

        if (containsObject(obj) && labels_ != 0) {

            for (auto& l : binRel_.getSuccessors(segIdMap_[obj])) {
                if (labels_->containsRank(l)) {
                    candCnts[labels_->unrank(l)]++;
                }
            }
//            binRel_.addSuccessorCountsOf(segIdMap_[obj], candCnts, *labels_);

        }

    }

    unsigned long countPairs() {
        return binRel_.countLinks();
    }

    unsigned long getNumRows() {
        return binRel_.getNumRows();
    }

    unsigned long getNumCols() {
        return binRel_.getNumCols();
    }

    void removeLabel(const L& lab) {
        if (labels_ != 0) {
            labels_->remove(lab);
        }
    }


private:
//    BasicK2Tree<bool> binRel_;
//    HybridK2Tree<bool> binRel_;
    KrKcTree<bool> binRel_;
//    UnevenKrKcTree<bool> binRel_;
//    AddingTree<numSeqs_t, S> binRel_;
//    UnevenKrKcOrMiniTree<bool> binRel_;
    std::map<StringIteratorPair, numSeqs_t, lessStringIteratorPair> segIdMap_;
    S* labels_;

};

template<typename S>
class RowTreeBinaryRelation {

    typedef StringIteratorPair O;
    typedef numSeqs_t L;

public:
    RowTreeBinaryRelation() {
        // nothing to do
    }

    RowTreeBinaryRelation(RelationPrecursor& ir, S& labels, const SuccinctConfig& sc) {

        // TODO optimise ( = no full copy), e.g. sorting + equal_range  + sub-copy OR ctor from iterators
#if 0
        std::vector<RelationList> rows(ir.mapping.size());
        for (auto& p : ir.pairs) {
            rows[p.first].push_back(p.second);
        }
        std::vector<std::pair<numSeqs_t, numSeqs_t>>().swap(ir.pairs);

        for (auto i = 0; i < rows.size(); i++) {

            binRel_.emplace_back(rows[i], sc.succK);
//            binRel_.emplace_back(rows[i], sc.succUK, sc.succUH, sc.succLK);
            RelationList().swap(rows[i]);

        }
#else
        std::sort(ir.pairs.begin(), ir.pairs.end());

        auto begin = ir.pairs.begin();
        for (numSeqs_t i = 0; i < ir.mapping.size(); i++) {

            auto end = std::upper_bound(begin, ir.pairs.end(), i,
                                        [](numSeqs_t val, const std::pair<numSeqs_t, numSeqs_t>& pair) {
                                            return val < pair.first;
                                        });
            binRel_.emplace_back(begin, end, sc.succK);
//            binRel_.emplace_back(begin, end, sc.succUK, sc.succUH, sc.succLK);
            begin = end;

        }
#endif

        segIdMap_.swap(ir.mapping);
        labels_ = &labels;

    }

    RowTreeBinaryRelation& operator=(const RowTreeBinaryRelation& other) {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        binRel_ = other.binRel_;
        segIdMap_ = other.segIdMap_;
        labels_ = other.labels_;

        return *this;

    }

    ~RowTreeBinaryRelation() {
        // nothing to do
    }

    bool areRelated(const O& obj, const L& lab) {
        return containsLabel(lab) && containsObject(obj) && binRel_[segIdMap_[obj]].isNotNull(lab);
    }

    bool containsLabel(const L& lab) {
        return (labels_ != 0) && labels_->contains(lab);
    }

    bool containsObject(const O& obj) {
        return segIdMap_.find(obj) != segIdMap_.end();
    }

    std::vector<L> getLabelsOf(const O& obj) {

        std::vector<L> res;
        if (containsObject(obj) && labels_ != 0) {

            auto tmp = binRel_[segIdMap_[obj]].getAllPositions();
            res.reserve(tmp.size());

            for (auto& l : tmp) {
                if (labels_->containsRank(l)) {
                    res.push_back(labels_->unrank(l));
                }
            }

        }

        return res;

    }

    void addLabelCountsOf(const O& obj, std::unordered_map<numSeqs_t, lenSeqs_t>& candCnts) {

        if (containsObject(obj) && labels_ != 0) {

            for (auto& l : binRel_[segIdMap_[obj]].getAllPositions()) {
                if (labels_->containsRank(l)) {
                    candCnts[labels_->unrank(l)]++;
                }
            }

        }

    }

    unsigned long countPairs() {

        unsigned long sum = 0;
        for (unsigned long i = 0; i < binRel_.size(); i++) {
            sum += binRel_[i].countElements();
        }

        return sum;

    }

    unsigned long getNumRows() {
        return binRel_.size();
    }

    unsigned long getNumCols() {

        unsigned long max = 0;
        for (unsigned long i = 0; i < binRel_.size(); i++) {
            max = std::max(max, binRel_[i].getLength());
        }

        return max;

    }

    void removeLabel(const L& lab) {
        if (labels_ != 0) {
            labels_->remove(lab);
        }
    }


private:
    std::vector<BasicRowTree<bool>> binRel_;
//    std::vector<HybridRowTree<bool>> binRel_;
    std::map<StringIteratorPair, numSeqs_t, lessStringIteratorPair> segIdMap_;
    S* labels_;

};

template<typename S>
class MiniRowTreeBinaryRelation {

    typedef StringIteratorPair O;
    typedef numSeqs_t L;

public:
    MiniRowTreeBinaryRelation() {
        // nothing to do
    }

    MiniRowTreeBinaryRelation(RelationPrecursor& ir, S& labels, const SuccinctConfig& sc) {

        // TODO optimise ( = no full copy), e.g. sorting + equal_range  + sub-copy OR ctor from iterators
#if 0
        std::vector<RelationList> rows(ir.mapping.size());
        for (auto& p : ir.pairs) {
            rows[p.first].push_back(p.second);
        }
        std::vector<std::pair<numSeqs_t, numSeqs_t>>().swap(ir.pairs);

        for (auto i = 0; i < rows.size(); i++) {

            if (rows[i].size() > 5) {
                binRel_.push_back(new BasicRowTree<bool>(rows[i], 10));
//                binRel_.push_back(new HybridRowTree<bool>(rows[i], 7, 2, 3));
            } else {
                binRel_.push_back(new MiniRowTree<bool>(rows[i]));
            }
            RelationList().swap(rows[i]);

        }
#else
        std::sort(ir.pairs.begin(), ir.pairs.end());

        auto begin = ir.pairs.begin();
        for (numSeqs_t i = 0; i < ir.mapping.size(); i++) {

            auto end = std::upper_bound(begin, ir.pairs.end(), i,
                                        [](numSeqs_t val, const std::pair<numSeqs_t, numSeqs_t>& pair) {
                                            return val < pair.first;
                                        });

            if (end - begin > sc.succMB) {
//                binRel_.push_back(new BasicRowTree<bool>(begin, end, sc.succK));
                binRel_.push_back(new HybridRowTree<bool>(begin, end, sc.succUK, sc.succUH, sc.succLK));
            } else {
                binRel_.push_back(new MiniRowTree<bool>(begin, end));
            }

            begin = end;

        }
#endif
        segIdMap_.swap(ir.mapping);
        labels_ = &labels;

    }

    MiniRowTreeBinaryRelation& operator=(const MiniRowTreeBinaryRelation& other) {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        for (auto i = 0; i < binRel_.size(); i++) {
            delete binRel_[i];
        }

        binRel_.clear();
        for (auto i = 0; i < other.binRel_.size(); i++) {
            binRel_.push_back(other.binRel_[i]->clone());
        }

        segIdMap_ = other.segIdMap_;
        labels_ = other.labels_;

        return *this;

    }

    ~MiniRowTreeBinaryRelation() {

        for (auto i = 0; i < binRel_.size(); i++) {
            delete binRel_[i];
        }

    }

    bool areRelated(const O& obj, const L& lab) {
        return containsLabel(lab) && containsObject(obj) && binRel_[segIdMap_[obj]]->isNotNull(lab);
    }

    bool containsLabel(const L& lab) {
        return (labels_ != 0) && labels_->contains(lab);
    }

    bool containsObject(const O& obj) {
        return segIdMap_.find(obj) != segIdMap_.end();
    }

    std::vector<L> getLabelsOf(const O& obj) {

        std::vector<L> res;
        if (containsObject(obj) && labels_ != 0) {

            auto tmp = binRel_[segIdMap_[obj]]->getAllPositions();
            res.reserve(tmp.size());

            for (auto &l : tmp) {
                if (labels_->containsRank(l)) {
                    res.push_back(labels_->unrank(l));
                }
            }

        }

        return res;

    }

    void addLabelCountsOf(const O& obj, std::unordered_map<numSeqs_t, lenSeqs_t>& candCnts) {

        if (containsObject(obj) && labels_ != 0) {

            for (auto& l : binRel_[segIdMap_[obj]]->getAllPositions()) {
                if (labels_->containsRank(l)) {
                    candCnts[labels_->unrank(l)]++;
                }
            }

        }

    }

    unsigned long countPairs() {

        unsigned long sum = 0;
        for (unsigned long i = 0; i < binRel_.size(); i++) {
            sum += binRel_[i]->countElements();
        }

        return sum;

    }

    unsigned long getNumRows() {
        return binRel_.size();
    }

    unsigned long getNumCols() {

        unsigned long max = 0;
        for (unsigned long i = 0; i < binRel_.size(); i++) {
            max = std::max(max, binRel_[i]->getLength());
        }

        return max;

    }

    void removeLabel(const L& lab) {
        if (labels_ != 0) {
            labels_->remove(lab);
        }
    }


private:
    std::vector<RowTree<bool>*> binRel_;
    std::map<StringIteratorPair, numSeqs_t, lessStringIteratorPair> segIdMap_;
    S* labels_;

};

template<typename S>
class MiniK2TreeBinaryRelation {

    typedef StringIteratorPair O;
    typedef numSeqs_t L;

public:
    MiniK2TreeBinaryRelation() {

        binRel_ = 0;
        labels_ = 0;

    }

    MiniK2TreeBinaryRelation(RelationPrecursor& ir, S& labels, const SuccinctConfig& sc) {


        if (ir.pairs.size() > sc.succMB) {
//            binRel_ = new BasicK2Tree<bool>(ir.pairs, sc.succK);
//            binRel_ = new HybridK2Tree<bool>(ir.pairs, sc.succUK, sc.succUH, sc.succLK);
//            binRel_ = new KrKcTree<bool>(ir.pairs, sc.succKR, sc.succKC);
            binRel_ = new UnevenKrKcTree<bool>(ir.pairs, sc.succKR, sc.succKC);
        } else{
            binRel_ = new MiniK2Tree<bool>(ir.pairs);
        }
        segIdMap_.swap(ir.mapping);
        labels_ = &labels;

    }

    MiniK2TreeBinaryRelation& operator=(const MiniK2TreeBinaryRelation& other) {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        delete binRel_;

        binRel_ = (other.binRel_ != 0) ? other.binRel_->clone() : 0;
        segIdMap_ = other.segIdMap_;
        labels_ = other.labels_;

        return *this;

    }

    ~MiniK2TreeBinaryRelation() {
        delete binRel_;
    }

    bool areRelated(const O& obj, const L& lab) {
        return containsLabel(lab) && containsObject(obj) && binRel_->areRelated(segIdMap_[obj], lab);
    }

    bool containsLabel(const L& lab) {
        return (labels_ != 0) && labels_->contains(lab);
    }

    bool containsObject(const O& obj) {
        return segIdMap_.find(obj) != segIdMap_.end();
    }

    std::vector<L> getLabelsOf(const O& obj) {

        std::vector<L> res;
        if (containsObject(obj) && labels_ != 0) {

            auto tmp = binRel_->getSuccessors(segIdMap_[obj]);
            res.reserve(tmp.size());

            for (auto& l : tmp) {
                if (labels_->containsRank(l)) {
                    res.push_back(labels_->unrank(l));
                }
            }

        }

        return res;

    }

    void addLabelCountsOf(const O& obj, std::unordered_map<numSeqs_t, lenSeqs_t>& candCnts) {

        if (containsObject(obj) && labels_ != 0) {

            for (auto& l : binRel_->getSuccessors(segIdMap_[obj])) {
                if (labels_->containsRank(l)) {
                    candCnts[labels_->unrank(l)]++;
                }
            }

        }

    }

    unsigned long countPairs() {
        return binRel_->countLinks();
    }

    unsigned long getNumRows() {
        return binRel_->getNumRows();
    }

    unsigned long getNumCols() {
        return binRel_->getNumCols();
    }

    void removeLabel(const L& lab) {
        if (labels_ != 0) {
            labels_->remove(lab);
        }
    }


private:
    K2Tree<bool>* binRel_;
    std::map<StringIteratorPair, numSeqs_t, lessStringIteratorPair> segIdMap_;
    S* labels_;

};

}

#endif //GEFAST_RELATIONSUCCINCT_HPP
