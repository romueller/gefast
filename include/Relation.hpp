/*
 * SCT-PJ
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

#ifndef SCT_PJ_RELATION_HPP
#define SCT_PJ_RELATION_HPP

#include <algorithm>
#include <mutex>
#include <unordered_map>
#include <vector>
#include <unordered_set>

#include "Base.hpp"


namespace SCT_PJ {

/*
 * Implementation of a dynamic binary relation between objects of type O and
 * labels of type L.
 * STL-based (unordered_map, vector).
 * Prefers object-based operations.
 */
template<typename O, typename L>
class SimpleBinaryRelation {

public:
    SimpleBinaryRelation() {
        // nothing to do
    }

    ~SimpleBinaryRelation() {
        // nothing to do
    }

    bool areRelated(const O& obj, const L& lab) {

        bool val = false;
        auto keyIter = binRel_.find(obj);

        if (keyIter != binRel_.end()) {

            auto& labels = keyIter->second;
            val = (std::find(labels.begin(), labels.end(), lab) != labels.end());

        }

        return val;

    }

    bool containsLabel(const L& lab) {

        bool notFound = true;

        for (auto iter = binRel_.begin(); notFound && iter != binRel_.end(); iter++) {
            notFound = std::find(iter->second.begin(), iter->second.end(), lab) == iter->second.end();
        }

        return !notFound;

    }

    bool containsObject(const O& obj) {
        return binRel_.find(obj) != binRel_.end();
    }

    std::vector<L> getLabelsOf(const O& obj) {

        std::vector<L> labels;
        auto keyIter = binRel_.find(obj);

        if (keyIter != binRel_.end()) {
            labels = keyIter->second;
        }

        return labels;

    }

    std::vector<O> getObjectsOf(const L& lab) {

        std::vector<O> objects;

        for (auto iter = binRel_.begin(); iter != binRel_.end(); iter++) {

            if (std::find(iter->second.begin(), iter->second.end(), lab) != iter->second.end()) {
                objects.push_back(iter->first);
            }

        }

        return objects;

    }


    unsigned long countPairs() {

        unsigned long sum = 0;

        for (auto iter = binRel_.begin(); iter != binRel_.end(); iter++) {
            sum += iter->second.size();
        }

        return sum;

    }

    unsigned long countLabelsOf(const O& obj) {

        unsigned long numLabels = 0;
        auto keyIter = binRel_.find(obj);

        if (keyIter != binRel_.end()) {
            numLabels = keyIter->second.size();
        }

        return numLabels;

    }

    unsigned long countObjectsOf(const L& lab) {

        unsigned long sum = 0;

        for (auto iter = binRel_.begin(); iter != binRel_.end(); iter++) {
            sum += std::find(iter->second.begin(), iter->second.end(), lab) != iter->second.end();
        }

        return sum;

    }


    void add(const O& obj, const L& lab) {

        // at() not necessary because creating new entry for s should be done anyway if it is not already included
        auto& labels = binRel_[obj];
        auto iter = std::find(labels.begin(), labels.end(), lab);

        if (iter == labels.end()) { // add label only if it is not already there
            labels.push_back(lab);
        }

    }

    void remove(const O& obj, const L& lab) {

        auto keyIter = binRel_.find(obj);

        if (keyIter != binRel_.end()) {

            auto& labels = keyIter->second;
            auto valueIter = std::find(labels.begin(), labels.end(), lab);

            if (valueIter != labels.end()) {
                labels.erase(valueIter);
            }

        }

    }


    void addObject(const O& obj) {
        binRel_[obj];
    }

    void removeObject(const O& obj) {
        binRel_.erase(obj);
    }

    void addLabel(const L& lab) {
        // nothing to do
    }

    void removeLabel(const L& lab) {

        for (auto iter = binRel_.begin(); iter != binRel_.end(); iter++) {

            auto pos = std::find(iter->second.begin(), iter->second.end(), lab);

            if (pos != iter->second.end()) {
                iter->second.erase(pos);
            }

        }

    }

    void join(const SimpleBinaryRelation<O, L>& sbr) {

        for (auto iter = sbr.binRel_.begin(); iter != sbr.binRel_.end(); iter++) {

            auto& labs = binRel_[iter->first];

            for (auto labIter = iter->second.begin(); labIter != iter->second.end(); labIter++) {

                if (std::find(labs.begin(), labs.end(), *labIter) == labs.end()) { // add label only if it is not already there
                    labs.push_back(*labIter);
                }

            }
        }

    }


private:
    std::unordered_map<O, std::vector<L>> binRel_;

};

template<typename O, typename L>
class SimpleBinaryRelationSet {

public:
    SimpleBinaryRelationSet() {
        // nothing to do
    }

    ~SimpleBinaryRelationSet() {
        // nothing to do
    }

    bool areRelated(const O& obj, const L& lab) {

        auto keyIter = binRel_.find(obj);

        return (keyIter != binRel_.end()) && (keyIter->second.find(lab) != keyIter->second.end());

    }

    bool containsLabel(const L& lab) {

        bool notFound = true;

        for (auto iter = binRel_.begin(); notFound && iter != binRel_.end(); iter++) {
            notFound = iter->second.find(lab) == iter->second.end();
        }

        return !notFound;

    }

    bool containsObject(const O& obj) {
        return binRel_.find(obj) != binRel_.end();
    }

    std::vector<L> getLabelsOf(const O& obj) {

        std::vector<L> res;
        auto keyIter = binRel_.find(obj);

        if (keyIter != binRel_.end()) {
            res = std::vector<L>(keyIter->second.begin(), keyIter->second.end());
        }

        return res;

    }

    std::vector<O> getObjectsOf(const L& lab) {

        std::vector<O> objects;

        for (auto iter = binRel_.begin(); iter != binRel_.end(); iter++) {

            if (iter->second.find(lab) != iter->second.end()) {
                objects.push_back(iter->first);
            }

        }

        return objects;

    }


    unsigned long countPairs() {

        unsigned long sum = 0;

        for (auto iter = binRel_.begin(); iter != binRel_.end(); iter++) {
            sum += iter->second.size();
        }

        return sum;

    }

    unsigned long countLabelsOf(const O& obj) {

        unsigned long numLabels = 0;
        auto keyIter = binRel_.find(obj);

        if (keyIter != binRel_.end()) {
            numLabels = keyIter->second.size();
        }

        return numLabels;

    }

    unsigned long countObjectsOf(const L& lab) {

        unsigned long sum = 0;

        for (auto iter = binRel_.begin(); iter != binRel_.end(); iter++) {
            sum += iter->second.find(lab) != iter->second.end();
        }

        return sum;

    }


    void add(const O& obj, const L& lab) {
        binRel_[obj].insert(lab);
    }

    void remove(const O& obj, const L& lab) {

        auto keyIter = binRel_.find(obj);

        if (keyIter != binRel_.end()) {
            keyIter->second.erase(lab);
        }

    }


    void addObject(const O& obj) {
        binRel_[obj];
    }

    void removeObject(const O& obj) {
        binRel_.erase(obj);
    }

    void addLabel(const L& lab) {
        // nothing to do
    }

    void removeLabel(const L& lab) {

        for (auto iter = binRel_.begin(); iter != binRel_.end(); iter++) {
            iter->second.erase(lab);
        }

    }

    void join(const SimpleBinaryRelationSet<O, L>& sbr) {

        for (auto iter = sbr.binRel_.begin(); iter != sbr.binRel_.end(); iter++) {

            auto& labs = binRel_[iter->first];
            labs.insert(iter->second.begin(), iter->second.end());

        }

    }


private:
    std::unordered_map<O, std::unordered_set<L>> binRel_;

};


/**
 * Variant of SimpleBinaryRelation(Set) with an optimisations for column-based queries.
 * Instead of iterating over all rows, each column (label) can be processed by following vertical links
 * pointing to the next row containing the label.
 */
template<typename O, typename L>
class TwoLinkBinaryRelation {

    struct Entry;

    typedef std::vector<Entry> LinkedRow;

    struct Entry {

        L value;
        std::pair<const O, LinkedRow>* next;

        Entry() {

        }

        Entry(L v, std::pair<const O, LinkedRow>* n) {

            value = v;
            next = n;

        }

        Entry(L v) {

            value = v;
            next = 0;

        }

        bool operator==(const Entry& e) {
            return value == e.value;
        }

    };

public:
    TwoLinkBinaryRelation() {
        // nothing to do
    }

    ~TwoLinkBinaryRelation() {
        // nothing to do
    }

    bool areRelated(const O& obj, const L& lab) {

        bool val = false;
        auto keyIter = binRel_.find(obj);

        if (keyIter != binRel_.end()) {

            auto& labels = keyIter->second;
            val = (std::find(labels.begin(), labels.end(), lab) != labels.end());

        }

        return val;

    }

    bool containsLabel(const L& lab) {
        return labels_.find(lab) != labels_.end();
    }

    bool containsObject(const O& obj) {
        return binRel_.find(obj) != binRel_.end();
    }

    std::vector<L> getLabelsOf(const O& obj) {

        std::vector<L> res;
        auto keyIter = binRel_.find(obj);

        if (keyIter != binRel_.end()) {

            for (auto iter = keyIter->second.begin(); iter != keyIter->second.end(); iter++) {
                res.push_back(iter->value);
            }

        }

        return res;

    }

    std::vector<O> getObjectsOf(const L& lab) {

        std::vector<O> objects;
        auto labIter = labels_.find(lab);
        auto entry = (labIter != labels_.end()) ? Entry(lab, labIter->second) : Entry(lab);

        while (entry.next != 0) {

            objects.push_back(entry.next->first);
            auto& nextRow = binRel_[entry.next->first];
            entry = *(std::find(nextRow.begin(), nextRow.end(), entry));

        }

        return objects;

    }


    unsigned long countPairs() {

        unsigned long sum = 0;

        for (auto iter = binRel_.begin(); iter != binRel_.end(); iter++) {
            sum += iter->second.size();
        }

        return sum;

    }

    unsigned long countLabelsOf(const O& obj) {

        unsigned long numLabels = 0;
        auto keyIter = binRel_.find(obj);

        if (keyIter != binRel_.end()) {
            numLabels = keyIter->second.size();
        }

        return numLabels;

    }

    unsigned long countObjectsOf(const L& lab) {

        unsigned long sum = 0;
        auto labIter = labels_.find(lab);
        auto entry = (labIter != labels_.end()) ? Entry(lab, labIter->second) : Entry(lab);

        while (entry.next != 0) {

            sum++;
            auto& nextRow = binRel_[entry.next->first];
            entry = *(std::find(nextRow.begin(), nextRow.end(), entry));

        }

        return sum;

    }


    void add(const O& obj, const L& lab) {

        auto labIter = labels_.find(lab);
        auto& row = binRel_[obj];

        if (labIter != labels_.end()) { // add to already existing column

            auto e = Entry(lab, labIter->second);
            auto iter = std::find(row.begin(), row.end(), e);

            if (iter == row.end()) { // only "redirect iterator" if insertion is actually done
                labels_[lab] = &(*(binRel_.find(obj)));
                row.push_back(e);
            }

        } else { // open new column

            auto e = Entry(lab);
            labels_[lab] = &(*(binRel_.find(obj)));
            row.push_back(e);

        }

    }

    void remove(const O& obj, const L& lab) {

        auto keyIter = binRel_.find(obj);
        auto labIter = labels_.find(lab);

        if (keyIter != binRel_.end() && labIter != labels_.end()) {

            auto entry = Entry(lab, labIter->second);
            auto prevEntry = entry;

            while (entry.next != 0) {

                if (entry.next->first == obj) {

                    auto& nextRow = binRel_[entry.next->first];
                    auto iter = std::find(nextRow.begin(), nextRow.end(), entry);

                    if (prevEntry.next->first == entry.next->first) { // remove first element of vertical chain
                        labels_[lab] = iter->next;
                    } else {

                        auto& prevRow = binRel_[prevEntry.next->first];
                        auto valueIter = std::find(prevRow.begin(), prevRow.end(), entry);
                        if (valueIter != prevRow.end()) {
                            *valueIter = Entry(lab, iter->next);
                        }

                    }

                    binRel_[entry.next->first].erase(iter);

                    return;

                }

                prevEntry = entry;
                entry = *(std::find(binRel_[entry.next->first].begin(), binRel_[entry.next->first].end(), entry));

            }

        }

    }


    void addObject(const O& obj) {
        binRel_[obj];
    }

    void removeObject(const O& obj) {
        binRel_.erase(obj);
    }

    void addLabel(const L& lab) {
        // nothing to do?
    }

    void removeLabel(const L& lab) {

        auto labIter = labels_.find(lab);

        if (labIter != labels_.end()) {

            auto entry = Entry(lab, labIter->second);
            labels_.erase(labIter);

            while (entry.next != 0) {

                auto iter = std::find(binRel_[entry.next->first].begin(), binRel_[entry.next->first].end(), entry);
                auto oldNext = entry.next->first;

                entry = *iter;
                binRel_[oldNext].erase(iter);

            }

        }

    }

    void join(const TwoLinkBinaryRelation& sbr) {

        for (auto iter = sbr.binRel_.begin(); iter != sbr.binRel_.end(); iter++) {

            for (auto labIter = iter->second.begin(); labIter != iter->second.end(); labIter++) {
                add(iter->first, labIter->value);
            }

        }

    }

private:
    std::unordered_map<L, std::pair<const O, LinkedRow>*> labels_;
    std::unordered_map<O, LinkedRow> binRel_;

};

// Note: Currently, this data structure allows only types with std::hash and std::equal_to methods
// for the labels (template parameter L).
template<typename O, typename L>
class TwoLinkBinaryRelationSet {

    struct Entry;

    struct EntryHash {

        std::hash<L> hf;

        size_t operator()(const Entry& e) const {
            return hf(e.value);
        }

    };

    struct EntryEqual {

        std::equal_to<L> ef;

        bool operator()(const Entry& e1, const Entry& e2) const {
            return ef(e1.value, e2.value);
        }

    };

    typedef std::unordered_set<Entry, EntryHash, EntryEqual> LinkedRow;

    struct Entry {

        L value;
        std::pair<const O, LinkedRow>* next;

        Entry() {

        }

        Entry(L v, std::pair<const O, LinkedRow>* n) {

            value = v;
            next = n;

        }

        Entry(L v) {

            value = v;
            next = 0;

        }

    };

public:
    TwoLinkBinaryRelationSet() {
        // nothing to do
    }

    ~TwoLinkBinaryRelationSet() {
        // nothing to do
    }

    bool areRelated(const O& obj, const L& lab) {

        auto keyIter = binRel_.find(obj);

        return (keyIter != binRel_.end()) && (keyIter->second.find(Entry(lab)) != keyIter->second.end());

    }

    bool containsLabel(const L& lab) {
        return labels_.find(lab) != labels_.end();
    }

    bool containsObject(const O& obj) {
        return binRel_.find(obj) != binRel_.end();
    }

    std::vector<L> getLabelsOf(const O& obj) {

        std::vector<L> res;
        auto keyIter = binRel_.find(obj);

        if (keyIter != binRel_.end()) {

            for (auto iter = keyIter->second.begin(); iter != keyIter->second.end(); iter++) {
                res.push_back(iter->value);
            }

        }

        return res;

    }

    std::vector<O> getObjectsOf(const L& lab) {

        std::vector<O> objects;
        auto labIter = labels_.find(lab);
        auto entry = (labIter != labels_.end()) ? Entry(lab, labIter->second) : Entry(lab);

        while (entry.next != 0) {

            objects.push_back(entry.next->first);
            entry = *(binRel_[entry.next->first].find(entry));

        }

        return objects;

    }


    unsigned long countPairs() {

        unsigned long sum = 0;

        for (auto iter = binRel_.begin(); iter != binRel_.end(); iter++) {
            sum += iter->second.size();
        }

        return sum;

    }

    unsigned long countLabelsOf(const O& obj) {

        unsigned long numLabels = 0;
        auto keyIter = binRel_.find(obj);

        if (keyIter != binRel_.end()) {
            numLabels = keyIter->second.size();
        }

        return numLabels;

    }

    unsigned long countObjectsOf(const L& lab) {

        unsigned long sum = 0;
        auto labIter = labels_.find(lab);
        auto entry = (labIter != labels_.end()) ? Entry(lab, labIter->second) : Entry(lab);

        while (entry.next != 0) {

            sum++;
            entry = *(binRel_[entry.next->first].find(entry));

        }

        return sum;

    }


    void add(const O& obj, const L& lab) {

        auto labIter = labels_.find(lab);

        if (labIter != labels_.end()) { // add to already existing column

            if (binRel_[obj].insert(Entry(lab, labIter->second)).second) { // only "redirect iterator" if insertion is actually done
                labels_[lab] = &(*(binRel_.find(obj)));
            }

        } else { // open new column

            binRel_[obj].insert(Entry(lab));
            labels_[lab] = &(*(binRel_.find(obj)));

        }

    }

    void remove(const O& obj, const L& lab) {

        auto keyIter = binRel_.find(obj);
        auto labIter = labels_.find(lab);

        if (keyIter != binRel_.end() && labIter != labels_.end()) {

            auto entry = Entry(lab, labIter->second);
            auto prevEntry = entry;

            while (entry.next != 0) {

                if (entry.next->first == obj) {

                    auto iter = binRel_[entry.next->first].find(entry);

                    if (prevEntry.next->first == entry.next->first) { // remove first element of vertical chain
                        labels_[lab] = iter->next;
                    } else {

                        binRel_[prevEntry.next->first].erase(entry);
                        binRel_[prevEntry.next->first].insert(Entry(lab, iter->next));

                    }

                    binRel_[entry.next->first].erase(iter);

                    return;

                }

                prevEntry = entry;
                entry = *(binRel_[entry.next->first].find(entry));

            }

        }

    }


    void addObject(const O& obj) {
        binRel_[obj];
    }

    void removeObject(const O& obj) {
        binRel_.erase(obj);
    }

    void addLabel(const L& lab) {
        // nothing to do?
    }

    void removeLabel(const L& lab) {

        auto labIter = labels_.find(lab);

        if (labIter != labels_.end()) {

            auto entry = Entry(lab, labIter->second);
            labels_.erase(labIter);

            while (entry.next != 0) {

                auto iter = binRel_[entry.next->first].find(entry);
                auto oldNext = entry.next->first;

                entry = *iter;
                binRel_[oldNext].erase(iter);

            }

        }

    }

    void join(const TwoLinkBinaryRelationSet& sbr) {

        for (auto iter = sbr.binRel_.begin(); iter != sbr.binRel_.end(); iter++) {

            for (auto labIter = iter->second.begin(); labIter != iter->second.end(); labIter++) {
                add(iter->first, labIter->value);
            }

        }

    }

private:
    std::unordered_map<L, std::pair<const O, LinkedRow>*> labels_;
    std::unordered_map<O, LinkedRow> binRel_;

};


// maps sequence segments to amplicon 'ids' (represented by their indices within the AmpliconPool)
//typedef SimpleBinaryRelation<std::string, numSeqs_t> InvertedIndex; // use with iterative index
typedef TwoLinkBinaryRelation<std::string, numSeqs_t> InvertedIndex; // use with full index


/*
 * Simple representation of matches (pairs of similar amplicons).
 * Stores the matches grouped by their (edit) distance;
 */
template<typename K, typename T>
class SimpleMatches {

public:
    SimpleMatches() {
        // nothing to do
    }

    ~SimpleMatches() {
        // nothing to do
    }


    // return whether the pair (obj,lab) is recorded as a pair
    bool contains(const T& t1, const T& t2) {

        bool found = false;

        for (auto iter = matches_.begin(); iter != matches_.end() && !found; iter++) {
            found = iter->second.areRelated(t1, t2);
        }

        return found;

    }

    // add a match
    void add(const T& t1, const T& t2, const K& d) {

        matches_[d].add(t1, t2);
        matches_[d].add(t2, t1); // only for symmetry

    }

    // remove (a) match
    void remove(const T& t1, const T& t2) {

        for (auto iter = matches_.begin(); iter != matches_.end(); iter++) {
            iter->second.remove(t1, t2);
            iter->second.remove(t2, t1); // only for symmetry
        }

    }


    // return all matches of the specified element
    std::vector<T> getMatchesOf(const T& t) {

        std::vector<T> res;
        std::vector<T> m;

        for (auto iter = matches_.begin(); iter != matches_.end(); iter++) {

            m = iter->second.getLabelsOf(t);
            res.insert(res.end(), m.begin(), m.end());

        }

        return res;

    }


    // return whether the pair (obj,lab) is recorded as a pair (1st component) and, if yes, what is their distance (2nd component)
    std::pair<bool, K> containsAugmented(const T& t1, const T& t2) {

        bool found = false;
        K d = 0;

        for (auto iter = matches_.begin(); iter != matches_.end() && !found; iter++) {
            found = iter->second.areRelated(t1, t2);
            d = iter->first;
        }

        return std::make_pair(found, d);

    }

    // return all matches (and their distances) of the specified element
    std::vector<std::pair<T, K>> getMatchesOfAugmented(const T& t) {

        std::vector<std::pair<T, K>> res;
        std::vector<T> m;

        for (auto iter = matches_.begin(); iter != matches_.end(); iter++) {

            m = iter->second.getLabelsOf(t);

            for (auto pIter = m.begin(); pIter != m.end(); pIter++) {
                res.push_back(std::make_pair(*pIter, iter->first));
            }

        }

        return res;

    }

    // count the total number of matches
    unsigned long countMatches() {

        unsigned long sum = 0;

        for (auto iter = matches_.begin(); iter != matches_.end(); iter++) {
            sum += iter->second.countPairs();
        }

//        return sum;
        return sum / 2; // only for symmetry

    }

    // count the number of matches of the specified element
    unsigned long countMatchesOf(const T& t) {

        unsigned long sum = 0;

        for (auto iter = matches_.begin(); iter != matches_.end(); iter++) {
            sum += iter->second.countLabelsOf(t);
        }

        return sum;

    }

    void join(const SimpleMatches<T, T>& sm) {

        for (auto iter = sm.matches_.begin(); iter != sm.matches_.end(); iter++) {
            matches_[iter->first].join(iter->second);
        }

    }


    bool syncContains(const T& t1, const T& t2) {

        std::lock_guard<std::mutex> lock(mtx_);
        return contains(t1, t2);

    }

    void syncAdd(const T& t1, const T& t2, const K& d) {

        std::lock_guard<std::mutex> lock(mtx_);
        add(t1, t2, d);

    }

    void syncRemove(const T& t1, const T& t2) {

        std::lock_guard<std::mutex> lock(mtx_);
        remove(t1, t2);

    }

    std::vector<T> syncGetMatchesOf(const T& t) {

        std::lock_guard<std::mutex> lock(mtx_);
        return getMatchesOf(t);

    }

    std::pair<bool, K> syncContainsAugmented(const T& t1, const T& t2) {

        std::lock_guard<std::mutex> lock(mtx_);
        return containsAugmented(t1, t2);

    }

    std::vector<std::pair<T, K>> syncGetMatchesOfAugmented(const T& t) {

        std::lock_guard<std::mutex> lock(mtx_);
        return getMatchesOfAugmented(t);

    }

    unsigned long syncCountMatches() {

        std::lock_guard<std::mutex> lock(mtx_);
        return countMatches();

    }

    unsigned long syncCountMatchesOf(const T& t) {

        std::lock_guard<std::mutex> lock(mtx_);
        return countMatchesOf(t);

    }

    void syncJoin(const SimpleMatches<T, T>& sm) { // assumption: no concurrent accesses on sm

        std::lock_guard<std::mutex> lock(mtx_);
        join(sm);

    }


private:
    std::map<K, SimpleBinaryRelation<T, T>> matches_;
    std::mutex mtx_;

};


// minimum interface: contains, add, remove
typedef SimpleMatches<lenSeqs_t, numSeqs_t> Matches;

}

#endif //SCT_PJ_RELATION_HPP