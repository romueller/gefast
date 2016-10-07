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

#ifndef SCT_PJ_RELATION_HPP
#define SCT_PJ_RELATION_HPP

#include <algorithm>
#include <mutex>
#include <unordered_map>
#include <vector>

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

            auto labels = keyIter->second;
            auto valueIter = std::find(labels.begin(), labels.end(), lab);

            val = (valueIter != labels.end());

        }

        return val;

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
//        auto labels = &(binRel_[obj]);
//        auto iter = std::find(labels->begin(), labels->end(), lab);
//
//        if (iter == labels->end()) { // add label only if it is not already there
//            labels->push_back(lab);
//        }
//
        std::vector<L>& labels = binRel_[obj];
        auto iter = std::find(labels.begin(), labels.end(), lab);

        if (iter == labels.end()) { // add label only if it is not already there
            labels.push_back(lab);
        }

    }

    void remove(const O& obj, const L& lab) {

        auto keyIter = binRel_.find(obj);

        if (keyIter != binRel_.end()) {

            auto labels = keyIter->second;
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
        //nothing to do
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
    std::unordered_map<O,std::vector<L>> binRel_;

};


// maps sequence segments to amplicon 'ids' (represented by their indices within the AmpliconPool)
typedef SimpleBinaryRelation<std::string, numSeqs_t> InvertedIndex;


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
        matches_[d].add(t2, t1); //NEW for symmetry

    }

    // remove (a) match
    void remove(const T& t1, const T& t2) {

        for (auto iter = matches_.begin(); iter != matches_.end(); iter++) {
            iter->second.remove(t1, t2);
            iter->second.remove(t2, t1); //NEW for symmetry
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
        return sum / 2; //NEW for symmetry

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

    void syncJoin(const SimpleMatches<T, T>& sm) {//assumption: no concurrent accesses on sm

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