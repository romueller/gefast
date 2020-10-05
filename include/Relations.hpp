/*
 * GeFaST
 *
 * Copyright (C) 2016 - 2020 Robert Mueller
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

#ifndef GEFAST_RELATIONS_HPP
#define GEFAST_RELATIONS_HPP

#include <algorithm>

#include "Utility.hpp"

namespace GeFaST {

    /*
     * Variant of a binary relation tailored to the case when there is only one object per label.
     * Instead of iterating over all rows, each column (label) can be processed by following the vertical link
     * pointing to the row containing the label.
     *
     * NOTE: The labels of the added elements have to be monotonically increasing,
     * so that the labels_ vector is sorted without any further action.
     */
    template<typename O, typename L, typename H = std::hash<O>, typename P = std::equal_to<O>>
    class LabelLinkBinaryRelation {

    public:
        LabelLinkBinaryRelation() = default;

        ~LabelLinkBinaryRelation() = default;

        LabelLinkBinaryRelation(const LabelLinkBinaryRelation& other) : bin_rel_(other.bin_rel_) { // copy constructor

            // iterate over the already sorted labels in the other instance
            // and find the one related object in the copy of the binary relation in the new instance
            // to obtain the pointer into the copied binary relation
            // and insert the new pair into labels_ of the new instance
            for (auto iter = other.labels_.begin(); iter != other.labels_.end(); iter++) {

                L lab = iter->first;
                const O obj = iter->second->first;

                labels_.emplace_back(lab, &(*(bin_rel_.find(obj))));

            }

        }

        LabelLinkBinaryRelation(LabelLinkBinaryRelation&& other) = default; // move constructor

        LabelLinkBinaryRelation& operator=(const LabelLinkBinaryRelation& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            // nothing to do

            // copy new resources
            bin_rel_ = other.bin_rel_;

            // iterate over the already sorted labels in the other instance
            // and find the one related object in the copy of the binary relation in the new instance
            // to obtain the pointer into the copied binary relation
            // and insert the new pair into labels_ of the new instance
            for (auto iter = other.labels_.begin(); iter != other.labels_.end(); iter++) {

                L lab = iter->first;
                const O obj = iter->second->first;

                labels_.emplace_back(lab, &(*(bin_rel_.find(obj))));

            }


        }

        LabelLinkBinaryRelation& operator=(LabelLinkBinaryRelation&& other) = default; // move assignment operator

        bool are_related(const O& obj, const L& lab) {

            auto lab_iter = std::lower_bound(labels_.begin(), labels_.end(), lab, cmp_);

            return (lab_iter != labels_.end()) && (lab_iter->first == lab)
                    && (lab_iter->second != 0) && (obj == lab_iter->second->first);

        }

        bool contains_label(const L& lab) {

            auto iter = std::lower_bound(labels_.begin(), labels_.end(), lab, cmp_);

            return (iter != labels_.end()) && (iter->first == lab) && (iter->second != 0);

        }

        bool contains_object(const O& obj) {
            return bin_rel_.find(obj) != bin_rel_.end();
        }

        std::vector<L> get_labels_of(const O& obj) {

            auto key_iter = bin_rel_.find(obj);

            return (key_iter != bin_rel_.end()) ? key_iter->second : std::vector<L>();

        }

        void add_label_counts_of(const O& obj, std::vector<numSeqs_t>& cand_cnts) {

            auto key_iter = bin_rel_.find(obj);

            if (key_iter != bin_rel_.end()) {
                cand_cnts.insert(cand_cnts.end(), key_iter->second.begin(), key_iter->second.end());
            }

        }

        std::vector<O> get_objects_of(const L& lab) {

            std::vector<O> objects;
            auto lab_iter = std::lower_bound(labels_.begin(), labels_.end(), lab, cmp_);
            if ((lab_iter != labels_.end()) && (lab_iter->first == lab) && (lab_iter->second != 0)) {
                objects.push_back(lab_iter->second->first);
            }

            return objects;

        }


        unsigned long count_pairs() {

            unsigned long sum = 0;

            for (auto iter = bin_rel_.begin(); iter != bin_rel_.end(); iter++) {
                sum += iter->second.size();
            }

            return sum;

        }

        unsigned long count_labels_of(const O& obj) {

            auto key_iter = bin_rel_.find(obj);

            return (key_iter != bin_rel_.end()) ? key_iter->second.size() : 0;

        }

        unsigned long count_objects_of(const L& lab) {

            auto lab_iter = std::lower_bound(labels_.begin(), labels_.end(), lab, cmp_);

            return ((lab_iter != labels_.end()) && (lab_iter->first == lab) && (lab_iter->second != 0));

        }


        void add(const O& obj, const L& lab) {

            bin_rel_[obj].push_back(lab);
            labels_.emplace_back(lab, &(*(bin_rel_.find(obj))));

        }

        void remove(const O& obj, const L& lab) {

            auto lab_iter = std::lower_bound(labels_.begin(), labels_.end(), lab, cmp_);

            if ((lab_iter != labels_.end()) && (lab_iter->first == lab)
                && (lab_iter->second != 0) && (obj == lab_iter->second->first)) {

                auto iter = std::find(lab_iter->second->second.begin(), lab_iter->second->second.end(), lab);
                lab_iter->second->second.erase(iter);
                lab_iter->second = 0;

            }

        }

        void remove_label(const L& lab) {

            auto lab_iter = std::lower_bound(labels_.begin(), labels_.end(), lab, cmp_);

            if ((lab_iter != labels_.end()) && (lab_iter->first == lab) && (lab_iter->second != 0)) {

                if (lab_iter->second != 0) {

                    auto iter = std::find(lab_iter->second->second.begin(), lab_iter->second->second.end(), lab);
                    lab_iter->second->second.erase(iter);
                    lab_iter->second = 0;

                }

//            labels_.erase(lab_iter);

            }

        }

    private:
        std::vector<std::pair<L, std::pair<const O, std::vector<L>>*>> labels_; // contained labels with pointer to row of occurrence
        std::unordered_map<O, std::vector<L>, H, P> bin_rel_; // binary relation between objects of type O and labels of type L

        struct CompareLabels {

            bool operator() (const std::pair<L, std::pair<const O, std::vector<L>>*>& lhs, const L rhs) {
                return lhs.first < rhs;
            }

        } cmp_;

    };

}

#endif //GEFAST_RELATIONS_HPP
