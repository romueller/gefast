/*
 * GeFaST
 *
 * Copyright (C) 2016 - 2021 Robert Mueller
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

#include "../../include/space/FlexibleRelations.hpp"

namespace GeFaST {

    /* === RankedAscendingLabels === */

    // to uphold sortedness while also allowing entries to be marked as removed without extra space,
    // removed entries receive a negative sign and comparisons use the absolute value
    struct CompareAbsoluteValues {

        bool operator() (const long long lhs, const long long rhs) {
            return std::llabs(lhs) < std::llabs(rhs);
        }

    } cmp_abs_values;

    RankedAscendingLabels::RankedAscendingLabels() {

        labels_ = nullptr;
        size_ = 0;
        capa_= 0;

    }

    RankedAscendingLabels::RankedAscendingLabels(const numSeqs_t capacity) {

        labels_ = new long long[capacity];
        size_ = 0;
        capa_= capacity;

    }

    RankedAscendingLabels::RankedAscendingLabels(std::vector<numSeqs_t>& labels) {

        capa_ = labels.size();
        size_ = labels.size();

        labels_ = new long long[labels.size()];
        for (numSeqs_t i = 0; i < size_; i++) {
            labels_[i] = static_cast<long long>(labels[i]) + 1;
        }

    }

    RankedAscendingLabels::~RankedAscendingLabels() {
        delete[] labels_;
    }

    RankedAscendingLabels::RankedAscendingLabels(const RankedAscendingLabels& other) {

        size_ = other.size_;
        capa_ = other.capa_;

        labels_ = new long long[capa_];
        for (numSeqs_t i = 0; i < size_; i++) {
            labels_[i] = other.labels_[i];
        }

    }

    RankedAscendingLabels::RankedAscendingLabels(RankedAscendingLabels&& other) noexcept {

        size_ = other.size_;
        capa_ = other.capa_;

        labels_ = other.labels_;

    }

    RankedAscendingLabels& RankedAscendingLabels::operator=(const RankedAscendingLabels& other) {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] labels_;

        // copy new resources
        size_ = other.size_;
        capa_ = other.capa_;

        labels_ = new long long[capa_];
        for (numSeqs_t i = 0; i < size_; i++) {
            labels_[i] = other.labels_[i];
        }

        return *this;

    }

    RankedAscendingLabels& RankedAscendingLabels::operator=(RankedAscendingLabels&& other) noexcept {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // remove old resources
        delete[] labels_;

        // copy / transfer new resources
        size_ = other.size_;
        capa_ = other.capa_;

        labels_ = other.labels_; other.labels_ = nullptr;

        return *this;

    }

    RankedAscendingLabels* RankedAscendingLabels::clone() const {
        return new RankedAscendingLabels(*this);
    }

    numSeqs_t RankedAscendingLabels::add(const numSeqs_t lab) {

        labels_[size_] = static_cast<long long>(lab) + 1;

        return size_++;

    }

    numSeqs_t RankedAscendingLabels::unrank(const numSeqs_t r) const {
        return std::llabs(labels_[r]) - 1;
    }

    bool RankedAscendingLabels::contains(const numSeqs_t lab) const {

        auto iter = std::lower_bound(labels_, labels_ + size_, lab + 1, cmp_abs_values);

        return (iter != (labels_ + size_)) && (lab + 1 == *iter);

    }

    bool RankedAscendingLabels::contains_rank(const numSeqs_t r) const {
        return (r < size_) && (labels_[r] > 0);
    }

    numSeqs_t RankedAscendingLabels::remove(const numSeqs_t lab) {

        auto iter = std::lower_bound(labels_, labels_ + size_, lab + 1, cmp_abs_values);

        // currently no check (iter != (labels_ + size_)) needed, because only initially existing values are removed
        *iter = -std::llabs(*iter);

        return iter - labels_;

    }

    numSeqs_t RankedAscendingLabels::size() const {
        return size_;
    }

    size_t RankedAscendingLabels::size_in_bytes() const {
        return sizeof(RankedAscendingLabels) // RankedAscendingLabels itself
            + sizeof(long long) * capa_; // labels_
    }

}
