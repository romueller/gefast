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

#include "../../../include/space/k2trees/Utility.hpp"

namespace GeFaST {

    bool is_all_zero(const std::vector<bool>& v) {
        return std::find(v.begin(), v.end(), true) == v.end();
    }


    NaiveBitVector::NaiveBitVector(size_t size) : NaiveBitVector(size, false) {
        // nothing else to do
    }

    NaiveBitVector::NaiveBitVector(size_t size, bool value) : bits_(std::vector<bool>(size, value)) {
        // nothing else to do
    }

    NaiveBitVector::NaiveBitVector(const std::vector<bool>& vec) : bits_(vec) {
        // nothing else to do
    }

    size_t NaiveBitVector::size() const {
        return bits_.size();
    }

    bool NaiveBitVector::empty() const {
        return bits_.empty();
    }

    void NaiveBitVector::set(size_t i, bool val) {
        bits_[i] = val;
    }

    bool NaiveBitVector::get(size_t i) const {
        return bits_[i];
    }

    bool NaiveBitVector::operator[](size_t i) const {
        return bits_[i];
    }

    void NaiveBitVector::append(const std::vector<bool>& bits) {
        for (bool bit : bits) bits_.push_back(bit);
    }

    size_t NaiveBitVector::set(size_t i, const std::vector<bool>& bits) {

        for (bool bit : bits) bits_[i++] = bit;
        return i;

    }

    bool NaiveBitVector::is_zero(size_t i, size_t len) {

        bool res = true;
        for (auto j = 0; res && j < len; j++) res &= !bits_[i + j];
        return res;

    }

    size_t NaiveBitVector::size_in_bytes() const {
        return sizeof(NaiveBitVector) // NaiveBitVector itself
            + bits_.capacity() / 8; // bits_, approximation
    }

}