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

#ifndef DACS_SPLITDACS_HPP
#define DACS_SPLITDACS_HPP

#include "SplitDacsBuilder.hpp"
#include "Utility.hpp"

namespace GeFaST {

    /*
     * Implementation of Directly Addressable Codes (DACs), a memory-efficient representation of
     * integer sequences using variable length coding.
     *
     * Based on the original DACs implementation:
     * DACs: Bringing Direct Access to Variable-Length Codes (Brisaboa et al., 2013)
     * http://lbd.udc.es/research/DACS/
     *
     * Template parameters (and extension to original implementation):
     *  - V = type of values to be stored (no longer fixed / restricted to 32 bits)
     *  - S = type underlying internal bit vectors (defines the 'word length', no longer fixed / restricted to 32 bits)
     *  - MAX_LEVELS = maximum allowed number of levels in the DACs representation (see below for more information)
     *  - FACTOR_RANK = parameter 'factor' of the rank data structure (controls its time / space trade-off)
     *
     * Our implementation covers the optimisation of the chunk lengths with and without a restricted number of levels.
     * Setting MAX_LEVELS to 0, the number of levels is optimised without limitation. This optimises the memory
     * consumption at the expense of time efficiency.
     *
     * Furthermore, our implementation simplifies the methods at several points and adds a second method to
     * extract consecutive elements of sequence (which is able to start from an index other than 0).
     *
     * Note that indices are zero-based in our implementation.
     *
     * Assumes that all integers are non-negative.
     */
    template<typename V, typename S, usmall_t MAX_LEVELS, usmall_t FACTOR_RANK>
    class SplitDacs {

    public:
        /*
         * Default constructor for "empty sequence" of length 0.
         */
        SplitDacs() {

            length_ = 0;
            num_levels_ = 0;
            accum_level_sizes_ = level_ranks_ = nullptr;
            extra_bits_ranks_ = nullptr;
            levels_ = nullptr;
            chunk_lengths_ = nullptr;
            min_values_ = nullptr;

        }

        // build SplitDacs directly from the integer sequence
        SplitDacs(const V* values, ularge_t length) {

            length_ = length;

            // determine number of levels and their chunk lengths
            Chunks chunks = (MAX_LEVELS == 0) ?
                            optimise_level_chunks<V, FACTOR_RANK>(values, length_)
                          : optimise_level_chunks_limited<V, FACTOR_RANK>(values, length_, MAX_LEVELS);
            chunk_lengths_ = chunks.chunk_lengths;
            num_levels_ = chunks.num_levels;

            // for each level, determine the minimum value that makes a chunk on the level necessary
            min_values_ = new V[num_levels_];
            accum_level_sizes_ = new ularge_t[num_levels_ + 1];
            V value = 0;
            V mult_val = 1;
            for (usmall_t i = 0; i < num_levels_; i++) {

                mult_val *= static_cast<V>(1) << chunk_lengths_[i];
                min_values_[i] = value;
                value += mult_val;

                accum_level_sizes_[i] = 0;

            }
            accum_level_sizes_[num_levels_] = 0;

            // for each level, determine how many chunks will be required
            for (ularge_t i = 0; i < length_; i++) {

                value = values[i];
                for (usmall_t j = 0; j < num_levels_; j++) {
                    accum_level_sizes_[j + 1] += (value >= min_values_[j]);
                }

            }

            // determine the total number of bits required for the chunks on all levels,
            // initialise the levels and accumulate the level sizes
            levels_ = new S*[num_levels_];
            for (usmall_t i = 0; i < num_levels_; i++) {

                levels_[i] = new S[(chunk_lengths_[i] * accum_level_sizes_[i + 1] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                accum_level_sizes_[i + 1] += accum_level_sizes_[i];

            }

            // set up cursors into the levels for filling tem
            ularge_t* level_cursors = new ularge_t[num_levels_]; // for each level, position (in bits) where to write the next chunk
            ularge_t* chunk_cursors = new ularge_t[num_levels_]; // for each level, position (in chunks) of the next chunk to be written

            for (usmall_t j = 0; j < num_levels_; j++) {

                level_cursors[j] = 0;
                chunk_cursors[j] = 0;

            }

            // initialise bit vector for Vbyte extra bits
            ularge_t len_extra_bits = accum_level_sizes_[num_levels_ - 1] + 1; // + 1 for first chunk on last level
            S* extra_bits = new S[(len_extra_bits + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
            for (ularge_t i = 0; i < (len_extra_bits + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); i++) {
                extra_bits[i] = 0;
            }

            // number of possibles values for a chunk on the different levels (used to write the chunks below)
            V* base = new V[num_levels_];
            for (usmall_t i = 0; i < num_levels_; i++) {
                base[i] = static_cast<V>(1) << chunk_lengths_[i];
            }

            // fill the levels and extra bits
            for (ularge_t i = 0; i < length_; i++) {

                value = values[i];

                for (ssmall_t j = num_levels_ - 1; j >= 0; j--) {

                    if (value >= min_values_[j]) { // find the highest required level for the current value

                        V new_value = value - min_values_[j]; // reduce the value

                        for (usmall_t k = 0; k <= j; k++) { // write the chunks in the respective levels

                            write_bits<V, S>(levels_[k], level_cursors[k], chunk_lengths_[k], new_value % base[k]);
                            level_cursors[k] += chunk_lengths_[k];
                            chunk_cursors[k]++;

                            new_value /= base[k];

                        }

                        // set the extra bit for the highest required level (unless it is the maximum level)
                        if (j < num_levels_ - 1) {
                            set_bit<S>(extra_bits, accum_level_sizes_[j] + chunk_cursors[j] - 1);
                        }

                        j = 0; // works as break only if j is a signed integer

                    }

                }

            }

            set_bit<S>(extra_bits, len_extra_bits - 1); // write the extra bit for the first chunk on the last level

            delete[] base;
            delete[] chunk_cursors;
            delete[] level_cursors;

            // set up rank data structures (receives ownership of extra_bits)
            extra_bits_ranks_ = new BitRank<S, FACTOR_RANK>(extra_bits, len_extra_bits, true);

            level_ranks_ = new ularge_t[num_levels_];
            level_ranks_[0] = 0;
            for (usmall_t j = 1; j < num_levels_; j++) {
                level_ranks_[j] = extra_bits_ranks_->rank(accum_level_sizes_[j] - 1);
            }

        }

        /*
         * Build SplitDacs by wrapping up its given members.
         */
        SplitDacs(ularge_t length, S** levels, usmall_t num_levels,
                  ularge_t* accum_level_sizes, BitRank<S, FACTOR_RANK>* extra_bits_ranks, ularge_t* level_ranks,
                  usmall_t* chunk_lengths, V* min_values) {

            length_ = length;
            levels_ = levels;
            num_levels_ = num_levels;
            accum_level_sizes_ = accum_level_sizes;
            extra_bits_ranks_ = extra_bits_ranks;
            level_ranks_ = level_ranks;
            chunk_lengths_ = chunk_lengths;
            min_values_ = min_values;

        }

        /*
         * Build SplitDacs from another instance and rearrange the elements according to vector 'index',
         * i.e. the i-th element of the new SplitDacs is other[index[i] - offset].
         */
        SplitDacs(SplitDacs& other, const std::vector<ularge_t>& index, V offset = 0) {

            /* copy members that do not change by reordering the elements */

            length_ = other.length_;
            num_levels_ = other.num_levels_;

            chunk_lengths_ = new usmall_t[num_levels_];
            for (usmall_t i = 0; i < num_levels_; i++) {
                chunk_lengths_[i] = other.chunk_lengths_[i];
            }

            accum_level_sizes_ = new ularge_t[num_levels_ + 1];
            for (usmall_t i = 0; i <= num_levels_; i++) {
                accum_level_sizes_[i] = other.accum_level_sizes_[i];
            }

            min_values_ = new V[num_levels_];
            for (usmall_t i = 0; i < num_levels_; i++) {
                min_values_[i] = other.min_values_[i];
            }


            /* fill the new SplitDacs with the rearranged elements */

            levels_ = new S*[num_levels_];
            for (usmall_t i = 0; i < num_levels_; i++) {
                levels_[i] = new S[(chunk_lengths_[i] * (accum_level_sizes_[i + 1] - accum_level_sizes_[i]) + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
            }

            // set up cursors into the levels for filling tem
            ularge_t* level_cursors = new ularge_t[num_levels_]; // for each level, position (in bits) where to write the next chunk
            ularge_t* chunk_cursors = new ularge_t[num_levels_]; // for each level, position (in chunks) of the next chunk to be written

            for (usmall_t j = 0; j < num_levels_; j++) {

                level_cursors[j] = 0;
                chunk_cursors[j] = 0;

            }

            // initialise bit vector for Vbyte extra bits
            ularge_t len_extra_bits = accum_level_sizes_[num_levels_ - 1] + 1; // + 1 for first chunk on last level
            S* extra_bits = new S[(len_extra_bits + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
            for (ularge_t i = 0; i < (len_extra_bits + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); i++) {
                extra_bits[i] = 0;
            }

            // number of possibles values for a chunk on the different levels (used to write the chunks below)
            V* base = new V[num_levels_];
            for (usmall_t i = 0; i < num_levels_; i++) {
                base[i] = static_cast<V>(1) << chunk_lengths_[i];
            }

            for (ularge_t i = 0; i < length_; i++) {

                V value = other[index[i] - offset];

                for (ssmall_t j = num_levels_ - 1; j >= 0; j--) {

                    if (value >= min_values_[j]) { // find the highest required level for the current value

                        V new_value = value - min_values_[j]; // reduce the value

                        for (usmall_t k = 0; k <= j; k++) { // write the chunks in the respective levels

                            write_bits<V, S>(levels_[k], level_cursors[k], chunk_lengths_[k], new_value % base[k]);
                            level_cursors[k] += chunk_lengths_[k];
                            chunk_cursors[k]++;

                            new_value /= base[k];

                        }

                        // set the extra bit for the highest required level (unless it is the maximum level)
                        if (j < num_levels_ - 1) {
                            set_bit<S>(extra_bits, accum_level_sizes_[j] + chunk_cursors[j] - 1);
                        }

                        j = 0; // works as break only if j is a signed integer

                    }

                }

            }

            set_bit<S>(extra_bits, len_extra_bits - 1); // write the extra bit for the first chunk on the last level

            delete[] base;
            delete[] chunk_cursors;
            delete[] level_cursors;

            // set up rank data structures (receives ownership of extra_bits)
            extra_bits_ranks_ = new BitRank<S, FACTOR_RANK>(extra_bits, len_extra_bits, true);

            level_ranks_ = new ularge_t[num_levels_];
            level_ranks_[0] = 0;
            for (usmall_t j = 1; j < num_levels_; j++) {
                level_ranks_[j] = extra_bits_ranks_->rank(accum_level_sizes_[j] - 1);
            }

        }

        /*
         * Builds SplitDacs instance from concatenation of multiple other SplitDacs instance.
         * The sequences are concatenated from left to right.
         * No reordering happens within each sequence.
         */
        explicit SplitDacs(const std::vector<SplitDacs<V, S, MAX_LEVELS, FACTOR_RANK>>& others) : SplitDacs() {

            SplitDacsDistributionBuilder<V, S, MAX_LEVELS, FACTOR_RANK>* dacs_builder;

            {

                std::map<ularge_t, ularge_t> counts;
                ularge_t total_length = 0;
                for (ularge_t i = 0; i < others.size(); i++) {

                    total_length += others[i].length();
                    for (ularge_t j = 0; j < others[i].length(); j++) {
                        counts[others[i][j]]++;
                    }

                }

                dacs_builder = new SplitDacsDistributionBuilder<V, S, MAX_LEVELS, FACTOR_RANK>(counts);

            }

            for (ularge_t i = 0; i < others.size(); i++) {
                for (ularge_t j = 0; j < others[i].length(); j++) {
                    dacs_builder->append(others[i][j]);
                }
            }

            auto dacs = dacs_builder->finalise();
            delete dacs_builder;

            *this = std::move(*dacs);
            delete dacs;

        }

        /*
         * Alternative constructor from multiple SplitDacs instances using the provided chunk length and capacity.
         */
        SplitDacs(const std::vector<SplitDacs<V, S, MAX_LEVELS, FACTOR_RANK>>& others, usmall_t chunk_length, ularge_t capacity, bool opt = true) : SplitDacs() {

            SplitDacsDirectBuilder<V, S, MAX_LEVELS, FACTOR_RANK> dacs_builder(chunk_length, capacity);
            for (ularge_t i = 0; i < others.size(); i++) {
                for (ularge_t j = 0; j < others[i].length(); j++) {
                    dacs_builder.append(others[i][j]);
                }
            }
            auto dacs = dacs_builder.finalise(opt);

            *this = std::move(*dacs);
            delete dacs;

        }

        ~SplitDacs() {

            for (usmall_t i = 0; i < num_levels_; i++) {
                delete[] levels_[i];
            }
            delete[] levels_;
            delete[] accum_level_sizes_;
            delete extra_bits_ranks_;
            delete[] level_ranks_;
            delete[] chunk_lengths_;
            delete[] min_values_;

        }

        SplitDacs(const SplitDacs& other) { // copy constructor

            length_ = other.length_;
            num_levels_ = other.num_levels_;

            accum_level_sizes_ = new ularge_t[num_levels_ + 1];
            std::copy(other.accum_level_sizes_, other.accum_level_sizes_ + num_levels_ + 1, accum_level_sizes_);

            extra_bits_ranks_ = other.extra_bits_ranks_->clone();

            level_ranks_ = new ularge_t[num_levels_];
            std::copy(other.level_ranks_, other.level_ranks_ + num_levels_, level_ranks_);

            chunk_lengths_ = new usmall_t[num_levels_];
            std::copy(other.chunk_lengths_, other.chunk_lengths_ + num_levels_, chunk_lengths_);

            min_values_ = new V[num_levels_];
            std::copy(other.min_values_, other.min_values_ + num_levels_, min_values_);

            levels_ = new S*[num_levels_];
            for (usmall_t i = 0; i < num_levels_; i++) {

                levels_[i] = new S[(chunk_lengths_[i] * (accum_level_sizes_[i + 1] - accum_level_sizes_[i]) + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                std::copy(other.levels_[i], other.levels_[i] + (chunk_lengths_[i] * (accum_level_sizes_[i + 1] - accum_level_sizes_[i]) + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), levels_[i]);

            }

        }

        SplitDacs(SplitDacs&& other) noexcept { // move constructor

            length_ = other.length_;
            num_levels_ = other.num_levels_; other.num_levels_ = 0;

            accum_level_sizes_ = other.accum_level_sizes_; other.accum_level_sizes_ = nullptr;
            extra_bits_ranks_ = other.extra_bits_ranks_; other.extra_bits_ranks_ = nullptr;
            level_ranks_ = other.level_ranks_; other.level_ranks_ = nullptr;
            chunk_lengths_ = other.chunk_lengths_; other.chunk_lengths_ = nullptr;
            min_values_ = other.min_values_; other.min_values_ = nullptr;
            levels_ = other.levels_; other.levels_ = nullptr;

        }

        SplitDacs& operator=(const SplitDacs& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            for (usmall_t i = 0; i < num_levels_; i++) {
                delete[] levels_[i];
            }
            delete[] levels_;
            delete[] accum_level_sizes_;
            delete extra_bits_ranks_;
            delete[] level_ranks_;
            delete[] chunk_lengths_;
            delete[] min_values_;

            // copy new resources
            length_ = other.length_;
            num_levels_ = other.num_levels_;

            accum_level_sizes_ = new ularge_t[num_levels_ + 1];
            std::copy(other.accum_level_sizes_, other.accum_level_sizes_ + num_levels_ + 1, accum_level_sizes_);

            extra_bits_ranks_ = other.extra_bits_ranks_->clone();

            level_ranks_ = new ularge_t[num_levels_];
            std::copy(other.level_ranks_, other.level_ranks_ + num_levels_, level_ranks_);

            chunk_lengths_ = new usmall_t[num_levels_];
            std::copy(other.chunk_lengths_, other.chunk_lengths_ + num_levels_, chunk_lengths_);

            min_values_ = new V[num_levels_];
            std::copy(other.min_values_, other.min_values_ + num_levels_, min_values_);

            levels_ = new S*[num_levels_];
            for (usmall_t i = 0; i < num_levels_; i++) {

                levels_[i] = new S[(chunk_lengths_[i] * (accum_level_sizes_[i + 1] - accum_level_sizes_[i]) + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                std::copy(other.levels_[i], other.levels_[i] + (chunk_lengths_[i] * (accum_level_sizes_[i + 1] - accum_level_sizes_[i]) + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), levels_[i]);

            }

            return *this;

        }

        SplitDacs& operator=(SplitDacs&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            for (usmall_t i = 0; i < num_levels_; i++) {
                delete[] levels_[i];
            }
            delete[] levels_;
            delete[] accum_level_sizes_;
            delete extra_bits_ranks_;
            delete[] level_ranks_;
            delete[] chunk_lengths_;
            delete[] min_values_;

            // copy / transfer new resources
            length_ = other.length_;
            num_levels_ = other.num_levels_; other.num_levels_ = 0;

            accum_level_sizes_ = other.accum_level_sizes_; other.accum_level_sizes_ = nullptr;
            extra_bits_ranks_ = other.extra_bits_ranks_; other.extra_bits_ranks_ = nullptr;
            level_ranks_ = other.level_ranks_; other.level_ranks_ = nullptr;
            chunk_lengths_ = other.chunk_lengths_; other.chunk_lengths_ = nullptr;
            min_values_ = other.min_values_; other.min_values_ = nullptr;
            levels_ = other.levels_; other.levels_ = nullptr;

            return *this;

        }

        SplitDacs* clone() const { // deep-copy clone method
            return new SplitDacs(*this);
        }

        /*
         * Extract the i-th integer of the sequence (index is zero-based).
         */
        V operator[](ularge_t i) const {

            // read first chunk of representation
            // (simplified by the convention that chunkLevelIndex[0] = bitLevelIndex[0] = 0)
            V partial_res = read_bits<V, S>(levels_[0], i * chunk_lengths_[0], chunk_lengths_[0]);

            // iterate over the higher levels (if necessary) to collect and combine the chunks of the representation
            ularge_t pos = i;
            ularge_t mult = chunk_lengths_[0];
            usmall_t j = 1;
            for (; j < num_levels_ && !extra_bits_ranks_->is_set(pos); j++) {

                // relative index in current level =
                //  relative rank on previous level
                //  - number of predecessors (in the sequence) having their last chunk on the previous level
                i -= (extra_bits_ranks_->rank(accum_level_sizes_[j - 1] + i) - level_ranks_[j - 1]);

                // read the chunk of the current level and add it to the partial result already computed
                partial_res += read_bits<V, S>(levels_[j], i * chunk_lengths_[j], chunk_lengths_[j]) << mult;

                // update the position of the next chunk (in the partial result) and of the extra bit
                mult += chunk_lengths_[j];
                pos = accum_level_sizes_[j] + i;

            }

            // undo the reduction
            return partial_res + min_values_[j - 1];

        }

        /*
         * Extract the first 'num' integers of the sequence.
         */
        V* extract(ularge_t num) const {

            if (num > length_ || length_ == 0) return 0;

            // initialise cursors to the first bit and chunk of each level
            ularge_t* level_cursors = new ularge_t[num_levels_];
            ularge_t* chunk_cursors = new ularge_t[num_levels_];
            for (usmall_t j = 0; j < num_levels_; j++) {

                level_cursors[j] = 0;
                chunk_cursors[j] = 0;

            }

            // decompress the integers one by one (while using and updating the cursors)
            V* list = new V[num];
            for (ularge_t i = 0; i < num; i++) {

                // read first chunk of representation
                V partial_res = read_bits<V, S>(levels_[0], level_cursors[0], chunk_lengths_[0]);

                ularge_t mult = chunk_lengths_[0];
                level_cursors[0] += chunk_lengths_[0];
                ularge_t p = chunk_cursors[0]++;

                // iterate over the higher levels (if necessary) to collect and combine the chunks of the representation
                usmall_t j = 1;
                for (; j < num_levels_ && !extra_bits_ranks_->is_set(p); j++) {

                    // read the chunk of the current level and add it to the partial result already computed
                    partial_res += read_bits<V, S>(levels_[j], level_cursors[j], chunk_lengths_[j]) << mult;

                    // update the position of the next chunk (in the partial result) and update cursors on the current level
                    mult += chunk_lengths_[j];
                    level_cursors[j] += chunk_lengths_[j];
                    p = accum_level_sizes_[j] + chunk_cursors[j]++;

                }

                // undo the reduction & save in result list
                list[i] = partial_res + min_values_[j - 1];

            }

            delete[] chunk_cursors;
            delete[] level_cursors;

            return list;

        }

        /*
         * Extract 'num' integers starting at index 'begin'.
         */
        V* extract_interval(ularge_t begin, ularge_t num) const {

            if ((begin + num) > length_ || length_ == 0) return 0;

            // initialise cursors to the first bit and chunk of each level
            ularge_t* level_cursors = new ularge_t[num_levels_];
            ularge_t* chunk_cursors = new ularge_t[num_levels_];

            level_cursors[0] = begin * chunk_lengths_[0];
            chunk_cursors[0] = begin;

            for (ularge_t i = begin, j = 1, pos = begin; j < num_levels_; j++) {

                // if there is no higher level for the current index, move to ("behind") the right-most predecessor
                // that has such a level no find the initial cursor positions on that level
                // + bitget(...) compensates for the additional 1-bit when there is no higher level for the current index
                i = (i + extra_bits_ranks_->is_set(pos)) - (extra_bits_ranks_->rank(accum_level_sizes_[j - 1] + i) - level_ranks_[j - 1]);

                level_cursors[j] = i * chunk_lengths_[j];
                chunk_cursors[j] = i;

                pos = accum_level_sizes_[j] + i;

            }

            // decompress the integers one by one (while using and updating the cursors)
            V* list = new V[num];
            for (ularge_t i = begin, pos = 0; i < begin + num; i++, pos++) {

                // read first chunk of representation
                V partial_res = read_bits<V, S>(levels_[0], level_cursors[0], chunk_lengths_[0]);

                ularge_t mult = chunk_lengths_[0];
                level_cursors[0] += chunk_lengths_[0];
                ularge_t p = chunk_cursors[0]++;

                // iterate over the higher levels (if necessary) to collect and combine the chunks of the representation
                usmall_t j = 1;
                for (; j < num_levels_ && !extra_bits_ranks_->is_set(p); j++) {

                    // read the chunk of the current level and add it to the partial result already computed
                    partial_res += read_bits<V, S>(levels_[j], level_cursors[j], chunk_lengths_[j]) << mult;

                    // update the position of the next chunk (in the partial result) and update cursors on the current level
                    mult += chunk_lengths_[j];
                    level_cursors[j] += chunk_lengths_[j];
                    p = accum_level_sizes_[j] + chunk_cursors[j]++;

                }

                // undo the reduction & save in result list
                list[pos] = partial_res + min_values_[j - 1];

            }

            delete[] chunk_cursors;
            delete[] level_cursors;

            return list;

        }

        /*
         * Determine the memory consumption (in bytes).
         */
        ularge_t size_in_bytes() const {

            ularge_t levels_content = 0;
            for (usmall_t i = 0; i < num_levels_; i++) {
                levels_content += (chunk_lengths_[i] * (accum_level_sizes_[i + 1] - accum_level_sizes_[i]) + (sizeof(S) * 8) - 1) / (sizeof(S) * 8);
            }

            return sizeof(SplitDacs) // SplitDacs itself
                + sizeof(S*) * num_levels_ // levels_
                + sizeof(ularge_t) * (num_levels_ + 1) // accum_level_sizes_
                + (extra_bits_ranks_ ? extra_bits_ranks_->size_in_bytes() : 0) // extra_bits_ranks_
                + sizeof(ularge_t) * num_levels_ // level_ranks_
                + sizeof(usmall_t) * num_levels_ // chunk_lengths_
                + sizeof(V) * num_levels_ // min_values_
                + sizeof(S) * levels_content; // actual content of levels_

        }

        /*
         * Return the length of the sequence.
         */
        ularge_t length() const {
            return length_;
        }

        /*
         * Restructure the SplitDacs instance to use the optimal chunk lengths and number of levels.
         */
        void optimise() {

            if (length_ == 0) return;

            Chunks chunks = (MAX_LEVELS == 0) ?
                            optimise_level_chunks<SplitDacs<V, S, MAX_LEVELS, FACTOR_RANK>, V, FACTOR_RANK>(*this, length_)
                          : optimise_level_chunks_limited<SplitDacs<V, S, MAX_LEVELS, FACTOR_RANK>, V, FACTOR_RANK>(*this, length_, MAX_LEVELS);
            usmall_t* new_chunk_lengths = chunks.chunk_lengths;
            usmall_t new_num_levels = chunks.num_levels;

            // check whether the optimisation will actually change the structure
            bool eq = (new_num_levels == num_levels_);
            for (usmall_t i = 0; eq && i < new_num_levels; i++) {
                eq = (new_chunk_lengths[i] == chunk_lengths_[i]);
            }
            if (eq) {

                delete[] new_chunk_lengths;
                return;

            }


            // for each level, determine the minimum value that makes a chunk on the level necessary
            V* new_min_values = new V[new_num_levels];
            ularge_t* new_accum_level_sizes = new ularge_t[new_num_levels + 1];
            V value = 0;
            V mult_val = 1;
            for (usmall_t i = 0; i < new_num_levels; i++) {

                mult_val *= static_cast<V>(1) << new_chunk_lengths[i];
                new_min_values[i] = value;
                value += mult_val;

                new_accum_level_sizes[i] = 0;

            }
            new_accum_level_sizes[new_num_levels] = 0;

            // for each level, determine how many chunks will be required
            for (ularge_t i = 0; i < length_; i++) {

                value = (*this)[i];
                for (usmall_t j = 0; j < new_num_levels; j++) {
                    new_accum_level_sizes[j + 1] += (value >= new_min_values[j]);
                }

            }

            // determine the total number of bits required for the chunks on all levels,
            // initialise the levels and accumulate the level sizes
            S** new_levels = new S*[new_num_levels];
            for (usmall_t i = 0; i < new_num_levels; i++) {

                new_levels[i] = new S[(new_chunk_lengths[i] * new_accum_level_sizes[i + 1] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                new_accum_level_sizes[i + 1] += new_accum_level_sizes[i];

            }

            // set up cursors into the levels for filling tem
            ularge_t* level_cursors = new ularge_t[new_num_levels]; // for each level, position (in bits) where to write the next chunk
            ularge_t* chunk_cursors = new ularge_t[new_num_levels]; // for each level, position (in chunks) of the next chunk to be written

            for (usmall_t j = 0; j < new_num_levels; j++) {

                level_cursors[j] = 0;
                chunk_cursors[j] = 0;

            }

            // initialise bit vector for Vbyte extra bits
            ularge_t len_extra_bits = new_accum_level_sizes[new_num_levels - 1] + 1; // + 1 for first chunk on last level
            S* extra_bits = new S[(len_extra_bits + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
            for (ularge_t i = 0; i < (len_extra_bits + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); i++) {
                extra_bits[i] = 0;
            }

            // number of possibles values for a chunk on the different levels (used to write the chunks below)
            V* base = new V[new_num_levels];
            for (usmall_t i = 0; i < new_num_levels; i++) {
                base[i] = static_cast<V>(1) << new_chunk_lengths[i];
            }

            // fill the levels and extra bits
            for (ularge_t i = 0; i < length_; i++) {

                value = (*this)[i];

                for (ssmall_t j = new_num_levels - 1; j >= 0; j--) {

                    if (value >= new_min_values[j]) { // find the highest required level for the current value

                        V new_value = value - new_min_values[j]; // reduce the value

                        for (usmall_t k = 0; k <= j; k++) { // write the chunks in the respective levels

                            write_bits<V, S>(new_levels[k], level_cursors[k], new_chunk_lengths[k], new_value % base[k]);
                            level_cursors[k] += new_chunk_lengths[k];
                            chunk_cursors[k]++;

                            new_value /= base[k];

                        }

                        // set the extra bit for the highest required level (unless it is the maximum level)
                        if (j < new_num_levels - 1) {
                            set_bit<S>(extra_bits, new_accum_level_sizes[j] + chunk_cursors[j] - 1);
                        }

                        j = 0; // works as break only if j is a signed integer

                    }

                }

            }

            set_bit<S>(extra_bits, len_extra_bits - 1); // write the extra bit for the first chunk on the last level

            delete[] base;
            delete[] chunk_cursors;
            delete[] level_cursors;

            // remove old contents
            for (usmall_t i = 0; i < num_levels_; i++) {
                delete[] levels_[i];
            }
            delete[] levels_;

            delete[] accum_level_sizes_;
            delete extra_bits_ranks_;
            delete[] level_ranks_;
            delete[] chunk_lengths_;
            delete[] min_values_;

            // put new contents in place
            num_levels_ = new_num_levels;

            levels_ = new_levels;
            accum_level_sizes_ = new_accum_level_sizes;
            chunk_lengths_ = new_chunk_lengths;
            min_values_ = new_min_values;

            // set up rank data structures (receives ownership of extra_bits)
            extra_bits_ranks_ = new BitRank<S, FACTOR_RANK>(extra_bits, len_extra_bits, true);

            level_ranks_ = new ularge_t[num_levels_];
            level_ranks_[0] = 0;
            for (usmall_t j = 1; j < num_levels_; j++) {
                level_ranks_[j] = extra_bits_ranks_->rank(accum_level_sizes_[j] - 1);
            }

        }

    private:
        ularge_t length_; // number of represented integers

        S** levels_; // concatenation of the chunks on all levels
        usmall_t num_levels_; // number of levels / streams
        ularge_t* accum_level_sizes_; // sum of the sizes of the levels (in chunks) before each level

        BitRank<S, FACTOR_RANK>* extra_bits_ranks_; // rank data structure on the Vbyte extra bits of the chunks
        ularge_t* level_ranks_; // number of 1-bits in the extra bits before each level (i.e. number of integers not reaching the level)

        usmall_t* chunk_lengths_; // chunk length (in bits) for each level
        V* min_values_; // for each level, minimum value that leads to a chunk on that level



    public:
        /*
         * Forward, read-only iterator to access subsequent elements.
         *
         * Provides efficient access to multiple elements without the need to store
         * all of them in decompressed from (unlike extract(Interval)).
         */
        class Reader {

        public:
            explicit Reader(const SplitDacs* d, ularge_t begin = 0) {

                dacs_ = d;

                pos_ = begin;

                level_cursors_ = new ularge_t[dacs_->num_levels_];
                chunk_cursors_ = new ularge_t[dacs_->num_levels_];

                if (begin == 0) {

                    for (usmall_t i = 0; i < dacs_->num_levels_; i++) {
                        level_cursors_[i] = chunk_cursors_[i] = 0;
                    }

                } else {

                    level_cursors_[0] = begin * dacs_->chunk_lengths_[0];
                    chunk_cursors_[0] = begin;
                    ularge_t p = begin;
                    ularge_t i = begin;

                    for (usmall_t j = 1; j < dacs_->num_levels_; j++) {

                        // if there is no higher level for the current index, move to ("behind") the right-most predecessor
                        // that has such a level no find the initial cursor positions on that level
                        // + bitget(...) compensates for the additional 1-bit when there is no higher level for the current index
                        i = (i + dacs_->extra_bits_ranks_->is_set(p)) - (dacs_->extra_bits_ranks_->rank(dacs_->accum_level_sizes_[j - 1] + i) - dacs_->level_ranks_[j - 1]);

                        level_cursors_[j] = i * dacs_->chunk_lengths_[j];
                        chunk_cursors_[j] = i;

                        p = dacs_->accum_level_sizes_[j] + i;

                    }

                }

            }

            ~Reader() {

                delete[] level_cursors_;
                delete[] chunk_cursors_;

            }

            Reader(const Reader& other) { // copy constructor

                dacs_ = other.dacs_;
                pos_ = other.pos_;

                level_cursors_ = new ularge_t[dacs_->num_levels_];
                std::copy(other.level_cursors_, other.level_cursors_ + dacs_->num_levels_, level_cursors_);

                chunk_cursors_ = new ularge_t[dacs_->num_levels_];
                std::copy(other.chunk_cursors_, other.chunk_cursors_ + dacs_->num_levels_, chunk_cursors_);

            }

            Reader(Reader&& other) noexcept { // move constructor

                dacs_ = other.dacs_;
                pos_ = other.pos_;

                level_cursors_ = other.level_cursors_; other.level_cursors_ = nullptr;
                chunk_cursors_ = other.chunk_cursors_; other.chunk_cursors_ = nullptr;

            }

            Reader& operator=(const Reader& other) { // copy assignment operator

                // check for self-assignment
                if (&other == this) {
                    return *this;
                }

                // release old resources
                delete[] level_cursors_;
                delete[] chunk_cursors_;

                // copy new resources
                dacs_ = other.dacs_;
                pos_ = other.pos_;

                level_cursors_ = new ularge_t[dacs_->num_levels_];
                std::copy(other.level_cursors_, other.level_cursors_ + dacs_->num_levels_, level_cursors_);

                chunk_cursors_ = new ularge_t[dacs_->num_levels_];
                std::copy(other.chunk_cursors_, other.chunk_cursors_ + dacs_->num_levels_, chunk_cursors_);

                return *this;

            }

            Reader& operator=(Reader&& other) noexcept { // move assignment operator

                // check for self-assignment
                if (&other == this) {
                    return *this;
                }

                // release old resources
                delete[] level_cursors_;
                delete[] chunk_cursors_;

                // copy / transfer new resources
                dacs_ = other.dacs_;
                pos_ = other.pos_;

                level_cursors_ = other.level_cursors_; other.level_cursors_ = nullptr;
                chunk_cursors_ = other.chunk_cursors_; other.chunk_cursors_ = nullptr;

                return *this;

            }

            Reader* clone() const { // deep-copy clone method
                return new Reader(*this);
            }

            /*
             * Check whether there is another element pointed to by the reader
             * (i.e. it has not reached the end of the sequence).
             */
            bool has_next() const {
                return pos_ < dacs_->length_;
            }

            /*
             * Provide the next element of the sequence and move the reader one position forward.
             * Does not check whether there is a next element (see has_next()).
             */
            V next() {

                usmall_t j;

                // read first chunk of representation
                V partial_res = read_bits<V, S>(dacs_->levels_[0], level_cursors_[0], dacs_->chunk_lengths_[0]);

                ularge_t mult = dacs_->chunk_lengths_[0];
                level_cursors_[0] += dacs_->chunk_lengths_[0];
                ularge_t p = chunk_cursors_[0]++;

                // iterate over the higher levels (if necessary) to collect and combine the chunks of the representation
                for (j = 1; j < dacs_->num_levels_ && !dacs_->extra_bits_ranks_->is_set(p); j++) {

                    // read the chunk of the current level and add it to the partial result already computed
                    partial_res += read_bits<V, S>(dacs_->levels_[j], level_cursors_[j], dacs_->chunk_lengths_[j]) << mult;

                    // update the position of the next chunk (in the partial result) and update cursors on the current level
                    mult += dacs_->chunk_lengths_[j];
                    level_cursors_[j] += dacs_->chunk_lengths_[j];
                    p = dacs_->accum_level_sizes_[j] + chunk_cursors_[j]++;

                }

                pos_++;

                // undo the reduction & save in result list
                return partial_res + dacs_->min_values_[j - 1];

            }

            /*
             * Reposition the reader on the same Dacs instance.
             */
            void reset(ularge_t begin = 0) {

                pos_ = begin;

                if (begin == 0) {

                    for (usmall_t i = 0; i < dacs_->num_levels_; i++) {
                        level_cursors_[i] = chunk_cursors_[i] = 0;
                    }

                } else {

                    level_cursors_[0] = begin * dacs_->chunk_lengths_[0];
                    chunk_cursors_[0] = begin;
                    ularge_t p = begin;
                    ularge_t i = begin;

                    for (usmall_t j = 1; j < dacs_->num_levels_; j++) {

                        // if there is no higher level for the current index, move to ("behind") the right-most predecessor
                        // that has such a level no find the initial cursor positions on that level
                        // + bitget(...) compensates for the additional 1-bit when there is no higher level for the current index
                        i = (i + dacs_->extra_bits_ranks_->is_set(p)) - (dacs_->extra_bits_ranks_->rank(dacs_->accum_level_sizes_[j - 1] + i) - dacs_->level_ranks_[j - 1]);

                        level_cursors_[j] = i * dacs_->chunk_lengths_[j];
                        chunk_cursors_[j] = i;

                        p = dacs_->accum_level_sizes_[j] + i;

                    }

                }

            }

        private:
            const SplitDacs* dacs_; // the SplitDacs instance over which the reader iterates
            ularge_t pos_; // current position of the reader in the Dacs instance
            ularge_t* level_cursors_; // level cursors into the Dacs instance for accessing the elements
            ularge_t* chunk_cursors_; // chunk cursors into the Dacs instance for accessing the elements

        };

        /*
         * Create another reader for the same SplitDacs instance starting at provided position.
         */
        Reader reader(ularge_t begin = 0) const {
            return Reader(this, begin);
        }

    };

}

#endif //DACS_SPLITDACS_HPP
