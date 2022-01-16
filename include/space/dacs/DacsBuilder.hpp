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

#ifndef DACS_DACSBUILDER_HPP
#define DACS_DACSBUILDER_HPP

#include "Utility.hpp"

namespace GeFaST {

    template<typename V, typename S, usmall_t MAX_LEVELS = 0, usmall_t FACTOR_RANK = 20>
    class Dacs;


    /*
     * Builder for a Dacs instance.
     *
     * The builder is initialised via a distribution of the elements to be added later on
     * in order to set up the internals correctly.
     * The elements are then added (appended) one by one to the sequence.
     * Once all elements are added, the construction process is completed by calling
     * finalise() once. This returns the finished Dacs.
     *
     * After finalising a construction, the builder has to be reset before it can be reused.
     */
    template<typename V, typename S, usmall_t MAX_LEVELS = 0, usmall_t FACTOR_RANK = 20>
    class DacsBuilder {

    public:
        /*
         * Initialise the builder via the distribution of the elements
         * (provided as a collection of pairs (value, number of occurrences).
         */
        explicit DacsBuilder(const std::map<ularge_t, ularge_t>& counts) {
            init(counts);
        }

        ~DacsBuilder() {
            clear();
        }


        DacsBuilder(const DacsBuilder& other) { // copy constructor

            length_ = other.length_;
            total_len_levels_ = other.total_len_levels_;
            num_levels_ = other.num_levels_;
            len_extra_bits_ = other.len_extra_bits_;
            finalised_ = other.finalised_;

            levels_ = new S[(total_len_levels_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
            std::copy(other.levels_, other.levels_ + (total_len_levels_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), levels_);

            chunk_level_index_ = new ularge_t[num_levels_ + 1];
            std::copy(other.chunk_level_index_, other.chunk_level_index_ + num_levels_ + 1, chunk_level_index_);

            bit_level_index_ = new ularge_t[num_levels_];
            std::copy(other.bit_level_index_, other.bit_level_index_ + num_levels_, bit_level_index_);

            chunk_lengths_ = new usmall_t[num_levels_];
            std::copy(other.chunk_lengths_, other.chunk_lengths_ + num_levels_, chunk_lengths_);

            min_values_ = new V[num_levels_];
            std::copy(other.min_values_, other.min_values_ + num_levels_, min_values_);

            level_cursors_ = new ularge_t[num_levels_];
            std::copy(other.level_cursors_, other.level_cursors_ + num_levels_, level_cursors_);

            chunk_cursors_ = new ularge_t[num_levels_];
            std::copy(other.chunk_cursors_, other.chunk_cursors_ + num_levels_, chunk_cursors_);

            extra_bits_ = new S[(len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
            std::copy(other.extra_bits_, other.extra_bits_ + (len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), extra_bits_);

            base_ = new V[num_levels_];
            std::copy(other.base_, other.base_ + num_levels_, base_);

        }

        DacsBuilder(DacsBuilder&& other) noexcept { // move constructor

            length_ = other.length_;
            total_len_levels_ = other.total_len_levels_;
            num_levels_ = other.num_levels_;
            len_extra_bits_ = other.len_extra_bits_;
            finalised_ = other.finalised_;

            levels_ = other.levels_; other.levels_ = nullptr;
            chunk_level_index_ = other.chunk_level_index_; other.chunk_level_index_ = nullptr;
            bit_level_index_ = other.bit_level_index_; other.bit_level_index_ = nullptr;
            chunk_lengths_ = other.chunk_lengths_; other.chunk_lengths_ = nullptr;
            min_values_ = other.min_values_; other.min_values_ = nullptr;
            level_cursors_ = other.level_cursors_; other.level_cursors_ = nullptr;
            chunk_cursors_ = other.chunk_cursors_; other.chunk_cursors_ = nullptr;
            extra_bits_ = other.extra_bits_; other.extra_bits_ = nullptr;
            base_ = other.base_; other.base_ = nullptr;

        }

        DacsBuilder& operator=(const DacsBuilder& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] levels_;
            delete[] chunk_level_index_;
            delete[] bit_level_index_;
            delete[] chunk_lengths_;
            delete[] min_values_;
            delete[] level_cursors_;
            delete[] chunk_cursors_;
            delete[] extra_bits_;
            delete[] base_;

            // copy new resources
            length_ = other.length_;
            total_len_levels_ = other.total_len_levels_;
            num_levels_ = other.num_levels_;
            len_extra_bits_ = other.len_extra_bits_;
            finalised_ = other.finalised_;

            levels_ = new S[(total_len_levels_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
            std::copy(other.levels_, other.levels_ + (total_len_levels_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), levels_);

            chunk_level_index_ = new ularge_t[num_levels_ + 1];
            std::copy(other.chunk_level_index_, other.chunk_level_index_ + num_levels_ + 1, chunk_level_index_);

            bit_level_index_ = new ularge_t[num_levels_];
            std::copy(other.bit_level_index_, other.bit_level_index_ + num_levels_, bit_level_index_);

            chunk_lengths_ = new usmall_t[num_levels_];
            std::copy(other.chunk_lengths_, other.chunk_lengths_ + num_levels_, chunk_lengths_);

            min_values_ = new V[num_levels_];
            std::copy(other.min_values_, other.min_values_ + num_levels_, min_values_);

            level_cursors_ = new ularge_t[num_levels_];
            std::copy(other.level_cursors_, other.level_cursors_ + num_levels_, level_cursors_);

            chunk_cursors_ = new ularge_t[num_levels_];
            std::copy(other.chunk_cursors_, other.chunk_cursors_ + num_levels_, chunk_cursors_);

            extra_bits_ = new S[(len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
            std::copy(other.extra_bits_, other.extra_bits_ + (len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), extra_bits_);

            base_ = new V[num_levels_];
            std::copy(other.base_, other.base_ + num_levels_, base_);

            return *this;

        }

        DacsBuilder& operator=(DacsBuilder&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] levels_;
            delete[] chunk_level_index_;
            delete[] bit_level_index_;
            delete[] chunk_lengths_;
            delete[] min_values_;
            delete[] level_cursors_;
            delete[] chunk_cursors_;
            delete[] extra_bits_;
            delete[] base_;

            // copy / transfer new resources
            length_ = other.length_;
            total_len_levels_ = other.total_len_levels_;
            num_levels_ = other.num_levels_;
            len_extra_bits_ = other.len_extra_bits_;
            finalised_ = other.finalised_;

            levels_ = other.levels_; other.levels_ = nullptr;
            chunk_level_index_ = other.chunk_level_index_; other.chunk_level_index_ = nullptr;
            bit_level_index_ = other.bit_level_index_; other.bit_level_index_ = nullptr;
            chunk_lengths_ = other.chunk_lengths_; other.chunk_lengths_ = nullptr;
            min_values_ = other.min_values_; other.min_values_ = nullptr;
            level_cursors_ = other.level_cursors_; other.level_cursors_ = nullptr;
            chunk_cursors_ = other.chunk_cursors_; other.chunk_cursors_ = nullptr;
            extra_bits_ = other.extra_bits_; other.extra_bits_ = nullptr;
            base_ = other.base_; other.base_ = nullptr;

            return *this;

        }

        DacsBuilder* clone() const { // deep-copy clone method
            return new DacsBuilder(*this);
        }

        /*
         * Add the value to the end of the Dacs representation of the sequence.
         */
        void append(V val) {

            for (ssmall_t j = num_levels_ - 1; j >= 0; j--) {

                if (val >= min_values_[j]) { // find the highest required level for the current value

                    V new_value = val - min_values_[j]; // reduce the value

                    for (usmall_t k = 0; k <= j; k++) { // write the chunks in the respective levels

                        write_bits<V, S>(levels_, level_cursors_[k], chunk_lengths_[k], new_value % base_[k]);
                        level_cursors_[k] += chunk_lengths_[k];
                        chunk_cursors_[k]++;

                        new_value /= base_[k];

                    }

                    // set the extra bit for the highest required level (unless it is the maximum level)
                    if (j < num_levels_ - 1) {
                        set_bit<S>(extra_bits_, chunk_cursors_[j] - 1);
                    }

                    j = 0; // works as break only if j is a signed integer

                }

            }

        }

        /*
         * Finalise the Dacs representation (after all elements of the sequence have been added to it)
         * and construct the actual Dacs instance (repeated calls return the null pointer).
         */
        Dacs<V, S, MAX_LEVELS, FACTOR_RANK>* finalise(bool opt = true) {

            if (finalised_) {
                return nullptr;
            } else {

                set_bit<S>(extra_bits_, len_extra_bits_ - 1); // write the extra bit for the first chunk on the last level

                delete[] base_; base_ = nullptr;
                delete[] chunk_cursors_; chunk_cursors_ = nullptr;
                delete[] level_cursors_; level_cursors_ = nullptr;

                // set up rank data structures (receives ownership of extra_bits_)
                BitRank<S, FACTOR_RANK>* extra_bits_ranks = new BitRank<S, FACTOR_RANK>(extra_bits_, len_extra_bits_, true);
                extra_bits_ = nullptr;

                ularge_t* level_ranks = new ularge_t[num_levels_];
                level_ranks[0] = 0;
                for (usmall_t j = 1; j < num_levels_; j++) {
                    level_ranks[j] = extra_bits_ranks->rank(chunk_level_index_[j] - 1);
                }

                auto dacs = new Dacs<V, S, MAX_LEVELS, FACTOR_RANK>(
                        length_,
                        levels_,
                        total_len_levels_,
                        num_levels_,
                        chunk_level_index_,
                        bit_level_index_,
                        extra_bits_ranks, // takes ownership
                        level_ranks, // takes ownership
                        chunk_lengths_,
                        min_values_
                );

                levels_ = nullptr;
                chunk_level_index_ = nullptr;
                bit_level_index_ = nullptr;
                chunk_lengths_ = nullptr;
                min_values_ = nullptr;
                extra_bits_ranks = nullptr;
                level_ranks = nullptr;

                finalised_ = true;

                return dacs;

            }

        }

        /*
         * Delete the current contents of the builder and prepare it for building another Dacs instance.
         */
        void reset(const std::map<ularge_t, ularge_t>& counts) {

            clear();
            init(counts);

        }

        /*
         * Determine the memory consumption (in bytes).
         */
        ularge_t size_in_bytes() const {
            return sizeof(DacsBuilder) // DacsBuilder itself
                + sizeof(S) * ((total_len_levels_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)) // levels_
                + sizeof(ularge_t) * (num_levels_ + 1) // chunk_level_index_
                + sizeof(ularge_t) * num_levels_ // bit_level_index_
                + sizeof(usmall_t) * num_levels_ // chunk_lengths_
                + sizeof(V) * num_levels_ // min_values_
                + sizeof(ularge_t) * num_levels_ // level_cursors_
                + sizeof(ularge_t) * num_levels_ // chunk_cursors_
                + sizeof(S) * ((len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)) // extra_bits_
                + sizeof(V) * num_levels_; // base_
        }

    private:
        // members of Dacs instance to be built (see documentation in Dacs)
        ularge_t length_;
        S* levels_;
        ularge_t total_len_levels_;
        usmall_t num_levels_;
        ularge_t* chunk_level_index_;
        ularge_t* bit_level_index_;
        usmall_t* chunk_lengths_;
        V* min_values_;

        // temporary helper members needed during the construction
        ularge_t* level_cursors_; // position (in bits) where to write the next chunk (per level)
        ularge_t* chunk_cursors_; // position (in chunks) of the next chunk to be written (per level)
        ularge_t len_extra_bits_; // length of bit vector for Vbyte extra bits
        S* extra_bits_; // bit vector for Vbyte extra bits
        V* base_; // number of possibles values on the different levels (based on chunk length)

        bool finalised_; // flag indicating that the builder is already been used to construct the Dacs instance

        /*
         * Prepare the builder for adding the elements of a sequence following a distribution
         * provided as a collection of pairs (value, number of occurrences).
         */
        void init(const std::map<ularge_t, ularge_t>& counts) {

            length_ = 0;
            for (auto& c : counts) {
                length_ += c.second;
            }

            // determine number of levels and their chunk lengths
            Chunks chunks = (MAX_LEVELS == 0) ?
                            optimise_level_chunks<V, FACTOR_RANK>(counts)
                          : optimise_level_chunks_limited<V, FACTOR_RANK>(counts, MAX_LEVELS);
            chunk_lengths_ = chunks.chunk_lengths;
            num_levels_ = chunks.num_levels;

            // for each level, determine the minimum value that makes a chunk on the level necessary
            min_values_ = new V[num_levels_];
            ularge_t* level_sizes = new ularge_t[num_levels_];
            V value = 0;
            V mult_val = 1;
            for (usmall_t i = 0; i < num_levels_; i++) {

                mult_val *= static_cast<V>(1) << chunk_lengths_[i];
                min_values_[i] = value;
                value += mult_val;

                level_sizes[i] = 0;

            }

            // for each level, determine how many chunks will be required
            for (auto& c : counts) {

                value = c.first;
                for (usmall_t j = 0; j < num_levels_; j++) {
                    level_sizes[j] += c.second * (value >= min_values_[j]);
                }

            }

            // determine the total number of bits required for the chunks on all levels
            total_len_levels_ = 0;
            for (usmall_t i = 0; i < num_levels_; i++) {
                total_len_levels_ += chunk_lengths_[i] * level_sizes[i];
            }

            // determine the positions of the first chunk / bit per level
            chunk_level_index_ = new ularge_t[num_levels_ + 1];
            bit_level_index_ = new ularge_t[num_levels_];
            level_cursors_ = new ularge_t[num_levels_]; // for each level, position (in bits) where to write the next chunk
            chunk_cursors_ = new ularge_t[num_levels_]; // for each level, position (in chunks) of the next chunk to be written

            bit_level_index_[0] = level_cursors_[0] = 0;
            chunk_level_index_[0] = chunk_cursors_[0] = 0;
            for (usmall_t j = 1; j < num_levels_; j++) {

                bit_level_index_[j] = level_cursors_[j] = bit_level_index_[j - 1] + level_sizes[j - 1] * chunk_lengths_[j - 1];
                chunk_level_index_[j] = chunk_cursors_[j] = chunk_level_index_[j - 1] + level_sizes[j - 1];

            }
            chunk_level_index_[num_levels_] = chunk_level_index_[num_levels_ - 1] + level_cursors_[num_levels_ - 1];

            delete[] level_sizes;

            // initialise bit vector for Vbyte extra bits
            len_extra_bits_ = chunk_level_index_[num_levels_ - 1] + 1; // + 1 for first chunk on last level
            extra_bits_ = new S[(len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
            for (ularge_t i = 0; i < (len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); i++) {
                extra_bits_[i] = 0;
            }

            // number of possibles values for a chunk on the different levels (used to write the chunks below)
            base_ = new V[num_levels_];
            for (usmall_t i = 0; i < num_levels_; i++) {
                base_[i] = static_cast<V>(1) << chunk_lengths_[i];
            }

            levels_ = new S[(total_len_levels_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];

            finalised_ = false;

        }

        /*
         * Delete the current contents of the builder (when destructing or resetting it).
         */
        void clear() {

            delete[] levels_;
            delete[] chunk_level_index_;
            delete[] bit_level_index_;
            delete[] chunk_lengths_;
            delete[] min_values_;

            delete[] level_cursors_;
            delete[] chunk_cursors_;
            delete[] extra_bits_;
            delete[] base_;

        }

    };

}

#endif //DACS_DACSBUILDER_HPP
