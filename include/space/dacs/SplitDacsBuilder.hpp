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

#ifndef DACS_SPLITDACSBUILDER_HPP
#define DACS_SPLITDACSBUILDER_HPP

#include "Utility.hpp"

namespace GeFaST {

    template<typename V, typename S, usmall_t MAX_LEVELS = 0, usmall_t FACTOR_RANK = 20>
    class SplitDacs;


    /*
     * Builder for a SplitDacs instance.
     *
     * The builder is initialised via a distribution of the elements to be added later on
     * in order to set up the internals correctly.
     * The elements are then added (appended) one by one to the sequence.
     * Once all elements are added, the construction process is completed by calling
     * finalise() once. This returns the finished SplitDacs.
     *
     * After finalising a construction, the builder has to be reset before it can be reused.
     */
    template<typename V, typename S, usmall_t MAX_LEVELS = 0, usmall_t FACTOR_RANK = 20>
    class SplitDacsDistributionBuilder {

    public:
        /*
         * Initialise the builder via the distribution of the elements
         * (provided as a collection of pairs (value, number of occurrences).
         */
        explicit SplitDacsDistributionBuilder(const std::map<ularge_t, ularge_t>& counts) {
            init(counts);
        }

        /*
         * Initialise the builder based on the given configuration.
         */
        explicit SplitDacsDistributionBuilder(SplitDacsParameters<V, S, MAX_LEVELS, FACTOR_RANK>& conf) {

            length_ = conf.length;
            num_levels_ = conf.num_levels;

            chunk_lengths_ = new usmall_t[num_levels_];
            std::copy(conf.chunk_lengths, conf.chunk_lengths + num_levels_, chunk_lengths_);

            min_values_ = new V[num_levels_];
            std::copy(conf.min_values, conf.min_values + num_levels_, min_values_);

            accum_level_sizes_ = new ularge_t[num_levels_ + 1];
            std::copy(conf.accum_level_sizes, conf.accum_level_sizes + num_levels_ + 1, accum_level_sizes_);

            // determine the total number of bits required for the chunks on all levels,
            // initialise the levels and accumulate the level sizes
            levels_ = new S*[num_levels_];
            for (usmall_t i = 0; i < num_levels_; i++) {

                levels_[i] = new S[(chunk_lengths_[i] * accum_level_sizes_[i + 1] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                accum_level_sizes_[i + 1] += accum_level_sizes_[i];

            }

            // set up cursors into the levels for filling tem
            level_cursors_ = new ularge_t[num_levels_]; // for each level, position (in bits) where to write the next chunk
            chunk_cursors_ = new ularge_t[num_levels_]; // for each level, position (in chunks) of the next chunk to be written

            for (usmall_t j = 0; j < num_levels_; j++) {

                level_cursors_[j] = 0;
                chunk_cursors_[j] = 0;

            }

            // initialise bit vector for Vbyte extra bits
            len_extra_bits_ = accum_level_sizes_[num_levels_ - 1] + 1; // + 1 for first chunk on last level
            extra_bits_ = new S[(len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
            for (ularge_t i = 0; i < (len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); i++) {
                extra_bits_[i] = 0;
            }

            // number of possibles values for a chunk on the different levels (used to write the chunks below)
            base_ = new V[num_levels_];
            for (usmall_t i = 0; i < num_levels_; i++) {
                base_[i] = static_cast<V>(1) << chunk_lengths_[i];
            }

            finalised_ = false;

        }

        ~SplitDacsDistributionBuilder() {
            clear();
        }

        SplitDacsDistributionBuilder(const SplitDacsDistributionBuilder& other) { // copy constructor

            length_ = other.length_;
            num_levels_ = other.num_levels_;
            len_extra_bits_ = other.len_extra_bits_;
            finalised_ = other.finalised_;

            accum_level_sizes_ = new ularge_t[num_levels_ + 1];
            std::copy(other.accum_level_sizes_, other.accum_level_sizes_ + num_levels_ + 1, accum_level_sizes_);

            chunk_lengths_ = new usmall_t[num_levels_];
            std::copy(other.chunk_lengths_, other.chunk_lengths_ + num_levels_, chunk_lengths_);

            levels_ = new S*[num_levels_];
            for (usmall_t i = 0; i < num_levels_; i++) {

                ularge_t len = (chunk_lengths_[i] * (accum_level_sizes_[i + 1] - accum_level_sizes_[i]) + (sizeof(S) * 8) - 1) / (sizeof(S) * 8);
                levels_[i] = new S[len];
                std::copy(other.levels_[i], other.levels_[i] + len, levels_[i]);

            }

            min_values_ = new V[num_levels_];
            std::copy(other.min_values_, other.min_values_ + num_levels_, min_values_);

            level_cursors_ = new ularge_t[num_levels_];
            std::copy(other.level_cursors_, other.level_cursors_ + num_levels_, level_cursors_);

            chunk_cursors_ = new ularge_t[num_levels_];
            std::copy(other.chunk_cursors_, other.chunk_cursors_ + num_levels_, chunk_cursors_);

            ularge_t len = (len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8);
            extra_bits_ = new S[len];
            std::copy(other.extra_bits_, other.extra_bits_ + len, extra_bits_);

            base_ = new V[num_levels_];
            std::copy(other.base_, other.base_ + num_levels_, base_);

        }

        SplitDacsDistributionBuilder(SplitDacsDistributionBuilder&& other) noexcept { // move constructor

            length_ = other.length_;
            num_levels_ = other.num_levels_; other.num_levels_ = 0;
            len_extra_bits_ = other.len_extra_bits_;
            finalised_ = other.finalised_;

            accum_level_sizes_ = other.accum_level_sizes_; other.accum_level_sizes_ = nullptr;
            chunk_lengths_ = other.chunk_lengths_; other.chunk_lengths_ = nullptr;
            levels_ = other.levels_; other.levels_ = nullptr;
            min_values_ = other.min_values_; other.min_values_ = nullptr;
            level_cursors_ = other.level_cursors_; other.level_cursors_ = nullptr;
            chunk_cursors_ = other.chunk_cursors_; other.chunk_cursors_ = nullptr;
            extra_bits_ = other.extra_bits_; other.extra_bits_ = nullptr;
            base_ = other.base_; other.base_ = nullptr;

        }

        SplitDacsDistributionBuilder& operator=(const SplitDacsDistributionBuilder& other) { // copy assignment operator

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
            delete[] chunk_lengths_;
            delete[] min_values_;
            delete[] level_cursors_;
            delete[] chunk_cursors_;
            delete[] extra_bits_;
            delete[] base_;

            // copy new resources
            length_ = other.length_;
            num_levels_ = other.num_levels_;
            len_extra_bits_ = other.len_extra_bits_;
            finalised_ = other.finalised_;

            accum_level_sizes_ = new ularge_t[num_levels_ + 1];
            std::copy(other.accum_level_sizes_, other.accum_level_sizes_ + num_levels_ + 1, accum_level_sizes_);

            chunk_lengths_ = new usmall_t[num_levels_];
            std::copy(other.chunk_lengths_, other.chunk_lengths_ + num_levels_, chunk_lengths_);

            levels_ = new S*[num_levels_];
            for (usmall_t i = 0; i < num_levels_; i++) {

                ularge_t len = (chunk_lengths_[i] * (accum_level_sizes_[i + 1] - accum_level_sizes_[i]) + (sizeof(S) * 8) - 1) / (sizeof(S) * 8);
                levels_[i] = new S[len];
                std::copy(other.levels_[i], other.levels_[i] + len, levels_[i]);

            }

            min_values_ = new V[num_levels_];
            std::copy(other.min_values_, other.min_values_ + num_levels_, min_values_);

            level_cursors_ = new ularge_t[num_levels_];
            std::copy(other.level_cursors_, other.level_cursors_ + num_levels_, level_cursors_);

            chunk_cursors_ = new ularge_t[num_levels_];
            std::copy(other.chunk_cursors_, other.chunk_cursors_ + num_levels_, chunk_cursors_);

            ularge_t len = (len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8);
            extra_bits_ = new S[len];
            std::copy(other.extra_bits_, other.extra_bits_ + len, extra_bits_);

            base_ = new V[num_levels_];
            std::copy(other.base_, other.base_ + num_levels_, base_);

            return *this;

        }

        SplitDacsDistributionBuilder& operator=(SplitDacsDistributionBuilder&& other) noexcept { // move assignment operator

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
            delete[] chunk_lengths_;
            delete[] min_values_;
            delete[] level_cursors_;
            delete[] chunk_cursors_;
            delete[] extra_bits_;
            delete[] base_;

            // copy / transfer new resources
            length_ = other.length_;
            num_levels_ = other.num_levels_; other.num_levels_ = 0;
            len_extra_bits_ = other.len_extra_bits_;
            finalised_ = other.finalised_;

            accum_level_sizes_ = other.accum_level_sizes_; other.accum_level_sizes_ = nullptr;
            chunk_lengths_ = other.chunk_lengths_; other.chunk_lengths_ = nullptr;
            levels_ = other.levels_; other.levels_ = nullptr;
            min_values_ = other.min_values_; other.min_values_ = nullptr;
            level_cursors_ = other.level_cursors_; other.level_cursors_ = nullptr;
            chunk_cursors_ = other.chunk_cursors_; other.chunk_cursors_ = nullptr;
            extra_bits_ = other.extra_bits_; other.extra_bits_ = nullptr;
            base_ = other.base_; other.base_ = nullptr;

            return *this;

        }

        SplitDacsDistributionBuilder* clone() const { // deep-copy clone method
            return new SplitDacsDistributionBuilder(*this);
        }

        /*
         * Add the value to the end of the SplitDacs representation of the sequence.
         */
        void append(V val) {

            for (ssmall_t j = num_levels_ - 1; j >= 0; j--) {

                if (val >= min_values_[j]) { // find the highest required level for the current value

                    V new_value = val - min_values_[j]; // reduce the value

                    for (usmall_t k = 0; k <= j; k++) { // write the chunks in the respective levels

                        write_bits<V, S>(levels_[k], level_cursors_[k], chunk_lengths_[k], new_value % base_[k]);
                        level_cursors_[k] += chunk_lengths_[k];
                        chunk_cursors_[k]++;

                        new_value /= base_[k];

                    }

                    // set the extra bit for the highest required level (unless it is the maximum level)
                    if (j < num_levels_ - 1) {
                        set_bit<S>(extra_bits_, accum_level_sizes_[j] + chunk_cursors_[j] - 1);
                    }

                    j = 0; // works as break only if j is a signed integer

                }

            }

        }

        /*
         * Finalise the SplitDacs representation (after all elements of the sequence have been added to it)
         * and construct the actual SplitDacs instance (repeated calls return the null pointer).
         */
        SplitDacs<V, S, MAX_LEVELS, FACTOR_RANK>* finalise() {

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
                    level_ranks[j] = extra_bits_ranks->rank(accum_level_sizes_[j] - 1);
                }

                auto dacs = new SplitDacs<V, S, MAX_LEVELS, FACTOR_RANK>(
                        length_,
                        levels_,
                        num_levels_,
                        accum_level_sizes_,
                        extra_bits_ranks, // takes ownership
                        level_ranks, // takes ownership
                        chunk_lengths_,
                        min_values_
                );

                num_levels_ = 0;
                levels_ = nullptr;
                accum_level_sizes_ = nullptr;
                chunk_lengths_ = nullptr;
                min_values_ = nullptr;
                extra_bits_ranks = nullptr;
                level_ranks = nullptr;

                finalised_ = true;

                return dacs;

            }

        }

        /*
         * Delete the current contents of the builder and prepare it for building another SplitDacs instance.
         */
        void reset(const std::map<ularge_t, ularge_t>& counts) {

            clear();
            init(counts);

        }

        /*
         * Determine memory consumption (in bytes).
         */
        ularge_t size_in_bytes() const {

            ularge_t contents = 0;
            for (usmall_t i = 0; i < num_levels_; i++) {
                contents += (chunk_lengths_[i] * (accum_level_sizes_[i + 1] - accum_level_sizes_[i]) + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); // levels_
            }

            return sizeof(SplitDacsDistributionBuilder) // SplitDacsDistributionBuilder itself
                + sizeof(S*) * num_levels_ // levels_
                + sizeof(ularge_t) * (num_levels_ + 1) // accum_level_sizes_
                + sizeof(usmall_t) * num_levels_ // chunk_lengths_
                + sizeof(V) * num_levels_ // min_values_
                + sizeof(ularge_t) * num_levels_ // level_cursors_
                + sizeof(ularge_t) * num_levels_ // chunk_cursors_
                + sizeof(S) * ((len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)) // extra_bits_
                + sizeof(V) * num_levels_ // base_
                + sizeof(S) * contents; // actual contents of levels_

        }

    private:
        // members of SplitDacs instance to manage persistently (see documentation in SplitDacs)
        ularge_t length_;
        S** levels_;
        usmall_t num_levels_;
        ularge_t* accum_level_sizes_;
        usmall_t* chunk_lengths_;
        V* min_values_;

        // temporary helper members needed during the construction
        ularge_t* level_cursors_; // position (in bits) where to write the next chunk (per level)
        ularge_t* chunk_cursors_; // position (in chunks) of the next chunk to be written (per level)
        ularge_t len_extra_bits_; // length of bit vector for Vbyte extra bits
        S* extra_bits_; // bit vector for Vbyte extra bits
        V* base_; // number of possibles values on the different levels (based on chunk lengths)

        bool finalised_; // flag indicating that the builder is already been used to construct the SplitDacs

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
            for (auto& c : counts) {

                value = c.first;
                for (usmall_t j = 0; j < num_levels_; j++) {
                    accum_level_sizes_[j + 1] += c.second * (value >= min_values_[j]);
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
            level_cursors_ = new ularge_t[num_levels_]; // for each level, position (in bits) where to write the next chunk
            chunk_cursors_ = new ularge_t[num_levels_]; // for each level, position (in chunks) of the next chunk to be written

            for (usmall_t j = 0; j < num_levels_; j++) {

                level_cursors_[j] = 0;
                chunk_cursors_[j] = 0;

            }

            // initialise bit vector for Vbyte extra bits
            len_extra_bits_ = accum_level_sizes_[num_levels_ - 1] + 1; // + 1 for first chunk on last level
            extra_bits_ = new S[(len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
            for (ularge_t i = 0; i < (len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); i++) {
                extra_bits_[i] = 0;
            }

            // number of possibles values for a chunk on the different levels (used to write the chunks below)
            base_ = new V[num_levels_];
            for (usmall_t i = 0; i < num_levels_; i++) {
                base_[i] = static_cast<V>(1) << chunk_lengths_[i];
            }

            finalised_ = false;

        }

        /*
         * Delete the current contents of the builder (when destructing or resetting it).
         */
        void clear() {

            for (usmall_t i = 0; i < num_levels_; i++) {
                delete[] levels_[i];
            }
            delete[] levels_;
            delete[] accum_level_sizes_;
            delete[] chunk_lengths_;
            delete[] min_values_;

            delete[] level_cursors_;
            delete[] chunk_cursors_;
            delete[] extra_bits_;
            delete[] base_;

        }

    };

    /*
     * Builder for a SplitDacs instance.
     *
     * The builder is initialised with the chunk length used for all levels and
     * the initial capacity of the first level.
     * The elements are then added (appended) one by one to the sequence.
     * This adds levels and enlarges existing levels as needed.
     * Once all elements are added, the construction process is completed by calling
     * finalise() once. This returns the finished SplitDacs - after the
     * number of levels and chunk lengths have been optimised (optional).
     *
     * After finalising a construction, the builder has to be reset before it can be reused.
     */
    template<typename V, typename S, usmall_t MAX_LEVELS = 0, usmall_t FACTOR_RANK = 20>
    class SplitDacsDirectBuilder {

    public:
        /*
         * Initialise the builder with the provided chunk length and an initial capacity.
         */
        SplitDacsDirectBuilder(usmall_t chunk_length, ularge_t capacity, float level_factor, float extend_factor) {
            init(chunk_length, capacity, level_factor, extend_factor);
        }

        ~SplitDacsDirectBuilder() {
            clear();
        }


        SplitDacsDirectBuilder(const SplitDacsDirectBuilder& other) { // copy constructor

            length_ = other.length_;
            num_levels_ = other.num_levels_;
            len_extra_bits_ = other.len_extra_bits_;
            chunk_length_ = other.chunk_length_;
            base_ = other.base_;
            levels_limit_ = other.levels_limit_;
            finalised_ = other.finalised_;
            level_factor_ = other.level_factor_;
            extend_factor_ = other.extend_factor_;

            min_values_ = new V[levels_limit_];
            std::copy(other.min_values_, other.min_values_ + levels_limit_, min_values_);

            level_cursors_ = new ularge_t[levels_limit_];
            std::copy(other.level_cursors_, other.level_cursors_ + levels_limit_, level_cursors_);

            chunk_cursors_ = new ularge_t[levels_limit_];
            std::copy(other.chunk_cursors_, other.chunk_cursors_ + levels_limit_, chunk_cursors_);

            level_capacities_ = new ularge_t[levels_limit_];
            std::copy(other.level_capacities_, other.level_capacities_ + levels_limit_, level_capacities_);

            levels_ = new S*[levels_limit_];
            level_extra_bits_ = new S*[levels_limit_];
            for (usmall_t i = 0; i < num_levels_; i++) {

                levels_[i] = new S[(chunk_length_ * level_capacities_[i] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                std::copy(other.levels_[i], other.levels_[i] + (chunk_length_ * level_capacities_[i] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), levels_[i]);

                level_extra_bits_[i] = new S[(level_capacities_[i] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                std::copy(other.level_extra_bits_[i], other.level_extra_bits_[i] + (level_capacities_[i] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), level_extra_bits_[i]);

            }
            for (usmall_t i = num_levels_; i < levels_limit_; i++) {
                levels_[i] = level_extra_bits_[i] = 0;
            }

        }

        SplitDacsDirectBuilder(SplitDacsDirectBuilder&& other) noexcept { // move constructor

            length_ = other.length_;
            num_levels_ = other.num_levels_;
            len_extra_bits_ = other.len_extra_bits_;
            chunk_length_ = other.chunk_length_;
            base_ = other.base_;
            levels_limit_ = other.levels_limit_;
            finalised_ = other.finalised_;
            level_factor_ = other.level_factor_;
            extend_factor_ = other.extend_factor_;

            min_values_ = other.min_values_; other.min_values_ = nullptr;
            level_cursors_ = other.level_cursors_; other.level_cursors_ = nullptr;
            chunk_cursors_ = other.chunk_cursors_; other.chunk_cursors_ = nullptr;
            level_capacities_ = other.level_capacities_; other.level_capacities_ = nullptr;
            levels_ = other.levels_; other.levels_ = nullptr;
            level_extra_bits_ = other.level_extra_bits_; other.level_extra_bits_ = nullptr;

        }

        SplitDacsDirectBuilder& operator=(const SplitDacsDirectBuilder& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            for (usmall_t i = 0; i < levels_limit_; i++) {
                delete[] levels_[i];
            }
            delete[] levels_;
            delete[] min_values_;

            delete[] level_cursors_;
            delete[] chunk_cursors_;

            delete[] level_capacities_;
            for (usmall_t i = 0; i < levels_limit_; i++) {
                delete[] level_extra_bits_[i];
            }
            delete[] level_extra_bits_;

            // copy new resources
            length_ = other.length_;
            num_levels_ = other.num_levels_;
            len_extra_bits_ = other.len_extra_bits_;
            chunk_length_ = other.chunk_length_;
            base_ = other.base_;
            levels_limit_ = other.levels_limit_;
            finalised_ = other.finalised_;
            level_factor_ = other.level_factor_;
            extend_factor_ = other.extend_factor_;

            min_values_ = new V[levels_limit_];
            std::copy(other.min_values_, other.min_values_ + levels_limit_, min_values_);

            level_cursors_ = new ularge_t[levels_limit_];
            std::copy(other.level_cursors_, other.level_cursors_ + levels_limit_, level_cursors_);

            chunk_cursors_ = new ularge_t[levels_limit_];
            std::copy(other.chunk_cursors_, other.chunk_cursors_ + levels_limit_, chunk_cursors_);

            level_capacities_ = new ularge_t[levels_limit_];
            std::copy(other.level_capacities_, other.level_capacities_ + levels_limit_, level_capacities_);

            levels_ = new S*[levels_limit_];
            level_extra_bits_ = new S*[levels_limit_];
            for (usmall_t i = 0; i < num_levels_; i++) {

                levels_[i] = new S[(chunk_length_ * level_capacities_[i] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                std::copy(other.levels_[i], other.levels_[i] + (chunk_length_ * level_capacities_[i] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), levels_[i]);

                level_extra_bits_[i] = new S[(level_capacities_[i] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                std::copy(other.level_extra_bits_[i], other.level_extra_bits_[i] + (level_capacities_[i] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), level_extra_bits_[i]);

            }
            for (usmall_t i = num_levels_; i < levels_limit_; i++) {
                levels_[i] = level_extra_bits_[i] = 0;
            }

            return *this;

        }

        SplitDacsDirectBuilder& operator=(SplitDacsDirectBuilder&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            for (usmall_t i = 0; i < levels_limit_; i++) {
                delete[] levels_[i];
            }
            delete[] levels_;
            delete[] min_values_;

            delete[] level_cursors_;
            delete[] chunk_cursors_;

            delete[] level_capacities_;
            for (usmall_t i = 0; i < levels_limit_; i++) {
                delete[] level_extra_bits_[i];
            }
            delete[] level_extra_bits_;

            // copy / transfer new resources
            length_ = other.length_;
            num_levels_ = other.num_levels_;
            len_extra_bits_ = other.len_extra_bits_;
            chunk_length_ = other.chunk_length_;
            base_ = other.base_;
            levels_limit_ = other.levels_limit_;
            finalised_ = other.finalised_;
            level_factor_ = other.level_factor_;
            extend_factor_ = other.extend_factor_;

            min_values_ = other.min_values_; other.min_values_ = nullptr;
            level_cursors_ = other.level_cursors_; other.level_cursors_ = nullptr;
            chunk_cursors_ = other.chunk_cursors_; other.chunk_cursors_ = nullptr;
            level_capacities_ = other.level_capacities_; other.level_capacities_ = nullptr;
            levels_ = other.levels_; other.levels_ = nullptr;
            level_extra_bits_ = other.level_extra_bits_; other.level_extra_bits_ = nullptr;

            return *this;

        }

        SplitDacsDirectBuilder* clone() const { // deep-copy clone method
            return new SplitDacsDirectBuilder(*this);
        }

        /*
         * Add the value to the end of the SplitDacs representation of the sequence.
         */
        void append(V val) {

            while (num_levels_ < levels_limit_ && val >= min_values_[num_levels_]) {

                // new level has half the capacity of the previous level, but at least a capacity of 1
                level_capacities_[num_levels_] = std::max(std::llroundf(level_capacities_[num_levels_ - 1] * level_factor_), 1ll);

                // create the new level
                levels_[num_levels_] = new S[(chunk_length_ * level_capacities_[num_levels_] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];

                // create extra bits for the new level
                level_extra_bits_[num_levels_] = new S[(level_capacities_[num_levels_] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                for (ularge_t i = 0; i < (level_capacities_[num_levels_] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); i++) {
                    level_extra_bits_[num_levels_][i] = 0;
                }

                num_levels_++;

            }

            for (ssmall_t j = num_levels_ - 1; j >= 0; j--) {

                if (val >= min_values_[j]) { // find the highest required level for the current value

                    V new_value = val - min_values_[j]; // reduce the value

                    for (usmall_t k = 0; k <= j; k++) { // write the chunks in the respective levels

                        if (chunk_cursors_[k] == level_capacities_[k]) { // level is full, extend first

                            // extend level itself
                            ularge_t new_level_capa = (k == 0) ? std::max(level_capacities_[k] + 1, static_cast<ularge_t>(std::llroundf(level_capacities_[k] * extend_factor_)))
                                    : std::min(level_capacities_[k - 1], std::max(level_capacities_[k] + 1, static_cast<ularge_t>(std::llroundf(level_capacities_[k] * extend_factor_))));
                            S* tmp = new S[(chunk_length_ * new_level_capa + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                            std::copy(levels_[k], levels_[k] + (chunk_length_ * level_capacities_[k] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), tmp);
                            delete[] levels_[k];
                            levels_[k] = tmp;

                            // extend extra bits of the level
                            tmp = new S[(new_level_capa + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                            for (ularge_t i = 0; i < (new_level_capa + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); i++) {
                                tmp[i] = 0;
                            }
                            std::copy(level_extra_bits_[k], level_extra_bits_[k] + (level_capacities_[k] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), tmp);
                            delete[] level_extra_bits_[k];
                            level_extra_bits_[k] = tmp;

                            level_capacities_[k] = new_level_capa;

                        }

                        len_extra_bits_++;
                        write_bits<V, S>(levels_[k], level_cursors_[k], chunk_length_, new_value % base_);
                        level_cursors_[k] += chunk_length_;
                        chunk_cursors_[k]++;

                        new_value /= base_;

                    }

                    // set the extra bit always (currently highest level might not be the final maximum level)
                    set_bit<S>(level_extra_bits_[j], chunk_cursors_[j] - 1);

                    j = 0; // works as break only if j is a signed integer

                }

            }

            length_++;

        }

        /*
         * Finalise the SplitDacs representation (after all elements of the sequence have been added to it),
         * optionally optimise the structure, and construct the actual SplitDacs instance (repeated calls return the null pointer).
         */
        SplitDacs<V, S, MAX_LEVELS, FACTOR_RANK>* finalise(bool opt = true) {

            if (finalised_) {
                return nullptr;
            } else {

                delete[] level_capacities_; level_capacities_ = nullptr;
                delete[] level_cursors_; level_cursors_ = nullptr;

                // determine accumulated level sizes
                ularge_t* accum_level_sizes = new ularge_t[num_levels_ + 1];
                accum_level_sizes[0] = 0;
                for (usmall_t j = 1; j <= num_levels_; j++) {
                    accum_level_sizes[j] = accum_level_sizes[j - 1] + chunk_cursors_[j - 1];
                }

                // combine the level-wise extra bits
                S* extra_bits = new S[(len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                for (ularge_t i = 0; i < (len_extra_bits_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); i++) {
                    extra_bits[i] = 0;
                }

                for (usmall_t i = 0; i < num_levels_; i++) {

                    ularge_t pos = accum_level_sizes[i];
                    for (ularge_t j = 0; j < chunk_cursors_[i] / (sizeof(S) * 8); j++, pos += (sizeof(S) * 8)) {
                        write_bits<S, S>(extra_bits, pos, sizeof(S) * 8, level_extra_bits_[i][j]);
                    }
                    if (chunk_cursors_[i] % (sizeof(S) * 8) != 0) {
                        write_bits<S, S>(extra_bits, pos, chunk_cursors_[i] % (sizeof(S) * 8),
                                         read_bits<S, S>(level_extra_bits_[i], pos - accum_level_sizes[i],
                                                         chunk_cursors_[i] % (sizeof(S) * 8)));
                    }

                }

                set_bit<S>(extra_bits, len_extra_bits_ - 1); // write the extra bit for the first chunk on the last level

                delete[] chunk_cursors_; chunk_cursors_ = nullptr;
                for (usmall_t i = 0; i < levels_limit_; i++) {
                    delete[] level_extra_bits_[i];
                }
                delete[] level_extra_bits_; level_extra_bits_ = nullptr;


                // set up rank data structures (receives ownership of extra_bits)
                BitRank<S, FACTOR_RANK>* extra_bits_ranks = new BitRank<S, FACTOR_RANK>(extra_bits, len_extra_bits_, true);
                extra_bits = nullptr;

                ularge_t* level_ranks = new ularge_t[num_levels_];
                level_ranks[0] = 0;
                for (usmall_t j = 1; j < num_levels_; j++) {
                    level_ranks[j] = extra_bits_ranks->rank(accum_level_sizes[j] - 1);
                }


                // expand the chunk length into an array handed over to the constructor
                usmall_t* chunk_lengths = new usmall_t[num_levels_];
                for (usmall_t j = 0; j < num_levels_; j++) {
                    chunk_lengths[j] = chunk_length_;
                }

                auto dacs = new SplitDacs<V, S, MAX_LEVELS, FACTOR_RANK>(
                        length_,
                        levels_,
                        num_levels_,
                        accum_level_sizes, // takes ownership
                        extra_bits_ranks, // takes ownership
                        level_ranks, // takes ownership
                        chunk_lengths, // takes ownership
                        min_values_
                );

                num_levels_ = 0;
                levels_limit_ = 0;
                levels_ = nullptr;
                min_values_ = nullptr;
                accum_level_sizes = nullptr;
                extra_bits_ranks = nullptr;
                level_ranks = nullptr;
                chunk_lengths = nullptr;

                if (opt) dacs->optimise();

                finalised_ = true;

                return dacs;

            }

        }

        /*
         * Delete the current contents of the builder and prepare it for building another SplitDacs instance.
         */
        void reset(usmall_t chunk_length, ularge_t capacity, float level_factor, float extend_factor) {

            clear();
            init(chunk_length, capacity, level_factor, extend_factor);

        }

        /*
         * Determine the memory consumption (in bytes).
         */
        ularge_t size_in_bytes() const {

            ularge_t contents = 0;
            for (usmall_t i = 0; i < num_levels_; i++) {

                contents += (chunk_length_ * level_capacities_[i] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); // levels_
                contents += (level_capacities_[i] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); // level_extra_bits_

            }

            return sizeof(SplitDacsDirectBuilder) // SplitDacsDirectBuilder itself
                + sizeof(S*) * levels_limit_ // levels_
                + sizeof(V) * levels_limit_ // min_values_
                + sizeof(ularge_t) * levels_limit_ // level_cursors_
                + sizeof(ularge_t) * levels_limit_ // chunk_cursors_
                + sizeof(ularge_t) * levels_limit_ // level_capacities_
                + sizeof(S*) * levels_limit_ // level_extra_bits_
                + sizeof(S) * contents; // actual contents of levels_ and level_extra_bits_

        }

    private:
        // members of SplitDacs instance to be built (see documentation in SplitDacs)
        ularge_t length_;
        S** levels_;
        usmall_t num_levels_;
        V* min_values_;

        // temporary helper members needed during the construction
        ularge_t* level_cursors_; // position (in bits) where to write the next chunk (per level)
        ularge_t* chunk_cursors_; // position (in chunks) of the next chunk to be written (per level)
        ularge_t len_extra_bits_; // length of bit vector for Vbyte extra bits
        usmall_t chunk_length_; // chunk length (in bits) used on all levels
        V base_; // number of possibles values on each level (based on the chunk length)

        ularge_t* level_capacities_; // maximum number of chunks for each level
        S** level_extra_bits_; // extra bits for each level
        usmall_t levels_limit_; // maximum number of levels that can exist for current chunk length and type V

        bool finalised_; // flag indicating that the builder is already been used to construct the SplitDacs

        float level_factor_; // factor for size of new level
        float extend_factor_; // factor for extending a full level

        /*
         * Initialise the builder with the provided chunk length and an initial capacity.
         */
        void init(usmall_t chunk_length, ularge_t capacity, float level_factor, float extend_factor) {

            length_ = 0;
            level_factor_ = level_factor;
            extend_factor_ = extend_factor;

            num_levels_ = 1;
            levels_limit_ = (bits(std::numeric_limits<V>::max()) + chunk_length - 1) / chunk_length;

            // determine number of levels and their chunk lengths
            chunk_length_ = chunk_length;

            // for each level, determine the minimum value that makes a chunk on the level necessary
            min_values_ = new V[levels_limit_];
            level_capacities_ = new ularge_t[levels_limit_];
            V value = 0;
            V mult_val = 1;
            for (usmall_t i = 0; i < levels_limit_; i++) {

                mult_val *= static_cast<V>(1) << chunk_length_;
                min_values_[i] = value;
                value += mult_val;

                level_capacities_[i] = 0;

            }
            level_capacities_[0] = capacity;

            // initially, only the first level exists, higher levels are added later as needed
            // keep a bit vector per level for Vbyte extra bits
            levels_ = new S*[levels_limit_];
            levels_[0] = new S[(chunk_length_ * level_capacities_[0] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];

            level_extra_bits_ = new S*[levels_limit_];
            level_extra_bits_[0] = new S[(level_capacities_[0] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
            len_extra_bits_ = 0;
            for (ularge_t i = 0; i < (level_capacities_[0] + (sizeof(S) * 8) - 1) / (sizeof(S) * 8); i++) {
                level_extra_bits_[0][i] = 0;
            }

            for (usmall_t i = 1; i < levels_limit_; i++) {

                levels_[i] = 0;
                level_extra_bits_[i] = 0;

            }

            // set up cursors into the levels for filling tem
            level_cursors_ = new ularge_t[levels_limit_]; // for each level, position (in bits) where to write the next chunk
            chunk_cursors_ = new ularge_t[levels_limit_]; // for each level, position (in chunks) of the next chunk to be written

            for (usmall_t j = 0; j < levels_limit_; j++) {

                level_cursors_[j] = 0;
                chunk_cursors_[j] = 0;

            }

            // number of possibles values for a chunk on the different levels (used to write the chunks below)
            base_ = static_cast<V>(1) << chunk_length_;

            finalised_ = false;

        }

        /*
         * Delete the current contents of the builder (when destructing or resetting it).
         */
        void clear() {

            for (usmall_t i = 0; i < levels_limit_; i++) {
                delete[] levels_[i];
            }
            delete[] levels_;
            delete[] min_values_;

            delete[] level_cursors_;
            delete[] chunk_cursors_;

            delete[] level_capacities_;
            for (usmall_t i = 0; i < levels_limit_; i++) {
                delete[] level_extra_bits_[i];
            }
            delete[] level_extra_bits_;

        }

    };

}

#endif //DACS_SPLITDACSBUILDER_HPP
