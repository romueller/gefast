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

#ifndef DACS_UTILITY_HPP
#define DACS_UTILITY_HPP

#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

namespace GeFaST {

    /* data types */

    // unsigned integers
    typedef unsigned long ularge_t;
    typedef unsigned short usmall_t;

    // signed integers
    typedef signed long slarge_t;
    typedef signed short ssmall_t;

    // result of DACs optimisation
    struct Chunks {

        Chunks(usmall_t* c, usmall_t n);

        usmall_t* chunk_lengths = nullptr; // chunk lengths on the different levels (not owned by the Chunks instance)
        usmall_t num_levels = 0; // number of levels

    };


    /*
     * Diverse methods for bit vector manipulations.
     *
     * Based on similar methods from the original DACs implementation:
     * DACs: Bringing Direct Access to Variable-Length Codes (Brisaboa et al., 2013)
     * http://lbd.udc.es/research/DACS/
     *
     * However, bits are now always written consecutively - unlike the original methods,
     * which break a bit pattern into two non-consecutive parts in case it overlaps with
     * two words of the bit vector representation.
     */

    /*
     * Determine the number of bits needed to represent the value in 'val'.
     */
    template<typename V>
    usmall_t bits(V val) {

        usmall_t num_bits = 0;
        while (val > 0) {

            num_bits++;
            val >>= 1;

        }

        return num_bits;

    }


    /*
     * Retrieve bit arr[pos].
     */
    template<typename S>
    inline bool get_bit(S* arr, ularge_t pos) {
        return (arr[pos / (sizeof(S) * 8)] >> ((sizeof(S) * 8) - 1 - pos % (sizeof(S) * 8))) & 1;
    }

    /*
     * Set arr[pos] to 1.
     */
    template<typename S>
    inline void set_bit(S* arr, ularge_t pos) {
        arr[pos / (sizeof(S) * 8)] |= static_cast<S>(1) << ((sizeof(S) * 8) - 1 - pos % (sizeof(S) * 8));
    }

    /*
     * Set arr[pos] to 0.
     */
    template<typename S>
    inline void reset_bit(S* arr, ularge_t pos) {
        arr[pos / (sizeof(S) * 8)] &= ~(static_cast<S>(1) << ((sizeof(S) * 8) - 1 - pos % (sizeof(S) * 8)));
    }


    /*
     * Retrieve the bits arr[p : p + len - 1] as an integer.
     *
     * Assumptions: indices not out of bounds, len <= (sizeof(S) * 8)
     */
    template<typename V, typename S>
    V read_bits(S* arr, ularge_t p, ularge_t len) {

        // initialise reading position
        arr += p / (sizeof(S) * 8);
        p %= (sizeof(S) * 8);

        V res;

        if (p + len <= (sizeof(S) * 8)) {
            res = (*arr >> ((sizeof(S) * 8) - (p + len))) & ((static_cast<S>(1) << len) - 1);
        } else {

            len = p + len - (sizeof(S) * 8); // number of overhanging bits

            res = ((*arr & ((static_cast<S>(1) << ((sizeof(S) * 8) - p)) - 1)) << len) // read non-overhanging bits and shift them to the right position in the result
                   | ((*(arr + 1) >> ((sizeof(S) * 8) - len)) & ((static_cast<S>(1) << len) - 1)); // read overhanging bits

        }

        return res;

    }

    /*
     * Write arr[p : p + len - 1] = val.
     *
     * Assumptions: indices not out of bounds, len <= (sizeof(S) * 8)
     */
    template<typename V, typename S>
    void write_bits(S* arr, ularge_t p, ularge_t len, V val) {

        // initialise writing position
        arr += p / (sizeof(S) * 8);
        p %= (sizeof(S) * 8);

        // write bits into one or two (neighbouring) words
        if (p + len <= (sizeof(S) * 8)) {
            *arr = (*arr & ~(((static_cast<S>(1) << len) - 1) << ((sizeof(S) * 8) - (p + len)))) | (static_cast<S>(val) << ((sizeof(S) * 8) - (p + len)));
        } else {

            len = p + len - (sizeof(S) * 8); // number of overhanging bits

            *arr = (*arr & ~((static_cast<S>(1) << ((sizeof(S) * 8) - p)) - 1)) | (static_cast<S>(val) >> len); // write non-overhanging bits

            arr++;

            *arr = (*arr & ((static_cast<S>(1) << ((sizeof(S) * 8) - len)) - 1)) | (static_cast<S>(val) << ((sizeof(S) * 8) - len)); // write overhanging bits

        }

    }

    /*
     * Insert len bits at position p in arr.
     *
     * Assumptions: indices not out of bounds, len <= (sizeof(S) * 8), capacity of arr is not exceeded by insertion
     */
    template<typename V, typename S>
    void insert_bits(S* arr, ularge_t arr_len, V val, ularge_t p, ularge_t len) {

        S* end = arr + arr_len;

        // initialise insertion position
        arr += p / (sizeof(S) * 8);
        p %= (sizeof(S) * 8);

        S tmp;
        S local_val = static_cast<S>(val);

        if (p + len <= (sizeof(S) * 8)) { // inserted bits completely within one word

            // first affected word: save bits that are pushed into the next word, insert val at position p
            // ternary operators avoid undefined behaviour (in case the right operand is >= width of S)
            tmp = (*arr << ((sizeof(S) * 8) - len)); // bits pushed into next word
            *arr = ((p) ? (*arr >> ((sizeof(S) * 8) - p)) << ((sizeof(S) * 8) - p) : 0) // kept bits before insertion position
                   | (local_val << ((sizeof(S) * 8) - len - p)) // bits of val
                   | ((p + len >= (sizeof(S) * 8)) ? 0 : ((*arr << p) >> (p + len))); // new last bits
            local_val = tmp;

        } else { // inserted bits overlap with next word

            // first affected word: save bits that are pushed into the next word, insert first bits of val at position p,
            // and combine remaining bits of val with saved bits to new value to be inserted
            tmp = *arr << p; tmp >>= p; // bits pushed into next word; note: doing it in one command is somehow buggy (no effect)
            *arr = (*arr >> ((sizeof(S) * 8) - p)) << ((sizeof(S) * 8) - p) // kept bits before insertion position
                   | (local_val >> (len - ((sizeof(S) * 8) - p))); // first bits of val
            local_val = ((local_val << ((sizeof(S) * 8) - p)) | tmp) << ((sizeof(S) * 8) - len); // last bits of val and saved bits

        }

        // subsequent words: save bits that are in turn pushed into the next word, insert previous value at position 0
        for (arr++; arr != end; arr++) {

            tmp = (*arr << ((sizeof(S) * 8) - len));
            *arr = (*arr >> len) | local_val;
            local_val = tmp;

        }

    }



    /*
     * Population-count methods as described in
     * "Practical implementation of rank and select queries"
     * (González et al., 2005)
     */


    // lookup table for one byte
    const std::array<unsigned char, 256> popcount_tab = {
            0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
            1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
            1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
            2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
            1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
            2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
            2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
            3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
    };

    // general population-count implementation for B-bit words (B = multiple of 8)
    template<usmall_t B>
    struct PopCounter {
        static usmall_t popcount(unsigned long long x) {

            usmall_t sum = 0;
            for (usmall_t shift = 0; shift < B; shift += 8) {
                sum += popcount_tab[(x >> shift) & 0xff];
            }

            return sum;

        };
    };

    // specialised population-count implementation for 8-bit words
    template<>
    struct PopCounter<8> {
        static usmall_t popcount(unsigned long long x) {
            return popcount_tab[x & 0xff];
        }
    };

    // specialised population-count implementation for 16-bit words
    template<>
    struct PopCounter<16> {
        static usmall_t popcount(unsigned long long x) {
            return popcount_tab[x & 0xff] + popcount_tab[(x >> 8) & 0xff];
        }
    };

    // specialised population-count implementation for 32-bit words
    template<>
    struct PopCounter<32> {
        static usmall_t popcount(unsigned long long x) {
            return popcount_tab[x & 0xff] + popcount_tab[(x >> 8) & 0xff]
                   + popcount_tab[(x >> 16) & 0xff] + popcount_tab[(x >> 24) & 0xff];
        }
    };

    // specialised population-count implementation for 64-bit words
    template<>
    struct PopCounter<64> {
        static usmall_t popcount(unsigned long long x) {
            return popcount_tab[x & 0xff] + popcount_tab[(x >> 8) & 0xff]
                   + popcount_tab[(x >> 16) & 0xff] + popcount_tab[(x >> 24) & 0xff]
                   + popcount_tab[(x >> 32) & 0xff] + popcount_tab[(x >> 40) & 0xff]
                   + popcount_tab[(x >> 48) & 0xff] + popcount_tab[(x >> 56) & 0xff];
        }
    };

    /*
     * Count the 1-bits in arr[0 : p[.
     *
     * Assumption: index not out of bounds
     */
    template<typename S>
    ularge_t count_bits(S* arr, ularge_t p) {

        ularge_t sum = 0;
        for (ularge_t i = 0; i < p / (sizeof(S) * 8); i++, arr++) {
            sum += PopCounter<sizeof(S) * 8>::popcount(*arr);
        }
        sum += PopCounter<sizeof(S) * 8>::popcount((*arr >> 1) & ~((static_cast<S>(1) << ((sizeof(S) * 8) - 1 - (p % (sizeof(S) * 8)))) - 1));

        return sum;

    }



    /*
     * Methods for optimisation of DACs.
     *
     * Based on similar methods from the original DACs implementation:
     * DACs: Bringing Direct Access to Variable-Length Codes (Brisaboa et al., 2013)
     * http://lbd.udc.es/research/DACS/
     */

    /*
     * The DACs representation consists of multiple levels, each handling chunks from a different position of
     * the represented integers. The length of the chunks on each level can be chosen separately.
     *
     * This method optimises the number of levels and the chunk lengths without a restriction on the number of levels.
     *
     * Template parameter 'FACTOR' is used to incorporate the contribution of the rank data structure to the overall
     * memory consumption.
     */
    template<typename V, usmall_t FACTOR>
    Chunks optimise_level_chunks(const V* values, ularge_t len) {

        V max_int = 0;
        for (ularge_t i = 0; i < len; i++) {
            max_int = std::max(max_int, values[i]);
        }

        usmall_t num_bits = bits(max_int);

        // special case: all values in the sequence are zero
        if (num_bits == 0) {

            usmall_t num_levels = 1;
            usmall_t* chunk_lengths = new usmall_t[num_levels];
            chunk_lengths[0] = 0;
            return Chunks(chunk_lengths, num_levels);

        }

        // count the number of integers < 2^i for 0 <= i <= numBits
        ularge_t* accum_cnts = new ularge_t[num_bits + 1];
        for (usmall_t i = 0; i <= num_bits; i++) {
            accum_cnts[i] = 0;
        }
        for (ularge_t i = 0; i < len; i++) {
            for (usmall_t j = 1; j < num_bits && values[i] < (static_cast<V>(1) << (num_bits - j)); j++) {
                accum_cnts[num_bits - j]++;
            }
        }
        accum_cnts[num_bits] = len; // in addition, store total number of elements as the last entry

        /* dynamic-programming scheme */

        ularge_t* s = new ularge_t[num_bits]; // optimal size of the encoding for the subproblems
        usmall_t* l = new usmall_t[num_bits]; // optimal number of levels for the subproblems
        usmall_t* b = new usmall_t[num_bits]; // size of the blocks for the first level of the optimal subdivision for the subproblems

        s[num_bits - 1] = 0;
        l[num_bits - 1] = 0;
        b[num_bits - 1] = 0;

        ularge_t inf = std::numeric_limits<ularge_t>::max();
        for (ssmall_t t = num_bits - 1; t >= 0; t--) { // consider subproblem t

            ularge_t min_size = inf;
            usmall_t min_pos = num_bits - 1;

            for (usmall_t i = num_bits - 1; i > t; i--) { // check combination with solution of subproblem i > t

                ularge_t current_size = s[i] + (accum_cnts[num_bits] - accum_cnts[t]) * (i - t + 1) + (accum_cnts[num_bits] - accum_cnts[t]) / FACTOR;

                if (min_size > current_size) {

                    min_size = current_size;
                    min_pos = i;

                }

            }

            // check whether creating a single level is better
            if (min_size < ((accum_cnts[num_bits] - accum_cnts[t]) * (num_bits - t))) {

                s[t] = min_size;
                l[t] = l[min_pos] + static_cast<usmall_t>(1);
                b[t] = min_pos - t;

            } else {

                s[t] = (accum_cnts[num_bits] - accum_cnts[t]) * (num_bits - t);
                l[t] = 1;
                b[t] = num_bits - t;

            }

        }

        // construct optimal solution from the arrays l and b
        usmall_t num_levels = l[0];
        usmall_t* chunk_lengths = new usmall_t[num_levels];
        for (ularge_t k = 0, t = 0; k < num_levels; k++, t += b[t]) {
            chunk_lengths[k] = b[t];
        }

        delete[] b;
        delete[] l;
        delete[] s;
        delete[] accum_cnts;

        return Chunks(chunk_lengths, num_levels);

    }

    /*
     * Same functionality as optimise_level_chunks(V* values, ...) above, but the values are handed over via
     * a container of type C, which has to offer the access operation ([], zero-based).
     */
    template<typename C, typename V, usmall_t FACTOR>
    Chunks optimise_level_chunks(const C& values, ularge_t len) {

        V max_int = 0;
        for (ularge_t i = 0; i < len; i++) {
            max_int = std::max(max_int, values[i]);
        }

        usmall_t num_bits = bits(max_int);

        // special case: all values in the sequence are zero
        if (num_bits == 0) {

            usmall_t num_levels = 1;
            usmall_t* chunk_lengths = new usmall_t[num_levels];
            chunk_lengths[0] = 0;
            return Chunks(chunk_lengths, num_levels);

        }

        // count the number of integers < 2^i for 0 <= i <= numBits
        ularge_t* accum_cnts = new ularge_t[num_bits + 1];
        for (usmall_t i = 0; i <= num_bits; i++) {
            accum_cnts[i] = 0;
        }
        for (ularge_t i = 0; i < len; i++) {
            for (usmall_t j = 1; j < num_bits && values[i] < (static_cast<V>(1) << (num_bits - j)); j++) {
                accum_cnts[num_bits - j]++;
            }
        }
        accum_cnts[num_bits] = len; // in addition, store total number of elements as the last entry

        /* dynamic-programming scheme */

        ularge_t* s = new ularge_t[num_bits]; // optimal size of the encoding for the subproblems
        usmall_t* l = new usmall_t[num_bits]; // optimal number of levels for the subproblems
        usmall_t* b = new usmall_t[num_bits]; // size of the blocks for the first level of the optimal subdivision for the subproblems

        s[num_bits - 1] = 0;
        l[num_bits - 1] = 0;
        b[num_bits - 1] = 0;

        ularge_t inf = std::numeric_limits<ularge_t>::max();
        for (ssmall_t t = num_bits - 1; t >= 0; t--) { // consider subproblem t

            ularge_t min_size = inf;
            usmall_t min_pos = num_bits - 1;

            for (usmall_t i = num_bits - 1; i > t; i--) { // check combination with solution of subproblem i > t

                ularge_t current_size = s[i] + (accum_cnts[num_bits] - accum_cnts[t]) * (i - t + 1) + (accum_cnts[num_bits] - accum_cnts[t]) / FACTOR;

                if (min_size > current_size) {

                    min_size = current_size;
                    min_pos = i;

                }

            }

            // check whether creating a single level is better
            if (min_size < ((accum_cnts[num_bits] - accum_cnts[t]) * (num_bits - t))) {

                s[t] = min_size;
                l[t] = l[min_pos] + static_cast<usmall_t>(1);
                b[t] = min_pos - t;

            } else {

                s[t] = (accum_cnts[num_bits] - accum_cnts[t]) * (num_bits - t);
                l[t] = 1;
                b[t] = num_bits - t;

            }

        }

        // construct optimal solution from the arrays l and b
        usmall_t num_levels = l[0];
        usmall_t* chunk_lengths = new usmall_t[num_levels];
        for (ularge_t k = 0, t = 0; k < num_levels; k++, t += b[t]) {
            chunk_lengths[k] = b[t];
        }

        delete[] b;
        delete[] l;
        delete[] s;
        delete[] accum_cnts;

        return Chunks(chunk_lengths, num_levels);

    }

    /*
     * Same functionality as optimise_level_chunks(V* values, ...) above, but the values are handed over
     * as a distribution (each occurring value is associated with the number of times it occurs).
     */
    template<typename V, usmall_t FACTOR>
    Chunks optimise_level_chunks(const std::map<ularge_t, ularge_t>& counts) {

        V max_int = counts.rbegin()->first;
        usmall_t num_bits = bits(max_int);

        // special case: all values in the sequence are zero
        if (num_bits == 0) {

            usmall_t num_levels = 1;
            usmall_t* chunk_lengths = new usmall_t[num_levels];
            chunk_lengths[0] = 0;
            return Chunks(chunk_lengths, num_levels);

        }

        // count the number of integers < 2^i for 0 <= i <= numBits
        ularge_t* accum_cnts = new ularge_t[num_bits + 1];
        ularge_t accum_value = 0;
        ularge_t cnt_bits = 1;

        accum_cnts[0] = 0;
        for (ularge_t i = 0; i <= max_int; i++) {

            if (i == (static_cast<ularge_t>(1) << cnt_bits)) { // next power of 2 reached -> store accumulated count

                accum_cnts[cnt_bits] = accum_value;
                cnt_bits++;

            }

            auto iter = counts.find(i);
            accum_value += (iter != counts.end()) ? iter->second : 0;

        }
        accum_cnts[cnt_bits] = accum_value; // in addition, store total number of elements as the last entry


        /* dynamic-programming scheme */

        ularge_t* s = new ularge_t[num_bits]; // optimal size of the encoding for the subproblems
        usmall_t* l = new usmall_t[num_bits]; // optimal number of levels for the subproblems
        usmall_t* b = new usmall_t[num_bits]; // size of the blocks for the first level of the optimal subdivision for the subproblems

        s[num_bits - 1] = 0;
        l[num_bits - 1] = 0;
        b[num_bits - 1] = 0;

        ularge_t inf = std::numeric_limits<ularge_t>::max();
        for (ssmall_t t = num_bits - 1; t >= 0; t--) { // consider subproblem t

            ularge_t min_size = inf;
            usmall_t min_pos = num_bits - 1;

            for (usmall_t i = num_bits - 1; i > t; i--) { // check combination with solution of subproblem i > t

                ularge_t current_size = s[i] + (accum_cnts[num_bits] - accum_cnts[t]) * (i - t + 1) + (accum_cnts[num_bits] - accum_cnts[t]) / FACTOR;

                if (min_size > current_size) {

                    min_size = current_size;
                    min_pos = i;

                }

            }

            // check whether creating a single level is better
            if (min_size < ((accum_cnts[num_bits] - accum_cnts[t]) * (num_bits - t))) {

                s[t] = min_size;
                l[t] = l[min_pos] + static_cast<usmall_t>(1);
                b[t] = min_pos - t;

            } else {

                s[t] = (accum_cnts[num_bits] - accum_cnts[t]) * (num_bits - t);
                l[t] = 1;
                b[t] = num_bits - t;

            }

        }

        // construct optimal solution from the arrays l and b
        usmall_t num_levels = l[0];
        usmall_t* chunk_lengths = new usmall_t[num_levels];
        for (ularge_t k = 0, t = 0; k < num_levels; k++, t += b[t]) {
            chunk_lengths[k] = b[t];
        }

        delete[] b;
        delete[] l;
        delete[] s;
        delete[] accum_cnts;

        return Chunks(chunk_lengths, num_levels);

    }

    /*
     * Helper method for optimise_level_chunks_limited containing the actual dynamic-programming computation
     * obeying the maximum number of levels.
     *
     * Performs the computation iteratively (instead of recursively as in the original DACs implementation).
     */
    template<usmall_t FACTOR>
    void optimise_for_limit(ularge_t** opt_sizes, usmall_t** opt_levels, usmall_t** opt_first_blocks, usmall_t max_levels, usmall_t num_bits, const ularge_t* accum_cnts) {

        // subproblems t, using only one level
        for (usmall_t t = 0; t < num_bits; t++) {

            opt_sizes[t][0] = (accum_cnts[num_bits] - accum_cnts[t]) * (num_bits - t);
            opt_levels[t][0] = 1;
            opt_first_blocks[t][0] = num_bits - t;

        }

        ularge_t inf = std::numeric_limits<ularge_t>::max();

        // subproblems t, using > 1 level
        for (usmall_t v = 1; v < max_levels; v++) {

            for (usmall_t t = 0; t < num_bits; t++) { // consider subproblem t

                ularge_t min_size = inf;
                usmall_t min_pos = num_bits - 1;

                for (usmall_t i = num_bits - 1; i > t; i--) { // check combination with solution of subproblem i > t and one level less

                    ularge_t current_size = opt_sizes[i][v - 1] + (accum_cnts[num_bits] - accum_cnts[t]) * (i - t + 1) + (accum_cnts[num_bits] - accum_cnts[t]) / FACTOR;

                    if (min_size > current_size) {

                        min_size = current_size;
                        min_pos = i;

                    }

                }

                // check whether creating a single level is better
                if (min_size < ((accum_cnts[num_bits] - accum_cnts[t]) * (num_bits - t))) {

                    opt_sizes[t][v] = min_size;
                    opt_levels[t][v] = opt_levels[min_pos][v - 1] + static_cast<usmall_t>(1);
                    opt_first_blocks[t][v] = min_pos - t;

                } else {

                    opt_sizes[t][v] = (accum_cnts[num_bits] - accum_cnts[t]) * (num_bits - t);
                    opt_levels[t][v] = 1;
                    opt_first_blocks[t][v] = num_bits - t;

                }

            }

        }

    }

    /*
     * The DACs representation consists of multiple levels, each handling chunks from a different position of
     * the represented integers. The length of the chunks on each level can be chosen separately.
     *
     * This method optimises the number of levels and the chunk lengths while not using more than 'maxLevels' levels.
     *
     * Template parameter 'FACTOR' is used to incorporate the contribution of the rank data structure to the overall
     * memory consumption.
     */
    template<typename V, usmall_t FACTOR>
    Chunks optimise_level_chunks_limited(const V* values, ularge_t len, usmall_t max_levels) {

        V max_int = 0;
        for (ularge_t i = 0; i < len; i++) {
            max_int = std::max(max_int, values[i]);
        }

        usmall_t num_bits = bits(max_int);

        // special case: all values in the sequence are zero
        if (num_bits == 0) {

            usmall_t num_levels = 1;
            usmall_t* chunk_lengths = new usmall_t[num_levels];
            chunk_lengths[0] = 0;
            return Chunks(chunk_lengths, num_levels);

        }

        // count the number of integers < 2^i for 0 <= i <= numBits
        ularge_t* accum_cnts = new ularge_t[num_bits + 1];
        for (usmall_t i = 0; i <= num_bits; i++) {
            accum_cnts[i] = 0;
        }
        for (ularge_t i = 0; i < len; i++) {
            for (usmall_t j = 1; j < num_bits && values[i] < (static_cast<V>(1) << (num_bits - j)); j++) {
                accum_cnts[num_bits - j]++;
            }
        }
        accum_cnts[num_bits] = len; // in addition, store total number of elements as the last entry

        /* dynamic-programming scheme */

        ularge_t** s = new ularge_t*[num_bits];
        usmall_t** l = new usmall_t*[num_bits];
        usmall_t** b = new usmall_t*[num_bits];

        for (usmall_t t = 0; t < num_bits; t++) {

            s[t] = new ularge_t[max_levels];
            l[t] = new usmall_t[max_levels];
            b[t] = new usmall_t[max_levels];

        }

        optimise_for_limit<FACTOR>(s, l, b, max_levels, num_bits, accum_cnts);

        // construct solution
        usmall_t num_levels = l[0][max_levels - 1];

        usmall_t* chunk_lengths = new usmall_t[num_levels];
        for (usmall_t t = 0, k = 0; k < num_levels; t += b[t][max_levels - 1 - k], k++) {
            chunk_lengths[k] = b[t][max_levels - 1 - k];
        }

        for (usmall_t t = 0; t < num_bits; t++) {

            delete[] s[t];
            delete[] l[t];
            delete[] b[t];

        }
        delete[] b;
        delete[] l;
        delete[] s;
        delete[] accum_cnts;

        return Chunks(chunk_lengths, num_levels);

    }

    /*
     * Same functionality as optimise_level_chunks_limited(V* values, ...) above, but the values are handed over via
     * a container of type C, which has to offer the access operation ([], zero-based).
     */
    template<typename C, typename V, usmall_t FACTOR>
    Chunks optimise_level_chunks_limited(const C& values, ularge_t len, usmall_t max_levels) {

        V max_int = 0;
        for (ularge_t i = 0; i < len; i++) {
            max_int = std::max(max_int, values[i]);
        }

        usmall_t num_bits = bits(max_int);

        // special case: all values in the sequence are zero
        if (num_bits == 0) {

            usmall_t num_levels = 1;
            usmall_t* chunk_lengths = new usmall_t[num_levels];
            chunk_lengths[0] = 0;
            return Chunks(chunk_lengths, num_levels);

        }

        // count the number of integers < 2^i for 0 <= i <= numBits
        ularge_t* accum_cnts = new ularge_t[num_bits + 1];
        for (usmall_t i = 0; i <= num_bits; i++) {
            accum_cnts[i] = 0;
        }
        for (ularge_t i = 0; i < len; i++) {
            for (usmall_t j = 1; j < num_bits && values[i] < (static_cast<V>(1) << (num_bits - j)); j++) {
                accum_cnts[num_bits - j]++;
            }
        }
        accum_cnts[num_bits] = len; // in addition, store total number of elements as the last entry

        /* dynamic-programming scheme */

        ularge_t** s = new ularge_t*[num_bits];
        usmall_t** l = new usmall_t*[num_bits];
        usmall_t** b = new usmall_t*[num_bits];

        for (usmall_t t = 0; t < num_bits; t++) {

            s[t] = new ularge_t[max_levels];
            l[t] = new usmall_t[max_levels];
            b[t] = new usmall_t[max_levels];

        }

        optimise_for_limit<FACTOR>(s, l, b, max_levels, num_bits, accum_cnts);

        // construct solution
        usmall_t num_levels = l[0][max_levels - 1];

        usmall_t* chunk_lengths = new usmall_t[num_levels];
        for (usmall_t t = 0, k = 0; k < num_levels; t += b[t][max_levels - 1 - k], k++) {
            chunk_lengths[k] = b[t][max_levels - 1 - k];
        }

        for (usmall_t t = 0; t < num_bits; t++) {

            delete[] s[t];
            delete[] l[t];
            delete[] b[t];

        }
        delete[] b;
        delete[] l;
        delete[] s;
        delete[] accum_cnts;

        return Chunks(chunk_lengths, num_levels);

    }

    /*
     * Same functionality as optimise_level_chunks_limited(V* values, ...) above, but the values are handed over
     * as a distribution (each occurring value is associated with the number of times it occurs).
     */
    template<typename V, usmall_t FACTOR>
    Chunks optimise_level_chunks_limited(const std::map<ularge_t, ularge_t>& counts, usmall_t max_levels) {

        V max_int = counts.rbegin()->first;
        usmall_t num_bits = bits(max_int);

        // special case: all values in the sequence are zero
        if (num_bits == 0) {

            usmall_t num_levels = 1;
            usmall_t* chunk_lengths = new usmall_t[num_levels];
            chunk_lengths[0] = 0;
            return Chunks(chunk_lengths, num_levels);

        }

        // count the number of integers < 2^i for 0 <= i <= numBits
        ularge_t* accum_cnts = new ularge_t[num_bits + 1];
        ularge_t accum_value = 0;
        ularge_t cnt_bits = 1;

        accum_cnts[0] = 0;
        for (ularge_t i = 0; i <= max_int; i++) {

            if (i == (static_cast<ularge_t>(1) << cnt_bits)) { // next power of 2 reached -> store accumulated count

                accum_cnts[cnt_bits] = accum_value;
                cnt_bits++;

            }

            auto iter = counts.find(i);
            accum_value += (iter != counts.end()) ? iter->second : 0;

        }
        accum_cnts[cnt_bits] = accum_value; // in addition, store total number of elements as the last entry


        /* dynamic-programming scheme */

        ularge_t** s = new ularge_t*[num_bits];
        usmall_t** l = new usmall_t*[num_bits];
        usmall_t** b = new usmall_t*[num_bits];

        for (usmall_t t = 0; t < num_bits; t++) {

            s[t] = new ularge_t[max_levels];
            l[t] = new usmall_t[max_levels];
            b[t] = new usmall_t[max_levels];

        }

        optimise_for_limit<FACTOR>(s, l, b, max_levels, num_bits, accum_cnts);

        // construct solution
        usmall_t num_levels = l[0][max_levels - 1];

        usmall_t* chunk_lengths = new usmall_t[num_levels];
        for (usmall_t t = 0, k = 0; k < num_levels; t += b[t][max_levels - 1 - k], k++) {
            chunk_lengths[k] = b[t][max_levels - 1 - k];
        }

        for (usmall_t t = 0; t < num_bits; t++) {

            delete[] s[t];
            delete[] l[t];
            delete[] b[t];

        }
        delete[] b;
        delete[] l;
        delete[] s;
        delete[] accum_cnts;

        return Chunks(chunk_lengths, num_levels);

    }

    /*
     * Helper data structure for recording the sizes of integers to determine the optimal chunk lengths
     * for a DACs representation (without building the DACs).
     */
    template<typename V, usmall_t FACTOR>
    struct ChunkOptimiser {

        /*
         * Initialise empty integer sequence and prepare counting the integer sizes.
         */
        ChunkOptimiser() {

            // member accumCnts_ is a fixed-size array
            for (size_t i = 0; i <= sizeof(V) * 8; i++) {
                accum_cnts_[i] = 0;
            }

            size_ = 0;
            max_val_ = 0;

        }

        /*
         * Add integer 'val' ('inc' times) to the sequence.
         */
        void add(V val, ularge_t inc = 1) {

            for (usmall_t j = 1; j < sizeof(V) * 8 && val < (static_cast<usmall_t>(1) << (sizeof(V) * 8 - j)); j++) {
                accum_cnts_[sizeof(V) * 8 - j] += inc;
            }

            size_ += inc;
            max_val_ = std::max(max_val_, static_cast<ularge_t>(val));

        }

        /*
         * Determine the optimal chunk lengths based on the recording integers.
         */
        Chunks optimise() {

            usmall_t num_bits = bits(max_val_);

            // special case: all values in the sequence are zero
            if (num_bits == 0) {

                usmall_t num_levels = 1;
                usmall_t* chunk_lengths = new usmall_t[num_levels];
                chunk_lengths[0] = 0;
                return Chunks(chunk_lengths, num_levels);

            }

            accum_cnts_[num_bits] = size_; // in addition, store total number of elements as the last entry

            /* dynamic-programming scheme */

            ularge_t* s = new ularge_t[num_bits]; // optimal size of the encoding for the subproblems
            usmall_t* l = new usmall_t[num_bits]; // optimal number of levels for the subproblems
            usmall_t* b = new usmall_t[num_bits]; // size of the blocks for the first level of the optimal subdivision for the subproblems

            s[num_bits - 1] = 0;
            l[num_bits - 1] = 0;
            b[num_bits - 1] = 0;

            ularge_t inf = std::numeric_limits<ularge_t>::max();
            for (ssmall_t t = num_bits - 1; t >= 0; t--) { // consider subproblem t

                ularge_t min_size = inf;
                usmall_t min_pos = num_bits - 1;

                for (usmall_t i = num_bits - 1; i > t; i--) { // check combination with solution of subproblem i > t

                    ularge_t current_size = s[i] + (accum_cnts_[num_bits] - accum_cnts_[t]) * (i - t + 1) + (accum_cnts_[num_bits] - accum_cnts_[t]) / FACTOR;

                    if (min_size > current_size) {

                        min_size = current_size;
                        min_pos = i;

                    }

                }

                // check whether creating a single level is better
                if (min_size < ((accum_cnts_[num_bits] - accum_cnts_[t]) * (num_bits - t))) {

                    s[t] = min_size;
                    l[t] = l[min_pos] + 1;
                    b[t] = min_pos - t;

                } else {

                    s[t] = (accum_cnts_[num_bits] - accum_cnts_[t]) * (num_bits - t);
                    l[t] = 1;
                    b[t] = num_bits - t;

                }

            }

            // construct optimal solution from the arrays l and b
            usmall_t num_levels = l[0];
            usmall_t* chunk_lengths = new usmall_t[num_levels];
            for (ularge_t k = 0, t = 0; k < num_levels; k++, t += b[t]) {
                chunk_lengths[k] = b[t];
            }

            accum_cnts_[num_bits] = 0;

            return Chunks(chunk_lengths, num_levels);

        }

        ularge_t accum_cnts_[sizeof(V) * 8 + 1]; // counts of integers < 2^i for 0 <= i <= num_bits in the sequence
        ularge_t size_; // length of the integer sequence
        ularge_t max_val_; // largest value in the integer sequence

    };

    /*
     * Helper data structure providing the configuration / major parameters of a SplitDacs
     * determined from a collection of counts to initialise the SplitDacs instance appropriately.
     */
    template<typename V, typename S, usmall_t MAX_LEVELS = 0, usmall_t FACTOR_RANK = 20>
    struct SplitDacsParameters {

        ularge_t length; // number of represented integers
        usmall_t num_levels; // number of levels / streams
        usmall_t* chunk_lengths; // chunk length (in bits) for each level
        V* min_values; // for each level, minimum value that leads to a chunk on that level
        ularge_t* accum_level_sizes; // sum of the sizes of the levels (in chunks) before each level

        /*
         * Determine the major parameters from the distribution of the elements
         * (provided as a collection of pairs (value, number of occurrences).
         */
        explicit SplitDacsParameters(const std::map<ularge_t, ularge_t>& counts) {

            ularge_t i;

            length = 0;
            for (auto& c : counts) {
                length += c.second;
            }

            // determine number of levels and their chunk lengths
            Chunks chunks = (MAX_LEVELS == 0) ?
                            optimise_level_chunks<V, FACTOR_RANK>(counts)
                          : optimise_level_chunks_limited<V, FACTOR_RANK>(counts, MAX_LEVELS);
            chunk_lengths = chunks.chunk_lengths;
            num_levels = chunks.num_levels;

            // for each level, determine the minimum value that makes a chunk on the level necessary
            min_values = new V[num_levels];
            accum_level_sizes = new ularge_t[num_levels + 1];
            V value = 0;
            V mult_val = 1;
            for (i = 0; i < num_levels; i++) {

                mult_val *= static_cast<V>(1) << chunk_lengths[i];
                min_values[i] = value;
                value += mult_val;

                accum_level_sizes[i] = 0;

            }
            accum_level_sizes[num_levels] = 0;

            // for each level, determine how many chunks will be required
            for (auto& c : counts) {

                value = c.first;
                for (usmall_t j = 0; j < num_levels; j++) {
                    accum_level_sizes[j + 1] += c.second * (value >= min_values[j]);
                }

            }

        }

        ~SplitDacsParameters() {

            delete[] chunk_lengths;
            delete[] min_values;
            delete[] accum_level_sizes;

        }

        SplitDacsParameters(const SplitDacsParameters& other) { // copy constructor

            length = other.length;
            num_levels = other.num_levels;

            chunk_lengths = new usmall_t[num_levels];
            min_values = new V[num_levels];
            accum_level_sizes = new ularge_t[num_levels + 1];
            for (auto i = 0; i < num_levels; i++) {

                chunk_lengths[i] = other.chunk_lengths[i];
                min_values[i] = other.min_values[i];
                accum_level_sizes[i] = other.accum_level_sizes[i];

            }
            accum_level_sizes[num_levels] = other.accum_level_sizes[num_levels];

        }

        SplitDacsParameters(SplitDacsParameters&& other) noexcept { // move constructor

            length = other.length;
            num_levels = other.num_levels;

            chunk_lengths = other.chunk_lengths; other.chunk_lengths = nullptr;
            min_values = other.min_values; other.min_values = nullptr;
            accum_level_sizes = other.accum_level_sizes; other.accum_level_sizes = nullptr;

        }

        SplitDacsParameters& operator=(const SplitDacsParameters& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] chunk_lengths;
            delete[] min_values;
            delete[] accum_level_sizes;

            // copy new resources
            length = other.length;
            num_levels = other.num_levels;

            chunk_lengths = new usmall_t[num_levels];
            min_values = new V[num_levels];
            accum_level_sizes = new ularge_t[num_levels + 1];
            for (auto i = 0; i < num_levels; i++) {

                chunk_lengths[i] = other.chunk_lengths[i];
                min_values[i] = other.min_values[i];
                accum_level_sizes[i] = other.accum_level_sizes[i];

            }
            accum_level_sizes[num_levels] = other.accum_level_sizes[num_levels];

            return *this;

        }

        SplitDacsParameters& operator=(SplitDacsParameters&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] chunk_lengths;
            delete[] min_values;
            delete[] accum_level_sizes;

            // copy / transfer new resources
            length = other.length;
            num_levels = other.num_levels;

            chunk_lengths = other.chunk_lengths; other.chunk_lengths = nullptr;
            min_values = other.min_values; other.min_values = nullptr;
            accum_level_sizes = other.accum_level_sizes; other.accum_level_sizes = nullptr;

            return *this;

        }

        SplitDacsParameters* clone() const { // deep-copy clone method
            return new SplitDacsParameters(*this);
        }

    };


    /*
     * Bit vector data structure.
     *
     * The bit vector consists of blocks of sizeof(S) * 8 bits.
     */
    template<typename S>
    class BitVector {

    public:
        // create a bit vector of length 'size' with all bits set to 'value'
        BitVector(size_t size, bool value) {

            size_ = size;
            alloc_size_ = (size_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8);
            bits_ = new S[alloc_size_];

            S val = value? -1 : 0;
            for (size_t i = 0; i < alloc_size_; i++) {
                bits_[i] = val;
            }

        }

        // create a bit vector of length 'size' with all bits set to false
        explicit BitVector(size_t size = 0) : BitVector(size, false) {
            // nothing else to do
        }

        // create a bit vector from std::vector<bool>
        explicit BitVector(const std::vector<bool>& vec) {

            size_ = vec.size();
            alloc_size_ = (size_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8);
            bits_ = new S[alloc_size_];
            for (size_t i = 0; i < alloc_size_; i++) {
                bits_[i] = 0;
            }

            for (size_t i = 0; i < size_; i++) {
                if (vec[i]) set_bit(bits_, i);
            }

        }

        ~BitVector() {
            delete[] bits_;
        }

        BitVector(const BitVector& other) { // copy constructor

            size_ = other.size_;
            alloc_size_ = other.alloc_size_;

            bits_ = new S[alloc_size_];
            std::copy(other.bits_, other.bits_ + alloc_size_, bits_);

        }

        BitVector(BitVector&& other) noexcept { // move constructor

            size_ = other.size_;
            alloc_size_ = other.alloc_size_;

            bits_ = other.bits_; other.bits_ = nullptr;

        }

        BitVector& operator=(const BitVector& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] bits_;

            // copy new resources
            size_ = other.size_;
            alloc_size_ = other.alloc_size_;

            bits_ = new S[alloc_size_];
            std::copy(other.bits_, other.bits_ + alloc_size_, bits_);

            return *this;

        }

        BitVector& operator=(BitVector&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] bits_;

            // copy / transfer new resources
            size_ = other.size_;
            alloc_size_ = other.alloc_size_;

            bits_ = other.bits_; other.bits_ = nullptr;

            return *this;

        }

        BitVector* clone() const { // deep-copy clone method
            return new BitVector<S>(*this);
        }

        // returns the length of the bit vector in bits
        size_t size() const {
            return size_;
        }

        // checks whether the length is zero
        bool empty() const {
            return size_ == 0;
        }

        // sets the i-th bit to 'val'
        void set(size_t i, bool val) {

            if (val) {
                set_bit(bits_, i);
            } else {
                reset_bit(bits_, i);
            }

        }

        // returns the value of the i-th bit
        bool get(size_t i) const {
            return get_bit(bits_, i);
        }

        // returns the value of the i-th bit
        bool operator[](size_t i) const {
            return get_bit(bits_, i);
        }

        // appends the bits of the provided bit vector at the end
        void append(const std::vector<bool>& bits) {

            size_t new_size = size_ + bits.size();
            size_t new_alloc_size = (new_size + (sizeof(S) * 8) - 1) / (sizeof(S) * 8);
            if (new_alloc_size > alloc_size_) {

                auto tmp = new S[new_alloc_size];
                std::copy(bits_, bits_ + alloc_size_, tmp);
                delete[] bits_;
                bits_ = tmp;
                alloc_size_ = new_alloc_size;

            }

            for (auto b : bits) {

                if (b) {
                    set_bit(bits_, size_);
                } else {
                    reset_bit(bits_, size_);
                }
                size_++;

            }


        }

        // overwrites the bits of the bit vector (starting at position i) with the bits
        // of the provided bit vector (does not check whether the bounds)
        size_t set(size_t i, const std::vector<bool>& bits) {

            for (auto b : bits) {

                if (b) {
                    set_bit(bits_, i);
                } else {
                    reset_bit(bits_, i);
                }
                i++;

            }

            return i;

        }

        // checks whether all bits with indices i to i + len - 1 are false
        bool is_zero(size_t i, size_t len) {

            size_t offset = i / (sizeof(S) * 8);
            size_t rel_i = i % (sizeof(S) * 8);

            return (rel_i == 0) ? (count_bits<S>(bits_ + offset, len) == 0) : ((count_bits<S>(bits_ + offset, rel_i + len) - count_bits<S>(bits_ + offset, rel_i)) == 0);

        }

        // determines memory consumption (in bytes)
        size_t size_in_bytes() const {
            return sizeof(BitVector<S>) // BitVector itself
                + sizeof(S) * alloc_size_; // bits_
        }


        S* get_bit_array() const {
            return bits_;
        }

    private:
        S* bits_; // actual bit vector
        size_t size_; // size of the represented bit vector in bits
        size_t alloc_size_; // size of the allocated array used to store the bits

    };


    /*
     * Rank data structure.
     *
     * The bit vector consists of blocks of sizeof(S) * 8 bits, which are grouped into superblocks
     * of 'FACTOR' blocks each.
     *
     * The number of ones in a block is computed using a precomputed popcount table.
     * In addition, we store the number of 1-bits before each superblock.
     *
     * Based on 'bitRankW32Int' from the original DACs implementation:
     * DACs: Bringing Direct Access to Variable-Length Codes (Brisaboa et al., 2013)
     * http://lbd.udc.es/research/DACS/
     */
    template<typename S, usmall_t FACTOR>
    class BitRank {

    public:
        BitRank() {

            data_ = nullptr;
            owner_ = false;
            length_ = 0;
            factor_ = 0;
            s_ = 0;
            super_blocks_ = nullptr;

        }

        /*
         * Construct rank data structure for the given bit vector.
         * The length is given in bits and the flag indicates
         * whether the rank data structure takes ownership
         * of the bit vector.
         */
        BitRank(S* bit_array, ularge_t length, bool owner) {

            data_ = bit_array;
            length_ = length;
            owner_ = owner;

            factor_ = (FACTOR == 0) ? bits(length_ - 1) : FACTOR;
            s_ = (sizeof(S) * 8) * factor_;

            // precompute the total number of ones before each superblock
            ularge_t num_super = (length_ + s_ - 1) / s_;
            super_blocks_ = new S[num_super];

            super_blocks_[0] = 0;
            for (ularge_t j = 1; j < num_super; j++) {

                super_blocks_[j] = super_blocks_[j - 1];
                super_blocks_[j] += count_in_interval((j - 1) * factor_, factor_);

            }

        }
        BitRank(const BitVector<S>& bit_vector) {

            data_ = bit_vector.get_bit_array();
            length_ = bit_vector.size();
            owner_ = false;

            factor_ = (FACTOR == 0) ? bits(length_ - 1) : FACTOR;
            s_ = (sizeof(S) * 8) * factor_;

            // precompute the total number of ones before each superblock
            ularge_t num_super = (length_ + s_ - 1) / s_;
            super_blocks_ = new S[num_super];

            super_blocks_[0] = 0;
            for (ularge_t j = 1; j < num_super; j++) {

                super_blocks_[j] = super_blocks_[j - 1];
                super_blocks_[j] += count_in_interval((j - 1) * factor_, factor_);

            }

        }

        ~BitRank() {

            delete[] super_blocks_;
            if (owner_) delete[] data_;

        }

        BitRank(const BitRank& other) { // copy constructor

            length_ = other.length_;
            owner_ = other.owner_;
            factor_ = other.factor_;
            s_ = other.s_;

            if (owner_) {

                data_ = new S[(length_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                std::copy(other.data_, other.data_ + (length_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), data_);

            } else {
                data_ = other.data_;
            }

            super_blocks_ = new S[(length_ + s_ - 1) / s_];
            std::copy(other.super_blocks_, other.super_blocks_ + (length_ + s_ - 1) / s_, super_blocks_);

        }

        BitRank(BitRank&& other) noexcept { // move constructor

            length_ = other.length_;
            owner_ = other.owner_;
            factor_ = other.factor_;
            s_ = other.s_;

            data_ = other.data_; other.data_ = nullptr;
            super_blocks_ = other.super_blocks_; other.super_blocks_ = nullptr;

        }

        BitRank& operator=(const BitRank& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            if (owner_) delete[] data_;
            delete[] super_blocks_;

            // copy new resources
            length_ = other.length_;
            owner_ = other.owner_;
            factor_ = other.factor_;
            s_ = other.s_;

            if (owner_) {

                data_ = new S[(length_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                std::copy(other.data_, other.data_ + (length_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), data_);

            } else {
                data_ = other.data_;
            }

            super_blocks_ = new S[(length_ + s_ - 1) / s_];
            std::copy(other.super_blocks_, other.super_blocks_ + (length_ + s_ - 1) / s_, super_blocks_);

            return *this;

        }

        BitRank& operator=(BitRank&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            if (owner_) delete[] data_;
            delete[] super_blocks_;

            // copy / transfer new resources
            length_ = other.length_;
            owner_ = other.owner_;
            factor_ = other.factor_;
            s_ = other.s_;

            data_ = other.data_; other.data_ = nullptr;
            super_blocks_ = other.super_blocks_; other.super_blocks_ = nullptr;

            return *this;

        }

        BitRank* clone() const { // deep-copy clone method
            return new BitRank<S, FACTOR>(*this);
        }

        /*
         * Determine the number of 1-bits up to and including index i (zero-based)
         *  = number of 1-bits before the superblock containing i
         *    + number of 1-bits in complete blocks before the block containing i (within the superblock containing i)
         *    + number of 1-bits in the block containing i up to and including index i
         */
        ularge_t rank(ularge_t i) const {

            ularge_t res = super_blocks_[i / s_];
            for (ularge_t a = (i / s_) * factor_; a < i / (sizeof(S) * 8); a++) {
                res += PopCounter<sizeof(S) * 8>::popcount(data_[a]);
            }
            res += PopCounter<sizeof(S) * 8>::popcount(data_[i / (sizeof(S) * 8)] & ~((static_cast<S>(1) << ((sizeof(S) * 8) - 1 - (i & ((sizeof(S) * 8) - 1)))) - 1));

            return res;

        }

        /*
         * Check whether the i-th bit is 1.
         */
        inline bool is_set(ularge_t i) const {
            return get_bit<S>(data_, i);
        }

        /*
         * Return the length of the underlying bit vector.
         */
        inline ularge_t size() const {
            return length_;
        }

        /*
         * Determine the memory consumption (in bytes).
         */
        ularge_t size_in_bytes() const {
            return sizeof(BitRank<S, FACTOR>) // BitRank itself
                + sizeof(S) * ((length_ + s_ - 1) / s_) // super_blocks_
                + (owner_ ? (sizeof(S) * ((length_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8))) : 0); // data_ (if owner)
        }


    private:
        S* data_; // pointer to the bit vector on which rank is computed
        bool owner_; // flag indicating whether the instance owns data_
        ularge_t length_; // length (in bits) of the underlying bit vector data_

        usmall_t factor_; // number of blocks per superblock (allows time / space trade-off)
        usmall_t s_; // length (in bits) of a superblock

        S* super_blocks_; // superblock array (stores total number of 1s before each superblock)


        /*
         * Return the number of 1s in the interval of blocks
         */
        ularge_t count_in_interval(ularge_t first_block, ularge_t num_blocks) const {

            ularge_t rank = 0;

            for (ularge_t i = first_block; (i < first_block + num_blocks); i++) {
                rank += PopCounter<sizeof(S) * 8>::popcount(data_[i]);
            }

            return rank;

        }

    };


    /*
     * Rank data structure.
     *
     * Similar to the one above but for cases assuming 1-based indexing.
     *
     * The bit vector consists of blocks of sizeof(S) * 8 bits, which are grouped into superblocks
     * of 'FACTOR' blocks each.
     *
     * The number of ones in a block is computed using a precomputed popcount table.
     * In addition, we store the number of 1-bits before each superblock.
     *
     * Based on 'bitRankW32Int' from the original DACs implementation:
     * DACs: Bringing Direct Access to Variable-Length Codes (Brisaboa et al., 2013)
     * http://lbd.udc.es/research/DACS/
     */
    template<typename S, usmall_t FACTOR>
    class BitRankOne {

    public:
        BitRankOne() {

            data_ = nullptr;
            perm_ = false;
            length_ = 0;
            factor_ = 0;
            s_ = 0;
            super_blocks_ = nullptr;

        }

        /*
         * Construct rank data structure for the given bit vector.
         * The length is given in bits and the flag indicates
         * whether a copy of the bit vector should be created for permanent internal use.
         */
        BitRankOne(S* bit_array, ularge_t length, bool perm) {

            length_ = length;
            perm_ = perm;
            if (perm_) {

                data_ = new S[(length_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                std::copy(bit_array, bit_array + (length_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), data_);

            } else {
                data_ = bit_array;
            }

            factor_ = (FACTOR == 0) ? bits(length_ - 1) : FACTOR;
            s_ = (sizeof(S) * 8) * factor_;

            // precompute the total number of ones before each superblock
            ularge_t num_super = std::max(1ul, (length_ + s_ - 1) / s_);
            super_blocks_ = new S[num_super];

            super_blocks_[0] = 0;
            for (ularge_t j = 1; j < num_super; j++) {

                super_blocks_[j] = super_blocks_[j - 1];
                super_blocks_[j] += count_in_interval((j - 1) * factor_, factor_);

            }

        }
        BitRankOne(const BitVector<S>& bit_vector, bool perm) {

            length_ = bit_vector.size();
            perm_ = perm;
            if (perm) {

                auto tmp = bit_vector.get_bit_array();
                data_ = new S[(length_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                std::copy(tmp, tmp + (length_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), data_);

            } else {
                data_ = bit_vector.get_bit_array();
            }

            factor_ = (FACTOR == 0) ? bits(length_ - 1) : FACTOR;
            s_ = (sizeof(S) * 8) * factor_;

            // precompute the total number of ones before each superblock
            ularge_t num_super = std::max(1ul, (length_ + s_ - 1) / s_);
            super_blocks_ = new S[num_super];

            super_blocks_[0] = 0;
            for (ularge_t j = 1; j < num_super; j++) {

                super_blocks_[j] = super_blocks_[j - 1];
                super_blocks_[j] += count_in_interval((j - 1) * factor_, factor_);

            }

        }

        ~BitRankOne() {

            delete[] super_blocks_;
            if (perm_) delete[] data_;

        }

        BitRankOne(const BitRankOne& other) { // copy constructor

            length_ = other.length_;
            perm_ = other.perm_;
            factor_ = other.factor_;
            s_ = other.s_;

            if (perm_) {

                data_ = new S[(length_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                std::copy(other.data_, other.data_ + (length_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), data_);

            } else {
                data_ = other.data_;
            }

            super_blocks_ = new S[std::max(1ul, (length_ + s_ - 1) / s_)];
            std::copy(other.super_blocks_, other.super_blocks_ + std::max(1ul, (length_ + s_ - 1) / s_), super_blocks_);

        }

        BitRankOne(BitRankOne&& other) noexcept { // move constructor

            length_ = other.length_; other.length_ = 0;
            perm_ = other.perm_; other.perm_ = false;
            factor_ = other.factor_;
            s_ = other.s_;

            data_ = other.data_; other.data_ = nullptr;
            super_blocks_ = other.super_blocks_; other.super_blocks_ = nullptr;

        }

        BitRankOne& operator=(const BitRankOne& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            if (perm_) delete[] data_;
            delete[] super_blocks_;

            // copy new resources
            length_ = other.length_;
            perm_ = other.perm_;
            factor_ = other.factor_;
            s_ = other.s_;

            if (perm_) {

                data_ = new S[(length_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8)];
                std::copy(other.data_, other.data_ + (length_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8), data_);

            } else {
                data_ = other.data_;
            }

            super_blocks_ = new S[std::max(1ul, (length_ + s_ - 1) / s_)];
            std::copy(other.super_blocks_, other.super_blocks_ + std::max(1ul, (length_ + s_ - 1) / s_), super_blocks_);

            return *this;

        }

        BitRankOne& operator=(BitRankOne&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            if (perm_) delete[] data_;
            delete[] super_blocks_;

            // copy / transfer new resources
            length_ = other.length_; other.length_ = 0;
            perm_ = other.perm_; other.perm_ = false;
            factor_ = other.factor_;
            s_ = other.s_;

            data_ = other.data_; other.data_ = nullptr;
            super_blocks_ = other.super_blocks_; other.super_blocks_ = nullptr;

            return *this;

        }

        BitRankOne* clone() const { // deep-copy clone method
            return new BitRankOne<S, FACTOR>(*this);
        }

        /*
         * Determine the number of 1-bits up to and including index i (zero-based)
         *  = number of 1-bits before the superblock containing i
         *    + number of 1-bits in complete blocks before the block containing i (within the superblock containing i)
         *    + number of 1-bits in the block containing i up to and including index i
         */
        ularge_t rank(ularge_t i) const {

            i--;
            ularge_t res = super_blocks_[i / s_];
            for (ularge_t a = (i / s_) * factor_; a < i / (sizeof(S) * 8); a++) {
                res += PopCounter<sizeof(S) * 8>::popcount(data_[a]);
            }
            res += PopCounter<sizeof(S) * 8>::popcount(data_[i / (sizeof(S) * 8)] & ~((static_cast<S>(1) << ((sizeof(S) * 8) - 1 - (i & ((sizeof(S) * 8) - 1)))) - 1));

            return res;

        }

        /*
         * Check whether the i-th bit is 1.
         */
        inline bool is_set(ularge_t i) const {
            return get_bit<S>(data_, i);
        }

        /*
         * Return the length of the underlying bit vector.
         */
        inline ularge_t size() const {
            return length_;
        }

        /*
         * Determine the memory consumption (in bytes).
         */
        ularge_t size_in_bytes() const {
            return sizeof(BitRankOne<S, FACTOR>) // BitRankOne itself
                + sizeof(S) * ((length_ + s_ - 1) / s_) // super_blocks_
                + (perm_ ? (sizeof(S) * ((length_ + (sizeof(S) * 8) - 1) / (sizeof(S) * 8))) : 0); // data_ (if owner)
        }


    private:
        S* data_; // pointer to the bit vector on which rank is computed
        bool perm_; // flag indicating whether data_ a permanent copy of the input bit vector
        ularge_t length_; // length (in bits) of the underlying bit vector data_

        usmall_t factor_; // number of blocks per superblock (allows time / space trade-off)
        usmall_t s_; // length (in bits) of a superblock

        S* super_blocks_; // superblock array (stores total number of 1s before each superblock)




        /*
         * Return the number of 1s in the interval of blocks
         */
        ularge_t count_in_interval(ularge_t first_block, ularge_t num_blocks) const {

            ularge_t rank = 0;

            for (ularge_t i = first_block; (i < first_block + num_blocks); i++) {
                rank += PopCounter<sizeof(S) * 8>::popcount(data_[i]);
            }

            return rank;

        }

    };

}

#endif //DACS_UTILITY_HPP
