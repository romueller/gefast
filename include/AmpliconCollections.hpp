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

#ifndef GEFAST_AMPLICONCOLLECTIONS_HPP
#define GEFAST_AMPLICONCOLLECTIONS_HPP

#include <algorithm>
#include <cstring>
#include <fstream>
#include <map>
#include <numeric>
#include <sstream>
#include <vector>

#include "Base.hpp"
#include "Utility.hpp"

namespace GeFaST {

    /*
     * Simple version of an amplicon collection.
     *
     * Identifiers and strings are stored as STL strings in two STL vectors.
     * Abundance values are stored in a third STL vector.
     */
    class SimpleAmpliconCollection : public AmpliconCollection {

    public:
        SimpleAmpliconCollection* clone() const override; // deep-copy clone method


        void emplace_back(const Defline& dl, const std::string& seq, const ExtraInfo& extra) override;

        numSeqs_t size() const override;

        inline const char* id(const numSeqs_t i) const override;

        inline const char* seq(const numSeqs_t i) const override;

        inline lenSeqs_t len(const numSeqs_t i) const override;

        inline numSeqs_t ab(const numSeqs_t i) const override;

        inline const char* quals(const numSeqs_t i) const override;

        /*
         * No q-gram profiles are stored.
         * Return 0 in order to avoid falsely rejecting the amplicon pair when checking the q-gram difference.
         */
        unsigned long qgram_diff(const numSeqs_t i, const numSeqs_t j) const override;

        /*
         * Sort the amplicons by abundance (decreasing) and
         * use the lexicographical rank (ascending) of the identifiers as the tie-breaker.
         */
        void sort();

        /*
         * Append the amplicons in the stored order in FASTA format to the specified file.
         */
        void print(const std::string& file, const Configuration& config) const;

    protected:
        std::vector<std::string> ids_; // amplicon identifiers
        std::vector<std::string> seqs_; // amplicon sequences
        std::vector<numSeqs_t> abunds_; // amplicon abundances

    };


    /*
     * More (space-)efficient version of an amplicon collection.
     *
     * Identifiers and sequences are stored as C-strings in two large C-style arrays.
     * The individual strings are accessed via stored pointers into the C-style arrays.
     * Lengths and abundance values are also stored in C-style arrays.
     */
    class ArrayAmpliconCollection : public AmpliconCollection {

    public:
        ArrayAmpliconCollection(const numSeqs_t capacity, const unsigned long long capacity_headers,
                const unsigned long long capacity_sequences);

        virtual ~ArrayAmpliconCollection();

        ArrayAmpliconCollection(const ArrayAmpliconCollection& other); // copy constructor

        ArrayAmpliconCollection(ArrayAmpliconCollection&& other) noexcept; // move constructor

        ArrayAmpliconCollection& operator=(const ArrayAmpliconCollection& other); // copy assignment operator

        ArrayAmpliconCollection& operator=(ArrayAmpliconCollection&& other) noexcept; // move assignment operator

        ArrayAmpliconCollection* clone() const override; // deep-copy clone method


        void emplace_back(const Defline& dl, const std::string& seq, const ExtraInfo& extra) override;

        numSeqs_t size() const override;

        inline const char* id(const numSeqs_t i) const override;

        inline const char* seq(const numSeqs_t i) const override;

        inline lenSeqs_t len(const numSeqs_t i) const override;

        inline numSeqs_t ab(const numSeqs_t i) const override;

        inline const char* quals(const numSeqs_t i) const override;

        /*
         * No q-gram profiles are stored.
         * Return 0 in order to avoid falsely rejecting the amplicon pair when checking the q-gram difference.
         */
        unsigned long qgram_diff(const numSeqs_t i, const numSeqs_t j) const override;

        /*
         * Sort the amplicons by abundance (decreasing) and
         * use the lexicographical rank (ascending) of the identifiers as the tie-breaker.
         */
        void sort();

        /*
         * Append amplicon to the collection, extending its capacity if necessary.
         *
         * Ideally, the amplicon collection is initialised with sufficient space (see the constructor)
         * but if the number of amplicons gets too large or the space for the identifiers resp. sequences
         * becomes insufficient, the capacity of the array is doubled.
         */
        void extend_arrays(const Defline& dl, const std::string& seq);

        /*
         * Append the amplicons in the stored order in FASTA format to the specified file.
         */
        void print(const std::string& file, const Configuration& config) const;

    protected:
        ArrayAmpliconCollection() = default;

        char* headers_; // overall identifiers array, each string ends with a \0
        char* next_header_; // position at which the next identifier would be inserted
        unsigned long long capacity_headers_; // capacity of headers_

        char* sequences_; // overall sequences array, each string ends with a \0
        char* next_sequence_; // position at which the next sequence would be inserted
        unsigned long long capacity_sequences_; // capacity of sequences_

        char** header_pointers_; // array of pointers into headers_ marking the start positions of the individual identifiers
        char** seq_pointers_; // array of pointers into sequences_ marking the start positions of the individual sequences
        numSeqs_t* abundances_; // array of abundances
        lenSeqs_t* lengths_; // array of lengths
        numSeqs_t size_; // current number of amplicons in the collection
        numSeqs_t capacity_; // capacity of the collection (in terms of the number of amplicons)

    };


    /*
     * Array-style version of an amplicon collection with quality scores.
     *
     * Identifiers, sequences and encoded quality scores are stored as C-strings in two large C-style arrays.
     * The individual strings are accessed via stored pointers into the C-style arrays.
     * The individual quality-score strings are accessed by reusing the offsets provided by the sequence pointers.
     * Lengths and abundance values are also stored in C-style arrays.
     */
    class ArrayQualityAmpliconCollection : public AmpliconCollection {

    public:
        ArrayQualityAmpliconCollection();

        ArrayQualityAmpliconCollection(const numSeqs_t capacity, const unsigned long long capacity_headers,
                const unsigned long long capacity_sequences);

        ~ArrayQualityAmpliconCollection();

        ArrayQualityAmpliconCollection(const ArrayQualityAmpliconCollection& other); // copy constructor

        ArrayQualityAmpliconCollection(ArrayQualityAmpliconCollection&& other) noexcept; // move constructor

        ArrayQualityAmpliconCollection& operator=(const ArrayQualityAmpliconCollection& other); // copy assignment operator

        ArrayQualityAmpliconCollection& operator=(ArrayQualityAmpliconCollection&& other) noexcept; // move assignment operator

        ArrayQualityAmpliconCollection* clone() const override; // deep-copy clone method


        void emplace_back(const Defline& dl, const std::string& seq, const ExtraInfo& extra) override;

        numSeqs_t size() const override;

        inline const char* id(const numSeqs_t i) const override;

        inline const char* seq(const numSeqs_t i) const override;

        inline lenSeqs_t len(const numSeqs_t i) const override;

        inline numSeqs_t ab(const numSeqs_t i) const override;

        inline const char* quals(const numSeqs_t i) const override;

        /*
         * Sort the amplicons by abundance (decreasing) and
         * use the lexicographical rank (ascending) of the identifiers as the tie-breaker.
         */
        void sort();

        /*
         * Append amplicon to the collection, extending its capacity if necessary.
         *
         * Ideally, the amplicon collection is initialised with sufficient space (see the constructor)
         * but if the number of amplicons gets too large or the space for the identifiers resp. sequences
         * becomes insufficient, the capacity of the array is doubled.
         */
        void extend_arrays(const Defline& dl, const std::string& seq);

        /*
         * Append the amplicons in the stored order in FASTQ format to the specified file.
         */
        void print(const std::string& file, const Configuration& config) const;

    protected:
        char* headers_; // overall identifiers array, each string ends with a \0
        char* next_header_; // position at which the next identifier would be inserted
        unsigned long long capacity_headers_; // capacity of headers_

        char* sequences_; // overall sequences array, each string ends with a \0
        char* next_sequence_; // position at which the next sequence would be inserted
        unsigned long long capacity_sequences_; // capacity of sequences_

        char* qualities_; // overall quality-scores array, each string ends with a \0
        char* next_qualities_; // position at which the next quality-score sequence would be inserted

        char** header_pointers_; // array of pointers into headers_ marking the start positions of the individual identifiers
        char** seq_pointers_; // array of pointers into sequences_ marking the start positions of the individual sequences
        numSeqs_t* abundances_; // array of abundances
        lenSeqs_t* lengths_; // array of lengths
        numSeqs_t size_; // current number of amplicons in the collection
        numSeqs_t capacity_; // capacity of the collection (in terms of the number of amplicons)

    };


    /*
     * More (space-)efficient version of an amplicon collection,
     * also storing a q-gram profile for each amplicon.
     *
     * Identifiers and sequences are stored as C-strings in two large C-style arrays.
     * The individual strings are accessed via stored pointers into the C-style arrays.
     * Lengths and abundance values are also stored in C-style arrays.
     */
    class ArrayQgramAmpliconCollection : public AmpliconCollection {

    public:
        ArrayQgramAmpliconCollection(const numSeqs_t capacity, const unsigned long long capacity_headers,
                const unsigned long long capacity_sequences);

        virtual ~ArrayQgramAmpliconCollection();

        ArrayQgramAmpliconCollection(const ArrayQgramAmpliconCollection& other); // copy constructor

        ArrayQgramAmpliconCollection(ArrayQgramAmpliconCollection&& other) noexcept; // move constructor

        ArrayQgramAmpliconCollection& operator=(const ArrayQgramAmpliconCollection& other); // copy assignment operator

        ArrayQgramAmpliconCollection& operator=(ArrayQgramAmpliconCollection&& other) noexcept; // move assignment operator

        ArrayQgramAmpliconCollection* clone() const override; // deep-copy clone method


        void emplace_back(const Defline& dl, const std::string& seq, const ExtraInfo& extra) override;

        numSeqs_t size() const override;

        inline const char* id(const numSeqs_t i) const override;

        inline const char* seq(const numSeqs_t i) const override;

        inline lenSeqs_t len(const numSeqs_t i) const override;

        inline numSeqs_t ab(const numSeqs_t i) const override;

        inline const char* quals(const numSeqs_t i) const override;

        /*
         * Compute lower bound on the q-gram distance between two amplicons.
         *
         * Based on the ideas from "Approximate string-matching with q-grams and maximal matches"
         * (Ukkonen, 1992, https://core.ac.uk/download/pdf/82613520.pdf),
         * as implemented in Swarm (https://github.com/torognes/swarm).
         *
         * Only considers whether a q-gram has an even (including 0) or odd number of occurrences,
         * instead of using actual counts of the q-grams.
         */
        unsigned long qgram_diff(const numSeqs_t i, const numSeqs_t j) const override;

        /*
         * Sort the amplicons by abundance (decreasing) and
         * use the lexicographical rank (ascending) of the identifiers as the tie-breaker.
         */
        void sort();

        /*
         * Append amplicon to the collection, extending its capacity if necessary.
         *
         * Ideally, the amplicon collection is initialised with sufficient space (see the constructor)
         * but if the number of amplicons gets too large or the space for the identifiers resp. sequences
         * becomes insufficient, the capacity of the array is doubled.
         */
        void extend_arrays(const Defline& dl, const std::string& seq);

        /*
         * Append the amplicons in the stored order in FASTA format to the specified file.
         */
        void print(const std::string& file, const Configuration& config) const;

    protected:
        ArrayQgramAmpliconCollection() = default;

        char* headers_; // overall identifiers array, each string ends with a \0
        char* next_header_; // position at which the next identifier would be inserted
        unsigned long long capacity_headers_; // capacity of headers_

        char* sequences_; // overall sequences array, each string ends with a \0
        char* next_sequence_; // position at which the next sequence would be inserted
        unsigned long long capacity_sequences_; // capacity of sequences_

        char** header_pointers_; // array of pointers into headers_ marking the start positions of the individual identifiers
        char** seq_pointers_; // array of pointers into sequences_ marking the start positions of the individual sequences
        numSeqs_t* abundances_; // array of abundances
        lenSeqs_t* lengths_; // array of lengths
        numSeqs_t size_; // current number of amplicons in the collection
        numSeqs_t capacity_; // capacity of the collection (in terms of the number of amplicons)

        unsigned char* qgrams_; // array of q-gram profiles
        numSeqs_t* qgram_offsets_; // offsets to start positions of the individual q-gram profiles in qgrams_

    };

    /*
     * More (space-)efficient version of an amplicon collection with quality scores,
     * also storing a q-gram profile for each amplicon.
     *
     * Identifiers and sequences are stored as C-strings in two large C-style arrays.
     * The individual strings are accessed via stored pointers into the C-style arrays.
     * The individual quality-score strings are accessed by reusing the offsets provided by the sequence pointers.
     * Lengths and abundance values are also stored in C-style arrays.
     */
    class ArrayQgramQualityAmpliconCollection : public AmpliconCollection {

    public:
        ArrayQgramQualityAmpliconCollection(const numSeqs_t capacity, const unsigned long long capacity_headers,
                const unsigned long long capacity_sequences);

        virtual ~ArrayQgramQualityAmpliconCollection();

        ArrayQgramQualityAmpliconCollection(const ArrayQgramQualityAmpliconCollection& other); // copy constructor

        ArrayQgramQualityAmpliconCollection(ArrayQgramQualityAmpliconCollection&& other) noexcept; // move constructor

        ArrayQgramQualityAmpliconCollection& operator=(const ArrayQgramQualityAmpliconCollection& other); // copy assignment operator

        ArrayQgramQualityAmpliconCollection& operator=(ArrayQgramQualityAmpliconCollection&& other) noexcept; // move assignment operator

        ArrayQgramQualityAmpliconCollection* clone() const override; // deep-copy clone method


        void emplace_back(const Defline& dl, const std::string& seq, const ExtraInfo& extra) override;

        numSeqs_t size() const override;

        inline const char* id(const numSeqs_t i) const override;

        inline const char* seq(const numSeqs_t i) const override;

        inline lenSeqs_t len(const numSeqs_t i) const override;

        inline numSeqs_t ab(const numSeqs_t i) const override;

        inline const char* quals(const numSeqs_t i) const override;

        /*
         * Sort the amplicons by abundance (decreasing) and
         * use the lexicographical rank (ascending) of the identifiers as the tie-breaker.
         */
        void sort();

        /*
         * Append amplicon to the collection, extending its capacity if necessary.
         *
         * Ideally, the amplicon collection is initialised with sufficient space (see the constructor)
         * but if the number of amplicons gets too large or the space for the identifiers resp. sequences
         * becomes insufficient, the capacity of the array is doubled.
         */
        void extend_arrays(const Defline& dl, const std::string& seq);

        /*
         * Append the amplicons in the stored order in FASTA format to the specified file.
         */
        void print(const std::string& file, const Configuration& config) const;

        /*
         * Compute lower bound on the q-gram distance between two amplicons.
         *
         * Based on the ideas from "Approximate string-matching with q-grams and maximal matches"
         * (Ukkonen, 1992, https://core.ac.uk/download/pdf/82613520.pdf),
         * as implemented in Swarm (https://github.com/torognes/swarm).
         *
         * Only considers whether a q-gram has an even (including 0) or odd number of occurrences,
         * instead of using actual counts of the q-grams.
         */
        unsigned long qgram_diff(const numSeqs_t i, const numSeqs_t j) const override;

    protected:
        ArrayQgramQualityAmpliconCollection() = default;

        char* headers_; // overall identifiers array, each string ends with a \0
        char* next_header_; // position at which the next identifier would be inserted
        unsigned long long capacity_headers_; // capacity of headers_

        char* sequences_; // overall sequences array, each string ends with a \0
        char* next_sequence_; // position at which the next sequence would be inserted
        unsigned long long capacity_sequences_; // capacity of sequences_

        char* qualities_; // overall strings (headers, sequences) array, each string ends with a \0
        char* next_qualities_; // position at which the next string would be inserted

        char** header_pointers_; // array of pointers into headers_ marking the start positions of the individual identifiers
        char** seq_pointers_; // array of pointers into sequences_ marking the start positions of the individual sequences
        numSeqs_t* abundances_; // array of abundances
        lenSeqs_t* lengths_; // array of lengths
        numSeqs_t size_; // current number of amplicons in the collection
        numSeqs_t capacity_; // capacity of the collection (in terms of the number of amplicons)

        unsigned char* qgrams_; // array of q-gram profiles
        numSeqs_t* qgram_offsets_; // offsets to start positions of the individual q-gram profiles in qgrams_

    };

}

#endif //GEFAST_AMPLICONCOLLECTIONS_HPP
