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

#ifndef GEFAST_BASICS_HPP
#define GEFAST_BASICS_HPP

#include <iostream>
#include <sstream>
#include <string>

#include "../Base.hpp"

namespace GeFaST {

    /*
     * Abstract class for collectively storing strings contained in an AmpliconCollection.
     *
     * Derived classes are used by FlexibleAmpliconCollection to represent headers and sequences.
     */
    struct StringComponent {

        virtual ~StringComponent() = default;

        virtual StringComponent* clone() const = 0;

        /*
         * Add another string to the end of the collection.
         */
        virtual void push_back(const std::string& s) = 0;

        /*
         * Get access to the string stored at position i through a SequenceWrapper.
         *
         * Derived classes are not required to store (or to be able to determine) the lengths of the strings.
         * Thus, the length of the sequence should be provided as well.
         */
        virtual std::unique_ptr<SequenceWrapper> get(const numSeqs_t i, const lenSeqs_t len) const = 0;

        /*
         * Retrieve the string stored at position i.
         *
         * Derived classes are not required to store (or to be able to determine) the lengths of the strings.
         * Thus, the length of the sequence should be provided as well.
         */
        virtual std::string get_str(const numSeqs_t i, const lenSeqs_t len) const = 0;

        /*
         * Check whether collection can store the provided string and, if necessary, allocate additional space.
         */
        virtual void extend(const std::string& s) = 0;

        /*
         * Filling the collection with strings and retrieving them might be two separate phases
         * for the derived classes. Once all strings have been added, finalise() should be called once.
         */
        virtual void finalise() = 0;

        /*
         * Rearrange the inserted strings based on the provided permutation.
         * The amplicons in an AmpliconCollection are sorted based on their abundance.
         * In the FlexibleAmpliconCollection, the different information are stored in separate containers,
         * which have to be rearranged individually using the same permutation.
         */
        virtual void sort(const std::vector<numSeqs_t>& permutation) = 0;

        /*
         * Determine the memory consumption (in bytes).
         */
        virtual size_t size_in_bytes() const = 0;

        /*
         * Show memory profile.
         */
        virtual void show_memory() const = 0;

    };

    /*
     * Abstract class for collectively storing numeric values contained in an AmpliconCollection.
     *
     * Derived classes are used by FlexibleAmpliconCollection to represent lengths and abundances.
     */
    template<typename V>
    struct NumericComponent {

        virtual ~NumericComponent() = default;

        virtual NumericComponent* clone() const = 0;

        /*
         * Add another value to the end of the collection.
         */
        virtual void push_back(const V v) = 0;

        /*
         * Retrieve the value stored at position i.
         */
        virtual V get(const numSeqs_t i) const = 0;

        /*
         * Check whether collection can store another value and, if necessary, allocate additional space.
         */
        virtual void extend() = 0;

        /*
         * Filling the collection with values and retrieving them might be two separate phases
         * for the derived classes. Once all values have been added, finalise() should be called once.
         */
        virtual void finalise() = 0;

        /*
         * Rearrange the inserted values based on the provided permutation.
         * The amplicons in an AmpliconCollection are sorted based on their abundance.
         * In the FlexibleAmpliconCollection, the different information are stored in separate containers,
         * which have to be rearranged individually using the same permutation.
         */
        virtual void sort(const std::vector<numSeqs_t>& permutation) = 0;

        /*
         * Determine the memory consumption (in bytes).
         */
        virtual size_t size_in_bytes() const = 0;

        /*
         * Show memory profile.
         */
        virtual void show_memory() const = 0;

    };

    /*
     * Abstract class for collectively storing q-gram profiles contained in an AmpliconCollection.
     *
     * Derived classes are used by FlexibleAmpliconCollection to represent the information underlying the q-gram filter.
     */
    struct QgramComponent {

        virtual ~QgramComponent() = default;

        virtual QgramComponent* clone() const = 0;

        /*
         * Determine the q-gram profiles for the provided string and add the information to the end of the collection.
         */
        virtual void push_back(const std::string& seq) = 0;

        /*
         * Get access to the q-gram profiles stored at position i.
         */
        virtual unsigned char* get(const numSeqs_t i) const = 0;

        /*
         * Check whether collection can store another q-gram profiles and, if necessary, allocate additional space.
         */
        virtual void extend() = 0;

        /*
         * Filling the collection with values and retrieving them might be two separate phases
         * for the derived classes. Once all values have been added, finalise() should be called once.
         */
        virtual void finalise() = 0;

        /*
         * Rearrange the inserted values based on the provided permutation.
         * The amplicons in an AmpliconCollection are sorted based on their abundance.
         * In the FlexibleAmpliconCollection, the different information are stored in separate containers,
         * which have to be rearranged individually using the same permutation.
         */
        virtual void sort(const std::vector<numSeqs_t>& permutation) = 0;

        /*
         * Determine the memory consumption (in bytes).
         */
        virtual size_t size_in_bytes() const = 0;

        /*
         * Show memory profile.
         */
        virtual void show_memory() const = 0;

    };



    /* === Approximation of memory consumption of STL containers === */

    /*
     * Approximate the memory consumption (in bytes) of an std::vector of POD elements.
     */
    template<typename T>
    size_t vec_pod_size_in_bytes(const std::vector<T>& vec) {
        return sizeof(std::vector<T>) // vector itself
            + sizeof(T) * vec.capacity(); // memory of elements
    }

    /*
     * Approximate the memory consumption (in bytes) of an std::vector of POD elements.
     */
    size_t vec_bool_size_in_bytes(const std::vector<bool>& vec);

    /*
     * Approximate the memory consumption (in bytes) of the std::vector of std::string instances.
     * Assumes that small-string optimisation does not occur (as, e.g., with GCC 4.9.2).
     */
    size_t vec_str_size_in_bytes(const std::vector<std::string>& vec);

    /*
     * Approximate the memory consumption (in bytes) of an std::map<std::string, unsigned char> instance.
     * Assumes the _Rb_tree implementation of GCC 4.9.2.
     * Assumes that small-string optimisation does not occur (as, e.g., with GCC 4.9.2).
     *
     * Inspired by https://stackoverflow.com/questions/720507/how-can-i-estimate-memory-usage-of-stdmap/720520.
     */
    size_t map_str_uchar_size_in_bytes(const std::map<std::string, unsigned char>& map);

    /*
     * Approximate the memory consumption (in bytes) of an std::map<K, V> instance with both K and V being POD types.
     * Assumes the _Rb_tree implementation of GCC 4.9.2.
     *
     * Inspired by https://stackoverflow.com/questions/720507/how-can-i-estimate-memory-usage-of-stdmap/720520.
     */
    template<typename K, typename V>
    size_t map_pod_size_in_bytes(const std::map<K, V>& map) {

        // internal typedefs of the map implementation (red-black trees)
        typedef K key_type;
        typedef V mapped_type;
        typedef std::pair<const key_type, mapped_type> value_type;
        typedef std::less<key_type> key_compare;
        typedef std::allocator<std::pair<const key_type, mapped_type>> allocator_type;
        typedef typename __gnu_cxx::__alloc_traits<allocator_type>::template rebind<value_type>::other pair_alloc_type;
        typedef std::_Rb_tree<key_type, value_type, std::_Select1st<value_type>, key_compare, pair_alloc_type> rep_type;

        return sizeof(rep_type) // std::map contains only member of type rep_type
            + sizeof(std::_Rb_tree_node<value_type>) * map.size(); // tree nodes (colour, pointers, key-value pair)

    }

    /*
     * Approximate the memory consumption (in bytes) of an std::unordered_map<K, V> instance with both K and V being POD types.
     * Assumes the _Hashtable implementation of GCC 4.9.2.
     */
    template<typename K, typename V>
    size_t umap_pod_size_in_bytes(const std::unordered_map<K, V>& umap) {

        size_t allocated_nodes = 0;
        for (auto i = 0; i < umap.bucket_count(); i++) {
            allocated_nodes += umap.bucket_size(i);
        }

        typedef K key_type;
        typedef V value_type;
        typedef std::hash<key_type> hash_type;
        typedef std::equal_to<key_type> pred_type;
        typedef std::allocator<std::pair<const key_type, value_type>> alloc_type;
        typedef std::__umap_traits<std::__cache_default<key_type, hash_type>::value> traits_type;
        typedef std::_Hashtable<key_type, std::pair<const key_type, value_type>, alloc_type,
            std::__detail::_Select1st, pred_type, hash_type, std::__detail::_Mod_range_hashing,
            std::__detail::_Default_ranged_hash, std::__detail::_Prime_rehash_policy,
            traits_type> hash_table_type;
        typedef std::__detail::_Hash_node<std::pair<const key_type, value_type>, traits_type::__hash_cached::value> node_type;

        return sizeof(hash_table_type) // std::unordered_map contains only member of type hash_table_type
            + sizeof(node_type*) * umap.bucket_count() // buckets (pointer to node before first element)
            + sizeof(node_type) * allocated_nodes; // nodes (at least next pointer and value)

    }

    /*
     * Approximate the memory consumption (in bytes) of a sip_map_type instance as defined in (Static)SpaceSegmentRelations.
     * Assumes the _Hashtable implementation of GCC 4.9.2.
     */
    size_t sip_map_single_size_in_bytes(std::unordered_map<std::unique_ptr<SubstringWrapper>, numSeqs_t, hashSubstringWrapper, equalSubstringWrapper>& umap);

    /*
     * Approximate the memory consumption (in bytes) of a sip_map_type instance as defined in SegmentRelations.
     * Assumes the _Hashtable implementation of GCC 4.9.2.
     */
    size_t sip_map_multiple_size_in_bytes(std::unordered_map<std::unique_ptr<SubstringWrapper>, std::vector<numSeqs_t>, hashSubstringWrapper, equalSubstringWrapper>& umap);

    /*
     * Approximate the memory consumption (in bytes) of an std::set of POD elements.
     * Assumes the _Rb_tree implementation of GCC 4.9.2.
     *
     * Inspired by https://stackoverflow.com/questions/720507/how-can-i-estimate-memory-usage-of-stdmap/720520.
     */
    template<typename T>
    size_t set_pod_size_in_bytes(const std::set<T>& set) {

        // internal typedefs of the set implementation (red-black trees)
        typedef T key_type;
        typedef T value_type;
        typedef std::less<key_type> key_compare;
        typedef std::less<key_type> value_compare;
        typedef std::allocator<key_type> allocator_type;
        typedef typename __gnu_cxx::__alloc_traits<allocator_type>::template rebind<key_type>::other key_alloc_type;
        typedef std::_Rb_tree<key_type, value_type, std::_Identity<value_type>, key_compare, key_alloc_type> rep_type;

        return sizeof(rep_type) // std::set contains only member of type rep_type
            + sizeof(std::_Rb_tree_node<value_type>) * set.size(); // tree nodes (colour, pointers, value)

    }

    /*
     * Approximate the memory consumption (in bytes) of an std::set<std::unique_ptr<SequenceWrapper>, lessSequenceWrapper> instance
     * as used by auxiliary data structures.
     * Assumes the _Rb_tree implementation of GCC 4.9.2.
     *
     * Inspired by https://stackoverflow.com/questions/720507/how-can-i-estimate-memory-usage-of-stdmap/720520.
     */
    size_t different_seqs_size_in_bytes(const std::set<std::unique_ptr<SequenceWrapper>, lessSequenceWrapper>& set);

    /*
     * Approximate the memory consumption (in bytes) of an std::map<std::pair<lenSeqs_t, lenSeqs_t>, Substrings*> instance
     * as used by auxiliary data structures.
     * Assumes the _Rb_tree implementation of GCC 4.9.2.
     *
     * Inspired by https://stackoverflow.com/questions/720507/how-can-i-estimate-memory-usage-of-stdmap/720520.
     */
    // size_t substrs_archive_size_in_bytes(const std::map<std::pair<lenSeqs_t, lenSeqs_t>, Substrings*>& map);
    // in AuxiliaryData.hpp to avoid #include problems

}

#endif //GEFAST_BASICS_HPP
