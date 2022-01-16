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

#include "../../include/space/Basics.hpp"

namespace GeFaST {

    size_t vec_bool_size_in_bytes(const std::vector<bool>& vec) {
        return sizeof(std::vector<bool>) // vector itself
            + (vec.capacity() + 7) / 8; // memory of elements
    }

    size_t vec_str_size_in_bytes(const std::vector<std::string>& vec) {

        size_t total_capa = 0; // total allocated capacity of all strings (including null terminators)
        for (auto i = 0; i < vec.size(); i++) {
            total_capa += vec[i].capacity() + 1; // capacity() does not account for null terminator
        }

        return sizeof(std::vector<std::string>) // vector itself
               + sizeof(std::string) * vec.capacity() // basic memory of the strings
               + sizeof(char) * total_capa; // actual string contents

    }

    size_t map_str_uchar_size_in_bytes(const std::map<std::string, unsigned char>& map) {

        // internal typedefs of the map implementation (red-black trees)
        typedef std::string key_type;
        typedef unsigned char mapped_type;
        typedef std::pair<const key_type, mapped_type> value_type;
        typedef std::less<key_type> key_compare;
        typedef std::allocator<std::pair<const key_type, mapped_type>> allocator_type;
        typedef typename __gnu_cxx::__alloc_traits<allocator_type>::template rebind<value_type>::other pair_alloc_type;
        typedef std::_Rb_tree<key_type, value_type, std::_Select1st<value_type>, key_compare, pair_alloc_type> rep_type;

        size_t total_capa = 0; // total allocated capacity of all strings (including null terminators)
        for (auto& kv : map) {
            total_capa += kv.first.capacity() + 1; // capacity() does not account for null terminator
        }

        return sizeof(rep_type) // std::map contains only member of type rep_type
               + sizeof(std::_Rb_tree_node<value_type>) * map.size() // tree nodes (colour, pointers, key-value pair)
               + sizeof(char) * total_capa; // actual string contents

    }

    size_t sip_map_single_size_in_bytes(std::unordered_map<std::unique_ptr<SubstringWrapper>, numSeqs_t, hashSubstringWrapper, equalSubstringWrapper>& umap) {

        size_t allocated_nodes = 0;
        for (auto i = 0; i < umap.bucket_count(); i++) {
            allocated_nodes += umap.bucket_size(i);
        }
        size_t allocated_wrappers = 0;
        for (auto& kv : umap) {
            allocated_wrappers += kv.first->size_in_bytes();
        }

        typedef std::unique_ptr<SubstringWrapper> key_type;
        typedef numSeqs_t value_type;
        typedef hashSubstringWrapper hash_type;
        typedef equalSubstringWrapper pred_type;
        typedef std::allocator<std::pair<const key_type, value_type>> alloc_type;
        typedef std::__umap_traits<std::__cache_default<key_type, hash_type>::value> traits_type;
        typedef std::_Hashtable<key_type, std::pair<const key_type, value_type>, alloc_type,
                std::__detail::_Select1st, pred_type, hash_type, std::__detail::_Mod_range_hashing,
                std::__detail::_Default_ranged_hash, std::__detail::_Prime_rehash_policy,
                traits_type> hash_table_type;
        typedef std::__detail::_Hash_node<std::pair<const key_type, value_type>, traits_type::__hash_cached::value> node_type;

        return sizeof(hash_table_type) // std::unordered_map contains only member of type hash_table_type
               + sizeof(node_type*) * umap.bucket_count() // buckets (pointer to node before first element)
               + sizeof(node_type) * allocated_nodes // nodes (at least next pointer and value)
               + allocated_wrappers; // memory of the SubstringWrappers

    }

    size_t sip_map_multiple_size_in_bytes(std::unordered_map<std::unique_ptr<SubstringWrapper>, std::vector<numSeqs_t>, hashSubstringWrapper, equalSubstringWrapper>& umap) {

        size_t allocated_nodes = 0;
        for (auto i = 0; i < umap.bucket_count(); i++) {
            allocated_nodes += umap.bucket_size(i);
        }
        size_t allocated_wrappers = 0;
        for (auto& kv : umap) {
            allocated_wrappers += kv.first->size_in_bytes();
        }
        size_t allocated_values = 0;
        for (auto& kv : umap) {
            allocated_values += vec_pod_size_in_bytes(kv.second) - sizeof(std::vector<numSeqs_t>);
        }

        typedef std::unique_ptr<SubstringWrapper> key_type;
        typedef std::vector<numSeqs_t> value_type;
        typedef hashSubstringWrapper hash_type;
        typedef equalSubstringWrapper pred_type;
        typedef std::allocator<std::pair<const key_type, value_type>> alloc_type;
        typedef std::__umap_traits<std::__cache_default<key_type, hash_type>::value> traits_type;
        typedef std::_Hashtable<key_type, std::pair<const key_type, value_type>, alloc_type,
                std::__detail::_Select1st, pred_type, hash_type, std::__detail::_Mod_range_hashing,
                std::__detail::_Default_ranged_hash, std::__detail::_Prime_rehash_policy,
                traits_type> hash_table_type;
        typedef std::__detail::_Hash_node<std::pair<const key_type, value_type>, traits_type::__hash_cached::value> node_type;

        return sizeof(hash_table_type) // std::unordered_map contains only member of type hash_table_type
               + sizeof(node_type*) * umap.bucket_count() // buckets (pointer to node before first element)
               + sizeof(node_type) * allocated_nodes // nodes (at least next pointer and value)
               + allocated_wrappers // memory of the SubstringWrappers
               + allocated_values; // memory of the vector contents

    }

    size_t different_seqs_size_in_bytes(const std::set<std::unique_ptr<SequenceWrapper>, lessSequenceWrapper>& set) {

        // internal typedefs of the set implementation (red-black trees)
        typedef std::unique_ptr<SequenceWrapper> key_type;
        typedef std::unique_ptr<SequenceWrapper> value_type;
        typedef lessSequenceWrapper key_compare;
        typedef lessSequenceWrapper value_compare;
        typedef std::allocator<key_type> allocator_type;
        typedef typename __gnu_cxx::__alloc_traits<allocator_type>::template rebind<key_type>::other key_alloc_type;
        typedef std::_Rb_tree<key_type, value_type, std::_Identity<value_type>, key_compare, key_alloc_type> rep_type;

        size_t allocated_wrappers = 0;
        for (auto& v : set) {
            allocated_wrappers += v->size_in_bytes();
        }

        return sizeof(rep_type) // std::set contains only member of type rep_type
            + sizeof(std::_Rb_tree_node<value_type>) * set.size() // tree nodes (colour, pointers, value)
            + allocated_wrappers; // memory of the SequenceWrapper

    }

}
