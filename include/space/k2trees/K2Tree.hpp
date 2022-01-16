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

#ifndef K2TREES_K2TREE_HPP
#define K2TREES_K2TREE_HPP

#include "Utility.hpp"

namespace GeFaST {

    /*
     * Representation of a (weighted / valued) relation.
     *
     * A K2Tree with n rows and m columns describes a relation R
     * between the sets [0 : n - 1] and [0 : m - 1].
     * The weights / values of an entry in the relation are of type V
     * and one value is designated the null value ("unrelated").
     * The sets (and, thus, the used indices) are of type I.
     *
     * The data structure is static (with the exception of the set_null() method).
     *
     * Originally adapted from:
     * Brisaboa, N. R., Ladra, S., & Navarro, G. (2014).
     * Compact representation of Web graphs with extended functionality.
     * Information Systems, 39(1), 152–174.
     * http://doi.org/10.1016/j.is.2013.08.003
     *
     * Refined class K2Tree<E> (https://github.com/romueller/k2trees/blob/master/K2Tree.hpp).
     */
    template<typename V, typename I = size_t>
    class K2Tree {

    public:
        // weights / values
        typedef V value_t;
        typedef I index_t;

        typedef std::vector<value_t> values_t;
        typedef std::vector<index_t> indices_t;

        // input type ((adjacency) matrix)
        typedef std::vector<std::vector<bool>> matrix_t;
        typedef std::vector<values_t> value_matrix_t;

        // input type (single (adjacency) list of pairs describing column and value)
        typedef std::vector<index_t> row_t;
        typedef std::vector<row_t> lists_t;

        typedef std::vector<ValueRowPosition<value_t, index_t>> value_row_t;
        typedef std::vector<value_row_t> value_lists_t;

        // output type (collection of positions (row and column number) of relation entries)
        // input / output type ((adjacency) pairs describing row, column and value)
        typedef std::vector<Position<index_t>> positions_t;
        typedef std::vector<ValuePosition<value_t, index_t>> value_positions_t;


        K2Tree() = default;
        virtual ~K2Tree() = default;

        K2Tree(const K2Tree& other) = default;
        K2Tree(K2Tree&& other) noexcept = default;
        K2Tree& operator=(const K2Tree& other) = default;
        K2Tree& operator=(K2Tree&& other) noexcept = default;

        // creates a deep copy
        virtual K2Tree* clone() const = 0;


        // returns the number of rows of the relation (n)
        virtual index_t get_num_rows() const = 0;

        // returns the number of columns of the relation (m)
        virtual index_t get_num_cols() const = 0;


        // checks whether (i,j) is in R
        virtual bool is_not_null(index_t i, index_t j) = 0;

        // sets the value of the pair (i,j) to null, i.e. removes it from the relation
        virtual void set_null(index_t i, index_t j, bool deep) = 0;
        virtual void set_null_deep(index_t i, index_t j) = 0;
        virtual void set_null_shallow(index_t i, index_t j) = 0;

        // sets the value of all pairs (.,j) to null, i.e. removes them from the relation
        virtual void set_null_column(index_t j, bool deep) = 0;
        virtual void set_null_column_deep(index_t j) = 0;
        virtual void set_null_column_shallow(index_t j) = 0;

        // returns the value of (i,j), if the pair is in R, null otherwise
        virtual value_t get_value(index_t i, index_t j) = 0;

        // returns the smallest column number j such that (i,j) is in R, or a value >= n if no such pairs exists
        virtual index_t get_first_successor(index_t i) = 0;

        // returns the values of all pairs in R whose first component is i
        virtual values_t get_successor_values(index_t i) = 0;

        // returns the column numbers of all pairs in R whose first component is i
        virtual indices_t get_successor_positions(index_t i) = 0;

        // returns all value positions in R whose first component is i
        virtual value_positions_t get_successor_value_positions(index_t i) = 0;

        // returns the values of all pairs in R whose second component is j
        virtual values_t get_predecessor_values(index_t j) = 0;

        // returns the row numbers of all pairs in R whose second component is j
        virtual indices_t get_predecessor_positions(index_t j) = 0;

        // returns all value positions in R whose second component is j
        virtual value_positions_t get_predecessor_value_positions(index_t j) = 0;

        // returns the values of all pairs (i,j) in R with i1 <= i <= i2 and j1 <= j <= j2
        virtual values_t get_values_in_range(index_t i1, index_t i2, index_t j1, index_t j2) = 0;

        // returns the positions of all pairs (i,j) in R with i1 <= i <= i2 and j1 <= j <= j2
        virtual positions_t get_positions_in_range(index_t i1, index_t i2, index_t j1, index_t j2) = 0;

        // returns the positions and values of all pairs (i,j) in R with i1 <= i <= i2 and j1 <= j <= j2
        virtual value_positions_t get_value_positions_in_range(index_t i1, index_t i2, index_t j1, index_t j2) = 0;

        // returns the values of all pairs in R
        virtual values_t get_all_values() = 0;

        // returns the positions of all pairs in R
        virtual positions_t get_all_positions() = 0;

        // returns all value positions in R
        virtual value_positions_t get_all_value_positions() = 0;

        // checks whether R contains a pair (i,j) with i1 <= i <= i2 and j1 <= j <= j2
        virtual bool contains_value(index_t i1, index_t i2, index_t j1, index_t j2) = 0;

        // returns the number of pairs in R
        virtual size_t count_values() = 0;


        // determines memory consumption (in bytes)
        virtual size_t size_in_bytes() const = 0;

        // prints the parameters (and contents) of the K2Tree
        virtual void print(bool all) = 0;


        /*
         * Method aliases using "relation nomenclature" (similar to the names proposed by Brisaboa et al.)
         */

        // alias of is_not_null()
        virtual bool are_related(index_t i, index_t j) = 0;

        // alias of get_successor_positions()
        virtual indices_t get_successors(index_t i) = 0;

        // alias of get_predecessor_positions()
        virtual indices_t get_predecessors(index_t j) = 0;

        // alias of get_positions_in_range()
        virtual positions_t get_range(index_t i1, index_t i2, index_t j1, index_t j2) = 0;

        // alias of contains_value()
        virtual bool contains_link(index_t i1, index_t i2, index_t j1, index_t j2) = 0;

        // alias of count_values()
        virtual size_t count_links() = 0;

    };

}

#endif //K2TREES_K2TREE_HPP
