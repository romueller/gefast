/*
 * SCT-PJ
 *
 * Copyright (C) 2016 Robert Mueller
 *
 * TODO add licence text (e.g. GNU (Affero) General Public License)
 *
 * Contact: Robert Mueller <romueller@techfak.uni-bielefeld.de>
 * Faculty of Technology, Bielefeld University,
 * PO box 100131, DE-33501 Bielefeld, Germany
 */

#ifndef SCT_PJ_BASE_HPP
#define SCT_PJ_BASE_HPP

#include <map>
#include <string>
#include <vector>

namespace SCT_PJ {

// type for everything related to counts of amplicons / sequences
typedef unsigned long numSeqs_t;

// type for everything related to the length of the sequence of an amplicon
typedef unsigned long lenSeqs_t;



// =====================================================
//              Data type for single amplicons
// =====================================================

struct Amplicon {

    std::string id;
    std::string seq;
    numSeqs_t abundance;

    Amplicon() {

        id = "";
        seq = "";
        abundance = 0;

    }

    Amplicon(std::string i, std::string s, numSeqs_t a) {

        id = i;
        seq = s;
        abundance = a;

    }

};

// comparer structures for Amplicon structures
struct AmpliconCompareAlph;
struct AmpliconCompareLen;
struct AmpliconCompareAbund;
struct AmpliconSeqEqual;


// =====================================================
//           Data types for multiple amplicons
// =====================================================

// (unspecific) type for referring to multiple amplicons
typedef std::vector<Amplicon> AmpliconCollection;


// comparer structure
// intended use: sort indices [1:n] according to the respective abundances of the amplicons in the referenced AmpliconCollection
struct IndexCompareAbund {

    const AmpliconCollection& ac;

    IndexCompareAbund(const AmpliconCollection& coll) : ac(coll) {
        // nothing more to do
    }

    bool operator()(numSeqs_t a, numSeqs_t b) {
        return ac[a].abundance > ac[b].abundance;
    }

};


// multiple amplicon collections
// intended use: see pool() method of LengthGroups
class AmpliconPools {

public:
    AmpliconPools();

    ~AmpliconPools();

    // append specified AmpliconCollection as a new pool
    void add(AmpliconCollection* ac);

    // return pointer to pool with the specified index (or null pointer if i is too large)
    AmpliconCollection* get(const lenSeqs_t i) const;

    lenSeqs_t numPools() const;

    // return total number of amplicons in all pools
    numSeqs_t numAmplicons() const;

private:
    std::vector<AmpliconCollection*> pools_;

};


// description of a subset of an AmpliconCollection suitable for segment filtering
// Forward filtering (beginIndex <= beginMatch < end, reading left to right):
// Amplicons with an index in [beginMatch, end) make up the actual content of the subpool, i.e. these are the amplicons to be filtered.
// The amplicons with an index in [beginIndex, beginMatch) (if any) are used for indexing only,
// in order to store the same information in the inverted indices as if we would process the whole pool at once.
//
// Backward filtering (beginMatch < beginIndex <= end, reading right to left):
// Amplicons with an index in [beginIndex, end) are used for indexing only, while those with an index in [beginMatch, beginIndex)
// have to be filtered. Here, starting at index end-1 and going down to beginIndex (inclusive) we are only indexing.
// Then, starting at beginIndex-1 and going down to beginMatch (inclusive) we are filtering and indexing.
struct Subpool {

    numSeqs_t beginIndex; // inclusive
    numSeqs_t beginMatch; // inclusive
    numSeqs_t end; // exclusive

    Subpool() {

        beginIndex = 0;
        beginMatch = 0;
        end = 0;

    }

    Subpool(numSeqs_t bi, numSeqs_t bm, numSeqs_t e) {

        beginIndex = bi;
        beginMatch = bm;
        end = e;

    }

};

std::vector<Subpool> getSubpoolBoundaries(const AmpliconCollection& ac, const numSeqs_t num, const lenSeqs_t threshold);


std::vector<Subpool> getSubpoolBoundariesBackward(const AmpliconCollection& ac, const numSeqs_t num, const lenSeqs_t threshold);



// amplicons grouped by (and accessed by) their (sequence) length
// intended use: during preprocessing (for more efficiency)
class LengthGroups {

public:
    LengthGroups();

    LengthGroups(const AmpliconCollection& amplicons);

    // add specified amplicon in the corresponding group based on its length (opens up new group if necessary)
    void add(const Amplicon& ampl);

    // return reference to all length groups
    std::map<lenSeqs_t, AmpliconCollection>& getGroups();

    /*
     * 	Split the total amplicon collection at 'breaks' into 'pools'.
     *	A 'break' (conceptually) lies between two neighboured length groups whose difference
     *	w.r.t. the length of the contained amplicon sequences is strictly larger than the specified threshold.
     *	There is also a break before the first and after the last length group.
     *	Then, a pool is a collection of amplicons created by combining all the length groups between two
     *	neighboured breaks. The amplicons within a pool are sorted by length and (then) alphabet in ascending order.
     *	Empties the LengthGroups structure.
     */
    AmpliconPools* pool(const lenSeqs_t threshold);

private:
    std::map<lenSeqs_t, AmpliconCollection> groups_;

    // sort every group separately (alphabetically, ascending)
    void sortInGroupsByLength();

};


// =====================================================
//                          Misc
// =====================================================

// pair of amplicon 'ids', amplicons are potentially similar (have passed the filter, but not yet verified)
typedef std::pair<numSeqs_t, numSeqs_t> Candidate;

/*
* Collection of (inverted) indices for the segment filter.
*
* The inverted indices are arranged in a grid where
* the columns correspond to segments and
* the rows correspond to sequence lengths.
*
* During the execution of the segment filter,
* 'older' rows can be removed once they correspond to
* sequences too short to be able to provide candidates for
* the current (and future) sequences.
*
*/
template<typename T>
class RollingIndices {
public:
    typedef std::vector<T> Row;


    RollingIndices(lenSeqs_t t, lenSeqs_t w, bool f) {

        threshold_ = t;
        width_ = w;
        forward_ = f;

        empty_ = T();
        emptyRow_ = Row(0);

    }


    // return the indices for the specified length
    Row& getIndicesRow(const lenSeqs_t len) {

        auto iter = indices_.find(len);

        if (iter != indices_.end()) {
            return iter->second;
        } else {
            return emptyRow_;
        }

    }


    // return the index corresponding to the specified length and segment
    T& getIndex(const lenSeqs_t len, const lenSeqs_t i) {

        if (i >= width_) return empty_;

        auto iter = indices_.find(len);

        if (iter != indices_.end()) {
            return (iter->second)[i];
        } else {
            return empty_;
        }

    }


    // add new row (and remove then outdated rows)
    void roll(const lenSeqs_t len) {

        if (indices_.find(len) == indices_.end()) {

            indices_[len] = Row(width_);
            shrink(len);

        }


    }


    // remove outdated rows
    void shrink(const lenSeqs_t cur) {

        if (forward_) {

            auto bound = indices_.lower_bound(cur - threshold_);
            indices_.erase(indices_.begin(), bound);

        } else {

            auto bound = indices_.upper_bound(cur + threshold_);
            indices_.erase(bound, indices_.end());

        }

    }


private:

    lenSeqs_t threshold_; // maximal number of rows
    lenSeqs_t width_; // number of columns / segments per row

    std::map<lenSeqs_t, Row> indices_; // indices grid

    T empty_; // empty (dummy) index returned for out-of-bounds queries
    Row emptyRow_; // empty (dummy) row returned for out-of-bounds queries

    bool forward_; // flag indicating of rolling forwards (increasingly larger lengths are 'inserted') or backwards (shorter lengths are 'inserted')

};


}


#endif //SCT_PJ_BASE_HPP