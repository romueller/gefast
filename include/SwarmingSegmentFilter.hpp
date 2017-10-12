/*
 * GeFaST
 *
 * Copyright (C) 2016 - 2017 Robert Mueller
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

#ifndef GEFAST_SWARMINGSEGMENTFILTER_HPP
#define GEFAST_SWARMINGSEGMENTFILTER_HPP

#include <thread>

#include "Base.hpp"
#include "Buffer.hpp"
#include "Relation.hpp"
#include "VerificationGotoh.hpp"
#include "SwarmClustering.hpp"

#define CHILDREN_FINDER 1

namespace GeFaST {
namespace SegmentFilter {

#if SUCCINCT

#include "RelationSuccinct.hpp"
typedef K2TreeBinaryRelation<RankedAscendingLabels> SwarmingInvertedIndex;
typedef SharingRollingIndices<RankedAscendingLabels, SwarmingInvertedIndex> SwarmingIndices;

#else

typedef LabelLinkBinaryRelation<StringIteratorPair, numSeqs_t, hashStringIteratorPair, equalStringIteratorPair> SwarmingInvertedIndex;
typedef RollingIndices<SwarmingInvertedIndex> SwarmingIndices;

#endif

/*
 * Looks up the segments of the amplicon in the inverted indices and makes a tally of the found candidates.
 */
void addCandCnts(const Amplicon& amplicon, lenSeqs_t childLen, lenSeqs_t numSegments, std::vector<numSeqs_t>& candCnts, SwarmingIndices& indices,
                 std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive);

/*
 * Applies forward filter and verifies candidates by computing the bounded edit distance.
 */
void verifyCands(const Amplicon& amplicon, const AmpliconCollection& ac, std::vector<numSeqs_t>& candCnts,
                 std::vector<std::pair<numSeqs_t, lenSeqs_t>>& matches, const SwarmClustering::SwarmConfig& sc,
                 lenSeqs_t* M, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP);

/*
 * Applies forward + pipelined backward filtering and verifies candidates by computing the bounded edit distance.
 */
void verifyCandsTwoWay(const Amplicon& amplicon, std::vector<std::string>& segmentStrs, const AmpliconCollection& ac,
                       std::vector<numSeqs_t>& candCnts, std::vector<Substrings>& candSubstrs,
                       std::vector<std::pair<numSeqs_t, lenSeqs_t>>& matches, const SwarmClustering::SwarmConfig& sc,
                       lenSeqs_t* M, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP);



/*
 * Encapsulates the computation of children of given amplicons (see getChildren[TwoWay] methods below) in one thread.
 */
class ChildrenFinder {

public:
    ChildrenFinder(const AmpliconCollection& ac, SwarmingIndices& indices,
                   std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive,
                   const SwarmClustering::SwarmConfig& sc, lenSeqs_t* M, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP);

    std::vector<std::pair<numSeqs_t, lenSeqs_t>> getChildren(const numSeqs_t id);

    void getChildren(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children);

    std::vector<std::pair<numSeqs_t, lenSeqs_t>> getChildrenTwoWay(const numSeqs_t id);

    void getChildrenTwoWay(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children);

private:
    const AmpliconCollection& ac_;
    SwarmingIndices& indices_;
    std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive_;
    const SwarmClustering::SwarmConfig& sc_;

    // arrays for verification computations
    lenSeqs_t* M_;
    val_t* D_;
    val_t* P_;
    lenSeqs_t* cntDiffs_;
    lenSeqs_t* cntDiffsP_;

};

/*
 * Encapsulates the computation of children of given amplicons (see getChildren[TwoWay] methods below) in multiple threads.
 */
class ParallelChildrenFinder {

public:
    ParallelChildrenFinder(const AmpliconCollection& ac, SwarmingIndices& indices,
                           std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive,
                           const lenSeqs_t width, const SwarmClustering::SwarmConfig& sc);

    ~ParallelChildrenFinder();

    std::vector<std::pair<numSeqs_t, lenSeqs_t>> getChildren(const numSeqs_t id);

    void getChildren(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children);

    std::vector<std::pair<numSeqs_t, lenSeqs_t>> getChildrenTwoWay(const numSeqs_t id);

    void getChildrenTwoWay(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children);

private:
    numSeqs_t sendCandsToVerification(const numSeqs_t id, const Amplicon& amplicon, std::vector<numSeqs_t>& candCnts);
    numSeqs_t sendCandsToVerificationTwoWay(const numSeqs_t id, const Amplicon& amplicon, std::vector<std::string>& segmentStrs,
                                            std::vector<numSeqs_t>& candCnts, std::vector<Substrings>& candSubstrs);

    void verify(std::vector<std::pair<numSeqs_t, lenSeqs_t>>& matches, Buffer<Candidate>& buf, lenSeqs_t width);

    void verifyGotoh(std::vector<std::pair<numSeqs_t, lenSeqs_t>>& matches, Buffer<Candidate>& buf, lenSeqs_t width);

    const AmpliconCollection& ac_;
    SwarmingIndices& indices_;
    std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive_;
    const SwarmClustering::SwarmConfig& sc_;

    // interface to the verification threads
    RotatingBuffers<Candidate> cbs_;
    std::vector<std::pair<numSeqs_t, lenSeqs_t>> matches_;
    std::vector<std::thread> verifierThreads_;

    std::mutex mtxMatches_;
    std::mutex mtxDiscard_;
    numSeqs_t numDiscarded_;

};

/*
 * Determine the children (= similar amplicons belonging to the next generation) of the given amplicon using a one-way filter.
 * Employs forward resp. backward filtering depending on the relative lengths of the amplicons.
 */
std::vector<std::pair<numSeqs_t, lenSeqs_t>> getChildren(const numSeqs_t id, const AmpliconCollection& ac, SwarmingIndices& indices,
                                                         std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive,
                                                         lenSeqs_t* M, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP,
                                                         const SwarmClustering::SwarmConfig& sc);
void getChildren(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children, const AmpliconCollection& ac, SwarmingIndices& indices,
                 std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive,
                 lenSeqs_t* M, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP,
                 const SwarmClustering::SwarmConfig& sc);

/*
 * Determine the children (= similar amplicons belonging to the next generation) of the given amplicon using a two-way filter.
 * Employs forward-backward resp. backward-forward filtering depending on the relative lengths of the amplicons.
 */
std::vector<std::pair<numSeqs_t, lenSeqs_t>> getChildrenTwoWay(const numSeqs_t id, const AmpliconCollection& ac, SwarmingIndices& indices,
                                                               std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive,
                                                               lenSeqs_t* M, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP,
                                                               const SwarmClustering::SwarmConfig& sc);
void getChildrenTwoWay(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children, const AmpliconCollection& ac, SwarmingIndices& indices,
                       std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive,
                       lenSeqs_t* M, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP,
                       const SwarmClustering::SwarmConfig& sc);


/*
 * Fill the inverted indices and the substrings archive.
 */
void prepareIndices(const AmpliconCollection& ac, SwarmingIndices& indices,
                    std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<Substrings>>>& substrsArchive,
                    const SwarmClustering::SwarmConfig& sc);

/*
 * Determine OTUs (swarms) like Swarm by using a segment filter.
 *
 * Indexes all amplicons and directly determines the swarms by iteratively adding amplicons
 * (and, thus, emptying the amplicon pool).
 * Delegates the verification of candidates to other threads.
 *
 * Implementation of the clustering strategy proposed in:
 * Mahé et al. (2015), Swarm v2: highly-scalable and high-resolution amplicon clustering
 */
void swarmFilter(const AmpliconCollection& ac, std::vector<SwarmClustering::Otu*>& otus, const SwarmClustering::SwarmConfig& sc);

/*
 * Determine OTUs (swarms) like Swarm by using a segment filter.
 *
 * Indexes all amplicons and directly determines the swarms by iteratively adding amplicons
 * (and, thus, emptying the amplicon pool).
 * Performs the verification of candidates itself.
 *
 * Implementation of the clustering strategy proposed in:
 * Mahé et al. (2015), Swarm v2: highly-scalable and high-resolution amplicon clustering
 */
void swarmFilterDirectly(const AmpliconCollection& ac, std::vector<SwarmClustering::Otu*>& otus, const SwarmClustering::SwarmConfig& sc);

}
}

#endif //GEFAST_SWARMINGSEGMENTFILTER_HPP
