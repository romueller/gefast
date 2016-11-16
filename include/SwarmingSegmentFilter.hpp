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

#ifndef SCT_PJ_SWARMINGSEGMENTFILTER_HPP
#define SCT_PJ_SWARMINGSEGMENTFILTER_HPP

#include <thread>

#include "Base.hpp"
#include "Buffer.hpp"
#include "Relation.hpp"
#include "VerificationGotoh.hpp"
#include "SwarmClustering.hpp"

#define CHILDREN_FINDER 1

namespace SCT_PJ {
namespace SegmentFilter {

    /**
     * Encapsulates the computation of children of given amplicons (see getChildren[TwoWay] methods below)
     * in one thread.
     */
    class ChildrenFinder {

    public:
        ChildrenFinder(const AmpliconCollection& ac, RollingIndices<InvertedIndex>& indices, std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<SegmentFilter::Substrings>>>& substrsArchive, const SwarmClustering::SwarmConfig& sc, lenSeqs_t* M, Verification::val_t* D, Verification::val_t* P, Verification::val_t* Q, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP, lenSeqs_t* cntDiffsQ);

        std::vector<std::pair<numSeqs_t, lenSeqs_t>> getChildren(const numSeqs_t id);

        void getChildren(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children);

        std::vector<std::pair<numSeqs_t, lenSeqs_t>> getChildrenTwoWay(const numSeqs_t id);

        void getChildrenTwoWay(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children);

    private:
        const AmpliconCollection& ac_;
        RollingIndices<InvertedIndex>& indices_;
        std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<SegmentFilter::Substrings>>>& substrsArchive_;
        const SwarmClustering::SwarmConfig& sc_;

        lenSeqs_t* M_;
        Verification::val_t* D_;
        Verification::val_t* P_;
        Verification::val_t* Q_;
        lenSeqs_t* cntDiffs_;
        lenSeqs_t* cntDiffsP_;
        lenSeqs_t* cntDiffsQ_;

    };

    /**
     * Encapsulates the computation of children of given amplicons (see getChildren[TwoWay] methods below)
     * in multiple threads.
     */
    class ParallelChildrenFinder {

    public:
        ParallelChildrenFinder(const AmpliconCollection& ac, RollingIndices<InvertedIndex>& indices, std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<SegmentFilter::Substrings>>>& substrsArchive, const lenSeqs_t width, const SwarmClustering::SwarmConfig& sc);

        ~ParallelChildrenFinder();

        std::vector<std::pair<numSeqs_t, lenSeqs_t>> getChildren(const numSeqs_t id);

        void getChildren(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children);

        std::vector<std::pair<numSeqs_t, lenSeqs_t>> getChildrenTwoWay(const numSeqs_t id);

        void getChildrenTwoWay(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children);

    private:
        void verify(std::vector<std::pair<numSeqs_t, lenSeqs_t>>& matches, Buffer<Candidate>& buf, lenSeqs_t width);

        void verifyGotoh(std::vector<std::pair<numSeqs_t, lenSeqs_t>>& matches, Buffer<Candidate>& buf, lenSeqs_t width);

        const AmpliconCollection& ac_;
        RollingIndices<InvertedIndex>& indices_;
        std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<SegmentFilter::Substrings>>>& substrsArchive_;
        const SwarmClustering::SwarmConfig& sc_;

        RotatingBuffers<Candidate> cbs_;
        std::vector<std::pair<numSeqs_t, lenSeqs_t>> matches_;
        std::vector<std::thread> verifierThreads_;

        std::mutex mtxMatches_;
        std::mutex mtxDiscard_;
        numSeqs_t numDiscarded_;

    };

    /**
     * Determine the children (= similar amplicons belonging to the next generation) of the given amplicon using a one-way filter.
     * Employs forward resp. backward filtering depending on the relative lengths of the amplicons.
     */
    std::vector<std::pair<numSeqs_t, lenSeqs_t>> getChildren(const numSeqs_t id, const AmpliconCollection& ac, RollingIndices<InvertedIndex>& indices, std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<SegmentFilter::Substrings>>>& substrsArchive, lenSeqs_t M[], Verification::val_t D[], Verification::val_t P[], Verification::val_t Q[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[], const SwarmClustering::SwarmConfig& sc);
    void getChildren(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children, const AmpliconCollection& ac, RollingIndices<InvertedIndex>& indices, std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<SegmentFilter::Substrings>>>& substrsArchive, lenSeqs_t M[], Verification::val_t D[], Verification::val_t P[], Verification::val_t Q[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[], const SwarmClustering::SwarmConfig& sc);

    /**
     * Determine the children (= similar amplicons belonging to the next generation) of the given amplicon using a two-way filter.
     * Employs forward-backward resp. backward-forward filtering depending on the relative lengths of the amplicons.
     */
    std::vector<std::pair<numSeqs_t, lenSeqs_t>> getChildrenTwoWay(const numSeqs_t id, const AmpliconCollection& ac, RollingIndices<InvertedIndex>& indices, std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<SegmentFilter::Substrings>>>& substrsArchive, lenSeqs_t M[], Verification::val_t D[], Verification::val_t P[], Verification::val_t Q[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[], const SwarmClustering::SwarmConfig& sc);
    void getChildrenTwoWay(const numSeqs_t id, std::vector<std::pair<numSeqs_t, lenSeqs_t>>& children, const AmpliconCollection& ac, RollingIndices<InvertedIndex>& indices, std::unordered_map<lenSeqs_t, std::unordered_map<lenSeqs_t, std::vector<SegmentFilter::Substrings>>>& substrsArchive, lenSeqs_t M[], Verification::val_t D[], Verification::val_t P[], Verification::val_t Q[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[], const SwarmClustering::SwarmConfig& sc);

    /**
     * Determine OTUs (swarms) like swarm by using a segment filter.
     *
     * Indexes all amplicons and directly determines the swarms by iteratively adding amplicons
     * (and, thus, emptying the amplicon pool).
     * Delegates the verification of candidates to other threads.
     */
    void swarmFilter(const AmpliconCollection& ac, std::vector<SwarmClustering::Otu*>& otus, const SwarmClustering::SwarmConfig& sc);

    /**
     * Determine OTUs (swarms) like swarm by using a segment filter.
     *
     * Indexes all amplicons and directly determines the swarms by iteratively adding amplicons
     * (and, thus, emptying the amplicon pool).
     * Performs the verification of candidates itself.
     */
    void swarmFilterDirectly(const AmpliconCollection& ac, std::vector<SwarmClustering::Otu*>& otus, const SwarmClustering::SwarmConfig& sc);

}
}

#endif //SCT_PJ_SWARMINGSEGMENTFILTER_HPP
