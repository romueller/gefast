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

#ifndef GEFAST_SWARMCLUSTERING_HPP
#define GEFAST_SWARMCLUSTERING_HPP

#include "Base.hpp"
#include "Relation.hpp"
#include "SegmentFilter.hpp"
#include "Verification.hpp"
#include "VerificationGotoh.hpp"

#define PRINT_INTERNAL_MODIFIED 0

#define FASTIDIOUS_PARALLEL_POOL 1
#define FASTIDIOUS_PARALLEL_CHECK 1

namespace GeFaST {
namespace SwarmClustering {

/**
 * Configuration parameters of the swarm procedure.
 */
struct SwarmConfig {

    // output flags, files & separators
    bool outInternals; // -i
    std::string oFileInternals;

    bool outOtus; // -o
    bool outMothur; // -r
    std::string oFileOtus;

    bool outStatistics; // -s
    std::string oFileStatistics;

    bool outSeeds; // -w
    std::string oFileSeeds;

    bool outUclust; // -u
    std::string oFileUclust;

    char sepInternals = '\t';
    char sepOtus = ' ';
    char sepStatistics = '\t';
    char sepMothur = ','; // for mothur-compatible output
    char sepUclust = '\t';
    std::string sepMothurOtu = "\t"; // for mothur-compatible output
    std::string sepAbundance;

    // clustering options
    bool noOtuBreaking = false; // -n
    bool dereplicate = false; // -d 0
    bool fastidious = false; // -f
    numSeqs_t boundary; // -b
    lenSeqs_t fastidiousThreshold;

    // scoring
    bool useScore = false;
    Verification::Scoring scoring; // -m, -p, -g, -e

    // misc
    lenSeqs_t threshold;
    lenSeqs_t extraSegs;

    bool filterTwoWay = false;

    unsigned long numExplorers = 1;
    unsigned long numThreadsPerExplorer = 1;
    unsigned long numGrafters = 1;

    unsigned long fastidiousCheckingMode = 0;
    unsigned long numThreadsPerCheck = 1;

};


/**
 * Representation of a member (amplicon) of an OTU, also describing its "position" within the OTU.
 * Serves also as a "point" in the amplicon space during the exploration.
 *
 * Stored information:
 *  - id: integer id of the amplicon within its pool (AmpliconCollection)
 *  - parentId: integer id of the parent amplicon within its pool (AmpliconCollection) OR the pool id (only if gen = 0)
 *  - parentDist: (edit) distance between the amplicon and its parent
 *  - gen: generation number of the amplicon
 *  - rad: cumulated differences between the OTU seed and the amplicon
 */
struct OtuEntry {

    numSeqs_t id;
    numSeqs_t parentId;
    lenSeqs_t parentDist;
    lenSeqs_t gen;
    lenSeqs_t rad;

    OtuEntry() {
        id = parentId = parentDist = gen = rad = 0;
    }

    OtuEntry(numSeqs_t i, numSeqs_t pi, lenSeqs_t pd, lenSeqs_t g, lenSeqs_t r) {

        id = i;
        parentId = pi;
        parentDist = pd;
        gen = g;
        rad = r;

    }

};

/**
 * Representation of one OTU.
 *
 * Stored information:
 *  [- id: OTU number (positive integer; OTUs are labelled in the order of occurrence)]
 *  - numUniqueSequences: number of distinct actual sequences in the cluster
 *  - mass: total abundance of all sequences in the cluster
 *  - seedId: sequence identifier of the OTU seed
 *  - seedAbundance: abundance (copy number) of OTU seed
 *  - numSingletons: number of sequences with abundance one in the OTU
 *  - maxGen: maximum number of generations (i.e. number of iterations before the OTU reached its natural limit
 *  - maxRad: maximum radius of the OTU (i.e. accumulated differences between the seed the furthermost sequence in the OTU)
 *  - members: "structure" of the OTU
 */
struct Otu {

//    numSeqs_t id;
    numSeqs_t numUniqueSequences;
    numSeqs_t mass;
    numSeqs_t seedId;
    numSeqs_t seedAbundance;
    numSeqs_t numSingletons;
    lenSeqs_t maxGen;
    lenSeqs_t maxRad;
    std::vector<OtuEntry> members; // always entry for seed at index 0 (with itself as its parent)
    bool attached; // only for fastidious

    numSeqs_t poolId;

    Otu(/*numSeqs_t i, */numSeqs_t sid, numSeqs_t scn) {

//        id = i;
        numUniqueSequences = 0;
        mass = 0;
        seedId = sid;
        seedAbundance = scn;
        numSingletons = 0;
        maxGen = 0;
        maxRad = 0;
        members = std::vector<OtuEntry>(0);
        poolId = 0;
        attached = false;

    }

};

/**
 * Inverted indices for applying the segment filter in the fastidious clustering phase
 * Maps sequence substrings onto OTU members.
 */
typedef SimpleBinaryRelation<std::string, OtuEntry*> InvertedIndexFastidious;
typedef SimpleBinaryRelation<std::string, std::pair<Otu*, OtuEntry*>> InvertedIndexFastidious2;

/**
 * Representation of a grafting candidate.
 * Parent and (potential) child amplicon are represented by the OTU member of their respective OTU.
 * If the grafting takes place, the OTU of the child is supposed to be grafted upon the parent's OTU.
 */
struct GraftCandidate {

    Otu* parentOtu;
    OtuEntry* parentMember;

    Otu* childOtu;
    OtuEntry* childMember;

    GraftCandidate() {

        parentOtu = childOtu = 0;
        parentMember = childMember = 0;

    }

    GraftCandidate(Otu* po, OtuEntry* pm, Otu* co, OtuEntry* cm) {

        parentOtu = po;
        parentMember = pm;
        childOtu = co;
        childMember = cm;

    }

};

/**
 * Representation of several pairs of potentially similar amplicons.
 * The first component of each pair of amplicon 'ids' (like the candidate pairs during matching) is 'parent',
 * while each entry in the 'children' vector is the second component of one pair.
 *
 * Additional information on the parent are stored in order to change the grafting candidate information
 * of the children when necessary.
 */
struct CandidateFastidious {

    numSeqs_t parent;
    Otu* parentOtu;
    OtuEntry* parentMember;

    std::vector<numSeqs_t> children;

    CandidateFastidious() {

        parent = 0;
        parentOtu =  0;
        parentMember = 0;

    }

    CandidateFastidious(numSeqs_t p, Otu* po, OtuEntry* pm) {

        parent = p;
        parentOtu = po;
        parentMember = pm;

    }

};

/*
 * Diverse comparer structures for the two swarm clustering phases.
 *
 * They can make use of the input order (rank, ascending) of amplicons in order to break ties and to mimic swarm's behaviour.
 * The usage can be toggled on / off through setting INPUT_RANK (in Base.hpp) to 1 resp. 0.
 */
// Sort indices [1:n] according to the respective abundances (descending) of the amplicons in the referenced AmpliconCollection
struct CompareIndicesAbund {

    const AmpliconCollection& ac;

    CompareIndicesAbund(const AmpliconCollection& coll) : ac(coll) {
        // nothing more to do
    }

    bool operator()(numSeqs_t a, numSeqs_t b) {
        return (ac[a].abundance > ac[b].abundance)
#if INPUT_RANK
            || ((ac[a].abundance == ac[b].abundance) && (ac[a].rank < ac[b].rank))
#endif
                ;
    }

};

// Sort OTU members according to the respective abundance (descending) of the amplicons in the referenced AmpliconCollection
struct CompareOtuEntriesAbund {

    const AmpliconCollection& ac;

    CompareOtuEntriesAbund(const AmpliconCollection& coll) : ac(coll) {
        // nothing more to do
    }

    bool operator()(const OtuEntry& a, const OtuEntry& b) {
        return (ac[a.id].abundance > ac[b.id].abundance)
#if INPUT_RANK
            || ((ac[a.id].abundance == ac[b.id].abundance) && (ac[a.id].rank < ac[b.id].rank))
#endif
                ;
    }

};

// Sort OTUs by their mass (descending)
struct CompareOtusMass {

    const AmpliconPools& pools;

    CompareOtusMass(const AmpliconPools& ap) : pools(ap) {
        // nothing more to do
    }

    bool operator()(const Otu* a, const Otu* b) {
        return (a->mass > b->mass)
#if INPUT_RANK
            || ((a->mass == b->mass) && ((*pools.get(a->poolId))[a->seedId].rank < (*pools.get(b->poolId))[b->seedId].rank))
#endif
                ;
    }
};

// Sort OTUs by the abundance of their seeds (descending)
struct CompareOtusSeedAbund {

    const AmpliconPools& pools;

    CompareOtusSeedAbund(const AmpliconPools& ap) : pools(ap) {
        // nothing more to do
    }

    bool operator()(const Otu* a, const Otu* b) {
        return (a->seedAbundance > b->seedAbundance)
#if INPUT_RANK
            || ((a->seedAbundance == b->seedAbundance) && ((*pools.get(a->poolId))[a->seedId].rank < (*pools.get(b->poolId))[b->seedId].rank))
#endif
                ;
    }
};

// Sort grafting candidates by the abundances of the parents and, if necessary, break ties through the abundances of the children (both descending)
struct CompareGraftCandidatesAbund {

    const AmpliconPools& pools;

    CompareGraftCandidatesAbund(const AmpliconPools& p) : pools(p) {
        // nothing more to do
    }

#if INPUT_RANK
    inline bool compareMember(const AmpliconCollection& poolA, numSeqs_t memberIdA, const AmpliconCollection& poolB, numSeqs_t memberIdB) {
        return (poolA[memberIdA].abundance > poolB[memberIdB].abundance) || ((poolA[memberIdA].abundance == poolB[memberIdB].abundance) && (poolA[memberIdA].rank < poolB[memberIdB].rank));
    }
#endif

    bool operator()(const GraftCandidate& a, const GraftCandidate& b) {

        AmpliconCollection& parentPoolA = *pools.get(a.parentOtu->poolId);
        AmpliconCollection& parentPoolB = *pools.get(b.parentOtu->poolId);
        AmpliconCollection& childPoolA = *pools.get(a.childOtu->poolId);
        AmpliconCollection& childPoolB = *pools.get(b.childOtu->poolId);

        return
#if INPUT_RANK
         compareMember(parentPoolA, a.parentMember->id, parentPoolB, b.parentMember->id)
                || (parentPoolA[a.parentMember->id].rank == parentPoolB[b.parentMember->id].rank && compareMember(childPoolA, a.childMember->id, childPoolB, b.childMember->id));
#else
        (parentPoolA[a.parentMember->id].abundance > parentPoolB[b.parentMember->id].abundance)
                || ((parentPoolA[a.parentMember->id].abundance == parentPoolB[b.parentMember->id].abundance) && (childPoolA[a.childMember->id].abundance > childPoolB[b.childMember->id].abundance));
#endif

    }

};



/**
 * Determine all OTUs for the given amplicons by exploring the possible links (matches).
 * The found OTUs are returned via the referenced OTU vector.
 */
void explorePool(const AmpliconCollection& ac, Matches& matches, std::vector<Otu*>& otus, const SwarmConfig& sc);


/* Implementation 1 of fastidious clustering */

/**
 * Index the amplicons of the given (light) OTU and prepares the child information of grafting candidate entries.
 * Potentially reuses already computed information on segment positions through segmentsArchive.
 */
void fastidiousIndexOtu(RollingIndices<InvertedIndexFastidious>& indices, std::unordered_map<lenSeqs_t, SegmentFilter::Segments>& segmentsArchive, const AmpliconCollection& ac, Otu& otu, std::vector<GraftCandidate>& graftCands, const SwarmConfig& sc);

/**
 * Verify the potentially similar amplicons arriving at a candidate buffer and,
 * if appropriate, change the grafting candidate information of the child amplicons.
 * Amplicons are similar if their edit distance is below the given threshold.
 */
void verifyFastidious(const AmpliconPools& pools, const AmpliconCollection& acOtus, const AmpliconCollection& acIndices, std::vector<GraftCandidate>& graftCands, Buffer<CandidateFastidious>& buf, const lenSeqs_t width, const lenSeqs_t t, std::mutex& mtx);

/**
 * Verify the potentially similar amplicons arriving at a candidate buffer and,
 * if appropriate, change the grafting candidate information of the child amplicons.
 * Amplicons are similar if the number of differences (mismatches, insertions, deletions)
 * in the optimal alignment for the given scoring function is below the given threshold.
 */
void verifyGotohFastidious(const AmpliconPools& pools, const AmpliconCollection& acOtus, const AmpliconCollection& acIndices, std::vector<GraftCandidate>& graftCands, Buffer<CandidateFastidious>& buf, const lenSeqs_t width, const lenSeqs_t t, const Verification::Scoring& scoring, std::mutex& mtx);

/**
 * Apply a (forward) segment filter on the amplicons from the heavy OTUs of the current pool using the indexed amplicons of light OTUs.
 * Determines the parent information of the grafting candidates.
 *
 * The method with the suffix 'Directly' verifies the candidates itself directly when they occur and does not hand them over to verifier threads through a buffer.
 */
void fastidiousCheckOtus(RotatingBuffers<CandidateFastidious>& cbs, const std::vector<Otu*>& otus, const AmpliconCollection& acOtus, RollingIndices<InvertedIndexFastidious>& indices, const AmpliconCollection& acIndices, std::vector<GraftCandidate>& graftCands, const SwarmConfig& sc);
void fastidiousCheckOtusDirectly(const AmpliconPools& pools, const std::vector<Otu*>& otus, const AmpliconCollection& acOtus, RollingIndices<InvertedIndexFastidious>& indices, const AmpliconCollection& acIndices, std::vector<GraftCandidate>& graftCands, const lenSeqs_t width, std::mutex& graftCandsMtx, const SwarmConfig& sc);

/**
 * Check for grafting candidates using a segment filter and multiple verifier threads.
 * Looks for grafting candidates for amplicons from 'acIndices' among the amplicons from 'acOtus'.
 */
void checkAndVerify(const AmpliconPools& pools, const std::vector<Otu*>& otus, const AmpliconCollection& acOtus, RollingIndices<InvertedIndexFastidious>& indices, const AmpliconCollection& acIndices, std::vector<GraftCandidate>& graftCands, const lenSeqs_t width, std::mutex& mtx, const SwarmConfig& sc);

/**
 * Determine the grafting candidates of the amplicons from all pools.
 */
void determineGrafts(const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, std::vector<GraftCandidate>& allGraftCands, const numSeqs_t p, std::mutex& allGraftCandsMtx, const SwarmConfig& sc);

/**
 *  Graft light OTUs onto heavy OTUs by "simulating virtual amplicons".
 *  To this end, index the amplicons from the light OTUs for the segment filter with a doubled threshold.
 *  Then, search for matches of the amplicons from the heavy OTUs among the indexed amplicons.
 *  Based on the resulting grafting candidates, light OTUs are grafted upon heavy ones:
 *   - The (final) grafting partner of an amplicon from a light OTU, is a matching amplicon with the highest abundance.
 *   - Each light OTU can be grafted upon at most one heavy OTU (even though there can be grafting candidates for several amplicons of the light OTU).
 *   - Grafting candidates with a higher parent amplicon abundance (and, for ties, higher child amplicon abundance) have a higher priority.
 */
void graftOtus(numSeqs_t& maxSize, numSeqs_t& numOtus, const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, const SwarmConfig& sc); //TODO? other segment filters (forward-backward etc.)

#if 0
/* Implementation 2 of fastidious clustering */

/**
 * Store child information of grafting candidate entries for amplicons of all OTUs of the current pool.
 */
void prepareGraftInfos(const numSeqs_t poolSize, const std::vector<Otu*>& otus, std::vector<GraftCandidate>& curGraftCands, std::vector<GraftCandidate>& nextGraftCands, const SwarmConfig& sc);

/**
 * Shifts the indexing window so that afterwards all amplicons with a length in [len - sc.fastidiousThreshold : len + sc.fastidiousThreshold] are indexed.
 * This may advances the indexing iterator beyond the pool containing the amplicons we are currently filtering.
 * In this case, the vector containing child information of grafting candidates for the next pool is (and has to be) initialised,
 * because these information are also used to decide whether an amplicon comes from a light or a heavy OTU.
 * However, the indexing iterator can only move ahead one pool, because the pools are separated by "gaps" (in terms of sequence length) of at least sc.threshold
 * and doubling the threshold for the fastidious stage thus cannot lead to matches beyond the directly neighbouring pools.
 */
AmpliconCollection::iterator shiftIndexWindow(RollingIndices<InvertedIndexFastidious2>& indices, const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, const numSeqs_t poolIndex, const AmpliconCollection::iterator indexIter, std::vector<GraftCandidate>& curGraftCands, std::vector<GraftCandidate>& nextGraftCands, const lenSeqs_t len, const bool forerunning, const SwarmConfig& sc);

/**
 *  Graft light OTUs onto heavy OTUs by "simulating virtual amplicons".
 *  To this end, we iterate over the amplicons of all pools (in the order of increasing length).
 *  For fastidious threshold T and an amplicon of length L, we have to consider all other amplicons with a length in [L - T : L + T].
 *  While iterating, we therefore keep the relevant amplicons in the inverted indices by shifting an indexing window appropriately whenever moving on to longer amplicons.
 *  Here, we are indexing amplicons from heavy OTUs and filter amplicons from light OTUs.
 *  Based on the resulting grafting candidates, light OTUs are grafted upon heavy ones:
 *   - The (final) grafting partner of an amplicon from a light OTU, is a matching amplicon with the highest abundance.
 *   - Each light OTU can be grafted upon at most one heavy OTU (even though there can be grafting candidates for several amplicons of the light OTU).
 *   - Grafting candidates with a higher parent amplicon abundance (and, for ties, higher child amplicon abundance) have a higher priority.
 */
void graftOtus2(numSeqs_t& maxSize, numSeqs_t& numOtus, const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, const SwarmConfig& sc);
#endif

/**
 * Determine overall statistics, start fastidious clustering phase (if requested) and output the results.
 */
void processOtus(const AmpliconPools& pools, std::vector<std::vector<Otu*>>& otus, const SwarmConfig& sc);

/**
 * Cluster amplicons according to swarm's iterative strategy (based on the provided matching information) and generates the requested outputs.
 * Supports also swarm's dereplication and fastidious clustering options.
 *
 * Version one uses the iterative segment filter (as proposed by Li et al.) and determines the OTUs from all matches.
 * Version two uses a "full index" version of the segment filter and directly determines the OTUs (like swarm).
 */
void cluster(const AmpliconPools& pools, std::vector<Matches*>& allMatches, const SwarmConfig& sc);
void cluster(const AmpliconPools& pools, const SwarmConfig& sc);


/**
 * Write the links of the given OTUs to file (corresponds to output of swarm's option -i).
 * Each line contains one link represented through the
 * (1) amplicon id of the parent,
 * (2) amplicon id of the child,
 * (3) distance,
 * (4) OTU id, and
 * (5) generation number of the child amplicon,
 * separated by the given separator.
 */
void outputInternalStructures(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const char sep);

/**
 * Writes the members of the given OTUs to file (corresponds to output of swarm's option -o).
 * Each line contains the members of one OTU represented through amplicon id and abundance.
 * The members are separated via sep, while id and abundance are separated via sepAbundance.
 */
void outputOtus(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const char sep, const std::string sepAbundance);

/**
 * Writes the members of the given OTUs to file (corresponds to output of swarm's option -o with -r).
 * On a single line, the members of all OTUs are represented through amplicon id and abundance.
 * The members of one OTU and the OTUs themselves are separated via sep resp. sepOtu, while id and abundance are separated via sepAbundance.
 */
void outputOtusMothur(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const lenSeqs_t threshold, const numSeqs_t numOtusAdjusted, const char sep, const std::string sepOtu, const std::string sepAbundance);

/**
 * Write the statistics of the given OTUs to file (corresponds to output of swarm's option -s).
 * Each line contains the statistics of one OTU represented through the
 * (1) number of unique sequences,
 * (2) mass of the OTU,
 * (3) amplicon id of the seed,
 * (4) abundance of the seed amplicon,
 * (5) number of singletons,
 * (6) number of iterations before the OTU reached its natural limits, and
 * (7) number of cumulated differences between the seed and the furthermost amplicon,
 * separated by the given separator.
 */
void outputStatistics(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const char sep);

/**
 * Write the seeds of the given OTUs to file (corresponds to output of swarm's option -w).
 * Each seed comprises two lines.
 * The first line contains the amplicon id of the seed preceded by '>' and followed by a separator and the mass of the OTU.
 * The second line describes the sequence of the seed amplicon.
 */
void outputSeeds(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const std::string sepAbundance);

/**
 * Writes the clustering results in a uclust-like format to file (corresponds to output of swarm's option -u).
 */
void outputUclust(const std::string oFile, const AmpliconPools& pools, const std::vector<Otu*>& otus, const SwarmConfig& sc);

/**
 * Write the requested dereplication outputs to file (SwarmConfig stores information on which are requested).
 * The outputs correspond to swarm's outputs (-i, -o, -s, -w, -u) when running it with -d 0.
 *
 * Minor differences to non-dereplicating output options:
 * -i: generation of child amplicons is also (always) 0
 * -o: <no differences>
 * -s: first column describes the number of amplicons in the OTU (not the number of unique sequences)
 * -w: <no differences>
 * -u: <no differences>
 */
void outputDereplicate(const AmpliconPools& pools, const std::vector<Otu*>& otus, const SwarmConfig& sc);
}
}

#endif //GEFAST_SWARMCLUSTERING_HPP