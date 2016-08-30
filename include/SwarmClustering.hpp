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

#ifndef SCT_PJ_SWARMCLUSTERING_HPP
#define SCT_PJ_SWARMCLUSTERING_HPP

#include "Base.hpp"
#include "Relation.hpp"


namespace SCT_PJ {
namespace SwarmClustering {

struct SwarmConfig {

    bool outInternals; //-i
    std::string oFileInternals;

    bool outOtus; //-o
    std::string oFileOtus;

    bool outStatistics; //-s
    std::string oFileStatistics;

    bool outSeeds; //-w
    std::string oFileSeeds;

    char sepInternals = '\t';
    char sepOtus = ' ';
    char sepStatistics = '\t';
    char sepAbundance = ';';

    bool noOtuBreaking = false; //-n

    bool dereplicate = false; //-d 0

    unsigned long numExplorers = 1;

};

/**
 * Representation of one OTU.
 *
 * Stored information:
 *  - id: OTU number (positive integer; OTUs are labelled in the order of occurrence)
 *  - numUniqueSequences: number of distinct actual sequences in the cluster
 *  - totalCopyNumber: total copy number of all sequences in the cluster
 *  - seedId: sequence identifier of the OTU seed
 *  - seedCopyNumber: abundance (copy number) of OTU seed
 *  - numSingletons: number of sequences with abundance one in the OTU
 *  - maxGen: maximum number of generations (i.e. number of iterations before the OTU reached its natural limit
 *  - maxRad: maximum radius of the OTU (i.e. accumulated differences between the seed the furthermost sequence in the OTU)
 *  - links: "structure" of the OTU
 */
struct Otu {

    numSeqs_t id;
    numSeqs_t numUniqueSequences;
    numSeqs_t totalCopyNumber;
    numSeqs_t seedId;
    numSeqs_t seedCopyNumber;
    numSeqs_t numSingletons;
    lenSeqs_t maxGen;
    lenSeqs_t maxRad;
    SimpleMatches<lenSeqs_t, numSeqs_t> links;

    numSeqs_t poolId; //only for dereplication

    Otu() {

        id = 0;
        numUniqueSequences = 0;
        totalCopyNumber = 0;
        seedId = 0;
        seedCopyNumber = 0;
        numSingletons = 0;
        maxGen = 0;
        maxRad = 0;
        links = SimpleMatches<lenSeqs_t, numSeqs_t>();
        poolId = 0;

    }

    Otu(numSeqs_t i, numSeqs_t sid, numSeqs_t scn) {

        id = i;
        numUniqueSequences = 0;
        totalCopyNumber = 0;
        seedId = sid;
        seedCopyNumber = scn;
        numSingletons = 0;
        maxGen = 0;
        maxRad = 0;
        links = SimpleMatches<lenSeqs_t, numSeqs_t>();
        poolId = 0;

    }

};

// "point" in the amplicon space during the exploration (characterised by the amplicon id as well as the distance to / number of generations from the root
struct OtuEntry {

    numSeqs_t id;
    lenSeqs_t gen;
    lenSeqs_t rad;

    OtuEntry() {
        id = gen = rad = 0;
    }

    OtuEntry(numSeqs_t i, lenSeqs_t g, lenSeqs_t r) {
        id = i;
        gen = g;
        rad = r;
    }

};

// comparer structures
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

// intended use: sort content of a deque (representing a generation in an OTU) according to the respective abundance of the amplicons in the referenced AmpliconCollection
struct DequeCompareAbund {

    const AmpliconCollection& ac;

    DequeCompareAbund(const AmpliconCollection& coll) : ac(coll) {
        // nothing more to do
    }

    bool operator()(OtuEntry& a, OtuEntry& b) {
        return ac[a.id].abundance > ac[b.id].abundance;
    }

};

// intended use: sort OTUs by their mass ( = total abundance)
struct OtusCompareAbund {
    bool operator()(Otu* a, Otu* b) {
        return a->totalCopyNumber > b->totalCopyNumber;
    }
};

void explorePool(const AmpliconCollection& ac, Matches& matches, std::vector<Otu*>& otus, const SwarmConfig& sc);


//option -i, --internal-structure output of swarm
void outputInternalStructures(const std::string oFile, const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, const char sep = '\t');

//option -o, --output-file output of swarm
void outputOtus(const std::string oFile, const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, const char sep = ' ', const char sepAbundance = ';');

//option -s, --statistics-file output of swarm
void outputStatistics(const std::string oFile, const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, const char sep = '\t');

//option -w, --seeds output of swarm
void outputSeeds(const std::string oFile, const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, const char sepAbundance = '\t');


void outputDereplicate(const AmpliconPools& pools, const std::vector<Otu*>& otus, const SwarmConfig& sc);


void exploreThenOutput(const AmpliconPools& pools, std::vector<Matches*>& allMatches, const SwarmConfig& sc);


void exploreAndOutput(const AmpliconPools& pools, std::vector<Matches*>& allMatches, const SwarmConfig& sc);

}
}

#endif //SCT_PJ_SWARMCLUSTERING_HPP