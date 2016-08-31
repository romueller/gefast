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

#include "../include/SwarmClustering.hpp"

#include <fstream>
#include <queue>
#include <sstream>
#include <thread>
#include <unordered_set>

namespace SCT_PJ {

void SwarmClustering::explorePool(const AmpliconCollection& ac, Matches& matches, std::vector<Otu*>& otus, const SwarmConfig& sc) {

    // determine order of amplicons based on abundance (descending) without invalidating the integer (position) ids of the amplicons
    std::vector<numSeqs_t> index(ac.size());
    std::iota(std::begin(index), std::end(index), 0);
    std::sort(index.begin(), index.end(), IndexCompareAbund(ac));

    Otu* curOtu = 0;
    numSeqs_t numOtu = 0;
    std::vector<bool> visited(ac.size(), false);

    std::deque<OtuEntry> otuRim;
    OtuEntry curSeed, newSeed;
    bool unique;
    std::unordered_set<std::string> nonUniques;
    std::vector<std::pair<numSeqs_t, lenSeqs_t>> next;
    lenSeqs_t lastGen;

    // open new OTU for the amplicon with the highest abundance that is not yet included in an OTU
    for (auto seedIter = index.begin(); seedIter != index.end(); seedIter++) {

        if (!visited[*seedIter]) { // visited amplicons are already included in an OTU

            /* (a) Initialise new OTU with seed */
            curOtu = new Otu(++numOtu, *seedIter, ac[*seedIter].abundance);

            newSeed.id = *seedIter;
            newSeed.gen = 0;
            newSeed.rad = 0;
            otuRim.push_back(newSeed);

            visited[*seedIter] = true;
            nonUniques.clear();

            lastGen = 0;


            /* (b) BFS through 'match space' */
            while (!otuRim.empty()) { // expand current OTU until no further similar amplicons can be added

                if (lastGen != otuRim.front().gen) {
                    std::sort(otuRim.begin(), otuRim.end(), DequeCompareAbund(ac));
                }

                // Get next OTU (sub)seed
                curSeed = otuRim.front();
                otuRim.pop_front();

                unique = true;

                // Update OTU information
                curOtu->totalCopyNumber += ac[curSeed.id].abundance;
                curOtu->numSingletons += (ac[curSeed.id].abundance == 1);

                if (curSeed.gen > curOtu->maxGen) curOtu->maxGen = curSeed.gen;
                if (curSeed.rad > curOtu->maxRad) curOtu->maxRad = curSeed.rad;

                // Consider yet unseen (unvisited) amplicons to continue exploring
                // An amplicon is marked as 'visited' as soon as it occurs the first time as matching partner
                // in order to prevent the algorithm from queueing it more than once from different amplicons.
                next = matches.getMatchesOfAugmented(curSeed.id);
                for (auto matchIter = next.begin(); matchIter != next.end(); matchIter++) {

                    unique &= (matchIter->second != 0);

                    if (!visited[matchIter->first] && (sc.noOtuBreaking || ac[matchIter->first].abundance <= ac[curSeed.id].abundance)) {

                        curOtu->links.add(curSeed.id, matchIter->first, matchIter->second);

                        newSeed.id = matchIter->first;
                        newSeed.gen = curSeed.gen + 1;
                        newSeed.rad = curSeed.rad + matchIter->second;
                        otuRim.push_back(newSeed);
                        visited[matchIter->first] = true;

                    }
                }

                unique = unique || sc.dereplicate || nonUniques.insert(ac[curSeed.id].seq).second;

                curOtu->numUniqueSequences += unique;

                lastGen = curSeed.gen;

            }

            /* (c) Close the no longer extendable OTU */
            otus.push_back(curOtu);

        }

    }

}

void SwarmClustering::outputInternalStructures(const std::string oFile, const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, const char sep) {

    std::ofstream oStream(oFile);
    std::stringstream sStream;

    AmpliconCollection* ac = 0;
    std::queue<std::pair<numSeqs_t, numSeqs_t>> otuRim;
    std::vector<bool> visited;
    std::vector<std::pair<numSeqs_t, lenSeqs_t>> partners;
    std::pair<numSeqs_t, numSeqs_t> curSeed;

    for (numSeqs_t p = 0; p < pools.numPools(); p++) {

        ac = pools.get(p);
        visited = std::vector<bool>(ac->size(), false);

        for (auto otuIter = otus[p].begin(); otuIter != otus[p].end(); otuIter++) {

            curSeed = std::make_pair((*otuIter)->seedId, 0);
            otuRim.push(curSeed);
            visited[curSeed.first] = true;

            while (!otuRim.empty()) {

                curSeed = otuRim.front();
                otuRim.pop();

                partners = (*otuIter)->links.getMatchesOfAugmented(curSeed.first);
                for (auto partnerIter = partners.begin(); partnerIter != partners.end(); partnerIter++) {

                    if (!visited[partnerIter->first]) {

                        sStream << (*ac)[curSeed.first].id << sep << (*ac)[partnerIter->first].id << sep << partnerIter->second << sep << (*otuIter)->id << sep << (curSeed.second + 1) << std::endl;
                        oStream << sStream.rdbuf();

                        otuRim.push(std::make_pair(partnerIter->first, curSeed.second + 1));
                        visited[partnerIter->first] = true;

                    }

                }

            }

        }

    }

    oStream.close();

}

void SwarmClustering::outputOtus(const std::string oFile, const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, const char sep, const char sepAbundance) {

    std::ofstream oStream(oFile);
    std::stringstream sStream;

    AmpliconCollection* ac = 0;
    std::queue<numSeqs_t> otuRim;
    std::vector<bool> visited;
    std::vector<numSeqs_t> partners;
    numSeqs_t curSeed;

    for (numSeqs_t p = 0; p < pools.numPools(); p++) {

        ac = pools.get(p);
        visited = std::vector<bool>(ac->size(), false);

        for (auto otuIter = otus[p].begin(); otuIter != otus[p].end(); otuIter++) {

            curSeed = (*otuIter)->seedId;
            sStream << (*ac)[curSeed].id << sepAbundance << (*ac)[curSeed].abundance;

            otuRim.push(curSeed);
            visited[curSeed] = true;

            while (!otuRim.empty()) {

                curSeed = otuRim.front();
                otuRim.pop();

                partners = (*otuIter)->links.getMatchesOf(curSeed);
                for (auto partnerIter = partners.begin(); partnerIter != partners.end(); partnerIter++) {

                    if (!visited[*partnerIter]) {

                        sStream << sep << (*ac)[*partnerIter].id << sepAbundance << (*ac)[*partnerIter].abundance;

                        otuRim.push(*partnerIter);
                        visited[*partnerIter] = true;

                    }

                }

            }

            sStream << std::endl;
            oStream << sStream.rdbuf();

        }

    }

    oStream.close();

}

void SwarmClustering::outputStatistics(const std::string oFile, const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, const char sep) {

    std::ofstream oStream(oFile);
    std::stringstream sStream;

    AmpliconCollection* ac = 0;

    for (numSeqs_t p = 0; p < pools.numPools(); p++) {

        ac = pools.get(p);

        for (auto otuIter = otus[p].begin(); otuIter != otus[p].end(); otuIter++) {

            sStream << (*otuIter)->numUniqueSequences << sep << (*otuIter)->totalCopyNumber << sep << (*ac)[(*otuIter)->seedId].id << sep << (*otuIter)->seedCopyNumber << sep << (*otuIter)->numSingletons << sep << (*otuIter)->maxGen << sep << (*otuIter)->maxRad << std::endl;
            oStream << sStream.rdbuf();

        }

    }

    oStream.close();

}

void SwarmClustering::outputSeeds(const std::string oFile, const AmpliconPools& pools, const std::vector<std::vector<Otu*>>& otus, const char sepAbundance) {

    std::ofstream oStream(oFile);
    std::stringstream sStream;

    AmpliconCollection* ac = 0;

    for (numSeqs_t p = 0; p < pools.numPools(); p++) {

        ac = pools.get(p);

        for (auto otuIter = otus[p].begin(); otuIter != otus[p].end(); otuIter++) {

            sStream << ">" << (*ac)[(*otuIter)->seedId].id << sepAbundance << (*otuIter)->totalCopyNumber << std::endl << (*ac)[(*otuIter)->seedId].seq << std::endl;
            oStream << sStream.rdbuf();

        }

    }

    oStream.close();

}


void SwarmClustering::outputDereplicate(const AmpliconPools& pools, const std::vector<Otu*>& otus, const SwarmConfig& sc) {

    std::ofstream oStreamInternals;
    std::ofstream oStreamOtus;
    std::ofstream oStreamStatistics;
    std::ofstream oStreamSeeds;
    std::stringstream sStreamInternals;
    std::stringstream sStreamOtus;
    std::stringstream sStreamStatistics;
    std::stringstream sStreamSeeds;

    if (sc.outInternals) oStreamInternals.open(sc.oFileInternals);
    if (sc.outOtus) oStreamOtus.open(sc.oFileOtus);
    if (sc.outStatistics) oStreamStatistics.open(sc.oFileStatistics);
    if (sc.outSeeds) oStreamSeeds.open(sc.oFileSeeds);

    std::vector<numSeqs_t> partners;

    for (auto i = 0; i < otus.size(); i++) {

        Otu& otu = *(otus[i]);
        AmpliconCollection& ac = *(pools.get(otu.poolId));

        if (sc.outInternals) {

            partners = otu.links.getMatchesOf(otu.seedId);
            for (auto partnerIter = partners.begin(); partnerIter != partners.end(); partnerIter++) {

                sStreamInternals << ac[otu.seedId].id << sc.sepInternals << ac[*partnerIter].id << sc.sepInternals << 0 << sc.sepInternals << (i + 1) << sc.sepInternals << 0 << std::endl;
                oStreamInternals << sStreamInternals.rdbuf();
                sStreamInternals.str(std::string());

            }

        }

        if (sc.outOtus) {

            sStreamOtus << ac[otu.seedId].id << sc.sepAbundance << ac[otu.seedId].abundance;

            partners = otu.links.getMatchesOf(otu.seedId);
            for (auto partnerIter = partners.begin(); partnerIter != partners.end(); partnerIter++) {
                sStreamOtus << sc.sepOtus << ac[*partnerIter].id << sc.sepAbundance << ac[*partnerIter].abundance;
            }

            sStreamOtus << std::endl;
            oStreamOtus << sStreamOtus.rdbuf();
            sStreamOtus.str(std::string());

        }

        if (sc.outStatistics) {

            sStreamStatistics << otu.numUniqueSequences << sc.sepStatistics << otu.totalCopyNumber << sc.sepStatistics << ac[otu.seedId].id << sc.sepStatistics << otu.seedCopyNumber << sc.sepStatistics << otu.numSingletons << sc.sepStatistics << 0 << sc.sepStatistics << 0 << std::endl;
            oStreamStatistics << sStreamStatistics.rdbuf();
            sStreamStatistics.str(std::string());

        }

        if (sc.outSeeds) {

            sStreamSeeds << ">" << ac[otu.seedId].id << sc.sepAbundance << otu.totalCopyNumber << std::endl << ac[otu.seedId].seq << std::endl;
            oStreamSeeds << sStreamSeeds.rdbuf();
            sStreamSeeds.str(std::string());

        }

    }

    if (sc.outInternals) oStreamInternals.close();
    if (sc.outOtus) oStreamOtus.close();
    if (sc.outStatistics) oStreamStatistics.close();
    if (sc.outSeeds) oStreamSeeds.close();

}


void SwarmClustering::exploreThenOutput(const AmpliconPools& pools, std::vector<Matches*>& allMatches, const SwarmConfig& sc) {

    std::vector<std::vector<Otu*>> otus(pools.numPools());
    std::thread explorers[sc.numExplorers];
    unsigned long r = 0;
    for (; r + sc.numExplorers <= pools.numPools(); r += sc.numExplorers) {

        for (unsigned long e = 0; e < sc.numExplorers; e++) {
            explorers[e] = std::thread(&SwarmClustering::explorePool, std::ref(*(pools.get(r + e))), std::ref(*(allMatches[r + e])), std::ref(otus[r + e]), std::ref(sc));
        }
        for (unsigned long e = 0; e < sc.numExplorers; e++) {
            explorers[e].join();
        }

    }

    for (unsigned long e = 0; e < pools.numPools() % sc.numExplorers; e++) {
        explorers[e] = std::thread(&SwarmClustering::explorePool, std::ref(*(pools.get(r + e))), std::ref(*(allMatches[r + e])), std::ref(otus[r + e]), std::ref(sc));
    }
    for (unsigned long e = 0; e < pools.numPools() % sc.numExplorers; e++) {
        explorers[e].join();
    }

    // make OTU IDs unique over all pools (so far IDs start at 1 in each pool) and add pool IDs
    numSeqs_t cnt = otus[0].size();
    for (numSeqs_t p = 1; p < pools.numPools(); p++) {

        for (auto otuIter = otus[p].begin(); otuIter != otus[p].end(); otuIter++) {

            (*otuIter)->id += cnt;
            (*otuIter)->poolId = p;

        }

        cnt += otus[p].size();

    }

    if (sc.dereplicate) {

        std::vector<Otu*> flattened(cnt);
        std::vector<Otu*>::iterator iter = flattened.begin();
        for (auto p = 0; p < pools.numPools(); p++) {

            iter = std::move(otus[p].begin(), otus[p].end(), iter);
            otus[p] = std::vector<Otu*>();

        }

        std::sort(flattened.begin(), flattened.end(), OtusCompareAbund());

        outputDereplicate(pools, flattened, sc);

        for (auto otuIter = flattened.begin(); otuIter != flattened.end(); otuIter++) {
                delete *otuIter;
        }

    } else {

        if (sc.outInternals) outputInternalStructures(sc.oFileInternals, pools, otus, sc.sepInternals);
        if (sc.outOtus) outputOtus(sc.oFileOtus, pools, otus, sc.sepOtus, sc.sepAbundance);
        if (sc.outStatistics) outputStatistics(sc.oFileStatistics, pools, otus, sc.sepStatistics);
        if (sc.outSeeds) outputSeeds(sc.oFileSeeds, pools, otus, sc.sepAbundance);

        for (auto pIter = otus.begin(); pIter != otus.end(); pIter++) {
            for (auto otuIter = pIter->begin(); otuIter != pIter->end(); otuIter++) {
                delete *otuIter;
            }
        }

    }

}

void SwarmClustering::exploreAndOutput(const AmpliconPools& pools, std::vector<Matches*>& allMatches, const SwarmConfig& sc) {

    std::ofstream oStreamInternals;
    std::ofstream oStreamOtus;
    std::ofstream oStreamStatistics;
    std::ofstream oStreamSeeds;
    std::stringstream sStreamInternals;
    std::stringstream sStreamOtus;
    std::stringstream sStreamStatistics;
    std::stringstream sStreamSeeds;

    if (!sc.dereplicate && sc.outInternals) oStreamInternals.open(sc.oFileInternals);
    if (!sc.dereplicate && sc.outOtus) oStreamOtus.open(sc.oFileOtus);
    if (!sc.dereplicate && sc.outStatistics) oStreamStatistics.open(sc.oFileStatistics);
    if (!sc.dereplicate && sc.outSeeds) oStreamSeeds.open(sc.oFileSeeds);


    std::vector<numSeqs_t> index;
    Otu* curOtu;
    std::vector<Otu*> otus;
    numSeqs_t numOtu = 0;
    std::vector<bool> visited;
    std::unordered_set<std::string> nonUniques;
    std::deque<OtuEntry> otuRim;
    OtuEntry curSeed, newSeed;
    bool unique;
    std::vector<std::pair<numSeqs_t, lenSeqs_t>> next;
    lenSeqs_t lastGen;


    for (numSeqs_t p = 0; p < pools.numPools(); p++) {

        AmpliconCollection& ac = *(pools.get(p));
        Matches& matches = *(allMatches[p]);

        // determine order of amplicons based on abundance (descending) without invalidating the integer (position) ids of the amplicons
        index = std::vector<numSeqs_t>(ac.size());//std::vector<numSeqs_t> index(ac.size());
        std::iota(std::begin(index), std::end(index), 0);
        std::sort(index.begin(), index.end(), IndexCompareAbund(ac));

        visited = std::vector<bool>(ac.size(), false);

        // open new OTU for the amplicon with the highest abundance that is not yet included in an OTU
        for (auto seedIter = index.begin(); seedIter != index.end(); seedIter++) {

            if (!visited[*seedIter]) { // visited amplicons are already included in an OTU

                /* (a) Initialise new OTU with seed */
                curOtu = new Otu(++numOtu, *seedIter, ac[*seedIter].abundance);
                curOtu->poolId = p;

                if (!sc.dereplicate && sc.outOtus) sStreamOtus << ac[*seedIter].id << sc.sepAbundance << ac[*seedIter].abundance;

                newSeed.id = *seedIter;
                newSeed.gen = 0;
                newSeed.rad = 0;
                otuRim.push_back(newSeed);

                visited[*seedIter] = true;
                nonUniques.clear();

                lastGen = 0;


                /* (b) BFS through 'match space' */
                while (!otuRim.empty()) { // expand current OTU until no further similar amplicons can be added

                    if (lastGen != otuRim.front().gen) {
                        std::sort(otuRim.begin(), otuRim.end(), DequeCompareAbund(ac));
                    }

                    // Get next OTU (sub)seed
                    curSeed = otuRim.front();
                    otuRim.pop_front();

                    unique = true;

                    // Update OTU information
                    curOtu->totalCopyNumber += ac[curSeed.id].abundance;
                    curOtu->numSingletons += (ac[curSeed.id].abundance == 1);

                    if (curSeed.gen > curOtu->maxGen) curOtu->maxGen = curSeed.gen;
                    if (curSeed.rad > curOtu->maxRad) curOtu->maxRad = curSeed.rad;

                    // Consider yet unseen (unvisited) amplicons to continue exploring
                    // An amplicon is marked as 'visited' as soon as it occurs the first time as matching partner
                    // in order to prevent the algorithm from queueing it more than once from different amplicons.
                    next = matches.getMatchesOfAugmented(curSeed.id);
                    for (auto matchIter = next.begin(); matchIter != next.end(); matchIter++) {

                        unique &= (matchIter->second != 0);

                        if (!visited[matchIter->first] && (sc.noOtuBreaking || ac[matchIter->first].abundance <= ac[curSeed.id].abundance)) {

                            if (sc.dereplicate) curOtu->links.add(curSeed.id, matchIter->first, matchIter->second);

                            newSeed.id = matchIter->first;
                            newSeed.gen = curSeed.gen + 1;
                            newSeed.rad = curSeed.rad + matchIter->second;
                            otuRim.push_back(newSeed);
                            visited[matchIter->first] = true;

                            if (!sc.dereplicate && sc.outInternals) {

                                sStreamInternals << ac[curSeed.id].id << sc.sepInternals << ac[matchIter->first].id << sc.sepInternals << matchIter->second << sc.sepInternals << curOtu->id << sc.sepInternals << (curSeed.gen + 1) << std::endl;
                                oStreamInternals << sStreamInternals.rdbuf();
                                sStreamInternals.str(std::string());

                            }
                            if (!sc.dereplicate && sc.outOtus) sStreamOtus << sc.sepOtus << ac[matchIter->first].id << sc.sepAbundance << ac[matchIter->first].abundance;

                        }

                    }

                    unique = unique || sc.dereplicate || nonUniques.insert(ac[curSeed.id].seq).second;

                    curOtu->numUniqueSequences += unique;

                    lastGen = curSeed.gen;

                }

                /* (c) Close the no longer extendable OTU */
                if (!sc.dereplicate && sc.outOtus) {

                    sStreamOtus << std::endl;
                    oStreamOtus << sStreamOtus.rdbuf();
                    sStreamOtus.str(std::string());

                }
                if (!sc.dereplicate && sc.outStatistics) {

                    sStreamStatistics << curOtu->numUniqueSequences << sc.sepStatistics << curOtu->totalCopyNumber << sc.sepStatistics << ac[curOtu->seedId].id << sc.sepStatistics << curOtu->seedCopyNumber << sc.sepStatistics << curOtu->numSingletons << sc.sepStatistics << curOtu->maxGen << sc.sepStatistics << curOtu->maxRad << std::endl;
                    oStreamStatistics << sStreamStatistics.rdbuf();
                    sStreamStatistics.str(std::string());

                }
                if (!sc.dereplicate && sc.outSeeds) {

                    sStreamSeeds << ">" << ac[curOtu->seedId].id << sc.sepAbundance << curOtu->totalCopyNumber << std::endl << ac[curOtu->seedId].seq << std::endl;
                    oStreamSeeds << sStreamSeeds.rdbuf();
                    sStreamSeeds.str(std::string());

                }

                if (sc.dereplicate) {
                    otus.push_back(curOtu);
                } else {
                    delete curOtu;
                }

            }

        }

    }

    if (sc.dereplicate) {

        std::sort(otus.begin(), otus.end(), OtusCompareAbund());

        outputDereplicate(pools, otus, sc);

    }

    if (!sc.dereplicate && sc.outInternals) oStreamInternals.close();
    if (!sc.dereplicate && sc.outOtus) oStreamOtus.close();
    if (!sc.dereplicate && sc.outStatistics) oStreamStatistics.close();
    if (!sc.dereplicate && sc.outSeeds) oStreamSeeds.close();

}



}