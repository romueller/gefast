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

namespace SCT_PJ {

void SwarmClustering::explorePool(const AmpliconCollection& ac, Matches& matches, std::vector<Otu*>& otus) {

    // determine order of amplicons based on abundance (descending) without invalidating the integer (position) ids of the amplicons
    std::vector<numSeqs_t> index(ac.size());
    std::iota(std::begin(index), std::end(index), 0);
    std::sort(index.begin(), index.end(), IndexCompareAbund(ac));

    Otu* curOtu = 0;
    numSeqs_t numOtu = 0;
    std::vector<bool> visited(ac.size(), false);

    std::queue<OtuEntry> otuRim;
    OtuEntry curSeed, newSeed;
    bool unique;
    std::vector<std::pair<numSeqs_t, lenSeqs_t>> next;

    // open new OTU for the amplicon with the highest abundance that is not yet included in an OTU
    for (auto seedIter = index.begin(); seedIter != index.end(); seedIter++) {

        if (!visited[*seedIter]) { // visited amplicons are already included in an OTU

            /* (a) Initialise new OTU with seed */
            curOtu = new Otu(++numOtu, *seedIter, ac[*seedIter].abundance);

            newSeed.id = *seedIter;
            newSeed.gen = 0;
            newSeed.rad = 0;
            otuRim.push(newSeed);

            visited[*seedIter] = true;


            /* (b) BFS through 'match space' */
            while (!otuRim.empty()) { // expand current OTU until no further similar amplicons can be added

                // Get next OTU (sub)seed
                curSeed = otuRim.front();
                otuRim.pop();

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

                    if (!visited[matchIter->first]) {

                        curOtu->links.add(curSeed.id, matchIter->first, matchIter->second);

                        newSeed.id = matchIter->first;
                        newSeed.gen = curSeed.gen + 1;
                        newSeed.rad = curSeed.rad + matchIter->second;
                        otuRim.push(newSeed);
                        visited[matchIter->first] = true;

                    }
                }

                curOtu->numUniqueSequences += unique;

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


void SwarmClustering::exploreAndOutput(const AmpliconPools& pools, std::vector<Matches*>& allMatches, /*std::vector<Otu*>& otus,*/ const SwarmConfig& sc) {

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


    std::vector<numSeqs_t> index;
    Otu* curOtu = 0;
    numSeqs_t numOtu = 0;
    std::vector<bool> visited;
    std::queue<OtuEntry> otuRim;
    OtuEntry curSeed, newSeed;
    bool unique;
    std::vector<std::pair<numSeqs_t, lenSeqs_t>> next;


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

                if (sc.outOtus) sStreamOtus << ac[*seedIter].id << sc.sepAbundance << ac[*seedIter].abundance;

                newSeed.id = *seedIter;
                newSeed.gen = 0;
                newSeed.rad = 0;
                otuRim.push(newSeed);

                visited[*seedIter] = true;


                /* (b) BFS through 'match space' */
                while (!otuRim.empty()) { // expand current OTU until no further similar amplicons can be added

                    // Get next OTU (sub)seed
                    curSeed = otuRim.front();
                    otuRim.pop();

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

                        if (!visited[matchIter->first]) {

                            //curOtu->links.add(curSeed.id, matchIter->first, matchIter->second);

                            newSeed.id = matchIter->first;
                            newSeed.gen = curSeed.gen + 1;
                            newSeed.rad = curSeed.rad + matchIter->second;
                            otuRim.push(newSeed);
                            visited[matchIter->first] = true;

                            if (sc.outInternals) {

                                sStreamInternals << ac[curSeed.id].id << sc.sepInternals << ac[matchIter->first].id << sc.sepInternals << matchIter->second << sc.sepInternals << curOtu->id << sc.sepInternals << (curSeed.gen + 1) << std::endl;
                                oStreamInternals << sStreamInternals.rdbuf();

                            }
                            if (sc.outOtus) sStreamOtus << sc.sepOtus << ac[matchIter->first].id << sc.sepAbundance << ac[matchIter->first].abundance;

                        }
                    }

                    curOtu->numUniqueSequences += unique;

                }

                /* (c) Close the no longer extendable OTU */
                //otus.push_back(curOtu);

                if (sc.outOtus) {

                    sStreamOtus << std::endl;
                    oStreamOtus << sStreamOtus.rdbuf();

                }
                if (sc.outStatistics) {

                    sStreamStatistics << curOtu->numUniqueSequences << sc.sepStatistics << curOtu->totalCopyNumber << sc.sepStatistics << ac[curOtu->seedId].id << sc.sepStatistics << curOtu->seedCopyNumber << sc.sepStatistics << curOtu->numSingletons << sc.sepStatistics << curOtu->maxGen << sc.sepStatistics << curOtu->maxRad << std::endl;
                    oStreamStatistics << sStreamStatistics.rdbuf();

                }
                if (sc.outSeeds) {

                    sStreamSeeds << ">" << ac[curOtu->seedId].id << sc.sepAbundance << curOtu->totalCopyNumber << std::endl << ac[curOtu->seedId].seq << std::endl;
                    oStreamSeeds << sStreamSeeds.rdbuf();

                }

            }

        }

    }

    if (sc.outInternals) oStreamInternals.close();
    if (sc.outOtus) oStreamOtus.close();
    if (sc.outStatistics) oStreamStatistics.close();
    if (sc.outSeeds) oStreamSeeds.close();

}



}