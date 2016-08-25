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

#include <ctime>
#include <iostream>
#include <vector>
#include <thread>

#include "include/Base.hpp"
#include "include/Preprocessor.hpp"
#include "include/Relation.hpp"
#include "include/SwarmClustering.hpp"
#include "include/Utility.hpp"
#include "include/Worker.hpp"


namespace SCT_PJ {

// Note: no internal dereplication (unlike swarm), TODO? incorporate
int run(int argc, const char* argv[]) {

    /* Bootstrapping */
    Config<std::string> c = getConfiguration(argc, argv);


    // if no list file is specified with -f / --files, then the first arguments not starting with a dash are assumed to be the input files
    std::vector<std::string> files;
    if (!c.peek(FILE_LIST)) {

        for (int i = 1; i < argc; i++) {

            if (argv[i][0] == '-') break;

            files.push_back(argv[i]);

        }

    } else {
        files = readFileList(c.get(FILE_LIST));
    }


    if (files.size() == 0) {

        std::cerr << "ERROR: No input files specified." << std::endl;
        return 1;

    }
    if (!c.peek(MATCHES_OUTPUT_FILE) && !c.peek(SWARM_OUTPUT_INTERNAL) && !c.peek(SWARM_OUTPUT_OTUS) && !c.peek(SWARM_OUTPUT_STATISTICS) && !c.peek(SWARM_OUTPUT_SEEDS)) {

        std::cerr << "ERROR: No output file specified." << std::endl;
        return 1;

    }


    //TODO remove or redirect to logger
    std::cout << "===== Configuration =====" << std::endl;
    c.print(std::cout);
    std::cout << "=========================" << std::endl;
    std::string jobName = c.get(NAME);
    std::replace(jobName.begin(), jobName.end(), ':', '-');
    if (c.peek(INFO_FOLDER)) writeJobParameters(c.get(INFO_FOLDER) + jobName + ".txt", c, files);


    /* Preprocessing */
    std::cout << "Preprocessing..." << std::endl;
    auto prepStart = clock();
    auto pools = Preprocessor::run(c, files);
    auto prepStop = clock();
    std::cout << "Preprocessing: " << (prepStop - prepStart) / CLOCKS_PER_SEC << std::endl;

    /* Matching */
    std::cout << "Matching..." << std::endl;
    auto matchStart = clock();
    unsigned long numWorkers = std::stoul(c.get(NUM_WORKERS));
    unsigned long numThreadsPerWorkers = std::stoul(c.get(NUM_THREADS_PER_WORKER));
    int mode = std::stoi(c.get(SEGMENT_FILTER));

    lenSeqs_t threshold = std::stoul(c.get(THRESHOLD));
    lenSeqs_t numExtraSegments = std::stoul(c.get(NUM_EXTRA_SEGMENTS));

    std::vector<Matches*> allMatches(pools->numPools()); // one Matches instance per pool
    std::vector<Subpool> subpools; // subpools of one pool
    std::thread workers[numWorkers]; // at most numWorkers per pool

    for (numSeqs_t p = 0; p < pools->numPools(); p++) { std::cout << "Pool #" << p << std::endl; //TODO remove print

        allMatches[p] = new Matches();
        if (mode % 2 == 0) {
            subpools = getSubpoolBoundaries(*(pools->get(p)), numWorkers, threshold);
        } else {
            subpools = getSubpoolBoundariesBackward(*(pools->get(p)), numWorkers, threshold);
        }

        for (unsigned long w = 0; w < subpools.size(); w++) {
            workers[w] = std::thread(&Worker::run, Worker(*(pools->get(p)), subpools[w], *(allMatches[p]), numThreadsPerWorkers, mode), threshold, numExtraSegments);
        }

        for (unsigned long w = 0; w < subpools.size(); w++) {
            workers[w].join();
        }

    }
    auto matchStop = clock();
    std::cout << "Matching: " << (matchStop - matchStart) / CLOCKS_PER_SEC << std::endl;

    /* Postprocessing */
    std::cout << "Postprocessing..." << std::endl;
    auto postpStart = clock();
    if (c.peek(MATCHES_OUTPUT_FILE)) {

        std::cout << "Writing matches..." << std::endl;
        SegmentFilter::writeMatches(c.get(MATCHES_OUTPUT_FILE), *pools, allMatches);

    }

//    // parallel computation of swarm clusters & subsequent generation of outputs
//    std::cout << "Finding OTUs..." << std::endl;
//    std::vector<std::vector<SwarmClustering::Otu*>> otus(pools->numPools());
//    unsigned long numExplorers = numWorkers * numThreadsPerWorkers;
//    std::thread explorers[numExplorers];
//    unsigned long r = 0;
//    for (; r + numExplorers <= pools->numPools(); r += numExplorers) {
//
//        for (unsigned long e = 0; e < numExplorers; e++) {
//            explorers[e] = std::thread(&SwarmClustering::explorePool, std::ref(*(pools->get(r + e))), std::ref(*(allMatches[r + e])), std::ref(otus[r + e]));
//        }
//        for (unsigned long e = 0; e < numExplorers; e++) {
//            explorers[e].join();
//        }
//
//    }
//
//    for (unsigned long e = 0; e < pools->numPools() % numExplorers; e++) {
//        explorers[e] = std::thread(&SwarmClustering::explorePool, std::ref(*(pools->get(r + e))), std::ref(*(allMatches[r + e])), std::ref(otus[r + e]));
//    }
//    for (unsigned long e = 0; e < pools->numPools() % numExplorers; e++) {
//        explorers[e].join();
//    }
//
//    // make OTU IDs unique over all pools (so far IDs start at 1 in each pool)
//    numSeqs_t cnt = otus[0].size();
//    for (numSeqs_t p = 1; p < pools->numPools(); p++) {
//
//        for (auto otuIter = otus[p].begin(); otuIter != otus[p].end(); otuIter++) {
//            (*otuIter)->id += cnt;
//        }
//
//        cnt += otus[p].size();
//
//    }
//
//    // swarm outputs
//    if (c.peek(SWARM_OUTPUT_INTERNAL)) SwarmClustering::outputInternalStructures(c.get(SWARM_OUTPUT_INTERNAL), *pools, otus);
//    if (c.peek(SWARM_OUTPUT_OTUS)) SwarmClustering::outputOtus(c.peek(SWARM_OUTPUT_OTUS), *pools, otus);
//    if (c.peek(SWARM_OUTPUT_STATISTICS)) SwarmClustering::outputStatistics(c.get(SWARM_OUTPUT_STATISTICS), *pools, otus);
//    if (c.peek(SWARM_OUTPUT_SEEDS)) SwarmClustering::outputSeeds(c.get(SWARM_OUTPUT_SEEDS), *pools, otus);

    // interleaved computation of swarm clusters & generation of outputs
    SwarmClustering::SwarmConfig sc;
    sc.outInternals = c.peek(SWARM_OUTPUT_INTERNAL);
    sc.outOtus = c.peek(SWARM_OUTPUT_OTUS);
    sc.outStatistics = c.peek(SWARM_OUTPUT_STATISTICS);
    sc.outSeeds = c.peek(SWARM_OUTPUT_SEEDS);
    if (sc.outInternals) sc.oFileInternals = c.get(SWARM_OUTPUT_INTERNAL);
    if (sc.outOtus) sc.oFileOtus = c.get(SWARM_OUTPUT_OTUS);
    if (sc.outStatistics) sc.oFileStatistics = c.get(SWARM_OUTPUT_STATISTICS);
    if (sc.outSeeds) sc.oFileSeeds = c.get(SWARM_OUTPUT_SEEDS);

    SwarmClustering::exploreAndOutput(*pools, allMatches, sc);
    auto postpStop = clock();
    std::cout << "Postprocessing: " << (postpStop - postpStart) / CLOCKS_PER_SEC << std::endl;


    /* Cleaning up */
    std::cout << "Cleaning up..." << std::endl;
    for (auto iter = allMatches.begin(); iter != allMatches.end(); iter++) {
        delete *iter;
    }
    delete pools;

    return 0;

}

}


int main(int argc, const char* argv[]) {

    SCT_PJ::run(argc, argv);

}
