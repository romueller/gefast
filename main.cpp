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

    SwarmClustering::SwarmConfig sc;
    sc.outInternals = c.peek(SWARM_OUTPUT_INTERNAL);
    sc.outOtus = c.peek(SWARM_OUTPUT_OTUS);
    sc.outStatistics = c.peek(SWARM_OUTPUT_STATISTICS);
    sc.outSeeds = c.peek(SWARM_OUTPUT_SEEDS);
    if (sc.outInternals) sc.oFileInternals = c.get(SWARM_OUTPUT_INTERNAL);
    if (sc.outOtus) sc.oFileOtus = c.get(SWARM_OUTPUT_OTUS);
    if (sc.outStatistics) sc.oFileStatistics = c.get(SWARM_OUTPUT_STATISTICS);
    if (sc.outSeeds) sc.oFileSeeds = c.get(SWARM_OUTPUT_SEEDS);
    if (c.peek(SWARM_NO_OTU_BREAKING)) sc.noOtuBreaking = (c.get(SWARM_NO_OTU_BREAKING) != "0");
    if (c.peek(SWARM_DEREPLICATE)) sc.dereplicate = (c.get(SWARM_DEREPLICATE) == "1");
    sc.sepAbundance = c.get(SEPARATOR_ABUNDANCE)[0];
    sc.extraSegs = std::stoul(c.get(NUM_EXTRA_SEGMENTS));

    if (sc.dereplicate) { // dereplication uses matching with distance 0
        c.set(THRESHOLD, "0");
    }

    sc.threshold = std::stoul(c.get(THRESHOLD));

    if (sc.dereplicate || (c.get(THRESHOLD) == "0")) { // fastidious option pointless when dereplicating or matching with distance 0

        c.set(SWARM_FASTIDIOUS, "0");
        sc.fastidious = false;

    }

    if (c.peek(SWARM_FASTIDIOUS)) sc.fastidious = (c.get(SWARM_FASTIDIOUS) == "1");
    if (sc.fastidious) sc.boundary = std::stoul(c.get(SWARM_BOUNDARY));


    //TODO remove or redirect to logger
    std::cout << "===== Configuration =====" << std::endl;
    c.print(std::cout);
    std::cout << "=========================" << std::endl;
    std::string jobName = c.get(NAME);
    std::replace(jobName.begin(), jobName.end(), ':', '-');
    if (c.peek(INFO_FOLDER)) writeJobParameters(c.get(INFO_FOLDER) + jobName + ".txt", c, files);


    /* Preprocessing */
    std::cout << "Preprocessing..." << std::flush;
//    auto prepStart = clock();
    auto pools = Preprocessor::run(c, files);
//    auto prepStop = clock();
    std::cout << "DONE" << std::endl;
//    std::cout << "Preprocessing: " << (prepStop - prepStart) / CLOCKS_PER_SEC << std::endl;

    /* Matching */
    std::cout << "Matching..." << std::flush;
//    auto matchStart = clock();
    unsigned long numWorkers = std::stoul(c.get(NUM_WORKERS));
    unsigned long numThreadsPerWorkers = std::stoul(c.get(NUM_THREADS_PER_WORKER));
    int mode = std::stoi(c.get(SEGMENT_FILTER));

    lenSeqs_t threshold = std::stoul(c.get(THRESHOLD));
    lenSeqs_t numExtraSegments = std::stoul(c.get(NUM_EXTRA_SEGMENTS));

    std::vector<Matches*> allMatches(pools->numPools()); // one Matches instance per pool
    std::vector<Subpool> subpools; // subpools of one pool
    std::thread workers[numWorkers]; // at most numWorkers per pool

    for (numSeqs_t p = 0; p < pools->numPools(); p++) {// std::cout << "Pool #" << p << std::endl; //TODO remove print

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
//    auto matchStop = clock();
    std::cout << "DONE" << std::endl;
//    std::cout << "Matching: " << (matchStop - matchStart) / CLOCKS_PER_SEC << std::endl;

    /* Postprocessing */
//    std::cout << "Postprocessing..." << std::endl;
//    auto postpStart = clock();
    if (c.peek(MATCHES_OUTPUT_FILE)) {

        std::cout << "Writing matches..." << std::flush;
        SegmentFilter::writeMatches(c.get(MATCHES_OUTPUT_FILE), *pools, allMatches);
        std::cout << "DONE" << std::endl;

    }

    std::cout << "Swarming results..." << std::endl;
    // parallel computation of swarm clusters & subsequent generation of outputs
    SwarmClustering::cluster(*pools, allMatches, sc);
    std::cout << "Swarming results...DONE" << std::endl;
//    auto postpStop = clock();
//    std::cout << "Postprocessing: " << (postpStop - postpStart) / CLOCKS_PER_SEC << std::endl;


    /* Cleaning up */
    std::cout << "Cleaning up..." << std::flush;
    for (auto iter = allMatches.begin(); iter != allMatches.end(); iter++) {
        delete *iter;
    }
    delete pools;
    std::cout << "DONE" << std::endl;

    return 0;

}

}


int main(int argc, const char* argv[]) {

    SCT_PJ::run(argc, argv);

}
