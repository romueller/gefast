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

#include <ctime>
#include <iostream>
#include <vector>
#include <thread>

#include "include/Base.hpp"
#include "include/Preprocessor.hpp"
#include "include/Relation.hpp"
#include "include/SIMD.hpp"
#include "include/SwarmClustering.hpp"
#include "include/Utility.hpp"
#include "include/Worker.hpp"







namespace GeFaST {

int run(int argc, const char* argv[]) {

    /* Bootstrapping */
    Config<std::string> c = getConfiguration(argc, argv);
#if QGRAM_FILTER
    cpu_features_detect();
#endif


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
    if (!((c.get(PREPROCESSING_ONLY) == "1") || c.peek(MATCHES_OUTPUT_FILE) || c.peek(SWARM_OUTPUT_INTERNAL) || c.peek(SWARM_OUTPUT_OTUS) || c.peek(SWARM_OUTPUT_STATISTICS) || c.peek(SWARM_OUTPUT_SEEDS) || c.peek(SWARM_OUTPUT_UCLUST))) {

        std::cerr << "ERROR: No output file specified." << std::endl;
        return 1;

    }

    SwarmClustering::SwarmConfig sc;
    sc.outInternals = c.peek(SWARM_OUTPUT_INTERNAL);
    sc.outOtus = c.peek(SWARM_OUTPUT_OTUS);
    sc.outMothur = (c.peek(SWARM_MOTHUR) && (c.get(SWARM_MOTHUR) == "1"));
    sc.outStatistics = c.peek(SWARM_OUTPUT_STATISTICS);
    sc.outSeeds = c.peek(SWARM_OUTPUT_SEEDS);
    sc.outUclust = c.peek(SWARM_OUTPUT_UCLUST);
    if (sc.outInternals) sc.oFileInternals = c.get(SWARM_OUTPUT_INTERNAL);
    if (sc.outOtus) sc.oFileOtus = c.get(SWARM_OUTPUT_OTUS);
    if (sc.outStatistics) sc.oFileStatistics = c.get(SWARM_OUTPUT_STATISTICS);
    if (sc.outSeeds) sc.oFileSeeds = c.get(SWARM_OUTPUT_SEEDS);
    if (sc.outUclust) sc.oFileUclust = c.get(SWARM_OUTPUT_UCLUST);
    sc.noOtuBreaking = (c.get(SWARM_NO_OTU_BREAKING) != "0");
    sc.dereplicate = (c.get(SWARM_DEREPLICATE) == "1");
    sc.sepAbundance = c.get(SEPARATOR_ABUNDANCE);
    sc.extraSegs = std::stoul(c.get(NUM_EXTRA_SEGMENTS));
    sc.filterTwoWay = (std::stoul(c.get(SEGMENT_FILTER)) == 2) || (std::stoul(c.get(SEGMENT_FILTER)) == 3);
    sc.numExplorers = std::stoul(c.get(SWARM_NUM_EXPLORERS));
    sc.numThreadsPerExplorer = std::stoul(c.get(SWARM_NUM_THREADS_PER_CHECK));
    sc.numGrafters = std::stoul(c.get(SWARM_NUM_GRAFTERS));
    sc.fastidiousCheckingMode = std::stoul(c.get(SWARM_FASTIDIOUS_CHECKING_MODE));
    sc.numThreadsPerCheck = std::stoul(c.get(SWARM_NUM_THREADS_PER_CHECK));
    sc.threshold = std::stoul(c.get(THRESHOLD));

    if (sc.threshold <= 0 && !sc.dereplicate) { // check for feasible threshold unless dereplication is chosen

        std::cerr << "ERROR: Only positive thresholds are allowed." << std::endl;
        return 1;

    }

    if (c.get(SWARM_FASTIDIOUS_THRESHOLD) == "0") { // "default" corresponds to swarm postulating one virtual linking amplicon

        c.set(SWARM_FASTIDIOUS_THRESHOLD, std::to_string(2 * sc.threshold));
        sc.fastidiousThreshold = 2 * sc.threshold;

    } else {
        sc.fastidiousThreshold = std::stoul(c.get(SWARM_FASTIDIOUS_THRESHOLD));
    }

    if (sc.dereplicate) { // fastidious option pointless when dereplicating or matching with distance 0

        c.set(SWARM_FASTIDIOUS, "0");
        sc.fastidious = false;

        c.set(USE_SCORE, "0");
        sc.useScore = false;

    }

    sc.fastidious = (c.get(SWARM_FASTIDIOUS) == "1");
    sc.boundary = std::stoul(c.get(SWARM_BOUNDARY));

    sc.useScore = (c.get(USE_SCORE) == "1");
    sc.scoring = Verification::Scoring(std::stoull(c.get(SWARM_MATCH_REWARD)), std::stoll(c.get(SWARM_MISMATCH_PENALTY)), std::stoll(c.get(SWARM_GAP_OPENING_PENALTY)), std::stoll(c.get(SWARM_GAP_EXTENSION_PENALTY)));

    //TODO remove or redirect to logger
    std::cout << "===== Configuration =====" << std::endl;
    c.print(std::cout);
    std::cout << "=========================" << std::endl;
    std::string jobName = c.get(NAME);
    std::replace(jobName.begin(), jobName.end(), ':', '-');
    if (c.peek(INFO_FOLDER)) writeJobParameters(c.get(INFO_FOLDER) + jobName + ".txt", c, files);


    /* Preprocessing */
    std::cout << "Preprocessing..." << std::flush;
    auto pools = Preprocessor::run(c, files);
#if 1 // full index, direct swarming
    AmpliconCollection* ac;
    for (numSeqs_t p = 0; p < pools->numPools(); p++) {

        ac = pools->get(p);
        std::sort(ac->begin(), ac->end(),
             [](const Amplicon& amplA, const Amplicon& amplB) {
                 return (amplA.abundance > amplB.abundance) || ((amplA.abundance == amplB.abundance) && (amplA.seq < amplB.seq));
             }
        );

    }

    std::cout << "DONE" << std::endl;

    if (c.get(PREPROCESSING_ONLY) == "1") {

        std::cout << "Cleaning up..." << std::flush;
        delete pools;
        std::cout << "DONE" << std::endl;

        return 0;

    }

    if (sc.dereplicate) { /* Dereplication */

        std::cout << "Dereplicating..." << std::endl;
        SwarmClustering::dereplicate(*pools, sc);
        std::cout << "Dereplicating...DONE" << std::endl;

    } else { /* Swarming */

        std::cout << "Swarming results..." << std::endl;
        // parallel computation of swarm clusters & subsequent generation of outputs
        SwarmClustering::cluster(*pools, sc);
        std::cout << "Swarming results...DONE" << std::endl;

    }

    /* Cleaning up */
    std::cout << "Cleaning up..." << std::flush;
    delete pools;
    std::cout << "DONE" << std::endl;
#else // iterative index, swarming via matches
    std::cout << "DONE" << std::endl;

    AmpliconCollection* ac;
    for (numSeqs_t p = 0; p < pools->numPools(); p++) {

        ac = pools->get(p);
        std::sort(ac->begin(), ac->end(),
                  [](const Amplicon& amplA, const Amplicon& amplB) {
                      return (amplA.len < amplB.len) || ((amplA.len == amplB.len) && (strcmp(amplA.seq, amplB.seq) < 0));
                  }
        );

    }

    if (c.get(PREPROCESSING_ONLY) == "1") {

        std::cout << "Cleaning up..." << std::flush;
        delete pools;
        std::cout << "DONE" << std::endl;

        return 0;

    }

    /* Matching */
    std::cout << "Matching..." << std::flush;
    unsigned long numWorkers = std::stoul(c.get(NUM_WORKERS));
    int mode = std::stoi(c.get(SEGMENT_FILTER));

    lenSeqs_t threshold = std::stoul(c.get(THRESHOLD));
    lenSeqs_t numExtraSegments = std::stoul(c.get(NUM_EXTRA_SEGMENTS));

    std::vector<Matches*> allMatches(pools->numPools()); // one Matches instance per pool
    std::vector<Subpool> subpools; // subpools of one pool
    std::thread workers[numWorkers]; // at most numWorkers per pool

    for (numSeqs_t p = 0; p < pools->numPools(); p++) {

        allMatches[p] = new Matches();
        if (mode % 2 == 0) {
            subpools = getSubpoolBoundaries(*(pools->get(p)), numWorkers, threshold);
        } else {
            subpools = getSubpoolBoundariesBackward(*(pools->get(p)), numWorkers, threshold);
        }

        for (unsigned long w = 0; w < subpools.size(); w++) {
            workers[w] = std::thread(&Worker::run, Worker(*(pools->get(p)), subpools[w], *(allMatches[p]), c), threshold, numExtraSegments);
        }

        for (unsigned long w = 0; w < subpools.size(); w++) {
            workers[w].join();
        }

    }

    std::cout << "DONE" << std::endl;

    /* Postprocessing */
    if (c.peek(MATCHES_OUTPUT_FILE)) {

        std::cout << "Writing matches..." << std::flush;
        SegmentFilter::writeMatches(c.get(MATCHES_OUTPUT_FILE), *pools, allMatches);
        std::cout << "DONE" << std::endl;

    }

    std::cout << "Swarming results..." << std::endl;
    // parallel computation of swarm clusters & subsequent generation of outputs
    SwarmClustering::cluster(*pools, allMatches, sc);
    std::cout << "Swarming results...DONE" << std::endl;


    /* Cleaning up */
    std::cout << "Cleaning up..." << std::flush;
    for (auto iter = allMatches.begin(); iter != allMatches.end(); iter++) {
        delete *iter;
    }
    delete pools;
    std::cout << "DONE" << std::endl;
#endif

    return 0;

}

}


int main(int argc, const char* argv[]) {

    GeFaST::run(argc, argv);

}
