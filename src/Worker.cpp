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

#include "../include/Worker.hpp"


namespace GeFaST {

Worker::Worker(AmpliconCollection& ac, Subpool& sp, Matches& matches, Config<std::string>& config) : ac_(ac), sp_(sp), matches_(matches), config_(config) {
    // nothing else to do
}

Worker::~Worker() {
    // nothing to do
}

void Worker::run(const lenSeqs_t threshold, const lenSeqs_t numExtraSegments) {

    unsigned long numThreads = std::stoul(config_.get(NUM_THREADS_PER_WORKER));
    int mode = std::stoi(config_.get(SEGMENT_FILTER));
    bool useScore = (config_.get(USE_SCORE) == "1");
    Verification::Scoring scoring(std::stoull(config_.get(SWARM_MATCH_REWARD)), std::stoll(config_.get(SWARM_MISMATCH_PENALTY)), std::stoll(config_.get(SWARM_GAP_OPENING_PENALTY)), std::stoll(config_.get(SWARM_GAP_EXTENSION_PENALTY)));

    if (numThreads == 1) { // completely sequential

        Matches localMatches;

#if 0

        RotatingBuffers<Candidate> cbs = RotatingBuffers<Candidate>(1);
        SegmentFilter::filter(ac_, sp_, cbs, threshold, numExtraSegments, mode);
        cbs.close();

        if (useScore) {
            Verification::verifyGotoh(ac_, localMatches, cbs.getBuffer(0), ac_.back().seq.length() + 1, threshold, scoring);
        } else {
            Verification::verify(ac_, localMatches, cbs.getBuffer(0), ac_.back().seq.length() + 1, threshold);
        }

#else

        SegmentFilter::filterDirectly(ac_, sp_, localMatches, threshold, numExtraSegments, mode, useScore, scoring);

#endif

        matches_.syncJoin(localMatches);

    } else { // worker thread + (numThreads - 1) verifier threads

        RotatingBuffers<Candidate> cbs = RotatingBuffers<Candidate>(numThreads - 1);

        Matches* matchesPerThread[numThreads - 1];
        for (unsigned long v = 0; v < numThreads - 1; v++) {
            matchesPerThread[v] = new Matches();
        }

        std::thread verifierThreads[numThreads - 1];
        for (unsigned long v = 0; v < numThreads - 1; v++) {
            verifierThreads[v] = useScore ?
                                   std::thread(&Verification::verifyGotoh, std::ref(ac_), std::ref(*(matchesPerThread[v])), std::ref(cbs.getBuffer(v)), ac_.back().seq.length() + 1, threshold, std::ref(scoring))
                                 : std::thread(&Verification::verify, std::ref(ac_), std::ref(*(matchesPerThread[v])), std::ref(cbs.getBuffer(v)), ac_.back().seq.length() + 1, threshold);
        }

        SegmentFilter::filter(ac_, sp_, cbs, threshold, numExtraSegments, mode);
        cbs.close();

        for (unsigned long v = 0; v < numThreads - 1; v++) {

            verifierThreads[v].join();

            matches_.syncJoin(*(matchesPerThread[v]));
            delete matchesPerThread[v];

        }

    }

}

}