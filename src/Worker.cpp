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

#include "../include/Worker.hpp"


namespace SCT_PJ {

Worker::Worker(AmpliconCollection& ac, Subpool& sp, Matches& matches, unsigned long numThreads, int mode) : ac_(ac), sp_(sp), matches_(matches) {

    mode_ = mode;
    numThreads_ = numThreads;

}

Worker::~Worker() {
    // nothing to do
}

void Worker::run(const lenSeqs_t threshold, const lenSeqs_t numExtraSegments) {

    RotatingBuffers<Candidate> cbs;

    if (numThreads_ == 1) { // completey sequential

        cbs = RotatingBuffers<Candidate>(1);
        SegmentFilter::filter(ac_, sp_, cbs, threshold, numExtraSegments, mode_);
        cbs.close();

        Verification::verify(ac_, matches_, cbs.getBuffer(0), threshold);

    } else { // worker thread + (numThreads - 1) verifier threads

        cbs = RotatingBuffers<Candidate>(numThreads_ - 1);

        Matches* matchesPerThread[numThreads_ - 1];
        matchesPerThread[0] = &matches_;
        for (unsigned long v = 1; v < numThreads_ - 1; v++) {
            matchesPerThread[v] = new Matches();
        }

        std::thread verifierThreads[numThreads_ - 1];
        for (unsigned long v = 0; v < numThreads_ - 1; v++) {
            verifierThreads[v] = std::thread(&Verification::verify, std::ref(ac_), std::ref(*(matchesPerThread[v])), std::ref(cbs.getBuffer(v)), threshold);
        }

        SegmentFilter::filter(ac_, sp_, cbs, threshold, numExtraSegments, mode_);
        cbs.close();

        for (unsigned long v = 0; v < numThreads_ - 1; v++) {

            verifierThreads[v].join();

            if (v > 0) {

                matches_.join(*(matchesPerThread[v]));
                delete matchesPerThread[v];

            }

        }

    }

}

}