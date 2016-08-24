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

#ifndef SCT_PJ_WORKER_HPP
#define SCT_PJ_WORKER_HPP

#include <thread>

#include "Base.hpp"
#include "Relation.hpp"
#include "SegmentFilter.hpp"
#include "Verification.hpp"


namespace SCT_PJ {

class Worker {

public:
    Worker(AmpliconCollection& ac, Subpool& sp, Matches& matches, unsigned long numThreads, int mode);

    ~Worker();

    void run(const lenSeqs_t threshold, const lenSeqs_t numExtraSegments);


private:
    AmpliconCollection& ac_;
    Subpool& sp_;
    Matches& matches_;
    unsigned long numThreads_;
    int mode_;

};

}

#endif //SCT_PJ_WORKER_HPP