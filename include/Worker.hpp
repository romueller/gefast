/*
 * SCT-PJ
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

#ifndef SCT_PJ_WORKER_HPP
#define SCT_PJ_WORKER_HPP

#include <thread>

#include "Base.hpp"
#include "Relation.hpp"
#include "SegmentFilter.hpp"
#include "Utility.hpp"
#include "Verification.hpp"


namespace SCT_PJ {

class Worker {

public:
    Worker(AmpliconCollection& ac, Subpool& sp, Matches& matches, Config<std::string>& config);

    ~Worker();

    void run(const lenSeqs_t threshold, const lenSeqs_t numExtraSegments);


private:
    AmpliconCollection& ac_;
    Subpool& sp_;
    Matches& matches_;
    Config<std::string>& config_;

};

}

#endif //SCT_PJ_WORKER_HPP