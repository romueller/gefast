/*
    SWARM

    Copyright (C) 2012-2017 Torbjorn Rognes and Frederic Mahe

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include "../../include/SIMD.hpp"

namespace GeFaST {
namespace SimdVerification {

    void fatal(const char *msg) {
        fprintf(stderr, "\nError: %s\n", msg);
        exit(1);
    }

    void *xmalloc(size_t size) {
        const size_t alignment = 16;
        void *t = NULL;
        if (posix_memalign(&t, alignment, size))
            fatal("Unable to allocate enough memory.");
        if (!t)
            fatal("Unable to allocate enough memory.");
        return t;
    }

}
}