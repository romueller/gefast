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

#ifdef __APPLE__
#include <sys/sysctl.h>
#else
#include <sys/sysinfo.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifdef __SSSE3__
#include <tmmintrin.h>
#endif

#include "Base.hpp"

#ifndef GEFAST_SIMD_HPP
#define GEFAST_SIMD_HPP

namespace GeFaST {

// Adapted from:
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

// =====================================================
//              CPU feature detection
// =====================================================

extern long mmx_present;
extern long sse_present;
extern long sse2_present;
extern long sse3_present;
extern long ssse3_present;
extern long sse41_present;
extern long sse42_present;
extern long popcnt_present;
extern long avx_present;
extern long avx2_present;

#define cpuid(f1, f2, a, b, c, d)                                       \
__asm__ __volatile__ ("cpuid"                                         \
                    : "=a" (a), "=b" (b), "=c" (c), "=d" (d)        \
                    : "a" (f1), "c" (f2));

void cpu_features_detect();



// =====================================================
//                  q-gram filter
// =====================================================

#if QGRAM_FILTER

#define popcnt_asm(x, y)                                         \
__asm__ __volatile__ ("popcnt %1,%0" : "=r"(y) : "r"(x));

inline unsigned long popcount(unsigned long x) {
    unsigned long y;
    popcnt_asm(x, y);
    return y;
}

unsigned long popcount_128(__m128i x);

unsigned long compareqgramvectors_128(const unsigned char * a, const unsigned char * b);

unsigned long compareqgramvectors_64(const unsigned char * a, const unsigned char * b);

unsigned long compareqgramvectors(const unsigned char * a, const unsigned char * b);

inline unsigned long qgram_diff(const Amplicon &a, const Amplicon &b) {
    unsigned long diffqgrams = compareqgramvectors(a.qGramVector, b.qGramVector);
    unsigned long mindiff = (diffqgrams + 2 * QGRAMLENGTH - 1) / (2 * QGRAMLENGTH);
    return mindiff;
}

#endif

}

#endif //GEFAST_SIMD_HPP
