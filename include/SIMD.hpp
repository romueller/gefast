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



// =====================================================
//                  verification
// =====================================================

namespace SimdVerification {

/* constants */

#ifndef LINE_MAX
#define LINE_MAX 2048
#endif

#ifndef MIN
#define MIN(x, y) ((x)<(y)?(x):(y))
#endif

    /* structures and data types */

    typedef unsigned short WORD;
    typedef unsigned char BYTE;

    struct queryinfo {
        long len;
        char *seq;
    };

    typedef struct queryinfo queryinfo_t;

    /* common data */

    extern unsigned long threads;
    extern long penalty_gapextend;
    extern long penalty_gapopen;
    extern long penalty_mismatch;
    extern unsigned long diff_saturation;

    extern unsigned char *score_matrix_8;
    extern unsigned short *score_matrix_16;

    extern char sym_nt[];

    extern unsigned long longestdbsequence;

    extern queryinfo_t query;


    /* functions in util.cc */

    void fatal(const char *msg);

    void *xmalloc(size_t size);


    /* functions in ssse3.cc */

    void dprofile_shuffle8(BYTE *dprofile,
                           BYTE *score_matrix,
                           BYTE *dseq_byte);

    void dprofile_shuffle16(WORD *dprofile,
                            WORD *score_matrix,
                            BYTE *dseq_byte);


    /* functions in search8.cc */

    void search8_reduced(BYTE **q_start,
                         BYTE gap_open_penalty,
                         BYTE gap_extend_penalty,
                         BYTE *score_matrix,
                         BYTE *dprofile,
                         BYTE *hearray,
                         unsigned long sequences,
                         numSeqs_t *seqnos,
                         lenSeqs_t *diffs,
                         unsigned long qlen,
                         unsigned long dirbuffersize,
                         unsigned long *dirbuffer,
                         AmpliconCollection& ac);


    /* functions in search16.cc */

    void search16_reduced(WORD **q_start,
                          WORD gap_open_penalty,
                          WORD gap_extend_penalty,
                          WORD *score_matrix,
                          WORD *dprofile,
                          WORD *hearray,
                          unsigned long sequences,
                          numSeqs_t *seqnos,
                          lenSeqs_t *diffs,
                          unsigned long qlen,
                          unsigned long dirbuffersize,
                          unsigned long *dirbuffer,
                          AmpliconCollection& ac);


    /* functions in matrix.cc */

    void score_matrix_init();

    void score_matrix_free();


    /* functions in SIMD.cpp */

    void swarm_simd_init(int numThreads);

    void swarm_simd_exit();

    void swarm_simd_verify(AmpliconCollection& ac, Amplicon& query_ampl, std::vector<numSeqs_t>* index_seqs, std::vector<lenSeqs_t>* diffs, long bits);

    std::vector<std::pair<numSeqs_t, lenSeqs_t>> computeDiffs(AmpliconCollection& ac, Amplicon& query_ampl, std::vector<numSeqs_t>& index_seqs);
    std::vector<std::pair<numSeqs_t, lenSeqs_t>> computeDiffsReduce(AmpliconCollection& ac, Amplicon& query_ampl, std::vector<numSeqs_t>& index_seqs, lenSeqs_t t);

    void args_init(lenSeqs_t, int penMismatch, int penOpen, int penExtend);

}

}

#endif //GEFAST_SIMD_HPP
