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

#include "../include/SIMD.hpp"

namespace GeFaST {

long mmx_present;
long sse_present;
long sse2_present;
long sse3_present;
long ssse3_present;
long sse41_present;
long sse42_present;
long popcnt_present;
long avx_present;
long avx2_present;

void cpu_features_detect() {
    unsigned int a, b, c, d;

    cpuid(0, 0, a, b, c, d);
    unsigned int maxlevel = a & 0xff;

    if (maxlevel >= 1) {
        cpuid(1, 0, a, b, c, d);
        mmx_present = (d >> 23) & 1;
        sse_present = (d >> 25) & 1;
        sse2_present = (d >> 26) & 1;
        sse3_present = (c >> 0) & 1;
        ssse3_present = (c >> 9) & 1;
        sse41_present = (c >> 19) & 1;
        sse42_present = (c >> 20) & 1;
        popcnt_present = (c >> 23) & 1;
        avx_present = (c >> 28) & 1;

        if (maxlevel >= 7) {
            cpuid(7, 0, a, b, c, d);
            avx2_present = (b >> 5) & 1;
        }
    }
}


#if QGRAM_FILTER

unsigned long popcount_128(__m128i x) {
    __m128i mask1 = _mm_set_epi8(0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55,
                                 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55);

    __m128i mask2 = _mm_set_epi8(0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33,
                                 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33);

    __m128i mask4 = _mm_set_epi8(0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f,
                                 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f);

    __m128i zero = _mm_setzero_si128();

    /* add together 2 bits: 0+1, 2+3, 3+4, ... 126+127 */

    __m128i a = _mm_srli_epi64(x, 1);
    __m128i b = _mm_and_si128(x, mask1);
    __m128i c = _mm_and_si128(a, mask1);
    __m128i d = _mm_add_epi64(b, c);

    /* add together 4 bits: (0+1)+(2+3), ... (124+125)+(126+127) */

    __m128i e = _mm_srli_epi64(d, 2);
    __m128i f = _mm_and_si128(d, mask2);
    __m128i g = _mm_and_si128(e, mask2);
    __m128i h = _mm_add_epi64(f, g);

    /* add together 8 bits: (0..3)+(4..7), ... (120..123)+(124..127) */

    __m128i i = _mm_srli_epi64(h, 4);
    __m128i j = _mm_add_epi64(h, i);
    __m128i k = _mm_and_si128(j, mask4);

    /* add together 8 bytes: (0..63) and (64..127) */

    __m128i l = _mm_sad_epu8(k, zero);

    /* add together 64-bit values into final 128 bit value */

    __m128i m = _mm_srli_si128(l, 8);
    __m128i n = _mm_add_epi64(m, l);

    /* return low 64 bits: return value is always in range 0 to 128 */

    unsigned long o = (unsigned long) _mm_movepi64_pi64(n);

    return o;
}

unsigned long compareqgramvectors_128(const unsigned char * a, const unsigned char * b) {
    /* Count number of different bits */
    /* Uses SSE2 but not POPCNT instruction */
    /* input MUST be 16-byte aligned */

    __m128i * ap = (__m128i *) a;
    __m128i * bp = (__m128i *) b;
    unsigned long count = 0;

    while ((unsigned char*)ap < a + QGRAMVECTORBYTES)
        count += popcount_128(_mm_xor_si128(*ap++, *bp++));

    return count;
}

unsigned long compareqgramvectors_64(const unsigned char * a, const unsigned char * b) {
    /* Count number of different bits */
    /* Uses the POPCNT instruction, requires CPU with this feature */

    unsigned long *ap = (unsigned long*)a;
    unsigned long *bp = (unsigned long*)b;
    unsigned long count = 0;

    while ((unsigned char*) ap < a + QGRAMVECTORBYTES)
        count += popcount(*ap++ ^ *bp++);

    return count;
}

unsigned long compareqgramvectors(const unsigned char * a, const unsigned char * b) {
    if (popcnt_present)
        return compareqgramvectors_64(a,b);
    else
        return compareqgramvectors_128(a,b);
}

#endif

}