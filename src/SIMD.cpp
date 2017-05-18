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

namespace SimdVerification {

    /* swarm.cc (adapted) */

    unsigned long threads;
    long penalty_gapextend;
    long penalty_gapopen;
    long penalty_mismatch;
    unsigned long diff_saturation;

    char sym_nt[] = "-acgt                           ";

    void args_init(lenSeqs_t maxLen, int penMismatch, int penOpen, int penExtend) {

        cpu_features_detect();

        longestdbsequence = maxLen;

        penalty_mismatch = penMismatch;
        penalty_gapopen = penOpen;
        penalty_gapextend = penExtend;

        diff_saturation = MIN(255 / penalty_mismatch, 255 / (penalty_gapopen + penalty_gapextend));

        score_matrix_init();

    }



    /* scan.cc (adapted) */

    static pthread_attr_t attr;

    static struct thread_info_s {

        /* generic thread info */
        pthread_t pthread;
        pthread_mutex_t workmutex;
        pthread_cond_t workcond;
        int work;

    } * ti;

    pthread_mutex_t workmutex = PTHREAD_MUTEX_INITIALIZER;

    queryinfo_t query;

    struct search_data {

        BYTE **qtable;
        WORD **qtable_w;

        BYTE *dprofile;
        WORD *dprofile_w;

        BYTE *hearray;

        unsigned long *dir_array;

        unsigned long target_count;
        unsigned long target_index;

    };

    struct search_data *sd;

    unsigned long master_next;
    unsigned long master_length;

    unsigned long remainingchunks;

    AmpliconCollection* master_pool;
    numSeqs_t *master_targets;
    lenSeqs_t *master_diffs;
    int master_bits;

    unsigned long longestdbsequence;
    unsigned long dirbufferbytes;

    void swarm_simd_search_alloc(struct search_data *sdp) {

        dirbufferbytes = 8 * longestdbsequence * ((longestdbsequence + 3) / 4) * 4;
        sdp->qtable = (BYTE **) xmalloc(longestdbsequence * sizeof(BYTE *));
        sdp->qtable_w = (WORD **) xmalloc(longestdbsequence * sizeof(WORD *));
        sdp->dprofile = (BYTE *) xmalloc(4 * 16 * 32);
        sdp->dprofile_w = (WORD *) xmalloc(4 * 2 * 8 * 32);
        sdp->hearray = (BYTE *) xmalloc(longestdbsequence * 32);
        sdp->dir_array = (unsigned long *) xmalloc(dirbufferbytes);

        memset(sdp->hearray, 0, longestdbsequence * 32);
        memset(sdp->dir_array, 0, dirbufferbytes);

    }

    void swarm_simd_search_free(struct search_data *sdp) {

        free(sdp->qtable);
        free(sdp->qtable_w);
        free(sdp->dprofile);
        free(sdp->dprofile_w);
        free(sdp->hearray);
        free(sdp->dir_array);

    }

    void swarm_simd_search_init(struct search_data *sdp) {

        for (long i = 0; i < query.len; i++) {
            sdp->qtable[i] = sdp->dprofile + 64 * query.seq[i];
            sdp->qtable_w[i] = sdp->dprofile_w + 32 * query.seq[i];
        }

    }


    void swarm_simd_search_chunk(struct search_data *sdp, long bits) {

        if (sdp->target_count == 0)
            return;

        if (bits == 16)
            search16_reduced(sdp->qtable_w,
                             penalty_gapopen,
                             penalty_gapextend,
                             (WORD *) score_matrix_16,
                             sdp->dprofile_w,
                             (WORD *) sdp->hearray,
                             sdp->target_count,
                             master_targets + sdp->target_index,
                             master_diffs + sdp->target_index,
                             query.len,
                             dirbufferbytes / 8,
                             sdp->dir_array,
                             *master_pool);
        else
            search8_reduced(sdp->qtable,
                            penalty_gapopen,
                            penalty_gapextend,
                            (BYTE *) score_matrix_8,
                            sdp->dprofile,
                            sdp->hearray,
                            sdp->target_count,
                            master_targets + sdp->target_index,
                            master_diffs + sdp->target_index,
                            query.len,
                            dirbufferbytes / 8,
                            sdp->dir_array,
                            *master_pool);

    }

    int swarm_simd_search_getwork(unsigned long *countref, unsigned long *firstref) {

        // * countref = how many sequences to search
        // * firstref = index into master_targets/scores/diffs where thread should start

        unsigned long status = 0;

        pthread_mutex_lock(&workmutex);

        if (master_next < master_length) {
            unsigned long chunksize =
                    ((master_length - master_next + remainingchunks - 1) / remainingchunks);

            *countref = chunksize;
            *firstref = master_next;

            master_next += chunksize;
            remainingchunks--;
            status = 1;
        }

        pthread_mutex_unlock(&workmutex);

        return status;

    }

    void verification_thread_worker_core(int t) {

        swarm_simd_search_init(sd + t);
        while (swarm_simd_search_getwork(&sd[t].target_count, &sd[t].target_index))
            swarm_simd_search_chunk(sd + t, master_bits);

    }

    void * verification_thread_worker(void *vp) {

        long t = (long) vp;
        struct thread_info_s *tip = ti + t;

        pthread_mutex_lock(&tip->workmutex);

        /* loop until signalled to quit */
        while (tip->work >= 0) {
            /* wait for work available */
            pthread_cond_wait(&tip->workcond, &tip->workmutex);
            if (tip->work > 0) {
                verification_thread_worker_core(t);
                tip->work = 0;
                pthread_cond_signal(&tip->workcond);
            }
        }
        pthread_mutex_unlock(&tip->workmutex);
        return 0;

    }

    void swarm_simd_init(int numThreads) {

        threads = numThreads;

        sd = (struct search_data *) xmalloc(sizeof(search_data) * threads);

        for (unsigned long t = 0; t < threads; t++)
            swarm_simd_search_alloc(sd + t);

        /* start threads */

        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

        /* allocate memory for thread info */
        ti = (struct thread_info_s *) xmalloc(threads *
                                              sizeof(struct thread_info_s));

        /* init and create worker threads */
        for (unsigned long t = 0; t < threads; t++) {
            struct thread_info_s *tip = ti + t;
            tip->work = 0;
            pthread_mutex_init(&tip->workmutex, NULL);
            pthread_cond_init(&tip->workcond, NULL);
            if (pthread_create(&tip->pthread, &attr, verification_thread_worker, (void *) (long) t))
                fatal("Cannot create thread");
        }

    }

    void swarm_simd_exit() {

        /* finish and clean up worker threads */

        for (unsigned long t = 0; t < threads; t++) {
            struct thread_info_s *tip = ti + t;

            /* tell worker to quit */
            pthread_mutex_lock(&tip->workmutex);
            tip->work = -1;
            pthread_cond_signal(&tip->workcond);
            pthread_mutex_unlock(&tip->workmutex);

            /* wait for worker to quit */
            if (pthread_join(tip->pthread, NULL))
                fatal("Cannot join thread");

            pthread_cond_destroy(&tip->workcond);
            pthread_mutex_destroy(&tip->workmutex);
        }

        free(ti);
        pthread_attr_destroy(&attr);

        for (unsigned long t = 0; t < threads; t++)
            swarm_simd_search_free(sd + t);
        free(sd);

    }

    void swarm_simd_verify(AmpliconCollection& ac, Amplicon& query_ampl, std::vector<numSeqs_t>* index_seqs, std::vector<lenSeqs_t>* diffs, long bits) {

        query.seq = (char*) query_ampl.seq.data();
        query.len = query_ampl.seq.length();

        master_next = 0;
        master_pool = &ac;
        master_length = index_seqs->size();
        master_targets = index_seqs->data();
        master_diffs = diffs->data();
        master_bits = bits;

        unsigned long thr = threads;

        if (bits == 8) {
            if (master_length <= (unsigned long) (15 * thr))
                thr = (master_length + 15) / 16;
        } else {
            if (master_length <= (unsigned long) (7 * thr))
                thr = (master_length + 7) / 8;
        }

        remainingchunks = thr;

        if (thr == 1) {
            verification_thread_worker_core(0);
        } else {
            /* wake up threads */
            for (unsigned long t = 0; t < thr; t++) {
                struct thread_info_s *tip = ti + t;
                pthread_mutex_lock(&tip->workmutex);
                tip->work = 1;
                pthread_cond_signal(&tip->workcond);
                pthread_mutex_unlock(&tip->workmutex);
            }

            /* wait for threads to finish their work */
            for (unsigned long t = 0; t < thr; t++) {
                struct thread_info_s *tip = ti + t;
                pthread_mutex_lock(&tip->workmutex);
                while (tip->work > 0)
                    pthread_cond_wait(&tip->workcond, &tip->workmutex);
                pthread_mutex_unlock(&tip->workmutex);
            }
        }

    }

    std::vector<std::pair<numSeqs_t, lenSeqs_t>> computeDiffs(AmpliconCollection& ac, Amplicon& query_ampl, std::vector<numSeqs_t>& index_seqs) {

        std::vector<lenSeqs_t> diffs(index_seqs.size());

        swarm_simd_verify(ac, query_ampl, &index_seqs, &diffs, 16);

        std::vector<std::pair<numSeqs_t, lenSeqs_t>> res(diffs.size());
        for (auto i = 0; i < index_seqs.size(); i++) {
            res[i] = std::make_pair(index_seqs[i], diffs[i]);
        }

        return res;

    }

    std::vector<std::pair<numSeqs_t, lenSeqs_t>> computeDiffsReduce(AmpliconCollection& ac, Amplicon& query_ampl, std::vector<numSeqs_t>& index_seqs, lenSeqs_t t) {

        std::vector<lenSeqs_t> diffs(index_seqs.size());

        swarm_simd_verify(ac, query_ampl, &index_seqs, &diffs, (t <= diff_saturation) ? 8 : 16);

        std::vector<std::pair<numSeqs_t, lenSeqs_t>> res;
        for (auto i = 0; i < index_seqs.size(); i++) {
            if (diffs[i] <= t) {
                res.push_back(std::make_pair(index_seqs[i], diffs[i]));
            }
        }

        return res;

    }

}

}