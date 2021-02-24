/*
 * GeFaST
 *
 * Copyright (C) 2016 - 2021 Robert Mueller
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

#include <sstream>

#include "../include/Distances.hpp"

namespace GeFaST {

    /* === ScoringFunction === */

    ScoringFunction::ScoringFunction(int match, int mismatch, int gap_open, int gap_extend) {

        // convert scores
        penalty_mismatch = 2 * (match - mismatch);
        penalty_open = -2 * gap_open;
        penalty_extend = match - 2 * gap_extend;

        auto penalty_factor = gcd(gcd(penalty_mismatch, penalty_open), penalty_extend);

        penalty_mismatch /= penalty_factor;
        penalty_open /= penalty_factor;
        penalty_extend /= penalty_factor;

    }


    /* === AlignmentInformation === */

    AlignmentInformation::AlignmentInformation(const std::string& c, const lenSeqs_t l, const lenSeqs_t n/*, const lenSeqs_t s*/) {

        cigar = c;
        length = l;
        num_diffs = n;
//        score = s;

    }

    // "backtracking flags" for Gotoh's algorithm using three matrices
    const char DIAGONAL_IN_D = 1;
    const char JUMP_TO_Q = 2;
    const char JUMP_TO_P = 4;
    const char UP_TO_D = 8;
    const char UP_IN_P = 16;
    const char LEFT_TO_D = 32;
    const char LEFT_IN_Q = 64;

    const char IN_D = 0;
    const char IN_P = 1;
    const char IN_Q = 2;

    AlignmentInformation compute_gotoh_cigar_row(const char* s, const lenSeqs_t len_s, const char* t,
            const lenSeqs_t len_t, const ScoringFunction& scoring, lenSeqs_t* d_row, lenSeqs_t* p_row, char* backtrack) {

        lenSeqs_t width = len_t + 1;

        // initialise first row
        d_row[0] = scoring.penalty_open;
        for (lenSeqs_t j = 1; j <= len_t; j++) {

            d_row[j] = d_row[j - 1] + scoring.penalty_extend;
            p_row[j] = POS_INF;

        }
        d_row[0] = 0;

        // compute remaining rows
        lenSeqs_t match, val_q, from_d, from_pq, min_val;
        char tmp;

        for (lenSeqs_t i = 1; i <= len_s; i++) {

            // handle left end
            d_row[0] = scoring.penalty_open + i * scoring.penalty_extend;
            match = (i > 1) * (d_row[0] - scoring.penalty_extend);
            val_q = POS_INF;

            // fill remaining row
            for (lenSeqs_t j = 1; j <= len_t; j++) {

                // array P
                from_d = d_row[j] + scoring.penalty_open + scoring.penalty_extend;
                from_pq = p_row[j] + scoring.penalty_extend;

                if (from_d <= from_pq) {

                    p_row[j] = from_d;
                    tmp = UP_TO_D;

                } else {

                    p_row[j] = from_pq;
                    tmp = UP_IN_P;

                }

                // array Q
                from_d = d_row[j - 1] + scoring.penalty_open + scoring.penalty_extend;
                from_pq = val_q + scoring.penalty_extend;

                if (from_d <= from_pq) {

                    val_q = from_d;
                    tmp |= LEFT_TO_D;

                } else {

                    val_q = from_pq;
                    tmp |= LEFT_IN_Q;

                }

                // arrays D & backtrack
                min_val = match + (s[i - 1] != t[j - 1]) * scoring.penalty_mismatch;
                backtrack[i * width + j] = DIAGONAL_IN_D;
                if (p_row[j] < min_val) {

                    min_val = p_row[j];
                    backtrack[i * width + j] = JUMP_TO_P;

                }
                if (val_q <= min_val){

                    min_val = val_q;
                    backtrack[i * width + j] = JUMP_TO_Q;

                }

                backtrack[i * width + j] |= tmp;

                match = d_row[j];
                d_row[j] = min_val;

            }

        }

        // backtracking (priorities: left > diagonal > up)
        std::vector<std::pair<char, lenSeqs_t >> cigarSegs = {std::make_pair('N',0)};
        lenSeqs_t len = 0;
        lenSeqs_t num_diffs = 0;
        lenSeqs_t i = len_s;
        lenSeqs_t j = len_t;

        while (i != 0 && j != 0) {

            if ((cigarSegs.back().first == 'I') && (backtrack[i * width + j + 1] & LEFT_IN_Q)) {

                len++;
                num_diffs++;
                cigarSegs.back().second++;

                j--;

            } else if ((cigarSegs.back().first == 'D') && (backtrack[(i + 1) * width + j] & UP_IN_P)) {

                len++;
                num_diffs++;
                cigarSegs.back().second++;

                i--;

            } else if (backtrack[i * width + j] & JUMP_TO_Q) {

                len++;
                num_diffs++;
                if (cigarSegs.back().first == 'I') {
                    cigarSegs.back().second++;
                } else {
                    cigarSegs.emplace_back('I', 1);
                }

                j--;

            } else if (backtrack[i * width + j] & DIAGONAL_IN_D) {

                len++;
                num_diffs += (s[i - 1] != t[j - 1]);
                if (cigarSegs.back().first == 'M') {
                    cigarSegs.back().second++;
                } else {
                    cigarSegs.emplace_back('M', 1);
                }

                i--;
                j--;

            } else { // backtrack[i][j] & JUMP_TO_P

                len++;
                num_diffs++;
                if (cigarSegs.back().first == 'D') {
                    cigarSegs.back().second++;
                } else {
                    cigarSegs.emplace_back('D', 1);
                }

                i--;

            }

        }

        while (i > 0) {

            len++;
            num_diffs++;
            if (cigarSegs.back().first == 'D') {
                cigarSegs.back().second++;
            } else {
                cigarSegs.emplace_back('D', 1);
            }

            i--;

        }

        while (j > 0) {

            len++;
            num_diffs++;
            if (cigarSegs.back().first == 'I') {
                cigarSegs.back().second++;
            } else {
                cigarSegs.emplace_back('I', 1);
            }

            j--;

        }

        std::stringstream str_stream;
        for (auto p = cigarSegs.size() - 1; p > 0; p--) {

            if (cigarSegs[p].second > 1) str_stream << cigarSegs[p].second;
            str_stream << cigarSegs[p].first;

        }

        return AlignmentInformation(str_stream.str(), len, num_diffs);

    }


    /* === BoundedLevenshteinDistance === */

    BoundedLevenshteinDistance::BoundedLevenshteinDistance(lenSeqs_t max_seq_len, lenSeqs_t threshold) {

        row_length_ = max_seq_len + 1;
        matrix_row_ = new lenSeqs_t[row_length_];
        threshold_ = threshold;

    }

    BoundedLevenshteinDistance::~BoundedLevenshteinDistance() {
        delete[] matrix_row_;
    }

    BoundedLevenshteinDistance::BoundedLevenshteinDistance(const BoundedLevenshteinDistance& other) : // copy constructor
            threshold_(other.threshold_), row_length_(other.row_length_) {

        matrix_row_ = new lenSeqs_t[row_length_];
        for (lenSeqs_t i = 0; i < row_length_; i++) {
            matrix_row_[i] = other.matrix_row_[i];
        }

    }

    BoundedLevenshteinDistance* BoundedLevenshteinDistance::clone() const {
        return new BoundedLevenshteinDistance(*this);
    }

    dist_t BoundedLevenshteinDistance::distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) {
        return compute_length_aware_row(ac.seq(i), ac.len(i), ac.seq(j), ac.len(j));
    }

    lenSeqs_t BoundedLevenshteinDistance::compute_length_aware_row(const char* s, const lenSeqs_t len_s,
            const char* t, const lenSeqs_t len_t) {

        // long computation not necessary if lengths differ too much
        if (((len_s > len_t) ? (len_s - len_t) : (len_t - len_s)) > threshold_) {
            return threshold_ + 1;
        }

        if (threshold_ == 0) {
            return lenSeqs_t(s != t);
        }


        const char* shorter = (len_s < len_t) ? s : t;
        lenSeqs_t len_shorter = std::min(len_s, len_t);
        const char* longer = (len_s >= len_t) ? s : t;
        lenSeqs_t len_longer = std::max(len_s, len_t);
        lenSeqs_t diff = len_longer - len_shorter;


        // instead of storing the full matrix M, only a single row R and two helper variables are used
        // (assume we want to compute the entry M[i][j]):
        // - match: stores the score M[i - 1, j - 1] from the preceding row, its position in R was overwritten in the previous step
        // - tmp: temporarily contains newly computed score for M[i][j] in order to move M[i - 1][j] (stored in R[j]) to match
        //   before reassigning the position
        lenSeqs_t match, tmp;

        // (mis)match is only possibility when we have to consider only one diagonal
        // [happens only if (a) bound = delta = 0, or (b) bound = 1 and delta = 0, but (a) is already covered above]
        if ((threshold_ - diff) / 2 == 0 && (threshold_ + diff) / 2 == 0) {

            lenSeqs_t diffs = 0;
            for (auto i = 0; diffs <= threshold_ && i < len_shorter; i++) {
                diffs += (shorter[i] != longer[i]);
            }

            return diffs;

        }

        // initialise necessary sections of first row
        for (lenSeqs_t j = 0; j <= (threshold_ + diff) / 2 && j <= len_longer; j++) {
            matrix_row_[j] = j;
        }

        lenSeqs_t j;
        bool early;

        // compute sections of remaining rows
        for (lenSeqs_t i = 1; i <= len_shorter; i++) {

            early = true; // early termination flag

            j = 1 + (i > (threshold_ - diff) / 2) * (i - (threshold_ - diff) / 2 - 1); // same as j = max(1, i - (threshold_ - diff) / 2) with signed integers
            match = matrix_row_[j - 1];
            matrix_row_[j - 1] = (i <= (threshold_ - diff) / 2) ? i : POS_INF; // handle left end to avoid case distinction
            if (i + (threshold_ + diff) / 2 <= len_longer) { // handle right end to avoid case distinction
                matrix_row_[i + (threshold_ + diff) / 2] = POS_INF;
            }

            for (; j <= i + (threshold_ + diff) / 2 && j <= len_longer; j++) {

                tmp = std::min({
                                       match + (shorter[i - 1] != longer[j - 1]), // (mis)match
                                       matrix_row_[j] + 1, // deletion
                                       matrix_row_[j - 1] + 1 // insertion
                               });

                match = matrix_row_[j];
                matrix_row_[j] = tmp;

                early &= ((matrix_row_[j] + llabs((long long)diff + (long long)i - (long long)j)) > threshold_); // improved e.t.

            }

            if (early) { // computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
                return threshold_ + 1;
            }


        }

        return (matrix_row_[len_longer] > threshold_) ? (threshold_ + 1) : matrix_row_[len_longer];

    }



    /* === BoundedScoreLevenshteinDistance === */

    BoundedScoreLevenshteinDistance::BoundedScoreLevenshteinDistance(lenSeqs_t max_seq_len, lenSeqs_t threshold, int match,
            int mismatch, int gap_open, int gap_extend) :
            scoring_(ScoringFunction(match, mismatch, gap_open, gap_extend)) {

        row_length_ = max_seq_len + 1;
        d_row_ = new lenSeqs_t[row_length_];
        p_row_ = new lenSeqs_t[row_length_];
        cnt_diffs_row_ = new lenSeqs_t[row_length_];
        cnt_diffs_p_row_ = new lenSeqs_t[row_length_];
        threshold_ = threshold;

    }

    BoundedScoreLevenshteinDistance::~BoundedScoreLevenshteinDistance() {

        delete[] cnt_diffs_p_row_;
        delete[] cnt_diffs_row_;
        delete[] p_row_;
        delete[] d_row_;

    }

    BoundedScoreLevenshteinDistance::BoundedScoreLevenshteinDistance(const BoundedScoreLevenshteinDistance& other) : // copy constructor
            scoring_(other.scoring_), threshold_(other.threshold_), row_length_(other.row_length_) {

        d_row_ = new lenSeqs_t[row_length_];
        p_row_ = new lenSeqs_t[row_length_];
        cnt_diffs_row_ = new lenSeqs_t[row_length_];
        cnt_diffs_p_row_ = new lenSeqs_t[row_length_];
        for (lenSeqs_t i = 0; i < row_length_; i++) {

            d_row_[i] = other.d_row_[i];
            p_row_[i] = other.p_row_[i];
            cnt_diffs_row_[i] = other.cnt_diffs_row_[i];
            cnt_diffs_p_row_[i] = other.cnt_diffs_p_row_[i];

        }

    }

    BoundedScoreLevenshteinDistance* BoundedScoreLevenshteinDistance::clone() const {
        return new BoundedScoreLevenshteinDistance(*this);
    }

    dist_t BoundedScoreLevenshteinDistance::distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) {
        return compute_gotoh_length_aware_early_row(ac.seq(i), ac.len(i), ac.seq(j), ac.len(j));
    }

    lenSeqs_t BoundedScoreLevenshteinDistance::compute_gotoh_length_aware_early_row(const char* s, const lenSeqs_t len_s,
            const char* t, const lenSeqs_t len_t) {

        // long computation not necessary if lengths differ too much
        if (((len_s > len_t) ? (len_s - len_t) : (len_t - len_s)) > threshold_) {
            return threshold_ + 1;
        }

        if (threshold_ == 0) {
            return lenSeqs_t(s != t);
        }

        const char* shorter = (len_s < len_t) ? s : t;
        lenSeqs_t len_shorter = std::min(len_s, len_t);
        const char* longer = (len_s >= len_t) ? s : t;
        lenSeqs_t len_longer = std::max(len_s, len_t);
        lenSeqs_t delta = len_longer - len_shorter;

        // (mis)match is only possibility when we have to consider only one diagonal
        // [happens only if (a) bound = delta = 0, or (b) bound = 1 and delta = 0, but (a) is already covered above]
        if ((threshold_ - delta) / 2 == 0 && (threshold_ + delta) / 2 == 0) {

            lenSeqs_t diffs = 0;
            for (auto i = 0; diffs <= threshold_ && i < len_shorter; i++) {
                diffs += (shorter[i] != longer[i]);
            }

            return diffs;

        }

        // initialise necessary section of first row
        d_row_[0] = scoring_.penalty_open;
        for (lenSeqs_t j = 1; j <= (threshold_ + delta) / 2 && j <= len_longer; j++) {

            d_row_[j] = d_row_[j - 1] + scoring_.penalty_extend;
            p_row_[j] = POS_INF;
            cnt_diffs_row_[j] = j;

        }
        d_row_[0] = 0;


        // similar to BoundedLevenshteinDistance, instead of storing the three full (score) matrices of Gotoh's scheme,
        // only a single row of D and P (Rd resp. Rp) and a few helper variables are used (assume we want to compute the entry D[i][j]):
        // - match: stores the score D[i - 1, j - 1] from the preceding row, its position in Rd was overwritten in the previous step
        // - min_val: temporarily contains newly computed score for D[i][j] in order to move D[i - 1][j] (stored in R[j]) to match
        //   before reassigning the position
        // - val_q: single value representing the entry to the left of the current row of matrix Q
        // - from_d: helper variable temporarily storing score obtained using a "alignment path" coming from matrix D
        // - from_pq: helper variable temporarily storing score obtained using a "alignment path" coming from matrix P (resp. Q)
        // since we are interested in the number of edit operations leading to the alignment score,
        // there are corresponding rows and helper variables counting the operations (differences)
        lenSeqs_t match, min_val, val_q, from_d, from_pq;
        lenSeqs_t j, diff, min_val_diff, diffs_q = 0;
        bool early;

        // compute remaining rows
        for (lenSeqs_t i = 1; i <= len_shorter; i++) {

            // handle left end
            d_row_[0] = scoring_.penalty_open + i * scoring_.penalty_extend;
            match = (i > (threshold_ - delta) / 2 + 1) ? d_row_[i - (threshold_ - delta) / 2 - 1] :
                    (((1 < i) && (i <= (threshold_ - delta) / 2 + 1)) * (d_row_[0] - scoring_.penalty_extend));
            val_q = POS_INF;
            diff = (i > (threshold_ - delta) / 2 + 1) ? cnt_diffs_row_[i - (threshold_ - delta) / 2 - 1] : (i - 1);
            cnt_diffs_row_[0] = i;

            early = true; // early termination flag

            // fill remaining row
            j = 1 + (i > (threshold_ - delta) / 2) * (i - (threshold_ - delta) / 2 - 1); // same as j = max(1, i - (threshold_ - delta) / 2) with signed integers
            d_row_[j - 1] = POS_INF;
            if (i + (threshold_ + delta) / 2 <= len_longer) {
                d_row_[i + (threshold_ + delta) / 2] = p_row_[i + (threshold_ + delta) / 2] = POS_INF;
            }

            for (; j <= i + (threshold_ + delta) / 2 && j <= len_longer; j++) {

                // arrays P & cntDiffsP
                from_d = d_row_[j] + scoring_.penalty_open + scoring_.penalty_extend;
                from_pq = p_row_[j] + scoring_.penalty_extend;

                if (from_d <= from_pq) {

                    p_row_[j] = from_d;
                    cnt_diffs_p_row_[j] = cnt_diffs_row_[j] + 1;

                } else {

                    p_row_[j] = from_pq;
                    cnt_diffs_p_row_[j]++;

                }

                // arrays Q & cntDiffsQ
                from_d = d_row_[j - 1] + scoring_.penalty_open + scoring_.penalty_extend;
                from_pq = val_q + scoring_.penalty_extend;
                if (from_d <= from_pq) {

                    val_q = from_d;
                    diffs_q = cnt_diffs_row_[j - 1] + 1;

                } else {

                    val_q = from_pq;
                    diffs_q++;

                }

                // arrays D & cntDiffs
                min_val = (match + (shorter[i - 1] != longer[j - 1]) * scoring_.penalty_mismatch);
                min_val_diff = diff + (shorter[i - 1] != longer[j - 1]);
                if (p_row_[j] < min_val) {

                    min_val = p_row_[j];
                    min_val_diff = cnt_diffs_p_row_[j];

                }
                if (val_q <= min_val){

                    min_val = val_q;
                    min_val_diff = diffs_q;

                }

                match = d_row_[j];
                d_row_[j] = min_val;

                diff = cnt_diffs_row_[j];
                cnt_diffs_row_[j] = min_val_diff;

                early &= ((cnt_diffs_row_[j] + llabs((long long)delta + (long long)i - (long long)j)) > threshold_); // improved e.t.

            }

            if (early) {// computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
                return threshold_ + 1;
            }

        }

        return (cnt_diffs_row_[len_longer] > threshold_) ? (threshold_ + 1) : cnt_diffs_row_[len_longer];

    }



    /* === QgramBoundedLevenshteinDistance === */

    QgramBoundedLevenshteinDistance::QgramBoundedLevenshteinDistance(Distance* m, lenSeqs_t threshold) {

        dist_fun_ = m;
        threshold_ = threshold;

    }

    QgramBoundedLevenshteinDistance::~QgramBoundedLevenshteinDistance() {
        delete dist_fun_;
    }

    QgramBoundedLevenshteinDistance::QgramBoundedLevenshteinDistance(const QgramBoundedLevenshteinDistance& other) : // copy constructor
            dist_fun_(other.dist_fun_->clone()), threshold_(other.threshold_) {

        // nothing else to do

    }

    QgramBoundedLevenshteinDistance* QgramBoundedLevenshteinDistance::clone() const {
        return new QgramBoundedLevenshteinDistance(*this);
    }

    dist_t QgramBoundedLevenshteinDistance::distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) {
        return (ac.qgram_diff(i, j) <= threshold_) ? dist_fun_->distance(ac, i, j) : (threshold_ + 1);
    }



    /* === BoundedScoreDistance === */

    BoundedScoreDistance::BoundedScoreDistance(lenSeqs_t max_seq_len, lenSeqs_t threshold, int match, int mismatch,
            int gap_open, int gap_extend) :
            scoring_(ScoringFunction(match, mismatch, gap_open, gap_extend)) {

        row_length_ = max_seq_len + 1;
        d_row_ = new lenSeqs_t[row_length_];
        p_row_ = new lenSeqs_t[row_length_];
        threshold_ = threshold;

    }

    BoundedScoreDistance::~BoundedScoreDistance() {

        delete[] p_row_;
        delete[] d_row_;

    }

    BoundedScoreDistance::BoundedScoreDistance(const BoundedScoreDistance& other) : // copy constructor
            scoring_(other.scoring_), threshold_(other.threshold_), row_length_(other.row_length_) {

        d_row_ = new lenSeqs_t[row_length_];
        p_row_ = new lenSeqs_t[row_length_];
        for (lenSeqs_t i = 0; i < row_length_; i++) {

            d_row_[i] = other.d_row_[i];
            p_row_[i] = other.p_row_[i];

        }

    }

    BoundedScoreDistance* BoundedScoreDistance::clone() const {
        return new BoundedScoreDistance(*this);
    }

    dist_t BoundedScoreDistance::distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) {
        return compute_bounded_gotoh_score(ac.seq(i), ac.len(i), ac.seq(j), ac.len(j));
    }

    lenSeqs_t BoundedScoreDistance::compute_bounded_gotoh_score(const char* s, const lenSeqs_t len_s,
            const char* t, const lenSeqs_t len_t) {

        const char* shorter = (len_s < len_t) ? s : t;
        lenSeqs_t len_shorter = std::min(len_s, len_t);
        const char* longer = (len_s >= len_t) ? s : t;
        lenSeqs_t len_longer = std::max(len_s, len_t);

        // initialise necessary section of first row
        d_row_[0] = scoring_.penalty_open;
        for (lenSeqs_t j = 1; j <= len_longer; j++) {

            d_row_[j] = d_row_[j - 1] + scoring_.penalty_extend;
            p_row_[j] = POS_INF;

        }
        d_row_[0] = 0;


        // similar to BoundedLevenshteinDistance, instead of storing the three full (score) matrices of Gotoh's scheme,
        // only a single row of D and P (Rd resp. Rp) and a few helper variables are used (assume we want to compute the entry D[i][j]):
        // - match: stores the score D[i - 1, j - 1] from the preceding row, its position in Rd was overwritten in the previous step
        // - tmp: temporarily contains newly computed score for D[i][j] in order to move D[i - 1][j] (stored in R[j]) to match
        //   before reassigning the position
        // - val_q: single value representing the entry to the left of the current row of matrix Q
        // the scheme is a bit simpler than the one of BoundedScoreLevenshteinDistance,
        // because we are only interested in the score (and not in the number of edit operations)
        lenSeqs_t match, val_q = POS_INF;
        bool early;

        // compute remaining rows
        for (lenSeqs_t i = 1; i <= len_shorter; i++) {

            // handle left end
            match = d_row_[0];
            d_row_[0] = scoring_.penalty_open + i * scoring_.penalty_extend;
            val_q = POS_INF;

            early = true; // early termination flag

            // fill remaining row
            for (lenSeqs_t j = 1; j <= len_longer; j++) {

                // array P
                p_row_[j] = std::min(d_row_[j] + scoring_.penalty_open + scoring_.penalty_extend, p_row_[j] + scoring_.penalty_extend);

                // array Q
                val_q = std::min(d_row_[j - 1] + scoring_.penalty_open + scoring_.penalty_extend, val_q + scoring_.penalty_extend);

                // array D
                lenSeqs_t tmp = (match + (shorter[i - 1] != longer[j - 1]) * scoring_.penalty_mismatch);
                match = d_row_[j];
                d_row_[j] = std::min({p_row_[j], val_q, tmp});

                early &= (d_row_[j] > threshold_);

            }

            if (early) {// computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
                return threshold_ + 1;
            }

        }

        return std::min({d_row_[len_longer], p_row_[len_longer], val_q});

    }


    /* === BoundedBandedScoreDistance == */

    BoundedBandedScoreDistance::BoundedBandedScoreDistance(lenSeqs_t max_seq_len, lenSeqs_t threshold, lenSeqs_t b, int match,
            int mismatch, int gap_open, int gap_extend) :
            scoring_(ScoringFunction(match, mismatch, gap_open, gap_extend)) {

        row_length_ = max_seq_len + 1;
        d_row_ = new lenSeqs_t[row_length_];
        p_row_ = new lenSeqs_t[row_length_];
        threshold_ = threshold;
        bands_per_side_ = b;

    }

    BoundedBandedScoreDistance::~BoundedBandedScoreDistance() {

        delete[] p_row_;
        delete[] d_row_;

    }

    BoundedBandedScoreDistance::BoundedBandedScoreDistance(const BoundedBandedScoreDistance& other) : // copy constructor
            scoring_(other.scoring_), threshold_(other.threshold_), bands_per_side_(other.bands_per_side_),
            row_length_(other.row_length_) {

        d_row_ = new lenSeqs_t[row_length_];
        p_row_ = new lenSeqs_t[row_length_];
        for (lenSeqs_t i = 0; i < row_length_; i++) {

            d_row_[i] = other.d_row_[i];
            p_row_[i] = other.p_row_[i];

        }

    }

    BoundedBandedScoreDistance* BoundedBandedScoreDistance::clone() const {
        return new BoundedBandedScoreDistance(*this);
    }

    dist_t BoundedBandedScoreDistance::distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) {
        return compute_bounded_banded_gotoh_score(ac.seq(i), ac.len(i), ac.seq(j), ac.len(j));
    }

    lenSeqs_t BoundedBandedScoreDistance::compute_bounded_banded_gotoh_score(const char* s, const lenSeqs_t len_s,
            const char* t, const lenSeqs_t len_t) {

        const char* shorter = (len_s < len_t) ? s : t;
        lenSeqs_t len_shorter = std::min(len_s, len_t);
        const char* longer = (len_s >= len_t) ? s : t;
        lenSeqs_t len_longer = std::max(len_s, len_t);

        // initialise necessary section of first row
        d_row_[0] = scoring_.penalty_open;
        for (lenSeqs_t j = 1; j <= bands_per_side_ && j <= len_longer; j++) {

            d_row_[j] = d_row_[j - 1] + scoring_.penalty_extend;
            p_row_[j] = POS_INF;

        }
        d_row_[0] = 0;


        // similar to BoundedLevenshteinDistance, instead of storing the three full (score) matrices of Gotoh's scheme,
        // only a single row of D and P (Rd resp. Rp) and a few helper variables are used (assume we want to compute the entry D[i][j]):
        // - match: stores the score D[i - 1, j - 1] from the preceding row, its position in Rd was overwritten in the previous step
        // - tmp: temporarily contains newly computed score for D[i][j] in order to move D[i - 1][j] (stored in R[j]) to match
        //   before reassigning the position
        // - val_q: single value representing the entry to the left of the current row of matrix Q
        // the scheme is a bit simpler than the one of BoundedScoreLevenshteinDistance,
        // because we are only interested in the score (and not in the number of edit operations)
        lenSeqs_t match, val_q = POS_INF;
        lenSeqs_t j;
        bool early;

        // compute remaining rows
        for (lenSeqs_t i = 1; i <= len_shorter; i++) {

            // handle left end
            d_row_[0] = scoring_.penalty_open + i * scoring_.penalty_extend;
            match = (i > bands_per_side_ + 1) ? d_row_[i - bands_per_side_ - 1] :
                    (((1 < i) && (i <= bands_per_side_ + 1)) * (d_row_[0] - scoring_.penalty_extend));
            val_q = POS_INF;

            early = true; // early termination flag

            // fill remaining row
            j = 1 + (i > bands_per_side_) * (i - bands_per_side_ - 1); // same as j = max(1, i - bands_per_side_) with signed integers
            d_row_[j - 1] = POS_INF;
            if (i + bands_per_side_ <= len_longer) {
                d_row_[i + bands_per_side_] = p_row_[i + bands_per_side_] = POS_INF;
            }

            for (; j <= i + bands_per_side_ && j <= len_longer; j++) {

                // array P
                p_row_[j] = std::min(d_row_[j] + scoring_.penalty_open + scoring_.penalty_extend, p_row_[j] + scoring_.penalty_extend);

                // array Q
                val_q = std::min(d_row_[j - 1] + scoring_.penalty_open + scoring_.penalty_extend, val_q + scoring_.penalty_extend);

                // array D
                lenSeqs_t tmp = (match + (shorter[i - 1] != longer[j - 1]) * scoring_.penalty_mismatch);
                match = d_row_[j];
                d_row_[j] = std::min({p_row_[j], val_q, tmp});

                early &= (d_row_[j] > threshold_);

            }

            if (early) {// computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
                return threshold_ + 1;
            }

        }

        return std::min({d_row_[len_longer], p_row_[len_longer], val_q});

    }


    /* === QgramBoundedScoreDistance === */

    QgramBoundedScoreDistance::QgramBoundedScoreDistance(Distance* m, lenSeqs_t threshold, lenSeqs_t max_operations) {

        dist_fun_ = m;
        threshold_ = threshold;
        max_operations_ = max_operations;

    }

    QgramBoundedScoreDistance::QgramBoundedScoreDistance(Distance* m, lenSeqs_t threshold, int match, int mismatch,
                                                         int gap_open, int gap_extend) {

        dist_fun_ = m;
        threshold_ = threshold;

        ScoringFunction sf(match, mismatch, gap_open, gap_extend);
        max_operations_ = threshold_ / std::min({sf.penalty_mismatch, sf.penalty_open, sf.penalty_extend});

    }

    QgramBoundedScoreDistance::~QgramBoundedScoreDistance() {
        delete dist_fun_;
    }

    QgramBoundedScoreDistance::QgramBoundedScoreDistance(const QgramBoundedScoreDistance& other) : // copy constructor
            dist_fun_(other.dist_fun_->clone()), threshold_(other.threshold_), max_operations_(other.max_operations_) {

        // nothing else to do

    }

    QgramBoundedScoreDistance* QgramBoundedScoreDistance::clone() const {
        return new QgramBoundedScoreDistance(*this);
    }

    dist_t QgramBoundedScoreDistance::distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) {
        return (ac.qgram_diff(i, j) <= max_operations_) ? dist_fun_->distance(ac, i, j) : (threshold_ + 1);
    }

}