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

#ifndef GEFAST_DISTANCES_HPP
#define GEFAST_DISTANCES_HPP

#include <limits>

#include "AmpliconCollections.hpp"
#include "Base.hpp"

namespace GeFaST {

    const lenSeqs_t POS_INF = INT16_MAX;
    const float POS_INF_FLOAT = std::numeric_limits<float>::max();


    /*
     * Representation of scoring function with affine gap costs.
     * Converted towards minimisation with a fixed match reward of 0.
     * The conversion is based on the ideas presented by Smith et al.
     * (https://www.ncbi.nlm.nih.gov/pubmed/7334527) in section "The Measures" (p. 39).
     * Similar to Swarm, the penalties are scaled down to smaller integer penalties (if possible).
     */
    struct ScoringFunction {

        ScoringFunction(int match, int mismatch, int gap_open, int gap_extend);

        lenSeqs_t penalty_mismatch;
        lenSeqs_t penalty_open;
        lenSeqs_t penalty_extend;

    };

    /*
     * Representation of an alignment in CIGAR format.
     *
     * Score information is currently not needed but could be easily included.
     */
    struct AlignmentInformation {

        std::string cigar; // alignment representation in CIGAR format
        lenSeqs_t length; // length of alignment (number of columns)
        lenSeqs_t num_diffs; // number of edit operations involved
//    lenSeqs_t score; // score of alignment

        AlignmentInformation(const std::string& c, const lenSeqs_t l, const lenSeqs_t n/*, const lenSeqs_t s*/);

    };

    /*
     * Compute alignment and perform backtracking to determine the CIGAR representation.
     */
    AlignmentInformation compute_gotoh_cigar_row(const char* s, const lenSeqs_t len_s, const char* t,
            const lenSeqs_t len_t, const ScoringFunction& scoring, lenSeqs_t* d_row, lenSeqs_t* p_row, char* backtrack);



    /*
     * Computes the bounded Levenshtein distance.
     * Uses only linear space as the alignment itself is not of interest.
     * Returns threshold + 1 (and, if possible, terminates early) when the distance is larger than the threshold.
     */
    class BoundedLevenshteinDistance : public Distance {

    public:
        BoundedLevenshteinDistance(lenSeqs_t max_seq_len, lenSeqs_t threshold);

        virtual ~BoundedLevenshteinDistance();

        BoundedLevenshteinDistance* clone() const override; // deep-copy clone method

        dist_t distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override;

    protected:
        BoundedLevenshteinDistance(const BoundedLevenshteinDistance& other); // copy constructor

        BoundedLevenshteinDistance(BoundedLevenshteinDistance&& other) = delete; // move constructor

        BoundedLevenshteinDistance& operator=(const BoundedLevenshteinDistance& other) = delete; // copy assignment operator

        BoundedLevenshteinDistance& operator=(BoundedLevenshteinDistance&& other) = delete; // move assignment operator

        /*
         * Compute the bounded Levenshtein distance.
         * The number of diagonals depends on the threshold and the length difference of s and t.
         * If the Levenshtein distance between s and t is larger than the threshold,
         * the returned distance is threshold + 1.
         */
        lenSeqs_t compute_length_aware_row(const char* s, const lenSeqs_t len_s, const char* t, const lenSeqs_t len_t);

        lenSeqs_t* matrix_row_; // single row of the DP-matrix
        lenSeqs_t threshold_; // distance threshold
        lenSeqs_t row_length_; // length of matrix row

    };


    /*
     * Computes a score-based bounded Levenshtein distance.
     * The distance between two amplicons is the number of edit operations
     * in their optimal alignment (according to the given scoring function).
     * Uses only linear space as the alignment itself is not of interest.
     * Returns threshold + 1 (and, if possible, terminates early) when the distance is larger than the threshold.
     */
    class BoundedScoreLevenshteinDistance : public Distance {

    public:
        BoundedScoreLevenshteinDistance(lenSeqs_t max_seq_len, lenSeqs_t threshold, int match, int mismatch, int gap_open, int gap_extend);

        ~BoundedScoreLevenshteinDistance();

        BoundedScoreLevenshteinDistance* clone() const override; // deep-copy clone method

        dist_t distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override;

    protected:
        BoundedScoreLevenshteinDistance(const BoundedScoreLevenshteinDistance& other); // copy constructor

        BoundedScoreLevenshteinDistance(BoundedScoreLevenshteinDistance&& other) = delete; // move constructor

        BoundedScoreLevenshteinDistance& operator=(const BoundedScoreLevenshteinDistance& other) = delete; // copy assignment operator

        BoundedScoreLevenshteinDistance& operator=(BoundedScoreLevenshteinDistance&& other) = delete; // move assignment operator

        /*
         * Compute a score-based bounded Levenshtein distance.
         * The number of diagonals depends on the threshold and the length difference of s and t.
         * If the distance between s and t is larger than the threshold,
         * the returned distance is threshold + 1.
         *
         * The method is a one-row adaption of Gotoh's algorithm using ideas similar to the ones
         * for the one-row version of the Needleman-Wunsch algorithm (score only).
         */
        lenSeqs_t compute_gotoh_length_aware_early_row(const char* s, const lenSeqs_t len_s, const char* t, const lenSeqs_t len_t);

        lenSeqs_t* d_row_; // one row of matrix D (cost of alignment of prefixes)
        lenSeqs_t* p_row_; // one row of matrix P (cost of alignment of prefixes that ends with gap in t)
        // row of matrix Q (for alignment of prefixes that ends with gap in s) not necessary (handled as one value of a row)
        lenSeqs_t* cnt_diffs_row_; // number of differences corresponding to alignments based on D
        lenSeqs_t* cnt_diffs_p_row_; // number of differences corresponding to alignments based on P
        ScoringFunction scoring_; // used scoring function
        lenSeqs_t threshold_; // distance threshold
        lenSeqs_t row_length_; // length of matrix row

    };


    /*
     * Computes the provided version of the bounded Levenshtein distance.
     * Tries to avoid the more expensive alignment computation by first computing
     * (a lower bound on) the q-gram distance (which is also a lower bound on the Levenshtein distance).
     */
    class QgramBoundedLevenshteinDistance : public Distance {

    public:
        QgramBoundedLevenshteinDistance(Distance* m, lenSeqs_t threshold);

        virtual ~QgramBoundedLevenshteinDistance();

        QgramBoundedLevenshteinDistance* clone() const override; // deep-copy clone method

        dist_t distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override;

    protected:
        QgramBoundedLevenshteinDistance(const QgramBoundedLevenshteinDistance& other); // copy constructor

        QgramBoundedLevenshteinDistance(QgramBoundedLevenshteinDistance&& other) = delete; // move constructor

        QgramBoundedLevenshteinDistance& operator=(const QgramBoundedLevenshteinDistance& other) = delete; // copy assignment operator

        QgramBoundedLevenshteinDistance& operator=(QgramBoundedLevenshteinDistance&& other) = delete; // move assignment operator

    private:
        Distance* dist_fun_; // distance function used when the amplicon pair is not rejected by the q-gram test
        lenSeqs_t threshold_; // distance threshold

    };


    /*
     * Computes the (bounded) score of an optimal alignment between two amplicons as their distance.
     * Uses only linear space as the alignment itself is not of interest.
     * Returns threshold + 1 (and, if possible, terminates early) when the score is larger than the threshold.
     */
    class BoundedScoreDistance : public Distance {

    public:
        BoundedScoreDistance(lenSeqs_t max_seq_len, lenSeqs_t threshold, int match, int mismatch, int gap_open, int gap_extend);

        virtual ~BoundedScoreDistance();

        BoundedScoreDistance* clone() const override; // deep-copy clone method

        dist_t distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override;

    protected:
        BoundedScoreDistance(const BoundedScoreDistance& other); // copy constructor

        BoundedScoreDistance(BoundedScoreDistance&& other) = delete; // move constructor

        BoundedScoreDistance& operator=(const BoundedScoreDistance& other) = delete; // copy assignment operator

        BoundedScoreDistance& operator=(BoundedScoreDistance&& other) = delete; // move assignment operator

        /*
         * Compute the (bounded) score of an optimal alignment.
         * If the score of an optimal alignment between s and t is larger than the threshold,
         * the returned score is threshold + 1.
         */
        lenSeqs_t compute_bounded_gotoh_score(const char* s, const lenSeqs_t len_s, const char* t, const lenSeqs_t len_t);

        lenSeqs_t* d_row_; // one row of matrix D (cost of alignment of prefixes)
        lenSeqs_t* p_row_; // one row of matrix P (cost of alignment of prefixes that ends with gap in t)
        ScoringFunction scoring_; // used scoring function
        lenSeqs_t threshold_; // score threshold
        lenSeqs_t row_length_; // length of matrix row

    };

    /*
     * Computes the (bounded & banded) score of an optimal alignment between two amplicons as their distance.
     * Speeds up the computation by considering only a limited number of diagonals besides the main diagonal
     * (even though this could prevent the computation from finding the overall optimal alignment).
     * Uses only linear space as the alignment itself is not of interest.
     * Returns threshold + 1 (and, if possible, terminates early) when the score is larger than the threshold.
     */
    class BoundedBandedScoreDistance : public Distance {

    public:
        BoundedBandedScoreDistance(lenSeqs_t max_seq_len, lenSeqs_t threshold, lenSeqs_t b, int match, int mismatch,
                int gap_open, int gap_extend);

        virtual ~BoundedBandedScoreDistance();

        BoundedBandedScoreDistance* clone() const override; // deep-copy clone method

        dist_t distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override;

    protected:
        BoundedBandedScoreDistance(const BoundedBandedScoreDistance& other); // copy constructor

        BoundedBandedScoreDistance(BoundedBandedScoreDistance&& other) = delete; // move constructor

        BoundedBandedScoreDistance& operator=(const BoundedBandedScoreDistance& other) = delete; // copy assignment operator

        BoundedBandedScoreDistance& operator=(BoundedBandedScoreDistance&& other) = delete; // move assignment operator

        /*
         * Compute the (bounded & banded) score of an optimal alignment.
         * If the score of an optimal alignment between s and t is larger than the threshold,
         * the returned score is threshold + 1.
         * Restricts the alignment to 2 * bands_per_side + 1 diagonals of the DP-matrix.
         */
        lenSeqs_t compute_bounded_banded_gotoh_score(const char* s, const lenSeqs_t len_s, const char* t, const lenSeqs_t len_t);

        lenSeqs_t* d_row_; // one row of matrix D (cost of alignment of prefixes)
        lenSeqs_t* p_row_; // one row of matrix P (cost of alignment of prefixes that ends with gap in t)
        ScoringFunction scoring_; // used scoring function
        lenSeqs_t threshold_; // score threshold
        lenSeqs_t bands_per_side_; // number of diagonals (bands) per side of the main diagonal considered
        lenSeqs_t row_length_; // length of matrix row

    };


    /*
     * Computes the given version of the (bounded) score of an optimal alignment between two amplicons.
     * Tries to avoid the more expensive alignment computation by first computing a lower bound
     * on the q-gram distance (which is also a lower bound on the Levenshtein distance)
     * and comparing this lower bound times the minimum edit-operation cost with the threshold.
     */
    class QgramBoundedScoreDistance : public Distance {

    public:
        QgramBoundedScoreDistance(Distance* m, lenSeqs_t threshold, lenSeqs_t max_operations);

        QgramBoundedScoreDistance(Distance* m, lenSeqs_t threshold, int match, int mismatch, int gap_open, int gap_extend);

        virtual ~QgramBoundedScoreDistance();

        QgramBoundedScoreDistance* clone() const override; // deep-copy clone method

        dist_t distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override;

    protected:
        QgramBoundedScoreDistance(const QgramBoundedScoreDistance& other); // copy constructor

        QgramBoundedScoreDistance(QgramBoundedScoreDistance&& other) = delete; // move constructor

        QgramBoundedScoreDistance& operator=(const QgramBoundedScoreDistance& other) = delete; // copy assignment operator

        QgramBoundedScoreDistance& operator=(QgramBoundedScoreDistance&& other) = delete; // move assignment operator

    private:
        Distance* dist_fun_; // distance function used when the amplicon pair is not rejected by the q-gram test
        lenSeqs_t threshold_; // score threshold
        lenSeqs_t max_operations_; // operations / q-gram threshold

    };



    /*
     * Computes the Manhattan distance between the features of two amplicons.
     *
     * Can only be used with FeatureAmpliconCollection.
     */
    class ManhattanDistance : public FeatureDistance {

    public:
        ManhattanDistance* clone() const override; // deep-copy clone method

        dist_t distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override;

        /*
         * Compute the Manhattan distance between the features.
         */
        dist_t distance(const FeatureAmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override;

    };


    /*
     * Computes the Euclidean distance between the features of two amplicons.
     *
     * Can only be used with FeatureAmpliconCollection.
     */
    class EuclideanDistance : public FeatureDistance {

    public:
        EuclideanDistance* clone() const override; // deep-copy clone method

        dist_t distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override;

        /*
         * Compute the Euclidean distance between the features.
         */
        dist_t distance(const FeatureAmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override;

    };


    /*
     * Computes the cosine distance between the features of two amplicons as 1 - (cosine similarity).
     *
     * Can only be used with FeatureAmpliconCollection.
     */
    class CosineDistance : public FeatureDistance {

    public:
        CosineDistance* clone() const override; // deep-copy clone method

        dist_t distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override;

        /*
         * Compute the cosine distance between the features.
         */
        dist_t distance(const FeatureAmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override;

    };

    /*
     * Computes the Pearson distance between the features of two amplicons as 1 - (Pearson correlation).
     *
     * Can only be used with FeatureAmpliconCollection.
     */
    class PearsonDistance : public FeatureDistance {

    public:
        PearsonDistance* clone() const override; // deep-copy clone method

        dist_t distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override;

        /*
         * Compute the Pearson distance between the features.
         */
        dist_t distance(const FeatureAmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override;

    };

}

#endif //GEFAST_DISTANCES_HPP
