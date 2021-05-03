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

#ifndef GEFAST_QUALITYWEIGHTEDDISTANCES_HPP
#define GEFAST_QUALITYWEIGHTEDDISTANCES_HPP

#include "Base.hpp"
#include "Distances.hpp"

namespace GeFaST {

    /* ===== Boosting options ===== */

    /*
     * Boosting function corresponding to the identity function.
     */
    struct NoEffectBooster {

        NoEffectBooster() = default;

        explicit NoEffectBooster(const Configuration& config) {
            // nothing to do
        }

        inline float boost(float in) const {
            return in;
        }

    };

    /*
     * Multiplicative boosting function, multiplying the probability by a constant factor and capping the value at 1.
     */
    struct MultBooster {

        MultBooster() = default;

        explicit MultBooster(const Configuration& config) {
            if (config.misc.find("mult_factor") != config.misc.end()) c_ = std::stof(config.misc.find("mult_factor")->second);
        }

        explicit MultBooster(float c) : c_(c) {
            // nothing else to do
        }

        inline float boost(float in) const {
            return std::min(1.0f, c_ * in);
        }

    private:
        float c_ = 1.0; // constant factor

    };

    /*
     * Root boosting function, extracting the root of the requested degree d
     * and shifting the result by a specified amount s towards 0.
     * Using the "full shift" option, the roots are shifted in a way
     * such that the boosted error probability for the highest quality score in the used encoding is 0.
     */
    struct RootBooster {

        RootBooster() = default;

        explicit RootBooster(const Configuration& config) {

            if (config.misc.find("root_degree") != config.misc.end()) d_ = std::stof(config.misc.find("root_degree")->second);
            if (config.misc.find("root_shift") != config.misc.end()) {

                auto val = config.misc.find("root_shift")->second;
                if (val == "full") {

                    auto qe = config.build_quality_encoding();
                    lenSeqs_t num_scores = qe->get_accepted_scores().size();
                    delete qe;

                    s_ = powf(powf(10.0, -static_cast<float>(num_scores - 1) / 10.0f), 1.0f / d_);

                } else {
                    s_ = std::stof(val);
                }

            }

        }

        explicit RootBooster(float d) : RootBooster(d, 0.0) {
            // nothing else to do
        }

        RootBooster(float d, float s) : d_(d), s_(s) {
            // nothing else to do
        }

        inline float boost(float in) const {
            return std::max(0.0f, powf(in, 1.0f / d_) - s_);
        }

    private:
        float d_ = 2.0; // degree of root
        float s_ = 0.0; // shift towards zero

    };

    /*
     * Boosting function which switches from exponential to linear decay after a specified quality level.
     */
    struct LinearBooster {

        LinearBooster() = default;

        explicit LinearBooster(const Configuration& config) {

            if (config.misc.find("linear_start") != config.misc.end()) start_ = std::stof(config.misc.find("linear_start")->second);
            if (config.misc.find("linear_levels") != config.misc.end()) {
                levels_ = std::stof(config.misc.find("linear_levels")->second);
            } else {

                auto qe = config.build_quality_encoding();
                levels_ = qe->get_accepted_scores().size();
                delete qe;

            }

            start_err_prob_ = powf(10.0, -start_ / 10.0f);
            steps_ = levels_ - 1 - start_;

        }

        LinearBooster(float s, float l) : start_(s), levels_(l) {

            start_err_prob_ = powf(10.0, -start_ / 10.0f);
            steps_ = levels_ - 1 - start_;

        }

        inline float boost(float in) const {

            float score = roundf(-10 * log10f(in));
            return (score <= start_) ? in : (start_err_prob_ * (steps_ - (score - start_)) / steps_);

        }

    private:
        float start_ = 0.0; // quality level at which the linear part starts
        float start_err_prob_ = 1.0; // error probability corresponding to quality level start_
        float steps_ = 41.0; // number of levels in the linear part
        float levels_ = 42.0; // number of levels in the used quality encoding (default value corresponds to Illumina 1.8+)

    };



    /* ===== Quality-weighting mechanisms =====
     *
     * Each *Scores class provides methods for computing (resp. looking up)
     * the (weighted) scores of matches, mismatches, gap-openings and gap-extensions.
     * These are used in the Distance implementations subsequently provided.
     */

    /*
     * Provides the unweighted scores of the operations, completely ignoring the quality information.
     */
    class UnweightedScores {

    public:
        UnweightedScores(const QualityEncoding<>& qe, const Configuration& config) {

            ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
            mismatch_ = sf.penalty_mismatch;
            gap_open_ = sf.penalty_open;
            gap_extend_ = sf.penalty_extend;

        }

        inline float match(char qual_a, char qual_b) { return 0; }
        inline float mismatch(char qual_a, char qual_b) { return mismatch_; }
        inline float substitute(char a, char qual_a, char b, char qual_b) { return (a == b) ? 0 : mismatch_; }
        inline float gap_open(char qual_a) { return gap_open_; }
        inline float gap_extend(char qual_a) { return gap_extend_; }

        inline float lookup_match(char qual_a, char qual_b) { return 0; }
        inline float lookup_mismatch(char qual_a, char qual_b) { return mismatch_; }
        inline float lookup_substitute(char a, char qual_a, char b, char qual_b) { return (a == b) ? 0 : mismatch_; }
        inline float lookup_gap_open(char qual_a) { return gap_open_; }
        inline float lookup_gap_extend(char qual_a) { return gap_extend_; }

    private:
        float mismatch_; // mismatch penalty
        float gap_open_; // gap-opening penalty
        float gap_extend_; // gap-extension penalty

    };


    /* === Clement et al., The GNUMAP algorithm: unbiased probabilistic mapping of oligonucleotides from next-generation sequencing === */

    /*
     * Provides quality-weighted scores using the mechanism from:
     * "The GNUMAP algorithm: unbiased probabilistic mapping of oligonucleotides from next-generation sequencing"
     * (Clement et al., 2010, https://doi.org/10.1093/bioinformatics/btp614)
     *
     * Affects only substitutions (modifies the mismatch penalty).
     *
     * The costs of matches are not quality-weighted (see ClementScoresWithMatches for a version also weighting matches).
     */
    template<typename Bi, typename Bo>
    class ClementScores {

    public:
        ClementScores(const QualityEncoding<>& qe, const Configuration& config) : bInner_(config), bOuter_(config), qe_(qe) {

            ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
            mismatch_ = sf.penalty_mismatch;
            gap_open_ = sf.penalty_open;
            gap_extend_ = sf.penalty_extend;

            alphabet_size_ = (config.alphabet.empty()) ? 4 : config.alphabet.length();

            // precompute the weighted mismatch scores
            auto chars = qe_.get_accepted_scores();
            num_scores_ = chars.size();
            score_shift_ = qe_.get_offset() + qe_.get_from();
            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_mismatch_scores_[index] = mismatch(a, b);

                }
            }

        }

        virtual ~ClementScores() {
            delete[] weighted_mismatch_scores_;
        }

        ClementScores(const ClementScores& other) : // copy constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                mismatch_(other.mismatch_), gap_open_(other.gap_open_), gap_extend_(other.gap_extend_),
                score_shift_(other.score_shift_), num_scores_(other.num_scores_) {

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

        }

        ClementScores(ClementScores&& other) noexcept : // move constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                mismatch_(other.mismatch_), gap_open_(other.gap_open_), gap_extend_(other.gap_extend_),
                weighted_mismatch_scores_(other.weighted_mismatch_scores_), score_shift_(other.score_shift_),
                num_scores_(other.num_scores_) {

            other.weighted_mismatch_scores_ = nullptr;

        }

        ClementScores& operator=(const ClementScores& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

            return *this;

        }

        ClementScores& operator=(ClementScores&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = other.weighted_mismatch_scores_;
            other.weighted_mismatch_scores_ = nullptr;

            return *this;

        }


        virtual inline float match(char qual_a, char qual_b) {
            return 0;
        }

        inline float mismatch(char qual_a, char qual_b) {

            float prob_a = 1 - bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float prob_b = 1 - bInner_.boost(qe_.encoded_quality_to_probability(qual_b));

            return mismatch_ * bOuter_.boost(1 + (2 - prob_a - prob_b - alphabet_size_ * (1 - prob_a * prob_b)) / ((alphabet_size_ - 1) * (alphabet_size_ - 1)));

        }

        inline float substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? match(qual_a, qual_b) : mismatch(qual_a, qual_b);
        }

        inline float gap_open(char qual_a) { return gap_open_; }
        inline float gap_extend(char qual_a) { return gap_extend_; }


        virtual inline float lookup_match(char qual_a, char qual_b) {
            return 0;
        }

        inline float lookup_mismatch(char qual_a, char qual_b) {
            return weighted_mismatch_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

        inline float lookup_substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? lookup_match(qual_a, qual_b) : lookup_mismatch(qual_a, qual_b);
        }

        inline float lookup_gap_open(char qual_a) { return gap_open_; }
        inline float lookup_gap_extend(char qual_a) { return gap_extend_; }

    protected:
        Bi bInner_; // inner boosting function (for individual error probabilities)
        Bo bOuter_; // outer boosting function (for combined error probabilities)
        QualityEncoding<> qe_; // encoding of the quality score
        lenSeqs_t alphabet_size_; // size of the used sequence alphabet

        float mismatch_; // mismatch penalty
        float gap_open_; // gap-opening penalty
        float gap_extend_; // gap-extension penalty

        float* weighted_mismatch_scores_; // matrix of weighted mismatch scores as 1-dimensional C-style array
        char score_shift_; // (offset + from) of the quality encoding to normalise the quality scores
        unsigned long num_scores_; // number of levels in the quality encoding

    };

    template<typename Bi, typename Bo>
    class ClementScoresWithMatches : public ClementScores<Bi, Bo> {

        using ClementScores<Bi, Bo>::bInner_;
        using ClementScores<Bi, Bo>::bOuter_;
        using ClementScores<Bi, Bo>::qe_;
        using ClementScores<Bi, Bo>::alphabet_size_;
        using ClementScores<Bi, Bo>::mismatch_;
        using ClementScores<Bi, Bo>::score_shift_;
        using ClementScores<Bi, Bo>::num_scores_;

        using ClementScores<Bi, Bo>::mismatch;
        using ClementScores<Bi, Bo>::lookup_mismatch;

    public:
        ClementScoresWithMatches(const QualityEncoding<>& qe, const Configuration& config) : ClementScores<Bi, Bo>(qe, config) {

            // also precompute the weighted match scores
            auto chars = qe_.get_accepted_scores();
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_match_scores_[index] = match(a, b);

                }
            }

        }

        ~ClementScoresWithMatches() {
            delete[] weighted_match_scores_;
        }

        ClementScoresWithMatches(const ClementScoresWithMatches& other) : // copy constructor
                ClementScores<Bi, Bo>(other) {

            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

        }

        ClementScoresWithMatches(ClementScoresWithMatches&& other) noexcept : // move constructor
                ClementScores<Bi, Bo>(std::move(other)), weighted_match_scores_(other.weighted_match_scores_) {

            other.weighted_match_scores_ = nullptr;

        }

        ClementScoresWithMatches& operator=(const ClementScoresWithMatches& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            ClementScores<Bi, Bo>::operator=(other);

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

            return *this;

        }

        ClementScoresWithMatches& operator=(ClementScoresWithMatches&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            ClementScores<Bi, Bo>::operator=(std::move(other));

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = other.weighted_match_scores_;
            other.weighted_match_scores_ = nullptr;

            return *this;

        }


        inline float match(char qual_a, char qual_b) override {

            float prob_a = 1 - bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float prob_b = 1 - bInner_.boost(qe_.encoded_quality_to_probability(qual_b));

            return mismatch_ * bOuter_.boost((prob_a + prob_b - alphabet_size_ * prob_a * prob_b + alphabet_size_ - 2) / (alphabet_size_ - 1));

        }

        inline float lookup_match(char qual_a, char qual_b) override {
            return weighted_match_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

    protected:
        float* weighted_match_scores_; // matrix of weighted match scores as 1-dimensional C-style array

    };


    /* === Malde, The effect of sequence quality on sequence alignment === */

    // Version A: alternative substitution matrix

    /*
     * Provides quality-weighted scores using the mechanism from:
     * "The effect of sequence quality on sequence alignment"
     * (Malde, 2008, https://doi.org/10.1093/bioinformatics/btn052)
     *
     * Affects only substitutions (ignores scoring function, only considers the error rates).
     *
     * The costs of matches are not quality-weighted (see MaldeScoresWithMatchesA for a version also weighting matches).
     */
    template<typename Bi, typename Bo>
    class MaldeScoresA {

    public:
        MaldeScoresA(const QualityEncoding<>& qe, const Configuration& config) : bInner_(config), bOuter_(config), qe_(qe) {

            alphabet_size_ = (config.alphabet.empty()) ? 4 : config.alphabet.length();
            scaling_factor_ = 1.0f / std::log(2.0f);
            if (config.misc.find("scaling_factor") != config.misc.end()) scaling_factor_ = std::stof(config.misc.find("scaling_factor")->second);

            // compute contents of alternative scoring matrix only considering the error probabilities
            // and record the minimum and maximum (mis)match scores
            auto chars = qe_.get_accepted_scores();
            num_scores_ = chars.size();
            score_shift_ = qe_.get_offset() + qe_.get_from();
            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            min_match_ = std::numeric_limits<float>::infinity();
            max_match_ = -std::numeric_limits<float>::infinity();
            min_mismatch_ = std::numeric_limits<float>::infinity();
            max_mismatch_ = -std::numeric_limits<float>::infinity();
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);

                    auto tmp = match_helper(a, b);
                    min_match_ = std::min(min_match_, tmp);
                    max_match_ = std::max(max_match_, tmp);

                    weighted_mismatch_scores_[index] = mismatch_helper(a, b);
                    min_mismatch_ = std::min(min_mismatch_, weighted_mismatch_scores_[index]);
                    max_mismatch_ = std::max(max_mismatch_, weighted_mismatch_scores_[index]);

                }
            }

            // using the minimum and maximum scores, map the score range of the matrix
            // to the score range of the minimising scoring function (only for mismatch)
            ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
            mismatch_ = sf.penalty_mismatch;
            gap_open_ = sf.penalty_open;
            gap_extend_ = sf.penalty_extend;
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_mismatch_scores_[index] = (max_mismatch_ - weighted_mismatch_scores_[index]) / (max_mismatch_ - min_mismatch_) * mismatch_;

                }
            }

        }

        virtual ~MaldeScoresA() {
            delete[] weighted_mismatch_scores_;
        }

        MaldeScoresA(const MaldeScoresA& other) : // copy constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                scaling_factor_(other.scaling_factor_), mismatch_(other.mismatch_), gap_open_(other.gap_open_),
                gap_extend_(other.gap_extend_), min_match_(other.min_match_), max_match_(other.max_match_),
                min_mismatch_(other.min_mismatch_), max_mismatch_(other.max_mismatch_),
                score_shift_(other.score_shift_), num_scores_(other.num_scores_) {

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

        }

        MaldeScoresA(MaldeScoresA&& other) noexcept : // move constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                scaling_factor_(other.scaling_factor_), mismatch_(other.mismatch_), gap_open_(other.gap_open_),
                gap_extend_(other.gap_extend_), min_match_(other.min_match_), max_match_(other.max_match_),
                min_mismatch_(other.min_mismatch_), max_mismatch_(other.max_mismatch_),
                weighted_mismatch_scores_(other.weighted_mismatch_scores_), score_shift_(other.score_shift_),
                num_scores_(other.num_scores_) {

            other.weighted_mismatch_scores_ = nullptr;

        }

        MaldeScoresA& operator=(const MaldeScoresA& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            scaling_factor_ = other.scaling_factor_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            min_match_ = other.min_match_;
            max_match_ = other.max_match_;
            min_mismatch_ = other.min_mismatch_;
            max_mismatch_ = other.max_mismatch_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

            return *this;

        }

        MaldeScoresA& operator=(MaldeScoresA&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            scaling_factor_ = other.scaling_factor_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            min_match_ = other.min_match_;
            max_match_ = other.max_match_;
            min_mismatch_ = other.min_mismatch_;
            max_mismatch_ = other.max_mismatch_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = other.weighted_mismatch_scores_;
            other.weighted_mismatch_scores_ = nullptr;

            return *this;

        }


        virtual inline float match(char qual_a, char qual_b) {
            return 0;
        }

        inline float mismatch(char qual_a, char qual_b) {
            return (max_mismatch_ - mismatch_helper(qual_a, qual_b)) / (max_mismatch_ - min_mismatch_) * mismatch_;
        }

        inline float substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? match(qual_a, qual_b) : mismatch(qual_a, qual_b);
        }

        inline float gap_open(char qual_a) { return gap_open_; }
        inline float gap_extend(char qual_a) { return gap_extend_; }


        virtual inline float lookup_match(char qual_a, char qual_b) {
            return 0;
        }

        inline float lookup_mismatch(char qual_a, char qual_b) {
            return weighted_mismatch_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

        inline float lookup_substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? lookup_match(qual_a, qual_b) : lookup_mismatch(qual_a, qual_b);
        }

        inline float lookup_gap_open(char qual_a) { return gap_open_; }
        inline float lookup_gap_extend(char qual_a) { return gap_extend_; }

    protected:
        MaldeScoresA() = default; // only for constructors of derived classes

        // helper function encapsulating the computation of the alternative match score
        inline float match_helper(char qual_a, char qual_b) {

            float err_prob_a = bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float err_prob_b = bInner_.boost(qe_.encoded_quality_to_probability(qual_b));
            float comb_err = bOuter_.boost(err_prob_a + err_prob_b - (1.0f + 1.0f / (alphabet_size_ - 1)) * err_prob_a * err_prob_b);

            return scaling_factor_ * std::log((1.0f - comb_err) * alphabet_size_); //TODO? boosting after difference / multiplication with alphabet size better?

        }

        // helper function encapsulating the computation of the alternative mismatch score
        inline float mismatch_helper(char qual_a, char qual_b) {

            float err_prob_a = bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float err_prob_b = bInner_.boost(qe_.encoded_quality_to_probability(qual_b));
            float comb_err = bOuter_.boost(err_prob_a + err_prob_b - (1.0f + 1.0f / (alphabet_size_ - 1)) * err_prob_a * err_prob_b);

            return scaling_factor_ * std::log(comb_err / (alphabet_size_ - 1.0f) * alphabet_size_ ); //TODO? boosting after difference / multiplication with alphabet size better?

        }


        Bi bInner_; // inner boosting function (for individual error probabilities)
        Bo bOuter_; // outer boosting function (for combined error probabilities)
        QualityEncoding<> qe_; // encoding of the quality score
        lenSeqs_t alphabet_size_; // size of the used sequence alphabet
        float scaling_factor_;

        float mismatch_; // mismatch penalty
        float gap_open_; // gap-opening penalty
        float gap_extend_; // gap-extension penalty
        float min_match_; // minimum match score in the alternative scoring matrix
        float max_match_; // maximum match score in the alternative scoring matrix
        float min_mismatch_; // minimum mismatch score in the alternative scoring matrix
        float max_mismatch_; // maximum mismatch score in the alternative scoring matrix

        float* weighted_mismatch_scores_; // matrix of weighted mismatch scores as 1-dimensional C-style array
        char score_shift_; // (offset + from) of the quality encoding to normalise the quality scores
        unsigned long num_scores_; // number of levels in the quality encoding

    };

    template<typename Bi, typename Bo>
    class MaldeScoresWithMatchesA : public MaldeScoresA<Bi, Bo> {

        using MaldeScoresA<Bi, Bo>::bInner_;
        using MaldeScoresA<Bi, Bo>::bOuter_;
        using MaldeScoresA<Bi, Bo>::qe_;
        using MaldeScoresA<Bi, Bo>::alphabet_size_;
        using MaldeScoresA<Bi, Bo>::scaling_factor_;
        using MaldeScoresA<Bi, Bo>::mismatch_;
        using MaldeScoresA<Bi, Bo>::gap_open_;
        using MaldeScoresA<Bi, Bo>::gap_extend_;
        using MaldeScoresA<Bi, Bo>::min_match_;
        using MaldeScoresA<Bi, Bo>::max_match_;
        using MaldeScoresA<Bi, Bo>::min_mismatch_;
        using MaldeScoresA<Bi, Bo>::max_mismatch_;
        using MaldeScoresA<Bi, Bo>::weighted_mismatch_scores_;
        using MaldeScoresA<Bi, Bo>::score_shift_;
        using MaldeScoresA<Bi, Bo>::num_scores_;

        using MaldeScoresA<Bi, Bo>::match_helper;
        using MaldeScoresA<Bi, Bo>::mismatch;
        using MaldeScoresA<Bi, Bo>::mismatch_helper;
        using MaldeScoresA<Bi, Bo>::lookup_mismatch;

    public:
        MaldeScoresWithMatchesA(const QualityEncoding<>& qe, const Configuration& config) {

            bInner_ = Bi(config);
            bOuter_ = Bo(config);
            qe_ = qe;

            alphabet_size_ = (config.alphabet.empty()) ? 4 : config.alphabet.length();
            scaling_factor_ = 1.0f / std::log(2.0f);
            if (config.misc.find("scaling_factor") != config.misc.end()) scaling_factor_ = std::stof(config.misc.find("scaling_factor")->second);

            // compute contents of alternative scoring matrix only considering the error probabilities
            // and record the minimum and maximum (mis)match scores
            auto chars = qe_.get_accepted_scores();
            num_scores_ = chars.size();
            score_shift_ = qe_.get_offset() + qe_.get_from();
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            min_match_ = std::numeric_limits<float>::infinity();
            max_match_ = -std::numeric_limits<float>::infinity();
            min_mismatch_ = std::numeric_limits<float>::infinity();
            max_mismatch_ = -std::numeric_limits<float>::infinity();
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);

                    weighted_match_scores_[index] = match_helper(a, b);
                    weighted_mismatch_scores_[index] = mismatch_helper(a, b);

                    min_match_ = std::min(min_match_, weighted_match_scores_[index]);
                    max_match_ = std::max(max_match_, weighted_match_scores_[index]);
                    min_mismatch_ = std::min(min_mismatch_, weighted_mismatch_scores_[index]);
                    max_mismatch_ = std::max(max_mismatch_, weighted_mismatch_scores_[index]);

                }
            }

            // using the minimum and maximum scores, map the score range of the matrix
            // to the score range of the minimising scoring function (for match and mismatch)
            ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
            mismatch_ = sf.penalty_mismatch;
            gap_open_ = sf.penalty_open;
            gap_extend_ = sf.penalty_extend;
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);

                    weighted_match_scores_[index] = (max_match_ - weighted_match_scores_[index]) / (max_match_ - min_match_) * mismatch_;
                    weighted_mismatch_scores_[index] = (max_mismatch_ - weighted_mismatch_scores_[index]) / (max_mismatch_ - min_mismatch_) * mismatch_;

                }
            }

        }

        ~MaldeScoresWithMatchesA() {
            delete[] weighted_match_scores_;
        }

        MaldeScoresWithMatchesA(const MaldeScoresWithMatchesA& other) : // copy constructor
                MaldeScoresA<Bi, Bo>(other) {

            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

        }

        MaldeScoresWithMatchesA(MaldeScoresWithMatchesA&& other) noexcept : // move constructor
                MaldeScoresA<Bi, Bo>(std::move(other)), weighted_match_scores_(other.weighted_match_scores_) {

            other.weighted_match_scores_ = nullptr;

        }

        MaldeScoresWithMatchesA& operator=(const MaldeScoresWithMatchesA& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            MaldeScoresA<Bi, Bo>::operator=(other);

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

            return *this;

        }

        MaldeScoresWithMatchesA& operator=(MaldeScoresWithMatchesA&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            MaldeScoresA<Bi, Bo>::operator=(std::move(other));

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = other.weighted_match_scores_;
            other.weighted_match_scores_ = nullptr;

            return *this;

        }


        inline float match(char qual_a, char qual_b) override {
            return (max_match_ - match_helper(qual_a, qual_b)) / (max_match_ - min_match_) * mismatch_;
        }

        inline float lookup_match(char qual_a, char qual_b) override {
            return weighted_match_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

    protected:
        float* weighted_match_scores_; // matrix of weighted match scores as 1-dimensional C-style array

    };


    // Version B: use combined error simply as weight

    /*
     * Provides quality-weighted scores using the mechanism from:
     * "The effect of sequence quality on sequence alignment"
     * (Malde, 2008, https://doi.org/10.1093/bioinformatics/btn052)
     *
     * Affects only substitutions (weights the mismatch penalty with the error rates).
     *
     * The costs of matches are not quality-weighted (see MaldeScoresWithMatchesB for a version also weighting matches).
     */
    template<typename Bi, typename Bo>
    class MaldeScoresB {

    public:
        MaldeScoresB(const QualityEncoding<>& qe, const Configuration& config) : bInner_(config), bOuter_(config), qe_(qe) {

            ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
            mismatch_ = sf.penalty_mismatch;
            gap_open_ = sf.penalty_open;
            gap_extend_ = sf.penalty_extend;

            alphabet_size_ = (config.alphabet.empty()) ? 4 : config.alphabet.length();

            // fill matrix of weighted mismatch scores by weighting the mismatch penalty with the combined error probability
            auto chars = qe_.get_accepted_scores();
            num_scores_ = chars.size();
            score_shift_ = qe_.get_offset() + qe_.get_from();
            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_mismatch_scores_[index] = mismatch(a, b);

                }
            }

        }

        virtual ~MaldeScoresB() {
            delete[] weighted_mismatch_scores_;
        }

        MaldeScoresB(const MaldeScoresB& other) : // copy constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                mismatch_(other.mismatch_), gap_open_(other.gap_open_), gap_extend_(other.gap_extend_),
                score_shift_(other.score_shift_), num_scores_(other.num_scores_) {

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

        }

        MaldeScoresB(MaldeScoresB&& other) noexcept : // move constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                mismatch_(other.mismatch_), gap_open_(other.gap_open_), gap_extend_(other.gap_extend_),
                weighted_mismatch_scores_(other.eighted_mismatch_scores_), score_shift_(other.score_shift_),
                num_scores_(other.num_scores_) {

            other.weighted_mismatch_scores_ = nullptr;

        }

        MaldeScoresB& operator=(const MaldeScoresB& other)  { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

            return *this;

        }

        MaldeScoresB& operator=(MaldeScoresB&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = other.weighted_mismatch_scores_;
            other.weighted_mismatch_scores_ = nullptr;

            return *this;

        }


        virtual inline float match(char qual_a, char qual_b) {
            return 0;
        }

        inline float mismatch(char qual_a, char qual_b) {

            float err_prob_a = bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float err_prob_b = bInner_.boost(qe_.encoded_quality_to_probability(qual_b));
            float comb_err = bOuter_.boost(err_prob_a + err_prob_b - (1.0f + 1.0f / (alphabet_size_ - 1)) * err_prob_a * err_prob_b);

            return (1 - comb_err / (alphabet_size_ - 1)) * mismatch_;
        }

        inline float substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? match(qual_a, qual_b) : mismatch(qual_a, qual_b);
        }

        inline float gap_open(char qual_a) { return gap_open_; }
        inline float gap_extend(char qual_a) { return gap_extend_; }


        virtual inline float lookup_match(char qual_a, char qual_b) {
            return 0;
        }

        inline float lookup_mismatch(char qual_a, char qual_b) {
            return weighted_mismatch_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

        inline float lookup_substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? lookup_match(qual_a, qual_b) : lookup_mismatch(qual_a, qual_b);
        }

        inline float lookup_gap_open(char qual_a) { return gap_open_; }
        inline float lookup_gap_extend(char qual_a) { return gap_extend_; }

    protected:
        Bi bInner_; // inner boosting function (for individual error probabilities)
        Bo bOuter_; // outer boosting function (for combined error probabilities)
        QualityEncoding<> qe_; // encoding of the quality score
        lenSeqs_t alphabet_size_; // size of the used sequence alphabet

        float mismatch_; // mismatch penalty
        float gap_open_; // gap-opening penalty
        float gap_extend_; // gap-extension penalty

        float* weighted_mismatch_scores_; // matrix of weighted mismatch scores as 1-dimensional C-style array
        char score_shift_; // (offset + from) of the quality encoding to normalise the quality scores
        unsigned long num_scores_; // number of levels in the quality encoding

    };

    template<typename Bi, typename Bo>
    class MaldeScoresWithMatchesB : public MaldeScoresB<Bi, Bo> {

        using MaldeScoresB<Bi, Bo>::bInner_;
        using MaldeScoresB<Bi, Bo>::bOuter_;
        using MaldeScoresB<Bi, Bo>::qe_;
        using MaldeScoresB<Bi, Bo>::alphabet_size_;
        using MaldeScoresB<Bi, Bo>::mismatch_;
        using MaldeScoresB<Bi, Bo>::score_shift_;
        using MaldeScoresB<Bi, Bo>::num_scores_;

        using MaldeScoresB<Bi, Bo>::mismatch;
        using MaldeScoresB<Bi, Bo>::lookup_mismatch;

    public:
        MaldeScoresWithMatchesB(const QualityEncoding<>& qe, const Configuration& config) : MaldeScoresB<Bi, Bo>(qe, config) {

            // also fill matrix of weighted match scores by weighting the mismatch penalty with the combined error probability
            auto chars = qe_.get_accepted_scores();
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_match_scores_[index] = match(a, b);

                }
            }

        }

        ~MaldeScoresWithMatchesB() {
            delete[] weighted_match_scores_;
        }

        MaldeScoresWithMatchesB(const MaldeScoresWithMatchesB& other) : // copy constructor
                MaldeScoresB<Bi, Bo>(other) {

            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

        }

        MaldeScoresWithMatchesB(MaldeScoresWithMatchesB&& other) noexcept : // move constructor
                MaldeScoresB<Bi, Bo>(std::move(other)), weighted_match_scores_(other.weighted_match_scores_) {

            other.weighted_match_scores_ = nullptr;

        }

        MaldeScoresWithMatchesB& operator=(const MaldeScoresWithMatchesB& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            MaldeScoresB<Bi, Bo>::operator=(other);

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

            return *this;

        }

        MaldeScoresWithMatchesB& operator=(MaldeScoresWithMatchesB&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            MaldeScoresB<Bi, Bo>::operator=(std::move(other));

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = other.weighted_match_scores_;
            other.weighted_match_scores_ = nullptr;

            return *this;

        }


        inline float match(char qual_a, char qual_b) override {

            float err_prob_a = bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float err_prob_b = bInner_.boost(qe_.encoded_quality_to_probability(qual_b));
            float comb_err = bOuter_.boost(err_prob_a + err_prob_b - (1.0f + 1.0f / (alphabet_size_ - 1)) * err_prob_a * err_prob_b);

            return comb_err * mismatch_;

        }

        inline float lookup_match(char qual_a, char qual_b) override {
            return weighted_match_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

    protected:
        float* weighted_match_scores_; // matrix of weighted match scores as 1-dimensional C-style array

    };


    // Version C: use (combined) error simply as weight (also for gaps)

    /*
     * Provides quality-weighted scores using the mechanism from:
     * "The effect of sequence quality on sequence alignment"
     * (Malde, 2008, https://doi.org/10.1093/bioinformatics/btn052)
     *
     * Affects all edit operations (weights the all penalties with the error rates).
     *
     * The costs of matches are not quality-weighted (see MaldeScoresWithMatchesC for a version also weighting matches).
     */
    template<typename Bi, typename Bo>
    class MaldeScoresC {

    public:
        MaldeScoresC(const QualityEncoding<>& qe, const Configuration& config) : bInner_(config), bOuter_(config), qe_(qe) {

            ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
            mismatch_ = sf.penalty_mismatch;
            gap_open_ = sf.penalty_open;
            gap_extend_ = sf.penalty_extend;

            alphabet_size_ = (config.alphabet.empty()) ? 4 : config.alphabet.length();

            // fill matrix of weighted mismatch scores by weighting the mismatch penalty with the combined error probability
            // and fill matrices of weighted gap costs by weighting the gap penalties with the error probability
            auto chars = qe_.get_accepted_scores();
            num_scores_ = chars.size();
            score_shift_ = qe_.get_offset() + qe_.get_from();
            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            weighted_gap_open_scores_ = new float[num_scores_];
            weighted_gap_extend_scores_ = new float[num_scores_];
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_mismatch_scores_[index] = mismatch(a, b);

                }

                weighted_gap_open_scores_[a - score_shift_] = gap_open(a);
                weighted_gap_extend_scores_[a - score_shift_] = gap_extend(a);

            }

        }

        virtual ~MaldeScoresC() {

            delete[] weighted_gap_extend_scores_;
            delete[] weighted_gap_open_scores_;
            delete[] weighted_mismatch_scores_;

        }

        MaldeScoresC(const MaldeScoresC& other) : // copy constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                mismatch_(other.mismatch_), gap_open_(other.gap_open_), gap_extend_(other.gap_extend_),
                score_shift_(other.score_shift_), num_scores_(other.num_scores_) {

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

            weighted_gap_open_scores_ = new float[num_scores_];
            weighted_gap_extend_scores_ = new float[num_scores_];
            for (unsigned long i = 0; i < num_scores_; i++) {

                weighted_gap_open_scores_[i] = other.weighted_gap_open_scores_[i];
                weighted_gap_extend_scores_[i] = other.weighted_gap_extend_scores_[i];

            }

        }

        MaldeScoresC(MaldeScoresC&& other) noexcept : // move constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                mismatch_(other.mismatch_), gap_open_(other.gap_open_), gap_extend_(other.gap_extend_),
                weighted_mismatch_scores_(other.weighted_mismatch_scores_),
                weighted_gap_open_scores_(other.weighted_gap_open_scores_),
                weighted_gap_extend_scores_(other.weighted_gap_extend_scores_),
                score_shift_(other.score_shift_), num_scores_(other.num_scores_) {

            other.weighted_mismatch_scores_ = nullptr;
            other.weighted_gap_open_scores_ = nullptr;
            other.weighted_gap_extend_scores_ = nullptr;

        }

        MaldeScoresC& operator=(const MaldeScoresC& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;
            delete[] weighted_gap_open_scores_;
            delete[] weighted_gap_extend_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

            weighted_gap_open_scores_ = new float[num_scores_];
            weighted_gap_extend_scores_ = new float[num_scores_];
            for (unsigned long i = 0; i < num_scores_; i++) {

                weighted_gap_open_scores_[i] = other.weighted_gap_open_scores_[i];
                weighted_gap_extend_scores_[i] = other.weighted_gap_extend_scores_[i];

            }

            return *this;

        }

        MaldeScoresC& operator=(MaldeScoresC&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;
            delete[] weighted_gap_open_scores_;
            delete[] weighted_gap_extend_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = other.weighted_mismatch_scores_;
            other.weighted_mismatch_scores_ = nullptr;
            weighted_gap_open_scores_ = other.weighted_gap_open_scores_;
            other.weighted_gap_open_scores_ = nullptr;
            weighted_gap_extend_scores_ = other.weighted_gap_extend_scores_;
            other.weighted_gap_extend_scores_ = nullptr;

            return *this;

        }


        virtual inline float match(char qual_a, char qual_b) {
            return 0;
        }

        inline float mismatch(char qual_a, char qual_b) {

            float err_prob_a = bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float err_prob_b = bInner_.boost(qe_.encoded_quality_to_probability(qual_b));
            float comb_err = bOuter_.boost(err_prob_a + err_prob_b - (1.0f + 1.0f / (alphabet_size_ - 1)) * err_prob_a * err_prob_b);

            return (1 - comb_err / (alphabet_size_ - 1)) * mismatch_;

        }

        inline float substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? match(qual_a, qual_b) : mismatch(qual_a, qual_b);
        }

        inline float gap_open(char qual_a) {
            return gap_open_ * (1 - bOuter_.boost(bInner_.boost(qe_.encoded_quality_to_probability(qual_a))));
        }

        inline float gap_extend(char qual_a) {
            return gap_extend_ * (1 - bOuter_.boost(bInner_.boost(qe_.encoded_quality_to_probability(qual_a))));
        }


        virtual inline float lookup_match(char qual_a, char qual_b) {
            return 0;
        }

        inline float lookup_mismatch(char qual_a, char qual_b) {
            return weighted_mismatch_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

        inline float lookup_substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? lookup_match(qual_a, qual_b) : lookup_mismatch(qual_a, qual_b);
        }

        inline float lookup_gap_open(char qual_a) {
            return weighted_gap_open_scores_[qual_a - score_shift_];
        }

        inline float lookup_gap_extend(char qual_a) {
            return weighted_gap_extend_scores_[qual_a - score_shift_];
        }

    protected:
        Bi bInner_; // inner boosting function (for individual error probabilities)
        Bo bOuter_; // outer boosting function (for combined error probabilities)
        QualityEncoding<> qe_; // encoding of the quality score
        lenSeqs_t alphabet_size_; // size of the used sequence alphabet

        float mismatch_; // mismatch penalty
        float gap_open_; // gap-opening penalty
        float gap_extend_; // gap-extension penalty

        float* weighted_mismatch_scores_; // matrix of weighted mismatch scores as 1-dimensional C-style array
        float* weighted_gap_open_scores_; // matrix of weighted gap-opening scores as 1-dimensional C-style array
        float* weighted_gap_extend_scores_; // matrix of weighted gap-extension scores as 1-dimensional C-style array
        char score_shift_; // (offset + from) of the quality encoding to normalise the quality scores
        unsigned long num_scores_; // number of levels in the quality encoding

    };

    template<typename Bi, typename Bo>
    class MaldeScoresWithMatchesC : public MaldeScoresC<Bi, Bo> {

        using MaldeScoresC<Bi, Bo>::bInner_;
        using MaldeScoresC<Bi, Bo>::bOuter_;
        using MaldeScoresC<Bi, Bo>::qe_;
        using MaldeScoresC<Bi, Bo>::alphabet_size_;
        using MaldeScoresC<Bi, Bo>::mismatch_;
        using MaldeScoresC<Bi, Bo>::score_shift_;
        using MaldeScoresC<Bi, Bo>::num_scores_;

        using MaldeScoresC<Bi, Bo>::mismatch;
        using MaldeScoresC<Bi, Bo>::lookup_mismatch;

    public:
        MaldeScoresWithMatchesC(const QualityEncoding<>& qe, const Configuration& config) : MaldeScoresC<Bi, Bo>(qe, config) {

            // also fill matrix of weighted match scores by weighting the mismatch penalty with the combined error probability
            auto chars = qe_.get_accepted_scores();
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_match_scores_[index] = match(a, b);

                }
            }

        }

        ~MaldeScoresWithMatchesC() {
            delete[] weighted_match_scores_;
        }

        MaldeScoresWithMatchesC(const MaldeScoresWithMatchesC& other) : // copy constructor
                MaldeScoresC<Bi, Bo>(other) {

            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

        }

        MaldeScoresWithMatchesC(MaldeScoresWithMatchesC&& other) noexcept : // move constructor
                MaldeScoresC<Bi, Bo>(std::move(other)), weighted_match_scores_(other.weighted_match_scores_) {

            other.weighted_match_scores_ = nullptr;

        }

        MaldeScoresWithMatchesC& operator=(const MaldeScoresWithMatchesC& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            MaldeScoresC<Bi, Bo>::operator=(other);

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

            return *this;

        }

        MaldeScoresWithMatchesC& operator=(MaldeScoresWithMatchesC&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            MaldeScoresC<Bi, Bo>::operator=(std::move(other));

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = other.weighted_match_scores_;
            other.weighted_match_scores_ = nullptr;

            return *this;

        }


        inline float match(char qual_a, char qual_b) override {

            float err_prob_a = bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float err_prob_b = bInner_.boost(qe_.encoded_quality_to_probability(qual_b));
            float comb_err = bOuter_.boost(err_prob_a + err_prob_b - (1.0f + 1.0f / (alphabet_size_ - 1)) * err_prob_a * err_prob_b);

            return comb_err * mismatch_;

        }

        inline float lookup_match(char qual_a, char qual_b) override {
            return weighted_match_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

    protected:
        float* weighted_match_scores_; // matrix of weighted match scores as 1-dimensional C-style array

    };


    /* === Frith et al., Incorporating sequence quality data into alignment improves DNA read mapping === */

    /*
     * Provides quality-weighted scores using the mechanism from:
     * "Incorporating sequence quality data into alignment improves DNA read mapping"
     * (Frith et al., 2010, https://doi.org/10.1093/nar/gkq010)
     *
     * Affects only substitutions (modifies the mismatch penalty).
     *
     * The costs of matches are not quality-weighted (see FrithScoresWithMatches for a version also weighting matches).
     */
    template<typename Bi, typename Bo>
    class FrithScores {

    public:
        FrithScores(const QualityEncoding<>& qe, const Configuration& config) : bInner_(config), bOuter_(config), qe_(qe) {

            ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
            mismatch_ = sf.penalty_mismatch;
            gap_open_ = sf.penalty_open;
            gap_extend_ = sf.penalty_extend;

            alphabet_size_ = (config.alphabet.empty()) ? 4 : config.alphabet.length();
            scaling_factor_ = 5.22113;
            if (config.misc.find("scaling_factor") != config.misc.end()) scaling_factor_ = std::stof(config.misc.find("scaling_factor")->second);

            // compute contents of modified scoring matrix by incorporating the error probabilities
            // and record the minimum and maximum (mis)match scores
            auto chars = qe_.get_accepted_scores();
            num_scores_ = chars.size();
            score_shift_ = qe_.get_offset() + qe_.get_from();
            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            orig_match_ratio_ = std::exp(config.match_reward / scaling_factor_);
            orig_mismatch_ratio_ = std::exp(config.mismatch_penalty / scaling_factor_);
            min_match_ = std::numeric_limits<float>::infinity();
            max_match_ = -std::numeric_limits<float>::infinity();
            min_mismatch_ = std::numeric_limits<float>::infinity();
            max_mismatch_ = -std::numeric_limits<float>::infinity();
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);

                    auto tmp = match_helper(a, b);
                    min_match_ = std::min(min_match_, tmp);
                    max_match_ = std::max(max_match_, tmp);

                    weighted_mismatch_scores_[index] = mismatch_helper(a, b);
                    min_mismatch_ = std::min(min_mismatch_, weighted_mismatch_scores_[index]);
                    max_mismatch_ = std::max(max_mismatch_, weighted_mismatch_scores_[index]);

                }
            }

            // using the minimum and maximum scores, map the score range of the matrix
            // to the score range of the minimising scoring function (only for mismatch)
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_mismatch_scores_[index] = (max_mismatch_ - weighted_mismatch_scores_[index]) / (max_mismatch_ - min_mismatch_) * mismatch_;

                }
            }

        }

        virtual ~FrithScores() {
            delete[] weighted_mismatch_scores_;
        }

        FrithScores(const FrithScores& other) : // copy constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                scaling_factor_(other.scaling_factor_), mismatch_(other.mismatch_), gap_open_(other.gap_open_),
                gap_extend_(other.gap_extend_), orig_match_ratio_(other.orig_match_ratio_),
                orig_mismatch_ratio_(other.orig_mismatch_ratio_), min_match_(other.min_match_),
                max_match_(other.max_match_), min_mismatch_(other.min_mismatch_), max_mismatch_(other.max_mismatch_),
                score_shift_(other.score_shift_), num_scores_(other.num_scores_) {

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

        }

        FrithScores(FrithScores&& other) noexcept : // move constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                scaling_factor_(other.scaling_factor_), mismatch_(other.mismatch_), gap_open_(other.gap_open_),
                gap_extend_(other.gap_extend_), orig_match_ratio_(other.orig_match_ratio_),
                orig_mismatch_ratio_(other.orig_mismatch_ratio_), min_match_(other.min_match_),
                max_match_(other.max_match_), min_mismatch_(other.min_mismatch_), max_mismatch_(other.max_mismatch_),
                weighted_mismatch_scores_(other.weighted_mismatch_scores_),
                score_shift_(other.score_shift_), num_scores_(other.num_scores_) {

            other.weighted_mismatch_scores_ = nullptr;

        }

        FrithScores& operator=(const FrithScores& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            scaling_factor_ = other.scaling_factor_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            orig_match_ratio_ = other.orig_match_ratio_;
            orig_mismatch_ratio_ = other.orig_mismatch_ratio_;
            min_match_ = other.min_match_;
            max_match_ = other.max_match_;
            min_mismatch_ = other.min_mismatch_;
            max_mismatch_ = other.max_mismatch_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

            return *this;

        }

        FrithScores& operator=(FrithScores&& other) noexcept  { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            scaling_factor_ = other.scaling_factor_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            orig_match_ratio_ = other.orig_match_ratio_;
            orig_mismatch_ratio_ = other.orig_mismatch_ratio_;
            min_match_ = other.min_match_;
            max_match_ = other.max_match_;
            min_mismatch_ = other.min_mismatch_;
            max_mismatch_ = other.max_mismatch_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = other.weighted_mismatch_scores_;
            other.weighted_mismatch_scores_ = nullptr;

            return *this;

        }


        virtual inline float match(char qual_a, char qual_b) {
            return 0;
        }

        inline float mismatch(char qual_a, char qual_b) {
            return (max_mismatch_ - mismatch_helper(qual_a, qual_b)) / (max_mismatch_ - min_mismatch_) * mismatch_;
        }

        inline float substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? match(qual_a, qual_b) : mismatch(qual_a, qual_b);
        }

        inline float gap_open(char qual_a) { return gap_open_; }
        inline float gap_extend(char qual_a) { return gap_extend_; }


        virtual inline float lookup_match(char qual_a, char qual_b) {
            return 0;
        }

        inline float lookup_mismatch(char qual_a, char qual_b) {
            return weighted_mismatch_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

        inline float lookup_substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? lookup_match(qual_a, qual_b) : lookup_mismatch(qual_a, qual_b);
        }

        inline float lookup_gap_open(char qual_a) { return gap_open_; }
        inline float lookup_gap_extend(char qual_a) { return gap_extend_; }

    protected:
        FrithScores() = default; // only for constructors of derived classes

        // helper function encapsulating the incorporation of the error probabilties into the match score
        inline float match_helper(char qual_a, char qual_b) {

            float prob_a = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float prob_b = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_b));
            float prob_match = bOuter_.boost((1.0f - prob_a - prob_b + alphabet_size_ * prob_a * prob_b) / (alphabet_size_ - 1));

            return scaling_factor_ * std::log(prob_match * orig_match_ratio_ + (1 - prob_match) * orig_mismatch_ratio_);

        }

        // helper function encapsulating the incorporation of the error probabilties into the mismatch score
        inline float mismatch_helper(char qual_a, char qual_b) {

            float prob_a = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float prob_b = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_b));
            float prob_match = bOuter_.boost((prob_a + prob_b - prob_a * prob_b * alphabet_size_ + alphabet_size_ - 2) / ((alphabet_size_ - 1) * (alphabet_size_ - 1)));

            return scaling_factor_ * std::log(prob_match * orig_match_ratio_ + (1 - prob_match) * orig_mismatch_ratio_);

        }

        Bi bInner_; // inner boosting function (for individual error probabilities)
        Bo bOuter_; // outer boosting function (for combined error probabilities)
        QualityEncoding<> qe_; // encoding of the quality score
        lenSeqs_t alphabet_size_; // size of the used sequence alphabet
        float scaling_factor_;

        float mismatch_; // mismatch penalty
        float gap_open_; // gap-opening penalty
        float gap_extend_; // gap-extension penalty
        float orig_match_ratio_; // likelihood ratio of match
        float orig_mismatch_ratio_; // likelihood ratio of mismatch
        float min_match_; // minimum match score in the alternative scoring matrix
        float max_match_; // maximum match score in the alternative scoring matrix
        float min_mismatch_; // minimum mismatch score in the alternative scoring matrix
        float max_mismatch_; // maximum mismatch score in the alternative scoring matrix

        float* weighted_mismatch_scores_; // matrix of weighted mismatch scores as 1-dimensional C-style array
        char score_shift_; // (offset + from) of the quality encoding to normalise the quality scores
        unsigned long num_scores_; // number of levels in the quality encoding

    };

    template<typename Bi, typename Bo>
    class FrithScoresWithMatches : public FrithScores<Bi, Bo> {

        using FrithScores<Bi, Bo>::bInner_;
        using FrithScores<Bi, Bo>::bOuter_;
        using FrithScores<Bi, Bo>::qe_;
        using FrithScores<Bi, Bo>::alphabet_size_;
        using FrithScores<Bi, Bo>::scaling_factor_;
        using FrithScores<Bi, Bo>::mismatch_;
        using FrithScores<Bi, Bo>::gap_open_;
        using FrithScores<Bi, Bo>::gap_extend_;
        using FrithScores<Bi, Bo>::orig_match_ratio_;
        using FrithScores<Bi, Bo>::orig_mismatch_ratio_;
        using FrithScores<Bi, Bo>::min_match_;
        using FrithScores<Bi, Bo>::max_match_;
        using FrithScores<Bi, Bo>::min_mismatch_;
        using FrithScores<Bi, Bo>::max_mismatch_;
        using FrithScores<Bi, Bo>::weighted_mismatch_scores_;
        using FrithScores<Bi, Bo>::score_shift_;
        using FrithScores<Bi, Bo>::num_scores_;

        using FrithScores<Bi, Bo>::match_helper;
        using FrithScores<Bi, Bo>::mismatch;
        using FrithScores<Bi, Bo>::mismatch_helper;
        using FrithScores<Bi, Bo>::lookup_mismatch;

    public:
        FrithScoresWithMatches(const QualityEncoding<>& qe, const Configuration& config) {

            bInner_ = Bi(config);
            bOuter_ = Bo(config);
            qe_ = qe;

            ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
            mismatch_ = sf.penalty_mismatch;
            gap_open_ = sf.penalty_open;
            gap_extend_ = sf.penalty_extend;

            alphabet_size_ = (config.alphabet.empty()) ? 4 : config.alphabet.length();
            scaling_factor_ = 5.22113;
            if (config.misc.find("scaling_factor") != config.misc.end()) scaling_factor_ = std::stof(config.misc.find("scaling_factor")->second);

            // compute contents of modified scoring matrix by incorporating the error probabilities
            // and record the minimum and maximum (mis)match scores
            auto chars = qe_.get_accepted_scores();
            num_scores_ = chars.size();
            score_shift_ = qe_.get_offset() + qe_.get_from();
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            orig_match_ratio_ = std::exp(config.match_reward / scaling_factor_);
            orig_mismatch_ratio_ = std::exp(config.mismatch_penalty / scaling_factor_);
            min_match_ = std::numeric_limits<float>::infinity();
            max_match_ = -std::numeric_limits<float>::infinity();
            min_mismatch_ = std::numeric_limits<float>::infinity();
            max_mismatch_ = -std::numeric_limits<float>::infinity();
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);

                    weighted_match_scores_[index] = match_helper(a, b);
                    weighted_mismatch_scores_[index] = mismatch_helper(a, b);

                    min_match_ = std::min(min_match_, weighted_match_scores_[index]);
                    max_match_ = std::max(max_match_, weighted_match_scores_[index]);
                    min_mismatch_ = std::min(min_mismatch_, weighted_mismatch_scores_[index]);
                    max_mismatch_ = std::max(max_mismatch_, weighted_mismatch_scores_[index]);

                }
            }

            // using the minimum and maximum scores, map the score range of the matrix
            // to the score range of the minimising scoring function (match and mismatch)
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);

                    weighted_match_scores_[index] = (max_match_ - weighted_match_scores_[index]) / (max_match_ - min_match_) * mismatch_;
                    weighted_mismatch_scores_[index] = (max_mismatch_ - weighted_mismatch_scores_[index]) / (max_mismatch_ - min_mismatch_) * mismatch_;

                }
            }

        }

        ~FrithScoresWithMatches() {
            delete[] weighted_match_scores_;
        }

        FrithScoresWithMatches(const FrithScoresWithMatches& other) : // copy constructor
                FrithScores<Bi, Bo>(other) {

            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

        }

        FrithScoresWithMatches(FrithScoresWithMatches&& other) noexcept : // move constructor
                FrithScores<Bi, Bo>(std::move(other)), weighted_match_scores_(other.weighted_match_scores_) {

            other.weighted_match_scores_ = nullptr;

        }

        FrithScoresWithMatches& operator=(const FrithScoresWithMatches& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            FrithScores<Bi, Bo>::operator=(other);

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

            return *this;

        }

        FrithScoresWithMatches& operator=(FrithScoresWithMatches&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            FrithScores<Bi, Bo>::operator=(std::move(other));

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = other.weighted_match_scores_;
            other.weighted_match_scores_ = nullptr;

            return *this;

        }


        inline float match(char qual_a, char qual_b) override {
            return (max_match_ - match_helper(qual_a, qual_b)) / (max_match_ - min_match_) * mismatch_;
        }

        inline float lookup_match(char qual_a, char qual_b) override {
            return weighted_match_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

    protected:
        float* weighted_match_scores_; // matrix of weighted match scores as 1-dimensional C-style array

    };


    /* === Kim et al., A DNA sequence alignment algorithm using quality information and a fuzzy inference method === */

    // Version A: approximated (as proposed in paper)

    /*
     * Provides quality-weighted scores using the mechanism from:
     * "A DNA sequence alignment algorithm using quality information and a fuzzy inference method"
     * (Kim et al., 2008, https://doi.org/10.1016/j.pnsc.2007.12.011)
     *
     * Affects all edit operations (modifies the penalties through a probability-weighted linear combination).
     *
     * The costs of matches are not quality-weighted (see KimScoresWithMatchesA for a version also weighting matches).
     */
    template<typename Bi, typename Bo>
    class KimScoresA {

    public:
        KimScoresA(const QualityEncoding<>& qe, const Configuration& config) : bInner_(config), bOuter_(config), qe_(qe) {

            ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
            mismatch_ = sf.penalty_mismatch;
            gap_open_ = sf.penalty_open;
            gap_extend_ = sf.penalty_extend;

            alphabet_size_ = (config.alphabet.empty()) ? 4 : config.alphabet.length();

            // precompute the weighted mismatch, gap-opening and gap-extension scores
            auto chars = qe_.get_accepted_scores();
            num_scores_ = chars.size();
            score_shift_ = qe_.get_offset() + qe_.get_from();
            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            weighted_gap_open_scores_ = new float[num_scores_];
            weighted_gap_extend_scores_ = new float[num_scores_];
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_mismatch_scores_[index] = mismatch(a, b);

                }

                weighted_gap_open_scores_[a - score_shift_] = gap_open(a);
                weighted_gap_extend_scores_[a - score_shift_] = gap_extend(a);

            }

        }

        virtual ~KimScoresA() {

            delete[] weighted_gap_extend_scores_;
            delete[] weighted_gap_open_scores_;
            delete[] weighted_mismatch_scores_;

        }

        KimScoresA(const KimScoresA& other) : // copy constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                mismatch_(other.mismatch_), gap_open_(other.gap_open_), gap_extend_(other.gap_extend_),
                score_shift_(other.score_shift_), num_scores_(other.num_scores_) {

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

            weighted_gap_open_scores_ = new float[num_scores_];
            weighted_gap_extend_scores_ = new float[num_scores_];
            for (unsigned long i = 0; i < num_scores_; i++) {

                weighted_gap_open_scores_[i] = other.weighted_gap_open_scores_[i];
                weighted_gap_extend_scores_[i] = other.weighted_gap_extend_scores_[i];

            }

        }

        KimScoresA(KimScoresA&& other) noexcept  : // move constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                mismatch_(other.mismatch_), gap_open_(other.gap_open_), gap_extend_(other.gap_extend_),
                weighted_mismatch_scores_(other.weighted_mismatch_scores_),
                weighted_gap_open_scores_(other.weighted_gap_open_scores_),
                weighted_gap_extend_scores_(other.weighted_gap_extend_scores_),
                score_shift_(other.score_shift_), num_scores_(other.num_scores_) {

            other.weighted_mismatch_scores_ = nullptr;
            other.weighted_gap_open_scores_ = nullptr;
            other.weighted_gap_extend_scores_ = nullptr;

        }

        KimScoresA& operator=(const KimScoresA& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;
            delete[] weighted_gap_open_scores_;
            delete[] weighted_gap_extend_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

            weighted_gap_open_scores_ = new float[num_scores_];
            weighted_gap_extend_scores_ = new float[num_scores_];
            for (unsigned long i = 0; i < num_scores_; i++) {

                weighted_gap_open_scores_[i] = other.weighted_gap_open_scores_[i];
                weighted_gap_extend_scores_[i] = other.weighted_gap_extend_scores_[i];

            }

            return *this;

        }

        KimScoresA& operator=(KimScoresA&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;
            delete[] weighted_gap_open_scores_;
            delete[] weighted_gap_extend_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = other.weighted_mismatch_scores_;
            other.weighted_mismatch_scores_ = nullptr;
            weighted_gap_open_scores_ = other.weighted_gap_open_scores_;
            other.weighted_gap_open_scores_ = nullptr;
            weighted_gap_extend_scores_ = other.weighted_gap_extend_scores_;
            other.weighted_gap_extend_scores_ = nullptr;

            return *this;

        }


        virtual inline float match(char qual_a, char qual_b) {
            return 0;
        }

        inline float mismatch(char qual_a, char qual_b) {

            float prob_a = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float prob_b = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_b));

            //float prob_match = bOuter_.boost((1 - prob_a) * prob_b / alphabet_size_ + prob_a * (1 - prob_b) / alphabet_size_);
            float prob_mismatch = bOuter_.boost(prob_a * prob_b + (alphabet_size_ - 2) * (1 - prob_a) * prob_b / alphabet_size_ + (alphabet_size_ - 2) * prob_a * (1 - prob_b) / alphabet_size_);
            float prob_gap = bOuter_.boost((1 - prob_a) * prob_b / alphabet_size_ + prob_a * (1 - prob_b) / alphabet_size_);

            return /* 0 * prob_match + */mismatch_ * prob_mismatch + (gap_open_ + gap_extend_) * prob_gap;

        }

        inline float substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? match(qual_a, qual_b) : mismatch(qual_a, qual_b);
        }

        inline float gap_open(char qual_a) {
            return gap_open_ * bOuter_.boost(1 - bInner_.boost(qe_.encoded_quality_to_probability(qual_a)) / alphabet_size_);
        }

        inline float gap_extend(char qual_a) {
            return gap_extend_ * bOuter_.boost(1 - bInner_.boost(qe_.encoded_quality_to_probability(qual_a)) / alphabet_size_);
        }


        virtual inline float lookup_match(char qual_a, char qual_b) {
            return 0;
        }

        inline float lookup_mismatch(char qual_a, char qual_b) {
            return weighted_mismatch_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

        inline float lookup_substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? lookup_match(qual_a, qual_b) : lookup_mismatch(qual_a, qual_b);
        }

        inline float lookup_gap_open(char qual_a) {
            return weighted_gap_open_scores_[qual_a - score_shift_];
        }

        inline float lookup_gap_extend(char qual_a) {
            return weighted_gap_extend_scores_[qual_a - score_shift_];
        }

    protected:
        Bi bInner_; // inner boosting function (for individual error probabilities)
        Bo bOuter_; // outer boosting function (for combined error probabilities)
        QualityEncoding<> qe_; // encoding of the quality score
        lenSeqs_t alphabet_size_; // size of the used sequence alphabet

        float mismatch_; // mismatch penalty
        float gap_open_; // gap-opening penalty
        float gap_extend_; // gap-extension penalty

        float* weighted_mismatch_scores_; // matrix of weighted mismatch scores as 1-dimensional C-style array
        float* weighted_gap_open_scores_; // matrix of weighted gap-opening scores as 1-dimensional C-style array
        float* weighted_gap_extend_scores_; // matrix of weighted gap-extension scores as 1-dimensional C-style array
        char score_shift_; // (offset + from) of the quality encoding to normalise the quality scores
        unsigned long num_scores_; // number of levels in the quality encoding

    };

    template<typename Bi, typename Bo>
    class KimScoresWithMatchesA : public KimScoresA<Bi, Bo> {

        using KimScoresA<Bi, Bo>::bInner_;
        using KimScoresA<Bi, Bo>::bOuter_;
        using KimScoresA<Bi, Bo>::qe_;
        using KimScoresA<Bi, Bo>::alphabet_size_;
        using KimScoresA<Bi, Bo>::mismatch_;
        using KimScoresA<Bi, Bo>::gap_open_;
        using KimScoresA<Bi, Bo>::gap_extend_;
        using KimScoresA<Bi, Bo>::score_shift_;
        using KimScoresA<Bi, Bo>::num_scores_;

        using KimScoresA<Bi, Bo>::mismatch;
        using KimScoresA<Bi, Bo>::lookup_mismatch;

    public:
        KimScoresWithMatchesA(const QualityEncoding<>& qe, const Configuration& config) : KimScoresA<Bi, Bo>(qe, config) {

            // also precompute the weighted match scores
            auto chars = qe_.get_accepted_scores();
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_match_scores_[index] = match(a, b);

                }
            }

        }

        ~KimScoresWithMatchesA() {
            delete[] weighted_match_scores_;
        }

        KimScoresWithMatchesA(const KimScoresWithMatchesA& other) : // copy constructor
                KimScoresA<Bi, Bo>(other) {

            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

        }

        KimScoresWithMatchesA(KimScoresWithMatchesA&& other) noexcept : // move constructor
                KimScoresA<Bi, Bo>(std::move(other)), weighted_match_scores_(other.weighted_match_scores_) {

            other.weighted_match_scores_ = nullptr;

        }

        KimScoresWithMatchesA& operator=(const KimScoresWithMatchesA& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            KimScoresA<Bi, Bo>::operator=(other);

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

            return *this;

        }

        KimScoresWithMatchesA& operator=(KimScoresWithMatchesA&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            KimScoresA<Bi, Bo>::operator=(std::move(other));

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = other.weighted_match_scores_;
            other.weighted_match_scores_ = nullptr;

            return *this;

        }


        inline float match(char qual_a, char qual_b) override {

            float prob_a = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float prob_b = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_b));

            //float prob_match = bOuter_.boost(prob_a * prob_b);
            float prob_mismatch = bOuter_.boost((alphabet_size_ - 1) * (1 - prob_a) * prob_b / alphabet_size_ + (alphabet_size_ - 1) * prob_a * (1 - prob_b) / alphabet_size_);
            float prob_gap = bOuter_.boost((1 - prob_a) * prob_b / alphabet_size_ + prob_a * (1 - prob_b) / alphabet_size_);

            return /* 0 * prob_match + */mismatch_ * prob_mismatch + (gap_open_ + gap_extend_) * prob_gap;

        }

        inline float lookup_match(char qual_a, char qual_b) override {
            return weighted_match_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

    protected:
        float* weighted_match_scores_; // matrix of weighted match scores as 1-dimensional C-style array

    };


    // Version B: not approximated

    /*
     * Provides quality-weighted scores using the mechanism from:
     * "A DNA sequence alignment algorithm using quality information and a fuzzy inference method"
     * (Kim et al., 2008, https://doi.org/10.1016/j.pnsc.2007.12.011)
     *
     * Affects all edit operations (modifies the penalties through a probability-weighted linear combination).
     *
     * The costs of matches are not quality-weighted (see KimScoresWithMatchesB for a version also weighting matches).
     */
    template<typename Bi, typename Bo>
    class KimScoresB {

    public:
        KimScoresB(const QualityEncoding<>& qe, const Configuration& config) : bInner_(config), bOuter_(config), qe_(qe) {

            ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
            mismatch_ = sf.penalty_mismatch;
            gap_open_ = sf.penalty_open;
            gap_extend_ = sf.penalty_extend;

            alphabet_size_ = (config.alphabet.empty()) ? 4 : config.alphabet.length();

            // precompute the weighted mismatch, gap-opening and gap-extension scores
            auto chars = qe_.get_accepted_scores();
            num_scores_ = chars.size();
            score_shift_ = qe_.get_offset() + qe_.get_from();
            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            weighted_gap_open_scores_ = new float[num_scores_];
            weighted_gap_extend_scores_ = new float[num_scores_];
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_mismatch_scores_[index] = mismatch(a, b);

                }

                weighted_gap_open_scores_[a - score_shift_] = gap_open(a);
                weighted_gap_extend_scores_[a - score_shift_] = gap_extend(a);

            }

        }

        virtual ~KimScoresB() {

            delete[] weighted_gap_extend_scores_;
            delete[] weighted_gap_open_scores_;
            delete[] weighted_mismatch_scores_;

        }

        KimScoresB(const KimScoresB& other) : // copy constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                mismatch_(other.mismatch_), gap_open_(other.gap_open_), gap_extend_(other.gap_extend_),
                score_shift_(other.score_shift_), num_scores_(other.num_scores_) {

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

            weighted_gap_open_scores_ = new float[num_scores_];
            weighted_gap_extend_scores_ = new float[num_scores_];
            for (unsigned long i = 0; i < num_scores_; i++) {

                weighted_gap_open_scores_[i] = other.weighted_gap_open_scores_[i];
                weighted_gap_extend_scores_[i] = other.weighted_gap_extend_scores_[i];

            }

        }

        KimScoresB(KimScoresB&& other) noexcept  : // move constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                mismatch_(other.mismatch_), gap_open_(other.gap_open_), gap_extend_(other.gap_extend_),
                weighted_mismatch_scores_(other.weighted_mismatch_scores_),
                weighted_gap_open_scores_(other.weighted_gap_open_scores_),
                weighted_gap_extend_scores_(other.weighted_gap_extend_scores_),
                score_shift_(other.score_shift_), num_scores_(other.num_scores_) {

            other.weighted_mismatch_scores_ = nullptr;
            other.weighted_gap_open_scores_ = nullptr;
            other.weighted_gap_extend_scores_ = nullptr;

        }

        KimScoresB& operator=(const KimScoresB& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;
            delete[] weighted_gap_open_scores_;
            delete[] weighted_gap_extend_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

            weighted_gap_open_scores_ = new float[num_scores_];
            weighted_gap_extend_scores_ = new float[num_scores_];
            for (unsigned long i = 0; i < num_scores_; i++) {

                weighted_gap_open_scores_[i] = other.weighted_gap_open_scores_[i];
                weighted_gap_extend_scores_[i] = other.weighted_gap_extend_scores_[i];

            }

            return *this;

        }

        KimScoresB& operator=(KimScoresB&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;
            delete[] weighted_gap_open_scores_;
            delete[] weighted_gap_extend_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = other.weighted_mismatch_scores_;
            other.weighted_mismatch_scores_ = nullptr;
            weighted_gap_open_scores_ = other.weighted_gap_open_scores_;
            other.weighted_gap_open_scores_ = nullptr;
            weighted_gap_extend_scores_ = other.weighted_gap_extend_scores_;
            other.weighted_gap_extend_scores_ = nullptr;

            return *this;

        }


        virtual inline float match(char qual_a, char qual_b) {
            return 0;
        }

        inline float mismatch(char qual_a, char qual_b) {

            float prob_a = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float prob_b = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_b));

            //float prob_match = bOuter_.boost((1 - prob_a) * prob_b / alphabet_size_ + prob_a * (1 - prob_b) / alphabet_size_ + (alphabet_size_ - 2) * (1 - prob_a) * (1 - prob_b) / (alphabet_size_ * alphabet_size_));
            float prob_mismatch = bOuter_.boost(prob_a * prob_b + (alphabet_size_ - 2) * (1 - prob_a) * prob_b / alphabet_size_ + (alphabet_size_ - 2) * prob_a * (1 - prob_b) / alphabet_size_ + ((alphabet_size_ - 1) * (alphabet_size_ - 1) - alphabet_size_ + 2) * (1 - prob_a) * (1 - prob_b) / (alphabet_size_ * alphabet_size_));
            float prob_gap = bOuter_.boost((1 - prob_a) * prob_b / alphabet_size_ + prob_a * (1 - prob_b) / alphabet_size_ + 2 * (alphabet_size_ - 1) * (1 - prob_a) * (1 - prob_b) / (alphabet_size_ * alphabet_size_));

            return /* 0 * prob_match + */mismatch_ * prob_mismatch + (gap_open_ + gap_extend_) * prob_gap;

        }

        inline float substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? match(qual_a, qual_b) : mismatch(qual_a, qual_b);
        }

        inline float gap_open(char qual_a) {
            return gap_open_ * bOuter_.boost(1 - bInner_.boost(qe_.encoded_quality_to_probability(qual_a)) / alphabet_size_);
        }

        inline float gap_extend(char qual_a) {
            return gap_extend_ * bOuter_.boost(1 - bInner_.boost(qe_.encoded_quality_to_probability(qual_a)) / alphabet_size_);
        }


        virtual inline float lookup_match(char qual_a, char qual_b) {
            return 0;
        }

        inline float lookup_mismatch(char qual_a, char qual_b) {
            return weighted_mismatch_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

        inline float lookup_substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? lookup_match(qual_a, qual_b) : lookup_mismatch(qual_a, qual_b);
        }

        inline float lookup_gap_open(char qual_a) {
            return weighted_gap_open_scores_[qual_a - score_shift_];
        }

        inline float lookup_gap_extend(char qual_a) {
            return weighted_gap_extend_scores_[qual_a - score_shift_];
        }

    protected:
        Bi bInner_; // inner boosting function (for individual error probabilities)
        Bo bOuter_; // outer boosting function (for combined error probabilities)
        QualityEncoding<> qe_; // encoding of the quality score
        lenSeqs_t alphabet_size_; // size of the used sequence alphabet

        float mismatch_; // mismatch penalty
        float gap_open_; // gap-opening penalty
        float gap_extend_; // gap-extension penalty

        float* weighted_mismatch_scores_; // matrix of weighted mismatch scores as 1-dimensional C-style array
        float* weighted_gap_open_scores_; // matrix of weighted gap-opening scores as 1-dimensional C-style array
        float* weighted_gap_extend_scores_; // matrix of weighted gap-extension scores as 1-dimensional C-style array
        char score_shift_; // (offset + from) of the quality encoding to normalise the quality scores
        unsigned long num_scores_; // number of levels in the quality encoding

    };

    template<typename Bi, typename Bo>
    class KimScoresWithMatchesB : public KimScoresB<Bi, Bo> {

        using KimScoresB<Bi, Bo>::bInner_;
        using KimScoresB<Bi, Bo>::bOuter_;
        using KimScoresB<Bi, Bo>::qe_;
        using KimScoresB<Bi, Bo>::alphabet_size_;
        using KimScoresB<Bi, Bo>::mismatch_;
        using KimScoresB<Bi, Bo>::gap_open_;
        using KimScoresB<Bi, Bo>::gap_extend_;
        using KimScoresB<Bi, Bo>::score_shift_;
        using KimScoresB<Bi, Bo>::num_scores_;

        using KimScoresB<Bi, Bo>::mismatch;
        using KimScoresB<Bi, Bo>::lookup_mismatch;

    public:
        KimScoresWithMatchesB(const QualityEncoding<>& qe, const Configuration& config) : KimScoresB<Bi, Bo>(qe, config) {

            // also precompute the weighted match scores
            auto chars = qe_.get_accepted_scores();
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_match_scores_[index] = match(a, b);

                }
            }

        }

        ~KimScoresWithMatchesB() {
            delete[] weighted_match_scores_;
        }

        KimScoresWithMatchesB(const KimScoresWithMatchesB& other) : // copy constructor
                KimScoresB<Bi, Bo>(other) {

            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

        }

        KimScoresWithMatchesB(KimScoresWithMatchesB&& other) noexcept : // move constructor
                KimScoresB<Bi, Bo>(std::move(other)), weighted_match_scores_(other.weighted_match_scores_) {

            other.weighted_match_scores_ = nullptr;

        }

        KimScoresWithMatchesB& operator=(const KimScoresWithMatchesB& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            KimScoresB<Bi, Bo>::operator=(other);

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

            return *this;

        }

        KimScoresWithMatchesB& operator=(KimScoresWithMatchesB&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            KimScoresB<Bi, Bo>::operator=(std::move(other));

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = other.weighted_match_scores_;
            other.weighted_match_scores_ = nullptr;

            return *this;

        }


        inline float match(char qual_a, char qual_b) override {

            float prob_a = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float prob_b = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_b));

            //float prob_match = bOuter_.boost(prob_a * prob_b + (alphabet_size_ - 1) * (1 - prob_a) * (1 - prob_b) / (alphabet_size_ * alphabet_size_));
            float prob_mismatch = bOuter_.boost((alphabet_size_ - 1) * (1 - prob_a) * prob_b / alphabet_size_ + (alphabet_size_ - 1) * prob_a * (1 - prob_b) / alphabet_size_ + (alphabet_size_ - 1) * (alphabet_size_ - 2) * (1 - prob_a) * (1 - prob_b) / (alphabet_size_ * alphabet_size_));
            float prob_gap = bOuter_.boost((1 - prob_a) * prob_b / alphabet_size_ + prob_a * (1 - prob_b) / alphabet_size_ + 2 * (alphabet_size_ - 1) * (1 - prob_a) * (1 - prob_b) / (alphabet_size_ * alphabet_size_));

            return /* 0 * prob_match + */mismatch_ * prob_mismatch + (gap_open_ + gap_extend_) * prob_gap;

        }

        inline float lookup_match(char qual_a, char qual_b) override {
            return weighted_match_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

    protected:
        float* weighted_match_scores_; // matrix of weighted match scores as 1-dimensional C-style array

    };


    /* === Own idea: match reward and mismatch / gap penalties converge with decreasing quality === */

    // Version A: convergence to balance cost

    /*
     * Provides quality-weighted scores using the following mechanism:
     * Matches receive an increasingly large penalty with decreasing quality, while the penalty of mismatches and gaps
     * is increasingly reduced with decreasing quality. This "convergence" might reach a common value (balance cost)
     * at which matches, mismatches and gaps have the same cost (when the quality is 0).
     *
     * Affects all edit operations (modifies the penalties).
     *
     * The costs of matches are not quality-weighted (see ConvergeScoresWithMatchesA for a version also weighting matches).
     */
    template<typename Bi, typename Bo>
    class ConvergeScoresA {

    public:
        ConvergeScoresA(const QualityEncoding<>& qe, const Configuration& config) : bInner_(config), bOuter_(config), qe_(qe) {

            ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
            mismatch_ = sf.penalty_mismatch;
            gap_open_ = sf.penalty_open;
            gap_extend_ = sf.penalty_extend;

            alphabet_size_ = (config.alphabet.empty()) ? 4 : config.alphabet.length();

            // determine how much the cost of the each operation can change before reaching the balance cost
            float balance_factor = 0.5;
            if (config.misc.find("balance_factor") != config.misc.end()) balance_factor = std::stof(config.misc.find("balance_factor")->second);
            float converged = std::min({mismatch_, gap_open_, gap_extend_}) * balance_factor;
            max_change_match_ = converged;
            max_change_mismatch_ = mismatch_ - converged;
            max_change_gap_open_ = gap_open_ - converged;
            max_change_gap_extend_ = gap_extend_ - converged;

            // precompute the weighted mismatch, gap-opening and gap-extension scores
            auto chars = qe_.get_accepted_scores();
            num_scores_ = chars.size();
            score_shift_ = qe_.get_offset() + qe_.get_from();
            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            weighted_gap_open_scores_ = new float[num_scores_];
            weighted_gap_extend_scores_ = new float[num_scores_];
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_mismatch_scores_[index] = mismatch(a, b);

                }

                weighted_gap_open_scores_[a - score_shift_] = gap_open(a);
                weighted_gap_extend_scores_[a - score_shift_] = gap_extend(a);

            }

        }

        virtual ~ConvergeScoresA() {

            delete[] weighted_gap_extend_scores_;
            delete[] weighted_gap_open_scores_;
            delete[] weighted_mismatch_scores_;

        }

        ConvergeScoresA(const ConvergeScoresA& other) : // copy constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                max_change_match_(other.max_change_match_), max_change_mismatch_(other.max_change_mismatch_),
                max_change_gap_open_(other.max_change_gap_open_), max_change_gap_extend_(other.max_change_gap_extend_),
                mismatch_(other.mismatch_), gap_open_(other.gap_open_), gap_extend_(other.gap_extend_),
                score_shift_(other.score_shift_), num_scores_(other.num_scores_) {

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

            weighted_gap_open_scores_ = new float[num_scores_];
            weighted_gap_extend_scores_ = new float[num_scores_];
            for (unsigned long i = 0; i < num_scores_; i++) {

                weighted_gap_open_scores_[i] = other.weighted_gap_open_scores_[i];
                weighted_gap_extend_scores_[i] = other.weighted_gap_extend_scores_[i];

            }

        }

        ConvergeScoresA(ConvergeScoresA&& other) noexcept : // move constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                max_change_match_(other.max_change_match_), max_change_mismatch_(other.max_change_mismatch_),
                max_change_gap_open_(other.max_change_gap_open_), max_change_gap_extend_(other.max_change_gap_extend_),
                mismatch_(other.mismatch_), gap_open_(other.gap_open_), gap_extend_(other.gap_extend_),
                weighted_mismatch_scores_(other.weighted_mismatch_scores_),
                weighted_gap_open_scores_(other.weighted_gap_open_scores_),
                weighted_gap_extend_scores_(other.weighted_gap_extend_scores_),
                score_shift_(other.score_shift_), num_scores_(other.num_scores_) {

            other.weighted_mismatch_scores_ = nullptr;
            other.weighted_gap_open_scores_ = nullptr;
            other.weighted_gap_extend_scores_ = nullptr;

        }

        ConvergeScoresA& operator=(const ConvergeScoresA& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;
            delete[] weighted_gap_open_scores_;
            delete[] weighted_gap_extend_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            max_change_match_ = other.max_change_match_;
            max_change_mismatch_ = other.max_change_mismatch_;
            max_change_gap_open_ = other.max_change_gap_open_;
            max_change_gap_extend_ = other.max_change_gap_extend_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

            weighted_gap_open_scores_ = new float[num_scores_];
            weighted_gap_extend_scores_ = new float[num_scores_];
            for (unsigned long i = 0; i < num_scores_; i++) {

                weighted_gap_open_scores_[i] = other.weighted_gap_open_scores_[i];
                weighted_gap_extend_scores_[i] = other.weighted_gap_extend_scores_[i];

            }

            return *this;

        }

        ConvergeScoresA& operator=(ConvergeScoresA&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;
            delete[] weighted_gap_open_scores_;
            delete[] weighted_gap_extend_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            max_change_match_ = other.max_change_match_;
            max_change_mismatch_ = other.max_change_mismatch_;
            max_change_gap_open_ = other.max_change_gap_open_;
            max_change_gap_extend_ = other.max_change_gap_extend_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = other.weighted_mismatch_scores_;
            other.weighted_mismatch_scores_ = nullptr;
            weighted_gap_open_scores_ = other.weighted_gap_open_scores_;
            other.weighted_gap_open_scores_ = nullptr;
            weighted_gap_extend_scores_ = other.weighted_gap_extend_scores_;
            other.weighted_gap_extend_scores_ = nullptr;

            return *this;

        }


        virtual inline float match(char qual_a, char qual_b) {
            return 0;
        }

        inline float mismatch(char qual_a, char qual_b) {

            float prob_a = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float prob_b = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_b));

            return mismatch_ - max_change_mismatch_ * bOuter_.boost((alphabet_size_ - 2 + prob_a + prob_b - prob_a * prob_b * alphabet_size_) / ((alphabet_size_ - 1) * (alphabet_size_ - 1)));

        }

        inline float substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? match(qual_a, qual_b) : mismatch(qual_a, qual_b);
        }

        inline float gap_open(char qual_a) {
            return gap_open_ - max_change_gap_open_ * bOuter_.boost(bInner_.boost(qe_.encoded_quality_to_probability(qual_a)));
        }

        inline float gap_extend(char qual_a) {
            return gap_extend_ - max_change_gap_extend_ * bOuter_.boost(bInner_.boost(qe_.encoded_quality_to_probability(qual_a)));
        }


        virtual inline float lookup_match(char qual_a, char qual_b) {
            return 0;
        }

        inline float lookup_mismatch(char qual_a, char qual_b) {
            return weighted_mismatch_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

        inline float lookup_substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? lookup_match(qual_a, qual_b) : lookup_mismatch(qual_a, qual_b);
        }

        inline float lookup_gap_open(char qual_a) {
            return weighted_gap_open_scores_[qual_a - score_shift_];
        }

        inline float lookup_gap_extend(char qual_a) {
            return weighted_gap_extend_scores_[qual_a - score_shift_];
        }

    protected:
        Bi bInner_; // inner boosting function (for individual error probabilities)
        Bo bOuter_; // outer boosting function (for combined error probabilities)
        QualityEncoding<> qe_; // encoding of the quality score
        lenSeqs_t alphabet_size_; // size of the used sequence alphabet
        float max_change_match_; // difference between unweighted match cost and balance cost
        float max_change_mismatch_; // difference between unweighted mismatch cost and balance cost
        float max_change_gap_open_; // difference between unweighted gap-opening cost and balance cost
        float max_change_gap_extend_; // difference between unweighted gap-extension cost and balance cost

        float mismatch_; // mismatch penalty
        float gap_open_; // gap-opening penalty
        float gap_extend_; // gap-extension penalty

        float* weighted_mismatch_scores_; // matrix of weighted mismatch scores as 1-dimensional C-style array
        float* weighted_gap_open_scores_; // matrix of weighted gap-opening scores as 1-dimensional C-style array
        float* weighted_gap_extend_scores_; // matrix of weighted gap-extension scores as 1-dimensional C-style array
        char score_shift_; // (offset + from) of the quality encoding to normalise the quality scores
        unsigned long num_scores_; // number of levels in the quality encoding

    };

    template<typename Bi, typename Bo>
    class ConvergeScoresWithMatchesA : public ConvergeScoresA<Bi, Bo> {

        using ConvergeScoresA<Bi, Bo>::bInner_;
        using ConvergeScoresA<Bi, Bo>::bOuter_;
        using ConvergeScoresA<Bi, Bo>::qe_;
        using ConvergeScoresA<Bi, Bo>::alphabet_size_;
        using ConvergeScoresA<Bi, Bo>::max_change_match_;
        using ConvergeScoresA<Bi, Bo>::score_shift_;
        using ConvergeScoresA<Bi, Bo>::num_scores_;

        using ConvergeScoresA<Bi, Bo>::mismatch;
        using ConvergeScoresA<Bi, Bo>::lookup_mismatch;

    public:
        ConvergeScoresWithMatchesA(const QualityEncoding<>& qe, const Configuration& config) : ConvergeScoresA<Bi, Bo>(qe, config) {

            // also precompute the weighted match scores
            auto chars = qe_.get_accepted_scores();
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_match_scores_[index] = match(a, b);

                }
            }

        }

        ~ConvergeScoresWithMatchesA() {
            delete[] weighted_match_scores_;
        }

        ConvergeScoresWithMatchesA(const ConvergeScoresWithMatchesA& other) : // copy constructor
                ConvergeScoresA<Bi, Bo>(other) {

            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

        }

        ConvergeScoresWithMatchesA(ConvergeScoresWithMatchesA&& other) noexcept : // move constructor
                ConvergeScoresA<Bi, Bo>(std::move(other)), weighted_match_scores_(other.weighted_match_scores_) {

            other.weighted_match_scores_ = nullptr;

        }

        ConvergeScoresWithMatchesA& operator=(const ConvergeScoresWithMatchesA& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            ConvergeScoresA<Bi, Bo>::operator=(other);

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

            return *this;

        }

        ConvergeScoresWithMatchesA& operator=(ConvergeScoresWithMatchesA&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            ConvergeScoresA<Bi, Bo>::operator=(std::move(other));

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = other.weighted_match_scores_;
            other.weighted_match_scores_ = nullptr;

            return *this;

        }


        inline float match(char qual_a, char qual_b) override {

            float prob_a = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float prob_b = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_b));

            return max_change_match_ * bOuter_.boost((alphabet_size_ - 2 + prob_a + prob_b - prob_a * prob_b * alphabet_size_) / (alphabet_size_ - 1));

        }

        inline float lookup_match(char qual_a, char qual_b) override {
            return weighted_match_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

    protected:
        float* weighted_match_scores_; // matrix of weighted match scores as 1-dimensional C-style array

    };


    // Version B: limited convergence per operation

    /*
     * Provides quality-weighted scores using the following mechanism:
     * Matches receive an increasingly large penalty with decreasing quality, while the penalty of mismatches and gaps
     * is increasingly reduced with decreasing quality. This "convergence" is limited by an per-operation
     * maximum modification of the costs.
     *
     * Affects all edit operations (modifies the penalties).
     *
     * The costs of matches are not quality-weighted (see ConvergeScoresWithMatchesB for a version also weighting matches).
     */
    template<typename Bi, typename Bo>
    class ConvergeScoresB {

    public:
        ConvergeScoresB(const QualityEncoding<>& qe, const Configuration& config) : bInner_(config), bOuter_(config), qe_(qe) {

            ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
            mismatch_ = sf.penalty_mismatch;
            gap_open_ = sf.penalty_open;
            gap_extend_ = sf.penalty_extend;

            alphabet_size_ = (config.alphabet.empty()) ? 4 : config.alphabet.length();

            // determine how much each operation can change at most due to quality-weighting
            max_change_match_ = mismatch_ / 2.0f;
            max_change_mismatch_ = mismatch_ / 2.0f;
            max_change_gap_open_ = gap_open_ / 2.0f;
            max_change_gap_extend_ = gap_extend_ / 2.0f;
            if (config.misc.find("max_change_match") != config.misc.end()) max_change_match_ = std::stof(config.misc.find("max_change_match")->second);
            if (config.misc.find("max_change_mismatch") != config.misc.end()) max_change_mismatch_ = std::stof(config.misc.find("max_change_mismatch")->second);
            if (config.misc.find("max_change_gap_open") != config.misc.end()) max_change_gap_open_ = std::stof(config.misc.find("max_change_gap_open")->second);
            if (config.misc.find("max_change_gap_extend") != config.misc.end()) max_change_gap_extend_ = std::stof(config.misc.find("max_change_gap_extend")->second);

            // precompute the weighted mismatch, gap-opening and gap-extension scores
            auto chars = qe_.get_accepted_scores();
            num_scores_ = chars.size();
            score_shift_ = qe_.get_offset() + qe_.get_from();
            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            weighted_gap_open_scores_ = new float[num_scores_];
            weighted_gap_extend_scores_ = new float[num_scores_];
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_mismatch_scores_[index] = mismatch(a, b);

                }

                weighted_gap_open_scores_[a - score_shift_] = gap_open(a);
                weighted_gap_extend_scores_[a - score_shift_] = gap_extend(a);

            }

        }

        virtual ~ConvergeScoresB() {

            delete[] weighted_gap_extend_scores_;
            delete[] weighted_gap_open_scores_;
            delete[] weighted_mismatch_scores_;

        }

        ConvergeScoresB(const ConvergeScoresB& other) : // copy constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                max_change_match_(other.max_change_match_), max_change_mismatch_(other.max_change_mismatch_),
                max_change_gap_open_(other.max_change_gap_open_), max_change_gap_extend_(other.max_change_gap_extend_),
                mismatch_(other.mismatch_), gap_open_(other.gap_open_), gap_extend_(other.gap_extend_),
                score_shift_(other.score_shift_), num_scores_(other.num_scores_) {

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

            weighted_gap_open_scores_ = new float[num_scores_];
            weighted_gap_extend_scores_ = new float[num_scores_];
            for (unsigned long i = 0; i < num_scores_; i++) {

                weighted_gap_open_scores_[i] = other.weighted_gap_open_scores_[i];
                weighted_gap_extend_scores_[i] = other.weighted_gap_extend_scores_[i];

            }

        }

        ConvergeScoresB(ConvergeScoresB&& other) noexcept : // move constructor
                bInner_(other.bInner_), bOuter_(other.bOuter_), qe_(other.qe_), alphabet_size_(other.alphabet_size_),
                max_change_match_(other.max_change_match_), max_change_mismatch_(other.max_change_mismatch_),
                max_change_gap_open_(other.max_change_gap_open_), max_change_gap_extend_(other.max_change_gap_extend_),
                mismatch_(other.mismatch_), gap_open_(other.gap_open_), gap_extend_(other.gap_extend_),
                weighted_mismatch_scores_(other.weighted_mismatch_scores_),
                weighted_gap_open_scores_(other.weighted_gap_open_scores_),
                weighted_gap_extend_scores_(other.weighted_gap_extend_scores_),
                score_shift_(other.score_shift_), num_scores_(other.num_scores_) {

            other.weighted_mismatch_scores_ = nullptr;
            other.weighted_gap_open_scores_ = nullptr;
            other.weighted_gap_extend_scores_ = nullptr;

        }

        ConvergeScoresB& operator=(const ConvergeScoresB& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;
            delete[] weighted_gap_open_scores_;
            delete[] weighted_gap_extend_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            max_change_match_ = other.max_change_match_;
            max_change_mismatch_ = other.max_change_mismatch_;
            max_change_gap_open_ = other.max_change_gap_open_;
            max_change_gap_extend_ = other.max_change_gap_extend_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_mismatch_scores_[i] = other.weighted_mismatch_scores_[i];
            }

            weighted_gap_open_scores_ = new float[num_scores_];
            weighted_gap_extend_scores_ = new float[num_scores_];
            for (unsigned long i = 0; i < num_scores_; i++) {

                weighted_gap_open_scores_[i] = other.weighted_gap_open_scores_[i];
                weighted_gap_extend_scores_[i] = other.weighted_gap_extend_scores_[i];

            }

            return *this;

        }

        ConvergeScoresB& operator=(ConvergeScoresB&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            // release old resources
            delete[] weighted_mismatch_scores_;
            delete[] weighted_gap_open_scores_;
            delete[] weighted_gap_extend_scores_;

            // copy new resources
            bInner_ = other.bInner_;
            bOuter_ = other.bOuter_;
            qe_ = other.qe_;
            alphabet_size_ = other.alphabet_size_;
            max_change_match_ = other.max_change_match_;
            max_change_mismatch_ = other.max_change_mismatch_;
            max_change_gap_open_ = other.max_change_gap_open_;
            max_change_gap_extend_ = other.max_change_gap_extend_;
            mismatch_ = other.mismatch_;
            gap_open_ = other.gap_open_;
            gap_extend_ = other.gap_extend_;
            score_shift_ = other.score_shift_;
            num_scores_ = other.num_scores_;

            weighted_mismatch_scores_ = other.weighted_mismatch_scores_;
            other.weighted_mismatch_scores_ = nullptr;
            weighted_gap_open_scores_ = other.weighted_gap_open_scores_;
            other.weighted_gap_open_scores_ = nullptr;
            weighted_gap_extend_scores_ = other.weighted_gap_extend_scores_;
            other.weighted_gap_extend_scores_ = nullptr;

            return *this;

        }


        virtual inline float match(char qual_a, char qual_b) {
            return 0;
        }

        inline float mismatch(char qual_a, char qual_b) {

            float prob_a = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float prob_b = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_b));

            return mismatch_ - max_change_mismatch_ * bOuter_.boost((alphabet_size_ - 2 + prob_a + prob_b - prob_a * prob_b * alphabet_size_) / ((alphabet_size_ - 1) * (alphabet_size_ - 1)));

        }

        inline float substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? match(qual_a, qual_b) : mismatch(qual_a, qual_b);
        }

        inline float gap_open(char qual_a) {
            return gap_open_ - max_change_gap_open_ * bOuter_.boost(bInner_.boost(qe_.encoded_quality_to_probability(qual_a)));
        }

        inline float gap_extend(char qual_a) {
            return gap_extend_ - max_change_gap_extend_ * bOuter_.boost(bInner_.boost(qe_.encoded_quality_to_probability(qual_a)));
        }


        virtual inline float lookup_match(char qual_a, char qual_b) {
            return 0;
        }

        inline float lookup_mismatch(char qual_a, char qual_b) {
            return weighted_mismatch_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

        inline float lookup_substitute(char a, char qual_a, char b, char qual_b) {
            return (a == b) ? lookup_match(qual_a, qual_b) : lookup_mismatch(qual_a, qual_b);
        }

        inline float lookup_gap_open(char qual_a) {
            return weighted_gap_open_scores_[qual_a - score_shift_];
        }

        inline float lookup_gap_extend(char qual_a) {
            return weighted_gap_extend_scores_[qual_a - score_shift_];
        }

    protected:
        Bi bInner_; // inner boosting function (for individual error probabilities)
        Bo bOuter_; // outer boosting function (for combined error probabilities)
        QualityEncoding<> qe_; // encoding of the quality score
        lenSeqs_t alphabet_size_; // size of the used sequence alphabet
        float max_change_match_; // maximum change of match cost due to quality-weighting
        float max_change_mismatch_; // maximum change of mismatch cost due to quality-weighting
        float max_change_gap_open_; // maximum change of gap-opening cost due to quality-weighting
        float max_change_gap_extend_; // maximum change of gap-extension cost due to quality-weighting

        float mismatch_; // mismatch penalty
        float gap_open_; // gap-opening penalty
        float gap_extend_; // gap-extension penalty

        float* weighted_mismatch_scores_; // matrix of weighted mismatch scores as 1-dimensional C-style array
        float* weighted_gap_open_scores_; // matrix of weighted gap-opening scores as 1-dimensional C-style array
        float* weighted_gap_extend_scores_; // matrix of weighted gap-extension scores as 1-dimensional C-style array
        char score_shift_; // (offset + from) of the quality encoding to normalise the quality scores
        unsigned long num_scores_; // number of levels in the quality encoding

    };

    template<typename Bi, typename Bo>
    class ConvergeScoresWithMatchesB : public ConvergeScoresB<Bi, Bo> {

        using ConvergeScoresB<Bi, Bo>::bInner_;
        using ConvergeScoresB<Bi, Bo>::bOuter_;
        using ConvergeScoresB<Bi, Bo>::qe_;
        using ConvergeScoresB<Bi, Bo>::alphabet_size_;
        using ConvergeScoresB<Bi, Bo>::max_change_match_;
        using ConvergeScoresB<Bi, Bo>::score_shift_;
        using ConvergeScoresB<Bi, Bo>::num_scores_;

        using ConvergeScoresB<Bi, Bo>::mismatch;
        using ConvergeScoresB<Bi, Bo>::lookup_mismatch;

    public:
        ConvergeScoresWithMatchesB(const QualityEncoding<>& qe, const Configuration& config) : ConvergeScoresB<Bi, Bo>(qe, config) {

            // also precompute the weighted match scores
            auto chars = qe_.get_accepted_scores();
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (auto a : chars) {
                for (auto b : chars) {

                    auto index = (a - score_shift_) * num_scores_ + (b - score_shift_);
                    weighted_match_scores_[index] = match(a, b);

                }
            }

        }

        ~ConvergeScoresWithMatchesB() {
            delete[] weighted_match_scores_;
        }

        ConvergeScoresWithMatchesB(const ConvergeScoresWithMatchesB& other) : // copy constructor
                ConvergeScoresB<Bi, Bo>(other) {

            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

        }

        ConvergeScoresWithMatchesB(ConvergeScoresWithMatchesB&& other) noexcept : // move constructor
                ConvergeScoresB<Bi, Bo>(std::move(other)), weighted_match_scores_(other.weighted_match_scores_) {

            other.weighted_match_scores_ = nullptr;

        }

        ConvergeScoresWithMatchesB& operator=(const ConvergeScoresWithMatchesB& other) { // copy assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            ConvergeScoresB<Bi, Bo>::operator=(other);

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = new float[num_scores_ * num_scores_];
            for (unsigned long i = 0; i < num_scores_ * num_scores_; i++) {
                weighted_match_scores_[i] = other.weighted_match_scores_[i];
            }

            return *this;

        }

        ConvergeScoresWithMatchesB& operator=(ConvergeScoresWithMatchesB&& other) noexcept { // move assignment operator

            // check for self-assignment
            if (&other == this) {
                return *this;
            }

            ConvergeScoresB<Bi, Bo>::operator=(std::move(other));

            // release old resources
            delete[] weighted_match_scores_;

            // copy new resources
            weighted_match_scores_ = other.weighted_match_scores_;
            other.weighted_match_scores_ = nullptr;

            return *this;

        }


        inline float match(char qual_a, char qual_b) override {

            float prob_a = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_a));
            float prob_b = 1.0f - bInner_.boost(qe_.encoded_quality_to_probability(qual_b));

            return max_change_match_ * bOuter_.boost((alphabet_size_ - 2 + prob_a + prob_b - prob_a * prob_b * alphabet_size_) / (alphabet_size_ - 1));

        }

        inline float lookup_match(char qual_a, char qual_b) override {
            return weighted_match_scores_[(qual_a - score_shift_) * num_scores_ + (qual_b - score_shift_)];
        }

    protected:
        float* weighted_match_scores_; // matrix of weighted match scores as 1-dimensional C-style array

    };



    /* ===== Alignment methods ===== */

    /*
     * Computes the bounded number of edit operations in a quality-weighted, optimal alignment
     * using the specified quality-weighting mechanism (template parameter S).
     *
     * Terminates early when the final number of edit operations for each entry in the current row
     * is guaranteed to exceed the threshold.
     */
    template<typename S>
    class BoundedQualityLevenshteinDistance : public Distance {

    public:
        BoundedQualityLevenshteinDistance(lenSeqs_t max_seq_len, lenSeqs_t threshold, const QualityEncoding<>& qe,
                const Configuration& config) : scores_(qe, config) {

            row_length_ = max_seq_len + 1;
            d_row_ = new float[row_length_];
            p_row_ = new float[row_length_];
            cnt_diffs_row_ = new lenSeqs_t[row_length_];
            cnt_diffs_p_row_ = new lenSeqs_t[row_length_];
            threshold_ = threshold;

        }

        virtual ~BoundedQualityLevenshteinDistance() {

            delete[] cnt_diffs_p_row_;
            delete[] cnt_diffs_row_;
            delete[] p_row_;
            delete[] d_row_;

        }

        BoundedQualityLevenshteinDistance* clone() const override { // deep-copy clone method
            return new BoundedQualityLevenshteinDistance(*this);
        }

        dist_t distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override {
            return compute_weighted_gotoh(ac.seq(i), ac.len(i), ac.quals(i), ac.seq(j), ac.len(j), ac.quals(j));
        }

    protected:
        BoundedQualityLevenshteinDistance(const BoundedQualityLevenshteinDistance& other) : // copy constructor
                threshold_(other.threshold_), row_length_(other.row_length_), scores_(other.scores_) {

            d_row_ = new dist_t[row_length_];
            p_row_ = new dist_t[row_length_];
            cnt_diffs_row_ = new lenSeqs_t[row_length_];
            cnt_diffs_p_row_ = new lenSeqs_t[row_length_];
            for (lenSeqs_t i = 0; i < row_length_; i++) {

                d_row_[i] = other.d_row_[i];
                p_row_[i] = other.p_row_[i];
                cnt_diffs_row_[i] = other.cnt_diffs_row_[i];
                cnt_diffs_p_row_[i] = other.cnt_diffs_p_row_[i];

            }

        }

        BoundedQualityLevenshteinDistance(BoundedQualityLevenshteinDistance&& other) = delete; // move constructor

        BoundedQualityLevenshteinDistance& operator=(const BoundedQualityLevenshteinDistance& other) = delete; // copy assignment operator

        BoundedQualityLevenshteinDistance& operator=(BoundedQualityLevenshteinDistance&& other) = delete; // move assignment operator


        virtual dist_t compute_weighted_gotoh(const char* s, const lenSeqs_t len_s, const char* quals_s,
                const char* t, const lenSeqs_t len_t, const char* quals_t) {

            // long computation not necessary if lengths differ too much
            if (((len_s > len_t) ? (len_s - len_t) : (len_t - len_s)) > threshold_) {
                return threshold_ + 1;
            }

            if (threshold_ == 0) {
                return lenSeqs_t(s != t);
            }

            const char* shorter = (len_s < len_t) ? s : t;
            const char* quals_shorter = (len_s < len_t) ? quals_s : quals_t;
            lenSeqs_t len_shorter = std::min(len_s, len_t);
            const char* longer = (len_s >= len_t) ? s : t;
            const char* quals_longer = (len_s >= len_t) ? quals_s : quals_t;
            lenSeqs_t len_longer = std::max(len_s, len_t);
            lenSeqs_t delta = len_longer - len_shorter;

            // (mis)match is only possibility when we have to consider only one diagonal
            // [happens only if (a) threshold = delta = 0, or (b) threshold = 1 and delta = 0, but (a) is already covered above]
            if ((threshold_ - delta) / 2 == 0 && (threshold_ + delta) / 2 == 0) {

                lenSeqs_t diffs = 0;
                for (auto i = 0; diffs <= threshold_ && i < len_shorter; i++) {
                    diffs += (shorter[i] != longer[i]);
                }

                return diffs;

            }

            // initialise necessary section of first row
            d_row_[0] = scores_.lookup_gap_open(quals_longer[0]);
            for (lenSeqs_t j = 1; j <= (threshold_ + delta) / 2 && j <= len_longer; j++) {

                d_row_[j] = d_row_[j - 1] + scores_.lookup_gap_extend(quals_longer[j - 1]);
                p_row_[j] = POS_INF_FLOAT;
                cnt_diffs_row_[j] = j;

            }
            d_row_[0] = 0;


            // the helper variables are used as in BoundedScoreLevenshteinDistance
            dist_t match, min_val, val_q, from_d, from_pq, tmp_gap;
            lenSeqs_t j, diff, min_val_diff, diffs_q = 0;
            bool early;

            // compute remaining rows
            tmp_gap = scores_.lookup_gap_open(quals_shorter[0]);
            for (lenSeqs_t i = 1; i <= len_shorter; i++) {

                // handle left end
                tmp_gap += scores_.lookup_gap_extend(quals_shorter[i - 1]);
                d_row_[0] = tmp_gap;
                match = (i > (threshold_ - delta) / 2 + 1) ? d_row_[i - (threshold_ - delta) / 2 - 1] :
                        (((1 < i) && (i <= (threshold_ - delta) / 2 + 1)) * (d_row_[0] - scores_.lookup_gap_extend(quals_shorter[i - 1])));
                val_q = POS_INF_FLOAT;
                diff = (i > (threshold_ - delta) / 2 + 1) ? cnt_diffs_row_[i - (threshold_ - delta) / 2 - 1] : (i - 1);
                cnt_diffs_row_[0] = i;

                early = true; // early termination flag

                // fill remaining row
                j = 1 + (i > (threshold_ - delta) / 2) * (i - (threshold_ - delta) / 2 - 1); // same as starting from j = max(1, i - (threshold - delta) / 2) with signed integers
                d_row_[j - 1] = POS_INF_FLOAT;
                if (i + (threshold_ + delta) / 2 <= len_longer) {
                    d_row_[i + (threshold_ + delta) / 2] = p_row_[i + (threshold_ + delta) / 2] = POS_INF_FLOAT;
                }

                for (; j <= i + (threshold_ + delta) / 2 && j <= len_longer; j++) {

                    // arrays P & cntDiffsP
                    from_d = d_row_[j] + scores_.lookup_gap_open(quals_shorter[i - 1]) + scores_.lookup_gap_extend(quals_shorter[i - 1]);
                    from_pq = p_row_[j] + scores_.lookup_gap_extend(quals_shorter[i - 1]);

                    if (from_d <= from_pq) {

                        p_row_[j] = from_d;
                        cnt_diffs_p_row_[j] = cnt_diffs_row_[j] + 1;

                    } else {

                        p_row_[j] = from_pq;
                        cnt_diffs_p_row_[j]++;

                    }

                    // arrays Q & cntDiffsQ
                    from_d = d_row_[j - 1] + scores_.lookup_gap_open(quals_longer[j - 1]) + scores_.lookup_gap_extend(quals_longer[j - 1]);
                    from_pq = val_q + scores_.lookup_gap_extend(quals_longer[j - 1]);
                    if (from_d <= from_pq) {

                        val_q = from_d;
                        diffs_q = cnt_diffs_row_[j - 1] + 1;

                    } else {

                        val_q = from_pq;
                        diffs_q++;

                    }

                    // arrays D & cntDiffs
                    min_val = match + scores_.lookup_substitute(shorter[i - 1], quals_shorter[i - 1], longer[j - 1], quals_longer[j - 1]);
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

                if (early) {// computation can be terminated early if computed row contains only values > threshold (because values are monotonically increasing)
                    return threshold_ + 1;
                }

            }

            return (cnt_diffs_row_[len_longer] > threshold_) ? (threshold_ + 1) : cnt_diffs_row_[len_longer];

        }


        float* d_row_; // one row of matrix D (cost of alignment of prefixes)
        float* p_row_; // one row of matrix P (cost of alignment of prefixes that ends with gap in t)
        lenSeqs_t* cnt_diffs_row_; // number of differences corresponding to alignments based on D
        lenSeqs_t* cnt_diffs_p_row_; // number of differences corresponding to alignments based on P
        lenSeqs_t threshold_; // distance threshold
        lenSeqs_t row_length_; // length of matrix row

        S scores_; // (quality-weighted) alignment score cost function

    };

    /*
     * Computes the bounded score of a quality-weighted, optimal alignment using the specified
     * quality-weighting mechanism (template parameter S).
     *
     * Terminates early when all scores in the current row are already larger than the threshold.
     */
    template<typename S>
    class BoundedQualityAlignmentScore : public Distance {

    public:
        BoundedQualityAlignmentScore(lenSeqs_t max_seq_len, lenSeqs_t threshold, const QualityEncoding<>& qe,
                const Configuration& config) : scores_(qe, config) {

            row_length_ = max_seq_len + 1;
            d_row_ = new float[row_length_];
            p_row_ = new float[row_length_];
            threshold_ = threshold;

        }

        virtual ~BoundedQualityAlignmentScore() {

            delete[] p_row_;
            delete[] d_row_;

        }

        BoundedQualityAlignmentScore* clone() const override { // deep-copy clone method
            return new BoundedQualityAlignmentScore(*this);
        }

        dist_t distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override {
            return compute_weighted_gotoh(ac.seq(i), ac.len(i), ac.quals(i), ac.seq(j), ac.len(j), ac.quals(j));
        }

    protected:
        BoundedQualityAlignmentScore(const BoundedQualityAlignmentScore& other) : // copy constructor
                threshold_(other.threshold_), row_length_(other.row_length_), scores_(other.scores_) {

            d_row_ = new dist_t[row_length_];
            p_row_ = new dist_t[row_length_];
            for (lenSeqs_t i = 0; i < row_length_; i++) {

                d_row_[i] = other.d_row_[i];
                p_row_[i] = other.p_row_[i];

            }

        }

        BoundedQualityAlignmentScore(BoundedQualityAlignmentScore&& other) = delete; // move constructor

        BoundedQualityAlignmentScore& operator=(const BoundedQualityAlignmentScore& other) = delete; // copy assignment operator

        BoundedQualityAlignmentScore& operator=(BoundedQualityAlignmentScore&& other) = delete; // move assignment operator


        virtual dist_t compute_weighted_gotoh(const char* s, const lenSeqs_t len_s, const char* quals_s,
                const char* t, const lenSeqs_t len_t, const char* quals_t) {

            const char* shorter = (len_s < len_t) ? s : t;
            const char* quals_shorter = (len_s < len_t) ? quals_s : quals_t;
            lenSeqs_t len_shorter = std::min(len_s, len_t);
            const char* longer = (len_s >= len_t) ? s : t;
            const char* quals_longer = (len_s >= len_t) ? quals_s : quals_t;
            lenSeqs_t len_longer = std::max(len_s, len_t);

            // initialise necessary section of first row
            d_row_[0] = scores_.lookup_gap_open(quals_longer[0]);
            for (lenSeqs_t j = 1; j <= len_longer; j++) {

                d_row_[j] = d_row_[j - 1] + scores_.lookup_gap_extend(quals_longer[j - 1]);
                p_row_[j] = POS_INF_FLOAT;

            }
            d_row_[0] = 0;

            // the helper variables are used as in BoundedScoreDistance
            dist_t match, val_q = POS_INF_FLOAT, tmp_gap;
            bool early;

            // compute remaining rows
            tmp_gap = scores_.lookup_gap_open(quals_shorter[0]);
            for (lenSeqs_t i = 1; i <= len_shorter; i++) {

                // handle left end
                match = d_row_[0];
                tmp_gap += scores_.lookup_gap_extend(quals_shorter[i - 1]);
                d_row_[0] = tmp_gap;
                val_q = POS_INF_FLOAT;

                early = true; // early termination flag

                // fill remaining row
                for (lenSeqs_t j = 1; j <= len_longer; j++) {

                    // array P
                    p_row_[j] = std::min(d_row_[j] + scores_.lookup_gap_open(quals_shorter[i - 1]) + scores_.lookup_gap_extend(quals_shorter[i - 1]),
                                         p_row_[j] + scores_.lookup_gap_extend(quals_shorter[i - 1]));

                    // array Q
                    val_q = std::min(d_row_[j - 1] + scores_.lookup_gap_open(quals_longer[j - 1]) + scores_.lookup_gap_extend(quals_longer[j - 1]),
                                     val_q + scores_.lookup_gap_extend(quals_longer[j - 1]));

                    // array D
                    dist_t tmp = match + scores_.lookup_substitute(shorter[i - 1], quals_shorter[i - 1], longer[j - 1], quals_longer[j - 1]);
                    match = d_row_[j];
                    d_row_[j] = std::min({p_row_[j], val_q, tmp});

                    early &= (d_row_[j] > threshold_);

                }

                if (early) {// computation can be terminated early if computed row contains only values > threshold (because values are monotonically increasing)
                    return threshold_ + 1;
                }

            }

            return std::min({d_row_[len_longer], p_row_[len_longer], val_q});

        }


        float* d_row_; // one row of matrix D (cost of alignment of prefixes)
        float* p_row_; // one row of matrix P (cost of alignment of prefixes that ends with gap in t)
        lenSeqs_t threshold_; // distance threshold
        lenSeqs_t row_length_; // length of matrix row

        S scores_; // (quality-weighted) alignment score cost function

    };

    /*
     * Computes the bounded score of a quality-weighted, optimal alignment using the specified
     * quality-weighting mechanism (template parameter S).
     * The dynamic-programming matrix is limited to a specified number of bands on each side
     * of the main diagonal.
     *
     * Terminates early when all scores in the current row are already larger than the threshold.
     */
    template<typename S>
    class BoundedBandedQualityAlignmentScore : public Distance {

    public:
        BoundedBandedQualityAlignmentScore(lenSeqs_t max_seq_len, lenSeqs_t threshold, lenSeqs_t b,
                const QualityEncoding<>& qe, const Configuration& config) : scores_(qe, config) {

            row_length_ = max_seq_len + 1;
            d_row_ = new float[row_length_];
            p_row_ = new float[row_length_];
            threshold_ = threshold;
            bands_per_side_ = b;

        }

        virtual ~BoundedBandedQualityAlignmentScore() {

            delete[] p_row_;
            delete[] d_row_;

        }

        BoundedBandedQualityAlignmentScore* clone() const override { // deep-copy clone method
            return new BoundedBandedQualityAlignmentScore(*this);
        }

        dist_t distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) override {
            return compute_weighted_gotoh(ac.seq(i), ac.len(i), ac.quals(i), ac.seq(j), ac.len(j), ac.quals(j));
        }

    protected:
        BoundedBandedQualityAlignmentScore(const BoundedBandedQualityAlignmentScore& other) : // copy constructor
                threshold_(other.threshold_), bands_per_side_(other.bands_per_side_), row_length_(other.row_length_), scores_(other.scores_) {

            d_row_ = new dist_t[row_length_];
            p_row_ = new dist_t[row_length_];
            for (lenSeqs_t i = 0; i < row_length_; i++) {

                d_row_[i] = other.d_row_[i];
                p_row_[i] = other.p_row_[i];

            }

        }

        BoundedBandedQualityAlignmentScore(BoundedBandedQualityAlignmentScore&& other) = delete; // move constructor

        BoundedBandedQualityAlignmentScore& operator=(const BoundedBandedQualityAlignmentScore& other) = delete; // copy assignment operator

        BoundedBandedQualityAlignmentScore& operator=(BoundedBandedQualityAlignmentScore&& other) = delete; // move assignment operator


        virtual dist_t compute_weighted_gotoh(const char* s, const lenSeqs_t len_s, const char* quals_s,
                const char* t, const lenSeqs_t len_t, const char* quals_t) {

            const char* shorter = (len_s < len_t) ? s : t;
            const char* quals_shorter = (len_s < len_t) ? quals_s : quals_t;
            lenSeqs_t len_shorter = std::min(len_s, len_t);
            const char* longer = (len_s >= len_t) ? s : t;
            const char* quals_longer = (len_s >= len_t) ? quals_s : quals_t;
            lenSeqs_t len_longer = std::max(len_s, len_t);

            // initialise necessary section of first row
            d_row_[0] = scores_.lookup_gap_open(quals_longer[0]);
            for (lenSeqs_t j = 1; j <= bands_per_side_ && j <= len_longer; j++) {

                d_row_[j] = d_row_[j - 1] + scores_.lookup_gap_extend(quals_longer[j - 1]);
                p_row_[j] = POS_INF_FLOAT;

            }
            d_row_[0] = 0;

            // the helper variables are used as in BoundedBandedScoreDistance
            dist_t match, val_q = POS_INF_FLOAT, tmp_gap;
            lenSeqs_t j;
            bool early;

            // compute remaining rows
            tmp_gap = scores_.lookup_gap_open(quals_shorter[0]);
            for (lenSeqs_t i = 1; i <= len_shorter; i++) {

                // handle left end
                match = (i > bands_per_side_ + 1) ? d_row_[i - bands_per_side_ - 1] : (((1 < i) && (i <= bands_per_side_ + 1)) * tmp_gap);
                tmp_gap += scores_.lookup_gap_extend(quals_shorter[i - 1]);
                d_row_[0] = tmp_gap;
                val_q = POS_INF_FLOAT;

                early = true; // early termination flag

                // fill remaining row
                j = 1 + (i > bands_per_side_) * (i - bands_per_side_ - 1); // same as starting from j = max(1, i - bound) with signed integers
                d_row_[j - 1] = POS_INF_FLOAT;
                if (i + bands_per_side_ <= len_longer) {
                    d_row_[i + bands_per_side_] = p_row_[i + bands_per_side_] = POS_INF_FLOAT;
                }
                for (; j <= i + bands_per_side_ && j <= len_longer; j++) {

                    // array P
                    p_row_[j] = std::min(d_row_[j] + scores_.lookup_gap_open(quals_shorter[i - 1]) + scores_.lookup_gap_extend(quals_shorter[i - 1]),
                                         p_row_[j] + scores_.lookup_gap_extend(quals_shorter[i - 1]));

                    // array Q
                    val_q = std::min(d_row_[j - 1] + scores_.lookup_gap_open(quals_longer[j - 1]) + scores_.lookup_gap_extend(quals_longer[j - 1]),
                                     val_q + scores_.lookup_gap_extend(quals_longer[j - 1]));

                    // array D
                    dist_t tmp = match + scores_.lookup_substitute(shorter[i - 1], quals_shorter[i - 1], longer[j - 1], quals_longer[j - 1]);
                    match = d_row_[j];
                    d_row_[j] = std::min({p_row_[j], val_q, tmp});

                    early &= (d_row_[j] > threshold_);

                }

                if (early) {// computation can be terminated early if computed row contains only values > threshold (because values are monotonically increasing)
                    return threshold_ + 1;
                }

            }

            return std::min({d_row_[len_longer], p_row_[len_longer], val_q});

        }


        float* d_row_; // one row of matrix D (cost of alignment of prefixes)
        float* p_row_; // one row of matrix P (cost of alignment of prefixes that ends with gap in t)
        lenSeqs_t threshold_; // distance threshold
        lenSeqs_t row_length_; // length of matrix row
        lenSeqs_t bands_per_side_; // number of diagonals (bands) per side of the main diagonal considered

        S scores_; // (quality-weighted) alignment score cost function

    };


    /*
     * Generator methods for constructing an alignment method in the demanded setup:
     *  - template parameter A = alignment method (e.g. BoundedQualityLevenshteinDistance)
     *  - template parameter W = quality-weighting mechanism (with / without matches; e.g. FrithScores)
     *  - template parameter C = configuration type (e.g. QualityAlignmentScoreConfiguration)
     *
     * Equips the aligner and quality-weighting with the "type" and "position" of the boosting method:
     * "type" = e.g. multiplicative boosting, "position" = inner or outer boosting
     */

    // no externally provided limitation of the bands of the dynamic-programming scheme
    template<template<class> class A, template<class, class> class W, typename C>
    Distance* get_qw_distance(lenSeqs_t max_seq_len, lenSeqs_t t, const C &config) {

        if (config.inner_boost) {
            switch (config.opt_boosting) {
                case BO_NO_EFFECT: return new A<W<NoEffectBooster, NoEffectBooster>>(max_seq_len, t, config.qe_map, config);
                case BO_MULT: return new A<W<MultBooster, NoEffectBooster>>(max_seq_len, t, config.qe_map, config);
                case BO_ROOT: return new A<W<RootBooster, NoEffectBooster>>(max_seq_len, t, config.qe_map, config);
                case BO_LINEAR: return new A<W<LinearBooster, NoEffectBooster>>(max_seq_len, t, config.qe_map, config);
                default: return nullptr;
            }
        } else {
            switch (config.opt_boosting) {
                case BO_NO_EFFECT: return new A<W<NoEffectBooster, NoEffectBooster>>(max_seq_len, t, config.qe_map, config);
                case BO_MULT: return new A<W<NoEffectBooster, MultBooster>>(max_seq_len, t, config.qe_map, config);
                case BO_ROOT: return new A<W<NoEffectBooster, RootBooster>>(max_seq_len, t, config.qe_map, config);
                case BO_LINEAR: return new A<W<NoEffectBooster, LinearBooster>>(max_seq_len, t, config.qe_map, config);
                default: return nullptr;
            }
        }

    }

    // limits the dynamic-programming scheme to b bands on each side of the main diagonal
    template<template<class> class A, template<class, class> class W, typename C>
    Distance* get_bqw_distance(lenSeqs_t max_seq_len, lenSeqs_t t, long long b, const C &config) {

        if (config.inner_boost) {
            switch (config.opt_boosting) {
                case BO_NO_EFFECT: return new A<W<NoEffectBooster, NoEffectBooster>>(max_seq_len, t, b, config.qe_map, config);
                case BO_MULT: return new A<W<MultBooster, NoEffectBooster>>(max_seq_len, t, b, config.qe_map, config);
                case BO_ROOT: return new A<W<RootBooster, NoEffectBooster>>(max_seq_len, t, b, config.qe_map, config);
                case BO_LINEAR: return new A<W<LinearBooster, NoEffectBooster>>(max_seq_len, t, b, config.qe_map, config);
                default: return nullptr;
            }
        } else {
            switch (config.opt_boosting) {
                case BO_NO_EFFECT: return new A<W<NoEffectBooster, NoEffectBooster>>(max_seq_len, t, b, config.qe_map, config);
                case BO_MULT: return new A<W<NoEffectBooster, MultBooster>>(max_seq_len, t, b, config.qe_map, config);
                case BO_ROOT: return new A<W<NoEffectBooster, RootBooster>>(max_seq_len, t, b, config.qe_map, config);
                case BO_LINEAR: return new A<W<NoEffectBooster, LinearBooster>>(max_seq_len, t, b, config.qe_map, config);
                default: return nullptr;
            }
        }

    }

}

#endif //GEFAST_QUALITYWEIGHTEDDISTANCES_HPP
