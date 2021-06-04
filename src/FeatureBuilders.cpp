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

#include <iostream>

#include "../include/FeatureBuilders.hpp"

namespace GeFaST {

    /* === WordCompBuilder == */

    WordCompBuilder::WordCompBuilder(const std::string& alphabet, short k) :
            alphabet_size_(alphabet.size()), k_(k) {

        short r = 0;
        for (char a : alphabet) {

            nt_map_[a] = r;
            r += 1;

        }

        // number of features = number of different k-mers
        num_features_ = 1;
        for (auto i = 0; i < k_; i++) {
            num_features_ *= alphabet_size_;
        }

    }

    WordCompBuilder* WordCompBuilder::clone() const {
        return new WordCompBuilder(*this);
    }

    size_t WordCompBuilder::num_features() {
        return num_features_;
    }

    std::vector<feat_t> WordCompBuilder::get_features(const std::string& seq) {

        std::vector<feat_t> features(num_features_);
        compute_features(seq, features.data());
        return features;

    }

    void WordCompBuilder::get_features(const std::string& seq, feat_t* features) {
        compute_features(seq, features);
    }

    void WordCompBuilder::get_features(const std::string& seq, std::vector<feat_t>& features) {
        compute_features(seq, features.data());
    }

    void WordCompBuilder::compute_features(const std::string& seq, feat_t* counts) {

        std::fill(counts, counts + num_features_, 0);

        // count the number of occurrences of the different k-mers
        size_t kmer = 0;
        lenSeqs_t i = 0;
        for (;(i < k_ - 1) && (i < seq.length()); i++) {
            kmer = kmer * alphabet_size_ + nt_map_[seq[i]];
        }
        for (; i < seq.length(); i++) {
            kmer = (kmer * alphabet_size_ + nt_map_[seq[i]]) % num_features_;
            counts[kmer] += 1;

        }

    }


    /* === CpfBuilder == */

    std::unordered_map<std::string, short> CpfBuilder::ry_map =
            {{"AA", 0}, {"AC", 1}, {"AG", 0}, {"AT", 1}, {"CA", 2}, {"CC", 3}, {"CG", 2}, {"CT", 3},
             {"GA", 0}, {"GC", 1}, {"GG", 0}, {"GT", 1}, {"TA", 2}, {"TC", 3}, {"TG", 2}, {"TT", 3},
             {"AU", 1}, {"CU", 3}, {"GU", 1}, {"UA", 2}, {"UC", 3}, {"UG", 2}, {"UU", 3}};
    std::unordered_map<std::string, short> CpfBuilder::mk_map =
            {{"AA", 0}, {"AC", 0}, {"AG", 1}, {"AT", 1}, {"CA", 0}, {"CC", 0}, {"CG", 1}, {"CT", 1},
             {"GA", 2}, {"GC", 2}, {"GG", 3}, {"GT", 3}, {"TA", 2}, {"TC", 2}, {"TG", 3}, {"TT", 3},
             {"AU", 1}, {"CU", 1}, {"GU", 3}, {"UA", 2}, {"UC", 2}, {"UG", 3}, {"UU", 3}};
    std::unordered_map<std::string, short> CpfBuilder::ws_map =
            {{"AA", 0}, {"AC", 1}, {"AG", 1}, {"AT", 0}, {"CA", 2}, {"CC", 3}, {"CG", 3}, {"CT", 2},
             {"GA", 2}, {"GC", 3}, {"GG", 3}, {"GT", 2}, {"TA", 0}, {"TC", 1}, {"TG", 1}, {"TT", 0},
             {"AU", 0}, {"CU", 2}, {"GU", 2}, {"UA", 0}, {"UC", 1}, {"UG", 1}, {"UU", 0}};

    CpfBuilder::CpfBuilder(lenSeqs_t max_len) : max_len_(max_len) {

        psums_ry_ = new feat_t[4 * max_len_];
        psums_mk_ = new feat_t[4 * max_len_];
        psums_ws_ = new feat_t[4 * max_len_];
        lasts_ry_ = new feat_t[4];
        lasts_mk_ = new feat_t[4];
        lasts_ws_ = new feat_t[4];

    }

    CpfBuilder::~CpfBuilder() {

        delete[] lasts_ws_;
        delete[] lasts_mk_;
        delete[] lasts_ry_;
        delete[] psums_ws_;
        delete[] psums_mk_;
        delete[] psums_ry_;

    }

    CpfBuilder::CpfBuilder(const CpfBuilder& other) { // copy constructor

        num_features_ = other.num_features_;
        max_len_ = other.max_len_;

        psums_ry_ = new feat_t[4 * max_len_];
        memcpy(psums_ry_, other.psums_ry_, sizeof(feat_t) * 4 * max_len_);
        psums_mk_ = new feat_t[4 * max_len_];
        memcpy(psums_mk_, other.psums_mk_, sizeof(feat_t) * 4 * max_len_);
        psums_ws_ = new feat_t[4 * max_len_];
        memcpy(psums_ws_, other.psums_ws_, sizeof(feat_t) * 4 * max_len_);
        lasts_ry_ = new feat_t[4];
        memcpy(lasts_ry_, other.lasts_ry_, sizeof(feat_t) * 4);
        lasts_mk_ = new feat_t[4];
        memcpy(lasts_mk_, other.lasts_mk_, sizeof(feat_t) * 4);
        lasts_ws_ = new feat_t[4];
        memcpy(lasts_ws_, other.lasts_ws_, sizeof(feat_t) * 4);

    }

    CpfBuilder::CpfBuilder(CpfBuilder&& other) noexcept { // move constructor

        num_features_ = other.num_features_;
        max_len_ = other.max_len_;

        psums_ry_ = other.psums_ry_; other.psums_ry_ = nullptr;
        psums_mk_ = other.psums_mk_; other.psums_mk_ = nullptr;
        psums_ws_ = other.psums_ws_; other.psums_ws_ = nullptr;
        lasts_ry_ = other.lasts_ry_; other.lasts_ry_ = nullptr;
        lasts_mk_ = other.lasts_mk_; other.lasts_mk_ = nullptr;
        lasts_ws_ = other.lasts_ws_; other.lasts_ws_ = nullptr;

    }

    CpfBuilder& CpfBuilder::operator=(const CpfBuilder& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] lasts_ws_;
        delete[] lasts_mk_;
        delete[] lasts_ry_;
        delete[] psums_ws_;
        delete[] psums_mk_;
        delete[] psums_ry_;

        // copy new resources
        num_features_ = other.num_features_;
        max_len_ = other.max_len_;

        psums_ry_ = new feat_t[4 * max_len_];
        memcpy(psums_ry_, other.psums_ry_, sizeof(feat_t) * 4 * max_len_);
        psums_mk_ = new feat_t[4 * max_len_];
        memcpy(psums_mk_, other.psums_mk_, sizeof(feat_t) * 4 * max_len_);
        psums_ws_ = new feat_t[4 * max_len_];
        memcpy(psums_ws_, other.psums_ws_, sizeof(feat_t) * 4 * max_len_);
        lasts_ry_ = new feat_t[4];
        memcpy(lasts_ry_, other.lasts_ry_, sizeof(feat_t) * 4);
        lasts_mk_ = new feat_t[4];
        memcpy(lasts_mk_, other.lasts_mk_, sizeof(feat_t) * 4);
        lasts_ws_ = new feat_t[4];
        memcpy(lasts_ws_, other.lasts_ws_, sizeof(feat_t) * 4);

        return *this;

    }

    CpfBuilder& CpfBuilder::operator=(CpfBuilder&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] lasts_ws_;
        delete[] lasts_mk_;
        delete[] lasts_ry_;
        delete[] psums_ws_;
        delete[] psums_mk_;
        delete[] psums_ry_;

        // copy new resources
        num_features_ = other.num_features_;
        max_len_ = other.max_len_;

        psums_ry_ = other.psums_ry_; other.psums_ry_ = nullptr;
        psums_mk_ = other.psums_mk_; other.psums_mk_ = nullptr;
        psums_ws_ = other.psums_ws_; other.psums_ws_ = nullptr;
        lasts_ry_ = other.lasts_ry_; other.lasts_ry_ = nullptr;
        lasts_mk_ = other.lasts_mk_; other.lasts_mk_ = nullptr;
        lasts_ws_ = other.lasts_ws_; other.lasts_ws_ = nullptr;

        return *this;

    }

    CpfBuilder* CpfBuilder::clone() const {
        return new CpfBuilder(*this);
    }

    size_t CpfBuilder::num_features() {
        return num_features_;
    }

    std::vector<feat_t> CpfBuilder::get_features(const std::string& seq) {

        std::vector<feat_t> features(num_features_);
        compute_features(seq, features.data());
        return features;

    }

    void CpfBuilder::get_features(const std::string& seq, feat_t* features) {
        compute_features(seq, features);
    }

    void CpfBuilder::get_features(const std::string& seq, std::vector<feat_t>& features) {
        compute_features(seq, features.data());
    }

    void CpfBuilder::compute_features(const std::string& seq, feat_t* entropies) {

        // compute partial sums of local frequencies
        std::fill(psums_ry_, psums_ry_ + 4 * seq.length(), 0);
        std::fill(psums_mk_, psums_mk_ + 4 * seq.length(), 0);
        std::fill(psums_ws_, psums_ws_ + 4 * seq.length(), 0);

        std::fill(lasts_ry_, lasts_ry_ + 4, -1);
        std::fill(lasts_mk_, lasts_mk_ + 4, -1);
        std::fill(lasts_ws_, lasts_ws_ + 4, -1);

        for (lenSeqs_t j = 0; j < seq.length() - 1; j++) {

            std::string word = seq.substr(j, 2);

            // The four LF sequences of each categorisation are concatenated
            // and accessed via different offsets (via mapped).
            // The indicator sequences are not constructed explicitly but used implicitly via mapped.

            short mapped = ry_map[word];
            psums_ry_[mapped * seq.length() + j] = static_cast<feat_t>(1.0) / (j - lasts_ry_[mapped]);
            lasts_ry_[mapped] = j;

            mapped = mk_map[word];
            psums_mk_[mapped * seq.length() + j] = static_cast<feat_t>(1.0) / (j - lasts_mk_[mapped]);
            lasts_mk_[mapped] = j;

            mapped = ws_map[word];
            psums_ws_[mapped * seq.length() + j] = static_cast<feat_t>(1.0) / (j - lasts_ws_[mapped]);
            lasts_ws_[mapped] = j;

        }

        for (lenSeqs_t i = 0; i < 4; i++) {
            for (lenSeqs_t s = i * seq.length(), j = 1; j < seq.length(); j++) {

                psums_ry_[s + j] += psums_ry_[s + j - 1];
                psums_mk_[s + j] += psums_mk_[s + j - 1];
                psums_ws_[s + j] += psums_ws_[s + j - 1];

            }
        }

        // turn partial sums into discrete probabilities
        for (lenSeqs_t i = 0; i < 4; i++) {

            feat_t total = 0;
            for (lenSeqs_t s = i * seq.length(), j = 0; j < seq.length(); j++) {
                total += psums_ry_[s + j];
            }
            if (total > 0) {
                for (lenSeqs_t s = i * seq.length(), j = 0; j < seq.length(); j++) {
                    psums_ry_[s + j] /= total;
                }
            }

            total = 0;
            for (lenSeqs_t s = i * seq.length(), j = 0; j < seq.length(); j++) {
                total += psums_mk_[s + j];
            }
            if (total > 0) {
                for (lenSeqs_t s = i * seq.length(), j = 0; j < seq.length(); j++) {
                    psums_mk_[s + j] /= total;
                }
            }

            total = 0;
            for (lenSeqs_t s = i * seq.length(), j = 0; j < seq.length(); j++) {
                total += psums_ws_[s + j];
            }
            if (total > 0) {
                for (lenSeqs_t s = i * seq.length(), j = 0; j < seq.length(); j++) {
                    psums_ws_[s + j] /= total;
                }
            }

        }

        // compute LF-based entropy values
        for (lenSeqs_t i = 0; i < 4; i++) {

            feat_t sum_ry = 0;
            feat_t sum_mk = 0;
            feat_t sum_ws = 0;
            for (lenSeqs_t s = i * seq.length(), j = 0; j < seq.length(); j++) {

                if (psums_ry_[s + j] > 0) sum_ry += psums_ry_[s + j] * log2(psums_ry_[s + j]);
                if (psums_mk_[s + j] > 0) sum_mk += psums_mk_[s + j] * log2(psums_mk_[s + j]);
                if (psums_ws_[s + j] > 0) sum_ws += psums_ws_[s + j] * log2(psums_ws_[s + j]);

            }

            entropies[i] = -sum_ry;
            entropies[4 + i] = -sum_mk;
            entropies[8 + i] = -sum_ws;

        }

    }


    /* === DetBuilder == */

    DetBuilder::DetBuilder(const std::string& alphabet, feat_t alpha, lenSeqs_t max_length) :
            num_features_(alphabet.size() * alphabet.size()), weights_(max_length), alpha_(alpha) {

        short r = 0;
        for (char a : alphabet) {
            for (char b : alphabet) {

                pair_map_[{a, b}] = r;
                r += 1;

            }
        }

        for (lenSeqs_t d = 1; d < max_length; d++) {
            weights_[d] = std::pow(d, -alpha_);
        }

    }

    DetBuilder::DetBuilder(const std::string& alphabet, lenSeqs_t pref_dist, lenSeqs_t max_length) :
            num_features_(alphabet.size() * alphabet.size()), weights_(max_length),
            alpha_(std::log(static_cast<feat_t>(10)) / std::log(static_cast<feat_t>(pref_dist))) {

        short r = 0;
        for (char a : alphabet) {
            for (char b : alphabet) {

                pair_map_[{a, b}] = r;
                r += 1;

            }
        }

        for (lenSeqs_t d = 1; d < max_length; d++) {
            weights_[d] = std::pow(d, -alpha_);
        }

    }

    DetBuilder* DetBuilder::clone() const {
        return new DetBuilder(*this);
    }

    size_t DetBuilder::num_features() {
        return num_features_;
    }

    std::vector<feat_t> DetBuilder::get_features(const std::string& seq) {

        std::vector<feat_t> arc_weights(num_features_);
        compute_features(seq, arc_weights.data());
        return arc_weights;

    }

    void DetBuilder::get_features(const std::string& seq, feat_t* features) {
        compute_features(seq, features);
    }

    void DetBuilder::get_features(const std::string& seq, std::vector<feat_t>& features) {
        compute_features(seq, features.data());
    }

    void DetBuilder::compute_features(const std::string& seq, feat_t* arc_weights) {

        std::fill(arc_weights, arc_weights + num_features_, 0);

        // compute aggregated arc weights in simplified graph
        for (lenSeqs_t d = 1; d < seq.length(); d++) {

            feat_t weight = weights_[d];

            for (lenSeqs_t i = 0; i < seq.length() - d; i++) {
                arc_weights[pair_map_[{seq[i], seq[i + d]}]] += weight;
            }

        }

    }


    /* === BbcBuilder == */

    BbcBuilder::BbcBuilder(const std::string& alphabet, short max_dist) :
            alphabet_size_(alphabet.size()), max_dist_(max_dist) {

        short r = 0;
        for (char a : alphabet) {

            nt_map_[a] = r;
            r += 1;

        }

        num_features_ = alphabet_size_ * alphabet_size_;

        single_probs_ = new feat_t[alphabet_size_];
        single_cnts_ = new feat_t[alphabet_size_];
        ind_probs_ = new feat_t[num_features_];
        joint_probs_ = new feat_t[num_features_];

    }

    BbcBuilder::~BbcBuilder() {

        delete[] joint_probs_;
        delete[] ind_probs_;
        delete[] single_cnts_;
        delete[] single_probs_;

    }

    BbcBuilder::BbcBuilder(const BbcBuilder& other) { // copy constructor

        num_features_ = other.num_features_;
        nt_map_ = other.nt_map_;
        max_dist_ = other.max_dist_;
        alphabet_size_ = other.alphabet_size_;

        single_probs_ = new feat_t[alphabet_size_];
        memcpy(single_probs_, other.single_probs_, sizeof(feat_t) * alphabet_size_);
        single_cnts_ = new feat_t[alphabet_size_];
        memcpy(single_cnts_, other.single_cnts_, sizeof(feat_t) * alphabet_size_);
        ind_probs_ = new feat_t[num_features_];
        memcpy(ind_probs_, other.ind_probs_, sizeof(feat_t) * num_features_);
        joint_probs_ = new feat_t[num_features_];
        memcpy(joint_probs_, other.joint_probs_, sizeof(feat_t) * num_features_);

    }

    BbcBuilder::BbcBuilder(BbcBuilder&& other) noexcept { // move constructor

        num_features_ = other.num_features_;
        nt_map_ = other.nt_map_;
        max_dist_ = other.max_dist_;
        alphabet_size_ = other.alphabet_size_;

        single_probs_ = other.single_probs_; other.single_probs_ = nullptr;
        single_cnts_ = other.single_cnts_; other.single_cnts_ = nullptr;
        ind_probs_ = other.ind_probs_; other.ind_probs_ = nullptr;
        joint_probs_ = other.joint_probs_; other.joint_probs_ = nullptr;

    }

    BbcBuilder& BbcBuilder::operator=(const BbcBuilder& other) { // copy assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] joint_probs_;
        delete[] ind_probs_;
        delete[] single_cnts_;
        delete[] single_probs_;

        // copy new resources
        num_features_ = other.num_features_;
        nt_map_ = other.nt_map_;
        max_dist_ = other.max_dist_;
        alphabet_size_ = other.alphabet_size_;

        single_probs_ = new feat_t[alphabet_size_];
        memcpy(single_probs_, other.single_probs_, sizeof(feat_t) * alphabet_size_);
        single_cnts_ = new feat_t[alphabet_size_];
        memcpy(single_cnts_, other.single_cnts_, sizeof(feat_t) * alphabet_size_);
        ind_probs_ = new feat_t[num_features_];
        memcpy(ind_probs_, other.ind_probs_, sizeof(feat_t) * num_features_);
        joint_probs_ = new feat_t[num_features_];
        memcpy(joint_probs_, other.joint_probs_, sizeof(feat_t) * num_features_);

        return *this;

    }

    BbcBuilder& BbcBuilder::operator=(BbcBuilder&& other) noexcept { // move assignment operator

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        // release old resources
        delete[] joint_probs_;
        delete[] ind_probs_;
        delete[] single_cnts_;
        delete[] single_probs_;

        // copy new resources
        num_features_ = other.num_features_;
        nt_map_ = other.nt_map_;
        max_dist_ = other.max_dist_;
        alphabet_size_ = other.alphabet_size_;

        single_probs_ = other.single_probs_; other.single_probs_ = nullptr;
        single_cnts_ = other.single_cnts_; other.single_cnts_ = nullptr;
        ind_probs_ = other.ind_probs_; other.ind_probs_ = nullptr;
        joint_probs_ = other.joint_probs_; other.joint_probs_ = nullptr;

        return *this;

    }

    BbcBuilder* BbcBuilder::clone() const {
        return new BbcBuilder(*this);
    }

    size_t BbcBuilder::num_features() {
        return num_features_;
    }

    std::vector<feat_t> BbcBuilder::get_features(const std::string& seq) {

        std::vector<feat_t> bbc_values(num_features_);
        compute_features(seq, bbc_values.data());
        return bbc_values;

    }

    void BbcBuilder::get_features(const std::string& seq, feat_t* features) {
        compute_features(seq, features);
    }

    void BbcBuilder::get_features(const std::string& seq, std::vector<feat_t>& features) {
        compute_features(seq, features.data());
    }

    void BbcBuilder::compute_features(const std::string& seq, feat_t* bbc_values) {

        std::fill(bbc_values, bbc_values + num_features_, 0);

        // compute probabilities of individual nucleotides in the sequence
        // and joint probabilities under the assumption of statistical independence
        std::fill(single_probs_, single_probs_ + alphabet_size_, 0.0);
        feat_t inc = static_cast<feat_t>(1.0) / seq.length();
        for (char c : seq) {
            single_probs_[nt_map_[c]] += inc;
        }
        for (auto i = 0; i < alphabet_size_; i++) {
            for (auto j = 0; j < alphabet_size_; j++) {
                ind_probs_[i * alphabet_size_ + j] = single_probs_[i] * single_probs_[j];
            }
        }

        // compute BBC features by iterating over the requested range of distances
        for (auto k = 1; k <= max_dist_ && k < seq.length(); k++) {

            // joint probabilities at distance k
            std::fill(single_cnts_, single_cnts_ + alphabet_size_, 0.0);
            std::fill(joint_probs_, joint_probs_ + num_features_, 0.0);
            for (auto i = 0; i < seq.length() - k; i++) {

                single_cnts_[nt_map_[seq[i]]] += 1;
                joint_probs_[nt_map_[seq[i]] * alphabet_size_ + nt_map_[seq[i + k]]] += 1;

            }
            for (auto i = 0; i < alphabet_size_; i++) {
                for (auto j = 0; j < alphabet_size_; j++) {
                    if (single_cnts_[i] > 0) joint_probs_[i * alphabet_size_ + j] *= single_probs_[i] / single_cnts_[i];
                }
            }

            for (size_t i = 0; i < num_features_; i++) {

                feat_t dev = joint_probs_[i] - ind_probs_[i]; // deviation from statistical independence
                // bbc_values[i] += dev + (dev * dev) / (2 * ind_probs_[i] * ind_probs_[i]); // paper contains error in development of equation 3
                if (ind_probs_[i] > 0) bbc_values[i] += dev + (dev * dev) / (2 * ind_probs_[i]);

            }

        }

        for (size_t i = 0; i < num_features_; i++) {
            bbc_values[i] /= std::log(2.0);
        }

    }


    /* === TwoDimWalkBuilder == */

    TwoDimWalkBuilder::TwoDimWalkBuilder(feat_t alpha, feat_t beta) {

        feat_t sqrt_beta = std::sqrt(beta);

        nt_map_x_['A'] = alpha;
        nt_map_x_['C'] = sqrt_beta;
        nt_map_x_['G'] = sqrt_beta;
        nt_map_x_['T'] = alpha;
        nt_map_x_['U'] = alpha;
        nt_map_y_['A'] = -sqrt_beta;
        nt_map_y_['C'] = alpha;
        nt_map_y_['G'] = -alpha;
        nt_map_y_['T'] = sqrt_beta;
        nt_map_y_['U'] = sqrt_beta;

    }

    TwoDimWalkBuilder* TwoDimWalkBuilder::clone() const {
        return new TwoDimWalkBuilder(*this);
    }

    size_t TwoDimWalkBuilder::num_features() {
        return num_features_;
    }

    std::vector<feat_t> TwoDimWalkBuilder::get_features(const std::string& seq) {

        std::vector<feat_t> features(num_features_);
        compute_features(seq, features.data());
        return features;

    }

    void TwoDimWalkBuilder::get_features(const std::string& seq, feat_t* features) {
        compute_features(seq, features);
    }

    void TwoDimWalkBuilder::get_features(const std::string& seq, std::vector<feat_t>& features) {
        compute_features(seq, features.data());
    }

    void TwoDimWalkBuilder::compute_features(const std::string& seq, feat_t* features) {

        feat_t prev_x = 0.0;
        feat_t prev_y = 0.0;
        feat_t sum_x = 0.0;
        feat_t sum_y = 0.0;

        for (char c : seq) {

            prev_x += nt_map_x_[c];
            prev_y += nt_map_y_[c];

            sum_x += prev_x;
            sum_y += prev_y;

        }

        features[0] = sum_x / seq.length();
        features[1] = sum_y / seq.length();

    }


    /* === CgrBuilder == */

    CgrBuilder::CgrBuilder() {

        nt_map_x_['A'] = 0;
        nt_map_x_['C'] = 0;
        nt_map_x_['G'] = 1;
        nt_map_x_['T'] = 1;
        nt_map_x_['U'] = 1;
        nt_map_y_['A'] = 0;
        nt_map_y_['C'] = 1;
        nt_map_y_['G'] = 1;
        nt_map_y_['T'] = 0;
        nt_map_y_['U'] = 0;

    }

    CgrBuilder* CgrBuilder::clone() const {
        return new CgrBuilder(*this);
    }

    size_t CgrBuilder::num_features() {
        return num_features_;
    }

    std::vector<feat_t> CgrBuilder::get_features(const std::string& seq) {

        std::vector<feat_t> features(num_features_);
        compute_features(seq, features.data());
        return features;

    }

    void CgrBuilder::get_features(const std::string& seq, feat_t* features) {
        compute_features(seq, features);
    }

    void CgrBuilder::get_features(const std::string& seq, std::vector<feat_t>& features) {
        compute_features(seq, features.data());
    }

    void CgrBuilder::compute_features(const std::string& seq, feat_t* features) {

        feat_t prev_x = 0.5;
        feat_t prev_y = 0.5;
        feat_t sum_x = 0.0;
        feat_t sum_y = 0.0;
        feat_t sq_sum = 0.0;

        for (char c : seq) {

            prev_x = (prev_x + nt_map_x_[c]) / static_cast<feat_t>(2);
            prev_y = (prev_y + nt_map_y_[c]) / static_cast<feat_t>(2);

            sum_x += prev_x;
            sum_y += prev_y;

            sq_sum += prev_x * prev_x + prev_y * prev_y;

        }

        features[0] = sum_x / seq.length(); // average x-coordinate
        features[1] = sum_y / seq.length(); // average y-coordinate
        features[2] = std::sqrt((sq_sum - seq.length() * (features[0] * features[0] + features[1] * features[1]))
                / (seq.length() - static_cast<feat_t>(1.0))); // mean deviation

    }


    /* === MultiCgrBuilder == */

    MultiCgrBuilder::MultiCgrBuilder() {

        nt_maps_x_[0]['A'] = 0; nt_maps_x_[0]['C'] = 0; nt_maps_x_[0]['G'] = 1; nt_maps_x_[0]['T'] = 1; nt_maps_x_[0]['U'] = 1; // CGR-RY
        nt_maps_x_[1]['A'] = 0; nt_maps_x_[1]['C'] = 1; nt_maps_x_[1]['G'] = 0; nt_maps_x_[1]['T'] = 1; nt_maps_x_[1]['U'] = 1; // CGR-RY
        nt_maps_x_[2]['A'] = 0; nt_maps_x_[2]['C'] = 0; nt_maps_x_[2]['G'] = 1; nt_maps_x_[2]['T'] = 1; nt_maps_x_[2]['U'] = 1; // CGR-WS

        nt_maps_y_[0]['A'] = 0; nt_maps_y_[0]['C'] = 1; nt_maps_y_[0]['G'] = 1; nt_maps_y_[0]['T'] = 0; nt_maps_y_[0]['U'] = 0; // CGR-RY
        nt_maps_y_[1]['A'] = 0; nt_maps_y_[1]['C'] = 1; nt_maps_y_[1]['G'] = 1; nt_maps_y_[1]['T'] = 0; nt_maps_y_[1]['U'] = 0; // CGR-RY
        nt_maps_y_[2]['A'] = 0; nt_maps_y_[2]['C'] = 1; nt_maps_y_[2]['G'] = 0; nt_maps_y_[2]['T'] = 1; nt_maps_y_[2]['U'] = 1; // CGR-WS

    }

    MultiCgrBuilder* MultiCgrBuilder::clone() const {
        return new MultiCgrBuilder(*this);
    }

    size_t MultiCgrBuilder::num_features() {
        return num_features_;
    }

    std::vector<feat_t> MultiCgrBuilder::get_features(const std::string& seq) {

        std::vector<feat_t> features(num_features_);
        compute_features(seq, features.data());
        return features;

    }

    void MultiCgrBuilder::get_features(const std::string& seq, feat_t* features) {
        compute_features(seq, features);
    }

    void MultiCgrBuilder::get_features(const std::string& seq, std::vector<feat_t>& features) {
        compute_features(seq, features.data());
    }

    void MultiCgrBuilder::compute_features(const std::string& seq, feat_t* features) {

        for (auto r = 0; r < 3; r++) {

            feat_t prev_x = 0.5;
            feat_t prev_y = 0.5;
            feat_t sum_x = 0.0;
            feat_t sum_y = 0.0;
            feat_t sq_sum = 0.0;

            auto& map_x = nt_maps_x_[r];
            auto& map_y = nt_maps_y_[r];

            for (char c : seq) {

                prev_x = (prev_x + map_x[c]) / static_cast<feat_t>(2);
                prev_y = (prev_y + map_y[c]) / static_cast<feat_t>(2);

                sum_x += prev_x;
                sum_y += prev_y;

                sq_sum += prev_x * prev_x + prev_y * prev_y;

            }

            features[3 * r] = sum_x / seq.length(); // average x-coordinate
            features[3 * r + 1] = sum_y / seq.length(); // average y-coordinate
            features[3 * r + 2] = std::sqrt((sq_sum - seq.length() * (features[3 * r] * features[3 * r] + features[3 * r + 1] * features[3 * r + 1]))
                    / (seq.length() - static_cast<feat_t>(1.0))); // mean deviation

        }

    }


    /* === ThreeDimCgrBuilder == */

    ThreeDimCgrBuilder::ThreeDimCgrBuilder() {

        nt_maps_x_[0]['A'] = 1; nt_maps_x_[0]['C'] = 1; nt_maps_x_[0]['G'] = 0; nt_maps_x_[0]['T'] = 0; nt_maps_x_[0]['U'] = 0; // A-CGT characteristic curve
        nt_maps_x_[1]['A'] = 1; nt_maps_x_[1]['C'] = 1; nt_maps_x_[1]['G'] = 0; nt_maps_x_[1]['T'] = 0; nt_maps_x_[1]['U'] = 0; // C-AGT characteristic curve
        nt_maps_x_[2]['A'] = 1; nt_maps_x_[2]['C'] = 0; nt_maps_x_[2]['G'] = 1; nt_maps_x_[2]['T'] = 0; nt_maps_x_[2]['U'] = 0; // G-ACT characteristic curve
        nt_maps_x_[3]['A'] = 1; nt_maps_x_[3]['C'] = 0; nt_maps_x_[3]['G'] = 0; nt_maps_x_[3]['T'] = 1; nt_maps_x_[3]['U'] = 1; // T-ACG characteristic curve

        nt_maps_y_[0]['A'] = 1; nt_maps_y_[0]['C'] = 0; nt_maps_y_[0]['G'] = 1; nt_maps_y_[0]['T'] = 0; nt_maps_y_[0]['U'] = 0; // A-CGT characteristic curve
        nt_maps_y_[1]['A'] = 0; nt_maps_y_[1]['C'] = 1; nt_maps_y_[1]['G'] = 1; nt_maps_y_[1]['T'] = 0; nt_maps_y_[1]['U'] = 0; // C-AGT characteristic curve
        nt_maps_y_[2]['A'] = 0; nt_maps_y_[2]['C'] = 1; nt_maps_y_[2]['G'] = 1; nt_maps_y_[2]['T'] = 0; nt_maps_y_[2]['U'] = 0; // G-ACT characteristic curve
        nt_maps_y_[3]['A'] = 0; nt_maps_y_[3]['C'] = 1; nt_maps_y_[3]['G'] = 0; nt_maps_y_[3]['T'] = 1; nt_maps_y_[3]['U'] = 1; // T-ACG characteristic curve

        nt_maps_z_[0]['A'] = 1; nt_maps_z_[0]['C'] = 0; nt_maps_z_[0]['G'] = 0; nt_maps_z_[0]['T'] = 1; nt_maps_z_[0]['U'] = 1; // A-CGT characteristic curve
        nt_maps_z_[1]['A'] = 0; nt_maps_z_[1]['C'] = 1; nt_maps_z_[1]['G'] = 0; nt_maps_z_[1]['T'] = 1; nt_maps_z_[1]['U'] = 1; // C-AGT characteristic curve
        nt_maps_z_[2]['A'] = 0; nt_maps_z_[2]['C'] = 0; nt_maps_z_[2]['G'] = 1; nt_maps_z_[2]['T'] = 1; nt_maps_z_[2]['U'] = 1; // G-ACT characteristic curve
        nt_maps_z_[3]['A'] = 0; nt_maps_z_[3]['C'] = 0; nt_maps_z_[3]['G'] = 1; nt_maps_z_[3]['T'] = 1; nt_maps_z_[3]['U'] = 1; // T-ACG characteristic curve

    }

    ThreeDimCgrBuilder* ThreeDimCgrBuilder::clone() const {
        return new ThreeDimCgrBuilder(*this);
    }

    size_t ThreeDimCgrBuilder::num_features() {
        return num_features_;
    }

    std::vector<feat_t> ThreeDimCgrBuilder::get_features(const std::string& seq) {

        std::vector<feat_t> features(num_features_);
        compute_features(seq, features.data());
        return features;

    }

    void ThreeDimCgrBuilder::get_features(const std::string& seq, feat_t* features) {
        compute_features(seq, features);
    }

    void ThreeDimCgrBuilder::get_features(const std::string& seq, std::vector<feat_t>& features) {
        compute_features(seq, features.data());
    }

    void ThreeDimCgrBuilder::compute_features(const std::string& seq, feat_t* features) {

        for (auto r = 0; r < 4; r++) {

            feat_t prev_x = 0.5;
            feat_t prev_y = 0.5;
            feat_t prev_z = 0.5;
            feat_t sum_x = 0.0;
            feat_t sum_y = 0.0;
            feat_t sum_z = 0.0;
            feat_t sq_sum = 0.0;

            auto& map_x = nt_maps_x_[r];
            auto& map_y = nt_maps_y_[r];
            auto& map_z = nt_maps_z_[r];

            for (char c : seq) {

                prev_x = (prev_x + map_x[c]) / static_cast<feat_t>(2);
                prev_y = (prev_y + map_y[c]) / static_cast<feat_t>(2);
                prev_z = (prev_z + map_z[c]) / static_cast<feat_t>(2);

                sum_x += prev_x;
                sum_y += prev_y;
                sum_z += prev_z;

                sq_sum += prev_x * prev_x + prev_y * prev_y + prev_z * prev_z;

            }

            features[4 * r] = sum_x / seq.length(); // average x-coordinate
            features[4 * r + 1] = sum_y / seq.length(); // average y-coordinate
            features[4 * r + 2] = sum_z / seq.length(); // average z-coordinate
            features[4 * r + 3] = std::sqrt((sq_sum - seq.length()
                    * (features[4 * r] * features[4 * r] + features[4 * r + 1] * features[4 * r + 1] + features[4 * r + 2] * features[4 * r + 2]))
                            / (seq.length() - static_cast<feat_t>(1.0))); // mean deviation

        }

    }

}