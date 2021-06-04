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

#ifndef GEFAST_FEATUREBUILDERS_HPP
#define GEFAST_FEATUREBUILDERS_HPP

#include "Base.hpp"
#include "Utility.hpp"

namespace GeFaST {

    /*
     * Word-composition vector (wcv) representation
     * TODO? absolute vs relative counts
     */
    class WordCompBuilder : public FeatureBuilder {

    public:
        /*
         * Prepares the counting process by establishing a map between the alphabet and integers
         * and by determining the number of features based on the alphabet size and the word length k.
         */
        WordCompBuilder(const std::string& alphabet, short k);

        WordCompBuilder(const WordCompBuilder& other) = default; // copy constructor

        WordCompBuilder(WordCompBuilder&& other) = default; // move constructor

        WordCompBuilder& operator=(const WordCompBuilder& other) = default; // copy assignment operator

        WordCompBuilder& operator=(WordCompBuilder&& other) = default; // move assignment operator

        WordCompBuilder* clone() const override; // deep-copy clone method


        size_t num_features() override;

        std::vector<feat_t> get_features(const std::string& seq) override;
        void get_features(const std::string& seq, std::vector<feat_t>& features) override;
        void get_features(const std::string& seq, feat_t* features) override;

    private:
        /*
         * Count the number of occurrences of the different words (k-mers) in the given string.
         */
        inline void compute_features(const std::string& seq, feat_t* counts);

        size_t num_features_; // number of features, depends on word size and alphabet size
        std::array<short, 256> nt_map_; // map translating alphabet letters to integers
        short k_; // word size / k-mer length
        size_t alphabet_size_;

    };


    /*
     * Category-Position-Frequency (cpf) representation
     *
     * Based on:
     * "An improved alignment-free model for dna sequence similarity metric"
     * (Bao et al., 2014, https://doi.org/10.1186/1471-2105-15-321)
     */
    class CpfBuilder : public FeatureBuilder {

    public:
        /*
         * Prepares the feature computation by initialising reused auxiliary arrays
         * based on the maximum sequence length.
         */
        CpfBuilder(lenSeqs_t max_len);

        ~CpfBuilder();

        CpfBuilder(const CpfBuilder& other); // copy constructor

        CpfBuilder(CpfBuilder&& other) noexcept; // move constructor

        CpfBuilder& operator=(const CpfBuilder& other); // copy assignment operator

        CpfBuilder& operator=(CpfBuilder&& other) noexcept; // move assignment operator

        CpfBuilder* clone() const override; // deep-copy clone method


        size_t num_features() override;

        std::vector<feat_t> get_features(const std::string& seq) override;
        void get_features(const std::string& seq, std::vector<feat_t>& features) override;
        void get_features(const std::string& seq, feat_t* features) override;

    private:
        /*
         * Compute the CPF representation of the given string.
         */
        inline void compute_features(const std::string& seq, feat_t* entropies);

        size_t num_features_ = 12; // number of features
        lenSeqs_t max_len_; // maximum expected sequence length
        feat_t* psums_ry_; // auxiliary array for partial sums related to RY-categorisation
        feat_t* psums_mk_; // auxiliary array for partial sums related to MK-categorisation
        feat_t* psums_ws_; // auxiliary array for partial sums related to WS-categorisation
        feat_t* lasts_ry_; // auxiliary array for last 1 in indicator sequences related to RY-categorisation
        feat_t* lasts_mk_; // auxiliary array for last 1 in indicator sequences related to MK-categorisation
        feat_t* lasts_ws_; // auxiliary array for last 1 in indicator sequences related to WS-categorisation

        // Dinucleotide mappings based on the chemical properties of the nucleotides.
        // There are two options per property, leading to 4 combinations per property.
        static std::unordered_map<std::string, short> ry_map; // purine R = {A,G} vs pyrimidine Y = {C,T}
        static std::unordered_map<std::string, short> mk_map; // amino M = {A,C} vs carbonyl / keto K = {G,T}
        static std::unordered_map<std::string, short> ws_map; // weak hydrogen bond W = {A,T} vs strong hydrogen bond S = {C,G}

    };


    /*
     * Directed-Euler-Tour (det) representation
     *
     * Based on:
     * "A novel model for DNA sequence similarity analysis based on graph theory"
     * (Qi et al., 2011, http://dx.doi.org/10.4137/EBO.S7364)
     */
    class DetBuilder : public FeatureBuilder {

    public:
        /*
         * Prepares the feature computation by establishing a map between the alphabet and integers
         * and by determining the number of features based on the alphabet size.
         * Also sets the used decay rate alpha (directly or derived from the preference distance).
         */
        DetBuilder(const std::string& alphabet, feat_t alpha, lenSeqs_t max_length);
        DetBuilder(const std::string& alphabet, lenSeqs_t pref_dist, lenSeqs_t max_length);

        DetBuilder(const DetBuilder& other) = default; // copy constructor

        DetBuilder(DetBuilder&& other) = default; // move constructor

        DetBuilder& operator=(const DetBuilder& other) = default; // copy assignment operator

        DetBuilder& operator=(DetBuilder&& other) = default; // move assignment operator

        DetBuilder* clone() const override; // deep-copy clone method


        size_t num_features() override;

        std::vector<feat_t> get_features(const std::string& seq) override;
        void get_features(const std::string& seq, std::vector<feat_t>& features) override;
        void get_features(const std::string& seq, feat_t* features) override;

    private:
        /*
         * Compute the DET representation of the given string.
         */
        inline void compute_features(const std::string& seq, feat_t* arc_weights);

        size_t num_features_; // number of features, depends on alphabet size
        std::unordered_map<std::string, short> pair_map_; // map translating letter pairs to integers
        feat_t alpha_; // decay rate

        std::vector<feat_t> weights_;

    };


    /*
     * Base-base correlation (bbc) representation
     *
     * Based on:
     * "A novel feature-based method for whole genome phylogenetic analysis without alignment"
     * (Liu et al., 2008, https://doi.org/10.1016/j.bbrc.2008.01.070)
     */
    class BbcBuilder : public FeatureBuilder {

    public:
        /*
         * Prepares the feature computation by establishing a map between the alphabet and integers,
         * initialising resued auxiliary arrays and by determining the number of features based on the alphabet size.
         * Also sets the maximum considered distance.
         */
        BbcBuilder(const std::string& alphabet, short max_dist);

        ~BbcBuilder();

        BbcBuilder(const BbcBuilder& other); // copy constructor

        BbcBuilder(BbcBuilder&& other) noexcept; // move constructor

        BbcBuilder& operator=(const BbcBuilder& other); // copy assignment operator

        BbcBuilder& operator=(BbcBuilder&& other) noexcept; // move assignment operator

        BbcBuilder* clone() const override; // deep-copy clone method


        size_t num_features() override;

        std::vector<feat_t> get_features(const std::string& seq) override;
        void get_features(const std::string& seq, std::vector<feat_t>& features) override;
        void get_features(const std::string& seq, feat_t* features) override;

    private:
        /*
         * Compute the BBC representation of the given string.
         */
        inline void compute_features(const std::string& seq, feat_t* arc_weights);

        size_t num_features_; // number of features, depends on alphabet size
        feat_t* single_probs_; // auxiliary array for probabilities of nucleotides
        feat_t* single_cnts_; // auxiliary array for nucleotide counts while considering a single distance
        feat_t* ind_probs_; // auxiliary array for joint probabilities under the assumption of statistical independence
        feat_t* joint_probs_; // auxiliary array for observed joint probabilities
        std::array<short, 256> nt_map_; // map translating alphabet letters to integers
        short max_dist_; // maximum considered distance between nucleotides
        short alphabet_size_;

    };


    /*
     * Two-dimensional walk (2d-w) representation
     *
     * Based on:
     * "A 2D graphical representation of DNA sequence"
     * (Liao, 2005, https://doi.org/10.1016/j.cplett.2004.11.059)
     */
    class TwoDimWalkBuilder : public FeatureBuilder {

    public:
        /*
         * Prepares the feature computation by establishing a map between the alphabet
         * and movements along the axes based on alpha and beta.
         */
        TwoDimWalkBuilder(feat_t alpha, feat_t beta);

        TwoDimWalkBuilder(const TwoDimWalkBuilder& other) = default; // copy constructor

        TwoDimWalkBuilder(TwoDimWalkBuilder&& other) = default; // move constructor

        TwoDimWalkBuilder& operator=(const TwoDimWalkBuilder& other) = default; // copy assignment operator

        TwoDimWalkBuilder& operator=(TwoDimWalkBuilder&& other) = default; // move assignment operator

        TwoDimWalkBuilder* clone() const override; // deep-copy clone method


        size_t num_features() override;

        std::vector<feat_t> get_features(const std::string& seq) override;
        void get_features(const std::string& seq, std::vector<feat_t>& features) override;
        void get_features(const std::string& seq, feat_t* features) override;

    private:
        /*
         * Compute the two-dimensional walk representation of the given string.
         */
        inline void compute_features(const std::string& seq, feat_t* features);

        size_t num_features_ = 2; // number of features
        std::array<feat_t, 256> nt_map_x_; // map translating letters into movements along the x-axis
        std::array<feat_t, 256> nt_map_y_; // map translating letters into movements along the y-axis

    };


    /*
     * Chaos Game Representation (cgr) representation
     *
     * Based on:
     * "Analysis of genomic sequences by Chaos Game Representation"
     * (Almeida et al., 2001, https://doi.org/10.1093/bioinformatics/17.5.429)
     * and
     * "Similarities of DNA sequences based on 3D Chaos Game Representation"
     * (Huang & Shi, 2010, https://doi.org/10.1109/BMEI.2010.5639720)
     *
     * Fixes the wrong sign in Equation 1 in Almeida et al. in order to actually move half the distance
     * between the previous position and the corner corresponding to the current symbol.
     *
     * In contrast to Almeida et al., we do not extract a frequency matrix from the CGR
     * but compute its geometric centre and mean deviation to obtain a representation
     * similar to the one in Huang & Shi.
     */
    class CgrBuilder : public FeatureBuilder {

    public:
        CgrBuilder();

        CgrBuilder(const CgrBuilder& other) = default; // copy constructor

        CgrBuilder(CgrBuilder&& other) = default; // move constructor

        CgrBuilder& operator=(const CgrBuilder& other) = default; // copy assignment operator

        CgrBuilder& operator=(CgrBuilder&& other) = default; // move assignment operator

        CgrBuilder* clone() const override; // deep-copy clone method


        size_t num_features() override;

        std::vector<feat_t> get_features(const std::string& seq) override;
        void get_features(const std::string& seq, std::vector<feat_t>& features) override;
        void get_features(const std::string& seq, feat_t* features) override;

    private:
        /*
         * Compute the CGR representation of the given string.
         */
        inline void compute_features(const std::string& seq, feat_t* features);

        size_t num_features_ = 3; // number of features
        std::array<feat_t, 256> nt_map_x_; // map translating letters into x-axis coordinates
        std::array<feat_t, 256> nt_map_y_; // map translating letters into y-axis coordinates

    };


    /*
     * Multipe CGR (mcgr) representation
     *
     * Based on:
     * "Analysis of similarity / dissimilarity of DNA sequences based on Chaos Game Representation"
     * (Deng & Luan, 2013, http://dx.doi.org/10.1155/2013/926519)
     * and
     * "Similarities of DNA sequences based on 3D Chaos Game Representation"
     * (Huang & Shi, 2010, https://doi.org/10.1109/BMEI.2010.5639720)
     *
     * Uses the idea to create three representations based on the chemical properties
     * of the nucleotides (purine / pyrimidine, amino / keto, weak / strong) from Deng & Luan,
     * but does not combine the coordinates into walk sequences to compute Hurst exponents from them.
     * Instead, geometric centres and mean deviations are computed to obtain a representation as in Huang & Shi.
     */
    class MultiCgrBuilder : public FeatureBuilder {

    public:
        MultiCgrBuilder();

        MultiCgrBuilder(const MultiCgrBuilder& other) = default; // copy constructor

        MultiCgrBuilder(MultiCgrBuilder&& other) = default; // move constructor

        MultiCgrBuilder& operator=(const MultiCgrBuilder& other) = default; // copy assignment operator

        MultiCgrBuilder& operator=(MultiCgrBuilder&& other) = default; // move assignment operator

        MultiCgrBuilder* clone() const override; // deep-copy clone method


        size_t num_features() override;

        std::vector<feat_t> get_features(const std::string& seq) override;
        void get_features(const std::string& seq, std::vector<feat_t>& features) override;
        void get_features(const std::string& seq, feat_t* features) override;

    private:
        /*
         * Compute the MCGR representation of the given string.
         */
        inline void compute_features(const std::string& seq, feat_t* features);

        size_t num_features_ = 9; // number of features
        std::array<std::array<feat_t, 256>, 3> nt_maps_x_; // maps translating letters into x-axis coordinates
        std::array<std::array<feat_t, 256>, 3> nt_maps_y_; // maps translating letters into y-axis coordinates

    };


    /*
     * Three-dimensional CGR (3d-cgr) representation
     *
     * Based on:
     * "Similarities of DNA sequences based on 3D Chaos Game Representation"
     * (Huang & Shi, 2010, https://doi.org/10.1109/BMEI.2010.5639720)
     *
     * Uses the proposed idea to compute geometric centres and mean deviations to obtain a representation.
     */
    class ThreeDimCgrBuilder : public FeatureBuilder {

    public:
        ThreeDimCgrBuilder();

        ThreeDimCgrBuilder(const ThreeDimCgrBuilder& other) = default; // copy constructor

        ThreeDimCgrBuilder(ThreeDimCgrBuilder&& other) = default; // move constructor

        ThreeDimCgrBuilder& operator=(const ThreeDimCgrBuilder& other) = default; // copy assignment operator

        ThreeDimCgrBuilder& operator=(ThreeDimCgrBuilder&& other) = default; // move assignment operator

        ThreeDimCgrBuilder* clone() const override; // deep-copy clone method


        size_t num_features() override;

        std::vector<feat_t> get_features(const std::string& seq) override;
        void get_features(const std::string& seq, std::vector<feat_t>& features) override;
        void get_features(const std::string& seq, feat_t* features) override;

    private:
        /*
         * Compute the three-dimensional CGR representation of the given string.
         */
        inline void compute_features(const std::string& seq, feat_t* features);

        size_t num_features_ = 16; // number of features
        std::array<std::array<feat_t, 256>, 4> nt_maps_x_; // maps translating letters into x-axis coordinates
        std::array<std::array<feat_t, 256>, 4> nt_maps_y_; // maps translating letters into y-axis coordinates
        std::array<std::array<feat_t, 256>, 4> nt_maps_z_; // maps translating letters into z-axis coordinates

    };

}

#endif //GEFAST_FEATUREBUILDERS_HPP
