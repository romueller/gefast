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

#ifndef GEFAST_CLUSTERERS_HPP
#define GEFAST_CLUSTERERS_HPP

#include "Base.hpp"

namespace GeFaST {

    /*
     * Implements Swarm's iterative amplicon clustering strategy in an abstract way.
     * The used distance function and auxiliary data structures depend on the mode / configuration.
     */
    class ClassicSwarmer : public Clusterer {

    public:
        ClassicSwarmer* clone() const override; // deep-copy clone method

        /*
         * Cluster amplicons according to Swarm's iterative strategy.
         */
        SwarmStorage* cluster(const AmpliconStorage& amplicon_storage, const Configuration& config) override;

    };

    /*
     * Dereplicates amplicons based on their sequence.
     */
    class Dereplicator : public Clusterer {

    public:
        Dereplicator* clone() const override; // deep-copy clone method

        /*
         * Generate swarms such that all amplicons in each swarm have the same sequence.
         */
        SwarmStorage* cluster(const AmpliconStorage& amplicon_storage, const Configuration& config) override;

    };

    /*
     * Clusterer that does nothing and returns an empty SwarmStorage instance.
     *
     * Used in preprocessing-only (prep) mode.
     */
    class IdleClusterer : public Clusterer {

    public:
        IdleClusterer* clone() const override; // deep-copy clone method

        /*
         * Do not generate any swarm.
         */
        SwarmStorage* cluster(const AmpliconStorage& amplicon_storage, const Configuration& config) override;

    };



    /* ===== Consistency-guided refinement methods =====
     *
     * The following clusterers are adaptations of the DADA2 implementation:
     * "DADA2: High-resolution sample inference from Illumina amplicon data"
     * (Callahan et al., 2016, https://doi.org/10.1038/nmeth.3869)
     *
     * Implementation obtained from: https://github.com/benjjneb/dada2 (especially from cluster.cpp and pval.cpp)
     */

    /*
     * Basic reimplementation of DADA2's split-and-shuffle strategy as an actual clustering method.
     */
    class Dada2Clusterer : public Clusterer {

    public:
        Dada2Clusterer() = default;

        Dada2Clusterer(const Configuration& config);

        Dada2Clusterer* clone() const override; // deep-copy clone method

        /*
         * Cluster amplicons according to DADA2's strategy and create swarms from the results.
         */
        SwarmStorage* cluster(const AmpliconStorage& amplicon_storage, const Configuration& config) override;

        /*
         * Apply DADA2's split-and-shuffle strategy to given group of amplicons and create swarms from the results.
         * Used by some of the consistency-guided cluster refiners as one of the processing options
         * for unattached light swarms.
         */
        static void add_swarms(const AmpliconCollection& ac, std::vector<numSeqs_t>& ampl_ids, Swarms& swarms, const Configuration& config,
                               bool use_k_mers, double k_dist_cutoff, int band_size, double omega_a,
                               double omega_p, int max_clust, double min_fold,
                               int min_hamming, int min_abund, int homo_gap, bool gapless, bool greedy, std::string matrix_file);

    private:
        typedef Dada2Utility::RawSequence RawSequence;
        typedef Dada2Utility::Comparison Comparison;
        typedef Dada2Utility::Sub Sub;

        // configuration of clustering
        bool use_k_mers = true;         // flag indicating whether a k-mer distance screen precedes the alignment computation
        double k_dist_cutoff = 0.42;    // k-mer distance threshold (experimentally calibrated by Callahan et al. for DADA2)
        int band_size = 16;             // number of bands per side in banded Needleman-Wunsch alignment computation
        double omega_a = 1e-40;         // threshold for calling an amplicon significantly overabundant and, thus, inconsistent
        double omega_p = 1e-4;          // relaxed threshold for unique sequences with prior evidence of existence
        int max_clust = 0;              // maximum number of partitions (ignored when set to 0)
        double min_fold = 1.0;          // minimum fold-overabundance for sequences to form new partitions
        int min_hamming = 1;            // minimum Hamming-separation for sequences to form new partitions
        int min_abund = 1;              // minimum abundance for unique sequences to form new partitions
        int homo_gap = 0;               // gap cost in homopolymer regions, default of 0 indicates that those gaps are treated as normal gaps
        bool gapless = true;            // flag indicating whether using the k-mer sequence (instead of counts) and gapless alignments (only relevant when use_k_mers is true)
        bool greedy = true;             // flag indicating whether the greedy speedup is activated
        std::string matrix_file = "";   // path to file containing the error matrix


        // Auxiliary representation of clusters / partitions during the refinement
        // containing additional information on its members.
        struct Cluster {

            char* seq;                      // representative sequence for the cluster
            RawSequence* centre;            // representative raw for the cluster (corresponds to seq)
            unsigned int size;              // number of raws in this cluster
            unsigned int mass;              // sum of abundances of raws in this cluster
            unsigned int i;                 // the cluster number in the total clustering
            RawSequence** raws;             // array of pointers to raws of cluster members
            unsigned int max_raws;          // number of member pointers currently allocated for in raws
            bool update_e;                  // set to true when consensus changes and when raws are shuffled
            bool check_locks;               // set to true when raws should be checked for locking to this cluster
            std::vector<Comparison> comps;  // all comparisons with raws that could potentially join this cluster

            Cluster(unsigned int num);

            ~Cluster();

            // add a raw to this cluster, update partition characteristics and return index to the raw
            unsigned int add_raw(RawSequence* raw);

            // iterate over the raws in the cluster and update mass and size
            void census();

            // calculate and assign the centre raw (most abundant one) and the cluster sequence (the one of the centre)
            void assign_centre();

            // find the abundance p-value from the raw in the cluster
            double get_p_a(RawSequence* raw);

            /* Remove the raw from the cluster, update partition characteristics and return pointer to that raw.
             *
             * Sets update_e flag.
             * Called in shuffle() method of ClusterCollection.
             *
             * Important note: Reorders the list of raws and, thus, can break indices.
             */
            RawSequence* pop_raw(unsigned int r);

        };

        // holds all clusters (the full partitioning)
        struct ClusterCollection {

            unsigned int num_clusters;  // current number of clusters in the clustering
            unsigned int num_raws;      // number of raws in the clustering
            unsigned int total_mass;    // sum of abundances of the raws in the clustering
            unsigned int max_clusters;  // number of cluster pointers currently allocated for in clusters
            double omega_a;             // significance threshold for abundance p-value
            double omega_p;             // significance threshold for prior-abundance p-value
            RawSequence** raws;         // pointers to raw sequences (no ownership)
            Cluster** clusters;         // pointers to clusters (ownership)

            int* tvec;                  // auxiliary array storing nucleotide transitions (-> compute_lambda(...))
            unsigned int* qind;         // auxiliary array storing quality scores (-> compute_lambda(...))


            ClusterCollection(RawSequence** raws, unsigned int raws_size, double om_a, double om_p);

            ~ClusterCollection();

            // add a new, empty cluster and return the index of the added cluster in the clustering
            unsigned int add_cluster();


            // calculate the abundance p-value for each raw in the clustering
            void p_update(bool greedy);

            // find the minimum p-value and, if significant, create a new cluster by moving the corresponding raw there,
            // return the index of the new cluster (resp. 0 if no new cluster was added)
            int bud(double min_fold, int min_hamming, int min_abund);

            // move each sequence to the cluster that produces the highest expected number of it (only centres are stationary)
            bool shuffle();

            // perform alignments and compute lambda for all raws to its current cluster
            void compare(unsigned int i, std::vector<std::vector<double>>& error_matrix,
                         int match, int mismatch, int gap_pen, int homo_gap_pen,
                         bool use_k_mers, double k_dist_cutoff, int band_size, bool gapless, bool greedy);

        };

        // central step of clustering phase executing the overall split-and-shuffle process
        static ClusterCollection* run_dada(RawSequence** raws, int num_raw, std::vector<std::vector<double>>& error_matrix,
                                          int match, int mismatch, int gap_pen, int homo_gap_pen, bool use_k_mers, double k_dist_cutoff, int band_size,
                                          double omega_a, double omega_p, int max_clust, double min_fold, int min_hamming, int min_abund,
                                          bool gapless, bool greedy);

    };

    /*
     * Implements Swarm's iterative amplicon clustering strategy in an abstract way
     * and coupled with a consistency check involving an abundance p-value inspired by DADA2.
     * The used distance function and auxiliary data structures depend on the mode / configuration.
     */
    class ConsistentClassicSwarmer : public Clusterer {

    public:
        ConsistentClassicSwarmer() = default;

        ConsistentClassicSwarmer(const Configuration& config);

        ConsistentClassicSwarmer* clone() const override; // deep-copy clone method

        /*
         * Cluster amplicons according to Swarm's iterative strategy as in ClassicSwarmer,
         * extended by a check of an abundance p-value similar to DADA2 (but in reverse direction).
         * For a subseed and a candidate amplicon falling within the clustering threshold (that would be
         * added to the current swarm immediately in the original swarming approach), it is now checked whether
         * the probability of deriving as many (or more) copies of the candidate from the subseed is significantly low.
         * In this case, the link between the subseed and the candidate is not established.
         */
        SwarmStorage* cluster(const AmpliconStorage& amplicon_storage, const Configuration& config) override;

    private:
        // configuration of clustering
        bool use_k_mers = true;         // flag indicating whether a k-mer distance screen precedes the alignment computation
        double k_dist_cutoff = 0.42;    // k-mer distance threshold (experimentally calibrated by Callahan et al. for DADA2)
        int band_size = 16;             // number of bands per side in banded Needleman-Wunsch alignment computation
        double omega_a = 1e-40;         // threshold for calling an amplicon significantly overabundant and, thus, inconsistent
        double min_fold = 1.0;          // minimum fold-overabundance for sequences to form new partitions
        int min_hamming = 1;            // minimum Hamming-separation for sequences to form new partitions
        bool gapless = true;            // flag indicating whether using the k-mer sequence (instead of counts) and gapless alignments (only relevant when use_k_mers is true)
        std::string matrix_file = "";   // path to file containing the error matrix

    };

    /*
     * Implements Swarm's iterative amplicon clustering strategy but instead of using a distance threshold,
     * an abundance p-value significance threshold is used to decide whether an amplicon is added.
     */
    class ConsistencySwarmer : public Clusterer {

    public:
        ConsistencySwarmer() = default;

        ConsistencySwarmer(const Configuration& config);

        ConsistencySwarmer* clone() const override; // deep-copy clone method

        /*
         * Cluster amplicons according to Swarm's iterative strategy but using an abundance p-value similar
         * to DADA2 (but in reverse direction). For a subseed and a candidate amplicon, it is checked whether
         * the probability of deriving as many (or more) copies of the candidate from the subseed is significantly low.
         * In this case, the link between the subseed and the candidate is not established.
         * Since singletons are never significant and, thus, would be always added to the first cluster,
         * they remain unswarmed and are rather grafted to the clusters in the refinement phase
         * using one of the DADA2-inspired methods.
         */
        SwarmStorage* cluster(const AmpliconStorage& amplicon_storage, const Configuration& config) override;

    private:
        // configuration of clustering
        bool use_k_mers = true;         // flag indicating whether a k-mer distance screen precedes the alignment computation
        double k_dist_cutoff = 0.42;    // k-mer distance threshold (experimentally calibrated by Callahan et al. for DADA2)
        int band_size = 16;             // number of bands per side in banded Needleman-Wunsch alignment computation
        double omega_a = 1e-40;         // threshold for calling an amplicon significantly overabundant and, thus, inconsistent
        double min_fold = 1.0;          // minimum fold-overabundance for sequences to form new partitions
        int min_hamming = 1;            // minimum Hamming-separation for sequences to form new partitions
        bool gapless = true;            // flag indicating whether using the k-mer sequence (instead of counts) and gapless alignments (only relevant when use_k_mers is true)
        std::string matrix_file = "";   // path to file containing the error matrix

    };

}

#endif //GEFAST_CLUSTERERS_HPP
