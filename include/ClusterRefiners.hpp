/*
 * GeFaST
 *
 * Copyright (C) 2016 - 2020 Robert Mueller
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

#ifndef GEFAST_CLUSTERREFINERS_HPP
#define GEFAST_CLUSTERREFINERS_HPP

#include "Base.hpp"

namespace GeFaST {

    /*
     * Representation of a grafting candidate.
     * Parent and (potential) child amplicon are represented by the member of their respective swarm.
     * If the grafting takes place, the swarm of the child is supposed to be grafted upon the parent's swarm.
     */
    struct GraftCandidate {

        Swarm* parent_swarm; // parent swarm (no ownership)
        numSeqs_t parent_id; // swarm-internal id

        Swarm* child_swarm; // child swarm (no ownership)
        numSeqs_t child_id; // swarm-internal id

        dist_t distance; // distance between parent amplicon and child amplicon

        GraftCandidate(); // empty constructor for null initialisation

        /*
         * Initialise with parent swarm, parent amplicon, child swarm, child amplicon
         * and the distance between the two amplicons.
         */
        GraftCandidate(Swarm* ps, numSeqs_t p, Swarm* cs, numSeqs_t c, dist_t d);

    };


    /*
     * Sort grafting candidates by the abundances of the parents and, if necessary, break ties through
     * the abundances of the children (both descending).
     * Use the lexicographical order of the ids of the amplicons as the tie-breaker for the abundance comparison.
     */
    struct CompareGraftCandidatesAbund {

        explicit CompareGraftCandidatesAbund(const AmpliconCollection& a);

        inline bool compare(const numSeqs_t ampl_a, const numSeqs_t ampl_b);

        bool operator()(const GraftCandidate& a, const GraftCandidate& b);

        const AmpliconCollection& ac; // reference to the collection from which the compared amplicons originate

    };


    /*
     * Implements a generalised version of Swarm's fastidious cluster refinement strategy in an abstract way.
     * The used distance function and auxiliary data structures depend on the mode / configuration.
     */
    class FastidiousRefiner : public ClusterRefiner {

    public:
        FastidiousRefiner* clone() const override; // deep-copy clone method

        /*
         * Graft light swarms onto heavy ones by enlarging the search radius (threshold).
         * To this end, index the amplicons from light swarms.
         * Then, search for matches between the amplicons from heavy swarms and the indexed amplicons.
         * Based on the resulting grafting-link candidates, light swarms are grafted upon heavy ones:
         *  - The (final) grafting partner of an amplicon from a light swarm is a matching amplicon with the highest abundance.
         *  - Each light swarm can be grafted upon at most one heavy swarm (even though there can be grafting candidates
         *    for several amplicons of the light swarm.
         *  - Grafting candidates with a higher parent amplicon abundance (an, for ties, higher child amplicon abundance)
         *    have a higher priority.
         */
        void refine(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage, const Configuration& config) override;

    private:
        /*
         * Compare the current candidate for the parent amplicon of a grafting link with another amplicon.
         * Return true iff ampl_a has a higher abundance or, in case of a tie, has the lexicographically smaller identifier.
         */
        inline bool compare_candidates(const AmpliconCollection& ac, const numSeqs_t ampl_a, const numSeqs_t ampl_b);

        /*
         * Refine the swarms of a particular amplicon pool as described above.
         */
        void refine_pool(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage, const numSeqs_t pool_id,
                         numSeqs_t& light_s, numSeqs_t& light_a, const Configuration& config);

    };


    /*
     * Implements an iterative version of Swarm's fastidious cluster refinement strategy in an abstract way.
     * The used distance function and auxiliary data structures depend on the mode / configuration.
     *
     * Performs several rounds of the classic approach with increasing threshold.
     */
    class IterativeFastidiousRefiner : public ClusterRefiner {

    public:
        IterativeFastidiousRefiner* clone() const override; // deep-copy clone method

        /*
         * Graft light swarms onto heavy ones in several rounds by successively increasing the threshold.
         * The individual rounds behave as described for FastidiousRefiner.
         * Light swarms already attached to a heavy one for some smaller threshold are not considered in
         * subsequent rounds.
         */
        void refine(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage, const Configuration& config) override;

    private:
        /*
         * Compare the current candidate for the parent amplicon of a grafting link with another amplicon.
         * Return true iff ampl_a has a higher abundance or, in case of a tie, has the lexicographically smaller identifier.
         */
        inline bool compare_candidates(const AmpliconCollection& ac, const numSeqs_t ampl_a, const numSeqs_t ampl_b);

        /*
         * Refine the swarms of a particular amplicon pool as described above using the given threshold.
         */
        void refine_pool(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage, const numSeqs_t pool_id,
                         const dist_t threshold, numSeqs_t& light_s, numSeqs_t& light_a, const Configuration& config);

    };


    /*
     * Refiner that does nothing.
     *
     * Used in preprocessing-only (prep) mode.
     */
    class IdleRefiner : public ClusterRefiner {

    public:
        IdleRefiner* clone() const override; // deep-copy clone method

        /*
         * Do not refine the collection swarms.
         */
        void refine(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage, const Configuration& config) override;

    };



    /* ===== Consistency-guided refinement methods =====
     *
     * The following refiners are adaptations of the DADA2 implementation:
     * "DADA2: High-resolution sample inference from Illumina amplicon data"
     * (Callahan et al., 2016, https://doi.org/10.1038/nmeth.3869)
     *
     * Implementation obtained from: https://github.com/benjjneb/dada2 (especially from cluster.cpp and pval.cpp)
     */

    /*
     * Consistency-guided cluster refinement considering light swarms as indivisible units.
     *
     * Options for further processing of light swarms unattached after the refinement:
     *  - 1 = keep the swarms without further changes
     *  - 2 = discard them
     *  - 3 = combine them into a single (star-shaped) swarm
     *  - 4 = combine and repartition them using DADA2's split-and-shuffle strategy
     *
     * By default, a limited form of greediness is activated:
     * If the mass of a light swarm is larger than the abundance of the seed of a heavy swarm,
     * the light swarm cannot be grafted onto that heavy swarm.
     */
    class LightSwarmAppender : public ClusterRefiner {

    public:
        LightSwarmAppender() = default;

        LightSwarmAppender(const Configuration& config);

        LightSwarmAppender* clone() const override; // deep-copy clone method

        /*
         * Graft light swarms onto heavy ones.
         * For each light swarm, its seed is compared to the seeds of the heavy swarms and checked for consistency.
         * Depending on the configuration, light swarms still unattached after refinement are further processed.
         * are kept without further changes, discarded
         * combined into a single swarm, or combined and repartitioned using DADA2's split-and-shuffle strategy.
         */
        void refine(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage, const Configuration& config) override;

    private:
        typedef Dada2Utility::RawSequence RawSequence;
        typedef Dada2Utility::Sub Sub;

        // configuration of refinement
        bool use_k_mers = false;        // flag indicating whether a k-mer distance screen precedes the alignment computation
        double k_dist_cutoff = 0.42;    // k-mer distance threshold (experimentally calibrated by Callahan et al. for DADA2)
        int band_size = 16;             // number of bands per side in banded Needleman-Wunsch alignment computation
        double omega_a = 1e-40;         // threshold for calling the light swarm significantly overabundant and, thus, inconsistent
        double omega_p = 1e-4;          // relaxed threshold for unique sequences with prior evidence of existence
        int max_clust = 0;              // maximum number of partitions (ignored when set to 0)
        double min_fold = 1.0;          // minimum fold-overabundance for sequences to form new partitions
        int min_hamming = 1;            // minimum Hamming-separation for sequences to form new partitions
        int min_abund = 1;              // minimum abundance for unique sequences to form new partitions
        int homo_gap = 0;               // gap cost in homopolymer regions, default of 0 indicates that those gaps are treated as normal gaps
        bool gapless = true;            // flag indicating whether using the k-mer sequence (instead of counts) and gapless alignments (only relevant when use_k_mers is true)
        bool greedy = true;             // flag indicating whether the greedy speedup is activated
        std::string matrix_file = "";   // path to file containing the error matrix
        int light_opt = 1;              // processing option for unattached light clusters

        /*
         * Auxiliary representation of target clusters (heavy swarms) during the refinement.
         * Contains additional information on the swarm and its seed to speed up the process
         * of finding the grafting targets.
         */
        struct TargetCluster {

            numSeqs_t sid; // pool-internal swarm id of this swarm
            RawSequence seed; // seed representation with additional information

            TargetCluster() = default;

            // initialise using swarm id, the corresponding amplicon collection, the seed integer id
            TargetCluster(numSeqs_t s, const AmpliconCollection& ac, numSeqs_t a, bool use_k_mers, QualityEncoding<>& qe);

        };

        // find the abundance p-value for the light swarm (represented by its seed)
        double get_p_a(numSeqs_t target_mass, RawSequence& ampl);

        // determine whether the (corrected) abundance p-value is significant
        bool is_significant(RawSequence& src, numSeqs_t target_mass, numSeqs_t num_amplicons);

    };


    /*
     * Consistency-guided cluster refinement considering the amplicons from light swarms independently.
     *
     * Options for further processing of light swarms unattached after the refinement:
     *  - 1 = keep the swarms without further changes
     *  - 2 = discard them
     *  - 3 = combine them into a single (star-shaped) swarm
     *  - 4 = combine and repartition them using DADA2's split-and-shuffle strategy
     *
     * By default, a limited form of greediness is activated:
     * If the abundance of an amplicon from a light swarm is larger than the abundance of the seed of a heavy swarm,
     * the light amplicon cannot be moved to that heavy swarm.
     */
    class LightSwarmResolver : public ClusterRefiner {//light swarms dismantled into separate entities that can be grafted independently

    public:
        LightSwarmResolver() = default;

        LightSwarmResolver(const Configuration& config);

        LightSwarmResolver* clone() const override; // deep-copy clone method

        /*
         * Graft light swarms onto heavy ones.
         * Each amplicon from a light swarm is compared to the seeds of the heavy swarms and checked for consistency.
         * If there is no consistent heavy swarm, the amplicon forms a new unattached light swarm.
         * Depending on the configuration, unattached light swarms are kept without further changes, discarded
         * combined into a single swarm, or combined and repartitioned using DADA2's split-and-shuffle strategy.
         */
        void refine(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage, const Configuration& config) override;

    private:
        typedef Dada2Utility::RawSequence RawSequence;
        typedef Dada2Utility::Sub Sub;

        // configuration of refinement
        bool use_k_mers = false;        // flag indicating whether a k-mer distance screen precedes the alignment computation
        double k_dist_cutoff = 0.42;    // k-mer distance threshold (experimentally calibrated by Callahan et al. for DADA2)
        int band_size = 16;             // number of bands per side in banded Needleman-Wunsch alignment computation
        double omega_a = 1e-40;         // threshold for calling the light swarm significantly overabundant and, thus, inconsistent
        double omega_p = 1e-4;          // relaxed threshold for unique sequences with prior evidence of existence
        int max_clust = 0;              // maximum number of partitions (ignored when set to 0)
        double min_fold = 1.0;          // minimum fold-overabundance for sequences to form new partitions
        int min_hamming = 1;            // minimum Hamming-separation for sequences to form new partitions
        int min_abund = 1;              // minimum abundance for unique sequences to form new partitions
        int homo_gap = 0;               // gap cost in homopolymer regions, default of 0 indicates that those gaps are treated as normal gaps
        bool gapless = true;            // flag indicating whether using the k-mer sequence (instead of counts) and gapless alignments (only relevant when use_k_mers is true)
        bool greedy = true;             // flag indicating whether the greedy speedup is activated
        std::string matrix_file = "";   // path to file containing the error matrix
        int light_opt = 1;              // processing option for unattached light clusters

        /*
         * Auxiliary representation of target clusters (heavy swarms) during the refinement.
         * Contains additional information on the swarm and its seed to speed up the process
         * of finding the grafting targets.
         */
        struct TargetCluster {

            numSeqs_t sid; // pool-internal swarm id of this swarm
            RawSequence seed; // seed representation with additional information

            TargetCluster() = default;

            // initialise using swarm id, the corresponding amplicon collection, the seed integer id
            TargetCluster(numSeqs_t s, const AmpliconCollection& ac, numSeqs_t a, bool use_k_mers, QualityEncoding<>& qe);

        };

        // find the abundance p-value for the amplicon from a light swarm
        double get_p_a(numSeqs_t target_mass, RawSequence& ampl);

        // determine whether the (corrected) abundance p-value is significant
        bool is_significant(RawSequence& src, numSeqs_t target_mass, numSeqs_t num_amplicons);

    };


    /*
     * Consistency-guided cluster refinement considering light swarms as loosely coupled groups that can be modified.
     *
     * By default, a limited form of greediness is activated:
     * If the abundance of an amplicon from a light swarm is larger than the abundance of the seed of a heavy swarm,
     * the light amplicon cannot be moved to that heavy swarm.
     */
    class LightSwarmShuffler : public ClusterRefiner {

    public:
        LightSwarmShuffler() = default;

        LightSwarmShuffler(const Configuration& config);

        LightSwarmShuffler* clone() const override; // deep-copy clone method

        /*
         * Graft light swarms onto heavy ones.
         * Amplicons from light swarms are again processed individually, but the light swarm is not disassembled right away.
         * Each amplicon from a light swarm is compared to the seeds of the heavy swarms and checked for consistency.
         * The seed of a light swarm can only be moved to another swarm when all other members have already been moved.
         * The remaining members of a light swarm (if any) are rearranged into a new star-shaped swarm with the seed as its centre.
         */
        void refine(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage, const Configuration& config) override;

    private:
        typedef Dada2Utility::RawSequence RawSequence;
        typedef Dada2Utility::Sub Sub;

        // configuration of refinement
        bool use_k_mers = false;        // flag indicating whether a k-mer distance screen precedes the alignment computation
        double k_dist_cutoff = 0.42;    // k-mer distance threshold (experimentally calibrated by Callahan et al. for DADA2)
        int band_size = 16;             // number of bands per side in banded Needleman-Wunsch alignment computation
        double omega_a = 1e-40;         // threshold for calling the light swarm significantly overabundant and, thus, inconsistent
        double min_fold = 1.0;          // minimum fold-overabundance for sequences to form new partitions
        int min_hamming = 1;            // minimum Hamming-separation for sequences to form new partitions
        int min_abund = 1;              // minimum abundance for unique sequences to form new partitions
        int homo_gap = 0;               // gap cost in homopolymer regions, default of 0 indicates that those gaps are treated as normal gaps
        bool gapless = true;            // flag indicating whether using the k-mer sequence (instead of counts) and gapless alignments (only relevant when use_k_mers is true)
        bool greedy = true;             // flag indicating whether the greedy speedup is activated
        std::string matrix_file = "";   // path to file containing the error matrix

        /*
         * Auxiliary representation of target clusters (heavy swarms) during the refinement.
         * Contains additional information on the swarm and its seed to speed up the process
         * of finding the grafting targets.
         */
        struct TargetCluster {

            numSeqs_t sid; // pool-internal swarm id of this swarm
            RawSequence seed; // seed representation with additional information

            TargetCluster() = default;

            // initialise using swarm id, the corresponding amplicon collection, the seed integer id
            TargetCluster(numSeqs_t s, const AmpliconCollection& ac, numSeqs_t a, bool use_k_mers, QualityEncoding<>& qe);

        };

        // find the abundance p-value for the amplicon from a light swarm
        double get_p_a(numSeqs_t target_mass, RawSequence& ampl);

        // determine whether the (corrected) abundance p-value is significant
        bool is_significant(RawSequence& src, numSeqs_t target_mass, numSeqs_t num_amplicons);

    };


    /*
     * Adaptation of DADA2's split-and-shuffle strategy to a shuffle-only refinement of the initial clustering.
     * Instead of starting the process with a single partition (cluster) containing all amplicons,
     * the process is initialised with a partitioning corresponding to the swarms determined in the clustering phase.
     * In contrast to the other consistency-guided refiners, there is no distinction between light and heavy swarms.
     */
    class SwarmShuffler : public ClusterRefiner {

    public:
        SwarmShuffler() = default;

        SwarmShuffler(const Configuration& config);

        SwarmShuffler* clone() const override; // deep-copy clone method

        /*
         * Refine all swarms. Amplicons are again processed individually and shuffled between the partitions
         * until overall consistency or the maximum number of rounds is reached.
         */
        void refine(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage, const Configuration& config) override;

    private:
        typedef Dada2Utility::RawSequence RawSequence;
        typedef Dada2Utility::Comparison Comparison;
        typedef Dada2Utility::Sub Sub;

        // configuration of refinement
        bool use_k_mers = false;        // flag indicating whether a k-mer distance screen precedes the alignment computation
        double k_dist_cutoff = 0.42;    // k-mer distance threshold (experimentally calibrated by Callahan et al. for DADA2)
        int band_size = 16;             // number of bands per side in banded Needleman-Wunsch alignment computation
        double omega_a = 1e-40;         // threshold for calling the light swarm significantly overabundant and, thus, inconsistent
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

            // find the abundance p-value from the raw in the cluster
            double get_p_a(RawSequence* raw);

            // determine whether the (corrected) abundance p-value of the raw is significant
            bool is_significant(RawSequence& src, numSeqs_t num_amplicons, double min_fold, int min_hamming, int min_abund, double omega_a);

            /*
             * Remove the raw from the cluster, update partition characteristics and return pointer to that raw.
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


            ClusterCollection(RawSequence** raws, Swarms& initial_swarms, unsigned int raws_size, double om_a, double om_p);

            ~ClusterCollection();


            // calculate the abundance p-value for each raw in the clustering
            void p_update(bool greedy);

            // move each sequence to the cluster that produces the highest expected number of it (only centres are stationary)
            bool shuffle(double min_fold, int min_hamming, int min_abund);

            // perform alignments and compute lambda for all raws to the specified cluster
            void compare(unsigned int i, std::vector<std::vector<double>>& error_matrix,
                         int match, int mismatch, int gap_pen, int homo_gap_pen,
                         bool use_k_mers, double k_dist_cutoff, int band_size, bool gapless, bool greedy);

            // perform alignments and compute lambda for all raws to its current cluster
            void compare(std::vector<std::vector<double>>& error_matrix,
                         int match, int mismatch, int gap_pen, int homo_gap_pen,
                         bool use_k_mers, double k_dist_cutoff, int band_size, bool gapless, bool greedy);

        };

        // central step of refinement phase executing the overall shuffling
        ClusterCollection* shuffle(RawSequence** raws, Swarms& initial_swarms, int num_raw,
                std::vector<std::vector<double>>& error_matrix, int match, int mismatch, int gap_pen, int homo_gap_pen);

    };

}

#endif //GEFAST_CLUSTERREFINERS_HPP
