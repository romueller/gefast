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

#ifndef GEFAST_BASE_HPP
#define GEFAST_BASE_HPP

#include <algorithm>
#include <cstring>
#include <vector>
#include <set>
#include <emmintrin.h>

#include "Utility.hpp"

namespace GeFaST {

    /* === Base classes of main data structures === */

    class Preprocessor;
    class Clusterer;
    class ClusterRefiner;
    class OutputGenerator;

    class AmpliconStorage;
    class Swarm;
    class SwarmStorage;
    class Distance;
    class AuxiliaryData;


    /*
     * Abstract configuration class.
     *
     * Manages the mode-independent configuration parameters like input / output files and thresholds.
     *
     * The configuration is the central structure governing the execution of GeFaST and stores all the necessary
     * configuration parameters. In addition, a specific configuration also serves as an intermediary factory
     * for the key components (e.g. Preprocessor, Clusterer) that are configured appropriately
     * and injected into the clustering framework in order to run it in a specific mode.
     */
    struct Configuration {

        virtual ~Configuration() = default;

        virtual Configuration* clone() const = 0; // deep-copy clone method


        /*
         * Build specific preprocessor based on the configuration.
         */
        virtual Preprocessor* build_preprocessor() const = 0;

        /*
         * Build specific quality encoding based on the configuration.
         */
        virtual QualityEncoding<>* build_quality_encoding() const = 0;

        /*
         * Build specific clusterer (clustering phase) based on the configuration.
         */
        virtual Clusterer* build_clusterer() const = 0;

        /*
        * Build specific cluster refiner (refinement phase) based on the configuration.
        */
        virtual ClusterRefiner* build_cluster_refiner() const = 0;

        /*
         * Build specific output generator based on the configuration.
         */
        virtual OutputGenerator* build_output_generator() const = 0;


        /*
         * Build specific storage data structure for the amplicons based on the configuration and data statistics.
         */
        virtual AmpliconStorage* build_amplicon_storage(const DataStatistics<>& ds) const = 0;

        /*
         * Build specific storage data structure for the swarms.
         */
        virtual SwarmStorage* build_swarm_storage(const AmpliconStorage& amplicon_storage) const = 0;

        /*
         * Build specific distance function based on the configuration and, potentially, the amplicons.
         */
        virtual Distance* build_distance_function(const AmpliconStorage &amplicon_storage, const dist_t threshold) const = 0;

        /*
         * Build specific auxiliary data structures for a pool based on the configuration and the pool.
         */
        virtual AuxiliaryData* build_auxiliary_data(const AmpliconStorage& amplicon_storage, const numSeqs_t pool_id,
                const dist_t threshold) const = 0;

        /*
         * Build specific auxiliary data structures for a pool during refinement phase
         * based on the configuration and the pool.
         */
        virtual AuxiliaryData* build_auxiliary_data(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                const numSeqs_t pool_id, const dist_t threshold) const = 0;


        /*
         * Print configuration on given output stream.
         */
        virtual void print(std::ostream& stream) const;


        /* Configuration parameters */

        std::string version = VERSION; // version number of the program

        // input & output files
        std::vector<std::string> input_files; // names of all input files
        std::string configuration_file; // config file used to load (parts of) the configuration

        std::string output_internal; // name of the output file corresponding to Swarm's output option -i (internal structures)
        std::string output_otus; // name of the output file corresponding to Swarm's output option -o (OTUs / clusters / swarms)
        std::string output_statistics; // name of the output file corresponding to Swarm's output option -s (statistics file)
        std::string output_seeds; // name of the output file corresponding to Swarm's output option -w (seeds)
        std::string output_uclust; // name of the output file corresponding to Swarm's output option -u (uclust)

        std::string separator; // separator symbol (string) between ID and abundance in a FASTA header line
        bool mothur; // Boolean flag indicating demand for output format compatible with mothur, corresponds to Swarm's option -r

        char sep_internals = SEP_INTERNALS; // separator between columns in output_internal
        char sep_otus = SEP_OTUS; // separator between amplicon identifiers in output_otus
        char sep_statistics = SEP_STATISTICS; // separator between columns in output_statistics
        char sep_uclust = SEP_UCLUST; // separator between columns in output_uclust
        char sep_mothur = SEP_MOTHUR; // for mothur-compatible output
        char sep_mothur_otu = SEP_MOTHUR_OTU; // for mothur-compatible output

        // preprocessing
        bool filter_alphabet; // flag for the alphabet filter
        int filter_length; // flag for the length filter (0 = no filtering, 1 = only maximum threshold, 2 = only minimum threshold, 3 = both)
        int filter_abundance; // flag for the abundance filter (0 = no filtering, 1 = only maximum threshold, 2 = only minimum threshold, 3 = both)
        std::string alphabet; // allowed alphabet for the amplicon sequences
        lenSeqs_t min_length; // minimum allowed sequence length
        lenSeqs_t max_length; // maximum allowed sequence length
        numSeqs_t min_abundance; // minimum sequence abundance
        numSeqs_t max_abundance; // maximum sequence abundance
        std::string keep_preprocessed; // file to store preprocessed inputs

        // swarm clustering
        std::string mode; // mode of clustering
        dist_t main_threshold; // clustering threshold for the clustering phase
        dist_t refinement_threshold; // clustering threshold for the refinement phase
        std::vector<dist_t> iterative_refinement_thresholds; // clustering thresholds for the refinement phase (iterative cluster refinement only)
        numSeqs_t boundary; // minimum mass of a heavy cluster, used only in the refinement phase
        bool break_swarms; // Boolean flag indicating usage of swarm breaking, equivalent of Swarm's option -n

        // scoring function
        int match_reward; // reward for a nucleotide match
        int mismatch_penalty; // penalty for a nucleotide mismatch
        int gap_opening_penalty; // penalty for opening a gap
        int gap_extension_penalty; // penalty for extending a gap

        // miscellaneous (provided via a string-to-string map)
        std::map<std::string, std::string> misc;

    protected:
        Configuration() = default;

        Configuration(const Configuration& other) = default; // copy constructor

        Configuration(Configuration&& other) = default; // move constructor

        Configuration& operator=(const Configuration& other) = default; // copy assignment operator

        Configuration& operator=(Configuration&& other) = default; // move assignment operator

        /*
         * Parse the refinement thresholds given in list format (comma-separated list of thresholds)
         * or in range format (<first>:<last>:<increment>).
         * Sets the iterative refinement thresholds (iterative_refinement_thresholds)
         * and picks the last one as the single refinement threshold (refinement_threshold).
         */
        void parse_and_set_refinement_thresholds(std::string arg);

        /*
         * Parse the list of miscellaneous configuration parameters.
         * The $-separated list consists of key-value pairs of the format <key>:<value>.
         */
        void parse_misc(std::string arg);

        /*
         * Set default values, add basic configuration from file and apply command-line parameters.
         *
         * Assumed syntax in configuration file:
         *  - Line comments are allowed and start with #.
         *  - Every comment is written in its own line.
         *  - Empty lines are allowed.
         *  - Every configuration parameter is written in its own line.
         *  - A line containing a configuration parameter must have the following form: <key>=<value>
         *
         *  The command-line parameters have the highest priority and overwrite the configuration read from file.
         */
        void set_general_parameters(int argc, const char* argv[]);

        // data structure options
        PreprocessorOption opt_preprocessor; // preprocessor used in this run
        QualityEncodingOption opt_quality_encoding; // expected quality-score encoding (e.g. in FASTQ files)
        ClustererOption opt_clusterer; // method for clustering phase
        ClusterRefinerOption opt_refiner; // method for optional refinement phase
        OutputGeneratorOption opt_output_generator; // determines available outputs and their format
        AmpliconStorageOption opt_amplicon_storage; // representation of the amplicons
        AmpliconCollectionOption opt_amplicon_collection; // amplicon pool
        SwarmStorageOption opt_swarm_storage; // representation of swarms (clusters)
        AuxiliaryDataOption opt_auxiliary_data; // auxiliary data (structures) supporting the clustering phase
        RefinementAuxiliaryDataOption opt_refinement_auxiliary_data; // auxiliary data (structures) supporting the refinement phase
        DistanceOption opt_distance; // computation of distance between amplicons

    };


    /*
     * Abstract preprocessor class.
     *
     * Reads amplicons from a list of input files into an AmpliconStorage object.
     * Preprocesses the amplicons according to the given configuration.
     */
    class Preprocessor {

    public:
        virtual ~Preprocessor() = default;

        virtual Preprocessor* clone() const = 0; // deep-copy clone method


        /*
         * Perform the overall preprocessing step from input files to (filtered, preprocessed) amplicons.
         */
        virtual AmpliconStorage* preprocess_inputs(const std::vector<std::string>& files, const Configuration& conf) = 0;

    protected:
        Preprocessor() = default;

        Preprocessor(const Preprocessor& other) = default; // copy constructor

        Preprocessor(Preprocessor&& other) = default; // move constructor

        Preprocessor& operator=(const Preprocessor& other) = default; // copy assignment operator

        Preprocessor& operator=(Preprocessor&& other) = default; // move assignment operator

    };


    /*
     * Abstract class for handling the clustering phase.
     *
     * Performs the initial clustering of the given amplicons into swarms (clusters).
     * Parameters like the used thresholds are obtained from the configuration.
     */
    class Clusterer {

    public:
        virtual ~Clusterer() = default;

        virtual Clusterer* clone() const = 0; // deep-copy clone method


        /*
         * Perform the overall clustering phase.
         */
        virtual SwarmStorage* cluster(const AmpliconStorage& amplicon_storage, const Configuration& config) = 0;

    protected:
        Clusterer() = default;

        Clusterer(const Clusterer& other) = default; // copy constructor

        Clusterer(Clusterer&& other) = default; // move constructor

        Clusterer& operator=(const Clusterer& other) = default; // copy assignment operator

        Clusterer& operator=(Clusterer&& other) = default; // move assignment operator

    };


    /*
     * Abstract class for handling the optional refinement phase.
     *
     * Refines the given swarms according to the configuration.
     */
    class ClusterRefiner {

    public:
        virtual ~ClusterRefiner() = default;

        virtual ClusterRefiner* clone() const = 0; // deep-copy clone method


        /*
         * Perform the overall refinement phase.
         */
        virtual void refine(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage,
                const Configuration& config) = 0;

    protected:
        ClusterRefiner() = default;

        ClusterRefiner(const ClusterRefiner& other) = default; // copy constructor

        ClusterRefiner(ClusterRefiner&& other) = default; // move constructor

        ClusterRefiner& operator=(const ClusterRefiner& other) = default; // copy assignment operator

        ClusterRefiner& operator=(ClusterRefiner&& other) = default; // move assignment operator

    };


    /*
     * Abstract class for handling the output generation.
     *
     * Determines the available output files and their format.
     */
    class OutputGenerator {

    public:
        virtual ~OutputGenerator() = default;

        virtual OutputGenerator* clone() const = 0; // deep-copy clone method


        /*
         * Perform the overall output generation according to the configuration.
         */
        virtual void output_results(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                const Configuration& config) = 0;

    protected:
        OutputGenerator() = default;

        OutputGenerator(const OutputGenerator& other) = default; // copy constructor

        OutputGenerator(OutputGenerator&& other) = default; // move constructor

        OutputGenerator& operator=(const OutputGenerator& other) = default; // copy assignment operator

        OutputGenerator& operator=(OutputGenerator&& other) = default; // move assignment operator

    };

#if 1
    // Parts of the code related to q-grams adapted from:
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

        https://github.com/torognes/swarm/
    */
    #define cpuid(f1, f2, a, b, c, d)                                       \
    __asm__ __volatile__ ("cpuid"                                         \
                        : "=a" (a), "=b" (b), "=c" (c), "=d" (d)        \
                        : "a" (f1), "c" (f2));

    #define popcnt_asm(x, y)                                         \
    __asm__ __volatile__ ("popcnt %1,%0" : "=r"(y) : "r"(x));

    // mapping of nucleotides onto integers (a/A -> 1, c/C -> 2, g/G -> 3, t/T/u/U -> 4)
    extern char acgtu_map[256];

    extern long popcnt_present;

    extern long QGRAM_LENGTH;
    extern long QGRAM_VECTOR_BYTES;

    /*
     * See also the following borrowed methods in class QgramComparer:
     *  unsigned long popcount(unsigned long x)
     *  void cpu_features_detect()
     *  unsigned long popcount_128(__m128i x)
     *  unsigned long compare_qgram_vectors_128(const unsigned char * a, const unsigned char * b)
     *  unsigned long compare_qgram_vectors_64(const unsigned char * a, const unsigned char * b)
     *  unsigned long compare_qgram_vectors(const unsigned char * a, const unsigned char * b)
     */

    /*
     * Abstract class for storing a collection of amplicons with q-gram profiles.
     *
     * Provides additional methods for handling the q-gram profiles.
     */
    struct QgramComparer {

        virtual ~QgramComparer() = default;

        /*
         * Compute lower bound on the q-gram distance between two amplicons.
         *
         * Based on the ideas from "Approximate string-matching with q-grams and maximal matches"
         * (Ukkonen, 1992, https://core.ac.uk/download/pdf/82613520.pdf),
         * as implemented in Swarm (https://github.com/torognes/swarm).
         *
         * Only considers whether a q-gram has an even (including 0) or odd number of occurrences,
         * instead of using actual counts of the q-grams.
         */
        virtual unsigned long qgram_diff(const numSeqs_t i, const numSeqs_t j) const = 0;

    protected:
        QgramComparer() = default;

        QgramComparer(const QgramComparer& other) = default; // copy constructor

        QgramComparer(QgramComparer&& other) = default; // move constructor

        QgramComparer& operator=(const QgramComparer& other) = default; // copy assignment operator

        QgramComparer& operator=(QgramComparer&& other) = default; // move assignment operator


        /*
         * POPCOUNT- & SIMD-based methods for quick comparison of q-gram profiles.
         *
         * Support the computation of the lower bound on the q-gram distance.
         *
         * Borrowed from Swarm (https://github.com/torognes/swarm).
         */

        inline unsigned long popcount(unsigned long x) const;

        void cpu_features_detect();

        unsigned long popcount_128(__m128i x) const;

        unsigned long compare_qgram_vectors_128(const unsigned char* a, const unsigned char* b) const;

        unsigned long compare_qgram_vectors_64(const unsigned char* a, const unsigned char* b) const;

        unsigned long compare_qgram_vectors(const unsigned char* a, const unsigned char* b) const;

    };
#endif

    /*
     * Abstract class for storing a collection of amplicons (optionally including quality scores).
     *
     * Provides access to (the characteristics of) the amplicons.
     * While working on a single amplicon collection, the amplicons are accessed to by their position
     * in the collection (also referred to as "pool-internal integer id").
     */
    struct AmpliconCollection : public QgramComparer {

        virtual ~AmpliconCollection() = default;

        virtual AmpliconCollection* clone() const = 0; // deep-copy clone method


        /*
         * Append amplicon to the collection.
         *
         * Identifier and abundance are obtained from the Defline instance.
         * Additional information (specific to the derived class) can be obtained from the ExtraInfo instance.
         */
        virtual void emplace_back(const Defline& dl, const std::string& seq, const ExtraInfo& extra) = 0;

        /*
         * Return the number of amplicons in the collection.
         */
        virtual numSeqs_t size() const = 0;

        /*
         * Return the respective characteristic of the i-th amplicon of the collection.
         */
        virtual const char* id(const numSeqs_t i) const = 0; // identifier
        virtual const char* seq(const numSeqs_t i) const = 0; // sequence
        virtual lenSeqs_t len(const numSeqs_t i) const = 0; // length
        virtual numSeqs_t ab(const numSeqs_t i) const = 0; // abundance
        virtual const char* quals(const numSeqs_t i) const = 0; // returns nullptr when there are no quality scores available

        /*
         * Compute q-gram distance between two amplicons.
         */
        unsigned long qgram_diff(const numSeqs_t i, const numSeqs_t j) const override;

    protected:
        AmpliconCollection() = default;

        AmpliconCollection(const AmpliconCollection& other) = default; // copy constructor

        AmpliconCollection(AmpliconCollection&& other) = default; // move constructor

        AmpliconCollection& operator=(const AmpliconCollection& other) = default; // copy assignment operator

        AmpliconCollection& operator=(AmpliconCollection&& other) = default; // move assignment operator

    };


    /*
     * Abstract amplicon storage class.
     *
     * AmpliconStorage manages the amplicons to be clustered.
     * The amplicons are stored in one or more pools, which are expected to be disconnected,
     * i.e. there can be no links between amplicons from different pools (during all clustering phases)
     * based on the clustering (and refinement) threshold.
     *
     * The actual storage of the amplicons in the pool(s) is subject to the underlying AmpliconCollection class.
     */
    class AmpliconStorage {

    public:
        virtual ~AmpliconStorage() = default;

        virtual AmpliconStorage* clone() const = 0; // deep-copy clone method


        /*
         * Add an amplicon to the storage.
         *
         * ExtraInfo contains additional information (if any) that accompany the amplicon (e.g. FASTQ quality scores).
         */
        virtual void add_amplicon(const Defline& dl, const std::string& seq, const ExtraInfo& extra) = 0;

        /*
         * Retrieve i-th pool from the storage.
         */
        virtual AmpliconCollection& get_pool(const numSeqs_t i) = 0;

        /*
         * Retrieve i-th pool from the storage.
         */
        virtual const AmpliconCollection& get_pool(const numSeqs_t i) const = 0;

        /*
         * Return the number of pools in the storage.
         */
        virtual numSeqs_t num_pools() const = 0;

        /*
         * Return the total number of amplicons in all pools.
         */
        virtual numSeqs_t num_amplicons() const = 0;

        /*
         * Return the number of amplicons in the i-th pool.
         */
        virtual numSeqs_t size_pool(numSeqs_t i) const = 0;

        /*
         * Return the maximum sequence length of all amplicons in the storage.
         */
        virtual lenSeqs_t max_length() const = 0;

        /*
         * Finalise the storage, i.e. prepare it for use in the clustering phases.
         *
         * Should be called once after the last amplicon is added and before data is read from
         * the storage for the first time.
         */
        virtual void finalise() = 0;

        /*
         * Print the amplicons to the specified file.
         */
        virtual void print(const std::string& file, const Configuration& config) const = 0;

    protected:
        AmpliconStorage() = default;

        AmpliconStorage(const AmpliconStorage& other) = default; // copy constructor

        AmpliconStorage(AmpliconStorage&& other) = default; // move constructor

        AmpliconStorage& operator=(const AmpliconStorage& other) = default; // copy assignment operator

        AmpliconStorage& operator=(AmpliconStorage&& other) = default; // move assignment operator

    };


    /*
     * Abstract class representing a collection of possible grafting links
     * used during the refinement phase.
     *
     * For a grafting link, 'parent' refers to the amplicon / swarm to which
     * the 'child' swarm is added (via the 'child' amplicon).
     */
    class GraftingInfo {

    public:
        virtual ~GraftingInfo() = default;

        virtual GraftingInfo* clone() const = 0; // deep-copy clone method


        /*
         * Return pool-internal integer id of the parent amplicon of the i-th grafting link in this collection.
         */
        virtual numSeqs_t get_parent(const numSeqs_t i) const = 0;

        /*
         * Return pool-internal integer id of the child amplicon of the i-th grafting link in this collection.
         */
        virtual numSeqs_t get_child(const numSeqs_t i) const = 0;

        /*
         * Return a pointer to the swarm of the child amplicon of the i-th grafting link in this collection.
         */
        virtual Swarm* get_child_swarm(const numSeqs_t i) const = 0;

        /*
         * Return the distance between the parent and child amplicon of the i-th grafting link in this collection.
         */
        virtual dist_t get_distance(const numSeqs_t i) const = 0;

        /*
         * Add grafting link between parent amplicon p and child amplicon c (from swarm cs) with distance d.
         */
        virtual void add_link(const numSeqs_t p, const numSeqs_t c, Swarm* cs, const dist_t d) = 0;

        /*
         * Return the number of grafting links in this collection.
         */
        virtual numSeqs_t num_links() const = 0;

    protected:
        GraftingInfo() = default;

        GraftingInfo(const GraftingInfo& other) = default; // copy constructor

        GraftingInfo(GraftingInfo&& other) = default; // move constructor

        GraftingInfo& operator=(const GraftingInfo& other) = default; // copy assignment operator

        GraftingInfo& operator=(GraftingInfo&& other) = default; // move assignment operator

    };


    /*
     * Simple representation of a collection of grafting links.
     *
     * The different components of the grafting links are stored in separate STL vectors.
     */
    class SimpleGraftingInfo : public GraftingInfo {

    public:
        SimpleGraftingInfo() = default;

        virtual ~SimpleGraftingInfo() = default;

        SimpleGraftingInfo(const SimpleGraftingInfo& other) = default; // copy constructor

        SimpleGraftingInfo(SimpleGraftingInfo&& other) = default; // move constructor

        SimpleGraftingInfo& operator=(const SimpleGraftingInfo& other) = default; // copy assignment operator

        SimpleGraftingInfo& operator=(SimpleGraftingInfo&& other) = default; // move assignment operator

        SimpleGraftingInfo* clone() const override; // deep-copy clone method


        numSeqs_t get_parent(const numSeqs_t i) const override;

        numSeqs_t get_child(const numSeqs_t i) const override;

        Swarm* get_child_swarm(const numSeqs_t i) const override;

        dist_t get_distance(const numSeqs_t i) const override;

        void add_link(const numSeqs_t p, const numSeqs_t c, Swarm* cs, const dist_t d) override;

        numSeqs_t num_links() const override;

    protected:
        std::vector<numSeqs_t> parents_; // pool-internal integer id of parent amplicon
        std::vector<numSeqs_t> children_; // pool-internal integer id of child amplicon
        std::vector<Swarm*> child_swarms_; // swarm the child amplicon belongs to
        std::vector<dist_t> distances_; // distance between parent and child amplicon

    };


    /*
     * Abstract class representing a member (entry) of a swarm.
     *
     * Provides access to the pool-internal ids of the member, its parent in the swarm,
     * the distance to the parent and its generation and radius in the swarm.
     */
    class SwarmEntry {

    public:
        virtual ~SwarmEntry() = default;

        virtual SwarmEntry* clone() const = 0; // deep-copy clone method


        /*
         * Return pool-internal integer id of the member.
         */
        virtual numSeqs_t member() const = 0;

        /*
         * Return the pool-internal integer id of the parent of the member in the swarm.
         */
        virtual numSeqs_t parent() const = 0;

        /*
         * Return the generation of the member in the swarm.
         */
        virtual numSeqs_t gen() const = 0;

        /*
         * Return the radius of the member in the swarm,
         * i.e. the sum of the distances on the links between the member and the seed of the swarm.
         */
        virtual dist_t rad() const = 0;

        /*
         * Return the distance between the member and its parent in the swarm.
         */
        virtual dist_t parent_dist() const = 0;

    protected:
        SwarmEntry() = default;

        SwarmEntry(const SwarmEntry& other) = default; // copy constructor

        SwarmEntry(SwarmEntry&& other) = default; // move constructor

        SwarmEntry& operator=(const SwarmEntry& other) = default; // copy assignment operator

        SwarmEntry& operator=(SwarmEntry&& other) = default; // move assignment operator

    };


    /*
     * Abstract class representing a swarm (cluster).
     *
     * Provides access to its members and their 'place' in it (parent, generation etc.).
     */
    class Swarm {

    public:
        virtual ~Swarm() = default;

        virtual Swarm* clone() const = 0; // deep-copy clone method


        /*
         * Return pool-internal integer id of the seed of the swarm.
         */
        virtual numSeqs_t seed() const = 0;

        /*
         * Construct representation of the seed of the swarm.
         * The second method reuses the provided SwarmEntry instance instead of allocating an additional one.
         */
        virtual SwarmEntry* seed_entry() const = 0;
        virtual SwarmEntry* seed_entry(SwarmEntry* entry) const = 0;

        /*
         * Return the pool-internal integer id of the i-th member of the swarm.
         */
        virtual numSeqs_t get(const numSeqs_t i) = 0;

        /*
         * Construct representation of the i-th member of the swarm.
         * The second method reuses the provided SwarmEntry instance instead of allocating an additional one.
         */
        virtual SwarmEntry* get_entry(const numSeqs_t i) = 0;
        virtual SwarmEntry* get_entry(const numSeqs_t i, SwarmEntry* entry) const = 0;

        /*
         * Return the pool-internal integer id of the i-th member of the swarm.
         */
        virtual numSeqs_t member(const numSeqs_t i) const = 0;

        /*
         * Return the pool-internal integer id of the i-th member of the swarm.
         */
        virtual numSeqs_t parent(const numSeqs_t i) const = 0;

        /*
         * Return the generation of the i-th member of the swarm.
         */
        virtual numSeqs_t gen(const numSeqs_t i) const = 0;

        /*
         * Return the radius of the i-th member of the swarm.
         */
        virtual dist_t rad(const numSeqs_t i) const = 0;

        /*
         * Return the distance of the i-th member to its parent in the swarm.
         */
        virtual dist_t parent_dist(const numSeqs_t i) const = 0;

        /*
         * Return the number of members in the swarm.
         * Does not include amplicons from attached swarms (see total_size()).
         */
        virtual numSeqs_t size() const = 0;

        /*
         * Return the number of different sequences found in the swarm.
         * This number can differ from the size of the swarm when the amplicons
         * were not dereplicated before clustering them.
         * Does not include amplicons from attached swarms (see total_num_different()).
         */
        virtual numSeqs_t num_different() const = 0;

        /*
         * Return the number of singletons (i.e. amplicons with abundance one) in the swarm.
         * Does not include amplicons from attached swarms (see total_num_singletons()).
         */
        virtual numSeqs_t num_singletons() const = 0;

        /*
         * Return the mass of the swarm (i.e. the sum of the abundances of all amplicons).
         * Does not include amplicons from attached swarms (see total_mass()).
         */
        virtual numSeqs_t mass() const = 0;

        /*
         * Sort the generation of amplicons of the swarm to be explored next.
         * That generation starts with the p-th member and ends with the currently last member.
         */
        virtual void sort_next_generation(const AmpliconCollection& ac, const numSeqs_t p) = 0;

        /*
         * Append new member m (generation g, radius r) with parent p (both via their pool-internal integer ids),
         * distance dist to p and abundance ab to the swarm.
         * Flag new_seq indicates whether the sequence of the amplicon has already been encountered
         * (influences num_different()).
         */
        virtual void append(const numSeqs_t m, const numSeqs_t g, const dist_t r, const numSeqs_t p, const dist_t dist,
                const numSeqs_t ab, const bool new_seq) = 0;

        /*
         * Append new member m (pool-internal integer id) with parent p (p-th member of swarm),
         * distance dist to p and abundance ab to the swarm.
         * Flag new_seq indicates whether the sequence of the amplicon has already been encountered
         * (influences num_different()).
         */
        virtual void append(const numSeqs_t m, const numSeqs_t p, const dist_t dist, const numSeqs_t ab, const bool new_seq) = 0;

        /*
         * Attach swarm cs to this swarm via a link between amplicon p from this swarm and amplicon c from cs.
         * The distance between p and c is d.
         */
        virtual void attach(const numSeqs_t p, const numSeqs_t c, Swarm* cs, const dist_t d) = 0;

        /*
         * Indicate whether this swarm is attached to (i.e. is the child of) another swarm.
         */
        virtual bool is_attached() const = 0;

        /*
         * Mark this swarm as attached.
         */
        virtual void mark_as_attached() = 0;

        /*
         * Return the number of members in the swarm (including attached swarms).
         */
        virtual numSeqs_t total_size() const = 0;

        /*
         * Return the number of different sequences in the swarm (including attached swarms).
         */
        virtual numSeqs_t total_num_different() const = 0;

        /*
         * Return the number of singletons in the swarm (including attached swarms).
         */
        virtual numSeqs_t total_num_singletons() const = 0;

        /*
         * Return the mass of the swarm (including attached swarms).
         */
        virtual numSeqs_t total_mass() const = 0;

        /*
         * Return the maximum radius of the swarm.
         * Does not consider attached swarms.
         */
        virtual dist_t max_rad() const = 0;

        /*
         * Return the maximum generation of the swarm
         * (i.e. number of iterations before the swarm reached its natural limit).
         * Does not consider attached swarms.
         */
        virtual lenSeqs_t max_gen() const = 0;

        /*
         * Retrieve the grafting links with this swarm as the parent.
         */
        virtual GraftingInfo* get_grafting_info() const = 0;

    protected:
        Swarm() = default;

        Swarm(const Swarm& other) = default; // copy constructor

        Swarm(Swarm&& other) = default; // move constructor

        Swarm& operator=(const Swarm& other) = default; // copy assignment operator

        Swarm& operator=(Swarm&& other) = default; // move assignment operator

    };


    /*
     * Abstract class representing a collection of swarms.
     */
    class Swarms {

    public:
        virtual ~Swarms() = default;

        virtual Swarms* clone() const = 0; // deep-copy clone method


        /*
         * Add new swarm to the collection and initialise it with the seed
         * (pool-internal integer id s, abundance ab).
         */
        virtual Swarm& initialise_cluster(const numSeqs_t s, const numSeqs_t ab) = 0;

        /*
         * Return the number of swarms in the collection.
         */
        virtual numSeqs_t size() const  = 0;

        /*
         * Return the i-th swarm in the collection.
         */
        virtual Swarm& get(const numSeqs_t i) = 0;
        virtual const Swarm& get(const numSeqs_t i) const = 0;

        /*
         * Discard the swarms in the collection, leaving behind an empty collection.
         */
        virtual void clear() = 0;

    protected:
        Swarms() = default;

        Swarms(const Swarms& other) = default; // copy constructor

        Swarms(Swarms&& other) = default; // move constructor

        Swarms& operator=(const Swarms& other) = default; // copy assignment operator

        Swarms& operator=(Swarms&& other) = default; // move assignment operator

    };


    /*
     * Abstract swarm storage class.
     *
     * SwarmStorage stores and manages swarms during clustering, refinement and output generation.
     * The swarms are stored per pool.
     *
     * The actual storage of the swarms is subject to the underlying Swarms class.
     */
    class SwarmStorage {

    public:
        virtual ~SwarmStorage() = default;

        virtual SwarmStorage* clone() const = 0; // deep-copy clone method


        /*
         * Return the swarms of the p-th pool.
         */
        virtual Swarms& get_swarms(const numSeqs_t p) = 0;
        virtual const Swarms& get_swarms(const numSeqs_t p) const = 0;

        /*
         * Return the number of pools.
         */
        virtual numSeqs_t num_pools() const = 0;

        /*
         * Return the number of swarms in all pools.
         */
        virtual numSeqs_t num_swarms() const = 0;

        /*
         * Return the number of swarms in the p-th pool.
         */
        virtual numSeqs_t num_swarms(const numSeqs_t p) const = 0;

    protected:
        SwarmStorage() = default;

        SwarmStorage(const SwarmStorage& other) = default; // copy constructor

        SwarmStorage(SwarmStorage&& other) = default; // move constructor

        SwarmStorage& operator=(const SwarmStorage& other) = default; // copy assignment operator

        SwarmStorage& operator=(SwarmStorage&& other) = default; // move assignment operator

    };


    /*
     * Method provider for generating an order of swarms across pools.
     *
     * The returned order consists of pairs of pointers to a Swarm instance (second member)
     * and the AmpliconCollection (first member), from which it originates.
     */
    struct SwarmOrderer {

        /*
         * Sort the swarms by their mass (descending)
         * and use the lexicographical rank (ascending) of the seeds as the tie-breaker.
         */
        static std::vector<std::pair<const AmpliconCollection*, const Swarm*>> order_by_mass(
                const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage);

        /*
         * Sort the swarms by the abundance of their seeds (descending)
         * and use the lexicographical rank (ascending) of the seeds as the tie-breaker.
         */
        static std::vector<std::pair<const AmpliconCollection*, const Swarm*>> order_by_seed_abundance(
                const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage);

    };


    /*
     * Abstract distance-function class.
     *
     * Provides the method to compute the distance between two amplicons.
     */
    class Distance {

    public:
        virtual ~Distance() = default;

        virtual Distance* clone() const = 0; // deep-copy clone method


        /*
         * Compute distance between the two given amplicons.
         */
        virtual dist_t distance(const AmpliconCollection& ac, const numSeqs_t i, const numSeqs_t j) = 0;

    protected:
        Distance() = default;

        Distance(const Distance& other) = default; // copy constructor

        Distance(Distance&& other) = default; // move constructor

        Distance& operator=(const Distance& other) = default; // copy assignment operator

        Distance& operator=(Distance&& other) = default; // move assignment operator

    };


    /*
     * Abstract class describing auxiliary data for the clustering phases.
     *
     * Auxiliary data supports the efficient identification of partners of an amplicon.
     */
    class AuxiliaryData {

    public:
        /*
         * Representation of a partner of an amplicon providing the necessary information
         * to build a link in the swarm.
         */
        struct Partner {

            numSeqs_t id; // pool-internal integer id of the partner
            dist_t dist; // distance of the partner to the currently explored amplicon

            Partner(numSeqs_t i, dist_t d) : id(i), dist(d) {
                // nothing else to do
            }

        };

        virtual ~AuxiliaryData() = default;

        virtual AuxiliaryData* clone() const = 0; // deep-copy clone method


        /*
         * Mark the amplicon with the given pool-internal integer id as swarmed,
         * such that it is not considered as a partner of other amplicons as well.
         */
        virtual void tick_off_amplicon(const numSeqs_t ampl_id) = 0;

        /*
         * Record occurrence of an amplicon.
         * This is used, e.g., to determine whether the sequence of the amplicon is new when adding it to a swarm.
         * Returns true iff the sequence is new.
         */
        virtual bool record_amplicon(const numSeqs_t ampl_id) = 0;

        /*
         * Reset the records on already observed amplicons.
         */
        virtual void clear_amplicon_records() = 0;

        /*
         * Determine the partners of the amplicon with the given pool-internal integer id.
         */
        virtual std::vector<Partner> find_partners(const numSeqs_t ampl_id) = 0;

    protected:
        AuxiliaryData() = default;

        AuxiliaryData(const AuxiliaryData& other) = default; // copy constructor

        AuxiliaryData(AuxiliaryData&& other) = default; // move constructor

        AuxiliaryData& operator=(const AuxiliaryData& other) = default; // copy assignment operator

        AuxiliaryData& operator=(AuxiliaryData&& other) = default; // move assignment operator

    };


    /* === Components related to alignment-free methods === */

    /*
     * Abstract class for storing a collection of amplicons with further features.
     *
     * Provides additional access to these features.
     */
    struct FeatureAmpliconCollection : public AmpliconCollection {

        virtual ~FeatureAmpliconCollection() = default;

        FeatureAmpliconCollection* clone() const override = 0; // deep-copy clone method


        /*
         * Return the features of the i-th amplicon of the collection.
         */
        virtual const feat_t* features(const numSeqs_t i) const = 0;

        /*
         * Return the features of all amplicons of the collection.
         */
        virtual const feat_t* all_features() const = 0;

        /*
         * Return the number of features per amplicon.
         */
        virtual size_t num_features() const = 0;

    protected:
        FeatureAmpliconCollection() = default;

        FeatureAmpliconCollection(const FeatureAmpliconCollection& other) = default; // copy constructor

        FeatureAmpliconCollection(FeatureAmpliconCollection&& other) = default; // move constructor

        FeatureAmpliconCollection& operator=(const FeatureAmpliconCollection& other) = default; // copy assignment operator

        FeatureAmpliconCollection& operator=(FeatureAmpliconCollection&& other) = default; // move assignment operator

    };

    /*
     * Abstract amplicon storage class for amplicons with further features.
     *
     * Essentially a covariant version of AmpliconStorage.
     */
    class FeatureAmpliconStorage : public AmpliconStorage {

    public:
        virtual ~FeatureAmpliconStorage() = default;

        FeatureAmpliconStorage* clone() const override = 0; // deep-copy clone method

        /*
         * Retrieve i-th pool from the storage.
         */
        virtual FeatureAmpliconCollection& get_pool(const numSeqs_t i) override = 0;

        /*
         * Retrieve i-th pool from the storage.
         */
        virtual const FeatureAmpliconCollection& get_pool(const numSeqs_t i) const override = 0;

    protected:
        FeatureAmpliconStorage() = default;

        FeatureAmpliconStorage(const FeatureAmpliconStorage& other) = default; // copy constructor

        FeatureAmpliconStorage(FeatureAmpliconStorage&& other) = default; // move constructor

        FeatureAmpliconStorage& operator=(const FeatureAmpliconStorage& other) = default; // copy assignment operator

        FeatureAmpliconStorage& operator=(FeatureAmpliconStorage&& other) = default; // move assignment operator

    };


    /*
     * Abstract class for determining the feature representation (as a "vector") of an amplicon.
     */
    class FeatureBuilder {

    public:
        virtual ~FeatureBuilder() = default;

        virtual FeatureBuilder* clone() const = 0; // deep-copy clone method

        /*
         * Return the number of features per amplicon.
         */
        virtual size_t num_features() = 0;

        /*
         * Determine the feature representation of the provided amplicon (sequence).
         *
         * The last two versions avoid constructing a new "vector" for each call.
         */
        virtual std::vector<feat_t> get_features(const std::string& seq) = 0;
        virtual void get_features(const std::string& seq, std::vector<feat_t>& features) = 0;
        virtual void get_features(const std::string& seq, feat_t* features) = 0;

    };


    /*
     * Abstract distance-function class when the distance is based on feature representations.
     *
     * Essentially a covariant version of Distance.
     */
    class FeatureDistance : public Distance {

    public:
        FeatureDistance* clone() const override = 0;

        virtual dist_t distance(const FeatureAmpliconCollection& fac, const numSeqs_t i, const numSeqs_t j) = 0;

    };

}

#endif //GEFAST_BASE_HPP
