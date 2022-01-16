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

#ifndef GEFAST_SPACELEVENSHTEINMODE_HPP
#define GEFAST_SPACELEVENSHTEINMODE_HPP

#include "../Base.hpp"

namespace GeFaST {

    /*
     * Chunk lengths and other parameters used to configure the DACs sequences
     * used in the amplicon and swarm representations.
     */
    struct ChunkLengths {

        usmall_t ids; // chunk length of member and parent ids
        usmall_t dists; // chunk length of parent distances
        usmall_t gens; // chunk length of generation numbers
        usmall_t rads; // chunk length of radii

        usmall_t starts; // chunk length of swarm start positions
        usmall_t sizes; // chunk length of swarm sizes
        usmall_t masses; // chunk length of swarm masses
        usmall_t differents; // chunk length of number of different sequences in swarms
        usmall_t singletons; // chunk length of number of singletons in swarms
        usmall_t max_rads; // chunk length of maximum radii of swarms
        float init_rel_capa; // initial capacity of DACs sequences for the swarm information
        float level_factor; // factor for size of new level, i.e. ((size of level) = (size of previous level) * (level_factor))
        float extend_factor; // factor for extending a full level, i.e. ((new size) = (current size) * (extend_factor))

        bool opt; // flag indicating whether the chunk lengths and the number of levels should be optimised during finalisation

        usmall_t headers_chunk; // chunk length of cursors in header component of FlexibleAmpliconCollection
        float headers_level; // level factor of cursors in header component of FlexibleAmpliconCollection
        float headers_extend; // extend factor of cursors in header component of FlexibleAmpliconCollection
        bool headers_opt; // optimisation flag of cursors in header component of FlexibleAmpliconCollection

        usmall_t prefixes_chunk; // chunk length of prefix numbers in header component of FlexibleAmpliconCollection (PrefixDacsStringArray only)
        float prefixes_level; // level factor of prefix numbers in header component of FlexibleAmpliconCollection (PrefixDacsStringArray only)
        float prefixes_extend; // extend factor of prefix numbers in header component of FlexibleAmpliconCollection (PrefixDacsStringArray only)
        bool prefixes_opt; // optimisation flag of prefix numbers in header component of FlexibleAmpliconCollection (PrefixDacsStringArray only)

        usmall_t seqs_chunk; // chunk length of cursors in sequence component of FlexibleAmpliconCollection
        float seqs_level; // level factor of cursors in sequence component of FlexibleAmpliconCollection
        float seqs_extend; // extend factor of cursors in sequence component of FlexibleAmpliconCollection
        bool seqs_opt; // optimisation flag of cursors in sequence component of FlexibleAmpliconCollection

        usmall_t qgram_chunk; // chunk length of cursors in q-gram component of FlexibleAmpliconCollection
        float qgram_level; // level factor of cursors in q-gram component of FlexibleAmpliconCollection
        float qgram_extend; // extend factor of cursors in q-gram component of FlexibleAmpliconCollection
        bool qgram_opt; // optimisation flag of cursors in q-gram component of FlexibleAmpliconCollection

        /*
         * Convenience method for determining initial capacity using init_rel capa and the provided capacity.
         */
        numSeqs_t get_capa(numSeqs_t capa) const;

    };

    /*
     * Parameters of the k^2-trees used to represent binary relations.
     */
    struct K2TreeParameters {

        size_t k = 0; // arity of k^2-tree
        size_t kr = 0; // row arity of k^2-tree
        size_t kc = 0; // column arity of k^2-tree
        size_t uk = 0; // arity in the upper part of the k^2-tree
        size_t uh = 0; // height of the upper part of the k^2-tree
        size_t lk = 0; // arity in the lower part of the k^2-tree
        size_t mb = 0; // boundary for decision between NaiveK2Tree (<=) and RectangularK2Tree (>)
        bool static_rels = true; // flag indicating whether elements to be removed are only marked as such
        bool deep = true; // flag indicating whether set_null_column() is executed in deep or shallow form

    };


    /*
     * Configuration of special mode for investigating the trade-off between runtime and memory consumption
     * when using space-efficient data structures in GeFaST's Levenshtein mode (space).
     */
    struct SpaceLevenshteinConfiguration : public Configuration {

        SpaceLevenshteinConfiguration(int argc, const char* argv[]);

        SpaceLevenshteinConfiguration* clone() const override; // deep-copy clone method


        Preprocessor* build_preprocessor() const override;

        QualityEncoding<>* build_quality_encoding() const override;

        Clusterer* build_clusterer() const override;

        ClusterRefiner* build_cluster_refiner() const override;

        OutputGenerator* build_output_generator() const override;

        AmpliconStorage* build_amplicon_storage(const DataStatistics<>& ds) const override;

        SwarmStorage* build_swarm_storage(const AmpliconStorage& amplicon_storage) const override;

        AuxiliaryData* build_auxiliary_data(const AmpliconStorage& amplicon_storage, const numSeqs_t pool_id,
                const dist_t threshold) const override;

        AuxiliaryData* build_auxiliary_data(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                const numSeqs_t pool_id, const dist_t threshold) const override;

        Distance* build_distance_function(const AmpliconStorage &amplicon_storage, const dist_t threshold) const override;


        void print(std::ostream& stream) const override;

        /*
         * Check whether the configuration obeys the restrictions of the mode.
         */
        void check() const;


        /* Mode-specific configuration parameters */

        lenSeqs_t num_extra_segments; // parameter of the pigeonhole principle (segment filter)
        bool use_score; // flag indicating whether to use an actual scoring function (not the edit distance)
        bool use_qgrams; // flag indicating whether alignment computations are preceded with a q-gram-based filtering step

        unsigned char prefix_len; // length of prefixes replaced by numbers (HuffmanStringArray)

        ChunkLengths chunk_lengths; // parameters for DACs

        K2TreeParameters tree_parameters; // parameters for k^2-trees used in (Static)SpaceSegmentRelations

    };

}

#endif //GEFAST_SPACELEVENSHTEINMODE_HPP
