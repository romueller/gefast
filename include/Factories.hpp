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

#ifndef GEFAST_FACTORIES_HPP
#define GEFAST_FACTORIES_HPP

#include <map>

#include "Base.hpp"
#include "modes/AlignmentScoreMode.hpp"
#include "modes/LevenshteinMode.hpp"
#include "modes/QualityAlignmentScoreMode.hpp"
#include "modes/QualityLevenshteinMode.hpp"

namespace GeFaST {

    /*
     * General note:
     * All create(...) methods return nullptr when an invalid option is provided
     * or when the factory could not build the requested object.
     */


    /*
     * Factory class for creating Configuration instances.
     *
     * The type of the configuration depends on the mode in which the user wants to run GeFaST.
     */
    class ConfigurationFactory {

    public:
        /*
         * Create a configuration based on the given command-line arguments.
         * The first argument has to be the name of the mode in which to run GeFaST.
         */
        static Configuration* create(int argc, const char* argv[]);

    };


    /*
     * Factory class for creating AmpliconStorage instances.
     *
     * The type of the amplicon storage depends on the provided AmpliconStorageOption (see Options.hpp).
     */
    class AmpliconStorageFactory {

    public:
        /*
         * Create an amplicon storage which is empty but prepared for the insertion of amplicons corresponding
         * to the information given by the DataStatistics instance.
         */
        static AmpliconStorage* create(const AmpliconStorageOption opt, const DataStatistics<>& ds, const Configuration& config);

    };


    /*
     * Factory class for creating Distance instances.
     *
     * The type of the distance function depends on the provided DistanceOption (see Options.hpp) and,
     * to ensure the compatibility in some cases, the type of the configuration.
     */
    class DistanceFactory {

    public:
        /*
         * Create a distance function for computing the (bounded) distance between amplicons.
         */
        static Distance* create(const DistanceOption opt, const AmpliconStorage& amplicon_storage, const dist_t threshold,
                const Configuration& config);

        /*
         * Create a distance function for computing the (bounded) distance between amplicon as the number of edit operations between them.
         */
        static Distance* create(const DistanceOption opt, const AmpliconStorage& amplicon_storage, const dist_t threshold,
                const LevenshteinConfiguration& config);

        /*
         * Create a distance function for computing the (bounded) distance between amplicons as the score of an optimal alignment between them.
         */
        static Distance* create(const DistanceOption opt, const AmpliconStorage& amplicon_storage, const dist_t threshold,
                const AlignmentScoreConfiguration& config);

        /*
         * Create a distance function for computing the (bounded) distance between amplicons as the quality-weighted score
         * of an optimal alignment between them.
         */
        static Distance* create(const DistanceOption opt, const AmpliconStorage& amplicon_storage, const dist_t threshold,
                const QualityAlignmentScoreConfiguration& config);

        /*
         * Create a distance function for computing the (bounded) distance between amplicons as the number of edit operations
         * in an optimal quality-weighted alignment between them.
         */
        static Distance* create(const DistanceOption opt, const AmpliconStorage& amplicon_storage, const dist_t threshold,
                const QualityLevenshteinConfiguration& config);

    };


    /*
     * Factory class for creating Preprocessor instances.
     *
     * The type of the preprocessor depends on the provided PreprocessorOption (see Options.hpp).
     */
    class PreprocessorFactory {

    public:
        /*
         * Create and configure a preprocessor, which uses the given QualityEncoding instance.
         */
        static Preprocessor* create(const PreprocessorOption opt, const QualityEncoding<>& qe, const Configuration& config);

        /*
         * Create and configure a preprocessor, which will be equipped with a new QualityEncoding instance.
         */
        static Preprocessor* create(const PreprocessorOption opt_prep, const QualityEncodingOption opt_qe, const Configuration& config);

    };


    /*
     * Factory class for creating Clusterer instances.
     *
     * The type of the clusterer depends on the provided ClustererOption (see Options.hpp).
     */
    class ClustererFactory {

    public:
        /*
         * Create a clusterer governing the clustering phase, i.e. how the amplicons
         * are grouped into swarms (clusters) initially.
         */
        static Clusterer* create(const ClustererOption opt, const Configuration& config);

    };


    /*
     * Factory class for creating ClusterRefiner instances.
     *
     * The type of the cluster refiner depends on the provided ClusterRefinerOption (see Options.hpp).
     */
    class ClusterRefinerFactory {

    public:
        /*
         * Create a cluster refiner governing the refinement phase,
         * allowing for adjustments of the swarms found in the clustering phase.
         */
        static ClusterRefiner* create(const ClusterRefinerOption opt, const Configuration& config);

    };


    /*
     * Factory class for creating OutputGenerator instances.
     *
     * The type of the output generator depends on the provided OutputGeneratorOption (see Options.hpp).
     */
    class OutputGeneratorFactory {

    public:
        /*
         * Create a output generator governing the content of the generated outputs.
         */
        static OutputGenerator* create(const OutputGeneratorOption opt, const Configuration& config);

    };


    /*
     * Factory class for creating SwarmStorage instances.
     *
     * The type of the swarm storage depends on the provided SwarmStorageOption (see Options.hpp).
     */
    class SwarmStorageFactory {

    public:
        /*
         * Create a swarm storage which is empty but prepared for storing swarms for the given amplicon storage.
         */
        static SwarmStorage* create(const SwarmStorageOption opt, const AmpliconStorage& amplicon_storage, const Configuration& config);

    };


    /*
     * Factory class for creating AuxiliaryData instances.
     *
     * The type of the auxiliary data depends on the provided AuxiliaryDataOption resp. RefinementAuxiliaryDataOption (see Options.hpp),
     * and, thus, phase in which the auxiliary data is used.
     * Groups of types that require additional parameters (which are only available from some configuration types) might be generated
     * in separate create(...) methods.
     */
    class AuxiliaryDataFactory {

    public:
        /*
         * Create the auxiliary data for a specific amplicon pool and threshold suitable for the clustering phase.
         */
        static AuxiliaryData* create(const AuxiliaryDataOption opt, const AmpliconStorage& amplicon_storage,
                const numSeqs_t pool_id, const dist_t threshold, const Configuration& config);

        /*
         * Create the auxiliary data employing a segment filter for a specific amplicon pool and threshold
         * suitable for the clustering phase.
         */
        static AuxiliaryData* create(const AuxiliaryDataOption opt, const AmpliconStorage& amplicon_storage,
                const numSeqs_t pool_id, const dist_t threshold, const lenSeqs_t num_extra_segments, const Configuration& config);

        /*
         * Create the auxiliary data for a specific amplicon pool and threshold suitable for the refinement phase.
         */
        static AuxiliaryData* create(const RefinementAuxiliaryDataOption opt, const AmpliconStorage& amplicon_storage,
                const SwarmStorage& swarm_storage, const numSeqs_t pool_id, const dist_t threshold, const Configuration& config);

        /*
         * Create the auxiliary data employing a segment filter for a specific amplicon pool and threshold
         * suitable for the refinement phase.
         */
        static AuxiliaryData* create(const RefinementAuxiliaryDataOption opt, const AmpliconStorage& amplicon_storage,
                const SwarmStorage& swarm_storage, const numSeqs_t pool_id, const dist_t threshold,
                const lenSeqs_t num_extra_segments, const Configuration& config);

    };

}

#endif //GEFAST_FACTORIES_HPP
