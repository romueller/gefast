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

#include <fstream>

#include "../../include/modes/ConsistencyMode.hpp"

#include "../../include/Factories.hpp"
#include "../../include/Distances.hpp"

namespace GeFaST {

    ConsistencyConfiguration::ConsistencyConfiguration(int argc, const char* argv[]) {

        set_general_parameters(argc, argv);

        /* Set default values and ensure a proper configuration */

        main_threshold = 0.0;

        opt_clusterer = (opt_clusterer == CL_IDLE) ? CL_CONS_SWARMER : opt_clusterer;
        opt_refiner = (opt_refiner == CR_IDLE) ? CR_CONS_LSA : opt_refiner;
        opt_output_generator = OG_CLASSIC;


        /* Read mode-specific configuration from file */

        std::ifstream in_stream(configuration_file);
        std::string line;
        unsigned long delim_pos;

        while (std::getline(in_stream, line).good()) {

            if (line.empty() || line.front() == '#') continue;

            delim_pos = line.find('=');
            std::string key = line.substr(0, delim_pos);

            // currently nothing to do

        }

        check();

        qe_map = QualityEncoding<>(opt_quality_encoding);


        /* mode-specific options */

        // opt_preprocessor: user-defined, handled in set_general_parameters()
        // opt_quality_encoding: user-defined, handled in set_general_parameters()
        // opt_clusterer: user-defined, handled above and in set_general_parameters()
        // opt_refiner: user-defined, handled above and in set_general_parameters()
        // opt_output_generator: user-defined, handled above and in set_general_parameters()
        // opt_amplicon_collection: given by opt_amplicon_storage

        opt_swarm_storage = SS_SIMPLE_PER_POOL;
        opt_amplicon_storage = AS_PREPARED_POOLS_QUALITY;
        opt_auxiliary_data = AD_NAIVE_AUXILIARY;
        opt_refinement_auxiliary_data = RD_NAIVE_AUXILIARY;
        opt_distance = DT_UNKNOWN;

    }

    ConsistencyConfiguration* ConsistencyConfiguration::clone() const {
        return new ConsistencyConfiguration(*this);
    }

    Preprocessor* ConsistencyConfiguration::build_preprocessor() const {
        return PreprocessorFactory::create(opt_preprocessor, opt_quality_encoding, *this);
    }

    QualityEncoding<>* ConsistencyConfiguration::build_quality_encoding() const {
        return new QualityEncoding<>(opt_quality_encoding);
    }

    Clusterer* ConsistencyConfiguration::build_clusterer() const {
        return ClustererFactory::create(opt_clusterer, *this);
    }

    ClusterRefiner* ConsistencyConfiguration::build_cluster_refiner() const {
        return ClusterRefinerFactory::create(opt_refiner, *this);
    }

    OutputGenerator* ConsistencyConfiguration::build_output_generator() const {
        return OutputGeneratorFactory::create(opt_output_generator, *this);
    }


    AmpliconStorage* ConsistencyConfiguration::build_amplicon_storage(const DataStatistics<>& ds) const {
        return AmpliconStorageFactory::create(opt_amplicon_storage, ds, *this);
    }

    SwarmStorage* ConsistencyConfiguration::build_swarm_storage(const AmpliconStorage& amplicon_storage) const {
        return SwarmStorageFactory::create(opt_swarm_storage, amplicon_storage, *this);
    }

    AuxiliaryData* ConsistencyConfiguration::build_auxiliary_data(const AmpliconStorage& amplicon_storage,
            const numSeqs_t pool_id, const dist_t threshold) const {
        return AuxiliaryDataFactory::create(opt_auxiliary_data, amplicon_storage, pool_id, threshold, 1, *this);
    }

    AuxiliaryData* ConsistencyConfiguration::build_auxiliary_data(const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, const numSeqs_t pool_id, const dist_t threshold) const {

        return AuxiliaryDataFactory::create(opt_refinement_auxiliary_data, amplicon_storage, swarm_storage,
                pool_id, threshold, 1, *this);

    }

    Distance* ConsistencyConfiguration::build_distance_function(const AmpliconStorage &amplicon_storage,
                                                          const dist_t threshold) const {
        return DistanceFactory::create(opt_distance, amplicon_storage, threshold, *this);
    }


    void ConsistencyConfiguration::print(std::ostream& stream) const {

        Configuration::print(stream);

        ScoringFunction sf(match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty);

        stream << ">>> Mode-specific" << std::endl;
        stream << "Match reward: " << match_reward << " (transformed: 0)" << std::endl;
        stream << "Mismatch penalty: " << mismatch_penalty << " (transformed: " << sf.penalty_mismatch << ")" << std::endl;
        stream << "Gap-opening penalty: " << gap_opening_penalty << " (transformed: " << sf.penalty_open << ")" << std::endl;
        stream << "Gap-extension penalty: " << gap_extension_penalty << " (transformed: " << sf.penalty_extend << ")" << std::endl;
        stream << "Quality encoding: " << get_quality_encoding_label(opt_quality_encoding) << std::endl;

        stream << std::endl;

    }

    void ConsistencyConfiguration::check() const {

        std::string mode = get_configuration_label(CF_CONSISTENCY);
        std::set<PreprocessorOption> valid_pp = {PP_FASTQ};
        std::set<ClustererOption > valid_cl = {CL_CONS_SWARMER, CL_DADA2};
        std::set<ClusterRefinerOption > valid_cr = {CR_CONS_LSA, CR_CONS_LSR, CR_CONS_LSS, CR_CONS_SS};
        std::set<OutputGeneratorOption> valid_og = {OG_CLASSIC};
        std::set<QualityEncodingOption> valid_qe = {QE_SANGER, QE_ILLUMINA13, QE_ILLUMINA15, QE_ILLUMINA18};
//        std::set<DistanceOption> valid_dt = {};

        if (valid_pp.find(opt_preprocessor) == valid_pp.end()){
            std::cerr << "ERROR: Mode '" << mode << "' supports only the FASTQ preprocessor." << std::endl;
            abort();
        }

        if (valid_cl.find(opt_clusterer) == valid_cl.end()) {
            std::cerr << "ERROR: Invalid 'clusterer' option for mode '" << mode << "'." << std::endl;
            abort();
        }

        if (refinement_threshold != 0 && valid_cr.find(opt_refiner) == valid_cr.end()) {
            std::cerr << "ERROR: Invalid 'refiner' option for mode '" << mode << "'." << std::endl;
            abort();
        }

        if (valid_og.find(opt_output_generator) == valid_og.end()) {
            std::cerr << "ERROR: Invalid 'output_generator' option for mode '" << mode << "'." << std::endl;
            abort();
        }

        if (valid_qe.find(opt_quality_encoding) == valid_qe.end()) {
            std::cerr << "ERROR: Invalid 'quality_encoding' option for mode '" << mode << "'." << std::endl;
            abort();
        }

//        if (valid_dt.find(opt_distance) == valid_dt.end()) {
//            std::cerr << "WARNING: Mode '' ..." << std::endl;
//            std::cerr << "ERROR: Mode '' ..." << std::endl;
//            abort();
//        }

    }

}
