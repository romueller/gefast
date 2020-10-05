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

#include "../../include/modes/QualityLevenshteinMode.hpp"

#include "../../include/Factories.hpp"
#include "../../include/Distances.hpp"

namespace GeFaST {

    QualityLevenshteinConfiguration::QualityLevenshteinConfiguration(int argc, const char* argv[]) {

        set_general_parameters(argc, argv);

        /* Set default values */

        if (main_threshold == 0) main_threshold = 1;
        num_extra_segments = 1;
        two_way_segment_filter = false;
        use_qgrams = false;

        opt_clusterer = (opt_clusterer == CL_IDLE) ? CL_CLASSIC_SWARMER : opt_clusterer;
        opt_refiner = (opt_refiner == CR_IDLE) ? CR_CLASSIC_FASTIDIOUS : opt_refiner;
        opt_output_generator = OG_CLASSIC;

        opt_distance = DT_QW_FRITH;
        unweighted_matches = false;
        inner_boost = true;
        opt_boosting = BO_NO_EFFECT;


        /* Read mode-specific configuration from file */

        std::ifstream in_stream(configuration_file);
        std::string line;
        unsigned long delim_pos;

        while (std::getline(in_stream, line).good()) {

            if (line.empty() || line.front() == '#') continue;

            delim_pos = line.find('=');
            std::string key = line.substr(0, delim_pos);

            if (key == "num_extra_segments") num_extra_segments = std::stoul(line.substr(delim_pos + 1));
            if (key == "two_way_segment_filter") two_way_segment_filter = (line.substr(delim_pos + 1) == "1");
            if (key == "use_qgrams") use_qgrams = (line.substr(delim_pos + 1) == "1");

            if (key == "distance") opt_distance = get_distance_option(line.substr(delim_pos + 1));

        }

        // quality-weighting and boosting parameters from misc
        auto iter = misc.begin();
        if ((iter = misc.find("unweighted_matches")) != misc.end()) {
            unweighted_matches = iter->second == "true";
        }
        if ((iter = misc.find("inner_boost")) != misc.end()) {
            inner_boost = iter->second == "true";
        }
        if ((iter = misc.find("boosting_method")) != misc.end()) {
            opt_boosting = get_boosting_option(iter->second);
        }

        qe_map = QualityEncoding<>(opt_quality_encoding);

        check();


        /* mode-specific options */

        // opt_preprocessor: user-defined, handled in set_general_parameters()
        // opt_quality_encoding: user-defined, handled in set_general_parameters()
        // opt_clusterer: user-defined, handled above and in set_general_parameters()
        // opt_refiner: user-defined, handled above and in set_general_parameters()
        // opt_output_generator: user-defined, handled above and in set_general_parameters()
        // opt_distance: user-defined, handled above
        // opt_boosting: user-defined, handled above
        // opt_amplicon_collection: given by opt_amplicon_storage

        opt_amplicon_storage = AS_PREPARED_LENGTH_POOLS_QUALITY;
        opt_swarm_storage = SS_SIMPLE_PER_POOL;
        opt_auxiliary_data = (two_way_segment_filter) ? AD_2W_SEGMENT_FILTER : AD_SEGMENT_FILTER;
        opt_refinement_auxiliary_data = (two_way_segment_filter) ? RD_2W_SEGMENT_FILTER : RD_SEGMENT_FILTER;

        if (use_qgrams) {
            opt_amplicon_storage = AS_PREPARED_QGRAM_LENGTH_POOLS_QUALITY;
        }

    }

    QualityLevenshteinConfiguration* QualityLevenshteinConfiguration::clone() const {
        return new QualityLevenshteinConfiguration(*this);
    }

    Preprocessor* QualityLevenshteinConfiguration::build_preprocessor() const {
        return PreprocessorFactory::create(opt_preprocessor, qe_map, *this);
    }

    QualityEncoding<>* QualityLevenshteinConfiguration::build_quality_encoding() const {
        return new QualityEncoding<>(opt_quality_encoding);
    }

    Clusterer* QualityLevenshteinConfiguration::build_clusterer() const {
        return ClustererFactory::create(opt_clusterer, *this);
    }

    ClusterRefiner* QualityLevenshteinConfiguration::build_cluster_refiner() const {
        return ClusterRefinerFactory::create(opt_refiner, *this);
    }

    OutputGenerator* QualityLevenshteinConfiguration::build_output_generator() const {
        return OutputGeneratorFactory::create(opt_output_generator, *this);
    }

    AmpliconStorage* QualityLevenshteinConfiguration::build_amplicon_storage(const DataStatistics<>& ds) const {
        return AmpliconStorageFactory::create(opt_amplicon_storage, ds, *this);
    }

    SwarmStorage* QualityLevenshteinConfiguration::build_swarm_storage(const AmpliconStorage& amplicon_storage) const {
        return SwarmStorageFactory::create(opt_swarm_storage, amplicon_storage, *this);
    }

    AuxiliaryData* QualityLevenshteinConfiguration::build_auxiliary_data(const AmpliconStorage& amplicon_storage,
            const numSeqs_t pool_id, const dist_t threshold) const {
        return AuxiliaryDataFactory::create(opt_auxiliary_data, amplicon_storage, pool_id, threshold, num_extra_segments, *this);
    }

    AuxiliaryData* QualityLevenshteinConfiguration::build_auxiliary_data(const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, const numSeqs_t pool_id, const dist_t threshold) const {

        return AuxiliaryDataFactory::create(opt_refinement_auxiliary_data, amplicon_storage, swarm_storage,
                pool_id, threshold, num_extra_segments, *this);

    }

    Distance* QualityLevenshteinConfiguration::build_distance_function(const AmpliconStorage &amplicon_storage,
                                                                       const dist_t threshold) const {
        return DistanceFactory::create(opt_distance, amplicon_storage, threshold, *this);
    }


    void QualityLevenshteinConfiguration::print(std::ostream& stream) const {

        Configuration::print(stream);

        ScoringFunction sf(match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty);

        stream << ">>> Mode-specific" << std::endl;
        stream << "Extra segments: " << num_extra_segments << std::endl;
        stream << "Segment filter: " << ((two_way_segment_filter) ? "two-way" : "one-way") << std::endl;
        stream << "Match reward: " << match_reward << " (transformed: 0)" << std::endl;
        stream << "Mismatch penalty: " << mismatch_penalty << " (transformed: " << sf.penalty_mismatch << ")" << std::endl;
        stream << "Gap-opening penalty: " << gap_opening_penalty << " (transformed: " << sf.penalty_open << ")" << std::endl;
        stream << "Gap-extension penalty: " << gap_extension_penalty << " (transformed: " << sf.penalty_extend << ")" << std::endl;
        stream << "Quality encoding: " << get_quality_encoding_label(opt_quality_encoding) << std::endl;
        if (use_qgrams) stream << "Prefiltering via q-grams: activated" << std::endl;

        stream << std::endl;

    }


    void QualityLevenshteinConfiguration::check() const {

        std::string mode = get_configuration_label(CF_QUALITY_LEVENSHTEIN);
        std::set<PreprocessorOption> valid_pp = {PP_FASTQ};
        std::set<ClustererOption > valid_cl = {CL_CLASSIC_SWARMER, CL_CONS_CLASSIC};
        std::set<ClusterRefinerOption > valid_cr = {CR_CLASSIC_FASTIDIOUS, CR_ITERATIVE_FASTIDIOUS, CR_CONS_LSA,
                                                    CR_CONS_LSR, CR_CONS_LSS, CR_CONS_SS};
        std::set<OutputGeneratorOption> valid_og = {OG_CLASSIC};
        std::set<QualityEncodingOption> valid_qe = {QE_SANGER, QE_ILLUMINA13, QE_ILLUMINA15, QE_ILLUMINA18};
        std::set<DistanceOption> valid_dt = {DT_QW_CLEMENT, DT_QW_MALDE_A, DT_QW_MALDE_B, DT_QW_MALDE_C, DT_QW_FRITH,
                                             DT_QW_KIM_A, DT_QW_KIM_B, DT_QW_CONVERGE_A, DT_QW_CONVERGE_B};
        std::set<BoostingOption> valid_bo = {BO_NO_EFFECT, BO_LINEAR, BO_MULT, BO_ROOT};

        if (valid_pp.find(opt_preprocessor) == valid_pp.end()){
            std::cerr << "ERROR: Mode '" << mode << "' supports only the FASTQ preprocessor." << std::endl;
            abort();
        }

        if (valid_cl.find(opt_clusterer) == valid_cl.end()) {
            std::cerr << "ERROR: Invalid 'clusterer' option for mode '" << mode << "'." << std::endl;
            abort();
        }

        if (valid_cr.find(opt_refiner) == valid_cr.end()) {
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

        if (valid_dt.find(opt_distance) == valid_dt.end()) {
            std::cerr << "ERROR: Invalid 'distance' option for mode '" << mode << "'." << std::endl;
            abort();
        }

        if (valid_bo.find(opt_boosting) == valid_bo.end()) {
            std::cerr << "ERROR: Invalid 'boosting_method' option for mode '" << mode << "'." << std::endl;
            abort();
        }

    }

}