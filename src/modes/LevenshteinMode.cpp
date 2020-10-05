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

#include "../../include/modes/LevenshteinMode.hpp"

#include "../../include/Factories.hpp"
#include "../../include/Distances.hpp"

namespace GeFaST {

    LevenshteinConfiguration::LevenshteinConfiguration(int argc, const char* argv[]) {

        set_general_parameters(argc, argv);

        /* Set default values */

        if (main_threshold == 0) main_threshold = 1;
        num_extra_segments = 1;
        two_way_segment_filter = false;
        use_score = false;
        use_qgrams = false;

        opt_clusterer = (opt_clusterer == CL_IDLE) ? CL_CLASSIC_SWARMER : opt_clusterer;
        opt_refiner = (opt_refiner == CR_IDLE) ? CR_CLASSIC_FASTIDIOUS : opt_refiner;
        opt_output_generator = OG_CLASSIC;


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
            if (key == "use_score") use_score = (line.substr(delim_pos + 1) == "1");
            if (key == "use_qgrams") use_qgrams = (line.substr(delim_pos + 1) == "1");

        }

        check();


        /* mode-specific data structures and settings */

        // opt_preprocessor: user-defined, handled in set_general_parameters()
        // opt_quality_encoding: user-defined, handled in set_general_parameters()
        // opt_clusterer: user-defined, handled above and in set_general_parameters()
        // opt_refiner: user-defined, handled above and in set_general_parameters()
        // opt_output_generator: user-defined, handled above and in set_general_parameters()
        // opt_amplicon_collection: given by opt_amplicon_storage

        opt_amplicon_storage = AS_PREPARED_LENGTH_POOLS;
        opt_swarm_storage = SS_SIMPLE_PER_POOL;
        opt_auxiliary_data = (two_way_segment_filter) ? AD_2W_SEGMENT_FILTER : AD_SEGMENT_FILTER;
        opt_refinement_auxiliary_data = (two_way_segment_filter) ? RD_2W_SEGMENT_FILTER : RD_SEGMENT_FILTER;
        opt_distance = (use_score) ? DT_BOUNDED_SCORE_LEVENSHTEIN : DT_BOUNDED_LEVENSHTEIN;

        bool qual_required = false;
        if (opt_refiner == CR_CONS_LSA || opt_refiner == CR_CONS_LSR || opt_refiner == CR_CONS_LSS ||
            opt_refiner == CR_CONS_SS || opt_clusterer == CL_CONS_CLASSIC) {
            qual_required = true;
        }

        if (use_qgrams) {
            opt_amplicon_storage = AS_PREPARED_QGRAM_LENGTH_POOLS;
        }

        if (qual_required) {

            if (opt_amplicon_storage == AS_PREPARED_POOLS) opt_amplicon_storage = AS_PREPARED_POOLS_QUALITY;
            if (opt_amplicon_storage == AS_PREPARED_LENGTH_POOLS) opt_amplicon_storage = AS_PREPARED_LENGTH_POOLS_QUALITY;
            if (opt_amplicon_storage == AS_PREPARED_QGRAM_LENGTH_POOLS) opt_amplicon_storage = AS_PREPARED_QGRAM_LENGTH_POOLS_QUALITY;

        }

    }

    LevenshteinConfiguration* LevenshteinConfiguration::clone() const {
        return new LevenshteinConfiguration(*this);
    }

    Preprocessor* LevenshteinConfiguration::build_preprocessor() const {
        return PreprocessorFactory::create(opt_preprocessor, opt_quality_encoding, *this);
    }

    QualityEncoding<>* LevenshteinConfiguration::build_quality_encoding() const {
        return new QualityEncoding<>(opt_quality_encoding);
    }

    Clusterer* LevenshteinConfiguration::build_clusterer() const {
        return ClustererFactory::create(opt_clusterer, *this);
    }

    ClusterRefiner* LevenshteinConfiguration::build_cluster_refiner() const {
        return ClusterRefinerFactory::create(opt_refiner, *this);
    }

    OutputGenerator* LevenshteinConfiguration::build_output_generator() const {
        return OutputGeneratorFactory::create(opt_output_generator, *this);
    }

    AmpliconStorage* LevenshteinConfiguration::build_amplicon_storage(const DataStatistics<>& ds) const {
        return AmpliconStorageFactory::create(opt_amplicon_storage, ds, *this);
    }

    SwarmStorage* LevenshteinConfiguration::build_swarm_storage(const AmpliconStorage& amplicon_storage) const {
        return SwarmStorageFactory::create(opt_swarm_storage, amplicon_storage, *this);
    }

    AuxiliaryData* LevenshteinConfiguration::build_auxiliary_data(const AmpliconStorage& amplicon_storage,
            const numSeqs_t pool_id, const dist_t threshold) const {
        return AuxiliaryDataFactory::create(opt_auxiliary_data, amplicon_storage, pool_id, threshold, num_extra_segments, *this);
    }

    AuxiliaryData* LevenshteinConfiguration::build_auxiliary_data(const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, const numSeqs_t pool_id, const dist_t threshold) const {

        return AuxiliaryDataFactory::create(opt_refinement_auxiliary_data, amplicon_storage, swarm_storage,
                pool_id, threshold, num_extra_segments, *this);

    }

    Distance* LevenshteinConfiguration::build_distance_function(const AmpliconStorage &amplicon_storage,
                                                                const dist_t threshold) const {
        return DistanceFactory::create(opt_distance, amplicon_storage, threshold, *this);
    }


    void LevenshteinConfiguration::print(std::ostream& stream) const {

        Configuration::print(stream);

        ScoringFunction sf(match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty);

        stream << ">>> Mode-specific" << std::endl;
        stream << "Extra segments: " << num_extra_segments << std::endl;
        stream << "Segment filter: " << ((two_way_segment_filter) ? "two-way" : "one-way") << std::endl;
        stream << "Use scoring function: " << use_score << std::endl;
        if (use_score) {

            stream << "Match reward: " << match_reward << " (transformed: 0)" << std::endl;
            stream << "Mismatch penalty: " << mismatch_penalty << " (transformed: " << sf.penalty_mismatch << ")" << std::endl;
            stream << "Gap-opening penalty: " << gap_opening_penalty << " (transformed: " << sf.penalty_open << ")" << std::endl;
            stream << "Gap-extension penalty: " << gap_extension_penalty << " (transformed: " << sf.penalty_extend << ")" << std::endl;

        }
        if (use_qgrams) stream << "Prefiltering via q-grams: activated" << std::endl;

        stream << std::endl;

    }


    void LevenshteinConfiguration::check() const {

        std::string mode = get_configuration_label(CF_LEVENSHTEIN);
        std::set<PreprocessorOption> valid_pp = {PP_FASTA, PP_FASTQ};
        std::set<ClustererOption > valid_cl = {CL_CLASSIC_SWARMER, CL_CONS_CLASSIC};
        std::set<ClusterRefinerOption > valid_cr = {CR_CLASSIC_FASTIDIOUS, CR_ITERATIVE_FASTIDIOUS, CR_CONS_LSA,
                                                    CR_CONS_LSR, CR_CONS_LSS, CR_CONS_SS};
        std::set<OutputGeneratorOption> valid_og = {OG_CLASSIC};
        std::set<QualityEncodingOption> valid_qe = {QE_SANGER, QE_ILLUMINA13, QE_ILLUMINA15, QE_ILLUMINA18};
//        std::set<DistanceOption> valid_dt = {};

        if (valid_pp.find(opt_preprocessor) == valid_pp.end()){
            std::cerr << "ERROR: Mode '" << mode << "' supports only the FASTA and FASTQ preprocessors." << std::endl;
            abort();
        }

        std::set<ClustererOption> qcl = {CL_CONS_CLASSIC};
        std::set<ClusterRefinerOption> qcr = {CR_CONS_LSA, CR_CONS_LSR, CR_CONS_LSS, CR_CONS_SS};
        if ((qcl.find(opt_clusterer) != qcl.end() || qcr.find(opt_refiner) != qcr.end()) && opt_preprocessor != PP_FASTQ) {
            std::cerr << "ERROR: The selected clusterer / refiner options require FASTQ files and the corresponding FASTQ preprocessor." << std::endl;
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

        if (opt_preprocessor == PP_FASTQ && valid_qe.find(opt_quality_encoding) == valid_qe.end()) {
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