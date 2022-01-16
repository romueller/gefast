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

#include "../../include/modes/SpaceLevenshteinMode.hpp"

#include "../../include/Factories.hpp"
#include "../../include/Distances.hpp"

namespace GeFaST {

    numSeqs_t ChunkLengths::get_capa(numSeqs_t capa) const {
        return std::max(1ll, std::llroundf(init_rel_capa * capa));
    }

    SpaceLevenshteinConfiguration::SpaceLevenshteinConfiguration(int argc, const char* argv[]) {

        set_general_parameters(argc, argv);

        /* Set default values */

        if (main_threshold == 0) main_threshold = 1;
        num_extra_segments = 1;
        use_score = false;
        use_qgrams = false;

        prefix_len = 0;

        chunk_lengths.ids = 8;
        chunk_lengths.dists = 8;
        chunk_lengths.gens = 8;
        chunk_lengths.rads = 8;
        chunk_lengths.starts = 8;
        chunk_lengths.sizes = 8;
        chunk_lengths.masses = 8;
        chunk_lengths.differents = 8;
        chunk_lengths.singletons = 8;
        chunk_lengths.max_rads = 8;
        chunk_lengths.init_rel_capa = 0.3;
        chunk_lengths.level_factor = 0.5;
        chunk_lengths.extend_factor = 2.0;
        chunk_lengths.opt = true;

        chunk_lengths.headers_chunk = 8;
        chunk_lengths.headers_level = 0.5;
        chunk_lengths.headers_extend = 2.0;
        chunk_lengths.headers_opt = true;

        chunk_lengths.prefixes_chunk = 8;
        chunk_lengths.prefixes_level = 0.5;
        chunk_lengths.prefixes_extend = 2.0;
        chunk_lengths.prefixes_opt = true;

        chunk_lengths.seqs_chunk = 8;
        chunk_lengths.seqs_level = 0.5;
        chunk_lengths.seqs_extend = 2.0;
        chunk_lengths.seqs_opt = true;

        chunk_lengths.qgram_chunk = 8;
        chunk_lengths.qgram_level = 0.5;
        chunk_lengths.qgram_extend = 2.0;
        chunk_lengths.qgram_opt = true;

        // TODO? improve default handling by postponing it until the selected k^2-tree implementation is known
        tree_parameters.k = 2;
        tree_parameters.kr = 2;
        tree_parameters.kc = 3;
        tree_parameters.uk = 2;
        tree_parameters.uh = 3;
        tree_parameters.lk = 4;
        tree_parameters.mb = 20;
        tree_parameters.static_rels = true;
        tree_parameters.deep = false;

        opt_amplicon_storage = AS_UNKNOWN;
        opt_clusterer = (opt_clusterer == CL_IDLE) ? CL_CLASSIC_SWARMER : opt_clusterer;
        opt_refiner = (opt_refiner == CR_IDLE) ? CR_CLASSIC_FASTIDIOUS : opt_refiner;
        opt_output_generator = OG_CLASSIC;
        opt_auxiliary_data = AD_SEGMENT_FILTER;
        opt_refinement_auxiliary_data = RD_SEGMENT_FILTER;


        /* Read mode-specific configuration from file */

        std::ifstream in_stream(configuration_file);
        std::string line;
        unsigned long delim_pos;

        while (std::getline(in_stream, line).good()) {

            if (line.empty() || line.front() == '#') continue;

            delim_pos = line.find('=');
            std::string key = line.substr(0, delim_pos);

            if (key == "num_extra_segments") num_extra_segments = std::stoul(line.substr(delim_pos + 1));
            if (key == "use_score") use_score = (line.substr(delim_pos + 1) == "1");
            if (key == "use_qgrams") use_qgrams = (line.substr(delim_pos + 1) == "1");
            if (key == "amplicon_collection") opt_amplicon_collection = get_amplicon_collection_option(line.substr(delim_pos + 1));
            if (key == "amplicon_storage") opt_amplicon_storage = get_amplicon_storage_option(line.substr(delim_pos + 1));
            if (key == "auxiliary_data") opt_auxiliary_data = get_auxiliary_data_option(line.substr(delim_pos + 1));
            if (key == "refinement_auxiliary_data") opt_refinement_auxiliary_data = get_refinement_auxiliary_data_option(line.substr(delim_pos + 1));

            if (key == "prefix_length") prefix_len = std::stoul(line.substr(delim_pos + 1));

            if (key == "chunk_length_ids") chunk_lengths.ids = std::stoul(line.substr(delim_pos + 1));
            if (key == "chunk_length_dists") chunk_lengths.dists = std::stoul(line.substr(delim_pos + 1));
            if (key == "chunk_length_gens") chunk_lengths.gens = std::stoul(line.substr(delim_pos + 1));
            if (key == "chunk_length_rads") chunk_lengths.rads = std::stoul(line.substr(delim_pos + 1));
            if (key == "chunk_length_starts") chunk_lengths.starts = std::stoul(line.substr(delim_pos + 1));
            if (key == "chunk_length_sizes") chunk_lengths.sizes = std::stoul(line.substr(delim_pos + 1));
            if (key == "chunk_length_masses") chunk_lengths.masses = std::stoul(line.substr(delim_pos + 1));
            if (key == "chunk_length_differents") chunk_lengths.differents = std::stoul(line.substr(delim_pos + 1));
            if (key == "chunk_length_singletons") chunk_lengths.singletons = std::stoul(line.substr(delim_pos + 1));
            if (key == "chunk_length_max_rads") chunk_lengths.max_rads = std::stoul(line.substr(delim_pos + 1));
            if (key == "init_rel_capa") chunk_lengths.init_rel_capa = std::stof(line.substr(delim_pos + 1));
            if (key == "level_factor") chunk_lengths.level_factor = std::stof(line.substr(delim_pos + 1));
            if (key == "extend_factor") chunk_lengths.extend_factor = std::stof(line.substr(delim_pos + 1));
            if (key == "optimise_dacs") chunk_lengths.opt = (line.substr(delim_pos + 1) == "1");

            if (key == "headers_chunk") chunk_lengths.headers_chunk = std::stoul(line.substr(delim_pos + 1));
            if (key == "headers_level") chunk_lengths.headers_level = std::stof(line.substr(delim_pos + 1));
            if (key == "headers_extend") chunk_lengths.headers_extend = std::stof(line.substr(delim_pos + 1));
            if (key == "headers_opt") chunk_lengths.headers_opt = (line.substr(delim_pos + 1) == "1");

            if (key == "prefixes_chunk") chunk_lengths.prefixes_chunk = std::stoul(line.substr(delim_pos + 1));
            if (key == "prefixes_level") chunk_lengths.prefixes_level = std::stof(line.substr(delim_pos + 1));
            if (key == "prefixes_extend") chunk_lengths.prefixes_extend = std::stof(line.substr(delim_pos + 1));
            if (key == "prefixes_opt") chunk_lengths.prefixes_opt = (line.substr(delim_pos + 1) == "1");

            if (key == "seqs_chunk") chunk_lengths.seqs_chunk = std::stoul(line.substr(delim_pos + 1));
            if (key == "seqs_level") chunk_lengths.seqs_level = std::stof(line.substr(delim_pos + 1));
            if (key == "seqs_extend") chunk_lengths.seqs_extend = std::stof(line.substr(delim_pos + 1));
            if (key == "seqs_opt") chunk_lengths.seqs_opt = (line.substr(delim_pos + 1) == "1");

            if (key == "qgram_chunk") chunk_lengths.qgram_chunk = std::stoul(line.substr(delim_pos + 1));
            if (key == "qgram_level") chunk_lengths.qgram_level = std::stof(line.substr(delim_pos + 1));
            if (key == "qgram_extend") chunk_lengths.qgram_extend = std::stof(line.substr(delim_pos + 1));
            if (key == "qgram_opt") chunk_lengths.qgram_opt = (line.substr(delim_pos + 1) == "1");

            if (key == "tree_parameter_k") tree_parameters.k = std::stoul(line.substr(delim_pos + 1));
            if (key == "tree_parameter_kr") tree_parameters.kr = std::stoul(line.substr(delim_pos + 1));
            if (key == "tree_parameter_kc") tree_parameters.kc = std::stoul(line.substr(delim_pos + 1));
            if (key == "tree_parameter_uk") tree_parameters.uk = std::stoul(line.substr(delim_pos + 1));
            if (key == "tree_parameter_uh") tree_parameters.uh = std::stoul(line.substr(delim_pos + 1));
            if (key == "tree_parameter_lk") tree_parameters.lk = std::stoul(line.substr(delim_pos + 1));
            if (key == "tree_parameter_mb") tree_parameters.mb = std::stoul(line.substr(delim_pos + 1));
            if (key == "tree_parameter_static_rels") tree_parameters.static_rels = (line.substr(delim_pos + 1) == "1");
            if (key == "tree_parameter_deep") tree_parameters.deep = (line.substr(delim_pos + 1) == "1");

        }

        check();


        /* mode-specific data structures and settings */

        // opt_preprocessor: user-defined, handled in set_general_parameters()
        // opt_quality_encoding: user-defined, handled in set_general_parameters()
        // opt_clusterer: user-defined, handled above and in set_general_parameters()
        // opt_refiner: user-defined, handled above and in set_general_parameters()
        // opt_output_generator: user-defined, handled above and in set_general_parameters()
        // opt_amplicon_collection: used-defined
        // opt_amplicon_storage: used-defined
        // opt_swarm_storage: user-defined, handled in set_general_parameters()

        if (opt_amplicon_storage == AS_PREPARED_QGRAM_LENGTH_POOLS) opt_amplicon_collection = AC_ARRAY_QGRAM;

        opt_distance = (use_score) ? DT_BOUNDED_SCORE_LEVENSHTEIN : DT_BOUNDED_LEVENSHTEIN;

    }

    SpaceLevenshteinConfiguration* SpaceLevenshteinConfiguration::clone() const {
        return new SpaceLevenshteinConfiguration(*this);
    }

    Preprocessor* SpaceLevenshteinConfiguration::build_preprocessor() const {
        return PreprocessorFactory::create(opt_preprocessor, opt_quality_encoding, *this);
    }

    QualityEncoding<>* SpaceLevenshteinConfiguration::build_quality_encoding() const {
        return new QualityEncoding<>(opt_quality_encoding);
    }

    Clusterer* SpaceLevenshteinConfiguration::build_clusterer() const {
        return ClustererFactory::create(opt_clusterer, *this);
    }

    ClusterRefiner* SpaceLevenshteinConfiguration::build_cluster_refiner() const {
        return ClusterRefinerFactory::create(opt_refiner, *this);
    }

    OutputGenerator* SpaceLevenshteinConfiguration::build_output_generator() const {
        return OutputGeneratorFactory::create(opt_output_generator, *this);
    }

    AmpliconStorage* SpaceLevenshteinConfiguration::build_amplicon_storage(const DataStatistics<>& ds) const {
        return AmpliconStorageFactory::create(opt_amplicon_storage, opt_amplicon_collection, ds, *this);
    }

    SwarmStorage* SpaceLevenshteinConfiguration::build_swarm_storage(const AmpliconStorage& amplicon_storage) const {
        return SwarmStorageFactory::create(opt_swarm_storage, amplicon_storage, *this);
    }

    AuxiliaryData* SpaceLevenshteinConfiguration::build_auxiliary_data(const AmpliconStorage& amplicon_storage,
            const numSeqs_t pool_id, const dist_t threshold) const {
        return AuxiliaryDataFactory::create(opt_auxiliary_data, amplicon_storage, pool_id, threshold, num_extra_segments, *this);
    }

    AuxiliaryData* SpaceLevenshteinConfiguration::build_auxiliary_data(const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, const numSeqs_t pool_id, const dist_t threshold) const {

        return AuxiliaryDataFactory::create(opt_refinement_auxiliary_data, amplicon_storage, swarm_storage,
                pool_id, threshold, num_extra_segments, *this);

    }

    Distance* SpaceLevenshteinConfiguration::build_distance_function(const AmpliconStorage &amplicon_storage,
                                                                const dist_t threshold) const {
        return DistanceFactory::create(opt_distance, amplicon_storage, threshold, *this);
    }


    void SpaceLevenshteinConfiguration::print(std::ostream& stream) const {

        Configuration::print(stream);

        ScoringFunction sf(match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty);

        stream << ">>> Mode-specific" << std::endl;
        stream << "Extra segments: " << num_extra_segments << std::endl;
        stream << "Segment filter: one-way" << std::endl;
        stream << "Use scoring function: " << use_score << std::endl;
        if (use_score) {

            stream << "Match reward: " << match_reward << " (transformed: 0)" << std::endl;
            stream << "Mismatch penalty: " << mismatch_penalty << " (transformed: " << sf.penalty_mismatch << ")" << std::endl;
            stream << "Gap-opening penalty: " << gap_opening_penalty << " (transformed: " << sf.penalty_open << ")" << std::endl;
            stream << "Gap-extension penalty: " << gap_extension_penalty << " (transformed: " << sf.penalty_extend << ")" << std::endl;

        }
        if (use_qgrams) stream << "Prefiltering via q-grams: activated" << std::endl;
        stream << "Swarm storage: " << get_swarm_storage_label(opt_swarm_storage) << std::endl;
        stream << "Amplicon storage: " << get_amplicon_storage_label(opt_amplicon_storage) << std::endl;
        stream << "Amplicon collection: " << get_amplicon_collection_label(opt_amplicon_collection) << std::endl;
        if (opt_amplicon_collection == AC_SPACE_ID_PREF) stream << "Prefix length (headers): " << static_cast<int>(prefix_len) << std::endl;
        if (opt_amplicon_collection == AC_SPACE_ID_PREF_DACS) {
            stream << "Prefix length: " << static_cast<int>(prefix_len) << std::endl;
            stream << "Cursor parameters (headers): " << chunk_lengths.headers_chunk << " / " << chunk_lengths.headers_level
                << " / " << chunk_lengths.headers_extend << " / " << chunk_lengths.headers_opt << std::endl;
            stream << "Prefix parameters (headers): " << chunk_lengths.prefixes_chunk << " / " << chunk_lengths.prefixes_level
                << " / " << chunk_lengths.prefixes_extend << " / " << chunk_lengths.prefixes_opt << std::endl;
        }
        if (opt_amplicon_collection == AC_SPACE_ID_HUFF_DACS) {
            stream << "Cursor parameters (headers): " << chunk_lengths.headers_chunk << " / " << chunk_lengths.headers_level
                << " / " << chunk_lengths.headers_extend << " / " << chunk_lengths.headers_opt << std::endl;
        }
        if (opt_amplicon_collection == AC_SPACE_SEQ_DACS) {
            stream << "Cursor parameters (sequences): " << chunk_lengths.seqs_chunk << " / " << chunk_lengths.seqs_level
                << " / " << chunk_lengths.seqs_extend << " / " << chunk_lengths.seqs_opt << std::endl;
        }
        if (opt_amplicon_collection == AC_SPACE_QGRAM_DACS) {
            stream << "Cursor parameters (q-grams): " << chunk_lengths.qgram_chunk << " / " << chunk_lengths.qgram_level
                << " / " << chunk_lengths.qgram_extend << " / " << chunk_lengths.qgram_opt << std::endl;
        }
        if (opt_amplicon_collection == AC_SPACE_FULL) {
            stream << "Prefix length: " << static_cast<int>(prefix_len) << std::endl;
            stream << "Cursor parameters (headers): " << chunk_lengths.headers_chunk << " / " << chunk_lengths.headers_level
                   << " / " << chunk_lengths.headers_extend << " / " << chunk_lengths.headers_opt << std::endl;
            stream << "Prefix parameters (headers): " << chunk_lengths.prefixes_chunk << " / " << chunk_lengths.prefixes_level
                   << " / " << chunk_lengths.prefixes_extend << " / " << chunk_lengths.prefixes_opt << std::endl;
            stream << "Cursor parameters (sequences): " << chunk_lengths.seqs_chunk << " / " << chunk_lengths.seqs_level
                   << " / " << chunk_lengths.seqs_extend << " / " << chunk_lengths.seqs_opt << std::endl;
            stream << "Cursor parameters (q-grams): " << chunk_lengths.qgram_chunk << " / " << chunk_lengths.qgram_level
                << " / " << chunk_lengths.qgram_extend << " / " << chunk_lengths.qgram_opt << std::endl;
        }
        stream << "Auxiliary data (clustering): " << get_auxiliary_data_label(opt_auxiliary_data) << std::endl;
        stream << "Auxiliary data (refinement): " << get_refinement_auxiliary_data_label(opt_refinement_auxiliary_data) << std::endl;

        stream << std::endl;
        if (opt_swarm_storage == SS_FORWARDING_DACS_PER_POOL) stream << "Chunk lengths (swarms): "
            << chunk_lengths.ids << " / " << chunk_lengths.dists << " / " << chunk_lengths.gens << " / " << chunk_lengths.rads << std::endl;
        if (opt_swarm_storage == SS_FORWARDING_FULL_DACS_PER_POOL) stream << "Chunk lengths (swarms): "
            << chunk_lengths.ids << " / " << chunk_lengths.dists << " / " << chunk_lengths.gens << " / " << chunk_lengths.rads << ", "
            << chunk_lengths.starts << " / " << chunk_lengths.sizes << " / " << chunk_lengths.masses << " / "
            << chunk_lengths.differents << " / " << chunk_lengths.singletons << " / " << chunk_lengths.max_rads << std::endl;
        if (opt_swarm_storage == SS_FORWARDING_FULL_DACS_PER_POOL) stream << "Initial relative capacity (swarms): " << chunk_lengths.init_rel_capa << std::endl;
        if (opt_swarm_storage == SS_FORWARDING_DACS_PER_POOL || opt_swarm_storage == SS_FORWARDING_FULL_DACS_PER_POOL) {
            stream << "Level factor (swarms): " << chunk_lengths.level_factor << std::endl;
            stream << "Extend factor (swarms): " << chunk_lengths.extend_factor << std::endl;
            stream << "Optimise DACs (swarms): " << chunk_lengths.opt << std::endl;
        }

        stream << std::endl;
        if (opt_auxiliary_data == AD_SPACE_SEGMENT_FILTER_BASIC) stream << "k^2-tree parameter k: " << tree_parameters.k << std::endl;
        if (opt_auxiliary_data == AD_SPACE_SEGMENT_FILTER_HYBRID) {
            stream << "k^2-tree parameter uk: " << tree_parameters.uk << std::endl;
            stream << "k^2-tree parameter uh: " << tree_parameters.uh << std::endl;
            stream << "k^2-tree parameter lk: " << tree_parameters.lk << std::endl;
        }
        if (opt_auxiliary_data == AD_SPACE_SEGMENT_FILTER_RECT || opt_auxiliary_data == AD_SPACE_SEGMENT_FILTER_PART) {
            stream << "k^2-tree parameter kr: " << tree_parameters.kr << std::endl;
            stream << "k^2-tree parameter kc: " << tree_parameters.kc << std::endl;
        }
        if (opt_auxiliary_data == AD_SPACE_SEGMENT_FILTER_VARI) {
            stream << "k^2-tree parameter kr: " << tree_parameters.kr << std::endl;
            stream << "k^2-tree parameter kc: " << tree_parameters.kc << std::endl;
            stream << "k^2-tree parameter mb: " << tree_parameters.mb << std::endl;
        }
        if (opt_auxiliary_data != AD_SEGMENT_FILTER || opt_refinement_auxiliary_data != RD_SEGMENT_FILTER) {
            stream << "Static k^2-trees: " << tree_parameters.static_rels << std::endl;
            stream << "Deep set_null: " << tree_parameters.deep << std::endl;
        }

        stream << std::endl;

    }


    void SpaceLevenshteinConfiguration::check() const {

        std::string mode = get_configuration_label(CF_LEVENSHTEIN);
        std::set<PreprocessorOption> valid_pp = {PP_FASTA, PP_FASTQ};
        std::set<ClustererOption > valid_cl = {CL_CLASSIC_SWARMER};
        std::set<ClusterRefinerOption > valid_cr = {CR_CLASSIC_FASTIDIOUS};
        std::set<OutputGeneratorOption> valid_og = {OG_CLASSIC};
        std::set<AuxiliaryDataOption> valid_ad = {AD_SEGMENT_FILTER, AD_SPACE_SEGMENT_FILTER_NAIVE, AD_SPACE_SEGMENT_FILTER_BASIC,
                                                  AD_SPACE_SEGMENT_FILTER_HYBRID, AD_SPACE_SEGMENT_FILTER_RECT,
                                                  AD_SPACE_SEGMENT_FILTER_PART, AD_SPACE_SEGMENT_FILTER_VARI};
        std::set<RefinementAuxiliaryDataOption> valid_rd = {RD_SEGMENT_FILTER, RD_SPACE_SEGMENT_FILTER_NAIVE, RD_SPACE_SEGMENT_FILTER_BASIC,
                                                            RD_SPACE_SEGMENT_FILTER_HYBRID, RD_SPACE_SEGMENT_FILTER_RECT,
                                                            RD_SPACE_SEGMENT_FILTER_PART, RD_SPACE_SEGMENT_FILTER_VARI};
        std::set<QualityEncodingOption> valid_qe = {QE_SANGER, QE_ILLUMINA13, QE_ILLUMINA15, QE_ILLUMINA18};
        std::set<AmpliconCollectionOption> valid_ac = {AC_UNKNOWN, AC_ARRAY_QGRAM,
                                                       AC_SPACE_REFERENCE, AC_SPACE_FULL,
                                                       AC_SPACE_ID_PREF, AC_SPACE_ID_HUFF, AC_SPACE_ID_PREF_DACS, AC_SPACE_ID_HUFF_DACS,
                                                       AC_SPACE_SEQ, AC_SPACE_SEQ_DACS,
                                                       AC_SPACE_AB, AC_SPACE_AB_BITS,
                                                       AC_SPACE_LEN, AC_SPACE_LEN_OFFSET,
                                                       AC_SPACE_QGRAM_DACS,
                                                       AC_SPACE_AB_SPLIT, AC_SPACE_AB_BITS_SPLIT, AC_SPACE_LEN_SPLIT, AC_SPACE_LEN_OFFSET_SPLIT};
        std::set<AmpliconStorageOption> valid_as = {AS_PREPARED_QGRAM_LENGTH_POOLS, AS_SPACE_FLEXIBLE};
        std::set<SwarmStorageOption> valid_ss = {SS_FORWARDING_PER_POOL, SS_FORWARDING_DACS_PER_POOL, SS_FORWARDING_FULL_DACS_PER_POOL};

        if (valid_pp.find(opt_preprocessor) == valid_pp.end()){
            std::cerr << "ERROR: Mode '" << mode << "' supports only the FASTA and FASTQ preprocessors." << std::endl;
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

        if (valid_ad.find(opt_auxiliary_data) == valid_ad.end()) {
            std::cerr << "ERROR: Invalid 'auxiliary_data' option for mode '" << mode << "'." << std::endl;
            abort();
        }

        if (valid_rd.find(opt_refinement_auxiliary_data) == valid_rd.end()) {
            std::cerr << "ERROR: Invalid 'refinement_auxiliary_data' option for mode '" << mode << "'." << std::endl;
            abort();
        }

        if (opt_preprocessor == PP_FASTQ && valid_qe.find(opt_quality_encoding) == valid_qe.end()) {
            std::cerr << "ERROR: Invalid 'quality_encoding' option for mode '" << mode << "'." << std::endl;
            abort();
        }

        if (valid_ac.find(opt_amplicon_collection) == valid_ac.end()) {
            std::cerr << "ERROR: Invalid 'amplicon_collection' option for mode '" << mode << "'." << std::endl;
            abort();
        }

        if (valid_as.find(opt_amplicon_storage) == valid_as.end()) {
            std::cerr << "ERROR: Invalid 'amplicon_storage' option for mode '" << mode << "'." << std::endl;
            abort();
        }

        if (valid_ss.find(opt_swarm_storage) == valid_ss.end()) {
            std::cerr << "ERROR: Invalid 'swarm_storage' option for mode '" << mode << "'." << std::endl;
            abort();
        }

        if (opt_amplicon_collection == AC_UNKNOWN && opt_amplicon_storage == AS_SPACE_FLEXIBLE) {
            std::cerr << "ERROR: The selected amplicon-storage option requires the specification of the amplicon-collection option." << std::endl;
            abort();
        }

        if (opt_amplicon_collection == AC_SPACE_ID_PREF && prefix_len == 0) {
            std::cerr << "ERROR: When using the PrefixStringArray, the 'prefix_length' option must be specified." << std::endl;
            abort();
        }

    }

}