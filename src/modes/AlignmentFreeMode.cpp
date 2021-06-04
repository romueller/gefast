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

#include <iomanip>

#include "../../include/modes/AlignmentFreeMode.hpp"

#include "../../include/Factories.hpp"

namespace GeFaST {

    AlignmentFreeConfiguration::AlignmentFreeConfiguration(int argc, const char* argv[]) {

        set_general_parameters(argc, argv);

        /* Set default values */

        opt_clusterer = (opt_clusterer == CL_IDLE) ? CL_CLASSIC_SWARMER : opt_clusterer;
        opt_refiner = (opt_refiner == CR_IDLE) ? CR_CLASSIC_FASTIDIOUS : opt_refiner;
        opt_output_generator = OG_CLASSIC;
        opt_auxiliary_data = AD_NAIVE_AUXILIARY;
        opt_refinement_auxiliary_data = RD_NAIVE_AUXILIARY;

        opt_feature = FB_UNKNOWN;
        opt_distance = DT_UNKNOWN;
        opt_amplicon_collection = AC_ARRAY;


        /* Read mode-specific configuration from file */

        std::ifstream in_stream(configuration_file);
        std::string line;
        unsigned long delim_pos;

        while (std::getline(in_stream, line).good()) {

            if (line.empty() || line.front() == '#') continue;

            delim_pos = line.find('=');
            std::string key = line.substr(0, delim_pos);

            if (key == "distance") opt_distance = get_distance_option(line.substr(delim_pos + 1));
            if (key == "representation") opt_feature = get_feature_builder_option(line.substr(delim_pos + 1));
            if (key == "amplicon_collection") opt_amplicon_collection = get_amplicon_collection_option(line.substr(delim_pos + 1));
            if (key == "auxiliary") {

                opt_auxiliary_data = get_auxiliary_data_option(line.substr(delim_pos + 1));
                opt_refinement_auxiliary_data = get_refinement_auxiliary_data_option(line.substr(delim_pos + 1));

            }

        }

        check(argc, argv);


        /* mode-specific options */

        // opt_preprocessor: user-defined, handled in set_general_parameters()
        // opt_quality_encoding: not relevant, some value is set in set_general_parameters()
        // opt_clusterer: user-defined, handled above and in set_general_parameters()
        // opt_refiner: user-defined, handled above and in set_general_parameters()
        // opt_output_generator: user-defined, handled above and in set_general_parameters()
        // opt_auxiliary_data: depends on user choices, handled above
        // opt_refinement_auxiliary_data: depends on user choices, handled above
        // opt_amplicon_collection: depends on user choice resp. given by opt_amplicon_storage
        // opt_swarm_storage: user-defined, handled in set_general_parameters()
        // opt_distance: depends on user choice, handled above
        // opt_feature: depends on user choice, handled above

        opt_amplicon_storage = AS_SIMPLE_FEATURES;

    }

    AlignmentFreeConfiguration* AlignmentFreeConfiguration::clone() const {// deep-copy clone method
        return new AlignmentFreeConfiguration(*this);
    }


    Preprocessor* AlignmentFreeConfiguration::build_preprocessor() const {
        return PreprocessorFactory::create(opt_preprocessor, opt_quality_encoding, *this);
    }

    QualityEncoding<>* AlignmentFreeConfiguration::build_quality_encoding() const {
        return new QualityEncoding<>(opt_quality_encoding);
    }

    Clusterer* AlignmentFreeConfiguration::build_clusterer() const {
        return ClustererFactory::create(opt_clusterer, *this);
    }

    ClusterRefiner* AlignmentFreeConfiguration::build_cluster_refiner() const {
        return ClusterRefinerFactory::create(opt_refiner, *this);
    }

    OutputGenerator* AlignmentFreeConfiguration::build_output_generator() const {
        return OutputGeneratorFactory::create(opt_output_generator, *this);
    }

    FeatureAmpliconStorage* AlignmentFreeConfiguration::build_amplicon_storage(const DataStatistics<>& ds) const {
        return AmpliconStorageFactory::create(opt_amplicon_storage, opt_feature, opt_amplicon_collection, ds, *this);
    }

    SwarmStorage* AlignmentFreeConfiguration::build_swarm_storage(const AmpliconStorage& amplicon_storage) const {
        return SwarmStorageFactory::create(opt_swarm_storage, amplicon_storage, *this);
    }

    AuxiliaryData* AlignmentFreeConfiguration::build_auxiliary_data(const AmpliconStorage& amplicon_storage, const numSeqs_t pool_id,
                                        const dist_t threshold) const {
        return AuxiliaryDataFactory::create(opt_auxiliary_data, amplicon_storage, pool_id, threshold, *this);
    }

    AuxiliaryData* AlignmentFreeConfiguration::build_auxiliary_data(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                                        const numSeqs_t pool_id, const dist_t threshold) const {
        return AuxiliaryDataFactory::create(opt_refinement_auxiliary_data, amplicon_storage, swarm_storage,
                                            pool_id, threshold, *this);
    }

    FeatureDistance* AlignmentFreeConfiguration::build_distance_function(const AmpliconStorage &amplicon_storage, const dist_t threshold) const {
        return DistanceFactory::create(opt_distance, amplicon_storage, threshold, *this);
    }


    void AlignmentFreeConfiguration::print(std::ostream& stream) const {

        Configuration::print(stream);

        stream << ">>> Mode-specific" << std::endl;
        stream << "Features: " << get_feature_builder_label(opt_feature) << std::endl;
//        if (opt_amplicon_collection != AC_UNKNOWN) stream << "Amplicon collection: " << get_amplicon_collection_label(opt_amplicon_collection) << std::endl;
        stream << "Auxiliary structure: " << get_auxiliary_data_label(opt_auxiliary_data) << std::endl;

        stream << std::endl;

    }

    DistanceOption AlignmentFreeConfiguration::get_opt_distance() const {
        return opt_distance;
    }




    void AlignmentFreeConfiguration::check(int argc, const char* argv[]) const {

        std::string mode = get_configuration_label(CF_ALIGNMENT_FREE);
        std::set<PreprocessorOption> valid_pp = {PP_FASTA, PP_FASTQ};
        std::set<ClustererOption > valid_cl = {CL_CLASSIC_SWARMER, CL_CONS_CLASSIC};
        std::set<ClusterRefinerOption > valid_cr = {CR_CLASSIC_FASTIDIOUS, CR_ITERATIVE_FASTIDIOUS, CR_CONS_LSA,
                                                    CR_CONS_LSR, CR_CONS_LSS, CR_CONS_SS};
        std::set<OutputGeneratorOption> valid_og = {OG_CLASSIC};
        std::set<QualityEncodingOption> valid_qe = {QE_SANGER, QE_ILLUMINA13, QE_ILLUMINA15, QE_ILLUMINA18};
        std::set<DistanceOption> valid_dt = {DT_MANHATTAN, DT_EUCLIDEAN, DT_COSINE, DT_PEARSON};
        std::set<FeatureBuilderOption> valid_fb = {FB_WCV, FB_CPF, FB_DET, FB_BBC, FB_2DW, FB_CGR, FB_MCGR, FB_3DCGR};
        std::set<AmpliconCollectionOption> valid_ac = {AC_UNKNOWN, AC_SIMPLE, AC_ARRAY};

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

        if (opt_distance == DT_UNKNOWN) {
            std::cerr << "ERROR: A 'distance' option has to be specified for mode '" << mode << "'." << std::endl;
            abort();
        }

        if (valid_dt.find(opt_distance) == valid_dt.end()) {
            std::cerr << "ERROR: Invalid 'distance' option for mode '" << mode << "'." << std::endl;
            abort();
        }

        if (valid_ac.find(opt_amplicon_collection) == valid_ac.end()) {
            std::cerr << "ERROR: Invalid 'amplicon_collection' option for mode '" << mode << "'." << std::endl;
            abort();
        }

        std::set<AuxiliaryDataOption> adc = {AD_KDTREE_AUXILIARY};
        std::set<RefinementAuxiliaryDataOption> adr = {RD_KDTREE_AUXILIARY};
        std::set<DistanceOption> ndt = {DT_MANHATTAN, DT_EUCLIDEAN};
        if ((adc.find(opt_auxiliary_data) != adc.end() || adr.find(opt_refinement_auxiliary_data) != adr.end()) && (ndt.find(opt_distance) == ndt.end())) {
            std::cerr << "ERROR: The selected 'distance' option is not compatible with the selected 'auxiliary' option." << std::endl;
            abort();
        }


        // exit with error when no threshold was specified
        bool threshold_specified = false;
        for (auto i = 0; i < argc && !threshold_specified; i++) {
            threshold_specified = (strcmp(argv[i], "-t") == 0) || (strcmp(argv[i], "--threshold") == 0);
        }

        std::ifstream in_stream(configuration_file);
        std::string line;
        unsigned long delim_pos;
        while (std::getline(in_stream, line).good()) {

            if (line.empty() || line.front() == '#') continue;

            delim_pos = line.find('=');
            std::string key = line.substr(0, delim_pos);

            threshold_specified |= (key == "threshold");

        }

        if (!threshold_specified) {
            std::cerr << "ERROR: No clustering 'threshold' was specified (mode '" << mode << "' does not provide a default)." << std::endl;
            abort();
        }

        // FeatureBuilder-related checks (choice, compatibility, existence of parameters)
        if (opt_feature == FB_UNKNOWN) {
            std::cerr << "ERROR: A 'representation' option has to be specified for mode '" << mode << "'." << std::endl;
            abort();
        }

        if (valid_fb.find(opt_feature) == valid_fb.end()) {
            std::cerr << "ERROR: Invalid 'representation' option for mode '" << mode << "'." << std::endl;
            abort();
        }
        std::string sorted_alphabet = alphabet;
        std::sort(sorted_alphabet.begin(), sorted_alphabet.end());

        switch (opt_feature) {

            case FB_WCV: {
                if (misc.find("wcv_word_length") == misc.end()) {
                    std::cerr << "ERROR: Feature representation " << get_feature_builder_label(opt_feature)
                              << " requires parameter 'wcv_word_length' (via misc)." << std::endl;
                    abort();
                }
                break;
            }

            case FB_DET: {
                if (misc.find("det_alpha") == misc.end() && misc.find("det_pref_dist") == misc.end()) {
                    std::cerr << "ERROR: Feature representation " << get_feature_builder_label(opt_feature)
                              << " requires parameter 'det_alpha' or 'det_pref_dist' (via misc)." << std::endl;
                    abort();
                }
                if (misc.find("det_alpha") != misc.end() && misc.find("det_pref_dist") != misc.end()) {
                    std::cerr << "ERROR: Only one of the two parameters 'det_alpha' or 'det_pref_dist' "
                              << "can be specified when using feature representation " << get_feature_builder_label(opt_feature) << std::endl;
                    abort();
                }
                break;
            }

            case FB_BBC: {
                if (misc.find("bbc_max_dist") == misc.end()) {
                    std::cerr << "ERROR: Feature representation " << get_feature_builder_label(opt_feature)
                              << " requires parameter 'bbc_max_dist' (via misc)." << std::endl;
                    abort();
                }
                break;
            }

            case FB_2DW: {
                if (misc.find("2dw_alpha") == misc.end() || misc.find("2dw_beta") == misc.end()) {
                    std::cerr << "ERROR: Feature representation " << get_feature_builder_label(opt_feature)
                              << " requires parameters '2dw_alpha' and '2dw_beta' (via misc)." << std::endl;
                    abort();
                }
                if (sorted_alphabet != "ACGT" && sorted_alphabet != "ACGU") {
                    std::cerr << "ERROR: Feature representation " << get_feature_builder_label(opt_feature)
                              << " can only be used with alphabet ACGT or ACGU." << std::endl;
                    abort();
                }
                break;
            }

            case FB_CPF: // fallthrough
            case FB_CGR: // fallthrough
            case FB_MCGR: // fallthrough
            case FB_3DCGR: {
                if (sorted_alphabet != "ACGT" && sorted_alphabet != "ACGU") {
                    std::cerr << "ERROR: Feature representation " << get_feature_builder_label(opt_feature)
                              << " can only be used with alphabet ACGT or ACGU." << std::endl;
                    abort();
                }
                break;
            }

        }

    }

}