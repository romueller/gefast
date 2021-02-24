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

#ifndef GEFAST_OPTIONS_HPP
#define GEFAST_OPTIONS_HPP

#include <map>
#include <string>

namespace GeFaST {

    /*
     * This file provides the options lists for the key data structures
     * underlying the different steps of GeFaST.
     *
     * Options are prefixed with their category:
     *  - Configuration             CF_
     *  - AmpliconStorage           AS_
     *  - AmpliconCollection        AC_
     *  - Distance                  DT_
     *  - Boosting                  BO_
     *  - Preprocessor              PP_
     *  - QualityEncoding           QE_
     *  - Clusterer                 CL_
     *  - ClusterRefiner            CR_
     *  - OutputGenerator           OG_
     *  - SwarmStorage              SS_
     *  - AuxiliaryData             AD_
     *  - RefinementAuxiliaryData   RD_
     *
     * Categories which allow a choice by the user or
     * which are printed as part of the configuration print
     * also have methods for converting between the enum and
     * a string label (and vice versa).
     * For category X, get_X_option(...) turns a string label into the
     * corresponding enum option (or a category-specific 'unknown' option),
     * while get_X_label(...) returns the string label for the enum option.
     */

    /* === Configuration === */

    enum ConfigurationOption {
        CF_UNKNOWN,
        CF_PREPROCESSING_ONLY,
        CF_DEREPLICATION,
        CF_LEVENSHTEIN,
        CF_QUALITY_LEVENSHTEIN,
        CF_ALIGNMENT_SCORE,
        CF_QUALITY_ALIGNMENT_SCORE,
        CF_CONSISTENCY
    };

    extern std::map<std::string, ConfigurationOption> config_map;

    ConfigurationOption get_configuration_option(const std::string& type);
    std::string get_configuration_label(const ConfigurationOption opt);


    /* === AmpliconStorage === */

    enum AmpliconStorageOption {
        AS_UNKNOWN,
        AS_SIMPLE_LENGTH_POOLS,
        AS_PREPARED_LENGTH_POOLS,
        AS_PREPARED_LENGTH_POOLS_QUALITY,
        AS_PREPARED_POOLS,
        AS_PREPARED_POOLS_QUALITY,
        AS_PREPARED_QGRAM_LENGTH_POOLS,
        AS_PREPARED_QGRAM_LENGTH_POOLS_QUALITY
    };

    extern std::map<std::string, AmpliconStorageOption> ampl_storage_map;

    AmpliconStorageOption get_amplicon_storage_option(const std::string& type);
    std::string get_amplicon_storage_label(const AmpliconStorageOption opt);


    /* === AmpliconCollection === */

    enum AmpliconCollectionOption {
        AC_UNKNOWN,
        AC_SIMPLE,
        AC_ARRAY,
        AC_ARRAY_QUALITY,
        AC_ARRAY_QGRAM,
        AC_ARRAY_QGRAM_QUALITY
    };

    extern std::map<std::string, AmpliconCollectionOption> ampl_collection_map;

    AmpliconCollectionOption get_amplicon_collection_option(const std::string& type);
    std::string get_amplicon_collection_label(const AmpliconCollectionOption opt);


    /* === Distance === */

    enum DistanceOption {
        DT_UNKNOWN,
        DT_BOUNDED_LEVENSHTEIN,
        DT_BOUNDED_SCORE_LEVENSHTEIN,
        DT_BOUNDED_SCORE,
        DT_BOUNDED_BANDED_SCORE,
        DT_QW_CLEMENT,
        DT_QW_MALDE_A,
        DT_QW_MALDE_B,
        DT_QW_MALDE_C,
        DT_QW_FRITH,
        DT_QW_KIM_A,
        DT_QW_KIM_B,
        DT_QW_CONVERGE_A,
        DT_QW_CONVERGE_B,
        DT_QWS_CLEMENT,
        DT_QWS_MALDE_A,
        DT_QWS_MALDE_B,
        DT_QWS_MALDE_C,
        DT_QWS_FRITH,
        DT_QWS_KIM_A,
        DT_QWS_KIM_B,
        DT_QWS_CONVERGE_A,
        DT_QWS_CONVERGE_B,
        DT_QWBS_CLEMENT,
        DT_QWBS_MALDE_A,
        DT_QWBS_MALDE_B,
        DT_QWBS_MALDE_C,
        DT_QWBS_FRITH,
        DT_QWBS_KIM_A,
        DT_QWBS_KIM_B,
        DT_QWBS_CONVERGE_A,
        DT_QWBS_CONVERGE_B
    };

    extern std::map<std::string, DistanceOption> distance_map;

    DistanceOption get_distance_option(const std::string &type);
    std::string get_distance_label(const DistanceOption opt);


    /* === Boosting === */

    enum BoostingOption {
        BO_UNKNOWN,
        BO_NO_EFFECT,
        BO_LINEAR,
        BO_MULT,
        BO_ROOT
    };

    extern std::map<std::string, BoostingOption> boosting_map;

    BoostingOption get_boosting_option(const std::string& type);
    std::string get_boosting_label(const BoostingOption opt);


    /* === Preprocessor === */

    enum PreprocessorOption {
        PP_UNKNOWN,
        PP_FASTA,
        PP_FASTQ
    };

    extern std::map<std::string, PreprocessorOption> prep_map;

    PreprocessorOption get_preprocessor_option(const std::string& type);
    std::string get_preprocessor_label(const PreprocessorOption opt);


    /* === QualityEncoding === */

    enum QualityEncodingOption {
        QE_UNKNOWN,
        QE_SANGER,
        QE_ILLUMINA13,
        QE_ILLUMINA15,
        QE_ILLUMINA18
    };

    extern std::map<std::string, QualityEncodingOption> enc_map;

    QualityEncodingOption get_quality_encoding_option(const std::string& type);
    std::string get_quality_encoding_label(const QualityEncodingOption opt);


    /* === Clusterer === */

    enum ClustererOption {
        CL_UNKNOWN,
        CL_IDLE,
        CL_CLASSIC_SWARMER,
        CL_DEREPLICATOR,
        CL_DADA2,
        CL_CONS_CLASSIC,
        CL_CONS_SWARMER
    };

    extern std::map<std::string, ClustererOption> clusterer_map;

    ClustererOption get_clusterer_option(const std::string& type);
    std::string get_clusterer_label(const ClustererOption opt);


    /* === ClusterRefiner === */

    enum ClusterRefinerOption {
        CR_UNKNOWN,
        CR_IDLE,
        CR_CLASSIC_FASTIDIOUS,
        CR_ITERATIVE_FASTIDIOUS,
        CR_CONS_LSA,
        CR_CONS_LSR,
        CR_CONS_LSS,
        CR_CONS_SS
    };

    extern std::map<std::string, ClusterRefinerOption> refiner_map;

    ClusterRefinerOption get_refiner_option(const std::string& type);
    std::string get_refiner_label(const ClusterRefinerOption opt);


    /* === OutputGenerator === */

    enum OutputGeneratorOption {
        OG_UNKNOWN,
        OG_IDLE,
        OG_CLASSIC,
        OG_DEREPLICATION
    };

    extern std::map<std::string, OutputGeneratorOption> output_generator_map;

    OutputGeneratorOption get_output_generator_option(const std::string& type);
    std::string get_output_generator_label(const OutputGeneratorOption opt);


    /* === SwarmStorage === */

    enum SwarmStorageOption {
        SS_UNKNOWN,
        SS_SIMPLE_PER_POOL,
        SS_SIMPLE_COMB_PER_POOL,
        SS_FORWARDING_PER_POOL,
        SS_FORWARDING_COMB_PER_POOL
    };

    extern std::map<std::string, SwarmStorageOption> swarm_storage_map;

    SwarmStorageOption get_swarm_storage_option(const std::string& type);
    std::string get_swarm_storage_label(const SwarmStorageOption opt);


    /* === AuxiliaryData === */

    enum AuxiliaryDataOption {
        AD_UNKNOWN,
        AD_NAIVE_AUXILIARY,
        AD_SEGMENT_FILTER,
        AD_2W_SEGMENT_FILTER,
        AD_SCORE_SEGMENT_FILTER
    };

    extern std::map<std::string, AuxiliaryDataOption> auxiliary_data_map;

    AuxiliaryDataOption get_auxiliary_data_option(const std::string& type);
    std::string get_auxiliary_data_label(const AuxiliaryDataOption opt);


    /* === RefinementAuxiliaryData === */

    enum RefinementAuxiliaryDataOption {
        RD_UNKNOWN,
        RD_NAIVE_AUXILIARY,
        RD_SEGMENT_FILTER,
        RD_2W_SEGMENT_FILTER,
        RD_SCORE_SEGMENT_FILTER
    };

    extern std::map<std::string, RefinementAuxiliaryDataOption> refinement_auxiliary_data_map;

    RefinementAuxiliaryDataOption get_refinement_auxiliary_data_option(const std::string &type);
    std::string get_refinement_auxiliary_data_label(const RefinementAuxiliaryDataOption opt);

}

#endif //GEFAST_OPTIONS_HPP
