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

#include <algorithm>

#include "../include/Options.hpp"

namespace GeFaST {

    /* === Configuration === */

    std::map<std::string, ConfigurationOption> config_map = {
            {"prep", CF_PREPROCESSING_ONLY},
            {"derep", CF_DEREPLICATION},
            {"lev", CF_LEVENSHTEIN},
            {"qlev", CF_QUALITY_LEVENSHTEIN},
            {"as", CF_ALIGNMENT_SCORE},
            {"qas", CF_QUALITY_ALIGNMENT_SCORE},
            {"cons", CF_CONSISTENCY}
    };

    ConfigurationOption get_configuration_option(const std::string& type) {
        return config_map[type];
    }

    std::string get_configuration_label(const ConfigurationOption opt) {
        return std::find_if(config_map.begin(), config_map.end(),
                [opt](const std::map<std::string, ConfigurationOption>::value_type& p) {return p.second == opt;})->first;
    }

    /* === AmpliconStorage === */

    std::map<std::string, AmpliconStorageOption> ampl_storage_map = {
            {"simple-length", AS_SIMPLE_LENGTH_POOLS},
            {"prepared-length", AS_PREPARED_LENGTH_POOLS},
            {"prepared-length-quality", AS_PREPARED_LENGTH_POOLS_QUALITY},
            {"prepared", AS_PREPARED_POOLS},
            {"prepared-quality", AS_PREPARED_POOLS_QUALITY},
            {"prepared-qgram-length", AS_PREPARED_QGRAM_LENGTH_POOLS},
            {"prepared-qgram-length-quality", AS_PREPARED_QGRAM_LENGTH_POOLS_QUALITY}
    };

    AmpliconStorageOption get_amplicon_storage_option(const std::string& type) {
        return ampl_storage_map[type];
    }

    std::string get_amplicon_storage_label(const AmpliconStorageOption opt) {
        return std::find_if(ampl_storage_map.begin(), ampl_storage_map.end(),
                [opt](const std::map<std::string, AmpliconStorageOption>::value_type& p) {return p.second == opt;})->first;
    }

    /* === AmpliconCollection === */

    std::map<std::string, AmpliconCollectionOption> ampl_collection_map = {
            {"simple", AC_SIMPLE},
            {"array", AC_ARRAY},
            {"array-quality", AC_ARRAY_QUALITY},
            {"array-qgram", AC_ARRAY_QGRAM},
            {"array-qgram-quality", AC_ARRAY_QGRAM_QUALITY}
    };

    AmpliconCollectionOption get_amplicon_collection_option(const std::string& type) {
        return ampl_collection_map[type];
    }

    std::string get_amplicon_collection_label(const AmpliconCollectionOption opt) {
        return std::find_if(ampl_collection_map.begin(), ampl_collection_map.end(),
                [opt](const std::map<std::string, AmpliconCollectionOption>::value_type& p) {return p.second == opt;})->first;
    }


    /* === Distance === */

    std::map<std::string, DistanceOption> distance_map = {
            {"none", DT_UNKNOWN},
            {"bounded-levenshtein", DT_BOUNDED_LEVENSHTEIN},
            {"bounded-score-levenshtein", DT_BOUNDED_SCORE_LEVENSHTEIN},
            {"bounded-score", DT_BOUNDED_SCORE},
            {"bounded-banded-score", DT_BOUNDED_BANDED_SCORE},
            {"clement", DT_QW_CLEMENT},
            {"malde-a", DT_QW_MALDE_A},
            {"malde-b", DT_QW_MALDE_B},
            {"malde-c", DT_QW_MALDE_C},
            {"frith", DT_QW_FRITH},
            {"kim-a", DT_QW_KIM_A},
            {"kim-b", DT_QW_KIM_B},
            {"converge-a", DT_QW_CONVERGE_A},
            {"converge-b", DT_QW_CONVERGE_B},
            {"score-clement", DT_QWS_CLEMENT},
            {"score-malde-a", DT_QWS_MALDE_A},
            {"score-malde-b", DT_QWS_MALDE_B},
            {"score-malde-c", DT_QWS_MALDE_C},
            {"score-frith", DT_QWS_FRITH},
            {"score-kim-a", DT_QWS_KIM_A},
            {"score-kim-b", DT_QWS_KIM_B},
            {"score-converge-a", DT_QWS_CONVERGE_A},
            {"score-converge-b", DT_QWS_CONVERGE_B},
            {"banded-score-clement", DT_QWBS_CLEMENT},
            {"banded-score-malde-a", DT_QWBS_MALDE_A},
            {"banded-score-malde-b", DT_QWBS_MALDE_B},
            {"banded-score-malde-c", DT_QWBS_MALDE_C},
            {"banded-score-frith", DT_QWBS_FRITH},
            {"banded-score-kim-a", DT_QWBS_KIM_A},
            {"banded-score-kim-b", DT_QWBS_KIM_B},
            {"banded-score-converge-a", DT_QWBS_CONVERGE_A},
            {"banded-score-converge-b", DT_QWBS_CONVERGE_B}
    };

    DistanceOption get_distance_option(const std::string &type) {
        return distance_map[type];
    }

    std::string get_distance_label(const DistanceOption opt) {
        return std::find_if(distance_map.begin(), distance_map.end(),
                [opt](const std::map<std::string, DistanceOption>::value_type& p) {return p.second == opt;})->first;
    }


    /* === Boosting === */

    std::map<std::string, BoostingOption> boosting_map = {
            {"none", BO_NO_EFFECT},
            {"linear", BO_LINEAR},
            {"mult", BO_MULT},
            {"root", BO_ROOT}
    };

    BoostingOption get_boosting_option(const std::string& type) {
        return boosting_map[type];
    }

    std::string get_boosting_label(const BoostingOption opt) {
        return std::find_if(boosting_map.begin(), boosting_map.end(),
                [opt](const std::map<std::string, BoostingOption>::value_type& p) {return p.second == opt;})->first;
    }


    /* === Preprocessor === */

    std::map<std::string, PreprocessorOption> prep_map = {
            {"fasta", PP_FASTA},
            {"fastq", PP_FASTQ},
    };

    PreprocessorOption get_preprocessor_option(const std::string& type) {
        return prep_map[type];
    }

    std::string get_preprocessor_label(const PreprocessorOption opt) {
        return std::find_if(prep_map.begin(), prep_map.end(),
                [opt](const std::map<std::string, PreprocessorOption>::value_type& p) {return p.second == opt;})->first;
    }


    /* === QualityEncoding === */

    std::map<std::string, QualityEncodingOption> enc_map = {
            {"sanger", QE_SANGER},
            {"illumina1.3", QE_ILLUMINA13},
            {"illumina1.5", QE_ILLUMINA15},
            {"illumina1.8", QE_ILLUMINA18}
    };

    QualityEncodingOption get_quality_encoding_option(const std::string& type) {
        return enc_map[type];
    }

    std::string get_quality_encoding_label(const QualityEncodingOption opt) {
        return std::find_if(enc_map.begin(), enc_map.end(),
                [opt](const std::map<std::string, QualityEncodingOption>::value_type& p) {return p.second == opt;})->first;
    }


    /* === Clusterer === */

    std::map<std::string, ClustererOption> clusterer_map = {
            {"idle", CL_IDLE},
            {"classic", CL_CLASSIC_SWARMER},
            {"dereplicator", CL_DEREPLICATOR},
            {"dada2", CL_DADA2},
            {"cons-classic", CL_CONS_CLASSIC},
            {"cons-swarmer", CL_CONS_SWARMER}
    };

    ClustererOption get_clusterer_option(const std::string& type) {
        return clusterer_map[type];
    }

    std::string get_clusterer_label(const ClustererOption opt) {
        return std::find_if(clusterer_map.begin(), clusterer_map.end(),
                [opt](const std::map<std::string, ClustererOption>::value_type& p) {return p.second == opt;})->first;
    }


    /* === ClusterRefiner === */

    std::map<std::string, ClusterRefinerOption> refiner_map = {
            {"idle", CR_IDLE},
            {"classic", CR_CLASSIC_FASTIDIOUS},
            {"iterative", CR_ITERATIVE_FASTIDIOUS},
            {"cons-lsa", CR_CONS_LSA},
            {"cons-lsr", CR_CONS_LSR},
            {"cons-lss", CR_CONS_LSS},
            {"cons-ss", CR_CONS_SS}
    };

    ClusterRefinerOption get_refiner_option(const std::string& type) {
        return refiner_map[type];
    }

    std::string get_refiner_label(const ClusterRefinerOption opt) {
        return std::find_if(refiner_map.begin(), refiner_map.end(),
                [opt](const std::map<std::string, ClusterRefinerOption>::value_type& p) {return p.second == opt;})->first;
    }


    /* === OutputGenerator === */

    std::map<std::string, OutputGeneratorOption> output_generator_map = {
            {"idle", OG_IDLE},
            {"classic", OG_CLASSIC},
            {"dereplication", OG_DEREPLICATION}
    };

    OutputGeneratorOption get_output_generator_option(const std::string& type) {
        return output_generator_map[type];
    }

    std::string get_output_generator_label(const OutputGeneratorOption opt) {
        return std::find_if(output_generator_map.begin(), output_generator_map.end(),
                [opt](const std::map<std::string, OutputGeneratorOption>::value_type& p) {return p.second == opt;})->first;
    }


    /* === SwarmStorage === */

    std::map<std::string, SwarmStorageOption> swarm_storage_map = {
            {"simple", SS_SIMPLE_PER_POOL}
    };

    SwarmStorageOption get_swarm_storage_option(const std::string& type) {
        return swarm_storage_map[type];
    }

    std::string get_swarm_storage_label(const SwarmStorageOption opt) {
        return std::find_if(swarm_storage_map.begin(), swarm_storage_map.end(),
                [opt](const std::map<std::string, SwarmStorageOption>::value_type& p) {return p.second == opt;})->first;
    }


    /* === AuxiliaryData === */

    std::map<std::string, AuxiliaryDataOption> auxiliary_data_map = {
            {"naive", AD_NAIVE_AUXILIARY},
            {"one-way", AD_SEGMENT_FILTER},
            {"two-way", AD_2W_SEGMENT_FILTER},
            {"score", AD_SCORE_SEGMENT_FILTER}
    };

    AuxiliaryDataOption get_auxiliary_data_option(const std::string& type) {
        return auxiliary_data_map[type];
    }

    std::string get_auxiliary_data_label(const AuxiliaryDataOption opt) {
        return std::find_if(auxiliary_data_map.begin(), auxiliary_data_map.end(),
                [opt](const std::map<std::string, AuxiliaryDataOption>::value_type& p) {return p.second == opt;})->first;
    }


    /* === RefinementAuxiliaryData === */

    std::map<std::string, RefinementAuxiliaryDataOption> refinement_auxiliary_data_map = {
            {"naive", RD_NAIVE_AUXILIARY},
            {"one-way", RD_SEGMENT_FILTER},
            {"two-way", RD_2W_SEGMENT_FILTER},
            {"score", RD_SCORE_SEGMENT_FILTER}
    };

    RefinementAuxiliaryDataOption get_refinement_auxiliary_data_option(const std::string &type) {
        return refinement_auxiliary_data_map[type];
    }

    std::string get_refinement_auxiliary_data_label(const RefinementAuxiliaryDataOption opt) {
        return std::find_if(refinement_auxiliary_data_map.begin(), refinement_auxiliary_data_map.end(),
                [opt](const std::map<std::string, RefinementAuxiliaryDataOption>::value_type& p) {return p.second == opt;})->first;
    }

}