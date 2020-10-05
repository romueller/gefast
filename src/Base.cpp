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
#include <iostream>
#include <fstream>
#include <map>

#include "../include/Base.hpp"

namespace GeFaST {

    /* === Configuration === */

    void Configuration::print(std::ostream& stream) const {

        stream << "===== Configuration =====\n" << std::endl;

        stream << ">>> Input / output" << std::endl;
        stream << "Input file(s): ";
        if (!input_files.empty()) {

            stream << input_files.front();
            for (auto i = 1; i < input_files.size(); i++) stream << ", " << input_files[i];
            stream << std::endl;

        }
        stream << "Configuration file: " << configuration_file << std::endl;
        if (!output_internal.empty()) stream << "Output file (internal structures): " << output_internal << std::endl;
        if (!output_otus.empty()) stream << "Output file (clusters): " << output_otus << std::endl;
        if (!output_statistics.empty()) stream << "Output file (statistics): " << output_statistics << std::endl;
        if (!output_seeds.empty()) stream << "Output file (seeds): " << output_seeds << std::endl;
        if (!output_uclust.empty()) stream << "Output file (UCLUST): " << output_uclust << std::endl;
        if (!keep_preprocessed.empty()) stream << "Output file (preprocessing results): " << keep_preprocessed << std::endl;
        stream << "Output generator: " << get_output_generator_label(opt_output_generator) << std::endl;
        stream << "Abundance separator: " << separator << std::endl;
        stream << "Mothur-compatible output format: " << mothur << std::endl;
        stream << std::endl;

        stream << ">>> Preprocessing" << std::endl;
        stream << "Preprocessor type: " << get_preprocessor_label(opt_preprocessor) << std::endl;
        if (filter_alphabet) stream << "Filter by alphabet: " << alphabet << std::endl;
        if (filter_length) stream << "Filter by length: " << filter_alphabet << std::endl;
        if (min_length > 0) stream << "   Minimum length: " << min_length << std::endl;
        if (max_length > 0) stream << "   Maximum length: " << max_length << std::endl;
        if (filter_abundance) stream << "Filter by abundance: " << filter_abundance << std::endl;
        if (min_abundance > 0) stream << "   Minimum abundance: " << min_abundance << std::endl;
        if (max_abundance > 0) stream << "   Maximum abundance: " << max_abundance << std::endl;
        stream << std::endl;

        stream << ">>> Swarm clustering" << std::endl;
        stream << "Mode: " << mode << std::endl;
        stream << "Clusterer: " << get_clusterer_label(opt_clusterer) << std::endl;
        if (refinement_threshold > 0) stream << "Cluster refiner: " << get_refiner_label(opt_refiner) << std::endl;
        stream << "Distance function: " << get_distance_label(opt_distance) << std::endl;
        stream << "Threshold: " << main_threshold << std::endl;
        if (refinement_threshold > 0) stream << "Refinement threshold: " << refinement_threshold << std::endl;
        if (iterative_refinement_thresholds.size() > 1) {

            stream << "Iterative refinement thresholds: " << iterative_refinement_thresholds.front();
            for (auto i = 1; i < iterative_refinement_thresholds.size(); i++) stream << ", " << iterative_refinement_thresholds[i];
            stream << std::endl;

        }
        if (refinement_threshold > 0) stream << "Mass boundary: " << boundary << std::endl;
        stream << "Swarm breaking: " << break_swarms << std::endl;
        stream << std::endl;

        if (!misc.empty()) {

            stream << ">>> Configuration parameters provided via 'misc'" << std::endl;
            for (auto iter = misc.begin(); iter != misc.end(); iter++) {
                stream << iter->first << ": " << iter->second << std::endl;
            }
            stream << std::endl;

        }

    }

    // helper method for parse_and_set_refinement_thresholds(...) to translate relative threshold
    dist_t get_absolute_threshold(const std::string& str, const dist_t base_threshold) {

        dist_t res = base_threshold;

        switch (str[0]) {

            case '+': {
                res += std::stof(str.substr(1));
                break;
            }

            case '-': {
                res -= std::stof(str.substr(1));
                break;
            }

            case '*': {
                res *= std::stof(str.substr(1));
                break;
            }

            case '/': {
                res /= std::stof(str.substr(1));
                break;
            }

            default: {
                res = std::stof(str.substr(0));
                break;
            }

        }

        return res;

    }

    void Configuration::parse_and_set_refinement_thresholds(std::string arg) {

        if (arg.find(':') != std::string::npos) { // range notation <from>:<to>:<increment>

            size_t pos = arg.find(':');
            dist_t first = get_absolute_threshold(arg.substr(0, pos), main_threshold);
            arg.erase(0, pos + 1);
            pos = arg.find(':');
            dist_t last = get_absolute_threshold(arg.substr(0, pos), main_threshold);
            arg.erase(0, pos + 1);
            dist_t inc = std::stof(arg);

            for (dist_t t = first; t <= last + inc / 2.0; t += inc) {
                iterative_refinement_thresholds.push_back(t);
            }

            refinement_threshold = iterative_refinement_thresholds.back();

        } else { // list notation <t>[,<t>]*

            size_t pos = 0;
            do {

                pos = arg.find(',');
                iterative_refinement_thresholds.push_back(get_absolute_threshold(arg.substr(0, pos), main_threshold));
                arg.erase(0, pos + 1);

            } while (pos != std::string::npos);

            refinement_threshold = iterative_refinement_thresholds.back();

        }

    }

    void Configuration::parse_misc(std::string arg) {

        while (arg.length() > 0) {

            std::string kv = arg.substr(0, arg.find('$'));
            size_t pos = arg.find(':');
            misc[kv.substr(0, pos)] = kv.substr(pos + 1);
            arg.erase(0, kv.length() + 1);

        }

    }

    void Configuration::set_general_parameters(int argc, const char* argv[]) {

        bool read_from_list_file = false;
        mode = argv[1];
        std::string input_parameter = argv[2];
        configuration_file = argv[3];


        /* Set default values */

        alphabet = "ACGT";
        main_threshold = 0;
        input_files = std::vector<std::string>();
        output_internal = "";
        output_otus = "";
        output_statistics = "";
        output_seeds = "";
        output_uclust = "";
        break_swarms = true;
        refinement_threshold = 0;
        iterative_refinement_thresholds.resize(0);
        boundary = 3;
        mothur = false;
        min_length = 0;
        max_length = 0;
        min_abundance = 0;
        max_abundance = 0;
        separator = '_';
        filter_alphabet = true;
        filter_length = false;
        filter_abundance = false;
        opt_preprocessor = PP_FASTA;
        opt_quality_encoding = QE_SANGER;
        opt_clusterer = CL_IDLE;
        opt_refiner = CR_IDLE;
        opt_output_generator = OG_IDLE;
        opt_amplicon_collection = AC_UNKNOWN;
        opt_distance = DT_UNKNOWN;
        match_reward = 5;
        mismatch_penalty = -4;
        gap_opening_penalty = -12;
        gap_extension_penalty = -4;
        keep_preprocessed = "";


        /* Add basic configuration from file */

        std::ifstream in_stream(configuration_file);
        std::string line;
        unsigned long delim_pos;
        std::string refinement_arg;

        if (!in_stream.good()) {
            std::cerr << "WARNING: Configuration file '" << configuration_file << "' not opened correctly. No parameters are read from it." << std::endl;
        }

        while (std::getline(in_stream, line).good()) {

            if (line.empty() || line.front() == '#') continue;

            delim_pos = line.find('=');
            std::string key = line.substr(0, delim_pos);

            if (key == "alphabet") alphabet = line.substr(delim_pos + 1);
            if (key == "threshold") main_threshold = std::stof(line.substr(delim_pos + 1));
            if (key == "output_internal") output_internal = line.substr(delim_pos + 1);
            if (key == "output_otus") output_otus = line.substr(delim_pos + 1);
            if (key == "output_statistics") output_statistics = line.substr(delim_pos + 1);
            if (key == "output_seeds") output_seeds = line.substr(delim_pos + 1);
            if (key == "output_uclust") output_uclust = line.substr(delim_pos + 1);
            if (key == "break_swarms") break_swarms = (line.substr(delim_pos + 1) == "1");
            if (key == "refinement_threshold") refinement_arg = line.substr(delim_pos + 1);
            if (key == "boundary") boundary = std::stoul(line.substr(delim_pos + 1));
            if (key == "mothur") mothur = (line.substr(delim_pos + 1) == "1");
            if (key == "min_length") min_length = std::stoul(line.substr(delim_pos + 1));
            if (key == "max_length") max_length = std::stoul(line.substr(delim_pos + 1));
            if (key == "min_abundance") min_abundance = std::stoul(line.substr(delim_pos + 1));
            if (key == "max_abundance") max_abundance = std::stoul(line.substr(delim_pos + 1));
            if (key == "sep_abundance") separator = line.substr(delim_pos + 1);
            if (key == "preprocessor") opt_preprocessor = get_preprocessor_option(line.substr(delim_pos + 1));
            if (key == "clusterer") opt_clusterer = get_clusterer_option(line.substr(delim_pos + 1));
            if (key == "refiner") opt_refiner = get_refiner_option(line.substr(delim_pos + 1));
            if (key == "output_generator") opt_output_generator = get_output_generator_option(line.substr(delim_pos + 1));
            if (key == "quality_encoding") opt_quality_encoding = get_quality_encoding_option(line.substr(delim_pos + 1));
            if (key == "match_reward") match_reward = std::stoi(line.substr(delim_pos + 1));
            if (key == "mismatch_penalty") mismatch_penalty = std::stoi(line.substr(delim_pos + 1));
            if (key == "gap_opening_penalty") gap_opening_penalty = std::stoi(line.substr(delim_pos + 1));
            if (key == "gap_extension_penalty") gap_extension_penalty = std::stoi(line.substr(delim_pos + 1));
            if (key == "keep_preprocessed") keep_preprocessed = line.substr(delim_pos + 1);
            if (key == "misc") parse_misc(line.substr(delim_pos + 1));

        }

        in_stream.close();


        /* Add configuration from command line */

        // set-up of parameter variables
        std::map<std::string, short> parameter_keys;

        // category 1: parameters with abbreviation and value
        parameter_keys["-i"] = 1;  parameter_keys["--output_internal"] = 2;
        parameter_keys["-o"] = 3;  parameter_keys["--output_otus"] = 4;
        parameter_keys["-s"] = 5;  parameter_keys["--output_statistics"] = 6;
        parameter_keys["-w"] = 7;  parameter_keys["--output_seeds"] = 8;
        parameter_keys["-u"] = 9;  parameter_keys["--output_uclust"] = 10;
        parameter_keys["-t"] = 11; parameter_keys["--threshold"] = 12;
        parameter_keys["-r"] = 13; parameter_keys["--refinement_threshold"] = 14;
        parameter_keys["-b"] = 15; parameter_keys["--boundary"] = 16;
        parameter_keys["-a"] = 17; parameter_keys["--alphabet"] = 18;
        parameter_keys["-m"] = 19; parameter_keys["--match_reward"] = 20;
        parameter_keys["-p"] = 21; parameter_keys["--mismatch_penalty"] = 22;
        parameter_keys["-g"] = 23; parameter_keys["--gap_opening_penalty"] = 24;
        parameter_keys["-e"] = 25; parameter_keys["--gap_extension_penalty"] = 26;
        parameter_keys["-k"] = 27; parameter_keys["--keep_preprocessed"] = 28;
        parameter_keys["-n"] = 29; parameter_keys["--break_swarms"] = 30;



        // category 2: parameters with value but without abbreviation
        parameter_keys["--min_length"] = 101;
        parameter_keys["--max_length"] = 102;
        parameter_keys["--sep_abundance"] = 103;
        parameter_keys["--preprocessor"] = 104;
        parameter_keys["--clusterer"] = 105;
        parameter_keys["--refiner"] = 106;
        parameter_keys["--output_generator"] = 107;
        parameter_keys["--quality_encoding"] = 108;
        parameter_keys["--min_abundance"] = 109;
        parameter_keys["--max_abundance"] = 110;

        // category 3: parameters with abbreviation but without value
//        parameter_keys["-"] = 1001; parameter_keys["--"] = 1002;

        // category 4: parameters without abbreviation and value
        parameter_keys["--list_file"] = 1101;
        parameter_keys["--mothur"] = 1102;


        // read remaining configuration from command line (potentially overwriting values from the config file)
        // first arguments / parameters are fixed (program name, mode, input file list and configuration file)
        // -> start at index 4
        for (int i = 4; i < argc; i++) {

            // one-part parameters (categories 3 and 4)
            switch (parameter_keys[argv[i]]) {

//                case 1001: // fallthrough - to --
//                case 1002:
//                    ...
//                    continue;


                case 1101:
                    read_from_list_file = true;
                    continue;

                case 1102:
                    mothur = true;
                    continue;

                default:
                    // do nothing
                    break;

            }

            // two-part parameters (categories 1 and 2)
            if (i + 1 != argc) {
                switch (parameter_keys[argv[i]]) {

                    case 1: // fallthrough -i to --output_internal
                    case 2:
                        output_internal = argv[++i];
                        break;

                    case 3: // fallthrough -o to --output_otus
                    case 4:
                        output_otus = argv[++i];
                        break;

                    case 5: // fallthrough -s to --output_statistics
                    case 6:
                        output_statistics = argv[++i];
                        break;

                    case 7: // fallthrough -w to --output_seeds
                    case 8:
                        output_seeds = argv[++i];
                        break;

                    case 9: // fallthrough -u to --uclust
                    case 10:
                        output_uclust = argv[++i];
                        break;

                    case 11: // fallthrough -t to --threshold
                    case 12:
                        main_threshold = std::stof(argv[++i]);
                        break;

                    case 13: // fallthrough -r to --refinement_threshold
                    case 14:
                        refinement_arg = argv[++i];
                        break;

                    case 15: // fallthrough -b to --boundary
                    case 16:
                        boundary = std::stoul(argv[++i]);
                        break;

                    case 17: // fallthrough -a to --alphabet
                    case 18:
                        alphabet = argv[++i];
                        break;

                    case 19: // fallthrough -m to --match_reward
                    case 20:
                        match_reward = std::stoi(argv[++i]);
                        break;

                    case 21: // fallthrough -p to --mismatch_penalty
                    case 22:
                        mismatch_penalty = std::stoi(argv[++i]);
                        break;

                    case 23: // fallthrough -g to --gap_opening_penalty
                    case 24:
                        gap_opening_penalty = std::stoi(argv[++i]);
                        break;

                    case 25: // fallthrough -e to --gap_extension_penalty
                    case 26:
                        gap_extension_penalty = std::stoi(argv[++i]);
                        break;

                    case 27: // fallthrough -k to --keep_preprocessed
                    case 28:
                        keep_preprocessed = argv[++i];
                        break;

                    case 29: // fallthrough -n to --break_swarms
                    case 30:
                        break_swarms = std::stoi(argv[++i]) == 1;
                        continue;


                    case 101:
                        min_length = std::stoul(argv[++i]);
                        break;

                    case 102:
                        max_length = std::stoul(argv[++i]);
                        break;

                    case 103:
                        separator = argv[++i];
                        break;

                    case 104:
                        opt_preprocessor = get_preprocessor_option(argv[++i]);
                        break;

                    case 105:
                        opt_clusterer = get_clusterer_option(argv[++i]);
                        break;

                    case 106:
                        opt_refiner = get_refiner_option(argv[++i]);
                        break;

                    case 107:
                        opt_output_generator = get_output_generator_option(argv[++i]);
                        break;

                    case 108:
                        opt_quality_encoding = get_quality_encoding_option(argv[++i]);
                        break;

                    case 109:
                        min_abundance = std::stoul(argv[++i]);
                        break;

                    case 110:
                        max_abundance = std::stoul(argv[++i]);
                        break;

                    default:
                        std::cout << "Unknown parameter: " << argv[i] << " (is ignored)" << std::endl;
                        break;

                }

            }

        }

        // further handling of some parameters
        filter_alphabet = !alphabet.empty();
        filter_length = 0 + (max_length != 0) + 2 * (min_length != 0);
        filter_abundance = 0 + (max_abundance != 0) + 2 * (min_abundance != 0);

        if (read_from_list_file) {
            input_files = read_list_file(input_parameter);
        } else {

            size_t pos = 0;
            do {

                pos = input_parameter.find(',');
                input_files.push_back(input_parameter.substr(0, pos));
                input_parameter.erase(0, pos + 1);

            } while (pos != std::string::npos);

        }

        // refinement thresholds are handled after reading the configuration file and command-line arguments
        // to ensure that the threshold parameter has already been set per configuration file or command line
        if (!refinement_arg.empty()) {
            parse_and_set_refinement_thresholds(refinement_arg);
        }

    }


    /* === AmpliconCollection === */

    unsigned long AmpliconCollection::qgram_diff(const numSeqs_t i, const numSeqs_t j) const {
        return 0; // return difference of 0 when q-gram information are not available
    }


    /* === QgramAmpliconCollection === */

    char acgtu_map[256] =
            {
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, 4, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, 4, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
            };

    long popcnt_present;

    long QGRAM_LENGTH = 5;
    long QGRAM_VECTOR_BYTES = (1 << (2 * QGRAM_LENGTH)) / 8;

    unsigned long QgramComparer::popcount(unsigned long x) const {

        unsigned long y;
        popcnt_asm(x, y);
        return y;

    }

    void QgramComparer::cpu_features_detect() {

        unsigned int a, b, c, d;

        cpuid(0, 0, a, b, c, d);
        unsigned int maxlevel = a & 0xff;

        if (maxlevel >= 1) {

            cpuid(1, 0, a, b, c, d);
            popcnt_present = (c >> 23) & 1;

        }

    }

    unsigned long QgramComparer::popcount_128(__m128i x) const {

        __m128i mask1 = _mm_set_epi8(0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55,
                                     0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55);

        __m128i mask2 = _mm_set_epi8(0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33,
                                     0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33);

        __m128i mask4 = _mm_set_epi8(0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f,
                                     0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f);

        __m128i zero = _mm_setzero_si128();

        /* add together 2 bits: 0+1, 2+3, 3+4, ... 126+127 */

        __m128i a = _mm_srli_epi64(x, 1);
        __m128i b = _mm_and_si128(x, mask1);
        __m128i c = _mm_and_si128(a, mask1);
        __m128i d = _mm_add_epi64(b, c);

        /* add together 4 bits: (0+1)+(2+3), ... (124+125)+(126+127) */

        __m128i e = _mm_srli_epi64(d, 2);
        __m128i f = _mm_and_si128(d, mask2);
        __m128i g = _mm_and_si128(e, mask2);
        __m128i h = _mm_add_epi64(f, g);

        /* add together 8 bits: (0..3)+(4..7), ... (120..123)+(124..127) */

        __m128i i = _mm_srli_epi64(h, 4);
        __m128i j = _mm_add_epi64(h, i);
        __m128i k = _mm_and_si128(j, mask4);

        /* add together 8 bytes: (0..63) and (64..127) */

        __m128i l = _mm_sad_epu8(k, zero);

        /* add together 64-bit values into final 128 bit value */

        __m128i m = _mm_srli_si128(l, 8);
        __m128i n = _mm_add_epi64(m, l);

        /* return low 64 bits: return value is always in range 0 to 128 */

        unsigned long o = (unsigned long) _mm_movepi64_pi64(n);

        return o;

    }

    unsigned long QgramComparer::compare_qgram_vectors_128(const unsigned char* a, const unsigned char* b) const {

        /* Count number of different bits */
        /* Uses SSE2 but not POPCNT instruction */
        /* input MUST be 16-byte aligned */

        __m128i* ap = (__m128i*) a;
        __m128i* bp = (__m128i*) b;
        unsigned long count = 0;

        while ((unsigned char*)ap < a + QGRAM_VECTOR_BYTES)
            count += popcount_128(_mm_xor_si128(*ap++, *bp++));

        return count;

    }

    unsigned long QgramComparer::compare_qgram_vectors_64(const unsigned char* a, const unsigned char* b) const {

        /* Count number of different bits */
        /* Uses the POPCNT instruction, requires CPU with this feature */

        unsigned long* ap = (unsigned long*)a;
        unsigned long* bp = (unsigned long*)b;
        unsigned long count = 0;

        while ((unsigned char*) ap < a + QGRAM_VECTOR_BYTES)
            count += popcount(*ap++ ^ *bp++);

        return count;

    }

    unsigned long QgramComparer::compare_qgram_vectors(const unsigned char* a, const unsigned char* b) const {

        if (popcnt_present)
            return compare_qgram_vectors_64(a,b);
        else
            return compare_qgram_vectors_128(a,b);

    }


    /* === SimpleGraftingInfo === */

    SimpleGraftingInfo* SimpleGraftingInfo::clone() const {
        return new SimpleGraftingInfo(*this);
    }

    numSeqs_t SimpleGraftingInfo::get_parent(const numSeqs_t i) const {
        return parents_[i];
    }

    numSeqs_t SimpleGraftingInfo::get_child(const numSeqs_t i) const {
        return children_[i];
    }

    Swarm* SimpleGraftingInfo::get_child_swarm(const numSeqs_t i) const {
        return child_swarms_[i];
    }

    dist_t SimpleGraftingInfo::get_distance(const numSeqs_t i) const {
        return distances_[i];
    }

    void SimpleGraftingInfo::add_link(const numSeqs_t p, const numSeqs_t c, Swarm* cs, const dist_t d) {

        parents_.push_back(p);
        children_.push_back(c);
        child_swarms_.push_back(cs);
        distances_.push_back(d);

    }

    numSeqs_t SimpleGraftingInfo::num_links() const {
        return parents_.size();
    }


    /* === SwarmOrderer === */

    std::vector<std::pair<const AmpliconCollection*, const Swarm*>> SwarmOrderer::order_by_mass(
            const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage) {

        std::vector<std::pair<const AmpliconCollection*, const Swarm*>> order;
        for (numSeqs_t pid = 0; pid < swarm_storage.num_pools(); pid++) {

            auto& swarms_in_pool = swarm_storage.get_swarms(pid);
            for (numSeqs_t sid = 0; sid < swarms_in_pool.size(); sid++) {
                order.emplace_back(&(amplicon_storage.get_pool(pid)), &(swarms_in_pool.get(sid)));
            }

        }

        std::sort(order.begin(), order.end(),
                  [](const std::pair<const AmpliconCollection *, const Swarm *> &a,
                     const std::pair<const AmpliconCollection *, const Swarm *> &b) {
                      return (a.second->total_mass() > b.second->total_mass())
                             || ((a.second->total_mass() == b.second->total_mass())
                                 && ((a.first->ab(a.second->seed()) > b.first->ab(b.second->seed()))
                                     || ((a.first->ab(a.second->seed()) == b.first->ab(b.second->seed()))
                                         && (strcmp(a.first->id(a.second->seed()), b.first->id(b.second->seed())) < 0)))
                             );
                  });

        return order;

    }

    std::vector<std::pair<const AmpliconCollection*, const Swarm*>> SwarmOrderer::order_by_seed_abundance(
            const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage) {

        std::vector<std::pair<const AmpliconCollection*, const Swarm*>> order;
        for (numSeqs_t pid = 0; pid < swarm_storage.num_pools(); pid++) {
            auto& swarms_in_pool = swarm_storage.get_swarms(pid);

            for (numSeqs_t sid = 0; sid < swarms_in_pool.size(); sid++) {
                order.emplace_back(&(amplicon_storage.get_pool(pid)), &(swarms_in_pool.get(sid)));
            }

        }

        std::sort(order.begin(), order.end(),
                  [](const std::pair<const AmpliconCollection *, const Swarm *> &a,
                     const std::pair<const AmpliconCollection *, const Swarm *> &b) {
                      return (a.first->ab(a.second->seed()) > b.first->ab(b.second->seed()))
                             || ((a.first->ab(a.second->seed()) > b.first->ab(b.second->seed()))
                                 && (strcmp(a.first->id(a.second->seed()), b.first->id(b.second->seed())) < 0)
                             );
                  });

        return order;

    }

}