/*
 * GeFaST
 *
 * Copyright (C) 2016 - 2017 Robert Mueller
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

#include "../include/Utility.hpp"


namespace GeFaST {

std::vector<std::string> readFileList(const std::string listFile) {

    std::vector<std::string> files;

    std::ifstream iStream(listFile);
    if (!iStream.good()) {
        std::cerr << "ERROR: File '" << listFile << "' not opened correctly." << std::endl;
        return files;
    }

    std::string file;

    while (std::getline(iStream, file).good()) {

        if (file.empty() || file[0] == ';') continue; // skip empty and comment lines (begin with ';')

        files.push_back(file);

    }

    return files;

}

Config<std::string> getConfiguration(int argc, const char* argv[]) {


    /* Set-up of parameter variables */
    std::map<std::string, short> parameters;
    parameters["-a"] = 1;
    parameters["--alphabet"] = 2;
    parameters["-e"] = 3;
    parameters["--extra-segments"] = 4;
    parameters["-n"] = 5;
    parameters["--name"] = 6;
    parameters["-t"] = 7;
    parameters["--threshold"] = 8;
    parameters["-c"] = 9;
    parameters["--config"] = 10;
    parameters["-i"] = 11;
    parameters["--input-files"] = 12;
    parameters["-s"] = 13;
    parameters["--seg-filter"] = 14;
    parameters["-o"] = 15;
    parameters["--output-file"] = 16;

    parameters["-si"] = 101;
    parameters["--swarm-internal"] = 102;
    parameters["-so"] = 103;
    parameters["--swarm-output"] = 104;
    parameters["-ss"] = 105;
    parameters["--swarm-statistics"] = 106;
    parameters["-sw"] = 107;
    parameters["--swarm-seeds"] = 108;
    parameters["-sn"] = 109;
    parameters["--swarm-no-otu-breaking"] = 110;
    parameters["-sd"] = 111;
    parameters["--swarm-dereplicate"] = 112;
    parameters["-sf"] = 113;
    parameters["--swarm-fastidious"] = 114;
    parameters["-sb"] = 115;
    parameters["--swarm-boundary"] = 116;
    parameters["-sm"] = 117;
    parameters["--swarm-match-reward"] = 118;
    parameters["-sp"] = 119;
    parameters["--swarm-mismatch-penalty"] = 120;
    parameters["-sg"] = 121;
    parameters["--swarm-gap-opening-penalty"] = 122;
    parameters["-se"] = 123;
    parameters["--swarm-gap-extension-penalty"] = 124;
    parameters["-sr"] = 125;
    parameters["--swarm-mothur"] = 126;
    parameters["-su"] = 127;
    parameters["--swarm-uclust"] = 128;

    parameters["--min-length"] = 1001;
    parameters["--max-length"] = 1002;
//    parameters["--per-worker"] = 1003;
//    parameters["--workers"] = 1004;
    parameters["--info-folder"] = 1005;
    parameters["--sep-abundance"] = 1006;
    parameters["--use-score"] = 1007;
    parameters["--preprocessing-only"] = 1008;

    parameters["--swarm-fastidious-checking-mode"] = 1101;
    parameters["--swarm-num-explorers"] = 1102;
    parameters["--swarm-num-grafters"] = 1103;
    parameters["--swarm-num-threads-per-check"] = 1104;
    parameters["--swarm-fastidious-threshold"] = 1105;


    std::string
            name = DEFAULT_JOB_NAME, // 5,6
            configFile = DEFAULT_CONFIG_FILE // 9,10
    ;

    unsigned long long val;
    long long signedVal;

    Config<std::string> config;
    config.set(VERSION, "0.7.0");

    /* Set default values */

    config.set(FILTER_ALPHABET, "0");
    config.set(FILTER_LENGTH, "0");
    config.set(NUM_EXTRA_SEGMENTS, "1");
//    config.set(NUM_THREADS_PER_WORKER, "1");
//    config.set(NUM_WORKERS, "1");
    config.set(PREPROCESSING_ONLY, "0");
    config.set(SEGMENT_FILTER, "0"); //TODO set "best" version of segment filter as default
    config.set(SEPARATOR_ABUNDANCE, "_");
    config.set(THRESHOLD, "1");
    config.set(USE_SCORE, "0");

    config.set(SWARM_BOUNDARY, "3");
    config.set(SWARM_DEREPLICATE, "0");
    config.set(SWARM_FASTIDIOUS, "0");
    config.set(SWARM_FASTIDIOUS_CHECKING_MODE, "0");
    config.set(SWARM_GAP_EXTENSION_PENALTY, "-4");
    config.set(SWARM_GAP_OPENING_PENALTY, "-12");
    config.set(SWARM_MATCH_REWARD, "5");
    config.set(SWARM_MISMATCH_PENALTY, "-4");
    config.set(SWARM_NO_OTU_BREAKING, "0");
    config.set(SWARM_NUM_EXPLORERS, "1");
    config.set(SWARM_NUM_GRAFTERS, "1");
    config.set(SWARM_NUM_THREADS_PER_CHECK, "1");
    config.set(SWARM_FASTIDIOUS_THRESHOLD, "0");

    /* Determine parameter values */

    // a) determine used config file and load basic configuration
    // b) determine position of first option (first arg starting with '-')
    int firstOpt = 0;
    for (int i = 1; i < argc; i++) {

        if (parameters[argv[i]] == 9 || parameters[argv[i]] == 10) {
            configFile = argv[++i];
        }

        if ((firstOpt == 0) && (argv[i][0] == '-')) {
            firstOpt = i;
        }

    }

    config.read(configFile, true);

    // read remaining configuration from command line (potentially overwriting values from the config file)
    // converting back and forth between string and integer values is used to "detect" (some) meaningless inputs
    for (int i = firstOpt; i < argc; i++) {

        // one-part parameters
        switch (parameters[argv[i]]) {
            case 109: // fallthrough -sn to --swarm-no-otu-breaking
            case 110:
                config.set(SWARM_NO_OTU_BREAKING, "1");
                continue;

            case 111: // fallthrough -sd to --swarm-dereplicate
            case 112:
                config.set(SWARM_DEREPLICATE, "1");
                continue;

            case 113: // fallthrough -sf to --swarm-fastidious
            case 114:
                config.set(SWARM_FASTIDIOUS, "1");
                continue;

            case 125: // fallthrough -sr to --swarm-mothur
            case 126:
                config.set(SWARM_MOTHUR, "1");
                continue;

            case 1007:
                config.set(USE_SCORE, "1");
                continue;

            case 1008:
                config.set(PREPROCESSING_ONLY, "1");
                continue;

            default:
                // do nothing
                break;
        }

        // two-part parameters
        if (i + 1 != argc) {
            switch (parameters[argv[i]]) {
                case 1: // fallthrough -a to --alphabet
                case 2:
                    config.set(FILTER_ALPHABET, std::to_string(true));
                    config.set(ALPHABET, argv[++i]);
                    break;

                case 3: // fallthrough -e to --extra-segments
                case 4:
                    val = std::stoul(argv[++i]);
                    config.set(NUM_EXTRA_SEGMENTS, std::to_string(val));
                    break;

                case 5: // fallthrough -n to --name
                case 6:
                    name = argv[++i]; // further handling after reading all parameters
                    break;

                case 7: // fallthrough -t to --threshold
                case 8:
                    val = std::stoul(argv[++i]);
                    config.set(THRESHOLD, std::to_string(val));
                    break;

                case 9: // fallthrough -c to --config
                case 10:
                    config.set(CONFIG_FILE, argv[++i]);
                    break;

                case 11: // fallthrough -i to --input-files
                case 12:
                    config.set(FILE_LIST, argv[++i]);
                    break;

                case 13: // fallthrough -s to --seg-filter
                case 14:
                    val = std::stoul(argv[++i]);
                    config.set(SEGMENT_FILTER, std::to_string(val));
                    break;

                case 15: // fallthrough -o to --output-file
                case 16:
                    config.set(MATCHES_OUTPUT_FILE, argv[++i]);
                    break;


                case 101: // fallthrough -si to --swarm-internal
                case 102:
                    config.set(SWARM_OUTPUT_INTERNAL, argv[++i]);
                    break;

                case 103: // fallthrough -so to --swarm-output
                case 104:
                    config.set(SWARM_OUTPUT_OTUS, argv[++i]);
                    break;

                case 105: // fallthrough -ss to --swarm-statistics
                case 106:
                    config.set(SWARM_OUTPUT_STATISTICS, argv[++i]);
                    break;

                case 107: // fallthrough -sw to --swarm-seeds
                case 108:
                    config.set(SWARM_OUTPUT_SEEDS, argv[++i]);
                    break;

                case 115: // fallthrough -sb to --swarm-boundary
                case 116:
                    val = std::stoul(argv[++i]);
                    config.set(SWARM_BOUNDARY, std::to_string(val));
                    break;

                case 117: // fallthrough -sm to --swarm-match-reward
                case 118:
                    val = std::stoul(argv[++i]);
                    config.set(SWARM_MATCH_REWARD, std::to_string(val));
                    break;

                case 119: // fallthrough -sp to --swarm-mismatch-penalty
                case 120:
                    signedVal = std::stoll(argv[++i]);
                    config.set(SWARM_MISMATCH_PENALTY, std::to_string(signedVal));
                    break;

                case 121: // fallthrough -sg to --swarm-gap-opening-penalty
                case 122:
                    signedVal = std::stoll(argv[++i]);
                    config.set(SWARM_GAP_OPENING_PENALTY, std::to_string(signedVal));
                    break;

                case 123: // fallthrough -se to --swarm-gap-extension-penalty
                case 124:
                    signedVal = std::stoll(argv[++i]);
                    config.set(SWARM_GAP_EXTENSION_PENALTY, std::to_string(signedVal));
                    break;

                case 127: // fallthrough -su to --swarm-uclust
                case 128:
                    config.set(SWARM_OUTPUT_UCLUST, argv[++i]);
                    break;

                // parameters without abbreviation
                case 1001:
                    val = std::stoul(argv[++i]); // further handling after reading all parameters
                    config.set(MIN_LENGTH, std::to_string(val));
                    break;

                case 1002:
                    val = std::stoul(argv[++i]); // further handling after reading all parameters
                    config.set(MAX_LENGTH, std::to_string(val));
                    break;

//                case 1003:
//                    val = std::stoul(argv[++i]);
//                    config.set(NUM_THREADS_PER_WORKER, std::to_string(val));
//                    break;

//                case 1004:
//                    val = std::stoul(argv[++i]);
//                    config.set(NUM_WORKERS, std::to_string(val));
//                    break;

                case 1005:
                    config.set(INFO_FOLDER, argv[++i]);
                    break;

                case 1006:
                    config.set(SEPARATOR_ABUNDANCE, argv[++i]);
                    break;

                case 1101:
                    val = std::stoul(argv[++i]);
                    config.set(SWARM_FASTIDIOUS_CHECKING_MODE, std::to_string(val));
                    break;

                case 1102:
                    val = std::stoul(argv[++i]);
                    config.set(SWARM_NUM_EXPLORERS, std::to_string(val));
                    break;

                case 1103:
                    val = std::stoul(argv[++i]);
                    config.set(SWARM_NUM_GRAFTERS, std::to_string(val));
                    break;

                case 1104:
                    val = std::stoul(argv[++i]);
                    config.set(SWARM_NUM_THREADS_PER_CHECK, std::to_string(val));
                    break;

                case 1105:
                    val = std::stoul(argv[++i]);
                    config.set(SWARM_FASTIDIOUS_THRESHOLD, std::to_string(val));
                    break;

                default:
                    std::cout << "Unknown parameter: " << argv[i] << " (is ignored)" << std::endl;
                    break;
            }
        }

    }


    // further handling of some parameters
    config.set(FILTER_LENGTH, std::to_string(
            0
            + config.peek(MAX_LENGTH)
            + 2 * config.peek(MIN_LENGTH)
    ));

    if (name == DEFAULT_JOB_NAME) { // append current time / date to default name

        time_t rawTime;
        time(&rawTime);
        std::string curTime = ctime(&rawTime);
        std::replace(curTime.begin(), curTime.end(), ' ', '_');
        name = name + "_" + curTime.substr(0, curTime.length() - 1);

    }
    config.set(NAME, name);

    return config;

}


void writeJobParameters(const std::string oFile, const Config<std::string>& conf, const std::vector<std::string>& inputFiles) {

    std::ofstream oStream(oFile);

    oStream << "### Configuration parameters (from file & command line):" << std::endl;
    conf.print(oStream);
    oStream << std::endl << std::endl;


    oStream << "### Input data files:" << std::endl;
    for (auto iter = inputFiles.begin(); iter != inputFiles.end(); iter++) {
        oStream << *iter << std::endl;
    }

    oStream.close();
}

unsigned long long gcd(unsigned long long a, unsigned long long b) {
    return (b == 0) ? a : gcd(b, a % b);
}

}