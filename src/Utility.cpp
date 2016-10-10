/*
 * SCT-PJ
 *
 * Copyright (C) 2016 Robert Mueller
 *
 * TODO add licence text (e.g. GNU (Affero) General Public License)
 *
 * Contact: Robert Mueller <romueller@techfak.uni-bielefeld.de>
 * Faculty of Technology, Bielefeld University,
 * PO box 100131, DE-33501 Bielefeld, Germany
 */

#include "../include/Utility.hpp"


namespace SCT_PJ {

std::vector<std::string> readFileList(const std::string listFile) {

    std::vector<std::string> files;

    std::ifstream iStream(listFile);
    if (!iStream.good()) {
        std::cerr << "ERROR: File '" << listFile << "' not opened correctly." << std::endl;
        return files;
    }

    std::string file;

    while (std::getline(iStream, file).good()) {

        if (file.empty() || file[0] == ';') continue; //skip empty and comment lines (begin with ';')

        files.push_back(file);

    }

    return files;

}

Config<std::string> getConfiguration(int argc, const char* argv[]) {


    /* Set-up of parameter variables */
    std::map<std::string, short> parameters;
    parameters["-a"] = 1;
    parameters["--alphabet"] = 2; // optional
    parameters["-e"] = 3;
    parameters["--extra-segments"] = 4; // optional (has default)
    parameters["-n"] = 5;
    parameters["--name"] = 6; // optional (has default)
    parameters["-t"] = 7;
    parameters["--threshold"] = 8; // optional (has default)
    parameters["-c"] = 9;
    parameters["--config"] = 10; // optional (has default)
    parameters["-i"] = 11;
    parameters["--input-files"] = 12; // (optional)
    parameters["-s"] = 13;
    parameters["--seg-filter"] = 14; // optional (has default)
    parameters["-o"] = 15;
    parameters["--output-file"] = 16; // (optional)

    parameters["-si"] = 101;
    parameters["--swarm-internal"] = 102; // (optional)
    parameters["-so"] = 103;
    parameters["--swarm-output"] = 104; // (optional)
    parameters["-ss"] = 105;
    parameters["--swarm-statistics"] = 106; // (optional)
    parameters["-sw"] = 107;
    parameters["--swarm-seeds"] = 108; // (optional)
    parameters["-sn"] = 109;
    parameters["--swarm-no-otu-breaking"] = 110; // (optional)
    parameters["-sd"] = 111;
    parameters["--swarm-dereplicate"] = 112; // optional
    parameters["-sf"] = 113;
    parameters["--swarm-fastidious"] = 114; // optional
    parameters["-sb"] = 115;
    parameters["--swarm-boundary"] = 116; // optional

    parameters["--min-length"] = 1001; // optional
    parameters["--max-length"] = 1002; // optional
    parameters["--per-worker"] = 1003; // optional (has default)
    parameters["--workers"] = 1004; // optional (has default)
    parameters["--info-folder"] = 1005; // optional
    parameters["--sep-abundance"] = 1006; // optional (has default)

    parameters["--swarm-fastidious-checking-mode"] = 1101; // optional (has default)
    parameters["--swarm-num-explorers"] = 1102; // optional (has default)
    parameters["--swarm-num-grafters"] = 1103; // optional (has default)
    parameters["--swarm-num-threads-per-check"] = 1104; // optional (has default)


    std::string
            name = DEFAULT_JOB_NAME, // 5,6
            configFile = DEFAULT_CONFIG_FILE // 9,10
    ;

    unsigned long long val;

    Config<std::string> config;
    config.set(VERSION, "1.0");

    /* Set default values */

    config.set(FILTER_ALPHABET, "0");
    config.set(FILTER_LENGTH, "0");
    config.set(NUM_EXTRA_SEGMENTS, "1");
    config.set(NUM_THREADS_PER_WORKER, "1");
    config.set(NUM_WORKERS, "1");
    config.set(SEGMENT_FILTER, "0"); //TODO set "best" version of segment filter as default
    config.set(SEPARATOR_ABUNDANCE, "_");
    config.set(THRESHOLD, "1");

    config.set(SWARM_BOUNDARY, "3");
    config.set(SWARM_DEREPLICATE, "0");
    config.set(SWARM_FASTIDIOUS, "0");
    config.set(SWARM_FASTIDIOUS_CHECKING_MODE, "0");
    config.set(SWARM_NO_OTU_BREAKING, "0");
    config.set(SWARM_NUM_EXPLORERS, "1");
    config.set(SWARM_NUM_GRAFTERS, "1");
    config.set(SWARM_NUM_THREADS_PER_CHECK, "1");

    /* Determine parameter values */

    // determine used config file and load basic configuration
    for (int i = 1; i < argc; i++) {

        if (parameters[argv[i]] == 9 || parameters[argv[i]] == 10) {

            configFile = argv[++i];
            break;

        }

    }

    config.read(configFile, true);

    // read remaining configuration from command line (potentially overwriting values from the config file)
    // converting back and forth between string and integer values is used to "detect" (some) meaningless inputs
    for (int i = 1; i < argc; i++) {

        // one-part parameters
        switch (parameters[argv[i]]) {
            case 109: //fallthrough -sn to --swarm-no-otu-breaking
            case 110:
                config.set(SWARM_NO_OTU_BREAKING, "1");
                break;

            case 111: //fallthrough -sd to --swarm-dereplicate
            case 112:
                config.set(SWARM_DEREPLICATE, "1");
                break;

            case 113: //fallthrough -sf to --swarm-fastidious
            case 114:
                config.set(SWARM_FASTIDIOUS, "1");
                break;

            default:
                // do nothing
                break;
        }

        // two-part parameters
        if (i + 1 != argc) {
            switch (parameters[argv[i]]) {
                case 1: //fallthrough -a to --alphabet
                case 2:
                    config.set(FILTER_ALPHABET, std::to_string(true));
                    config.set(ALPHABET, argv[++i]);
                    break;

                case 3: //fallthrough -e to --extra-segments
                case 4:
                    val = std::stoul(argv[++i]);
                    config.set(NUM_EXTRA_SEGMENTS, std::to_string(val));
                    break;

                case 5: //fallthrough -n to --name
                case 6:
                    name = argv[++i]; // further handling after reading all parameters
                    break;

                case 7: //fallthrough -t to --threshold
                case 8:
                    val = std::stoul(argv[++i]);
                    config.set(THRESHOLD, std::to_string(val));
                    break;

                case 11: //fallthrough -i to --input-files
                case 12:
                    config.set(FILE_LIST, argv[++i]);
                    break;

                case 13: //fallthrough -s to --seg-filter
                case 14:
                    val = std::stoul(argv[++i]);
                    config.set(SEGMENT_FILTER, std::to_string(val));
                    break;

                case 15: //fallthrough -o to --output-file
                case 16:
                    config.set(MATCHES_OUTPUT_FILE, argv[++i]);
                    break;


                case 101: //fallthrough -si to --swarm-internal
                case 102:
                    config.set(SWARM_OUTPUT_INTERNAL, argv[++i]);
                    break;

                case 103: //fallthrough -so to --swarm-output
                case 104:
                    config.set(SWARM_OUTPUT_OTUS, argv[++i]);
                    break;

                case 105: //fallthrough -ss to --swarm-statistics
                case 106:
                    config.set(SWARM_OUTPUT_STATISTICS, argv[++i]);
                    break;

                case 107: //fallthrough -sw to --swarm-seeds
                case 108:
                    config.set(SWARM_OUTPUT_SEEDS, argv[++i]);
                    break;

                case 115: //fallthrough -sb to --swarm-boundary
                case 116:
                    val = std::stoul(argv[++i]);
                    config.set(SWARM_BOUNDARY, std::to_string(val));
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

                case 1003:
                    val = std::stoul(argv[++i]);
                    config.set(NUM_THREADS_PER_WORKER, std::to_string(val));
                    break;

                case 1004:
                    val = std::stoul(argv[++i]);
                    config.set(NUM_WORKERS, std::to_string(val));
                    break;

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

}