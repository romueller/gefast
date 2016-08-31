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

    parameters["--min-length"] = 1001; // optional
    parameters["--max-length"] = 1002; // optional
    parameters["--per-worker"] = 1003; // optional (has default)
    parameters["--workers"] = 1004; // optional (has default)
    parameters["--info-file"] = 1005; // optional
    parameters["--sep-abundance"] = 1006; // optional (has default)


    std::string
            alphabet = "", // 1,2
            name = DEFAULT_JOB_NAME, // 5,6
            configFile = "", // 9,10
            fileList = "", // 11,12
            outputFile = "", // 15,16
            oFileSwarmInternal = "", // 101,102
            oFileSwarmOutput = "", // 103,104
            oFileSwarmStatistics = "", // 105,106
            oFileSwarmSeeds = "", // 107,108
            infoFile = "", // 1005
            sepAbundance = ";" // 1006 //TODO choose "best" default (underscore as in swarm?)
    ;

    unsigned long long
            extraSegs = 1, // 3,4
            threshold = 1, // 7,8
            segFilter = 0, // 13,14 //TODO set "best" version of segment filter as default
            minLength = 0, // 1001
            maxLength = 0, // 1002
            numThreadsPerWorker = 1, // 1003
            numWorkers = 1 // 1004
    ;

    bool
            filterAlphabet = false, // 1,2
            filterMinLength = false, // 1001
            filterMaxLength = false, // 1002
            flagExtraSegs = false, // (3,4)
            flagThreshold = false, // (7,8)
            flagSegFilter = false, // (13,14)
            flagNumThreadsPerWorker = false, // (1003)
            flagNumWorkers = false, // (1004)
            flagDereplicate = false, // 111,112
            flagSepAbundance = false // 1006
    ;

    int flagNoOtuBreaking = -1; // (109,110)


    /* Determine parameter values */

    for (int i = 1; i < argc; i++) {

        // one-part parameters
        switch (parameters[argv[i]]) {
            case 111: //fallthrough -sd to -swarm-dereplicate
            case 112:
                flagDereplicate = true;
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
                    alphabet = argv[++i];
                    filterAlphabet = true;
                    break;

                case 3: //fallthrough -e to --extra-segments
                case 4:
                    extraSegs = std::stoul(argv[++i]);
                    flagExtraSegs = true;
                    break;

                case 5: //fallthrough -n to --name
                case 6:
                    name = argv[++i];
                    break;

                case 7: //fallthrough -t to --threshold
                case 8:
                    threshold = std::stoul(argv[++i]);
                    flagThreshold = true;
                    break;

                case 9: //fallthrough -c to --config
                case 10:
                    configFile = argv[++i];
                    break;

                case 11: //fallthrough -i to --input-files
                case 12:
                    fileList = argv[++i];
                    break;

                case 13: //fallthrough -s to --seg-filter
                case 14:
                    segFilter = std::stoul(argv[++i]);
                    flagSegFilter = true;
                    break;

                case 15: //fallthrough -o to --output-file
                case 16:
                    outputFile = argv[++i];
                    break;


                case 101: //fallthrough -si to --swarm-internal
                case 102:
                    oFileSwarmInternal = argv[++i];
                    break;

                case 103: //fallthrough -so to --swarm-output
                case 104:
                    oFileSwarmOutput = argv[++i];
                    break;

                case 105: //fallthrough -ss to --swarm-statistics
                case 106:
                    oFileSwarmStatistics = argv[++i];
                    break;


                case 107: //fallthrough -sw to --swarm-seeds
                case 108:
                    oFileSwarmSeeds = argv[++i];
                    break;

                case 109: //fallthrough -sn to --swarm-no-otu-breaking
                case 110:
                    flagNoOtuBreaking = std::stoi(argv[++i]);
                    break;

                // parameters without abbreviation
                case 1001:
                    minLength = std::stoul(argv[++i]);
                    filterMinLength = true;
                    break;

                case 1002:
                    maxLength = std::stoul(argv[++i]);
                    filterMaxLength = true;
                    break;

                case 1003:
                    numThreadsPerWorker = std::stoul(argv[++i]);
                    flagNumThreadsPerWorker = true;
                    break;

                case 1004:
                    numWorkers = std::stoul(argv[++i]);
                    flagNumWorkers = true;
                    break;

                case 1005:
                    infoFile = std::stoul(argv[++i]);
                    break;

                case 1006:
                    sepAbundance = argv[++i];
                    flagSepAbundance = true;
                    break;

                default:
                    std::cout << "Unknown parameter: " << argv[i] << " (is ignored)" << std::endl;
                    break;
            }
        }

    }


    /* Populate configuration */
    // load (basic) parameters from config file
    Config<std::string> config = Config<std::string>((configFile == "") ? DEFAULT_CONFIG_FILE : configFile);

    // overwrite with / add command-line inputs
    config.set(FILTER_ALPHABET, std::to_string(filterAlphabet || config.peek(ALPHABET)));
    if (alphabet != "") config.set(ALPHABET, alphabet);

    if (fileList != "") config.set(FILE_LIST, fileList);

    config.set(FILTER_LENGTH, std::to_string(
            0
            + (filterMaxLength || config.peek(MAX_LENGTH))
            + 2 * (filterMinLength || config.peek(MIN_LENGTH))
    ));
    if (maxLength != 0) config.set(MAX_LENGTH, std::to_string(maxLength));
    if (minLength != 0) config.set(MIN_LENGTH, std::to_string(minLength));

    if (name == DEFAULT_JOB_NAME) { // append current time / date to default name

        time_t rawTime;
        time(&rawTime);
        std::string curTime = ctime(&rawTime);
        std::replace(curTime.begin(), curTime.end(), ' ', '_');
        name = name + "_" + curTime.substr(0, curTime.length() - 1);

    }
    config.set(NAME, name);

    if (flagExtraSegs || !config.peek(NUM_EXTRA_SEGMENTS)) config.set(NUM_EXTRA_SEGMENTS, std::to_string(extraSegs));

    if (flagNumThreadsPerWorker || !config.peek(NUM_THREADS_PER_WORKER)) config.set(NUM_THREADS_PER_WORKER, std::to_string(numThreadsPerWorker));

    if (flagNumWorkers || !config.peek(NUM_WORKERS)) config.set(NUM_WORKERS, std::to_string(numWorkers));

    if (flagThreshold || !config.peek(THRESHOLD)) config.set(THRESHOLD, std::to_string(threshold));

    if (flagSegFilter || !config.peek(SEGMENT_FILTER)) config.set(SEGMENT_FILTER, std::to_string(segFilter));

    if (outputFile != "") config.set(MATCHES_OUTPUT_FILE, outputFile);


    if (oFileSwarmInternal != "") config.set(SWARM_OUTPUT_INTERNAL, oFileSwarmInternal);

    if (oFileSwarmOutput != "") config.set(SWARM_OUTPUT_OTUS, oFileSwarmOutput);

    if (oFileSwarmStatistics != "") config.set(SWARM_OUTPUT_STATISTICS, oFileSwarmStatistics);

    if (oFileSwarmSeeds != "") config.set(SWARM_OUTPUT_SEEDS, oFileSwarmSeeds);

    if (flagNoOtuBreaking != -1) config.set(SWARM_NO_OTU_BREAKING, std::to_string(flagNoOtuBreaking));

    if (flagDereplicate) config.set(SWARM_DEREPLICATE, "1");

    if (flagSepAbundance || !config.peek(SEPARATOR_ABUNDANCE)) config.set(SEPARATOR_ABUNDANCE, sepAbundance);

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