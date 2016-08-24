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

#ifndef SCT_PJ_UTILITY_HPP
#define SCT_PJ_UTILITY_HPP

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>


namespace SCT_PJ {

const std::string DEFAULT_CONFIG_FILE = "default.conf";
const std::string DEFAULT_JOB_NAME = "SCT-PJ";


// list of configuration parameters (changes should be mirrored in method initParamNames() of Config)
enum ConfigParameters {
    ALPHABET,               // allowed alphabet for the amplicon sequences
    FILE_LIST,              // file containing list of input file names
    FILTER_ALPHABET,        // flag for the alphabet filter
    FILTER_LENGTH,          // flag for the length filter
    FILTER_REGEX,           // flag for the regex filter
    INFO_FOLDER,            // name of the folder for storing files showing the configuration information of the current job
    MATCHES_OUTPUT_FILE,    // name of the output file containing all matches found
    MAX_LENGTH,             // maximal sequence length
    MIN_LENGTH,             // minimal sequence length
    NAME,                   // name of the job to be executed
    NUM_EXTRA_SEGMENTS,     // parameter of the pigeonhole principle (segment filter)
    NUM_THREADS_PER_WORKER, // number of parallel threads employed by each work
    NUM_WORKERS,            // number of parallel workers
    SEGMENT_FILTER,         // mode of the segment filter (forward, backward, forward-backward, backward-forward)
    SWARM_OUTPUT_INTERNAL,  // name of the output file corresponding to swarm's output option -i (internal structures)
    SWARM_OUTPUT_OTUS,      // name of the output file corresponding to swarm's output option -o (OTUs)
    SWARM_OUTPUT_STATISTICS,// name of the output file corresponding to swarm's output option -s (statistics file)
    SWARM_OUTPUT_SEEDS,     // name of the output file corresponding to swarm's output option -w (seeds)
    THRESHOLD,              // edit distance threshold for the clustering
    VERSION                 // version number of the program
};


// from http://stackoverflow.com/questions/18837857/cant-use-enum-class-as-unordered-map-key
struct EnumClassHash {
    template<typename T>
    std::size_t operator()(T t) const {
        return static_cast<std::size_t>(t);
    }
};

/**
 * Manages the configuration parameters.
 * Every parameter consists of a key and an associated value.
 * All values are stored as the same type V (e.g. string) and have to be converted
 * into the actual type at the points of use (if necessary).
 */
template<typename V>
class Config {

public:
    Config() {

        initParamNames();

    }

    /**
     * Assumed syntax:
     *  - Line comments are allowed and start with #.
     *  - Every comment is written in its own line.
     *  - Empty lines are allowed.
     *  - Every configuration parameter is written in its own line.
     *  - A line containing a configuration parameter must have
     *    the following form: <key>=<value>
     */
    Config(const std::string file) {

        initParamNames();

        std::ifstream iStream(file);
        std::string line;
        unsigned long delimPos;

        while (std::getline(iStream, line).good()) {

            if (line.empty() || line.front() == '#') continue;

            delimPos = line.find('=');
            conf_[paramNames_[line.substr(0, delimPos)]] = line.substr(delimPos + 1);

        }

        iStream.close();

    }

    void set(const ConfigParameters key, const V value) {
        conf_[key] = value;
    }

    V get(const ConfigParameters key) const { // throws exception if key is not contained
        return conf_.at(key);
    }

    // return true if something is stored for the specified key
    bool peek(const ConfigParameters key) const {

        return conf_.find(key) != conf_.end();

    }

    void print(std::ostream& stream) const {

        for (auto iter = paramNames_.begin(); iter != paramNames_.end(); iter++) {

            if (peek(iter->second)) {
                stream << iter->first << " = " << conf_.at(iter->second) << std::endl;
            }

        }

    }


private:
    std::unordered_map<ConfigParameters, V, EnumClassHash> conf_;

    std::unordered_map<std::string, ConfigParameters> paramNames_;

    void initParamNames() {

        paramNames_ = std::unordered_map<std::string, ConfigParameters>(
                {
                        {"ALPHABET",               ALPHABET},
                        {"FILE_LIST",              FILE_LIST},
                        {"FILTER_ALPHABET",        FILTER_ALPHABET},
                        {"FILTER_LENGTH",          FILTER_LENGTH},
                        {"FILTER_REGEX",           FILTER_REGEX},
                        {"INFO_FOLDER",            INFO_FOLDER},
                        {"MAX_LENGTH",             MAX_LENGTH},
                        {"MIN_LENGTH",             MIN_LENGTH},
                        {"NAME",                   NAME},
                        {"NUM_EXTRA_SEGMENTS",     NUM_EXTRA_SEGMENTS},
                        {"NUM_THREADS_PER_WORKER", NUM_THREADS_PER_WORKER},
                        {"NUM_WORKERS",            NUM_WORKERS},
                        {"MATCHES_OUTPUT_FILE",    MATCHES_OUTPUT_FILE},
                        {"SEGMENT_FILTER",         SEGMENT_FILTER},
                        {"SWARM_OUTPUT_INTERNAL",  SWARM_OUTPUT_INTERNAL},
                        {"SWARM_OUTPUT_OTUS",      SWARM_OUTPUT_OTUS},
                        {"SWARM_OUTPUT_STATISTICS",SWARM_OUTPUT_STATISTICS},
                        {"SWARM_OUTPUT_SEEDS",     SWARM_OUTPUT_SEEDS},
                        {"THRESHOLD",              THRESHOLD},
                        {"VERSION",                VERSION}
                }
        );

    }

};


// read file names from a specified file
// each line consists of a single file name (including its path)
// empty lines and comment lines (starting with ';') are allowed
std::vector<std::string> readFileList(const std::string listFile);


// process program arguments to get the configuration (parameters)
// command line overwrites information read from config file
Config<std::string> getConfiguration(int argc, const char* argv[]);


// write job parameters to file
void writeJobParameters(const std::string oFile, const Config<std::string>& conf, const std::vector<std::string>& inputFiles);

}

#endif //SCT_PJ_UTILITY_HPP