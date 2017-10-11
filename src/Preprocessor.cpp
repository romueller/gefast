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

#include <fstream>
#include <limits>

#include "../include/Preprocessor.hpp"


namespace GeFaST {

Preprocessor::Defline Preprocessor::parseDescriptionLine(const std::string& defline, const std::string sep) {

    auto pos = defline.find(' ');
    numSeqs_t abundance = 1;
    std::string firstPart = defline.substr(1, pos - 1);
    std::string secondPart = defline.substr(pos + 1);

    pos = firstPart.find(sep);

    if (pos != std::string::npos) {
        abundance = std::stoul(firstPart.substr(pos + sep.size()));
        firstPart = firstPart.substr(0, pos);
    }

    return Defline(firstPart, secondPart, abundance);

}


namespace Preprocessor {
    char convert[128] = { // upper-case for A-Z, everything else remains the same
            0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  15,
            16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
            32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
            48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
            64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
            80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,
            96,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
            80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90, 123, 124, 125, 126, 127
    };
}

void Preprocessor::upperCase(std::string& s) {

    //std::transform(s.begin(), s.end(), s.begin(), ::toupper); // safer, more general (but slower) version
    for (auto i = 0; i < s.size(); i++) {
        s[i] = convert[s[i]];
    }

}

bool Preprocessor::checkSequence(const std::string& seq, const std::string& alphabet, lenSeqs_t minLen, lenSeqs_t maxLen, bool flagAlph, int flagLen) {

    bool valid = true;

    if (flagAlph) {
        valid = seq.find_first_not_of(alphabet) == std::string::npos;
    }

    switch (flagLen) {
        case 1: { // 01 = only max
            valid = valid && seq.length() <= maxLen;
            break;
        }

        case 2: { // 10 = only min
            valid = valid && minLen <= seq.length();
            break;
        }

        case 3: { // 11 = min & max
            valid = valid && minLen <= seq.length() && seq.length() <= maxLen;
            break;
        }

        default: {
            // do nothing
        }
    }

    return valid;

}

unsigned long long Preprocessor::analyseInput(const Config<std::string>& conf, std::map<lenSeqs_t, numSeqs_t>& counts, const std::string fileName, const std::string sep) {

    std::ifstream iStream(fileName);
    if (!iStream.good()) {

        std::cerr << "ERROR: File '" << fileName << "' not opened correctly. No sequences are read from it." << std::endl;
        return 0;

    }

    // set up filters (depending on configuration)
    lenSeqs_t minLength = 0;
    lenSeqs_t maxLength = std::numeric_limits<lenSeqs_t>::max();
    std::string alphabet;
    unsigned long long totalLength = 0;

    int flagLength = std::stoi(conf.get(FILTER_LENGTH));

#if QGRAM_FILTER

    bool flagAlph = true;
    alphabet = "ACGTU";

#else

    bool flagAlph = bool(std::stoi(conf.get(FILTER_ALPHABET)));
    if (flagAlph) {
        alphabet = conf.get(ALPHABET);
    }

#endif

    switch (flagLength) {
        case 1: { // 01 = only max
            maxLength = std::stoul(conf.get(MAX_LENGTH));
            break;
        }

        case 2: { // 10 = only min
            minLength = std::stoul(conf.get(MIN_LENGTH));
            break;
        }

        case 3: { // 11 = min & max
            minLength = std::stoul(conf.get(MIN_LENGTH));
            maxLength = std::stoul(conf.get(MAX_LENGTH));
            break;
        }

        default: {
            // do nothing
        }
    }


    Defline dl;
    std::string line, seq;
    bool first = true;

    while (std::getline(iStream, line).good()) {

        if (line.empty() || line[0] == ';') continue; // skip empty and comment lines (begin with ';')

        if (line[0] == '>') { // header line

            if (first) { // first entry found (no previous entry to finish), simply parse header

                dl = parseDescriptionLine(line, sep);
                first = false;

            } else { // finish and store previous entry, then collect parse header of new entry

                upperCase(seq);

                if (checkSequence(seq, alphabet, minLength, maxLength, flagAlph, flagLength)) {

                    counts[seq.length()]++;
                    totalLength += dl.id.length() + seq.length() + 2;

                }

                seq.clear();
                dl = parseDescriptionLine(line, sep);

            }

        } else { // still the same entry, continue to collect sequence
            seq += line;
        }
    }

    if (!first) { // ensures that last entry (if any) is written to file.

        upperCase(seq);

        if (checkSequence(seq, alphabet, minLength, maxLength, flagAlph, flagLength)) {

            counts[seq.length()]++;
            totalLength += dl.id.length() + seq.length() + 2;

        }

    }

    return totalLength;

}

// adds contents of specified file (if readable) to the given AmpliconPools object
// amplicons are assigned to the respective pools established by prior calls of analyseInput()
// "normalises" to lower-case letters, filters according to configuration
void Preprocessor::appendInput(const Config<std::string>& conf, AmpliconPools& pools, std::map<lenSeqs_t, numSeqs_t>& poolMap, const std::string fileName, const std::string sep) {

    std::ifstream iStream(fileName);
    if (!iStream.good()) {

        std::cerr << "ERROR: File '" << fileName << "' not opened correctly. No sequences are read from it." << std::endl;
        return;

    }

    // set up filters (depending on configuration)
    lenSeqs_t minLength = 0;
    lenSeqs_t maxLength = std::numeric_limits<lenSeqs_t>::max();
    std::string alphabet;

    int flagLength = std::stoi(conf.get(FILTER_LENGTH));

#if QGRAM_FILTER

    bool flagAlph = true;
    alphabet = "ACGTU";

#else

    bool flagAlph = bool(std::stoi(conf.get(FILTER_ALPHABET)));
    if (flagAlph) {
        alphabet = conf.get(ALPHABET);
    }

#endif

    switch (flagLength) {
        case 1: { // 01 = only max
            maxLength = std::stoul(conf.get(MAX_LENGTH));
            break;
        }

        case 2: { // 10 = only min
            minLength = std::stoul(conf.get(MIN_LENGTH));
            break;
        }

        case 3: { // 11 = min & max
            minLength = std::stoul(conf.get(MIN_LENGTH));
            maxLength = std::stoul(conf.get(MAX_LENGTH));
            break;
        }

        default: {
            // do nothing
        }
    }


    Defline dl;
    std::string line, seq;
    bool first = true;


    while (std::getline(iStream, line).good()) {

        if (line.empty() || line[0] == ';') continue; // skip empty and comment lines (begin with ';')

        if (line[0] == '>') { // header line

            if (first) { // first entry found (no previous entry to finish), simply parse header

                dl = parseDescriptionLine(line, sep);
                first = false;

            } else { // finish and store previous entry, then collect parse header of new entry

                upperCase(seq);

                if (checkSequence(seq, alphabet, minLength, maxLength, flagAlph, flagLength)) {
                    pools.add(poolMap[seq.length()], dl.id, seq, dl.abundance);
                }

                seq.clear();
                dl = parseDescriptionLine(line, sep);

            }

        } else { // still the same entry, continue to collect sequence
            seq += line;
        }
    }

    if (!first) { // ensures that last entry (if any) is written to file.

        upperCase(seq);

        if (checkSequence(seq, alphabet, minLength, maxLength, flagAlph, flagLength)) {
            pools.add(poolMap[seq.length()], dl.id, seq, dl.abundance);
        }

    }

}


AmpliconPools* Preprocessor::run(const Config<std::string>& conf, const std::vector<std::string>& fileNames) {

    std::string sep = conf.get(SEPARATOR_ABUNDANCE);

    unsigned long long totalLength = 0;
    std::map<lenSeqs_t, numSeqs_t> counts;
    for (auto iter = fileNames.begin(); iter != fileNames.end(); iter++) {

        std::cout << "Analysing inputs: file " << (iter - fileNames.begin()) + 1 << " / " << fileNames.size() << "\r" << std::flush;
        totalLength += analyseInput(conf, counts, *iter, sep);

    }
    std::cout << std::endl;

    AmpliconPools* pools = new AmpliconPools(counts, totalLength, std::stoul(conf.get(THRESHOLD)));

    for (auto iter = fileNames.begin(); iter != fileNames.end(); iter++) {

        std::cout << "Reading inputs: file " << (iter - fileNames.begin()) + 1 << " / " << fileNames.size() << "\r" << std::flush;
        appendInput(conf, *pools, counts, *iter, sep);

    }
    std::cout << std::endl;

    for (numSeqs_t p = 0; p < pools->numPools(); p++) {

        std::cout << "Sorting amplicons: pool " << p + 1 << " / " << pools->numPools() << "\r" << std::flush;

        auto ac = pools->get(p);
        std::sort(ac->begin(), ac->end(),
                  [](const Amplicon& amplA, const Amplicon& amplB) {
                      return (amplA.abundance > amplB.abundance) || ((amplA.abundance == amplB.abundance) && (strcmp(amplA.id, amplB.id) < 0));
                  }
        );

    }
    std::cout << std::endl;

    return pools;

}

}