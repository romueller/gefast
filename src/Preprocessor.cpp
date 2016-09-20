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

#include <fstream>

#include "../include/Preprocessor.hpp"


namespace SCT_PJ {

Preprocessor::Defline Preprocessor::parseDescriptionLine(const std::string defline, const char sep) {

    auto pos = defline.find(' ');
    numSeqs_t abundance = 1;
    std::string firstPart = defline.substr(1, pos - 1);
    std::string secondPart = defline.substr(pos + 1);

    pos = firstPart.find(sep);

    if (pos != std::string::npos) {
        abundance = std::stoul(firstPart.substr(pos + 1));
        firstPart = firstPart.substr(0, pos);
    }

    return Defline(firstPart, secondPart, abundance);

}


namespace Preprocessor {
    char convert[128] = { // lower-case for A-Z, everything else remains the same
            0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  15,
            16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
            32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
            48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
            64,  97,  98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
            112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 91,  92,  93,  94,  95,
            96,  97,  98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
            112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127
    };
}

void Preprocessor::lowerCase(std::string& s) {

    //std::transform(s.begin(), s.end(), s.begin(), ::tolower); // safer, more general (but slower) version
    for (auto i = 0; i < s.size(); i++) {
        s[i] = convert[s[i]];
    }

}

//adds contents of specified file (if readable) to the given LengthGroups object, returns pointer to this object
// "normalises" to lower-case letters
LengthGroups& Preprocessor::appendInput(LengthGroups& ampls, const std::string fileName, const char sep) {

    std::ifstream iStream(fileName);
    if (!iStream.good()) {

        std::cerr << "ERROR: File '" << fileName << "' not opened correctly. No sequences are read from it."
                  << std::endl;
        return ampls;

    }

#if INPUT_RANK
    numSeqs_t prevSize = ampls.size();
#endif
    numSeqs_t cnt = 0;
    Defline dl;
    std::string line, seq;


    while (std::getline(iStream, line).good()) {

        if (line.empty() || line[0] == ';') continue; //skip empty and comment lines (begin with ';')

        if (line[0] == '>') { //header line

            cnt++;

            if (cnt == 1) { //first entry found (no previous entry to finish), simply parse header

                dl = parseDescriptionLine(line, sep);

            } else { //finish and store previous entry, then collect parse header of new entry

                lowerCase(seq);
#if INPUT_RANK
                ampls.add(Amplicon(dl.id, seq, dl.abundance, prevSize + cnt));
#else
                ampls.add(Amplicon(dl.id, seq, dl.abundance));
#endif

                seq.clear();
                dl = parseDescriptionLine(line, sep);

            }

        } else { //still the same entry, continue to collect sequence

            seq += line;

        }
    }

    if (cnt > 0) { //ensures that last entry (if any) is written to file.

        lowerCase(seq);
#if INPUT_RANK
        ampls.add(Amplicon(dl.id, seq, dl.abundance, prevSize + cnt));
#else
        ampls.add(Amplicon(dl.id, seq, dl.abundance));
#endif

    }

    return ampls;


}


LengthGroups& Preprocessor::filterMinLength(LengthGroups& lg, const numSeqs_t minLen) {

    auto groups = lg.getGroups();
    auto endIter = groups.begin();
    for (; endIter != groups.end() && endIter->first < minLen; endIter++) {}

    groups.erase(groups.begin(), endIter);

    return lg;

}


LengthGroups& Preprocessor::filterMaxLength(LengthGroups& lg, const numSeqs_t maxLen) {

    auto groups = lg.getGroups();
    auto beginIter = groups.begin();
    for (; beginIter != groups.end() && beginIter->first <= maxLen; beginIter++) {}

    groups.erase(beginIter, groups.end());

    return lg;

}


LengthGroups& Preprocessor::filterMinMaxLength(LengthGroups& lg, const numSeqs_t minLen, const numSeqs_t maxLen) {

    return filterMaxLength(filterMinLength(lg, minLen), maxLen);

}


// bei Benutzung von regex_match() Problem mit erster Gruppe, die Sequenzen länger als 3000 Zeichen enthält -- warum?
// regex_search() mit negiertem Alphabet bedeutend langsamer
LengthGroups& Preprocessor::filterAlphabet(LengthGroups& lg, const std::string& alph) {

    auto groups = lg.getGroups();

    // remove (in each length group) all amplicons with unsuitable sequences
    for (auto groupIter = groups.begin(); groupIter != groups.end(); groupIter++) {

        groupIter->second.erase(
                std::remove_if(
                        groupIter->second.begin(),
                        groupIter->second.end(),
                        [alph](Amplicon& a) {
                            return a.seq.find_first_not_of(alph) != std::string::npos;
                        }),
                groupIter->second.end()
        );

    }

    // clear LengthGroups from empty groups
    for (auto iter = groups.begin(); iter != groups.end();) {
        if (iter->second.size() == 0) {
            iter = groups.erase(iter);
        } else {
            iter++;
        }
    }

    return lg;

}


// bei Benutzung von regex_match() Problem mit erster Gruppe, die Sequenzen länger als 3000 Zeichen enthält -- warum?
// e.g. filterRegex(*ampls, std::regex("[^aAcCgGtTuU]"));
LengthGroups& Preprocessor::filterRegex(LengthGroups& lg, const std::regex& expr) {

    auto groups = lg.getGroups();

    // remove (in each length group) all amplicons with unsuitable sequences
    for (auto groupIter = groups.begin(); groupIter != groups.end(); groupIter++) {

        groupIter->second.erase(
                std::remove_if(
                        groupIter->second.begin(),
                        groupIter->second.end(),
                        [expr](Amplicon& a) {
                            return !std::regex_search(a.seq, expr);
                        }),
                groupIter->second.end()
        );

    }

    // clear LengthGroups from empty groups
    for (auto iter = groups.begin(); iter != groups.end();) {
        if (iter->second.size() == 0) {
            iter = groups.erase(iter);
        } else {
            iter++;
        }
    }

    return lg;

}


AmpliconPools* Preprocessor::run(const Config<std::string>& conf, const std::vector<std::string>& fileNames) {

    // collect amplicons from all input files
    LengthGroups* ampls = new LengthGroups();
    char sep = conf.get(SEPARATOR_ABUNDANCE)[0];
    for (auto iter = fileNames.begin(); iter != fileNames.end(); iter++) {
        appendInput(*ampls, *iter, sep);
    }


    // apply filters (depending on configuration)
    lenSeqs_t minLength, maxLength;
    std::string alphabet;

    int flagLength = std::stoi(conf.get(FILTER_LENGTH));
    bool flagAlph = bool(std::stoi(conf.get(FILTER_ALPHABET)));

    switch (flagLength) {
        case 1: { //01 = only max
            maxLength = std::stoul(conf.get(MAX_LENGTH));
            Preprocessor::filterMaxLength(*ampls, maxLength);
            break;
        }

        case 2: { //10 = only min
            minLength = std::stoul(conf.get(MIN_LENGTH));
            Preprocessor::filterMinLength(*ampls, minLength);
            break;
        }

        case 3: { //11 = min & max
            minLength = std::stoul(conf.get(MIN_LENGTH));
            maxLength = std::stoul(conf.get(MAX_LENGTH));
            Preprocessor::filterMinMaxLength(*ampls, minLength, maxLength);
            break;
        }

        default: {
            //do nothing
        }
    }

    if (flagAlph) {

        alphabet = conf.get(ALPHABET);
        Preprocessor::filterAlphabet(*ampls, alphabet);

    }


    // pool length groups into amplicon collections
    lenSeqs_t threshold = std::stoul(conf.get(THRESHOLD));
    AmpliconPools* pools = ampls->pool(threshold);

    delete ampls;

    return pools;

}

}