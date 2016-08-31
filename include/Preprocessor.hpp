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

#ifndef SCT_PJ_PREPROCESSOR_HPP
#define SCT_PJ_PREPROCESSOR_HPP

#include <regex>

#include "Base.hpp"
#include "Utility.hpp"


namespace SCT_PJ {
namespace Preprocessor {

    struct Defline {// representation of a header / defline of a FASTA file entry

        std::string id;
        std::string information;
        numSeqs_t abundance;

        Defline() {

            id = "";
            information = "";
            abundance = 0;

        }

        Defline(std::string i, std::string info, numSeqs_t a) {

            id = i;
            information = info;
            abundance = a;

        }

    };


    Defline parseDescriptionLine(const std::string defline, const char sep);

    void lowerCase(std::string& s);

    LengthGroups& appendInput(LengthGroups& lg, const std::string fileName, const char sep);

    LengthGroups& filterMinLength(LengthGroups& lg, const numSeqs_t minLen);

    LengthGroups& filterMaxLength(LengthGroups& lg, const numSeqs_t maxLen);

    LengthGroups& filterMinMaxLength(LengthGroups& lg, const numSeqs_t minLen, const numSeqs_t maxLen);

    LengthGroups& filterAlphabet(LengthGroups& lg, const std::string& alph);

    LengthGroups& filterRegex(LengthGroups& lg, const std::regex& expr);

    AmpliconPools* run(const Config<std::string>& conf, const std::vector<std::string>& fileNames);

}
}


#endif //SCT_PJ_PREPROCESSOR_HPP
