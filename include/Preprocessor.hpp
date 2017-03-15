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

#ifndef GEFAST_PREPROCESSOR_HPP
#define GEFAST_PREPROCESSOR_HPP

#include <regex>

#include "Base.hpp"
#include "Utility.hpp"


namespace GeFaST {
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


    Defline parseDescriptionLine(const std::string defline, const std::string sep);

    void lowerCase(std::string& s);

    LengthGroups& appendInput(LengthGroups& lg, const std::string fileName, const std::string sep);

    LengthGroups& filterMinLength(LengthGroups& lg, const numSeqs_t minLen);

    LengthGroups& filterMaxLength(LengthGroups& lg, const numSeqs_t maxLen);

    LengthGroups& filterMinMaxLength(LengthGroups& lg, const numSeqs_t minLen, const numSeqs_t maxLen);

    LengthGroups& filterAlphabet(LengthGroups& lg, const std::string& alph);

    LengthGroups& filterRegex(LengthGroups& lg, const std::regex& expr);

    AmpliconPools* run(const Config<std::string>& conf, const std::vector<std::string>& fileNames);

}
}


#endif //GEFAST_PREPROCESSOR_HPP
