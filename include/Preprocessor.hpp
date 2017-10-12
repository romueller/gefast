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


    // splits description line into actual header, abundance value and additional information (if any)
    Defline parseDescriptionLine(const std::string& defline, const std::string sep);

    // sets sequence in upper case
    void upperCase(std::string& s);

    // checks whether the given sequence satisfies the given filter criteria
    bool checkSequence(const std::string& seq, const std::string& alphabet, lenSeqs_t minLen, lenSeqs_t maxLen,
                       bool flagAlph, int flagLen);

    /*
     * Scans the given input file for sequences that satisfy the preprocessor's filters.
     * The sequences passing the filters are counted and, according to their length, recorded in the counts mapping.
     * For passing sequences, also the combined length of identifier and sequence (+2 for terminating \0's) is summed up
     * and finally returned.
     */
    unsigned long long analyseInput(const Config<std::string>& conf, std::map<lenSeqs_t, numSeqs_t>& counts,
                                    const std::string fileName, const std::string sep);

    /*
     * Rereads the analysed input and inserts the suitable amplicons into the pools.
     * poolMap assigns each sequence length the index of the pool into which it is to be inserted.
     */
    void appendInput(const Config<std::string>& conf, AmpliconPools& pools, std::map<lenSeqs_t, numSeqs_t>& poolMap,
                     const std::string fileName, const std::string sep);

    /*
     * Manages the overall preprocessing step.
     *
     * First, all input files are analysed to determine
     *  (a) how many sequences of which lengths will be retained after the preprocessing, and
     *  (b) how much characters are necessary to store the sequences and their identifiers as C-strings in one long char array.
     * Second, the amplicon pools are initialised based on the results of above analysis (see the constructor of
     * AmpliconPools for important details).
     * Third, all input files are read a second time to get the actual amplicons and fill the amplicon pools.
     * Finally, sort the amplicons within each pool by abundance (using the lexicographical order of the headers as the tie-breaker).
     */
    AmpliconPools* run(const Config<std::string>& conf, const std::vector<std::string>& fileNames);

}
}


#endif //GEFAST_PREPROCESSOR_HPP
