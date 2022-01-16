/*
 * GeFaST
 *
 * Copyright (C) 2016 - 2021 Robert Mueller
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

#include "../include/Preprocessors.hpp"

namespace GeFaST {

    /* === FastaPreprocessor === */

    void FastaPreprocessor::analyse_input_file(const std::string& file, DataStatistics<>& ds, const Configuration& config) {

        std::string sep = config.separator;

        std::ifstream in_stream(file);
        if (!in_stream.good()) {

            std::cerr << " -- WARNING: File '" << file << "' not opened correctly." << std::endl;
            return;

        }

        // set up filters (depending on configuration)
        lenSeqs_t min_length = 0;
        lenSeqs_t max_length = std::numeric_limits<lenSeqs_t>::max();
        switch (config.filter_length) {
            case 1: { // 01 = only max
                max_length = config.max_length;
                break;
            }

            case 2: { // 10 = only min
                min_length = config.min_length;
                break;
            }

            case 3: { // 11 = min & max
                min_length = config.min_length;
                max_length = config.max_length;
                break;
            }

            default: {
                // do nothing
            }
        }

        numSeqs_t min_abundance = 0;
        numSeqs_t max_abundance = std::numeric_limits<numSeqs_t>::max();
        switch (config.filter_abundance) {
            case 1: { // 01 = only max
                max_abundance = config.max_abundance;
                break;
            }

            case 2: { // 10 = only min
                min_abundance = config.min_abundance;
                break;
            }

            case 3: { // 11 = min & max
                min_abundance = config.min_abundance;
                max_abundance = config.max_abundance;
                break;
            }

            default: {
                // do nothing
            }
        }


        Defline dl;
        std::string line, seq;
        bool first = true;

        while (std::getline(in_stream, line).good()) {

            if (line.empty() || line[0] == ';') continue; // skip empty and comment lines (begin with ';')

            if (line[0] == '>') { // header line

                if (first) { // first entry found (no previous entry to finish), simply parse header

                    dl = Defline::parse_description_line(line, sep);
                    first = false;

                } else { // finish and store previous entry, then collect parse header of new entry

                    upper_case(seq);

                    if (check_sequence(seq, config.alphabet, dl.abundance, min_length, max_length, min_abundance, max_abundance)) {
                        ds.record(dl.id, seq, dl.abundance);
                    }

                    seq.clear();
                    dl = Defline::parse_description_line(line, sep);

                }

            } else { // still the same entry, continue to collect sequence
                seq += line;
            }
        }

        if (!first) { // ensures that last entry (if any) is written to file.

            upper_case(seq);

            if (check_sequence(seq, config.alphabet, dl.abundance, min_length, max_length, min_abundance, max_abundance)) {
                ds.record(dl.id, seq, dl.abundance);
            }

        }

    }

    void FastaPreprocessor::read_input_file(const std::string& file, AmpliconStorage& amplicon_storage, const Configuration& config) {

        std::ifstream in_stream(file);
        if (!in_stream.good()) {

            std::cerr << " -- WARNING: File '" << file << "' not opened correctly. No sequences are read from it." << std::endl;
            return;

        }

        // set up filters (depending on configuration)
        lenSeqs_t min_length = 0;
        lenSeqs_t max_length = std::numeric_limits<lenSeqs_t>::max();
        switch (config.filter_length) {
            case 1: { // 01 = only max
                max_length = config.max_length;
                break;
            }

            case 2: { // 10 = only min
                min_length = config.min_length;
                break;
            }

            case 3: { // 11 = min & max
                min_length = config.min_length;
                max_length = config.max_length;
                break;
            }

            default: {
                // do nothing
            }
        }

        numSeqs_t min_abundance = 0;
        numSeqs_t max_abundance = std::numeric_limits<numSeqs_t>::max();
        switch (config.filter_abundance) {
            case 1: { // 01 = only max
                max_abundance = config.max_abundance;
                break;
            }

            case 2: { // 10 = only min
                min_abundance = config.min_abundance;
                break;
            }

            case 3: { // 11 = min & max
                min_abundance = config.min_abundance;
                max_abundance = config.max_abundance;
                break;
            }

            default: {
                // do nothing
            }
        }


        Defline dl;
        ExtraInfo ei; // stays empty for FASTA sequences
        std::string line, seq;
        bool first = true;

        while (std::getline(in_stream, line).good()) {

            if (line.empty() || line[0] == ';') continue; // skip empty and comment lines (begin with ';')

            if (line[0] == '>') { // header line

                if (first) { // first entry found (no previous entry to finish), simply parse header

                    dl = Defline::parse_description_line(line, config.separator);
                    first = false;

                } else { // finish and store previous entry, then collect parse header of new entry

                    upper_case(seq);

                    if (check_sequence(seq, config.alphabet, dl.abundance, min_length, max_length, min_abundance, max_abundance)) {
                        amplicon_storage.add_amplicon(dl, seq, ei);
                    }

                    seq.clear();

                    dl = Defline::parse_description_line(line, config.separator);

                }

            } else { // still the same entry, continue to collect sequence
                seq += line;
            }
        }

        if (!first) { // ensures that last entry (if any) is written to file.

            upper_case(seq);

            if (check_sequence(seq, config.alphabet, dl.abundance, min_length, max_length, min_abundance, max_abundance)) {
                amplicon_storage.add_amplicon(dl, seq, ei);
            }

        }

    }

    AmpliconStorage* FastaPreprocessor::preprocess_inputs(const std::vector<std::string>& files, const Configuration& config) {

        AmpliconStorage* amplicon_storage = nullptr;

        {
            DataStatistics<> ds;

            std::cout << " -- Scanning input files..." << std::endl;
            for (auto& f : files) {
                analyse_input_file(f, ds, config);
            }

            amplicon_storage = config.build_amplicon_storage(ds);

            std::cout << " -- Scan completed." << std::endl;

        }

        std::cout << " -- Reading input files..." << std::endl;
        for (auto& f : files) {
            read_input_file(f, *amplicon_storage, config);
        }

        amplicon_storage->finalise();

        std::cout << " -- Construction of amplicon storage completed." << std::endl;

        return amplicon_storage;

    }

    FastaPreprocessor* FastaPreprocessor::clone() const {
        return new FastaPreprocessor(*this);
    }


    /* === MaxReadRepresentation === */

    void MaxReadRepresentation::initialise(const Defline& dl, const std::string& seq, const std::string& quality_scores, const QualityEncoding<>& qe) {

        id = dl.id;
        abundance = dl.abundance;
        max_scores = quality_scores;

    }

    void MaxReadRepresentation::add_duplicate(const numSeqs_t ab, const std::string& quality_scores, const QualityEncoding<>& qe) {

        abundance += ab;
        for (auto i = 0; i < quality_scores.length(); i++) {
            max_scores[i] = std::max(max_scores[i], quality_scores[i]);
        }

    }

    std::string MaxReadRepresentation::get_scores(const QualityEncoding<>& qe) const {
        return max_scores;
    }

    /* === MeanReadRepresentation === */

    void MeanReadRepresentation::initialise(const Defline& dl, const std::string& seq, const std::string& quality_scores, const QualityEncoding<>& qe) {

        id = dl.id;
        abundance = dl.abundance;
        abundances.push_back(dl.abundance);
        scores.emplace_back(quality_scores);

    }

    void MeanReadRepresentation::add_duplicate(const numSeqs_t ab, const std::string& quality_scores, const QualityEncoding<>& qe) {

        abundance += ab;
        abundances.push_back(ab);
        scores.emplace_back(quality_scores);

    }

    std::string MeanReadRepresentation::get_scores(const QualityEncoding<>& qe) const {

        std::vector<float> probability_sums = std::vector<float>(scores[0].length(), 0);
        for (auto d = 0; d < abundances.size(); d++) {
            for (auto i = 0; i < scores[d].length(); i++) {
                probability_sums[i] += abundances[d] * qe.encoded_quality_to_probability(scores[d][i]);
            }
        }

        char max_score = qe.get_accepted_scores().back();
        std::string mean_scores(probability_sums.size(), ' ');
        for (auto i = 0; i < probability_sums.size(); i++) {
            // limit to the score to the encoding range (as a precaution against imprecision and rounding errors)
            mean_scores[i] = std::min(max_score, qe.probability_to_encoded_quality(probability_sums[i] / abundance));
        }

        return mean_scores;

    }

    /* === ModeReadRepresentation === */

    void ModeReadRepresentation::initialise(const Defline& dl, const std::string& seq, const std::string& quality_scores, const QualityEncoding<>& qe) {

        id = dl.id;
        abundance = dl.abundance;
        score_tallies = std::vector<std::map<char, numSeqs_t>>(quality_scores.length());
        for (auto i = 0; i < quality_scores.length(); i++) {
            score_tallies[i][quality_scores[i]] = abundance;
        }

    }

    void ModeReadRepresentation::add_duplicate(const numSeqs_t ab, const std::string& quality_scores, const QualityEncoding<>& qe) {

        abundance += ab;
        for (auto i = 0; i < quality_scores.size(); i++) {
            score_tallies[i][quality_scores[i]] += ab;
        }

    }

    std::string ModeReadRepresentation::get_scores(const QualityEncoding<>& qe) const {

        std::string scores(score_tallies.size(), ' ');
        for (auto i = 0; i < score_tallies.size(); i++) {
            scores[i] = std::max_element(
                    score_tallies[i].begin(),
                    score_tallies[i].end(),
                    [](const std::map<char, numSeqs_t>::value_type& a, const std::map<char, numSeqs_t>::value_type& b) {
                        return a.second < b.second;
                    })->first;
        }

        return scores;

    }

    /* === MedianReadRepresentation === */

    void MedianReadRepresentation::initialise(const Defline& dl, const std::string& seq, const std::string& quality_scores, const QualityEncoding<>& qe) {

        id = dl.id;
        abundance = dl.abundance;
        score_tallies = std::vector<std::map<char, numSeqs_t>>(quality_scores.length());
        for (auto i = 0; i < quality_scores.length(); i++) {
            score_tallies[i][quality_scores[i]] = abundance;
        }

    }

    void MedianReadRepresentation::add_duplicate(const numSeqs_t ab, const std::string& quality_scores, const QualityEncoding<>& qe) {

        abundance += ab;
        for (auto i = 0; i < quality_scores.size(); i++) {
            score_tallies[i][quality_scores[i]] += ab;
        }

    }

    std::string MedianReadRepresentation::get_scores(const QualityEncoding<>& qe) const {

        std::string scores(score_tallies.size(), ' ');

        if (abundance == 1) { // copy the single existing value

            for (auto i = 0; i < score_tallies.size(); i++) {
                scores[i] = score_tallies[i].begin()->first;
            }

        } else { // determine median from value-count pairs

            for (auto i = 0; i < score_tallies.size(); i++) {

                numSeqs_t pos = abundance / 2;
                numSeqs_t cnt = 0;
                for (auto& kv : score_tallies[i]) {

                    cnt += kv.second;
                    if (cnt > pos) {

                        scores[i] = kv.first;
                        break;

                    }

                }

            }

        }

        return scores;

    }

    /* === ProductReadRepresentation === */

    void ProductReadRepresentation::initialise(const Defline& dl, const std::string& seq, const std::string& quality_scores, const QualityEncoding<>& qe) {

        id = dl.id;
        abundance = dl.abundance;
        probability_products = std::vector<float>(quality_scores.length());
        for (auto i = 0; i < quality_scores.length(); i++) {
            probability_products[i] = std::pow(qe.encoded_quality_to_probability(quality_scores[i]), abundance);
        }

    }

    void ProductReadRepresentation::add_duplicate(const numSeqs_t ab, const std::string& quality_scores, const QualityEncoding<>& qe) {

        abundance += ab;
        for (auto i = 0; i < quality_scores.size(); i++) {
            probability_products[i] *= std::pow(qe.encoded_quality_to_probability(quality_scores[i]), ab);
        }

    }

    std::string ProductReadRepresentation::get_scores(const QualityEncoding<>& qe) const {

        float min_err_prob = qe.encoded_quality_to_probability(qe.get_accepted_scores().back());
        std::string scores(probability_products.size(), ' ');
        for (auto i = 0; i < probability_products.size(); i++) {
            scores[i] = qe.probability_to_encoded_quality(std::max(min_err_prob, probability_products[i]));
        }

        return scores;

    }

}