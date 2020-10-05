/*
 * GeFaST
 *
 * Copyright (C) 2016 - 2020 Robert Mueller
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
                        ds.record(dl.id.length(), seq.length());
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
                ds.record(dl.id.length(), seq.length());
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


    /* === FastqPreprocessor === */

    FastqPreprocessor::ReadRepresentation::ReadRepresentation(const Defline& dl, const std::string& seq,
            const std::vector<float>& quality_scores) {

        id = dl.id;
        abundance = dl.abundance;
        probability_sums = std::vector<float>(seq.length(), 0.0f);
        for (auto i = 0; i < probability_sums.size(); i++) {
            probability_sums[i] += abundance * quality_scores[i];
        }

    }

    void FastqPreprocessor::ReadRepresentation::add_duplicate(const numSeqs_t ab, const std::vector<float>& quality_scores) {

        abundance += ab;
        for (auto i = 0; i < probability_sums.size(); i++) {
            probability_sums[i] += ab * quality_scores[i];
        }

    }


    void FastqPreprocessor::read_input_file(const std::string& file, const Configuration& config, DataStatistics<>& ds,
            std::map<std::string, ReadRepresentation>& different_reads) {

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


        Defline dl;
        std::string line, seq;

        while (std::getline(in_stream, line).good()) {

            if (line.empty() || line[0] == ';') continue; // skip empty and comment lines (begin with ';')

            if (line[0] == '@') { // header line

                dl = Defline::parse_description_line(line, config.separator);

                std::getline(in_stream, seq); // read sequence
                upper_case(seq);

                std::getline(in_stream, line); // read '+' line (ignored)
                std::getline(in_stream, line); // read quality scores

                if (check_sequence(seq, config.alphabet, dl.abundance, min_length, max_length, min_abundance, max_abundance)) {

                    std::vector<float> quality_scores(seq.length(), 0.0f); // error probability per position corresponding to quality score
                    for (auto i = 0; i < seq.length(); i++) {
                        quality_scores[i] += quality_encoding_.encoded_quality_to_probability(line[i]);
                    }

                    auto iter = different_reads.find(seq);
                    if (iter != different_reads.end()) {
                        iter->second.add_duplicate(dl.abundance, quality_scores);
                    } else {

                        different_reads[seq] = ReadRepresentation(dl, seq, quality_scores);
                        ds.record(dl.id.length(), seq.length());

                    }

                }

            }

        }

    }

    FastqPreprocessor::FastqPreprocessor(const QualityEncoding<>& qe) : quality_encoding_(qe) {
        // nothing else to do
    }

    FastqPreprocessor::FastqPreprocessor(const QualityEncodingOption opt) : quality_encoding_(QualityEncoding<>(opt)) {
        // nothing else to do
    }

    FastqPreprocessor* FastqPreprocessor::clone() const {
        return new FastqPreprocessor(*this);
    }

    AmpliconStorage* FastqPreprocessor::preprocess_inputs(const std::vector<std::string>& files, const Configuration& config) {

        AmpliconStorage* amplicon_storage = nullptr;
        std::map<std::string, ReadRepresentation> different_reads;

        {
            DataStatistics<> ds;

            std::cout << " -- Scanning input files..." << std::endl;
            for (auto& f : files) {
                read_input_file(f, config, ds, different_reads);
            }

            if (config.filter_abundance > 0) {
                filter_by_abundance(different_reads, ds, config);
            }

            amplicon_storage = config.build_amplicon_storage(ds);

            std::cout << " -- Scan completed." << std::endl;

        }

        std::cout << " -- Inserting amplicons..." << std::endl;
        for (auto iter = different_reads.begin(); iter != different_reads.end(); iter++) {

            auto& rr = iter->second;
            ExtraInfoQuality ei;
            for (auto i = 0; i < rr.probability_sums.size(); i++) {
                ei.quality_scores.push_back(quality_encoding_.probability_to_encoded_quality(rr.probability_sums[i] / rr.abundance));
            }
            Defline dl(rr.id, "", rr.abundance);
            amplicon_storage->add_amplicon(dl, iter->first, ei);

        }

        amplicon_storage->finalise();

        std::cout << " -- Construction of amplicon storage completed." << std::endl;

        return amplicon_storage;

    }

    void FastqPreprocessor::filter_by_abundance(std::map<std::string, ReadRepresentation>& different_reads,
            DataStatistics<>& ds, const Configuration& config) {

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

        for (auto iter = different_reads.begin(); iter != different_reads.end(); ) {

            if (iter->second.abundance < min_abundance || iter->second.abundance > max_abundance) {

                ds.unrecord(iter->second.id.length(), iter->first.length());
                iter = different_reads.erase(iter);

            } else {
                iter++;
            }

        }

    }

}