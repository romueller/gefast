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

#include <fstream>
#include <iomanip>
#include <sstream>

#include "../include/Distances.hpp"
#include "../include/OutputGenerators.hpp"

namespace GeFaST {

    /* === ClassicOutputGenerator === */

    void ClassicOutputGenerator::output_internal_structures(const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, std::vector<std::pair<const AmpliconCollection*, const Swarm*>>& swarm_order,
            const Configuration& config) {

        std::ofstream out_stream(config.output_internal);
        std::stringstream str_stream;

        numSeqs_t swarm_id = 0;
        for (auto& s : swarm_order) {

            const AmpliconCollection& ac = *(s.first);
            const Swarm* swarm = s.second;
            const Swarm* main_swarm = s.second;

            if (!swarm->is_attached()) {

                swarm_id++;
                GraftingInfo* graftings = swarm->get_grafting_info();
                numSeqs_t next = 0;

                do {

                    // swarm-internal links (1st round = main swarm, further rounds = grafted swarms)
                    for (numSeqs_t i = 1; i < swarm->size(); i++) {

                        str_stream << ac.id(swarm->parent(i)) << config.sep_internals
                            << ac.id(swarm->member(i)) << config.sep_internals
                            << swarm->parent_dist(i) << config.sep_internals
                            << swarm_id << config.sep_internals << swarm->gen(i) << std::endl;
                        out_stream << str_stream.rdbuf();
                        str_stream.str(std::string());

                    }

                    // if there are grafted swarms, first print the grafting link and then the internal links (above)
                    if (graftings != nullptr && next < graftings->num_links()) {

                        str_stream << ac.id(main_swarm->member(graftings->get_parent(next))) << config.sep_internals
                            << ac.id(graftings->get_child_swarm(next)->member(graftings->get_child(next))) << config.sep_internals
                            << graftings->get_distance(next) << config.sep_internals
                            << swarm_id << config.sep_internals << (main_swarm->gen(graftings->get_parent(next)) + 1) << std::endl;

                        swarm = graftings->get_child_swarm(next);

                        next++;

                    } else {
                        swarm = nullptr;
                    }

                } while (swarm != nullptr);

            }

        }

        out_stream.close();

        std::cout << " -- Internal-structures output created." << std::endl;

    }

    void ClassicOutputGenerator::output_otus(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
            std::vector<std::pair<const AmpliconCollection*, const Swarm*>>& swarm_order, const Configuration& config) {

        std::ofstream out_stream(config.output_otus);
        std::stringstream str_stream;

        for (auto& s : swarm_order) {

            const AmpliconCollection& ac = *(s.first);
            const Swarm* swarm = s.second;

            if (!swarm->is_attached()) {

                GraftingInfo* graftings = swarm->get_grafting_info();
                numSeqs_t next = 0;

                do {

                    // swarm members (1st round = main swarm, further rounds = grafted swarms)
                    str_stream << ac.id(swarm->seed()) << config.separator << ac.ab(swarm->seed());

                    for (numSeqs_t i = 1; i < swarm->size(); i++) {
                        str_stream << config.sep_otus << ac.id(swarm->member(i)) << config.separator << ac.ab(swarm->member(i));
                    }

                    str_stream << std::flush;
                    out_stream << str_stream.rdbuf();
                    str_stream.str(std::string());

                    // if there are grafted swarms, print additional separator
                    if (graftings != nullptr && next < graftings->num_links()) {

                        str_stream << config.sep_otus;

                        swarm = graftings->get_child_swarm(next);

                        next++;

                    } else {
                        swarm = nullptr;
                    }

                } while (swarm != nullptr);

                out_stream << std::endl;

            }

        }

        out_stream.close();

        std::cout << " -- Clusters output created." << std::endl;

    }

    void ClassicOutputGenerator::output_otus_mothur(const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, std::vector<std::pair<const AmpliconCollection*, const Swarm*>>& swarm_order,
            const Configuration& config) {

        std::ofstream out_stream(config.output_otus);
        std::stringstream str_stream;

        numSeqs_t num_swarms_adjusted = swarm_storage.num_swarms();
        for (auto& s : swarm_order) {
            num_swarms_adjusted -= s.second->is_attached();
        }
        out_stream << "swarm_" << config.main_threshold << '\t' << num_swarms_adjusted << config.sep_mothur_otu;


        for (auto& s : swarm_order) {

            const AmpliconCollection& ac = *(s.first);
            const Swarm* swarm = s.second;

            if (!swarm->is_attached()) {

                GraftingInfo* graftings = swarm->get_grafting_info();
                numSeqs_t next = 0;

                do {

                    // swarm members (1st round = main swarm, further rounds = grafted swarms)
                    str_stream << ac.id(swarm->seed()) << config.separator << ac.ab(swarm->seed());

                    for (numSeqs_t i = 1; i < swarm->size(); i++) {
                        str_stream << config.sep_mothur << ac.id(swarm->member(i)) << config.separator << ac.ab(swarm->member(i));
                    }

                    out_stream << str_stream.rdbuf();
                    str_stream.str(std::string());

                    // if there are grafted swarms, print additional separator
                    if (graftings != nullptr && next < graftings->num_links()) {

                        str_stream << config.sep_mothur;

                        swarm = graftings->get_child_swarm(next);

                        next++;

                    } else {
                        swarm = nullptr;
                    }

                } while (swarm != nullptr);

            }

        }

        out_stream << std::endl;

        out_stream.close();

        std::cout << " -- Clusters output created." << std::endl;

    }

    void ClassicOutputGenerator::output_statistics(const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, std::vector<std::pair<const AmpliconCollection*, const Swarm*>>& swarm_order,
            const Configuration& config) {

        std::ofstream out_stream(config.output_statistics);
        std::stringstream str_stream;

        for (auto& s : swarm_order) {

            const AmpliconCollection& ac = *(s.first);
            const Swarm& swarm = *(s.second);

            if (!swarm.is_attached()) {

                str_stream << swarm.total_num_different() << config.sep_statistics << swarm.total_mass() << config.sep_statistics
                    << ac.id(swarm.seed()) << config.sep_statistics << ac.ab(swarm.seed()) << config.sep_statistics
                    << swarm.total_num_singletons() << config.sep_statistics << swarm.max_gen() << config.sep_statistics
                    << swarm.max_rad() << std::endl;
                out_stream << str_stream.rdbuf();
                str_stream.str(std::string());

            }

        }

        out_stream.close();

        std::cout << " -- Statistics output created." << std::endl;

    }

    void ClassicOutputGenerator::output_seeds(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
            std::vector<std::pair<const AmpliconCollection*, const Swarm*>>& swarm_order, const Configuration& config) {

        std::ofstream out_stream(config.output_seeds);
        std::stringstream str_stream;

        for (auto& s : swarm_order) {

            const AmpliconCollection& ac = *(s.first);
            const Swarm& swarm = *(s.second);

            if (!swarm.is_attached()) {

                str_stream << '>' << ac.id(swarm.seed()) << config.separator << swarm.mass() << std::endl
                    << ac.seq(swarm.seed()) << std::endl;
                out_stream << str_stream.rdbuf();
                str_stream.str(std::string());

            }

        }

        out_stream.close();

        std::cout << " -- Seeds output created." << std::endl;

    }

    void ClassicOutputGenerator::output_uclust(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
            std::vector<std::pair<const AmpliconCollection*, const Swarm*>>& swarm_order, const Configuration& config) {

        std::ofstream out_stream(config.output_uclust);
        std::stringstream str_stream;
        str_stream << std::fixed << std::setprecision(1);

        ScoringFunction sf(config.match_reward, config.mismatch_penalty, config.gap_opening_penalty, config.gap_extension_penalty);
        lenSeqs_t d_row[amplicon_storage.max_length() + 1];
        lenSeqs_t p_row[amplicon_storage.max_length() + 1];
        char backtrack[(amplicon_storage.max_length() + 1) * (amplicon_storage.max_length() + 1)];

        numSeqs_t swarm_id = 0;
        for (auto& s : swarm_order) {

            const AmpliconCollection& ac = *(s.first);
            const Swarm* swarm = s.second;

            if (!swarm->is_attached()) {

                numSeqs_t seed = swarm->seed();
                GraftingInfo* graftings = swarm->get_grafting_info();
                numSeqs_t next = 0;

                str_stream << 'C' << config.sep_uclust << swarm_id << config.sep_uclust << swarm->total_size() << config.sep_uclust
                    << '*' << config.sep_uclust << '*' << config.sep_uclust << '*' << config.sep_uclust << '*' << config.sep_uclust
                    << '*' << config.sep_uclust << ac.id(swarm->seed()) << config.separator << ac.ab(swarm->seed()) << config.sep_uclust
                    << '*' << '\n';
                str_stream << 'S' << config.sep_uclust << swarm_id << config.sep_uclust << ac.len(swarm->seed()) << config.sep_uclust
                    << '*' << config.sep_uclust << '*' << config.sep_uclust << '*' << config.sep_uclust << '*' << config.sep_uclust
                    << '*' << config.sep_uclust << ac.id(swarm->seed()) << config.separator << ac.ab(swarm->seed()) << config.sep_uclust
                    << '*' << '\n';
                out_stream << str_stream.rdbuf();
                str_stream.str(std::string());

                do {

                    // H lines (1st round = main swarm, further rounds = grafted swarms)
                    for (numSeqs_t i = 1; i < swarm->size(); i++) {

                        auto ai = compute_gotoh_cigar_row(ac.seq(seed), ac.len(seed), ac.seq(swarm->member(i)),
                                ac.len(swarm->member(i)), sf, d_row, p_row, backtrack);

                        str_stream << 'H' << config.sep_uclust << swarm_id << config.sep_uclust
                            << ac.len(swarm->member(i)) << config.sep_uclust << (100.0 * (ai.length - ai.num_diffs) / ai.length)
                            << config.sep_uclust << '+' << config.sep_uclust << '0' << config.sep_uclust << '0' << config.sep_uclust
                            << ((ai.num_diffs == 0) ? "=" : ai.cigar) << config.sep_uclust
                            << ac.id(swarm->member(i)) << config.separator << ac.ab(swarm->member(i)) << config.sep_uclust
                            << ac.id(seed) << config.separator << ac.ab(seed) << '\n';
                        out_stream << str_stream.rdbuf();
                        str_stream.str(std::string());

                    }

                    // if there are grafted swarms, print H lines for their members but no C or S line (above)
                    if (graftings != nullptr && next < graftings->num_links()) {

                        swarm = graftings->get_child_swarm(next);

                        next++;

                    } else {
                        swarm = nullptr;
                    }

                } while (swarm != nullptr);

                swarm_id++;

            }

        }

        out_stream.close();

        std::cout << " -- UCLUST output created." << std::endl;

    }

    ClassicOutputGenerator* ClassicOutputGenerator::clone() const {
        return new ClassicOutputGenerator(*this);
    }

    void ClassicOutputGenerator::output_results(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
            const Configuration& config) {

        auto swarm_order = SwarmOrderer::order_by_seed_abundance(amplicon_storage, swarm_storage);

        if (!config.output_internal.empty()) output_internal_structures(amplicon_storage, swarm_storage, swarm_order, config);
        if (!config.output_otus.empty()) {
            if (config.mothur) {
                output_otus_mothur(amplicon_storage, swarm_storage, swarm_order, config);
            } else {
                output_otus(amplicon_storage, swarm_storage, swarm_order, config);
            }
        }
        if (!config.output_statistics.empty()) output_statistics(amplicon_storage, swarm_storage, swarm_order, config);
        if (!config.output_seeds.empty()) output_seeds(amplicon_storage, swarm_storage, swarm_order, config);
        if (!config.output_uclust.empty()) output_uclust(amplicon_storage, swarm_storage, swarm_order, config);

    }


    /* === DereplicationOutputGenerator === */

    DereplicationOutputGenerator* DereplicationOutputGenerator::clone() const {
        return new DereplicationOutputGenerator(*this);
    }

    void DereplicationOutputGenerator::output_dereplicate(const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, std::vector<std::pair<const AmpliconCollection*, const Swarm*>>& swarm_order,
            const Configuration& config) {

        std::ofstream out_internals, out_otus, out_statistics, out_seeds, out_uclust;
        std::stringstream str_internals, str_otus, str_statistics, str_seeds, str_uclust;
        str_uclust << std::fixed << std::setprecision(1);

        if (!config.output_internal.empty()) out_internals.open(config.output_internal);
        if (!config.output_otus.empty()) out_otus.open(config.output_otus);
        if (!config.output_statistics.empty()) out_statistics.open(config.output_statistics);
        if (!config.output_seeds.empty()) out_seeds.open(config.output_seeds);
        if (!config.output_uclust.empty()) out_uclust.open(config.output_uclust);

        if (!config.output_otus.empty() && config.mothur) out_otus << "swarm_" << config.main_threshold << '\t' << swarm_storage.num_swarms();

        numSeqs_t swarm_id = 0;
        for (auto& s : swarm_order) {

            const AmpliconCollection& ac = *(s.first);
            const Swarm& swarm = *(s.second);

            if (!config.output_internal.empty()) {

                for (numSeqs_t i = 1; i < swarm.size(); i++) {

                    str_internals << ac.id(swarm.seed()) << config.sep_internals << ac.id(swarm.member(i)) << config.sep_internals
                        << 0 << config.sep_internals << (swarm_id + 1) << config.sep_internals << 0 << std::endl;
                    out_internals << str_internals.rdbuf();
                    str_internals.str(std::string());

                }

            }

            if (!config.output_otus.empty()) {

                if (config.mothur) {

                    str_otus << config.sep_mothur_otu << ac.id(swarm.seed()) << config.separator << ac.ab(swarm.seed());

                    for (numSeqs_t i = 1; i < swarm.size(); i++) {
                        str_otus << config.sep_mothur << ac.id(swarm.member(i)) << config.separator << ac.ab(swarm.member(i));
                    }

                    out_otus << str_otus.rdbuf();
                    str_otus.str(std::string());

                } else {

                    str_otus << ac.id(swarm.seed()) << config.separator << ac.ab(swarm.seed());

                    for (numSeqs_t i = 1; i < swarm.size(); i++) {
                        str_otus << config.sep_otus << ac.id(swarm.member(i)) << config.separator << ac.ab(swarm.member(i));
                    }

                    str_otus << std::endl;
                    out_otus << str_otus.rdbuf();
                    str_otus.str(std::string());

                }

            }

            if (!config.output_statistics.empty()) {

                str_statistics << swarm.num_different() << config.sep_statistics << swarm.mass() << config.sep_statistics
                    << ac.id(swarm.seed()) << config.sep_statistics << ac.ab(swarm.seed()) << config.sep_statistics
                    << swarm.num_singletons() << config.sep_statistics << 0 << config.sep_statistics << 0 << std::endl;
                out_statistics << str_statistics.rdbuf();
                str_statistics.str(std::string());

            }

            if (!config.output_seeds.empty()) {

                str_seeds << '>' << ac.id(swarm.seed()) << config.separator << swarm.mass() << std::endl <<
                    ac.seq(swarm.seed()) << std::endl;
                out_seeds << str_seeds.rdbuf();
                str_seeds.str(std::string());

            }

            if (!config.output_uclust.empty()) {

                str_uclust << 'C' << config.sep_uclust << swarm_id << config.sep_uclust << swarm.size() << config.sep_uclust
                    << '*' << config.sep_uclust << '*' << config.sep_uclust << '*' << config.sep_uclust << '*' << config.sep_uclust
                    << '*' << config.sep_uclust << ac.id(swarm.seed()) << config.separator << ac.ab(swarm.seed()) << config.sep_uclust
                    << '*' << '\n';
                str_uclust << 'S' << config.sep_uclust << swarm_id << config.sep_uclust << ac.len(swarm.seed()) << config.sep_uclust
                    << '*' << config.sep_uclust << '*' << config.sep_uclust << '*' << config.sep_uclust << '*' << config.sep_uclust
                    << '*' << config.sep_uclust << ac.id(swarm.seed()) << config.separator << ac.ab(swarm.seed()) << config.sep_uclust
                    << '*' << '\n';
                out_uclust << str_uclust.rdbuf() << std::flush;
                str_uclust.str(std::string());

                for (numSeqs_t i = 1; i < swarm.size(); i++) {

                    str_uclust << 'H' << config.sep_uclust << swarm_id << config.sep_uclust
                        << ac.len(swarm.member(i)) << config.sep_uclust << "100.0" << config.sep_uclust
                        << '+' << config.sep_uclust << '0' << config.sep_uclust << '0' << config.sep_uclust
                        << '=' << config.sep_uclust << ac.id(swarm.member(i)) << config.separator
                        << ac.ab(swarm.member(i)) << config.sep_uclust << ac.id(swarm.seed()) << config.separator
                        << ac.ab(swarm.seed()) << '\n';

                    out_uclust << str_uclust.rdbuf() << std::flush;
                    str_uclust.str(std::string());

                }


            }

        }

        if (!config.output_otus.empty() && config.mothur) out_otus << std::endl;

        if (!config.output_internal.empty()) out_internals.close();
        if (!config.output_otus.empty()) out_otus.close();
        if (!config.output_statistics.empty()) out_statistics.close();
        if (!config.output_seeds.empty()) out_seeds.close();
        if (!config.output_uclust.empty()) out_uclust.close();

        std::cout << " -- Outputs created." << std::endl;

    }

    void DereplicationOutputGenerator::output_results(const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, const Configuration& config) {

        auto swarm_order = SwarmOrderer::order_by_mass(amplicon_storage, swarm_storage);

        output_dereplicate(amplicon_storage, swarm_storage, swarm_order, config);

    }


    /* === IdleOutputGenerator === */

    IdleOutputGenerator* IdleOutputGenerator::clone() const {
        return new IdleOutputGenerator(*this);
    }

    void IdleOutputGenerator::output_results(const AmpliconStorage& amplicon_storage,
            const SwarmStorage& swarm_storage, const Configuration& config) {
        // do nothing
    }

}