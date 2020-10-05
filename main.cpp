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

#include "include/Factories.hpp"

namespace GeFaST {

    /*
     * Controls the overall workflow of GeFaST:
     * 0) read / set up configuration
     * 1) reading and preprocessing input data
     * 2) clustering amplicons
     * 3) cluster refinement (optional)
     * 4) output generation
     */
    int run(int argc, const char* argv[]) {

        print_information();


        /* ===== Configuration ===== */

        if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)) {

            print_help();
            return 0;

        }

        if (argc < 4) {

            std::cerr << "ERROR: Not enough arguments." << std::endl;
            print_help();
            return 1;

        }

        Configuration* config = ConfigurationFactory::create(argc, argv);

        if (config->input_files.empty()) {

            std::cerr << "ERROR: No input file(s) specified." << std::endl;
            print_help();
            delete config;
            return 1;

        }
        if ((config->output_internal + config->output_otus + config->output_statistics + config->output_seeds +
                config->output_uclust + config->keep_preprocessed).empty()) {

            std::cerr << "ERROR: No output file(s) specified." << std::endl;
            print_help();
            delete config;
            return 1;

        }

        config->print(std::cout);

        std::cout << "===== Computation =====\n" << std::endl;

        /* ===== Preprocessing ===== */

        std::cout << "Preprocessing inputs..." << std::endl;
        Preprocessor* preprocessor = config->build_preprocessor();
        AmpliconStorage* amplicon_storage = preprocessor->preprocess_inputs(config->input_files, *config);

        if (!config->keep_preprocessed.empty()) {
            amplicon_storage->print(config->keep_preprocessed, *config);
        }


        if (amplicon_storage->num_amplicons() > 0) {

            std::cout << " >> Number of amplicons: " << amplicon_storage->num_amplicons() << std::endl;

            /* ===== Clustering ===== */

            Clusterer* clusterer = config->build_clusterer();
            SwarmStorage* swarm_storage = clusterer->cluster(*amplicon_storage, *config);

            std::cout << " >> Number of swarms: " << swarm_storage->num_swarms() << std::endl;

            /* ===== Cluster refinement ===== */

            if (config->refinement_threshold > 0) {

                std::cout << "Refining swarms..." << std::endl;
                ClusterRefiner* refiner = config->build_cluster_refiner();
                refiner->refine(*amplicon_storage, *swarm_storage, *config);
                delete refiner;

                // determine number of amplicons from swarm storage (amplicons not added to / removed from swarms
                // during refinement cannot be figured from the amplicon storage)
                numSeqs_t num_ampl = 0;
                for (auto p = 0; p < swarm_storage->num_pools(); p++) {
                    auto& sws = swarm_storage->get_swarms(p);
                    for (auto s = 0; s < sws.size(); s++) {
                        auto& sw = sws.get(s);
                        if (!sw.is_attached()) num_ampl += sw.total_size();
                    }
                }
                std::cout << " >> Number of amplicons after refinement: " << num_ampl << std::endl;
                std::cout << " >> Number of swarms after refinement: " << swarm_storage->num_swarms() << std::endl;

            }


            /* ===== Output generation ===== */

            std::cout << "Generating outputs..." << std::endl;
            OutputGenerator* generator = config->build_output_generator();
            generator->output_results(*amplicon_storage, *swarm_storage, *config);

            delete generator;
            delete clusterer;
            delete swarm_storage;

        } else {
            std::cout << " >> No amplicons remained after the preprocessing, i.e. there is nothing to cluster." << std::endl;
        }


        /* ===== Cleaning up ===== */

        std::cout << "Cleaning up..." << std::endl;
        delete amplicon_storage;
        delete preprocessor;
        delete config;
        std::cout << "Computation finished." << std::endl;

        return 0;

    }

}


int main(int argc, const char* argv[]) {
    GeFaST::run(argc, argv);
}
