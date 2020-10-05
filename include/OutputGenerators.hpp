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

#ifndef GEFAST_OUTPUTGENERATORS_HPP
#define GEFAST_OUTPUTGENERATORS_HPP

#include "Base.hpp"

namespace GeFaST {

    /*
     * Output generator mimicking the behaviour of Swarm.
     * The output consists of up to five different files, which can be selected independently.
     * The specifics of the different files are described above the respective method.
     */
    class ClassicOutputGenerator : public OutputGenerator {

    public:
        ClassicOutputGenerator* clone() const override; // deep-copy clone method

        /*
         * Output the determined swarms as requested in the configuration.
         */
        void output_results(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                const Configuration& config) override;

    private:
        /*
         * Write the links of the given swarms to file (corresponds to output of Swarm's option -i).
         * Each line contains one link represented through the
         * (1) amplicon id of the parent,
         * (2) amplicon id of the child,
         * (3) distance,
         * (4) cluster id, and
         * (5) generation number of the child amplicon,
         * separated by the given separator.
         */
        void output_internal_structures(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                std::vector<std::pair<const AmpliconCollection*, const Swarm*>>& swarm_order, const Configuration& config);

        /*
         * Write the members of the given clusters to file (corresponds to output of Swarm's option -o).
         * Each line contains the members of one cluster represented through amplicon id and abundance.
         * The members are separated via sep, while id and abundance are separated via sepAbundance.
         */
        void output_otus(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                std::vector<std::pair<const AmpliconCollection*, const Swarm*>>& swarm_order, const Configuration& config);

        /*
         * Write the members of the given clusters to file (corresponds to output of Swarm's option -o with -r).
         * On a single line, the members of all clusters are represented through amplicon id and abundance.
         * The members of one cluster and the clusters themselves are separated via sep resp. sepOtu, while id and abundance are separated via sepAbundance.
         */
        void output_otus_mothur(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                std::vector<std::pair<const AmpliconCollection*, const Swarm*>>& swarm_order, const Configuration& config);

        /*
         * Write the statistics of the given clusters to file (corresponds to output of Swarm's option -s).
         * Each line contains the statistics of one cluster represented through the
         * (1) number of unique sequences,
         * (2) mass of the cluster,
         * (3) amplicon id of the seed,
         * (4) abundance of the seed amplicon,
         * (5) number of singletons,
         * (6) number of iterations before the cluster reached its natural limits, and
         * (7) number of cumulated differences between the seed and the furthermost amplicon,
         * separated by the given separator.
         * The statistics (6) and (7) are not affected by the refinement phase.
         */
        void output_statistics(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                std::vector<std::pair<const AmpliconCollection*, const Swarm*>>& swarm_order, const Configuration& config);

        /*
         * Write the seeds of the given clusters to file (corresponds to output of Swarm's option -w).
         * Each seed comprises two lines.
         * The first line contains the amplicon id of the seed preceded by '>' and followed by a separator and the mass of the cluster.
         * The second line describes the sequence of the seed amplicon.
         */
        void output_seeds(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                std::vector<std::pair<const AmpliconCollection*, const Swarm*>>& swarm_order, const Configuration& config);

        /*
         * Write the clustering results in a uclust-like format to file (corresponds to output of Swarm's option -u).
         */
        void output_uclust(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                std::vector<std::pair<const AmpliconCollection*, const Swarm*>>& swarm_order, const Configuration& config);

    };

    /*
     * Output generator mimicking the behaviour of Swarm when dereplicating.
     * The output consists of up to five different files, which can be selected independently.
     * The specifics of the different files are described above the respective method.
     */
    class DereplicationOutputGenerator : public OutputGenerator {

    public:
        DereplicationOutputGenerator* clone() const override; // deep-copy clone method

        /*
         * Output the determined swarms as requested in the configuration.
         */
        void output_results(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                const Configuration& config) override;

    private:
        /*
         * Write the requested dereplication outputs to file (SwarmConfig stores information on which are requested).
         * The outputs correspond to Swarm's outputs (-i, -o, -s, -w, -u) when running it with -d 0.
         *
         * Minor differences to non-dereplicating output options:
         * -i: generation of child amplicons is also (always) 0
         * -o: <no differences>
         * -s: first column describes the number of amplicons in the cluster (not the number of unique sequences)
         * -w: <no differences>
         * -u: <no differences>
         */
        void output_dereplicate(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                std::vector<std::pair<const AmpliconCollection*, const Swarm*>>& swarm_order, const Configuration& config);

    };

    /*
     * Output generator that does nothing.
     *
     * Used in preprocessing-only (prep) mode.
     */
    class IdleOutputGenerator : public OutputGenerator {

    public:
        IdleOutputGenerator* clone() const override; // deep-copy clone method

        /*
         * Do not generate any output.
         */
        void output_results(const AmpliconStorage& amplicon_storage, const SwarmStorage& swarm_storage,
                const Configuration& config) override;

    };

}

#endif //GEFAST_OUTPUTGENERATORS_HPP
