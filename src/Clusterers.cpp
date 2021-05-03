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

#include <cstring>
#include <fstream>

#include "../include/Clusterers.hpp"

namespace GeFaST {

    ClassicSwarmer* ClassicSwarmer::clone() const {
        return new ClassicSwarmer(*this);
    }

    SwarmStorage* ClassicSwarmer::cluster(const AmpliconStorage& amplicon_storage, const Configuration& config) {

        std::cout << "Clustering amplicons into swarms..." << std::endl;

        SwarmStorage* swarm_storage = config.build_swarm_storage(amplicon_storage);

        for (numSeqs_t p = 0; p < amplicon_storage.num_pools(); p++) {

            auto& ac = amplicon_storage.get_pool(p);
            auto& swarms = swarm_storage->get_swarms(p);

            AuxiliaryData* aux = config.build_auxiliary_data(amplicon_storage, p, config.main_threshold);
            std::vector<bool> swarmed(ac.size(), false); // swarmed amplicons are already included in a cluster

            // open new swarm for the amplicon with the highest abundance that is not yet included in a swarm
            // (corresponds to next unswarmed amplicon due to sorting of amplicons in the pool)
            for (numSeqs_t seed = 0; seed < ac.size(); seed++) {

                if (!swarmed[seed]) {

                    /* (a) Initialise new swarm with seed */
                    auto& cur_swarm = swarms.initialise_cluster(seed, ac.ab(seed));
                    swarmed[seed] = true;
                    aux->record_amplicon(seed);
                    aux->tick_off_amplicon(seed);
                    SwarmEntry* cur_seed = cur_swarm.seed_entry();

                    lenSeqs_t last_gen = 0;

                    /* (b) Proceed breadth-first through the 'amplicon space' */
                    numSeqs_t pos = 0;
                    while (pos < cur_swarm.size()) { // expand current swarm until no further similar amplicons can be added

                        if (last_gen != cur_swarm.gen(pos)) { // sort upcoming generation to work through it by decreasing abundance

                            cur_swarm.sort_next_generation(ac, pos);
                            aux->clear_amplicon_records();

                        }

                        // get next swarm (sub)seed
                        cur_seed = cur_swarm.get_entry(pos, cur_seed);
                        numSeqs_t cur_seed_gen = cur_seed->gen();
                        dist_t cur_seed_rad = cur_seed->rad();
                        numSeqs_t cur_seed_mem = cur_seed->member();

                        // Consider yet unseen (unswarmed) amplicons to continue the exploration.
                        // An amplicon is marked as 'swarmed' as soon as it occurs the first time as a partner
                        // in order to prevent the algorithm from queueing it more than once coming from different amplicons.
                        for (auto& partner : aux->find_partners(cur_seed_mem)) {

                            // append new member
                            cur_swarm.append(
                                    partner.id, // pool-internal integer id of new member
                                    cur_seed_gen + 1, // generation number of new member
                                    cur_seed_rad + partner.dist, // radius of new member
                                    cur_seed_mem, // pool-internal integer id of parent
                                    partner.dist, // distance to parent
                                    ac.ab(partner.id), // abundance of new member
                                    (partner.dist != 0) && aux->record_amplicon(partner.id) // flag indicating whether new member is also a new sequence
                                    );

                            swarmed[partner.id] = true;
                            aux->tick_off_amplicon(partner.id);

                        }

                        last_gen = cur_seed_gen;
                        pos++;

                    }

                    /* (c) Close the no longer extendable swarm */
                    aux->clear_amplicon_records();
                    delete cur_seed;

                }

            }

            delete aux;

        }

        return swarm_storage;

    }


    Dereplicator* Dereplicator::clone() const {
        return new Dereplicator(*this);
    }

    SwarmStorage* Dereplicator::cluster(const AmpliconStorage& amplicon_storage, const Configuration& config) {

        struct lessCharArray {
            bool operator()(const char* lhs, const char* rhs) const {
                return strcmp(lhs, rhs) < 0;
            }
        };

        std::cout << "Dereplicating amplicons..." << std::endl;

        SwarmStorage* swarm_storage = config.build_swarm_storage(amplicon_storage);

        for (numSeqs_t p = 0; p < amplicon_storage.num_pools(); p++) {

            auto& ac = amplicon_storage.get_pool(p);
            std::map<const char*, std::vector<numSeqs_t>, lessCharArray> groups;

            // group amplicons based on their sequence
            for (numSeqs_t i = 0; i < ac.size(); i++) {
                groups[ac.seq(i)].push_back(i);
            }

            auto& swarms = swarm_storage->get_swarms(p);

            // output each group as a swarm
            for (auto& g : groups) {

                auto& members = g.second;
                numSeqs_t seed = members[0];

                auto& cur_swarm = swarms.initialise_cluster(seed, ac.ab(seed));
                for (auto i = 1; i < members.size(); i++) {
                    cur_swarm.append(members[i], 0, 0, seed, 0, ac.ab(members[i]), true);
                }

            }

        }

        return swarm_storage;

    }


    IdleClusterer* IdleClusterer::clone() const {
        return new IdleClusterer(*this);
    }

    SwarmStorage* IdleClusterer::cluster(const AmpliconStorage& amplicon_storage, const Configuration& config) {
        return config.build_swarm_storage(amplicon_storage);
    }


    Dada2Clusterer::Dada2Clusterer(const Configuration& config) {

        auto& misc = config.misc;
        auto iter = misc.begin();

        if ((iter = misc.find("use_k_mers")) != misc.end()) {
            use_k_mers = iter->second == "true";
        }
        if ((iter = misc.find("k_dist_cutoff")) != misc.end()) {
            k_dist_cutoff = std::stod(iter->second);
        }
        if ((iter = misc.find("band_size")) != misc.end()) {
            band_size = std::stoi(iter->second);
        }
        if ((iter = misc.find("omega_a")) != misc.end()) {
            omega_a = std::stod(iter->second);
        }
        if ((iter = misc.find("omega_p")) != misc.end()) {
            omega_p = std::stod(iter->second);
        }
        if ((iter = misc.find("max_clust")) != misc.end()) {
            max_clust = std::stoi(iter->second);
        }
        if ((iter = misc.find("min_fold")) != misc.end()) {
            min_fold = std::stod(iter->second);
        }
        if ((iter = misc.find("min_hamming")) != misc.end()) {
            min_hamming = std::stoi(iter->second);
        }
        if ((iter = misc.find("min_abund")) != misc.end()) {
            min_abund = std::stoi(iter->second);
        }
        if ((iter = misc.find("homo_gap")) != misc.end()) {
            homo_gap = std::stoi(iter->second);
        }
        if ((iter = misc.find("gapless")) != misc.end()) {
            gapless = iter->second == "true";
        }
        if ((iter = misc.find("greedy")) != misc.end()) {
            greedy = iter->second == "true";
        }
        if ((iter = misc.find("error_matrix")) != misc.end()) {
            matrix_file = iter->second;
        }

    }

    Dada2Clusterer* Dada2Clusterer::clone() const {
        return new Dada2Clusterer(*this);
    }

    SwarmStorage* Dada2Clusterer::cluster(const AmpliconStorage& amplicon_storage, const Configuration& config) {

        std::cout << "Clustering amplicons into swarms..." << std::endl;

        homo_gap = (homo_gap == 0) ? (config.gap_opening_penalty + config.gap_extension_penalty) : homo_gap;
        lenSeqs_t alphabet_size = (config.alphabet.empty()) ? 4 : config.alphabet.size();

        // set up error matrix
        QualityEncoding<>* qe = config.build_quality_encoding();
        char score_shift = qe->get_offset() + qe->get_from(); // subtracted from quality scores before accessing the error matrix
        std::vector<std::vector<double>> error_matrix;
        if (!Dada2Utility::fill_error_matrix(error_matrix, matrix_file, qe->get_accepted_scores().size())) {

            for (auto c1 = 0; c1 < alphabet_size; c1++) {
                for (auto c2 = 0; c2 < alphabet_size; c2++) {
                    error_matrix.emplace_back();
                    for (auto q : qe->get_accepted_scores()) {
                        auto err_prob = qe->encoded_quality_to_probability(q);
                        error_matrix[c1 * alphabet_size + c2].push_back((c1 == c2) ? (1 - err_prob) : (err_prob / (alphabet_size - 1)));
                    }
                }
            }

            std::cout << " -- Error matrix filled automatically based on the quality-encoding scheme." << std::endl;

        }

        auto swarm_storage = config.build_swarm_storage(amplicon_storage);

        for (numSeqs_t p = 0; p < amplicon_storage.num_pools(); p++) {

            auto& ac = amplicon_storage.get_pool(p);
            auto& swarms = swarm_storage->get_swarms(p);

            lenSeqs_t max_len = ac.len(0);
            lenSeqs_t min_len = ac.len(0);
            for (numSeqs_t i = 1; i < ac.size(); i++) {

                max_len = std::max(max_len, ac.len(i));
                min_len = std::min(min_len, ac.len(i));

            }

            if (min_len <= Dada2Utility::k_mer_size) throw std::runtime_error("Input sequences must all be longer than the k-mer size (" + std::to_string(Dada2Utility::k_mer_size) + ").");

            // construct a raw for each input sequence, store in raws[index]
            char seq[max_len + 1];
            uint8_t qual[max_len];
            RawSequence** raws = new RawSequence*[ac.size()];
            for (numSeqs_t i = 0; i < ac.size(); i++) {

                Dada2Utility::nt2int(seq, ac.seq(i));

                auto quals = ac.quals(i);
                for (lenSeqs_t pos = 0; pos < ac.len(i); pos++) {
                    qual[pos] = (uint8_t)(quals[pos] - score_shift);
                }
                raws[i] = new RawSequence(seq, ac.len(i), qual, ac.ab(i), false);

                raws[i]->index = i;

            }

            uint8_t* k8;
            uint16_t* k16;
            uint16_t* kord;

            if (use_k_mers) {

                // add uint8_t k-mer index in contiguous memory block
                size_t n_kmer = 1 << (2 * Dada2Utility::k_mer_size);
                k8 = new uint8_t[ac.size() * n_kmer];
                for (numSeqs_t index = 0; index < ac.size(); index++) {

                    RawSequence* raw = raws[index];
                    raw->kmer8 = &k8[index * n_kmer];
                    Dada2Utility::assign_kmer8(raw->kmer8, raw->seq, Dada2Utility::k_mer_size);

                }

                // add uint16_t k-mer index in contiguous memory block
                k16 = new uint16_t[ac.size() * n_kmer];
                for (numSeqs_t index = 0; index < ac.size(); index++) {

                    RawSequence* raw = raws[index];
                    raw->kmer = &k16[index * n_kmer];
                    Dada2Utility::assign_kmer(raw->kmer, raw->seq, Dada2Utility::k_mer_size);

                }

                // add uint16_t ordered k-mer record in contiguous memory block
                kord = new uint16_t[ac.size() * max_len];
                for (numSeqs_t index = 0; index < ac.size(); index++) {

                    RawSequence* raw = raws[index];
                    raw->kord = &kord[index * max_len];
                    Dada2Utility::assign_kmer_order(raw->kord, raw->seq, Dada2Utility::k_mer_size);

                }

            }

            ClusterCollection* clusters = run_dada(raws, ac.size(), error_matrix, config.match_reward,
                    config.mismatch_penalty, config.gap_opening_penalty + config.gap_extension_penalty, homo_gap,
                    use_k_mers, k_dist_cutoff, band_size, omega_a, omega_p, max_clust, min_fold, min_hamming, min_abund,
                    gapless, greedy);

            for (auto i = 0; i < clusters->num_clusters; i++) {

                Cluster* c = clusters->clusters[i];
                auto& new_swarm = swarms.initialise_cluster(c->centre->index, c->centre->abundance);
                for (auto j = 0; j < c->size; j++) {
                    if (c->raws[j]->index != c->centre->index) new_swarm.append(c->raws[j]->index, 1, config.main_threshold + 1, c->centre->index, config.main_threshold + 1, c->raws[j]->abundance, true);
                }

            }

            // free memory
            delete clusters;
            for (numSeqs_t index = 0; index < ac.size(); index++) {
                delete raws[index];
            }
            delete[] raws;
            if (use_k_mers) {

                delete[] k8;
                delete[] k16;
                delete[] kord;

            }

        }

        delete qe;

        return swarm_storage;

    }

    void Dada2Clusterer::add_swarms(const AmpliconCollection& ac, std::vector<numSeqs_t>& ampl_ids, Swarms& swarms, const Configuration& config,
                                    bool use_k_mers, double k_dist_cutoff, int band_size, double omega_a,
                                    double omega_p, int max_clust, double min_fold, int min_hamming, int min_abund,
                                    int homo_gap, bool gapless, bool greedy, std::string matrix_file) {

        homo_gap = (homo_gap == 0) ? (config.gap_opening_penalty + config.gap_extension_penalty) : homo_gap;
        lenSeqs_t alphabet_size = (config.alphabet.empty()) ? 4 : config.alphabet.size();

        if (ampl_ids.empty()) return; // nothing to do

        lenSeqs_t max_len = ac.len(ampl_ids[0]);
        lenSeqs_t min_len = ac.len(ampl_ids[0]);

        for (numSeqs_t i = 1; i < ampl_ids.size(); i++) {

            max_len = std::max(max_len, ac.len(ampl_ids[i]));
            min_len = std::min(min_len, ac.len(ampl_ids[i]));

        }
        if (min_len <= Dada2Utility::k_mer_size) throw std::runtime_error("Input sequences must all be longer than the k-mer size (" + std::to_string(Dada2Utility::k_mer_size) + ").");


        // set up error matrix
        QualityEncoding<>* qe = config.build_quality_encoding();
        char score_shift = qe->get_offset() + qe->get_from(); // subtracted from quality scores before accessing the error matrix
        std::vector<std::vector<double>> error_matrix;
        if (!Dada2Utility::fill_error_matrix(error_matrix, matrix_file, qe->get_accepted_scores().size())) {

            for (auto c1 = 0; c1 < alphabet_size; c1++) {
                for (auto c2 = 0; c2 < alphabet_size; c2++) {
                    error_matrix.emplace_back();
                    for (auto q : qe->get_accepted_scores()) {
                        auto err_prob = qe->encoded_quality_to_probability(q);
                        error_matrix[c1 * alphabet_size + c2].push_back((c1 == c2) ? (1 - err_prob) : (err_prob / (alphabet_size - 1)));
                    }
                }
            }

        }


        char seq[max_len + 1];
        uint8_t qual[max_len];
        RawSequence** raws = new RawSequence*[ampl_ids.size()];

        // construct a raw for each input sequence, store in raws[index]
        for (numSeqs_t i = 0; i < ampl_ids.size(); i++) {

            Dada2Utility::nt2int(seq, ac.seq(ampl_ids[i]));

            auto quals = ac.quals(ampl_ids[i]);
            for (unsigned int pos = 0; pos < ac.len(ampl_ids[i]); pos++) {
                qual[pos] = (uint8_t)(quals[pos] - score_shift);
            }
            raws[i] = new RawSequence(seq, ac.len(ampl_ids[i]), qual, ac.ab(ampl_ids[i]), false);

            raws[i]->index = i;

        }

        uint8_t* k8;
        uint16_t* k16;
        uint16_t* kord;

        if (use_k_mers) {

            // add uint8_t k-mer index in contiguous memory block
            size_t n_kmer = 1 << (2 * Dada2Utility::k_mer_size);
            k8 = new uint8_t[ampl_ids.size() * n_kmer];
            for (numSeqs_t index = 0; index < ampl_ids.size(); index++) {

                RawSequence* raw = raws[index];
                raw->kmer8 = &k8[index * n_kmer];
                Dada2Utility::assign_kmer8(raw->kmer8, raw->seq, Dada2Utility::k_mer_size);

            }

            // add uint16_t k-mer index in contiguous memory block
            k16 = new uint16_t[ampl_ids.size() * n_kmer];
            for (numSeqs_t index = 0; index < ampl_ids.size(); index++) {

                RawSequence* raw = raws[index];
                raw->kmer = &k16[index * n_kmer];
                Dada2Utility::assign_kmer(raw->kmer, raw->seq, Dada2Utility::k_mer_size);

            }

            // add uint16_t ordered k-mer record in contiguous memory block
            kord = new uint16_t[ampl_ids.size() * max_len];
            for (numSeqs_t index = 0; index < ampl_ids.size(); index++) {

                RawSequence* raw = raws[index];
                raw->kord = &kord[index * max_len];
                Dada2Utility::assign_kmer_order(raw->kord, raw->seq, Dada2Utility::k_mer_size);

            }

        }

        ClusterCollection* clusters = run_dada(raws, ampl_ids.size(), error_matrix, config.match_reward,
                config.mismatch_penalty, config.gap_opening_penalty + config.gap_extension_penalty, homo_gap,
                use_k_mers, k_dist_cutoff, band_size, omega_a, omega_p, max_clust, min_fold, min_hamming, min_abund,
                gapless, greedy);


        for (auto i = 0; i < clusters->num_clusters; i++) {

            Cluster* c = clusters->clusters[i];
            Swarm& new_swarm = swarms.initialise_cluster(ampl_ids[c->centre->index], c->centre->abundance);
            for (auto j = 0; j < c->size; j++) {
                if (c->raws[j]->index != c->centre->index) new_swarm.append(ampl_ids[c->raws[j]->index], 1, config.main_threshold + 1, c->centre->index, config.main_threshold + 1, c->raws[j]->abundance, true);
            }

        }

        // free memory
        delete clusters;
        for (numSeqs_t index = 0; index < ampl_ids.size(); index++) {
            delete raws[index];
        }
        delete[] raws;
        if (use_k_mers) {

            delete[] k8;
            delete[] k16;
            delete[] kord;

        }

        delete qe;

    }

    Dada2Clusterer::Cluster::Cluster(unsigned int num) {

        seq = new char[1];
        seq[0] = '\0';
        centre = nullptr;
        size = 0;
        mass = 0;
        i = num;
        raws = new RawSequence*[Dada2Utility::raw_buffer];
        max_raws = Dada2Utility::raw_buffer;

        update_e = true;
        check_locks = true;

    }

    Dada2Clusterer::Cluster::~Cluster() {

        delete[] seq;
        delete[] raws;

    }

    unsigned int Dada2Clusterer::Cluster::add_raw(RawSequence* raw) {

        // allocate more space if needed
        if (size >= max_raws) { // extend RawSequence* buffer

            RawSequence** tmp = new RawSequence*[max_raws + Dada2Utility::raw_buffer];
            std::memcpy(tmp, raws, max_raws * sizeof(RawSequence*));
            delete[] raws;
            raws = tmp;
            max_raws += Dada2Utility::raw_buffer;

        }

        // add raw and update mass / size
        raws[size] = raw;
        mass += raw->abundance;
        update_e = true;

        return size++;

    }

    void Dada2Clusterer::Cluster::census() {

        unsigned int new_mass = 0, new_nraw = 0;

        for (unsigned int r = 0; r < size; r++) {

            new_mass += raws[r]->abundance;
            new_nraw++;

        }

        if (new_mass != mass) {
            update_e = true;
        }

        mass = new_mass;
        size = new_nraw;

    }

    void Dada2Clusterer::Cluster::assign_centre() {

        unsigned int max_reads = 0;

        // assign the raw with the most reads as the centre
        centre = nullptr;
        for (unsigned int r = 0; r < size; r++) {

            raws[r]->lock = false; // unlock everything as centre is changing
            if (raws[r]->abundance > max_reads) {

                centre = raws[r];
                max_reads = centre->abundance;

            }

        }

        // assign centre sequence to seq and flag check_locks
        if (centre) {

            delete[] seq;
            seq = new char[centre->length + 1];
            strcpy(seq, centre->seq);

        }
        check_locks = true;

    }

    double Dada2Clusterer::Cluster::get_p_a(RawSequence* raw) {

        double pval;

        if (raw->abundance == 1 && !raw->prior) { // singleton -> no abundance p-value
            pval = 1.0;
        }
        else if (raw->comp.hamming == 0) { // cluster centre (or no mismatch to centre)
            pval = 1.0;
        }
        else if (raw->comp.lambda == 0) { // zero expected reads of this raw
            pval = 0.0;
        }
        else { // calculate abundance p-value

            // exp_reads is the expected number of reads for this raw
            double exp_reads = raw->comp.lambda * mass;
            pval = Dada2Utility::calc_p_a(raw->abundance, exp_reads, raw->prior);

        }

        return pval;

    }

    Dada2Clusterer::RawSequence* Dada2Clusterer::Cluster::pop_raw(unsigned int r) {

        RawSequence* pop;

        if (r < size) {

            pop = raws[r];
            raws[r] = raws[size - 1]; // popped raw is replaced by last raw
            raws[size - 1] = nullptr;
            size--;
            mass -= pop->abundance;
            update_e = true;

        } else {
            throw std::runtime_error("Container Error (Cluster): Tried to pop out-of-range raw.");
        }

        return pop;

    }

    Dada2Clusterer::ClusterCollection::ClusterCollection(RawSequence** raw_sequences, unsigned int raws_size, double om_a, double om_p) {

        // allocate memory
        clusters = new Cluster*[Dada2Utility::clust_buffer];
        max_clusters = Dada2Utility::clust_buffer;

        // initialise basic values
        num_clusters = 0;
        total_mass = 0;
        num_raws = raws_size;
        omega_a = om_a;
        omega_p = om_p;

        // initialise array of pointers to raws (allocated outside the construction of the cluster collection)
        raws = raw_sequences;
        unsigned int max_len = 0;
        for (unsigned int index = 0; index < num_raws; index++) {

            raws[index]->index = index;
            total_mass += raws[index]->abundance;
            max_len = std::max(max_len, raws[index]->length);

        }

        tvec = new int[max_len];
        qind = new unsigned int[max_len];

        // initialise with one cluster containing all the raws
        add_cluster();
        for (unsigned int index = 0; index < num_raws; index++) {
            clusters[0]->add_raw(raws[index]);
        }

        clusters[0]->census();
        clusters[0]->assign_centre(); // makes cluster centre sequence

    }

    Dada2Clusterer::ClusterCollection::~ClusterCollection() {

        delete[] qind;
        delete[] tvec;
        for (int i = 0; i < num_clusters; i++) {
            delete clusters[i];
        }
        delete[] clusters;

    }

    unsigned int Dada2Clusterer::ClusterCollection::add_cluster() {

        // allocate more space if needed
        if (num_clusters >= max_clusters) { // extend Cluster* buffer

            Cluster** tmp = new Cluster*[max_clusters + Dada2Utility::clust_buffer];
            std::memcpy(tmp, clusters, max_clusters * sizeof(Cluster*));
            delete[] clusters;
            clusters = tmp;
            max_clusters += Dada2Utility::clust_buffer;

        }

        // add cluster and update num_clusters
        clusters[num_clusters] = new Cluster(num_clusters);

        return num_clusters++;

    }

    void Dada2Clusterer::ClusterCollection::p_update(bool greedy) {

        for (unsigned int i = 0; i < num_clusters; i++) {

            Cluster* cluster = clusters[i];
            if (cluster->update_e) {

                for (unsigned int r = 0; r < cluster->size; r++) {

                    RawSequence* raw_sequence = cluster->raws[r];
                    raw_sequence->p = cluster->get_p_a(raw_sequence);

                }
                cluster->update_e = false;

            }

            if (greedy && cluster->check_locks) { // lock raw if its abundance is less than the expected abundance from just the centre

                for (unsigned int r = 0; r < cluster->size; r++) {

                    RawSequence* raw_sequence = cluster->raws[r];
                    double exp_reads_centre = clusters[i]->centre->abundance * raw_sequence->comp.lambda;
                    if (exp_reads_centre > raw_sequence->abundance) {
                        raw_sequence->lock = true;
                    }
                    if (raw_sequence == cluster->centre) {
                        raw_sequence->lock = true;
                    }

                }
                cluster->check_locks = false; // locking can only happen first time around

            }

        }

    }

    int Dada2Clusterer::ClusterCollection::bud(double min_fold, int min_hamming, int min_abund) {

        int min_cluster = -1, min_raw_index = -1, min_cluster_prior = -1, min_raw_index_prior = -1; // negative initialisation

        // assumption: complete alignment / p-value calculations were performed in the initial cluster
        RawSequence* min_raw = clusters[0]->centre;
        RawSequence* min_raw_prior = clusters[0]->centre;

        // find i, r indices and value of minimum p-value
        for (int i = 0; i < num_clusters; i++) {

            for (int r = 1; r < clusters[i]->size; r++) { // skip centre at r = 0

                RawSequence* raw = clusters[i]->raws[r];
                if (raw->abundance < min_abund) {
                    continue;
                }
                int hamming = raw->comp.hamming;
                double lambda = raw->comp.lambda;

                // calculate the fold over-abundance and the Hamming distance to this raw
                if (hamming >= min_hamming) { // only those passing the Hamming / fold screens can be budded

                    if (min_fold <= 1 || ((double)raw->abundance) >= min_fold * lambda * clusters[i]->mass) {

                        if ((raw->p < min_raw->p) || ((raw->p == min_raw->p && raw->abundance > min_raw->abundance))) { // most significant

                            min_cluster = i; min_raw_index = r;
                            min_raw = raw;

                        }
                        if (raw->prior && ((raw->p < min_raw_prior->p) || (raw->p == min_raw_prior->p && raw->abundance > min_raw_prior->abundance))) { // most significant

                            min_cluster_prior = i; min_raw_index_prior = r;
                            min_raw_prior = raw;

                        }

                    }

                }

            }

        }

        // Bonferroni correct the abundance p-value by the number of raws and compare to omega_a
        // (quite conservative, although probably unimportant given the abundance model issues)
        double p_a = min_raw->p * num_raws;
        double p_p = min_raw_prior->p;
        if (p_a < omega_a && min_cluster >= 0) {  // a significant abundance p-value

            RawSequence* raw = clusters[min_cluster]->pop_raw(min_raw_index);
            int i = add_cluster();

            // add raw to new cluster
            clusters[i]->add_raw(raw);
            clusters[i]->assign_centre();
            return i;

        } else if (p_p < omega_p && min_cluster_prior >= 0) { // a significant prior-abundance p-value

            RawSequence* raw = clusters[min_cluster_prior]->pop_raw(min_raw_index_prior);
            int i = add_cluster();

            // add raw to new cluster
            clusters[i]->add_raw(raw);
            clusters[i]->assign_centre();
            return i;

        }

        return 0;

    }

    bool Dada2Clusterer::ClusterCollection::shuffle() {

        double* emax = new double[num_raws]; // highest expected abundance for each raw
        Comparison** compmax = new Comparison*[num_raws]; // pointer to comparison leading to emax-value for each raw

        // initialise emax / compmax off of cluster 0
        // comparisons to all raws exist in cluster 0 (in index order)
        for (unsigned int index = 0; index < num_raws; index++) {

            compmax[index] = &clusters[0]->comps[index];
            emax[index] = compmax[index]->lambda * clusters[0]->mass;

        }

        // iterate over comparisons and find comparison with best expectation for each raw
        for (unsigned int i = 1; i < num_clusters; i++) {

            for (unsigned int cind = 0; cind < clusters[i]->comps.size(); cind++) {

                Comparison* comp = &clusters[i]->comps[cind];
                unsigned int index = comp->index;
                double e = comp->lambda * clusters[i]->mass;
                if (e > emax[index]) { // better expectation

                    compmax[index] = comp;
                    emax[index] = e;

                }

            }

        }

        bool shuffled = false;

        // iterate over raws: if best i different than current, then move
        for (unsigned int i = 0; i < num_clusters; i++) {

            // important: iterate backwards due to pop_raw(...)
            for (int r = clusters[i]->size - 1; r >= 0; r--) {

                RawSequence* raw_sequence = clusters[i]->raws[r];

                // if a better cluster was found, move the raw to the new cluster
                if (compmax[raw_sequence->index]->i != i) {

                    if (raw_sequence->index == clusters[i]->centre->index) {  // check if centre
                        continue;
                    }

                    // move raw and assign raw the comparison of its new cluster
                    clusters[i]->pop_raw(r);
                    clusters[compmax[raw_sequence->index]->i]->add_raw(raw_sequence);
                    raw_sequence->comp = *compmax[raw_sequence->index];
                    shuffled = true;

                }

            }

        }

        delete[] compmax;
        delete[] emax;

        return shuffled;

    }

    void Dada2Clusterer::ClusterCollection::compare(unsigned int i, std::vector<std::vector<double>>& error_matrix,
                                                 int match, int mismatch, int gap_pen, int homo_gap_pen, bool use_k_mers,
                                                 double k_dist_cutoff, int band_size, bool gapless, bool greedy) {

        // align all raws to this sequence and compute corresponding lambda
        unsigned int centre_reads = clusters[i]->centre->abundance;

        for (unsigned int index = 0; index < num_raws; index++) {

            RawSequence* raw_sequence = raws[index];
            Sub* sub;

            // get sub object
            if (greedy && (raw_sequence->abundance > centre_reads)) {
                sub = nullptr;
            } else if (greedy && raw_sequence->lock) { // was (raw->p > 0.5)
                sub = nullptr;
            } else {
                sub = new Sub(clusters[i]->centre, raw_sequence, match, mismatch, gap_pen, homo_gap_pen, use_k_mers,
                              k_dist_cutoff, band_size, gapless);
            }

            // calculate lambda for that sub
            double lambda = Dada2Utility::compute_lambda(raw_sequence->seq, raw_sequence->length, raw_sequence->qual, sub, error_matrix, tvec, qind);

            // store comparison if potentially useful
            if (lambda * total_mass > raw_sequence->max_exp) { // this cluster could attract this raw

                if (lambda * clusters[i]->centre->abundance > raw_sequence->max_exp) { // better max_exp
                    raw_sequence->max_exp = lambda * clusters[i]->centre->abundance;
                }

                Comparison comp;
                comp.i = i;
                comp.index = index;
                comp.lambda = lambda;
                comp.hamming = sub->nsubs;

                clusters[i]->comps.push_back(comp);
                if (i == 0 || raw_sequence == clusters[i]->centre) { // update on init (i = 0) or if the centre (as b_bud doesn't update raw->comp)
                    raw_sequence->comp = comp;
                }

            }

            delete sub;

        }


    }

    Dada2Clusterer::ClusterCollection* Dada2Clusterer::run_dada(RawSequence** raws, int num_raw,
            std::vector<std::vector<double>>& error_matrix, int match, int mismatch, int gap_pen, int homo_gap_pen,
            bool use_k_mers, double k_dist_cutoff, int band_size, double omega_a, double omega_p, int max_clust,
            double min_fold, int min_hamming, int min_abund, bool gapless, bool greedy) {

        auto bb = new ClusterCollection(raws, num_raw, omega_a, omega_p); // new cluster with all sequences in one cluster

        // everyone gets aligned within the initial cluster, no k-mer screen
        bb->compare(0, error_matrix, match, mismatch, gap_pen, homo_gap_pen, false, 1.0, band_size, gapless, greedy);
        bb->p_update(greedy); // calculates abundance p-value for each raw in its cluster (consensuses)

        if (max_clust < 1) {
            max_clust = bb->num_raws;
        }

        int newi = 0;
        int nshuffle = 0;
        bool shuffled;
        while ((bb->num_clusters < max_clust) && (newi = bb->bud(min_fold, min_hamming, min_abund))) {

            bb->compare(newi, error_matrix, match, mismatch, gap_pen, homo_gap_pen, use_k_mers, k_dist_cutoff,
                        band_size, gapless, greedy);

            // keep shuffling and updating until no more shuffles
            nshuffle = 0;
            do {
                shuffled = bb->shuffle();
            } while (shuffled && ++nshuffle < Dada2Utility::max_shuffle);
            if (nshuffle >= Dada2Utility::max_shuffle) {
                printf("Warning: Reached maximum (%i) shuffles.\n", Dada2Utility::max_shuffle);
            }

            bb->p_update(greedy);

        }

        return bb;

    }


    ConsistentClassicSwarmer::ConsistentClassicSwarmer(const Configuration& config) {

        auto& misc = config.misc;
        auto iter = misc.begin();

        if ((iter = misc.find("use_k_mers")) != misc.end()) {
            use_k_mers = iter->second == "true";
        }
        if ((iter = misc.find("k_dist_cutoff")) != misc.end()) {
            k_dist_cutoff = std::stod(iter->second);
        }
        if ((iter = misc.find("band_size")) != misc.end()) {
            band_size = std::stoi(iter->second);
        }
        if ((iter = misc.find("omega_a")) != misc.end()) {
            omega_a = std::stod(iter->second);
        }
        if ((iter = misc.find("min_fold")) != misc.end()) {
            min_fold = std::stod(iter->second);
        }
        if ((iter = misc.find("min_hamming")) != misc.end()) {
            min_hamming = std::stoi(iter->second);
        }
        if ((iter = misc.find("gapless")) != misc.end()) {
            gapless = iter->second == "true";
        }
        if ((iter = misc.find("error_matrix")) != misc.end()) {
            matrix_file = iter->second;
        }

    }

    ConsistentClassicSwarmer* ConsistentClassicSwarmer::clone() const {
        return new ConsistentClassicSwarmer(*this);
    }

    SwarmStorage* ConsistentClassicSwarmer::cluster(const AmpliconStorage& amplicon_storage, const Configuration& config) {

        std::cout << "Clustering amplicons into swarms..." << std::endl;

        SwarmStorage* swarm_storage = config.build_swarm_storage(amplicon_storage);

        // set up error matrix & temporary arrays
        lenSeqs_t alphabet_size = (config.alphabet.empty()) ? 4 : config.alphabet.size();
        QualityEncoding<>* qe = config.build_quality_encoding();
        char score_shift = qe->get_offset() + qe->get_from(); // subtracted from quality scores before accessing the error matrix
        std::vector<std::vector<double>> error_matrix;
        if (!Dada2Utility::fill_error_matrix(error_matrix, matrix_file, qe->get_accepted_scores().size())) {

            for (auto c1 = 0; c1 < alphabet_size; c1++) {
                for (auto c2 = 0; c2 < alphabet_size; c2++) {
                    error_matrix.emplace_back();
                    for (auto q : qe->get_accepted_scores()) {
                        auto err_prob = qe->encoded_quality_to_probability(q);
                        error_matrix[c1 * alphabet_size + c2].push_back((c1 == c2) ? (1 - err_prob) : (err_prob / (alphabet_size - 1)));
                    }
                }
            }

            std::cout << " -- Error matrix filled automatically based on the quality-encoding scheme." << std::endl;

        }
        int* tvec = new int[amplicon_storage.max_length()];
        unsigned int* qind = new unsigned int[amplicon_storage.max_length()];
        char* seq = new char[amplicon_storage.max_length() + 1];
        uint8_t* qual = new uint8_t[amplicon_storage.max_length()];

        size_t n_kmer = 1 << (2 * Dada2Utility::k_mer_size);
        auto cur_seed_kmer8 = new uint8_t[n_kmer];
        auto cur_seed_kmer = new uint16_t[n_kmer];
        auto cur_seed_kord = new uint16_t[amplicon_storage.max_length()];
        auto cur_partner_kmer8 = new uint8_t[n_kmer];
        auto cur_partner_kmer = new uint16_t[n_kmer];
        auto cur_partner_kord = new uint16_t[amplicon_storage.max_length()];

        for (numSeqs_t p = 0; p < amplicon_storage.num_pools(); p++) {

            auto& ac = amplicon_storage.get_pool(p);
            auto& swarms = swarm_storage->get_swarms(p);

            AuxiliaryData* aux = config.build_auxiliary_data(amplicon_storage, p, config.main_threshold);
            std::vector<bool> swarmed(ac.size(), false); // swarmed amplicons are already included in a cluster

            // open new swarm for the amplicon with the highest abundance that is not yet included in a swarm
            // (corresponds to next unswarmed amplicon due to sorting of amplicons in the pool)
            for (numSeqs_t seed = 0; seed < ac.size(); seed++) {

                if (!swarmed[seed]) {

                    /* (a) Initialise new swarm with seed */
                    auto& cur_swarm = swarms.initialise_cluster(seed, ac.ab(seed));
                    swarmed[seed] = true;
                    aux->record_amplicon(seed);
                    aux->tick_off_amplicon(seed);
                    SwarmEntry* cur_seed = cur_swarm.seed_entry();

                    lenSeqs_t last_gen = 0;

                    /* (b) Proceed breadth-first through the 'amplicon space' */
                    numSeqs_t pos = 0;
                    while (pos < cur_swarm.size()) { // expand current swarm until no further similar amplicons can be added

                        if (last_gen != cur_swarm.gen(pos)) { // sort upcoming generation to work through it by decreasing abundance

                            cur_swarm.sort_next_generation(ac, pos);
                            aux->clear_amplicon_records();

                        }

                        // get next swarm (sub)seed
                        cur_seed = cur_swarm.get_entry(pos, cur_seed);
                        numSeqs_t cur_seed_gen = cur_seed->gen();
                        dist_t cur_seed_rad = cur_seed->rad();
                        numSeqs_t cur_seed_mem = cur_seed->member();

                        // construct a raw for the current (sub)seed
                        Dada2Utility::nt2int(seq, ac.seq(cur_seed_mem));
                        auto quals = ac.quals(cur_seed_mem);
                        for (lenSeqs_t ppos = 0; ppos < ac.len(cur_seed_mem); ppos++) {
                            qual[ppos] = (uint8_t)(quals[ppos] - score_shift);
                        }
                        Dada2Utility::RawSequence cur_seed_raw(seq, ac.len(cur_seed_mem), qual, ac.ab(cur_seed_mem), false);
                        if (use_k_mers) {

                            cur_seed_raw.kmer8 = cur_seed_kmer8;
                            Dada2Utility::assign_kmer8(cur_seed_raw.kmer8, cur_seed_raw.seq, Dada2Utility::k_mer_size);
                            cur_seed_raw.kmer = cur_seed_kmer;
                            Dada2Utility::assign_kmer(cur_seed_raw.kmer, cur_seed_raw.seq, Dada2Utility::k_mer_size);
                            cur_seed_raw.kord = cur_seed_kord;
                            Dada2Utility::assign_kmer_order(cur_seed_raw.kord, cur_seed_raw.seq, Dada2Utility::k_mer_size);

                        }


                        // Consider yet unseen (unswarmed) amplicons to continue the exploration.
                        // An amplicon is marked as 'swarmed' as soon as it occurs the first time as a partner
                        // in order to prevent the algorithm from queueing it more than once coming from different amplicons.
                        for (auto& partner : aux->find_partners(cur_seed_mem)) {

                            Dada2Utility::nt2int(seq, ac.seq(partner.id));
                            quals = ac.quals(partner.id);
                            for (lenSeqs_t ppos = 0; ppos < ac.len(partner.id); ppos++) {
                                qual[ppos] = (uint8_t)(quals[ppos] - score_shift);
                            }
                            Dada2Utility::RawSequence cur_partner_raw(seq, ac.len(partner.id), qual, ac.ab(partner.id), false);
                            if (use_k_mers) {

                                cur_partner_raw.kmer8 = cur_partner_kmer8;
                                Dada2Utility::assign_kmer8(cur_partner_raw.kmer8, cur_partner_raw.seq, Dada2Utility::k_mer_size);
                                cur_partner_raw.kmer = cur_partner_kmer;
                                Dada2Utility::assign_kmer(cur_partner_raw.kmer, cur_partner_raw.seq, Dada2Utility::k_mer_size);
                                cur_partner_raw.kord = cur_partner_kord;
                                Dada2Utility::assign_kmer_order(cur_partner_raw.kord, cur_partner_raw.seq, Dada2Utility::k_mer_size);

                            }

                            Dada2Utility::Sub sub(&cur_seed_raw, &cur_partner_raw,
                                    config.match_reward, config.mismatch_penalty,
                                    config.gap_opening_penalty + config.gap_extension_penalty,
                                    config.gap_opening_penalty + config.gap_extension_penalty,
                                    use_k_mers, k_dist_cutoff, band_size, gapless);

                            // calculate lambda for that sub
                            double lambda = Dada2Utility::compute_lambda(cur_partner_raw.seq, cur_partner_raw.length,
                                    cur_partner_raw.qual, &sub, error_matrix, tvec, qind);

                            // calculate the fold over-abundance and the Hamming distance to this raw
                            if (sub.nsubs >= min_hamming && (min_fold <= 1 || ac.ab(partner.id) >= min_fold * lambda * ac.ab(cur_seed_mem))) { // only those passing the Hamming / fold screens can be budded

                                double p_value;
                                if (ac.ab(partner.id) == 1 || sub.nsubs == 0) {
                                    p_value = 1.0;
                                } else if (lambda == 0) {
                                    p_value = 0.0;
                                } else {

                                    // exp_reads is the expected number of reads for this raw
                                    double exp_reads = lambda * ac.ab(cur_seed_mem);
                                    p_value = Dada2Utility::calc_p_a(ac.ab(partner.id), exp_reads, false);

                                }

                                // Bonferroni correct the abundance p-value by the number of raws and compare to omega_a
                                // (quite conservative, although probably unimportant given the abundance model issues)
                                double p_a = p_value * amplicon_storage.num_amplicons();
                                if (p_a >= omega_a) {  // a non-significant abundance p-value -> would not bud (i.e. stay in cluster) in DADA2

                                    cur_swarm.append(
                                            partner.id, // pool-internal integer id of new member
                                            cur_seed_gen + 1, // generation number of new member
                                            cur_seed_rad + partner.dist, // radius of new member
                                            cur_seed_mem, // pool-internal integer id of parent
                                            partner.dist, // distance to parent
                                            ac.ab(partner.id), // abundance of new member
                                            (partner.dist != 0) && aux->record_amplicon(partner.id) // flag indicating whether new member is also a new sequence
                                            );

                                    swarmed[partner.id] = true;
                                    aux->tick_off_amplicon(partner.id);

                                }

                            }

                        }

                        last_gen = cur_seed_gen;
                        pos++;

                    }

                    /* (c) Close the no longer extendable swarm */
                    aux->clear_amplicon_records();
                    delete cur_seed;

                }

            }

            delete aux;

        }

        delete[] cur_partner_kord;
        delete[] cur_partner_kmer;
        delete[] cur_partner_kmer8;
        delete[] cur_seed_kord;
        delete[] cur_seed_kmer;
        delete[] cur_seed_kmer8;

        delete[] qual;
        delete[] seq;
        delete[] qind;
        delete[] tvec;

        delete qe;

        return swarm_storage;

    }


    ConsistencySwarmer::ConsistencySwarmer(const Configuration& config) {

        auto& misc = config.misc;
        auto iter = misc.begin();

        if ((iter = misc.find("use_k_mers")) != misc.end()) {
            use_k_mers = iter->second == "true";
        }
        if ((iter = misc.find("k_dist_cutoff")) != misc.end()) {
            k_dist_cutoff = std::stod(iter->second);
        }
        if ((iter = misc.find("band_size")) != misc.end()) {
            band_size = std::stoi(iter->second);
        }
        if ((iter = misc.find("omega_a")) != misc.end()) {
            omega_a = std::stod(iter->second);
        }
        if ((iter = misc.find("min_fold")) != misc.end()) {
            min_fold = std::stod(iter->second);
        }
        if ((iter = misc.find("min_hamming")) != misc.end()) {
            min_hamming = std::stoi(iter->second);
        }
        if ((iter = misc.find("gapless")) != misc.end()) {
            gapless = iter->second == "true";
        }
        if ((iter = misc.find("error_matrix")) != misc.end()) {
            matrix_file = iter->second;
        }

    }


    ConsistencySwarmer* ConsistencySwarmer::clone() const {
        return new ConsistencySwarmer(*this);
    }

    SwarmStorage* ConsistencySwarmer::cluster(const AmpliconStorage& amplicon_storage, const Configuration& config) {

        std::cout << "Clustering amplicons into swarms..." << std::endl;

        SwarmStorage* swarm_storage = config.build_swarm_storage(amplicon_storage);

        // set up error matrix & temporary arrays
        lenSeqs_t alphabet_size = (config.alphabet.empty()) ? 4 : config.alphabet.size();
        QualityEncoding<>* qe = config.build_quality_encoding();
        char score_shift = qe->get_offset() + qe->get_from(); // subtracted from quality scores before accessing the error matrix
        std::vector<std::vector<double>> error_matrix;
        if (!Dada2Utility::fill_error_matrix(error_matrix, matrix_file, qe->get_accepted_scores().size())) {

            for (auto c1 = 0; c1 < alphabet_size; c1++) {
                for (auto c2 = 0; c2 < alphabet_size; c2++) {
                    error_matrix.emplace_back();
                    for (auto q : qe->get_accepted_scores()) {
                        auto err_prob = qe->encoded_quality_to_probability(q);
                        error_matrix[c1 * alphabet_size + c2].push_back((c1 == c2) ? (1 - err_prob) : (err_prob / (alphabet_size - 1)));
                    }
                }
            }

            std::cout << " -- Error matrix filled automatically based on the quality-encoding scheme." << std::endl;

        }
        int* tvec = new int[amplicon_storage.max_length()];
        unsigned int* qind = new unsigned int[amplicon_storage.max_length()];
        char* seq = new char[amplicon_storage.max_length() + 1];
        uint8_t* qual = new uint8_t[amplicon_storage.max_length()];

        size_t n_kmer = 1 << (2 * Dada2Utility::k_mer_size);
        auto cur_seed_kmer8 = new uint8_t[n_kmer];
        auto cur_seed_kmer = new uint16_t[n_kmer];
        auto cur_seed_kord = new uint16_t[amplicon_storage.max_length()];
        auto cur_partner_kmer8 = new uint8_t[n_kmer];
        auto cur_partner_kmer = new uint16_t[n_kmer];
        auto cur_partner_kord = new uint16_t[amplicon_storage.max_length()];

        for (numSeqs_t p = 0; p < amplicon_storage.num_pools(); p++) {

            auto& ac = amplicon_storage.get_pool(p);
            auto& swarms = swarm_storage->get_swarms(p);

            std::vector<bool> swarmed(ac.size(), false); // swarmed amplicons are already included in a cluster
            std::set<numSeqs_t> unswarmed; // pool-internal integer ids of non-singleton amplicons not yet included in a swarm
            for (numSeqs_t i = 0; i < ac.size(); i++) {
                if (ac.ab(i) > 1) unswarmed.insert(i);
            }

            // open new swarm for the amplicon with the highest abundance that is not yet included in a swarm
            // (corresponds to next unswarmed amplicon due to sorting of amplicons in the pool)
            for (numSeqs_t seed = 0; seed < ac.size(); seed++) {

                if (!swarmed[seed]) {

                    /* (a) Initialise new swarm with seed */
                    auto& cur_swarm = swarms.initialise_cluster(seed, ac.ab(seed));
                    swarmed[seed] = true;
                    auto iii = unswarmed.erase(seed);
                    SwarmEntry* cur_seed = cur_swarm.seed_entry();

                    lenSeqs_t last_gen = 0;

                    /* (b) Proceed breadth-first through the 'amplicon space' */
                    numSeqs_t pos = 0;
                    while (pos < cur_swarm.size()) { // expand current swarm until no further similar amplicons can be added

                        if (last_gen != cur_swarm.gen(pos)) { // sort upcoming generation to work through it by decreasing abundance
                            cur_swarm.sort_next_generation(ac, pos);
                        }

                        // get next swarm (sub)seed
                        cur_seed = cur_swarm.get_entry(pos, cur_seed);
                        numSeqs_t cur_seed_gen = cur_seed->gen();
                        dist_t cur_seed_rad = cur_seed->rad();
                        numSeqs_t cur_seed_mem = cur_seed->member();

                        // construct a raw for the current (sub)seed
                        Dada2Utility::nt2int(seq, ac.seq(cur_seed_mem));
                        auto quals = ac.quals(cur_seed_mem);
                        for (lenSeqs_t ppos = 0; ppos < ac.len(cur_seed_mem); ppos++) {
                            qual[ppos] = (uint8_t)(quals[ppos] - score_shift);
                        }
                        Dada2Utility::RawSequence cur_seed_raw(seq, ac.len(cur_seed_mem), qual, ac.ab(cur_seed_mem), false);
                        if (use_k_mers) {

                            cur_seed_raw.kmer8 = cur_seed_kmer8;
                            Dada2Utility::assign_kmer8(cur_seed_raw.kmer8, cur_seed_raw.seq, Dada2Utility::k_mer_size);
                            cur_seed_raw.kmer = cur_seed_kmer;
                            Dada2Utility::assign_kmer(cur_seed_raw.kmer, cur_seed_raw.seq, Dada2Utility::k_mer_size);
                            cur_seed_raw.kord = cur_seed_kord;
                            Dada2Utility::assign_kmer_order(cur_seed_raw.kord, cur_seed_raw.seq, Dada2Utility::k_mer_size);

                        }


                        // Consider yet unseen (unswarmed) amplicons to continue the exploration.
                        // An amplicon is marked as 'swarmed' as soon as it occurs the first time as a partner
                        // in order to prevent the algorithm from queueing it more than once coming from different amplicons.
                        for (auto iter = unswarmed.begin(); iter != unswarmed.end(); ) {

                            numSeqs_t partner_id = *iter;

                            if (ac.ab(partner_id) <= ac.ab(cur_seed_mem)) {

                                Dada2Utility::nt2int(seq, ac.seq(partner_id));
                                quals = ac.quals(partner_id);
                                for (lenSeqs_t ppos = 0; ppos < ac.len(partner_id); ppos++) {
                                    qual[ppos] = (uint8_t)(quals[ppos] - score_shift);
                                }
                                Dada2Utility::RawSequence cur_partner_raw(seq, ac.len(partner_id), qual, ac.ab(partner_id), false);
                                if (use_k_mers) {

                                    cur_partner_raw.kmer8 = cur_partner_kmer8;
                                    Dada2Utility::assign_kmer8(cur_partner_raw.kmer8, cur_partner_raw.seq, Dada2Utility::k_mer_size);
                                    cur_partner_raw.kmer = cur_partner_kmer;
                                    Dada2Utility::assign_kmer(cur_partner_raw.kmer, cur_partner_raw.seq, Dada2Utility::k_mer_size);
                                    cur_partner_raw.kord = cur_partner_kord;
                                    Dada2Utility::assign_kmer_order(cur_partner_raw.kord, cur_partner_raw.seq, Dada2Utility::k_mer_size);

                                }

                                Dada2Utility::Sub sub(&cur_seed_raw, &cur_partner_raw,
                                        config.match_reward, config.mismatch_penalty, config.gap_opening_penalty + config.gap_extension_penalty, config.gap_opening_penalty + config.gap_extension_penalty,
                                        use_k_mers, k_dist_cutoff, band_size, gapless);

                                // calculate lambda for that sub
                                double lambda = Dada2Utility::compute_lambda(cur_partner_raw.seq, cur_partner_raw.length, cur_partner_raw.qual, &sub, error_matrix, tvec, qind);

                                // calculate the fold over-abundance and the Hamming distance to this raw
                                if (sub.nsubs >= min_hamming && (min_fold <= 1 || ac.ab(partner_id) >= min_fold * lambda * ac.ab(cur_seed_mem))) { // only those passing the Hamming / fold screens can be budded

                                    double p_value;
                                    if (ac.ab(partner_id) == 1 || sub.nsubs == 0) {
                                        p_value = 1.0;
                                    } else if (lambda == 0) {
                                        p_value = 0.0;
                                    } else {

                                        // exp_reads is the expected number of reads for this raw
                                        double exp_reads = lambda * ac.ab(cur_seed_mem);
                                        p_value = Dada2Utility::calc_p_a(ac.ab(partner_id), exp_reads, false);

                                    }

                                    // Bonferroni correct the abundance p-value by the number of raws and compare to omega_a
                                    // (quite conservative, although probably unimportant given the abundance model issues)
                                    double p_a = p_value * amplicon_storage.num_amplicons();
                                    if (p_a >= omega_a) {  // a non-significant abundance p-value -> would not bud (i.e. stay in cluster) in DADA2

                                        cur_swarm.append(
                                                partner_id, // pool-internal integer id of new member
                                                cur_seed_gen + 1, // generation number of new member
                                                cur_seed_rad + sub.nsubs, // radius of new member
                                                cur_seed_mem, // pool-internal integer id of parent
                                                sub.nsubs, // distance to parent
                                                ac.ab(partner_id), // abundance of new member
                                                true // flag indicating whether new member is also a new sequence
                                                );

                                        swarmed[partner_id] = true;

                                    }

                                }

                            }

                            if (swarmed[partner_id]) {
                                iter = unswarmed.erase(iter);
                            } else {
                                iter++;
                            }

                        }

                        last_gen = cur_seed_gen;
                        pos++;

                    }

                    /* (c) Close the no longer extendable swarm */
                    delete cur_seed;

                }

            }

        }

        delete[] cur_partner_kord;
        delete[] cur_partner_kmer;
        delete[] cur_partner_kmer8;
        delete[] cur_seed_kord;
        delete[] cur_seed_kmer;
        delete[] cur_seed_kmer8;

        delete[] qual;
        delete[] seq;
        delete[] qind;
        delete[] tvec;

        delete qe;

        return swarm_storage;

    }

}