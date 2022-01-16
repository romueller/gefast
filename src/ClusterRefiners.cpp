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

#include "../include/Clusterers.hpp"
#include "../include/ClusterRefiners.hpp"

namespace GeFaST {

    /* === GraftCandidate === */

    GraftCandidate::GraftCandidate() {

        parent_swarm = child_swarm = nullptr;
        parent_id = child_id = 0;
        distance = 0;

    }

    GraftCandidate::GraftCandidate(Swarm* ps, numSeqs_t p, Swarm* cs, numSeqs_t c, dist_t d) {

        parent_swarm = ps;
        parent_id = p;
        child_swarm = cs;
        child_id = c;
        distance = d;

    }


    /* === CompareGraftCandidatesAbund === */

    CompareGraftCandidatesAbund::CompareGraftCandidatesAbund(const AmpliconCollection& a) : ac(a) {
        // nothing else to do
    }

    bool CompareGraftCandidatesAbund::compare(const numSeqs_t ampl_a, const numSeqs_t ampl_b) {
        return (ac.ab(ampl_a) > ac.ab(ampl_b))
            || ((ac.ab(ampl_a) == ac.ab(ampl_b)) && (*ac.id(ampl_a) < *ac.id(ampl_b)));
    }

    bool CompareGraftCandidatesAbund::operator()(const GraftCandidate& a, const GraftCandidate& b) {
        return compare(a.parent_swarm->member(a.parent_id), b.parent_swarm->member(b.parent_id))
            || ((*ac.id(a.parent_swarm->member(a.parent_id)) == *ac.id(b.parent_swarm->member(b.parent_id)))
                && compare(a.child_swarm->member(a.child_id), b.child_swarm->member(b.child_id)));
    }


    /* === FastidiousRefiner === */

    FastidiousRefiner* FastidiousRefiner::clone() const {
        return new FastidiousRefiner(*this);
    }

    void FastidiousRefiner::refine(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage,
            const Configuration& config) {

        numSeqs_t out_cnt_s = 0;
        numSeqs_t out_cnt_a = 0;
        for (auto p = 0; p < swarm_storage.num_pools(); p++) {
            auto& out_sws = swarm_storage.get_swarms(p);
            for (auto s = 0; s < out_sws.size(); s++) {
                if (out_sws.get(s).mass() < config.boundary) out_cnt_s++;
                if (out_sws.get(s).mass() < config.boundary) out_cnt_a += out_sws.get(s).size();
            }
        }
        std::cout << " -- Light swarms before refinement: " << out_cnt_s << " (of " << swarm_storage.num_swarms() << ")" << std::endl;
        std::cout << " -- Amplicons in light swarms before refinement: " << out_cnt_a << " (of " << amplicon_storage.num_amplicons() << ")" << std::endl;

        // refine swarms in a pool when there are both light and heavy swarms
        for (numSeqs_t p = 0; p < amplicon_storage.num_pools(); p++) {

            bool light = false;
            bool heavy = false;
            auto& swarms = swarm_storage.get_swarms(p);

            for (numSeqs_t s = 0; !(light && heavy) && s < swarms.size(); s++) {

                light |= (swarms.get(s).mass() < config.boundary);
                heavy |= (swarms.get(s).mass() >= config.boundary);

            }

            if (light && heavy) {
                refine_pool(amplicon_storage, swarm_storage, p, out_cnt_s, out_cnt_a, config);
            }

        }

        std::cout << " -- Light swarms after refinement: " << out_cnt_s << " (of " << swarm_storage.num_swarms() << ")" << std::endl;
        std::cout << " -- Amplicons in light swarms after refinement: " << out_cnt_a << " (of " << amplicon_storage.num_amplicons() << ")" << std::endl;

    }

    bool FastidiousRefiner::compare_candidates(const AmpliconCollection& ac,
            const numSeqs_t ampl_a, const numSeqs_t ampl_b) {

        return (ac.ab(ampl_a) > ac.ab(ampl_b))
            || ((ac.ab(ampl_a) == ac.ab(ampl_b)) && (*ac.id(ampl_a) < *ac.id(ampl_b)));

    }

    void FastidiousRefiner::refine_pool(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage,
            const numSeqs_t pool_id, numSeqs_t& light_s, numSeqs_t& light_a, const Configuration& config) {

        auto& ac = amplicon_storage.get_pool(pool_id);
        auto& swarms = swarm_storage.get_swarms(pool_id);

        std::vector<GraftCandidate> graft_cands(ac.size()); // initially, graft candidates for all amplicons of the pool are "empty"
        for (numSeqs_t s = 0; s < swarms.size(); s++) {

            Swarm& swarm = swarms.get(s);
            for (numSeqs_t i = 0; i < swarm.size(); i++) {

                graft_cands[swarm.member(i)].child_swarm = &swarm;
                graft_cands[swarm.member(i)].child_id = i;

            }

        }

        AuxiliaryData* aux = config.build_auxiliary_data(amplicon_storage, swarm_storage, pool_id, config.refinement_threshold);

        // iterate over the (heavy) swarms and determine the graft candidates
        for (numSeqs_t s = 0; s < swarms.size(); s++) {

            Swarm& swarm = swarms.get(s);

            if (swarm.mass() >= config.boundary) {

                for (numSeqs_t i = 0; i < swarm.size(); i++) {

                    numSeqs_t ampl_id = swarm.member(i);
                    auto partners = aux->find_partners(ampl_id);

                    for (auto& p : partners) {

                        // add new graft candidate or replace the old candidate when the new one has a higher priority
                        if (((graft_cands[p.id].parent_swarm == nullptr)
                            || compare_candidates(ac, ampl_id, graft_cands[p.id].parent_swarm->member(graft_cands[p.id].parent_id)))) {

                            graft_cands[p.id].parent_swarm = &swarm;
                            graft_cands[p.id].parent_id = i;
                            graft_cands[p.id].distance = p.dist;

                        }

                    }
                }

            }

        }

        delete aux;


        // sort all actual grafting candidates and perform graftings
        auto new_end = std::remove_if(
                graft_cands.begin(),
                graft_cands.end(),
                [](GraftCandidate& gc) {
                    return gc.parent_swarm == nullptr;
                });

        std::sort(graft_cands.begin(), new_end, CompareGraftCandidatesAbund(ac));

        for (auto graft_iter = graft_cands.begin(); graft_iter != new_end; graft_iter++) {

            if (!(graft_iter->child_swarm->is_attached())) {

                Swarm* parent_swarm = graft_iter->parent_swarm;
                Swarm* child_swarm = graft_iter->child_swarm;

                // "attach the seed of the light swarm to the tail of the heavy swarm"
                // cluster entries are moved unchanged (i.e. entry of seed of attached swarm
                // stays 'incomplete', grafting 'link' is not recorded)
                parent_swarm->attach(graft_iter->parent_id, graft_iter->child_id, child_swarm, graft_iter->distance);

                light_s--;
                light_a -= child_swarm->size();

            }

        }

    }


    /* === IterativeFastidiousRefiner === */

    IterativeFastidiousRefiner* IterativeFastidiousRefiner::clone() const {
        return new IterativeFastidiousRefiner(*this);
    }

    void IterativeFastidiousRefiner::refine(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage,
            const Configuration& config) {

        numSeqs_t out_cnt_s = 0;
        numSeqs_t out_cnt_a = 0;
        for (auto p = 0; p < swarm_storage.num_pools(); p++) {
            auto& out_sws = swarm_storage.get_swarms(p);
            for (auto s = 0; s < out_sws.size(); s++) {
                if (out_sws.get(s).mass() < config.boundary) out_cnt_s++;
                if (out_sws.get(s).mass() < config.boundary) out_cnt_a += out_sws.get(s).size();
            }
        }
        std::cout << " -- Light swarms before refinement: " << out_cnt_s << " (of " << swarm_storage.num_swarms() << ")" << std::endl;
        std::cout << " -- Amplicons in light swarms before refinement: " << out_cnt_a << " (of " << amplicon_storage.num_amplicons() << ")" << std::endl;

        for (auto t : config.iterative_refinement_thresholds) {

            // refine swarms in a pool when there are both light and heavy swarms
            for (numSeqs_t p = 0; p < amplicon_storage.num_pools(); p++) {

                bool light = false;
                bool heavy = false;
                auto& swarms = swarm_storage.get_swarms(p);

                for (numSeqs_t s = 0; !(light && heavy) && s < swarms.size(); s++) {

                    light |= (swarms.get(s).mass() < config.boundary && !swarms.get(s).is_attached());
                    heavy |= (swarms.get(s).mass() >= config.boundary);

                }

                if (light && heavy) {
                    refine_pool(amplicon_storage, swarm_storage, p, t, out_cnt_s, out_cnt_a, config);
                }

            }

            std::cout << " -- Light swarms after refinement with threshold " << t << " : "
                      << out_cnt_s << " (of " << swarm_storage.num_swarms() << ")" << std::endl;
            std::cout << " -- Amplicons in light swarms after refinement with threshold " << t << ": "
                      << out_cnt_a << " (of " << amplicon_storage.num_amplicons() << ")" << std::endl;

        }

    }

    bool IterativeFastidiousRefiner::compare_candidates(const AmpliconCollection& ac,
            const numSeqs_t ampl_a, const numSeqs_t ampl_b) {

        return (ac.ab(ampl_a) > ac.ab(ampl_b))
            || ((ac.ab(ampl_a) == ac.ab(ampl_b)) && (*ac.id(ampl_a) < *ac.id(ampl_b)));

    }

    void IterativeFastidiousRefiner::refine_pool(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage,
            const numSeqs_t pool_id, const dist_t threshold, numSeqs_t& light_s, numSeqs_t& light_a, const Configuration& config) {

        auto& ac = amplicon_storage.get_pool(pool_id);
        auto& swarms = swarm_storage.get_swarms(pool_id);

        std::vector<GraftCandidate> graft_cands(ac.size()); // initially, graft candidates for all amplicons of the pool are "empty"
        for (numSeqs_t s = 0; s < swarms.size(); s++) {

            Swarm& swarm = swarms.get(s);
            for (numSeqs_t i = 0; i < swarm.size(); i++) {

                graft_cands[swarm.member(i)].child_swarm = &swarm;
                graft_cands[swarm.member(i)].child_id = i;

            }

        }

        AuxiliaryData* aux = config.build_auxiliary_data(amplicon_storage, swarm_storage, pool_id, threshold);

        // iterate over the (heavy) swarms and determine the graft candidates
        for (numSeqs_t s = 0; s < swarms.size(); s++) {

            Swarm& swarm = swarms.get(s);

            if (swarm.mass() >= config.boundary) {

                for (numSeqs_t i = 0; i < swarm.size(); i++) {

                    numSeqs_t ampl_id = swarm.member(i);
                    auto partners = aux->find_partners(ampl_id);

                    for (auto& p : partners) {

                        // add new graft candidate or replace the old candidate when the new one has a higher priority
                        if (((graft_cands[p.id].parent_swarm == nullptr)
                            || compare_candidates(ac, ampl_id, graft_cands[p.id].parent_swarm->member(graft_cands[p.id].parent_id)))) {

                            graft_cands[p.id].parent_swarm = &swarm;
                            graft_cands[p.id].parent_id = i;
                            graft_cands[p.id].distance = p.dist;

                        }

                    }
                }

            }

        }

        delete aux;


        // Sort all actual grafting candidates and perform graftings
        auto new_end = std::remove_if(
                graft_cands.begin(),
                graft_cands.end(),
                [](GraftCandidate& gc) {
                    return gc.parent_swarm == nullptr;
                });

        std::sort(graft_cands.begin(), new_end, CompareGraftCandidatesAbund(ac));

        for (auto graft_iter = graft_cands.begin(); graft_iter != new_end; graft_iter++) {

            if (!(graft_iter->child_swarm->is_attached())) {

                Swarm* parent_swarm = graft_iter->parent_swarm;
                Swarm* child_swarm = graft_iter->child_swarm;

                // "attach the seed of the light swarm to the tail of the heavy swarm"
                // cluster entries are moved unchanged (i.e. entry of seed of attached swarm
                // stays 'incomplete', grafting 'link' is not recorded)
                parent_swarm->attach(graft_iter->parent_id, graft_iter->child_id, child_swarm, graft_iter->distance);

                light_s--;
                light_a -= child_swarm->size();

            }

        }

    }


    /* === IdleRefiner === */

    IdleRefiner* IdleRefiner::clone() const {
        return new IdleRefiner(*this);
    }

    void IdleRefiner::refine(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage, const Configuration& config) {
        // do nothing
    }


    /* === LightSwarmAppender === */

    LightSwarmAppender::LightSwarmAppender(const Configuration& config) {

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
        if ((iter = misc.find("light_opt")) != misc.end()) {
            light_opt = std::stoi(iter->second);
        }

    }

    LightSwarmAppender* LightSwarmAppender::clone() const {
        return new LightSwarmAppender(*this);
    }

    void LightSwarmAppender::refine(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage, const Configuration& config) {

        char tmp_seq[amplicon_storage.max_length()];
        uint8_t tmp_qual[amplicon_storage.max_length()];
        int tvec[amplicon_storage.max_length()];
        unsigned int qind[amplicon_storage.max_length()];
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

            std::cout << " -- Error matrix filled automatically based on the quality-encoding scheme..." << std::endl;

        }

        // first use: number of light swarms / amplicons before refinement
        numSeqs_t out_cnt_s = 0;
        numSeqs_t out_cnt_a = 0;
        for (auto p = 0; p < swarm_storage.num_pools(); p++) {
            auto& out_sws = swarm_storage.get_swarms(p);
            for (auto s = 0; s < out_sws.size(); s++) {
                if (out_sws.get(s).mass() < config.boundary) out_cnt_s++;
                if (out_sws.get(s).mass() < config.boundary) out_cnt_a += out_sws.get(s).size();
            }
        }

        std::cout << " -- Number of light swarms before refinement: " << out_cnt_s << " (of " << swarm_storage.num_swarms() << ")" << std::endl;
        std::cout << " -- Number of amplicons in light swarms before refinement: " << out_cnt_a << " (of " << amplicon_storage.num_amplicons() << ")" << std::endl;

        // second use: number of ungrafted light swarms / amplicons
        out_cnt_s = 0;
        out_cnt_a = 0;

        for (numSeqs_t p = 0; p < amplicon_storage.num_pools(); p++) {

            auto& ac = amplicon_storage.get_pool(p);
            auto& swarms = swarm_storage.get_swarms(p);

            // build target clusters from heavy swarms (for access to mass / size, properties of seed)
            std::vector<TargetCluster> target_clusters;
            std::vector<numSeqs_t> light_swarms;
            numSeqs_t total_abund = 0;
            numSeqs_t total_num = ac.size();
            for (numSeqs_t s = 0; s < swarms.size(); s++) {

                auto& swarm = swarms.get(s);

                if (swarm.mass() >= config.boundary) {
                    target_clusters.emplace_back(s, ac, swarm.seed(), use_k_mers, *qe);
                } else {
                    light_swarms.emplace_back(s);
                }

                total_abund += swarm.mass();

            }

            std::vector<numSeqs_t> unattached_lights; // swarm ids of light swarms that were not grafted to a heavy swarm

            // for each light swarm:
            // - compare light seed (using the mass of the light swarm as the abundance) with centres of target clusters
            // - move whole light swarm to target cluster with highest E_minmax if not significant
            for (auto src_sid : light_swarms) {

                auto& src_swarm = swarms.get(src_sid);

                std::string src_seq = ac.seq_str(src_swarm.seed());
                Dada2Utility::nt2int(tmp_seq, src_seq.c_str());
                auto quals = ac.quals(src_swarm.seed());
                for (lenSeqs_t pos = 0; pos < ac.len(src_swarm.seed()); pos++) {
                    tmp_qual[pos] = (uint8_t)(quals[pos] - score_shift);
                }
                RawSequence src_seed(tmp_seq, ac.len(src_swarm.seed()), tmp_qual, ac.ab(src_swarm.seed()), false);
                src_seed.abundance = src_swarm.mass(); // use mass instead of abundance of seed as whole swarm is considered as a unit

                double max_exp = 0.0; // minimum expected abundance of the light swarm's seed, maximised over all target clusters (E_minmax in DADA2)
                numSeqs_t max_sid = target_clusters.size(); // swarm id of the target cluster corresponding to max_exp

                for (auto& tc : target_clusters) {

                    auto& target_swarm = swarms.get(tc.sid);
                    numSeqs_t target_seed_abund = ac.ab(target_swarm.seed());

                    Dada2Utility::Sub* sub;
                    if (greedy && (src_seed.abundance > target_seed_abund)) {
                        sub = nullptr;
                    } else {
                        sub = new Sub(&tc.seed, &src_seed, config.match_reward, config.mismatch_penalty,
                                      config.gap_opening_penalty + config.gap_extension_penalty, homo_gap, use_k_mers,
                                      k_dist_cutoff, band_size, gapless);
                    }

                    // calculate lambda for that sub
                    double lambda = Dada2Utility::compute_lambda(src_seed.seq, src_seed.length, src_seed.qual, sub, error_matrix, tvec, qind);

                    // store comparison if potentially useful
                    if (lambda * total_abund > max_exp) { // cluster could attract this raw //TODO? total_abund = mass of target cluster OR mass of all clusters (as in DADA2)?

                        if (lambda * target_seed_abund > max_exp) { // better E_minmax -> set

                            max_exp = lambda * target_seed_abund;
                            max_sid = tc.sid;

                            src_seed.comp.lambda = lambda;
                            src_seed.comp.hamming = sub->nsubs;

                            src_seed.p = get_p_a(target_swarm.mass(), src_seed);

                        }

                    }

                    delete sub;

                }

                if (max_sid < target_clusters.size() && !is_significant(src_seed, swarms.get(max_sid).mass(), total_num)) {
                    swarms.get(max_sid).attach(0, 0, &src_swarm, config.main_threshold + 1);
                } else {
                    unattached_lights.emplace_back(src_sid);
                    out_cnt_a += src_swarm.size();
                    out_cnt_s++;
                }

            }

            // processing of unattached light swarms
            std::vector<numSeqs_t> ungrafted_amplicons;
            switch (light_opt) {

                case 1:
                    // keep remaining light swarms as separate clusters
                    // nothing to do
                    break;

                case 2:
                    // discard the remaining light swarms completely
                    for (auto sid : unattached_lights) {
                        swarms.get(sid).mark_as_attached(); // "lazy delete" (will be ignored during output generation)
                    }
                    break;

                case 3:
                    // collect members of remaining light swarms in a single "unstructured" cluster
                    for (auto sid : unattached_lights) {
                        auto& swarm = swarms.get(sid);
                        for (auto i = 0; i < swarms.get(sid).size(); i++) {
                            ungrafted_amplicons.emplace_back(swarm.member(i));
                        }
                        swarm.mark_as_attached(); // "lazy delete" (will be ignored during output generation)
                    }
                    if (!ungrafted_amplicons.empty()) {
                        Swarm& new_swarm = swarms.initialise_cluster(ungrafted_amplicons[0], ac.ab(ungrafted_amplicons[0]));
                        for (auto i = 1; i < ungrafted_amplicons.size(); i++) {
                            new_swarm.append(ungrafted_amplicons[i], 1, config.main_threshold + 1, ungrafted_amplicons[0], config.main_threshold + 1, ac.ab(ungrafted_amplicons[i]), true);
                        }
                    }
                    break;

                case 4:
                    // collect members of remaining light swarms and perform a separate DADA2-like clustering on them
                    for (auto sid : unattached_lights) {
                        auto& swarm = swarms.get(sid);
                        for (auto i = 0; i < swarms.get(sid).size(); i++) {
                            ungrafted_amplicons.emplace_back(swarm.member(i));
                        }
                        swarm.mark_as_attached(); // "lazy delete" (will be ignored during output generation)
                    }
                    if (!ungrafted_amplicons.empty()) {
                        Dada2Clusterer::add_swarms(ac, ungrafted_amplicons, swarms, config, use_k_mers, k_dist_cutoff, band_size, omega_a,
                                                 omega_p, max_clust, min_fold, min_hamming, min_abund, homo_gap, gapless, greedy, matrix_file);
                    }
                    break;

                default:
                    std::cerr << " -- WARNING: Invalid light_opt option. Remaining light swarms are kept as separate clusters." << std::endl;

            }

        }

        std::cout << " -- Ungrafted light swarms: " << out_cnt_s << std::endl;
        std::cout << " -- Amplicons in ungrafted light swarms: " << out_cnt_a << std::endl;
        std::cout << " -- Processing option for ungrafted light swarms: " << light_opt << std::endl;

        // third use: number of light swarms / amplicons after refinement
        out_cnt_s = 0;
        out_cnt_a = 0;
        for (auto p = 0; p < swarm_storage.num_pools(); p++) {
            auto& out_sws = swarm_storage.get_swarms(p);
            for (auto s = 0; s < out_sws.size(); s++) {
                auto& sw = out_sws.get(s);
                if (sw.mass() < config.boundary && !sw.is_attached()) out_cnt_s++;
                if (sw.mass() < config.boundary && !sw.is_attached()) out_cnt_a += sw.size();
            }
        }

        std::cout << " -- Light swarms after refinement: " << out_cnt_s << " (of " << swarm_storage.num_swarms() << ")" << std::endl;
        std::cout << " -- Amplicons in light swarms after refinement: " << out_cnt_a << " (of " << amplicon_storage.num_amplicons() << ")" << std::endl;

        delete qe;

    }


    LightSwarmAppender::TargetCluster::TargetCluster(numSeqs_t s, const AmpliconCollection& ac, numSeqs_t a, bool use_k_mers, QualityEncoding<>& qe) : sid(s), seed() {

        seed.length = ac.len(a);
        seed.abundance = ac.ab(a);

        seed.seq = nullptr;
        seed.seq = new char[seed.length + 1];
        auto seq = ac.seq_str(a);
        Dada2Utility::nt2int(seed.seq, seq.c_str());

        seed.qual = new uint8_t[seed.length];
        auto q = ac.quals(a);
        char score_shift = qe.get_offset() + qe.get_from(); // subtracted from quality scores before accessing the error matrix
        for (auto i = 0; i < seed.length; i++) {
            seed.qual[i] = (uint8_t)(q[i] - score_shift); //QAO: subtract 'from'
        }

        if (use_k_mers) {

            // Add uint8_t kmer index in contiguous memory block
            size_t n_kmer = 1 << (2 * Dada2Utility::k_mer_size);
            seed.kmer8 = new uint8_t[n_kmer];
            Dada2Utility::assign_kmer8(seed.kmer8, seq.c_str(), Dada2Utility::k_mer_size);

            // Add uint16_t kmer index in contiguous memory block
            seed.kmer = new uint16_t[n_kmer];
            Dada2Utility::assign_kmer(seed.kmer, seq.c_str(), Dada2Utility::k_mer_size);

            // Add uint16_t ordered kmer record in contiguous memory block
            seed.kord = new uint16_t[seed.length];
            Dada2Utility::assign_kmer_order(seed.kord, seq.c_str(), Dada2Utility::k_mer_size);

        }

        seed.index = a;

    }


    double LightSwarmAppender::get_p_a(numSeqs_t target_mass, RawSequence& ampl) {

        double pval;

        if (ampl.abundance == 1 && !ampl.prior) { // singleton -> no abundance p-value
            pval = 1.0;
        }
        else if (ampl.comp.hamming == 0) { // cluster centre (or no mismatch to centre)
            pval = 1.0;
        }
        else if (ampl.comp.lambda == 0) { // zero expected reads of this raw
            pval = 0.0;
        }
        else { // calculate abundance p-value

            // exp_reads is the expected number of reads for this raw
            double exp_reads = ampl.comp.lambda * target_mass;
            pval = Dada2Utility::calc_p_a(ampl.abundance, exp_reads, ampl.prior);

        }

        return pval;

    }

    bool LightSwarmAppender::is_significant(RawSequence& src, numSeqs_t target_mass, numSeqs_t num_amplicons) {

        bool is_sig = false;

        // calculate the fold over-abundance and the Hamming distance to this raw
        if ((src.abundance >= min_abund) && (src.comp.hamming >= min_hamming) && (min_fold <= 1 || ((double)src.abundance) >= min_fold * src.comp.lambda * target_mass)) { // only those passing the Hamming/fold screens can be budded

            // Bonferroni correct the abundance p-value by the number of raws and compare to omega_a
            // (quite conservative, although probably unimportant given the abundance model issues)
            double p_a = src.p * num_amplicons;
            is_sig = (p_a < omega_a); // significant abundance p-value

        }

        return is_sig;

    }


    /* === LightSwarmResolver === */

    LightSwarmResolver::LightSwarmResolver(const Configuration& config) {

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
        if ((iter = misc.find("light_opt")) != misc.end()) {
            light_opt = std::stoi(iter->second);
        }

    }

    LightSwarmResolver* LightSwarmResolver::clone() const {
        return new LightSwarmResolver(*this);
    }

    void LightSwarmResolver::refine(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage, const Configuration& config) {

        char tmp_seq[amplicon_storage.max_length()];
        uint8_t tmp_qual[amplicon_storage.max_length()];
        int tvec[amplicon_storage.max_length()];
        unsigned int qind[amplicon_storage.max_length()];
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

            std::cout << " -- Error matrix filled automatically based on the quality-encoding scheme..." << std::endl;

        }

        // first use: number of light swarms / amplicons before refinement
        numSeqs_t out_cnt_s = 0;
        numSeqs_t out_cnt_a = 0;
        for (auto p = 0; p < swarm_storage.num_pools(); p++) {
            auto& out_sws = swarm_storage.get_swarms(p);
            for (auto s = 0; s < out_sws.size(); s++) {
                if (out_sws.get(s).mass() < config.boundary) out_cnt_s++;
                if (out_sws.get(s).mass() < config.boundary) out_cnt_a += out_sws.get(s).size();
            }
        }

        std::cout << " -- Number of light swarms before refinement: " << out_cnt_s << " (of " << swarm_storage.num_swarms() << ")" << std::endl;
        std::cout << " -- Number of amplicons in light swarms before refinement: " << out_cnt_a << " (of " << amplicon_storage.num_amplicons() << ")" << std::endl;

        // second use: number of ungrafted light swarms / amplicons
        out_cnt_s = 0;
        out_cnt_a = 0;

        for (numSeqs_t p = 0; p < amplicon_storage.num_pools(); p++) {

            auto& ac = amplicon_storage.get_pool(p);
            auto& swarms = swarm_storage.get_swarms(p);

            // build target clusters from heavy swarms (for access to mass / size, properties of seed)
            std::vector<TargetCluster> target_clusters;
            std::vector<numSeqs_t> light_swarms;
            numSeqs_t total_abund = 0;
            numSeqs_t total_num = ac.size();
            for (numSeqs_t s = 0; s < swarms.size(); s++) {

                auto& swarm = swarms.get(s);

                if (swarm.mass() >= config.boundary) {
                    target_clusters.emplace_back(s, ac, swarm.seed(), use_k_mers, *qe);
                } else {
                    light_swarms.emplace_back(s);
                }

                total_abund += swarm.mass();

            }

            std::vector<numSeqs_t> ungrafted_amplicons; // ids of amplicons from light swarms that were not grafted to a heavy swarm

            // for each light amplicon:
            // - compare it (using its abundance) with centres of target clusters
            // - move light amplicon (independently from other potential amplicons of the corresponding light swarm) to target cluster with highest E_minmax if not significant
            for (auto src_sid : light_swarms) {

                auto& src_swarm = swarms.get(src_sid);
                for (numSeqs_t m = 0; m < src_swarm.size(); m++) {

                    numSeqs_t aid = src_swarm.member(m);

                    std::string src_seq = ac.seq_str(aid);
                    Dada2Utility::nt2int(tmp_seq, src_seq.c_str());
                    auto quals = ac.quals(aid);
                    for (lenSeqs_t pos = 0; pos < ac.len(aid); pos++) {
                        tmp_qual[pos] = (uint8_t)(quals[pos] - score_shift);
                    }
                    RawSequence src_ampl(tmp_seq, ac.len(aid), tmp_qual, ac.ab(aid), false);

                    double max_exp = 0.0; // minimum expected abundance of the light swarm's seed, maximised over all target clusters (E_minmax in DADA2)
                    numSeqs_t max_sid = target_clusters.size(); // swarm id of the target cluster corresponding to max_exp

                    for (auto& tc : target_clusters) {

                        auto& target_swarm = swarms.get(tc.sid);
                        numSeqs_t target_seed_abund = ac.ab(target_swarm.seed());

                        Dada2Utility::Sub* sub;
                        if (greedy && (src_ampl.abundance > target_seed_abund)) {
                            sub = nullptr;
                        } else {
                            sub = new Sub(&tc.seed, &src_ampl, config.match_reward, config.mismatch_penalty,
                                          config.gap_opening_penalty + config.gap_extension_penalty, homo_gap, use_k_mers,
                                          k_dist_cutoff, band_size, gapless);
                        }

                        // calculate lambda for that sub
                        double lambda = Dada2Utility::compute_lambda(src_ampl.seq, src_ampl.length, src_ampl.qual, sub, error_matrix, tvec, qind);

                        // store comparison if potentially useful
                        if (lambda * total_abund > max_exp) { // cluster could attract this raw //TODO? total_abund = mass of target cluster OR mass of all clusters (as in DADA2)?

                            if (lambda * target_seed_abund > max_exp) { // better E_minmax -> set

                                max_exp = lambda * target_seed_abund;
                                max_sid = tc.sid;

                                src_ampl.comp.lambda = lambda;
                                src_ampl.comp.hamming = sub->nsubs;

                                src_ampl.p = get_p_a(target_swarm.mass(), src_ampl);

                            }

                        }

                        delete sub;

                    }

                    if (max_sid < target_clusters.size() && !is_significant(src_ampl, swarms.get(max_sid).mass(), total_num)) {

                        auto& target_swarm = swarms.get(max_sid);
                        Swarm& new_swarm = swarms.initialise_cluster(aid, ac.ab(aid));
                        target_swarm.attach(0, 0, &new_swarm, config.main_threshold + 1);

                    } else {
                        ungrafted_amplicons.emplace_back(aid);
                        out_cnt_a++;
                        out_cnt_s++;
                    }

                }

                src_swarm.mark_as_attached(); // "lazy delete" (will be ignored during output generation)

            }

            // processing of ungrafted light-swarm amplicons
            switch (light_opt) {

                case 1:
                    // each remaining light-swarm amplicon becomes its own swarm
                    for (auto aid : ungrafted_amplicons) {
                        swarms.initialise_cluster(aid, ac.ab(aid));
                    }
                    break;

                case 2:
                    // discard the remaining light-swarm amplicons completely by ignoring the previously determined light swarms
                    // (those members grafted to heavy swarms are handled above)
                    // nothing to do (as light swarms are already marked as attached above)
                    break;

                case 3:
                    // collect remaining light-swarm members in a single "unstructured" cluster
                    if (!ungrafted_amplicons.empty()) {
                        Swarm& new_swarm = swarms.initialise_cluster(ungrafted_amplicons[0], ac.ab(ungrafted_amplicons[0]));
                        for (auto i = 1; i < ungrafted_amplicons.size(); i++) {
                            new_swarm.append(ungrafted_amplicons[i], 1, config.main_threshold + 1, ungrafted_amplicons[0], config.main_threshold + 1, ac.ab(ungrafted_amplicons[i]), true);
                        }
                    }
                    break;

                case 4:
                    // collect remaining light-swarm members and perform a separate DADA2-like clustering on them
                    if (!ungrafted_amplicons.empty()) {
                        Dada2Clusterer::add_swarms(ac, ungrafted_amplicons, swarms, config, use_k_mers, k_dist_cutoff, band_size, omega_a,
                                                 omega_p, max_clust, min_fold, min_hamming, min_abund, homo_gap, gapless, greedy, matrix_file);
                    }
                    break;

                default:
                    std::cerr << " -- WARNING: Invalid light_opt option. Remaining light swarms are kept as separate clusters." << std::endl;
                    // same as light_opt = 1
                    for (auto aid : ungrafted_amplicons) {
                        swarms.initialise_cluster(aid, ac.ab(aid));
                    }

            }

        }

        std::cout << " -- Ungrafted light swarms: " << out_cnt_s << std::endl;
        std::cout << " -- Amplicons in ungrafted light swarms: " << out_cnt_a << std::endl;
        std::cout << " -- Processing option for ungrafted light swarms: " << light_opt << std::endl;

        // third use: number of light swarms / amplicons after refinement
        out_cnt_s = 0;
        out_cnt_a = 0;
        for (auto p = 0; p < swarm_storage.num_pools(); p++) {
            auto& out_sws = swarm_storage.get_swarms(p);
            for (auto s = 0; s < out_sws.size(); s++) {
                auto& sw = out_sws.get(s);
                if (sw.mass() < config.boundary && !sw.is_attached()) out_cnt_s++;
                if (sw.mass() < config.boundary && !sw.is_attached()) out_cnt_a += sw.size();
            }
        }

        std::cout << " -- Light swarms after refinement: " << out_cnt_s << " (of " << swarm_storage.num_swarms() << ")" << std::endl;
        std::cout << " -- Amplicons in light swarms after refinement: " << out_cnt_a << " (of " << amplicon_storage.num_amplicons() << ")" << std::endl;

        delete qe;

    }


    LightSwarmResolver::TargetCluster::TargetCluster(numSeqs_t s, const AmpliconCollection& ac, numSeqs_t a, bool use_k_mers, QualityEncoding<>& qe) : sid(s), seed() {

        seed.length = ac.len(a);
        seed.abundance = ac.ab(a);

        seed.seq = nullptr;
        seed.seq = new char[seed.length + 1];
        auto seq = ac.seq_str(a);
        Dada2Utility::nt2int(seed.seq, seq.c_str());

        seed.qual = new uint8_t[seed.length];
        auto q = ac.quals(a);
        char score_shift = qe.get_offset() + qe.get_from(); // subtracted from quality scores before accessing the error matrix
        for (auto i = 0; i < seed.length; i++) {
            seed.qual[i] = (uint8_t)(q[i] - score_shift);
        }

        if (use_k_mers) {

            // Add uint8_t kmer index in contiguous memory block
            size_t n_kmer = 1 << (2 * Dada2Utility::k_mer_size);
            seed.kmer8 = new uint8_t[n_kmer];
            Dada2Utility::assign_kmer8(seed.kmer8, seq.c_str(), Dada2Utility::k_mer_size);

            // Add uint16_t kmer index in contiguous memory block
            seed.kmer = new uint16_t[n_kmer];
            Dada2Utility::assign_kmer(seed.kmer, seq.c_str(), Dada2Utility::k_mer_size);

            // Add uint16_t ordered kmer record in contiguous memory block
            seed.kord = new uint16_t[seed.length];
            Dada2Utility::assign_kmer_order(seed.kord, seq.c_str(), Dada2Utility::k_mer_size);

        }

        seed.index = a;

    }


    double LightSwarmResolver::get_p_a(numSeqs_t target_mass, RawSequence& ampl) {

        double pval;

        if (ampl.abundance == 1 && !ampl.prior) { // singleton -> no abundance p-value
            pval = 1.0;
        }
        else if (ampl.comp.hamming == 0) { // cluster centre (or no mismatch to centre)
            pval = 1.0;
        }
        else if (ampl.comp.lambda == 0) { // zero expected reads of this raw
            pval = 0.0;
        }
        else { // calculate abundance p-value

            // exp_reads is the expected number of reads for this raw
            double exp_reads = ampl.comp.lambda * target_mass;
            pval = Dada2Utility::calc_p_a(ampl.abundance, exp_reads, ampl.prior);

        }

        return pval;

    }

    bool LightSwarmResolver::is_significant(RawSequence& src, numSeqs_t target_mass, numSeqs_t num_amplicons) {

        bool is_sig = false;

        // calculate the fold over-abundance and the Hamming distance to this raw
        if ((src.abundance >= min_abund) && (src.comp.hamming >= min_hamming) && (min_fold <= 1 || ((double)src.abundance) >= min_fold * src.comp.lambda * target_mass)) { // only those passing the Hamming/fold screens can be budded

            // Bonferroni correct the abundance p-value by the number of raws and compare to omega_a
            // (quite conservative, although probably unimportant given the abundance model issues)
            double p_a = src.p * num_amplicons;
            is_sig = (p_a < omega_a); // significant abundance p-value

        }

        return is_sig;

    }


    /* === LightSwarmShuffler === */

    LightSwarmShuffler::LightSwarmShuffler(const Configuration& config) {

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

    LightSwarmShuffler* LightSwarmShuffler::clone() const {
        return new LightSwarmShuffler(*this);
    }

    void LightSwarmShuffler::refine(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage, const Configuration& config) {

        char tmp_seq[amplicon_storage.max_length()];
        uint8_t tmp_qual[amplicon_storage.max_length()];
        int tvec[amplicon_storage.max_length()];
        unsigned int qind[amplicon_storage.max_length()];
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

            std::cout << " -- Error matrix filled automatically based on the quality-encoding scheme..." << std::endl;

        }

        // first use: number of light swarms / amplicons before refinement
        numSeqs_t out_cnt_s = 0;
        numSeqs_t out_cnt_a = 0;
        for (auto p = 0; p < swarm_storage.num_pools(); p++) {
            auto& out_sws = swarm_storage.get_swarms(p);
            for (auto s = 0; s < out_sws.size(); s++) {
                if (out_sws.get(s).mass() < config.boundary) out_cnt_s++;
                if (out_sws.get(s).mass() < config.boundary) out_cnt_a += out_sws.get(s).size();
            }
        }

        std::cout << " -- Number of light swarms before refinement: " << out_cnt_s << " (of " << swarm_storage.num_swarms() << ")" << std::endl;
        std::cout << " -- Number of amplicons in light swarms before refinement: " << out_cnt_a << " (of " << amplicon_storage.num_amplicons() << ")" << std::endl;

        // second use: number of ungrafted light swarms / amplicons = number of light swarms / amplicons after refinement
        out_cnt_s = 0;
        out_cnt_a = 0;

        for (numSeqs_t p = 0; p < amplicon_storage.num_pools(); p++) {

            auto& ac = amplicon_storage.get_pool(p);
            auto& swarms = swarm_storage.get_swarms(p);

            // build target clusters from heavy swarms (for access to mass / size, properties of seed)
            std::vector<TargetCluster> target_clusters;
            std::vector<numSeqs_t> light_swarms;
            numSeqs_t total_abund = 0;
            numSeqs_t total_num = ac.size();
            for (numSeqs_t s = 0; s < swarms.size(); s++) {

                auto& swarm = swarms.get(s);

                if (swarm.mass() >= config.boundary) {
                    target_clusters.emplace_back(s, ac, swarm.seed(), use_k_mers, *qe);
                } else {
                    light_swarms.emplace_back(s);
                }

                total_abund += swarm.mass();

            }


            // for each light amplicon:
            // - compare it (using its abundance) with centres of target clusters
            // - move light amplicon (independently from other potential amplicons of the corresponding light swarm) to target cluster with highest E_minmax if not significant
            for (auto src_sid : light_swarms) {

                auto& src_swarm = swarms.get(src_sid);
                std::vector<numSeqs_t> remaining_members; // swarm-internal member ids

                for (numSeqs_t m = src_swarm.size(); m > 0; m--) {

                    numSeqs_t aid = src_swarm.member(m - 1);

                    std::string src_seq = ac.seq_str(aid);
                    Dada2Utility::nt2int(tmp_seq, src_seq.c_str());
                    auto quals = ac.quals(aid);
                    for (lenSeqs_t pos = 0; pos < ac.len(aid); pos++) {
                        tmp_qual[pos] = (uint8_t)(quals[pos] - score_shift);
                    }
                    RawSequence src_ampl(tmp_seq, ac.len(aid), tmp_qual, ac.ab(aid), false);

                    double max_exp = 0.0; // minimum expected abundance of the light swarm's seed, maximised over all target clusters (E_minmax in DADA2)
                    numSeqs_t max_sid = target_clusters.size(); // swarm id of the target cluster corresponding to max_exp

                    for (auto& tc : target_clusters) {

                        auto& target_swarm = swarms.get(tc.sid);
                        numSeqs_t target_seed_abund = ac.ab(target_swarm.seed());

                        Dada2Utility::Sub* sub;
                        if (greedy && (src_ampl.abundance > target_seed_abund)) {
                            sub = nullptr;
                        } else {
                            sub = new Sub(&tc.seed, &src_ampl, config.match_reward, config.mismatch_penalty,
                                          config.gap_opening_penalty + config.gap_extension_penalty, homo_gap, use_k_mers,
                                          k_dist_cutoff, band_size, gapless);
                        }

                        // calculate lambda for that sub
                        double lambda = Dada2Utility::compute_lambda(src_ampl.seq, src_ampl.length, src_ampl.qual, sub, error_matrix, tvec, qind);

                        // store comparison if potentially useful
                        if (lambda * total_abund > max_exp) { // cluster could attract this raw //TODO? total_abund = mass of target cluster OR mass of all clusters (as in DADA2)?

                            if (lambda * target_seed_abund > max_exp) { // Better E_minmax -> set

                                max_exp = lambda * target_seed_abund;
                                max_sid = tc.sid;

                                src_ampl.comp.lambda = lambda;
                                src_ampl.comp.hamming = sub->nsubs;

                                src_ampl.p = get_p_a(target_swarm.mass(), src_ampl);

                            }

                        }

                        delete sub;

                    }

                    if ((aid != src_swarm.seed() || remaining_members.empty()) && max_sid < target_clusters.size() && !is_significant(src_ampl, swarms.get(max_sid).mass(), total_num)) {

                        auto& target_swarm = swarms.get(max_sid);
                        Swarm& new_swarm = swarms.initialise_cluster(aid, ac.ab(aid));
                        target_swarm.attach(0, 0, &new_swarm, config.main_threshold + 1);

                    } else {
                        remaining_members.emplace_back(m - 1);
                    }

                }

                // add new reduced (and "restructured") light swarm if any member left, delete old version (always)
                if (!remaining_members.empty()) {

                    auto& new_swarm = swarms.initialise_cluster(src_swarm.seed(), ac.ab(src_swarm.seed()));
                    for (auto iter = remaining_members.rbegin() + 1; iter != remaining_members.rend(); iter++) {
                        new_swarm.append(src_swarm.member(*iter), 1, config.main_threshold + 1, src_swarm.seed(),
                                config.main_threshold + 1, ac.ab(src_swarm.member(*iter)), true);
                    }

                    out_cnt_a += remaining_members.size();
                    out_cnt_s++;

                }
                src_swarm.mark_as_attached(); // "lazy delete" (will be ignored during output generation)

            }

            // processing of ungrafted light-swarm amplicons is fixed: keep remaining light swarms as they are
            // nothing to do

        }

        std::cout << " -- Light swarms after refinement: " << out_cnt_s << " (of " << swarm_storage.num_swarms() << ")" << std::endl;
        std::cout << " -- Amplicons in light swarms after refinement: " << out_cnt_a << " (of " << amplicon_storage.num_amplicons() << ")" << std::endl;

        delete qe;

    }


    LightSwarmShuffler::TargetCluster::TargetCluster(numSeqs_t s, const AmpliconCollection& ac, numSeqs_t a, bool use_k_mers, QualityEncoding<>& qe) : sid(s), seed() {

        seed.length = ac.len(a);
        seed.abundance = ac.ab(a);

        seed.seq = nullptr;
        seed.seq = new char[seed.length + 1];
        auto seq = ac.seq_str(a);
        Dada2Utility::nt2int(seed.seq, seq.c_str());

        seed.qual = new uint8_t[seed.length];
        auto q = ac.quals(a);
        char score_shift = qe.get_offset() + qe.get_from(); // subtracted from quality scores before accessing the error matrix
        for (auto i = 0; i < seed.length; i++) {
            seed.qual[i] = (uint8_t)(q[i] - score_shift);
        }

        if (use_k_mers) {

            // Add uint8_t kmer index in contiguous memory block
            size_t n_kmer = 1 << (2 * Dada2Utility::k_mer_size);
            seed.kmer8 = new uint8_t[n_kmer];
            Dada2Utility::assign_kmer8(seed.kmer8, seq.c_str(), Dada2Utility::k_mer_size);

            // Add uint16_t kmer index in contiguous memory block
            seed.kmer = new uint16_t[n_kmer];
            Dada2Utility::assign_kmer(seed.kmer, seq.c_str(), Dada2Utility::k_mer_size);

            // Add uint16_t ordered kmer record in contiguous memory block
            seed.kord = new uint16_t[seed.length];
            Dada2Utility::assign_kmer_order(seed.kord, seq.c_str(), Dada2Utility::k_mer_size);

        }

        seed.index = a;

    }


    double LightSwarmShuffler::get_p_a(numSeqs_t target_mass, RawSequence& ampl) {

        double pval;

        if (ampl.abundance == 1 && !ampl.prior) { // singleton -> no abundance p-value
            pval = 1.0;
        }
        else if (ampl.comp.hamming == 0) { // cluster centre (or no mismatch to centre)
            pval = 1.0;
        }
        else if (ampl.comp.lambda == 0) { // zero expected reads of this raw
            pval = 0.0;
        }
        else { // calculate abundance p-value

            // exp_reads is the expected number of reads for this raw
            double exp_reads = ampl.comp.lambda * target_mass;
            pval = Dada2Utility::calc_p_a(ampl.abundance, exp_reads, ampl.prior);

        }

        return pval;

    }

    bool LightSwarmShuffler::is_significant(RawSequence& src, numSeqs_t target_mass, numSeqs_t num_amplicons) {

        bool is_sig = false;

        // calculate the fold over-abundance and the Hamming distance to this raw
        if ((src.abundance >= min_abund) && (src.comp.hamming >= min_hamming) && (min_fold <= 1 || ((double)src.abundance) >= min_fold * src.comp.lambda * target_mass)) { // only those passing the Hamming/fold screens can be budded

            // Bonferroni correct the abundance p-value by the number of raws and compare to omega_a
            // (quite conservative, although probably unimportant given the abundance model issues)
            double p_a = src.p * num_amplicons;
            is_sig = (p_a < omega_a); // significant abundance p-value

        }

        return is_sig;

    }


    /* === SwarmShuffler === */

    SwarmShuffler::SwarmShuffler(const Configuration& config) {

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

    SwarmShuffler* SwarmShuffler::clone() const {
        return new SwarmShuffler(*this);
    }

    void SwarmShuffler::refine(const AmpliconStorage& amplicon_storage, SwarmStorage& swarm_storage, const Configuration& config) {

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

            std::cout << " -- Error matrix filled automatically based on the quality-encoding scheme..." << std::endl;

        }

        // first use: number of light swarms / amplicons before refinement
        numSeqs_t out_cnt_s = 0;
        numSeqs_t out_cnt_a = 0;
        for (auto p = 0; p < swarm_storage.num_pools(); p++) {
            auto& out_sws = swarm_storage.get_swarms(p);
            for (auto s = 0; s < out_sws.size(); s++) {
                if (out_sws.get(s).mass() < config.boundary) out_cnt_s++;
                if (out_sws.get(s).mass() < config.boundary) out_cnt_a += out_sws.get(s).size();
            }
        }

        std::cout << " -- Number of light swarms before refinement: " << out_cnt_s << " (of " << swarm_storage.num_swarms() << ")" << std::endl;
        std::cout << " -- Number of amplicons in light swarms before refinement: " << out_cnt_a << " (of " << amplicon_storage.num_amplicons() << ")" << std::endl;

        for (numSeqs_t p = 0; p < amplicon_storage.num_pools(); p++) {

            auto& ac = amplicon_storage.get_pool(p);
            auto& swarms = swarm_storage.get_swarms(p);

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

                std::string src_seq = ac.seq_str(i);
                Dada2Utility::nt2int(seq, src_seq.c_str());

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

            ClusterCollection* clusters = shuffle(raws, swarms, ac.size(), error_matrix, config.match_reward,
                                                   config.mismatch_penalty, config.gap_opening_penalty + config.gap_extension_penalty, homo_gap);

            // delete old swarms of the current pool and build new ones from the DADA2 shuffling
            swarms.clear();
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

        // second use: number of light swarms / amplicons after refinement
        out_cnt_s = 0;
        out_cnt_a = 0;
        for (auto p = 0; p < swarm_storage.num_pools(); p++) {
            auto& out_sws = swarm_storage.get_swarms(p);
            for (auto s = 0; s < out_sws.size(); s++) {
                auto& sw = out_sws.get(s);
                if (sw.mass() < config.boundary && !sw.is_attached()) out_cnt_s++;
                if (sw.mass() < config.boundary && !sw.is_attached()) out_cnt_a += sw.size();
            }
        }

        std::cout << " -- Light swarms after refinement: " << out_cnt_s << " (of " << swarm_storage.num_swarms() << ")" << std::endl;
        std::cout << " -- Amplicons in light swarms after refinement: " << out_cnt_a << " (of " << amplicon_storage.num_amplicons() << ")" << std::endl;

        delete qe;

    }

    SwarmShuffler::Cluster::Cluster(unsigned int num) {

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

    SwarmShuffler::Cluster::~Cluster() {

        delete[] seq;
        delete[] raws;

    }

    unsigned int SwarmShuffler::Cluster::add_raw(RawSequence* raw) {

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

    double SwarmShuffler::Cluster::get_p_a(RawSequence* raw) {

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

    bool SwarmShuffler::Cluster::is_significant(RawSequence& src, numSeqs_t num_amplicons, double min_fold, int min_hamming, int min_abund, double omega_a) {

        bool is_sig = false;

        // calculate the fold over-abundance and the Hamming distance to this raw
        if ((src.abundance >= min_abund) && (src.comp.hamming >= min_hamming) && (min_fold <= 1 || ((double)src.abundance) >= min_fold * src.comp.lambda * mass)) { // only those passing the Hamming/fold screens can be budded

            // Bonferroni correct the abundance p-value by the number of raws and compare to omega_a
            // (quite conservative, although probably unimportant given the abundance model issues)
            double p_a = src.p * num_amplicons;
            is_sig = (p_a < omega_a); // significant abundance p-value

        }

        return is_sig;

    }

    SwarmShuffler::RawSequence* SwarmShuffler::Cluster::pop_raw(unsigned int r) {

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

    SwarmShuffler::ClusterCollection::ClusterCollection(RawSequence** raw_sequences, Swarms& initial_swarms, unsigned int raws_size, double om_a, double om_p) {

        // allocate memory
        clusters = new Cluster*[initial_swarms.size()];

        // initialise basic values
        num_clusters = initial_swarms.size();
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

        // initialise with given clusters
        for (numSeqs_t s = 0; s < initial_swarms.size(); s++) {

            auto& swarm = initial_swarms.get(s);

            clusters[s] = new Cluster(s);
            for (numSeqs_t i = 0; i < swarm.size(); i++) {
                clusters[s]->add_raw(raws[swarm.member(i)]);
            }

            // assign the raw with the most reads as the centre
            clusters[s]->centre = clusters[s]->raws[0];
            delete[] clusters[s]->seq;
            clusters[s]->seq = new char[clusters[s]->centre->length + 1];
            strcpy(clusters[s]->seq, clusters[s]->centre->seq);

        }

    }

    SwarmShuffler::ClusterCollection::~ClusterCollection() {

        delete[] qind;
        delete[] tvec;
        for (int i = 0; i < num_clusters; i++) {
            delete clusters[i];
        }
        delete[] clusters;

    }

    void SwarmShuffler::ClusterCollection::p_update(bool greedy) {

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
                    double exp_reads_centre = cluster->centre->abundance * raw_sequence->comp.lambda;
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

    bool SwarmShuffler::ClusterCollection::shuffle(double min_fold, int min_hamming, int min_abund) {

        double* emax = new double[num_raws]; // highest expected abundance for each raw
        Comparison** compmax = new Comparison*[num_raws]; // pointer to comparison leading to emax-value for each raw

        // initialise emax / compmax off of cluster 0
        // comparisons to all raws exist in cluster 0 (in index order)
        for (numSeqs_t c = 0; c < num_clusters; c++) {

            for (numSeqs_t i = 0; i < clusters[c]->size; i++) {

                numSeqs_t index = clusters[c]->raws[i]->index;
                compmax[index] = &clusters[c]->comps[i];
                emax[index] = compmax[index]->lambda * clusters[c]->mass;

            }

        }

        // iterate over comparisons and find comparison with best expectation for each raw
        for (unsigned int i = 0; i < num_clusters; i++) {

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
                if (compmax[raw_sequence->index]->i != i && !clusters[compmax[raw_sequence->index]->i]->is_significant(*raw_sequence, num_raws, min_fold, min_hamming, min_abund, omega_a)) {

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

    void SwarmShuffler::ClusterCollection::compare(unsigned int i, std::vector<std::vector<double>>& error_matrix,
            int match, int mismatch, int gap_pen, int homo_gap_pen,
            bool use_k_mers, double k_dist_cutoff, int band_size, bool gapless, bool greedy) {

        // align all raws to this sequence and compute corresponding lambda
        unsigned int centre_reads = clusters[i]->centre->abundance;

        for (unsigned int index = 0; index < num_raws; index++) {

            RawSequence* raw_sequence = raws[index];
            Sub* sub;

            //TODO? speed-up: check whether raw is already in cluster c and continue to next raw in this case

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
                if (/*i == 0 || */raw_sequence == clusters[i]->centre) { // update /*on init (i = 0) or*/ if the centre
                    raw_sequence->comp = comp;
                }

            }

            delete sub;

        }

    }

    void SwarmShuffler::ClusterCollection::compare(std::vector<std::vector<double>>& error_matrix,
            int match, int mismatch, int gap_pen, int homo_gap_pen,
            bool use_k_mers, double k_dist_cutoff, int band_size, bool gapless, bool greedy) {

        for (numSeqs_t c = 0; c < num_clusters; c++) {

            Cluster* cluster = clusters[c];

            // align all raws to this sequence and compute corresponding lambda
            unsigned int centre_reads = cluster->centre->abundance;

            for (unsigned int index = 0; index < cluster->size; index++) {

                RawSequence* raw_sequence = cluster->raws[index];
                Sub* sub;

                // get sub object
                if (greedy && (raw_sequence->abundance > centre_reads)) {
                    sub = nullptr;
                } else if (greedy && raw_sequence->lock) { // was (raw->p > 0.5)
                    sub = nullptr;
                } else {
                    sub = new Sub(cluster->centre, raw_sequence, match, mismatch, gap_pen, homo_gap_pen, use_k_mers,
                                  k_dist_cutoff, band_size, gapless);
                }

                // calculate lambda for that sub
                double lambda = Dada2Utility::compute_lambda(raw_sequence->seq, raw_sequence->length, raw_sequence->qual, sub, error_matrix, tvec, qind);

                // store comparison if potentially useful
                if (lambda * total_mass > raw_sequence->max_exp) { // this cluster could attract this raw

                    if (lambda * cluster->centre->abundance > raw_sequence->max_exp) { // better max_exp
                        raw_sequence->max_exp = lambda * cluster->centre->abundance;
                    }

                    Comparison comp;
                    comp.i = c;
                    comp.index = raw_sequence->index;
                    comp.lambda = lambda;
                    comp.hamming = sub->nsubs;

                    cluster->comps.push_back(comp);
                    raw_sequence->comp = comp;

                }

                delete sub;

            }

        }

    }

    SwarmShuffler::ClusterCollection* SwarmShuffler::shuffle(RawSequence** raws, Swarms& initial_swarms, int num_raw,
            std::vector<std::vector<double>>& error_matrix, int match, int mismatch, int gap_pen, int homo_gap_pen) {

        auto bb = new ClusterCollection(raws, initial_swarms, num_raw, omega_a, omega_p); // new cluster with all sequences in one cluster

        // everyone gets aligned within its initial cluster, no k-mer screen
        bb->compare(error_matrix, match, mismatch, gap_pen, homo_gap_pen, false, 1.0, band_size, gapless, greedy);
        bb->p_update(greedy); // calculates abundance p-value for each raw in its cluster (consensuses)

        // keep shuffling and updating until no more shuffles
        int nshuffle = 0;
        bool shuffled;
        do {

            for (numSeqs_t c = 0; c < bb->num_clusters; c++) {
                bb->compare(c, error_matrix, match, mismatch, gap_pen, homo_gap_pen, use_k_mers, k_dist_cutoff, band_size, gapless, greedy);
            }
            shuffled = bb->shuffle(min_fold, min_hamming, min_abund);

        } while (shuffled && ++nshuffle < Dada2Utility::max_shuffle);
        if (nshuffle >= Dada2Utility::max_shuffle) {
            printf("Warning: Reached maximum (%i) shuffles.\n", Dada2Utility::max_shuffle);
        }

        bb->p_update(greedy);

        return bb;

    }


}