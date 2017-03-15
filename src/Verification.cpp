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

#include "../include/Verification.hpp"


namespace GeFaST {

// ===== Exact computation =====

lenSeqs_t Verification::computeClassicFull(const std::string& s, const std::string& t) {

    lenSeqs_t M[s.size() + 1][t.size() + 1];

    // initialise first column and row
    for (lenSeqs_t i = 0; i <= s.size(); i++) {
        M[i][0] = i;
    }
    for (lenSeqs_t j = 1; j <= t.size(); j++) {
        M[0][j] = j;
    }

    // compute remaining rows
    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        for (lenSeqs_t j = 1; j <= t.size(); j++) {

            M[i][j] = std::min({
                                       M[i - 1][j - 1] + (s[i - 1] != t[j - 1]), // (mis)match
                                       M[i - 1][j] + 1, // deletion
                                       M[i][j - 1] + 1 // insertion
                               });

        }

    }

    return M[s.size()][t.size()];

}


lenSeqs_t Verification::computeClassicRow(const std::string& s, const std::string& t, lenSeqs_t* M) {

//    lenSeqs_t M[t.size() + 1];
    lenSeqs_t match, tmp;

    // initialise first row
    for (lenSeqs_t j = 0; j <= t.size(); j++) {
        M[j] = j;
    }

    // compute remaining rows
    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        match = i - 1; // value of the epsilon column in the current row
        M[0] = i;

        for (lenSeqs_t j = 1; j <= t.size(); j++) {

            tmp = std::min({
                                   match + (s[i - 1] != t[j - 1]), // (mis)match
                                   M[j] + 1, // deletion
                                   M[j - 1] + 1 // insertion
                           });
            match = M[j];
            M[j] = tmp;

        }

    }

    return M[t.size()];

}


// ===== Bounded computation =====

lenSeqs_t Verification::mapColIndex(const lenSeqs_t j, const lenSeqs_t i, const lenSeqs_t bound) {
    return j - (i > bound) * (i - bound);
}


// rationale for condition of the first case of the entry calculation: avoid underflow (due to unsigned integer type) by basically "disabling" this case for rows where it should not occur

lenSeqs_t Verification::computeBoundedFull(const std::string& s, const std::string& t, const lenSeqs_t bound) {

    // long computation not necessary if lengths differ too much
    if (((s.size() > t.size()) ? (s.size() - t.size()) : (t.size() - s.size())) > bound) {
        return bound + 1;
    }

    if (bound == 0) {
        return lenSeqs_t(s != t);
    }


    lenSeqs_t M[s.size() + 1][t.size() + 1];

    // initialise necessary sections of first column and row
    for (lenSeqs_t i = 0; i <= bound && i <= s.size(); i++) {
        M[i][0] = i;
    }
    for (lenSeqs_t j = 1; j <= bound && j <= t.size(); j++) {
        M[0][j] = j;
    }


    bool early;

    // compute sections of remaining rows
    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        early = true; // early termination flag

        for (lenSeqs_t j = 1 + (i > bound) * (i - bound - 1); j <= i + bound && j <= t.size(); j++) { // same as starting from j = max(1, i - bound) with signed integers

            if (j == ((i > bound) * (i - bound))) {

                M[i][j] = std::min({
                                           M[i - 1][j - 1] + (s[i - 1] != t[j - 1]), // (mis)match
                                           M[i - 1][j] + 1 // deletion
                                   });

            } else if (j == (i + bound)) {

                M[i][j] = std::min({
                                           M[i - 1][j - 1] + (s[i - 1] != t[j - 1]), // (mis)match
                                           M[i][j - 1] + 1 // insertion
                                   });

            } else {

                M[i][j] = std::min({
                                           M[i - 1][j - 1] + (s[i - 1] != t[j - 1]), // (mis)match
                                           M[i - 1][j] + 1, // deletion
                                           M[i][j - 1] + 1 // insertion
                                   });

            }

            early &= (M[i][j] > bound);

        }

        if (early) { // computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }

    }

    return (M[s.size()][t.size()] > bound) ? (bound + 1) : M[s.size()][t.size()];

}

lenSeqs_t Verification::computeBoundedRow(const std::string& s, const std::string& t, const lenSeqs_t bound, lenSeqs_t* M) {

    // long computation not necessary if lengths differ too much
    if (((s.size() > t.size()) ? (s.size() - t.size()) : (t.size() - s.size())) > bound) {
        return bound + 1;
    }

    if (bound == 0) {
        return lenSeqs_t(s != t);
    }


//    lenSeqs_t M[t.size() + 1];
    lenSeqs_t match, tmp;

    // initialise necessary section of first row
    for (lenSeqs_t j = 0; j <= bound && j <= t.size(); j++) {
        M[j] = j;
    }

    bool early;

    // compute sections of remaining rows
    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        early = true; // early termination flag

        match = (i <= bound + 1) * (i - 1)  + (i > bound + 1) * M[i - bound - 1]; // same as: match = (i <= bound + 1) ? (i - 1) : M[i - bound - 1]
        M[0] = i;

        for (lenSeqs_t j = 1 + (i > bound) * (i - bound - 1); j <= i + bound && j <= t.size(); j++) { // same as starting from j = max(1, i - bound) with signed integers

            if (j == ((i > bound) * (i - bound))) {

                tmp = std::min({
                                       match + (s[i - 1] != t[j - 1]), // (mis)match
                                       M[j] + 1 // deletion
                               });

            } else if (j == (i + bound)) {

                tmp = std::min({
                                       match + (s[i - 1] != t[j - 1]), // (mis)match
                                       M[j - 1] + 1 // insertion
                               });

            } else {

                tmp = std::min({
                                       match + (s[i - 1] != t[j - 1]), // (mis)match
                                       M[j] + 1, // deletion
                                       M[j - 1] + 1 // insertion
                               });

            }

            match = M[j];
            M[j] = tmp;

            early &= (M[j] > bound);

        }

        if (early) { // computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }

    }

    return (M[t.size()] > bound) ? (bound + 1) : M[t.size()];

}

lenSeqs_t Verification::computeBoundedFullSlim(const std::string& s, const std::string& t, const lenSeqs_t bound) {

    // long computation not necessary if lengths differ too much
    if (((s.size() > t.size()) ? (s.size() - t.size()) : (t.size() - s.size())) > bound) {
        return bound + 1;
    }

    if (bound == 0) {
        return lenSeqs_t(s != t);
    }


    lenSeqs_t M[s.size() + 1][std::min(t.size() + 1, 2 * bound + 1)];

    // initialise necessary sections of first column and row
    for (lenSeqs_t i = 0; i <= bound && i <= s.size(); i++) {
        M[i][0] = i;
    }
    for (lenSeqs_t j = 1; j <= bound && j <= t.size(); j++) {
        M[0][j] = j;
    }


    bool early;

    // compute sections of remaining rows
    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        early = true; // early termination flag

        for (lenSeqs_t j = 1 + (i > bound) * (i - bound - 1); j <= i + bound && j <= t.size(); j++) { // same as starting from j = max(1, i - bound) with signed integers

            if (j == ((i > bound) * (i - bound))) {

                M[i][mapColIndex(j, i, bound)] = std::min({
                                                                M[i - 1][mapColIndex(j - 1, i - 1, bound)] + (s[i - 1] != t[j - 1]), // (mis)match
                                                                M[i - 1][mapColIndex(j, i - 1, bound)] + 1 // deletion
                                                        });

            } else if (j == (i + bound)) {

                M[i][mapColIndex(j, i, bound)] = std::min({
                                                                M[i - 1][mapColIndex(j - 1, i - 1, bound)] + (s[i - 1] != t[j - 1]), // (mis)match
                                                                M[i][mapColIndex(j, i, bound) - 1] + 1 // insertion
                                                        });

            } else {

                M[i][mapColIndex(j, i, bound)] = std::min({
                                                                M[i - 1][mapColIndex(j - 1, i - 1, bound)] + (s[i - 1] != t[j - 1]), // (mis)match
                                                                M[i - 1][mapColIndex(j, i - 1, bound)] + 1, // deletion
                                                                M[i][mapColIndex(j, i, bound) - 1] + 1 // insertion
                                                        });

            }

            early &= (M[i][mapColIndex(j, i, bound)] > bound);

        }
        if (early) { // computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }

    }

    return (M[s.size()][t.size() - (s.size() > bound) * (s.size() - bound)] > bound) ? (bound + 1) : M[s.size()][t.size() - (s.size() > bound) * (s.size() - bound)];

}

lenSeqs_t Verification::computeBoundedRowSlim(const std::string& s, const std::string& t, const lenSeqs_t bound, lenSeqs_t* M) {

    // long computation not necessary if lengths differ too much
    if (((s.size() > t.size()) ? (s.size() - t.size()) : (t.size() - s.size())) > bound) {
        return bound + 1;
    }

    if (bound == 0) {
        return lenSeqs_t(s != t);
    }


//    lenSeqs_t M[std::min(t.size() + 1, 2 * bound + 1)];
    lenSeqs_t match, tmp;

    // initialise necessary section of first row
    for (lenSeqs_t j = 0; j <= bound && j <= t.size(); j++) {
        M[j] = j;
    }

    bool early;

    // compute sections of remaining rows
    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        early = true; // early termination flag

        match = (i <= bound) * (i - 1)  + (i > bound) * M[mapColIndex(i - bound - 1, i - 1, bound)]; // same as: match = (i <= bound) ? (i - 1) : M[mapColIndex(i - bound - 1, i - 1, bound)]
        M[0] = i;


        for (lenSeqs_t j = 1 + (i > bound) * (i - bound - 1); j <= i + bound && j <= t.size(); j++) { // same as starting from j = max(1, i - bound) with signed integers

            if (j == ((i > bound) * (i - bound))) {

                tmp = std::min({
                                       match + (s[i - 1] != t[j - 1]), // (mis)match
                                       M[mapColIndex(j, i - 1, bound)] + 1 // deletion
                               });

            } else if (j == (i + bound)) {

                tmp = std::min({
                                       match + (s[i - 1] != t[j - 1]), // (mis)match
                                       M[mapColIndex(j - 1, i, bound)] + 1 // insertion
                               });

            } else {

                tmp = std::min({
                                       match + (s[i - 1] != t[j - 1]), // (mis)match
                                       M[mapColIndex(j, i - 1, bound)] + 1, // deletion
                                       M[mapColIndex(j - 1, i, bound)] + 1 // insertion
                               });

            }

            match = M[mapColIndex(j, i - 1, bound)];
            M[mapColIndex(j, i, bound)] = tmp;



            early &= (M[mapColIndex(j, i, bound)] > bound);

        }

        if (early) { // computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }

    }

    return (M[mapColIndex(t.size(), s.size(), bound)] > bound) ? (bound + 1) : M[mapColIndex(t.size(), s.size(), bound)];

}


// ===== Length-aware bounded computation =====


// rationale for condition of the first case of the entry calculation: avoid underflow (due to unsigned integer type) by basically "disabling" this case for rows where it should not occur

lenSeqs_t Verification::computeLengthAwareFull(const std::string& s, const std::string& t, const lenSeqs_t bound) {

    // long computation not necessary if lengths differ too much
    if (((s.size() > t.size()) ? (s.size() - t.size()) : (t.size() - s.size())) > bound) {
        return bound + 1;
    }

    if (bound == 0) {
        return lenSeqs_t(s != t);
    }


    std::string shorter = (s.size() < t.size()) ? s : t;
    std::string longer = (s.size() >= t.size()) ? s : t;
    lenSeqs_t diff = longer.size() - shorter.size();

    lenSeqs_t M[shorter.size() + 1][longer.size() + 1];

    // initialise necessary sections of first column and row
    for (lenSeqs_t i = 0; i <= (bound - diff) / 2 && i <= shorter.size(); i++) {
        M[i][0] = i;
    }
    for (lenSeqs_t j = 1; j <= (bound + diff) / 2 && j <= longer.size(); j++) {
        M[0][j] = j;
    }


    bool early;

    // compute sections of remaining rows
    for (lenSeqs_t i = 1; i <= shorter.size(); i++) {

        early = true; // early termination flag

        for (lenSeqs_t j = 1 + (i > (bound - diff) / 2) * (i - (bound - diff) / 2 - 1); j <= i + (bound + diff) / 2 && j <= longer.size(); j++) { // same as starting from j = max(1, i - (bound - diff) / 2) with signed integers

            if ((bound - diff) / 2 == 0 && (bound + diff) / 2 == 0) {

                M[i][j] = M[i - 1][j - 1] + (shorter[i - 1] != longer[j - 1]); // (mis)match is only possibility since we have to consider only one diagonal

            } else if (j == ((i > (bound - diff) / 2) * (i - (bound - diff) / 2))) {

                M[i][j] = std::min({
                                           M[i - 1][j - 1] + (shorter[i - 1] != longer[j - 1]), // (mis)match
                                           M[i - 1][j] + 1 // deletion
                                   });

            } else if (j == (i + (bound + diff) / 2)) {

                M[i][j] = std::min({
                                           M[i - 1][j - 1] + (shorter[i - 1] != longer[j - 1]), // (mis)match
                                           M[i][j - 1] + 1 // insertion
                                   });

            } else {

                M[i][j] = std::min({
                                           M[i - 1][j - 1] + (shorter[i - 1] != longer[j - 1]), // (mis)match
                                           M[i - 1][j] + 1, // deletion
                                           M[i][j - 1] + 1 // insertion
                                   });

            }

            early &= ((M[i][j] + ((diff + i >= j) * (diff + i - j) + (diff + i < j) * (j - diff - i))) > bound); // improved e.t.

        }

        if (early) { // computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }

    }

    return (M[shorter.size()][longer.size()] > bound) ? (bound + 1) : M[shorter.size()][longer.size()];

}

lenSeqs_t Verification::computeLengthAwareRow(const std::string& s, const std::string& t, const lenSeqs_t bound, lenSeqs_t* M) {

    // long computation not necessary if lengths differ too much
    if (((s.size() > t.size()) ? (s.size() - t.size()) : (t.size() - s.size())) > bound) {
        return bound + 1;
    }

    if (bound == 0) {
        return lenSeqs_t(s != t);
    }


    std::string shorter = s;
    std::string longer = t;
    if (shorter.size() > longer.size()) shorter.swap(longer);
    lenSeqs_t diff = longer.size() - shorter.size();


    lenSeqs_t match, tmp;

    // (mis)match is only possibility when we have to consider only one diagonal [happens only if (a) bound = delta = 0, or (b) bound = 1 and delta = 0, but (a) is already covered above]
    if ((bound - diff) / 2 == 0 && (bound + diff) / 2 == 0) {

        lenSeqs_t diffs = 0;
        for (auto i = 0; diffs <= bound && i < shorter.size(); i++) {
            diffs += (shorter[i] != longer[i]);
        }

        return diffs;

    }

    // initialise necessary sections of first row
    for (lenSeqs_t j = 0; j <= (bound + diff) / 2 && j <= longer.size(); j++) {
        M[j] = j;
    }

    lenSeqs_t j;
    bool early;

    // compute sections of remaining rows
    for (lenSeqs_t i = 1; i <= shorter.size(); i++) {

        early = true; // early termination flag

        j = 1 + (i > (bound - diff) / 2) * (i - (bound - diff) / 2 - 1);
        match = M[j - 1];
        M[j - 1] = (i <= (bound - diff) / 2) ? i : POS_INF; // handle left end to avoid case distinction
        if (i + (bound + diff) / 2 <= longer.size()) { // handle right end to avoid case distinction
            M[i + (bound + diff) / 2] = POS_INF;
        }

        for (; j <= i + (bound + diff) / 2 && j <= longer.size(); j++) { // same as starting from j = max(1, i - (bound - diff) / 2) with signed integers

            tmp = std::min({
                                   match + (shorter[i - 1] != longer[j - 1]), // (mis)match
                                   M[j] + 1, // deletion
                                   M[j - 1] + 1 // insertion
                           });

            match = M[j];
            M[j] = tmp;

            early &= ((M[j] + llabs((long long)diff + (long long)i - (long long)j)) > bound); // improved e.t.

        }

        if (early) { // computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }


    }

    return (M[longer.size()] > bound) ? (bound + 1) : M[longer.size()];

}

lenSeqs_t Verification::computeLengthAwareFullSlim(const std::string& s, const std::string& t, const lenSeqs_t bound) {

    // long computation not necessary if lengths differ too much
    if (((s.size() > t.size()) ? (s.size() - t.size()) : (t.size() - s.size())) > bound) {
        return bound + 1;
    }

    if (bound == 0) {
        return lenSeqs_t(s != t);
    }


    std::string shorter = (s.size() < t.size()) ? s : t;
    std::string longer = (s.size() >= t.size()) ? s : t;
    lenSeqs_t diff = longer.size() - shorter.size();

    lenSeqs_t M[shorter.size() + 1][std::min(longer.size() + 1, (bound - diff) / 2 + 1 + (bound + diff) / 2)];

    // initialise necessary sections of first column and row
    for (lenSeqs_t i = 0; i <= (bound - diff) / 2 && i <= shorter.size(); i++) {
        M[i][0] = i;
    }
    for (lenSeqs_t j = 1; j <= (bound + diff) / 2 && j <= longer.size(); j++) {
        M[0][j] = j;
    }


    bool early;

    // compute sections of remaining rows
    for (lenSeqs_t i = 1; i <= shorter.size(); i++) {

        early = true; // early termination flag

        for (lenSeqs_t j = 1 + (i > (bound - diff) / 2) * (i - (bound - diff) / 2 - 1); j <= i + (bound + diff) / 2 && j <= longer.size(); j++) { // same as starting from j = max(1, i - (bound - diff) / 2) with signed integers

            if ((bound - diff) / 2 == 0 && (bound + diff) / 2 == 0) {

                M[i][mapColIndex(j, i, (bound - diff) / 2)] = M[i - 1][mapColIndex(j - 1, i - 1, (bound - diff) / 2)] + (shorter[i - 1] != longer[j - 1]); // (mis)match is only possibility since we have to consider only one diagonal

            } else if (j == ((i > (bound - diff) / 2) * (i - (bound - diff) / 2))) {

                M[i][mapColIndex(j, i, (bound - diff) / 2)] = std::min({
                                                                           M[i - 1][mapColIndex(j - 1, i - 1, (bound - diff) / 2)] + (shorter[i - 1] != longer[j - 1]), // (mis)match
                                                                           M[i - 1][mapColIndex(j, i - 1, (bound - diff) / 2)] + 1 // deletion
                                                                   });

            } else if (j == (i + (bound + diff) / 2)) {

                M[i][mapColIndex(j, i, (bound - diff) / 2)] = std::min({
                                                                           M[i - 1][mapColIndex(j - 1, i - 1, (bound - diff) / 2)] + (shorter[i - 1] != longer[j - 1]), // (mis)match
                                                                           M[i][mapColIndex(j - 1, i, (bound - diff) / 2)] + 1 // insertion
                                                                   });

            } else {

                M[i][mapColIndex(j, i, (bound - diff) / 2)] = std::min({
                                                                           M[i - 1][mapColIndex(j - 1, i - 1, (bound - diff) / 2)] + (shorter[i - 1] != longer[j - 1]), // (mis)match
                                                                           M[i - 1][mapColIndex(j, i - 1, (bound - diff) / 2)] + 1, // deletion
                                                                           M[i][mapColIndex(j - 1, i, (bound - diff) / 2)] + 1 // insertion
                                                                   });

            }

            early &= ((M[i][mapColIndex(j, i, (bound - diff) / 2)] + ((diff + i >= j) * (diff + i - j) + (diff + i < j) * (j - diff - i))) > bound); // improved e.t.

        }
        if (early) { // computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }

    }

    return (M[shorter.size()][mapColIndex(longer.size(), shorter.size(), (bound - diff) / 2)] > bound) ? (bound + 1) : M[shorter.size()][mapColIndex(longer.size(), shorter.size(), (bound - diff) / 2)];

}

lenSeqs_t Verification::computeLengthAwareRowSlim(const std::string& s, const std::string& t, const lenSeqs_t bound, lenSeqs_t* M) {

    // long computation not necessary if lengths differ too much
    if (((s.size() > t.size()) ? (s.size() - t.size()) : (t.size() - s.size())) > bound) {
        return bound + 1;
    }

    if (bound == 0) {
        return lenSeqs_t(s != t);
    }


    std::string shorter = (s.size() < t.size()) ? s : t;
    std::string longer = (s.size() >= t.size()) ? s : t;
    lenSeqs_t diff = longer.size() - shorter.size();


//    lenSeqs_t M[std::min(longer.size() + 1, (bound - diff) / 2 + 1 + (bound + diff) / 2)];
    lenSeqs_t match, tmp;

    // initialise necessary sections of first row
    for (lenSeqs_t j = 0; j <= (bound + diff) / 2 && j <= longer.size(); j++) {
        M[j] = j;
    }

    bool early;

    // compute sections of remaining rows
    for (lenSeqs_t i = 1; i <= shorter.size(); i++) {

        early = true; // early termination flag

        match = (i <= (bound - diff) / 2) * (i - 1) + (i > (bound - diff) / 2) * M[mapColIndex(i - (bound - diff) / 2 - 1, i - 1, (bound - diff) / 2)]; // same as match = (i <= (bound - diff) / 2) ? (i - 1) : M[mapColIndex(i - (bound - diff) / 2 - 1, i - 1, (bound - diff) / 2)]
        M[0] = i;

        for (lenSeqs_t j = 1 + (i > (bound - diff) / 2) * (i - (bound - diff) / 2 - 1); j <= i + (bound + diff) / 2 && j <= longer.size(); j++) { // same as starting from j = max(1, i - (bound - diff) / 2) with signed integers

            if ((bound - diff) / 2 == 0 && (bound + diff) / 2 == 0) {

                tmp = match + (shorter[i - 1] != longer[j - 1]); // (mis)match is only possibility since we have to consider only one diagonal

            } else if (j == ((i > (bound - diff) / 2) * (i - (bound - diff) / 2))) {

                tmp = std::min({
                                       match + (shorter[i - 1] != longer[j - 1]), // (mis)match
                                       M[mapColIndex(j, i - 1, (bound - diff) / 2)] + 1 // deletion
                               });

            } else if (j == (i + (bound + diff) / 2)) {

                tmp = std::min({
                                       match + (shorter[i - 1] != longer[j - 1]), // (mis)match
                                       M[mapColIndex(j - 1, i, (bound - diff) / 2)] + 1 // insertion
                               });

            } else {

                tmp = std::min({
                                       match + (shorter[i - 1] != longer[j - 1]), // (mis)match
                                       M[mapColIndex(j, i - 1, (bound - diff) / 2)] + 1, // deletion
                                       M[mapColIndex(j - 1, i, (bound - diff) / 2)] + 1 // insertion
                               });

            }

            match = M[mapColIndex(j, i - 1, (bound - diff) / 2)];
            M[mapColIndex(j, i, (bound - diff) / 2)] = tmp;

            early &= ((M[mapColIndex(j, i, (bound - diff) / 2)] + ((diff + i >= j) * (diff + i - j) + (diff + i < j) * (j - diff - i))) > bound); // improved e.t.

        }

        if (early) { // computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }

    }

    return (M[mapColIndex(longer.size(), shorter.size(), (bound - diff) / 2)] > bound) ? (bound + 1) : M[mapColIndex(longer.size(), shorter.size(), (bound - diff) / 2)];

}


void Verification::verify(const AmpliconCollection& ac, Matches& mat, Buffer<Candidate>& buf, lenSeqs_t width, lenSeqs_t t) {

    Candidate c;
    Buffer<Candidate> localBuffer;
    lenSeqs_t M[width]; // reusable DP-matrix (wide enough for all possible calculations for this AmpliconCollection)

    while (!buf.isClosed() || buf.syncSize() > 0) {

        buf.syncSwapContents(localBuffer);

        while (localBuffer.size() > 0) {

            c = localBuffer.pop();

            if (!mat.contains(c.first, c.second)) {

                lenSeqs_t d = computeLengthAwareRow(ac[c.first].seq, ac[c.second].seq, t, M);

                if (d <= t) {
                    mat.add(c.first, c.second, d);
                }

            }

        }

    }

}

}