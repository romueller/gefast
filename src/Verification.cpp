/*
 * SCT-PJ
 *
 * Copyright (C) 2016 Robert Mueller
 *
 * TODO add licence text (e.g. GNU (Affero) General Public License)
 *
 * Contact: Robert Mueller <romueller@techfak.uni-bielefeld.de>
 * Faculty of Technology, Bielefeld University,
 * PO box 100131, DE-33501 Bielefeld, Germany
 */

#include "../include/Verification.hpp"


namespace SCT_PJ {

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


lenSeqs_t Verification::computeClassicRow(const std::string& s, const std::string& t) {

    lenSeqs_t M[t.size() + 1];
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

lenSeqs_t Verification::computeBoundedRow(const std::string& s, const std::string& t, const lenSeqs_t bound) {

    // long computation not necessary if lengths differ too much
    if (((s.size() > t.size()) ? (s.size() - t.size()) : (t.size() - s.size())) > bound) {
        return bound + 1;
    }

    if (bound == 0) {
        return lenSeqs_t(s != t);
    }


    lenSeqs_t M[t.size() + 1];
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

lenSeqs_t Verification::computeBoundedRowSlim(const std::string& s, const std::string& t, const lenSeqs_t bound) {

    // long computation not necessary if lengths differ too much
    if (((s.size() > t.size()) ? (s.size() - t.size()) : (t.size() - s.size())) > bound) {
        return bound + 1;
    }

    if (bound == 0) {
        return lenSeqs_t(s != t);
    }


    lenSeqs_t M[std::min(t.size() + 1, 2 * bound + 1)];
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

lenSeqs_t Verification::computeLengthAwareRow(const std::string& s, const std::string& t, const lenSeqs_t bound) {

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


    lenSeqs_t M[longer.size() + 1];
    lenSeqs_t match, tmp;

    // initialise necessary sections of first row
    for (lenSeqs_t j = 0; j <= (bound + diff) / 2 && j <= longer.size(); j++) {
        M[j] = j;
    }

    bool early;

    // compute sections of remaining rows
    for (lenSeqs_t i = 1; i <= shorter.size(); i++) {

        early = true; // early termination flag

        match = (i <= (bound - diff) / 2 + 1) * (i - 1) + (i > (bound - diff) / 2 + 1) * M[i - (bound - diff) / 2 - 1]; // same as match = (i <= (bound - diff) / 2 + 1) ? (i - 1) : M[i - (bound - diff) / 2 - 1]
        M[0] = i;

        for (lenSeqs_t j = 1 + (i > (bound - diff) / 2) * (i - (bound - diff) / 2 - 1); j <= i + (bound + diff) / 2 && j <= longer.size(); j++) { // same as starting from j = max(1, i - (bound - diff) / 2) with signed integers

            if ((bound - diff) / 2 == 0 && (bound + diff) / 2 == 0) {

                tmp = match + (shorter[i - 1] != longer[j - 1]); // (mis)match is only possibility since we have to consider only one diagonal

            } else if (j == ((i > (bound - diff) / 2) * (i - (bound - diff) / 2))) {

                tmp = std::min({
                                       match + (shorter[i - 1] != longer[j - 1]), // (mis)match
                                       M[j] + 1 // deletion
                               });

            } else if (j == (i + (bound + diff) / 2)) {

                tmp = std::min({
                                       match + (shorter[i - 1] != longer[j - 1]), // (mis)match
                                       M[j - 1] + 1 // insertion
                               });

            } else {

                tmp = std::min({
                                       match + (shorter[i - 1] != longer[j - 1]), // (mis)match
                                       M[j] + 1, // deletion
                                       M[j - 1] + 1 // insertion
                               });

            }

            match = M[j];
            M[j] = tmp;

            early &= ((M[j] + ((diff + i >= j) * (diff + i - j) + (diff + i < j) * (j - diff - i))) > bound); // improved e.t.

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

lenSeqs_t Verification::computeLengthAwareRowSlim(const std::string& s, const std::string& t, const lenSeqs_t bound) {

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


    lenSeqs_t M[std::min(longer.size() + 1, (bound - diff) / 2 + 1 + (bound + diff) / 2)];
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


void Verification::verify(const AmpliconCollection& ac, Matches& mat, Buffer<Candidate>& buf, lenSeqs_t t) {

    numSeqs_t matches = 0; //TODO remove
    Candidate c;

    while (!buf.getFlag() || buf.size() > 0) {

        if (buf.size() > 0) {

            c = buf.pop();

            if (!mat.contains(c.first, c.second)) {

                lenSeqs_t d = computeBoundedRow(ac[c.first].seq, ac[c.second].seq, t); //TODO choose "best" implementation

                if (d <= t) { matches++; //TODO remove
                    mat.add(c.first, c.second, d);
                }

            }

        }

    }

    std::cout << "#matches = " << matches << std::endl; //TODO remove

}

}