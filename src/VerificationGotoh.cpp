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

#include <algorithm>
#include <cstring>
#include <iostream>
#include <sstream>

#include "../include/VerificationGotoh.hpp"


namespace GeFaST {

val_t Verification::computeGotohScoreFull(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring) {

    val_t D[lenS + 1][lenT + 1];
    val_t P[lenS + 1][lenT + 1];
    val_t Q[lenS + 1][lenT + 1];

    // initialisation
    D[0][0] = scoring.penOpen;
    for (lenSeqs_t i = 1; i <= lenS; i++) {

        D[i][0] = D[i - 1][0] + scoring.penExtend;
        Q[i][0] = POS_INF;

    }
    for (lenSeqs_t j = 1; j <= lenT; j++) {

        D[0][j] = D[0][j - 1] + scoring.penExtend;
        P[0][j] = POS_INF;

    }
    D[0][0] = 0;

    // fill matrix
    for (lenSeqs_t i = 1; i <= lenS; i++) {

        for (lenSeqs_t j = 1; j <= lenT; j++) {

            P[i][j] = std::min({
                                       D[i - 1][j] + scoring.penOpen + scoring.penExtend,
                                       P[i - 1][j] + scoring.penExtend
                               });

            Q[i][j] = std::min({
                                       D[i][j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[i][j - 1] + scoring.penExtend
                               });

            D[i][j] = std::min({
                                       D[i - 1][j - 1] + (s[i - 1] != t[j - 1]) * scoring.penMismatch,
                                       P[i][j],
                                       Q[i][j]
                               });

        }

    }

    return D[lenS][lenT];

}

std::pair<val_t, std::string> Verification::computeGotohAlignmentFull(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring) {

    val_t D[lenS + 1][lenT + 1];
    val_t P[lenS + 1][lenT + 1];
    val_t Q[lenS + 1][lenT + 1];
    char BT[lenS + 1][lenT + 1];

    // initialisation
    D[0][0] = scoring.penOpen;
    for (lenSeqs_t i = 1; i <= lenS; i++) {

        D[i][0] = D[i - 1][0] + scoring.penExtend;
        BT[i][0] = UP_IN_P;
        Q[i][0] = POS_INF;

    }
    for (lenSeqs_t j = 1; j <= lenT; j++) {

        D[0][j] = D[0][j - 1] + scoring.penExtend;
        BT[0][j] = LEFT_IN_Q;
        P[0][j] = POS_INF;

    }
    D[0][0] = 0;
    BT[1][0] = UP_TO_D;
    BT[0][1] = LEFT_TO_D;

    // fill matrix
    for (lenSeqs_t i = 1; i <= lenS; i++) {

        for (lenSeqs_t j = 1; j <= lenT; j++) {

            P[i][j] = std::min({
                                       D[i - 1][j] + scoring.penOpen + scoring.penExtend,
                                       P[i - 1][j] + scoring.penExtend
                               });

            Q[i][j] = std::min({
                                       D[i][j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[i][j - 1] + scoring.penExtend
                               });

            D[i][j] = std::min({
                                       D[i - 1][j - 1] + (s[i - 1] != t[j - 1]) * scoring.penMismatch,
                                       P[i][j],
                                       Q[i][j]
                               });

            BT[i][j] = (D[i - 1][j - 1] + (s[i - 1] != t[j - 1]) * scoring.penMismatch == D[i][j])
                       | ((Q[i][j] == D[i][j]) << 1)
                       | ((P[i][j] == D[i][j]) << 2)
                       | (((D[i - 1][j] + scoring.penOpen + scoring.penExtend) == P[i][j]) << 3)
                       | (((P[i - 1][j] + scoring.penExtend) == P[i][j]) << 4)
                       | (((D[i][j - 1] + scoring.penOpen + scoring.penExtend) == Q[i][j]) << 5)
                       | (((Q[i][j - 1] + scoring.penExtend) == Q[i][j]) << 6);

        }

    }

    // backtracking (priorities: left > up > diagonal)
    std::string path = "";
    lenSeqs_t i = lenS;
    lenSeqs_t j = lenT;

    long MATRIX = 0;

    while (i != 0 && j != 0) {

        switch (MATRIX) {

            case IN_D: {

                if (BT[i][j] & JUMP_TO_Q) {

                    MATRIX = IN_Q;

                } else if (BT[i][j] & JUMP_TO_P) {

                    MATRIX = IN_P;

                } else if (BT[i][j] & DIAGONAL_IN_D) {

                    path.push_back('\\');
                    i--;
                    j--;

                }

                break;

            }

            case IN_P: {

                /*if (BT[i][j] & UP_IN_P) {

                } else*/ if (BT[i][j] & UP_TO_D) {

                    MATRIX = IN_D;

                }

                path.push_back('|');
                i--;

                break;

            }

            case IN_Q: {

                /*if (BT[i][j] & LEFT_IN_Q) {

                } else*/ if (BT[i][j] & LEFT_TO_D) {

                    MATRIX = IN_D;
                }

                path.push_back('-');
                j--;

                break;

            }

            default: {
                // do nothing (should not occur)
            }

        }

    }

    std::reverse(path.begin(), path.end());

    return std::make_pair(D[lenS][lenT], path);

}


lenSeqs_t Verification::computeGotohFull(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring) {

    val_t D[lenS + 1][lenT + 1];
    val_t P[lenS + 1][lenT + 1];
    val_t Q[lenS + 1][lenT + 1];
    char BT[lenS + 1][lenT + 1];
    lenSeqs_t cntDiffs[lenS + 1][lenT + 1];
    lenSeqs_t cntDiffsP[lenS + 1][lenT + 1];
    lenSeqs_t cntDiffsQ[lenS + 1][lenT + 1];

    // initialisation
    D[0][0] = scoring.penOpen;
    for (lenSeqs_t i = 1; i <= lenS; i++) {

        D[i][0] = D[i - 1][0] + scoring.penExtend;
        BT[i][0] = UP_IN_P;
        Q[i][0] = POS_INF;
        cntDiffs[i][0] = i;

    }
    for (lenSeqs_t j = 1; j <= lenT; j++) {

        D[0][j] = D[0][j - 1] + scoring.penExtend;
        BT[0][j] = LEFT_IN_Q;
        P[0][j] = POS_INF;
        cntDiffs[0][j] = j;

    }
    D[0][0] = 0;
    BT[1][0] = UP_TO_D;
    BT[0][1] = LEFT_TO_D;



    // fill matrix
    for (lenSeqs_t i = 1; i <= lenS; i++) {

        for (lenSeqs_t j = 1; j <= lenT; j++) {

            P[i][j] = std::min({
                                       D[i - 1][j] + scoring.penOpen + scoring.penExtend,
                                       P[i - 1][j] + scoring.penExtend
                               });

            Q[i][j] = std::min({
                                       D[i][j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[i][j - 1] + scoring.penExtend
                               });

            D[i][j] = std::min({
                                       D[i - 1][j - 1] + (s[i - 1] != t[j - 1]) * scoring.penMismatch,
                                       P[i][j],
                                       Q[i][j]
                               });

            BT[i][j] = (D[i - 1][j - 1] + (s[i - 1] != t[j - 1]) * scoring.penMismatch == D[i][j])
                       | ((Q[i][j] == D[i][j]) << 1)
                       | ((P[i][j] == D[i][j]) << 2)
                       | (((D[i - 1][j] + scoring.penOpen + scoring.penExtend) == P[i][j]) << 3)
                       | (((P[i - 1][j] + scoring.penExtend) == P[i][j]) << 4)
                       | (((D[i][j - 1] + scoring.penOpen + scoring.penExtend) == Q[i][j]) << 5)
                       | (((Q[i][j - 1] + scoring.penExtend) == Q[i][j]) << 6);


            if (BT[i][j] & UP_TO_D) {
                cntDiffsP[i][j] = cntDiffs[i - 1][j] + 1;
            } else {
                cntDiffsP[i][j] = cntDiffsP[i - 1][j] + 1;
            }

            if (BT[i][j] & LEFT_TO_D) {
                cntDiffsQ[i][j] = cntDiffs[i][j - 1] + 1;
            } else {
                cntDiffsQ[i][j] = cntDiffsQ[i][j - 1] + 1;
            }

            if (BT[i][j] & JUMP_TO_Q) {
                cntDiffs[i][j] = cntDiffsQ[i][j];
            } else if (BT[i][j] & JUMP_TO_P) {
                cntDiffs[i][j] = cntDiffsP[i][j];
            } else if (BT[i][j] & DIAGONAL_IN_D) {
                cntDiffs[i][j] = cntDiffs[i - 1][j - 1] + (s[i - 1] != t[j - 1]);
            }

        }

    }

#if 0
    std::cout << s << " vs. " << t << std::endl;
    std::cout << "===========D=============" << std::endl;
    for (auto i = 0; i <= lenS; i++) {
        for (auto j = 0; j <= lenT; j++) {
            std::cout << D[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "========================" << std::endl;
    std::cout << "===========P=============" << std::endl;
    for (auto i = 0; i <= lenS; i++) {
        for (auto j = 0; j <= lenT; j++) {
            std::cout << P[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "========================" << std::endl;
    std::cout << "===========Q=============" << std::endl;
    for (auto i = 0; i <= lenS; i++) {
        for (auto j = 0; j <= lenT; j++) {
            std::cout << Q[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "========================" << std::endl;
    std::cout << "===========BT=============" << std::endl;
    for (auto i = 0; i <= lenS; i++) {
        for (auto j = 0; j <= lenT; j++) {
            std::cout << (BT[i][j] & JUMP_TO_Q) << (BT[i][j] & JUMP_TO_P) << (BT[i][j] & DIAGONAL_IN_D) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "========================" << std::endl;
    std::cout << "=============cntDiffs===========" << std::endl;
    for (auto i = 0; i <= lenS; i++) {
        for (auto j = 0; j <= lenT; j++) {
            std::cout << cntDiffs[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "========================" << std::endl;
    std::cout << "=============cntDiffsP===========" << std::endl;
    for (auto i = 0; i <= lenS; i++) {
        for (auto j = 0; j <= lenT; j++) {
            std::cout << cntDiffsP[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "========================" << std::endl;
    std::cout << "=============cntDiffsQ===========" << std::endl;
    for (auto i = 0; i <= lenS; i++) {
        for (auto j = 0; j <= lenT; j++) {
            std::cout << cntDiffsQ[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "========================" << std::endl;
#endif

    return cntDiffs[lenS][lenT];

}

lenSeqs_t Verification::computeGotohRow(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring, val_t* D, val_t* P, val_t* Q, char* BT, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP, lenSeqs_t* cntDiffsQ) {

    // initialise first row
    D[0] = scoring.penOpen;
    for (lenSeqs_t j = 1; j <= lenT; j++) {

        D[j] = D[j - 1] + scoring.penExtend;
        BT[j] = LEFT_IN_Q;
        P[j] = POS_INF;
        cntDiffs[j] = j;

    }
    D[0] = 0;
    BT[1] = LEFT_TO_D;

    // compute remaining rows
    val_t match, tmpD, tmpP;
    lenSeqs_t diff, tmp;

    for (lenSeqs_t i = 1; i <= lenS; i++) {

        // handle left end
        D[0] = scoring.penOpen + i * scoring.penExtend;
        match = (i > 1) * (D[0] - scoring.penExtend);
        Q[0] = POS_INF;
        BT[0] = (i == 1) ? UP_TO_D : UP_IN_P;
        diff = i - 1;
        cntDiffs[0] = i;

        // fill remaining row
        for (lenSeqs_t j = 1; j <= lenT; j++) {

            // arrays D, P, Q, BT
            tmpP = std::min({
                                       D[j] + scoring.penOpen + scoring.penExtend,
                                       P[j] + scoring.penExtend
                            });

            Q[j] = std::min({
                                       D[j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[j - 1] + scoring.penExtend
                                   });

            tmpD = std::min({
                                       match + (s[i - 1] != t[j - 1]) * scoring.penMismatch,
                                       tmpP,
                                       Q[j]
                           });

            BT[j] = (match + (s[i - 1] != t[j - 1]) * scoring.penMismatch == tmpD)
                           | ((Q[j] == tmpD) << 1)
                           | ((tmpP == tmpD) << 2)
                           | (((D[j] + scoring.penOpen + scoring.penExtend) == tmpP) << 3)
                           | (((P[j] + scoring.penExtend) == tmpP) << 4)
                           | (((D[j - 1] + scoring.penOpen + scoring.penExtend) == Q[j]) << 5)
                           | (((Q[j - 1] + scoring.penExtend) == Q[j]) << 6);

            P[j] = tmpP;

            match = D[j];
            D[j] = tmpD;


            // arrays cntDiffs, cntDiffsP, cntDiffsQ
            if (BT[j] & UP_TO_D) {
                cntDiffsP[j] = cntDiffs[j] + 1;
            } else {
                cntDiffsP[j] = cntDiffsP[j] + 1;
            }

            if (BT[j] & LEFT_TO_D) {
                cntDiffsQ[j] = cntDiffs[j - 1] + 1;
            } else {
                cntDiffsQ[j] = cntDiffsQ[j - 1] + 1;
            }

            if (BT[j] & JUMP_TO_Q) {

                diff = cntDiffs[j];
                cntDiffs[j] = cntDiffsQ[j];

            } else if (BT[j] & JUMP_TO_P) {

                diff = cntDiffs[j];
                cntDiffs[j] = cntDiffsP[j];

            } else if (BT[j] & DIAGONAL_IN_D) {

                tmp = diff + (s[i - 1] != t[j - 1]);
                diff = cntDiffs[j];
                cntDiffs[j] = tmp;

            }

        }

    }

    return cntDiffs[lenT];

}

lenSeqs_t Verification::computeGotohEarlyRow(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const lenSeqs_t bound, const Scoring& scoring, val_t* D, val_t* P, val_t* Q, char* BT, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP, lenSeqs_t* cntDiffsQ) {

    // initialise first row
    D[0] = scoring.penOpen;
    for (lenSeqs_t j = 1; j <= lenT; j++) {

        D[j] = D[j - 1] + scoring.penExtend;
        BT[j] = LEFT_IN_Q;
        P[j] = POS_INF;
        cntDiffs[j] = j;

    }
    D[0] = 0;
    BT[1] = LEFT_TO_D;

    // compute remaining rows
    val_t match, tmpD, tmpP;
    lenSeqs_t diff, tmp;
    bool early;

    for (lenSeqs_t i = 1; i <= lenS; i++) {

        // handle left end
        D[0] = scoring.penOpen + i * scoring.penExtend;
        match = (i > 1) * (D[0] - scoring.penExtend);
        Q[0] = POS_INF;
        BT[0] = (i == 1) ? UP_TO_D : UP_IN_P;
        diff = i - 1;
        cntDiffs[0] = i;

        early = true; // early termination flag

        // fill remaining row
        for (lenSeqs_t j = 1; j <= lenT; j++) {

            // arrays D, P, Q, BT
            tmpP = std::min({
                                       D[j] + scoring.penOpen + scoring.penExtend,
                                       P[j] + scoring.penExtend
                            });

            Q[j] = std::min({
                                       D[j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[j - 1] + scoring.penExtend
                                   });

            tmpD = std::min({
                                       match + (s[i - 1] != t[j - 1]) * scoring.penMismatch,
                                       tmpP,
                                       Q[j]
                           });

            BT[j] = (match + (s[i - 1] != t[j - 1]) * scoring.penMismatch == tmpD)
                           | ((Q[j] == tmpD) << 1)
                           | ((tmpP == tmpD) << 2)
                           | (((D[j] + scoring.penOpen + scoring.penExtend) == tmpP) << 3)
                           | (((P[j] + scoring.penExtend) == tmpP) << 4)
                           | (((D[j - 1] + scoring.penOpen + scoring.penExtend) == Q[j]) << 5)
                           | (((Q[j - 1] + scoring.penExtend) == Q[j]) << 6);

            P[j] = tmpP;

            match = D[j];
            D[j] = tmpD;


            // arrays cntDiffs, cntDiffsP, cntDiffsQ
            if (BT[j] & UP_TO_D) {
                cntDiffsP[j] = cntDiffs[j] + 1;
            } else {
                cntDiffsP[j] = cntDiffsP[j] + 1;
            }

            if (BT[j] & LEFT_TO_D) {
                cntDiffsQ[j] = cntDiffs[j - 1] + 1;
            } else {
                cntDiffsQ[j] = cntDiffsQ[j - 1] + 1;
            }

            if (BT[j] & JUMP_TO_Q) {

                diff = cntDiffs[j];
                cntDiffs[j] = cntDiffsQ[j];

            } else if (BT[j] & JUMP_TO_P) {

                diff = cntDiffs[j];
                cntDiffs[j] = cntDiffsP[j];

            } else if (BT[j] & DIAGONAL_IN_D) {

                tmp = diff + (s[i - 1] != t[j - 1]);
                diff = cntDiffs[j];
                cntDiffs[j] = tmp;

            }

            early &= (cntDiffs[j] > bound);

        }

        if (early) {// computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }

    }

    return (cntDiffs[lenT] > bound) ? (bound + 1) : cntDiffs[lenT];

}


lenSeqs_t Verification::computeGotohLengthAwareEarlyRow(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const lenSeqs_t bound, const Scoring& scoring, val_t* D, val_t* P, lenSeqs_t* cntDiffs, lenSeqs_t* cntDiffsP) {

    // long computation not necessary if lengths differ too much
    if (((lenS > lenT) ? (lenS - lenT) : (lenT - lenS)) > bound) {
        return bound + 1;
    }

    if (bound == 0) {
        return lenSeqs_t(s != t);
    }

    const char* shorter = (lenS < lenT) ? s : t;
    lenSeqs_t lenShorter = std::min(lenS, lenT);
    const char* longer = (lenS >= lenT) ? s : t;
    lenSeqs_t lenLonger = std::max(lenS, lenT);
    lenSeqs_t delta = lenLonger - lenShorter;

    // (mis)match is only possibility when we have to consider only one diagonal [happens only if (a) bound = delta = 0, or (b) bound = 1 and delta = 0, but (a) is already covered above]
    if ((bound - delta) / 2 == 0 && (bound + delta) / 2 == 0) {

        lenSeqs_t diffs = 0;
        for (auto i = 0; diffs <= bound && i < lenShorter; i++) {
            diffs += (shorter[i] != longer[i]);
        }

        return diffs;

    }

    // initialise necessary section of first row
    D[0] = scoring.penOpen;
    for (lenSeqs_t j = 1; j <= (bound + delta) / 2 && j <= lenLonger; j++) {

        D[j] = D[j - 1] + scoring.penExtend;
        P[j] = POS_INF;
        cntDiffs[j] = j;

    }
    D[0] = 0;

    // compute remaining rows
    val_t match, minVal, valQ, fromD, fromPQ;
    lenSeqs_t j, diff, minValDiff, diffsQ;
    bool early;

    for (lenSeqs_t i = 1; i <= lenShorter; i++) {

        // handle left end
        D[0] = scoring.penOpen + i * scoring.penExtend;
        match = ((1 < i) && (i <= (bound - delta) / 2 + 1)) * (D[0] - scoring.penExtend) + (i > (bound - delta) / 2 + 1) * D[i - (bound - delta) / 2 - 1];
        valQ = POS_INF;
        diff = (i <= (bound - delta) / 2 + 1) * (i - 1) + (i > (bound - delta) / 2 + 1) * cntDiffs[i - (bound - delta) / 2 - 1];
        cntDiffs[0] = i;

        early = true; // early termination flag

        // fill remaining row
        j = 1 + (i > (bound - delta) / 2) * (i - (bound - delta) / 2 - 1); // same as starting from j = max(1, i - (bound - delta) / 2) with signed integers
        D[j - 1] = POS_INF;
        if (i + (bound + delta) / 2 <= lenLonger) {
            D[i + (bound + delta) / 2] = P[i + (bound + delta) / 2] = POS_INF;
        }

        for (; j <= i + (bound + delta) / 2 && j <= lenLonger; j++) {

            // arrays P & cntDiffsP
            fromD = D[j] + scoring.penOpen + scoring.penExtend;
            fromPQ = P[j] + scoring.penExtend;

            if (fromD <= fromPQ) {

                P[j] = fromD;
                cntDiffsP[j] = cntDiffs[j] + 1;

            } else {

                P[j] = fromPQ;
                cntDiffsP[j]++;

            }

            // arrays Q & cntDiffsQ
            fromD = D[j - 1] + scoring.penOpen + scoring.penExtend;
            fromPQ = valQ + scoring.penExtend;
            if (fromD <= fromPQ) {

                valQ = fromD;
                diffsQ = cntDiffs[j - 1] + 1;

            } else {

                valQ = fromPQ;
                diffsQ++;

            }

            // arrays D & cntDiffs
            minVal = (match + (shorter[i - 1] != longer[j - 1]) * scoring.penMismatch);
            minValDiff = diff + (shorter[i - 1] != longer[j - 1]);
            if (P[j] < minVal) {

                minVal = P[j];
                minValDiff = cntDiffsP[j];

            }
            if (valQ <= minVal){

                minVal = valQ;
                minValDiff = diffsQ;

            }

            match = D[j];
            D[j] = minVal;

            diff = cntDiffs[j];
            cntDiffs[j] = minValDiff;

            early &= ((cntDiffs[j] + llabs((long long)delta + (long long)i - (long long)j)) > bound); // improved e.t.

        }

        if (early) {// computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }

    }

    return (cntDiffs[lenLonger] > bound) ? (bound + 1) : cntDiffs[lenLonger];

}

void Verification::verifyGotoh(const AmpliconCollection& ac, Matches& mat, Buffer<Candidate>& buf, lenSeqs_t width, lenSeqs_t t, const Scoring& scoring) {

    Candidate c;
    Buffer<Candidate> localBuffer;
    val_t D[width];
    val_t P[width];
    lenSeqs_t cntDiffs[width];
    lenSeqs_t cntDiffsP[width];

    while (!buf.isClosed() || buf.syncSize() > 0) {

        buf.syncSwapContents(localBuffer);

        while (localBuffer.size() > 0) {

            c = localBuffer.pop();

            if (!mat.contains(c.first, c.second)) {

                lenSeqs_t d = computeGotohLengthAwareEarlyRow(ac[c.first].seq, ac[c.first].len, ac[c.second].seq, ac[c.second].len, t, scoring, D, P, cntDiffs, cntDiffsP);

                if (d <= t) {
                    mat.add(c.first, c.second, d);
                }

            }

        }

    }

}



Verification::AlignmentInformation Verification::computeGotohCigarFull(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring) {

    val_t D[lenS + 1][lenT + 1];
    val_t P[lenS + 1][lenT + 1];
    val_t Q[lenS + 1][lenT + 1];
    char BT[lenS + 1][lenT + 1];

    // initialisation
    D[0][0] = scoring.penOpen;
    for (lenSeqs_t i = 1; i <= lenS; i++) {

        D[i][0] = D[i - 1][0] + scoring.penExtend;
        BT[i][0] = UP_IN_P;
        Q[i][0] = POS_INF;

    }
    for (lenSeqs_t j = 1; j <= lenT; j++) {

        D[0][j] = D[0][j - 1] + scoring.penExtend;
        BT[0][j] = LEFT_IN_Q;
        P[0][j] = POS_INF;

    }
    D[0][0] = 0;
    BT[1][0] = UP_TO_D;
    BT[0][1] = LEFT_TO_D;

    // fill matrix
    val_t fromD, fromPQ, minVal;
    char tmp;

    for (lenSeqs_t i = 1; i <= lenS; i++) {

        for (lenSeqs_t j = 1; j <= lenT; j++) {

            // array P
            fromD = D[i - 1][j] + scoring.penOpen + scoring.penExtend;
            fromPQ = P[i - 1][j] + scoring.penExtend;

            if (fromD <= fromPQ) {

                P[i][j] = fromD;
                tmp = UP_TO_D;

            } else {

                P[i][j] = fromPQ;
                tmp = UP_IN_P;

            }

            // array Q
            fromD = D[i][j - 1] + scoring.penOpen + scoring.penExtend;
            fromPQ = Q[i][j - 1] + scoring.penExtend;

            if (fromD <= fromPQ) {

                Q[i][j] = fromD;
                tmp |= LEFT_TO_D;

            } else {

                Q[i][j] = fromPQ;
                tmp |= LEFT_IN_Q;

            }

            // arrays D & BT
            minVal = D[i - 1][j - 1] + (s[i - 1] != t[j - 1]) * scoring.penMismatch;
            BT[i][j] = DIAGONAL_IN_D;
            if (P[i][j] < minVal) {

                minVal = P[i][j];
                BT[i][j] = JUMP_TO_P;

            }
            if (Q[i][j] <= minVal){

                minVal = Q[i][j];
                BT[i][j] = JUMP_TO_Q;

            }

            BT[i][j] |= tmp;

            D[i][j] = minVal;

        }

    }

    // backtracking (priorities: left > diagonal > up)
    std::vector<std::pair<char, lenSeqs_t >> cigarSegs = {std::make_pair('N',0)};
    lenSeqs_t len = 0;
    lenSeqs_t numDiffs = 0;
    lenSeqs_t i = lenS;
    lenSeqs_t j = lenT;

    while (i != 0 && j != 0) {

        if ((cigarSegs.back().first == 'I') && (BT[i][j + 1] & LEFT_IN_Q)) {

            len++;
            numDiffs++;
            cigarSegs.back().second++;

            j--;

        } else if ((cigarSegs.back().first == 'D') && (BT[i + 1][j] & UP_IN_P)) {

            len++;
            numDiffs++;
            cigarSegs.back().second++;

            i--;

        } else if (BT[i][j] & JUMP_TO_Q) {

            len++;
            numDiffs++;
            if (cigarSegs.back().first == 'I') {
                cigarSegs.back().second++;
            } else {
                cigarSegs.push_back(std::make_pair('I', 1));
            }

            j--;

        } else if (BT[i][j] & DIAGONAL_IN_D) {

            len++;
            numDiffs += (s[i - 1] != t[j - 1]);
            if (cigarSegs.back().first == 'M') {
                cigarSegs.back().second++;
            } else {
                cigarSegs.push_back(std::make_pair('M', 1));
            }

            i--;
            j--;

        } else { // BT[i][j] & JUMP_TO_P

            len++;
            numDiffs++;
            if (cigarSegs.back().first == 'D') {
                cigarSegs.back().second++;
            } else {
                cigarSegs.push_back(std::make_pair('D', 1));
            }

            i--;

        }

    }

    while (i > 0) {

        len++;
        numDiffs++;
        if (cigarSegs.back().first == 'D') {
            cigarSegs.back().second++;
        } else {
            cigarSegs.push_back(std::make_pair('D', 1));
        }

        i--;

    }

    while (j > 0) {

        len++;
        numDiffs++;
        if (cigarSegs.back().first == 'I') {
            cigarSegs.back().second++;
        } else {
            cigarSegs.push_back(std::make_pair('I', 1));
        }

        j--;

    }

    std::stringstream sStream;
    for (auto p = cigarSegs.size() - 1; p > 0; p--) {

        if (cigarSegs[p].second > 1) sStream << cigarSegs[p].second;
        sStream << cigarSegs[p].first;

    }

    return AlignmentInformation(sStream.str(), len, numDiffs);

}

Verification::AlignmentInformation Verification::computeGotohCigarFull1(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring, val_t* D, val_t* P, val_t* Q, char* BT) {

    lenSeqs_t width = lenT + 1;

    // initialisation
    D[0] = scoring.penOpen;
    for (lenSeqs_t i = 1; i <= lenS; i++) {

        D[i * width] = D[(i - 1) * width] + scoring.penExtend;
        BT[i * width] = UP_IN_P;
        Q[i * width] = POS_INF;

    }
    for (lenSeqs_t j = 1; j <= lenT; j++) {

        D[j] = D[j - 1] + scoring.penExtend;
        BT[j] = LEFT_IN_Q;
        P[j] = POS_INF;

    }
    D[0] = 0;
    BT[1 * width] = UP_TO_D;
    BT[1] = LEFT_TO_D;

    // fill matrix
    val_t fromD, fromPQ, minVal;
    char tmp;

    for (lenSeqs_t i = 1; i <= lenS; i++) {

        for (lenSeqs_t j = 1; j <= lenT; j++) {

            // array P
            fromD = D[(i - 1) * width + j] + scoring.penOpen + scoring.penExtend;
            fromPQ = P[(i - 1) * width + j] + scoring.penExtend;

            if (fromD <= fromPQ) {

                P[i * width + j] = fromD;
                tmp = UP_TO_D;

            } else {

                P[i * width + j] = fromPQ;
                tmp = UP_IN_P;

            }

            // array Q
            fromD = D[i * width + j - 1] + scoring.penOpen + scoring.penExtend;
            fromPQ = Q[i * width + j - 1] + scoring.penExtend;

            if (fromD <= fromPQ) {

                Q[i * width + j] = fromD;
                tmp |= LEFT_TO_D;

            } else {

                Q[i * width + j] = fromPQ;
                tmp |= LEFT_IN_Q;

            }

            // arrays D & BT
            minVal = D[(i - 1) * width + j - 1] + (s[i - 1] != t[j - 1]) * scoring.penMismatch;
            BT[i * width + j] = DIAGONAL_IN_D;
            if (P[i * width + j] < minVal) {

                minVal = P[i * width + j];
                BT[i * width + j] = JUMP_TO_P;

            }
            if (Q[i * width + j] <= minVal){

                minVal = Q[i * width + j];
                BT[i * width + j] = JUMP_TO_Q;

            }

            BT[i * width + j] |= tmp;

            D[i * width + j] = minVal;

        }

    }

    // backtracking (priorities: left > diagonal > up)
    std::vector<std::pair<char, lenSeqs_t >> cigarSegs = {std::make_pair('N',0)};
    lenSeqs_t len = 0;
    lenSeqs_t numDiffs = 0;
    lenSeqs_t i = lenS;
    lenSeqs_t j = lenT;

    while (i != 0 && j != 0) {

        if ((cigarSegs.back().first == 'I') && (BT[i * width  + j + 1] & LEFT_IN_Q)) {

            len++;
            numDiffs++;
            cigarSegs.back().second++;

            j--;

        } else if ((cigarSegs.back().first == 'D') && (BT[(i + 1) * width + j] & UP_IN_P)) {

            len++;
            numDiffs++;
            cigarSegs.back().second++;

            i--;

        } else if (BT[i * width + j] & JUMP_TO_Q) {

            len++;
            numDiffs++;
            if (cigarSegs.back().first == 'I') {
                cigarSegs.back().second++;
            } else {
                cigarSegs.push_back(std::make_pair('I', 1));
            }

            j--;

        } else if (BT[i * width + j] & DIAGONAL_IN_D) {

            len++;
            numDiffs += (s[i - 1] != t[j - 1]);
            if (cigarSegs.back().first == 'M') {
                cigarSegs.back().second++;
            } else {
                cigarSegs.push_back(std::make_pair('M', 1));
            }

            i--;
            j--;

        } else { // BT[i * width + j] & JUMP_TO_P

            len++;
            numDiffs++;
            if (cigarSegs.back().first == 'D') {
                cigarSegs.back().second++;
            } else {
                cigarSegs.push_back(std::make_pair('D', 1));
            }

            i--;

        }

    }

    while (i > 0) {

        len++;
        numDiffs++;
        if (cigarSegs.back().first == 'D') {
            cigarSegs.back().second++;
        } else {
            cigarSegs.push_back(std::make_pair('D', 1));
        }

        i--;

    }

    while (j > 0) {

        len++;
        numDiffs++;
        if (cigarSegs.back().first == 'I') {
            cigarSegs.back().second++;
        } else {
            cigarSegs.push_back(std::make_pair('I', 1));
        }

        j--;

    }

    std::stringstream sStream;
    for (auto p = cigarSegs.size() - 1; p > 0; p--) {

        if (cigarSegs[p].second > 1) sStream << cigarSegs[p].second;
        sStream << cigarSegs[p].first;

    }

    return AlignmentInformation(sStream.str(), len, numDiffs);

}


Verification::AlignmentInformation Verification::computeGotohCigarRow(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring) {

    val_t D[lenT + 1];
    val_t P[lenT + 1];
    char BT[lenS + 1][lenT + 1];

    // initialise first row
    D[0] = scoring.penOpen;
    for (lenSeqs_t j = 1; j <= lenT; j++) {

        D[j] = D[j - 1] + scoring.penExtend;
        P[j] = POS_INF;

    }
    D[0] = 0;

    // compute remaining rows
    val_t match, valQ, fromD, fromPQ, minVal;
    char tmp;

    for (lenSeqs_t i = 1; i <= lenS; i++) {

        // handle left end
        D[0] = scoring.penOpen + i * scoring.penExtend;
        match = (i > 1) * (D[0] - scoring.penExtend);
        valQ = POS_INF;

        // fill remaining row
        for (lenSeqs_t j = 1; j <= lenT; j++) {

            // array P
            fromD = D[j] + scoring.penOpen + scoring.penExtend;
            fromPQ = P[j] + scoring.penExtend;

            if (fromD <= fromPQ) {

                P[j] = fromD;
                tmp = UP_TO_D;

            } else {

                P[j] = fromPQ;
                tmp = UP_IN_P;

            }

            // array Q
            fromD = D[j - 1] + scoring.penOpen + scoring.penExtend;
            fromPQ = valQ + scoring.penExtend;

            if (fromD <= fromPQ) {

                valQ = fromD;
                tmp |= LEFT_TO_D;

            } else {

                valQ = fromPQ;
                tmp |= LEFT_IN_Q;

            }

            // arrays D & BT
            minVal = match + (s[i - 1] != t[j - 1]) * scoring.penMismatch;
            BT[i][j] = DIAGONAL_IN_D;
            if (P[j] < minVal) {

                minVal = P[j];
                BT[i][j] = JUMP_TO_P;

            }
            if (valQ <= minVal){

                minVal = valQ;
                BT[i][j] = JUMP_TO_Q;

            }

            BT[i][j] |= tmp;

            match = D[j];
            D[j] = minVal;

        }

    }

    // backtracking (priorities: left > diagonal > up)
    std::vector<std::pair<char, lenSeqs_t >> cigarSegs = {std::make_pair('N',0)};
    lenSeqs_t len = 0;
    lenSeqs_t numDiffs = 0;
    lenSeqs_t i = lenS;
    lenSeqs_t j = lenT;

    while (i != 0 && j != 0) {

        if ((cigarSegs.back().first == 'I') && (BT[i][j + 1] & LEFT_IN_Q)) {

            len++;
            numDiffs++;
            cigarSegs.back().second++;

            j--;

        } else if ((cigarSegs.back().first == 'D') && (BT[i + 1][j] & UP_IN_P)) {

            len++;
            numDiffs++;
            cigarSegs.back().second++;

            i--;

        } else if (BT[i][j] & JUMP_TO_Q) {

            len++;
            numDiffs++;
            if (cigarSegs.back().first == 'I') {
                cigarSegs.back().second++;
            } else {
                cigarSegs.push_back(std::make_pair('I', 1));
            }

            j--;

        } else if (BT[i][j] & DIAGONAL_IN_D) {

            len++;
            numDiffs += (s[i - 1] != t[j - 1]);
            if (cigarSegs.back().first == 'M') {
                cigarSegs.back().second++;
            } else {
                cigarSegs.push_back(std::make_pair('M', 1));
            }

            i--;
            j--;

        } else { // BT[i][j] & JUMP_TO_P

            len++;
            numDiffs++;
            if (cigarSegs.back().first == 'D') {
                cigarSegs.back().second++;
            } else {
                cigarSegs.push_back(std::make_pair('D', 1));
            }

            i--;

        }

    }

    while (i > 0) {

        len++;
        numDiffs++;
        if (cigarSegs.back().first == 'D') {
            cigarSegs.back().second++;
        } else {
            cigarSegs.push_back(std::make_pair('D', 1));
        }

        i--;

    }

    while (j > 0) {

        len++;
        numDiffs++;
        if (cigarSegs.back().first == 'I') {
            cigarSegs.back().second++;
        } else {
            cigarSegs.push_back(std::make_pair('I', 1));
        }

        j--;

    }

    std::stringstream sStream;
    for (auto p = cigarSegs.size() - 1; p > 0; p--) {

        if (cigarSegs[p].second > 1) sStream << cigarSegs[p].second;
        sStream << cigarSegs[p].first;

    }

    return AlignmentInformation(sStream.str(), len, numDiffs);

}

Verification::AlignmentInformation Verification::computeGotohCigarRow1(const char* s, const lenSeqs_t lenS, const char* t, const lenSeqs_t lenT, const Scoring& scoring, val_t* D, val_t* P, char* BT) {

    lenSeqs_t width = lenT + 1;

    // initialise first row
    D[0] = scoring.penOpen;
    for (lenSeqs_t j = 1; j <= lenT; j++) {

        D[j] = D[j - 1] + scoring.penExtend;
        P[j] = POS_INF;

    }
    D[0] = 0;

    // compute remaining rows
    val_t match, valQ, fromD, fromPQ, minVal;
    char tmp;

    for (lenSeqs_t i = 1; i <= lenS; i++) {

        // handle left end
        D[0] = scoring.penOpen + i * scoring.penExtend;
        match = (i > 1) * (D[0] - scoring.penExtend);
        valQ = POS_INF;

        // fill remaining row
        for (lenSeqs_t j = 1; j <= lenT; j++) {

            // array P
            fromD = D[j] + scoring.penOpen + scoring.penExtend;
            fromPQ = P[j] + scoring.penExtend;

            if (fromD <= fromPQ) {

                P[j] = fromD;
                tmp = UP_TO_D;

            } else {

                P[j] = fromPQ;
                tmp = UP_IN_P;

            }

            // array Q
            fromD = D[j - 1] + scoring.penOpen + scoring.penExtend;
            fromPQ = valQ + scoring.penExtend;

            if (fromD <= fromPQ) {

                valQ = fromD;
                tmp |= LEFT_TO_D;

            } else {

                valQ = fromPQ;
                tmp |= LEFT_IN_Q;

            }

            // arrays D & BT
            minVal = match + (s[i - 1] != t[j - 1]) * scoring.penMismatch;
            BT[i * width + j] = DIAGONAL_IN_D;
            if (P[j] < minVal) {

                minVal = P[j];
                BT[i * width + j] = JUMP_TO_P;

            }
            if (valQ <= minVal){

                minVal = valQ;
                BT[i * width + j] = JUMP_TO_Q;

            }

            BT[i * width + j] |= tmp;

            match = D[j];
            D[j] = minVal;

        }

    }

    // backtracking (priorities: left > diagonal > up)
    std::vector<std::pair<char, lenSeqs_t >> cigarSegs = {std::make_pair('N',0)};
    lenSeqs_t len = 0;
    lenSeqs_t numDiffs = 0;
    lenSeqs_t i = lenS;
    lenSeqs_t j = lenT;

    while (i != 0 && j != 0) {

        if ((cigarSegs.back().first == 'I') && (BT[i * width + j + 1] & LEFT_IN_Q)) {

            len++;
            numDiffs++;
            cigarSegs.back().second++;

            j--;

        } else if ((cigarSegs.back().first == 'D') && (BT[(i + 1) * width + j] & UP_IN_P)) {

            len++;
            numDiffs++;
            cigarSegs.back().second++;

            i--;

        } else if (BT[i * width + j] & JUMP_TO_Q) {

            len++;
            numDiffs++;
            if (cigarSegs.back().first == 'I') {
                cigarSegs.back().second++;
            } else {
                cigarSegs.push_back(std::make_pair('I', 1));
            }

            j--;

        } else if (BT[i * width + j] & DIAGONAL_IN_D) {

            len++;
            numDiffs += (s[i - 1] != t[j - 1]);
            if (cigarSegs.back().first == 'M') {
                cigarSegs.back().second++;
            } else {
                cigarSegs.push_back(std::make_pair('M', 1));
            }

            i--;
            j--;

        } else { // BT[i][j] & JUMP_TO_P

            len++;
            numDiffs++;
            if (cigarSegs.back().first == 'D') {
                cigarSegs.back().second++;
            } else {
                cigarSegs.push_back(std::make_pair('D', 1));
            }

            i--;

        }

    }

    while (i > 0) {

        len++;
        numDiffs++;
        if (cigarSegs.back().first == 'D') {
            cigarSegs.back().second++;
        } else {
            cigarSegs.push_back(std::make_pair('D', 1));
        }

        i--;

    }

    while (j > 0) {

        len++;
        numDiffs++;
        if (cigarSegs.back().first == 'I') {
            cigarSegs.back().second++;
        } else {
            cigarSegs.push_back(std::make_pair('I', 1));
        }

        j--;

    }

    std::stringstream sStream;
    for (auto p = cigarSegs.size() - 1; p > 0; p--) {

        if (cigarSegs[p].second > 1) sStream << cigarSegs[p].second;
        sStream << cigarSegs[p].first;

    }

    return AlignmentInformation(sStream.str(), len, numDiffs);

}

}