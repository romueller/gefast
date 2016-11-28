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

#include <algorithm>
#include <cstring>
#include <iostream>
#include <sstream>

#include "../include/VerificationGotoh.hpp"


namespace SCT_PJ {

Verification::val_t Verification::computeGotohScoreFull(const std::string &s, const std::string &t, const Scoring& scoring) {

    ValArray D(s.size() + 1, std::vector<val_t>(t.size() + 1));
    ValArray P(s.size() + 1, std::vector<val_t>(t.size() + 1));
    ValArray Q(s.size() + 1, std::vector<val_t>(t.size() + 1));

    // initialisation
    D[0][0] = scoring.penOpen;
    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        D[i][0] = D[i - 1][0] + scoring.penExtend;
        Q[i][0] = NEG_INF;

    }
    for (lenSeqs_t j = 1; j <= t.size(); j++) {

        D[0][j] = D[0][j - 1] + scoring.penExtend;
        P[0][j] = NEG_INF;

    }
    D[0][0] = 0;

    // fill matrix
    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        for (lenSeqs_t j = 1; j <= t.size(); j++) {

            P[i][j] = std::max({
                                       D[i - 1][j] + scoring.penOpen + scoring.penExtend,
                                       P[i - 1][j] + scoring.penExtend
                               });

            Q[i][j] = std::max({
                                       D[i][j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[i][j - 1] + scoring.penExtend
                               });

            D[i][j] = std::max({
                                       D[i - 1][j - 1] + ((s[i - 1] == t[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       P[i][j],
                                       Q[i][j]
                               });

        }

    }

    return D[s.size()][t.size()];

}

std::pair<Verification::val_t, std::string> Verification::computeGotohAlignmentFull(const std::string &s, const std::string &t, const Scoring& scoring) {

    ValArray D(s.size() + 1, std::vector<val_t>(t.size() + 1));
    ValArray P(s.size() + 1, std::vector<val_t>(t.size() + 1));
    ValArray Q(s.size() + 1, std::vector<val_t>(t.size() + 1));
    BacktrackingArray BT(s.size() + 1, std::vector<char>(t.size() + 1));

    // initialisation
    D[0][0] = scoring.penOpen;
    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        D[i][0] = D[i - 1][0] + scoring.penExtend;
        BT[i][0] = UP_IN_P;
        Q[i][0] = NEG_INF;

    }
    for (lenSeqs_t j = 1; j <= t.size(); j++) {

        D[0][j] = D[0][j - 1] + scoring.penExtend;
        BT[0][j] = LEFT_IN_Q;
        P[0][j] = NEG_INF;

    }
    D[0][0] = 0;
    BT[1][0] = UP_TO_D;
    BT[0][1] = LEFT_TO_D;

    // fill matrix
    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        for (lenSeqs_t j = 1; j <= t.size(); j++) {

            P[i][j] = std::max({
                                       D[i - 1][j] + scoring.penOpen + scoring.penExtend,
                                       P[i - 1][j] + scoring.penExtend
                               });

            Q[i][j] = std::max({
                                       D[i][j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[i][j - 1] + scoring.penExtend
                               });

            D[i][j] = std::max({
                                       D[i - 1][j - 1] + ((s[i - 1] == t[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       P[i][j],
                                       Q[i][j]
                               });

            BT[i][j] = (D[i - 1][j - 1] + ((s[i - 1] == t[j - 1]) ? scoring.rewardMatch : scoring.penMismatch) == D[i][j])
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
    lenSeqs_t i = s.size();
    lenSeqs_t j = t.size();

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

    return std::make_pair(D[s.size()][t.size()], path);

}


lenSeqs_t Verification::computeGotohFull(const std::string &s, const std::string &t, const Scoring& scoring) {

    ValArray D(s.size() + 1, std::vector<val_t>(t.size() + 1));
    ValArray P(s.size() + 1, std::vector<val_t>(t.size() + 1));
    ValArray Q(s.size() + 1, std::vector<val_t>(t.size() + 1));
    BacktrackingArray BT(s.size() + 1, std::vector<char>(t.size() + 1));
    CountArray cntDiffs(s.size() + 1, std::vector<lenSeqs_t>(t.size() + 1));
    CountArray cntDiffsP(s.size() + 1, std::vector<lenSeqs_t>(t.size() + 1));
    CountArray cntDiffsQ(s.size() + 1, std::vector<lenSeqs_t>(t.size() + 1));

    // initialisation
    D[0][0] = scoring.penOpen;
    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        D[i][0] = D[i - 1][0] + scoring.penExtend;
        BT[i][0] = UP_IN_P;
        Q[i][0] = NEG_INF;
        cntDiffs[i][0] = i;

    }
    for (lenSeqs_t j = 1; j <= t.size(); j++) {

        D[0][j] = D[0][j - 1] + scoring.penExtend;
        BT[0][j] = LEFT_IN_Q;
        P[0][j] = NEG_INF;
        cntDiffs[0][j] = j;

    }
    D[0][0] = 0;
    BT[1][0] = UP_TO_D;
    BT[0][1] = LEFT_TO_D;



    // fill matrix
    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        for (lenSeqs_t j = 1; j <= t.size(); j++) {

            P[i][j] = std::max({
                                       D[i - 1][j] + scoring.penOpen + scoring.penExtend,
                                       P[i - 1][j] + scoring.penExtend
                               });

            Q[i][j] = std::max({
                                       D[i][j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[i][j - 1] + scoring.penExtend
                               });

            D[i][j] = std::max({
                                       D[i - 1][j - 1] + ((s[i - 1] == t[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       P[i][j],
                                       Q[i][j]
                               });

            BT[i][j] = (D[i - 1][j - 1] + ((s[i - 1] == t[j - 1]) ? scoring.rewardMatch : scoring.penMismatch) == D[i][j])
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
    for (auto i = 0; i <= s.size(); i++) {
        for (auto j = 0; j <= t.size(); j++) {
            std::cout << D[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "========================" << std::endl;
    std::cout << "===========P=============" << std::endl;
    for (auto i = 0; i <= s.size(); i++) {
        for (auto j = 0; j <= t.size(); j++) {
            std::cout << P[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "========================" << std::endl;
    std::cout << "===========Q=============" << std::endl;
    for (auto i = 0; i <= s.size(); i++) {
        for (auto j = 0; j <= t.size(); j++) {
            std::cout << Q[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "========================" << std::endl;
    std::cout << "===========BT=============" << std::endl;
    for (auto i = 0; i <= s.size(); i++) {
        for (auto j = 0; j <= t.size(); j++) {
            std::cout << (BT[i][j] & JUMP_TO_Q) << (BT[i][j] & JUMP_TO_P) << (BT[i][j] & DIAGONAL_IN_D) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "========================" << std::endl;
    std::cout << "=============cntDiffs===========" << std::endl;
    for (auto i = 0; i <= s.size(); i++) {
        for (auto j = 0; j <= t.size(); j++) {
            std::cout << cntDiffs[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "========================" << std::endl;
    std::cout << "=============cntDiffsP===========" << std::endl;
    for (auto i = 0; i <= s.size(); i++) {
        for (auto j = 0; j <= t.size(); j++) {
            std::cout << cntDiffsP[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "========================" << std::endl;
    std::cout << "=============cntDiffsQ===========" << std::endl;
    for (auto i = 0; i <= s.size(); i++) {
        for (auto j = 0; j <= t.size(); j++) {
            std::cout << cntDiffsQ[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "========================" << std::endl;
    if (s ==
        "agctccaatagcgaatactaattctgctgcagttaaaaagctcgtagttgaacttttggtaggacgtgcccagccttagagtgaactctctttgtgttctgtgtgtacgtttggccatactctgcctttgcaaaagggtagtctttcactgtaaaaaaattagggtgtttcaagcaaaagatttttggaatatattagtatgggatgataagataggctctgagtgctattttgttggtttgtacatttagagaatgattaacagggacagttgggggtattcatatttacatgtcagaggtgaaattcttgggatttttggaaagatgaacgactgcgaaggcatttaccaaggatgttttca" &&
        t ==
        "agctccaatagcgaatactaattctgctgcagttaaaaagctcgtagttgaacttttggtaggacgtgcccagccttagagtgaactctctttgtgttctgtgtgtacgtttggccatactctgcctttgcaaaagggtagtctttcactgtaaaaaaattagggtgtttcaagcaaaagatttttggaatatattagtatgggatgataagataggctctgagtgctattttgttggtttgtacatttagagaatgattaacagggacagttgggggtattcatatttacatgtcagaggtgaaattcttggattttggaaagatgaacgactgcgaaggcatttaccaaggatgttttca" ||
        t ==
        "agctccaatagcgaatactaattctgctgcagttaaaaagctcgtagttgaacttttggtaggacgtgcccagccttagagtgaactctctttgtgttctgtgtgtacgtttggccatactctgcctttgcaaaagggtagtctttcactgtaaaaaaattagggtgtttcaagcaaaagatttttggaatatattagtatgggatgataagataggctctgagtgctattttgttggtttgtacatttagagaatgattaacagggacagttgggggtattcatatttacatgtcagaggtgaaattcttgggatttttggaaagatgaacgactgcgaaggcatttaccaaggatgttttca" &&
        s ==
        "agctccaatagcgaatactaattctgctgcagttaaaaagctcgtagttgaacttttggtaggacgtgcccagccttagagtgaactctctttgtgttctgtgtgtacgtttggccatactctgcctttgcaaaagggtagtctttcactgtaaaaaaattagggtgtttcaagcaaaagatttttggaatatattagtatgggatgataagataggctctgagtgctattttgttggtttgtacatttagagaatgattaacagggacagttgggggtattcatatttacatgtcagaggtgaaattcttggattttggaaagatgaacgactgcgaaggcatttaccaaggatgttttca")
        std::cout << "GotohFull: s = " << s << "; t = " << t << std::endl << "D = " << D[s.size()][t.size()]
                  << "; P = " << P[s.size()][t.size()] << "; Q = " << Q[s.size()][t.size()] << "; cntDiffs = "
                  << cntDiffs[s.size()][t.size()] << std::endl;
#endif

    return cntDiffs[s.size()][t.size()];

}

lenSeqs_t Verification::computeGotohRow(const std::string &s, const std::string &t, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], char BT[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]) {

    // initialise first row
    D[0] = scoring.penOpen;
    for (lenSeqs_t j = 1; j <= t.size(); j++) {

        D[j] = D[j - 1] + scoring.penExtend;
        BT[j] = LEFT_IN_Q;
        P[j] = NEG_INF;
        cntDiffs[j] = j;

    }
    D[0] = 0;
    BT[1] = LEFT_TO_D;

    // compute remaining rows
    val_t match, tmpD, tmpP;
    lenSeqs_t diff, tmp;

    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        // handle left end
        D[0] = scoring.penOpen + i * scoring.penExtend;
        match = D[0] - scoring.penExtend - (i == 1) * scoring.penOpen;
        Q[0] = NEG_INF;
        BT[0] = (i == 1) ? UP_TO_D : UP_IN_P;
        diff = i - 1;
        cntDiffs[0] = i;

        // fill remaining row
        for (lenSeqs_t j = 1; j <= t.size(); j++) {

            // arrays D, P, Q, BT
            tmpP = std::max({
                                       D[j] + scoring.penOpen + scoring.penExtend,
                                       P[j] + scoring.penExtend
                            });

            Q[j] = std::max({
                                       D[j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[j - 1] + scoring.penExtend
                                   });

            tmpD = std::max({
                                       match + ((s[i - 1] == t[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       tmpP,
                                       Q[j]
                           });

            BT[j] = (match + ((s[i - 1] == t[j - 1]) ? scoring.rewardMatch : scoring.penMismatch) == tmpD)
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

    return cntDiffs[t.size()];

}

lenSeqs_t Verification::computeGotohEarlyRow(const std::string &s, const std::string &t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], char BT[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]) {

    // initialise first row
    D[0] = scoring.penOpen;
    for (lenSeqs_t j = 1; j <= t.size(); j++) {

        D[j] = D[j - 1] + scoring.penExtend;
        BT[j] = LEFT_IN_Q;
        P[j] = NEG_INF;
        cntDiffs[j] = j;

    }
    D[0] = 0;
    BT[1] = LEFT_TO_D;

    // compute remaining rows
    val_t match, tmpD, tmpP;
    lenSeqs_t diff, tmp;
    bool early;

    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        //handle left end
        D[0] = scoring.penOpen + i * scoring.penExtend;
        match = D[0] - scoring.penExtend - (i == 1) * scoring.penOpen;
        Q[0] = NEG_INF;
        BT[0] = (i == 1) ? UP_TO_D : UP_IN_P;
        diff = i - 1;
        cntDiffs[0] = i;

        early = true; // early termination flag

        // fill remaining row
        for (lenSeqs_t j = 1; j <= t.size(); j++) {

            // arrays D, P, Q, BT
            tmpP = std::max({
                                       D[j] + scoring.penOpen + scoring.penExtend,
                                       P[j] + scoring.penExtend
                            });

            Q[j] = std::max({
                                       D[j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[j - 1] + scoring.penExtend
                                   });

            tmpD = std::max({
                                       match + ((s[i - 1] == t[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       tmpP,
                                       Q[j]
                           });

            BT[j] = (match + ((s[i - 1] == t[j - 1]) ? scoring.rewardMatch : scoring.penMismatch) == tmpD)
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

    return (cntDiffs[t.size()] > bound) ? (bound + 1) : cntDiffs[t.size()];

}

lenSeqs_t Verification::computeGotohBoundedEarlyRow(const std::string &s, const std::string &t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], char BT[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]) {

    // long computation not necessary if lengths differ too much
    if (((s.size() > t.size()) ? (s.size() - t.size()) : (t.size() - s.size())) > bound) {
        return bound + 1;
    }

    if (bound == 0) {
        return lenSeqs_t(s != t);
    }

    // initialise necessary section of first row
    D[0] = scoring.penOpen;
    for (lenSeqs_t j = 1; j <= bound && j <= t.size(); j++) {

        D[j] = D[j - 1] + scoring.penExtend;
        BT[j] = LEFT_IN_Q;
        P[j] = NEG_INF;
        cntDiffs[j] = j;

    }
    D[0] = 0;
    BT[1] = LEFT_TO_D;

    // compute remaining rows
    val_t match, tmpD, tmpP;
    lenSeqs_t diff, tmp;
    bool early;

    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        // handle left end
        D[0] = scoring.penOpen + i * scoring.penExtend;
        match = (i <= bound + 1) * (D[0] - scoring.penExtend - (i == 1) * scoring.penOpen) + (i > bound + 1) * D[i - bound - 1];
        Q[0] = NEG_INF;
        BT[0] = (i == 1) ? UP_TO_D : UP_IN_P;
        diff = (i <= bound + 1) * (i - 1) + (i > bound + 1) * cntDiffs[i - bound - 1];
        cntDiffs[0] = i;

        early = true; // early termination flag

        // fill remaining row
        for (lenSeqs_t j = 1 + (i > bound) * (i - bound - 1); j <= i + bound && j <= t.size(); j++) { // same as starting from j = max(1, i - bound) with signed integers

            // arrays D, P, Q, BT
            if (j == ((i > bound) * (i - bound))) {

                tmpP = std::max({
                                        D[j] + scoring.penOpen + scoring.penExtend,
                                        P[j] + scoring.penExtend
                                });

                Q[j] = NEG_INF;

                tmpD = std::max({
                                       match + ((s[i - 1] == t[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       tmpP
                               });

            } else if (j == (i + bound)) {

                tmpP = NEG_INF;

                Q[j] = std::max({
                                       D[j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[j - 1] + scoring.penExtend
                                });

                tmpD = std::max({
                                       match + ((s[i - 1] == t[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       Q[j]
                               });

            } else {

                tmpP = std::max({
                                        D[j] + scoring.penOpen + scoring.penExtend,
                                        P[j] + scoring.penExtend
                                });

                Q[j] = std::max({
                                        D[j - 1] + scoring.penOpen + scoring.penExtend,
                                        Q[j - 1] + scoring.penExtend
                                });

                tmpD = std::max({
                                       match + ((s[i - 1] == t[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       tmpP,
                                       Q[j]
                               });

            }

            BT[j] = (match + ((s[i - 1] == t[j - 1]) ? scoring.rewardMatch : scoring.penMismatch) == tmpD)
                    | (j != ((i > bound) * (i - bound))) * ((Q[j] == tmpD) << 1)
                    | (j != (i + bound)) * ((tmpP == tmpD) << 2)
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

    return (cntDiffs[t.size()] > bound) ? (bound + 1) : cntDiffs[t.size()];

}

lenSeqs_t Verification::computeGotohLengthAwareEarlyRow(const std::string &s, const std::string &t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], char BT[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]) {

    // long computation not necessary if lengths differ too much
    if (((s.size() > t.size()) ? (s.size() - t.size()) : (t.size() - s.size())) > bound) {
        return bound + 1;
    }

    if (bound == 0) {
        return lenSeqs_t(s != t);
    }

    std::string shorter = s;
    std::string longer = t;
    if (shorter.size() > longer.size()) {
        shorter.swap(longer);
    }
    lenSeqs_t delta = longer.size() - shorter.size();

    // initialise necessary section of first row
    D[0] = scoring.penOpen;
    for (lenSeqs_t j = 1; j <= (bound + delta) / 2 && j <= longer.size(); j++) {

        D[j] = D[j - 1] + scoring.penExtend;
        BT[j] = LEFT_IN_Q;
        P[j] = NEG_INF;
        cntDiffs[j] = j;

    }
    D[0] = 0;
    BT[1] = LEFT_TO_D;

    // compute remaining rows
    val_t match, tmpD, tmpP = NEG_INF;
    lenSeqs_t diff, tmp;
    bool early;

    for (lenSeqs_t i = 1; i <= shorter.size(); i++) {

        // handle left end
        D[0] = scoring.penOpen + i * scoring.penExtend;
        match = (i <= (bound - delta) / 2 + 1) * (D[0] - scoring.penExtend - (i == 1) * scoring.penOpen) + (i > (bound - delta) / 2 + 1) * D[i - (bound - delta) / 2 - 1];
        Q[0] = NEG_INF;
        BT[0] = (i == 1) ? UP_TO_D : UP_IN_P;
        diff = (i <= (bound - delta) / 2 + 1) * (i - 1) + (i > (bound - delta) / 2 + 1) * cntDiffs[i - (bound - delta) / 2 - 1];
        cntDiffs[0] = i;

        early = true; // early termination flag

        // fill remaining row
        for (lenSeqs_t j = 1 + (i > (bound - delta) / 2) * (i - (bound - delta) / 2 - 1); j <= i + (bound + delta) / 2 && j <= longer.size(); j++) { // same as starting from j = max(1, i - (bound - delta) / 2) with signed integers

            // arrays D, P, Q, BT
            if ((bound - delta) / 2 == 0 && (bound + delta) / 2 == 0) {

                tmpD = match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch); // (mis)match is only possibility since we have to consider only one diagonal

            } else if (j == ((i > (bound - delta) / 2) * (i - (bound - delta) / 2))) {

                tmpP = std::max({
                                        D[j] + scoring.penOpen + scoring.penExtend,
                                        P[j] + scoring.penExtend
                                });

                Q[j] = NEG_INF;

                tmpD = std::max({
                                       match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       tmpP
                               });

            } else if (j == (i + (bound + delta) / 2)) {

                tmpP = NEG_INF;

                Q[j] = std::max({
                                        D[j - 1] + scoring.penOpen + scoring.penExtend,
                                        Q[j - 1] + scoring.penExtend
                                });

                tmpD = std::max({
                                       match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       Q[j]
                               });

            } else {

                tmpP = std::max({
                                        D[j] + scoring.penOpen + scoring.penExtend,
                                        P[j] + scoring.penExtend
                                });

                Q[j] = std::max({
                                        D[j - 1] + scoring.penOpen + scoring.penExtend,
                                        Q[j - 1] + scoring.penExtend
                                });

                tmpD = std::max({
                                       match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       tmpP,
                                       Q[j]
                               });

            }

            BT[j] = (match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch) == tmpD)
                           | ((bound - delta) / 2 != 0 || (bound + delta) / 2 != 0) * (j != ((i > (bound - delta) / 2) * (i - (bound - delta) / 2))) * ((Q[j] == tmpD) << 1)
                           | ((bound - delta) / 2 != 0 || (bound + delta) / 2 != 0) * (j != (i + (bound + delta) / 2)) * ((tmpP == tmpD) << 2)
                           | ((bound - delta) / 2 != 0 || (bound + delta) / 2 != 0) * (((D[j] + scoring.penOpen + scoring.penExtend) == tmpP) << 3)
                           | ((bound - delta) / 2 != 0 || (bound + delta) / 2 != 0) * (((P[j] + scoring.penExtend) == tmpP) << 4)
                           | ((bound - delta) / 2 != 0 || (bound + delta) / 2 != 0) * (((D[j - 1] + scoring.penOpen + scoring.penExtend) == Q[j]) << 5)
                           | ((bound - delta) / 2 != 0 || (bound + delta) / 2 != 0) * (((Q[j - 1] + scoring.penExtend) == Q[j]) << 6);

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

                tmp = diff + (shorter[i - 1] != longer[j - 1]);
                diff = cntDiffs[j];
                cntDiffs[j] = tmp;

            }

            early &= ((cntDiffs[j] + ((delta + i >= j) * (delta + i - j) + (delta + i < j) * (j - delta - i))) > bound); // improved e.t.

        }

        if (early) {// computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }

    }

    return (cntDiffs[longer.size()] > bound) ? (bound + 1) : cntDiffs[longer.size()];

}

lenSeqs_t Verification::computeGotohLengthAwareEarlyRow2(const std::string &s, const std::string &t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], char BT[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]) {

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
    lenSeqs_t delta = longer.size() - shorter.size();

    // (mis)match is the only possibility when we have to consider only one diagonal [happens only if (a) bound = delta = 0, or (b) bound = 1 and delta = 0, but (a) is already covered above]
    if ((bound - delta) / 2 == 0 && (bound + delta) / 2 == 0) {

        lenSeqs_t diffs = 0;
        for (auto i = 0; diffs <= bound && i < shorter.size(); i++) {
            diffs += (shorter[i] != longer[i]);
        }

        return diffs;

    }

    // initialise necessary section of first row
    D[0] = scoring.penOpen;
    for (lenSeqs_t j = 1; j <= (bound + delta) / 2 && j <= longer.size(); j++) {

        D[j] = D[j - 1] + scoring.penExtend;
        BT[j] = LEFT_IN_Q;
        P[j] = NEG_INF;
        cntDiffs[j] = j;

    }
    D[0] = 0;
    BT[1] = LEFT_TO_D;

    // compute remaining rows
    val_t match, tmpD, tmpP;
    lenSeqs_t diff, tmp;
    bool early;

    for (lenSeqs_t i = 1; i <= shorter.size(); i++) {

        // handle left end
        D[0] = scoring.penOpen + i * scoring.penExtend;
        match = (i <= (bound - delta) / 2 + 1) * (D[0] - scoring.penExtend - (i == 1) * scoring.penOpen) + (i > (bound - delta) / 2 + 1) * D[i - (bound - delta) / 2 - 1];
        Q[0] = NEG_INF;
        BT[0] = (i == 1) ? UP_TO_D : UP_IN_P;
        diff = (i <= (bound - delta) / 2 + 1) * (i - 1) + (i > (bound - delta) / 2 + 1) * cntDiffs[i - (bound - delta) / 2 - 1];
        cntDiffs[0] = i;

        early = true; // early termination flag

        // fill remaining row
        for (lenSeqs_t j = 1 + (i > (bound - delta) / 2) * (i - (bound - delta) / 2 - 1); j <= i + (bound + delta) / 2 && j <= longer.size(); j++) { // same as starting from j = max(1, i - (bound - delta) / 2) with signed integers

            // arrays D, P, Q, BT
            if (j == ((i > (bound - delta) / 2) * (i - (bound - delta) / 2))) {

                tmpP = std::max({
                                       D[j] + scoring.penOpen + scoring.penExtend,
                                       P[j] + scoring.penExtend
                                });

                Q[j] = NEG_INF;

                tmpD = std::max({
                                       match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       tmpP
                               });

            } else if (j == (i + (bound + delta) / 2)) {

                tmpP = NEG_INF;

                Q[j] = std::max({
                                       D[j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[j - 1] + scoring.penExtend
                                       });

                tmpD = std::max({
                                       match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       Q[j]
                               });

            } else {

                tmpP = std::max({
                                       D[j] + scoring.penOpen + scoring.penExtend,
                                       P[j] + scoring.penExtend
                                });

                Q[j] = std::max({
                                       D[j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[j - 1] + scoring.penExtend
                                       });

                tmpD = std::max({
                                       match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       tmpP,
                                       Q[j]
                               });

            }

            BT[j] = (match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch) == tmpD)
                           | (j != ((i > (bound - delta) / 2) * (i - (bound - delta) / 2))) * ((Q[j] == tmpD) << 1)
                           | (j != (i + (bound + delta) / 2)) * ((tmpP == tmpD) << 2)
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

                tmp = diff + (shorter[i - 1] != longer[j - 1]);
                diff = cntDiffs[j];
                cntDiffs[j] = tmp;

            }

            early &= ((cntDiffs[j] + ((delta + i >= j) * (delta + i - j) + (delta + i < j) * (j - delta - i))) > bound); // improved e.t.

        }

        if (early) {// computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }

    }

    return (cntDiffs[longer.size()] > bound) ? (bound + 1) : cntDiffs[longer.size()];

}

lenSeqs_t Verification::computeGotohLengthAwareEarlyRow3(const std::string &s, const std::string &t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]) {

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
    lenSeqs_t delta = longer.size() - shorter.size();

    // (mis)match is only possibility when we have to consider only one diagonal [happens only if (a) bound = delta = 0, or (b) bound = 1 and delta = 0, but (a) is already covered above]
    if ((bound - delta) / 2 == 0 && (bound + delta) / 2 == 0) {

        lenSeqs_t diffs = 0;
        for (auto i = 0; diffs <= bound && i < shorter.size(); i++) {
            diffs += (shorter[i] != longer[i]);
        }

        return diffs;

    }

    // initialise necessary section of first row
    D[0] = scoring.penOpen;
    for (lenSeqs_t j = 1; j <= (bound + delta) / 2 && j <= longer.size(); j++) {

        D[j] = D[j - 1] + scoring.penExtend;
        P[j] = NEG_INF;
        cntDiffs[j] = j;

    }
    D[0] = 0;

    // compute remaining rows
    val_t match, tmpD, tmpP;
    lenSeqs_t diff, tmp;
    bool early;

    for (lenSeqs_t i = 1; i <= shorter.size(); i++) {

        // handle left end
        D[0] = scoring.penOpen + i * scoring.penExtend;
        match = (i <= (bound - delta) / 2 + 1) * (D[0] - scoring.penExtend - (i == 1) * scoring.penOpen) + (i > (bound - delta) / 2 + 1) * D[i - (bound - delta) / 2 - 1];
        Q[0] = NEG_INF;
        diff = (i <= (bound - delta) / 2 + 1) * (i - 1) + (i > (bound - delta) / 2 + 1) * cntDiffs[i - (bound - delta) / 2 - 1];
        cntDiffs[0] = i;

        early = true; // early termination flag

        // fill remaining row
        for (lenSeqs_t j = 1 + (i > (bound - delta) / 2) * (i - (bound - delta) / 2 - 1); j <= i + (bound + delta) / 2 && j <= longer.size(); j++) { // same as starting from j = max(1, i - (bound - delta) / 2) with signed integers

            // arrays D, P, Q, BT
            if (j == ((i > (bound - delta) / 2) * (i - (bound - delta) / 2))) {

                tmpP = std::max({
                                       D[j] + scoring.penOpen + scoring.penExtend,
                                       P[j] + scoring.penExtend
                                });

                Q[j] = NEG_INF;

                tmpD = std::max({
                                       match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       tmpP
                                });

            } else if (j == (i + (bound + delta) / 2)) {

                tmpP = NEG_INF;

                Q[j] = std::max({
                                       D[j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[j - 1] + scoring.penExtend
                                       });

                tmpD = std::max({
                                       match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       Q[j]
                               });

            } else {

                tmpP = std::max({
                                       D[j] + scoring.penOpen + scoring.penExtend,
                                       P[j] + scoring.penExtend
                                });

                Q[j] = std::max({
                                       D[j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[j - 1] + scoring.penExtend
                                       });

                tmpD = std::max({
                                       match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       tmpP,
                                       Q[j]
                               });

            }

            // arrays cntDiffs, cntDiffsP, cntDiffsQ
            if (tmpP == (D[j] + scoring.penOpen + scoring.penExtend)) {
                cntDiffsP[j] = cntDiffs[j] + 1;
            } else {
                cntDiffsP[j] = cntDiffsP[j] + 1;
            }

            if (Q[j] ==(D[j - 1] + scoring.penOpen + scoring.penExtend)) {
                cntDiffsQ[j] = cntDiffs[j - 1] + 1;
            } else {
                cntDiffsQ[j] = cntDiffsQ[j - 1] + 1;
            }

            if ((j != ((i > (bound - delta) / 2) * (i - (bound - delta) / 2))) && (tmpD == Q[j])) {

                diff = cntDiffs[j];
                cntDiffs[j] = cntDiffsQ[j];

            } else if ((j != (i + (bound + delta) / 2)) && (tmpD == tmpP)) {

                diff = cntDiffs[j];
                cntDiffs[j] = cntDiffsP[j];

            } else if (tmpD == match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch)) {

                tmp = diff + (shorter[i - 1] != longer[j - 1]);
                diff = cntDiffs[j];
                cntDiffs[j] = tmp;

            }

            P[j] = tmpP;

            match = D[j];
            D[j] = tmpD;

            early &= ((cntDiffs[j] + ((delta + i >= j) * (delta + i - j) + (delta + i < j) * (j - delta - i))) > bound); // improved e.t.

        }

        if (early) {// computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }

    }

    return (cntDiffs[longer.size()] > bound) ? (bound + 1) : cntDiffs[longer.size()];

}

lenSeqs_t Verification::computeGotohLengthAwareEarlyRow4(const std::string &s, const std::string &t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]) {

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
    lenSeqs_t delta = longer.size() - shorter.size();

    // (mis)match is only possibility when we have to consider only one diagonal [happens only if (a) bound = delta = 0, or (b) bound = 1 and delta = 0, but (a) is already covered above]
    if ((bound - delta) / 2 == 0 && (bound + delta) / 2 == 0) {

        lenSeqs_t diffs = 0;
        for (auto i = 0; diffs <= bound && i < shorter.size(); i++) {
            diffs += (shorter[i] != longer[i]);
        }

        return diffs;

    }

    // initialise necessary section of first row
    D[0] = scoring.penOpen;
    for (lenSeqs_t j = 1; j <= (bound + delta) / 2 && j <= longer.size(); j++) {

        D[j] = D[j - 1] + scoring.penExtend;
        P[j] = NEG_INF;
        cntDiffs[j] = j;

    }
    D[0] = 0;

    // compute remaining rows
    val_t match, tmpD, tmpP;
    lenSeqs_t j, diff, tmp;
    bool early;

    for (lenSeqs_t i = 1; i <= shorter.size(); i++) {

        // handle left end
        D[0] = scoring.penOpen + i * scoring.penExtend;
        match = (i <= (bound - delta) / 2 + 1) * (D[0] - scoring.penExtend - (i == 1) * scoring.penOpen) + (i > (bound - delta) / 2 + 1) * D[i - (bound - delta) / 2 - 1];
        Q[0] = NEG_INF;
        diff = (i <= (bound - delta) / 2 + 1) * (i - 1) + (i > (bound - delta) / 2 + 1) * cntDiffs[i - (bound - delta) / 2 - 1];
        cntDiffs[0] = i;

        early = true; // early termination flag

        // fill remaining row
        j = 1 + (i > (bound - delta) / 2) * (i - (bound - delta) / 2 - 1);
        D[j - 1] = Q[j - 1] = NEG_INF;
        if (i + (bound + delta) / 2 <= longer.size()) {D[i + (bound + delta) / 2] = P[i + (bound + delta) / 2] = NEG_INF;}
        for (; j <= i + (bound + delta) / 2 && j <= longer.size(); j++) { // same as starting from j = max(1, i - (bound - delta) / 2) with signed integers

            // arrays D, P, Q, BT
            tmpP = std::max({
                                           D[j] + scoring.penOpen + scoring.penExtend,
                                           P[j] + scoring.penExtend
                            });

            Q[j] = std::max({
                                           D[j - 1] + scoring.penOpen + scoring.penExtend,
                                           Q[j - 1] + scoring.penExtend
                            });

            tmpD = std::max({
                                   match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                   tmpP,
                                   Q[j]
                           });


            // arrays cntDiffs, cntDiffsP, cntDiffsQ
            if (tmpP == (D[j] + scoring.penOpen + scoring.penExtend)) {
                cntDiffsP[j] = cntDiffs[j] + 1;
            } else {
                cntDiffsP[j] = cntDiffsP[j] + 1;
            }

            if (Q[j] == (D[j - 1] + scoring.penOpen + scoring.penExtend)) {
                cntDiffsQ[j] = cntDiffs[j - 1] + 1;
            } else {
                cntDiffsQ[j] = cntDiffsQ[j - 1] + 1;
            }

            if (tmpD == Q[j]) {

                diff = cntDiffs[j];
                cntDiffs[j] = cntDiffsQ[j];

            } else if (tmpD == tmpP) {

                diff = cntDiffs[j];
                cntDiffs[j] = cntDiffsP[j];

            } else if (tmpD == match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch)) {

                tmp = diff + (shorter[i - 1] != longer[j - 1]);
                diff = cntDiffs[j];
                cntDiffs[j] = tmp;

            }

            P[j] = tmpP;

            match = D[j];
            D[j] = tmpD;

            early &= ((cntDiffs[j] + labs((long) delta + (long) i - (long) j)) > bound); // improved e.t.

        }

        if (early) {// computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }

    }

    return (cntDiffs[longer.size()] > bound) ? (bound + 1) : cntDiffs[longer.size()];

}

lenSeqs_t Verification::computeGotohLengthAwareEarlyRow5(const std::string &s, const std::string &t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]) {

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
    lenSeqs_t delta = longer.size() - shorter.size();

    // (mis)match is only possibility when we have to consider only one diagonal [happens only if (a) bound = delta = 0, or (b) bound = 1 and delta = 0, but (a) is already covered above]
    if ((bound - delta) / 2 == 0 && (bound + delta) / 2 == 0) {

        lenSeqs_t diffs = 0;
        for (auto i = 0; diffs <= bound && i < shorter.size(); i++) {
            diffs += (shorter[i] != longer[i]);
        }

        return diffs;

    }

    lenSeqs_t bl = (bound - delta) / 2;
    lenSeqs_t br = (bound + delta) / 2;

    // initialise necessary section of first row
    D[0] = scoring.penOpen;
    for (lenSeqs_t j = 1; j <= br && j <= longer.size(); j++) {

        D[j] = D[j - 1] + scoring.penExtend;
        P[j] = NEG_INF;
        cntDiffs[j] = j;

    }
    D[0] = 0;

    // compute remaining rows
    val_t match, tmpD, tmpP;
    lenSeqs_t diff, tmp;
    bool early;

    for (lenSeqs_t i = 1; i <= shorter.size(); i++) {

        // handle left end
        D[0] = scoring.penOpen + i * scoring.penExtend;
        match = (i <= bl + 1) * (D[0] - scoring.penExtend - (i == 1) * scoring.penOpen) + (i > bl + 1) * D[i - bl - 1];
        Q[0] = NEG_INF;
        diff = (i <= bl + 1) * (i - 1) + (i > bl + 1) * cntDiffs[i - bl - 1];
        cntDiffs[0] = i;

        early = true; // early termination flag

        // fill remaining row
        lenSeqs_t j = 1 + (i > bl) * (i - bl - 1);
        D[j - 1] = Q[j - 1] = NEG_INF;
        if (i + br <= longer.size()) {D[i + br] = P[i + br] = NEG_INF;}
        for (; j <= i + br && j <= longer.size(); j++) { // same as starting from j = max(1, i - bl) with signed integers

            // arrays D, P, Q, BT
            tmpP = std::max({
                                           D[j] + scoring.penOpen + scoring.penExtend,
                                           P[j] + scoring.penExtend
                            });

            Q[j] = std::max({
                                           D[j - 1] + scoring.penOpen + scoring.penExtend,
                                           Q[j - 1] + scoring.penExtend
                            });

            tmpD = std::max({
                                   match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                   tmpP,
                                   Q[j]
                           });


            // arrays cntDiffs, cntDiffsP, cntDiffsQ
            if (tmpP == (D[j] + scoring.penOpen + scoring.penExtend)) {
                cntDiffsP[j] = cntDiffs[j] + 1;
            } else {
                cntDiffsP[j] = cntDiffsP[j] + 1;
            }

            if (Q[j] == (D[j - 1] + scoring.penOpen + scoring.penExtend)) {
                cntDiffsQ[j] = cntDiffs[j - 1] + 1;
            } else {
                cntDiffsQ[j] = cntDiffsQ[j - 1] + 1;
            }

            if (tmpD == Q[j]) {

                diff = cntDiffs[j];
                cntDiffs[j] = cntDiffsQ[j];

            } else if (tmpD == tmpP) {

                diff = cntDiffs[j];
                cntDiffs[j] = cntDiffsP[j];

            } else if (tmpD == match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch)) {

                tmp = diff + (shorter[i - 1] != longer[j - 1]);
                diff = cntDiffs[j];
                cntDiffs[j] = tmp;

            }

            P[j] = tmpP;

            match = D[j];
            D[j] = tmpD;

            early &= ((cntDiffs[j] + labs((long) delta + (long) i - (long) j)) > bound); // improved e.t.

        }

        if (early) {// computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }

    }

    return (cntDiffs[longer.size()] > bound) ? (bound + 1) : cntDiffs[longer.size()];

}

lenSeqs_t Verification::computeGotohLengthAwareEarlyRow6(const std::string &s, const std::string &t, const lenSeqs_t bound, const Scoring& scoring, val_t D[], val_t P[], val_t Q[], lenSeqs_t cntDiffs[], lenSeqs_t cntDiffsP[], lenSeqs_t cntDiffsQ[]) {

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
    lenSeqs_t delta = longer.size() - shorter.size();

    // (mis)match is only possibility when we have to consider only one diagonal [happens only if (a) bound = delta = 0, or (b) bound = 1 and delta = 0, but (a) is already covered above]
    if ((bound - delta) / 2 == 0 && (bound + delta) / 2 == 0) {

        lenSeqs_t diffs = 0;
        for (auto i = 0; diffs <= bound && i < shorter.size(); i++) {
            diffs += (shorter[i] != longer[i]);
        }

        return diffs;

    }

    // initialise necessary section of first row
    D[0] = scoring.penOpen;
    for (lenSeqs_t j = 1; j <= (bound + delta) / 2 && j <= longer.size(); j++) {

        D[j] = D[j - 1] + scoring.penExtend;
        P[j] = NEG_INF;
        cntDiffs[j] = j;

    }
    D[0] = 0;

    // compute remaining rows
    val_t match, maxVal;
    lenSeqs_t j, diff, maxDiff;
    bool early;

    for (lenSeqs_t i = 1; i <= shorter.size(); i++) {

        // handle left end
        D[0] = scoring.penOpen + i * scoring.penExtend;
        match = (i <= (bound - delta) / 2 + 1) * (D[0] - scoring.penExtend - (i == 1) * scoring.penOpen) + (i > (bound - delta) / 2 + 1) * D[i - (bound - delta) / 2 - 1];
        Q[0] = NEG_INF;
        diff = (i <= (bound - delta) / 2 + 1) * (i - 1) + (i > (bound - delta) / 2 + 1) * cntDiffs[i - (bound - delta) / 2 - 1];
        cntDiffs[0] = i;

        early = true; // early termination flag

        // fill remaining row
        j = 1 + (i > (bound - delta) / 2) * (i - (bound - delta) / 2 - 1);
        D[j - 1] = Q[j - 1] = NEG_INF;
        if (i + (bound + delta) / 2 <= longer.size()) {
            D[i + (bound + delta) / 2] = P[i + (bound + delta) / 2] = NEG_INF;
        }
        for (; j <= i + (bound + delta) / 2 && j <= longer.size(); j++) { // same as starting from j = max(1, i - (bound - delta) / 2) with signed integers

            // arrays D, P, Q, cntDiffs, cntDiffsP, cntDiffsQ
            if ((D[j] + scoring.penOpen + scoring.penExtend) >= (P[j] + scoring.penExtend)) {

                P[j] = D[j] + scoring.penOpen + scoring.penExtend;
                cntDiffsP[j] = cntDiffs[j] + 1;

            } else {

                P[j] = P[j] + scoring.penExtend;
                cntDiffsP[j] = cntDiffsP[j] + 1;

            }

            if ((D[j - 1] + scoring.penOpen + scoring.penExtend) >= (Q[j - 1] + scoring.penExtend)) {

                Q[j] = D[j - 1] + scoring.penOpen + scoring.penExtend;
                cntDiffsQ[j] = cntDiffs[j - 1] + 1;

            } else {

                Q[j] = Q[j - 1] + scoring.penExtend;
                cntDiffsQ[j] = cntDiffsQ[j - 1] + 1;

            }

            maxVal = (match + ((shorter[i - 1] == longer[j - 1]) ? scoring.rewardMatch : scoring.penMismatch));
            maxDiff = diff + (shorter[i - 1] != longer[j - 1]);
            if (P[j] >= maxVal) {

                maxVal = P[j];
                maxDiff = cntDiffsP[j];

            }
            if (Q[j] >= maxVal){

                maxVal = Q[j];
                maxDiff = cntDiffsQ[j];

            }

            match = D[j];
            D[j] = maxVal;

            diff = cntDiffs[j];
            cntDiffs[j] = maxDiff;

            early &= ((cntDiffs[j] + labs((long) delta + (long) i - (long) j)) > bound); // improved e.t.

        }

        if (early) {// computation can be terminated early if computed row contains only values > bound (because values are monotonically increasing)
            return bound + 1;
        }

    }

    return (cntDiffs[longer.size()] > bound) ? (bound + 1) : cntDiffs[longer.size()];

}

void Verification::verifyGotoh(const AmpliconCollection& ac, Matches& mat, Buffer<Candidate>& buf, lenSeqs_t width, lenSeqs_t t, const Scoring& scoring) {

    Candidate c;
    Buffer<Candidate> localBuffer;
    Verification::val_t D[width];
    Verification::val_t P[width];
    Verification::val_t Q[width];
    lenSeqs_t cntDiffs[width];
    lenSeqs_t cntDiffsP[width];
    lenSeqs_t cntDiffsQ[width];

    while (!buf.isClosed() || buf.syncSize() > 0) {

        buf.syncSwapContents(localBuffer);

        while (localBuffer.size() > 0) {

            c = localBuffer.pop();

            if (!mat.contains(c.first, c.second)) {

                lenSeqs_t d = computeGotohLengthAwareEarlyRow6(ac[c.first].seq, ac[c.second].seq, t, scoring, D, P, Q, cntDiffs, cntDiffsP, cntDiffsQ); //TODO choose "best" implementation

                if (d <= t) {
                    mat.add(c.first, c.second, d);
                }

            }

        }

    }

}



Verification::AlignmentInformation Verification::computeGotohCigarFull(const std::string &s, const std::string &t, const Scoring& scoring) {

    ValArray D(s.size() + 1, std::vector<val_t>(t.size() + 1));
    ValArray P(s.size() + 1, std::vector<val_t>(t.size() + 1));
    ValArray Q(s.size() + 1, std::vector<val_t>(t.size() + 1));
    BacktrackingArray BT(s.size() + 1, std::vector<char>(t.size() + 1));

    std::vector<std::pair<char, lenSeqs_t >> cigarSegs;
    lenSeqs_t len = 0;
    lenSeqs_t numDiffs = 0;

    // initialisation
    D[0][0] = scoring.penOpen;
    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        D[i][0] = D[i - 1][0] + scoring.penExtend;
        BT[i][0] = UP_IN_P;
        Q[i][0] = NEG_INF;

    }
    for (lenSeqs_t j = 1; j <= t.size(); j++) {

        D[0][j] = D[0][j - 1] + scoring.penExtend;
        BT[0][j] = LEFT_IN_Q;
        P[0][j] = NEG_INF;

    }
    D[0][0] = 0;
    BT[1][0] = UP_TO_D;
    BT[0][1] = LEFT_TO_D;

    // fill matrix
    for (lenSeqs_t i = 1; i <= s.size(); i++) {

        for (lenSeqs_t j = 1; j <= t.size(); j++) {

            P[i][j] = std::max({
                                       D[i - 1][j] + scoring.penOpen + scoring.penExtend,
                                       P[i - 1][j] + scoring.penExtend
                               });

            Q[i][j] = std::max({
                                       D[i][j - 1] + scoring.penOpen + scoring.penExtend,
                                       Q[i][j - 1] + scoring.penExtend
                               });

            D[i][j] = std::max({
                                       D[i - 1][j - 1] + ((s[i - 1] == t[j - 1]) ? scoring.rewardMatch : scoring.penMismatch),
                                       P[i][j],
                                       Q[i][j]
                               });

            BT[i][j] = (D[i - 1][j - 1] + ((s[i - 1] == t[j - 1]) ? scoring.rewardMatch : scoring.penMismatch) == D[i][j])
                       | ((Q[i][j] == D[i][j]) << 1)
                       | ((P[i][j] == D[i][j]) << 2)
                       | (((D[i - 1][j] + scoring.penOpen + scoring.penExtend) == P[i][j]) << 3)
                       | (((P[i - 1][j] + scoring.penExtend) == P[i][j]) << 4)
                       | (((D[i][j - 1] + scoring.penOpen + scoring.penExtend) == Q[i][j]) << 5)
                       | (((Q[i][j - 1] + scoring.penExtend) == Q[i][j]) << 6);

        }

    }

    // backtracking (priorities: left > up > diagonal)
    lenSeqs_t i = s.size();
    lenSeqs_t j = t.size();

    long MATRIX = 0;

    while (i != 0 && j != 0) {

        switch (MATRIX) {

            case IN_D: {

                if (BT[i][j] & JUMP_TO_Q) {

                    MATRIX = IN_Q;

                } else if (BT[i][j] & JUMP_TO_P) {

                    MATRIX = IN_P;

                } else if (BT[i][j] & DIAGONAL_IN_D) {

                    len++;
                    numDiffs += (s[i - 1] != t[j - 1]);
                    if ((cigarSegs.size() > 0) && (cigarSegs.back().first == 'M')) {
                        cigarSegs.back().second++;
                    } else {
                        cigarSegs.push_back(std::make_pair('M', 1));
                    }

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

                len++;
                numDiffs++;
                if ((cigarSegs.size() > 0) && (cigarSegs.back().first == 'D')) {
                    cigarSegs.back().second++;
                } else {
                    cigarSegs.push_back(std::make_pair('D', 1));
                }

                i--;

                break;

            }

            case IN_Q: {

                /*if (BT[i][j] & LEFT_IN_Q) {

                } else*/ if (BT[i][j] & LEFT_TO_D) {

                    MATRIX = IN_D;
                }

                len++;
                numDiffs++;
                if ((cigarSegs.size() > 0) && (cigarSegs.back().first == 'I')) {
                    cigarSegs.back().second++;
                } else {
                    cigarSegs.push_back(std::make_pair('I', 1));
                }

                j--;

                break;

            }

            default: {
                // do nothing (should not occur)
            }

        }

    }

    while (i > 0) {

        len++;
        numDiffs++;
        if ((cigarSegs.size() > 0) && (cigarSegs.back().first == 'd')) {
            cigarSegs.back().second++;
        } else {
            cigarSegs.push_back(std::make_pair('d', 1));
        }

        i--;

    }

    while (j > 0) {

        len++;
        numDiffs++;
        if ((cigarSegs.size() > 0) && (cigarSegs.back().first == 'i')) {
            cigarSegs.back().second++;
        } else {
            cigarSegs.push_back(std::make_pair('i', 1));
        }

        j--;

    }

    std::stringstream sStream;
    for (auto iter = cigarSegs.rbegin(); iter != cigarSegs.rend(); iter++) {

        if (iter->second > 1) sStream << iter->second;
        sStream << iter->first;

    }

    return AlignmentInformation(sStream.str(), len, numDiffs);

}


// Slightly adapted nw(...) code taken from swarm.
// TODO Reimplement properly or change above method computeGotohCigarFull(...).
void pushop(char newop, char ** cigarendp, char * op, int * count)
{
    if (newop == *op)
        (*count)++;
    else
    {
        *--*cigarendp = *op;
        if (*count > 1)
        {
            char buf[25];
            int len = sprintf(buf, "%d", *count);
            *cigarendp -= len;
            strncpy(*cigarendp, buf, len);
        }
        *op = newop;
        *count = 1;
    }
}

void finishop(char ** cigarendp, char * op, int * count)
{
    if ((op) && (count))
    {
        *--*cigarendp = *op;
        if (*count > 1)
        {
            char buf[25];
            int len = sprintf(buf, "%d", *count);
            *cigarendp -= len;
            strncpy(*cigarendp, buf, len);
        }
        *op = 0;
        *count = 0;
    }
}

Verification::AlignmentInformation Verification::computeGotohCigarSwarm(const std::string &s, const std::string &t, const Scoring& scoring) {
//begin - "external"
    unsigned char maskup      = 1;
    unsigned char maskleft    = 2;
    unsigned char maskextup   = 4;
    unsigned char maskextleft = 8;

    unsigned long gapopen = 12;
    unsigned long gapextend = 4;
    unsigned long penalty_mismatch = 9;



    char map_nt[256] =
            {
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1,  1, -1,  2, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1,  4,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1,  1, -1,  2, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1,  4,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
            };
    std::vector<long> score_matrix(32 * 32);
    for (int a = 0; a < 16; a++) {
        for (int b = 0; b < 16; b++) {
            score_matrix[(a << 5) + b] = (((a == b) && (a > 0) && (b > 0)) ? 0 : penalty_mismatch);
        }
    }


    long longestamplicon = std::max(s.length(), t.length());
    unsigned char * dir = new unsigned char[longestamplicon * longestamplicon];
    unsigned long * hearray = new unsigned long[2 * longestamplicon];

//end - "external"

    /* dir must point to at least qlen*dlen bytes of allocated memory
    hearray must point to at least 2*qlen longs of allocated memory (8*qlen bytes) */

    long n, e;
    long unsigned *hep;

    long qlen = s.length();
    long dlen = t.length();

    memset(dir, 0, qlen*dlen);

    long i, j;

    for(i=0; i<qlen; i++)
    {
        hearray[2*i]   = 1 * gapopen + (i+1) * gapextend; // H (N)
        hearray[2*i+1] = 2 * gapopen + (i+2) * gapextend; // E
    }

    for(j=0; j<dlen; j++)
    {
        hep = &(hearray[0]);
        long f = 2 * gapopen + (j+2) * gapextend;
        long h = (j == 0) ? 0 : (gapopen + j * gapextend);

        for(i=0; i<qlen; i++)
        {
            long index = qlen*j+i;

            n = *hep;
            e = *(hep+1);
            h += score_matrix[(map_nt[t[j]]<<5) + map_nt[s[i]]]; //h += score_matrix[(t[j]<<5) + s[i]];


            dir[index] |= (f < h ? maskup : 0);
            h = std::min(h, f);
            h = std::min(h, e);
            dir[index] |= (e == h ? maskleft : 0);

            *hep = h;

            h += gapopen + gapextend;
            e += gapextend;
            f += gapextend;

            dir[index] |= (f < h ? maskextup : 0);
            dir[index] |= (e < h ? maskextleft : 0);
            f = std::min(h,f);
            e = std::min(h,e);

            *(hep+1) = e;
            h = n;
            hep += 2;
        }
    }

    long dist = hearray[2*qlen-2];

    /* backtrack: count differences and save alignment in cigar string */

    long score = 0;
    long alength = 0;
    long matches = 0;

    char * cigar = new char[qlen + dlen + 1];//(char *) xmalloc(qlen + dlen + 1);
    char * cigarend = cigar+qlen+dlen+1;


    char op = 0;
    int count = 0;
    *(--cigarend) = 0;

    i = qlen;
    j = dlen;

    while ((i>0) && (j>0))
    {
        int d = dir[qlen*(j-1)+(i-1)];

        alength++;

        if ((op == 'I') && (d & maskextleft))
        {
            score += gapextend;
            j--;
            pushop('I', &cigarend, &op, &count);
        }
        else if ((op == 'D') && (d & maskextup))
        {
            score += gapextend;
            i--;
            pushop('D', &cigarend, &op, &count);
        }
        else if (d & maskleft)
        {
            score += gapextend;
            if (op != 'I')
                score += gapopen;
            j--;
            pushop('I', &cigarend, &op, &count);
        }
        else if (d & maskup)
        {
            score += gapextend;
            if (op != 'D')
                score +=gapopen;
            i--;
            pushop('D', &cigarend, &op, &count);
        }
        else
        {
            score += score_matrix[(map_nt[t[j-1]] << 5) + map_nt[s[i-1]]]; //score += score_matrix[(t[j-1] << 5) + s[i-1]];
            if (s[i-1] == t[j-1])
                matches++;
            i--;
            j--;
            pushop('M', &cigarend, &op, &count);
        }
    }

    while(i>0)
    {
        alength++;
        score += gapextend;
        if (op != 'D')
            score += gapopen;
        i--;
        pushop('D', &cigarend, &op, &count);
    }

    while(j>0)
    {
        alength++;
        score += gapextend;
        if (op != 'I')
            score += gapopen;
        j--;
        pushop('I', &cigarend, &op, &count);

    }

    finishop(&cigarend, &op, &count);


    /* move and reallocate cigar */

    long cigarlength = cigar+qlen+dlen-cigarend;
//    memmove(cigar, cigarend, cigarlength+1);
//    cigar = (char*) xrealloc(cigar, cigarlength+1);
    char* finalCigar = new char[cigarlength + 1];
    memmove(finalCigar, cigarend, cigarlength + 1);

    std::string cigStr = std::string(finalCigar);

    delete[] cigar;
    delete[] finalCigar;
    delete[] dir;
    delete[] hearray;

    return AlignmentInformation(cigStr, alength, alength - matches);

}

}