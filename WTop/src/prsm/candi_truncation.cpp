//
// Created by LYC on 2025/4/3.
//

#include "candi_truncation.h"
using namespace std;

int candi_trucation::getCutN() const {
    return cut_n;
}
void candi_trucation::setCutN(int cutN) {
    cut_n = cutN;
}
int candi_trucation::getCutC() const {
    return cut_c;
}
void candi_trucation::setCutC(int cutC) {
    cut_c = cutC;
}
double candi_trucation::getUnkmass() const {
    return unkmass;
}
void candi_trucation::setUnkmass(double unkmass) {
    candi_trucation::unkmass = unkmass;
}
int candi_trucation::getMatchPeaks() const {
    return match_peaks;
}
void candi_trucation::setMatchPeaks(int matchPeaks) {
    match_peaks = matchPeaks;
}
double candi_trucation::getVarptmmass() const {
    return varptmmass;
}
void candi_trucation::setVarptmmass(double varptmmass) {
    candi_trucation::varptmmass = varptmmass;
}
candi_trucation::candi_trucation(int cutN, int cutC, double unkmass, double varptmmass, int matchPeaks) :
    cut_n(cutN), cut_c(cutC), unkmass(unkmass), varptmmass(varptmmass), match_peaks(matchPeaks) {
}
