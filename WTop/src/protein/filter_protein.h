//
// Created by wenzhong on 2023/3/21.
//

#ifndef WTOP_FILTER_PROTEIN_H
#define WTOP_FILTER_PROTEIN_H

#include "protein.hpp"

class filter_protein {
public:
    filter_protein(protein *proteinPtr) :
        proteinPtr(proteinPtr) {
    }

    filter_protein(protein *proteinPtr, double ptmMassGap) :
        proteinPtr(proteinPtr), ptmMassGap(ptmMassGap) {
    }

    void setFilterTagLength(int filterTagLength) {
        filter_tag_length = filterTagLength;
    }

    void setFilterMatchIons(int filterMatchIons) {
        filter_match_ions = filterMatchIons;
    }

    void setFilterTagNum(int filterTagNum) {
        filter_tag_num = filterTagNum;
    }

    void setComplementIonsInProtein(double complementIonsInProtein) {
        filter_protein::complementIonsInProtein = complementIonsInProtein;
    }

    protein *getProteinPtr() {
        return proteinPtr;
    }

    int getFilterTagLength() const {
        return filter_tag_length;
    }

    int getFilterMatchIons() const {
        return filter_match_ions;
    }

    int getFilterTagNum() const {
        return filter_tag_num;
    }

    double getComplementIonsInProtein() const {
        return complementIonsInProtein;
    }

    bool isLogLeastScope() const {
        return logLeastScope;
    }

    void setLogLeastScope(bool logLeastScope) {
        filter_protein::logLeastScope = logLeastScope;
    }

    void setOnePtmMatchPeaks(int onePtmMatchPeaks) {
        one_ptm_match_peaks = onePtmMatchPeaks;
    }

    int getOnePtmMatchPeaks() const {
        return one_ptm_match_peaks;
    }

    double getMassGap() const {
        return mass_gap;
    }

    void setMassGap(double massGap) {
        mass_gap = massGap;
    }

    double getPtmMassGap() const {
        return ptmMassGap;
    }

    protein_ptr proteinPtr;

    string tagN, tagC;
    int filter_tag_length = 0; //最长tag的长度
    double tagMassError = 0.0;

    int filter_match_ions = 0; //完全匹配离子的个数
    int filter_tag_num = 0;    //匹配上tag的个数
    double complementIonsInProtein = 0.0;
    bool logLeastScope; //是否满足范围

    double mass_gap = 0.0;
    int one_ptm_match_peaks = 0;

    double ptmMassGap = 0.0;

private:
};

typedef filter_protein *filter_protein_ptr;

#endif //WTOP_FILTER_PROTEIN_H
