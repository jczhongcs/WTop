//
// Created by Administrator on 2023/3/31.
//

#include "ptm_filter.hpp"
#include <algorithm>
#include "../mylib/data_steam.hpp"
#include "unordered_set"
#include "../util/ppm_value.hpp"

bool isPtmConMass(vector<double> &theoMassVec, vector<double> &ptmCMassVec,
                  double precursorMass, double &ptmCMass,
                  double maxPtmMass, double minPtmMass, double minPPMValue) {
    int startIndex = mylib::data_stream::txpb_binarey_search_ex(theoMassVec, theoMassVec.size(), precursorMass + maxPtmMass);
    if (startIndex == -1 && startIndex > theoMassVec.size()) {
        return false;
    }
    double massGap = 0;
    int ptmIndex = -1;
    ppm_value_ptr ppmValuePtr = std::make_shared<ppm_value>(minPPMValue);
    for (int i = startIndex; i >= 0; i--) {
        massGap = precursorMass - theoMassVec[i];
        if (massGap < minPtmMass) {
            break;
        }
        ptmIndex = mylib::data_stream::txpb_binarey_search_ex(ptmCMassVec, ptmCMassVec.size(), massGap);
        if (ptmIndex == -1) {
            continue;
        }
        if (ppmValuePtr->calculate_ppm_value(ptmCMassVec[ptmIndex], massGap)) {
            ptmCMass = ptmCMassVec[ptmIndex];
            return true;
        }
    }
//    delete ppmValuePtr;
    return false;
}

vector<filter_protein> ptm_filter::ptmFilterProcess(const vector<msalign_ptr> &msVec, const vector<protein_ptr> &prVec,
                                                    prsm::modify_ptr modPtr, prsm::config_argument_ptr configArgumentPtr) {
    vector<double> ptmComMassVec = modPtr->modifyMassSeq;
    vector<filter_protein> proFilterVec;
    for (msalign_ptr msalignPtr : msVec) {
        for (protein_ptr proteinPtr : prVec) {
            double ptmCMass = 0.0;
            if (isPtmConMass(proteinPtr->getTheoryMassCBy(), ptmComMassVec,
                             msalignPtr->getPrecursorMass(), ptmCMass,
                             configArgumentPtr->ptmMassMax, configArgumentPtr->ptmMassMin,
                             configArgumentPtr->min_ppm_value)) {
                filter_protein filterProtein(proteinPtr, ptmCMass);
                proFilterVec.push_back(filterProtein);
            }
        }
    }
    return proFilterVec;
    //    ofstream out("ptm_filter_protein.txt");
    //    for (filter_protein filterProtein : proFilterVec) {
    //        out << filterProtein.getProteinPtr()->getProteinTitle()
    //        <<" " << filterProtein.getPtmMassGap() << endl ;
    //    }
    //    out.close() ;
}