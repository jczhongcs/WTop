//
// Created by wenzhong on 2023/3/22.
//

#ifndef WTOP_PRSM_TERMINAL_TRUNCATION_HPP
#define WTOP_PRSM_TERMINAL_TRUNCATION_HPP

#include "candi_truncation.h"
#include "../merge/modify.hpp"
#include "../protein/protein.hpp"
#include "../ms/msalign.hpp"
#include "../merge/config_argument.hpp"
#include "../util/ppm_value.hpp"
#include "../mylib/speculate_lib.h"
#include "algorithm"

class prsm_terminal_truncation {
public:
    void judge_terminal_truncation_By(prsm::modify_ptr mp,
                                      protein_ptr proteinPtr,
                                      msalign_ptr msalignPtr,
                                      const prsm::config_argument_ptr &prsmArg,
                                      int &cut_n, int cut_c);

    int match_peaks_search_ppm(const vector<double> &mono_mass,
                               const vector<double> &theo_mass,
                               double para_ppm);

    int match_peaks_search_value(const vector<double> &mono_mass,
                                 const vector<double> &theo_mass,
                                 double value);

    void judge_terminal_truncation_Cz(prsm::modify_ptr mp, protein_ptr proteinPtr, msalign_ptr msalignPtr,
                                      prsm::config_argument_ptr const &prsmArg, int &cut_n, int cut_c);


    void test_NC_By(prsm::modify_ptr mp,protein_ptr proteinPtr,msalign_ptr msalignPtr,const prsm::config_argument_ptr &prsmArg,map<double,int> &modCutN, map<double,int> &modCutC,map<double,vector<double>> &modUnk);

    void test_NC_Cz(prsm::modify_ptr mp,protein_ptr proteinPtr,msalign_ptr msalignPtr,const prsm::config_argument_ptr &prsmArg,map<double,int> &modCutN, map<double,int> &modCutC,map<double,vector<double>> &modUnk);

    void locate_truncation_by(prsm::modify_ptr mp,protein_ptr proteinPtr,msalign_ptr msalignPtr,const prsm::config_argument_ptr &prsmArg,map<double,int> &modCutN, map<double,int> &modCutC,map<double,vector<double>> &modUnk);

    void locate_truncation_cz(prsm::modify_ptr mp,protein_ptr proteinPtr,msalign_ptr msalignPtr,const prsm::config_argument_ptr &prsmArg,map<double,int> &modCutN, map<double,int> &modCutC,map<double,vector<double>> &modUnk);

    void locate_truncation_by(prsm::modify_ptr mp,protein_ptr proteinPtr,msalign_ptr msalignPtr,const prsm::config_argument_ptr &prsmArg,map<double,candi_trucation_vec> &ptm_candiTru_map);

    void locate_truncation_cz(prsm::modify_ptr mp,protein_ptr proteinPtr,msalign_ptr msalignPtr,const prsm::config_argument_ptr &prsmArg,map<double,candi_trucation_vec> &ptm_candiTru_map);

};

typedef prsm_terminal_truncation* prsm_terminal_truncation_ptr;


#endif //WTOP_PRSM_TERMINAL_TRUNCATION_HPP
