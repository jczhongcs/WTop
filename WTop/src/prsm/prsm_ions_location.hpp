//
// Created by wenzhong on 2023/3/22.
//

#ifndef WTOP_PRSM_IONS_LOCATION_HPP
#define WTOP_PRSM_IONS_LOCATION_HPP

#include "../ms/msalign.hpp"
#include "../protein/protein.hpp"
#include "../merge/modify.hpp"
#include "../merge/config_argument.hpp"
#include "../merge/temp_prsm_argument.hpp"
#include "../mylib//etd_process.hpp"

class prsm_ions_location {
public:

    void ions_location_match_N(vector<node> &n_ions, vector<node> &c_ions,
                               vector<node> &merg_ions_n,vector<node> &merg_ions_c,
                               vector<modi> &mod, const vector<double> &modifyTable, int lenth, double modifier_mass);

    void ions_location_match_C(vector<node> &n_ions, vector<node> &c_ions,
                               vector<node> &merg_ions_n,vector<node> &merg_ions_c,
                          vector<modi> &mod, const vector<double> &modifyTable, int lenth, double modifier_mass);

    void ions_location_match_N_v2(vector<node> &n_ions, vector<node> &c_ions, vector<node> &merg_ions_n,
                                  vector<node> &merg_ions_c, vector<modi> &mod, const vector<double> &modifyTable,
                                  int lenth, double modifier_mass);

    void ions_location_match_C_v2(vector<node> &n_ions, vector<node> &c_ions, vector<node> &merg_ions_n,
                                  vector<node> &merg_ions_c,
                                  vector<modi> &mod, const vector<double> &modifyTable, int lenth,
                                  double modifier_mass);

    void
    confirm_final_ions(vector<node> &n_ions_n, vector<node> &n_ions_c, vector<modi> &n_ptms, vector<node> &c_ions_n,
                       vector<node> &c_ions_c, vector<modi> &c_ptms, prsm::temp_prsm_argument_ptr tempPrsmArgumentPtr);
};
typedef prsm_ions_location* prsm_ions_location_ptr;

#endif //WTOP_PRSM_IONS_LOCATION_HPP
