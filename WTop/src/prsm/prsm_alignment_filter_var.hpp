//
// Created by Administrator on 2023/3/24.
//

#ifndef WTOP_PRSM_ALIGNMENT_FILTER_VAR_HPP
#define WTOP_PRSM_ALIGNMENT_FILTER_VAR_HPP

#include "prsm_alignment_filter.hpp"

class prsm_alignment_filter_var : public prsm_alignment_filter{
public:
    void index_ions_mass(vector<double> &sub, vector<double> &index,
                         vector<double> &modifier, double &modifier_mass);

};
typedef prsm_alignment_filter_var* prsm_alignment_filter_var_ptr ;

#endif //WTOP_PRSM_ALIGNMENT_FILTER_VAR_HPP
