//
// Created by wenzhong on 2023/3/22.
//

#ifndef WTOP_PRSM_SCOPE_ALIGNMENT_HPP
#define WTOP_PRSM_SCOPE_ALIGNMENT_HPP

#include <vector>
#include <set>
#include <math.h>
#include <algorithm>

#include "../mylib/data_steam.hpp"

using namespace std;
class prsm_scope_alignment {
public:
    void scopeAlignment(const vector<double> &reference, const vector<double> &peer,
                        vector<std::pair<long, long>> &alignment,
                        vector<double> &modifier, const double &modifier_mass, double maxMass);
};
typedef prsm_scope_alignment* prsm_scope_alignment_ptr;


#endif //WTOP_PRSM_SCOPE_ALIGNMENT_HPP
