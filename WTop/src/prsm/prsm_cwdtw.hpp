//
// Created by wenzhong on 2023/3/22.
//

#ifndef WTOP_PRSM_CWDTW_HPP
#define WTOP_PRSM_CWDTW_HPP

#include <vector>
#include <math.h>
#include "../mylib/fun.h"

using namespace std;
class prsm_cwdtw {
public:
    double cwdtw_algorithm(vector<double> &reference_1, vector<double> &peer_1,
                           vector<std::pair<long, long>> &alignment,
                           vector<double> &ref_zscore, vector<double> &peer_zscore);
};
typedef prsm_cwdtw* prsm_cwdtw_ptr;


#endif //WTOP_PRSM_CWDTW_HPP
