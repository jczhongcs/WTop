//
// Created by Administrator on 2023/3/31.
//

#ifndef WTOP_PTM_FILTER_HPP
#define WTOP_PTM_FILTER_HPP

#include "../merge/modify.hpp"
#include "../ms/msalign.hpp"
#include "../protein/protein.hpp"
#include "../protein/filter_protein.h"
#include "../merge/config_argument.hpp"

using namespace std;

class ptm_filter {
public:
    vector<filter_protein> ptmFilterProcess(const vector<msalign_ptr> & msVec, const vector<protein_ptr> & prVec,
                          prsm::modify_ptr modPtr,prsm::config_argument_ptr configArgumentPtr);
};
typedef ptm_filter* ptm_filter_ptr ;

#endif //WTOP_PTM_FILTER_HPP
