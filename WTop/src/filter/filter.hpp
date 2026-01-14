//
// Created by wenzhong on 2023/3/21.
//

#ifndef WTOP_FILTER_HPP
#define WTOP_FILTER_HPP

#include <thread>
#include <algorithm>
#include <chrono>

#include "../ms/msalign_tag.hpp"
#include "../ms/complement_ions.hpp"

#include "../merge/config_argument.hpp"

#include "ms_protein_filter.hpp"

#include "ptm_filter.hpp"
#include "../util/string_utils.h"


using namespace std::chrono; // 使用命名空间

class filter {
public:
    void filtrationProcess(const vector<msalign_ptr> &msContainer,
                           const vector<protein_ptr> &prContainer,
                           const prsm::modify_ptr modifyPtr,
                           const prsm::config_argument_ptr configArgumentPtr);

    void tagFilterProcess(const vector<msalign_ptr> &msContainer,
                          const vector<protein_ptr> &prContainer,
                          const prsm::modify_ptr modifyPtr,
                          const prsm::config_argument_ptr configArgumentPtr);
    void getToppicFilterResult(const vector<msalign_ptr> &msContainer,
                               const vector<protein_ptr> &prContainer,
                               const prsm::modify_ptr modifyPtr,
                               const prsm::config_argument_ptr configArgumentPtr);

private:
};
typedef filter *filter_ptr;

#endif //WTOP_FILTER_HPP
