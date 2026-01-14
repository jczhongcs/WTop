//
// Created by Administrator on 2023/3/23.
//

#ifndef WTOP_PRSM_VARIBLE_PTMS_HPP
#define WTOP_PRSM_VARIBLE_PTMS_HPP

#include "prsm_unknown_ptms.hpp"

class prsm_varible_ptms : public prsm_unknown_ptms{
public:

    prsm_varible_ptms() {

        if (this->getPrsmAlignmentFilterPtr() != NULL) {
            delete this->getPrsmAlignmentFilterPtr();
            this->setPrsmAlignmentFilterPtr(new prsm_alignment_filter_var());
        }

    }

    void search_ptms_(vector<msalign_ptr> & msalign_ptrs,
                      prsm::modify_ptr modifyPtr,
                      prsm::config_argument_ptr configArgumentPtr,
                      const string & outFileName);

    void searchPtms_V2(vector<msalign_ptr> &msalign_ptrs,
                       prsm::modify_ptr modifyPtr,
                       prsm::config_argument_ptr configArgumentPtr,
                       const string &outFileName);

    void terminal_trunction(protein_ptr proteinPtr, msalign_ptr msalignPtr, prsm::config_argument_ptr prsmArg,
                            const vector<double> &ptmMass,
                            const double &precursorMass, map<int, int> &cutMap);

    int cur_ms_num = 1;

};
typedef prsm_varible_ptms* prsm_varible_ptms_ptr ;

#endif //WTOP_PRSM_VARIBLE_PTMS_HPP
