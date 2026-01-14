//
// Created by wenzhong on 2023/3/15.
//

#ifndef DTW_WTOP_PROCESSOR_MANAGEMENT_H
#define DTW_WTOP_PROCESSOR_MANAGEMENT_H

#include <vector>
#include <iostream>

#include "msalign_processor.hpp"
#include "protein_processor.hpp"
#include "config_argument.hpp"

#include "../util/string_utils.h"
#include "../ms/msalign.hpp"

namespace prsm {
    namespace processor_management {

        void save_best_prsm_information(protein_processor_ptr protein,
                                        msalign_processor_ptr msalign,
                                        temp_prsm_argument_ptr temp_arg);

        void initializerPtmsPtr(prsm::modify_ptr ptmsPtr, prsm::modify_ptr one_mod_ptr,vector<double> &one_ptm_table);

        void init_msalign_container(const string & ms_file_name,
                                    vector<msalign_ptr> &ms_ptr_container);

        void init_protein_container(const string & protein_file,
                                    vector<prsm::protein_processor_ptr> & prContainer);

        void writeResultNews(msalign_processor_ptr msIt,
                             string Result_file, modify_ptr modptr);

        void test_out(char pI, temp_prsm_argument_ptr arg, msalign_processor_ptr spact);

        void writeIonsMessage(msalign_processor_ptr msIt,
                              config_argument_ptr configArgumentPtr,
                              string msSrcPath);


    }
}
#endif //DTW_WTOP_PROCESSOR_MANAGEMENT_H
