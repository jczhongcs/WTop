#ifndef PRSM_PROCESS_HPP
#define PRSM_PROCESS_HPP

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <algorithm>
#include <thread>
#include <mutex>
#include "temp_prsm_argument.hpp"
#include "cluster_msalign_protein.hpp"
#include "unknown_ptm_search.h"
#include "processor_management.h"
#include "config_argument.hpp"

#include "../mylib/io.h"
#include "../mylib/speculate_lib.h"
#include "../mylib/data_steam.hpp"
#include "../mylib/etd_process.hpp"
#include "../mylib/terminal_truncation.h"
#include "../mylib/calculate_lib.h"
#include "../protein/protein.hpp"
#include "../ms/msalign.hpp"

#include "../prsm/prsm_varible_ptms.hpp"

#include "../filter/filter.hpp"
#include "../filter/ptm_filter.hpp"

using namespace std;

namespace prsm {

    class WTop_Process {
        public:

        void wtop_process(const string &msalign_file,
                          const string & var_ptm_file,
                          const string & var_ptm_result_file, const string & unknown_ptm_result_file,
                          const vector<protein_ptr> & protein_container);

        void prsm_process_CID(Modify * modptr, protein_processor_ptr pro,
                              msalign_processor_ptr spact,
                              temp_prsm_argument_ptr arg, config_argument_ptr prsmArg);

        void prsm_process_ETD(Modify * modptr,
                              protein_processor_ptr pro,
                              msalign_processor_ptr spact,
                              temp_prsm_argument_ptr temp_arg, config_argument_ptr config_arg);


        void process_var_ptm_prsm(vector<msalign_processor_ptr> & msalign_container,
                                  prsm::modify_ptr modptr,
                                  prsm::config_argument_ptr & config_arg,
                                  const string &out_file_name);

        void process_one_ptm_prsm(modify_ptr one_mod_ptr,
                                  vector<msalign_processor_ptr> &msalign_container,
                                  config_argument_ptr config_arg,
                                  const string &variablePTMsFileOutPath, const string & outFileName);



        vector<double> modifier_table;

        vector<double> one_ptm_table;

//        vector<protein_processor_ptr> protein_container;
        //参数设置
//        double h20 = 18.01056 ;     //h20
//        double core = 1.9919 ;
//        double ptmMass = 500 ; // 输入修饰质量
        char pI ;    //记录n端标识

        int test_code_out = 0 ; // 1 输出测试
    };
    typedef std::shared_ptr<WTop_Process> wtop_process_share_ptr;
    typedef WTop_Process * wtop_process_ptr;
}
#endif