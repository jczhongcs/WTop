//
// Created by 碎雨粘霓裳 on 2022/4/29.
//

#ifndef DTW_TOPZ_UNKNOW_PTM_SEARCH_H
#define DTW_TOPZ_UNKNOW_PTM_SEARCH_H

#include <iostream>
#include <vector>
#include <algorithm>     // 算法头文件，提供迭代器
#include <fstream>
#include <iomanip>
#include <map>           // STL
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <string>
#include <unistd.h>
#include <cmath>

#include "wtop_prsm_processor.hpp"
#include "unknown_ptm_utils.h"
#include "config_argument.hpp"

#include "../protein/protein.hpp"

namespace prsm {
    class mspr{
    public:
        msalign_processor_ptr it;
        vector<protein_processor_ptr> prContainer ;   //保存过滤后的蛋白
        int match_size = 0 ;
    };
    typedef mspr* msprP ;

    void test_filter_result(msalign_processor_ptr it) ;
    void one_ptm_filter(vector<msalign_processor_ptr> &msContainer, vector<protein_processor_ptr> &prContainer,
                        vector<double> &ptmMassSeq, vector<msprP> &unknow_ptm_filter,
                        modify_ptr one_mod_ptr, const string & outFileName, config_argument_ptr prsmArg);

    void unknow_ptms_prsm_process_CID(prsm::modify_ptr modptr, protein_processor_ptr pro, msalign_processor_ptr spact, config_argument_ptr prsmArg,
                                      map<int,int> & iiMap, const string &variablePTMsFileOutPath);

    void unknow_ptms_prsm_process_ETD(prsm::modify_ptr mp, protein_processor_ptr pro, msalign_processor_ptr spact, config_argument_ptr prsmArg,
                                      map<int,int> & iiMap, const string &variablePTMsFileOutPath) ;

    void initial_Unknown_ptmTable(temp_prsm_argument_ptr arg, vector<double> &modifierTable,
                                  map<double,string> &mmap,
                                  double ptmMass, string ptm, double bestUnknown) ;

    void prsm_process_unknown(prsm::modify_ptr modptr, protein_processor_ptr pro, msalign_processor_ptr spact, temp_prsm_argument_ptr arg,
                              vector<double> &mass, config_argument_ptr prsmArg);

    void writeResult(msalign_processor_ptr msIt, string Result_file, modify_ptr modptr);

    double judge_Unknown_Mass(double PtmMass, double prMass, double unk_mass,
                              temp_prsm_argument_ptr arg, msalign_processor_ptr spact, protein_processor_ptr pro);

    void makeTagByMonoMass(protein_processor_ptr prIt, msalign_processor_ptr msIt, map<char, char> &ccMap, set<string> &set);

//    void filtrationProcess(const vector<msalign_ptr> &msContainer,
//                           const vector<protein_ptr> &prContainer,
//                           const config_argument_ptr prsmArg);

    void one_ptm_search_process(modify_ptr one_mod_ptr, vector<msalign_processor_ptr> &msContainer, config_argument_ptr prsmArg,
                                const string &variablePTMsFileOutPath, string & outFileName);

    void mergeOnePTMPrSMsProcess(Modify *mp,
                                 protein_processor_ptr pro,
                                 msalign_processor_ptr spact,
                                 config_argument_ptr prsmArg,
                                 const string &variablePTMsFileOutPath);
}


#endif //DTW_TOPZ_UNKNOW_PTM_SEARCH_H
