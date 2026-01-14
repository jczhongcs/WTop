#ifndef THREAD_MANAGEMENT_HPP
#define THREAD_MANAGEMENT_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <mutex>
#include "data_steam.hpp"

#include "../merge/wtop_prsm_processor.hpp"
#include "../ms/msalign.hpp"
#include "../protein/protein.hpp"

using namespace std;


namespace mylib {

    namespace thread_management {

        static vector<vector<protein_ptr>> separate_protein_ptr ;
        static int init_protein_process_num = 0 ;
        static double protein_total_number = 0 ;

        static std::mutex mutex_;

        void multi_thread_startor(string &protein_file, string &msalign_file,
                                  string &var_ptm_combination_file,
                                  int thread_num, vector<string> & prsm_out_file);

        void separate_msalign_file(string &msalign_file,
                                   vector<string> &separate_msalign_file, int separate_num);

        void separate_protein_file(string & protein_file,
                                   vector<string> & separate_protein_file_container, int separate_num);

        void init_multi_protein_thread_start(const string &separate_protein_file, const int separate_container_index);

        void init_protein_process_rate_thread();

        void thread_wtop_porcess_start(const string &spactrue_file_name,
                                       const string &variablePTMsFileOutPath, const string &ptmOutFileName,
                                       const string &ptmOutFileName2,
                                       const vector<protein_ptr> &protein_ptr_container);
    }
}
#endif