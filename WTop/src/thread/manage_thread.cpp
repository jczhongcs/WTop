//
// Created by Administrator on 2023/3/18.
//

#include "manage_thread.hpp"

void manage_thread::start_prsm_threads(std::vector<string> &prsm_out_file)
{
    int thread_num = initConfigPtr->read_config_ptr->getThreadNum();
    std::string protein_file = initConfigPtr->read_config_ptr->getProteinFile();
    std::set<string> msalign_file_set = initConfigPtr->read_config_ptr->getMsalignFileSet();
    std::string var_ptm_combination_file = initConfigPtr->read_config_ptr->getVarPtmCombinationOutFile();

    cout << "thread number - " << thread_num << endl;
    for (string msalign_file : msalign_file_set)
    {
        if (thread_num > 0)
        {
            mylib::thread_management::multi_thread_startor(protein_file, msalign_file,
                                                           var_ptm_combination_file,
                                                           thread_num, prsm_out_file);
        }
        else
        {
            cout << "thread number wrong" << endl;
        }
    }

    return;
}