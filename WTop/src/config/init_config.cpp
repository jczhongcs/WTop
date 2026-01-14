//
// Created by wenzhong on 2023/3/17.
//

#include "init_config.hpp"

void init_config::init_combination_ptm()
{
    ptm_combination_sharePtr ptmCombinationPtr = std::make_shared<ptm_combination>();
    ptmCombinationPtr->init_var_ptm_combination(read_config_ptr->getPtmFile(),
                                                read_config_ptr->getVarPtmCombinationOutFile(),
                                                6,500) ;
//    delete ptmCombinationPtr ;
}