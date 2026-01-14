//
// Created by wenzhong on 2023/3/17.
//

#ifndef WTOP_INIT_CONFIG_HPP
#define WTOP_INIT_CONFIG_HPP

#include "read_config.hpp"
#include "../ptm/ptm_combination.hpp"

class init_config {

public:
    read_config_Ptr read_config_ptr ;

    ~init_config() {
        delete read_config_ptr ;
    }

    init_config() {
        read_config_ptr = new read_config();
    }

    void init_combination_ptm();
};
typedef shared_ptr<init_config> init_config_sharePtr;
typedef init_config* init_config_Ptr;


#endif //WTOP_INIT_CONFIG_HPP
