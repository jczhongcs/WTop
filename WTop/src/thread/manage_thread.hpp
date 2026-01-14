//
// Created by Administrator on 2023/3/18.
//

#ifndef WTOP_MANAGE_THREAD_HPP
#define WTOP_MANAGE_THREAD_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <mutex>

#include "../ms/msalign.hpp"
#include "../mylib/thread_management.hpp"
#include "../config/init_config.hpp"

using namespace std;

class manage_thread {

public:
    void start_prsm_threads(std::vector<string> & prsm_out_file);

    manage_thread(init_config_Ptr initConfigPtr) { this->initConfigPtr = initConfigPtr;}

private:
    init_config_Ptr initConfigPtr ;
};
typedef manage_thread* manage_thread_Ptr ;


#endif //WTOP_MANAGE_THREAD_HPP
