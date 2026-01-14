//
// Created by wenzhong on 2023/3/17.
//

#ifndef WTOP_PTM_COMBINATION_HPP
#define WTOP_PTM_COMBINATION_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
#include <map>
#include <fstream>
#include <iomanip>
#include <memory>

#include "../util/string_utils.h"

using namespace std;
class ptm_combination {
public:
    void init_var_ptm_combination(const string & variableFileName, string outPath, int maxSize, double maxMass);

};
typedef shared_ptr<ptm_combination> ptm_combination_sharePtr;
typedef ptm_combination* ptm_combination_Ptr ;


#endif //WTOP_PTM_COMBINATION_HPP
