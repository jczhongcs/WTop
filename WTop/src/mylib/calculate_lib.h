//
// Created by 碎雨粘霓裳 on 2022/9/4.
//

#ifndef DTW_TOPZ_SPUTILS_H
#define DTW_TOPZ_SPUTILS_H

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

using namespace std;

namespace mylib {
    namespace calculate_lib {
        /**
         * 判断理论和precursormass是否相等
         * @param tMass
         * @param pMass
         * @param hMass
         */
        bool calculate_ppm_value_is_low(double tMass, double pMass, double hMass, double ppm) ;
    }
}


#endif //DTW_TOPZ_SPUTILS_H
