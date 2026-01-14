//
// Created by wenzhong on 2023/1/5.
//

#ifndef DTW_WTOP_SERVICE_H
#define DTW_WTOP_SERVICE_H


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <algorithm>

#include "Argument.h"

using namespace std;

namespace mylib {
    class Service {
    public:
        //依据质量排序
        static void MassSortArgument(vector<mylib::Argument> &aContainer);

        //依据分数排序
        static void ScoreSortArgument(vector<mylib::Argument> &aContainer);

        static int calculateFdr(vector<mylib::Argument> &aContainer);

        //选出fdr < value 并且 不含DECOY 的进入bContainer
        static void fdrChoice(vector<mylib::Argument> &aContainer, vector<mylib::Argument> &bContainer,
                              double value, int minNub, double filterArg);

        static void getProteinForm(vector<mylib::Argument> &bContainer,
                                   vector<mylib::Argument> &cContainer);

        static void getUnknownProteinForm(vector<mylib::Argument> &bContainer, vector<mylib::Argument> &cContainer);
        static map<string,string> getScansScoreMap(vector<mylib::Argument> &Container,
                                                   map<string,mylib::Argument> &m);
    };
}


#endif //DTW_WTOP_SERVICE_H
