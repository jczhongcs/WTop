//
// Created by Administrator on 2022/9/22.
//

#ifndef DTW_TOPZ_UNKNOWNPTMUTILS_H
#define DTW_TOPZ_UNKNOWNPTMUTILS_H

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


#include "../mylib/data_steam.hpp"
//#include "protein_processor.hpp"

using namespace std;

namespace unknown {
    namespace ptm_utils{
        //用来过滤时，获取N端的连续离子对。
        int getSeqIonsPair(vector<node> & t, vector<node> & m) ;

        //根据monoMass，快速制造tag
        void getTagByMonoMass(vector<double>& monoMass,map<char,double> & cMass,set<string> &set,map<char,char> &replaceMap);

        //从mono中寻找出互补离子
        void getComplementIonsMap(double precursorMass,const vector<double> &mono,map<int,int> &indexMap);

        //寻找到最大和最小截断
        void getCutMinAndMax(vector<double> & cT,double precursorMass,
                             double scopeValue,
                             int &cutMin,int &cutMax);


        //得到该蛋白质形式下ptm mass对应的unknown mass 用map保存
        void getUnknownMassByPtmMass(const vector<double> & cT,const vector<double> & ptmMass,
                                     double precursorMass,int location,map<double,double> &ptmUnknown,double scopeValue);

        //确定互补离子的C端离子
        bool getComplementIonsCTIndex(double m1,double m2,int &m1Index,int &m2Index,int cutMax,
                                      const vector<double> & cT,const vector<double> & nT,
                                      vector<int> &trueModIndex,
                                      const map<double,double> & ptmUnknown,double nTCutMass,double scopeValue);

//        void getConPro(vector<prsm::protein_processor_ptr> & p1, vector<prsm::protein_processor_ptr> & p2, vector<prsm::protein_processor_ptr> & con);
    }
}


#endif //DTW_TOPZ_UNKNOWNPTMUTILS_H
