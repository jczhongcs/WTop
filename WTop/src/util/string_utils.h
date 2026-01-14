//
// Created by wenzhong on 2022/11/15.
//

#ifndef DTW_TOPZ_STRUTILS_H
#define DTW_TOPZ_STRUTILS_H

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <vector>

#include <math.h>
#include <sstream>

#include "Argument.h"       //用于处理输出文件
#include "Service.h"        //库
#include "WriteResultUtils.h"

using namespace std ;

namespace utils_string_util{

    void processProteinFom(string &rawSeq,vector<string> & ptms,map<char,string> & charPTMsMap,map<int,int> &startEnd,string &pfSeq);

    void processProteinFomUnkMass(string &rawSeq,vector<string> & ptms,map<char,string> ptmsMap,map<int,int> &startEnd,string &pfSeq,string unk_mass);

    void processCutLocation(string &seq,int cutN,int cutC);

    void txtFileToCsvFile(vector<string> & outFile, map<string,string> & ssMap);

    string convertToString(double d) ;

    void Stringsplit(string str, const char split,vector<string>& rst);

    vector<string> split(const string& str, const string& delim);

    string cutProName(const string& str);

    string getProID(const string& str);


}
//处理输出的蛋白质形


#endif //DTW_TOPZ_STRUTILS_H
