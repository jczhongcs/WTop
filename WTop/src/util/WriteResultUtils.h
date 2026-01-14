//
// Created by wenzhong on 2023/1/5.
//

#ifndef DTW_WTOP_WRITERESULTUTILS_H
#define DTW_WTOP_WRITERESULTUTILS_H

#include <iostream>

#include "Argument.h"
#include "Service.h"

void processResult(ofstream &out, vector<mylib::Argument> &aContainer);

//输出初始文件
void WriteArgumentResult(ofstream &out, vector<mylib::Argument> &aContainer);

vector<mylib::Argument> processOneFileResult(string inFileName , string outFileName, map<string,string> & ssMap);


#endif //DTW_WTOP_WRITERESULTUTILS_H
