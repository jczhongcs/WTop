//
// Created by LYC on 2025/3/29.
//

#ifndef WTOP_XMLUTIL_H
#define WTOP_XMLUTIL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "tinyxml.h"
#include "string"
#include <cstdio>
#include <unistd.h>
#include <filesystem>

namespace fs = std::filesystem;
using namespace std;
void CombineAllXmlToCsv(vector<string> evalueFileVec);
void outXmlForWTopMLS();

#endif //WTOP_XMLUTIL_H
