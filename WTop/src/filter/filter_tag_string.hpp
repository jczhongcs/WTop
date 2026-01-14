//
// Created by wenzhong on 2023/4/13.
//

#ifndef WTOP_FILTER_TAG_STRING_HPP
#define WTOP_FILTER_TAG_STRING_HPP

#include <vector>
#include <iostream>
#include <map>
#include <string>

using namespace std ;

class filter_tag_string {
public:

private:
    int index ;
    string tag ;
    int matchSeqIndex ;

    map<int,string> tagsMap ; // match tag
};


#endif //WTOP_FILTER_TAG_STRING_HPP
