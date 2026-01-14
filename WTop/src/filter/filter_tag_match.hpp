//
// Created by Administrator on 2023/4/12.
//

#ifndef WTOP_FILTER_TAG_MATCH_HPP
#define WTOP_FILTER_TAG_MATCH_HPP

#include <iostream>
#include <map>
#include <vector>

using namespace std;

class filter_tag_match {

public:

     map<int, double> &getIndexMsMap()  {
        return indexMsMap;
    }

     multimap<int, string> &getIndexTagMMap()  {
        return indexTagMMap;
    }

     map<string, vector<int>> &getIndexTagBeginMap()  {
        return indexTagBeginMap;
    }

    void insertIndexTagBeginMap(string index, vector<int> &beginVec) {
        indexTagBeginMap.insert(make_pair(index,beginVec));
    }

    void insertIndexTagMMap(int mtIndex,string tag) {
        indexTagMMap.insert(make_pair(mtIndex,tag));
    }

    void insertIndexMsMap(int index, double monoMass) {
        indexMsMap.insert(make_pair(index,monoMass));
    }

    void outMsgTest(){
        for (auto it = indexTagMMap.begin(); it != indexTagMMap.end(); ++it) {
            vector<int> tagMatchSeq = indexTagBeginMap[to_string(it->first) + it->second] ;
            cout << indexMsMap[it->first] << " -> " << it->first
            << " -> " << it->second ;
            for (int tagMatchIndex : tagMatchSeq) {
                cout << " -> " << tagMatchIndex ;
            }
            cout << endl ;
        }
    }
private:
    map<int,double> indexMsMap ;
    multimap<int,string> indexTagMMap ;
    map<string,vector<int>> indexTagBeginMap;  // mono index + tag -> tag match index Vec
};
typedef filter_tag_match* filter_tag_match_ptr;

#endif //WTOP_FILTER_TAG_MATCH_HPP
