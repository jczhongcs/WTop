//
// Created by Administrator on 2023/3/30.
//

#ifndef WTOP_FILTER_TAG_MSG_HPP
#define WTOP_FILTER_TAG_MSG_HPP

#include <iostream>
#include <string>
#include <vector>
#include <map>

using namespace std;

class filter_tag_msg {
public:
    void insertTagMap(int fromIndex, string tag) {
        tagMap.insert(make_pair(fromIndex, tag));
    }

    void insertMassErrorMap(string tag, double massError) {
        massErrorMap.insert(make_pair(tag, massError));
    }

    void setCurIndex(int index) {
        curIndex = index;
    }

    int getCurIndex() const {
        return curIndex;
    }

    const multimap<int, string> &getTagMap() const {
        return tagMap;
    }

    const map<string, double> &getMassErrorMap() const {
        return massErrorMap;
    }

    void filterTagByProteinSeq(string &proteinSeq, string &revProSeq);

    void outProteinSeqMap(vector<double> &monoMass);

    void outTheMapMsg(vector<double> &monoMass);

    int curIndex = 0;
    multimap<int, string> tagMap;
    map<string, double> massErrorMap;
    multimap<int, string> seqTagMapN, seqTagMapC;

private:
};
typedef filter_tag_msg *filter_tag_msg_ptr;

#endif //WTOP_FILTER_TAG_MSG_HPP
