//
// Created by wenzhong on 2023/3/21.
//

#ifndef WTOP_MSALIGN_TAG_HPP
#define WTOP_MSALIGN_TAG_HPP

#include "msalign.hpp"

using namespace std;

class msalign_tag {
public:
    msalign_tag(){

    };

    msalign_tag(const msalign_ptr msalignPtr, map<char, double> &acidMap) {
        getTagByMonoMassV2(msalignPtr, acidMap);
    }

    void getTagByMonoMass(const vector<double> &monoMass, map<char, double> &cMass,
                          set<string> &set, map<char, char> &replaceMap);

    void replaceCMap(map<char, double> &cMass, map<char, char> &rMap, map<int, char> &iCMap);

    void getTagByMonoAndIntChar(vector<double> &monoRaw, map<int, char> &massChar, set<string> &set);

    void getTagByMonoAndIntCharV2(vector<double> &monoRaw,
                                  map<int, char> &standAcidMap,
                                  set<string> &tagSet,
                                  map<char, double> &acidMap,
                                  vector<double> &acidMassVec);

    void getTagByMonoMassV2(const msalign_ptr &msalignPtr, map<char, double> &acidMap);

    vector<string> getTagsVec() {
        return tagsVec;
    }

    map<char, char> getAcidReplaceMap() {
        return acidReplaceMap;
    }

    map<int, filter_tag_msg_ptr> getTagMap() {
        return tagMap;
    }

    ~msalign_tag() {
        for (auto it = tagMap.begin(); it != tagMap.end(); ++it) {
            delete it->second;
        }
    }

    map<string, double> tagErrorMap;

    map<char, char> acidReplaceMap;      // acid change map (raw acid , change acid)
    map<int, filter_tag_msg_ptr> tagMap; // ms tag msg map
    vector<string> tagsVec;              // tagVec sort by taglength
    double tagMassGap = 0.1;

    int minTagLen = 2;

    double H = 1.007276;

private:
};
typedef msalign_tag *msalign_tag_ptr;

#endif // WTOP_MSALIGN_TAG_HPP
