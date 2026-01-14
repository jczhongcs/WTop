//
// Created by wenzhong on 2023/3/21.
//

#include "msalign_tag.hpp"
#include <set>
// 质量相近的氨基酸进行统一替换
void msalign_tag::replaceCMap(map<char, double> &cMass, map<char, char> &rMap, map<int, char> &iCMap) {
    vector<double> pMass;
    // 先将质量相同的替换
    map<double, int> newCMap;
    for (auto it = cMass.begin(); it != cMass.end(); ++it) {
        if (!newCMap.count(it->second))
            newCMap.insert(make_pair(it->second, it->first));
        else {
            char c2 = newCMap[it->second];
            rMap.insert(make_pair(it->first, c2));
        }
    }
    //    for( auto it = newCMap.begin();it != newCMap.end(); ++it) {
    //        cout <<" first = " << it->first<< " second = " << it->second <<endl ;
    //    }
    //    cout  <<"=================================" << endl;

    for (auto it = newCMap.begin(); it != newCMap.end(); ++it) {
        int temp = it->first;
        char c = it->second;
        if (iCMap.count(temp) || iCMap.count(temp - 1) || iCMap.count(temp + 1)) {
            char c2;
            if (iCMap.count(temp)) {
                c2 = iCMap[temp];
            } else if (iCMap.count(temp - 1)) {
                c2 = iCMap[temp - 1];
            } else {
                c2 = iCMap[temp + 1];
            }
            rMap.insert(make_pair(c, c2));
            if (!iCMap.count(temp)) {
                iCMap.insert(make_pair(temp, c2));
            }
            if (!iCMap.count(temp - 1)) {
                iCMap.insert(make_pair(temp - 1, c2));
            }
            if (!iCMap.count(temp + 1)) {
                iCMap.insert(make_pair(temp + 1, c2));
            }
        } else {
            iCMap.insert(make_pair(temp - 1, c));
            iCMap.insert(make_pair(temp, c));
            iCMap.insert(make_pair(temp + 1, c));
        }
    }
    // 测试输出
    //    for( auto it = rMap.begin();it != rMap.end(); ++it) {
    //        cout <<" first = " << it->first<< " second = " << it->second <<endl ;
    //    }
    //    cout <<"=================================" << endl;
    //    for( auto it = iCMap.begin();it != iCMap.end(); ++it) {
    //        cout <<" first = " << it->first<< " second = " << it->second <<endl ;
    //    }
}

/**
 * 根据mono的质量寻找tag ，并保存到set中
 * @param mono
 * @param massChar
 */
void msalign_tag::getTagByMonoAndIntChar(vector<double> &monoRaw, map<int, char> &massChar, set<string> &set) {
    map<int, string> intStringMap; // 前一个为下标位置，后一个为下标对应的tag;
    // 针对可能存在的多个，将string改进为vector<string>
    map<int, vector<string>> intVectorMap; // 前一个为下标位置，后一个为下标对应的Vector tag;
    vector<double> mono;
    for (int i = 0; i < monoRaw.size(); ++i) {
        if (i > 1 && (int)monoRaw[i] == (int)monoRaw[i - 1]) {
            continue;
        }
        mono.push_back(monoRaw[i - 1]);
    }
    // 进行tag的判定
    int minMass = massChar.begin()->first;
    int maxMass = 0;
    for (auto it = massChar.begin(); it != massChar.end(); ++it) {
        if (it->first > maxMass) {
            maxMass = it->first;
        }
    }
    for (int i = 1; i < mono.size(); ++i) {
        for (int j = i - 1; j >= 0; j--) {
            int temp = mono[i] - mono[j];
            string ts;
            if (!(temp >= minMass && temp <= maxMass)) {
                continue;
            }
            if (massChar.count(temp)) {
                char c = massChar[temp];     // 得到当前氨基酸
                ts += c;                     // 加上当前氨基酸
                if (intStringMap.count(j)) { // 如果被减的值，拥有了氨基酸序列,则进行拼接
                    string s = intStringMap[j];
                    ts += s;
                    intStringMap.insert(make_pair(i, ts));
                    // cout << " mono[" << i << "] = " << mono[i] << " mono[" << j << "]=  " << mono[j]
                    //      << " temp = " << temp
                    //      << " ts = " << ts
                    //      << endl;
                } else { // 直接加入
                    intStringMap.insert(make_pair(i, ts));
                    // cout << " mono[" << i << "] = " << mono[i] << " mono[" << j << "]=  " << mono[j]
                    //      << " temp = " << temp
                    //      << " ts = " << ts
                    //      << endl;
                }
                break;
            }
        }
    }

    for (auto it = intStringMap.begin(); it != intStringMap.end(); ++it) {
        set.insert(it->second);
    }

    // for (string s : set) {
    //     cout << " s = " << s << endl;
    // }
}

void msalign_tag::getTagByMonoAndIntCharV2(vector<double> &monoRaw,
                                           map<int, char> &standAcidMap,
                                           set<string> &tagsSet,
                                           map<char, double> &acidMap,
                                           vector<double> &acidMassVec) {
    double minAcidMass = standAcidMap.begin()->first;
    double maxAcidMass = 0;
    for (map<int, char>::iterator it = standAcidMap.begin(); it != standAcidMap.end(); ++it) {
        maxAcidMass = max((double)it->first, maxAcidMass);
    }

    // start to calculate tag
    for (int i = 1; i < monoRaw.size(); ++i) {
        filter_tag_msg_ptr filterTagMsgPtr = new filter_tag_msg();
        std::set<string> tagSet;

        for (int j = i - 1; j >= 0; --j) {
            int acidMass = monoRaw[i] - monoRaw[j]; // round((monoRaw[i] - monoRaw[j]) * 100) / 100;
            if (acidMass > maxAcidMass) {           //(acidMass < minAcidMass || acidMass > maxAcidMass) {
                break;
            }
            if (!standAcidMap.count(acidMass)) {
                continue;
            }
            char acid = standAcidMap[acidMass];
            double massClose = mylib::data_stream::findClose(acidMassVec, acidMassVec.size(), monoRaw[i] - monoRaw[j]);
            double massError = fabs(monoRaw[i] - monoRaw[j] - massClose);
            if (tagMap.count(j)) {
                multimap<int, string> iTagMap = tagMap[j]->getTagMap(); // key -> i , value -> tag
                map<string, double> massErrorMap = tagMap[j]->getMassErrorMap();
                for (auto it = iTagMap.begin(); it != iTagMap.end(); ++it) {
                    string acidStr;
                    acidStr += it->second;
                    double avgMassError = (massError + massErrorMap[to_string(it->first) + it->second]) / (double)(acidStr.size() + 1);
                    acidStr.push_back(acid);
                    if (tagSet.count(acidStr) != 1) {
                        filterTagMsgPtr->tagMap.insert(make_pair(j, acidStr));
                        filterTagMsgPtr->massErrorMap.insert(
                            make_pair(to_string(j) + acidStr, massError + massErrorMap[to_string(it->first) + it->second])); //insertMassErrorMap(to_string(j) + acidStr, massError + massErrorMap[to_string(it->first) + it->second]);
                        tagSet.insert(acidStr);
                    }
                }
                string acidStr;
                acidStr.push_back(acid);
                filterTagMsgPtr->tagMap.insert(make_pair(j, acidStr));                              //insertTagMap(j, acidStr);
                filterTagMsgPtr->massErrorMap.insert(make_pair(to_string(j) + acidStr, massError)); //insertMassErrorMap(to_string(j) + acidStr, massError);
            } else {
                string acidStr;
                acidStr.push_back(acid);
                if (tagSet.count(acidStr) != 1) {
                    filterTagMsgPtr->tagMap.insert(make_pair(j, acidStr));
                    filterTagMsgPtr->massErrorMap.insert(make_pair(to_string(j) + acidStr, massError)); //insertMassErrorMap(to_string(j) + acidStr, massError);
                    tagSet.insert(acidStr);
                }
            }
        }

        if (filterTagMsgPtr->getTagMap().size() != 0) {
            filterTagMsgPtr->setCurIndex(i);
            tagMap.insert(make_pair(i, filterTagMsgPtr));
            // filterTagMsgPtr->outTheMapMsg(monoRaw);
        } else {
            delete filterTagMsgPtr;
        }
    }
}

void msalign_tag::getTagByMonoMass(const vector<double> &monoMass, map<char, double> &cMass,
                                   set<string> &set, map<char, char> &replaceMap) {
    vector<double> mono = monoMass;
    vector<double> pMass;
    for (auto it = cMass.begin(); it != cMass.end(); ++it) {
        pMass.push_back(it->second);
    }
    sort(mono.begin(), mono.end());
    sort(pMass.begin(), pMass.end());

    // map<char,char> replaceMap; //前一个为被替换的氨基酸，后一个为替换后的氨基酸
    map<int, char> massChar; // 前一个为标准化之后的数据，后一个为对应的氨基酸
    // 质量相近的氨基酸进行统一替换,并得到替换的氨基酸map和一个int，char的map
    msalign_tag::replaceCMap(cMass, replaceMap, massChar); // 得到替换的氨基酸map
    // set<string> set ; //保存所有的tag
    getTagByMonoAndIntChar(mono, massChar, set); // 得到tag

    //    getTagByMonoAndIntCharV2(mono,massChar,set);
}

bool isDuplicate(double a, double b) {
    return std::abs(a - b) < 1.5;
}

vector<int> removeDupMonoTag(vector<double> &mono, map<int, filter_tag_msg_ptr> &tagMap) {
    map<int, filter_tag_msg_ptr> resultTagMap;
    vector<int> indexKey;
    map<int, int> indexTagSizeMap;
    for (map<int, filter_tag_msg_ptr>::iterator it = tagMap.begin(); it != tagMap.end(); ++it) {
        indexKey.push_back(it->first);
        indexTagSizeMap.insert(make_pair(it->first, it->second->getTagMap().size()));
    }

    vector<int> saveIndex;
    for (int i = 0; i < indexKey.size(); ++i) {
        if (indexKey[i] == -1) {
            continue;
        }

        // record cur dup mono index
        map<int, int> dupIndexJIndexMap;
        dupIndexJIndexMap.insert(make_pair(indexKey[i], i));
        for (int j = i + 1; j < indexKey.size(); ++j) {
            if (isDuplicate(mono[indexKey[i]], mono[indexKey[j]])) {
                dupIndexJIndexMap.insert(make_pair(indexKey[j], j));
            } else {
                break;
            }
        }

        // select best mono index by tag size
        int bestIndex = 0;
        int bestTagSize = 0;
        for (map<int, int>::iterator it = dupIndexJIndexMap.begin(); it != dupIndexJIndexMap.end(); ++it) {
            if (indexTagSizeMap[it->first] > bestTagSize) {
                bestIndex = it->first;
                bestTagSize = indexTagSizeMap[it->first];
            }
        }

        // put -1 in to indexKey
        for (map<int, int>::iterator it = dupIndexJIndexMap.begin(); it != dupIndexJIndexMap.end(); ++it) {
            if (it->first != bestIndex) {
                //                int j = dupIndexJIndexMap[it->first];
                indexKey[dupIndexJIndexMap[it->first]] = -1;
            }
        }

        saveIndex.push_back(bestIndex);
    }

    return saveIndex;
}

void acidPreProcess(map<char, double> &acidMap, map<char, char> &replaceAcidMap, map<double, char> &standAcidMap) {
    double H = 1.01;
    // standardization acid mass , precision is 0.00 . and record same acid
    for (map<char, double>::iterator it = acidMap.begin(); it != acidMap.end(); ++it) {
        double acidMass = round(it->second * 100) / 100;
        if (replaceAcidMap.count(it->first)) {
            standAcidMap.insert(make_pair(acidMass, replaceAcidMap[it->first]));
        } else {
            standAcidMap.insert(make_pair(acidMass, it->first));
        }
    }

    // replace the scope acid mass
    map<double, char> extendAcidMap;
    for (map<double, char>::iterator it = standAcidMap.begin(); it != standAcidMap.end(); ++it) {
        extendAcidMap.insert(make_pair(it->first - H, it->second));
        extendAcidMap.insert(make_pair(it->first, it->second));
        extendAcidMap.insert(make_pair(it->first + H, it->second));
    }

    for (map<double, char>::iterator it = extendAcidMap.begin(); it != extendAcidMap.end(); ++it) {
        cout << it->first << " " << it->second << endl;
    }
    standAcidMap = extendAcidMap;
}

bool compareTagLength(string a, string b) {
    return a.length() > b.length();
}

void msalign_tag::getTagByMonoMassV2(const msalign_ptr &msalignPtr,
                                     map<char, double> &acidMap) {
    vector<double> mono = msalignPtr->getIonsMassContainer();
    vector<double> pMass;
    for (map<char, double>::iterator it = acidMap.begin(); it != acidMap.end(); ++it) {
        pMass.push_back(it->second);
    }
    sort(mono.begin(), mono.end());
    sort(pMass.begin(), pMass.end());

    map<double, char> standAcidMap;
    map<int, char> massChar; // 前一个为标准化之后的数据，后一个为对应的氨基酸
    // standardization acid mass , precision is 0.00 . and record same acid
    msalign_tag::replaceCMap(acidMap, acidReplaceMap, massChar); // 得到替换的氨基酸map

    // acidPreProcess(cMass,replaceMap,standAcidMap);  (acid pre process . deprecate !! )

    // get acid tag (key -> indexI , value -> tagMap(key -> indexj , value -> tag) )
    // map<int, filter_tag_msg_ptr> coarseTagMap = getTagByMonoAndIntCharV2(mono, massChar, tagsSet);

    // // remove duplicate tag
    // vector<int> saveIndex = removeDupMonoTag(mono, coarseTagMap);

    // // get final tag map
    // for (int i = 0; i < saveIndex.size(); ++i)
    // {
    //     tagMap.insert(make_pair(saveIndex[i], coarseTagMap[saveIndex[i]]));
    // }

    vector<double> acidMassVec;
    for (map<char, double>::iterator acidIt = acidMap.begin(); acidIt != acidMap.end(); ++acidIt) {
        acidMassVec.push_back(acidIt->second + H);
        acidMassVec.push_back(acidIt->second);
        acidMassVec.push_back(acidIt->second - H);
    }
    sort(acidMassVec.begin(), acidMassVec.end());

    //get tags
    set<string> tagsSet;
    getTagByMonoAndIntCharV2(mono, massChar, tagsSet, acidMap, acidMassVec);

    // get all tag in tagsSet
    for (auto it = tagMap.begin(); it != tagMap.end(); ++it) {
        for (auto mapIt = it->second->getTagMap().begin(); mapIt != it->second->getTagMap().end(); ++mapIt) {
            map<string, double> errorMap = it->second->getMassErrorMap();
            if (mapIt->second.size() > minTagLen && !tagsSet.count(mapIt->second)) {
                tagsSet.insert(mapIt->second);
                if (tagErrorMap.count(mapIt->second)) {
                    if (tagErrorMap[mapIt->second] > errorMap[to_string(mapIt->first) + mapIt->second]) {
                        tagErrorMap[mapIt->second] = errorMap[to_string(mapIt->first) + mapIt->second];
                    }
                } else {
                    tagErrorMap.insert(make_pair(mapIt->second, errorMap[to_string(mapIt->first) + mapIt->second]));
                }
            }
        }
    }

    for (string t : tagsSet) {
        tagsVec.push_back(t);
    }

    // sort by tag length
    sort(tagsVec.begin(), tagsVec.end(), compareTagLength);

    // // free malloc
    // for (map<int, filter_tag_msg_ptr>::iterator it = coarseTagMap.begin(); it != coarseTagMap.end(); ++it)
    // {
    //     if (tagMap.count(it->first) == 0)
    //     {
    //         delete it->second;
    //     }
    // }
}
