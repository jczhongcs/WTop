//
// Created by wenzhong on 2023/3/21.
//

#include "ms_protein_filter.hpp"

vector<int> findSubsequences(const string &A, const string &B) {
    vector<int> result;

    int n = A.size();
    int m = B.size();

    for (int i = 0; i <= n - m; ++i) {
        if (A.substr(i, m) == B) {
            result.push_back(i);
        }
    }

    return result;
}

void getTagString(const protein_ptr prIt, const msalign_ptr msIt, string &proteinTSeq,
                  map<int, vector<string>> &indexTagMMapN, map<int, vector<string>> &indexTagMMapC) {
    string proteinSeq = proteinTSeq;
    multimap<int, string> indexTagNMMap, indexTagCMMap;
    for (auto it = indexTagMMapN.begin(); it != indexTagMMapN.end(); ++it) {
        for (string s : it->second) {
            indexTagNMMap.insert(make_pair(it->first, s));
        }
    }

    for (auto it = indexTagMMapC.begin(); it != indexTagMMapC.end(); ++it) {
        for (string s : it->second) {
            indexTagCMMap.insert(make_pair(it->first, s));
        }
    }

    protein_ptr proteinPtr = new protein();
    for (multimap<int, string>::iterator it = indexTagNMMap.begin(); it != indexTagNMMap.end(); ++it) {
        set<string> judgeSet;
        for (multimap<int, string>::iterator itRef = indexTagNMMap.begin(); itRef != indexTagNMMap.end(); ++itRef) {
            if (it->first >= itRef->first || judgeSet.count(itRef->second)) {
                continue;
            }
            judgeSet.insert(itRef->second);

            vector<int> itTagLocVec = findSubsequences(proteinSeq, it->second);
            vector<int> itRefTagLocVec = findSubsequences(proteinSeq, itRef->second);

            for (int itTagLoc : itTagLocVec) {
                for (int itRefTagLoc : itRefTagLocVec) {
                    if (itTagLoc >= itRefTagLoc - it->second.size()) {
                        continue;
                    }
                    string subStr = prIt->getProteinSequence().substr(itTagLoc + it->second.size(), itRefTagLoc + itRef->second.size() - (itTagLoc + it->second.size()));
                    double subStrMass = proteinPtr->getAcidSeqMassN(subStr);

                    if (subStr.size() >= 2 && abs(msIt->getIonsMassContainer()[itRef->first] - msIt->getIonsMassContainer()[it->first] - subStrMass) < 2) {
                        cout << msIt->getIonsMassContainer()[itRef->first]
                             << "(" << itRef->first << ")"
                             << " - " << msIt->getIonsMassContainer()[it->first]
                             << "(" << it->first << ")"
                             << " -> " << msIt->getIonsMassContainer()[itRef->first] - msIt->getIonsMassContainer()[it->first]
                             << " -> " << itRef->second << " -> " << it->second
                             << " -> " << subStr
                             << " -> " << subStrMass
                             << " -> " << msIt->getIonsMassContainer()[itRef->first] - msIt->getIonsMassContainer()[it->first] - subStrMass << endl;
                    }
                }
            }
        }
    }

    string revProteinSeq = string(proteinSeq.rbegin(), proteinSeq.rend());
    string revRawProSeq = string(prIt->getProteinSequence().rbegin(), prIt->getProteinSequence().rend());
    for (multimap<int, string>::iterator it = indexTagCMMap.begin(); it != indexTagCMMap.end(); ++it) {
        set<string> judgeSet;
        for (multimap<int, string>::iterator itRef = indexTagCMMap.begin(); itRef != indexTagCMMap.end(); ++itRef) {
            if (it->first >= itRef->first || judgeSet.count(itRef->second)) {
                continue;
            }
            judgeSet.insert(itRef->second);

            vector<int> itTagLocVec = findSubsequences(revProteinSeq, it->second);
            vector<int> itRefTagLocVec = findSubsequences(revProteinSeq, itRef->second);

            for (int itTagLoc : itTagLocVec) {
                for (int itRefTagLoc : itRefTagLocVec) {
                    if (itTagLoc >= itRefTagLoc - it->second.size()) {
                        continue;
                    }
                    string subStr = revRawProSeq.substr(itTagLoc + it->second.size(), itRefTagLoc + itRef->second.size() - (itTagLoc + it->second.size()));
                    double subStrMass = proteinPtr->getAcidSeqMassN(subStr);

                    if (subStr.size() >= 2 && abs(msIt->getIonsMassContainer()[itRef->first] - msIt->getIonsMassContainer()[it->first] - subStrMass) < 2) {
                        cout << msIt->getIonsMassContainer()[itRef->first]
                             << "(" << itRef->first << ")"
                             << " - " << msIt->getIonsMassContainer()[it->first]
                             << "(" << it->first << ")"
                             << " -> " << msIt->getIonsMassContainer()[itRef->first] - msIt->getIonsMassContainer()[it->first]
                             << " -> " << itRef->second << " -> " << it->second
                             << " -> " << subStr
                             << " -> " << subStrMass
                             << " -> " << msIt->getIonsMassContainer()[itRef->first] - msIt->getIonsMassContainer()[it->first] - subStrMass << endl;
                    }
                }
            }
        }
    }

    delete proteinPtr;
}

int ms_protein_filter::getFilterTagMatch(const protein_ptr prIt, const msalign_ptr msIt,
                                         const map<char, char> &ccMap, map<int, filter_tag_msg_ptr> &tagMap,
                                         vector<filter_tag_match_ptr> &filterTagMatchNVec,
                                         vector<filter_tag_match_ptr> &filterTagMatchCVec) {
    // get protein sequence
    string proSeq = prIt->getProteinSequence();
    // replace acid by ccMap
    //    cout <<"\n"<< proSeq <<endl;
    for (auto it = ccMap.begin(); it != ccMap.end(); ++it) {
        //        cout << it->first << " translate -> " << it->second << endl ;
        replace(proSeq.begin(), proSeq.end(), it->first, it->second);
    }
    string revProSeq = string(proSeq.rbegin(), proSeq.rend());
    //    cout << proSeq << endl ;

    // get all tag
    map<int, vector<string>> indexTagMMapN;
    map<int, vector<string>> indexTagMMapC;

    for (map<int, filter_tag_msg_ptr>::iterator it = tagMap.begin(); it != tagMap.end(); ++it) {
        filter_tag_msg_ptr filterTagMsgPtr = it->second;
        for (map<int, string>::const_iterator fit = filterTagMsgPtr->getTagMap().begin();
             fit != filterTagMsgPtr->getTagMap().end(); ++fit) {
            if (fit->second.size() < minTagLength) {
                continue;
            }
            if (proSeq.find(fit->second) != std::string::npos) {
                if (indexTagMMapN.count(it->first)) {
                    indexTagMMapN[it->first].push_back(fit->second);
                } else {
                    vector<string> tagVec;
                    tagVec.push_back(fit->second);
                    indexTagMMapN.insert(make_pair(it->first, tagVec));
                }
            } else if (revProSeq.find(fit->second) != std::string::npos) {
                if (indexTagMMapC.count(it->first)) {
                    indexTagMMapC[it->first].push_back(fit->second);
                } else {
                    vector<string> tagVec;
                    tagVec.push_back(fit->second);
                    indexTagMMapC.insert(make_pair(it->first, tagVec));
                }
            }
        }
    }

    int maxTagLengthN = 0, maxTagLengthC = 0;
    string maxTagN, maxTagC;
    multimap<int, string> indexTagNMMap, indexTagCMMap;
    for (auto it = indexTagMMapN.begin(); it != indexTagMMapN.end(); ++it) {
        for (string s : it->second) {
            if (s.size() > maxTagLengthN) {
                maxTagN = s;
                maxTagLengthN = s.size();
            }
            indexTagNMMap.insert(make_pair(it->first, s));
        }
    }

    for (auto it = indexTagMMapC.begin(); it != indexTagMMapC.end(); ++it) {
        for (string s : it->second) {
            if (s.size() > maxTagLengthC) {
                maxTagC = s;
                maxTagLengthC = s.size();
            }
            indexTagCMMap.insert(make_pair(it->first, s));
        }
    }

    //    getTagString(prIt,msIt,proSeq,indexTagMMapN,indexTagMMapC);
    return maxTagLengthN + maxTagLengthC;

    // filter the same index sub tag
    for (auto it = indexTagMMapN.begin(); it != indexTagMMapN.end(); ++it) {
        vector<string> tagVec = it->second;
        vector<string> removeSubTag;
        for (int i = 0; i < tagVec.size(); ++i) {
            bool isSubTag = false;
            for (int j = 0; j < tagVec.size(); ++j) {
                if (i != j && tagVec[j].find(tagVec[i]) != std::string::npos) {
                    isSubTag = true;
                    break;
                }
            }
            if (!isSubTag) {
                removeSubTag.push_back(tagVec[i]);
            }
        }
        indexTagMMapN[it->first] = removeSubTag;
    }

    for (auto it = indexTagMMapC.begin(); it != indexTagMMapC.end(); ++it) {
        vector<string> tagVec = it->second;
        vector<string> removeSubTag;
        for (int i = 0; i < tagVec.size(); ++i) {
            bool isSubTag = false;
            for (int j = 0; j < tagVec.size(); ++j) {
                if (i != j && tagVec[j].find(tagVec[i]) != std::string::npos) {
                    isSubTag = true;
                    break;
                }
            }
            if (!isSubTag) {
                removeSubTag.push_back(tagVec[i]);
            }
        }
        indexTagMMapC[it->first] = removeSubTag;
    }

    // filter the diff index sub tag
    multimap<int, vector<string>> fIndexTagMMapN, fIndexTagMMapC;
    for (auto it = indexTagMMapN.begin(); it != indexTagMMapN.end(); ++it) {
        vector<string> tempTagVec;
        for (string curIndexTag : it->second) {
            bool isDiffIndexTag = false;
            for (auto itNext = indexTagMMapN.begin(); itNext != indexTagMMapN.end(); ++itNext) {
                if (itNext == it) {
                    continue;
                }
                if (isDiffIndexTag == true) {
                    break;
                }
                for (string diffIndexTag : itNext->second) {
                    if (diffIndexTag.size() <= curIndexTag.size()) {
                        continue;
                    }
                    if (diffIndexTag.find(curIndexTag) != std::string::npos) {
                        isDiffIndexTag = true;
                        break;
                    }
                }
            }
            if (isDiffIndexTag == true) {
                continue;
            } else {
                tempTagVec.push_back(curIndexTag);
            }
        }
        if (tempTagVec.size() > 0) {
            fIndexTagMMapN.insert(make_pair(it->first, tempTagVec));
        }
    }

    for (auto it = indexTagMMapC.begin(); it != indexTagMMapC.end(); ++it) {
        vector<string> tempTagVec;
        for (string curIndexTag : it->second) {
            bool isDiffIndexTag = false;
            for (auto itNext = indexTagMMapC.begin(); itNext != indexTagMMapC.end(); ++itNext) {
                if (itNext == it) {
                    continue;
                }
                if (isDiffIndexTag == true) {
                    break;
                }
                for (string diffIndexTag : itNext->second) {
                    if (diffIndexTag.size() <= curIndexTag.size()) {
                        continue;
                    }
                    if (diffIndexTag.find(curIndexTag) != std::string::npos) {
                        isDiffIndexTag = true;
                        break;
                    }
                }
            }
            if (isDiffIndexTag == true) {
                continue;
            } else {
                tempTagVec.push_back(curIndexTag);
            }
        }
        if (tempTagVec.size() > 0) {
            fIndexTagMMapC.insert(make_pair(it->first, tempTagVec));
        }
    }

    // start analysis the tag mass gap
    for (auto it = fIndexTagMMapN.begin(); it != fIndexTagMMapN.end(); ++it) {
        filter_tag_match_ptr filterTagMatch = new filter_tag_match();
        filterTagMatch->insertIndexMsMap(it->first, msIt->getIonsMassContainer()[it->first]);
        for (string tag : it->second) {
            vector<int> tagMatchIndexVec = findSubsequences(proSeq, tag);
            filterTagMatch->insertIndexTagBeginMap(to_string(it->first) + tag, tagMatchIndexVec);
            filterTagMatch->insertIndexTagMMap(it->first, tag);
        }
        filterTagMatchNVec.push_back(filterTagMatch);
    }
    for (auto it = fIndexTagMMapC.begin(); it != fIndexTagMMapC.end(); ++it) {
        filter_tag_match_ptr filterTagMatch = new filter_tag_match();
        filterTagMatch->insertIndexMsMap(it->first, msIt->getIonsMassContainer()[it->first]);
        for (string tag : it->second) {
            vector<int> tagMatchIndexVec = findSubsequences(string(proSeq.rbegin(), proSeq.rend()), tag);
            filterTagMatch->insertIndexTagBeginMap(to_string(it->first) + tag, tagMatchIndexVec);
            filterTagMatch->insertIndexTagMMap(it->first, tag);
        }
        filterTagMatchCVec.push_back(filterTagMatch);
    }
}

void analysisTagMatch(vector<filter_tag_match_ptr> &filterTagMatchNVec,
                      vector<filter_tag_match_ptr> &filterTagMatchCVec,
                      const vector<double> &msMass, const string &rawProteinSeq) {
    cout << "----------- N Terminal Tag is : ------------- " << endl;
    for (filter_tag_match_ptr filterTagMatch : filterTagMatchNVec) {
        filterTagMatch->outMsgTest();
    }
    //    cout << "----------- C Terminal Tag is : ------------- " << endl ;
    //    for (filter_tag_match_ptr filterTagMatch : filterTagMatchCVec) {
    //        filterTagMatch->outMsgTest();
    //    }

    multimap<int, int> indexTagLocMap;
    map<string, string> indexTagMap;
    for (filter_tag_match_ptr filterTagMatchPtr : filterTagMatchNVec) {
        for (auto it = filterTagMatchPtr->getIndexTagMMap().begin(); it != filterTagMatchPtr->getIndexTagMMap().end(); ++it) {
            vector<int> tagMatchSeq = filterTagMatchPtr->getIndexTagBeginMap()[to_string(it->first) + it->second];
            for (int tagMatchIndex : tagMatchSeq) {
                indexTagLocMap.insert(make_pair(it->first, tagMatchIndex));
                indexTagMap.insert(make_pair((to_string(it->first) + to_string(tagMatchIndex)), it->second));
            }
        }
    }

    cout << "======================== " << endl;
    protein_ptr proteinPtr = new protein();
    for (multimap<int, int>::iterator it = indexTagLocMap.begin(); it != indexTagLocMap.end(); ++it) {
        for (multimap<int, int>::iterator itRef = indexTagLocMap.begin(); itRef != indexTagLocMap.end(); ++itRef) {
            string tag = indexTagMap[to_string(it->first) + to_string(it->second)];
            if (itRef->first <= it->first || it->second >= itRef->second) {
                continue;
            }
            double massGap = msMass[itRef->first] - msMass[it->first];
            string behindSeq = rawProteinSeq.substr(it->second + tag.size(), itRef->second - (it->second + tag.size()));
            double acidSeqMass = proteinPtr->getAcidSeqMassN(behindSeq);
            if (fabs(acidSeqMass - massGap) > 3) {
                continue;
            } else {
                cout << behindSeq << endl;
                cout << it->first << " -> " << tag << "(" << it->second << ")"
                     << " -> " << itRef->first
                     << " -> " << indexTagMap[to_string(itRef->first) + to_string(itRef->second)] << "(" << itRef->second << ")"
                     << " -> " << massGap << " -> " << acidSeqMass
                     << " -> " << acidSeqMass - massGap << endl;
            }
        }
    }

    delete proteinPtr;
}

int compIonsTagFind(const protein_ptr prIt,
                    const msalign_ptr msIt,
                    const map<char, char> &ccMap,
                    map<int, filter_tag_msg_ptr> &tagMap,
                    map<int, int> &compIonsMap) {
    string proSeq = prIt->getProteinSequence(); // get protein sequence
    string revProSeq;                           // reverse protein sequence
    int tagMinSize = 2;
    // replace acid by ccMap
    for (auto it = ccMap.begin(); it != ccMap.end(); ++it) {
        replace(proSeq.begin(), proSeq.end(), it->first, it->second);
    }
    revProSeq = string(proSeq.rbegin(), proSeq.rend());
    cout << "\n"
         << proSeq << "\n"
         << revProSeq << endl;

    multimap<int, string> seqTagMapN, seqTagMapC; // match N tag , C tag

    for (map<int, filter_tag_msg_ptr>::iterator it = tagMap.begin(); it != tagMap.end(); ++it) {
        for (multimap<int, string>::const_iterator tagIt = it->second->getTagMap().begin();
             tagIt != it->second->getTagMap().end(); ++tagIt) {
            if (tagIt->second.size() < 2) {
                continue;
            }
            if (proSeq.find(tagIt->second) != std::string::npos) {
                seqTagMapN.insert(make_pair(tagIt->first, tagIt->second));
            }
            if (revProSeq.find(tagIt->second) != std::string::npos) {
                seqTagMapC.insert(make_pair(tagIt->first, tagIt->second));
            }
        }
    }

    for (map<int, filter_tag_msg_ptr>::iterator it = tagMap.begin(); it != tagMap.end(); ++it) {
        it->second->outProteinSeqMap(msIt->getIonsMassContainer());
    }

    return 0; //
}

string ms_protein_filter::tagsMatchSequence(const protein_ptr &prIt, const msalign_ptr &msIt,
                                            const complement_ions_ptr &compIonsPtr,
                                            const msalign_tag_ptr &msalignTagPtr,
                                            filter_protein_ptr &filterProPtr) {
    string proSeq = prIt->getProteinSequence();
    map<char, char> acidReplaceMap = msalignTagPtr->getAcidReplaceMap();
    for (map<char, char>::iterator it = acidReplaceMap.begin(); it != acidReplaceMap.end(); ++it) {
        replace(proSeq.begin(), proSeq.end(), it->first, it->second);
    }
    string revProSeq = string(proSeq.rbegin(), proSeq.rend());

    for (string tag : msalignTagPtr->tagsVec) {
        if (filterProPtr->tagN.length() != 0 && filterProPtr->tagC.length() != 0) {
            break;
        }
        if (proSeq.find(tag) != string::npos && filterProPtr->tagN.length() == 0) {
            filterProPtr->tagN = tag;
            filterProPtr->tagMassError += msalignTagPtr->tagErrorMap[tag] / filterProPtr->tagN.length();
        } else if (revProSeq.find(tag) != string::npos && filterProPtr->tagC.length() == 0) {
            filterProPtr->tagC = tag;
            filterProPtr->tagMassError += msalignTagPtr->tagErrorMap[tag] / filterProPtr->tagC.length();
        }
    }
    filterProPtr->filter_tag_length = filterProPtr->tagN.length() + filterProPtr->tagC.length();
    return "0";
    //    cout << proSeq << endl ;

    // vector<filter_tag_match_ptr> filterTagMatchNVec, filterTagMatchCVec;

    // if ms has complement ions , it's to find tag seq ;
    // compIonsTagFind(prIt, msIt, ccMap, tagMap,compIonsMap) ;
    // return getFilterTagMatch(prIt, msIt, ccMap, tagMap, filterTagMatchNVec, filterTagMatchCVec);

    // analysis tag match
    //    analysisTagMatch(filterTagMatchNVec,filterTagMatchCVec,
    //                     msIt->getIonsMassContainer(),prIt->getProteinSequence());
}

/**
 * 根据ccMap 变化蛋白，并且寻找 匹配的tag
 * @param prIt
 * @param msIt
 * @param ccMap     key-变化前的氨基酸 ，value- 变化后的氨基酸
 * @param set   所有的tag
 * @param match_tag_length     记录多少个tag被匹配
 * @param filterMass    由于是补充过滤并且是按照tag ，过滤大质量蛋白
 * @return  返回的是匹配tag的最大长度
 */
int ms_protein_filter::findBestProteinByTag(const protein_ptr prIt, const msalign_ptr msIt,
                                            const map<char, char> &ccMap, const set<string> &set,
                                            const double filterMass) {
    int match_tag_length = 0;
    if (fabs(filterMass) > 2000) {
        return 0;
    }
    // 得到该条蛋白序列
    string proSeq = prIt->getProteinSequence();
    // 替换
    for (auto it = ccMap.begin(); it != ccMap.end(); ++it) {
        //        cout << it->first << " " << it->second << endl ;
        replace(proSeq.begin(), proSeq.end(), it->first, it->second);
    }
    //    cout << proSeq << endl ;
    //    cout << " ==== " << endl ;
    //    for (string s : set) {
    //        cout << s << endl ;
    //    }
    int count = 0;
    int maxTagSize = 0;
    for (string str : set) {
        if (str.size() < 2)
            continue; // 如果只有1个氨基酸，则直接过滤
        string rev = string(str.rbegin(), str.rend());
        if (proSeq.find(str) != -1) {
            count++;
            maxTagSize = maxTagSize > str.size() ? maxTagSize : str.size();
            //            cout << "n terminal match seq : " << str << endl ;
        } else if (proSeq.find(rev) != -1) {
            count++;
            maxTagSize = maxTagSize > str.size() ? maxTagSize : str.size();
            //            cout << "c terminal match seq : " << str << endl ;
        }
    }
    match_tag_length = maxTagSize;
    return maxTagSize;
}

int ms_protein_filter::zero_match_peaks_search(const vector<double> &mono_mass,
                                               const vector<double> &theo_mass,
                                               int &one_ptm_match_peaks) {
    int count = 0;
    vector<int> seq_id;
    vector<double> short_1;
    vector<double> long_1;
    if (mono_mass.size() > theo_mass.size()) {
        long_1 = mono_mass;
        short_1 = theo_mass;
    } else {
        short_1 = mono_mass;
        long_1 = theo_mass;
    }
    ppm_value_ptr ppmValuePtr = std::make_shared<ppm_value>(15.0);
    for (int i = 0; i < short_1.size(); i++) { // 获取匹配个数
        int n = mylib::speculate_lib::txpb_binarey_search_ex(long_1, long_1.size(), short_1[i]);
        if (ppmValuePtr->calculate_ppm_value(short_1[i], long_1[n])) {
            count++;
        }
    }
//    delete ppmValuePtr;

    if (count < 2) {
        /**
         * 进行N端的搜寻
         * 1、判定最大截断位置
         * 2、添加修饰峰
         * 3、返回取匹配最优的
         */
        //             int maxCount = getMaxCutLocETD(prIt,msIt,mdIt,10);
        //             seq_match_peaks = maxCount;
        return 0;
    }
    int max_seq_peaks = 0;
    for (int i = 0; i < seq_id.size(); ++i) {
        int c = 1;
        // 1 3 4 5 9
        for (int j = i + 1; j < seq_id.size(); ++j) {
            if (seq_id[j] - seq_id[j - 1] < 1) {
                c++;
            } else if (c > max_seq_peaks) {
                max_seq_peaks = c;
                i = j - 1;
                break;
            }
        }
    }
    one_ptm_match_peaks = max_seq_peaks;
    return count;
}

bool ms_protein_filter::precursor_mass_scope_search(const vector<double> &mono_mass, double &mass_gap,
                                                    double precursor_mass, double scope) {
    if (fabs(mono_mass.back() - precursor_mass) < scope) {
        mass_gap = fabs(mono_mass.back() - precursor_mass);
        return true;
    }
    return false;
}