//
// Created by Administrator on 2023/3/30.
//

#include "filter_tag_msg.hpp"
#include <fstream>
#include <iomanip>

void filter_tag_msg::filterTagByProteinSeq(string &proteinSeq, string &revProSeq) {
    for (auto it = tagMap.begin(); it != tagMap.end(); ++it) {
        if (it->second.size() < 2) {
            continue;
        }
        if (proteinSeq.find(it->second) != std::string::npos) {
            seqTagMapN.insert(make_pair(it->first, it->second));
        }
        if (revProSeq.find(it->second) != std::string::npos) {
            seqTagMapC.insert(make_pair(it->first, it->second));
        }
    }
}

void filter_tag_msg::outProteinSeqMap(vector<double> &monoMass) {
    if (seqTagMapN.size() == 0 && seqTagMapC.size() == 0) {
        return;
    }
    cout << "index = " << curIndex << " -> " << monoMass[curIndex] << endl;
    cout << " ============N:============== " << endl;
    for (multimap<int, string>::iterator it = seqTagMapN.begin(); it != seqTagMapN.end(); ++it) {
        cout << it->first << " -> " << monoMass[it->first]
             << " sub -> " << monoMass[curIndex] - monoMass[it->first]
             << " -> " << it->second << endl;
    }
    cout << " ============C:============== " << endl;
    for (multimap<int, string>::iterator it = seqTagMapC.begin(); it != seqTagMapC.end(); ++it) {
        cout << it->first << " -> " << monoMass[it->first]
             << " sub -> " << monoMass[curIndex] - monoMass[it->first]
             << " -> " << it->second << endl;
    }
    cout << endl;
}

void filter_tag_msg::outTheMapMsg(vector<double> &monoMass) {
    ofstream out("filter-tag-msg.txt", ios::out);
    out << "index = " << curIndex << " -> " << monoMass[curIndex] << endl;
    for (multimap<int, string>::iterator it = tagMap.begin(); it != tagMap.end(); ++it) {
        out << fixed << setprecision(6)
            << it->first
            << " -> " << monoMass[it->first]
            << " sub -> " << monoMass[curIndex] - monoMass[it->first]
            << " tag -> " << it->second
            << " massError -> " << massErrorMap[to_string(it->first) + it->second]
            << " avgMassError -> " << massErrorMap[to_string(it->first) + it->second] / (double)it->second.length()
            << endl;
    }
    out << endl;
    out.close();
}
