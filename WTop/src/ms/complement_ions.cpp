//
// Created by wenzhong on 2023/3/21.
//

#include "complement_ions.hpp"


// 得到MS中的互补离子index
void complement_ions::getComplementIonsMap(const msalign_ptr & msalignPtr)
{
//    vector<double> msVec;
//    for (int i = 0; i < mono.size(); ++i) {
//        if (i == 0) {
//            msVec.push_back(mono[i]);
//        }
//        else if ( msVec.size() > 0 && fabs(mono[i] - msVec[msVec.size()-1]) > 2) {
//            msVec.push_back(mono[i]);
//        }
//    }
    vector<double> monoVec = msalignPtr->getIonsMassContainer();
    double precursorMass = msalignPtr->getPrecursorMass();
    set<int> s ; //用于去重复
    for (int i = 0 ; i < monoVec.size(); ++i) {
        double t = precursorMass - monoVec[i];
        int index = mylib::data_stream::txpb_binarey_search_ex(monoVec, monoVec.size(), t);
        if ( index >= 0 && !s.count(index) && fabs(t - monoVec[index]) < massError) {
            indexMap.insert(make_pair(i,index));
            indexSet.insert(i);
            indexSet.insert(index);
            s.insert(i);
        }
    }


}