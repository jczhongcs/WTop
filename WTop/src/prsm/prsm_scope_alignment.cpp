//
// Created by wenzhong on 2023/3/22.
//

#include "prsm_scope_alignment.hpp"


void prsm_scope_alignment::scopeAlignment(const vector<double> &reference, const vector<double> &peer,
                                          std::vector<std::pair<long, long> > &alignment,
                                          std::vector<double> &modifier,
                                          const double &modifier_mass, double maxMass) {
    vector<pair<long, long>> addAlignment;
    int sex = 5 ; //5da 常熟项
    set<long> removeRef;
    for (int i = 0 ; i < alignment.size(); ++i) {
        if (removeRef.count(alignment[i].first))    {  // 如果一个理论值已经对准过，直接进行下一个
            continue;
        } else {
            removeRef.insert(alignment[i].first);
        }
        set<long> removePeer;
        removePeer.insert(alignment[i].second);
        double theoMass = reference[alignment[i].first];
        //向上边查找
        for(long left = alignment[i].second - 1; left >= 0 ; --left) {
            double monoMass = peer[left] ;
            if (fabs(monoMass - theoMass) > maxMass) {
                break;
            }
            if ( removePeer.count(left) )  {
                continue ;
            }
            double monoSubTheo = monoMass - theoMass ;
            if (fabs(monoSubTheo) < sex) {   //1、当两者相差不大
                addAlignment.push_back(make_pair(alignment[i].first,left));
            } else if (fabs(monoSubTheo - modifier_mass) < sex) {   //2.1 当是 modifier_mass 的值，直接加入
                addAlignment.push_back(make_pair(alignment[i].first,left));
            } else {    //2.2 当是modifier_mass的子值，判断是否有可达modifier_mass的质量，如果有，则加入，没有，则不加入
                double subMass = modifier_mass - monoSubTheo ;
                double searchIndex = mylib::data_stream::txpb_binarey_search_ex(modifier, modifier.size(), subMass);
                if (searchIndex < 0 ) continue;
                double temp =   modifier[searchIndex]  - subMass ;
                if (fabs(temp) < sex) {
                    addAlignment.push_back(make_pair(alignment[i].first,left));
                }
            }
            removePeer.insert(left);
        }
        //向下边查找
        for(long right = alignment[i].second+1; right < peer.size() ; ++right) {  //向下边查找 指向mono的指针
            double monoMass = peer[right] ;
            if (fabs(monoMass - theoMass) > maxMass) {
                break;
            }
            if ( removePeer.count(right) )  {
                continue ;
            }
            double monoSubTheo = monoMass - theoMass ;
            if (fabs(monoSubTheo) < sex) {   //1、当两者相差不大
                addAlignment.push_back(make_pair(alignment[i].first,right));
            } else if (fabs(monoSubTheo - modifier_mass) < sex) {   //2.1 当是 modifier_mass 的值，直接加入
                addAlignment.push_back(make_pair(alignment[i].first,right));
            } else {    //2.2 当是modifier_mass的子值，判断是否有可达modifier_mass的质量，如果有，则加入，没有，则不加入
                double subMass = modifier_mass - monoSubTheo ;
                double searchIndex = mylib::data_stream::txpb_binarey_search_ex(modifier, modifier.size(), subMass);
                if (searchIndex < 0 ) continue;
                double temp =  modifier[searchIndex] - subMass ;
                if (fabs(temp) < sex) {
                    addAlignment.push_back(make_pair(alignment[i].first,right));
                }
            }
            removePeer.insert(right);
        }
    }
    for(int i = 0 ; i < addAlignment.size() ; ++i) {
        alignment.push_back(addAlignment[i]);
    }
    std::sort(alignment.begin(), alignment.end(), mylib::data_stream::cmp1);
    alignment.erase(unique(alignment.begin(), alignment.end()), alignment.end());
}
