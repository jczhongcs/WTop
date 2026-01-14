//
// Created by wenzhong on 2023/3/22.
//

#include "prsm_alignment_filter.hpp"


int prsm_alignment_filter::get_mass_sub(const vector<double> &reference, const vector<double> &reference_y,
                                const vector<double> &peer,
                                vector<double> &sub, vector<double> &sub_y,
                                const std::vector<std::pair<long, long> > &alignment,
                                const std::vector<std::pair<long, long> > &alignment_y) {
    int i = 0;
    for (i = 0; i < alignment.size(); i++) {
        sub.push_back((peer[alignment[i].second] - reference[alignment[i].first]));
    }
    for (i = 0; i < alignment_y.size(); i++) {
        sub_y.push_back(peer[alignment_y[i].second] - reference_y[alignment_y[i].first]);
    }
    return 0;
}

int find_key_prsm(double key, vector<double> &mod) {
    for (int i = 0; i < mod.size(); i++) {
        if ((key >= mod[i] - 2) && (key <= mod[i] + 2)) {
            return i;
        }
    }
    return -1;
}

int remove_false(double index, vector<double> &modifier, const double &modifier_mass)        //过滤掉不可能的值
{
    int i = 0, j = 0;
    if (index > modifier_mass - 2 && index < modifier_mass + 2) {
        return 1;
    }
    for (j = 0; j < modifier.size(); ++j) {
        if ((index + modifier[j] >= modifier_mass - 2) && (index + modifier[j] <= modifier_mass + 2)) {
            return 1;
        }
        if (index + modifier[j] > modifier_mass) {
            return 0;
        }
    }
    return 0;
}

void prsm_alignment_filter::index_ions(vector<double> &sub, vector<double> &index_y,
                                       std::vector<double> &modifier,
                                       double &modifier_mass)    //定位修饰偏差值
{
    int i = 0;
    for (i = 0; i < sub.size(); ++i) {
        //判断该值小于修饰值总质量，同时为一定修饰组合的质量
        if ((sub[i] <= modifier_mass + 5) && (find_key_prsm(sub[i], modifier) != -1)) {
            index_y.push_back(sub[i]);
        } else if (fabs(sub[i]) <= 2) {      //是否该值为1以下？
            index_y.push_back(sub[i]);
        } else {                         //既不是质量偏移峰，也不是正确的峰
            index_y.push_back(INT_MAX);
        }
    }
    //------------------过滤一下----------------//
    //有些离子的组合，无法到修饰质量 . 例如 14*6 = 74 + ？ = 98
    for (i = 0; i < index_y.size(); ++i) {
        int temp;
        if (index_y[i] != -1 && index_y[i] < modifier_mass - 2) {
            temp = remove_false(index_y[i], modifier, modifier_mass);
            if (temp == 0) {
                index_y[i] = INT_MAX;
            }
        }
    }
}

void prsm_alignment_filter::index_ions_mass(vector<double> &sub, vector<double> &index,
                                          std::vector<double> &modifier, double &modifier_mass)
{
    double sex = 5 ;
    for (int i = 0; i < sub.size(); ++i) {
        if (fabs(sub[i]) < sex) {
            index.push_back(sub[i]);
        } else if (fabs(sub[i] - modifier_mass) < 5) {
            index.push_back(sub[i]);
        } else {
            double searchIndex = mylib::data_stream::txpb_binarey_search_ex(modifier, modifier.size(), sub[i]);
            if (searchIndex < 0 || modifier.size() == 0) continue;
            double temp =  modifier[searchIndex] - sub[i] ;
            if (fabs(temp) < sex) {
                index.push_back(sub[i]);
            } else {
                index.push_back(INT_MAX) ;
            }
        }
    }
}

void prsm_alignment_filter::get_ppm_unknown(vector<double> &reference, vector<double> &peer,
                                              std::vector<std::pair<long, long> > &alignment,
                                              std::vector<double> &index, std::vector<double> &modifier, std::vector<double> &ppm,
                                              double &modifier_mass) {
    int i = 0;
    for (i = 0; i < index.size(); i++) {
        int searchIndex = mylib::data_stream::txpb_binarey_search_ex(modifier, modifier.size(), index[i]);   //寻找一个最接近index的修饰质量
        if (searchIndex < 0){
            continue;
        }
        double massShift = modifier[searchIndex] ;
        double theoryMass = reference[alignment[i].first];
        double monoMass = peer[alignment[i].second] ;
        double minPPM = INT_MAX ;
        if (mylib::etd_process::_PPM_H(theoryMass + massShift, monoMass, minPPM)) {
            ppm.push_back(minPPM);
        } else {
            ppm.push_back(INT_MAX);
        }
    }
}

void prsm_alignment_filter::get_ions_peaks(vector<node> &ions,vector<pair<long, long> > &alignment,
                                           vector<double> &index,vector<double> &ppm, double minPPMValue)                //ions取出所有的离子
{
    for (int i = 0; i < alignment.size(); i++)      //分别取出y离子质量偏移为0 ，和带修饰的离子。
    {
        if (ppm[i] != INT_MAX) {
            if (ppm[i] <= minPPMValue && index[i] != INT_MAX) {
                node t;
                t.thoe_id = alignment[i].first;
                t.mono_id = alignment[i].second;
                t.index = index[i];
                t.ppm = ppm[i];
                ions.push_back(t);
            }
        }
    }
}
