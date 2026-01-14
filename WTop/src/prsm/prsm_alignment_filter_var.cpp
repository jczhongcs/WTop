//
// Created by Administrator on 2023/3/24.
//

#include "prsm_alignment_filter_var.hpp"


void prsm_alignment_filter_var::index_ions_mass(vector<double> &sub, vector<double> &index,
                                                vector<double> &modifier, double &modifier_mass)
{
    int i = 0;
    for (i = 0; i < sub.size(); ++i) {
        //判断该值小于修饰值总质量，同时为一定修饰组合的质量
        if ((sub[i] <= modifier_mass + 5) && (mylib::data_stream::find_key(sub[i], modifier) != -1)) {
            index.push_back(sub[i]);
        } else if (fabs(sub[i]) <= 2) {      //是否该值为1以下？
            index.push_back(sub[i]);
        } else {                         //既不是质量偏移峰，也不是正确的峰
            index.push_back(INT_MAX);
        }
    }
    //------------------过滤一下----------------//
    //有些离子的组合，无法到修饰质量 . 例如 14*6 = 74 + ？ = 98
    for (i = 0; i < index.size(); ++i) {
        int temp = 0;
        if (index[i] != -1 && index[i] < modifier_mass - 2) {
            temp = mylib::data_stream::remove_false(index[i], modifier, modifier_mass);
            if (temp == 0) {
                index[i] = INT_MAX;
            }
        }
    }
}