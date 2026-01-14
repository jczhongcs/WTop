//
// Created by wenzhong on 2022/10/18.
//
/**
 * 该类用于存放prsm的参数设置
 */
#ifndef DTW_TOPZ_INITPRSMARGUMENT_HPP
#define DTW_TOPZ_INITPRSMARGUMENT_HPP

#include <memory>
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <vector>
#include <set>

using namespace std;

namespace prsm {

    class config_argument {
    public:
        const double h20 = 18.01056 ;     //h20
        const double core = 1.9919 ;
        const double H = 1.007276 ;
        const double scope_h = 1.5 ; //mass 掉H少H的范围设置
        const double min_value = 0.000001 ;  //用于判断double是否为0
        const double min_ppm_value = 15 ;  // ppmValue

        int selectPrSMsArg = 1 ; // 1标识进行未知修饰鉴定，0标识鉴定已知修饰鉴定

        double sub_etd = h20 - core ;
//        double searchNumb = 5000 ;       //搜寻库的数量
        double min_precursor_mass = 300;       //前提质量最小

        int minMonoPeaks = 5 ;       // 最小质谱峰数
        double ptmMassMax = 500 ;       //已知修饰的最大质量
        double ptmMassMin = -500 ;       //已知修饰的最大质量

        set<string> BYIonsActions {"HCD","CID"};    //BY离子的action
        set<string> CZIonsActions {"ETD"};      //CZ离子的action

        int unknownPTMNumber = 1 ; //未知修饰个数，支持0,1
        double unknownPTMsMassMax = 500 ;   // 未知修饰的质量最大质量
        double unknownPTMsMassMin = -500 ;   //未知修饰质量最小质量
        double unknownPTMsFilterProteinMassScope = 3000 ; //未知修饰第三步过滤的范围大小

        int unknownPTMsFilterCountSize = 2 ; //未知修饰过滤第一步和第二步的最小长度条件
        int unknownFilterOneAndTwoTopSelect = 24 ;  //未知修饰过滤第一步和第二步的top选取个数
        int unknownFilterThreeTopSelect = 100 ;     //未知修饰第三步的top选取个数

        double cutSizeMin = 1 ; //H2a,h2b,h3 ,h4 = 10 , blod = 1;
        // 未知修饰的质量在(-0.015 - 0.015)中浮动  对应10000da的ppm15
        double unknownPTMsMassFlog = -0.15 ;
        double unknownPTMsMassFlogChangeMass = 0.01 ;   //每次变化 0.001
        int unknownPTMsMassFlogFloatSize = 30 ; // 变化30次(-0.15~0.15)

    };
    typedef std::shared_ptr<config_argument> config_arg_share_ptr;
    typedef config_argument* config_argument_ptr;
}


#endif //DTW_TOPZ_INITPRSMARGUMENT_HPP
