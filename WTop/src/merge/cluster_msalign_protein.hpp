#ifndef SAVEMSPR_HPP__
#define SAVEMSPR_HPP__

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <algorithm>

#include "msalign_processor.hpp"
#include "temp_prsm_argument.hpp"
#include "config_argument.hpp"

#include "../mylib/terminal_truncation.h"
#include "../mylib/calculate_lib.h"
using namespace std ;
/*
 * 主要用于第一步质谱与蛋白质数据库聚类的函数
 * */
class cutMass {
    public : 
    int Ncut = 0 ;
    int Ccut = 0 ; 
    int SumCut = 0 ;
    int log = 0 ;   //  记录蛋白质指针下标
};
typedef std::shared_ptr<cutMass> cutMassPtr ; 
typedef cutMass* cutMassPtrP; 

namespace prsm {

    //质谱聚类，取整，整数相同的为同一类，搜同样的截断后选取的蛋白质序列库
    class dataSet { 
        public:
        std::vector<protein_processor_ptr> protein ; //保存蛋白质数据集
        std::vector<cutMassPtrP> cutPtms ;     //保存截断信息

        int prNumber = 0 ; 
    };
    typedef std::shared_ptr<dataSet> dataSetPtr; 
    typedef dataSet* dataSetPtrP ; 

    void sameConbine(const vector<msalign_processor_ptr> &msContainer, const vector<protein_processor_ptr> &prContainer,
                     const vector<double> &ptmMassSeq,
                     map<long,dataSetPtrP> &CIDmslogMap, map<long,dataSetPtr> &ETDmslogMap, double addMassCID, double addMassETD); //----------------聚

    void spectProteinFormReWriteCID(protein_processor_ptr pro, const vector<double> &ptmMass ,
                                    const double & precursorMass , dataSetPtrP p, double addMass); //CID 分析

    void spectProteinFormReWriteETD(protein_processor_ptr pro, const vector<double> &ptmMass ,
                                    const double & precursorMass , dataSetPtrP p, double addMass); //ETD 分析

    void spectProteinFormReWriteETD_1(protein_processor_share_ptr pro, const vector<double> &ptmMass ,
                                      const double & precursorMass , dataSetPtr p, double addMass);
    void mslogProduceProteinData(map<long,dataSetPtr> &CIDmslogMap, const vector<protein_processor_share_ptr> &prContainer, const map<double,string>& modifyTable, const vector<double> & PtmsSeq);

    void sortdataSetContiner(map<long,dataSetPtr> &CIDmslogMap);

    //20221114
    void judgeProteinFormCID(protein_processor_ptr pro, const vector<double> &ptmMass, const double &precursorMass, config_argument_ptr &prsmArg,
                             map<int,int> &cutMap);
    void judgeProteinFormETD(protein_processor_ptr pro, const vector<double> &ptmMass, const double &precursorMass, config_argument_ptr &prsmArg,
                             map<int,int> &cutMap);

    void DetermineProteinForm(protein_processor_ptr pro, msalign_processor_ptr sp, config_argument_ptr prsmArg, const vector<double> &ptmMass, const double &precursorMass,
                              map<int,int> &cutMap);

}
#endif