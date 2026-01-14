#ifndef INITSPAC_HPP
#define INITSPAC_HPP

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include "modify.hpp"
#include "protein_processor.hpp"

#include "../mylib/speculate_lib.h"
#include "../mylib/data_steam.hpp"

#include "../util/string_utils.h"

using namespace std;

namespace prsm{

    class Msalign{
    public:
    modify_ptr unknownModPtr ;
    std::vector<protein_processor_ptr> candidate_protein_container; //过滤蛋白容器

    std::string Id ;
    std::string scans ;
    std::string retention_time ;
    std::string activation;
    std::string activation_sub ;
    std::string ms_one_id;
    std::string ms_one_scans;
    std::string precursor_mz;
    std::string precursor_charge;
    std::string precursor_mass;
    double precursor_mass_double = 0 ;
    std::string precursor_intensity;

    std::vector<double> ions_mass_container;

    std::map<string,string> ms_msg_map;

    map<double,double> keyMsValueAbRadio;  //key - mass , value - Radio
    map<int,int> compIonsIndexMap; //存放互补离子下标
    map<string,string> unknownPTMMassMap ; //存放未知修饰的map，key PTM or Unknown , value mass .

    protein_processor_ptr matchProtein ; //最佳匹配的蛋白
    vector<node> ionsMessageContainerN ; //离子信息容器
    vector<node> ionsMessageContainerC ; //离子信息容器
    vector<double> theoryMassN ;
    vector<double> theoryMassC ;

    double best_score = 0;  //最佳分值
    double best_peak = 0;
    double theory_proteoform_mass ;
    
    double modify_mass;
    
    string match_protein_name;
    string match_protein_seq;
    vector<modi>  best_ptm;

    int cut_c_n = 0;    //N cut
    int cut_n_c = 0;    //C cut

    double mut_s = 0 ;
    long fragment_ions = 0 ;

    long bestIonsPairNub = 0 ;

    double huBuCount = 0 ;
    double seqScore = 0 ;

    int seqScoreCount = 0 ;
    int matchFragmentIonsSize = 0 ; //匹配碎片离子数

    void Init_Msalign_Data(ifstream &ms_file);

    void findComplementaryIons();

    void initCompIonsIndexMap() ;

    int getComplementIonsSize() {
        return compIonsIndexMap.size();
    }

    };
    typedef std::shared_ptr<Msalign> msalign_processor_share_ptr;
    typedef Msalign* msalign_processor_ptr ;
}
#endif