#ifndef INIT_PROTE_HPP
#define INIT_PROTE_HPP

#include <memory>
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <vector>
#include <map>

#include "temp_prsm_argument.hpp"

using namespace std;

namespace prsm{

    class Protein{
        public:


        string protein_name;     //蛋白质名称
        string protein_sequence;      //蛋白质序列
        vector<double> theoryMassC_By;      //理论质量
        vector<double> theoryMassN_By;      //理论质量

        vector<double> theoryMassC_Cz;      //理论质量
        vector<double> theoryMassN_Cz;      //理论质量

        vector<double> seqMass ;

        Protein();

        Protein(const string &proteinName, const string &proteinSequence);

        void initTitleSeq(const string &title,const string &seq);
        void initTheoryMass() ;

        void Init_Sequence(ifstream &prote_file_ptr,int &first_flag);

        void Init_Sequence_re(ifstream &prote_file_ptr);

        void get_seq_mass(string seq , vector<double> & seq_mass) ;

        //得到N端截断个数的质量
        double getNCutMass(int cutSize);

        std::map<char,double> getMap();

        std::string next_title;
        std::string first_title;
        std::string sequence;

        int matchTagSize = 0 ;  //匹配tag的数量
        int maxTagLen = 0 ; //匹配的tag最大长度
        int one_ptm_match_peaks = 0 ;
        int one_ptm_seq_match_peaks = 0 ;

        double complementIonsInProtein = 0 ; //占比

        double leastSub = 0 ; //用于过滤时范围差

    };
    typedef std::shared_ptr<Protein> protein_processor_share_ptr;
    typedef Protein* protein_processor_ptr ;
}
#endif