#ifndef PRSM_ARGUMENT_HPP
#define PRSM_ARGUMENT_HPP

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "protein_processor.hpp"
#include "modify.hpp"

#include "../mylib/io.h"
#include "../mylib/data_steam.hpp"
#include "../mylib/speculate_lib.h"
#include "../mylib/terminal_truncation.h"

#include "../console/opts.h"

using namespace std;

namespace prsm{
    
    class temp_prsm_argument{
        public:

        vector<c::arg::cutLocationPtr> cutPtrVec;
        int cut_location_n = 0;           //记录末端截断位置
        int cut_location_c = 0;
        double var_mod_mass = 0;         //修饰质量
        double pr_spec_pr_sub = 0;          //理论protemass与实际之差
        double prsm_score;
        int seq_score_count;

        double complement_ions_number = 0 ;

        int matchFragmentIonsSize = 0  ; //匹配碎片离子个数
        std::vector<double> modifier_table;     //修饰表值记录

        vector<double> theory_ions_mass_c;  //y离子
        vector<double> theory_ions_mass_n;  //b离子

        std::vector<std::pair<long,long> > alignment_ions_n; //得到b lons对准数据
        std::vector<std::pair<long,long> > alignment_ions_c; //得到y lons对准数据
        // std::vector<double> pre_reference;
        // std::vector<double> pre_peer;
        vector<double> ref_score_n;       //理论b离子得分
        vector<double> peer_score_n;       //实际b离子得分
        vector<double> ref_score_c;       //理论y离子得分
        vector<double> peer_score_c;       //实际y离子得分

        std::vector<double> sub_n;        //保存b离子相减数据
        std::vector<double> sub_c;       //保存y离子相减数据
        std::vector<double> filter_ions_n;    //保存过滤后的b ions数据
        std::vector<double> filter_ions_c;    //保存过滤后的y ions数据

        std::vector<double> ppm_value_n;
        std::vector<double> ppm_value_c;

        std::vector<double> msMass ; 

        vector<modi> ptms;
        std::vector<locat_ions> location;
        vector<node> ions_n;        //n端离子（b ，c 离子）
        vector<node> ions_c;        //c端离子（y ，z 离子）

        vector<node> merge_ions_c;
        vector<node> merge_ions_n;

        int peaks_match = 0;

        int length = 0;

        // vector<locat_ions> location; 

        double score_c = 0, score_n = 0, score = 0;


        double mut_x = 0 ;

        double true_mass = 0 ;

        long pairIonsNub = 0 ;

        ///调试打分输出特征参数
        double feature1=0;
        double feature2=0;
        double feature3=0;
        double feature4=0;
        double feature5=0;
        double feature6=0;
        double feature7=0;
        double feature8=0;
        double feature9=0;
        double feature10=0;
        double feature11=0;
        double feature12=0;
        double feature13=0;
        double feature14=0;
        double feature15=0;
        double feature16=0;
        double feature17=0;
        double feature18=0;
        double feature19=0;
        double feature20=0;
        double feature21=0;
        double feature22=0;
        double feature23=0;
        double feature24=0;
        double feature25=0;
        double feature26=0;
        double feature27=0;
        double feature28=0;
        double feature29=0;
        double feature30=0;
        double feature31=0;
        double feature32=0;
        double feature33=0;
        double feature34=0;
        double feature35=0;
        double feature36=0;
        double feature37=0;
        double feature38=0;
        double feature39=0;
        double feature40=0;
        double feature41=0;
        double feature42=0;
        double feature43=0;
        double feature44=0;
        double feature45=0;
        double feature46=0;
        double feature47=0;
        double feature48=0;
        double feature49=0;
        double feature50=0;



        ///调试参数结束

        void get_true_modMass ( double mass) {
            true_mass = mass ;
        }
        //----------------------------------------------------------------21.11.1 wzreplace 


        int get_lenth(int seq_length)
        {
            length = seq_length;
            return length;
        }

        int MaxSize = 0 ;
        void get_match_peaks()
        {
            peaks_match = ions_n.size() + ions_c.size();
        }

        void get_match_peaks_ETD()
        {
            peaks_match = ions_n.size() + ions_c.size();
        }

        void get_score(double mux,double h ,int subCor,double seqScore,double ionsA)
        {
            score = (score_c + score_n) * sqrt(h) * seqScore * ionsA;
//            cout << " score = " << score
//            <<" score_N_C = " <<score_N_C
//            << " score_C_N = " << score_C_N
//            <<" seqScore = " << seqScore
//            <<" ionsA = " << ionsA
//            <<endl;
//            if (fabs(subCor) > 1) {
//                score /= subCor;
//            }
            score *= subCor;
        }


        void get_ncter_score(double mux,double h ,double subCor,double seqScore,double ionsA)
        {
            score = (score_c + score_n) * sqrt(h) * seqScore * ionsA;
//            cout << " score = " << score
//            <<" score_N_C = " <<score_N_C
//            << " score_C_N = " << score_C_N
//            <<" seqScore = " << seqScore
//            <<" ionsA = " << ionsA
//            <<endl;
//            if (fabs(subCor) > 1) {
//                score /= subCor;
//            }
//            cout<<"before_score:"<<score<<endl;
            score *= subCor;
//            cout<<"after_score:"<<score<<endl;
        }





        void get_score_newlength(double mux,double h ,int subCor,double seqScore,double ionsA,double ion_len_pp)
        {
            score = (score_c + score_n) * sqrt(h) * seqScore * ionsA;
//            cout << " score = " << score
//            <<" score_N_C = " <<score_N_C
//            << " score_C_N = " << score_C_N
//            <<" seqScore = " << seqScore
//            <<" ionsA = " << ionsA
//            <<endl;
            if (fabs(subCor) > 1) {
                score /= subCor;
            }
                score *= (ion_len_pp+1);
        }



        void GetMaxSize() {
            MaxSize = alignment_ions_n.size();
        }

        void SetSize();

        void InitArgument(int Ncut,int Ccut ,const vector<double> & TN,const vector<double> & TC) {
            cut_location_n = Ncut ;
            cut_location_c = Ccut ;
            theory_ions_mass_c = TC ;
            theory_ions_mass_n = TN ;
        }

    };
    typedef std::shared_ptr<temp_prsm_argument> temp_prsm_argument_share_ptr;
    typedef temp_prsm_argument* temp_prsm_argument_ptr ;
}
#endif