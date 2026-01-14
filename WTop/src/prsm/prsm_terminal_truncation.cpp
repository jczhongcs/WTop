//
// Created by wenzhong on 2023/3/22.
//

#include "prsm_terminal_truncation.hpp"

int prsm_terminal_truncation::match_peaks_search_ppm(const vector<double> &mono_mass,
                                                     const vector<double> &theo_mass,
                                                     double para_ppm) {
    int count = 0;
    //获取匹配个数
    vector<double> short_1 = theo_mass;
    vector<double> long_1 = mono_mass;
    if (short_1.size() == 0 || long_1.size() == 0) {
        return 0;
    }
    ppm_value_ptr ppmValuePtr = std::make_shared<ppm_value>(para_ppm);
    for (int i = 0; i < short_1.size(); i++) {
        int n = mylib::speculate_lib::txpb_binarey_search_ex(long_1, long_1.size(), short_1[i]);
        if (n >=0 && ppmValuePtr->calculate_ppm_value(short_1[i],long_1[n])) {
            count++;
        }
    }
//    delete ppmValuePtr ;
    return count;
}


int prsm_terminal_truncation::match_peaks_search_value(const vector<double> &mono_mass,
                             const vector<double> &theo_mass,
                             double value) {
    int count = 0;
    //获取匹配个数
    vector<double> short_1 = theo_mass;
    vector<double> long_1 = mono_mass;

    if (short_1.size() == 0 || long_1.size() == 0) {
        return 0;
    }
    for (int i = 0; i < short_1.size(); i++) {
        int n = mylib::speculate_lib::txpb_binarey_search_ex(long_1, long_1.size(), short_1[i]);
        if (fabs((short_1[i]) - (long_1[n])) < value) {
            count++;
        }
    }
    return count;
}





//判断未知修饰最合适的截断
void prsm_terminal_truncation::judge_terminal_truncation_By(prsm::modify_ptr mp,
                                                      protein_ptr proteinPtr,
                                                      msalign_ptr msalignPtr,
                                                      const prsm::config_argument_ptr &prsmArg,
                                                      int &cut_n, int cut_c)
{
    vector<double> mono_mass = msalignPtr->getIonsMassContainer();    //临时记录质谱质量
    double max_ptm_mass = 0;        //记录下最大修饰质量
    double min_ptm_mass = 0;        //记录下最小修饰质量
    double precursor_mass = msalignPtr->getPrecursorMass();

    vector<double> ptm_mass;    //修饰的质量，方便直接取值
    map<double, string> cmap = (*mp).mapi;  //保存修饰map
    vector<double> seq_mass; //记录氨基酸序列质量
    proteinPtr->get_seq_mass(seq_mass); //记录氨基酸质量序列

    //得到修饰质量向量
    for (map<double, string>::iterator it_mp = cmap.begin(); it_mp != cmap.end(); ++it_mp) {
        if (max_ptm_mass < it_mp->first) {
            max_ptm_mass = it_mp->first;
        }
        if (min_ptm_mass > it_mp->first) {
            min_ptm_mass = it_mp->first;
        }
        ptm_mass.push_back(it_mp->first);
    }   //for end
    int cut_n_best = -1;
    double max_theory_mass = proteinPtr->getTheoryMassCBy().back();
    int peaks_match = 0;    //记录下匹配的总峰个数

    //遍历寻找最优的蛋白质形,判断最优截断位置
    for (int co = 0; co < seq_mass.size(); ++co) {
        if (precursor_mass - max_theory_mass > max_ptm_mass + prsmArg->unknownPTMsMassMax) {
            break;
        }
        if (precursor_mass - max_theory_mass < prsmArg->unknownPTMsMassMin + min_ptm_mass) {
            max_theory_mass -= seq_mass[co];
            continue;
        }
        max_theory_mass -= seq_mass[co];
        vector<double> theory_ptm_unknow;   //保存截断后的理论峰值
        //遍历可能的修饰与未知修饰，并且添加修饰峰与未知修饰峰
        //当前截断下，有N种修饰，就有N+1种蛋白质形,取出当前截断下,1条最合适的蛋白质形,并保存
        for (map<double, string>::iterator it_mp = mp->mapi.begin(); it_mp != mp->mapi.end(); ++it_mp) {
//            cout<<"co = " << co << "   " <<it_mp->first << " " << it_mp->second << endl ;
            if (fabs(it_mp->first) < prsmArg->min_value) {
                continue;
            }
            double max_theory_mass_1 = proteinPtr->getTheoryMassCBy().back();
            theory_ptm_unknow.assign(proteinPtr->getTheoryMassNBy().begin() + co, proteinPtr->getTheoryMassNBy().end());
//            cout << theory_ptm_unknow.size() << endl ;
            vector<double> theory_with_ptm;
            vector<double> theory_with_unk;
            vector<double> theory_with_ptm_unk;
            double sub_mass = 0;    //记录截断的氨基酸质量
            for (int i = 0; i < co; ++i) {
                sub_mass += seq_mass[i];
            }
            max_theory_mass_1 -= sub_mass;
            double unk_mass = precursor_mass - max_theory_mass_1 - it_mp->first;    //记录未知修饰质量
            //N端减去
            for (int i = 0; i < theory_ptm_unknow.size(); ++i) {
                theory_ptm_unknow[i] -= sub_mass;
            }
            //峰值序列添加
            for (int iu = 0; iu < theory_ptm_unknow.size(); ++iu) {
                theory_with_unk.push_back(theory_ptm_unknow[iu] + unk_mass);
                if (it_mp->first != 0) {
                    theory_with_ptm.push_back(theory_ptm_unknow[iu] + it_mp->first);
                    theory_with_ptm_unk.push_back(theory_ptm_unknow[iu] + it_mp->first + unk_mass);
                }
            }
            int match_peaks = match_peaks_search_ppm(mono_mass, theory_ptm_unknow,prsmArg->min_ppm_value);
            int ptm_peaks = match_peaks_search_ppm(mono_mass, theory_with_ptm,prsmArg->min_ppm_value);
//                int unk_peaks = match_peaks_search_fa(mono_seq_mass, theory_with_unk, mono_seq_mass.size());
            int unk_peaks = 0 ;
            int ptm_unk_peaks = match_peaks_search_value(mono_mass, theory_with_ptm_unk,0.2);
            int all_peaks = match_peaks + ptm_peaks + unk_peaks + ptm_unk_peaks;

            if (fabs(it_mp->first) < prsmArg->min_value) {
                if (unk_peaks > 0) {
                    if (peaks_match < all_peaks) {
                        peaks_match = all_peaks;
                        cut_n_best = co;
                    }
                }
            } else if (fabs(it_mp->first) > prsmArg->min_value) {
                if (peaks_match < all_peaks) {
                    peaks_match = all_peaks;
                    cut_n_best = co;
                }
            }
//            cout<<endl<<"co:"<<co<<"  unk_mass:"<<unk_mass<<"  match_peak:"<<match_peaks
//                 <<"  ptm_peak:"<<ptm_peaks<<"  ptm_unk_peak:"<<ptm_unk_peaks<<"  all_peaks:"<<all_peaks<<endl;
        }   //遍历修饰
    }  //遍历截断
    cut_n = cut_n_best; //纪录当前最佳截断
}


void prsm_terminal_truncation::judge_terminal_truncation_Cz(prsm::modify_ptr mp,
                                  protein_ptr proteinPtr,
                                  msalign_ptr msalignPtr,
                                  const prsm::config_argument_ptr &prsmArg,
                                  int &cut_n, int cut_c) {
    vector<double> mono_seq_mass = msalignPtr->getIonsMassContainer();    //临时记录质谱质量
    double max_ptm_mass = 0;        //记录下最大修饰质量
    double min_ptm_mass = 0;        //记录下最小修饰质量
//        double subMassETD = 18.01056 - 1.9919;
    double subMassETD = prsmArg->sub_etd;
    double precursor_mass = msalignPtr->getPrecursorMass();
    vector<double> ptm_mass;    //修饰的质量，方便直接取值
    map<double, string> cmap = (*mp).mapi;  //保存修饰map
    vector<double> seq_mass;    //记录氨基酸序列质量
    proteinPtr->get_seq_mass(seq_mass); //记录氨基酸质量序列
    //得到修饰质量向量和最大、最小修饰质量
    for (map<double, string>::iterator it_mp = cmap.begin(); it_mp != cmap.end(); ++it_mp) {
        if (max_ptm_mass < it_mp->first) {
            max_ptm_mass = it_mp->first;
        }
        if (min_ptm_mass > it_mp->first) {
            min_ptm_mass = it_mp->first;
        }
        ptm_mass.push_back(it_mp->first);
    }
    //以插峰的方式,遍历寻找最优的蛋白质形,判断最优截断位置
    int cut_n_best = -1;     //临时变量
    int peaks_match_temp = 0;    //记录最优的匹配的总峰个数
    double max_theory_mass = proteinPtr->getTheoryMassCCz().back();    //从C端判断
    for (int co = 0; co < seq_mass.size(); ++co) {
//            double max_theory_mass_1 = pro->theoryMassC_Cz.back();
        if (precursor_mass - max_theory_mass > max_ptm_mass + prsmArg->unknownPTMsMassMax) {
            break;
        }
        if (precursor_mass - max_theory_mass <  min_ptm_mass + prsmArg->unknownPTMsMassMin) {
            max_theory_mass -= seq_mass[co];
            continue;
        }
        max_theory_mass -= seq_mass[co];     //最大理论质量每次截断后减去对应的氨基酸质量
        vector<double> theory_ptm_unknow;   //保存截断后的理论峰值
        /**
         * 遍历可能的修饰与未知修饰，并且添加修饰峰与未知修饰峰
         * 当前截断下，有N种修饰，就有N+1种蛋白质形,取出当前截断下,1条最合适的蛋白质形,并保存
         */
        for (map<double, string>::iterator it_mp = mp->mapi.begin(); it_mp != mp->mapi.end(); ++it_mp) {
            if (fabs(it_mp->first) < 0.000001) {
                continue;
            }
            double max_theory_mass_1 = proteinPtr->getTheoryMassCCz().back();
            theory_ptm_unknow.assign(proteinPtr->getTheoryMassNCz().begin() + co, proteinPtr->getTheoryMassNCz().end());  // 初始化N端理论峰序列
            vector<double> theory_with_ptm;     //插入ptm之后的序列
            vector<double> theory_with_unk;     //插入未知ptm之后的序列
            vector<double> theory_with_ptm_unk;     //插入混合ptm之后的序列
            double sub_mass = 0;    //记录截断的氨基酸质量
            for (int i = 0; i < co; ++i) {      //得到当前截断下的氨基酸质量,co是当前截断的个数
                sub_mass += seq_mass[i];
            }
            max_theory_mass_1 -= sub_mass;
//                cout << " max_theory_mass_1 = " <<max_theory_mass_1 <<endl;
            double unk_mass = precursor_mass - max_theory_mass_1 - it_mp->first -(18.0106 - 1.9919);    //记录未知修饰质量
            if (unk_mass > prsmArg->unknownPTMsMassMax || unk_mass < prsmArg->unknownPTMsMassMin) {
                continue;
            }
            for (int i = 0; i < theory_ptm_unknow.size() && fabs(sub_mass) > 0.00001; ++i) {    //N端理论峰减去截断质量
                theory_ptm_unknow[i] -= sub_mass;
            }
            for (int iu = 0; iu < theory_ptm_unknow.size(); ++iu) {      //峰值序列添加
                theory_with_unk.push_back(theory_ptm_unknow[iu] + unk_mass);        //未知修饰峰
                if (fabs(it_mp->first) > 0.00001) {     //当只有未知修饰时
                    theory_with_ptm.push_back(theory_ptm_unknow[iu] + it_mp->first);    //已知修饰峰
                    theory_with_ptm_unk.push_back(theory_ptm_unknow[iu] + it_mp->first + unk_mass);     //合并修饰峰
                }
            }
            int match_peaks = match_peaks_search_ppm(mono_seq_mass, theory_ptm_unknow, prsmArg->min_ppm_value);   //判断是否是掉H多H
            int ptm_peaks = match_peaks_search_ppm(mono_seq_mass, theory_with_ptm, prsmArg->min_ppm_value);       //判断带修饰的峰
//                int unk_peaks = match_peaks_search_fa(mono_seq_mass, theory_with_unk, mono_seq_mass.size());    //判断是否在1.5范围内
            int ptm_unk_peaks = match_peaks_search_value(mono_seq_mass, theory_with_ptm_unk,0.2);
            int unk_peaks = 0 ;     //由于截断的质量会添加至未知修饰的质量上，故暂时不考虑只有未知修饰的峰个数
//                int ptm_unk_peaks = 0 ;
            int all_peaks_temp = match_peaks + ptm_peaks + unk_peaks + ptm_unk_peaks;
//                cout <<" co = " << co
//                 << " unknownMss = " << unk_mass
//                 << " ptmMass = " <<it_mp->first
//                 <<  " unknownMss + ptmMass = " << unk_mass + it_mp->first
//                 <<endl;
//                cout << " match_peaks = " << match_peaks
//                << " ptm_peaks = "<< ptm_peaks
//                << " unk_peaks = " <<unk_peaks
//                <<" ptm_unk_peaks = " << ptm_unk_peaks
//                << " all_peaks_temp = " <<all_peaks_temp
//                <<endl<<endl;
            if (fabs(it_mp->first) < prsmArg->min_value) {
                if (unk_peaks > 0) {
                    if (all_peaks_temp > peaks_match_temp) {
                        peaks_match_temp = all_peaks_temp;
                        cut_n_best = co;
                    }
                }
            } else if (fabs(it_mp->first) >prsmArg->min_value) {
                if (all_peaks_temp > peaks_match_temp) {
                    peaks_match_temp = all_peaks_temp;
                    cut_n_best = co;
                }
//                    if (ptm_unk_peaks > 0) {
//                        if (all_peaks_temp > peaks_match_temp) {
//                            peaks_match_temp = all_peaks_temp;
//                            cut_n_best = co;
//                        }
//                    }
            }
        }   //遍历修饰
    }  //遍历截断

    cut_n = cut_n_best; //纪录当前最佳截断

}



///双端截断 by lyc
void prsm_terminal_truncation::test_NC_By(prsm::modify_ptr mp,
                                                     protein_ptr proteinPtr,
                                                     msalign_ptr msalignPtr,
                                                     const prsm::config_argument_ptr &prsmArg,
                                                     map<double,int> &modCutN, map<double,int> &modCutC,map<double,vector<double>> &modUnk)
{

   vector<double> mono_mass = msalignPtr->getIonsMassContainer();    //临时记录质谱质量
    double max_ptm_mass = 0;        //记录下最大修饰质量
    double min_ptm_mass = 0;        //记录下最小修饰质量
    double precursor_mass = msalignPtr->getPrecursorMass();
    vector<double> ptm_mass;    //修饰的质量，方便直接取值
    map<double, string> cmap = (*mp).mapi;  //保存修饰map
    vector<double> seq_mass; //记录氨基酸序列质量
    proteinPtr->get_seq_mass(seq_mass); //记录氨基酸质量序列


    map<double,int> modPeaks;//存放已知修饰对应的当前最多匹配峰数
    map<double,int> curModUnk;//存放已知修饰对应的当前最多匹配峰数时，对应的未知修饰质量
    for(map<double, string>::iterator it_mp = cmap.begin(); it_mp != cmap.end(); ++it_mp){
        modPeaks.insert(pair<double,int>(it_mp->first,0));
        curModUnk.insert(pair<double,int>(it_mp->first,500));
    }

    //得到修饰质量向量
    for (map<double, string>::iterator it_mp = cmap.begin(); it_mp != cmap.end(); ++it_mp) {
        if (max_ptm_mass < it_mp->first) {
            max_ptm_mass = it_mp->first;
        }
        if (min_ptm_mass > it_mp->first) {
            min_ptm_mass = it_mp->first;
        }
        ptm_mass.push_back(it_mp->first);
    }   //for end


    int minCutLength,maxCutLength;
    if((precursor_mass + prsmArg->unknownPTMsMassMin + min_ptm_mass) <= 0)
        minCutLength = 0;
    else{
        minCutLength = (precursor_mass + prsmArg->unknownPTMsMassMin + min_ptm_mass) / 186.079354;
        if((minCutLength - 1) > 0)
            minCutLength--;
    }
    if(seq_mass.size() < minCutLength)
        return;

    maxCutLength = (precursor_mass + prsmArg->unknownPTMsMassMax + max_ptm_mass)/57.02146;
    maxCutLength++;
    if(maxCutLength > seq_mass.size())
        maxCutLength = seq_mass.size();


    int n = seq_mass.size();

    ///遍历子串开头，即0到seqlenth-min
    ///然后动态规划得到每个固定开头，长度min到max的各个字串的峰数
    ///因为只考虑b离子峰数，故有动态规划递推关系，只判定最后一个最大理论峰，是否存在于质谱峰中
    for(int cutStart = 0; cutStart <= seq_mass.size() - minCutLength ; cutStart++){
        if(cutStart == seq_mass.size())
            continue;
        int dplength = maxCutLength - minCutLength + 1;
        if(cutStart + maxCutLength >= seq_mass.size())
            dplength = seq_mass.size() - cutStart - minCutLength + 1;
        if(dplength <= 0) dplength = 1;
        for (map<double, string>::iterator it_mp = mp->mapi.begin(); it_mp != mp->mapi.end(); ++it_mp) {
            for(int cutEnd = (cutStart + minCutLength - 1) ; cutEnd < (cutStart + minCutLength + dplength - 1) ; cutEnd++) {
//                cout<<"cutStart:"<<cutStart<<"  cutEnd:"<<cutEnd<<"  itmp:"<<it_mp->second<<endl;

                std::unique_ptr<ppm_value> ppmValuePtr(new ppm_value(15));
                ///初始化n端peak_select_nter

                vector<double> temp_nter_theory_ptm_unknow;
                vector<double> temp_nter_theory_with_ptm;
                vector<double> temp_nter_theory_with_unk;
                vector<double> temp_nter_theory_with_unk_ptm;

                temp_nter_theory_ptm_unknow.assign(proteinPtr->getTheoryMassNBy().begin() + cutStart, proteinPtr->getTheoryMassNBy().begin() + cutEnd + 1);
                int cuted_seq_length = cutEnd - cutStart + 1;


                ///计算unk需要用到蛋白序列质量，所以用C端的getMass
                double peptide_mass = proteinPtr->getTheoryMassCBy()[seq_mass.size()-cutStart-1];
                double peptide_sub_mass = 0;
                for (int i = seq_mass.size() - 1; i > cutEnd; --i) {
                    peptide_sub_mass += seq_mass[i];
                }
                peptide_mass -= peptide_sub_mass;

                double nter_sub_mass = 0;
                for(int i = 0;i < cutStart;i++){
                    nter_sub_mass += seq_mass[i];
                }

                for (int i = 0; i < temp_nter_theory_ptm_unknow.size(); ++i) {
                    temp_nter_theory_ptm_unknow[i] -= nter_sub_mass;
                }

                ///初始化unk修饰表
                double h_value = 1.007276;
                double ori_unk_mass = 0;
                ori_unk_mass =  precursor_mass - peptide_mass - it_mp->first;
                if(fabs(ori_unk_mass) > prsmArg->unknownPTMsMassMax)
                    continue;
                if(fabs(ori_unk_mass + it_mp->first) < h_value)
                    continue;
                double unk_mass_dH = ori_unk_mass + h_value;
                double unk_mass_sH = ori_unk_mass - h_value;
                vector<double> unkMassVec;
                unkMassVec.push_back(ori_unk_mass);
                unkMassVec.push_back(unk_mass_dH);
                unkMassVec.push_back(unk_mass_sH);
                for(int um = 0;um<unkMassVec.size();um++){
                    if(fabs(unkMassVec[um])<0.115)
                        unkMassVec[um] = 0;
                }
                int cur_tru_max_peaks = 0;
                double cur_unk;
                for(double unk_mass:unkMassVec){
                    if(fabs(unk_mass) > prsmArg->unknownPTMsMassMax + 2) {
                        continue;
                    }
                    for(int i = 0; i < temp_nter_theory_ptm_unknow.size(); ++i){
                        if (fabs(it_mp->first) > prsmArg->min_value) {
                            temp_nter_theory_with_ptm.push_back(temp_nter_theory_ptm_unknow[i] + it_mp->first);
                        }else{
                            temp_nter_theory_with_ptm.push_back(0);
                        }
                        if(fabs(unk_mass) > 0.1) {
                            temp_nter_theory_with_unk.push_back(temp_nter_theory_ptm_unknow[i] + unk_mass);
                        }else{
                            temp_nter_theory_with_unk.push_back(0);
                        }
                        if(fabs(it_mp->first) > prsmArg->min_value && fabs(unk_mass) > 0.1){
                            temp_nter_theory_with_unk_ptm.push_back(temp_nter_theory_ptm_unknow[i] + unk_mass + it_mp->first);
                        }
                        else{
                            temp_nter_theory_with_unk_ptm.push_back(0);
                        }
                    }
                    ///得到C端3种情况下的匹配峰数
                    vector<double> temp_cter_theory_ptm_unknow;
                    vector<double> temp_cter_theory_with_ptm;
                    vector<double> temp_cter_theory_with_unk;
                    vector<double> temp_cter_theory_with_unk_ptm;

                    int pro_seq_size = seq_mass.size();
                    for(int i = pro_seq_size -cutStart-1;i >= pro_seq_size-cutEnd-1;i--){
                        temp_cter_theory_ptm_unknow.push_back(proteinPtr->getTheoryMassCBy()[i]);
                    }


                    ///理论峰减去截断质量
                    double cter_sub_mass = 0;
                    for(int i = seq_mass.size() - 1 ; i > cutEnd; i--)
                        cter_sub_mass += seq_mass[i];
                    for(int i = 0;i < temp_cter_theory_ptm_unknow.size();i++){
                        temp_cter_theory_ptm_unknow[i] -= cter_sub_mass;
                    }


                    for(int i = 0; i < temp_cter_theory_ptm_unknow.size(); ++i){
                        if (fabs(it_mp->first) > prsmArg->min_value) {
                            temp_cter_theory_with_ptm.push_back(temp_cter_theory_ptm_unknow[i] + it_mp->first);
                        }else{
                            temp_cter_theory_with_ptm.push_back(0);
                        }
                        if(fabs(unk_mass) > 0.1) {
                            temp_cter_theory_with_unk.push_back(temp_cter_theory_ptm_unknow[i] + unk_mass);
                        }else{
                            temp_cter_theory_with_unk.push_back(0);
                        }
                        if(fabs(it_mp->first) > prsmArg->min_value && fabs(unk_mass) > 0.1){
                            temp_cter_theory_with_unk_ptm.push_back(temp_cter_theory_ptm_unknow[i] + unk_mass + it_mp->first);
                        }
                        else{
                            temp_cter_theory_with_unk_ptm.push_back(0);
                        }
                    }
                    std::vector<std::vector<int>> peak_select(5, std::vector<int>(cuted_seq_length, 0));

                    int final_peaks=0;
                    for(int i=0; i<cuted_seq_length;i++){
                        ///step1:N端搜完全匹配的峰
                        int ncom_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_nter_theory_ptm_unknow[i]);
                        if (ncom_index >=0 && ppmValuePtr->calculate_ppm_value(temp_nter_theory_ptm_unknow[i],mono_mass[ncom_index])) {
                            final_peaks++;
                        }
                        ///step2：N端搜带ptm的峰
                        int nptm_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_nter_theory_with_ptm[i]);
                        if (nptm_index >=0 && ppmValuePtr->calculate_ppm_value(temp_nter_theory_with_ptm[i],mono_mass[nptm_index])) {
                            final_peaks++;
                        }
                        ///step3：N端搜带unk的峰
                        int nunk_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_nter_theory_with_unk[i]);
                        if (nunk_index >=0 && ppmValuePtr->calculate_ppm_value(temp_nter_theory_with_unk[i],mono_mass[nunk_index])) {
                            final_peaks++;
                        }
                        ///step4：N端搜带unk和ptm的峰
                        int nup_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_nter_theory_with_unk_ptm[i]);
                        if (nup_index >=0 && ppmValuePtr->calculate_ppm_value(temp_nter_theory_with_unk_ptm[i],mono_mass[nup_index])) {
                            final_peaks++;
                        }

                        ///step1:C端搜完全匹配的峰
                        int ccom_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_cter_theory_ptm_unknow[i]);
                        if (ccom_index >=0 && ppmValuePtr->calculate_ppm_value(temp_cter_theory_ptm_unknow[i],mono_mass[ccom_index])) {
                            final_peaks++;
                        }
                        ///step2：C端搜带ptm的峰
                        int cptm_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_cter_theory_with_ptm[i]);
                        if (cptm_index >=0 && ppmValuePtr->calculate_ppm_value(temp_cter_theory_with_ptm[i],mono_mass[cptm_index])) {
                            final_peaks++;
                        }

                        ///step3：C端搜带unk的峰
                        int cunk_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_cter_theory_with_unk[i]);
                        if (cunk_index >=0 && ppmValuePtr->calculate_ppm_value(temp_cter_theory_with_unk[i],mono_mass[cunk_index])) {
                            final_peaks++;
                        }
                        ///step4：C端搜带unk和ptm的峰
                        int cup_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_cter_theory_with_unk_ptm[i]);
                        if (cup_index >=0 && ppmValuePtr->calculate_ppm_value(temp_cter_theory_with_unk_ptm[i],mono_mass[cup_index])) {
                            final_peaks++;
                        }
                    }

                    if(final_peaks > cur_tru_max_peaks) {
                        cur_tru_max_peaks = final_peaks;
                        cur_unk = unk_mass;
                    }
                    ///test out
//                    int first = cutStart + 1;
//                    int last = cutEnd +1;
//                    bool test_out_tru_peak = true;
//                    if(test_out_tru_peak &&it_mp->first ==0){
//                        cout<<"mod:"<<it_mp->first<<" unk:"<<unk_mass<<" first:"<<first<<" last:"<<last<<" peaks:"<<final_peaks<<endl;
//                    }
                    ///end test
                }//for unkMass

                ///每个cut_end比较结果，峰数相同保留未知修饰最小的
                bool ex_result = false;
                if(cur_tru_max_peaks < 5)
                    continue;
                double temp_th;
                int th;
                int cur_cutN = cutStart;
                int cur_cutC = seq_mass.size() - cutEnd  - 1;
                if(modCutC[it_mp->first] == cur_cutC){
                    temp_th = fabs(cur_unk-curModUnk[it_mp->first])/30.0;
                    th = ceil(temp_th) + 1;
                    if(cur_tru_max_peaks - modPeaks[it_mp->first] > th)
                        ex_result = true;
                }
                else if(modCutN[it_mp->first] == cur_cutN){
                    temp_th = fabs(cur_unk-curModUnk[it_mp->first])/30.0;
                    th = ceil(temp_th) + 1;
                    if(cur_tru_max_peaks - modPeaks[it_mp->first] > th)
                        ex_result = true;
                }
                else {

                    if((fabs(cur_unk) - fabs(curModUnk[it_mp->first]))>0.0)
                        temp_th = (fabs(cur_unk) - fabs(curModUnk[it_mp->first]))/30.0;
                    else
                        temp_th = (fabs(cur_unk) - fabs(curModUnk[it_mp->first]))/70.0;
                    if(temp_th > 0){
                        th = ceil(temp_th) + 1;
                    } else if(temp_th == 0){
                        th = 0;
                    } else{
                        th = floor(temp_th) + 1;
                    }
                    if(cur_tru_max_peaks - modPeaks[it_mp->first] > th)
                        ex_result = true;
                }
                if(modCutN[it_mp->first]==-1 && modCutC[it_mp->first]==0)
                    ex_result = true;
                if(ex_result){
                    modPeaks[it_mp->first] = cur_tru_max_peaks;
                    curModUnk[it_mp->first] = cur_unk;
                    modUnk[it_mp->first] = unkMassVec;
                    modCutN[it_mp->first] = cutStart;
                    modCutC[it_mp->first] = seq_mass.size() - cutEnd  - 1;
                    ///test out
//                    if(fabs(it_mp->first-79.966)<0.01 || fabs(it_mp->first-14.02)<0.01){
//                        cout << "itmp:" << it_mp->first << " peak:" << cur_tru_max_peaks << " modCutN:" << modCutN[it_mp->first]
//                             << " modCutC:" << modCutC[it_mp->first] <<" unk:"<<unkMassVec[0] <<endl;
//                    }
                }

            }// for cutEnd
        }// for it_mp 遍历修饰组合
    }// for cutStart

}

void prsm_terminal_truncation::test_NC_Cz(prsm::modify_ptr mp,
                                          protein_ptr proteinPtr,
                                          msalign_ptr msalignPtr,
                                          const prsm::config_argument_ptr &prsmArg,
                                          map<double,int> &modCutN, map<double,int> &modCutC,map<double,vector<double>> &modUnk)
{
    vector<double> mono_mass = msalignPtr->getIonsMassContainer();    //临时记录质谱质量
    double max_ptm_mass = 0;        //记录下最大修饰质量
    double min_ptm_mass = 0;        //记录下最小修饰质量
    double precursor_mass = msalignPtr->getPrecursorMass();

    vector<double> ptm_mass;    //修饰的质量，方便直接取值
    map<double, string> cmap = (*mp).mapi;  //保存修饰map
    vector<double> seq_mass; //记录氨基酸序列质量
    proteinPtr->get_seq_mass(seq_mass); //记录氨基酸质量序列



    map<double,int> modPeaks;//存放已知修饰对应的当前最多匹配峰数
    map<double,int> curModUnk;//存放已知修饰对应的当前最多匹配峰数时，对应的未知修饰质量
    for(map<double, string>::iterator it_mp = cmap.begin(); it_mp != cmap.end(); ++it_mp){
        modPeaks.insert(pair<double,int>(it_mp->first,0));
        curModUnk.insert(pair<double,int>(it_mp->first,500));
    }


    //得到修饰质量向量
    for (map<double, string>::iterator it_mp = cmap.begin(); it_mp != cmap.end(); ++it_mp) {
        if (max_ptm_mass < it_mp->first) {
            max_ptm_mass = it_mp->first;
        }
        if (min_ptm_mass > it_mp->first) {
            min_ptm_mass = it_mp->first;
        }
        ptm_mass.push_back(it_mp->first);
    }   //for end
    int cut_n_best = -1,cut_c_best = 0;

//    int peaks_match = 0;    //记录下匹配的总峰个数


    int minCutLength,maxCutLength;

    if((precursor_mass + prsmArg->unknownPTMsMassMin + min_ptm_mass) <= 0)
        minCutLength = 0;
    else{
        minCutLength = (precursor_mass + prsmArg->unknownPTMsMassMin + min_ptm_mass) / 186.079354;
        if((minCutLength - 1) > 0)
            minCutLength--;
    }
    if(seq_mass.size() < minCutLength)
        return;

    maxCutLength = (precursor_mass + prsmArg->unknownPTMsMassMax + max_ptm_mass)/57.02146;
    maxCutLength++;
    if(maxCutLength > seq_mass.size())
        maxCutLength = seq_mass.size();


    int n = seq_mass.size();

    ///遍历子串开头，即0到seqlenth-min
    ///然后动态规划得到每个固定开头，长度min到max的各个字串的峰数
    ///因为只考虑b离子峰数，故有动态规划递推关系，只判定最后一个最大理论峰，是否存在于质谱峰中
    for(int cutStart = 0; cutStart <= seq_mass.size() - minCutLength ; cutStart++){
        if(cutStart == seq_mass.size())
            continue;
        int dplength = maxCutLength - minCutLength + 1;
        if(cutStart + maxCutLength >= seq_mass.size())
            dplength = seq_mass.size() - cutStart - minCutLength + 1;
        if(dplength <= 0) dplength = 1;
        for (map<double, string>::iterator it_mp = mp->mapi.begin(); it_mp != mp->mapi.end(); ++it_mp) {
            for(int cutEnd = (cutStart + minCutLength - 1) ; cutEnd < (cutStart + minCutLength + dplength - 1) ; cutEnd++) {
                std::unique_ptr<ppm_value> ppmValuePtr(new ppm_value(15));
                ///初始化n端peak_select_nter

                vector<double> temp_nter_theory_ptm_unknow;
                vector<double> temp_nter_theory_with_ptm;
                vector<double> temp_nter_theory_with_unk;
                vector<double> temp_nter_theory_with_unk_ptm;

                temp_nter_theory_ptm_unknow.assign(proteinPtr->getTheoryMassNCz().begin() + cutStart, proteinPtr->getTheoryMassNCz().begin() + cutEnd + 1);
                int cuted_seq_length = cutEnd - cutStart + 1;


                cout<<"getBy.size:"<<proteinPtr->getTheoryMassCBy().size()<<"  seq_mass.size():"<<seq_mass.size()
                    <<"  cutStart:"<<cutStart<<"  seq_mass.size()-cutStart-1:"<<seq_mass.size()-cutStart-1<<endl;

                ///计算unk需要用到蛋白序列质量，所以用C端的getMass
                double peptide_mass = proteinPtr->getTheoryMassCCz()[seq_mass.size()-cutStart-1];
                double peptide_sub_mass = 0;
                for (int i = seq_mass.size() - 1; i > cutEnd; --i) {
                    peptide_sub_mass += seq_mass[i];
                }
                peptide_mass -= peptide_sub_mass;

                double nter_sub_mass = 0;
                for(int i = 0;i < cutStart;i++){
                    nter_sub_mass += seq_mass[i];
                }

                for (int i = 0; i < temp_nter_theory_ptm_unknow.size(); ++i) {
                    temp_nter_theory_ptm_unknow[i] -= nter_sub_mass;
                }

                ///初始化unk修饰表
                double h_value = 1.007276;
                double ori_unk_mass = 0;
                ori_unk_mass =  precursor_mass - peptide_mass - it_mp->first;
                if(fabs(ori_unk_mass) > prsmArg->unknownPTMsMassMax)
                    continue;
                if(fabs(ori_unk_mass + it_mp->first) < h_value)
                    continue;
                double unk_mass_dH = ori_unk_mass + h_value;
                double unk_mass_sH = ori_unk_mass - h_value;
                vector<double> unkMassVec;
                unkMassVec.push_back(ori_unk_mass);
                unkMassVec.push_back(unk_mass_dH);
                unkMassVec.push_back(unk_mass_sH);
                for(int um = 0;um<unkMassVec.size();um++){
                    if(fabs(unkMassVec[um])<0.115)
                        unkMassVec[um] = 0;
                }
                int cur_tru_max_peaks = 0;
                double cur_unk;
                for(double unk_mass:unkMassVec){
                    if(fabs(unk_mass) > prsmArg->unknownPTMsMassMax + 2) {
                        continue;
                    }
                    for(int i = 0; i < temp_nter_theory_ptm_unknow.size(); ++i){
                        if (fabs(it_mp->first) > prsmArg->min_value) {
                            temp_nter_theory_with_ptm.push_back(temp_nter_theory_ptm_unknow[i] + it_mp->first);
                        }else{
                            temp_nter_theory_with_ptm.push_back(0);
                        }
                        if(fabs(unk_mass) > 0.1) {
                            temp_nter_theory_with_unk.push_back(temp_nter_theory_ptm_unknow[i] + unk_mass);
                        }else{
                            temp_nter_theory_with_unk.push_back(0);
                        }
                        if(fabs(it_mp->first) > prsmArg->min_value && fabs(unk_mass) > 0.1){
                            temp_nter_theory_with_unk_ptm.push_back(temp_nter_theory_ptm_unknow[i] + unk_mass + it_mp->first);
                        }
                        else{
                            temp_nter_theory_with_unk_ptm.push_back(0);
                        }
                    }

                    ///得到C端3种情况下的匹配峰数

                    vector<double> temp_cter_theory_ptm_unknow;
                    vector<double> temp_cter_theory_with_ptm;
                    vector<double> temp_cter_theory_with_unk;
                    vector<double> temp_cter_theory_with_unk_ptm;

                    int pro_seq_size = seq_mass.size();
                    for(int i = pro_seq_size -cutStart-1;i >= pro_seq_size-cutEnd-1;i--){
                        temp_cter_theory_ptm_unknow.push_back(proteinPtr->getTheoryMassCCz()[i]);
                    }


                    ///理论峰减去截断质量
                    double cter_sub_mass = 0;
                    for(int i = seq_mass.size() - 1 ; i > cutEnd; i--)
                        cter_sub_mass += seq_mass[i];
                    for(int i = 0;i < temp_cter_theory_ptm_unknow.size();i++){
                        temp_cter_theory_ptm_unknow[i] -= cter_sub_mass;
                    }


                    for(int i = 0; i < temp_cter_theory_ptm_unknow.size(); ++i){
                        if (fabs(it_mp->first) > prsmArg->min_value) {
                            temp_cter_theory_with_ptm.push_back(temp_cter_theory_ptm_unknow[i] + it_mp->first);
                        }else{
                            temp_cter_theory_with_ptm.push_back(0);
                        }
                        if(fabs(unk_mass) > 0.1) {
                            temp_cter_theory_with_unk.push_back(temp_cter_theory_ptm_unknow[i] + unk_mass);
                        }else{
                            temp_cter_theory_with_unk.push_back(0);
                        }
                        if(fabs(it_mp->first) > prsmArg->min_value && fabs(unk_mass) > 0.1){
                            temp_cter_theory_with_unk_ptm.push_back(temp_cter_theory_ptm_unknow[i] + unk_mass + it_mp->first);
                        }
                        else{
                            temp_cter_theory_with_unk_ptm.push_back(0);
                        }
                    }
                    std::vector<std::vector<int>> peak_select(5, std::vector<int>(cuted_seq_length, 0));

                    int final_peaks=0;
                    for(int i=0; i<cuted_seq_length;i++){
                        ///step1:N端搜完全匹配的峰
                        int ncom_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_nter_theory_ptm_unknow[i]);
                        if (ncom_index >=0 && ppmValuePtr->calculate_ppm_value(temp_nter_theory_ptm_unknow[i],mono_mass[ncom_index])) {
                            final_peaks++;
                        }
                        ///step2：N端搜带ptm的峰
                        int nptm_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_nter_theory_with_ptm[i]);
                        if (nptm_index >=0 && ppmValuePtr->calculate_ppm_value(temp_nter_theory_with_ptm[i],mono_mass[nptm_index])) {
                            final_peaks++;
                        }
                        ///step3：N端搜带unk的峰
                        int nunk_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_nter_theory_with_unk[i]);
                        if (nunk_index >=0 && ppmValuePtr->calculate_ppm_value(temp_nter_theory_with_unk[i],mono_mass[nunk_index])) {
                            final_peaks++;
                        }
                        ///step4：N端搜带unk和ptm的峰
                        int nup_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_nter_theory_with_unk_ptm[i]);
                        if (nup_index >=0 && ppmValuePtr->calculate_ppm_value(temp_nter_theory_with_unk_ptm[i],mono_mass[nup_index])) {
                            final_peaks++;
                        }

                        //step1:C端搜完全匹配的峰
                        int ccom_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_cter_theory_ptm_unknow[i]);
                        if (ccom_index >=0 && ppmValuePtr->calculate_ppm_value(temp_cter_theory_ptm_unknow[i],mono_mass[ccom_index])) {
                            final_peaks++;
                        }
                        //step2：C端搜带ptm的峰
                        int cptm_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_cter_theory_with_ptm[i]);
                        if (cptm_index >=0 && ppmValuePtr->calculate_ppm_value(temp_cter_theory_with_ptm[i],mono_mass[cptm_index])) {
                            final_peaks++;
                        }

                        //step3：C端搜带unk的峰
                        int cunk_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_cter_theory_with_unk[i]);
                        if (cunk_index >=0 && ppmValuePtr->calculate_ppm_value(temp_cter_theory_with_unk[i],mono_mass[cunk_index])) {
                            final_peaks++;
                        }
                        //step4：C端搜带unk和ptm的峰
                        int cup_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), temp_cter_theory_with_unk_ptm[i]);
                        if (cup_index >=0 && ppmValuePtr->calculate_ppm_value(temp_cter_theory_with_unk_ptm[i],mono_mass[cup_index])) {
                            final_peaks++;
                        }

                    }



                    if(final_peaks > cur_tru_max_peaks) {
                        cur_tru_max_peaks = final_peaks;
                        cur_unk = unk_mass;
                    }
                    ///test out
//                    int first = cutStart + 1;
//                    int last = cutEnd +1;
//                    bool test_out_tru_peak = true;
//                    if(test_out_tru_peak && last==76 && fabs(it_mp->first-42.01)<0.01){
//                        cout<<"mod:"<<it_mp->first<<" unk:"<<unk_mass<<" first:"<<first<<" last:"<<last<<" peaks:"<<final_peaks<<endl;
//                    }
                    ///end test
                }//for unkMass

                ///每个cut_end比较结果，峰数相同保留未知修饰最小的
                bool ex_result = false;
                if(cur_tru_max_peaks < 5)
                    continue;
                double temp_th;
                int th;
                int cur_cutN = cutStart;
                int cur_cutC = seq_mass.size() - cutEnd  - 1;
                if(modCutC[it_mp->first] == cur_cutC){
                    temp_th = fabs(cur_unk-curModUnk[it_mp->first])/30.0;
                    th = ceil(temp_th) + 1;
                    if(cur_tru_max_peaks - modPeaks[it_mp->first] > th)
                        ex_result = true;
                }
                else if(modCutN[it_mp->first] == cur_cutN){
                    temp_th = fabs(cur_unk-curModUnk[it_mp->first])/30.0;
                    th = ceil(temp_th) + 1;
                    if(cur_tru_max_peaks - modPeaks[it_mp->first] > th)
                        ex_result = true;
                }
                else {
                    temp_th = (fabs(cur_unk) - fabs(curModUnk[it_mp->first]))/30.0;
                    if(temp_th > 0){
                        th = ceil(temp_th) + 1;
                    } else if(temp_th == 0){
                        th = 0;
                    } else{
                        th = floor(temp_th) - 1;
                    }
                    if(cur_tru_max_peaks - modPeaks[it_mp->first] > th)
                        ex_result = true;
                }
                if(modCutN[it_mp->first]==-1 && modCutC[it_mp->first]==0)
                    ex_result = true;
                if(ex_result){
                    modPeaks[it_mp->first] = cur_tru_max_peaks;
                    curModUnk[it_mp->first] = cur_unk;
                    modUnk[it_mp->first] = unkMassVec;
                    modCutN[it_mp->first] = cutStart;
                    modCutC[it_mp->first] = seq_mass.size() - cutEnd  - 1;
                    ///test out
//                    if(fabs(it_mp->first-79.966)<0.01 || fabs(it_mp->first-14.02)<0.01){
//                        cout << "itmp:" << it_mp->first << " peak:" << cur_tru_max_peaks << " modCutN:" << modCutN[it_mp->first]
//                             << " modCutC:" << modCutC[it_mp->first] <<" unk:"<<unkMassVec[0] <<endl;
//                    }
                }

            }// for cutEnd
        }// for it_mp 遍历修饰组合
    }// for cutStart


}



///双端截断 by lyc
void prsm_terminal_truncation::locate_truncation_by(prsm::modify_ptr mp,
                                                    protein_ptr proteinPtr,
                                                    msalign_ptr msalignPtr,
                                                    const prsm::config_argument_ptr &prsmArg,
                                                    map<double,candi_trucation_vec> &ptm_candiTru_map)
{

    ///初始化
    vector<double> mono_mass = msalignPtr->getIonsMassContainer();
    double max_ptm_mass = 0;
    double min_ptm_mass = 0;
    double precursor_mass = msalignPtr->getPrecursorMass();
    vector<double> ptm_mass;
    map<double, string> cmap = (*mp).mapi;



    ///初始化理论质量
    vector<double> seq_mass;
    proteinPtr->get_seq_mass(seq_mass);
    vector<double> prefix_mass(seq_mass.size()+1, 0), suffix_mass(seq_mass.size()+1, 0);
    for (int i = 0; i < seq_mass.size(); ++i) prefix_mass[i+1] = prefix_mass[i] + seq_mass[i];
    for (int i = seq_mass.size()-1; i >= 0; --i) suffix_mass[i] = suffix_mass[i+1] + seq_mass[i];




    ///初始化已知修饰
    for (map<double, string>::iterator it_mp = cmap.begin(); it_mp != cmap.end(); ++it_mp) {
        if (max_ptm_mass < it_mp->first) {max_ptm_mass = it_mp->first;}
        if (min_ptm_mass > it_mp->first) {min_ptm_mass = it_mp->first;}
        ptm_mass.push_back(it_mp->first);
    }



    ///初始化截断长度
    int minCutLength,maxCutLength;
    if((precursor_mass + prsmArg->unknownPTMsMassMin + min_ptm_mass) <= 0)
        minCutLength = 0;
    else{
        minCutLength = (precursor_mass + prsmArg->unknownPTMsMassMin + min_ptm_mass) / 186.079354;
        if((minCutLength) > 1)  minCutLength--;
    }
    if(minCutLength < 1) minCutLength = 1;
    if(seq_mass.size() < minCutLength)  return;
    maxCutLength = (precursor_mass + prsmArg->unknownPTMsMassMax + max_ptm_mass)/57.02146 + 1;
    if(maxCutLength > seq_mass.size())  maxCutLength = seq_mass.size();



    ///初始化中间结果
    map<double,int> modFrg;//存放已知修饰对应的当前最多匹配峰数
    for(map<double, string>::iterator it_mp = cmap.begin(); it_mp != cmap.end(); ++it_mp){
        modFrg.insert(pair<double,int>(it_mp->first,0));
    }



    ///遍历截断位置
    int n = seq_mass.size();
    for(int cutStart = 0; cutStart <= seq_mass.size() - minCutLength ; cutStart++){
        ///边界条件
        if(cutStart == seq_mass.size())
            continue;
        int dplength = maxCutLength - minCutLength + 1;
        if(cutStart + maxCutLength >= seq_mass.size())
            dplength = seq_mass.size() - cutStart - minCutLength + 1;
        if(dplength <= 0) dplength = 1;


        ///初始化N端截断后的理论质量，B离子直接用氨基酸质量累加，Y离子多加水
        vector<double> Nter_mass;
        for(int i=cutStart; i<n; ++i)
            Nter_mass.push_back(proteinPtr->getTheoryMassNBy()[i] - prefix_mass[cutStart]);


        ///对每种已知修饰遍历C端截断
        for (map<double, string>::iterator it_mp = mp->mapi.begin(); it_mp != mp->mapi.end(); ++it_mp) {



            ///对每种修饰遍历C端截断位置
            for(int cutEnd = (cutStart + minCutLength - 1) ; cutEnd < (cutStart + minCutLength + dplength - 1) ; cutEnd++) {
                std::vector<double> Nter_mass_ori(Nter_mass);
                vector<double> Nter_mass_ptm(Nter_mass.size(),0.0);
                for(int i=0; i<Nter_mass.size(); ++i) Nter_mass_ptm[i] = Nter_mass[i] + it_mp->first;


                ///计算unk需要用到蛋白序列质量，所以用C端的getMass
                double peptide_mass = proteinPtr->getTheoryMassCBy()[n-cutStart-1];
                double peptide_sub_mass =suffix_mass[cutEnd+1];
                peptide_mass -= peptide_sub_mass;
                double nter_sub_mass = 0;
                if(cutStart != 0) nter_sub_mass = prefix_mass[cutStart-1];


                ///初始化unk
                double h_value = 1.007276;
                double unk_mass = 0;
                unk_mass =  precursor_mass - peptide_mass - it_mp->first;
                if(fabs(unk_mass) > prsmArg->unknownPTMsMassMax)    continue;
                if(fabs(unk_mass)<1.2)  unk_mass = 0;
                if(fabs(unk_mass) > prsmArg->unknownPTMsMassMax + 2)    continue;
                if(fabs(unk_mass + it_mp->first) < h_value && fabs(unk_mass) > 0.01 )     continue;


                ///根据C端截断位置，初始化N端理论质量
                std::unique_ptr<ppm_value> ppmValuePtr(new ppm_value(15));
                vector<double> Nter_mass_unk;
                vector<double> Nter_mass_unk_ptm;
                if(cutEnd < n - 1) {
                    Nter_mass_ori.erase(Nter_mass_ori.begin() + cutEnd - cutStart + 1, Nter_mass_ori.end());
                    Nter_mass_ptm.erase(Nter_mass_ptm.begin() + cutEnd - cutStart + 1, Nter_mass_ptm.end());
                }
                for(int i = 0; i < Nter_mass_ori.size(); ++i){
                    if(fabs(unk_mass) > 0.1) {
                        Nter_mass_unk.push_back(Nter_mass_ori[i] + unk_mass);
                    }else{
                        Nter_mass_unk.push_back(Nter_mass_ori[i]);
                    }
                    if(fabs(it_mp->first) > prsmArg->min_value || fabs(unk_mass) > 0.1){
                        Nter_mass_unk_ptm.push_back(Nter_mass_ori[i] + unk_mass + it_mp->first);
                    }
                    else{
                        Nter_mass_unk_ptm.push_back(Nter_mass_ori[i]);
                    }
                }


                ///根据C端截断位置，初始化C端理论质量
                vector<double> Cter_mass_ori;
                vector<double> Cter_mass_ptm;
                vector<double> Cter_mass_unk;
                vector<double> Cter_mass_unk_ptm;
                const auto& theoryMassCBy = proteinPtr->getTheoryMassCBy();
                Cter_mass_ori.reserve(cutEnd - cutStart + 1);
                Cter_mass_ori.insert(Cter_mass_ori.end(),theoryMassCBy.rbegin() + cutStart,theoryMassCBy.rbegin() + cutEnd + 1);
                double cter_sub_mass = suffix_mass[cutEnd+1];
                for(int i = 0;i < Cter_mass_ori.size();i++){
                    Cter_mass_ori[i] -= cter_sub_mass;
                }
                for(int i = 0; i < Cter_mass_ori.size(); ++i){
                    if (fabs(it_mp->first) > prsmArg->min_value) {
                        Cter_mass_ptm.push_back(Cter_mass_ori[i] + it_mp->first);
                    }else{
                        Cter_mass_ptm.push_back(Cter_mass_ori[i]);
                    }
                    if(fabs(unk_mass) > 0.1) {
                        Cter_mass_unk.push_back(Cter_mass_ori[i] + unk_mass);
                    }else{
                        Cter_mass_unk.push_back(Cter_mass_ori[i]);
                    }
                    if(fabs(it_mp->first) > prsmArg->min_value || fabs(unk_mass) > 0.1){
                        Cter_mass_unk_ptm.push_back(Cter_mass_ori[i] + unk_mass + it_mp->first);
                    }
                    else{
                        Cter_mass_unk_ptm.push_back(Cter_mass_ori[i]);
                    }
                }



                ///初始化匹配质量数
                int cuted_seq_length = cutEnd - cutStart + 1;
                vector<double> Nter_ori_count(cuted_seq_length,0);
                vector<double> Nter_ptm_count(cuted_seq_length,0);
                vector<double> Nter_unk_count(cuted_seq_length,0);
                vector<double> Nter_unk_ptm_count(cuted_seq_length,0);
                vector<double> Cter_ori_count(cuted_seq_length,0);
                vector<double> Cter_ptm_count(cuted_seq_length,0);
                vector<double> Cter_unk_count(cuted_seq_length,0);
                vector<double> Cter_unk_ptm_count(cuted_seq_length,0);



                ///不同修饰情况的匹配数量
                for(int i=0; i<cuted_seq_length;i++){
                    int ncom_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Nter_mass_ori[i]);
                    if (ncom_index >=0 && ppmValuePtr->calculate_ppm_value(Nter_mass_ori[i],mono_mass[ncom_index])) {
                        if(i==0)
                            Nter_ori_count[0] = 1;
                        else
                            Nter_ori_count[i] = Nter_ori_count[i-1]+1;
                    }else{
                        if(i==0)    Nter_ori_count[0] = 0;
                        else    Nter_ori_count[i] = Nter_ori_count[i-1];
                    }
                    int nptm_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Nter_mass_ptm[i]);
                    if (nptm_index >=0 && ppmValuePtr->calculate_ppm_value(Nter_mass_ptm[i],mono_mass[nptm_index])) {
                        if(i==0)
                            Nter_ptm_count[0] = 1;
                        else{
                            Nter_ptm_count[i] = Nter_ptm_count[i-1]+1;
                        }
                    }else{
                        if(i==0)    Nter_ptm_count[0] = 0;
                        else    Nter_ptm_count[i] = Nter_ptm_count[i-1];
                    }
                    int nunk_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Nter_mass_unk[i]);
                    if (nunk_index >=0 && ppmValuePtr->calculate_ppm_value(Nter_mass_unk[i],mono_mass[nunk_index])) {
                        if(i==0)
                            Nter_unk_count[0] = 1;
                        else{
                            Nter_unk_count[i] = Nter_unk_count[i-1]+1;
                        }
                    }else{
                        if(i==0)    Nter_unk_count[0] = 0;
                        else    Nter_unk_count[i] = Nter_unk_count[i-1];
                    }
                    int nunk_ptm_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Nter_mass_unk_ptm[i]);
                    if (nunk_ptm_index >=0 && ppmValuePtr->calculate_ppm_value(Nter_mass_unk_ptm[i],mono_mass[nunk_ptm_index])) {
                        if(i==0)
                            Nter_unk_ptm_count[0] = 1;
                        else{
                            Nter_unk_ptm_count[i] = Nter_unk_ptm_count[i-1]+1;
                        }
                    }else{
                        if(i==0)    Nter_unk_ptm_count[0] = 0;
                        else    Nter_unk_ptm_count[i] = Nter_unk_ptm_count[i-1];
                    }
                    int ccom_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Cter_mass_ori[cuted_seq_length-1-i]);
                    if (ccom_index >=0 && ppmValuePtr->calculate_ppm_value(Cter_mass_ori[cuted_seq_length-1-i],mono_mass[ccom_index])) {
                        if(i==0)
                            Cter_ori_count[cuted_seq_length-1] = 1;
                        else{
                            Cter_ori_count[cuted_seq_length-1-i] = Cter_ori_count[cuted_seq_length-1-i+1]+1;
                        }
                    }else{
                        if(i==0)    Cter_ori_count[cuted_seq_length-1] = 0;
                        else     Cter_ori_count[cuted_seq_length-1-i] = Cter_ori_count[cuted_seq_length-1-i+1];
                    }
                    int cptm_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Cter_mass_ptm[cuted_seq_length-1-i]);
                    if (cptm_index >=0 && ppmValuePtr->calculate_ppm_value(Cter_mass_ptm[cuted_seq_length-1-i],mono_mass[cptm_index])) {
                        if(i==0)
                            Cter_ptm_count[cuted_seq_length-1] = 1;
                        else{
                            Cter_ptm_count[cuted_seq_length-1-i] = Cter_ptm_count[cuted_seq_length-1-i+1]+1;
                        }
                    }else{
                        if(i==0)    Cter_ptm_count[cuted_seq_length-1] = 0;
                        else     Cter_ptm_count[cuted_seq_length-1-i] = Cter_ptm_count[cuted_seq_length-1-i+1];
                    }
                    int cunk_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Cter_mass_unk[cuted_seq_length-1-i]);
                    if (cunk_index >=0 && ppmValuePtr->calculate_ppm_value(Cter_mass_unk[cuted_seq_length-1-i],mono_mass[cunk_index])) {
                        if(i==0)
                            Cter_unk_count[cuted_seq_length-1] = 1;
                        else{
                            Cter_unk_count[cuted_seq_length-1-i] = Cter_unk_count[cuted_seq_length-1-i+1]+1;
                        }
                    }else{
                        if(i==0)    Cter_unk_count[cuted_seq_length-1] = 0;
                        else     Cter_unk_count[cuted_seq_length-1-i] = Cter_unk_count[cuted_seq_length-1-i+1];
                    }
                    int cunk_ptm_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Cter_mass_unk_ptm[cuted_seq_length-1-i]);
                    if (cunk_ptm_index >=0 && ppmValuePtr->calculate_ppm_value(Cter_mass_unk_ptm[cuted_seq_length-1-i],mono_mass[cunk_ptm_index])) {
                        if(i==0)
                            Cter_unk_ptm_count[cuted_seq_length-1] = 1;
                        else{
                            Cter_unk_ptm_count[cuted_seq_length-1-i] = Cter_unk_ptm_count[cuted_seq_length-1-i+1]+1;
                        }
                    }else{
                        if(i==0)    Cter_unk_ptm_count[cuted_seq_length-1] = 0;
                        else     Cter_unk_ptm_count[cuted_seq_length-1-i] = Cter_unk_ptm_count[cuted_seq_length-1-i+1];
                    }
                }




                ///进行k均分假设
                int k = 20;
                int span = cuted_seq_length/k;
//                if(span <= 0)
//                    span = 1;
                int ptm_loc= 0,unk_loc = 0;
                int frg = 0,pre = 0,suf = 0,match_peak = 0;
                int cur_max_frg = 0,cur_max_pre = 0,cur_max_suf = 0,cur_max_match_peak = 0;
                while (ptm_loc <= cuted_seq_length && span != 0){
                    unk_loc = 0;
                    while(unk_loc <= cuted_seq_length){
                        ///判定部分
                        if(unk_loc == ptm_loc) {
                            if(unk_loc == cuted_seq_length-1)   break;
                            unk_loc += span;
                            if(unk_loc > cuted_seq_length-1)    unk_loc = cuted_seq_length-1;
                            continue;
                        }


                        ///计算匹配数
                        if(unk_loc < ptm_loc){
                            frg = Nter_unk_count[ptm_loc] - Nter_unk_count[unk_loc] + Cter_ptm_count[unk_loc] - Cter_ptm_count[ptm_loc];
                            pre = Nter_ori_count[unk_loc];
                            suf = Cter_ori_count[ptm_loc];
                            match_peak = Nter_ori_count[unk_loc] + Nter_unk_count[ptm_loc] - Nter_unk_count[unk_loc] + Nter_unk_ptm_count[Nter_unk_ptm_count.size()-1] -Nter_unk_ptm_count[ptm_loc]
                                         + Cter_ori_count[ptm_loc] + Cter_ptm_count[unk_loc] - Cter_ptm_count[ptm_loc] + Cter_unk_ptm_count[0] - Cter_unk_ptm_count[unk_loc];
                        }
                        else{
                            frg = Nter_ptm_count[unk_loc] - Nter_ptm_count[ptm_loc] + Cter_ptm_count[ptm_loc] - Cter_ptm_count[unk_loc];
                            pre = Nter_ori_count[ptm_loc];
                            suf = Cter_ori_count[unk_loc];
                            match_peak = Nter_ori_count[ptm_loc] + Nter_ptm_count[unk_loc] - Nter_ptm_count[ptm_loc] + Nter_unk_ptm_count[Nter_unk_ptm_count.size()-1] -Nter_unk_ptm_count[unk_loc]
                                         +Cter_ori_count[unk_loc] + Cter_unk_count[ptm_loc] - Cter_unk_count[unk_loc] + Cter_unk_ptm_count[0] - Cter_unk_ptm_count[ptm_loc];
                        }

                        if(frg > cur_max_frg){cur_max_frg = frg;}
                        if(pre > cur_max_pre){cur_max_pre = pre;}
                        if(suf > cur_max_suf){cur_max_suf = suf;}
                        if(match_peak > cur_max_match_peak){cur_max_match_peak = match_peak;}
                        ///结束条件
                        if(unk_loc == cuted_seq_length-1)   break;
                        unk_loc += span;
                        if(unk_loc > cuted_seq_length-1)    unk_loc = cuted_seq_length-1;
                    }
                    ///结束条件
                    if(ptm_loc == cuted_seq_length-1)
                        break;
                    ptm_loc += span;
                    if(ptm_loc > cuted_seq_length-1)
                        ptm_loc = cuted_seq_length-1;
                }
                if(cur_max_frg > ptm_candiTru_map[it_mp->first][0]->getMatchPeaks()) {
                    ptm_candiTru_map[it_mp->first][0]->setMatchPeaks(cur_max_frg);
                    ptm_candiTru_map[it_mp->first][0]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][0]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][0]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][0]->setCutC(seq_mass.size() - cutEnd  - 1);
                }else if(cur_max_frg == ptm_candiTru_map[it_mp->first][0]->getMatchPeaks() && unk_mass == 0){
                    ptm_candiTru_map[it_mp->first][0]->setMatchPeaks(cur_max_frg);
                    ptm_candiTru_map[it_mp->first][0]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][0]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][0]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][0]->setCutC(seq_mass.size() - cutEnd  - 1);
                }
                if(cur_max_pre > ptm_candiTru_map[it_mp->first][1]->getMatchPeaks()) {
                    ptm_candiTru_map[it_mp->first][1]->setMatchPeaks(cur_max_pre);
                    ptm_candiTru_map[it_mp->first][1]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][1]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][1]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][1]->setCutC(seq_mass.size() - cutEnd  - 1);
                }else if(cur_max_pre == ptm_candiTru_map[it_mp->first][1]->getMatchPeaks() && unk_mass == 0){
                    ptm_candiTru_map[it_mp->first][1]->setMatchPeaks(cur_max_pre);
                    ptm_candiTru_map[it_mp->first][1]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][1]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][1]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][1]->setCutC(seq_mass.size() - cutEnd  - 1);
                }
                if(cur_max_suf > ptm_candiTru_map[it_mp->first][2]->getMatchPeaks()) {
                    ptm_candiTru_map[it_mp->first][2]->setMatchPeaks(cur_max_suf);
                    ptm_candiTru_map[it_mp->first][2]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][2]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][2]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][2]->setCutC(seq_mass.size() - cutEnd  - 1);
                } else if(cur_max_suf == ptm_candiTru_map[it_mp->first][2]->getMatchPeaks() && unk_mass == 0){
                    ptm_candiTru_map[it_mp->first][2]->setMatchPeaks(cur_max_suf);
                    ptm_candiTru_map[it_mp->first][2]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][2]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][2]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][2]->setCutC(seq_mass.size() - cutEnd  - 1);
                }
                if(cur_max_match_peak > ptm_candiTru_map[it_mp->first][3]->getMatchPeaks()) {
                    ptm_candiTru_map[it_mp->first][3]->setMatchPeaks(cur_max_match_peak);
                    ptm_candiTru_map[it_mp->first][3]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][3]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][3]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][3]->setCutC(seq_mass.size() - cutEnd  - 1);
                }else if(cur_max_suf == ptm_candiTru_map[it_mp->first][2]->getMatchPeaks() && unk_mass == 0){
                    ptm_candiTru_map[it_mp->first][3]->setMatchPeaks(cur_max_match_peak);
                    ptm_candiTru_map[it_mp->first][3]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][3]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][3]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][3]->setCutC(seq_mass.size() - cutEnd  - 1);
                }
                Nter_mass_ori.clear();
                Nter_mass_ptm.clear();
                Nter_mass_unk.clear();
                Nter_mass_unk_ptm.clear();
                Cter_mass_ori.clear();
                Cter_mass_ptm.clear();
                Cter_mass_unk.clear();
                Cter_mass_unk_ptm.clear();
                Nter_ori_count.clear();
                Nter_ptm_count.clear();
                Nter_unk_count.clear();
                Nter_unk_ptm_count.clear();
                Cter_ori_count.clear();
                Cter_ptm_count.clear();
                Cter_unk_count.clear();
                Cter_unk_ptm_count.clear();

            }// for cutEnd
        }// for it_mp 遍历修饰组合
    }// for cutStart

}



void prsm_terminal_truncation::locate_truncation_cz(prsm::modify_ptr mp,
                                                    protein_ptr proteinPtr,
                                                    msalign_ptr msalignPtr,
                                                    const prsm::config_argument_ptr &prsmArg,
                                                    map<double,candi_trucation_vec> &ptm_candiTru_map)
{

    ///初始化
    vector<double> mono_mass = msalignPtr->getIonsMassContainer();
    double max_ptm_mass = 0;
    double min_ptm_mass = 0;
    double precursor_mass = msalignPtr->getPrecursorMass();
    vector<double> ptm_mass;
    map<double, string> cmap = (*mp).mapi;



    ///初始化理论质量
    vector<double> seq_mass;
    proteinPtr->get_seq_mass(seq_mass);
    vector<double> prefix_mass(seq_mass.size()+1, 0), suffix_mass(seq_mass.size()+1, 0);
    for (int i = 0; i < seq_mass.size(); ++i) prefix_mass[i+1] = prefix_mass[i] + seq_mass[i];
    for (int i = seq_mass.size()-1; i >= 0; --i) suffix_mass[i] = suffix_mass[i+1] + seq_mass[i];




    ///初始化已知修饰
    for (map<double, string>::iterator it_mp = cmap.begin(); it_mp != cmap.end(); ++it_mp) {
        if (max_ptm_mass < it_mp->first) {max_ptm_mass = it_mp->first;}
        if (min_ptm_mass > it_mp->first) {min_ptm_mass = it_mp->first;}
        ptm_mass.push_back(it_mp->first);
    }



    ///初始化截断长度
    int minCutLength,maxCutLength;
    if((precursor_mass + prsmArg->unknownPTMsMassMin + min_ptm_mass) <= 0)
        minCutLength = 0;
    else{
        minCutLength = (precursor_mass + prsmArg->unknownPTMsMassMin + min_ptm_mass) / 186.079354;
        if((minCutLength) > 1)  minCutLength--;
    }
    if(minCutLength < 1) minCutLength = 1;
    if(seq_mass.size() < minCutLength)  return;
    maxCutLength = (precursor_mass + prsmArg->unknownPTMsMassMax + max_ptm_mass)/57.02146 + 1;
    if(maxCutLength > seq_mass.size())  maxCutLength = seq_mass.size();



    ///初始化中间结果
    map<double,int> modFrg;//存放已知修饰对应的当前最多匹配峰数
    for(map<double, string>::iterator it_mp = cmap.begin(); it_mp != cmap.end(); ++it_mp){
        modFrg.insert(pair<double,int>(it_mp->first,0));
    }



    ///遍历截断位置
    int n = seq_mass.size();
    for(int cutStart = 0; cutStart <= seq_mass.size() - minCutLength ; cutStart++){
        ///边界条件
        if(cutStart == seq_mass.size())
            continue;
        int dplength = maxCutLength - minCutLength + 1;
        if(cutStart + maxCutLength >= seq_mass.size())
            dplength = seq_mass.size() - cutStart - minCutLength + 1;
        if(dplength <= 0) dplength = 1;


        ///初始化N端截断后的理论质量，B离子直接用氨基酸质量累加，Y离子多加水
        vector<double> Nter_mass;
        for(int i=cutStart; i<n; ++i)
            Nter_mass.push_back(proteinPtr->getTheoryMassNCz()[i] - prefix_mass[cutStart]);


        ///对每种已知修饰遍历C端截断
        for (map<double, string>::iterator it_mp = mp->mapi.begin(); it_mp != mp->mapi.end(); ++it_mp) {



            ///对每种修饰遍历C端截断位置
            for(int cutEnd = (cutStart + minCutLength - 1) ; cutEnd < (cutStart + minCutLength + dplength - 1) ; cutEnd++) {
                std::vector<double> Nter_mass_ori(Nter_mass);
                ///初始化N端截断后的带PTM的质量
                vector<double> Nter_mass_ptm(Nter_mass.size(),0.0);
                for(int i=0; i<Nter_mass.size(); ++i) Nter_mass_ptm[i] = Nter_mass[i] + it_mp->first;


                ///计算unk需要用到蛋白序列质量，所以用C端的getMass
                double peptide_mass = proteinPtr->getTheoryMassCCz()[n-cutStart-1];
                double peptide_sub_mass =suffix_mass[cutEnd+1];
                peptide_mass -= peptide_sub_mass;
                double nter_sub_mass = 0;
                if(cutStart != 0) nter_sub_mass = prefix_mass[cutStart-1];


                ///初始化unk
                double h_value = 1.007276;
                double unk_mass = 0;
                unk_mass =  precursor_mass - peptide_mass - it_mp->first;
                if(fabs(unk_mass) > prsmArg->unknownPTMsMassMax)    continue;
                if(fabs(unk_mass)<1.2)  unk_mass = 0;
                if(fabs(unk_mass) > prsmArg->unknownPTMsMassMax + 2)    continue;
                if(fabs(unk_mass + it_mp->first) < h_value && fabs(unk_mass) > 0.01 )     continue;


                ///根据C端截断位置，初始化N端理论质量
                std::unique_ptr<ppm_value> ppmValuePtr(new ppm_value(15));
                vector<double> Nter_mass_unk;
                vector<double> Nter_mass_unk_ptm;
                if(cutEnd < n - 1) {
                    Nter_mass_ori.erase(Nter_mass_ori.begin() + cutEnd - cutStart + 1, Nter_mass_ori.end());
                    Nter_mass_ptm.erase(Nter_mass_ptm.begin() + cutEnd - cutStart + 1, Nter_mass_ptm.end());
                }
                for(int i = 0; i < Nter_mass_ori.size(); ++i){
                    if(fabs(unk_mass) > 0.1) {
                        Nter_mass_unk.push_back(Nter_mass_ori[i] + unk_mass);
                    }else{
                        Nter_mass_unk.push_back(Nter_mass_ori[i]);
                    }
                    if(fabs(it_mp->first) > prsmArg->min_value || fabs(unk_mass) > 0.1){
                        Nter_mass_unk_ptm.push_back(Nter_mass_ori[i] + unk_mass + it_mp->first);
                    }
                    else{
                        Nter_mass_unk_ptm.push_back(Nter_mass_ori[i]);
                    }
                }


                ///根据C端截断位置，初始化C端理论质量
                vector<double> Cter_mass_ori;
                vector<double> Cter_mass_ptm;
                vector<double> Cter_mass_unk;
                vector<double> Cter_mass_unk_ptm;
                const auto& theoryMassCCz = proteinPtr->getTheoryMassCCz();
                Cter_mass_ori.reserve(cutEnd - cutStart + 1);
                Cter_mass_ori.insert(Cter_mass_ori.end(),theoryMassCCz.rbegin() + cutStart,theoryMassCCz.rbegin() + cutEnd + 1);
                double cter_sub_mass = suffix_mass[cutEnd+1];
                for(int i = 0;i < Cter_mass_ori.size();i++){
                    Cter_mass_ori[i] -= cter_sub_mass;
                }
                for(int i = 0; i < Cter_mass_ori.size(); ++i){
                    if (fabs(it_mp->first) > prsmArg->min_value) {
                        Cter_mass_ptm.push_back(Cter_mass_ori[i] + it_mp->first);
                    }else{
                        Cter_mass_ptm.push_back(Cter_mass_ori[i]);
                    }
                    if(fabs(unk_mass) > 0.1) {
                        Cter_mass_unk.push_back(Cter_mass_ori[i] + unk_mass);
                    }else{
                        Cter_mass_unk.push_back(Cter_mass_ori[i]);
                    }
                    if(fabs(it_mp->first) > prsmArg->min_value || fabs(unk_mass) > 0.1){
                        Cter_mass_unk_ptm.push_back(Cter_mass_ori[i] + unk_mass + it_mp->first);
                    }
                    else{
                        Cter_mass_unk_ptm.push_back(Cter_mass_ori[i]);
                    }
                }



                ///初始化匹配质量数
                int cuted_seq_length = cutEnd - cutStart + 1;
                vector<double> Nter_ori_count(cuted_seq_length,0);
                vector<double> Nter_ptm_count(cuted_seq_length,0);
                vector<double> Nter_unk_count(cuted_seq_length,0);
                vector<double> Nter_unk_ptm_count(cuted_seq_length,0);
                vector<double> Cter_ori_count(cuted_seq_length,0);
                vector<double> Cter_ptm_count(cuted_seq_length,0);
                vector<double> Cter_unk_count(cuted_seq_length,0);
                vector<double> Cter_unk_ptm_count(cuted_seq_length,0);



                ///不同修饰情况的匹配数量
                for(int i=0; i<cuted_seq_length;i++){
                    int ncom_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Nter_mass_ori[i]);
                    if (ncom_index >=0 && ppmValuePtr->calculate_ppm_value(Nter_mass_ori[i],mono_mass[ncom_index])) {
                        if(i==0)
                            Nter_ori_count[0] = 1;
                        else
                            Nter_ori_count[i] = Nter_ori_count[i-1]+1;
                    }else{
                        if(i==0)    Nter_ori_count[0] = 0;
                        else    Nter_ori_count[i] = Nter_ori_count[i-1];
                    }
                    int nptm_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Nter_mass_ptm[i]);
                    if (nptm_index >=0 && ppmValuePtr->calculate_ppm_value(Nter_mass_ptm[i],mono_mass[nptm_index])) {
                        if(i==0)
                            Nter_ptm_count[0] = 1;
                        else{
                            Nter_ptm_count[i] = Nter_ptm_count[i-1]+1;
                        }
                    }else{
                        if(i==0)    Nter_ptm_count[0] = 0;
                        else    Nter_ptm_count[i] = Nter_ptm_count[i-1];
                    }
                    int nunk_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Nter_mass_unk[i]);
                    if (nunk_index >=0 && ppmValuePtr->calculate_ppm_value(Nter_mass_unk[i],mono_mass[nunk_index])) {
                        if(i==0)
                            Nter_unk_count[0] = 1;
                        else{
                            Nter_unk_count[i] = Nter_unk_count[i-1]+1;
                        }
                    }else{
                        if(i==0)    Nter_unk_count[0] = 0;
                        else    Nter_unk_count[i] = Nter_unk_count[i-1];
                    }
                    int nunk_ptm_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Nter_mass_unk_ptm[i]);
                    if (nunk_ptm_index >=0 && ppmValuePtr->calculate_ppm_value(Nter_mass_unk_ptm[i],mono_mass[nunk_ptm_index])) {
                        if(i==0)
                            Nter_unk_ptm_count[0] = 1;
                        else{
                            Nter_unk_ptm_count[i] = Nter_unk_ptm_count[i-1]+1;
                        }
                    }else{
                        if(i==0)    Nter_unk_ptm_count[0] = 0;
                        else    Nter_unk_ptm_count[i] = Nter_unk_ptm_count[i-1];
                    }
                    int ccom_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Cter_mass_ori[cuted_seq_length-1-i]);
                    if (ccom_index >=0 && ppmValuePtr->calculate_ppm_value(Cter_mass_ori[cuted_seq_length-1-i],mono_mass[ccom_index])) {
                        if(i==0)
                            Cter_ori_count[cuted_seq_length-1] = 1;
                        else{
                            Cter_ori_count[cuted_seq_length-1-i] = Cter_ori_count[cuted_seq_length-1-i+1]+1;
                        }
                    }else{
                        if(i==0)    Cter_ori_count[cuted_seq_length-1] = 0;
                        else     Cter_ori_count[cuted_seq_length-1-i] = Cter_ori_count[cuted_seq_length-1-i+1];
                    }
                    int cptm_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Cter_mass_ptm[cuted_seq_length-1-i]);
                    if (cptm_index >=0 && ppmValuePtr->calculate_ppm_value(Cter_mass_ptm[cuted_seq_length-1-i],mono_mass[cptm_index])) {
                        if(i==0)
                            Cter_ptm_count[cuted_seq_length-1] = 1;
                        else{
                            Cter_ptm_count[cuted_seq_length-1-i] = Cter_ptm_count[cuted_seq_length-1-i+1]+1;
                        }
                    }else{
                        if(i==0)    Cter_ptm_count[cuted_seq_length-1] = 0;
                        else     Cter_ptm_count[cuted_seq_length-1-i] = Cter_ptm_count[cuted_seq_length-1-i+1];
                    }
                    int cunk_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Cter_mass_unk[cuted_seq_length-1-i]);
                    if (cunk_index >=0 && ppmValuePtr->calculate_ppm_value(Cter_mass_unk[cuted_seq_length-1-i],mono_mass[cunk_index])) {
                        if(i==0)
                            Cter_unk_count[cuted_seq_length-1] = 1;
                        else{
                            Cter_unk_count[cuted_seq_length-1-i] = Cter_unk_count[cuted_seq_length-1-i+1]+1;
                        }
                    }else{
                        if(i==0)    Cter_unk_count[cuted_seq_length-1] = 0;
                        else     Cter_unk_count[cuted_seq_length-1-i] = Cter_unk_count[cuted_seq_length-1-i+1];
                    }
                    int cunk_ptm_index = mylib::speculate_lib::txpb_binarey_search_ex(mono_mass, mono_mass.size(), Cter_mass_unk_ptm[cuted_seq_length-1-i]);
                    if (cunk_ptm_index >=0 && ppmValuePtr->calculate_ppm_value(Cter_mass_unk_ptm[cuted_seq_length-1-i],mono_mass[cunk_ptm_index])) {
                        if(i==0)
                            Cter_unk_ptm_count[cuted_seq_length-1] = 1;
                        else{
                            Cter_unk_ptm_count[cuted_seq_length-1-i] = Cter_unk_ptm_count[cuted_seq_length-1-i+1]+1;
                        }
                    }else{
                        if(i==0)    Cter_unk_ptm_count[cuted_seq_length-1] = 0;
                        else     Cter_unk_ptm_count[cuted_seq_length-1-i] = Cter_unk_ptm_count[cuted_seq_length-1-i+1];
                    }
                }



                ///进行k均分假设
                int k = 20;
                int span = cuted_seq_length/k;
                int ptm_loc= 0,unk_loc = 0;
                int frg = 0,pre = 0,suf = 0,match_peak = 0;
                int cur_max_frg = 0,cur_max_pre = 0,cur_max_suf = 0,cur_max_match_peak = 0;
                while (ptm_loc <= cuted_seq_length && span != 0){
                    unk_loc = 0;
                    while(unk_loc <= cuted_seq_length){
                        ///判定部分
                        if(unk_loc == ptm_loc) {
                            if(unk_loc == cuted_seq_length-1)   break;
                            unk_loc += span;
                            if(unk_loc > cuted_seq_length-1)    unk_loc = cuted_seq_length-1;
                            continue;
                        }


                        ///计算匹配数
                        if(unk_loc < ptm_loc){
                            frg = Nter_unk_count[ptm_loc] - Nter_unk_count[unk_loc] + Cter_ptm_count[unk_loc] - Cter_ptm_count[ptm_loc];
                            pre = Nter_ori_count[unk_loc];
                            suf = Cter_ori_count[ptm_loc];
                            match_peak = Nter_ori_count[unk_loc] + Nter_unk_count[ptm_loc] - Nter_unk_count[unk_loc] + Nter_unk_ptm_count[Nter_unk_ptm_count.size()-1] -Nter_unk_ptm_count[ptm_loc]
                                         + Cter_ori_count[ptm_loc] + Cter_ptm_count[unk_loc] - Cter_ptm_count[ptm_loc] + Cter_unk_ptm_count[0] - Cter_unk_ptm_count[unk_loc];
                        }
                        else{
                            frg = Nter_ptm_count[unk_loc] - Nter_ptm_count[ptm_loc] + Cter_ptm_count[ptm_loc] - Cter_ptm_count[unk_loc];
                            pre = Nter_ori_count[ptm_loc];
                            suf = Cter_ori_count[unk_loc];
                            match_peak = Nter_ori_count[ptm_loc] + Nter_ptm_count[unk_loc] - Nter_ptm_count[ptm_loc] + Nter_unk_ptm_count[Nter_unk_ptm_count.size()-1] -Nter_unk_ptm_count[unk_loc]
                                         +Cter_ori_count[unk_loc] + Cter_unk_count[ptm_loc] - Cter_unk_count[unk_loc] + Cter_unk_ptm_count[0] - Cter_unk_ptm_count[ptm_loc];
                        }

                        if(frg > cur_max_frg){cur_max_frg = frg;}
                        if(pre > cur_max_pre){cur_max_pre = pre;}
                        if(suf > cur_max_suf){cur_max_suf = suf;}
                        if(match_peak > cur_max_match_peak){cur_max_match_peak = match_peak;}
                        ///结束条件
                        if(unk_loc == cuted_seq_length-1)   break;
                        unk_loc += span;
                        if(unk_loc > cuted_seq_length-1)    unk_loc = cuted_seq_length-1;
                    }
                    ///结束条件
                    if(ptm_loc == cuted_seq_length-1)
                        break;
                    ptm_loc += span;
                    if(ptm_loc > cuted_seq_length-1)
                        ptm_loc = cuted_seq_length-1;
                }
                if(cur_max_frg > ptm_candiTru_map[it_mp->first][0]->getMatchPeaks()) {
                    ptm_candiTru_map[it_mp->first][0]->setMatchPeaks(cur_max_frg);
                    ptm_candiTru_map[it_mp->first][0]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][0]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][0]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][0]->setCutC(seq_mass.size() - cutEnd  - 1);
                }else if(cur_max_frg == ptm_candiTru_map[it_mp->first][0]->getMatchPeaks() && unk_mass == 0){
                    ptm_candiTru_map[it_mp->first][0]->setMatchPeaks(cur_max_frg);
                    ptm_candiTru_map[it_mp->first][0]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][0]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][0]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][0]->setCutC(seq_mass.size() - cutEnd  - 1);
                }
                if(cur_max_pre > ptm_candiTru_map[it_mp->first][1]->getMatchPeaks()) {
                    ptm_candiTru_map[it_mp->first][1]->setMatchPeaks(cur_max_pre);
                    ptm_candiTru_map[it_mp->first][1]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][1]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][1]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][1]->setCutC(seq_mass.size() - cutEnd  - 1);
                }else if(cur_max_pre == ptm_candiTru_map[it_mp->first][1]->getMatchPeaks() && unk_mass == 0){
                    ptm_candiTru_map[it_mp->first][1]->setMatchPeaks(cur_max_pre);
                    ptm_candiTru_map[it_mp->first][1]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][1]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][1]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][1]->setCutC(seq_mass.size() - cutEnd  - 1);
                }
                if(cur_max_suf > ptm_candiTru_map[it_mp->first][2]->getMatchPeaks()) {
                    ptm_candiTru_map[it_mp->first][2]->setMatchPeaks(cur_max_suf);
                    ptm_candiTru_map[it_mp->first][2]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][2]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][2]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][2]->setCutC(seq_mass.size() - cutEnd  - 1);
                } else if(cur_max_suf == ptm_candiTru_map[it_mp->first][2]->getMatchPeaks() && unk_mass == 0){
                    ptm_candiTru_map[it_mp->first][2]->setMatchPeaks(cur_max_suf);
                    ptm_candiTru_map[it_mp->first][2]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][2]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][2]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][2]->setCutC(seq_mass.size() - cutEnd  - 1);
                }
                if(cur_max_match_peak > ptm_candiTru_map[it_mp->first][3]->getMatchPeaks()) {
                    ptm_candiTru_map[it_mp->first][3]->setMatchPeaks(cur_max_match_peak);
                    ptm_candiTru_map[it_mp->first][3]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][3]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][3]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][3]->setCutC(seq_mass.size() - cutEnd  - 1);
                }else if(cur_max_suf == ptm_candiTru_map[it_mp->first][2]->getMatchPeaks() && unk_mass == 0){
                    ptm_candiTru_map[it_mp->first][3]->setMatchPeaks(cur_max_match_peak);
                    ptm_candiTru_map[it_mp->first][3]->setUnkmass(unk_mass);
                    ptm_candiTru_map[it_mp->first][3]->setVarptmmass(it_mp->first);
                    ptm_candiTru_map[it_mp->first][3]->setCutN(cutStart);
                    ptm_candiTru_map[it_mp->first][3]->setCutC(seq_mass.size() - cutEnd  - 1);
                }
                Nter_mass_ori.clear();
                Nter_mass_ptm.clear();
                Nter_mass_unk.clear();
                Nter_mass_unk_ptm.clear();
                Cter_mass_ori.clear();
                Cter_mass_ptm.clear();
                Cter_mass_unk.clear();
                Cter_mass_unk_ptm.clear();
                Nter_ori_count.clear();
                Nter_ptm_count.clear();
                Nter_unk_count.clear();
                Nter_unk_ptm_count.clear();
                Cter_ori_count.clear();
                Cter_ptm_count.clear();
                Cter_unk_count.clear();
                Cter_unk_ptm_count.clear();
            }// for cutEnd
        }// for it_mp 遍历修饰组合
    }// for cutStart

}