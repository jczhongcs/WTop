//
// Created by wenzhong on 2023/3/22.
//

#include "prsm_unknown_ptms.hpp"
#include "../merge//wtop_prsm_processor.hpp"
#include "../merge//processor_management.h"
#include "candidate_prsm.hpp"


bool candidate_prsm_sort(candidate_prsm_ptr a,candidate_prsm_ptr b){
    return a->getPrsmBestScore()>=b->getPrsmBestScore();
}



//只有未知修饰的情况
void unknown_match_ptmZero(const vector<double> &mono_mass, const vector<double> &theo_mass,
                           vector<double> &deviSeq, double unknownPtmMass, int searchSize) {
    //获取匹配个数
    vector<double> theorySeq_UnknownPtm;
    vector<double> mono;
    mono.assign(mono_mass.begin() + searchSize, mono_mass.end());
    //添加峰值
    for (int i = 0; i < theo_mass.size(); ++i) {
        theorySeq_UnknownPtm.push_back(theo_mass[i] + unknownPtmMass);
    } // for end

    vector<double> short_1;
    vector<double> long_1;
    if (mono.size() > theorySeq_UnknownPtm.size()) {
        long_1 = mono;
        short_1 = theorySeq_UnknownPtm;
    } else {
        short_1 = mono;
        long_1 = theorySeq_UnknownPtm;
    }
    for (int i = 0; i < short_1.size(); i++) {
        int n = mylib::speculate_lib::txpb_binarey_search_ex(long_1, long_1.size(), short_1[i]);
        if (mylib::data_stream::_ppm_H(short_1[i], long_1[n])) {
            deviSeq.push_back(short_1[i] - long_1[n]);
        }
    }
    theorySeq_UnknownPtm.clear();
    mono.clear();

}

//未知修饰和已知修饰的情况
void unknown_match_ptmNum(const vector<double> &mono_mass,
                          const vector<double> &theo_mass,
                          vector<double> &deviSeq,
                          double ptmMass, double unknownPtmMass, int searchSize) {
    vector<double> theorySeq_UnknownPtm;
    vector<double> theorySeq_Ptm_UnknownPtm;
    vector<double> mono;
    mono.assign(mono_mass.begin() + searchSize, mono_mass.end());
    //添加峰值
    for (int i = 0; i < theo_mass.size(); ++i) {
        //            theorySeq_Ptm.push_back(theo_mass[i] + ptmMass) ;
        theorySeq_UnknownPtm.push_back(theo_mass[i] + unknownPtmMass);
        theorySeq_Ptm_UnknownPtm.push_back(theo_mass[i] + ptmMass + unknownPtmMass);
    } // for end

    vector<double> short_1;
    vector<double> long_1;
    if (mono_mass.size() > theorySeq_UnknownPtm.size()) {
        long_1 = mono;
        short_1 = theorySeq_UnknownPtm;

    } else {
        short_1 = mono;
        long_1 = theorySeq_UnknownPtm;
    }
    for (int i = 0; i < short_1.size(); i++) {
        int n = mylib::speculate_lib::txpb_binarey_search_ex(long_1, long_1.size(), short_1[i]);
        if (mylib::data_stream::_ppm_H(short_1[i], long_1[n])) {
            deviSeq.push_back(short_1[i] - long_1[n]);
        }
    }

    if (mono_mass.size() > theorySeq_Ptm_UnknownPtm.size()) {
        long_1 = mono;
        short_1 = theorySeq_Ptm_UnknownPtm;

    } else {
        short_1 = mono;
        long_1 = theorySeq_Ptm_UnknownPtm;
    }

    for (int i = 0; i < short_1.size(); i++) {
        int n = mylib::speculate_lib::txpb_binarey_search_ex(long_1, long_1.size(), short_1[i]);
        if (mylib::data_stream::_ppm_H(short_1[i], long_1[n])) {
            deviSeq.push_back(short_1[i] - long_1[n]);
        }
    }
    theorySeq_UnknownPtm.clear();
    theorySeq_Ptm_UnknownPtm.clear();
}

int judge_best_unknown_mass(msalign_ptr msalignPtr,
                            prsm::temp_prsm_argument_ptr tempPrsmArgumentPtr,
                            double ptmMass, double unknown_mass, int searchN, int searchC) {
    vector<double> n_devi_seq;
    vector<double> c_devi_seq;
    if (ptmMass == 0) {
        unknown_match_ptmZero(msalignPtr->getIonsMassContainer(),
                              tempPrsmArgumentPtr->theory_ions_mass_c,
                              c_devi_seq, unknown_mass, searchC);
        unknown_match_ptmZero(msalignPtr->getIonsMassContainer(),
                              tempPrsmArgumentPtr->theory_ions_mass_n,
                              n_devi_seq, unknown_mass, searchN);
    } else {
        unknown_match_ptmNum(msalignPtr->getIonsMassContainer(),
                             tempPrsmArgumentPtr->theory_ions_mass_c,
                             c_devi_seq, ptmMass, unknown_mass, searchC);
        unknown_match_ptmNum(msalignPtr->getIonsMassContainer(),
                             tempPrsmArgumentPtr->theory_ions_mass_n,
                             n_devi_seq, ptmMass, unknown_mass, searchN);
    }
    return n_devi_seq.size() + c_devi_seq.size();
}

int get_unknown_mass_searchSize(const vector<double> &mono_mass,
                                const vector<double> &theo_mass) {
    int count = 0;
    //获取匹配个数
    vector<double> short_1;
    vector<double> long_1;
    char s_log;
    if (mono_mass.size() > theo_mass.size()) {
        long_1 = mono_mass;
        short_1 = theo_mass;
        s_log = 't';
    } else {
        short_1 = mono_mass;
        long_1 = theo_mass;
        s_log = 'm';
    }
    int n = 0;
    vector<int> mIndex;
    for (int i = 0; i < short_1.size(); i++) {
        n = mylib::speculate_lib::txpb_binarey_search_ex(long_1, long_1.size(), short_1[i]);
        if (mylib::data_stream::_ppm_H(short_1[i], long_1[n])) {
            count++;
            if (s_log == 'm') {
                mIndex.push_back(i);
            } else {
                mIndex.push_back(n);
            }
        }
    }

    int searchSize = mIndex.size();
    if (searchSize < 5) {
        return 0;
    }
    return mIndex[searchSize-1];
}

//2023-1-8 final merge
double determine_unknown_mass(double ptm_mass, double precursor_mass, double unk_mass,
                              prsm::temp_prsm_argument_ptr tempPrsmArgumentPtr,
                              msalign_ptr msalignPtr,
                              protein_ptr proteinPtr,
                              prsm::config_argument_ptr prsmArg) {
    vector<double> unknownMassSeq;
    // 前提质量判断的时候
    double flog = prsmArg->unknownPTMsMassFlog;
    for (int i = 0; i < prsmArg->unknownPTMsMassFlogFloatSize; i++) {
        unknownMassSeq.push_back(unk_mass + flog);
        flog += prsmArg->unknownPTMsMassFlogChangeMass;
    }

    int searchSize_C = 0, searchSize_N = 0;
//    if (prsmArg->BYIonsActions.count(msalignPtr->getActivationSub())) {
//        searchSize_C = get_unknown_mass_searchSize(msalignPtr->getIonsMassContainer(), proteinPtr->getTheoryMassCBy());
//        searchSize_N = get_unknown_mass_searchSize(msalignPtr->getIonsMassContainer(), proteinPtr->getTheoryMassNBy());
//    }
//    if (prsmArg->CZIonsActions.count(msalignPtr->getActivationSub())) {
//        searchSize_C = get_unknown_mass_searchSize(msalignPtr->getIonsMassContainer(), proteinPtr->getTheoryMassCCz());
//        searchSize_N = get_unknown_mass_searchSize(msalignPtr->getIonsMassContainer(), proteinPtr->getTheoryMassNCz());
//    }
    searchSize_C = 0;
    searchSize_N = 0;





    //建立理论峰值,判断最佳合适
    double bestCount = 0;
    double bestUnknown = unk_mass; //真正的修饰
    ppm_value_ptr ppmValuePtr = std::make_shared<ppm_value>(prsmArg->min_ppm_value);
    for (double unknownMass : unknownMassSeq) { //确定真正的未知修饰
        double ppm = 0;
        bool boolean = false;
        if (prsmArg->BYIonsActions.count(msalignPtr->getActivationSub())) {
            boolean = ppmValuePtr->calculate_ppm_with_H(unknownMass + ptm_mass + tempPrsmArgumentPtr->theory_ions_mass_c.back(),
                                                       precursor_mass);
        }
        if (prsmArg->CZIonsActions.count(msalignPtr->getActivationSub())) {
            boolean = ppmValuePtr->calculate_ppm_with_H(unknownMass + ptm_mass + tempPrsmArgumentPtr->theory_ions_mass_c.back() + prsmArg->sub_etd, precursor_mass);
        }
        if (boolean) {
            int count = judge_best_unknown_mass(msalignPtr, tempPrsmArgumentPtr,
                                                ptm_mass, unknownMass,
                                                searchSize_N, searchSize_C);
            if (count > bestCount) {
                bestCount = count;
                bestUnknown = unknownMass;
            }

        }
    }

//    delete ppmValuePtr;

    if (fabs(bestUnknown) < 0.1) {
        return 0;
    }
    return bestUnknown;
}

void prsm_unknown_ptms::terminal_trunc_mass(vector<double> &theoryMass, vector<double> &theoryMassC,
                                            protein_ptr pro, int cutSizeN, int cutSizeC) {
    vector<double> seqMass;
    pro->get_seq_mass(seqMass);
    double mass = 0;
    for (int i = 0; i < cutSizeN; ++i) {
        mass += seqMass[i];
    }
    for (int i = 0; i < theoryMass.size(); ++i) {
        theoryMass[i] -= mass;
    }
    double mass2 = 0;
    int t = seqMass.size() - 1;
    for (int i = 0; i < cutSizeC; ++i) {
        mass2 += seqMass[t - i];
    }
    for (int i = 0; i < theoryMassC.size(); ++i) {
        theoryMassC[i] -= mass2;
    }
}

void init_unknown_ptmTable(prsm::temp_prsm_argument_ptr arg,
                           vector<double> &modifierTable, map<double, string> &mmap,
                           double ptmMass, string ptm, double bestUnknown) {
    modifierTable.clear();
    mmap.clear();
    mmap.insert(make_pair(0, "0"));
    if (fabs(ptmMass) <= 1e-15) {
        //只有未知修饰
        mmap.insert(make_pair(bestUnknown, "*"));
        arg->var_mod_mass = bestUnknown;
    } else {
        //插入已知修饰，未知修饰，未知+已知修饰组合
        mmap.insert(make_pair(ptmMass, ptm));
        mmap.insert(make_pair(bestUnknown, "*"));
        mmap.insert(make_pair(bestUnknown + ptmMass, "*" + ptm + ""));
        arg->var_mod_mass = ptmMass + bestUnknown;
    }
    for (map<double, string>::iterator it = mmap.begin(); it != mmap.end(); ++it) {
        modifierTable.push_back(it->first);
    }
}

void prsm_unknown_ptms::storagePrSMsMsg(protein_ptr proteinPtr, msalign_ptr msalignPtr,
                                        prsm::temp_prsm_argument_ptr &tempPrsmArgumentPtr) {
    msalignPtr->setPrsmBestScore(tempPrsmArgumentPtr->score);
    msalignPtr->setPrsmMatchPeaks(tempPrsmArgumentPtr->peaks_match);
    msalignPtr->setPrsmCutN(tempPrsmArgumentPtr->cut_location_n);
    msalignPtr->setPrsmCutC(tempPrsmArgumentPtr->cut_location_c);
    msalignPtr->setPrsmPtm(tempPrsmArgumentPtr->ptms);
    msalignPtr->setPrsmMatchProteinName(proteinPtr->getProteinName());
    msalignPtr->setPrsmMatchProteinSeq(proteinPtr->getProteinSequence());
    msalignPtr->setPrsmMutX(tempPrsmArgumentPtr->mut_x);
    msalignPtr->setPrsmIonsPairNum(tempPrsmArgumentPtr->pairIonsNub);
    msalignPtr->setPrsmComplementsIonsNum(tempPrsmArgumentPtr->complement_ions_number);
    msalignPtr->setPrsmSeqScore(tempPrsmArgumentPtr->prsm_score);
    msalignPtr->setPrsmSeqScoreCount(tempPrsmArgumentPtr->seq_score_count);
    msalignPtr->setPrsmMatchFragmentIonsNum(tempPrsmArgumentPtr->matchFragmentIonsSize);
    msalignPtr->c_ter = tempPrsmArgumentPtr->ions_c;
    msalignPtr->n_ter = tempPrsmArgumentPtr->ions_n;

    ///调试打分
    msalignPtr->setFeature1(tempPrsmArgumentPtr->feature1);
    msalignPtr->setFeature2(tempPrsmArgumentPtr->feature2);
    msalignPtr->setFeature3(tempPrsmArgumentPtr->feature3);
    msalignPtr->setFeature4(tempPrsmArgumentPtr->feature4);
    msalignPtr->setFeature5(tempPrsmArgumentPtr->feature5);
    msalignPtr->setFeature6(tempPrsmArgumentPtr->feature6);
    msalignPtr->setFeature7(tempPrsmArgumentPtr->feature7);
    msalignPtr->setFeature8(tempPrsmArgumentPtr->feature8);
    msalignPtr->setFeature9(tempPrsmArgumentPtr->feature9);
    msalignPtr->setFeature10(tempPrsmArgumentPtr->feature10);
    msalignPtr->setFeature11(tempPrsmArgumentPtr->feature11);
    msalignPtr->setFeature12(tempPrsmArgumentPtr->feature12);
    msalignPtr->setFeature13(tempPrsmArgumentPtr->feature13);
    msalignPtr->setFeature14(tempPrsmArgumentPtr->feature14);
    msalignPtr->setFeature15(tempPrsmArgumentPtr->feature15);
    msalignPtr->setFeature16(tempPrsmArgumentPtr->feature16);
    msalignPtr->setFeature17(tempPrsmArgumentPtr->feature17);
    msalignPtr->setFeature18(tempPrsmArgumentPtr->feature18);
    msalignPtr->setFeature19(tempPrsmArgumentPtr->feature19);
    msalignPtr->setFeature20(tempPrsmArgumentPtr->feature20);
    msalignPtr->setFeature21(tempPrsmArgumentPtr->feature21);
    msalignPtr->setFeature22(tempPrsmArgumentPtr->feature22);
    msalignPtr->setFeature23(tempPrsmArgumentPtr->feature23);
    msalignPtr->setFeature24(tempPrsmArgumentPtr->feature24);
    msalignPtr->setFeature25(tempPrsmArgumentPtr->feature25);
    msalignPtr->setFeature26(tempPrsmArgumentPtr->feature26);
    msalignPtr->setFeature27(tempPrsmArgumentPtr->feature27);
    msalignPtr->setFeature28(tempPrsmArgumentPtr->feature28);
    msalignPtr->setFeature29(tempPrsmArgumentPtr->feature29);
    msalignPtr->setFeature30(tempPrsmArgumentPtr->feature30);
    msalignPtr->setFeature31(tempPrsmArgumentPtr->feature31);
    msalignPtr->setFeature32(tempPrsmArgumentPtr->feature32);
    msalignPtr->setFeature33(tempPrsmArgumentPtr->feature33);
    msalignPtr->setFeature34(tempPrsmArgumentPtr->feature34);
    msalignPtr->setFeature35(tempPrsmArgumentPtr->feature35);
    msalignPtr->setFeature36(tempPrsmArgumentPtr->feature36);
    msalignPtr->setFeature37(tempPrsmArgumentPtr->feature37);
    msalignPtr->setFeature38(tempPrsmArgumentPtr->feature38);
    msalignPtr->setFeature39(tempPrsmArgumentPtr->feature39);
    msalignPtr->setFeature40(tempPrsmArgumentPtr->feature40);
    msalignPtr->setFeature41(tempPrsmArgumentPtr->feature41);
    msalignPtr->setFeature42(tempPrsmArgumentPtr->feature42);
    msalignPtr->setFeature43(tempPrsmArgumentPtr->feature43);
    msalignPtr->setFeature44(tempPrsmArgumentPtr->feature44);
    msalignPtr->setFeature45(tempPrsmArgumentPtr->feature45);
    msalignPtr->setFeature46(tempPrsmArgumentPtr->feature46);
    msalignPtr->setFeature47(tempPrsmArgumentPtr->feature47);
    msalignPtr->setFeature48(tempPrsmArgumentPtr->feature48);
    msalignPtr->setFeature49(tempPrsmArgumentPtr->feature49);
    msalignPtr->setFeature50(tempPrsmArgumentPtr->feature50);

    //test set byions 23.7.12




    //    spact->prsm_best_score = arg->score;
    //    spact->best_peak = arg->peaks_match;
    //    spact->cut_c_n = arg->cut_location_n;
    //    spact->cut_n_c = arg->cut_location_c;
    //    spact->best_ptm = arg->ptms;
    //    spact->match_protein_name = pro->protein_name;
    //    spact->match_protein_seq = pro->protein_sequence;
    //    spact->mut_s = arg->mut_x;
    //    spact->bestIonsPairNub = arg->pairIonsNub;
    //    spact->huBuCount = arg->complement_ions_number;
    //    spact->seqScore = arg->prsm_score;
    //    spact->seqScoreCount = arg->seq_score_count;
    //    spact->matchFragmentIonsSize = arg->matchFragmentIonsSize;
}

void prsm_unknown_ptms::search_ptms_(prsm::modify_ptr one_mod_ptr,
                                     vector<msalign_ptr> &msalign_ptrs,
                                     prsm::config_argument_ptr configArgumentPtr,
                                     const string &variablePTMsFileOutPath,
                                     const string &outFIleName) {
    for (vector<msalign_ptr>::iterator ms_it = msalign_ptrs.begin(); ms_it != msalign_ptrs.end(); ++ms_it, ++cur_ms_size) {
        if ((*ms_it)->getPrecursorMass() < configArgumentPtr->min_precursor_mass) {
            continue;
        }
        if (cur_ms_size % 10 == 0) {
            cout << "\rcurrent thread id : " << std::this_thread::get_id() << " . "
                 << "Number of PrSMs : " << cur_ms_size << " th ."
                 << "Progress bar : " << (double)cur_ms_size / (double)msalign_ptrs.size() * 100 << " % ";
            std::cout.flush();
        }
       ///原鉴定，只输出一条质谱对应一条蛋白变体
        vector<protein_ptr>::const_iterator pr_it = (*ms_it)->getCandidateProteins().begin();   //候选蛋白库起始位置
        vector<protein_ptr>::const_iterator pr_it_end = (*ms_it)->getCandidateProteins().end(); //候选蛋白库终止位置

        cout<<endl<<"unk sp:"<<(*ms_it)->getId()<<"-candi_protein_size:"<<(*ms_it)->getCandidateProteins().size()<<endl;

        for (; pr_it != pr_it_end; ++pr_it) {
//            cout<<endl<<"sp:"<<(*ms_it)->getId()<<"  pro:"<<(*pr_it)->getProteinTitle()<<endl;
            prsm_unknown_ptms::two_unknown_ptm_Process(one_mod_ptr, (*pr_it), (*ms_it),
                                                       configArgumentPtr, variablePTMsFileOutPath);
        }
//        cout<<endl<<"sp:"<<(*ms_it)->getId()<<" prsm end"<<endl;

        prsm::modify_ptr resultMPtr = (*ms_it)->getPrsmUnknownModPtr();
        ///输出PrSM结果
//        prsmMsgWritePtr->writeResultNews((*ms_it), outFIleName, resultMPtr);


        ///原鉴定结束.记得打开最后的delete resultMptr



        ///test_out_candidate_prsm
        bool out_candidate_prsm = false;
        if(out_candidate_prsm){
            //sort((*ms_it)->candi_prsms.begin(),(*ms_it)->candi_prsms.end(),candidate_prsm_sort);
            ofstream candidate_out_file;
            candidate_out_file.open("D:\\LYC\\wtop\\WTop\\bin\\candidate_prsm.txt");
            for(auto cp:(*ms_it)->candi_prsms){
//                candidate_out_file<<"###begin###\n";
//                candidate_out_file<<"prsm_best_score:"<<cp->getPrsmBestScore()<<endl;
//                candidate_out_file<<"prsm_match_peaks:"<<cp->getPrsmMatchPeaks()<<endl;
//                candidate_out_file<<"prsm_spectrum_id:"<<(*ms_it)->getId()<<endl;
//                candidate_out_file<<"prsm_match_protein_name:"<<cp->getPrsmMatchProteinName()<<endl;
//                candidate_out_file<<"prsm_cut_n:"<<cp->getPrsmCutN()<<endl;
//                candidate_out_file<<"prsm_cut_c"<<cp->getPrsmCutC()<<endl;
//                candidate_out_file<<"###end###\n\n";
///cout
                cout<<"###begin###\n";
                cout<<"prsm_best_score:"<<cp->getPrsmBestScore()<<endl;
                cout<<"prsm_match_peaks:"<<cp->getPrsmMatchPeaks()<<endl;
                cout<<"prsm_spectrum_id:"<<(*ms_it)->getId()<<endl;
                cout<<"prsm_match_protein_name:"<<cp->getPrsmMatchProteinName()<<endl;
                cout<<"prsm_cut_n:"<<cp->getPrsmCutN()<<endl;
                cout<<"prsm_cut_c"<<cp->getPrsmCutC()<<endl;
                cout<<"###end###\n\n";
            }
            candidate_out_file.close();
        }



//        delete (*ms_it);
//        delete resultMPtr;
    }

    ///合并Wtop_prsm.test_combine,并清除临时文件


}

void prsm_unknown_ptms::two_unknown_ptm_Process(prsm::modify_ptr mp,
                                                protein_ptr proteinPtr,
                                                msalign_ptr msalignPtr,
                                                prsm::config_argument_ptr prsmArg,
                                                const string &variablePTMsFileOutPath) {

    msalignPtrVec evalueMsVec;

    //先给理论峰值,添加所有修饰峰值
    //最后定位蛋白质形
    //减前端截断，选取最为合适的蛋白质形分析
    //此处cut_n表示n端截断的截止索引，即[0,cut_n]截去，若cut_n == -1，则无n端截断
    int cut_n = -1, cut_c = 0;
    prsm_terminal_truncation_ptr prsmTerminalTruncationPtr = new prsm_terminal_truncation();
    bool judegeCter = true;


    map<double,candi_trucation_vec> ptm_candiTru_map;
    if(proteinPtr->getProteinSequence().size() <= 1000 && judegeCter) {
        for(auto& it_mp : mp->mapi) {
            candi_trucation_vec newVec;
            for(int i=0; i<4; i++) {
                newVec.push_back(make_shared<candi_trucation>(-1,0,0.0,0.0,0));
            }
            ptm_candiTru_map.emplace(it_mp.first, std::move(newVec));
        }
        ///测试推测NC两端截
        if (prsmArg->BYIonsActions.count(msalignPtr->getActivationSub())) {prsmTerminalTruncationPtr->locate_truncation_by(mp, proteinPtr, msalignPtr,prsmArg,ptm_candiTru_map);}
        if (prsmArg->CZIonsActions.count(msalignPtr->getActivationSub())) {prsmTerminalTruncationPtr->locate_truncation_cz(mp, proteinPtr, msalignPtr,prsmArg,ptm_candiTru_map);}

    }
    else{
        ///原只推断N端截断
        if (prsmArg->BYIonsActions.count(msalignPtr->getActivationSub())) {prsmTerminalTruncationPtr->judge_terminal_truncation_By(mp, proteinPtr, msalignPtr,prsmArg, cut_n, cut_c);}
        if (prsmArg->CZIonsActions.count(msalignPtr->getActivationSub())) {prsmTerminalTruncationPtr->judge_terminal_truncation_Cz(mp, proteinPtr, msalignPtr,prsmArg, cut_n, cut_c);}
        if (cut_n == -1 && cut_c == 0) { return; }
        ///原截断推断结束
        for(auto& it_mp : mp->mapi) {
            candi_trucation_vec newVec;
            for(int i=0; i<1; i++) {
                newVec.push_back(make_shared<candi_trucation>(cut_n,cut_c,0.0,0.0,0));
            }
            ptm_candiTru_map.emplace(it_mp.first, std::move(newVec));
        }
    }

    delete prsmTerminalTruncationPtr;


    //遍历每一个修饰
    for (map<double, string>::iterator it_mp = mp->mapi.begin(); it_mp != mp->mapi.end(); ++it_mp) {
        candi_trucation_vec cur_var_tru_vec = ptm_candiTru_map[it_mp->first];
        for(auto cur_var_tru:cur_var_tru_vec){
            prsm::temp_prsm_argument_ptr tempPrsmArgPtr = new prsm::temp_prsm_argument();
            if(proteinPtr->getProteinSequence().size() <= 1000 && judegeCter){
                tempPrsmArgPtr->cut_location_n = cur_var_tru->getCutN();
                tempPrsmArgPtr->cut_location_c = cur_var_tru->getCutC();
                cut_n = cur_var_tru->getCutN();
                cut_c = cur_var_tru->getCutC();
                if (cut_n == -1 && cut_c == 0)  ///匹配峰数小于5直接跳过
                    continue;
            }
            else{
                tempPrsmArgPtr->cut_location_n = cut_n;
                tempPrsmArgPtr->cut_location_c = cut_c;
            }
            // 初始化理论峰序列
            if (prsmArg->BYIonsActions.count(msalignPtr->getActivationSub())) {
                tempPrsmArgPtr->theory_ions_mass_n.assign(proteinPtr->getTheoryMassNBy().begin() + cut_n,
                                                          proteinPtr->getTheoryMassNBy().end() - cut_c);
                tempPrsmArgPtr->theory_ions_mass_c.assign(proteinPtr->getTheoryMassCBy().begin() + cut_c,
                                                          proteinPtr->getTheoryMassCBy().end() - cut_n);
            }
            if (prsmArg->CZIonsActions.count(msalignPtr->getActivationSub())) {
                tempPrsmArgPtr->theory_ions_mass_n.assign(proteinPtr->getTheoryMassNCz().begin() + cut_n,
                                                          proteinPtr->getTheoryMassNCz().end() - cut_c);
                tempPrsmArgPtr->theory_ions_mass_c.assign(proteinPtr->getTheoryMassCCz().begin() + cut_c,
                                                          proteinPtr->getTheoryMassCCz().end() - cut_n);
            }
            if(tempPrsmArgPtr->theory_ions_mass_n.size()<=5 || tempPrsmArgPtr->theory_ions_mass_c.size()<=5){
                delete tempPrsmArgPtr;
                continue;
            }

            // 减去截断质量
            terminal_trunc_mass(tempPrsmArgPtr->theory_ions_mass_n,tempPrsmArgPtr->theory_ions_mass_c,proteinPtr, cut_n, cut_c);

            double ms_precursor_mass = msalignPtr->getPrecursorMass();

            double unk_mass;
            if(proteinPtr->getProteinSequence().size() <= 1000 && judegeCter){
                unk_mass = cur_var_tru->getUnkmass();
            }
            else{
                if (prsmArg->BYIonsActions.count(msalignPtr->getActivationSub())) {
                    unk_mass = ms_precursor_mass - tempPrsmArgPtr->theory_ions_mass_c.back() - it_mp->first;
                }
                if (prsmArg->CZIonsActions.count(msalignPtr->getActivationSub())) {
                    unk_mass = ms_precursor_mass - tempPrsmArgPtr->theory_ions_mass_c.back() - it_mp->first - prsmArg->sub_etd;
                }
            }
            double bestUnknown = determine_unknown_mass(it_mp->first, ms_precursor_mass, unk_mass,
                                                        tempPrsmArgPtr, msalignPtr, proteinPtr, prsmArg);
            if (fabs(bestUnknown) > prsmArg->unknownPTMsMassMax + 2) {
                delete tempPrsmArgPtr;
                continue;
            }

            vector<double> modifierTable; //获取当前修饰所有质量
            map<double, string> mmap;     //double - 组合质量 , string - 组合修饰
            // 初始化mmap 和 modifierTable
            ///modifierTable：vector<修饰质量>，mmap：<修饰质量，修饰标记符号>
            init_unknown_ptmTable(tempPrsmArgPtr, modifierTable, mmap,
                                  it_mp->first, it_mp->second, bestUnknown);
            prsm::modify_ptr temp_mod_ptr = std::make_shared<prsm::Modify>(variablePTMsFileOutPath); //这里参数无影响
            temp_mod_ptr->Init_modify();
            temp_mod_ptr->mapi = mmap;
            temp_mod_ptr->massStringToStringMass();
            sort(modifierTable.begin(), modifierTable.end());

            map<string, string> ssMap;
            ssMap["PTM_Mass"] = to_string(it_mp->first);
            ssMap["Unknown_PTM_Mass"] = to_string(bestUnknown);

            // PrSMs 鉴定流程
            prsm_unknown_ptms::prsm_unk_process(temp_mod_ptr, proteinPtr, msalignPtr,
                                                tempPrsmArgPtr, modifierTable, prsmArg,it_mp->first,bestUnknown);

            ///用于存储各种修饰的Top1，共mapi.size()个，用于Evalue候选。同时需打开下方读写xml代码块
            bool out_for_evalue = true;
            if(out_for_evalue){
                msalign_ptr curModMsalign = std::make_shared<msalign>(*msalignPtr);
                prsm_unknown_ptms::storagePrSMsMsg(proteinPtr, curModMsalign, tempPrsmArgPtr);
                curModMsalign->setPrsmPtmsMass(bestUnknown + it_mp->first);
                curModMsalign->setPrsmUnknownPtmMassMap(ssMap);
                prsm::modify_ptr tmp = std::make_shared<prsm::Modify>(*temp_mod_ptr);
                curModMsalign->setPrsmUnknownModPtr(tmp);
                if (prsmArg->BYIonsActions.count(curModMsalign->getActivationSub())) {
                    curModMsalign->setPrsmTheoryProteoformMass(curModMsalign->getPrsmPtmsMass() + tempPrsmArgPtr->theory_ions_mass_c.back());
                }
                if (prsmArg->CZIonsActions.count(curModMsalign->getActivationSub())) {
                    curModMsalign->setPrsmTheoryProteoformMass(curModMsalign->getPrsmPtmsMass() + tempPrsmArgPtr->theory_ions_mass_c.back() + prsmArg->sub_etd);
                }
                evalueMsVec.push_back(std::make_shared<msalign>(*curModMsalign));
//            itmpCount++;
//                delete curModMsalign;
            }
            ///end test


            if ((tempPrsmArgPtr->score - msalignPtr->getPrsmBestScore()) > 0.000001) {
                // 交换公共部分的最佳PrSM result
                prsm_unknown_ptms::storagePrSMsMsg(proteinPtr, msalignPtr, tempPrsmArgPtr);
                msalignPtr->setPrsmPtmsMass(bestUnknown + it_mp->first);
                msalignPtr->setPrsmUnknownPtmMassMap(ssMap);
//                delete msalignPtr->getPrsmUnknownModPtr();
                msalignPtr->setPrsmUnknownModPtr(temp_mod_ptr);
                if (prsmArg->BYIonsActions.count(msalignPtr->getActivationSub())) {
                    msalignPtr->setPrsmTheoryProteoformMass(msalignPtr->getPrsmPtmsMass() + tempPrsmArgPtr->theory_ions_mass_c.back());
                }
                if (prsmArg->CZIonsActions.count(msalignPtr->getActivationSub())) {
                    msalignPtr->setPrsmTheoryProteoformMass(msalignPtr->getPrsmPtmsMass() + tempPrsmArgPtr->theory_ions_mass_c.back() + prsmArg->sub_etd);
                }
            }
//            else {
//                delete temp_mod_ptr;
//            }
            delete tempPrsmArgPtr;

        }
        cur_var_tru_vec.clear();
    }//for it_mp



    ptm_candiTru_map.clear();
    ///与out_for_evalue变量同步
    if(evalueMsVec.size() != 0){
        string xmlFileName = "wtop_prsm.test_combined";
        prsm_msg_write::OutXmlToEvalue(evalueMsVec, xmlFileName);
//        for(auto x:evalueMsVec){
//            delete x->getPrsmUnknownModPtr();
//            delete x;
//        }
        evalueMsVec.clear();
    }
}

void prsm_unknown_ptms::prsm_unk_process(prsm::modify_ptr modptr,
                                     protein_ptr proteinPtr,
                                     msalign_ptr msalignPtr,
                                     prsm::temp_prsm_argument_ptr tempPrsmArgumentPtr,
                                     vector<double> &ptms_mass,
                                     prsm::config_argument_ptr prsmArg,double var_ptm_mass,double unk_mass) {


    //Coarse alignment
    //    cout << " cwdtw - start " << endl ;
    prsmCwdtwPtr->cwdtw_algorithm(tempPrsmArgumentPtr->theory_ions_mass_n, msalignPtr->getIonsMassContainer(),
                                  tempPrsmArgumentPtr->alignment_ions_n, tempPrsmArgumentPtr->ref_score_n,
                                  tempPrsmArgumentPtr->peer_score_n);
    prsmCwdtwPtr->cwdtw_algorithm(tempPrsmArgumentPtr->theory_ions_mass_c, msalignPtr->getIonsMassContainer(),
                                  tempPrsmArgumentPtr->alignment_ions_c, tempPrsmArgumentPtr->ref_score_c,
                                  tempPrsmArgumentPtr->peer_score_c);
    //    cout << " cwdtw - end " << endl ;

    //    cout << " scopeAlignment - start " << endl ;
    //tempPrsmArgumentPtr->var_mod_mass为已知修饰+未知修饰质量和
    prsmScopeAlignmentPtr->scopeAlignment(tempPrsmArgumentPtr->theory_ions_mass_n,
                                          msalignPtr->getIonsMassContainer(),
                                          tempPrsmArgumentPtr->alignment_ions_n,
                                          ptms_mass, tempPrsmArgumentPtr->var_mod_mass, prsmArg->ptmMassMax);
    prsmScopeAlignmentPtr->scopeAlignment(tempPrsmArgumentPtr->theory_ions_mass_c,
                                          msalignPtr->getIonsMassContainer(),
                                          tempPrsmArgumentPtr->alignment_ions_c,
                                          ptms_mass, tempPrsmArgumentPtr->var_mod_mass, prsmArg->ptmMassMax);
    //    cout << " scopeAlignment - end " << endl ;

    prsmAlignmentFilterPtr->get_mass_sub(tempPrsmArgumentPtr->theory_ions_mass_n, tempPrsmArgumentPtr->theory_ions_mass_c,
                                         msalignPtr->getIonsMassContainer(),
                                         tempPrsmArgumentPtr->sub_n, tempPrsmArgumentPtr->sub_c,
                                         tempPrsmArgumentPtr->alignment_ions_n, tempPrsmArgumentPtr->alignment_ions_c);

    prsmAlignmentFilterPtr->index_ions_mass(tempPrsmArgumentPtr->sub_n, tempPrsmArgumentPtr->filter_ions_n,
                                            ptms_mass, tempPrsmArgumentPtr->var_mod_mass);
    prsmAlignmentFilterPtr->index_ions_mass(tempPrsmArgumentPtr->sub_c, tempPrsmArgumentPtr->filter_ions_c,
                                            ptms_mass, tempPrsmArgumentPtr->var_mod_mass);

    prsmAlignmentFilterPtr->get_ppm_unknown(tempPrsmArgumentPtr->theory_ions_mass_n, msalignPtr->getIonsMassContainer(),
                                            tempPrsmArgumentPtr->alignment_ions_n, tempPrsmArgumentPtr->filter_ions_n,
                                            ptms_mass, tempPrsmArgumentPtr->ppm_value_n, tempPrsmArgumentPtr->var_mod_mass);
    prsmAlignmentFilterPtr->get_ppm_unknown(tempPrsmArgumentPtr->theory_ions_mass_c, msalignPtr->getIonsMassContainer(),
                                            tempPrsmArgumentPtr->alignment_ions_c, tempPrsmArgumentPtr->filter_ions_c,
                                            ptms_mass, tempPrsmArgumentPtr->ppm_value_c, tempPrsmArgumentPtr->var_mod_mass);

    prsmAlignmentFilterPtr->get_ions_peaks(tempPrsmArgumentPtr->merge_ions_n, tempPrsmArgumentPtr->alignment_ions_n,
                                           tempPrsmArgumentPtr->filter_ions_n, tempPrsmArgumentPtr->ppm_value_n,
                                           prsmArg->min_ppm_value);
    prsmAlignmentFilterPtr->get_ions_peaks(tempPrsmArgumentPtr->merge_ions_c, tempPrsmArgumentPtr->alignment_ions_c,
                                           tempPrsmArgumentPtr->filter_ions_c, tempPrsmArgumentPtr->ppm_value_c,
                                           prsmArg->min_ppm_value);

    vector<node> N_ions_c_n, N_ions_n_c;
    vector<modi> N_Ptms;
    vector<node> C_ions_c_n, C_ions_n_c;
    vector<modi> C_Ptms;
    prsmIonsLocationPtr->ions_location_match_N(N_ions_c_n, N_ions_n_c, tempPrsmArgumentPtr->merge_ions_n,
                                               tempPrsmArgumentPtr->merge_ions_c, N_Ptms, ptms_mass,
                                               tempPrsmArgumentPtr->get_lenth(tempPrsmArgumentPtr->theory_ions_mass_n.size()),
                                               tempPrsmArgumentPtr->var_mod_mass);

    prsmIonsLocationPtr->ions_location_match_C(C_ions_c_n, C_ions_n_c, tempPrsmArgumentPtr->merge_ions_n,
                                               tempPrsmArgumentPtr->merge_ions_c,C_Ptms, ptms_mass,
                                               tempPrsmArgumentPtr->get_lenth(tempPrsmArgumentPtr->theory_ions_mass_n.size()),
                                               tempPrsmArgumentPtr->var_mod_mass);



    prsmIonsLocationPtr->confirm_final_ions(N_ions_c_n, N_ions_n_c, N_Ptms,
                                            C_ions_c_n, C_ions_n_c, C_Ptms,
                                            tempPrsmArgumentPtr);
    if (tempPrsmArgumentPtr->ions_n.size() == 0 && tempPrsmArgumentPtr->ions_c.size() == 0) {
        return;
    }

    /// test by zhou
//    cout << " ============ n ============ " << endl ;
//    for (node n : tempPrsmArgumentPtr->ions_n) {
//            cout <<"n.thoe_id:"<< n.thoe_id << " "
//                 <<"n.mono_id:"<< n.mono_id << " "<<endl
//                 <<"theory_ions_mass:"<< tempPrsmArgumentPtr->theory_ions_mass_n[n.thoe_id] << " "
//                 <<"spec_Ions_Mass:"<<setprecision(10)<< msalignPtr->getIonsMassContainer()[n.mono_id] << " "
//                 <<"index:"<< n.index << " "
//                 <<"ppm:"<< n.ppm << endl ;
//        }
//        cout << " ============ c ============ " << endl ;
//        for (node n : tempPrsmArgumentPtr->ions_c) {
//            cout <<"n.thoe_id:"<< n.thoe_id << " "
//                 <<"n.mono_id:"<< n.mono_id << " "<<endl
//                 <<"theory_ions_mass:"<< tempPrsmArgumentPtr->theory_ions_mass_c[n.thoe_id] << " "
//                 <<"spec_Ions_Mass:"<<setprecision(10)<< msalignPtr->getIonsMassContainer()[n.mono_id] << " "
//                 <<"index:"<< n.index << " "
//                 <<"ppm:"<< n.ppm << endl ;
//        }



    //计算mut_x系数(无用)
    tempPrsmArgumentPtr->mut_x = prsmScorePtr->analysis_mutx(tempPrsmArgumentPtr->ions_n, tempPrsmArgumentPtr->ions_c,
                                                             msalignPtr->getIonsMassContainer(),
                                                             ptms_mass, tempPrsmArgumentPtr->var_mod_mass);

    //计算互补离子个数
    tempPrsmArgumentPtr->complement_ions_number = prsmScorePtr->get_complementIons_score(tempPrsmArgumentPtr->ions_n, tempPrsmArgumentPtr->ions_c,
                                                                                         tempPrsmArgumentPtr->theory_ions_mass_c.size())
                                                  + 1;

    //计算连续的离子个数，返回的不是个数，而是经过计算的seqScore
    double seqScore = prsmScorePtr->get_seqIons_score(tempPrsmArgumentPtr->ions_n, tempPrsmArgumentPtr->ions_c,
                                                      tempPrsmArgumentPtr->seq_score_count);

    tempPrsmArgumentPtr->pairIonsNub = prsmScorePtr->get_SeqIons_pair(tempPrsmArgumentPtr->ions_n,
                                                                      tempPrsmArgumentPtr->ions_c);



    ///打分

    tempPrsmArgumentPtr->prsm_score = sqrt(((double)tempPrsmArgumentPtr->pairIonsNub * 2)
                                           / (double)msalignPtr->getIonsMassContainer().size());                          //计算连续离子得分
    tempPrsmArgumentPtr->matchFragmentIonsSize += mylib::data_stream::get_fragment_ions_nub(tempPrsmArgumentPtr->ions_n); //得到匹配离子数
    tempPrsmArgumentPtr->matchFragmentIonsSize += mylib::data_stream::get_fragment_ions_nub(tempPrsmArgumentPtr->ions_c); //得到匹配离子数

    tempPrsmArgumentPtr->prsm_score += ((double)tempPrsmArgumentPtr->matchFragmentIonsSize
                                        / (double)msalignPtr->getIonsMassContainer().size());


    prsmScorePtr->get_score(tempPrsmArgumentPtr->ions_n, tempPrsmArgumentPtr->ref_score_n, ptms_mass,
                            tempPrsmArgumentPtr->peer_score_n, tempPrsmArgumentPtr->score_n, 0.1,
                            tempPrsmArgumentPtr->var_mod_mass);
    prsmScorePtr->get_score(tempPrsmArgumentPtr->ions_c, tempPrsmArgumentPtr->ref_score_c, ptms_mass,
                            tempPrsmArgumentPtr->peer_score_c, tempPrsmArgumentPtr->score_c, 0.1,
                            tempPrsmArgumentPtr->var_mod_mass);

    //验证ptm
    double mass = 0;
    for (int xi = 0; xi < tempPrsmArgumentPtr->ptms.size(); ++xi) {
        mass += modptr->mapStrMass[modptr->analysis(tempPrsmArgumentPtr->ptms[xi].mod_mass)];
    }
    tempPrsmArgumentPtr->get_true_modMass(mass);
    if (fabs(mass - tempPrsmArgumentPtr->var_mod_mass) > 3) {
        return;
    }
    if (mass == 0 && tempPrsmArgumentPtr->var_mod_mass != 0) {
        return;
    }
    double ions_peaks_proportion = ((double)tempPrsmArgumentPtr->ions_n.size() + tempPrsmArgumentPtr->ions_c.size())
                                   / (double)msalignPtr->getIonsMassContainer().size() * 10; //计算峰值占比


    //去重鉴定离子数
    int dedup_lenn=0,dedup_lenc=0;

    int last_ion;
    if(tempPrsmArgumentPtr->ions_n.size() !=0) {
        for (int i = 0; i < tempPrsmArgumentPtr->ions_n.size(); i++) {
            if (i == 0) {
                last_ion = tempPrsmArgumentPtr->ions_n[i].thoe_id;
                dedup_lenn++;
            } else {
                if (tempPrsmArgumentPtr->ions_n[i].thoe_id != last_ion) {
                    dedup_lenn++;
                }
                last_ion = tempPrsmArgumentPtr->ions_n[i].thoe_id;
            }
        }
    }
    if(tempPrsmArgumentPtr->ions_c.size() !=0) {
        for (int i = 0; i < tempPrsmArgumentPtr->ions_c.size(); i++) {
            if (i == 0) {
                last_ion = tempPrsmArgumentPtr->ions_c[i].thoe_id;
                dedup_lenc++;
            } else {
                if (tempPrsmArgumentPtr->ions_c[i].thoe_id != last_ion) {
                    dedup_lenc++;
                }
                last_ion = tempPrsmArgumentPtr->ions_c[i].thoe_id;
            }
        }
    }



    ///截断后比例
    double truPro = (double)(proteinPtr->getProteinSequence().size()-tempPrsmArgumentPtr->cut_location_n - tempPrsmArgumentPtr->cut_location_c)/(double)proteinPtr->getProteinSequence().size();

    double tempPeakPro = ((double)tempPrsmArgumentPtr->ions_n.size() + (double)tempPrsmArgumentPtr->ions_c.size()+dedup_lenc+dedup_lenn)/(double)msalignPtr->getIonsMassContainer().size();
    tempPeakPro = exp(tempPeakPro);


    int min_tru = max(min(tempPrsmArgumentPtr->cut_location_n,tempPrsmArgumentPtr->cut_location_c),5);
    double peakionpro=(double)(tempPrsmArgumentPtr->ions_n.size()+tempPrsmArgumentPtr->ions_c.size()+dedup_lenc+dedup_lenn)/(double)proteinPtr->getProteinSequence().size();

    double subCoe = truPro * tempPeakPro * peakionpro/(double)min_tru;


    tempPrsmArgumentPtr->get_ncter_score(tempPrsmArgumentPtr->mut_x, tempPrsmArgumentPtr->complement_ions_number, subCoe,
                                   tempPrsmArgumentPtr->prsm_score, ions_peaks_proportion); // 计算最终得分

    ///测试打分
    //完全匹配离子数
    int com_match_n = 0,com_match_c = 0;
    for(auto ion:tempPrsmArgumentPtr->ions_n){
        if(fabs(ion.index) < 0.2)
            com_match_n++;
    }
    for(auto ion:tempPrsmArgumentPtr->ions_c){
        if(fabs(ion.index) < 0.2)
            com_match_c++;
    }
    //去重连续离子数
    int seq_match_n = 1,seq_match_c = 1;
    seq_match_n = prsmScorePtr->get_seq_count(tempPrsmArgumentPtr->ions_n);
    seq_match_c = prsmScorePtr->get_seq_count(tempPrsmArgumentPtr->ions_c);
    //总鉴定离子数
    int lenn =  tempPrsmArgumentPtr->theory_ions_mass_n.size();
    int lenc =  tempPrsmArgumentPtr->theory_ions_mass_c.size();
    double ion_len_pp =  ((double)tempPrsmArgumentPtr->ions_n.size() + tempPrsmArgumentPtr->ions_c.size())/lenn;
    //离子匹配误差加和
    double b_index_sum=0,y_index_sum=0,index_sum=0;
    double b_avg_index=0,y_avg_index=0,avg_index=0;
    for(int i = 0;i < tempPrsmArgumentPtr->ions_n.size();i++){
        b_index_sum += fabs(tempPrsmArgumentPtr->ions_n[i].index);
    }
    for(int i = 0;i < tempPrsmArgumentPtr->ions_c.size();i++){
        y_index_sum += fabs(tempPrsmArgumentPtr->ions_c[i].index);
    }
    index_sum = b_index_sum + y_index_sum;
    if(tempPrsmArgumentPtr->ions_n.size()!=0)
        b_avg_index = b_index_sum / (double)tempPrsmArgumentPtr->ions_n.size();
    if(tempPrsmArgumentPtr->ions_c.size()!=0)
        y_avg_index = y_index_sum / (double)tempPrsmArgumentPtr->ions_c.size();
    if( (tempPrsmArgumentPtr->ions_c.size()+tempPrsmArgumentPtr->ions_n.size())  != 0)
        avg_index = index_sum /(double)(tempPrsmArgumentPtr->ions_n.size()+tempPrsmArgumentPtr->ions_c.size());
    //被修饰匹配峰数
    int nter_mod_match_peak = tempPrsmArgumentPtr->ions_n.size() - com_match_n;
    int cter_mod_match_peak = tempPrsmArgumentPtr->ions_c.size() - com_match_c;
    //平均ppm
    double nter_avg_ppm,cter_avg_ppm,avg_ppm;
    double n_ppm_sum = 0,c_ppm_sum = 0,ppm_sum = 0;
    for(auto temp_ion:tempPrsmArgumentPtr->ions_n){
        n_ppm_sum += fabs(temp_ion.ppm);
        ppm_sum += fabs(temp_ion.ppm);
    }
    nter_avg_ppm = n_ppm_sum/(double)tempPrsmArgumentPtr->ions_n.size();
    for(auto temp_ion:tempPrsmArgumentPtr->ions_c){
        c_ppm_sum += fabs(temp_ion.ppm);
        ppm_sum += fabs(temp_ion.ppm);
    }
    cter_avg_ppm = c_ppm_sum/(double)tempPrsmArgumentPtr->ions_c.size();
    avg_ppm = ppm_sum/(double)(tempPrsmArgumentPtr->ions_n.size()+tempPrsmArgumentPtr->ions_c.size());
    //蛋白变体长度、原始蛋白序列长度
    int cutProLen = proteinPtr->getProteinSequence().size()-tempPrsmArgumentPtr->cut_location_n - tempPrsmArgumentPtr->cut_location_c;
    int oriProLen = proteinPtr->getProteinSequence().size();
    ///保存特征
    ///feature1：N端完全匹配离子数
    tempPrsmArgumentPtr->feature1 = com_match_n;
    ///feature2：C端完全匹配离子数
    tempPrsmArgumentPtr->feature2 = com_match_c;
    ///feature3：互补离子数
    tempPrsmArgumentPtr->feature3 = tempPrsmArgumentPtr->complement_ions_number;
    ///feature4:N端去重鉴定离子数
    tempPrsmArgumentPtr->feature4 = dedup_lenn;
    ///feature5:C端去重鉴定离子数
    tempPrsmArgumentPtr->feature5 = dedup_lenc;
    ///feature6:去重鉴定离子个数/蛋白变体长度
    tempPrsmArgumentPtr->feature6 = (double)(dedup_lenc + dedup_lenn)/(double)cutProLen;
    ///feature7:N端去重鉴定离子个数/蛋白变体长度
    tempPrsmArgumentPtr->feature7 = (double)(dedup_lenn)/(double)cutProLen;
    ///feature8:C端去重鉴定离子个数/蛋白变体长度
    tempPrsmArgumentPtr->feature8 = (double)(dedup_lenc)/(double)cutProLen;
    ///feature9:n端连续匹配离子个数（自写）
    tempPrsmArgumentPtr->feature9 = seq_match_n;
    ///feature10:c端连续匹配离子个数（自写）
    tempPrsmArgumentPtr->feature10 = seq_match_c;
    ///feature11:N端匹配峰数
    tempPrsmArgumentPtr->feature11 = tempPrsmArgumentPtr->ions_n.size();
    ///feature12:C端匹配峰数
    tempPrsmArgumentPtr->feature12 = tempPrsmArgumentPtr->ions_c.size();
    ///feature13:匹配峰数
    tempPrsmArgumentPtr->feature13 = tempPrsmArgumentPtr->peaks_match;
    ///feature14:离子对数量
    tempPrsmArgumentPtr->feature14 = tempPrsmArgumentPtr->pairIonsNub;
    ///feature15:完全匹配离子数
    tempPrsmArgumentPtr->feature15 = com_match_n + com_match_c;
    ///feature16:去重鉴定离子总数
    tempPrsmArgumentPtr->feature16 = dedup_lenn + dedup_lenc;
    ///feature17:连续匹配离子个数（自写）
    tempPrsmArgumentPtr->feature17 = seq_match_n + seq_match_c;
    ///feature18:N端被修饰匹配峰数
    tempPrsmArgumentPtr->feature18 = nter_mod_match_peak;
    ///feature19:C端被修饰匹配峰数
    tempPrsmArgumentPtr->feature19 = cter_mod_match_peak;
    ///feature20:被修饰匹配峰和
    tempPrsmArgumentPtr->feature20 = nter_mod_match_peak + cter_mod_match_peak;
    ///feature21:N端匹配峰平均ppm
    tempPrsmArgumentPtr->feature21 = nter_avg_ppm;
    ///feature22:C端匹配峰平均ppm
    tempPrsmArgumentPtr->feature22 = cter_avg_ppm;
    ///feature23:匹配峰平均ppm
    tempPrsmArgumentPtr->feature23 = avg_ppm;
    ///feature24:b离子丰度得分
    tempPrsmArgumentPtr->feature24 = prsmScorePtr->get_abundance_score(tempPrsmArgumentPtr->ions_n
        ,msalignPtr->getIonsMassContainer(),msalignPtr->getKeyMsValueAbRadio());
    ///feature25:y离子丰度得分
    tempPrsmArgumentPtr->feature25 = prsmScorePtr->get_abundance_score(tempPrsmArgumentPtr->ions_c
        ,msalignPtr->getIonsMassContainer(),msalignPtr->getKeyMsValueAbRadio());
    ///feature26:by离子丰度得分和
    tempPrsmArgumentPtr->feature26 = tempPrsmArgumentPtr->feature24 + tempPrsmArgumentPtr->feature25;
    ///feature27:b离子丰度基于前体质量动态得分
    tempPrsmArgumentPtr->feature27 = prsmScorePtr->get_dyn_abu_score(tempPrsmArgumentPtr->ions_n
        ,msalignPtr->getIonsMassContainer(),msalignPtr->getKeyMsValueAbRadio(),msalignPtr->getPrecursorMass());
    ///feature28:y离子丰度基于前体质量动态得分
    tempPrsmArgumentPtr->feature28 = prsmScorePtr->get_dyn_abu_score(tempPrsmArgumentPtr->ions_c
        ,msalignPtr->getIonsMassContainer(),msalignPtr->getKeyMsValueAbRadio(),msalignPtr->getPrecursorMass());
    ///feature29:by离子丰度基于前体质量动态得分和
    tempPrsmArgumentPtr->feature29 = tempPrsmArgumentPtr->feature27 + tempPrsmArgumentPtr->feature28;
    ///feature30:by平均离子丰度得分和
    tempPrsmArgumentPtr->feature30 = prsmScorePtr->get_avg_abundance_score(tempPrsmArgumentPtr->ions_n
        ,tempPrsmArgumentPtr->ions_c,msalignPtr->getIonsMassContainer(),msalignPtr->getKeyMsValueAbRadio());
    ///feature31:by平均离子丰度基于前体质量动态得分
    tempPrsmArgumentPtr->feature31 = prsmScorePtr->get_avg_dyn_abu_score(tempPrsmArgumentPtr->ions_n
        ,tempPrsmArgumentPtr->ions_c,msalignPtr->getIonsMassContainer(),msalignPtr->getKeyMsValueAbRadio()
        ,msalignPtr->getPrecursorMass());
    ///feature32:N端匹配峰数/蛋白变体长度
    tempPrsmArgumentPtr->feature32 = (double)tempPrsmArgumentPtr->ions_n.size()/(double)cutProLen;
    ///feature33:C端匹配峰数/蛋白变体长度
    tempPrsmArgumentPtr->feature33 = (double)tempPrsmArgumentPtr->ions_c.size()/(double)cutProLen;
    ///feature34:匹配峰数/蛋白变体长度
    tempPrsmArgumentPtr->feature34 = (double)(tempPrsmArgumentPtr->ions_n.size() + tempPrsmArgumentPtr->ions_c.size())/(double)cutProLen;
    ///feature35:匹配峰数/质谱峰总数
    tempPrsmArgumentPtr->feature35 = (double)(tempPrsmArgumentPtr->ions_n.size() + tempPrsmArgumentPtr->ions_c.size())/(double)msalignPtr->getIonsMassContainer().size();
    ///feature36:匹配离子数/质谱峰总数
    tempPrsmArgumentPtr->feature36 = (double)(dedup_lenn + dedup_lenc)/(double)msalignPtr->getIonsMassContainer().size();
    ///feature37:完全匹配离子数/质谱峰总数
    tempPrsmArgumentPtr->feature37 = (double)(com_match_n + com_match_c)/(double)msalignPtr->getIonsMassContainer().size();
    ///feature38:被修饰匹配峰数/质谱峰总数
    tempPrsmArgumentPtr->feature38 = (double)(nter_mod_match_peak + cter_mod_match_peak)/(double)msalignPtr->getIonsMassContainer().size();
    ///feature39:连续匹配离子数/质谱峰总数
    tempPrsmArgumentPtr->feature39 = (double)(seq_match_n + seq_match_c)/(double)msalignPtr->getIonsMassContainer().size();
    ///feature40:N端截断长度/原蛋白序列长度
    tempPrsmArgumentPtr->feature40 = (double)tempPrsmArgumentPtr->cut_location_n/(double)oriProLen;
    ///feature41:C端截断长度/原蛋白序列长度
    tempPrsmArgumentPtr->feature41 = (double)tempPrsmArgumentPtr->cut_location_c/(double)oriProLen;
    ///feature42:MIN（N端截断长度，C端截断长度）/原蛋白序列长度
    tempPrsmArgumentPtr->feature42 = (double)min(tempPrsmArgumentPtr->cut_location_n,tempPrsmArgumentPtr->cut_location_c)/(double)oriProLen;
    ///feature43:蛋白变体长度/原蛋白序列长度
    tempPrsmArgumentPtr->feature43 = (double)cutProLen/(double)oriProLen;
    ///feature44:N端C端匹配峰比值
    tempPrsmArgumentPtr->feature44 = (double)min(tempPrsmArgumentPtr->ions_n.size(),tempPrsmArgumentPtr->ions_c.size())/(double)max(tempPrsmArgumentPtr->ions_n.size(),tempPrsmArgumentPtr->ions_c.size());
    ///feature45:N端C端匹配离子数比值
    tempPrsmArgumentPtr->feature45 = (double)min(dedup_lenn,dedup_lenc)/(double)max(dedup_lenn,dedup_lenc);
    ///feature46:N端C端完全匹配离子数比值
    tempPrsmArgumentPtr->feature46 = (double)min(com_match_n,com_match_c)/(double)max(com_match_n,com_match_c);
    ///feature47:N端C端被修饰匹配峰数比值
    tempPrsmArgumentPtr->feature47 = (double)min(nter_mod_match_peak,cter_mod_match_peak)/(double)max(nter_mod_match_peak,cter_mod_match_peak);
    ///feature48:N端C端平均ppm比值
    tempPrsmArgumentPtr->feature48 = (double)min(nter_avg_ppm,cter_avg_ppm)/(double)max(nter_avg_ppm,cter_avg_ppm);
    ///feature49:完全匹配离子数/蛋白变体长度
    tempPrsmArgumentPtr->feature49 = (double)(com_match_n + com_match_c)/(double)cutProLen;
    ///feature50:连续匹配离子数/蛋白变体长度
    tempPrsmArgumentPtr->feature50 = (double)(seq_match_n + seq_match_c)/(double)cutProLen;

}

int prsm_unknown_ptms::getCurMsSize() {
    return cur_ms_size;
}

prsm_cwdtw *prsm_unknown_ptms::getPrsmCwdtwPtr() {
    return prsmCwdtwPtr;
}

prsm_scope_alignment *prsm_unknown_ptms::getPrsmScopeAlignmentPtr() {
    return prsmScopeAlignmentPtr;
}

prsm_alignment_filter *prsm_unknown_ptms::getPrsmAlignmentFilterPtr() {
    return prsmAlignmentFilterPtr;
}

prsm_ions_location *prsm_unknown_ptms::getPrsmIonsLocationPtr() {
    return prsmIonsLocationPtr;
}

prsm_score *prsm_unknown_ptms::getPrsmScorePtr() {
    return prsmScorePtr;
}

prsm_msg_write *prsm_unknown_ptms::getPrsmMsgWritePtr() {
    return prsmMsgWritePtr;
}

void prsm_unknown_ptms::setCurMsSize(int curMsSize) {
    cur_ms_size = curMsSize;
}

void prsm_unknown_ptms::setPrsmCwdtwPtr(prsm_cwdtw *prsmCwdtwPtr) {
    prsm_unknown_ptms::prsmCwdtwPtr = prsmCwdtwPtr;
}

void prsm_unknown_ptms::setPrsmScopeAlignmentPtr(prsm_scope_alignment *prsmScopeAlignmentPtr) {
    prsm_unknown_ptms::prsmScopeAlignmentPtr = prsmScopeAlignmentPtr;
}

void prsm_unknown_ptms::setPrsmAlignmentFilterPtr(prsm_alignment_filter *prsmAlignmentFilterPtr) {
    prsm_unknown_ptms::prsmAlignmentFilterPtr = prsmAlignmentFilterPtr;
}

void prsm_unknown_ptms::setPrsmIonsLocationPtr(prsm_ions_location *prsmIonsLocationPtr) {
    prsm_unknown_ptms::prsmIonsLocationPtr = prsmIonsLocationPtr;
}

void prsm_unknown_ptms::setPrsmScorePtr(prsm_score *prsmScorePtr) {
    prsm_unknown_ptms::prsmScorePtr = prsmScorePtr;
}

void prsm_unknown_ptms::setPrsmMsgWritePtr(prsm_msg_write *prsmMsgWritePtr) {
    prsm_unknown_ptms::prsmMsgWritePtr = prsmMsgWritePtr;
}

void prsm_unknown_ptms::prsm_process_var(prsm::modify_ptr modptr,
                                         protein_ptr proteinPtr,
                                         msalign_ptr msalignPtr,
                                         prsm::temp_prsm_argument_ptr tempPrsmArgumentPtr,
                                         vector<double> &ptms_mass,
                                         prsm::config_argument_ptr prsmArg) {
    //Coarse alignment
    //    cout << " cwdtw - start " << endl ;
    prsmCwdtwPtr->cwdtw_algorithm(tempPrsmArgumentPtr->theory_ions_mass_n, msalignPtr->getIonsMassContainer(),
                                  tempPrsmArgumentPtr->alignment_ions_n, tempPrsmArgumentPtr->ref_score_n,
                                  tempPrsmArgumentPtr->peer_score_n);
    prsmCwdtwPtr->cwdtw_algorithm(tempPrsmArgumentPtr->theory_ions_mass_c, msalignPtr->getIonsMassContainer(),
                                  tempPrsmArgumentPtr->alignment_ions_c, tempPrsmArgumentPtr->ref_score_c,
                                  tempPrsmArgumentPtr->peer_score_c);
    //    cout << " cwdtw - end " << endl ;

    //    cout << " scopeAlignment - start " << endl ;
    prsmScopeAlignmentPtr->scopeAlignment(tempPrsmArgumentPtr->theory_ions_mass_n,
                                          msalignPtr->getIonsMassContainer(),
                                          tempPrsmArgumentPtr->alignment_ions_n,
                                          ptms_mass, tempPrsmArgumentPtr->var_mod_mass, prsmArg->ptmMassMax);
    prsmScopeAlignmentPtr->scopeAlignment(tempPrsmArgumentPtr->theory_ions_mass_c,
                                          msalignPtr->getIonsMassContainer(),
                                          tempPrsmArgumentPtr->alignment_ions_c,
                                          ptms_mass, tempPrsmArgumentPtr->var_mod_mass, prsmArg->ptmMassMax);
    //    cout << " scopeAlignment - end " << endl ;

    prsmAlignmentFilterPtr->get_mass_sub(tempPrsmArgumentPtr->theory_ions_mass_n, tempPrsmArgumentPtr->theory_ions_mass_c,
                                         msalignPtr->getIonsMassContainer(),
                                         tempPrsmArgumentPtr->sub_n, tempPrsmArgumentPtr->sub_c,
                                         tempPrsmArgumentPtr->alignment_ions_n, tempPrsmArgumentPtr->alignment_ions_c);

    //change
    prsmAlignmentFilterPtr->index_ions(tempPrsmArgumentPtr->sub_n, tempPrsmArgumentPtr->filter_ions_n,
                                       ptms_mass, tempPrsmArgumentPtr->var_mod_mass);

    prsmAlignmentFilterPtr->index_ions(tempPrsmArgumentPtr->sub_c, tempPrsmArgumentPtr->filter_ions_c,
                                       ptms_mass, tempPrsmArgumentPtr->var_mod_mass);

    prsmAlignmentFilterPtr->get_ppm_unknown(tempPrsmArgumentPtr->theory_ions_mass_n, msalignPtr->getIonsMassContainer(),
                                            tempPrsmArgumentPtr->alignment_ions_n, tempPrsmArgumentPtr->filter_ions_n,
                                            ptms_mass, tempPrsmArgumentPtr->ppm_value_n, tempPrsmArgumentPtr->var_mod_mass);
    prsmAlignmentFilterPtr->get_ppm_unknown(tempPrsmArgumentPtr->theory_ions_mass_c, msalignPtr->getIonsMassContainer(),
                                            tempPrsmArgumentPtr->alignment_ions_c, tempPrsmArgumentPtr->filter_ions_c,
                                            ptms_mass, tempPrsmArgumentPtr->ppm_value_c, tempPrsmArgumentPtr->var_mod_mass);

    prsmAlignmentFilterPtr->get_ions_peaks(tempPrsmArgumentPtr->merge_ions_n, tempPrsmArgumentPtr->alignment_ions_n,
                                           tempPrsmArgumentPtr->filter_ions_n, tempPrsmArgumentPtr->ppm_value_n,
                                           prsmArg->min_ppm_value);
    prsmAlignmentFilterPtr->get_ions_peaks(tempPrsmArgumentPtr->merge_ions_c, tempPrsmArgumentPtr->alignment_ions_c,
                                           tempPrsmArgumentPtr->filter_ions_c, tempPrsmArgumentPtr->ppm_value_c,
                                           prsmArg->min_ppm_value);

    vector<node> N_ions_c_n, N_ions_n_c;
    vector<modi> N_Ptms;
    vector<node> C_ions_c_n, C_ions_n_c;
    vector<modi> C_Ptms;
    //change
    prsmIonsLocationPtr->ions_location_match_N_v2(N_ions_c_n, N_ions_n_c, tempPrsmArgumentPtr->merge_ions_n,
                                                  tempPrsmArgumentPtr->merge_ions_c, N_Ptms, ptms_mass,
                                                  tempPrsmArgumentPtr->get_lenth(tempPrsmArgumentPtr->theory_ions_mass_n.size()),
                                                  tempPrsmArgumentPtr->var_mod_mass);

    prsmIonsLocationPtr->ions_location_match_C_v2(C_ions_c_n, C_ions_n_c, tempPrsmArgumentPtr->merge_ions_n,
                                                  tempPrsmArgumentPtr->merge_ions_c,
                                                  C_Ptms, ptms_mass,
                                                  tempPrsmArgumentPtr->get_lenth(tempPrsmArgumentPtr->theory_ions_mass_n.size()),
                                                  tempPrsmArgumentPtr->var_mod_mass);

    prsmIonsLocationPtr->confirm_final_ions(N_ions_c_n, N_ions_n_c, N_Ptms,
                                            C_ions_c_n, C_ions_n_c, C_Ptms,
                                            tempPrsmArgumentPtr);

    if (tempPrsmArgumentPtr->ions_n.size() == 0 && tempPrsmArgumentPtr->ions_c.size() == 0) {
        return;
    }

    if(N_Ptms.size()>3 || C_Ptms.size()>3){
        return;
    }

    //计算mut_x系数(无用)
    tempPrsmArgumentPtr->mut_x =
        prsmScorePtr->analysis_mutx(tempPrsmArgumentPtr->ions_n,
                                    tempPrsmArgumentPtr->ions_c,
                                    msalignPtr->getIonsMassContainer(),
                                    ptms_mass, tempPrsmArgumentPtr->var_mod_mass);

    //计算互补离子个数
    tempPrsmArgumentPtr->complement_ions_number =
        prsmScorePtr->get_complementIons_score(tempPrsmArgumentPtr->ions_n,
                                               tempPrsmArgumentPtr->ions_c,
                                               tempPrsmArgumentPtr->theory_ions_mass_c.size()) + 1;

    //计算连续的离子个数，返回的不是个数，而是经过计算的seqScore
    double seqScore = prsmScorePtr->get_seqIons_score(tempPrsmArgumentPtr->ions_n, tempPrsmArgumentPtr->ions_c,
                                                      tempPrsmArgumentPtr->seq_score_count);

    tempPrsmArgumentPtr->pairIonsNub =
        prsmScorePtr->get_SeqIons_pair(tempPrsmArgumentPtr->ions_n,
                                       tempPrsmArgumentPtr->ions_c);

    tempPrsmArgumentPtr->prsm_score = sqrt(((double)tempPrsmArgumentPtr->pairIonsNub * 2)
                                           / (double)msalignPtr->getIonsMassContainer().size());                          //计算连续离子得分
    tempPrsmArgumentPtr->matchFragmentIonsSize += mylib::data_stream::get_fragment_ions_nub(tempPrsmArgumentPtr->ions_n); //得到匹配离子数
    tempPrsmArgumentPtr->matchFragmentIonsSize += mylib::data_stream::get_fragment_ions_nub(tempPrsmArgumentPtr->ions_c); //得到匹配离子数

    tempPrsmArgumentPtr->prsm_score += ((double)tempPrsmArgumentPtr->matchFragmentIonsSize
                                        / (double)msalignPtr->getIonsMassContainer().size());

    prsmScorePtr->get_score(tempPrsmArgumentPtr->ions_n, tempPrsmArgumentPtr->ref_score_n, ptms_mass,
                            tempPrsmArgumentPtr->peer_score_n, tempPrsmArgumentPtr->score_n, 0.1,
                            tempPrsmArgumentPtr->var_mod_mass);
    prsmScorePtr->get_score(tempPrsmArgumentPtr->ions_c, tempPrsmArgumentPtr->ref_score_c, ptms_mass,
                            tempPrsmArgumentPtr->peer_score_c, tempPrsmArgumentPtr->score_c, 0.1,
                            tempPrsmArgumentPtr->var_mod_mass);

    //验证ptm
    double mass = 0;
    for (int xi = 0; xi < tempPrsmArgumentPtr->ptms.size(); ++xi) {
        mass += modptr->mapStrMass[modptr->analysis(tempPrsmArgumentPtr->ptms[xi].mod_mass)];
    }
    tempPrsmArgumentPtr->get_true_modMass(mass);
    if (fabs(mass - tempPrsmArgumentPtr->var_mod_mass) > 3) {
        return;
    }
    if (mass == 0 && tempPrsmArgumentPtr->var_mod_mass != 0) {
        return;
    }
    double ions_peaks_proportion = ((double)tempPrsmArgumentPtr->ions_n.size() + tempPrsmArgumentPtr->ions_c.size())
                                   / (double)msalignPtr->getIonsMassContainer().size() * 10; //计算峰值占比

    int subCoe = (tempPrsmArgumentPtr->cut_location_n / 10) + (tempPrsmArgumentPtr->cut_location_c / 10); //计算截断总数

    tempPrsmArgumentPtr->get_score(tempPrsmArgumentPtr->mut_x, tempPrsmArgumentPtr->complement_ions_number, subCoe,
                                   tempPrsmArgumentPtr->prsm_score, ions_peaks_proportion); // 计算最终得分

    ///完全匹配离子数
    int com_match_n = 0,com_match_c = 0;
    for(auto ion:tempPrsmArgumentPtr->ions_n){
        if(fabs(ion.index) < 0.2)
            com_match_n++;
    }
    for(auto ion:tempPrsmArgumentPtr->ions_c){
        if(fabs(ion.index) < 0.2)
            com_match_c++;
    }
    ///去重连续离子数
    int seq_match_n = 1,seq_match_c = 1;
    seq_match_n = prsmScorePtr->get_seq_count(tempPrsmArgumentPtr->ions_n);
    seq_match_c = prsmScorePtr->get_seq_count(tempPrsmArgumentPtr->ions_c);
    ///总鉴定离子数
    int lenn =  tempPrsmArgumentPtr->theory_ions_mass_n.size();
    int lenc =  tempPrsmArgumentPtr->theory_ions_mass_c.size();
    double ion_len_pp =  ((double)tempPrsmArgumentPtr->ions_n.size() + tempPrsmArgumentPtr->ions_c.size())/lenn;
    ///离子匹配误差加和
    double b_index_sum=0,y_index_sum=0,index_sum=0;
    double b_avg_index=0,y_avg_index=0,avg_index=0;
    for(int i = 0;i < tempPrsmArgumentPtr->ions_n.size();i++){
        b_index_sum += fabs(tempPrsmArgumentPtr->ions_n[i].index);
    }
    for(int i = 0;i < tempPrsmArgumentPtr->ions_c.size();i++){
        y_index_sum += fabs(tempPrsmArgumentPtr->ions_c[i].index);
    }
    index_sum = b_index_sum + y_index_sum;
    if(tempPrsmArgumentPtr->ions_n.size()!=0)
        b_avg_index = b_index_sum / (double)tempPrsmArgumentPtr->ions_n.size();
    if(tempPrsmArgumentPtr->ions_c.size()!=0)
        y_avg_index = y_index_sum / (double)tempPrsmArgumentPtr->ions_c.size();
    if( (tempPrsmArgumentPtr->ions_c.size()+tempPrsmArgumentPtr->ions_n.size())  != 0)
        avg_index = index_sum /(double)(tempPrsmArgumentPtr->ions_n.size()+tempPrsmArgumentPtr->ions_c.size());
    ///被修饰匹配峰数
    int nter_mod_match_peak = tempPrsmArgumentPtr->ions_n.size() - com_match_n;
    int cter_mod_match_peak = tempPrsmArgumentPtr->ions_c.size() - com_match_c;
    ///平均ppm
    double nter_avg_ppm,cter_avg_ppm,avg_ppm;
    double n_ppm_sum = 0,c_ppm_sum = 0,ppm_sum = 0;
    for(auto temp_ion:tempPrsmArgumentPtr->ions_n){
        n_ppm_sum += fabs(temp_ion.ppm);
        ppm_sum += fabs(temp_ion.ppm);
    }
    nter_avg_ppm = n_ppm_sum/(double)tempPrsmArgumentPtr->ions_n.size();
    for(auto temp_ion:tempPrsmArgumentPtr->ions_c){
        c_ppm_sum += fabs(temp_ion.ppm);
        ppm_sum += fabs(temp_ion.ppm);
    }
    cter_avg_ppm = c_ppm_sum/(double)tempPrsmArgumentPtr->ions_c.size();
    avg_ppm = ppm_sum/(double)(tempPrsmArgumentPtr->ions_n.size()+tempPrsmArgumentPtr->ions_c.size());
    ///蛋白变体长度、原始蛋白序列长度
    int cutProLen = proteinPtr->getProteinSequence().size()-tempPrsmArgumentPtr->cut_location_n - tempPrsmArgumentPtr->cut_location_c;
    int oriProLen = proteinPtr->getProteinSequence().size();
    ///去重鉴定离子数
    int dedup_lenn=0,dedup_lenc=0;

    int last_ion;
    if(tempPrsmArgumentPtr->ions_n.size() !=0) {
        for (int i = 0; i < tempPrsmArgumentPtr->ions_n.size(); i++) {
            if (i == 0) {
                last_ion = tempPrsmArgumentPtr->ions_n[i].thoe_id;
                dedup_lenn++;
            } else {
                if (tempPrsmArgumentPtr->ions_n[i].thoe_id != last_ion) {
                    dedup_lenn++;
                }
                last_ion = tempPrsmArgumentPtr->ions_n[i].thoe_id;
            }
        }
    }
    if(tempPrsmArgumentPtr->ions_c.size() !=0) {
        for (int i = 0; i < tempPrsmArgumentPtr->ions_c.size(); i++) {
            if (i == 0) {
                last_ion = tempPrsmArgumentPtr->ions_c[i].thoe_id;
                dedup_lenc++;
            } else {
                if (tempPrsmArgumentPtr->ions_c[i].thoe_id != last_ion) {
                    dedup_lenc++;
                }
                last_ion = tempPrsmArgumentPtr->ions_c[i].thoe_id;
            }
        }
    }



///保存特征
    ///feature1：N端完全匹配离子数
    tempPrsmArgumentPtr->feature1 = com_match_n;
    ///feature2：C端完全匹配离子数
    tempPrsmArgumentPtr->feature2 = com_match_c;
    ///feature3：互补离子数
    tempPrsmArgumentPtr->feature3 = tempPrsmArgumentPtr->complement_ions_number;
    ///feature4:N端去重鉴定离子数
    tempPrsmArgumentPtr->feature4 = dedup_lenn;
    ///feature5:C端去重鉴定离子数
    tempPrsmArgumentPtr->feature5 = dedup_lenc;
    ///feature6:去重鉴定离子个数/蛋白变体长度
    tempPrsmArgumentPtr->feature6 = (double)(dedup_lenc + dedup_lenn)/(double)cutProLen;
    ///feature7:N端去重鉴定离子个数/蛋白变体长度
    tempPrsmArgumentPtr->feature7 = (double)(dedup_lenn)/(double)cutProLen;
    ///feature8:C端去重鉴定离子个数/蛋白变体长度
    tempPrsmArgumentPtr->feature8 = (double)(dedup_lenc)/(double)cutProLen;
    ///feature9:n端连续匹配离子个数（自写）
    tempPrsmArgumentPtr->feature9 = seq_match_n;
    ///feature10:c端连续匹配离子个数（自写）
    tempPrsmArgumentPtr->feature10 = seq_match_c;
    ///feature11:N端匹配峰数
    tempPrsmArgumentPtr->feature11 = tempPrsmArgumentPtr->ions_n.size();
    ///feature12:C端匹配峰数
    tempPrsmArgumentPtr->feature12 = tempPrsmArgumentPtr->ions_c.size();
    ///feature13:匹配峰数
    tempPrsmArgumentPtr->feature13 = tempPrsmArgumentPtr->peaks_match;
    ///feature14:离子对数量
    tempPrsmArgumentPtr->feature14 = tempPrsmArgumentPtr->pairIonsNub;
    ///feature15:完全匹配离子数
    tempPrsmArgumentPtr->feature15 = com_match_n + com_match_c;
    ///feature16:去重鉴定离子总数
    tempPrsmArgumentPtr->feature16 = dedup_lenn + dedup_lenc;
    ///feature17:连续匹配离子个数（自写）
    tempPrsmArgumentPtr->feature17 = seq_match_n + seq_match_c;
    ///feature18:N端被修饰匹配峰数
    tempPrsmArgumentPtr->feature18 = nter_mod_match_peak;
    ///feature19:C端被修饰匹配峰数
    tempPrsmArgumentPtr->feature19 = cter_mod_match_peak;
    ///feature20:被修饰匹配峰和
    tempPrsmArgumentPtr->feature20 = nter_mod_match_peak + cter_mod_match_peak;
    ///feature21:N端匹配峰平均ppm
    tempPrsmArgumentPtr->feature21 = nter_avg_ppm;
    ///feature22:C端匹配峰平均ppm
    tempPrsmArgumentPtr->feature22 = cter_avg_ppm;
    ///feature23:匹配峰平均ppm
    tempPrsmArgumentPtr->feature23 = avg_ppm;
    ///feature24:b离子丰度得分
    tempPrsmArgumentPtr->feature24 = prsmScorePtr->get_abundance_score(tempPrsmArgumentPtr->ions_n
        ,msalignPtr->getIonsMassContainer(),msalignPtr->getKeyMsValueAbRadio());
    ///feature25:y离子丰度得分
    tempPrsmArgumentPtr->feature25 = prsmScorePtr->get_abundance_score(tempPrsmArgumentPtr->ions_c
        ,msalignPtr->getIonsMassContainer(),msalignPtr->getKeyMsValueAbRadio());
    ///feature26:by离子丰度得分和
    tempPrsmArgumentPtr->feature26 = tempPrsmArgumentPtr->feature24 + tempPrsmArgumentPtr->feature25;
    ///feature27:b离子丰度基于前体质量动态得分
    tempPrsmArgumentPtr->feature27 = prsmScorePtr->get_dyn_abu_score(tempPrsmArgumentPtr->ions_n
        ,msalignPtr->getIonsMassContainer(),msalignPtr->getKeyMsValueAbRadio(),msalignPtr->getPrecursorMass());
    ///feature28:y离子丰度基于前体质量动态得分
    tempPrsmArgumentPtr->feature28 = prsmScorePtr->get_dyn_abu_score(tempPrsmArgumentPtr->ions_c
        ,msalignPtr->getIonsMassContainer(),msalignPtr->getKeyMsValueAbRadio(),msalignPtr->getPrecursorMass());
    ///feature29:by离子丰度基于前体质量动态得分和
    tempPrsmArgumentPtr->feature29 = tempPrsmArgumentPtr->feature27 + tempPrsmArgumentPtr->feature28;
    ///feature30:by平均离子丰度得分和
    tempPrsmArgumentPtr->feature30 = prsmScorePtr->get_avg_abundance_score(tempPrsmArgumentPtr->ions_n
        ,tempPrsmArgumentPtr->ions_c,msalignPtr->getIonsMassContainer(),msalignPtr->getKeyMsValueAbRadio());
    ///feature31:by平均离子丰度基于前体质量动态得分
    tempPrsmArgumentPtr->feature31 = prsmScorePtr->get_avg_dyn_abu_score(tempPrsmArgumentPtr->ions_n
        ,tempPrsmArgumentPtr->ions_c,msalignPtr->getIonsMassContainer(),msalignPtr->getKeyMsValueAbRadio()
        ,msalignPtr->getPrecursorMass());
    ///feature32:N端匹配峰数/蛋白变体长度
    tempPrsmArgumentPtr->feature32 = (double)tempPrsmArgumentPtr->ions_n.size()/(double)cutProLen;
    ///feature33:C端匹配峰数/蛋白变体长度
    tempPrsmArgumentPtr->feature33 = (double)tempPrsmArgumentPtr->ions_c.size()/(double)cutProLen;
    ///feature34:匹配峰数/蛋白变体长度
    tempPrsmArgumentPtr->feature34 = (double)(tempPrsmArgumentPtr->ions_n.size() + tempPrsmArgumentPtr->ions_c.size())/(double)cutProLen;
    ///feature35:匹配峰数/质谱峰总数
    tempPrsmArgumentPtr->feature35 = (double)(tempPrsmArgumentPtr->ions_n.size() + tempPrsmArgumentPtr->ions_c.size())/(double)msalignPtr->getIonsMassContainer().size();
    ///feature36:匹配离子数/质谱峰总数
    tempPrsmArgumentPtr->feature36 = (double)(dedup_lenn + dedup_lenc)/(double)msalignPtr->getIonsMassContainer().size();
    ///feature37:完全匹配离子数/质谱峰总数
    tempPrsmArgumentPtr->feature37 = (double)(com_match_n + com_match_c)/(double)msalignPtr->getIonsMassContainer().size();
    ///feature38:被修饰匹配峰数/质谱峰总数
    tempPrsmArgumentPtr->feature38 = (double)(nter_mod_match_peak + cter_mod_match_peak)/(double)msalignPtr->getIonsMassContainer().size();
    ///feature39:连续匹配离子数/质谱峰总数
    tempPrsmArgumentPtr->feature39 = (double)(seq_match_n + seq_match_c)/(double)msalignPtr->getIonsMassContainer().size();
    ///feature40:N端截断长度/原蛋白序列长度
    tempPrsmArgumentPtr->feature40 = (double)tempPrsmArgumentPtr->cut_location_n/(double)oriProLen;
    ///feature41:C端截断长度/原蛋白序列长度
    tempPrsmArgumentPtr->feature41 = (double)tempPrsmArgumentPtr->cut_location_c/(double)oriProLen;
    ///feature42:MIN（N端截断长度，C端截断长度）/原蛋白序列长度
    tempPrsmArgumentPtr->feature42 = (double)min(tempPrsmArgumentPtr->cut_location_n,tempPrsmArgumentPtr->cut_location_c)/(double)oriProLen;
    ///feature43:蛋白变体长度/原蛋白序列长度
    tempPrsmArgumentPtr->feature43 = (double)cutProLen/(double)oriProLen;
    ///feature44:N端C端匹配峰比值
    tempPrsmArgumentPtr->feature44 = (double)min(tempPrsmArgumentPtr->ions_n.size(),tempPrsmArgumentPtr->ions_c.size())/(double)max(tempPrsmArgumentPtr->ions_n.size(),tempPrsmArgumentPtr->ions_c.size());
    ///feature45:N端C端匹配离子数比值
    tempPrsmArgumentPtr->feature45 = (double)min(dedup_lenn,dedup_lenc)/(double)max(dedup_lenn,dedup_lenc);
    ///feature46:N端C端完全匹配离子数比值
    tempPrsmArgumentPtr->feature46 = (double)min(com_match_n,com_match_c)/(double)max(com_match_n,com_match_c);
    ///feature47:N端C端被修饰匹配峰数比值
    tempPrsmArgumentPtr->feature47 = (double)min(nter_mod_match_peak,cter_mod_match_peak)/(double)max(nter_mod_match_peak,cter_mod_match_peak);
    ///feature48:N端C端平均ppm比值
    tempPrsmArgumentPtr->feature48 = (double)min(nter_avg_ppm,cter_avg_ppm)/(double)max(nter_avg_ppm,cter_avg_ppm);
    ///feature49:完全匹配离子数/蛋白变体长度
    tempPrsmArgumentPtr->feature49 = (double)(com_match_n + com_match_c)/(double)cutProLen;
    ///feature50:连续匹配离子数/蛋白变体长度
    tempPrsmArgumentPtr->feature50 = (double)(seq_match_n + seq_match_c)/(double)cutProLen;



}