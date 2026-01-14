//
// Created by Administrator on 2023/3/23.
//

#include "prsm_varible_ptms.hpp"
#include "../mylib/calculate_lib.h"

bool select_cut_sortV2(c::arg::cutLocationPtrP a, c::arg::cutLocationPtrP b) {
    return a->cutLocationC + a->cutLocationN < b->cutLocationC + b->cutLocationN;
}


bool spectSortFunVV2(c::arg::cutLocationPtrP a, c::arg::cutLocationPtrP b) {
    return a->cutLocationC+a->cutLocationN < b->cutLocationC+b->cutLocationN;
}



void prsm_varible_ptms::terminal_trunction(protein_ptr proteinPtr, msalign_ptr msalignPtr,
                                           prsm::config_argument_ptr prsmArg,
                                           const vector<double> &ptmMass, const double &precursorMass,
                                           map<int, int> &cutMap) {
    if (proteinPtr->getTheoryMassCBy().empty() || proteinPtr->getTheoryMassCCz().empty() || ptmMass.empty()) { //条件判断
        return;
    }
    double minPrecursorMass = precursorMass - prsmArg->ptmMassMax;
    double maxPrecursorMass = precursorMass + prsmArg->ptmMassMax;
    int leftRange = 0;
    int seqLog = proteinPtr->getProteinSequence().length() - 1;
    double subMass = 0;
    double hMass = prsmArg->H;
    double ppm = prsmArg->min_ppm_value;
    vector<c::arg::cutLocationPtrP> cutPtrVec;
    vector<double> TheroMass1; // = pro->theoryMassC_Cz;

    ppm_value_ptr ppmValuePtr = std::make_shared<ppm_value>(prsmArg->min_ppm_value);

    if (prsmArg->BYIonsActions.count(msalignPtr->getActivationSub())) {
        TheroMass1 = proteinPtr->getTheoryMassCBy();
    }
    if (prsmArg->CZIonsActions.count(msalignPtr->getActivationSub())) {
        TheroMass1 = proteinPtr->getTheoryMassCCz();
    }

    int rightRange = TheroMass1.size();
    int n = 0;
    for (int ji = 0; ji < TheroMass1.size(); ++ji) { //控制C端截取
        if (TheroMass1.back() < minPrecursorMass) {      //如果理论最大值 < precursorMass - 500.无论如何添加ptmMass都无法到达
            break;
        }
        for (int xi = leftRange + 1; xi < TheroMass1.size(); ++xi) { //控制N端截取
            if (TheroMass1[xi] < minPrecursorMass) {                 //如果小于precursormass - 500 ，continue
                leftRange = xi;
                continue;
            }
            if (TheroMass1[xi] > maxPrecursorMass) { //如果大于precursormass + 500 , break ;
                rightRange = xi;
                break;
            }
            if (prsmArg->BYIonsActions.count(msalignPtr->getActivationSub())) {
                n = mylib::speculate_lib::txpb_binarey_search_ex(ptmMass, ptmMass.size(), precursorMass - TheroMass1[xi]);
            }
            if (prsmArg->CZIonsActions.count(msalignPtr->getActivationSub())) {
                n = mylib::speculate_lib::txpb_binarey_search_ex(ptmMass, ptmMass.size(), precursorMass - prsmArg->sub_etd - TheroMass1[xi]);
            }

            if (n < 0) {
                continue;
            }

            double tMass = 0; // TheroMass1[xi] + subMassETD + ptmMass[n] ;
            if (prsmArg->BYIonsActions.count(msalignPtr->getActivationSub())) {
                tMass = TheroMass1[xi] + ptmMass[n];
            }
            if (prsmArg->CZIonsActions.count(msalignPtr->getActivationSub())) {
                tMass = TheroMass1[xi] + prsmArg->sub_etd + ptmMass[n];
            }
            //进行三个H范围的判断
            if (ppmValuePtr->calculate_ppm_value(tMass, precursorMass)) {
                c::arg::cutLocationPtrP cp = new c::arg::cutLocation();
                cp->cutLocationC = ji;
                cp->cutLocationN = TheroMass1.size() - 1 - xi;
                cp->variablePtmsMasss = ptmMass[n];
                if(proteinPtr->getTheoryMassNBy().size()-(cp->cutLocationC+cp->cutLocationN)>7)
                    cutPtrVec.push_back(cp);
            }
        } //后端截断结束
        if (mapi.count(proteinPtr->getProteinSequence()[seqLog])) {
            subMass = mapi[proteinPtr->getProteinSequence()[seqLog--]];
        } else {
            seqLog--;
        }
        for (int ui = leftRange; ui < TheroMass1.size(); ++ui) { //截去前端
            TheroMass1[ui] -= subMass;
        }
    } //前段截断结束

    sort(cutPtrVec.begin(), cutPtrVec.end(), select_cut_sortV2);

    if (cutPtrVec.empty()) {
        return;
    } else {
        cutMap.insert(make_pair(cutPtrVec.front()->cutLocationN, cutPtrVec.front()->cutLocationC));
        int sumCut = cutPtrVec.front()->cutLocationN + cutPtrVec.front()->cutLocationC;
        for (int xi = 1; xi < cutPtrVec.size(); xi++) {
            if (sumCut == cutPtrVec[xi]->cutLocationC + cutPtrVec[xi]->cutLocationN) {
                cutMap.insert(make_pair(cutPtrVec[xi]->cutLocationN, cutPtrVec[xi]->cutLocationC));
            } else {
                break;
            }
        }
    }
    //释放空间
    for (int xi = 0; xi < cutPtrVec.size(); xi++) {
        delete cutPtrVec[xi];
    }
//    delete ppmValuePtr;
}


bool preProcessTempArg(prsm::config_argument_ptr &configArgumentPtr, prsm::temp_prsm_argument_ptr &tempPrsmArgumentPtr,
                       msalign_ptr &msalignPtr, protein_ptr &proteinPtr, prsm::modify_ptr modifyPtr,
                       int cutN, int cutC) {
    int initLog = 0, var_ptm_mass_index = -1;

    if (configArgumentPtr->BYIonsActions.count(msalignPtr->getActivationSub())) {
        tempPrsmArgumentPtr->InitArgument(cutN, cutC, proteinPtr->getTheoryMassNBy(), proteinPtr->getTheoryMassCBy());
        initLog = mylib::data_stream::init_theory_ions_mass_(tempPrsmArgumentPtr->theory_ions_mass_n,
                                                             tempPrsmArgumentPtr->cut_location_n,
                                                             tempPrsmArgumentPtr->cut_location_c, 0);
        initLog = mylib::data_stream::init_theory_ions_mass_(tempPrsmArgumentPtr->theory_ions_mass_c,
                                                             tempPrsmArgumentPtr->cut_location_c,
                                                             tempPrsmArgumentPtr->cut_location_n,
                                                             configArgumentPtr->h20);
        var_ptm_mass_index = mylib::data_stream::txpb_binarey_search_ex(modifyPtr->modifyMassSeq,
                                                                        modifyPtr->modifyMassSeq.size(),
                                                                        msalignPtr->getPrecursorMass() - tempPrsmArgumentPtr->theory_ions_mass_c.back());
    }
    if (configArgumentPtr->CZIonsActions.count(msalignPtr->getActivationSub())) {
        tempPrsmArgumentPtr->InitArgument(cutN, cutC, proteinPtr->getTheoryMassNCz(), proteinPtr->getTheoryMassCCz());
        initLog = mylib::data_stream::init_theory_ions_mass_(tempPrsmArgumentPtr->theory_ions_mass_n,
                                                             tempPrsmArgumentPtr->cut_location_n,
                                                             tempPrsmArgumentPtr->cut_location_c, 17.026);
        initLog = mylib::data_stream::init_theory_ions_mass_(tempPrsmArgumentPtr->theory_ions_mass_c,
                                                             tempPrsmArgumentPtr->cut_location_c,
                                                             tempPrsmArgumentPtr->cut_location_n,
                                                             configArgumentPtr->core);
        var_ptm_mass_index = mylib::data_stream::txpb_binarey_search_ex(modifyPtr->modifyMassSeq, modifyPtr->modifyMassSeq.size(),
                                                                        msalignPtr->getPrecursorMass() - configArgumentPtr->sub_etd - tempPrsmArgumentPtr->theory_ions_mass_c.back());
    }

    tempPrsmArgumentPtr->var_mod_mass = modifyPtr->modifyMassSeq[var_ptm_mass_index]; // 可变修饰的质量
    double theory_proteoform_mass = 0.0;                                              // 记录当下理论蛋白质变体质量
    if (configArgumentPtr->BYIonsActions.count(msalignPtr->getActivationSub())) {
        theory_proteoform_mass = tempPrsmArgumentPtr->theory_ions_mass_c.back() + modifyPtr->modifyMassSeq[var_ptm_mass_index];
    }
    if (configArgumentPtr->CZIonsActions.count(msalignPtr->getActivationSub())) {
        theory_proteoform_mass = tempPrsmArgumentPtr->theory_ions_mass_c.back()
                                 + modifyPtr->modifyMassSeq[var_ptm_mass_index] + configArgumentPtr->sub_etd;
    }


//    if (!mylib::calculate_lib::calculate_ppm_value_is_low(theory_proteoform_mass,
//                                                          msalignPtr->getPrecursorMass(),
//                                                          configArgumentPtr->H,
//                                                          configArgumentPtr->min_ppm_value)) {
//        return false;
//    }
    return true;
}

void prsm_varible_ptms::search_ptms_(vector<msalign_ptr> &msalign_ptrs,
                                     prsm::modify_ptr modifyPtr,
                                     prsm::config_argument_ptr configArgumentPtr,
                                     const string &outFileName) {

    //对每一条ms进行鉴定
    // ms 在候选蛋白质库中的每一套蛋白（n个）进行鉴定，先判断当前蛋白最佳的几种（k个）截断形式。时间复杂度为O(n*k)
    // 最后，选取当前蛋白质下打分最高的形。并选出候选库中，最优的蛋白质变体形式。
    vector<double> select_ptm_mass;
    int max_ptm_num = 3;
    for (const auto& entry : modifyPtr->mapi) {
        // entry.first 是 double, entry.second 是 string
        if (entry.second.length() <= max_ptm_num) {
            // 如果string的长度小于等于max_ptm_num，将对应的double添加到数组中
            select_ptm_mass.push_back(entry.first);
        }
    }


    for (auto msIt = msalign_ptrs.begin(); msIt != msalign_ptrs.end(); msIt++, ++cur_ms_num) {
        msalignPtrVec evalueMsVec;
        if ((*msIt)->getPrecursorMass() < configArgumentPtr->min_precursor_mass) {
            continue;
        }
        if (cur_ms_num % 100 == 0) {
            cout << "\rCurrent thread id : " << std::this_thread::get_id()
                 << " . number of PrSMs : " << cur_ms_num << " th ";
            fflush(stdout);
        }
        cout<<endl<<"var sp:"<<(*msIt)->getId()<<"-candi_protein_size:"<<(*msIt)->getCandidateProteins().size()<<endl;
        for (auto prIt = (*msIt)->candidate_proteins.begin(); prIt != (*msIt)->candidate_proteins.end(); ++prIt) {
//            cout<<endl<<" sp:"<<(*msIt)->getId()<<"  pro:"<<(*prIt)->getProteinTitle()<<" start"<<endl;
            map<int, int> cutMap; // 截断map
            prsm_varible_ptms::terminal_trunction((*prIt), (*msIt), configArgumentPtr, select_ptm_mass,
                                                  (*msIt)->getPrecursorMass(), cutMap); //探测最优的几种截断形式


            for (map<int, int>::iterator cutMapIt = cutMap.begin(); cutMapIt != cutMap.end(); ++cutMapIt) {
                prsm::temp_prsm_argument_ptr tempPrsmArgPtr = new prsm::temp_prsm_argument(); // 临时参数指针
                if (!preProcessTempArg(configArgumentPtr, tempPrsmArgPtr,
                                       (*msIt), (*prIt), modifyPtr, cutMapIt->first, cutMapIt->second)) {
                    continue;
                }
                this->prsm_process_var(modifyPtr, (*prIt), (*msIt), tempPrsmArgPtr,modifyPtr->modifyMassSeq, configArgumentPtr); // 鉴定
                ///保存xml格式中间结果
                bool out_for_evalue = true;
                if(out_for_evalue){
                    msalign_ptr curModMsalign = std::make_shared<msalign>(**msIt);
                    prsm_unknown_ptms::storagePrSMsMsg(*prIt, curModMsalign, tempPrsmArgPtr);
                    curModMsalign->setPrsmPtmsMass(tempPrsmArgPtr->var_mod_mass);
                    prsm::modify_ptr tmp = std::make_shared<prsm::Modify>(*modifyPtr);
                    curModMsalign->setPrsmUnknownModPtr(tmp);
                    if (configArgumentPtr->BYIonsActions.count(curModMsalign->getActivationSub())) {
                        curModMsalign->setPrsmTheoryProteoformMass(curModMsalign->getPrsmPtmsMass() + tempPrsmArgPtr->theory_ions_mass_c.back());
                    }
                    if (configArgumentPtr->CZIonsActions.count(curModMsalign->getActivationSub())) {
                        curModMsalign->setPrsmTheoryProteoformMass((curModMsalign->getPrsmPtmsMass() + tempPrsmArgPtr->theory_ions_mass_c.back()) + configArgumentPtr->sub_etd);
                    }
                    evalueMsVec.push_back(make_shared<msalign>(*curModMsalign));
//                    delete curModMsalign;
                }
                ///end

                if ((tempPrsmArgPtr->prsm_score - (*msIt)->getPrsmBestScore()) > 0.000001) {      //当前的蛋白质形鉴定完成，保存最优解
                    msalign_ptr msalignPtr = (*msIt);
                    prsm_unknown_ptms::storagePrSMsMsg((*prIt), (*msIt), tempPrsmArgPtr);
                    msalignPtr->setPrsmPtmsMass(tempPrsmArgPtr->var_mod_mass);
                    if (configArgumentPtr->BYIonsActions.count(msalignPtr->getActivationSub())) {
                        msalignPtr->setPrsmTheoryProteoformMass(msalignPtr->getPrsmPtmsMass() + tempPrsmArgPtr->theory_ions_mass_c.back());
                    }
                    if (configArgumentPtr->CZIonsActions.count(msalignPtr->getActivationSub())) {
                        msalignPtr->setPrsmTheoryProteoformMass((msalignPtr->getPrsmPtmsMass() + tempPrsmArgPtr->theory_ions_mass_c.back()) + configArgumentPtr->sub_etd);
                    }
                }
                //释放空间
                delete tempPrsmArgPtr;
            }
//            cout<<endl<<" sp:"<<(*msIt)->getId()<<"  pro:"<<(*prIt)->getProteinTitle()<<" end"<<endl;
        }

//        cout<<endl<<"sp:"<<(*msIt)->getId()<<"-candi_protein_size:"<<(*msIt)->getCandidateProteins().size()<<" "<<evalueMsVec.size()<<"   end"<<endl;
        //保存鉴定出的最佳蛋白质变体
//        this->getPrsmMsgWritePtr()->writeResultNews((*msIt), outFileName, modifyPtr);

        ///out_for_rfscore
        if(evalueMsVec.size() != 0){
            string xmlFileName = "wtop_prsm.test_combined";
            prsm_msg_write::OutXmlToEvalue(evalueMsVec, xmlFileName);
//            for(auto x:evalueMsVec){
//                delete x->getPrsmUnknownModPtr();
//                delete x;
//            }
            evalueMsVec.clear();
        }
    }
    select_ptm_mass.clear();
}

