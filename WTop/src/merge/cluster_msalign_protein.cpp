#include "cluster_msalign_protein.hpp"


bool spectSortFun1(c::arg::cutLocationPtrP a, c::arg::cutLocationPtrP b) {
    return fabs(a->variablePtmsMasss) < fabs(b->variablePtmsMasss);
}

bool spectSortFun2(c::arg::cutLocationPtrP a, c::arg::cutLocationPtrP b) {
    return a->cutLocationC+a->cutLocationN < b->cutLocationC+b->cutLocationN;
}

bool cmpCutPtms(cutMassPtrP a, cutMassPtrP b) {
    return a->SumCut < b->SumCut;
}

namespace prsm {
    /**
     * 由质谱质量过滤
     * @param msContainer 质谱序列
     * @param prContainer 蛋白序列
     * @param ptmMassSeq  已知修饰质量序列
     * @param CIDmslogMap CID过滤结果
     * @param ETDmslogMap ETD过滤结果
     * @param addMassCID
     * @param addMassETD
     */
    void sameConbine(const vector<msalign_processor_ptr> &msContainer, const vector<protein_processor_ptr> &prContainer,
                     const vector<double> &ptmMassSeq,
                     map<long, dataSetPtrP> &CIDmslogMap, map<long, dataSetPtr> &ETDmslogMap, double addMassCID,
                     double addMassETD) {
        double msSize = msContainer.size();
        double hasDeal = 0;
        vector<msalign_processor_ptr>::const_iterator it = msContainer.begin();
        while (it != msContainer.end()) {       //
            hasDeal += 1;
            cout << "\rMsDataSet has deal   " << hasDeal / msSize * 100 << " %\n";
            fflush(stdout);
            if ((*it)->ions_mass_container.size() < 20) {
                it++;
                continue;
            }
            long mass = (long) ((*it)->precursor_mass_double * 100); // 小数点往后移二位
            if (!CIDmslogMap.count(mass)) {       // 如果该质量无
                dataSetPtrP p = new dataSet();    //创建数据集指针
                CIDmslogMap[mass] = p;                        //添加进指针
                double prItPtrSize = prContainer.size();        //遍历蛋白质数据库
                double prHasDeal = 0;          //记录进度
                vector<protein_processor_ptr>::const_iterator prItPtr = prContainer.begin();
                for ( ; prItPtr != prContainer.end(); prItPtr++) {
                    prHasDeal += 1;
                    cout << "\rMs Protein has deal  " << prHasDeal / prItPtrSize * 100 << " % ";
                    fflush(stdout);
                    if ((*it)->activation.find("HCD") != (*it)->activation.npos ||
                        (*it)->activation.find("CID") != (*it)->activation.npos) {
                        spectProteinFormReWriteCID((*prItPtr), ptmMassSeq, (*it)->precursor_mass_double, p, addMassCID);
                    } else if ((*it)->activation.find("ETD") != (*it)->activation.npos) {
                        spectProteinFormReWriteETD((*prItPtr), ptmMassSeq, (*it)->precursor_mass_double, p, addMassETD);
                    }
                }
                sort(p->cutPtms.begin(), p->cutPtms.end(), cmpCutPtms);
            }
            it++;
        }
        // cout<<"Cid = "<<CIDmslogMap.size()<<endl;
        // cout<<"Etd = "<<ETDmslogMap.size()<<endl;
        // for (map<long,dataSetPtr>::iterator it = CIDmslogMap.begin(); it != CIDmslogMap.end(); ++it) {
        //    cout<<"size = "<<it->second->protein.size()<<endl;
        // }
        // cout<<mslogMap.size()<<endl;
    }

    /**
     * subMass = H2o-1.9919 ;
     * 最终Adjust mass = 理论最大值 - 截断质量 + PtmMass + subMass = PrecursorMass ;
     * PtmMass = PrecursorMass - 理论最大值 + 截断质量 - subMass ;
     * 从C端开始截
     * @param pro
     * @param ptmMass
     * @param precursorMass
     * @param p
     * @param addMass
     */
    void spectProteinFormReWriteCID(protein_processor_ptr pro, const vector<double> &ptmMass,
                                    const double &precursorMass, dataSetPtrP p, double addMass) {

        if (pro->theoryMassC_By.empty() || ptmMass.empty()) {    //条件判断
            return;
        }

        vector<double> TheroMass1 = pro->theoryMassC_By; //理论质量由C端产生
        double minPrecursorMass = precursorMass - 500; //产生范围
        double maxPrecursorMass = precursorMass + 500;
        int leftRange = 0;
        int rightRange = TheroMass1.size();
        int seqLog = pro->protein_sequence.length() - 1;  //C端
        double subMass = 0;
        vector<c::arg::cutLocationPtrP> cutPtrVec; //用于接受该条蛋白中，所有可能的截断情况,最终,取截断质量相等的蛋白质形
        //产生截断向量序列
        for (int ji = 0; ji < TheroMass1.size()/2; ++ji) {   //控制前截,N端截断的个数
            if (TheroMass1.back() < minPrecursorMass) { //如果理论最大值<precursormass - 500,无论如何添加PTMS质量都无法到达
                break;
            }
            for (int xi = leftRange + 1; xi < TheroMass1.size(); ++xi) { //xi控制后截
                if (TheroMass1[xi] < minPrecursorMass) {    //如果小于precursormass - 500,继续
                    leftRange = xi;
                    continue;
                }
                if (TheroMass1[xi] > maxPrecursorMass) { //如果大于precursormass + 500,break ;
                    rightRange = xi;
                    break;
                }
                //查找当前C端截断情况下是否有合适的PTMs的偏移
                int n = mylib::speculate_lib::txpb_binarey_search_ex(ptmMass, ptmMass.size(),
                                                        precursorMass - TheroMass1[xi]);
                if (n < 0 ) {
                    continue;
                }
                double ppm = (TheroMass1[xi] + ptmMass[n] - precursorMass) / (TheroMass1[xi] + ptmMass[n]) * 1000000;
                double ppm_DeH = (TheroMass1[xi] + ptmMass[n]-1.007276 - precursorMass) / (TheroMass1[xi] + ptmMass[n]-1.007276) * 1000000;
                double ppm_AeH = (TheroMass1[xi] + ptmMass[n]+1.007276 - precursorMass) / (TheroMass1[xi] + ptmMass[n]+1.007276) * 1000000;
                if (fabs(ppm) < 15 || fabs(ppm_DeH) < 15 || fabs(ppm_AeH) < 15 ) {   //如果ppm<15
                    c::arg::cutLocationPtrP cp = new c::arg::cutLocation();
                    cp->cutLocationC = ji;
                    cp->cutLocationN = TheroMass1.size() - 1 - xi;
                    // cp->ppm = ppm ;
                    // cp->precursorMass = precursorMass ;
                    // // cp->ptmsSeq = it->second ;
                    // cp->theroMass = TheroMass1[xi];
                    cp->variablePtmsMasss = ptmMass[n];
                    cutPtrVec.push_back(cp);
                }
            } //后端截断结束
            if (mapi.count(pro->protein_sequence[seqLog])) {
                subMass = mapi[pro->protein_sequence[seqLog--]];
            } else {
                seqLog--;
            }
            for (int ui = leftRange; ui < TheroMass1.size(); ++ui) {   //截去前端质量
                TheroMass1[ui] -= subMass;
            }
        }//前段截断结束
        sort(cutPtrVec.begin(), cutPtrVec.end(), spectSortFun2);
        if (cutPtrVec.empty()) {
            return;
        } else {
            cutMassPtrP t = new cutMass();
            t->Ncut = cutPtrVec.front()->cutLocationN;
            t->Ccut = cutPtrVec.front()->cutLocationC;
            t->SumCut = cutPtrVec.front()->cutLocationN + cutPtrVec.front()->cutLocationC;
            t->log = p->prNumber++;
            p->cutPtms.push_back(t);       //记录截断信息
            p->protein.push_back(pro);      //记录蛋白信息
            int sumCut = t->SumCut;
            for (int xi = 1; xi < cutPtrVec.size(); xi++) {
                if (sumCut == cutPtrVec[xi]->cutLocationC + cutPtrVec[xi]->cutLocationN) {
                    cutMassPtrP tp = new cutMass();
                    tp->Ncut = cutPtrVec[xi]->cutLocationN;
                    tp->Ccut = cutPtrVec[xi]->cutLocationC;
                    tp->SumCut = cutPtrVec[xi]->cutLocationC + cutPtrVec[xi]->cutLocationN;
                    tp->log = p->prNumber++;
                    p->cutPtms.push_back(tp);
                    p->protein.push_back(pro);
                } else {
                    break ;
                }
            }
        }
//        cout << "可能的蛋白质形式个数为 : " << cutPtrVec.size() << endl;
//        for (int xi = 0; xi < cutPtrVec.size(); ++xi) {
//            cout << "n端截断为 : " << cutPtrVec[xi]->cutLocationN << " c端截断为 : " << cutPtrVec[xi]->cutLocationC << endl;
//            cout << "ppm = " << cutPtrVec[xi]->ppm << endl;
//            cout << "Ptms Mass = " << cutPtrVec[xi]->variablePtmsMasss << endl;
//            cout << "ProteinformMass = " << cutPtrVec[xi]->theroMass + cutPtrVec[xi]->variablePtmsMasss << endl;
//            cout << "PrecursorMass = " << cutPtrVec[xi]->precursorMass << endl;
//            cout << "#END\n\n";
//        }
        for (int xi = 0; xi < cutPtrVec.size(); xi++) {     //释放空间
            delete cutPtrVec[xi];
        }
    }

    /**
     * subMass = H2o-1.9919 ;
     * 最终Adjust mass = 理论最大值 - 截断质量 + PtmMass + subMass = PrecursorMass ;
     * PtmMass = PrecursorMass - 理论最大值 + 截断质量 - subMass ;
     * 从C端开始截
     * @param pro
     * @param ptmMass
     * @param precursorMass
     * @param p
     * @param addMass
     */
    void spectProteinFormReWriteETD(protein_processor_ptr pro, const vector<double> &ptmMass,
                                    const double &precursorMass, dataSetPtrP p, double addMass) {
        if (pro->theoryMassC_Cz.empty() || ptmMass.empty()) {    //条件判断
            EX_TRACE("TheroMass or modifyTable is null ~\n");
            return;
        }
        double subMassETD = 18.01056 - 1.9919;      //
        vector<double> TheroMass1 = pro->theoryMassC_Cz;
        double minPrecursorMass = precursorMass - 500; //范围
        double maxPrecursorMass = precursorMass + 500;
        int leftRange = 0;
        int rightRange = TheroMass1.size();
        int seqLog = pro->protein_sequence.length() - 1;
        double subMass = 0;
        double hMass = 1.007276 ;
        double ppm = 15 ;
        vector<c::arg::cutLocationPtrP> cutPtrVec;
        int n = 0 ;
        for (int ji = 0; ji < TheroMass1.size()/2; ++ji) {   //控制C端截取
            if (TheroMass1.back() < minPrecursorMass) {      //如果理论最大值 < precursorMass - 500.无论如何添加ptmMass都无法到达
                break;
            }
            for (int xi = leftRange + 1; xi < TheroMass1.size(); ++xi) {      //控制N端截取
                if (TheroMass1[xi] < minPrecursorMass) {    //如果小于precursormass - 500 ，continue
                    leftRange = xi;
                    continue;
                }
                if (TheroMass1[xi] > maxPrecursorMass) { //如果大于precursormass + 500 , break ;
                    rightRange = xi;
                    break;
                }

                n = mylib::speculate_lib::txpb_binarey_search_ex(ptmMass, ptmMass.size(),
                                                    precursorMass - subMassETD - TheroMass1[xi]);
                if (n < 0) {
                    continue;
                }

//                cout << " N cut = " << TheroMass1.size() - 1 - xi
//                << " C cut = " << ji
//                << " ptmMass = " << ptmMass[n]
//                << endl ;

                //进行三个H范围的判断
                double tMass = TheroMass1[xi] + subMassETD + ptmMass[n] ;
                if (mylib::calculate_lib::calculate_ppm_value_is_low(tMass, precursorMass, hMass, ppm)) {
//                    cout << " N cut = " << TheroMass1.size() - 1 - xi
//                    << " C cut = " << ji
//                    << " ptmMass = " << ptmMass[n]
//                    << " tMass = " << tMass
//                    << " pMass = " << precursorMass
//                    << endl ;
                    c::arg::cutLocationPtrP cp = new c::arg::cutLocation();
                    cp->cutLocationC = ji;
                    cp->cutLocationN = TheroMass1.size() - 1 - xi;
                    // cp->ppm = ppm ;
                    // cp->precursorMass = precursorMass ;
                    // cp->ptmsSeq = it->second ;
                    // cp->theroMass = TheroMass1[xi];
                    cp->variablePtmsMasss = ptmMass[n];
                    cutPtrVec.push_back(cp);
                    // }
                }
            }  //后端截断结束
            if (mapi.count(pro->protein_sequence[seqLog])) {
                subMass = mapi[pro->protein_sequence[seqLog--]];
            } else {
                seqLog--;
            }
            for (int ui = leftRange; ui < TheroMass1.size(); ++ui) {   //截去前端
                TheroMass1[ui] -= subMass;
            }
        }//前段截断结束
        //截断越少认为该蛋白质形最优
        sort(cutPtrVec.begin(), cutPtrVec.end(), spectSortFun2);
        if (cutPtrVec.empty()) {
            return;
        } else {
            cutMassPtrP t = new cutMass();
            t->Ncut = cutPtrVec.front()->cutLocationN;
            t->Ccut = cutPtrVec.front()->cutLocationC;
            t->SumCut = cutPtrVec.front()->cutLocationN + cutPtrVec.front()->cutLocationC;
            t->log = p->prNumber++;
            p->cutPtms.push_back(t);       //记录截断信息
            p->protein.push_back(pro);      //记录蛋白信息
            int sumCut = t->SumCut;
            //
            for (int xi = 1; xi < cutPtrVec.size(); xi++) {
                if (sumCut == cutPtrVec[xi]->cutLocationC + cutPtrVec[xi]->cutLocationN) {
                    cutMassPtrP tp = new cutMass();
                    tp->Ncut = cutPtrVec[xi]->cutLocationN;
                    tp->Ccut = cutPtrVec[xi]->cutLocationC;
                    tp->SumCut = cutPtrVec[xi]->cutLocationC + cutPtrVec[xi]->cutLocationN;
                    tp->log = p->prNumber++;
                    p->cutPtms.push_back(tp);
                    p->protein.push_back(pro);
                } else {
                    break ;
                }
            }
        }
//        cout << "可能的蛋白质形式个数为 : " << cutPtrVec.size() << endl;
//        for (int xi = 0; xi < cutPtrVec.size(); ++xi) {
//            cout << "n端截断为 : " << cutPtrVec[xi]->cutLocationN << " c端截断为 : " << cutPtrVec[xi]->cutLocationC << endl;
//            cout << "ppm = " << cutPtrVec[xi]->ppm << endl;
//            cout << "Ptms Mass = " << cutPtrVec[xi]->variablePtmsMasss << endl;
//            cout << "ProteinformMass = " << cutPtrVec[xi]->theroMass + cutPtrVec[xi]->variablePtmsMasss << endl;
//            cout << "PrecursorMass = " << cutPtrVec[xi]->precursorMass << endl;
//            cout << "#END\n\n";
//        }
        //释放空间
        for (int xi = 0; xi < cutPtrVec.size(); xi++) {
            delete cutPtrVec[xi];
        }
    }

    // 2023-1-8 合并By Cz
    void DetermineProteinForm(protein_processor_ptr pro, msalign_processor_ptr sp, config_argument_ptr prsmArg, const vector<double> &ptmMass, const double &precursorMass,
                              map<int,int> &cutMap)
    {
        if (pro->theoryMassC_By.empty() || pro->theoryMassC_Cz.empty() || ptmMass.empty()) {    //条件判断
            return;
        }
        double minPrecursorMass = precursorMass - prsmArg->ptmMassMax;
        double maxPrecursorMass = precursorMass + prsmArg->ptmMassMax;
        int leftRange = 0;
        int seqLog = pro->protein_sequence.length() - 1;
        double subMass = 0;
        double hMass = prsmArg->H ;
        double ppm = prsmArg->min_ppm_value ;
        vector<c::arg::cutLocationPtrP> cutPtrVec;
        vector<double> TheroMass1 ;// = pro->theoryMassC_Cz;
        if (prsmArg->BYIonsActions.count(sp->activation_sub)) {
            TheroMass1 = pro->theoryMassC_By ;
        }
        if (prsmArg->CZIonsActions.count(sp->activation_sub)) {
            TheroMass1 = pro->theoryMassC_Cz ;
        }
        int rightRange = TheroMass1.size();

        int n = 0 ;
        for (int ji = 0; ji < TheroMass1.size() / 2; ++ji) {   //控制C端截取
            if (TheroMass1.back() < minPrecursorMass) {      //如果理论最大值 < precursorMass - 500.无论如何添加ptmMass都无法到达
                break;
            }
            for (int xi = leftRange + 1; xi < TheroMass1.size(); ++xi) {      //控制N端截取
                if (TheroMass1[xi] < minPrecursorMass) {    //如果小于precursormass - 500 ，continue
                    leftRange = xi;
                    continue;
                }
                if (TheroMass1[xi] > maxPrecursorMass) { //如果大于precursormass + 500 , break ;
                    rightRange = xi;
                    break;
                }
                if (prsmArg->BYIonsActions.count(sp->activation_sub)) {
                    n = mylib::speculate_lib::txpb_binarey_search_ex(ptmMass, ptmMass.size(), precursorMass - TheroMass1[xi]);                }
                if (prsmArg->CZIonsActions.count(sp->activation_sub)) {
                    n = mylib::speculate_lib::txpb_binarey_search_ex(ptmMass, ptmMass.size(), precursorMass - prsmArg->sub_etd - TheroMass1[xi]);
                }

                if (n < 0) {
                    continue;
                }

                double tMass = 0 ; // TheroMass1[xi] + subMassETD + ptmMass[n] ;
                if (prsmArg->BYIonsActions.count(sp->activation_sub)) {
                    tMass = TheroMass1[xi] + ptmMass[n] ;
                }
                if (prsmArg->CZIonsActions.count(sp->activation_sub)) {
                    tMass = TheroMass1[xi] + prsmArg->sub_etd + ptmMass[n] ;
                }
                //进行三个H范围的判断
                if (mylib::calculate_lib::calculate_ppm_value_is_low(tMass, precursorMass, hMass, ppm)) {
                    c::arg::cutLocationPtrP cp = new c::arg::cutLocation();
                    cp->cutLocationC = ji;
                    cp->cutLocationN = TheroMass1.size() - 1 - xi;
                    cp->variablePtmsMasss = ptmMass[n];
                    cutPtrVec.push_back(cp);
                }
            }  //后端截断结束
            if (mapi.count(pro->protein_sequence[seqLog])) {
                subMass = mapi[pro->protein_sequence[seqLog--]];
            } else {
                seqLog--;
            }
            for (int ui = leftRange; ui < TheroMass1.size(); ++ui) {   //截去前端
                TheroMass1[ui] -= subMass;
            }
        }  //前段截断结束
        sort(cutPtrVec.begin(), cutPtrVec.end(), spectSortFun2);
        if (cutPtrVec.empty()) {
            return ;
        } else {
            cutMap.insert(make_pair(cutPtrVec.front()->cutLocationN,cutPtrVec.front()->cutLocationC));
            int sumCut = cutPtrVec.front()->cutLocationN + cutPtrVec.front()->cutLocationC;
            for (int xi = 1; xi < cutPtrVec.size(); xi++)
            {
                if (sumCut == cutPtrVec[xi]->cutLocationC + cutPtrVec[xi]->cutLocationN)
                {
                    cutMap.insert(make_pair(cutPtrVec[xi]->cutLocationN,cutPtrVec[xi]->cutLocationC));
                } else {
                    break ;
                }
            }
        }
        //释放空间
        for (int xi = 0; xi < cutPtrVec.size(); xi++) {
            delete cutPtrVec[xi];
        }
    }

    /**
     * 20221114 推测合适的蛋白质形
     * @param pro
     * @param ptmMass
     * @param precursorMass
     * @param prsmArg
     * @param addMass
     */
    void judgeProteinFormCID(protein_processor_ptr pro, const vector<double> &ptmMass, const double &precursorMass, config_argument_ptr &prsmArg,
                             map<int,int> &cutMap)
    {
        if (pro->theoryMassC_By.empty() || ptmMass.empty()) {    //条件判断
            return;
        }
        vector<double> TheroMass1 = pro->theoryMassC_By; //理论质量由C端产生
        double minPrecursorMass = precursorMass - prsmArg->ptmMassMax; //产生范围
        double maxPrecursorMass = precursorMass + prsmArg->ptmMassMax;
        int leftRange = 0;
        int rightRange = TheroMass1.size();
        int seqLog = pro->protein_sequence.length() - 1;  //C端
        double subMass = 0;
        vector<c::arg::cutLocationPtrP> cutPtrVec; //用于接受该条蛋白中，所有可能的截断情况,最终,取截断质量相等的蛋白质形
        //产生截断向量序列
        for (int ji = 0; ji < TheroMass1.size()/2; ++ji) {   //控制前截,N端截断的个数
            if (TheroMass1.back() < minPrecursorMass) { //如果理论最大值<precursormass - 500,无论如何添加PTMS质量都无法到达
                break;
            }
            for (int xi = leftRange + 1; xi < TheroMass1.size(); ++xi) { //xi控制后截
                if (TheroMass1[xi] < minPrecursorMass) {    //如果小于precursormass - 500,继续
                    leftRange = xi;
                    continue;
                }
                if (TheroMass1[xi] > maxPrecursorMass) { //如果大于precursormass + 500,break ;
                    rightRange = xi;
                    break;
                }
                //查找当前C端截断情况下是否有合适的PTMs的偏移
                int n = mylib::speculate_lib::txpb_binarey_search_ex(ptmMass, ptmMass.size(), precursorMass - TheroMass1[xi]);
                if (n < 0 ) {
                    continue;
                }
                double ppm = (TheroMass1[xi] + ptmMass[n] - precursorMass) / (TheroMass1[xi] + ptmMass[n]) * 1000000;
                double ppm_DeH = (TheroMass1[xi] + ptmMass[n] - prsmArg->H - precursorMass) / (TheroMass1[xi] + ptmMass[n]-prsmArg->H) * 1000000;
                double ppm_AeH = (TheroMass1[xi] + ptmMass[n] + prsmArg->H - precursorMass) / (TheroMass1[xi] + ptmMass[n]+prsmArg->H) * 1000000;
                if (fabs(ppm) < prsmArg->min_ppm_value || fabs(ppm_DeH) < prsmArg->min_ppm_value || fabs(ppm_AeH) < prsmArg->min_ppm_value )
                {   //如果ppm<15
                    c::arg::cutLocationPtrP cp = new c::arg::cutLocation();
                    cp->cutLocationC = ji;
                    cp->cutLocationN = TheroMass1.size() - 1 - xi;
                    // cp->ppm = ppm ;
                    // cp->precursorMass = precursorMass ;
                    // // cp->ptmsSeq = it->second ;
                    // cp->theroMass = TheroMass1[xi];
                    cp->variablePtmsMasss = ptmMass[n];
                    cutPtrVec.push_back(cp);
                }
            } //后端截断结束
            if (mapi.count(pro->protein_sequence[seqLog])) {
                subMass = mapi[pro->protein_sequence[seqLog--]];
            } else {
                seqLog--;
            }
            for (int ui = leftRange; ui < TheroMass1.size(); ++ui) {   //截去前端质量
                TheroMass1[ui] -= subMass;
            }
        }
        // 截断 排序
        sort(cutPtrVec.begin(), cutPtrVec.end(), spectSortFun2);
//        map<int,int> cutMap ; //key - N截  value - C截
        if (cutPtrVec.empty())
        {
            return ;
        } else {
            cutMap.insert(make_pair(cutPtrVec.front()->cutLocationN,cutPtrVec.front()->cutLocationC));
            int sumCut = cutPtrVec.front()->cutLocationN + cutPtrVec.front()->cutLocationC;
            for (int xi = 1; xi < cutPtrVec.size(); xi++)
            {
                if (sumCut == cutPtrVec[xi]->cutLocationC + cutPtrVec[xi]->cutLocationN)
                {
                    cutMap.insert(make_pair(cutPtrVec[xi]->cutLocationN,cutPtrVec[xi]->cutLocationC));
                } else {
                    break ;
                }
            }
        }
//        cout << "可能的蛋白质形式个数为 : " << cutPtrVec.size() << endl;
//        for (int xi = 0; xi < cutPtrVec.size(); ++xi) {
//            cout << "n端截断为 : " << cutPtrVec[xi]->cutLocationN << " c端截断为 : " << cutPtrVec[xi]->cutLocationC << endl;
//            cout << "ppm = " << cutPtrVec[xi]->ppm << endl;
//            cout << "Ptms Mass = " << cutPtrVec[xi]->variablePtmsMasss << endl;
//            cout << "ProteinformMass = " << cutPtrVec[xi]->theroMass + cutPtrVec[xi]->variablePtmsMasss << endl;
//            cout << "PrecursorMass = " << cutPtrVec[xi]->precursorMass << endl;
//            cout << "#END\n\n";
//        }
        for (int xi = 0; xi < cutPtrVec.size(); xi++) {     //释放空间
            delete cutPtrVec[xi];
        }
    }


    void judgeProteinFormETD(protein_processor_ptr pro, const vector<double> &ptmMass, const double &precursorMass, config_argument_ptr &prsmArg,
                             map<int,int> &cutMap)
     {
        if (pro->theoryMassC_Cz.empty() || ptmMass.empty()) {    //条件判断
            return;
        }
//        double subMassETD = 18.01056 - 1.9919;      //
        double subMassETD = prsmArg->h20-prsmArg->core;
        vector<double> TheroMass1 = pro->theoryMassC_Cz;
        double minPrecursorMass = precursorMass - prsmArg->ptmMassMax; //范围
        double maxPrecursorMass = precursorMass + prsmArg->ptmMassMax;
        int leftRange = 0;
        int rightRange = TheroMass1.size();
        int seqLog = pro->protein_sequence.length() - 1;
        double subMass = 0;
        double hMass = prsmArg->H ;
        double ppm = prsmArg->min_ppm_value ;
        vector<c::arg::cutLocationPtrP> cutPtrVec;
        int n = 0 ;
        for (int ji = 0; ji < TheroMass1.size()/2; ++ji) {   //控制C端截取
            if (TheroMass1.back() < minPrecursorMass) {      //如果理论最大值 < precursorMass - 500.无论如何添加ptmMass都无法到达
                break;
            }
            for (int xi = leftRange + 1; xi < TheroMass1.size(); ++xi) {      //控制N端截取
                if (TheroMass1[xi] < minPrecursorMass) {    //如果小于precursormass - 500 ，continue
                    leftRange = xi;
                    continue;
                }
                if (TheroMass1[xi] > maxPrecursorMass) { //如果大于precursormass + 500 , break ;
                    rightRange = xi;
                    break;
                }
                n = mylib::speculate_lib::txpb_binarey_search_ex(ptmMass, ptmMass.size(),
                                                    precursorMass - subMassETD - TheroMass1[xi]);
                if (n < 0) {
                    continue;
                }

//                cout << " N cut = " << TheroMass1.size() - 1 - xi
//                << " C cut = " << ji
//                << " ptmMass = " << ptmMass[n]
//                << endl ;

                //进行三个H范围的判断
                double tMass = TheroMass1[xi] + subMassETD + ptmMass[n] ;
                if (mylib::calculate_lib::calculate_ppm_value_is_low(tMass, precursorMass, hMass, ppm)) {
//                    cout << " N cut = " << TheroMass1.size() - 1 - xi
//                    << " C cut = " << ji
//                    << " ptmMass = " << ptmMass[n]
//                    << " tMass = " << tMass
//                    << " pMass = " << precursorMass
//                    << endl ;
                    c::arg::cutLocationPtrP cp = new c::arg::cutLocation();
                    cp->cutLocationC = ji;
                    cp->cutLocationN = TheroMass1.size() - 1 - xi;
                    // cp->ppm = ppm ;
                    // cp->precursorMass = precursorMass ;
                    // cp->ptmsSeq = it->second ;
                    // cp->theroMass = TheroMass1[xi];
                    cp->variablePtmsMasss = ptmMass[n];
                    cutPtrVec.push_back(cp);
                    // }
                }
            }  //后端截断结束
            if (mapi.count(pro->protein_sequence[seqLog])) {
                subMass = mapi[pro->protein_sequence[seqLog--]];
            } else {
                seqLog--;
            }
            for (int ui = leftRange; ui < TheroMass1.size(); ++ui) {   //截去前端
                TheroMass1[ui] -= subMass;
            }
        }//前段截断结束
        //截断越少认为该蛋白质形最优
        sort(cutPtrVec.begin(), cutPtrVec.end(), spectSortFun2);
        if (cutPtrVec.empty())
        {
            return ;
        } else {
            cutMap.insert(make_pair(cutPtrVec.front()->cutLocationN,cutPtrVec.front()->cutLocationC));
            int sumCut = cutPtrVec.front()->cutLocationN + cutPtrVec.front()->cutLocationC;
            for (int xi = 1; xi < cutPtrVec.size(); xi++)
            {
                if (sumCut == cutPtrVec[xi]->cutLocationC + cutPtrVec[xi]->cutLocationN)
                {
                    cutMap.insert(make_pair(cutPtrVec[xi]->cutLocationN,cutPtrVec[xi]->cutLocationC));
                } else {
                    break ;
                }
            }
        }
//        cout << "可能的蛋白质形式个数为 : " << cutPtrVec.size() << endl;
//        for (int xi = 0; xi < cutPtrVec.size(); ++xi) {
//            cout << "n端截断为 : " << cutPtrVec[xi]->cutLocationN << " c端截断为 : " << cutPtrVec[xi]->cutLocationC << endl;
//            cout << "ppm = " << cutPtrVec[xi]->ppm << endl;
//            cout << "Ptms Mass = " << cutPtrVec[xi]->variablePtmsMasss << endl;
//            cout << "ProteinformMass = " << cutPtrVec[xi]->theroMass + cutPtrVec[xi]->variablePtmsMasss << endl;
//            cout << "PrecursorMass = " << cutPtrVec[xi]->precursorMass << endl;
//            cout << "#END\n\n";
//        }
        //释放空间
        for (int xi = 0; xi < cutPtrVec.size(); xi++) {
            delete cutPtrVec[xi];
        }
    }

}