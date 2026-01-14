//
// Created by 碎雨粘霓裳 on 2022/4/29.
//

#include "unknown_ptm_search.h"

namespace prsm {

    /**
     * 将结果写入文件
     * @param msIt 质谱
     * @param Result_file 结果文件名
     * @param modptr 修饰 (用于解析质量偏移)
     */
    void writeResult(msalign_processor_ptr msIt, string Result_file, modify_ptr modptr) {
        //输出匹配结果入txt文件
        if (msIt->match_protein_name.size() > 0) {
            std::ofstream out(Result_file, ios::out | ios::app);
            // std::ofstream out_pro("pro_title.txt",ios::out|ios::app);
            // out_pro<<Ms->BestSeqTitle<<endl;
            double t_mass = atof(
                    msIt->precursor_mass.substr(msIt->precursor_mass.find('=') + 1, msIt->precursor_mass.size()).c_str());
            out << "#BEGIN" << endl;
            out << "#" << msIt->Id
                << " , #" << msIt->scans
                << ", #Types : " << msIt->activation.substr(msIt->activation.find('=') + 1, msIt->activation.size())
                << ", #Protero : " << msIt->match_protein_name << "." << endl;
            out << "Match Peaks = " << msIt->best_peak
                << ",score = " << msIt->best_score << endl;
            out << "CutLocationN = " << msIt->cut_c_n
                << " ,CutLocationC = " << msIt->cut_n_c << endl;
            out << "Precursor mass = " << setiosflags(ios::fixed) << setprecision(4)
                << msIt->precursor_mass.substr(msIt->precursor_mass.find('=') + 1, msIt->precursor_mass.size())
                << " ,Proteoform mass = " << msIt->theory_proteoform_mass << endl;
            out << "PPM = " << (t_mass - msIt->theory_proteoform_mass) / t_mass * 1000000 << endl;

            //处理seq
            string be;
            if (msIt->cut_c_n != 0) {
                msIt->match_protein_seq.insert(msIt->cut_c_n, ".");
            }
            if (msIt->cut_n_c != 0) {
                msIt->match_protein_seq.insert(msIt->match_protein_seq.size() - msIt->cut_n_c, ".");
            }
            out << msIt->match_protein_seq << endl;
            out << "#Mod BEGIN" << endl;
            double pMass = 0;
            for (int j = 0; j < msIt->best_ptm.size(); j++) {       //输出修饰
                string ptms = modptr->analysis(msIt->best_ptm[j].mod_mass);
                if (ptms != "0" && msIt->modify_mass != 0) {
                    out << "start: " << (int) msIt->best_ptm[j].first << "-" << (int) msIt->best_ptm[j].second
                        << ";mass: " << msIt->best_ptm[j].mod_mass << " ; Mod# " << ptms << endl;
                    pMass += msIt->best_ptm[j].mod_mass;
                }
            }
            out << "#Mod END" << endl;
            out << "Spact ModifyMass = " << msIt->modify_mass << "\nTrue ModifyMass = " << pMass << endl;
            out << "All peaks = " << msIt->ions_mass_container.size() << endl;
            out << "match peaks = " << (int) msIt->best_peak << endl;
            out << "Ms complementIons =  " << msIt->compIonsIndexMap.size() << endl;
            out << "matched fragment ions = " << msIt->matchFragmentIonsSize << endl;
            out << "mut_x = " << msIt->mut_s << endl;
            out << "ionsPair = " << msIt->bestIonsPairNub << endl;
            out << "huBuIonsCount = " << (int) msIt->huBuCount - 1 << endl;
            out << "seqIonsCount = " << msIt->seqScoreCount << endl;
            out << "seqScore = " << msIt->seqScore << endl;
            out << "#END\n" << endl;
            out.close();
        }   //输出答案结束
    }

    //将结果写入文件，仅供测试数据使用
    void WriteSequenceAlignment(const char *output,
                                const std::vector<double> &reference_orig, const std::vector<double> &peer_orig,
                                vector<pair<long, long> > &alignment, vector<double> &sub,
                                vector<double> &index, vector<double> &ppm,
                                vector<double> &RefeZscore, vector<double> &PeerZscore) {
        vector<std::string> tmp_rec;
        std::ostringstream title;
        title << setw(10) << "Thero ID" << setw(10) << " Mono ID" << "  | " << setw(15) << " Thero Mass" << setw(15)
              << "Mono Mass" << "  | " << setw(10) << "Mono-Thero" << " ｜ " << setw(10) << "Index" << " ｜ " << setw(10)
              << "Ppm" << " ｜ ";
        tmp_rec.push_back(title.str());
        for (long i = 0; i < alignment.size(); i++) {
            //----- output to string ----//
            std::ostringstream o;
            o << setw(10) << alignment[i].first << " " << setw(10) << alignment[i].second << " | ";
            o << setw(15) << setiosflags(ios::fixed) << setprecision(4) << reference_orig[alignment[i].first] << " "
              << setw(15) << setiosflags(ios::fixed) << setprecision(4) << peer_orig[alignment[i].second] << " | ";
            o << setw(10) << setiosflags(ios::fixed) << setprecision(4) << sub[i] << " ｜ ";
            if (index[i] != INT_MAX) {
                if (fabs(ppm[i]) < 15) {
                    // o<<setw(10)<<alignment[i].first<<" "<<setw(10)<<alignment[i].second<<" | ";
                    // o<<setw(15)<<setiosflags(ios::fixed)<<setprecision(4)<<reference_orig[alignment[i].first+start]<<" "<<setw(15)<<setiosflags(ios::fixed)<<setprecision(4)<<peer_orig[alignment[i].second]<<" | ";
                    o << setw(10) << setiosflags(ios::fixed) << setprecision(4) << index[i] << " ｜ ";
                    o << setw(10) << setiosflags(ios::fixed) << setprecision(4) << ppm[i] << " ｜ ";
                    o << setw(10) << setiosflags(ios::fixed) << setprecision(4) << RefeZscore[alignment[i].first]
                      << " ｜ ";
                    o << setw(10) << setiosflags(ios::fixed) << setprecision(4) << PeerZscore[alignment[i].second]
                      << " ｜ ";
                    // std::string s=o.str();
                    // tmp_rec.push_back(s);
                }
            }
            //----- record string -----//
            std::string s = o.str();
            tmp_rec.push_back(s);
        }
        //----- output to file ------//
        FILE *fp = fopen(output, "wb");
        for (long i = 0; i < (long) tmp_rec.size(); i++) fprintf(fp, "%s\n", tmp_rec[i].c_str());
        fclose(fp);
    }

    void test_out(char pI, prsm::temp_prsm_argument_ptr arg, prsm::msalign_processor_ptr spact) {
        //            if(arg->varible_mod_mass > (double)(500/2)) {
        //                int count = 0 ;
        //                count += g::proc::find_error_match_ions(modifier_table,C_ions_c_n);
        //                count += g::proc::find_error_match_ions(modifier_table,C_ions_n_c);
        //                count += g::proc::find_error_match_ions(modifier_table,N_ions_c_n);
        //                count += g::proc::find_error_match_ions(modifier_table,N_ions_n_c);
        //                cout<<" error match ions = "<< count<<endl;
        //            }
        cout << "\n选取的是 : " << pI << " 端 !" << endl;
        cout << "N端选取的离子为 : \n";
        for (int pi = 0; pi < arg->ions_n.size(); pi++) {
            cout << " T : " << arg->ions_n[pi].thoe_id
                 << " M : " << arg->ions_n[pi].mono_id
                 << " mono : " << spact->ions_mass_container[arg->ions_n[pi].mono_id]
                 << " radio : " << spact->keyMsValueAbRadio[spact->ions_mass_container[arg->ions_n[pi].mono_id]]
                 << " idex :" << arg->ions_n[pi].index
                 << endl;
        }
        cout << "C端选取的离子为 : \n";
        for (int pi = 0; pi < arg->ions_c.size(); pi++) {
            cout << " T : " << arg->ions_c[pi].thoe_id
                 << " M : " << arg->ions_c[pi].mono_id
                 << " mono : " << spact->ions_mass_container[arg->ions_c[pi].mono_id]
                 << " radio : " << spact->keyMsValueAbRadio[spact->ions_mass_container[arg->ions_c[pi].mono_id]]
                 << " idex :" << arg->ions_c[pi].index
                 << endl;
        }
        //将部分结果写入文件，仅供测试使用
        WriteSequenceAlignment("../wtop_process/data/unknownTest/test-result-c-n.txt", arg->theory_ions_mass_n, spact->ions_mass_container,
                               arg->alignment_ions_n, arg->sub_n, arg->filter_ions_n, arg->ppm_value_n,
                               arg->ref_score_n, arg->peer_score_n);
        WriteSequenceAlignment("../wtop_process/data/unknownTest/test-result-n-c.txt", arg->theory_ions_mass_c, spact->ions_mass_container,
                               arg->alignment_ions_c, arg->sub_c, arg->filter_ions_c, arg->ppm_value_c, arg->ref_score_c,
                               arg->peer_score_c);
        cout << " cut n = " << arg->cut_location_n
             << " cut c = " << arg->cut_location_c << endl;
        cout << "score = " << arg->score
             << ", match Peaks = " << arg->peaks_match
             << " , all peaks = " << spact->ions_mass_container.size()
             << " mut-x = " << arg->mut_x
             << endl;
    }


    bool sort_one_ptm_protein(protein_processor_ptr a, protein_processor_ptr b) {
        return a->one_ptm_match_peaks > b->one_ptm_match_peaks;
    }

    bool sort_one_ptm_protein_2(protein_processor_ptr a, protein_processor_ptr b) {
        if (a->matchTagSize != b->matchTagSize) {
            return a->matchTagSize > b->matchTagSize;
        } else {
            return a->maxTagLen > b->maxTagLen;
        }
    }

    bool sort_one_ptm_protein_3(protein_processor_ptr a, protein_processor_ptr b) {
        return a->complementIonsInProtein > b->complementIonsInProtein;
    }

    bool sort_one_ptm_protein_4(protein_processor_ptr a, protein_processor_ptr b) {
        return a->leastSub < b->leastSub;
    }

    /**
     * 根据ccMap 变化蛋白，并且寻找 匹配的tag
     * @param prIt
     * @param msIt
     * @param ccMap     key-变化前的氨基酸 ，value- 变化后的氨基酸
     * @param set   所有的tag
     * @param m     记录多少个tag被匹配
     * @param filterMass    由于是补充过滤并且是按照tag ，过滤大质量蛋白
     * @return  返回的是匹配tag的最大长度
     */
    int findBestProteinByTag(protein_processor_ptr prIt,
                             msalign_processor_ptr msIt,
                             map<char, char> ccMap,
                             set<string> set,
                             int &m, double filterMass)
    {
        if (fabs(filterMass) > 2000) {
            return 0;
        }
        //得到该条蛋白序列
        string proSeq = prIt->protein_sequence;
        //替换
        for (auto it = ccMap.begin(); it != ccMap.end(); ++it) {
            replace(proSeq.begin(), proSeq.end(), it->first, it->second);
        }
        int count = 0;
        int maxTagSize = 0;
        for (string s : set) {
            if (s.size() < 2) continue;     //如果只有1个氨基酸，则直接过滤
            string rev = string(s.rbegin(), s.rend());
            if (proSeq.find(s) != -1) {
                count++;
                maxTagSize = maxTagSize > s.size() ? maxTagSize : s.size();
            } else if (proSeq.find(rev) != -1) {
                count++;
                maxTagSize = maxTagSize > s.size() ? maxTagSize : s.size();
            }
        }
        m = count;
        return maxTagSize;
    }

    /**
     * 过滤寻找峰值策略
     * @param mono_mass
     * @param theo_mass
     * @param prIt
     * @param msIt
     * @param mdIt
     * @param search_size
     * @param seq_match_peaks
     * @return
     */
    int head_match_peaks_search(const vector<double> &mono_mass, const vector<double> &theo_mass,
                                protein_processor_ptr prIt, msalign_processor_ptr msIt,
                                 int &seq_match_peaks)
    {
        int count = 0;
        vector<int> seq_id;
        vector<double> short_1;
        vector<double> long_1;
        if (mono_mass.size() > theo_mass.size()) {
            long_1 = mono_mass;
            short_1 = theo_mass;
        } else {
            short_1 = mono_mass;
            long_1 = theo_mass;
        }

        for (int i = 0; i < short_1.size(); i++) {      //获取匹配个数
            int n = mylib::speculate_lib::txpb_binarey_search_ex(long_1, long_1.size(), short_1[i]);
            double minPPM = 0;
            if (mylib::data_stream::_ppm_H(short_1[i], long_1[n])) {
                count++;
            }
        }
        if (count < 2) {
            /**
             * 进行N端的搜寻
             * 1、判定最大截断位置
             * 2、添加修饰峰
             * 3、返回取匹配最优的
             */
//             int maxCount = getMaxCutLocETD(prIt,msIt,mdIt,10);
//             seq_match_peaks = maxCount;
            return seq_match_peaks;
        }
        int max_seq_peaks = 0;
        for (int i = 0; i < seq_id.size(); ++i) {
            int c = 1;
            //1 3 4 5 9
            for (int j = i + 1; j < seq_id.size(); ++j) {
                if (seq_id[j] - seq_id[j - 1] < 1) {
                    c++;
                } else if (c > max_seq_peaks) {
                    max_seq_peaks = c;
                    i = j - 1;
                    break;
                }
            }
        }
        seq_match_peaks = max_seq_peaks;
        return count;
    }

    bool selectMassLeast1000(const vector<double> &cT, double precursorMass, double &leastSub, double scopeValue) {
        if (fabs(cT.back() - precursorMass) < scopeValue) {
            leastSub = fabs(cT.back() - precursorMass);
            return true;
        }
        return false;
    }

    /**
     * 过滤处理函数
     * @param msContainer   质谱
     * @param prContainer   蛋白
     * @param one_mod_ptr   修饰
     * @param prsmArg       参数
     */
    void filtrationProcess(const vector<msalign_ptr> &msContainer,
                           const vector<protein_ptr> &prContainer,
                           const config_argument_ptr prsmArg)
    {
//        int count = 1 ;
//        double msSize = msContainer.size();
//        vector<msalign_ptr>::const_iterator ms_it = msContainer.begin() ;
//        for ( ; ms_it != msContainer.end(); ++ms_it, ++count)
//        {
//            std::cout << "\rCurrent thread id : "<<std::this_thread::get_id()<<" . "
//                <<"Number of Ms Filter : " << count <<" th ";
//                std::cout.flush();
//            if ((*ms_it)->getIonsMassContainer().size() < prsmArg->minMonoPeaks
//              || (*ms_it)->getPrecursorMass() < prsmArg->min_precursor_mass) //质谱峰条件
//            {
//                continue;
//            }
//            vector<protein_ptr> match_ions_filters;    //保存由zero 匹配的蛋白，根据匹配个数，最后取top
//            vector<protein_ptr> match_tah_filters;  //保存由tag 匹配的蛋白，最后取top
//            vector<protein_ptr> complement_ions_filter;  //保存由互补离子 匹配的蛋白，最后取top
//            vector<protein_ptr> con;
//            vector<protein_ptr> precursor_mass_filters;  //保存前提质量 匹配的蛋白，最后取top
//            set<string> tagSet;    // 保存所有的tag
//            map<char, char> ccMap;   //保存变化过得氨基酸，前一个为变化的，后一个为变化后的
//            map<int, int> mmMap;    //key为互补离子下标，value 为互补离子下标
//
//            int comIonsSize = 0;  // 互补离子查找到的个数
//
//            if (prContainer.size() == 0) {
//                cout << " error protein database " << endl ;
//                exit(23) ;
//            }
//            vector<protein_ptr>::const_iterator prit = prContainer.begin();
//            map<char,double> acid_map = (*prit)->get_acid_Map();
//            unknown::ptm_utils::getTagByMonoMass((*ms_it)->ions_mass_container,
//                                                 acid_map, tagSet, ccMap);
//
////            makeTagByMonoMass((*prit), (*it), ccMap, tagSet);// 根据该条质谱得到tag 和转变后的ccMap
//
//            //得到质谱互补离子对的map
//            unknown::ptm_utils::getComplementIonsMap((*ms_it)->precursor_mass_double, (*ms_it)->ions_mass_container, mmMap);
//
//            //判定蛋白中互补的占比
//            //            double complementIonsInProtein = isHaveComplementIons((*prit),(*it),one_mod_ptr,mmMap,comIonsSize);
//            //开始过滤
//            for ( ; prit != prContainer.end(); ++prit)
//            {
//                if ((*prit)->getTheoryMassCBy().size() == 0 || (*prit)->getTheoryMassCCz().size() == 0)
//                {
//                    continue;
//                }
//                //这里过滤主要考虑匹配离子的个数
//                int zero_match_peaks = 0;   //完全匹配离子的个数
//                int matchTagSize = 0;      //匹配上tag的个数
//                int maxTagLen = 0;         //最长tag的长度
//                double complementIonsInProtein = 0;
//                bool logLeastScope;   //是否满足范围
//                double leastSub = 0;   //范围内的差值
//                if (prsmArg->BYIonsActions.count((*ms_it)->activation_sub)) {       //如果是BY离子Action
//                    double precursor_mass_gap = (*ms_it)->precursor_mass_double - (*prit)->getTheoryMassCBy().back();
//                    if (precursor_mass_gap > prsmArg->ptmMassMax + prsmArg->unknownPTMsMassMax) {   //如果该蛋白质量太小
//                        continue;
//                    }
//                    // 过滤第一步，根据tag ， 寻找到有几个匹配的tag（matchTagSize） 和 匹配上tag的最大长度（返回值）
//                    maxTagLen = findBestProteinByTag((*prit), (*ms_it), ccMap, tagSet, matchTagSize, precursor_mass_gap);
//
//                    // 过滤第二步，完全匹配
//                    zero_match_peaks = head_match_peaks_search((*ms_it)->ions_mass_container, (*prit)->theoryMassC_By,
//                                                               (*prit), (*ms_it),
//                                                               (*prit)->one_ptm_seq_match_peaks);
//                    logLeastScope = selectMassLeast1000((*prit)->theoryMassC_By,
//                                                        (*ms_it)->precursor_mass_double,
//                                                        (*prit)->leastSub,
//                                                        prsmArg->unknownPTMsFilterProteinMassScope);
//                }
//                if (prsmArg->CZIonsActions.count((*ms_it)->activation_sub)) {   //如果是CZ离子Action
//
//                    double fMass = (*ms_it)->precursor_mass_double - (*prit)->theoryMassC_Cz.back();
//                    if (fMass > prsmArg->ptmMassMax + prsmArg->unknownPTMsMassMax) {   //如果说该蛋白质量太小
//                        continue;
//                    }
//
//                    // 根据tag ， 寻找到有几个匹配的tag（matchTagSize） 和 匹配上tag的最大长度（返回值）
//                    maxTagLen = findBestProteinByTag((*prit), (*ms_it), ccMap, tagSet, matchTagSize, fMass);
//                    zero_match_peaks = head_match_peaks_search((*ms_it)->ions_mass_container, (*prit)->theoryMassC_Cz,
//                                                               (*prit), (*ms_it),
//                                                               (*prit)->one_ptm_seq_match_peaks);
//                    logLeastScope = selectMassLeast1000((*prit)->theoryMassC_Cz, (*ms_it)->precursor_mass_double, (*prit)->leastSub,
//                                                        prsmArg->unknownPTMsFilterProteinMassScope);
//                }
//                //第一步过滤结果保存
//                if (zero_match_peaks > prsmArg->unknownPTMsFilterCountSize) {
////                    (*prit)->one_ptm_match_peaks = zero_match_peaks;
//                    //                    (*prit)->complementIonsInProtein = isHaveComplementIons((*prit),(*it),one_mod_ptr,mmMap,comIonsSize);
//                    match_ions_filters.push_back((*prit));
//                }
//                //第二步过滤结果保存
//                if (maxTagLen > prsmArg->unknownPTMsFilterCountSize) {
////                    (*prit)->matchTagSize = matchTagSize;  //匹配tag数
//                    //                    (*prit)->complementIonsInProtein = isHaveComplementIons((*prit),(*it),one_mod_ptr,mmMap,comIonsSize);
////                    (*prit)->maxTagLen = maxTagLen;
//                    match_tah_filters.push_back((*prit));
//                }
//                //取第一步和第二步交叠
//                if (zero_match_peaks > prsmArg->unknownPTMsFilterCountSize && maxTagLen > prsmArg->unknownPTMsFilterCountSize) {
//                    complement_ions_filter.push_back((*prit));
//                }
//                //第三步过滤结果保存
//                if (logLeastScope) {
//                    precursor_mass_filters.push_back((*prit));
//                }
//            }
//
//            //选取top过滤蛋白
//            int topSelect = prsmArg->unknownFilterOneAndTwoTopSelect;
//            int leastTop = prsmArg->unknownFilterThreeTopSelect;
//            sort(match_ions_filters.begin(), match_ions_filters.end(), sort_one_ptm_protein);   //根据 匹配数量排序
//            sort(match_tah_filters.begin(), match_tah_filters.end(), sort_one_ptm_protein_2);
//            sort(complement_ions_filter.begin(), complement_ions_filter.end(), sort_one_ptm_protein);
//            sort(precursor_mass_filters.begin(), precursor_mass_filters.end(), sort_one_ptm_protein_4);
//            // selectFilterProtein(px,pTag,p);
//            set<protein_processor_ptr> setP;      //用于去重，并且确定最终的过滤蛋白
//            for (int i = 0; i < match_ions_filters.size() && i < topSelect; ++i) {    //取 count top
//                setP.insert(match_ions_filters[i]);
//            }
//            for (int i = 0; i < topSelect && i < match_tah_filters.size(); ++i) {   //取 tag top
//                setP.insert(match_tah_filters[i]);
//            }
//            for (int i = 0; i < topSelect && i < complement_ions_filter.size(); ++i) {   //取 pCom top 上两者
//                setP.insert(complement_ions_filter[i]);
//            }
//            for (int i = 0; i < precursor_mass_filters.size() && i < leastTop; ++i) {    //取 ms top
//                setP.insert(precursor_mass_filters[i]);
//            }
//            for (set<protein_processor_ptr>::iterator itp = setP.begin(); itp != setP.end(); ++itp)
//            {
//                (*ms_it)->candidate_protein_container.push_back((*itp));
//            }
//        }
    }


    int match_peaks_search(const vector<double> &mono_mass, const vector<double> &theo_mass,int search_size) {
        int count = 0;
        //获取匹配个数
        vector<double> short_1 = theo_mass;
        vector<double> long_1 = mono_mass;
//        if (mono_mass.size() > theo_mass.size()) {
//            long_1 = mono_mass;
//            short_1 = theo_mass;
//        } else {
//            short_1 = mono_mass;
//            long_1 = theo_mass;
//        }
        for (int i = 0; i < short_1.size(); i++) {
            int n = mylib::speculate_lib::txpb_binarey_search_ex(long_1, long_1.size(), short_1[i]);
            if (mylib::data_stream::_ppm_H(short_1[i], long_1[n])) {
                count++;
            }
        }
        return count;
    }
//    co = 1 all_peaks_temp = 25 unknownMss = 140.0369 ptmMass = 14.0157
//    match_peaks = 0 ptm_peaks = 0 unk_peaks = 8 ptm_unk_peaks = 17
//    co = 5 all_peaks_temp = 26 unknownMss = 497.2129 ptmMass = 14.0157
//    match_peaks = 0 ptm_peaks = 1 unk_peaks = 8 ptm_unk_peaks = 17
    int match_peaks_search_fa(const vector<double> &mono_mass, const vector<double> &theo_mass,
                              int search_size) {
        int count = 0;
        //获取匹配个数
        vector<double> short_1 = theo_mass;
        vector<double> long_1 = mono_mass;
//        if (mono_mass.size() > theo_mass.size()) {
//            long_1 = mono_mass;
//            short_1 = theo_mass;
//        } else {
//            short_1 = mono_mass;
//            long_1 = theo_mass;
//        }
        if (short_1.size() == 0 || long_1.size() == 0) {
            return 0;
        }
        for (int i = 0; i < short_1.size(); i++) {
            int n = mylib::speculate_lib::txpb_binarey_search_ex(long_1, long_1.size(), short_1[i]);
            if (fabs((short_1[i]) - (long_1[n])) < 1.5) {
                count++;
            }
        }
        return count;
    }


    int getMatchPeakSize(vector<double> &a, vector<double> &b) {
        int count;
        vector<double> short_1;
        vector<double> long_1;
        if (a.size() > b.size()) {
            long_1 = a;
            short_1 = b;
        } else {
            short_1 = b;
            long_1 = a;
        }
        for (int i = 0; i < short_1.size(); i++) {      //获取匹配个数
            int n = mylib::speculate_lib::txpb_binarey_search_ex(long_1, long_1.size(), short_1[i]);
            if (mylib::data_stream::_ppm_H(short_1[i], long_1[n])) {
                count++;
            }
        }
        return count;
    }

    /**
     * 当前提质量  截断质量 修饰质量 未知修饰质量 肽段质量
     * @param prIt
     * @param msIt
     * @param mdIt
     * @param maxPtmMass
     */
    int getMaxCutLocETD(protein_processor_ptr prIt, msalign_processor_ptr msIt, modify_ptr mdIt, int maxCutSize) {
        vector<double> seqMass;
        prIt->get_seq_mass(prIt->protein_sequence, seqMass);
        string pre = msIt->precursor_mass;
        double msPreMass = atof(pre.substr(pre.find('=') + 1, pre.size()).c_str());
        double maxTheoryMass;
        vector<double> TheoryMass;
        if (msIt->activation.find("HCD") != msIt->activation.npos ||
            msIt->activation.find("CID") != msIt->activation.npos) {
            TheoryMass = prIt->theoryMassN_By;
            maxTheoryMass = prIt->theoryMassN_By.back();
        }
        if (msIt->activation.find("ETD") != msIt->activation.npos) {
            TheoryMass = prIt->theoryMassN_Cz;
            maxTheoryMass = prIt->theoryMassN_Cz.back();
        }
        if (fabs(msPreMass - maxTheoryMass) > 1000) {
            return 0;
        }
        vector<double> PtmMass;
        for (map<double, string>::iterator it = mdIt->mapi.begin(); it != mdIt->mapi.end(); ++it) {
            PtmMass.push_back(it->first);
        }
        /**
         * 循环遍历maxCutSize次
         * 添加峰值寻找
         */
        double subMass;
        int bestCount = 0;
        for (int i = 0; i < maxCutSize && i < prIt->protein_sequence.size(); ++i) {
            if (i > 1) {
                subMass = seqMass[i - 1];
                for (int j = i; j < TheoryMass.size(); ++j) {
                    TheoryMass[j] -= subMass;
                }
            }
            vector<double> TheoryPtmMass;
            for (int j = 0; j < PtmMass.size(); ++j) {
                for (int t = i; t < TheoryMass.size(); ++t) {
                    TheoryPtmMass.push_back(TheoryMass[t] + PtmMass[j]);
                }
            }
            int count = getMatchPeakSize(TheoryPtmMass, msIt->ions_mass_container);
            if (count > bestCount) {
                bestCount = count;
            }
        }
        return bestCount;
    }





    void test_filter_result(msalign_processor_ptr it)
    {
        int op = 0;
        ofstream out_test("one_ptm_filter_result.txt", ios::app);
        out_test << it->scans << endl;
        vector<protein_processor_ptr>::iterator itp = it->candidate_protein_container.begin();
        for (; itp != it->candidate_protein_container.end(); ++itp) {
            //                if ( (*itp)->one_ptm_match_peaks > 3 && (*itp)->one_ptm_seq_match_peaks > 1 ) {
            string name = (*itp)->protein_name.substr(0, (*itp)->protein_name.find_first_of(" "));
            bool temp;
//            if (temp) {
//                out_test << name << " , count = " << (*itp)->one_ptm_match_peaks
//                         << " , seq_match peaks = " << (*itp)->one_ptm_seq_match_peaks
//                         << " , matchTagSize = " << (*itp)->matchTagSize
//                         << " , maxTagLen = " << (*itp)->maxTagLen
//                         << ", complementIonsInProtein = " << (*itp)->complementIonsInProtein << endl;
//                break;
//            }
            out_test << name << " , count = " << (*itp)->one_ptm_match_peaks
                     << " , seq_match peaks = " << (*itp)->one_ptm_seq_match_peaks
                     << " , matchTagSize = " << (*itp)->matchTagSize
                     << " , maxTagLen = " << (*itp)->maxTagLen
                     << ", complementIonsInProtein = " << (*itp)->complementIonsInProtein << endl;
        }
        out_test << endl;
        out_test.close();
    }

    /**
     * 判断是否有互补离子的实现
     * @param prIt
     * @param msIt
     * @param iiMap
     * @param nT
     * @param cT
     * @param onePtmMass
     */
    double judgeComplementIons(protein_processor_ptr prIt, msalign_processor_ptr msIt, map<int, int> &iiMap,
                               const int cutMax, const int cutMin,
                               const vector<double> &nT, const vector<double> &cT,
                               const vector<double> &onePtmMass) {
        string prSeq = prIt->protein_sequence;
        vector<double> mono = msIt->ions_mass_container;
        set<int> haveCompIons; //记录的key表示含有
        int startIndex = mono.size() - cutMax; //起始位置
        int endIndex = mono.size() - cutMin;   //最终位置
        int cutSize = startIndex - 1; //每次截断之后要更新N端的理论
        double prMass = msIt->precursor_mass_double; //蛋白质前体质量
        map<double, double> ptmUnknown; //key 为 已知修饰质量，value为未知修饰质量
        int proSeqLen = prIt->protein_sequence.size();
        int complementIonsSize = 0;
//        for (map<int,int>::iterator itMap = iiMap.begin();itMap!=iiMap.end();++itMap) {
//            cout << " first = " << mono[itMap->first] << " second = " << mono[itMap->second] <<endl;
//        }
        //判定每一对互补离子，只要判定一对，其余对皆按此
        for (map<int, int>::iterator itMap = iiMap.begin(); itMap != iiMap.end(); ++itMap) {
            double m1 = mono[itMap->first];         // 互补峰1 未知NC
            double m2 = mono[itMap->second];        // 互补峰2 未知NC
            //步骤：
            // 1、先查找该形式下 C 端是否存在值。存在 则判定下一步，超过截断（未知修饰超过500da）则直接滤掉
            // 2、从N端找，此时确定的是C端离子右边的修饰质量，以及N端离子的位置，截断则需要去遍历判定（相减即可）
            for (int j = cutMax; j != cutMin; j--) { //每个蛋白质形
                int m1Index = -1;
                int m2Index = -1;
                unknown::ptm_utils::getUnknownMassByPtmMass(cT, onePtmMass, prMass, j, ptmUnknown,
                                                            500);//该蛋白质形下，每个已知ptm的未知修饰质量getUnknownMassByPtmMass
                if (ptmUnknown.size() == 0) {   //如果没有满足条件的未知修饰质量 直接退出当前的蛋白质形
                    break;
                }
                //得到该蛋白形下截断的质量
                double nTCutMass = prIt->getNCutMass(nT.size() - j - 1);
                //步骤一：先确定C端的离子是m1 还是 m2 , 在cT中确定一个m1或m2如果，则可以知道该蛋白形下另外一个
//                bool loop = u::proc::getComplementIonsCTIndex(m1,m2,m1Index,m2Index,j+1,cT,nT,ptmUnknown,nTCutMass,1.5);
//                if (loop) {
//                    complementIonsSize ++ ;
//                    cout << "cutSize = " << nT.size() - j - 1
//                    << " m1 = " << m1 << " m2 = "<< m2
//                    << endl ;
//                    break ;
//                }
            }
        }
//        cout << " iiMap.size = " << iiMap.size() << " complementIonsSize = " << complementIonsSize <<endl;
        return (double) complementIonsSize / (double) iiMap.size() * 100;
    }

    /**
     * 根据质谱得到tag
     * @param prIt
     * @param msIt
     * @param ccMap key-变化前的氨基酸 ，value- 变化后的氨基酸
     * @param set 所有的tag
     * @return
     */
    void makeTagByMonoMass(msalign_processor_ptr msIt,
                          map<char, char> &ccMap,map<char,double> & acid_map ,
                          set<string> &set) {
        unknown::ptm_utils::getTagByMonoMass(msIt->ions_mass_container, acid_map, set, ccMap);
    }

    //如果有互补离子，过滤蛋白理论序列中必须要有，否则直接干掉
    double isHaveComplementIons(protein_processor_ptr prIt, msalign_processor_ptr msIt,
                                modify_ptr oneModPtm, map<int, int> &iiMap, int comIonsSize) {
        if (iiMap.size() == 0) {
            return 0;
        }
        vector<double> mono = msIt->ions_mass_container;
        vector<double> nT;
        vector<double> cT;
        vector<double> ptmMass;
        oneModPtm->getModMassSeq(ptmMass);  //得到可变ptm质量序列
        //寻找到最大和最小截断数
        int cutMin = 0;
        int cutMax = 0;
        //init理论峰
        if (msIt->activation.find("HCD") != msIt->activation.npos ||
            msIt->activation.find("CID") != msIt->activation.npos) {
            nT = prIt->theoryMassN_By;
            cT = prIt->theoryMassC_By;
            //传入c端序列并判断
            unknown::ptm_utils::getCutMinAndMax(cT, msIt->precursor_mass_double, 1000, cutMin, cutMax);
        }
        //unk_mass = prMass - arg->thero_mass_N_C.back() - it_mp->first - subMassETD;
        if (msIt->activation.find("ETD") != msIt->activation.npos) {
            double subMassETD = 18.01056 - 1.9919;
            nT = prIt->theoryMassN_Cz;
            cT = prIt->theoryMassC_Cz;
            for (int i = 0; i < cT.size(); ++i) {
                cT[i] -= subMassETD;
            }
            unknown::ptm_utils::getCutMinAndMax(cT, msIt->precursor_mass_double, 1000, cutMin, cutMax);
            for (int i = 0; i < cT.size(); ++i) {
                cT[i] += subMassETD;
            }
            msIt->precursor_mass_double -= subMassETD;
        }
        if (cutMax == 0 && cutMin == 0) {
            comIonsSize = 0;
            return 0;
        }
//        cout << " min = " << cutMin <<" max = "<< cutMax <<endl;
//        cout << " seqSize = " << prIt->theoryMassC_By.size()<<endl;
        //判断是否存在（思路：截断，后n+1个蛋白质形，每个蛋白质形进行判断（2个修饰）
        //2个修饰进行遍历，在右边的遍历，并且在右边分为N端和C端的遍历
        double t = judgeComplementIons(prIt, msIt, iiMap, cutMax, cutMin, nT, cT, ptmMass);
        return t;
    }


    /**
     * 根据匹配离子数，匹配tag数，动态取过滤后的protein
     * @param one
     * @param two
     * @param pc
     */
    void selectFilterProtein(vector<protein_processor_ptr> &one, vector<protein_processor_ptr> &two, msprP &pc) {
        int oneMinCount = 0; //记录最小的匹配
        int selectPro = 100; //选取top20
        if (one.size() > selectPro) {
            oneMinCount = one[selectPro - 1]->one_ptm_match_peaks;
            for (protein_processor_ptr p : one) {
                if (p->one_ptm_match_peaks < oneMinCount) {
                    break;
                }
                pc->prContainer.push_back(p);
            }
        } else {
            for (protein_processor_ptr p : one) {
                pc->prContainer.push_back(p);
            }
        }
        int selectTwo = 100;
        if (two.size() > selectTwo) {
            oneMinCount = two[selectTwo - 1]->matchTagSize;
            for (protein_processor_ptr p : two) {
                if (p->matchTagSize < oneMinCount) {
                    break;
                }
                pc->prContainer.push_back(p);
            }
        } else {
            for (protein_processor_ptr p : two) {
                pc->prContainer.push_back(p);
            }
        }
    }



    int judgeComplementIonsByCut(protein_processor_ptr prIt, msalign_processor_ptr msIt, const map<int, int> &iiMap,
                                 const int cutMax, const int cutMin,
                                 const vector<double> &nT, const vector<double> &cT,
                                 const vector<double> &onePtmMass, double &modMass) {
        string prSeq = prIt->protein_sequence;   //蛋白序列
        vector<double> mono = msIt->ions_mass_container;  //质谱质量
        double prMass = msIt->precursor_mass_double; //蛋白质前体质量
        map<int, vector<int> > iVmap; //key 为互补的iiMap的key ,后一个为可行的截断数
        map<int, map<int, vector<int> > > mIdCutModId; //key对应mId, value (key - Cut数, value - ModIndex)
        int proSeqLen = prIt->protein_sequence.size();
        int complementIonsSize = 0;
        //判定每一对互补离子，只要判定一对，其余对皆按此
        int logM = 0;
        for (map<int, int>::const_iterator itMap = iiMap.begin(); itMap != iiMap.end(); ++itMap) {
            double m1 = mono[itMap->first];         // 互补峰1 未知NC
            double m2 = mono[itMap->second];        // 互补峰2 未知NC
//            cout << " m1 = " <<m1 << " m2 = "<< m2 <<endl;
            //步骤：
            // 1、先查找该形式下 C 端是否存在值。存在 则判定下一步，超过截断（未知修饰超过500da）则直接滤掉
            // 2、从N端找，此时确定的是C端离子右边的修饰质量，以及N端离子的位置，截断则需要去遍历判定（相减即可）
            for (int j = cutMax; j != cutMin; j--) { //每个蛋白质形
                int m1Index = -1;
                int m2Index = -1;
                map<double, double> ptmUnknown; //key 为 已知修饰质量，value为未知修饰质量
                vector<int> trueModIndex;      //该对互补离子，当前截断下可行的修饰下标
                unknown::ptm_utils::getUnknownMassByPtmMass(cT, onePtmMass, prMass, j, ptmUnknown,
                                                            500);//该蛋白质形下，每个已知ptm的未知修饰质量getUnknownMassByPtmMass
                if (ptmUnknown.size() == 0) {   //如果没有满足条件的未知修饰质量 直接退出当前的蛋白质形
                    continue;
                }

                double nTCutMass = prIt->getNCutMass(nT.size() - j - 1);    //得到该蛋白形下截断的质量
                // 步骤一：先确定C端的离子是m1 还是 m2 , 在cT中确定一个m1或m2如果，则可以知道该蛋白形下另外一个
                bool loop = unknown::ptm_utils::getComplementIonsCTIndex(m1, m2, m1Index, m2Index, j + 1, cT, nT, trueModIndex,
                                                                         ptmUnknown, nTCutMass, 1.2);
                if (loop) { //说明当前蛋白质形存在该互补离子
//                    if (iVmap.count(itMap->first)) {
//                        iVmap[itMap->first].push_back(nT.size() - j - 1);
//                    } else {
//                        vector<int> c ;
//                        c.push_back(nT.size() - j - 1);
//                        iVmap.insert(make_pair(itMap->first,c));
//                    }
                    if (mIdCutModId.count(itMap->first)) {
                        mIdCutModId[itMap->first].insert(make_pair(nT.size() - j - 1, trueModIndex));
                    } else {
                        map<int, vector<int> > t;
                        t.insert(make_pair(nT.size() - j - 1, trueModIndex));
                        mIdCutModId.insert(make_pair(itMap->first, t));
                    }
                }
            }
        }
        int cm = 0;
        int mLen = onePtmMass.size() - 1;
        for (auto iit = mIdCutModId.begin(); iit != mIdCutModId.end(); iit++) {
//            cout << "mId = " << iit->first  ;
            map<int, vector<int>> t = iit->second;
            for (auto ic = t.begin(); ic != t.end(); ++ic) {
                vector<int> v = ic->second;
                if (ic->first >= cm) {
                    cm = ic->first + 1;
                }
//                cout <<" cut = " << ic->first ;
//                cout << " modIndex = " ;
//                for( int i : v) {
//                    cout  << i <<" " ;
//                }
//                cout << endl;
            }
        }
        //得到最佳截断
        int bestCut = -1;
        int rw[cm][mLen];  //将截断，匹配构造成二维数组，在二位数组中寻求最大截断
        memset(rw, 0, sizeof(rw));
        for (auto iit = mIdCutModId.begin(); iit != mIdCutModId.end(); iit++) {
            map<int, vector<int>> t = iit->second;
            for (auto ic = t.begin(); ic != t.end(); ++ic) {
                vector<int> v = ic->second;
                int m = ic->first;
                for (int i : v) {
                    rw[m][i] += 1;
                }
            }
        }
        int tempMax = 0;
        int bestModIndex = -1;
        for (int i = 0; i < cm; ++i) {
            for (int j = 0; j < mLen; ++j) {
                if (rw[i][j] > tempMax) {
                    tempMax = rw[i][j];
                    bestCut = i;
                    bestModIndex = j;
                }
//                cout << rw[i][j] << " " ;
            }
//            cout << endl ;
        }
        int sc = 0;
        for (int i = 0; i < onePtmMass.size(); ++i) {
            if (fabs(onePtmMass[i]) < 0.00001) {
                continue;
            }
            if (sc == bestModIndex) {
                modMass = onePtmMass[i];
                break;
            }
            sc++;
        }
        return bestCut;
//        for( auto it = iVmap.begin() ; it!= iVmap.end(); ++it) {
//            cout << " first = "<< it->first ;
//            for (int d : it->second) {
//                cout <<"  "<< d ;
//            }
//            cout <<endl;
//        }
//        map<int,int> ccMap ;    //key-截断数 value 出现的次数
//            for( auto it = iVmap.begin() ; it!= iVmap.end(); ++it) {
//                for (int d : it->second) {
//                    if (ccMap.count(d)) {
//                        ccMap[d]++ ;
//                    } else {
//                        ccMap.insert(make_pair(d,1));
//                    }
//                }
//            }
//            int bestCut = -1 ;
//            int temp = 0 ;
//            for( auto it = ccMap.begin(); it != ccMap.end(); ++it) {
//        //            cout << it->first << " " <<it->second <<endl;
//                if( it->second > temp) {
//                    bestCut = it->first ;
//                    temp = it->second ;
//                }
//            }
    }

    //22-10-8日更新，判断最合适的N端截断
    //思路为，如果有互补离子，从互补离子下手，否则就是judge_best_n_cut
    void judgeBestNCutSize(protein_processor_ptr prIt, msalign_processor_ptr msIt, modify_ptr oneModPtm, const map<int, int> &iiMap,
                           const config_argument_ptr &prsmArg, double &mod, int &cutN) {
        vector<double> mono = msIt->ions_mass_container;
        vector<double> nT;
        vector<double> cT;
        vector<double> ptmMass;
//        double subMassETD = 18.01056 - 1.9919;
        double subMassETD = prsmArg->sub_etd;
        oneModPtm->getModMassSeq(ptmMass);  //得到可变ptm质量序列
        int cutMin = 0;     //最大截断位置
        int cutMax = 0;     //最小截断位置
        //init理论峰
        if (prsmArg->BYIonsActions.count(msIt->activation_sub)) {
            nT = prIt->theoryMassN_By;
            cT = prIt->theoryMassC_By;
            //传入c端序列并判断 //获取最大截断和最小截断
            unknown::ptm_utils::getCutMinAndMax(cT, msIt->precursor_mass_double, 1000, cutMin, cutMax);
        }
        //unk_mass = prMass - arg->thero_mass_N_C.back() - it_mp->first - subMassETD;
        if (prsmArg->CZIonsActions.count(msIt->activation_sub)) {
            nT = prIt->theoryMassN_Cz;
            cT = prIt->theoryMassC_Cz;
            for (int i = 0; i < cT.size(); ++i) {
                cT[i] -= subMassETD;
            }
            unknown::ptm_utils::getCutMinAndMax(cT, msIt->precursor_mass_double, 1000, cutMin, cutMax);
            for (int i = 0; i < cT.size(); ++i) {
                cT[i] += subMassETD;
            }
            msIt->precursor_mass_double -= subMassETD;
        }
        if (cutMin == 0 && cutMax == 0) {
            cutN = -1;
            return;
        }
//        if (prIt->proteinName.find(">sp|P62805|H4_HUMAN") != prIt->proteinName.npos) {
//            cout << " max = " << cutMax << " min = " << cutMin <<endl;
//            for( double d : cT) {
//                cout << d <<endl;
//            }
//            cout <<"d["<<cutMax<<"] = " << cT[cutMax]
//            <<"d["<<cutMin<<"] = " << cT[cutMin]
//            <<endl;
//            cout <<" prMass = "<< msIt->Pr_mass <<endl;
//        }
//        cout << " seqSize = " << prIt->theoryMassC_By.size()<<endl;
        //判断是否存在（思路：截断，后n+1个蛋白质形，每个蛋白质形进行判断（2个修饰）
        //2个修饰进行遍历，在右边的遍历，并且在右边分为N端和C端的遍历
        cutN = judgeComplementIonsByCut(prIt, msIt, iiMap, cutMax, cutMin, nT, cT, ptmMass, mod);
        if (prsmArg->CZIonsActions.count(msIt->activation_sub)) {
            msIt->precursor_mass_double += subMassETD;
        }
        if (cutN < 0) {
            cutN = -1;
        }
    }

    //判断未知修饰最合适的截断
    void judge_best_n_cut_CID(prsm::modify_ptr mp, protein_processor_ptr pro, msalign_processor_ptr spact, const config_argument_ptr &prsmArg,
                              int &cut_n, int cut_c) {
        vector<double> mono_seq_mass = spact->ions_mass_container;    //临时记录质谱质量
        double max_ptm_mass = 0;        //记录下最大修饰质量
        double min_ptm_mass = 0;        //记录下最小修饰质量
        double precursor_mass = atof(spact->precursor_mass.substr(spact->precursor_mass.find("=") + 1,
                                                                  spact->precursor_mass.size() + 1).c_str());
        vector<double> ptm_mass;    //修饰的质量，方便直接取值
        map<double, string> cmap = (*mp).mapi;  //保存修饰map
        vector<double> seq_mass; //记录氨基酸序列质量
        pro->get_seq_mass(pro->protein_sequence, seq_mass); //记录氨基酸质量序列
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
        double max_theory_mass = pro->theoryMassC_By.back();
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
                if (fabs(it_mp->first) < prsmArg->min_value) {
                    continue;
                }
                double max_theory_mass_1 = pro->theoryMassC_By.back();
                theory_ptm_unknow.assign(pro->theoryMassN_By.begin() + co, pro->theoryMassN_By.end());
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
                int match_peaks = match_peaks_search(mono_seq_mass, theory_ptm_unknow, mono_seq_mass.size());
                int ptm_peaks = match_peaks_search(mono_seq_mass, theory_with_ptm, mono_seq_mass.size());
//                int unk_peaks = match_peaks_search_fa(mono_seq_mass, theory_with_unk, mono_seq_mass.size());
                int unk_peaks = 0 ;
                int ptm_unk_peaks = match_peaks_search_fa(mono_seq_mass, theory_with_ptm_unk, mono_seq_mass.size());
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
            }   //遍历修饰
        }  //遍历截断
//        cout <<"name" << pro->proteinName<<endl;
        cut_n = cut_n_best; //纪录当前最佳截断
    }

    //判断未知修饰最合适的截断
    /**
     * 判断未知修饰最合适的截断
     * @param mp
     * @param pro
     * @param spact
     * @param cut_n
     * @param cut_c
     */
    void judge_best_n_cut_ETD(prsm::modify_ptr mp, protein_processor_ptr pro, msalign_processor_ptr spact, const config_argument_ptr &prsmArg,
                              int &cut_n, int cut_c) {
        vector<double> mono_seq_mass = spact->ions_mass_container;    //临时记录质谱质量
        double max_ptm_mass = 0;        //记录下最大修饰质量
        double min_ptm_mass = 0;        //记录下最小修饰质量
//        double subMassETD = 18.01056 - 1.9919;
        double subMassETD = prsmArg->sub_etd;
        double precursor_mass = atof(spact->precursor_mass.substr(spact->precursor_mass.find("=") + 1,
                                                                  spact->precursor_mass.size() + 1).c_str()) - subMassETD;
        vector<double> ptm_mass;    //修饰的质量，方便直接取值
        map<double, string> cmap = (*mp).mapi;  //保存修饰map
        vector<double> seq_mass;    //记录氨基酸序列质量
        pro->get_seq_mass(pro->protein_sequence, seq_mass); //记录氨基酸质量序列
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
        double max_theory_mass = pro->theoryMassC_Cz.back();    //从C端判断
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
                double max_theory_mass_1 = pro->theoryMassC_Cz.back();
                theory_ptm_unknow.assign(pro->theoryMassN_Cz.begin() + co, pro->theoryMassN_Cz.end());  // 初始化N端理论峰序列
                vector<double> theory_with_ptm;     //插入ptm之后的序列
                vector<double> theory_with_unk;     //插入未知ptm之后的序列
                vector<double> theory_with_ptm_unk;     //插入混合ptm之后的序列
                double sub_mass = 0;    //记录截断的氨基酸质量
                for (int i = 0; i < co; ++i) {      //得到当前截断下的氨基酸质量,co是当前截断的个数
                    sub_mass += seq_mass[i];
                }
                max_theory_mass_1 -= sub_mass;
//                cout << " max_theory_mass_1 = " <<max_theory_mass_1 <<endl;
                double unk_mass = precursor_mass - max_theory_mass_1 - it_mp->first;    //记录未知修饰质量
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
                int match_peaks = match_peaks_search(mono_seq_mass, theory_ptm_unknow, mono_seq_mass.size());   //判断是否是掉H多H
                int ptm_peaks = match_peaks_search(mono_seq_mass, theory_with_ptm, mono_seq_mass.size());       //判断带修饰的峰
//                int unk_peaks = match_peaks_search_fa(mono_seq_mass, theory_with_unk, mono_seq_mass.size());    //判断是否在1.5范围内
                int ptm_unk_peaks = match_peaks_search_fa(mono_seq_mass, theory_with_ptm_unk, mono_seq_mass.size());
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

//        cout <<"name : " << pro->proteinName<< " :::: best cut = " << cut_n<< endl;
    }


    //只有未知修饰的情况
    void unknown_match_search_ptmZero(const vector<double> &mono_mass, const vector<double> &theo_mass,
                                      vector<double> &deviSeq, double unknownPtmMass, int searchSize) {
        //获取匹配个数
        vector<double> theorySeq_UnknownPtm = theo_mass;
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
    }

    //未知修饰和已知修饰的情况
    void unknown_match_search_ptmNum(const vector<double> &mono_mass, const vector<double> &theo_mass,
                                     vector<double> &deviSeq,
                                     double ptmMass, double unknownPtmMass, int searchSize) {
        vector<double> theorySeq_UnknownPtm = theo_mass;
        vector<double> theorySeq_Ptm_UnknownPtm = theo_mass;
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

    }

    int judge_best_unknown_ptm_mass(msalign_processor_ptr spact, temp_prsm_argument_ptr arg,
                                    double ptmMass, double unknown_mass,
                                    int searchN, int searchC) {
        vector<double> n_devi_seq;
        vector<double> c_devi_seq;
        if (ptmMass == 0) {
            unknown_match_search_ptmZero(spact->ions_mass_container, arg->theory_ions_mass_c,
                                         c_devi_seq, unknown_mass, searchC);
            unknown_match_search_ptmZero(spact->ions_mass_container, arg->theory_ions_mass_n,
                                         n_devi_seq, unknown_mass, searchN);
        } else {
            unknown_match_search_ptmNum(spact->ions_mass_container, arg->theory_ions_mass_c, c_devi_seq,
                                        ptmMass, unknown_mass, searchC);
            unknown_match_search_ptmNum(spact->ions_mass_container, arg->theory_ions_mass_n, n_devi_seq,
                                        ptmMass, unknown_mass, searchN);
        }
        return n_devi_seq.size() + c_devi_seq.size();
//        cout <<" n size = " << n_devi_seq.size() << " c size = " << c_devi_seq.size()<<endl ;
//        for(int i = 0 ; i < n_devi_seq.size() ; ++i ) {
//            cout << n_devi_seq[i] << endl ;
//        }
//        cout <<" ------------------------- " <<endl ;
//        for(int i = 0 ; i < c_devi_seq.size() ; ++i ) {
//            cout << c_devi_seq[i] << endl ;
//        }
    }


    int get_unknown_Ptm_SearchSize(const vector<double> &mono_mass, const vector<double> &theo_mass) {
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

        int searchSize = mIndex.size() - 5;
        if (searchSize < 0) {
            return 0;
        }
        return mIndex[searchSize];
    }

    void initTheoryMass(vector<double> &theoryMass, vector<double> &theoryMassC,
                        protein_processor_ptr pro, int cutSizeN, int cutSizeC) {
        vector<double> seqMass;
        pro->get_seq_mass(pro->protein_sequence, seqMass);
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
            mass2 += seqMass[t];
        }
        for (int i = 0; i < theoryMassC.size(); ++i) {
            theoryMassC[i] -= mass2;
        }
    }

    void initial_Unknown_ptmTable(temp_prsm_argument_ptr arg, vector<double> &modifierTable, map<double, string> &mmap,
                                  double ptmMass, string ptm, double bestUnknown) {
        modifierTable.clear();
        mmap.clear();
        mmap.insert(make_pair(0, "0"));
        if (fabs(ptmMass) <= 1e-15) {
            mmap.insert(make_pair(bestUnknown, "*"));
            arg->var_mod_mass = bestUnknown;
        } else {
            mmap.insert(make_pair(ptmMass, ptm));
            mmap.insert(make_pair(bestUnknown, "*"));
            mmap.insert(make_pair(bestUnknown + ptmMass, "*" + ptm + ""));
            arg->var_mod_mass = ptmMass + bestUnknown;
        }
        for (map<double, string>::iterator it = mmap.begin(); it != mmap.end(); ++it) {
            modifierTable.push_back(it->first);
        }

    }

    //2023-1-8 final merge
    double DetermineTheMostAppropriateUnknownModifierMass(double PtmMass, double prMass, double unk_mass, temp_prsm_argument_ptr arg,
                                                          msalign_processor_ptr spact, protein_processor_ptr pro, config_argument_ptr prsmArg)
    {
        vector<double> unknownMassSeq;
        // 前提质量判断的时候
        double flog = prsmArg->unknownPTMsMassFlog;
        for (int i = 0; i < prsmArg->unknownPTMsMassFlogFloatSize; i++) {
            unknownMassSeq.push_back(unk_mass + flog);
            flog += prsmArg->unknownPTMsMassFlogChangeMass;
        }
        //确定搜索范围是从倒数第五个完全匹配的离子开始
        int searchSize_C = 0 , searchSize_N = 0 ;
//        if (prsmArg->BYIonsActions.count(spact->actionSub)) {
            searchSize_C = get_unknown_Ptm_SearchSize(spact->ions_mass_container, pro->theoryMassC_By);
            searchSize_N = get_unknown_Ptm_SearchSize(spact->ions_mass_container, pro->theoryMassN_By);
//        }
//        if (prsmArg->CZIonsActions.count(spact->actionSub)) {
//            searchSize_C = get_unknown_Ptm_SearchSize(spact->mono_mass, pro->theoryMassC_By);
//            searchSize_N = get_unknown_Ptm_SearchSize(spact->mono_mass, pro->theoryMassN_By);
//        }

        //建立理论峰值,判断最佳合适
        double bestCount = 0;
        double bestUnknown = 0;        //真正的修饰
        for (double unknownMass : unknownMassSeq) {  //确定真正的未知修饰
            double ppm = 0;
            bool boolean = false ;
            if (prsmArg->BYIonsActions.count(spact->activation_sub)) {
                boolean = mylib::etd_process::_PPM_H(unknownMass + PtmMass + arg->theory_ions_mass_c.back(), prMass, ppm) ;
            }
            if (prsmArg->CZIonsActions.count(spact->activation_sub)) {
                boolean = mylib::etd_process::_PPM_H(unknownMass + PtmMass + arg->theory_ions_mass_c.back() + prsmArg->sub_etd, prMass, ppm);
            }
            if (boolean) {
                int count = judge_best_unknown_ptm_mass(spact, arg, PtmMass, unknownMass, searchSize_N,searchSize_C);
                if (count > bestCount) {
                    bestCount = count;
                    bestUnknown = unknownMass;
                }
            }
        }
        if (fabs(bestUnknown) < 0.1) {
            return 0;
        }
        return bestUnknown;
    }
    /**
     * 判断最合适的未知修饰质量
     * @param PtmMass
     * @param arg
     * @param spact
     * @param pro
     * @return
     */
    double judge_Unknown_Mass(double PtmMass, double prMass, double unk_mass, temp_prsm_argument_ptr arg, msalign_processor_ptr spact,
                              protein_processor_ptr pro) {
        vector<double> unknownMassSeq;
        // 前提质量判断的时候
        double flog = -1.5;
        for (int i = 0; i < 30; i++) {     // -1.5 - 1.5 之间的未知修饰范围
            unknownMassSeq.push_back(unk_mass + flog);
            flog += 0.1;
        }
        //确定搜索范围是从倒数第五个完全匹配的离子开始
        int searchSize_C = get_unknown_Ptm_SearchSize(spact->ions_mass_container, pro->theoryMassC_By);
        int searchSize_N = get_unknown_Ptm_SearchSize(spact->ions_mass_container, pro->theoryMassN_By);
        //建立理论峰值,判断最佳合适
        double bestCount = 0;
        double bestUnknown = 0;        //真正的修饰
        for (double unknownMass : unknownMassSeq) {  //确定真正的未知修饰
            double ppm = 0;
            if (mylib::etd_process::_PPM_H(unknownMass + PtmMass + arg->theory_ions_mass_c.back(), prMass, ppm)) {
                int count = judge_best_unknown_ptm_mass(spact, arg, PtmMass, unknownMass, searchSize_N,
                                                        searchSize_C);
                if (count > bestCount) {
                    bestCount = count;
                    bestUnknown = unknownMass;
                }
            }
        }
//        cout << " best Unknown mass = " << bestUnknown << endl;
        if (fabs(bestUnknown) < 0.1) {
            return 0;
        }
        return bestUnknown;
    }

    double judge_Unknown_Mass_ETD(double PtmMass, double prMass, double unk_mass, temp_prsm_argument_ptr arg, msalign_processor_ptr spact,
                                  protein_processor_ptr pro) {
        vector<double> unknownMassSeq;
        // 前提质量判断的时候
        double flog = -1.5;
        for (int i = 0; i < 30; i++) {     // -1.5 - 1.5 之间的未知修饰范围
            unknownMassSeq.push_back(unk_mass + flog);
            flog += 0.1;
        }
        //确定搜索范围是从倒数第五个完全匹配的离子开始
        int searchSize_C = get_unknown_Ptm_SearchSize(spact->ions_mass_container, pro->theoryMassC_Cz);
        int searchSize_N = get_unknown_Ptm_SearchSize(spact->ions_mass_container, pro->theoryMassN_Cz);
        //建立理论峰值,判断最佳合适
        double bestCount = 0;
        double bestUnknown = 0;        //真正的修饰
        for (double unknownMass : unknownMassSeq) {  //确定真正的未知修饰
            double ppm = 0;
            if (mylib::etd_process::_PPM_H(unknownMass + PtmMass + arg->theory_ions_mass_c.back() + 18.01056 - 1.9919, prMass, ppm)) {
//                cout <<" unknownMass = " << unknownMass
//                << " ppm = " << ppm <<endl;
                int count = judge_best_unknown_ptm_mass(spact, arg, PtmMass, unknownMass, searchSize_N,searchSize_C);
                if (count > bestCount) {
                    bestCount = count;
                    bestUnknown = unknownMass;
                }
            }
        }
//        cout << " best Unknown mass = " << bestUnknown << endl;
        if (fabs(bestUnknown) < 0.1) {
            return 0;
        }
        return bestUnknown;
    }

    /**
     * 没有过滤流程的prsm
     * @param one_mod_ptr
     * @param msContainer
     * @param prsmArg
     * @param iiMap
     * @param outFileName
     */
    void  one_ptm_search_process(modify_ptr one_mod_ptr, vector<msalign_processor_ptr> &msContainer, config_argument_ptr prsmArg,
                                 const string &variablePTMsFileOutPath, string & outFileName)
    {
        vector<msalign_processor_ptr>::iterator msIt = msContainer.begin() ;
        for ( ; msIt != msContainer.end(); ++msIt) {
            vector<protein_processor_ptr>::iterator pr_it_begin = (*msIt)->candidate_protein_container.begin();  //过滤蛋白库起始位置
            vector<protein_processor_ptr>::iterator pr_it_end = (*msIt)->candidate_protein_container.end();     //过滤蛋白库终止位置
            for (; pr_it_begin != pr_it_end; ++pr_it_begin) {
                if (prsmArg->BYIonsActions.count((*msIt)->activation_sub)) {
                    unknow_ptms_prsm_process_CID(one_mod_ptr, (*pr_it_begin),(*msIt), prsmArg, (*msIt)->compIonsIndexMap,
                                                 variablePTMsFileOutPath);    //修饰，蛋白，质谱
                }
                if (prsmArg->CZIonsActions.count((*msIt)->activation_sub)) {
                    unknow_ptms_prsm_process_ETD(one_mod_ptr, (*pr_it_begin), (*msIt), prsmArg, (*msIt)->compIonsIndexMap,
                                                 variablePTMsFileOutPath);
                }
            }
            //输出prsm结果
            modify_ptr resultMPtr = (*msIt)->unknownModPtr;
            writeResult((*msIt), outFileName, resultMPtr);
//            delete resultMPtr;
        }
    }
    /**
     * 2023-1-6 deprecate
     * 独立出未知ptm识别修饰结果,与已知ptm识别结果分离
     * 主要流程:
     * 1、每条质谱与输入蛋白质库进行过滤，得出结果，取前20条。
     * 2、将质谱与前20条蛋白进行鉴定，最终得出prsm结果
     * @param msContainer
     * @param prContainer
     * @param ptmMassSeq
     * @param one_ptm_filter_result
     */
    void one_ptm_filter(vector<msalign_processor_ptr> &msContainer, vector<protein_processor_ptr> &prContainer, vector<double> &ptmMassSeq,
                        vector<msprP> &one_ptm_filter_result, modify_ptr one_mod_ptr, const string &outFileName,
                        config_argument_ptr prsmArg) {
        double max_ptm_mass = 0;
        if (ptmMassSeq.size() != 0) {
            sort(ptmMassSeq.begin(), ptmMassSeq.end());
            max_ptm_mass = ptmMassSeq.back();
        }
        /**
         * 过滤逻辑
         * 遍历所有的质谱和蛋白进行对比
         */
        vector<msalign_processor_ptr>::const_iterator it = msContainer.begin();
        double prSize = prContainer.size();
        int count = 1;
        int msSize = msContainer.size();
        //        string log = "4285";
        //        bool co = false;
        for (; it != msContainer.end(); ++it, ++count)
        {
            if (count % 5 == 0)
            {
                cout << "\rone unknown has deal  " << (double) count / (double) msSize * 100 << " % ";
            }
            //            if ((*it)->Scans.find(log) != (*it)->Scans.npos) {
            //                co = true ;
            //            }
            //            if (!co) {
            //                continue;
            //            }
            int search_size = (*it)->ions_mass_container.size(); //搜索大小
            msprP p = new mspr();
            vector<protein_processor_ptr> px;    //保存由zero 匹配的蛋白，根据匹配个数，最后取top20
            vector<protein_processor_ptr> pTag; //保存由tag 匹配的蛋白，最后取top20
            vector<protein_processor_ptr> pCom; //保存由互补离子 匹配的蛋白，最后取top20
            vector<protein_processor_ptr> con;
            vector<protein_processor_ptr> least1000;
            p->it = (*it);
            if ((*it)->ions_mass_container.size() < prsmArg->minMonoPeaks || (*it)->precursor_mass_double < prsmArg->min_precursor_mass) { //质谱峰条件
                continue;
            }
            vector<protein_processor_ptr>::iterator prit = prContainer.begin();
            double count_s = 0;     //用于计数
            set<string> tagSet;    // 保存所有的tag
            map<char, char> ccMap;   //保存变化过得氨基酸，前一个为变化的，后一个为变化后的
            map<int, int> mmMap;    //key为互补离子下标，value 为互补离子下标
            int comIonsSize = 0;  // 互补离子查找到的个数
//            makeTagByMonoMass((*prit), (*it), ccMap, tagSet);// 根据该条质谱得到tag 和转变后的ccMap
            //得到质谱互补离子对的map
            unknown::ptm_utils::getComplementIonsMap((*it)->precursor_mass_double, (*it)->ions_mass_container, mmMap);
            //判定蛋白中互补的占比
            //            double complementIonsInProtein = isHaveComplementIons((*prit),(*it),one_mod_ptr,mmMap,comIonsSize);
            //开始过滤
            for (; prit != prContainer.end(); ++prit) {
                count_s++;
                if ((*prit)->theoryMassC_By.size() == 0) {
                    continue;
                }
                //这里过滤主要考虑匹配离子的个数
                int zero_match_peaks = 0;   //完全匹配离子的个数
                int matchTagSize = 0;      //匹配上tag的个数
                int maxTagLen = 0;         //最长tag的长度
                double complementIonsInProtein = 0;
                bool logLeastScope;   //是否满足范围
                double leastSub = 0;   //范围内的差值
                if (prsmArg->BYIonsActions.count((*it)->activation_sub)) {       //如果是BY离子Action
                    double fMass = (*it)->precursor_mass_double - (*prit)->theoryMassC_By.back();
                    if (fMass > max_ptm_mass + prsmArg->unknownPTMsMassMax) {   //如果该蛋白质量太小
                        continue;
                    }
                    // 过滤第一步，根据tag ， 寻找到有几个匹配的tag（matchTagSize） 和 匹配上tag的最大长度（返回值）
                    maxTagLen = findBestProteinByTag((*prit), (*it), ccMap, tagSet, matchTagSize, fMass);
                    // 过滤第二步，完全匹配
                    zero_match_peaks = head_match_peaks_search((*it)->ions_mass_container, (*prit)->theoryMassC_By,
                                                               (*prit), (*it), (*prit)->one_ptm_seq_match_peaks);
                    //过滤第三步，质量添加
                    logLeastScope = selectMassLeast1000((*prit)->theoryMassC_By, (*it)->precursor_mass_double, (*prit)->leastSub,
                                                        prsmArg->unknownPTMsFilterProteinMassScope);
                }
                if (prsmArg->CZIonsActions.count((*it)->activation_sub)) {   //如果是CZ离子Action
                    double fMass = (*it)->precursor_mass_double - (*prit)->theoryMassC_Cz.back();
                    if (fMass > max_ptm_mass + prsmArg->unknownPTMsMassMax) {   //如果说该蛋白质量太小
                        continue;
                    }
                    // 根据tag ， 寻找到有几个匹配的tag（matchTagSize） 和 匹配上tag的最大长度（返回值）
                    maxTagLen = findBestProteinByTag((*prit), (*it), ccMap, tagSet, matchTagSize, fMass);
                    zero_match_peaks = head_match_peaks_search((*it)->ions_mass_container, (*prit)->theoryMassC_Cz,
                                                               (*prit), (*it), (*prit)->one_ptm_seq_match_peaks);
                    logLeastScope = selectMassLeast1000((*prit)->theoryMassC_Cz, (*it)->precursor_mass_double, (*prit)->leastSub,
                                                        prsmArg->unknownPTMsFilterProteinMassScope);
                }
                // 3 step 过滤保存结果
                //第一步过滤结果保存
                if (zero_match_peaks > prsmArg->unknownPTMsFilterCountSize) {
                    (*prit)->one_ptm_match_peaks = zero_match_peaks;
                    //                    (*prit)->complementIonsInProtein = isHaveComplementIons((*prit),(*it),one_mod_ptr,mmMap,comIonsSize);
                    px.push_back((*prit));
                }
                //第二步过滤结果保存
                if (maxTagLen > prsmArg->unknownPTMsFilterCountSize) {
                    (*prit)->matchTagSize = matchTagSize;  //匹配tag数
                    //                    (*prit)->complementIonsInProtein = isHaveComplementIons((*prit),(*it),one_mod_ptr,mmMap,comIonsSize);
                    (*prit)->maxTagLen = maxTagLen;
                    pTag.push_back((*prit));
                }
                //取第一步和第二步交叠
                if (zero_match_peaks > prsmArg->unknownPTMsFilterCountSize &&
                    maxTagLen > prsmArg->unknownPTMsFilterCountSize) {
                    pCom.push_back((*prit));
                }
                //第三步过滤结果保存
                if (logLeastScope) {
                    least1000.push_back((*prit));
                }
            }
            //选取top过滤蛋白
            int topSelect = prsmArg->unknownFilterOneAndTwoTopSelect;
            int leastTop = prsmArg->unknownFilterThreeTopSelect;
            sort(px.begin(), px.end(), sort_one_ptm_protein);   //根据 匹配数量排序
            sort(pTag.begin(), pTag.end(), sort_one_ptm_protein_2);
            sort(pCom.begin(), pCom.end(), sort_one_ptm_protein);
            sort(least1000.begin(), least1000.end(), sort_one_ptm_protein_4);
            // selectFilterProtein(px,pTag,p);
            set<protein_processor_ptr> setP;      //用于去重，并且确定最终的过滤蛋白
            for (int i = 0; i < px.size() && i < topSelect; ++i) {    //取 count top
                setP.insert(px[i]);
            }
            for (int i = 0; i < topSelect && i < pTag.size(); ++i) {   //取 tag top
                setP.insert(pTag[i]);
            }
            for (int i = 0; i < topSelect && i < pCom.size(); ++i) {   //取 pCom top 上两者
                setP.insert(pCom[i]);
            }
            for (int i = 0; i < least1000.size() && i < leastTop; ++i) {    //取 ms top
                setP.insert(least1000[i]);
            }
            for (set<protein_processor_ptr>::iterator it = setP.begin(); it != setP.end(); ++it) {  //
                p->prContainer.push_back((*it));
            }
//            cout << (*it)->Scans << " 过滤蛋白个数 : " << p->prContainer.size() << endl;
//            test_filter_result(it, p);   //测试输出蛋白质过滤结果
//             continue ;
//            one_ptm_filter_result.push_back(p);     //过滤蛋白指针用vector保存，
            //开始prsm
            vector<protein_processor_ptr>::iterator pr_it_begin = p->prContainer.begin();  //过滤蛋白库起始位置
            vector<protein_processor_ptr>::iterator pr_it_end = p->prContainer.end();     //过滤蛋白库终止位置
            for (; pr_it_begin != pr_it_end; ++pr_it_begin) {
                if (prsmArg->BYIonsActions.count(p->it->activation_sub)) {
//                    unknow_ptms_prsm_process_CID(one_mod_ptr, (*pr_it_begin), p->it, prsmArg, mmMap);    //修饰，蛋白，质谱
                }
                if (prsmArg->CZIonsActions.count(p->it->activation_sub)) {
//                    unknow_ptms_prsm_process_ETD(one_mod_ptr, (*pr_it_begin), p->it, prsmArg, mmMap);
                }
            }
            //输出prsm结果
            modify_ptr resultMPtr = p->it->unknownModPtr;
            writeResult(p->it, outFileName, resultMPtr);
//            delete resultMPtr;
        }   // ms end
    }


    void changeBestPrSMs(protein_processor_ptr &pro, msalign_processor_ptr &spact, temp_prsm_argument_ptr &arg)
    {
        spact->best_score = arg->score;
        spact->best_peak = arg->peaks_match;
        spact->match_protein_name = pro->protein_name;
        spact->cut_c_n = arg->cut_location_n;
        spact->cut_n_c = arg->cut_location_c;
        spact->best_ptm = arg->ptms;
        spact->match_protein_seq = pro->protein_sequence;
        spact->mut_s = arg->mut_x;
        spact->bestIonsPairNub = arg->pairIonsNub;
        spact->huBuCount = arg->complement_ions_number;
        spact->seqScore = arg->prsm_score;
        spact->seqScoreCount = arg->seq_score_count;
        spact->matchFragmentIonsSize = arg->matchFragmentIonsSize;
    }

    // 2023-1-7 更新将流程合并
    void mergeOnePTMPrSMsProcess(prsm::modify_ptr mp,
                                 protein_processor_ptr pro,
                                 msalign_processor_ptr spact,
                                 config_argument_ptr prsmArg,
                                 const string &variablePTMsFileOutPath)
    {
        //先给理论峰值,添加所有修饰峰值
        //最后定位蛋白质形
        //减前端截断，选取最为合适的蛋白质形分析
        int cut_n = -1, cut_c = 0;
        if (cut_n == -1)
        {
            if (prsmArg->BYIonsActions.count(spact->activation_sub)) {
                judge_best_n_cut_CID(mp, pro, spact, prsmArg, cut_n, cut_c);
            }
            if (prsmArg->CZIonsActions.count(spact->activation_sub)) {
                judge_best_n_cut_ETD(mp, pro, spact, prsmArg, cut_n, cut_c);
            }
        }

        if (cut_n == -1 && cut_c == 0)
        {
            return;
        }
        //遍历每一个修饰
        for (map<double, string>::iterator it_mp = mp->mapi.begin(); it_mp != mp->mapi.end(); ++it_mp)
        {
            temp_prsm_argument_ptr arg = new temp_prsm_argument();
            arg->cut_location_n = cut_n;
            arg->cut_location_c = cut_c;
            // 初始化理论峰序列
            if (prsmArg->BYIonsActions.count(spact->activation_sub)) {
                arg->theory_ions_mass_n.assign(pro->theoryMassN_By.begin() + cut_n, pro->theoryMassN_By.end() - cut_c);
                arg->theory_ions_mass_c.assign(pro->theoryMassC_By.begin() + cut_c, pro->theoryMassC_By.end() - cut_n);
            }
            if (prsmArg->CZIonsActions.count(spact->activation_sub)) {
                arg->theory_ions_mass_n.assign(pro->theoryMassN_Cz.begin() + cut_n, pro->theoryMassN_Cz.end() - cut_c);
                arg->theory_ions_mass_c.assign(pro->theoryMassC_Cz.begin() + cut_c, pro->theoryMassC_Cz.end() - cut_n);
            }
            // 减去截断质量
            initTheoryMass(arg->theory_ions_mass_n, arg->theory_ions_mass_c, pro, cut_n, cut_c);
            string precursorMassStr = spact->precursor_mass.substr(spact->precursor_mass.find("=") + 1);
            double prMass = atof(precursorMassStr.c_str());
            double unk_mass = 0.0 , bestUnknown = 0.0;
            // 计算未知修饰质量 2022-6-6 update
            if (prsmArg->BYIonsActions.count(spact->activation_sub) ){
                unk_mass = prMass - arg->theory_ions_mass_c.back() - it_mp->first;
//                bestUnknown = judge_Unknown_Mass(it_mp->first, prMass, unk_mass, arg, spact, pro);
            }
            if (prsmArg->CZIonsActions.count(spact->activation_sub) ) {
                unk_mass = prMass - arg->theory_ions_mass_c.back() - it_mp->first - prsmArg->sub_etd;
//                bestUnknown = judge_Unknown_Mass_ETD(it_mp->first, prMass, unk_mass, arg, spact, pro);
            }
            //探测最佳未知修饰质量
            bestUnknown = DetermineTheMostAppropriateUnknownModifierMass(it_mp->first, prMass, unk_mass, arg, spact, pro,prsmArg);
            if (fabs(unk_mass - bestUnknown) > prsmArg->scope_h) {
                continue;
            }
            vector<double> modifierTable;       //获取当前修饰所有质量
            map<double,string> mmap;           //double - 组合质量 , string - 组合修饰
            // 初始化mmap 和 modifierTable
            initial_Unknown_ptmTable(arg, modifierTable, mmap, it_mp->first, it_mp->second, bestUnknown);
            prsm::modify_ptr tmp = std::make_shared<prsm::Modify>(variablePTMsFileOutPath);   //这里参数无影响
            tmp->Init_modify();
            tmp->mapi = mmap;
            tmp->massStringToStringMass();
            sort(modifierTable.begin(), modifierTable.end());

            map<string,string> ssMap;
            ssMap["PTM_Mass"] = to_string(it_mp->first);
            ssMap["Unknown_PTM_Mass"] = to_string(bestUnknown) ;
            // PrSMs 鉴定流程
            prsm_process_unknown(tmp, pro, spact, arg, modifierTable,prsmArg);
            if ((arg->score - spact->best_score) > 0.000001 ) {
                // 交换公共部分的最佳PrSM result
                changeBestPrSMs(pro,spact,arg);
                spact->modify_mass = bestUnknown + it_mp->first;
                spact->unknownPTMMassMap = ssMap ;
//                delete spact->unknownModPtr;
                spact->unknownModPtr = tmp;
                if (prsmArg->BYIonsActions.count(spact->activation_sub)) {
                    spact->theory_proteoform_mass = spact->modify_mass + arg->theory_ions_mass_c.back();
                }
                if (prsmArg->CZIonsActions.count(spact->activation_sub)) {
                    spact->theory_proteoform_mass = spact->modify_mass + arg->theory_ions_mass_c.back() + prsmArg->sub_etd;
                }
            }
//            else {
//                delete tmp ;
//            }
            delete arg;
        }
    }

    //CID process  deprecate
    void unknow_ptms_prsm_process_CID(prsm::modify_ptr mp, protein_processor_ptr pro, msalign_processor_ptr spact, config_argument_ptr prsmArg,
                                      map<int, int> &iiMap, const string &variablePTMsFileOutPath) {
        //先给理论峰值,添加所有修饰峰值
        //最后定位蛋白质形
        //减前端截断，选取最为合适的蛋白质形分析
        int cut_n = -1, cut_c = 0;
        double mod = 0;
//        if (iiMap.size() > 0) {
//            judgeBestNCutSize(pro, spact, mp, iiMap, prsmArg, mod, cut_n);
//            cout << " 1 cut = " << cut_n << endl;
//        }
        if (cut_n == -1) {
            judge_best_n_cut_CID(mp, pro, spact, prsmArg, cut_n, cut_c);
        }
        if (cut_n == -1) {
            return;
        }
        int po = 0;
        vector<modify_ptr> mP;
        int best_mP = 0;
        int temp = 0;

        /**
         * 最外层循环是修饰个数的循环， 循环遍历
         * 对齐初始化工作  1、有N个修饰，那么循环N+1种情况，选出最合适的情况
         * 修饰 + 0 ， 还应加上未知的质量修饰
         */
        for (map<double, string>::iterator it_mp = mp->mapi.begin(); it_mp != mp->mapi.end(); ++it_mp, ++temp) {
            //初始化理论峰
            temp_prsm_argument_ptr arg = new temp_prsm_argument();
            arg->cut_location_n = cut_n;
            arg->cut_location_c = cut_c;
            arg->theory_ions_mass_n.assign(pro->theoryMassN_By.begin() + cut_n, pro->theoryMassN_By.end() - cut_c);
            arg->theory_ions_mass_c.assign(pro->theoryMassC_By.begin() + cut_c, pro->theoryMassC_By.end() - cut_n);
            initTheoryMass(arg->theory_ions_mass_n, arg->theory_ions_mass_c, pro, cut_n, cut_c);
            //推测最优未知修饰 2022-6-6 update
            string precursorMassStr= spact->precursor_mass.substr(spact->precursor_mass.find("=") + 1, spact->precursor_mass.size() + 1);
            double prMass = atof(precursorMassStr.c_str());
            double unk_mass = prMass - arg->theory_ions_mass_c.back() - it_mp->first;   //未知修饰质量
            double bestUnknown = judge_Unknown_Mass(it_mp->first, prMass, unk_mass, arg, spact, pro);
            if (fabs(unk_mass - bestUnknown) > prsmArg->scope_h) {
                continue;
            }
            vector<double> modifierTable;
            map<double, string> mmap;
            //初始化修饰Map
            initial_Unknown_ptmTable(arg, modifierTable, mmap, it_mp->first, it_mp->second, bestUnknown);
            modify_ptr tmp = std::make_shared<prsm::Modify>(variablePTMsFileOutPath);   //这里参数无影响
            tmp->mapi = mmap;
            tmp->massStringToStringMass();
            sort(modifierTable.begin(), modifierTable.end());
            //prsm 接口
            prsm_process_unknown(tmp, pro, spact, arg, modifierTable,prsmArg);
            if ((arg->score > spact->best_score)) {
                //交换公共部分的最佳PrSM
                changeBestPrSMs(pro,spact,arg);
                spact->theory_proteoform_mass = spact->modify_mass + arg->theory_ions_mass_c.back();
                spact->modify_mass = bestUnknown + it_mp->first;
                best_mP = temp;
                spact->unknownModPtr = tmp;
            }
            delete arg;
        }
        for (int i = 0; i < mP.size(); ++i) {
            if (i == best_mP) {
                continue;
            }
//            delete mP[i];
        }
        //map转化为vec
    }



    //ETD
    /**
     * deprecate
     * 鉴定主要流程:
     * 1、选出最为合适的蛋白质形，也就是进行截断的判定  判定接口 judge_best_n_cut
     * 2、初始化理论峰，进行最佳未知修饰质量的判定 judge_Unknown_Mass
     * @param mp
     * @param pro
     * @param spact
     */
    void unknow_ptms_prsm_process_ETD(modify_ptr mp, protein_processor_ptr pro, msalign_processor_ptr spact, config_argument_ptr prsmArg,
                                      map<int, int> &iiMap, const string &variablePTMsFileOutPath) {
        //先给理论峰值,添加所有修饰峰值
        //然后依次减前端截断，选取出最为合适的蛋白质形分析
        //最后定位蛋白质形
        int cut_n = -1, cut_c = 0;
        //减前端截断，选取出最为合适的蛋白质形分析
        double mod = 0;

        if (cut_n == -1) {
            judge_best_n_cut_ETD(mp, pro, spact, prsmArg, cut_n, cut_c);
        }

        if (cut_n == -1) {
            return;
        }
        int po = 0;
        vector<modify_ptr> mP;
        int best_mP = 0;
        int temp = 0;

        /**
         * 最外层循环是修饰个数的循环， 循环遍历
         * 对齐初始化工作  1、有N个修饰，那么循环N+1种情况，选出最合适的情况
         * 修饰 + 0 ， 还应加上未知的质量修饰
         */
         for (map<double, string>::iterator it_mp = mp->mapi.begin(); it_mp != mp->mapi.end(); ++it_mp, ++temp) {
            double subMassETD = prsmArg->sub_etd;//18.01056 - 1.9919;
            temp_prsm_argument_ptr arg = new temp_prsm_argument();
            arg->cut_location_n = cut_n;
            arg->cut_location_c = cut_c;
            arg->theory_ions_mass_n.assign(pro->theoryMassN_Cz.begin() + cut_n, pro->theoryMassN_Cz.end() - cut_c);
            arg->theory_ions_mass_c.assign(pro->theoryMassC_Cz.begin() + cut_c, pro->theoryMassC_Cz.end() - cut_n);
            //初始化理论质量
            initTheoryMass(arg->theory_ions_mass_n, arg->theory_ions_mass_c, pro, cut_n, cut_c);

            string precursorMassStr = spact->precursor_mass.substr(spact->precursor_mass.find("=") + 1, spact->precursor_mass.size() + 1);
            double prMass = atof(precursorMassStr.c_str());
            double unk_mass = prMass - arg->theory_ions_mass_c.back() - it_mp->first - subMassETD;  //未知修饰质量
            //求出最佳的未知修饰
            double bestUnknown = judge_Unknown_Mass_ETD(it_mp->first, prMass, unk_mass, arg, spact, pro);

            if (fabs(unk_mass - bestUnknown) > 1.5) {
                continue;
            }

            vector<double> modifierTable;   //修饰质量
            map<double, string> mmap; // 当前状态下的修饰map
            initial_Unknown_ptmTable(arg, modifierTable, mmap,it_mp->first, it_mp->second,bestUnknown);
            //初始化未知修饰  得到mmp参数
            modify_ptr tmp = std::make_shared<prsm::Modify>(variablePTMsFileOutPath);   //这里参数无影响
            tmp->mapi = mmap;
            tmp->massStringToStringMass();
            sort(modifierTable.begin(), modifierTable.end());
            //处理 Prsm 入口
            prsm_process_unknown(tmp, pro, spact, arg, modifierTable,prsmArg);
            if ((arg->score > spact->best_score)) {
                //交换公共部分的最佳PrSM
                changeBestPrSMs(pro,spact,arg);
                spact->modify_mass = bestUnknown + it_mp->first;
                spact->theory_proteoform_mass = spact->modify_mass + arg->theory_ions_mass_c.back() + 18.01056 - 1.9919;
                best_mP = temp;
                spact->unknownModPtr = tmp;
            }
            delete arg;
        }
        for (int i = 0; i < mP.size(); ++i) {
            if (i == best_mP) {
                continue;
            }
//            delete mP[i];
        }
    }

    /**
     * 鉴定prsm
     * @param modptr 修饰
     * @param pro   蛋白
     * @param spact 质谱
     * @param arg   临时参数
     * @param m     修饰质量序列
     */
    void prsm_process_unknown(prsm::modify_ptr modptr, protein_processor_ptr pro, msalign_processor_ptr spact, temp_prsm_argument_ptr arg, vector<double> &m,
                              config_argument_ptr prsmArg)
    {
        char p;
//        EX_TRACE("------------------CWDTW-alignment BEGIN------------------");

        mylib::data_stream::cwdtw_map(arg->theory_ions_mass_n, spact->ions_mass_container, arg->alignment_ions_n, arg->ref_score_n,
                                      arg->peer_score_n);
        mylib::data_stream::cwdtw_map(arg->theory_ions_mass_c, spact->ions_mass_container, arg->alignment_ions_c, arg->ref_score_c,
                                      arg->peer_score_c);

//         cout<<"------------------CWDTW-alignment END------------------"<<endl;
//         cout<<"------------------Scope_Align START------------------"<<endl;
        // 2022.9.19 更新  ， 获取所有可能的对齐数据
        mylib::data_stream::scopeAlignReWrite(arg->theory_ions_mass_n, spact->ions_mass_container, arg->alignment_ions_n, m, arg->var_mod_mass,
                                              prsmArg->ptmMassMax);
        mylib::data_stream::scopeAlignReWrite(arg->theory_ions_mass_c, spact->ions_mass_container, arg->alignment_ions_c, m, arg->var_mod_mass,
                                              prsmArg->ptmMassMax);
//        g::proc::Scope_Align(arg->thero_mass_C_N, spact->mono_mass,
//                             arg->alignment_C_N,
//                             m, arg->varible_mod_mass);
//        g::proc::Scope_Align(arg->thero_mass_N_C, spact->mono_mass,
//                             arg->alignment_N_C,
//                             m, arg->varible_mod_mass);
//        EX_TRACE("------------------Scope_Align END   ------------------");

        //没有做任何限定，单纯就是得出sub
        mylib::data_stream::Get_sub(arg->theory_ions_mass_n, arg->theory_ions_mass_c, spact->ions_mass_container, arg->sub_n, arg->sub_c,
                                    arg->alignment_ions_n, arg->alignment_ions_c);

        // 对所有的 sub 进行判定选取
        mylib::etd_process::indexIonsReWrite(arg->sub_n, arg->filter_ions_n, m,
                                             arg->var_mod_mass);      // (进行小部分数据过滤)(带有部分噪声离子数据)
        mylib::etd_process::indexIonsReWrite(arg->sub_c, arg->filter_ions_c, m,
                                             arg->var_mod_mass);      // (带有部分噪声离子数据)index 保存过滤数据

//        e::proc::index_ions(arg->sub_C_N, arg->index_C_N,
//                            m, arg->varible_mod_mass);      // (进行小部分数据过滤)(带有部分噪声离子数据)
//        e::proc::index_ions(arg->sub_N_C, arg->index_N_C,
//                            m,arg->varible_mod_mass);      // (带有部分噪声离子数据)index 保存过滤数据

//         cout<<"------------------GET    PPM    START------------------"<<endl;
        mylib::etd_process::GetPpmUnknownReWrite(arg->theory_ions_mass_n, spact->ions_mass_container, arg->alignment_ions_n, arg->filter_ions_n,
                                                 m, arg->ppm_value_n, arg->var_mod_mass);
        mylib::etd_process::GetPpmUnknownReWrite(arg->theory_ions_mass_c, spact->ions_mass_container, arg->alignment_ions_c, arg->filter_ions_c,
                                                 m, arg->ppm_value_c, arg->var_mod_mass);
//         cout<<"------------------GET    PPM    END------------------"<<endl;
        mylib::etd_process::Get_all_ions_peaks(arg->merge_ions_c, arg->alignment_ions_c, arg->filter_ions_c, arg->ppm_value_c,
                                               prsmArg->min_ppm_value);
        mylib::etd_process::Get_all_ions_peaks(arg->merge_ions_n, arg->alignment_ions_n, arg->filter_ions_n, arg->ppm_value_n,
                                               prsmArg->min_ppm_value);

        vector<node> N_ions_c_n;
        vector<node> N_ions_n_c;
        vector<modi> N_Ptms;
        vector<node> C_ions_c_n;
        vector<node> C_ions_n_c;
        vector<modi> C_Ptms;
        //定位未知修饰和已知修饰未知
        //C端选取结果
        mylib::etd_process::ions_location_match_N(N_ions_c_n, N_ions_n_c, arg->merge_ions_n, arg->merge_ions_c,
                                                  N_Ptms, m, arg->get_lenth(arg->theory_ions_mass_n.size()), arg->var_mod_mass);
        mylib::etd_process::ions_location_match_C(C_ions_c_n, C_ions_n_c, arg->merge_ions_n, arg->merge_ions_c,
                                                  C_Ptms, m, arg->get_lenth(arg->theory_ions_mass_n.size()), arg->var_mod_mass);

        if (arg->var_mod_mass != 0) {
            if (N_ions_n_c.size() + N_ions_c_n.size() < C_ions_n_c.size() + C_ions_c_n.size()) {
                arg->ions_n = C_ions_c_n;
                arg->ions_c = C_ions_n_c;
                arg->ptms = C_Ptms;
                p = 'C';
            } else {
                arg->ions_n = N_ions_c_n;
                arg->ions_c = N_ions_n_c;
                arg->ptms = N_Ptms;
                p = 'N';
            }
        } else if (arg->var_mod_mass == 0) {
            arg->ions_n = N_ions_c_n;
            arg->ions_c = C_ions_n_c;
            p = 'Z';
        }
        arg->get_match_peaks();
        if (arg->ions_n.size() == 0 && arg->ions_c.size() == 0) {
            return;
        }
//        test_out(p, arg, spact);      //测试输出
        //计算mut_x系数无用
        double mut_x = mylib::data_stream::analysis_all_ions(arg->ions_n, arg->ions_c, spact->ions_mass_container,
                                                             m, arg->var_mod_mass);
        arg->mut_x = mut_x;    //记录mut_x
        //计算互补离子个数
        double huBuIonsC = mylib::data_stream::getHuBuIonsScore(arg->ions_n, arg->ions_c, arg->theory_ions_mass_c.size()) + 1;
        arg->complement_ions_number = huBuIonsC;
        //计算连续的离子个数，返回的不是个数，而是经过计算的seqScore
        double seqScore = mylib::data_stream::getSeqIonsScore(arg->ions_n, arg->ions_c, arg->seq_score_count);
//        arg->seqScore = seqScore;   //返回的不是个数，而是经过计算的seqScore
//        9-22 得到连续离子对个数
        arg->pairIonsNub = unknown::ptm_utils::getSeqIonsPair(arg->ions_n, arg->ions_c);
        //2022-9-22 18:50 更改    10-20 当互补离子对为0 时， 得分也为 0
        arg->prsm_score = sqrt(((double) arg->pairIonsNub * 2) / (double) spact->ions_mass_container.size());    //计算连续离子得分
        arg->matchFragmentIonsSize += mylib::data_stream::get_fragment_ions_nub(arg->ions_n);    //得到匹配离子数
        arg->matchFragmentIonsSize += mylib::data_stream::get_fragment_ions_nub(arg->ions_c);    //得到匹配离子数
        //22-10-20 添加
        arg->prsm_score += ((double)arg->matchFragmentIonsSize / (double)spact->ions_mass_container.size());
        //2022-9-28 更改
//        double ionsSize = arg->ions_N_C.size() + arg->ions_C_N.size() ;
//        arg->seqScore =  sqrt( ((double)arg->pairIonsNub*2) / ionsSize ) + 0.1 ;
//        计算基础得分   10-13修正 添加了匹配碎片离子计算
        mylib::data_stream::Get_SCORE_Un(arg->ions_n, arg->ref_score_n, m, arg->peer_score_n, arg->score_n,
                                         0.1, arg->var_mod_mass);
        mylib::data_stream::Get_SCORE_Un(arg->ions_c, arg->ref_score_c, m, arg->peer_score_c, arg->score_c,
                                         0.1, arg->var_mod_mass);
        //验证ptm
        double mass = 0;
        for (int xi = 0; xi < arg->ptms.size(); ++xi) {
            mass += modptr->mapStrMass[modptr->analysis(arg->ptms[xi].mod_mass)];
        }
        arg->get_true_modMass(mass);
        if (fabs(mass - arg->var_mod_mass) > 3) {
            return;
        }
        if (mass == 0 && arg->var_mod_mass != 0) {
            return;
        }
        double ionsA = ((double) arg->ions_n.size() + arg->ions_c.size()) / (double) spact->ions_mass_container.size() *
                       10;    //计算峰值占比
                       //subCoe = /10 H2a
        int subCoe = (arg->cut_location_n / 10) + (arg->cut_location_c / 10); //计算截断总数
        arg->get_score(mut_x, huBuIonsC, subCoe, arg->prsm_score, ionsA); // 计算最终得分

        //测试输出函数
//        test_out(p, arg, spact);
    }


}

#include "unknown_ptm_search.h"
