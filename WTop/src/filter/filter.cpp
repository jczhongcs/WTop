//
// Created by wenzhong on 2023/3/21.
//

#include "filter.hpp"

bool match_ions_sort(filter_protein_ptr a, filter_protein_ptr b) {
    return a->getFilterMatchIons() > b->getFilterMatchIons();
}

bool match_tag_length_sort(filter_protein_ptr a, filter_protein_ptr b) {
    return a->filter_tag_length > b->filter_tag_length;
}

bool massErrorSort(filter_protein_ptr a, filter_protein_ptr b) {
    return a->tagMassError < b->tagMassError;
}

bool one_ptm_search_peaks_sort(filter_protein_ptr a, filter_protein_ptr b) {
    return a->getOnePtmMatchPeaks() > b->getOnePtmMatchPeaks();
}

bool mass_gap_sort(filter_protein_ptr a, filter_protein_ptr b) {
    return a->getMassGap() < b->getMassGap();
}

void selectTopCandidateProtein(const msalign_ptr &msalignPtr,
                               vector<filter_protein_ptr> &filterProPtrVec,
                               int topSelect) {
    vector<filter_protein_ptr> finalProPtrVec;
    for (int i = 0; i < filterProPtrVec.size(); i++) {
        if (finalProPtrVec.size() > topSelect) {
            break;
        }
        vector<filter_protein_ptr> tempVec;
        tempVec.push_back(filterProPtrVec[i]);

        for (int j = i + 1; j < filterProPtrVec.size(); ++j) {
            if (filterProPtrVec[j]->filter_tag_length == filterProPtrVec[i]->filter_tag_length) {
                tempVec.push_back(filterProPtrVec[j]);
            } else {
                break;
            }
        }
        sort(tempVec.begin(), tempVec.end(), massErrorSort);
        for (int j = 0; j < tempVec.size(); ++j) {
            if (finalProPtrVec.size() > topSelect) {
                break;
            }
            finalProPtrVec.push_back(tempVec[j]);
        }
        i += tempVec.size() - 1;
    }

    //insert in to candidate protein
    for (filter_protein_ptr te : finalProPtrVec) {
        msalignPtr->insert_protein_ptr(te->proteinPtr);
    }
}

void filter::filtrationProcess(const vector<msalign_ptr> &msContainer,
                               const vector<protein_ptr> &prContainer,
                               const prsm::modify_ptr modifyPtr,
                               const prsm::config_argument_ptr configArgumentPtr) {
    int count = 1;
    map<char, double> acidMap;
    msalign_tag_ptr msalignTagPtr = new msalign_tag();
    if (prContainer.size() > 0) {
        acidMap = prContainer[0]->get_acid_Map();
    } else {
        cout << "protein database error, filter candidate protein failed." << endl;
        exit(10);
    }

    //    ptm_filter_ptr ptmFilterPtr = new ptm_filter();
    //    vector<filter_protein> filterProteinVec =
    //            ptmFilterPtr->ptmFilterProcess(msContainer,prContainer,modifyPtr,configArgumentPtr);
    //    delete ptmFilterPtr ;



//    const std::string ion_file_name = "ion_filter_result.txt";
//    const std::string tag_file_name = "tag_filter_result.txt";
    // 一次性打开文件，保持打开状态（ios::app 确保追加写入，文件不存在则自动创建）
//    std::ofstream ion_file(ion_file_name, std::ios::out | std::ios::app);
//    std::ofstream tag_file(tag_file_name, std::ios::out | std::ios::app);



    for (vector<msalign_ptr>::const_iterator ms_it = msContainer.begin(); ms_it != msContainer.end(); ++ms_it, ++count) {
        std::cout << std::flush << "\rCurrent thread id : " << std::this_thread::get_id() << " . "
                  << "Number of Ms Filter : " << count << " th ";

        if ((*ms_it)->getIonsMassContainer().size() < configArgumentPtr->minMonoPeaks
            || (*ms_it)->getPrecursorMass() < configArgumentPtr->min_precursor_mass) {
            continue;
        }

        vector<filter_protein_ptr> match_ions_filters; // 保存由zero 匹配的蛋白，根据匹配个数，最后取top
        vector<filter_protein_ptr> match_tag_filters;  // 保存由tag 匹配的蛋白，最后取top
        //        vector<filter_protein_ptr> complement_ions_filter;  //保存由互补离子 匹配的蛋白，最后取top

        vector<filter_protein_ptr> overlap_ions_tag_filters;
        vector<filter_protein_ptr> precursor_mass_filters; // 保存前提质量 匹配的蛋白，最后取top
        set<filter_protein_ptr> final_candidate_protein;

        set<string> tagSet;    // store ms tag
        map<char, char> ccMap; // 保存变化过得氨基酸，前一个为变化的，后一个为变化后的
        map<int, int> mmMap;   // complement ions index

        if (prContainer.size() == 0) {
            cout << " error protein database " << endl;
            exit(7);
        }
        // get tag
        msalignTagPtr->getTagByMonoMass((*ms_it)->getIonsMassContainer(),
                                        acidMap, tagSet, ccMap);

        complement_ions_ptr complementIonsPtr = new complement_ions();
        // get complement ions index mmMap
        complementIonsPtr->getComplementIonsMap((*ms_it));

        // 判定蛋白中互补的占比
        //             double complementIonsInProtein = isHaveComplementIons((*prit),(*it),one_mod_ptr,mmMap,comIonsSize);
        // 开始过滤
        for (vector<protein_ptr>::const_iterator prit = prContainer.begin(); prit != prContainer.end(); ++prit) {
            filter_protein_ptr filterProteinPtr = new filter_protein((*prit));

            if ((*prit)->getTheoryMassCBy().size() == 0 || (*prit)->getTheoryMassCCz().size() == 0) {
                continue;
            }

            if (configArgumentPtr->BYIonsActions.count((*ms_it)->getActivationSub())) {
                double precursor_mass_gap = (*ms_it)->getPrecursorMass() - (*prit)->getTheoryMassCBy().back();
                if (precursor_mass_gap > (configArgumentPtr->ptmMassMax + configArgumentPtr->unknownPTMsMassMax)) {
                    continue;
                }
                ms_protein_filter_ptr filterPtr = new ms_protein_filter();
                filterProteinPtr->setFilterTagLength(filterPtr->findBestProteinByTag((*prit), (*ms_it), ccMap, tagSet, precursor_mass_gap));

                int one_ptm_match_peaks = 0;
                filterProteinPtr->setFilterMatchIons(
                    filterPtr->zero_match_peaks_search((*ms_it)->getIonsMassContainer(),
                                                       (*prit)->getTheoryMassCBy(),
                                                       one_ptm_match_peaks));
                filterProteinPtr->setOnePtmMatchPeaks(one_ptm_match_peaks);

                double mass_gap = 0.0;
                filterProteinPtr->setLogLeastScope(
                    filterPtr->precursor_mass_scope_search((*prit)->getTheoryMassCBy(),
                                                           mass_gap, (*ms_it)->getPrecursorMass(),
                                                           configArgumentPtr->unknownPTMsFilterProteinMassScope));
                filterProteinPtr->setMassGap(mass_gap);
                delete filterPtr;
            }

            if (configArgumentPtr->CZIonsActions.count((*ms_it)->getActivationSub())) { // 如果是CZ离子Action
                double precursor_mass_gap = (*ms_it)->getPrecursorMass() - (*prit)->getTheoryMassCCz().back();
                if (precursor_mass_gap > configArgumentPtr->ptmMassMax + configArgumentPtr->unknownPTMsMassMax) { // 如果说该蛋白质量太小
                    continue;
                }

                ms_protein_filter_ptr filterPtr = new ms_protein_filter();
                filterProteinPtr->setFilterTagLength(
                    filterPtr->findBestProteinByTag((*prit), (*ms_it), ccMap,
                                                    tagSet, precursor_mass_gap));

                int one_ptm_match_peaks = 0;
                filterProteinPtr->setFilterMatchIons(
                    filterPtr->zero_match_peaks_search((*ms_it)->getIonsMassContainer(),
                                                       (*prit)->getTheoryMassCCz(),
                                                       one_ptm_match_peaks));
                filterProteinPtr->setOnePtmMatchPeaks(one_ptm_match_peaks);

                double mass_gap = 0.0;
                filterProteinPtr->setLogLeastScope(
                    filterPtr->precursor_mass_scope_search((*prit)->getTheoryMassCCz(),
                                                           mass_gap, (*ms_it)->getPrecursorMass(),
                                                           configArgumentPtr->unknownPTMsFilterProteinMassScope));
                filterProteinPtr->setMassGap(mass_gap);

                delete filterPtr;
            }

            // 第一步过滤结果保存
            bool whether_release = true;
            if (filterProteinPtr->getFilterMatchIons() > configArgumentPtr->unknownPTMsFilterCountSize) {
                whether_release = false;
                match_ions_filters.push_back(filterProteinPtr);
            }
            // 第二步过滤结果保存
            if (filterProteinPtr->getFilterTagLength() > configArgumentPtr->unknownPTMsFilterCountSize) {
                whether_release = false;
                match_tag_filters.push_back(filterProteinPtr);
            }
            // 取第一步和第二步交叠
            if (filterProteinPtr->getFilterMatchIons() > configArgumentPtr->unknownPTMsFilterCountSize && filterProteinPtr->getFilterTagLength() > configArgumentPtr->unknownPTMsFilterCountSize) {
                whether_release = false;
                overlap_ions_tag_filters.push_back(filterProteinPtr);
            }
            // 第三步过滤结果保存
            if (filterProteinPtr->isLogLeastScope()) {
                whether_release = false;
                precursor_mass_filters.push_back(filterProteinPtr);
            }
            if (whether_release) {
                delete filterProteinPtr;
            }
        } // for protein final

        // 选取top过滤蛋白
        int topSelect = configArgumentPtr->unknownFilterOneAndTwoTopSelect;
        int leastTop = configArgumentPtr->unknownFilterThreeTopSelect;

        sort(match_ions_filters.begin(), match_ions_filters.end(), match_ions_sort); // 根据 匹配数量排序
        sort(match_tag_filters.begin(), match_tag_filters.end(), match_tag_length_sort);
        sort(overlap_ions_tag_filters.begin(), overlap_ions_tag_filters.end(), one_ptm_search_peaks_sort);
        sort(precursor_mass_filters.begin(), precursor_mass_filters.end(), mass_gap_sort);

        // ofstream out("match_ions_filter.txt", ios::out);
        // for (int i = 0; i < match_ions_filters.size(); ++i) { // 取 count top
        //     filter_protein_ptr filterProteinPtr = match_ions_filters[i];
        //     out << filterProteinPtr->getProteinPtr()->getProteinTitle()
        //         << " ions : " << filterProteinPtr->getFilterMatchIons() << endl;
        // }
        // ofstream out2("match_tag_filter.txt", ios::out);
        // for (filter_protein_ptr filterProteinPtr : match_tag_filters) {
        //     out2 << filterProteinPtr->getProteinPtr()->getProteinTitle()
        //          << " tag : " << filterProteinPtr->getFilterTagLength() << endl;
        // }
        // ofstream out3("match_tag_ions_filter.txt", ios::out);
        // for (filter_protein_ptr filterProteinPtr : overlap_ions_tag_filters) {
        //     out3 << filterProteinPtr->getProteinPtr()->getProteinTitle() << endl;
        // }

        for (int i = 0; i < topSelect && i < match_ions_filters.size(); ++i) { // 取 count top
            final_candidate_protein.insert(match_ions_filters[i]);
            size_t equal_pos = (*ms_it)->getId().find('=');
            std::string id_num_str = (*ms_it)->getId().substr(equal_pos + 1);
//            ion_file<<id_num_str<<" "<<match_ions_filters[i]->getProteinPtr()->getProteinTitle()<<endl;

        }
        for (int i = 0; i < topSelect && i < match_tag_filters.size(); ++i) { // 取 tag top
            final_candidate_protein.insert(match_tag_filters[i]);
            size_t equal_pos = (*ms_it)->getId().find('=');
            std::string id_num_str = (*ms_it)->getId().substr(equal_pos + 1);
//            tag_file<<id_num_str<<" "<<match_tag_filters[i]->getProteinPtr()->getProteinTitle()<<endl;
        }
        for (int i = 0; i < topSelect && i < overlap_ions_tag_filters.size(); ++i) { // 取 pCom top 上两者
            final_candidate_protein.insert(overlap_ions_tag_filters[i]);
        }
        for (int i = 0; i < precursor_mass_filters.size() && i < leastTop; ++i) { // 取 ms top
            final_candidate_protein.insert(precursor_mass_filters[i]);
        }

        //        get candidate protein
        for (filter_protein_ptr filterProteinPtr : final_candidate_protein) {
            (*ms_it)->insert_protein_ptr(filterProteinPtr->getProteinPtr());
        }

        // free filterProteinPtr
        set<filter_protein_ptr> all_filter_protein_ptr;
        for (int i = 0; i < match_ions_filters.size(); ++i) { // 取 count top
            all_filter_protein_ptr.insert(match_ions_filters[i]);
        }
        for (int i = 0; i < match_tag_filters.size(); ++i) { // 取 tag top
            all_filter_protein_ptr.insert(match_tag_filters[i]);
        }
        for (int i = 0; i < overlap_ions_tag_filters.size(); ++i) { // 取 pCom top 上两者
            all_filter_protein_ptr.insert(overlap_ions_tag_filters[i]);
        }
        for (int i = 0; i < precursor_mass_filters.size(); ++i) { // 取 ms top
            all_filter_protein_ptr.insert(precursor_mass_filters[i]);
        }

        for (filter_protein_ptr filterProteinPtr : all_filter_protein_ptr) {
            delete filterProteinPtr;
        }
        delete complementIonsPtr;

        vector<filter_protein_ptr>().swap(match_ions_filters);
        vector<filter_protein_ptr>().swap(match_tag_filters);
        vector<filter_protein_ptr>().swap(overlap_ions_tag_filters);
        vector<filter_protein_ptr>().swap(precursor_mass_filters);
        set<filter_protein_ptr>().swap(all_filter_protein_ptr);
        set<filter_protein_ptr>().swap(final_candidate_protein);

        //        cout << (*ms_it)->getScans() << " candidate protein number : " << (*ms_it)->getCandidateProteins().size() << endl;
    }

//    ion_file.flush();
//    ion_file.close();
//    tag_file.flush();
//    tag_file.close();

    delete msalignTagPtr;
}


void filter::getToppicFilterResult(const vector<msalign_ptr> &msContainer,
                               const vector<protein_ptr> &prContainer,
                               const prsm::modify_ptr modifyPtr,
                               const prsm::config_argument_ptr configArgumentPtr) {



    map<string,set<string>> flResultMap;
    ifstream in("total_result.txt", ios::in);
    string buff = " ";
    int i = 0 ;
        while(in.good()){
            getline(in,buff) ;
            if (!buff.empty() && buff.back() == 13) {
                buff.pop_back();
            }
            if (buff.find(" ") != buff.npos)
            {
                vector<string> str;
                utils_string_util::Stringsplit(buff,' ',str) ;

                if (str.size() != 2) {
                    cout<<"tempstr.size=="<<str.size()<<endl;
                    cout << " toppic result read error ! " << endl ; 
                    exit(1);
                }
                if (flResultMap.count(str[0])) {
                    flResultMap[str[0]].insert(str[1]);
                } else {
                    set<string> tempSet ; 
                    tempSet.insert(str[1]);
                    flResultMap.insert(make_pair(str[0],tempSet));
                }
            }
        }
    // cout<<endl;
    // for(auto pair : flResultMap) {
    //     cout << "id : " << pair.first << endl; 
    //     for (auto proName : pair.second) {
    //         cout <<"proName:" <<proName << endl ; 
    //     }
    // }

    int count = 1;
//    for (vector<msalign_ptr>::const_iterator ms_it = msContainer.begin(); ms_it != msContainer.end(); ++ms_it, ++count) {
//        for(auto pair:flResultMap){
//            std::string ms_it_Id= (*ms_it)->getId();
//            //cout<<endl<<"ms_it_ID="<<ms_it_Id<<"###"<<"pair.first="<<"ID="+pair.first<<endl;
//            if(strcmp((*ms_it)->getId().c_str(),("ID="+pair.first).c_str()) == 0){
//              //  cout<<endl<<"ms_it_ID="<<ms_it_Id<<"###"<<"pair.first="<<"ID="+pair.first<<endl;
//                for(auto proName:pair.second){
//                    for(auto pr_it : prContainer){
//
//                        string str = pr_it->getProteinName();
//                        string wtopProName = utils_string_util::cutProName(str);
//                        string pName = ">" +proName;
//                //        cout<<endl<<"wtop_pro_name="<<wtopProName<<endl<<"proName="<<pName<<endl;
//                        if (strcmp(wtopProName.c_str(),pName.c_str())==0)
//                        {
//                  //          cout<<endl<<"wtop_pro_name="<<wtopProName<<endl<<"proName="<<pName<<endl;
//                            (*ms_it)->insert_protein_ptr(pr_it);
//                        }
//                    }
//                }
//            }
//        }
//    }

    map<string,protein_ptr> getProPtrMap;
    for(protein_ptr pro_ptr :prContainer){
        getProPtrMap.insert(make_pair(utils_string_util::cutProName(pro_ptr->getProteinName()),pro_ptr));
    }


    for (vector<msalign_ptr>::const_iterator ms_it = msContainer.begin(); ms_it != msContainer.end(); ++ms_it, ++count) {
        string cutID =  utils_string_util::getProID((*ms_it)->getId());

        set<string> proset = flResultMap[cutID];
        //cout<<endl<<"###(*ms_it)->getId()="<<(*ms_it)->getId()<<"####proset.size()="<<proset.size()<<"#######"<<endl;


        if(proset.size()==0) {
            continue;
        } else{
           // cout<<endl<<"ms_it_ID="<<(*ms_it)->getId()<<"###"<<"pair.first="<<"ID="+(*ms_it)->getId()<<endl;
            for(auto proName:proset){

                string pName = ">" +proName;
                    //cout<<endl<<"wtop_pro_name="<<wtopProName<<endl<<"proName="<<pName<<endl;
                    (*ms_it)->insert_protein_ptr(getProPtrMap[pName]);
            }

        }
    }
}

















/**
 * tag推测
 * @param msContainer
 * @param prContainer
 * @param modifyPtr
 * @param configArgumentPtr
 */
void filter::tagFilterProcess(const vector<msalign_ptr> &msContainer,
                              const vector<protein_ptr> &prContainer,
                              const prsm::modify_ptr modifyPtr,
                              const prsm::config_argument_ptr configArgumentPtr) {
    int count = 0;
    map<char, double> acidMap; // acid

    if (prContainer.size() > 0) {
        acidMap = prContainer[0]->get_acid_Map();
    } else {
        cout << "protein database error, filter candidate protein failed." << endl;
        exit(10);
    }

    ofstream outt("get_tag_time.txt", ios::app);
    for (vector<msalign_ptr>::const_iterator msIt = msContainer.begin(); msIt != msContainer.end(); ++msIt) {
        std::cout << std::flush << " \r "
                  << "Scans : " << (*msIt)->getScans()
                  << " Current Thread Id : " << std::this_thread::get_id()
                  << " -> Number Of Ms Filter : " << ++count << " th ";

        if ((*msIt)->getIonsMassContainer().size() < configArgumentPtr->minMonoPeaks
            || (*msIt)->getPrecursorMass() < configArgumentPtr->min_precursor_mass) {
            continue;
        }
        auto start = system_clock::now();                                  // 记录开始时间
        msalign_tag_ptr msalignTagPtr = new msalign_tag((*msIt), acidMap); // ms align tag processor
        auto end = system_clock::now();                                    // 记录开始时间
        auto duration = duration_cast<milliseconds>(end - start);          // 计算时间差
        outt << "Scans : " << (*msIt)->getScans() << " ==>"
             << "Time taken by function: " << duration.count() << " milliseconds" << endl;
        continue;
        complement_ions_ptr complementIonsPtr = new complement_ions();      // complementIons Ptr
        ms_protein_filter_ptr msProteinFilterPtr = new ms_protein_filter(); // filter processer
        vector<filter_protein_ptr> filterProPtrVec;                         // final filter protein ptr
        for (vector<protein_ptr>::const_iterator prIt = prContainer.begin(); prIt != prContainer.end(); ++prIt) {
            filter_protein_ptr filterProteinPtr = new filter_protein(*(prIt));
            // msProteinFilterPtr->tagsMatchSequence((*prIt), (*msIt), complementIonsPtr, msalignTagPtr, filterProteinPtr);
            if (filterProteinPtr->tagMassError < 0.1
                && (filterProteinPtr->tagN.length() != 0 || filterProteinPtr->tagC.length() != 0)) {
                filterProPtrVec.push_back(filterProteinPtr);
            } else {
                delete filterProteinPtr;
            }
        }

        sort(filterProPtrVec.begin(), filterProPtrVec.end(), match_tag_length_sort);
        selectTopCandidateProtein((*msIt), filterProPtrVec, configArgumentPtr->unknownFilterOneAndTwoTopSelect);
        // free malloc
        for (filter_protein_ptr te : filterProPtrVec) {
            delete te;
        }

        delete msProteinFilterPtr;
        delete complementIonsPtr;
        delete msalignTagPtr;
    }
}
