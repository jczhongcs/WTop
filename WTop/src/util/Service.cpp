//
// Created by 碎雨粘霓裳 on 2021/12/7.
//

#include "Service.h"

namespace mylib {

    bool cmpMass(mylib::Argument a, mylib::Argument b) {
        return a.pPrecursorMass < b.pPrecursorMass;
    }

    void Service::MassSortArgument(vector<mylib::Argument> &aContainer) {
        sort(aContainer.begin(), aContainer.end(), cmpMass);
    }

    bool cmpScore(mylib::Argument a, mylib::Argument b) {
        return atof(a.score.c_str()) > atof(b.score.c_str());
    }

    void Service::ScoreSortArgument(vector<mylib::Argument> &aContainer) {
        sort(aContainer.begin(), aContainer.end(), cmpScore);
    }

    //计算fdr
    int Service::calculateFdr(vector<mylib::Argument> &aContainer) {
        double fdrNub = 0;
        double p = 1;
        for (int xi = 0; xi < aContainer.size(); xi++, p++) {
            if (aContainer[xi].proteinAccession.find(">DECOY") != aContainer[xi].proteinAccession.npos) {
                fdrNub++;
            }
            if (fdrNub >= 1) {
                aContainer[xi].fdrScore = (fdrNub / p);
            } else {
                aContainer[xi].fdrScore = 0;
            }
        }
        return fdrNub;
    }

    //去掉只有1条Prsm的，fdr>vakue,decoy,Prsm的个数大于 minNub.
    void Service::fdrChoice(vector<mylib::Argument> &aContainer, vector<mylib::Argument> &bContainer, double value,
                            int minNub, double filterArg) {
        std::map<string, int> mapi;
        vector<mylib::Argument> cContainer;
        vector<mylib::Argument> resultContainer;
        //去掉fdr >0.01 ,Decoy    用cContainer
        for (int xi = 0; xi < aContainer.size(); xi++) {
            if (aContainer[xi].fdrScore > value)
                continue;
            if (aContainer[xi].proteinAccession.find(">DECOY") != aContainer[xi].proteinAccession.npos) {
                continue;
            }
            cContainer.push_back(aContainer[xi]);
        }

        for (int xi = 0; xi < cContainer.size(); xi++) {
            if (!mapi.count(cContainer[xi].proteinAccession)) {
                mapi[cContainer[xi].proteinAccession] = 1;
            } else {
                mapi[cContainer[xi].proteinAccession]++;
            }
        }

        for (int xi = 0; xi < cContainer.size(); xi++) {
            if (mapi[cContainer[xi].proteinAccession] <= minNub) {
                continue;
            }
            bContainer.push_back(cContainer[xi]);
        }
        double loop_1 = 0;   // 1条Prsm
        double all_Prsm_Size = bContainer.size();
        double loop_2 = 0; // 2条Prsm
        double loop_3 = 0; // 3条Prsm
        double loop_4 = 0; // 3条Prsm
        vector<string> killProtein;
//        cout << "过滤prsm = 1 之前的 Prsm Number = " << all_Prsm_Size << endl;
        //string 为蛋白，int为Prsm的个数
        for (map<string, int>::iterator it = mapi.begin(); it != mapi.end(); ++it) {
//            cout<<"Protein = "<<it->first<<" size = "<<it->second<<endl;
            if (it->second > 3) {
                loop_4 += it->second;
                continue;
            }
            if (it->second == 1) {
                killProtein.push_back(it->first);
                loop_1 += it->second;
            }
            if (it->second == 2) {
                loop_2 += it->second;
            }
            if (it->second == 3) {
                loop_3 += it->second;
            }
        }
//        cout << "Prsm = 1 ,size = " << loop_1 << " 占比" << loop_1 / all_Prsm_Size * 100 << " % " << endl;
//        cout << "Prsm = 2 ,size = " << loop_2 << " 占比" << loop_2 / all_Prsm_Size * 100 << " % " << endl;
//        cout << "Prsm = 3 ,size = " << loop_3 << " 占比" << loop_3 / all_Prsm_Size * 100 << " % " << endl;
//        cout << "Prsm > 3 ,size = " << loop_4 << " 占比" << loop_4 / all_Prsm_Size * 100 << " % " << endl;

        //得到最终Prsm
        if ((loop_1 / all_Prsm_Size) <= filterArg) {
            for (int i = 0; i < killProtein.size(); ++i) {
                mapi[killProtein[i]] = 0;
            }
        }
        for (int xi = 0; xi < cContainer.size(); xi++) {
            if (mapi[cContainer[xi].proteinAccession] <= minNub) {
                continue;
            }
            resultContainer.push_back(cContainer[xi]);
        }
        bContainer = resultContainer;
//        cout << "Prsm Number = " << bContainer.size() << endl;
    }

    /**
     * 处理已知修饰的Protenform
     * @param bContainer fdr 之后的PrSMs
     * @param cContainer protenform
     */
    void Service::getProteinForm(vector<mylib::Argument> &bContainer, vector<mylib::Argument> &cContainer) {
        multimap<string,mylib::Argument> mmp;
        int count = 0 ;

        for (int i = 0; i < bContainer.size(); i++)
        {
            if (i == 0 && bContainer.size() != 0)
            {
                mmp.insert(make_pair(bContainer[0].proteinAccession,bContainer[0]));
                continue;
            }
            mylib::Argument comp = bContainer[i] ;
            count++;
            bool boolean = false ;
            for (auto iter = mmp.begin(); iter != mmp.end(); ++iter)
            {
                mylib::Argument refer = iter->second ;
                if (comp.proteinAccession == refer.proteinAccession)
                {
                    if (comp.firstResidue == refer.firstResidue && comp.lastResidue == refer.lastResidue )
                    {
                        if (fabs(comp.pPrecursorMass - refer.pPrecursorMass)<0.1)
                        {
                            if (fabs(comp.trueModMass - refer.trueModMass) < 1)
                            {
                                boolean = true;
                                break ;
                            }
                        }
                    }
                }
            }
            if (boolean == false)
            {
                mmp.insert(make_pair(comp.proteinAccession,comp));
            }
        }
        for (auto iter = mmp.begin(); iter != mmp.end(); ++iter)
        {
            cContainer.push_back(iter->second);
        }
        cout << cContainer.size()<<endl;
    }

    
    void Service::getUnknownProteinForm(vector<mylib::Argument> &bContainer, vector<mylib::Argument> &cContainer) {
        multimap<string,mylib::Argument> mmp;
        int count = 0 ;

        for (int i = 0; i < bContainer.size(); i++)
        {
            if (i == 0 && bContainer.size() != 0)
            {
                mmp.insert(make_pair(bContainer[0].proteinAccession,bContainer[0]));
                continue;
            }
            mylib::Argument comp = bContainer[i] ;
            count++;
            bool boolean = false ;
            for (auto iter = mmp.begin(); iter != mmp.end(); ++iter)
            {
                mylib::Argument refer = iter->second ;
                if (comp.proteinAccession == refer.proteinAccession)
                {
                    if (comp.firstResidue == refer.firstResidue && comp.lastResidue == refer.lastResidue )
                    {
                        if (fabs(comp.pPrecursorMass - refer.pPrecursorMass)<1)
                        {
                            double compPTMMass = stod(comp.ptmMass);
                            double compUnknownPTMMass = stod(comp.unknownPtmMass);
                            double referPTMMass = stod(refer.ptmMass);
                            double referUnknownPTMMass = stod(refer.unknownPtmMass);
                            if (fabs(compPTMMass - referPTMMass) < 1 && fabs(compUnknownPTMMass- referUnknownPTMMass)<1)
                            {
                                boolean = true;
                                break ;
                            }
                        }
                    }
                }
            }
            if (boolean == false)
            {
                mmp.insert(make_pair(comp.proteinAccession,comp));
            }
        }
        for (auto iter = mmp.begin(); iter != mmp.end(); ++iter)
        {
            cContainer.push_back(iter->second);
        }
        cout << cContainer.size()<<endl;
    }


    map<string,string> Service::getScansScoreMap(vector<mylib::Argument> &Container, map<string,mylib::Argument> &m)
    {
        map<string,string> result;
        for (auto wrapIter = Container.begin(); wrapIter != Container.end(); ++wrapIter) {
            result.insert(make_pair(wrapIter->Scans,wrapIter->score));
            m.insert(make_pair(wrapIter->Scans,(*wrapIter)));
        }
        return result;
    }
}
