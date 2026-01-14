//
// Created by Administrator on 2022/9/22.
//

#include "unknown_ptm_utils.h"


/**
 * 得到选取离子中，离子对的个数
 * @param t
 * @param m
 * @return
 */
int unknown::ptm_utils::getSeqIonsPair(vector<node> & t, vector<node> & m)
{
    int pairSize = 0 ;
    for (int i = 1 ; i < t.size(); ++i) {
        if (t[i].thoe_id - t[i-1].thoe_id == 1) {
            pairSize++;
            i += 1 ;
        }
    }

    for (int i = 1 ; i < m.size(); ++i) {
        if (m[i].thoe_id - m[i-1].thoe_id == 1) {
            pairSize++;
            i += 1 ;
        }
    }
    return pairSize ;
}

//质量相近的氨基酸进行统一替换
void replaceCMap(map<char,double> & cMass,map<char,char> &rMap,map<int,char> &iCMap)
{
    vector<double> pMass ;
    //先将质量相同的替换
    map<double,int> newCMap;
    for( auto it = cMass.begin();it != cMass.end(); ++it) {
        if (!newCMap.count(it->second))
        newCMap.insert(make_pair(it->second,it->first));
        else {
            char c2 = newCMap[it->second];
            rMap.insert(make_pair(it->first,c2));
        }
    }
//    for( auto it = newCMap.begin();it != newCMap.end(); ++it) {
//        cout <<" first = " << it->first<< " second = " << it->second <<endl ;
//    }
//    cout  <<"=================================" << endl;

    for(auto it = newCMap.begin() ; it != newCMap.end(); ++it) {
        int temp = it->first ;
        char c = it->second ;
        if (iCMap.count(temp) || iCMap.count(temp-1) || iCMap.count(temp+1)) {
           char c2 ;
           if(iCMap.count(temp)) {
               c2 = iCMap[temp];
           } else if (iCMap.count(temp-1)){
               c2 = iCMap[temp-1];
           } else {
               c2 = iCMap[temp+1];
           }
           rMap.insert(make_pair(c,c2));
            if (!iCMap.count(temp)){
                iCMap.insert(make_pair(temp,c2));
            }
           if (!iCMap.count(temp-1)){
               iCMap.insert(make_pair(temp-1,c2));
           }
           if (!iCMap.count(temp+1)){
                iCMap.insert(make_pair(temp+1,c2));
           }
        } else {
            iCMap.insert(make_pair(temp-1,c));
            iCMap.insert(make_pair(temp,c));
            iCMap.insert(make_pair(temp+1,c));
        }
    }
    //测试输出
//    for( auto it = rMap.begin();it != rMap.end(); ++it) {
//        cout <<" first = " << it->first<< " second = " << it->second <<endl ;
//    }
//    cout <<"=================================" << endl;
//    for( auto it = iCMap.begin();it != iCMap.end(); ++it) {
//        cout <<" first = " << it->first<< " second = " << it->second <<endl ;
//    }
}
/**
 * 根据mono的质量寻找tag ，并保存到set中
 * @param mono
 * @param massChar
 */
void getTagByMonoAndIntChar(vector<double>& monoRaw,map<int,char> & massChar,set<string> &set)
{
    map<int,string> intStringMap;//前一个为下标位置，后一个为下标对应的tag;
    //针对可能存在的多个，将string改进为vector<string>
    map<int,vector<string>> intVectorMap;//前一个为下标位置，后一个为下标对应的Vector tag;
    vector<double> mono ;
    for (int i = 0 ; i < monoRaw.size(); ++i) {
        if ( i > 1 && (int)monoRaw[i] == (int) monoRaw[i-1]){
            continue ;
        }
        mono.push_back(monoRaw[i-1]);
    }
    //进行tag的判定
    int minMass = massChar.begin()->first;
    int maxMass ;
    for(auto it = massChar.begin(); it!=massChar.end();++it) {

        if (it->first > maxMass) {
            maxMass = it->first;
        }
    }
    for(int i = 1 ; i < mono.size() ; ++i) {
        for (int j = i - 1; j >= 0 ; j--) {
            int temp = mono[i] - mono[j];
            string ts ;
            if (!(temp >= minMass && temp <= maxMass)) {
                continue;
            }
            if (massChar.count(temp)) {
                char c = massChar[temp] ;   //得到当前氨基酸
                ts += c ;                   //加上当前氨基酸
                if (intStringMap.count(j)) {        //如果被减的值，拥有了氨基酸序列,则进行拼接
                    string s = intStringMap[j];
                    ts += s ;
                    intStringMap.insert(make_pair(i,ts));
//                    cout << " mono[" <<i <<"] = " << mono[i] << " mono["<<j<<"]=  " <<mono[j]
//                         <<" temp = " << temp
//                         << " ts = " << ts
//                         <<endl;
                } else {        //直接加入
                    intStringMap.insert(make_pair(i,ts));
//                    cout << " mono[" <<i <<"] = " << mono[i] << " mono["<<j<<"]=  " <<mono[j]
//                         <<" temp = " << temp
//                         << " ts = " << ts
//                         <<endl;
                }
                break ;
            }
        }
    }

    for (auto it = intStringMap.begin(); it != intStringMap.end(); ++it) {
        set.insert(it->second);
    }
//    for(string s : set) {
//        cout <<" s = "<<s <<endl;
//    }
}
/**
 * 根据monoMass 制造tag
 * @param monoMass
 * @param cMass
 *
 * @return
 */
void unknown::ptm_utils::getTagByMonoMass(vector<double> &monoMass, map<char,double> & cMass,
                                          set<string> &set, map<char,char> &replaceMap)
{
    vector<double> mono = monoMass ;
    vector<double> pMass ;
    for(auto it = cMass.begin() ; it != cMass.end(); ++it) {
        pMass.push_back(it->second);
    }
    sort(mono.begin(),mono.end());
    sort(pMass.begin(),pMass.end());

    //map<char,char> replaceMap; //前一个为被替换的氨基酸，后一个为替换后的氨基酸
    map<int,char> massChar ; //前一个为标准化之后的数据，后一个为对应的氨基酸
    //质量相近的氨基酸进行统一替换,并得到替换的氨基酸map和一个int，char的map
    replaceCMap(cMass,replaceMap,massChar); //得到替换的氨基酸map
    //set<string> set ; //保存所有的tag
    getTagByMonoAndIntChar(mono,massChar,set);  //得到tag

}

/**
 * 得到互补map key-下标，value-下标
 * @param precursorMass 前提质量。
 * @param mono          质谱
 * @param indexMap
 */
void unknown::ptm_utils::getComplementIonsMap(double precursorMass, const vector<double> &mono, map<int,int> &indexMap)
{
    set<int> s ; //用于去重复
    for (int i = 0 ; i < mono.size(); ++i) {
        double t = precursorMass - mono[i];
        int index = mylib::data_stream::txpb_binarey_search_ex(mono, mono.size(), t);
        if ( index>=0 && !s.count(index) && fabs(t - mono[index]) < 2) {
            indexMap.insert(make_pair(i,index));
            s.insert(i);
        }
    }
}


void unknown::ptm_utils::getCutMinAndMax(vector<double> & cT, double precursorMass, double scopeValue, int &cutMin, int &cutMax)
{
    double minMass = precursorMass - scopeValue ;
    double maxMass = precursorMass + scopeValue ;
    double seqMaxMass = cT[cT.size()-1];
    if (seqMaxMass < minMass) {
        cutMin = 0 ;
        cutMax = 0 ;
        return ;
    }
    cutMin = mylib::data_stream::txpb_binarey_search_ex(cT, cT.size(), minMass);
    cutMax = mylib::data_stream::txpb_binarey_search_ex(cT, cT.size(), maxMass);

}


/**
 * 得到该蛋白质形式下ptm mass对应的unknown mass 用map保存
 * @param cT    C端理论峰
 * @param ptmMass   ptm质量序列
 * @param precursorMass 蛋白前体质量
 * @param location  当前理论峰的位置
 * @param ptmUnknown    key - ptmMass .value - unknown mass
 */
void unknown::ptm_utils::getUnknownMassByPtmMass(const vector<double> & cT, const vector<double> & ptmMass,
                                                 double precursorMass, int location, map<double,double> &ptmUnknown,
                                                 double scopeValue)
{
    ptmUnknown.clear();
    double locationMass = 0 ;
    if (location >= 0 && location < cT.size()) {
        locationMass = cT[location];        //得到肽段质量
    }
    for (double ptm : ptmMass) {    //得到每个已知修饰的质量
        if (fabs(ptm) < 0.00001) {
            continue;
        }
        double allPtmMass = precursorMass - locationMass ; //得到总的理论修饰质量
        double unknownMass = allPtmMass - ptm ;     //得到未知修饰的质量
        if (fabs(unknownMass) > scopeValue) {
            continue;
        } else {
            ptmUnknown.insert(make_pair(ptm,unknownMass));
        }
    }
}

//判断N端的
bool judgeNTIonsIsRight(const vector<double> &nT,const vector<double> &cT,
                        double nMono,int cutCount,int cutMax,int complementLoc,double nTCutMass,double pMass)
{
    //得到当前蛋白质形的互补质量 nT ->value 截断个数为cutCount . complementLoc = index ,index + 1
    // index1 + 1 + index2 + 1 = cutMax -> index2 = cutMax - index1 - 2 + cutCount;
    int index1 = complementLoc ;
    int index2 = cutMax - index1 - 2 + cutCount ;
    double complementNTIons = 0 ;
    if (index2 >= 0 && index2 < nT.size()) {
        complementNTIons = nT[index2] - nTCutMass;
    }
//    cout << " complementNTIons  = " << complementNTIons <<endl;
//    cout << " pMass = " << pMass <<endl;
//    cout << " complementNTIons + pMass = " << complementNTIons + pMass<<endl;
//    cout << " nTCutMass = " << nTCutMass << endl;
//    cout << " nMono = " << nMono <<endl;
//    cout << endl;
    if (fabs(nMono - complementNTIons - pMass) < 1.2) {
        return true;
    }
    return false ;
}

/**
 *
 * @param m1    互补离子1
 * @param m2    互补离子2
 * @param m1Index   互补离子1的下标
 * @param m2Index   互补离子2的下标
 * @param cT    c端理论峰
 * @param ptmUnknown    key为已知ptm修饰，value为未知修饰
 * @param scopeValue
 * @return
 */
bool unknown::ptm_utils::getComplementIonsCTIndex(double m1, double m2, int &m1Index, int &m2Index, int cutMax,
                                                  const vector<double> & cT, const vector<double> & nT,
                                                  vector<int> &trueModIndex,
                                                  const map<double,double> & ptmUnknown, double nTCutMass,
                                                  double scopeValue)
{
    int count = 0 ;
    //如果带修饰,判断是否是当前的蛋白质形带
//    cout << " m1 = " << m1 << " m2 = " << m2 << " cutSize = " << nT.size() - cutMax << " cutMax = " << cutMax <<endl;
    int index = 0 ;
    for (auto it = ptmUnknown.begin(); it != ptmUnknown.end(); ++it,index++)        //遍历未知的修饰
    {
        int m1Loc = -1;
        int m2Loc = -1;
        int m1Log = -1 ; // 0 表示不带 修饰 1表示带ptm修饰 2表示带未知修饰  3表示带混合修饰
        int m2Log = -1 ;
        double ptmMass = it->first;
        double unknownMass = it->second ;
        double allMass = ptmMass + unknownMass ;
        m1Loc = mylib::data_stream::txpb_binarey_search_ex(cT, cT.size(), m1 - ptmMass);
        if (m1Loc >= 0 && fabs(m1-ptmMass - cT[m1Loc]) < scopeValue) {
            m1Index = m1Loc ;
            m1Log = 1;
        }
        m1Loc = mylib::data_stream::txpb_binarey_search_ex(cT, cT.size(), m1 - unknownMass);
        if (m1Loc >= 0 && fabs(m1-unknownMass - cT[m1Loc]) < scopeValue) {
            m1Index = m1Loc ;
            m1Loc = 2;
        }
//        m1Loc = g::proc::txpb_binarey_search_ex(cT,cT.size(),m1-allMass);
//        if (m1Loc >= 0 && fabs(m1-allMass - cT[m1Loc]) < scopeValue) {
//            m1Index = m1Loc ;
//            m1Loc = 3 ;
//        }
        //判定另外一个离子,正确性
        if (m1Index != -1) {    //说明m1 是C端离子
//            cout << " m1Index = " << m1Index  << " m1Log = " << m1Log
//            <<" ontPtmMass = " << ptmMass
//            << "unknownMass = " << unknownMass
//            << "all mass = " <<allMass
//            << "findIonsMass = " << cT[m1Index]
//            << endl;
//            cout <<endl;
            if (m1Log == 1) {       //如果是带ptm修饰
                bool t = judgeNTIonsIsRight(nT,cT,m2,nT.size()-cutMax,cutMax,m1Index,nTCutMass,unknownMass);
                if (t) {
                    trueModIndex.push_back(index);
                    count ++ ;
                    continue;
                }
            } else if (m1Log == 2) {    //如果是带未知修饰
                bool t = judgeNTIonsIsRight(nT,cT,m2,nT.size()-cutMax,cutMax,m1Index,nTCutMass,ptmMass);
                if (t) {
                    trueModIndex.push_back(index);
                    count ++ ;
                    continue ;
                }
            } else if (m1Log == 3) {   //如果是带混合修饰
                bool t = judgeNTIonsIsRight(nT,cT,m2,nT.size()-cutMax,cutMax,m1Index,nTCutMass,0);
                if (t) {
                    trueModIndex.push_back(index);
                    count ++ ;
                    continue ;
                }
            }
        }

        m2Loc = mylib::data_stream::txpb_binarey_search_ex(cT, cT.size(), m2 - ptmMass);
        if (m2Loc >= 0 && fabs(m2 - ptmMass - cT[m2Loc]) < scopeValue) {
            m2Index = m2Loc ;
            m2Log = 1;
        }
        m2Loc = mylib::data_stream::txpb_binarey_search_ex(cT, cT.size(), m2 - unknownMass);
        if (m2Loc >= 0 && fabs(m2 - unknownMass - cT[m2Loc]) < scopeValue) {
            m2Index = m2Loc ;
            m2Log = 2;
        }
//        m2Loc = g::proc::txpb_binarey_search_ex(cT,cT.size(),m2-allMass);
//        if (m2Loc >= 0 && fabs(m2 - allMass - cT[m2Loc]) < scopeValue) {
//            m2Index = m2Loc ;
//            m2Log = 3;
//        }

        //判定另外一个离子,正确性
        if (m2Index != -1) {    //说明m2 是C端离子
//            cout << " m2Index = " << m2Index  << " m2Log = " << m2Log
//            <<" ontPtmMass = " << ptmMass
//            << " unknownMass = " << unknownMass
//            << " all mass = " <<allMass
//            << " findIonsMass = " << cT[m2Index]
//            << endl;
//            cout <<endl;
            if (m2Log == 1) {       //如果是带ptm修饰
                bool t = judgeNTIonsIsRight(nT,cT,m1,nT.size()-cutMax,cutMax,m2Index,nTCutMass,unknownMass);
                if (t) {
                    trueModIndex.push_back(index);
                    count ++ ;
                    continue ;
                }
            } else if (m2Log == 2) {    //如果是带未知修饰
                bool t = judgeNTIonsIsRight(nT,cT,m1,nT.size()-cutMax,cutMax,m2Index,nTCutMass,ptmMass);
                if (t) {
                    trueModIndex.push_back(index);
                    count ++ ;
                    continue ;
                }
            } else if (m2Log == 3) {   //如果是带混合修饰
                bool t = judgeNTIonsIsRight(nT,cT,m1,nT.size()-cutMax,cutMax,m2Index,nTCutMass,0);
                if (t) {
                    trueModIndex.push_back(index);
                    count ++ ;
                    continue ;
                }
            }
        }
    }

    return trueModIndex.size() > 0 ? true : false;
}


//void unknown::ptm_utils::getConPro(vector<prsm::protein_processor_ptr> & p1, vector<prsm::protein_processor_ptr> & p2, vector<prsm::protein_processor_ptr> & con)
//{
//    map<prsm::protein_processor_ptr,int> map1 ;
//    for (int i = 0; i < p1.size(); ++i) {
//        map1.insert(make_pair(p1[i],1));
//    }
//    for (int i = 0; i < p2.size(); ++i) {
//        if (map1.count(p2[i])){
//            map1[p2[i]]++ ;
//        }
//    }
//    for (auto it = map1.begin(); it != map1.end(); ++it) {
//        if (it->second == 2) {
//            con.push_back(it->first);
//        }
//    }
//}