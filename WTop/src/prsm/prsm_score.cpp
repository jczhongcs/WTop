//
// Created by wenzhong on 2023/3/22.
//

#include "prsm_score.hpp"


double prsm_score::analysis_mutx(const std::vector<node> &arr_ions, const std::vector<node> &arr_ions_2,
                                             const std::vector<double> &mono,
                                             const std::vector<double> &arr_mod, double mod_mass) {
    std::vector<node> all = arr_ions;
    int count_zero = 0;
    int count_error = 0;
    int count_match = 0;
    for (int i = 0; i < arr_ions_2.size(); ++i) {
        all.push_back(arr_ions_2[i]);
    }
    for (int i = 0; i < all.size(); ++i) {
        if (fabs(all[i].index) < 1.1) {     //首端匹配的离子
            count_zero++;
        }
        if (fabs(all[i].index - mod_mass) < 2) {  //末端匹配的离子
            count_match++;
        }
        if (!mylib::data_stream::if_is_error_ions(arr_mod, all[i])) {    //匹配不对齐的离子
            count_error++;
        }
    }
    double p = 0;
    double all_size = all.size();
    double mono_size = mono.size();
    if (all_size / mono_size * 100 > 10 && all_size / mono_size * 100 <= 15) {
        p = 1;
    } else if (all_size / mono_size * 100 > 15 && all_size / mono_size * 100 <= 30) {
        p = 2;
    } else if (all_size / mono_size * 100 > 30) {
        p = 3;
    }
    //如果满足count_zero > 0 , count_match > 0 ; 最终得分翻倍。
//    if (fabs(mod_mass)<=eps && count_zero > 0 && count_match > 0) {
//        return 2 + p ;
//    } else
    if (count_zero > 0 && count_match > 0) {
        return 2 + p;
    } else if (count_zero > 0 || count_match > 0) {
        return 1 + p;
    }
    return 0.5 + p;
}

double prsm_score::get_complementIons_score(vector<node> &ionsC, vector<node> &ionsN, int length) {
    double score = 0;
    set<long> setC;
    for(int i = 0; i <ionsC.size(); ++i) {
        setC.insert(ionsC[i].thoe_id+1);
    }
    int count = 0 ;
    for (int i = 0; i < ionsN.size(); ++i) {
        if (setC.count(length - (ionsN[i].thoe_id + 1) )) {
            count ++ ;
        }
    }
    return count ;
}

double prsm_score::get_seqIons_score(vector<node> &ionsC, vector<node> &ionsN, int &count) {
    double sumSeqIons = 1 ;
    double scoreC = 0 ;
    double scoreN = 0 ;
    for (int i = 1 ; i < ionsC.size(); ++i) {
        if (ionsC[i].thoe_id - ionsC[i-1].thoe_id <= 1) {
            sumSeqIons++;
        }
    }
    for (int i = 1 ; i < ionsN.size(); ++i) {
        if (ionsN[i].thoe_id - ionsN[i-1].thoe_id <= 1) {
            sumSeqIons++;
        }
    }
    count = sumSeqIons ;
    double cof = sumSeqIons / (double)(ionsC.size()+ionsN.size()) ;
    return sqrt(cof * 100) ;
}

int prsm_score::get_SeqIons_pair(vector<node> & t, vector<node> & m)
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

void prsm_score::get_score(vector<node> &ions_r, vector<double> &refe_score,
                           const vector<double> &mod_arr, vector<double> &peer_score,
                           double &score, double base_score, double mod_mass) {
    vector<node> ions;
    double count = 0;
    double punish_nub = 0;
    //如果修饰质量太大，可能误差就大，故适当降低
    if (fabs(mod_mass) > 500.0 / 2.0) {
        for (int pi = 0; pi < ions_r.size(); ++pi) {
            if (!mylib::data_stream::if_is_error_ions(mod_arr, ions_r[pi])) {
                count++;
                continue;
            }
            ions.push_back(ions_r[pi]);
        }
        punish_nub = 1 - (count / ions.size());
    } else {
        ions = ions_r;
        punish_nub = 1;
    }
    //如果没有匹配离子，则得分为 0
    if (ions.size() == 0) {
        score = 0;
        return;
    }
    int i = 0;
    double p = 0.0;    //标识修饰离子个数
    double z = 0.0;    //标识完全对准个数
    vector<double> mutScore(ions.size(), 0);  //扩大系数
    vector<int> matchIndex(ions.size(), 0);  //标识是否为完全匹配
    vector<int> matchSeq(ions.size(), 0);   //标识是否连续
    map<int, int> rm_dup;
    //标识出完全匹配的离子和修饰离子matchIndex
    //得到完全匹配离子的扩大系数和修饰离子的扩大系数
    for (i = 0; i < ions.size(); i++) {
        if (!rm_dup.count(ions[i].thoe_id)) {
            rm_dup[ions[i].thoe_id] = 1;
        }
        if (fabs(ions[i].index) < 1.1) {
            matchIndex[i] = 1;
            ++z;
        } else {
            ++p;
        }
    }
    //基础得分系数
    p = 1 - (p / (ions.size())) + base_score;        //缩小修饰得分比例
    z = 1 + (z / (ions.size())) + base_score;        //扩大完全对准得分比例
    //标识出连续匹配的离子matchSeq
    for (i = 0; i < ions.size() - 1; i++) {
        if (fabs(ions[i].thoe_id - ions[i + 1].thoe_id) == 1) {
            matchSeq[i] = 1;
            matchSeq[i + 1] = 1;
        }
    }
    //    p = 1 - (p / (ions.size())) + 0.1;        //缩小修饰得分比例
    //    z = 1 + (z / (ions.size())) + 0.1;        //扩大完全对准得分比例
    //当只有完全匹配的离子时，判定是否存在连续匹配离子
    if ((1 - (p / (ions.size()))) == 0) {
        int lo = 0;
        for (i = 0; i < matchSeq.size(); i++) {
            if (matchSeq[i] == 1) {
                lo = 1;
                break;
            }
        }
        if (!lo) {
            z -= (z / (ions.size()));
        }
    }
    //设置连续匹配系数mutScore
    for (i = 0; i < ions.size(); i++) {
        if (matchSeq[i] == 0) {
            if (matchIndex[i]) {
                mutScore[i] += z;
            } else {
                mutScore[i] += p;
            }
        } else if (matchSeq[i] == 1) {
            double log = 0;
            int first = i;
            while (i < ions.size() && matchSeq[i] == 1 ) {  //记下多少个连续匹配的
                log++;
                i++;
            }
            double logs = log / 10.0;
            //            double logs = 0 ;
            while (first != i) {
                if (matchIndex[first]) {
                    mutScore[first] += (z + logs);
                } else {
                    mutScore[first] += (p + logs);
                }
                first++;
            }
            if (i != ions.size()) {
                i--;
            }
        }
    }

//            cout << "完全匹配峰" << endl;
//            for (i = 0; i < matchIndex.size(); ++i) {
//                cout << " " << matchIndex[i];
//            }
//            cout << endl;
//            cout << "匹配连续峰" << endl;
//            for (i = 0; i < matchSeq.size(); ++i) {
//                cout << " " << matchSeq[i];
//            }
//            cout << endl;
//            cout << "系数" << endl;
//            for (i = 0; i < mutScore.size(); ++i) {
//                cout << " " << mutScore[i];
//            }
//            cout << endl;

    //22-10-12添加模块，如果存在修饰，那么该修饰有几个峰值10个峰值为1 ， 1个峰值为0.1 。 2个0.2
    double base = 0.1 ;
    vector<double> countBase ; //计算有几个相同的
    for (i = 0; i < ions.size(); ++i) {
        double hasSeqBase = 0.1 ;
        int j = i + 1;
        set<int> dup ;
        dup.insert(ions[i].thoe_id);
        for (; j < ions.size(); ++j) {
            if ( dup.count(ions[j].thoe_id)) {
                continue ;
            }
            if (fabs(ions[j].index - ions[i].index) > 1.5) {
                break ;
            }
            hasSeqBase += base ;
            dup.insert(ions[j].thoe_id);
        }
        while (i != j) {
            countBase.push_back(hasSeqBase);
            i++;
        }
        i = j - 1;
    }
//    for (double c : countBase) {
//        cout << " c = " << c << endl;
//    }
//
//    cout << " ==== " << endl;
    for (i = 0; i < ions.size(); ++i) {
        if (rm_dup[ions[i].thoe_id]) { //refer 和peer的分都算
            rm_dup[ions[i].thoe_id] = 0;
            score += (fabs(refe_score[ions[i].thoe_id]) + fabs(peer_score[ions[i].mono_id]) * mutScore[i]) * countBase[i];
        } else {    //只算peer的分
            score += fabs(peer_score[ions[i].mono_id]) * mutScore[i] * countBase[i];
        }
    }
    score *= punish_nub;
    //    cout<<"score = "<<score<<endl;
}

int prsm_score::get_seq_count(vector<node> &ions_r) {
    int seq_match_r = 0;
    if(ions_r.size() <=1) {
        return 1;
    }

    for (int i = 0; i<(ions_r.size()-2); i++) {

        //cout<<endl<<"ions_r[i].thoe_id="<<ions_r[i].thoe_id<<"   ions_r[i+1].thoe_id="<<ions_r[i+1].thoe_id<<endl;
        if (fabs(ions_r[i].thoe_id - ions_r[i + 1].thoe_id) == 1) {
            seq_match_r++;
        }
    }
    return seq_match_r;
}

double prsm_score::get_abundance_score(vector<node> &ions_n,const vector<double> &IonsMass,const map<double,double> &keyMsValueAbRadio){
    double score=0;
    for(node ion : ions_n){
        vector<double> abunds;
        double ion_abundence = keyMsValueAbRadio.at(IonsMass[ion.mono_id]);

//        cout<<"ion_mass:"<<IonsMass[ion.mono_id]<<"   abundence:"<<ion_abundence<<endl;
        abunds.push_back(ion_abundence);
        for(int i =0;i<IonsMass.size();i++){
            if(fabs(IonsMass[i]-IonsMass[ion.mono_id])<500)
                abunds.push_back(keyMsValueAbRadio.at(IonsMass[i]));
        }
        std::sort(abunds.begin(), abunds.end(), std::greater<double>());
        // 使用 std::find() 查找给定数字在排序后的向量中的位置
        auto it = std::find(abunds.begin(), abunds.end(), ion_abundence);

        int rank = -1;
        // 如果找到了目标数字，则计算排名并返回
        if (it != abunds.end()) {
            rank = std::distance(abunds.begin(), it) + 1;
        } else{
            continue;
        }

        double abun_score = max(log(20/rank),0.0);
        score += abun_score;
    }
    return score;
}

double prsm_score::get_dyn_abu_score(vector<node> &ions_n,const vector<double> &IonsMass,const map<double,double> &keyMsValueAbRadio,double precursor_mass){
    double score=0;
    double jug_mass = precursor_mass/(double)70;
    for(node ion : ions_n){
        vector<double> abunds;
        double ion_abundence = keyMsValueAbRadio.at(IonsMass[ion.mono_id]);

//        cout<<"ion_mass:"<<IonsMass[ion.mono_id]<<"   abundence:"<<ion_abundence<<endl;
        abunds.push_back(ion_abundence);
        for(int i =0;i<IonsMass.size();i++){
            if(fabs(IonsMass[i]-IonsMass[ion.mono_id]) < jug_mass)
                abunds.push_back(keyMsValueAbRadio.at(IonsMass[i]));
        }
        std::sort(abunds.begin(), abunds.end(), std::greater<double>());
        // 使用 std::find() 查找给定数字在排序后的向量中的位置
        auto it = std::find(abunds.begin(), abunds.end(), ion_abundence);

        int rank = -1;
        // 如果找到了目标数字，则计算排名并返回
        if (it != abunds.end()) {
            rank = std::distance(abunds.begin(), it) + 1;
        } else{
            continue;
        }

        double abun_score = max(log(20/rank),0.0);
        score += abun_score;
    }
    return score;
}

double prsm_score::get_avg_abundance_score(vector<node> &ions_n,vector<node> &ions_c,const vector<double> &IonsMass,const map<double,double> &keyMsValueAbRadio){
    double score = 0,count = 0;
    vector<node> all_ions;

    for(node in:ions_n){
        all_ions.push_back(in);
    }
    for(node ic:ions_c){
        all_ions.push_back(ic);
    }
    for(node ion : all_ions){
        vector<double> abunds;
        double ion_abundence = keyMsValueAbRadio.at(IonsMass[ion.mono_id]);

//        cout<<"ion_mass:"<<IonsMass[ion.mono_id]<<"   abundence:"<<ion_abundence<<endl;
        abunds.push_back(ion_abundence);
        for(int i =0;i<IonsMass.size();i++){
            if(fabs(IonsMass[i]-IonsMass[ion.mono_id])<500)
                abunds.push_back(keyMsValueAbRadio.at(IonsMass[i]));
        }
        std::sort(abunds.begin(), abunds.end(), std::greater<double>());
        // 使用 std::find() 查找给定数字在排序后的向量中的位置
        auto it = std::find(abunds.begin(), abunds.end(), ion_abundence);

        int rank = -1;
        // 如果找到了目标数字，则计算排名并返回
        if (it != abunds.end()) {
            rank = std::distance(abunds.begin(), it) + 1;
        } else{
            continue;
        }

        double abun_score = max(log(20/rank),0.0);
        if(abun_score > 0 ) count++;
        score += abun_score;
    }
    if(count == 0)  return 0;

    score /= (double)count;
    return score;
}

double prsm_score::get_avg_dyn_abu_score(vector<node> &ions_n,vector<node> &ions_c,const vector<double> &IonsMass,const map<double,double> &keyMsValueAbRadio,double precursor_mass){
    double score = 0,count = 0;

    double jug_mass = precursor_mass/(double)70;
    vector<node> all_ions;

    for(node in:ions_n){
        all_ions.push_back(in);
    }
    for(node ic:ions_c){
        all_ions.push_back(ic);
    }

    for(node ion : all_ions){
        vector<double> abunds;
        double ion_abundence = keyMsValueAbRadio.at(IonsMass[ion.mono_id]);

//        cout<<"ion_mass:"<<IonsMass[ion.mono_id]<<"   abundence:"<<ion_abundence<<endl;
        abunds.push_back(ion_abundence);
        for(int i =0;i<IonsMass.size();i++){
            if(fabs(IonsMass[i]-IonsMass[ion.mono_id]) < jug_mass)
                abunds.push_back(keyMsValueAbRadio.at(IonsMass[i]));
        }

        std::sort(abunds.begin(), abunds.end(), std::greater<double>());
        // 使用 std::find() 查找给定数字在排序后的向量中的位置
        auto it = std::find(abunds.begin(), abunds.end(), ion_abundence);

        int rank = -1;
        // 如果找到了目标数字，则计算排名并返回
        if (it != abunds.end()) {
            rank = std::distance(abunds.begin(), it) + 1;
        } else{
            continue;
        }

        double abun_score = max(log(20/rank),0.0);
        if(abun_score > 0) count++;
        score += abun_score;
    }
    if(count == 0)  return 0;


    score /= (double)count;
    return score;
}


int prsm_score::get_diffcharge_score(vector<node> &ions_n){
    int diffcharge_score = 0;
    map<int,int> keyThoeidValueCount;
    for(node ion : ions_n){
        if(keyThoeidValueCount.find(ion.thoe_id) == keyThoeidValueCount.end())
            keyThoeidValueCount.insert(std::make_pair(ion.thoe_id, 0));
        keyThoeidValueCount[ion.thoe_id]++;
    }

    for(auto m:keyThoeidValueCount) {
        if (m.second > 1)
            diffcharge_score += m.second;
    }

    return diffcharge_score;
}