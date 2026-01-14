//
// Created by wenzhong on 2023/3/22.
//

#include "prsm_ions_location.hpp"


void n_terminal_transaction_c_terminal(vector<modi> &N_Terminal_Inter, vector<modi> &C_Terminal_Inter, int log) {
    log -= 2 ;
    for (int xi = N_Terminal_Inter.size() - 1; xi >= 0; xi--)
    {
        modi t;
        t.first = log - N_Terminal_Inter[xi].second ;
        t.second = log - N_Terminal_Inter[xi].first ;
        t.mod_mass = N_Terminal_Inter[xi].mod_mass;
        C_Terminal_Inter.push_back(t);
    }
}

void get_Cterminal_ions_and_change_ptms_tnter(vector<modi> &Terminal_Inter, vector<node> &merg_ions,std::map<long, long> &mapDup,
                                               vector<node> &get_ions,vector<modi> &ptmsTerminal, double ptm_mass)
{
    if (Terminal_Inter.empty())
    {
        return;
    }
    map<long, long> mapD;       //用于去重
    bool firstInt = false;      // 是否是第一个区间
    //遍历每一个可能的离子,取可能的离子
    for (int xi = 0; xi < merg_ions.size(); xi++)
    {
        long first = merg_ions[xi].thoe_id;
        long second = merg_ions[xi].mono_id;
        if (mapDup.count(second)) //如果该质谱峰已有匹配项,过滤
        {
            continue;
        }
        double maxMass = 0;    //记录区间的修饰质量
        int pLoc = 0;        //pLoc [0,Terminal_Inter.size() - 1 ],标识在第几个修饰区间
        //得出该理论峰的区间位置以及修饰质量限值
        for (int ji = 0; ji < Terminal_Inter.size(); ji++)
        {
            if (first <= Terminal_Inter[ji].first)  //实现选取峰值时，前开后闭的选取原则
            {
                break;
            }
            maxMass += Terminal_Inter[ji].mod_mass; //如果在下一个区间，则加上下一个区间的修饰质量
            pLoc = ji;
        }
//        int curInterFirst = Terminal_Inter[pLoc].first ;    //当前区间的起始位置
        if (pLoc == 0) //如果在第一个区间内
        {
            if (fabs(merg_ions[xi].index) < 2 && !mapD.count(merg_ions[xi].thoe_id) && !firstInt) //可以取 0
            {
                get_ions.push_back(merg_ions[xi]);        //将该离子选入
                mapD[merg_ions[xi].thoe_id] = merg_ions[xi].mono_id;    //将理论峰与质谱峰添加进map 去重
                //更新修饰区间
                if (merg_ions[xi].thoe_id > ptmsTerminal[pLoc].first && merg_ions[xi].thoe_id < ptmsTerminal[pLoc].second)
                {
                    ptmsTerminal[pLoc].first = merg_ions[xi].thoe_id;
                }
//                if (merg_ions[xi].thoe_id > Terminal_Inter[pLoc].first &&  merg_ions[xi].thoe_id < Terminal_Inter[pLoc].second)
//                {
//                    Terminal_Inter[pLoc].first = merg_ions[xi].thoe_id;
//                }
            } else if (fabs(merg_ions[xi].index - maxMass) < 3) {       //如果是修饰离子
                get_ions.push_back(merg_ions[xi]);      //取该离子
                mapD[merg_ions[xi].thoe_id] = merg_ions[xi].mono_id;    //去重
                firstInt = true;        //不能够再取 0 离子
                //更新修饰区间
                if (merg_ions[xi].thoe_id > ptmsTerminal[pLoc].first && merg_ions[xi].thoe_id < ptmsTerminal[pLoc].second)
                {
                    ptmsTerminal[pLoc].first = merg_ions[xi].thoe_id;
                }
//                if (merg_ions[xi].thoe_id > Terminal_Inter[pLoc].first && merg_ions[xi].thoe_id < Terminal_Inter[pLoc].second)
//                {
//                    Terminal_Inter[pLoc].first = merg_ions[xi].thoe_id;
//                }
            }
        } else if (fabs(merg_ions[xi].index - maxMass) < 3) {    //如果在后续区间，且为修饰离子
            get_ions.push_back(merg_ions[xi]);
            mapD[merg_ions[xi].thoe_id] = merg_ions[xi].mono_id;
            //则更新修饰区间
            if (merg_ions[xi].thoe_id > ptmsTerminal[pLoc].first && merg_ions[xi].thoe_id < ptmsTerminal[pLoc].second)
            {
                ptmsTerminal[pLoc].first = merg_ions[xi].thoe_id;
            }
//            if (merg_ions[xi].thoe_id > Terminal_Inter[pLoc].first && merg_ions[xi].thoe_id < Terminal_Inter[pLoc].second)
//            {
//                Terminal_Inter[pLoc].first = merg_ions[xi].thoe_id;
//            }
        }
    }//For merg_ions END
}

void get_Other_terminal_ions(vector<modi> &Terminal_Inter, vector<node> &merg_ions,
                             std::map<long, long> &mapDup,
                             vector<node> &get_ions,vector<modi> &ptmsTerminal, double ptm_mass)
{
    if (Terminal_Inter.empty())
    {
        return;
    }
    map<long, long> dupMap;       //用于去重
    bool firstInt = false;      // 是否是第一个区间
    //遍历每一个可能的离子,取可能的离子
    for (int xi = 0; xi < merg_ions.size(); xi++)
    {
        long first = merg_ions[xi].thoe_id;
        long second = merg_ions[xi].mono_id;
        if (mapDup.count(second)) //如果该质谱峰已有匹配项,过滤
        {
            continue;
        }
        double maxMass = 0;    //记录区间的修饰质量
        int pLoc = -1 ;        //pLoc [0,Terminal_Inter.size() - 1 ],标识在第几个修饰区间
        //得出该理论峰的区间位置以及修饰质量限值
        for (int ji = 0; ji < Terminal_Inter.size(); ji++) {
            if (first > Terminal_Inter[ji].first && first <= Terminal_Inter[ji].second) {
                pLoc = ji ;
                break ;
            }
        }
        for (int i = 0; i < Terminal_Inter.size(); i++) {
            if (first <= Terminal_Inter[i].first) {
                break ;
            }
            maxMass += Terminal_Inter[i].mod_mass ;
        }
//        int curInterFirst = Terminal_Inter[pLoc].first ;    //当前区间的起始位置
        if (pLoc == 0) //如果在第一个区间内
        {
            if (fabs(merg_ions[xi].index) < 2 && !dupMap.count(merg_ions[xi].thoe_id) && !firstInt) //可以取 0
            {
                get_ions.push_back(merg_ions[xi]);        //将该离子选入
                dupMap[merg_ions[xi].thoe_id] = merg_ions[xi].mono_id;    //将理论峰与质谱峰添加进map 去重
                //更新修饰区间
                if (merg_ions[xi].thoe_id > ptmsTerminal[pLoc].first && merg_ions[xi].thoe_id < ptmsTerminal[pLoc].second)
                {
                    ptmsTerminal[pLoc].first = merg_ions[xi].thoe_id;
                }
//                if (merg_ions[xi].thoe_id > Terminal_Inter[pLoc].first &&  merg_ions[xi].thoe_id < Terminal_Inter[pLoc].second)
//                {
//                    Terminal_Inter[pLoc].first = merg_ions[xi].thoe_id;
//                }
            } else if (fabs(merg_ions[xi].index - maxMass) < 3) {       //如果是修饰离子
                get_ions.push_back(merg_ions[xi]);      //取该离子
                dupMap[merg_ions[xi].thoe_id] = merg_ions[xi].mono_id;    //去重
                firstInt = true;        //不能够再取 0 离子
                //更新修饰区间
                if (merg_ions[xi].thoe_id > ptmsTerminal[pLoc].first && merg_ions[xi].thoe_id < ptmsTerminal[pLoc].second)
                {
                    ptmsTerminal[pLoc].first = merg_ions[xi].thoe_id;
                }
//                if (merg_ions[xi].thoe_id > Terminal_Inter[pLoc].first && merg_ions[xi].thoe_id < Terminal_Inter[pLoc].second)
//                {
//                    Terminal_Inter[pLoc].first = merg_ions[xi].thoe_id;
//                }
            }
        } else if (fabs(merg_ions[xi].index - maxMass) < 3) {    //如果在后续区间，且为修饰离子
            get_ions.push_back(merg_ions[xi]);
            dupMap[merg_ions[xi].thoe_id] = merg_ions[xi].mono_id;
            //则更新修饰区间
            if(pLoc == -1)
                pLoc++;
            if (merg_ions[xi].thoe_id > ptmsTerminal[pLoc].first && merg_ions[xi].thoe_id < ptmsTerminal[pLoc].second)
            {
                ptmsTerminal[pLoc].first = merg_ions[xi].thoe_id;
            }
//            if (merg_ions[xi].thoe_id > Terminal_Inter[pLoc].first && merg_ions[xi].thoe_id < Terminal_Inter[pLoc].second)
//            {
//                Terminal_Inter[pLoc].first = merg_ions[xi].thoe_id;
//            }
        }
    }//For merg_ions END
}


void prsm_ions_location::ions_location_match_N(vector<node> &n_ions, vector<node> &c_ions,
                                               vector<node> &merg_ions_n,
                                               vector<node> &merg_ions_c, vector<modi> &mod,
                                               const vector<double> &modifyTable,
                                               int lenth, double modifier_mass)
{
    vector<modi> N_Terminal_Inter;    //计算N端修饰区间
    vector<modi> C_Terminal_Inter;    //计算C端修饰区间
    std::map<long, long> mapDup;      //first保存质谱峰下标，second保存理论峰下标，用于去重
    int log = lenth - 1;
    //最长路径算法
    mylib::ions_path_Ptr long_ions_path_ptr = new mylib::ions_path();
    long_ions_path_ptr->initializeData(merg_ions_n, merg_ions_c, modifyTable, n_ions, lenth);
    delete long_ions_path_ptr;
    if (n_ions.size() == 0)
    {
        return ;
    }

    /// 计算修饰区间 （node -> modi 的转变)
    mylib::etd_process::index_Unknown_Ptm(n_ions, N_Terminal_Inter, modifyTable);
    /// 2022.9.2 添加 -- 进行修饰质量的判断，假设缺少的修饰质量在后半段
    mylib::etd_process::assumMissingMass(n_ions, N_Terminal_Inter, modifyTable, modifier_mass, lenth);

    ///  转化为另外一端可取离子区间
    n_terminal_transaction_c_terminal(N_Terminal_Inter, C_Terminal_Inter, lenth);

    ///  根据另外一端修饰区间取离子，并且改变修饰区间，区间的取值是前开后闭
    vector<modi> ptmsTerminal = C_Terminal_Inter ;  //变动后的区间
    get_Other_terminal_ions(C_Terminal_Inter, merg_ions_c, mapDup, c_ions, ptmsTerminal,modifier_mass);

    /// 由最终的 terminal-Inter 得出最终修饰区间
    log = lenth - 2 ;
    for (int xi = ptmsTerminal.size() - 1; xi >= 0; xi--) {
        modi t;
        t.first = log - ptmsTerminal[xi].second;
        t.second = log - ptmsTerminal[xi].first;
        t.mod_mass = ptmsTerminal[xi].mod_mass;
        mod.push_back(t);
    }

    return;
}

void prsm_ions_location::ions_location_match_C(vector<node> &n_ions, vector<node> &c_ions,
                                               vector<node> &merg_ions_n, vector<node> &merg_ions_c,
                                               vector<modi> &mod, const vector<double> &modifyTable, int lenth, double modifier_mass)
{
    int i = 0, j = 0;
    vector<modi> N_Terminal_Inter;    //计算N端修饰区间
    vector<modi> C_Terminal_Inter;    //计算C端修饰区间
    std::map<long, long> mapDup;      //first保存质谱峰下标，second保存理论峰下标，用于去重
    int log = lenth - 1;
    /**
     * 求出最长路径
     */
    mylib::ions_path_Ptr lpCP = new mylib::ions_path();
    lpCP->initializeData(merg_ions_c, merg_ions_n, modifyTable, c_ions, lenth);
    delete lpCP;
    //主要用于利用n端离子制造修饰区间( node -> modi)
    if (c_ions.size() == 0) {
        return ;
    }
    /**
     * 得出相应的修饰区间
     */
    mylib::etd_process::index_Unknown_Ptm(c_ions, C_Terminal_Inter, modifyTable);
//    for (node t : c_ions) {
//        cout
//        <<" f : "<< t.thoe_id
//        << " m : " << t.mono_id
//        << " mass : " << t.index <<endl;
//    }
    /**
     * 2022.9.2 添加 -- 进行修饰质量的判断，假设缺少的修饰质量在后半段
     */
    mylib::etd_process::assumMissingMass(c_ions, C_Terminal_Inter, modifyTable, modifier_mass, lenth);
//    cout << " C mod inter : \n" ;
//    for (i = 0; i < C_Terminal_Inter.size(); i++) {
//        cout << " C_varible : first = " << C_Terminal_Inter[i].first << " second = " << C_Terminal_Inter[i].second
//             << " mass = " << C_Terminal_Inter[i].mod_mass << endl;
//    }

    //由该端修饰区间转化为另一端修饰区间
    n_terminal_transaction_c_terminal(C_Terminal_Inter, N_Terminal_Inter, lenth);
//        cout <<" =========== "<<endl;
//     for (i = 0 ; i < N_Terminal_Inter.size() ; i++) {
//     	cout
//     	<<" N_varible : first = "<<N_Terminal_Inter[i].first
//     	<<" second = "<<N_Terminal_Inter[i].second
//     	<<" mass = " <<N_Terminal_Inter[i].mod_mass
//     	<<endl;
//     }
    //根据修饰区间，另外一端所有离子，取出另外一端的离子
    vector<modi> ptmsTerminal  = N_Terminal_Inter;
    get_Other_terminal_ions(N_Terminal_Inter, merg_ions_n, mapDup, n_ions, ptmsTerminal,modifier_mass);

//    get_Cterminal_ions_and_change_ptms_tnter(N_Terminal_Inter, merg_ions_n, mapDup, n_ions, ptmsTerminal,modifier_mass);

//    for (i = 0 ; i < ptmsTerminal.size() ; i++) {
//        cout
//                <<" ptmsTerminal : first = "<<ptmsTerminal[i].first
//                <<" second = "<<ptmsTerminal[i].second
//                <<" mass = " <<ptmsTerminal[i].mod_mass
//                <<endl;
//    }
//    cout << modifier_mass <<endl;
//         cout<<"N ions : \n";
//     for(int xi = 0 ; xi < n_ions.size() ; xi++ ) {
//         cout<<"thero id = "<<n_ions[xi].thoe_id<<", mono id = "<<n_ions[xi].mono_id<<" index = "<<n_ions[xi].index<<endl;
//     }
//    cout<<"C ions : \n";
//    for(int xi = 0 ; xi < c_ions.size() ; xi++ ) {
//        cout<<"thero id = "<<c_ions[xi].thoe_id<<", mono id = "<<c_ions[xi].mono_id<<" index = "<<c_ions[xi].index<<endl;
//    }
    //-------------------------------由最终的N-terminal-Inter 得出最终修饰区间
    mod = ptmsTerminal;
    return;
}


void prsm_ions_location::confirm_final_ions(vector<node> &n_ions_n, vector<node> &n_ions_c,vector<modi> & n_ptms,
                                            vector<node> &c_ions_n, vector<node> &c_ions_c,vector<modi> & c_ptms,
                                            prsm::temp_prsm_argument_ptr tempPrsmArgumentPtr)
{
    ///原函数
//    if (tempPrsmArgumentPtr->var_mod_mass != 0) {
//        if (n_ions_n.size() + n_ions_c.size() < c_ions_n.size() + c_ions_c.size()) {
//            tempPrsmArgumentPtr->ions_n = c_ions_n;
//            tempPrsmArgumentPtr->ions_c = c_ions_c;
//            tempPrsmArgumentPtr->ptms = c_ptms;
//        } else {
//            tempPrsmArgumentPtr->ions_n = n_ions_n;
//            tempPrsmArgumentPtr->ions_c = n_ions_c;
//            tempPrsmArgumentPtr->ptms = n_ptms;
//        }
//    } else if (tempPrsmArgumentPtr->var_mod_mass == 0) {
//        tempPrsmArgumentPtr->ions_n = n_ions_n;
//        tempPrsmArgumentPtr->ions_c = c_ions_c;
//    }
//    tempPrsmArgumentPtr->get_match_peaks();
    ///原函数end

    ///test by lyc 2024.04.10
//    cout<<"###n_ptms###"<<endl;
//    for(auto x:n_ptms){
//        cout<<x.first<<"  "<<x.second<<"  "<<x.mod_mass<<endl;
//    }
//    cout<<"###c_ptms###"<<endl;
//    for(auto x:c_ptms){
//        cout<<x.first<<"  "<<x.second<<"  "<<x.mod_mass<<endl;
//    }


    double th = 1.1;
    int pro_lenth = tempPrsmArgumentPtr->theory_ions_mass_n.size();
    vector<double> N_bIons_mod(pro_lenth,0);
    //Bi离子所对应的修饰值为N_bIons_mod[i-1]，即：i == thoe_id
    for(auto x:n_ptms){
        for(int j = x.first+1;j<N_bIons_mod.size();j++){
            N_bIons_mod[j] += x.mod_mass;
        }
    }
    //Yi离子所对应的修饰值为N_yIons_mod[i-1]，即：i == thoe_id
    //n_ptm表示之前函数由n端判定修饰得出的ptm
    vector<double> N_yIons_mod(pro_lenth,0);
    for(auto x:n_ptms){
        for(int j = pro_lenth-x.second-1;j<N_yIons_mod.size();j++){
            N_yIons_mod[j] += x.mod_mass;
        }
    }
    ///过滤违规被修饰离子
    vector<node> sel_n_ions_n,sel_n_ions_c;
    for(auto x:n_ions_n){
        if(fabs(x.index - N_bIons_mod[x.thoe_id])<th)
            sel_n_ions_n.push_back(x);
    }
    for(auto x:n_ions_c){
        if(fabs(x.index - N_yIons_mod[x.thoe_id])<th)
            sel_n_ions_c.push_back(x);
    }

    //由C端判定修饰时b离子对应的修饰值
    vector<double> C_bIons_mod(pro_lenth,0);
    for(auto x:c_ptms){
        for(int j = x.first+1;j<C_bIons_mod.size();j++){
            C_bIons_mod[j] += x.mod_mass;
        }
    }
    vector<double> C_yIons_mod(pro_lenth,0);
    for(auto x:c_ptms){
        for(int j = pro_lenth-x.second-1;j<C_yIons_mod.size();j++){
            C_yIons_mod[j] += x.mod_mass;
        }
    }
    //过滤违规被修饰离子
    vector<node> sel_c_ions_n,sel_c_ions_c;
    for(auto x:c_ions_n){
        if(fabs(x.index - C_bIons_mod[x.thoe_id])<th)
            sel_c_ions_n.push_back(x);
    }
    for(auto x:c_ions_c){
        if(fabs(x.index - C_yIons_mod[x.thoe_id])<th)
            sel_c_ions_c.push_back(x);
    }


    if (tempPrsmArgumentPtr->var_mod_mass != 0) {
        if (sel_n_ions_n.size() + sel_n_ions_c.size() < sel_c_ions_n.size() + sel_c_ions_c.size()) {
            tempPrsmArgumentPtr->ions_n = sel_c_ions_n;
            tempPrsmArgumentPtr->ions_c = sel_c_ions_c;
            tempPrsmArgumentPtr->ptms = c_ptms;
        } else {
            tempPrsmArgumentPtr->ions_n = sel_n_ions_n;
            tempPrsmArgumentPtr->ions_c = sel_n_ions_c;
            tempPrsmArgumentPtr->ptms = n_ptms;
        }
    } else if (tempPrsmArgumentPtr->var_mod_mass == 0) {
        tempPrsmArgumentPtr->ions_n = sel_n_ions_n;
        tempPrsmArgumentPtr->ions_c = sel_c_ions_c;
    }
    tempPrsmArgumentPtr->get_match_peaks();


    N_bIons_mod.clear();
    N_yIons_mod.clear();
    C_bIons_mod.clear();
    C_yIons_mod.clear();
    sel_c_ions_n.clear();
    sel_c_ions_c.clear();
    sel_n_ions_n.clear();
    sel_n_ions_c.clear();
}


void Get_C_Terminal_Ions_And_Change_Ptms_Inter_v2(vector<modi> &Terminal_Inter, vector<node> &merg_ions,std::map<long, long> &mapDup,
                                                  vector<node> &get_ions,vector<modi> &ptmsTerminal, double ptm_mass)
{
    if (Terminal_Inter.empty())
    {
        return;
    }
    map<long, long> mapD;       //用于去重
    bool firstInt = false;      // 是否是第一个区间
    //遍历每一个可能的离子,取可能的离子
    for (int xi = 0; xi < merg_ions.size(); xi++)
    {
        long first = merg_ions[xi].thoe_id;
        long second = merg_ions[xi].mono_id;
        if (mapDup.count(second)) //如果该质谱峰已有匹配项,过滤
        {
            continue;
        }
        double maxMass = 0;    //记录区间的修饰质量
        int pLoc = 0;        //pLoc [0,Terminal_Inter.size() - 1 ],标识在第几个修饰区间
        //得出该理论峰的区间位置以及修饰质量限值
        for (int ji = 0; ji < Terminal_Inter.size(); ji++)
        {
            if (first <= Terminal_Inter[ji].first)  //实现选取峰值时，前开后闭的选取原则
            {
                break;
            }
            maxMass += Terminal_Inter[ji].mod_mass; //如果在下一个区间，则加上下一个区间的修饰质量
            pLoc = ji;
        }
//        int curInterFirst = Terminal_Inter[pLoc].first ;    //当前区间的起始位置
        if (pLoc == 0) //如果在第一个区间内
        {
            if (fabs(merg_ions[xi].index) < 2 && !mapD.count(merg_ions[xi].thoe_id) && !firstInt) //可以取 0
            {
                get_ions.push_back(merg_ions[xi]);        //将该离子选入
                mapD[merg_ions[xi].thoe_id] = merg_ions[xi].mono_id;    //将理论峰与质谱峰添加进map 去重
                //更新修饰区间
                if (merg_ions[xi].thoe_id > ptmsTerminal[pLoc].first && merg_ions[xi].thoe_id < ptmsTerminal[pLoc].second)
                {
                    ptmsTerminal[pLoc].first = merg_ions[xi].thoe_id;
                }
//                if (merg_ions[xi].thoe_id > Terminal_Inter[pLoc].first &&  merg_ions[xi].thoe_id < Terminal_Inter[pLoc].second)
//                {
//                    Terminal_Inter[pLoc].first = merg_ions[xi].thoe_id;
//                }
            } else if (fabs(merg_ions[xi].index - maxMass) < 3) {       //如果是修饰离子
                get_ions.push_back(merg_ions[xi]);      //取该离子
                mapD[merg_ions[xi].thoe_id] = merg_ions[xi].mono_id;    //去重
                firstInt = true;        //不能够再取 0 离子
                //更新修饰区间
                if (merg_ions[xi].thoe_id > ptmsTerminal[pLoc].first && merg_ions[xi].thoe_id < ptmsTerminal[pLoc].second)
                {
                    ptmsTerminal[pLoc].first = merg_ions[xi].thoe_id;
                }
//                if (merg_ions[xi].thoe_id > Terminal_Inter[pLoc].first && merg_ions[xi].thoe_id < Terminal_Inter[pLoc].second)
//                {
//                    Terminal_Inter[pLoc].first = merg_ions[xi].thoe_id;
//                }
            }
        } else if (fabs(merg_ions[xi].index - maxMass) < 3) {    //如果在后续区间，且为修饰离子
            get_ions.push_back(merg_ions[xi]);
            mapD[merg_ions[xi].thoe_id] = merg_ions[xi].mono_id;
            //则更新修饰区间
            if (merg_ions[xi].thoe_id > ptmsTerminal[pLoc].first && merg_ions[xi].thoe_id < ptmsTerminal[pLoc].second)
            {
                ptmsTerminal[pLoc].first = merg_ions[xi].thoe_id;
            }
//            if (merg_ions[xi].thoe_id > Terminal_Inter[pLoc].first && merg_ions[xi].thoe_id < Terminal_Inter[pLoc].second)
//            {
//                Terminal_Inter[pLoc].first = merg_ions[xi].thoe_id;
//            }
        }
    }//For merg_ions END
}

void prsm_ions_location::ions_location_match_N_v2(vector<node> &n_ions, vector<node> &c_ions,vector<node> &merg_ions_n,
                                                  vector<node> &merg_ions_c,vector<modi> &mod,const vector<double> &modifyTable,
                                                  int lenth, double modifier_mass)
{
    vector<modi> N_Terminal_Inter;    //计算N端修饰区间
    vector<modi> C_Terminal_Inter;    //计算C端修饰区间
    std::map<long, long> mapDup;      //first保存质谱峰下标，second保存理论峰下标，用于去重
    int log = lenth - 1;
    //最长路径算法
    mylib::ions_path_Ptr lpCP = new mylib::ions_path();
    lpCP->initializeData(merg_ions_n, merg_ions_c, modifyTable, n_ions, lenth);
    delete lpCP;
    if (n_ions.size() == 0)
    {
        return ;
    }
    // 计算修饰区间 （node -> modi 的转变)
    mylib::etd_process::index_Unknown_Ptm(n_ions,N_Terminal_Inter, modifyTable);
    // 2022.9.2 添加 -- 进行修饰质量的判断，假设缺少的修饰质量在后半段
    mylib::etd_process::assumMissingMass(n_ions,N_Terminal_Inter,modifyTable,modifier_mass,lenth);
    //2022-11-20 修正修饰
//  转化为另外一端可取离子区间
    n_terminal_transaction_c_terminal(N_Terminal_Inter, C_Terminal_Inter, lenth);
//  根据另外一端修饰区间取离子，并且改变修饰区间，区间的取值是前开后闭
    vector<modi> ptmsTerminal = C_Terminal_Inter ;  //变动后的区间
    Get_C_Terminal_Ions_And_Change_Ptms_Inter_v2(C_Terminal_Inter, merg_ions_c, mapDup, c_ions, ptmsTerminal,modifier_mass);
    //由最终的 terminal-Inter 得出最终修饰区间
    log = lenth - 2 ;
    for (int xi = ptmsTerminal.size() - 1; xi >= 0; xi--) {
        modi t;
        t.first = log - ptmsTerminal[xi].second;
        t.second = log - ptmsTerminal[xi].first;
        t.mod_mass = ptmsTerminal[xi].mod_mass;
        mod.push_back(t);
    }
    return;
}

void prsm_ions_location::ions_location_match_C_v2(vector<node> &n_ions, vector<node> &c_ions,vector<node> &merg_ions_n, vector<node> &merg_ions_c,
                                                  vector<modi> &mod,const vector<double> &modifyTable,int lenth, double modifier_mass)
{
    int i = 0, j = 0;
    vector<modi> N_Terminal_Inter;    //计算N端修饰区间
    vector<modi> C_Terminal_Inter;    //计算C端修饰区间
    std::map<long, long> mapDup;      //first保存质谱峰下标，second保存理论峰下标，用于去重
    int log = lenth - 1;
    //求出最长路径
    mylib::ions_path_Ptr lpCP = new mylib::ions_path();
    lpCP->initializeData(merg_ions_c, merg_ions_n, modifyTable, c_ions, lenth);
    delete lpCP;
    //主要用于利用n端离子制造修饰区间( node -> modi)
    if (c_ions.size() == 0) {
        return ;
    }
    //得出相应的修饰区间
    mylib::etd_process::index_Unknown_Ptm(c_ions,C_Terminal_Inter, modifyTable);
    //2022.9.2 添加 -- 进行修饰质量的判断，假设缺少的修饰质量在后半段
    mylib::etd_process::assumMissingMass(c_ions,C_Terminal_Inter,modifyTable,modifier_mass,lenth);
    //由该端修饰区间转化为另一端修饰区间
    n_terminal_transaction_c_terminal(C_Terminal_Inter, N_Terminal_Inter, lenth);
    //根据修饰区间，另外一端所有离子，取出另外一端的离子
    vector<modi> ptmsTerminal  = N_Terminal_Inter;
    Get_C_Terminal_Ions_And_Change_Ptms_Inter_v2(N_Terminal_Inter, merg_ions_n, mapDup, n_ions, ptmsTerminal,modifier_mass);

    //-------------------------------由最终的N-terminal-Inter 得出最终修饰区间
    mod = ptmsTerminal;
    return;
}