#include "etd_process.hpp"


void mylib::etd_process::index_ions_Unknown(vector<double> &sub, vector<double> &index_y,
                                            std::vector<double> &modifier,
                                            double &modifier_mass)    //定位修饰偏差值
{
    int i = 0;
    for (i = 0; i < sub.size(); ++i) {
        //判断该值小于修饰值总质量，同时为一定修饰组合的质量
        if ((sub[i] <= modifier_mass + 5) && (mylib::data_stream::find_key(sub[i], modifier) != -1)) {
            index_y.push_back(sub[i]);
            cout << " index = " << sub[i] << endl;
        } else if (fabs(sub[i]) <= 2) {      //是否该值为1以下？
            index_y.push_back(sub[i]);
        } else {                         //既不是质量偏移峰，也不是正确的峰
            index_y.push_back(INT_MAX);
        }
    }
    //------------------过滤一下----------------//
    //有些离子的组合，无法到修饰质量 . 例如 14*6 = 74 + ？ = 98
    for (i = 0; i < index_y.size(); ++i) {
        int temp;
        if (index_y[i] != -1 && index_y[i] < modifier_mass - 2) {
            temp = mylib::data_stream::remove_false(index_y[i], modifier, modifier_mass);
            if (temp == 0) {
                index_y[i] = INT_MAX;
            }
        }
    }
}

/**
 * 进行偏差的过滤，主要的思路是：先判定该值
 * @param sub
 * @param index_y
 * @param modifier
 * @param modifier_mass
 */

void mylib::etd_process::index_ions(vector<double> &sub, vector<double> &index_y,
                                    std::vector<double> &modifier,
                                    double &modifier_mass)    //定位修饰偏差值
{
    int i = 0;
    for (i = 0; i < sub.size(); ++i) {
        //判断该值小于修饰值总质量，同时为一定修饰组合的质量
        if ((sub[i] <= modifier_mass + 5) && (mylib::data_stream::find_key(sub[i], modifier) != -1)) {
            index_y.push_back(sub[i]);
        } else if (fabs(sub[i]) <= 2) {      //是否该值为1以下？
            index_y.push_back(sub[i]);
        } else {                         //既不是质量偏移峰，也不是正确的峰
            index_y.push_back(INT_MAX);
        }
    }
    //------------------过滤一下----------------//
    //有些离子的组合，无法到修饰质量 . 例如 14*6 = 74 + ？ = 98
    for (i = 0; i < index_y.size(); ++i) {
        int temp;
        if (index_y[i] != -1 && index_y[i] < modifier_mass - 2) {
            temp = mylib::data_stream::remove_false(index_y[i], modifier, modifier_mass);
            if (temp == 0) {
                index_y[i] = INT_MAX;
            }
        }
    }
}

/**
 * 9.20 重写了满足条件的定位离子峰的方法
 *           //1、当两者相差不大
            //2、如果相差较大，是否满足ptmMassShift
            //2.1 当是 modifier_mass 的值，直接加入
            //2.2 当是 modifier_mass 的子值，判断是否有可达modifier_mass的质量，如果有，则加入，没有，则不加入
 * @param sub
 * @param index
 * @param modifier
 * @param modifier_mass
 */
void mylib::etd_process::indexIonsReWrite(vector<double> &sub, vector<double> &index,
                                          std::vector<double> &modifier, double &modifier_mass)
{
    double sex = 5 ;
    for (int i = 0; i < sub.size(); ++i) {
        if (fabs(sub[i]) < sex) {
            index.push_back(sub[i]);
        } else if (fabs(sub[i] - modifier_mass) < 5) {
            index.push_back(sub[i]);
        } else {
            double searchIndex = mylib::data_stream::txpb_binarey_search_ex(modifier, modifier.size(), sub[i]);
            if (searchIndex < 0 ) continue;
            double temp =  modifier[searchIndex] - sub[i] ;
            if (fabs(temp) < sex) {
                index.push_back(sub[i]);
            } else {
                index.push_back(INT_MAX) ;
            }
        }
    }
}

bool mylib::etd_process::_PPM_H(double t, double m, double &minPPM) {
    ///原ppm考虑多氢少氢
//    double uP = (t - m)
//                / (t) * 1000000;
//    double uP_uH = (t - 1.007276 - m)
//                   / (t - 1.007276) * 1000000;
//    double uP_dH = (t + 1.007276 - m)
//                   / (t + 1.007276) * 1000000;
//    if (fabs(uP) < 15 || fabs(uP_uH) < 15 || fabs(uP_dH) < 15) {
//        minPPM = fabs(uP) > fabs(uP_uH) ? (fabs(uP_uH) > fabs(uP_dH) ? uP_dH : uP_uH) : fabs(fabs(uP) > fabs(uP_dH) ? uP_dH : uP);
//        return true;
//    }
//    return false;

    //不考虑多氢少氢，by lyc
//    double uP = (t - m)
//                / (t) * 1000000;
//    if (fabs(uP) < 15 ) {
//        minPPM = fabs(uP);
//        return true;
//    }
//    return false;

    ///考虑同位素峰
    double uP = (t - m)
                / (t) * 1000000;
    double uP_uH = (t - 1.0024 - m)
                   / (t - 1.0024) * 1000000;
    double uP_dH = (t + 1.0024 - m)
                   / (t + 1.0024) * 1000000;
    if (fabs(uP) < 15 || fabs(uP_uH) < 15 || fabs(uP_dH) < 15) {
        minPPM = fabs(uP) > fabs(uP_uH) ? (fabs(uP_uH) > fabs(uP_dH) ? uP_dH : uP_uH) : fabs(fabs(uP) > fabs(uP_dH) ? uP_dH : uP);
        return true;
    }
    return false;

}


void mylib::etd_process::GetPpmUnknownReWrite(vector<double> &reference, vector<double> &peer,
                                              std::vector<std::pair<long, long> > &alignment,
                                              std::vector<double> &index, std::vector<double> &modifier, std::vector<double> &ppm,
                                              double &modifier_mass) {
    int i = 0;
    for (i = 0; i < index.size(); i++) {
        int searchIndex = mylib::data_stream::txpb_binarey_search_ex(modifier, modifier.size(), index[i]);   //寻找一个最接近index的修饰质量
        if (searchIndex < 0){
            continue;
        }
        double massShift = modifier[searchIndex] ;
        double theoryMass = reference[alignment[i].first];
        double monoMass = peer[alignment[i].second] ;
        double minPPM = INT_MAX ;
        if (mylib::etd_process::_PPM_H(theoryMass + massShift, monoMass, minPPM)) {
            ppm.push_back(minPPM);
        } else {
            ppm.push_back(INT_MAX);
        }
    }
}


void mylib::etd_process::Get_ppm(vector<double> &reference, vector<double> &peer, std::vector<std::pair<long, long> > &alignment,
                                 std::vector<double> &index, std::vector<double> &modifier, std::vector<double> &ppm,
                                 double &modifier_mass) {
    int i = 0;
    double tol_ppm = 0.00003;   //  扩大一点点
    for (i = 0; i < index.size(); i++) {
        if (index[i] != INT_MAX && index[i] > 12)            //index为修饰值时
        {
            int k = mylib::data_stream::find_value(index[i], modifier, modifier_mass);
            if (k != -1) {
                double sub = peer[alignment[i].second] - reference[alignment[i].first];
                double tol_mass = (reference[alignment[i].first] + modifier[k]) * tol_ppm;        //容错质量
                //当差值为正，即质谱位移往右偏差(两种情况，1.差一丢丢 2.容错设计。3.1.0024的设计)
                if (sub > modifier[k]) {
                    double pm = fabs((reference[alignment[i].first] + modifier[k] + tol_mass) -
                                     (peer[alignment[i].second])) /
                                (reference[alignment[i].first] + modifier[k] + tol_mass) * 1000000;
                    double pm2 = fabs((reference[alignment[i].first] + modifier[k]) - (peer[alignment[i].second])) /
                                 (reference[alignment[i].first] + modifier[k]) * 1000000;
                    double pm3 =
                            fabs((reference[alignment[i].first] + modifier[k] + 1.0024) - (peer[alignment[i].second])) /
                            (reference[alignment[i].first] + modifier[k]) * 1000000;
                    double t = pm > pm2 ? (pm2 > pm3 ? pm3 : pm2) : (pm > pm3 ? pm3 : pm);
                    ppm.push_back(t);
                } else {       //当差值为负时，质量向左偏。（可能少H，可能仅仅只是偏差一点）
                    double pm = fabs((reference[alignment[i].first] + modifier[k] + tol_mass) -
                                     (peer[alignment[i].second] + 1.007276)) /
                                (reference[alignment[i].first] + modifier[k] + tol_mass) * 1000000;    //少H的容错设计
                    double pm0 = fabs((reference[alignment[i].first] + modifier[k]) -
                                      (peer[alignment[i].second] + 1.007276)) /
                                 (reference[alignment[i].first] + modifier[k]) * 1000000;    //少H
                    double pm1 = fabs((reference[alignment[i].first] + modifier[k]) - (peer[alignment[i].second])) /
                                 (reference[alignment[i].first] + modifier[k]) * 1000000;        //偏差一点点
                    double pm2 = fabs((reference[alignment[i].first] + modifier[k] + tol_mass) -
                                      (peer[alignment[i].second])) /
                                 (reference[alignment[i].first] + modifier[k] + tol_mass) * 1000000;    //容错设计
                    double t = pm < pm1 ? (pm < pm2 ? pm : pm2) : (pm1 < pm2 ? pm1 : pm2);
                    if (pm0 > t) {
                        ppm.push_back(t);
                    } else {
                        ppm.push_back(pm0);
                    }
                }
            } else {
                ppm.push_back(INT_MAX);
            }
        } else if (index[i] != INT_MAX && fabs(index[i]) < 2) {        //当index值为小于1时
            if (index[i] > 0)                //当质量为右偏移，那么只进行容错的设计
            {
                double tol_mass = reference[alignment[i].first] * tol_ppm;
                double pm = fabs(reference[alignment[i].first] - peer[alignment[i].second]) /
                            (reference[alignment[i].first]) * 1000000;
                double pm_h = fabs(reference[alignment[i].first] + 1.0024 - peer[alignment[i].second]) /
                              (reference[alignment[i].first]) * 1000000;
                double pm_tol_add = fabs(reference[alignment[i].first] + tol_mass - peer[alignment[i].second]) /
                                    (reference[alignment[i].first]) * 1000000;
                double pm_tol_sub = fabs(reference[alignment[i].first] - tol_mass - peer[alignment[i].second]) /
                                    (reference[alignment[i].first]) * 1000000;
                double t = pm_tol_add > pm_tol_sub ? (pm_tol_sub > pm ? pm : pm_tol_sub) : (pm_tol_add > pm ? pm
                                                                                                            : pm_tol_add);
                if (pm_h > t)
                    ppm.push_back(t);
                else {
                    ppm.push_back(pm_h);
                }
            }
            if (index[i] < 0) {
                double tol_mass = reference[alignment[i].first] * tol_ppm;
                double pm = fabs((reference[alignment[i].first] + tol_mass) - (peer[alignment[i].second] + 1.007276)) /
                            (reference[alignment[i].first] + tol_mass) * 1000000;    //少H的容错设计
                double pm0 = fabs((reference[alignment[i].first]) - (peer[alignment[i].second] + 1.007276)) /
                             (reference[alignment[i].first]) * 1000000;    //少H
                double pm1 = fabs((reference[alignment[i].first]) - (peer[alignment[i].second])) /
                             (reference[alignment[i].first]) * 1000000;        //偏差一点点
                double pm2 = fabs((reference[alignment[i].first] + tol_mass) - (peer[alignment[i].second])) /
                             (reference[alignment[i].first] + tol_mass) * 1000000;    //容错设计
                double t = pm < pm1 ? (pm < pm2 ? pm : pm2) : (pm1 < pm2 ? pm1 : pm2);
                if (pm0 > t) {
                    ppm.push_back(t);
                } else {
                    ppm.push_back(pm0);
                }
            }
        } else {
            ppm.push_back(INT_MAX);
        }
    }
}

//参数1，保存最终离子，参数2保存所有的可行离子，参数三，为修饰离子
void Greed_Get_Ions_strategy(vector<node> &c_n_ions, vector<node> &merg_ions_c_n,
                             vector<node> &mod_ions)    //--------前两个个为输入，后一个为输出
{
    vector<node>::iterator it = merg_ions_c_n.begin();
    vector<node> mod_ions_re;    //--------用于临时保存
    //----取出N端完全匹配的离子
    while (fabs(it->index) < 2 && it != merg_ions_c_n.end()) {
        c_n_ions.push_back(*it);
        it++;
    }
    //----贪婪策略取修饰值
    for (vector<node>::iterator itk = it; itk != merg_ions_c_n.end(); itk++) {
        for (vector<node>::iterator it_m = itk; it_m != merg_ions_c_n.end(); it_m++) {
            //取出第一个修饰离子
            if (mod_ions.size() == 0) {
                mod_ions.push_back(*it_m);
                continue;
            }
            //取后面的修饰值
            if (it_m->thoe_id != mod_ions.back().thoe_id) {
                if (it_m->index + 3 < mod_ions.back().index) {
                    continue;    //太小了
                } else {            //取得值必然大于或等于最后一个修饰
                    mod_ions.push_back(*it_m);
                }
            } else if (it_m->thoe_id == mod_ions.back().thoe_id) {            //取同位素或者中间多的修饰
                if (it_m->index + 3 < mod_ions.back().index) {    //小于
                    continue;
                } else if (fabs(it_m->index - mod_ions.back().index) < 2) {        //等于
                    mod_ions.push_back(*it_m);
                } else if (mod_ions.size() != 1) {    //大于
                    if (mod_ions.size() >= 2) {
                        if (fabs(mod_ions[mod_ions.size() - 2].index - mod_ions[mod_ions.size() - 1].index) < 3) {
                            continue;
                        } else {
                            mod_ions.pop_back();
                        }
                    }
                    mod_ions.push_back(*it_m);
                }
            }
        }
        if (mod_ions.size() > mod_ions_re.size()) {
            mod_ions_re = mod_ions;
        }
        mod_ions.clear();
    }

    for (int xi = 0; xi < mod_ions_re.size(); ++xi) {
        c_n_ions.push_back(mod_ions_re[xi]);
    }

    mod_ions_re.swap(mod_ions);
}

/**
 * 取对应的互补离子峰为
 * @param N_Terminal_Inter
 * @param C_Terminal_Inter
 * @param log
 */
void N_Terminal_Inter_transaction_C_Terminal_Inter(vector<modi> &N_Terminal_Inter, vector<modi> &C_Terminal_Inter, int log) {
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

//得到另外一端的离子，并且改变Ptms区间
//根据 teriminal_ inter , 将merg_ions 中的离子选入 get_ions
void Get_C_Terminal_Ions_And_Change_Ptms_Inter(vector<modi> &Terminal_Inter, vector<node> &merg_ions,std::map<long, long> &mapDup,
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

bool testModTerminal(const vector<modi> & terminal,const vector<double> & massTable,
                     double modifyMass)
{
    double m = 0;
    for (int i = 0; i < terminal.size(); i++) {
        m += terminal[i].mod_mass;
    }
    if ( fabs(m - modifyMass) > 2) {
        return true ;
    }
    return false ;
}

/**
 * 主要进行缺失修饰质量的假设设定
 * @param N_Terminal_Inter
 * @param modifyTable
 * @param mMass
 * @param seqLen
 */
void mylib::etd_process::assumMissingMass(vector<node> &ions, vector<modi> &N_Terminal_Inter,
                                          const std::vector<double> &modifyTable, double mMass , int seqLen) {
    double tMass = 0 ;
    int first = 0 ;
    for (int i = ions.size()-1; i >= 0 ; i--) {
        first = ions[i].thoe_id;
        break ;
    }
    for (int i = 0; i < N_Terminal_Inter.size(); ++i) {
        tMass += N_Terminal_Inter[i].mod_mass ;
    }
    double missingMass = mMass - tMass ;
    int index = mylib::speculate_lib::txpb_binarey_search_ex(modifyTable, modifyTable.size(), mMass - tMass) ;
    //如果存在修饰
    if ( fabs(modifyTable[index]) > 0.000001 )
    {
        if (first == seqLen - 1)
        {
            //去除最后的一个离子
            vector<node> tempIons ;
            for (node t : ions) {
                if (t.thoe_id != seqLen - 1) {
                    tempIons.push_back(t) ;
                    first = t.thoe_id;
                } else {
                    break;
                }
            }
            ions = tempIons ;
            //判断当前区间是否存在
            int lastTerminalSecond = -1;
            int indexT = -1 ;
            for (int i = N_Terminal_Inter.size()-1; i >=0 ; --i) {
                lastTerminalSecond = N_Terminal_Inter[i].second ;
                indexT = i ;
                break ;
            }
            if (lastTerminalSecond == seqLen - 1) {       //如果区间存在
                N_Terminal_Inter[indexT].mod_mass += modifyTable[index];
            } else {
                modi m ;
                m.first = first ;
                m.second = seqLen-1 ;
                m.mod_mass = modifyTable[index] ;
                N_Terminal_Inter.push_back(m) ;
            }
        } else {
            modi m ;
            m.first = first ;
            m.second = seqLen-1 ;
            m.mod_mass = modifyTable[index] ;
            N_Terminal_Inter.push_back(m) ;
        }
    }
}

void mylib::etd_process::ions_location_match_N(vector<node> &n_ions, vector<node> &c_ions, vector<node> &merg_ions_n,
                                               vector<node> &merg_ions_c, vector<modi> &mod, const vector<double> &modifyTable,
                                               int lenth, double modifier_mass)
{
    vector<modi> N_Terminal_Inter;    //计算N端修饰区间
    vector<modi> C_Terminal_Inter;    //计算C端修饰区间
    std::map<long, long> mapDup;      //first保存质谱峰下标，second保存理论峰下标，用于去重
    int log = lenth - 1;
    //最长路径算法
    mylib::ions_path_Ptr lpN = new mylib::ions_path();
    lpN->initializeData(merg_ions_n, merg_ions_c, modifyTable, n_ions, lenth);
    delete lpN;
    if (n_ions.size() == 0)
    {
        return ;
    }
//    e::proc::N_terminal_Get_Varible_Interval(c_n_ions, N_Terminal_Inter, mapDup, modify_table, modifier_mass, log);
    // 计算修饰区间 （node -> modi 的转变)
    mylib::etd_process::index_Unknown_Ptm(n_ions, N_Terminal_Inter, modifyTable);
    // 2022.9.2 添加 -- 进行修饰质量的判断，假设缺少的修饰质量在后半段
    mylib::etd_process::assumMissingMass(n_ions, N_Terminal_Inter, modifyTable, modifier_mass, lenth);
    //2022-11-20 修正修饰
//    for (int i = 0 ; i < N_Terminal_Inter.size(); ++i)
//    {
//        int second = N_Terminal_Inter[i].first ;
//        if (i > 0 && second == N_Terminal_Inter[i-1].second) {
//            N_Terminal_Inter[i].first += 1 ;
//        }
//    }
//    cout <<" =========== "<<endl;
//    for (int i = 0; i < N_Terminal_Inter.size(); i++) {
//        cout << " N_varible : first = " << N_Terminal_Inter[i].first
//        << " second = " << N_Terminal_Inter[i].second
//        << " mass = " << N_Terminal_Inter[i].mod_mass << endl;
//    }
//  转化为另外一端可取离子区间
    N_Terminal_Inter_transaction_C_Terminal_Inter(N_Terminal_Inter, C_Terminal_Inter, lenth);
//     cout<<" Transaction Interval "<<endl;
//     for (int i = 0 ; i < C_Terminal_Inter.size() ; i++) {
//     	cout<<" C_varible : first = "<<C_Terminal_Inter[i].first<<" second = "<<C_Terminal_Inter[i].second<<" mass = " <<C_Terminal_Inter[i].mod_mass<<endl;
//     }
//  根据另外一端修饰区间取离子，并且改变修饰区间，区间的取值是前开后闭
    vector<modi> ptmsTerminal = C_Terminal_Inter ;  //变动后的区间
    Get_C_Terminal_Ions_And_Change_Ptms_Inter(C_Terminal_Inter, merg_ions_c, mapDup, c_ions, ptmsTerminal,modifier_mass);
//     cout<<"C ions : \n";
//     for(int xi = 0 ; xi < c_ions.size() ; xi++ ) {
//         cout<<"thero id = "<<c_ions[xi].thoe_id<<", mono id = "<<c_ions[xi].mono_id<<" index = "<<c_ions[xi].index<<endl;
//     }
//    for (int i = 0 ; i < ptmsTerminal.size() ; i++) {
//        cout<<" C_varible : first = "<<ptmsTerminal[i].first<<" second = "<<ptmsTerminal[i].second<<" mass = " <<ptmsTerminal[i].mod_mass<<endl;
//    }
    /**
     * 由最终的 terminal-Inter 得出最终修饰区间
     */
     log = lenth - 2 ;
    for (int xi = ptmsTerminal.size() - 1; xi >= 0; xi--) {
        modi t;
        t.first = log - ptmsTerminal[xi].second;
        t.second = log - ptmsTerminal[xi].first;
        t.mod_mass = ptmsTerminal[xi].mod_mass;
        mod.push_back(t);
    }
//    for (modi t : mod) {
//        cout <<" first = " << t.first
//        <<" second = "<< t.second
//        << " mass = " << t.mod_mass
//        <<endl;
//    }
    return;
}

bool judge_Peaks(const vector<double> &mass, double pre, double next) {
    int n = mylib::data_stream::txpb_binarey_search_ex(mass, mass.size(), next - pre);
    int n_dH = mylib::data_stream::txpb_binarey_search_ex(mass, mass.size(), next - 1.007276 - pre);
    int n_aH = mylib::data_stream::txpb_binarey_search_ex(mass, mass.size(), next + 1.007276 - pre);
    if (n != -1 && n_aH != -1 && n_aH != -1) {
        if (mass[n] == 0 || mass[n_dH] == 0 || mass[n_aH] == 0) {
            return true;
        }
    }
    return false;
}


/**
 *
 * @param longIons_N
 * @param longIons_C
 * @param mod
 * @param PtmMass
 * @param onePtmMass
 */
void mylib::etd_process::index_Unknown_Ptm(vector<node> &longIons_N, vector<modi> &mod, const vector<double> &PtmMass) {
    vector<modi> modIndex_N;
    vector<int> indexIons;

    int px = 0;
    int index;
     if (longIons_N.size() == 0) {
         return ;
     }
     indexIons.push_back(px);
    // 先标识出不同的修饰
    for (int i = 0, j = 1; i < longIons_N.size() && j < longIons_N.size(); ++i, ++j) {
        if (judge_Peaks(PtmMass, longIons_N[i].index, longIons_N[j].index)) {
            indexIons.push_back(px);
            continue;
        }
        px++;
        indexIons.push_back(px);
    }
    // 粗糙定位修饰
    for (int i = 0; i < indexIons.size(); ++i)
    {
        index = mylib::data_stream::txpb_binarey_search_ex(PtmMass, PtmMass.size(), longIons_N[i].index);
        if (i == 0 && index != -1 && PtmMass[index] != 0)
        {
            modi t;
            t.first = -1;
            t.second = longIons_N[i].thoe_id;
            t.mod_mass = PtmMass[index];
            modIndex_N.push_back(t);
        }
        if (i > 0 && indexIons[i] != indexIons[i - 1])
        {
            index = mylib::data_stream::txpb_binarey_search_ex(PtmMass, PtmMass.size(), longIons_N[i].index - longIons_N[i - 1].index);
            modi t;
            t.first = longIons_N[i - 1].thoe_id;
            t.second = longIons_N[i].thoe_id;
            t.mod_mass = PtmMass[index];
            modIndex_N.push_back(t);
        }
    }
    mod = modIndex_N;
    return;
}

/**
 *
 * @param n_ions   选出n端的离子
 * @param c_ions   选出C端的离子
 * @param merg_ions_n n端所有符合条件的离子
 * @param merg_ions_c c端所有符合条件的离子
 * @param alignment_n n端对齐
 * @param alignment_c c端对齐
 * @param reference 理论峰
 * @param peer      质谱峰
 * @param index_c_n n端差值
 * @param index_n_c c端差值
 * @param ppm_c_n   n端ppm值
 * @param ppm_n_c   c端ppm值
 * @param mod       最终的修饰区间
 * @param modify_table 修饰质量向量
 * @param lenth     长度
 * @param modifier_mass  总修饰质量
 */
void mylib::etd_process::unknown_Ptm_Search(vector<node> &n_ions, vector<node> &c_ions,
                                            vector<node> &merg_ions_n, vector<node> &merg_ions_c,
                                            vector<pair<long, long> > &alignment_n,
                                            vector<pair<long, long> > &alignment_c,
                                            vector<double> &reference, vector<double> &peer,
                                            vector<double> &index_c_n, vector<double> &index_n_c,
                                            vector<double> &ppm_c_n, vector<double> &ppm_n_c,
                                            vector<modi> &mod,
                                            const vector<double> &modify_table,
                                            const map<double, string> &modMap,
                                            int lenth, double modifier_mass) {
    int i = 0, j = 0;
    vector<modi> N_Terminal_Inter;    //计算N端修饰区间
    vector<modi> C_Terminal_Inter;    //计算C端修饰区间
    std::map<long, long> mapDup;      //first保存质谱峰下标，second保存理论峰下标，用于去重
    int log = lenth - 1;
    /**
     * 0、先获取未知修饰和已知修饰
     */
    vector<double> onePtmMass;
    for (map<double, string>::const_iterator it = modMap.begin(); it != modMap.end(); it++) {
        if (it->second.size() == 1 && it->second == "Unknown" && it->second != "0") {
            onePtmMass.push_back(it->first);
        }
    }
    /**
     * 1、选出c端最长路径 , 那么在其中定位已知修饰
     * 1.1 输入为 X 端所有可能的离子,将最长选择路径离子放入 ions 中
     * 1.2 以 N 端 为参照，先尝试确定修饰未知
     */
    mylib::ions_path_Ptr lpCP = new mylib::ions_path();
    lpCP->initializeData(merg_ions_c, merg_ions_n, modify_table, c_ions, lenth);
    delete lpCP;
    mylib::ions_path_Ptr lpCP_2 = new mylib::ions_path();
    lpCP->initializeData(merg_ions_n, merg_ions_c, modify_table, n_ions, lenth);
    delete lpCP;
    vector<modi> modN;
    index_Unknown_Ptm(n_ions, modN, modify_table);
    vector<modi> modC;
    N_Terminal_Inter_transaction_C_Terminal_Inter(modN, modC, log);
    for (i = 0; i < modC.size(); i++) {
        cout << " modC  : first = " << modC[i].first << " second = " << modC[i].second << " mass = " << modC[i].mod_mass
             << endl;
    }
    c_ions.clear();
//    Get_C_Terminal_Ions_And_Change_Ptms_Inter(modC, merg_ions_c, mapDup, c_ions, modifier_mass);
    cout << "C ions : \n";
    for (int xi = 0; xi < c_ions.size(); xi++) {
        cout << "thero id = " << c_ions[xi].thoe_id << ", mono id = " << c_ions[xi].mono_id << " index = "
             << c_ions[xi].index << endl;
    }
    mod = modC;
    for (i = 0; i < mod.size(); i++) {
        cout << " mod : first = " << mod[i].first << " second = " << mod[i].second << " mass = " << mod[i].mod_mass
             << endl;
    }
    /**
     * 2、定位已知修饰
     * 根据c_ions 选出最佳修饰位置
     */

    return;
}

void mylib::etd_process::ions_location_match_C(vector<node> &n_ions, vector<node> &c_ions, vector<node> &merg_ions_n, vector<node> &merg_ions_c,
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
     * modifyTable：修饰质量（double）数组
     */
    mylib::etd_process::index_Unknown_Ptm(c_ions, C_Terminal_Inter, modifyTable);
//    for (node t : c_ions) {
//        cout << "mod mass = " << modifier_mass
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
    N_Terminal_Inter_transaction_C_Terminal_Inter(C_Terminal_Inter, N_Terminal_Inter, lenth);
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
    Get_C_Terminal_Ions_And_Change_Ptms_Inter(N_Terminal_Inter, merg_ions_n, mapDup, n_ions, ptmsTerminal,modifier_mass);

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


void mylib::etd_process::N_terminal_Get_Varible_Interval(vector<node> &ions, vector<modi> &Interval,
                                                         map<long, long> &mapDup,
                                                         const vector<double> &modify_table,
                                                         double PtmsMass, int lenth) {
    //lenth为截断后的序列长度
    int ji = 0;
    double mass = 0;
    if (ions.size() == 0) {
        return;
    }
    if (ions.size() == 1) {
        int n = mylib::data_stream::txpb_binarey_search_ex(modify_table, modify_table.size(), ions[0].index - 0);
        if (n == -1 || modify_table[n] == 0) {
            return;
        }
        modi t;
        t.first = 0;
        t.second = ions[0].thoe_id;
        t.mod_mass = modify_table[n];
        Interval.push_back(t);
        mapDup[ions[0].mono_id] = ions[0].thoe_id;        //用来记录下质谱的碎片下标
        return;
    }
    int fms_log = mylib::data_stream::txpb_binarey_search_ex(modify_table, modify_table.size(), ions[0].index);
    if (modify_table[fms_log] != 0) {
        modi t;
        t.first = 0;
        t.second = ions[0].thoe_id;
        t.mod_mass = modify_table[fms_log];
        Interval.push_back(t);
        mass += modify_table[fms_log];
        mapDup[ions[0].mono_id] = ions[0].thoe_id;        //用来记录下质谱的碎片下标
    }
    // 离子与离子之间要么存在修饰质量偏差，要么不存在修饰质量偏差。
    vector<double> indexIons;

    for (int xi = 0; xi < ions.size() - 1; xi++) {
        ji = xi + 1;
        int n = mylib::data_stream::txpb_binarey_search_ex(modify_table, modify_table.size(), ions[ji].index - ions[xi].index);
        if (n == -1 || modify_table[n] == 0) {    // binarySearch Failed
            continue;
        }

        modi t;
        t.first = ions[xi].thoe_id;
        t.second = ions[ji].thoe_id;
        t.mod_mass = modify_table[n];
        mass += modify_table[n];
        Interval.push_back(t);
        mapDup[ions[xi].mono_id] = ions[xi].thoe_id;        //用来记录下质谱的碎片下标
    }
    // 末端存在修饰的情况
    int n = mylib::data_stream::txpb_binarey_search_ex(modify_table, modify_table.size(), PtmsMass - mass);
    if (n == -1 || modify_table[n] == 0) {
        return;
    }
    modi t_end;
    if (!Interval.empty()) {        // 如果修饰不为空
        t_end.first = Interval.back().second;
        t_end.second = lenth;
        t_end.mod_mass = modify_table[n];
    } else {                        // 如果为空
        t_end.first = 0;
        t_end.second = lenth;
        t_end.mod_mass = modify_table[n];
    }
    Interval.push_back(t_end);
}

// struct node{
//     double thoe_id; //理论序号
//     double mono_id; //质谱序号
//     double index;   //差值
//     double ppm;     //ppm值
// };

// struct modi{        //存储修饰信息
//     double first;
//     double second;
//     double mod_mass;
// };

void mylib::etd_process::Get_mod_ions_loca(vector<node> &merg_ions, vector<node> &merg_ions_other, vector<node> &ions,
                                           vector<modi> &modf,
                                           std::vector<std::pair<long, long> > &alignment, double &modifier_mass, int &lenth)
//在所有离子中取出修饰离子，并入到merg_ions中，取出来的修饰区间用modf保存

//---------注意：该函数只能运用于原数据中，必然包含修饰离子的情况-------------//
{
    int i = 0;
    node flag = ions[0];
    vector<node> temp;   //将初步取的离子数值保存。
    temp.push_back(flag);
    for (i = 1; i < ions.size(); i++)      //贪婪取值(往下搜寻，当有相等的加入，没有则加入大值)
    {
        if (ions[i].index > flag.index - 1.0024 && ions[i].index < flag.index + 1.0024) {
            temp.push_back(ions[i]);
            flag = ions[i];
            continue;
        }
        if (ions[i].index > flag.index && ions[i].thoe_id != flag.thoe_id) {
            temp.push_back(ions[i]);
            flag = ions[i];
            continue;
        }
    }
    if (merg_ions.size() != 0)                //定位第一个修饰的区间
    {
        modi r;
        r.first = merg_ions.back().thoe_id;
        r.second = temp.front().thoe_id;
        r.mod_mass = temp.front().index;
        modf.push_back(r);
    } else {
        modi r;
        r.first = 0;
        r.second = temp.front().thoe_id;
        r.mod_mass = temp.front().index;
        modf.push_back(r);
    }
    for (i = 1; i < temp.size(); i++) {
        double sub_ms = temp[i].index - temp[i - 1].index;
        if (sub_ms > 12)//相差一个修饰
        {
            modi r;
            r.first = temp[i - 1].thoe_id;
            r.second = temp[i].thoe_id;
            r.mod_mass = sub_ms;
            modf.push_back(r);
        } else if (fabs(sub_ms) < 2 && modf.size() > 0)  //如果相等
        {
            // modf[modf.size()-1].second = temp[i].thoe_id;
            // continue;
        }
        if (i == temp.size() - 1 && modf.size() > 0) {
            if (temp[i].index + 2 < modifier_mass)   //判断尾部是否存在修饰质量
            {
                modi r;
                r.first = modf[modf.size() - 1].second;
                if (merg_ions_other.size() != 0) {
                    r.second = lenth - merg_ions_other.back().thoe_id;
                } else {
                    r.second = alignment.back().first;
                }
                r.mod_mass = modifier_mass - temp[i].index;
                modf.push_back(r);
            }
        }
    }

    for (i = 0; i < temp.size(); i++)          //得到带修饰的y ions
    {
        merg_ions.push_back(temp[i]);
    }
    // for(i = 0 ; i < arr.size();i++)
    // {
    //     cout<<"arr first = "<<arr[i].first<<" arr second = "<<arr[i].second<<" arr index = "<<arr[i].mod_mass<<endl;
    // }
    // for(i = 0 ; i < modf.size();i++)
    // {
    //     cout<<"modf first = "<<modf[i].first<<" modf second = "<<modf[i].second<<" modf index = "<<modf[i].mod_mass<<endl;
    // }   
    // cout<<"NEXT "<<endl; 
}

void mylib::etd_process::Get_all_ions(vector<node> &modf_ions, vector<node> &modf_ions_other,
                                      vector<node> &ions,
                                      vector<pair<long, long> > &alignment,
                                      vector<double> &reference, vector<double> &peer,
                                      vector<double> &index,
                                      vector<double> &ppm,
                                      int &lenth)                //选出质量位移偏移为0的离子和修饰离子。	ions保存偏差为0左右的离子，modf_ions保存修饰离子
{
    node loca;
    int i = 0;
    for (i = 0; i < alignment.size(); i++)      //分别取出y离子质量偏移为0 ，和带修饰的离子。
    {
        if (modf_ions_other.size() != 0)   //检查修饰离子是否在另外一端0离子之后
        {
            if (alignment[i].first >= (lenth - modf_ions_other.back().thoe_id))
                break;
        }
        if (ppm[i] < 15 && index[i] != -1) {
            if (fabs(index[i]) < 0.015 && ppm[i] < 15)        //该y离子必然正确
            {
                if (i > 0 && alignment[i - 1].second != alignment[i].second && ppm[i - 1] < 15 &&
                    index[i - 1] != -1)//选取同位素
                {
                    double sub = fabs(reference[alignment[i - 1].second] - reference[alignment[i].second]);
                    if (sub < 1.007276 && ions[ions.size() - 1].mono_id != alignment[i - 1].second) {
                        loca.thoe_id = alignment[i - 1].first;
                        loca.mono_id = alignment[i - 1].second;
                        loca.index = index[i - 1];
                        loca.ppm = ppm[i - 1];
                        ions.push_back(loca);
                    }
                }
                loca.thoe_id = alignment[i].first;
                loca.mono_id = alignment[i].second;
                loca.index = index[i];
                loca.ppm = ppm[i];
                ions.push_back(loca);
                if (modf_ions.size() != 0)         //如果该离子必然正确，则之前的修饰离子清空。
                {
                    modf_ions.clear();
                }
            } else if (fabs(index[i]) > 0.015 && ppm[i] < 15 && fabs(index[i]) < 2)    //该离子可能正确，取决于上一条的信息
            {
                if (ions.size() != 0) {                                        //判断ions[-1]是否合法
                    if (alignment[i].first - 1 == ions[ions.size() - 1].thoe_id ||
                        modf_ions.size() <= 1)    //1.若是该离子不在修饰堆中，2.若在修饰堆中,同位素峰
                    {
                        if (i > 0 && alignment[i - 1].second != alignment[i].second && ppm[i - 1] < 15 &&
                            index[i - 1] != INT_MAX)    //选去同位素峰
                        {
                            double sub = fabs(reference[alignment[i - 1].second] - reference[alignment[i].second]);
                            if (sub < 1.007276 && ions[ions.size() - 1].mono_id != alignment[i - 1].second) {
                                loca.thoe_id = alignment[i - 1].first;
                                loca.mono_id = alignment[i - 1].second;
                                loca.index = index[i - 1];
                                loca.ppm = ppm[i - 1];
                                ions.push_back(loca);
                            }
                        }
                        loca.thoe_id = alignment[i].first;
                        loca.mono_id = alignment[i].second;
                        loca.index = index[i];
                        loca.ppm = ppm[i];
                        ions.push_back(loca);
                        if (modf_ions.size() != 0) {
                            modf_ions.clear();
                        }
                    } else { continue; }
                } else {
                    loca.thoe_id = alignment[i].first;
                    loca.mono_id = alignment[i].second;
                    loca.index = index[i];
                    loca.ppm = ppm[i];
                    ions.push_back(loca);
                    if (modf_ions.size() != 0) {
                        modf_ions.clear();
                    }
                }
            } else {                                //该值为修饰值
                if (ions.size() > 0)                    //判断究竟有没有确定的离子。
                {
                    if (alignment[i].first != ions[ions.size() - 1].thoe_id)    //如果有确定的离子则往下取修饰值
                    {
                        loca.thoe_id = alignment[i].first;
                        loca.mono_id = alignment[i].second;
                        loca.index = index[i];
                        loca.ppm = ppm[i];
                        modf_ions.push_back(loca);
                    }
                } else {
                    loca.thoe_id = alignment[i].first;
                    loca.mono_id = alignment[i].second;
                    loca.index = index[i];
                    loca.ppm = ppm[i];
                    modf_ions.push_back(loca);
                }
            }
        }
    }
    //--------2021.5.24-----------//

}

void mylib::etd_process::Get_all_ions_peaks(vector<node> &ions, vector<pair<long, long> > &alignment, vector<double> &index,
                                            vector<double> &ppm, double minPPMValue)                //ions取出所有的离子
{
    node loca;
    int i = 0;
    for (i = 0; i < alignment.size(); i++)      //分别取出y离子质量偏移为0 ，和带修饰的离子。
    {
        if (ppm[i] != INT_MAX) {
            if (ppm[i] <= minPPMValue && index[i] != INT_MAX) {
                node t;
                t.thoe_id = alignment[i].first;
                t.mono_id = alignment[i].second;
                t.index = index[i];
                t.ppm = ppm[i];
                ions.push_back(t);
            }
        }
    }
}


void mylib::etd_process::loca_ions(vector<node> &n_ions, vector<node> &c_ions, vector<locat_ions> &location,
                                   vector<double> peer) {
    int i = 0;
    for (i = 0; i < peer.size(); i++) {
        locat_ions t;
        t.c = "x";
        t.n = "x";
        location.push_back(t);
    }
    for (i = 0; i < c_ions.size(); i++) {
        location[c_ions[i].mono_id].c = "Z" + to_string(c_ions[i].thoe_id + 1);
    }
    for (i = 0; i < n_ions.size(); i++) {
        location[n_ions[i].mono_id].n = "C" + to_string(n_ions[i].thoe_id + 1);
    }

}


