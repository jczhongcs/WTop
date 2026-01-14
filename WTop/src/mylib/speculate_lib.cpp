#include "speculate_lib.h"

using namespace std;

std::map<char, double> mapi =
        {
                {'G', 57.02146},
                {'A', 71.03711},
                {'S', 87.03202},
                {'P', 97.05276},
                {'V', 99.06841},
                {'T', 101.04768},
                {'C', 103.00918},
                {'I', 113.08406},
                {'L', 113.08406},
                {'N', 114.04293},
                {'D', 115.02694},
                {'Q', 128.05858},
                {'K', 128.095},
                {'E', 129.0426},
                {'M', 131.0405},
                {'H', 137.05891},
                {'F', 147.06841},
                {'R', 156.10111},
                {'Y', 163.06333},
                {'W', 186.079354},
        };

//推测末端截断以及所带修饰质量
int mylib::speculate_lib::spect(double &c, double key, vector<double> &mod, char a[], int &cut_laction) {
    std::map<char, double>::iterator iter = mapi.begin();
    int temp = 0;
    double k;
    int i, j = 0, ac = 0;
    if (fabs(c - key) < 1) {
        return 0;
    } else {
        while (c > key) {
            iter = mapi.find(a[temp]);
            c -= iter->second;
            temp++;
        }
        double sub = key - c;  //记录差值
        k = sub - mod[0];
        while (ac != 1 && (temp != strlen(a) - 1)) {
            for (i = 1; i < mod.size(); i++) {
                if (sub - mod[i] > 0 && mod[i] != 0) {
                    if (k > sub - mod[i]) {
                        k = sub - mod[i];
                        j = i;
                    }
                    if (fabs(k) < 1) {
                        ac = 1;
                    }
                } else break;
            }
            if (ac != 1) {
                iter = mapi.find(a[temp]);
                c -= iter->second;
                sub = key - c;
                k = sub - mod[0];
                temp++;
            }
        }
        if (temp == strlen(a) - 1) {
            return -1;
        }
        cut_laction = temp;
        if (cut_laction > 30) {
            return -1;
        }
    }
    return j;
}

//推测末端截断以及所带修饰质量
int mylib::speculate_lib::spect_s(double &ThreoMass, double MonoMass,
                                  vector<double> &mod, const string &seq, int &cut_laction) {
    std::map<char, double>::iterator iter = mapi.begin();
    int temp = 0;
    double k;
    int i, j = 0, ac = 0;
    if (fabs(ThreoMass - MonoMass) < 1) {
        return 0;
    } else {
        while (ThreoMass > MonoMass) {
            iter = mapi.find(seq[temp]);
            ThreoMass -= iter->second;
            temp++;
        }
        double sub = MonoMass - ThreoMass;  //记录差值
        k = sub - mod[0];
        while (ac != 1 && (temp != seq.length() - 1)) {
            for (i = 1; i < mod.size(); i++) {
                if (sub - mod[i] > 0 && mod[i] != 0) {
                    if (k > sub - mod[i]) {
                        k = sub - mod[i];
                        j = i;
                    }
                    if (fabs(k) < 1) {
                        ac = 1;
                    }
                } else break;
            }
            if (ac != 1) {
                iter = mapi.find(seq[temp]);
                ThreoMass -= iter->second;
                sub = MonoMass - ThreoMass;
                k = sub - mod[0];
                temp++;
            }
        }
        if (temp == seq.size() - 1) {
            return -1;
        }
        cut_laction = temp;
        if (cut_laction > 30) {
            return -1;
        }
    }
    return j;
}

int mylib::speculate_lib::spect_ETD(double &c, double key, vector<double> &mod, char a[], int &cut_laction_C_N,
                                    int &cut_laction_N_C)             //推测末端截断以及所带修饰质量
{
    std::map<char, double>::iterator iter = mapi.begin();
    int temp = 0;
    int temp1 = 0;
    int last = strlen(a) - 1;
    double k;
    int i, j = 0, ac = 0;
    if (fabs(c - key) < 1) {
        cut_laction_C_N = 0;
        cut_laction_N_C = 0;
        return 0;
    } else {
        iter = mapi.find(a[temp]);
        c -= iter->second;
        temp++;
        // cout<<"max = "<<c<<endl;
        double sub = key - c + 1.9919 - 18.01056;  //记录差值
        // cout<<"sub = "<<sub<<endl;
        k = sub - mod[0];
        // cout<<"k = "<<k <<endl;
        while (ac != 1 && (last > strlen(a) - 6)) {
            for (i = 1; i < mod.size(); i++) {
                if (sub - mod[i] > 0 && mod[i] != 0) {
                    if (k > sub - mod[i]) {
                        k = sub - mod[i];
                        j = i;
                    }
                    if (fabs(k) < 1) {
                        ac = 1;
                    }
                } else break;
            }
            if (ac != 1) {
                iter = mapi.find(a[last]);
                c -= iter->second;
                sub = key - c + 1.9919 - 18.01056;
                // cout<<"fdsfsfsf = "<<sub <<endl;
                k = sub - mod[0];
                last--;
                temp1++;
            }
        }
        cut_laction_C_N = temp;
        cut_laction_N_C = temp1;
        // cout<<"c_n = "<<cut_laction_C_N<< "  n_C = "<<cut_laction_N_C<<endl;
        if (ac == 0) {
            return -1;
        }
    }
    return j;
}

int mylib::speculate_lib::spect_ETD_s(double &c, double key, vector<double> &mod,
                                      const string &a, int &cut_laction_C_N, int &cut_laction_N_C)             //推测末端截断以及所带修饰质量
{
    std::map<char, double>::iterator iter = mapi.begin();
    int temp = 0;
    int temp1 = 0;
    int last = a.length() - 1;
    double k;
    int i, j = 0, ac = 0;
    if (fabs(c - key) < 1) {
        cut_laction_C_N = 0;
        cut_laction_N_C = 0;
        return 0;
    } else {
        while (temp != last && ac != 1) {
            iter = mapi.find(a[temp]);
            c -= iter->second;
            temp++;
            // cout<<"max = "<<c<<endl;
            double sub = key - c + 1.9919 - 18.01056;  //记录差值
            // cout<<"sub = "<<sub<<endl;
            k = sub - mod[0];
            // cout<<"k = "<<k <<endl;
            while (ac != 1 && (last > a.length() - 6)) {
                for (i = 1; i < mod.size(); i++) {
                    if (sub - mod[i] > 0 && mod[i] != 0) {
                        if (k > sub - mod[i]) {
                            k = sub - mod[i];
                            j = i;
                        }
                        if (fabs(k) < 1) {
                            ac = 1;
                        }
                    } else break;
                }
                if (ac != 1) {
                    iter = mapi.find(a[last]);
                    c -= iter->second;
                    sub = key - c + 1.9919 - 18.01056;
                    // cout<<"fdsfsfsf = "<<sub <<endl;
                    k = sub - mod[0];
                    last--;
                    temp1++;
                }
            }
        }
        cut_laction_C_N = temp;
        cut_laction_N_C = temp1;
        if (ac == 0) {
            return -1;
        }
    }
    return j;
}

int mylib::speculate_lib::Set_Data(std::vector<double> &data, std::vector<double> &data_y,
                                   char reference[], double Pr_mass, vector<double> &mod, double &modifier_mass, double &pr_sub,
                                   int &cut_locat) {
    int n = strlen(reference);
    int i = 0;
    std::map<char, double>::iterator iter = mapi.begin();
    //制造y离子理论谱
    int temp = 0;
    // cout<<reference[0]<<endl;
    // cout<<"make"<< "mapi.size = "<<mapi.size()<<" n = "<<n<<endl;
    for (i = n - 1, iter = mapi.begin(); i >= 0; i--) {
        iter = mapi.find(reference[i]);
        if (i == n - 1) {
            data_y.push_back(iter->second + 19.017836 - 1.007276);//(18.01056)
            temp++;
            // cout<<"temp = "<<temp<<endl;
        } else if (iter != mapi.end()) {
            data_y.push_back(iter->second + data_y[temp - 1]);
            temp++;
        } else if (iter == mapi.end()) {
            data_y.push_back(data_y[temp - 1]);
            temp++;
        }
    }
    // cout<<"data_y.size() = "<<data_y.size()<<endl;
    double max = data_y[data_y.size() - 1];   //记录整条肽段的质量
    cout << "max = " << max << endl;
    if (fabs(max - Pr_mass) > 500) {
        // cout<<"no ! mass error "<<endl;
        return 0;
    }
    int cut_location = 0;                  //记录截断个数
    int modif_location = mylib::speculate_lib::spect(max, Pr_mass, mod, reference, cut_location);

    if (modif_location == -1) {
        //cout<<"no ! error "<<endl;
        return 0;
    }
    modifier_mass = mod[modif_location];        //修饰质量
    //----------------验证一下----------------//
    cout << "Modification   mass = " << setiosflags(ios::fixed) << setprecision(5) << mod[modif_location] << endl;
    cout << "Throe mass Pre mass = " << setiosflags(ios::fixed) << setprecision(5) << mod[modif_location] + max << endl;
    cout << " Mono mass Pre mass = " << setiosflags(ios::fixed) << setprecision(5) << Pr_mass << endl;
    cout << "末端截断的个数 = " << cut_location << endl;
    pr_sub = mod[modif_location] + max - Pr_mass;
    i = cut_location;
    cut_locat = cut_location;
    while (i > 0)          //截去末端截断
    {
        data_y.pop_back();
        i--;
    }
    n = data_y.size();
    //-------------------制造b离子理论谱----------------------//(测试)
    for (i = n - 2; i >= 0; i--) {
        data.push_back(max - data_y[i]);
    }
    data.push_back(data_y[data_y.size() - 1]);
    n = data.size();
    return n;
}

int mylib::speculate_lib::Set_Data_s(std::vector<double> &data, std::vector<double> &data_y,
                                     const string &seq, double &Pr_mass, vector<double> &mod,
                                     double &modifier_mass, double &pr_sub, int &cut_locat) {
    int n = seq.length();
    // cout<<seq.size()<<" lenth = " <<seq.length()<<endl;
    int i = 0;
    std::map<char, double>::iterator iter = mapi.begin();
    //制造y离子理论谱
    int temp = 0;
    // cout<<reference[0]<<endl;
    // cout<<"make"<< "mapi.size = "<<mapi.size()<<" n = "<<n<<endl;
    for (i = n - 1, iter = mapi.begin(); i >= 0; i--) {
        iter = mapi.find(seq[i]);
        if (i == n - 1) {
            data_y.push_back(iter->second + 19.017836 - 1.007276);//(18.01056)
            temp++;
            // cout<<"temp = "<<temp<<endl;
        } else if (iter != mapi.end() && data_y.size() > 0) {
            data_y.push_back(iter->second + data_y[temp - 1]);
            temp++;
        } else if (iter == mapi.end() && data_y.size() > 0) {
            data_y.push_back(data_y[temp - 1]);
            temp++;
        }
    }
    // cout<<"data_y.size() = "<<data_y.size()<<endl;
    double max = data_y[data_y.size() - 1];   //记录整条肽段的质量
    // cout<<"y 理论值最大值 = "<<max<<"质谱中该肽段的质量 = "<<Pr_mass<<" 相差的质量 = " <<Pr_mass - max<<endl;

    if (fabs(max - Pr_mass) > 500) {
        // cout<<"no ! mass error "<<endl;
        return 0;
    }
    int cut_location = 0;                  //记录截断个数
    //modif_location接受到的是一个修饰值所在数组的序号
    //SPECT_s返回-1表示推断失败。返回0，应该表示无修饰的情况。其他数则为修饰总质量所在数组的下标。
    int modif_location = mylib::speculate_lib::spect_s(max, Pr_mass, mod, seq, cut_location);
    // cout<<"cutLocation = "<<modif_location<<endl;
    if (modif_location == -1) {
        //cout<<"no ! error "<<endl;
        return 0;
    }
    modifier_mass = mod[modif_location];        //修饰质量
    //----------------验证一下----------------//
    // cout<<"Modification   mass = "<<setiosflags(ios::fixed)<<setprecision(5)<< mod[modif_location]<<endl;
    // cout<<"Throe mass Pre mass = "<<setiosflags(ios::fixed)<<setprecision(5)<< mod[modif_location] + max<<endl;
    // cout<<" Mono mass Pre mass = " <<setiosflags(ios::fixed)<<setprecision(5)<<Pr_mass<<endl;
    // cout<<"末端截断的个数 = " << cut_location<<endl;
    pr_sub = mod[modif_location] + max - Pr_mass;
    i = cut_location;
    cut_locat = cut_location;
    while (i > 0)          //截去末端截断
    {
        data_y.pop_back();
        i--;
    }
    // cout<<"截断后data_y.size() = "<<data_y.size()<<endl;
    n = data_y.size();
    //-------------------制造b离子理论谱----------------------//(测试)
    for (i = n - 2; i >= 0; i--) {
        data.push_back(max - data_y[i]);
    }
    data.push_back(data_y[data_y.size() - 1] - 18.01056);
    n = data.size();
    // cout<<"b 理论峰的个数为: "<<data.size() <<" y 理论峰的个数为 : "<<data_y.size()<<endl;
    // cout<<"b max = "<<data.back()<<" y max = "<<data_y.back()<<endl;
    return n;
}

int mylib::speculate_lib::Set_Data_ETD(std::vector<double> &data, std::vector<double> &data_y,
                                       char reference[], double Pr_mass, vector<double> &mod, double &modifier_mass, double &pr_sub,
                                       int &cut_locat_c_n, int &cut_locat_n_c) {
    int n = strlen(reference);
    // cout<<"序列总长度为 : "<<n<<endl;
    int i = 0;
    std::map<char, double>::iterator iter;
    //制造Z离子理论谱
    int temp = 0;
    for (i = n - 1, iter = mapi.begin(); i >= 0; i--) {
        iter = mapi.find(reference[i]);
        if (i == n - 1) {
            data_y.push_back(iter->second + 1.9919);//(18.01056)
            temp++;
        } else if (iter != mapi.end()) {
            data_y.push_back(iter->second + data_y[temp - 1]);
            temp++;
        } else if (iter == mapi.end()) {
            data_y.push_back(data_y[temp - 1]);
            temp++;
        }
    }
    //制造z离子理论谱结束
    double max = data_y[data_y.size() - 1];   //记录整条肽段的质量
    if (fabs(max - Pr_mass) > 500) {
        // cout<<"no ! mass error "<<endl;
        return 0;
    }
    int cut_location_C_N = 0;
    int cut_location_N_C = 0;                  //记录截断个数
    int modif_location = mylib::speculate_lib::spect_ETD(max, Pr_mass, mod, reference, cut_location_C_N, cut_location_N_C);

    if (modif_location == -1) {
        // cout<<"no ! error "<<endl;
        return 0;
    }
    modifier_mass = mod[modif_location];        //修饰质量
    max += 18.01056 + 1.007276 - 1.9919;
    //----------------验证一下----------------//
    // cout<<"Modification   mass = "<<setiosflags(ios::fixed)<<setprecision(5)<< mod[modif_location]<<endl;
    // cout<<"Throe mass Pre mass = "<<setiosflags(ios::fixed)<<setprecision(5)<< mod[modif_location] + max <<endl;
    // cout<<"Mono  mass Pre mass = " <<setiosflags(ios::fixed)<<setprecision(5)<<Pr_mass<<endl;
    // cout<<"C 端截断的个数 = " << cut_location_C_N<<" , N 端截断的个数 = " << cut_location_N_C<<endl;
    pr_sub = mod[modif_location] + max - Pr_mass;
    cut_locat_n_c = cut_location_N_C;
    cut_locat_c_n = cut_location_C_N;
    i = cut_locat_c_n;
    while (i > 0)          //截去C端截断
    {
        data_y.pop_back();
        i--;
    }
    if (cut_location_N_C != 0)    //截去N端截断,此时，之后的每一个点都要截去最后一段肽段的质量
    {
        iter = mapi.begin();
        n = strlen(reference) - 1;
        double cut_mass = 0;
        vector<double>::iterator it = data_y.begin();
        i = 0;
        for (it = data_y.begin(); it != data_y.end();) {
            if (i < cut_location_N_C) {
                it = data_y.erase(it);
                iter = mapi.find(reference[n]);
                cut_mass += iter->second;
                // cout<<"cut_mass = "<<cut_mass<<endl;
                i++;
                n--;
            } else {
                (*it) = (*it) - cut_mass;
                it++;
            }
        }
    }
    n = data_y.size();
    // cout<<"data_y.size() = "<<data_y.size()<<endl;
    // for(i = 0 ; i < data_y.size() ; i++)
    // {
    //     cout<<setiosflags(ios::fixed)<<setprecision(5)<<data_y[i]<<endl;
    // }
    //-------------------制造c离子理论谱----------------------//(测试)
    for (i = n - 2; i >= 0; i--) {
        data.push_back(max - data_y[i]);
    }
    // data.push_back(data_y[data_y.size()-1]);
    n = data.size();
    return n;
}

bool spectSortFun(c::arg::cutLocationPtr a, c::arg::cutLocationPtr b) {
    return fabs(a->variablePtmsMasss) < fabs(b->variablePtmsMasss);
}

int mylib::speculate_lib::Set_Data_ETD_s(std::vector<double> &data, std::vector<double> &data_y,
                                         const string &seq, double Pr_mass, map<double, string> &modifyTable,
                                         vector<c::arg::cutLocationPtr> &cutPtrVec) {
    int n = seq.length();
    // cout<<"序列总长度为 : "<<n<<endl;
    int i = 0;
    std::map<char, double>::iterator iter;
    //制造Z离子理论谱
    int temp = 0;
    for (i = n - 1, iter = mapi.begin(); i >= 0; i--) {
        iter = mapi.find(seq[i]);
        if (i == n - 1) {
            data_y.push_back(iter->second + 1.9919);//(18.01056)
            temp++;
        } else if (iter != mapi.end()) {
            data_y.push_back(iter->second + data_y[temp - 1]);
            temp++;
        } else if (iter == mapi.end()) {
            data_y.push_back(data_y[temp - 1]);
            temp++;
        }
    }
    //制造z离子理论谱结束
    //----------------制造c离子
    for (i = data_y.size() - 2; i >= 0; i--) {
        data.push_back(data_y.back() - data_y[i]);
    }
    data.push_back(data_y.back() - 1.9919);
    //----------------制造c离子结束
    // d::deal::SpectReWrite(data_y,mod,modifyTable,Pr_mass,cutPtrVec);
    // cout<<"可能的形式为:"<<cutPtrVec.size()<<endl;
    mylib::speculate_lib::spectProteinForm(data_y, modifyTable, Pr_mass, cutPtrVec);
    sort(cutPtrVec.begin(), cutPtrVec.end(), spectSortFun);
    return 0;
}

//函数返回INT_MAX或者修饰表下标的值，表示推测修饰成功，或失败。
int mylib::speculate_lib::SpectReWrite(vector<double> &TheroMass, vector<double> &mod, const map<double, string> &modifyTable,
                                       const double &MonoMass, vector<c::arg::cutLocationPtr> &cutPtrVec) {
    // vector<c::arg::cutLocationPtr> cutPtrVec;   //存放截断位置的指针
    vector<int> temp(TheroMass.size(), -1);      //存放截断
    for (int i = TheroMass.size() - 1; i >= 0; i--) {     //------------------------遍历截断向量T
        if (TheroMass.size() - 1 - i > 30) { //--------------------
            break;
        }
        if (fabs(TheroMass[i] - MonoMass) > 500) {      //--------------------如果相差过大，则当前截断不成立
            continue;
        }
        for (map<double, string>::const_iterator j = modifyTable.begin();
             j != modifyTable.end(); j++) {          //--------------------对于向量T每个元素Ti都需要查询修饰表的质量
            double tMass = TheroMass[i] + j->first - MonoMass;
            double sMass = MonoMass + j->first - TheroMass[i];
            if (fabs(tMass) < 2 || fabs(sMass) < 2) {    //如果理论修饰得到一定的Ptms偏移在precurmass的一定范围内
                if (fabs(tMass) / (TheroMass[i] + j->first) * 1000000 <= 15 ||
                    fabs(sMass) / (MonoMass + j->first) * 1000000 <= 15) {
                    c::arg::cutLocationPtr cutPtr = std::make_shared<c::arg::cutLocation>();   //----------------创建可能的截断点
                    if (j->first <= 0) {
                        cutPtr->cutLocationN = TheroMass.size() - 1 - i;   //--------------------保存截断位置
                        cutPtr->variablePtmsMasss = j->first;         //--------------------保存可变修饰质量偏移
                        cutPtr->ptmsSeq = j->second;                //------------------------保存当前可能的修饰
                    } else if (j->first > 0 && fabs(tMass) <
                                               2) {     //------------------------如果修饰为正，那么一定是理论修饰+可变PTMs质量偏移才能接近PRESCORE_MASS
                        cutPtr->cutLocationN = TheroMass.size() - 1 - i;
                        cutPtr->variablePtmsMasss = j->first;
                        cutPtr->ptmsSeq = j->second;
                    }
                    cutPtr->ppm = (tMass) / (TheroMass[i] + j->first) * 1000000;         //----------------保存ppm
                    cutPtr->theroMass = TheroMass[i];       //-------------------------------保存当前截断理论质量
                    cutPtr->precursorMass = MonoMass;       //----------------------------------保存PrecursorMass
                    cutPtrVec.push_back(cutPtr);
                }
            }
        }
    }
    // for(int i = 0 ; i < cutPtrVec.size() ;i++) {            //--------------------输出推测结果。测试使用
    //     cout<<"variblePtmsMass = "<<cutPtrVec[i]->variablePtmsMasss<<" Ptms = "<<cutPtrVec[i]->ptmsSeq<<endl;;
    //     cout<<"ProteinformMass = "<<cutPtrVec[i]->theroMass+cutPtrVec[i]->variablePtmsMasss<<endl;
    //     cout<<"PrecursorMass = "<<cutPtrVec[i]->precursorMass<<endl;
    //     cout<<"PPM = "<<cutPtrVec[i]->ppm<<endl;
    //     cout<<"cutSize = "<<cutPtrVec[i]->cutLocationN<<endl;
    //     cout<<endl;
    // }
    return 0;
}


int mylib::speculate_lib::SetDataHCD(std::vector<double> &data, std::vector<double> &data_y,
                                     const string &seq, double &Pr_mass,
                                     vector<double> &mod, const map<double, string> &modifyTable,
                                     vector<c::arg::cutLocationPtr> &cutPtrVec, double addMass) {
    int n = seq.length();
    int i = 0;
    std::map<char, double>::iterator iter = mapi.begin();
    // d::deal::SpectReWrite(data_y,mod,modifyTable,Pr_mass,cutPtrVec);
    // cout<<"可能的形式为:"<<cutPtrVec.size()<<endl;
    if (cutPtrVec.empty()) {
        // d::deal::spectProteinFormReWrite(data_y ,modifyTable,Pr_mass ,cutPtrVec,addMass,                                                                                                                                                                                                                                             seq);
    }
    sort(cutPtrVec.begin(), cutPtrVec.end(), spectSortFun);
    // cout<<"可能的形式为:"<<cutPtrVec.size()<<endl;
    return n;
}

void mylib::speculate_lib::spectProteinForm(vector<double> &TheroMass, const map<double, string> &modifyTable,
                                            const double &precursorMass, vector<c::arg::cutLocationPtr> &cutPtrVec) {
    if (TheroMass.empty() || modifyTable.empty()) {    //--------------------条件判断
        return;
    }
    //----------------产生范围
    vector<double> TheroMass1 = TheroMass;
    // double minPrecursorMass = precursorMass - 500 ; 
    // double maxPrecursorMass = precursorMass + 500 ; 
    //----------------产生截断向量序列
    for (int ji = 0; ji < TheroMass1.size(); ++ji) {   //--------------------控制前截
        for (int xi = TheroMass1.size() - 1; xi >= ji; --xi) {      //--------------------控制后截
            if (fabs(TheroMass1[xi] - precursorMass) > 500) {       //--------------------如果超过500，遍历下一值
                continue;
            } else {    //----------------对于该值进行修饰探测
                for (map<double, string>::const_iterator it = modifyTable.begin(); it != modifyTable.end(); ++it) {
                    if (fabs(TheroMass1[xi] + it->first - precursorMass) > 1) {
                        continue;
                    }
                    double ppm = (TheroMass1[xi] + it->first - precursorMass) / (TheroMass1[xi] + it->first) * 1000000;
                    if (fabs(ppm) < 15) {
                        c::arg::cutLocationPtr cp = std::make_shared<c::arg::cutLocation>();
                        cp->cutLocationC = ji;
                        cp->cutLocationN = TheroMass1.size() - 1 - xi;
                        cp->ppm = ppm;
                        cp->precursorMass = precursorMass;
                        cp->ptmsSeq = it->second;
                        cp->theroMass = TheroMass1[xi];
                        cp->variablePtmsMasss = it->first;
                        cutPtrVec.push_back(cp);
                    }
                }
            }
        }//----------------后端截断结束
        for (int ui = ji + 1; ui < TheroMass1.size(); ++ui) {   //--------------------截去前端
            TheroMass1[ui] -= TheroMass1[ji];
        }
        if (TheroMass1.back() + 500 - precursorMass < 0) {     //-----如果截取到理论最大值与质谱前体质量相差500,推断结束
            break;
        }
    }//--------------------前段截断结束
    // cout<<"可能的蛋白质形式个数为 : " <<cutPtrVec.size() <<endl;
    // for(int xi = 0 ; xi < cutPtrVec.size() ; ++xi ) {
    //     cout<<"n端截断为 : "<<cutPtrVec[xi]->cutLocationN<<" c端截断为 : " <<cutPtrVec[xi]->cutLocationC<<endl;
    //     cout<<"ppm = "<<cutPtrVec[xi]->ppm<<endl;
    //     cout<<"Ptms Mass = "<<cutPtrVec[xi]->variablePtmsMasss<<endl;
    //     cout<<"ProteinformMass = "<<cutPtrVec[xi]->theroMass + cutPtrVec[xi]->variablePtmsMasss <<endl;
    //     cout<<"PrecursorMass = "<<cutPtrVec[xi]->precursorMass<<endl;
    // }
}


void findProteinSeq(const vector<double> &theroMass, vector<double> &ms, int cutLocStart, int cutLocEnd, double addMass,
                    int &score) {
    if (theroMass.empty()) {
        return;
    }
    int length = theroMass.size();
    // cout<<"theroMasssize = "<<theroMass.size()<<"  cutsizestart = "<<cutLocStart<<" cutend = "<<cutLocEnd<<" cutSum = "<<cutLocStart+cutLocEnd<<endl;
    vector<double> theroMassT;
    vector<double>::const_iterator itStart = theroMass.begin() + cutLocStart;
    vector<double>::const_iterator itEnd = theroMass.end() - cutLocEnd;
    theroMassT.assign(itStart, itEnd);
    if (cutLocStart >= 1) {
        double subMass = theroMass[cutLocStart - 1] - addMass;    //----------------前截质量
        for (int i = 0; i <= theroMassT.size() / 5; i++) {        //----------------更新质量
            theroMassT[i] -= subMass;
        }
    }
    for (int xi = 0; xi < theroMassT.size(); xi++) {
        if (xi > theroMassT.size() / 5) {
            break;
        }
        int n = mylib::speculate_lib::txpb_binarey_search_ex(ms, ms.size() / 2, theroMassT[xi]);
        if (n < 0) {
            continue;
        } else if (fabs((theroMassT[xi] - ms[n]) / theroMassT[xi] * 1000000) < 15) {
            // cout<<"Thero = "<<theroMassT[xi]<<" ms = "<<ms[n]<<endl;
            score++;
        }
    }
}

void mylib::speculate_lib::cutMassSeq(vector<double> &theroMass, int cutLocStart, int cutLocEnd, double addMass) {
    if (theroMass.empty()) {
        return;
    }
    vector<double> theroMassT;
    if (cutLocStart > theroMass.size() - cutLocEnd) {
        return;
    }
    vector<double>::const_iterator itStart = theroMass.begin() + cutLocStart;
    vector<double>::const_iterator itEnd = theroMass.end() - cutLocEnd;
    theroMassT.assign(itStart, itEnd);
    if (cutLocStart >= 1) {
        double subMass = theroMassT[cutLocStart - 1] - addMass;
        for (int i = 0; i < theroMassT.size(); i++) {
            theroMassT[i] -= subMass;
        }
    }
    theroMass = theroMassT;
}


void mylib::speculate_lib::spectProteinFormReWrite(const vector<double> &TheroMass, const vector<double> &TheroMassN,
                                                   vector<double> &ms,
                                                   const map<double, string> &modifyTable, vector<double> &ptmMass,
                                                   const double &precursorMass, vector<c::arg::cutLocationPtr> &cutPtrVec,
                                                   double addMass, const string &prtSeq) {
    if (TheroMass.empty() || modifyTable.empty()) {    //--------------------条件判断
        return;
    }
    //----------------产生范围
    vector<double> TheroMass1 = TheroMass;
    double minPrecursorMass = precursorMass - 500;
    double maxPrecursorMass = precursorMass + 500;
    int leftRange = 0;
    int rightRange = TheroMass1.size();
    int seqLog = prtSeq.length() - 1;
    double subMass = 0;
    //----------------产生截断向量序列
    for (int ji = 0; ji < TheroMass1.size(); ++ji) {   //--------------------控制前截
        if (TheroMass1.back() < minPrecursorMass) { //----------------如果理论最大值  < precursormass - 500  。无论如何添加PTMS质量都无法到达
            break;
        }
        for (int xi = leftRange + 1; xi < TheroMass1.size(); ++xi) {      //--------------------控制后截
            if (TheroMass1[xi] < minPrecursorMass) {    //----------------如果小于precursormass - 500 ，continue
                leftRange = xi;
                continue;
            }
            if (TheroMass1[xi] > maxPrecursorMass) { //----------------如果大于precursormass + 500 , break ;
                rightRange = xi;
                break;
            }
            int n = mylib::speculate_lib::txpb_binarey_search_ex(ptmMass, ptmMass.size(), precursorMass - TheroMass1[xi]);
            // cout<<"TheroMass1 = "<<TheroMass1[xi]<<"Ther+modmass = "<<TheroMass1[xi]+ptmMass[n]<<" cutEnd = "<<TheroMass1.size() - 1 - xi<<endl;
            if (n < 0) {
                continue;
            }
            if (fabs(TheroMass1[xi] + ptmMass[n] - precursorMass) > 1) {  //----------------控制precuremass与蛋白质形式质量 < 1
                continue;
            }
            double ppm = (TheroMass1[xi] + ptmMass[n] - precursorMass) / (TheroMass1[xi] + ptmMass[n]) * 1000000;
            if (fabs(ppm) < 15) {   //--------------------如果ppm < 15 ，判断是否有完全相等的离子
                // int score = 0 ;   
                //     findProteinSeq(TheroMass,ms,ji,TheroMass1.size() - 1 - xi,addMass,score);
                // if(score == 0 ) {
                //     findProteinSeq(TheroMassN,ms,TheroMass1.size()-1-xi,ji,addMass,score);
                // }
                // // cout<<"cutStart = "<<ji<<" cutEnd = "<<TheroMass1.size() - 1 - xi<<"score = "<<score<<endl;
                // if( score >= 1 ) {
                c::arg::cutLocationPtr cp = std::make_shared<c::arg::cutLocation>();
                cp->cutLocationC = ji;
                cp->cutLocationN = TheroMass1.size() - 1 - xi;
                cp->ppm = ppm;
                // cp->precursorMass = precursorMass ; 
                // cp->ptmsSeq = it->second ; 
                // cp->theroMass = TheroMass1[xi]; 
                cp->variablePtmsMasss = ptmMass[n];
                cutPtrVec.push_back(cp);
                // }
            }
        }//----------------后端截断结束
        if (mapi.count(prtSeq[seqLog])) {
            subMass = mapi[prtSeq[seqLog--]];
        } else {
            seqLog--;
        }
        for (int ui = leftRange; ui < TheroMass1.size(); ++ui) {   //--------------------截去前端
            TheroMass1[ui] -= subMass;
        }
    }//--------------------前段截断结束
    sort(cutPtrVec.begin(), cutPtrVec.end(), spectSortFun);
    // cout<<"可能的蛋白质形式个数为 : " <<cutPtrVec.size() <<endl;
    // for(int xi = 0 ; xi < cutPtrVec.size() ; ++xi ) {
    //     cout<<"n端截断为 : "<<cutPtrVec[xi]->cutLocationN<<" c端截断为 : " <<cutPtrVec[xi]->cutLocationC<<endl;
    //     cout<<"ppm = "<<cutPtrVec[xi]->ppm<<endl;
    //     cout<<"Ptms Mass = "<<cutPtrVec[xi]->variablePtmsMasss<<endl;
    //     cout<<"ProteinformMass = "<<cutPtrVec[xi]->theroMass + cutPtrVec[xi]->variablePtmsMasss <<endl;
    //     cout<<"PrecursorMass = "<<cutPtrVec[xi]->precursorMass<<endl;
    //     cout<<"#END\n\n";
    // }
}

int mylib::speculate_lib::txpb_binarey_search_ex(const vector<double> &pstArray, int iLength, double ullTime) {
//    int middle = 0;
//    int low = 0, high = iLength - 1;
//    if (pstArray.empty() || iLength < 1) {
//        printf("%s %d error !\n", __func__, __LINE__);
//        return -1;
//    }
//
//    if (iLength == 1) {
//        return iLength - 1;
//    }
//
//    if (ullTime <= pstArray[0]) {
//        return 0;
//    } else if (ullTime >= pstArray[iLength - 1]) {
//        return iLength - 1;
//    }
//
//    while (low <= high) {
//
//
//        middle = (low + high) / 2;
//        if (abs(pstArray[middle + 1] - ullTime) > abs(pstArray[middle] - ullTime)) {
//            high = middle - 1;
//        } else {
//            low = middle + 1;
//        }
//    }
//    return (abs(pstArray[middle + 1] - ullTime) > abs(pstArray[middle] - ullTime)) ? middle : (middle + 1);
//
//
//
    ///by lyc


    if (pstArray.empty() || iLength < 1) {
        printf("%s %d error !\n", __func__, __LINE__);
        return -1;
    }

    if (iLength == 1) {
        return iLength - 1;
    }

    if (ullTime <= pstArray[0]) {
        return 0;
    } else if (ullTime >= pstArray[iLength - 1]) {
        return iLength - 1;
    }

    int left = 0;
    int right = iLength - 1;

    if(ullTime <= pstArray[left])
        return left;
    if(ullTime >= pstArray[right])
        return right;


    while(right - left > 1)        //while循环里面的r - l > 1 这个限定条件很重要，当r yu l 的差值为1 的时候就应该停止循环了，在分别ar[l] 与 ar[r] 两个值讨论看那个值更合适
    {
        int mid = (left + right) / 2;
        if(ullTime >= pstArray[mid])
            left = mid;	  //在这里 l 直接等于 mid ，否则会出错
        else
            right = mid;	  // 同理
    }
    long long int x = abs(pstArray[left] - ullTime);
    long long int y = abs(pstArray[right] - ullTime);
    return x <= y ? left : right;

}


void mylib::speculate_lib::spectProteinFormETD(const vector<double> &TheroMass, const vector<double> &TheroMassN,
                                               vector<double> &ms,
                                               const map<double, string> &modifyTable, vector<double> &ptmMass,
                                               const double &precursorMass, vector<c::arg::cutLocationPtr> &cutPtrVec,
                                               double addMass, const string &prtSeq) {
    if (TheroMass.empty() || modifyTable.empty()) {    //--------------------条件判断
        EX_TRACE("TheroMass or modifyTable is null ~\n");
        return;
    }
    //----------------产生范围
    double subMassETD = 18.01056 - 1.9919;
    vector<double> TheroMass1 = TheroMass;
    double minPrecursorMass = precursorMass - 500;
    double maxPrecursorMass = precursorMass + 500;
    int leftRange = 0;
    int rightRange = TheroMass1.size();
    int seqLog = prtSeq.length() - 1;
    double subMass = 0;
    //----------------产生截断向量序列
    for (int ji = 0; ji < TheroMass1.size(); ++ji) {   //--------------------控制前截
        if (TheroMass1.back() < minPrecursorMass) { //----------------如果理论最大值  < precursormass - 500  。无论如何添加PTMS质量都无法到达
            break;
        }
        for (int xi = leftRange + 1; xi < TheroMass1.size(); ++xi) {      //--------------------控制后截
            if (TheroMass1[xi] < minPrecursorMass) {    //----------------如果小于precursormass - 500 ，continue
                leftRange = xi;
                continue;
            }
            if (TheroMass1[xi] > maxPrecursorMass) { //----------------如果大于precursormass + 500 , break ;
                rightRange = xi;
                break;
            }
            int n = mylib::speculate_lib::txpb_binarey_search_ex(ptmMass, ptmMass.size(),
                                                    precursorMass - subMassETD - TheroMass1[xi]);
            // cout<<"TheroMass1 = "<<TheroMass1[xi] <<" Ther+modmass = "<<TheroMass1[xi]+ptmMass[n] <<" cutEnd = "<<TheroMass1.size() - 1 - xi
            // <<" ModMass = "<<ptmMass[n]<<endl;
            if (n < 0) {
                continue;
            }
            if (fabs(TheroMass1[xi] + subMassETD + ptmMass[n] - precursorMass) >
                1) {  //----------------控制precuremass与蛋白质形式质量 < 1
                continue;
            }
            // cout<<"TMass = "<<TheroMass1[xi] + subMassETD + ptmMass[n]<<endl;
            // cout<<"cutStart = "<<ji<<" cutEnd = "<<TheroMass1.size() - 1 - xi<<"mass = "<<TheroMass1[xi] + ptmMass[n] <<endl;
            double ppm = (TheroMass1[xi] + subMassETD + ptmMass[n] - precursorMass) /
                         (TheroMass1[xi] + subMassETD + ptmMass[n]) * 1000000;
            if (fabs(ppm) < 15) {   //--------------------如果ppm < 15 ，判断是否有完全相等的离子
                // cout<<"cutStart = "<<ji<<" cutEnd = "<<TheroMass1.size() - 1 - xi<<"mass = "<<TheroMass1[xi] + ptmMass[n] +subMassETD<<" precursorMass = "<<precursorMass<<endl;
                // int score = 0 ;   
                //     findProteinSeq(TheroMass,ms,ji,TheroMass1.size() - 1 - xi,addMass,score);
                // if(score <= 1 ) {
                //     findProteinSeq(TheroMassN,ms,TheroMass1.size()-1-xi,ji,addMass,score);
                // }
                // cout<<"cutStart = "<<ji<<" cutEnd = "<<TheroMass1.size() - 1 - xi<<"score = "<<score<<endl;
                // if( score > 1 ) {
                c::arg::cutLocationPtr cp = std::make_shared<c::arg::cutLocation>();
                cp->cutLocationC = ji;
                cp->cutLocationN = TheroMass1.size() - 1 - xi;
                cp->ppm = ppm;
                // cp->precursorMass = precursorMass ; 
                // cp->ptmsSeq = it->second ; 
                cp->theroMass = TheroMass1[xi];
                cp->variablePtmsMasss = ptmMass[n];
                cutPtrVec.push_back(cp);
                // }
            }
        }//----------------后端截断结束
        if (mapi.count(prtSeq[seqLog])) {
            subMass = mapi[prtSeq[seqLog--]];
        } else {
            seqLog--;
        }
        for (int ui = leftRange; ui < TheroMass1.size(); ++ui) {   //--------------------截去前端
            TheroMass1[ui] -= subMass;
        }
    }//--------------------前段截断结束
    sort(cutPtrVec.begin(), cutPtrVec.end(), spectSortFun);
    // cout<<"可能的蛋白质形式个数为 : " <<cutPtrVec.size() <<endl;
    // for(int xi = 0 ; xi < cutPtrVec.size() ; ++xi ) {
    //     cout<<"n端截断为 : "<<cutPtrVec[xi]->cutLocationN<<" c端截断为 : " <<cutPtrVec[xi]->cutLocationC<<endl;
    //     cout<<"ppm = "<<cutPtrVec[xi]->ppm<<endl;
    //     cout<<"Ptms Mass = "<<cutPtrVec[xi]->variablePtmsMasss<<endl;
    //     cout<<"ProteinformMass = "<<cutPtrVec[xi]->theroMass + cutPtrVec[xi]->variablePtmsMasss + subMassETD<<endl;
    //     cout<<"PrecursorMass = "<<precursorMass<<endl;
    //     cout<<"#END\n\n";
    // }
}  
