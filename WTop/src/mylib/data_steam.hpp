#ifndef SETDATA_HPP
#define SETDATA_HPP

#include <iostream>
#include <vector>
#include <algorithm> // 算法头文件，提供迭代器
#include <fstream>
#include <iomanip>
#include <map> // STL
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <string>
#include <unistd.h>
#include <cmath>
#include <set>
#include <climits>

#include "../console/opts.h"
#include "../mylib/io.h"
#include "../mylib/fun.h"
#include "../mylib/cwdtw.h"

using namespace std;

// 贪婪策略取修饰值
struct node
{
    int thoe_id;  // 理论序号
    int mono_id;  // 质谱序号
    double index; // 差值
    double ppm;   // ppm值
};

struct modi
{ // 存储修饰信息
    double first;
    double second;
    double mod_mass;
};

struct locat_ions
{
    string n;
    string c;
};

namespace mylib
{

    namespace data_stream
    {
        int setdata(std::vector<double> &data,
                    std::vector<char> &reference);

        void rtrim_s(string &str);

        void rtrim(char *str);

        double cwdtw_map(vector<double> &reference, vector<double> &peer,
                         std::vector<std::pair<long, long>> &alignment,
                         vector<double> &ref_zscore, vector<double> &peer_zscore);

        int Get_sub(const vector<double> &reference, const vector<double> &reference_y,
                    const vector<double> &peer,
                    vector<double> &sub, vector<double> &sub_y,
                    const std::vector<std::pair<long, long>> &alignment,
                    const std::vector<std::pair<long, long>> &alignment_y);

        void WriteSequenceAlignment_peer_lons(const char *output, const std::vector<double> &peer_orig,
                                              std::vector<locat_ions> lons);

        void Read_massdata(const char a[], std::vector<double> &peer, double &Pr_mass);

        void MergeSort(std::vector<double> &arr, int l, int r);

        void Merge(std::vector<double> &arr, int l, int q, int r);

        int removeDuplicate(std::vector<double> &a, int l);

        void Scope_Align(const vector<double> &reference, const vector<double> &peer,
                         std::vector<std::pair<long, long>> &alignment,
                         std::vector<double> &modifier, const double &modifier_mass);

        void align_modification_down(vector<double> reference, vector<double> peer,
                                     std::vector<std::pair<long, long>> &alignment, std::vector<double> &modifier,
                                     double &modifier_mass);

        void align_modification_top(vector<double> reference, vector<double> peer,
                                    std::vector<std::pair<long, long>> &alignment, std::vector<double> &modifier,
                                    double &modifier_mass);

        int find_key(double key, vector<double> &mod);

        void index_ions(const vector<double> &sub, vector<double> &index_y,
                        std::vector<double> &modifier,
                        double &modifier_mass);

        int remove_false(double index, vector<double> &modifier, const double &modifier_mass);

        void align_top_down(vector<double> reference, vector<double> peer,
                            std::vector<std::pair<long, long>> &alignment, std::vector<double> &modifier,
                            double &modifier_mass);

        bool cmp1(pair<int, int> a, pair<int, int> b);

        void Get_ppm(const vector<double> &reference, const vector<double> &peer,
                     const std::vector<std::pair<long, long>> &alignment,
                     const std::vector<double> &index,
                     std::vector<double> &modifier, std::vector<double> &ppm,
                     double &modifier_mass);

        int find_value(double key, std::vector<double> &modifier, double &modifier_mass);

        void loca_ions(vector<node> &n_ions, vector<node> &c_ions, vector<locat_ions> &location,
                       vector<double> peer);

        void ions_location(vector<node> &IonsCN, vector<node> &IonsNC,
                           vector<node> &MergIonsCN, vector<node> &MergIonsNC,
                           const vector<pair<long, long>> &AlignmentCN, const vector<pair<long, long>> &AlignmentNC,
                           const vector<double> &ReferenceCN, const vector<double> &ReferenceNC,
                           const vector<double> &peer,
                           const vector<double> &IndexCN, const vector<double> &IndexNC,
                           const vector<double> &ppmCN, const vector<double> &ppmNC,
                           vector<modi> &Ptms, const vector<double> &modify_table,
                           int lenth, double ModifyMass);

        void Process_mod_ions(vector<node> &c_ions, vector<node> &merg_ions, vector<modi> &modf,
                              std::vector<std::pair<long, long>> &alignment, double modifier_mass);

        void Get_mod_ions_loca(vector<node> &c_ions, vector<node> &merg_ions, vector<modi> &mod,
                               const std::vector<std::pair<long, long>> &alignment,
                               const vector<double> &modify_table,
                               double modifier_mass);

        void Get_SCORE(vector<node> &ions, vector<double> &refe_score, const vector<double> &mod_arr,
                       vector<double> &peer_score, double &score, double base_score, double mod_mass);

        void Get_all_ions(vector<node> &c_ions, vector<node> &c_modf_ions,
                          const vector<pair<long, long>> &alignment_c,
                          const vector<double> &reference, const vector<double> &peer,
                          const vector<double> &index_c,
                          const vector<double> &ppm_c);

        void Last_ptm(vector<modi> &mod, double modifier);

        void Get_other_ions(const vector<pair<long, long>> &alignment, vector<modi> &mod,
                            vector<node> &ions,
                            const vector<double> &index, const vector<double> &ppm,
                            int lenth);

        void GetOtherIonsZero(const vector<pair<long, long>> &alignment, vector<node> &ions,
                              vector<double> &index, vector<double> &ppm);

        void IfDyHe(vector<node> &MergIonsNC, vector<modi> &Dyhe, double &ModifyMass);

        void AssumePtmInLast(vector<node> &MergIons, const vector<double> &Reference, vector<modi> &Ptms,
                             double ModifyMass);

        void PtmDyHeFind(vector<node> &MergIons, vector<modi> &Dyhe, double ModifyMass);

        void AnalysisCNIons(vector<node> &MergIonsCN, vector<node> &Ions, vector<modi> &Ptms, vector<double> &Reference,
                            double ModifyMass);

        void ModifyMassBiger0(vector<node> &IonsCN, vector<node> &IonsNC,
                              vector<node> &MergIonsCN, vector<node> &MergIonsNC,
                              const vector<pair<long, long>> &AlignmentCN,
                              const vector<pair<long, long>> &AlignmentNC,
                              const vector<double> &ReferenceCN, const vector<double> &ReferenceNC,
                              const vector<double> &peer,
                              const vector<double> &IndexCN, const vector<double> &IndexNC,
                              const vector<double> &ppmCN, const vector<double> &ppmNC,
                              vector<modi> &Ptms,
                              const vector<double> &modify_table,
                              int lenth, double ModifyMass);

        void ModifyMassSmaller0(vector<node> &IonsCN, vector<node> &IonsNC,
                                vector<node> &MergIonsCN, vector<node> &MergIonsNC,
                                const vector<pair<long, long>> &AlignmentCN,
                                const vector<pair<long, long>> &AlignmentNC,
                                const vector<double> &ReferenceCN, const vector<double> &ReferenceNC,
                                const vector<double> &peer,
                                const vector<double> &IndexCN, const vector<double> &IndexNC,
                                const vector<double> &ppmCN, const vector<double> &ppmNC,
                                vector<modi> &Ptms,
                                int lenth, double ModifyMass);

        void GetAnotherIons(const vector<pair<long, long>> &alignment, vector<modi> &mod,
                            const vector<node> &Mergions, vector<node> &Ions,
                            const vector<double> &index, const vector<double> &ppm,
                            int lenth);

        int init_theory_ions_mass_(vector<double> &theory_ions_mass,
                                   int cut_n, int cut_c, double addMass);

        bool verifyModify(string &raw, string &Ptms);

        bool searchMatchIons(vector<pair<long, long>> &alignmentN, vector<pair<long, long>> &alignmentC,
                             vector<double> &referenceN, vector<double> &referenceC, vector<double> &peer);

        void filterAlignment(vector<pair<long, long>> &alignment, vector<double> &sub,
                             vector<double> &index, vector<double> &ppm);

        int txpb_binarey_search_ex(const vector<double> &pstArray, int iLength, double ullTime);

        void scopeAlignemnt(const vector<double> &reference, const vector<double> &peer,
                            std::vector<std::pair<long, long>> &alignment,
                            std::vector<double> &modifier, const double &modifier_mass);

        long get_fragment_ions_nub(const std::vector<node> &all_ions);

        int get_modification_table(const string &mod_file_name, vector<string> &mod);

        int find_match_ions(const std::vector<double> &arr_mod, const std::vector<node> &arr_ions);

        int find_error_match_ions(const std::vector<double> &arr_mod, const std::vector<node> &arr_ions);

        bool if_is_error_ions(const std::vector<double> &arr_mod, const node &p);

        // 分析离子组成
        double analysis_all_ions(const std::vector<node> &arr_ions, const std::vector<node> &arr_ions_2,
                                 const std::vector<double> &mono, const std::vector<double> &arr_mod, double mod_mass);

        bool _ppm_H(double t, double m);
        bool _ppm(double t, double m);
        bool _ppm_H_re(double t, double m);

        void Scope_Align_ReWrite(const vector<double> &reference, const vector<double> &peer,
                                 std::vector<std::pair<long, long>> &alignment,
                                 std::vector<double> &modifier, const double &modifier_mass);

        void scopeAlignReWrite(const vector<double> &reference, const vector<double> &peer,
                               std::vector<std::pair<long, long>> &alignment,
                               std::vector<double> &modifier, const double &modifier_mass, double maxMass);

        double getHuBuIonsScore(vector<node> &ionsC, vector<node> &ionsN, int length);

        bool _ppm_H_getMinPPM(double t, double m, double &minPPM);

        double getSeqIonsScore(vector<node> &ionsC, vector<node> &ionsN, int &count);

        void Get_SCORE_Un(vector<node> &ions_r, vector<double> &refe_score, const vector<double> &mod_arr,
                          vector<double> &peer_score,
                          double &score, double base_score,
                          double mod_mass);

        double findClose(vector<double> arr, int n, double target);

    }

}

#endif