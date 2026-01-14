#ifndef Garg_H_
#define Garg_H_

#include <iostream>
#include <vector>
#include <map>           // STL
#include <stdio.h>
#include <string>
#include <string.h>
#include <math.h>
#include <iomanip>
#include <utility>

#include "data_steam.hpp"
#include "terminal_truncation.h"

using namespace std;

extern std::map<char, double> mapi;

namespace mylib {
    namespace speculate_lib {
        int spect(double &c, double key, vector<double> &mod, char a[], int &cut_laction);             //推测末端截断以及所带修饰质量

        int spect_s(double &c, double key, vector<double> &mod, const string &a,
                    int &cut_laction);             //推测末端截断以及所带修饰质量

        int spect_ETD(double &c, double key, vector<double> &mod, char a[], int &cut_laction_C_N,
                      int &cut_laction_N_C);             //推测末端截断以及所带修饰质量

        int spect_ETD_s(double &c, double key, vector<double> &mod, const string &a, int &cut_laction_C_N,
                        int &cut_laction_N_C);             //推测末端截断以及所带修饰质量

        int Set_Data(std::vector<double> &data, std::vector<double> &data_y,
                     char reference[], double Pr_mass, vector<double> &mod, double &modifier_mass, double &pr_sub,
                     int &cut_locat);

        int Set_Data_s(std::vector<double> &data, std::vector<double> &data_y,
                       const string &seq, double &Pr_mass, vector<double> &mod, double &modifier_mass, double &pr_sub,
                       int &cut_locat);

        int Set_Data_ETD(std::vector<double> &data, std::vector<double> &data_y,
                         char reference[], double Pr_mass, vector<double> &mod, double &modifier_mass, double &pr_sub,
                         int &cut_locat_c_n, int &cut_locat_n_c);

        int Set_Data_ETD_s(std::vector<double> &data, std::vector<double> &data_y,
                           const string &seq, double Pr_mass, map<double, string> &modifyTable,
                           vector<c::arg::cutLocationPtr> &cutPtrVec);

        //HCD
        int SpectReWrite(vector<double> &TheroMass, vector<double> &mod, const map<double, string> &modifyTable,
                         const double &MonoMass, vector<c::arg::cutLocationPtr> &cutPtrVec);

        int SetDataHCD(std::vector<double> &data, std::vector<double> &data_y,
                       const string &seq, double &Pr_mass, vector<double> &mod, const map<double, string> &modifyTable,
                       vector<c::arg::cutLocationPtr> &cutPtrVec, double addMass);

        //-------推测蛋白质形
        void spectProteinForm(vector<double> &TheroMass, const map<double, string> &modifyTable,
                              const double &precursorMass, vector<c::arg::cutLocationPtr> &cutPtrVec);


        void
        spectProteinFormReWrite(const vector<double> &TheroMass, const vector<double> &TheroMassN, vector<double> &ms,
                                const map<double, string> &modifyTable, vector<double> &ptmMass,
                                const double &precursorMass, vector<c::arg::cutLocationPtr> &cutPtrVec, double addMass,
                                const string &prtSeq);

        int txpb_binarey_search_ex(const vector<double> &pstArray, int iLength, double ullTime);

        void cutMassSeq(vector<double> &theroMass, int cutLocStart, int cutLocEnd, double addMass);

        void spectProteinFormETD(const vector<double> &TheroMass, const vector<double> &TheroMassN,
                                 vector<double> &ms,
                                 const map<double, string> &modifyTable, vector<double> &ptmMass,
                                 const double &precursorMass, vector<c::arg::cutLocationPtr> &cutPtrVec, double addMass,
                                 const string &prtSeq);

    }
}
#endif