//
// Created by wenzhong on 2023/3/22.
//

#ifndef WTOP_PRSM_ALIGNMENT_FILTER_HPP
#define WTOP_PRSM_ALIGNMENT_FILTER_HPP

#include <vector>
#include <math.h>

#include "../mylib/data_steam.hpp"
#include "../mylib/etd_process.hpp"

using namespace std;

class prsm_alignment_filter {
public:
    int get_mass_sub(const vector<double> &reference, const vector<double> &reference_y, const vector<double> &peer,
                     vector<double> &sub, vector<double> &sub_y, const vector<std::pair<long, long>> &alignment,
                     const vector<std::pair<long, long>> &alignment_y);

    void index_ions(vector<double> &sub, vector<double> &index_y,
                    std::vector<double> &modifier,
                    double &modifier_mass);

    void index_ions_mass(vector<double> &sub, vector<double> &index, vector<double> &modifier, double &modifier_mass);

    void get_ppm_unknown(vector<double> &reference, vector<double> &peer, vector<std::pair<long, long>> &alignment,
                         vector<double> &index, vector<double> &modifier, vector<double> &ppm, double &modifier_mass);

    void get_ions_peaks(vector<node> &ions, vector<pair<long, long>> &alignment, vector<double> &index, vector<double> &ppm,
                   double minPPMValue);
};
typedef prsm_alignment_filter* prsm_alignment_filter_ptr ;

#endif //WTOP_PRSM_ALIGNMENT_FILTER_HPP
