//
// Created by wenzhong on 2023/3/21.
//

#include "ppm_value.hpp"


bool ppm_value::calculate_ppm_value(double theory_value, double actual_value)
{
    ///原计算ppm，考虑多氢少氢
//    double minPPM = 0.0 ;
//    double mut_value = 1000000.0;
//    double h_value = 1.007276 ;
//    double PPM_re = (theory_value - actual_value)
//                    / (theory_value) * mut_value;
//    double PPM_DeH_re = (theory_value - h_value - actual_value)
//                        / (theory_value - h_value) * mut_value;
//    double PPM_AeH_re = (theory_value + h_value - actual_value)
//                        / (theory_value + h_value) * mut_value;
//    minPPM = fabs(PPM_re) > fabs(PPM_DeH_re) ? PPM_DeH_re : PPM_re;
//    minPPM = fabs(minPPM) > fabs(PPM_AeH_re) ? PPM_AeH_re : minPPM;
//    if (fabs(PPM_re) < min_ppm_value || fabs(PPM_DeH_re) < min_ppm_value || fabs(PPM_AeH_re) < min_ppm_value) {
//        return true;
//    }
//    return false;
    ///end

    ///计算ppm，by lyc，不计算多氢少氢
//
//    double mut_value = 1000000.0;
//
//    double PPM_re = (theory_value - actual_value)
//                    / (theory_value) * mut_value;
//
//    if (fabs(PPM_re) < min_ppm_value) {
//        return true;
//    }
//    return false;

    ///end test


    ///计算ppm，考虑同位素峰
    double mut_value = 1000000.0;
    double mass_gap = 1.0024 ;
    double PPM_re = (theory_value - actual_value)
                    / (theory_value) * mut_value;
    double PPM_DMass_re = (theory_value - mass_gap - actual_value)
                        / (theory_value - mass_gap) * mut_value;
    double PPM_AMass_re = (theory_value + mass_gap - actual_value)
                        / (theory_value + mass_gap) * mut_value;

    if (fabs(PPM_re) < min_ppm_value || fabs(PPM_DMass_re) < min_ppm_value || fabs(PPM_AMass_re) < min_ppm_value) {
        return true;
    }
    return false;


}

bool ppm_value::calculate_ppm_with_H(double theory_value, double actual_value)
{
    ///原计算ppm，考虑多氢少氢
//    double minPPM = 0.0 ;
//    double mut_value = 1000000.0;
//    double h_value = 1.007276 ;
//    double PPM_re = (theory_value - actual_value)
//                    / (theory_value) * mut_value;
//    double PPM_DeH_re = (theory_value - h_value - actual_value)
//                        / (theory_value - h_value) * mut_value;
//    double PPM_AeH_re = (theory_value + h_value - actual_value)
//                        / (theory_value + h_value) * mut_value;
//    minPPM = fabs(PPM_re) > fabs(PPM_DeH_re) ? PPM_DeH_re : PPM_re;
//    minPPM = fabs(minPPM) > fabs(PPM_AeH_re) ? PPM_AeH_re : minPPM;
//    if (fabs(PPM_re) < min_ppm_value || fabs(PPM_DeH_re) < min_ppm_value || fabs(PPM_AeH_re) < min_ppm_value) {
//        return true;
//    }
//    return false;
    ///end
    ///计算ppm，by lyc，不计算多氢少氢
//
//    double mut_value = 1000000.0;
//
//    double PPM_re = (theory_value - actual_value)
//                    / (theory_value) * mut_value;
//
//    if (fabs(PPM_re) < min_ppm_value) {
//        return true;
//    }
//    return false;

    ///end test

    ///计算ppm，考虑同位素峰
    double mut_value = 1000000.0;
    double mass_gap = 1.0024 ;
    double PPM_re = (theory_value - actual_value)
                    / (theory_value) * mut_value;
    double PPM_DMass_re = (theory_value - mass_gap - actual_value)
                          / (theory_value - mass_gap) * mut_value;
    double PPM_AMass_re = (theory_value + mass_gap - actual_value)
                          / (theory_value + mass_gap) * mut_value;

    if (fabs(PPM_re) < min_ppm_value || fabs(PPM_DMass_re) < min_ppm_value || fabs(PPM_AMass_re) < min_ppm_value) {
        return true;
    }
    return false;


}