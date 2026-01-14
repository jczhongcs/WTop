//
// Created by 碎雨粘霓裳 on 2022/9/4.
//

#include "calculate_lib.h"

bool mylib::calculate_lib::calculate_ppm_value_is_low(double tMass, double pMass, double hMass, double sPpm) {
    double ppm = (tMass - pMass) / (tMass) * 1000000;
    double ppm_DeH = (tMass - hMass - pMass) / (tMass - hMass) * 1000000;
    double ppm_AeH = (tMass + hMass - pMass) / (tMass + hMass) * 1000000;
    if (fabs(ppm) <= ppm || fabs(ppm_DeH) <= ppm || fabs(ppm_AeH) <= ppm) {
        return true ;
    }
    return false ;
}