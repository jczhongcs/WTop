//
// Created by wenzhong on 2023/3/21.
//

#ifndef WTOP_PPM_VALUE_HPP
#define WTOP_PPM_VALUE_HPP

#include <math.h>
#include "memory"
class ppm_value {
public:

    ppm_value(double para_value) {
        min_ppm_value = para_value ;
    }
    bool calculate_ppm_value(double theory_value, double actual_value);
    bool calculate_ppm_with_H(double theory_value, double actual_value);

private:
    double min_ppm_value ;
};
typedef std::shared_ptr<ppm_value> ppm_value_ptr;

#endif //WTOP_PPM_VALUE_HPP
