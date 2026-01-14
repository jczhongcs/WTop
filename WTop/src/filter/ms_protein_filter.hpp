//
// Created by wenzhong on 2023/3/21.
//

#ifndef WTOP_MS_PROTEIN_FILTER_HPP
#define WTOP_MS_PROTEIN_FILTER_HPP

#include <iostream>
#include <string>
#include <vector>

#include "filter_tag_match.hpp"
#include "filter_tag_string.hpp"

#include "../ms/complement_ions.hpp"
#include "../ms/msalign.hpp"
#include "../ms/msalign_tag.hpp"

#include "../protein/filter_protein.h"

#include "../mylib/speculate_lib.h"
#include "../protein/protein.hpp"
#include "../util/ppm_value.hpp"

class ms_protein_filter {
public:
    int findBestProteinByTag(const protein_ptr prIt, const msalign_ptr msIt,
                             const map<char, char> &ccMap, const set<string> &set,
                             const double filterMass);

    int zero_match_peaks_search(const vector<double> &mono_mass,
                                const vector<double> &theo_mass,
                                int &one_ptm_match_peaks);

    bool precursor_mass_scope_search(const vector<double> &mono_mass,
                                     double &mass_gap, double precursor_mass,
                                     double scope);

    int getFilterTagMatch(const protein_ptr prIt, const msalign_ptr msIt,
                          const map<char, char> &ccMap,
                          map<int, filter_tag_msg_ptr> &tagMap,
                          vector<filter_tag_match_ptr> &filterTagMatchNVec,
                          vector<filter_tag_match_ptr> &filterTagMatchCVec);

    string tagsMatchSequence(const protein_ptr &prIt, const msalign_ptr &msIt,
                             const complement_ions_ptr &compIonsPtr,
                             const msalign_tag_ptr &msalignTagPtr,
                             filter_protein_ptr &filterProPtr);

private:
    int minTagLength = 3;

    double minMassGap = 500;
};
typedef ms_protein_filter *ms_protein_filter_ptr;

#endif // WTOP_MS_PROTEIN_FILTER_HPP
