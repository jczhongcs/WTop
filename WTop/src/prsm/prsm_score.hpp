//
// Created by wenzhong on 2023/3/22.
//

#ifndef WTOP_PRSM_SCORE_HPP
#define WTOP_PRSM_SCORE_HPP

#include "../ms/msalign.hpp"
#include "../protein/protein.hpp"
#include "../merge/modify.hpp"
#include "algorithm"

class prsm_score {
public:
    double analysis_mutx(const vector<node> &arr_ions, const vector<node> &arr_ions_2, const vector<double> &mono,
                         const vector<double> &arr_mod, double mod_mass);

    double get_complementIons_score(vector<node> &ionsC, vector<node> &ionsN, int length);

    double get_seqIons_score(vector<node> &ionsC, vector<node> &ionsN, int &count);

    int get_SeqIons_pair(vector<node> &t, vector<node> &m);

    void
    get_score(vector<node> &ions_r, vector<double> &refe_score, const vector<double> &mod_arr,
              vector<double> &peer_score,
              double &score, double base_score, double mod_mass);
    int get_seq_count(vector<node> &ions_r);


    double get_abundance_score(vector<node> &ions_n,const vector<double> &IonsMass,const map<double,double> &keyMsValueAbRadio);
    double get_avg_abundance_score(vector<node> &ions_n,vector<node> &ions_c,const vector<double> &IonsMass,const map<double,double> &keyMsValueAbRadio);
    double get_dyn_abu_score(vector<node> &ions_n,const vector<double> &IonsMass,const map<double,double> &keyMsValueAbRadio,double precursor_mass);
    double get_avg_dyn_abu_score(vector<node> &ions_n,vector<node> &ions_c,const vector<double> &IonsMass,const map<double,double> &keyMsValueAbRadio,double precursor_mass);
    int get_diffcharge_score(vector<node> &ions_n);
};

typedef prsm_score* prsm_score_ptr;


#endif //WTOP_PRSM_SCORE_HPP
