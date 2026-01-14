#ifndef ETD_PROCESS_HPP
#define ETD_PROCESS_HPP

#include <iostream>
#include <vector>
#include <math.h>

#include "data_steam.hpp"
#include "ions_path.h"

using namespace std;

namespace mylib {
    namespace etd_process {
        void ions_location_n_c(vector<node> &c_n_ions, vector<node> &n_c_ions, vector<node> &merg_ions,
                               vector<pair<long, long> > &alignment_n, vector<pair<long, long> > &alignment_c,
                               vector<double> &reference, vector<double> &peer,
                               vector<double> &index_n, vector<double> &index_c,
                               vector<double> &ppm_n, vector<double> &ppm_c, vector<modi> &mod,
                               int lenth, double modifier_mass);

        void ions_location_c_n(vector<node> &n_ions, vector<node> &c_ions, vector<node> &merg_ions,
                               vector<pair<long, long> > &alignment_n, vector<pair<long, long> > &alignment_c,
                               vector<double> &reference, vector<double> &peer,
                               vector<double> &index_n, vector<double> &index_c,
                               vector<double> &ppm_n, vector<double> &ppm_c, vector<modi> &mod,
                               int lenth, double modifier_mass);

        void Get_all_ions(vector<node> &modf_ions, vector<node> &modf_ions_other,
                          vector<node> &ions,
                          vector<pair<long, long> > &alignment,
                          vector<double> &reference, vector<double> &peer,
                          vector<double> &index,
                          vector<double> &ppm, int &lenth);

        void Get_mod_ions_loca(vector<node> &merg_ions, vector<node> &merg_ions_other,
                               vector<node> &ions, vector<modi> &modf,
                               std::vector<std::pair<long, long> > &alignment, double &modifier_mass, int &lenth);

        void loca_ions(vector<node> &n_ions, vector<node> &c_ions, vector<locat_ions> &location,
                       vector<double> peer);

        void ions_location_match_N(vector<node> &c_n_ions, vector<node> &n_c_ions,
                                   vector<node> &merg_ions_c_n, vector<node> &merg_ions_n_c,vector<modi> &mod,
                                   const vector<double> &modify_table,int lenth, double modifier_mass);

        void ions_location_match_C(vector<node> &c_n_ions, vector<node> &n_c_ions,
                                   vector<node> &merg_ions_c_n, vector<node> &merg_ions_n_c,vector<modi> &mod,
                                   const vector<double> &modify_table,
                                   int lenth, double modifier_mass);

        void Get_ppm(vector<double> &reference, vector<double> &peer, std::vector<std::pair<long, long> > &alignment,
                     std::vector<double> &index, std::vector<double> &modifier, std::vector<double> &ppm,
                     double &modifier_mass);

        void Get_ppm_n_c(vector<double> reference, vector<double> peer, std::vector<std::pair<long, long> > &alignment,
                         std::vector<double> &index, std::vector<double> &modifier, std::vector<double> &ppm,
                         double &modifier_mass);

        void Get_ppm_c_n(vector<double> reference, vector<double> peer, std::vector<std::pair<long, long> > &alignment,
                         std::vector<double> &index, std::vector<double> &modifier, std::vector<double> &ppm,
                         double &modifier_mass);

        void index_ions(vector<double> &sub, vector<double> &index_y,
                        std::vector<double> &modifier,
                        double &modifier_mass);

        void Get_all_ions_peaks(vector<node> &ions,vector<pair<long, long> > &alignment,vector<double> &index,
                vector<double> &ppm,double minPPMValue);


        void unknown_Ptm_Search(vector<node> &c_n_ions, vector<node> &n_c_ions,
                                         vector<node> &merg_ions_c_n, vector<node> &merg_ions_n_c,
                                         vector<pair<long, long> > &alignment_c_n,
                                         vector<pair<long, long> > &alignment_n_c,
                                         vector<double> &reference, vector<double> &peer,
                                         vector<double> &index_c_n, vector<double> &index_n_c,
                                         vector<double> &ppm_c_n, vector<double> &ppm_n_c,
                                         vector<modi> &mod,
                                         const vector<double> &modify_table,
                                         const map<double,string> &modMap,
                                         int lenth, double modifier_mass);

        void index_ions_Unknown(vector<double> &sub, vector<double> &index_y,
                                         std::vector<double> &modifier,
                                         double &modifier_mass);
        void getAllPerfectMatchIons(std::vector<std::pair<long, long> > &alignment, vector<double> &sub,
                                    vector<double> &index, vector<double> &ppm);

        void N_terminal_Get_Varible_Interval(vector<node> &ions, vector<modi> &Interval, map<long, long> &mapDup,
                                             const vector<double> &modify_table, double PtmsMass, int lenth);

        void Get_ppm_Unknown_(vector<double> &reference, vector<double> &peer, std::vector<std::pair<long, long> > &alignment,
                                      std::vector<double> &index, std::vector<double> &modifier, std::vector<double> &ppm,
                                      double &modifier_mass);

        bool _PPM_H(double t, double m,double &minPPM);

        void index_Unknown_Ptm(vector<node> &longIons_N,vector<modi> &mod,const vector<double> &PtmMass);

        /**
         * 进行推断修饰质量少的位置设定
         * @param mod
         * @param modifier
         * @param mMass
         */
        void assumMissingMass(vector<node> &longIons_N,vector<modi> &mod,
                              const std::vector<double> & modifier,double mMass,int len);


        void GetPpmUnknownReWrite(vector<double> &reference, vector<double> &peer,
                                           std::vector<std::pair<long, long> > &alignment,
                                           std::vector<double> &index, std::vector<double> &modifier, std::vector<double> &ppm,
                                           double &modifier_mass);

        void indexIonsReWrite(vector<double> &sub, vector<double> &index_y,
                                       std::vector<double> &modifier,
                                       double &modifier_mass);
    }
}
#endif