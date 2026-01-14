#include "cwdtw.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <memory.h>
#include <climits>
#include <cfloat>
#include "fun.h"

using namespace std;

void mylib::cwdtw::CWTAnalysis(
        const std::vector<double> &raw,
        std::vector<std::vector<double> > &output,
        double scale0, double dscale, long npyr) {
    const double *sigs = &raw[0];        //sst_nino3.dat
    cwt_object wt;

    long N = raw.size();
    double dt = 1;//2;		//sample rate	>  maybe we should use 2?

    wt = cwt_init("dog", 2.0, N, dt, npyr);    //"morlet", "dog", "paul"
    setCWTScales(wt, scale0, dscale, "pow", 2.0);
    cwt(wt, sigs);

    output.resize(npyr);
    for (long k = npyr; k--;) {
        long idx = npyr - k - 1;
        output[idx].resize(raw.size());
        long offset = k * raw.size();
        for (long i = 0; i < output[idx].size(); i++) {
            output[idx][i] = wt->output[i + offset].re;

        }
    }
    cwt_free(wt);
}


//----------------- boundary generation for constrained dynamic time warping (cDTW) -------------//
void mylib::cwdtw::BoundGeneration(
        std::vector<std::pair<long, long> > &cosali,
        long neib, std::vector<std::pair<long, long> > &bound, int mode,
        int RENMIN_or_SHENG) {

//if(mode!=-1) //-> mode = -1 means Renmin mode
    vector<pair<long, long> > cosali_ = cosali;
    if (mode == 1) //-> use partial-diagonol alignment
    {

        //-------- generate partial-diagonol alignment -------//
        vector<pair<long, long> > cosali2;
        cosali2.push_back(cosali[0]);
        for (long i = 1; i < cosali.size(); i++) {
            if (cosali[i].first != cosali2[cosali2.size() - 1].first) {
                cosali2.push_back(cosali[i]);
            } else {
                cosali2[cosali2.size() - 1].second = cosali[i].second;
            }
        }

        cosali_.clear();
        for (long i = 1; i < cosali2.size(); i++) {
            long pre_idx = cosali2[i - 1].first, cur_idx = cosali2[i].first;
            long pre_anchor = cosali2[i - 1].second, cur_anchor = cosali2[i].second;
            double anchor_diff = (cur_anchor - pre_anchor) / (cur_idx - pre_idx);
            for (long k = pre_idx, count = 0; k < cur_idx; k++, count++) {
                long mid = pre_anchor + (long) (count * anchor_diff);  //assign point relationship
                if (mid < pre_anchor)mid = pre_anchor;
                if (mid > cur_anchor)mid = cur_anchor;
                cosali_.push_back(make_pair(k, mid));
            }
        }
        cosali_.push_back(cosali2[cosali2.size() - 1]);
    }

    //---------- use block bound ------------//start
    //-> get signal length
    long moln1 = cosali_[cosali_.size() - 1].first + 1;
    long moln2 = cosali_[cosali_.size() - 1].second + 1;
    //-> renmin align to sheng style
    std::vector<std::pair<long, long> > align_sheng;
    mylib::fun::Renmin_To_Sheng_align(moln1, moln2, cosali_, align_sheng);
    //-> get bound in sheng style
    std::vector<std::pair<long, long> > bound_sheng;
    mylib::fun::From_Align_Get_Bound(moln1, moln2, align_sheng, bound_sheng, neib);
    //-> transfer bound to renmin style
    mylib::fun::Sheng_To_Renmin_bound(moln1, moln2, bound_sheng, bound);
    //----- maybe useful -----//
    bound[0].first = 0;
    bound[bound.size() - 1].first = bound[bound.size() - 2].first;
    bound[bound.size() - 1].second = cosali[cosali.size() - 1].second;


    if (RENMIN_or_SHENG == 1) //-> use Sheng's bound definition
    {
        bound = bound_sheng;
    }


    return;

}

//====================== continous wavelet dynamic time warping (cwDTW) ========================//
void mylib::cwdtw::MultiLevel_WaveletDTW(
        std::vector<double> &in1, std::vector<double> &in2,
        std::vector<std::vector<double> > &sig1,
        std::vector<std::vector<double> > &sig2,
        std::vector<std::pair<long, long> > &alignment,
        long radius, int test, int mode,
        double *totaldiff) {
    double tdiff;
    std::vector<std::pair<long, double> > sig1peaks, sig2peaks;
    double length1 = sig1[0].size(), length2 = sig2[0].size();

    long tot_size = sig1.size();
    for (long k = 0; k < tot_size; k++) {
        //------ peakpick CWT signal -------------//
        mylib::fun::PeakPick(sig1[k], sig1peaks);
        mylib::fun::PeakPick(sig2[k], sig2peaks);
//		std::cout<<sig1peaks.size()<<"\t"<<sig2peaks.size()<<std::endl;
        std::vector<double> peak1(sig1peaks.size());
        std::vector<double> peak2(sig2peaks.size());


        //-------- use peak picking result ---------//
        for (long i = sig1peaks.size(); i--;) {
            peak1[i] = sig1peaks[i].second;
        }
        for (long i = sig2peaks.size(); i--;) {
            peak2[i] = sig2peaks[i].second;
        }


//================================ peak average =====================//start
        if (test == 2)  //-> peak_ave
        {
            //-------- use peak average reso --------------//(just for comparison)
            for (long i = 1; i < sig1peaks.size() - 1; i++) {
                long start = (sig1peaks[i - 1].first + sig1peaks[i].first) / 2;
                long end = (sig1peaks[i].first + sig1peaks[i + 1].first) / 2;
                double val = 0;
                long cc;
                for (long j = start; j <= end; j++) {
                    val += in1[j];
                    cc++;
                }
                peak1[i] = val;
            }
            for (long i = 1; i < sig2peaks.size() - 1; i++) {
                long start = (sig2peaks[i - 1].first + sig2peaks[i].first) / 2;
                long end = (sig2peaks[i].first + sig2peaks[i + 1].first) / 2;
                double val = 0;
                long cc;
                for (long j = start; j <= end; j++) {
                    val += in2[j];
                    cc++;
                }
                peak2[i] = val;
            }
        }
//================================ peak average =====================//over


//================================ equal average =====================//start
        if (test == 1)  //-> equal_ave
        {
            //------ using Equally Averaging ------//(just for comparison)
            //--> proc peak1
            long cur1 = 0;
            long num1 = (long) (1.0 * in1.size() / sig1peaks.size());
            for (long i = 0; i < sig1peaks.size() - 1; i++) {
                double val = 0;
                for (long j = 0; j < num1; j++) {
                    val += in1[cur1];
                    cur1++;
                }
                val /= num1;
                sig1peaks[i].first = i * num1;
//		sig1peaks[i].first=i*num1+num1/2 < in1.size()-1? i*num1+num1/2:in1.size()-1 ;
//		sig1peaks[i].second=val;
                peak1[i] = val;
            }
            {
                long i = sig1peaks.size() - 1;
                long cc = 0;
                double val = 0;
                for (long j = cur1; j < in1.size(); j++) {
                    val += in1[j];
                    cc++;
                }
                val /= cc;
                sig1peaks[i].first = in1.size() - 1;
//		sig1peaks[i].first=i*num1+num1/2 < in1.size()-1? i*num1+num1/2:in1.size()-1 ;
//		sig1peaks[i].second=val;
                peak1[i] = val;
            }
            //--> proc peak2
            long cur2 = 0;
            long num2 = (long) (1.0 * in2.size() / sig2peaks.size());
            for (long i = 0; i < sig2peaks.size() - 1; i++) {
                double val = 0;
                for (long j = 0; j < num2; j++) {
                    val += in2[cur2];
                    cur2++;
                }
                val /= num2;
                sig2peaks[i].first = i * num2;
//		sig2peaks[i].first=i*num2+num2/2 < in2.size()-1? i*num2+num2/2:in2.size()-1 ;
//		sig2peaks[i].second=val;
                peak2[i] = val;
            }
            {
                long i = sig2peaks.size() - 1;
                long cc = 0;
                double val = 0;
                for (long j = cur2; j < in2.size(); j++) {
                    val += in2[j];
                    cc++;
                }
                val /= cc;
                sig2peaks[i].first = in2.size() - 1;
//		sig2peaks[i].first=i*num2+num2/2 < in2.size()-1? i*num2+num2/2:in2.size()-1 ;
//		sig2peaks[i].second=val;
                peak2[i] = val;
            }
        }
//================================ equal average =====================//over

        //----- apply DTW or cDTW dependent on k-th level -------//
        if (k == 0) {
            tdiff = mylib::fun::DynamicTimeWarping(peak1, peak2, alignment);
        } else {
            //----- ReMapIndex_partI (map ground level upon k-th level) -----//
            long c = 0;
            for (long i = 0; i < alignment.size(); i++) {
                while (sig1peaks[c].first < alignment[i].first) {
                    c++;
                }
                alignment[i].first = c;
            }

            long d = 0;
            for (long i = 0; i < alignment.size(); i++) {
                while (sig2peaks[d].first < alignment[i].second) {
                    d++;
                }
                alignment[i].second = d;
            }
            //----- cDWT (constrained DWT) -------//
            std::vector<std::pair<long, long> > bound;
            long neib = radius * pow(2, tot_size - k); //adaptive radius
            BoundGeneration(alignment, neib, bound, mode);
            tdiff = mylib::fun::BoundDynamicTimeWarping(peak1, peak2, bound, alignment);
        }
        //----- ReMapIndex_partII (map k-th level back to ground level) -----//
        for (long i = alignment.size(); i--;) {
            alignment[i].first = sig1peaks[alignment[i].first].first;
            alignment[i].second = sig2peaks[alignment[i].second].first;
        }

    }

    if (totaldiff) {
        *totaldiff = tdiff;
    }
}