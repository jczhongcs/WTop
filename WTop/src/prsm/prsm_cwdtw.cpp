//
// Created by wenzhong on 2023/3/22.
//

#include "prsm_cwdtw.hpp"

double prsm_cwdtw::cwdtw_algorithm(vector<double> &reference_1, vector<double> &peer_1,
                                     std::vector<std::pair<long, long> > &alignment,
                                     vector<double> &ref_zscore, vector<double> &peer_zscore) {
    struct options {
        //-> required parameter
        char output[65532];
        //-> key parameter
        int radius;
        int level;
        float scale0;
        //-> vice parameter
        int verbose;
        int test;
        int mode;
    };
    options opts;
    opts.radius = 50;
    opts.level = 1;
    opts.scale0 = sqrt(2);
    opts.verbose = 0;       //-> [0] no verbose; 1 verbose
    opts.test = 0;       //-> [0] not use test mode; 1 equal_ave, 2 peak_ave, 3 Fast_DTW
    opts.mode = 0;       //-> [0] block bound; 1 diagonol bound


    //======================= START Procedure ===================================//
    vector<double> reference = reference_1;
    std::vector<double> reference_orig = reference;

    vector<double> peer;

    for (int i = 0; i < peer_1.size(); i++) {
        peer.push_back(peer_1[i]);
    }
    std::vector<double> peer_orig = peer;
    //----- length check ------//
    int swap = 0;
    if (reference.size() > peer.size()) {
        std::vector<double> tmp = peer;
        peer = reference;
        reference = tmp;
        tmp = ref_zscore;
        ref_zscore = peer_zscore;
        peer_zscore = tmp;
        swap = 1;
    }

    //==================================================//
    //------3. process initial input signals ----------//
    //----- 3.1 Zscore normaliza on both signals -----//
    double avg, dev;

    //-------3.1 test consol Zsocrenormalize --------------//
    std::vector<double> consol;
    int i = 0;

    for (i = 0; i < peer.size(); i++) {
        consol.push_back(peer[i]);
    }
    for (i = 0; i < reference.size(); i++) {
        consol.push_back(reference[i]);
    }
    mylib::fun::ZScoreNormalize(consol, &avg, &dev);
    for (i = 0; i < peer.size(); i++) {
        peer[i] = consol[i];
    }
    for (int j = 0; i < consol.size(); i++, j++) {
        reference[j] = consol[i];
    }
    ref_zscore = reference;
    peer_zscore = peer;
    //----- 3.2 calculate length ratio between input signals -----//
    double alpha = (double) peer.size() / reference.size();


    //====================================================//
    //----- 4. continous wavelet transform --------------//
    std::vector<std::vector<double> > rcwt, pcwt;

    if (opts.verbose == 1) {
        cout << "CWT Analysis...\n" ;
    }

    long npyr = opts.level;          // default: 3
    double scale0 = opts.scale0;    // default: sqrt(2)
    double dscale = 1;              // default: 1

    mylib::cwdtw::CWTAnalysis(reference, rcwt, scale0, dscale, npyr);
    mylib::cwdtw::CWTAnalysis(peer, pcwt, scale0 * alpha, dscale, npyr);

    //------ 4.1 Zscore normaliza on both CWT signals -----//
    //if multiscale is used, pyr logical should be added.
    for (long i = 0; i < npyr; i++) {
        mylib::fun::ZScoreNormalize(rcwt[i], &avg, &dev);
        mylib::fun::ZScoreNormalize(pcwt[i], &avg, &dev);
    }

    //============================================//
    //------ 5. multi-level WaveletDTW ----------//
    std::vector<std::pair<long, long> > cosali;
    double tdiff;

    if (opts.verbose == 1) {
        EX_TRACE("Coarse Alignment...\n");
    }
    mylib::cwdtw::MultiLevel_WaveletDTW(reference, peer, rcwt, pcwt, cosali, opts.radius, opts.test, opts.mode, &tdiff);
    if (opts.verbose == 1) {
        EX_TRACE("Average Deviation (%.1lf/%ld=%.3lf)\n", tdiff, cosali.size(), tdiff / cosali.size());
    }

    //------ 5.1 generate final boundary -------//
    std::vector<std::pair<long, long> > bound;
    mylib::cwdtw::BoundGeneration(cosali, opts.radius, bound, opts.mode);

    //------ 5.2 generate final alignment via cDTW ------//
    // std::vector<std::pair<long,long> > alignment;
    tdiff = mylib::fun::BoundDynamicTimeWarping(reference, peer, bound, alignment);
    // fprintf(stderr,"%lf %d %lf\n",tdiff,alignment.size(),tdiff/alignment.size());

    double c1 = tdiff / alignment.size();

//=================================================//
    //------ 6. output final alignment to file -------//
    // if(output!="")
    // {
    // 	WriteSequenceAlignment(output, reference_orig, peer_orig, reference, peer, alignment, swap,tdiff);
    // }
    //----- exit -----//
//=================================================//
    if (swap == 1) {
        for (i = 0; i < alignment.size(); i++) {
            long t;
            t = alignment[i].first;
            alignment[i].first = alignment[i].second;
            alignment[i].second = t;
        }
        std::vector<double> tmp = peer;
        peer = reference;
        reference = tmp;
        tmp = ref_zscore;
        ref_zscore = peer_zscore;
        peer_zscore = tmp;
        swap = 1;
    }
    return tdiff;
}