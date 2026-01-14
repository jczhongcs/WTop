//
// Created by LYC on 2025/4/3.
//

#ifndef WTOP_CANDI_TRUNCATION_H
#define WTOP_CANDI_TRUNCATION_H
#include <memory>
#include <vector>


using namespace std;
class candi_trucation{


public:
    int getCutN() const;
    void setCutN(int cutN);
    int getCutC() const;
    void setCutC(int cutC);
    double getUnkmass() const;
    void setUnkmass(double unkmass);
    int getMatchPeaks() const;
    void setMatchPeaks(int matchPeaks);
    double getVarptmmass() const;
    void setVarptmmass(double varptmmass);
    candi_trucation(int cutN, int cutC, double unkmass, double varptmmass, int matchPeaks);

private:
    int cut_n;
    int cut_c;
    double unkmass;
    double varptmmass;
    int match_peaks;
};

typedef shared_ptr<candi_trucation> candi_trucation_ptr;
typedef vector<candi_trucation_ptr> candi_trucation_vec;

#endif //WTOP_CANDI_TRUNCATION_H