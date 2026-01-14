//
// Created by Administrator on 2023/3/18.
//

#ifndef WTOP_PROTEIN_HPP
#define WTOP_PROTEIN_HPP


#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;
class protein {
public:

    protein(){

    };
    protein(const string &proteinName, const string &proteinSequence);

    void initTheoryMass() ;

    double getAcidSeqMassN(string & acidSeq);

    void get_seq_mass(vector<double> & seq_mass) ;

    //得到N端截断个数的质量
    double getNCutMass(int cutSize);

    std::map<char,double> get_acid_Map();

    const string &getProteinName() const;

    const string &getProteinSequence() const;

     vector<double> &getTheoryMassCBy() ;

    const vector<double> &getTheoryMassNBy() const;

    const vector<double> &getTheoryMassCCz() const;

    const vector<double> &getTheoryMassNCz() const;

    const vector<double> &getSeqMass() const;

    string getProteinTitle() {
        return protein_name.substr(0,protein_name.find_first_of(" "));
    }

private:
    string protein_name;     //蛋白质名称
    string protein_sequence;      //蛋白质序列

    vector<double> theoryMassC_By;      //理论质量
    vector<double> theoryMassN_By;      //理论质量

    vector<double> theoryMassC_Cz;      //理论质量
    vector<double> theoryMassN_Cz;      //理论质量

    vector<double> seqMass ;

};
typedef protein* protein_ptr ;

#endif //WTOP_PROTEIN_HPP
