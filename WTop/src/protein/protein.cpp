//
// Created by Administrator on 2023/3/18.
//

#include "protein.hpp"

std::map<char, double> mapMass=
        {
                {'G',57.02146},
                {'A',71.03711},
                {'S',87.03202},
                {'P',97.05276},
                {'V',99.06841},
                {'T',101.04768},
                {'C',103.00918},
                {'I',113.08406},
                {'L',113.08406},
                {'N',114.04293},
                {'D',115.02694},
                {'Q',128.05858},
                {'K',128.095},
                {'E',129.0426},
                {'M',131.0405},
                {'H',137.05891},
                {'F',147.06841},
                {'R',156.10111},
                {'Y',163.06333},
                {'W',186.079354},
                {'U',168.96420},
                {'X',0.0},
                {'O',150.953636},
        };

double protein::getAcidSeqMassN(string & acidSeq)
{
    double mass = 0;
    for (char c : acidSeq) {
        mass += mapMass[c];
    }
    return mass ;
}
double protein::getNCutMass(int cutSize)
{
    string seq = protein_sequence ;
    double mass = 0 ;
    for(int i = 0 ; i < cutSize ; ++i) {
        mass += mapMass[protein_sequence[i]];
    }
    return mass ;
}


map<char,double> protein::get_acid_Map() {
    return mapMass;
}

void protein::initTheoryMass() {
    double massAdd = 0.0;
    //C端 Y ions and ZIons
    for(int xi = protein_sequence.size() - 1; xi >= 0 ; xi--)  {
        if(mapMass.count(protein_sequence[xi])) {     //C端构造理论峰
            massAdd += mapMass[protein_sequence[xi]];
            theoryMassC_By.push_back(massAdd + 18.0106);
            theoryMassC_Cz.push_back(massAdd + 1.9919);
        }
    }
    massAdd = 0.0;
    double czAddMass = 17.0265;
    //N 端 B ions and C ions
    for(int xi = 0 ; xi < protein_sequence.size(); ++xi ) {   //N端构造理论峰
        if(mapMass.count(protein_sequence[xi])) {
            massAdd += mapMass[protein_sequence[xi]];
            theoryMassN_By.push_back(massAdd);
            theoryMassN_Cz.push_back(massAdd);
        }
    }
    for(int xi = 0 ; xi < theoryMassN_Cz.size(); ++xi ) { //
        theoryMassN_Cz[xi] += czAddMass;
    }
}

void protein::get_seq_mass(vector<double> & seq_mass)
{
    for(int i = 0 ; i < protein_sequence.size() ;++i) {
        seq_mass.push_back(mapMass[protein_sequence[i]]);
    }
}

protein::protein(const string &proteinName, const string &proteinSequence) : protein_name(proteinName),
                                                                             protein_sequence(proteinSequence) {
//    if (protein_sequence.find("X") != std::string::npos) {
//        protein_sequence.erase(std::remove(protein_sequence.begin(), protein_sequence.end(),'X'),
//                               protein_sequence.end());
//    }
    if (protein_sequence.find("*") != std::string::npos) {
        protein_sequence.erase(std::remove(protein_sequence.begin(), protein_sequence.end(),'*'),
                               protein_sequence.end());
    }
    protein::initTheoryMass( );
}

const string &protein::getProteinName() const {
    return protein_name;
}

const string &protein::getProteinSequence() const {
    return protein_sequence;
}

 vector<double> &protein::getTheoryMassCBy()  {
    return theoryMassC_By;
}

const vector<double> &protein::getTheoryMassNBy() const {
    return theoryMassN_By;
}

const vector<double> &protein::getTheoryMassCCz() const {
    return theoryMassC_Cz;
}

const vector<double> &protein::getTheoryMassNCz() const {
    return theoryMassN_Cz;
}

const vector<double> &protein::getSeqMass() const {
    return seqMass;
}
