#include "protein_processor.hpp"

std::map <char, double> map_mass=
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
    };

namespace prsm{

    double Protein::getNCutMass(int cutSize)
    {
        string seq = protein_sequence ;
        double mass = 0 ;
        for(int i = 0 ; i < cutSize ; ++i) {
            mass += map_mass[protein_sequence[i]];
        }
        return mass ;
    }

    void Protein::Init_Sequence(ifstream &prote_file_ptr,int &first_flag)
    {
    std::string line;
    int flag = 0;
    if(first_flag==0)
    {
        getline(prote_file_ptr,first_title);
        mylib::data_stream::rtrim_s(first_title);
        first_flag++;
    }
    else
    {
        first_title = next_title;
    }
    while(!prote_file_ptr.eof() && flag<1) {    //每次只读取一条蛋白质序列

        getline(prote_file_ptr,line);
        mylib::data_stream::rtrim_s(line);
        if (line[0] == '>') 
        {
            next_title = line;
            flag++;
        }
        else if (line[0] != '>'){
            sequence += line;
        }
    }
    }

    void Protein::initTitleSeq(const string &title,const string &seq) {
        this->protein_name = title;
        this->protein_sequence = seq ;
    }

    map<char,double> Protein::getMap() {
        return map_mass;
    }

    void Protein::initTheoryMass() {
        double massAdd = 0.0;
        //C端 Y ions and ZIons
        for(int xi = protein_sequence.size() - 1; xi >= 0 ; xi--)  {
            if(map_mass.count(protein_sequence[xi])) {     //C端构造理论峰
                massAdd += map_mass[protein_sequence[xi]];
                theoryMassC_By.push_back(massAdd + 18.01056);
                theoryMassC_Cz.push_back(massAdd + 1.9919);
            }
        }
        massAdd = 0.0;
        double czAddMass = 17.026;
        //N 端 B ions and C ions
        for(int xi = 0 ; xi < protein_sequence.size(); ++xi ) {   //N端构造理论峰
            if(map_mass.count(protein_sequence[xi])) {
                massAdd += map_mass[protein_sequence[xi]];
                theoryMassN_By.push_back(massAdd);
                theoryMassN_Cz.push_back(massAdd);
            }
        }
        for(int xi = 0 ; xi < theoryMassN_Cz.size(); ++xi ) { //
            theoryMassN_Cz[xi] += czAddMass;
        }
    }

    void Protein::get_seq_mass(string seq , vector<double> & seq_mass)
    {
        if (seq.size() == 0 ) {
            return ;
        }
        for(int i = 0 ; i < seq.size() ;++i) {
            seq_mass.push_back(map_mass[seq[i]]);
        }
    }

    Protein::Protein(const string &proteinName, const string &proteinSequence) : protein_name(proteinName),
                                                                                 protein_sequence(proteinSequence) {

        if (protein_sequence.find("X")!=std::string::npos) {
            protein_sequence.erase(std::remove(protein_sequence.begin(), protein_sequence.end(), 'X'), protein_sequence.end());
        }
        Protein::initTheoryMass();
    }
}