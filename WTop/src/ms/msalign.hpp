//
// Created by Administrator on 2023/3/18.
//

#ifndef WTOP_MSALIGN_HPP
#define WTOP_MSALIGN_HPP

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <memory>

#include "../prsm/candidate_prsm.hpp"

#include "../mylib/data_steam.hpp"
#include "../util/string_utils.h"

#include "../protein/protein.hpp"
#include "../filter/filter_tag_msg.hpp"

#include "../merge/modify.hpp"


using namespace std;

class msalign {
public:
    void init_ms_msg(ifstream &ms_file);

    void insert_ions_msg_container(std::string &msg) {
        ions_msg_container.push_back(msg);
    }

    void initCompIonsIndexMap();

    const std::map<std::string, std::string> &getMsgMap() const;

    const std::vector<std::string> &getIonsMsgContainer() const;

    std::vector<double> &getIonsMassContainer();

    void setPrecursorMass(double precursorMass);

    double getPrecursorMass() const;

    const map<int, int> &getCompIonsIndexMap() const;

    const map<double, double> &getKeyMsValueAbRadio() const;

    const map<double, double> &getKeyMsValueCharge() const;

    const map<std::string, std::string> &getMsMsgMap() const;

    const string &getId() const;

    const string &getScans() const;

    const string &getRetentionTime() const;

    const string &getActivation() const;

    const string &getActivationSub() const;

    const string &getMsOneId() const;

    const string &getMsOneScans() const;

    const string &getPrecursorMz() const;

    const string &getPrecursorCharge() const;

    const string &getPrecursorMassStr() const;

    const string &getPrecursorIntensity() const;

    const vector<protein_ptr> &getCandidateProteins() const;

    void insert_protein_ptr(protein_ptr proteinPtr) {
        candidate_proteins.push_back(proteinPtr);
    }

    double getPrsmBestScore() const;

    void setPrsmBestScore(double prsmBestScore);

    int getPrsmMatchPeaks() const;

    const string &getPrsmMatchProteinName() const;

    string &getPrsmMatchProteinSeq();

    int getPrsmCutN() const;

    int getPrsmCutC() const;

    const vector<modi> &getPrsmPtm() const;

    double getPrsmMutX() const;

    int getPrsmIonsPairNum() const;

    int getPrsmComplementsIonsNum() const;

    double getPrsmSeqScore() const;

    int getPrsmSeqScoreCount() const;

    int getPrsmMatchFragmentIonsNum() const;

    void setPrsmMatchPeaks(int prsmMatchPeaks);

    void setPrsmMatchProteinName(const string &prsmMatchProteinName);

    void setPrsmMatchProteinSeq(const string &prsmMatchProteinSeq);

    void setPrsmCutN(int prsmCutN);

    void setPrsmCutC(int prsmCutC);

    void setPrsmPtm(const vector<modi> &prsmPtm);

    void setPrsmMutX(double prsmMutX);

    void setPrsmIonsPairNum(int prsmIonsPairNum);

    void setPrsmComplementsIonsNum(int prsmComplementsIonsNum);

    void setPrsmSeqScore(double prsmSeqScore);

    void setPrsmSeqScoreCount(int prsmSeqScoreCount);

    void setPrsmMatchFragmentIonsNum(int prsmMatchFragmentIonsNum);

    double getPrsmPtmsMass() const;

    const map<string, string> &getPrsmUnknownPtmMassMap() const;

    prsm::modify_ptr getPrsmUnknownModPtr();

    double getPrsmTheoryProteoformMass() const;

    void setPrsmPtmsMass(double prsmPtmsMass);

    void setPrsmUnknownPtmMassMap(const map<string, string> &prsmUnknownPtmMassMap);

    void setPrsmUnknownModPtr(prsm::modify_ptr prsmUnknownModPtr);

    void setPrsmTheoryProteoformMass(double prsmTheoryProteoformMass);

    std::vector<protein_ptr> candidate_proteins;

    std::vector<node> n_ter;
    std::vector<node> c_ter;
    //test_outcandi_prsm
    vector<candidate_prsm_ptr> candi_prsms;

private:
    std::map<int, int> compIonsIndexMap;           //存放互补离子下标
    std::map<double, double> keyMsValueAbRadio;    // ions mass - radio
    std::map<std::string, std::string> ms_msg_map; // title msg
    std::vector<std::string> ions_msg_container;   //
    std::vector<double> ions_mass_container;
    std::map<double, double> keyMsValueCharge;    // ions mass - Charge

    std::string id; //"ID=XX"
    std::string scans;
    std::string retention_time;
    std::string activation;
    std::string activation_sub;
    std::string ms_one_id;
    std::string ms_one_scans;
    std::string precursor_mz;
    std::string precursor_charge;
    std::string precursor_mass_str;
    double precursor_mass;
    std::string precursor_intensity;

    double prsm_best_score = 0;
    int prsm_match_peaks;
    std::string prsm_match_protein_name;
    std::string prsm_match_protein_seq;
    int prsm_cut_n;
    int prsm_cut_c;
    std::vector<modi> prsm_ptm;
    double prsm_mut_x;
    int prsm_ionsPair_num;
    int prsm_complementsIons_num;
    double prsm_seqScore;
    int prsm_seqScore_count;
    int prsm_match_fragmentIons_num;
    double prsm_ptms_mass;
    map<string, string> prsm_unknown_ptmMass_map;

    prsm::modify_ptr prsm_unknownModPtr;
    double prsm_theory_proteoform_mass;


    ///调试打分特征变量
    double feature1;
    double feature2;
    double feature3;
    double feature4;
    double feature5;
    double feature6;
    double feature7;
    double feature8;
    double feature9;
    double feature10;
    double feature11;
    double feature12;
    double feature13;
    double feature14;
    double feature15;
    double feature16;
    double feature17;
    double feature18;
    double feature19;
    double feature20;
    double feature21;
    double feature22;
    double feature23;
    double feature24;
    double feature25;
    double feature26;
    double feature27;
    double feature28;
    double feature29;
    double feature30;
    double feature31;
    double feature32;
    double feature33;
    double feature34;
    double feature35;
    double feature36;
    double feature37;
    double feature38;
    double feature39;
    double feature40;
    double feature41;
    double feature42;
    double feature43;
    double feature44;
    double feature45;
    double feature46;
    double feature47;
    double feature48;
    double feature49;
    double feature50;



public:
    double getFeature1() const;
    void setFeature1(double feature1);
    double getFeature2() const;
    void setFeature2(double feature2);
    double getFeature3() const;
    void setFeature3(double feature3);
    double getFeature4() const;
    void setFeature4(double feature4);
    double getFeature5() const;
    void setFeature5(double feature5);
    double getFeature6() const;
    void setFeature6(double feature6);
    double getFeature7() const;
    void setFeature7(double feature7);
    double getFeature8() const;
    void setFeature8(double feature8);
    double getFeature9() const;
    void setFeature9(double feature9);
    double getFeature10() const;
    void setFeature10(double feature10);
    double getFeature11() const;
    void setFeature11(double feature11);
    double getFeature12() const;
    void setFeature12(double feature12);
    double getFeature13() const;
    void setFeature13(double feature13);
    double getFeature14() const;
    void setFeature14(double feature14);
    double getFeature15() const;
    void setFeature15(double feature15);
    double getFeature16() const;
    void setFeature16(double feature16);
    double getFeature17() const;
    void setFeature17(double feature17);
    double getFeature18() const;
    void setFeature18(double feature18);
    double getFeature19() const;
    void setFeature19(double feature19);
    double getFeature20() const;
    void setFeature20(double feature20);
    double getFeature21() const;
    void setFeature21(double feature21);
    double getFeature22() const;
    void setFeature22(double feature22);
    double getFeature23() const;
    void setFeature23(double feature23);
    double getFeature24() const;
    void setFeature24(double feature24);
    double getFeature25() const;
    void setFeature25(double feature25);
    double getFeature26() const;
    void setFeature26(double feature26);
    double getFeature27() const;
    void setFeature27(double feature27);
    double getFeature28() const;
    void setFeature28(double feature28);
    double getFeature29() const;
    void setFeature29(double feature29);
    double getFeature30() const;
    void setFeature30(double feature30);
    double getFeature31() const;
    void setFeature31(double feature31);
    double getFeature32() const;
    void setFeature32(double feature32);
    double getFeature33() const;
    void setFeature33(double feature33);
    double getFeature34() const;
    void setFeature34(double feature34);
    double getFeature35() const;
    void setFeature35(double feature35);
    double getFeature36() const;
    void setFeature36(double feature36);
    double getFeature37() const;
    void setFeature37(double feature37);
    double getFeature38() const;
    void setFeature38(double feature38);
    double getFeature39() const;
    void setFeature39(double feature39);
    double getFeature40() const;
    void setFeature40(double feature40);
    double getFeature41() const;
    void setFeature41(double feature41);
    double getFeature42() const;
    void setFeature42(double feature42);
    double getFeature43() const;
    void setFeature43(double feature43);
    double getFeature44() const;
    void setFeature44(double feature44);
    double getFeature45() const;
    void setFeature45(double feature45);
    double getFeature46() const;
    void setFeature46(double feature46);
    double getFeature47() const;
    void setFeature47(double feature47);
    double getFeature48() const;
    void setFeature48(double feature48);
    double getFeature49() const;
    void setFeature49(double feature49);
    double getFeature50() const;
    void setFeature50(double feature50);


};

typedef shared_ptr<msalign> msalign_ptr;
typedef vector<msalign_ptr> msalignPtrVec;

#endif //WTOP_MSALIGN_HPP
