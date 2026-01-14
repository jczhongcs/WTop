//
// Created by LYC on 2023/7/14.
//
#ifndef WTOP_CANDIDATE_PRSM_HPP
#define WTOP_CANDIDATE_PRSM_HPP

#include <string>
#include <vector>
#include "../merge/modify.hpp"
#include "../mylib/data_steam.hpp"





class candidate_prsm{
private:


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

public:
    double getPrsmBestScore() const;
    void setPrsmBestScore(double prsmBestScore);
    int getPrsmMatchPeaks() const;
    void setPrsmMatchPeaks(int prsmMatchPeaks);
    const std::string &getPrsmMatchProteinName() const;
    void setPrsmMatchProteinName(const std::string &prsmMatchProteinName);
    const std::string &getPrsmMatchProteinSeq() const;
    void setPrsmMatchProteinSeq(const std::string &prsmMatchProteinSeq);
    int getPrsmCutN() const;
    void setPrsmCutN(int prsmCutN);
    int getPrsmCutC() const;
    void setPrsmCutC(int prsmCutC);
    const std::vector<modi> &getPrsmPtm() const;
    void setPrsmPtm(const std::vector<modi> &prsmPtm);
    double getPrsmMutX() const;
    void setPrsmMutX(double prsmMutX);
    int getPrsmIonsPairNum() const;
    void setPrsmIonsPairNum(int prsmIonsPairNum);
    int getPrsmComplementsIonsNum() const;
    void setPrsmComplementsIonsNum(int prsmComplementsIonsNum);
    double getPrsmSeqScore() const;
    void setPrsmSeqScore(double prsmSeqScore);
    int getPrsmSeqScoreCount() const;
    void setPrsmSeqScoreCount(int prsmSeqScoreCount);
    int getPrsmMatchFragmentIonsNum() const;
    void setPrsmMatchFragmentIonsNum(int prsmMatchFragmentIonsNum);
    double getPrsmPtmsMass() const;
    void setPrsmPtmsMass(double prsmPtmsMass);
    const map<string, string> &getPrsmUnknownPtmMassMap() const;
    void setPrsmUnknownPtmMassMap(const map<string, string> &prsmUnknownPtmMassMap);
    const prsm::modify_ptr &getPrsmUnknownModPtr() const;
    void setPrsmUnknownModPtr(const prsm::modify_ptr &prsmUnknownModPtr);
    double getPrsmTheoryProteoformMass() const;
    void setPrsmTheoryProteoformMass(double prsmTheoryProteoformMass);
};
typedef candidate_prsm *candidate_prsm_ptr;

#endif //WTOP_CANDIDATE_PRSM_HPP