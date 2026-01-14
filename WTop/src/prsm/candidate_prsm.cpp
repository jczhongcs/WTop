#include "candidate_prsm.hpp"
//
// Created by LYC on 2023/7/14.
//

#include "candidate_prsm.hpp"

double candidate_prsm::getPrsmBestScore() const {
    return prsm_best_score;
}
void candidate_prsm::setPrsmBestScore(double prsmBestScore) {
    prsm_best_score = prsmBestScore;
}
int candidate_prsm::getPrsmMatchPeaks() const {
    return prsm_match_peaks;
}
void candidate_prsm::setPrsmMatchPeaks(int prsmMatchPeaks) {
    prsm_match_peaks = prsmMatchPeaks;
}
const std::string &candidate_prsm::getPrsmMatchProteinName() const {
    return prsm_match_protein_name;
}
void candidate_prsm::setPrsmMatchProteinName(const std::string &prsmMatchProteinName) {
    prsm_match_protein_name = prsmMatchProteinName;
}
const std::string &candidate_prsm::getPrsmMatchProteinSeq() const {
    return prsm_match_protein_seq;
}
void candidate_prsm::setPrsmMatchProteinSeq(const std::string &prsmMatchProteinSeq) {
    prsm_match_protein_seq = prsmMatchProteinSeq;
}
int candidate_prsm::getPrsmCutN() const {
    return prsm_cut_n;
}
void candidate_prsm::setPrsmCutN(int prsmCutN) {
    prsm_cut_n = prsmCutN;
}
int candidate_prsm::getPrsmCutC() const {
    return prsm_cut_c;
}
void candidate_prsm::setPrsmCutC(int prsmCutC) {
    prsm_cut_c = prsmCutC;
}
const std::vector<modi> &candidate_prsm::getPrsmPtm() const {
    return prsm_ptm;
}
void candidate_prsm::setPrsmPtm(const std::vector<modi> &prsmPtm) {
    prsm_ptm = prsmPtm;
}
double candidate_prsm::getPrsmMutX() const {
    return prsm_mut_x;
}
void candidate_prsm::setPrsmMutX(double prsmMutX) {
    prsm_mut_x = prsmMutX;
}
int candidate_prsm::getPrsmIonsPairNum() const {
    return prsm_ionsPair_num;
}
void candidate_prsm::setPrsmIonsPairNum(int prsmIonsPairNum) {
    prsm_ionsPair_num = prsmIonsPairNum;
}
int candidate_prsm::getPrsmComplementsIonsNum() const {
    return prsm_complementsIons_num;
}
void candidate_prsm::setPrsmComplementsIonsNum(int prsmComplementsIonsNum) {
    prsm_complementsIons_num = prsmComplementsIonsNum;
}
double candidate_prsm::getPrsmSeqScore() const {
    return prsm_seqScore;
}
void candidate_prsm::setPrsmSeqScore(double prsmSeqScore) {
    prsm_seqScore = prsmSeqScore;
}
int candidate_prsm::getPrsmSeqScoreCount() const {
    return prsm_seqScore_count;
}
void candidate_prsm::setPrsmSeqScoreCount(int prsmSeqScoreCount) {
    prsm_seqScore_count = prsmSeqScoreCount;
}
int candidate_prsm::getPrsmMatchFragmentIonsNum() const {
    return prsm_match_fragmentIons_num;
}
void candidate_prsm::setPrsmMatchFragmentIonsNum(int prsmMatchFragmentIonsNum) {
    prsm_match_fragmentIons_num = prsmMatchFragmentIonsNum;
}
double candidate_prsm::getPrsmPtmsMass() const {
    return prsm_ptms_mass;
}
void candidate_prsm::setPrsmPtmsMass(double prsmPtmsMass) {
    prsm_ptms_mass = prsmPtmsMass;
}
const map<string, string> &candidate_prsm::getPrsmUnknownPtmMassMap() const {
    return prsm_unknown_ptmMass_map;
}
void candidate_prsm::setPrsmUnknownPtmMassMap(const map<string, string> &prsmUnknownPtmMassMap) {
    prsm_unknown_ptmMass_map = prsmUnknownPtmMassMap;
}
const prsm::modify_ptr &candidate_prsm::getPrsmUnknownModPtr() const {
    return prsm_unknownModPtr;
}
void candidate_prsm::setPrsmUnknownModPtr(const prsm::modify_ptr &prsmUnknownModPtr) {
    prsm_unknownModPtr = prsmUnknownModPtr;
}
double candidate_prsm::getPrsmTheoryProteoformMass() const {
    return prsm_theory_proteoform_mass;
}
void candidate_prsm::setPrsmTheoryProteoformMass(double prsmTheoryProteoformMass) {
    prsm_theory_proteoform_mass = prsmTheoryProteoformMass;
}
