//
// Created by Administrator on 2023/3/18.
//

#include "msalign.hpp"

const std::map<std::string, std::string> &msalign::getMsgMap() const {
    return ms_msg_map;
}

const std::vector<std::string> &msalign::getIonsMsgContainer() const {
    return ions_msg_container;
}

std::vector<double> &msalign::getIonsMassContainer()  {
    return ions_mass_container;
}

void msalign::setPrecursorMass(double precursorMass) {
    precursor_mass = precursorMass;
}

double msalign::getPrecursorMass() const {
    return precursor_mass;
}


void msalign::init_ms_msg(std::ifstream &ms_file) {
    if (!ms_file.good()) {
        std::cerr << "ms file read - failed" << std::endl;
        std::exit(1);
    }

    std::string buff;
    ms_msg_map.clear();                  // 清空 map，避免之前内容残留
    ions_mass_container.clear();        // 清空向量
    keyMsValueAbRadio.clear();
    keyMsValueCharge.clear();

    while (std::getline(ms_file, buff)) {
        mylib::data_stream::rtrim_s(buff);  // 清理尾部空白

        if (buff == "END IONS") {
            break;
        }

        size_t pos = buff.find('=');
        if (pos != std::string::npos) {
            // 使用 emplace 避免不必要的复制
            ms_msg_map.emplace(buff.substr(0, pos), buff.substr(pos + 1));
        } else {
            std::vector<std::string> s_str;
            utils_string_util::Stringsplit(buff, '\t', s_str);

            if (s_str.size() == 3) {
                try {
                    double mass = std::stod(s_str[0]);
                    ions_mass_container.emplace_back(mass);
                    keyMsValueAbRadio.emplace(mass, std::stod(s_str[1]));
                    keyMsValueCharge.emplace(mass, std::stod(s_str[2]));
                } catch (const std::invalid_argument &e) {
                    std::cerr << "Invalid number in MS line: " << buff << std::endl;
                }
            }
        }
    }

    // 使用 find 检查键是否存在，防止异常或不完整数据崩溃
    auto safe_get = [&](const std::string &key) -> std::string {
      auto it = ms_msg_map.find(key);
      return (it != ms_msg_map.end()) ? it->second : "";
    };

    id = "ID=" + safe_get("ID");
    scans = "SCANS=" + safe_get("SCANS");
    retention_time = "RETENTION_TIME=" + safe_get("RETENTION_TIME");
    activation = "ACTIVATION=" + safe_get("ACTIVATION");
    activation_sub = safe_get("ACTIVATION");
    ms_one_id = "MS_ONE_ID=" + safe_get("MS_ONE_ID");
    ms_one_scans = "MS_ONE_SCAN=" + safe_get("MS_ONE_SCAN");
    precursor_mz = "PRECURSOR_MZ=" + safe_get("PRECURSOR_MZ");
    precursor_charge = "PRECURSOR_CHARGE=" + safe_get("PRECURSOR_CHARGE");
    precursor_mass_str = "PRECURSOR_MASS=" + safe_get("PRECURSOR_MASS");

    try {
        precursor_mass = std::stod(safe_get("PRECURSOR_MASS"));
    } catch (...) {
        precursor_mass = 0.0;
    }

    precursor_intensity = "PRECURSOR_INTENSITY=" + safe_get("PRECURSOR_INTENSITY");

    std::sort(ions_mass_container.begin(), ions_mass_container.end());
}




void msalign::initCompIonsIndexMap()
{
    set<int> s ; //用于去重复
    for (int i = 0 ; i < ions_mass_container.size(); ++i) {
        double t = precursor_mass - ions_mass_container[i];
        int index = mylib::data_stream::txpb_binarey_search_ex(ions_mass_container, ions_mass_container.size(), t);
        if ( index>=0 && !s.count(index) && fabs(t - ions_mass_container[index]) < 2)
        {
            compIonsIndexMap.insert(make_pair(i,index));
            s.insert(i);
        }
    }
}

const map<int, int> &msalign::getCompIonsIndexMap() const {
    return compIonsIndexMap;
}

const map<double, double> &msalign::getKeyMsValueAbRadio() const {
    return keyMsValueAbRadio;
}

const map<double, double> &msalign::getKeyMsValueCharge() const {
    return keyMsValueCharge;
}

const map<std::string, std::string> &msalign::getMsMsgMap() const {
    return ms_msg_map;
}

const string &msalign::getId() const {
    return id;
}

const string &msalign::getScans() const {
    return scans;
}

const string &msalign::getRetentionTime() const {
    return retention_time;
}

const string &msalign::getActivation() const {
    return activation;
}

const string &msalign::getActivationSub() const {
    return activation_sub;
}

const string &msalign::getMsOneId() const {
    return ms_one_id;
}

const string &msalign::getMsOneScans() const {
    return ms_one_scans;
}

const string &msalign::getPrecursorMz() const {
    return precursor_mz;
}

const string &msalign::getPrecursorCharge() const {
    return precursor_charge;
}

const string &msalign::getPrecursorMassStr() const {
    return precursor_mass_str;
}

const string &msalign::getPrecursorIntensity() const {
    return precursor_intensity;
}

const vector<protein_ptr> &msalign::getCandidateProteins() const {
    return candidate_proteins;
}

double msalign::getPrsmBestScore() const {
    return prsm_best_score;
}

void msalign::setPrsmBestScore(double prsmBestScore) {
    prsm_best_score = prsmBestScore;
}

int msalign::getPrsmMatchPeaks() const {
    return prsm_match_peaks;
}

const string &msalign::getPrsmMatchProteinName() const {
    return prsm_match_protein_name;
}

string &msalign::getPrsmMatchProteinSeq()  {
    return prsm_match_protein_seq;
}

int msalign::getPrsmCutN() const {
    return prsm_cut_n;
}

int msalign::getPrsmCutC() const {
    return prsm_cut_c;
}

const vector<modi> &msalign::getPrsmPtm() const {
    return prsm_ptm;
}

double msalign::getPrsmMutX() const {
    return prsm_mut_x;
}

int msalign::getPrsmIonsPairNum() const {
    return prsm_ionsPair_num;
}

int msalign::getPrsmComplementsIonsNum() const {
    return prsm_complementsIons_num;
}

double msalign::getPrsmSeqScore() const {
    return prsm_seqScore;
}

int msalign::getPrsmSeqScoreCount() const {
    return prsm_seqScore_count;
}

int msalign::getPrsmMatchFragmentIonsNum() const {
    return prsm_match_fragmentIons_num;
}

void msalign::setPrsmMatchPeaks(int prsmMatchPeaks) {
    prsm_match_peaks = prsmMatchPeaks;
}

void msalign::setPrsmMatchProteinName(const string &prsmMatchProteinName) {
    prsm_match_protein_name = prsmMatchProteinName;
}

void msalign::setPrsmMatchProteinSeq(const string &prsmMatchProteinSeq) {
    prsm_match_protein_seq = prsmMatchProteinSeq;
}

void msalign::setPrsmCutN(int prsmCutN) {
    prsm_cut_n = prsmCutN;
}

void msalign::setPrsmCutC(int prsmCutC) {
    prsm_cut_c = prsmCutC;
}

void msalign::setPrsmPtm(const vector<modi> &prsmPtm) {
    prsm_ptm = prsmPtm;
}

void msalign::setPrsmMutX(double prsmMutX) {
    prsm_mut_x = prsmMutX;
}

void msalign::setPrsmIonsPairNum(int prsmIonsPairNum) {
    prsm_ionsPair_num = prsmIonsPairNum;
}

void msalign::setPrsmComplementsIonsNum(int prsmComplementsIonsNum) {
    prsm_complementsIons_num = prsmComplementsIonsNum;
}

void msalign::setPrsmSeqScore(double prsmSeqScore) {
    prsm_seqScore = prsmSeqScore;
}

void msalign::setPrsmSeqScoreCount(int prsmSeqScoreCount) {
    prsm_seqScore_count = prsmSeqScoreCount;
}

void msalign::setPrsmMatchFragmentIonsNum(int prsmMatchFragmentIonsNum) {
    prsm_match_fragmentIons_num = prsmMatchFragmentIonsNum;
}

double msalign::getPrsmPtmsMass() const {
    return prsm_ptms_mass;
}

const map<string, string> &msalign::getPrsmUnknownPtmMassMap() const {
    return prsm_unknown_ptmMass_map;
}

prsm::modify_ptr msalign::getPrsmUnknownModPtr() {
    return prsm_unknownModPtr;
}

double msalign::getPrsmTheoryProteoformMass() const {
    return prsm_theory_proteoform_mass;
}

void msalign::setPrsmPtmsMass(double prsmPtmsMass) {
    prsm_ptms_mass = prsmPtmsMass;
}

void msalign::setPrsmUnknownPtmMassMap(const map<string, string> &prsmUnknownPtmMassMap) {
    prsm_unknown_ptmMass_map = prsmUnknownPtmMassMap;
}

void msalign::setPrsmUnknownModPtr(prsm::modify_ptr prsmUnknownModPtr) {
    prsm_unknownModPtr = prsmUnknownModPtr;
}

void msalign::setPrsmTheoryProteoformMass(double prsmTheoryProteoformMass) {
    prsm_theory_proteoform_mass = prsmTheoryProteoformMass;
}
double msalign::getFeature1() const {
    return feature1;
}
void msalign::setFeature1(double feature1) {
    msalign::feature1 = feature1;
}
double msalign::getFeature2() const {
    return feature2;
}
void msalign::setFeature2(double feature2) {
    msalign::feature2 = feature2;
}
double msalign::getFeature3() const {
    return feature3;
}
void msalign::setFeature3(double feature3) {
    msalign::feature3 = feature3;
}
double msalign::getFeature4() const {
    return feature4;
}
void msalign::setFeature4(double feature4) {
    msalign::feature4 = feature4;
}
double msalign::getFeature5() const {
    return feature5;
}
void msalign::setFeature5(double feature5) {
    msalign::feature5 = feature5;
}
double msalign::getFeature6() const {
    return feature6;
}
void msalign::setFeature6(double feature6) {
    msalign::feature6 = feature6;
}
double msalign::getFeature7() const {
    return feature7;
}
void msalign::setFeature7(double feature7) {
    msalign::feature7 = feature7;
}
double msalign::getFeature8() const {
    return feature8;
}
void msalign::setFeature8(double feature8) {
    msalign::feature8 = feature8;
}
double msalign::getFeature9() const {
    return feature9;
}
void msalign::setFeature9(double feature9) {
    msalign::feature9 = feature9;
}
double msalign::getFeature10() const {
    return feature10;
}
void msalign::setFeature10(double feature10) {
    msalign::feature10 = feature10;
}
double msalign::getFeature11() const {
    return feature11;
}
void msalign::setFeature11(double feature11) {
    msalign::feature11 = feature11;
}
double msalign::getFeature12() const {
    return feature12;
}
void msalign::setFeature12(double feature12) {
    msalign::feature12 = feature12;
}
double msalign::getFeature13() const {
    return feature13;
}
void msalign::setFeature13(double feature13) {
    msalign::feature13 = feature13;
}
double msalign::getFeature14() const {
    return feature14;
}
void msalign::setFeature14(double feature14) {
    msalign::feature14 = feature14;
}
double msalign::getFeature15() const {
    return feature15;
}
void msalign::setFeature15(double feature15) {
    msalign::feature15 = feature15;
}
double msalign::getFeature16() const {
    return feature16;
}
void msalign::setFeature16(double feature16) {
    msalign::feature16 = feature16;
}
double msalign::getFeature17() const {
    return feature17;
}
void msalign::setFeature17(double feature17) {
    msalign::feature17 = feature17;
}
double msalign::getFeature18() const {
    return feature18;
}
void msalign::setFeature18(double feature18) {
    msalign::feature18 = feature18;
}
double msalign::getFeature19() const {
    return feature19;
}
void msalign::setFeature19(double feature19) {
    msalign::feature19 = feature19;
}
double msalign::getFeature20() const {
    return feature20;
}
void msalign::setFeature20(double feature20) {
    msalign::feature20 = feature20;
}
double msalign::getFeature21() const {
    return feature21;
}
void msalign::setFeature21(double feature21) {
    msalign::feature21 = feature21;
}
double msalign::getFeature22() const {
    return feature22;
}
void msalign::setFeature22(double feature22) {
    msalign::feature22 = feature22;
}
double msalign::getFeature23() const {
    return feature23;
}
void msalign::setFeature23(double feature23) {
    msalign::feature23 = feature23;
}
double msalign::getFeature24() const {
    return feature24;
}
void msalign::setFeature24(double feature24) {
    msalign::feature24 = feature24;
}
double msalign::getFeature25() const {
    return feature25;
}
void msalign::setFeature25(double feature25) {
    msalign::feature25 = feature25;
}
double msalign::getFeature26() const {
    return feature26;
}
void msalign::setFeature26(double feature26) {
    msalign::feature26 = feature26;
}
double msalign::getFeature27() const {
    return feature27;
}
void msalign::setFeature27(double feature27) {
    msalign::feature27 = feature27;
}
double msalign::getFeature28() const {
    return feature28;
}
void msalign::setFeature28(double feature28) {
    msalign::feature28 = feature28;
}
double msalign::getFeature29() const {
    return feature29;
}
void msalign::setFeature29(double feature29) {
    msalign::feature29 = feature29;
}
double msalign::getFeature30() const {
    return feature30;
}
void msalign::setFeature30(double feature30) {
    msalign::feature30 = feature30;
}
double msalign::getFeature31() const {
    return feature31;
}
void msalign::setFeature31(double feature31) {
    msalign::feature31 = feature31;
}
double msalign::getFeature32() const {
    return feature32;
}
void msalign::setFeature32(double feature32) {
    msalign::feature32 = feature32;
}
double msalign::getFeature33() const {
    return feature33;
}
void msalign::setFeature33(double feature33) {
    msalign::feature33 = feature33;
}
double msalign::getFeature34() const {
    return feature34;
}
void msalign::setFeature34(double feature34) {
    msalign::feature34 = feature34;
}
double msalign::getFeature35() const {
    return feature35;
}
void msalign::setFeature35(double feature35) {
    msalign::feature35 = feature35;
}
double msalign::getFeature36() const {
    return feature36;
}
void msalign::setFeature36(double feature36) {
    msalign::feature36 = feature36;
}
double msalign::getFeature37() const {
    return feature37;
}
void msalign::setFeature37(double feature37) {
    msalign::feature37 = feature37;
}
double msalign::getFeature38() const {
    return feature38;
}
void msalign::setFeature38(double feature38) {
    msalign::feature38 = feature38;
}
double msalign::getFeature39() const {
    return feature39;
}
void msalign::setFeature39(double feature39) {
    msalign::feature39 = feature39;
}
double msalign::getFeature40() const {
    return feature40;
}
void msalign::setFeature40(double feature40) {
    msalign::feature40 = feature40;
}
double msalign::getFeature41() const {
    return feature41;
}
void msalign::setFeature41(double feature41) {
    msalign::feature41 = feature41;
}
double msalign::getFeature42() const {
    return feature42;
}
void msalign::setFeature42(double feature42) {
    msalign::feature42 = feature42;
}
double msalign::getFeature43() const {
    return feature43;
}
void msalign::setFeature43(double feature43) {
    msalign::feature43 = feature43;
}
double msalign::getFeature44() const {
    return feature44;
}
void msalign::setFeature44(double feature44) {
    msalign::feature44 = feature44;
}
double msalign::getFeature45() const {
    return feature45;
}
void msalign::setFeature45(double feature45) {
    msalign::feature45 = feature45;
}
double msalign::getFeature46() const {
    return feature46;
}
void msalign::setFeature46(double feature46) {
    msalign::feature46 = feature46;
}
double msalign::getFeature47() const {
    return feature47;
}
void msalign::setFeature47(double feature47) {
    msalign::feature47 = feature47;
}
double msalign::getFeature48() const {
    return feature48;
}
void msalign::setFeature48(double feature48) {
    msalign::feature48 = feature48;
}
double msalign::getFeature49() const {
    return feature49;
}
void msalign::setFeature49(double feature49) {
    msalign::feature49 = feature49;
}
double msalign::getFeature50() const {
    return feature50;
}
void msalign::setFeature50(double feature50) {
    msalign::feature50 = feature50;
}
