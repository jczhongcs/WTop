//
// Created by wenzhong on 2023/3/22.
//

#include <mutex>
#include "prsm_msg_write.hpp"
#include "tinyxml.h"
#include "prsm_store_ptm.h"
#include <sys/file.h>

#include <string>
#include <thread>

std::string prsm_msg_write::generateSubstring(const std::string& str) {
    size_t lastDotPos = str.rfind('.');
    if (lastDotPos != std::string::npos) {
        // 返回最后一个'.'之前的子串
        return str.substr(0, lastDotPos+1);
    }
    return str;  // 如果没有找到'.'，返回原始字符串
}


std::string prsm_msg_write::doubleToStringWithPrecision(double value, int precision) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}


void prsm_msg_write::AddPrsm(TiXmlElement* prsm_list,msalign_ptr msalignPtr,prsm::modify_ptr modptr){
    utils_string_util::processCutLocation(msalignPtr->getPrsmMatchProteinSeq(),
                                          msalignPtr->getPrsmCutN(), msalignPtr->getPrsmCutC());
    TiXmlElement* prsm = new TiXmlElement("prsm");
    prsm_list->LinkEndChild(prsm);
    //prsm基本信息
    prsm->LinkEndChild(new TiXmlElement("file_name"))->LinkEndChild(new TiXmlText("Wtop_Prsms"));
    prsm->LinkEndChild(new TiXmlElement("prsm_id"))->LinkEndChild(new TiXmlText("-1"));
    string spectrum_id = msalignPtr->getId().substr(3);
    prsm->LinkEndChild(new TiXmlElement("spectrum_id"))->LinkEndChild(new TiXmlText(spectrum_id.c_str()));
    string spectrum_scan = msalignPtr->getScans().substr(6);
    prsm->LinkEndChild(new TiXmlElement("spectrum_scan"))->LinkEndChild(new TiXmlText(spectrum_scan.c_str()));
    prsm->LinkEndChild(new TiXmlElement("precursor_id"))->LinkEndChild(new TiXmlText("0"));
    prsm->LinkEndChild(new TiXmlElement("sample_feature_id"))->LinkEndChild(new TiXmlText("-1"));
    prsm->LinkEndChild(new TiXmlElement("sample_feature_inte"))->LinkEndChild(new TiXmlText("-1.0"));
    prsm->LinkEndChild(new TiXmlElement("frac_feature_score"))->LinkEndChild(new TiXmlText("-1000.000"));
    prsm->LinkEndChild(new TiXmlElement("spectrum_number"))->LinkEndChild(new TiXmlText("1"));

    string ori_prec_mass = doubleToStringWithPrecision(msalignPtr->getPrecursorMass(), 4);
    prsm->LinkEndChild(new TiXmlElement("ori_prec_mass"))->LinkEndChild(new TiXmlText(ori_prec_mass.c_str()));
    string adjusted_prec_mass = doubleToStringWithPrecision(msalignPtr->getPrsmTheoryProteoformMass(),4);
    prsm->LinkEndChild(new TiXmlElement("adjusted_prec_mass"))->LinkEndChild(new TiXmlText(adjusted_prec_mass.c_str()));

    prsm->LinkEndChild(new TiXmlElement("fdr"))->LinkEndChild(new TiXmlText("-1.000000"));
    prsm->LinkEndChild(new TiXmlElement("proteoform_fdr"))->LinkEndChild(new TiXmlText("-1.00"));

    string match_peak_num = to_string(msalignPtr->getPrsmMatchPeaks());
    prsm->LinkEndChild(new TiXmlElement("match_peak_num"))->LinkEndChild(new TiXmlText(match_peak_num.c_str()));
    string match_fragment_num = to_string(msalignPtr->getPrsmMatchFragmentIonsNum());
    prsm->LinkEndChild(new TiXmlElement("match_fragment_num"))->LinkEndChild(new TiXmlText(match_fragment_num.c_str()));

    prsm->LinkEndChild(new TiXmlElement("norm_match_fragment_num"))->LinkEndChild(new TiXmlText("0"));
    prsm->LinkEndChild(new TiXmlElement("fraction_feature_time_apex"))->LinkEndChild(new TiXmlText("-1.00"));


    // Proteoform节点和嵌套
    TiXmlElement* proteoform = new TiXmlElement("proteoform");
    prsm->LinkEndChild(proteoform);

    // Fasta序列节点
    TiXmlElement* fastaSeq = new TiXmlElement("fasta_seq");
    proteoform->LinkEndChild(fastaSeq);
    istringstream iss(msalignPtr->getPrsmMatchProteinName().substr(1));
    string seq_name, seq_desc;
    // 以第一个空格为分隔符，将字符串分为两部分
    getline(iss, seq_name, ' ');
    getline(iss, seq_desc);
    fastaSeq->LinkEndChild(new TiXmlElement("seq_name"))->LinkEndChild(new TiXmlText(seq_name.c_str()));
    fastaSeq->LinkEndChild(new TiXmlElement("seq_desc"))->LinkEndChild(new TiXmlText(seq_desc.c_str()));

    //prot_mod节点
    TiXmlElement* prot_mod = new TiXmlElement("prot_mod");
    proteoform->LinkEndChild(prot_mod);
    prot_mod->LinkEndChild(new TiXmlElement("name"))->LinkEndChild(new TiXmlText("NONE"));


    string start_pos = to_string(msalignPtr->getPrsmCutN());
    proteoform->LinkEndChild(new TiXmlElement("start_pos"))->LinkEndChild(new TiXmlText(start_pos.c_str()));
    int end_pos_int = msalignPtr->getPrsmCutN() + msalignPtr->getPrsmMatchProteinSeq().size() -3;
    string end_pos = to_string(end_pos_int);
    proteoform->LinkEndChild(new TiXmlElement("end_pos"))->LinkEndChild(new TiXmlText(end_pos.c_str()));

    proteoform->LinkEndChild(new TiXmlElement("proteo_cluster_id"))->LinkEndChild(new TiXmlText("-1"));
    proteoform->LinkEndChild(new TiXmlElement("prot_id"))->LinkEndChild(new TiXmlText("-1"));



    StorePtmPtrVec storePtmPtrVec;
    map<char,string> charUniMap = {
        {'A',"34"},
        {'B',"36"},
        {'C',"1"},
        {'D',"37"},
        {'E',"21"},
        {'*',"-1"}
    };

    map<int,int> startEnd;  //修饰所在的序列区间
    vector<string> varPtmsContainer ;  //保存修饰
    map<string,string> ptmMassMap = msalignPtr->getPrsmUnknownPtmMassMap();
    for (int j = 0; j < msalignPtr->getPrsmPtm().size(); j++) {
        string ptms = modptr->analysis(msalignPtr->getPrsmPtm()[j].mod_mass);
        if (ptms != "0" && msalignPtr->getPrsmPtmsMass() != 0) {
          int leftBpPos = (int) msalignPtr->getPrsmPtm()[j].first + 1;
          int rightBpPos = (int) msalignPtr->getPrsmPtm()[j].second + 1;
          for(char x:ptms){
              if(x != '*'){//已知修饰
                  StorePtmPtr storePtm = make_shared<StorePtm>(atoi(charUniMap[x].c_str()),x,modptr->charPTMsMap[x]
                                                               ,leftBpPos,rightBpPos,modptr->charDtoSPtmsMap[x],false);
                  storePtmPtrVec.push_back(storePtm);
              }else{//未知修饰
                  StorePtmPtr storePtm = make_shared<StorePtm>(atoi(charUniMap[x].c_str()),x,"Unexpected"
                      ,leftBpPos,rightBpPos,ptmMassMap["Unknown_PTM_Mass"], true);
                  storePtmPtrVec.push_back(storePtm);
              }
          }


          int protein_len = 0;
          int dotcout = 0;
          for(auto x:msalignPtr->getPrsmMatchProteinSeq()){
              if(x == '.')
                  dotcout++;
              if(dotcout == 1)
                  protein_len++;
          }
          protein_len--;//减去第一个点的长度


          int ptm_left = msalignPtr->getPrsmPtm()[j].first;
          int ptm_right = msalignPtr->getPrsmPtm()[j].second;

          if(ptm_right > msalignPtr->getPrsmMatchProteinSeq().size() -1){
              ptm_right = msalignPtr->getPrsmMatchProteinSeq().size() -1;
          }
          if(ptm_left == ptm_right){
              ptm_left--;
          }
          if(ptm_left < -1){
                ptm_left = -1;
          }

          varPtmsContainer.push_back(ptms) ;
          startEnd.insert(make_pair(ptm_left, ptm_right));
        }
    }

    int variable_ptm_num = 0,unexpected_ptm_num = 0;
    for(auto sPtm:storePtmPtrVec){
        if(sPtm->isUnknowPtm){
            unexpected_ptm_num++;
        }else{
            variable_ptm_num++;
        }
    }

    proteoform->LinkEndChild(new TiXmlElement("variable_ptm_num"))->LinkEndChild(new TiXmlText(to_string(variable_ptm_num).c_str()));
    proteoform->LinkEndChild(new TiXmlElement("unexpected_ptm_num"))->LinkEndChild(new TiXmlText(to_string(unexpected_ptm_num).c_str()));



    string proteo_match_seq = "";
    utils_string_util::processProteinFomUnkMass(msalignPtr->getPrsmMatchProteinSeq(), varPtmsContainer,
                                         modptr->charPTMsMap, startEnd, proteo_match_seq,ptmMassMap["Unknown_PTM_Mass"]);

    string cuted_str = generateSubstring(proteo_match_seq);

    proteoform->LinkEndChild(new TiXmlElement("proteo_match_seq"))->LinkEndChild(new TiXmlText(cuted_str.c_str()));


    // Mass Shift节点和嵌套
    TiXmlElement* massShiftList = new TiXmlElement("mass_shift_list");
    proteoform->LinkEndChild(massShiftList);

    map<char,string> amoCharNameMap = {
        {'X',"None"},
        {'A',"Alanine"},
        {'R',"Arginine"},
        {'N',"Asparagine"},
        {'D',"Aspartic_acid"},
        {'C',"Cysteine"},
        {'E',"Glutamic_acid"},
        {'Q',"Glutamine"},
        {'G',"Glycine"},
        {'H',"Histidine"},
        {'I',"Isoleucine"},
        {'L',"Leucine"},
        {'K',"Lysine"},
        {'M',"Methionine"},
        {'F',"Phenylalanine"},
        {'P',"Proline"},
        {'S',"Serine"},
        {'T',"Threonine"},
        {'W',"Tryptophan"},
        {'Y',"Tyrosine"},
        {'V',"Valine"},
        {'U',"Selenocysteine"},
        {'O',"Ornithine"},
        {'C',"Cysteine"}
    };

    for (int j = 0; j < msalignPtr->getPrsmPtm().size(); j++) {
        string ptms = modptr->analysis(msalignPtr->getPrsmPtm()[j].mod_mass);
        if (ptms != "0" && msalignPtr->getPrsmPtmsMass() != 0) {
            int leftBpPos = (int) msalignPtr->getPrsmPtm()[j].first + 1;
            int rightBpPos = (int) msalignPtr->getPrsmPtm()[j].second + 1;
            if(leftBpPos == rightBpPos)
                leftBpPos--;
            if(leftBpPos < 0 )
                leftBpPos = 0;

            TiXmlElement* massShift = new TiXmlElement("mass_shift");
            massShiftList->LinkEndChild(massShift);
            massShift->LinkEndChild(new TiXmlElement("left_bp_pos"))->LinkEndChild(new TiXmlText(to_string(leftBpPos).c_str()));
            massShift->LinkEndChild(new TiXmlElement("right_bp_pos"))->LinkEndChild(new TiXmlText(to_string(rightBpPos).c_str()));
            double shift = 0;
            for(auto single_ptm:ptms){
                if (single_ptm != '*') {
                    shift += atof(modptr->charDtoSPtmsMap[single_ptm].c_str());
                } else{
                    shift += stold(ptmMassMap["Unknown_PTM_Mass"]);
                }
            }
            massShift->LinkEndChild(new TiXmlElement("shift"))->LinkEndChild(new TiXmlText(to_string(shift).c_str()));


            TiXmlElement* alteration_list = new TiXmlElement("alteration_list");
            massShift->LinkEndChild(alteration_list);
            //区间相同的修饰为同一个mass_shift,区间相同不同类型的修饰放alteration_list下
            for(char x:ptms){
                TiXmlElement* alteration = new TiXmlElement("alteration");
                alteration_list->LinkEndChild(alteration);

                if(x != '*'){//已知修饰
                    alteration->LinkEndChild(new TiXmlElement("left_bp_pos"))->LinkEndChild(new TiXmlText(to_string(leftBpPos).c_str()));
                    alteration->LinkEndChild(new TiXmlElement("right_bp_pos"))->LinkEndChild(new TiXmlText(to_string(leftBpPos+1).c_str()));
                    TiXmlElement* alter_type = new TiXmlElement("alter_type");
                    alteration->LinkEndChild(alter_type);
                    alter_type->LinkEndChild(new TiXmlElement("name"))->LinkEndChild(new TiXmlText("Variable"));
                    alteration->LinkEndChild(new TiXmlElement("mass"))->LinkEndChild(new TiXmlText(modptr->charDtoSPtmsMap[x].c_str()));

                    TiXmlElement* mod = new TiXmlElement("mod");
                    alteration->LinkEndChild(mod);
                    TiXmlElement* ori_residue = new TiXmlElement("ori_residue");
                    mod->LinkEndChild(ori_residue);
                    TiXmlElement* amino_acid = new TiXmlElement("amino_acid");
                    ori_residue->LinkEndChild(amino_acid);
                    char alter_amino = msalignPtr->getPrsmMatchProteinSeq()[leftBpPos+1];
                    amino_acid->LinkEndChild(new TiXmlElement("name"))->LinkEndChild(new TiXmlText(amoCharNameMap[alter_amino].c_str()));
                    TiXmlElement* ptm = new TiXmlElement("ptm");
                    ori_residue->LinkEndChild(ptm);
                    ptm->LinkEndChild(new TiXmlElement("abbreviation"))->LinkEndChild(new TiXmlText("No PTM"));
                    ptm->LinkEndChild(new TiXmlElement("unimod"))->LinkEndChild(new TiXmlText("-1"));

                    TiXmlElement* mod_residue = new TiXmlElement("mod_residue");
                    mod->LinkEndChild(mod_residue);
                    TiXmlElement* amino_acid1 = new TiXmlElement("amino_acid");
                    mod_residue->LinkEndChild(amino_acid1);
                    amino_acid1->LinkEndChild(new TiXmlElement("name"))->LinkEndChild(new TiXmlText(amoCharNameMap[alter_amino].c_str()));
                    TiXmlElement* ptm1 = new TiXmlElement("ptm");
                    mod_residue->LinkEndChild(ptm1);
                    ptm1->LinkEndChild(new TiXmlElement("abbreviation"))->LinkEndChild(new TiXmlText(modptr->charPTMsMap[x].c_str()));
                    ptm1->LinkEndChild(new TiXmlElement("unimod"))->LinkEndChild(new TiXmlText(charUniMap[x].c_str()));
                }else{//未知修饰
                    alteration->LinkEndChild(new TiXmlElement("left_bp_pos"))->LinkEndChild(new TiXmlText(to_string(leftBpPos).c_str()));
                    alteration->LinkEndChild(new TiXmlElement("right_bp_pos"))->LinkEndChild(new TiXmlText(to_string(rightBpPos).c_str()));
                    TiXmlElement* alter_type = new TiXmlElement("alter_type");
                    alteration->LinkEndChild(alter_type);
                    alter_type->LinkEndChild(new TiXmlElement("name"))->LinkEndChild(new TiXmlText("Unexpected"));
                    alteration->LinkEndChild(new TiXmlElement("mass"))->LinkEndChild(new TiXmlText(ptmMassMap["Unknown_PTM_Mass"].c_str()));

                    TiXmlElement* mod = new TiXmlElement("mod");
                    alteration->LinkEndChild(mod);
                    TiXmlElement* ori_residue = new TiXmlElement("ori_residue");
                    mod->LinkEndChild(ori_residue);
                    TiXmlElement* amino_acid = new TiXmlElement("amino_acid");
                    ori_residue->LinkEndChild(amino_acid);
                    char alter_amino = msalignPtr->getPrsmMatchProteinSeq()[leftBpPos+1];
                    amino_acid->LinkEndChild(new TiXmlElement("name"))->LinkEndChild(new TiXmlText("None"));
                    TiXmlElement* ptm = new TiXmlElement("ptm");
                    ori_residue->LinkEndChild(ptm);
                    ptm->LinkEndChild(new TiXmlElement("abbreviation"))->LinkEndChild(new TiXmlText("No PTM"));
                    ptm->LinkEndChild(new TiXmlElement("unimod"))->LinkEndChild(new TiXmlText("-1"));

                    TiXmlElement* mod_residue = new TiXmlElement("mod_residue");
                    mod->LinkEndChild(mod_residue);
                    TiXmlElement* amino_acid1 = new TiXmlElement("amino_acid");
                    mod_residue->LinkEndChild(amino_acid1);
                    amino_acid1->LinkEndChild(new TiXmlElement("name"))->LinkEndChild(new TiXmlText("None"));
                    TiXmlElement* ptm1 = new TiXmlElement("ptm");
                    mod_residue->LinkEndChild(ptm1);
                    ptm1->LinkEndChild(new TiXmlElement("abbreviation"))->LinkEndChild(new TiXmlText("No PTM"));
                    ptm1->LinkEndChild(new TiXmlElement("unimod"))->LinkEndChild(new TiXmlText("-1"));
                }
            }
        }
    }

    prsm->LinkEndChild(new TiXmlElement("OriProLen"))->LinkEndChild(new TiXmlText(to_string( msalignPtr->getPrsmMatchProteinSeq().size()).c_str()));
    ///特征输出到xml文件
    prsm->LinkEndChild(new TiXmlElement("feature1"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature1()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature2"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature2()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature3"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature3()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature4"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature4()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature5"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature5()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature6"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature6()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature7"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature7()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature8"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature8()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature9"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature9()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature10"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature10()).c_str()));

    prsm->LinkEndChild(new TiXmlElement("feature11"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature11()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature12"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature12()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature13"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature13()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature14"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature14()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature15"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature15()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature16"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature16()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature17"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature17()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature18"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature18()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature19"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature19()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature20"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature20()).c_str()));

    prsm->LinkEndChild(new TiXmlElement("feature21"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature21()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature22"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature22()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature23"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature23()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature24"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature24()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature25"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature25()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature26"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature26()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature27"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature27()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature28"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature28()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature29"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature29()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature30"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature30()).c_str()));

    prsm->LinkEndChild(new TiXmlElement("feature31"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature31()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature32"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature32()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature33"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature33()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature34"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature34()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature35"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature35()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature36"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature36()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature37"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature37()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature38"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature38()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature39"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature39()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature40"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature40()).c_str()));

    prsm->LinkEndChild(new TiXmlElement("feature41"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature41()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature42"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature42()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature43"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature43()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature44"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature44()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature45"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature45()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature46"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature46()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature47"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature47()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature48"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature48()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature49"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature49()).c_str()));
    prsm->LinkEndChild(new TiXmlElement("feature50"))->LinkEndChild(new TiXmlText(to_string(msalignPtr->getFeature50()).c_str()));


}

void prsm_msg_write::OutXmlToEvalue( msalignPtrVec msalignptrVec,
                                    string xmlFileName){
//    std::ostringstream oss;
//    oss << std::this_thread::get_id();
//    std::string stringID = oss.str();
//    xmlFileName = xmlFileName + "_" + stringID;
    string ms_id = msalignptrVec[0]->getId();
    ms_id = ms_id.substr(3);
    xmlFileName = xmlFileName + "_" + ms_id;

    FILE* file = fopen(xmlFileName.c_str(), "rb");
    bool xmlExist;
    if (file) {
        xmlExist = true;
        fclose(file);
    }else{
        xmlExist = false;
    }
    if(xmlExist){
        //xml文件已存在，添加新prsm结点
        TiXmlDocument doc(xmlFileName.c_str());
        doc.LoadFile();
        TiXmlElement* prsm_list = doc.FirstChildElement("prsm_list");
        if (prsm_list) {
            // prsm_list node exists, add new child nodes here
            for(auto msalignPtr:msalignptrVec){
                if (msalignPtr->getPrsmMatchProteinName().size() > 0) {
                    prsm::modify_ptr modptr = msalignPtr->getPrsmUnknownModPtr();
                    AddPrsm(prsm_list, msalignPtr, modptr);
                }
            }
            // Save the modified XML document
            doc.SaveFile(xmlFileName.c_str());
        }else{
            cout<<endl<<"xml error"<<endl;
            exit(718);
        }
    }
    else{
        //新建xml文件
        TiXmlDocument doc;
        TiXmlDeclaration* declaration = new TiXmlDeclaration("1.0", "UTF-8", "");
        doc.LinkEndChild(declaration);
        TiXmlElement* prsm_list = new TiXmlElement("prsm_list");
        doc.LinkEndChild(prsm_list);
        //添加prsm节点
        for(auto msalignPtr:msalignptrVec){
            if (msalignPtr->getPrsmMatchProteinName().size() > 0) {
                prsm::modify_ptr modptr = msalignPtr->getPrsmUnknownModPtr();
                AddPrsm(prsm_list, msalignPtr, modptr);
            }
        }
        doc.SaveFile(xmlFileName.c_str());
    }



}


void prsm_msg_write::writeResultNews(msalign_ptr msalignPtr,
                                     string Result_file,
                                     prsm::modify_ptr modptr)
{

    if (msalignPtr->getPrsmMatchProteinName().size() > 0)
    {
        std::ofstream out(Result_file, ios::out | ios::app);
        double precursorMassDouble = msalignPtr->getPrecursorMass(); //atof(msIt->Pre_mass.substr(msIt->Pre_mass.find('=') + 1, msIt->Pre_mass.size()).c_str());
        out << "#BEGIN" << endl;
        out << "#" << msalignPtr->getId()
            << " , #" << msalignPtr->getScans()
            << ", #Types : " << msalignPtr->getActivationSub()
            << ", #Protero : " << msalignPtr->getPrsmMatchProteinName() << "." << endl;
        out << "Match Peaks = " << msalignPtr->getPrsmMatchPeaks()
            << ",score = " << msalignPtr->getPrsmBestScore() << endl;
        out << "CutLocationN = " << msalignPtr->getPrsmCutN()
            << " ,CutLocationC = " << msalignPtr->getPrsmCutC() << endl;
        out << "Precursor mass = " << setiosflags(ios::fixed) << setprecision(4)
            << msalignPtr->getPrecursorMass()
            << " ,Proteoform mass = " << msalignPtr->getPrsmTheoryProteoformMass() << endl;
        out << "PPM = " << (precursorMassDouble - msalignPtr->getPrsmTheoryProteoformMass()) / precursorMassDouble * 1000000 << endl;
        out << "Protein length = " << msalignPtr->getPrsmMatchProteinSeq().size()<< endl;
        //处理截断
        utils_string_util::processCutLocation(msalignPtr->getPrsmMatchProteinSeq(),
                                              msalignPtr->getPrsmCutN(), msalignPtr->getPrsmCutC());

        //输出修饰
        out << "#Mod BEGIN" << endl;
        double pMass = 0;
        map<int,int> startEnd;  //修饰所在的序列区间
        vector<string> varPtmsContainer ;  //保存修饰
        for (int j = 0; j < msalignPtr->getPrsmPtm().size(); j++) {
            string ptms = modptr->analysis(msalignPtr->getPrsmPtm()[j].mod_mass);
            if (ptms != "0" && msalignPtr->getPrsmPtmsMass() != 0) {
                out << "start: " << (int) msalignPtr->getPrsmPtm()[j].first << "-" << (int) msalignPtr->getPrsmPtm()[j].second
                    << ";mass: " << msalignPtr->getPrsmPtm()[j].mod_mass
                    << " ; Mod# " << ptms << endl;
                pMass += msalignPtr->getPrsmPtm()[j].mod_mass;
                varPtmsContainer.push_back(ptms) ;
                startEnd.insert(make_pair(msalignPtr->getPrsmPtm()[j].first, msalignPtr->getPrsmPtm()[j].second));
            }
        }
        out << "#Mod END" << endl;
        for (auto it = msalignPtr->getPrsmUnknownPtmMassMap().begin(); it !=  msalignPtr->getPrsmUnknownPtmMassMap().end(); ++it)
        {
            out << it->first <<" = "<<it->second<<endl;
        }
        string pfSeq = "";  //蛋白质形序列
//        for(auto it = modptr->charPTMsMap.begin(); it != modptr->charPTMsMap.end();++it)
//        {
//            cout <<it->first <<" "<< it->second <<endl;
//        }
        //处理蛋白质形
        utils_string_util::processProteinFom(msalignPtr->getPrsmMatchProteinSeq(), varPtmsContainer,
                                             modptr->charPTMsMap, startEnd, pfSeq);
        out << "RawProteinSeq# " << msalignPtr->getPrsmMatchProteinSeq() << endl;
        out <<"ProteoForm# "<< pfSeq <<endl;
        //其余数据输出
        out << "Spact ModifyMass = " << msalignPtr->getPrsmPtmsMass() << "\nTrue ModifyMass = " << pMass << endl;
        out << "All peaks = " << msalignPtr->getIonsMassContainer().size() << endl;
        out << "match peaks = " << (int) msalignPtr->getPrsmMatchPeaks() << endl;
        out << "Ms complementIons =  " << msalignPtr->getCompIonsIndexMap().size() << endl;
        out << "matched fragment ions = " << msalignPtr->getPrsmMatchFragmentIonsNum() << endl;
        out << "mut_x = " << msalignPtr->getPrsmMutX() << endl;
        out << "ionsPair = " << msalignPtr->getPrsmIonsPairNum() << endl;
        out << "huBuIonsCount = " << (int) msalignPtr->getPrsmComplementsIonsNum() - 1 << endl;
        out << "seqIonsCount = " << msalignPtr->getPrsmSeqScoreCount() << endl;
        out << "seqScore = " << msalignPtr->getPrsmSeqScore() << endl;

        ///调试打分
        out << "feature01 = " << msalignPtr->getFeature1() << endl;
        out << "feature02 = " << msalignPtr->getFeature2() << endl;
        out << "feature03 = " << msalignPtr->getFeature3() << endl;
        out << "feature04 = " << msalignPtr->getFeature4() << endl;
        out << "feature05 = " << msalignPtr->getFeature5() << endl;
        out << "feature06 = " << msalignPtr->getFeature6() << endl;
        out << "feature07 = " << msalignPtr->getFeature7() << endl;
        out << "feature08 = " << msalignPtr->getFeature8() << endl;
        out << "feature09 = " << msalignPtr->getFeature9() << endl;
        out << "feature10 = " << msalignPtr->getFeature10() << endl;
        out << "feature11 = " << msalignPtr->getFeature11() << endl;
        out << "feature12 = " << msalignPtr->getFeature12() << endl;
        out << "feature13 = " << msalignPtr->getFeature13() << endl;
        out << "feature14 = " << msalignPtr->getFeature14() << endl;
        out << "feature15 = " << msalignPtr->getFeature15() << endl;
        out << "feature16 = " << msalignPtr->getFeature16() << endl;
        out << "feature17 = " << msalignPtr->getFeature17() << endl;
        out << "feature18 = " << msalignPtr->getFeature18() << endl;
        out << "feature19 = " << msalignPtr->getFeature19() << endl;
        out << "feature20 = " << msalignPtr->getFeature20() << endl;
        out << "feature21 = " << msalignPtr->getFeature21() << endl;
        out << "feature22 = " << msalignPtr->getFeature22() << endl;
        out << "feature23 = " << msalignPtr->getFeature23() << endl;
        out << "feature24 = " << msalignPtr->getFeature24() << endl;
        out << "feature25 = " << msalignPtr->getFeature25() << endl;
        out << "feature26 = " << msalignPtr->getFeature26() << endl;
        out << "feature27 = " << msalignPtr->getFeature27() << endl;
        out << "feature28 = " << msalignPtr->getFeature28() << endl;
        out << "feature29 = " << msalignPtr->getFeature29() << endl;
        out << "feature30 = " << msalignPtr->getFeature30() << endl;
        out << "feature31 = " << msalignPtr->getFeature31() << endl;
        out << "feature32 = " << msalignPtr->getFeature32() << endl;
        out << "feature33 = " << msalignPtr->getFeature33() << endl;
        out << "feature34 = " << msalignPtr->getFeature34() << endl;
        out << "feature35 = " << msalignPtr->getFeature35() << endl;
        out << "feature36 = " << msalignPtr->getFeature36() << endl;
        out << "feature37 = " << msalignPtr->getFeature37() << endl;
        out << "feature38 = " << msalignPtr->getFeature38() << endl;
        out << "feature39 = " << msalignPtr->getFeature39() << endl;
        out << "feature40 = " << msalignPtr->getFeature40() << endl;
        out << "feature41 = " << msalignPtr->getFeature41() << endl;
        out << "feature42 = " << msalignPtr->getFeature42() << endl;
        out << "feature43 = " << msalignPtr->getFeature43() << endl;
        out << "feature44 = " << msalignPtr->getFeature44() << endl;
        out << "feature45 = " << msalignPtr->getFeature45() << endl;
        out << "feature46 = " << msalignPtr->getFeature46() << endl;
        out << "feature47 = " << msalignPtr->getFeature47() << endl;
        out << "feature48 = " << msalignPtr->getFeature48() << endl;
        out << "feature49 = " << msalignPtr->getFeature49() << endl;
        out << "feature50 = " << msalignPtr->getFeature50() << endl;

        out << "#END\n" << endl;
        out.close();
    }   //输出答案结束

}

