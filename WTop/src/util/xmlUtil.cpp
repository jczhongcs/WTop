//
// Created by LYC on 2025/3/29.
//

#include "xmlUtil.h"

using namespace std;

struct PrsmNode {
    int id;
    TiXmlElement* element;

    PrsmNode(int _id, TiXmlElement* _element) : id(_id), element(_element) {}
};

// 比较函数用于排序
bool compareById(const PrsmNode& a, const PrsmNode& b) {
    return a.id < b.id;
}


void sortPrsm(string combinedPrsmFile){

    string filename = combinedPrsmFile;
    TiXmlDocument doc(filename.c_str());
    bool loadOkay = doc.LoadFile();
    string out_filename = "out_spid_evalue.csv";

    if (loadOkay) {
        std::cout << "XML 文件加载成功" << std::endl;

        // 获取根元素
        TiXmlElement* prsm_list = doc.FirstChildElement();
        if (prsm_list) {
            std::vector<PrsmNode> prsmNodes;
            for (TiXmlElement* prsm = prsm_list->FirstChildElement("prsm"); prsm != NULL; prsm = prsm->NextSiblingElement("prsm"))
            {
                // 获取 evalue 节点的文本值
                TiXmlElement* spectrum_id = prsm->FirstChildElement("spectrum_id");
                int id = stoi(spectrum_id->GetText());
                prsmNodes.emplace_back(id, prsm);
            }
            sort(prsmNodes.begin(), prsmNodes.end(), compareById);



            TiXmlDocument doc;
            TiXmlDeclaration* declaration = new TiXmlDeclaration("1.0", "UTF-8", "");
            doc.LinkEndChild(declaration);
            TiXmlElement* new_prsm_list = new TiXmlElement("prsm_list");
            doc.LinkEndChild(new_prsm_list);
            for(int i=0;i<prsmNodes.size();i++){
                TiXmlElement* prsmCopy = new TiXmlElement(*prsmNodes[i].element);
                new_prsm_list->LinkEndChild(prsmCopy);
            }


            //添加prsm节点
            doc.SaveFile("sorted_prsm.wtop_combined");




        } else {
            std::cerr << "未找到 prsm_list 节点" << std::endl;
        }
    } else {
        std::cerr << "无法加载 XML 文件" << std::endl;
    }

}

void evalueToCsv(string evalueFile){

    string filename = evalueFile;
    TiXmlDocument doc(filename.c_str());
    bool loadOkay = doc.LoadFile();
    string out_filename = "out_spid_evalue.csv";
    std::ofstream outputFile;
    outputFile.open(out_filename,ios::out|ios::trunc);
    if (loadOkay) {
        std::cout << "XML 文件加载成功" << std::endl;

        // 获取根元素
        TiXmlElement* prsm_list = doc.FirstChildElement();
        if (prsm_list) {
            outputFile <<"Data file name"<<","<<"Prsm ID"<<","<<"Spectrum ID"<<","<<"Fragmentation"<<","<<"Scan(s)"<<","
                       <<"Retention time"<<","<<"#peaks"<<","<<"Charge"<<","<<"Precursor mass"<<","<<"Adjusted precursor mass"<<","
                       <<"Proteoform ID"<<","<<"Feature intensity"<<","<<"Feature score"<<","<<"Feature apex time"<<","<<"#Protein hits"<<","
                       <<"Protein accession"<<","<<"Protein description"<<","<<"First residue"<<","<<"Last residue"<<","<<"Special amino acids"<<","
                       <<"Database protein sequence"<<","<<"Proteoform"<<","<<"Proteoform mass"<<","<<"evalue"<<","<<"match peak"<<","<<"match fragment"<<std::endl;

            for (TiXmlElement* prsm = prsm_list->FirstChildElement("prsm"); prsm != NULL; prsm = prsm->NextSiblingElement("prsm")) {
                // 获取 evalue 节点的文本值
                TiXmlElement* spectrum_id = prsm->FirstChildElement("spectrum_id");
                TiXmlElement* spectrum_scan = prsm->FirstChildElement("spectrum_scan");
                TiXmlElement* ori_prec_mass = prsm->FirstChildElement("ori_prec_mass");
                TiXmlElement* match_peak_num = prsm->FirstChildElement("match_peak_num");
                TiXmlElement* match_fragment_num = prsm->FirstChildElement("match_fragment_num");

                TiXmlElement* proteoform =  prsm->FirstChildElement("proteoform");
                TiXmlElement* fasta_seq =  proteoform->FirstChildElement("fasta_seq");
                TiXmlElement* start_pos =  proteoform->FirstChildElement("start_pos");
                TiXmlElement* end_pos =  proteoform->FirstChildElement("end_pos");
                TiXmlElement* seq_name =  fasta_seq->FirstChildElement("seq_name");
                TiXmlElement* seq_desc =  fasta_seq->FirstChildElement("seq_desc");
                TiXmlElement* adjusted_prec_mass =  prsm->FirstChildElement("adjusted_prec_mass");

                TiXmlElement* proteo_match_seq = proteoform->FirstChildElement("proteo_match_seq");

                TiXmlElement* extreme_value = prsm->FirstChildElement("extreme_value");
                TiXmlElement* e_value = extreme_value->FirstChildElement("e_value");

                outputFile <<"ms2.msalign"<<","<<"-1"<<","<<spectrum_id->GetText()<<","<<"CID"<<","<<spectrum_scan->GetText()<<","
                           <<"-"<<","<<"-"<<","<<"-"<<","<<ori_prec_mass->GetText()<<","<<adjusted_prec_mass->GetText()<<","
                           <<"0"<<","<<"-"<<","<<"-1"<<","<<"-1"<<","<<"1"<<","
                           <<seq_name->GetText()<<","<<"-"<<","<<start_pos->GetText()<<","<<end_pos->GetText()<<","<<"-"<<","
                           <<"-"<<","<<proteo_match_seq->GetText()<<","<<adjusted_prec_mass->GetText()<<","<<e_value->GetText()
                           <<","<<match_peak_num->GetText()<<","<<match_fragment_num->GetText()<<std::endl;




            }
        } else {
            std::cerr << "未找到 prsm_list 节点" << std::endl;
        }
    } else {
        std::cerr << "无法加载 XML 文件" << std::endl;
    }
    outputFile.close();

}


void PrsmXmlToCsv(string evalueFile){

    string filename = evalueFile;
    TiXmlDocument doc(filename.c_str());
    bool loadOkay = doc.LoadFile();
    string out_filename = "out_all_prsm.csv";
    std::ofstream outputFile;
    outputFile.open(out_filename,ios::out|ios::trunc);
    if (loadOkay) {
        std::cout << "XML 文件加载成功" << std::endl;

        // 获取根元素
        TiXmlElement* prsm_list = doc.FirstChildElement();
        if (prsm_list) {
            outputFile <<"Spectrum ID"<<","<<"Scan(s)"<<","<<"Precursor mass"<<","<<"Adjusted precursor mass"<<","
                       <<"Protein accession"<<","<<"First residue"<<","<<"Last residue"<<","
                       <<"Proteoform"<<","<<"match peak"<<","<<"match fragment"<<","
                       <<"feature1 "<<","
                       <<"feature2 "<<","
                       <<"feature3 "<<","
                       <<"feature4 "<<","
                       <<"feature5 "<<","
                       <<"feature6 "<<","
                       <<"feature7 "<<","
                       <<"feature8 "<<","
                       <<"feature9 "<<","
                       <<"feature10"<<","
                       <<"feature11"<<","
                       <<"feature12"<<","
                       <<"feature13"<<","
                       <<"feature14"<<","
                       <<"feature15"<<","
                       <<"feature16"<<","
                       <<"feature17"<<","
                       <<"feature18"<<","
                       <<"feature19"<<","
                       <<"feature20"<<","
                       <<"feature21"<<","
                       <<"feature22"<<","
                       <<"feature23"<<","
                       <<"feature24"<<","
                       <<"feature25"<<","
                       <<"feature26"<<","
                       <<"feature27"<<","
                       <<"feature28"<<","
                       <<"feature29"<<","
                       <<"feature30"<<","
                       <<"feature31"<<","
                       <<"feature32"<<","
                       <<"feature33"<<","
                       <<"feature34"<<","
                       <<"feature35"<<","
                       <<"feature36"<<","
                       <<"feature37"<<","
                       <<"feature38"<<","
                       <<"feature39"<<","
                       <<"feature40"<<","
                       <<"feature41"<<","
                       <<"feature42"<<","
                       <<"feature43"<<","
                       <<"feature44"<<","
                       <<"feature45"<<","
                       <<"feature46"<<","
                       <<"feature47"<<","
                       <<"feature48"<<","
                       <<"feature49"<<","
                       <<"feature50"<<std::endl;

            for (TiXmlElement* prsm = prsm_list->FirstChildElement("prsm"); prsm != NULL; prsm = prsm->NextSiblingElement("prsm")) {
                // 获取 evalue 节点的文本值
                TiXmlElement* spectrum_id = prsm->FirstChildElement("spectrum_id");
                TiXmlElement* spectrum_scan = prsm->FirstChildElement("spectrum_scan");
                TiXmlElement* ori_prec_mass = prsm->FirstChildElement("ori_prec_mass");
                TiXmlElement* match_peak_num = prsm->FirstChildElement("match_peak_num");
                TiXmlElement* match_fragment_num = prsm->FirstChildElement("match_fragment_num");

                TiXmlElement* proteoform =  prsm->FirstChildElement("proteoform");
                TiXmlElement* fasta_seq =  proteoform->FirstChildElement("fasta_seq");
                TiXmlElement* start_pos =  proteoform->FirstChildElement("start_pos");
                TiXmlElement* end_pos =  proteoform->FirstChildElement("end_pos");
                TiXmlElement* seq_name =  fasta_seq->FirstChildElement("seq_name");
                TiXmlElement* adjusted_prec_mass =  prsm->FirstChildElement("adjusted_prec_mass");
                TiXmlElement* proteo_match_seq = proteoform->FirstChildElement("proteo_match_seq");


                TiXmlElement* feature1 = prsm->FirstChildElement("feature1");
                TiXmlElement* feature2 = prsm->FirstChildElement("feature2");
                TiXmlElement* feature3 = prsm->FirstChildElement("feature3");
                TiXmlElement* feature4 = prsm->FirstChildElement("feature4");
                TiXmlElement* feature5 = prsm->FirstChildElement("feature5");
                TiXmlElement* feature6 = prsm->FirstChildElement("feature6");
                TiXmlElement* feature7 = prsm->FirstChildElement("feature7");
                TiXmlElement* feature8 = prsm->FirstChildElement("feature8");
                TiXmlElement* feature9 = prsm->FirstChildElement("feature9");
                TiXmlElement* feature10 = prsm->FirstChildElement("feature10");

                TiXmlElement* feature11 = prsm->FirstChildElement("feature11");
                TiXmlElement* feature12 = prsm->FirstChildElement("feature12");
                TiXmlElement* feature13 = prsm->FirstChildElement("feature13");
                TiXmlElement* feature14 = prsm->FirstChildElement("feature14");
                TiXmlElement* feature15 = prsm->FirstChildElement("feature15");
                TiXmlElement* feature16 = prsm->FirstChildElement("feature16");
                TiXmlElement* feature17 = prsm->FirstChildElement("feature17");
                TiXmlElement* feature18 = prsm->FirstChildElement("feature18");
                TiXmlElement* feature19 = prsm->FirstChildElement("feature19");
                TiXmlElement* feature20 = prsm->FirstChildElement("feature20");

                TiXmlElement* feature21 = prsm->FirstChildElement("feature21");
                TiXmlElement* feature22 = prsm->FirstChildElement("feature22");
                TiXmlElement* feature23 = prsm->FirstChildElement("feature23");
                TiXmlElement* feature24 = prsm->FirstChildElement("feature24");
                TiXmlElement* feature25 = prsm->FirstChildElement("feature25");
                TiXmlElement* feature26 = prsm->FirstChildElement("feature26");
                TiXmlElement* feature27 = prsm->FirstChildElement("feature27");
                TiXmlElement* feature28 = prsm->FirstChildElement("feature28");
                TiXmlElement* feature29 = prsm->FirstChildElement("feature29");
                TiXmlElement* feature30 = prsm->FirstChildElement("feature30");

                TiXmlElement* feature31 = prsm->FirstChildElement("feature31");
                TiXmlElement* feature32 = prsm->FirstChildElement("feature32");
                TiXmlElement* feature33 = prsm->FirstChildElement("feature33");
                TiXmlElement* feature34 = prsm->FirstChildElement("feature34");
                TiXmlElement* feature35 = prsm->FirstChildElement("feature35");
                TiXmlElement* feature36 = prsm->FirstChildElement("feature36");
                TiXmlElement* feature37 = prsm->FirstChildElement("feature37");
                TiXmlElement* feature38 = prsm->FirstChildElement("feature38");
                TiXmlElement* feature39 = prsm->FirstChildElement("feature39");
                TiXmlElement* feature40 = prsm->FirstChildElement("feature40");

                TiXmlElement* feature41 = prsm->FirstChildElement("feature41");
                TiXmlElement* feature42 = prsm->FirstChildElement("feature42");
                TiXmlElement* feature43 = prsm->FirstChildElement("feature43");
                TiXmlElement* feature44 = prsm->FirstChildElement("feature44");
                TiXmlElement* feature45 = prsm->FirstChildElement("feature45");
                TiXmlElement* feature46 = prsm->FirstChildElement("feature46");
                TiXmlElement* feature47 = prsm->FirstChildElement("feature47");
                TiXmlElement* feature48 = prsm->FirstChildElement("feature48");
                TiXmlElement* feature49 = prsm->FirstChildElement("feature49");
                TiXmlElement* feature50 = prsm->FirstChildElement("feature50");



                outputFile <<spectrum_id->GetText()<<","<<spectrum_scan->GetText()<<","<<ori_prec_mass->GetText()<<","
                           <<adjusted_prec_mass->GetText()<<","<<seq_name->GetText()<<","<<start_pos->GetText()<<","
                           <<end_pos->GetText()<<","<<proteo_match_seq->GetText()<<","<<match_peak_num->GetText()<<","
                           <<match_fragment_num->GetText()<<","
                           <<feature1 ->GetText()<<","
                           <<feature2 ->GetText()<<","
                           <<feature3 ->GetText()<<","
                           <<feature4 ->GetText()<<","
                           <<feature5 ->GetText()<<","
                           <<feature6 ->GetText()<<","
                           <<feature7 ->GetText()<<","
                           <<feature8 ->GetText()<<","
                           <<feature9 ->GetText()<<","
                           <<feature10 ->GetText()<<","
                           <<feature11 ->GetText()<<","
                           <<feature12 ->GetText()<<","
                           <<feature13 ->GetText()<<","
                           <<feature14 ->GetText()<<","
                           <<feature15 ->GetText()<<","
                           <<feature16 ->GetText()<<","
                           <<feature17 ->GetText()<<","
                           <<feature18 ->GetText()<<","
                           <<feature19 ->GetText()<<","
                           <<feature20 ->GetText()<<","
                           <<feature21 ->GetText()<<","
                           <<feature22 ->GetText()<<","
                           <<feature23 ->GetText()<<","
                           <<feature24 ->GetText()<<","
                           <<feature25 ->GetText()<<","
                           <<feature26 ->GetText()<<","
                           <<feature27 ->GetText()<<","
                           <<feature28 ->GetText()<<","
                           <<feature29 ->GetText()<<","
                           <<feature30 ->GetText()<<","
                           <<feature31 ->GetText()<<","
                           <<feature32 ->GetText()<<","
                           <<feature33 ->GetText()<<","
                           <<feature34 ->GetText()<<","
                           <<feature35 ->GetText()<<","
                           <<feature36 ->GetText()<<","
                           <<feature37 ->GetText()<<","
                           <<feature38 ->GetText()<<","
                           <<feature39 ->GetText()<<","
                           <<feature40 ->GetText()<<","
                           <<feature41 ->GetText()<<","
                           <<feature42 ->GetText()<<","
                           <<feature43 ->GetText()<<","
                           <<feature44 ->GetText()<<","
                           <<feature45 ->GetText()<<","
                           <<feature46 ->GetText()<<","
                           <<feature47 ->GetText()<<","
                           <<feature48 ->GetText()<<","
                           <<feature49 ->GetText()<<","
                           <<feature50 ->GetText()<<std::endl;

            }
        } else {
            std::cerr << "未找到 prsm_list 节点" << std::endl;
        }
    } else {
        std::cerr << "无法加载 XML 文件" << std::endl;
    }
    outputFile.close();

}


void CombineAllXmlToCsv(vector<string> evalueFileVec){

    string out_filename = "out_all_prsm.csv";
    std::ofstream outputFile;
    outputFile.open(out_filename,ios::out|ios::trunc);
    if(outputFile){
        outputFile <<"Spectrum ID"<<","<<"Scan(s)"<<","<<"Precursor mass"<<","<<"Adjusted precursor mass"<<","
                   <<"Protein accession"<<","<<"First residue"<<","<<"Last residue"<<","
                   <<"Proteoform"<<","<<"match peak"<<","<<"match fragment"<<","
                   <<"feature1 "<<","
                   <<"feature2 "<<","
                   <<"feature3 "<<","
                   <<"feature4 "<<","
                   <<"feature5 "<<","
                   <<"feature6 "<<","
                   <<"feature7 "<<","
                   <<"feature8 "<<","
                   <<"feature9 "<<","
                   <<"feature10"<<","
                   <<"feature11"<<","
                   <<"feature12"<<","
                   <<"feature13"<<","
                   <<"feature14"<<","
                   <<"feature15"<<","
                   <<"feature16"<<","
                   <<"feature17"<<","
                   <<"feature18"<<","
                   <<"feature19"<<","
                   <<"feature20"<<","
                   <<"feature21"<<","
                   <<"feature22"<<","
                   <<"feature23"<<","
                   <<"feature24"<<","
                   <<"feature25"<<","
                   <<"feature26"<<","
                   <<"feature27"<<","
                   <<"feature28"<<","
                   <<"feature29"<<","
                   <<"feature30"<<","
                   <<"feature31"<<","
                   <<"feature32"<<","
                   <<"feature33"<<","
                   <<"feature34"<<","
                   <<"feature35"<<","
                   <<"feature36"<<","
                   <<"feature37"<<","
                   <<"feature38"<<","
                   <<"feature39"<<","
                   <<"feature40"<<","
                   <<"feature41"<<","
                   <<"feature42"<<","
                   <<"feature43"<<","
                   <<"feature44"<<","
                   <<"feature45"<<","
                   <<"feature46"<<","
                   <<"feature47"<<","
                   <<"feature48"<<","
                   <<"feature49"<<","
                   <<"feature50"<<std::endl;

    }
    for(string filename:evalueFileVec){
        TiXmlDocument doc(filename.c_str());
        bool loadOkay = doc.LoadFile();
        if (loadOkay) {
            std::cout<<filename<< " 文件加载成功" << std::endl;
            // 获取根元素
            TiXmlElement* prsm_list = doc.FirstChildElement();
            if (prsm_list) {
                for (TiXmlElement* prsm = prsm_list->FirstChildElement("prsm"); prsm != NULL; prsm = prsm->NextSiblingElement("prsm")) {
                    // 获取 evalue 节点的文本值
                    TiXmlElement* spectrum_id = prsm->FirstChildElement("spectrum_id");
                    TiXmlElement* spectrum_scan = prsm->FirstChildElement("spectrum_scan");
                    TiXmlElement* ori_prec_mass = prsm->FirstChildElement("ori_prec_mass");
                    TiXmlElement* match_peak_num = prsm->FirstChildElement("match_peak_num");
                    TiXmlElement* match_fragment_num = prsm->FirstChildElement("match_fragment_num");

                    TiXmlElement* proteoform =  prsm->FirstChildElement("proteoform");
                    TiXmlElement* fasta_seq =  proteoform->FirstChildElement("fasta_seq");
                    TiXmlElement* start_pos =  proteoform->FirstChildElement("start_pos");
                    TiXmlElement* end_pos =  proteoform->FirstChildElement("end_pos");
                    TiXmlElement* seq_name =  fasta_seq->FirstChildElement("seq_name");
                    TiXmlElement* adjusted_prec_mass =  prsm->FirstChildElement("adjusted_prec_mass");
                    TiXmlElement* proteo_match_seq = proteoform->FirstChildElement("proteo_match_seq");


                    TiXmlElement* feature1 = prsm->FirstChildElement("feature1");
                    TiXmlElement* feature2 = prsm->FirstChildElement("feature2");
                    TiXmlElement* feature3 = prsm->FirstChildElement("feature3");
                    TiXmlElement* feature4 = prsm->FirstChildElement("feature4");
                    TiXmlElement* feature5 = prsm->FirstChildElement("feature5");
                    TiXmlElement* feature6 = prsm->FirstChildElement("feature6");
                    TiXmlElement* feature7 = prsm->FirstChildElement("feature7");
                    TiXmlElement* feature8 = prsm->FirstChildElement("feature8");
                    TiXmlElement* feature9 = prsm->FirstChildElement("feature9");
                    TiXmlElement* feature10 = prsm->FirstChildElement("feature10");

                    TiXmlElement* feature11 = prsm->FirstChildElement("feature11");
                    TiXmlElement* feature12 = prsm->FirstChildElement("feature12");
                    TiXmlElement* feature13 = prsm->FirstChildElement("feature13");
                    TiXmlElement* feature14 = prsm->FirstChildElement("feature14");
                    TiXmlElement* feature15 = prsm->FirstChildElement("feature15");
                    TiXmlElement* feature16 = prsm->FirstChildElement("feature16");
                    TiXmlElement* feature17 = prsm->FirstChildElement("feature17");
                    TiXmlElement* feature18 = prsm->FirstChildElement("feature18");
                    TiXmlElement* feature19 = prsm->FirstChildElement("feature19");
                    TiXmlElement* feature20 = prsm->FirstChildElement("feature20");

                    TiXmlElement* feature21 = prsm->FirstChildElement("feature21");
                    TiXmlElement* feature22 = prsm->FirstChildElement("feature22");
                    TiXmlElement* feature23 = prsm->FirstChildElement("feature23");
                    TiXmlElement* feature24 = prsm->FirstChildElement("feature24");
                    TiXmlElement* feature25 = prsm->FirstChildElement("feature25");
                    TiXmlElement* feature26 = prsm->FirstChildElement("feature26");
                    TiXmlElement* feature27 = prsm->FirstChildElement("feature27");
                    TiXmlElement* feature28 = prsm->FirstChildElement("feature28");
                    TiXmlElement* feature29 = prsm->FirstChildElement("feature29");
                    TiXmlElement* feature30 = prsm->FirstChildElement("feature30");

                    TiXmlElement* feature31 = prsm->FirstChildElement("feature31");
                    TiXmlElement* feature32 = prsm->FirstChildElement("feature32");
                    TiXmlElement* feature33 = prsm->FirstChildElement("feature33");
                    TiXmlElement* feature34 = prsm->FirstChildElement("feature34");
                    TiXmlElement* feature35 = prsm->FirstChildElement("feature35");
                    TiXmlElement* feature36 = prsm->FirstChildElement("feature36");
                    TiXmlElement* feature37 = prsm->FirstChildElement("feature37");
                    TiXmlElement* feature38 = prsm->FirstChildElement("feature38");
                    TiXmlElement* feature39 = prsm->FirstChildElement("feature39");
                    TiXmlElement* feature40 = prsm->FirstChildElement("feature40");

                    TiXmlElement* feature41 = prsm->FirstChildElement("feature41");
                    TiXmlElement* feature42 = prsm->FirstChildElement("feature42");
                    TiXmlElement* feature43 = prsm->FirstChildElement("feature43");
                    TiXmlElement* feature44 = prsm->FirstChildElement("feature44");
                    TiXmlElement* feature45 = prsm->FirstChildElement("feature45");
                    TiXmlElement* feature46 = prsm->FirstChildElement("feature46");
                    TiXmlElement* feature47 = prsm->FirstChildElement("feature47");
                    TiXmlElement* feature48 = prsm->FirstChildElement("feature48");
                    TiXmlElement* feature49 = prsm->FirstChildElement("feature49");
                    TiXmlElement* feature50 = prsm->FirstChildElement("feature50");



                    outputFile <<spectrum_id->GetText()<<","<<spectrum_scan->GetText()<<","<<ori_prec_mass->GetText()<<","
                               <<adjusted_prec_mass->GetText()<<","<<seq_name->GetText()<<","<<start_pos->GetText()<<","
                               <<end_pos->GetText()<<","<<proteo_match_seq->GetText()<<","<<match_peak_num->GetText()<<","
                               <<match_fragment_num->GetText()<<","
                               <<feature1 ->GetText()<<","
                               <<feature2 ->GetText()<<","
                               <<feature3 ->GetText()<<","
                               <<feature4 ->GetText()<<","
                               <<feature5 ->GetText()<<","
                               <<feature6 ->GetText()<<","
                               <<feature7 ->GetText()<<","
                               <<feature8 ->GetText()<<","
                               <<feature9 ->GetText()<<","
                               <<feature10 ->GetText()<<","
                               <<feature11 ->GetText()<<","
                               <<feature12 ->GetText()<<","
                               <<feature13 ->GetText()<<","
                               <<feature14 ->GetText()<<","
                               <<feature15 ->GetText()<<","
                               <<feature16 ->GetText()<<","
                               <<feature17 ->GetText()<<","
                               <<feature18 ->GetText()<<","
                               <<feature19 ->GetText()<<","
                               <<feature20 ->GetText()<<","
                               <<feature21 ->GetText()<<","
                               <<feature22 ->GetText()<<","
                               <<feature23 ->GetText()<<","
                               <<feature24 ->GetText()<<","
                               <<feature25 ->GetText()<<","
                               <<feature26 ->GetText()<<","
                               <<feature27 ->GetText()<<","
                               <<feature28 ->GetText()<<","
                               <<feature29 ->GetText()<<","
                               <<feature30 ->GetText()<<","
                               <<feature31 ->GetText()<<","
                               <<feature32 ->GetText()<<","
                               <<feature33 ->GetText()<<","
                               <<feature34 ->GetText()<<","
                               <<feature35 ->GetText()<<","
                               <<feature36 ->GetText()<<","
                               <<feature37 ->GetText()<<","
                               <<feature38 ->GetText()<<","
                               <<feature39 ->GetText()<<","
                               <<feature40 ->GetText()<<","
                               <<feature41 ->GetText()<<","
                               <<feature42 ->GetText()<<","
                               <<feature43 ->GetText()<<","
                               <<feature44 ->GetText()<<","
                               <<feature45 ->GetText()<<","
                               <<feature46 ->GetText()<<","
                               <<feature47 ->GetText()<<","
                               <<feature48 ->GetText()<<","
                               <<feature49 ->GetText()<<","
                               <<feature50 ->GetText()<<std::endl;

                }
                std::remove(filename.c_str());
            } else {
                std::cerr<<filename << " 未找到 prsm_list 节点" << std::endl;
            }
        } else {
            std::cerr<<filename << " 无法加载 XML 文件" << std::endl;
        }
        doc.Clear();
    }


    outputFile.close();

}


void combined_xml(vector<string> combinedPrsmFiles,string combinedPrsmFile){
    string filename = combinedPrsmFile;

    TiXmlDocument doc(filename.c_str());


    FILE* file = fopen(filename.c_str(), "rb");
    bool xmlExist;
    if (file) {
        xmlExist = true;
        fclose(file);
    }else{
        xmlExist = false;
    }
    if(xmlExist){
        //xml文件已存在，添加新prsm结点
        TiXmlDocument doc(filename.c_str());
        doc.LoadFile();
        TiXmlElement* prsm_list = doc.FirstChildElement("prsm_list");
        if (prsm_list) {
            for(string tempFileName:combinedPrsmFiles){
                TiXmlDocument tempDoc(tempFileName.c_str());
                if(tempDoc.LoadFile()){
                    TiXmlElement* temp_prsm_list = tempDoc.FirstChildElement();
                    for (TiXmlElement* tempPrsm = temp_prsm_list->FirstChildElement("prsm"); tempPrsm != NULL; tempPrsm = tempPrsm->NextSiblingElement("prsm"))
                    {
                        TiXmlElement* prsmCopy = new TiXmlElement(*tempPrsm);
                        prsm_list->LinkEndChild(prsmCopy);
                    }
                }else{
                    cout<<endl<<"待合并文件错误"<<endl;
                    exit(718);
                }

                tempDoc.Clear();
                cout<<tempFileName<<" finish"<<endl;
                std::remove(tempFileName.c_str());
            }
            doc.SaveFile(filename.c_str());
            doc.Clear();
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
        for(string tempFileName:combinedPrsmFiles){
            TiXmlDocument tempDoc(tempFileName.c_str());
            if(tempDoc.LoadFile()){
                TiXmlElement* temp_prsm_list = tempDoc.FirstChildElement();
                for (TiXmlElement* tempPrsm = temp_prsm_list->FirstChildElement("prsm"); tempPrsm != NULL; tempPrsm = tempPrsm->NextSiblingElement("prsm"))
                {
                    TiXmlElement* prsmCopy = new TiXmlElement(*tempPrsm);
                    prsm_list->LinkEndChild(prsmCopy);
                }
            }else{
                cout<<endl<<"待合并文件错误"<<endl;
                exit(718);
            }

            tempDoc.Clear();
            std::remove(tempFileName.c_str());
            cout<<tempFileName<<" finish"<<endl;
        }

        doc.SaveFile(filename.c_str());
        doc.Clear();
    }

    fclose(file);

}


void conditional_combined_xml(vector<string> combinedPrsmFiles,string combinedPrsmFile){
    string filename = combinedPrsmFile;

    TiXmlDocument doc(filename.c_str());


    FILE* file = fopen(filename.c_str(), "rb");
    bool xmlExist;
    if (file) {
        xmlExist = true;
        fclose(file);
    }else{
        xmlExist = false;
    }
    if(xmlExist){
        //xml文件已存在，添加新prsm结点
        TiXmlDocument doc(filename.c_str());
        doc.LoadFile();
        TiXmlElement* prsm_list = doc.FirstChildElement("prsm_list");
        if (prsm_list) {
            for(string tempFileName:combinedPrsmFiles){
                if(access(tempFileName.c_str(),F_OK)!=0)
                    continue;
                TiXmlDocument tempDoc(tempFileName.c_str());
                if(tempDoc.LoadFile()){
                    TiXmlElement* temp_prsm_list = tempDoc.FirstChildElement();
                    for (TiXmlElement* tempPrsm = temp_prsm_list->FirstChildElement("prsm"); tempPrsm != NULL; tempPrsm = tempPrsm->NextSiblingElement("prsm"))
                    {
                        TiXmlElement* prsmCopy = new TiXmlElement(*tempPrsm);
                        TiXmlElement* proteoform = prsmCopy->FirstChildElement("proteoform");
                        TiXmlElement* match_peak_num = prsmCopy->FirstChildElement("match_peak_num");

                        TiXmlElement* variable_ptm_num = proteoform->FirstChildElement("variable_ptm_num");
                        TiXmlElement* unexpected_ptm_num = proteoform->FirstChildElement("unexpected_ptm_num");
                        TiXmlElement* start_pos = proteoform->FirstChildElement("start_pos");
                        TiXmlElement* end_pos = proteoform->FirstChildElement("end_pos");

                        int ptm_num = stoi(variable_ptm_num->GetText());
                        int unk_num = stoi(unexpected_ptm_num->GetText());
                        int start_num = stoi(start_pos->GetText());
                        int end_num = stoi(end_pos->GetText());
                        int peak_num = stoi(match_peak_num->GetText());

                        bool isAdd= true;

                        if(peak_num < 5){
                            isAdd = false;
                        }

                        if(ptm_num + unk_num>= 5) {
                            isAdd = false;
                        }
                        if(start_num>=1000 ||end_num>=1000)
                            isAdd = false;
                        if(isAdd)
                            prsm_list->LinkEndChild(prsmCopy);
                    }
                }else{
                    cout<<endl<<"待合并文件错误"<<endl;
                    exit(718);
                }

                tempDoc.Clear();
                cout<<tempFileName<<" finish"<<endl;
                std::remove(tempFileName.c_str());
            }
            doc.SaveFile(filename.c_str());
            doc.Clear();
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
        for(string tempFileName:combinedPrsmFiles){
            if(access(tempFileName.c_str(),F_OK)!=0)
                continue;
            TiXmlDocument tempDoc(tempFileName.c_str());
            if(tempDoc.LoadFile()){
                TiXmlElement* temp_prsm_list = tempDoc.FirstChildElement();
                for (TiXmlElement* tempPrsm = temp_prsm_list->FirstChildElement("prsm"); tempPrsm != NULL; tempPrsm = tempPrsm->NextSiblingElement("prsm"))
                {
                    TiXmlElement* prsmCopy = new TiXmlElement(*tempPrsm);
                    TiXmlElement* proteoform = prsmCopy->FirstChildElement("proteoform");
                    TiXmlElement* match_peak_num = prsmCopy->FirstChildElement("match_peak_num");

                    TiXmlElement* variable_ptm_num = proteoform->FirstChildElement("variable_ptm_num");
                    TiXmlElement* unexpected_ptm_num = proteoform->FirstChildElement("unexpected_ptm_num");
                    TiXmlElement* start_pos = proteoform->FirstChildElement("start_pos");
                    TiXmlElement* end_pos = proteoform->FirstChildElement("end_pos");

                    int ptm_num = stoi(variable_ptm_num->GetText());
                    int unk_num = stoi(unexpected_ptm_num->GetText());
                    int start_num = stoi(start_pos->GetText());
                    int end_num = stoi(end_pos->GetText());
                    int peak_num = stoi(match_peak_num->GetText());

                    bool isAdd= true;

                    if(peak_num < 5){
                        isAdd = false;
                    }

                    if(ptm_num + unk_num>= 5) {
                        isAdd = false;
                    }
                    if(start_num>=1000 ||end_num>=1000)
                        isAdd = false;
                    if(isAdd)
                        prsm_list->LinkEndChild(prsmCopy);
                }
            }else{
                cout<<endl<<"待合并文件错误"<<endl;
                exit(718);
            }

            tempDoc.Clear();
            std::remove(tempFileName.c_str());
            cout<<tempFileName<<" finish"<<endl;
        }

        doc.SaveFile(filename.c_str());
        doc.Clear();
    }

    fclose(file);

}


void select_combined_xml(string combinedPrsmFile,string thread_id){
    string filename = combinedPrsmFile;
    string tempFileName = filename + "_" +thread_id;
    TiXmlDocument doc(filename.c_str());
    TiXmlDocument docTemp(tempFileName.c_str());

    FILE* file = fopen(filename.c_str(), "rb");
    bool xmlExist;
    if (file) {
        xmlExist = true;
        fclose(file);
    }else{
        xmlExist = false;
    }
    if(xmlExist){
        //xml文件已存在，添加新prsm结点
        TiXmlDocument doc(filename.c_str());
        doc.LoadFile();
        TiXmlElement* prsm_list = doc.FirstChildElement("prsm_list");
        if (prsm_list) {
            TiXmlDocument tempDoc(tempFileName.c_str());
            if(tempDoc.LoadFile()){
                TiXmlElement* temp_prsm_list = tempDoc.FirstChildElement();
                for (TiXmlElement* tempPrsm = temp_prsm_list->FirstChildElement("prsm"); tempPrsm != NULL; tempPrsm = tempPrsm->NextSiblingElement("prsm"))
                {
                    TiXmlElement* proteoform = tempPrsm->FirstChildElement("proteoform");
                    TiXmlElement* variable_ptm_num = proteoform->FirstChildElement("variable_ptm_num");
                    TiXmlElement* unexpected_ptm_num = proteoform->FirstChildElement("unexpected_ptm_num");
                    int ptm_num = stoi(variable_ptm_num->GetText());
                    int unk_num = stoi(unexpected_ptm_num->GetText());
                    if(ptm_num <= 5 && unk_num<= 5) {
                        TiXmlElement *prsmCopy = new TiXmlElement(*tempPrsm);
                        prsm_list->LinkEndChild(prsmCopy);
                    }
                    std::remove(tempFileName.c_str());
                }
            }else{
                cout<<endl<<"待合并文件错误"<<endl;
                exit(718);
            }
            doc.SaveFile(filename.c_str());
            doc.Clear();
            tempDoc.Clear();

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
        TiXmlDocument tempDoc(tempFileName.c_str());
        if(tempDoc.LoadFile()){
            TiXmlElement* temp_prsm_list = tempDoc.FirstChildElement();
            for (TiXmlElement* tempPrsm = temp_prsm_list->FirstChildElement("prsm"); tempPrsm != NULL; tempPrsm = tempPrsm->NextSiblingElement("prsm"))
            {
                TiXmlElement* proteoform = tempPrsm->FirstChildElement("proteoform");
                TiXmlElement* variable_ptm_num = proteoform->FirstChildElement("variable_ptm_num");
                TiXmlElement* unexpected_ptm_num = proteoform->FirstChildElement("unexpected_ptm_num");
                int ptm_num = stoi(variable_ptm_num->GetText());
                int unk_num = stoi(unexpected_ptm_num->GetText());
                if(ptm_num<=5 && unk_num<= 5) {
                    TiXmlElement *prsmCopy = new TiXmlElement(*tempPrsm);
                    prsm_list->LinkEndChild(prsmCopy);
                }
                std::remove(tempFileName.c_str());
            }
        }else{
            cout<<endl<<"待合并文件错误"<<endl;
            exit(718);
        }
        doc.SaveFile(filename.c_str());
        doc.Clear();
        tempDoc.Clear();
    }

    fclose(file);

}


void getToppicFilterResult(string toppicCombinedFile,string output_filename){
    fstream f;
    f.open(output_filename,ios::out);
    //输入你想写入的内容

    TiXmlDocument tempDoc(toppicCombinedFile.c_str());
    if(tempDoc.LoadFile()){
        TiXmlElement* temp_prsm_list = tempDoc.FirstChildElement();
        for (TiXmlElement* tempPrsm = temp_prsm_list->FirstChildElement("prsm"); tempPrsm != NULL; tempPrsm = tempPrsm->NextSiblingElement("prsm"))
        {
            TiXmlElement* spectrum_id = tempPrsm->FirstChildElement("spectrum_id");
            TiXmlElement* proteoform = tempPrsm->FirstChildElement("proteoform");
            TiXmlElement* fasta_seq = proteoform->FirstChildElement("fasta_seq");
            TiXmlElement* seq_name = fasta_seq->FirstChildElement("seq_name");

            f<<spectrum_id->GetText()<<" "<<seq_name->GetText()<<endl;

        }
    }else{
        cout<<"加载Toppic_combined文件失败"<<endl;
    }
    f.close();
}

void outXmlForWTopMLS(){
    string keyword = "wtop_prsm.test_combined";
    vector<string> files;
    string directory = "./";  // 目标文件夹（可修改为具体路径）
    // 遍历当前目录，查找匹配的文件名
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.is_regular_file()) {
            string filename = entry.path().filename().string();
            if (filename.find(keyword) != string::npos) {
                files.push_back(entry.path().string()); // 记录完整路径
            }
        }
    }
    // 调用函数处理文件列表
    CombineAllXmlToCsv(files);
    cout << endl << "Conditional combine XML finish" << endl;
}



struct ProteinEntry {
    std::string name;
    size_t length;
};

// 解析FASTA文件
std::vector<ProteinEntry> parseFasta(const std::string& filename) {
    std::vector<ProteinEntry> proteins;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return proteins;
    }

    std::string line, name, sequence;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') { // 标题行
            if (!name.empty()) {
                proteins.push_back({name, sequence.length()});
                sequence.clear();
            }
            size_t spacePos = line.find(' ');
            name = (spacePos != std::string::npos) ? line.substr(1, spacePos - 1) : line.substr(1); // 去掉 '>' 并截取第一个空格前的内容
        } else {
            sequence += line; // 连接序列
        }
    }
    if (!name.empty()) {
        proteins.push_back({name, sequence.length()});
    }

    return proteins;
}


// 输出到CSV文件
void exportToCSV(const std::vector<ProteinEntry>& proteins, const std::string& outputFilename) {
    std::ofstream file(outputFilename);
    if (!file.is_open()) {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }
    file << "Protein Name,Sequence Length\n";
    for (const auto& protein : proteins) {
        file << protein.name << "," << protein.length << "\n";
    }
}

//读取打分
//int main() {
//
////    getToppicFilterResult("Experiment_4_620_F5_01_ms2.toppic_combined","f5_01_total_result.txt");
//
////
//
////    string filename = "D:\\LYC\\wtopMLS\\528\\F8_1\\copy_file\\wtop_prsm.test_combined";
////    vector<string> files;
////    for(int i=0;i<=1910;i++){
////        string tempfile = filename + "_" + to_string(i);
////        files.push_back(tempfile);
////    }
////    conditional_combined_xml(files,filename);
////    cout<<endl<<"conditional combin xml finish"<<endl;
//
//
//
////    sortPrsm("wtop_prsm.test_combined");
//
////    evalueToCsv("wtop_prsm.test_combined");
//
//    string filename = "wtop_prsm.test_combined";
//    vector<string> files;
//    for(int i=990;i<=991;i++){
//        string tempfile = filename + "_" + to_string(i);
//        files.push_back(tempfile);
//    }
//    CombineAllXmlToCsv(files);
//    cout<<endl<<"conditional combin xml finish"<<endl;
//
//
//
//    return 0;
//}




///统计序列
//int main() {
//    std::string inputFilename = "uniprotkb_UP000000589_AND_model_organis_2024_04_24_target_decoy.fasta";
//    std::string outputFilename = "proteins.csv";
//
//    std::vector<ProteinEntry> proteins = parseFasta(inputFilename);
//    exportToCSV(proteins, outputFilename);
//
//    std::cout << "Exported to " << outputFilename << " successfully!" << std::endl;
//    return 0;
//}