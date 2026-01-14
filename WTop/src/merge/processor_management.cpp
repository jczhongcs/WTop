//
// Created by wenzhong on 2023/3/15.
//

#include "processor_management.h"

namespace prsm {
    namespace processor_management {
        void save_best_prsm_information(protein_processor_ptr protein,
                                        msalign_processor_ptr msalign,
                                        temp_prsm_argument_ptr temp_arg) {
            msalign->best_score = temp_arg->score; //分值
            msalign->best_peak = temp_arg->peaks_match;    //匹配峰数
            msalign->matchProtein = protein; //匹配的蛋白
            msalign->match_protein_name = protein->protein_name;
            msalign->match_protein_seq = protein->protein_sequence;
            msalign->cut_c_n = temp_arg->cut_location_n; // n cut size
            msalign->cut_n_c = temp_arg->cut_location_c;  //c cut size
            msalign->best_ptm = temp_arg->ptms;        // var ptm ptr
            msalign->modify_mass = temp_arg->var_mod_mass;;    // ptms mass
            msalign->mut_s = temp_arg->mut_x;  // mut s
            msalign->bestIonsPairNub = temp_arg->pairIonsNub;
            msalign->huBuCount = temp_arg->complement_ions_number;
            msalign->seqScore = temp_arg->prsm_score;
            msalign->seqScoreCount = temp_arg->seq_score_count;
            msalign->matchFragmentIonsSize = temp_arg->matchFragmentIonsSize;
            msalign->ionsMessageContainerN = temp_arg->ions_n;
            msalign->ionsMessageContainerC = temp_arg->ions_c;
            msalign->theoryMassN = temp_arg->theory_ions_mass_n;
            msalign->theoryMassC = temp_arg->theory_ions_mass_c;
        }


        /**
 * 初始化修饰
 * @param ptmsPtr
 * @param one_mod_ptr
 * @param modifier_table
 * @param one_ptm_table
 */
        void initializerPtmsPtr(prsm::modify_ptr ptmsPtr, prsm::modify_ptr one_mod_ptr,vector<double> &one_ptm_table) {
            ptmsPtr->Init_modify();                  //初始化修饰表
            std::map<double, string> modifyTable = ptmsPtr->getTable();       //修饰表保存    double ,string
            map<double, string> one_ptm_map;
            for (std::map<double, string>::iterator mit = modifyTable.begin(); mit != modifyTable.end(); mit++) {
                if (mit->second.size() == 1) {       //保存一个修饰的质量
                    one_ptm_map[mit->first] = mit->second;
                    one_ptm_table.push_back(mit->first);
                }
            }
            one_mod_ptr->get_mod_Map(one_ptm_map);
            if (ptmsPtr->modifyMassSeq.size() == 0) {
                EX_TRACE("The modification_table_file can't find !\n");
                exit(1);
            }
        }

        /**
        * 读取质谱数据
        * @param msFileName
        * @param msContainer
        */
        void init_msalign_container(const string & ms_file_name,
                                    vector<msalign_ptr> &msContainer)
        {
            std::ifstream inMsFile(ms_file_name, ios::in);
            if (!inMsFile.good()) {
                EX_TRACE("Can't Find The Ms File !\n");
                exit(1);
            }
            while (!inMsFile.eof()) {       //读取质谱数据
                string buff;
                getline(inMsFile, buff);
                if (buff.find("BEGIN") == 0) {  //读到标志位，开始往下读取质谱数据
                    msalign_ptr msalignPtr = std::make_shared<msalign>();  //创建质谱指针.初始化质谱数据
                    msalignPtr->init_ms_msg(inMsFile); //初始化数据
                    if(msalignPtr->getIonsMassContainer().size() > 5) {
                        msContainer.push_back(msalignPtr);
                    }
                }
            }
            //初始化质谱互补离子对下标
            for(int i = 0 ; i < msContainer.size(); ++i) {
                msContainer[i]->initCompIonsIndexMap();
            }


        }


        /**
        * 初始化标准蛋白质数据库
        * @param protein_file
        * @param prContainer
        */
        void init_protein_container(const string & protein_file,
                                    vector<prsm::protein_processor_ptr> &prContainer)
        {
            double cur_pro_num = 0;     //记录蛋白数量
            char *cProteinFileName = (char *) malloc((protein_file.length() + 1) * sizeof(char));
            strcpy(cProteinFileName, protein_file.c_str());
            map<string, string> proteinTileSeq;      //蛋白质文件的title和sequence
            mylib::io::get_protein_name_seq_map(cProteinFileName, proteinTileSeq);        //用map保存每条蛋白质(first->name,second->seq.)
            for (map<string, string>::iterator pro_it = proteinTileSeq.begin();
            pro_it != proteinTileSeq.end(); pro_it++) {
                prsm::protein_processor_ptr p = new prsm::Protein(pro_it->first,pro_it->second);
                prContainer.push_back(p);   //添加进蛋白库

                cout << "\rProtein Has Inited " << setiosflags(ios::fixed) << setprecision(4)
                     << ++cur_pro_num / proteinTileSeq.size() * 100 << " %";
                fflush(stdout);

            }
            free(cProteinFileName);
        }

//20221114,将PrSM结果写入文本
        void writeResultNews(msalign_processor_ptr msIt,
                             string Result_file,
                             modify_ptr modptr)
        {
            if (msIt->match_protein_name.size() > 0)
            {
                std::ofstream out(Result_file, ios::out | ios::app);
                double t_mass = msIt->precursor_mass_double; //atof(msIt->Pre_mass.substr(msIt->Pre_mass.find('=') + 1, msIt->Pre_mass.size()).c_str());
                out << "#BEGIN" << endl;
                out << "#" << msIt->Id
                    << " , #" << msIt->scans
                    << ", #Types : " << msIt->activation.substr(msIt->activation.find('=') + 1, msIt->activation.size())
                    << ", #Protero : " << msIt->match_protein_name << "." << endl;
                out << "Match Peaks = " << msIt->best_peak
                    << ",score = " << msIt->best_score << endl;
                out << "CutLocationN = " << msIt->cut_c_n
                    << " ,CutLocationC = " << msIt->cut_n_c << endl;
                out << "Precursor mass = " << setiosflags(ios::fixed) << setprecision(4)
                    << msIt->precursor_mass.substr(msIt->precursor_mass.find('=') + 1, msIt->precursor_mass.size())
                    << " ,Proteoform mass = " << msIt->theory_proteoform_mass << endl;
                out << "PPM = " << (t_mass - msIt->theory_proteoform_mass) / t_mass * 1000000 << endl;

                //处理截断
                utils_string_util::processCutLocation(msIt->match_protein_seq, msIt->cut_c_n, msIt->cut_n_c);
                //输出修饰
                out << "#Mod BEGIN" << endl;
                double pMass = 0;
                map<int,int> startEnd;  //修饰所在的序列区间
                vector<string> varPtmsContainer ;  //保存修饰
                for (int j = 0; j < msIt->best_ptm.size(); j++) {
                    string ptms = modptr->analysis(msIt->best_ptm[j].mod_mass);
                    if (ptms != "0" && msIt->modify_mass != 0) {
                        out << "start: " << (int) msIt->best_ptm[j].first << "-" << (int) msIt->best_ptm[j].second
                            << ";mass: " << msIt->best_ptm[j].mod_mass
                            << " ; Mod# " << ptms << endl;
                        pMass += msIt->best_ptm[j].mod_mass;
                        varPtmsContainer.push_back(ptms) ;
                        startEnd.insert(make_pair(msIt->best_ptm[j].first,msIt->best_ptm[j].second));
                    }
                }
                out << "#Mod END" << endl;
                for (auto it = msIt->unknownPTMMassMap.begin(); it != msIt->unknownPTMMassMap.end();++it)
                {
                    out << it->first <<" = "<<it->second<<endl;
                }
                string pfSeq = "";  //蛋白质形序列
//        for(auto it = modptr->charPTMsMap.begin(); it != modptr->charPTMsMap.end();++it)
//        {
//            cout <<it->first <<" "<< it->second <<endl;
//        }
                //处理蛋白质形
                utils_string_util::processProteinFom(msIt->match_protein_seq, varPtmsContainer, modptr->charPTMsMap, startEnd, pfSeq);
                out << "RawProteinSeq# " << msIt->match_protein_seq << endl;
                out <<"ProteoForm# "<< pfSeq <<endl;
                //其余数据输出
                out << "Spact ModifyMass = " << msIt->modify_mass << "\nTrue ModifyMass = " << pMass << endl;
                out << "All peaks = " << msIt->ions_mass_container.size() << endl;
                out << "match peaks = " << (int) msIt->best_peak << endl;
                out << "Ms complementIons =  " << msIt->compIonsIndexMap.size() << endl;
                out << "matched fragment ions = " << msIt->matchFragmentIonsSize << endl;
                out << "mut_x = " << msIt->mut_s << endl;
                out << "ionsPair = " << msIt->bestIonsPairNub << endl;
                out << "huBuIonsCount = " << (int) msIt->huBuCount - 1 << endl;
                out << "seqIonsCount = " << msIt->seqScoreCount << endl;
                out << "seqScore = " << msIt->seqScore << endl;
                out << "#END\n" << endl;
                out.close();
            }   //输出答案结束
        }


        //将结果写入文件，仅供测试数据使用
        void WriteSequenceAlignment(const char *output,const std::vector<double> &reference_orig, const std::vector<double> &peer_orig,
                                    vector<pair<long, long> > &alignment, vector<double> &sub,vector<double> &index, vector<double> &ppm,
                                    vector<double> &RefeZscore, vector<double> &PeerZscore)
        {
            vector<std::string> tmp_rec;
            std::ostringstream title;
            title << setw(10) << "Thero ID" << setw(10) << " Mono ID" << "  | " << setw(15) << " Thero Mass" << setw(15)
                  << "Mono Mass" << "  | " << setw(10) << "Mono-Thero" << " ｜ " << setw(10) << "Index" << " ｜ " << setw(10)
                  << "Ppm" << " ｜ ";
            tmp_rec.push_back(title.str());
            for (long i = 0; i < alignment.size(); i++) {
                //----- output to string ----//
                std::ostringstream o;
                o << setw(10) << alignment[i].first << " " << setw(10) << alignment[i].second << " | ";
                o << setw(15) << setiosflags(ios::fixed) << setprecision(4) << reference_orig[alignment[i].first] << " "
                  << setw(15) << setiosflags(ios::fixed) << setprecision(4) << peer_orig[alignment[i].second] << " | ";
                o << setw(10) << setiosflags(ios::fixed) << setprecision(4) << sub[i] << " ｜ ";
                if (index[i] != INT_MAX) {
                    if (fabs(ppm[i]) < 15) {
                        // o<<setw(10)<<alignment[i].firstwe<<" "<<setw(10)<<alignment[i].second<<" | ";
                        // o<<setw(15)<<setiosflags(ios::fixed)<<setprecision(4)<<reference_orig[alignment[i].first+start]<<" "<<setw(15)<<setiosflags(ios::fixed)<<setprecision(4)<<peer_orig[alignment[i].second]<<" | ";
                        o << setw(10) << setiosflags(ios::fixed) << setprecision(4) << index[i] << " ｜ ";
                        o << setw(10) << setiosflags(ios::fixed) << setprecision(4) << ppm[i] << " ｜ ";
                        o << setw(10) << setiosflags(ios::fixed) << setprecision(4) << RefeZscore[alignment[i].first] << " ｜ ";
                        o << setw(10) << setiosflags(ios::fixed) << setprecision(4) << PeerZscore[alignment[i].second] << " ｜ ";
                        // std::string s=o.str();
                        // tmp_rec.push_back(s);
                    }
                }
                //----- record string -----//
                std::string s = o.str();
                tmp_rec.push_back(s);
            }
            //----- output to file ------//
            FILE *fp = fopen(output, "wb");
            for (long i = 0; i < (long) tmp_rec.size(); i++) fprintf(fp, "%s\n", tmp_rec[i].c_str());
            fclose(fp);
        }

        /**
        * 适用于一条质谱对应一条蛋白的测试输出
        * @param pI 标识选取端
        * @param arg 参数
        * @param spact 质谱
        */
        void test_out(char pI, temp_prsm_argument_ptr arg, msalign_processor_ptr spact)
        {
            cout << "\n选取的是 : " << pI << " 端 !" << endl;
            cout << "N端选取的离子为 : \n";
            for (int pi = 0; pi < arg->ions_n.size(); pi++) {
                cout << " T : " << arg->ions_n[pi].thoe_id << " M : " << arg->ions_n[pi].mono_id << " idex :"
                     << arg->ions_n[pi].index << endl;
            }
            cout << "C端选取的离子为 : \n";
            for (int pi = 0; pi < arg->ions_c.size(); pi++) {
                cout << " T : " << arg->ions_c[pi].thoe_id << " M : " << arg->ions_c[pi].mono_id << " idex :"
                     << arg->ions_c[pi].index << endl;
            }
            //将部分结果写入文件，仅供测试使用
            WriteSequenceAlignment("../wtop_process/data/test-result-c-n.txt", arg->theory_ions_mass_n, spact->ions_mass_container,
                                   arg->alignment_ions_n, arg->sub_n, arg->filter_ions_n, arg->ppm_value_n,
                                   arg->ref_score_n, arg->peer_score_n);
            WriteSequenceAlignment("../wtop_process/data/test-result-n-c.txt", arg->theory_ions_mass_c, spact->ions_mass_container,
                                   arg->alignment_ions_c, arg->sub_c, arg->filter_ions_c, arg->ppm_value_c, arg->ref_score_c,
                                   arg->peer_score_c);
            cout << "score = " << arg->score << ", match Peaks = " << arg->peaks_match << " , all peaks = "
                 << spact->ions_mass_container.size() << endl;
        }


        void writeIonsMessage(msalign_processor_ptr msIt,
                              config_argument_ptr configArgumentPtr,
                              string msSrcPath)
        {
            // string srcPath = msSrcPath.substr(0,msSrcPath.find_last_of("/"));
            // srcPath += "\\PrSMs" ; //输出路径
            // string scans = msIt->Scans.substr(msIt->Scans.find_first_of("=")+1);    //Scams
            // string outFile = srcPath +"\\"+ scans + ".txt" ; //输出文件名
            // if (_access(srcPath.c_str(),0) == -1) {
            //     _mkdir(srcPath.c_str());
            // }
            // std::ofstream out(outFile,ios::out) ;

            // vector<node> nodeN = msIt->ionsMessageContainerN ;
            // vector<node> nodeC = msIt->ionsMessageContainerC ;
            // vector<double> theoryN = msIt->theoryMassN;
            // vector<double> theoryC = msIt->theoryMassC;

            // out <<msIt->Scans<<endl ;
            // out <<msIt->Activation<<endl ;
            // out<<"<<N Terminal Ions Begin>>"<<endl;
            // for (int i = 0; i < nodeN.size(); i++) {
            //     out << setiosflags(ios::fixed) << setprecision(2)
            //     <<setw(10) <<"theoryId : "<<setw(5)<< nodeN[i].thoe_id
            //     <<setw(10)<<" msId : "<<setw(5)<< nodeN[i].mono_id
            //     <<setw(10)<<" Theoretical mass :"<<setw(10)<< theoryN[nodeN[i].thoe_id]
            //     <<setw(10)<< " Mono mass :"<<setw(10)<< msIt->mono_mass[nodeN[i].mono_id]
            //     <<setw(10)<<" Sub :"<<setw(5)<< nodeN[i].index
            //     <<setw(10)<< " PPM :" <<setw(5)<< nodeN[i].ppm <<endl;
            // }
            // out<<"<<N Terminal Ions END>>"<<endl;
            // out<<"<<C Terminal Ions Begin>>"<<endl;
            // for (int i = 0; i < nodeC.size(); ++i) {
            //     out << setiosflags(ios::fixed) << setprecision(2)
            //         <<setw(10) <<"theoryId : "<<setw(5)<< nodeC[i].thoe_id
            //         <<setw(10)<<" msId : "<<setw(5)<< nodeC[i].mono_id
            //         <<setw(10)<<" Theoretical mass :"<<setw(10)<< theoryC[nodeC[i].thoe_id]
            //         <<setw(10)<< " Mono mass :"<<setw(10)<< msIt->mono_mass[nodeC[i].mono_id]
            //         <<setw(10)<<" Sub :"<<setw(5)<< nodeC[i].index
            //         <<setw(10)<< " PPM :" <<setw(5)<< nodeC[i].ppm <<endl;
            // }
            // out<<"<<C Terminal Ions End>>"<<endl;
            // out.close();
        }

    }
}
