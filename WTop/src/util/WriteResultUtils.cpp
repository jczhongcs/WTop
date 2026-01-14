//
// Created by wenzhong on 2022/11/22.
//

#include "WriteResultUtils.h"

/**
 * 将g::Argument容器写入文件流中
 * @param out 输出流
 * @param aContainer 结果容器
 */
void WriteArgumentResult(ofstream &out, vector<mylib::Argument> &aContainer) {
    //输出表头
    out << "Spectrum ID\t"
        << "Fragmentation\t"
        << "Scan(s)\t"
        << "Fragment Ions Size\t"
        << "Match Peaks\t"
        << "All Peaks\t"
        << "M/A*100\t"
        << "Mut_X\t"
        << "Precursor mass\t"
        << "Adjusted precursor mass\t"
        << "Protein accession\t"
        << "Protein description\t"
        << "First residue\t"
        << "Last residue\t"
        << "#variablePTMs\t"
        << "Proteoform\t"
        << "Score\t"
        << "Fdr\t"
        << "Protein\t"
        << "Prsm\t"
        << "Choice\t"
        << "ionsPair\t"
        << "huBuIonsCount\t"
        << "seqIonsCount\t"
        << "seqScore\t"
        << "Protein Length\t"
        << "feature1\t"
        << "feature2\t"
        << "feature3\t"
        << "feature4\t"
        << "feature5\t"
        << "feature6\t"
        << "feature7\t"
        << "feature8\t"
        << "feature9\t"
        << "feature10\t"
        << "feature11\t"
        << "feature12\t"
        << "feature13\t"
        << "feature14\t"
        << "feature15\t"
        << "feature16\t"
        << "feature17\t"
        << "feature18\t"
        << "feature19\t"
        << "feature20\t"
        << "feature21\t"
        << "feature22\t"
        << "feature23\t"
        << "feature24\t"
        << "feature25\t"
        << "feature26\t"
        << "feature27\t"
        << "feature28\t"
        << "feature29\t"
        << "feature30\t"
        << "feature31\t"
        << "feature32\t"
        << "feature33\t"
        << "feature34\t"
        << "feature35\t"
        << "feature36\t"
        << "feature37\t"
        << "feature38\t"
        << "feature39\t"
        << "feature40\t"
        << "feature41\t"
        << "feature42\t"
        << "feature43\t"
        << "feature44\t"
        << "feature45\t"
        << "feature46\t"
        << "feature47\t"
        << "feature48\t"
        << "feature49\t"
        << "feature50\t"
        << endl;

    vector<mylib::Argument>::iterator it = aContainer.begin();
    //输出每条结果
    for (; it != aContainer.end(); it++) {
        out << it->spectrumID
            << "\t" << it->action
            << "\t" << it->Scans
            << "\t" << it->fragmentIonsSize
            << "\t" << it->peaks
            << "\t" << it->all_peaks
            << "\t" << (double) (atof(it->peaks.c_str()) / (double) it->all_peaks) * 100
            << "\t" << it->mut_x
            << "\t" << it->precursorMass
            << "\t" << it->adjustedPrecursorMass
            << "\t" << it->proteinAccession
            << "\t" << it->proteinDescription
            << "\t" << it->firstResidue
            << "\t" << it->lastResidue
            << "\t" << it->variblePTMPbuff
            << "\t" << it->proteoform
            << "\t" << it->score
            << "\t" << it->fdrScore
            << "\t" << it->proteinAccession.substr(1)
            << "\t" << it->Scans + it->proteinAccession.substr(1)
            << "\t" << it->chiose
            << "\t" << it->ionsPair
            << "\t" << it->huBuIonsCount
            << "\t" << it->seqIonsCount
            << "\t" << it->seqScore
            << "\t" << it->proteinLength
            << "\t" << it->feature1
            << "\t" << it->feature2
            << "\t" << it->feature3
            << "\t" << it->feature4
            << "\t" << it->feature5
            << "\t" << it->feature6
            << "\t" << it->feature7
            << "\t" << it->feature8
            << "\t" << it->feature9
            << "\t" << it->feature10
            << "\t" << it->feature11
            << "\t" << it->feature12
            << "\t" << it->feature13
            << "\t" << it->feature14
            << "\t" << it->feature15
            << "\t" << it->feature16
            << "\t" << it->feature17
            << "\t" << it->feature18
            << "\t" << it->feature19
            << "\t" << it->feature20
            << "\t" << it->feature21
            << "\t" << it->feature22
            << "\t" << it->feature23
            << "\t" << it->feature24
            << "\t" << it->feature25
            << "\t" << it->feature26
            << "\t" << it->feature27
            << "\t" << it->feature28
            << "\t" << it->feature29
            << "\t" << it->feature30
            << "\t" << it->feature31
            << "\t" << it->feature32
            << "\t" << it->feature33
            << "\t" << it->feature34
            << "\t" << it->feature35
            << "\t" << it->feature36
            << "\t" << it->feature37
            << "\t" << it->feature38
            << "\t" << it->feature39
            << "\t" << it->feature40
            << "\t" << it->feature41
            << "\t" << it->feature42
            << "\t" << it->feature43
            << "\t" << it->feature44
            << "\t" << it->feature45
            << "\t" << it->feature46
            << "\t" << it->feature47
            << "\t" << it->feature48
            << "\t" << it->feature49
            << "\t" << it->feature50

            << endl;
    }
}

//Data file name
//Prsm ID\t
//Spectrum ID\t
//Fragmentation\t
//Scan(s)\t
//Retention time\t
//#peaks	\t
//Charge\t
//Precursor mass\t
//Adjusted precursor mass\t
//Proteoform ID\t
//Feature intensity\t
//Protein accession\t
//Protein description\t
//First residue\t
//Last residue\t
//Proteoform\t
//#unexpected modifications\t
//#variable PTMs\t
//#matched peaks\t
//#matched fragment ions\t
//P-value	E-value	Q-value (spectral FDR)	Proteoform FDR

void processResult(ofstream &out, vector<mylib::Argument> &aContainer)
{
    //输出表头
    out<< "Data file name\t" << "Prsm ID\t"<<"Spectrum ID\t"<<"Fragmentation\t"<<"Scan(s)\t"<<"Retention time\t"
       <<"#peaks\t"<<"Charge\t"<<"Precursor mass\t"<<"Adjusted precursor mass\t"<<"Proteoform ID\t"<<"Feature intensity\t"
       <<"Protein accession\t"<<"Protein description\t"<<"First residue\t"<<"Last residue\t"<<"Proteoform\t"
       <<"#unexpected modifications\t"<<"#variable PTMs\t"<<"#matched peaks\t"<<"#matched fragment ions\t"
       <<"Score\t"<<"FDR\t"<<endl;
    vector<mylib::Argument>::iterator it = aContainer.begin();
    //输出每条结果
    int i = 0 ;
    for (vector<mylib::Argument>::iterator it = aContainer.begin(); it != aContainer.end(); it++)
    {
        out <<"-"
            << "\t" << it->prsmId
            << "\t" << it->spectrumID
            << "\t" << it->fragmentation
            << "\t" << it->Scans
            << "\t" << it->RetentionTime
            << "\t" << it->all_peaks
            << "\t" << it->charge
            << "\t" << it->precursorMass
            << "\t" << it->adjustedPrecursorMass
            << "\t" << it->ProteoformID
            << "\t" << it->FeatureIntensity
            << "\t" << it->proteinAccession.substr(it->proteinAccession.find_first_of(">")+1)
            << "\t" << it->proteinDescription
            << "\t" << it->firstResidue
            << "\t" << it->lastResidue
            << "\t" << it->proteoform
            << "\t" << it->unexpectedModifications
            << "\t" << it->variablePTMsString
            << "\t" << it->peaks
            << "\t" << it->fragmentIonsSize
            << "\t" << it->score
            << "\t" << it->fdrScore
            <<endl;
    }
}

//Data file name,Prsm ID,Spectrum ID,Fragmentation,Scan(s),Retention time,
// #peaks,Charge,Precursor mass,Adjusted precursor mass,Proteoform ID,Feature intensity,
// Protein accession,Protein description,First residue,Last residue,Proteoform,
// #unexpected modifications,MIScore,#variable PTMs,#matched peaks,#matched fragment ions
// ,P-value,E-value,Q-value (spectral FDR),Proteoform FDR
//"Data file name,Prsm ID,Spectrum ID,Fragmentation,Scan(s),Retention time,
// #peaks,Charge,Precursor mass,Adjusted precursor mass,Proteoform ID,Feature intensity,
// Protein accession,Protein description,First residue,Last residue,Proteoform,
// #unexpected modifications,MIScore,#variable PTMs,#matched peaks,#matched fragment ions,
// P-value,E-value,Q-value (spectral FDR),Proteoform FDR",

void processResultCSV(ofstream &out, vector<mylib::Argument> &aContainer, map<string,string> & ssMap)
{
    for(map<string,string>::iterator it = ssMap.begin(); it!= ssMap.end(); ++it)
    {
        out <<it->first<<","<<it->second<<endl;
    }
    //输出表头
    out<< "Data file name," << "Prsm ID,"<<"Spectrum ID,"<<"Fragmentation,"<<"Scan(s),"<<"Retention time,"
       <<"#peaks,"<<"Charge,"<<"Precursor mass,"<<"Adjusted precursor mass,"<<"Proteoform ID,"<<"Feature intensity,"
       <<"Protein accession,"<<"Protein description,"<<"First residue,"<<"Last residue,"<<"Proteoform,"
       <<"#unexpected modifications,"<<"MIScore,"<<"#variable PTMs,"<<"#matched peaks,"<<"#matched fragment ions,"
       <<"P-value,"<<"E-value,"<<"Q-value (spectral FDR),"<<"Proteoform FDR"<<endl;
    vector<mylib::Argument>::iterator it = aContainer.begin();
    //输出每条结果
    int i = 0 ;
    for (vector<mylib::Argument>::iterator it = aContainer.begin(); it != aContainer.end(); it++)
    {
        out <<"-"
            << "," << it->prsmId
            << "," << it->spectrumID
            << "," << it->fragmentation
            << "," << it->Scans
            << "," << it->RetentionTime
            << "," << it->all_peaks
            << "," << it->charge
            << "," << it->precursorMass
            << "," << it->adjustedPrecursorMass
            << "," << ++i
            << "," << it->FeatureIntensity
            << "," << it->proteinAccession.substr(it->proteinAccession.find_first_of(">")+1)
            << "," << it->proteinDescription
            << "," << it->firstResidueInt
            << "," << it->lastResidueInt
            << "," << it->proteoform
            << "," << it->unexpectedModifications
            <<","<<"-"
            << "," << it->variablePTMsString
            << "," << it->peaks
            << "," << it->fragmentIonsSize
            <<","<<"-"
            <<","<<"-"
            << "," << it->score
            << "," << it->fdrScore
            <<endl;
    }
}


void processProteoform(vector<mylib::Argument> &bContainer)
{

}
/**
 * 将txt文本处理为xlsx文件
 * @param inFileName
 * @param outFileName
 * @return
 */
vector<mylib::Argument> processOneFileResult(string inFileName , string outFileName , map<string,string> & ssMap)
{

    //输出文件名字
    //prsm输出文件名字
//    string outFileName = "2DLC_H2A_ms2_4_15.xlsx";
    int prsmMinNub = 0;    //prsm > 1
    double fdr = 0.01;     //fdr < 0.01
    double filterArg = 0; //Prsm = 1 的占比小于Arg则过滤
    string preName = outFileName.substr(0, outFileName.find_last_of('.'));
    //proteinform输出文件名

    string outPrFormFilename = preName + "_Raw_ProteinForm.xlsx";
    string outAllFieName = preName + "_Raw_Decoy_PrSMs.xlsx" ;

    string outPrSMsTop = preName+"_WTop_PrSM.csv";       
    string outProtenFormTop = preName+"_WTop_ProteinForm.csv";

    const char *cinFileName = inFileName.c_str();
    const char *coutFileName = outFileName.c_str();
    const char *coutPrFileName = outPrFormFilename.c_str();
    const char *coutAllFileName = outAllFieName.c_str();

    const char *coutPrSMsTop = outPrSMsTop.c_str();
    const char *coutProtenFormsTop = outProtenFormTop.c_str();

    ifstream in(cinFileName);    // 输入文件名`
//    remove(coutFileName);
//    ofstream out(coutFileName);   // 输出Prsm文件名
//    ofstream out_prForm(coutPrFileName);   // 输出ProteinForm文件名
//    ofstream out_all(coutAllFileName);   // 输出ProteinForm文件名
//
//    ofstream outTop(coutPrSMsTop);
//    ofstream outTopProtenForm(coutProtenFormsTop);
    if (!in)
    {
        cout << "in file name is null " << endl;
    }
    cout << "in file name : " << inFileName << "\nout file name : " << outFileName << endl;
    string buff;   // 用于接受读入的字符串
    bool loop = false; // 用于判断该段的内容为一条质谱的内容输入
    vector<mylib::Argument> aContainer; //所有结果
    vector<mylib::Argument> bContainer; //-fdr-decoy
    vector<mylib::Argument> cContainer; //-fdr-decoy
    int prsmId = 1 ;
    while (in.good()) {      //用于判断文件是否读取完毕
        buff.clear();
        getline(in, buff);          //输入该行内容
        while (buff.find("#BEGIN") != buff.npos) {    // 结果标识位
            loop = true;
            mylib::Argument p;          //一条结果的存储对象
            do {                    //循环读取结果
                buff.clear();       //清空缓存
                getline(in, buff);
                if (buff.find("fragment") != buff.npos)
                {
                    p.processFragmentIons(buff);
                }
                //得到spectrumID,Scans,action,proteinAccession,proteinDescription
                if (buff.find("#ID") != buff.npos) {
                    p.processID(buff);
                }
                //得到peaks和Score
                if (buff.find("Peaks") != buff.npos) {
                    p.processPeaksAndScore(buff);
                }
                //得到截断信息
                if (buff.find("CutLocationN") != buff.npos) {
                    p.processCut(buff);
                }
                //得到前缀质量
                if (buff.find("Precursor mass") != buff.npos) {
                    p.processMass(buff);
                }
                //得到variblePTMs
                if (buff.find("Mod#") != buff.npos) {
                    p.processVarible(buff);
                }
                if (buff.find("All") != buff.npos && buff.find("peaks") != buff.npos) {
                    p.process_all_peaks(buff);
                }
                if (buff.find("mut_x") != buff.npos ) {
                    p.process_mut_x(buff);
                }
                if (buff.find("ionsPair") != buff.npos) {
                    p.processIonsPair(buff);
                }
                if (buff.find("huBuIonsCount") != buff.npos) {
                    p.processHuBuIonsCount(buff);
                }
                if (buff.find("seqIonsCount") != buff.npos) {
                    p.processSeqIonsCount(buff);
                }
                if (buff.find("seqScore") != buff.npos) {
                    p.processSeqScore(buff);
                }
                if (buff.find("ProteoForm#")!=buff.npos) {
                    p.getProteoform(buff) ;
                }
                if (buff.find("RawProteinSeq#")!=buff.npos) {
                    p.processLastResidue(buff);
                }
                if (buff.find("True ModifyMass") != buff.npos) {
                    p.getTrueModMass(buff);
                }
                if (buff.find("Protein length") != buff.npos) {
                    p.processProteinLength(buff);
                }


                ///调试打分
                if (buff.find("feature01") != buff.npos) {
                    p.processFeature1(buff);
                }
                if (buff.find("feature02") != buff.npos) {
                    p.processFeature2(buff);
                }
                if (buff.find("feature03") != buff.npos) {
                    p.processFeature3(buff);
                }
                if (buff.find("feature04") != buff.npos) {
                    p.processFeature4(buff);
                }
                if (buff.find("feature05") != buff.npos) {
                    p.processFeature5(buff);
                }
                if (buff.find("feature06") != buff.npos) {
                    p.processFeature6(buff);
                }
                if (buff.find("feature07") != buff.npos) {
                    p.processFeature7(buff);
                }
                if (buff.find("feature08") != buff.npos) {
                    p.processFeature8(buff);
                }
                if (buff.find("feature09") != buff.npos) {
                    p.processFeature9(buff);
                }
                if (buff.find("feature10") != buff.npos) {
                    p.processFeature10(buff);
                }
                if (buff.find("feature11") != buff.npos) {
                    p.processFeature11(buff);
                }
                if (buff.find("feature12") != buff.npos) {
                    p.processFeature12(buff);
                }
                if (buff.find("feature13") != buff.npos) {
                    p.processFeature13(buff);
                }
                if (buff.find("feature14") != buff.npos) {
                    p.processFeature14(buff);
                }
                if (buff.find("feature15") != buff.npos) {
                    p.processFeature15(buff);
                }
                if (buff.find("feature16") != buff.npos) {
                    p.processFeature16(buff);
                }
                if (buff.find("feature17") != buff.npos) {
                    p.processFeature17(buff);
                }
                if (buff.find("feature18") != buff.npos) {
                    p.processFeature18(buff);
                }
                if (buff.find("feature19") != buff.npos) {
                    p.processFeature19(buff);
                }
                if (buff.find("feature20") != buff.npos) {
                    p.processFeature20(buff);
                }
                if (buff.find("feature21") != buff.npos) {
                    p.processFeature21(buff);
                }
                if (buff.find("feature22") != buff.npos) {
                    p.processFeature22(buff);
                }
                if (buff.find("feature23") != buff.npos) {
                    p.processFeature23(buff);
                }
                if (buff.find("feature24") != buff.npos) {
                    p.processFeature24(buff);
                }
                if (buff.find("feature25") != buff.npos) {
                    p.processFeature25(buff);
                }
                if (buff.find("feature26") != buff.npos) {
                    p.processFeature26(buff);
                }
                if (buff.find("feature27") != buff.npos) {
                    p.processFeature27(buff);
                }
                if (buff.find("feature28") != buff.npos) {
                    p.processFeature28(buff);
                }
                if (buff.find("feature29") != buff.npos) {
                    p.processFeature29(buff);
                }
                if (buff.find("feature30") != buff.npos) {
                    p.processFeature30(buff);
                }
                if (buff.find("feature31") != buff.npos) {
                    p.processFeature31(buff);
                }
                if (buff.find("feature32") != buff.npos) {
                    p.processFeature32(buff);
                }
                if (buff.find("feature33") != buff.npos) {
                    p.processFeature33(buff);
                }
                if (buff.find("feature34") != buff.npos) {
                    p.processFeature34(buff);
                }
                if (buff.find("feature35") != buff.npos) {
                    p.processFeature35(buff);
                }
                if (buff.find("feature36") != buff.npos) {
                    p.processFeature36(buff);
                }
                if (buff.find("feature37") != buff.npos) {
                    p.processFeature37(buff);
                }
                if (buff.find("feature38") != buff.npos) {
                    p.processFeature38(buff);
                }
                if (buff.find("feature39") != buff.npos) {
                    p.processFeature39(buff);
                }
                if (buff.find("feature40") != buff.npos) {
                    p.processFeature40(buff);
                }
                if (buff.find("feature41") != buff.npos) {
                    p.processFeature41(buff);
                }
                if (buff.find("feature42") != buff.npos) {
                    p.processFeature42(buff);
                }
                if (buff.find("feature43") != buff.npos) {
                    p.processFeature43(buff);
                }
                if (buff.find("feature44") != buff.npos) {
                    p.processFeature44(buff);
                }
                if (buff.find("feature45") != buff.npos) {
                    p.processFeature45(buff);
                }
                if (buff.find("feature46") != buff.npos) {
                    p.processFeature46(buff);
                }
                if (buff.find("feature47") != buff.npos) {
                    p.processFeature47(buff);
                }
                if (buff.find("feature48") != buff.npos) {
                    p.processFeature48(buff);
                }
                if (buff.find("feature49") != buff.npos) {
                    p.processFeature49(buff);
                }
                if (buff.find("feature50") != buff.npos) {
                    p.processFeature50(buff);
                }
                ///调试结束
                if (buff.find("Unknown_PTM") != std::string::npos) {
                    p.getUnknwonPtmMass(buff);
                } else if (buff.find("PTM_Mass") != std::string::npos) {
                    p.getPTMMass(buff);
                }
                if (buff.find("#END") != buff.npos) {
                    loop = false;
                }
            } while (loop);
            p.prsmId = prsmId++ ;
            aContainer.push_back(p);        //结果容器接受结果
        }
    }

    cout << "End of processing out file ..." << endl;

    //原始数据Score排序
    mylib::Service::ScoreSortArgument(aContainer);

    //排序后计算每条的fdr
    int numberFdr = mylib::Service::calculateFdr(aContainer);   //得到fdr 并返回fdr的个数
//    cout<<"Decoy的个数为: "<<numberFdr<<endl ;
    //去掉fdr和Decoy并且去掉只有1条Prsm的占比低于filterArg -> bContainer
    mylib::Service::fdrChoice(aContainer, bContainer, fdr, prsmMinNub, filterArg);

    //根据质量排序
    mylib::Service::MassSortArgument(bContainer);
    //得到最终的蛋白质形
    if (inFileName.find("unknown")!=std::string::npos) {
        mylib::Service::getUnknownProteinForm(bContainer, cContainer);
    } else {
        mylib::Service::getProteinForm(bContainer, cContainer);
    }



    ofstream outTop(coutPrSMsTop);
    ofstream outTopProtenForm(coutProtenFormsTop);

     ofstream out(coutFileName);   // 输出Prsm文件名
     ofstream out_prForm(coutPrFileName);   // 输出ProteinForm文件名
     ofstream out_all(coutAllFileName);   // 输出ProteinForm文件名

     //输出xlsx
     WriteArgumentResult(out, bContainer);   //fdr 之后的prsm
     out.close();
     WriteArgumentResult(out_prForm, cContainer); //proteoForm->cContainer
     out_prForm.close();
     WriteArgumentResult(out_all, aContainer);   //输出所有 prsm
     out_all.close();

    //输出topmg表格形式
    processResultCSV(outTop,bContainer,ssMap);
    processResultCSV(outTopProtenForm,cContainer,ssMap);
    in.close();
//    out.close();
    return aContainer;
}
