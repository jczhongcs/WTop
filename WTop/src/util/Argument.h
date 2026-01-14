//
// Created by wenzhong on 2023/1/5.
//

#ifndef DTW_WTOP_ARGUMENT_H
#define DTW_WTOP_ARGUMENT_H


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

/*
 * 对结果文件进行操作
 * 1、文件对象，并获取对象信息
 * */

namespace mylib {
    class Argument {
    public:
        int prsmId ; //prsmId
        string spectrumID;      //spectrumId
        string fragmentation ; //CID or ETD
        string Scans;           //Scans(s)
        string RetentionTime = "-" ; //Retention time
        string action;   //(Action)
        int all_peaks  = 0 ;    //#peaks
        string charge="-" ; //charge
        string precursorMass;   //Precursor mass
        double pPrecursorMass = 0; // double pre
        string adjustedPrecursorMass;   // Adjusted precursor mass
        double aPrecursorMass = 0;
        string ProteoformID ="-"; //Proteoform ID
        string FeatureIntensity= "-";//Feature intensity
        string proteinAccession;        //Protein accession
        string proteinDescription;      //Protein description
        string firstResidue;            //First residue
        int firstResidueInt;
        string lastResidue;             //Last residue
        int lastResidueInt ;
        string proteoform;              //Proteoform
        string unexpectedModifications="0"; //#unexpected modifications = "0";
        string MIScore = "-";
        string variablePTMsString = "0" ;
        int matchPeaks ;      //#matched peaks peaks
        string fragmentIonsSize ; //#matched fragment ions
        string P_value = "-";
        string E_value = "-" ;

        string ptmMass;
        string unknownPtmMass ; 
        double trueModMass = 0 ;

        string peaks;

        int variablePTMs = 0;
        int matchedPeaks = 0;
        string score = "";
        double score_db = 0 ;
        string variblePTMPbuff;
        double fdrScore = 0 ;

        //得到spectrumID,Scans,action,proteinAccession,proteinDescription
        double mut_x ;
        string chiose  ;  //

        //2022-9-22添加
        string ionsPair ;
        string huBuIonsCount ;
        string seqIonsCount ;
        string seqScore ;

        string proteinLength;
        ///调试打分特征
        string feature1;
        string feature2;
        string feature3;
        string feature4;
        string feature5;
        string feature6;
        string feature7;
        string feature8;
        string feature9;
        string feature10;
        string feature11;
        string feature12;
        string feature13;
        string feature14;
        string feature15;
        string feature16;
        string feature17;
        string feature18;
        string feature19;
        string feature20;
        string feature21;
        string feature22;
        string feature23;
        string feature24;
        string feature25;
        string feature26;
        string feature27;
        string feature28;
        string feature29;
        string feature30;
        string feature31;
        string feature32;
        string feature33;
        string feature34;
        string feature35;
        string feature36;
        string feature37;
        string feature38;
        string feature39;
        string feature40;
        string feature41;
        string feature42;
        string feature43;
        string feature44;
        string feature45;
        string feature46;
        string feature47;
        string feature48;
        string feature49;
        string feature50;

        void processID(string buff);

        //得到peaks和Score
        void processPeaksAndScore(string buff) ;

        //得到mass信息
        void processMass(string buff);

        //得到截断
        void processCut(string buff);

        void processVarible(string buff) ;

        void process_all_peaks(const string &buff);

        void process_mut_x(const string &buff) ;

        map<string,string> getScansScoreMap();

        void processIonsPair(string basicString);
        void processHuBuIonsCount(string basicString);
        void processSeqIonsCount(string basicString);
        void processSeqScore(string basicString);

        void processFragmentIons(string basicString);

        void getProteoform (string basicString);

        void processLastResidue(string basicString);

        void processProteinLength(string basicString);

        ///调试打分
        void processFeature1(string basicString);
        void processFeature2(string basicString);
        void processFeature3(string basicString);
        void processFeature4(string basicString);
        void processFeature5(string basicString);
        void processFeature6(string basicString);
        void processFeature7(string basicString);
        void processFeature8(string basicString);
        void processFeature9(string basicString);
        void processFeature10(string basicString);
        void processFeature11(string basicString);
        void processFeature12(string basicString);
        void processFeature13(string basicString);
        void processFeature14(string basicString);
        void processFeature15(string basicString);
        void processFeature16(string basicString);
        void processFeature17(string basicString);
        void processFeature18(string basicString);
        void processFeature19(string basicString);
        void processFeature20(string basicString);
        void processFeature21(string basicString);
        void processFeature22(string basicString);
        void processFeature23(string basicString);
        void processFeature24(string basicString);
        void processFeature25(string basicString);
        void processFeature26(string basicString);
        void processFeature27(string basicString);
        void processFeature28(string basicString);
        void processFeature29(string basicString);
        void processFeature30(string basicString);
        void processFeature31(string basicString);
        void processFeature32(string basicString);
        void processFeature33(string basicString);
        void processFeature34(string basicString);
        void processFeature35(string basicString);
        void processFeature36(string basicString);
        void processFeature37(string basicString);
        void processFeature38(string basicString);
        void processFeature39(string basicString);
        void processFeature40(string basicString);
        void processFeature41(string basicString);
        void processFeature42(string basicString);
        void processFeature43(string basicString);
        void processFeature44(string basicString);
        void processFeature45(string basicString);
        void processFeature46(string basicString);
        void processFeature47(string basicString);
        void processFeature48(string basicString);
        void processFeature49(string basicString);
        void processFeature50(string basicString);


        void getTrueModMass(string basicString)
        {
            string mass = basicString.substr(basicString.find("=")+2);
            trueModMass = std::stod(mass);
        }

        void getPTMMass(string & ptmString)
        {
            ptmMass = ptmString.substr(ptmString.find("=")+2) ; 
        }

        void getUnknwonPtmMass(string & unknownPTMString)
        {
            unknownPtmMass = unknownPTMString.substr(unknownPTMString.find("=")+2);
        }

    };
}

#endif //DTW_WTOP_ARGUMENT_H
