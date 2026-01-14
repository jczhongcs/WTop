//
// Created by 碎雨粘霓裳 on 2021/12/7.
//

#include "Argument.h"

namespace mylib {


    //得到spectrumID,Scans,action,proteinAccession,proteinDescription
    void Argument::processID(string buff) {
        int idBegin = buff.find("#ID") + 4;     //id的起始位置
        int scansBegin = buff.find("#SCANS") + 7;
        int actionBegin = buff.find("#Types") + 9;
        string pro;
        // 获取蛋白信息
        if (buff.find(">sp") != buff.npos) {
            pro = buff.substr(buff.find(">sp"), buff.size() - buff.find(">sp") - 1);
        } else if (buff.find(">DECOY") != buff.npos) {
            pro = buff.substr(buff.find(">DECOY"), buff.size() - buff.find(">DECOY") - 1);
        } else if (buff.find(">tr") != buff.npos) {
            pro = buff.substr(buff.find(">tr"), buff.size() - buff.find(">tr") - 1);

        }
        int cutProLoc = pro.find_first_of(' '); //蛋白质前缀分界点
        //得到蛋白质前缀proteinAccession和蛋白质描述proteinDescription
        proteinAccession = pro.substr(0, cutProLoc);
        proteinDescription = pro.substr(cutProLoc + 1, buff.size() - cutProLoc - 2);
        for(int i = 0 ; i < proteinDescription.size();++i)
        {
            if (proteinDescription[i]==',') {
                proteinDescription[i] = ' ';
            }
        }
        //得到ID
        while (buff[idBegin] != ' ' || buff[idBegin] != ',') {
            this->spectrumID = this->spectrumID + buff[idBegin];
            idBegin++;
            if (buff[idBegin] == ' ' || buff[idBegin] == ',') {
                break;
            }
        }
        //得到scans
        while (buff[scansBegin] != ' ' || buff[scansBegin] != ',') {
            this->Scans = this->Scans + buff[scansBegin];
            scansBegin++;
            if (buff[scansBegin] == ' ' || buff[scansBegin] == ',') {
                break;
            }
        }
        //得到action
        while (buff[actionBegin] != ' ' || buff[actionBegin] != ',') {
            this->action = this->action + buff[actionBegin];
            this->fragmentation = action ;
            actionBegin++;
            if (buff[actionBegin] == ' ' || buff[actionBegin] == ',') {
                break;
            }
        }
    }


    //得到peaks和Score
    void Argument::processPeaksAndScore(string buff) {
        int peaksBegin = buff.find("Peaks") + 8;
        int scoreBegin = buff.find_last_of("=") + 2;
        //得到peaks
        while (buff[peaksBegin] != ' ' || buff[peaksBegin] != ',' && scoreBegin < buff.size()) {
            this->peaks = this->peaks + buff[peaksBegin];
            peaksBegin++;
            if (buff[peaksBegin] == ' ' || buff[peaksBegin] == ',') {
                break;
            }
        }
        //得到score
        this->score = buff.substr(scoreBegin);

        if (this->score.find('\n') != this->score.npos) {
            this->score.erase((std::remove(this->score.begin(), this->score.end(), '\n'),
                    this->score.end()));
        }

    }

    //得到mass信息
    void Argument::processMass(string buff) {
        int pBegin = buff.find_first_of("=") + 2;
        int aBegin = buff.find_last_of("=") + 2;
        //得到PreMass
        while (buff[pBegin] != ' ' || buff[pBegin] != ',' && pBegin < buff.size()) {
            this->precursorMass = this->precursorMass + buff[pBegin];
            pBegin++;
            if (buff[pBegin] == ' ' || buff[pBegin] == ',') {
                break;
            }
        }
        this->adjustedPrecursorMass = buff.substr(aBegin);
        this->adjustedPrecursorMass.pop_back();
        this->pPrecursorMass = atof(this->precursorMass.c_str());
        this->aPrecursorMass = atof(this->adjustedPrecursorMass.c_str());
    }

    //得到截断
    void Argument::processCut(string buff) {
        int equalFirst = buff.find_first_of("=") + 2;
        int equalLast = buff.find_last_of("=") + 2;
        while (buff[equalFirst] != ' ' || buff[equalFirst] != ',' && equalFirst < buff.size()) {
            this->firstResidue = this->firstResidue + buff[equalFirst];
            equalFirst++;
            if (buff[equalFirst] == ' ' || buff[equalFirst] == ',') {
                break;
            }
        }
        this->lastResidue = buff.substr(equalLast);
    }

    /**
     * 处理修饰
     * @param buff
     */
    void Argument::processVarible(string buff) {
        string pBuff;
        int dFlag = buff.find_last_of("#") + 2;   //定位MOD
        pBuff = buff.substr(dFlag);     //获取修饰
        variblePTMPbuff += pBuff + " & ";       //保存所有修饰
        variablePTMs += pBuff.length();      //取得修饰个数
    }

    void Argument::processFragmentIons(string basicString)
    {
        string pBuff ;
        pBuff = basicString.substr(basicString.find("=") + 2 );
        fragmentIonsSize = pBuff ;
    }

    void Argument::process_all_peaks(const string &buff)
    {
        int dFlag = buff.find_last_of("=") + 2 ;
        string pBuff = buff.substr(dFlag) ;
        all_peaks = atoi(pBuff.c_str());
    }

    void Argument::process_mut_x(const string &buff) {
        int log = buff.find("=")+2 ;
        string nub = buff.substr(log);
        this->mut_x = atof(nub.c_str()) ;
    }

    void Argument::processIonsPair(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->ionsPair = t ;
    }

    void Argument::processHuBuIonsCount(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->huBuIonsCount = t ;
    }
    void Argument::processSeqIonsCount(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->seqIonsCount = t ;
    }
    void Argument::processSeqScore(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->seqScore = t ;
    }

    void Argument::getProteoform (string basicString)
    {
        string s = basicString.substr(basicString.find_first_of("#")+2) ;
        this->proteoform = s ;
    }

    void Argument::processLastResidue(string basicString)
    {
        string str = basicString.substr(basicString.find("#")+2);
        int size = 0;
        for(int i = 0 ; i < str.size();++i) {
            if(str[i]!='.') {
                size++;
            }
        }
        this->lastResidueInt = (stoi(proteinLength) - stoi(this->lastResidue)) ;


        this->firstResidueInt = stoi(this->firstResidue) + 1;

    }
    void Argument::processProteinLength(string basicString){
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->proteinLength = t ;
    }

    ///调试打分  2023-8-14
    void Argument::processFeature1(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature1 = t ;
    }

    void Argument::processFeature2(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature2 = t ;
    }
    void Argument::processFeature3(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature3 = t ;
    }
    void Argument::processFeature4(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature4 = t ;
    }
    void Argument::processFeature5(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature5 = t ;
    }
    void Argument::processFeature6(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature6 = t ;
    }
    void Argument::processFeature7(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature7 = t ;
    }
    void Argument::processFeature8(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature8 = t ;
    }
    void Argument::processFeature9(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature9 = t ;
    }
    void Argument::processFeature10(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature10 = t ;
    }
    void Argument::processFeature11(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature11 = t ;
    }
    void Argument::processFeature12(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature12 = t ;
    }
    void Argument::processFeature13(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature13 = t ;
    }
    void Argument::processFeature14(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature14 = t ;
    }
    void Argument::processFeature15(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature15 = t ;
    }
    void Argument::processFeature16(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature16 = t ;
    }
    void Argument::processFeature17(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature17 = t ;
    }
    void Argument::processFeature18(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature18 = t ;
    }
    void Argument::processFeature19(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature19 = t ;
    }
    void Argument::processFeature20(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature20 = t ;
    }
    void Argument::processFeature21(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature21 = t ;
    }
    void Argument::processFeature22(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature22 = t ;
    }
    void Argument::processFeature23(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature23 = t ;
    }
    void Argument::processFeature24(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature24 = t ;
    }
    void Argument::processFeature25(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature25 = t ;
    }
    void Argument::processFeature26(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature26 = t ;
    }
    void Argument::processFeature27(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature27 = t ;
    }
    void Argument::processFeature28(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature28 = t ;
    }
    void Argument::processFeature29(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature29 = t ;
    }
    void Argument::processFeature30(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature30 = t ;
    }
    void Argument::processFeature31(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature31 = t ;
    }
    void Argument::processFeature32(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature32 = t ;
    }
    void Argument::processFeature33(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature33 = t ;
    }
    void Argument::processFeature34(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature34 = t ;
    }
    void Argument::processFeature35(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature35 = t ;
    }
    void Argument::processFeature36(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature36 = t ;
    }
    void Argument::processFeature37(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature37 = t ;
    }
    void Argument::processFeature38(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature38 = t ;
    }
    void Argument::processFeature39(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature39 = t ;
    }
    void Argument::processFeature40(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature40 = t ;
    }
    void Argument::processFeature41(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature41 = t ;
    }
    void Argument::processFeature42(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature42 = t ;
    }
    void Argument::processFeature43(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature43 = t ;
    }
    void Argument::processFeature44(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature44 = t ;
    }
    void Argument::processFeature45(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature45 = t ;
    }
    void Argument::processFeature46(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature46 = t ;
    }
    void Argument::processFeature47(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature47 = t ;
    }
    void Argument::processFeature48(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature48 = t ;
    }
    void Argument::processFeature49(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature49 = t ;
    }
    void Argument::processFeature50(string basicString) {
        int index = basicString.find("=") + 2;
        string t =  basicString.substr(index);
        this->feature50 = t ;
    }

    ///测试结束
}
