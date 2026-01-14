//
// Created by wenzhong on 2022/11/15.
//

#include "string_utils.h"
#include <cstring>
//map<char,string> ptmsMap = {
//        {'A', "Acetyl"},
//        {'M', "Methyl"},
//        {'D', "Dimethyl"},
//        {'P', "Phospho"},
//};
namespace utils_string_util{

    void Stringsplit(string str, const char split,vector<string>& rst)
    {
        istringstream iss(str);	// 输入流
        string token;			// 接收缓冲区
        while (getline(iss, token, split))	// 以split为分隔符
        {
            rst.push_back(token) ;
        }
    }

    vector<string> split(const string& str, const string& delim) {  
        vector<string> res;  
        if ("" == str) return res;  
        //先将要切割的字符串从string类型转换为char*类型  
        char * strs = new char[str.length() + 1];  
        strcpy(strs, str.c_str());  
    
        char * d = new char[delim.length() + 1];  
        strcpy(d, delim.c_str());  
    
        char *p = strtok(strs, d);  
        while (p) {  
            string s = p;  //分割得到的字符串转换为string类型  
            res.push_back(s); //存入结果数组  
            p = strtok(NULL, d);  
        }  
        return res;  
    }  





    void processProteinFom(string &rawSeq,vector<string> & ptms,map<char,string> & ptmsMap,map<int,int> &startEnd,string &pfSeq)
    {
        ptmsMap.insert(make_pair('*',"UNKNOWN"));
        int firstPos = rawSeq.find_first_of(".") ;
        int lastPos = rawSeq.find_last_of(".");
        int posLen = lastPos - firstPos - 1 ;
        string preStr = rawSeq.substr(0,firstPos+1);
        string backStr = rawSeq.substr(lastPos);
        string formSeq =rawSeq.substr(rawSeq.find_first_of(".")+1,posLen);
        set<char> s = {'(',')','[',']'};
        string leftChar = "(";
        string rightChar = ")";
        string leftChar2 = "[";
        string rightChar2 = "]";
        map<int,char> charInt;  //key - seqChar , value - loc
        int ptmsLog = 0 ;
        //解析ptms
        for(int i = 0; i < ptms.size(); ++i) {
            string tempStr = ptms[i] ;
            string resultStr = "";
            int sSize = 0;
            while (sSize != tempStr.size()) {
                char c = tempStr[sSize];
                if (sSize == tempStr.size()-1 && ptmsMap.count(c)) {
                    resultStr += ptmsMap[c];
                } else if (ptmsMap.count(c)) {
                    resultStr += ptmsMap[c] + ";";
                }
                sSize++;
            }
            ptms[i] = resultStr ;
        }
        //计算蛋白质形
        for (map<int,int>::iterator it = startEnd.begin(); it != startEnd.end(); ++it) {
            int left = it->first + 1;
            int right = it->second;
//        if (it->first == 0) {
//             left = it->first + 1 ;
//             right = it->second + 1;
//        } else if (it->first == it->second){
//             left = it->first + 1 ;
//             right = it->second + 1;
//        } else {
//            left = it->first + 2 ;
//            right = it->second + 1;
//        }
//        cout <<left << " --- "<< right <<endl;
            int curCharLoc = -1 ;
            for (int i = 0 ; i < formSeq.size();++i) {
                if (s.count(formSeq[i])) {
                    continue;
                } else {
                    curCharLoc++;
                    if (curCharLoc == left) {
                        formSeq.insert(i,leftChar);
                        ++i;
                    }
                    if (curCharLoc == right) {
                        formSeq.insert(i+1,rightChar);
                        i++;
                        formSeq.insert(i+1,leftChar2);
                        i++;
                        formSeq.insert(i+1,rightChar2);
                        i++;
                        break ;
                    }
                }
            }
        }
        //插入修饰
        for (int i = 0; i < ptms.size(); ++i) {
            if (formSeq.find("[]") != formSeq.npos) {
                int loc = formSeq.find("[]") + 1 ;
                formSeq.insert(loc,ptms[i]);
            }
        }
        //前缀+中缀+后缀
        pfSeq = preStr + formSeq + backStr ;
    }

    void processProteinFomUnkMass(string &rawSeq,vector<string> & ptms,map<char,string> ptmsMap,map<int,int> &startEnd,string &pfSeq,string unk_mass)
    {   ptmsMap['*'] = unk_mass;
        ptmsMap.insert(make_pair('*',"UNKNOWN"));
        int firstPos = rawSeq.find_first_of(".") ;
        int lastPos = rawSeq.find_last_of(".");
        int posLen = lastPos - firstPos - 1 ;
        string preStr = rawSeq.substr(0,firstPos+1);
        string backStr = rawSeq.substr(lastPos);
        string formSeq =rawSeq.substr(rawSeq.find_first_of(".")+1,posLen);
        set<char> s = {'(',')','[',']'};
        string leftChar = "(";
        string rightChar = ")";
        string leftChar2 = "[";
        string rightChar2 = "]";
        map<int,char> charInt;  //key - seqChar , value - loc
        int ptmsLog = 0 ;
        //解析ptms
        for(int i = 0; i < ptms.size(); ++i) {
            string tempStr = ptms[i] ;
            string resultStr = "";
            int sSize = 0;
            while (sSize != tempStr.size()) {
                char c = tempStr[sSize];
                if (sSize == tempStr.size()-1 && ptmsMap.count(c)) {
                    resultStr += ptmsMap[c];
                } else if (ptmsMap.count(c)) {
                    resultStr += ptmsMap[c] + ";";
                }
                sSize++;
            }
            ptms[i] = resultStr ;
        }
        //计算蛋白质形
        for (map<int,int>::iterator it = startEnd.begin(); it != startEnd.end(); ++it) {
            int left = it->first + 1;
            int right = it->second;
//        if (it->first == 0) {
//             left = it->first + 1 ;
//             right = it->second + 1;
//        } else if (it->first == it->second){
//             left = it->first + 1 ;
//             right = it->second + 1;
//        } else {
//            left = it->first + 2 ;
//            right = it->second + 1;
//        }
//        cout <<left << " --- "<< right <<endl;
            int curCharLoc = -1 ;
            for (int i = 0 ; i < formSeq.size();++i) {
                if (s.count(formSeq[i])) {
                    continue;
                } else {
                    curCharLoc++;
                    if (curCharLoc == left) {
                        formSeq.insert(i,leftChar);
                        ++i;
                    }
                    if (curCharLoc == right) {
                        formSeq.insert(i+1,rightChar);
                        i++;
                        formSeq.insert(i+1,leftChar2);
                        i++;
                        formSeq.insert(i+1,rightChar2);
                        i++;
                        break ;
                    }
                }
            }
        }
        //插入修饰
        for (int i = 0; i < ptms.size(); ++i) {
            if (formSeq.find("[]") != formSeq.npos) {
                int loc = formSeq.find("[]") + 1 ;
                formSeq.insert(loc,ptms[i]);
            }
        }
        //前缀+中缀+后缀
        pfSeq = preStr + formSeq + backStr ;

    }

    void processCutLocation(string &seq,int cutN,int cutC)
    {
        if (cutN != 0) {
            seq.insert(cutN, ".");
            ///避免变体序列过长，溢出excel单元格  test:2023-08-21
            seq = seq.substr(cutN);
        } else {
            seq.insert(0, ".");
        }
        if (cutC != 0) {
            seq.insert(seq.size() - cutC, ".");
            seq = seq.substr(0,seq.size() - cutC);
        } else {
            seq.insert(seq.size() , ".");
        }
    }


    vector<mylib::Argument> conbineResult(vector<mylib::Argument> & known, vector<mylib::Argument> &unknown)
    {
        cout <<"开始合并未知修饰和已知修饰文件 其中已知prsm size = "<< known.size()
             <<"  未知prsm size = " << unknown.size()<<endl;
        map<string,mylib::Argument> kScansArg ;
        map<string,mylib::Argument> UnkScansArg ;
        map<string,string> knownScansScore = mylib::Service::getScansScoreMap(known, kScansArg);
        map<string,string> UnknownScansScore = mylib::Service::getScansScoreMap(unknown, UnkScansArg);
        vector<mylib::Argument> result ;
        cout<< "已知PTMs prsm 个数 : "<< knownScansScore.size()
            <<" , 未知PTMs prsm 个数 : " << UnknownScansScore.size()
            <<endl;
        return result;
    }

    void txtFileToCsvFile(vector<string> & outFile, map<string,string> & ssMap)
    {
        vector<string> csvOutFileName ;
        for (string inFileName : outFile)
        {
            string outFileName = inFileName.substr(0,inFileName.find_last_of(".") ) + ".xlsx";

            if (inFileName.find("unknown") != inFileName.npos)
            {
                //处理未知修饰结果
                vector<mylib::Argument> allUnKnownPrsm  = processOneFileResult(inFileName, outFileName, ssMap);
            } else {
                vector<mylib::Argument> allKnownPrsm = processOneFileResult(inFileName, outFileName, ssMap);

            }

            //将allKnown和allUnknown合并
//    vector<g::Argument> conbinePrsm = conbineResult(allKnownPrsm,allUnKnownPrsm);
//    cout << "合并后的prsm size = " << conbinePrsm.size()<<endl;
//    s::Service::ScoreSortArgument(conbinePrsm);
//
//    numberFdr = s::Service::calculateFdr(conbinePrsm);
//
//    vector<g::Argument> bContainer ;
//
//    s::Service::fdrChoice(conbinePrsm, bContainer, fdr, prsmMinNub, filterArg);

//    WriteArgumentResult(outConbine, conbinePrsm);
//    WriteArgumentResult(outFdrConbine, bContainer);
//    outConbine.close();
//    outFdrConbine.close();

        }

    }

    string convertToString(double d) {
        ostringstream os;
        if (os << d)
            return os.str();
        return "invalid conversion";
    }

    string cutProName(const string& str){
        std::string word;
        for (const char& c : str) {
            if (c == ' ') {
                break;
            }
            word += c;
        }
        return word;

    }

    string getProID(const string& str){
        std::string word;
        for (const char& c : str) {
            if (c >= '0' && c <= '9') {
                word += c;
            }
            else {
                continue;
            }
        }
        return word;

    }





}