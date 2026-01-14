//
// Created by wenzhong on 2023/3/17.
//

#include "ptm_combination.hpp"
#include <climits>

void ptm_loop(int n,vector<double> & arr, int maxSize,vector<int> & arrT,double maxMass,vector<vector<int>>& ans)
{
    if(n == arr.size())  //如果n = arr.size 说明最后一层
        return;
    for(int i = 0; i < maxSize; i++)
    {
        arrT[n] = i ;
        if (n == arr.size() - 1 )
        {
            double sum = 0 ;
            for(int j = 0 ; j < arr.size() ; ++j)
            {
                sum += arrT[j] * arr[j] ;
            }
            if (sum <= maxMass) {
                ans.push_back(arrT);

//                for(int j = 0 ; j < arr.size() ; ++j)
//                {
//                    out << "w" << " : " << arrT[j] <<" ";
//                }
//                out << " sum = " << sum <<endl;
            }
        }
        ptm_loop(n+1,arr,maxSize,arrT,maxMass,ans);
    }
}

void ptm_doubleArrayRemoveRepeat(vector<double>& array,vector<string> & arrStr)
{
    int length = array.size();
    for (int i = 0; i < length; ++i)
    {
        if (array[i] == INT_MAX)
        {
            continue;
        }
        for (int j = i + 1; j < length; ++j)
        {
            if (fabs(array[i] - array[j]) <= 0.00001)
            {
                if (arrStr[i].length() > arrStr[j].length())
                {
                    arrStr[i] = arrStr[j];
                }
                array[j] = INT_MAX;
            }
        }
    }
}

/**
 * 初始化修饰
 * @param variableFileName
 * @param maxSize
 * @param maxMass
 */
void ptm_combination::init_var_ptm_combination(const string & variableFileName,
                                               string outPath, int maxSize, double maxMass)
{
    ifstream in(variableFileName, ios::in);
    ofstream out(outPath,ios::out) ;
    string buff = "";
    map<double,string> ptmMassMap ;
    vector<char> logPTMs ;      //标识修饰log
    map<char,string> logPTMsMap ;
    map<char,double> logMassMap ;
    int i = 0 ;
    while(in.good())
    {
        getline(in,buff) ;
        if (buff.find(",") != buff.npos)
        {
            vector<string> str;
            utils_string_util::Stringsplit(buff,',',str) ;

            ptmMassMap.insert(make_pair(atof(str[1].c_str()),str[0]));
        }
    }
    //利用字母标识输入最多42个输入修饰
    for (int i = 65 ; i < 65 + 26; ++i)
    {
        logPTMs.push_back(i) ;
    }
    for (int i = 97 ; i < 97 + 26; ++i)
    {
        logPTMs.push_back(i) ;
    }
    if (ptmMassMap.size() > 42)
    {
        cout << "variable PTMs size need <= 42 " << endl;
        exit(42);
    }
    i = 0 ;
    string logStr = "" ;
    vector<double> arrMass ;
    vector<vector<int> > ans ; //存放组合数
    for (auto it = ptmMassMap.begin(); it != ptmMassMap.end() && i < logPTMs.size() ; ++it,++i)
    {
        logPTMsMap.insert(make_pair(logPTMs[i],it->second));    //log - PTMs string
        logMassMap.insert(make_pair(logPTMs[i],it->first));     // log - Mass double
        logStr += logPTMs[i] ;
        arrMass.push_back(it->first);
    }
    vector<int>  arrLog(arrMass.size(),0);
    ptm_loop(0,arrMass,maxSize,arrLog,maxMass,ans);

    map<double,string>  combinationMassMap ;
    vector<string> combinationStr ;
    vector<double> combinationMass;
    for(vector<int> v :ans)
    {
        string con = "";
        double conMass = 0.0;
        for (i = 0 ; i < v.size(); ++i)
        {
            char curChar = logStr[i] ;
            int countNum = v[i] ;
            while(countNum > 0 )
            {
                con += curChar ;
                conMass += logMassMap[curChar];
                countNum -= 1 ;
            }
        }
        combinationMass.push_back(conMass);
        if (fabs(conMass - 0.0) < 0.000001)
        {
            combinationStr.push_back("0");
        } else {
            combinationStr.push_back(con);

        }
    }
    ptm_doubleArrayRemoveRepeat(combinationMass,combinationStr);
    for (i = 0 ; i < combinationMass.size() ; ++i)
    {
        if (combinationMass[i] == INT_MAX)
        {
            continue;
        }
        combinationMassMap.insert(make_pair(combinationMass[i],combinationStr[i]));
    }
    for(auto it = combinationMassMap.begin(); it != combinationMassMap.end();++it)
    {
        out <<setiosflags(ios::fixed) << setprecision(4)<< it->first <<" "<< it->second<<endl;
    }
    out << "PTMs Map" <<endl;
    for (auto it = logPTMsMap.begin(); it != logPTMsMap.end(); ++it)
    {
        out << it->first << " " << it->second << " " << logMassMap[it->first]<<endl;
    }
    out.close();
    in.close() ;

}
