#include <map>
#include "io.h"

bool mylib::io::Readdata(const char *name, std::vector<double> &signals) {
    std::ifstream in(name);
    if (!in.good()) {
        return false;
    }

    while (in.good()) {
        double item;
        in >> item;
        if (in.fail()) {
            break;
        }
        signals.push_back(item);
    }
    in.close();
    return true;
}

void mylib::io::quickSort(std::vector<double> &arr, int begin, int end) {
    //如果区间不只一个数
    if (begin < end) {
        double temp = arr[begin]; //将区间的第一个数作为基准数
        int i = begin; //从左到右进行查找时的“指针”，指示当前左位置
        int j = end; //从右到左进行查找时的“指针”，指示当前右位置
        //不重复遍历
        while (i < j) {
            //当右边的数大于基准数时，略过，继续向左查找
            //不满足条件时跳出循环，此时的j对应的元素是小于基准元素的
            while (i < j && arr[j] > temp)
                j--;
            //将右边小于等于基准元素的数填入右边相应位置
            arr[i] = arr[j];
            //当左边的数小于等于基准数时，略过，继续向右查找
            //(重复的基准元素集合到左区间)
            //不满足条件时跳出循环，此时的i对应的元素是大于等于基准元素的
            while (i < j && arr[i] <= temp)
                i++;
            //将左边大于基准元素的数填入左边相应位置
            arr[j] = arr[i];
        }
        //将基准元素填入相应位置
        arr[i] = temp;
        //此时的i即为基准元素的位置
        //对基准元素的左边子区间进行相似的快速排序
        quickSort(arr, begin, i - 1);
        //对基准元素的右边子区间进行相似的快速排序
        quickSort(arr, i + 1, end);
    }
        //如果区间只有一个数，则返回
    else
        return;
}


void mylib::io::quickSort_ll(std::vector<double> &arr, int begin, int end, std::vector<std::vector<double> > &kk) {
    //如果区间不只一个数
    if (begin < end) {
        double temp = arr[begin]; //将区间的第一个数作为基准数
        int i = begin; //从左到右进行查找时的“指针”，指示当前左位置
        int j = end; //从右到左进行查找时的“指针”，指示当前右位置

        int a = kk[1][begin];
        int a1 = kk[2][begin];
        int a2 = kk[3][begin];
        int a3 = kk[4][begin];
        int a4 = kk[5][begin];

        int order = kk[6][begin];
        int pre = kk[7][begin];
        while (i < j) {

            while (i < j && arr[j] > temp)
                j--;

            arr[i] = arr[j];
            kk[1][i] = kk[1][j];
            kk[2][i] = kk[2][j];
            kk[3][i] = kk[3][j];
            kk[4][i] = kk[4][j];
            kk[5][i] = kk[5][j];
            kk[6][i] = kk[6][j];
            kk[7][i] = kk[7][j];

            while (i < j && arr[i] <= temp)
                i++;

            arr[j] = arr[i];
            kk[1][j] = kk[1][i];
            kk[2][j] = kk[2][i];
            kk[3][j] = kk[3][i];
            kk[4][j] = kk[4][i];
            kk[5][j] = kk[5][i];
            kk[6][j] = kk[6][i];
            kk[7][j] = kk[7][i];
        }

        arr[i] = temp;
        kk[1][i] = a;
        kk[2][i] = a1;
        kk[3][i] = a2;
        kk[4][i] = a3;
        kk[5][i] = a4;
        kk[6][i] = order;
        kk[7][i] = pre;

        quickSort_ll(arr, begin, i - 1, kk);

        quickSort_ll(arr, i + 1, end, kk);
    }
        //如果区间只有一个数，则返回
    else
        return;
}

void mylib::io::get_protein_name_seq_map(char *protein_file,
                                         std::map<std::string, std::string> &pro_name_seq_map) {
    int i, a, seq_len, seq_number = 0;
//    printf("%s\n", proteinFileName);
    mylib::data_stream::rtrim(protein_file);
    ifstream in(protein_file, ios::in);
    if (!in) {
        printf("Can't open %s for reading.\n", protein_file);
        return;
    }
    string title;
    string seq;
    string buff;
    vector<string> t;
    vector<string> s;
    while (in.good()) {
//        fgets(T1, LINE, in);
        getline(in, buff);
        mylib::data_stream::rtrim_s(buff);
        if (buff[0] != '>') {
            mylib::data_stream::rtrim_s(buff);
            seq += buff;
        } else {
            seq_number++;
            if (seq_number > 1) {
//               cout<<seq<<endl;
                s.emplace_back(seq);
            }
//            printf("%s", T1);
//            cout<<buff<<endl;
            t.emplace_back(buff);
            seq.clear();
        }
    }
//    printf("%s\n", fa);
//    cout<<seq<<endl;
    s.emplace_back(seq);
    for (int xi = 0; xi < t.size(); ++xi) {
        pro_name_seq_map[t[xi]] = s[xi];
    }
    in.close();
}

void Stringsplit(string str, const char split,vector<string>& rst)
{
    istringstream iss(str);	// 输入流
    string token;			// 接收缓冲区
    while (getline(iss, token, split))	// 以split为分隔符
    {
        rst.push_back(token) ;
    }

}


/**
 * 递归制造输入修饰组合函数
 * @param n 当前第几个修饰
 * @param arr   修饰质量arr
 * @param maxSize   每个修饰出现最大个数
 * @param arrT  修饰出现的个数
 * @param maxMass   组合最大质量
 */

void loop(int n,vector<double> & arr, int maxSize,vector<int> & arrT,double maxMass,vector<vector<int>>& ans)
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
        loop(n+1,arr,maxSize,arrT,maxMass,ans);
    }
}

void doubleArrayRemoveRepeat(vector<double>& array,vector<string> & arrStr)
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
void mylib::io::init_var_ptm_combination(const string & variableFileName, string outPath, int maxSize, double maxMass)
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
            Stringsplit(buff,',',str) ;
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
        EX_TRACE("variable PTMs size need <= 42 ");
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
    loop(0,arrMass,maxSize,arrLog,maxMass,ans);

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
    doubleArrayRemoveRepeat(combinationMass,combinationStr);
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
