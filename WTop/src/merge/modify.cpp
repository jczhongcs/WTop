#include "modify.hpp"

#include <cmath>

namespace prsm {

    void Stringsplit(string str, const char split,vector<string>& rst)
    {
        istringstream iss(str);	// 输入流
        string token;			// 接收缓冲区
        while (getline(iss, token, split))	// 以split为分隔符
        {
            rst.push_back(token) ;
        }

    }

    void Modify::Init_modify() {
        ifstream in;
        in.open(ptmsMapFile,ios::in);
        if (!in) {
            cout << "打开修饰表文件失败!" << endl;
            exit(20);
        }
        string t;
        int temp = 0 ;
        while (getline(in, t)) {
            if (t.find("PTMs") != t.npos) {
                temp = 1 ;
                continue;
            }
            if (!temp) {
                char data[30];
                int i = 0;
                for (i = 0; t[i] != ' '; ++i) {
                    data[i] = t[i];
                }
                t = t.substr(i + 1, t.size() - 1);
                double a = strtod(data, NULL);
                mapi[a] = t;
                mapStrMass[t] = a;
            } else {
                vector<string> logPTMs ;
                Stringsplit(t,' ',logPTMs);
                charPTMsMap.insert(make_pair(logPTMs[0].at(0),logPTMs[1]));
                charDtoSPtmsMap.insert(make_pair(logPTMs[0].at(0),logPTMs[2]));
            }
        }

        for(auto it = mapi.begin(); it != mapi.end(); ++it) {
            modifyMassSeq.push_back(it->first);
        }

    }


    string Modify::analysis(const double &mass) {
        auto it = mapi.begin();
        double temp = 65535;
        string ans;
        while (it != mapi.end()) {
            double t = fabs(mass - it->first);
            if (temp >= t) {
                temp = t;
                ans = it->second;
            }
            if (it->first - mass >= 100)
                break;
            ++it;
        }
        return ans;
    }

}