#ifndef MODIFY_HPP
#define MODIFY_HPP

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <memory>
#include <vector>
#include <sstream>

using namespace std;

namespace prsm{
    class Modify{
        public:
        map<double,string> mapi;  //double - 组合质量, string - 组合修饰

        map<string,double> mapStrMass;  // string - 组合修饰 double - 组合质量
        string ptmsMapFile ;        //输入的修饰文件
        vector<double> modifyMassSeq ;  //double - 组合修饰

        map<char,string> charPTMsMap ;  // char - 唯一标识 string - 修饰全名
        map<char,string> charDtoSPtmsMap; //char - 唯一标识 string ptm质量
        void Init_modify();     //读取数据


        Modify(const string & ptmsMassMapFile) : ptmsMapFile(ptmsMassMapFile) {

        };

        void get_mod_Map(map<double,string> a ) {
            this->mapi = a ;
        };

        //获取double - string map
        map<double,string> getTable() { return mapi; }

        string analysis(const double & mass);   //由修饰质量转换为修饰

        void massStringToStringMass()
        {
            for(map<double,string>::iterator it = mapi.begin(); it!= mapi.end(); ++it) {
                mapStrMass.insert(make_pair(it->second,it->first)) ;
            }
        }


        void getModMassSeq(std::vector<double> &ptmMass) {
            for(auto it = mapi.begin(); it != mapi.end(); ++it) {
                ptmMass.push_back(it->first);
            }
        }

    };
    typedef shared_ptr<Modify> modify_ptr;
}
#endif