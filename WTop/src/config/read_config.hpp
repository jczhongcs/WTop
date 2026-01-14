//
// Created by wenzhong on 2023/3/17.
//

#ifndef WTOP_READ_CONFIG_HPP
#define WTOP_READ_CONFIG_HPP

#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <memory>

#include "../console/opts.h"

using namespace std;

class read_config {
public:
    map<string,string> readConfigFile();

    const string &getProteinFile() const;

    const set<string> &getMsalignFileSet() const;

    const string &getPtmFile() const;

    int getThreadNum() const;

    bool get_input_config_arg(options opts);

    const string &getVarPtmCombinationOutFile() const;

    map<string,string> &get_arg_map() { return arg_map;}

    void insert_arg_map(string & key ,string & value) {
        arg_map.insert(make_pair(key,value));
    }
private:
    string protein_file ;
    set<string> msalign_file_set ;
    string ptm_file;
    std::string var_ptm_combination_out_file;
    int thread_num;

    map<string,string> arg_map ;

    bool judge_input_file();
};
typedef shared_ptr<read_config> read_config_sharePtr;
typedef read_config* read_config_Ptr;


#endif //WTOP_READ_CONFIG_HPP
