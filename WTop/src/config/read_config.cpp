//
// Created by wenzhong on 2023/3/17.
//

#include "read_config.hpp"

map<string, string> read_config::readConfigFile() {
    map<string, string> config_arg_map;
    string config_file = "config_argument";
    ifstream in(config_file, ios::in);
    string buff;
    int msalign_file_num = 0;

    while (getline(in, buff)) {
        if (buff.size() == 0 || buff[0] == '#') {
            continue;
        }
        string argument_name = buff.substr(0, buff.find("="));
        string argument_val = buff.substr(buff.find("=") + 1);
        if (!argument_val.empty() && argument_val.back() == 13) {
            argument_val.pop_back();
        }
        if (argument_name.find("protein_file") != std::string::npos) {
            protein_file = argument_val;
        } else if (argument_name.find("msalign_file") != std::string::npos) {
            argument_name += +"_" + to_string(++msalign_file_num);
            msalign_file_set.insert(argument_val);
        } else if (argument_name.find("ptm_file") != std::string::npos) {
            ptm_file = argument_val;
        } else if (argument_name.find("thread_num") != std::string::npos) {
            thread_num = stoi(argument_val);
        }
        config_arg_map.insert(make_pair(argument_name, argument_val));
    }

    return config_arg_map;
}

bool read_config::judge_input_file() {
    if (protein_file.empty() && protein_file.find(".fasta") == std::string::npos) {
        cout << "no Protein File Path " << endl;
        return false;
    }
    if (msalign_file_set.empty()) {
        cout << "no msalign File Path " << endl;
        return false;
    }
    return true;
}
bool read_config::get_input_config_arg(struct options opts) {
    string style_ = opts.style;
    //command input
    if (style_.find("command") != std::string::npos) {
        protein_file = opts.input;          //蛋白质文件
        msalign_file_set.insert(opts.peer); //质谱文件
        ptm_file = opts.output;             //输出文件
        thread_num = opts.test;

        arg_map.insert(make_pair("protein_file:", protein_file));
        arg_map.insert(make_pair("msalign_file:", opts.peer));
        arg_map.insert(make_pair("ptm_file:", ptm_file));
        arg_map.insert(make_pair("thread_num:", to_string(thread_num)));
    } else if (style_.find("config") != std::string::npos) { // config input
        arg_map = read_config::readConfigFile();
    }

    //judge input file is right
    if (ptm_file.size() > 0 && ptm_file.find(".txt") != std::string::npos) {
        var_ptm_combination_out_file = ptm_file.substr(0, ptm_file.find_last_of(".")) + "_out_ptm_map.txt";
    }
    return judge_input_file();
}

const string &read_config::getProteinFile() const {
    return protein_file;
}

const set<string> &read_config::getMsalignFileSet() const {
    return msalign_file_set;
}

const string &read_config::getPtmFile() const {
    return ptm_file;
}

int read_config::getThreadNum() const {
    return thread_num;
}

const string &read_config::getVarPtmCombinationOutFile() const {
    return var_ptm_combination_out_file;
}
