#include "wtop_prsm_processor.hpp"

namespace prsm {

/**
     * 处理prsm流程
     * @param msalign_file  质谱文件名
     * @param protein_file   蛋白文件名
     * @param var_ptm_result_file    可变修饰输出文件名
     * @param unknown_ptm_result_file   未知可变修饰输出文件名
     */
void WTop_Process::wtop_process(const string &msalign_file,
                                const string &var_ptm_file,
                                const string &var_ptm_result_file,
                                const string &unknown_ptm_result_file,
                                const vector<protein_ptr> &protein_container) // 匹配程序入口
{
    // 初始化参数
    config_argument_ptr config_prsm_arg = new config_argument(); // 创建prsm参数

    // init ptm
    cout << "variable ptms init - start" << endl;
    modify_ptr modptr = std::make_shared<prsm::Modify>(var_ptm_file);      // 创建修饰表指针.0表示跑H3数据集修饰，1表示跑WHIM2数据集
    modify_ptr one_mod_ptr = std::make_shared<prsm::Modify>(var_ptm_file); // 0 ,1 无所谓
    processor_management::initializerPtmsPtr(modptr, one_mod_ptr, one_ptm_table);
    cout << "variable ptms init - end" << endl;


    cout << "mass spectrum init - start" << endl;
    vector<msalign_ptr> msalign_ptr_container;
    processor_management::init_msalign_container(msalign_file, msalign_ptr_container);
    remove(msalign_file.c_str());
    cout << "mass spectrum init - end" << endl;


    cout << "filter out candidate proteins - start" << endl;
    filter_ptr filterPtr = new filter();
    //filterPtr->tagFilterProcess(msalign_ptr_container, protein_container, modptr, config_prsm_arg);


///    test wtop result use this function
//    filterPtr->filtrationProcess(msalign_ptr_container, protein_container, modptr, config_prsm_arg);


    filterPtr->getToppicFilterResult(msalign_ptr_container, protein_container, modptr, config_prsm_arg);





    
    delete filterPtr;
    cout << "filter out candidate proteins - end" << endl;

    // 未知修饰识别
    config_prsm_arg->selectPrSMsArg = 2;
    if (config_prsm_arg->selectPrSMsArg == 1) {
        cout << "one unknown ptm search - start" << endl;

        prsm_unknown_ptms_ptr prsmUnknownPtmsPtr = new prsm_unknown_ptms();
        prsmUnknownPtmsPtr->search_ptms_(one_mod_ptr, msalign_ptr_container, config_prsm_arg,
                                         var_ptm_file, unknown_ptm_result_file);
        delete prsmUnknownPtmsPtr;
        cout << "one unknown ptm search - end" << endl;
    }

    // 已知修饰识别(不识别未知修饰)
    if (config_prsm_arg->selectPrSMsArg == 0) {
        cout << "var ptms search - start" << endl;
        prsm_varible_ptms_ptr prsmVariblePtmsPtr = new prsm_varible_ptms();
        prsmVariblePtmsPtr->search_ptms_(msalign_ptr_container, modptr,
                                         config_prsm_arg, var_ptm_result_file);
        delete prsmVariblePtmsPtr;
        cout << "var ptms search - end" << endl;
    }

    if (config_prsm_arg->selectPrSMsArg == 2) {
        cout << "var ptms search - start" << endl;
        prsm_varible_ptms_ptr prsmVariblePtmsPtr = new prsm_varible_ptms();
        prsmVariblePtmsPtr->search_ptms_(msalign_ptr_container, modptr,
                                         config_prsm_arg, var_ptm_result_file);
        delete prsmVariblePtmsPtr;
        cout << "var ptms search - end" << endl;
        cout << "one unknown ptm search - start" << endl;
        prsm_unknown_ptms_ptr prsmUnknownPtmsPtr = new prsm_unknown_ptms();
        prsmUnknownPtmsPtr->search_ptms_(one_mod_ptr, msalign_ptr_container, config_prsm_arg,
                                         var_ptm_file, unknown_ptm_result_file);
        delete prsmUnknownPtmsPtr;
        cout << "one unknown ptm search - end" << endl;

    }


    // 释放空间
    delete config_prsm_arg;
//    delete modptr;
//    delete one_mod_ptr;
    return;
}
// 搜库流程


// Process

} // namespace prsm
