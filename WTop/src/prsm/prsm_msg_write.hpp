//
// Created by wenzhong on 2023/3/22.
//

#ifndef WTOP_PRSM_MSG_WRITE_HPP
#define WTOP_PRSM_MSG_WRITE_HPP

#include <tinyxml.h>
#include "../ms/msalign.hpp"
#include "../protein/protein.hpp"
#include "../merge/modify.hpp"

class prsm_msg_write {
public:


    static void writeResultNews(msalign_ptr msalignPtr,
                                         string Result_file,
                                         prsm::modify_ptr modptr);


    static void OutXmlToEvalue( msalignPtrVec msalignptrVec,
                                         string xmlFileName);

    static void AddPrsm(TiXmlElement* prsm_list,msalign_ptr msalignPtr,prsm::modify_ptr modptr);

    static std::string doubleToStringWithPrecision(double value, int precision);

    static std::string generateSubstring(const std::string& str);
};
typedef prsm_msg_write* prsm_msg_write_ptr;

#endif //WTOP_PRSM_MSG_WRITE_HPP
