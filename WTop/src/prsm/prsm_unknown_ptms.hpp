//
// Created by wenzhong on 2023/3/22.
//

#ifndef WTOP_PRSM_UNKNOWN_PTMS_HPP
#define WTOP_PRSM_UNKNOWN_PTMS_HPP

#include <thread>
#include <vector>

#include "prsm_terminal_truncation.hpp"
#include "prsm_cwdtw.hpp"
#include "prsm_scope_alignment.hpp"
#include "prsm_alignment_filter.hpp"
#include "prsm_ions_location.hpp"
#include "prsm_score.hpp"
#include "prsm_msg_write.hpp"
#include "prsm_alignment_filter_var.hpp"
#include "candi_truncation.h"

#include "../ms/msalign.hpp"
#include "../protein/protein.hpp"
#include "../merge/modify.hpp"
#include "../merge/config_argument.hpp"
#include "../merge/temp_prsm_argument.hpp"

#include "../mylib//etd_process.hpp"


using namespace std;

class prsm_unknown_ptms {
public:

    prsm_unknown_ptms() {
        prsmCwdtwPtr = new prsm_cwdtw();
        prsmScopeAlignmentPtr = new prsm_scope_alignment();
        prsmAlignmentFilterPtr = new prsm_alignment_filter();
        prsmIonsLocationPtr = new prsm_ions_location();
        prsmScorePtr = new prsm_score();
        prsmMsgWritePtr = new prsm_msg_write();
    }

    void terminal_trunc_mass(vector<double> &theoryMass, vector<double> &theoryMassC,
                                                protein_ptr pro, int cutSizeN, int cutSizeC);

    void search_ptms_(prsm::modify_ptr one_mod_ptr,
                      vector<msalign_ptr> & msalign_ptrs,
                      prsm::config_argument_ptr configArgumentPtr,
                      const string & variablePTMsFileOutPath,
                      const string & outFIleName);

    void two_unknown_ptm_Process(prsm::modify_ptr mp, protein_ptr proteinPtr, msalign_ptr msalignPtr,
                                 prsm::config_argument_ptr prsmArg, const string &variablePTMsFileOutPath);

    void prsm_unk_process(prsm::modify_ptr modptr, protein_ptr proteinPtr, msalign_ptr msalignPtr,
                      prsm::temp_prsm_argument_ptr tempPrsmArgumentPtr, vector<double> &ptms_mass,
                      prsm::config_argument_ptr prsmArg,double var_ptm_mass,double unk_mass);

    void prsm_process_var(prsm::modify_ptr modptr,
                          protein_ptr proteinPtr,
                          msalign_ptr msalignPtr,
                          prsm::temp_prsm_argument_ptr tempPrsmArgumentPtr,
                          vector<double> &ptms_mass,
                          prsm::config_argument_ptr prsmArg);

    void storagePrSMsMsg(protein_ptr proteinPtr, msalign_ptr msalignPtr, prsm::temp_prsm_argument_ptr &tempPrsmArgumentPtr);

    int getCurMsSize() ;

     prsm_cwdtw *getPrsmCwdtwPtr() ;

     prsm_scope_alignment *getPrsmScopeAlignmentPtr() ;

     prsm_alignment_filter *getPrsmAlignmentFilterPtr() ;

     prsm_ions_location *getPrsmIonsLocationPtr() ;

     prsm_score *getPrsmScorePtr() ;

     prsm_msg_write *getPrsmMsgWritePtr() ;

    void setCurMsSize(int curMsSize);

    void setPrsmCwdtwPtr( prsm_cwdtw *prsmCwdtwPtr);

    void setPrsmScopeAlignmentPtr( prsm_scope_alignment *prsmScopeAlignmentPtr);

    void setPrsmAlignmentFilterPtr( prsm_alignment_filter *prsmAlignmentFilterPtr);

    void setPrsmIonsLocationPtr( prsm_ions_location *prsmIonsLocationPtr);

    void setPrsmScorePtr( prsm_score *prsmScorePtr);

    void setPrsmMsgWritePtr( prsm_msg_write *prsmMsgWritePtr);

    ~prsm_unknown_ptms() {
        delete prsmCwdtwPtr ;
        delete prsmScopeAlignmentPtr;
        delete prsmAlignmentFilterPtr;
        delete prsmIonsLocationPtr;
        delete prsmScorePtr;
        delete prsmMsgWritePtr;
    }

private:
    int cur_ms_size = 1 ;
    prsm_cwdtw_ptr prsmCwdtwPtr ;
    prsm_scope_alignment_ptr prsmScopeAlignmentPtr ;
    prsm_alignment_filter_ptr prsmAlignmentFilterPtr;
    prsm_ions_location_ptr prsmIonsLocationPtr ;
    prsm_score_ptr prsmScorePtr;
    prsm_msg_write_ptr prsmMsgWritePtr;
    msalignPtrVec evalueMsVec;
};
typedef prsm_unknown_ptms* prsm_unknown_ptms_ptr;

#endif //WTOP_PRSM_UNKNOWN_PTMS_HPP
