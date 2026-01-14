//
// Created by wenzhong on 2023/3/21.
//

#ifndef WTOP_COMPLEMENT_IONS_HPP
#define WTOP_COMPLEMENT_IONS_HPP

#include "msalign.hpp"

using namespace std;
class complement_ions
{
public:

    complement_ions(){

    };
    
    complement_ions(const msalign_ptr &msalignPtr) {
        getComplementIonsMap(msalignPtr);
    }

    void getComplementIonsMap(const msalign_ptr &msalignPtr);

    map<int, int> getIndexMap()
    {
        return indexMap;
    };

    set<int> getMember()
    {
        return indexSet;
    }

private:
    map<int, int> indexMap;
    set<int> indexSet;
    double massError = 2.0 ; 
};
typedef complement_ions *complement_ions_ptr;

#endif // WTOP_COMPLEMENT_IONS_HPP
