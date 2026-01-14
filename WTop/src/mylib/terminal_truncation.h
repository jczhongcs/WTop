#ifndef CUTLOCATION_H__
#define CUTLOCATION_H__

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <memory>

using namespace std;

namespace c {
    namespace arg {

    class cutLocation {
        public:
        int cutLocationC = 0 ;
        int cutLocationN = 0 ;

        double theroMass = 0 ;
        double precursorMass = 0 ;
        double variablePtmsMasss = 0 ;
        double ppm = 0 ; 
        string ptmsSeq;
    };
    typedef shared_ptr<cutLocation> cutLocationPtr;
    typedef cutLocation* cutLocationPtrP;
    }
}
#endif 