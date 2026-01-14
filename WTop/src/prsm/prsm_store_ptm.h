//
// Created by LYC on 2024/3/20.
//

#ifndef WTOP_PRSM_STORE_PTM_H
#define WTOP_PRSM_STORE_PTM_H

#endif //WTOP_PRSM_STORE_PTM_H
class StorePtm{
public:
    StorePtm(int unimod, char uniChar, const string &ptmName, int leftBpPos, int rightBpPos, const string &stringMass, bool isUnknowPtm) :
        unimod(unimod), uniChar(uniChar), ptmName(ptmName), left_bp_pos(leftBpPos), right_bp_pos(rightBpPos), stringMass(stringMass), isUnknowPtm(isUnknowPtm) {
    }
    int unimod;//toppic中可变修饰的编码，未知修饰则为-1
    char uniChar;//唯一标识
    string ptmName;//修饰全名
    int left_bp_pos;//所修饰的起始位置（与toppic规则相同）
    int right_bp_pos;//所修饰的结束位置
    string stringMass;//修饰质量
    bool isUnknowPtm;//用于区分var_ptm和unk_ptm

    int getUnimod() const {
        return unimod;
    }
    void setUnimod(int unimod) {
        StorePtm::unimod = unimod;
    }
    char getUniChar() const {
        return uniChar;
    }
    void setUniChar(char uniChar) {
        StorePtm::uniChar = uniChar;
    }
    const string &getPtmName() const {
        return ptmName;
    }
    void setPtmName(const string &ptmName) {
        StorePtm::ptmName = ptmName;
    }
    int getLeftBpPos() const {
        return left_bp_pos;
    }
    void setLeftBpPos(int leftBpPos) {
        left_bp_pos = leftBpPos;
    }
    int getRightBpPos() const {
        return right_bp_pos;
    }
    void setRightBpPos(int rightBpPos) {
        right_bp_pos = rightBpPos;
    }
    const string &getStringMass() const {
        return stringMass;
    }
    void setStringMass(const string &stringMass) {
        StorePtm::stringMass = stringMass;
    }
    bool isUnknowPtm1() const {
        return isUnknowPtm;
    }
    void setIsUnknowPtm(bool isUnknowPtm) {
        StorePtm::isUnknowPtm = isUnknowPtm;
    }
};

typedef std::shared_ptr<StorePtm> StorePtmPtr;
typedef std::vector<StorePtmPtr> StorePtmPtrVec;