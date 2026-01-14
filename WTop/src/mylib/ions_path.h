#ifndef LONGPATH_H__
#define LONGPATH_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <map>
#include <set>
#include "data_steam.hpp"
#include "speculate_lib.h"

using namespace std;

struct vw {
    int v;
    int w;
};
const int len = 9024;

namespace mylib {
    class ions_path {
    public :
        vector<vw> G[len];
        int dp[len];
        bool visit[len];
        int choice[len];//存路径

        int nodeNumber = 0;
        int edgeNumber = 0;
        int start;
        int weight = 1;

        //需要让每个节点的下标唯一
        long tp = 100;
        std::map<long, int> preMut; //-------------first保存原始节点下标，second保存出现重复的次数。防止出现环。每次重复出现，将该节点下标乘以10000,在添加入图中.
        std::map<long, long> backPre; //-------------first用来保存变化后的节点，second保存由原始节点下标

        /** 1.将对齐数据转化为图
         *  2.以C端为参考。
         *  3.
         */
        void initializeData(vector<node> &N_All_Ions, vector<node> &C_All_Ions, const vector<double> &PtmsMassSeq,
                            vector<node> &resultIons, int seqLenth);

        //----------------------------------------------------------------
        int DP(int i);

        void printPath(int i, map<long, long> &ans);

        void process(vector<node> &ModifyIons, vector<double> &ptmsMass);

    };

    typedef std::shared_ptr<ions_path> ions_path_share_Ptr;
    typedef ions_path* ions_path_Ptr;
}

#endif
