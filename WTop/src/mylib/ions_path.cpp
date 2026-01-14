#include "ions_path.h"

namespace mylib {

    int ions_path::DP(int i) {
        if (visit[i])
            return dp[i];
        visit[i] = true;
        for (int j = 0; j < G[i].size(); j++) {
            if (dp[i] < DP(G[i][j].v) + G[i][j].w) {//逆拓扑求最长路径
                dp[i] = DP(G[i][j].v) + G[i][j].w;
                choice[i] = G[i][j].v;//存后继节点
            }
        }
        return dp[i];
    }

    void ions_path::printPath(int i, map<long, long> &ans) {
        // cout<<i;
        ans[i] = 1;
        while (choice[i] != -1) {
            // cout<<"->"<<choice[i];
            ans[choice[i]] = 1;
            i = choice[i];
        }
    }

    void ions_path::process(vector<node> &ModifyIons, vector<double> &ptmsMass) {
        int n, m, start;
        // cin>>n>>m>>start;
        n = nodeNumber;
        m = edgeNumber;
        fill(dp, dp + len, 0);
        fill(visit, visit + len, false);
        fill(choice, choice + len, -1);
        for (int i = 0; i < len; i++) {
            choice[i] = -1;
        }
        for (int i = 0; i < m; i++) {
            int x, y, w;
            cin >> x >> y >> w;       //------x为起始节点，y为目标节点，w为权值
            vw t;
            t.v = y;
            t.w = w;
            G[x].push_back(t);
        }
        // cout<<DP(start)<<endl;
        // printPath(start);
        return;
    }


    void ions_path::initializeData(vector<node> &N_All_Ions, vector<node> &C_All_Ions,
                                   const vector<double> &PtmsMassSeq, vector<node> &resultIons,
                                   int seqLenth) {
        int n, m, start = -1;
        n = 0;
        m = edgeNumber;
        int log = seqLenth - 2;  // 利用log - 2 - 下标  = 互补峰下标
        fill(dp, dp + len, 0);
        fill(visit, visit + len, false);
        fill(choice, choice + len, -1);
        if (N_All_Ions.empty()) {
            return;
        }
        int max = N_All_Ions.back().thoe_id + 1;
        //数据转变，每一个节点唯一话
        vector<node> ionsT = N_All_Ions;
        map<long, long> preMut;    //前一个保存出现的节点，后一个保存出现的倍数(-1,-2,-3)
        map<long, long> backPre;    //first
        map<long, set<long> > mapSet;
//         for(int i = 0 ; i < ionsT.size() ; i++ ) {
//         	cout<<"thero = "<<ionsT[i].thoe_id<<" mono id = "<<ionsT[i].mono_id<<" mass = "<<ionsT[i].index<<endl;
//         }
        /**
         * 构造图
         */
        for (int xi = 0; xi < ionsT.size(); xi++) {
            n++;
            if (!preMut.count(ionsT[xi].thoe_id)) {
                preMut[ionsT[xi].thoe_id] = max;
            } else {
                long t = ionsT[xi].thoe_id;
                ionsT[xi].thoe_id = ionsT[xi].thoe_id + preMut[ionsT[xi].thoe_id];        //如果有，将该位点扩大
                if (preMut[t] + max <= len) {            //如果该位点可行
                    preMut[t] = preMut[t] + max;        //将该倍数扩大
                } else {
                    continue;
                }
                backPre[ionsT[xi].thoe_id] = t;
                mapSet[t].insert(ionsT[xi].thoe_id);
                mapSet[t].insert(t);
            }
        }
//         for(int i = 0 ; i < ionsT.size() ; i++ ) {
//         	cout<<"thero = "<<ionsT[i].thoe_id<<" mono id = "<<ionsT[i].mono_id<<" mass = "<<ionsT[i].index<<endl;
//         }
//         for(auto it = mapSet.begin() ; it != mapSet.end() ; it++ ) {
//         	cout<<"it->first = "<<it->first <<endl;
//         	auto it2 = it->second.begin();
//         	while(it2 != it->second.end() ) {
//         		cout<<" " <<(*it2)<<" ";
//         		it2 ++ ;
//         	}
//         	cout<<endl ;
//         }
        //------------------------
        for (int xi = 0; xi < ionsT.size(); xi++) {
            for (int ji = xi + 1; ji < ionsT.size(); ji++) {
                int n = mylib::data_stream::txpb_binarey_search_ex(PtmsMassSeq, PtmsMassSeq.size(),
                                                        ionsT[ji].index - ionsT[xi].index);
                if (n < 0)  continue;
                if (fabs( (ionsT[ji].index - ionsT[xi].index) - PtmsMassSeq[n]) < 3) {    //满足相等或存在修饰
                    if (!backPre.count(ionsT[ji].thoe_id)) {    // 该节点不是任何点变化而来
                        vw t;
                        t.v = ionsT[ji].thoe_id;
                        t.w = 1;
                        G[ionsT[xi].thoe_id].push_back(t);
                    } else {                //该节点由其他点变化而来
                        // mapSet[backPre[ionsT[ji].thoe_id]] 代表原始点所变化点产生的集合
                        if (!mapSet[backPre[ionsT[ji].thoe_id]].count(ionsT[xi].thoe_id)) { //代表如果前驱节点不为集合中的点,直接添加节点
                            vw t;
                            t.v = ionsT[ji].thoe_id;
                            t.w = 1;
                            G[ionsT[xi].thoe_id].push_back(t);
                        } else { // 如果前驱节点是集合中的点,相差不大才有边
                            if (fabs(ionsT[ji].index - ionsT[xi].index) < 3) {
                                vw t;
                                t.v = ionsT[ji].thoe_id;
                                t.w = 1;
                                G[ionsT[xi].thoe_id].push_back(t);
                            }
                        }
                    }
                }
            }   // for end
        }   // for end
        map<long, int> startMap;
        int maxLong = 0;
        int bestStart = 0;
        for (int xi = 0; xi < ionsT.size(); xi++) {
            int pLong = 0;
            fill(dp, dp + len, 0);
            fill(visit, visit + len, false);
            fill(choice, choice + len, -1);
            if (!startMap.count(start)) {
                start = ionsT[xi].thoe_id;
                pLong = DP(start);
                if (pLong > maxLong) {
                    maxLong = pLong;
                    bestStart = start;
                }
            }
        }
        // fill(dp,dp+len,0);
        // fill(visit,visit+len,false);
        // fill(choice,choice+len,-1);
        DP(bestStart);
        map<long, long> ans;
        printPath(bestStart, ans);
        //------------------------解析
        for (int xi = 0; xi < ionsT.size(); xi++) {
            if (ans.count(ionsT[xi].thoe_id)) {
                if (backPre.count(ionsT[xi].thoe_id)) {
                    ionsT[xi].thoe_id = backPre[ionsT[xi].thoe_id];
                }
                resultIons.push_back(ionsT[xi]);
            }
        }
        // for (int xi = 0 ; xi < resultIons.size() ; xi++ ) {
        // 	cout<<" theroId = "<<resultIons[xi].thoe_id<<" ,monoId = "<<resultIons[xi].mono_id<<" ,mass = "<<resultIons[xi].index<<endl;
        // }
        return;
    }

}
