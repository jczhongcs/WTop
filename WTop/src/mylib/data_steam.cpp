#include "data_steam.hpp"

// #include "../process/opts.h"
// #include "../mylib/io.h"
// #include "../mylib/fun.h"
// #include "../mylib/cwdtw.h"

using namespace std;

#define eps 1e-10

void mylib::data_stream::rtrim(char *str)
{
    char *s;
    s = str + strlen(str);
    while (--s >= str)
    {
        if (!isspace(*s))
            break;
        *s = 0;
    }
}

void mylib::data_stream::rtrim_s(string &str)
{
    char *s = (char *)malloc((str.length() + 1) * sizeof(char));
    strcpy(s, str.c_str());
    rtrim(s);
    str = s;
    free(s);
}

void WriteSequenceAlignment(const char *output,
                            const std::vector<double> &reference_orig, const std::vector<double> &peer_orig,
                            const std::vector<double> &reference, const std::vector<double> &peer,
                            vector<pair<long, long>> &alignment, int swap, double tdiff)
{
    vector<std::string> tmp_rec;
    double diff;

    for (long i = 0; i < alignment.size(); i++)
    {
        //----- output to string ----//
        std::ostringstream o;
        diff = std::fabs(reference[alignment[i].first] - peer[alignment[i].second]);
        if (swap == 1)
        {
            o << setw(10) << alignment[i].second + 1 << " " << setw(10) << alignment[i].first + 1 << " | ";

            o << setw(10) << fabs(reference_orig[alignment[i].second] - peer_orig[alignment[i].first]) << " | ";

            o << setw(15) << reference_orig[alignment[i].second] << " " << setw(15) << peer_orig[alignment[i].first]
              << " | ";
            o << setw(15) << peer[alignment[i].second] << " " << setw(15) << reference[alignment[i].first];
        }
        else
        {
            o << setw(10) << alignment[i].first << " " << setw(10) << alignment[i].second << " | ";

            o << setw(10) << fabs(reference_orig[alignment[i].first] - peer_orig[alignment[i].second]) << " | ";

            o << setw(15) << reference_orig[alignment[i].first] << " " << setw(15) << peer_orig[alignment[i].second]
              << " | ";
            o << setw(15) << reference[alignment[i].first] << " " << setw(15) << peer[alignment[i].second];
        }
        o << "          diff:" << setw(15) << diff;
        //----- record string -----//
        std::string s = o.str();
        tmp_rec.push_back(s);
    }
    //----- output to file ------//

    FILE *fp = fopen(output, "wb");
    for (long i = 0; i < (long)tmp_rec.size(); i++)
        fprintf(fp, "%s\n", tmp_rec[i].c_str());
    // fprintf(fp,"| tdiff = %lf | alignment.size() = %d | tdiff/alignment.size() = %lf |\n",tdiff ,alignment.size(),tdiff/alignment.size());
    fclose(fp);
}

//
double mylib::data_stream::cwdtw_map(vector<double> &reference_1, vector<double> &peer_1,
                                     std::vector<std::pair<long, long>> &alignment,
                                     vector<double> &ref_zscore, vector<double> &peer_zscore)
{
    struct options
    {
        //-> required parameter
        char output[65532];
        //-> key parameter
        int radius;
        int level;
        float scale0;
        //-> vice parameter
        int verbose;
        int test;
        int mode;
    };
    options opts;
    opts.radius = 50;
    opts.level = 1;
    opts.scale0 = sqrt(2);
    opts.verbose = 0; //-> [0] no verbose; 1 verbose
    opts.test = 0;    //-> [0] not use test mode; 1 equal_ave, 2 peak_ave, 3 Fast_DTW
    opts.mode = 0;    //-> [0] block bound; 1 diagonol bound

    //======================= START Procedure ===================================//
    vector<double> reference = reference_1;
    std::vector<double> reference_orig = reference;

    vector<double> peer;

    for (int i = 0; i < peer_1.size(); i++)
    {
        peer.push_back(peer_1[i]);
    }
    std::vector<double> peer_orig = peer;
    //----- length check ------//
    int swap = 0;
    if (reference.size() > peer.size())
    {
        std::vector<double> tmp = peer;
        peer = reference;
        reference = tmp;
        tmp = ref_zscore;
        ref_zscore = peer_zscore;
        peer_zscore = tmp;
        swap = 1;
    }

    //==================================================//
    //------3. process initial input signals ----------//
    //----- 3.1 Zscore normaliza on both signals -----//
    double avg, dev;

    //-------3.1 test consol Zsocrenormalize --------------//
    std::vector<double> consol;
    int i = 0;

    for (i = 0; i < peer.size(); i++)
    {
        consol.push_back(peer[i]);
    }
    for (i = 0; i < reference.size(); i++)
    {
        consol.push_back(reference[i]);
    }
    mylib::fun::ZScoreNormalize(consol, &avg, &dev);
    for (i = 0; i < peer.size(); i++)
    {
        peer[i] = consol[i];
    }
    for (int j = 0; i < consol.size(); i++, j++)
    {
        reference[j] = consol[i];
    }
    ref_zscore = reference;
    peer_zscore = peer;
    //----- 3.2 calculate length ratio between input signals -----//
    double alpha = (double)peer.size() / reference.size();

    //====================================================//
    //----- 4. continous wavelet transform --------------//
    std::vector<std::vector<double>> rcwt, pcwt;

    if (opts.verbose == 1)
    {
        EX_TRACE("CWT Analysis...\n");
    }

    long npyr = opts.level;      // default: 3
    double scale0 = opts.scale0; // default: sqrt(2)
    double dscale = 1;           // default: 1

    mylib::cwdtw::CWTAnalysis(reference, rcwt, scale0, dscale, npyr);
    mylib::cwdtw::CWTAnalysis(peer, pcwt, scale0 * alpha, dscale, npyr);

    //------ 4.1 Zscore normaliza on both CWT signals -----//
    // if multiscale is used, pyr logical should be added.
    for (long i = 0; i < npyr; i++)
    {
        mylib::fun::ZScoreNormalize(rcwt[i], &avg, &dev);
        mylib::fun::ZScoreNormalize(pcwt[i], &avg, &dev);
    }

    //============================================//
    //------ 5. multi-level WaveletDTW ----------//
    std::vector<std::pair<long, long>> cosali;
    double tdiff;

    if (opts.verbose == 1)
    {
        EX_TRACE("Coarse Alignment...\n");
    }
    mylib::cwdtw::MultiLevel_WaveletDTW(reference, peer, rcwt, pcwt, cosali, opts.radius, opts.test, opts.mode, &tdiff);
    if (opts.verbose == 1)
    {
        EX_TRACE("Average Deviation (%.1lf/%ld=%.3lf)\n", tdiff, cosali.size(), tdiff / cosali.size());
    }

    //------ 5.1 generate final boundary -------//
    std::vector<std::pair<long, long>> bound;
    mylib::cwdtw::BoundGeneration(cosali, opts.radius, bound, opts.mode);

    //------ 5.2 generate final alignment via cDTW ------//
    // std::vector<std::pair<long,long> > alignment;
    tdiff = mylib::fun::BoundDynamicTimeWarping(reference, peer, bound, alignment);
    // fprintf(stderr,"%lf %d %lf\n",tdiff,alignment.size(),tdiff/alignment.size());

    double c1 = tdiff / alignment.size();

    //=================================================//
    //------ 6. output final alignment to file -------//
    // if(output!="")
    // {
    // 	WriteSequenceAlignment(output, reference_orig, peer_orig, reference, peer, alignment, swap,tdiff);
    // }
    //----- exit -----//
    //=================================================//
    if (swap == 1)
    {
        for (i = 0; i < alignment.size(); i++)
        {
            long t;
            t = alignment[i].first;
            alignment[i].first = alignment[i].second;
            alignment[i].second = t;
        }
        std::vector<double> tmp = peer;
        peer = reference;
        reference = tmp;
        tmp = ref_zscore;
        ref_zscore = peer_zscore;
        peer_zscore = tmp;
        swap = 1;
    }
    return tdiff;
}

/**
 * 得到所有的sub值，从sub中选取合适的峰值
 * @param reference
 * @param reference_y
 * @param peer
 * @param sub
 * @param sub_y
 * @param alignment
 * @param alignment_y
 * @return
 */
int mylib::data_stream::Get_sub(const vector<double> &reference, const vector<double> &reference_y,
                                const vector<double> &peer,
                                vector<double> &sub, vector<double> &sub_y,
                                const std::vector<std::pair<long, long>> &alignment,
                                const std::vector<std::pair<long, long>> &alignment_y)
{
    int i = 0;
    for (i = 0; i < alignment.size(); i++)
    {
        sub.push_back((peer[alignment[i].second] - reference[alignment[i].first]));
    }
    for (i = 0; i < alignment_y.size(); i++)
    {
        sub_y.push_back(peer[alignment_y[i].second] - reference_y[alignment_y[i].first]);
    }
    return 0;
}

void mylib::data_stream::WriteSequenceAlignment_peer_lons(const char *output,
                                                          const std::vector<double> &peer_orig,
                                                          std::vector<locat_ions> location)
{
    vector<std::string> tmp_rec;
    for (long i = 0; i < peer_orig.size(); i++)
    {
        //----- output to string ----//
        std::ostringstream o;
        o << setw(10) << setiosflags(ios::fixed) << setprecision(4) << peer_orig[i] << " | ";
        if (location[i].c != "x" || location[i].n != "x")
            o << setw(10) << location[i].c << " | " << setw(10) << location[i].n << " | ";
        //----- record string -----//
        std::string s = o.str();
        tmp_rec.push_back(s);
    }
    //----- output to file ------//
    FILE *fp = fopen(output, "wb");
    for (long i = 0; i < (long)tmp_rec.size(); i++)
        fprintf(fp, "%s\n", tmp_rec[i].c_str());
    fclose(fp);
}

void mylib::data_stream::Scope_Align_ReWrite(const vector<double> &reference, const vector<double> &peer,
                                             std::vector<std::pair<long, long>> &alignment,
                                             std::vector<double> &modifier, const double &modifier_mass)
{
    int i = 0, j;
    vector<pair<long, long>>::iterator iter = alignment.begin();
    vector<pair<long, long>>::iterator iter_temp = alignment.begin();
    int n = alignment.size();
    vector<pair<long, long>> s;
    vector<pair<long, long>> m;
    int dup = 1;
    double sValue;
    int index;
    // 对于第一个数据的范围搜索
    for (i = 0; i < 10 && i < reference.size(); i++)
    {
        if (mylib::data_stream::find_key(peer[0] - reference[i], modifier) != -1)
        {
            if (mylib::data_stream::remove_false(peer[0] - reference[i], modifier, modifier_mass))
            {
                m.push_back(make_pair(i, 0));
            }
        }
    }
    // 对于最后一个数据进行的范围搜索
    for (i = reference.size() - 1; i >= reference.size() - 10 && i >= 0; --i)
    {
        if (mylib::data_stream::find_key(peer.back() - reference[i], modifier) != -1)
        {
            if (mylib::data_stream::remove_false(peer.back() - reference[i], modifier, modifier_mass))
            {
                m.push_back(make_pair(i, peer.size() - 1));
            }
        }
    }
    for (i = 1; i < n - 1; i++)
    {
        double min = peer[alignment[i].second] - modifier_mass - 10;
        double max = peer[alignment[i].second] + (modifier_mass + 10);
        if (peer[alignment[i].second] != peer[alignment[i + 1].second])
        {
            for (j = alignment[i - dup].first; j > 0; j--)
            {
                if (reference[j] > min)
                {
                    double key = peer[alignment[i].second] - reference[j];
                    if (fabs(key) < 1 && (reference[j] != reference[alignment[i].first]))
                    {
                        s.push_back(make_pair(j, alignment[i].second));
                    }
                    if (mylib::data_stream::find_key(key, modifier) != -1 && (reference[j] != reference[alignment[i].first]))
                    {
                        if (mylib::data_stream::remove_false(key, modifier, modifier_mass))
                        {
                            m.push_back(make_pair(j, alignment[i].second));
                        }
                    }
                }
                else
                {
                    break;
                }
            }
            for (j = alignment[i + 1].first; j < reference.size(); j++)
            {
                if (reference[j] < max)
                {
                    double key = peer[alignment[i].second] - reference[j];
                    if (fabs(key) < 1 && (reference[j] != reference[alignment[i].first]))
                    {
                        s.push_back(make_pair(j, alignment[i].second));
                    }
                    if (mylib::data_stream::find_key(key, modifier) != -1 && (reference[j] != reference[alignment[i].first]))
                    {
                        if (mylib::data_stream::remove_false(key, modifier, modifier_mass))
                        {
                            m.push_back(make_pair(j, alignment[i].second));
                        }
                    }
                }
                else
                {
                    break;
                }
            }
            dup = 1;
        }
        else
        {
            dup++;
        }
    }
    for (i = 0; i < s.size(); i++)
    {
        alignment.push_back(s[i]);
    }
    for (i = 0; i < m.size(); i++)
    {
        alignment.push_back(m[i]);
    }
    std::sort(alignment.begin(), alignment.end(), mylib::data_stream::cmp1);
    alignment.erase(unique(alignment.begin(), alignment.end()), alignment.end());
}

/**
 * 2022-9-19
 * @param reference
 * @param peer
 * @param alignment
 * @param modifier
 * @param modifier_mass
 *          //满足下述条件可以认为是需要添加的数据对齐对
            //1、当两者相差不大
            //2、如果相差较大，是否满足ptmMassShift
            //2.1 当是 modifier_mass 的值，直接加入
            //2.2 当是 modifier_mass 的子值，判断是否有可达modifier_mass的质量，如果有，则加入，没有，则不加入
            //首先去重复。
 */
void mylib::data_stream::scopeAlignReWrite(const vector<double> &reference, const vector<double> &peer,
                                           std::vector<std::pair<long, long>> &alignment,
                                           std::vector<double> &modifier, const double &modifier_mass, double maxMass)
{
    vector<pair<long, long>> addAlignment;
    int sex = 5; // 5da 常熟项
    set<long> removeRef;
    for (int i = 0; i < alignment.size(); ++i)
    {
        if (removeRef.count(alignment[i].first))
        { // 如果一个理论值已经对准过，直接进行下一个
            continue;
        }
        else
        {
            removeRef.insert(alignment[i].first);
        }
        set<long> removePeer;
        removePeer.insert(alignment[i].second);
        double theoMass = reference[alignment[i].first];
        // 向上边查找
        for (long left = alignment[i].second - 1; left >= 0; --left)
        {
            double monoMass = peer[left];
            if (fabs(monoMass - theoMass) > maxMass)
            {
                break;
            }
            if (removePeer.count(left))
            {
                continue;
            }
            double monoSubTheo = monoMass - theoMass;
            if (fabs(monoSubTheo) < sex)
            { // 1、当两者相差不大
                addAlignment.push_back(make_pair(alignment[i].first, left));
            }
            else if (fabs(monoSubTheo - modifier_mass) < sex)
            { // 2.1 当是 modifier_mass 的值，直接加入
                addAlignment.push_back(make_pair(alignment[i].first, left));
            }
            else
            { // 2.2 当是modifier_mass的子值，判断是否有可达modifier_mass的质量，如果有，则加入，没有，则不加入
                double subMass = modifier_mass - monoSubTheo;
                double searchIndex = mylib::data_stream::txpb_binarey_search_ex(modifier, modifier.size(), subMass);
                if (searchIndex < 0)
                    continue;
                double temp = modifier[searchIndex] - subMass;
                if (fabs(temp) < sex)
                {
                    addAlignment.push_back(make_pair(alignment[i].first, left));
                }
            }
            removePeer.insert(left);
        }
        // 向下边查找
        for (long right = alignment[i].second + 1; right < peer.size(); ++right)
        { // 向下边查找 指向mono的指针
            double monoMass = peer[right];
            if (fabs(monoMass - theoMass) > maxMass)
            {
                break;
            }
            if (removePeer.count(right))
            {
                continue;
            }
            double monoSubTheo = monoMass - theoMass;
            if (fabs(monoSubTheo) < sex)
            { // 1、当两者相差不大
                addAlignment.push_back(make_pair(alignment[i].first, right));
            }
            else if (fabs(monoSubTheo - modifier_mass) < sex)
            { // 2.1 当是 modifier_mass 的值，直接加入
                addAlignment.push_back(make_pair(alignment[i].first, right));
            }
            else
            { // 2.2 当是modifier_mass的子值，判断是否有可达modifier_mass的质量，如果有，则加入，没有，则不加入
                double subMass = modifier_mass - monoSubTheo;
                double searchIndex = mylib::data_stream::txpb_binarey_search_ex(modifier, modifier.size(), subMass);
                if (searchIndex < 0)
                    continue;
                double temp = modifier[searchIndex] - subMass;
                if (fabs(temp) < sex)
                {
                    addAlignment.push_back(make_pair(alignment[i].first, right));
                }
            }
            removePeer.insert(right);
        }
    }
    for (int i = 0; i < addAlignment.size(); ++i)
    {
        alignment.push_back(addAlignment[i]);
    }
    std::sort(alignment.begin(), alignment.end(), mylib::data_stream::cmp1);
    alignment.erase(unique(alignment.begin(), alignment.end()), alignment.end());
}

void mylib::data_stream::Scope_Align(const vector<double> &reference, const vector<double> &peer,
                                     std::vector<std::pair<long, long>> &alignment,
                                     std::vector<double> &modifier, const double &modifier_mass)
{
    int i = 0, j;
    vector<pair<long, long>>::iterator iter = alignment.begin();
    vector<pair<long, long>>::iterator iter_temp = alignment.begin();
    int n = alignment.size();
    vector<pair<long, long>> s;
    vector<pair<long, long>> m;
    int dup = 1;
    // 对于第一个数据的范围搜索
    for (i = 0; i < 10 && i < reference.size(); i++)
    {
        if (mylib::data_stream::find_key(peer[0] - reference[i], modifier) != -1)
        {
            if (mylib::data_stream::remove_false(peer[0] - reference[i], modifier, modifier_mass))
            {
                m.push_back(make_pair(i, 0));
            }
        }
    }
    // 对于最后一个数据进行的范围搜索
    for (i = reference.size() - 1; i >= reference.size() - 10 && i >= 0; --i)
    {
        if (mylib::data_stream::find_key(peer.back() - reference[i], modifier) != -1)
        {
            if (mylib::data_stream::remove_false(peer.back() - reference[i], modifier, modifier_mass))
            {
                m.push_back(make_pair(i, peer.size() - 1));
            }
        }
    }

    for (i = 1; i < n - 1; i++)
    {
        double min = peer[alignment[i].second] - modifier_mass - 10;
        double max = peer[alignment[i].second] + (modifier_mass + 10);
        if (peer[alignment[i].second] != peer[alignment[i + 1].second])
        {
            for (j = alignment[i - dup].first; j > 0; j--)
            {
                if (reference[j] > min)
                {
                    double key = peer[alignment[i].second] - reference[j];
                    if (fabs(key) < 1 && (reference[j] != reference[alignment[i].first]))
                    {
                        s.push_back(make_pair(j, alignment[i].second));
                    }
                    /**
                     * find_key
                     */
                    if (mylib::data_stream::find_key(key, modifier) != -1 && (reference[j] != reference[alignment[i].first]))
                    {
                        if (mylib::data_stream::remove_false(key, modifier, modifier_mass))
                        {
                            m.push_back(make_pair(j, alignment[i].second));
                        }
                    }
                }
                else
                {
                    break;
                }
            }
            for (j = alignment[i + 1].first; j < reference.size(); j++)
            {
                if (reference[j] < max)
                {
                    double key = peer[alignment[i].second] - reference[j];
                    if (fabs(key) < 1 && (reference[j] != reference[alignment[i].first]))
                    {
                        s.push_back(make_pair(j, alignment[i].second));
                    }
                    if (mylib::data_stream::find_key(key, modifier) != -1 && (reference[j] != reference[alignment[i].first]))
                    {
                        if (mylib::data_stream::remove_false(key, modifier, modifier_mass))
                        {
                            m.push_back(make_pair(j, alignment[i].second));
                        }
                    }
                }
                else
                {
                    break;
                }
            }
            dup = 1;
        }
        else
        {
            dup++;
        }
    }

    for (i = 0; i < s.size(); i++)
    {
        alignment.push_back(s[i]);
    }
    for (i = 0; i < m.size(); i++)
    {
        alignment.push_back(m[i]);
    }
    std::sort(alignment.begin(), alignment.end(), mylib::data_stream::cmp1);
    alignment.erase(unique(alignment.begin(), alignment.end()), alignment.end());
}

void mylib::data_stream::align_modification_down(vector<double> reference, vector<double> peer,
                                                 std::vector<std::pair<long, long>> &alignment, std::vector<double> &modifier,
                                                 double &modifier_mass)
{ // 向下范围搜{
    int i = 0;
    vector<pair<long, long>>::iterator iter = alignment.begin();
    vector<pair<long, long>>::iterator iter_temp = alignment.begin();
    int r_size = reference.size() - 1;
    int p_size = peer.size() - 1;
    for (i = 0; i < alignment.size(); ++i)
    {
        iter_temp++;
        if (mylib::data_stream::find_key(peer[alignment[i].second] - reference[alignment[i].first], modifier) == -1 &&
            peer[alignment[i].second] - reference[alignment[i].first] > 2)
        {
            int j = i + 1;
            iter = iter_temp;
            while (j < alignment.size() && peer[alignment[i].second] - reference[alignment[j].first] > 0 &&
                   peer[alignment[i].second] - reference[alignment[j].first] <= modifier_mass + 1.007276) // 往下搜索{
            {
                iter++;
                if (mylib::data_stream::find_key(peer[alignment[i].second] - reference[alignment[j].first], modifier) != -1)
                {
                    std::pair<long, long> s = make_pair(alignment[j].first, alignment[i].second);
                    alignment.insert(iter, s);
                    i++;
                    // cout<<s.first<<"  " <<s.second<<endl;
                    break;
                }
                j++;
            }
        }
    }
}

void mylib::data_stream::align_modification_top(vector<double> reference, vector<double> peer,
                                                std::vector<std::pair<long, long>> &alignment, std::vector<double> &modifier,
                                                double &modifier_mass) // 向上范围搜
{
    int i = 0;
    vector<pair<long, long>>::iterator iter = alignment.begin();
    vector<pair<long, long>>::iterator iter_temp = alignment.begin();
    int r_size = reference.size() - 1;
    int p_size = peer.size() - 1;
    int j = 0;

    for (i = 0; i < alignment.size(); i++)
    {
        if ((peer[alignment[i].second] - peer[alignment[i - 1].second] >
             0.3)) // g::proc::find_key(peer[alignment[i].second]-reference[alignment[i].first],modifier)==-1
        {
            int j = i - 1;
            iter = iter_temp; // 指向当前元素
            iter--;
            while (j >= 0 && (peer[alignment[i].second] - reference[alignment[j].first] <=
                              modifier_mass + 1.007276)) // (peer[alignment[i].second]-reference[alignment[j].first] >0)
            {
                std::vector<std::pair<long, long>> fo;
                if (mylib::data_stream::find_key(peer[alignment[i].second] - reference[alignment[j].first], modifier) != -1 &&
                    (reference[alignment[j].first] != reference[alignment[i].first]))
                {
                    std::pair<long, long> s = make_pair(alignment[j].first, alignment[i].second);
                    alignment.insert(iter, s);
                    break;
                }
                j--;
            }
        }
        iter_temp++;
    }
}

int mylib::data_stream::find_key(double key, vector<double> &mod)
{
    for (int i = 0; i < mod.size(); i++)
    {
        if ((key >= mod[i] - 2) && (key <= mod[i] + 2))
        {
            return i;
        }
    }
    return -1;
}

void mylib::data_stream::Read_massdata(const char file_name[], std::vector<double> &peer, double &Pr_mass)
{ // 读取质谱数据
    std::string flag;
    std::string Id;
    std::string Scans;
    std::string Ret_time;
    std::string Activation;
    std::string Ms_one_id;
    std::string Ms_one_scan;
    std::string Pre_mz;
    std::string Pre_charge;
    std::string Pre_mass;
    std::string Pre_intensity;
    std::string Initial_data;
    ifstream ifs;
    std::string buf;
    ifs.open(file_name, ios::in);
    if (!ifs.is_open())
    {
        cout << "Error opening file";
        exit(1);
    }
    else
    {
        getline(ifs, flag);
        if (flag != "END IONS")
        {
            getline(ifs, Id);
            getline(ifs, Scans);
            getline(ifs, Ret_time);
            getline(ifs, Activation);
            getline(ifs, Ms_one_id);
            getline(ifs, Ms_one_scan);
            getline(ifs, Pre_mz);
            getline(ifs, Pre_charge);
            getline(ifs, Pre_mass);
            getline(ifs, Pre_intensity);
        }
        getline(ifs, flag);
        while (flag != "END IONS")
        {
            char data[30];
            for (int i = 0; flag[i] != '\t'; i++)
            {
                data[i] = flag[i];
            }
            peer.push_back(strtod(data, NULL));
            getline(ifs, flag);
        }
        ifs.close();
    }

    char c[100];
    for (int i = 0; Pre_mass[i] != '\0'; i++)
    {
        if (Pre_mass[i] == '=')
        {
            int k = 0;
            for (int j = i + 1; Pre_mass[j] != '\0'; j++, k++)
            {
                c[k] = Pre_mass[j];
            }
            break;
        }
    }
    Pr_mass = strtod(c, NULL);
    sort(peer.begin(), peer.end());
    // MergeSort(peer,0,peer.size()-1);
}

//----------------排序--------------//
void mylib::data_stream::Merge(vector<double> &arr, int l, int q, int r)
{
    int n = r - l + 1; // 临时数组存合并后的有序序列
    double *tmp = new double[n];
    int i = 0;
    double left = l;
    double right = q + 1;
    while (left <= q && right <= r)
        tmp[i++] = arr[left] <= arr[right] ? arr[left++] : arr[right++];
    while (left <= q)
        tmp[i++] = arr[left++];
    while (right <= r)
        tmp[i++] = arr[right++];
    for (int j = 0; j < n; ++j)
        arr[l + j] = tmp[j];
    delete[] tmp; // 删掉堆区的内存
}

void mylib::data_stream::MergeSort(vector<double> &arr, int l, int r)
{
    if (l == r)
        return; // 递归基是让数组中的每个数单独成为长度为1的区间
    double q = (l + r) / 2;
    MergeSort(arr, l, q);
    MergeSort(arr, q + 1, r);
    Merge(arr, l, q, r);
}

int mylib::data_stream::removeDuplicate(vector<double> &A, int n)
{
    for (int i = 0; i < n; i++)
    {
        double x = A[i]; // 暂存变量
        for (int j = i + 1; j < n; j++)
        {
            while (x == A[j])
            {
                for (int l = j; l < n; l++)
                {
                    A[l] = A[l + 1];
                }
                n--;
            }
        }
    }
    return n;
}

void mylib::data_stream::index_ions(const vector<double> &sub, vector<double> &index_y,
                                    std::vector<double> &modifier,
                                    double &modifier_mass) // 定位修饰偏差值
{
    int i = 0;
    for (i = 0; i < sub.size(); ++i)
    {
        // 判断该值小于修饰值总质量，同时为一定修饰组合的质量
        if ((sub[i] <= modifier_mass + 2) && (mylib::data_stream::find_key(sub[i], modifier) != -1))
        {
            index_y.push_back(sub[i]);
        }
        else if (fabs(sub[i]) <= 2) // 是否该值为2以下？
        {
            index_y.push_back(sub[i]);
        }
        else // 既不是质量偏移峰，也不是正确的峰
        {
            index_y.push_back(INT_MAX);
        }
    }
    //------------------过滤一下----------------//		//有些离子的组合，无法到修饰质量例如14*6 = 74 + ？ = 98
    for (i = 0; i < index_y.size(); ++i)
    {
        int temp;
        if (index_y[i] != INT_MAX && index_y[i] < modifier_mass - 2)
        {
            temp = mylib::data_stream::remove_false(index_y[i], modifier, modifier_mass);
            if (temp == 0)
            {
                index_y[i] = INT_MAX;
            }
        }
    }
    //
    if (modifier_mass < 0)
    {
        for (i = 0; i < sub.size(); ++i)
        {
            if (sub[i] < 0)
                if (fabs(sub[i]) <= fabs(modifier_mass) + 0.5)
                {
                    index_y[i] = sub[i];
                }
        }
    }
}

int mylib::data_stream::remove_false(double index, vector<double> &modifier, const double &modifier_mass) // 过滤掉不可能的值
{
    int i = 0, j = 0;
    if (index > modifier_mass - 2 && index < modifier_mass + 2)
    {
        return 1;
    }
    for (j = 0; j < modifier.size(); ++j)
    {
        if ((index + modifier[j] >= modifier_mass - 2) && (index + modifier[j] <= modifier_mass + 2))
        {
            return 1;
        }
        if (index + modifier[j] > modifier_mass)
        {
            return 0;
        }
    }
    return 0;
}

void mylib::data_stream::align_top_down(vector<double> reference, vector<double> peer,
                                        std::vector<std::pair<long, long>> &alignment, std::vector<double> &modifier,
                                        double &modifier_mass) // 范围搜索
{
    int i = 0;
    vector<pair<long, long>>::iterator iter = alignment.begin();
    vector<pair<long, long>>::iterator iter_temp = alignment.begin();
    int r_size = reference.size() - 1;
    int p_size = peer.size() - 1;
    int j = 0;
    int lenth = alignment.size() - 1;
    for (i = 0; i < lenth - 2; i++)
    { // 往下搜索，最后一个数，不做搜索
        if (alignment[i].second != alignment[i + 1].second)
        {
            int rflag = alignment[i].first + 1;
            while (peer[alignment[i].second] > reference[rflag])
            {
                if (mylib::data_stream::find_key(peer[alignment[i].second] - reference[rflag], modifier) != -1 &&
                    peer[alignment[i].second] - reference[rflag] <= modifier_mass + 2)
                {
                    if (mylib::data_stream::remove_false(peer[alignment[i].second] - reference[rflag], modifier, modifier_mass))
                    {
                        alignment.push_back(make_pair(rflag, alignment[i].second));
                    }
                }
                rflag++;
                if (peer[alignment[i].second] < reference[rflag])
                {
                    break;
                }
            }
        }
    }

    for (; lenth > 0; lenth--)
    { // 向上搜索
        if (alignment[lenth].second != alignment[lenth - 1].second)
        {
            int lflag = alignment[lenth].first - 1;
            while (peer[alignment[lenth].second] > reference[lflag])
            {
                if (mylib::data_stream::find_key(peer[alignment[lenth].second] - reference[lflag], modifier) != -1 &&
                    peer[alignment[lenth].second] - reference[lflag] <= modifier_mass + 2)
                {
                    if (mylib::data_stream::remove_false(peer[alignment[lenth].second] - reference[lflag], modifier,
                                                         modifier_mass))
                    {
                        alignment.push_back(make_pair(lflag, alignment[lenth].second));
                    }
                }
                lflag--;
                if (peer[alignment[lenth].second] - reference[lflag] > modifier_mass + 2)
                {
                    break;
                }
            }
        }
    }
    std::sort(alignment.begin(), alignment.end(), mylib::data_stream::cmp1);
}

bool mylib::data_stream::cmp1(pair<int, int> a, pair<int, int> b)
{
    if (a.first != b.first)
        return a.first < b.first;
    if (a.first == b.first)
        return a.second < b.second;
}

void mylib::data_stream::Get_ppm(const vector<double> &reference, const vector<double> &peer,
                                 const std::vector<std::pair<long, long>> &alignment,
                                 const std::vector<double> &index,
                                 std::vector<double> &modifier, std::vector<double> &ppm,
                                 double &modifier_mass)
{
    int i = 0;
    double tol_ppm = 0.00003;
    for (i = 0; i < index.size(); i++)
    {
        if (index[i] != INT_MAX && fabs(index[i]) > 2)
        { // index为修饰值时
            int k = mylib::data_stream::find_value(index[i], modifier, modifier_mass);
            if (k != -1)
            {
                double sub = peer[alignment[i].second] - reference[alignment[i].first];
                double tol_mass = (reference[alignment[i].first] + modifier[k]) * tol_ppm; // 容错值
                if (sub >
                    modifier[k])
                { // 当差值为正，即质谱位移往右偏差(两种情况，1.差一丢丢 2.容错设计。3.1.0024的设计)
                    double pm = fabs((reference[alignment[i].first] + modifier[k] + tol_mass) -
                                     (peer[alignment[i].second])) /
                                (reference[alignment[i].first] + modifier[k] + tol_mass) * 1000000;
                    double pm2 = fabs((reference[alignment[i].first] + modifier[k]) - (peer[alignment[i].second])) /
                                 (reference[alignment[i].first] + modifier[k]) * 1000000;
                    double pm3 =
                        fabs((reference[alignment[i].first] + modifier[k] + 1.0024) - (peer[alignment[i].second])) /
                        (reference[alignment[i].first] + modifier[k]) * 1000000;
                    double t = pm > pm2 ? (pm2 > pm3 ? pm3 : pm2) : (pm > pm3 ? pm3 : pm);
                    ppm.push_back(t);
                }
                else
                { // 当差值为负时，质量向左偏。（可能少H，可能仅仅只是偏差一点）
                    double pm = fabs((reference[alignment[i].first] + modifier[k] + tol_mass) -
                                     (peer[alignment[i].second] + 1.007276)) /
                                (reference[alignment[i].first] + modifier[k] + tol_mass) * 1000000; // 少H的容错设计
                    double pm0 = fabs((reference[alignment[i].first] + modifier[k]) -
                                      (peer[alignment[i].second] + 1.007276)) /
                                 (reference[alignment[i].first] + modifier[k]) * 1000000; // 少H
                    double pm1 = fabs((reference[alignment[i].first] + modifier[k]) - (peer[alignment[i].second])) /
                                 (reference[alignment[i].first] + modifier[k]) * 1000000; // 偏差一点点
                    double pm2 = fabs((reference[alignment[i].first] + modifier[k] + tol_mass) -
                                      (peer[alignment[i].second])) /
                                 (reference[alignment[i].first] + modifier[k] + tol_mass) * 1000000; // 容错设计
                    double t = pm < pm1 ? (pm < pm2 ? pm : pm2) : (pm1 < pm2 ? pm1 : pm2);
                    if (pm0 > t)
                    {
                        ppm.push_back(t);
                    }
                    else
                    {
                        ppm.push_back(pm0);
                    }
                }
            }
            else
            {
                ppm.push_back(INT_MAX);
            }
        }
        else if (index[i] != INT_MAX && fabs(index[i]) < 2)
        { // 当index值为小于2时
            if (index[i] > 0)
            { // 当质量为右偏移，那么只进行容错的设计,还要考虑多H的设计
                double tol_mass = reference[alignment[i].first] * tol_ppm;
                double pm = fabs(reference[alignment[i].first] - peer[alignment[i].second]) /
                            (reference[alignment[i].first]) * 1000000;
                double pmAdd_H = fabs(reference[alignment[i].first] + 1.007276 - peer[alignment[i].second]) /
                                 (reference[alignment[i].first] + 1.007276) * 1000000;
                double pm_tol_add = fabs(reference[alignment[i].first] + tol_mass - peer[alignment[i].second]) /
                                    (reference[alignment[i].first]) * 1000000;
                double pm_tol_sub = fabs(reference[alignment[i].first] - tol_mass - peer[alignment[i].second]) /
                                    (reference[alignment[i].first]) * 1000000;
                double t = pm_tol_add > pm_tol_sub ? (pm_tol_sub > pm ? pm : pm_tol_sub) : (pm_tol_add > pm ? pm : pm_tol_add);
                double t1 = t > pmAdd_H ? pmAdd_H : t;
                ppm.push_back(t1);
            }
            if (index[i] < 0)
            {
                double tol_mass = reference[alignment[i].first] * tol_ppm;
                double pm = fabs((reference[alignment[i].first] + tol_mass) - (peer[alignment[i].second] + 1.007276)) /
                            (reference[alignment[i].first] + tol_mass) * 1000000; // 少H的容错设计
                double pm0 = fabs((reference[alignment[i].first]) - (peer[alignment[i].second] + 1.007276)) /
                             (reference[alignment[i].first]) * 1000000; // 少H
                double pm1 = fabs((reference[alignment[i].first]) - (peer[alignment[i].second])) /
                             (reference[alignment[i].first]) * 1000000; // 偏差一点点
                double pm2 = fabs((reference[alignment[i].first] + tol_mass) - (peer[alignment[i].second])) /
                             (reference[alignment[i].first] + tol_mass) * 1000000; // 容错设计
                double t = pm < pm1 ? (pm < pm2 ? pm : pm2) : (pm1 < pm2 ? pm1 : pm2);
                if (pm0 > t)
                {
                    ppm.push_back(t);
                }
                else
                {
                    ppm.push_back(pm0);
                }
            }
        }
        else
        {
            ppm.push_back(INT_MAX);
        }
    }
    if (modifier_mass < 0)
    {
        for (i = 0; i < index.size(); ++i)
        {
            if (ppm[i] < 15)
            {
                continue;
            }
            if (fabs(index[i]) > 0.9 && index[i] < 0 && fabs(index[i]) < fabs(modifier_mass) + 0.5)
            {
                double tol_ppm = 0.00003;
                double DeH = 1.007825;
                double tol_mass = reference[alignment[i].first] * tol_ppm;
                double pm =
                    (reference[alignment[i].first] - (peer[alignment[i].second] + round(fabs(index[i])) * DeH)) /
                    (reference[alignment[i].first]) * 1000000;
                // cout<<"index = "<<index[i]<<" round = "<<round(fabs(index[i])) * DeH<<"  pm = "<<pm<<endl;
                double pm2 = (reference[alignment[i].first] - (peer[alignment[i].second] + round(fabs(index[i])) * DeH +
                                                               tol_ppm * reference[alignment[i].first])) /
                             (reference[alignment[i].first]) * 1000000;
                double t = fabs(pm) > fabs(pm2) ? fabs(pm2) : fabs(pm);
                ppm[i] = t;
            }
        }
    }
}

int mylib::data_stream::find_value(double key, std::vector<double> &modifier, double &modifier_mass)
{ // 寻找修饰值
    int i = 0;
    double temp = key - modifier[0];
    int flag = 0;
    if (fabs(key) < 2)
    {
        return -1;
    }
    for (i = 1; i < modifier.size(); ++i)
    {
        if (temp > fabs(key - modifier[i]))
        {
            temp = fabs(key - modifier[i]);
            flag = i;
        }
        if (fabs(modifier[i] - modifier_mass) < 0.0000001)
        {
            return flag;
        }
    }
    return -1;
}

void mylib::data_stream::loca_ions(vector<node> &IonsCN, vector<node> &MergIonsNC,
                                   vector<locat_ions> &location,
                                   vector<double> peer)
{
    int i = 0;
    for (i = 0; i < peer.size(); i++)
    {
        locat_ions t;
        t.c = "x";
        t.n = "x";
        location.push_back(t);
    }
    for (i = 0; i < MergIonsNC.size(); i++)
    {
        location[MergIonsNC[i].mono_id].c = "Y" + to_string(MergIonsNC[i].thoe_id + 1);
    }
    for (i = 0; i < IonsCN.size(); i++)
    {
        location[IonsCN[i].mono_id].n = "B" + to_string(IonsCN[i].thoe_id + 1);
    }
}

void mylib::data_stream::ModifyMassBiger0(vector<node> &IonsCN, vector<node> &IonsNC,
                                          vector<node> &MergIonsCN, vector<node> &MergIonsNC,
                                          const vector<pair<long, long>> &AlignmentCN,
                                          const vector<pair<long, long>> &AlignmentNC,
                                          const vector<double> &ReferenceCN, const vector<double> &ReferenceNC,
                                          const vector<double> &peer,
                                          const vector<double> &IndexCN, const vector<double> &IndexNC,
                                          const vector<double> &ppmCN, const vector<double> &ppmNC,
                                          vector<modi> &Ptms,
                                          const vector<double> &modify_table,
                                          int lenth, double ModifyMass)
{
    // 1.以NC端作为参照，搜索所有可行值
    mylib::data_stream::Get_all_ions(MergIonsNC, IonsNC, AlignmentNC, ReferenceNC, peer, IndexNC, ppmNC);
    if (IonsNC.size() != 0)
    {
        // NC端找到了修饰离子
        // 2.利用修饰离子，获取修饰区间(参考NC端得出最终的MergIonsNC ， 并且得出Ptms的修饰区间)
        mylib::data_stream::Get_mod_ions_loca(IonsNC, MergIonsNC, Ptms, AlignmentNC, modify_table, ModifyMass);
        std::map<int, bool> ItMap;
        for (int i = 0; i < MergIonsNC.size(); ++i)
        {
            ItMap[MergIonsNC[i].mono_id] = true;
        }
        if (Ptms.size() == 0)
        {
            // cout<<"NC端有修饰离子,修饰区间为0"<<endl;
            return;
        }
        else
        {
            // cout<<"寻找另一段离子"<<endl;
            // 3.利用修饰区间寻找另一端离子,取出所有可行CN Ions离子
            for (int i = 0; i < AlignmentCN.size(); ++i)
            {
                if (ppmCN[i] < 15 && !ItMap.count(AlignmentCN[i].second))
                {
                    node t;
                    t.thoe_id = AlignmentCN[i].first;
                    t.mono_id = AlignmentCN[i].second;
                    t.index = IndexCN[i];
                    t.ppm = ppmCN[i];
                    MergIonsCN.push_back(t);
                }
            }
            if (MergIonsCN.size() == 0)
                return;
            // g::proc::Get_other_ions(AlignmentCN,Ptms,IonsCN,IndexCN,ppmCN,lenth);
            mylib::data_stream::GetAnotherIons(AlignmentCN, Ptms, MergIonsCN, IonsCN, IndexCN, ppmCN, lenth);
        }
    }
    else
    {
        // NC端找不到修饰离子（有可能修饰在最后一截）
        if (MergIonsNC.size() > 0)
        {
            // 1.首先判定掉H行为,若掉N个H则 ModifyMass = (ModifyMass + 1.0078 * N),并且将所得 ModifyMass 制在最后一段肽段上
            vector<modi> Dyhe;
            mylib::data_stream::IfDyHe(MergIonsNC, Dyhe, ModifyMass);
            for (int l = 0; l < Dyhe.size(); ++l)
            {
                Ptms.push_back(Dyhe[l]);
            }
            // 2.将修饰区间安在最一段区间,若安置不上则Ptms清空.
            // 只有当NC端存在修饰区间时,才寻找另外一端离子。
            mylib::data_stream::AssumePtmInLast(MergIonsNC, ReferenceNC, Ptms, ModifyMass);
            std::map<int, bool> ItMap;
            for (int i = 0; i < MergIonsNC.size(); ++i)
            {
                ItMap[MergIonsNC[i].mono_id] = true;
            }
            if (Ptms.size() != 0)
            {
                for (int i = 0; i < AlignmentCN.size(); ++i)
                {
                    if (ppmCN[i] < 15 && !ItMap.count(AlignmentCN[i].second))
                    {
                        node t;
                        t.thoe_id = AlignmentCN[i].first;
                        t.mono_id = AlignmentCN[i].second;
                        t.index = IndexCN[i];
                        t.ppm = ppmCN[i];
                        MergIonsCN.push_back(t);
                    }
                }
                if (MergIonsCN.size() == 0)
                    return;
                // g::proc::Get_other_ions(AlignmentCN,Ptms,IonsCN,IndexCN,ppmCN,lenth);
                mylib::data_stream::GetAnotherIons(AlignmentCN, Ptms, MergIonsCN, IonsCN, IndexCN, ppmCN, lenth);
            }
            else
            {
                return;
            }
        }
        else
        { // if(MergIonsNC.size() = 0 && IonsNC.size == 0)
            // IonsNC.clear();
            // 这里的逻辑是，直接参考CN端，取峰，取修饰
            mylib::data_stream::Get_all_ions(IonsCN, MergIonsCN, AlignmentCN, ReferenceCN, peer, IndexCN, ppmCN);
            if (IonsCN.size() == 0)
                return;
            if (MergIonsCN.size() > 0)
            {
                mylib::data_stream::Get_mod_ions_loca(MergIonsCN, IonsCN, Ptms, AlignmentCN, modify_table, ModifyMass);
            }
            else
            {
                mylib::data_stream::AssumePtmInLast(IonsCN, ReferenceCN, Ptms, ModifyMass);
            }
        }
    }
}

bool ModifyMassSmallerCmp(node a, node b)
{
    return a.thoe_id < b.thoe_id;
}

static bool myfunc(node i, node j)
{
    return i.thoe_id == j.thoe_id;
}

void mylib::data_stream::ModifyMassSmaller0(vector<node> &IonsCN, vector<node> &IonsNC,
                                            vector<node> &MergIonsCN, vector<node> &MergIonsNC,
                                            const vector<pair<long, long>> &AlignmentCN,
                                            const vector<pair<long, long>> &AlignmentNC,
                                            const vector<double> &ReferenceCN, const vector<double> &ReferenceNC,
                                            const vector<double> &peer,
                                            const vector<double> &IndexCN, const vector<double> &IndexNC,
                                            const vector<double> &ppmCN, const vector<double> &ppmNC,
                                            vector<modi> &Ptms,
                                            int lenth, double ModifyMass)
{
    int i = 0;
    vector<node> merge;
    // 1.取出两端离子
    // 这里只考虑到掉H的情况-------------------------需要修改-------------------//21.8.25
    for (i = 0; i < AlignmentCN.size(); ++i)
    {
        if (fabs(IndexCN[i]) < fabs(ModifyMass) && ppmCN[i] < 15)
        {
            node t;
            t.thoe_id = AlignmentCN[i].first;
            t.mono_id = AlignmentCN[i].second;
            t.index = IndexCN[i];
            t.ppm = ppmCN[i];
            IonsCN.push_back(t);
            if (fabs(IndexCN[i]) > 0.9 && IndexCN[i] < 0)
            {
                merge.push_back(t);
            }
        }
        else
        {
            continue;
        }
    }
    for (i = 0; i < AlignmentNC.size(); ++i)
    {
        if (fabs(IndexNC[i]) < fabs(ModifyMass) && ppmNC[i] < 15)
        {
            node t;
            t.thoe_id = AlignmentNC[i].first;
            t.mono_id = AlignmentNC[i].second;
            t.index = IndexNC[i];
            t.ppm = ppmNC[i];
            MergIonsNC.push_back(t);
            if (fabs(IndexNC[i]) > 0.9 && IndexNC[i] < 0)
            {
                t.thoe_id = lenth - AlignmentNC[i].first;
                merge.push_back(t);
            }
        }
        else
        {
            continue;
        }
    }
    if (merge.size() == 0)
    {
        modi Temp;
        if (IonsCN.size() != 0)
        {
            merge = IonsCN;
            Temp.first = IonsCN.back().thoe_id;
        }
        else
        {
            Temp.first = 0;
        }
        for (i = 0; i < MergIonsNC.size(); ++i)
        {
            node t;
            t.thoe_id = lenth - MergIonsNC[i].thoe_id;
            t.mono_id = MergIonsNC[i].mono_id;
            t.index = MergIonsNC[i].index;
            t.ppm = MergIonsNC[i].ppm;
            merge.push_back(t);
        }
        sort(merge.begin(), merge.end(), ModifyMassSmallerCmp);
        merge.erase(unique(merge.begin(), merge.end(), myfunc), merge.end());
        for (i = 0; i < merge.size(); ++i)
        {
            if (merge[i].thoe_id > Temp.first)
            {
                Temp.second = merge[i].thoe_id;
                Temp.mod_mass = ModifyMass;
                Ptms.push_back(Temp);
                break;
            }
        }
    }
    else
    {
        sort(merge.begin(), merge.end(), ModifyMassSmallerCmp);
        merge.erase(unique(merge.begin(), merge.end(), myfunc), merge.end());
        // 2.定位修饰区间
        int MaxDyHe = round(fabs(ModifyMass));
        double H = 1.007825;
        double mass = 0;
        for (i = 0; i < merge.size(); ++i)
        {
            if (Ptms.size() == 0)
            {
                modi t;
                t.first = 0;
                t.second = merge[i].thoe_id;
                t.mod_mass = round(fabs(merge[i].index)) * H * -1;
                mass += t.mod_mass;
                Ptms.push_back(t);
            }
            else if (fabs(merge[i].index) - fabs(mass) > 0.9)
            {
                modi t;
                t.first = merge[i - 1].thoe_id;
                t.second = merge[i].thoe_id;
                t.mod_mass = round(fabs(fabs(merge[i].index) - fabs(mass))) * H * -1;
                mass += t.mod_mass;
                Ptms.push_back(t);
            }
            // if(fabs(mass) - fabs(ModifyMass) < 0.5) break;
        }
        if (fabs(ModifyMass) - fabs(mass) > 0.9)
        {
            modi t;
            t.first = merge.back().thoe_id;
            t.second = ReferenceCN.size() - 1;
            t.mod_mass = round(fabs(fabs(ModifyMass) - fabs(mass))) * H * -1;
            mass += t.mod_mass;
            Ptms.push_back(t);
        }
    }
}

// 2021.7.24 星期六待完善,已完善
void mylib::data_stream::ions_location(vector<node> &IonsCN, vector<node> &IonsNC,
                                       vector<node> &MergIonsCN, vector<node> &MergIonsNC,
                                       const vector<pair<long, long>> &AlignmentCN, const vector<pair<long, long>> &AlignmentNC,
                                       const vector<double> &ReferenceCN, const vector<double> &ReferenceNC,
                                       const vector<double> &peer,
                                       const vector<double> &IndexCN, const vector<double> &IndexNC,
                                       const vector<double> &ppmCN, const vector<double> &ppmNC,
                                       vector<modi> &Ptms, const vector<double> &modify_table,
                                       int lenth, double ModifyMass)
{
    int i = 0;
    if (ModifyMass > 0)
    { //-----------ModifyMass > 0
        // cout<<"ModifyMassBiger0 "<<endl;
        mylib::data_stream::ModifyMassBiger0(IonsCN, IonsNC,
                                             MergIonsCN, MergIonsNC, AlignmentCN, AlignmentNC,
                                             ReferenceCN, ReferenceNC, peer,
                                             IndexCN, IndexNC, ppmCN, ppmNC, Ptms,
                                             modify_table,
                                             lenth, ModifyMass);
    }
    else
    { // modifymass <= 0
        if (ModifyMass == 0)
        { //-----------ModifyMass = 0
            // cout<<"ModifyMass = 0  "<<endl;
            mylib::data_stream::Get_all_ions(MergIonsNC, IonsNC, AlignmentNC, ReferenceNC, peer, IndexNC, ppmNC);
            mylib::data_stream::Get_all_ions(IonsCN, MergIonsCN, AlignmentCN, ReferenceCN, peer, IndexCN, ppmCN);
        }
        else if (ModifyMass < 0)
        { //-----------ModifyMass < 0
            // cout<<"ModifyMass < 0  "<<endl;
            mylib::data_stream::ModifyMassSmaller0(IonsCN, IonsNC,
                                                   MergIonsCN, MergIonsNC, AlignmentCN, AlignmentNC,
                                                   ReferenceCN, ReferenceNC, peer,
                                                   IndexCN, IndexNC, ppmCN, ppmNC, Ptms, lenth, ModifyMass);
        }
    } // else End .
}

void mylib::data_stream::AnalysisCNIons(vector<node> &MergIonsCN, vector<node> &Ions, vector<modi> &Ptms,
                                        vector<double> &Reference, double ModifyMass)
{
    if (MergIonsCN.size() != 0)
    {
        vector<node> Copy;
        double mass = 0;
        for (int i = 0; i < MergIonsCN.size(); ++i)
        {
            if (i == 0)
            {
                Copy.push_back(MergIonsCN[i]);
            }
            else
            {
                if (MergIonsCN[i].index - Copy.back().index > 12)
                {
                    Copy.push_back(MergIonsCN[i]);
                }
            }
        }
        for (int i = 0; i < Copy.size(); ++i)
        {
            if (Ptms.size() == 0)
            {
                modi t;
                t.first = Ions.back().thoe_id;
                t.second = Copy[i].thoe_id;
                t.mod_mass = Copy[i].index;
                Ptms.push_back(t);
            }
            else
            {
                modi t;
                t.first = Ptms.back().first;
                t.second = Copy[i].thoe_id;
                t.mod_mass = Copy[i].index;
                Ptms.push_back(t);
            }
        }
    }
    else
    {
        modi t;
        t.first = Ions.back().thoe_id;
        t.second = Reference.size() - 1;
        t.mod_mass = ModifyMass;
        Ptms.push_back(t);
    }
}

// 将修饰固定在最后一段
void mylib::data_stream::AssumePtmInLast(vector<node> &MergIons, const vector<double> &Reference,
                                         vector<modi> &Ptms, double ModifyMass)
{
    modi temp;
    if (MergIons[MergIons.size() - 1].thoe_id == Reference.size() - 1)
    {
        // cout<<"最后一端的修饰无法固定。"<<endl;
        Ptms.clear();
        return;
    }
    if (MergIons.size() < 0)
        return;
    temp.first = MergIons[MergIons.size() - 1].thoe_id;
    temp.second = Reference.size() - 1;
    temp.mod_mass = ModifyMass;
    Ptms.push_back(temp);
    // cout<<"temp.first = "<<temp.first<<" temp.scond = "<<temp.second <<"temp.modimass = "<<temp.mod_mass<<endl;
}

// 寻找少H的修饰区间
void mylib::data_stream::IfDyHe(vector<node> &MergIonsNC, vector<modi> &Dyhe, double &ModifyMass)
{
    double first;
    double last;
    double He = 1.007825;
    double Temp = 0;

    for (int i = 1; i < MergIonsNC.size(); ++i)
    {
        // 太小的值直接过滤
        if (fabs(MergIonsNC[i].index) <= 0.1 || MergIonsNC[i].index > 0)
            continue;
        // 累计掉H.那如果一次性出现3个H？
        int Ze = round(fabs(MergIonsNC[i].index)) > 1 ? round(fabs(MergIonsNC[i].index)) : 1;
        // 首先，判定该修饰值是否满足少H的条件（fabs (fabs(MergIonsNC[i].index) - Ze * He ) < 0.1）
        // 其次判定，是否要比我之前少H的修饰量大fabs (fabs(MergIonsNC[i].index) - Temp )>0.8。
        if (fabs(fabs(MergIonsNC[i].index) - Ze * He) < 0.1 && fabs(MergIonsNC[i].index) - Temp > 0.8)
        {
            if (Dyhe.size() == 0)
            {
                modi t;
                t.first = MergIonsNC[i - 1].thoe_id;
                t.second = MergIonsNC[i].thoe_id;
                t.mod_mass = -1 * He * Ze;
                Dyhe.push_back(t);
                Temp += Ze * He;
            }
            else
            {
                modi t;
                t.first = Dyhe.back().first;
                t.second = MergIonsNC[i].thoe_id;
                t.mod_mass = -1 * (He * Ze - Temp);
                Dyhe.push_back(t);
                Temp = Ze * He;
            }
        }
    }
    // 将少的修饰值添加回总质量
    ModifyMass += Temp;
    // cout<<"ModifyMass = "<<ModifyMass<<endl;
    // cout<<"少H的查找结果如下:"<<endl;
    // for(int i = 0 ; i < Dyhe.size() ; ++ i )
    // {
    // 	cout<<Dyhe[i].first<<"  "<<Dyhe[i].second<<" "<<Dyhe[i].mod_mass<<endl;
    // }
}

// 寻找另外一端的离子（差值为O）
void mylib::data_stream::GetOtherIonsZero(const vector<pair<long, long>> &alignment, vector<node> &ions,
                                          vector<double> &index, vector<double> &ppm)
{
    for (int i = 0; i < alignment.size(); ++i)
    {
        if (ppm[i] <= 15)
        {
            node t;
            t.thoe_id = alignment[i].first;
            t.mono_id = alignment[i].second;
            t.index = index[i];
            t.ppm = ppm[i];
            ions.push_back(t);
        }
    }
}

void mylib::data_stream::Get_other_ions(const vector<pair<long, long>> &alignment,
                                        vector<modi> &mod,
                                        vector<node> &ions,
                                        const vector<double> &index, const vector<double> &ppm,
                                        int lenth)
{
    vector<modi> copy = mod;
    cout << copy.size() << "  " << copy.capacity() << endl;
    // cout<<"修饰区间:\n";
    // for(int i = 0 ; i < mod.size() ; ++i )
    // {
    // 	cout<<mod[i].first <<" "<<mod[i].second<<" " << mod[i].mod_mass<<endl;
    // }
    reverse(copy.begin(), copy.end());
    // cout<<"修饰区间recovey:\n";
    for (int i = 0; i < copy.size(); i++)
    {
        copy[i].first += copy[i].second;
        copy[i].second = copy[i].first - copy[i].second;
        copy[i].first -= copy[i].second;

        copy[i].first = lenth - copy[i].first;
        copy[i].second = lenth - copy[i].second;
        // cout<<copy[i].first<<" "<<copy[i].second <<" " << copy[i].mod_mass<<endl;
    }
    int flag = 0;
    int maxsize = copy.size();
    double mass = copy.front().mod_mass;
    map<long, int> mapi;
    for (int i = 0; i < alignment.size(); ++i)
    {
        if (alignment[i].first < copy[flag].second)
        {
            // 如果在第一个区间内为 0 离子
            if (flag == 0 && fabs(index[i]) <= 2 && ppm[i] <= 15)
            {
                // 如果之前都没有离子取值
                if (ions.size() == 0 && !mapi.count(alignment[i].second))
                {
                    node t;
                    t.index = index[i];
                    t.thoe_id = alignment[i].first;
                    t.mono_id = alignment[i].second;
                    t.ppm = ppm[i];
                    ions.push_back(t);
                    mapi[alignment[i].second] = 1;
                    copy[flag].first = alignment[i].first + 1;
                }
                else if (fabs(ions.back().index) <= 2 && !mapi.count(alignment[i].second))
                {
                    // 如果之前取了离子值，看是否是相同类型的离子
                    node t;
                    t.index = index[i];
                    t.thoe_id = alignment[i].first;
                    t.mono_id = alignment[i].second;
                    t.ppm = ppm[i];
                    ions.push_back(t);
                    mapi[alignment[i].second] = 1;
                    copy[flag].first = alignment[i].first + 1;
                }
                // 如果不是0离子
            }
            else if (!mapi.count(alignment[i].second) && fabs(index[i] - mass) <= 4 && ppm[i] <= 15 &&
                     alignment[i].first >= copy[flag].first)
            {
                node t;
                t.index = index[i];
                t.thoe_id = alignment[i].first;
                t.mono_id = alignment[i].second;
                t.ppm = ppm[i];
                ions.push_back(t);
                mapi[alignment[i].second] = 1;
                copy[flag].first = alignment[i].first;
            }
        }
        else if (flag < maxsize)
        { // 也就是说进入到下一个区间
            ++flag;
            mass += copy[flag].mod_mass;
        }
        else if (flag == maxsize)
        {
            if (alignment[i].first > copy.back().second && !mapi.count(alignment[i].second))
            {
                if (fabs(index[i] - mass) < 4 && ppm[i] <= 15)
                {
                    node t;
                    t.index = index[i];
                    t.thoe_id = alignment[i].first;
                    t.mono_id = alignment[i].second;
                    t.ppm = ppm[i];
                    ions.push_back(t);
                    mapi[alignment[i].second] = 1;
                }
            }
        }
    }
    // cout<<" b ions size = "<<ions.size()<<endl;
    // cout<<"---------"<<endl;
    // for(int i = 0 ; i < copy.size() ; i++)
    // {
    // 	cout<<copy[i].first<<" "<<copy[i].second <<" " << copy[i].mod_mass<<endl;
    // }
    mod.swap(copy);
}

bool JugeModIsRight(double mass, const vector<double> &mod)
{
    for (int i = 0; i < mod.size(); ++i)
    {
        if (fabs(mass - mod[i]) < 1)
        {
            return true;
        }
    }
    return false;
}

// 取修饰值的策略
void mylib::data_stream::Get_mod_ions_loca(vector<node> &Ions, vector<node> &MergeIons, vector<modi> &mod,
                                           const std::vector<std::pair<long, long>> &alignment,
                                           const vector<double> &modify_table,
                                           double modifier_mass)
{
    int i = 0;
    vector<node> temp_best; // 将最优的离子取值保存
    node flag;
    vector<node> temp; // 将初步取的离子数值保存。去除同大小的修饰
    vector<modi> modf;
    // 循环取最优值
    for (int j = 0; j < Ions.size(); j++)
    {
        // cout<<" "<<flag.thoe_id<<" "<<flag.mono_id<<" "<<flag.index<<endl;
        // node flag = Ions[j];
        // vector<node>  temp;   //将初步取的离子数值保存。去除同大小的修饰
        // vector<modi>  modf;
        flag = Ions[j];
        temp.push_back(flag);
        for (i = j + 1; i < Ions.size(); i++) // 贪婪取值(往下搜寻，当有相等的加入，没有则加入大值)
        {
            if (fabs(Ions[i].index - flag.index) < 2)
            {
                temp.push_back(Ions[i]);
                flag = Ions[i];
                continue;
            }
            // if(Ions[i].index > flag.index && fabs(Ions[i].index - flag.index) > 10 &&Ions[i].thoe_id!=flag.thoe_id) {
            if (Ions[i].thoe_id != flag.thoe_id && Ions[i].index > flag.index)
            {
                double mass = Ions[i].index - flag.index;
                // cout<<"mass = "<<mass<<" "<<Ions[i].index<<" "<<flag.index<<endl;
                if (JugeModIsRight(mass, modify_table))
                {
                    // cout<<Ions[i].index<<" "<<flag.index<<endl;
                    temp.push_back(Ions[i]);
                    flag = Ions[i];
                }
                continue;
            }
        }
        // 定位第一个修饰的区间
        if (MergeIons.size() != 0)
        {
            modi r;
            r.first = MergeIons.back().thoe_id;
            r.second = temp.front().thoe_id;
            r.mod_mass = temp.front().index;
            modf.push_back(r);
        }
        else
        {
            modi r;
            r.first = 0;
            r.second = temp.front().thoe_id;
            r.mod_mass = temp.front().index;
            modf.push_back(r);
        }
        for (i = 1; i < temp.size(); i++)
        {
            double sub_ms = temp[i].index - temp[i - 1].index;
            // if(sub_ms > 12 ) {//相差一个修饰
            if (JugeModIsRight(sub_ms, modify_table) && fabs(sub_ms) > 0.5)
            {
                modi r;
                r.first = temp[i - 1].thoe_id;
                r.second = temp[i].thoe_id;
                r.mod_mass = sub_ms;
                modf.push_back(r);
            }
            else if (fabs(sub_ms) < 2 && modf.size() > 0)
            {
                // modf[modf.size()-1].second = temp[i].thoe_id;
            }
            // 在最后一段肽段定位修饰
            if (i == temp.size() - 1 && modf.size() > 0)
            {
                if (temp[i].index + 2 < modifier_mass)
                {
                    modi r;
                    r.first = modf[modf.size() - 1].second;
                    r.second = alignment[alignment.size() - 1].first;
                    r.mod_mass = modifier_mass - temp[i].index;
                    modf.push_back(r);
                }
            }
        }

        if (temp.size() > temp_best.size())
        {
            temp_best = temp;
            mod = modf;
        }
        temp.clear();
        modf.clear();
    } // 寻找最优结束

    for (i = 0; i < temp_best.size(); i++) // 得到带修饰的y ions
    {
        MergeIons.push_back(temp_best[i]);
    }
    // for(i = 0 ; i<temp_best.size();i++)          //得到带修饰的y ions
    // {
    //     cout<<"修饰离子 : "<<temp_best[i].thoe_id<<" "<<temp_best[i].mono_id<<" "<<temp_best[i].index<<endl;
    // }
    // cout<<"================================"<<endl;
    // for(i = 0 ; i <mod.size() ; i++) {
    // 	cout<<" Frist = "<<mod[i].first<<" Second = "<<mod[i].second<<" mass = "<<mod[i].mod_mass<<endl;
    // }
}

// void g::proc::Process_mod_ions(vector<node> & c_ions,vector<node> & merg_ions,vector<modi> & modf,
// 								std::vector<std::pair<long,long> >& alignment,double modifier_mass)   //贪婪策略取值
// {
//     int i = 0 ;
//     node flag = c_ions[0];
//     vector<node>  temp;   //将初步取的离子数值保存。
//     temp.push_back(flag);
//     for(i = 1 ; i<c_ions.size() ; i++)      //贪婪取值(往下搜寻，当有相等的加入，没有则加入大值)
//     {
//         if(c_ions[i].index>flag.index-2&&c_ions[i].index<flag.index+2)
//         {
//             temp.push_back(c_ions[i]);
//             flag = c_ions[i];
//             continue;
//         }
//         if(c_ions[i].index>flag.index&&c_ions[i].thoe_id!=flag.thoe_id)
//         {
//             temp.push_back(c_ions[i]);
//             flag = c_ions[i];
//             continue;
//         }
//     }
// 	cout<<"----temp_ions-------"<<endl;
// 	for(i=0;i<temp.size();i++)
// 	{
// 		cout<<temp[i].thoe_id<<" || "<<temp[i].mono_id<<endl;
// 	}
// 	cout<<"----temp_ions-------"<<endl;
//     if(merg_ions.size()!=0)				//定位第一个修饰的区间
//     {
//         modi r ;
//         r.first= merg_ions.back().thoe_id;
//         r.second = temp.front().thoe_id;
//         r.mod_mass = temp.front().index;
//         modf.push_back(r);
//     }
//     else
//     {
//         modi r ;
//         r.first= 0;
//         r.second = temp.front().thoe_id;
//         r.mod_mass = temp.front().index;
//         modf.push_back(r);
//     }
//     for(i = 1 ; i < temp.size() ; i++)      //取修饰值区间
//     {
//         if(temp[i].index - temp[i-1].index >12)		//相差一个修饰
//         {
//             modi r ;
//             r.first = temp[i-1].thoe_id;
//             r.second = temp[i].thoe_id;
//             r.mod_mass = temp[i].index - temp[i-1].index;
//             modf.push_back(r);
//         }
// 		if(fabs(temp[i].index - temp[i-1].index)<2)		//相差一个H？
// 		{
// 			if(modf.size()>=1)
// 			{
// 				modf[i-1].second = temp[i].thoe_id;
// 			}
// 		}
//     }

//     for(i = 0 ; i<temp.size();i++)          //得到带修饰的y ions
//     {
//         merg_ions.push_back(temp[i]);
//     }
// }

// 打分策略
// 1、得到该段蛋白序列的完全匹配峰和修饰匹配峰
// 2、占比 z,p
// 3、取出连续峰和不连续峰
// 4、每段连续的峰值，按照1.0,1.1,1.2,1.3,依例添加其分值。完全匹配的则1.0,1.2,1.4，添加
// void g::proc::Get_SCORE(vector<node> &ions, vector<double> &refe_score,
//                         vector<double> &peer_score, double &score,double base_score) {
//     int i = 0;
//     double p = 0.0;    //-标识修饰离子个数
//     double z = 0.0;    //->标识完全对准个数
//     for (i = 0; i < ions.size(); ++i) {
//         if (fabs(ions[i].index) < 0.1) {
//             ++z;
//         } else {
//             ++p;
//         }
//     }
//     p = 1 - (p / (ions.size()));        //-----缩小修饰得分比例
//     z = 1 + (z / (ions.size()));        //-----扩大完全对准得分比例
//     for (i = 0; i < ions.size(); i++) {
//         double t = 0;
//         if (fabs(ions[i].index) < 0.1) {
//             t = fabs(fabs(refe_score[ions[i].thoe_id]) + fabs(peer_score[ions[i].mono_id])) * (z + 0.1);
//         } else {
//             t = fabs(fabs(refe_score[ions[i].thoe_id]) + fabs(peer_score[ions[i].mono_id])) * (p + 0.1);
//         }
//         score += t;
//     }
//     if (p == 1) {
//         score /= 4;
//     }
//     // cout<<"ions size = "<<ions.size()<<endl;
//     // cout<<"score = "<<score<<endl;
// }

// 打分策略news
// 1、得到该段蛋白序列的完全匹配峰和修饰匹配峰
// 2、占比 z,p
// 3、取出连续峰和不连续峰
// 4、每段连续的峰值，按照1.0,1.1,1.2,1.3,依例添加其分值。完全匹配的则1.0,1.2,1.4，添加
void mylib::data_stream::Get_SCORE(vector<node> &ions_r, vector<double> &refe_score, const vector<double> &mod_arr,
                                   vector<double> &peer_score,
                                   double &score, double base_score,
                                   double mod_mass)
{
    vector<node> ions;
    double count = 0;
    double punish_nub = 0;
    if (fabs(mod_mass) > 500.0 / 2.0)
    { // 如果修饰质量太大，可能误差就大，故适当降低
        for (int pi = 0; pi < ions_r.size(); ++pi)
        {
            if (!mylib::data_stream::if_is_error_ions(mod_arr, ions_r[pi]))
            {
                count++;
                continue;
            }
            ions.push_back(ions_r[pi]);
        }
        punish_nub = 1 - (count / ions.size());
    }
    else
    {
        ions = ions_r;
        punish_nub = 1;
    }

    if (ions.size() == 0)
    {
        score = 0;
        return;
    }
    int i = 0;
    double p = 0.0;                          // 标识修饰离子个数
    double z = 0.0;                          // 标识完全对准个数
    vector<double> mutScore(ions.size(), 0); // 扩大系数
    vector<int> matchIndex(ions.size(), 0);  // 标识是否为完全匹配
    vector<int> matchSeq(ions.size(), 0);    // 标识是否连续
    map<int, int> rm_dup;
    // 标识出完全匹配的离子和修饰离子matchIndex
    // 得到完全匹配离子的扩大系数和修饰离子的扩大系数
    for (i = 0; i < ions.size(); i++)
    {
        if (!rm_dup.count(ions[i].thoe_id))
        {
            rm_dup[ions[i].thoe_id] = 1;
        }
        if (fabs(ions[i].index) < 1.1)
        {
            matchIndex[i] = 1;
            ++z;
        }
        else
        {
            ++p;
        }
    }
    // 基础得分系数
    p = 1 - (p / (ions.size())) + base_score; // 缩小修饰得分比例
    z = 1 + (z / (ions.size())) + base_score; // 扩大完全对准得分比例
    // 标识出连续匹配的离子matchSeq
    for (i = 0; i < ions.size() - 1; i++)
    {
        if (fabs(ions[i].thoe_id - ions[i + 1].thoe_id) == 1)
        {
            matchSeq[i] = 1;
            matchSeq[i + 1] = 1;
        }
    }
    //    p = 1 - (p / (ions.size())) + 0.1;        //缩小修饰得分比例
    //    z = 1 + (z / (ions.size())) + 0.1;        //扩大完全对准得分比例
    // 当只有完全匹配的离子时，判定是否存在连续匹配离子
    if ((1 - (p / (ions.size()))) == 0)
    {
        int lo = 0;
        for (i = 0; i < matchSeq.size(); i++)
        {
            if (matchSeq[i] == 1)
            {
                lo = 1;
                break;
            }
        }
        if (!lo)
        {
            z -= (z / (ions.size()));
        }
    }
    // 设置连续匹配系数mutScore
    for (i = 0; i < ions.size(); i++)
    {
        if (matchSeq[i] == 0)
        {
            if (matchIndex[i])
            {
                mutScore[i] += z;
            }
            else
            {
                mutScore[i] += p;
            }
        }
        else if (matchSeq[i] == 1)
        {
            double log = 0;
            int first = i;
            while (matchSeq[i] == 1 && i < ions.size())
            { // 记下多少个连续匹配的
                log++;
                i++;
            }
            double logs = log / 10.0;
            //            double logs = 0 ;
            while (first != i)
            {
                if (matchIndex[first])
                {
                    mutScore[first] += (z + logs);
                }
                else
                {
                    mutScore[first] += (p + logs);
                }
                first++;
            }
            if (i != ions.size())
            {
                i--;
            }
        }
    }

    //        cout << "完全匹配峰" << endl;
    //        for (i = 0; i < matchIndex.size(); ++i) {
    //            cout << " " << matchIndex[i];
    //        }
    //        cout << endl;
    //        cout << "匹配连续峰" << endl;
    //        for (i = 0; i < matchSeq.size(); ++i) {
    //            cout << " " << matchSeq[i];
    //        }
    //        cout << endl;
    //        cout << "系数" << endl;
    //        for (i = 0; i < mutScore.size(); ++i) {
    //            cout << " " << mutScore[i];
    //        }
    //        cout << endl;
    for (i = 0; i < ions.size(); ++i)
    {
        if (rm_dup[ions[i].thoe_id])
        {
            rm_dup[ions[i].thoe_id] = 0;
            score += (fabs(refe_score[ions[i].thoe_id]) + fabs(peer_score[ions[i].mono_id]) * mutScore[i]);
        }
        else
        {
            score += fabs(peer_score[ions[i].mono_id]) * mutScore[i];
        }
    }
    score *= punish_nub;
    //    cout<<"score = "<<score<<endl;
}

// void g::proc::Get_SCORE(vector<node> &ions, vector<double> &refe_score,
//                         vector<double> &peer_score, double &score,double base_score) {
//     if (ions.size() == 0 ) {
//         score = 0 ;
//         return ;
//     }
//     int i = 0;
//     double p = 0.0;    //-标识修饰离子个数
//     double z = 0.0;    //->标识完全对准个数
//     vector<double> mutScore(ions.size(),0);  //扩大系数
//     vector<int> matchIndex(ions.size(),0); //标识是否为完全匹配
//     vector<int> matchSeq(ions.size(),0);   //标识是否连续
//     map<int,double> rm_dup_t ;    //前一个为下标，后一个为分值系数
//     map<int,double> rm_dup_m ;    //前一个为下标，后一个为分值系数
//     //记下选取离子的位置下标
//     for(i = 0 ; i < ions.size() ; i++) {
//         //标识完全匹配的离子
//         if (fabs(ions[i].index) < 1.1) {
//             matchIndex[i]=1;
//             ++z ; //完全匹配峰的个数
//         } else {
//             ++p ; //
//         }
//         //标识唯一理论离子 , rm_dup_t,size 是离子的个数
//         if (!rm_dup_t.count(ions[i].thoe_id)) {
//             rm_dup_t[ions[i].thoe_id] = 0 ;
//         }
//         //标识唯一质谱离子
//         if (!rm_dup_m.count(ions[i].mono_id)) {
//             rm_dup_m[ions[i].mono_id] = 0 ;
//         }
//     }
//     //标识出连续匹配的离子matchSeq
//     for (i = 0; i < ions.size() - 1; i++) {
//         if (fabs(ions[i].thoe_id - ions[i+1].thoe_id ) == 1 ) {
//             matchSeq[i] = 1 ;
//             matchSeq[i+1] = 1 ;
//         }
//     }
//     //基础得分系数
//     p = 1 - (p / (ions.size())) + base_score;        //缩小修饰得分比例
//     z = 1 + (z / (ions.size())) + base_score;        //扩大完全对准得分比例
//     //    p = 1 - (p / (ions.size())) + 0.1;        //缩小修饰得分比例
//     //    z = 1 + (z / (ions.size())) + 0.1;        //扩大完全对准得分比例
//     //当只有完全匹配的离子时，判定是否存在连续匹配离子
//     if ( (1 - (p / (ions.size()))) == 0) {
//         int lo = 0 ;
//         for(i = 0; i < matchSeq.size();i++) {
//             if (matchSeq[i]==1) {
//                 lo = 1;
//                 break ;
//             }
//         }
//         if (!lo) {
//             z -= (z / (ions.size()));
//         }
//     }
//     //设置连续匹配系数mutScore
//     for (i = 0; i < ions.size(); i++ ) {
//         if (matchSeq[i] == 0) {
//             if (matchIndex[i]) {
//                 mutScore[i] += z ;
//             } else {
//                 mutScore[i] += p ;
//             }
//         } else if ( matchSeq[i] == 1) {
//             double log = 0 ;
//             int first = i ;
//             while (matchSeq[i] == 1 && i < ions.size()) {  //记下多少个连续匹配的
//                 log++;
//                 i++;
//             }
//             //            double logs = log/10.0 ;
//             double logs = base_score ;
//             while (first != i ) {
//                 if (matchIndex[first]) {
//                     mutScore[first] += (z+logs) ;
//                 } else {
//                     mutScore[first] += (p+logs) ;
//                 }
//                 first++;
//             }
//             if (i != ions.size()) {
//                 i--;
//             }
//         }
//     }
//
//     //    cout << "完全匹配峰" << endl;
//     //    for (i = 0; i < matchIndex.size(); ++i) {
//     //        cout << " " << matchIndex[i];
//     //    }
//     //    cout << endl;
//     //    cout << "匹配连续峰" << endl;
//     //    for (i = 0; i < matchSeq.size(); ++i) {
//     //        cout << " " << matchSeq[i];
//     //    }
//     //    cout << endl;
//     //    cout << "系数" << endl;
//     //    for (i = 0; i < mutScore.size(); ++i) {
//     //        cout << " " << mutScore[i];
//     //    }
//     //    cout << endl;
//     for (i = 0 ; i <ions.size(); ++i) {
//         score += (fabs(refe_score[ions[i].thoe_id]) + fabs(peer_score[ions[i].mono_id])*mutScore[i]) ;
//     }
//     //    cout<<"score = "<<score<<endl;
// }

// 选出质量位移偏移为0的离子和修饰离子。	ions保存偏差为0左右的离子，modf_ions保存修饰离子
void mylib::data_stream::Get_all_ions(vector<node> &ions, vector<node> &ModifyIons,
                                      const vector<pair<long, long>> &alignment,
                                      const vector<double> &reference, const vector<double> &peer,
                                      const vector<double> &index,
                                      const vector<double> &ppm)
{
    node loca;
    int i = 0;
    // 分别取出y离子质量偏移为0 ，和带修饰的离子。
    for (i = 0; i < alignment.size(); i++)
    {
        // cout<<"i = " <<i<<endl;
        if (ppm[i] < 15 && index[i] != INT_MAX)
        {
            // 该离子可能正确
            if (fabs(index[i]) < 0.02)
            {
                // 此处看是在修饰离子中偶然出现的，则过滤掉.
                if (ModifyIons.size() > 10)
                {
                    continue;
                }
                if (i - 1 >= 0 && alignment[i - 1].second != alignment[i].second && ppm[i - 1] < 15 &&
                    index[i - 1] != INT_MAX) // 选取同位素
                {
                    double sub = fabs(reference[alignment[i - 1].second] - reference[alignment[i].second]);
                    if (ions.size() > 0)
                    {
                        if (sub < 1.007276 && ions[ions.size() - 1].mono_id != alignment[i - 1].second)
                        {
                            loca.thoe_id = alignment[i - 1].first;
                            loca.mono_id = alignment[i - 1].second;
                            loca.index = index[i - 1];
                            loca.ppm = ppm[i - 1];
                            ions.push_back(loca);
                        }
                    }
                    else if (sub < 1.007276)
                    {
                        loca.thoe_id = alignment[i - 1].first;
                        loca.mono_id = alignment[i - 1].second;
                        loca.index = index[i - 1];
                        loca.ppm = ppm[i - 1];
                        ions.push_back(loca);
                    }
                }
                loca.thoe_id = alignment[i].first;
                loca.mono_id = alignment[i].second;
                loca.index = index[i];
                loca.ppm = ppm[i];
                ions.push_back(loca);
                if (ModifyIons.size() != 0) // 此时之前的修饰值过滤掉
                {
                    ModifyIons.clear();
                }
            }
            else if (fabs(index[i]) > 0.02 && fabs(index[i]) < 2) // 该离子可能正确，取决于上一条的信息
            {
                if (ions.size() > 0)
                { // 判断c_ions[-1]是否合法
                    if (alignment[i].first - 1 == ions[ions.size() - 1].thoe_id ||
                        ModifyIons.size() <= 1) // 1.若是该离子不在修饰堆中，2.若在修饰堆中,同位素峰
                    {
                        if (i - 1 >= 0 && alignment[i - 1].second != alignment[i].second && ppm[i - 1] < 15 &&
                            index[i - 1] != INT_MAX) // 选去同位素峰
                        {
                            double sub = fabs(reference[alignment[i - 1].second] - reference[alignment[i].second]);
                            if (ions.size() > 0)
                            {
                                if (sub < 1.007276 && ions[ions.size() - 1].mono_id != alignment[i - 1].second)
                                {
                                    loca.thoe_id = alignment[i - 1].first;
                                    loca.mono_id = alignment[i - 1].second;
                                    loca.index = index[i - 1];
                                    loca.ppm = ppm[i - 1];
                                    ions.push_back(loca);
                                }
                            }
                            else if (sub < 1.007276)
                            {
                                loca.thoe_id = alignment[i - 1].first;
                                loca.mono_id = alignment[i - 1].second;
                                loca.index = index[i - 1];
                                loca.ppm = ppm[i - 1];
                                ions.push_back(loca);
                            }
                        }
                        loca.thoe_id = alignment[i].first;
                        loca.mono_id = alignment[i].second;
                        loca.index = index[i];
                        loca.ppm = ppm[i];
                        ions.push_back(loca);
                        if (ModifyIons.size() != 0)
                        {
                            ModifyIons.clear();
                        }
                    }
                    else
                    {
                        continue;
                    }
                }
                else
                {
                    loca.thoe_id = alignment[i].first;
                    loca.mono_id = alignment[i].second;
                    loca.index = index[i];
                    loca.ppm = ppm[i];
                    ions.push_back(loca);
                    if (ModifyIons.size() != 0)
                    {
                        ModifyIons.clear();
                    }
                }
            }
            else
            {                        // 该值为修饰值
                if (ions.size() > 0) // 判断究竟有没有确定的离子。
                {
                    if (alignment[i].first != ions[ions.size() - 1].thoe_id) // 如果有确定的离子则往下取修饰值
                    {
                        loca.thoe_id = alignment[i].first;
                        loca.mono_id = alignment[i].second;
                        loca.index = index[i];
                        loca.ppm = ppm[i];
                        ModifyIons.push_back(loca);
                    }
                }
                else
                {
                    loca.thoe_id = alignment[i].first;
                    loca.mono_id = alignment[i].second;
                    loca.index = index[i];
                    loca.ppm = ppm[i];
                    ModifyIons.push_back(loca);
                }
            }
        }
    }
    // for -- end
}

// Ions是需要得到的离子ptms以NC端作为参照
void mylib::data_stream::GetAnotherIons(const vector<pair<long, long>> &alignment, vector<modi> &mod,
                                        const vector<node> &MergIons, vector<node> &Ions,
                                        const vector<double> &index, const vector<double> &ppm,
                                        int lenth)
{
    // for(int i = 0; i < mod.size(); ++i ) {
    // 	cout<<mod[i].first<<" "<<mod[i].second <<" " << mod[i].mod_mass<<endl;
    // }
    reverse(mod.begin(), mod.end());
    // cout<<"翻转后:\n"<<endl;
    for (int i = 0; i < mod.size(); ++i)
    {
        double t = mod[i].first;
        mod[i].first = mod[i].second;
        mod[i].second = t;

        mod[i].first = lenth - mod[i].first;
        mod[i].second = lenth - mod[i].second;
        // cout<<mod[i].first<<" "<<mod[i].second <<" " << mod[i].mod_mass<<endl;
    }
    int k = 0;
    int i = 0;
    // 第一个区间单独拿出来
    for (k = 0; k < MergIons.size(); ++k)
    {
        // 首先判断是否小于第一个区间的Second位置
        if (MergIons[k].thoe_id <= mod[0].second)
        {
            // 如果为0 则 OK ，添加进来，并且移动frist 位置 。
            if (fabs(MergIons[k].index) < 1)
            {
                Ions.push_back(MergIons[k]);
                mod[0].first = MergIons[k].thoe_id;
                // 如果相近则OK，并且退出该循环
            }
            else if (fabs(fabs(fabs(MergIons[k].index) - fabs(mod[0].mod_mass)) < 2))
            {
                Ions.push_back(MergIons[k]);
                mod[0].first = MergIons[k].thoe_id;
                break;
            }
        }
        else if (MergIons[k].thoe_id > mod[0].second)
        {
            break;
        }
    }
    for (i = k; i < MergIons.size(); ++i)
    {
        double mass = mod[0].mod_mass;
        for (int j = 1; j < mod.size(); ++j)
        {
            if (MergIons[i].thoe_id >= mod[j].first)
            {
                mass += mod[j].mod_mass;
            }
            if (MergIons[i].thoe_id > mod[j - 1].second && MergIons[i].thoe_id < mod[j].first)
            {
                if (fabs(fabs(MergIons[i].index) - fabs(mass)) < 4)
                {
                    Ions.push_back(MergIons[i]);
                }
                continue;
            }
            if (MergIons[i].thoe_id > mod[j].first && MergIons[i].thoe_id < mod[j].second)
            {
                if (fabs(fabs(MergIons[i].index) - fabs(mass)) < 4)
                {
                    Ions.push_back(MergIons[i]);
                    mod[j].first = MergIons[i].thoe_id;
                }
                continue;
            }
        }
        if (fabs(fabs(MergIons[i].index) - fabs(mass)) < 4)
        {
            Ions.push_back(MergIons[i]);
        }
    }
    // for(i = 0 ; i < Ions.size(); i++) {
    // 	cout<<"F = "<<Ions[i].thoe_id<<" S = "<<Ions[i].mono_id<<" Index = "<<Ions[i].index<<endl;
    // }
}

int mylib::data_stream::init_theory_ions_mass_(vector<double> &theory_ions_mass,
                                               int cut_n, int cut_c,
                                               double addMass) // 初始化当前蛋白质形式质量
{
    if (theory_ions_mass.empty())
    {
        return 0;
    }
    if (cut_n > theory_ions_mass.size() - cut_c - 1)
    {
        return 0;
    }
    if (cut_n > 0)
    { //-----------如果c端截断
        double subMass = theory_ions_mass[cut_n - 1] - addMass;
        theory_ions_mass.erase(theory_ions_mass.begin(), theory_ions_mass.begin() + cut_n); //-----y离子截去前面离子
        for (int i = 0; i < theory_ions_mass.size(); ++i)
        { //-----对于每一个y离子进行更新
            theory_ions_mass[i] -= subMass;
        }
    }
    if (cut_c > 0)
    {
        theory_ions_mass.erase(theory_ions_mass.end() - cut_c, theory_ions_mass.end()); //-----b离子截去末端
    }
    return 1;
    //----------------------------------------------------------------
}

bool mylib::data_stream::verifyModify(string &raw, string &Ptms)
{

    string cmpRaw;
    string cmpPtms;
    for (char ch : raw)
    {
        if (ch == 'D')
        {
            cmpRaw += "MM";
        }
        else
        {
            cmpRaw.push_back(ch);
        }
    }
    for (char ch : Ptms)
    {
        if (ch == 'D')
        {
            cmpPtms += "MM";
        }
        else
        {
            cmpPtms.push_back(ch);
        }
    }
    sort(cmpRaw.begin(), cmpRaw.end());
    sort(cmpPtms.begin(), cmpPtms.end());
    // cout<<"cmpRaw = "<<cmpRaw<< " cmpPtms = "<<cmpPtms<<endl;
    if (cmpRaw == cmpPtms)
    {
        return true;
    }
    return false;
}

bool mylib::data_stream::searchMatchIons(vector<pair<long, long>> &alignmentN, vector<pair<long, long>> &alignmentC,
                                         vector<double> &referenceN, vector<double> &referenceC, vector<double> &peer)
{
    for (int xi = 0; xi < alignmentN.size(); ++xi)
    {
        double ppm = fabs(referenceN[alignmentN[xi].first] - peer[alignmentN[xi].second]) /
                     referenceN[alignmentN[xi].first] * 1000000;
        if (ppm < 15)
        {
            return true;
        }
    }
    for (int xi = 0; xi < alignmentC.size(); ++xi)
    {
        double ppm = fabs(referenceC[alignmentC[xi].first] - peer[alignmentC[xi].second]) /
                     referenceC[alignmentC[xi].first] * 1000000;
        if (ppm < 15)
        {
            return true;
        }
    }
    return false;
}

void mylib::data_stream::filterAlignment(vector<pair<long, long>> &alignment, vector<double> &sub,
                                         vector<double> &index, vector<double> &ppm)
{
    vector<pair<long, long>> alignment_T;
    vector<double> sub_T;
    vector<double> index_T;
    vector<double> ppm_T;
    for (int xi = 0; xi < alignment.size(); xi++)
    {
        if (ppm[xi] < 15)
        {
            alignment_T.push_back(alignment[xi]);
            sub_T.push_back(sub[xi]);
            index_T.push_back(index[xi]);
            ppm_T.push_back(ppm[xi]);
        }
    }
    alignment.swap(alignment_T);
    sub.swap(sub_T);
    index.swap(index_T);
    ppm.swap(ppm_T);
}

int mylib::data_stream::txpb_binarey_search_ex(const vector<double> &pstArray, int iLength, double ullTime)
{
//    int middle = 0;
//    int low = 0, high = iLength - 1;
//    if (pstArray.empty() || iLength < 1)
//    {
//        printf("%s %d error !\n", __func__, __LINE__);
//        return -1;
//    }
//
//    if (iLength == 1)
//    {
//        return iLength - 1;
//    }
//
//    if (ullTime <= pstArray[0])
//    {
//        return 0;
//    }
//    else if (ullTime >= pstArray[iLength - 1])
//    {
//        return iLength - 1;
//    }
//
//    while (low <= high)
//    {
//        middle = (low + high) / 2;
//        if (abs(pstArray[middle + 1] - ullTime) > abs(pstArray[middle] - ullTime))
//        {
//            high = middle - 1;
//        }
//        else
//        {
//            low = middle + 1;
//        }
//    }
//    return (abs(pstArray[middle + 1] - ullTime) > abs(pstArray[middle] - ullTime)) ? middle : (middle + 1);



    ///by lyc
    if (pstArray.empty() || iLength < 1)
    {
        printf("%s %d error !\n", __func__, __LINE__);
        return -1;
    }

    if (iLength == 1)
    {
        return iLength - 1;
    }

    if (ullTime <= pstArray[0])
    {
        return 0;
    }
    else if (ullTime >= pstArray[iLength - 1])
    {
        return iLength - 1;
    }

    int left = 0;
    int right = iLength - 1;


    if(ullTime <= pstArray[left])
        return left;
    if(ullTime >= pstArray[right])
        return right;


    while(right - left > 1)        //while循环里面的r - l > 1 这个限定条件很重要，当r yu l 的差值为1 的时候就应该停止循环了，在分别ar[l] 与 ar[r] 两个值讨论看那个值更合适
    {
        int mid = (left + right) / 2;
        if(ullTime >= pstArray[mid])
            left = mid;	  //在这里 l 直接等于 mid ，否则会出错
        else
            right = mid;	  // 同理
    }
    long long int x = abs(pstArray[left] - ullTime);
    long long int y = abs(pstArray[right] - ullTime);
    return x <= y ? left : right;


}

long mylib::data_stream::get_fragment_ions_nub(const std::vector<node> &all_ions)
{
    map<int, int> remove;
    for (int xi = 0; xi < all_ions.size(); ++xi)
    {
        if (remove.count(all_ions[xi].thoe_id))
        {
            continue;
        }
        //        remove[all_ions[xi].thoe_id] = all_ions[xi].mono_id;
        remove.insert(make_pair(all_ions[xi].thoe_id, 1));
    }
    return remove.size();
}

int mylib::data_stream::get_modification_table(const string &mod_file_name, vector<string> &mod)
{
    ifstream in(mod_file_name, ios::in);
    ofstream out_massage("data/modification/mod_table_massage.txt");
    map<string, double> mod_mass_map;
    if (!in.good())
    {
        return 0;
    }
    string buff;
    while (in.good())
    {
        getline(in, buff);
        if (buff.size() <= 10 || buff[0] == '#')
        {
            continue;
        }
        mod.push_back(buff);
        string mod_name = buff.substr(0, buff.find_first_of(','));
        string next_ = buff.substr(buff.find_first_of(',') + 1);
        double mod_mass = atof(next_.substr(0, next_.find_first_of(',')).c_str());
        mod_mass_map[mod_name] = mod_mass;
    }
    // 输出日志信息
    int mod_nub = mod_mass_map.size();
    out_massage << "Mod number = " << mod_nub << endl;
    for (map<string, double>::iterator it = mod_mass_map.begin(); it != mod_mass_map.end(); it++)
    {
        out_massage << "mod name : " << setiosflags(ios::fixed) << setprecision(5) << it->first << " , mod mass = "
                    << it->second << endl;
    }
    time_t now_time = time(NULL);
    // 获取本地时间
    tm *t_tm = localtime(&now_time);
    // 转换为年月日星期时分秒结果，如图：
    out_massage << "local time is : " << asctime(t_tm) << endl;
    // 输出日志结束
    // 组合修饰
    int n = mod_nub * 5; // 25
    int k = 1;
    //    Solution s; //组合接口
    map<string, double> combine; // id_id_id_id , mass
    map<int, string> id_mod;
    map<int, double> id_mass;
    double max_mass = 0;
    char max_id;
    for (map<string, double>::iterator it = mod_mass_map.begin(); it != mod_mass_map.end(); it++)
    {
        id_mod[k] = it->first;
        id_mass[k++] = it->second;
        if (max_mass < it->second)
        {
            max_mass = it->second;
            max_id = k + '0';
        }
    }
    int maxSize = 500 / max_mass;
    string out_str; // 退出边界
    for (int i = 0; i < maxSize; ++i)
    {
        out_str += max_id;
    }
    bool c = true;
    map<string, double> con_copy;
    while (c)
    {
    }
    //
    return 1;
}

int mylib::data_stream::find_match_ions(const std::vector<double> &arr_mod, const std::vector<node> &arr_ions)
{
    int log = 0;
    for (int i = 0; i < arr_ions.size(); ++i)
    {
        int n = mylib::data_stream::txpb_binarey_search_ex(arr_mod, arr_mod.size(), arr_ions[i].index);
        if (fabs(arr_mod[n] - arr_ions[i].index) > 0.1)
        {
            continue;
        }
        log++;
    }
    return log;
}

int mylib::data_stream::find_error_match_ions(const std::vector<double> &arr_mod, const std::vector<node> &arr_ions)
{
    int count = 0;
    for (int pi = 0; pi < arr_ions.size(); ++pi)
    {
        int cc = mylib::data_stream::txpb_binarey_search_ex(arr_mod, arr_mod.size(), arr_ions[pi].index);
        if (arr_mod[cc] == 0 || fabs(arr_mod[cc] - arr_ions[pi].index) <= 0.1)
        {
            continue;
        }
        count++;
    }
    return count;
}

bool mylib::data_stream::if_is_error_ions(const std::vector<double> &arr_mod, const node &p)
{
    int cc = mylib::data_stream::txpb_binarey_search_ex(arr_mod, arr_mod.size(), p.index);
    if (arr_mod[cc] == 0 || fabs(arr_mod[cc] - p.index) <= 0.1)
    {
        return true;
    }
    return false;
}

double mylib::data_stream::analysis_all_ions(const std::vector<node> &arr_ions, const std::vector<node> &arr_ions_2,
                                             const std::vector<double> &mono,
                                             const std::vector<double> &arr_mod, double mod_mass)
{
    std::vector<node> all = arr_ions;
    int count_zero = 0;
    int count_error = 0;
    int count_match = 0;
    for (int i = 0; i < arr_ions_2.size(); ++i)
    {
        all.push_back(arr_ions_2[i]);
    }
    for (int i = 0; i < all.size(); ++i)
    {
        if (fabs(all[i].index) < 1.1)
        { // 首端匹配的离子
            count_zero++;
        }
        if (fabs(all[i].index - mod_mass) < 2)
        { // 末端匹配的离子
            count_match++;
        }
        if (!mylib::data_stream::if_is_error_ions(arr_mod, all[i]))
        { // 匹配不对齐的离子
            count_error++;
        }
    }
    double p = 0;
    double all_size = all.size();
    double mono_size = mono.size();
    if (all_size / mono_size * 100 > 10 && all_size / mono_size * 100 <= 15)
    {
        p = 1;
    }
    else if (all_size / mono_size * 100 > 15 && all_size / mono_size * 100 <= 30)
    {
        p = 2;
    }
    else if (all_size / mono_size * 100 > 30)
    {
        p = 3;
    }
    // 如果满足count_zero > 0 , count_match > 0 ; 最终得分翻倍。
    //    if (fabs(mod_mass)<=eps && count_zero > 0 && count_match > 0) {
    //        return 2 + p ;
    //    } else
    if (count_zero > 0 && count_match > 0)
    {
        return 2 + p;
    }
    else if (count_zero > 0 || count_match > 0)
    {
        return 1 + p;
    }
    return 0.5 + p;
}
// 5 标识识别离子占比高于30 , 并且有完全识别的离子和最后推测出的修饰离子
// 4 标识占比在 15 - 30 之间,并且有有完全识别的离子和最后推测出的修饰离子
// 3 标识识别离子占比在10 - 15 之间 ,并且有并且有有完全识别的离子和最后推测出的修饰离子
// 1.5
// 2.5
// 3.5

bool mylib::data_stream::_ppm_H(double t, double m)
{
//    double PPM_re = (t - m) / (t)*1000000;
//    double PPM_DeH_re = (t - 1.007276 - m) / (t - 1.007276) * 1000000;
//    double PPM_AeH_re = (t + 1.007276 - m) / (t + 1.007276) * 1000000;
//    if (fabs(PPM_re) < 15 || fabs(PPM_DeH_re) < 15 || fabs(PPM_AeH_re) < 15)
//    {
//        return true;
//    }
//    return false;


    ///by lyc ppm不计算多氢少氢
//    double PPM_re = (t - m) / (t);
//    PPM_re *= 1000000;
//     if (fabs(PPM_re) < 15)
//    {
//        return true;
//    }
//    return false;

    ///考虑同位素峰
    double PPM_re = (t - m) / (t)*1000000;
    double PPM_DeH_re = (t - 1.0024 - m) / (t - 1.0024) * 1000000;
    double PPM_AeH_re = (t + 1.0024 - m) / (t + 1.0024) * 1000000;
    if (fabs(PPM_re) < 15 || fabs(PPM_DeH_re) < 15 || fabs(PPM_AeH_re) < 15)
    {
        return true;
    }
    return false;

}

bool mylib::data_stream::_ppm_H_getMinPPM(double t, double m, double &minPPM)
{
    double PPM_re = (t - m) / (t)*1000000;
    double PPM_DeH_re = (t - 1.007276 - m) / (t - 1.007276) * 1000000;
    double PPM_AeH_re = (t + 1.007276 - m) / (t + 1.007276) * 1000000;
    minPPM = fabs(PPM_re) > fabs(PPM_DeH_re) ? PPM_DeH_re : PPM_re;
    minPPM = fabs(minPPM) > fabs(PPM_AeH_re) ? PPM_AeH_re : minPPM;
    if (fabs(PPM_re) < 15 || fabs(PPM_DeH_re) < 15 || fabs(PPM_AeH_re) < 15)
    {
        return true;
    }
    return false;
}

bool mylib::data_stream::_ppm_H_re(double t, double m)
{
    double PPM_re = (t - m) / (t)*1000000;
    double PPM_DeH_re = (t - 1.007276 - m) / (t - 1.007276) * 1000000;
    double PPM_AeH_re = (t + 1.007276 - m) / (t + 1.007276) * 1000000;
    if (PPM_re < 15 || PPM_DeH_re < 15 || PPM_AeH_re < 15)
    {
        return true;
    }
    return false;
}

bool mylib::data_stream::_ppm(double t, double m)
{
    double PPM_re = (t - m) / (t)*1000000;
    if (PPM_re <= 15)
    {
        return true;
    }
    return false;
}

/**
 * 用于计算互补离子系数得分
 * @param ionsC
 * @param ionsN
 * @param length
 * @return
 */
double mylib::data_stream::getHuBuIonsScore(vector<node> &ionsC, vector<node> &ionsN, int length)
{
    double score = 0;
    set<long> setC;
    for (int i = 0; i < ionsC.size(); ++i)
    {
        setC.insert(ionsC[i].thoe_id + 1);
    }
    int count = 0;
    for (int i = 0; i < ionsN.size(); ++i)
    {
        if (setC.count(length - (ionsN[i].thoe_id + 1)))
        {
            count++;
        }
    }
    return count;
}

double mylib::data_stream::getSeqIonsScore(vector<node> &ionsC, vector<node> &ionsN, int &count)
{
    double sumSeqIons = 1;
    double scoreC = 0;
    double scoreN = 0;
    for (int i = 1; i < ionsC.size(); ++i)
    {
        if (ionsC[i].thoe_id - ionsC[i - 1].thoe_id <= 1)
        {
            sumSeqIons++;
        }
    }
    for (int i = 1; i < ionsN.size(); ++i)
    {
        if (ionsN[i].thoe_id - ionsN[i - 1].thoe_id <= 1)
        {
            sumSeqIons++;
        }
    }
    count = sumSeqIons;
    double cof = sumSeqIons / (double)(ionsC.size() + ionsN.size());
    return sqrt(cof * 100);
}

/**
 * 计算选出离子的得分
 * @param ions_r
 * @param refe_score
 * @param mod_arr
 * @param peer_score
 * @param score
 * @param base_score
 * @param mod_mass
 */
void mylib::data_stream::Get_SCORE_Un(vector<node> &ions_r, vector<double> &refe_score,
                                      const vector<double> &mod_arr, vector<double> &peer_score,
                                      double &score, double base_score, double mod_mass)
{
    vector<node> ions;
    double count = 0;
    double punish_nub = 0;
    // 如果修饰质量太大，可能误差就大，故适当降低
    if (fabs(mod_mass) > 500.0 / 2.0)
    {
        for (int pi = 0; pi < ions_r.size(); ++pi)
        {
            if (!mylib::data_stream::if_is_error_ions(mod_arr, ions_r[pi]))
            {
                count++;
                continue;
            }
            ions.push_back(ions_r[pi]);
        }
        punish_nub = 1 - (count / ions.size());
    }
    else
    {
        ions = ions_r;
        punish_nub = 1;
    }
    // 如果没有匹配离子，则得分为 0
    if (ions.size() == 0)
    {
        score = 0;
        return;
    }
    int i = 0;
    double p = 0.0;                          // 标识修饰离子个数
    double z = 0.0;                          // 标识完全对准个数
    vector<double> mutScore(ions.size(), 0); // 扩大系数
    vector<int> matchIndex(ions.size(), 0);  // 标识是否为完全匹配
    vector<int> matchSeq(ions.size(), 0);    // 标识是否连续
    map<int, int> rm_dup;
    // 标识出完全匹配的离子和修饰离子matchIndex
    // 得到完全匹配离子的扩大系数和修饰离子的扩大系数
    for (i = 0; i < ions.size(); i++)
    {
        if (!rm_dup.count(ions[i].thoe_id))
        {
            rm_dup[ions[i].thoe_id] = 1;
        }
        if (fabs(ions[i].index) < 1.1)
        {
            matchIndex[i] = 1;
            ++z;
        }
        else
        {
            ++p;
        }
    }
    // 基础得分系数
    p = 1 - (p / (ions.size())) + base_score; // 缩小修饰得分比例
    z = 1 + (z / (ions.size())) + base_score; // 扩大完全对准得分比例
    // 标识出连续匹配的离子matchSeq
    for (i = 0; i < ions.size() - 1; i++)
    {
        if (fabs(ions[i].thoe_id - ions[i + 1].thoe_id) == 1)
        {
            matchSeq[i] = 1;
            matchSeq[i + 1] = 1;
        }
    }
    //    p = 1 - (p / (ions.size())) + 0.1;        //缩小修饰得分比例
    //    z = 1 + (z / (ions.size())) + 0.1;        //扩大完全对准得分比例
    // 当只有完全匹配的离子时，判定是否存在连续匹配离子
    if ((1 - (p / (ions.size()))) == 0)
    {
        int lo = 0;
        for (i = 0; i < matchSeq.size(); i++)
        {
            if (matchSeq[i] == 1)
            {
                lo = 1;
                break;
            }
        }
        if (!lo)
        {
            z -= (z / (ions.size()));
        }
    }
    // 设置连续匹配系数mutScore
    for (i = 0; i < ions.size(); i++)
    {
        if (matchSeq[i] == 0)
        {
            if (matchIndex[i])
            {
                mutScore[i] += z;
            }
            else
            {
                mutScore[i] += p;
            }
        }
        else if (matchSeq[i] == 1)
        {
            double log = 0;
            int first = i;
            while (matchSeq[i] == 1 && i < ions.size())
            { // 记下多少个连续匹配的
                log++;
                i++;
            }
            double logs = log / 10.0;
            //            double logs = 0 ;
            while (first != i)
            {
                if (matchIndex[first])
                {
                    mutScore[first] += (z + logs);
                }
                else
                {
                    mutScore[first] += (p + logs);
                }
                first++;
            }
            if (i != ions.size())
            {
                i--;
            }
        }
    }

    //            cout << "完全匹配峰" << endl;
    //            for (i = 0; i < matchIndex.size(); ++i) {
    //                cout << " " << matchIndex[i];
    //            }
    //            cout << endl;
    //            cout << "匹配连续峰" << endl;
    //            for (i = 0; i < matchSeq.size(); ++i) {
    //                cout << " " << matchSeq[i];
    //            }
    //            cout << endl;
    //            cout << "系数" << endl;
    //            for (i = 0; i < mutScore.size(); ++i) {
    //                cout << " " << mutScore[i];
    //            }
    //            cout << endl;

    // 22-10-12添加模块，如果存在修饰，那么该修饰有几个峰值10个峰值为1 ， 1个峰值为0.1 。 2个0.2
    double base = 0.1;
    vector<double> countBase; // 计算有几个相同的
    for (i = 0; i < ions.size(); ++i)
    {
        double hasSeqBase = 0.1;
        int j = i + 1;
        set<int> dup;
        dup.insert(ions[i].thoe_id);
        for (; j < ions.size(); ++j)
        {
            if (dup.count(ions[j].thoe_id))
            {
                continue;
            }
            if (fabs(ions[j].index - ions[i].index) > 1.5)
            {
                break;
            }
            hasSeqBase += base;
            dup.insert(ions[j].thoe_id);
        }
        while (i != j)
        {
            countBase.push_back(hasSeqBase);
            i++;
        }
        i = j - 1;
    }
    //    for (double c : countBase) {
    //        cout << " c = " << c << endl;
    //    }
    //
    //    cout << " ==== " << endl;
    for (i = 0; i < ions.size(); ++i)
    {
        if (rm_dup[ions[i].thoe_id])
        { // refer 和peer的分都算
            rm_dup[ions[i].thoe_id] = 0;
            score += (fabs(refe_score[ions[i].thoe_id]) + fabs(peer_score[ions[i].mono_id]) * mutScore[i]) * countBase[i];
        }
        else
        { // 只算peer的分
            score += fabs(peer_score[ions[i].mono_id]) * mutScore[i] * countBase[i];
        }
    }
    score *= punish_nub;
    //    cout<<"score = "<<score<<endl;
}

double mylib::data_stream::findClose(vector<double> arr, int n, double target)
{
    if (arr.size() == 0) {
        return 0 ; 
    }
    
    if (target <= arr[0])
    { // 特殊处理小于等于第一个数
        return arr[0];
    }
    if (target >= arr[n - 1])
    { // 特殊处理大于等于最后一个数
        return arr[n - 1];
    }

    int left = 0, right = n - 1;
    while (left <= right)
    { // 二分查找
        int mid = (left + right) / 2;
        if (arr[mid] == target)
        { // 直接找到目标值
            return arr[mid];
        }
        if (arr[mid] < target)
        { // 在右半部分查找
            left = mid + 1;
        }
        else
        { // 在左半部分查找
            right = mid - 1;
        }
    }

    // left指向左半部分的最后一个数，right指向右半部分的第一个数
    if (abs(arr[left] - target) < abs(arr[right] - target))
    {
        return arr[left]; // 返回最接近的值
    }
    else
    {
        return arr[right];
    }
}
