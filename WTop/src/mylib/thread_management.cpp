#include "thread_management.hpp"

bool cmpMs(msalign a, msalign b)
{
    return a.getPrecursorMass() < b.getPrecursorMass();
}

// 多线程程序函数调用入口
void mylib::thread_management::thread_wtop_porcess_start(const string &spactrue_file_name,
                                                         const string &variablePTMsFileOutPath,
                                                         const string &ptmOutFileName, const string &ptmOutFileName2,
                                                         const vector<protein_ptr> &protein_container)
{
    prsm::wtop_process_share_ptr proc = std::make_shared<prsm::WTop_Process>();
    proc->wtop_process(spactrue_file_name,
                       variablePTMsFileOutPath,
                       ptmOutFileName, ptmOutFileName2, protein_container);
    proc = nullptr;
}

void mylib::thread_management::separate_msalign_file(string &msalign_file,
                                                     vector<string> &separate_msalign_file,
                                                     int separate_num)
{

    if (msalign_file.find("msalign") == std::string::npos)
    {
        EX_TRACE("spectrum file Name Error!");
        exit(1);
    }

    ifstream in(msalign_file, ios::in);
    string buff;
    vector<msalign> ms_container; // 保存该文件数据
    while (getline(in, buff))
    { // 读取所有质谱数据
        mylib::data_stream::rtrim_s(buff);
        if (buff == "BEGIN IONS")
        {
            msalign ms;
            ms.insert_ions_msg_container(buff);
            while (buff != "END IONS")
            {
                if (in.good())
                {
                    getline(in, buff);
                    mylib::data_stream::rtrim_s(buff);
                    if (buff.find("PRECURSOR_MASS") != std::string::npos)
                    {
                        string value_str = buff.substr(buff.find("=") + 1);
                        ms.setPrecursorMass(stod(value_str));
                    }
                }
                ms.insert_ions_msg_container(buff);
            }
            ms_container.push_back(ms);
        }
    }
    in.close();

    int mid_n = ms_container.size() / separate_num;

    sort(ms_container.begin(), ms_container.end(), cmpMs);
    vector<vector<msalign>> eachFileMs(separate_num);

    // 均分质谱数量
    int eachIndex = 0;
    for (int i = 0; i < ms_container.size(); ++i)
    {
        eachIndex = i % separate_num;
        eachFileMs[eachIndex].push_back(ms_container[i]);
    }

    // 得到每个质谱文件名
    string sFileName = "";
    for (int i = 1; i <= separate_num; ++i)
    {
        sFileName = msalign_file.substr(0, msalign_file.find("msalign") - 1) + "_" + to_string(i) + ".msalign";
        separate_msalign_file.push_back(sFileName);
    }

    // 输出到对应文件
    int total_number = 0;
    for (int i = 0; i < eachFileMs.size(); ++i)
    {
        //    cout << "ms number : " << eachFileMs[i].size() << " file path : " << separate_msalign_file[i] << endl;
        ofstream out(separate_msalign_file[i], ios::out);
        vector<msalign> for_ms_con = eachFileMs[i];
        total_number += for_ms_con.size();
        for (msalign for_ms : for_ms_con)
        {
            for (string f_str : for_ms.getIonsMsgContainer())
            {
                out << f_str << endl;
            }
            out << endl;
        }
        out.close();
    }

    cout << "ms number - " << ms_container.size() << endl;
    cout << "separate ms total number - " << total_number << endl;
}

void mylib::thread_management::init_multi_protein_thread_start(const string &separate_protein_file,
                                                               const int separate_container_index)
{
    map<string, string> protein_name_seq_map;
    char *temp_protein_file = (char *)malloc((separate_protein_file.size() + 1) * sizeof(char));
    strcpy(temp_protein_file, separate_protein_file.c_str());

    mylib::io::get_protein_name_seq_map(temp_protein_file, protein_name_seq_map);

    for (map<string, string>::const_iterator it = protein_name_seq_map.begin(); it != protein_name_seq_map.end(); ++it)
    {
        protein_ptr proteinPtr = new protein(it->first, it->second);
        separate_protein_ptr[separate_container_index].push_back(proteinPtr); // prsm :: protein_container ;

        //init protein processor 
        std::lock_guard<std::mutex> lock(mutex_);
        ++init_protein_process_num;
        if (init_protein_process_num != 0 && (int)init_protein_process_num % 1000 == 0)
        {
            cout << std::flush << "\r"
                 << "init protein - process " << setprecision(4)
                 << (double)init_protein_process_num / protein_total_number * 100 << " % ";
        }
    }

    free(temp_protein_file);
}

void mylib::thread_management::init_protein_process_rate_thread()
{
    while (init_protein_process_num != protein_total_number)
    {
        if (init_protein_process_num != 0 && (int)init_protein_process_num % 500 == 0)
        {
            cout << std::flush << "\r"
                 << "init protein - process " << setprecision(4)
                 << (double)init_protein_process_num / protein_total_number * 100 << " % ";
        }
    }
}

void mylib::thread_management::multi_thread_startor(string &protein_file, string &msalign_file,
                                                    string &var_ptm_combination_file,
                                                    int thread_num, vector<string> &prsm_out_file)
{
    vector<string> separate_msalign_files; // 多线程质谱分离文件
    vector<string> separate_protein_files;
    vector<string> separate_var_ptm_prsm_out_file; // 多线程已知ptms输出文件
    vector<string> separate_one_ptm_prsm_out_file; // 多线程未知ptms输出文件
    vector<thread> wtop_thread_container;          // wtop process all thread
    vector<thread> init_protein_thread_container;  // protein file process thread
    // separate msalign files and protein files
    mylib::thread_management::separate_msalign_file(msalign_file,
                                                    separate_msalign_files, thread_num); // 分层质谱文件，并得出每个文件的输出文件

    mylib::thread_management::separate_protein_file(protein_file,
                                                    separate_protein_files, thread_num);

    cout << "protein_database have " << protein_total_number << " proteins " << endl;

    // init each separate protein file to mylib::thread_management::separate_protein_ptr_container
    for (int i = 0; i < separate_protein_files.size(); ++i)
    {
        // if (i == 0)
        // {
        //     init_protein_thread_container.push_back(thread(init_protein_process_rate_thread));
        // }
        vector<protein_ptr> protein_ptr_container;
        separate_protein_ptr.push_back(protein_ptr_container);
        init_protein_thread_container.push_back(thread(mylib::thread_management::init_multi_protein_thread_start,
                                                       separate_protein_files[i], i));
    }
    // wait each init protein thread end
    for (int i = 0; i < init_protein_thread_container.size(); ++i)
    {
        init_protein_thread_container[i].join();
    }
    // clear the separate protein files
    for (int i = 0; i < separate_protein_files.size(); ++i)
    {
        remove(separate_protein_files[i].c_str());
    }
    // combination all protein_ptr to combination_protein_container
    vector<protein_ptr> com_protein_ptr_container;
    for (vector<protein_ptr> protein_ptr_container : separate_protein_ptr)
    {
        for (int i = 0; i < protein_ptr_container.size(); ++i)
        {
            com_protein_ptr_container.push_back(protein_ptr_container[i]);
        }
    }
    separate_protein_ptr.clear(); // free

    //test out protein seq max length and avg length
    //test by luo

//    long max_length = 0,sum_length = 0;
//    for(protein_ptr prp : com_protein_ptr_container){
//        sum_length += prp->getProteinSequence().size();
//        if(prp->getProteinSequence().size() > max_length)
//            max_length = prp->getProteinSequence().size();
//    }
//    cout<<endl<<"######################################"<<endl;
//    cout<<"max_protein_sequence_length:"<<max_length<<"   avg_protein_sequence_length:"<<sum_length/com_protein_ptr_container.size()<<endl;
//    cout<<"######################################"<<endl;
    //end test



    // every prsm thread - start
    for (int i = 0; i < thread_num; ++i)
    {
        string cur_thread_var_prsm_out_file =
            separate_msalign_files[i].substr(0, separate_msalign_files[i].find("msalign") - 1) + "_result" + ".txt"; // 每个线程输出的结果
        string cur_thread_one_prsm_out_file =
            separate_msalign_files[i].substr(0, separate_msalign_files[i].find("msalign") - 1) + "_unknown_result" + ".txt"; // 每个线程输出的结果

        wtop_thread_container.push_back(thread(mylib::thread_management::thread_wtop_porcess_start,
                                               separate_msalign_files[i],
                                               var_ptm_combination_file,
                                               cur_thread_var_prsm_out_file,
                                               cur_thread_one_prsm_out_file,
                                               com_protein_ptr_container));

        separate_var_ptm_prsm_out_file.push_back(cur_thread_var_prsm_out_file);
        separate_one_ptm_prsm_out_file.push_back(cur_thread_one_prsm_out_file);
    }
    // wait every thread - end
    for (int i = 0; i < thread_num; ++i)
    {
        wtop_thread_container[i].join();
    }
    // free protein memory
    for (int i = 0; i < com_protein_ptr_container.size(); ++i)
    {
        delete com_protein_ptr_container[i];
    }
    // merge var ptms prsm result file
    string outFileName = "";
    for (int i = 0; i < separate_var_ptm_prsm_out_file.size(); ++i)
    {
        outFileName = msalign_file.substr(0, msalign_file.find(".msalign")) + "_merge_result.txt";
        ifstream in(separate_var_ptm_prsm_out_file[i], ios::in);
        ofstream out(outFileName, ios::app);
        string buff;
        while (getline(in, buff))
        {
            out << buff << endl;
        }
        in.close();
        out.close();
    }
    prsm_out_file.push_back(outFileName);

    // merge unknown ptm prsm result file
    for (int i = 0; i < separate_one_ptm_prsm_out_file.size(); ++i)
    { // 合并未知修饰结果
        outFileName = msalign_file.substr(0, msalign_file.find(".msalign")) + "_unknown_merge_result.txt";
        ifstream in(separate_one_ptm_prsm_out_file[i], ios::in);
        ofstream out(outFileName, ios::app);
        string buff;
        while (getline(in, buff))
        {
            out << buff << endl;
        }
        in.close();
        out.close();
    }
    prsm_out_file.push_back(outFileName);

    // remove files
    for (int i = 0; i < separate_msalign_files.size(); ++i)
    {
        remove(separate_var_ptm_prsm_out_file[i].c_str());
        remove(separate_one_ptm_prsm_out_file[i].c_str());
    }
}

void mylib::thread_management::separate_protein_file(string &protein_file,
                                                     vector<string> &separate_protein_file_container,
                                                     int separate_num)
{
    // read all protein
    std::map<std::string, std::string> protein_name_seq_map;
    char *char_protein_file = (char *)malloc((protein_file.length() + 1) * sizeof(char));
    strcpy(char_protein_file, protein_file.c_str());
    mylib::io::get_protein_name_seq_map(char_protein_file, protein_name_seq_map);

    // separate protein to other file
    std::string protein_file_pre_ = protein_file.substr(0, protein_file.find("."));
    std::vector<std::vector<pair<string, string>>> separate_file_protein_container;
    for (int i = 0; i < separate_num; ++i)
    {
        vector<pair<string, string>> temp;
        separate_file_protein_container.push_back(temp);
    }
    int cur_index = 0;
    for (map<string, string>::iterator map_it = protein_name_seq_map.begin(); map_it != protein_name_seq_map.end(); ++map_it, ++cur_index)
    {
        separate_file_protein_container[cur_index % separate_num].push_back(make_pair(map_it->first, map_it->second));
    }

    // out to separate file
    int protein_num = 0;
    for (int i = 0; i < separate_num; ++i)
    {
        string separate_protein_file = protein_file_pre_ + "_" + to_string(i) + ".fasta";
        ofstream out(separate_protein_file, ios::out);
        protein_num += separate_file_protein_container[i].size();
        for (int j = 0; j < separate_file_protein_container[i].size(); ++j)
        {
            out << separate_file_protein_container[i][j].first << "\n"
                << separate_file_protein_container[i][j].second << endl;
        }
        out.close();

        separate_protein_file_container.push_back(separate_protein_file);
    }

    protein_total_number = protein_num;
    free(char_protein_file);
    //    cout << "raw protein num is " << protein_name_seq_map.size() << " separate protein all num is " << protein_num << endl ;
}
