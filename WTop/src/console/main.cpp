#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>
#include <string>
#include <time.h>
#include <thread>

#include <fstream>
#include <cstring>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "opts.h"

#include "../config/init_config.hpp"
#include "../thread/manage_thread.hpp"
#include "../util/xmlUtil.h"


int main(int argc, char *argv[]) {
    clock_t start_time = clock();
    time_t t2 = time(0);
    char startTime[32], endTime[32];
    struct options opts;

    strftime(startTime, sizeof(startTime), "%Y-%m-%d %H:%M:%S", localtime(&t2));

    if (GetOpts(argc, argv, &opts) < 0) {
        EX_TRACE("**WRONG INPUT!**\n");
        return -1;
    }

    // init config_argument
    init_config_Ptr init_config_ptr = new init_config();
    bool input_is_right = init_config_ptr->read_config_ptr->get_input_config_arg(opts);

    if (!input_is_right) {
        EX_TRACE("**WRONG INPUT!**\n");
        return -1;
    }

    // ptm combination
    init_config_ptr->init_combination_ptm();

    // ptm search
    std::vector<string> prsm_out_file; // 记录输出文件(txt -> csv)

    manage_thread_Ptr manageThreadPtr = new manage_thread(init_config_ptr);
    manageThreadPtr->start_prsm_threads(prsm_out_file);
    delete manageThreadPtr;

    time_t t = time(0);
    strftime(endTime, sizeof(endTime), "%Y-%m-%d %H:%M:%S", localtime(&t));
    string startTimeStr(startTime);
    string endTimeStr(endTime);
    clock_t end_time = clock();


    cout << "Start Text File To CSV File " << endl ;
//    utils_string_util::txtFileToCsvFile(prsm_out_file, init_config_ptr->read_config_ptr->get_arg_map());  //txt -> csv
    outXmlForWTopMLS();
    cout << "End Text File To CSV File " << endl ;

    cout << "Running time is: "
         << static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC
         << " s " << endl;

    delete init_config_ptr;
    return 0;
}
