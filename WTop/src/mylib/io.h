#ifndef IO_H__
#define IO_H__

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
#include <map>
#include <fstream>

#include "data_steam.hpp"

using namespace std;

namespace mylib {
    namespace io {

        bool Readdata(const char *name, std::vector<double> &signals);

        bool ReadSignalSequence_int(const char *name, std::vector<int> &signals);

        bool WriteSignalSequence(const char *name, const std::vector<double> &signals);

        bool WriteSignalSequence_int(const char *name, const std::vector<int> &signals);

//---- write signal with name ----//
        bool WriteSignalSequence_withName(const char *name,
                                          const std::vector<double> &signals,
                                          const std::vector<std::string> &kmer_rec);

        bool WriteSignalSequence_int_withName(const char *name,
                                              const std::vector<int> &signals,
                                              const std::vector<std::string> &kmer_rec);

        void quickSort(std::vector<double> &arr, int begin, int end);

        void quickSort_ll(std::vector<double> &arr, int begin, int end, std::vector<std::vector<double> > &kk);

        void get_protein_name_seq_map(char *protein_file,
                                      std::map<std::string, std::string> &pro_name_seq_map);

        void init_var_ptm_combination(const string & variableFileName, string outPath, int maxSize, double maxMass);
    }

}

#endif
