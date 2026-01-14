#include "msalign_processor.hpp"

using namespace std;
namespace prsm {

    void Msalign::initCompIonsIndexMap()
    {
        set<int> s ; //用于去重复
        for (int i = 0 ; i < ions_mass_container.size(); ++i) {
            double t = precursor_mass_double - ions_mass_container[i];
            int index = mylib::data_stream::txpb_binarey_search_ex(ions_mass_container, ions_mass_container.size(), t);
            if ( index>=0 && !s.count(index) && fabs(t - ions_mass_container[index]) < 2)
            {
                compIonsIndexMap.insert(make_pair(i,index));
                s.insert(i);
            }
        }
    }

    void Msalign::Init_Msalign_Data(ifstream &ms_file) {
        string buff;
        if (!ms_file.good()) {
            cout << "Prsm End!" << endl;
            exit(1);
        }
        while (getline(ms_file,buff)) {
            mylib::data_stream::rtrim_s(buff);
            if (buff == "END IONS") {
                break;
            }
            if (buff.find("=") != std::string::npos) {
                string title = buff.substr(0,buff.find("="));
                string value = buff.substr(buff.find("=")+1);
                ms_msg_map.insert(make_pair(title,value));
            } else {
                vector<string> s_str;
                utils_string_util::Stringsplit(buff,'\t',s_str);
                if (s_str.size() == 3) {
                    ions_mass_container.push_back(stod(s_str[0]));
                    keyMsValueAbRadio.insert(make_pair(stod(s_str[0]), stod(s_str[1])));
                }
            }
        }
        Id = "ID=" + ms_msg_map["ID"];
        scans="SCANS=" + ms_msg_map["SCANS"];
        retention_time = "RETENTION_TIME=" + ms_msg_map["RETENTION_TIME"];
        activation = "ACTIVATION=" + ms_msg_map["ACTIVATION"];
        activation_sub = ms_msg_map["ACTIVATION"];
        ms_one_id = "MS_ONE_ID=" + ms_msg_map["MS_ONE_ID"];
        ms_one_scans = "MS_ONE_SCAN=" + ms_msg_map["MS_ONE_SCAN"];
        precursor_mz = "PRECURSOR_MZ=" + ms_msg_map["PRECURSOR_MZ"];
        precursor_charge = "PRECURSOR_CHARGE=" + ms_msg_map["PRECURSOR_CHARGE"];
        precursor_mass = "PRECURSOR_MASS=" + ms_msg_map["PRECURSOR_MASS"];
        precursor_mass_double = stod(ms_msg_map["PRECURSOR_MASS"]);
        precursor_intensity = "PRECURSOR_INTENSITY=" + ms_msg_map["PRECURSOR_INTENSITY"];

        sort(ions_mass_container.begin(), ions_mass_container.end());
        // g::proc::MergeSort(mono_mass,0,mono_mass.size());
    }


    void Msalign::findComplementaryIons() {
//        map<int, int> mp;
//        // vector<double> peer1 = this->mono_mass ;
//        Pr_mass -= 18.01056;
//        // cout<<"prmass = "<<Pr_mass<<endl;
//        for (int xi = 0; xi < this->mono_mass.size(); ++xi) {
//            if (mp.count(xi)) {
//                continue;
//            }
//            int n = d::deal::txpb_binarey_search_ex(this->mono_mass, this->mono_mass.size(),
//                                                    Pr_mass - this->mono_mass[xi]);
//            if (fabs((this->mono_mass[n] + this->mono_mass[xi] - Pr_mass) / (this->mono_mass[n] + this->mono_mass[xi]) *
//                     1000000) < 15) {
//                this->cpmlementaryIons.push_back(make_pair(xi, n));
//                // cout<<"xi ="<<mono_mass[xi]<<" n = "<<mono_mass[n]<<endl;
//                // cout<<"xi+n = "<<mono_mass[xi] + mono_mass[n]<<endl;
//                cmplNumb++;
//                mp[xi] = 1;
//                mp[n] = 1;
//            }
//        }
    }

}