#ifndef OPTS_H__
#define OPTS_H__

#include <iostream>
#include <sstream>
#include <cassert>
extern "C" {
#include <getopt.h>
}
#include "../util/exception.h"

struct options {
	//-> required parameter
	char input[65532];
	char peer[65532];
	char output[65532];
    char testoutb[65532]; //x
    char testouty[65532];
	//-> key parameter
	int radius;
	int level;
	float scale0;
	//-> vice parameter
	int verbose;
	int test = -1;

	char style[65532];

	int mode;
};

inline int GetOpts(int argc, char **argv, options* opts_) {

    static struct option longopts[] = {
	{ "help",            no_argument,            NULL,              'h' },
	{ "input",     required_argument,      NULL,    			  'i' },
	{ "peer",    	     required_argument,      NULL,      	'p' },
	{ "output",      	required_argument,      NULL,              'o' },
    { "testoutb",      	required_argument,      NULL,              'c' },
    { "testouty",      	required_argument,      NULL,              'x' },
	{ "radius",      required_argument,      NULL,              'r' },
	{ "level",      required_argument,      NULL,              'l' },
	{ "scale",      required_argument,      NULL,              's' },
	{ "verbose",         no_argument,            NULL,              'v' },
	{ "test",            no_argument,            NULL,              't' },
	{ "mode",            no_argument,            NULL,              'm' },
    { "get",            no_argument,            NULL,              'g' },
    { NULL,              0,                      NULL,               0 }
    };

    if((argc !=5 && argc != 7 && argc != 9 && argc != 11
    && argc != 13 && argc != 15 && argc != 17 && argc != 19 )
    && argc >= 3 || (argc == 2 && argv[1][0] != '-' && argv[1][1] != 'h') || argc == 1) {

	EX_TRACE("\n");
	EX_TRACE("_  ______  ____  ___ ___ ____ \n");
	EX_TRACE(" _  |  -___----|-- |   _   _ \n");
	EX_TRACE("____|__ | __| __| __| __|   | \n");
	EX_TRACE("____|___|___|___|___|___|___|\n");
	EX_TRACE("\n");
	EX_TRACE("version v1.0 (2020-09-01) \n");
	EX_TRACE("-------------------------------------------------------------\n");
	EX_TRACE("Arg Option :\n");
    EX_TRACE("[-i Fasta File Name] [-p Mass spectrum File Name] (Must)\n");
    EX_TRACE("[-t Number of Threads][-o Variable PTMs File Name] (Optional)\n");
	EX_TRACE("-------------------------------------------------------------\n");
	EX_TRACE("        for more detailed description, type '-h'  \n");
        return -1;
    }

    int ch;
    while((ch = getopt_long(argc, argv, "hi:p:o:c:x:r:l:s:v:t:m:g", longopts, NULL))!= -1) {
        switch(ch) {
        case '?':
        {
            EX_TRACE("Invalid option '%s'.", argv[optind-1]);
            return -1;
        }

        case ':':
        {
            EX_TRACE("Missing option argument for '%s'.", argv[optind-1]);
            return -1;
        }

        case 'h':
        {
            EX_TRACE("----------- WTop ---------- \n"
                     "version v1.0 (2020-9-1) \n"
                     "-------------------------------------------------------------\n"
                     "required:\n"
                     "[-i Fasta File Name] [-p Mass spectrum File Name]\n"
                     "optional:\n"
                     "[-t Number of Threads][-o Variable PTMs File Name]"
                     "-------------------------------------------------------------\n"
                     "**** required: ******\n"
                     "Fasta File Name: Human Protein sequence ...(.fasta);\n"
                     "Mass spectrum File Name:  Mass spectrum file ... (.msalign);\n"
                     "OUTPUT: Do not specified,  It's automatically generated according to the mass spectrum file;  \n"
                     "**** other parameters: ******\n"
                     "Number of Threads :  The number of threads (default 1)\n"
                     "Variable PTMs file Name : The variable ptms files (default null);\n"
                    );
            return -1;
        }

        case 'i':
        {
            std::istringstream iss(optarg);
            iss >> opts_->input;
            if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;

        case 'p':
        {
            std::istringstream iss(optarg);
            iss >> opts_->peer;
            if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;

        case 'o':
        {
            std::istringstream iss(optarg);
            iss >> opts_->output;
            if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;
       case 'c':
        {
            std::istringstream iss(optarg);
            iss >> opts_->testoutb;
            if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
       case 'x':
        {
            std::istringstream iss(optarg);
            iss >> opts_->style;
            if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;

        case 't':
        {
             std::istringstream iss(optarg);
             iss >> opts_->test;
             if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
             }
        }
        break;

        case 'g':
        {
            printf("444\n");
            std::istringstream iss(optarg);
            iss >> opts_->style;
            printf("555\n");
            if(iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;

        case 0:
            break;

        default:
            assert(false);
        }
    }
    return 1;
}

#endif

