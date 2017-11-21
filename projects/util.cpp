//
//  util.cpp
//  flash
//
//  Created by Jose Luis Rodriguez on 11/15/17.
//  Copyright © 2017 Jose Luis Rodriguez. All rights reserved.
//

// icpc -qopenmp -qopt-report util.cpp
// install_name_tool -change @rpath/libiomp5.dylib /opt/intel/compilers_and_libraries_2018.1.126/mac/compiler/lib/libiomp5.dylib a.out

#include "util.hpp"
#include <iostream>
#include <fstream>
#include <regex>

using namespace std;

void help()
{
    fprintf(stderr,"nbody3 --help|-h --nparticles|-n --nsteps|-s --stepsize|-t\n");
}

vector<string> openFix(string file_name, vector<string> fixdata, int length) {
    ifstream fixfile;
    fixfile.open(file_name);
    
    if(!fixfile.is_open()) {
        cout << "The read file could not be opened. Check the location.\t";
        return fixdata;
    }
    
    //cnt = 100000;
    for ( int i = 0; i < length; ++i )
        getline(fixfile, fixdata[i]);
    fixfile.close();
    
    return fixdata;
};



vector<string> tagSearch(vector<string> fixdata,
                         string fixtag,
                         int length,
                         vector<string> search) {
    
    fixtag = "\x01" + fixtag + "=";
    char *buf = new char[fixtag.length() + 1];
    strcpy(buf, fixtag.c_str());
    
#pragma omp parallel
    {
        vector<string> vec_private;
#pragma omp for nowait //fill vec_private in parallel
        for(int i=0; i<length; i++) {
            //char buf[] = fixtag;
            vec_private.push_back(fixdata[i].substr(fixdata[i].find(buf)+4,8));
            //vec_private.push_back(fixdata[i]);
        }
#pragma omp critical
        search.insert(search.end(), vec_private.begin(), vec_private.end());
        
    }
    
    return search;
};

int main (int argc, char* argv[])
{
    string in_file;
    
    for (int i = 1; i < argc; ++i)
    {
#define check_index(i,str) \
if ((i) >= argc) \
{ fprintf(stderr,"Missing 2nd argument for %s\n", str); return 1; }
        
        if ( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0)
        {
            help();
            return 1;
        }
        else if (strcmp(argv[i],"--file") == 0 || strcmp(argv[i],"-f") == 0)
        {
            check_index(i+1,"--file|-f");
            i++;
            if (typeid(*argv[i]) == typeid(string))
                in_file = string(argv[i]) ;
        }
        else
        {
            fprintf(stderr,"Unknown option %s\n", argv[i]);
            help();
            return 1;
        }
    }
    
    int length = 100000;
    vector<string> fixdata(length);
    vector<string> dates;
    
    //fprintf(stderr,"Number Objects = %d\n", n);
    in_file = "/Users/jlroo/cme/data/2010/raw/XCME";
    fprintf(stderr,"Path to file   = %s\n", in_file.c_str());
    
    fixdata = openFix(in_file,fixdata,length);
    string fixtag = "52";
    dates = tagSearch(fixdata,fixtag, length, dates);

    /*
    for(vector<string>::iterator it = dates.begin(); it != dates.end(); ++it) {
        cout << *it << endl;
    }
    */
    
    return 0;
}
