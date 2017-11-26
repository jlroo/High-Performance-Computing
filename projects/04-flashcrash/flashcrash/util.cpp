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
#include <iostream> // std::cout
#include <fstream>  // std::ifstream
#include <unordered_set>

using namespace std;

void help()
{
    fprintf(stderr,"nbody3 --help|-h --nparticles|-n --nsteps|-s --stepsize|-t\n");
}


/**
 map with private histogram
 critical
 update share map with values of private map
 
 **/

void tagSearch(const vector<string> &fixdata,
                         string fixtag,
                         int length,
                         vector<string> &search) {
    
    fixtag = "\x01" + fixtag + "=";
    char *buf = new char[fixtag.length() + 1];
    strcpy(buf, fixtag.c_str());
    
#pragma omp parallel
    {
        vector<string> vec_private;
#pragma omp for nowait //fill vec_private in parallel
        for(int i=0; i<length; i++) {
            //char buf[] = fixtag;
            string smatch = fixdata[i].substr(fixdata[i].find(buf)+4,8);
            vec_private.push_back(smatch);
        }
#pragma omp critical
        search.insert(search.end(), vec_private.begin(), vec_private.end());
    }
};


void Counter(vector<string> data, int *count){
    int nsize = (int) data.size();
    count = new int [nsize];

#pragma omp parallel for default(shared)
    for (int i = 0; i < nsize; ++i)
    {
    #pragma omp atomic update
    count[ stoi(data[i])] ++;
    }
}


// Read the file all at once.
int read_all (const char *filename)
{
    std::ifstream is (filename, std::ifstream::binary);
    if (is)
    {
        // get length of file:
        is.seekg (0, is.end);
        const int length = is.tellg();
        is.seekg (0, is.beg);
        
        char * buffer = new char [length];
        
        std::cout << "Reading " << length << " characters... ";
        // read data as a single block:
        is.read (buffer,length);
        
        if (is)
            std::cout << "all characters read successfully." << std::endl;
        else
            std::cerr << "error: only " << is.gcount() << " could be read" << std::endl;
        is.close();
        
        // ...buffer contains the entire file...
        cout<< buffer;
        
        delete[] buffer;
    }
    
    return 0;
}


int read_batches (const char *filename)
{
    std::ifstream is (filename, std::ifstream::binary);
    if (not(is))
    {
        std::cerr << "Error opening file " << filename << std::endl;
        return 1;
    }
    
    std::string str;
    while (is.good())
    {
        const size_t chunkSize = 1024;
        char buffer[chunkSize];
        
        std::cout << "Reading chunk: ";
        // read data as a single block:
        size_t nread = is.read (buffer, chunkSize).gcount();
        
        str.append( buffer, nread );
        
        if (is)
            std::cout << nread << " characters read successfully." << std::endl;
        else
            std::cerr << " only " << is.gcount() << " could be read" << std::endl;
    }
    
    if (is.eof())
        std::cout << "File sucessful read with " << str.length() << std::endl;
    
    is.close();
    
    return 0;
}

/*

int main (int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "filename not specified\n";
        return 1;
    }
    
    read_all( argv[1] );
    read_batches( argv[1] );
    
    return 0;
}
*/

int main (int argc, char* argv[])
{
    
    /*

    if (argc < 2)
    {
        std::cerr << "filename not specified\n";
        return 1;
    }
    
    */
    
    //read_all( argv[1] );
    //read_batches( argv[1] );
    
    //fprintf(stderr,"Number Objects = %d\n", n);
    //fprintf(stderr,"Path to file   = %s\n", in_file.c_str());
    string in_file = "/Users/jlroo/cme/data/2010/raw/XCME";
    read_all( in_file.c_str() );
    //fixdata = openFix(in_file,fixdata,length);
    string fixtag = "52";
    //dates = tagSearch(fixdata,fixtag, dates);

    
    return 0;
}

