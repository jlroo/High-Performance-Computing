//
//  util.hpp
//  flash
//
//  Created by Jose Luis Rodriguez on 11/15/17.
//  Copyright © 2017 Jose Luis Rodriguez. All rights reserved.
//

#ifndef util_hpp
#define util_hpp

#include <stdio.h>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;


bool read_fix(string filename, vector<string> &data) {
    
    string line;
    ifstream fixfile;
    fixfile.open(filename);
    
    if(!fixfile.is_open()) {
        cout << "The read file could not be opened. Check the location.\t";
        return 0;
    }
    
    int i=0;
    while(std::getline(fixfile, line))
    {
        data.insert(data.begin() + i, line);
        ++i;
    }
    
    fixfile.close();
    return 0;
};


// Read the file all at once.
int read_buffer(const char *filename, vector<char> &buff)
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
        //cout<< buffer;
        std::copy(buffer, buffer + length, std::back_inserter(buff));
        
        delete[] buffer;
    }
    return 0;
};


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
};



#endif /* util_hpp */
