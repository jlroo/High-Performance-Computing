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
#include <string.h>

using namespace std;


bool read_fix(const char *filename, vector<string> &data) {
    
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

// Fills lps[] for given patttern pat[0..M-1]
void computeLPSArray(const char *pat, size_t M, int *lps)
{
    // length of the previous longest prefix suffix
    int len = 0;
    
    lps[0] = 0; // lps[0] is always 0
    
    // the loop calculates lps[i] for i = 1 to M-1
    int i = 1;
    while (i < M)
    {
        if (pat[i] == pat[len])
        {
            len++;
            lps[i] = len;
            i++;
        }
        else // (pat[i] != pat[len])
        {
            // This is tricky. Consider the example.
            // AAACAAAA and i = 7. The idea is similar
            // to search step.
            if (len != 0)
            {
                len = lps[len-1];
                
                // Also, note that we do not increment
                // i here
            }
            else // if (len == 0)
            {
                lps[i] = 0;
                i++;
            }
        }
    }
};


// Prints occurrences of txt[] in pat[]
void KMPSearch(const char *pat, const char *txt, vector<int> &end_lines, size_t N)
{
    size_t M = strlen(pat);
    //size_t N = strlen(txt);
    
    // create lps[] that will hold the longest prefix suffix
    // values for pattern
    int lps[M];
    
    // Preprocess the pattern (calculate lps[] array)
    computeLPSArray(pat, M, lps);
    
    int i = 0; // index for txt[]
    int j = 0; // index for pat[]
    while (i < N)
    {
        if (pat[j] == txt[i])
        {
            j++;
            i++;
        }
        
        if (j == M)
        {
            //printf("Found pattern at index %d n", i-j);
            end_lines.push_back(1+i-j);
            j = lps[j-1];
        }
        
        // mismatch after j matches
        else if (i < N && pat[j] != txt[i])
        {
            // Do not match lps[0..lps[j-1]] characters,
            // they will match anyway
            if (j != 0)
                j = lps[j-1];
            else
                i = i+1;
        }
    }
};

// Read the file all at once.
char * read_buffer(const char *filename, size_t &length)
{
    char * buffer = NULL;
    std::ifstream is (filename, std::ifstream::binary);
    if (is)
    {
        // get length of file:
        is.seekg (0, is.end);
        length = is.tellg();
        is.seekg (0, is.beg);
        
        //char * buffer = new char [length];
        buffer = new char [length];
        
        std::cout << "Reading " << length << " characters... ";
        
        // read data as a single block:
        is.read (buffer,length);
        
        if (is)
        std::cout << "all characters read successfully." << std::endl;
        else
        std::cerr << "error: only " << is.gcount() << " could be read" << std::endl;
        is.close();
        
    }
    return buffer;
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
