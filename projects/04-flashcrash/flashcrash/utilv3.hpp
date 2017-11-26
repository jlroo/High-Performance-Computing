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

/*
bool openFix(string file_name, vector<string> & fixdata) {
    ifstream fixfile;
    fixfile.open(file_name);
    string line;
    
    if(!fixfile.is_open()) {
        cout << "The read file could not be opened. Check the location.\t";
        return 0;
    }
    
    while(fixfile.good()) {
        getline(fixfile, line);
        char * message = new char [line.size()];
        fixdata.push_back(message);
    }
    fixfile.clear();
    fixfile.close();
    return 1;
};
*/


#endif /* util_hpp */
