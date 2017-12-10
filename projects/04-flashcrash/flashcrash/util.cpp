
/**
  util.cpp
  flash

  Created by Jose Luis Rodriguez on 11/15/17.
  Copyright © 2017 Jose Luis Rodriguez. All rights reserved.
  
  icpc -qopenmp -qopt-report util.cpp
  install_name_tool -change @rpath/libiomp5.dylib /opt/intel/compilers_and_libraries_2018.1.126/mac/compiler/lib/libiomp5.dylib a.out

**/

#include "util.hpp"
#include <iostream> 
#include <fstream>
#include <map>


using namespace std;

void help()
{
    fprintf(stderr,"\t--help | -h       : Print help message.\n");
    fprintf(stderr,"\t--start_tag | -b  : FixTag start, default (\x01).\n");
    fprintf(stderr,"\t--end_tag | -e    : FixTag end, default ('=').\n");
    fprintf(stderr,"\t--num_start | -n  : Index # for the start of the string.\n");
    fprintf(stderr,"\t--num_end | -m    : Index # for the end of the string.\n");
    fprintf(stderr,"\t--search | -r     : Search for an string in the fixdata.\n");
}


void tagSearch(const vector<string> &data,
               vector<string> &search,
               string fixtag,
               int num_start,
               int num_end,
               size_t nsize) {
#pragma omp parallel
    {
        vector<string> vec_private;
#pragma omp for nowait
        for(int i=0; i<nsize; i++) {
            string smatch = data[i].substr(data[i].find(fixtag.c_str())+num_start,num_end);
            vec_private.push_back(smatch);
        }
#pragma omp critical
        search.insert(search.end(), vec_private.begin(), vec_private.end());
    }
}


void dateVolume(vector<string> &data, map<string, int> &count_dates, size_t nsize){
#pragma omp parallel for default(shared)
    for (int i = 0; i < nsize; ++i)
    {
    #pragma omp atomic update
        count_dates[data[i]]+=1;
    }
}



int main (int argc, char* argv[])
{
    
    string path;
    string fixtag;
    int num_start = 0;
    int num_end = 0;
    const char * tag_search= NULL;
    const char * tag_start = NULL;
    const char * tag_end = NULL;
    
    //string path = "/Users/jlroo/cme/data/2010/XCME";
    //string  path = "/work/05191/jlroo/stampede2/data/01/XCME_MD_ES_20100104_20100108";
    //const char * tag_start = "\x01";
    //const char * tag_end = "=";
    //fixtag = "52";
    //int num_start = 4;
    //int num_end = 8;
    
    vector<string> data;
    vector<string> search;
    map<string, int> volume;
    int week_volume = 0;
    
    for (int i = 1; i < argc; ++i)
    {
#define check_index(i,str) \
if ((i) >= argc) \
{ fprintf(stderr,"Missing 2nd argument for %s\n", str); return 1; }
        
        
        fprintf(stderr,"\t--help | -h       : Print help message.\n");
        fprintf(stderr,"\t--start_tag | -b  : FixTag start, default (\x01).\n");
        fprintf(stderr,"\t--end_tag | -e    : FixTag end, default ('=').\n");
        fprintf(stderr,"\t--num_start | -n  : Index # for the start of the string.\n");
        fprintf(stderr,"\t--num_end | -m    : Index # for the end of the string.\n");
        fprintf(stderr,"\t--search | -r     : Search for an string in the fixdata.\n");
        
        if ( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0)
        {
            help();
            return 1;
        }
        else if (strcmp(argv[i],"--path") == 0 || strcmp(argv[i],"-p") == 0)
        {
            path = argv[i];
        }
        else if (strcmp(argv[i],"--tag") == 0 || strcmp(argv[i],"-t") == 0)
        {
            fixtag = argv[i];
        }
        else if (strcmp(argv[i],"--tag_start") == 0 || strcmp(argv[i],"-s") == 0)
        {
            tag_start = argv[i];
        }
        else if (strcmp(argv[i],"--tag_end") == 0 || strcmp(argv[i],"-d") == 0)
        {
            tag_end = argv[i];
        }
        else if (strcmp(argv[i],"--num_start") == 0 || strcmp(argv[i],"-n") == 0)
        {
            num_start = atoi( argv[i] );
        }
        else if (strcmp(argv[i],"--num_end") == 0 || strcmp(argv[i],"-m") == 0)
        {
            num_end = atoi( argv[i] );
        }
        else if (strcmp(argv[i],"--search") == 0 || strcmp(argv[i],"-r") == 0)
        {
            tag_search = argv[i];
        }
        else
        {
            fprintf(stderr,"Unknown option %s\n", argv[i]);
            help();
            return 1;
        }
    }
    
    if (tag_start == NULL || tag_end == NULL )  {
        const char * tag_start = "\x01";
        const char * tag_end = "=";
        fixtag = tag_start + fixtag + tag_end;
    }else if (tag_search != NULL){
        fixtag = tag_search;
    }
    else{
        fixtag = tag_start + fixtag + tag_end;
    }

    read_fix(path.c_str(), data);
    size_t n = data.size();
    tagSearch(data, search, fixtag,num_start,num_end, n);

    dateVolume(search, volume, n);
    
    for (auto& iter : volume)  {
        std::cout << "volume[" << iter.first << "] = " << iter.second <<std::endl;
        week_volume +=iter.second;
    }
    std::cout << "Total Volume = " << week_volume <<std::endl;
    
    return 0;
}


