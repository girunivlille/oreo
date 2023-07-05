#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <unistd.h>
#include <climits>
#include <algorithm>
#include <map>
#include <iomanip>
#include <list>
#include <stack>

using namespace std;

void copy_the_file(const string& reads_file){
    string readsfile_copy;
    readsfile_copy = "copie-"+reads_file;
    fstream reads(reads_file, ios::in);
    fstream copy_file(readsfile_copy, ios::out);
    int compteur=0;
    string ajout;
    if(reads && copy_file){
        string line;
        uint32_t total_count = 0;
        while(!reads.eof()){
            getline(reads, line);
            if(line.length()>0 and line[0]=='>'){
                //la ligne est un id 
                ajout = ">."+to_string(compteur);
                copy_file << ajout << endl;
                compteur++;
            }
            else if(line.length()>0){
                copy_file << line << endl;
            }
        }
    }
    else{
        if(!reads){
            cerr << "Error opening the reads file "+reads_file+" (copy_the_file)." << endl;
        }
        if(!copy_file){
            cerr << "Error opening the copy file "+readsfile_copy+" (copy_the_file)." << endl;
        }
    }
}

int main(int argc, char *argv[])
{
    string readsfile;
    if (argc!=2){
        cerr << "Mauvais nombre d'arguments" << endl;
        return EXIT_FAILURE;
    }
    else{
        readsfile = argv[1];
    }

    copy_the_file(readsfile);
    return 0;
}