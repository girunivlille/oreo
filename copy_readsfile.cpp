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
#include <cstdint>

using namespace std;

vector<string> split_string(const string& str, char separator) {
    //Separates the string with the separator and puts the resulting words in a vector of strings
    vector<string> result;
    istringstream iss(str);
    string word;

    while (getline(iss, word, separator)) {
        result.push_back(word);
    }

    return result;
}

void copy_the_file(const string& reads_file,const string& readsfile_copy){
    // string readsfile_copy="";
    vector<string> repo_files = split_string(reads_file, '/');
    // readsfile_copy = "copie-"+repo_files[repo_files.size()-1];
    fstream reads(reads_file, ios::in);
    fstream copy_file(readsfile_copy, ios::out);
    int compteur=0;
    string ajout;
    if(reads && copy_file){
        string line;
        while(!reads.eof()){
            getline(reads, line);
            if(line.length()>0 and line[0]=='>'){
                //la ligne est un id
                ajout = ">"+to_string(compteur);
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
    string readsfile_copy;
    if (argc!=3){
        cerr << "Mauvais nombre d'arguments" << endl;
        return EXIT_FAILURE;
    }
    else{
        readsfile = argv[1];
        readsfile_copy = argv[2];
    }

    copy_the_file(readsfile,readsfile_copy);
    return 0;
}

