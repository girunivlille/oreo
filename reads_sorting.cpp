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
#include <queue>
#include <iomanip>
#include <stack>
#include "unordered_dense.h"
#include "zstr/src/zstr.hpp"

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

vector<uint64_t> read_line_position(const string& filename){
    //Notes in a vector the position of each begining of line in a fasta file
    vector<uint64_t> read_line_pos={0};
    fstream file(filename, ios::in);
    //zstr::ifstream file(filename, ios::in);
    if(file){
        string line;
        uint64_t total_count = 0;
        while(!file.eof()){
            getline(file, line);
            total_count += line.length();
            read_line_pos.push_back(total_count+1);
        }
        return read_line_pos;
    }
    else{
        cerr << "Error opening the file (read_line_position)." << endl;
    }
    return {};
}

vector<uint64_t> contig_line_position(const string& gfa_file){
    //Notes in a vector the position of the begining of each contig in a gfa file
    vector<uint64_t> contig_line_pos={0};
    fstream gfa(gfa_file, ios::in);
    if(gfa){
        string line;
        uint64_t total_count = 0;
        getline(gfa, line);
        total_count += line.length()+1;
        while(!gfa.eof()){
            getline(gfa, line);
            vector<string> words_of_line = split_string(line, '\t');
            if(words_of_line.size()>0 && words_of_line[0]=="S"){
                contig_line_pos.push_back(total_count);
            }
            total_count += line.length()+1;
        }
        return contig_line_pos;
    }
    else{
        cerr << "Error opening the file (read_line_position)." << endl;
    }
    return {};
}

string get_read_sequence(const vector<uint64_t>& read_line_pos, const string& filename, uint64_t read_id){
    //Returns in a string the 'read_id'th sequence of the fasta file
    auto pos_seq = read_line_pos[(read_id*2)+1];
    auto pos_entete_suivant = read_line_pos[(read_id*2)+2];

    //zstr::ifstream file_in(filename, ios::in);
    fstream file_in(filename, ios::in);
    char line_seq[100000];
    if(file_in){
        file_in.seekg((pos_seq+read_id*2), file_in.beg);
        file_in.read(line_seq, ((pos_entete_suivant+read_id)-(pos_seq+read_id)));
        line_seq[pos_entete_suivant-pos_seq+2] = 0;
        file_in.close();
        return line_seq;
    }
    else{
        cerr << "Error opening the input file (get_read_sequence)." << endl;
    }
    return "";
}

string get_read_header(const vector<uint64_t>& read_line_pos, const string& filename, uint64_t read_id){
    //Returns in a string the 'read_id'th header of the fasta file
    auto pos_seq = read_line_pos[(read_id*2)+1]+(read_id*2);
    auto pos_header = read_line_pos[(read_id)*2]+(read_id*2-1);
    if(read_id==0){
        pos_header=0;
        pos_seq=read_line_pos[1];
    }
    fstream file_in(filename, ios::in);
    //zstr::ifstream file_in(filename, ios::in);
    char line_seq[100000];
    if(file_in){
        file_in.seekg(pos_header, file_in.beg);
        file_in.read(line_seq, (pos_seq-pos_header));
        line_seq[pos_seq-pos_header-1] = 0;
        file_in.close();
        return line_seq;
    }
    else{
        cerr << "Error opening the input file (get_read_header)." << endl;
    }
    return "";
}

vector<int> order_ctgs(const string& gfa_file, int algo){
    //Creates a vector containing the numbers of the contigs in the right order (depth-first search if algo=1, breadth-first search if algo=2)
    if(algo!=1 && algo!=2){
        cerr << "Wrong algorithm number (order_ctgs). Must be 1 or 2" << endl;
    }
    //Fills a map with the graph : some contigs are linked with one or more other contigs, the info is found in the gfa file 
    ankerl::unordered_dense::map<int, vector<int>> links;
    int nb_contigs=0;

    typedef pair <int, vector<int>> Int_vector;
    string ctg1, ctg2;
    int num_ctg1, num_ctg2=0;
    char c;
    int i;
    vector<int>::iterator it;
    fstream gfa(gfa_file, ios::in);
    if(gfa){
        string line;
        while(!gfa.eof()){
            getline(gfa, line);
            vector<string> words_of_line = split_string(line, '\t');
            if(words_of_line.size()>0 && words_of_line[0]=="L"){
                //Retrieve the couple of contigs
                ctg1="";
                i=0;
                c=words_of_line[1][i];
                while(c=='u' || c=='t' || c=='g' || c=='0'){
                    i++;
                    c=words_of_line[1][i];
                }
                while(c!='l' && c!='c'){
                    ctg1+=c;
                    i++;
                    c=words_of_line[1][i];
                }
                num_ctg1 = stoi(ctg1);
                num_ctg1-=1;

                ctg2="";
                i=0;
                c=words_of_line[3][i];
                while(c=='u' || c=='t' || c=='g' || c=='0'){
                    i++;
                    c=words_of_line[3][i];
                }
                while(c!='l'){
                    ctg2+=c;
                    i++;
                    c=words_of_line[3][i];
                }
                num_ctg2 = stoi(ctg2);
                num_ctg2-=1;

                //Insert the couple in the map
                if(num_ctg1!=num_ctg2){
                    if(links.find(num_ctg1)!=links.end()){
                        //The first contig is a key so we add the second to his list if it isn't already
                        it=links[num_ctg1].begin();
                        while(it!=links[num_ctg1].end() && *it!=num_ctg2){
                            it++;
                        }
                        if(it==links[num_ctg1].end()){
                            links[num_ctg1].push_back(num_ctg2);
                        }
                    }
                    else{
                        //The first contig isn't a key so we add it to the map with a list containing the second contig
                        links.insert(Int_vector(num_ctg1, {num_ctg2}));
                    }
                }
            }
            if(words_of_line.size()>0 && words_of_line[0]=="S"){
                nb_contigs++;
            } 
        }
    }
    //Depth-first search :
    //The vector seen stores a 1 if the contig is already in the order, 0 otherwise
    //The vector order stores the resulting order of contigs
    vector<int> order;
    vector<int> seen(nb_contigs, 0);
    if(algo==1){
        stack<int> pilegraph;

        for (auto& elem : links){
            if(seen[elem.first]==0){
                seen[elem.first]=1;
                pilegraph.push(elem.first);
                while(!pilegraph.empty()){
                    num_ctg2 = pilegraph.top();
                    order.push_back(num_ctg2);
                    pilegraph.pop();
                    for(vector<int>::iterator ite = links[num_ctg2].begin(); ite!=links[num_ctg2].end(); ite++){
                        if(seen[*ite]==0){
                            seen[*ite]=1;
                            pilegraph.push(*ite);
                        }
                    }
                }
            }
        }
    }
    else{
        queue<int> pilegraph;

        for (auto& elem : links){
            if(seen[elem.first]==0){
                seen[elem.first]=1;
                pilegraph.push(elem.first);
                while(!pilegraph.empty()){
                    num_ctg2 = pilegraph.front();
                    order.push_back(num_ctg2);
                    pilegraph.pop();
                    for(vector<int>::iterator ite = links[num_ctg2].begin(); ite!=links[num_ctg2].end(); ite++){
                        if(seen[*ite]==0){
                            seen[*ite]=1;
                            pilegraph.push(*ite);
                        }
                    }
                }
            }
        }
    }

    return(order);
}

ankerl::unordered_dense::map<int, int> unmapped_reads_miniasm(int nb_reads, const string& gfa_file){
    //Returns a map containing the number of the reads that are not mapped on a contig in the gfa file
    typedef pair <int, int> Int_Pair;
    //The vector reads notes if a read is seen (1) or not (0)
    vector<int> reads(nb_reads, 0);
    fstream gfa(gfa_file, ios::in);
    if(gfa){
        string line;
        while(!gfa.eof()){
            getline(gfa, line);
            vector<string> words_of_line = split_string(line, '\t');
            if(words_of_line.size()>0 && words_of_line[0]=="a"){
                //The line begins with an "a" so it corresponds to a read, update of the vector reads
                string num_read = split_string(split_string(words_of_line[3], ':')[0], '.')[split_string(split_string(words_of_line[3], ':')[0], '.').size()-1]; 
                reads[stoul(num_read)] = 1;
            }
        }
    }
    //We browse the vector and write in a map the number of unmapped reads
    ankerl::unordered_dense::map<int, int> unmapped;
    for(int i=0; i<reads.size(); i++){
        if(reads[i]==0){
            unmapped.insert(Int_Pair(i, 0));
        }
    }
    cout << "reads non mappés : "+ to_string(unmapped.size()) << endl;
    return(unmapped);
}

ankerl::unordered_dense::map<int, vector<int>> link_unmapped_gz(int nb_reads, const string& gfa_file, const string& compressed_paf_file){
    //Returns a map that links a mapped read to a vector of unmapped reads according to the alignment in the paf file
    typedef pair <int, int> Int_Pair;
    typedef pair <int, vector<int>> Int_vector;
    //The map unmapped contains the reads that are not mapped and not associated with an other read
    ankerl::unordered_dense::map<int, int> unmapped = unmapped_reads_miniasm(nb_reads, gfa_file);
    //The map newly_mapped contains the reads that are not mapped in the gfa file but already associated with another read
    ankerl::unordered_dense::map<int, int> newly_mapped;
    //The map links contains the associations of mapped reads with a vector of unmapped reads
    ankerl::unordered_dense::map<int, vector<int>> links;
    int num_current_read;
    int current_read_align;
    zstr::ifstream paf(compressed_paf_file);
    if(paf){
        //Reads the paf file containing alignments between reads
        string line;
        getline(paf, line);
        while(!paf.eof()){
            //For each alignment, if the first read is unmapped and unassociated with an other, we look at the read aligned with it
            //If the read aligned with it is unmapped (associated or not with an other), we associate the 2 reads and marked it in 'newly_mapped'
            //Otherwise we take a look at other alignments until we find a mapped read to associate our unmapped read with
            vector<string> words_of_line = split_string(line, '\t');
            if(words_of_line.size()>0){
                num_current_read = stoi(split_string(words_of_line[0], '.')[split_string(words_of_line[0], '.').size()-1]);
                if(unmapped.find(num_current_read)!=unmapped.end()){
                    //The read is unmapped and unassociated
                    //We associate it with the first mapped read we find
                    current_read_align = stoi(split_string(words_of_line[5], '.')[split_string(words_of_line[5], '.').size()-1]); 
                    while(words_of_line.size()>0 && (unmapped.find(stoi(split_string(words_of_line[5], '.')[split_string(words_of_line[5], '.').size()-1]))!=unmapped.end() || newly_mapped.find(stoi(split_string(words_of_line[5], '.')[split_string(words_of_line[5], '.').size()-1]))!=newly_mapped.end()) && num_current_read == stoi(split_string(words_of_line[0], '.')[split_string(words_of_line[0], '.').size()-1]) && !paf.eof()){
                        getline(paf, line);
                        words_of_line = split_string(line, '\t');
                        if(words_of_line.size()>0 && num_current_read == stoi(split_string(words_of_line[0], '.')[split_string(words_of_line[0], '.').size()-1])){
                            current_read_align = stoi(split_string(words_of_line[5], '.')[split_string(words_of_line[5], '.').size()-1]);
                        }
                    }
                    //We reached the end of the bloc or we found a mapped read
                    if(unmapped.count(current_read_align)==0 && newly_mapped.count(current_read_align)==0){
                        //If the read is mapped, we associate the couple of reads in 'links'
                        if(links.find(current_read_align)==links.end()){
                            links.insert(Int_vector(current_read_align, {num_current_read}));
                        }
                        else{
                            links[current_read_align].push_back(num_current_read);
                        }
                        unmapped.erase(num_current_read);
                        newly_mapped.insert(Int_Pair(num_current_read, 0));
                    }
                    
                }
                while(words_of_line.size()>0 && num_current_read == stoi(split_string(words_of_line[0], '.')[split_string(words_of_line[0], '.').size()-1]) && !paf.eof()){
                    getline(paf, line);
                    words_of_line = split_string(line, '\t');
                }
            }
        }
    }
    //Writes the last reads that are not mapped in the map at the index -1
    links.insert(Int_vector(-1, {}));
    for (auto& elem : unmapped){
        links[-1].push_back(elem.first);
    }
    return links;
}

vector<int> unincluded_ctgs(int nb_ctgs, vector<int> ctgs_order){
    //Returns a vector containing the numbers of the contigs that are not in the contigs order (because they are not linked in the graph)
    vector<int> ctgs(nb_ctgs, 0);
    vector<int> lasting_ctgs;

    for(int i:ctgs_order){
        ctgs[i]=1;
    }
    for(int i=0; i<ctgs.size(); i++){
        if(ctgs[i]==0){
            lasting_ctgs.push_back(i);
        }
    }

    cout << "solo ctgs "+to_string(lasting_ctgs.size()) << endl;
    
    return(lasting_ctgs);
}

vector<int> shuffle_vector(vector<int> vector_to_shuffle){
    srand((unsigned) time(NULL));
    int index=0;
    for(int i = 0; i<vector_to_shuffle.size(); i++){
        index=rand()%vector_to_shuffle.size();
        int tmp = vector_to_shuffle[index];
        vector_to_shuffle[index]=vector_to_shuffle[i];
        vector_to_shuffle[i] = tmp;
    }
    return(vector_to_shuffle);
}

void ctgs_and_unmapped_sorting_compressed(const string& reads_file, const string& gfa_file, const string& out_file, const string& paf_file, int ctgs_sort){
    //Writes a new file containing the same reads as 'reads_file', in a different order : 
    //the reads are in the same order as they are mapped to contigs, contigs are put in a depth-first search order, and the reads that are not in any contigs are brougth together to similar reads in the contigs
    typedef pair <int, int> Int_Pair;
    vector<uint64_t> read_line_pos = read_line_position(reads_file);
    int nb_reads = read_line_pos.size()/2-1;
    ankerl::unordered_dense::map<int, vector<int>> unmapped_align = link_unmapped_gz(nb_reads, gfa_file, paf_file);
    
    vector<uint64_t> contig_line_pos = contig_line_position(gfa_file);
    int nb_contigs = contig_line_pos.size();
    int mapped=0;
    int non_map = 0;
    int reads_ecrits=0;
    vector<int> ctgs(nb_contigs);
    for(int i=0; i<nb_contigs; i++){
        ctgs[i]=i;
    }
    //Ordering the contigs with the chosen algorithm (in ctgs_sort)
    vector<int> ctgs_sorted;
    if(ctgs_sort==0){
        //Random
        ctgs_sorted = shuffle_vector(ctgs);
    }
    else if(ctgs_sort==1){
        //depth-first
        ctgs_sorted = order_ctgs(gfa_file, 1);
    }
    else{
        //breadth-first
        ctgs_sorted = order_ctgs(gfa_file, 2);
    }
    
    cout << "number of contigs : "+to_string(nb_contigs) << endl;
    vector<int> lasting_ctgs = unincluded_ctgs(nb_contigs, ctgs_sorted);
    fstream gfa(gfa_file, ios::in);
    fstream out(out_file, ios::out);
    //zstr::ifstream reads(reads_file, ios::in);
    fstream reads(reads_file, ios::in);
    if(gfa && out && reads){
        string line;
        //Writes the contigs that are in the order
        for (int i: ctgs_sorted){
            gfa.clear();
            gfa.seekg(contig_line_pos[i], gfa.beg);
            getline(gfa, line);
            //Line begining with "S" that is not interesting, we copy only the ones that start with "a"
            getline(gfa, line);
            vector<string> words_of_line = split_string(line, '\t');
            while(words_of_line.size()>0 && (words_of_line[0]=="a" || words_of_line[0]=="L")){
                if(words_of_line[0]=="a"){
                    string num_read = split_string(split_string(words_of_line[3], ':')[0], '.')[1];
                    out << get_read_header(read_line_pos, reads_file, stoul(num_read))+"\n"+get_read_sequence(read_line_pos, reads_file, stoul(num_read)) << endl;
                    reads_ecrits++;
                    //Insertion of unmapped reads linked to this read
                    if(unmapped_align.find(stoi(num_read))!=unmapped_align.end()){
                        for(vector<int>::iterator it = unmapped_align[stoi(num_read)].begin(); it!=unmapped_align[stoi(num_read)].end(); it++){
                            out << get_read_header(read_line_pos, reads_file, *it)+"\n"+get_read_sequence(read_line_pos, reads_file, *it) << endl;
                            reads_ecrits++;
                            mapped++;
                        }
                        unmapped_align.erase(stoi(num_read));
                    }
                }
                getline(gfa, line);
                words_of_line = split_string(line, '\t');
            }
        }
        //Writes the contigs that are alone
        for (int i: lasting_ctgs){
            gfa.clear();
            gfa.seekg(contig_line_pos[i], gfa.beg);
            getline(gfa, line);
            //Line begining with "S" that is not interesting, we copy only the ones that start with "a"
            getline(gfa, line);
            vector<string> words_of_line = split_string(line, '\t');
            while(words_of_line.size()>0 && (words_of_line[0]=="a" || words_of_line[0]=="L")){
                if(words_of_line[0]=="a"){
                   string num_read = split_string(split_string(words_of_line[3], ':')[0], '.')[1];
                    out << get_read_header(read_line_pos, reads_file, stoul(num_read))+"\n"+get_read_sequence(read_line_pos, reads_file, stoul(num_read)) << endl;
                    reads_ecrits++;
                    //Insertion of unmapped reads linked to this read
                    if(unmapped_align.find(stoi(num_read))!=unmapped_align.end()){
                        for(vector<int>::iterator it = unmapped_align[stoi(num_read)].begin(); it!=unmapped_align[stoi(num_read)].end(); it++){
                            out << get_read_header(read_line_pos, reads_file, *it)+"\n"+get_read_sequence(read_line_pos, reads_file, *it) << endl;
                            reads_ecrits++;
                            mapped++;
                        }
                        unmapped_align.erase(stoi(num_read));
                    }
                }
                getline(gfa, line);
                words_of_line = split_string(line, '\t'); 
            }
        }
        
        //Write the reads that are still alone
        for(vector<int>::iterator it = unmapped_align[-1].begin(); it!=unmapped_align[-1].end(); it++){
            out << get_read_header(read_line_pos, reads_file, *it)+"\n"+get_read_sequence(read_line_pos, reads_file, *it) << endl;
            reads_ecrits++;
            non_map++;
        }
        cout << "reads mappés par la fonction : "+to_string(mapped) << endl;
        cout << "reads toujours non mappés : "+to_string(non_map) << endl;
        cout << "reads ecrits : "+to_string(reads_ecrits)<<endl;
    }
}

int main(int argc, char *argv[])
{
    string readsfile;
    string gfafile;
    string outfile;
    string paffile;
    int ctgs_sort=1;
    if (argc<5){
        cerr << "Usage : ./reads_sorting [READS FILE] [GFA FILE] [OUT FILE] [PAF FILE] [ctgs sorting algo]" << endl;
        return EXIT_FAILURE;
    }
    else{
        readsfile = argv[1];
        gfafile = argv[2];
        outfile = argv[3];
        paffile = argv[4];
        if(argc==6){
            //Option : ctgs sort 0:random order, 1:depth-first search, 2:breadth-first search. Default=1
            ctgs_sort = stoi(argv[5]);
        }
    }

    ctgs_and_unmapped_sorting_compressed(readsfile, gfafile, outfile, paffile, ctgs_sort);

    return 0;
}
