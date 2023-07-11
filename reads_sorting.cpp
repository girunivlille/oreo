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
#include <stack>
#include "unordered_dense.h"
#include "zstr/src/zstr.hpp"

using namespace std;

vector<string> split_string(const string& str, char separator) {
    vector<string> result;
    istringstream iss(str);
    string word;

    while (getline(iss, word, separator)) {
        result.push_back(word);
    }

    return result;
}

vector<uint32_t> read_line_position(const string& filename){
    vector<uint32_t> read_line_pos={0};
    fstream file(filename, ios::in);
    if(file){
        string line;
        uint32_t total_count = 0;
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

vector<uint32_t> contig_line_position(const string& gfa_file){
    vector<uint32_t> contig_line_pos={0};
    fstream gfa(gfa_file, ios::in);
    if(gfa){
        string line;
        uint32_t total_count = 0;
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

string get_read_sequence(const vector<uint32_t>& read_line_pos, const string& filename, uint32_t read_id){
    auto pos_seq = read_line_pos[(read_id*2)+1];
    auto pos_entete_suivant = read_line_pos[(read_id*2)+2];

    fstream file_in(filename, ios::in);
    char line_seq[100000];
    if(file_in){
        file_in.seekg((pos_seq+read_id*2), file_in.beg);
        file_in.read(line_seq, ((pos_entete_suivant+read_id)-(pos_seq+read_id)));
        //line_seq[((pos_entete_suivant+read_id*2)-(pos_seq+read_id))] = 0;
        line_seq[pos_entete_suivant-pos_seq+2] = 0;
        file_in.close();
        return line_seq;
    }
    else{
        cerr << "Error opening the input file (get_read_sequence)." << endl;
    }
    return "";
}

string get_read_header(const vector<uint32_t>& read_line_pos, const string& filename, uint32_t read_id){
    auto pos_seq = read_line_pos[(read_id*2)+1]+(read_id*2);
    auto pos_header = read_line_pos[(read_id)*2]+(read_id*2-1);
    if(read_id==0){
        pos_header=0;
        pos_seq=read_line_pos[1];
    }
    fstream file_in(filename, ios::in);
    char line_seq[100000];
    if(file_in){
        file_in.seekg(pos_header, file_in.beg);
        file_in.read(line_seq, (pos_seq-pos_header));
        //line_seq[((pos_entete_suivant+read_id*2)-(pos_seq+read_id))] = 0;
        line_seq[pos_seq-pos_header-1] = 0;
        file_in.close();
        //cout << line_seq << endl;
        return line_seq;
    }
    else{
        cerr << "Error opening the input file (get_read_header)." << endl;
    }
    return "";
}

vector<int> order_ctgs(const string& gfa_file){
    //parcourir le gfa et créer une map int, vector<int> (à éléments uniques) pour stocker les liens entre contigs
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
                //récupération du couple de contigs 
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

                //insertion du couple de contigs dans links (seulement si le contig2 n'est pas déjà dans la liste du ctg1 et que ctg1!=ctg2)
                if(num_ctg1!=num_ctg2){
                    if(links.find(num_ctg1)!=links.end()){
                        //ctg1 est trouvé donc on insère ctg2 dans sa liste s'il n'y est pas déjà
                        it=links[num_ctg1].begin();
                        while(it!=links[num_ctg1].end() && *it!=num_ctg2){
                            it++;
                        }
                        if(it==links[num_ctg1].end()){
                            //on a pas trouvé ctg2 donc on peut l'ajouter dans la liste
                            links[num_ctg1].push_back(num_ctg2);
                        }
                    }
                    else{
                        //ctg1 n'est pas trouvé donc on le créé dans links avec ctg2 dans sa liste
                        links.insert(Int_vector(num_ctg1, {num_ctg2}));
                    }
                }
            }
            if(words_of_line.size()>0 && words_of_line[0]=="S"){
                nb_contigs++;
            } 
        }
    }
    //Affichage de la map de liens du graphe
    /*for (auto& elem : links){
        cout << "first "+to_string(elem.first+1) << endl;
        for(vector<int>::iterator iter = elem.second.begin(); iter!=elem.second.end(); iter++){
            cout << *iter+1 << endl;
        }
    }*/
    //créer un vecteur<int> de la taille nb_ctgs qui dit si on a déjà classé le contig
    //créer une pile
    //créer un vecteur<int> qui renvoie le tri des contigs 
    vector<int> seen(nb_contigs, 0);
    stack<int> pilegraph;
    vector<int> order;

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

    /*for(int i: order){
        cout << to_string(i+1) << " ";
    }*/

    return(order);
}

ankerl::unordered_dense::map<int, int> unmapped_reads_miniasm(int nb_reads, const string& gfa_file){
    typedef pair <int, int> Int_Pair;
    //renvoie les numéros des reads qui n'ont pas été intégrés dans le gfa, dans une map
    
    vector<int> reads(nb_reads, 0);
    fstream gfa(gfa_file, ios::in);
    if(gfa){
        string line;
        while(!gfa.eof()){
            getline(gfa, line);
            vector<string> words_of_line = split_string(line, '\t');
            if(words_of_line.size()>0 && words_of_line[0]=="a"){
                //words_of_line[7][2]
                //string num_read = words_of_line[4];
                string num_read = split_string(split_string(words_of_line[3], ':')[0], '.')[split_string(split_string(words_of_line[3], ':')[0], '.').size()-1]; 
                reads[stoul(num_read)] = 1;
            }
        }
    }
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
    //renvoie une map qui lie un read mappé à un read non mappé
    typedef pair <int, int> Int_Pair;
    typedef pair <int, vector<int>> Int_vector;
    ankerl::unordered_dense::map<int, int> unmapped = unmapped_reads_miniasm(nb_reads, gfa_file);
    ankerl::unordered_dense::map<int, int> newly_mapped;
    ankerl::unordered_dense::map<int, vector<int>> links;
    int num_current_read;
    int current_read_align;
    zstr::ifstream paf(compressed_paf_file);
    //fstream paf(compressed_paf_file, ios::in);
    if(paf){
        string line;
        getline(paf, line);
        while(!paf.eof()){
            vector<string> words_of_line = split_string(line, '\t');
            if(words_of_line.size()>0){
                num_current_read = stoi(split_string(words_of_line[0], '.')[split_string(words_of_line[0], '.').size()-1]);
                if(unmapped.find(num_current_read)!=unmapped.end()){
                    //le read est non mappé
                    //on l'associe au premier read aligné qui est mappé qu'on trouve
                    current_read_align = stoi(split_string(words_of_line[5], '.')[split_string(words_of_line[5], '.').size()-1]); 
                    while(words_of_line.size()>0 && (unmapped.find(stoi(split_string(words_of_line[5], '.')[split_string(words_of_line[5], '.').size()-1]))!=unmapped.end() || newly_mapped.find(stoi(split_string(words_of_line[5], '.')[split_string(words_of_line[5], '.').size()-1]))!=newly_mapped.end()) && num_current_read == stoi(split_string(words_of_line[0], '.')[split_string(words_of_line[0], '.').size()-1]) && !paf.eof()){
                        getline(paf, line);
                        words_of_line = split_string(line, '\t');
                        if(words_of_line.size()>0 && num_current_read == stoi(split_string(words_of_line[0], '.')[split_string(words_of_line[0], '.').size()-1])){
                            current_read_align = stoi(split_string(words_of_line[5], '.')[split_string(words_of_line[5], '.').size()-1]);
                        }
                    }
                    //si on est arrivé à la fin du bloc sans trouver un read mappé à associer à notre read mappé
                    //càd current_read_align est non mappé malgré tous nos efforts :( alors on ne fait rien 
                    //si current read est mappé càd il n'est ni dans unmapped ni dans newly_mapped
                    if(unmapped.find(current_read_align)==unmapped.end() && newly_mapped.find(current_read_align)==newly_mapped.end()){
                        if(links.find(current_read_align)==links.end()){
                            //le read n'était pas déjà associé à un non mappé
                            links.insert(Int_List(current_read_align, {num_current_read}));
                        }
                        else{
                            //le read était déjà dans la map
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
    /*cout << "links ---------"<<endl;
    for (auto& elem : links){
        cout << "cle "+to_string(elem.first) << endl;
        for(vector<int>::iterator it = elem.second.begin(); it!=elem.second.end(); it++){
            cout << to_string(*it) +" ";
        }
        cout << "\n";
    }
    cout << "unmapped ----------"<<endl;
    for (auto& elem : unmapped){
        cout << to_string(elem.first) +" ";
    }
    cout << "\nnewly_mapped -------------"<<endl;
    for (auto& elem : newly_mapped){
        cout << to_string(elem.first) +" ";
    }
    cout <<"\n";*/

    return links;
}

vector<int> unincluded_ctgs(int nb_ctgs, vector<int> ctgs_order){
    //renvoie les numéros des contigs qui n'ont pas été intégrés dans l'ordre des contigs, dans un vecteur
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

void ctgs_and_unmapped_sorting_compressed(const string& reads_file, const string& gfa_file, const string& out_file, const string& paf_file){
    //les contigs sont triés aléatoirement et chaque read restant est inséré après le read sur lequel il a le meilleur alignement (d'après le paf)
    typedef pair <int, int> Int_Pair;
    vector<uint32_t> read_line_pos = read_line_position(reads_file);
    int nb_reads = read_line_pos.size()/2-1;
    ankerl::unordered_dense::map<int, vector<int>> unmapped_align = link_unmapped_gz(nb_reads, gfa_file, paf_file);
    
    vector<uint32_t> contig_line_pos = contig_line_position(gfa_file);
    int nb_contigs = contig_line_pos.size();
    int mapped=0;
    int non_map = 0;
    vector<int> ctgs(nb_contigs);
    for(int i=0; i<nb_contigs; i++){
        ctgs[i]=i;
    }
    vector<int> ctgs_sorted = order_ctgs(gfa_file);
    cout << "number of contigs : "+to_string(nb_contigs) << endl;
    vector<int> lasting_ctgs = unincluded_ctgs(nb_contigs, ctgs_sorted);
    fstream gfa(gfa_file, ios::in);
    fstream out(out_file, ios::out);
    fstream reads(reads_file, ios::in);
    if(gfa && out && reads){
        string line;
        for (int i: ctgs_sorted){
            gfa.clear();
            gfa.seekg(contig_line_pos[i], gfa.beg);
            getline(gfa, line);
            //techniquement cette ligne c'est le "S" donc poubelle
            getline(gfa, line);
            vector<string> words_of_line = split_string(line, '\t');
            while(words_of_line.size()>0 && words_of_line[0]=="a"){
                string num_read = split_string(split_string(words_of_line[3], ':')[0], '.')[1];
                //string num_read = words_of_line[4];
                out << get_read_header(read_line_pos, reads_file, stoul(num_read))+"\n"+get_read_sequence(read_line_pos, reads_file, stoul(num_read)) << endl;
                //s'il y a des reads non mappés à insérer ici, on les insère
                if(unmapped_align.find(stoi(num_read))!=unmapped_align.end()){
                    for(vector<int>::iterator it = unmapped_align[stoi(num_read)].begin(); it!=unmapped_align[stoi(num_read)].end(); it++){
                        out << get_read_header(read_line_pos, reads_file, *it)+"\n"+get_read_sequence(read_line_pos, reads_file, *it) << endl;
                        mapped++;
                    }
                    unmapped_align[stoi(num_read)] = {};
                }
                getline(gfa, line);
                words_of_line = split_string(line, '\t');
            }
        }
        for (int i: lasting_ctgs){
            gfa.clear();
            gfa.seekg(contig_line_pos[i], gfa.beg);
            getline(gfa, line);
            //techniquement cette ligne c'est le "S" donc poubelle
            getline(gfa, line);
            vector<string> words_of_line = split_string(line, '\t');
            while(words_of_line.size()>0 && words_of_line[0]=="a"){
                string num_read = split_string(split_string(words_of_line[3], ':')[0], '.')[1];
                //string num_read = words_of_line[4];
                out << get_read_header(read_line_pos, reads_file, stoul(num_read))+"\n"+get_read_sequence(read_line_pos, reads_file, stoul(num_read)) << endl;
                //s'il y a des reads non mappés à insérer ici, on les insère
                if(unmapped_align.find(stoi(num_read))!=unmapped_align.end()){
                    for(vector<int>::iterator it = unmapped_align[stoi(num_read)].begin(); it!=unmapped_align[stoi(num_read)].end(); it++){
                        out << get_read_header(read_line_pos, reads_file, *it)+"\n"+get_read_sequence(read_line_pos, reads_file, *it) << endl;
                        mapped++;
                    }
                    unmapped_align[stoi(num_read)] = {};
                }
                getline(gfa, line);
                words_of_line = split_string(line, '\t');
            }
        }
        for (auto& elem : unmapped_align){
            for(vector<int>::iterator it = elem.second.begin(); it!=elem.second.end(); it++){
                out << get_read_header(read_line_pos, reads_file, *it)+"\n"+get_read_sequence(read_line_pos, reads_file, *it) << endl;
                non_map++;
            }
        }
        cout << "reads mappés par la fonction : "+to_string(mapped) << endl;
        cout << "reads toujours non mappés : "+to_string(non_map) << endl;
    }
}

int main(int argc, char *argv[])
{
    string readsfile;
    string gfafile;
    string outfile;
    string paffile;
    if (argc!=5){
        cerr << "Mauvais nombre d'arguments" << endl;
        return EXIT_FAILURE;
    }
    else{
        readsfile = argv[1];
        gfafile = argv[2];
        outfile = argv[3];
        paffile = argv[4];
    }

    ctgs_and_unmapped_sorting_compressed(readsfile, gfafile, outfile, paffile);

    /*for(uint32_t i = 0; i<11; i++){
        get_read_header(read_line_position("test_reads2.fa"), "test_reads2.fa", i);
    }*/
    


    return 0;
}