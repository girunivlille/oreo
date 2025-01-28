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
#include <regex>
#include <algorithm>
#include "unordered_dense.h"

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

string reverse_complement(const string& sequence){
    string result;

    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        switch (*it) {
            case 'A': result += 'T'; break;
            case 'T': result += 'A'; break;
            case 'C': result += 'G'; break;
            case 'G': result += 'C'; break;
            case 'a': result += 't'; break;
            case 't': result += 'a'; break;
            case 'c': result += 'g'; break;
            case 'g': result += 'c'; break;
            default: result += *it; break;
        }
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
            total_count += line.length()+1;
            if (line!=""){
                read_line_pos.push_back(total_count); 
            }
        }
        //read_line_pos.pop_back(); <- finalement utile pour récupérer la taille de la dernière séquence
        return read_line_pos;
    }
    else{
        cerr << "Error opening the file (read_line_position)." << endl;
    }
    return {};
}

vector<uint64_t> contig_line_position(const string& paf_file, int nb_contigs, int rev_comp){
    //Notes in a vector the position of the begining of each contig in a paf file
    uint64_t nb_contigs_positions = (rev_comp == 0) ? nb_contigs * 2 : nb_contigs;
    vector<uint64_t> contig_line_pos(nb_contigs_positions, 0);
    fstream paf(paf_file, ios::in);
    if(paf){
        string line;
        uint64_t total_count = 0;
        uint64_t ctg_nb = 0;
        string ctg = "";
        std::regex pattern(R"(utg0+(\d+)[lc])");
        std::smatch match;
        while(!paf.eof()){
            getline(paf, line);
            vector<string> words_of_line = split_string(line, '\t');
            if(words_of_line.size()>0 && words_of_line[5]!=ctg){
                ctg = words_of_line[5];
                if (std::regex_match(words_of_line[5], match, pattern)){
                    ctg_nb = std::stoi(match[1])-1;
                }
                if (rev_comp==0 && words_of_line[4]=="-"){
                    contig_line_pos[ctg_nb+nb_contigs] = total_count;
                }
                else{
                    contig_line_pos[ctg_nb] = total_count;
                }
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

string get_read_part(const vector<uint64_t>& read_line_pos, const string& filename, const string& format, uint64_t read_num, const string& part){
    //Returns in a string the 'read_id'th part (header, sequence, separator, qscore) of the fasta/fastq file
    int block_size = 2;
    if (format=="fastq"){
        block_size = 4;
    }

    int part_num = 0;
    if (part=="sequence"){
        part_num = 1;
    }
    else if (part=="separator"){
        part_num = 2;
    }
    else if (part=="qscore"){
        part_num = 3;
    }

    auto pos_part = read_line_pos[(read_num)*block_size+part_num];
    auto pos_next_line = read_line_pos[(read_num*block_size)+1+part_num];

    fstream file_in(filename, ios::in);
    char line_seq[5000000];
    if(file_in){
        file_in.seekg(pos_part, file_in.beg);
        file_in.read(line_seq, (pos_next_line-pos_part));
        line_seq[pos_next_line-pos_part-1] = 0;
        file_in.close();
        return line_seq;
    }
    else{
        cerr << "Error opening the input file (get_read_header)." << endl;
    }
    return "";
}

vector<int> order_ctgs(const string& gfa_links_file, int nb_contigs, int algo){
    //Creates a vector containing the numbers of the contigs in the right order (depth-first search if algo=1, breadth-first search if algo=2)
    if(algo!=1 && algo!=2){
        cerr << "Wrong algorithm number (order_ctgs). Must be 1 or 2" << endl;
    }
    //Fills a map with the graph : some contigs are linked with one or more other contigs, the info is found in the gfa file 
    ankerl::unordered_dense::map<int, vector<int>> links;

    typedef pair <int, vector<int>> Int_vector;
    string ctg1, ctg2;
    int num_ctg1, num_ctg2=0;
    char c;
    int i;
    std::regex pattern(R"(utg0+(\d+)[lc])");
    std::smatch match1, match2;
    vector<int>::iterator it;
    fstream gfa(gfa_links_file, ios::in);
    if(gfa){
        string line;
        while(!gfa.eof()){
            getline(gfa, line);
            vector<string> words_of_line = split_string(line, '\t');
            if(words_of_line.size()>0){
                //Retrieve the couple of contigs
                if (std::regex_match(words_of_line[1], match1, pattern)){
                    num_ctg1 = std::stoi(match1[1])-1;
                }
                if (std::regex_match(words_of_line[3], match2, pattern)){
                    num_ctg2 = std::stoi(match2[1])-1;
                }

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

void write_read_in_file(int read_num, fstream& out, const string& reads_file, int rev_comp, const string& strand, const string& format, const vector<uint64_t>& read_line_pos){
    //écrire le read selon son numéro (champ 0) (si fastq écrire 4 lignes, si rev_comp et orientation = - (champ 4) retourner la seq)
    if (rev_comp==1 && strand=="-"){
        out << get_read_part(read_line_pos, reads_file, format, read_num, "header")+"\n"+reverse_complement(get_read_part(read_line_pos, reads_file, format, read_num, "sequence")) << endl;
    }
    else{
        out << get_read_part(read_line_pos, reads_file, format, read_num, "header")+"\n"+get_read_part(read_line_pos, reads_file, format, read_num, "sequence") << endl;
    }
    if (format=="fastq"){
        if (rev_comp==1 && strand=="-"){
            string qscore = get_read_part(read_line_pos, reads_file, format, read_num, "qscore");
            reverse(qscore.begin(), qscore.end());
            out << get_read_part(read_line_pos, reads_file, format, read_num, "separator")+"\n"+qscore << endl;
        }
        else{
            out << get_read_part(read_line_pos, reads_file, format, read_num, "separator")+"\n"+get_read_part(read_line_pos, reads_file, format, read_num, "qscore") << endl;
        }
    }
}

void write_all_reads_from_ctg(fstream& out, fstream& paf, const string& reads_file, int rev_comp, const string& format, const vector<uint64_t>& read_line_pos, const vector<uint64_t>& contig_line_pos, int ctg_written_num, int& mapped, int& reads_ecrits, vector<int>& written_reads, vector<int>& reads_in_ctg, string strand_to_write, int index_to_add_if_strand_minus){
    regex pattern(R"(utg0+(\d+)[lc])");
    string line;
    string strand;
    int ctg_num;
    int read_num;

    paf.clear();
    paf.seekg(contig_line_pos[index_to_add_if_strand_minus + ctg_written_num], paf.beg);
    // lire les lignes et écrire les reads tant que numéro de contig (champ 5) == i
    getline(paf, line);
    vector<string> words_of_line = split_string(line, '\t');
    //extraire le numéro de contig
    ctg_num = -1;
    smatch match;
    if (regex_match(words_of_line[5], match, pattern)){
        ctg_num = std::stoi(match[1])-1;
    }
    while (!paf.eof() && ctg_num==ctg_written_num && (strand_to_write=="both" || (strand_to_write!="both" && words_of_line[4]==strand_to_write))){
        read_num = stoi(words_of_line[0]);
        strand = words_of_line[4];
        // si le read n'a pas encore été écrit
        if (written_reads[read_num]==0){
            //écrire le read selon son numéro (champ 0) (si fastq écrire 4 lignes, si rev_comp et orientation = - (champ 4) retourner la seq)
            write_read_in_file(read_num, out, reads_file, rev_comp, strand, format, read_line_pos);
            //noter quelque part que le read a été écrit et augmenter le compteur du nombre de reads écrits
            written_reads[read_num] = 1;
            mapped++;
            reads_ecrits++;
            reads_in_ctg[ctg_num]++;
        }
        //lire la ligne suivante
        getline(paf, line);
        words_of_line = split_string(line, '\t');
        if (words_of_line.size()>0 && regex_match(words_of_line[5], match, pattern)){
            ctg_num = std::stoi(match[1])-1;
        }
    }
}

void ctgs_and_unmapped_sorting_compressed(const string& reads_file, const string& format, const string& gfa_links_file, int nb_ctgs, const string& out_file, const string& paf_file, const string& log_file, int rev_comp, int ctgs_sort){
    //Writes a new file containing the same reads as 'reads_file', in a different order : 
    //the reads are in the same order as they are mapped to contigs, contigs are put in a depth-first search order, and the reads that are not in any contigs are brougth together at the end of the file

    vector<uint64_t> read_line_pos = read_line_position(reads_file);
    vector<uint64_t> contig_line_pos = contig_line_position(paf_file, nb_ctgs, rev_comp);


    int nb_reads = (read_line_pos.size() - 1) / (format == "fasta" ? 2 : 4);

    //Ordering the contigs with the chosen algorithm (in ctgs_sort)
    vector<int> ctgs(nb_ctgs);
    for(int i=0; i<nb_ctgs; i++){
        ctgs[i]=i;
    }
    vector<int> ctgs_sorted;
    if(ctgs_sort==0){
        //Random
        ctgs_sorted = shuffle_vector(ctgs);
    }
    else if(ctgs_sort==1){
        //depth-first
        ctgs_sorted = order_ctgs(gfa_links_file, nb_ctgs, 1);
    }
    else{
        //breadth-first
        ctgs_sorted = order_ctgs(gfa_links_file, nb_ctgs, 2);
    }
    
    vector<int> lasting_ctgs = unincluded_ctgs(nb_ctgs, ctgs_sorted);

    int mapped=0;
    int non_map = 0;
    int reads_ecrits=0;

    vector<int> written_reads(nb_reads, 0);
    vector<int> reads_in_ctg(nb_ctgs, 0);

    std::regex pattern(R"(utg0+(\d+)[lc])");

    fstream paf(paf_file, ios::in);
    fstream out(out_file, ios::out);
    fstream reads(reads_file, ios::in);
    if(paf && out && reads){
        //Writes the contigs that are in the order
        if (rev_comp==0){
            // Write the strand +
            for (int i: ctgs_sorted){
                write_all_reads_from_ctg(out, paf, reads_file, rev_comp, format, read_line_pos, contig_line_pos, i, mapped, reads_ecrits, written_reads, reads_in_ctg, "+", 0);
            }
            // Write the strand -
            for (int i: ctgs_sorted){
                write_all_reads_from_ctg(out, paf, reads_file, rev_comp, format, read_line_pos, contig_line_pos, i, mapped, reads_ecrits, written_reads, reads_in_ctg, "-", nb_ctgs);
            }
        }
        else{
            for (int i: ctgs_sorted){
                write_all_reads_from_ctg(out, paf, reads_file, rev_comp, format, read_line_pos, contig_line_pos, i, mapped, reads_ecrits, written_reads, reads_in_ctg, "both", 0);
            }
        }

        //Writes the contigs that are alone
        if (rev_comp==0){
            // Write the strand +
            for (int i: lasting_ctgs){
                write_all_reads_from_ctg(out, paf, reads_file, rev_comp, format, read_line_pos, contig_line_pos, i, mapped, reads_ecrits, written_reads, reads_in_ctg, "+", 0);
            }
            // Write the strand -
            for (int i: lasting_ctgs){
                write_all_reads_from_ctg(out, paf, reads_file, rev_comp, format, read_line_pos, contig_line_pos, i, mapped, reads_ecrits, written_reads, reads_in_ctg, "-", nb_ctgs);
            }
        }
        else{
            for (int i: lasting_ctgs){
                write_all_reads_from_ctg(out, paf, reads_file, rev_comp, format, read_line_pos, contig_line_pos, i, mapped, reads_ecrits, written_reads, reads_in_ctg, "both", 0);
            }
        }

        //Write the reads that are still alone
        for (int i=0; i<nb_reads; i++){
            if (written_reads[i]==0){
                out << get_read_part(read_line_pos, reads_file, format, i, "header")+"\n"+get_read_part(read_line_pos, reads_file, format, i, "sequence") << endl;
                if (format=="fastq"){
                    out << get_read_part(read_line_pos, reads_file, format, i, "separator")+"\n"+get_read_part(read_line_pos, reads_file, format, i, "qscore") << endl;
                }
                non_map++;
                reads_ecrits++;
            }
        }

        fstream log(log_file, ios::out | ios::app);
        log << "reads mappés par la fonction : "+to_string(mapped) << endl;
        log << "reads toujours non mappés : "+to_string(non_map) << endl;
        log << "reads ecrits : "+to_string(reads_ecrits)<<endl;

        cout << "reads mappés par la fonction : "+to_string(mapped) << endl;
        cout << "reads toujours non mappés : "+to_string(non_map) << endl;
        cout << "reads ecrits : "+to_string(reads_ecrits)<<endl;
    }
}

int main(int argc, char *argv[])
{
    string readsfile;
    string format;
    string gfafile;
    string gfa_links;
    int nb_ctgs;
    string outfile;
    string paffile;
    string logfile;
    int rev_comp;
    int ctgs_sort=1;
    if (argc<9){
        cerr << "Usage : ./reads_sorting [READS FILE] [FORMAT] [GFA LINKS FILE] [NB OF CTGS] [OUT FILE] [PAF FILE] [LOG FILE] [REVERSE COMPLEMENT] [ctgs sorting algo]" << endl;
        return EXIT_FAILURE;
    }
    else{
        readsfile = argv[1];
        format = argv[2];
        gfa_links = argv[3];
        nb_ctgs = stoi(argv[4]);
        outfile = argv[5];
        paffile = argv[6];
        logfile = argv[7];
        rev_comp = stoi(argv[8]);
        if(argc==10){
            //Option : ctgs sort 0:random order, 1:depth-first search, 2:breadth-first search. Default=1
            ctgs_sort = stoi(argv[9]);
        }
    }

    ctgs_and_unmapped_sorting_compressed(readsfile, format, gfa_links, nb_ctgs, outfile, paffile, logfile, rev_comp, ctgs_sort);

    return 0;
}

