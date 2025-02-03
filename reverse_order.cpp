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
        return read_line_pos;
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

void write_read_in_file(int read_num, fstream& out, const string& reads_file, int rev_comp, const string& strand, const string& format, const vector<uint64_t>& read_line_pos, fstream& reverse_order){
    //écrire le read selon son numéro (champ 0) (si fastq écrire 4 lignes, si rev_comp et orientation = - (champ 4) retourner la seq)
    if (rev_comp==1 && strand=="-"){
        out << get_read_part(read_line_pos, reads_file, format, read_num, "header")+"\n"+reverse_complement(get_read_part(read_line_pos, reads_file, format, read_num, "sequence")) << endl;
        reverse_order << to_string(read_num)+" -" << endl;
    }
    else{
        out << get_read_part(read_line_pos, reads_file, format, read_num, "header")+"\n"+get_read_part(read_line_pos, reads_file, format, read_num, "sequence") << endl;
        reverse_order << to_string(read_num)+" +" << endl;
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

void reverse_file_order(const string& sorted_read_file, const string& reverse_order_file, const string& output_file, const string& format){
    vector<uint64_t> read_line_pos = read_line_position(sorted_read_file);
    int nb_reads = (read_line_pos.size() - 1) / (format == "fasta" ? 2 : 4);

    vector<int> reads_permutation(nb_reads, 0);
    vector<char> reads_strands(nb_reads, '+');

    // Retrieve read order
    fstream file(reverse_order_file, ios::in);
    if(file){
        string line;
        uint64_t total_count = 0;

        int read_nb;
        char strand;
        
        while(!file.eof()){
            getline(file, line);
            if (line!=""){
                istringstream flux(line);
                flux >> read_nb >> strand;
                reads_permutation[read_nb] = total_count;
                reads_strands[read_nb] = strand;
                total_count++;
            }
        }
    }

    // Write reads in said order
    fstream outfile(output_file, ios::out);
    if(outfile){
        int i = 0;
        while (i<nb_reads){
            if (reads_strands[i]=='-'){
                outfile << get_read_part(read_line_pos, sorted_read_file, format, reads_permutation[i], "header")+"\n"+reverse_complement(get_read_part(read_line_pos, sorted_read_file, format, reads_permutation[i], "sequence")) << endl;
                if (format=="fastq"){
                    string qscore = get_read_part(read_line_pos, sorted_read_file, format, reads_permutation[i], "qscore");
                    outfile << get_read_part(read_line_pos, sorted_read_file, format, reads_permutation[i], "separator")+"\n"+qscore << endl;
                }
            }
            else{
                outfile << get_read_part(read_line_pos, sorted_read_file, format, reads_permutation[i], "header")+"\n"+get_read_part(read_line_pos, sorted_read_file, format, reads_permutation[i], "sequence") << endl;
                if (format=="fastq"){
                    outfile << get_read_part(read_line_pos, sorted_read_file, format, reads_permutation[i], "separator")+"\n"+get_read_part(read_line_pos, sorted_read_file, format, reads_permutation[i], "qscore") << endl;
                }
            }
            i++;
        }
    }

}


int main(int argc, char *argv[])
{
    string sorted_read_file;
    string format;
    string reverse_order_file;
    string output_file;
    if (argc<3){
        cerr << "Usage : ./reverse_order [SORTED READS FILE] [FORMAT] [REVERSE ORDER FILE] [OUTPUT FILE]" << endl;
        return EXIT_FAILURE;
    }
    else{
        sorted_read_file = argv[1];
        format = argv[2];
        reverse_order_file = argv[3];
        output_file = argv[4];
    }

    reverse_file_order(sorted_read_file, reverse_order_file, output_file, format);

    return 0;
}