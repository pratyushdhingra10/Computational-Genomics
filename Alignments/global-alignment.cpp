#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <limits>
#include <algorithm>

using namespace std;

struct DP_cell{
    double Sscore;
    double Dscore;
    double Iscore;
};

// Function to trim the spaces at the beginning and the end of the line.
string trim(const string& str)
{
    size_t first = str.find_first_not_of(' ');
    
    if (string::npos == first)
    {
        return str;
    }

    size_t last = str.find_last_not_of(' ');
    
    return str.substr(first, (last - first + 1));
}

// Reads the sequence input from the given filename and returns the sequences as a vector.
vector<string> read_input_file(string filename){
    
    fstream sequence_file;

    sequence_file.open(filename,ios::in);

    if(!sequence_file.is_open()){
        cout << "Something went wrong. The file couldnot open.";
        exit(1);
    }

    string line{};
    vector<string> sequences;
    
    bool new_sequence;
    string sequence_segment{};

    int blank_lines = 0;

    while(blank_lines < 2){
        if(!getline(sequence_file, line )){
            blank_lines+=1;

            continue;
        }

        blank_lines = 0;

        if (!line.empty() && line[0] == '>' && !sequence_segment.empty() ){
            sequences.push_back(sequence_segment);
            sequence_segment = "";
        }

        if (!line.empty() && line[0] != '>'){
            sequence_segment += trim(line);
        }

    };

    if(!sequence_segment.empty()){
        sequences.push_back(sequence_segment);
        sequence_segment = "";
    }

    return sequences;
};

class GlobalAlignment{
    private:
        string sequence_1;
        string sequence_2;

        int match_score;

        signed int gap_score;
        signed int mismatch_score;
        signed int affinity_score;

        vector<vector<DP_cell>> table;

        // Gives the maximum value for the provided three values. Used to find max while filling the 
        // Dp table.
        double getMaxValue(double value_1 , double value_2 ,double value_3){
            return max(max(value_1,value_2),value_3);
        }

        /* Returns the type with the corresponding max values. Used to find which is the next case i.e D, I ,S while
        retracing. */
        char getMaxValueType(vector<double> comparables, vector<char> types){
            int max_index = max_element(comparables.begin(), comparables.end()) - comparables.begin();
    
            return types[max_index];
        }

        /* Function to calculate the S score while filling the table. */
        double getSscore(vector<vector<DP_cell>>table, int index_i , int index_j){
            double prev_Sscore = table[index_i - 1][index_j -1].Sscore;
            double prev_Dscore = table[index_i - 1][index_j -1].Dscore;
            double prev_Iscore = table[index_i - 1][index_j -1].Iscore;

            double score = getMaxValue(prev_Dscore,prev_Iscore,prev_Sscore);

            if(sequence_1.at(index_i -1) == sequence_2.at(index_j-1)){
                return score + this->match_score;
            }

            return score + this->mismatch_score;
        }

        /* Function to find the previous best type if the current best type is S. Used for retracing. */
        char getSscorePreviousType(vector<vector<DP_cell>>table, int index_i , int index_j){
            double prev_Sscore = table[index_i - 1][index_j -1].Sscore;
            double prev_Dscore = table[index_i - 1][index_j -1].Dscore;
            double prev_Iscore = table[index_i - 1][index_j -1].Iscore;

            vector<char> types{'d','i','s'};
            vector<double> values{prev_Dscore, prev_Iscore, prev_Sscore};

            char type = getMaxValueType(values,types);

            return type;
        }

        /* Function to calculate the D score while filling the table. */
        double getDscore(vector<vector<DP_cell>>table, int index_i , int index_j){
            double prev_Sscore = table[index_i - 1][index_j].Sscore;
            double prev_Dscore = table[index_i - 1][index_j].Dscore;
            double prev_Iscore = table[index_i - 1][index_j].Iscore;

            double score = getMaxValue(prev_Dscore + this->gap_score,prev_Iscore + this->affinity_score + this->gap_score,prev_Sscore + this->affinity_score + this->gap_score);

            return score;
        }

        /* Function to find the previous best type if the current best type is D. Used for retracing. */        
        char getDscorePreviousType(vector<vector<DP_cell>>table, int index_i , int index_j){
            double prev_Sscore = table[index_i - 1][index_j].Sscore;
            double prev_Dscore = table[index_i - 1][index_j].Dscore;
            double prev_Iscore = table[index_i - 1][index_j].Iscore;

            vector<char> types{'d','i','s'};
            vector<double> values{prev_Dscore + this->gap_score, prev_Iscore + this->affinity_score + this->gap_score, prev_Sscore + this->affinity_score + this->gap_score};


            char type = getMaxValueType(values, types);

            return type;
        }
        
        /* Function to calculate the I score while filling the table. */
        double getIscore(vector<vector<DP_cell>>table, int index_i , int index_j){
            double prev_Sscore = table[index_i][index_j - 1].Sscore;
            double prev_Dscore = table[index_i][index_j - 1].Dscore;
            double prev_Iscore = table[index_i][index_j - 1].Iscore;


            double score = getMaxValue(prev_Iscore + this->gap_score, prev_Dscore + this->affinity_score + this->gap_score, prev_Sscore + this->affinity_score + this->gap_score);

            return score;
        }

        /* Function to find the previous best type if the current best type is I. Used for retracing. */   
        char getIscorePreviousType(vector<vector<DP_cell>>table, int index_i , int index_j){
            double prev_Sscore = table[index_i][index_j - 1].Sscore;
            double prev_Dscore = table[index_i][index_j - 1].Dscore;
            double prev_Iscore = table[index_i][index_j - 1].Iscore;

            vector<char> types{'d','i','s'};
            vector<double> values{prev_Dscore + this->affinity_score + this->gap_score,prev_Iscore + this->gap_score ,prev_Sscore + this->affinity_score + this->gap_score};
            
            char type = getMaxValueType(values,types);

            return type;
        }

        /* Function gives the best previous case for the retracing path. */
        char getRetraceCaseType(vector<vector<DP_cell>> table , int index_i , int index_j ,char current_case_type){
            char selected_type;

            if(current_case_type == 'd'){
                return getDscorePreviousType(table,index_i,index_j);
            }
            if(current_case_type == 'i'){
                return getIscorePreviousType(table,index_i,index_j);
            }

            return getSscorePreviousType(table,index_i,index_j);
        }

        /* Function that performs the retrace tasks. It is used a recursive fucntion. */
        void retrace(vector<vector<DP_cell>> table, int index_i , int index_j, vector<char> *aligned_str_1 , vector<char> *aligned_str_2, char cell_max_type){
            if(table[index_i][index_j].Dscore == 0 && table[index_i][index_j].Sscore == 0 && table[index_i][index_j].Iscore == 0 ){
                return;
            }

            char next_case_type;
            
            if(index_i == this->sequence_1.length() && index_j == this->sequence_2.length()){
                vector<char> types{'d','i','s'};
                vector<double> values{table[index_i][index_j].Dscore,table[index_i][index_j].Iscore,table[index_i][index_j].Sscore};

                cell_max_type = getMaxValueType(values,types);
            }

            next_case_type = getRetraceCaseType(table ,index_i, index_j , cell_max_type);
            
            if(cell_max_type == 'd'){
                
                /* Using index_i - 1 here because sequence starts at the 1 index in the table rather than 0. So to fetch the
                correct character from string we need to subtract the current index by 1. */
                aligned_str_1->push_back(this->sequence_1.at(index_i-1));
                aligned_str_2->push_back('-');

                return retrace(table,index_i -1 ,index_j, aligned_str_1,aligned_str_2,next_case_type);

            }else if(cell_max_type == 'i'){

                aligned_str_1->push_back('-');
                aligned_str_2->push_back(this->sequence_2.at(index_j-1));

                return retrace(table,index_i ,index_j - 1, aligned_str_1,aligned_str_2,next_case_type);
            }else{

                aligned_str_1->push_back(this->sequence_1.at(index_i-1));
                aligned_str_2->push_back(this->sequence_2.at(index_j-1));

                return retrace(table,index_i - 1,index_j - 1, aligned_str_1,aligned_str_2,next_case_type);
            }

            return;
            
        }

    public:
        GlobalAlignment(string sequence_1 , string sequence_2, int match_score , signed int mismatch_score , signed int gap_score, signed int affinity_score ){
            this->sequence_1 = sequence_1;
            this->sequence_2 = sequence_2;
            this->match_score = match_score;
            this->mismatch_score = mismatch_score;
            this->gap_score = gap_score;
            this->affinity_score = affinity_score;
            vector<vector<DP_cell>> table(sequence_1.length() +1, vector<DP_cell>(sequence_2.length() +1));

            this->table = table;
        }

        /* Function that fills the DP table in the feed forward method. */
        vector<vector<DP_cell>> fillAndGetTable(){
            
            double infinity = numeric_limits<double>::infinity();

            this->table[0][0].Sscore = 0;
            this->table[0][0].Dscore = 0;
            this->table[0][0].Iscore = 0;

            for (int index_i = 1 ; index_i <= this->sequence_1.length() ;  index_i++ ){
                this->table[index_i][0].Sscore = infinity * -1;
                this->table[index_i][0].Iscore = infinity * -1;
                this->table[index_i][0].Dscore =  this->affinity_score + index_i * this->gap_score;
            }

            for (int index_j = 1 ; index_j <= this->sequence_2.length() ;  index_j++ ){
                this->table[0][index_j].Sscore = infinity * -1;
                this->table[0][index_j].Iscore = this->affinity_score + index_j * this->gap_score;
                this->table[0][index_j].Dscore =  infinity * -1;
            }

            for (int index_1 = 2; index_1 <= this->sequence_1.length(); index_1++){
                for(int index_2 =2 ; index_2 <= this->sequence_2.length(); index_2++){
                    this->table[index_1][index_2].Sscore = getSscore(table, index_1,index_2);
                    this->table[index_1][index_2].Dscore = getDscore(table, index_1,index_2);
                    this->table[index_1][index_2].Iscore = getIscore(table, index_1,index_2);
                }
            }
            cout << this->affinity_score << "\n";
            return table;
        }

        /* Retraces the alignment from the end of the table and displays the output in the desired format. */
        void retraceAndGetAlignment(){

            int end_index_i = this->sequence_1.length();
            int end_index_j = this->sequence_2.length();

            vector<char> aligned_string_1;
            vector<char> aligned_string_2;

            char type;

            retrace(this->table , end_index_i, end_index_j, &aligned_string_1, &aligned_string_2, type);

            vector<char> match_signs;

            for(int i = aligned_string_1.size() - 1; i>= 0 ; i-- ){
                if(aligned_string_1[i] ==aligned_string_2[i]){
                    match_signs.push_back('|');
                }else{
                    match_signs.push_back(' ');
                }
            }

            for(int i = aligned_string_1.size() - 1; i>= 0 ; i-- ){
                cout << aligned_string_1[i];
            }

            cout <<endl;

            for(int i = 0; i < match_signs.size() ; i++ ){
                cout << match_signs[i];
            }

            cout <<endl;

            for(int i = aligned_string_2.size() - 1; i>= 0 ; i-- ){
                cout << aligned_string_2[i];
            }

        }
};

int main()
{
    vector<string> sequences = read_input_file("sequence.fasta");

    string sequence_1 = sequences[0];
    string sequence_2 = sequences[1];

    GlobalAlignment alignment_g(sequence_1,sequence_2, 1 , -2, -1,-5);

    vector<vector<DP_cell>> alignment_table = alignment_g.fillAndGetTable();

    alignment_g.retraceAndGetAlignment();

    return 0;
}
