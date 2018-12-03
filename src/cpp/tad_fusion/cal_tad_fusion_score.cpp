//#: Date         :05/05/2017
//#: Author       :Linh Huynh
//#: Version      :1.0.0 

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

using namespace std;

typedef vector<int> IdList;
typedef vector<double> ValueList;
typedef vector<IdList> IdMatrix;
typedef vector<ValueList> ValueMatrix; 
typedef vector<string> StringList;
typedef vector<StringList> StringMatrix;

typedef map<string,int> IdMap;

#define STR_EQ        0
#define ROW_DELIM     '\n'
#define FIELD_DELIM   '\t'
#define COMMA_DELIM   ','

// Auxiliary functions
string num_to_string (int);
double max(double, double);
double min(double, double);
StringList split(const string&s, char delim);
// Common data structure functions
void resize_and_fill(ValueList&, int size, int value);                  // size, value
void resize_and_fill(ValueMatrix&, int size, int value);                // Squre matrix only
void read_a_table_file (string filename, char row_delim_char, char field_delim_char, 
                        char comment_char, int skipping_header_line_num, 
                        const IdList& col_list, StringMatrix& str_mat); // index of 1st col = 0
void read_a_value_list (string filename, int col, ValueList& val_list); // index of 1st col = 0
void write_a_value_list (string filename, ValueList&);
void sort_by_value (const ValueList& l, ValueList& sorted_index, ValueList& rank, bool is_asc);

int main (int argc, char** argv) {
  string del_filename, model_dir, out_filename;
  int window_size;
  double delta;

  // Parse the parameters
  if (argc % 2 == 0) {
    cout << "ERROR: Number of parameters must be even" << endl;
    return -1;
  }
  for (int i = 1; i < argc; i = i + 2) {
    string option(argv[i]);
    //cout << option << "\t" << argv[i+1] << endl;
    if (option.compare("-f") == STR_EQ)
      del_filename = argv[i+1];
    else if (option.compare("-md") == STR_EQ)
      model_dir = argv[i+1];
    else if (option.compare("-w") == STR_EQ)
      window_size = atoi(argv[i+1]);
    else if (option.compare("-d") == STR_EQ)
      delta = atof(argv[i+1]);
    else if (option.compare("-o") == STR_EQ)
      out_filename = argv[i+1];
    else {
      cout << "ERROR: Can not recognize the option " << option << endl;
      return -1;
    }
  }
  // Read the deletion file
  StringMatrix del_str_mat;
  read_a_table_file (del_filename, ROW_DELIM, FIELD_DELIM, '#', 0, {0,1,2}, del_str_mat);
  ValueList tad_fusion_score_list(del_str_mat.size(), 0);
  //for (int chr = 1; chr <= 23; chr++) {
  for (int chr = 1; chr <= 23; chr++) {
    string chr_str = "chr" + ((chr < 23)? num_to_string(chr):"X");
    cout << chr_str << endl;
    StringMatrix para_str_mat;
    read_a_table_file (model_dir + "/" + chr_str + ".model", ROW_DELIM, FIELD_DELIM,'#', 0, {0,1,2,3}, para_str_mat);
    double resolution = atof(para_str_mat[0][0].c_str());
    int chr_length = para_str_mat.size() - 1;  // skip the first line
    ValueList alpha(chr_length, 0), 
              beta_1(chr_length, 0), 
              beta_2(chr_length, 0),
              insulation(chr_length, 0); 
    for (int para = 0; para < chr_length; para++) {
      alpha[para] = atof(para_str_mat[para + 1][0].c_str());
      beta_1[para] = atof(para_str_mat[para + 1][1].c_str());
      beta_2[para] = atof(para_str_mat[para + 1][2].c_str());
      insulation[para] = atof(para_str_mat[para + 1][3].c_str());
    }
    for (int del = 0; del < del_str_mat.size(); del++) {
      if (del_str_mat[del][0].compare(chr_str) == STR_EQ) {
        int start_bin = round(atoi(del_str_mat[del][1].c_str())/resolution),
            end_bin = round(atoi(del_str_mat[del][2].c_str())/resolution),
            start_window = ((start_bin > window_size)? (start_bin - window_size):0),
            end_window = ((end_bin + window_size < (para_str_mat.size() - 1))? (end_bin + window_size):(para_str_mat.size() - 1));
        for (int i = start_window; i <= start_bin; i++) {
          for (int j = end_bin; j <= end_window; j++) {
            double all_insulation = 0, del_insulation = 0;
            for (int k = i; k <= j; k++) {
              all_insulation += insulation[k];
              if (k >= start_bin && k <= end_bin)
                del_insulation += insulation[k];
            }
            double before_del_val = exp(((beta_2[i] + beta_2[j])/2)*log(j - i) - all_insulation),
                   after_del_val = exp(((beta_2[i] + beta_2[j])/2)*log(j - i - (end_bin - start_bin)) - all_insulation + del_insulation);
            if (before_del_val < delta && after_del_val > delta)
              tad_fusion_score_list[del]++;
            //cout << i << "\t" << j << "\t" 
            //     << beta_2[i] << "\t" << beta_2[j] << "\t"
            //     << all_insulation << "\t" << del_insulation << "\t"
            //     << before_del_val << "\t" << after_del_val << endl;
          }
          //break;
        }
      }
    }
  }
  ofstream out_file(out_filename.c_str());
  for (int del = 0; del < del_str_mat.size(); del++) {
    if (del > 0)
      out_file << endl;
    out_file << del_str_mat[del][0] << "\t" 
             << del_str_mat[del][1] << "\t" 
             << del_str_mat[del][2] << "\t" 
             << tad_fusion_score_list[del];
  }
  out_file.close();
  cout << "Complete!" << endl;
}


void read_a_table_file (string filename, char row_delim_char, char field_delim_char, char comment_char, int skipping_header_line_num, const IdList& col_list, StringMatrix& str_mat) {
  str_mat.clear();
  ifstream in_file(filename.c_str());
  if (in_file.is_open()) {		
    cout << "Reading file " << filename << endl;
    int line_num = 0;
    for (string row; getline(in_file, row, row_delim_char); ) {
      line_num++;
      if (line_num > skipping_header_line_num) {
        istringstream ss(row);
        StringList record;
        for (string field; getline(ss, field, field_delim_char); )
        record.push_back(field);
        if (!record.empty() && !record[0].empty() && record[0][0] != comment_char) {
          StringList tmp;
          for (int i = 0; i < col_list.size(); i++)
            if (record.size() > col_list[i])
              tmp.push_back(record[col_list[i]]);
            else
              cout << "\tline " << line_num << " has only " << record.size() << " fields while the required index is " << col_list[i] << endl;
          if (tmp.size() == col_list.size())
              str_mat.push_back(tmp);
        }
      }
    }
    cout << "\t#lines " << line_num << "\t#records " << str_mat.size() << endl;
    in_file.close();	
  }
  else {
    cout << "Can not open file " << filename << endl;
    return;
  }
}

void read_a_value_list (string filename, int col, ValueList& val_list) {
  IdList col_list(1, col);
  StringMatrix str_mat;
  read_a_table_file (filename, ROW_DELIM, FIELD_DELIM, '#', false, col_list, str_mat);
  val_list.clear();
  for (int i = 0; i < str_mat.size(); i++)
    val_list.push_back(strtod(str_mat[i][0].c_str(), NULL));
}

void write_a_value_list (string filename, ValueList& val_list) {
  ofstream out_file(filename.c_str());
  for (int i = 0; i < val_list.size(); i++)
    out_file << i << "\t" << val_list[i] << endl;
  out_file.close();
}

void sort_by_value (const ValueList& l, ValueList& sorted_index, ValueList& rank, bool is_asc) {
  ValueList x;
  for (int i = 0; i < l.size(); i++)
    x.push_back(l[i]*(is_asc? 1:-1));
  int n = x.size();
  sorted_index.resize(n);
  std::size_t tmp(0);
  std::generate(std::begin(sorted_index), std::end(sorted_index), [&]{ return tmp++; });
  std::sort(std::begin(sorted_index), std::end(sorted_index), [&](int i, int j) { // Sort min -> max
    return (x[i] < x[j]);
  });
  rank.resize(n);
  for (int i = 0; i < n; i++)
    rank[sorted_index[i]] = i + 1;
}

StringList split(const string&s, char delim) {
  stringstream ss(s);
  string item;
  StringList elems;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

string num_to_string (int n) {
  stringstream ss;
  ss << n;
  return ss.str();
}

double max(double a, double b) {
  return (a > b)? a:b;
}

double min(double a, double b) {
  return (a < b)? a:b;
}

void resize_and_fill(ValueList& x, int size, int val) {
  if (x.size() != size)
    x.resize(size);
  for (int i = 0; i < size; i++)
    x[i] = val;
}
void resize_and_fill(ValueMatrix& x, int size, int val) {
  if (x.size() != size)
    x.resize(size);
  for (int i = 0; i < size; i++)
    resize_and_fill(x[i], size, val);
}

