//#: Date         :05/05/2017
//#: Author       :Linh Huynh
//#: Version      :1.0.0 

#include <ilcplex/cplex.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <math.h>

using namespace std;

typedef vector<int> IdList;
typedef vector<double> ValueList;
typedef vector<IdList> IdMatrix;
typedef vector<ValueList> ValueMatrix; 
typedef vector<string> StringList;

typedef map<string,int> IdMap;

#define STR_EQ 0
#define ROW_DELIM '\n'
#define FIELD_DELIM '\t'
#define COMMA_DELIM ','

#define CHR_NUM			23
#define RESOLUTION		5000

#define ZERO_ENTRY		1e-1

#define ARROW_HEAD		1
#define CATCH			2
#define HICSEG			3
#define INSULATION_SCORE	4
#define ARMATUS			5

#define FLOYD_SHORTEST_PATH	1
#define APPROX_SHORTEST_PATH	0

#define MAX_LENGTH_OBJ			100
#define SHORTEST_PATH_WINDOW_SIZE	50

struct Distribution {
	void add_value (double x) {
		val_list.push_back(x);
		if (val_list.size() == 1)
			min_val = max_val = x;
		else {
			if (x < min_val)
				min_val = x;
			if (x > max_val)
				max_val = x;
		}
	}
	void add_a_value_list (const ValueList& val_list) {
		for (int i = 0; i < val_list.size(); i++)
			add_value(val_list[i]);
	}
	void make_dist (int bin_num) {
		freq.resize(bin_num);
		for (int i = 0; i < bin_num; i++)
			freq[i] = 0;
		if (val_list.empty()) {
			cout << "Error: Make distribution from an empty list!" << endl;
			return;
		}
		if (max_val == min_val) {
			mean = min_val;
			std = 0;
			return;
		}
		mean = 0;
		for (int i = 0; i < val_list.size(); i++) {
			int l = floor((val_list[i] - min_val)*bin_num/(max_val - min_val));
			if (l == bin_num)
				l--;
			freq[l] += 1.0/val_list.size();
			mean += val_list[i]/val_list.size();
		}
		std = 0;
		for (int i = 0; i < val_list.size(); i++) 
			std += (val_list[i] - mean)*(val_list[i] - mean);
		std = sqrt(std/val_list.size());
		//val_list.clear();
	}
	double get_prob (double x) const {
		double prob = 0;
		if (x >= min_val || x <= max_val) {
			if (min_val == max_val) {
				cout << "Warning: Distribution has only one value!" << min_val << endl;
				return 1;
			}
			int l =  floor((x - min_val)*freq.size()/(max_val - min_val));
			if (l == freq.size())
				l--;
			return freq[l];
		}
		return prob;
	}
	double get_cdf (double x) const {
		double count = 0;
		for (int i = 0; i < val_list.size(); i++)
			if (val_list[i] <= x)
				count++;
		return ((count > 0)? count/val_list.size(): 0);
	}
	void print (string filename) {
		ofstream out_file(filename.c_str());
		out_file << min_val - 0.5*(max_val - min_val)/freq.size() << "\t" << 0 << endl;
		for (int i = 0; i < freq.size(); i++)
			out_file << min_val + (0.5 + i)*(max_val - min_val)/freq.size() << "\t" << freq[i]*100 << endl;
		out_file << max_val + 0.5*(max_val - min_val)/freq.size() << "\t" << 0 << endl;
		out_file.close();
	}
	double min_val, max_val, mean, std;
	ValueList freq, val_list;
};

struct DistributionList {
	void reset (int max_length) {
		dist_list.resize(max_length + 1);
	}
	void add_value (int length, double x) {
		if (length >= dist_list.size())
			cout << "Error: length = " << length << " while max length = " << dist_list.size() - 1 << endl;
		else
			dist_list[length].add_value(x);
	}
	void make_dist(int bin_num) {
		for (int i = 0; i < dist_list.size(); i++) {
			dist_list[i].make_dist(bin_num);
		}
	}

	double get_prob (int length, double x) const {
		double prob = 0;
		if (length >= dist_list.size()) 
			cout << "Error: length = " << length << " while max length = " << dist_list.size() - 1 << endl;
		else 
			prob = dist_list[length].get_prob(x);
		return prob;
	}
	vector<Distribution> dist_list;
};

typedef vector<vector<Distribution> > DistributionMatrix;

struct DistributionSummary {
	DistributionSummary() {
		size = 0;
		mean = std = 0;
	}
	void add_value (double x) {
		double new_mean = (size*mean + x)/(size + 1);
		double tmp = (x*x + size*(std*std + mean*mean))/(size+1) - new_mean*new_mean;
		if (tmp > 1e-10)
			std = sqrt(tmp);
		else
			std = 0;
		mean = new_mean;
		size++;
	}
	void merge (const DistributionSummary& s) {
		int new_size = size + s.size;
		double new_mean = (size*mean + s.size*s.mean)/new_size;
		std = sqrt((size*(mean*mean + std*std) + s.size*(s.mean*s.mean + s.std*s.std))/new_size - new_mean*new_mean);
		mean = new_mean;
		size += s.size;
	}
	int size;
	double mean, std;
};

typedef vector<DistributionSummary> DistributionSummaryList;
typedef vector<DistributionSummaryList> DistributionSummaryMatrix;

struct Deletion {
	Deletion () {
		left = right = 0;
	}
	Deletion (int l, int r) {
		left = l;
		right = r;
	}	
	int left, right;
};
typedef vector<Deletion> DeletionList;
typedef vector<DeletionList> DeletionMatrix;

struct TAD {
	TAD () {
		left = right = 0;
	}
	TAD (int l, int r) {
		left = l;
		right = r;
	}
	int left, right;
};
typedef vector<TAD> TADList;
typedef vector<TADList> TADMatrix;

struct TADBoundary {
	TADBoundary () {
		left = right = score = -1;
	}
	TADBoundary (int l, int r, double s) {
		left = l;
		right = r;
		score = s;
	}
	int left, right;
	double score;
};
typedef vector<TADBoundary> TADBoundaryList;
typedef vector<TADBoundaryList> TADBoundaryMatrix;

struct DeletionInfo {
	DeletionInfo(int c, int d) {
		chr = c;
		del = d;
	};
	DeletionInfo() {
		chr = del = -1;
	};
	int chr;
	int del;
};
typedef vector<DeletionInfo> DeletionInfoList;

struct HiCDataModel {
	int get_right() const {
		return left + alpha.size() - 1;
	}
	ValueList alpha;
	double beta;
	ValueList resistance;
	int left;
};
typedef vector<HiCDataModel> HiCDataModelList;

void write_hi_c_data_model_list(string filename, const HiCDataModelList& hi_c_data_model_list) {
	ofstream out_file(filename.c_str());
	for (int i = 0; i < hi_c_data_model_list.size(); i++) {
		out_file << "#" << endl;
		out_file << hi_c_data_model_list[i].left << "\t" << hi_c_data_model_list[i].beta << endl;
		for (int k = 0; k < hi_c_data_model_list[i].alpha.size(); k++)
			out_file << hi_c_data_model_list[i].alpha[k] << "\t" << hi_c_data_model_list[i].resistance[k] << endl;
	}
	out_file.close();
}

void read_hi_c_data_model_list(string filename, HiCDataModelList& hi_c_data_model_list) {
	ifstream in_file(filename.c_str());
	if (in_file.is_open()) {		
		int line_num = 0, status = -1;
		for (string row; getline(in_file, row, ROW_DELIM); ) {
			line_num++;
			istringstream ss(row);
			StringList record;
			for (string field; getline(ss, field, FIELD_DELIM); )
				record.push_back(field);
			if (!record.empty() && !record[0].empty()) {
				if (record[0][0] == '#') {
					status = 0;
					HiCDataModel tmp;
					hi_c_data_model_list.push_back(tmp);
				}
				else if (record.size() == 2) {
					int l = hi_c_data_model_list.size() - 1;
					if (status == 0) {
						hi_c_data_model_list[l].left = strtod(record[0].c_str(), NULL);
						hi_c_data_model_list[l].beta = strtod(record[1].c_str(), NULL);
						status = 1;
					}
					else {
						hi_c_data_model_list[l].alpha.push_back(strtod(record[0].c_str(), NULL));
						hi_c_data_model_list[l].resistance.push_back(strtod(record[1].c_str(), NULL));
					}
				}
				else 
					cout << "ERROR: Can not read the HiC data model file at line " << line_num << ", there are " << record.size() << " != 2 lines!" << endl;			
			}
		}
		cout << "File " << filename << ", #lines " << line_num << endl;
		in_file.close();	
	}
	else {
		cout << "Can not open file " << filename << endl;
		return;
	}
}

string num_to_string (int n) {
	stringstream ss;
	ss << n;
	return ss.str();
}

float max(double a, double b) {
	return (a > b)? a:b;
}

float min(double a, double b) {
        return (a < b)? a:b;
}

void resize_and_fill(ValueList& x, int n, int a) {
	if (x.size() != n)
		x.resize(n);
	for (int i = 0; i < n; i++)
		x[i] = a;
}

void read_a_value_list (string filename, ValueList& val_list) {
	val_list.clear();
	ifstream in_file(filename.c_str());
	if (in_file.is_open()) {		
		int line_num = 0;
		for (string row; getline(in_file, row, ROW_DELIM); ) {
			line_num++;
			istringstream ss(row);
			StringList record;
			for (string field; getline(ss, field, FIELD_DELIM); )
				record.push_back(field);
			if (!record.empty() && !record[0].empty() && record[0][0] != '#') {
				if (record.size() >= 1) 
					val_list.push_back(strtod(record[0].c_str(), NULL));
			}
		}
		cout << "File " << filename << ", #lines " << line_num << endl;
		in_file.close();	
	}
	else {
		cout << "Can not open file " << filename << endl;
		return;
	}
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

void find_TAD_boundary_insulation_index (const ValueMatrix& mat, int resolution, int window_size, ValueList& index_list) {
	ValueList tmp(mat.size(), 0), cell_num(mat.size(), 0);
	tmp[0] = 0;
	cell_num[0] = 0;
	index_list.clear();
	for (int i = 1; i < mat.size(); i++) {
		int col_num = 0, row_num = 0;
		double col = 0, row = 0;
		for (int j = i - 2 - window_size; j <= i - 2; j++) {
			if (j >= 0) {
				col += mat[j][i];
				col_num++;
			}
		}
		for (int j = i + 1; j <= i + 1 + window_size; j++) {
			if (j < mat.size()) {
				row += mat[i-1][j];
				row_num++;
			}
		}
		tmp[i] = tmp[i-1] - col + row;
		cell_num[i] = cell_num[i-1] - col_num + row_num;
		index_list.push_back((cell_num[i] > 0)?tmp[i]/cell_num[i]:0);
	}	
}

void read_deletion_mutation_file (string filename, int length_threshold, DeletionMatrix& del_matrix) {
	IdMap chr_name_to_id;
	del_matrix.resize(23);
	for (int chr = 1; chr <= 23; chr++) {
		string chr_name = "chr" + ((chr < 23)? num_to_string(chr):"X");
		chr_name_to_id[chr_name] = chr - 1;
		del_matrix[chr - 1].clear();
	}

	int del_num = 0;
	ifstream in_file(filename.c_str());
	if (in_file.is_open()) {		
		int line_num = 0;
		for (string row; getline(in_file, row, ROW_DELIM); ) {
			line_num++;
			istringstream ss(row);
			StringList record;
			for (string field; getline(ss, field, FIELD_DELIM); )
				record.push_back(field);
			if (!record.empty() && !record[0].empty() && record[0][0] != '#') {
				if (record.size() >= 3) {
					if (chr_name_to_id.find(record[0]) != chr_name_to_id.end()) {
						Deletion del_tmp(atoi(record[1].c_str()), atoi(record[2].c_str()));
						if (del_tmp.right - del_tmp.left >= length_threshold) 
							del_matrix[chr_name_to_id[record[0]]].push_back(del_tmp);						
					}
					else
						cout << "ERROR: Unknown chromosome of the deletion " << record[0] << "\t" << record[1] << "\t" << record[2] << endl;
					//cout << record[0] << "\t" << record[1] << "\t" << record[2] << endl;
					del_num++;								
				}
			}
		}
		cout << "File " << filename << ", #lines " << line_num << endl;
		cout << "#deletion " << del_num << endl;
		in_file.close();		
	}
	else {
		cout << "Can not open file " << filename << endl;
		return;
	}
}

void read_CTCF (string filename, DeletionMatrix& CTCF_matrix) {
	IdMap chr_name_to_id;
	CTCF_matrix.resize(23);
	for (int chr = 1; chr <= 23; chr++) {
		string chr_name = "chr" + ((chr < 23)? num_to_string(chr):"X");
		chr_name_to_id[chr_name] = chr - 1;
		CTCF_matrix[chr - 1].clear();
	}

	int site_num = 0;
	ifstream in_file(filename.c_str());
	if (in_file.is_open()) {		
		int line_num = 0;
		for (string row; getline(in_file, row, ROW_DELIM); ) {
			line_num++;
			istringstream ss(row);
			StringList record;
			for (string field; getline(ss, field, FIELD_DELIM); )
				record.push_back(field);
			if (!record.empty() && !record[0].empty() && record[0][0] != '#') {
				if (record.size() >= 3) {
					if (chr_name_to_id.find(record[0]) != chr_name_to_id.end()) {
						Deletion del_tmp(atoi(record[1].c_str()), atoi(record[2].c_str()));
						if (del_tmp.right >= del_tmp.left)
							CTCF_matrix[chr_name_to_id[record[0]]].push_back(del_tmp);
					}
					else
						cout << "ERROR: Unknown chromosome of the CTCF site " << record[0] << "\t" << record[1] << "\t" << record[2] << endl;
					site_num++;							
				}
			}
		}
		cout << "File " << filename << ", #lines " << line_num << endl;
		cout << "#site " << site_num << endl;
		in_file.close();
	}
	else {
		cout << "Can not open file " << filename << endl;
		return;
	}
}

void read_hi_c_data (string filename, int resolution, ValueMatrix& hi_c_matrix) {
	IdList src, dest;
	ValueList val_list;
	int pos_max = -1;

	ifstream in_file(filename.c_str());
	if (in_file.is_open()) {		
		int line_num = 0;
		for (string row; getline(in_file, row, ROW_DELIM); ) {
			line_num++;
			istringstream ss(row);
			StringList record;
			for (string field; getline(ss, field, FIELD_DELIM); )
				record.push_back(field);
			if (!record.empty() && !record[0].empty() && record[0][0] != '#') {
				if (record.size() >= 3) {
					int x = atoi(record[0].c_str()), y = atoi(record[1].c_str());
					if (pos_max < x)
						pos_max = x;
					if (pos_max < y)
						pos_max = y;
					src.push_back(x);
					dest.push_back(y);
					val_list.push_back(atof(record[2].c_str()));				
				}
			}
		}
		cout << "File " << filename << ", #lines " << line_num << endl;
		in_file.close();
	}
	else {
		cout << "Can not open file " << filename << endl;
		return;
	}

	int size  = floor(pos_max/resolution) + 1;
	cout << "Matrix size " << size << endl;
	hi_c_matrix.resize(size);
	for (int i = 0; i < size; i++) {
		hi_c_matrix[i].resize(size);
		for (int j = 0; j < size; j++)
			hi_c_matrix[i][j] = 0;
	}
	for (int i = 0; i < src.size(); i++) {
		int x = floor(src[i]/resolution);
		int y = floor(dest[i]/resolution);
		if (hi_c_matrix[x][y] > 0 || hi_c_matrix[y][x] > 0)
			cout << "ERROR: matrix at " << x << "\t" << y << " has value " << hi_c_matrix[x][y] << endl;
		else {			
			hi_c_matrix[x][y] = val_list[i];
			hi_c_matrix[y][x] = val_list[i];
		}
	}
}

void read_TAD_data_with_Juice_format (const string& filename, TADMatrix& tad_mat, int min_size) {
	IdMap chr_name_to_id;
	tad_mat.resize(23);
	for (int chr = 1; chr <= 23; chr++) {
		string chr_name = ((chr < 23)? num_to_string(chr):"X");
		chr_name_to_id[chr_name] = chr - 1;
		tad_mat[chr - 1].clear();
	}
	
	ifstream in_file(filename.c_str());
	if (in_file.is_open()) {
		ValueList left_list, right_list;
		int line_num = 0;
		for (string row; getline(in_file, row, ROW_DELIM); ) {
			line_num++;
			istringstream ss(row);
			StringList record;
			for (string field; getline(ss, field, FIELD_DELIM); )
				record.push_back(field);
			if (!record.empty() && !record[0].empty() && record[0][0] != '#') {
				if (record.size() >= 4) {
					if (record[0].compare(record[3]) == STR_EQ) { // the first header line will be skipped by this
						if (chr_name_to_id.find(record[0]) != chr_name_to_id.end()) {
							TAD tad_tmp (atoi(record[1].c_str()), atoi(record[2].c_str()));
							if (tad_tmp.right - tad_tmp.left >= min_size)
								tad_mat[chr_name_to_id[record[0]]].push_back(tad_tmp);
						}
						else
							cout << "ERROR: Unknown chromosome of the TAD " << record[0] << "\t" << record[1] << "\t" << record[2] << endl;
					}
				}
			}
		}
		cout << "File " << filename << ", #lines " << line_num << endl;
		in_file.close();
	}
	else
		cout << "Can not open file " << filename << endl;
}

void read_TAD_boundary_with_4_col_format (const string& filename, TADBoundaryMatrix& tad_boundary_mat) {
	IdMap chr_name_to_id;
	tad_boundary_mat.resize(23);
	for (int chr = 1; chr <= 23; chr++) {
		string chr_name = ((chr < 23)? num_to_string(chr):"X");
		chr_name_to_id[chr_name] = chr - 1;
		tad_boundary_mat[chr - 1].clear();
	}

	ifstream in_file(filename.c_str());
	if (in_file.is_open()) {
		ValueList left_list, right_list;
		int line_num = 0;
		for (string row; getline(in_file, row, ROW_DELIM); ) {
			line_num++;
			istringstream ss(row);
			StringList record;
			for (string field; getline(ss, field, FIELD_DELIM); )
				record.push_back(field);
			if (!record.empty() && !record[0].empty() && record[0][0] != '#') {
				if (record.size() >= 4) {
					if (chr_name_to_id.find(record[0]) != chr_name_to_id.end()) {
						TADBoundary tad_boundary_tmp (atoi(record[1].c_str()), atoi(record[2].c_str()), atof(record[3].c_str()));
						tad_boundary_mat[chr_name_to_id[record[0]]].push_back(tad_boundary_tmp);
					}
					else
						cout << "ERROR: Unknown chromosome of the TAD boundary " << record[0] << "\t" << record[1] << "\t" << record[2] << "\t" << record[3] << endl;
				}
			}
		}
		cout << "File " << filename << ", #lines " << line_num << endl;
		in_file.close();
	}
	else
		cout << "Can not open file " << filename << endl;
}

void print_tad_list (string filename, const TADList& tad_list, int chr) {
	ofstream out_file (filename.c_str());
	out_file << "chr1\tx1\tx2\tchr2\ty1\ty2\tcolor" << endl;
	for (int i = 0; i < tad_list.size(); i++)
		out_file << chr << "\t" << tad_list[i].left << "\t" << tad_list[i].right << "\t" << chr << "\t" << tad_list[i].left << "\t" << tad_list[i].right << "\t0,0,225" << endl; 
	out_file.close();
}

void print_boundary_list (string filename, int chr_num, const TADBoundaryList& boundary_list) {
	ofstream boundary_file (filename.c_str());
	boundary_file << "variableStep chrom=chr" + num_to_string(chr_num) << endl;
	for (int i = 0; i < boundary_list.size(); i++) {
		boundary_file << round((boundary_list[i].left + boundary_list[i].right)/2) << "\t" << boundary_list[i].score << endl;
	}
	boundary_file.close();
}

void print_wig_file (string filename, int chr_num, int left, int right, int resolution, const ValueList& score_list) {
	if (right - left + 1 < score_list.size())
		cout << "Error: Exoprt wig file error!" << endl;
	ofstream out_file (filename.c_str());
	out_file << "variableStep chrom=chr" + num_to_string(chr_num) << endl;
	for (int i = left; i <= right; i++)
		out_file << i*resolution << "\t" << score_list[i - left] << endl;
	out_file.close();
}

void print_data_file (string filename, int left, int right, int resolution, const ValueList& score_list, const IdList& left_TAD_list, const IdList& right_TAD_list) {
	if (right - left + 1 < score_list.size())
		cout << "Error: Exoprt wig file error!" << endl;
	ofstream out_file (filename.c_str());
	for (int i = left; i <= right; i++)
		out_file << i << "\t" << score_list[i - left] << "\t" << left_TAD_list[i - left] << "\t" << right_TAD_list[i- left] << "\t"
			<< i*resolution << "\t" << left_TAD_list[i - left]*resolution << "\t" << right_TAD_list[i- left]*resolution << endl;
	out_file.close();
}

void print_dist_list (const DistributionSummaryList& dsl) {
	for (int i = 0; i < dsl.size(); i++)
		cout << i << "\t" << dsl[i].mean << "\t" << dsl[i].std << endl;
}

double RSQ (const ValueList& x, const ValueList& y) { // x: pred, y: obs
	if (x.size() != y.size() || x.empty() || y.empty()) {
		cout << "Error: RSQ of two different vector sizes " << x.size() << " and " << y.size() << endl;
		return 0;
	}
	double x2 = 0, y2 = 0, xy = 0, x_mean, y_mean;
	for (int i = 0; i < x.size(); i++) {
		x_mean += ((x[i] > ZERO_ENTRY)? x[i]:ZERO_ENTRY);
		y_mean += ((y[i] > ZERO_ENTRY)? y[i]:ZERO_ENTRY);
	}
	x_mean /= x.size();
	y_mean /= y.size();
	for (int i = 0; i < x.size(); i++) {
		double w = ((x[i] > ZERO_ENTRY)? x[i]:ZERO_ENTRY);
		double z = ((y[i] > ZERO_ENTRY)? y[i]:ZERO_ENTRY);
		x2 += (w - x_mean)*(w - x_mean);
		y2 += (z - y_mean)*(z - y_mean);
		xy += (w - x_mean)*(z - y_mean);
	}
	//cout << x2 << "\t" << y2 << "\t" << xy << endl;
	if (x2 < 1e-9 && y2 < 1e-9)
		return 1;
	else if (x2 < 1e-9 || y2 < 1e-9)
		return 0;
	else
		return xy/sqrt(x2*y2);		
}

double average_log_ratio (const ValueList& pred, const ValueList& obs) {
	double s = 0, count = 0;
	for (int i = 0; i < pred.size(); i++) {
		double w = ((pred[i] > ZERO_ENTRY)? pred[i]:ZERO_ENTRY);
		double z = ((obs[i] > ZERO_ENTRY)? obs[i]:ZERO_ENTRY);
		//if (obs[i] > ZERO_ENTRY) {
			s += abs(log(w) - log(z));
			count++;
		//}		
	}
	if (count > 0)		// for non-zero only
		s = s/count;
	return s;
}

double L1(const ValueList& pred, const ValueList& obs) {
	double s = 0;
	for (int i = 0; i < pred.size(); i++)
		s += abs(pred[i] - obs[i]);
	if (pred.size() > 0)
		s = s/pred.size();
	return s;
}

double L2(const ValueList& pred, const ValueList& obs) {
	double s = 0;
	for (int i = 0; i < pred.size(); i++)
		s += pow(pred[i] - obs[i], 2);
	if (pred.size() > 0)
		s = s/pred.size();
	return s;
}

void line_LQ (const ValueList& x, const ValueList& y, double& a, double& b) {	// y = ax + b
	double x_sum = 0, x2 = 0, y_sum = 0, xy = 0;
	for (int i = 0; i < x.size(); i++) {
		x_sum += x[i];
		x2 += x[i]*x[i];
		y_sum += y[i];
		xy += x[i]*y[i];
	}
	a = (x.size()*xy - x_sum*y_sum)/(x.size()*x2 - x_sum*x_sum);
	b = (y_sum*x2 - x_sum*xy)/(x.size()*x2 - x_sum*x_sum);
}

void print_prediction (string prefix, const ValueMatrix& hi_c_mat, int left_1, int right_1, int left_2, int right_2, const ValueMatrix& pred_mat) {
	string filename = prefix + "_tab.dat";
	ofstream tab_file (filename.c_str());
	filename = prefix + ".dat";
	ofstream col_file (filename.c_str());
	int l1 = right_1 - left_1 + 1, l2 = right_2 - left_2 + 1;
	for (int i = 0; i < l1 + l2; i++) {
		for (int j = 0; j < l1 + l2; j++) {
			if (i <= j) {
				if (i < l1 && j < l1)
					tab_file << hi_c_mat[left_1 + i][left_1 + j] << "\t";
				else if (i >= l1 && j >= l1)
					tab_file << hi_c_mat[left_2 + i - l1][left_2 + j - l1] << "\t";
				else {
					tab_file << "(" << hi_c_mat[left_1 + i][left_2 + j - l1] << "|" << round(100*pred_mat[i][j - l1])/100 << ")\t";
					col_file << i + left_1 << "\t" << j - l1 + left_2 << "\t" << hi_c_mat[left_1 + i][left_2 + j - l1] << "\t" << round(100*pred_mat[i][j - l1])/100 << endl; 
				}
			}
			else 
				tab_file << "-1\t";			
		}
		tab_file << endl;
	}
	tab_file.close();
	col_file.close();
}

ValueList naive_prediction (const ValueMatrix& hi_c_mat, int left_1, int right_1, int left_2, int right_2, bool print_out) {
	int max_length = (right_1 - left_1 > right_2 - left_2)? (right_1 - left_1):(right_2 - left_2);
	vector<DistributionSummary> dist_sum_list(max_length + 1);
	for (int length = 1; length <= right_1 - left_1; length++)
		for (int pos = left_1; pos + length <= right_1; pos++)
			dist_sum_list[length].add_value(hi_c_mat[pos][pos + length]);
	for (int length = 1; length <= right_2 - left_2; length++)
		for (int pos = left_2; pos + length <= right_2; pos++)
			dist_sum_list[length].add_value(hi_c_mat[pos][pos + length]);
	
	ValueList x,y;
	for (int i = 1; i <= max_length; i++) {
		x.push_back(log(i));
		y.push_back(log((dist_sum_list[i].mean > 0)? dist_sum_list[i].mean:ZERO_ENTRY));
	}
	double a,b;
	line_LQ(x, y, a, b);

	//cout << "(a,b)" << a << "\t" << b << endl;

	//pred_mat.resize(right_1 - left_1 + 1);
	ValueList obs, pred;
	for (int i = left_1; i <= right_1; i++) {
		for (int j = left_2; j <= right_2; j++) {
			double obs_val = hi_c_mat[i][j];
			double log_ij = log(right_1 - i + j - left_2 + 1);
			double pred_val = exp(a*log_ij + b);
			//pred_mat[i - left_1].push_back(pred_val);
			obs.push_back(obs_val);
			pred.push_back(pred_val);
		}
	}
	if (print_out) {
		string file_name = "Debug/Naive_" + num_to_string(left_1) + "_" + num_to_string(right_1) + ".dat";
		ofstream out_file(file_name.c_str());
		for (int i = 0; i < pred.size(); i++)
			out_file << pred[i] << "\t" << obs[i] << endl;
		out_file.close();
	}

	ValueList tmp;
	tmp.push_back(RSQ(pred, obs));
	tmp.push_back(average_log_ratio(pred, obs));
	tmp.push_back(L1(pred, obs));
	tmp.push_back(L2(pred, obs));
	return tmp;
}

ValueList norm_naive_prediction (const ValueMatrix& hi_c_mat, int left_1, int right_1, int left_2, int right_2, ValueMatrix& pred_mat) {
	int max_length = (right_1 - left_1 > right_2 - left_2)? (right_1 - left_1):(right_2 - left_2);
	ValueList left_norm (right_1 - left_1 + 1), right_norm(right_2 - left_2 + 1);
	for (int i = left_1; i <= right_1; i++)
		left_norm[i - left_1] = (hi_c_mat[i][i] > 0)? hi_c_mat[i][i]:ZERO_ENTRY;
	for (int i = left_2; i <= right_2; i++)
		right_norm[i- left_2] = (hi_c_mat[i][i] > 0)? hi_c_mat[i][i]:ZERO_ENTRY;
	vector<DistributionSummary> dist_sum_list(max_length + 1);
	for (int length = 1; length <= right_1 - left_1; length++)
		for (int pos = left_1; pos + length <= right_1; pos++)
			dist_sum_list[length].add_value(hi_c_mat[pos][pos + length]/sqrt(left_norm[pos - left_1]*left_norm[pos + length - left_1]));
	for (int length = 1; length <= right_2 - left_2; length++)
		for (int pos = left_2; pos + length <= right_2; pos++)
			dist_sum_list[length].add_value(hi_c_mat[pos][pos + length]/sqrt(right_norm[pos - left_2]*right_norm[pos + length - left_2]));
	
	ValueList x,y;
	for (int i = 1; i <= max_length; i++) {
		x.push_back(log(i));
		y.push_back(log((dist_sum_list[i].mean > 0)? dist_sum_list[i].mean:ZERO_ENTRY));
	}
	double a,b;
	line_LQ(x, y, a, b);

	//cout << "(a,b)" << a << "\t" << b << endl;

	pred_mat.resize(right_1 - left_1 + 1);
	ValueList obs, pred;
	for (int i = left_1; i <= right_1; i++) {
		for (int j = left_2; j <= right_2; j++) {
			double obs_val = hi_c_mat[i][j];
			double log_ij = log(right_1 - i + j - left_2 + 1);
			double pred_val = exp(a*log_ij + b)*sqrt(left_norm[i - left_1]*right_norm[j - left_2]);
			pred_mat[i - left_1].push_back(pred_val);
			obs.push_back(obs_val);
			pred.push_back(pred_val);
		}
	}
	ValueList tmp;
	tmp.push_back(RSQ(pred, obs));
	tmp.push_back(average_log_ratio(pred, obs));
	return tmp;
}

void Floyd_APSP (const ValueMatrix& graph_dist_mat, ValueMatrix& shortest_dist_mat) {
	int n = graph_dist_mat.size();
	ValueMatrix shortest_dist_mat_tmp = graph_dist_mat;
	ValueMatrix prev_node_mat = graph_dist_mat;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			shortest_dist_mat_tmp[i][j] = graph_dist_mat[i][j];
			prev_node_mat[i][j] = -1;
		}
		shortest_dist_mat_tmp[i][i] = 0.0;
	}
	for (int k = 0; k < n; k++) {
		//cout << k << endl;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (shortest_dist_mat_tmp[i][k] + shortest_dist_mat_tmp[k][j] < shortest_dist_mat_tmp[i][j]) {
					shortest_dist_mat_tmp[i][j] = shortest_dist_mat_tmp[i][k] + shortest_dist_mat_tmp[k][j];
					//prev_node_mat[i][j] = k;
  				}
			}
		}
	}
	shortest_dist_mat.resize(n);
	for (int i = 0; i < n; i++) {
		shortest_dist_mat[i].resize(i);
		for (int j = 0; j < i; j++)
			shortest_dist_mat[i][j] = shortest_dist_mat_tmp[i][i - j - 1];	// shortest_dist_mat[i][j] = shortest path from i to i-j-1,
	}
}

void limited_shortest_path(const ValueMatrix& distance_mat, int min_length, int window, ValueMatrix& shortest_dist_mat) {
	int n = distance_mat.size();
	shortest_dist_mat.resize(n);	// shortest_dist_mat[i][j] = shortest path from i to i-j-1, j <= min_length, all path are in [i-j-w,i+w]
	ValueMatrix shortest_dist_mat_tmp(n);
	for (int k = 0; k < n; k++) {
		//if (k % 100 == 0)
		//	cout << k << endl;
		int l = (k >= min_length + 2*window)? (min_length + 2*window):k;
		resize_and_fill(shortest_dist_mat_tmp[k], l, 0);
		ValueList tmp(l, 0);
		for (int i = 1; i <= l; i++) {
			shortest_dist_mat_tmp[k][i - 1] = tmp[i - 1] = distance_mat[k][k - i];
			for (int j = 1; j <= l; j++)
				if (i != j) {
					//cout << k << "\t" << i << "\t" << j << endl;
					//cout << shortest_dist_mat[((j > i)?(k-i):(k-j))].size() << endl;
					double d_tmp = ((j > i)? shortest_dist_mat_tmp[k-i][j-i-1]:shortest_dist_mat_tmp[k-j][i-j-1]);
					if (tmp[i-1] > distance_mat[k][k-j] + d_tmp)
						tmp[i-1] = distance_mat[k][k-j] + d_tmp;
				}
		}
		for (int i = 1; i <= l; i++)
			for (int j = i + 1; j <= l; j++)
				if (shortest_dist_mat_tmp[k-i][j-i-1] > tmp[i-1] + tmp[j-1])
					shortest_dist_mat_tmp[k-i][j-i-1] = tmp[i-1] + tmp[j-1];
		IdList node_list;
		if (k >= window) 
			node_list.push_back(k - window);
		if (k == n - 1) {
			int i_min = ((k >= window - 1)? (k - window + 1):0);
			for (int i = i_min; i < n; i++)
				node_list.push_back(i);
		}	

		for (int i = 0; i < node_list.size(); i++) {
			int l_p = (node_list[i] >= min_length)? min_length : node_list[i];
			shortest_dist_mat[node_list[i]].resize(l_p);
			for (int j = 0; j < l_p; j++)
				shortest_dist_mat[node_list[i]][j] = shortest_dist_mat_tmp[node_list[i]][j];
		}
	}

	/*for (int i = 0; i < distance_mat.size(); i++) {
		for (int j = 0; j < distance_mat.size(); j++)
			cout << distance_mat[i][j] << " ";
		cout << endl;
	}
	for (int i = 0; i < shortest_dist_mat.size(); i++) {
		for (int j = 0; j < shortest_dist_mat[i].size(); j++)
			cout << shortest_dist_mat[i][j] << " ";
		cout << endl;
	}*/		
}

// Only interactions that have the length <= max_length_in_obj will contribute to the obj function
void fit_hi_c_data_model (const ValueMatrix& hi_c_mat, int left, int right, int shortest_path_mode, int max_length_in_obj, int shortest_path_window_size, HiCDataModel& hi_c_data_model) {
	//cout << hi_c_mat.size() << "\t" << left << "\t" << right << endl;
	int length = right - left + 1;
	ValueList norm_val(length);
	for (int i = 0; i < length; i++)
		norm_val[i] = ((hi_c_mat[left + i][left + i] > 0)? hi_c_mat[left + i][left + i]:ZERO_ENTRY);

	ValueMatrix distance_mat(length);
	for (int i = 0; i < length; i++) {
		distance_mat[i].resize(length);
		for (int j = 0; j < length; j++) {
			int pos_i = left + i, pos_j = left + j;
			distance_mat[i][j] = pow(((hi_c_mat[pos_i][pos_j] > 0)? hi_c_mat[pos_i][pos_j]:ZERO_ENTRY)/sqrt(norm_val[i]*norm_val[j]),-1);
			
		}
	}

	ValueMatrix shortest_dist_mat;
	if (shortest_path_mode == FLOYD_SHORTEST_PATH)
		Floyd_APSP(distance_mat, shortest_dist_mat);		
	else if (shortest_path_mode == APPROX_SHORTEST_PATH)
		limited_shortest_path(distance_mat, max_length_in_obj, shortest_path_window_size, shortest_dist_mat); 
	else
		cout << "Can not find the mode of finding the shortest path " << shortest_path_mode << endl;
	//for (int i = 0; i < distance_mat.size(); i++)
	//	cout << i << "\t" << shortest_dist_mat[i].size() << endl;
	//cout << "Complete finding all the shortest paths" << endl;

	CPXENVptr env = NULL;
	CPXLPptr lp = NULL;

	int alpha_size = length;
	int beta_size = 1;
	int e_size = ((length > max_length_in_obj)? ((length - max_length_in_obj)*max_length_in_obj + max_length_in_obj*(max_length_in_obj - 1)/2) : (length*(length - 1)/2));	// slack
	int r_size = length;															// resistance

	int var_num = alpha_size + beta_size + e_size + r_size; 
	int status = 0;
	env = CPXopenCPLEX (&status);
	/* Turn on output to the screen */
	status = CPXsetintparam (env, CPXPARAM_ScreenOutput, CPX_OFF);
	/* Turn on data checking */
	status = CPXsetintparam (env, CPXPARAM_Read_DataCheck, CPX_ON);
	//status = CPXsetintparam (env, CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
	/* Create the problem */
	lp = CPXcreateprob (env, &status, "linear_model");
	/* Populate the problem */
	/* Problem is minimization */
	status = CPXchgobjsen (env, lp, CPX_MIN);
	// For obj and bounds
	double* obj = (double*) malloc(var_num*sizeof(double));
	double* lb = (double*) malloc(var_num*sizeof(double));
	double* ub = (double*) malloc(var_num*sizeof(double));
	char** colname = (char**) malloc(var_num*sizeof(char*));
	
	int  current_id = 0;
	// alpha
	for (int i = 0; i < alpha_size; i++) {
		lb[current_id] = -CPX_INFBOUND;
		ub[current_id] = CPX_INFBOUND;
		//lb[current_id] = log(hi_c_mat[left + i][left + i]);
		//ub[current_id] = log(hi_c_mat[left + i][left + i] + 1e-10);

		colname[current_id] = (char*) malloc(50*sizeof(char));
		sprintf(colname[current_id], "a_%i\0", i);
		obj[current_id] = 0;
		current_id++;
	}
	// beta
	lb[current_id] = -2;
	ub[current_id] = -1;
	//lb[current_id] = -1.00001;
	//ub[current_id] = -1;
	
	colname[current_id] = (char*) malloc(50*sizeof(char));
	sprintf(colname[current_id], "beta\0");
	obj[current_id] = 0;
	current_id++;
	// slack
	int tmp = current_id;
	int cell_num_tmp = 0;
	for (int i = 0; i < length; i++) {
		int max_i = ((i + max_length_in_obj < length)? (i + max_length_in_obj):(length-1));
		for (int j = i + 1; j <= max_i; j++) {
			lb[current_id] = 0;
			ub[current_id] = CPX_INFBOUND;
			colname[current_id] = (char*) malloc(50*sizeof(char));
			sprintf(colname[current_id], "e_%i,%i\0", i, j);
			obj[current_id] = 1;
			current_id++;
			cell_num_tmp += (4 + j - i - 1);
		}
	}
	//cout << "#Var " << current_id - tmp << "\t" << e_size << endl;
	// resistance
	for (int i = 0; i < length; i++) {
		lb[current_id] = 0;
		ub[current_id] = CPX_INFBOUND;
		colname[current_id] = (char*) malloc(50*sizeof(char));
		sprintf(colname[current_id], "r_%i\0", i);
		obj[current_id] = 0;
		current_id++;
	}
	/* Obj and the bound */
	status = CPXnewcols (env, lp, var_num, obj, lb, ub, NULL, colname);
	if (status) {
		cout << "LINH: Fail to set up (bound and obj) the LP problem!" << endl;
		cout << "(left,right) = " << left << "\t" << right << endl;
	}
	//else 
	//	cout << "LINH: Add all variables successfully!" << endl;
	/* Add constraints */
	int row_num = 2*e_size;
	int cell_num = 2*cell_num_tmp;
	int* rmatbeg = (int*) malloc(row_num*sizeof(int));
	int* rmatind = (int*) malloc(cell_num*sizeof(int));
	double* rmatval = (double*) malloc(cell_num*sizeof(double));
	double* rhs = (double*) malloc(row_num*sizeof(double));
	char* sense = (char*) malloc(row_num*sizeof(char));
	char** rowname = (char**) malloc(row_num*sizeof(char*));
	
	int current_row_id = 0;
	int current_cell_id = 0;
	//cout << "Total: " << length << endl;
	for (int i = 0; i < length; i++) {
		int i_max = ((i + max_length_in_obj < length)? (i + max_length_in_obj):(length-1));
		//if (i > n - 10)
		//	cout << i << endl;
		for (int j = i + 1; j <= i_max; j++) {
			double logH = ((hi_c_mat[left + i][left + j] > 0)? log(hi_c_mat[left + i][left + j]) : log(ZERO_ENTRY));
			// (a_i + a_j)/2 + beta*log(d_ij) - sum(r_k) - e_ij <= log(H_ij)
			rmatbeg[current_row_id] = current_cell_id;
			sense[current_row_id] = 'L';
			rhs[current_row_id] = logH;
			rowname[current_row_id] = (char*) malloc(20*sizeof(char));
			sprintf(rowname[current_row_id], "e_%i,%i\0", i, j);
			rmatind[current_cell_id] = i;
			rmatval[current_cell_id] = 0.5;
			rmatind[current_cell_id + 1] = j;
			rmatval[current_cell_id + 1] = 0.5;
			rmatind[current_cell_id + 2] = alpha_size;
			rmatval[current_cell_id + 2] = log(shortest_dist_mat[j][j - i - 1]);
			rmatind[current_cell_id + 3] = alpha_size + 1 + round(current_row_id/2);
			rmatval[current_cell_id + 3] = -1;
			for (int k = i + 1; k < j; k++) {
				rmatind[current_cell_id + 3 + k - i] = alpha_size + 1 + e_size + k;
				rmatval[current_cell_id + 3 + k - i] = -1; 
			}
			current_row_id++;
			current_cell_id += 3 + (j-i);
			// -(a_i + a_j)/2 - beta*log(d_ij) + sum(r_k) - e_ij <= -log(H_ij)
			rmatbeg[current_row_id] = current_cell_id;
			sense[current_row_id] = 'L';
			rhs[current_row_id] = -logH;
			rowname[current_row_id] = (char*) malloc(20*sizeof(char));
			sprintf(rowname[current_row_id], "e_%i,%i\0", i, j);
			rmatind[current_cell_id] = i;
			rmatval[current_cell_id] = -0.5;
			rmatind[current_cell_id + 1] = j;
			rmatval[current_cell_id + 1] = -0.5;
			rmatind[current_cell_id + 2] = alpha_size;
			rmatval[current_cell_id + 2] = -log(shortest_dist_mat[j][j - i - 1]);
			rmatind[current_cell_id + 3] = alpha_size + 1 + round(current_row_id/2);
			rmatval[current_cell_id + 3] = -1;
			for (int k = i + 1; k < j; k++) {
				rmatind[current_cell_id + 3 + k - i] = alpha_size + 1 + e_size + k;
				rmatval[current_cell_id + 3 + k - i] = 1; 
			}
			current_row_id++;
			current_cell_id += 3 + (j-i);
		}
	}
	//cout << "Estimate: " << row_num << "\t" << cell_num << endl;
	//cout << "Iterate: " << current_row_id << "\t" << current_cell_id << endl;
	
	//status = CPXaddrows (env, lp, 0, row_num, cell_num, rhs, sense, rmatbeg, rmatind, rmatval, NULL, rowname);
	status = CPXaddrows (env, lp, 0, current_row_id, current_cell_id, rhs, sense, rmatbeg, rmatind, rmatval, NULL, rowname);
	if (status) {
		cout << "LINH: Fail to add constraints to LP problem!" << endl;
		/*for (int i = 0; i < var_num; i++)
			cout << i << "\t" << lb[i] << "\t" << ub[i] << "\t" << obj[i] << "\t" << colname[i] << endl;
		cout << "=============================" << endl;
		for (int i = 0; i < current_row_id; i++)
			cout << i << "\t" << "rmatbeg: " << rmatbeg[i] << "\tsense: " << sense[i] << "\trhs: " << rhs[i] << "\t" << rowname[i] << endl;
		for (int i = 0; i < current_cell_id; i++)
			cout << i << "\t" << "rmatind: " << rmatind[i] << "\t rmatval: " << rmatval[i] << endl;
		*/
	}
	else {
		//cout << "LINH: Add all constraints successfully!" << endl;
	}	
	/* Solve the problem */
	status = CPXlpopt (env, lp); 
	/* Get the actual size of the problem */
	int cur_numrows = CPXgetnumrows (env, lp);
	int cur_numcols = CPXgetnumcols (env, lp);
	double* x = (double*) malloc (cur_numcols*sizeof(double));
	double* slack = (double*) malloc (cur_numrows*sizeof(double));
	double* dj = (double*) malloc (cur_numcols*sizeof(double));
	double* pi = (double*) malloc (cur_numrows*sizeof(double));
	int solstat;
	double objval;
	/* Get the solution */
	//status = CPXwriteprob (env, lp, "linear_model.lp", NULL);
	status = CPXsolution (env, lp, &solstat, &objval, x, pi, slack, dj);
	if (status)
		cout << "LINH: Fail to solve the LP problem!" << endl;
	else {
		//cout << "LINH: Final obj value = " << objval << endl;
		//for (int i = 0; i < var_num; i++)
		//	cout << colname[i] << "\t" << x[i] << endl;
	}
	
	// Export the prediction
	hi_c_data_model.beta = x[alpha_size];
	resize_and_fill(hi_c_data_model.alpha, length, 0);
	resize_and_fill(hi_c_data_model.resistance, length, 0);
	for (int i = 0; i < length; i++) {
		hi_c_data_model.alpha[i] = x[i];
		hi_c_data_model.resistance[i] = x[alpha_size + 1 + e_size + i];
	}

	/* Free the problem */
	CPXfreeprob(env, &lp);
	/* Free the environment */
	CPXcloseCPLEX (&env);
	/* Free all memory */
	free(obj);
	free(ub);
	free(lb);
	for (int i = 0 ; i < var_num; i++)
		free(colname[i]);
	free(colname);
	free(x);
	free(slack);
	free(dj);
	free(pi);
	free(rmatbeg);
	free(rmatind);
	free(rmatval);
	free(rhs);
	free(sense);
	for (int i = 0; i < current_row_id; i++)
		free(rowname[i]);
	free(rowname);
}

void Dijkstra (const ValueMatrix& hi_c_mat, int left, int right, int src, ValueList& shortest_path_val_list) {
	int length = right - left + 1;
	ValueList norm_val(length);
	for (int i = 0; i < length; i++)
		norm_val[i] = ((hi_c_mat[left + i][left + i] > 0)? hi_c_mat[left + i][left + i]:ZERO_ENTRY);

	ValueMatrix graph_val(length);
	for (int i = 0; i < length; i++) {
		graph_val[i].resize(length);
		for (int j = 0; j < length; j++) {
			int pos_i = left + i, pos_j = left + j;
			graph_val[i][j] = pow(((hi_c_mat[pos_i][pos_j] > 0)? hi_c_mat[pos_i][pos_j]:ZERO_ENTRY)/sqrt(norm_val[i]*norm_val[j]),-1);
		}
	}
	src = src - left;
	IdList Q, is_in_Q, prev;
	ValueList dist;
	for (int i = 0; i < graph_val.size(); i++) {
		dist.push_back(1e10);
		Q.push_back(i);
		is_in_Q.push_back(1);
		prev.push_back(-1);	
	}
	dist[src] = 0;
	while (!Q.empty()) {
		IdList::iterator min_it = Q.begin();
		for (IdList::iterator it = Q.begin(); it != Q.end(); it++) {
			if (dist[(*it)] < dist[(*min_it)])
				min_it = it;
			//cout << (*it) << "\t" << dist[(*it)] << endl;
		}
		int u = (*min_it);
		Q.erase(min_it);
		is_in_Q[u] = 0;
		//cout << u << endl;
		for (int v = 0; v < graph_val.size(); v++) 
			if (is_in_Q[v] > 0) {
				double alt = dist[u] + graph_val[u][v];
				if (alt < dist[v]) {
					dist[v] = alt;
					prev[v] = u;
				}
		}
	}
	shortest_path_val_list = dist;
	//for (int i = 0; i < shortest_path_val_list.size(); i++)
	//	if (shortest_path_val_list[i] < 1e-10)
	//		cout << "\t Dij " << left << "\t" << right << "\t" << src << "\t" << i << "\t" << shortest_path_val_list[i] << endl;
	//for (int i = 0; i < graph_val.size(); i++) {
	//	for (int j = 0; j < graph_val.size(); j++)
	//		cout << graph_val[i][j] << "\t";
	//		cout << endl;
	//}
	//for (int i = 0; i < shortest_path_val_list.size(); i++)
	//	cout << shortest_path_val_list[i] << "\t";
	//cout << endl;
}

ValueList predict_3D (const ValueMatrix& hi_c_mat, int left_1, int right_1, int left_2, int right_2, int shortest_path_mode, int min_length, int window_size, bool print_out) {
	if (left_1 >= right_1 || left_2 >= right_2)
		cout << "ERROR: Can not predict because the left reigon or the right region is so small!" << endl;

	HiCDataModel hi_c_data_model_1, hi_c_data_model_2;
	fit_hi_c_data_model(hi_c_mat, left_1, right_1, shortest_path_mode, min_length, window_size, hi_c_data_model_1);
	fit_hi_c_data_model(hi_c_mat, left_2, right_2, shortest_path_mode, min_length, window_size, hi_c_data_model_2);
	
	ValueList shortest_path_list_1, shortest_path_list_2;
	Dijkstra(hi_c_mat, left_1, right_1, right_1, shortest_path_list_1);
	Dijkstra(hi_c_mat, left_2, right_2, left_2, shortest_path_list_2);

	// print the boundary score
	//cout << "Boundary" << endl;
	//for (int i = 0; i < left_length; i++)
	//	cout << (left_1 + i)*5000 << "\t" << boundary_score_left[i] << endl;
	//for (int i = 0; i < right_length; i++)
	//	cout << (left_2 + i)*5000 << "\t" << boundary_score_right[i] << endl;

	// Export the prediction
	int length_1 = right_1 - left_1 + 1, length_2 = right_2 - left_2 + 1;
	//pred_mat.resize(right_1 - left_1 + 1);
	ValueList obs, pred;
	for (int i = 0; i < length_1; i++) {
		//pred_mat[i].clear();
		for (int j = 0; j < length_2; j++) {
			int pos_i = left_1 + i, pos_j = left_2 + j;
			double obs_val = hi_c_mat[pos_i][pos_j];
			double unit_distance = (shortest_path_list_1[length_1 - 2] + shortest_path_list_2[1])/2;	// One unit here
			double new_distance = shortest_path_list_1[i] + shortest_path_list_2[j] + unit_distance;
			double log_pred = (hi_c_data_model_1.alpha[i] + hi_c_data_model_2.alpha[j])/2 + ((hi_c_data_model_1.beta + hi_c_data_model_2.beta)/2)*log(new_distance);
			for (int k = i + 1; k < right_1 - left_1; k++)
				log_pred -= hi_c_data_model_1.resistance[k];
			for (int k = 1; k < j; k++)
				log_pred -= hi_c_data_model_2.resistance[k];
			double pred_val = exp(log_pred);
			//cout << left_shortest_path_list[i] << "\t" << right_shortest_path_list[j] << "\t" << pred_val << endl;
			//pred_mat[i].push_back(pred_val);
			obs.push_back(obs_val);
			pred.push_back(pred_val);			
		}
		//cout << endl;
	}
	if (print_out) {
		string file_name = "Debug/New_" + num_to_string(left_1) + "_" + num_to_string(right_1) + ".dat";
		ofstream out_file(file_name.c_str());
		for (int i = 0; i < pred.size(); i++)
			out_file << pred[i] << "\t" << obs[i] << endl;
		out_file.close();
	}

	ValueList tmp;
	tmp.push_back(RSQ(pred, obs));
	tmp.push_back(average_log_ratio(pred, obs));
	tmp.push_back(L1(pred, obs));
	tmp.push_back(L2(pred, obs));
	return tmp;
}

bool predict_3D (const ValueMatrix& hi_c_mat, const HiCDataModelList& hi_c_data_model_list, int left_1, int right_1, int left_2, int right_2, ValueMatrix& pred_mat) {	
	if (left_1 >= right_1 || left_2 >= right_2)
		cout << "ERROR: Can not predict because the left reigon or the right region is so small!" << endl;

	for (int k = 0; k < hi_c_data_model_list.size(); k++) {
		if (left_1 >= hi_c_data_model_list[k].left && left_2 >= right_1 && right_2 <= hi_c_data_model_list[k].get_right()) {
			ValueList shortest_path_list_1, shortest_path_list_2;
			Dijkstra(hi_c_mat, left_1, right_1, right_1, shortest_path_list_1);
			Dijkstra(hi_c_mat, left_2, right_2, left_2, shortest_path_list_2);

			int domain_left = hi_c_data_model_list[k].left;
			int length_1 = right_1 - left_1 + 1, length_2 = right_2 - left_2 + 1;
			pred_mat.resize(right_1 - left_1 + 1);
			for (int i = 0; i < length_1; i++) {
				pred_mat[i].clear();
				for (int j = 0; j < length_2; j++) {
					int pos_i = left_1 + i, pos_j = left_2 + j;
					double unit_distance = (shortest_path_list_1[length_1 - 2] + shortest_path_list_2[1])/2;	// One unit here
					double new_distance = shortest_path_list_1[i] + shortest_path_list_2[j] + unit_distance;
					double log_pred = (hi_c_data_model_list[k].alpha[i + left_1 - domain_left] + hi_c_data_model_list[k].alpha[j + left_2 - domain_left])/2 + hi_c_data_model_list[k].beta*log(new_distance);
					for (int r = i + 1; r < right_1 - left_1; r++)
						log_pred -= hi_c_data_model_list[k].resistance[r + left_1 - domain_left];
					for (int r = 1; r < j; r++)
						log_pred -= hi_c_data_model_list[k].resistance[r + left_2 - domain_left];
					double pred_val = exp(log_pred);
					//cout << left_shortest_path_list[i] << "\t" << right_shortest_path_list[j] << "\t" << pred_val << endl;
					pred_mat[i].push_back(pred_val);
				}
			}
			return true;
		}
	}	
	cout << "Error: Can not find the Hi-C data model! " << left_1 << "\t" << right_1 << "\t" << left_2 << "\t" << right_2 << endl;
	int length_1 = right_1 - left_1 + 1, length_2 = right_2 - left_2 + 1;
	pred_mat.resize(length_1);
	for (int i = 0; i < length_1; i++) {
		pred_mat[i].clear();
		for (int j = 0; j < length_2; j++)
			pred_mat[i].push_back(ZERO_ENTRY);
	}
	return false;
}

void export_hi_c_model_list (const ValueMatrix& hi_c_mat, int window_size, int shift_length, int min_length, int shortest_path_window, string filename) {
	HiCDataModelList hi_c_data_model_list;
	int current_left = 0;
	while (current_left < hi_c_mat.size()) {
		HiCDataModel hi_c_data_model;
		int right_pos = (current_left + window_size < hi_c_mat.size())? (current_left + window_size): (hi_c_mat.size() - 1);
		cout << hi_c_mat.size() << "\t" << current_left << "\t" << right_pos << endl;
		fit_hi_c_data_model(hi_c_mat, current_left, right_pos, APPROX_SHORTEST_PATH, min_length, shortest_path_window, hi_c_data_model);
		hi_c_data_model.left = current_left;
		hi_c_data_model_list.push_back(hi_c_data_model);
		current_left += shift_length;	// Correct here for the final version
	}
	write_hi_c_data_model_list(filename, hi_c_data_model_list);
}

ValueList estimate_fusion_level (const ValueMatrix& hi_c_mat, const HiCDataModelList& hi_c_data_model_list, const DeletionList& del_list, int resolution, int window_size, const ValueList& significance_threshold_list, ValueMatrix& fusion_level_matrix) {
	fusion_level_matrix.resize(significance_threshold_list.size());
	ValueList total_fusion_list (significance_threshold_list.size(), 0);
	for (int i = 0; i < significance_threshold_list.size(); i++)
		resize_and_fill(fusion_level_matrix[i], del_list.size(), 0);
	for (int i = 0; i < del_list.size(); i++) {
		//if (i % 10 == 0)
		//	cout << i << endl;
		int del_left = floor(del_list[i].left/resolution), del_right = ceil(del_list[i].right/resolution) - 1;
		if (del_right >= hi_c_mat.size())
			cout << "Error: Deletion right = " << del_right << ", hi-c matrix size = " << hi_c_mat.size() << endl;
		int left_window_size = (del_left >= window_size)? window_size:del_left;
		int right_window_size = (del_right + window_size < hi_c_mat.size())? window_size: hi_c_mat.size() - 1 - del_right;
		ValueMatrix pred;
		bool is_complete = predict_3D(hi_c_mat, hi_c_data_model_list, del_left - left_window_size, del_left, del_right, del_right + right_window_size, pred);
		for (int l = 0; l <= left_window_size; l++)
			for (int r = 0; r <= right_window_size; r++) {
				double old_val = hi_c_mat[del_left - left_window_size + l][del_right + r], new_val = pred[l][r];
				for (int t = 0; t < significance_threshold_list.size(); t++)
					if (old_val < significance_threshold_list[t] && new_val > significance_threshold_list[t]) 
						fusion_level_matrix[t][i]++;
			}
		for (int t = 0; t < significance_threshold_list.size(); t++)
			total_fusion_list[t] += fusion_level_matrix[t][i];
	}
	return total_fusion_list;
}

double conventional_count_fusion (const DeletionList& del_list, const TADList& tad_list, ValueList& fusion_count_list, bool is_printed) {
	resize_and_fill(fusion_count_list, del_list.size(), 0);
	double count = 0;
	for (int k = 0; k < del_list.size(); k++) {
		bool is_covered_left = false, is_covered_right = false;
		for (int i = 0; i < tad_list.size(); i++)
			if (del_list[k].left >= tad_list[i].left && del_list[k].left <= tad_list[i].right && del_list[k].right >= tad_list[i].right) {
					is_covered_left = true;
					break;
			}
		for (int i = 0; i < tad_list.size(); i++)
			if (del_list[k].left <= tad_list[i].left && del_list[k].right >= tad_list[i].left && del_list[k].right <= tad_list[i].right) {
					is_covered_right = true;
					break;
			}
		if (is_covered_left && is_covered_right) {			
			count++;
			fusion_count_list[k] = 1;
			if (is_printed)
				cout << del_list[k].left << "\t" << del_list[k].right << endl;
		}
	}
	return count;
}

double conventional_count_fusion (const DeletionList& del_list, const TADBoundaryList& tad_boundary_list, ValueList& fusion_count_list, bool is_printed) {
	resize_and_fill(fusion_count_list, del_list.size(), 0);
	double count = 0;
	for (int k = 0; k < del_list.size(); k++) {
		for (int i = 0; i < tad_boundary_list.size(); i++)
			if ((del_list[k].left >= tad_boundary_list[i].left && del_list[k].left <= tad_boundary_list[i].right && del_list[k].right >= tad_boundary_list[i].right) 
				|| (del_list[k].right >= tad_boundary_list[i].left && del_list[k].right <= tad_boundary_list[i].right && del_list[k].left <= tad_boundary_list[i].left) 
				|| (del_list[k].left <= tad_boundary_list[i].left && del_list[k].right >= tad_boundary_list[i].right)) { 
					if (fusion_count_list[k] < tad_boundary_list[i].score)
						fusion_count_list[k] = tad_boundary_list[i].score;
			}
		count += fusion_count_list[k];
		if (is_printed)
			cout << del_list[k].left << "\t" << del_list[k].right << endl;
	}
	return count;
}

void sampling_deletion (const DeletionList& del_list, double chr_size, int window_size, DeletionList& sampling_del_list) {
	sampling_del_list.resize(del_list.size());
	for (int del = 0; del < del_list.size(); del++) {
		int l = del_list[del].right - del_list[del].left;
		double r = (rand() % 10001)/10000.0;
		sampling_del_list[del].left = round(r*(chr_size - l - 1 - 2*window_size)) + window_size;
		sampling_del_list[del].right = sampling_del_list[del].left + l;		
	}
}

void prediction_experiment(int chr, int bin_min, int bin_max, const IdList& window_size_list, bool print_out) {
	string raw_data_file_folder = "/share/hormozdiarilab/Data/HiC/GM12878_combined/5kb_resolution_intrachromosomal/chr";
	//string raw_data_file_folder = "/share/hormozdiarilab/Data/HiC/IMR90/5kb_resolution_intrachromosomal/chr";
	string chr_name = ((chr < 23)? num_to_string(chr):"X");
	string chr_raw_data_filename = raw_data_file_folder + chr_name +  "/MAPQGE30/chr" + chr_name + "_5kb.RAWobserved";
	ValueMatrix hi_c_matrix;
	read_hi_c_data(chr_raw_data_filename, RESOLUTION, hi_c_matrix);

		
	int max_window_size = window_size_list[window_size_list.size() - 1];
	if (bin_min < max_window_size)
		bin_min = max_window_size;
	if (bin_max > hi_c_matrix.size() - max_window_size - 1)
		bin_max = hi_c_matrix.size() - max_window_size - 1;	
	IdList bin_list;
	for (int i = bin_min; i <= bin_max; i = i + 5) 
		bin_list.push_back(i);

	DistributionSummaryMatrix naive_error_dist_mat(window_size_list.size()), new_error_dist_mat(window_size_list.size());
	for (int w = 0; w < window_size_list.size(); w++) {
		naive_error_dist_mat[w].resize(4);
		new_error_dist_mat[w].resize(4);
	}
	
	string out_filename = "Debug/error_comparison_" + num_to_string(chr) + "_" + num_to_string(bin_list[0]) + "_" + num_to_string(bin_list[bin_list.size() - 1]) + ".dat";
	ofstream error_out_file (out_filename.c_str());
	error_out_file << "Distance\tNaive\tNew" << endl;
	
	for (int i = 0; i < bin_list.size(); i++) {
		if (bin_list[i] % 100 == 0)
			cout << bin_list[i] << endl;		
		int bin = bin_list[i];
		error_out_file << bin;
		ValueMatrix naive_result(window_size_list.size()), new_result(window_size_list.size());
		for (int w = 0; w < window_size_list.size(); w++) { 
			ValueList naive_result = naive_prediction(hi_c_matrix, bin - window_size_list[w], bin, bin + 1, bin + window_size_list[w], print_out);
			ValueList new_result = predict_3D(hi_c_matrix, bin - window_size_list[w], bin, bin + 1, bin + window_size_list[w], 
				FLOYD_SHORTEST_PATH, MAX_LENGTH_OBJ, SHORTEST_PATH_WINDOW_SIZE, print_out);
			for (int k = 0; k < 4; k++) {
				naive_error_dist_mat[w][k].add_value(naive_result[k]);
				new_error_dist_mat[w][k].add_value(new_result[k]);
				error_out_file << "\t" << naive_result[k] << "\t" << new_result[k];
			}
		}
		error_out_file << endl;		
	}
	error_out_file.close();
	out_filename = "Debug/pred_comparison_summary_" + num_to_string(chr) + "_" + num_to_string(bin_list[0]) + "_" + num_to_string(bin_list[bin_list.size() - 1]) + ".dat";
	ofstream summary_out (out_filename.c_str());
	summary_out << "Naive(R)\tNew(R)\tNaive(log)\tNew(log)\tNaive(L1)\tNew(L1)\tNaive(L2)\tNew(L2)" << endl;
	for (int w = 0; w < window_size_list.size(); w++) {
		for (int k = 0; k < 4; k++)
			summary_out << ((k > 0)? "\t":"") 
				<< naive_error_dist_mat[w][k].mean << " +/- " << naive_error_dist_mat[w][k].std << "\t" 
				<< new_error_dist_mat[w][k].mean << "+/-" << new_error_dist_mat[w][k].std;
		summary_out << endl;
	}
	summary_out.close();
}

void export_all_hi_c_model () {
	int chr_min = 1;
	int chr_max = 23;
	int hi_c_model_window_size = round(6000000/RESOLUTION);
	// For generating Hi-C data model list
	string data_file_folder = "/share/hormozdiarilab/Data/HiC/GM12878_combined/5kb_resolution_intrachromosomal/chr" ;
	string model_file_folder = "/share/hormozdiarilab/Data/HiC/Model/GM12878/";
	for (int chr = chr_min; chr <= chr_max; chr++) { 
		string chr_name = ((chr < 23)? num_to_string(chr):"X");
		string chr_raw_data_filename = data_file_folder + chr_name +  "/MAPQGE30/chr" + chr_name + "_5kb.RAWobserved"; 
		ValueMatrix hi_c_mat;
		read_hi_c_data(chr_raw_data_filename, RESOLUTION, hi_c_mat);
		export_hi_c_model_list (hi_c_mat, hi_c_model_window_size, round(hi_c_model_window_size/2), MAX_LENGTH_OBJ, SHORTEST_PATH_WINDOW_SIZE, model_file_folder + "chr" + chr_name + ".model");
		cout << "Complete generating Hi-C model of chromosome " << chr_name << endl;
	}
}

void deletion_filter (DeletionMatrix& original_del_mat, const DeletionMatrix& filter_del_mat) {
	for (int chr = 0; chr < original_del_mat.size(); chr++) {
		//for (int l = 0; l < filter_del_mat[chr].size(); l++)
		//	cout << chr + 1 << "\t" << filter_del_mat[chr][l].left << "\t" << filter_del_mat[chr][l].right << endl; 
		DeletionList tmp;
		for (int del = 0; del < original_del_mat[chr].size(); del++) {
			bool is_repeated = false;
			for (int k = 0; k < filter_del_mat[chr].size(); k++)
				if (original_del_mat[chr][del].left == filter_del_mat[chr][k].left && original_del_mat[chr][del].right == filter_del_mat[chr][k].right) {
					is_repeated = true;
					break;
				}
			if (!is_repeated)
				tmp.push_back(original_del_mat[chr][del]);
		}
		//cout << chr << "\t" << original_del_mat[chr].size() << "\t" << tmp.size() << endl;
		original_del_mat[chr] = tmp;
	}
}

void deletion_filter (DeletionMatrix& original_del_mat, int max_del_length) {
	for (int chr = 0; chr < original_del_mat.size(); chr++) {
		DeletionList tmp;
		for (int del = 0; del < original_del_mat[chr].size(); del++) {
			if (original_del_mat[chr][del].right - original_del_mat[chr][del].left <= max_del_length)
				tmp.push_back(original_del_mat[chr][del]);
		}
		//cout << chr << "\t" << original_del_mat[chr].size() << "\t" << tmp.size() << endl;
		original_del_mat[chr] = tmp;
	}
}

void fusion_experiment (string exp_name, string del_filename, bool is_filtered, int max_del_length, const ValueList& contact_threshold_list, int final_contact_threshold_index, const ValueList& fusion_threshold_list, int sample_num) {
	int chr_min = 1;	// 1;
	int chr_max = 23; 	// 23;
	int min_TAD_size = 10*RESOLUTION;
	int min_del_size = 2*RESOLUTION;
	int window_size = round(250000/RESOLUTION);
	string model_file_folder = "/share/hormozdiarilab/Data/HiC/Model/GM12878/";
	string raw_data_file_folder = "/share/hormozdiarilab/Data/HiC/GM12878_combined/5kb_resolution_intrachromosomal/chr" ;

	if (final_contact_threshold_index < 0 || final_contact_threshold_index >= contact_threshold_list.size())
		cout << "ERROR: Not valid index of contact threshold " << final_contact_threshold_index << " vs size = " << contact_threshold_list.size() << endl;

	DeletionMatrix del_mat;
	read_deletion_mutation_file(del_filename, min_del_size, del_mat);

	if (is_filtered) {
		string GM12878_variant_filename = "/share/hormozdiarilab/Data/ValidatedVariants/1000G_SV_Phase3/NA12878_Del.BED";
		DeletionMatrix GM12878_del_mat;
		read_deletion_mutation_file(GM12878_variant_filename, min_del_size, GM12878_del_mat);
		deletion_filter(del_mat, GM12878_del_mat);
	}
	deletion_filter(del_mat, max_del_length);
	//for (int chr = 0; chr < 23; chr++)
	//	for (int del = 0; del < del_mat[chr].size(); del++)
	//		cout << chr + 1 << "\t" << del_mat[chr][del].left << "\t" << del_mat[chr][del].right << endl;

	string CTCF_filename = "/share/hormozdiarilab/Data/CTCF_Binding/UW_hg19_GM12878_CTCFBSDB.bed";
	DeletionMatrix CTCF_mat;
	read_CTCF(CTCF_filename, CTCF_mat);
	
	StringList tool_name_list = {"Arrowhead", "CaTCH", "Amaratus", "InsulationScore"};
	ValueList tool_type_list = {0, 0, 0, 1};	// 0: discrete, 1: continuous
	StringList TAD_info_filename_list = {"GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt", "CaTCH_TAD.dat", "Amaratus_TAD.dat", "insulation_score_TAD.dat"};
	vector<TADMatrix> TAD_mat_list(tool_name_list.size());
	vector<TADBoundaryMatrix> TAD_boundary_mat_list(tool_name_list.size());

	ValueMatrix other_observation_fusion_count_matrix(tool_name_list.size());			// tool_num x CHR_NUM
	ValueMatrix our_observation_fusion_level_matrix(CHR_NUM);					// CHR_NUM x contact_threshold
	ValueMatrix our_observation_fusion_count_matrix(CHR_NUM);					// CHR_NUM x fusion_threshold

	ValueList other_total_observation_fusion_level(tool_name_list.size(), 0);			// tool_num
	ValueList our_total_observation_fusion_level(contact_threshold_list.size(), 0);			// contact_threshold
	ValueList our_total_observation_fusion_count(fusion_threshold_list.size(), 0);			// fusion_threshold

	DistributionMatrix other_permutation_fusion_dist_matrix(tool_name_list.size());			// tool_num x CHR_NUM
	DistributionMatrix our_permutation_fusion_level_dist_matrix(CHR_NUM);				// CHR_NUM x contact_threshold
	DistributionMatrix our_permutation_fusion_count_dist_matrix(CHR_NUM);				// CHR_NUM x fusion_threshold

	ValueMatrix other_total_permutation_fusion_level(tool_name_list.size());			// tool_num x sample_num
	ValueMatrix our_total_permutation_fusion_level(sample_num);					// sample_num x contact_threshold
	ValueMatrix our_total_permutation_fusion_count(sample_num);					// sample_num x fusion_threshold

	srand (time(NULL));
	// Read TAD info of different tools
	string TAD_info_folder = "/share/hormozdiarilab/Data/HiC/TAD_Annotation/GM12878/";
	for (int tool = 0; tool < tool_type_list.size(); tool++) {			
		if (tool_type_list[tool] < 1) 
			read_TAD_data_with_Juice_format (TAD_info_folder + TAD_info_filename_list[tool], TAD_mat_list[tool], min_TAD_size);
		else 
			read_TAD_boundary_with_4_col_format (TAD_info_folder + TAD_info_filename_list[tool], TAD_boundary_mat_list[tool]);		
		other_observation_fusion_count_matrix[tool].resize(CHR_NUM, 0);
		other_permutation_fusion_dist_matrix[tool].resize(CHR_NUM);
		resize_and_fill(other_total_permutation_fusion_level[tool], sample_num, 0);		

	}
	for (int s = 0; s < sample_num; s++) {
		resize_and_fill(our_total_permutation_fusion_level[s], contact_threshold_list.size(), 0);
		resize_and_fill(our_total_permutation_fusion_count[s], fusion_threshold_list.size(), 0);
	}
	for (int chr = 0; chr < CHR_NUM; chr++) {
		our_observation_fusion_level_matrix[chr].resize(contact_threshold_list.size());
		our_observation_fusion_count_matrix[chr].resize(fusion_threshold_list.size());
		our_permutation_fusion_level_dist_matrix[chr].resize(contact_threshold_list.size());
		our_permutation_fusion_count_dist_matrix[chr].resize(fusion_threshold_list.size());
	}

	vector<ValueMatrix> other_observation_fusion_level_matrix_list(CHR_NUM);			// CHR_NUM x tool_num x del_num
	vector<ValueMatrix> our_observation_fusion_level_matrix_list(CHR_NUM);				// CHR_NUM x contact_threshold x del_num
	for (int chr = chr_min; chr <= chr_max; chr++) { 
		string chr_name = ((chr < 23)? num_to_string(chr):"X");
		string chr_raw_data_filename = raw_data_file_folder + chr_name +  "/MAPQGE30/chr" + chr_name + "_5kb.RAWobserved"; 
		ValueMatrix hi_c_mat;
		if (!del_mat[chr - 1].empty())
			read_hi_c_data(chr_raw_data_filename, RESOLUTION, hi_c_mat);		
		// Observation 
		HiCDataModelList hi_c_data_model_list;
		read_hi_c_data_model_list(model_file_folder + "chr" + chr_name + ".model", hi_c_data_model_list);
		cout << "Complete loading Hi-C model of chromosome " << chr_name << endl;

		other_observation_fusion_level_matrix_list[chr - 1].resize(tool_name_list.size());
		our_observation_fusion_level_matrix_list[chr - 1].resize(contact_threshold_list.size());	
		for (int tool = 0; tool < tool_name_list.size(); tool++) {
			if (tool_type_list[tool] < 1)
				other_observation_fusion_count_matrix[tool][chr - 1] = conventional_count_fusion (del_mat[chr - 1], TAD_mat_list[tool][chr - 1], other_observation_fusion_level_matrix_list[chr - 1][tool], false);
			else 
				other_observation_fusion_count_matrix[tool][chr - 1] = conventional_count_fusion (del_mat[chr - 1], TAD_boundary_mat_list[tool][chr - 1], other_observation_fusion_level_matrix_list[chr - 1][tool], false);
			// total
			other_total_observation_fusion_level[tool] += other_observation_fusion_count_matrix[tool][chr - 1];
		}
		our_observation_fusion_level_matrix[chr - 1] = estimate_fusion_level(hi_c_mat, hi_c_data_model_list, del_mat[chr - 1], RESOLUTION, window_size, contact_threshold_list, our_observation_fusion_level_matrix_list[chr - 1]);
		for (int t = 0; t < contact_threshold_list.size(); t++)
			our_total_observation_fusion_level[t] += our_observation_fusion_level_matrix[chr - 1][t];
		for (int t = 0; t < fusion_threshold_list.size(); t++) {
			our_observation_fusion_count_matrix[chr - 1][t] = 0;
			for (int del = 0; del < our_observation_fusion_level_matrix_list[chr - 1][final_contact_threshold_index].size(); del++)
				if (our_observation_fusion_level_matrix_list[chr - 1][final_contact_threshold_index][del] >= fusion_threshold_list[t])
					our_observation_fusion_count_matrix[chr - 1][t]++;
			our_total_observation_fusion_count[t] += our_observation_fusion_count_matrix[chr - 1][t];
		}
		// Permutated
		ValueList tmp;
		for (int sample = 0; sample < sample_num; sample++) {
			if (sample % 100 == 0)
				cout << sample << endl;
			DeletionList sampling_del_list;
			sampling_deletion (del_mat[chr - 1], hi_c_mat.size()*RESOLUTION, window_size*RESOLUTION, sampling_del_list);
			for (int tool = 0; tool < tool_name_list.size(); tool++) {
				double sample_fusion_level_tmp;
				if (tool_type_list[tool] < 1)
					sample_fusion_level_tmp = conventional_count_fusion (sampling_del_list, TAD_mat_list[tool][chr - 1], tmp, false);
				else 
					sample_fusion_level_tmp = conventional_count_fusion (sampling_del_list, TAD_boundary_mat_list[tool][chr - 1], tmp, false);
				other_permutation_fusion_dist_matrix[tool][chr - 1].add_value(sample_fusion_level_tmp);
				// total
				other_total_permutation_fusion_level[tool][sample] += sample_fusion_level_tmp;
			}
			ValueMatrix mat_tmp;
			ValueList sample_fusion_level_list_tmp = estimate_fusion_level(hi_c_mat, hi_c_data_model_list, sampling_del_list, RESOLUTION, window_size, contact_threshold_list, mat_tmp);
			for (int t = 0; t < contact_threshold_list.size(); t++) {
				our_permutation_fusion_level_dist_matrix[chr - 1][t].add_value(sample_fusion_level_list_tmp[t]);
				// total
				our_total_permutation_fusion_level[sample][t] += sample_fusion_level_list_tmp[t];
			}
			for (int t = 0; t < fusion_threshold_list.size(); t++) {
				int our_permutation_fusion_count_tmp = 0;
				for (int del = 0; del < mat_tmp[final_contact_threshold_index].size(); del++)
					if (mat_tmp[final_contact_threshold_index][del] >= fusion_threshold_list[t])
						our_permutation_fusion_count_tmp++;
				our_permutation_fusion_count_dist_matrix[chr - 1][t].add_value(our_permutation_fusion_count_tmp);
				// total
				our_total_permutation_fusion_count[sample][t] += our_permutation_fusion_count_tmp;
			}
		}
	}

	ValueMatrix CTCF_del_matrix(CHR_NUM);
	ValueList CTCF_del_num_list(CHR_NUM, 0);
	ValueList fusion_count_list(tool_name_list.size(), 0);
	//ValueList our_fusion_count_list(contact_threshold_list.size(), 0);
	DeletionInfoList del_info_list;
	ValueMatrix our_all_fusion_level_matrix(contact_threshold_list.size());
	double LOOCV_cut_off = 0.7;
	// Write the prediction of all tools	
	string filename_tmp = "Debug/" + exp_name + "_fusion_summary.dat";
	ofstream out_observation_fusion_summary_file(filename_tmp.c_str());
	out_observation_fusion_summary_file << "chr\tdel_left\tdel_right\tCTCF";
	for (int tool = 0; tool < tool_name_list.size(); tool++)
		out_observation_fusion_summary_file << "\t" << tool_name_list[tool];
	for (int t = 0; t < contact_threshold_list.size(); t++)
		out_observation_fusion_summary_file << "\t" << "New (" << contact_threshold_list[t] << ")";
	out_observation_fusion_summary_file << endl;
	for (int chr = chr_min; chr <= chr_max; chr++) {
		string chr_name = ((chr < 23)? num_to_string(chr):"X");
		resize_and_fill(CTCF_del_matrix[chr - 1], del_mat[chr - 1].size(), 0);
		for (int del = 0; del < del_mat[chr - 1].size(); del++) {
			for (int ctcf = 0; ctcf < CTCF_mat[chr - 1].size(); ctcf++)
				if (CTCF_mat[chr - 1][ctcf].left >= del_mat[chr - 1][del].left && CTCF_mat[chr - 1][ctcf].right <= del_mat[chr - 1][del].right) {
					CTCF_del_matrix[chr - 1][del]++;					
				}
			CTCF_del_num_list[chr - 1] += CTCF_del_matrix[chr - 1][del];
			out_observation_fusion_summary_file << chr_name << "\t" << del_mat[chr - 1][del].left << "\t" << del_mat[chr - 1][del].right << "\t" << CTCF_del_matrix[chr - 1][del];
			for (int tool = 0; tool < tool_name_list.size(); tool++) {
				out_observation_fusion_summary_file << "\t" << other_observation_fusion_level_matrix_list[chr - 1][tool][del];
				if (other_observation_fusion_level_matrix_list[chr - 1][tool][del] >= LOOCV_cut_off)
					fusion_count_list[tool]++;
			}
			for (int t = 0; t < contact_threshold_list.size(); t++) {
				out_observation_fusion_summary_file << "\t" << our_observation_fusion_level_matrix_list[chr - 1][t][del];
				//if (our_One_KG_fusion_level_matrix_list[chr - 1][t][del] >= LOOCV_cut_off)
				//	our_fusion_count_list[t]++;
				our_all_fusion_level_matrix[t].push_back(our_observation_fusion_level_matrix_list[chr - 1][t][del]);
			}
			out_observation_fusion_summary_file << endl;
			// prepare for LOOCV
			del_info_list.push_back(DeletionInfo(chr, del));
		}
	}
	out_observation_fusion_summary_file.close();

	// LOOCV
	ValueList rank_tmp;
	ValueMatrix sorted_index_mat (contact_threshold_list.size());
	for (int t = 0; t < contact_threshold_list.size(); t++)
		sort_by_value (our_all_fusion_level_matrix[t], sorted_index_mat[t], rank_tmp, false);
	//ValueList concensus_fusion_num(tool_name_list.size(), 0), 
	//	other_CTCF_fusion(tool_name_list.size(), 0), our_CTCF_fusion(tool_name_list.size(), 0), 
	//	other_tp(tool_name_list.size(), 0), our_tp(tool_name_list.size(), 0);
	filename_tmp = "Debug/" + exp_name + "_LOOCV.dat";
	ofstream loocv_file(filename_tmp.c_str());
	loocv_file << "#Tool\tConcensus\t#fusion\tCTCF_del\tCTCF_tot\tAgree";
	for (int t = 0; t < contact_threshold_list.size(); t++)
		loocv_file << "\tCTCF_del(" << contact_threshold_list[t] << ")\tCTCF_tot(" << contact_threshold_list[t] << ")\tAgree(" << contact_threshold_list[t] << ")";
	loocv_file << endl;
	for (int tool = 0; tool < tool_name_list.size(); tool++) {
		int concensus_fusion_num = 0, other_CTCF_fusion_del = 0, other_CTCF_fusion_total = 0, other_tp = 0;
		ValueMatrix concensus_fusion_matrix(CHR_NUM);
		for (int chr = chr_min; chr <= chr_max; chr++)
			resize_and_fill(concensus_fusion_matrix[chr - 1], del_mat[chr - 1].size(), 0);
		for (int chr = chr_min; chr <= chr_max; chr++) {
			for (int del = 0; del < del_mat[chr - 1].size(); del++) {
				for (int another_tool = 0; another_tool < tool_name_list.size(); another_tool++) {
					if (another_tool != tool)
						concensus_fusion_matrix[chr - 1][del] += other_observation_fusion_level_matrix_list[chr - 1][another_tool][del];
				}
				if (concensus_fusion_matrix[chr - 1][del] >= 1 + LOOCV_cut_off) {
					concensus_fusion_num++;
					if (other_observation_fusion_level_matrix_list[chr - 1][tool][del] >= LOOCV_cut_off)
						other_tp++;				
				}
				if (other_observation_fusion_level_matrix_list[chr - 1][tool][del] >= LOOCV_cut_off) {
					if (CTCF_del_matrix[chr - 1][del] > 0)
						other_CTCF_fusion_del++;
					other_CTCF_fusion_total += CTCF_del_matrix[chr - 1][del];
				}
			}
		}
		loocv_file << tool_name_list[tool] << "\t" << concensus_fusion_num << "\t" << fusion_count_list[tool] << "\t"
			<< other_CTCF_fusion_del << "\t" << other_CTCF_fusion_total << "\t" << other_tp;
		for (int t = 0; t < contact_threshold_list.size(); t++) {
			// Find the comparable fusion list from our tool
			int our_tp = 0, our_CTCF_fusion_del = 0, our_CTCF_fusion_total = 0;
			ValueMatrix our_fusion_matrix(CHR_NUM);
			for (int chr = chr_min; chr <= chr_max; chr++) 				
				resize_and_fill(our_fusion_matrix[chr - 1], del_mat[chr - 1].size(), 0);
			for (int k = 0; k < fusion_count_list[tool]; k++) 
				our_fusion_matrix[del_info_list[round(sorted_index_mat[t][k])].chr - 1][del_info_list[round(sorted_index_mat[t][k])].del] = 1;
			for (int chr = chr_min; chr <= chr_max; chr++) {
				for (int del = 0; del < del_mat[chr - 1].size(); del++) {
					if (our_fusion_matrix[chr - 1][del] >= LOOCV_cut_off) {
						if (concensus_fusion_matrix[chr - 1][del] >= 1 + LOOCV_cut_off) 
							our_tp++;
						if (CTCF_del_matrix[chr - 1][del] > 0)
							our_CTCF_fusion_del++;
						our_CTCF_fusion_total += CTCF_del_matrix[chr - 1][del];
					}
				}
			}
		 	loocv_file << "\t" << our_CTCF_fusion_del << "\t" << our_CTCF_fusion_total << "\t" << our_tp;
		}
		loocv_file << endl;
	}
	// p-value
	filename_tmp = "Debug/" + exp_name + "_p_value.dat";
	ofstream p_value_file(filename_tmp.c_str());
	for (int tool = 0; tool < tool_name_list.size(); tool++) {
		DistributionSummary tmp;
		for (int sample = 0; sample < sample_num; sample++)
			tmp.add_value(other_total_permutation_fusion_level[tool][sample]);
		p_value_file << tool_name_list[tool] << "\t" << other_total_observation_fusion_level[tool] << "\t"
			<< round(tmp.mean*100)/100 << " +/- " << round(tmp.std*100)/100 << "\t"  
			<< ((tmp.std > 0)? round(100*(other_total_observation_fusion_level[tool] - tmp.mean)/tmp.std)/100 : 0) << endl;
	}
	for (int t = 0; t < contact_threshold_list.size(); t++) {
		DistributionSummary tmp;
		for (int sample = 0; sample < sample_num; sample++)
			tmp.add_value(our_total_permutation_fusion_level[sample][t]);
		p_value_file << "Threshold = " << contact_threshold_list[t] << "\t" << our_total_observation_fusion_level[t] << "\t"
			<< round(tmp.mean) << " +/- " << round(tmp.std) << "\t"  
			<< ((tmp.std > 0)? round(100*(our_total_observation_fusion_level[t] - tmp.mean)/tmp.std)/100 : 0) << endl;
	}
	for (int t = 0; t < fusion_threshold_list.size(); t++) {
		DistributionSummary tmp;
		for (int sample = 0; sample < sample_num; sample++)
			tmp.add_value(our_total_permutation_fusion_count[sample][t]);
		p_value_file << "Threshold = " << fusion_threshold_list[t] << "\t" << our_total_observation_fusion_count[t] << "\t"
			<< round(tmp.mean) << " +/- " << round(tmp.std) << "\t"  
			<< ((tmp.std > 0)? round(100*(our_total_observation_fusion_count[t] - tmp.mean)/tmp.std)/100 : 0) << endl;	
	}
	p_value_file.close();
	// print the distribution
	filename_tmp = "Debug/" + exp_name + "_distribution.dat";
	ofstream dist_file(filename_tmp.c_str());
	dist_file << "#";
	for (int sample = -1; sample < sample_num; sample++) {
		for (int tool = 0; tool < tool_name_list.size(); tool++) {
			if (tool > 0)
				dist_file << "\t";
			if (sample < 0)
				dist_file << tool_name_list[tool];
			else
				dist_file << other_total_permutation_fusion_level[tool][sample];
		}
		for (int t = 0; t < contact_threshold_list.size(); t++) {
			if (sample < 0) 
				dist_file << "\t" << "Threshold = " << contact_threshold_list[t];
			else
				dist_file << "\t" << our_total_permutation_fusion_level[sample][t];
		}
		for (int t = 0; t < fusion_threshold_list.size(); t++) {
			if (sample < 0)
				dist_file << "\t" << "Threshold = " << fusion_threshold_list[t];
			else
				dist_file << "\t" << our_total_permutation_fusion_count[sample][t];
		}
		dist_file << endl;
	}
	dist_file.close();
}

void final_prediction_experiment() {
	IdList window_size_list = {25, 50, 100, 200, 400};
	prediction_experiment(21, 1, 100000, window_size_list, false);
}

void final_prediction_plot() {
	IdList window_size_list = {25, 50, 100};
	prediction_experiment(21, 4200, 4200, window_size_list, true);
}

void final_1KG_fusion_experiment(int sample_num) {
	ValueList contact_threshold_list;
	for (int i = 5; i <= 15; i++)
		contact_threshold_list.push_back(i);
	ValueList fusion_threshold_list;
	for (int i = 350; i <= 750; i = i + 50)
		fusion_threshold_list.push_back(i);
	string del_filename = "/share/hormozdiarilab/Data/ValidatedVariants/1000G_SV_Phase3/ALL.Del.Bed";
	fusion_experiment("GM12878_1KG", del_filename, true, 2500000, contact_threshold_list, 2, fusion_threshold_list, sample_num);
}

void final_disease_fusion_experiment(int sample_num) {
	ValueList contact_threshold_list;
	for (int i = 5; i <= 15; i++)
		contact_threshold_list.push_back(i);
	ValueList fusion_threshold_list;
	for (int i = 350; i <= 750; i = i + 50)
		fusion_threshold_list.push_back(i);
	string del_filename = "/share/hormozdiarilab/Data/ValidatedVariants/Disease/deletion_mutation.dat";
	fusion_experiment("GM12878_disease", del_filename, false, 2500000, contact_threshold_list, 2, fusion_threshold_list, sample_num);
}

int main() {
	//export_all_hi_c_model();
	//final_prediction_experiment();
	//final_prediction_plot();
	final_1KG_fusion_experiment(10000);
	//final_disease_fusion_experiment(100);
	cout << "Complete" << endl;
	return 1;
}
