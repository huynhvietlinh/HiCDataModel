//#: Date         :05/05/2017
//#: Author       :Linh Huynh
//#: Version      :1.0.0 

#include <ilcplex/cplex.h>
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

#define STR_EQ 		0
#define ROW_DELIM 	'\n'
#define FIELD_DELIM 	'\t'
#define COMMA_DELIM 	','

#define CHR_NUM				23
#define RESOLUTION			5000

#define ZERO				1e-1

#define ARROW_HEAD			1
#define CATCH				2
#define HICSEG				3
#define INSULATION_SCORE		4
#define ARMATUS				5

#define DELETION_FILE			0
#define DELETION_1KG_FILE 		1
#define CTCF_FILE			2
#define RAO_LOOP_FILE			3
#define JUICE_FILE			4
#define INSULATION_SCORE_FILE 		5
#define TAD_FUSION_SCORE_FILE		6

#define MAX_LENGTH_OBJ			100
#define FUSION_SCORE_WINDOW_SIZE	50

#define INSULATION_SCORE_THRESHOLD	0.7

struct ExplicitDistribution {
	ExplicitDistribution();	
	void add_value (int);
	double get_prob (int) const;
	void read_from_file (string);
	void print_to_file (string) const;
 	
	ValueList freq;
	double n = 0;
};

struct DistributionSummary {
	DistributionSummary();
	void add_value (double);
	void merge (const DistributionSummary&);

	int size;
	double mean, std;
};
typedef vector<DistributionSummary> DistributionSummaryList;
typedef vector<DistributionSummaryList> DistributionSummaryMatrix;

struct Sequence {
	Sequence ();
	Sequence (int l, int r, double s = 0);
	int left, right;
	ValueList score_list;
};
typedef vector<Sequence> SequenceList;
typedef vector<SequenceList> SequenceMatrix;

struct SequencePair {
	Sequence src, dest;
};
typedef vector<SequencePair> SequencePairList;
typedef vector<SequencePairList> SequencePairMatrix;
typedef vector<SequencePairList> SequencePairCube;

struct HiCDataModel {
	int get_right() const;
	double get_alpha(int bin_id) const;

	ValueList alpha;
	double beta;
	ValueList resistance;
	int left;		// bin, not base-pair
};
typedef vector<HiCDataModel> HiCDataModelList;

struct Genomic_Interaction_Distribution {
	void load(string filename, int round_factor_);
	double get_prob(double contact_freq) const;

	int round_factor;
	ValueList freq_distribution;
};

// Auxiliari functions
string num_to_string (int);
double max(double, double);
double min(double, double);
//
double PCC (const ValueList&, const ValueList&, double);  		// pred, obs, zero
double average_log_ratio (const ValueList&, const ValueList&, double);	
double L1(const ValueList&, const ValueList&);
double L2(const ValueList&, const ValueList&);
void line_LQ (const ValueList&, const ValueList&, double&, double&);	//x, y, a, b: y = ax + b
//
void resize_and_fill(ValueList&, int, int);				// size, value
void read_a_table_file (string filename, char row_delim_char, char field_delim_char, 
	char comment_char, int skipping_header_line_num, 
	const IdList& col_list, StringMatrix& str_mat);			// index of 1st col = 0
void read_a_value_list (string filename, int col, ValueList& val_list);	// index of 1st col = 0
void write_a_value_list (string filename, ValueList&);
void sort_by_value (const ValueList& l, ValueList& sorted_index, ValueList& rank, bool is_asc);
//
void Dijkstra (const ValueMatrix& hi_c_mat, int left, int right, int src, ValueList& shortest_path_val_list);
void Floyd_APSP (const ValueMatrix& graph_dist_mat, ValueMatrix& shortest_dist_mat);

// Hi-C processing functions
void read_sequence_list_file (string filename, int file_format, int min_length, int max_length, SequenceMatrix& seq_mat);
int count(const SequenceMatrix&); 
void deletion_filter (SequenceMatrix& original_del_mat, const SequenceMatrix& filter_del_mat);
void extract_1KG_without_GM12878 (string One_KG_deletion_filename, string GM12878_variant_filename, int min_length, int max_length);
//
void read_hi_c_data (string filename, int resolution, ValueMatrix& hi_c_matrix);
void read_5C_file (string filename, int header_line_num, int matrix_size,  ValueMatrix& hi_c_matrix);
//
void write_hi_c_data_model_list(string, const HiCDataModelList&);
void read_hi_c_data_model_list(string, HiCDataModelList&);
void export_hi_c_model_list (const ValueMatrix& hi_c_mat, int window_size, int shift_length, int max_length, string filename);
void export_all_hi_c_model (int chr_min, int chr_max, int max_length_obj, int window_length, int shift_length, string model_folder); 
// 
void write_TAD_fusion_score_primary (const SequenceMatrix&, string filename);
void write_TAD_fusion_score_permutation( const vector<SequenceMatrix>&, string filename);
void read_TAD_fusion_score_permutation (string filename, const SequenceMatrix&, vector<SequenceMatrix>&);
double find_fusion_score_sum(const SequenceMatrix& seq_mat, int score_col_index);

// Hi-C main functions
ValueList length_based_model_prediction (const ValueMatrix& hi_c_mat_wt, const ValueMatrix& hi_c_mat_mut, int left_1, int right_1, int left_2, int right_2, string print_filename = "");
//
void fit_hi_c_data_model (const ValueMatrix& hi_c_mat, int left, int right, int max_length_in_obj, HiCDataModel& hi_c_data_model, bool print_solver_status = false);
//
ValueList predict_HiC (const ValueMatrix& hi_c_mat_wt, const ValueMatrix& hi_c_mat_mut, const HiCDataModel& hi_c_data_model, const Genomic_Interaction_Distribution& genomic_interaction_dist, int left_1, int right_1, int left_2, int right_2, string print_prefix = "");
//
double predict_HiC (const ValueMatrix& hi_c_mat, const HiCDataModelList& hi_c_data_model_list, const Genomic_Interaction_Distribution& genomic_interaction_dist, int left_1, int right_1, int left_2, int right_2, string prefix_print = "");
//
ValueList predict_HiC (const ValueMatrix& hi_c_mat_wt, const ValueMatrix& hi_c_mat_mut, int left_1, int right_1, int left_2, int right_2, int max_length, const Genomic_Interaction_Distribution& genomic_interaction_dist, string print_prefix = "");
//
void print_distribution(int max_loop_length, int resolution);
//
void conventional_count_fusion (SequenceList& del_list, const SequenceList& tad_list, int score_index, bool is_printed = false); 		// overlap or not
void conventional_estimate_fusion (SequenceList& del_list, const SequenceList& tad_boundary_list, int score_index, bool is_printed = false);	// find the max boundary score
//
void generate_TAD_fusion_score(const SequenceMatrix& del_mat, string hi_c_folder, string model_folder, string TAD_annotation_folder, string CTCF_filename, const Genomic_Interaction_Distribution& genomic_interaction_dist, int window_size, int sample_num, vector<SequenceMatrix>& result); // sample_num <= 0: find the score of primary deletions

// Experiments
void fig_S1_robustness_experiment();
void print_prediction_case_study(int chr, int del_left, int del_right, int window_size, const Genomic_Interaction_Distribution& genomic_interaction_distribution, bool is_with_K562);
void HoxD_case_study(const Genomic_Interaction_Distribution& genomic_interaction_dist, int window_size);
void Firre_case_study(const Genomic_Interaction_Distribution& genomic_interaction_dist, int window_size);
void K562_prediction(string del_filename, const Genomic_Interaction_Distribution& genomic_interaction_dist, string out_filename);
void generate_TAD_fusion_score(string del_filename, int min_length, int max_length, string hi_c_folder, string model_folder, string TAD_annotation_folder, string CTCF_filename, const Genomic_Interaction_Distribution& genomic_interaction_dist, int window_size, int sample_num, string print_prefix); // sample_num <= 0: find the score of primary deletions
void tilling_experiment(string del_filename, string hi_c_folder, string model_folder, string TAD_annotation_folder, string CTCF_filename, const Genomic_Interaction_Distribution& genomic_interaction_dist, int window_size, string output_prefix);
double find_top_TAD_fusion_score_threshold(const string& primary_filename, int top);	// return a score at the top ranking position
void analyze_TAD_fusion_score (const string& primary_filename, const StringList& permutation_filename_list = {}, const ValueList& cut_off_list = {10}, double score_threshold = INSULATION_SCORE_THRESHOLD);
void compare_TAD_fusion_score_by_binning(string filename_1, string filename_2, double score_threshold, int bin_num, int bin_size, string out_file);

int main() {
	string CTCF_binding_filename = "/share/hormozdiarilab/Data/CTCF_Binding/UW_hg19_GM12878_CTCFBSDB.bed";
	string TAD_annotation_folder = "/share/hormozdiarilab/Data/HiC/TAD_Annotation/GM12878/";
	//string genomic_interaction_distribution_filename = "/share/hormozdiarilab/Data/HiC/Metric/top_restrict_contact_prob_chr_1_1000000_1000000.dat";
	//string genomic_interaction_distribution_filename = "/share/hormozdiarilab/Data/HiC/Metric/Rao_loop_contact_prob_chr_1_1_max_length_1000000_rounded_500.dat";
	string genomic_interaction_distribution_filename = "/share/hormozdiarilab/Data/HiC/Metric/Rao_loop_contact_prob_chr_1_23_max_length_1000000_rounded_500.dat";

	string GM12878_hi_c_intra_chr_data_folder = "/share/hormozdiarilab/Data/HiC/GM12878_combined/5kb_resolution_intrachromosomal/chr",
		K562_hi_c_intra_chr_data_folder = "",
		GM12878_ReModel_folder = "/share/hormozdiarilab/Data/HiC/Model/GM12878/ReModel/";
	string K562_del_filename = "/share/hormozdiarilab/Data/RECOMB_2018/old_K562_del.dat";

	//string Original_1KG_deletion_filename = "/share/hormozdiarilab/Data/ValidatedVariants/1000G_SV_Phase3/ALL.Del.Bed",
	//	GM12878_deletion_filename = "/share/hormozdiarilab/Data/ValidatedVariants/1000G_SV_Phase3/NA12878_Del.BED";

	string One_KG_deletion_without_GM12878_filename = "/share/hormozdiarilab/Data/ValidatedVariants/1000G_SV_Phase3/1KG_without_GM12878_10000_10000000.dat";
	string diseased_deletion_filename = "/share/hormozdiarilab/Data/RECOMB_2018/Diseased_del.dat";
	string GreatApe_deletion_filename = "/share/hormozdiarilab/Data/RECOMB_2018/GreatApe/GreatApe.Del.Hg19";

	string TCGA_ALL_filename = "/share/hormozdiarilab/Data/TCGA/TCGA_del.dat";
	string TCGA_LAML_filename = "/share/hormozdiarilab/Data/TCGA/TCGA_LAML_CNV.dat";
 	
	string DD_filename = "/share/hormozdiarilab/Data/ValidatedVariants/CNV_Morbidity_Map/Signature_Del_hg37.RemoveSegDups.Longer5kbp.Unique.bed";

	cout << "MAX_OBJ_LENGTH = " << MAX_LENGTH_OBJ << endl;
	cout << "FUSION_SCORE_WINDOW_SIZE = " << FUSION_SCORE_WINDOW_SIZE << endl;

	fig_S1_robustness_experiment();

	//export_all_hi_c_model (1,23);
	//print_distribution(1000000, 5000);

	Genomic_Interaction_Distribution genomic_interaction_dist;
	genomic_interaction_dist.load(genomic_interaction_distribution_filename, 500);

	//HoxD_case_study(genomic_interaction_dist, FUSION_SCORE_WINDOW_SIZE);
	//Firre_case_study(genomic_interaction_dist, FUSION_SCORE_WINDOW_SIZE);
	//K562_prediction(K562_del_filename, genomic_interaction_dist, "Debug/K562_stat.dat");

	//generate_TAD_fusion_score(One_KG_deletion_without_GM12878_filename, 1e4, 1e7, 
	//	GM12878_hi_c_intra_chr_data_folder, GM12878_ReModel_folder,	// old model 
	//	TAD_annotation_folder, CTCF_binding_filename, genomic_interaction_dist,
	//	FUSION_SCORE_WINDOW_SIZE, 0,"Debug/1KG_FinalModel");
	//generate_TAD_fusion_score(One_KG_deletion_without_GM12878_filename, 1e4, 1e7, 
	//	GM12878_hi_c_intra_chr_data_folder, GM12878_ReModel_folder,	// old model 
	//	TAD_annotation_folder, CTCF_binding_filename, genomic_interaction_dist,
	//	FUSION_SCORE_WINDOW_SIZE, 1500,"Debug/1KG_FinalModel");

	//generate_TAD_fusion_score(GreatApe_deletion_filename, 1e4, 1e7, 
	//	GM12878_hi_c_intra_chr_data_folder, GM12878_ReModel_folder,	// old model 
	//	TAD_annotation_folder, CTCF_binding_filename, genomic_interaction_dist,
	//	FUSION_SCORE_WINDOW_SIZE, 0,"Debug/GreatApe_FinalModel");
	//generate_TAD_fusion_score(GreatApe_deletion_filename, 1e4, 1e7, 
	//	GM12878_hi_c_intra_chr_data_folder, GM12878_ReModel_folder,	// old model 
	//	TAD_annotation_folder, CTCF_binding_filename, genomic_interaction_dist,
	//	FUSION_SCORE_WINDOW_SIZE, 500,"Debug/GreatApe_FinalModel");

	//generate_TAD_fusion_score(diseased_deletion_filename, 1e4, 1e7, 
	//	GM12878_hi_c_intra_chr_data_folder, GM12878_ReModel_folder, // old model 
	//	TAD_annotation_folder, CTCF_binding_filename, genomic_interaction_dist,
	//	FUSION_SCORE_WINDOW_SIZE, 0, "Debug/diseased_FinalModel");

	//tilling_experiment(diseased_deletion_filename, GM12878_hi_c_intra_chr_data_folder, GM12878_ReModel_folder, 
	//	TAD_annotation_folder, CTCF_binding_filename, genomic_interaction_dist, FUSION_SCORE_WINDOW_SIZE, "Debug/diseased_FinalModel");

	//generate_TAD_fusion_score(TCGA_ALL_filename, 1e4, 5e5,
	//	GM12878_hi_c_intra_chr_data_folder, GM12878_ReModel_folder, 
	//	TAD_annotation_folder, CTCF_binding_filename, genomic_interaction_dist,
	//	FUSION_SCORE_WINDOW_SIZE, 0,"Debug/TCGA_ALL_FinalModel");	// old model
	//generate_TAD_fusion_score(DD_filename, 1e4, 5e5,
	//	GM12878_hi_c_intra_chr_data_folder, GM12878_ReModel_folder, 
	//	TAD_annotation_folder, CTCF_binding_filename, genomic_interaction_dist,
	//	FUSION_SCORE_WINDOW_SIZE, 0,"Debug/DD_ReModel");	// old model

	//cout << "New top 100 " << find_top_TAD_fusion_score_threshold("Debug/1KG_new_model_primary.dat", 100) << endl; 
	//cout << "Old top 100 " << find_top_TAD_fusion_score_threshold("Debug/1KG_old_model_primary.dat", 100) << endl;
	
	//analyze_TAD_fusion_score("Debug/1KG_FinalModel_50_primary.dat");
	//analyze_TAD_fusion_score("Debug/1KG_FinalModel_50_primary.dat", 
	//	{"Permutation/1KG_FinalModel_50_permutation_5000.dat"}, 
	//		{10, 50, 100, 200, 300, 400, 500, 1000, 1500, 2000});
	//analyze_TAD_fusion_score("Debug/GreatApe_FinalModel_50_primary.dat", 
	//	{"Permutation/GreatApe_FinalModel_50_permutation_5000.dat"}, 
	//		{6, 9, 12, 15, 100, 500});

	//compare_TAD_fusion_score_by_binning("Debug/1KG_FinalModel_50_primary.dat", "Debug/TCGA_ALL_FinalModel_50_primary.dat", 
	//	find_top_TAD_fusion_score_threshold("Debug/1KG_FinalModel_50_primary.dat", 100),
	//	10, 50000, "Debug/1KG_vs_TCGA_ALL_FinalModel");
	
	//compare_TAD_fusion_score_by_binning("Debug/1KG_ReModel_50_primary.dat", "Debug/DD_ReModel_50_primary.dat", 
	//	find_top_TAD_fusion_score_threshold("Debug/1KG_ReModel_50_primary.dat", 100),
	//	10, 50000, "Debug/1KG_vs_DD_ReModel");

	//compare_TAD_fusion_score_by_binning("Debug/1KG_old_model_primary.dat", "Debug/GreatApe_old_model_primary.dat", 
	//	find_top_TAD_fusion_score_threshold("Debug/1KG_old_model_primary.dat", 100),
	//	10, 50000, "Debug/binning_1KG_vs_GreatApe_old_model.dat");

	cout << "Complete" << endl;
	return 1;
}

void tilling_experiment(string del_filename, string hi_c_folder, string model_folder, string tad_annotation_folder, string CTCF_filename, const Genomic_Interaction_Distribution& genomic_interaction_dist, int window_size, string output_prefix) {
	SequenceMatrix del_mat;
	read_sequence_list_file (del_filename, DELETION_FILE, 1e4, 1e7, del_mat);
	// Sampling deletions

	SequenceMatrix tilling_del_mat(del_mat.size());
	SequenceMatrix del_index_mat = del_mat;	// each of them contains the low bound index and the up bound index of a deletion

	int max_index = 0;
	for (int chr = 0; chr < del_mat.size(); chr++) {
		for (int del = 0; del < del_mat[chr].size(); del++) {
			del_index_mat[chr][del].left = tilling_del_mat[chr].size();			
			int current_left = del_mat[chr][del].left;
			while (current_left < del_mat[chr][del].right) {
				Sequence seq(current_left, current_left + 20000);
				tilling_del_mat[chr].push_back(seq);
				current_left += 1e4;
			}
			del_index_mat[chr][del].right = tilling_del_mat[chr].size() - 1;
			max_index = max(max_index, del_index_mat[chr][del].right - del_index_mat[chr][del].left + 1);
		}
	}
	vector<SequenceMatrix> result_del_cube;
	generate_TAD_fusion_score (tilling_del_mat, hi_c_folder, model_folder, tad_annotation_folder, CTCF_filename, genomic_interaction_dist, window_size, 0, result_del_cube);

	write_TAD_fusion_score_primary (result_del_cube[0], "Debug/tilling_tmp");

	int del_num = count(del_mat);
	ValueMatrix tilling_del_score_mat(del_num), CTCF_del_score_mat(del_num), without_CTCF_del_score_mat(del_num);
	int del_index = 0;
	for (int chr = 0; chr < del_mat.size(); chr++) {
		for (int del = 0; del < del_mat[chr].size(); del++) {
			for (int k = del_index_mat[chr][del].left; k <= del_index_mat[chr][del].right; k++) {
				double fusion_score = result_del_cube[0][chr][k].score_list[6];
				tilling_del_score_mat[del_index].push_back(fusion_score);	
				if (result_del_cube[0][chr][k].score_list[5] > 0)
					CTCF_del_score_mat[del_index].push_back(fusion_score);
				else
					without_CTCF_del_score_mat[del_index].push_back(fusion_score);
			}
			del_index++;
		}			
	}

	string tilling_filename = output_prefix + "_tilling_all.dat";
	ofstream tilling_file(tilling_filename.c_str());
	tilling_file << "#";
	for (int i = 0; i < max_index; i++) {
		tilling_file << endl;
		for (int j = 0; j < del_num; j++) {
			if (i < tilling_del_score_mat[j].size())
				tilling_file << tilling_del_score_mat[j][i];
			else
				tilling_file << ".";
			tilling_file << "\t";
		}
	}

	tilling_file.close();

	string out_filename = output_prefix + "_tilling_vs_CTCF.dat";
	ofstream out_file(out_filename.c_str());
	out_file << "#";
	for (int i = 0; i < max_index; i++) {
		out_file << endl;
		for (int j = 0; j < del_num; j++) {
			if (i < CTCF_del_score_mat[j].size())
				out_file << CTCF_del_score_mat[j][i];
			else
				out_file << ".";
			out_file << "\t";
		}
		for (int j = 0; j < del_num; j++) {
			if (i < without_CTCF_del_score_mat[j].size())
				out_file << without_CTCF_del_score_mat[j][i];
			else
				out_file << ".";
			out_file << "\t";
		}
	}
	out_file.close();
}

void HoxD_case_study(const Genomic_Interaction_Distribution& genomic_interaction_dist, int fusion_score_window_size) {
	ValueMatrix wt_mat, del_1_13_mat, del_attP_Rel5_mat;
	read_5C_file("/share/hormozdiarilab/Data/HiC/HoxD_Gene_Development_2017/PL_HiC_E12_Wt__chr2__20kb__raw.matrix", 6, 9106, wt_mat);
	read_5C_file("/share/hormozdiarilab/Data/HiC/HoxD_Gene_Development_2017/PL_HiC_E12_del1-13d9lac__chr2__20kb__raw.matrix_add_zero_lines.dat", 6, 9106, del_1_13_mat);	// origin: 9102
	read_5C_file("/share/hormozdiarilab/Data/HiC/HoxD_Gene_Development_2017/PL_HiC_E12_delattP-Rel5d9lac__chr2__20kb__raw.matrix_add_zero_lines.dat", 6, 9106, del_attP_Rel5_mat); // origin: 9089
	
	//for (int i = 500; i < 1000; i++)
	//	cout << wt_mat[i][i] << endl;
	// HoxD: 73,500,000 - 76,000,000
	// del_1_13_d9lac: 74,663,890 - 74,764,314
	int resolution = 20000;
	int left = round(73500000/resolution), right = round(76000000/resolution);
	int del_1_left = floor(74663890/resolution), del_1_right = floor(74764314/resolution);
	ValueList result = predict_HiC(wt_mat, del_1_13_mat, left, del_1_left, del_1_right, right, MAX_LENGTH_OBJ, genomic_interaction_dist, "Debug/Del_1_13_full");
	cout << "PCC = " << result[0] << ", L1 = " << result[2] << endl;
	// del_attP_Rel5_d9lac: 74,422,050 - 74,768,746
	int del_2_left = floor(74422050/resolution), del_2_right = floor(74768746/resolution) + 1;
	result = predict_HiC(wt_mat, del_attP_Rel5_mat, left, del_2_left, del_2_right, right, MAX_LENGTH_OBJ, genomic_interaction_dist, "Debug/Del_attP_Rel5_full");
	cout << "PCC = " << result[0] << ", L1 = " << result[2] << endl;

	// For scoring
	result = predict_HiC(wt_mat, del_1_13_mat, del_1_left - fusion_score_window_size, del_1_left, del_1_right, del_1_right + fusion_score_window_size, 
		MAX_LENGTH_OBJ, genomic_interaction_dist, "Debug/Del_1_13_score_only");
	cout << "Del_1_13, fusion score = " << result[4] << endl;

	result = predict_HiC(wt_mat, del_attP_Rel5_mat, del_2_left - fusion_score_window_size, del_2_left, del_2_right, del_2_right + fusion_score_window_size, 
		MAX_LENGTH_OBJ, genomic_interaction_dist, "Debug/Del_attP_Rel5_score_only");
	cout << "Del_attP_Rel5, fusion score = " << result[4] << endl;
}

void Firre_case_study(const Genomic_Interaction_Distribution& genomic_interaction_dist, int fusion_score_window_size) {
	ValueMatrix wt_mat_10kb, Firre_del_mat_10kb, wt_mat_40kb, Firre_del_mat_40kb; 

	read_5C_file("/share/hormozdiarilab/Data/HiC/Firre_Nat_Comm/Female_MEF_FirreWT_mm9_40000_iced_chrX_dense.addedHeaders.matrix", 19, 4167, wt_mat_40kb);
	read_5C_file("/share/hormozdiarilab/Data/HiC/Firre_Nat_Comm/Female_MEF_FirreKO_mm9_40000_iced_chrX_dense.addedHeaders.matrix", 19, 4167, Firre_del_mat_40kb);
	//read_5C_file("/share/hormozdiarilab/Data/HiC/Firre_Nat_Comm/Combined_femaleWT_10000_iced_chrX_dense.addedHeaders.matrix", 19, 16666, wt_mat_10kb);
	//read_5C_file("/share/hormozdiarilab/Data/HiC/Firre_Nat_Comm/Combined_femaleKO_10000_iced_chrX_dense.addedHeaders.matrix", 19, 16666, Firre_del_mat_10kb);

	//for (int i = 500; i < 1000; i++)
	//	cout << wt_mat[i][i] << endl;
	int resolution_40kb = 40000;
	int left = round(46500000/resolution_40kb), right = round(49500000/resolution_40kb);
	int del_left = round(47908463/resolution_40kb), del_right = round(47990293/resolution_40kb);
	ValueList result = predict_HiC(wt_mat_40kb, Firre_del_mat_40kb, left, del_left, del_right, right, MAX_LENGTH_OBJ, genomic_interaction_dist, "Debug/Firre_Del_40kb_full");
	cout << "40kb, PCC = " << result[0] << ", L1 = " << result[2] << endl;

	//int resolution_10kb = 10000;
	//int left = round(46500000/resolution_10kb), right = round(49500000/resolution_10kb);
	//int del_left = floor(47908463/resolution_10kb), del_right = floor(47990293/resolution_10kb) + 1;
	//ValueList result = predict_HiC(wt_mat_10kb, Firre_del_mat_10kb, left, del_left, del_right, right, MAX_LENGTH_OBJ, dummy_contact_prob_list, "Debug/Firre_Del_10kb");
	//cout << "10kb, PCC = " << result[0] << ", L1 = " << result[2] << endl;

	// For scoring
	result = predict_HiC(wt_mat_40kb, Firre_del_mat_40kb, del_left - fusion_score_window_size, del_left, del_right, del_right + fusion_score_window_size, 
		MAX_LENGTH_OBJ, genomic_interaction_dist, "Debug/Firre_Del_40kb_score_only");
	cout << "Firre deletion, fusion score = " << result[4] << endl;

}

void analyze_TAD_fusion_score (const string& primary_filename, const StringList& permutation_filename_list, const ValueList& top_cut_off_list, double score_threshold) {
	// For the primary
	SequenceMatrix primary_fusion_score_mat;
	read_sequence_list_file (primary_filename, TAD_FUSION_SCORE_FILE, 0, 1e9, primary_fusion_score_mat);

	StringList tool_name_list = {"Insulation_score", "Armatus          ", "CaTCH            ", "Arrowhead         "};
	
	vector<StringList> consensus_content (tool_name_list.size() + 1);
	for (int chr = 0; chr < 23; chr++) 
		for (int del = 0; del < primary_fusion_score_mat[chr].size(); del++) {
			int consensus = 0;
			for (int t = 0; t < tool_name_list.size(); t++)
				if (primary_fusion_score_mat[chr][del].score_list[t] >= score_threshold)
					consensus++;
			consensus_content[consensus].push_back(((chr <=21)? num_to_string(chr+1):"X") 
				+ "\t" + num_to_string(primary_fusion_score_mat[chr][del].left) 
				+ "\t" + num_to_string(primary_fusion_score_mat[chr][del].right) 
				+ "\t" + to_string(primary_fusion_score_mat[chr][del].score_list[6]));
		}

	ofstream consensus_file(primary_filename + "_consensus.dat");
	consensus_file << "#";
	for (int i = 0; i < consensus_content[0].size(); i++) {
		consensus_file << endl;
		for (int t = 0; t <= tool_name_list.size(); t++)
			if (i < consensus_content[t].size())
				consensus_file << "\t" << consensus_content[t][i];
			else
				consensus_file << "\t" << ".\t.\t.\t.";
	}
	consensus_file.close();
	
	//ValueList top_cut_off_list = {80, 90, 100, 110, 120};
	ValueList fusion_score_cut_off_list (top_cut_off_list.size(), 0);
	for (int i = 0; i < top_cut_off_list.size(); i++)
		fusion_score_cut_off_list[i] = find_top_TAD_fusion_score_threshold(primary_filename, top_cut_off_list[i]);

	ValueMatrix fusion_count_mat;
	ValueList fusion_score_sum_list;
	//vector<SequenceMatrix> permutation_fusion_score_cube;
	for (int file = 0; file < permutation_filename_list.size(); file++) {
		vector<SequenceMatrix> fusion_score_cube_tmp;
		read_TAD_fusion_score_permutation(permutation_filename_list[file], primary_fusion_score_mat, fusion_score_cube_tmp);
		for (int rep = 0; rep < fusion_score_cube_tmp.size(); rep++) {
			//cout << fusion_score_cube_tmp.size() << endl;
			ValueList fusion_count_list_tmp(top_cut_off_list.size(), 0);
			for (int chr = 0; chr < 23; chr++) {
				for (int del = 0; del < fusion_score_cube_tmp[rep][chr].size(); del++) {
					for (int k = 0; k < top_cut_off_list.size(); k++)
						if (fusion_score_cube_tmp[rep][chr][del].score_list[4] >= fusion_score_cut_off_list[k])
							fusion_count_list_tmp[k]++;
				}
			}
			fusion_count_mat.push_back(fusion_count_list_tmp);
			fusion_score_sum_list.push_back(find_fusion_score_sum(fusion_score_cube_tmp[rep], 4));
		}	
		//for (int rep = 0; rep < fusion_score_cube_tmp.size(); rep++)
		//	permutation_fusion_score_cube.push_back(fusion_score_cube_tmp[rep]);		
	}
	string summary_filename =  primary_filename + "_" + num_to_string(fusion_count_mat.size()) + "_permutation_summary.dat";
	ofstream out_file(summary_filename.c_str());
	out_file << "#";
	for (int rep = 0; rep <= fusion_count_mat.size(); rep++) {
		for (int i = 0; i < top_cut_off_list.size(); i++) {			
			if (rep == 0)
				out_file << "\t" << top_cut_off_list[i];
			else
				out_file << "\t" << fusion_count_mat[rep - 1][i];
		}
		if (rep == 0)
			out_file << "\t" << find_fusion_score_sum(primary_fusion_score_mat, 6);
		else
			out_file << "\t" << fusion_score_sum_list[rep - 1];
		out_file << endl;
	}
	out_file.close();
}

double find_top_TAD_fusion_score_threshold(const string& primary_filename, int top) {
	SequenceMatrix fusion_score_mat;
	read_sequence_list_file (primary_filename, TAD_FUSION_SCORE_FILE, 0, 1e9, fusion_score_mat);
	ValueList score_list_tmp;
	for (int chr = 0; chr < 23; chr++)
		for (int del = 0; del < fusion_score_mat[chr].size(); del++)
			score_list_tmp.push_back(fusion_score_mat[chr][del].score_list[6]);
	ValueList sorted_index_list, rank_tmp;
	sort_by_value (score_list_tmp, sorted_index_list, rank_tmp, false);
	return score_list_tmp[sorted_index_list[top]];
}

void make_TAD_fusion_score_binning (const SequenceMatrix& fusion_score_mat, double score_threshold, int bin_num, int bin_size, ValueMatrix& raw_score_mat, ValueList& del_num_list, ValueList& fusion_num_list, ValueList& percentage_list) {
	raw_score_mat.resize(bin_num);
	del_num_list.resize(bin_num, 0);
	fusion_num_list.resize(bin_num, 0);
	percentage_list.resize(bin_num, 0);
	for (int chr = 0; chr < 23; chr++) {
		for (int del = 0; del < fusion_score_mat[chr].size(); del++) {
			int del_length = fusion_score_mat[chr][del].right - fusion_score_mat[chr][del].left + 1;
			int bin = -1;
			for (int k = 0; k < bin_num; k++)
				if (del_length >= k*bin_size && del_length <= (k+1)*bin_size) {
					bin = k;
					break;
				}
			if (bin >= 0) {
				double fusion_score = fusion_score_mat[chr][del].score_list[6];
				raw_score_mat[bin].push_back(fusion_score);
				del_num_list[bin]++;
				if (fusion_score >= score_threshold)
					fusion_num_list[bin]++;
			}
		}
	}
	for (int bin = 0; bin < bin_num; bin++)
		percentage_list[bin] = ((del_num_list[bin] > 0)? (fusion_num_list[bin]/del_num_list[bin]):0);
}

void compare_TAD_fusion_score_by_binning(string filename_1, string filename_2, double score_threshold, int bin_num, int bin_size, string out_filename) {
	ValueMatrix raw_score_mat_1, raw_score_mat_2;
	ValueList del_num_list_1, del_num_list_2,
		fusion_num_list_1, fusion_num_list_2,
		percentage_list_1, percentage_list_2;
	SequenceMatrix fusion_score_mat_1, fusion_score_mat_2;
	read_sequence_list_file (filename_1, TAD_FUSION_SCORE_FILE, 0, 1e9, fusion_score_mat_1);
	read_sequence_list_file (filename_2, TAD_FUSION_SCORE_FILE, 0, 1e9, fusion_score_mat_2);
	make_TAD_fusion_score_binning (fusion_score_mat_1, score_threshold, bin_num, bin_size, raw_score_mat_1, del_num_list_1, fusion_num_list_1, percentage_list_1);
	make_TAD_fusion_score_binning (fusion_score_mat_2, score_threshold, bin_num, bin_size, raw_score_mat_2, del_num_list_2, fusion_num_list_2, percentage_list_2);

	
	int del_num = max(count(fusion_score_mat_1), count(fusion_score_mat_2));
	string raw_filename = out_filename + "_raw.dat";
	ofstream raw_file(raw_filename.c_str());
	raw_file << "#";
	for (int i = -1; i < del_num; i++) {
		if (i >= 0)
			raw_file << endl;
		for (int j = 0; j < bin_num; j++)
			if (i < 0)
				raw_file << j*bin_size << "\t";
			else {
				if (i < raw_score_mat_1[j].size())
					raw_file << raw_score_mat_1[j][i] << "\t";
				else
					raw_file << "." << "\t";
			}
		for (int j = 0; j < bin_num; j++)
			if (i < 0)
				raw_file << j*bin_size << "\t";
			else {
				if (i < raw_score_mat_2[j].size())
					raw_file << raw_score_mat_2[j][i] << "\t";
				else
					raw_file << "." << "\t";
			}

	}
	raw_file.close();

	string cut_off_filename = out_filename + "_cut_off.dat";
	ofstream cut_off_file(cut_off_filename.c_str());	
	for (int i = 0; i < bin_num; i++)
		cut_off_file << i*bin_size << "-" << (i+1)*bin_size 
			<< "\t" << del_num_list_1[i] << "\t" << fusion_num_list_1[i] << "\t" << percentage_list_1[i] 
			<< "\t" << del_num_list_2[i] << "\t" << fusion_num_list_2[i] << "\t" << percentage_list_2[i]
			<< endl;
	cut_off_file.close();
}

void print_prediction_case_study(int chr, int del_left, int del_right, int window_size, const Genomic_Interaction_Distribution& genomic_interaction_dist, bool is_with_K562) {
	string GM12878_intra_data_folder = "/share/hormozdiarilab/Data/HiC/GM12878_combined/5kb_resolution_intrachromosomal/chr";
	string K562_intra_data_folder = "/share/hormozdiarilab/Data/HiC/K562/K562_intrachromosomal/5kb_resolution_intrachromosomal/chr";
	string chr_name = ((chr < 23)? num_to_string(chr):"X");
	string GM_intra_filename = GM12878_intra_data_folder + chr_name +  "/MAPQGE30/chr" + chr_name + "_5kb.RAWobserved";
	ValueMatrix GM_hi_c_mat;
	read_hi_c_data(GM_intra_filename, RESOLUTION, GM_hi_c_mat);

	int del_left_bin = floor(del_left/RESOLUTION), 
		del_right_bin = floor(del_right/RESOLUTION);
	int left_1 = ((del_left_bin >= window_size)? (del_left_bin - window_size):0),
		right_2 = ((del_right_bin + window_size < GM_hi_c_mat.size())? (del_right_bin + window_size):(GM_hi_c_mat.size() - 1));
	printf("%d\t%d\t%d\t%d\n", left_1, del_left_bin - 1, del_right_bin + 1, right_2);
	
	string print_out_filename = "Debug/pred_vs_obs_chr_" + num_to_string(chr) + "_" + num_to_string(del_left) + "_" + num_to_string(del_right) + "_" + num_to_string(window_size) + ".dat";	
	if (is_with_K562) {
		string K_intra_filename = K562_intra_data_folder + chr_name +  "/MAPQGE30/chr" + chr_name + "_5kb.RAWobserved";
		ValueMatrix K_hi_c_mat;
		read_hi_c_data(K_intra_filename, RESOLUTION, K_hi_c_mat);
		ValueList new_result = predict_HiC(GM_hi_c_mat, K_hi_c_mat, left_1, del_left_bin - 1, del_right_bin + 1, right_2, 
			MAX_LENGTH_OBJ, genomic_interaction_dist, print_out_filename);
	}
	else
		ValueList new_result = predict_HiC(GM_hi_c_mat, GM_hi_c_mat, left_1, del_left_bin - 1, del_right_bin + 1, right_2, 
			MAX_LENGTH_OBJ, genomic_interaction_dist, print_out_filename);
}

void K562_prediction(string del_filename, const Genomic_Interaction_Distribution& genomic_interaction_dist, string out_filename) {
	// For deletions	
	SequenceMatrix del_mat;
	read_sequence_list_file (del_filename, DELETION_FILE, 5000, 10000000, del_mat);
	string GM12878_intra_data_folder = "/share/hormozdiarilab/Data/HiC/GM12878_combined/5kb_resolution_intrachromosomal/chr";
	string K562_intra_data_folder = "/share/hormozdiarilab/Data/HiC/K562/K562_intrachromosomal/5kb_resolution_intrachromosomal/chr";

	ValueList window_size_list = {50, 100, 200, 400};
	ValueMatrix chr_list(window_size_list.size()), 
		left_list(window_size_list.size()), 
		right_list(window_size_list.size()), 
		PCC_list(window_size_list.size()),
		L1_list(window_size_list.size()),
		max_list(window_size_list.size()),
		mean_list(window_size_list.size());

	for (int chr = 1; chr <= 23; chr++) {
		if (chr == 8)		// This deletion is from GM12878
			continue;
		string chr_name = ((chr < 23)? num_to_string(chr):"X");
		string GM_intra_filename = GM12878_intra_data_folder + chr_name +  "/MAPQGE30/chr" + chr_name + "_5kb.RAWobserved",
			K_intra_filename = K562_intra_data_folder + chr_name +  "/MAPQGE30/chr" + chr_name + "_5kb.RAWobserved"; 
		ValueMatrix GM_hi_c_mat, K_hi_c_mat;
		if (!del_mat[chr - 1].empty()) {
			read_hi_c_data(GM_intra_filename, RESOLUTION, GM_hi_c_mat);
			read_hi_c_data(K_intra_filename, RESOLUTION, K_hi_c_mat);
		}
		for (int del = 0; del < del_mat[chr - 1].size(); del++) {
			for (int w = 0; w < window_size_list.size(); w++) {
				int window_size = window_size_list[w];
				int del_left = floor(del_mat[chr - 1][del].left/RESOLUTION), 
					del_right = floor(del_mat[chr - 1][del].right/RESOLUTION);
				int left_1 = ((del_left >= window_size)? (del_left - window_size):0),
					right_2 = ((del_right + window_size < GM_hi_c_mat.size())? (del_right + window_size):(GM_hi_c_mat.size() - 1));
				//printf("%d\t%d\t%d\t%d\n", left_1, del_left - 1, del_right + 1, right_2);
	
				//ValueList GM_inter_list, K_inter_list;
				//for (int i = left_1; i <= del_left - 1; i++)
				//	for (int j = del_right + 1; j <= right_2; j++) {
				//		GM_inter_list.push_back(GM_hi_c_mat[i][j]);
				//		K_inter_list.push_back(K_hi_c_mat[i][j]);
				//	}
			
				//ValueList length_based_result = length_based_model_prediction(GM_hi_c_mat, K_hi_c_mat, left_1, del_left - 1, del_right + 1, right_2, ((window_size <= 100)? ("Debug/A1_chr_" + num_to_string(chr) + "_" + num_to_string(del_mat[chr-1][del].left) + "_" + num_to_string(del_mat[chr-1][del].right) + "_" + num_to_string(window_size) + ".dat"):""));
				//ValueList upgrade_length_based_result = predict_HiC(GM_hi_c_mat, K_hi_c_mat, left_1, del_left - 1, del_right + 1, right_2, MAX_LENGTH_OBJ, contact_prob_list, ((window_size <= 100)? ("Debug/A2_chr_" + num_to_string(chr) + "_" + num_to_string(del_mat[chr-1][del].left) + "_" + num_to_string(del_mat[chr-1][del].right) + "_" + num_to_string(window_size) + ".dat"):""));
;
				ValueList new_result = predict_HiC(GM_hi_c_mat, K_hi_c_mat, left_1, del_left - 1, del_right + 1, right_2, MAX_LENGTH_OBJ, genomic_interaction_dist, "DO_NOT_WRITE_A_FILE");
				chr_list[w].push_back(chr);
				left_list[w].push_back(del_mat[chr - 1][del].left);
				right_list[w].push_back(del_mat[chr - 1][del].right);
				PCC_list[w].push_back(new_result[0]);
				L1_list[w].push_back(new_result[2]);
				max_list[w].push_back(new_result[5]);
				mean_list[w].push_back(new_result[6]);
			}
		}	
	}
	ofstream out_file(out_filename.c_str());	
	for (int w = 0; w < window_size_list.size(); w++) {
		out_file << "Window size = " << window_size_list[w] << endl;
		out_file << "Chr\tDel_left\tDel_right\tDel_length\tPCC\tMax\tMean\tL1" << endl;
		for (int i = 0; i < chr_list[w].size(); i++)
			out_file << chr_list[w][i] << "\t" << left_list[w][i] << "\t" << right_list[w][i] 
				<< "\t" << right_list[w][i] - left_list[w][i] << "\t" << PCC_list[w][i]
				<< "\t" << max_list[w][i] << "\t" << mean_list[w][i]  << "\t" << L1_list[w][i] << endl; 
	}
	out_file.close();
}


void fig_S1_robustness_experiment() {
	ValueList max_interaction_length = {100, 150, 200, 250, 300};
	//ValueList max_training_window_length = {5e6, 5.5e6, 6e6, 6.5e6, 7e6};
	ValueList max_training_window_length = {4e6, 4.5e6, 5.5e6, 6e6};
	int remaining_length = 3e6;

	//for (int i = 0; i < max_interaction_length.size(); i++)
	//	export_all_hi_c_model (21, 21, max_interaction_length[i], 4e6, 1e6);
	//for (int i = 0; i < max_training_window_length.size(); i++)
	//	export_all_hi_c_model (21, 21, 100, max_training_window_length[i], max_training_window_length[i] - remaining_length, "Model_Robustness/");

	SequenceMatrix One_KG_del_mat;
	string One_KG_deletion_without_GM12878_filename = "/share/hormozdiarilab/Data/ValidatedVariants/1000G_SV_Phase3/1KG_without_GM12878_10000_10000000.dat";
	read_sequence_list_file (One_KG_deletion_without_GM12878_filename, DELETION_FILE, 5000, 1000000, One_KG_del_mat);

	int chr = 21;
	cout << "Robustness Exp: #Deletion = " << One_KG_del_mat[chr - 1].size() << endl;

	string GM12878_intra_data_folder = "/share/hormozdiarilab/Data/HiC/GM12878_combined/5kb_resolution_intrachromosomal/chr";
	string chr_name = ((chr < 23)? num_to_string(chr):"X");
	string GM_intra_filename = GM12878_intra_data_folder + chr_name + "/MAPQGE30/chr" + chr_name + "_5kb.RAWobserved";
	ValueMatrix hi_c_mat;
	read_hi_c_data(GM_intra_filename, RESOLUTION, hi_c_mat);


	ValueList score_window_length = {50, 100, 150, 200, 250, 300};
	ValueMatrix max_interaction_length_fusion_score_mat (max_interaction_length.size()),
		max_training_window_length_fusion_score_mat (max_training_window_length.size()),
		score_window_length_fusion_score_mat (score_window_length.size());
	

	string genomic_interaction_distribution_filename = "/share/hormozdiarilab/Data/HiC/Metric/Rao_loop_contact_prob_chr_1_23_max_length_1000000_rounded_500.dat";
	Genomic_Interaction_Distribution genomic_interaction_dist;
	genomic_interaction_dist.load(genomic_interaction_distribution_filename, 500);


	int window_size = 50;
	string model_folder = "Model_Robustness/";
	/*for (int i = 0; i < max_interaction_length.size(); i++) {
		HiCDataModelList hi_c_data_model_list;
		read_hi_c_data_model_list(model_folder + "chr" + chr_name + "_" + num_to_string(round(max_interaction_length[i])) + "_4000000_1000000.model", hi_c_data_model_list);
		for (int del = 0; del < One_KG_del_mat[chr - 1].size(); del++) {
			int left = floor(One_KG_del_mat[chr - 1][del].left/RESOLUTION), 
				right = ceil(One_KG_del_mat[chr - 1][del].right/RESOLUTION) - 1;
			if (right < hi_c_mat.size() - 5) {
				int left_boundary = ((left >= window_size)? (left - window_size):0),
					right_boundary = ((right + window_size < hi_c_mat.size())? (right + window_size):(hi_c_mat.size()-1));
				max_interaction_length_fusion_score_mat[i].push_back(predict_HiC(hi_c_mat, hi_c_data_model_list, 
					genomic_interaction_dist, left_boundary, left - 1, right + 1, right_boundary));	// corrected	
			}
		}
	}*/		
	for (int i = 0; i < max_training_window_length.size(); i++) {
		HiCDataModelList hi_c_data_model_list;
		StringList xxx = {"4000000_1000000", "4500000_1500000", "5500000_2500000", "6000000_3000000"};
		read_hi_c_data_model_list(model_folder + "chr" + chr_name + "_100_" + xxx[i] + ".model", hi_c_data_model_list);

		//read_hi_c_data_model_list(model_folder + "chr" + chr_name + "_100_" 
		//	+ num_to_string(round(max_training_window_length[i])) + "_" + num_to_string(round(max_training_window_length[i] - remaining_length))  + ".model", hi_c_data_model_list);
		for (int del = 0; del < One_KG_del_mat[chr - 1].size(); del++) {
			int left = floor(One_KG_del_mat[chr - 1][del].left/RESOLUTION), 
				right = ceil(One_KG_del_mat[chr - 1][del].right/RESOLUTION) - 1;
			if (right < hi_c_mat.size() - 5) {
				int left_boundary = ((left >= window_size)? (left - window_size):0),
					right_boundary = ((right + window_size < hi_c_mat.size())? (right + window_size):(hi_c_mat.size()-1));
				max_training_window_length_fusion_score_mat[i].push_back(predict_HiC(hi_c_mat, hi_c_data_model_list, 
					genomic_interaction_dist, left_boundary, left - 1, right + 1, right_boundary));	// corrected	
			}
		}
	}		
	/*for (int i = 0; i < score_window_length.size(); i++) {
		HiCDataModelList hi_c_data_model_list;
		read_hi_c_data_model_list(model_folder + "chr" + chr_name + "_100_7000000_3000000.model", hi_c_data_model_list);
		window_size = score_window_length[i];
		for (int del = 0; del < One_KG_del_mat[chr - 1].size(); del++) {
			int left = floor(One_KG_del_mat[chr - 1][del].left/RESOLUTION), 
				right = ceil(One_KG_del_mat[chr - 1][del].right/RESOLUTION) - 1;
			if (right < hi_c_mat.size() - 5) {
				int left_boundary = ((left >= window_size)? (left - window_size):0),
					right_boundary = ((right + window_size < hi_c_mat.size())? (right + window_size):(hi_c_mat.size()-1));
				score_window_length_fusion_score_mat[i].push_back(predict_HiC(hi_c_mat, hi_c_data_model_list, 
					genomic_interaction_dist, left_boundary, left - 1, right + 1, right_boundary));	// corrected	
			}
		}
	}*/		

	//cout << "Max interaction length" << endl;
	//for (int j = 0; j < max_interaction_length_fusion_score_mat[0].size(); j++) {
	//	for (int i = 0; i < max_interaction_length_fusion_score_mat.size(); i++)
	//		cout << "\t" << max_interaction_length_fusion_score_mat[i][j];
	//	cout << endl;
	//}
	//for (int i = 0; i < max_interaction_length.size(); i++) {
	//	cout << endl;
	//	for (int j = 0; j < max_interaction_length.size(); j++)
	//		cout_file << PCC(max_interaction_length_fusion_score_mat[i], max_interaction_length_fusion_score_mat[j], ZERO) << "\t";
	//}
	cout << "===============" << endl;
	cout << "Max training window length" << endl;
	for (int j = 0; j < max_training_window_length_fusion_score_mat[0].size(); j++) {
		for (int i = 0; i < max_training_window_length_fusion_score_mat.size(); i++)
			cout << "\t" << max_training_window_length_fusion_score_mat[i][j];
		cout << endl;
	}
	//for (int i = 0; i < max_training_window_length.size(); i++) {
	//	cout << endl;
	//	for (int j = 0; j < max_training_window_length.size(); j++)
	//		cout << PCC(max_training_window_length_fusion_score_mat[i], max_training_window_length_fusion_score_mat[j], ZERO) << "\t";
	//}	
	cout << "===============" << endl;
	//cout << "Score window size" << endl;
	//for (int j = 0; j < score_window_length_fusion_score_mat[0].size(); j++) {
	//	for (int i = 0; i < score_window_length_fusion_score_mat.size(); i++)
	//		cout << "\t" << score_window_length_fusion_score_mat[i][j];
	//	cout << endl;
	//}
	//for (int i = 0; i < score_window_length.size(); i++) {
	//	cout << endl;
	//	for (int j = 0; j < score_window_length.size(); j++)
	//		cout << PCC(score_window_length_fusion_score_mat[i], score_window_length_fusion_score_mat[j], ZERO) << "\t";
	//}
	//cout << "===============" << endl;

	ofstream out_file;
	/*out_file.open("Debug/robust_max_interaction_length.dat");
	out_file << "#";
	for (int i = 0; i < max_interaction_length.size(); i++) {
		for (int j = 0; j < max_interaction_length.size(); j++)
			out_file << endl << i + 1 << "\t" << j + 1 << "\t" 
				<< PCC(max_interaction_length_fusion_score_mat[i], max_interaction_length_fusion_score_mat[j], ZERO);
	}
	out_file.close();*/
	out_file.open("Debug/robust_training_window_size.dat");
	out_file << "#";
	for (int i = 0; i < max_training_window_length.size(); i++) {
		for (int j = 0; j < max_training_window_length.size(); j++)
			out_file << endl << i + 1 << "\t" << j + 1 << "\t"
				<< PCC(max_training_window_length_fusion_score_mat[i], max_training_window_length_fusion_score_mat[j], ZERO);
	}	
	out_file.close();
	/*out_file.open("Debug/robust_score_window_size.dat");
	out_file << "#";
	for (int i = 0; i < score_window_length.size(); i++) {
		for (int j = 0; j < score_window_length.size(); j++)
			out_file << endl << i + 1 << "\t" << j + 1 << "\t" 
				<< PCC(score_window_length_fusion_score_mat[i], score_window_length_fusion_score_mat[j], ZERO);
	}
	out_file.close();*/
}

void generate_TAD_fusion_score(const SequenceMatrix& del_mat, string hi_c_folder, string model_folder, string tad_annotation_folder, string CTCF_filename, const Genomic_Interaction_Distribution& genomic_interaction_dist, int window_size, int sample_num, vector<SequenceMatrix>& result_del_cube) {
	result_del_cube.clear();
	if (sample_num <= 0)
		result_del_cube.push_back(del_mat);
	else {
		double rand_resolution = 1e9;
		int rand_resolution_plus = round(rand_resolution + 1.0);
		int buffer_size = 100000;	// Minimum window region for estimating the TAD-fusion score
		for (int s = 0; s < sample_num; s++) {
			SequenceMatrix del_mat_tmp(23);
			ValueList chr_bin_length_list = {49847,
				48638,
				39585,
				38209,
				36178,
				34210,
				31826,
				29261,
				28223,
				27105,
				26990,
				26769,
				23022,
				21458,
				20505,
				18059,
				16240,
				15604,
				11824,
				12594,
				9623,
				10249,
				31052
			};
			for (int chr = 0; chr < 23; chr++) {
				del_mat_tmp[chr].resize(del_mat[chr].size());
				for (int del = 0; del < del_mat[chr].size(); del++) {
					int l = del_mat[chr][del].right - del_mat[chr][del].left;
					double r = (rand() % rand_resolution_plus)/rand_resolution;
					del_mat_tmp[chr][del].left = round(r*(chr_bin_length_list[chr]*RESOLUTION - l - 2*buffer_size)) + buffer_size;
					del_mat_tmp[chr][del].right = del_mat_tmp[chr][del].left + l;
				}			
			}
			result_del_cube.push_back(del_mat_tmp);
		}			
	}

	// Load other TAD callers
	string insulation_score_tad_filename = tad_annotation_folder + "insulation_score_TAD.dat",
		Armatus_tad_filename = tad_annotation_folder + "Amaratus_TAD.dat",
		CaTCH_tad_filename = tad_annotation_folder + "CaTCH_TAD.dat",
		Arrowhead_tad_filename = tad_annotation_folder + "GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt";

	SequenceMatrix insulation_score_tad_mat, Armatus_tad_mat, CaTCH_tad_mat, Arrowhead_tad_mat;
	int min_tad_length = 50000, max_tad_length = 1e7;
	read_sequence_list_file (insulation_score_tad_filename, INSULATION_SCORE_FILE, min_tad_length/10, max_tad_length, insulation_score_tad_mat);
	read_sequence_list_file (Armatus_tad_filename, JUICE_FILE, min_tad_length, max_tad_length, Armatus_tad_mat);
	read_sequence_list_file (CaTCH_tad_filename, JUICE_FILE, min_tad_length, max_tad_length, CaTCH_tad_mat);
	read_sequence_list_file (Arrowhead_tad_filename, JUICE_FILE, min_tad_length, max_tad_length, Arrowhead_tad_mat);

	// Load CFCF file
	SequenceMatrix CTCF_mat;
	read_sequence_list_file(CTCF_filename, CTCF_FILE, 0, 1e9, CTCF_mat);

	for (int chr = 1; chr <= 23; chr++) { 
		string chr_name = ((chr < 23)? num_to_string(chr):"X");
		string chr_raw_data_filename = hi_c_folder + chr_name +  "/MAPQGE30/chr" + chr_name + "_5kb.RAWobserved"; 
		ValueMatrix hi_c_mat;
		if (!del_mat[chr - 1].empty()) {
			read_hi_c_data(chr_raw_data_filename, RESOLUTION, hi_c_mat);
			HiCDataModelList hi_c_data_model_list;
			read_hi_c_data_model_list(model_folder + "chr" + chr_name + ".model", hi_c_data_model_list);
			cout << "Complete loading Hi-C model of chromosome " << chr_name << endl;

			for (int rep = 0; rep < result_del_cube.size(); rep++) {
				for (int del = 0; del < result_del_cube[rep][chr - 1].size(); del++) {
					result_del_cube[rep][chr - 1][del].score_list.resize(7);
					int left = floor(result_del_cube[rep][chr - 1][del].left/RESOLUTION), 
						right = ceil(result_del_cube[rep][chr - 1][del].right/RESOLUTION) - 1;
					if (right < hi_c_mat.size() - 5) {
						int left_boundary = ((left >= window_size)? (left - window_size):0),
						right_boundary = ((right + window_size < hi_c_mat.size())? (right + window_size):(hi_c_mat.size()-1));
						result_del_cube[rep][chr - 1][del].score_list[6] = predict_HiC(hi_c_mat, hi_c_data_model_list, genomic_interaction_dist, left_boundary, left - 1, right + 1, right_boundary);	// corrected
					}
					else
						result_del_cube[rep][chr - 1][del].score_list[6] = 0;
					//del_cube[rep][chr - 1][del].score_list[6] = predict_HiC(hi_c_mat, hi_c_data_model_list, contact_prob_list, left_boundary, left, right, right_boundary);
					//if (del_cube[rep][chr - 1][del].score_list[6] < 0)
					//	cout << "Negative score " << del_cube[rep][chr - 1][del].left << "\t" << del_cube[rep][chr - 1][del].right << endl;
	
				}
				conventional_estimate_fusion (result_del_cube[rep][chr-1], insulation_score_tad_mat[chr-1], 0);
				conventional_count_fusion (result_del_cube[rep][chr-1], Armatus_tad_mat[chr-1], 1);
				conventional_count_fusion (result_del_cube[rep][chr-1], CaTCH_tad_mat[chr-1], 2);
				conventional_count_fusion (result_del_cube[rep][chr-1], Arrowhead_tad_mat[chr-1], 3);
				for (int del = 0; del < result_del_cube[rep][chr - 1].size(); del++) {
					// Sum
					result_del_cube[rep][chr-1][del].score_list[4] = 0;
					for (int t = 0; t < 4; t++) 
						result_del_cube[rep][chr-1][del].score_list[4] += result_del_cube[rep][chr-1][del].score_list[t];
					result_del_cube[rep][chr-1][del].score_list[5] = 0;
					for (int ctcf = 0; ctcf < CTCF_mat[chr-1].size(); ctcf++)
						if (CTCF_mat[chr-1][ctcf].left >= result_del_cube[rep][chr-1][del].left 
							&& CTCF_mat[chr-1][ctcf].right <= result_del_cube[rep][chr-1][del].right)
								result_del_cube[rep][chr-1][del].score_list[5]++;
				}
		
			}		
		}
	}

}

void generate_TAD_fusion_score (string del_filename, int min_length, int max_length, string hi_c_folder, string model_folder, string tad_annotation_folder, string CTCF_filename, const Genomic_Interaction_Distribution& genomic_interaction_dist, int window_size, int sample_num, string print_prefix) {
	SequenceMatrix del_mat;
	read_sequence_list_file (del_filename, DELETION_FILE, min_length, max_length, del_mat);
	// Sampling deletions
	vector<SequenceMatrix> result_del_cube;
	generate_TAD_fusion_score (del_mat, hi_c_folder, model_folder, tad_annotation_folder, CTCF_filename, genomic_interaction_dist, window_size, sample_num, result_del_cube);
	if (sample_num <= 0) {
		string out_filename = print_prefix + "_" + num_to_string(window_size) + "_primary.dat";
		write_TAD_fusion_score_primary (result_del_cube[0], out_filename);
	}
	else {
		string out_filename = print_prefix + "_" + num_to_string(window_size) + "_permutation_" + num_to_string(sample_num) + ".dat";
		write_TAD_fusion_score_permutation(result_del_cube, out_filename);
	}
}

void write_TAD_fusion_score_primary (const SequenceMatrix& del_mat, string filename) {
	ofstream out_file(filename.c_str());
	out_file << "#Chr\tLeft\tRight\tInsulation-score\tArmatus\tCaTCH\tArrowhead\tAgreement\tCTCF\tFusion-score" << endl;
	for (int chr = 1; chr <= 23; chr++) {
		string chr_name = ((chr < 23)? num_to_string(chr):"X");
		for (int del = 0; del < del_mat[chr - 1].size(); del++) {
			Sequence tmp = del_mat[chr - 1][del];
			out_file << chr_name << "\t" << tmp.left << "\t" << tmp.right;
			for (int k = 0; k < tmp.score_list.size(); k++)
				out_file << "\t" << tmp.score_list[k];
			out_file << endl;
		}
	}
	out_file.close();
}

void write_TAD_fusion_score_permutation( const vector<SequenceMatrix>& del_cube, string filename) {
	ofstream out_file(filename.c_str());
	out_file << del_cube.size() << endl;
	for (int rep = 0; rep < del_cube.size(); rep++) {
		for (int chr = 0; chr < 23; chr++)
			for (int del = 0; del < del_cube[rep][chr].size(); del++) {
				out_file << del_cube[rep][chr][del].score_list[0];
				out_file << "\t" << del_cube[rep][chr][del].score_list[1];
				out_file << "\t" << del_cube[rep][chr][del].score_list[2];
				out_file << "\t" << del_cube[rep][chr][del].score_list[3];
				out_file << "\t" << del_cube[rep][chr][del].score_list[6];
				out_file << endl;
			}
	}
	out_file.close();	
}

void read_TAD_fusion_score_permutation (string filename, const SequenceMatrix& primary_fusion_score_mat, vector<SequenceMatrix>& fusion_score_cube) {
	ifstream in_file(filename.c_str());
	int sample_num;
	in_file >> sample_num;
	fusion_score_cube.resize(sample_num);
	for (int sample = 0; sample < sample_num; sample++) {
		fusion_score_cube[sample].resize(23);
		for (int chr = 0; chr < 23; chr++) {
			fusion_score_cube[sample][chr].resize(primary_fusion_score_mat[chr].size());
			for (int del = 0; del < primary_fusion_score_mat[chr].size(); del++) {
				fusion_score_cube[sample][chr][del].score_list.resize(5);
				for (int k = 0; k < 5; k++)
					in_file >> fusion_score_cube[sample][chr][del].score_list[k];
			}
		}
	}
	in_file.close();
}

double find_fusion_score_sum(const SequenceMatrix& del_mat, int score_col_index) {
	double sum = 0;
	for (int chr = 0; chr < 23; chr++)
		for (int del = 0; del < del_mat[chr].size(); del++)
			sum += del_mat[chr][del].score_list[score_col_index];
	return sum;
}


// Only interactions that have the length <= max_length_in_obj will contribute to the obj function
void fit_hi_c_data_model (const ValueMatrix& hi_c_mat, int left, int right, int max_length_in_obj, HiCDataModel& hi_c_data_model, bool print_solver_status) {
	//cout << hi_c_mat.size() << "\t" << left << "\t" << right << endl;
	CPXENVptr env = NULL;
	CPXLPptr lp = NULL;

	//int min_length_in_obj = 2;
	int length = right - left + 1;
	int alpha_size = length;
	int beta_size = 1;
	int e_size = ((length > max_length_in_obj)? ((length - max_length_in_obj)*max_length_in_obj + max_length_in_obj*(max_length_in_obj - 1)/2) : (length*(length - 1)/2));	// slack
	//e_size = e_size - (length - min_length_in_obj + 1)*(min_length_in_obj - 1) - (min_length_in_obj - 1)*(min_length_in_obj - 2)/2;
	int r_size = length;	// resistance

	int var_num = alpha_size + beta_size + e_size + r_size; 
	int status = 0;
	env = CPXopenCPLEX (&status);
	// Turn on output to the screen 
	status = CPXsetintparam (env, CPXPARAM_ScreenOutput, (print_solver_status? CPX_ON:CPX_OFF));
	// Turn on data checking 
	status = CPXsetintparam (env, CPXPARAM_Read_DataCheck, CPX_ON);
	//status = CPXsetintparam (env, CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
	// Create the problem 
	lp = CPXcreateprob (env, &status, "linear_model");
	// Populate the problem
	// Problem is minimization
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
	//lb[current_id] = -2;
	//ub[current_id] = -1;
	lb[current_id] = -2;
	ub[current_id] = -0.1;
	
	colname[current_id] = (char*) malloc(50*sizeof(char));
	sprintf(colname[current_id], "beta\0");
	obj[current_id] = 0;
	current_id++;
	// slack
	int tmp = current_id;
	int cell_num_tmp = 0;
	for (int i = 0; i < length; i++) {
		int max_i = ((i + max_length_in_obj < length)? (i + max_length_in_obj):(length-1));
		//for (int j = i + min_length_in_obj; j <= max_i; j++) {
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
	cout << "#Var " << current_id - tmp << "\t" << e_size << endl;
	// resistance
	for (int i = 0; i < length; i++) {
		lb[current_id] = 0;
		ub[current_id] = CPX_INFBOUND;
		colname[current_id] = (char*) malloc(50*sizeof(char));
		sprintf(colname[current_id], "r_%i\0", i);
		obj[current_id] = 0;
		current_id++;
	}
	// Obj and the bound
	status = CPXnewcols (env, lp, var_num, obj, lb, ub, NULL, colname);
	if (status) {
		cout << "LINH: Fail to set up (bound and obj) the LP problem!" << endl;
		cout << "(left,right) = " << left << "\t" << right << endl;
	}
	//else 
	//	cout << "LINH: Add all variables successfully!" << endl;
	// Add constraints
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
	cout << "Total: " << length << endl;
	for (int i = 0; i < length; i++) {
		int j_max = ((i + max_length_in_obj < length)? (i + max_length_in_obj):(length-1));
		//if (i > n - 10)
		//	cout << i << endl;
		//for (int j = i + min_length_in_obj; j <= j_max; j++) {
		for (int j = i + 1; j <= j_max; j++) {
			double logH = ((hi_c_mat[left + i][left + j] > 0)? log(hi_c_mat[left + i][left + j]) : log(ZERO));
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
			rmatval[current_cell_id + 2] = log(j - i);
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
			rmatval[current_cell_id + 2] = -log(j - i);
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
		//for (int i = 0; i < var_num; i++)
		//	cout << i << "\t" << lb[i] << "\t" << ub[i] << "\t" << obj[i] << "\t" << colname[i] << endl;
		//cout << "=============================" << endl;
		//for (int i = 0; i < current_row_id; i++)
		//	cout << i << "\t" << "rmatbeg: " << rmatbeg[i] << "\tsense: " << sense[i] << "\trhs: " << rhs[i] << "\t" << rowname[i] << endl;
		//for (int i = 0; i < current_cell_id; i++)
		//	cout << i << "\t" << "rmatind: " << rmatind[i] << "\t rmatval: " << rmatval[i] << endl;		
	}
	else {
		//cout << "LINH: Add all constraints successfully!" << endl;
	}	
	// Solve the problem
	time_t start, end;
	start = clock();

	status = CPXlpopt (env, lp); 
	
	end = clock();
	double total_time = (double)( end - start )/(double)CLOCKS_PER_SEC;
	if (print_solver_status)
		printf( "\t\t\tElapsed time : %0.3f \n", total_time );


	// Get the actual size of the problem
	int cur_numrows = CPXgetnumrows (env, lp);
	int cur_numcols = CPXgetnumcols (env, lp);
	double* x = (double*) malloc (cur_numcols*sizeof(double));
	double* slack = (double*) malloc (cur_numrows*sizeof(double));
	double* dj = (double*) malloc (cur_numcols*sizeof(double));
	double* pi = (double*) malloc (cur_numrows*sizeof(double));
	int solstat;
	double objval;
	// Get the solution
	//status = CPXwriteprob (env, lp, "linear_model.lp", NULL);
	status = CPXsolution (env, lp, &solstat, &objval, x, pi, slack, dj);
	if (status)
		cout << "LINH: Fail to solve the LP problem!" << endl;
	else {
		if (print_solver_status)
			cout << "LINH: Final obj value = " << objval << endl;
		//for (int i = 0; i < var_num; i++)
		//	cout << colname[i] << "\t" << x[i] << endl;
	}
	
	// Export the prediction
	hi_c_data_model.left = left;
	hi_c_data_model.beta = x[alpha_size];
	cout << "Beta_2 = " << x[alpha_size] << endl;
	resize_and_fill(hi_c_data_model.alpha, length, 0);
	resize_and_fill(hi_c_data_model.resistance, length, 0);
	for (int i = 0; i < length; i++) {
		hi_c_data_model.alpha[i] = x[i];
		hi_c_data_model.resistance[i] = x[alpha_size + 1 + e_size + i];
	}

	// Free the problem
	CPXfreeprob(env, &lp);
	// Free the environment
	CPXcloseCPLEX (&env);
	// Free all memory
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

ValueList predict_HiC (const ValueMatrix& hi_c_mat_wt, const ValueMatrix& hi_c_mat_mut, const HiCDataModel& hi_c_data_model, const Genomic_Interaction_Distribution& genomic_interaction_dist, int left_1, int right_1, int left_2, int right_2, string print_prefix) {	
	if (left_1 >= right_1 || left_2 >= right_2 || right_1 >= left_2) {
		cout << "ERROR: The left reigon & the right region must be large and not overlapped!" << endl;
		cout << "(left_1, right_1, left_2, right_2) = " << left_1 << "\t" << right_1 << "\t" << left_2 << "\t" << right_2 << endl;
	}

	int length = right_2 - left_1 + 1;
	// Export the prediction
	int index_shift = left_1 - hi_c_data_model.left;
	int length_1 = right_1 - left_1 + 1, length_2 = right_2 - left_2 + 1;
	ValueList pred, old, obs, norm_pred, norm_old;
	for (int i = 0; i < length_1; i++) {
		for (int j = 0; j < length_2; j++) {
			int pos_i = left_1 + i, pos_j = left_2 + j;
			int i_new = i, j_new = j + left_2 - left_1;		
			
			double alpha_i = hi_c_data_model.alpha[i_new + index_shift], alpha_j = hi_c_data_model.alpha[j_new + index_shift];
			double alpha_ij_2 = (alpha_i + alpha_j)/2;
			// Distance from bin i-th to bin j-th is j - i
			double log_pred = alpha_ij_2 + hi_c_data_model.beta*log(length_1 - i + j);	// new distance = length_1 - 1 - i + j + 1
			double log_old = alpha_ij_2 + hi_c_data_model.beta*log(length_1 - i + j + left_2 - right_1 - 1);
			//for (int k = i + 1; k < right_1 - left_1; k++)
			double del_resistance = 0;
			for (int k = i + 1; k <= right_1 - left_1; k++)
				del_resistance += hi_c_data_model.resistance[k + index_shift];
			//for (int k = 1; k < j; k++)
			for (int k = 0; k < j; k++)
				del_resistance += hi_c_data_model.resistance[left_2 - left_1 + k + index_shift];
			double total_resistance = del_resistance;
			for (int k = 1; k < left_2 - right_1; k++)
				total_resistance += hi_c_data_model.resistance[right_1 - left_1 + k + index_shift];

			pred.push_back(exp(log_pred - del_resistance));
			old.push_back(hi_c_mat_wt[pos_i][pos_j]);
			obs.push_back(hi_c_mat_mut[pos_i][pos_j]);
			if (alpha_i >= 0 && alpha_j >= 0) {
				norm_pred.push_back(exp(log_pred - del_resistance - alpha_ij_2));
				norm_old.push_back(exp(log_old - total_resistance - alpha_ij_2));
			}
			else {
				norm_pred.push_back(0);
				norm_old.push_back(0);
			}
			//cout << i << "\t" << j << "\t" << exp(alpha_i/2) << "\t" << exp(alpha_j/2) << endl;
		}		
	}
	double mean_exp = 0, max_exp = 0;
	if (!print_prefix.empty()) {
		// For plotting
		ofstream plot_file;
		bool is_write_file = (print_prefix.compare("DO_NOT_WRITE_A_FILE") != 0);
		if (is_write_file) {
			string plot_filename = print_prefix + "_plot_remodel.dat";
			plot_file.open(plot_filename.c_str());
		}
		int resolution = 20000;
		if (is_write_file)
			plot_file << "#" << left_1*resolution/1000000.0 << "\t" << right_2*resolution/1000000.0 << endl;
		double max_old = 0, sum_old = 0.0, max_pred = 0, sum_pred = 0.0, max_obs = 0, sum_obs = 0.0, n = 0;
		for (int i = left_1; i <= right_2; i++)
			for (int j = left_1; j <= right_2; j++) {				
				if (i != left_1 || j != left_1)
					if (is_write_file)
						plot_file << endl;
				if (is_write_file)
					plot_file << (i*resolution)/1000000.0 << "\t" << (j*resolution)/1000000.0 << "\t" << hi_c_mat_wt[i][j] << "\t" << hi_c_mat_mut[i][j];
				double pred_tmp;
				if ((i <= right_1 && j <= right_1) || (i >= left_2 && j >= left_2))
					pred_tmp = hi_c_mat_wt[i][j];
				else if ((i > right_1 && i < left_2) || (j > right_1 && j < left_2))
					pred_tmp = 0;
				else if (i <= right_1 && j >= left_2) {
					pred_tmp = pred[(i - left_1)*(right_2 - left_2 + 1) + j - left_2]; 
				
					n++;
					sum_old += hi_c_mat_wt[i][j];
					if (max_old < hi_c_mat_wt[i][j])
						max_old = hi_c_mat_wt[i][j];
					sum_pred += pred_tmp;
					if (max_pred < pred_tmp)
						max_pred = pred_tmp;
					sum_obs += hi_c_mat_mut[i][j];
					if (max_obs < hi_c_mat_mut[i][j])
						max_obs = hi_c_mat_mut[i][j];			
				}
				else if (i >= left_2 && j <= right_1) 
					pred_tmp = pred[(j - left_1)*(right_2 - left_2 + 1) + i - left_2];
				else 
					cout << "ERROR: Can not print for that case! " << i << "\t" << j 
						<< "\t" << left_1 << "\t" << right_1 << "\t" << left_2 << "\t" << right_2 << endl;
				if (is_write_file)
					plot_file << "\t" << pred_tmp << "\t" << pred_tmp - hi_c_mat_mut[i][j];
		
			}
		if (is_write_file) {
			plot_file << endl << "#\t\t" << max_old << "\t" << max_obs << "\t" << max_pred;
			plot_file << "#n = " << n << " vs " << norm_old.size() << " vs " << right_1 - left_1 + 1 << " x " << right_2 - left_2 + 1 << endl;
 			plot_file << endl << "#\t\t" << sum_old/n << "\t" << sum_obs/n << "\t" << sum_pred/n; 
			plot_file.close();
		}
		mean_exp = sum_obs/n;
		max_exp = max_obs; 
		// For debugging the score
		if (!genomic_interaction_dist.freq_distribution.empty() && is_write_file) {
			string score_filename = print_prefix + "_score.dat";
			ofstream score_file(score_filename.c_str());
			int index = 0;
			for (int i = 0; i < length_1; i++)
				for (int j = 0; j < length_2; j++) {
				double old_p = genomic_interaction_dist.get_prob(norm_old[index]),  
					new_p = genomic_interaction_dist.get_prob(norm_pred[index]);
					score_file << i + left_1 << "\t" << j + left_2 << "\t" << old[index] << "\t" << pred[index] 
						<< "\t" << norm_old[index] << "\t" << norm_pred[index]
						//<< "\t" << old_p << "\t" << ((new_p < old_p)? old_p:new_p) << endl;
						<< "\t" << old_p << "\t" << new_p << endl;
					index++;
				}
			score_file.close();
		}
	}

	double fusion_score = 0;
	if (!genomic_interaction_dist.freq_distribution.empty()) {
		for (int i = 0; i < norm_old.size(); i++) {
			double old_p = genomic_interaction_dist.get_prob(norm_old[i]),  
				new_p = genomic_interaction_dist.get_prob(norm_pred[i]);
			fusion_score += (new_p - old_p);	// corrected			
		}
	}

	ValueList tmp;
	tmp.push_back(PCC(pred, obs, ZERO));
	tmp.push_back(average_log_ratio(pred, obs, ZERO));
	tmp.push_back(L1(pred, obs));
	tmp.push_back(L2(pred, obs));
	tmp.push_back(fusion_score);
	tmp.push_back(max_exp);
	tmp.push_back(mean_exp);
	return tmp;
}

double predict_HiC (const ValueMatrix& hi_c_mat, const HiCDataModelList& hi_c_data_model_list, const Genomic_Interaction_Distribution& genomic_interaction_dist, int left_1, int right_1, int left_2, int right_2, string print_prefix) {
	if (left_1 >= right_1 || left_2 >= right_2) {
		cout << "ERROR: Can not predict because the left reigon or the right region is so small!" << endl;
		cout << left_1 << "\t" << right_1 << "\t" << left_2 << "\t" << right_2 << endl;
		return -1e6;
	}
	for (int k = 0; k < hi_c_data_model_list.size(); k++) {
		if (left_1 >= hi_c_data_model_list[k].left && left_2 >= right_1 && right_2 <= hi_c_data_model_list[k].get_right()) {
			ValueList result = predict_HiC (hi_c_mat, hi_c_mat, hi_c_data_model_list[k], genomic_interaction_dist, 
				left_1, right_1, left_2, right_2, print_prefix);	
			return result[4];
		}
	}	
	cout << "Error: Can not find the Hi-C data model! " << left_1 << "\t" << right_1 << "\t" << left_2 << "\t" << right_2 << endl;
	return -1e6;
}

ValueList predict_HiC (const ValueMatrix& hi_c_mat_wt, const ValueMatrix& hi_c_mat_mut, int left_1, int right_1, int left_2, int right_2, int max_length, const Genomic_Interaction_Distribution& genomic_interaction_dist, string print_prefix) {
	if (left_1 >= right_1 || left_2 >= right_2 || right_1 > left_2) {
		cout << "ERROR: The left reigon & the right region must be large and not overlapped!" << endl;
		cout << "(left_1, right_1, left_2, right_2) = " << left_1 << "\t" << right_1 << "\t" << left_2 << "\t" << right_2 << endl;
	}

	HiCDataModel hi_c_data_model;
	fit_hi_c_data_model(hi_c_mat_wt, left_1, right_2, max_length, hi_c_data_model, false);
	cout << "Finish training " << endl;

	return predict_HiC (hi_c_mat_wt, hi_c_mat_mut, hi_c_data_model, genomic_interaction_dist,
		left_1, right_1, left_2, right_2, print_prefix);
}

ValueList length_based_model_prediction (const ValueMatrix& hi_c_mat_wt, const ValueMatrix& hi_c_mat_mut, int left_1, int right_1, int left_2, int right_2, string print_filename) {
	int max_length = (right_1 - left_1 > right_2 - left_2)? (right_1 - left_1):(right_2 - left_2);
	vector<DistributionSummary> dist_sum_list(max_length + 1);
	for (int length = 1; length <= right_1 - left_1; length++)
		for (int pos = left_1; pos + length <= right_1; pos++)
			dist_sum_list[length].add_value(hi_c_mat_wt[pos][pos + length]);
	for (int length = 1; length <= right_2 - left_2; length++)
		for (int pos = left_2; pos + length <= right_2; pos++)
			dist_sum_list[length].add_value(hi_c_mat_wt[pos][pos + length]);
	
	ValueList x,y;
	for (int i = 1; i <= max_length; i++) {
		x.push_back(log(i));
		y.push_back(log((dist_sum_list[i].mean > 0)? dist_sum_list[i].mean:ZERO));
	}
	double a,b;
	line_LQ(x, y, a, b);

	//cout << "(a,b)" << a << "\t" << b << endl;

	ValueMatrix pred_mat(right_1 - left_1 + 1), obs_mat(right_1 - left_1 + 1);
	ValueList obs, pred;
	for (int i = left_1; i <= right_1; i++) {
		for (int j = left_2; j <= right_2; j++) {
			double obs_val = hi_c_mat_mut[i][j];
			double log_ij = log(right_1 - i + j - left_2 + 1);
			double pred_val = exp(a*log_ij + b);
			obs.push_back(obs_val);
			pred.push_back(pred_val);
			pred_mat[i - left_1].push_back(pred_val);
			obs_mat[i - left_1].push_back(obs_val);
		}
	}

	if (!print_filename.empty()) {
		ofstream out_file(print_filename.c_str());
		for (int i = 0; i < pred_mat.size(); i++)
			for (int j = 0; j < pred_mat[i].size(); j++)
				out_file << i << "\t" << j << "\t" << obs_mat[i][j] << "\t" << pred_mat[i][j] << endl;
		out_file.close();
	}

	//cout << "Naive list size " << pred.size() << "\t" << obs.size() << endl;
	ValueList tmp;
	tmp.push_back(PCC(pred, obs, ZERO));
	tmp.push_back(average_log_ratio(pred, obs, ZERO));
	tmp.push_back(L1(pred, obs));
	tmp.push_back(L2(pred, obs));
	return tmp;
}

void conventional_count_fusion (SequenceList& del_list, const SequenceList& tad_list, int score_index, bool is_printed) {
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
		del_list[k].score_list[score_index] = ((is_covered_left && is_covered_right)? 1:0);
		if (is_printed && is_covered_left && is_covered_right)
			cout << del_list[k].left << "\t" << del_list[k].right << endl;		
	}	
}

void conventional_estimate_fusion (SequenceList& del_list, const SequenceList& tad_boundary_list, int score_index, bool is_printed) {	
	for (int k = 0; k < del_list.size(); k++) {
		double fusion_score = 0;
		for (int i = 0; i < tad_boundary_list.size(); i++)
			if ((del_list[k].left >= tad_boundary_list[i].left && del_list[k].left <= tad_boundary_list[i].right && del_list[k].right >= tad_boundary_list[i].right) 
				|| (del_list[k].right >= tad_boundary_list[i].left && del_list[k].right <= tad_boundary_list[i].right && del_list[k].left <= tad_boundary_list[i].left) 
				|| (del_list[k].left <= tad_boundary_list[i].left && del_list[k].right >= tad_boundary_list[i].right)) { 
					if (fusion_score < tad_boundary_list[i].score_list[0])
						fusion_score = tad_boundary_list[i].score_list[0];
			}
		del_list[k].score_list[score_index] = fusion_score;
		if (is_printed && fusion_score > 0.5)
			cout << del_list[k].left << "\t" << del_list[k].right << endl;
	}	
}

void deletion_filter (SequenceMatrix& original_del_mat, const SequenceMatrix& filter_del_mat) {
	for (int chr = 0; chr < original_del_mat.size(); chr++) {
		//for (int l = 0; l < filter_del_mat[chr].size(); l++)
		//	cout << chr + 1 << "\t" << filter_del_mat[chr][l].left << "\t" << filter_del_mat[chr][l].right << endl; 
		SequenceList tmp;
		for (int del = 0; del < original_del_mat[chr].size(); del++) {
			bool is_repeated = false;
			for (int k = 0; k < filter_del_mat[chr].size(); k++) {
				//if (original_del_mat[chr][del].left == filter_del_mat[chr][k].left && original_del_mat[chr][del].right == filter_del_mat[chr][k].right) {
				//	is_repeated = true;
				//	break;
				//}
				int left_intersection = -1, right_intersection = -1;
				if (original_del_mat[chr][del].left >= filter_del_mat[chr][k].left && original_del_mat[chr][del].left <= filter_del_mat[chr][k].right)
					left_intersection = original_del_mat[chr][del].left;
				if (filter_del_mat[chr][k].left >= original_del_mat[chr][del].left && filter_del_mat[chr][k].left <= original_del_mat[chr][del].right)
					left_intersection = filter_del_mat[chr][k].left;
				if (original_del_mat[chr][del].right >= filter_del_mat[chr][k].left && original_del_mat[chr][del].right <= filter_del_mat[chr][k].right)
					right_intersection = original_del_mat[chr][del].right;
				if (filter_del_mat[chr][k].right >= original_del_mat[chr][del].left && filter_del_mat[chr][k].right <= original_del_mat[chr][del].right)
					right_intersection = filter_del_mat[chr][k].right;
				int intersection_length = right_intersection - left_intersection,
					original_del_length = original_del_mat[chr][del].right - original_del_mat[chr][del].left,
					filter_del_length = filter_del_mat[chr][k].right - filter_del_mat[chr][k].left;

				//if (original_del_mat[chr][del].left == filer_del_mat[chr][
				if (intersection_length + 1000 > original_del_length && intersection_length + 1000 > filter_del_length) {
					is_repeated = true;
					break;
				}
			}
			if (!is_repeated)
				tmp.push_back(original_del_mat[chr][del]);
		}
		//cout << chr << "\t" << original_del_mat[chr].size() << "\t" << tmp.size() << endl;
		original_del_mat[chr] = tmp;
	}
}

void deletion_filter (SequenceMatrix& original_del_mat, int max_del_length) {
	for (int chr = 0; chr < original_del_mat.size(); chr++) {
		SequenceList tmp;
		for (int del = 0; del < original_del_mat[chr].size(); del++) {
			if (original_del_mat[chr][del].right - original_del_mat[chr][del].left <= max_del_length)
				tmp.push_back(original_del_mat[chr][del]);
		}
		//cout << chr << "\t" << original_del_mat[chr].size() << "\t" << tmp.size() << endl;
		original_del_mat[chr] = tmp;
	}
}

void extract_1KG_without_GM12878 (string One_KG_deletion_filename, string GM12878_variant_filename, int min_length, int max_length) {	
	SequenceMatrix One_KG_del_mat, GM12878_del_mat;
	read_sequence_list_file (One_KG_deletion_filename, DELETION_1KG_FILE, 10000, 10000000, One_KG_del_mat);
	read_sequence_list_file (GM12878_variant_filename, DELETION_FILE, 10000, 10000000, GM12878_del_mat);
	cout << "All 1KG deletions " << count(One_KG_del_mat) << endl;
	cout << "GM12878 deletions " << count(GM12878_del_mat) << endl;
	deletion_filter(One_KG_del_mat, GM12878_del_mat);
	cout << "After removing GM12878 deletions " << count(One_KG_del_mat) << endl;
	string filename = "1KG_without_GM12878_" + num_to_string(min_length) + "_" + num_to_string(max_length) + ".dat";
	ofstream out_file(filename.c_str());
	for (int chr = 1; chr <= 23; chr++)
		for (int del = 0; del < One_KG_del_mat[chr - 1].size(); del++) {
			string chr_name = "chr" + ((chr < 23)? num_to_string(chr):"X");
			int length = One_KG_del_mat[chr - 1][del].right - One_KG_del_mat[chr - 1][del].left + 1;
			if (length >= min_length && length <= max_length)
	 			out_file << chr_name << "\t" << One_KG_del_mat[chr - 1][del].left << "\t" << One_KG_del_mat[chr - 1][del].right << endl;
		}
	out_file.close();
}

void print_distribution(int max_loop_length, int resolution) {
	SequenceMatrix loop_mat;
	string loop_filename = "/share/hormozdiarilab/Data/HiC/Interaction_Data/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt"; 
	//read_Rao_loop_file(loop_filename, max_loop_length, loop_mat);
	read_sequence_list_file (loop_filename, RAO_LOOP_FILE, 20*resolution, max_loop_length, loop_mat);

	string model_folder = "/share/hormozdiarilab/Data/HiC/Model/GM12878/ReModel/";
	string raw_data_file_folder = "/share/hormozdiarilab/Data/HiC/GM12878_combined/5kb_resolution_intrachromosomal/chr" ;	
	int max_loop_bin_length = round(max_loop_length/resolution);

	ExplicitDistribution genomic_contact_dist, all_dist;

	int chr_min = 1, chr_max = 23, round_factor = 500;
	for (int chr = chr_min; chr <= chr_max; chr++) {
		string chr_name = ((chr < 23)? num_to_string(chr):"X");
	
		HiCDataModelList hi_c_data_model_list;
		read_hi_c_data_model_list(model_folder + "chr" + chr_name + ".model", hi_c_data_model_list);

		string chr_raw_data_filename = raw_data_file_folder + chr_name +  "/MAPQGE30/chr" + chr_name + "_5kb.RAWobserved"; 

		ValueMatrix hi_c_mat, mark_mat;
		read_hi_c_data(chr_raw_data_filename, RESOLUTION, hi_c_mat);
		mark_mat.resize(hi_c_mat.size());
		for (int i = 0; i < hi_c_mat.size(); i++) 
			resize_and_fill(mark_mat[i], max_loop_bin_length + 1, 0);
	

		int loop_num = 0;
		ValueList covered_model_id_list(hi_c_mat.size(), -1);
		for (int loop = 0; loop < loop_mat[chr - 1].size(); loop++) {
			int left = round(loop_mat[chr - 1][loop].left/resolution);
			int right = round(loop_mat[chr - 1][loop].right/resolution);
			
			//cout << left << "\t" << right << endl;		

			int max_shift = -1, covered_model_id = -1;
			for (int k = 0; k < hi_c_data_model_list.size(); k++) {
				int l = left - hi_c_data_model_list[k].left, r = hi_c_data_model_list[k].get_right() - right;
				if (l >= 0 && r >= 0) {
					int min_shift_left_right = ((l >= r)? l:r);
					if (min_shift_left_right > max_shift) {
						covered_model_id = k;
						max_shift = min_shift_left_right;
					}
				}
			}
				
			if (covered_model_id >= 0) {
				loop_num++;
				for (int i = left; i <= right; i++) {
					if (i < hi_c_mat.size())
						covered_model_id_list[i] = covered_model_id;
					for (int j = i; j <= right; j++)
						if (mark_mat[i][j - i] < 1) 
							mark_mat[i][j - i] = 1;					
				}
			}		
		}	
		cout << "Chr " << chr << ", #loops = " << loop_mat[chr - 1].size() << ", #loops with length < " << max_loop_length << " & covered by the model: " << loop_num <<  endl;		
		for (int i = 0; i < hi_c_mat.size(); i++)
			for (int j = i; j < hi_c_mat.size(); j++)
				if (j - i <= max_loop_bin_length && covered_model_id_list[i] >= 0 && covered_model_id_list[j] >= 0) {
					double alpha_i = hi_c_data_model_list[covered_model_id_list[i]].get_alpha(i), 
						alpha_j = hi_c_data_model_list[covered_model_id_list[j]].get_alpha(j);
					double norm_freq = hi_c_mat[i][j]*exp(-(alpha_i + alpha_j)/2);
					
					//cout << norm_freq << "\t" << ((mark_mat[i][j - i] > 0)? 1:0) << endl;
					if (norm_freq < 0)
						cout << "ERROR: Norm contact freq can not be negative! " << norm_freq << endl;
					if (norm_freq < 20) {
						int rounded_freq = round(round_factor*norm_freq);
						all_dist.add_value(round(rounded_freq));
						if (mark_mat[i][j - i] > 0)
							genomic_contact_dist.add_value(rounded_freq);
					}
					else {
						if (mark_mat[i][j - i] < 1)
							cout << "WARNING: Norm contact freq of an non-interaction pair is so large! " << i*resolution << "\t" << j*resolution << "\t" << norm_freq;
					}
				}

	}
	string pos_fix = "_chr_" + num_to_string(chr_min) + "_" + num_to_string(chr_max) + "_max_length_" + num_to_string(max_loop_length) + "_rounded_" + num_to_string(round_factor) + ".dat";
	genomic_contact_dist.print_to_file("Debug/genomic_contact_distribution" + pos_fix);
	all_dist.print_to_file("Debug/all_distribution" + pos_fix);
	//top_restrict_all_dist.print_to_file("Debug/top_restrict_all_distribution_chr_" + num_to_string(chr) + ".dat");
		
	ValueList contact_prob_list;
	int pos = 0;	
	while (pos < all_dist.freq.size() && pos < genomic_contact_dist.freq.size() && all_dist.freq[pos] > 0) {
		double r = max(genomic_contact_dist.freq[pos]/all_dist.freq[pos], (contact_prob_list.empty()? 0:contact_prob_list.back()));
		contact_prob_list.push_back(r);
		pos++;
	};
	write_a_value_list ("Debug/Rao_loop_contact_prob" + pos_fix, contact_prob_list);
}

void export_all_hi_c_model (int chr_min, int chr_max, int max_length_obj, int window_length, int shift_length, string model_file_folder) { // Origin: window length = 4Mb, shift length = 1Mb, max_length_obj = 0.5Mb
	int hi_c_model_window_size = round(window_length/RESOLUTION);
	int shift_size = round(shift_length/RESOLUTION);
	// For generating Hi-C data model list
	string data_file_folder = "/share/hormozdiarilab/Data/HiC/GM12878_combined/5kb_resolution_intrachromosomal/chr";
	//string model_file_folder = "/share/hormozdiarilab/Data/HiC/Model/GM12878/";
	for (int chr = chr_min; chr <= chr_max; chr++) { 
		string chr_name = ((chr < 23)? num_to_string(chr):"X");
		string chr_raw_data_filename = data_file_folder + chr_name +  "/MAPQGE30/chr" + chr_name + "_5kb.RAWobserved"; 
		ValueMatrix hi_c_mat;
		read_hi_c_data(chr_raw_data_filename, RESOLUTION, hi_c_mat);
		export_hi_c_model_list (hi_c_mat, hi_c_model_window_size, shift_size, max_length_obj, 
			model_file_folder + "chr" + chr_name + "_" + to_string(max_length_obj) + "_" + to_string(window_length) + "_" + to_string(shift_length) + ".model");
		cout << "Complete generating Hi-C model of chromosome " << chr_name << endl;
	}
}

void export_hi_c_model_list (const ValueMatrix& hi_c_mat, int window_size, int shift_length, int max_length, string filename) {
	HiCDataModelList hi_c_data_model_list;
	int current_left = 0;
	while (current_left < hi_c_mat.size()) {
		HiCDataModel hi_c_data_model;
		int right_pos = (current_left + window_size < hi_c_mat.size())? (current_left + window_size): (hi_c_mat.size() - 1);
		cout << hi_c_mat.size() << "\t" << current_left << "\t" << right_pos << endl;
		fit_hi_c_data_model(hi_c_mat, current_left, right_pos, max_length, hi_c_data_model);
		hi_c_data_model.left = current_left;
		hi_c_data_model_list.push_back(hi_c_data_model);
		current_left += shift_length;	// Correct here for the final version
	}
	write_hi_c_data_model_list(filename, hi_c_data_model_list);
}

int HiCDataModel::get_right() const {
	return left + alpha.size() - 1;
}

double HiCDataModel::get_alpha(int bin_id) const {
	if (bin_id < left || bin_id - left >= alpha.size())
		cout << "ERROR: Try to access alpha from model with bin_id = " << bin_id << " while left = " << left << endl;
	return alpha.at(bin_id - left);
}

void Genomic_Interaction_Distribution::load(string filename, int round_factor_) {
	round_factor = round_factor_;
	read_a_value_list (filename, 1, freq_distribution);
}

double Genomic_Interaction_Distribution::get_prob(double contact_freq) const {
	int tmp = round(contact_freq*round_factor);
	if (tmp >= freq_distribution.size())
		return 1;
	else
		return freq_distribution[tmp];
}

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

void read_sequence_list_file (string filename, int file_format, int min_length, int max_length, SequenceMatrix& seq_mat) {
	IdList col_list;
	int skipping_header_line_num, chr_name_format;
	if (file_format == DELETION_FILE) {
		col_list = {0, 1, 2};
		skipping_header_line_num = 0;
		chr_name_format = 1;	// i.e. chr1
	}
	else if (file_format == DELETION_1KG_FILE) {
		col_list = {0, 1, 2};
		skipping_header_line_num = 1;
		chr_name_format = 1;	// i.e. chr1
	}
	else if (file_format == CTCF_FILE) {
		col_list = {0, 1, 2};
		skipping_header_line_num = 0;
		chr_name_format = 1;	// i.e. chr1
	}
	else if (file_format == RAO_LOOP_FILE) {
		col_list = {0, 1, 4};
		skipping_header_line_num = 1;
		chr_name_format = 0;	// i.e. 1
	}
	else if (file_format == JUICE_FILE) {
		col_list = {0, 1, 2};
		skipping_header_line_num = 1;
		chr_name_format = 0;	// i.e. 1
	}
	else if (file_format == INSULATION_SCORE_FILE) {
		col_list = {0, 1, 2, 3};
		skipping_header_line_num = 0;
		chr_name_format = 0;	// i.e. 1
	}
	else if (file_format == TAD_FUSION_SCORE_FILE) {
		col_list = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
		skipping_header_line_num = 0;
		chr_name_format = 0; 	// i.e. 1

	}
	else {
		cout << "Can not read the sequence list from file: Unkown file type " << file_format << endl;
		return;
	}

	seq_mat.resize(23);
	IdMap chr_name_to_id;
	for (int chr = 1; chr <= 23; chr++) {
		string chr_name;
		if (chr_name_format == 1)
			chr_name = "chr" + ((chr < 23)? num_to_string(chr):"X");
		else
			chr_name =  ((chr < 23)? num_to_string(chr):"X");
		chr_name_to_id[chr_name] = chr - 1;
		seq_mat[chr - 1].clear();
	}
	StringMatrix str_mat;
	read_a_table_file (filename, ROW_DELIM, FIELD_DELIM, '#', skipping_header_line_num, col_list, str_mat);
	int seq_num = 0;
	for (int i = 0; i < str_mat.size(); i++) {
		if (str_mat[i].size() >= 3 && chr_name_to_id.find(str_mat[i][0]) != chr_name_to_id.end()) {
			Sequence seq_tmp(atoi(str_mat[i][1].c_str()), atoi(str_mat[i][2].c_str()));
			int length = seq_tmp.right - seq_tmp.left;
			if (length >= min_length && length <= max_length) {
				if (col_list.size() > 3) {
					seq_tmp.score_list.resize(col_list.size() - 3);
					for (int k = 0; k < col_list.size() - 3; k++)
						seq_tmp.score_list[k] = stod(str_mat[i][k + 3].c_str());
				}
				seq_mat[chr_name_to_id[str_mat[i][0]]].push_back(seq_tmp);						
				seq_num++;
			}			
		}
		else
			cout << "\tERROR: Unknown chromosome of the sequence " << str_mat[i][0] << "\t" << str_mat[i][1] << "\t" << str_mat[i][2] << endl;
	}
	cout << "\t#Sequence: " << seq_num << endl;
}

int count(const SequenceMatrix& seq_mat) {
	int n = 0;
	for (int chr = 0; chr < 23; chr++)
		n += seq_mat[chr].size();
	return n;
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
void read_5C_file (string filename, int header_line_num, int matrix_size,  ValueMatrix& hi_c_matrix) {
	hi_c_matrix.resize(matrix_size);

	ifstream in_file(filename.c_str());
	if (in_file.is_open()) {		
		int line_num = 0;
		for (string row; getline(in_file, row, ROW_DELIM); ) {
			line_num++;
			istringstream ss(row);
			StringList record;
			for (string field; getline(ss, field, FIELD_DELIM); )
				record.push_back(field);
			if (line_num > header_line_num) {
				if (record.size() != matrix_size + 1) {
					cout << "ERROR: line " << line_num << " has " << record.size() << " elements while the matrix size = " << matrix_size << endl; 
					return;
				}
				hi_c_matrix[line_num - header_line_num - 1].clear();
				for (int i = 0; i < matrix_size; i++)
					hi_c_matrix[line_num - header_line_num - 1].push_back(atoi(record[i + 1].c_str()));			
			}
		}
		cout << "File " << filename << ", #lines " << line_num << endl;
		if (line_num != header_line_num + matrix_size)
			cout << "ERROR: Reading a 5C file with " << header_line_num << " header lines, matrix size = " << matrix_size << " with " << line_num << " lines!" << endl;
		in_file.close();
	}
	else {
		cout << "Can not open file " << filename << endl;
		return;
	}
}

Sequence::Sequence () {
	left = right = -1;
	score_list = {0};	
}
Sequence::Sequence (int l, int r, double s) {
	left = l;
	right = r;
	score_list.clear();
	score_list.push_back(s);
}	

DistributionSummary::DistributionSummary() {
	size = 0;
	mean = std = 0;
}
void DistributionSummary::add_value (double x) {
	double new_mean = (size*mean + x)/(size + 1);
	double tmp = (x*x + size*(std*std + mean*mean))/(size+1) - new_mean*new_mean;
	if (tmp > 1e-10)
		std = sqrt(tmp);
	else
		std = 0;
	mean = new_mean;
	size++;
}
void DistributionSummary::merge (const DistributionSummary& s) {
	int new_size = size + s.size;
	double new_mean = (size*mean + s.size*s.mean)/new_size;
	std = sqrt((size*(mean*mean + std*std) + s.size*(s.mean*s.mean + s.std*s.std))/new_size - new_mean*new_mean);
	mean = new_mean;
	size += s.size;
}


ExplicitDistribution::ExplicitDistribution() {
	freq.clear();
	n = 0;
}
void ExplicitDistribution::add_value (int x) {
	if (x < 0)
		cout << "Can not add a negative value " << x << " to a ExplicitDistribution!" << endl;
	else if (x < freq.size())
		freq[x]++;
	else {
		ValueList freq_tmp(x + 1, 0);
		for (int i = 0; i < freq.size(); i++)
			freq_tmp[i] = freq[i];
		freq_tmp[x] = 1;
		freq = freq_tmp;
	}
	n++;
}
double ExplicitDistribution::get_prob (int x) const {
	if (x >= 0 && x < freq.size())
		return freq[x]/n;
	else
		return 0;
}
void ExplicitDistribution::read_from_file (string filename) {
	IdList col_list = {0, 1};
	StringMatrix str_mat;
	read_a_table_file (filename, ROW_DELIM, FIELD_DELIM, '#', true, col_list, str_mat);
	freq.clear();
	n = 0;
	for (int i = 0; i < str_mat.size(); i++) {
		int x = atoi(str_mat[i][0].c_str());
		for ( ;x >= freq.size(); )
			freq.push_back(0);
		int f = atoi(str_mat[i][1].c_str());
		freq[x] = f;
		n += f;
	}
}
void ExplicitDistribution::print_to_file (string filename) const { 
	ofstream out_file (filename.c_str());
	for (int i = 0; i < freq.size(); i++)
		out_file << i << "\t" << freq[i] << endl;
	out_file.close();
}


double PCC (const ValueList& x, const ValueList& y, double zero) { // x: pred, y: obs
	if (x.size() != y.size() || x.empty() || y.empty()) {
		cout << "Error: PCC of two different vector sizes " << x.size() << " and " << y.size() << endl;
		return 0;
	}
	double x2 = 0, y2 = 0, xy = 0, x_mean, y_mean;
	for (int i = 0; i < x.size(); i++) {
		x_mean += ((x[i] > zero)? x[i]:zero);
		y_mean += ((y[i] > zero)? y[i]:zero);
	}
	x_mean /= x.size();
	y_mean /= y.size();
	for (int i = 0; i < x.size(); i++) {
		double w = ((x[i] > zero)? x[i]:zero);
		double z = ((y[i] > zero)? y[i]:zero);
		x2 += (w - x_mean)*(w - x_mean);
		y2 += (z - y_mean)*(z - y_mean);
		xy += (w - x_mean)*(z - y_mean);
	}
	cout << x_mean << "\t" << y_mean << "\t" << x2 << "\t" << y2 << "\t" << xy << endl;
	if (x2 < 1e-9 && y2 < 1e-9)
		return 1;
	else if (x2 < 1e-9 || y2 < 1e-9)
		return 0;
	else
		return xy/sqrt(x2*y2);		
}

double average_log_ratio (const ValueList& pred, const ValueList& obs, double zero) {
	double s = 0, count = 0;
	for (int i = 0; i < pred.size(); i++) {
		double w = ((pred[i] > zero)? pred[i]:zero);
		double z = ((obs[i] > zero)? obs[i]:zero);
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
							cout << "\tline " << line_num << " has only " << record.size() 
								<< " fields while the required index is " << col_list[i] << endl;
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

void resize_and_fill(ValueList& x, int n, int a) {
	if (x.size() != n)
		x.resize(n);
	for (int i = 0; i < n; i++)
		x[i] = a;
}

