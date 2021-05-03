#pragma once
#ifndef INIT_VISUAL_H
#define INIT_VISUAL_H

#define SPACE_MODE

#include "Library.h"

const double inf = 999999999999;
const double small_epsilon = 0.000000000000001;
const double pi = 3.14159265358979323846;

struct Gap
{
	double D_dist;
	double S_dist;
	double Q_dist;
	int num_of_elements;
};

struct statistics
{
	int n; //number of points in datasets
	double***out_visual;
	vector<double> bandwidth_vec; //a set of bandwidth values
	int n_row; //number of discrete region in the row, e.g. 100
	int n_col; //number of discrete region in the col, e.g. 100
	double row_L, row_U; //row region, e.g. [-2,2]
	double col_L, col_U; //col region e.g. [-2,2]
	double incr_row; //incremental row e.g. (2-(-2))/100=0.04
	double incr_col; //incremental col e.g. (2-(-2))/100=0.04
	int leafCapacity; //leaf capacity for the P set
	double**featureVector; //feature vector of all data points
	double**queryVector; //query vector
	char*outMatrixFileName_Prefix; //output file name
	int L; //number of bandwidth values

	//bandwidth generation
	double start_b; 
	double incr_b;

	//Used in indexing framework
	int cur_b_Num;
	int cur_r;
	int cur_c;
	double q_SquareNorm;
	double**query_boundary;

	//Used in Z-order method
	int ori_n;

	//Used in our method SAFE
	vector<int> range_result_idList;
	vector<Gap> Gap_vec;
	int num_of_elements_cum;
	double D_dist_cum;
	double S_dist_cum;
	double Q_dist_cum;

	//Used in our methods SAFE_all and SAFE_exp
	vector< vector<Gap> > Gap_vec_otf;
	vector< vector<double> > dist_vec_otf;
	vector<double> dist_vec;

	#ifdef SPACE_MODE
		double SAFE_all_space;
		double SAFE_exp_space;
		double tree_space;
		double feature_vector_space;
		double visual_space;
	#endif

	const int dim = 2; //dim = 2 (2d visualization)
	int method; //chosen method
	double epsilon; //epsilon value for epsilon-KVQ
	int kernel_type; //kernel type //0: Triangular kernel, 1: Epanechnikov kernel 2: Quartic kernel
};

void initStat(int argc, char**argv, statistics& stat);
void gen_bandwidth(statistics& stat);
void load_bandwidth(char*bandwidthFileName, statistics& stat);
void updateRegion(statistics& stat);
void initQuery(statistics& stat);
void extract_FeatureVector(char*fileName, statistics& stat);
void update_visual(statistics& stat);
void output_visual(statistics& stat);

#endif