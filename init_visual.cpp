#include "init_visual.h"

void initStat(int argc, char**argv, statistics& stat)
{
	char*dataFileName = (char*)argv[1];
	char*bandwidthFileName = (char*)argv[2];
	stat.outMatrixFileName_Prefix = (char*)argv[3];
	stat.method = atoi(argv[4]);
	stat.n_row = atoi(argv[5]);
	stat.n_col = atoi(argv[6]);
	stat.kernel_type = atoi(argv[7]);
	stat.L = atoi(argv[8]);
	stat.epsilon = atof(argv[9]);
	stat.leafCapacity = atoi(argv[10]);

	if (stat.method == 11) //data_sampling
		stat.ori_n = atoi(argv[11]);

	/*char*dataFileName = (char*)"../../../Datasets/Chicago/Chicago";
	char*bandwidthFileName = (char*)"../../../Datasets/Chicago/Chicago";
	stat.outMatrixFileName_Prefix = (char*)"./Results/Chicago_M10_test";
	stat.method = 10;
	stat.n_row = 64;
	stat.n_col = 64;
	stat.kernel_type = 1;
	stat.L = 10;
	//stat.start_b = 1;
	//stat.incr_b = 2;
	stat.epsilon = 0.01;
	stat.leafCapacity = 40;*/

	extract_FeatureVector(dataFileName, stat);
	updateRegion(stat);
	initQuery(stat);

	load_bandwidth(bandwidthFileName, stat);
	//gen_bandwidth(stat);

	stat.out_visual = new double**[stat.L];
	for (int b = 0; b < stat.L; b++)
	{
		stat.out_visual[b] = new double*[stat.n_row];
		for (int r = 0; r < stat.n_row; r++)
			stat.out_visual[b][r] = new double[stat.n_col];
	}
}

void gen_bandwidth(statistics& stat)
{
	for (int b = 0; b < stat.L; b++)
		stat.bandwidth_vec.push_back(stat.start_b + stat.incr_b*b);
}

void load_bandwidth(char*bandwidthFileName, statistics& stat)
{
	fstream bandwidth_file;
	stringstream bandwidth_fileName_ss;
	double bandwidth_value;

	bandwidth_fileName_ss << bandwidthFileName << "_bandwidth_values_" << stat.L;

	bandwidth_file.open(bandwidth_fileName_ss.str().c_str());
	if (bandwidth_file.is_open() == false)
	{
		cout << "Cannot open bandwidth file!" << endl;
		exit(0);
	}

	for (int b = 0; b < stat.L; b++)
	{
		bandwidth_file >> bandwidth_value;
		stat.bandwidth_vec.push_back(bandwidth_value);
	}

	bandwidth_file.close();
}

void updateRegion(statistics& stat)
{
	stat.row_L = inf;
	stat.row_U = -inf;
	stat.col_L = inf;
	stat.col_U = -inf;

	for (int i = 0; i < stat.n; i++)
	{
		if (stat.featureVector[i][0] < stat.row_L)
			stat.row_L = stat.featureVector[i][0];
		if (stat.featureVector[i][0] > stat.row_U)
			stat.row_U = stat.featureVector[i][0];

		if (stat.featureVector[i][1] < stat.col_L)
			stat.col_L = stat.featureVector[i][1];
		if (stat.featureVector[i][1] > stat.col_U)
			stat.col_U = stat.featureVector[i][1];
	}
}

void initQuery(statistics& stat)
{
	int total_q = stat.n_row*stat.n_col;
	double x_coord;
	double y_coord;
	stat.queryVector = new double*[stat.n_row*stat.n_col];

	if (stat.n_row != 1 || stat.n_col != 1)
	{
		stat.incr_row = (stat.row_U - stat.row_L) / (stat.n_row - 1);
		stat.incr_col = (stat.col_U - stat.col_L) / (stat.n_col - 1);
	}

	if (stat.n_row == 1)
		stat.incr_row = 0;
	if (stat.n_col == 1)
		stat.incr_col = 0;

	for (int q = 0; q < total_q; q++)
		stat.queryVector[q] = new double[stat.dim];

	for (int r = 0; r < stat.n_row; r++)
	{
		x_coord = stat.row_L + r * stat.incr_row;
		for (int c = 0; c < stat.n_col; c++)
		{
			y_coord = stat.col_L + c * stat.incr_col;
			stat.queryVector[r*stat.n_col + c][0] = x_coord;
			stat.queryVector[r*stat.n_col + c][1] = y_coord;
		}
	}
}

void extract_FeatureVector(char*fileName, statistics& stat)
{
	//load data to feature array
	fstream file;
	int n;
	int dim;

	file.open(fileName);
	if (file.is_open() == false)
	{
		cout << "Cannot Open File!" << endl;
		exit(1);
	}

	file >> n;
	file >> dim;

	stat.n = n;
	if (stat.dim != dim)
	{
		cout << "dimension is not matched!" << endl;
		exit(0);
	}

	stat.featureVector = new double*[n];
	for (int i = 0; i < n; i++)
		stat.featureVector[i] = new double[dim];

	for (int i = 0; i < n; i++)
		for (int d = 0; d < dim; d++)
			file >> stat.featureVector[i][d];

	file.close();
}

void update_visual(statistics& stat)
{
	for (int b = 0; b < stat.L; b++)
		for (int r = 0; r < stat.n_row; r++)
			for (int c = 0; c < stat.n_col; c++)
				stat.out_visual[b][r][c] *= ((double)stat.ori_n / stat.n);
}

void output_visual(statistics& stat)
{
	stringstream outMatrixFileName_ss;
	//char*outMatrixFileName;
	fstream outMatrixFile;

	for (int b = 0; b < stat.L; b++)
	{
		outMatrixFileName_ss << stat.outMatrixFileName_Prefix << "_b" << stat.bandwidth_vec[b];
		//outMatrixFileName = (char*)outMatrixFileName_ss.str().c_str();
		outMatrixFile.open(outMatrixFileName_ss.str().c_str(), ios::in | ios::out | ios::trunc);
		if (outMatrixFile.is_open() == false)
		{
			cout << "Cannot open output file!" << endl;
			exit(0);
		}

		for (int r = 0; r < stat.n_row; r++)
		{
			for (int c = 0; c < stat.n_col; c++)
				outMatrixFile << stat.out_visual[b][r][c] << " ";
			outMatrixFile << endl;
		}

		outMatrixFile.close();
		outMatrixFile.clear();
		outMatrixFileName_ss.str("");
	}
}