#include "SAFE_on_the_fly.h"

void init_SAFE_otf(statistics& stat)
{
	vector<Gap> temp_Gap_vec;
	vector<double> temp_double_vec;

	for (int r = 0; r < stat.n_row; r++)
	{
		for (int c = 0; c < stat.n_col; c++)
		{
			stat.Gap_vec_otf.push_back(temp_Gap_vec);
			stat.dist_vec_otf.push_back(temp_double_vec);
		}
	}
}

void aggregate_one_element(statistics& stat, int prev_element_num)
{
	Gap g;
	double sq_dist;
	int cur_element_num = prev_element_num + 1;
	int cur_q_pos= stat.cur_r*stat.n_col + stat.cur_c;

	if (prev_element_num == -1) //The first element
	{
		g.num_of_elements = 1;
		if (stat.kernel_type == 0) //Triangular kernel
			g.D_dist = stat.dist_vec[cur_element_num];
		if (stat.kernel_type == 1 || stat.kernel_type == 2) //Epanechnikov and Quartic kernels
		{
			sq_dist = stat.dist_vec[cur_element_num] * stat.dist_vec[cur_element_num];
			g.S_dist = sq_dist;

			if (stat.kernel_type == 2) //Quartic kernel
				g.Q_dist = sq_dist * sq_dist;
		}
	}
	else
	{
		g.num_of_elements = stat.Gap_vec_otf[cur_q_pos][prev_element_num].num_of_elements + 1;
		if (stat.kernel_type == 0) //Triangular kernel
			g.D_dist = stat.Gap_vec_otf[cur_q_pos][prev_element_num].D_dist + stat.dist_vec[cur_element_num];
		if (stat.kernel_type == 1 || stat.kernel_type == 2) //Epanechnikov and Quartic kernels
		{
			sq_dist = stat.dist_vec[cur_element_num] * stat.dist_vec[cur_element_num];
			g.S_dist = stat.Gap_vec_otf[cur_q_pos][prev_element_num].S_dist + sq_dist;

			if (stat.kernel_type == 2) //Quartic kernel
				g.Q_dist = stat.Gap_vec_otf[cur_q_pos][prev_element_num].Q_dist + sq_dist * sq_dist;
		}
	}

	stat.Gap_vec_otf[cur_q_pos].push_back(g);
}

void aggregate_exp(statistics& stat)
{
	vector<Gap> exp_gap_vec;
	int cur_q_pos = stat.cur_r*stat.n_col + stat.cur_c;
	int index = 0;
	int size = stat.Gap_vec_otf[cur_q_pos].size();

	if (size != 0)
		exp_gap_vec.push_back(stat.Gap_vec_otf[cur_q_pos][0]);
		
	while (index < size - 1)
	{
		index = min(2 * index + 1, size - 1);
		exp_gap_vec.push_back(stat.Gap_vec_otf[cur_q_pos][index]);
		stat.dist_vec_otf[cur_q_pos].push_back(stat.dist_vec[index]);
	}

	stat.Gap_vec_otf[cur_q_pos].clear();

	//Only insert those elements which have the exponential number in the index value.
	for (int i = 0; i < (int)exp_gap_vec.size(); i++)
		stat.Gap_vec_otf[cur_q_pos].push_back(exp_gap_vec[i]);

	stat.Gap_vec_otf[cur_q_pos].shrink_to_fit();
	
	exp_gap_vec.clear();
}

void aggregate_otf(statistics& stat)
{
	int num_of_elements = stat.dist_vec.size();
	int cur_q_pos = stat.cur_r*stat.n_col + stat.cur_c;

	sort(stat.dist_vec.begin(), stat.dist_vec.end());

	for (int i = 0; i < num_of_elements; i++)
	{
		aggregate_one_element(stat, i - 1);
		//SAFE_all (kd-tree and ball-tree)
		if (stat.method == 7 || stat.method == 8)
			stat.dist_vec_otf[cur_q_pos].push_back(stat.dist_vec[i]);
	}

	//SAFE_exp (kd-tree and ball-tree), removing the elements from the function aggregate_one_element
	if (stat.method == 9 || stat.method == 10)
		aggregate_exp(stat);

	stat.dist_vec.clear();
}

double KAF_b_otf(statistics& stat)
{
	int cur_q_pos = stat.cur_r*stat.n_col + stat.cur_c;
	double bandwidth = stat.bandwidth_vec[stat.cur_b_Num];
	double bandwidth_sq = bandwidth * bandwidth;
	vector<double>::iterator it;
	int pos; 
	double value;

	it = lower_bound(stat.dist_vec_otf[cur_q_pos].begin(), stat.dist_vec_otf[cur_q_pos].end(), bandwidth);
	pos = (it - stat.dist_vec_otf[cur_q_pos].begin()) - 1;

	if (pos == -1)
		return 0;

	if (stat.kernel_type == 0) //Triangular kernel
		value = stat.Gap_vec_otf[cur_q_pos][pos].num_of_elements - (stat.Gap_vec_otf[cur_q_pos][pos].D_dist / bandwidth);
	if (stat.kernel_type == 1) //Epanechnikov kernel
		value = stat.Gap_vec_otf[cur_q_pos][pos].num_of_elements - (stat.Gap_vec_otf[cur_q_pos][pos].S_dist / bandwidth_sq);
	if (stat.kernel_type == 2) //Quartic kernel
		value = stat.Gap_vec_otf[cur_q_pos][pos].num_of_elements - (2.0 / bandwidth_sq)*stat.Gap_vec_otf[cur_q_pos][pos].S_dist
		+ (stat.Gap_vec_otf[cur_q_pos][pos].Q_dist / (bandwidth_sq*bandwidth_sq));

	return value;
}

void SAFE_b_L_aggregate(statistics& stat, Tree& t)
{
	init_SAFE_otf(stat);

	for (int r = 0; r < stat.n_row; r++)
	{
		stat.cur_r = r;
		for (int c = 0; c < stat.n_col; c++)
		{
			stat.cur_c = c;
			stat.cur_b_Num = stat.L - 1;
			
			t.RQS(stat);

			#ifdef SPACE_MODE
				update_SAFE_all_space_iter(stat);
			#endif

			aggregate_otf(stat);
		}
	}
}

void SAFE_otf(statistics& stat, Tree& t)
{
	SAFE_b_L_aggregate(stat, t);

	for (int r = 0; r < stat.n_row; r++)
	{
		stat.cur_r = r;
		for (int c = 0; c < stat.n_col; c++)
		{
			stat.cur_c = c;

			for (int b = 0; b < stat.L; b++)
			{
				stat.cur_b_Num = b;
				stat.out_visual[b][r][c] = KAF_b_otf(stat);
			}
		}
	}

	#ifdef SPACE_MODE
		if (stat.method == 9 || stat.method == 10)
			obtain_SAFE_exp_space(stat);
	#endif
}

#ifdef SPACE_MODE
void init_space(statistics& stat)
{
	stat.SAFE_all_space = 0;
	stat.SAFE_exp_space = 0;
	stat.feature_vector_space = stat.n*stat.dim * sizeof(double);
	stat.tree_space = stat.n * 4 * sizeof(double);
	stat.visual_space = stat.n_row*stat.n_col * sizeof(double);
}

void update_SAFE_all_space_iter(statistics& stat)
{
	int num_of_elements = stat.dist_vec.size();
	stat.SAFE_all_space += num_of_elements * sizeof(double);
	if (stat.kernel_type == 0 || stat.kernel_type == 1)
		stat.SAFE_all_space += num_of_elements * (sizeof(double) + sizeof(int));
	if (stat.kernel_type == 2)
		stat.SAFE_all_space += num_of_elements * (2 * sizeof(double) + sizeof(int));
}

void obtain_SAFE_exp_space(statistics& stat)
{
	int index = 0;
	int num_of_elements;
	for (int r = 0; r < stat.n_row; r++)
	{
		for (int c = 0; c < stat.n_col; c++)
		{
			num_of_elements = stat.Gap_vec_otf[index].size();
			if (stat.kernel_type == 0 || stat.kernel_type == 1)
				stat.SAFE_exp_space += num_of_elements * (sizeof(double) + sizeof(int));
			if (stat.kernel_type == 2)
				stat.SAFE_exp_space += num_of_elements * (2 * sizeof(double) + sizeof(int));
		}
	}
}

void output_space(statistics& stat)
{
	cout << "Space consumption is: ";
	if (stat.method == 0)
		cout << (stat.feature_vector_space + stat.visual_space) << endl;

	if (stat.method >= 1 && stat.method <= 4)
		cout << (stat.feature_vector_space + stat.visual_space + stat.tree_space) << endl;

	if (stat.method >= 9 && stat.method <= 10)
	{
		cout << "feature vector space: "<< stat.feature_vector_space << endl;
		cout << "visual space: " << stat.visual_space << endl;
		cout << "Tree space: " << stat.tree_space << endl;
		cout << "SAFE_all: " << (stat.feature_vector_space + stat.visual_space + stat.tree_space + stat.SAFE_all_space) << endl;
		cout << "SAFE_exp: " << (stat.feature_vector_space + stat.visual_space + stat.tree_space + stat.SAFE_exp_space) << endl;
	}
}
#endif