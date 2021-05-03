#include "SAFE.h"

void init_SAFE(statistics& stat)
{
	Gap temp_gap;

	for (int g = 0; g < stat.L; g++)
		stat.Gap_vec.push_back(temp_gap);
}

double KAF_b(statistics& stat)
{
	double value;
	double bandwidth = stat.bandwidth_vec[stat.cur_b_Num];
	double bandwidth_sq = bandwidth * bandwidth;
	if (stat.kernel_type == 0) //Triangular kernel
		value = stat.num_of_elements_cum - (stat.D_dist_cum / bandwidth);
	if (stat.kernel_type == 1) //Epanechnikov kernel
		value = stat.num_of_elements_cum - (stat.S_dist_cum / bandwidth_sq);
	if (stat.kernel_type == 2) //Quartic kernel
		value = stat.num_of_elements_cum - (2.0 / bandwidth_sq)*stat.S_dist_cum + (stat.Q_dist_cum / (bandwidth_sq*bandwidth_sq));

	return value;
}

void share(statistics& stat)
{
	int id;
	int pos;
	double dist_value;
	double dist_value_sq;
	double*q = stat.queryVector[stat.cur_r*stat.n_col + stat.cur_c];
	int result_size = (int)stat.range_result_idList.size();

	for (int r = 0; r < result_size; r++)
	{
		id = stat.range_result_idList[r];
		dist_value = euclid_dist(q, stat.featureVector[id], stat.dim); 
		pos = lower_bound(stat.bandwidth_vec.begin(), stat.bandwidth_vec.end(), dist_value) - stat.bandwidth_vec.begin();
		stat.Gap_vec[pos].num_of_elements++;

		//Triangular kernel
		if (stat.kernel_type == 0)
			stat.Gap_vec[pos].D_dist += dist_value;

		//Epanechnikov and Quartic kernels
		if (stat.kernel_type == 1 || stat.kernel_type == 2)
		{
			dist_value_sq = dist_value * dist_value;
			stat.Gap_vec[pos].S_dist += dist_value_sq;

			if (stat.kernel_type == 2)
				stat.Gap_vec[pos].Q_dist += dist_value_sq * dist_value_sq;
		}
	}
}

void aggregate(statistics& stat)
{
	if (stat.L > 1)
	{
		stat.cur_b_Num = 0;
		stat.num_of_elements_cum = stat.Gap_vec[0].num_of_elements;

		if (stat.kernel_type == 0) //Triangular kernel
			stat.D_dist_cum = stat.Gap_vec[0].D_dist;

		if (stat.kernel_type == 1 || stat.kernel_type == 2) //Epanechnikov and Quartic kernels
		{
			stat.S_dist_cum = stat.Gap_vec[0].S_dist;
			if (stat.kernel_type == 2)
				stat.Q_dist_cum = stat.Gap_vec[0].Q_dist;
		}			

		stat.out_visual[0][stat.cur_r][stat.cur_c] = KAF_b(stat);
	}

	for (int b = 1; b < stat.L-1; b++)
	{
		stat.num_of_elements_cum += stat.Gap_vec[b].num_of_elements;

		if (stat.kernel_type == 0) //Triangular kernel
			stat.D_dist_cum += stat.Gap_vec[b].D_dist;

		if (stat.kernel_type == 1 || stat.kernel_type == 2) //Epanechnikov and Quartic kernels
		{
			stat.S_dist_cum += stat.Gap_vec[b].S_dist;
			if (stat.kernel_type == 2)
				stat.Q_dist_cum += stat.Gap_vec[b].Q_dist;
		}

		stat.cur_b_Num = b;
		stat.out_visual[b][stat.cur_r][stat.cur_c] = KAF_b(stat);
	}
}

void clear(statistics& stat)
{
	stat.range_result_idList.clear();
	for (int g = 0; g < (int)stat.Gap_vec.size(); g++)
	{
		stat.Gap_vec[g].num_of_elements = 0;
		stat.Gap_vec[g].D_dist = 0;
		stat.Gap_vec[g].S_dist = 0;
		stat.Gap_vec[g].Q_dist = 0;
	}
}

void SAFE(statistics& stat, Tree& t)
{
	init_SAFE(stat);

	for (int r = 0; r < stat.n_row; r++)
	{
		stat.cur_r = r;
		for (int c = 0; c < stat.n_col; c++)
		{
			stat.cur_c = c;
			stat.cur_b_Num = stat.L - 1;
			clear(stat);
			t.RQS(stat);

			share(stat);
			aggregate(stat);
		}
	}
}