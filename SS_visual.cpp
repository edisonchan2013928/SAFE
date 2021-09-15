#include "SS_visual.h"

double SCAN(double*q, double bandwidth, statistics& stat)
{
	double temp_value;
	double sq_dist_value;
	double dist_value;
	double incr_value = 0;

	for (int i = 0; i < stat.n; i++)
	{
		sq_dist_value = sq_euclid_dist(q, stat.featureVector[i], stat.dim);

		if (sq_dist_value > bandwidth*bandwidth)
			continue;

		if (stat.kernel_type == 0) //Triangular kernel
		{
			dist_value = sqrt(sq_dist_value);
			incr_value += 1.0 - (1.0 / bandwidth)*dist_value;
		}
		if (stat.kernel_type == 1) //Epanechnikov kernel
			incr_value += 1.0 - (1.0 / (bandwidth*bandwidth))*sq_dist_value;
		if (stat.kernel_type == 2) //Quartic kernel
		{
			temp_value = 1.0 - (1.0 / (bandwidth*bandwidth))*sq_dist_value;
			incr_value += temp_value * temp_value;
		}
	}

	return incr_value;
}

void KDE_visual(statistics& stat)
{
	double bandwidth;
	double KDE_value;
	for (int b = 0; b < stat.L; b++)
	{
		bandwidth = stat.bandwidth_vec[b];
		for (int r = 0; r < stat.n_row; r++)
		{
			for (int c = 0; c < stat.n_col; c++)
			{
				KDE_value = SCAN(stat.queryVector[r*stat.n_col + c], bandwidth, stat);
				stat.out_visual[b][r][c] = KDE_value;
			}
		}
	}
}

double refinement(Node*curNode, statistics& stat)
{
	double f_cur;
	double temp_value;
	double sq_dist_value;
	double dist_value;
	int id;
	double*q = stat.queryVector[stat.cur_r*stat.n_col + stat.cur_c];
	double bandwidth = stat.bandwidth_vec[stat.cur_b_Num];
	f_cur = 0;

	for (int i = 0; i < (int)curNode->idList.size(); i++)
	{
		id = curNode->idList[i];
		sq_dist_value = sq_euclid_dist(q, stat.featureVector[id], stat.dim);

		if (sq_dist_value > bandwidth*bandwidth)
			continue;

		//SAFE_kd and SAFE_ball
		if (stat.method == 5 || stat.method == 6 || stat.method == 12)
			stat.range_result_idList.push_back(id);

		//SAFE_all (kd-tree and ball-tree)
		if (stat.method == 7 || stat.method == 8 || stat.method == 9 || stat.method == 10)
		{
			dist_value = sqrt(sq_dist_value);
			stat.dist_vec.push_back(dist_value);
		}

		if (stat.kernel_type == 0) //Triangular kernel
		{
			dist_value = sqrt(sq_dist_value);
			f_cur += 1.0 - (1.0 / bandwidth)*dist_value;
		}
		if (stat.kernel_type == 1) //Epanechnikov kernel
			f_cur += 1.0 - (1.0 / (bandwidth*bandwidth))*sq_dist_value;
		if (stat.kernel_type == 2) //Quartic kernel
		{
			temp_value = 1.0 - (1.0 / (bandwidth*bandwidth))*sq_dist_value;
			f_cur += temp_value * temp_value;
		}
	}

	return f_cur;
}