#include "SAFE_parallel.h"

void initStat_parallel(int tid, statistics& stat, statistics*stat_thread_vec)
{
	stat_thread_vec[tid].method = stat.method;
	stat_thread_vec[tid].n_row = stat.n_row;
	stat_thread_vec[tid].n_col = stat.n_col;
	stat_thread_vec[tid].kernel_type = stat.kernel_type;
	stat_thread_vec[tid].L = stat.L;
	stat_thread_vec[tid].epsilon = stat.epsilon;
	stat_thread_vec[tid].leafCapacity = stat.leafCapacity;
	stat_thread_vec[tid].n = stat.n;

	stat_thread_vec[tid].row_L = stat.row_L;
	stat_thread_vec[tid].row_U = stat.row_U;
	stat_thread_vec[tid].col_L = stat.col_L;
	stat_thread_vec[tid].col_U = stat.col_U;
	stat_thread_vec[tid].incr_row = stat.incr_row;
	stat_thread_vec[tid].incr_col = stat.incr_col;
	stat_thread_vec[tid].start_b = stat.start_b;
	stat_thread_vec[tid].incr_b = stat.incr_b;

	stat_thread_vec[tid].featureVector = stat.featureVector;
	stat_thread_vec[tid].queryVector = stat.queryVector;
	stat_thread_vec[tid].out_visual = stat.out_visual;
	
	for (int b = 0; b < (int)stat.bandwidth_vec.size(); b++)
		stat_thread_vec[tid].bandwidth_vec.push_back(stat.bandwidth_vec[b]);
}

void SAFE_thread_parallel(int tid, statistics& stat, statistics*stat_thread_vec, Tree& t)
{
	int total_num_entries = stat.n_row * stat.n_col;
	int num_entries_per_thread = (int)ceil((double)(stat.n_row * stat.n_col) / stat.thread_num);
	int start_entry = tid * num_entries_per_thread;
	int end_entry = (tid + 1)*num_entries_per_thread;
	//int v;

	initStat_parallel(tid, stat, stat_thread_vec);
	init_SAFE(stat_thread_vec[tid]);

	for (int e = start_entry; e < min(end_entry, total_num_entries); e++)
	{
		stat_thread_vec[tid].cur_r = (int)floor((double)e / stat.n_col);
		stat_thread_vec[tid].cur_c = e - stat_thread_vec[tid].cur_r*stat.n_col;
		stat_thread_vec[tid].cur_b_Num = stat.L - 1;

		clear(stat_thread_vec[tid]);
		t.RQS(stat_thread_vec[tid]);

		//if (stat_thread_vec[tid].range_result_idList.size() > 0)
		//	v = 0;

		share(stat_thread_vec[tid]);
		aggregate(stat_thread_vec[tid]);
	}
}

void SAFE_parallel(statistics& stat, Tree& t)
{
	statistics*stat_thread_vec;
	thread*th_vec = new thread[stat.thread_num];
	stat_thread_vec = new statistics[stat.thread_num];

	for (int tid = 0; tid < stat.thread_num; tid++)
		th_vec[tid] = thread(SAFE_thread_parallel, tid, ref(stat), stat_thread_vec, ref(t));

	for (int tid = 0; tid < stat.thread_num; tid++)
		th_vec[tid].join();
}