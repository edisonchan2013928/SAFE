#include "alg_visual.h"

inline void clearHeap(PQ& pq)
{
	int heapSize = (int)pq.size();
	for (int h = 0; h < heapSize; h++)
		pq.pop();
}

void GBF_iter(Tree& tree, statistics& stat)
{
	static PQ pq;
	pqNode pq_entry;
	Node*curNode;
	double L, U;
	double f_cur;
	double val_R;
	
	Node*rootNode = tree.rootNode;
	double*cur_q = stat.queryVector[stat.cur_r*stat.n_col + stat.cur_c];
	stat.q_SquareNorm = computeSqNorm(cur_q, stat.dim);

	if (stat.method >= 1 && stat.method <= 2)
	{
		L = rootNode->LB(cur_q, stat);
		U = rootNode->UB(cur_q, stat);
	}

	pq_entry.node = rootNode;
	pq_entry.node_L = L;
	pq_entry.node_U = U;
	pq_entry.discrepancy = U - L;

	pq.push(pq_entry);

	while (pq.size() != 0)
	{
		if (validate_best(L, U, stat.epsilon, val_R) == true)
		{
			stat.out_visual[stat.cur_b_Num][stat.cur_r][stat.cur_c] = val_R;
			clearHeap(pq);
			return;
		}

		pq_entry = pq.top();
		pq.pop();

		L = L - pq_entry.node_L;
		U = U - pq_entry.node_U;

		curNode = pq_entry.node;

		//leaf Node
		if ((int)curNode->idList.size() <= tree.leafCapacity)
		{
			f_cur = refinement(curNode, stat);
			L = L + f_cur;
			U = U + f_cur;

			continue;
		}

		//Non-Leaf Node
		for (int c = 0; c < (int)curNode->childVector.size(); c++)
		{
			pq_entry.node_L = curNode->childVector[c]->LB(cur_q, stat);
			pq_entry.node_U = curNode->childVector[c]->UB(cur_q, stat);

			pq_entry.discrepancy = pq_entry.node_U - pq_entry.node_L;
			pq_entry.node = curNode->childVector[c];

			L = L + pq_entry.node_L;
			U = U + pq_entry.node_U;

			pq.push(pq_entry);
		}
	}

	stat.out_visual[stat.cur_b_Num][stat.cur_r][stat.cur_c] = L;
	clearHeap(pq);
}

void visual_Algorithm(statistics& stat)
{
	double run_time;
	//Different algorithms
	auto start_s = chrono::high_resolution_clock::now();

	//Preprocessing stage for indexing framework
	kdTree kd_Tree(stat.dim, stat.featureVector, stat.leafCapacity);
	ballTree ball_Tree(stat.dim, stat.featureVector, stat.leafCapacity);

	#ifdef SPACE_MODE
		init_space(stat);
	#endif

	if (stat.method == 1 || stat.method == 2 || stat.method == 3 || stat.method == 5 || stat.method == 7 || stat.method == 9 || stat.method == 11)
	{
		if (stat.method == 1 || stat.method == 3 || stat.method == 5 || stat.method == 7 || stat.method == 9 || stat.method == 11) //aKDE (LB_MBR and UB_MBR), RQS_kd, SAFE (kd-tree), SAFE_all (kd-tree), SAFE_exp (kd-tree)
			kd_Tree.rootNode = new kdNode();
		if (stat.method == 2) //QUAD (LB_QUAD and UB_QUAD)
			kd_Tree.rootNode = new kdQuadAugNode();

		if (stat.method == 3 || stat.method == 5 || stat.method == 11) //RQS_kd, SAFE_kd
			kd_Tree.init_RQS(stat);

		kd_Tree.build_kdTree(stat);
		kd_Tree.updateAugment((kdNode*)kd_Tree.rootNode);
	}
	if (stat.method == 4 || stat.method == 6 || stat.method == 8 || stat.method == 10) //RQS_ball, SAFE (ball-tree), SAFE_all (ball-tree), SAFE_exp (ball-tree)
	{
		ball_Tree.rootNode = new ballNode();
		ball_Tree.build_ballTree(stat);
	}
	
	//Online stage
	if (stat.method == 0) //SCAN method
		KDE_visual(stat);

	if ((stat.method >= 1 && stat.method <= 4) || stat.method == 11) //aKDE, QUAD, RQS_kd and RQS_ball methods
	{
		for (int b = 0; b < stat.L; b++)
		{
			stat.cur_b_Num = b;
			for (int r = 0; r < stat.n_row; r++)
			{
				stat.cur_r = r;
				for (int c = 0; c < stat.n_col; c++)
				{
					stat.cur_c = c;
					if (stat.method == 1 || stat.method == 2)
						GBF_iter(kd_Tree, stat);
					if (stat.method == 3 || stat.method == 11)
						kd_Tree.RQS(stat);
					if (stat.method == 4)
						ball_Tree.RQS(stat);
				}
			}
		}
	}

	if (stat.method == 5) //SAFE (kd-tree)
		SAFE(stat, kd_Tree);
	if (stat.method == 6) //SAFE (ball-tree)
		SAFE(stat, ball_Tree);
	if (stat.method == 7 || stat.method == 9) //SAFE_all and SAFE_exp (kd-tree)
		SAFE_otf(stat, kd_Tree);
	if (stat.method == 8 || stat.method == 10) //SAFE_all and SAFE_exp (ball-tree)
		SAFE_otf(stat, ball_Tree);
	if (stat.method == 11)
		update_visual(stat);

	auto end_s = chrono::high_resolution_clock::now();

	run_time = (chrono::duration_cast<chrono::nanoseconds>(end_s - start_s).count()) / 1000000000.0;
	std::cout << "method " << stat.method << ":" << run_time << endl;
	#ifdef SPACE_MODE
		output_space(stat);
	#endif
	output_visual(stat);
}

