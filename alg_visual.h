#pragma once
#ifndef ALG_VISUAL_H
#define ALG_VISUAL_H

#include "init_visual.h"
#include "SS_visual.h"
#include "kd_tree.h"
#include "ball_tree.h"
#include "SAFE.h"
#include "SAFE_on_the_fly.h"
#include "SAFE_parallel.h"
#include "Validation.h"

struct pqNode
{
	Node*node;
	double discrepancy;
	double node_L;
	double node_U;
};

//This is the maximum heap
struct comparePriority
{
	bool operator()(pqNode& p1, pqNode& p2)
	{
		return p1.discrepancy < p2.discrepancy;
	}
};

typedef priority_queue<pqNode, vector<pqNode>, comparePriority> PQ;

void visual_Algorithm(statistics& stat);

#endif