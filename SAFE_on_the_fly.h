#pragma once
#ifndef SAFE_ON_THE_FLY_H
#define SAFE_ON_THE_FLY_H

#include "init_visual.h"
#include "kd_tree.h"
#include "ball_tree.h"

void init_SAFE_otf(statistics& stat);
double KAF_b_otf(statistics& stat);
void SAFE_b_L_aggregate(statistics& stat, Tree& t);
void SAFE_otf(statistics& stat, Tree& t);

#ifdef SPACE_MODE
void init_space(statistics& stat);
void update_SAFE_all_space_iter(statistics& stat);
void obtain_SAFE_exp_space(statistics& stat);
void output_space(statistics& stat);
#endif

#endif