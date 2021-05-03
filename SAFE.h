#pragma once
#ifndef SAFE_H
#define SAFE_H

#include "init_visual.h"
#include "kd_tree.h"
#include "ball_tree.h"

void init_SAFE(statistics& stat);
double KAF_b(statistics& stat);
void share(statistics& stat);
void aggregate(statistics& stat);
void clear(statistics& stat);
void SAFE(statistics& stat, Tree& t);

#endif