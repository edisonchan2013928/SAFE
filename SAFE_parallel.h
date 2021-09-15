#pragma once
#ifndef SAFE_PARALLEL_H
#define SAFE_PARALLEL_H

#include "SAFE.h"

void initStat_parallel(int tid, statistics& stat, statistics*stat_thread_vec);
void SAFE_thread_parallel(int tid, statistics& stat, statistics*stat_thread_vec, Tree& t);
void SAFE_parallel(statistics& stat, Tree& t);

#endif