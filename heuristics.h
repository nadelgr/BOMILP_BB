
#include<stdlib.h>
#include<stdio.h>
#include <math.h>
#include<time.h>
#include "cplex.h"
#include "max_tree.h"
#include "biobjective_bb.h"
#include "callbacks.h"
#include "bb-bicriteria.h"

int binary_heuristic(CPXCENVptr env, CPXLPptr lp1, CPXLPptr lp2, const int *indices);
int mixed_heuristic(CPXCENVptr env);
