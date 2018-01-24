/* File created by Nathan Adelgren, Graduate Assisistant at Clemson University.
Started: 9/3/2014 
Finished: N/A
This work serves as part of my doctoral research under the advisement of Dr. Akshay Gupte. 
The goal of this program is to solve Biobjective Mixed-Integer Linear Programs using an MIP
based technique and employing the Parametric Simplex Algorithm (PSA) to find Pareto optimal
integer feasible line segments. (Small pieces of this code may be copied from a BB based 
code authored by Drs. Banu Soylu and Pietro Belotti. I collaborated on this work.)*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <time.h>
#include <pthread.h>

extern int exact_mips;
extern double time_limit;

struct split_pt //Points used for "objective branching"
{
	double f1;
	double f2;
	double *x;
	struct split_pt *prev;
  	struct split_pt *next;
};

typedef struct split_pt split_pt;

struct box
{
	double x_ub;
	double x_lb;
	double y_ub;
	double y_lb;
  	struct box *next;
};

struct slope_val
{
	double slope;
	double lb;
	double ub;
  	struct slope_val *next;
  	struct slope_val *prev;
};

int comparison_function4(const void *a, const void *b);

struct interval
{
	double lb;
	double ub;
	int dir;
	int miss_count;
	struct interval *next;
};

typedef struct box box;
typedef struct slope_val slope_val;
typedef struct interval interval;

struct store_it
{
    double ratio;
    int index;
};

//void add_binary_cut(CPXCENVptr env, CPXLPptr prob1, CPXLPptr prob2, double *x, int *indices, double *coefs);
