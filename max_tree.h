#ifndef MAX_TREE_H
#define MAX_TREE_H

#include<stdlib.h>
#include<stdio.h>
#include <math.h>
#include<time.h>
#include "cplex.h"

/* File created by Nathan Adelgren, Graduate Assisistant at Clemson University on 6/23/2013 (updated 1/13/14) as part of doctoral research
under the advisement of Drs. Pietro Belotti and Akshay Gupte. The goal of this program is to take singleton points and 
line segments as input, treat these input as possible solutions (in the objective space) of a bi-objective minimization 
problem, and create a binary tree that stores the Pareto set.*/

/* This file uses the same code as the file for minimization. All inserted points are "reflected" about some point that is
guaranteed to be north-east of all inserted points. The ideal point is a good candidate for such a point. Once reflected, 
points can be treated as if the problem is a minimization. Then, when reporting the Pareto set, points are simply 
"reflected" back to their original position.*/

struct node
{
  	int type;			
	int subtree_size;	
	double slope;
	double nw_x;
	double nw_y;
	double se_x;
	double se_y;
	struct node *parent;
  	struct node *left;	
  	struct node *right;
};

typedef struct node node;

struct nadir
{
	double nw_x; /* end1_z1 < end2_z1   */
  	double se_x;
  	double nw_y;
  	double se_y;
  	int type;       /* 1-isolated, 2-segment  */
  	double rank_val;
};

struct closest_nodes
{
	struct node *closest;
  	struct node *next;
};

typedef struct closest_nodes closest_nodes;

struct dual_bd
{
  	int seqnum1;			
	int seqnum2;	
	double x1;
	double y1;
	double x2;
	double y2;
	struct dual_bd *next;
  	struct dual_bd *prev;
};

typedef struct dual_bd dual_bd;

struct user_data
{
	double *x_ws;
	double *x1;
	double *x2;
	int ws_still_feas;
	int ob1_still_feas;
	int ob2_still_feas;
	CPXLPptr prob;
	double f1;
	double f2;
};

typedef struct user_data user_data;

node *find_right_leaf(node *leaf);
node *find_left_leaf(node *leaf);
node *find_rightmost_leaf(node *leaf);
node *find_leftmost_leaf(node *leaf);
node *find_deepest_leaf_right(node *leaf);

node *find_first_node_left_of_val(double val, double lower_val, node *cur_node);
closest_nodes *find_two_nodes_left_of_val(double val, double lower_val, node *cur_node);
node *find_first_node_right_of_val(double val, double upper_val, node *cur_node);
closest_nodes *find_two_nodes_right_of_val(double val, double upper_val, node *cur_node);
void reduce_box(double *x, double *y, node *leaf);

void insert(int type, double nw_x, double nw_y, double se_x, double se_y, double slope, struct node **leaf, struct node **leaf2);
void insert_db(int type, double nw_x, double nw_y, double se_x, double se_y, double slope, struct node **leaf, struct node **leaf2);
void adjust_ranges(node *leaf);
void scan_proper_tree(node *leaf);
void print_preorder(node *tree, FILE *file);
int compare_ranks(const void* px, const void* py);
double calculate_hd_dist(struct node *primal, struct node *dual, double best_val);
double calculate_max_proximal_hd_dist(struct node *primal, struct node *dual, double best_val);
double get_nadirs(node *n1, int starting, double best_val);
double get_length(node *n1, double val);
int find_separations(node *n1, int starting);
node *copy_tree(struct node *current, struct node *parent);
double get_hypervolume(node *n1, int starting, double current_val);

/**************** Added for use with Chebychev branching strategy 10/15/14 ********************/

double *find_leaf_in_bounds(double x_lb,double x_ub,double y_lb,double y_ub,node *cur_node);
int semi_mock_insert(int type, double nw_x, double nw_y, double se_x, double se_y, double slope, struct node **leaf);
int mock_insert(int type, double nw_x, double nw_y, double se_x, double se_y, double slope, struct node **leaf);

extern node *prev_node;

#endif
