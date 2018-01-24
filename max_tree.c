#include<stdlib.h>
#include<stdio.h>
#include <math.h>
#include<time.h>

#include "max_tree.h"
#include "user_set_parameters.h"

/* File created by Nathan Adelgren, Graduate Assisistant at Clemson University on 6/23/2013 (updated 1/13/14) as part of doctoral research
under the advisement of Drs. Pietro Belotti and Akshay Gupte. The goal of this program is to take singleton points and 
line segments as input, treat these input as possible solutions (in the objective space) of a bi-objective minimization 
problem, and create a binary tree that stores the Pareto set.*/

/* This file uses the same code as the file for minimization. All inserted points are "reflected" about some point that is
guaranteed to be north-east of all inserted points. The ideal point is a good candidate for such a point. Once reflected, 
points can be treated as if the problem is a minimization. Then, when reporting the Pareto set, points are simply 
"reflected" back to their original position.*/

#define BILLION  1000000000L;

FILE *file;		// Uncomment these lines if you wish to output MATLAB syntax
//static FILE *drawing;

double x_ideal;
double y_ideal;
int insert_to_potential_branch_tree = 0;

int temp = 0;
int res = 0;
double delta = 0.3;
double a = 0;
double b = 0;
double c = 0;
double d = 0;
int n = 0;
int m = 0;
//int p = 0;
int q = 0;
int printing = 0;

int counter = 0;
//int insert_counter = 0;

double new_nw_x;
double new_se_x;
double new_nw_y;
double new_se_y;
double x_intersect;
double y_intersect;

int temp_type;
double temp_nw_x;
double temp_se_x;
double temp_nw_y;
double temp_se_y;

clock_t start_time, finish_time;

struct node *root = NULL;

node *tree = NULL;
node *tree2 = NULL;
node *potential_branch_tree = NULL;

node *find_right_leaf(node *leaf);
node *find_left_leaf(node *leaf);
node *find_deepest_leaf_right(node *leaf);
void insert2(int type, double nw_x, double nw_y, double se_x, double se_y, double slope, struct node **leaf, struct node **leaf2);
void insert_db(int type, double nw_x, double nw_y, double se_x, double se_y, double slope, struct node **leaf, struct node **leaf2);
void adjust_ranges(node *leaf);
void scan_proper_tree(node *leaf);
node *copy_tree(struct node *current, struct node *parent);

int there_will_only_be_points = 0;
int integer_objective = 0;
double smallest_coef = 0.;
double integer_objective_lb = 0., extreme_x = 0., extreme_y = 0.;
/*int integer_bb = 1;*/

clock_t start_struct_timer;
clock_t finish_struct_timer;
double struct_time;

clock_t start_insert_timer;
clock_t finish_insert_timer;
double insert_time;

struct timespec start1, stop1;
double accum, struct_time2;
struct timespec start2, stop2;
double insert_time2;

/*int use_hausdorff = 0;*/
/*double epsilon = 0.0001;*/
double max_range = 0., x_range = 0., y_range = 0.;

int min(int x, int y)
{
  return ((x < y) ? x : y);
}

int max(int x, int y)
{
  return ((x > y) ? x : y);
}

void declare_ideals(double x_i, double y_i)
{
	x_ideal = x_i;
	y_ideal = y_i;
}

void open_files()
{
	file = fopen("pareto_output.txt", "w+");
}
void close_files()
{
	fclose(file);
}

void update_depth(node *leaf)	// This function can be used to update the depths and subtree sizes associated with each node within a certain path of the tree whenever a new node is added or a node has been removed.
{	
/*	printf("leaf: %lf,%lf to %lf,%lf\n",leaf->nw_x,leaf->nw_y,leaf->se_x,leaf->se_y);*/
/*	printf("tree: %lf,%lf to %lf,%lf\n",tree->nw_x,tree->nw_y,tree->se_x,tree->se_y);*/
	if(!leaf) return;
	if(leaf->right)
	{	
		if(leaf->left)
		{	
			leaf->subtree_size = leaf->left->subtree_size + leaf->right->subtree_size + 1;
		}
		else
		{
			leaf->subtree_size = leaf->right->subtree_size + 1;
		}
	}
	else if(leaf->left)
	{
		leaf->subtree_size = leaf->left->subtree_size + 1;
	}
	else
	{
		leaf->subtree_size = 1;
	}
	if(leaf && leaf != tree && leaf->parent)// && leaf->parent)
	{	
/*		if(printing)*/
/*		{*/
/*			printf("leaf size: %d\n",leaf->subtree_size);*/
/*			printf("tree size: %d\n",tree->subtree_size);*/
/*		}*/
		update_depth(leaf->parent);
	}
}

void print_preorder(node *n1, FILE *file)
{	
    	if (n1)
    	{	
   		//printf("%d \t %lf %lf %lf %lf, tree size: %d \n",tree->type, x_ideal-tree->nw_x, y_ideal-tree->nw_y, x_ideal-tree->se_x, y_ideal-tree->se_y, tree->subtree_size);
		// Uncomment the following lines in order to print MATLAB syntax for plotting the Pareto optimal solutions
		if(n1->type == 2)
		{
		  // the plotutil graph allows one to draw a scatter
		  // graph by specifying the x,y coordinates of points
		  // and segments. For segments:
		  //
		  //#m=1,S=3
		  //1 3
		  //2 4
		  //
		  // For points instead:
		  //
		  //#m=0,S=4
		  //2 5
		  //3 3
		  //
		  // this is now disabled as the bomip2graph.sh script does the same.


#ifdef DRAW_GRAPH
			fprintf(file,"#m=1,S=2\n%lf %lf\n%lf %lf\n\n",
				x_ideal-n1->nw_x,
				y_ideal-n1->nw_y,
				x_ideal-n1->se_x,
				y_ideal-n1->se_y);
#else
/*			fprintf(file,"%d: (%lf,%lf) to (%lf,%lf)\n",*/
/*				n1->type,*/
/*				x_ideal-n1->nw_x,*/
/*				y_ideal-n1->nw_y,*/
/*				x_ideal-n1->se_x,*/
/*				y_ideal-n1->se_y);*/
#endif
  	    		//fprintf(file,"plot([%.24lf,%.24lf],[%.24lf,%.24lf],'-ro');\n",x_ideal-n1->nw_x,x_ideal-n1->se_x,y_ideal-n1->nw_y,y_ideal-n1->se_y);
  	    		printf("plot([%lf,%lf],[%lf,%lf],'-ro');\n",x_ideal-n1->nw_x,x_ideal-n1->se_x,y_ideal-n1->nw_y,y_ideal-n1->se_y);
		}
		else 
		{
#ifdef DRAW_GRAPH
		  fprintf(file,"#m=0,S=4\n%lf %lf\n\n",
			  x_ideal-n1->nw_x,
			  y_ideal-n1->nw_y);
#else
/*		  fprintf(file,"%d: (%lf,%lf)\n",*/
/*			  n1->type,*/
/*			  x_ideal-n1->nw_x,*/
/*			  y_ideal-n1->nw_y);*/
#endif		
	//Also uncomment the following line in order to print the MATLAB syntax.
			//fprintf(file,"plot(%.24lf,%.24lf,'ro');\n",x_ideal-n1->nw_x,y_ideal-n1->nw_y);
			printf("plot(%lf,%lf,'ro');\n",x_ideal-n1->nw_x,y_ideal-n1->nw_y);
		}
		if(n1->left)
		{
			printf("left child:\n");
		  print_preorder(n1->left, file);
        	}
        	if(n1->right)
        	{
        		printf("right child:\n");
		  print_preorder(n1->right, file);
        	}
   	}
}

void print_inorder(node *n1, int color)
{	
    	if (n1)
    	{	
   		if(n1->left)
		{
			//printf("left child:\n");
		  print_inorder(n1->left,color);
        	}
		// Uncomment the following lines in order to print MATLAB syntax for plotting the Pareto optimal solutions
		if(n1->type == 2)
		{
  	    		if(color == 1) printf("plot([%lf,%lf],[%lf,%lf],'-ro');\n",x_ideal-n1->nw_x,x_ideal-n1->se_x,y_ideal-n1->nw_y,y_ideal-n1->se_y);
  	    		else if(color == 2) printf("plot([%lf,%lf],[%lf,%lf],'-bo');\n",x_ideal-n1->nw_x,x_ideal-n1->se_x,y_ideal-n1->nw_y,y_ideal-n1->se_y);
/*  	    		if(color == 1) printf("plot([%lf,%lf],[%lf,%lf],'-rx');\n",n1->nw_x,n1->se_x,n1->nw_y,n1->se_y);*/
/*  	    		else if(color == 2) printf("plot([%lf,%lf],[%lf,%lf],'-bo');\n",n1->nw_x,n1->se_x,n1->nw_y,n1->se_y);*/
		}
		else 
		{
			if(color == 1) printf("plot(%lf,%lf,'ro');\n",x_ideal-n1->nw_x,y_ideal-n1->nw_y);
			else if(color == 2) printf("plot(%lf,%lf,'bo');\n",x_ideal-n1->nw_x,y_ideal-n1->nw_y);
/*			if(color == 1) printf("plot(%lf,%lf,'rx');\n",n1->nw_x,n1->nw_y);*/
/*			else if(color == 2) printf("plot(%lf,%lf,'bo');\n",n1->nw_x,n1->nw_y);*/
		}
        	if(n1->right)
        	{
        		//printf("right child:\n");
		  print_inorder(n1->right,color);
        	}
   	}
}

void destroy_tree(struct node *leaf)
{	
  	if (leaf)
  	{
		if(leaf->left)
  		{
     			destroy_tree(leaf->left);
     		}
     		if(leaf->right)
     		{
      			destroy_tree(leaf->right);
      		}
		free( leaf ); 
		// leaf = NULL; // ineffectual
	}
}

// Uncomment the following function if you wish to output MATLAB syntax that plots the structure of the tree

/*void draw_tree(node *leaf, double pos1, double pos2, double inc) // Debugging procedure that prints MATLAB commands designed to draw the tree. The commands are stored in a text file: drawn_tree.txt ... when plotting, include the commands "figure();" and "hold on;"*/
/*{*/
/*	if(leaf->type == 1)*/
/*	{*/
/*		fprintf(drawing,"plot(%lf,%lf,'o');\n",pos1,pos2);*/
/*	}*/
/*	else*/
/*	{*/
/*		fprintf(drawing,"plot(%lf,%lf,'+');\n",pos1,pos2);*/
/*	}*/
/*	if(leaf->left)*/
/*	{*/
/*		fprintf(drawing,"plot([%lf,%lf],[%lf,%lf],'-');\n",pos1,pos1-inc,pos2,pos2-1);*/
/*		draw_tree(leaf->left,pos1-inc,pos2-1,inc*0.5);*/
/*	}*/
/*	if(leaf->right)*/
/*	{*/
/*		fprintf(drawing,"plot([%lf,%lf],[%lf,%lf],'-');\n",pos1,pos1+inc,pos2,pos2-1);*/
/*		draw_tree(leaf->right,pos1+inc,pos2-1,inc*0.5);*/
/*	}*/
/*} */

void scan_proper_tree(node *leaf) // Debugging procedure that makes sure only dummy nodes have children. Note that this should only be run before or after a user-implemented insertion into the tree (not an insertion generated by a recursion).
{
/*	printf("scanning node %lf,%lf to %lf,%lf\n",leaf->nw_x,leaf->nw_y,leaf->se_x,leaf->se_y);*/
	if(leaf->parent == leaf)
	{
		printf("error: leaf's parent is itself\n");
		exit(0);
	}
	if(leaf->left == leaf || leaf->right == leaf)
	{	
/*		printf("the leaf: %d: %lf,%lf to %lf,%lf\n",leaf->type,leaf->nw_x,leaf->nw_y,leaf->se_x,leaf->se_y);*/
		printf("error: leaf's child is itself\n");
		exit(0);
	}
	if(leaf->nw_x - leaf->se_x > .00001 )
	{
		printf("3 -----------\n");
/*		print_inorder(tree,1);*/
		printf("scanning node %lf,%lf to %lf,%lf\n",leaf->nw_x,leaf->nw_y,leaf->se_x,leaf->se_y);
		exit(0);
	}
	if(leaf->parent && (leaf->parent->left != leaf && leaf->parent->right != leaf))
	{
		printf("error: parent pointers are F-ed up\n");
		exit(0);
	}
	if(leaf->left)
	{
		if(leaf->left->parent != leaf)
		{
			printf("error: left child has incorrect parent pointer");
			exit(0);
		}	
		else if(leaf->left == leaf->parent)
		{	
/*			printf("1\n");*/
		}
		else
		{
			scan_proper_tree(leaf->left);
		}
	}
	else if(leaf->right)
	{
		if(leaf->right->parent != leaf)
		{
			printf("error: right child has incorrect parent pointer");
			exit(0);
		}
		else if(leaf->right == leaf->parent)
		{	
/*			printf("2\n");*/
		}	
		else
		{
			scan_proper_tree(leaf->right);
		}
	}
}

void Rebalance_Right1(struct node *leaf)
{	
/*	printf("rebalance right1\n");*/
	if(leaf->parent)
	{
		int t = 0;
		if(leaf == leaf->parent->right) leaf->parent->right = leaf->left;
		else {leaf->parent->left = leaf->left; t = 1;}
		if(leaf->left->right) leaf->left = leaf->left->right;
		else leaf->left = NULL;
		if(t == 0)
		{
			leaf->parent->right->right = leaf;
			leaf->parent->right->parent = leaf->parent;
			leaf->parent = leaf->parent->right;
		}
		else
		{
			leaf->parent->left->right = leaf;
			leaf->parent->left->parent = leaf->parent;
			leaf->parent = leaf->parent->left;
		}
		if(leaf->left) leaf->left->parent = leaf;
	}
	else
	{
		tree = leaf->left;
		leaf->parent = leaf->left;
		if(leaf->left->right) leaf->left = leaf->left->right;
		else leaf->left = NULL;
		leaf->parent->parent = NULL;
		leaf->parent->right = leaf;
		if(leaf->left) leaf->left->parent = leaf;
	}
	if(leaf->left) update_depth(leaf->left);
	else if(leaf->right) update_depth(leaf->right);
	else update_depth(leaf);
}

void Rebalance_Left1(struct node *leaf)
{	
/*	printf("rebalance left1\n");*/
	if(leaf->parent)
	{
		int t = 0;
		if(leaf == leaf->parent->left) leaf->parent->left = leaf->right;
		else {leaf->parent->right = leaf->right; t = 1;}
		if(leaf->right->left) leaf->right = leaf->right->left;
		else leaf->right = NULL;
		if(t == 0)
		{
			leaf->parent->left->left = leaf;
			leaf->parent->left->parent = leaf->parent;
			leaf->parent = leaf->parent->left;
		}
		else
		{
			leaf->parent->right->left = leaf;
			leaf->parent->right->parent = leaf->parent;
			leaf->parent = leaf->parent->right;
		}
		if(leaf->right) leaf->right->parent = leaf;
	}
	else
	{
		tree = leaf->right;
		leaf->parent = leaf->right;
		if(leaf->right->left) leaf->right = leaf->right->left;
		else leaf->right = NULL;
		leaf->parent->parent = NULL;
		leaf->parent->left = leaf;
		if(leaf->right) leaf->right->parent = leaf;
	}
	if(leaf->left) update_depth(leaf->left);
	else if(leaf->right) update_depth(leaf->right);
	else update_depth(leaf);
}

node *find_rightmost_leaf(struct node *leaf);
node *find_leftmost_leaf(struct node *leaf);
void delete_node(node *leaf);

void Rebalance_Left2(struct node *leaf, struct node *leaf2)
{	
/*	printf("rebalance left2\n");*/
	if(leaf->left)
	{
		node *temp_node = find_rightmost_leaf(leaf->left);
/*		printf("inserting 40\n");*/
		insert2(leaf->type,leaf->nw_x,leaf->nw_y,leaf->se_x,leaf->se_y,leaf->slope,&temp_node->right, &temp_node);
		if(temp_node->right) 
		{
			temp_node->right->parent = temp_node;
			update_depth(temp_node);
		}
	}
	else
	{
/*		printf("inserting 41\n");*/
		insert2(leaf->type,leaf->nw_x,leaf->nw_y,leaf->se_x,leaf->se_y,leaf->slope,&leaf->left, &leaf);
		if(leaf->left) 
		{
			leaf->left->parent = leaf;
			update_depth(leaf);
		}
	}
	leaf->type = leaf2->type;
	leaf->nw_x = leaf2->nw_x;
	leaf->nw_y = leaf2->nw_y;
	leaf->se_x = leaf2->se_x;
	leaf->se_y = leaf2->se_y;
	leaf->slope = leaf2->slope;
	delete_node(leaf2);
}

void Rebalance_Right2(struct node *leaf, struct node *leaf2)
{	
/*	printf("rebalance right2\n");*/
	if(leaf->right)
	{
		node *temp_node = find_leftmost_leaf(leaf->right);
/*		printf("inserting 42\n");*/
		if(printing) print_preorder(tree,NULL);
		insert2(leaf->type,leaf->nw_x,leaf->nw_y,leaf->se_x,leaf->se_y,leaf->slope,&temp_node->left, &temp_node);
		if(temp_node->left)
		{
			temp_node->left->parent = temp_node;
			update_depth(temp_node);
		}
	}
	else
	{
/*		printf("inserting 43\n");*/
		insert2(leaf->type,leaf->nw_x,leaf->nw_y,leaf->se_x,leaf->se_y,leaf->slope,&leaf->right, &leaf);
		if(leaf->right)
		{
			leaf->right->parent = leaf;
			update_depth(leaf);
		}
	}
	leaf->type = leaf2->type;
	leaf->nw_x = leaf2->nw_x;
	leaf->nw_y = leaf2->nw_y;
	leaf->se_x = leaf2->se_x;
	leaf->se_y = leaf2->se_y;
	leaf->slope = leaf2->slope;
	delete_node(leaf2);
}

node *find_rightmost_leaf(struct node *leaf)
{
	node *temp_node = leaf;
	while(temp_node->right)
	{
		temp_node = temp_node->right;
	}
	return temp_node;
}

node *find_leftmost_leaf(struct node *leaf)
{
	node *temp_node = leaf;
	while(temp_node->left)
	{
		temp_node = temp_node->left;
	}
	return temp_node;
}

void Rebalance(struct node *leaf);
//int counter1 = 0;
//int counter2 = 0;

void Rebalance2(struct node *leaf)
{	
	if(leaf->subtree_size > 2)
	{
		//if(leaf->left && leaf->left->subtree_size > (leaf->subtree_size - 1)/(2-delta))
		if(leaf->left && leaf->left->subtree_size > (leaf->subtree_size)/(2-delta))
		{	//print_preorder(leaf);
/*			if(leaf->right) printf("leaf subtree size: %d, left subtree size: %d, right subtree size: %d, cutoff value: %lf\n",leaf->subtree_size,leaf->left->subtree_size,leaf->right->subtree_size,(leaf->subtree_size - 1)/(2-delta));*/
			node *temp_node = find_rightmost_leaf(leaf->left);
			//if((leaf->left->left && leaf->left->left->subtree_size > (leaf->subtree_size - 1)/(2-delta) ) || temp_node == leaf->left)
			if((leaf->left->left && (leaf->left->left->subtree_size >= ((1-delta)*leaf->subtree_size)/(2-delta) - 1) && leaf->left->left->subtree_size <= (leaf->subtree_size)/(2-delta)) || temp_node == leaf->left)
			{	//counter1++;
				temp_node = leaf->left;
				Rebalance_Right1(leaf);
				Rebalance(temp_node);
			}
			else
			{	//counter2++;
				Rebalance_Right2(leaf,temp_node);
				Rebalance(leaf);
			}
		}
		//else if(leaf->right && leaf->right->subtree_size > (leaf->subtree_size - 1)/(2-delta))
		else if(leaf->right && leaf->right->subtree_size > (leaf->subtree_size)/(2-delta))
		{
			node *temp_node = find_leftmost_leaf(leaf->right);
			//if((leaf->right->right && leaf->right->right->subtree_size > (leaf->subtree_size - 1)/(2-delta) ) || temp_node == leaf->right)
			if((leaf->right->right && (leaf->right->right->subtree_size >= ((1-delta)*leaf->subtree_size)/(2-delta) - 1) && leaf->right->right->subtree_size <= (leaf->subtree_size)/(2-delta)) || temp_node == leaf->right)
			{	//counter1++;
				temp_node = leaf->right;
				Rebalance_Left1(leaf);
				Rebalance(temp_node);
			}
			else
			{	//counter2++;
				Rebalance_Left2(leaf,temp_node);
				Rebalance(leaf);
			}
		}
	}
}

void Rebalance(struct node *leaf)
{	
	Rebalance2(leaf);
	if(leaf->left && leaf->left->subtree_size > 2) Rebalance(leaf->left);
	if(leaf->right && leaf->right->subtree_size > 2) Rebalance(leaf->right);
}

int times_happened = 0;
int del_root = 0;

void delete_node(node *leaf)
{	
/*	printf("del_root: %d\n",del_root);*/
/*	printf("node to be deleted: %d: %lf,%lf to %lf,%lf\n",leaf->type,leaf->nw_x,leaf->nw_y,leaf->se_x,leaf->se_y);*/
/*	if(leaf)*/
/*	{*/

/*	if(!(leaf->left) && !(leaf->right))*/
/*	{*/
/*		printf("node to be deleted: %d: %lf,%lf to %lf,%lf\n",leaf->type,leaf->nw_x,leaf->nw_y,leaf->se_x,leaf->se_y);*/
/*		printf("no left or right. subtree size: %d\n",leaf->subtree_size);*/
/*		print_inorder(leaf,1);*/
/*	}*/

	if(leaf->left && (!leaf->right || (leaf->right && (leaf->left->subtree_size > leaf->right->subtree_size))))
	{	
/*		printf("replacing with a left node\n");*/
		node *temp_node2 = find_rightmost_leaf(leaf->left);
		leaf->type = temp_node2->type;
		leaf->nw_x = temp_node2->nw_x;
		leaf->nw_y = temp_node2->nw_y;
		leaf->se_x = temp_node2->se_x;
		leaf->se_y = temp_node2->se_y;
		leaf->slope = temp_node2->slope;
		delete_node(temp_node2);
	}
	else if(leaf->right)
	{	
/*		printf("replacing with a right node\n");*/
		node *temp_node2 = find_leftmost_leaf(leaf->right);
		leaf->type = temp_node2->type;
		leaf->nw_x = temp_node2->nw_x;
		leaf->nw_y = temp_node2->nw_y;
		leaf->se_x = temp_node2->se_x;
		leaf->se_y = temp_node2->se_y;
		leaf->slope = temp_node2->slope;
		delete_node(temp_node2);
	}
	else
	{	
/*		printf("else\n");*/
		node *temp_node2 = leaf->parent;
		int t = 0;
		if(!temp_node2) t = -1;
		if(leaf == tree || leaf == tree2) t = 1;
		if(leaf->parent && leaf->parent->left && leaf == leaf->parent->left) t = 2;
		else if(leaf->parent && leaf == leaf->parent->right) t = 3;
/*		printf("t:\t%d\n",t);*/
		if(leaf && leaf->nw_x < -1335.17 && leaf->nw_x > -1335.179) 
		{
			times_happened++;
/*			printf("th: %d\n",times_happened);*/
/*			if(times_happened == 4) printing = 1;*/
		}
/*		print_preorder(tree,NULL);*/
		free(leaf);
		leaf = NULL;
		if(t == 2) temp_node2->left = NULL;
		else if(t == 3) temp_node2->right = NULL;
		if(temp_node2) update_depth(temp_node2);
	}
/*	}*/
}

/*void do_nothing(int i)*/
/*{*/
/*}*/
/*node_to_delete *first1 = NULL;*/
/*node_to_delete *last1 = NULL;*/

/*void add_node_to_delete(node *nd)*/
/*{*/
/*	if(first1 == NULL)*/
/*	{*/
/*		first1 = (struct node_to_delete*) malloc( sizeof( struct node_to_delete ) );*/
/*		first1->current = nd;*/
/*		first1->next = NULL;*/
/*	}*/
/*}*/
//node *freed_node = NULL;
node *prev_node = NULL;

void clean_it(node *nd1, int dir)
{
	node *temp = NULL;
	if(dir != 0) temp = nd1->parent;
	if(nd1 && nd1->left)
	{
		clean_it(nd1->left,1);
	}
	if(nd1 && prev_node && nd1 != prev_node)
	{
		double xdiff = 1000.;
		double ydiff = 1000.;
		double slopediff = fabs(nd1->slope - prev_node->slope);
/*		printf("node1: %d (%lf,%lf),(%lf,%lf)\tprev_node: %d (%lf,%lf),(%lf,%lf)\n",nd1->type,nd1->nw_x,nd1->nw_y,nd1->se_x,nd1->se_y,prev_node->type,prev_node->nw_x,prev_node->nw_y,prev_node->se_x,prev_node->se_y);*/
		xdiff = nd1->nw_x - prev_node->se_x;
		ydiff = prev_node->se_y - nd1->nw_y;
		if(prev_node->type == 1)
		{
/*			if( xdiff < .00001 && ydiff < .00001 )*/
/*			{*/
/*				//freed_node = prev_node;*/
/*				delete_node(prev_node);*/
/*			}*/
			//else 
			if( xdiff < .01 && ydiff < .01 )
			{
				if( xdiff < .005 || ydiff < .005 )
				{
					nd1->type = 2;
					nd1->nw_x = prev_node->nw_x;
					nd1->nw_y = prev_node->nw_y;
					nd1->slope = (nd1->se_y-prev_node->nw_y)/(nd1->se_x-prev_node->nw_x);
					//freed_node = prev_node;
					delete_node(prev_node);
				}
			}
		}
		else if(prev_node->type == 2)
		{
			if( slopediff < .01 && xdiff < .01 && ydiff < .01 )
			{
				if( xdiff < .005 || ydiff < .005 )
				{
					nd1->type = 2;
					nd1->nw_x = prev_node->nw_x;
					nd1->nw_y = prev_node->nw_y;
					nd1->slope = (nd1->se_y-prev_node->nw_y)/(nd1->se_x-prev_node->nw_x);
					//freed_node = prev_node;
					delete_node(prev_node);
				}
			}
		}
	}
	if( dir == 0 ) nd1 = tree;
	else if( dir == 1 ) nd1 = temp->left;
	else if( dir == 2 ) nd1 = temp-> right;
	prev_node = nd1;
	if(nd1)
	{
		if(nd1->right) clean_it(nd1->right,2);
	}
}

int insert_counter = 0;
int insert_counter2 = 0;
int rebalance_count = 100;
int rebuild_count = 0;
int another_counter = 0;
/*int del_root = 0;*/

int semi_mock_insert_return_val = -1;

int semi_mock_insert2(int type, double nw_x, double nw_y, double se_x, double se_y, double slope, struct node **leaf) // This function mimics the actions of the insert function, takes an input point and determines whether or not that point should be added to the tree. If it should be, it is added and 1 is returned, otherwise 0 is returned.
{	
/*	printf("******************************************\n");*/
/*	print_preorder(tree,NULL);*/
/*	printf("******************************************\n");*/
/*	printf("------------------------------------------\n");*/
/*	print_inorder(tree);*/
/*	printf("------------------------------------------\n");*/
/*	printf("inserted data is: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",type,nw_x,nw_y,se_x,se_y);*/
/*	printf("or maybe it was: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",type,x_ideal-nw_x,y_ideal-nw_y,x_ideal-se_x,y_ideal-se_y);*/
/*	printf("slope: %lf\t calculated: %lf\n",slope,(nw_y - se_y)/(nw_x-se_x));*/
/*	node *temp_node3 = tree;*/
/*	if(*leaf){*/
/*	printf("comparing against: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",(*leaf)->type,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y);*/
/*	printf("or ...: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",(*leaf)->type,x_ideal-(*leaf)->nw_x,y_ideal-(*leaf)->nw_y,x_ideal-(*leaf)->se_x,y_ideal-(*leaf)->se_y);*/
/*	}*/
/*	if (*leaf && (*leaf)->type == 2 && (*leaf)->nw_x < 5.377857 && (*leaf)->nw_x > 5.37785 )*/
/*	{*/
/*		do_nothing(1);*/
/*	}*/
/*	if(insert_counter2 > 20) exit(0);*/
	semi_mock_insert_return_val = 0;
	if((*leaf) == tree || insert_to_potential_branch_tree == 1)
	{
		insert_to_potential_branch_tree = 0;
		if(type == 1) // Ensures points are dealt with appropriately
		{
			se_x = nw_x;
			se_y = nw_y;
		}
		if(del_root == 0)
		{	
			nw_x =  x_ideal - nw_x;
			se_x =  x_ideal - se_x;
			nw_y =  y_ideal - nw_y;
			se_y =  y_ideal - se_y;
			insert_counter2++;
		}
		else del_root = 0;
		if(type == 2) slope = (se_y-nw_y)/(se_x-nw_x);
	}
	//printf("after messing: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",type,nw_x,nw_y,se_x,se_y);
	if(nw_x > se_x || se_y > nw_y)
	{
		return 0;
	}
	else if(type == 2 && (double) nw_x == (double) se_x)
	{
		type = 1;
		nw_y = se_y;
	}
	else if(type == 2 && (double) nw_y == (double) se_y)
	{
		type = 1;
		se_x = nw_x;
	}
	if(*leaf && type == 1 && ((fabs(nw_x - (*leaf)->nw_x) < .000001 && fabs(nw_y - (*leaf)->nw_y) < .000001) || (fabs(nw_x - (*leaf)->se_x) < .000001 && fabs(nw_y - (*leaf)->se_y) < .000001)))
	{	//printf("point already in tree\n");
		return 0;
	}
	else if(*leaf && type == 2 && (fabs(nw_x - (*leaf)->nw_x) < .000001 && fabs(nw_y - (*leaf)->nw_y) < .000001 && fabs(se_x - (*leaf)->se_x) < .000001 && fabs(se_y - (*leaf)->se_y) < .000001))
	{	//printf("segment already in tree\n");
		return 0;
	}
	if(*leaf && *leaf == tree)
	{
		if(insert_counter > rebalance_count)
		{
			rebalance_count = round(1.00733*insert_counter);	//use 1.00733 if rebalancing only using periodic full tree, use 8 if combining with traversed path rebalance
			Rebalance(tree);
			another_counter++;
			prev_node = NULL;
/*			printf("before cleaning, size is %d\n",tree->subtree_size);*/
/*			start_time = clock();*/
/*			clean_it(tree,0);*/
/*			finish_time = clock();*/
/*			printf("after cleaning, size is %d and the time to do it was %lf\n",tree->subtree_size,(double)(finish_time - start_time) / CLOCKS_PER_SEC);*/
		}
	}
//	if(*leaf && (*leaf)->subtree_size > 2)
//	{
//		Rebalance2(*leaf);
//	}
/*	if(tree)						//These lines can be used for debugging*/
/*	{*/
/*		scan_proper_tree(tree);*/
/*	}*/
    	if( *leaf == NULL) // Node where we are trying to insert does not exist. This means that a new *leaf must be created.
   	{	//printf("creating new node: %d \t %.12lf,%.12lf  %.12lf,%.12lf\n",type,nw_x,nw_y,se_x,se_y);
      		*leaf = (struct node*) malloc( sizeof( struct node ) ); 
        	(*leaf)->type = type;
		(*leaf)->nw_x = nw_x;
		(*leaf)->se_x = se_x;
		(*leaf)->nw_y = nw_y;
		(*leaf)->se_y = se_y;
		(*leaf)->slope = slope;
		(*leaf)->subtree_size = 1;
		(*leaf)->parent = NULL;
        	(*leaf)->left = NULL;    
        	(*leaf)->right = NULL;
        	insert_counter++;
        	return 1;
	}
	else // The node being compared to does exist
	{	
		if(type == 1)
		{
			if((*leaf)->type == 1)
			{
				if(nw_x <= (*leaf)->nw_x && nw_y <= (*leaf)->se_y) // The point is dominated
				{	node *temp_node = NULL;
					if((*leaf)->left) temp_node = (*leaf)->left;
					else if((*leaf)->right) temp_node = (*leaf)->right;
					if((*leaf)->subtree_size == 1)
					{
						(*leaf)->type = type;
						(*leaf)->nw_x = nw_x;
						(*leaf)->se_x = se_x;
						(*leaf)->nw_y = nw_y;
						(*leaf)->se_y = se_y;
						(*leaf)->slope = slope;
						insert_counter++;
						return 1;
					}
					else
					{
						if(*leaf == tree) del_root = 1;
						delete_node(*leaf);
						//printf("inserting 1\n");
						semi_mock_insert_return_val = semi_mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf));
					}
				}
				else if(nw_x < (*leaf)->nw_x && nw_y >= (*leaf)->se_y) // The point needs inserted left
				{	//printf("inserting 2\n");
					semi_mock_insert_return_val = semi_mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left);
					if((*leaf)->left)
					{
						(*leaf)->left->parent = *leaf;
						if((*leaf)->left->subtree_size == 1) update_depth(*leaf);
					}
				}
				else if(nw_x >= (*leaf)->nw_x && nw_y < (*leaf)->se_y) // The point needs inserted right
				{	//printf("inserting 3\n");
					semi_mock_insert_return_val = semi_mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right);
					if((*leaf)->right)
					{
						(*leaf)->right->parent = *leaf;
						if((*leaf)->right->subtree_size == 1) update_depth(*leaf);
					}
				}
			}	
			else
			{
				if(nw_x <= (*leaf)->nw_x && nw_y >= (*leaf)->nw_y) // The point needs inserted left
				{	//printf("inserting 4\n");
					semi_mock_insert_return_val = semi_mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left);
					if((*leaf)->left)
					{
						(*leaf)->left->parent = *leaf;
						if((*leaf)->left->subtree_size == 1) update_depth(*leaf);
					}
				}
				else if(nw_x >= (*leaf)->se_x && nw_y <= (*leaf)->se_y) // The point needs inserted right
				{	//printf("inserting 5\n");
					semi_mock_insert_return_val = semi_mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right);
					if((*leaf)->right)
					{
						(*leaf)->right->parent = *leaf;
						if((*leaf)->right->subtree_size == 1) update_depth(*leaf);
					}
				}
				else
				{
					double x_proj = (*leaf)->se_x + (1./(*leaf)->slope)*(nw_y - (*leaf)->se_y);
					double y_proj = (*leaf)->se_y + (*leaf)->slope*(nw_x - (*leaf)->se_x);
					//printf("incoming: %lf,%lf to %lf,%lf\t current: %lf,%lf to %lf,%lf\t projection:%lf,%lf\n",nw_x,nw_y,se_x,se_y,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y,x_proj,y_proj);
					if(nw_x < x_proj || nw_y < y_proj) // The point is under (the extension of) the segment
					{
						if( nw_x <= (*leaf)->nw_x && x_proj >= (*leaf)->se_x ) // The entire segment is dominated
						{	node *temp_node = NULL;
							if((*leaf)->left) temp_node = (*leaf)->left;
							else if((*leaf)->right) temp_node = (*leaf)->right;
							if((*leaf)->subtree_size == 1)
							{
								(*leaf)->type = type;
								(*leaf)->nw_x = nw_x;
								(*leaf)->se_x = se_x;
								(*leaf)->nw_y = nw_y;
								(*leaf)->se_y = se_y;
								(*leaf)->slope = slope;
								insert_counter++;
								return 1;
							}
							else
							{
								if(*leaf == tree) del_root = 1;
								delete_node(*leaf);
								//if(printing) printf("inserting 6\n");
								semi_mock_insert_return_val = semi_mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf));
							}
						}
						else if( nw_x <= (*leaf)->nw_x ) // The left portion of the segment is dominated
						{
							(*leaf)->nw_x = x_proj;
							(*leaf)->nw_y = nw_y;
							//if(printing) printf("inserting 7\n");
							semi_mock_insert_return_val = semi_mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left);
							if((*leaf)->left)
							{
								(*leaf)->left->parent = *leaf;
								if((*leaf)->left->subtree_size == 1) update_depth(*leaf);
							}
						}
						else if( x_proj >= (*leaf)->se_x ) // The right portion of the segment is dominated
						{
							(*leaf)->se_y = (*leaf)->se_y + ((*leaf)->se_y - (*leaf)->nw_y)/((*leaf)->se_x - (*leaf)->nw_x)*(nw_x - (*leaf)->se_x);
							(*leaf)->se_x = nw_x;
							//printf("inserting 8\n");
							semi_mock_insert_return_val = semi_mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right);
							if((*leaf)->right)
							{
								(*leaf)->right->parent = *leaf;
								if((*leaf)->right->subtree_size == 1) update_depth(*leaf);
							}
						}
						else // The center portion of the segment is dominated
						{
							double save1 = (*leaf)->se_x;
							double save2 = (*leaf)->se_y;
							(*leaf)->se_y = (*leaf)->se_y + (*leaf)->slope*(nw_x - (*leaf)->se_x);
							(*leaf)->se_x = nw_x;
							//printf("inserting 8b\n");
							semi_mock_insert2((*leaf)->type,x_proj,nw_y,save1,save2,(*leaf)->slope,&(*leaf)->right);
							if((*leaf)->right)
							{	
								(*leaf)->right->parent = *leaf;

								if((*leaf)->right->subtree_size == 1) update_depth(*leaf);
							}
							//printf("inserting 9\n");
							semi_mock_insert_return_val = semi_mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right);
							if((*leaf)->right)
							{
								(*leaf)->right->parent = *leaf;
								if((*leaf)->right->subtree_size == 1) update_depth(*leaf);
							}
						}
					}
				}
			}
		}
	}
	return semi_mock_insert_return_val;
}

int mock_insert_return_val = -1, mock_insert_counter = 0;

int mock_insert(int type, double nw_x, double nw_y, double se_x, double se_y, double slope, struct node **leaf)
{
	int ret_val;
	
	start_struct_timer = clock();
	if( clock_gettime( CLOCK_REALTIME, &start2) == -1 ) {
      		perror( "clock gettime" );
      		exit( EXIT_FAILURE );
    		}
	ret_val = mock_insert2(type, nw_x, nw_y, se_x, se_y, slope, &(*leaf));
	finish_struct_timer = clock();
	if( clock_gettime( CLOCK_REALTIME, &stop2) == -1 ) {
      		perror( "clock gettime" );
      		exit( EXIT_FAILURE );
	    	}

	    	accum = ( stop2.tv_sec - start2.tv_sec )
		  	+ ( stop2.tv_nsec - start2.tv_nsec )
		    	/ BILLION;
	     	struct_time2 += accum;
        struct_time += (double)(finish_struct_timer - start_struct_timer) / CLOCKS_PER_SEC;
	return ret_val;
}

int mock_insert2(int type, double nw_x, double nw_y, double se_x, double se_y, double slope, struct node **leaf) // This function mimics the actions of the insert function, takes an input point and determines whether or not that point would be added to the tree if so desired. Will be used for a fathoming rule.
{	
/*	printing = 1;*/
	if(type == 2 && there_will_only_be_points && integer_bb)
	{
/*		printf("integer_objective_lb: %lf\n",integer_objective_lb);*/
		if(integer_objective == 1)
		{
			double multiplier = floor( (nw_x - integer_objective_lb)/smallest_coef );
			double upper_val = integer_objective_lb + multiplier*smallest_coef;
			int rval = 0, it = 0;
/*			printf("nw_x: %lf, se_x: %lf, upper_val: %lf\n",nw_x,se_x,upper_val);*/
			while(upper_val > se_x && rval == 0)
			{
				it++;
				rval = mock_insert2(1,upper_val,(nw_y-se_y)/(nw_x-se_x)*(upper_val-nw_x)+nw_y,0.,0.,0.,&tree);
/*				printf("plot(%lf,%lf,'mo');\n",upper_val,(nw_y-se_y)/(nw_x-se_x)*(upper_val-nw_x)+nw_y);*/
				upper_val -= smallest_coef;
				if(rval == 1) return rval;
			}
			multiplier = floor( (se_x - integer_objective_lb)/smallest_coef );
			upper_val = integer_objective_lb + multiplier*smallest_coef;
/*			printf("plot(%lf,%lf,'mo');\n",upper_val,se_y);*/
			rval = mock_insert2(1,upper_val,se_y,0.,0.,0.,&tree);
			if(rval == 1) return rval;
			else return 0;
		}
		else if(integer_objective == 2)
		{
			double multiplier = floor( (se_y - integer_objective_lb)/smallest_coef );
			double upper_val = integer_objective_lb + multiplier*smallest_coef;
			int rval = 0, it = 0;
/*			printf("nw_y: %lf, se_y: %lf, upper_val: %lf\n",nw_y,se_y,upper_val);*/
			while(upper_val > nw_y && rval == 0)
			{
				it++;
				rval = mock_insert2(1,(nw_x-se_x)/(nw_y-se_y)*(upper_val-nw_y)+nw_x,upper_val,0.,0.,0.,&tree);
/*				printf("plot(%lf,%lf,'mo');\n",(nw_x-se_x)/(nw_y-se_y)*(upper_val-nw_y)+nw_x,upper_val);*/
				upper_val -= smallest_coef;
				if(rval == 1) return rval;
			}
			multiplier = floor( (nw_y - integer_objective_lb)/smallest_coef );
			upper_val = integer_objective_lb + multiplier*smallest_coef;
/*			printf("plot(%lf,%lf,'mo');\n",upper_val,se_y);*/
			rval = mock_insert2(1,nw_x,upper_val,0.,0.,0.,&tree);
			if(rval == 1) return rval;
			else return 0;
		}
	}
	if(printing)
	{
	/*	printf("inserted data is: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",type,nw_x,nw_y,se_x,se_y);*/
	/*	printf("or maybe it was: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",type,x_ideal-nw_x,y_ideal-nw_y,x_ideal-se_x,y_ideal-se_y);*/
		if(type == 1)
		{
	/*		printf("_$^&__$^&__$^&__$^&__$^&__$^&__$^&_\n");*/
	/*		print_inorder(tree);*/
	/*		printf("_$^&__$^&__$^&__$^&__$^&__$^&__$^&_\n");*/
/*			printf("plot([%lf],[%lf],'-o');\n",nw_x,nw_y);*/
/*			printf("plot([%lf],[%lf],'-o');\n",x_ideal-nw_x,y_ideal-nw_y);*/
		}
		else
		{
/*			printf("plot([%lf,%lf],[%lf,%lf],'-o');\n",nw_x,se_x,nw_y,se_y);*/
/*			printf("plot([%lf,%lf],[%lf,%lf],'-o');\n",x_ideal-nw_x,x_ideal-se_x,y_ideal-nw_y,y_ideal-se_y);*/
		}
	}
	//printf("slope: %lf\t calculated: %lf\n",slope,(nw_y - se_y)/(nw_x-se_x));
	mock_insert_return_val = 0;
/*	node *temp_node3 = tree;*/
	if(printing && *leaf)
	{
	/*	printf("comparing against: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",(*leaf)->type,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y);*/
/*		if((*leaf)->type == 1) */
/*		{*/
/*			printf("plot([%lf],[%lf],'-go');\n",x_ideal-(*leaf)->nw_x,y_ideal-(*leaf)->nw_y);*/
/*			printf("plot([%lf],[%lf],'-go');\n",(*leaf)->nw_x,(*leaf)->nw_y);*/
/*		}*/
/*		else	*/
/*		{*/
/*			printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-(*leaf)->nw_x,x_ideal-(*leaf)->se_x,y_ideal-(*leaf)->nw_y,y_ideal-(*leaf)->se_y);*/
/*			printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",(*leaf)->nw_x,(*leaf)->se_x,(*leaf)->nw_y,(*leaf)->se_y);*/
/*		}*/
	}
/*	if (*leaf && (*leaf)->type == 2 && (*leaf)->nw_x < 5.377857 && (*leaf)->nw_x > 5.37785 )*/
/*	{*/
/*		do_nothing(1);*/
/*	}*/
	if((*leaf) == tree || insert_to_potential_branch_tree == 1)
	{
		insert_to_potential_branch_tree = 0;
		if(type == 1) // Ensures points are dealt with appropriately
		{
			se_x = nw_x;
			se_y = nw_y;
		}
		nw_x =  x_ideal - nw_x;
		se_x =  x_ideal - se_x;
		nw_y =  y_ideal - nw_y;
		se_y =  y_ideal - se_y;
		if(type == 2) slope = (se_y-nw_y)/(se_x-nw_x);
		mock_insert_counter++;
	}
	//printf("after messing: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",type,nw_x,nw_y,se_x,se_y);
	if(nw_x > se_x || se_y > nw_y)
	{
/*		if(printing)*/
/*		 printf("the segment is backwards\n");*/
		return 0;
	}
	else if(type == 2 && (double) nw_x == (double) se_x)
	{
		type = 1;
		nw_y = se_y;
	}
	else if(type == 2 && (double) nw_y == (double) se_y)
	{
		type = 1;
		se_x = nw_x;
	}
	if(*leaf && type == 1 && ((fabs(nw_x - (*leaf)->nw_x) < .00001 && fabs(nw_y - (*leaf)->nw_y) < .00001) || (fabs(nw_x - (*leaf)->se_x) < .00001 && fabs(nw_y - (*leaf)->se_y) < .00001)))
	{	
/*		printf("point already in tree\n");*/
		return 0;
	}
	else if(*leaf && type == 2 && (fabs(nw_x - (*leaf)->nw_x) < .00001 && fabs(nw_y - (*leaf)->nw_y) < .00001 && fabs(se_x - (*leaf)->se_x) < .00001 && fabs(se_y - (*leaf)->se_y) < .00001))
	{	
		if(printing) printf("segment already in tree\n");
		return 0;
	}
	if( !(*leaf) || (*leaf) == NULL) return 1;
/*	if(*leaf && *leaf == tree)*/
/*	{*/
/*		if(insert_counter > rebalance_count)*/
/*		{*/
/*			rebalance_count = round(1.00733*insert_counter);	//use 1.00733 if rebalancing only using periodic full tree, use 8 if combining with traversed path rebalance*/
/*			Rebalance(tree);*/
/*			another_counter++;*/
/*		}*/
/*	}*/
//	if(*leaf && (*leaf)->subtree_size > 2)
//	{
//		Rebalance2(*leaf);
//	}
/*	if(tree)						//These lines can be used for debugging*/
/*	{*/
/*		scan_proper_tree(tree);*/
/*	}*/
/*    	if( *leaf == NULL) // Node where we are trying to insert does not exist. This means that a new *leaf must be created.*/
/*   	{	//printf("creating new node: %d \t %.12lf,%.12lf  %.12lf,%.12lf\n",type,nw_x,nw_y,se_x,se_y);*/
/*      		*leaf = (struct node*) malloc( sizeof( struct node ) ); */
/*        	(*leaf)->type = type;*/
/*		(*leaf)->nw_x = nw_x;*/
/*		(*leaf)->se_x = se_x;*/
/*		(*leaf)->nw_y = nw_y;*/
/*		(*leaf)->se_y = se_y;*/
/*		(*leaf)->slope = slope;*/
/*		(*leaf)->subtree_size = 1;*/
/*		(*leaf)->parent = NULL;*/
/*        	(*leaf)->left = NULL;    */
/*        	(*leaf)->right = NULL;*/
/*        	insert_counter++;*/
/*	}*/
	double distance = 0.;
	if(*leaf != NULL) // The node being compared to does exist
	{	
		if(type == 1) // The inserted solution is a point.
		{
			if((*leaf)->type == 1) // The current solution is a point.
			{
				if(use_hausdorff && epsilon > 0.)
				{
					distance = sqrt( (nw_x - (*leaf)->nw_x)*(nw_x - (*leaf)->nw_x) + (nw_y - (*leaf)->nw_y)*(nw_y - (*leaf)->nw_y) );
/*					printf("distance of %lf btwn pts\n", distance);*/
					if(distance < epsilon*max_range) 
					{
/*						printf("distance of %lf btwn pts did not warrant continuing\n", distance);*/
						return 0;
					}
				}
				
				if(nw_x - (*leaf)->nw_x <= .00000001 && nw_y - (*leaf)->se_y <= .000000001 ) // The current point is dominated
				{	
/*					printf("cur pt dom\n");*/
					return 1; 
				}
				else if(nw_x < (*leaf)->nw_x && nw_y >= (*leaf)->se_y) // The point needs inserted left
				{	
/*					printf("going left\n");*/
					mock_insert_return_val = mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left);
				}
				else if(nw_x >= (*leaf)->nw_x && nw_y < (*leaf)->se_y) // The point needs inserted right
				{	
/*					printf("going right\n");*/
					mock_insert_return_val = mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right);
				}
			}	
			else // The current solution is a segment
			{
				if(use_hausdorff && epsilon > 0.)
				{
					distance = fmax( sqrt( (nw_x - (*leaf)->nw_x)*(nw_x - (*leaf)->nw_x) + (nw_y - (*leaf)->nw_y)*(nw_y - (*leaf)->nw_y) ),
						sqrt( (nw_x - (*leaf)->se_x)*(nw_x - (*leaf)->se_x) + (nw_y - (*leaf)->se_y)*(nw_y - (*leaf)->se_y) ) );
/*					printf("distance of %lf btwn pt and segment\n", distance);*/
					if(distance < epsilon*max_range)
					{
/*						printf("distance of %lf btwn pt and segment did not warrant continuing\n", distance);*/
						return 0;
					}
				}
				
				if(nw_x <= (*leaf)->nw_x && nw_y >= (*leaf)->nw_y) // The point needs inserted left
				{	
/*					printf("going left\n");*/
					mock_insert_return_val = mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left);
				}
				else if(nw_x >= (*leaf)->se_x && nw_y <= (*leaf)->se_y) // The point needs inserted right
				{	
/*					printf("going right\n");*/
					mock_insert_return_val = mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right);
				}
				else
				{
					double x_proj = (*leaf)->se_x + (1./(*leaf)->slope)*(nw_y - (*leaf)->se_y);
					double y_proj = (*leaf)->se_y + (*leaf)->slope*(nw_x - (*leaf)->se_x);
/*					printf("incoming: %lf,%lf to %lf,%lf\t current: %lf,%lf to %lf,%lf\t projection:%lf,%lf\n",nw_x,nw_y,se_x,se_y,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y,x_proj,y_proj);*/
					if(nw_x < x_proj - 0.00001 || nw_y < y_proj - 0.00001) // The point is under (the extension of) the segment
					{
/*						printf("pt under sgmt\n");*/
						if(use_hausdorff && epsilon > 0.)
						{
							double temp_nw_x = fmax(nw_x,(*leaf)->nw_x);
							double temp_nw_y = (*leaf)->nw_y;
							if(temp_nw_x != (*leaf)->nw_x) temp_nw_y = (*leaf)->se_y + (*leaf)->slope*(temp_nw_x - (*leaf)->se_x);
							double temp_se_y = fmax(se_y,(*leaf)->se_y);
							double temp_se_x = (*leaf)->se_x;
							if(temp_se_y != (*leaf)->se_y) temp_se_x = (*leaf)->se_x + (1./(*leaf)->slope)*(temp_se_y - (*leaf)->se_y);
							
/*							printf("plot([%lf],[%lf],'-ro');\n",temp_nw_x,temp_nw_y);*/
/*							printf("plot([%lf],[%lf],'-ro');\n",temp_se_x,temp_se_y);*/
							
							distance = fmax( sqrt( (nw_x - temp_nw_x)*(nw_x - temp_nw_x) + (nw_y - temp_nw_y)*(nw_y - temp_nw_y) ),
								sqrt( (nw_x - temp_se_x)*(nw_x - temp_se_x) + (nw_y - temp_se_y)*(nw_y - temp_se_y) ) );
/*							printf("distance of %lf btwn pt and segment\n", distance);*/
							if(distance < epsilon*max_range)
							{
		/*						printf("distance of %lf btwn pt and segment did not warrant continuing\n", distance);*/
								return 0;
							}
							else return 1;
						}
						else return 1;
					}
				}
			}
		}
		else // Input is a line segment
		{
			int returned_val = 0;
			if((*leaf)->type == 1) // Incoming segment compared against current point
			{
				if(use_hausdorff && epsilon > 0.)
				{
					distance = fmax( sqrt( (nw_x - (*leaf)->nw_x)*(nw_x - (*leaf)->nw_x) + (nw_y - (*leaf)->nw_y)*(nw_y - (*leaf)->nw_y) ),
						sqrt( (se_x - (*leaf)->nw_x)*(se_x - (*leaf)->nw_x) + (se_y - (*leaf)->nw_y)*(se_y - (*leaf)->nw_y) ) );
/*					printf("distance of %lf btwn pt and segment\n",distance);*/
				 	if(distance < epsilon*max_range)
				 	{
/*						printf("distance of %lf btwn pt and segment did not warrant continuing\n",distance);*/
						return 0;
					}
				}
			
				double y_proj = se_y + slope*((*leaf)->nw_x - se_x);
				double x_proj = se_x + (1./slope)*((*leaf)->nw_y - se_y);
				//printf("incoming: %lf,%lf to %lf,%lf\t current: %lf,%lf to %lf,%lf\n",nw_x,nw_y,se_x,se_y,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y);
				//printf("the projections are: %.12lf,%.12lf and %.12lf,%.12lf\n",x_proj,(*leaf)->nw_y,(*leaf)->nw_x,y_proj);
				//printf("slope: %lf\t calculated: %lf\n",slope,(se_y - nw_y)/(se_x - nw_x));
				if(se_x - (*leaf)->nw_x <= .0000001 && se_y - (*leaf)->nw_y >= -.00000001) // The segment needs inserted left
				{	//printf("inserting 10\n");
					if(printing) printf("location 1\n");
					returned_val = mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left);
					if(returned_val == 1) return 1;
					mock_insert_return_val = max(mock_insert_return_val,returned_val);
				}
				else if(nw_x - (*leaf)->nw_x >= -.00000001 && nw_y - (*leaf)->nw_y <= .000000001 ) // The segment needs inserted right
				{	//printf("inserting 11\n");
					if(printing) printf("location 2\n");
					returned_val = mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right);
					if(returned_val == 1) return 1;
					mock_insert_return_val = max(mock_insert_return_val,returned_val);
				}
				else if(x_proj <= ((*leaf)->nw_x + .0001) || y_proj <= ((*leaf)->nw_y + .0001)) // The point is dominated
				{	
					if(printing) printf("location 3\n");
					if(use_hausdorff && epsilon > 0.)
					{
						double temp_nw_y = fmin(nw_y,(*leaf)->nw_y);
						double temp_nw_x = nw_x;
						if(temp_nw_y != nw_y) temp_nw_x = se_x + (1./slope)*(temp_nw_y - se_y);
						double temp_se_x = fmin(se_x,(*leaf)->se_x);
						double temp_se_y = se_y;
						if(temp_se_x != se_x) temp_se_y = se_y + slope*(temp_se_x - se_x);
						
/*						printf("plot([%lf],[%lf],'-ro');\n",temp_nw_x,temp_nw_y);*/
/*						printf("plot([%lf],[%lf],'-ro');\n",temp_se_x,temp_se_y);*/
						
						distance = fmax( sqrt( (temp_nw_x - (*leaf)->nw_x)*(temp_nw_x - (*leaf)->nw_x) + 
							(temp_nw_y - (*leaf)->nw_y)*(temp_nw_y - (*leaf)->nw_y) ),
							sqrt( (temp_se_x - (*leaf)->nw_x)*(temp_se_x - (*leaf)->nw_x) + 
							(temp_se_y - (*leaf)->nw_y)*(temp_se_y - (*leaf)->nw_y) ) );
/*						printf("distance of %lf btwn pt and segment\n",distance);*/
					 	if(distance < epsilon*max_range)
					 	{
					 		int yes = 0;
					 		if(temp_nw_x > nw_x)
					 		{
					 			yes = 1;
						 		returned_val = mock_insert2(type,nw_x,nw_y,temp_nw_x,temp_nw_y,slope,&(*leaf)->left);
								if(returned_val == 1) return 1;
								mock_insert_return_val = max(mock_insert_return_val,returned_val);
							}
							if(temp_se_x < se_x)
							{
								yes = 1;
								returned_val = mock_insert2(type,temp_se_x,temp_se_y,se_x,se_y,slope,&(*leaf)->right);
								if(returned_val == 1) return 1;
								mock_insert_return_val = max(mock_insert_return_val,returned_val);
							}
							if(!yes) return 1;
						}
						else return 1;
					}
					else return 1;
				}
				else 
				{
/*					printf("location 4\n");*/
					if(nw_x != (*leaf)->nw_x)
					{
					returned_val = mock_insert2(type,nw_x,nw_y,(*leaf)->nw_x,y_proj,slope,&(*leaf)->left);
					if(returned_val == 1) return 1;
					mock_insert_return_val = max(mock_insert_return_val,returned_val);
					}
					if((*leaf)->nw_y != se_y)
					{
					returned_val = mock_insert2(type,x_proj,(*leaf)->nw_y,se_x,se_y,slope,&(*leaf)->right);
					if(returned_val == 1) return 1;
					mock_insert_return_val = max(mock_insert_return_val,returned_val);
					}
				}
			}	
			else // Incoming segment compared against current segment
			{	
				if(use_hausdorff && epsilon > 0.)
				{
/*					printf("plot([%lf,%lf],[%lf,%lf],'-o');\n",nw_x,se_x,nw_y,se_y);*/
/*					printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",(*leaf)->nw_x,(*leaf)->se_x,(*leaf)->nw_y,(*leaf)->se_y);*/
/*					printf("the supposed intersections of projections:\n");*/
				
					double dist1 = 0., dist2 = 0, dist3 = 0., dist4 = 0.;
					double x_intersect = ((1./(*leaf)->slope)*se_x+(*leaf)->slope*(*leaf)->se_x+se_y-(*leaf)->se_y)/
													((*leaf)->slope+(1./(*leaf)->slope)); // proj from se_x
					double y_intersect = (*leaf)->nw_y + (*leaf)->slope*(x_intersect - (*leaf)->nw_x);
					if(y_intersect != y_intersect) y_intersect = se_y + (-1./(*leaf)->slope)*(x_intersect - se_x);
					
/*					printf("plot([%lf],[%lf],'-ro');\n",x_intersect,y_intersect);*/
				
					if(x_intersect <= (*leaf)->se_x && x_intersect >= (*leaf)->nw_x) dist1 = 
						sqrt( (se_x - x_intersect)*(se_x - x_intersect) + (se_y - y_intersect)*(se_y - y_intersect) );
					else dist1 = fmin( sqrt( (se_x - (*leaf)->nw_x)*(se_x - (*leaf)->nw_x) + (se_y - (*leaf)->nw_y)*(se_y - (*leaf)->nw_y) ),
							sqrt( (se_x - (*leaf)->se_x)*(se_x - (*leaf)->se_x) + (se_y - (*leaf)->se_y)*(se_y - (*leaf)->se_y) ) );
						
					x_intersect = ((1./(*leaf)->slope)*nw_x+(*leaf)->slope*(*leaf)->se_x+nw_y-(*leaf)->se_y)/
													((*leaf)->slope+(1./(*leaf)->slope)); // proj from nw_x
					y_intersect = (*leaf)->nw_y + (*leaf)->slope*(x_intersect - (*leaf)->nw_x);
					if(y_intersect != y_intersect) y_intersect = nw_y + (-1./(*leaf)->slope)*(x_intersect - nw_x);
					
/*					printf("plot([%lf],[%lf],'-mo');\n",x_intersect,y_intersect);*/
				
					if(x_intersect <= (*leaf)->se_x && x_intersect >= (*leaf)->nw_x) dist2 = 
						sqrt( (nw_x - x_intersect)*(nw_x - x_intersect) + (nw_y - y_intersect)*(nw_y - y_intersect) );
					else dist2 = fmin( sqrt( (nw_x - (*leaf)->nw_x)*(nw_x - (*leaf)->nw_x) + (nw_y - (*leaf)->nw_y)*(nw_y - (*leaf)->nw_y) ),
							sqrt( (nw_x - (*leaf)->se_x)*(nw_x - (*leaf)->se_x) + (nw_y - (*leaf)->se_y)*(nw_y - (*leaf)->se_y) ) );
						
					x_intersect = ((slope)*nw_x+(1./slope)*(*leaf)->se_x-nw_y+(*leaf)->se_y)/((1./slope)+(slope)); // proj from leaf se_x
					y_intersect = (*leaf)->se_y + (-1./slope)*(x_intersect - (*leaf)->se_x);
					if(y_intersect != y_intersect) y_intersect = nw_y + (slope)*(x_intersect - nw_x);
					
/*					printf("plot([%lf],[%lf],'-ko');\n",x_intersect,y_intersect);*/
				
					if(x_intersect <= se_x && x_intersect >= nw_x) dist3 = sqrt( ((*leaf)->se_x - x_intersect)*((*leaf)->se_x - x_intersect) +
												((*leaf)->se_y - y_intersect)*((*leaf)->se_y - y_intersect) );
					else dist3 = fmin( sqrt( (nw_x - (*leaf)->se_x)*(nw_x - (*leaf)->se_x) + (nw_y - (*leaf)->se_y)*(nw_y - (*leaf)->se_y) ),
							sqrt( (se_x - (*leaf)->se_x)*(se_x - (*leaf)->se_x) + (se_y - (*leaf)->se_y)*(se_y - (*leaf)->se_y) ) );
						
					x_intersect = ((slope)*nw_x+(1./slope)*(*leaf)->nw_x-nw_y+(*leaf)->nw_y)/((1./slope)+(slope)); // proj from leaf nw_x
					y_intersect = (*leaf)->nw_y + (-1./slope)*(x_intersect - (*leaf)->nw_x);
					if(y_intersect != y_intersect) y_intersect = nw_y + (slope)*(x_intersect - nw_x);
					
/*					printf("plot([%lf],[%lf],'-o');\n",x_intersect,y_intersect);*/
				
					if(x_intersect <= se_x && x_intersect >= nw_x) dist4 = sqrt( ((*leaf)->nw_x - x_intersect)*((*leaf)->nw_x - x_intersect) +
												((*leaf)->nw_y - y_intersect)*((*leaf)->nw_y - y_intersect) );
					else dist4 = fmin( sqrt( (nw_x - (*leaf)->nw_x)*(nw_x - (*leaf)->nw_x) + (nw_y - (*leaf)->nw_y)*(nw_y - (*leaf)->nw_y) ),
							sqrt( (se_x - (*leaf)->nw_x)*(se_x - (*leaf)->nw_x) + (se_y - (*leaf)->nw_y)*(se_y - (*leaf)->nw_y) ) );
						
					dist1 = fmax(dist1,dist2);
					dist3 = fmax(dist3,dist4);
					distance = fmax(dist1,dist3);
					
/*					printf("the max distance: %lf\n",distance);*/
					
/*					printf("distance of %lf btwn segments\n", distance);*/
					if(distance < epsilon*max_range)
					{
/*						printf("distance of %lf btwn segments did not warrant continuing\n", distance);*/
						return 0;
					}
				}
			
				if(se_x <= (*leaf)->nw_x + .00001 && se_y >= (*leaf)->nw_y - .00001) // The segment needs inserted left
				{	
					if(printing) printf("location 5\n");
					//printf("inserting 15\n");
					returned_val = mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left);
					if(returned_val == 1) return 1;
					mock_insert_return_val = max(mock_insert_return_val,returned_val);
				}
				else if(nw_x - (*leaf)->se_x >= -.00000001 && nw_y - (*leaf)->se_y <= .00000001) // The segment needs inserted right
				{	//printf("incoming segment: %lf,%lf to %lf,%lf \t current segment: %lf,%lf to %lf,%lf\n",nw_x,nw_y,se_x,se_y,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y);
					//printf("inserting 16\n");
					if(printing) printf("location 6\n");
					returned_val = mock_insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right);
					if(returned_val == 1) return 1;
					mock_insert_return_val = max(mock_insert_return_val,returned_val);
				}
				else
				{	
					if(printing) printf("location 7\n");
					//printf("incoming segment: %.12lf,%.12lf to %.12lf,%.12lf \t current segment: %.12lf,%.12lf to %.12lf,%.12lf\n",nw_x,nw_y,se_x,se_y,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y);
/*					double x_intersect = ((*leaf)->nw_y*se_x*(*leaf)->se_x-(*leaf)->nw_y*nw_x*(*leaf)->se_x-nw_y*se_x*(*leaf)->se_x+nw_y*se_x*(*leaf)->nw_x+se_y*nw_x*(*leaf)->se_x-se_y*nw_x*(*leaf)->nw_x-(*leaf)->se_y*se_x*(*leaf)->nw_x+(*leaf)->se_y*nw_x*(*leaf)->nw_x)/(-(*leaf)->se_y*se_x+(*leaf)->se_y*nw_x+(*leaf)->nw_y*se_x-(*leaf)->nw_y*nw_x+se_y*(*leaf)->se_x-se_y*(*leaf)->nw_x-nw_y*(*leaf)->se_x+nw_y*(*leaf)->nw_x); // Find the intersection of (the possible extensions of) the line segments being compared*/
					double x_intersect = (slope*se_x-(*leaf)->slope*(*leaf)->se_x-se_y+(*leaf)->se_y)/(-(*leaf)->slope+slope);
					double y_intersect = (*leaf)->nw_y + (*leaf)->slope*(x_intersect - (*leaf)->nw_x);
					if(y_intersect != y_intersect) y_intersect = nw_y + slope*(x_intersect - nw_x);
					//printf("the intersection is: %.12lf,%.12lf\n",x_intersect,y_intersect);
					double nw_y_proj_c2i = se_y + slope*((*leaf)->nw_x - se_x); // projection of the endpoints of current segment onto incoming segment (hence, c2i)
					double se_y_proj_c2i = se_y + slope*((*leaf)->se_x - se_x);
					double nw_x_proj_c2i = se_x + (1./slope)*((*leaf)->nw_y - se_y);
					double se_x_proj_c2i = se_x + (1./slope)*((*leaf)->se_y - se_y);
					double nw_y_proj_i2c = (*leaf)->se_y + (*leaf)->slope*(nw_x - (*leaf)->se_x); // projection of the endpoints of incoming segment onto current segment (hence, i2c)
					double se_y_proj_i2c = (*leaf)->se_y + (*leaf)->slope*(se_x - (*leaf)->se_x);
					double nw_x_proj_i2c = (*leaf)->se_x + (1./(*leaf)->slope)*(nw_y - (*leaf)->se_y);
					double se_x_proj_i2c = (*leaf)->se_x + (1./(*leaf)->slope)*(se_y - (*leaf)->se_y);
					//printf("the projections: %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf\n",(*leaf)->nw_x,nw_y_proj_c2i,(*leaf)->se_x,se_y_proj_c2i,nw_x_proj_c2i,(*leaf)->nw_y,se_x_proj_c2i,(*leaf)->se_y,nw_x,nw_y_proj_i2c,se_x,se_y_proj_i2c,nw_x_proj_i2c,nw_y,se_x_proj_i2c,se_y);
					if( (*leaf)->nw_x <= x_intersect && x_intersect <= (*leaf)->se_x && nw_x <= x_intersect && x_intersect <= se_x ) // The segments intersect
					{
						if(printing) printf("location 8\n");
						if( fabs(nw_x - x_intersect) <= .0000001 || fabs(se_x - x_intersect) <= .0000001)
						{
/*							printf("the intersection is exactly an endpoint of the incoming line\n");*/
							goto THEY_DONT_INTERSECT;
						}
						if(use_hausdorff && epsilon > 0.)
						{
							if(nw_y >= nw_y_proj_i2c) // left side of incoming segment is above current segment
							{
								double temp_leaf_se_y = fmax((*leaf)->se_y,se_y);
								double temp_leaf_se_x = (*leaf)->se_x;
								if(temp_leaf_se_y != (*leaf)->se_y) temp_leaf_se_x = (*leaf)->se_x + 
									(1./(*leaf)->slope)*(temp_leaf_se_y - (*leaf)->se_y);
								double temp_se_x = fmin(se_x,(*leaf)->se_x);
								double temp_se_y = se_y;
								if(temp_se_x != se_x) temp_se_y = se_y + slope*(temp_se_x - se_x);
								
/*								printf("plot([%lf],[%lf],'-ro');\n",temp_leaf_se_x,temp_leaf_se_y);*/
/*								printf("plot([%lf],[%lf],'-ro');\n",temp_se_x,temp_se_y);*/
								
								double dist1 = 0., dist2 = 0.;
								double proj_x = 
									((1./(*leaf)->slope)*temp_se_x+(*leaf)->slope*(*leaf)->se_x+temp_se_y-(*leaf)->se_y)/
									((*leaf)->slope+(1./(*leaf)->slope)); // proj from se_x
								double proj_y = (*leaf)->nw_y + (*leaf)->slope*(proj_x - (*leaf)->nw_x);
								if(proj_y != proj_y) proj_y = temp_se_y + (-1./(*leaf)->slope)*(proj_x - temp_se_x);
				
								if(proj_x <= temp_leaf_se_x && proj_x >= (*leaf)->nw_x) dist1 = 
									sqrt( (temp_se_x - proj_x)*(temp_se_x - proj_x) + 
									(temp_se_y - proj_y)*(temp_se_y - proj_y) );
								else dist1 = sqrt( (temp_se_x - temp_leaf_se_x)*(temp_se_x - temp_leaf_se_x) + 
									(temp_se_y - temp_leaf_se_y)*(temp_se_y - temp_leaf_se_y) );
									
								proj_x = ((slope)*nw_x+(1./slope)*temp_leaf_se_x-nw_y+temp_leaf_se_y)/((1./slope)+(slope)); // proj from leaf se_x
								proj_y = temp_leaf_se_y + (-1./slope)*(proj_x - temp_leaf_se_x);
								if(proj_y != proj_y) proj_y = nw_y + (slope)*(proj_x - nw_x);
					
			/*					printf("plot([%lf],[%lf],'-ko');\n",x_intersect,y_intersect);*/
				
								if(proj_x <= temp_se_x && proj_x >= nw_x) dist2 = 
									sqrt( (temp_leaf_se_x - proj_x)*(temp_leaf_se_x - proj_x) +
									(temp_leaf_se_y - proj_y)*(temp_leaf_se_y - proj_y) );
								else dist2 = sqrt( (temp_se_x - temp_leaf_se_x)*(temp_se_x - temp_leaf_se_x) + 
									(temp_se_y - temp_leaf_se_y)*(temp_se_y - temp_leaf_se_y) );
									
								distance = fmax(dist1,dist2);
/*								printf("distance of %lf btwn segments\n", distance);*/
								if(distance < epsilon*max_range)
								{
									returned_val = mock_insert2(type,nw_x,nw_y,x_intersect,y_intersect,slope,&(*leaf)->left);
									if(returned_val == 1) return 1;
									mock_insert_return_val = max(mock_insert_return_val,returned_val);
									if(se_x > temp_se_x)
									{
										returned_val = mock_insert2(type,temp_se_x,temp_se_y,se_x,se_y,slope,
											&(*leaf)->right);
										if(returned_val == 1) return 1;
										mock_insert_return_val = max(mock_insert_return_val,returned_val);
									}
								}
								else return 1;
								
							}
							else // right side of incoming segment is above current segment
							{
								double temp_leaf_nw_x = fmax((*leaf)->nw_x,nw_x);
								double temp_leaf_nw_y = (*leaf)->nw_y;
								if(temp_leaf_nw_x != (*leaf)->nw_x) temp_leaf_nw_y = (*leaf)->se_y + 
									((*leaf)->slope)*(temp_leaf_nw_x - (*leaf)->se_x);
								double temp_nw_y = fmin(nw_y,(*leaf)->nw_y);
								double temp_nw_x = nw_x;
								if(temp_nw_y != nw_y) temp_nw_x = se_x + (1./slope)*(temp_nw_y - se_y);
								
/*								printf("plot([%lf],[%lf],'-ro');\n",temp_nw_x,temp_nw_y);*/
/*								printf("plot([%lf],[%lf],'-ro');\n",temp_leaf_nw_x,temp_leaf_nw_y);*/
								
								double dist1 = 0., dist2 = 0.;
								double proj_x = 
									((1./(*leaf)->slope)*temp_nw_x+(*leaf)->slope*(*leaf)->se_x+temp_nw_y-(*leaf)->se_y)/
									((*leaf)->slope+(1./(*leaf)->slope)); // proj from nw_x
								double proj_y = (*leaf)->nw_y + (*leaf)->slope*(proj_x - (*leaf)->nw_x);
								if(proj_y != proj_y) proj_y = temp_nw_y + (-1./(*leaf)->slope)*(proj_x - temp_nw_x);
				
								if(proj_x <= (*leaf)->se_x && proj_x >= temp_leaf_nw_x) dist1 = 
									sqrt( (temp_nw_x - proj_x)*(temp_nw_x - proj_x) + 
									(temp_nw_y - proj_y)*(temp_nw_y - proj_y) );
								else dist1 = sqrt( (temp_nw_x - temp_leaf_nw_x)*(temp_nw_x - temp_leaf_nw_x) + 
									(temp_nw_y - temp_leaf_nw_y)*(temp_nw_y - temp_leaf_nw_y) );
									
								proj_x = ((slope)*nw_x+(1./slope)*temp_leaf_nw_x-nw_y+temp_leaf_nw_y)/((1./slope)+(slope)); // proj from leaf nw_x
								proj_y = temp_leaf_nw_y + (-1./slope)*(proj_x - temp_leaf_nw_x);
								if(proj_y != proj_y) proj_y = nw_y + (slope)*(proj_x - nw_x);
					
			/*					printf("plot([%lf],[%lf],'-ko');\n",x_intersect,y_intersect);*/
				
								if(proj_x <= se_x && proj_x >= temp_nw_x) dist2 = 
									sqrt( (temp_leaf_nw_x - proj_x)*(temp_leaf_nw_x - proj_x) +
									(temp_leaf_nw_y - proj_y)*(temp_leaf_nw_y - proj_y) );
								else dist2 = sqrt( (temp_nw_x - temp_leaf_nw_x)*(temp_nw_x - temp_leaf_nw_x) + 
									(temp_nw_y - temp_leaf_nw_y)*(temp_nw_y - temp_leaf_nw_y) );
									
								distance = fmax(dist1,dist2);
/*								printf("distance of %lf btwn segments\n", distance);*/
								if(distance < epsilon*max_range)
								{
									returned_val = mock_insert2(type,x_intersect,y_intersect,se_x,se_y,slope,&(*leaf)->left);
									if(returned_val == 1) return 1;
									mock_insert_return_val = max(mock_insert_return_val,returned_val);
									if(nw_x < temp_nw_x)
									{
										returned_val = mock_insert2(type,nw_x,nw_y,temp_nw_x,temp_nw_y,slope,
											&(*leaf)->right);
										if(returned_val == 1) return 1;
										mock_insert_return_val = max(mock_insert_return_val,returned_val);
									}
								}
								else return 1;
							}
						}
						else return 1;
					}
					else // The segments don't intersect
					{
						THEY_DONT_INTERSECT:
						if( nw_x - (*leaf)->nw_x <= .0000001) // first point of incoming segment is left of current segment
						{
							if((*leaf)->nw_y - (double)nw_y_proj_c2i <= .00000001 ) // || (*leaf)->se_y <= (double)se_y_proj_c2i )
							{
								if(nw_y_proj_c2i == nw_y_proj_c2i)
								{
									if(printing) printf("location 9\n");
									//printf("inserting 27\n");
									returned_val = mock_insert2(type,nw_x,nw_y,
											(*leaf)->nw_x,(double)nw_y_proj_c2i,slope,&(*leaf)->left);
									if(returned_val == 1) return 1;
									mock_insert_return_val = max(mock_insert_return_val,returned_val);
								}
								if(se_x_proj_c2i == se_x_proj_c2i)
								{
									if(printing) printf("location 10\n");
									//printf("inserting 28\n");
									returned_val = mock_insert2(type,
											(double)se_x_proj_c2i,(*leaf)->se_y,se_x,se_y,slope,&(*leaf)->right);
									if(returned_val == 1) return 1;
									mock_insert_return_val = max(mock_insert_return_val,returned_val);
								}
							}
							else // right portion of segment dominates a portion of current segment
							{
								if(printing) printf("location 11\n");
								goto LOCATION12;
							}
						}
						else if( nw_y - (double)nw_y_proj_i2c <= .00000001 ) // first point of incoming segment lies below the current segment
						{
							LOCATION12:
							
							if(printing) printf("location 12\n");
							if(use_hausdorff && epsilon > 0.)
							{
								double dist1 = 0., dist2 = 0, dist3 = 0., dist4 = 0.;
								
								double temp_leaf_nw_x = fmax((*leaf)->nw_x,nw_x);
								double temp_leaf_nw_y = (*leaf)->nw_y;
								if(temp_leaf_nw_x != (*leaf)->nw_x) temp_leaf_nw_y = (*leaf)->se_y + 
									((*leaf)->slope)*(temp_leaf_nw_x - (*leaf)->se_x);
								double temp_nw_y = fmin(nw_y,(*leaf)->nw_y);
								double temp_nw_x = nw_x;
								if(temp_nw_y != nw_y) temp_nw_x = se_x + (1./slope)*(temp_nw_y - se_y);
								double temp_leaf_se_y = fmax((*leaf)->se_y,se_y);
								double temp_leaf_se_x = (*leaf)->se_x;
								if(temp_leaf_se_y != (*leaf)->se_y) temp_leaf_se_x = (*leaf)->se_x + 
									(1./(*leaf)->slope)*(temp_leaf_se_y - (*leaf)->se_y);
								double temp_se_x = fmin(se_x,(*leaf)->se_x);
								double temp_se_y = se_y;
								if(temp_se_x != se_x) temp_se_y = se_y + slope*(temp_se_x - se_x);
								
/*								printf("plot([%lf],[%lf],'-ro');\n",temp_nw_x,temp_nw_y);*/
/*								printf("plot([%lf],[%lf],'-ro');\n",temp_se_x,temp_se_y);*/
/*								printf("plot([%lf],[%lf],'-ro');\n",temp_leaf_nw_x,temp_leaf_nw_y);*/
/*								printf("plot([%lf],[%lf],'-ro');\n",temp_leaf_se_x,temp_leaf_se_y);*/
								
								double x_proj = 
									((1./(*leaf)->slope)*temp_se_x+(*leaf)->slope*temp_leaf_se_x+temp_se_y-temp_leaf_se_y)/
									((*leaf)->slope+(1./(*leaf)->slope)); // proj from se_x
								double y_proj = temp_leaf_nw_y + (*leaf)->slope*(x_proj - temp_leaf_nw_x);
								if(y_proj != y_proj) y_proj = temp_se_y + (-1./(*leaf)->slope)*(x_proj - temp_se_x);
				
								if(x_proj <= temp_leaf_se_x && x_proj >= temp_leaf_nw_x) dist1 = 
									sqrt( (temp_se_x - x_proj)*(temp_se_x - x_proj) + 
									(temp_se_y - y_proj)*(temp_se_y - y_proj) );
								else dist1 = fmin( sqrt( (temp_se_x - temp_leaf_nw_x)*(temp_se_x - temp_leaf_nw_x) + 
									(temp_se_y - temp_leaf_nw_y)*(temp_se_y - temp_leaf_nw_y) ),
									sqrt( (temp_se_x - temp_leaf_se_x)*(temp_se_x - temp_leaf_se_x) + 
									(temp_se_y - temp_leaf_se_y)*(temp_se_y - temp_leaf_se_y) ) );
						
								x_proj = ((1./(*leaf)->slope)*temp_nw_x+(*leaf)->slope*temp_leaf_se_x+temp_nw_y-temp_leaf_se_y)/
									((*leaf)->slope+(1./(*leaf)->slope)); // proj from nw_x
								y_proj = temp_leaf_nw_y + (*leaf)->slope*(x_proj - temp_leaf_nw_x);
								if(y_proj != y_proj) y_proj = temp_nw_y + (-1./(*leaf)->slope)*(x_proj - temp_nw_x);
				
								if(x_proj <= temp_leaf_se_x && x_proj >= temp_leaf_nw_x) dist2 = 
									sqrt( (temp_nw_x - x_proj)*(temp_nw_x - x_proj) + 
									(temp_nw_y - y_proj)*(temp_nw_y - y_proj) );
								else dist2 = fmin( sqrt( (temp_nw_x - temp_leaf_nw_x)*(temp_nw_x - temp_leaf_nw_x) +
									(temp_nw_y - temp_leaf_nw_y)*(temp_nw_y - temp_leaf_nw_y) ),
									sqrt( (temp_nw_x - temp_leaf_se_x)*(temp_nw_x - temp_leaf_se_x) + 
									(temp_nw_y - temp_leaf_se_y)*(temp_nw_y - temp_leaf_se_y) ) );
						
								x_proj = ((slope)*temp_nw_x+(1./slope)*temp_leaf_se_x-temp_nw_y+temp_leaf_se_y)/
									((1./slope)+(slope)); // proj from leaf se_x
								y_proj = temp_leaf_se_y + (-1./slope)*(x_proj - temp_leaf_se_x);
								if(y_proj != y_proj) y_proj = temp_nw_y + (slope)*(x_proj - temp_nw_x);
				
								if(x_proj <= temp_se_x && x_proj >= temp_nw_x) dist3 = 
									sqrt( (temp_leaf_se_x - x_proj)*(temp_leaf_se_x - x_proj) +
									(temp_leaf_se_y - y_proj)*(temp_leaf_se_y - y_proj) );
								else dist3 = fmin( sqrt( (temp_nw_x - temp_leaf_se_x)*(temp_nw_x - temp_leaf_se_x) + 
									(temp_nw_y - temp_leaf_se_y)*(temp_nw_y - temp_leaf_se_y) ),
									sqrt( (temp_se_x - temp_leaf_se_x)*(temp_se_x - temp_leaf_se_x) + 
									(temp_se_y - temp_leaf_se_y)*(temp_se_y - temp_leaf_se_y) ) );
						
								x_proj = ((slope)*temp_nw_x+(1./slope)*temp_leaf_nw_x-temp_nw_y+temp_leaf_nw_y)/
									((1./slope)+(slope)); // proj from leaf nw_x
								y_proj = temp_leaf_nw_y + (-1./slope)*(x_proj - temp_leaf_nw_x);
								if(y_proj != y_proj) y_proj = temp_nw_y + (slope)*(x_proj - temp_nw_x);
				
								if(x_proj <= temp_se_x && x_proj >= temp_nw_x) dist4 = 
									sqrt( (temp_leaf_nw_x - x_proj)*(temp_leaf_nw_x - x_proj) +
									(temp_leaf_nw_y - y_proj)*(temp_leaf_nw_y - y_proj) );
								else dist4 = fmin( sqrt( (temp_nw_x - temp_leaf_nw_x)*(temp_nw_x - temp_leaf_nw_x) + 
									(temp_nw_y - temp_leaf_nw_y)*(temp_nw_y - temp_leaf_nw_y) ),
									sqrt( (temp_se_x - temp_leaf_nw_x)*(temp_se_x - temp_leaf_nw_x) + 
									(temp_se_y - temp_leaf_nw_y)*(temp_se_y - temp_leaf_nw_y) ) );
						
								dist1 = fmax(dist1,dist2);
								dist3 = fmax(dist3,dist4);
								distance = fmax(dist1,dist3);
					
			/*					printf("the max distance: %lf\n",distance);*/
					
/*								printf("distance of %lf btwn segments\n", distance);*/
								if(distance < epsilon*max_range)
								{
									int yes = 0;
									if(se_x > temp_se_x)
									{
										yes = 1;
										returned_val = mock_insert2(type,temp_se_x,temp_se_y,se_x,se_y,slope,
											&(*leaf)->right);
										if(returned_val == 1) return 1;
										mock_insert_return_val = max(mock_insert_return_val,returned_val);
									}
									if(nw_x < temp_nw_x)
									{
										yes = 1;
										returned_val = mock_insert2(type,nw_x,nw_y,temp_nw_x,temp_nw_y,slope,
											&(*leaf)->right);
										if(returned_val == 1) return 1;
										mock_insert_return_val = max(mock_insert_return_val,returned_val);
									}
									if(!yes) return 1;
								}
								else return 1;
							}
							else return 1;
						}
						else // first point of incoming segment lies above or to the right of current segment
						{
							if(se_x_proj_c2i == se_x_proj_c2i)
							{
								if(printing) printf("location 13\n");
								//printf("inserting 33\n");
								returned_val = mock_insert2(type,(double)se_x_proj_c2i,
											(*leaf)->se_y,se_x,se_y,slope,&(*leaf)->right);
								if(returned_val == 1) return 1;
								mock_insert_return_val = max(mock_insert_return_val,returned_val);
							}
						}
					}
				}
			}
		}
	}
/*	printf("returning %d\n",mock_insert_return_val);*/
	return mock_insert_return_val;
}

/*int printing = 0;*/

int points_only = 0, its_been_only_points = 0;
node *empty_node2 = NULL;
int printing_ = 0, real_counter = 0;

void insert(int type, double nw_x, double nw_y, double se_x, double se_y, double slope, struct node **leaf, struct node **leaf2)
{
	start_insert_timer = clock();
	if( clock_gettime( CLOCK_REALTIME, &start1) == -1 ) {
      		perror( "clock gettime" );
      		exit( EXIT_FAILURE );
    		}
	insert2(type, nw_x, nw_y, se_x, se_y, slope, &(*leaf), &(*leaf2));
	finish_insert_timer = clock();
	if( clock_gettime( CLOCK_REALTIME, &stop1) == -1 ) {
      		perror( "clock gettime" );
      		exit( EXIT_FAILURE );
	    	}

	    	accum = ( stop1.tv_sec - start1.tv_sec )
		  	+ ( stop1.tv_nsec - start1.tv_nsec )
		    	/ BILLION;
	     	insert_time2 += accum;
        insert_time += (double)(finish_struct_timer - start_struct_timer) / CLOCKS_PER_SEC;
} 

void insert2(int type, double nw_x, double nw_y, double se_x, double se_y, double slope, struct node **leaf, struct node **leaf2) // This is the main function, takes an input point or segment and removes portions of the tree that are dominated and then inserts portions of the point/segment that are not dominated into the appropriate places within the tree.
{	
/*	printing = 0;*/
/*	printf("******************************************\n");*/
/*	print_preorder(tree,NULL);*/
/*	printf("******************************************\n");*/
/*	printf("------------------------------------------\n");*/
/*	print_inorder(tree);*/
/*	printf("------------------------------------------\n");*/

/*	if(*leaf == tree) printf("inserting at tree\n");*/
/*	else if(*leaf == tree2 && tree2 != NULL) */
/*	{*/
/*		printf("inserting at dual bd tree\n");	*/
/*	}*/
	
/*	if(*leaf2) printf("leaf2 exists\n");*/
/*	if(!(*leaf2)) printf("leaf2 doesn't exist\n");*/
/*	if(*leaf2){*/
/*	printf("leaf2: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",(*leaf2)->type,x_ideal-(*leaf2)->nw_x,y_ideal-(*leaf2)->nw_y,x_ideal-(*leaf2)->se_x,y_ideal-(*leaf2)->se_y);*/
/*	}*/
	
	if(printing)
	{
	printf("inserted data is: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",type,nw_x,nw_y,se_x,se_y);
	printf("or maybe it was: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",type,x_ideal-nw_x,y_ideal-nw_y,x_ideal-se_x,y_ideal-se_y);
	
/*	printf("slope: %lf\t calculated: %lf\n",slope,(nw_y - se_y)/(nw_x-se_x));*/
/*	node *temp_node3 = tree;*/
	if(*leaf){
	printf("comparing against: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",(*leaf)->type,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y);
	printf("or ...: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",(*leaf)->type,x_ideal-(*leaf)->nw_x,y_ideal-(*leaf)->nw_y,x_ideal-(*leaf)->se_x,y_ideal-(*leaf)->se_y);
	}
	}
/*	if (*leaf && (*leaf)->type == 2 && (*leaf)->nw_x < 5.377857 && (*leaf)->nw_x > 5.37785 )*/
/*	{*/
/*		do_nothing(1);*/
/*	}*/
/*	if(insert_counter2 > 20) exit(0);*/
	if((*leaf) == tree || insert_to_potential_branch_tree == 1)
	{
		insert_to_potential_branch_tree = 0;
		if(type == 1) // Ensures points are dealt with appropriately
		{
			se_x = nw_x;
			se_y = nw_y;
		}
		if(del_root == 0)
		{	
			nw_x =  x_ideal - nw_x;
			se_x =  x_ideal - se_x;
			nw_y =  y_ideal - nw_y;
			se_y =  y_ideal - se_y;
			insert_counter++;
			real_counter++;
		}
		else del_root = 0;
		if(type == 2) slope = (se_y-nw_y)/(se_x-nw_x);
	}
	//printf("after messing: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",type,nw_x,nw_y,se_x,se_y);
	if(nw_x > se_x || se_y > nw_y)
	{
		return;
	}
	else if(type == 2 && fabs((double) nw_x - (double) se_x) <= .00001)
	{
		type = 1;
		nw_y = se_y;
	}
	else if(type == 2 && fabs((double) nw_y - (double) se_y) <= .00001)
	{
		type = 1;
		se_x = nw_x;
	}
	if(*leaf && type == 1 && ((fabs(nw_x - (*leaf)->nw_x) < .000001 && fabs(nw_y - (*leaf)->nw_y) < .000001) || (fabs(nw_x - (*leaf)->se_x) < .000001 && fabs(nw_y - (*leaf)->se_y) < .000001)))
	{	//printf("point already in tree\n");
		return;
	}
	else if(*leaf && type == 2 && (fabs(nw_x - (*leaf)->nw_x) < .00001 && fabs(nw_y - (*leaf)->nw_y) < .00001 && fabs(se_x - (*leaf)->se_x) < .00001 && fabs(se_y - (*leaf)->se_y) < .00001))
	{	
/*		printf("segment already in tree\n");*/
		return;
	}
	if(*leaf && *leaf == tree)
	{
		if(insert_counter > rebalance_count) // && !show_progress)
		{
			rebalance_count = round(1.00733*insert_counter);	//use 1.00733 if rebalancing only using periodic full tree, use 8 if combining with traversed path rebalance
/*			printf("rebalancing tree\n");*/
			Rebalance(tree);
			another_counter++;
			prev_node = NULL;
/*			printf("before cleaning, size is %d\n",tree->subtree_size);*/
/*			start_time = clock();*/
/*			clean_it(tree,0);*/
/*			finish_time = clock();*/
/*			double time_to = (double)(finish_time - start_time) / CLOCKS_PER_SEC;*/
/*			printf("after cleaning, size is %d and the time to do it was %lf\n",tree->subtree_size,time_to);*/
		}
	}
//	if(*leaf && (*leaf)->subtree_size > 2)
//	{
//		Rebalance2(*leaf);
//	}
//	if(tree)						//These lines can be used for debugging
//	{
//		scan_proper_tree(tree);
//	}
    	if( *leaf == NULL) // Node where we are trying to insert does not exist. This means that a new *leaf must be created.
   	{	
/*   		printf("creating new node: %d \t %.12lf,%.12lf  %.12lf,%.12lf\n",type,nw_x,nw_y,se_x,se_y);*/
   	
   		if(type == 2) 
		{
			points_only = 0;
			its_been_only_points = 0;
		}
   	
      		*leaf = (struct node*) malloc( sizeof( struct node ) ); 
        	(*leaf)->type = type;
		(*leaf)->nw_x = nw_x;
		(*leaf)->se_x = se_x;
		(*leaf)->nw_y = nw_y;
		(*leaf)->se_y = se_y;
		(*leaf)->slope = slope;
		(*leaf)->subtree_size = 1;
		if(*leaf == tree || !(*leaf2) )
		{
/*			printf("setting parent to null\n");*/
/*			if(*leaf != tree)*/
/*			{	*/
/*				 printf("this is weird and should not be happening\n");*/
/*				 exit(0);*/
/*			}*/
			(*leaf)->parent = NULL;
		}
		else 
		{
/*			printf("setting parent to incoming node\n");*/
			(*leaf)->parent = *leaf2;
		}
        	(*leaf)->left = NULL;    
        	(*leaf)->right = NULL;
        	insert_counter++;
        	if(*leaf != tree && *leaf2) update_depth(*leaf2);
	}
	else // The node being compared to does exist
	{	
		if(type == 1)
		{
			if((*leaf)->type == 1)
			{
				if(nw_x <= (*leaf)->nw_x && nw_y <= (*leaf)->se_y) // The point is dominated
				{	node *temp_node = NULL;
					if((*leaf)->left) temp_node = (*leaf)->left;
					else if((*leaf)->right) temp_node = (*leaf)->right;
					if((*leaf)->subtree_size == 1)
					{
						(*leaf)->type = type;
						(*leaf)->nw_x = nw_x;
						(*leaf)->se_x = se_x;
						(*leaf)->nw_y = nw_y;
						(*leaf)->se_y = se_y;
						(*leaf)->slope = slope;
						insert_counter++;
						return;
					}
					else
					{
						if(*leaf == tree) del_root = 1;
/*						printf("(%d) deleting root\n",__LINE__);*/
						delete_node(*leaf);
						if(printing_) printf("inserting 1\n");
						insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf),&(*leaf));
					}
				}
				else if(nw_x < (*leaf)->nw_x && nw_y >= (*leaf)->se_y) // The point needs inserted left
				{	
					if(printing_) printf("inserting 2\n");
					insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,&(*leaf));
/*					if((*leaf)->left)*/
/*					{*/
/*						(*leaf)->left->parent = *leaf;*/
/*						if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
				else if(nw_x >= (*leaf)->nw_x && nw_y < (*leaf)->se_y) // The point needs inserted right
				{	
					if(printing_) printf("inserting 3\n");
					insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*					if((*leaf)->right)*/
/*					{*/
/*						(*leaf)->right->parent = *leaf;*/
/*						if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
			}	
			else
			{
				if(nw_x <= (*leaf)->nw_x && nw_y >= (*leaf)->nw_y) // The point needs inserted left
				{	
					if(printing_) printf("inserting 4\n");
					insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,&(*leaf));
/*					if((*leaf)->left)*/
/*					{*/
/*						(*leaf)->left->parent = *leaf;*/
/*						if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
				else if(nw_x >= (*leaf)->se_x && nw_y <= (*leaf)->se_y) // The point needs inserted right
				{	
					if(printing_) printf("inserting 5\n");
					insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*					if((*leaf)->right)*/
/*					{*/
/*						(*leaf)->right->parent = *leaf;*/
/*						if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
				else
				{
					double x_proj = (*leaf)->se_x + (1./(*leaf)->slope)*(nw_y - (*leaf)->se_y);
					double y_proj = (*leaf)->se_y + (*leaf)->slope*(nw_x - (*leaf)->se_x);
					//printf("incoming: %lf,%lf to %lf,%lf\t current: %lf,%lf to %lf,%lf\t projection:%lf,%lf\n",nw_x,nw_y,se_x,se_y,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y,x_proj,y_proj);
					if(nw_x < x_proj || nw_y < y_proj) // The point is under (the extension of) the segment
					{
						if( nw_x <= (*leaf)->nw_x && x_proj >= (*leaf)->se_x ) // The entire segment is dominated
						{	node *temp_node = NULL;
							if((*leaf)->left) temp_node = (*leaf)->left;
							else if((*leaf)->right) temp_node = (*leaf)->right;
							if((*leaf)->subtree_size == 1)
							{
								(*leaf)->type = type;
								(*leaf)->nw_x = nw_x;
								(*leaf)->se_x = se_x;
								(*leaf)->nw_y = nw_y;
								(*leaf)->se_y = se_y;
								(*leaf)->slope = slope;
								insert_counter++;
								return;
							}
							else
							{
								if(*leaf == tree) del_root = 1;
/*								printf("(%d) deleting root\n",__LINE__);*/
								delete_node(*leaf);
								if(printing_) printf("inserting 6\n");
								insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf),&(*leaf));
							}
						}
						else if( nw_x <= (*leaf)->nw_x ) // The left portion of the segment is dominated
						{
							(*leaf)->nw_x = x_proj;
							(*leaf)->nw_y = nw_y;
							if(printing_) printf("inserting 7\n");
							insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,&(*leaf));
/*							if((*leaf)->left)*/
/*							{*/
/*								(*leaf)->left->parent = *leaf;*/
/*								if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*							}*/
						}
						else if( x_proj >= (*leaf)->se_x ) // The right portion of the segment is dominated
						{
							(*leaf)->se_y = (*leaf)->se_y + ((*leaf)->se_y - (*leaf)->nw_y)/((*leaf)->se_x - (*leaf)->nw_x)*(nw_x - (*leaf)->se_x);
							(*leaf)->se_x = nw_x;
							if(printing_) printf("inserting 8\n");
							insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*							if((*leaf)->right)*/
/*							{*/
/*								(*leaf)->right->parent = *leaf;*/
/*								if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*							}*/
						}
						else // The center portion of the segment is dominated
						{
							double save1 = (*leaf)->se_x;
							double save2 = (*leaf)->se_y;
							(*leaf)->se_y = (*leaf)->se_y + (*leaf)->slope*(nw_x - (*leaf)->se_x);
							(*leaf)->se_x = nw_x;
							if(printing_) printf("inserting 8b\n");
							insert2((*leaf)->type,x_proj,nw_y,save1,save2,(*leaf)->slope,&(*leaf)->right,&(*leaf));
/*							if((*leaf)->right)*/
/*							{	*/
/*								(*leaf)->right->parent = *leaf;*/
/*								if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*							}*/
							if(printing_) printf("inserting 9\n");
							insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*							if((*leaf)->right)*/
/*							{*/
/*								(*leaf)->right->parent = *leaf;*/
/*								if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*							}*/
						}
					}
				}
			}
		}
		else // Input is a line segment
		{
			if((*leaf)->type == 1) // Incoming segment compared against current point
			{
				double y_proj = se_y + slope*((*leaf)->nw_x - se_x);
				double x_proj = se_x + (1./slope)*((*leaf)->nw_y - se_y);
				//printf("incoming: %lf,%lf to %lf,%lf\t current: %lf,%lf to %lf,%lf\n",nw_x,nw_y,se_x,se_y,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y);
				//printf("the projections are: %.12lf,%.12lf and %.12lf,%.12lf\n",x_proj,(*leaf)->nw_y,(*leaf)->nw_x,y_proj);
				//printf("slope: %lf\t calculated: %lf\n",slope,(se_y - nw_y)/(se_x - nw_x));
				if(se_x <= (*leaf)->nw_x && se_y >= (*leaf)->nw_y) // The segment needs inserted left
				{	if(printing_) printf("inserting 10\n");
					insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,&(*leaf));
/*					if((*leaf)->left)*/
/*					{*/
/*						(*leaf)->left->parent = *leaf;*/
/*						if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
				else if(nw_x >= (*leaf)->nw_x && nw_y <= (*leaf)->nw_y) // The segment needs inserted right
				{	
					if(printing_) printf("inserting 11\n");
					insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*					if((*leaf)->right)*/
/*					{*/
/*						(*leaf)->right->parent = *leaf;*/
/*						if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
				else if(x_proj <= ((*leaf)->nw_x + .0001) || y_proj <= ((*leaf)->nw_y + .0001)) // The point is dominated
				{	node *temp_node = NULL;
					if((*leaf)->left) temp_node = (*leaf)->left;
					else if((*leaf)->right) temp_node = (*leaf)->right;
					if((*leaf)->subtree_size == 1)
					{
						(*leaf)->type = type;
						(*leaf)->nw_x = nw_x;
						(*leaf)->se_x = se_x;
						(*leaf)->nw_y = nw_y;
						(*leaf)->se_y = se_y;
						(*leaf)->slope = slope;
						insert_counter++;
						return;
					}
					else
					{
						if(*leaf == tree) del_root = 1;
/*						printf("(%d) deleting root\n",__LINE__);*/
						delete_node(*leaf);
						if(printing_) printf("inserting 12\n");
						insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf),&(*leaf));
					}
				}
				else 
				{
					double save = (*leaf)->nw_y;
					if(printing_) printf("inserting 13\n");
					insert2(type,nw_x,nw_y,(*leaf)->nw_x,y_proj,slope,&(*leaf)->left,&(*leaf));
/*					if((*leaf)->left)*/
/*					{*/
/*						(*leaf)->left->parent = *leaf;*/
/*						if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
					if(printing_) printf("inserting 14\n");
					insert2(type,x_proj,save,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*					if((*leaf)->right)*/
/*					{*/
/*						(*leaf)->right->parent = *leaf;*/
/*						if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
			}	
			else // Incoming segment compared against current segment
			{	
				if(se_x <= (*leaf)->nw_x && se_y >= (*leaf)->nw_y) // The segment needs inserted left
				{	
					if(printing_) printf("inserting 15\n");
					insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,&(*leaf));
/*					if((*leaf)->left)*/
/*					{*/
/*						(*leaf)->left->parent = *leaf;*/
/*						if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
				else if(nw_x >= (*leaf)->se_x && nw_y <= (*leaf)->se_y) // The segment needs inserted right
				{	//printf("incoming segment: %lf,%lf to %lf,%lf \t current segment: %lf,%lf to %lf,%lf\n",nw_x,nw_y,se_x,se_y,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y);
					if(printing_) printf("inserting 16\n");
					insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*					if((*leaf)->right)*/
/*					{*/
/*						(*leaf)->right->parent = *leaf;*/
/*						if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
				else
				{	
/*					printf("incoming segment: %.12lf,%.12lf to %.12lf,%.12lf \t current segment: %.12lf,%.12lf to %.12lf,%.12lf\n",nw_x,nw_y,se_x,se_y,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y);*/
/*					double x_intersect = ((*leaf)->nw_y*se_x*(*leaf)->se_x-(*leaf)->nw_y*nw_x*(*leaf)->se_x-nw_y*se_x*(*leaf)->se_x+nw_y*se_x*(*leaf)->nw_x+se_y*nw_x*(*leaf)->se_x-se_y*nw_x*(*leaf)->nw_x-(*leaf)->se_y*se_x*(*leaf)->nw_x+(*leaf)->se_y*nw_x*(*leaf)->nw_x)/(-(*leaf)->se_y*se_x+(*leaf)->se_y*nw_x+(*leaf)->nw_y*se_x-(*leaf)->nw_y*nw_x+se_y*(*leaf)->se_x-se_y*(*leaf)->nw_x-nw_y*(*leaf)->se_x+nw_y*(*leaf)->nw_x); // Find the intersection of (the possible extensions of) the line segments being compared*/
					double x_intersect = (slope*se_x-(*leaf)->slope*(*leaf)->se_x-se_y+(*leaf)->se_y)/(-(*leaf)->slope+slope);
					double y_intersect = (*leaf)->nw_y + (*leaf)->slope*(x_intersect - (*leaf)->nw_x);
					if(y_intersect != y_intersect) y_intersect = nw_y + slope*(x_intersect - nw_x);
					if(printing) printf("the intersection is: %.12lf,%.12lf\n",x_intersect,y_intersect);
					if(printing) printf("plot(%lf,%lf,'go');\n",x_intersect,y_intersect);
					double nw_y_proj_c2i = se_y + slope*((*leaf)->nw_x - se_x); // projection of the endpoints of current segment onto incoming segment (hence, c2i)
					double se_y_proj_c2i = se_y + slope*((*leaf)->se_x - se_x);
					double nw_x_proj_c2i = se_x + (1./slope)*((*leaf)->nw_y - se_y);
					double se_x_proj_c2i = se_x + (1./slope)*((*leaf)->se_y - se_y);
					double nw_y_proj_i2c = (*leaf)->se_y + (*leaf)->slope*(nw_x - (*leaf)->se_x); // projection of the endpoints of incoming segment onto current segment (hence, i2c)
					double se_y_proj_i2c = (*leaf)->se_y + (*leaf)->slope*(se_x - (*leaf)->se_x);
					double nw_x_proj_i2c = (*leaf)->se_x + (1./(*leaf)->slope)*(nw_y - (*leaf)->se_y);
					double se_x_proj_i2c = (*leaf)->se_x + (1./(*leaf)->slope)*(se_y - (*leaf)->se_y);
					//printf("the projections: %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf\n",(*leaf)->nw_x,nw_y_proj_c2i,(*leaf)->se_x,se_y_proj_c2i,nw_x_proj_c2i,(*leaf)->nw_y,se_x_proj_c2i,(*leaf)->se_y,nw_x,nw_y_proj_i2c,se_x,se_y_proj_i2c,nw_x_proj_i2c,nw_y,se_x_proj_i2c,se_y);
					if( (*leaf)->nw_x - x_intersect <= .0000001 && x_intersect - (*leaf)->se_x <= .0000001 && nw_x - x_intersect <= .0000001 
						&& x_intersect - se_x <= .0000001 ) // The segments intersect
					{
						if( (double)nw_x_proj_c2i <= (*leaf)->nw_x || (double)nw_y_proj_c2i <= (*leaf)->nw_y ) // Left side of incoming segment is below current segment
						{	
							if( (*leaf)->nw_x < nw_x ) // Left most portion of current segment is non-dominated
							{
								double save1 = (*leaf)->se_y;
								double save2 = (*leaf)->se_x;
        							(*leaf)->se_x = nw_x;
								(*leaf)->se_y = (double)nw_y_proj_i2c;
								if(printing_) printf("inserting 17\n");
								insert2((*leaf)->type,x_intersect,y_intersect,save2,save1,(*leaf)->slope,&(*leaf)->right,&(*leaf));
/*								if((*leaf)->right)*/
/*								{*/
/*									(*leaf)->right->parent = *leaf;*/
/*									if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*								}*/
								if(printing_) printf("inserting 18\n");
								insert2(type,nw_x,nw_y,x_intersect,y_intersect,slope,&(*leaf)->right,&(*leaf));
/*								if((*leaf)->right)*/
/*								{*/
/*									(*leaf)->right->parent = *leaf;*/
/*									if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*								}*/
								if(se_x_proj_c2i == se_x_proj_c2i)
								{
									if(printing_) printf("inserting 19\n");
									insert2(type,(double)se_x_proj_c2i,save1,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*									if((*leaf)->right)*/
/*									{*/
/*										(*leaf)->right->parent = *leaf;*/
/*										if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*									}*/
								}
							}
							else
							{
								(*leaf)->nw_x = x_intersect;
								(*leaf)->nw_y = y_intersect;
								double save = (*leaf)->se_y;
								if(printing_) printf("inserting 20\n");
								insert2(type,nw_x,nw_y,x_intersect,y_intersect,slope,&(*leaf)->left,&(*leaf));
/*								if((*leaf)->left)*/
/*								{*/
/*									(*leaf)->left->parent = *leaf;*/
/*									if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*								}*/
								if(se_x_proj_c2i == se_x_proj_c2i)
								{
									if(printing_) printf("inserting 21\n");
									insert2(type,(double)se_x_proj_c2i,save,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*									if((*leaf)->right)*/
/*									{*/
/*										(*leaf)->right->parent = *leaf;*/
/*										if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*									}*/
								}
							}
						}
						else // Right side of incoming segment is below current segment
						{	//printf("right side of incoming segment is below current segment\n");
							if( (*leaf)->se_y < se_y ) // Right most portion of current segment is non-dominated
							{	//printf("Right most portion of current segment is non-dominated\n");
								double save1 = (*leaf)->se_x;
								double save2 = (*leaf)->se_y;
								(*leaf)->se_x = x_intersect;
								(*leaf)->se_y = y_intersect;
								double save = (*leaf)->nw_x;
								if(se_x_proj_i2c == se_x_proj_i2c)
								{
									if(printing_) printf("inserting 22\n");
									insert2((*leaf)->type,(double)se_x_proj_i2c,se_y,save1,save2,
										(*leaf)->slope,&(*leaf)->right,&(*leaf));
/*										if((*leaf)->right)*/
/*										{*/
/*											(*leaf)->right->parent = *leaf;*/
/*											if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*										}*/
								}
								if(printing_) printf("inserting 23\n");
								insert2(type,x_intersect,y_intersect,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*									if((*leaf)->right)*/
/*									{*/
/*										(*leaf)->right->parent = *leaf;*/
/*										if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*									}*/
								if(nw_y_proj_c2i == nw_y_proj_c2i)
								{
									if(printing_) printf("inserting 24\n");
									insert2(type,nw_x,nw_y,save,(double)nw_y_proj_c2i,slope,&(*leaf)->left,&(*leaf));
/*										if((*leaf)->left)*/
/*										{*/
/*											(*leaf)->left->parent = *leaf;*/
/*											if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*										}*/
								}
							}
							else
							{
								(*leaf)->se_x = x_intersect;
								(*leaf)->se_y = y_intersect;
								double save = (*leaf)->nw_x;
								if(printing_) printf("inserting 25\n");
								insert2(type,x_intersect,y_intersect,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*									if((*leaf)->right)*/
/*									{*/
/*										(*leaf)->right->parent = *leaf;*/
/*										if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*									}*/
								if(nw_y_proj_c2i == nw_y_proj_c2i)
								{
									if(printing_) printf("inserting 26\n");
									insert2(type,nw_x,nw_y,save,(double)nw_y_proj_c2i,slope,&(*leaf)->left,&(*leaf));
/*										if((*leaf)->left)*/
/*										{*/
/*											(*leaf)->left->parent = *leaf;*/
/*											if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*										}*/
								}
							}
						}
					}
					else // The segments don't intersect
					{
							if( nw_x <= (*leaf)->nw_x ) // first point of incoming segment is left of current segment
							{
								if((*leaf)->nw_y <= (double)nw_y_proj_c2i) // || (*leaf)->se_y <= (double)se_y_proj_c2i )
								{
									double save = (*leaf)->se_y;
									if(nw_y_proj_c2i == nw_y_proj_c2i)
									{
										if(printing_) printf("inserting 27\n");
										insert2(type,nw_x,nw_y,(*leaf)->nw_x,(double)nw_y_proj_c2i,slope,
											&(*leaf)->left,&(*leaf));
/*										if((*leaf)->left)*/
/*										{*/
/*											(*leaf)->left->parent = *leaf;*/
/*											if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*										}*/
									}
									if(se_x_proj_c2i == se_x_proj_c2i)
									{
										if(printing_) printf("inserting 28\n");
										insert2(type,(double)se_x_proj_c2i,save,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*										if((*leaf)->right)*/
/*										{*/
/*											(*leaf)->right->parent = *leaf;*/
/*											if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*										}*/
									}
								}
								else
								{
									if( se_y <= (*leaf)->se_y ) // incoming segment dominates current segment
									{	node *temp_node = NULL;
										if((*leaf)->left) temp_node = (*leaf)->left;
										else if((*leaf)->right) temp_node = (*leaf)->right;
										else temp_node = (*leaf)->parent;
										if((*leaf)->subtree_size == 1)
										{
											(*leaf)->type = type;
											(*leaf)->nw_x = nw_x;
											(*leaf)->se_x = se_x;
											(*leaf)->nw_y = nw_y;
											(*leaf)->se_y = se_y;
											(*leaf)->slope = slope;
											insert_counter++;
											return;
										}
										else
										{
											if(*leaf == tree) del_root = 1;
/*											printf("(%d) deleting root\n",__LINE__);*/
											delete_node(*leaf);
											if(printing_) printf("inserting 29\n");
											insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf),&(*leaf));
										}
									}
									else // upper portion of current segment dominated
									{
										(*leaf)->nw_x = (double)se_x_proj_i2c;
										(*leaf)->nw_y = se_y;
										if(printing_) printf("inserting 30\n");
										insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,&(*leaf));
/*										if((*leaf)->left)*/
/*										{*/
/*											(*leaf)->left->parent = *leaf;*/
/*											if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*										}*/
									}
								}
							}
							else if( nw_y <= (double)nw_y_proj_i2c ) // first point of incoming segment lies below the current segment
							{
								if( se_y <= (*leaf)->se_y ) // lower portion of current segment dominated
								{
									(*leaf)->se_x = nw_x;
									(*leaf)->se_y = (double)nw_y_proj_i2c;
									if(printing_) printf("inserting 31\n");
									insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*									if((*leaf)->right)*/
/*									{*/
/*										(*leaf)->right->parent = *leaf;*/
/*										if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*									}*/
								}
								else // center portion of current segment dominated
								{
									double save1 = (*leaf)->se_x;
									double save2 = (*leaf)->se_y;
        								(*leaf)->se_x = nw_x;
									(*leaf)->se_y = (double)nw_y_proj_i2c;
									if(se_x_proj_i2c == se_x_proj_i2c)
									{
										if(printing_) printf("inserting 31.5\n");
										insert2((*leaf)->type,(double)se_x_proj_i2c,se_y,save1,save2,(*leaf)->slope,
											&(*leaf)->right,&(*leaf));
/*										if((*leaf)->right)*/
/*										{*/
/*											(*leaf)->right->parent = *leaf;*/
/*											if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*										}*/
									}
									if(printing_) printf("inserting 32\n");
									insert2(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*									if((*leaf)->right)*/
/*									{*/
/*										(*leaf)->right->parent = *leaf;*/
/*										if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*									}*/
								}
							}
							else // first point of incoming segment lies above or to the right of current segment
							{
								if(se_x_proj_c2i == se_x_proj_c2i)
								{
									if(printing_) printf("inserting 33\n");
									insert2(type,(double)se_x_proj_c2i,(*leaf)->se_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*									if((*leaf)->right)*/
/*									{*/
/*										(*leaf)->right->parent = *leaf;*/
/*										if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*									}*/
								}
							}
						}
					}
				}
			}
		}
}
int cn = 0;

void insert_db(int type, double nw_x, double nw_y, double se_x, double se_y, double slope, struct node **leaf, struct node **leaf2) // This is the main function, takes an input point or segment and removes portions of the tree that are dominated and then inserts portions of the point/segment that are not dominated into the appropriate places within the tree.
{	
/*	cn++;*/
/*	if(cn > 12000) */
/*	{*/
/*		printf("------------------------------------------\n");*/
/*		print_inorder(tree2,2);*/
/*		printf("------------------------------------------\n");*/
/*		exit(0);*/
/*	}*/
/*	printing = 1;*/
/*	printf("******************************************\n");*/
/*	print_preorder(tree,NULL);*/
/*	printf("******************************************\n");*/
/*	printf("------------------------------------------\n");*/
/*	print_inorder(tree);*/
/*	printf("------------------------------------------\n");*/
	if(printing)
	{
		if(*leaf == tree) 
		{
/*			printf("inserted data is: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",type,nw_x,nw_y,se_x,se_y);*/
			printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",nw_x,se_x,nw_y,se_y);
		}
		else 
		{
/*			printf("or maybe it was: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",type,x_ideal-nw_x,y_ideal-nw_y,x_ideal-se_x,y_ideal-se_y);*/
			printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",x_ideal-nw_x,x_ideal-se_x,y_ideal-nw_y,y_ideal-se_y);
		}
	
/*	printf("slope: %lf\t calculated: %lf\n",slope,(nw_y - se_y)/(nw_x-se_x));*/
/*	node *temp_node3 = tree;*/
	if(*leaf){
/*	printf("comparing against: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",(*leaf)->type,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y);*/
/*	printf("or ...: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",(*leaf)->type,x_ideal-(*leaf)->nw_x,y_ideal-(*leaf)->nw_y,x_ideal-(*leaf)->se_x,y_ideal-(*leaf)->se_y);*/
		printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ideal-(*leaf)->nw_x,x_ideal-(*leaf)->se_x,y_ideal-(*leaf)->nw_y,y_ideal-(*leaf)->se_y);
	}
	}
/*	if (*leaf && (*leaf)->type == 2 && (*leaf)->nw_x < 5.377857 && (*leaf)->nw_x > 5.37785 )*/
/*	{*/
/*		do_nothing(1);*/
/*	}*/
/*	if(insert_counter2 > 20) exit(0);*/
	if((*leaf) == tree2 || insert_to_potential_branch_tree == 1)
	{
		insert_to_potential_branch_tree = 0;
		if(type == 1) // Ensures points are dealt with appropriately
		{
			se_x = nw_x;

			se_y = nw_y;
		}
		if(del_root == 0)
		{	
			nw_x =  x_ideal - nw_x;
			se_x =  x_ideal - se_x;
			nw_y =  y_ideal - nw_y;
			se_y =  y_ideal - se_y;
/*			insert_counter2++;*/
		}
		else del_root = 0;
		if(type == 2) slope = (se_y-nw_y)/(se_x-nw_x);
	}
	//printf("after messing: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",type,nw_x,nw_y,se_x,se_y);
	if(nw_x > se_x || se_y > nw_y)
	{
		return;
	}
	else if(type == 2 && fabs((double) nw_x - (double) se_x) <= .00001)
	{
		type = 1;
		nw_y = se_y;
	}
	else if(type == 2 && fabs((double) nw_y - (double) se_y) <= .00001)
	{
		type = 1;
		se_x = nw_x;
	}
	if(*leaf && type == 1 && ((fabs(nw_x - (*leaf)->nw_x) < .000001 && fabs(nw_y - (*leaf)->nw_y) < .000001) || (fabs(nw_x - (*leaf)->se_x) < .000001 && fabs(nw_y - (*leaf)->se_y) < .000001)))
	{	//printf("point already in tree\n");
		return;
	}
	else if(*leaf && type == 2 && (fabs(nw_x - (*leaf)->nw_x) < .00001 && fabs(nw_y - (*leaf)->nw_y) < .00001 && fabs(se_x - (*leaf)->se_x) < .00001 && fabs(se_y - (*leaf)->se_y) < .00001))
	{	
/*		printf("segment already in tree\n");*/
		return;
	}
    	if( *leaf == NULL) // Node where we are trying to insert does not exist. This means that a new *leaf must be created.
   	{	//printf("creating new node: %d \t %.12lf,%.12lf  %.12lf,%.12lf\n",type,nw_x,nw_y,se_x,se_y);
      		*leaf = (struct node*) malloc( sizeof( struct node ) ); 
        	(*leaf)->type = type;
		(*leaf)->nw_x = nw_x;
		(*leaf)->se_x = se_x;
		(*leaf)->nw_y = nw_y;
		(*leaf)->se_y = se_y;
		(*leaf)->slope = slope;
		(*leaf)->subtree_size = 1;
		if(*leaf == tree2 || !(*leaf2) )
		{
/*			printf("setting parent to null\n");*/
			(*leaf)->parent = NULL;
		}
		else 
		{
/*			printf("setting parent to incoming node\n");*/
			(*leaf)->parent = *leaf2;
		}
        	(*leaf)->left = NULL;    
        	(*leaf)->right = NULL;
/*        	insert_counter++;*/
        	if(*leaf != tree2) update_depth(*leaf2);
	}
	else // The node being compared to does exist
	{	
		if(type == 1)
		{
			if((*leaf)->type == 1)
			{
				if(nw_x <= (*leaf)->nw_x && nw_y <= (*leaf)->se_y) // The point is dominated
				{	node *temp_node = NULL;
					if((*leaf)->left) temp_node = (*leaf)->left;
					else if((*leaf)->right) temp_node = (*leaf)->right;
					if((*leaf)->subtree_size == 1)
					{
						(*leaf)->type = type;
						(*leaf)->nw_x = nw_x;
						(*leaf)->se_x = se_x;
						(*leaf)->nw_y = nw_y;
						(*leaf)->se_y = se_y;
						(*leaf)->slope = slope;
/*						insert_counter++;*/
						return;
					}
					else
					{
						if(*leaf == tree2) del_root = 1;
						delete_node(*leaf);
						if(printing) printf("inserting 1\n");
						insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf),&(*leaf));
					}
				}
				else if(nw_x < (*leaf)->nw_x && nw_y >= (*leaf)->se_y) // The point needs inserted left
				{	
					if(printing) printf("inserting 2\n");
					insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,&(*leaf));
/*					if((*leaf)->left)*/
/*					{*/
/*						(*leaf)->left->parent = *leaf;*/
/*						if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
				else if(nw_x >= (*leaf)->nw_x && nw_y < (*leaf)->se_y) // The point needs inserted right
				{	
					if(printing) printf("inserting 3\n");
					insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*					if((*leaf)->right)*/
/*					{*/
/*						(*leaf)->right->parent = *leaf;*/
/*						if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
			}	
			else
			{
				if(nw_x <= (*leaf)->nw_x && nw_y >= (*leaf)->nw_y) // The point needs inserted left
				{	
					if(printing) printf("inserting 4\n");
					insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,&(*leaf));
/*					if((*leaf)->left)*/
/*					{*/
/*						(*leaf)->left->parent = *leaf;*/
/*						if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
				else if(nw_x >= (*leaf)->se_x && nw_y <= (*leaf)->se_y) // The point needs inserted right
				{	
					if(printing) printf("inserting 5\n");
					insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*					if((*leaf)->right)*/
/*					{*/
/*						(*leaf)->right->parent = *leaf;*/
/*						if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
				else
				{
					double x_proj = (*leaf)->se_x + (1./(*leaf)->slope)*(nw_y - (*leaf)->se_y);
					double y_proj = (*leaf)->se_y + (*leaf)->slope*(nw_x - (*leaf)->se_x);
					//printf("incoming: %lf,%lf to %lf,%lf\t current: %lf,%lf to %lf,%lf\t projection:%lf,%lf\n",nw_x,nw_y,se_x,se_y,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y,x_proj,y_proj);
					if(nw_x < x_proj || nw_y < y_proj) // The point is under (the extension of) the segment
					{
						if( nw_x <= (*leaf)->nw_x && x_proj >= (*leaf)->se_x ) // The entire segment is dominated
						{	node *temp_node = NULL;
							if((*leaf)->left) temp_node = (*leaf)->left;
							else if((*leaf)->right) temp_node = (*leaf)->right;
							if((*leaf)->subtree_size == 1)
							{
								(*leaf)->type = type;
								(*leaf)->nw_x = nw_x;
								(*leaf)->se_x = se_x;
								(*leaf)->nw_y = nw_y;
								(*leaf)->se_y = se_y;
								(*leaf)->slope = slope;
/*								insert_counter++;*/
								return;
							}
							else
							{
								if(*leaf == tree2) del_root = 1;
								delete_node(*leaf);
								if(printing) printf("inserting 6\n");
								insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf),&(*leaf));
							}
						}
						else if( nw_x <= (*leaf)->nw_x ) // The left portion of the segment is dominated
						{
							(*leaf)->nw_x = x_proj;
							(*leaf)->nw_y = nw_y;
							if(printing) printf("inserting 7\n");
							insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,&(*leaf));
/*							if((*leaf)->left)*/
/*							{*/
/*								(*leaf)->left->parent = *leaf;*/
/*								if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*							}*/
						}
						else if( x_proj >= (*leaf)->se_x ) // The right portion of the segment is dominated
						{
							(*leaf)->se_y = (*leaf)->se_y + ((*leaf)->se_y - (*leaf)->nw_y)/((*leaf)->se_x - (*leaf)->nw_x)*(nw_x - (*leaf)->se_x);
							(*leaf)->se_x = nw_x;
							if(printing) printf("inserting 8\n");
							insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*							if((*leaf)->right)*/
/*							{*/
/*								(*leaf)->right->parent = *leaf;*/
/*								if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*							}*/
						}
						else // The center portion of the segment is dominated
						{
							double save1 = (*leaf)->se_x;
							double save2 = (*leaf)->se_y;
							(*leaf)->se_y = (*leaf)->se_y + (*leaf)->slope*(nw_x - (*leaf)->se_x);
							(*leaf)->se_x = nw_x;
							if(printing) printf("inserting 8b\n");
							insert_db((*leaf)->type,x_proj,nw_y,save1,save2,(*leaf)->slope,&(*leaf)->right,&(*leaf));
/*							if((*leaf)->right)*/
/*							{	*/
/*								(*leaf)->right->parent = *leaf;*/
/*								if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*							}*/
							if(printing) printf("inserting 9\n");
							insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*							if((*leaf)->right)*/
/*							{*/
/*								(*leaf)->right->parent = *leaf;*/
/*								if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*							}*/
						}
					}
				}
			}
		}
		else // Input is a line segment
		{
			if((*leaf)->type == 1) // Incoming segment compared against current point
			{
				double y_proj = se_y + slope*((*leaf)->nw_x - se_x);
				double x_proj = se_x + (1./slope)*((*leaf)->nw_y - se_y);
				//printf("incoming: %lf,%lf to %lf,%lf\t current: %lf,%lf to %lf,%lf\n",nw_x,nw_y,se_x,se_y,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y);
				//printf("the projections are: %.12lf,%.12lf and %.12lf,%.12lf\n",x_proj,(*leaf)->nw_y,(*leaf)->nw_x,y_proj);
				//printf("slope: %lf\t calculated: %lf\n",slope,(se_y - nw_y)/(se_x - nw_x));
				if(se_x <= (*leaf)->nw_x && se_y >= (*leaf)->nw_y) // The segment needs inserted left
				{	if(printing) printf("inserting 10\n");
					insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,&(*leaf));
/*					if((*leaf)->left)*/
/*					{*/
/*						(*leaf)->left->parent = *leaf;*/
/*						if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
				else if(nw_x >= (*leaf)->nw_x && nw_y <= (*leaf)->nw_y) // The segment needs inserted right
				{	
					if(printing) printf("inserting 11\n");
					insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*					if((*leaf)->right)*/
/*					{*/
/*						(*leaf)->right->parent = *leaf;*/
/*						if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
				else if(x_proj <= ((*leaf)->nw_x + .0001) || y_proj <= ((*leaf)->nw_y + .0001)) // The point is dominated
				{	node *temp_node = NULL;
					if((*leaf)->left) temp_node = (*leaf)->left;
					else if((*leaf)->right) temp_node = (*leaf)->right;
					if((*leaf)->subtree_size == 1)
					{
						(*leaf)->type = type;
						(*leaf)->nw_x = nw_x;
						(*leaf)->se_x = se_x;
						(*leaf)->nw_y = nw_y;
						(*leaf)->se_y = se_y;
						(*leaf)->slope = slope;
/*						insert_counter++;*/
						return;
					}
					else
					{
						if(*leaf == tree2) del_root = 1;
						delete_node(*leaf);
						if(printing) printf("inserting 12\n");
						insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf),&(*leaf));
					}
				}
				else 
				{
					double save = (*leaf)->nw_y;
					if(printing) printf("inserting 13\n");
					insert_db(type,nw_x,nw_y,(*leaf)->nw_x,y_proj,slope,&(*leaf)->left,&(*leaf));
/*					if((*leaf)->left)*/
/*					{*/
/*						(*leaf)->left->parent = *leaf;*/
/*						if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
					if(printing) printf("inserting 14\n");
					insert_db(type,x_proj,save,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*					if((*leaf)->right)*/
/*					{*/
/*						(*leaf)->right->parent = *leaf;*/
/*						if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
			}	
			else // Incoming segment compared against current segment
			{	
				if(se_x <= (*leaf)->nw_x && se_y >= (*leaf)->nw_y) // The segment needs inserted left
				{	
					if(printing) printf("inserting 15\n");
					insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,&(*leaf));
/*					if((*leaf)->left)*/
/*					{*/
/*						(*leaf)->left->parent = *leaf;*/
/*						if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
				else if(nw_x >= (*leaf)->se_x && nw_y <= (*leaf)->se_y) // The segment needs inserted right
				{	//printf("incoming segment: %lf,%lf to %lf,%lf \t current segment: %lf,%lf to %lf,%lf\n",nw_x,nw_y,se_x,se_y,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y);
					if(printing) printf("inserting 16\n");
					insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*					if((*leaf)->right)*/
/*					{*/
/*						(*leaf)->right->parent = *leaf;*/
/*						if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*					}*/
				}
				else
				{	
/*					printf("incoming segment: %.12lf,%.12lf to %.12lf,%.12lf \t current segment: %.12lf,%.12lf to %.12lf,%.12lf\n",nw_x,nw_y,se_x,se_y,(*leaf)->nw_x,(*leaf)->nw_y,(*leaf)->se_x,(*leaf)->se_y);*/
/*					double x_intersect = ((*leaf)->nw_y*se_x*(*leaf)->se_x-(*leaf)->nw_y*nw_x*(*leaf)->se_x-nw_y*se_x*(*leaf)->se_x+nw_y*se_x*(*leaf)->nw_x+se_y*nw_x*(*leaf)->se_x-se_y*nw_x*(*leaf)->nw_x-(*leaf)->se_y*se_x*(*leaf)->nw_x+(*leaf)->se_y*nw_x*(*leaf)->nw_x)/(-(*leaf)->se_y*se_x+(*leaf)->se_y*nw_x+(*leaf)->nw_y*se_x-(*leaf)->nw_y*nw_x+se_y*(*leaf)->se_x-se_y*(*leaf)->nw_x-nw_y*(*leaf)->se_x+nw_y*(*leaf)->nw_x); // Find the intersection of (the possible extensions of) the line segments being compared*/
					double x_intersect = (slope*se_x-(*leaf)->slope*(*leaf)->se_x-se_y+(*leaf)->se_y)/(-(*leaf)->slope+slope);
					double y_intersect = (*leaf)->nw_y + (*leaf)->slope*(x_intersect - (*leaf)->nw_x);
					if(y_intersect != y_intersect) y_intersect = nw_y + slope*(x_intersect - nw_x);
					if(printing) printf("the intersection is: %.12lf,%.12lf\n",x_intersect,y_intersect);
					if(printing) printf("plot(%lf,%lf,'go');\n",x_intersect,y_intersect);
					double nw_y_proj_c2i = se_y + slope*((*leaf)->nw_x - se_x); // projection of the endpoints of current segment onto incoming segment (hence, c2i)
					double se_y_proj_c2i = se_y + slope*((*leaf)->se_x - se_x);
					double nw_x_proj_c2i = se_x + (1./slope)*((*leaf)->nw_y - se_y);

					double se_x_proj_c2i = se_x + (1./slope)*((*leaf)->se_y - se_y);
					double nw_y_proj_i2c = (*leaf)->se_y + (*leaf)->slope*(nw_x - (*leaf)->se_x); // projection of the endpoints of incoming segment onto current segment (hence, i2c)
					double se_y_proj_i2c = (*leaf)->se_y + (*leaf)->slope*(se_x - (*leaf)->se_x);
					double nw_x_proj_i2c = (*leaf)->se_x + (1./(*leaf)->slope)*(nw_y - (*leaf)->se_y);
					double se_x_proj_i2c = (*leaf)->se_x + (1./(*leaf)->slope)*(se_y - (*leaf)->se_y);
					//printf("the projections: %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf,  %.12lf,%.12lf\n",(*leaf)->nw_x,nw_y_proj_c2i,(*leaf)->se_x,se_y_proj_c2i,nw_x_proj_c2i,(*leaf)->nw_y,se_x_proj_c2i,(*leaf)->se_y,nw_x,nw_y_proj_i2c,se_x,se_y_proj_i2c,nw_x_proj_i2c,nw_y,se_x_proj_i2c,se_y);
					if( (*leaf)->nw_x - x_intersect <= .0000001 && x_intersect - (*leaf)->se_x <= .0000001 && nw_x - x_intersect <= .0000001 
						&& x_intersect - se_x <= .0000001 ) // The segments intersect
					{
						if( (double)nw_x_proj_c2i <= (*leaf)->nw_x || (double)nw_y_proj_c2i <= (*leaf)->nw_y ) // Left side of incoming segment is below current segment
						{	
							if( (*leaf)->nw_x < nw_x ) // Left most portion of current segment is non-dominated
							{
								double save1 = (*leaf)->se_y;
								double save2 = (*leaf)->se_x;
        							(*leaf)->se_x = nw_x;
								(*leaf)->se_y = (double)nw_y_proj_i2c;
								if(printing) printf("inserting 17\n");
								insert_db((*leaf)->type,x_intersect,y_intersect,save2,save1,(*leaf)->slope,&(*leaf)->right,&(*leaf));
/*								if((*leaf)->right)*/
/*								{*/
/*									(*leaf)->right->parent = *leaf;*/
/*									if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*								}*/
								if(printing) printf("inserting 18\n");
								insert_db(type,nw_x,nw_y,x_intersect,y_intersect,slope,&(*leaf)->right,&(*leaf));
/*								if((*leaf)->right)*/
/*								{*/
/*									(*leaf)->right->parent = *leaf;*/
/*									if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*								}*/
								if(se_x_proj_c2i == se_x_proj_c2i)
								{
									if(printing) printf("inserting 19\n");
									insert_db(type,(double)se_x_proj_c2i,save1,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*									if((*leaf)->right)*/
/*									{*/
/*										(*leaf)->right->parent = *leaf;*/
/*										if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*									}*/
								}
							}
							else
							{
								(*leaf)->nw_x = x_intersect;
								(*leaf)->nw_y = y_intersect;
								double save = (*leaf)->se_y;
								if(printing) printf("inserting 20\n");
								insert_db(type,nw_x,nw_y,x_intersect,y_intersect,slope,&(*leaf)->left,&(*leaf));
/*								if((*leaf)->left)*/
/*								{*/
/*									(*leaf)->left->parent = *leaf;*/
/*									if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*								}*/
								if(se_x_proj_c2i == se_x_proj_c2i)
								{
									if(printing) printf("inserting 21\n");
									insert_db(type,(double)se_x_proj_c2i,save,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*									if((*leaf)->right)*/
/*									{*/
/*										(*leaf)->right->parent = *leaf;*/
/*										if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*									}*/
								}
							}
						}
						else // Right side of incoming segment is below current segment
						{	//printf("right side of incoming segment is below current segment\n");
							if( (*leaf)->se_y < se_y ) // Right most portion of current segment is non-dominated
							{	//printf("Right most portion of current segment is non-dominated\n");
								double save1 = (*leaf)->se_x;
								double save2 = (*leaf)->se_y;
								(*leaf)->se_x = x_intersect;
								(*leaf)->se_y = y_intersect;
								double save = (*leaf)->nw_x;
								if(se_x_proj_i2c == se_x_proj_i2c)
								{
									if(printing) printf("inserting 22\n");
									insert_db((*leaf)->type,(double)se_x_proj_i2c,se_y,save1,save2,
										(*leaf)->slope,&(*leaf)->right,&(*leaf));
/*										if((*leaf)->right)*/
/*										{*/
/*											(*leaf)->right->parent = *leaf;*/
/*											if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*										}*/
								}
								if(printing) printf("inserting 23\n");
								insert_db(type,x_intersect,y_intersect,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*									if((*leaf)->right)*/
/*									{*/
/*										(*leaf)->right->parent = *leaf;*/
/*										if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*									}*/
								if(nw_y_proj_c2i == nw_y_proj_c2i)
								{
									if(printing) printf("inserting 24\n");
									insert_db(type,nw_x,nw_y,save,(double)nw_y_proj_c2i,slope,&(*leaf)->left,&(*leaf));
/*										if((*leaf)->left)*/
/*										{*/
/*											(*leaf)->left->parent = *leaf;*/
/*											if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*										}*/
								}
							}
							else
							{
								(*leaf)->se_x = x_intersect;
								(*leaf)->se_y = y_intersect;
								double save = (*leaf)->nw_x;
								if(printing) printf("inserting 25\n");
								insert_db(type,x_intersect,y_intersect,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*									if((*leaf)->right)*/
/*									{*/
/*										(*leaf)->right->parent = *leaf;*/
/*										if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*									}*/
								if(nw_y_proj_c2i == nw_y_proj_c2i)
								{
									if(printing) printf("inserting 26\n");
									insert_db(type,nw_x,nw_y,save,(double)nw_y_proj_c2i,slope,&(*leaf)->left,&(*leaf));
/*										if((*leaf)->left)*/
/*										{*/
/*											(*leaf)->left->parent = *leaf;*/
/*											if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*										}*/
								}
							}
						}
					}
					else // The segments don't intersect
					{
							if( nw_x <= (*leaf)->nw_x ) // first point of incoming segment is left of current segment
							{
								if((*leaf)->nw_y <= (double)nw_y_proj_c2i) // || (*leaf)->se_y <= (double)se_y_proj_c2i )
								{
									double save = (*leaf)->se_y;
									if(nw_y_proj_c2i == nw_y_proj_c2i)
									{
										if(printing) printf("inserting 27\n");
										insert_db(type,nw_x,nw_y,(*leaf)->nw_x,(double)nw_y_proj_c2i,slope,
											&(*leaf)->left,&(*leaf));
/*										if((*leaf)->left)*/
/*										{*/
/*											(*leaf)->left->parent = *leaf;*/
/*											if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*										}*/
									}
									if(se_x_proj_c2i == se_x_proj_c2i)
									{
										if(printing) printf("inserting 28\n");
										insert_db(type,(double)se_x_proj_c2i,save,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*										if((*leaf)->right)*/
/*										{*/
/*											(*leaf)->right->parent = *leaf;*/
/*											if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*										}*/
									}
								}
								else
								{
									if( se_y <= (*leaf)->se_y ) // incoming segment dominates current segment
									{	node *temp_node = NULL;
										if((*leaf)->left) temp_node = (*leaf)->left;
										else if((*leaf)->right) temp_node = (*leaf)->right;
										else temp_node = (*leaf)->parent;
										if((*leaf)->subtree_size == 1)
										{
											(*leaf)->type = type;
											(*leaf)->nw_x = nw_x;
											(*leaf)->se_x = se_x;
											(*leaf)->nw_y = nw_y;
											(*leaf)->se_y = se_y;
											(*leaf)->slope = slope;
/*											insert_counter++;*/
											return;
										}
										else
										{
											if(*leaf == tree2) del_root = 1;
											delete_node(*leaf);
											if(printing) printf("inserting 29\n");
											insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf),&(*leaf));
										}
									}
									else // upper portion of current segment dominated
									{
										(*leaf)->nw_x = (double)se_x_proj_i2c;
										(*leaf)->nw_y = se_y;
										if(printing) printf("inserting 30\n");
										insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,&(*leaf));
/*										if((*leaf)->left)*/
/*										{*/
/*											(*leaf)->left->parent = *leaf;*/
/*											if((*leaf)->left->subtree_size == 1) update_depth(*leaf);*/
/*										}*/

									}
								}
							}
							else if( nw_y <= (double)nw_y_proj_i2c ) // first point of incoming segment lies below the current segment
							{
								if( se_y <= (*leaf)->se_y ) // lower portion of current segment dominated
								{
									(*leaf)->se_x = nw_x;
									(*leaf)->se_y = (double)nw_y_proj_i2c;
									if(printing) printf("inserting 31\n");
									insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*									if((*leaf)->right)*/
/*									{*/
/*										(*leaf)->right->parent = *leaf;*/
/*										if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*									}*/
								}
								else // center portion of current segment dominated
								{
									double save1 = (*leaf)->se_x;
									double save2 = (*leaf)->se_y;
        								(*leaf)->se_x = nw_x;
									(*leaf)->se_y = (double)nw_y_proj_i2c;
									if(se_x_proj_i2c == se_x_proj_i2c)
									{
										if(printing) printf("inserting 31.5\n");
										insert_db((*leaf)->type,(double)se_x_proj_i2c,se_y,save1,save2,(*leaf)->slope,
											&(*leaf)->right,&(*leaf));
/*										if((*leaf)->right)*/
/*										{*/
/*											(*leaf)->right->parent = *leaf;*/
/*											if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*										}*/
									}
									if(printing) printf("inserting 32\n");
									insert_db(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*									if((*leaf)->right)*/
/*									{*/
/*										(*leaf)->right->parent = *leaf;*/
/*										if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*									}*/
								}
							}
							else // first point of incoming segment lies above or to the right of current segment
							{
								if(se_x_proj_c2i == se_x_proj_c2i)
								{
									if(printing) printf("inserting 33\n");
									insert_db(type,(double)se_x_proj_c2i,(*leaf)->se_y,se_x,se_y,slope,&(*leaf)->right,&(*leaf));
/*									if((*leaf)->right)*/
/*									{*/
/*										(*leaf)->right->parent = *leaf;*/
/*										if((*leaf)->right->subtree_size == 1) update_depth(*leaf);*/
/*									}*/
								}
							}
						}
					}
				}
			}
		}
}

int done_val = -1;

double get_max_proximal_hd_dist(int type, double nw_x, double nw_y, double se_x, double se_y, double slope, struct node **leaf,	double best_val) 
{	
	double new_val = 0.;
/*	if(best_val == 0.)*/
/*	{*/
/*		printf("...\n");*/
/*		printf("...\n");*/
/*		printf("...\n");*/
/*	}*/
	if( *leaf == NULL) 
   	{   	
/*   		printf("returning %lf (%d)\n",best_val,__LINE__);*/
   		return best_val;
	}
	
	if((*leaf) && ((*leaf) == tree || insert_to_potential_branch_tree == 1))
	{
		insert_to_potential_branch_tree = 0;
		if(type == 1) // Ensures points are dealt with appropriately
		{
			se_x = nw_x;
			se_y = nw_y;
		}
		if(del_root == 0)
		{	
			nw_x =  x_ideal - nw_x;
			se_x =  x_ideal - se_x;
			nw_y =  y_ideal - nw_y;
			se_y =  y_ideal - se_y;
		}
		else del_root = 0;
		if(type == 2) slope = (se_y-nw_y)/(se_x-nw_x);
	}
	//printf("after messing: %d\t %.12lf,%.12lf - %.12lf,%.12lf\n",type,nw_x,nw_y,se_x,se_y);
	if(nw_x > se_x || se_y > nw_y)
	{
/*		printf("the segment is backwards\n");*/

		if((*leaf)->left && (*leaf)->left->subtree_size > done_val) 
		{
			new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,best_val);
			(*leaf)->left->subtree_size = done_val;
		}
		best_val = fmax(new_val,best_val);
		if((*leaf)->right && (*leaf)->right->subtree_size > done_val) 
		{
			new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,best_val);
			(*leaf)->right->subtree_size = done_val;
		}
		best_val = fmax(new_val,best_val);

/*		printf("returning %lf (%d)\n",0.,__LINE__);*/
		return best_val;
	}
	else if(type == 2 && fabs((double) nw_x - (double) se_x) <= .00001)
	{
		type = 1;
		nw_y = se_y;
	}
	else if(type == 2 && fabs((double) nw_y - (double) se_y) <= .00001)
	{
		type = 1;
		se_x = nw_x;
	}
	if(*leaf && type == 1 && ((fabs(nw_x - (*leaf)->nw_x) < .000001 && fabs(nw_y - (*leaf)->nw_y) < .000001) || (fabs(nw_x - (*leaf)->se_x) < .000001 && fabs(nw_y - (*leaf)->se_y) < .000001)))
	{	
/*		printf("point already in tree\n");*/

/*		if((*leaf)->left) new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,best_val);*/
/*		best_val = fmax(new_val,best_val);*/
/*		if((*leaf)->right) new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,best_val);*/
/*		best_val = fmax(new_val,best_val);		*/
		
/*		printf("returning %lf (%d)\n",best_val,__LINE__);*/
		return best_val;
	}
	else if(*leaf && type == 2 && (fabs(nw_x - (*leaf)->nw_x) < .00001 && fabs(nw_y - (*leaf)->nw_y) < .00001 && fabs(se_x - (*leaf)->se_x) < .00001 && fabs(se_y - (*leaf)->se_y) < .00001))
	{	
/*		printf("segment already in tree\n");*/
/*		if((*leaf)->left) new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,best_val);*/
/*		best_val = fmax(new_val,best_val);*/
/*		if((*leaf)->right) new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,best_val);*/
/*		best_val = fmax(new_val,best_val);*/

/*		printf("returning %lf (%d)\n",best_val,__LINE__);*/
		return best_val;
	}
	
	double distance = 0.;
    	if( *leaf == NULL) 
   	{   	
/*   		printf("returning %lf (%d)\n",best_val,__LINE__);*/
   		return best_val;
	}
	else // The node being compared to does exist
	{	
/*		printf("\tplot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-(*leaf)->nw_x,x_ideal-(*leaf)->se_x,y_ideal-(*leaf)->nw_y,y_ideal-(*leaf)->se_y);*/
/*		printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",x_ideal-nw_x,x_ideal-se_x,y_ideal-nw_y,y_ideal-se_y);*/
		if(type == 1) // The inserted solution is a point.
		{
			if((*leaf)->type == 1) // The current solution is a point.
			{
				distance = sqrt( (nw_x - (*leaf)->nw_x)*(nw_x - (*leaf)->nw_x) + (nw_y - (*leaf)->nw_y)*(nw_y - (*leaf)->nw_y) );
				
				if(nw_x - (*leaf)->nw_x <= .00000001 && nw_y - (*leaf)->se_y <= .000000001 ) // The current point is dominated
				{	
					distance = sqrt( (nw_x - (*leaf)->nw_x)*(nw_x - (*leaf)->nw_x) + (nw_y - (*leaf)->nw_y)*(nw_y - (*leaf)->nw_y) );
					best_val = fmax(best_val,distance);
					
					if((*leaf)->left && (*leaf)->left->subtree_size > done_val) 
					{
						new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,best_val);
						(*leaf)->left->subtree_size = done_val;
					}
					best_val = fmax(new_val,best_val);
					if((*leaf)->right && (*leaf)->right->subtree_size > done_val) 
					{
						new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,best_val);
						(*leaf)->right->subtree_size = done_val;
					}
					best_val = fmax(new_val,best_val);
					
/*					printf("returning %lf (%d)\n",best_val,__LINE__);*/
					return best_val; 
				}
				else if(nw_x < (*leaf)->nw_x && nw_y >= (*leaf)->se_y) // The point needs inserted left
				{	
/*					printf("going left\n");*/
					new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,best_val);
					best_val = fmax(new_val,best_val);
				}
				else if(nw_x >= (*leaf)->nw_x && nw_y < (*leaf)->se_y) // The point needs inserted right
				{	
/*					printf("going right\n");*/
					new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,best_val);
					best_val = fmax(new_val,best_val);
				}
			}	
			else // The current solution is a segment
			{
				if(nw_x <= (*leaf)->nw_x && nw_y >= (*leaf)->nw_y) // The point needs inserted left
				{	
/*					printf("going left\n");*/
					new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,best_val);
					best_val = fmax(new_val,best_val);
				}
				else if(nw_x >= (*leaf)->se_x && nw_y <= (*leaf)->se_y) // The point needs inserted right
				{	
/*					printf("going right\n");*/
					new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,best_val);
					best_val = fmax(best_val,new_val);
				}
				else
				{
					double x_proj = (*leaf)->se_x + (1./(*leaf)->slope)*(nw_y - (*leaf)->se_y);
					double y_proj = (*leaf)->se_y + (*leaf)->slope*(nw_x - (*leaf)->se_x);
					if(nw_x < x_proj - 0.00001 || nw_y < y_proj - 0.00001) // The point is under (the extension of) the segment
					{
						double temp_nw_x = fmax(nw_x,(*leaf)->nw_x);
						double temp_nw_y = (*leaf)->nw_y;
						if(temp_nw_x != (*leaf)->nw_x) temp_nw_y = (*leaf)->se_y + (*leaf)->slope*(temp_nw_x - (*leaf)->se_x);
						double temp_se_y = fmax(se_y,(*leaf)->se_y);
						double temp_se_x = (*leaf)->se_x;
						if(temp_se_y != (*leaf)->se_y) temp_se_x = (*leaf)->se_x + (1./(*leaf)->slope)*(temp_se_y - (*leaf)->se_y);
						
						distance = fmax( sqrt( (nw_x - temp_nw_x)*(nw_x - temp_nw_x) + (nw_y - temp_nw_y)*(nw_y - temp_nw_y) ),
							sqrt( (nw_x - temp_se_x)*(nw_x - temp_se_x) + (nw_y - temp_se_y)*(nw_y - temp_se_y) ) );
						
						best_val = fmax(best_val,distance);
						
						if((*leaf)->left && (*leaf)->left->subtree_size > done_val) 
						{
							new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,best_val);
							(*leaf)->left->subtree_size = done_val;
						}
						best_val = fmax(new_val,best_val);
						if((*leaf)->right && (*leaf)->right->subtree_size > done_val) 
						{
							new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,best_val);
							(*leaf)->right->subtree_size = done_val;
						}
						best_val = fmax(new_val,best_val);
						
/*						printf("returning %lf (%d)\n",best_val,__LINE__);*/
						return best_val;
					}
				}
			}
		}
		else // Input is a line segment
		{
			int returned_val = 0;
			if((*leaf)->type == 1) // Incoming segment compared against current point
			{
/*				printf("%d\n",__LINE__);*/
				double y_proj = se_y + slope*((*leaf)->nw_x - se_x);
				double x_proj = se_x + (1./slope)*((*leaf)->nw_y - se_y);
				
				if(se_x - (*leaf)->nw_x <= .0000001 && se_y - (*leaf)->nw_y >= -.00000001) // The segment needs inserted left
				{	//printf("inserting 10\n");
					if(printing) printf("location 1\n");
					new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,best_val);
					best_val = fmax(new_val,best_val);
				}
				else if(nw_x - (*leaf)->nw_x >= -.00000001 && nw_y - (*leaf)->nw_y <= .000000001 ) // The segment needs inserted right
				{	//printf("inserting 11\n");
					if(printing) printf("location 2\n");
					new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,best_val);
					best_val = fmax(new_val,best_val);
				}
				else if(x_proj <= ((*leaf)->nw_x + .0001) || y_proj <= ((*leaf)->nw_y + .0001)) // The point is dominated
				{	
					double temp_nw_y = fmin(nw_y,(*leaf)->nw_y);
					double temp_nw_x = nw_x;
					if(temp_nw_y != nw_y) temp_nw_x = se_x + (1./slope)*(temp_nw_y - se_y);
					double temp_se_x = fmin(se_x,(*leaf)->se_x);
					double temp_se_y = se_y;
					if(temp_se_x != se_x) temp_se_y = se_y + slope*(temp_se_x - se_x);
					
					distance = fmax( sqrt( (temp_nw_x - (*leaf)->nw_x)*(temp_nw_x - (*leaf)->nw_x) + 
						(temp_nw_y - (*leaf)->nw_y)*(temp_nw_y - (*leaf)->nw_y) ),
						sqrt( (temp_se_x - (*leaf)->nw_x)*(temp_se_x - (*leaf)->nw_x) + 
						(temp_se_y - (*leaf)->nw_y)*(temp_se_y - (*leaf)->nw_y) ) );
						
					best_val = fmax(best_val,distance);
						
			 		if(temp_nw_x > nw_x)
			 		{
				 		new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,temp_nw_x,temp_nw_y,slope,&(*leaf)->left,best_val);
						best_val = fmax(new_val,best_val);
					}
					if(temp_se_x < se_x)
					{
						new_val = get_max_proximal_hd_dist(type,temp_se_x,temp_se_y,se_x,se_y,slope,&(*leaf)->right,best_val);
						best_val = fmax(new_val,best_val);
					}
					
					if((*leaf)->left && (*leaf)->left->subtree_size > done_val) 
					{
						new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,best_val);
						(*leaf)->left->subtree_size = done_val;
					}
					best_val = fmax(new_val,best_val);
					if((*leaf)->right && (*leaf)->right->subtree_size > done_val) 
					{
						new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,best_val);
						(*leaf)->right->subtree_size = done_val;
					}
					best_val = fmax(new_val,best_val);
					
/*					printf("returning %lf (%d)\n",best_val,__LINE__);*/
					return best_val;
				}
				else 
				{
/*					printf("location 4\n");*/
					if(nw_x != (*leaf)->nw_x)
					{
						new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,(*leaf)->nw_x,y_proj,slope,&(*leaf)->left,best_val);
						best_val = fmax(new_val,best_val);
					}
					if((*leaf)->nw_y != se_y)
					{
						new_val = get_max_proximal_hd_dist(type,x_proj,(*leaf)->nw_y,se_x,se_y,slope,&(*leaf)->right,best_val);
						best_val = fmax(new_val,best_val);
					}
				}
			}	
			else // Incoming segment compared against current segment
			{	
				if(se_x <= (*leaf)->nw_x + .00001 && se_y >= (*leaf)->nw_y - .00001) // The segment needs inserted left
				{	
					if(printing) printf("location 5\n");
					//printf("inserting 15\n");
					new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,best_val);
					best_val = fmax(new_val,best_val);
				}
				else if(nw_x - (*leaf)->se_x >= -.00000001 && nw_y - (*leaf)->se_y <= .00000001) // The segment needs inserted right
				{	
					new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,best_val);
					best_val = fmax(new_val,best_val);
				}
				else
				{	
					double x_intersect = (slope*se_x-(*leaf)->slope*(*leaf)->se_x-se_y+(*leaf)->se_y)/(-(*leaf)->slope+slope);
					double y_intersect = (*leaf)->nw_y + (*leaf)->slope*(x_intersect - (*leaf)->nw_x);
					if(y_intersect != y_intersect) y_intersect = nw_y + slope*(x_intersect - nw_x);
					
					double nw_y_proj_c2i = se_y + slope*((*leaf)->nw_x - se_x); 
					double se_y_proj_c2i = se_y + slope*((*leaf)->se_x - se_x);
					double nw_x_proj_c2i = se_x + (1./slope)*((*leaf)->nw_y - se_y);
					double se_x_proj_c2i = se_x + (1./slope)*((*leaf)->se_y - se_y);
					double nw_y_proj_i2c = (*leaf)->se_y + (*leaf)->slope*(nw_x - (*leaf)->se_x); 
					double se_y_proj_i2c = (*leaf)->se_y + (*leaf)->slope*(se_x - (*leaf)->se_x);
					double nw_x_proj_i2c = (*leaf)->se_x + (1./(*leaf)->slope)*(nw_y - (*leaf)->se_y);
					double se_x_proj_i2c = (*leaf)->se_x + (1./(*leaf)->slope)*(se_y - (*leaf)->se_y);
					
					if( (*leaf)->nw_x <= x_intersect && x_intersect <= (*leaf)->se_x && nw_x <= x_intersect && x_intersect <= se_x ) 
					{
						if(printing) printf("location 8\n");
						if( fabs(nw_x - x_intersect) <= .0000001 || fabs(se_x - x_intersect) <= .0000001)

						{
/*							printf("the intersection is exactly an endpoint of the incoming line\n");*/
							goto THEY_DONT_INTERSECT2;
						}
						
						if(nw_y >= nw_y_proj_i2c) // left side of incoming segment is above current segment
						{
							double temp_leaf_se_y = fmax((*leaf)->se_y,se_y);
							double temp_leaf_se_x = (*leaf)->se_x;
							if(temp_leaf_se_y != (*leaf)->se_y) temp_leaf_se_x = (*leaf)->se_x + 
								(1./(*leaf)->slope)*(temp_leaf_se_y - (*leaf)->se_y);
							double temp_se_x = fmin(se_x,(*leaf)->se_x);
							double temp_se_y = se_y;
							if(temp_se_x != se_x) temp_se_y = se_y + slope*(temp_se_x - se_x);
							
							double dist1 = 0., dist2 = 0.;
							double proj_x = 
								((1./(*leaf)->slope)*temp_se_x+(*leaf)->slope*(*leaf)->se_x+temp_se_y-(*leaf)->se_y)/
								((*leaf)->slope+(1./(*leaf)->slope)); // proj from se_x
							double proj_y = (*leaf)->nw_y + (*leaf)->slope*(proj_x - (*leaf)->nw_x);
							if(proj_y != proj_y) proj_y = temp_se_y + (-1./(*leaf)->slope)*(proj_x - temp_se_x);
			
							if(proj_x <= temp_leaf_se_x && proj_x >= (*leaf)->nw_x) dist1 = 
								sqrt( (temp_se_x - proj_x)*(temp_se_x - proj_x) + 
								(temp_se_y - proj_y)*(temp_se_y - proj_y) );
							else dist1 = sqrt( (temp_se_x - temp_leaf_se_x)*(temp_se_x - temp_leaf_se_x) + 
								(temp_se_y - temp_leaf_se_y)*(temp_se_y - temp_leaf_se_y) );
								
							proj_x = ((slope)*nw_x+(1./slope)*temp_leaf_se_x-nw_y+temp_leaf_se_y)/((1./slope)+(slope)); 
							proj_y = temp_leaf_se_y + (-1./slope)*(proj_x - temp_leaf_se_x);
							if(proj_y != proj_y) proj_y = nw_y + (slope)*(proj_x - nw_x);
			
							if(proj_x <= temp_se_x && proj_x >= nw_x) dist2 = 
								sqrt( (temp_leaf_se_x - proj_x)*(temp_leaf_se_x - proj_x) +
								(temp_leaf_se_y - proj_y)*(temp_leaf_se_y - proj_y) );
							else dist2 = sqrt( (temp_se_x - temp_leaf_se_x)*(temp_se_x - temp_leaf_se_x) + 
								(temp_se_y - temp_leaf_se_y)*(temp_se_y - temp_leaf_se_y) );
								
							distance = fmax(dist1,dist2);
							best_val = fmax(best_val,distance);
							
							new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,x_intersect,y_intersect,slope,&(*leaf)->left,best_val);
							best_val = fmax(new_val,best_val);
							if(se_x > temp_se_x)
							{
								new_val = get_max_proximal_hd_dist(type,temp_se_x,temp_se_y,se_x,se_y,slope,
									&(*leaf)->right,best_val);
								best_val = fmax(new_val,best_val);
							}
							
							if((*leaf)->left && (*leaf)->left->subtree_size > done_val) 
							{
								new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,best_val);
								(*leaf)->left->subtree_size = done_val;
							}
							best_val = fmax(new_val,best_val);
							if((*leaf)->right && (*leaf)->right->subtree_size > done_val) 
							{
								new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,best_val);
								(*leaf)->right->subtree_size = done_val;
							}
							
/*							printf("returning %lf (%d)\n",best_val,__LINE__);*/
							return best_val;
						}
						else // right side of incoming segment is above current segment
						{
							double temp_leaf_nw_x = fmax((*leaf)->nw_x,nw_x);
							double temp_leaf_nw_y = (*leaf)->nw_y;
							if(temp_leaf_nw_x != (*leaf)->nw_x) temp_leaf_nw_y = (*leaf)->se_y + 
								((*leaf)->slope)*(temp_leaf_nw_x - (*leaf)->se_x);
							double temp_nw_y = fmin(nw_y,(*leaf)->nw_y);
							double temp_nw_x = nw_x;
							if(temp_nw_y != nw_y) temp_nw_x = se_x + (1./slope)*(temp_nw_y - se_y);
							
							double dist1 = 0., dist2 = 0.;
							double proj_x = 
								((1./(*leaf)->slope)*temp_nw_x+(*leaf)->slope*(*leaf)->se_x+temp_nw_y-(*leaf)->se_y)/
								((*leaf)->slope+(1./(*leaf)->slope)); // proj from nw_x
							double proj_y = (*leaf)->nw_y + (*leaf)->slope*(proj_x - (*leaf)->nw_x);
							if(proj_y != proj_y) proj_y = temp_nw_y + (-1./(*leaf)->slope)*(proj_x - temp_nw_x);
			
							if(proj_x <= (*leaf)->se_x && proj_x >= temp_leaf_nw_x) dist1 = 
								sqrt( (temp_nw_x - proj_x)*(temp_nw_x - proj_x) + 
								(temp_nw_y - proj_y)*(temp_nw_y - proj_y) );
							else dist1 = sqrt( (temp_nw_x - temp_leaf_nw_x)*(temp_nw_x - temp_leaf_nw_x) + 
								(temp_nw_y - temp_leaf_nw_y)*(temp_nw_y - temp_leaf_nw_y) );
								
							proj_x = ((slope)*nw_x+(1./slope)*temp_leaf_nw_x-nw_y+temp_leaf_nw_y)/((1./slope)+(slope)); 
							proj_y = temp_leaf_nw_y + (-1./slope)*(proj_x - temp_leaf_nw_x);
							if(proj_y != proj_y) proj_y = nw_y + (slope)*(proj_x - nw_x);
			
							if(proj_x <= se_x && proj_x >= temp_nw_x) dist2 = 
								sqrt( (temp_leaf_nw_x - proj_x)*(temp_leaf_nw_x - proj_x) +
								(temp_leaf_nw_y - proj_y)*(temp_leaf_nw_y - proj_y) );
							else dist2 = sqrt( (temp_nw_x - temp_leaf_nw_x)*(temp_nw_x - temp_leaf_nw_x) + 
								(temp_nw_y - temp_leaf_nw_y)*(temp_nw_y - temp_leaf_nw_y) );
								
							distance = fmax(dist1,dist2);
							
							best_val = fmax(distance,best_val);
							
							new_val = get_max_proximal_hd_dist(type,x_intersect,y_intersect,se_x,se_y,slope,&(*leaf)->left,best_val);
							best_val = fmax(new_val,best_val);
							if(nw_x < temp_nw_x)
							{
								new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,temp_nw_x,temp_nw_y,slope,
									&(*leaf)->right,best_val);
								best_val = fmax(new_val,best_val);
							}
							
							if((*leaf)->left && (*leaf)->left->subtree_size > done_val) 
							{
								new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,best_val);
								(*leaf)->left->subtree_size = done_val;
							}
							best_val = fmax(new_val,best_val);
							if((*leaf)->right && (*leaf)->right->subtree_size > done_val) 
							{
								new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,best_val);
								(*leaf)->right->subtree_size = done_val;
							}
							
/*							printf("returning %lf (%d)\n",best_val,__LINE__);*/
							return best_val;
						}
					}
					else // The segments don't intersect
					{
						THEY_DONT_INTERSECT2:
						if( nw_x - (*leaf)->nw_x <= .0000001) // first point of incoming segment is left of current segment
						{
							if((*leaf)->nw_y - (double)nw_y_proj_c2i <= .00000001 ) // || (*leaf)->se_y <= (double)se_y_proj_c2i )
							{
								if(nw_y_proj_c2i == nw_y_proj_c2i)
								{
									if(printing) printf("location 9\n");
									//printf("inserting 27\n");
									new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,
											(*leaf)->nw_x,(double)nw_y_proj_c2i,slope,&(*leaf)->left,best_val);
									best_val = fmax(new_val,best_val);
								}
								if(se_x_proj_c2i == se_x_proj_c2i)
								{
									if(printing) printf("location 10\n");
									//printf("inserting 28\n");
									new_val = get_max_proximal_hd_dist(type,(double)se_x_proj_c2i,
											(*leaf)->se_y,se_x,se_y,slope,&(*leaf)->right,best_val);
									best_val = fmax(new_val,best_val);
								}
							}
							else // right portion of segment dominates a portion of current segment
							{
								if(printing) printf("location 11\n");
								goto LOCATION12b;
							}
						}
						else if( nw_y - (double)nw_y_proj_i2c <= .00000001 ) 
						{
							LOCATION12b:
							;
							
/*							printf("leaf: %lf\n",(*leaf)->nw_x);*/
/*							if((*leaf)->left) printf("leaf_left: %lf\n",(*leaf)->left->nw_x);*/
/*							if((*leaf)->right) printf("leaf_right: %lf\n",(*leaf)->right->nw_x);*/
							
							double dist1 = 0., dist2 = 0, dist3 = 0., dist4 = 0.;
								
							double temp_leaf_nw_x = fmax((*leaf)->nw_x,nw_x);
							double temp_leaf_nw_y = (*leaf)->nw_y;
							if(temp_leaf_nw_x != (*leaf)->nw_x) temp_leaf_nw_y = (*leaf)->se_y + 
								((*leaf)->slope)*(temp_leaf_nw_x - (*leaf)->se_x);
							double temp_nw_y = fmin(nw_y,(*leaf)->nw_y);
							double temp_nw_x = nw_x;
							if(temp_nw_y != nw_y) temp_nw_x = se_x + (1./slope)*(temp_nw_y - se_y);
							double temp_leaf_se_y = fmax((*leaf)->se_y,se_y);
							double temp_leaf_se_x = (*leaf)->se_x;
							if(temp_leaf_se_y != (*leaf)->se_y) temp_leaf_se_x = (*leaf)->se_x + 
								(1./(*leaf)->slope)*(temp_leaf_se_y - (*leaf)->se_y);
							double temp_se_x = fmin(se_x,(*leaf)->se_x);
							double temp_se_y = se_y;
							if(temp_se_x != se_x) temp_se_y = se_y + slope*(temp_se_x - se_x);
							
							double x_proj = 
								((1./(*leaf)->slope)*temp_se_x+(*leaf)->slope*temp_leaf_se_x+temp_se_y-temp_leaf_se_y)/
								((*leaf)->slope+(1./(*leaf)->slope)); // proj from se_x
							double y_proj = temp_leaf_nw_y + (*leaf)->slope*(x_proj - temp_leaf_nw_x);
							if(y_proj != y_proj) y_proj = temp_se_y + (-1./(*leaf)->slope)*(x_proj - temp_se_x);
			
							if(x_proj <= temp_leaf_se_x && x_proj >= temp_leaf_nw_x) dist1 = 

								sqrt( (temp_se_x - x_proj)*(temp_se_x - x_proj) + 
								(temp_se_y - y_proj)*(temp_se_y - y_proj) );
							else dist1 = fmin( sqrt( (temp_se_x - temp_leaf_nw_x)*(temp_se_x - temp_leaf_nw_x) + 
								(temp_se_y - temp_leaf_nw_y)*(temp_se_y - temp_leaf_nw_y) ),
								sqrt( (temp_se_x - temp_leaf_se_x)*(temp_se_x - temp_leaf_se_x) + 
								(temp_se_y - temp_leaf_se_y)*(temp_se_y - temp_leaf_se_y) ) );
					
							x_proj = ((1./(*leaf)->slope)*temp_nw_x+(*leaf)->slope*temp_leaf_se_x+temp_nw_y-temp_leaf_se_y)/
								((*leaf)->slope+(1./(*leaf)->slope)); // proj from nw_x
							y_proj = temp_leaf_nw_y + (*leaf)->slope*(x_proj - temp_leaf_nw_x);
							if(y_proj != y_proj) y_proj = temp_nw_y + (-1./(*leaf)->slope)*(x_proj - temp_nw_x);
			
							if(x_proj <= temp_leaf_se_x && x_proj >= temp_leaf_nw_x) dist2 = 
								sqrt( (temp_nw_x - x_proj)*(temp_nw_x - x_proj) + 
								(temp_nw_y - y_proj)*(temp_nw_y - y_proj) );
							else dist2 = fmin( sqrt( (temp_nw_x - temp_leaf_nw_x)*(temp_nw_x - temp_leaf_nw_x) +
								(temp_nw_y - temp_leaf_nw_y)*(temp_nw_y - temp_leaf_nw_y) ),
								sqrt( (temp_nw_x - temp_leaf_se_x)*(temp_nw_x - temp_leaf_se_x) + 
								(temp_nw_y - temp_leaf_se_y)*(temp_nw_y - temp_leaf_se_y) ) );
					
							x_proj = ((slope)*temp_nw_x+(1./slope)*temp_leaf_se_x-temp_nw_y+temp_leaf_se_y)/
								((1./slope)+(slope)); // proj from leaf se_x
							y_proj = temp_leaf_se_y + (-1./slope)*(x_proj - temp_leaf_se_x);
							if(y_proj != y_proj) y_proj = temp_nw_y + (slope)*(x_proj - temp_nw_x);
			
							if(x_proj <= temp_se_x && x_proj >= temp_nw_x) dist3 = 
								sqrt( (temp_leaf_se_x - x_proj)*(temp_leaf_se_x - x_proj) +
								(temp_leaf_se_y - y_proj)*(temp_leaf_se_y - y_proj) );
							else dist3 = fmin( sqrt( (temp_nw_x - temp_leaf_se_x)*(temp_nw_x - temp_leaf_se_x) + 
								(temp_nw_y - temp_leaf_se_y)*(temp_nw_y - temp_leaf_se_y) ),
								sqrt( (temp_se_x - temp_leaf_se_x)*(temp_se_x - temp_leaf_se_x) + 
								(temp_se_y - temp_leaf_se_y)*(temp_se_y - temp_leaf_se_y) ) );
					
							x_proj = ((slope)*temp_nw_x+(1./slope)*temp_leaf_nw_x-temp_nw_y+temp_leaf_nw_y)/
								((1./slope)+(slope)); // proj from leaf nw_x
							y_proj = temp_leaf_nw_y + (-1./slope)*(x_proj - temp_leaf_nw_x);
							if(y_proj != y_proj) y_proj = temp_nw_y + (slope)*(x_proj - temp_nw_x);
			
							if(x_proj <= temp_se_x && x_proj >= temp_nw_x) dist4 = 
								sqrt( (temp_leaf_nw_x - x_proj)*(temp_leaf_nw_x - x_proj) +
								(temp_leaf_nw_y - y_proj)*(temp_leaf_nw_y - y_proj) );
							else dist4 = fmin( sqrt( (temp_nw_x - temp_leaf_nw_x)*(temp_nw_x - temp_leaf_nw_x) + 
								(temp_nw_y - temp_leaf_nw_y)*(temp_nw_y - temp_leaf_nw_y) ),
								sqrt( (temp_se_x - temp_leaf_nw_x)*(temp_se_x - temp_leaf_nw_x) + 
								(temp_se_y - temp_leaf_nw_y)*(temp_se_y - temp_leaf_nw_y) ) );
					
							dist1 = fmax(dist1,dist2);
							dist3 = fmax(dist3,dist4);
							distance = fmax(dist1,dist3);
				
							best_val = fmax(distance,best_val);
				
							if(se_x > temp_se_x)
							{
								new_val = get_max_proximal_hd_dist(type,temp_se_x,temp_se_y,se_x,se_y,slope,
									&(*leaf)->right,best_val);
								best_val = fmax(new_val,best_val);
							}
							if(nw_x < temp_nw_x)
							{
								new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,temp_nw_x,temp_nw_y,slope,
									&(*leaf)->right,best_val);
								best_val = fmax(new_val,best_val);
							}
							
							if((*leaf)->left && (*leaf)->left->subtree_size > done_val) 
							{
								new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,best_val);
								(*leaf)->left->subtree_size = done_val;
							}
							best_val = fmax(new_val,best_val);
							if((*leaf)->right && (*leaf)->right->subtree_size > done_val) 
							{
								new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,best_val);
								(*leaf)->right->subtree_size = done_val;
							}
							
/*							printf("returning %lf (%d)\n",best_val,__LINE__);*/
							return best_val;
						}
						else // first point of incoming segment lies above or to the right of current segment
						{
							if(se_x_proj_c2i == se_x_proj_c2i)
							{
								if(printing) printf("location 13\n");
								//printf("inserting 33\n");
								new_val = get_max_proximal_hd_dist(type,(double)se_x_proj_c2i,
											(*leaf)->se_y,se_x,se_y,slope,&(*leaf)->right,best_val);
								best_val = fmax(new_val,best_val);
							}
						}
					}
				}
			}
		}
	}
	if((*leaf)->left && (*leaf)->left->subtree_size > done_val) 
	{
		new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->left,best_val);
		(*leaf)->left->subtree_size = done_val;
	}
	best_val = fmax(new_val,best_val);
	if((*leaf)->right && (*leaf)->right->subtree_size > done_val) 
	{
		new_val = get_max_proximal_hd_dist(type,nw_x,nw_y,se_x,se_y,slope,&(*leaf)->right,best_val);
		(*leaf)->right->subtree_size = done_val;
	}
	best_val = fmax(new_val,best_val);
	
/*	printf("returning %lf (%d)\n",best_val,__LINE__);*/
	return best_val;
}

int get_tree_depth(node *leaf)
{
	if(leaf->left && leaf->right) return max(get_tree_depth(leaf->left) + 1,get_tree_depth(leaf->right) + 1);
	else if(leaf->left) return get_tree_depth(leaf->left) + 1;
	else if(leaf->right) return get_tree_depth(leaf->right) + 1;
	else return 0;
}

int get_num_nodes(node *leaf)
{
	if(leaf->left && leaf->right) return get_num_nodes(leaf->left)+get_num_nodes(leaf->right)+1;
	else if(leaf->left) return get_num_nodes(leaf->left) + 1;
	else if(leaf->right) return get_num_nodes(leaf->right) + 1;
	else return 1;
}

int get_num_inserts()
{
	return real_counter;
}

int get_num_mock_inserts()
{
	return mock_insert_counter;
}

node *temp_node4;
struct nadir *theta;
int cnt = 0;
int ck = 0;

void populate_theta(node *cur_node)
{
	if(cur_node)
	{
		if(cur_node->right) populate_theta(cur_node->right);
		if(cur_node->type == 2 || ck != 1)
		{	
			//theta[cnt] = cur_node;
			//ck = 0;
			if(cur_node->type == 1)
			{
				theta[cnt].type = 1;
				theta[cnt].nw_x = x_ideal - temp_node4->nw_x;
				theta[cnt].nw_y = y_ideal - cur_node->se_y;
				theta[cnt].se_x = theta[cnt].nw_x;
				theta[cnt].se_y = theta[cnt].nw_y;
				cnt++;
				ck = 0;
			}
			else
			{
				if(ck != 1 && fabs(cur_node->se_x - temp_node4->nw_x) > 1/10000. && fabs(cur_node->se_y - temp_node4->nw_y) > 1/10000.)
				{
					theta[cnt].type = 1;
					theta[cnt].nw_x = x_ideal - temp_node4->nw_x;
					theta[cnt].nw_y = y_ideal - cur_node->se_y;
					theta[cnt].se_x = theta[cnt].nw_x;
					theta[cnt].se_y = theta[cnt].nw_y;
					cnt++;
					ck = 0;
				}
				if(fabs(cur_node->se_x - cur_node->nw_x) <= .00001 || fabs(cur_node->se_y - cur_node->nw_y) <= .00001) theta[cnt].type = 1;
				else theta[cnt].type = 2;
				theta[cnt].nw_x = x_ideal - cur_node->se_x;
				theta[cnt].nw_y = y_ideal - cur_node->se_y;
				theta[cnt].se_x = x_ideal - cur_node->nw_x;
				theta[cnt].se_y = y_ideal - cur_node->nw_y;
				cnt++;
				ck = 0;
			}
			//printf("cur node: %d \t (%lf,%lf),(%lf,%lf)\n",cur_node->type,x_ideal-cur_node->nw_x,y_ideal-cur_node->nw_y,x_ideal-cur_node->se_x,y_ideal-cur_node->se_y);
			//printf("associated theta entry: %d \t (%lf,%lf),(%lf,%lf)\n",theta[cnt-1].type,theta[cnt-1].nw_x,theta[cnt-1].nw_y,theta[cnt-1].se_x,theta[cnt].se_y);
			//cnt++;
		}
		else ck = 0;
		temp_node4 = cur_node;
		if(cur_node->left) populate_theta(cur_node->left);
	}
}

//struct nadir *build_theta(int size)
void build_theta(int size)
{
	//node *temp_node4 = (struct node*) malloc( sizeof( struct node ) );
	//temp_node4 = NULL;
	ck = 1;
	cnt = 1;
	//free(th);
	theta = (struct nadir*) realloc(theta, (2*size + 2) * sizeof(struct nadir));
	populate_theta(tree);
	//free(temp_node4);
	//return th;
}

//void build_theta2(int size)
//{
//	//node *temp_node4 = (struct node*) malloc( sizeof( struct node ) );
//	//temp_node4 = NULL;
//	ck = 1;
//	cnt = 0;
//	//free(th);
//	theta = (struct nadir*) realloc(theta, (2*size + 2) * sizeof(struct nadir));
//	populate_theta(tree);
//	//free(temp_node4);
//	//return th;
//}

int compare_ranks(const void* px, const void* py)
{
	int x = *(const int*)px;
   	int y = *(const int*)py;
	if(theta[x].rank_val > theta[y].rank_val) return 1;
	else return -1;
}

/*int main()*/
/*{*/
/*	srand48(time(NULL));*/
/*  	srand(time(NULL));*/
/*	//srand(72); //started w 65 and 6 then 71 and 12*/
/*	//srand48(13);*/

/*	#define MAX_ITER 10/2*/
/*	#define N_TIMES 2*/

/*	double k = 1.;//MAX_ITER;*/
/*	int i = 0;*/
/*	int counter = 0;*/
/*	int max_counter = 100000;*/
/*	//file = fopen("pareto_output.txt", "w");	//Uncomment these lines to output MATLAB syntax*/
/*	//drawing = fopen("drawn_tree.txt", "w");*/

/*	//FILE *fp;*/
/*   */
/*    	//fp = fopen("full_output.txt", "w");*/

/*	while (counter < max_counter)*/
/*    	{*/
/*		n = round(1+(rand() % 100)/20.);*/
/*		m = n;*/
/*		if (counter == 401)*/
/*		{*/
/*			int yourmom = 16;*/
/*		}*/
/*		double points[2*n];*/
/*		a =(floor (drand48() * 1.0e9)) / 1.0e8;*/
/*		b = (10-a)*(10-a);*/
/*		c = a + (5 - k);*/
/*		points[0] = c;*/
/*		points[1] = b;*/
/*		n--;*/
/*		int j = 2;*/
/*		while( n >= 0 )*/
/*		{*/
/*			a = a + ((floor (drand48() * 1.0e9)) / 1.0e8)/(10/(10-a));*/
/*			b = (10-a)*(10-a);*/
/*			c = a + (5 - k);*/
/*			if ( c - points[j-2] < 0.1 )*/
/*			{*/
/*				m = j/2-1;*/
/*				break;*/
/*			}*/
/*			points[j] = c;*/
/*			points[j+1] = b;*/
/*			j = j+2;*/
/*			n--;*/
/*		}*/
/*		int r;*/
/*		for (r = 0; r < 2*m; r = r+2)*/
/*		{*/
/*			if (m == 1)*/
/*			{*/
/*				insert2(1,points[0],points[1],0,0,0.,&tree);*/
/*				counter++;*/
/*				//printf("The generated point is: %lf,%lf \n",points[0],points[1]);*/
/*				//fprintf(fp,"plot(%lf,%lf,'o');\n",points[0],points[1]);*/
/*			}*/
/*			else*/
/*			{*/
/*				insert2(2,fmax(points[r],points[r+2]),fmin(points[r+1],points[r+3]),fmin(points[r],points[r+2]),fmax(points[r+1],points[r+3]),(points[r+3]-points[r+1])/(points[r+2]-points[r]),&tree);*/
/*				counter++;*/
/*				//printf("The %d generated segment is: %lf,%lf to %lf,%lf \n",r/2+1,points[r],points[r+1],points[r+2],points[r+3]);*/
/*				//fprintf(fp,"plot([%lf,%lf],[%lf,%lf],'-o');\n",points[r],points[r+2],points[r+1],points[r+3]);*/
/*				if (m == 2) break;*/
/*			}*/
/*		}*/
/*		k = k + 2./(max_counter);*/
/*		//printf("k: %.16lf\n",k);*/
/*		i++;*/
/*    }*/
/*    	//printf("tree depth: %d\n",tree->depth);*/
/*    	//printf("rebalance count: %d\n",another_counter);*/
/*    	//printf("single rebalance count: %d\n",counter2);*/
/*    	//printf("powerful rebalance count: %d\n",counter1);*/
/*    	//print_preorder(tree);*/
/*	//draw_tree(tree,0,0,1);*/
/*	destroy_tree(tree);*/
/*	//fclose(file);			//Uncomment if you are using files to output MATLAB syntax*/
/*	//fclose(drawing);*/
/*	//fclose(fp);*/
/*}*/



/**************** Added for use with Chebychev branching strategy 10/15/14 ********************/
double *sol = NULL;
double solution[2] = {0.,0.};

double *find_leaf_in_bounds(double x_lb,double x_ub,double y_lb,double y_ub,node *cur_node)
{
	solution[0] = x_ideal;
	solution[1] = y_ideal;
	if(x_ub > cur_node->nw_x)
	{
		if(y_lb < cur_node->nw_y)
		{
			if(!cur_node->left && !cur_node->right)
			{
				if(x_lb < cur_node->nw_x)
				{
					if(y_ub > cur_node->nw_y)
					{
						solution[0] = cur_node->nw_x;
						solution[1] = cur_node->nw_y;
						return solution;
					}
					else return solution;
				}
				else return solution;
			}
			else 
			{
				if(cur_node->right) 
				{
					sol = find_leaf_in_bounds(x_lb,x_ub,y_lb,y_ub,cur_node->right);
					solution[0] = sol[0];
					solution[1] = sol[1];
				}
			}
		}
		else if(x_lb < cur_node->nw_x)
		{
			if(!cur_node->left && !cur_node->right) return solution;
			else 
			{
				if(cur_node->left)
				{
					sol = find_leaf_in_bounds(x_lb,x_ub,y_lb,y_ub,cur_node->left);
					solution[0] = sol[0];
					solution[1] = sol[1];
				}
			}
		}
	}
	else if(y_ub > cur_node->nw_y)
	{
		if(!cur_node->left && !cur_node->right)	return solution;
		else 
		{
			if(cur_node->left) 
			{
				sol = find_leaf_in_bounds(x_lb,x_ub,y_lb,y_ub,cur_node->left);
				solution[0] = sol[0];
				solution[1] = sol[1];
			}
		}
	}
	return solution;
}

node *find_first_node_left_of_val(double val, double lower_val, node *cur_node)
{
	node *temp_node = NULL;
/*	printf("current node: %lf,%lf - %lf,%lf \n",x_ideal - cur_node->se_x,y_ideal - cur_node->se_y,x_ideal - cur_node->nw_x,y_ideal - cur_node->nw_y);*/
	if( (x_ideal - cur_node->nw_x - val <= .0000001 || x_ideal - cur_node->se_x - val < -.0001) && y_ideal - cur_node->se_y - lower_val >= .0000001)
	{
		if(cur_node->left) temp_node = find_first_node_left_of_val(val, lower_val, cur_node->left);
	}
	else if(cur_node->right) temp_node = find_first_node_left_of_val(val, lower_val, cur_node->right);
/*	if(temp_node) printf("temp node: %lf,%lf - %lf,%lf \n",x_ideal-temp_node->se_x,y_ideal -temp_node->se_y,x_ideal -temp_node->nw_x,y_ideal -temp_node->nw_y);*/
	if(temp_node && (x_ideal - temp_node->nw_x - val <= .0000001 || x_ideal - temp_node->se_x - val < -.0001) && 
		y_ideal - temp_node->se_y - lower_val >= .0000001) 
	{
/*		printf("returning temp node\n");*/
		return temp_node;
	}
	else if( (x_ideal - cur_node->nw_x - val <= .0000001 || x_ideal - cur_node->se_x - val < -.0001) && y_ideal - cur_node->se_y - lower_val >= .0000001)
	{
/*		printf("returning current node\n");*/
		return cur_node;
	}
	else 
	{
/*		printf("returning NULL\n");*/
		return NULL;
	}
}

closest_nodes *find_two_nodes_left_of_val(double val, double lower_val, node *cur_node)
{
	node *closest = find_first_node_left_of_val(val, lower_val, cur_node);
	node *next = closest;
	if(closest)
	{
		if(closest->right) next = find_leftmost_leaf(closest->right);
		else 
		{
/*			printf("using parents for next\n");*/
	/*		print_preorder(tree,NULL);*/
			while(next->parent)
			{
				 if(next == next->parent->right) next = next->parent;
				 else break;
			}
			if(next->parent) next = next->parent;
		}
	}
	else return NULL;
	closest_nodes *ret_nodes = (struct closest_nodes*) malloc( sizeof( struct closest_nodes ) ); 
	ret_nodes->closest = closest;
	ret_nodes->next = next;
	return ret_nodes;
}

node *find_first_node_right_of_val(double val, double upper_val, node *cur_node)
{
	node *temp_node = NULL;
	if( (x_ideal - cur_node->se_x - val >= -.0000001 || x_ideal - cur_node->nw_x - val > .0001) && y_ideal - cur_node->nw_y - upper_val <= -.0000001)
	{
		if(cur_node->right) temp_node = find_first_node_right_of_val(val, upper_val, cur_node->right);
	}
	else if(cur_node->left) temp_node = find_first_node_right_of_val(val, upper_val, cur_node->left);
	if(temp_node && (x_ideal - temp_node->se_x - val >= -.0000001 || x_ideal - temp_node->nw_x - val > .0001) && 
		y_ideal - temp_node->nw_y - upper_val <= -.0000001) return temp_node;
	else if( (x_ideal - cur_node->se_x - val >= -.0000001 || x_ideal - cur_node->nw_x - val > .0001) && 
		y_ideal - cur_node->nw_y - upper_val <= -.0000001) 	return cur_node;
	else return NULL;
}

closest_nodes *find_two_nodes_right_of_val(double val, double upper_val, node *cur_node)
{
	node *closest = find_first_node_right_of_val(val, upper_val, cur_node);
	node *next = closest;
	if(closest)
	{
		if(closest->left) 
		{
/*			printf("using rightmost leaf of left node for next\n");*/
			next = find_rightmost_leaf(closest->left);
		}
		else 
		{
/*			printf("using parents for next\n");*/
	/*		print_preorder(tree,NULL);*/
			while(next->parent)
			{
				 if(next == next->parent->left) next = next->parent;
				 else break;
			}
			if(next->parent) next = next->parent;
		}
	}
	else return NULL;
	closest_nodes *ret_nodes = (struct closest_nodes*) malloc( sizeof( struct closest_nodes ) ); 
	ret_nodes->closest = closest;
	ret_nodes->next = next;
	return ret_nodes;
}

void reduce_box2(double *x, double *y, node *leaf);

void reduce_box(double *x, double *y, node *leaf)
{
	start_struct_timer = clock();
	if( clock_gettime( CLOCK_REALTIME, &start2) == -1 ) {
      		perror( "clock gettime" );
      		exit( EXIT_FAILURE );
    	}
	reduce_box2(x, y, leaf);
	finish_struct_timer = clock();
	if( clock_gettime( CLOCK_REALTIME, &stop2) == -1 ) {
      		perror( "clock gettime" );
      		exit( EXIT_FAILURE );
    	}

    	accum = ( stop2.tv_sec - start2.tv_sec )
          	+ ( stop2.tv_nsec - start2.tv_nsec )
            	/ BILLION;
     	struct_time2 += accum;
        struct_time += (double)(finish_struct_timer - start_struct_timer) / CLOCKS_PER_SEC;
}

void reduce_box2(double *x, double *y, node *leaf)
{
/*	printf("comparing box against: %d -- %lf,%lf to %lf,%lf\n",leaf->type,x_ideal-leaf->se_x,y_ideal-leaf->se_y,x_ideal-leaf->nw_x,y_ideal-leaf->nw_y);*/
	if(y_ideal - leaf->nw_y >= y[1])
	{
		if(x_ideal - leaf->nw_x >= x[0]) 
		{
			x[0] = x_ideal - leaf->nw_x + .001;
/*			printf("changing box x lb\n");*/
		}
		if(leaf->left) reduce_box2(x, y, leaf->left);
	}
	//else if(x_ideal - leaf->nw_x <= x[1] && y_ideal - leaf->nw_y >= y[0] && leaf->left) reduce_box2(x, y, leaf->left);
	//else if(leaf && leaf->left) reduce_box2(x, y, leaf->left);
	else if(x_ideal - leaf->nw_x <= x[1] && y_ideal - leaf->nw_y >= y[0] && leaf->right) 
	{
/*		printf("going right\n");*/
		reduce_box2(x, y, leaf->right);
	}
	
	if(x_ideal - leaf->se_x >= x[1])
	{
		if(y_ideal - leaf->se_y >= y[0]) 
		{
			y[0] = y_ideal - leaf->se_y;
/*			printf("changing box y lb\n");*/
		}
		if(leaf->right) reduce_box2(x, y, leaf->right);
	}
	//else if(x_ideal - leaf->se_x >= x[0] && y_ideal - leaf->se_y <= y[1] && leaf->right) reduce_box2(x, y, leaf->right);
	//else if(leaf && leaf->right) reduce_box2(x, y, leaf->right);
	else if(x_ideal - leaf->se_x >= x[0] && y_ideal - leaf->se_y <= y[1] && leaf->left) 
	{
/*		printf("going left\n");*/
		reduce_box2(x, y, leaf->left);
	}
}

dual_bd *first_dual_bd = NULL;
dual_bd *last_dual_bd = NULL;

void store_dual_bd(int seqnum1, int seqnum2, double x1, double y1, double x2, double y2)
{
/*	printf("adding a dual bound for seqnums %d and %d:\n",seqnum1,seqnum2);*/
/*	printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",x1,x2,y1,y2);*/
	dual_bd *new_dual_bd = (struct dual_bd*) malloc( sizeof( struct dual_bd ) ); 
	new_dual_bd->seqnum1 = seqnum1;
	new_dual_bd->seqnum2 = seqnum2;
	new_dual_bd->x1 = x1;
	new_dual_bd->x2 = x2;
	new_dual_bd->y1 = y1;
	new_dual_bd->y2 = y2;
	new_dual_bd->next = NULL;
	new_dual_bd->prev = NULL;
	if(!first_dual_bd)
	{
/*		printf("adding as first\n");*/
		first_dual_bd = new_dual_bd;
		first_dual_bd->prev = NULL;
		first_dual_bd->next = NULL;
	}
	else if(!last_dual_bd)
	{
/*		printf("adding as last\n");*/
		last_dual_bd = new_dual_bd;
		last_dual_bd->prev = first_dual_bd;
		first_dual_bd->next = last_dual_bd;
		last_dual_bd->next = NULL;
	}
	else
	{
/*		printf("replacing the last\n");*/
		last_dual_bd->next = new_dual_bd;
		new_dual_bd->prev = last_dual_bd;
		last_dual_bd = new_dual_bd;
	}
}

void remove_dual_bd(int seqnum)
{
/*	printf("removing dual bd for seqnum %d\n",seqnum);*/
	if(seqnum == 0) return;
	if(first_dual_bd && seqnum == first_dual_bd->seqnum1 || seqnum == first_dual_bd->seqnum2)
	{
/*		printf("the stored seqnums: %d and %d\n",first_dual_bd->seqnum1,first_dual_bd->seqnum2);*/
		if(seqnum == first_dual_bd->seqnum1) first_dual_bd->seqnum1 = -1;
		else first_dual_bd->seqnum2 = -1;
		if(first_dual_bd->seqnum1 == -1 && first_dual_bd->seqnum2 == -1)
		{
			dual_bd *temp = first_dual_bd;
			first_dual_bd = first_dual_bd->next;
			if(first_dual_bd) first_dual_bd->prev = NULL;
			if(first_dual_bd == last_dual_bd) last_dual_bd = NULL;
			free(temp);
/*			printf("the bound is actually being deleted (1)\n");*/
		}
	}
	else if(last_dual_bd && seqnum == last_dual_bd->seqnum1 || seqnum == last_dual_bd->seqnum2)
	{
		if(seqnum == last_dual_bd->seqnum1) last_dual_bd->seqnum1 = -1;
		else last_dual_bd->seqnum2 = -1;
		if(last_dual_bd->seqnum1 == -1 && last_dual_bd->seqnum2 == -1)
		{
			dual_bd *temp = last_dual_bd;
			if(last_dual_bd->prev != first_dual_bd) last_dual_bd = last_dual_bd->prev;
			else 
			{
				first_dual_bd->next = NULL;
				last_dual_bd = NULL;
			}
			if(last_dual_bd) last_dual_bd->next = NULL;
/*			if(first_dual_bd == last_dual_bd) first_dual_bd = NULL;*/
			free(temp);
/*			printf("the bound is actually being deleted (2)\n");*/
		}
	}
	else
	{
		dual_bd *temp = last_dual_bd->prev;
		while(seqnum != temp->seqnum1 && seqnum != temp->seqnum2) temp = temp->prev;
		if(seqnum == temp->seqnum1) temp->seqnum1 = -1;
		else temp->seqnum2 = -1;
		if(temp->seqnum1 == -1 && temp->seqnum2 == -1)
		{
			temp->prev->next = temp->next;
			temp->next->prev = temp->prev;
			free(temp);
/*			printf("the bound is actually being deleted (3)\n");*/
		}
	}
}

void print_current_dual_bd()
{
	dual_bd *temp = first_dual_bd;
	if(!temp && last_dual_bd) temp = last_dual_bd;
	if(temp)
	{
/*		printf("plot([%lf,%lf],[%lf,%lf],'-bo');\n",temp->x1,temp->x2,temp->y1,temp->y2);*/
		while(temp->next)
		{
			temp = temp->next;
/*			printf("plot([%lf,%lf],[%lf,%lf],'-bo');\n",temp->x1,temp->x2,temp->y1,temp->y2);*/
		}
	}
}

void create_dual_bd()
{
	dual_bd *temp = first_dual_bd;
	if(!temp && last_dual_bd) temp = last_dual_bd;
	if(temp)
	{
		insert2(2,temp->x1,temp->y1,temp->x2,temp->y2,(temp->y1 - temp->y2)/(temp->x1 - temp->x2),&tree,&empty_node2);
		while(temp->next)
		{
			insert2(2,temp->x1,temp->y1,temp->x2,temp->y2,(temp->y1 - temp->y2)/(temp->x1 - temp->x2),&tree,&empty_node2);
			temp = temp->next;
		}
	}
}

double calculate_hd_dist(struct node *primal, struct node *dual, double best_val)
{	
	double new_val = 0., distance = 0.;
  	if (dual)
  	{
/*  		printf("got a dual soln\n");*/
		if(dual->subtree_size > -1 && dual->left)
  		{
     			new_val = calculate_hd_dist(primal,dual->left,best_val);
     			best_val = fmax(best_val,new_val);
/*     			printf("after comparison, best value: %lf\n",best_val);*/
/*     			printf("\n");*/
     		}
     		if(dual->subtree_size > -1 && dual->right)
     		{
      			new_val = calculate_hd_dist(primal,dual->right,best_val);
      			best_val = fmax(best_val,new_val);
/*      			printf("after comparison, best value: %lf\n",best_val);*/
/*      			printf("\n");*/
      		}
      		dual->subtree_size = -1;
/*      		printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-dual->nw_x,x_ideal-dual->se_x,y_ideal-dual->nw_y,y_ideal-dual->se_y);      		*/
	      		
      		int insert_check = 0;
      		if(dual->type == 1) insert_check = mock_insert2(1,x_ideal-dual->nw_x,y_ideal-dual->nw_y,0.,0.,0.,&tree);
      		else insert_check = mock_insert2(2,x_ideal-dual->nw_x,y_ideal-dual->nw_y,x_ideal-dual->se_x,y_ideal-dual->se_y,
      			(dual->nw_y-dual->se_y)/(dual->nw_x-dual->se_x),&tree);
      		if(!insert_check) 
      		{
/*      			printf("returning %lf (%d)\n",best_val,__LINE__);*/
      			return (double) best_val;
      		}
      		if (primal)
		{
/*			printf("plot([%lf,%lf],[%lf,%lf],'-o');\n",x_ideal-primal->nw_x,x_ideal-primal->se_x,y_ideal-primal->nw_y,y_ideal-primal->se_y);  */
/*			printf("got a primal soln\n");*/
			if(dual->se_x < primal->nw_x)
			{
				if(primal->left)
				{
					new_val = calculate_hd_dist(primal->left,dual,best_val);
		     			if(best_val == 0.) best_val = new_val;
		     			else best_val = fmax(best_val,new_val);
	     			}
			}
			else if(dual->nw_x > primal->se_x)
			{
				if( primal->right)
				{
					new_val = calculate_hd_dist(primal->right,dual,best_val);
		      			if(best_val == 0.) best_val = new_val;
		     			else best_val = fmax(best_val,new_val);
				}
			}
			else
			{
				if(primal->left)
		  		{
		     			new_val = calculate_hd_dist(primal->left,dual,best_val);
		     			if(best_val == 0.) best_val = new_val;
		     			else best_val = fmax(best_val,new_val);
/*		     			printf("after comparison, best value: %lf\n",best_val);*/
		     		}
		     		if(primal->right)
		     		{
		      			new_val = calculate_hd_dist(primal->right,dual,best_val);
		      			if(best_val == 0.) best_val = new_val;
		     			else best_val = fmax(best_val,new_val);
/*		      			printf("after comparison, best value: %lf\n",best_val);*/
		      		}
		      		
		      		if(primal->type == 1) // The inserted solution is a point.
				{
					if(dual->type == 1) // The current solution is a point.
					{
						distance = sqrt( (primal->nw_x - dual->nw_x)*(primal->nw_x - dual->nw_x) + 
							(primal->nw_y - dual->nw_y)*(primal->nw_y - dual->nw_y) );
/*						printf("comparing two points\nThe distance: %lf\n",distance);*/
/*						*/
/*						printf("returning (%d)\n",__LINE__);*/
						return (double) distance;
					}
					else
					{
						double temp_nw_x = fmax(primal->nw_x,dual->nw_x);
						double temp_nw_y = dual->nw_y;
						if(temp_nw_x != dual->nw_x) temp_nw_y = dual->se_y + dual->slope*(temp_nw_x - dual->se_x);
						double temp_se_y = fmax(primal->se_y,dual->se_y);
						double temp_se_x = dual->se_x;
						if(temp_se_y != dual->se_y) temp_se_x = dual->se_x + (1./dual->slope)*(temp_se_y - dual->se_y);
					
						distance = fmax( sqrt( (primal->nw_x - temp_nw_x)*(primal->nw_x - temp_nw_x) + 
							(primal->nw_y - temp_nw_y)*(primal->nw_y - temp_nw_y) ), 
							sqrt( (primal->nw_x - temp_se_x)*(primal->nw_x - temp_se_x) + 
							(primal->nw_y - temp_se_y)*(primal->nw_y - temp_se_y) ) );
/*						printf("comparing a point and a segment\nThe distance: %lf\n",distance);*/

/*						printf("returning (%d)\n",__LINE__);*/
						return (double) distance;
					}
				}
				else
				{
					if(dual->type == 1)
					{
						distance = fmax( sqrt( (primal->nw_x - dual->nw_x)*(primal->nw_x - dual->nw_x) + 
							(primal->nw_y - dual->nw_y)*(primal->nw_y - dual->nw_y) ), 
							sqrt( (primal->se_x - dual->nw_x)*(primal->se_x - dual->nw_x) +
							(primal->se_y - dual->nw_y)*(primal->se_y - dual->nw_y) ) );
						
/*						printf("comparing a point and a segment\nThe distance: %lf\n",distance);*/
/*						printf("returning (%d)\n",__LINE__);*/
					 	return (double) distance;
				 	}
				 	else
				 	{
				 		double dist1 = 0., dist2 = 0, dist3 = 0., dist4 = 0.;
						double x_intersect = ((1./dual->slope)*primal->se_x+dual->slope*dual->se_x+primal->se_y-dual->se_y)/
														(dual->slope+(1./dual->slope)); // proj from se_x
						double y_intersect = dual->nw_y + dual->slope*(x_intersect - dual->nw_x);
						if(y_intersect != y_intersect) y_intersect = primal->se_y + (-1./dual->slope)*(x_intersect - primal->se_x);
				
						if(x_intersect <= dual->se_x && x_intersect >= dual->nw_x) dist1 = 
							sqrt( (primal->se_x - x_intersect)*(primal->se_x - x_intersect) + 
							(primal->se_y - y_intersect)*(primal->se_y - y_intersect) );

						else dist1 = fmin( sqrt( (primal->se_x - dual->nw_x)*(primal->se_x - dual->nw_x) + 
							(primal->se_y - dual->nw_y)*(primal->se_y - dual->nw_y) ), 
							sqrt( (primal->se_x - dual->se_x)*(primal->se_x - dual->se_x) + 
							(primal->se_y - dual->se_y)*(primal->se_y - dual->se_y) ) );
						
						x_intersect = ((1./dual->slope)*primal->nw_x+dual->slope*dual->se_x+primal->nw_y-dual->se_y)/
														(dual->slope+(1./dual->slope)); // proj from nw_x
						y_intersect = dual->nw_y + dual->slope*(x_intersect - dual->nw_x);
						if(y_intersect != y_intersect) y_intersect = primal->nw_y + (-1./dual->slope)*(x_intersect - primal->nw_x);
					
				
						if(x_intersect <= dual->se_x && x_intersect >= dual->nw_x) dist2 = 
							sqrt( (primal->nw_x - x_intersect)*(primal->nw_x - x_intersect) + 
							(primal->nw_y - y_intersect)*(primal->nw_y - y_intersect) );
						else dist2 = fmin( sqrt( (primal->nw_x - dual->nw_x)*(primal->nw_x - dual->nw_x) + 
							(primal->nw_y - dual->nw_y)*(primal->nw_y - dual->nw_y) ),
							sqrt( (primal->nw_x - dual->se_x)*(primal->nw_x - dual->se_x) + 
							(primal->nw_y - dual->se_y)*(primal->nw_y - dual->se_y) ) );
						
						x_intersect = ((primal->slope)*primal->nw_x+(1./primal->slope)*dual->se_x-primal->nw_y+dual->se_y)/
							((1./primal->slope)+(primal->slope));  
						y_intersect = dual->se_y + (-1./primal->slope)*(x_intersect - dual->se_x);
						if(y_intersect != y_intersect) y_intersect = primal->nw_y + (primal->slope)*(x_intersect - primal->nw_x);
					
				
						if(x_intersect <= primal->se_x && x_intersect >= primal->nw_x) dist3 = 
							sqrt( (dual->se_x - x_intersect)*(dual->se_x - x_intersect) + 
							(dual->se_y - y_intersect)*(dual->se_y - y_intersect) );
						else dist3 = fmin( sqrt( (primal->nw_x - dual->se_x)*(primal->nw_x - dual->se_x) + 
							(primal->nw_y - dual->se_y)*(primal->nw_y - dual->se_y) ),
							sqrt( (primal->se_x - dual->se_x)*(primal->se_x - dual->se_x) + 
							(primal->se_y - dual->se_y)*(primal->se_y - dual->se_y) ) );
						
						x_intersect = ((primal->slope)*primal->nw_x+(1./primal->slope)*dual->nw_x-primal->nw_y+dual->nw_y)/
							((1./primal->slope)+(primal->slope)); 
						y_intersect = dual->nw_y + (-1./primal->slope)*(x_intersect - dual->nw_x);
						if(y_intersect != y_intersect) y_intersect = primal->nw_y + (primal->slope)*(x_intersect - primal->nw_x);
					
				
						if(x_intersect <= primal->se_x && x_intersect >= primal->nw_x) dist4 = 
							sqrt( (dual->nw_x - x_intersect)*(dual->nw_x - x_intersect) +
							(dual->nw_y - y_intersect)*(dual->nw_y - y_intersect) );
						else dist4 = fmin( sqrt( (primal->nw_x - dual->nw_x)*(primal->nw_x - dual->nw_x) + 
							(primal->nw_y - dual->nw_y)*(primal->nw_y - dual->nw_y) ),
							sqrt( (primal->se_x - dual->nw_x)*(primal->se_x - dual->nw_x) + 
							(primal->se_y - dual->nw_y)*(primal->se_y - dual->nw_y) ) );
						
						dist1 = fmax(dist1,dist2);
						dist3 = fmax(dist3,dist4);
						distance = fmax(dist1,dist3);
					
/*						printf("returning (%d)\n",__LINE__);*/
/*						printf("comparing two segments\nThe distance: %lf\n",distance);*/
					
						return (double) distance;
				 	}
				}	
			}
		}
/*		printf("returning %lf (%d)\n",best_val,__LINE__);*/
/*		if(best_val == 0) best_val = max_range;*/
		return (double) best_val;
	}
/*	printf("returning %lf (%d)\n",best_val,__LINE__);*/
/*	if(best_val == 0) best_val = max_range;*/
	return (double) best_val;
}

int done_val2 = -2;
double calculate_max_proximal_hd_dist(struct node *primal, struct node *dual, double best_val)
{	
	double new_val = 0., distance = 0.;
  	if (dual)
  	{
/*  		printf("\tgot a dual soln\n");*/
		if(dual->subtree_size > done_val2 && dual->left)
  		{
     			new_val = calculate_max_proximal_hd_dist(primal,dual->left,best_val);
     			best_val = fmax(best_val,new_val);
     		}
     		if(dual->subtree_size > done_val2 && dual->right)
     		{
      			new_val = calculate_max_proximal_hd_dist(primal,dual->right,best_val);
      			best_val = fmax(best_val,new_val);
      		}
      		dual->subtree_size = done_val2;
/*      		printf("\tplot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-dual->nw_x,x_ideal-dual->se_x,y_ideal-dual->nw_y,y_ideal-dual->se_y);  */
/*      		printf("current best val: %lf\n",best_val);   */
      		done_val--; 		
	      		
      		int insert_check = 0;
      		if(dual->type == 1) insert_check = mock_insert2(1,x_ideal-dual->nw_x,y_ideal-dual->nw_y,0.,0.,0.,&tree);
      		else insert_check = mock_insert2(2,x_ideal-dual->nw_x,y_ideal-dual->nw_y,x_ideal-dual->se_x,y_ideal-dual->se_y,
      			(dual->nw_y-dual->se_y)/(dual->nw_x-dual->se_x),&tree);
      		if(!insert_check) 
      		{
/*      			printf("\treturning %lf (%d)\n",best_val,__LINE__);*/
      			return (double) best_val;
      		}
      		
      		if (primal) 
      		{
      			if(primal == tree) best_val = get_max_proximal_hd_dist(dual->type,x_ideal-dual->nw_x,y_ideal-dual->nw_y,x_ideal-dual->se_x,
      						y_ideal-dual->se_y,(dual->nw_y-dual->se_y)/(dual->nw_x-dual->se_x),&primal,best_val);
      			else best_val = get_max_proximal_hd_dist(dual->type,dual->nw_x,dual->nw_y,dual->se_x,dual->se_y,
      						(dual->nw_y-dual->se_y)/(dual->nw_x-dual->se_x),&primal,best_val);
      		}
/*      		printf("\treturning %lf (%d)\n",best_val,__LINE__);*/
		return (double) best_val;
	}
/*	printf("\treturning %lf (%d)\n",best_val,__LINE__);*/
	return (double) best_val;
}

double glbl_x, glbl_y;
double get_nadirs(node *n1, int starting, double best_val)
{	
	double new_val = 0.;
	if(starting)
	{
		glbl_x = n1->se_x;
		glbl_y = n1->se_y;
	}
    	if (n1)
    	{	
   		if(n1->left)
		{
			new_val = get_nadirs(n1->left,starting,best_val);
			best_val = fmax(new_val,best_val);
        	}
        	starting = 0;
		if(!starting) 
		{
/*			printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ideal-n1->nw_x,x_ideal-n1->se_x,y_ideal-n1->nw_y,y_ideal-n1->se_y);*/
/*			printf("plot(%lf,%lf,'go');\n",x_ideal-n1->nw_x,y_ideal-glbl_y);*/
			node *temp_node = (struct node*) malloc( sizeof( struct node ) );
			temp_node->type = 1; 
			temp_node->nw_x = n1->nw_x;
			temp_node->se_x = n1->nw_x;
			temp_node->nw_y = glbl_y;
			temp_node->se_y = glbl_y;
			temp_node->left = NULL;
			temp_node->right = NULL;
			temp_node->subtree_size = 0;
			done_val2--;
			double new_val = calculate_max_proximal_hd_dist(temp_node, tree2, best_val);
			best_val = fmax(new_val,best_val);
/*			printf("distance: %lf\n",best_val);*/
			glbl_x = n1->se_x;
			glbl_y = n1->se_y;
			free(temp_node);
		}
        	if(n1->right)
        	{
        		new_val = get_nadirs(n1->right,starting,best_val);
        		best_val = fmax(new_val,best_val);
        	}
   	}
   	starting = 0;
   	return best_val;
}

double get_length(node *n1, double val)
{
    	if (n1)
    	{	
		if(mock_insert2(n1->type,x_ideal-n1->nw_x,y_ideal-n1->nw_y,x_ideal-n1->se_x,y_ideal-n1->se_y,
      			(n1->nw_y-n1->se_y)/(n1->nw_x-n1->se_x),&tree) )
      		{
      			double length = sqrt( (n1->nw_x-n1->se_x)*(n1->nw_x-n1->se_x) + (n1->nw_y-n1->se_y)*(n1->nw_y-n1->se_y) );
			val += length;
		}
		
   		if(n1->left)
		{
			val = get_length(n1->left,val);
        	}
		
        	if(n1->right)
        	{
        		val = get_length(n1->right,val);
        	}
   	}
   	return val;
}

double *x_separators = NULL;
double *y_separators = NULL;
int *separator_dir = NULL;
int num_separators = 0, seps_been_malloced = 0, malloced_val = 0;

int find_separations(node *n1, int starting)
{	
	double new_val = 0.;
	if(starting)
	{
		glbl_x = n1->se_x;
		glbl_y = n1->se_y;
/*		points_only = 1;*/
/*		its_been_only_points = 1;*/
	}
    	if (n1)
    	{	
   		if(n1->left)
		{
			find_separations(n1->left,0);
        	}
		if(!starting) 
		{
			double x_sep = n1->nw_x - glbl_x;
			double y_sep = -n1->nw_y + glbl_y;
/*			printf("x sep: %lf, y_sep: %lf\n",x_sep,y_sep);*/
			
			double x_per = 100.*x_sep/x_range;
			double y_per = 100.*y_sep/y_range;
/*			printf("x per: %lf, y per: %lf\n",x_per,y_per);*/
			
			if((y_per > 1.5 && x_per > 5.) || (x_per > 1.5 && y_per > 5.))
			{
				if(!seps_been_malloced)
				{
					x_separators = (double*) malloc( 10*sizeof( double ) );
					y_separators = (double*) malloc( 10*sizeof( double ) );
					separator_dir = (int*) malloc( 10*sizeof( int ) );
					seps_been_malloced = 1;
					malloced_val = 10;
				}
				else if(num_separators >= malloced_val)
				{
					x_separators = (double*) realloc(x_separators, 2*num_separators*sizeof( double ) );
					y_separators = (double*) realloc(y_separators, 2*num_separators*sizeof( double ) );
					separator_dir = (int*) realloc(separator_dir, 2*num_separators*sizeof( int ) );
					malloced_val = 2*num_separators;
				}
				x_separators[num_separators] = x_ideal-n1->nw_x;
				y_separators[num_separators] = y_ideal-n1->nw_y;
				separator_dir[num_separators] = 0;
/*				printf("just created values for separator number %d\n",num_separators);*/
				num_separators++;
			}
/*			else if(y_per > 10.)*/
/*			{*/
/*				if(!seps_been_malloced)*/
/*				{*/
/*					x_separators = (double*) malloc( 10*sizeof( double ) );*/
/*					y_separators = (double*) malloc( 10*sizeof( double ) );*/
/*					separator_dir = (int*) malloc( 10*sizeof( int ) );*/
/*					seps_been_malloced = 1;*/
/*					malloced_val = 10;*/
/*				}*/
/*				else if(num_separators >= malloced_val)*/
/*				{*/
/*					x_separators = (double*) realloc(x_separators, 2*num_separators*sizeof( double ) );*/
/*					y_separators = (double*) realloc(y_separators, 2*num_separators*sizeof( double ) );*/
/*					separator_dir = (int*) realloc(separator_dir, 2*num_separators*sizeof( int ) );*/
/*					malloced_val = 2*num_separators;*/
/*				}*/
/*				x_separators[num_separators] = x_ideal-n1->nw_x;*/
/*				y_separators[num_separators] = y_ideal-n1->nw_y;*/
/*				separator_dir[num_separators] = 1;*/
/*				num_separators++;*/
/*			}*/
			
			if(n1->type != 1) 
			{
				points_only = 0;
				its_been_only_points = 0;
			}
			
			glbl_x = n1->se_x;
			glbl_y = n1->se_y;
		}
        	if(n1->right)
        	{
        		find_separations(n1->right,0);
        	}
   	}
   	starting = 0;
}

node *copy_tree(struct node *current, struct node *parent)
{
	node *copy = NULL;
	if(current)
	{
		copy = (struct node*) malloc( sizeof( struct node ) ); 
		if(current->left) copy->left = copy_tree(current->left,copy);
        	else copy->left = NULL;
        	copy->type = current->type;
		copy->nw_x = current->nw_x;
		copy->se_x = current->se_x;
		copy->nw_y = current->nw_y;
		copy->se_y = current->se_y;
		copy->slope = current->slope;
		copy->subtree_size = current->subtree_size;
		copy->parent = parent;    
        	if(current->right) copy->right = copy_tree(current->right,copy);
        	else copy->right = NULL;    
	}
	return copy;
}

node *last_used = NULL;

double get_hypervolume(node *n1, int starting, double current_val)
{	
	if (starting) last_used = NULL;
    	if (n1)
    	{	
   		if(n1->left)
		{
			//printf("left child:\n");
			current_val = fmax(current_val,get_hypervolume(n1->left,0,current_val));
        	}
		// Uncomment the following lines in order to print MATLAB syntax for plotting the Pareto optimal solutions
		if(last_used != NULL)
		{
			if(current_val <= 0.)
			{
				if(y_ideal-n1->nw_y > extreme_y)
				{
					current_val += fabs( (last_used->nw_y - (y_ideal-extreme_y))*(last_used->nw_x - (x_ideal-extreme_x)) );
/*					printf("fill([%lf,%lf,%lf,%lf],[%lf,%lf,%lf,%lf],'g');\n",x_ideal-last_used->nw_x,extreme_x,extreme_x,x_ideal-last_used->nw_x,y_ideal-last_used->nw_y,y_ideal-last_used->nw_y,extreme_y,extreme_y);*/
/*					printf("current_val: %lf\n",current_val);*/
				}
				if(last_used->type == 2) 
				{
					current_val += fabs(.5*( (last_used->se_x - (x_ideal-extreme_x)) + (last_used->nw_x - (x_ideal-extreme_x)) )*(last_used->se_y - last_used->nw_y));
/*					printf("fill([%lf,%lf,%lf,%lf],[%lf,%lf,%lf,%lf],'g');\n",x_ideal-last_used->nw_x,extreme_x,extreme_x,x_ideal-last_used->se_x,y_ideal-last_used->nw_y,y_ideal-last_used->nw_y,y_ideal-last_used->se_y,y_ideal-last_used->se_y);*/
/*					printf("current_val: %lf\n",current_val);*/
				}
			}
			if( fabs(n1->nw_y - last_used->se_y) > .000001) //&& fabs(n1->se_x - last_used->nw_x) > .000001)
			{
				current_val += fabs( (last_used->se_y - n1->nw_y)*(n1->nw_x - (x_ideal-extreme_x)) );
/*				printf("fill([%lf,%lf,%lf,%lf],[%lf,%lf,%lf,%lf],'g');\n",x_ideal-n1->nw_x,x_ideal-n1->nw_x,extreme_x,extreme_x,y_ideal-last_used->se_y,y_ideal-n1->nw_y,y_ideal-n1->nw_y,y_ideal-last_used->se_y);*/
/*				printf("current_val: %lf\n",current_val);*/
			}
			if(n1->type == 2) 
			{
				current_val += fabs(.5*( (n1->se_x - (x_ideal-extreme_x)) + (n1->nw_x - (x_ideal-extreme_x)) )*(n1->se_y - n1->nw_y));
/*				printf("fill([%lf,%lf,%lf,%lf],[%lf,%lf,%lf,%lf],'g');\n",x_ideal-n1->nw_x,extreme_x,extreme_x,x_ideal-n1->se_x,y_ideal-n1->nw_y,y_ideal-n1->nw_y,y_ideal-n1->se_y,y_ideal-n1->se_y);*/
/*				printf("current_val: %lf\n",current_val);*/
			}
/*			printf("plot([%lf,%lf],[%lf,%lf],'-ro');\n",n1->nw_x,n1->se_x,n1->nw_y,n1->se_y);*/
/*			printf("plot([%lf,%lf],[%lf,%lf],'-bo');\n",last_used->nw_x,last_used->se_x,last_used->nw_y,last_used->se_y);*/
/*			printf("NW_x: %lf\n",x_ideal-extreme_x);*/
/*			exit(0);*/
		}
		last_used = n1;
        	if(n1->right)
        	{
        		//printf("right child:\n");
		  	current_val = fmax(current_val,get_hypervolume(n1->right,0,current_val));
        	}
   	}
   	return current_val;
}
