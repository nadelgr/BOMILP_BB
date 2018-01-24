/* File created by Nathan Adelgren, Graduate Assisistant at Clemson University.
Started: 9/28/2014 
Finished: N/A
*/

#include<stdlib.h>
#include<stdio.h>
#include <math.h>
#include<time.h>

#include "cplex.h"
#include "callbacks.h"
#include "bb-bicriteria.h"
#include "max_tree.h"

int print_on = 0;
int num_integer = 10;

FILE *inserted_data2 = NULL;

branch_node *first = NULL;
branch_node *last = NULL;

pareto_branch_node *first_pareto = NULL;
pareto_branch_node *last_pareto = NULL;

/*node *potential_branch_tree = NULL;*/
int num_still_to_branch = 0;

double *ob_coef1=NULL;
double *ob_coef2=NULL;
double *weighted_coefs=NULL;
double x_ideal = 0.;
double y_ideal = 0.;
int node_iterations = 0;
char *xctype = NULL;

int cur_numcols ;//= -1;
int cur_numrows ;//= -1;
int sub_prob_numrows;
int original_num_rows = -1;
int cur_split_dir = 0;

int first_new_row_index = -1;

int prev_seqnum = -1;
double prev_slope = 0.;
int same_seq_flag = 0;
int next_seqnum = -1;
int need_to_fathom = 0;
int fathom_seqnum = -1;
int reduce_to_lower = 0;
int lower_seqnum = -1;
int reduce_to_left = 0;
int left_seqnum = -1;
int num_stored = 4;
/*int psa_complete_snum = 0;*/

int fathomed_by_dominated_local_ideal_pt = 0;
int fathomed_by_dominated_local_ideal_segment = 0;

int fathomed_by_dominated_search_region = 0;
int fathomed_by_dominated_local_ideal_pts = 0;
int fathomed_by_dominated_local_ideal_segments = 0;
int fathomed_by_1_dominated_pt_1_dominated_segment = 0;
int fathomed_by_PSA_completion = 0;
int fathomed_by_dominated_lb = 0;
int num_nodes_with_mips_solved = 0;

double sub_pr1_x_ub = 0;
double sub_pr1_y_lb = 0;
double sub_pr1_x_lb = 0;
double sub_pr1_y_ub = 0;
double sub_pr2_x_lb = 0;
double sub_pr2_y_ub = 0;
double sub_pr2_x_ub = 0;
double sub_pr2_y_lb = 0;

int reduced_subprob_seqnum = 0;
double reduced_subprob_x_ub = 0.;
double reduced_subprob_y_ub = 0.;
double reduced_subprob_x_lb = 0.;
double reduced_subprob_y_lb = 0.;

int psa_int_pivot = 0;
int num_x_frac_at_both = 0;
int frac_index = -1;
double frac_val = 0.;
double branch_seqnum = -1.;

int pareto_branching = 0;
int pareto_start = 50;
int pareto_limit = 500;
double prob_width = 0.;
int width_divisor = 100;

CPXLPptr lp_1 = NULL;
CPXLPptr lp_2 = NULL;
CPXLPptr mip1 = NULL;

double **stored_x = NULL;
int x_rotation;
int exact_mips = 0;

/*double time_limit = 60.;*/

//int insert_to_potential_branch_tree = 0;

clock_t start_heur, finish_heur;

int been_passed = 0;

int *integer_indices = NULL;
int *integer_var_scores = NULL;
int *integer_var_last_val = NULL;
int total_num_integer = 0;

int bds_reduced_by_dominated_ideal_pt = 0;
int bds_reduced_by_dominated_ideal_segment = 0;
int bds_reduced_by_dominated_lb = 0;
int bds_reduced_by_infeasibility = 0;
int bds_reduced = 0;

/*int only_infeasibility = 0;*/
/*int check_bound = 0;*/
int check_bound_in_nodeop = 1;
/*int bd_reduction_during_branching = 1;*/
int slope_scores_in_psa = 0;

CPXENVptr  env2=NULL;
CPXENVptr  env3=NULL;
CPXENVptr  env_global=NULL;

double multiplier = 1.;
int num_frac = 0;
double within_PSA_score = 5.;

user_data *userhandle_cutcallback_zero = NULL;
int current_userhandle_was_empty = 0;
/*int points_only = 0, its_been_only_points = 0;*/

/*double *start_x;*/
/*double *start_slack;*/

int *global_beg;
int *global_varindices; 
double *global_values; 
int *global_effortlevel;
int global_num_starts;
int global_startspace;

CPXLPptr   global_mip=NULL;
CPXENVptr  env_just_solve_mips=NULL;

int check_for_stopping_PSA_full = 0;
/*int check_for_early_PSA_reduce = 0;*/
/*int check_for_early_PSA_reduce_after_n_iter = 1;*/
/*int n_iter = 15;*/

int prev_numsolns = 0;
int times_to_run = 0;
int remove_dominated_middle = 1;

clock_t start_BB, finish_BB;
/*double max_time = 43200.;*/
/*double max_time = 2400.;*/
double duration_BB = 0.;

/*int limit_bd_reduction = 0;*/

double time_processing_nodes = 0.;
double time_branching = 0.;
double time_solving_mips = 0.;
double time_generating_cuts = 0.;
double time_tightening_variable_bounds = 0;

double ob1_offset = 0., ob2_offset = 0.;

double NW_extreme_x = 0.;
double NW_extreme_y = 0.;
double SE_extreme_x = 0.;
double SE_extreme_y = 0.;

int cut_prob_first_row_to_change_index = 0;
int cut_prob_second_row_to_change_index = 0;
int *cut_prob_bd_indices_disj1 = NULL;
int *cut_prob_bd_indices_disj2 = NULL;
int cut_prob_first_pi_index = 0;
int *cut_prob_pi_indices = NULL;
CPXLPptr cut_problem=NULL;

double *temp_x_r = NULL;
double *temp_x_l = NULL;

double time_generating_disjunction_cuts = 0.;
int number_disj_cuts = 0;
int number_pareto_branches = 0;

int *dominated_indices = NULL;
int *dominating_indices = NULL;
int num_dom_cols = 0;
double time_per = 0.;
int approximate_dual_bd = 0;
int quit_generating_dual_bd = 0;

double max_time_to_solve_a_mip = 0.;

/*double max_time_phase_2 = 20.;*/
/*double max_time_bd_red_b4_branching = 1.;*/

/*********************************************************************************************************************** 

	This function is used for decoding the scoring scheme used for determining the branching variable.
	
***********************************************************************************************************************/

double convert_it (int n) {
    if (n < 10) return .1;
    if (n < 100) return .01;
    if (n < 1000) return .001;
    if (n < 10000) return .0001;
    if (n < 100000) return .00001;
    if (n < 1000000) return .000001;
    if (n < 10000000) return .0000001;
    if (n < 100000000) return .00000001;
    if (n < 1000000000) return .000000001;
    return .0000000001;
}

void give_env(CPXENVptr env)
{
	env_global = env;
}

void PSA_all(CPXCENVptr env, CPXLPptr prob);
void chg_coefs(CPXCENVptr env, CPXLPptr prob, int *indices, double slope);

int stopped_early_last_time = 0, index_we_left_off_at = 0;

/*CPXLPptr create_cut_prob(CPXCENVptr env, CPXLPptr prob)*/
/*{*/
/*	int status = 0, i = 0, z = 0, nzcnt = 0, surplus = 0, cmatspace = 0;*/
/*	*/
/*	double *lb = (double *) malloc( (cur_numrows+2*cur_numcols)*sizeof(double) );*/
/*	double *ub = (double *) malloc( (cur_numrows+2*cur_numcols)*sizeof(double) );*/
/*	char *sense = (char *) malloc ( (cur_numrows+2*cur_numcols)*sizeof(char) );*/
/*	char *lu = (char *) malloc ( (cur_numcols+1)*sizeof(char) );*/
/*	int *rmatbeg = (int *) malloc ( 2*cur_numcols*sizeof(int) );*/
/*	int *rmatind = (int *) malloc ( 4*cur_numcols*sizeof(int) );*/
/*	double *rmatval = (double *) malloc ( 4*cur_numcols*sizeof(double) );*/
/*	double *rhs = (double *) malloc ( cur_numrows*sizeof(double) );*/
/*	double *zeros = (double *) calloc( cur_numrows, sizeof(double) );*/
/*	int *indices = (int *) malloc (fmax(cur_numrows,cur_numcols+1)*sizeof(int) );*/
/*	int *col_list = (int *) malloc (cur_numrows*sizeof(int) );*/
/*	*/
/*	CPXLPptr cut_prob = CPXcloneprob (env, prob, &status);*/
/*	if ( status ) {*/
/*   		printf ("%s(%d): Failed to copy prob, error code %d\n", __FILE__, __LINE__, status);*/
/*  	}*/
/*  	*/
/*  	status = ( CPXgetub (env, cut_prob, ub, 0, cur_numcols-1) || CPXgetlb (env, cut_prob, lb, 0, cur_numcols-1) );*/
/*  	if ( status ) {*/
/*   		printf ("%s(%d): Failed to get variable bounds, error code %d\n", __FILE__, __LINE__, status);*/
/*  	}*/
/*  	*/
/*  	status = CPXgetrhs (env, cut_prob, rhs, 0, cur_numrows-1);*/
/*  	if ( status ) {*/
/*   		printf ("%s(%d): Failed to get rhs, error code %d\n", __FILE__, __LINE__, status);*/
/*  	}*/
/*  	*/
/*  	for(i=0;i<cur_numrows;i++)*/
/*  	{*/
/*  		indices[i] = i;*/
/*  		col_list[i] = cur_numcols;*/
/*  	}*/
/*  	status = CPXchgrhs (env, cut_prob, cur_numrows, indices, zeros);*/
/*  	if ( status ) {*/
/*   		printf ("%s(%d): Failed to chg rhs, error code %d\n", __FILE__, __LINE__, status);*/
/*  	}*/
/*  	*/
/*  	for(i=0;i<2*cur_numcols;i++)*/
/*  	{*/
/*  		sense[i] = 'G';*/
/*  		rmatbeg[i] = 2*i;*/
/*  		if(i < cur_numcols) rmatind[2*i] = i;*/
/*  		else rmatind[2*i] = i - cur_numcols;*/
/*  		rmatind[2*i+1] = cur_numcols;*/
/*  		if(i < cur_numcols)*/
/*  		{*/
/*  			rmatval[2*i] = 1.;*/
/*  			rmatval[2*i+1] = ub[i];*/
/*  		}*/
/*  		else*/
/*  		{*/
/*  			rmatval[2*i] = -1.;*/
/*  			rmatval[2*i+1] = -lb[i-cur_numcols];*/
/*  		}*/
/*  	}*/
/*  	*/
/*  	status = CPXchgsense (env, cut_prob, cur_numrows, indices, sense);*/
/*  	*/
/*  	status = CPXaddrows (env, cut_prob, 1, 2*cur_numcols, 4*cur_numcols, NULL, sense, rmatbeg, rmatind, rmatval, NULL, NULL);*/
/*  	if ( status ) {*/
/*   		printf ("%s(%d): Failed to add rows, error code %d\n", __FILE__, __LINE__, status);*/
/*  	}*/
/*  	*/
/*  	status = CPXchgcoeflist (env, cut_prob, cur_numrows, indices, col_list, rhs);*/
/*  	if ( status ) {*/
/*   		printf ("%s(%d): Failed to change col coefs, error code %d\n", __FILE__, __LINE__, status);*/
/*  	}*/
/*  	*/
/*  	for(i=0;i<cur_numrows+2*cur_numcols;i++)*/
/*  	{*/
/*  		lb[i] = -CPX_INFBOUND;*/
/*  		ub[i] = CPX_INFBOUND;*/
/*  	}*/
/*  	if(cur_numcols + 1 > cur_numrows) for(i=cur_numrows;i<cur_numcols+1;i++) indices[i] = i;*/
/*  	*/
/*  	for(i=0;i<cur_numcols+1;i++) lu[i] = 'L';*/
/*  	*/
/*  	status = CPXchgbds (env, cut_prob, cur_numcols+1, indices, lu, lb);*/
/*  	if ( status ) {*/
/*   		printf ("%s(%d): Failed to change lbs, error code %d\n", __FILE__, __LINE__, status);*/
/*  	}*/
/*  	*/
/*  	for(i=0;i<cur_numcols+1;i++) lu[i] = 'U';*/
/*  	*/
/*  	status = CPXchgbds (env, cut_prob, cur_numcols+1, indices, lu, ub);*/
/*  	if ( status ) {*/
/*   		printf ("%s(%d): Failed to change ubs, error code %d\n", __FILE__, __LINE__, status);*/
/*  	}*/
/*  	*/
/*  	for(i=2*cur_numcols;i<cur_numrows+2*cur_numcols;i++) sense[i] = 'G';*/
/*  	*/
/*  	status = CPXaddrows (env, cut_prob, 0, cur_numrows+2*cur_numcols, 0, NULL, sense, &z, NULL, NULL, NULL, NULL);*/
/*  	if ( status ) {*/
/*   		printf ("%s(%d): Failed to add rows, error code %d\n", __FILE__, __LINE__, status);*/
/*  	}*/
/*  	*/
/*  	status = CPXgetcols (env, cut_prob, &nzcnt, NULL, NULL, NULL, 0, &surplus, 0, cur_numcols);*/
/*  	if ( status ) {*/
/*   		printf ("%s(%d): Failed to 'get cols', error code %d\n", __FILE__, __LINE__, status);*/
/*  	}*/
/*  	*/
/*  	cmatspace = -surplus;*/
/*  	int *cmatbeg = (int *) malloc ( (cur_numcols+1)*sizeof(int) );*/
/*	int *cmatind = (int *) malloc ( cmatspace*sizeof(int) );*/
/*	double *cmatval = (double *) malloc ( cmatspace*sizeof(double) );*/
/* */
/* 	status = CPXgetcols (env, cut_prob, &nzcnt, cmatbeg, cmatind, cmatval, cmatspace, &surplus, 0, cur_numcols);*/
/*  	if ( status ) {*/
/*   		printf ("%s(%d): Failed to actually get cols, error code %d\n", __FILE__, __LINE__, status);*/
/*  	}*/
/*  	*/
/*  	for(i=0;i<cmatspace;i++) cmatind[i] = cmatind[i]+ cur_numrows + 2*cur_numcols;*/
/*  	*/
/*  	status = CPXaddcols (env, cut_prob, cur_numcols+1, nzcnt, NULL, cmatbeg, cmatind, cmatval, lb, ub, NULL);*/
/*  	if ( status ) {*/
/*   		printf ("%s(%d): Failed to add cols, error code %d\n", __FILE__, __LINE__, status);*/
/*  	}*/
/*  	*/
/*  	for(i=0;i<cur_numcols;i++)*/
/*  	{*/
/*  		rmatind[2*i] = i;*/
/*  		rmatind[2*i+1] = i + cur_numcols + 1;*/
/*  		rmatval[2*i] = -1.;*/
/*  		rmatval[2*i+1] = -1.;*/
/*  	}*/
/*  	*/
/*  	status = CPXaddrows (env, cut_prob, 0, cur_numcols, 2*cur_numcols, NULL, NULL, rmatbeg, rmatind, rmatval, NULL, NULL);*/
/*  	if ( status ) {*/
/*   		printf ("%s(%d): Failed to add rows, error code %d\n", __FILE__, __LINE__, status);*/
/*  	}*/
/*  	*/
/*  	CPXchgobjsen (env, cut_prob, CPX_MIN);*/
/*  	*/
/*  	free(lb); */
/*  	free(ub);*/
/*  	free(sense);*/
/*  	free(rmatbeg);*/
/*  	free(rmatind);*/
/*  	free(rmatval);*/
/*  	free(rhs);*/
/*  	free(zeros);*/
/*  	free(indices);*/
/*  	free(col_list);*/
/*  	free(cmatbeg);*/
/*  	free(cmatind);*/
/*  	free(cmatval);*/
/*  	free(lu);*/
/*  	*/
/*  	return(cut_prob);*/
/*}*/

CPXLPptr create_cut_prob(CPXCENVptr env, CPXLPptr prob)
{
	int status = 0, i = 0, z = 0, nzcnt = 0, surplus = 0, cmatspace = 0;
	
	double *lb = (double *) malloc( (cur_numrows+2*cur_numcols)*sizeof(double) );
	double *ub = (double *) malloc( (cur_numrows+2*cur_numcols)*sizeof(double) );
	char *sense = (char *) malloc ( (cur_numrows+2*cur_numcols+4)*sizeof(char) );
/*	char *sense2 = (char *) malloc ( (cur_numrows+2*cur_numcols)*sizeof(char) );*/
	char *lu = (char *) malloc ( (cur_numcols+1)*sizeof(char) );
	int *rmatbeg = (int *) malloc ( 2*cur_numcols*sizeof(int) );
	int *rmatind = (int *) malloc ( 4*cur_numcols*sizeof(int) );
	double *rmatval = (double *) malloc ( 4*cur_numcols*sizeof(double) );
	double *rhs = (double *) malloc ( cur_numrows*sizeof(double) );
	double *zeros = (double *) calloc( cur_numrows, sizeof(double) );
	int *indices = (int *) malloc (fmax(cur_numrows,cur_numcols+1)*sizeof(int) );
	int *col_list = (int *) malloc (cur_numrows*sizeof(int) );
	
	CPXLPptr cut_prob = CPXcloneprob (env, prob, &status);
	if ( status ) {
   		printf ("%s(%d): Failed to copy prob, error code %d\n", __FILE__, __LINE__, status);
  	}
  	
  	status = ( CPXgetub (env, cut_prob, ub, 0, cur_numcols-1) || CPXgetlb (env, cut_prob, lb, 0, cur_numcols-1) );
  	if ( status ) {
   		printf ("%s(%d): Failed to get variable bounds, error code %d\n", __FILE__, __LINE__, status);
  	}
  	
  	status = CPXgetrhs (env, cut_prob, rhs, 0, cur_numrows-1);
  	if ( status ) {
   		printf ("%s(%d): Failed to get rhs, error code %d\n", __FILE__, __LINE__, status);
  	}
  	
  	for(i=0;i<cur_numrows;i++)
  	{
  		indices[i] = i;
  		col_list[i] = cur_numcols;
  	}
  	status = CPXchgrhs (env, cut_prob, cur_numrows, indices, zeros);
  	if ( status ) {
   		printf ("%s(%d): Failed to chg rhs, error code %d\n", __FILE__, __LINE__, status);
  	}
  	status = CPXgetsense (env, cut_prob, sense, 0, cur_numrows-1);
  	
  	for(i=0;i<2*cur_numcols;i++)
  	{
/*  		if(i<cur_numrows && sense[i] != 'E') sense[i] = 'G';*/
  		sense[i] = 'G';
  		rmatbeg[i] = 2*i;
  		if(i < cur_numcols) rmatind[2*i] = i;
  		else rmatind[2*i] = i - cur_numcols;
  		rmatind[2*i+1] = cur_numcols;
  		if(i < cur_numcols)
  		{
  			rmatval[2*i] = 1.;
  			rmatval[2*i+1] = ub[i];
  		}
  		else
  		{
  			rmatval[2*i] = -1.;
  			rmatval[2*i+1] = -lb[i-cur_numcols];
  		}
  	}
  	if(cur_numrows > 2*cur_numcols) 
  	{
  		for(i=2*cur_numcols;i<cur_numrows;i++) 
  		{
/*  			if(sense[i] != 'E') sense[i] = 'G';*/
  			sense[i] = 'G';
  		}
/*		for(i=2*cur_numcols;i<cur_numrows;i++) sense[i] = 'L';*/
  	}
  	
  	status = CPXchgsense (env, cut_prob, cur_numrows, indices, sense);
  	
/*  	status = CPXwriteprob (env, cut_prob, "myprob1.lp", "LP");*/
/*  	status = CPXwriteprob (env, prob, "myprob2.lp", "LP");*/
/*	exit(0);*/
  	
  	status = CPXaddrows (env, cut_prob, 1, 2*cur_numcols, 4*cur_numcols, NULL, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
  	if ( status ) {
   		printf ("%s(%d): Failed to add rows, error code %d\n", __FILE__, __LINE__, status);
  	}
  	
  	status = CPXchgcoeflist (env, cut_prob, cur_numrows, indices, col_list, rhs);
  	if ( status ) {
   		printf ("%s(%d): Failed to change col coefs, error code %d\n", __FILE__, __LINE__, status);
  	}
  	
/*  	status = CPXwriteprob (env, cut_prob, "myprob1.lp", "LP");*/
/*  	status = CPXwriteprob (env, prob, "myprob2.lp", "LP");*/
/*	exit(0);*/
  	
  	for(i=0;i<cur_numrows+2*cur_numcols;i++)
  	{
  		lb[i] = -CPX_INFBOUND;
  		ub[i] = CPX_INFBOUND;
  	}
/*  	lb[cur_numcols] = 0.;*/
/*  	ub[cur_numcols] = CPX_INFBOUND;*/
  	if(cur_numcols + 1 > cur_numrows) for(i=cur_numrows;i<cur_numcols+1;i++) indices[i] = i;
  	
  	for(i=0;i<cur_numcols+1;i++) lu[i] = 'L';
  	
  	status = CPXchgbds (env, cut_prob, cur_numcols+1, indices, lu, lb);
  	if ( status ) {
   		printf ("%s(%d): Failed to change lbs, error code %d\n", __FILE__, __LINE__, status);
  	}
  	
  	for(i=0;i<cur_numcols+1;i++) lu[i] = 'U';
  	
  	status = CPXchgbds (env, cut_prob, cur_numcols+1, indices, lu, ub);
  	if ( status ) {
   		printf ("%s(%d): Failed to change ubs, error code %d\n", __FILE__, __LINE__, status);
  	}
  	
  	for(i=2*cur_numcols;i<cur_numrows+2*cur_numcols+4;i++) sense[i] = 'G';
/*	for(i=2*cur_numcols;i<cur_numrows+2*cur_numcols;i++) sense[i] = 'L';*/


/*	for(i=0;i<cur_numrows+2*cur_numcols;i++) if(sense[i] != 'G') printf("sense_%d: %c\n",i,sense[i]);*/
/*	printf("%d\n", cur_numrows+2*cur_numcols);*/
  	
  	double zero[1] = {0.};
  	int *zeros_ = calloc(cur_numrows+2*cur_numcols, sizeof(int));
  	
/*  	status = CPXwriteprob (env, cut_prob, "myprob1.lp", "LP");*/
  	 
  	status = CPXaddrows (env, cut_prob, 0, cur_numrows+2*cur_numcols, 0, NULL, sense, zeros_, &z, zero, NULL, NULL);
  	if ( status ) {
   		printf ("%s(%d): Failed to add rows, error code %d\n", __FILE__, __LINE__, status);
  	}
  	
  	status = CPXgetcols (env, cut_prob, &nzcnt, NULL, NULL, NULL, 0, &surplus, 0, cur_numcols);
/*  	if ( status ) {*/
/*   		printf ("%s(%d): Failed to 'get cols', error code %d\n", __FILE__, __LINE__, status);*/
/*  	}*/
  	
  	cmatspace = -surplus;
  	int *cmatbeg = (int *) malloc ( (cur_numcols+1)*sizeof(int) );
	int *cmatind = (int *) malloc ( cmatspace*sizeof(int) );
	double *cmatval = (double *) malloc ( cmatspace*sizeof(double) );
 
 	status = CPXgetcols (env, cut_prob, &nzcnt, cmatbeg, cmatind, cmatval, cmatspace, &surplus, 0, cur_numcols);
  	if ( status ) {
   		printf ("%s(%d): Failed to actually get cols, error code %d\n", __FILE__, __LINE__, status);
  	}
  	
  	for(i=0;i<cmatspace;i++) cmatind[i] = cmatind[i]+ cur_numrows + 2*cur_numcols;
  	
  	status = CPXaddcols (env, cut_prob, cur_numcols+1, nzcnt, NULL, cmatbeg, cmatind, cmatval, lb, ub, NULL);
  	if ( status ) {
   		printf ("%s(%d): Failed to add cols, error code %d\n", __FILE__, __LINE__, status);
  	}
  	
  	for(i=0;i<cur_numcols;i++)
  	{
  		rmatind[2*i] = i;
  		rmatind[2*i+1] = i + cur_numcols + 1;
  		rmatval[2*i] = -1.;
  		rmatval[2*i+1] = -1.;
  	}
  	
  	status = CPXaddrows (env, cut_prob, 0, cur_numcols, 2*cur_numcols, NULL, NULL, rmatbeg, rmatind, rmatval, NULL, NULL);
  	if ( status ) {
   		printf ("%s(%d): Failed to add rows, error code %d\n", __FILE__, __LINE__, status);
  	}
  	
  	CPXchgobjsen (env, cut_prob, CPX_MIN);
/*	CPXchgobjsen (env, cut_prob, CPX_MAX);*/
  	
  	free(lb); 
  	free(ub);
  	free(sense);
  	free(rmatbeg);
  	free(rmatind);
  	free(rmatval);
  	free(rhs);
  	free(zeros);
  	free(indices);
  	free(col_list);
  	free(cmatbeg);
  	free(cmatind);
  	free(cmatval);
  	free(lu);
  	free(zeros_);
  	
  	return(cut_prob);
}

/*********************************************************************************************************************** 

	This function is used for reducing bounds after processing each node but before branching.
	
***********************************************************************************************************************/

int reduce_bds_before_branching(CPXENVptr env, CPXLPptr prob, double *bds_b1, double *bds_b2, int reducing_size, int branch_index, double branch_val, double slope)
{
	int status = 0;
	int i,j,k;
	double NW_objectives[2] = {0.,0.};
	double SE_objectives[2] = {0.,0.};
	double WS_objectives[2] = {0.,0.};
	
	int insert_check = 0;
	int PSA_reduce_val = 0;
	double temp_x,temp_y;
	int fix_it2 = 0;
	int num_bounds_reduced = 0;
	int extra = 0;
	int retval = 0;
	int counter = 0;
	
	int *indices = NULL;
    	indices = (int *) malloc (max(cur_numcols+2,cur_numrows+2)*sizeof(int));
	for(i=0;i<max(cur_numcols+2,cur_numrows+2);i++) indices[i] = i;
	
	double *x = (double *) malloc (cur_numcols*sizeof(double));
/*	printf("allocating\n");*/
	int *basis_col = (int *) malloc (cur_numcols*sizeof(int));
	int *basis_row = (int *) malloc (cur_numrows*sizeof(int));
	
	int lpstat;
	int lowest_index = 0;
	
	clock_t bd_red_start_time, bd_red_intermediate_time;
	double time_so_far = 0.;
	
	CPXLPptr reduced_lp = CPXcloneprob (env, prob, &status);
	if ( status ) {
   		printf ("%s(%d): CPXchgbds, Failed to copy prob, error code %d\n", __FILE__, __LINE__, status);
  	}
	if(reducing_size)
	{
		int ind[4] = {obj1_index,obj1_index,obj2_index,obj2_index};
		char lu[4] = {'L','U','L','U'};
		double bds[4] = {bds_b1[2*total_num_integer],bds_b1[2*total_num_integer+1],bds_b1[2*total_num_integer+2],bds_b1[2*total_num_integer+3]};
		status = CPXchgbds (env, reduced_lp, 4, ind, lu, bds);
		if ( status ) {
	   		printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
	  	}
	  	extra = 2;
	}
	CPXLPptr reduced_lp2 = CPXcloneprob (env, reduced_lp, &status);
	if ( status ) {
   		printf ("%s(%d): CPXchgbds, Failed to copy prob, error code %d\n", __FILE__, __LINE__, status);
  	}
  	CPXLPptr temp_lp = CPXcloneprob (env, reduced_lp, &status);
	if ( status ) {
   		printf ("%s(%d): CPXchgbds, Failed to copy prob, error code %d\n", __FILE__, __LINE__, status);
  	}
  	
  	char lu2[1] = {'U'};
  	double bd = (int) branch_val;
  	
  	START:
/*  	printf("changing ub of %d to %lf\n",branch_index,bd);*/

	bd_red_start_time = clock();

  	status = CPXchgbds (env, reduced_lp, 1, &branch_index, lu2, &bd);
	if ( status ) {
	   	printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
	}
/*	status = CPXwriteprob (env, reduced_lp, "myprob1.lp", "LP");*/
	status = CPXlpopt (env, reduced_lp);
	lpstat = CPXgetstat (env, reduced_lp);
/*  	printf("the status of the solve: %d\n",lpstat);*/
  	if(lpstat == 3)
  	{
/*  		printf("changing ub of x%d to %lf resulted in an infeasible problem.\n",branch_index,bd);*/
 		CPXLPptr reduced_lp = CPXcloneprob (env, temp_lp, &status); 		
 		j = 0;
 		for(i=0;i<total_num_integer;i++)
 		{
 			if(integer_indices[i] == branch_index) 
 			{
 				j = i;
 				break;
 			}
 		}
 		bds_b1[2*j] = bd+1;
 		bds_b2[2*j] = bd+1;
 		lu2[0] = 'L';
 		status = CPXchgbds (env, reduced_lp, 1, &branch_index, lu2, &bds_b1[2*j]);
 		status = CPXchgbds (env, reduced_lp2, 1, &branch_index, lu2, &bds_b1[2*j]);
 		status = CPXchgbds (env, temp_lp, 1, &branch_index, lu2, &bds_b1[2*j]);
 		lu2[0] = 'U';
 		counter = 0;
 		for(i=lowest_index;i<total_num_integer;i++)
 		{
 			if(bds_b1[2*i] != bds_b1[2*i+1])
 			{
/* 				printf("integer index: %d, branch index: %d\n",integer_indices[i],branch_index);*/
 				if(integer_indices[i] == branch_index) goto END;
 				branch_index = integer_indices[i];
 				branch_val = fmax(bds_b1[2*i],bds_b1[2*i+1])-.5;
 				bd = (int) branch_val;
 				lowest_index = i;
 				counter++;
 			}
 		}
 		if(!counter) goto END;
/* 		printf("changing branch variable to x%d.\n",branch_index);*/
  		goto START;
  	}
  	
  	lu2[0] = 'L';
  	bd = bd + 1;
/*  	printf("changing lb of %d to %lf\n",branch_index,bd);*/
  	status = CPXchgbds (env, reduced_lp2, 1, &branch_index, lu2, &bd);
	if ( status ) {
	   	printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
	}
/*	status = CPXwriteprob (env, presolve2_lp, "myprob1.lp", "LP");*/
	status = CPXlpopt (env, reduced_lp2);
	lpstat = CPXgetstat (env, reduced_lp2);
/*  	printf("the status of the solve: %d\n",lpstat);*/
  	if(lpstat == 3)
  	{
/*  		printf("changing lb of x%d to %lf resulted in an infeasible problem.\n",branch_index,bd);*/
 		CPXLPptr reduced_lp2 = CPXcloneprob (env, temp_lp, &status);
 		j = 0;
 		for(i=0;i<total_num_integer;i++)
 		{
 			if(integer_indices[i] == branch_index) 
 			{
 				j = i;
 				break;
 			}
 		}
 		bds_b1[2*j+1] = bd-1;
 		bds_b2[2*j+1] = bd-1;
 		lu2[0] = 'U';
 		status = CPXchgbds (env, reduced_lp, 1, &branch_index, lu2, &bds_b1[2*j+1]);
 		status = CPXchgbds (env, reduced_lp2, 1, &branch_index, lu2, &bds_b1[2*j+1]);
 		status = CPXchgbds (env, temp_lp, 1, &branch_index, lu2, &bds_b1[2*j+1]);
 		lu2[0] = 'L';
 		counter = 0;
 		for(i=lowest_index;i<total_num_integer;i++)
 		{
 			if(bds_b1[2*i] != bds_b1[2*i+1])
 			{
/* 				printf("integer index: %d, branch index: %d\n",integer_indices[i],branch_index);*/
 				if(integer_indices[i] == branch_index) goto END;
 				branch_index = integer_indices[i];
 				branch_val = fmax(bds_b1[2*i],bds_b1[2*i+1])-.5;
 				bd = (int) branch_val;
 				lowest_index = i;
 				counter++;
 			}
 		}
 		if(!counter) goto END;
/* 		printf("changing branch variable to x%d.\n",branch_index);*/
  		goto START;
  	}
  	
/*  	printf("making it here means that x%d will definitely be the branch variable\n",branch_index);*/
  	
  	j = 0;
 	for(i=0;i<total_num_integer;i++)
 	{
 		if(integer_indices[i] == branch_index) 
 		{
 			j = i;
 			break;
 		}
 	}
 	bd = (int) branch_val;
 	bds_b1[2*j+1] = bd;
 	bds_b2[2*j] = bd + 1;
 	lu2[0] = 'U';
 	status = CPXchgbds (env, reduced_lp, 1, &branch_index, lu2, &bds_b1[2*j+1]);
 	lu2[0] = 'L';
 	status = CPXchgbds (env, reduced_lp2, 1, &branch_index, lu2, &bds_b2[2*j]);
  	
  	char sym[1] = {'U'};
  	
  	double ran = (double) rand() / ( (double) RAND_MAX);
	int ten_percent = (int) (total_num_integer/10.);
	int starting_index = (int) (ran*((double) total_num_integer - (double) ten_percent));
/*	printf("ran: %lf, ten_percent: %d, starting_index: %d\n",ran,ten_percent,starting_index);*/

	if(!limit_bd_reduction)
	{
		starting_index = 0;
		ten_percent = total_num_integer;
	}
	else if(stopped_early_last_time) 
	{
		starting_index = index_we_left_off_at;
		ten_percent = total_num_integer - starting_index;
	}
	
	for(i=starting_index;i<starting_index+ten_percent;i++)
	{
		bd_red_intermediate_time = clock();	
		time_so_far = (double)(bd_red_intermediate_time - bd_red_start_time) / CLOCKS_PER_SEC;
		if(time_so_far > max_time_bd_red_b4_branching)
		{
/*			printf("Exceeded time limit of %lf when reducing variable bounds before branching. Stopping!\n", max_time_bd_red_b4_branching);*/
			stopped_early_last_time = 1;
			index_we_left_off_at = i;
			
			free(indices);
		    	free(x);
/*		    	printf("freeing\n");*/
			free(basis_col);
			free(basis_row);
			CPXfreeprob(env, &reduced_lp);
			CPXfreeprob(env, &reduced_lp2);
			CPXfreeprob(env, &temp_lp);
			
			return 1;
		} 
	
/*		printf("i: %d\n",i);*/
		fix_it2 = 0;
		if(bds_b1[2*i] != bds_b1[2*i+1] && bds_b1[2*i+1] - bds_b1[2*i] < 10.)
		{
			REPEAT_IT:
			sym[0] = 'U';
/*			printf("changing upper bound of x%d to %lf\n",integer_indices[i], bds_b1[2*i]);*/
    			status = CPXchgbds (env, reduced_lp, 1, &integer_indices[i], sym, &bds_b1[2*i]);
    			
/*    			printf("_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-\n");*/
/*			PSA_all(env,reduced_lp);*/
/*			printf("_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-\n");*/
			
/*			printf("slope: %lf\n",slope);*/
    			
			chg_coefs(env,reduced_lp,indices,slope);
			status = CPXlpopt (env, reduced_lp);
			lpstat = CPXgetstat (env, reduced_lp);
			if(lpstat == 3) 
	  		{
	  			fix_it2 = 1;
/*	  			printf("fixing by infeasibility\n");*/
	  			bds_reduced_by_infeasibility++;
	  			goto FIXING;
			}
			if(!only_infeasibility)
			{
				status = CPXgetx (env, reduced_lp, WS_objectives, obj1_index, obj2_index);
	  			if(status)
	  			{
					printf("Failed to get x-values from CPLEX. Status: %d Line: %d\n",status,__LINE__);
					exit(0);
	  			}
/*	  			printf("plot(%lf,%lf,'o');\n",WS_objectives[0],WS_objectives[1]);*/
	  			insert_check = mock_insert(1,WS_objectives[0],WS_objectives[1],0,0,0,&tree);
	  			if(!insert_check)
	  			{
/*	  				printf("WS pt dominated when fixing x%d to its lb\n",i);*/
/*			  		printf("plot(%lf,%lf,'o');\n",WS_objectives[0],WS_objectives[1]);*/
	  				chg_coefs(env,reduced_lp,indices,-.000000001);
/*	  				status = CPXwriteprob (env, reduced_lp, "myprob2.lp", "LP");*/
					status = CPXlpopt (env, reduced_lp);
/*					printf("(%d) status: %d\n",__LINE__,status);*/
/*					lpstat = CPXgetstat (env, reduced_lp);*/
/*					printf("(%d) lpstat: %d\n",__LINE__,lpstat);*/
					status = CPXgetx (env, reduced_lp, SE_objectives, obj1_index, obj2_index);
			  		if(status)
					{
/*						printf("Failed to get x-values from CPLEX. Line: %d\n",__LINE__);*/
			  		}
			  		insert_check = mock_insert(1,SE_objectives[0],SE_objectives[1],0,0,0,&tree);
			  		if(!insert_check)
	  				{
/*		  				printf("SE pt dominated when fixing x%d to its lb\n",i);*/
/*			  			printf("plot(%lf,%lf,'o');\n",SE_objectives[0],SE_objectives[1]);*/
			  			chg_coefs(env,reduced_lp,indices,-10000000.);
/*			  			status = CPXwriteprob (env, reduced_lp, "myprob3.lp", "LP");*/
						status = CPXlpopt (env, reduced_lp);
/*						printf("(%d) status: %d\n",__LINE__,status);*/
/*						lpstat = CPXgetstat (env, reduced_lp);*/
/*						printf("(%d) lpstat: %d\n",__LINE__,lpstat);*/
						status = CPXgetx (env, reduced_lp, NW_objectives, obj1_index, obj2_index);
			  			if(status)
			  			{
/*							printf("Failed to get x-values from CPLEX. Line: %d\n",__LINE__);*/
			  			}
			  			insert_check = mock_insert(1,NW_objectives[0],NW_objectives[1],0,0,0,&tree);
		  				if(!insert_check)
		  				{
/*		  					printf("NW pt also dominated when fixing x%d to its lb\n",i);*/
/*				  			printf("plot(%lf,%lf,'o');\n",NW_objectives[0],NW_objectives[1]);*/
				  			insert_check = mock_insert(1,NW_objectives[0],SE_objectives[1],0,0,0,&tree);
				  			if(!insert_check)
				  			{
/*				  				printf("local ideal pt dominated when fixing x%d to its lb\n",i);*/
/*				  				printf("plot(%lf,%lf,'go');\n",NW_objectives[0],SE_objectives[1]);*/
				  				fix_it2 = 1;
				  				bds_reduced_by_dominated_ideal_pt++;
				  			}
				  			else
				  			{
				  				insert_check = mock_insert(2,NW_objectives[0],(slope)*(NW_objectives[0]-WS_objectives[0])+
				  					WS_objectives[1], (1./slope)*(SE_objectives[1]-WS_objectives[1])+WS_objectives[0],
				  					SE_objectives[1],(slope),&tree);
								if(!insert_check)
								{
/*									printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",NW_objectives[0],*/
/*										(1./slope)*(SE_objectives[1]-WS_objectives[1])+WS_objectives[0],*/
/*										(slope)*(NW_objectives[0]-WS_objectives[0])+WS_objectives[1],SE_objectives[1]);*/
/*									printf("line: %d\n",__LINE__);*/
									fix_it2 = 1;
				  					bds_reduced_by_dominated_ideal_segment++;
								}
				  				else if(check_bound)
				  				{
									status = CPXgetx (env, reduced_lp, x, 0, cur_numcols-1);
									status = CPXgetbase (env, reduced_lp, basis_col, basis_row);
									PSA_reduce_val = PSA_reduce_right(env, reduced_lp, x, basis_col, basis_row, indices, -1);
/*									printf("PSA reduce val: %d\n",PSA_reduce_val);*/
									if(PSA_reduce_val == 2) 
									{
										fix_it2 = 1;
										bds_reduced_by_dominated_lb++;
									}
								}
				  			}
				  		}
			  		}
	  			}
  			}
  			FIXING:
  			sym[0] = 'U';
/*  			printf("changing upper bound of x%d back to %lf\n",integer_indices[i], bds_b1[2*i+1]);*/
    			status = CPXchgbds (env, reduced_lp, 1, &integer_indices[i], sym, &bds_b1[2*i+1]);
    			if(fix_it2 == 1)
    			{
/*    				printf("Changing lb of x%d to %d!\n",integer_indices[i],(int)(bds_b1[2*i]+1));*/
    				sym[0] = 'L';
    				double bound[1] = {bds_b1[2*i]+1.};
    				status = CPXchgbds (env, reduced_lp, 1, &integer_indices[i], sym, bound);
    				bds_b1[2*i]++;
    				num_bounds_reduced++;
    				bds_reduced++;
    				if(bds_b1[2*i] + 1 < bds_b1[2*i+1])
    				{
/*    					printf("x%d has had its lower bound increased by 1. Attempting to increase further\n",i);*/
    					fix_it2 = 0;
    					goto REPEAT_IT;
    				}
    			}
    			if(bds_b1[2*i] != bds_b1[2*i+1] && (!fix_it2 || xctype[i] != 'B') && bds_b1[2*i+1] - bds_b1[2*i] < 10.)
			{
				REPEAT_IT2:
				fix_it2 = 0;
				sym[0] = 'L';
/*				printf("changing lower bound of x%d to %lf\n",integer_indices[i], bds_b1[2*i+1]);*/
	    			status = CPXchgbds (env, reduced_lp, 1, &integer_indices[i], sym, &bds_b1[2*i+1]);
	    			
				chg_coefs(env,reduced_lp,indices,slope);
				status = CPXlpopt (env, reduced_lp);
				lpstat = CPXgetstat (env, reduced_lp);
/*				printf("status: %d\n",lpstat);*/
				if(lpstat == 3) 
	  			{
	  				fix_it2 = 1;
/*	  				printf("fixing by infeasibility\n");*/
	  				bds_reduced_by_infeasibility++;
	  				goto FIXING2;
	  			}
	  			if(!only_infeasibility)
	  			{
				status = CPXgetx (env, reduced_lp, WS_objectives, obj1_index, obj2_index);
	  			if(status)
	  			{
					printf("Failed to get x-values from CPLEX. Status: %d Line: %d\n",status,__LINE__);
					exit(0);
	  			}
/*	  			printf("plot(%lf,%lf,'o');\n",WS_objectives[0],WS_objectives[1]);*/
	  			insert_check = mock_insert(1,WS_objectives[0],WS_objectives[1],0,0,0,&tree);
	  			if(!insert_check)
	  			{
/*	  				printf("WS pt dominated when fixing x%d to its lb\n",i);*/
/*	  				printf("plot(%lf,%lf,'o');\n",WS_objectives[0],WS_objectives[1]);*/
	  				chg_coefs(env,reduced_lp,indices,-.000000001);
					status = CPXlpopt (env, reduced_lp);
					status = CPXgetx (env, reduced_lp, SE_objectives, obj1_index, obj2_index);
			  		if(status)
					{
/*						printf("Failed to get x-values from CPLEX. Line: %d\n",__LINE__);*/
			  		}
			  		insert_check = mock_insert(1,SE_objectives[0],SE_objectives[1],0,0,0,&tree);
			  		if(!insert_check)
	  				{
/*		  				printf("SE pt dominated when fixing x%d to its ub\n",i);*/
/*			  			printf("plot(%lf,%lf,'o');\n",SE_objectives[0],SE_objectives[1]);*/
			  			chg_coefs(env,reduced_lp,indices,-10000000.);
						status = CPXlpopt (env, reduced_lp);
						status = CPXgetx (env, reduced_lp, NW_objectives, obj1_index, obj2_index);
			  			if(status)
			  			{
/*							printf("Failed to get x-values from CPLEX. Line: %d\n",__LINE__);*/
			  			}
			  			insert_check = mock_insert(1,NW_objectives[0],NW_objectives[1],0,0,0,&tree);
		  				if(!insert_check)
		  				{
/*		  					printf("NW pt also dominated when fixing x%d to its ub\n",i);*/
/*				  			printf("plot(%lf,%lf,'o');\n",NW_objectives[0],NW_objectives[1]);*/
				  			insert_check = mock_insert(1,NW_objectives[0],SE_objectives[1],0,0,0,&tree);
				  			if(!insert_check)
				  			{
/*				  				printf("local ideal pt dominated when fixing x%d to its ub\n",i);*/
/*				  				printf("plot(%lf,%lf,'go');\n",NW_objectives[0],SE_objectives[1]);*/
				  				fix_it2 = 1;
				  				bds_reduced_by_dominated_ideal_pt++;
				  			}
				  			else
				  			{
								insert_check = mock_insert(2,NW_objectives[0],
									(slope)*(NW_objectives[0]-WS_objectives[0])+WS_objectives[1],
									(1./slope)*(SE_objectives[1]-WS_objectives[1])+WS_objectives[0],
									SE_objectives[1],(slope),&tree);
								if(!insert_check)
								{
/*									printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",NW_objectives[0],*/
/*										(1./slope)*(SE_objectives[1]-WS_objectives[1])+WS_objectives[0],*/
/*										(slope)*(NW_objectives[0]-WS_objectives[0])+WS_objectives[1],SE_objectives[1]);*/
/*									printf("line: %d\n",__LINE__);*/
									fix_it2 = 1;
				  					bds_reduced_by_dominated_ideal_segment++;
								}
				  				else if(check_bound)
				  				{
									status = CPXgetx (env, reduced_lp, x, 0, cur_numcols-1);
									status = CPXgetbase (env, reduced_lp, basis_col, basis_row);
									PSA_reduce_val = PSA_reduce_right(env, reduced_lp, x, basis_col, basis_row, indices, -1);
/*									printf("PSA reduce val: %d\n",PSA_reduce_val);*/
									if(PSA_reduce_val == 2) 
									{
										fix_it2 = 1;
										bds_reduced_by_dominated_lb++;
									}
								}
				  			}
				  		}
				  	}
	  			}
	  			}
	  			FIXING2:
	  			sym[0] = 'L';
/*	  			printf("changing lower bound of x%d back to %lf\n",integer_indices[i], bds_b1[2*i]);*/
	    			status = CPXchgbds (env, reduced_lp, 1, &integer_indices[i], sym, &bds_b1[2*i]);
	    			if(fix_it2 == 1)
	    			{
/*	    				printf("Changing ub of x%d to %d!\n",integer_indices[i],(int)(bds_b1[2*i+1]-1));*/
	    				sym[0] = 'U';
	    				double bound[1] = {bds_b1[2*i+1]-1.};
	    				status = CPXchgbds (env, reduced_lp, 1, &integer_indices[i], sym, bound);
	    				bds_b1[2*i+1]--;
	    				num_bounds_reduced++;
	    				bds_reduced++;
	    				if(bds_b1[2*i] + 1 < bds_b1[2*i+1])
	    				{
/*	    					printf("x%d has had its upper bound decreased by 1. Attempting to decrease further\n",i);*/
	    					fix_it2 = 0;
	    					goto REPEAT_IT2;
	    				}
	    			}
	    		}	
    		}
    	}
/*    	printf("done attempting extra fixings for branch1\n");*/
	int count = num_bounds_reduced;
/*	printf("number of bounds reduced for branch1: %d\n",num_bounds_reduced); 	*/
    	
    	for(i=starting_index;i<starting_index+ten_percent;i++)
	{
		bd_red_intermediate_time = clock();	
		time_so_far = (double)(bd_red_intermediate_time - bd_red_start_time) / CLOCKS_PER_SEC;
		if(time_so_far > max_time_bd_red_b4_branching)
		{
/*			printf("Exceeded time limit of %lf when reducing variable bounds before branching. Stopping!\n", max_time_bd_red_b4_branching);*/
			stopped_early_last_time = 1;
			index_we_left_off_at = i;
			
			free(indices);
		    	free(x);
/*		    	printf("freeing\n");*/
			free(basis_col);
			free(basis_row);
			CPXfreeprob(env, &reduced_lp);
			CPXfreeprob(env, &reduced_lp2);
			CPXfreeprob(env, &temp_lp);
			
			return 1;
		} 
		fix_it2 = 0;
		if(bds_b2[2*i] != bds_b2[2*i+1] && bds_b2[2*i+1] - bds_b2[2*i] < 10.)
		{
			REPEAT_IT3:
			sym[0] = 'U';
/*			printf("changing upper bound of x%d to %lf\n",integer_indices[i], bds_b2[2*i]);*/
    			status = CPXchgbds (env, reduced_lp2, 1, &integer_indices[i], sym, &bds_b2[2*i]);
    			
			chg_coefs(env,reduced_lp2,indices,slope);
			status = CPXlpopt (env, reduced_lp2);
			lpstat = CPXgetstat (env, reduced_lp2);
/*			printf("the status of the solve: %d\n",lpstat);*/
			if(lpstat == 3) 
	  		{
	  			fix_it2 = 1;
/*	  			printf("fixing by infeasibility\n");*/
	  			bds_reduced_by_infeasibility++;
	  			goto FIXING3;
			}
			if(!only_infeasibility)
			{
			status = CPXgetx (env, reduced_lp2, WS_objectives, obj1_index, obj2_index);
  			if(status)
  			{
				printf("Failed to get x-values from CPLEX. Status: %d Line: %d\n",status,__LINE__);
				exit(0);
  			}
/*  			printf("plot(%lf,%lf,'go');\n",WS_objectives[0],WS_objectives[1]);*/
  			insert_check = mock_insert(1,WS_objectives[0],WS_objectives[1],0,0,0,&tree);
  			if(!insert_check)
  			{
/*  				printf("WS pt dominated when fixing x%d to its lb\n",i);*/
/*  				printf("plot(%lf,%lf,'go');\n",WS_objectives[0],WS_objectives[1]);*/
  				chg_coefs(env,reduced_lp2,indices,-.000000001);
				status = CPXlpopt (env, reduced_lp2);
				status = CPXgetx (env, reduced_lp2, SE_objectives, obj1_index, obj2_index);
		  		if(status)
				{
					printf("Failed to get x-values from CPLEX. Line: %d\n",__LINE__);
		  		}
		  		insert_check = mock_insert(1,SE_objectives[0],SE_objectives[1],0,0,0,&tree);
		  		if(!insert_check)
  				{
/*	  				printf("SE pt dominated when fixing x%d to its lb\n",i);*/
/*		  			printf("plot(%lf,%lf,'o');\n",SE_objectives[0],SE_objectives[1]);*/
		  			chg_coefs(env,reduced_lp2,indices,-10000000.);
					status = CPXlpopt (env, reduced_lp2);
					status = CPXgetx (env, reduced_lp2, NW_objectives, obj1_index, obj2_index);
		  			if(status)
		  			{
/*						printf("Failed to get x-values from CPLEX. Line: %d\n",__LINE__);*/
		  			}
		  			insert_check = mock_insert(1,NW_objectives[0],NW_objectives[1],0,0,0,&tree);
	  				if(!insert_check)
	  				{
/*	  					printf("NW pt also dominated when fixing x%d to its lb\n",i);*/
/*			  			printf("plot(%lf,%lf,'o');\n",NW_objectives[0],NW_objectives[1]);*/
			  			insert_check = mock_insert(1,NW_objectives[0],SE_objectives[1],0,0,0,&tree);
			  			if(!insert_check)
			  			{
/*			  				printf("local ideal pt dominated when fixing x%d to its lb\n",i);*/
/*			  				printf("plot(%lf,%lf,'go');\n",NW_objectives[0],SE_objectives[1]);*/
			  				fix_it2 = 1;
			  				bds_reduced_by_dominated_ideal_pt++;
			  			}
			  			else
			  			{
			  				insert_check = mock_insert(2,NW_objectives[0],(slope)*(NW_objectives[0]-WS_objectives[0])+WS_objectives[1],
								(1./slope)*(SE_objectives[1]-WS_objectives[1])+WS_objectives[0],SE_objectives[1],(slope),&tree);
							if(!insert_check)
							{
/*								printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",NW_objectives[0],*/
/*									(1./slope)*(SE_objectives[1]-WS_objectives[1])+WS_objectives[0],*/
/*									(slope)*(NW_objectives[0]-WS_objectives[0])+WS_objectives[1],SE_objectives[1]);*/
/*								printf("line: %d\n",__LINE__);*/
								fix_it2 = 1;
			  					bds_reduced_by_dominated_ideal_segment++;
							}
			  				else if(check_bound)
			  				{
								status = CPXgetx (env, reduced_lp2, x, 0, cur_numcols-1);
								status = CPXgetbase (env, reduced_lp2, basis_col, basis_row);
								PSA_reduce_val = PSA_reduce_right(env, reduced_lp2, x, basis_col, basis_row, indices, -1);
/*								printf("PSA reduce val: %d\n",PSA_reduce_val);*/
								if(PSA_reduce_val == 2) 
								{
									fix_it2 = 1;
									bds_reduced_by_dominated_lb++;
								}
							}
			  			}
			  		}
			  	}
  			}
  			}
  			FIXING3:
  			sym[0] = 'U';
    			status = CPXchgbds (env, reduced_lp2, 1, &integer_indices[i], sym, &bds_b2[2*i+1]);
/*    			printf("changing upper bound of x%d back to %lf\n",integer_indices[i], bds_b2[2*i+1]);*/
    			if(fix_it2 == 1)
    			{
/*    				printf("Changing lb of x%d to %d!\n",integer_indices[i],(int)(bds_b2[2*i]+1));*/
    				sym[0] = 'L';
    				double bound[1] = {bds_b2[2*i]+1.};
    				status = CPXchgbds (env, reduced_lp2, 1, &integer_indices[i], sym, bound);
    				bds_b2[2*i]++;
    				num_bounds_reduced++;
    				bds_reduced++;
    				if(bds_b2[2*i] + 1 < bds_b2[2*i+1])
    				{
/*    					printf("x%d has had its lower bound increased by 1. Attempting to increase further\n",i);*/
    					fix_it2 = 0;
    					goto REPEAT_IT3;
    				}
    			}
    			if(bds_b2[2*i] != bds_b2[2*i+1] && (!fix_it2 || xctype[i] != 'B') && bds_b2[2*i+1] - bds_b2[2*i] < 10.)
			{
				REPEAT_IT4:
				fix_it2 = 0;
				sym[0] = 'L';
/*				printf("changing lower bound of x%d to %lf\n",integer_indices[i], bds_b2[2*i+1]);*/
	    			status = CPXchgbds (env, reduced_lp2, 1, &integer_indices[i], sym, &bds_b2[2*i+1]);
	    			
				chg_coefs(env,reduced_lp2,indices,slope);
				status = CPXlpopt (env, reduced_lp2);
				lpstat = CPXgetstat (env, reduced_lp2);
				if(lpstat == 3) 
		  		{
		  			fix_it2 = 1;
/*		  			printf("fixing by infeasibility\n");*/
		  			bds_reduced_by_infeasibility++;
		  			goto FIXING4;
				}
				if(!only_infeasibility)
				{
				status = CPXgetx (env, reduced_lp2, WS_objectives, obj1_index, obj2_index);
	  			if(status)
	  			{
					printf("Failed to get x-values from CPLEX. Status: %d Line: %d\n",status,__LINE__);
					exit(0);
	  			}
/*	  			printf("plot(%lf,%lf,'o');\n",WS_objectives[0],WS_objectives[1]);*/
	  			insert_check = mock_insert(1,WS_objectives[0],WS_objectives[1],0,0,0,&tree);
	  			if(!insert_check)
	  			{
/*	  				printf("WS pt dominated when fixing x%d to its ub\n",i);*/
/*	  				printf("plot(%lf,%lf,'o');\n",WS_objectives[0],WS_objectives[1]);*/
	  				chg_coefs(env,reduced_lp2,indices,-.000000001);
					status = CPXlpopt (env, reduced_lp2);
					status = CPXgetx (env, reduced_lp2, SE_objectives, obj1_index, obj2_index);
			  		if(status)
					{
/*						printf("Failed to get x-values from CPLEX. Line: %d\n",__LINE__);*/
			  		}
			  		insert_check = mock_insert(1,SE_objectives[0],SE_objectives[1],0,0,0,&tree);
			  		if(!insert_check)
	  				{
/*		  				printf("SE pt dominated when fixing x%d to its ub\n",i);*/
/*			  			printf("plot(%lf,%lf,'o');\n",SE_objectives[0],SE_objectives[1]);*/
			  			chg_coefs(env,reduced_lp2,indices,-10000000.);
						status = CPXlpopt (env, reduced_lp2);
						status = CPXgetx (env, reduced_lp2, NW_objectives, obj1_index, obj2_index);
			  			if(status)
			  			{
/*							printf("Failed to get x-values from CPLEX. Line: %d\n",__LINE__);*/
			  			}
			  			insert_check = mock_insert(1,NW_objectives[0],NW_objectives[1],0,0,0,&tree);
		  				if(!insert_check)
		  				{
/*		  					printf("NW pt also dominated when fixing x%d to its ub\n",i);*/
/*				  			printf("plot(%lf,%lf,'o');\n",NW_objectives[0],NW_objectives[1]);*/
				  			insert_check = mock_insert(1,NW_objectives[0],SE_objectives[1],0,0,0,&tree);
				  			if(!insert_check)
				  			{
/*				  				printf("local ideal pt dominated when fixing x%d to its ub\n",i);*/
/*				  				printf("plot(%lf,%lf,'go');\n",NW_objectives[0],SE_objectives[1]);*/
				  				fix_it2 = 1;
				  				bds_reduced_by_dominated_ideal_pt++;
				  			}
				  			else
				  			{
				  				insert_check = mock_insert(2,NW_objectives[0],
									(slope)*(NW_objectives[0]-WS_objectives[0])+WS_objectives[1],
									(1./slope)*(SE_objectives[1]-WS_objectives[1])+WS_objectives[0],
									SE_objectives[1],(slope),&tree);
								if(!insert_check)
								{
/*									printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",NW_objectives[0],*/
/*										(1./slope)*(SE_objectives[1]-WS_objectives[1])+WS_objectives[0],*/
/*										(slope)*(NW_objectives[0]-WS_objectives[0])+WS_objectives[1],SE_objectives[1]);*/
/*									printf("line: %d\n",__LINE__);*/
									fix_it2 = 1;
				  					bds_reduced_by_dominated_ideal_segment++;
								}
				  				else if(check_bound)
				  				{
									status = CPXgetx (env, reduced_lp2, x, 0, cur_numcols-1);
									status = CPXgetbase (env, reduced_lp2, basis_col, basis_row);
									PSA_reduce_val = PSA_reduce_right(env, reduced_lp2, x, basis_col, basis_row, indices, -1);
/*									printf("PSA reduce val: %d\n",PSA_reduce_val);*/
									if(PSA_reduce_val == 2) 
									{
										fix_it2 = 1;
										bds_reduced_by_dominated_lb++;
									}
								}
				  			}
				  		}
				  	}
	  			}
	  			}
	  			FIXING4:
	  			sym[0] = 'L';
	    			status = CPXchgbds (env, reduced_lp2, 1, &integer_indices[i], sym, &bds_b2[2*i]);
/*	    			printf("changing lower bound of x%d back to %lf\n",integer_indices[i], bds_b2[2*i]);*/
	    			if(fix_it2 == 1)
	    			{
/*	    				printf("Changing ub of x%d to %d!\n",integer_indices[i],(int)(bds_b2[2*i+1]-1));*/
	    				sym[0] = 'U';
	    				double bound[1] = {bds_b2[2*i+1]-1.};
	    				status = CPXchgbds (env, reduced_lp2, 1, &integer_indices[i], sym, bound);
	    				bds_b2[2*i+1]--;
	    				num_bounds_reduced++;
	    				bds_reduced++;
	    				if(bds_b2[2*i] + 1 < bds_b2[2*i+1])
	    				{
/*	    					printf("x%d has had its upper bound decreased by 1. Attempting to decrease further\n",i);*/
	    					fix_it2 = 0;
	    					goto REPEAT_IT4;
	    				}
	    			}
	    		}	
    		}
    	}
/*    	printf("done attempting extra fixings for branch2\n");*/
/*	printf("number of bounds reduced for branch2: %d\n",num_bounds_reduced-count);*/
    	retval = 1;
    	
    	END:
    	
    	stopped_early_last_time = 0;
	index_we_left_off_at = 0;
    	
    	free(indices);
    	free(x);
/*    	printf("freeing\n");*/
	free(basis_col);
	free(basis_row);
	CPXfreeprob(env, &reduced_lp);
	CPXfreeprob(env, &reduced_lp2);
	CPXfreeprob(env, &temp_lp);
    	
    	return retval;
}

void close_all_env()
{
	int status = 0;
/*	if(env_global!=NULL) status=CPXcloseCPLEX(&env_global);*/
	if(env2!=NULL) status=CPXcloseCPLEX(&env2);
	if(env3!=NULL) status=CPXcloseCPLEX(&env3);
}

/*********************************************************************************************************************** 

	The following function, as well as some of the others that follow, is designed to pass certain
	information from one C file to another. This should probably be done in future versions using
	external global variables.
	
***********************************************************************************************************************/

void pass_the_lps(CPXENVptr env, CPXLPptr prob1, CPXLPptr prob2)
{
	int status = 0;
	lp_1 = CPXcloneprob (env, prob1, &status);
/*	printf("status: %d\n",status);*/
	lp_2 = CPXcloneprob (env, prob2, &status);
/*	printf("status: %d\n",status);*/
        if(!been_passed) 
        {
        	give_env(env);
        	env2 = CPXopenCPLEX(&status);
		if(status)
		{
			printf("failed to open another copy of CPLEX\n");
			exit(0);
		}
		status = CPXsetdblparam (env2, CPX_PARAM_TILIM, time_limit);
		if ( status ) {
		    	printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);
		   	exit(0);
		}
		status = CPXsetintparam (env2, CPX_PARAM_THREADS, 1);
	    	if ( status ) {
	    		printf ("Failure to set threads to 1, error %d.\n",status);
	    		exit(0);
	    	}
	    	env3 = CPXopenCPLEX(&status);
		if(status)
		{
			printf("failed to open another copy of CPLEX\n");
			exit(0);
		}
		status = CPXsetdblparam (env3, CPX_PARAM_TILIM, time_limit);
		if ( status ) {
		    	printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);
		   	exit(0);
		}
		status = CPXsetintparam (env3, CPX_PARAM_THREADS, 1);
	    	if ( status ) {
	    		printf ("Failure to set threads to 1, error %d.\n",status);
	    		exit(0);
	    	}
	    	status = CPXsetintparam (env3, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	    	status = CPXsetintparam (env3, CPX_PARAM_PREIND, CPX_ON);
	    	status = CPXsetincumbentcallbackfunc 	(env3, userincumbent2,  NULL);
	    	status = CPXsetintparam (env3, CPX_PARAM_LANDPCUTS, 	3);
	  	  status = CPXsetintparam (env3, CPX_PARAM_CLIQUES, 	3);
		  status = CPXsetintparam (env3, CPX_PARAM_COVERS, 	3);
		  status = CPXsetintparam (env3, CPX_PARAM_DISJCUTS, 	3);
		  status = CPXsetintparam (env3, CPX_PARAM_FRACCUTS, 	3);
		  status = CPXsetintparam (env3, CPX_PARAM_GUBCOVERS, 	3);
		  status = CPXsetintparam (env3, CPX_PARAM_IMPLBD, 	3);
		  status = CPXsetintparam (env3, CPX_PARAM_MIRCUTS, 	3);
		  status = CPXsetintparam (env3, CPX_PARAM_ZEROHALFCUTS, 3);
		  status = CPXsetintparam (env3, CPX_PARAM_FLOWCOVERS, 	3);
		  status = CPXsetintparam (env3, CPX_PARAM_FLOWPATHS, 	3);
	}
	been_passed = 1;
	status = CPXchgprobtype (env, lp_1, CPXPROB_LP);
/*	printf("status: %d\n",status);*/
	status = CPXchgprobtype (env, lp_2, CPXPROB_LP);
/*	printf("status: %d\n",status);*/
}

void pass_a_mip(CPXLPptr prob1)
{
	int status = 0;
	mip1 = CPXcloneprob (env2, prob1, &status);
}

void free_probs()
{
	CPXfreeprob(env2, &mip1);
	int status = 0;
	status=CPXcloseCPLEX(&env2);
}

void free_coefs()
{
	free(ob_coef1);
	free(ob_coef2);
	free(weighted_coefs);
}

void provide_coefs(double *coef1, double *coef2, int num_cols)
{
	ob_coef1 = (double *) malloc ((num_cols+3)*sizeof(double));
	ob_coef2 = (double *) malloc ((num_cols+3)*sizeof(double));
	weighted_coefs = (double *) malloc ((num_cols+3)*sizeof(double));
	int i;
	for(i=0;i<num_cols+3;i++)
	{
		if(i<num_cols) 
		{
			ob_coef1[i] = coef1[i];
			ob_coef2[i] = coef2[i];
			weighted_coefs[i] = .5*ob_coef1[i]+.5*ob_coef2[i];
		}
		else
		{
			ob_coef1[i] = 0.;
			ob_coef2[i] = 0.;
			weighted_coefs[i] = 0.;
		}
	}
}

void free_xctype()
{
	free(xctype);
	free(integer_indices);
	free(integer_var_scores);
	free(integer_var_last_val);
}

int pure_binary = 0;
double multiplying_factor = 1.;

void provide_xctype(CPXENVptr env, char *types, int num_cols)
{
	xctype = (char *) malloc ((num_cols+3)*sizeof(char));
	integer_indices = (int *) malloc ((num_cols+3)*sizeof(int));
	integer_var_scores = (int *) malloc ((num_cols+3)*sizeof(int));
	integer_var_last_val = (int *) malloc ((num_cols+3)*sizeof(int));
	
	int status = 0;
/*	status = CPXwriteprob (env, lp_1, "myprob1.lp", "LP");*/
/*  	status = CPXwriteprob (env, lp_2, "myprob2.lp", "LP");*/
/*  	exit(0);*/
	
	double *ub = (double *) malloc ((num_cols)*sizeof(double));
	status = CPXgetub (env, lp_1, ub, 0, num_cols-1);
	if ( status ) {
		printf ("(%d) Failed to get ub's. Error code %d\n", __LINE__,status);
		exit(0);
	}
	int num_int_not_binary = 0;
	
	int i;
	int count = 0;
	total_num_integer = 0;
	for(i=0;i<num_cols+3;i++)
	{
		if(i<num_cols) 
		{
			xctype[i] = types[i];
/*			printf("xctype[%d] = %c, ub[%d] = %lf\n",i,xctype[i],i,ub[i]);*/
			if(xctype[i] == 'I' && ub[i] != 1.) num_int_not_binary++;
		}
		else xctype[i] = 'C';
		if( xctype[i] == 'B' || xctype[i] == 'I')
		{
/*			printf("coef: %lf, int of it: %d, difference: %lf\n",multiplying_factor*fabs(ob_coef1[i]), (int) floor(multiplying_factor*fabs(ob_coef1[i])), fabs(multiplying_factor*fabs(ob_coef1[i]) - (int) floor(multiplying_factor*fabs(ob_coef1[i]))));*/
	    		if(there_will_only_be_points && integer_objective == 1 && fabs(ob_coef1[i]) != 0.) 
	    			while( fabs(multiplying_factor*fabs(ob_coef1[i]) - floor(multiplying_factor*fabs(ob_coef1[i]))) > 0.000001 && 
	    				fabs(multiplying_factor*fabs(ob_coef1[i]) - floor(multiplying_factor*fabs(ob_coef1[i]))) < .999999)
		    		{
/*		    			printf("the coef was: %lf\n",ob_coef1[i]);*/
/*		   			printf("scaled: %lf\n",multiplying_factor*fabs(ob_coef1[i]));*/
		    			multiplying_factor = 10.*multiplying_factor;
		    		}
	    		else if(there_will_only_be_points && integer_objective == 2 && fabs(ob_coef2[i]) != 0.) 
	    			while( fabs(multiplying_factor*fabs(ob_coef2[i]) - floor(multiplying_factor*fabs(ob_coef2[i]))) > 0.000001 &&
	    				fabs(multiplying_factor*fabs(ob_coef2[i]) - floor(multiplying_factor*fabs(ob_coef2[i]))) < .999999)
		    		{
/*		    			printf("the coef was: %lf\n",ob_coef1[i]);*/
/*		   			printf("scaled: %lf\n",multiplying_factor*fabs(ob_coef1[i]));*/
		    			multiplying_factor = 10.*multiplying_factor;
		    		}
/*			printf("integer index: %d\n",i);*/
			integer_indices[total_num_integer] = i;
			total_num_integer++;
		}
		integer_var_scores[i] = 0;
		integer_var_last_val[i] = 0;
	}
	
	if(there_will_only_be_points) printf("In order to make all coefficients of integer variables integral, objective %d must be scaled by a factor of %lf.\n",
		integer_objective,multiplying_factor);
	if(num_int_not_binary == 0) pure_binary = 1;
	free(ub);
	printf("number integer vars: %d\n",total_num_integer);
	num_integer = total_num_integer;
}

static void free_and_null (char **ptr)
{
  if ( *ptr != NULL ) {
    free (*ptr);
    *ptr = NULL;
  }
}

/*********************************************************************************************************************** 

	This function is used for changing the weights on the coefficients of the objective functions 
	so that they correspond with a particular slope in the objective space.
	
***********************************************************************************************************************/

void chg_coefs(CPXCENVptr env, CPXLPptr prob, int *indices, double slope)
{
	int i,status;
	if(slope == 0.) 
	{
/*		printf("doing this\n");*/
		for(i=0;i<obj1_index;i++)
		{
			weighted_coefs[i] = ob_coef2[i];
		}
		weighted_coefs[obj1_index] = 0.;
		weighted_coefs[obj2_index] = 0.;
	}
	else if(slope < -10000000.)
	{
/*		printf("doing this2\n");*/
		for(i=0;i<obj1_index;i++)
		{
			weighted_coefs[i] = ob_coef1[i];
		}
		weighted_coefs[obj1_index] = 0.;
		weighted_coefs[obj2_index] = 0.;
	}
	else
	{
		for(i=0;i<obj1_index;i++) weighted_coefs[i] = 0.;
		weighted_coefs[obj1_index] = 1.;
		weighted_coefs[obj2_index] = -1./slope;
	}

/*	if(prob == lp_1 || prob == lp_2) status = CPXchgobj (env, prob, obj1_index, indices, weighted_coefs);*/
/*	else status = CPXchgobj (env, prob, obj1_index+2, indices, weighted_coefs);*/
	status = CPXchgobj (env, prob, obj1_index+2, indices, weighted_coefs);
	if ( status ) {
		printf ("Failed to change obj coef. Error code %d\n", status);
	}
/*	status = CPXwriteprob (env, prob, "myprob2.lp", "LP");*/
}

/*********************************************************************************************************************** 

	These functions are used for qsort.
	
***********************************************************************************************************************/

int comparison_function(const void *a, const void *b)
{
	double a_ = *(double *)a - (int) (*(double *)a);
	double b_ = *(double *)b - (int) (*(double *)b);
/*	printf("a_: %lf, b_:%lf\n",a_,b_);*/
	if(a_ < b_) return 1;
	else return 0;
}

int comparison_function3(const void *a, const void *b)
{
	double a_ = *(double *)a ;
	double b_ = *(double *)b ;
/*	printf("a_: %lf, b_:%lf\n",a_,b_);*/
	if(a_ < b_) return 1;
	else return 0;
}

/*********************************************************************************************************************** 

	This function is used when BB terminates early. It is called at each open node upon termination.
	It generates the dual bound associated with a node and adds it to a data structure for storage.
	
***********************************************************************************************************************/

node *empty_node = NULL;
clock_t build_st, build_fi;

void build_dual_bd_approximately(CPXCENVptr env, CPXLPptr prob1, double slp, int *indices)
{
/*	printf("attempting to build dual bd approximately\n");*/
	build_st = clock();
	double ti = 0.;
	
	double first_obj1_val, first_obj2_val;
/*	int *indices = (int *) malloc (cur_numcols * sizeof (int));*/
	
	/*************** Getting started ***************/
	
	int status, i;
/*	for(i=0;i<cur_numcols;i++) indices[i] = i;*/
	
  	chg_coefs(env,prob1,indices,slp);
	
	/*************** Solve LP ***************/
		
	status = CPXlpopt (env, prob1);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	int lpstat = CPXgetstat (env, prob1);
  	
  	if(lpstat == 12)
  	{
  		double upper_limit,lower_limit,what_it_should_be;
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_UpperObj, &upper_limit);
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_LowerObj, &lower_limit);
  		what_it_should_be = pow(10.,75.);
  		if(upper_limit < what_it_should_be || lower_limit > -what_it_should_be)
  		{
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_UpperObj, what_it_should_be);
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_LowerObj, -what_it_should_be);
  			printf("Warning: Within PSA_full encountered an issue where objective value limits were set to strange values. This was NOT a user error. Resetting these values to default and restarting PSA_full\n");
/*  			PSA_all(env,prob);*/
  			goto TERMINATE;
  		}
  	}
  	
  	double orig_obj_vals[2] = {0.,0.};
  	double prev_obj_vals[2] = {0.,0.};
  	double new_obj_vals[2] = {0.,0.};
  	double coef_lb = 0., coef_ub = 0., coef_lb_orig = 0.;

	status = CPXgetx (env, prob1, prev_obj_vals, obj1_index, obj2_index);
	orig_obj_vals[0] = prev_obj_vals[0];
	orig_obj_vals[1] = prev_obj_vals[1];
	
	int iter_cnt = 0, same_cnt = 0;
/*	while( fabs(prev_obj_vals[0] - obj_vals[0]) > .01 )*/
	for(i=0;i<5;i++)
	{
		iter_cnt++;
		build_fi = clock();
		ti = (double)(build_fi - build_st) / CLOCKS_PER_SEC;
/*		printf("iteration: %d\n",iter_cnt);*/
		if( ti > time_per || (ti > 5. && same_cnt > 25) || (iter_cnt > 100 && same_cnt >= iter_cnt -5))//if(iter_cnt > 1000000)
		{
			if(same_cnt <= 0 && ti > time_per)
			{
				printf("Warning: time to generate dual bound at an individual node has been exceeded. An approximate dual bound is being used for this node!\n");
			}
			printf("same count: %d\n",same_cnt);
/*			printf("%lf, %lf\n",prev_obj_vals[0],obj_vals[0]);*/
			printf("%d\n",__LINE__);
			break;
		}
		if(iter_cnt == 1) status = CPXobjsa(env,prob1,obj2_index,obj2_index,&coef_lb_orig,&coef_ub);
		else status = CPXobjsa(env,prob1,obj2_index,obj2_index,&coef_lb,&coef_ub);
		chg_coefs(env,prob1,indices, -1./(coef_ub + .00001));
		status = CPXlpopt (env, prob1);
		if ( status ) {
			printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
			goto TERMINATE;
		}
		lpstat = CPXgetstat (env, prob1);
		status = CPXgetx (env, prob1, new_obj_vals, obj1_index, obj2_index);
		
		if(new_obj_vals[0] < prev_obj_vals[0])
		{
/*			printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",prev_obj_vals[0],new_obj_vals[0],prev_obj_vals[1],new_obj_vals[1]);*/
			insert_db(2,prev_obj_vals[0],prev_obj_vals[1],new_obj_vals[0],new_obj_vals[1],
				(prev_obj_vals[1]-new_obj_vals[1])/(prev_obj_vals[0]-new_obj_vals[0]),&tree2,&empty_node);
			prev_obj_vals[0] = new_obj_vals[0];
			prev_obj_vals[1] = new_obj_vals[1];
			same_cnt = 0;
		}
		else same_cnt++;
/*		if( fabs(prev_obj_vals[0] - obj_vals[0]) < .1 && same_cnt >=5 ) break;*/
	}
	
	chg_coefs(env,prob1,indices,  -1./(coef_lb - .00001));
	prev_obj_vals[0] = orig_obj_vals[0];
	prev_obj_vals[1] = orig_obj_vals[1];
	
	iter_cnt = 0;
	same_cnt = 0;
	for(i=0;i<5;i++)
	{
		iter_cnt++;
		build_fi = clock();
		ti = (double)(build_fi - build_st) / CLOCKS_PER_SEC;
/*		printf("iteration: %d\n",iter_cnt);*/
		if( ti > time_per || (ti > 5. && same_cnt > 25) || (iter_cnt > 100 && same_cnt >= iter_cnt -5))//if(iter_cnt > 1000000)
		{
			if(same_cnt <= 0 && ti > time_per)
			{
				printf("Warning: time to generate dual bound at an individual node has been exceeded. An approximate dual bound is being used for this node!\n");
			}
			printf("same count: %d\n",same_cnt);
/*			printf("%lf, %lf\n",prev_obj_vals[0],obj_vals[0]);*/
			printf("%d\n",__LINE__);
			break;
		}
		status = CPXobjsa(env,prob1,obj2_index,obj2_index,&coef_lb,&coef_ub);
		chg_coefs(env,prob1,indices, -1./(coef_lb - .00001));
		status = CPXlpopt (env, prob1);
		if ( status ) {
			printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
			goto TERMINATE;
		}
		lpstat = CPXgetstat (env, prob1);
		status = CPXgetx (env, prob1, new_obj_vals, obj1_index, obj2_index);
		
		if(new_obj_vals[0] > prev_obj_vals[0])
		{
/*			printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",prev_obj_vals[0],new_obj_vals[0],prev_obj_vals[1],new_obj_vals[1]);*/
			insert_db(2,new_obj_vals[0],new_obj_vals[1],prev_obj_vals[0],prev_obj_vals[1],
				(prev_obj_vals[1]-new_obj_vals[1])/(prev_obj_vals[0]-new_obj_vals[0]),&tree2,&empty_node);
			prev_obj_vals[0] = new_obj_vals[0];
			prev_obj_vals[1] = new_obj_vals[1];
			same_cnt = 0;
		}
		else same_cnt++;
/*		if( fabs(prev_obj_vals[0] - obj_vals[0]) < .1 && same_cnt >=5 ) break;*/
	}
 	
 	TERMINATE:
 	;
}

void build_dual_bd(CPXCENVptr env, CPXLPptr prob)
{
/*	printf("attempting to build dual bd\n");*/
	build_st = clock();
	double ti = 0.;
	if( !prob ) prob = lp_1;
	
	double first_obj1_val, first_obj2_val;
	int *indices = (int *) malloc (cur_numcols * sizeof (int));
	
	/*************** Getting started ***************/
	
	int status, i;
	CPXLPptr prob1 = CPXcloneprob (env, prob, &status);
	for(i=0;i<cur_numcols;i++) indices[i] = i;
	
	double obj1_ub = 0., obj1_lb = 0., obj2_ub = 0., obj2_lb = 0.;
	
	status = CPXgetlb (env, prob1, &obj1_lb, obj1_index, obj1_index);
	status = CPXgetub (env, prob1, &obj1_ub, obj1_index, obj1_index);
	
	status = CPXgetlb (env, prob1, &obj2_lb, obj2_index, obj2_index);
	status = CPXgetub (env, prob1, &obj2_ub, obj2_index, obj2_index);
	
	if(max_time < .1 || max_nodes <= 1) //if(obj1_lb <= NW_extreme_x || (exploit_objective_gaps == 0 && obj2_ub  >= NW_extreme_y))
	{
		int chg_ind[2] = {obj1_index,obj2_index};
		char chg_lu[2] = {'L','U'};
		double chg_bd[2] = {-CPX_INFBOUND,CPX_INFBOUND};
		status = CPXchgbds (env, prob1, 2, chg_ind, chg_lu, chg_bd);
	}
/*	if(obj1_ub  >= SE_extreme_x || (exploit_objective_gaps == 0 && obj2_lb <= SE_extreme_y))*/
/*	{*/
/*		int chg_ind[2] = {obj1_index,obj2_index};*/
/*		char chg_lu[2] = {'U','L'};*/
/*		double chg_bd[2] = {CPX_INFBOUND,-CPX_INFBOUND};*/
/*		status = CPXchgbds (env, prob1, 2, chg_ind, chg_lu, chg_bd);*/
/*	}*/
	
	if(approximate_dual_bd) 
	{
		build_dual_bd_approximately(env, prob1, (obj2_ub-obj2_lb)/(obj1_lb-obj1_ub), indices);
		goto TERMINATE;
	}
	/*************** Build prob 2 ***************/
	
/*	START_DB_OVER:*/
/*	;*/
	
	CPXLPptr prob2 = CPXcloneprob (env, prob1, &status);
  	
  	chg_coefs(env,prob1,indices,-10000000.);
  	chg_coefs(env,prob2,indices,-.000000001);
	
	/*************** Solve LP1 ***************/
		
	status = CPXlpopt (env, prob1);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	int lpstat = CPXgetstat (env, prob1);
  	
  	if(lpstat == 12)
  	{
  		double upper_limit,lower_limit,what_it_should_be;
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_UpperObj, &upper_limit);
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_LowerObj, &lower_limit);
  		what_it_should_be = pow(10.,75.);
  		if(upper_limit < what_it_should_be || lower_limit > -what_it_should_be)
  		{
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_UpperObj, what_it_should_be);
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_LowerObj, -what_it_should_be);
  			printf("Warning: Within PSA_full encountered an issue where objective value limits were set to strange values. This was NOT a user error. Resetting these values to default and restarting PSA_full\n");
  			PSA_all(env,prob);
  			goto TERMINATE;
  		}
  	}
  	
  	double prev_obj_vals[2] = {0.,0.};
  	double new_obj_vals[2] = {0.,0.};
  	double obj_vals[2] = {0.,0.};
  	double coef_lb = 0., coef_ub = 0.;

	status = CPXgetx (env, prob1, prev_obj_vals, obj1_index, obj2_index);
	
	status = CPXlpopt (env, prob2);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	status = CPXgetx (env, prob2, obj_vals, obj1_index, obj2_index);
  	
  	if( fabs(obj_vals[0] - prev_obj_vals[0]) < .0001 || fabs(obj_vals[1] - prev_obj_vals[1]) < .0001 ) 
  	{
/*  		printf("%lf,%lf,%lf,%lf ... %lf,%lf\n",SE_extreme_x,SE_extreme_y,NW_extreme_x,NW_extreme_y,obj_vals[0],obj_vals[1]);*/
/*  		if(first_time && (obj_vals[0] >= SE_extreme_x || obj_vals[1] >= NW_extreme_y))*/
/*  		{*/
/*  			printf("changing and starting over\n");*/
/*  			int chg_ind[4] = {obj1_index,obj1_index,obj2_index,obj2_index};*/
/*			char chg_lu[4] = {'L','U','L','U'};*/
/*			double chg_bd[4] = {-CPX_INFBOUND,CPX_INFBOUND,-CPX_INFBOUND,CPX_INFBOUND};*/
/*			status = CPXchgbds (env, prob1, 4, chg_ind, chg_lu, chg_bd);*/
/*			first_time = 0;*/
/*			goto START_DB_OVER;*/
/*  		}*/
/*  		printf("we didn't start over, or its the 2nd time. Check: %d\n",1-first_time);*/
/*		printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",obj_vals[0],obj_vals[0],obj_vals[1],obj_vals[1]);*/
  		insert_db(2,obj_vals[0],obj_vals[1],obj_vals[0],obj_vals[1],0.,&tree2,&empty_node);
  		goto TERMINATE;
  	}
	
	int iter_cnt = 0, same_cnt = 0;
	while( fabs(prev_obj_vals[0] - obj_vals[0]) > .01 )
	{
		iter_cnt++;
		build_fi = clock();
		ti = (double)(build_fi - build_st) / CLOCKS_PER_SEC;
/*		printf("iteration: %d\n",iter_cnt);*/
/*		if( ti > time_per || (ti > 5. && same_cnt > 25) || (iter_cnt > 100 && same_cnt >= iter_cnt -5))//if(iter_cnt > 1000000)*/
/*		{*/
/*			if(same_cnt <= 0 && ti > time_per)*/
/*			{*/
/*				printf("Warning: time to generate dual bound at an individual node has been exceeded. An approximate dual bound is being used for this node!\n");*/
/*			}*/
/*			printf("same count: %d\n",same_cnt);*/
/*			printf("%lf, %lf\n",prev_obj_vals[0],obj_vals[0]);*/
/*			printf("%d\n",__LINE__);*/
/*			break;*/
/*		}*/
		status = CPXobjsa(env,prob1,obj2_index,obj2_index,&coef_lb,&coef_ub);
		chg_coefs(env,prob1,indices, -1./(coef_ub + .00001));
		status = CPXlpopt (env, prob1);
		if ( status ) {
			printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
			goto TERMINATE;
		}
		lpstat = CPXgetstat (env, prob1);
		status = CPXgetx (env, prob1, new_obj_vals, obj1_index, obj2_index);
		
		if(new_obj_vals[0] < prev_obj_vals[0])
		{
/*			printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",prev_obj_vals[0],new_obj_vals[0],prev_obj_vals[1],new_obj_vals[1]);*/
			insert_db(2,prev_obj_vals[0],prev_obj_vals[1],new_obj_vals[0],new_obj_vals[1],
				(prev_obj_vals[1]-new_obj_vals[1])/(prev_obj_vals[0]-new_obj_vals[0]),&tree2,&empty_node);
			prev_obj_vals[0] = new_obj_vals[0];
			prev_obj_vals[1] = new_obj_vals[1];
			same_cnt = 0;
		}
		else same_cnt++;
		if( fabs(prev_obj_vals[0] - obj_vals[0]) < .1 && same_cnt >=5 ) break;
	}
 	
 	TERMINATE:
 	
	free_and_null ((char **) &indices);
	CPXfreeprob(env,&prob1);
	CPXfreeprob(env,&prob2);
}

/*********************************************************************************************************************** 

	This function can be used to generate the dual bound associated with a node. It is really only used
	for debugging.
	
***********************************************************************************************************************/

void PSA_all(CPXCENVptr env, CPXLPptr prob)
{
	if( !prob ) prob = lp_1;
	
	double first_obj1_val, first_obj2_val;
	int *indices = (int *) malloc (cur_numcols * sizeof (int));
	
	/*************** Getting started ***************/
	
	int status, i;
	CPXLPptr prob1 = CPXcloneprob (env, prob, &status);
  	for(i=0;i<cur_numcols;i++) indices[i] = i;
	
	double obj1_ub = 0., obj1_lb = 0.;
	
	status = CPXgetlb (env, prob1, &obj1_lb, obj1_index, obj1_index);
	status = CPXgetub (env, prob1, &obj1_ub, obj1_index, obj1_index);
	
	if(obj1_lb <= NW_extreme_x || obj1_ub  >= SE_extreme_x)
	{
		int chg_ind[4] = {obj1_index,obj1_index,obj2_index,obj2_index};
		char chg_lu[4] = {'L','U','L','U'};
		double chg_bd[4] = {-CPX_INFBOUND,CPX_INFBOUND,-CPX_INFBOUND,CPX_INFBOUND};
		status = CPXchgbds (env, prob1, 4, chg_ind, chg_lu, chg_bd);
	}
	
	/*************** Build prob 2 ***************/
	
	CPXLPptr prob2 = CPXcloneprob (env, prob1, &status);
  	
  	chg_coefs(env,prob1,indices,-10000000.);
  	chg_coefs(env,prob2,indices,-.000000001);
/*  	status = CPXchgobj (env, prob2, numcols, indices, ob_coef2);*/
/*	if ( status ) {*/
/*		printf ("Failed to get change obj coef. Error code %d\n", status);*/
/*		goto TERMINATE;*/
/*	}*/
	
/*	status = CPXwriteprob (env, prob1, "myprob1.lp", "LP");*/
/*	status = CPXwriteprob (env, prob2, "myprob2.lp", "LP");*/
/*	exit(0);*/

/*	printf("problem 1 obj2_index coef: %lf\n",1./1000.);*/
	
	/*************** Solve LP1 ***************/
		
	status = CPXlpopt (env, prob1);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	int lpstat = CPXgetstat (env, prob1);
  	
  	if(lpstat == 12)
  	{
  		double upper_limit,lower_limit,what_it_should_be;
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_UpperObj, &upper_limit);
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_LowerObj, &lower_limit);
  		what_it_should_be = pow(10.,75.);
  		if(upper_limit < what_it_should_be || lower_limit > -what_it_should_be)
  		{
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_UpperObj, what_it_should_be);
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_LowerObj, -what_it_should_be);
  			printf("Warning: Within PSA_full encountered an issue where objective value limits were set to strange values. This was NOT a user error. Resetting these values to default and restarting PSA_full\n");
  			PSA_all(env,prob);
  			goto TERMINATE;
  		}
  	}
  	
  	double prev_obj_vals[2] = {0.,0.};
  	double new_obj_vals[2] = {0.,0.};
  	double obj_vals[2] = {0.,0.};
  	double coef_lb = 0., coef_ub = 0.;
/*  	printf("the status of the solve: %d\n",lpstat);*/

	status = CPXgetx (env, prob1, prev_obj_vals, obj1_index, obj2_index);
	
/*	printf("plot([%lf],[%lf],'go');\n",prev_obj_vals[0],prev_obj_vals[1]);*/
	
	status = CPXlpopt (env, prob2);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	status = CPXgetx (env, prob2, obj_vals, obj1_index, obj2_index);
  	
  	if( fabs(obj_vals[0] - prev_obj_vals[0]) > .0001 && fabs(obj_vals[1] - prev_obj_vals[1]) > .0001 )
  	{
/*  		printf("obj1 and obj2 solution are different. Keep going\n");*/
/*  		printf("plot([%lf],[%lf],'go');\n",obj_vals[0],obj_vals[1]);*/
  	}
  	else
  	{
/*  		printf("obj1 and obj2 solutions are the same. Stop\n");*/
  		goto TERMINATE;
  	}
	
	int iter_cnt = 0, same_cnt = 0;
	while( fabs(prev_obj_vals[0] - obj_vals[0]) > .01 )
	{
		iter_cnt++;
		if(iter_cnt > 100000)
		{
			printf("%lf, %lf\n",prev_obj_vals[0],obj_vals[0]);
			printf("%d\n",__LINE__);
			break;
			exit(0);
		}
		status = CPXobjsa(env,prob1,obj2_index,obj2_index,&coef_lb,&coef_ub);
/*		printf("bounds on coefficient: %lf to %lf\n",coef_lb,coef_ub);*/
		chg_coefs(env,prob1,indices, -1./(coef_ub + .00001));
		status = CPXlpopt (env, prob1);
		if ( status ) {
			printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
			goto TERMINATE;
		}
		lpstat = CPXgetstat (env, prob1);
/*	  	printf("the status of the solve: %d\n",lpstat);*/
		status = CPXgetx (env, prob1, new_obj_vals, obj1_index, obj2_index);
/*		printf("the next point might be: (%lf,%lf)\n",new_obj_vals[0],new_obj_vals[1]);*/
		
		if(new_obj_vals[0] < prev_obj_vals[0])
		{
			printf("plot([%lf,%lf],[%lf,%lf],'-o');\n",prev_obj_vals[0],new_obj_vals[0],prev_obj_vals[1],new_obj_vals[1]);
			prev_obj_vals[0] = new_obj_vals[0];
			prev_obj_vals[1] = new_obj_vals[1];
			same_cnt = 0;
		}
		else same_cnt++;
		if( fabs(prev_obj_vals[0] - obj_vals[0]) < .1 && same_cnt >=5 ) break;
	}
 	
 	TERMINATE:
 	
 	status = CPXsetintparam ( (CPXENVptr) env, CPX_PARAM_ITLIM, 2100000000);
  	if ( status ) {
    		printf ("%s(%d): CPXsetintparam, Failed to set iterations limit of SIMPLEX, error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
 	
	free_and_null ((char **) &indices);
	CPXfreeprob(env,&prob1);
	CPXfreeprob(env,&prob2);
}

/*********************************************************************************************************************** 

	This function is used for generating integer feasible line segments that are part of a node's dual 
	bound. The function exits when an integer variable changes value. (Works from right to left)
	
***********************************************************************************************************************/
clock_t PSA_st, PSA_fi;
double coef_add = .00001;

int PSA_iter_cnt = 0, PSA_left_iter_cnt = 0;

int PSA(CPXCENVptr env, CPXLPptr prob, double *x_orig, int *basis_col_info_orig, int *basis_row_info_orig, int snum, int same_seq,  CPXLPptr original_prob)
{
	int first_it = 1;
	int retval = 0;
	PSA_iter_cnt = 0;
	
	if(generate_disjunctive_cuts_from_obj_space_disjunctions && !temp_x_r) temp_x_r = (double *) malloc (cur_numcols*sizeof(double));
	
/*	printf("calling PSA from seqnum: %d\n",snum);*/
	/*************** Getting started ***************/
	
	int status;
	CPXLPptr prob1 = CPXcloneprob (env, prob, &status);
	CPXchgprobtype(env, prob1, CPXPROB_LP);
	
	/*************** Fix the integer variables ***************/
	
	char up[1] = {'U'};
	
	status = CPXchgbds (env, prob1, 1, &obj1_index, up, &x_orig[obj1_index]);
	if ( status ){
		 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
		 goto TERMINATE;}
	
	int i;
	int *indices;
	indices = (int *) malloc (cur_numcols * sizeof (int));
	double *x = (double *) malloc (cur_numcols * sizeof (double));
	double *x_new = (double *) malloc (cur_numcols * sizeof (double));
	int *basis_col_info = (int *) malloc (cur_numcols * sizeof (int));
	int *basis_row_info = (int *) malloc (cur_numrows * sizeof (int));
	double coef_lb = 0., coef_ub = 0.;
	
	for(i=0;i<cur_numcols;i++)
	{
		x[i] = x_orig[i];
		basis_col_info[i] = basis_col_info_orig[i];
	}
	for(i=0;i<cur_numrows;i++) basis_row_info[i] = basis_row_info_orig[i];
	
	
	/*************** Build prob 2 ***************/
	
	CPXLPptr prob2 = CPXcloneprob (env, prob1, &status);
  	for(i=0;i<cur_numcols;i++) indices[i] = i;
  	
  	chg_coefs(env,prob1,indices,-10000000.);
  	chg_coefs(env,prob2,indices,-.000000001);
	
	/*************** Solve LP1 ***************/
	
	status = CPXcopybase (env, prob1, basis_col_info, basis_row_info);
  	if ( status ) {
    		printf ("%s(%d): CPXcopybase, Failed to copy basis, error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	status = CPXcopybase (env, prob2, basis_col_info, basis_row_info);
  	if ( status ) {
    		printf ("%s(%d): CPXcopybase, Failed to copy basis, error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	status = CPXlpopt (env, prob1);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	int lpstat = CPXgetstat (env, prob1);
  	
  	if(lpstat == 3) 
  	{
  		retval = 2;	
  		goto TERMINATE;
  	}
  	
  	status = CPXgetx (env, prob1, x, 0, cur_numcols-1);
  	if ( status ) {
   		printf ("%s(%d): CPXgetx, Failed to get x-values, error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
	
/*	printf("plot([%lf],[%lf],'go');\n",x[obj1_index],x[obj2_index]);*/
	
	status = CPXlpopt (env, prob2);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	double prob2_obvals[2] = {0.,0.};
  	
  	status = CPXgetx (env, prob2, prob2_obvals, obj1_index, obj2_index);
  	
  	if( fabs(x[obj1_index] - prob2_obvals[0]) > .0001 && fabs(x[obj2_index] - prob2_obvals[1]) > .0001 )
  	{
/*  		printf("obj1 and obj2 solution are different. Keep going\n");*/
/*  		printf("plot([%lf],[%lf],'go');\n",prob2_obvals[0],prob2_obvals[1]);*/
  	}
  	else
  	{
/*  		printf("obj1 and obj2 solutions are the same. Stop\n");*/
/*		printf("plot([%lf],[%lf],'go');\n",x[obj1_index],x[obj2_index]);*/
  		insert(1,x[obj1_index],x[obj2_index],0.,0.,0.,&tree,&empty_node);
  		retval = 2;
  		goto TERMINATE;
  	}
  	
  	for(i=0;i<cur_numcols;i++) x_new[i] = x[i];
  	
  	double overall_coef_lb = 0., prev_coef_lb = -100000000000000000., add_val = .001;
	
	while( fabs(prob2_obvals[0] - x_new[obj1_index]) > .01 )
	{
		PSA_iter_cnt++;
		if(PSA_iter_cnt > 5000)
		{
			printf("%lf, %lf\n",prob2_obvals[0],x_new[obj1_index]);
/*			printf("%d\n",__LINE__);*/
			goto TERMINATE;
		}
		status = CPXobjsa(env,prob1,obj2_index,obj2_index,&coef_lb,&coef_ub);
/*		printf("bounds on coefficient: %lf to %lf\n",coef_lb,coef_ub);*/
		if(fabs(prev_coef_lb - coef_lb) < .00000001)
		{
			coef_ub += add_val;
			add_val = add_val*2;
/*			retval = 0;*/
/*			goto TERMINATE;*/
		}
		prev_coef_lb = coef_lb;
		if(coef_ub > 10000000000.)
		{
			retval = 2;
			goto TERMINATE;
		}
		chg_coefs(env,prob1,indices, -1./(coef_ub + .00001));
		status = CPXlpopt (env, prob1);
		if ( status ) {
			printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
			goto TERMINATE;
		}
		lpstat = CPXgetstat (env, prob1);
/*	  	printf("the status of the solve: %d\n",lpstat);*/
		status = CPXgetx (env, prob1, x_new, 0, cur_numcols-1);
/*		printf("the next point might be: (%lf,%lf)\n",new_obj_vals[0],new_obj_vals[1]);*/
		
		if(x_new[obj1_index] - x[obj1_index] < -.000001 )
		{
			overall_coef_lb = coef_lb;
			int num_diff = 0;
			for(i=0;i<total_num_integer;i++)
			{
/*				printf(" x%d: %lf\n",integer_indices[i],x[integer_indices[i]]);*/
/*				printf(" x%d: %lf\n",integer_indices[i],x_new[integer_indices[i]]);*/
				double diff = x_new[integer_indices[i]]-floor(x_new[integer_indices[i]]);
				if(diff >= .00001 && diff <= .99999) 
				{
/*					printf("from within PSA(2): since it became fractional, setting branch var to %d for seqnum %d\n",best_index,snum);*/
					branch_seqnum = snum;
/*					printf("changing frac index to %d (%d)\n",integer_indices[i],__LINE__);*/
/*					printf("plot(%lf,%lf,'co');\n",x_new[obj1_index],x_new[obj2_index]);*/
					frac_index = integer_indices[i];
					frac_val = x_new[integer_indices[i]];
					if(frac_scores[i] > 0.0001) 
					{
/*						printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
						frac_scores[i] += within_PSA_score;
/*						printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
					}
					else
					{
						frac_scores[i] = within_PSA_score + multiplier*integer_indices[i];
/*						printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
						num_frac++;
					}
					frac_values[integer_indices[i]] = frac_val;
/*					printf("num_frac: %d\n",num_frac);*/
					num_diff++;
/*					printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
				}
				else if( fabs(x[integer_indices[i]]-x_new[integer_indices[i]]) >= .00001)
				{
					branch_seqnum = snum;
/*					printf("changing frac index to %d (%d)\n",integer_indices[i],__LINE__);*/
/*					printf("plot(%lf,%lf,'mo');\n",x[obj1_index],x[obj2_index]);*/
/*					printf("plot(%lf,%lf,'mo');\n",x_new[obj1_index],x_new[obj2_index]);*/
					frac_index = integer_indices[i];
					frac_val = fmax(x[integer_indices[i]],x_new[integer_indices[i]])-.5;
					if(frac_scores[i] > 0.0001) 
					{
						frac_scores[i] += within_PSA_score;
/*						printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
					}
					else
					{
						frac_scores[i] = within_PSA_score + multiplier*integer_indices[i];
/*						printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
						num_frac++;
					}
					frac_values[integer_indices[i]] = frac_val;
/*					printf("num_frac: %d\n",num_frac);*/
					num_diff++;
/*					printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
				}
			}
			if(slope_scores_in_psa && num_diff > 1) 
			{
/*				printf("number of integer variables that changed: %d\n",num_diff);*/
				double *new_up = (double *) malloc (total_num_integer* sizeof (double));
				double *new_low = (double *) malloc (total_num_integer * sizeof (double));
				char *new_lu = (char *) malloc (total_num_integer * sizeof (char));
				double *slope_score = (double *) malloc (num_diff* sizeof (double));
				int counter = 0, changed_up_or_down = 0;
				
				for(i=0;i<total_num_integer;i++)
				{
					new_up[i] = x[integer_indices[i]];
					new_low[i] = x[integer_indices[i]];
				}
				
				CPXLPptr temp_lp = CPXcloneprob (env3, prob1, &status);
				
				for(i=0;i<total_num_integer;i++) new_lu[i] = 'L';
				status = CPXchgbds (env3, temp_lp, total_num_integer, integer_indices, new_lu, new_low);
		    		if ( status ){
		     			 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
		     			 goto TERMINATE;}
		     			 
		     		for(i=0;i<total_num_integer;i++) new_lu[i] = 'U';
		     		status = CPXchgbds (env3, temp_lp, total_num_integer, integer_indices, new_lu, new_up);
		    		if ( status ){
		     			 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
		     			 goto TERMINATE;}
				
				int j = 0;
				int change_indices[2] = {0,0};
				char change_lu[2] = {'L','U'};
				double change_bds[2] = {0.,0.};
				while(counter < num_diff)
				{
					for(i=j;i<total_num_integer;i++)
					{
						if( x[integer_indices[i]] < x_new[integer_indices[i]] )
						{
							changed_up_or_down = 0;
							change_indices[0] = integer_indices[i];
							change_indices[1] = integer_indices[i];
							change_bds[0] = x[integer_indices[i]] + 1.;
							change_bds[1] = x[integer_indices[i]] + 1.;
							j = i+1;
							break;
						}
						else if( x[integer_indices[i]] > x_new[integer_indices[i]] )
						{
							changed_up_or_down = 1;
							change_indices[0] = integer_indices[i];
							change_indices[1] = integer_indices[i];
							change_bds[0] = x[integer_indices[i]] - 1.;
							change_bds[1] = x[integer_indices[i]] - 1.;
							j = i+1;
							break;
						}
					}
					
					status = CPXchgbds (env3, temp_lp, 2, change_indices, change_lu, change_bds);
			    		if ( status ){
			     			 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
			     			 goto TERMINATE;}
			     		
			     		status = CPXlpopt(env3, temp_lp);
				  	if ( status ) {
				    		printf ("Failed to optimize MIP, error code %d\n", status);
				    		goto TERMINATE;
				  	}
				  	int lpstat = CPXgetstat (env3, temp_lp);
				  	
				  	double objvals[2] = {0.,0.};
				  	double slope = 0.;
				  	if(lpstat != 3)
				  	{
				  		status = CPXgetx (env3, temp_lp, objvals, obj1_index, obj2_index);
/*				  		printf("plot([%lf],[%lf],'o');\n",objvals[0],objvals[1]);*/
				  		slope = (objvals[1]-x[obj2_index])/(objvals[0]-x[obj1_index]);
				  	}
				  	if(slope < -99999) slope = -99999.;
				  	if(slope >= 0.) slope_score[counter] = j-1;
				  	else slope_score[counter] = j-1 - (-1. - (-100000. - slope)/(100000.));
				  	
				  	if(changed_up_or_down)
				  	{
				  		change_bds[0] += 1.;
				  		change_bds[1] += 1.;
				  	}
				  	else 
				  	{
				  		change_bds[1] -= 1.;
				  		change_bds[0] -= 1.;
				  	}
				  	status = CPXchgbds (env3, temp_lp, 2, change_indices, change_lu, change_bds);
			    		if ( status ){
			     			 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
			     			 goto TERMINATE;}
				  	
				  	counter++;
				}
				
/*				for(i=0;i<num_diff;i++) printf("slope_score%d: %lf\n",i,slope_score[i]);*/
		     		qsort(slope_score, num_diff, sizeof(double), comparison_function);
/*		     		for(i=0;i<num_diff;i++) printf("slope_score%d: %lf\n",i,slope_score[i]);*/
		     		
		     		int best_index_based_on_slope = (int) slope_score[0];	
/*		     		printf("best index: %d\tfrac score: %lf\n",best_index_based_on_slope,frac_scores[best_index_based_on_slope]);*/
		     		frac_scores[best_index_based_on_slope] += 2;
/*		     		printf("best index: %d\tfrac score: %lf\n",best_index_based_on_slope,frac_scores[best_index_based_on_slope]);*/
		     		
		     		if(slope_score) free_and_null ((char **) &slope_score);
		     		if(new_lu) free_and_null ((char **) &new_lu);
		     		if(new_up) free_and_null ((char **) &new_up);
		     		if(new_low) free_and_null ((char **) &new_low);
		     		
				goto TERMINATE;
			}
			if(num_diff) goto TERMINATE;
			
/*			printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x[obj1_index],x_new[obj1_index],x[obj2_index],x_new[obj2_index]);*/
			insert(2,x[obj1_index],x[obj2_index],x_new[obj1_index],x_new[obj2_index],
			      (x[obj1_index]-x_new[obj1_index])/(x[obj2_index]-x_new[obj2_index]),&tree,&empty_node);
			      
			for(i=0;i<obj2_index+1;i++) 
			{
				x[i] = x_new[i];
				if(generate_disjunctive_cuts_from_obj_space_disjunctions && PSA_iter_cnt == 1) temp_x_r[i] = x[i];
			}
			
			first_it = 0;
			sub_pr1_x_ub = x[obj1_index];
			sub_pr1_y_lb = x[obj2_index];
		}
	}
	
	if(snum != 0 && same_seq != 1) retval = 2;
 	
 	TERMINATE:
 	
 	if(first_it != 1) 
 	{
		if(retval != 2)
		{
			int add_check = mock_insert(1,x[obj1_index],x[obj2_index],0,0,0.,&tree);
			if(!add_check)
			{
/*				printf("from within PSA, last found value is dominated\n");*/
/*				printf("plot([%lf],[%lf],'ko');\n",x[obj1_index],x[obj2_index]);*/
				if(x[obj1_index] != x_new[obj1_index]) 
				{
/*					printf("checking other soln\n");*/
					add_check = mock_insert(2,x[obj1_index],x[obj2_index],x_new[obj1_index],x_new[obj2_index],
						(x[obj1_index]-x_new[obj1_index])/(x[obj2_index]-x_new[obj2_index]),&tree);
				}
				if(!add_check)
				{
/*					printf("running PSA reduce right to see if the remainder of the bound is dominated\n");*/
/*					printf("plot([%lf],[%lf],'ko');\n",x_new[obj1_index],x_new[obj2_index]);*/

					status = CPXgetbase (env, prob1, basis_col_info, basis_row_info);
				      	if ( status ) {
						printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
				      	}

					retval = PSA_reduce_right(env, prob1, x_new, basis_col_info, basis_row_info, indices, snum);
/*					printf("value of the return: %d\n",retval);*/
				}
			}
		}

  		if(retval != 2) retval = 1;
	}
 	
	if(x) free_and_null ((char **) &x);
	if(indices) free_and_null ((char **) &indices);
	if(x_new) free_and_null ((char **) &x_new);
	if(basis_col_info) free_and_null ((char **) &basis_col_info);
	if(basis_row_info) free_and_null ((char **) &basis_row_info);
	CPXfreeprob(env, &prob1);
	CPXfreeprob(env, &prob2);

	return(retval);
}

/*********************************************************************************************************************** 

	This function is used for generating integer feasible line segments that are part of a node's dual 
	bound. The function exits when an integer variable changes value. (Works from left to right)
	
***********************************************************************************************************************/

int PSA_left(CPXCENVptr env, CPXLPptr prob, double *x_orig, int *basis_col_info_orig, int *basis_row_info_orig, int snum, int same_seq,  CPXLPptr original_prob)
{
	int first_it = 1;
	int retval = 0;
/*	printf("calling PSA left from seqnum: %d\n",snum);*/
	int numcols = CPXgetnumcols (env, prob);
	int numrows = CPXgetnumrows (env, prob);
	
	if(generate_disjunctive_cuts_from_obj_space_disjunctions && !temp_x_l) temp_x_l = (double *) malloc (cur_numcols*sizeof(double));
	
	/*************** Getting started ***************/
	
	int status;
	PSA_left_iter_cnt = 0;
	CPXLPptr prob1 = CPXcloneprob (env, prob, &status);
	CPXchgprobtype(env, prob1, CPXPROB_LP);
	
	char up[1] = {'U'};
	
	status = CPXchgbds (env, prob1, 1, &obj2_index, up, &x_orig[obj2_index]);
	if ( status ){
		 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
		 goto TERMINATE;}
	
	/*************** Fix the integer variables ***************/
	
	int i;
	int *indices;
	indices = (int *) malloc (cur_numcols * sizeof (int));
	double *x = (double *) malloc (cur_numcols * sizeof (double));
	double *x_new = (double *) malloc (cur_numcols * sizeof (double));
	int *basis_col_info = (int *) malloc (cur_numcols * sizeof (int));
	int *basis_row_info = (int *) malloc (cur_numrows * sizeof (int));
	double coef_lb = 0., coef_ub = 0.;
	
	for(i=0;i<cur_numcols;i++)
	{
		x[i] = x_orig[i];
		basis_col_info[i] = basis_col_info_orig[i];
	}
	for(i=0;i<cur_numrows;i++) basis_row_info[i] = basis_row_info_orig[i];
	
	/*************** Build prob 2 ***************/
	
	CPXLPptr prob2 = CPXcloneprob (env, prob1, &status);
  	for(i=0;i<cur_numcols;i++) indices[i] = i;
  	
  	chg_coefs(env,prob2,indices,-10000000.);
  	chg_coefs(env,prob1,indices,-.000000001);
	
/*	if(snum == 359){*/
/*	status = CPXwriteprob (env, prob1, "myprob1.lp", "LP");*/
/*	status = CPXwriteprob (env, prob2, "myprob2.lp", "LP");*/
/*	}*/
	
	/*************** Solve LP1 ***************/
	
	status = CPXcopybase (env, prob1, basis_col_info, basis_row_info);
  	if ( status ) {
    		printf ("%s(%d): CPXcopybase, Failed to copy basis, error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	status = CPXlpopt (env, prob1);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	int lpstat = CPXgetstat (env, prob1);
	
	if(lpstat == 3) 
	{
		retval = 2;
		goto TERMINATE;
	}
  	
  	status = CPXgetx (env, prob1, x, 0, cur_numcols-1);
  	if ( status ) {
   		printf ("%s(%d): CPXgetx, Failed to get x-values, error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
/*  	printf("plot([%lf],[%lf],'go');\n",x[obj1_index],x[obj2_index]);*/
	
	/*************** Get info from LP2 ***************/
	
	status = CPXcopybase (env, prob2, basis_col_info, basis_row_info);
  	if ( status ) {
    		printf ("%s(%d): CPXcopybase, Failed to copy basis, error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}

  	status = CPXlpopt (env, prob2);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
 
 	double prob2_obvals[2] = {0.,0.};
  	
  	status = CPXgetx (env, prob2, prob2_obvals, obj1_index, obj2_index);
  	
/*  	printf("plot([%lf],[%lf],'go');\n",prob2_obvals[0],prob2_obvals[1]);*/
  	
  	if( fabs(x[obj1_index] - prob2_obvals[0]) > .0001 && fabs(x[obj2_index] - prob2_obvals[1]) > .0001 )
  	{
/*  		printf("obj1 and obj2 solution are different. Keep going\n");*/
/*  		printf("plot([%lf],[%lf],'go');\n",prob2_obvals[0],prob2_obvals[1]);*/
  	}
  	else
  	{
/*  		printf("obj1 and obj2 solutions are the same. Stop\n");*/
/*		printf("plot([%lf],[%lf],'go');\n",x[obj1_index],x[obj2_index]);*/
  		insert(1,x[obj1_index],x[obj2_index],0.,0.,0.,&tree,&empty_node);
  		retval = 2;
  		goto TERMINATE;
  	}
  	
  	for(i=0;i<cur_numcols;i++) x_new[i] = x[i];
  	
  	lpstat = 0;
  	double overall_coef_ub = 0., prev_coef_ub = 100000000000000000., add_val = .001;
	
	while( fabs(prob2_obvals[0] - x_new[obj1_index]) > .01 )
	{
		PSA_left_iter_cnt++;
		if(PSA_left_iter_cnt > 5000)
		{
			printf("%lf, %lf\n",prob2_obvals[0],x_new[obj1_index]);
/*			printf("%d\n",__LINE__);*/
			goto TERMINATE;
		}
		status = CPXobjsa(env,prob1,obj2_index,obj2_index,&coef_lb,&coef_ub);
/*		printf("bounds on coefficient: %lf to %lf\n",coef_lb,coef_ub);*/
		if(fabs(prev_coef_ub - coef_ub) < .00000001)
		{
			coef_lb -= add_val;
			add_val = add_val*2;
/*			retval = 0;*/
/*			goto TERMINATE;*/
		}
		prev_coef_ub = coef_ub;
		if(coef_lb < -10000000000.)
		{
			retval = 2;
			goto TERMINATE;
		}
		
		chg_coefs(env,prob1,indices, -1./(coef_lb - .00001));
		status = CPXlpopt (env, prob1);
		if ( status ) {
			printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
			goto TERMINATE;
		}
		lpstat = CPXgetstat (env, prob1);
/*	  	printf("the status of the solve: %d\n",lpstat);*/
		status = CPXgetx (env, prob1, x_new, 0, cur_numcols-1);
/*		printf("the next point might be: (%lf,%lf)\n",new_obj_vals[0],new_obj_vals[1]);*/
		
		if(x[obj1_index] - x_new[obj1_index] < -.000001 )
		{
			overall_coef_ub = coef_ub;
			int num_diff = 0;
			for(i=0;i<total_num_integer;i++)
			{
/*				printf(" x%d: %lf\n",integer_indices[i],x_new[integer_indices[i]]);*/
/*				printf(" x%d: %lf\n",integer_indices[i],x[integer_indices[i]]);*/
				double diff = x_new[integer_indices[i]]-floor(x_new[integer_indices[i]]);
				if(diff >= .00001 && diff <= .99999) 
				{
/*					printf("from within PSA(2): since it became fractional, setting branch var to %d for seqnum %d\n",best_index,snum);*/
					branch_seqnum = snum;
/*					printf("changing frac index to %d (%d)\n",integer_indices[i],__LINE__);*/
/*					printf("plot(%lf,%lf,'co');\n",x_new[obj1_index],x_new[obj2_index]);*/
					frac_index = integer_indices[i];
					frac_val = x_new[integer_indices[i]];
					if(frac_scores[i] > 0.0001) 
					{
						frac_scores[i] += within_PSA_score;
/*						printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
					}
					else
					{
						frac_scores[i] = within_PSA_score + multiplier*integer_indices[i];
/*						printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
						num_frac++;
					}
					frac_values[integer_indices[i]] = frac_val;
/*					printf("num_frac: %d\n",num_frac);*/
					num_diff++;
/*					printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
				}
				else if( fabs(x[integer_indices[i]]-x_new[integer_indices[i]]) >= .00001)
				{
					branch_seqnum = snum;
/*					printf("changing frac index to %d (%d)\n",integer_indices[i],__LINE__);*/
/*					printf("plot(%lf,%lf,'mo');\n",x[obj1_index],x[obj2_index]);*/
/*					printf("plot(%lf,%lf,'mo');\n",x_new[obj1_index],x_new[obj2_index]);*/
					frac_index = integer_indices[i];
					frac_val = fmax(x[integer_indices[i]],x_new[integer_indices[i]])-.5;
					if(frac_scores[i] > 0.0001) 
					{
						frac_scores[i] += within_PSA_score;
/*						printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
					}
					else
					{
						frac_scores[i] = within_PSA_score + multiplier*integer_indices[i];
/*						printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
						num_frac++;
					}
					frac_values[integer_indices[i]] = frac_val;
/*					printf("num_frac: %d\n",num_frac);*/
					num_diff++;
/*					printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
				}
			}
			if(slope_scores_in_psa && num_diff > 1) 
			{
/*				printf("number of integer variables that changed: %d\n",num_diff);*/
				double *new_up = (double *) malloc (total_num_integer* sizeof (double));
				double *new_low = (double *) malloc (total_num_integer * sizeof (double));
				char *new_lu = (char *) malloc (total_num_integer * sizeof (char));
				double *slope_score = (double *) malloc (num_diff* sizeof (double));
				int counter = 0, changed_up_or_down = 0;
				
				for(i=0;i<total_num_integer;i++)
				{
					new_up[i] = x[integer_indices[i]];
					new_low[i] = x[integer_indices[i]];
				}
				
				CPXLPptr temp_lp = CPXcloneprob (env3, prob1, &status);
				
				for(i=0;i<total_num_integer;i++) new_lu[i] = 'L';
				status = CPXchgbds (env3, temp_lp, total_num_integer, integer_indices, new_lu, new_low);
		    		if ( status ){
		     			 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
		     			 goto TERMINATE;}
		     			 
		     		for(i=0;i<total_num_integer;i++) new_lu[i] = 'U';
		     		status = CPXchgbds (env3, temp_lp, total_num_integer, integer_indices, new_lu, new_up);
		    		if ( status ){
		     			 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
		     			 goto TERMINATE;}
				
				int j = 0;
				int change_indices[2] = {0,0};
				char change_lu[2] = {'L','U'};
				double change_bds[2] = {0.,0.};
				while(counter < num_diff)
				{
					for(i=j;i<total_num_integer;i++)
					{
/*						printf("i: %d, j: %d\n",i,j);*/
						if( x[integer_indices[i]] < x_new[integer_indices[i]] )
						{
							changed_up_or_down = 0;
							change_indices[0] = integer_indices[i];
							change_indices[1] = integer_indices[i];
							change_bds[0] = x[integer_indices[i]] + 1.;
							change_bds[1] = x[integer_indices[i]] + 1.;
							j = i+1;
							break;
						}
						else if( x[integer_indices[i]] > x_new[integer_indices[i]] )
						{
							changed_up_or_down = 1;
							change_indices[0] = integer_indices[i];
							change_indices[1] = integer_indices[i];
							change_bds[0] = x[integer_indices[i]] - 1.;
							change_bds[1] = x[integer_indices[i]] - 1.;
							j = i+1;
							break;
						}
					}
					
					status = CPXchgbds (env3, temp_lp, 2, change_indices, change_lu, change_bds);
			    		if ( status ){
			     			 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
			     			 goto TERMINATE;}
			     		
/*			     		status = CPXwriteprob (env3, temp_lp, "lp3.lp", "LP");*/
			     		
			     		status = CPXlpopt(env3, temp_lp);
				  	if ( status ) {
				    		printf ("Failed to optimize MIP, error code %d\n", status);
				    		goto TERMINATE;
				  	}
				  	int lpstat = CPXgetstat (env3, temp_lp);
				  	
/*				  	printf("status of solve: %d\n",lpstat);*/
				  	
				  	double objvals[2] = {0.,0.};
				  	double slope = 0.;
				  	if(lpstat != 3)
				  	{
				  		status = CPXgetx (env3, temp_lp, objvals, obj1_index, obj2_index);
/*				  		printf("plot([%lf],[%lf],'o');\n",objvals[0],objvals[1]);*/
				  		slope = (objvals[1]-x[obj2_index])/(objvals[0]-x[obj1_index]);
/*				  		printf("slope: %lf\n",slope);*/
				  	}
				  	if(slope < -99999) slope = -99999.;
				  	if(slope >= 0.) slope_score[counter] = j-1;
				  	else slope_score[counter] = j-1 - (1. - (slope - .00001)/slope);
				  	
				  	if(changed_up_or_down)
				  	{
				  		change_bds[0] += 1.;
				  		change_bds[1] += 1.;
				  	}
				  	else 
				  	{
				  		change_bds[1] -= 1.;
				  		change_bds[0] -= 1.;
				  	}
				  	status = CPXchgbds (env3, temp_lp, 2, change_indices, change_lu, change_bds);
			    		if ( status ){
			     			 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
			     			 goto TERMINATE;}
				  	
				  	counter++;
				}
				
/*				for(i=0;i<num_diff;i++) printf("slope_score%d: %lf\n",i,slope_score[i]);*/
		     		qsort(slope_score, num_diff, sizeof(double), comparison_function);
/*		     		for(i=0;i<num_diff;i++) printf("slope_score%d: %lf\n",i,slope_score[i]);*/
		     		
		     		int best_index_based_on_slope = (int) slope_score[0];	
		     		frac_scores[best_index_based_on_slope] += 2;
/*		     		printf("best index: %d\tfrac score: %lf\n",best_index_based_on_slope,frac_scores[best_index_based_on_slope]);*/
		     		
		     		if(slope_score) free_and_null ((char **) &slope_score);
		     		if(new_lu) free_and_null ((char **) &new_lu);
		     		if(new_up) free_and_null ((char **) &new_up);
		     		if(new_low) free_and_null ((char **) &new_low);
		     		
				goto TERMINATE;
			}
			if(num_diff) goto TERMINATE;
			
/*			printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x[obj1_index],x_new[obj1_index],x[obj2_index],x_new[obj2_index]);*/
			insert(2,x_new[obj1_index],x_new[obj2_index],x[obj1_index],x[obj2_index],
			      (x[obj1_index]-x_new[obj1_index])/(x[obj2_index]-x_new[obj2_index]),&tree,&empty_node);
			      
			for(i=0;i<obj2_index+1;i++) 
			{
				x[i] = x_new[i];
				if(generate_disjunctive_cuts_from_obj_space_disjunctions && PSA_left_iter_cnt == 1) temp_x_l[i] = x[i];
			}
			
			first_it = 0;
			sub_pr2_x_lb = x[obj1_index];
  			sub_pr2_y_ub = x[obj2_index];
		}
	}
	
	if(snum != 0 && same_seq != 1) retval = 2;
 	
 	TERMINATE:
 	
 	if(first_it != 1) 
 	{
/*		printf("changing x_ub to %lf and y_lb to %lf\n",sub_pr1_x_ub,sub_pr1_y_lb);*/
		int ind[2] = {obj1_index,obj2_index};
		char lu[2] = {'L','U'};
		double bd[2] = {sub_pr2_x_lb,sub_pr2_y_ub};
		
		if(retval != 2)
		{
			int add_check = mock_insert(1,x[obj1_index],x[obj2_index],0,0,0.,&tree);
			if(!add_check)
			{
/*				printf("from within PSA left, last found value is dominated\n");*/
/*				printf("plot([%lf],[%lf],'ko');\n",x[obj1_index],x[obj2_index]);*/
				if(x[obj1_index] != x_new[obj1_index])
				{
/*					printf("checking other soln\n");*/
					add_check = mock_insert(2,x_new[obj1_index],x_new[obj2_index],x[obj1_index],x[obj2_index],
						(x[obj1_index]-x_new[obj1_index])/(x[obj2_index]-x_new[obj2_index]),&tree);
				}
				if(!add_check)
				{
/*					printf("running PSA reduce left to see if the remainder of the bound is dominated\n");*/
/*					printf("plot([%lf],[%lf],'ko');\n",x_new[obj1_index],x_new[obj2_index]);*/
					retval = PSA_reduce_left(env, prob1, x, basis_col_info, basis_row_info, indices, snum);
/*					printf("value of the return: %d\n",retval);*/
				}
			}
		}

  		if(retval != 2) retval = 1;
	}
 	
	if(x) free_and_null ((char **) &x);
	if(basis_col_info) free_and_null ((char **) &basis_col_info);
	if(basis_row_info) free_and_null ((char **) &basis_row_info);
	if(indices) free_and_null ((char **) &indices);
	if(x_new) free_and_null ((char **) &x_new);
	CPXfreeprob(env, &prob1);
	CPXfreeprob(env, &prob2);

	return(retval);
}

/*********************************************************************************************************************** 

	This function is used for generating integer feasible line segments that are part of the Pareto
	set of a slice problem. (Works from right to left)
	
***********************************************************************************************************************/

void PSA_full_right(CPXCENVptr env, CPXLPptr prob, double *x_orig, int *basis_col_info_orig, int *basis_row_info_orig)
{
/*	print_on = 1;*/
/*	printf("PSA_full right\n");*/

	int iter_cnt = 0;
	
	if( !prob ) prob = lp_1;
	
	double first_obj1_val, first_obj2_val;
	
	/*************** Getting started ***************/
	
	int status;
	CPXLPptr prob1 = CPXcloneprob (env, prob, &status);
	
	/*************** Fix the integer variables ***************/
	
	int i;
	char *b;
	double *integer_x;
	int *indices;
	b = (char *) malloc (total_num_integer*sizeof(char));
	integer_x = (double *) malloc (total_num_integer*sizeof(double));
	indices = (int *) malloc (cur_numcols * sizeof (int));
			
	for(i=0;i<total_num_integer;i++)
	{
		b[i] = 'B';
		integer_x[i] = x_orig[integer_indices[i]];
	}
	
/*	printf("changing bounds 1\n");*/
	status = CPXchgbds (env, prob1, total_num_integer, integer_indices, b, integer_x);
	
	char up[1] = {'U'};
	
	status = CPXchgbds (env, prob1, 1, &obj1_index, up, &x_orig[obj1_index]);
	if ( status ){
		 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
		 goto TERMINATE;}
	
	/*************** Build prob 2 ***************/
	
	CPXLPptr prob2 = CPXcloneprob (env, prob1, &status);
  	for(i=0;i<cur_numcols;i++) indices[i] = i;
  	
  	chg_coefs(env,prob1,indices,-10000000.);
  	chg_coefs(env,prob2,indices,-.000000001);
	
	/*************** Solve LP1 ***************/
		
	status = CPXlpopt (env, prob1);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	int lpstat = CPXgetstat (env, prob1);
  	
  	if(lpstat == 12)
  	{
  		double upper_limit,lower_limit,what_it_should_be;
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_UpperObj, &upper_limit);
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_LowerObj, &lower_limit);
  		what_it_should_be = pow(10.,75.);
  		if(upper_limit < what_it_should_be || lower_limit > -what_it_should_be)
  		{
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_UpperObj, what_it_should_be);
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_LowerObj, -what_it_should_be);
  			printf("Warning: Within PSA_full encountered an issue where objective value limits were set to strange values. This was NOT a user error. Resetting these values to default and restarting PSA_full\n");
  			PSA_full_right(env,NULL,x_orig,NULL,NULL);
  			goto TERMINATE;
  		}
  	}
  	else if(lpstat == 3) goto TERMINATE;
  	
  	double prev_obj_vals[2] = {0.,0.};
  	double new_obj_vals[2] = {0.,0.};
  	double obj_vals[2] = {0.,0.};
  	double coef_lb = 0., coef_ub = 0.;
/*  	printf("the status of the solve: %d\n",lpstat);*/

	status = CPXgetx (env, prob1, prev_obj_vals, obj1_index, obj2_index);
	
/*	printf("plot([%lf],[%lf],'go');\n",prev_obj_vals[0],prev_obj_vals[1]);*/
	
	status = CPXlpopt (env, prob2);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	status = CPXgetx (env, prob2, obj_vals, obj1_index, obj2_index);
  	
  	if( fabs(obj_vals[0] - prev_obj_vals[0]) > .0001 && fabs(obj_vals[1] - prev_obj_vals[1]) > .0001 )
  	{
/*  		printf("obj1 and obj2 solution are different. Keep going\n");*/
/*  		printf("plot([%lf],[%lf],'go');\n",obj_vals[0],obj_vals[1]);*/
  	}
  	else
  	{
/*  		printf("obj1 and obj2 solutions are the same. Stop\n");*/
/*		printf("plot([%lf],[%lf],'go');\n",prev_obj_vals[0],prev_obj_vals[1]);*/
  		insert(1,prev_obj_vals[0],prev_obj_vals[1],0.,0.,0.,&tree,&empty_node);
  		goto TERMINATE;
  	}
	
	int same_cnt = 0;
	double prev_coef_lb = -100000000000000000., add_val = .001;
	while( fabs(prev_obj_vals[0] - obj_vals[0]) > .01 )
	{
		iter_cnt++;
		if(iter_cnt > 5000)
		{
			printf("%lf, %lf\n",prev_obj_vals[0], obj_vals[0]);
			printf("same count: %d\n",same_cnt);
/*			printf("%d\n",__LINE__);*/
			goto TERMINATE;
		}
/*		printf("prev_obj_vals[0]: %lf, obj_vals[0]: %lf\n",prev_obj_vals[0],obj_vals[0]);*/
		status = CPXobjsa(env,prob1,obj2_index,obj2_index,&coef_lb,&coef_ub);
/*		printf("bounds on coefficient: %lf to %lf\n",coef_lb,coef_ub);*/
		if(coef_ub > 10000000000.) break;
		if(fabs(prev_coef_lb - coef_lb) < .00000001)
		{
			coef_ub += add_val;
			add_val = add_val*2;
/*			retval = 0;*/
/*			goto TERMINATE;*/
		}
		prev_coef_lb = coef_lb;
		chg_coefs(env,prob1,indices, -1./(coef_ub + .00001));
		status = CPXlpopt (env, prob1);
		if ( status ) {
			printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
			goto TERMINATE;
		}
		lpstat = CPXgetstat (env, prob1);
/*	  	printf("the status of the solve: %d\n",lpstat);*/
		status = CPXgetx (env, prob1, new_obj_vals, obj1_index, obj2_index);
/*		printf("the next point might be: (%lf,%lf)\n",new_obj_vals[0],new_obj_vals[1]);*/
		
		if(new_obj_vals[0] < prev_obj_vals[0])
		{
/*			printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",prev_obj_vals[0],new_obj_vals[0],prev_obj_vals[1],new_obj_vals[1]);*/
			if(check_for_stopping_PSA_full)
			{
				int insert_check = mock_insert(2,prev_obj_vals[0],prev_obj_vals[1],new_obj_vals[0],new_obj_vals[1],
				(prev_obj_vals[1]-new_obj_vals[1])/(prev_obj_vals[0]-new_obj_vals[0]),&tree);
				if(!insert_check) 
				{
/*					printf("Stopping PSA full right early\n");	*/
					goto TERMINATE;
				}
			}
			insert(2,prev_obj_vals[0],prev_obj_vals[1],new_obj_vals[0],new_obj_vals[1],
				(prev_obj_vals[1]-new_obj_vals[1])/(prev_obj_vals[0]-new_obj_vals[0]),&tree,&empty_node);
			prev_obj_vals[0] = new_obj_vals[0];
			prev_obj_vals[1] = new_obj_vals[1];
			same_cnt = 0;
		}
		else same_cnt++;
		if( fabs(prev_obj_vals[0] - obj_vals[0]) < .1 && same_cnt >=5 ) break;
	}
 	
 	TERMINATE:
 	
 	free_and_null ((char **) &b);
	free_and_null ((char **) &integer_x);
	free_and_null ((char **) &indices);
	CPXfreeprob(env,&prob1);
	CPXfreeprob(env,&prob2);
}

/*********************************************************************************************************************** 

	This function is used for generating integer feasible line segments that are part of the Pareto
	set of a slice problem. (Works from left to right)
	
***********************************************************************************************************************/

void PSA_full_left(CPXCENVptr env, CPXLPptr prob, double *x_orig, int *basis_col_info_orig, int *basis_row_info_orig)
{
/*	print_on = 1;*/
/*	printf("PSA_full left\n");*/

	int iter_cnt = 0;
	
	if( !prob ) prob = lp_1;
	
	double first_obj1_val, first_obj2_val;
	
	/*************** Getting started ***************/
	
	int status;
	CPXLPptr prob1 = CPXcloneprob (env, prob, &status);
	
	/*************** Fix the integer variables ***************/
	
	int i;
	char *b;
	double *integer_x;
	int *indices;
	b = (char *) malloc (total_num_integer*sizeof(char));
	integer_x = (double *) malloc (total_num_integer*sizeof(double));
	indices = (int *) malloc (cur_numcols * sizeof (int));
	
/*	printf("_________----------___________\n");*/
/*	for(i=0;i<cur_numcols;i++) printf("x_%d: %lf\n",i,x_orig[i]);*/
/*	printf("_________----------___________\n");*/
			
	for(i=0;i<total_num_integer;i++)
	{
		b[i] = 'B';
		integer_x[i] = x_orig[integer_indices[i]];
/*		printf("x%d: %lf\n",integer_indices[i],x_orig[integer_indices[i]]);*/
	}
	
/*	printf("changing bounds 1\n");*/
	status = CPXchgbds (env, prob1, total_num_integer, integer_indices, b, integer_x);
	
	char low[1] = {'L'};
	
	status = CPXchgbds (env, prob1, 1, &obj1_index, low, &x_orig[obj1_index]);
	if ( status ){
		 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
		 goto TERMINATE;}
	
	/*************** Build prob 2 ***************/
	
	CPXLPptr prob2 = CPXcloneprob (env, prob1, &status);
  	for(i=0;i<cur_numcols;i++) indices[i] = i;
  	
  	chg_coefs(env,prob1,indices,-.0000001);
  	chg_coefs(env,prob2,indices,-1000000000.);
	
	/*************** Solve LP1 ***************/
		
	status = CPXlpopt (env, prob1);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	int lpstat = CPXgetstat (env, prob1);
  	
  	if(lpstat == 12)
  	{
  		double upper_limit,lower_limit,what_it_should_be;
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_UpperObj, &upper_limit);
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_LowerObj, &lower_limit);
  		what_it_should_be = pow(10.,75.);
  		if(upper_limit < what_it_should_be || lower_limit > -what_it_should_be)
  		{
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_UpperObj, what_it_should_be);
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_LowerObj, -what_it_should_be);
  			printf("Warning: Within PSA_full encountered an issue where objective value limits were set to strange values. This was NOT a user error. Resetting these values to default and restarting PSA_full\n");
  			PSA_full_left(env,NULL,x_orig,NULL,NULL);
  			goto TERMINATE;
  		}
  	}
  	else if(lpstat == 3) goto TERMINATE;
  	
  	double prev_obj_vals[2] = {0.,0.};
  	double new_obj_vals[2] = {0.,0.};
  	double obj_vals[2] = {0.,0.};
  	double coef_lb = 0., coef_ub = 0.;
/*  	printf("the status of the solve: %d\n",lpstat);*/

	status = CPXgetx (env, prob1, prev_obj_vals, obj1_index, obj2_index);
	
/*	printf("plot([%lf],[%lf],'go');\n",prev_obj_vals[0],prev_obj_vals[1]);*/
	
	status = CPXlpopt (env, prob2);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	status = CPXgetx (env, prob2, obj_vals, obj1_index, obj2_index);
  	
  	if( fabs(obj_vals[0] - prev_obj_vals[0]) > .0001 && fabs(obj_vals[1] - prev_obj_vals[1]) > .0001 )
  	{
/*  		printf("obj1 and obj2 solution are different. Keep going\n");*/
/*  		printf("plot([%lf],[%lf],'go');\n",obj_vals[0],obj_vals[1]);*/
  	}
  	else
  	{
/*  		printf("obj1 and obj2 solutions are the same. Stop\n");*/
/*		printf("plot([%lf],[%lf],'go');\n",prev_obj_vals[0],prev_obj_vals[1]);*/
  		insert(1,prev_obj_vals[0],prev_obj_vals[1],0.,0.,0.,&tree,&empty_node);
  		goto TERMINATE;
  	}
	
	int same_cnt = 0;
	double prev_coef_ub = 100000000000000000., add_val = .001;
	while( fabs(prev_obj_vals[0] - obj_vals[0]) > .01 )
	{
		iter_cnt++;
		if(iter_cnt > 5000)
		{
			printf("%lf, %lf\n",prev_obj_vals[0], obj_vals[0]);
			printf("same count: %d\n",same_cnt);
/*			printf("it_cnt: %d, (%d)\n",iter_cnt,__LINE__);*/
			goto TERMINATE;
		}
/*		printf("prev_obj_vals[0]: %lf, obj_vals[0]: %lf\n",prev_obj_vals[0],obj_vals[0]);*/
		status = CPXobjsa(env,prob1,obj2_index,obj2_index,&coef_lb,&coef_ub);
/*		printf("bounds on coefficient: %lf to %lf\n",coef_lb,coef_ub);*/
		if(coef_lb < -10000000000.) break;
		if(fabs(prev_coef_ub - coef_ub) < .00000001)
		{
			coef_lb -= add_val;
			add_val = add_val*2;
/*			retval = 0;*/
/*			goto TERMINATE;*/
		}
		prev_coef_ub = coef_ub;
		chg_coefs(env,prob1,indices, -1./(coef_lb - .00001));
		status = CPXlpopt (env, prob1);
		if ( status ) {
			printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
			goto TERMINATE;
		}
		lpstat = CPXgetstat (env, prob1);
/*	  	printf("the status of the solve: %d\n",lpstat);*/
		status = CPXgetx (env, prob1, new_obj_vals, obj1_index, obj2_index);
/*		printf("the next point might be: (%lf,%lf)\n",new_obj_vals[0],new_obj_vals[1]);*/
		
		if(new_obj_vals[0] > prev_obj_vals[0])
		{
/*			printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",prev_obj_vals[0],new_obj_vals[0],prev_obj_vals[1],new_obj_vals[1]);*/
			if(check_for_stopping_PSA_full)
			{
				int insert_check = mock_insert(2,new_obj_vals[0],new_obj_vals[1],prev_obj_vals[0],prev_obj_vals[1],
					(prev_obj_vals[1]-new_obj_vals[1])/(prev_obj_vals[0]-new_obj_vals[0]),&tree);
				if(!insert_check) 
				{
/*					printf("Stopping PSA full left early\n");*/
					goto TERMINATE;
				}
			}
			insert(2,new_obj_vals[0],new_obj_vals[1],prev_obj_vals[0],prev_obj_vals[1],
				(prev_obj_vals[1]-new_obj_vals[1])/(prev_obj_vals[0]-new_obj_vals[0]),&tree,&empty_node);
			prev_obj_vals[0] = new_obj_vals[0];
			prev_obj_vals[1] = new_obj_vals[1];
			same_cnt = 0;
		}
		else same_cnt++;
		if( fabs(prev_obj_vals[0] - obj_vals[0]) < .1 && same_cnt >=5 ) break;
	}
 	
 	TERMINATE:
 	
 	free_and_null ((char **) &b);
	free_and_null ((char **) &integer_x);
	free_and_null ((char **) &indices);
	CPXfreeprob(env,&prob1);
	CPXfreeprob(env,&prob2);
}

void PSA_full(CPXCENVptr env, CPXLPptr prob, double *x_orig, int *basis_col_info_orig, int *basis_row_info_orig)
{
/*	printf("plot([%lf],[%lf],'go');\n",x_orig[obj1_index],x_orig[obj2_index]);*/
/*	printf("calling PSA full right\n");*/
	PSA_full_right(env, prob, x_orig, basis_col_info_orig, basis_row_info_orig);
/*	printf("calling PSA full left\n");*/
	PSA_full_left(env, prob, x_orig, basis_col_info_orig, basis_row_info_orig);
}

/*********************************************************************************************************************** 

	This function is unused. It was originally used for generating integer feasible line segments that are 
	part of the Pareto set of a slice problem.
	
***********************************************************************************************************************/

void PSA_full2(CPXCENVptr env, CPXLPptr prob, double *x_orig, int *basis_col_info_orig, int *basis_row_info_orig)
{
/*	print_on = 1;*/
/*	printf("PSA_full\n");*/

	int iter_cnt = 0;
	
	if( !prob ) prob = lp_1;
	
	double first_obj1_val, first_obj2_val;
	
	/*************** Getting started ***************/
	
	int status;
	CPXLPptr prob1 = CPXcloneprob (env, prob, &status);
	
	/*************** Fix the integer variables ***************/
	
	int i;
	char *b;
	double *integer_x;
	int *indices;
	b = (char *) malloc (total_num_integer*sizeof(char));
	integer_x = (double *) malloc (total_num_integer*sizeof(double));
	indices = (int *) malloc (cur_numcols * sizeof (int));
			
	for(i=0;i<total_num_integer;i++)
	{
		b[i] = 'B';
		integer_x[i] = x_orig[integer_indices[i]];
	}
	
/*	printf("changing bounds 1\n");*/
	status = CPXchgbds (env, prob1, total_num_integer, integer_indices, b, integer_x);
	
	/*************** Build prob 2 ***************/
	
	CPXLPptr prob2 = CPXcloneprob (env, prob1, &status);
  	for(i=0;i<cur_numcols;i++) indices[i] = i;
  	
  	chg_coefs(env,prob1,indices,-10000000.);
  	chg_coefs(env,prob2,indices,-.000000001);
/*  	status = CPXchgobj (env, prob2, numcols, indices, ob_coef2);*/
/*	if ( status ) {*/
/*		printf ("Failed to get change obj coef. Error code %d\n", status);*/
/*		goto TERMINATE;*/
/*	}*/
	
/*	status = CPXwriteprob (env, prob1, "myprob1.lp", "LP");*/
/*	status = CPXwriteprob (env, prob2, "myprob2.lp", "LP");*/
/*	exit(0);*/

/*	printf("problem 1 obj2_index coef: %lf\n",1./1000.);*/
	
	/*************** Solve LP1 ***************/
		
	status = CPXlpopt (env, prob1);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	int lpstat = CPXgetstat (env, prob1);
  	
  	if(lpstat == 12)
  	{
  		double upper_limit,lower_limit,what_it_should_be;
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_UpperObj, &upper_limit);
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_LowerObj, &lower_limit);
  		what_it_should_be = pow(10.,75.);
  		if(upper_limit < what_it_should_be || lower_limit > -what_it_should_be)
  		{
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_UpperObj, what_it_should_be);
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_LowerObj, -what_it_should_be);
  			printf("Warning: Within PSA_full encountered an issue where objective value limits were set to strange values. This was NOT a user error. Resetting these values to default and restarting PSA_full\n");
  			PSA_full(env,NULL,x_orig,NULL,NULL);
  			goto TERMINATE;
  		}
  	}
  	
  	double prev_obj_vals[2] = {0.,0.};
  	double new_obj_vals[2] = {0.,0.};
  	double obj_vals[2] = {0.,0.};
  	double coef_lb = 0., coef_ub = 0.;
/*  	printf("the status of the solve: %d\n",lpstat);*/

	status = CPXgetx (env, prob1, prev_obj_vals, obj1_index, obj2_index);
	
/*	printf("plot([%lf],[%lf],'go');\n",prev_obj_vals[0],prev_obj_vals[1]);*/
	
	status = CPXlpopt (env, prob2);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	status = CPXgetx (env, prob2, obj_vals, obj1_index, obj2_index);
  	
  	if( fabs(obj_vals[0] - prev_obj_vals[0]) > .0001 && fabs(obj_vals[1] - prev_obj_vals[1]) > .0001 )
  	{
/*  		printf("obj1 and obj2 solution are different. Keep going\n");*/
/*  		printf("plot([%lf],[%lf],'go');\n",obj_vals[0],obj_vals[1]);*/
  	}
  	else
  	{
/*  		printf("obj1 and obj2 solutions are the same. Stop\n");*/
/*		printf("plot([%lf],[%lf],'go');\n",prev_obj_vals[0],prev_obj_vals[1]);*/
  		insert(1,prev_obj_vals[0],prev_obj_vals[1],0.,0.,0.,&tree,&empty_node);
  		goto TERMINATE;
  	}
	
	int same_cnt = 0;
	while( fabs(prev_obj_vals[0] - obj_vals[0]) > .01 )
	{
		iter_cnt++;
		if(iter_cnt > 5000)
		{
			printf("%lf, %lf\n",prev_obj_vals[0], obj_vals[0]);
			printf("%d\n",__LINE__);
			exit(0);
		}
/*		printf("prev_obj_vals[0]: %lf, obj_vals[0]: %lf\n",prev_obj_vals[0],obj_vals[0]);*/
		status = CPXobjsa(env,prob1,obj2_index,obj2_index,&coef_lb,&coef_ub);
/*		printf("bounds on coefficient: %lf to %lf\n",coef_lb,coef_ub);*/
		chg_coefs(env,prob1,indices, -1./(coef_ub + .00001));
		status = CPXlpopt (env, prob1);
		if ( status ) {
			printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
			goto TERMINATE;
		}
		lpstat = CPXgetstat (env, prob1);
/*	  	printf("the status of the solve: %d\n",lpstat);*/
		status = CPXgetx (env, prob1, new_obj_vals, obj1_index, obj2_index);
/*		printf("the next point might be: (%lf,%lf)\n",new_obj_vals[0],new_obj_vals[1]);*/
		
		if(new_obj_vals[0] < prev_obj_vals[0])
		{
/*			printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",prev_obj_vals[0],new_obj_vals[0],prev_obj_vals[1],new_obj_vals[1]);*/
			insert(2,prev_obj_vals[0],prev_obj_vals[1],new_obj_vals[0],new_obj_vals[1],
				(prev_obj_vals[1]-new_obj_vals[1])/(prev_obj_vals[0]-new_obj_vals[0]),&tree,&empty_node);
			prev_obj_vals[0] = new_obj_vals[0];
			prev_obj_vals[1] = new_obj_vals[1];
			same_cnt = 0;
		}
		else same_cnt++;
		if( fabs(prev_obj_vals[0] - obj_vals[0]) < .1 && same_cnt >=5 ) break;
	}
 	
 	TERMINATE:
 	
/* 	printf("PSA_full, itcnt: %d\n",iter_cnt);*/
 	
/* 	status = CPXsetintparam ( (CPXENVptr) env, CPX_PARAM_ITLIM, 2100000000);*/
/*  	if ( status ) {*/
/*    		printf ("%s(%d): CPXsetintparam, Failed to set iterations limit of SIMPLEX, error code %d\n", __FILE__, __LINE__, status);*/
/*    		goto TERMINATE;*/
/*  	}*/
 	
 	free_and_null ((char **) &b);
	free_and_null ((char **) &integer_x);
	free_and_null ((char **) &indices);
	CPXfreeprob(env,&prob1);
	CPXfreeprob(env,&prob2);
}

/*********************************************************************************************************************** 

	This function is used for generating line segments that are portions of a node's dual bound and are dominated
	by the primal bound. Stops when a segment has been generated that is not dominated. (Works from right to left)
	
***********************************************************************************************************************/
int PSA_rr_iter_cnt = 0, PSA_rl_iter_cnt = 0;
int PSA_reduce_right(CPXCENVptr env, CPXLPptr prob, double *x_orig, int *basis_col_info_orig, int *basis_row_info_orig, int *indices, int snum)
{	
/*	printf("reducing right side\n");*/
	PSA_st = clock();
	PSA_rr_iter_cnt = 0;
	double ti = 0.;
	
	if(generate_disjunctive_cuts_from_obj_space_disjunctions && !temp_x_r) temp_x_r = (double *) malloc (cur_numcols*sizeof(double));

	/*************** Solve LP1 ***************/
	int ret_val = 0;
	int first_it = 1;
	int status = 0;
	
  	double	 *x = NULL;
  	int 	 *basis_col_info = NULL;
  	int 	 *basis_row_info = NULL;
  	
  	x = (double *) malloc (cur_numcols * sizeof (double));
  	double *x_new = (double *) malloc (cur_numcols * sizeof (double));
  	double coef_lb = 0., coef_ub = 0.;
	
	int i = 0;
	for(i=0;i<cur_numcols;i++)
	{
		x[i] = x_orig[i];
	}

	sub_pr1_x_ub = x[obj1_index];
	sub_pr1_y_lb = x[obj2_index];
  	
  	CPXLPptr prob1 = CPXcloneprob (env, prob, &status);
  	CPXchgprobtype(env, prob1, CPXPROB_LP);	
	CPXLPptr prob2 = CPXcloneprob (env, prob1, &status);
	
	char up[1] = {'U'};
	
	status = CPXchgbds (env, prob1, 1, &obj1_index, up, &x_orig[obj1_index]);
	if ( status ){
		 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
		 goto TERMINATE;}
  	
  	chg_coefs(env,prob1,indices,-10000000.);
  	chg_coefs(env,prob2,indices,-.000000001);
	
	/*************** Get info from LP2 ***************/
	
	status = CPXcopybase (env, prob1, basis_col_info_orig, basis_row_info_orig);
  	if ( status ) {
    		printf ("%s(%d): CPXcopybase, Failed to copy basis, error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	status = CPXlpopt (env, prob1);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	status = CPXcopybase (env, prob2, basis_col_info_orig, basis_row_info_orig);
  	if ( status ) {
    		printf ("%s(%d): CPXcopybase, Failed to copy basis, error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}

  	int lpstat = CPXgetstat (env, prob1);
  	
  	if(lpstat == 3)
  	{
  		ret_val = 2;
  		goto TERMINATE;
  	}
  	
  	status = CPXgetx (env, prob1, x, 0, cur_numcols-1);
  	if ( status ) {
   		printf ("%s(%d): CPXgetx, Failed to get x-values, error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
	
/*	printf("plot([%lf],[%lf],'go');\n",x[obj1_index],x[obj2_index]);*/
	
	status = CPXlpopt (env, prob2);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	double prob2_obvals[2] = {0.,0.};
  	
  	status = CPXgetx (env, prob2, prob2_obvals, obj1_index, obj2_index);
  	
  	if( fabs(x[obj1_index] - prob2_obvals[0]) > .0001 && fabs(x[obj2_index] - prob2_obvals[1]) > .0001 )
  	{
/*  		printf("obj1 and obj2 solution are different. Keep going\n");*/
/*  		printf("plot([%lf],[%lf],'go');\n",prob2_obvals[0],prob2_obvals[1]);*/
  	}
  	else
  	{
/*  		printf("obj1 and obj2 solutions are the same. Stop\n");*/
		ret_val = 2;
  		goto TERMINATE;
  	}
  	
  	double overall_coef_lb = 0.;
  	
  	for(i=0;i<cur_numcols;i++) x_new[i] = x[i];
	
	int same_cnt = 0;
	while( fabs(prob2_obvals[0] - x_new[obj1_index]) > .01 )
	{
		PSA_rr_iter_cnt++;
		PSA_fi = clock();
		ti = (double)(PSA_fi - PSA_st) / CLOCKS_PER_SEC;
		if( (ti > .5 && same_cnt > 25) || (PSA_rr_iter_cnt > 100 && same_cnt >= PSA_rr_iter_cnt -5))//if(PSA_rr_iter_cnt > 1000)
		{
			ti = (double)(PSA_fi - PSA_st) / CLOCKS_PER_SEC;
			printf("time in PSA reduce right: %lf\n",ti);
			printf("iterations: %d\n",PSA_rr_iter_cnt);
			printf("same count: %d\n",same_cnt);
			printf("%lf, %lf\n",prob2_obvals[0],x_new[obj1_index]);
			printf("%d\n",__LINE__);
			if(coef_add < .01) coef_add = coef_add*5.;
			goto TERMINATE;
		}
		status = CPXobjsa(env,prob1,obj2_index,obj2_index,&coef_lb,&coef_ub);
/*		printf("status: %d\n",status);*/
/*		printf("bounds on coefficient: %lf to %lf\n",coef_lb,coef_ub);*/
		if(coef_lb > 10000000000.)
		{
			ret_val = 2;
			goto TERMINATE;
		}
/*		printf("new coef: %lf\n",-1./(-1./(coef_ub + coef_add)));*/
		chg_coefs(env,prob1,indices, -1./(coef_ub + coef_add));
		status = CPXlpopt (env, prob1);
		if ( status ) {
			printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
			goto TERMINATE;
		}
		lpstat = CPXgetstat (env, prob1);
/*	  	printf("the status of the solve: %d\n",lpstat);*/
		status = CPXgetx (env, prob1, x_new, 0, cur_numcols-1);
/*		printf("the next point might be: (%lf,%lf)\n",new_obj_vals[0],new_obj_vals[1]);*/
		
		if(x_new[obj1_index] - x[obj1_index] < -.000001 )
		{
			overall_coef_lb = coef_lb;
/*			printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",x[obj1_index],x_new[obj1_index],x[obj2_index],x_new[obj2_index]);*/
			int keep_check = mock_insert(2,x[obj1_index],x[obj2_index],x_new[obj1_index],x_new[obj2_index],
				      			(x[obj1_index]-x_new[obj1_index])/(x[obj2_index]-x_new[obj2_index]),&tree);
/*			exit(0);*/
			if(keep_check)
			{
/*					printf("keep! Exiting ...\n");*/
				if(first_it == 1) goto TERMINATE;

				ret_val = 1;
/*					printf("num integer: %d\n",num_integer);*/
				for(i=0;i<num_integer;i++)
				{
/*						printf("x_%d: %lf\n",integer_indices[i],x[integer_indices[i]]);*/
					double diff = x[integer_indices[i]] - floor(x[integer_indices[i]]);
					int any_diff = 0;
					if( diff >= .00001 && diff <= .99999) 
					{
/*  						printf("changing within PSA reduce right\n");*/
						branch_seqnum = snum;
/*  						printf("changing frac index to %d (%d)\n",integer_indices[i],__LINE__);*/
/*  						printf("plot(%lf,%lf,'co');\n",x[obj1_index],x[obj2_index]);*/
						frac_index = integer_indices[i];
						frac_val = x[integer_indices[i]];
						if(frac_scores[i] > 0.0001) 
						{
							frac_scores[i] += within_PSA_score;
/*							printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
						}
						else
						{
							frac_scores[i] = within_PSA_score + multiplier*integer_indices[i];
/*							printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
							num_frac++;
						}
						frac_values[integer_indices[i]] = frac_val;
/*						printf("num_frac: %d\n",num_frac);*/
						any_diff = 1;
/*						printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
					}
					if(any_diff) goto TERMINATE;
				}
				goto TERMINATE;
			}
			else 
			{
				sub_pr1_x_ub = x_new[obj1_index];
				sub_pr1_y_lb = x_new[obj2_index];
				first_it = 0;
				if(check_for_early_PSA_reduce || (check_for_early_PSA_reduce_after_n_iter && PSA_rr_iter_cnt >= n_iter))
				{
					double far_y = (x[obj2_index]-x_new[obj2_index])/(x[obj1_index]-x_new[obj1_index])*(prob2_obvals[0]-x[obj1_index])+
						x[obj2_index];
					keep_check = mock_insert(2,x[obj1_index],x[obj2_index],prob2_obvals[0],far_y,
				      			(x[obj1_index]-x_new[obj1_index])/(x[obj2_index]-x_new[obj2_index]),&tree);
/*				      	printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",x[obj1_index],prob2_obvals[0],x[obj2_index],far_y);*/
				      	if(!keep_check) 
				      	{
/*				      		printf("Could exit PSA reduce right after iteration %d\n",iter_cnt);*/
						ret_val = 2;
						goto TERMINATE;
					}
				}
			}
			same_cnt = 0;
		}
		else same_cnt++;
		
		for(i=0;i<cur_numcols;i++) 
		{
			x[i] = x_new[i];
			if(generate_disjunctive_cuts_from_obj_space_disjunctions && PSA_rr_iter_cnt == 1) temp_x_r[i] = x[i];
		}
		
		if( fabs(prob2_obvals[0] - x_new[obj1_index]) < .1 && same_cnt >=5 ) break;
	}
  	
  	ret_val = 2;
 	
 	
 	TERMINATE:
 	
/* 	PSA_fi = clock();*/
/* 	ti = (double)(PSA_fi - PSA_st) / CLOCKS_PER_SEC;*/
/*	printf("time in PSA reduce right: %lf\n",ti);*/
 	
  	CPXfreeprob(env,&prob1);
  	CPXfreeprob(env,&prob2);
 	
  	free_and_null ((char **) &x);
  	free_and_null ((char **) &x_new);

	return ret_val;
}

int global_snum = 0;

/*********************************************************************************************************************** 

	This function is used for generating line segments that are portions of a node's dual bound and are dominated
	by the primal bound. Stops when a segment has been generated that is not dominated. (Works from left to right)
	
***********************************************************************************************************************/

int PSA_reduce_left(CPXCENVptr env, CPXLPptr prob, double *x_orig, int *basis_col_info_orig, int *basis_row_info_orig, int *indices)
{	
/*	printf("reducing left side\n");*/
	PSA_st = clock();
	PSA_rl_iter_cnt = 0;
	int retval = 0;
	double ti = 0.;
	
	if(generate_disjunctive_cuts_from_obj_space_disjunctions && !temp_x_l) temp_x_l = (double *) malloc (cur_numcols*sizeof(double));
	
	/*************** Solve LP1 ***************/
	int first_it = 1;
	int status = 0;
	
  	double	 *x = (double *) malloc (cur_numcols * sizeof (double));
  	double	 *x_new = (double *) malloc (cur_numcols * sizeof (double));
  	double coef_lb = 0., coef_ub = 0.;
  	
  	int i = 0;
  	for(i=0;i<cur_numcols;i++)
	{
		x[i] = x_orig[i];
	}

	sub_pr2_x_lb = x[obj1_index];
	sub_pr2_y_ub = x[obj2_index];
	
/*	printf("within reduce left, new bounds: x_lb: %lf, y_ub: %lf\n",sub_pr2_x_lb,sub_pr2_y_ub);*/
  	
  	CPXLPptr prob1 = CPXcloneprob (env, prob, &status);	
  	CPXchgprobtype(env, prob1, CPXPROB_LP);
	CPXLPptr prob2 = CPXcloneprob (env, prob1, &status);
	
	char up[1] = {'U'};
	
	status = CPXchgbds (env, prob1, 1, &obj2_index, up, &x_orig[obj2_index]);
	if ( status ){
		 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
		 goto TERMINATE;}
		 
	chg_coefs(env,prob2,indices,-10000000.);
  	chg_coefs(env,prob1,indices,-.000000001);
  	
	
	status = CPXcopybase (env, prob1, basis_col_info_orig, basis_row_info_orig);
  	if ( status ) {
    		printf ("%s(%d): CPXcopybase, Failed to copy basis, error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	status = CPXlpopt (env, prob1);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	int lpstat = CPXgetstat (env, prob1);
  	
  	if(lpstat == 3)
  	{
  		retval = 2;
  		goto TERMINATE;
  	}
  	
  	status = CPXgetx (env, prob1, x, 0, cur_numcols-1);
  	if ( status ) {
   		printf ("%s(%d): CPXgetx, Failed to get x-values, error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	
  	for(i=0;i<cur_numcols;i++)
	{
		x_new[i] = x[i];
	}

	
/*	printf("plot([%lf],[%lf],'go');\n",x[obj1_index],x[obj2_index]);*/
  	
  	status = CPXcopybase (env, prob2, basis_col_info_orig, basis_row_info_orig);
  	if ( status ) {
    		printf ("%s(%d): CPXcopybase, Failed to copy basis, error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}
  	status = CPXlpopt (env, prob2);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		goto TERMINATE;
  	}

	double prob2_obvals[2] = {0.,0.};
  	
  	status = CPXgetx (env, prob2, prob2_obvals, obj1_index, obj2_index);
  	
  	if( fabs(x[obj1_index] - prob2_obvals[0]) > .0001 && fabs(x[obj2_index] - prob2_obvals[1]) > .0001 )
  	{
/*  		printf("obj1 and obj2 solution are different. Keep going\n");*/
/*  		printf("plot([%lf],[%lf],'go');\n",prob2_obvals[0],prob2_obvals[1]);*/
  	}
  	else
  	{
/*  		printf("obj1 and obj2 solutions are the same. Stop\n");*/
		retval = 2;
  		goto TERMINATE;
  	}
  	
  	double overall_coef_ub = 0.;
  	
  	for(i=0;i<cur_numcols;i++) x_new[i] = x[i];
	
	int same_cnt = 0;
	while( fabs(prob2_obvals[0] - x_new[obj1_index]) > .01 )
	{
		PSA_rl_iter_cnt++;
		PSA_fi = clock();
		ti = (double)(PSA_fi - PSA_st) / CLOCKS_PER_SEC;
		if( (ti > .5 && same_cnt > 25) || (PSA_rl_iter_cnt > 100 && same_cnt >= PSA_rl_iter_cnt -5))//if(PSA_rl_iter_cnt > 1000)
		{
			ti = (double)(PSA_fi - PSA_st) / CLOCKS_PER_SEC;
			printf("time in PSA reduce left: %lf\n",ti);
			printf("iterations: %d\n",PSA_rl_iter_cnt);
			printf("same count: %d\n",same_cnt);
			printf("%lf, %lf\n",prob2_obvals[0],x_new[obj1_index]);
			printf("%d\n",__LINE__);
			if(coef_add < .01) coef_add = coef_add*5.;
			goto TERMINATE;
		}
		status = CPXobjsa(env,prob1,obj2_index,obj2_index,&coef_lb,&coef_ub);
/*		printf("bounds on coefficient: %lf to %lf\n",coef_lb,coef_ub);*/
		if(coef_lb < -10000000000.)
		{
			retval = 2;
			goto TERMINATE;
		}
/*		printf("new coef: %lf\n",-1./(-1./(coef_lb - coef_add)));*/
		chg_coefs(env,prob1,indices, -1./(coef_lb - coef_add));
		status = CPXlpopt (env, prob1);
		if ( status ) {
			printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
			goto TERMINATE;
		}
		lpstat = CPXgetstat (env, prob1);
/*	  	printf("the status of the solve: %d\n",lpstat);*/
		status = CPXgetx (env, prob1, x_new, 0, cur_numcols-1);
/*		printf("the next point might be: (%lf,%lf)\n",new_obj_vals[0],new_obj_vals[1]);*/
		
		if(x[obj1_index] - x_new[obj1_index] < -.000001 )
		{
			overall_coef_ub = coef_ub;
/*			printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",x[obj1_index],x_new[obj1_index],x[obj2_index],x_new[obj2_index]);*/
			int keep_check = mock_insert(2,x_new[obj1_index],x_new[obj2_index],x[obj1_index],x[obj2_index],
				      			(x[obj1_index]-x_new[obj1_index])/(x[obj2_index]-x_new[obj2_index]),&tree);
			if(keep_check)
			{
/*				printf("keep! Exiting ...\n");*/
				if(first_it == 1) goto TERMINATE;

				retval = 1;
/*				printf("num integer: %d\n",num_integer);*/
				for(i=0;i<num_integer;i++)
				{
/*					printf("x_%d: %lf\n",integer_indices[i],x[integer_indices[i]]);*/
					double diff = x[integer_indices[i]] - floor(x[integer_indices[i]]);
					int any_diff = 0;
					if( diff >= .00001 && diff <= .99999) 
					{
/*  						printf("changing frac index to %d (%d)\n",integer_indices[i],__LINE__);*/
/*  						printf("plot(%lf,%lf,'co');\n",x[obj1_index],x[obj2_index]);*/
						frac_index = integer_indices[i];
						frac_val = x[integer_indices[i]];
						if(frac_scores[i] > 0.0001) 
						{
							frac_scores[i] += within_PSA_score;
/*							printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
						}
						else
						{
							frac_scores[i] = within_PSA_score + multiplier*integer_indices[i];
/*							printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
							num_frac++;
						}
						frac_values[integer_indices[i]] = frac_val;
/*						printf("num_frac: %d\n",num_frac);*/
						any_diff = 1;
/*						printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
					}
					if(any_diff) goto TERMINATE;
				}
				goto TERMINATE;
			}
			else 
			{
				sub_pr2_x_lb = x_new[obj1_index];
				sub_pr2_y_ub = x_new[obj2_index];
				first_it = 0;
				if(check_for_early_PSA_reduce || (check_for_early_PSA_reduce_after_n_iter && PSA_rl_iter_cnt >= n_iter))
				{
					double far_y = (x[obj2_index]-x_new[obj2_index])/(x[obj1_index]-x_new[obj1_index])*(prob2_obvals[0]-x[obj1_index])+
						x[obj2_index];
					keep_check = mock_insert(2,prob2_obvals[0],far_y,x[obj1_index],x[obj2_index],
				      			(x[obj1_index]-x_new[obj1_index])/(x[obj2_index]-x_new[obj2_index]),&tree);
/*				      	printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",x[obj1_index],prob2_obvals[0],x[obj2_index],far_y);*/
				      	if(!keep_check) 
				      	{
/*				      		printf("Could exit PSA reduce left after iteration %d\n",iter_cnt);*/
				      		retval = 2;
				      		goto TERMINATE;
				      	}
				}
			}
			same_cnt = 0;
		}
		else same_cnt++;
		
		for(i=0;i<cur_numcols;i++) 
		{	
			x[i] = x_new[i];
			if(generate_disjunctive_cuts_from_obj_space_disjunctions && PSA_rl_iter_cnt == 1) temp_x_l[i] = x[i];
		}
		
		if( fabs(prob2_obvals[0] - x_new[obj1_index]) < .1 && same_cnt >=5 ) break;
	}
  	
  	retval = 2;
 	
 	TERMINATE:
 	
/* 	PSA_fi = clock();*/
/* 	ti = (double)(PSA_fi - PSA_st) / CLOCKS_PER_SEC;*/
/*	printf("time in PSA reduce left: %lf\n",ti);*/
 	
  	CPXfreeprob(env,&prob1);
  	CPXfreeprob(env,&prob2);
	
  	free_and_null ((char **) &x);
  	free_and_null ((char **) &x_new);

	return retval;
}

/*********************************************************************************************************************** 

	This function is unused.
	
***********************************************************************************************************************/

int PSA_reduce(CPXCENVptr env, CPXLPptr prob, int *indices, int snum)
{
/*	printf("PSA reduce, seqnum: %d\n",snum);*/
	int ret_val = 0;
	/*************** Getting started ***************/
	
	int status;
	CPXLPptr prob1 = CPXcloneprob (env, prob, &status);
	CPXLPptr prob2 = CPXcloneprob (env, prob, &status);
	
	/*************** Build probs ***************/
  	
  	status = CPXchgobj (env, prob1, obj1_index+2, indices, ob_coef1);
	if ( status ) {
		printf ("(%d) Failed to get change obj coef. Error code %d\n", __LINE__,status);
		goto TERMINATE;
	}
  	status = CPXchgobj (env, prob2, obj1_index+2, indices, ob_coef2);
	if ( status ) {
		printf ("(%d) Failed to get change obj coef. Error code %d\n", __LINE__,status);
		goto TERMINATE;
	}
 	
 	TERMINATE:
 	
 	CPXfreeprob(env,&prob1);
	CPXfreeprob(env,&prob2);
 	
 	return ret_val;
}

int first_time_printing_starting_number_int_feas = 0;
int turn_cplex_cuts_off = 1;
int ran_hueristic = 0;

/*********************************************************************************************************************** 

	This function is the genetic algorithm used for generating integer feasible solutions during preprocessing
	when general integer variables are involved in the problem.
	
***********************************************************************************************************************/

int mixed_heuristic(CPXCENVptr env)
{
/*	printf("in mixed heur\n");*/
	int status = 0;
	CPXLPptr temp_lp = CPXcloneprob (env, lp_2, &status);
	double ran1 = 0.;
	double ran2 = 0.;
	int ind1 = 0;
	int ind2 = 0;
	double obj_val1 = 0.;
	double obj_val2 = 0.;
	int added = 0;
	//int num_stored = 0;
	double objval1 = 0;
	double objval2 = 0;
	int lpstat = 0;
	
	int numcols = CPXgetnumcols (env, temp_lp);
	int numrows = CPXgetnumrows (env, temp_lp);
	
	char *b = NULL;
	double *low_bd = NULL;
	double *up_bd = NULL;
	double *temp_x = NULL;
	int *indices = NULL;
	int *cstat = NULL;
	int *rstat = NULL;
			
	//int *integer_indices;
	temp_x = (double *) malloc (numcols*sizeof(double));
	//integer_indices = (int *) malloc (cur_numcols*sizeof(int));
	low_bd = (double *) malloc (numcols*sizeof(double));
	up_bd = (double *) malloc (numcols*sizeof(double));
	indices = (int *) malloc (numcols*sizeof(int));
	cstat = (int *) malloc (numcols*sizeof(int));
	rstat = (int *) malloc (numrows*sizeof(int));
			
	int num_added = 0;
	int i = 0;
	
	for(i=0;i<numcols;i++) indices[i] = i;
	
	status = CPXgetlb (env, temp_lp, low_bd, 0, numcols-1);
	if ( status ) {
		printf ("(%d) Failed to get lb's. Error code %d\n", __LINE__,status);
		goto TERMINATE;
	}
	status = CPXgetub (env, temp_lp, up_bd, 0, numcols-1);
	if ( status ) {
		printf ("(%d) Failed to get ub's. Error code %d\n", __LINE__,status);
		goto TERMINATE;
	}
			
	b = (char *) malloc (total_num_integer*sizeof(char));
			
	for(i=0;i<total_num_integer;i++)
	{
		b[i] = 'B';
	}
	int j = 0;
	for(j=0;j<heur_limit;j++)
	{
		ran1 = (double) rand() / ( (double) RAND_MAX);
		ind1 = round((num_stored-1)*ran1);
		ran2 = (double) rand() / ( (double) RAND_MAX);
		ind2 = round((num_stored-1)*ran2);
		if(ind2 == ind1) ind2 = (ind2 + 1) % num_x_to_store;
		for(i=0;i<total_num_integer;i++)
		{
			double ran3 = (double) rand() / ( (double) RAND_MAX);
			temp_x[i] = round( ran3*stored_x[ind1][integer_indices[i]] + (1.-ran3)*stored_x[ind2][integer_indices[i]] );
		}
			
		status = CPXchgbds (env, temp_lp, total_num_integer, integer_indices, b, temp_x);
			
		status = CPXlpopt(env, temp_lp);
		lpstat = CPXgetstat (env, temp_lp);
  		if(lpstat != 2 && lpstat !=3)
  		{
  			status = CPXgetx (env, temp_lp, temp_x, 0, numcols-1);
  			status = CPXgetbase (env, temp_lp, cstat, rstat);
      			if ( status ) {
       				printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
      			}
  			obj_val1 = obj1_extra_val;
  			obj_val2 = obj2_extra_val;
  			for(i=0;i<obj1_index;i++)
  			{
  				obj_val1 = obj_val1 + temp_x[i]*ob_coef1[i];
  				obj_val2 = obj_val2 + temp_x[i]*ob_coef2[i];
  			}
  			added = mock_insert(1,obj_val1,obj_val2,0,0,0,&tree);
  			if(added)
  			{
  				for(i=0;i<obj1_index;i++)
		      		{
			      		stored_x[x_rotation][i] = temp_x[i];
			      	}
			      	x_rotation = (x_rotation + 1) % num_x_to_store;
				added = 0;
		      		num_added++;
		      		num_stored = min(num_stored+1,num_x_to_store);
		      		PSA_full(env, lp_1, temp_x, cstat, rstat);
			      	int beg[1] = {0};
			      	int effortlevel[1] = {CPX_MIPSTART_AUTO};
			      	double *values = (double *) malloc (total_num_integer*sizeof(double));
			      	for(i=0;i<total_num_integer;i++) values[i] = temp_x[integer_indices[i]];
			      	status = CPXaddmipstarts (env_just_solve_mips, global_mip, 1, total_num_integer, beg, integer_indices,
                           					values, effortlevel, NULL);
                           	free(values);
  			}
  			else
	  		{
	  			CPXLPptr temp_lp2 = CPXcloneprob (env, temp_lp, &status);
	  			status = CPXchgobj (env, temp_lp2, numcols, indices, weighted_coefs);
				if ( status ) {
					printf ("(%d) Failed to get change obj coef. Error code %d\n", __LINE__,status);
					goto TERMINATE;
				}
				status = CPXlpopt(env, temp_lp2);
				lpstat = CPXgetstat (env, temp_lp2);
	  			if(lpstat != 2 && lpstat !=3)
	  			{
	  				status = CPXgetx (env, temp_lp2, temp_x, 0, numcols-1);
	  				status = CPXgetbase (env, temp_lp2, cstat, rstat);
	      				if ( status ) {
	       					printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
	      				}
	  				obj_val1 = obj1_extra_val;
	  				obj_val2 = obj2_extra_val;
	  				for(i=0;i<obj1_index;i++)
	  				{
	  					obj_val1 += temp_x[i]*ob_coef1[i];
	  					obj_val2 += temp_x[i]*ob_coef2[i];
	  				}
	  				added = mock_insert(1,obj_val1,obj_val2,0,0,0,&tree);
	  				if(added)
	  				{
	  					for(i=0;i<obj1_index;i++)
			      			{
				      			stored_x[x_rotation][i] = temp_x[i];
				      		}
				      		x_rotation = (x_rotation + 1) % num_x_to_store;
						added = 0;
			      			num_added++;
			      			num_stored = min(num_stored+1,num_x_to_store);
			      			PSA_full(env, lp_1, temp_x, cstat, rstat);
				      		int beg[1] = {0};
					      	int effortlevel[1] = {CPX_MIPSTART_AUTO};
					      	double *values = (double *) malloc (total_num_integer*sizeof(double));
					      	for(i=0;i<total_num_integer;i++) values[i] = temp_x[integer_indices[i]];
					      	status = CPXaddmipstarts (env_just_solve_mips, global_mip, 1, total_num_integer, beg, integer_indices,
				           					values, effortlevel, NULL);
				           	free(values);
	  				}
	  				else
			  		{
			  			CPXLPptr temp_lp3 = CPXcloneprob (env, temp_lp, &status);
			  			status = CPXchgobj (env, temp_lp2, numcols, indices, ob_coef1);
						if ( status ) {
							printf ("(%d) Failed to get change obj coef. Error code %d\n", __LINE__,status);
							goto TERMINATE;
						}
						status = CPXlpopt(env, temp_lp3);
						lpstat = CPXgetstat (env, temp_lp3);
			  			if(lpstat != 2 && lpstat !=3)
			  			{
			  				status = CPXgetx (env, temp_lp3, temp_x, 0, numcols-1);
			  				status = CPXgetbase (env, temp_lp3, cstat, rstat);
			      				if ( status ) {
			       					printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
			      				}
			  				obj_val1 = obj1_extra_val;
			  				obj_val2 = obj2_extra_val;
			  				for(i=0;i<obj1_index;i++)
			  				{
			  					obj_val1 += temp_x[i]*ob_coef1[i];
			  					obj_val2 += temp_x[i]*ob_coef2[i];
			  				}
			  				added = mock_insert(1,obj_val1,obj_val2,0,0,0,&tree);
			  				if(added)
			  				{
			  					for(i=0;i<obj1_index;i++)
					      			{
						      			stored_x[x_rotation][i] = temp_x[i];
						      		}
						      		x_rotation = (x_rotation + 1) % num_x_to_store;
								added = 0;
					      			num_added++;
					      			num_stored = min(num_stored+1,num_x_to_store);
					      			PSA_full(env, lp_1, temp_x, cstat, rstat);
						      		int beg[1] = {0};
							      	int effortlevel[1] = {CPX_MIPSTART_AUTO};
							      	double *values = (double *) malloc (total_num_integer*sizeof(double));
							      	for(i=0;i<total_num_integer;i++) values[i] = temp_x[integer_indices[i]];
							      	status = CPXaddmipstarts (env_just_solve_mips, global_mip, 1, total_num_integer, beg,
							      		integer_indices, values, effortlevel, NULL);
							   	free(values);
			  				}
			  			}
			  			CPXfreeprob(env,&temp_lp3);
			  		}
	  			}
	  			CPXfreeprob(env,&temp_lp2);
	  		}
  		}
  			
  		ran1 = (double) rand() / ( (double) RAND_MAX);
  		ind1 = round((num_stored-1)*ran1);
  		if(stored_x[ind1][1] > 1000.) ind1 = (ind1 + 1) % num_x_to_store;
  		ran2 = (double) rand() / ( (double) RAND_MAX);
		ind2 = round(total_num_integer*ran2);
		if(ind2 == total_num_integer) ind2--;
			
		for(i=0;i<total_num_integer;i++)
		{
			temp_x[i] = stored_x[ind1][integer_indices[i]];
		}
		if(temp_x[ind2] < up_bd[integer_indices[ind2]]) temp_x[ind2]++;
		else if(temp_x[ind2] > low_bd[integer_indices[ind2]]) temp_x[ind2]--;
		
		status = CPXchgbds (env, temp_lp, total_num_integer, integer_indices, b, temp_x);
			
		status = CPXlpopt(env, temp_lp);
		lpstat = CPXgetstat (env, temp_lp);
  		if(lpstat != 2 && lpstat !=3)
  		{
  			status = CPXgetx (env, temp_lp, temp_x, 0, numcols-1);
  			status = CPXgetbase (env, temp_lp, cstat, rstat);
      			if ( status ) {
       				printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
      			}
  			obj_val1 = obj1_extra_val;
  			obj_val2 = obj2_extra_val;
  			for(i=0;i<obj1_index;i++)
  			{
  				obj_val1 = obj_val1 + temp_x[i]*ob_coef1[i];
  				obj_val2 = obj_val2 + temp_x[i]*ob_coef2[i];
  			}
  			added = mock_insert(1,obj_val1,obj_val2,0,0,0,&tree);
  			if(added)
  			{
  				for(i=0;i<obj1_index;i++)
			      	{
			      		stored_x[x_rotation][i] = temp_x[i];
		 		}
		      		x_rotation = (x_rotation + 1) % num_x_to_store;
			      	added = 0;
			      	num_added++;
			      	num_stored = min(num_stored+1,num_x_to_store);
				PSA_full(env, lp_1, temp_x, cstat, rstat);
			      	int beg[1] = {0};
			      	int effortlevel[1] = {CPX_MIPSTART_AUTO};
			      	double *values = (double *) malloc (total_num_integer*sizeof(double));
			      	for(i=0;i<total_num_integer;i++) values[i] = temp_x[integer_indices[i]];
			      	status = CPXaddmipstarts (env_just_solve_mips, global_mip, 1, total_num_integer, beg, integer_indices,
                           					values, effortlevel, NULL);
                           	free(values);
  			}
  			else
	  		{
	  			CPXLPptr temp_lp2 = CPXcloneprob (env, temp_lp, &status);
	  			status = CPXchgobj (env, temp_lp2, numcols, indices, weighted_coefs);
				if ( status ) {
					printf ("(%d) Failed to get change obj coef. Error code %d\n", __LINE__,status);
					goto TERMINATE;
				}
				status = CPXlpopt(env, temp_lp2);
				lpstat = CPXgetstat (env, temp_lp2);
	  			if(lpstat != 2 && lpstat !=3)
	  			{
	  				status = CPXgetx (env, temp_lp2, temp_x, 0, numcols-1);
	  				status = CPXgetbase (env, temp_lp2, cstat, rstat);
	      				if ( status ) {
	       					printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
	      				}
	  				obj_val1 = obj1_extra_val;
	  				obj_val2 = obj2_extra_val;
	  				for(i=0;i<obj1_index;i++)
	  				{
	  					obj_val1 += temp_x[i]*ob_coef1[i];
	  					obj_val2 += temp_x[i]*ob_coef2[i];
	  				}
	  				added = mock_insert(1,obj_val1,obj_val2,0,0,0,&tree);
	  				if(added)
	  				{
	  					for(i=0;i<obj1_index;i++)
			      			{
				      			stored_x[x_rotation][i] = temp_x[i];
				      		}
				      		x_rotation = (x_rotation + 1) % num_x_to_store;
						added = 0;
			      			num_added++;
			      			num_stored = min(num_stored+1,num_x_to_store);
			      			PSA_full(env, lp_1, temp_x, cstat, rstat);
				      		int beg[1] = {0};
					      	int effortlevel[1] = {CPX_MIPSTART_AUTO};
					      	double *values = (double *) malloc (total_num_integer*sizeof(double));
					      	for(i=0;i<total_num_integer;i++) values[i] = temp_x[integer_indices[i]];
					      	status = CPXaddmipstarts (env_just_solve_mips, global_mip, 1, total_num_integer, beg, integer_indices,
				           					values, effortlevel, NULL);
				           	free(values);
	  				}
	  				else
			  		{
			  			CPXLPptr temp_lp3 = CPXcloneprob (env, temp_lp, &status);
			  			status = CPXchgobj (env, temp_lp2, numcols, indices, ob_coef1);
						if ( status ) {
							printf ("(%d) Failed to get change obj coef. Error code %d\n", __LINE__,status);
							goto TERMINATE;
						}
						status = CPXlpopt(env, temp_lp3);
						lpstat = CPXgetstat (env, temp_lp3);
			  			if(lpstat != 2 && lpstat !=3)
			  			{
			  				status = CPXgetx (env, temp_lp3, temp_x, 0, numcols-1);
			  				status = CPXgetbase (env, temp_lp3, cstat, rstat);
			      				if ( status ) {
			       					printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
			      				}
			  				obj_val1 = obj1_extra_val;
			  				obj_val2 = obj2_extra_val;
			  				for(i=0;i<obj1_index;i++)
			  				{
			  					obj_val1 += temp_x[i]*ob_coef1[i];
			  					obj_val2 += temp_x[i]*ob_coef2[i];
			  				}
			  				added = mock_insert(1,obj_val1,obj_val2,0,0,0,&tree);
			  				if(added)
			  				{
			  					for(i=0;i<obj1_index;i++)
					      			{
						      			stored_x[x_rotation][i] = temp_x[i];
						      		}
						      		x_rotation = (x_rotation + 1) % num_x_to_store;
								added = 0;
					      			num_added++;
					      			num_stored = min(num_stored+1,num_x_to_store);
					      			PSA_full(env, lp_1, temp_x, cstat, rstat);
						      		int beg[1] = {0};
							      	int effortlevel[1] = {CPX_MIPSTART_AUTO};
							      	double *values = (double *) malloc (total_num_integer*sizeof(double));
							      	for(i=0;i<total_num_integer;i++) values[i] = temp_x[integer_indices[i]];
							      	status = CPXaddmipstarts (env_just_solve_mips, global_mip, 1, total_num_integer, beg,
							      		 integer_indices, values, effortlevel, NULL);
							   	free(values);
			  				}
			  			}
			  			CPXfreeprob(env,&temp_lp3);
			  		}
	  			}
	  			CPXfreeprob(env,&temp_lp2);
	  		}
  		}
  			
  		ran1 = (double) rand() / ( (double) RAND_MAX);
  		ind1 = round((num_stored-1)*ran1);
  		if(stored_x[ind1][1] > 1000.) ind1 = (ind1 + 1) % num_x_to_store;
			
		for(i=0;i<total_num_integer;i++)
		{
			temp_x[i] = stored_x[ind1][integer_indices[i]];
			ran2 = (double) rand() / ( (double) RAND_MAX);
			if(ran2 < .33) temp_x[i] = fmax(low_bd[integer_indices[i]],temp_x[i]-1);
			else if(ran2 < .66) temp_x[i] = fmin(up_bd[integer_indices[i]],temp_x[i]+1);
		}
		status = CPXchgbds (env, temp_lp, total_num_integer, integer_indices, b, temp_x);
			
		status = CPXlpopt(env, temp_lp);
		lpstat = CPXgetstat (env, temp_lp);
  		if(lpstat != 2 && lpstat !=3)
  		{
  			status = CPXgetx (env, temp_lp, temp_x, 0, numcols-1);
  			status = CPXgetbase (env, temp_lp, cstat, rstat);
      			if ( status ) {
       				printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
      			}
  			obj_val1 = obj1_extra_val;
  			obj_val2 = obj2_extra_val;
  			for(i=0;i<obj1_index;i++)
  			{
  				obj_val1 = obj_val1 + temp_x[i]*ob_coef1[i];
  				obj_val2 = obj_val2 + temp_x[i]*ob_coef2[i];
  			}
  			added = mock_insert(1,obj_val1,obj_val2,0,0,0,&tree);
  			if(added)
  			{
  				for(i=0;i<obj1_index;i++)
			      	{
					stored_x[x_rotation][i] = temp_x[i];
			      	}
			      	x_rotation = (x_rotation + 1) % num_x_to_store;
				added = 0;
			      	num_added++;
			      	num_stored = min(num_stored+1,num_x_to_store);
			      	PSA_full(env, lp_1, temp_x, cstat, rstat);
				int beg[1] = {0};
			      	int effortlevel[1] = {CPX_MIPSTART_AUTO};
			      	double *values = (double *) malloc (total_num_integer*sizeof(double));
			      	for(i=0;i<total_num_integer;i++) values[i] = temp_x[integer_indices[i]];
			      	status = CPXaddmipstarts (env_just_solve_mips, global_mip, 1, total_num_integer, beg, integer_indices,
                           					values, effortlevel, NULL);
                           	free(values);
  			}
  			else
	  		{
	  			CPXLPptr temp_lp2 = CPXcloneprob (env, temp_lp, &status);
	  			status = CPXchgobj (env, temp_lp2, numcols, indices, weighted_coefs);
				if ( status ) {
					printf ("(%d) Failed to get change obj coef. Error code %d\n", __LINE__,status);
					goto TERMINATE;
				}
				status = CPXlpopt(env, temp_lp2);
				lpstat = CPXgetstat (env, temp_lp2);
	  			if(lpstat != 2 && lpstat !=3)
	  			{
	  				status = CPXgetx (env, temp_lp2, temp_x, 0, numcols-1);
	  				status = CPXgetbase (env, temp_lp2, cstat, rstat);
	      				if ( status ) {
	       					printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
	      				}
	  				obj_val1 = obj1_extra_val;
	  				obj_val2 = obj2_extra_val;
	  				for(i=0;i<obj1_index;i++)
	  				{
	  					obj_val1 += temp_x[i]*ob_coef1[i];
	  					obj_val2 += temp_x[i]*ob_coef2[i];
	  				}
	  				added = mock_insert(1,obj_val1,obj_val2,0,0,0,&tree);
	  				if(added)
	  				{
	  					for(i=0;i<obj1_index;i++)
			      			{
				      			stored_x[x_rotation][i] = temp_x[i];
				      		}
				      		x_rotation = (x_rotation + 1) % num_x_to_store;
						added = 0;
			      			num_added++;
			      			num_stored = min(num_stored+1,num_x_to_store);
			      			PSA_full(env, lp_1, temp_x, cstat, rstat);
				      		int beg[1] = {0};
					      	int effortlevel[1] = {CPX_MIPSTART_AUTO};
					      	double *values = (double *) malloc (total_num_integer*sizeof(double));
					      	for(i=0;i<total_num_integer;i++) values[i] = temp_x[integer_indices[i]];
					      	status = CPXaddmipstarts (env_just_solve_mips, global_mip, 1, total_num_integer, beg, integer_indices,
				           					values, effortlevel, NULL);
				           	free(values);
	  				}
	  				else
			  		{
			  			CPXLPptr temp_lp3 = CPXcloneprob (env, temp_lp, &status);
			  			status = CPXchgobj (env, temp_lp2, numcols, indices, ob_coef1);
						if ( status ) {
							printf ("(%d) Failed to get change obj coef. Error code %d\n", __LINE__,status);
							goto TERMINATE;
						}
						status = CPXlpopt(env, temp_lp3);
						lpstat = CPXgetstat (env, temp_lp3);
			  			if(lpstat != 2 && lpstat !=3)
			  			{
			  				status = CPXgetx (env, temp_lp3, temp_x, 0, numcols-1);
			  				status = CPXgetbase (env, temp_lp3, cstat, rstat);
			      				if ( status ) {
			       					printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
			      				}
			  				obj_val1 = obj1_extra_val;
			  				obj_val2 = obj2_extra_val;
			  				for(i=0;i<obj1_index;i++)
			  				{
			  					obj_val1 += temp_x[i]*ob_coef1[i];
			  					obj_val2 += temp_x[i]*ob_coef2[i];
			  				}
			  				added = mock_insert(1,obj_val1,obj_val2,0,0,0,&tree);
			  				if(added)
			  				{
			  					for(i=0;i<obj1_index;i++)
					      			{
						      			stored_x[x_rotation][i] = temp_x[i];
						      		}
						      		x_rotation = (x_rotation + 1) % num_x_to_store;
								added = 0;
					      			num_added++;
					      			num_stored = min(num_stored+1,num_x_to_store);
					      			PSA_full(env, lp_1, temp_x, cstat, rstat);
						      		int beg[1] = {0};
							      	int effortlevel[1] = {CPX_MIPSTART_AUTO};
							      	double *values = (double *) malloc (total_num_integer*sizeof(double));
							      	for(i=0;i<total_num_integer;i++) values[i] = temp_x[integer_indices[i]];
							      	status = CPXaddmipstarts (env_just_solve_mips, global_mip, 1, total_num_integer, beg,
							      		 integer_indices, values, effortlevel, NULL);
							   	free(values);
			  				}
			  			}
			  			CPXfreeprob(env,&temp_lp3);
			  		}
	  			}
	  			CPXfreeprob(env,&temp_lp2);
	  		}
  		}
  	}
  	TERMINATE:
  	if(indices) free_and_null ((char **) &indices);
    	if(cstat) free_and_null ((char **) &cstat);
    	if(rstat) free_and_null ((char **) &rstat);
    	if(b)  free_and_null ((char **) &b);
    	if(low_bd)  free_and_null ((char **) &low_bd);
   	if(up_bd)  free_and_null ((char **) &up_bd);
    	if(temp_x)  free_and_null ((char **) &temp_x);
    	CPXfreeprob(env,&temp_lp);
    	
    	return num_added;
}

/*********************************************************************************************************************** 

	This function is the genetic algorithm used for generating integer feasible solutions during preprocessing
	when only binary variables are involved in the problem.
	
***********************************************************************************************************************/

int binary_heuristic(CPXCENVptr env, CPXLPptr lp1, CPXLPptr lp2)
{
/*	printf("in bin heur\n");*/
	int status = 0, i = 0;
	CPXLPptr temp_lp = CPXcloneprob (env, lp_2, &status);
	double ran1 = 0.;
	double ran2 = 0.;
	int ind1 = 0;
	int ind2 = 0;
	double obj_val1 = obj1_extra_val;
	double obj_val2 = obj2_extra_val;
	int added = 0;
	double objval1 = 0;
	double objval2 = 0;
	int lpstat = 0;
	
	int numcols = CPXgetnumcols (env, temp_lp);
	int numrows = CPXgetnumrows (env, temp_lp);
	
	int *indices = NULL;
	
	indices = (int *) malloc (numcols*sizeof(int));
	for(i=0;i<numcols;i++) indices[i] = i;
	
	char *b = NULL;
	double *low_bd = NULL;
	double *up_bd = NULL;
	double *temp_x = NULL;
	int *cstat = NULL;
	int *rstat = NULL;
			
	temp_x = (double *) malloc (numcols*sizeof(double));
	low_bd = (double *) malloc (numcols*sizeof(double));
	up_bd = (double *) malloc (numcols*sizeof(double));
	cstat = (int *) malloc (numcols*sizeof(int));
	rstat = (int *) malloc (5*numrows*sizeof(int));
			
	int num_added = 0;
	
	status = CPXgetlb (env, temp_lp, low_bd, 0, numcols-1);
	if ( status ) {
		printf ("(%d) Failed to get lb's. Error code %d\n", __LINE__,status);
		goto TERMINATE;
	}
	status = CPXgetub (env, temp_lp, up_bd, 0, numcols-1);
	if ( status ) {
		printf ("(%d) Failed to get ub's. Error code %d\n", __LINE__,status);
		goto TERMINATE;
	}
			
	b = (char *) malloc (num_integer*sizeof(char));
			
	for(i=0;i<num_integer;i++)
	{
		b[i] = 'B';
	}
	int j = 0;
	for(j=0;j<heur_limit;j++)
	{
		ran1 = (double) rand() / ( (double) RAND_MAX);
		ind1 = round((num_stored-1)*ran1);
		ran2 = (double) rand() / ( (double) RAND_MAX);
		ind2 = round((num_stored-1)*ran2);
		if(ind2 == ind1) ind2 = (ind2 + 1) % num_x_to_store;
		for(i=0;i<num_integer;i++)
		{
			double ran3 = (double) rand() / ( (double) RAND_MAX);
			if(ran3 < .5) temp_x[i] = stored_x[ind1][integer_indices[i]];
			else temp_x[i] = stored_x[ind2][integer_indices[i]];
		}
		status = CPXchgbds (env, temp_lp, num_integer, integer_indices, b, temp_x);
			
		status = CPXlpopt(env, temp_lp);
		lpstat = CPXgetstat (env, temp_lp);
  		if(lpstat != 2 && lpstat !=3)
  		{
  			status = CPXgetx (env, temp_lp, temp_x, 0, numcols-1);
  			status = CPXgetbase (env, temp_lp, cstat, rstat);
      			if ( status ) {
       				printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
      			}
  			obj_val1 = obj1_extra_val;
  			obj_val2 = obj2_extra_val;
  			for(i=0;i<obj1_index;i++)
  			{
  				obj_val1 = obj_val1 + temp_x[i]*ob_coef1[i];
  				obj_val2 = obj_val2 + temp_x[i]*ob_coef2[i];
  			}
  			added = mock_insert(1,obj_val1,obj_val2,0,0,0,&tree);
  			if(added)
  			{
  				for(i=0;i<obj1_index;i++)
		      		{
			      		stored_x[x_rotation][i] = temp_x[i];
			      	}
			      	x_rotation = (x_rotation + 1) % num_x_to_store;
				added = 0;
		      		num_added++;
		      		num_stored = min(num_stored+1,num_x_to_store);
		      		PSA_full(env, lp_1, temp_x, cstat, rstat);
				int beg[1] = {0};
			      	int effortlevel[1] = {CPX_MIPSTART_AUTO};
			      	double *values = (double *) malloc (total_num_integer*sizeof(double));
			      	for(i=0;i<total_num_integer;i++) values[i] = temp_x[integer_indices[i]];
			      	status = CPXaddmipstarts (env_just_solve_mips, global_mip, 1, total_num_integer, beg, integer_indices,
                           					values, effortlevel, NULL);
                           	free(values);
  			}
  			else
	  		{
	  			CPXLPptr temp_lp2 = CPXcloneprob (env, temp_lp, &status);
	  			status = CPXchgobj (env, temp_lp2, numcols, indices, weighted_coefs);
				if ( status ) {
					printf ("(%d) Failed to get change obj coef. Error code %d\n", __LINE__,status);
					goto TERMINATE;
				}
				status = CPXlpopt(env, temp_lp2);
				lpstat = CPXgetstat (env, temp_lp2);
	  			if(lpstat != 2 && lpstat !=3)
	  			{
	  				status = CPXgetx (env, temp_lp2, temp_x, 0, numcols-1);
	  				status = CPXgetbase (env, temp_lp2, cstat, rstat);
	      				if ( status ) {
	       					printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
	      				}
	  				obj_val1 = obj1_extra_val;
		  			obj_val2 = obj2_extra_val;
		  			for(i=0;i<obj1_index;i++)
		  			{
		  				obj_val1 = obj_val1 + temp_x[i]*ob_coef1[i];
		  				obj_val2 = obj_val2 + temp_x[i]*ob_coef2[i];
		  			}
	  				added = mock_insert(1,obj_val1,obj_val2,0,0,0,&tree);
	  				if(added)
	  				{
	  					for(i=0;i<obj1_index;i++)
			      			{
				      			stored_x[x_rotation][i] = temp_x[i];
				      		}
				      		x_rotation = (x_rotation + 1) % num_x_to_store;
						added = 0;
			      			num_added++;
			      			num_stored = min(num_stored+1,num_x_to_store);
			      			PSA_full(env, lp_1, temp_x, cstat, rstat);
						int beg[1] = {0};
					      	int effortlevel[1] = {CPX_MIPSTART_AUTO};
					      	double *values = (double *) malloc (total_num_integer*sizeof(double));
					      	for(i=0;i<total_num_integer;i++) values[i] = temp_x[integer_indices[i]];
					      	status = CPXaddmipstarts (env_just_solve_mips, global_mip, 1, total_num_integer, beg, integer_indices,
				           					values, effortlevel, NULL);
				           	free(values);
	  				}
	  			}
	  			CPXfreeprob(env,&temp_lp2);
	  		}
  		}
  			
  		ran1 = (double) rand() / ( (double) RAND_MAX);
  		ind1 = round((num_stored-1)*ran1);
  		if(stored_x[ind1][1] > 1000.) ind1 = (ind1 + 1) % num_x_to_store;
  		ran2 = (double) rand() / ( (double) RAND_MAX);
		ind2 = round(num_integer*ran2);
			
		for(i=0;i<num_integer;i++)
		{
			temp_x[i] = stored_x[ind1][integer_indices[i]];
		}
		temp_x[ind2] = ((int)temp_x[ind2] + 1) % 2;
		status = CPXchgbds (env, temp_lp, num_integer, integer_indices, b, temp_x);
			
		status = CPXlpopt(env, temp_lp);
		lpstat = CPXgetstat (env, temp_lp);
  		if(lpstat != 2 && lpstat !=3)
  		{
  			status = CPXgetx (env, temp_lp, temp_x, 0, numcols-1);
  			status = CPXgetbase (env, temp_lp, cstat, rstat);
      			if ( status ) {
       				printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
      			}
  			obj_val1 = obj1_extra_val;
  			obj_val2 = obj2_extra_val;
  			for(i=0;i<obj1_index;i++)
  			{
  				obj_val1 = obj_val1 + temp_x[i]*ob_coef1[i];
  				obj_val2 = obj_val2 + temp_x[i]*ob_coef2[i];
  			}
  			added = mock_insert(1,obj_val1,obj_val2,0,0,0,&tree);
  			if(added)
  			{
  				for(i=0;i<obj1_index;i++)
			      	{
			      		stored_x[x_rotation][i] = temp_x[i];
		 		}
		      		x_rotation = (x_rotation + 1) % num_x_to_store;
			      	added = 0;
			      	num_added++;
			      	num_stored = min(num_stored+1,num_x_to_store);
				PSA_full(env, lp_1, temp_x, cstat, rstat);
				int beg[1] = {0};
			      	int effortlevel[1] = {CPX_MIPSTART_AUTO};
			      	double *values = (double *) malloc (total_num_integer*sizeof(double));
			      	for(i=0;i<total_num_integer;i++) values[i] = temp_x[integer_indices[i]];
			      	status = CPXaddmipstarts (env_just_solve_mips, global_mip, 1, total_num_integer, beg, integer_indices,
                           					values, effortlevel, NULL);
                           	free(values);
  			}
  			else
	  		{
	  			CPXLPptr temp_lp2 = CPXcloneprob (env, temp_lp, &status);
	  			status = CPXchgobj (env, temp_lp2, numcols, indices, weighted_coefs);
				if ( status ) {
					printf ("Failed to get change obj coef. Error code %d\n", status);
					goto TERMINATE;
				}
				status = CPXlpopt(env, temp_lp2);
				lpstat = CPXgetstat (env, temp_lp2);
	  			if(lpstat != 2 && lpstat !=3)
	  			{
	  				status = CPXgetx (env, temp_lp2, temp_x, 0, numcols-1);
	  				status = CPXgetbase (env, temp_lp2, cstat, rstat);
	      				if ( status ) {
	       					printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
	      				}
	  				obj_val1 = obj1_extra_val;
		  			obj_val2 = obj2_extra_val;
		  			for(i=0;i<obj1_index;i++)
		  			{
		  				obj_val1 = obj_val1 + temp_x[i]*ob_coef1[i];
		  				obj_val2 = obj_val2 + temp_x[i]*ob_coef2[i];
		  			}
	  				added = mock_insert(1,obj_val1,obj_val2,0,0,0,&tree);
	  				if(added)
	  				{
	  					for(i=0;i<obj1_index;i++)
			      			{
				      			stored_x[x_rotation][i] = temp_x[i];
				      		}
				      		x_rotation = (x_rotation + 1) % num_x_to_store;
						added = 0;
			      			num_added++;
			      			num_stored = min(num_stored+1,num_x_to_store);
			      			PSA_full(env, lp_1, temp_x, cstat, rstat);
						int beg[1] = {0};
					      	int effortlevel[1] = {CPX_MIPSTART_AUTO};
					      	double *values = (double *) malloc (total_num_integer*sizeof(double));
					      	for(i=0;i<total_num_integer;i++) values[i] = temp_x[integer_indices[i]];
					      	status = CPXaddmipstarts (env_just_solve_mips, global_mip, 1, total_num_integer, beg, integer_indices,
				           					values, effortlevel, NULL);
				           	free(values);
	  				}
	  			}
	  			CPXfreeprob(env,&temp_lp2);
	  		}
  		}
  			
  		ran1 = (double) rand() / ( (double) RAND_MAX);
  		ind1 = round((num_stored-1)*ran1);
  		if(stored_x[ind1][1] > 1000.) ind1 = (ind1 + 1) % num_x_to_store;
			
		for(i=0;i<num_integer;i++)
		{
			temp_x[i] = stored_x[ind1][integer_indices[i]];
			ran2 = (double) rand() / ( (double) RAND_MAX);
			if(ran2 < .1) 
			{
				temp_x[i] = ( (int) temp_x[i] + 1) % 2;
			}
		}
			
		status = CPXchgbds (env, temp_lp, num_integer, integer_indices, b, temp_x);
			
		status = CPXlpopt(env, temp_lp);
		lpstat = CPXgetstat (env, temp_lp);
  		if(lpstat != 2 && lpstat !=3)
  		{
  			status = CPXgetx (env, temp_lp, temp_x, 0, numcols-1);
  			status = CPXgetbase (env, temp_lp, cstat, rstat);
      			if ( status ) {
       				printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
      			}
  			obj_val1 = obj1_extra_val;
  			obj_val2 = obj2_extra_val;
  			for(i=0;i<obj1_index;i++)
  			{
  				obj_val1 = obj_val1 + temp_x[i]*ob_coef1[i];
  				obj_val2 = obj_val2 + temp_x[i]*ob_coef2[i];
  			}
  			added = mock_insert(1,obj_val1,obj_val2,0,0,0,&tree);
  			if(added)
  			{
  				for(i=0;i<obj1_index;i++)
			      	{
					stored_x[x_rotation][i] = temp_x[i];
			      	}
			      	x_rotation = (x_rotation + 1) % num_x_to_store;
				added = 0;
			      	num_added++;
			      	num_stored = min(num_stored+1,num_x_to_store);
			      	PSA_full(env, lp_1, temp_x, cstat, rstat);
				int beg[1] = {0};
			      	int effortlevel[1] = {CPX_MIPSTART_AUTO};
			      	double *values = (double *) malloc (total_num_integer*sizeof(double));
			      	for(i=0;i<total_num_integer;i++) values[i] = temp_x[integer_indices[i]];
			      	status = CPXaddmipstarts (env_just_solve_mips, global_mip, 1, total_num_integer, beg, integer_indices,
                           					values, effortlevel, NULL);
                           	free(values);
  			}
  			else
	  		{
	  			CPXLPptr temp_lp2 = CPXcloneprob (env, temp_lp, &status);
	  			status = CPXchgobj (env, temp_lp2, numcols, indices, weighted_coefs);
				if ( status ) {
					printf ("(%d) Failed to get change obj coef. Error code %d\n", __LINE__,status);
					goto TERMINATE;
				}
				status = CPXlpopt(env, temp_lp2);
				lpstat = CPXgetstat (env, temp_lp2);
	  			if(lpstat != 2 && lpstat !=3)
	  			{
	  				status = CPXgetx (env, temp_lp2, temp_x, 0, numcols-1);
	  				status = CPXgetbase (env, temp_lp2, cstat, rstat);
	      				if ( status ) {
	       					printf ("CPXcopybase, Failed to copy basis, error code %d\n", status);
	      				}
	  				obj_val1 = obj1_extra_val;
		  			obj_val2 = obj2_extra_val;
		  			for(i=0;i<obj1_index;i++)
		  			{
		  				obj_val1 = obj_val1 + temp_x[i]*ob_coef1[i];
		  				obj_val2 = obj_val2 + temp_x[i]*ob_coef2[i];
		  			}
	  				added = mock_insert(1,obj_val1,obj_val2,0,0,0,&tree);
	  				if(added)
	  				{
	  					for(i=0;i<obj1_index;i++)
			      			{
				      			stored_x[x_rotation][i] = temp_x[i];
				      		}
				      		x_rotation = (x_rotation + 1) % num_x_to_store;
						added = 0;
			      			num_added++;
			      			num_stored = min(num_stored+1,num_x_to_store);
			      			PSA_full(env, lp_1, temp_x, cstat, rstat);
						int beg[1] = {0};
					      	int effortlevel[1] = {CPX_MIPSTART_AUTO};
					      	double *values = (double *) malloc (total_num_integer*sizeof(double));
					      	for(i=0;i<total_num_integer;i++) values[i] = temp_x[integer_indices[i]];
					      	status = CPXaddmipstarts (env_just_solve_mips, global_mip, 1, total_num_integer, beg, integer_indices,
				           					values, effortlevel, NULL);
				           	free(values);
	  				}
	  			}
	  			CPXfreeprob(env,&temp_lp2);
	  		}
  		}
  	}
  	TERMINATE:
  	
    	if(cstat) free_and_null ((char **) &cstat);
    	if(rstat) free_and_null ((char **) &rstat);
    	if(b)  free_and_null ((char **) &b);
    	if(low_bd)  free_and_null ((char **) &low_bd);
   	if(up_bd)  free_and_null ((char **) &up_bd);
    	if(temp_x)  free_and_null ((char **) &temp_x);
    	if(indices)  free_and_null ((char **) &indices);
    	CPXfreeprob(env,&temp_lp);
    	
    	return num_added;
}

user_data *userhandle_current = NULL;
user_data *userhandle_up = NULL;
user_data *userhandle_down = NULL;
double *x_ws = NULL;
double *x1 = NULL;
double *x2 = NULL;

double mip_slope = 0.;
int direction = 0;
int node_op_mips_solved = 0;
int num_iterations_since_last_add = 0;

int been_done = 0;
CPXCLPptr orig_mip = NULL;  
int seqnum2 = -1;
int seqnum1 = 0;
int branch_iterations = 0;
int printing_in_setbranch = 0;
int infeasible_in_cutcallback = 0;

clock_t start_time, finish_time, mip_start_time, mip_finish_time, bd_red_start_time, bd_red_finish_time;
double cumulative_time = 0.;
int printed_yet = 0, fathoming = 0, pareto = 0, bound_reduction = 0, last_cutcallback_seqnum = -1;
int ws_mip_opt = 0, ob1_mip_opt = 0, ob2_mip_opt = 0, mip_solved = 0, left_side_dom = 0, right_side_dom = 0, last_time_shown = 0, first_time_showing_progress = 1;
int break_early = 0, cut_problem_created = 0;
  
     /**************************************************************************************/
     /*	 _  _  ____  ____  ____    ____  ____  ____    ____  ____   __   __ _   ___  _  _  */
     /*	/ )( \/ ___)(  __)(  _ \  / ___)(  __)(_  _)  (  _ \(  _ \ / _\ (  ( \ / __)/ )( \ */
     /*	) \/ (\___ \ ) _)  )   /  \___ \ ) _)   )(     ) _ ( )   //    \/    /( (__ ) __ ( */
     /*	\____/(____/(____)(__\_)  (____/(____) (__)   (____/(__\_)\_/\_/\_)__) \___)\_)(_/ */
     /* 										   */
     /*  This is the callback used for branching. We also process each node in this        */
     /*  callback because of the ease of using the information discovered here in          */
     /*  branching.                                                                        */
     /**************************************************************************************/

int CPXPUBLIC
usersetbranch (CPXCENVptr   env,
               void         *cbdata,
               int          wherefrom,
               void         *cbhandle,
               int          brtype,
               int          sos,
               int          nodecnt,
               int          bdcnt,
               const int    *nodebeg,
               const int    *indices,
               const char   *lu,
               const double *bd,
               const double *nodeest,
               int          *useraction_p)
{

	/*************** Getting started *******************************************************/

	if(printing_in_setbranch) printf("********* in the branch callback *******\n");
/*	printf("indices: %d, %d\t nodeest: %lf, %lf\n",indices[0],indices[1],nodeest[0],nodeest[1]);*/
	start_time = clock();
	
	cumulative_time = (double)(start_time - start_BB) / CLOCKS_PER_SEC;
	
	int status = 0, i = 0, k = 0;

	branch_iterations++;
	print_on = 1;
   	
   	frac_index = -1;
   	frac_val = 0.;
 
   	char     varlu[1];
   	double   varbd[1];
   	int 	 vars[1];
   	char     varlu2[2];
   	double   varbd2[2];
   	int	 vars2[2];
   	char     varlu3[5];
   	double   varbd3[5];
   	int 	 vars3[5];
   	seqnum1 += 2;
   	seqnum2 += 2;
   	
   	CPXLPptr nodelp = NULL, nodelp2 = NULL, nodelp_copy = NULL, nodelp_copy2 = NULL, lp_ob1 = NULL, lp_ob2 = NULL, mip_ob1 = NULL, mip_ob2 = NULL;
	int lpstat = 0, depth = 0;
	int      j = -1, bestj = -1;
   	double   maxinf = -CPX_INFBOUND;
   	double   xj_inf = 0., objval = 0., best_bound = 0., projection = 0., diff = 0.;
   	clock_t start_mipsolve, finish_mipsolve;
   	
   	int seqnum = -1;
   	int add_check = 0, all_feas = 0, surplus = 0, infeas = 0, ws_lp_dom = 0, ob1_lp_dom = 0, ob2_lp_dom = 0; 
   	int ws_mip_ran = 0, ob1_mip_ran = 0, ob2_mip_ran = 0, ob1_lp_been_solved = 0, ob2_lp_been_solved = 0, ws_lp_int_feas = 0;
   	int ob1_lp_int_feas = 0, ob2_lp_int_feas = 0, right_pt_dom = 0, right_seg_dom = 0;
   	int ws_still_feas = 1, ob1_still_feas = 1, ob2_still_feas = 1, going_back = 0, ws_sol_interior = 0, mipstarts_copied = 0;
   	int from_pareto = 0, bound_reduced_from_PSA_reduce_left = 0, going_back_for_lp1 = 0, going_back_for_lp2 = 0, ws_lp_sol_on_top_bd = 0;
   	int ws_lp_sol_on_right_bd = 0;
   	
   	double endpoint1_x = 0., endpoint2_x = 0., endpoint1_y = 0., endpoint2_y = 0.;
   	
   	int PSA_right_check = 0, PSA_left_check = 0;
   	
	int *feas = NULL;
	char *low_up = NULL;
	double *lb_ = NULL;
	double *ub_ = NULL;
	double *bds_br1 = NULL;
	double *bds_br2 = NULL;
	double *var_ubs = NULL;
	double *var_lbs = NULL;
	int *br_ind = NULL;
	int *indexes = NULL;
	int *cstat = NULL;
	int *rstat = NULL;
	int *cstat_ws = NULL;
	int *rstat_ws = NULL;
	
	if(!x_ws) x_ws = (double *) malloc (cur_numcols*sizeof(double));
	if(!x1) x1 = (double *) malloc (cur_numcols*sizeof(double));
	if(!x2) x2 = (double *) malloc (cur_numcols*sizeof(double));
	
	/*************** Exit early if too much time has been used *************************/
	
	if(cumulative_time > max_time || branch_iterations > max_nodes || (show_progress && break_early))
	{
		show_progress = 0;
		if(cumulative_time > max_time + max_time_build_dual_bound && !approximate_dual_bd)
		{
			approximate_dual_bd = 1;
			printf("Maximimum time allowed for generating dual bound has been exceeded. Approximating bound from here on.\n");
		}
		if(cumulative_time > max_time + 2*max_time_build_dual_bound)
		{
			if(!quit_generating_dual_bd)
			{
				quit_generating_dual_bd = 1;
				printf("Double the maximimum time allowed for generating dual bound has been exceeded. Exitting!!\nWARNING: This indicates that the dual bound reported will most likely not be accurate and thus the quality of the primal solutions is UNKNOWN!!!");
			}
			fathoming = 1;
			*useraction_p = CPX_CALLBACK_SET;
			nodecnt = 0;
	  		goto TERMINATE;	
		}
		if(!printed_yet)
		{
			printed_yet = 1;
			int nodesleft = 0;
			
			duration_BB = cumulative_time;
			
			if ((status = CPXgetcallbackinfo (env, cbdata, wherefrom, 
			CPX_CALLBACK_INFO_NODES_LEFT, &nodesleft))) goto TERMINATE;
			
			if(cumulative_time > max_time) printf("Branch and Bound time has exceeded limit of %lf! Terminating after %lf seconds.\n",
				max_time,cumulative_time); 
			else if(branch_iterations > max_nodes) printf("Branch and Bound has processed the maximum number of nodes, %d! Terminating after %lf seconds.\n",max_nodes,cumulative_time);
			else if(show_progress && break_early) printf("Duality gap is below desired threshold of %lf! Terminating after %lf seconds.\n",
				duality_gap_limit,cumulative_time);
			printf("There are currently %d open nodes left, closing them and calculating global dual bound.\n",nodesleft); 
			fprintf(bb_results,"%lf\t",cumulative_time);
			fprintf(bb_results,"%d\t",nodesleft);
			fclose(bb_results);
			bb_results = fopen ("bb_results.txt", "a+");
			
			time_per = max_time_build_dual_bound/nodesleft;
		} 
		CPXLPptr nodelp = NULL;
		status = CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp);
		if ( status ) {
			printf ("CPXgetcallbacknodelp, Failed to catch nodelp, error code %d\n", status);
			goto TERMINATE;
		}
		printf("(%d) about to call build_dual_bd\n",__LINE__);
		build_dual_bd(env,nodelp);
		
		fathoming = 1;
		*useraction_p = CPX_CALLBACK_SET;
		nodecnt = 0;
  		goto TERMINATE;	
	}
	
	status = CPXgetcallbacknodeinfo(env,
                                 	cbdata,
                                 	wherefrom,
                                 	0,
                                 	CPX_CALLBACK_INFO_NODE_SEQNUM_LONG,
                                 	&seqnum);
   	
	if(last_cutcallback_seqnum == seqnum)
	{
		if(printing_in_setbranch) printf("node %d was processed during cutcallback, just branch\n",seqnum);
		
		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
							reduced_subprob_y_lb,reduced_subprob_y_lb);					
		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
						reduced_subprob_y_ub,reduced_subprob_y_ub);
		
		status = CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp);
		if ( status ) {
			printf ("CPXgetcallbacknodelp, Failed to catch nodelp, error code %d\n", status);
			goto TERMINATE;
		}
	   	
	   	cur_numcols = CPXgetnumcols (env, nodelp);
		cur_numrows = CPXgetnumrows (env, nodelp); 
	
		if(fathoming)
		{
			*useraction_p = CPX_CALLBACK_SET;
			nodecnt = 0;
	  		goto TERMINATE;
		}
	
		if(pareto)
		{
			pareto = 0;
			goto PARETO_BRANCH2;
		}
		goto BRANCHING2;
	}
	
	BEGINNING:
	
	if(printing_in_setbranch) printf("Processing node %d from branch callback\n",seqnum);
	
	ws_mip_opt = 0;
	ob1_mip_opt = 0;
   	ob2_mip_opt = 0;
   	mip_solved = 0;
   	bound_reduction = 0;
   	left_side_dom = 0;
   	right_side_dom = 0;
   	fathoming = 0;
   	cut_problem_created = 0;
	
	double ws_lp_objvals[2] = {0.,0.};
	double ob1_lp_objvals[2] = {0.,0.};
	double ob2_lp_objvals[2] = {0.,0.};
	double ws_mip_objvals[2] = {0.,0.};
	double ob1_mip_objvals[2] = {0.,0.};
	double ob2_mip_objvals[2] = {0.,0.};
	double lbs[2] = {0.,0.};
	double ubs[2] = {0.,0.};
	
	status = CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp);
	if ( status ) {
		printf ("CPXgetcallbacknodelp, Failed to catch nodelp, error code %d\n", status);
		goto TERMINATE;
	}
   	
   	cur_numcols = CPXgetnumcols (env, nodelp);
	cur_numrows = CPXgetnumrows (env, nodelp); 
	
	if(!x_ws) x_ws = (double *) malloc (cur_numcols*sizeof(double));
   	if(!x1) x1 = (double *) malloc (cur_numcols*sizeof(double));
   	if(!x2) x2 = (double *) malloc (cur_numcols*sizeof(double));
	feas = (int *) malloc (cur_numcols*sizeof(int));
	low_up = (char *) malloc ((2*total_num_integer+16)*sizeof(char));
	lb_ = (double *) malloc ((cur_numcols)*sizeof(double));
	ub_ = (double *) malloc ((cur_numcols)*sizeof(double));
	bds_br1 = (double *) malloc ((2*total_num_integer+16)*sizeof(double));
	bds_br2 = (double *) malloc ((2*total_num_integer+16)*sizeof(double));
	br_ind = (int *) malloc ((2*total_num_integer+16)*sizeof(int));
	indexes = (int *) malloc (cur_numcols*sizeof(int));
	cstat = (int *) malloc (cur_numcols*sizeof(int));
	rstat = (int *) malloc (cur_numrows*sizeof(int));
	cstat_ws = (int *) malloc (cur_numcols*sizeof(int));
	rstat_ws = (int *) malloc (cur_numrows*sizeof(int));
	for(i=0;i<cur_numcols;i++) indexes[i] = i;
	
/*	printf("(%d) mallocing uh's up & down\n",__LINE__);*/
	userhandle_up = (user_data*) malloc( sizeof( user_data ) );
/*	if(!userhandle_up->x_ws) */
	userhandle_up->x_ws = calloc ((cur_numcols+1),sizeof(double));
	userhandle_up->x1 = calloc (cur_numcols,sizeof(double));
	userhandle_up->x2 = calloc (cur_numcols,sizeof(double));
	userhandle_up->ws_still_feas = 0;
	userhandle_up->ob1_still_feas = 0;
	userhandle_up->ob2_still_feas = 0;
	userhandle_down = (user_data*) malloc( sizeof( user_data ) );
/*	if(!userhandle_down->x_ws) */
	userhandle_down->x_ws = calloc ((cur_numcols+1),sizeof(double));
	userhandle_down->x1 = calloc (cur_numcols,sizeof(double));
	userhandle_down->x2 = calloc (cur_numcols,sizeof(double));
	userhandle_down->ws_still_feas = 0;
	userhandle_down->ob1_still_feas = 0;
	userhandle_down->ob2_still_feas = 0;
	
	for(i=0;i<cur_numcols;i++)
	{
		userhandle_up->x_ws[i] = -100000.;
		userhandle_up->x1[i] = -100000.;
		userhandle_up->x2[i] = -100000.;
		userhandle_down->x_ws[i] = -100000.;
		userhandle_down->x1[i] = -100000.;
		userhandle_down->x2[i] = -100000.;
	}
	
	userhandle_up->prob = NULL;
	userhandle_down->prob = NULL;
	
	if(show_progress)
	{
		userhandle_up->prob = CPXcloneprob (env, nodelp, &status);
		if ( status ) {
			printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
			goto TERMINATE;
		}
		userhandle_down->prob = CPXcloneprob (env, nodelp, &status);
		if ( status ) {
			printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
			goto TERMINATE;
		}
	}
	
	for(i=0;i<total_num_integer;i++) frac_scores[i] = 0.;
   	num_frac = 0;
   	multiplier = convert_it(cur_numcols);
	
	nodelp2 = CPXcloneprob (env, nodelp, &status);
	if ( status ) {
		printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
		goto TERMINATE;
	}
	
	if(infeasible_in_cutcallback)
	{
		fathoming = 1;
		infeasible_in_cutcallback = 0;
		*useraction_p = CPX_CALLBACK_SET;
		nodecnt = 0;
  		goto TERMINATE;
	}
	
	status = CPXgetlb (env, nodelp2, lbs, obj1_index, obj2_index);
	if ( status ) {
		printf ("Failed to get lb's for objectives. Error code %d\n", status);
		goto TERMINATE;
	}
	
	status = CPXgetub (env, nodelp2, ubs, obj1_index, obj2_index);
	if ( status ) {
		printf ("Failed to get ub's for objectives. Error code %d\n", status);
		goto TERMINATE;
	}
	
	reduced_subprob_x_lb = lbs[0];
	reduced_subprob_y_ub = ubs[1];
	reduced_subprob_x_ub = ubs[0];
	reduced_subprob_y_lb = lbs[1];
	
	if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
	if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
	if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
						reduced_subprob_y_lb,reduced_subprob_y_lb);					
	if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
						reduced_subprob_y_ub,reduced_subprob_y_ub);
	
	if(! mock_insert(1,reduced_subprob_x_ub,reduced_subprob_y_ub,0,0,0,&tree) )
	{
		if(printing_in_setbranch) printf("the entire search region is dominated. Fathom.\n");
		fathoming = 1;
		*useraction_p = CPX_CALLBACK_SET;
		nodecnt = 0;
  		goto TERMINATE;
	}
	
	double slope = (ubs[1]-lbs[1])/(lbs[0]-ubs[0]);
	if(slope != slope) //slope < -1000. || slope > -.001) 
	{
		if(printing_in_setbranch) printf("setting slope to -1\n");
		slope = -1.;
	}
	if(printing_in_setbranch) printf("slope: %lf\n",slope);
   	
   	if(objective_space_fathoming)
   	{
	   	if (branch_iterations > -1) pareto_branching = 1;
		if (branch_iterations > 1000)
		{
			int nodesleft = 0;
			if ((status = CPXgetcallbackinfo (env, cbdata, wherefrom, 
			CPX_CALLBACK_INFO_NODES_LEFT, &nodesleft))) goto TERMINATE;
		
			if(nodesleft > 1000) 
			{
				pareto_branching = 0;
				if(printing_in_setbranch) printf("(%d) changing dom mid to 0\n",__LINE__);
				remove_dominated_middle = 0;
			}
			else pareto_branching = 1;
		}
	}
	
/*	printf("cplex index: %d, value: %lf\n",indices[0],bd[0]);*/
	
/*	if (branch_iterations > 100000000)*/
/*	{*/
/*		printf("**************reached iteration limit*****************\n");*/
/*		exit(0);*/
/*	}*/
   	
   	/* Get the objective value associated with this node. This is used for providing CPLEX with an estimated value of the solution an the next node. */
		
	status = CPXgetcallbacknodeobjval (env, cbdata, wherefrom, &objval);
  	if ( status ) {
      		fprintf (stdout, "Can't get node objective value.");
		goto TERMINATE;
	}
	
	status = CPXgetcallbacknodex (env, cbdata, wherefrom, x_ws, 0, cur_numcols-1);
	if ( status ) {
		printf ("CPXgetcallbacknodelp, Failed to get node x, error code %d\n", status);
		goto TERMINATE;
	}
	ws_lp_objvals[0] = x_ws[obj1_index];
	ws_lp_objvals[1] = x_ws[obj2_index];
	
	double proj1 = (slope)*(ubs[0]-x_ws[obj1_index])+x_ws[obj2_index];
	double proj2 = (1./slope)*(ubs[1]-x_ws[obj2_index])+x_ws[obj1_index];
	endpoint1_x = ubs[0];
	endpoint2_x = proj2;
	endpoint1_y = proj1;
	endpoint2_y = ubs[1];
	if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",x_ws[obj1_index],x_ws[obj2_index]);
	
/*	status = CPXwriteprob (env, nodelp, "myprob2.lp", "LP");*/
/*	exit(0);*/
	
	status = CPXgetcallbacknodeinfo(	env,
                                 		cbdata,
                                 		wherefrom,
                                 		0,
                               			CPX_CALLBACK_INFO_NODE_SEQNUM_LONG,
                               			&seqnum);
        if ( status ) {
		printf ("CPXgetcallbacknodeinfo, Failed to get seqnum, error code %d\n", status);
		goto TERMINATE;
	}
	
	if(printing_in_setbranch) printf("at branching seqnum is: %d\n",seqnum);
	
/*	if(seqnum == 0) printf("iterations: %d\n",branch_iterations);*/
	
/*	if( seqnum == 3)*/
/*	{*/
/*		status = CPXwriteprob (env, nodelp, "myprob2.lp", "LP");*/
/*		exit(0);*/
/*	}*/
	
	if(seqnum == 0 && objective_space_fathoming == 0) 
	{
		printf("(%d) changing dom mid to 0\n",__LINE__);
		remove_dominated_middle = 0;
	}
	
	if(printing_in_setbranch){
	printf("_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-\n");
	PSA_all(env,nodelp);
	printf("_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-\n");}
	
	status = CPXgetcallbacknodeinfo(	env,
                                 		cbdata,
                                 		wherefrom,
                                 		0,
                                 		CPX_CALLBACK_INFO_NODE_DEPTH,
                                 		&depth);
        if ( status ) {
		printf ("CPXgetcallbacknodeinfo, Failed to get node depth, error code %d\n", status);
		goto TERMINATE;
	}
	
	status = CPXgetcallbacknodeinfo(	env,
	                         		cbdata,
	                         		wherefrom,
	                         		0,
	                        		CPX_CALLBACK_INFO_NODE_USERHANDLE,
	                         		&userhandle_current);
	                         		
	var_ubs   = (double *) malloc (cur_numcols * sizeof (double));
	var_lbs   = (double *) malloc (cur_numcols * sizeof (double));
/*		double *x   = (double *) malloc (cur_numcols * sizeof (double));*/
	int ws_feas = 1, ob1_feas = 1, ob2_feas = 1;
	
	status = CPXgetlb (env, nodelp, var_lbs, 0, cur_numcols-1);
	if ( status ) {
		printf ("(%d) Failed to get lb's for objectives. Error code %d\n",__LINE__, status);
		exit(0);
		goto TERMINATE;
	}

	status = CPXgetub (env, nodelp, var_ubs, 0, cur_numcols-1);
	if ( status ) {
		printf ("(%d) Failed to get ub's for objectives. Error code %d\n",__LINE__, status);
		goto TERMINATE;
	}
	
/*	if(!userhandle_current) printf("userhandle_current doesn't exist\n");*/
		
	if(seqnum != 0 && userhandle_current)
	{
		int k = 0;
		if(!(userhandle_current->x_ws)) 
		{	
/*				printf("uh_cur x_ws didn't exist\n");*/
			ws_feas = 0;
		}
		if(!(userhandle_current->x1))
		{
/*				printf("uh_cur x1 didn't exist\n");*/
			ob1_feas = 0;
		}
		if(!(userhandle_current->x2))
		{
/*				printf("uh_cur x2 didn't exist\n");*/
			ob2_feas = 0;
		}
		while( (ws_feas || ob1_feas || ob2_feas) && k < cur_numcols )
		{
			
			if(ws_feas)
			{
/*					printf("k: %d, lb: %lf, x_ws: %lf, ub: %lf\n",k,var_lbs[k],userhandle_current->x_ws[k],var_ubs[k]);*/
				if(userhandle_current->x_ws[k] != userhandle_current->x_ws[k] || userhandle_current->x_ws[k] - var_lbs[k] < -.00001 
					|| userhandle_current->x_ws[k] - var_ubs[k] > .00001) ws_feas = 0;
			}
			if(ob1_feas)
			{
/*					printf("k: %d, lb: %lf, x1: %lf, ub: %lf\n",k,var_lbs[k],userhandle_current->x1[k],var_ubs[k]);*/
				if(userhandle_current->x1[k] != userhandle_current->x1[k] || userhandle_current->x1[k] - var_lbs[k] < -.00001 
					|| userhandle_current->x1[k] - var_ubs[k] > .00001) ob1_feas = 0;
			}
			if(ob2_feas)
			{
/*					printf("k: %d, lb: %lf, x2: %lf, ub: %lf\n",k,var_lbs[k],userhandle_current->x2[k],var_ubs[k]);*/
				if(userhandle_current->x2[k] != userhandle_current->x2[k] || userhandle_current->x2[k] - var_lbs[k] < -.00001 
					|| userhandle_current->x2[k] - var_ubs[k] > .00001) ob2_feas = 0;
			}
			k++;
		}
		if(ws_feas)
		{
				if(printing_in_setbranch) printf("previous ws mip soln still feasible here\n");
			userhandle_current->ws_still_feas = 1;
		}
		if(ob1_feas)
		{
				if(printing_in_setbranch) printf("previous ob1 mip soln still feasible here\n");
			userhandle_current->ob1_still_feas = 1;
		}
		if(ob2_feas)
		{
				if(printing_in_setbranch) printf("previous ob2 mip soln still feasible here\n");
			userhandle_current->ob2_still_feas = 1;
		}
	}
	
	if(printing_in_setbranch){
	printf("*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*\n");
	print_inorder(tree,1);
	printf("*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*\n");
	}
	
	*useraction_p = CPX_CALLBACK_DEFAULT;

	if(there_will_only_be_points) 
	{
		points_only = 0;
	}
	if(points_only && seqnum == 0)
	{
		clock_t current_time = clock();	
		cumulative_time = (double)(current_time - start_BB) / CLOCKS_PER_SEC;
		status = CPXsetdblparam (env_just_solve_mips, CPX_PARAM_TILIM, max_time - cumulative_time);
		if ( status ) {
		    	printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);
		   	exit(0);
		}
		goto SOLVE_OB2_MIP;
	}
	else
	{
		if(time_vs_node_lim)
		{
			status = CPXsetintparam (env_just_solve_mips, CPX_PARAM_NODELIM, 0);
		    	if ( status ) {
		    		printf ("Failure to set MIP node limit to 0, error %d.\n",status);
		    		exit(0);
		    	}
		}
		else
		{
			status = CPXsetdblparam (env_just_solve_mips, CPX_PARAM_TILIM, time_limit);
			if ( status ) {
			    	printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);
			   	exit(0);
			}
		}
	}

	nodelp_copy = CPXcloneprob (env, nodelp2, &status);
	if ( status ) {
		printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
		goto TERMINATE;
	}
	
	if(userhandle_current && userhandle_current->ws_still_feas) goto SOLVE_WS_MIP;
	
	chg_coefs(env, nodelp_copy, indexes, slope);
	
/*	if(x_ws[obj1_index] > reduced_subprob_x_lb && x_ws[obj2_index] < reduced_subprob_y_ub && */
/*				x_ws[obj1_index] < reduced_subprob_x_ub && x_ws[obj2_index] > reduced_subprob_y_lb) */
/*	{*/
/*		ws_sol_interior = 1;*/
/*		slope = initial_slope;*/
/*		goto CHECK_FEASIBILITY;*/
/*	}*/
	
	status = CPXlpopt (env, nodelp_copy);
 	if ( status ) {
   		printf ("%s(%d): CPXlpopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
		goto TERMINATE;
	}
	
	lpstat = CPXgetstat (env, nodelp_copy);
	if(lpstat == 3 || lpstat == 4)
	{
/*		printf("the ws lp is infeasible, fathom\n");*/
		fathoming = 1;
		*useraction_p = CPX_CALLBACK_SET;
		nodecnt = 0;
  		goto TERMINATE;
	}
	
  	status = CPXgetx (env, nodelp_copy, x_ws, 0, cur_numcols-1);
	if(status) 
	{
		printf ("(%d) CPXgetx, Failed to get x values, error code %d\n", __LINE__,status);
		goto TERMINATE;
	}
	ws_lp_objvals[0] = x_ws[obj1_index];
	ws_lp_objvals[1] = x_ws[obj2_index];
	endpoint1_x = x_ws[obj1_index];
	endpoint2_x = x_ws[obj1_index];
	endpoint1_y = x_ws[obj2_index];
	endpoint2_y = x_ws[obj2_index];
	
	if(x_ws[obj1_index] > reduced_subprob_x_lb && x_ws[obj2_index] < reduced_subprob_y_ub && 
				x_ws[obj1_index] < reduced_subprob_x_ub && x_ws[obj2_index] > reduced_subprob_y_lb) ws_sol_interior = 1;
	else
	{
		if(x_ws[obj2_index] >= reduced_subprob_y_ub) 
		{
/*			ob2_lp_been_solved = 1;*/
			ws_lp_sol_on_top_bd  = 1; 
		}
		if(x_ws[obj1_index] >= reduced_subprob_x_ub) 
		{
/*			ob1_lp_been_solved = 1;*/
			ws_lp_sol_on_right_bd = 1;
		}
	}
	
	/*************** Here we check to see if the solutions found at the parent node are still feasible at this one ***************************/
	
	CHECK_FEASIBILITY:
	
	if(control_node_selection)
	{
		userhandle_up->f1 = x_ws[obj1_index];
		userhandle_down->f1 = x_ws[obj1_index];
		userhandle_up->f2 = x_ws[obj2_index];
		userhandle_down->f2 = x_ws[obj2_index];
	}
	
	if(printing_in_setbranch) printf("cplex index: %d\n",indices[0]);
	
	for(i=0;i<total_num_integer;i++)
	{
		k = integer_indices[i];
	  	diff = x_ws[k] - floor(x_ws[k]);
/*		printf("diff%d: %lf\n",i,diff);*/
  		if( diff >= .00001 && diff <= .99999)
  		{
  			if(printing_in_setbranch) printf("changing frac index to %d (%d)\n",integer_indices[i],__LINE__);
/*  			printf("plot(%lf,%lf,'go');\n",x_ws[obj1_index],x_ws[obj2_index]);*/
	  		frac_index = k;
	  		frac_val = x_ws[k];
			if(frac_scores[i] > 0.0001) frac_scores[i] += 1.;
			else 
			{
				frac_scores[i] += 1. + multiplier*k;
				num_frac++;
			}
			frac_values[k] = frac_val;
		}
	}
	if(num_frac == 0)
	{
		if(printing_in_setbranch) printf("ws lp solution is integer feasible\n");
/*		for(i=0;i<total_num_integer;i++) printf("x%d: %lf\n",integer_indices[i],x_ws[integer_indices[i]]);*/
		all_feas = 0;
		
		ws_mip_opt = 1;
		ws_mip_objvals[0] = x_ws[obj1_index];
		ws_mip_objvals[1] = x_ws[obj2_index];
		
		if(!userhandle_up->x_ws) userhandle_up->x_ws = calloc ((cur_numcols+1),sizeof(double));
		if(!userhandle_down->x_ws) userhandle_down->x_ws = calloc ((cur_numcols+1),sizeof(double));
	
		for(i=0;i<cur_numcols;i++) 
		{
			userhandle_up->x_ws[i] = x_ws[i];
			userhandle_down->x_ws[i] = x_ws[i];
		}
		userhandle_up->x_ws[cur_numcols] = slope;
		userhandle_down->x_ws[cur_numcols] = slope;
		
		ws_lp_int_feas = 1;
	}
	
	if(indices[0] > integer_indices[0] && indices[0] < cur_numcols)
	{
		k = 0;
		for(i=0;i<total_num_integer;i++)
		{
			k = integer_indices[i];
			if(k == indices[0])
			{
				k = i;
				break;
			}
		}
	
/*		printf("cplex index: %d, cols: %d, k: %d, i: %d\n",indices[0],cur_numcols,k,i);*/
		if(frac_scores[k] > .0001) frac_scores[k] += 1.;
		else frac_scores[k] += 1. + multiplier*indices[0];
		frac_values[indices[0]] = x_ws[indices[0]];
		num_frac++;
	}
	
	/*************** We now start to process the weighted sum LP solution ****************************/

	if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",x_ws[obj1_index],x_ws[obj2_index]);
	add_check = mock_insert(1,x_ws[obj1_index],x_ws[obj2_index],0,0,0,&tree);
	if(add_check)
	{
		if(printing_in_setbranch) printf("weighted sum solution not dominated\n");
		if(ws_lp_int_feas)
		{
			if(printing_in_setbranch) printf("weighted sum solution is integer feasible\n");
			
			ws_mip_opt = 1;
			ws_mip_objvals[0] = x_ws[obj1_index];
			ws_mip_objvals[1] = x_ws[obj2_index];
			
			add_check = mock_insert(1,x_ws[obj1_index],x_ws[obj2_index],0,0,0,&tree);
			if(add_check)
			{
				if(printing_in_setbranch) printf("adding solution\n");
				for(i=0;i<cur_numcols;i++)
		      		{
		      			stored_x[x_rotation][i] = x_ws[i];
		      		}
		      		x_rotation = (x_rotation + 1) % num_x_to_store;
		      		add_check = 0;
		      		PSA_full(env,NULL,x_ws,NULL,NULL);
			}
			
			if(!ws_sol_interior) 
			{
				status = CPXgetbase (env, nodelp_copy, cstat_ws, rstat_ws);
				if(status) 
				{
					printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
					goto TERMINATE;
				}
			}
			else 
			{
				CPXgetbase (env, nodelp, cstat_ws, rstat_ws);
				if(status) 
				{
					printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
					goto TERMINATE;
				}
			}
			if(printing_in_setbranch) printf("status: %d\n",status);
			
	  		sub_pr1_x_ub = x_ws[obj1_index];
			sub_pr1_y_lb = x_ws[obj2_index];
			sub_pr2_x_lb = x_ws[obj1_index];
			sub_pr2_y_ub = x_ws[obj2_index];
		  	PSA_right_check = PSA(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, seqnum, 0, nodelp_copy);
		  	PSA_left_check = PSA_left(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, seqnum, 0, nodelp_copy);
		  	if(within_PSA_score > 1.) within_PSA_score -= 1.;
		  	
			if(PSA_right_check == 2 && PSA_left_check == 2)
	  		{
		  		if(printing_in_setbranch) printf("fathoming node %d for PSA completion (%d)\n",seqnum,__LINE__);
			  	fathomed_by_PSA_completion++;
				fathoming = 1;
				*useraction_p = CPX_CALLBACK_SET;
				nodecnt = 0;
		  		goto TERMINATE;
	  		}
		  	else if(PSA_right_check == 2)
		  	{
				if(printing_in_setbranch) printf("the left subproblem is now empty by PSA completion\n");
		  		reduced_subprob_x_lb = sub_pr2_x_lb;
		  		reduced_subprob_y_ub = sub_pr2_y_ub;
		  		reduced_subprob_x_ub = ubs[0];
				reduced_subprob_y_lb = lbs[1];
				if(printing_in_setbranch) printf("after changing:\n");
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_ub,reduced_subprob_y_ub);
				bound_reduction = 1;
				left_side_dom = 1;
		  		goto SOLVE_OB1_LP;
		  	}
		  	else if(PSA_left_check == 2)
			{
	  			if(printing_in_setbranch) printf("the right subproblem is now empty by PSA completion\n");
		  		reduced_subprob_x_lb = lbs[0];
		  		reduced_subprob_y_ub = ubs[1];
				reduced_subprob_x_ub = sub_pr1_x_ub;
	  			reduced_subprob_y_lb = sub_pr1_y_lb;
	  			if(printing_in_setbranch) printf("after changing:\n");
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_ub,reduced_subprob_y_ub);
	  			bound_reduction = 1;
	  			right_side_dom = 1;
		  		goto SOLVE_OB2_LP;
		  	}
		  	else if(pareto_branching && sub_pr2_x_lb > sub_pr1_x_ub)
			{
		  		goto PARETO_BRANCH;
		  	}
		  	else
		  	{
		  		if(printing_in_setbranch) printf("none of the above\n");
		  	}
			
			if(!ob1_lp_been_solved) goto SOLVE_OB1_LP;
			else if(!ob1_lp_int_feas && !ob1_mip_ran) goto SOLVE_OB1_MIP;
			else goto BRANCHING;
		}
		else
		{
			if(printing_in_setbranch) printf("weighted sum solution is not integer feasible\n");
			
			/****************************************************************************
				 This is a scheme for detecting if solutions are all singletons in the 
			         objective space. We then branch in the objective space 
			****************************************************************************/
			
/*			DO_IT_AGAIN:*/
/*			;*/
/*			closest_nodes *two_nodes = find_two_nodes_right_of_val(reduced_subprob_x_lb, reduced_subprob_y_ub, tree);*/
/*			if(!two_nodes) */
/*			{*/
/*				if(points_only || its_been_only_points) goto SOLVE_OB2_MIP;*/
/*				else goto AFTER_THIS1;*/
/*			}*/
/*			if(printing_in_setbranch) printf("the two nodes:\n");*/
/*			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->closest->nw_x,*/
/*				x_ideal-two_nodes->closest->se_x,y_ideal-two_nodes->closest->nw_y,y_ideal-two_nodes->closest->se_y);*/
/*			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->next->nw_x,x_ideal-two_nodes->next->se_x,*/
/*				y_ideal-two_nodes->next->nw_y,y_ideal-two_nodes->next->se_y);*/
/*			if(two_nodes->closest->type == 2 && fabs(two_nodes->closest->nw_x - two_nodes->closest->se_x) < .00000001 && */
/*				fabs(two_nodes->closest->nw_y - two_nodes->closest->se_y) < .00000001) two_nodes->closest->type = 1;*/
/*			if(two_nodes->next->type == 2 && fabs(two_nodes->next->nw_x - two_nodes->next->se_x) < .00000001 && */
/*				fabs(two_nodes->next->nw_y - two_nodes->next->se_y) < .00000001) two_nodes->next->type = 1;*/
/*			*/
/*			if(two_nodes->closest->type == 1 && ((fabs(two_nodes->closest->nw_x - two_nodes->next->nw_x) < .0000001 && */
/*				two_nodes->closest->nw_y - two_nodes->next->nw_y >= -.0000001) || */
/*				(fabs(two_nodes->closest->nw_x - two_nodes->next->se_x) < .0000001 && */
/*				two_nodes->closest->nw_y - two_nodes->next->se_y >= -.0000001)))*/
/*			{*/
/*				if(printing_in_setbranch) */
/*					printf("the closest node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
/*				delete_node(two_nodes->closest);*/
/*				free(two_nodes);*/
/*				goto DO_IT_AGAIN;*/
/*			}*/
/*			*/
/*			if(x_ideal - two_nodes->closest->nw_x - reduced_subprob_x_ub > .0001 || */
/*				y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_lb < -.0001 || */
/*				x_ideal - two_nodes->closest->nw_x - reduced_subprob_x_lb < -.0001 || */
/*				y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_ub > .0001)*/
/*			{*/
/*				if(printing_in_setbranch) printf("closest node is outside region\n");*/
/*				free(two_nodes);*/
/*				goto SOLVE_WS_MIP;*/
/*			}*/
/*			*/
/*			else if(two_nodes->next->type == 1 && ((fabs(two_nodes->next->nw_x - two_nodes->closest->nw_x) < .0000001 && */
/*				two_nodes->next->nw_y - two_nodes->closest->nw_y >= -.0000001) || */
/*				(fabs(two_nodes->next->nw_x - two_nodes->closest->se_x) < .0000001 && */
/*				two_nodes->next->nw_y - two_nodes->closest->se_y >= -.0000001)))*/
/*			{*/
/*				if(printing_in_setbranch) */
/*					printf("the next node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
/*				delete_node(two_nodes->next);*/
/*				free(two_nodes);*/
/*				goto DO_IT_AGAIN;*/
/*			}*/
/*			*/
/*			if(!there_will_only_be_points && two_nodes->closest->type == 1 && two_nodes->next->type == 1)*/
/*			{*/
/*				points_only = 1;*/
/*				its_been_only_points = 1;*/
/*				status = CPXsetdblparam (env_just_solve_mips, CPX_PARAM_TILIM, pow(10.,75.));*/
/*				if ( status ) {*/
/*				    	printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);*/
/*				   	exit(0);*/
/*				}*/
/*				if(x_ideal - two_nodes->next->se_x < x_ideal - two_nodes->closest->nw_x )*/
/*				{*/
/*					if(y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_lb > .0001)*/
/*					{*/
/*						printf("there was only one node right of the val. Pareto branch based on its location\n");*/
/*						sub_pr1_x_ub = reduced_subprob_x_ub;*/
/*						sub_pr1_y_lb = y_ideal - two_nodes->closest->nw_y;*/
/*					 	sub_pr2_x_lb = x_ideal - two_nodes->closest->nw_x + .001;*/
/*						sub_pr2_y_ub = y_ideal - two_nodes->closest->nw_y;*/
/*						free(two_nodes);*/
/*						goto PARETO_BRANCH;*/
/*					}*/
/*				}*/
/*				else if(x_ideal - two_nodes->next->se_x - reduced_subprob_x_ub <= .0000001 && */
/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001)*/
/*				{*/
/*					double ran = (double) rand() / ( (double) RAND_MAX);*/
/*					if(printing_in_setbranch) printf("left pt of next node is also inside the search region.\n");*/
/*					if(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y) < -.0001 && */
/*						y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001 && */
/*						y_ideal - two_nodes->next->se_y - reduced_subprob_y_ub <= .0000001)*/
/*					{*/
/*						if(printing_in_setbranch) printf("separation between y-values, split\n");*/
/*						printf("y separation is %lf percent of y_range\n", */
/*							100.*(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y))/y_range);*/
/*						sub_pr1_x_ub = reduced_subprob_x_ub;*/
/*						sub_pr1_y_lb = y_ideal - two_nodes->closest->nw_y;*/
/*					 	sub_pr2_x_lb = x_ideal - two_nodes->closest->nw_x + .001;*/
/*						sub_pr2_y_ub = y_ideal - two_nodes->closest->nw_y;*/
/*						free(two_nodes);*/
/*						goto PARETO_BRANCH;*/
/*					}*/
/*					else if(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x) > .0001)*/
/*					{*/
/*						if(printing_in_setbranch) printf("separation between x-values, split\n");*/
/*						printf("x separation is %lf percent of x_range\n", */
/*							100.*(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x))/x_range);*/
/*						exit(0);*/
/*						sub_pr1_x_ub = ((x_ideal - two_nodes->next->se_x) + (x_ideal - two_nodes->closest->nw_x))/2.;*/
/*						sub_pr1_y_lb = y_ideal - two_nodes->next->se_y + .001;//reduced_subprob_y_lb;*/
/*					 	sub_pr2_x_lb = ((x_ideal - two_nodes->next->se_x) + (x_ideal - two_nodes->closest->nw_x))/2.;*/
/*						sub_pr2_y_ub = reduced_subprob_y_ub;*/
/*						free(two_nodes);*/
/*						goto PARETO_BRANCH;*/
/*					}*/
/*				}*/
/*				else*/
/*				{*/
/*					printf("the 2nd node was out of the region. Go back to beginning\n");*/
/*					ob2_mip_ran = 0;*/
/*					ob2_mip_opt = 0;*/
/*					goto BEGINNING;*/
/*				}*/
/*			}*/
/*			else */
/*			{*/
/*				points_only = 0;*/
/*				its_been_only_points = 0;*/
/*				status = CPXsetdblparam (env_just_solve_mips, CPX_PARAM_TILIM, time_limit);*/
/*				if ( status ) {*/
/*				    	printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);*/
/*				   	exit(0);*/
/*				}*/
/*				*/
/*				goto AFTER_THIS1;*/
/*				*/
/*				if(x_ideal - two_nodes->next->se_x < x_ideal - two_nodes->closest->nw_x )*/
/*				{*/
/*					if(y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_lb > .0001)*/
/*					{*/
/*						goto AFTER_THIS1;*/
/*					}*/
/*				}*/
/*				else if(x_ideal - two_nodes->next->se_x - reduced_subprob_x_ub <= .0000001 && */
/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001)*/
/*				{*/
/*					double ran = (double) rand() / ( (double) RAND_MAX);*/
/*					if(printing_in_setbranch) printf("left pt of next node is also inside the search region.\n");*/
/*					if( (y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y) < -.0001 && */
/*						y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001 && */
/*						y_ideal - two_nodes->next->se_y - reduced_subprob_y_ub <= .0000001) && */
/*						(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x) > .0001) )*/
/*					{*/
/*						if(printing_in_setbranch) printf("separation between x and y-values, split\n");*/
/*						double y_sep = fabs(100.*(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y))/y_range);*/
/*						double x_sep = fabs(100.*(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x))/x_range);*/
/*						printf("y separation is %lf percent of y_range\n", */
/*							100.*(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y))/y_range);*/
/*						printf("x separation is %lf percent of x_range\n", */
/*							100.*(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x))/x_range);*/
/*						if(x_sep > 5. || y_sep > 5.)*/
/*						{*/
/*							printf("(%d) exploiting separation\n",__LINE__);*/
/*							sub_pr1_x_ub = reduced_subprob_x_ub;*/
/*							sub_pr1_y_lb = y_ideal - two_nodes->closest->nw_y;*/
/*						 	sub_pr2_x_lb = x_ideal - two_nodes->closest->nw_x + .001;*/
/*							sub_pr2_y_ub = y_ideal - two_nodes->closest->nw_y;*/
/*							free(two_nodes);*/
/*							goto PARETO_BRANCH;*/
/*						}*/
/*						else goto AFTER_THIS1;*/
/*					}*/
/*				}*/
/*			}*/
/*			AFTER_THIS1:*/
/*			*/
/*			free(two_nodes);*/
/*			goto SOLVE_WS_MIP;*/
			goto BRANCHING;
		}
	}
	else
	{
		if(printing_in_setbranch) printf("weighted sum solution dominated\n");

		ws_lp_dom = 1;
		
		if(ws_lp_int_feas)
		{
			if(!ws_sol_interior) 
			{
				status = CPXgetbase (env, nodelp_copy, cstat_ws, rstat_ws);
				if(status) 
				{
					printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
					goto TERMINATE;
				}
			}
			else 
			{
				CPXgetbase (env, nodelp, cstat_ws, rstat_ws);
				if(status) 
				{
					printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
					goto TERMINATE;
				}
			}
			if(printing_in_setbranch) printf("status: %d\n",status);
			
			sub_pr1_x_ub = x_ws[obj1_index];
			sub_pr1_y_lb = x_ws[obj2_index];
			sub_pr2_x_lb = x_ws[obj1_index];
			sub_pr2_y_ub = x_ws[obj2_index];
		  	PSA_right_check = PSA(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, seqnum, 0, nodelp_copy);
		  	PSA_left_check = PSA_left(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, seqnum, 0, nodelp_copy);
		  	
			if(PSA_right_check == 2 && PSA_left_check == 2)
	  		{
		  		if(printing_in_setbranch) printf("fathoming node %d for PSA completion (%d)\n",seqnum,__LINE__);
			  	fathomed_by_PSA_completion++;
				fathoming = 1;
				*useraction_p = CPX_CALLBACK_SET;
				nodecnt = 0;
		  		goto TERMINATE;
	  		}
		  	else if(PSA_right_check == 2)
		  	{
				if(printing_in_setbranch) printf("the left subproblem is now empty by PSA completion\n");
		  		reduced_subprob_x_lb = x_ws[obj1_index];
		  		if(there_will_only_be_points && integer_bb && integer_objective == 1) reduced_subprob_x_lb += smallest_coef;
		  		reduced_subprob_y_ub = x_ws[obj2_index];
		  		if(there_will_only_be_points && integer_bb && integer_objective == 2) reduced_subprob_y_ub -= smallest_coef;
		  		reduced_subprob_x_ub = ubs[0];
				reduced_subprob_y_lb = lbs[1];
				if(printing_in_setbranch) printf("after changing:\n");
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_ub,reduced_subprob_y_ub);
				bound_reduction = 1;
				left_side_dom = 1;
		  		goto SOLVE_OB1_LP;
		  	}
		  	else if(PSA_left_check == 2)
			{
	  			if(printing_in_setbranch) printf("the right subproblem is now empty by PSA completion\n");
		  		reduced_subprob_x_lb = lbs[0];
		  		reduced_subprob_y_ub = ubs[1];
				reduced_subprob_x_ub = x_ws[obj1_index];
				if(there_will_only_be_points && integer_bb && integer_objective == 1) reduced_subprob_x_ub -= smallest_coef;
	  			reduced_subprob_y_lb = x_ws[obj2_index];
	  			if(there_will_only_be_points && integer_bb && integer_objective == 2) reduced_subprob_y_lb += smallest_coef;
	  			if(printing_in_setbranch) printf("after changing:\n");
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_ub,reduced_subprob_y_ub);
	  			bound_reduction = 1;
	  			right_side_dom = 1;
		  		goto SOLVE_OB2_LP;
		  	}
		  	else if(pareto_branching && sub_pr2_x_lb > sub_pr1_x_ub)
			{
		  		goto PARETO_BRANCH;
		  	}
		  	else
		  	{
		  		if(printing_in_setbranch) printf("none of the above\n");
/*		  		for(i=0;i<total_num_integer;i++)*/
/*				{*/
/*					k = integer_indices[i];*/
/*					if(k == frac_index) */
/*					{*/
/*						k = i;*/
/*						break;*/
/*					}*/
/*				}*/
/*				frac_scores[k] += 3.;*/
/*				frac_values[frac_index] = frac_val;*/
		  	}
		}
		else if(remove_dominated_middle)
		{
			if(!ws_sol_interior) 
			{
				status = CPXgetbase (env, nodelp_copy, cstat_ws, rstat_ws);
				if(status) 
				{
					printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
					goto TERMINATE;
				}
			}
			else 
			{
				CPXgetbase (env, nodelp, cstat_ws, rstat_ws);
				if(status) 
				{
					printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
					goto TERMINATE;
				}
			}
			if(printing_in_setbranch) printf("status: %d\n",status);

			sub_pr1_x_ub = x_ws[obj1_index];
			sub_pr1_y_lb = x_ws[obj2_index];
			sub_pr2_x_lb = x_ws[obj1_index];
			sub_pr2_y_ub = x_ws[obj2_index];
		  	PSA_left_check = PSA_reduce_left(env, nodelp2, x_ws, cstat_ws, rstat_ws, indexes);
		  	PSA_right_check = PSA_reduce_right(env, nodelp2, x_ws, cstat_ws, rstat_ws, indexes, seqnum);
		  	
		  	if(PSA_left_check == 2 && PSA_right_check == 2)
		  	{
		  		fathomed_by_dominated_lb++;
		  		fathoming = 1;
		  		goto TERMINATE;
		  	}
			if(PSA_left_check == 2)
			{
	  			if(printing_in_setbranch) printf("the right subproblem is now empty by PSA completion\n");
				reduced_subprob_x_ub = x_ws[obj1_index];
	  			reduced_subprob_y_lb = x_ws[obj2_index];
	  			if(printing_in_setbranch) printf("after changing:\n");
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_ub,reduced_subprob_y_ub);
	  			bound_reduction = 1;
	  			right_side_dom = 1;
		  		goto SOLVE_OB2_LP;
		  	}
		  	if(PSA_right_check == 2)
			{
	  			if(printing_in_setbranch) printf("the right subproblem is now empty by PSA completion\n");
				reduced_subprob_x_lb = x_ws[obj1_index];
	  			reduced_subprob_y_ub = x_ws[obj2_index];
	  			if(printing_in_setbranch) printf("after changing:\n");
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_ub,reduced_subprob_y_ub);
	  			bound_reduction = 1;
	  			left_side_dom = 1;
		  		goto SOLVE_OB1_LP;
		  	}
		  	if(sub_pr2_x_lb > sub_pr1_x_ub)
			{
				if(pareto_branching) goto PARETO_BRANCH;
		  	}
		}
		
/*		DO_IT_AGAIN3:*/
/*		;*/
/*		closest_nodes *two_nodes = find_two_nodes_right_of_val(reduced_subprob_x_lb, reduced_subprob_y_ub, tree);*/
/*		if(!two_nodes) */
/*		{*/
/*			if(points_only || its_been_only_points) goto SOLVE_OB2_MIP;*/
/*			else goto AFTER_THIS2;*/
/*		}*/
/*		if(printing_in_setbranch) printf("the two nodes:\n");*/
/*		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->closest->nw_x,*/
/*			x_ideal-two_nodes->closest->se_x,y_ideal-two_nodes->closest->nw_y,y_ideal-two_nodes->closest->se_y);*/
/*		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->next->nw_x,x_ideal-two_nodes->next->se_x,*/
/*			y_ideal-two_nodes->next->nw_y,y_ideal-two_nodes->next->se_y);*/
/*		if(two_nodes->closest->type == 2 && fabs(two_nodes->closest->nw_x - two_nodes->closest->se_x) < .00000001 && */
/*			fabs(two_nodes->closest->nw_y - two_nodes->closest->se_y) < .00000001) two_nodes->closest->type = 1;*/
/*		if(two_nodes->next->type == 2 && fabs(two_nodes->next->nw_x - two_nodes->next->se_x) < .00000001 && */
/*			fabs(two_nodes->next->nw_y - two_nodes->next->se_y) < .00000001) two_nodes->next->type = 1;*/
/*		if(two_nodes->closest->type == 1 && ((fabs(two_nodes->closest->nw_x - two_nodes->next->nw_x) < .0000001 && */
/*			two_nodes->closest->nw_y - two_nodes->next->nw_y >= -.0000001) || */
/*			(fabs(two_nodes->closest->nw_x - two_nodes->next->se_x) < .0000001 && */
/*			two_nodes->closest->nw_y - two_nodes->next->se_y >= -.0000001)))*/
/*		{*/
/*			if(printing_in_setbranch) */
/*				printf("the closest node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
/*			delete_node(two_nodes->closest);*/
/*			free(two_nodes);*/
/*			goto DO_IT_AGAIN3;*/
/*		}*/
/*		*/
/*		if(x_ideal - two_nodes->closest->nw_x - reduced_subprob_x_ub > .0001 || */
/*			y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_lb < -.0001 || */
/*			x_ideal - two_nodes->closest->nw_x - reduced_subprob_x_lb < -.0001 || */
/*			y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_ub > .0001)*/
/*		{*/
/*			if(printing_in_setbranch) printf("closest node is outside region\n");*/
/*			free(two_nodes);*/
/*			goto SOLVE_WS_MIP;*/
/*		}*/
/*		*/
/*		else if(two_nodes->next->type == 1 && ((fabs(two_nodes->next->nw_x - two_nodes->closest->nw_x) < .0000001 && */
/*			two_nodes->next->nw_y - two_nodes->closest->nw_y >= -.0000001) || */
/*			(fabs(two_nodes->next->nw_x - two_nodes->closest->se_x) < .0000001 && */
/*			two_nodes->next->nw_y - two_nodes->closest->se_y >= -.0000001)))*/
/*		{*/
/*			if(printing_in_setbranch) */
/*				printf("the next node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
/*			delete_node(two_nodes->next);*/
/*			free(two_nodes);*/
/*			goto DO_IT_AGAIN3;*/
/*		}*/
/*		*/
/*		if(!there_will_only_be_points && two_nodes->closest->type == 1 && two_nodes->next->type == 1)*/
/*		{*/
/*			points_only = 1;*/
/*			its_been_only_points = 1;*/
/*			status = CPXsetdblparam (env_just_solve_mips, CPX_PARAM_TILIM, pow(10.,75.));*/
/*			if ( status ) {*/
/*			    	printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);*/
/*			   	exit(0);*/
/*			}*/
/*			if(x_ideal - two_nodes->next->se_x < x_ideal - two_nodes->closest->nw_x)*/
/*			{*/
/*				if(y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_lb > .0001)*/
/*				{*/
/*					printf("there was only one node right of the val. Pareto branch based on its location\n");*/
/*					sub_pr1_x_ub = reduced_subprob_x_ub;*/
/*					sub_pr1_y_lb = y_ideal - two_nodes->closest->nw_y;*/
/*				 	sub_pr2_x_lb = x_ideal - two_nodes->closest->nw_x + .001;*/
/*					sub_pr2_y_ub = y_ideal - two_nodes->closest->nw_y;*/
/*					free(two_nodes);*/
/*					goto PARETO_BRANCH;*/
/*				}*/
/*			}*/
/*			else if(x_ideal - two_nodes->next->se_x - reduced_subprob_x_ub <= .0000001 && */
/*				y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001)*/
/*			{*/
/*				if(printing_in_setbranch) printf("left pt of next node is also inside the search region.\n");*/
/*				if(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y) < -.0001 && */
/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001 && */
/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_ub <= .0000001)*/
/*				{*/
/*					if(printing_in_setbranch) printf("separation between y-values, split\n");*/
/*					printf("y separation is %lf percent of y_range\n", */
/*							100.*(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y))/y_range);*/
/*					sub_pr1_x_ub = reduced_subprob_x_ub;*/
/*					sub_pr1_y_lb = y_ideal - two_nodes->closest->nw_y;*/
/*				 	sub_pr2_x_lb = x_ideal - two_nodes->closest->nw_x + .001;*/
/*					sub_pr2_y_ub = y_ideal - two_nodes->closest->nw_y;*/
/*					free(two_nodes);*/
/*					goto PARETO_BRANCH;*/
/*				}*/
/*				else if(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x) > .0001)*/
/*				{*/
/*					if(printing_in_setbranch) printf("separation between x-values, split\n");*/
/*					printf("y separation is %lf percent of y_range\n", */
/*							100.*(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x))/x_range);*/
/*					sub_pr1_x_ub = ((x_ideal - two_nodes->next->se_x) + (x_ideal - two_nodes->closest->nw_x))/2.;*/
/*					sub_pr1_y_lb = y_ideal - two_nodes->next->se_y + .001;//reduced_subprob_y_lb;*/
/*				 	sub_pr2_x_lb = ((x_ideal - two_nodes->next->se_x) + (x_ideal - two_nodes->closest->nw_x))/2.;*/
/*					sub_pr2_y_ub = reduced_subprob_y_ub;*/
/*					free(two_nodes);*/
/*					goto PARETO_BRANCH;*/
/*				}*/
/*			}*/
/*		}*/
/*		else */
/*		{*/
/*			points_only = 0;*/
/*			its_been_only_points = 0;*/
/*			status = CPXsetdblparam (env_just_solve_mips, CPX_PARAM_TILIM, time_limit);*/
/*			if ( status ) {*/
/*			    	printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);*/
/*			   	exit(0);*/
/*			}*/
/*			*/
/*			goto AFTER_THIS2;*/
/*			*/
/*			if(x_ideal - two_nodes->next->se_x < x_ideal - two_nodes->closest->nw_x )*/
/*			{*/
/*				if(y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_lb > .0001)*/
/*				{*/
/*					goto AFTER_THIS2;*/
/*				}*/
/*			}*/
/*			else if(x_ideal - two_nodes->next->se_x - reduced_subprob_x_ub <= .0000001 && */
/*				y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001)*/
/*			{*/
/*					double ran = (double) rand() / ( (double) RAND_MAX);*/
/*				if(printing_in_setbranch) printf("left pt of next node is also inside the search region.\n");*/
/*				if( (y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y) < -.0001 && */
/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001 && */
/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_ub <= .0000001) && */
/*					(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x) > .0001) )*/
/*				{*/
/*					if(printing_in_setbranch) printf("separation between x and y-values, split\n");*/
/*					double y_sep = fabs(100.*(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y))/y_range);*/
/*					double x_sep = fabs(100.*(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x))/x_range);*/
/*					printf("y separation is %lf percent of y_range\n", */
/*						100.*(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y))/y_range);*/
/*					printf("x separation is %lf percent of x_range\n", */
/*						100.*(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x))/x_range);*/
/*					if(x_sep > 5. || y_sep > 5.)*/
/*					{*/
/*						printf("(%d) exploiting separation\n",__LINE__);*/
/*						sub_pr1_x_ub = reduced_subprob_x_ub;*/
/*						sub_pr1_y_lb = y_ideal - two_nodes->closest->nw_y;*/
/*					 	sub_pr2_x_lb = x_ideal - two_nodes->closest->nw_x + .001;*/
/*						sub_pr2_y_ub = y_ideal - two_nodes->closest->nw_y;*/
/*						free(two_nodes);*/
/*						goto PARETO_BRANCH;*/
/*					}*/
/*					else goto AFTER_THIS2;*/
/*				}*/
/*			}*/
/*		}*/
/*		AFTER_THIS2:*/
/*		*/
/*		free(two_nodes);*/
		goto SOLVE_OB1_LP;
		
	}
	
	goto SOLVE_OB1_LP;
	
	/*************** Here we solve the weighted sum single objective MIP if its determined that we should ******/

	SOLVE_WS_MIP:

/*	right_side_dom = 0;*/
/*	left_side_dom = 0;*/
	ws_mip_ran = 1;
	
	if(bound_reduction)
	{
		int ind4[4] = {obj1_index,obj1_index,obj2_index,obj2_index};
		char lu4[4] = {'L','U','L','U'};
		double bds4[4] = {reduced_subprob_x_lb,reduced_subprob_x_ub,reduced_subprob_y_lb,reduced_subprob_y_ub};
		status = CPXchgbds (env, nodelp2, 4, ind4, lu4, bds4);
		if ( status ) {
	   		printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
	  	}
	}
	
	nodelp_copy2 = CPXcloneprob (env_just_solve_mips, nodelp2, &status);
	if ( status ) {
		printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
		goto TERMINATE;
	}
	CPXchgprobtype(env_just_solve_mips, nodelp_copy2, CPXPROB_MILP);
	status = CPXchgctype(env_just_solve_mips, nodelp_copy2, cur_numcols, indexes, xctype);
	
/*	if(exact_mips) goto SKIP_THIS4;*/
	if(userhandle_current && userhandle_current->ws_still_feas)
	{
		if(printing_in_setbranch) printf("ws soln from parent is still feasible here\n");
		for(j=0;j<cur_numcols;j++) x_ws[j] = userhandle_current->x_ws[j];
		ws_mip_objvals[0] = x_ws[obj1_index];
		ws_mip_objvals[1] = x_ws[obj2_index];
		slope = userhandle_current->x_ws[cur_numcols];
		
		if(!userhandle_up->x_ws) userhandle_up->x_ws = calloc ((cur_numcols+1),sizeof(double));
		if(!userhandle_down->x_ws) userhandle_down->x_ws = calloc ((cur_numcols+1),sizeof(double));
	
		for(i=0;i<cur_numcols;i++) 
		{
			userhandle_up->x_ws[i] = x_ws[i];
			userhandle_down->x_ws[i] = x_ws[i];
		}
		userhandle_up->x_ws[cur_numcols] = slope;
		userhandle_down->x_ws[cur_numcols] = slope;
		if(control_node_selection)
		{
			userhandle_up->f1 = x_ws[obj1_index];
			userhandle_down->f1 = x_ws[obj1_index];
			userhandle_up->f2 = x_ws[obj2_index];
			userhandle_down->f2 = x_ws[obj2_index];
		}
	
		if(printing_in_setbranch) printf("plotting\n");
/*		if(printing_in_setbranch) */
/*		printf("plot(%lf,%lf,'mo');\n",x_ws[obj1_index],x_ws[obj2_index]);*/
	
		if(printing_in_setbranch) printf("the mip solution therefore stays optimal\n");
		ws_mip_opt = 1;
		mip_solved = 1;
		
		if(going_back) goto BRANCHING;
		
		if(!ob1_lp_been_solved) goto SOLVE_OB1_LP;
		else if(!ob1_lp_int_feas) goto SOLVE_OB1_MIP;
		else goto BRANCHING;
	}
	
/*	printf("slope before any changes: %lf\n",slope);*/
		
	if(bound_reduction) slope = (reduced_subprob_y_ub - reduced_subprob_y_lb)/(reduced_subprob_x_lb - reduced_subprob_x_ub);
		
	if(slope != prev_slope || same_seq_flag != 1)
	{
		if(printing_in_setbranch) printf("slope: %lf\n",slope);
		prev_slope = slope;
		chg_coefs(env_just_solve_mips, nodelp_copy2, indexes, slope);
		if(ubs[0]-lbs[0] < prob_width/width_divisor) pareto_branching = 0;
	}
	
	if(!mip_solved) 
	{
		if(printing_in_setbranch) printf("here: %d\n",__LINE__);
		mip_solved = 1;
		num_nodes_with_mips_solved++;
	}
	
	if(printing_in_setbranch) printf("optimizing mip\n");
	start_mipsolve = clock();
	
	int num_starts = CPXgetnummipstarts (env_just_solve_mips, global_mip);
	int nzcnt = 0, prev_numsols = 0;
		
    	if(num_starts > global_num_starts)
    	{
    		global_num_starts = num_starts;
    		global_startspace = cur_numcols*global_num_starts;
    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
		global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
		global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
		global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
    	}
	
	status = CPXgetmipstarts (env_just_solve_mips, global_mip, &nzcnt, global_beg, global_varindices, 
			   global_values, global_effortlevel, global_startspace,
			   &surplus, 0, num_starts-1);
			   
	status = CPXaddmipstarts (env_just_solve_mips, nodelp_copy2, num_starts, nzcnt, global_beg, global_varindices,
			   global_values, global_effortlevel, NULL);

	CPXmipopt (env_just_solve_mips, nodelp_copy2);

	int numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, nodelp_copy2);
	int numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, nodelp_copy2);

	num_starts = numsolns - prev_numsols + numrep;

	prev_numsols = numsolns;

    	if(num_starts > global_num_starts)
    	{
    		global_num_starts = num_starts;
    		global_startspace = cur_numcols*global_num_starts;
    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
		global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
		global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
		global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
    	}

	status = CPXgetmipstarts (env_just_solve_mips, nodelp_copy2, &nzcnt, global_beg, global_varindices, 
			   global_values, global_effortlevel, global_startspace,
			   &surplus, 0, num_starts-1);
			   
	status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,
			   global_values, global_effortlevel, NULL);
  	
  	finish_mipsolve = clock();
	double duration_mipsolve = (double)(finish_mipsolve- start_mipsolve) / CLOCKS_PER_SEC;
	if(duration_mipsolve > max_time_to_solve_a_mip) max_time_to_solve_a_mip = duration_mipsolve;
	time_solving_mips += duration_mipsolve;
	if(printing_in_setbranch) printf("time to solve mip at seqnum %d: %lf\n",seqnum,duration_mipsolve);
	
	cumulative_time = (double)(finish_mipsolve - start_BB) / CLOCKS_PER_SEC;
	if(cumulative_time > max_time) goto BRANCHING;
	
	int num_solns = CPXgetsolnpoolnumsolns (env_just_solve_mips, nodelp_copy2);
  	
  	int insert_check = 0;

	if(num_solns >= prev_numsolns) times_to_run = num_solns - prev_numsolns;
  	else times_to_run = num_solns;
  	prev_numsolns = num_solns;
  	
  	for(j=0;j<times_to_run;j++)
  	{
  		status = CPXgetsolnpoolx (env_just_solve_mips, nodelp_copy2, j, x_ws, 0, cur_numcols-1);
	      	insert_check = mock_insert(1,x_ws[obj1_index],x_ws[obj2_index],0,0,0,&tree);
	      	if(insert_check)
	      	{
	      		if(branch_iterations < 5) for(i=0;i<cur_numcols;i++)
	      		{
	      			stored_x[x_rotation][i] = x_ws[i];
	      		}
	      		x_rotation = (x_rotation + 1) % num_x_to_store;
	      		insert_check = 0;
	      		PSA_full(env_global,NULL,x_ws,NULL,NULL);
	      	}
      	}
  	
  	lpstat = CPXgetstat (env_just_solve_mips, nodelp_copy2);
  	if(printing_in_setbranch) printf("solve status: %d\n",lpstat);
  	
  	if(keep_solving_infeasible_MIPs) while(lpstat == 108)
	{
		CPXmipopt (env_just_solve_mips, nodelp_copy2);
	 	lpstat = CPXgetstat (env_just_solve_mips, nodelp_copy2);
	}

	if(lpstat == 103 || lpstat == 119)
	{
		if(printing_in_setbranch) printf("infeasible mip. Fathoming\n");
		fathoming = 1;
		*useraction_p = CPX_CALLBACK_SET;
		nodecnt = 0;
	  	goto TERMINATE;
	}
	else if(lpstat == 108)
	{
		goto BRANCHING;
		printf("Warning: Fathoming a node as infeasible because no feasible solution was found within time limit %lf s. This may cause incorrect solutions\n",time_limit);
		if(printing_in_setbranch) printf("infeasible mip. Fathoming\n");
		fathoming = 1;
		*useraction_p = CPX_CALLBACK_SET;
		nodecnt = 0;
	  	goto TERMINATE;
	}
  	
  	if(printing_in_setbranch) printf("getting x\n");
  	status = CPXgetx (env_just_solve_mips, nodelp_copy2, x_ws, 0, cur_numcols-1);
	if(status) 
	{
		printf ("(%d) CPXgetx, Failed to get x values, error code %d\n", __LINE__,status);
		printf("status of the solve: %d\n",lpstat);
		goto TERMINATE;
	}
	ws_mip_objvals[0] = x_ws[obj1_index];
	ws_mip_objvals[1] = x_ws[obj2_index];
	
	if(!userhandle_up->x_ws) userhandle_up->x_ws = calloc ( (cur_numcols+1),sizeof(double));
	if(!userhandle_down->x_ws) userhandle_down->x_ws = calloc ( (cur_numcols+1),sizeof(double));
	
	for(i=0;i<cur_numcols;i++) 
	{
		userhandle_up->x_ws[i] = x_ws[i];
		userhandle_down->x_ws[i] = x_ws[i];
	}
	userhandle_up->x_ws[cur_numcols] = slope;
	userhandle_down->x_ws[cur_numcols] = slope;
	
	if(x_ws[obj1_index] - reduced_subprob_x_ub >= -.0000001) 
	{
		right_side_dom = 1;
		if(reduced_subprob_y_lb < x_ws[obj2_index])
		{
			reduced_subprob_y_lb = x_ws[obj2_index];
			bound_reduction = 1;
		}
	}
	if(x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001) 
	{
		left_side_dom = 1;
		if(reduced_subprob_y_ub > x_ws[obj2_index])
		{
			reduced_subprob_y_ub = x_ws[obj2_index];
			bound_reduction = 1;
		}
	}

	if(printing_in_setbranch) printf("plotting\n");
/*	if(printing_in_setbranch) */
	if(printing_in_setbranch) printf("plot(%lf,%lf,'mo');\n",x_ws[obj1_index],x_ws[obj2_index]);

	if(lpstat == 101 || lpstat == 102)
	{
		if(printing_in_setbranch) printf("the mip solution was optimal\n");
		ws_mip_opt = 1;
		
		if(control_node_selection)
		{
			userhandle_up->f1 = x_ws[obj1_index];
			userhandle_down->f1 = x_ws[obj1_index];
			userhandle_up->f2 = x_ws[obj2_index];
			userhandle_down->f2 = x_ws[obj2_index];
		}
		
		add_check = mock_insert(1,x_ws[obj1_index],x_ws[obj2_index],0,0,0,&tree);
		if(add_check)
		{
			if(printing_in_setbranch) printf("adding solution\n");
			for(i=0;i<cur_numcols;i++)
	      		{
	      			stored_x[x_rotation][i] = x_ws[i];
	      		}
	      		x_rotation = (x_rotation + 1) % num_x_to_store;
	      		add_check = 0;
	      		PSA_full(env,NULL,x_ws,NULL,NULL);
		}
		else if(objective_space_fathoming)
		{
			int add_check_1 = mock_insert(1,x_ws[obj1_index]+.0001*x_range,x_ws[obj2_index],0,0,0,&tree);
			int add_check_2 = mock_insert(1,x_ws[obj1_index],x_ws[obj2_index]+.0001*y_range,0,0,0,&tree);
			if(!add_check_1 && !add_check_2)
			{
				sub_pr1_x_ub = x_ws[obj1_index];
				sub_pr1_y_lb = x_ws[obj2_index]+.0001*y_range;
			 	sub_pr2_x_lb = x_ws[obj1_index]+.0001*x_range;
				sub_pr2_y_ub = x_ws[obj2_index];
				goto PARETO_BRANCH;
			}
			else if(!add_check_1)
			{
				sub_pr1_x_ub = x_ws[obj1_index];
				sub_pr1_y_lb = x_ws[obj2_index];
			 	sub_pr2_x_lb = x_ws[obj1_index]+.0001*x_range;
				sub_pr2_y_ub = x_ws[obj2_index];
				goto PARETO_BRANCH;
			}
			else if(!add_check_2)
			{
				sub_pr1_x_ub = x_ws[obj1_index];
				sub_pr1_y_lb = x_ws[obj2_index]+.0001*y_range;
			 	sub_pr2_x_lb = x_ws[obj1_index];
				sub_pr2_y_ub = x_ws[obj2_index];
				goto PARETO_BRANCH;
			}
		}
		
		if(its_been_only_points)
		{
			sub_pr1_x_ub = x_ws[obj1_index];
			sub_pr1_y_lb = x_ws[obj2_index];
		 	sub_pr2_x_lb = x_ws[obj2_index] + .001;
			sub_pr2_y_ub = x_ws[obj2_index];
			goto PARETO_BRANCH;
		}
		
		if(x_ws[obj2_index] - reduced_subprob_y_ub >= -.0000001) reduced_subprob_x_lb = x_ws[obj1_index];
		if(x_ws[obj2_index] - reduced_subprob_y_lb <= .0000001) reduced_subprob_x_ub = x_ws[obj1_index];
		
		if(!exact_mips && x_ws[obj1_index] > reduced_subprob_x_lb && x_ws[obj2_index] < reduced_subprob_y_ub && 
			x_ws[obj1_index] < reduced_subprob_x_ub && x_ws[obj2_index] > reduced_subprob_y_lb)
		{
			if(printing_in_setbranch) printf("weighted sum solution in interior of search region\n");
			
			if(!its_been_only_points && !points_only)
			{
				if(going_back) goto BRANCHING;
		
				if(!ob1_lp_been_solved) goto SOLVE_OB1_LP;
				else if(!ob1_lp_int_feas && !ob1_mip_ran) goto SOLVE_OB1_MIP;
				else if(ob1_mip_opt) goto RECHECK_IDEALS;
				else goto BRANCHING;
			}
			
			sub_pr1_x_ub = x_ws[obj1_index];
			sub_pr1_y_lb = x_ws[obj2_index];
		 	sub_pr2_x_lb = x_ws[obj1_index] + .001;
			sub_pr2_y_ub = x_ws[obj2_index];
			add_check = mock_insert(1,x_ws[obj1_index],reduced_subprob_y_ub,0,0,0,&tree);
			if(!add_check)
			{
				if(printing_in_setbranch) printf("left side is dominated\n");
				left_side_dom = 1;
			}
			int add_check2 = mock_insert(1,reduced_subprob_x_ub,x_ws[obj2_index],0,0,0,&tree);
			if(!add_check2)
			{
				if(printing_in_setbranch) printf("right side is dominated\n");
				right_side_dom = 1;
			}
			if(left_side_dom && right_side_dom)
			{
				fathomed_by_dominated_local_ideal_pts++;
				fathoming = 1;
				*useraction_p = CPX_CALLBACK_SET;
				nodecnt = 0;
			  	goto TERMINATE;
			}
			else if(left_side_dom)
			{
				reduced_subprob_x_lb = x_ws[obj1_index];
				reduced_subprob_y_ub = x_ws[obj2_index];
				bound_reduction = 1;
			  	goto SOLVE_OB1_LP;
			}
			else if(right_side_dom)
			{
				reduced_subprob_x_ub = x_ws[obj1_index];
				reduced_subprob_y_lb = x_ws[obj2_index];
				bound_reduction = 1;
			  	goto SOLVE_OB2_LP;
			}
			else if(!ws_lp_dom && pareto_branching) goto PARETO_BRANCH;
		}
		else if(!exact_mips && !ws_lp_dom && pareto_branching && x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
		{
			if(printing_in_setbranch) printf("weighted sum solution on left boundary of search region\n");
/*			DO_IT_AGAIN1:*/
/*			;*/
/*			closest_nodes *two_nodes = find_two_nodes_right_of_val(reduced_subprob_x_lb, reduced_subprob_y_ub, tree);*/
/*			if(!two_nodes) */
/*			{*/
/*				if(points_only || its_been_only_points) goto SOLVE_OB2_MIP;*/
/*				else goto AFTER_THIS3;*/
/*			}*/
/*			if(printing_in_setbranch) printf("the two nodes:\n");*/
/*			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->closest->nw_x,*/
/*				x_ideal-two_nodes->closest->se_x,y_ideal-two_nodes->closest->nw_y,y_ideal-two_nodes->closest->se_y);*/
/*			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->next->nw_x,x_ideal-two_nodes->next->se_x,*/
/*				y_ideal-two_nodes->next->nw_y,y_ideal-two_nodes->next->se_y);*/
/*			if(two_nodes->closest->type == 2 && fabs(two_nodes->closest->nw_x - two_nodes->closest->se_x) < .00000001 && */
/*				fabs(two_nodes->closest->nw_y - two_nodes->closest->se_y) < .00000001) two_nodes->closest->type = 1;*/
/*			if(two_nodes->next->type == 2 && fabs(two_nodes->next->nw_x - two_nodes->next->se_x) < .00000001 && */
/*				fabs(two_nodes->next->nw_y - two_nodes->next->se_y) < .00000001) two_nodes->next->type = 1;*/
/*			*/
/*			if(two_nodes->closest->type == 1 && ((fabs(two_nodes->closest->nw_x - two_nodes->next->nw_x) < .0000001 && */
/*				two_nodes->closest->nw_y - two_nodes->next->nw_y >= -.0000001) || */
/*				(fabs(two_nodes->closest->nw_x - two_nodes->next->se_x) < .0000001 && */
/*				two_nodes->closest->nw_y - two_nodes->next->se_y >= -.0000001)))*/
/*			{*/
/*				if(printing_in_setbranch) */
/*					printf("the closest node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
/*				delete_node(two_nodes->closest);*/
/*				free(two_nodes);*/
/*				goto DO_IT_AGAIN1;*/
/*			}*/
/*			*/
/*			if(x_ideal - two_nodes->closest->nw_x - reduced_subprob_x_ub > .0001 || */
/*				y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_lb < -.0001 || */
/*				x_ideal - two_nodes->closest->nw_x - reduced_subprob_x_lb < -.0001 || */
/*				y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_ub > .0001)*/
/*			{*/
/*				if(printing_in_setbranch) printf("closest node is outside region\n");*/
/*				free(two_nodes);*/
/*				goto SOLVE_OB1_LP;*/
/*			}*/
/*			*/
/*			else if(two_nodes->next->type == 1 && ((fabs(two_nodes->next->nw_x - two_nodes->closest->nw_x) < .0000001 && */
/*				two_nodes->next->nw_y - two_nodes->closest->nw_y >= -.0000001) || */
/*				(fabs(two_nodes->next->nw_x - two_nodes->closest->se_x) < .0000001 && */
/*				two_nodes->next->nw_y - two_nodes->closest->se_y >= -.0000001)))*/
/*			{*/
/*				if(printing_in_setbranch) */
/*					printf("the next node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
/*				delete_node(two_nodes->next);*/
/*				free(two_nodes);*/
/*				goto DO_IT_AGAIN1;*/
/*			}*/
/*			if(x_ideal - two_nodes->next->se_x - reduced_subprob_x_ub <= .0000001)*/
/*			{*/
/*				if(printing_in_setbranch) printf("left pt of next node is also inside the search region.\n");*/
/*				if(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y) < -.0001 && */
/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001 && */
/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_ub <= .0000001)*/
/*				{*/
/*					if(printing_in_setbranch) printf("separation between y-values, split\n");*/
/*					printf("y separation is %lf percent of y_range\n", */
/*							100.*(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y))/y_range);*/
/*					sub_pr1_x_ub = reduced_subprob_x_ub;*/
/*					sub_pr1_y_lb = ((y_ideal - two_nodes->next->se_y) + (y_ideal - two_nodes->closest->nw_y))/2.;*/
/*				 	sub_pr2_x_lb = x_ideal - two_nodes->closest->nw_x + .001;//reduced_subprob_x_lb;*/
/*					sub_pr2_y_ub = ((y_ideal - two_nodes->next->se_y) + (y_ideal - two_nodes->closest->nw_y))/2.;*/
/*					free(two_nodes);*/
/*					goto PARETO_BRANCH;*/
/*				}*/
/*				else if(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x) > .0001)*/
/*				{*/
/*					if(printing_in_setbranch) printf("separation between x-values, split\n");*/
/*					printf("y separation is %lf percent of y_range\n", */
/*							100.*(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x))/x_range);*/
/*					sub_pr1_x_ub = ((x_ideal - two_nodes->next->se_x) + (x_ideal - two_nodes->closest->nw_x))/2.;*/
/*					sub_pr1_y_lb = y_ideal - two_nodes->next->se_y + .001;//reduced_subprob_y_lb;*/
/*				 	sub_pr2_x_lb = ((x_ideal - two_nodes->next->se_x) + (x_ideal - two_nodes->closest->nw_x))/2.;*/
/*					sub_pr2_y_ub = reduced_subprob_y_ub;*/
/*					free(two_nodes);*/
/*					goto PARETO_BRANCH;*/
/*				}*/
/*				else if(x_ideal - two_nodes->next->se_x - reduced_subprob_x_ub < -.001 && */
/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb > .001)*/
/*				{*/
/*					if(printing_in_setbranch) printf("left pt of next node is strictly inside the search region.\n");*/
/*					double ran = (double) rand() / ( (double) RAND_MAX);*/
/*					if(ran < .5)*/
/*					{*/
/*						sub_pr1_x_ub = (x_ideal - two_nodes->closest->nw_x);*/
/*						sub_pr1_y_lb = (y_ideal - two_nodes->closest->nw_y);//reduced_subprob_y_lb;*/
/*					 	sub_pr2_x_lb = (x_ideal - two_nodes->closest->nw_x);*/
/*						sub_pr2_y_ub = reduced_subprob_y_ub;*/
/*						free(two_nodes);*/
/*						goto PARETO_BRANCH;*/
/*					}*/
/*					else*/
/*					{*/
/*						sub_pr1_x_ub = reduced_subprob_x_ub;*/
/*						sub_pr1_y_lb = (y_ideal - two_nodes->closest->nw_y);*/
/*					 	sub_pr2_x_lb = (x_ideal - two_nodes->closest->nw_x);//reduced_subprob_x_lb;*/
/*						sub_pr2_y_ub = (y_ideal - two_nodes->closest->nw_y);*/
/*						free(two_nodes);*/
/*						goto PARETO_BRANCH;*/
/*					}*/
/*				}*/
/*			}*/
/*			AFTER_THIS3:*/
/*			*/
/*			free(two_nodes);*/
			goto SOLVE_OB1_LP;
			
		}
		else if(!exact_mips && !ws_lp_dom && pareto_branching && x_ws[obj1_index] - reduced_subprob_x_ub >= -.0000001)
		{
			if(printing_in_setbranch) printf("weighted sum solution on right boundary of search region\n");
/*			DO_IT_AGAIN2:*/
/*			;*/
/*			closest_nodes *two_nodes = find_two_nodes_left_of_val(reduced_subprob_x_ub, reduced_subprob_y_lb, tree);*/
/*			if(!two_nodes) */
/*			{*/
/*				if(points_only || its_been_only_points) goto SOLVE_OB2_MIP;*/
/*				else goto AFTER_THIS4;*/
/*			}*/
/*			if(printing_in_setbranch) printf("the two nodes:\n");*/
/*			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->closest->nw_x,*/
/*				x_ideal-two_nodes->closest->se_x,y_ideal-two_nodes->closest->nw_y,y_ideal-two_nodes->closest->se_y);*/
/*			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->next->nw_x,x_ideal-two_nodes->next->se_x,*/
/*				y_ideal-two_nodes->next->nw_y,y_ideal-two_nodes->next->se_y);*/
/*			*/
/*			if(two_nodes->closest->type == 2 && fabs(two_nodes->closest->nw_x - two_nodes->closest->se_x) < .00000001 && */
/*				fabs(two_nodes->closest->nw_y - two_nodes->closest->se_y) < .00000001) two_nodes->closest->type = 1;*/
/*			if(two_nodes->next->type == 2 && fabs(two_nodes->next->nw_x - two_nodes->next->se_x) < .00000001 && */
/*				fabs(two_nodes->next->nw_y - two_nodes->next->se_y) < .00000001) two_nodes->next->type = 1;*/
/*			if(two_nodes->closest->type == 1 && ((fabs(two_nodes->closest->nw_x - two_nodes->next->nw_x) < .0000001 && */
/*				two_nodes->closest->nw_y - two_nodes->next->nw_y >= -.0000001) || */
/*				(fabs(two_nodes->closest->nw_x - two_nodes->next->se_x) < .0000001 && */
/*				two_nodes->closest->nw_y - two_nodes->next->se_y >= -.0000001)))*/
/*			{*/
/*				if(printing_in_setbranch) */
/*					printf("the closest node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
/*				delete_node(two_nodes->closest);*/
/*				free(two_nodes);*/
/*				goto DO_IT_AGAIN2;*/
/*			}*/
/*			else if(two_nodes->next->type == 1 && ((fabs(two_nodes->next->nw_x - two_nodes->closest->nw_x) < .0000001 && */
/*				two_nodes->next->nw_y - two_nodes->closest->nw_y >= -.0000001) || */
/*				(fabs(two_nodes->next->nw_x - two_nodes->closest->se_x) < .0000001 && */
/*				two_nodes->next->nw_y - two_nodes->closest->se_y >= -.0000001)))*/
/*			{*/
/*				if(printing_in_setbranch) */
/*					printf("the next node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
/*				delete_node(two_nodes->next);*/
/*				free(two_nodes);*/
/*				goto DO_IT_AGAIN2;*/
/*			}*/
/*			*/
/*			if(x_ideal - two_nodes->closest->se_x - reduced_subprob_x_ub > .0001 || */
/*				y_ideal - two_nodes->closest->se_y - reduced_subprob_y_lb < -.0001 || */
/*				x_ideal - two_nodes->closest->se_x - reduced_subprob_x_lb < -.0001 || */
/*				y_ideal - two_nodes->closest->se_y - reduced_subprob_y_ub > .0001)*/
/*			{*/
/*				if(printing_in_setbranch) printf("closest node is outside region\n");*/
/*				free(two_nodes);*/
/*				goto SOLVE_OB2_LP;*/
/*			}*/
/*			*/
/*			if(x_ideal - two_nodes->next->nw_x - reduced_subprob_x_lb >= -.0000001 && */
/*				y_ideal - two_nodes->next->nw_y - reduced_subprob_y_ub <= .0000001 )*/
/*			{*/
/*				if(printing_in_setbranch) printf("right pt of next node is also inside the search region.\n");*/
/*				if(x_ideal - two_nodes->closest->se_x - (x_ideal - two_nodes->next->nw_x) > .0001)*/
/*				{*/
/*					sub_pr1_x_ub = ((x_ideal - two_nodes->closest->se_x) + (x_ideal - two_nodes->next->nw_x))/2.;*/
/*					sub_pr1_y_lb = y_ideal - two_nodes->closest->se_y;//reduced_subprob_y_lb;*/
/*				 	sub_pr2_x_lb = ((x_ideal - two_nodes->closest->se_x) + (x_ideal - two_nodes->next->nw_x))/2.;*/
/*					sub_pr2_y_ub = reduced_subprob_y_ub;*/
/*					free(two_nodes);*/
/*					goto PARETO_BRANCH;*/
/*				}*/
/*				else if(y_ideal - two_nodes->closest->se_y - (y_ideal - two_nodes->next->nw_y) < -.0001 && */
/*					y_ideal - two_nodes->closest->se_y - reduced_subprob_y_ub <= .0000001 && */
/*					y_ideal - two_nodes->closest->se_y - reduced_subprob_y_lb >= -.0000001)*/
/*				{*/
/*					sub_pr1_x_ub = reduced_subprob_x_ub;*/
/*					sub_pr1_y_lb = ((y_ideal - two_nodes->closest->se_y) + (y_ideal - two_nodes->next->nw_y))/2.;*/
/*				 	sub_pr2_x_lb = x_ideal - two_nodes->next->nw_x;//reduced_subprob_x_lb;*/
/*					sub_pr2_y_ub = ((y_ideal - two_nodes->closest->se_y) + (y_ideal - two_nodes->next->nw_y))/2.;*/
/*					free(two_nodes);*/
/*					goto PARETO_BRANCH;*/
/*				}*/
/*				else if(x_ideal - two_nodes->next->nw_x - reduced_subprob_x_lb > .001 && */
/*					y_ideal - two_nodes->next->nw_y - reduced_subprob_y_ub < -.001)*/
/*				{*/
/*					double ran = (double) rand() / ( (double) RAND_MAX);*/
/*					if(ran < .5)*/
/*					{*/
/*						sub_pr1_x_ub = (x_ideal - two_nodes->next->se_x);*/
/*						sub_pr1_y_lb = (y_ideal - two_nodes->next->se_y);//reduced_subprob_y_lb;*/
/*					 	sub_pr2_x_lb = (x_ideal - two_nodes->next->se_x);*/
/*						sub_pr2_y_ub = reduced_subprob_y_ub;*/
/*						free(two_nodes);*/
/*						goto PARETO_BRANCH;*/
/*					}*/
/*					else*/
/*					{*/
/*						sub_pr1_x_ub = reduced_subprob_x_ub;*/
/*						sub_pr1_y_lb = (y_ideal - two_nodes->next->se_y);*/
/*					 	sub_pr2_x_lb = (x_ideal - two_nodes->next->se_x);//reduced_subprob_x_lb;*/
/*						sub_pr2_y_ub = (y_ideal - two_nodes->next->se_y);*/
/*						free(two_nodes);*/
/*						goto PARETO_BRANCH;*/
/*					}*/
/*				}*/
/*			}*/
/*			AFTER_THIS4:*/
/*			*/
/*			free(two_nodes);*/
			goto SOLVE_OB2_LP;
		}
		
		if(going_back) goto BRANCHING;
		
		if(!ob1_lp_been_solved) goto SOLVE_OB1_LP;
		else if(!ob1_lp_int_feas && !ob1_mip_ran) goto SOLVE_OB1_MIP;
		else if(ob1_mip_opt) goto RECHECK_IDEALS;
		else goto BRANCHING;
	}
	else if(lpstat == 105 || lpstat == 107 || lpstat == 109 || lpstat == 111 || lpstat == 113 || lpstat == 116 || lpstat == 106)
	{
		if(printing_in_setbranch) printf("the mip solution was less than optimal\n");
		ws_mip_opt = 0;
		add_check = mock_insert(1,x_ws[obj1_index],x_ws[obj2_index],0,0,0,&tree);
		if(add_check)
		{
			if(printing_in_setbranch) printf("adding solution\n");
			for(i=0;i<cur_numcols;i++)
	      		{
	      			stored_x[x_rotation][i] = x_ws[i];
	      		}
	      		x_rotation = (x_rotation + 1) % num_x_to_store;
	      		add_check = 0;
	      		PSA_full(env,NULL,x_ws,NULL,NULL);
		}
		
		if(going_back) goto BRANCHING;
		
		status = CPXgetbestobjval (env_just_solve_mips, nodelp_copy2, &best_bound);
		projection = -ws_lp_objvals[1]+best_bound;
		if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",projection,ws_lp_objvals[1]);
/*		printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",-65.,best_bound-5./(-1.),(-1.)*(-best_bound-65),-5.);*/
		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",-65.,best_bound-5./(slope),(slope)*(-best_bound-65),-5.);
	
		add_check = mock_insert(1,projection,ws_lp_objvals[1],0,0,0,&tree);
		if(add_check)
		{
			if(printing_in_setbranch) printf("the projection onto the best objval level curve is not dominated, go to branching\n");
			goto BRANCHING;
		}
		else
		{
			if(printing_in_setbranch) printf("the projection onto the best objval level curve is dominated\n");
			if(printing_in_setbranch) printf("trying this: setting obj vals to be this projected point and checking domination\n");
			x_ws[obj1_index] = projection;
			x_ws[obj2_index] = ws_lp_objvals[1];
			if(!ob1_lp_been_solved) goto SOLVE_OB1_LP;
			else if(!ob1_lp_int_feas) goto SOLVE_OB1_MIP;
			else goto BRANCHING;
		}
	}
	else
	{
		printf("the status of the FAILED mipsolve: %d\n",lpstat);
		exit(0);
	}
	
	/*********** Here we solve and process the single objective LP solution associated with f_1 ************/
	
	SOLVE_OB1_LP:
	
/*	printf("arrived at solve of ob1 lp\n");*/

	if(userhandle_current && userhandle_current->ob1_still_feas) goto SOLVE_OB1_MIP;
	
	if(ws_lp_sol_on_right_bd) 
	{
		ob1_lp_been_solved = 1;
		goto SOLVE_OB2_LP;
	}
	
	lp_ob1 = CPXcloneprob (env, nodelp2, &status);
  	if ( status ) {
    		printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
		goto TERMINATE;
	}
	
	chg_coefs(env, lp_ob1, indexes,-10000000.);
	
  	ob1_lp_been_solved = 1;
  	status = CPXlpopt (env, lp_ob1);
 	if ( status ) {
   		printf ("%s(%d): CPXlpopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
		goto TERMINATE;
	}
	
  	status = CPXgetx (env, lp_ob1, x1, 0, cur_numcols-1);
	if(status) 
	{
		printf ("(%d) CPXgetx, Failed to get x values, error code %d\n", __LINE__,status);
		goto TERMINATE;
	}
	ob1_lp_objvals[0] = x1[obj1_index];
	ob1_lp_objvals[1] = x1[obj2_index];
	
/*	if(x1[obj2_index] > x_ws[obj2_index]) goto SOLVE_WS_MIP;*/
	
	if(ob1_lp_objvals[0] < reduced_subprob_x_ub)
	{
		reduced_subprob_x_ub = ob1_lp_objvals[0];
/*		reduced_subprob_y_lb = ob1_lp_objvals[1];*/
		if(printing_in_setbranch) printf("after changing:\n");
		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
				reduced_subprob_y_lb,reduced_subprob_y_ub);
		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
							reduced_subprob_y_lb,reduced_subprob_y_lb);					
		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
							reduced_subprob_y_ub,reduced_subprob_y_ub);
		bound_reduction = 1;
	}

	if(printing_in_setbranch) printf("plotting ob1 lp soln\n");
	if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",x1[obj1_index],x1[obj2_index]);
	
	all_feas = 1;
	for(i=0;i<total_num_integer;i++)
	{
		k = integer_indices[i];
	  	diff = x1[k] - floor(x1[k]);
/*		printf("diff%d: %lf\n",i,diff);*/
  		if( diff >= .00001 && diff <= .99999)
  		{
  			if(printing_in_setbranch) printf("changing frac index to %d (%d)\n",integer_indices[i],__LINE__);
/*  			printf("plot(%lf,%lf,'o');\n",x1[obj1_index],x1[obj2_index]);*/
	  		frac_index = k;
	  		frac_val = x1[k];
	  		frac_values[k] = frac_val;
	  		if(frac_scores[i] > 0.0001)
	  		{
		  		frac_scores[i] += 1.;
/*		  		num_frac++;*/
		  	}
		  	else
		  	{
		  		frac_scores[i] += multiplier*k;
	  			num_frac++;
		  	}
/*		  	printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
	  		all_feas = 0;
/*	  		printf("num_frac: %d\n",num_frac);*/
/*			printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
		}
	}
	if(all_feas == 1)
	{
		if(printing_in_setbranch) printf("ob1 lp solution is integer feasible\n");
		all_feas = 0;
		
		ob1_mip_opt = 1;
		ob1_mip_objvals[0] = x1[obj1_index];
		ob1_mip_objvals[1] = x1[obj2_index];
		
		if(!userhandle_up->x1) userhandle_up->x1 = calloc ((cur_numcols),sizeof(double));
		if(!userhandle_down->x1) userhandle_down->x1 = calloc ((cur_numcols),sizeof(double));
	
		for(i=0;i<cur_numcols;i++) 
		{
			userhandle_up->x1[i] = x1[i];
			userhandle_down->x1[i] = x1[i];
		}
		
		ob1_lp_int_feas = 1;
	}
		
	add_check = mock_insert(1,x1[obj1_index],x1[obj2_index],0,0,0,&tree);
	if(add_check)
	{
		if(printing_in_setbranch) printf("the ob1 lp solution is not dominated\n");
		if(ob1_lp_int_feas == 1)
		{
			add_check = mock_insert(1,x1[obj1_index],x1[obj2_index],0,0,0,&tree);
			if(add_check)
			{
				if(printing_in_setbranch) printf("adding solution\n");
				for(i=0;i<cur_numcols;i++)
		      		{
			      		stored_x[x_rotation][i] = x1[i];
			      	}
			      	x_rotation = (x_rotation + 1) % num_x_to_store;
		 		add_check = 0;
		     		if(check_for_stopping_PSA_full)
				{
					check_for_stopping_PSA_full = 0;
					PSA_full(env,NULL,x1,NULL,NULL);
					check_for_stopping_PSA_full = 1;
				}
				else PSA_full(env,NULL,x1,NULL,NULL);
			}
			
			add_check = mock_insert(1,x1[obj1_index],x_ws[obj2_index],0,0,0,&tree);
			if(add_check)
			{
				if(printing_in_setbranch) printf("the right partial ideal is not dominated\n");
/*				projection = (-1.)*(x1[obj1_index]-x_ws[obj1_index])+x_ws[obj2_index];*/
				projection = (slope)*(x1[obj1_index]-x_ws[obj1_index])+x_ws[obj2_index];
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],x1[obj1_index],x_ws[obj2_index],projection);
/*				add_check = mock_insert(2,x1[obj1_index],projection,x_ws[obj1_index],x_ws[obj2_index],-1.,&tree);*/
/*				printf("slope used: %lf\n",slope);*/
				add_check = mock_insert(2,x1[obj1_index],projection,x_ws[obj1_index],x_ws[obj2_index],slope,&tree);
				if(add_check)
				{
					endpoint1_x = x1[obj1_index];
					endpoint1_y = projection;
					if(printing_in_setbranch) printf("the right partial ideal segment is not dominated (%d)\n",__LINE__);
					if(ws_mip_ran)
					{
						if(printing_in_setbranch) printf("the ws lp soln was not dominated, so PSA will not produce a dominated dual bound.\n");
						goto BRANCHING;
					}
					else
					{
						if(printing_in_setbranch) printf("ws lp soln was dominated. Consider using PSA here to check for a dominated dual bound.\n");
						status = CPXgetbase (env, lp_ob1, cstat, rstat);
						if(printing_in_setbranch) printf("calling PSA\n");
						int reduce_val = PSA(env, lp_ob1, x1, cstat, rstat, seqnum, 0, lp_ob1);
		  				if(reduce_val == 2)
		  				{
		  					if(printing_in_setbranch) printf("fathoming node for completed PSA\n");
						  	fathomed_by_PSA_completion++;
							fathoming = 1;
							*useraction_p = CPX_CALLBACK_SET;
							nodecnt = 0;
						  	goto TERMINATE;
		  				}
						else 
						{
			  				if(reduce_val == 1)
					  		{
					  			reduced_subprob_x_ub = sub_pr1_x_ub;
				  				reduced_subprob_y_lb = sub_pr1_y_lb;
				  				if(printing_in_setbranch) printf("after changing:\n");
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
									reduced_subprob_x_lb,reduced_subprob_x_lb,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
									reduced_subprob_x_ub,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
									reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
									reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_ub,reduced_subprob_y_ub);
				  				bound_reduction = 1;
					  		}
							if(printing_in_setbranch) printf("PSA did not complete. Try solving MIPs.\n");
							if(!ws_lp_int_feas) goto SOLVE_WS_MIP;
							else goto BRANCHING;
						}
					}
				}
				else
				{
					if(printing_in_setbranch) printf("the right partial ideal segment is dominated. We can begin exploring ob2\n");
					right_pt_dom = 1;
					reduced_subprob_x_ub = x_ws[obj1_index];
					reduced_subprob_y_lb = x_ws[obj2_index];
					if(printing_in_setbranch) printf("after changing:\n");
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_ub,reduced_subprob_y_ub);
					bound_reduction = 1;
					
					if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
					{
						fathomed_by_dominated_local_ideal_segments++;
						if(printing_in_setbranch) printf("fathoming because left side already dominated\n");
						fathoming = 1;
						*useraction_p = CPX_CALLBACK_SET;
						nodecnt = 0;
					  	goto TERMINATE;
					}
					
					if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
					else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
					else goto BRANCHING;
				}
			}
			else
			{
				if(printing_in_setbranch) printf("the right partial ideal is dominated. We can begin exploring ob2\n");
				right_seg_dom = 1;
				reduced_subprob_x_ub = x_ws[obj1_index];
				reduced_subprob_y_lb = x_ws[obj2_index];
				if(printing_in_setbranch) printf("after changing:\n");
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_ub,reduced_subprob_y_ub);
				bound_reduction = 1;
				
				if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
				{
					fathomed_by_dominated_local_ideal_pts++;
					fathoming = 1;
					if(printing_in_setbranch) printf("fathoming because left side already dominated\n");
					*useraction_p = CPX_CALLBACK_SET;
					nodecnt = 0;
				  	goto TERMINATE;
				}
				
				if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
				else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
				else goto BRANCHING;
			}
		}
		else
		{
			
			if(!left_side_dom && ws_lp_dom && remove_dominated_middle)
			{
/*				printf("running this 1\n");*/
				if(!ws_sol_interior) 
				{
					status = CPXgetbase (env, nodelp_copy, cstat_ws, rstat_ws);
					if(status) 
					{
						printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
						goto TERMINATE;
					}
				}
				else 
				{
					CPXgetbase (env, nodelp, cstat_ws, rstat_ws);
					if(status) 
					{
						printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
						goto TERMINATE;
					}
				}
				if(printing_in_setbranch) printf("status: %d\n",status);
		
				sub_pr1_x_ub = x_ws[obj1_index];
				sub_pr1_y_lb = x_ws[obj2_index];
				sub_pr2_x_lb = x_ws[obj1_index];
				sub_pr2_y_ub = x_ws[obj2_index];
			  	PSA_left_check = PSA_reduce_left(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, indexes);
			  	
				if(PSA_left_check == 2)
				{
		  			if(printing_in_setbranch) printf("the right subproblem is now empty by PSA completion\n");
					reduced_subprob_x_ub = x_ws[obj1_index];
		  			reduced_subprob_y_lb = x_ws[obj2_index];
		  			if(printing_in_setbranch) printf("after changing:\n");
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_ub,reduced_subprob_y_ub);
		  			bound_reduction = 1;
		  			right_side_dom = 1;
			  		goto SOLVE_OB2_LP;
			  	}
			  	else if(pareto_branching && sub_pr2_x_lb > sub_pr1_x_ub)
				{
			  		goto PARETO_BRANCH;
			  	}
			}
			
			if(!ob1_mip_ran && !going_back_for_lp1) goto SOLVE_OB1_MIP;
			else if(!ob2_lp_been_solved) 
			{
				if(!ws_lp_sol_on_top_bd) goto SOLVE_OB2_LP;
				else
				{
					ob2_lp_been_solved = 1;
					goto BRANCHING;
				}
			}
			else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
			else goto BRANCHING;
		}
	}
	else
	{
		if(printing_in_setbranch) printf("the ob1 lp solution is dominated\n");
		ob1_lp_dom = 1;
		
		add_check = mock_insert(1,x1[obj1_index],x_ws[obj2_index],0,0,0,&tree);
		if(add_check)
		{
			if(printing_in_setbranch) printf("the right partial ideal is not dominated\n");
/*			projection = (-1.)*(x1[obj1_index]-x_ws[obj1_index])+x_ws[obj2_index];*/
			projection = (slope)*(x1[obj1_index]-x_ws[obj1_index])+x_ws[obj2_index];
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],x1[obj1_index],x_ws[obj2_index],projection);
/*			add_check = mock_insert(2,x1[obj1_index],projection,x_ws[obj1_index],x_ws[obj2_index],-1.,&tree);*/
/*			printf("slope used: %lf\n",slope);*/
			add_check = mock_insert(2,x1[obj1_index],projection,x_ws[obj1_index],x_ws[obj2_index],slope,&tree);
			if(add_check)
			{
				endpoint1_x = x1[obj1_index];
				endpoint1_y = projection;
				if(printing_in_setbranch) printf("the right partial ideal segment is not dominated\n");
				if(printing_in_setbranch) printf("consider solving mip(s) to reduce bounds\n"); 
/*				if(ws_mip_ran) */
/*				{*/
/*					if(printing_in_setbranch) printf("the ws lp soln was not dominated, PSA will not produce dominated dual bound.\n");*/
/*					if(!ob1_mip_ran && !ob1_lp_int_feas && !going_back_for_lp1) goto SOLVE_OB1_MIP;*/
/*					else goto BRANCHING;*/
/*				}*/
/*				else*/
				{
					status = CPXgetbase (env, lp_ob1, cstat, rstat);
					if(ws_lp_int_feas && ob1_lp_int_feas)
					{
						int reduce_val = PSA(env, lp_ob1, x1, cstat, rstat, seqnum, 0, lp_ob1);
		  				if(reduce_val == 2)
		  				{
		  					if(printing_in_setbranch) printf("fathoming node for completed PSA\n");
						  	fathomed_by_PSA_completion++;
							fathoming = 1;
							*useraction_p = CPX_CALLBACK_SET;
							nodecnt = 0;
						  	goto TERMINATE;
		  				}
		  				else if(reduce_val == 1)
				  		{
				  			reduced_subprob_x_ub = sub_pr1_x_ub;
			  				reduced_subprob_y_lb = sub_pr1_y_lb;
			  				if(printing_in_setbranch) printf("after changing:\n");
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
								reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
								reduced_subprob_x_ub,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
								reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
								reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_ub,reduced_subprob_y_ub);
			  				bound_reduction = 1;
				  		}
				  		goto BRANCHING;
					}
					else //if(ws_lp_dom)
					{
						if(printing_in_setbranch) printf("calling PSA reduce right\n");
						int PSA_reduce_right_val = PSA_reduce_right(env, lp_ob1, x1, cstat, rstat, indexes, seqnum);
		  				if(PSA_reduce_right_val == 2)
		  				{
		  					if(printing_in_setbranch) printf("fathoming node for dominated lower bound\n");
						  	fathomed_by_dominated_lb++;
							fathoming = 1;
							*useraction_p = CPX_CALLBACK_SET;
							nodecnt = 0;
						  	goto TERMINATE;
		  				}
							
						if(PSA_reduce_right_val == 1)
		  				{
		  					reduced_subprob_x_ub = sub_pr1_x_ub;
			  				reduced_subprob_y_lb = sub_pr1_y_lb;
			  				if(printing_in_setbranch) printf("after changing:\n");
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
								reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
								reduced_subprob_x_ub,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
								reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
								reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_ub,reduced_subprob_y_ub);
			  				bound_reduction = 1;
			  				if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
		  				}
		  				
						if(printing_in_setbranch) printf("lp dual bound not dominated. Try solving MIPs.\n");
						if(!ws_mip_ran) goto SOLVE_WS_MIP;
						else goto SOLVE_OB1_MIP;
/*						goto BRANCHING;*/
					}
					goto BRANCHING;
				}
			}
			else
			{
				if(printing_in_setbranch) printf("the right partial ideal segment is dominated. Move to consideration of ob2.\n");
				right_pt_dom = 1;
				reduced_subprob_x_ub = x_ws[obj1_index];
				reduced_subprob_y_lb = x_ws[obj2_index];
				if(printing_in_setbranch) printf("after changing:\n");
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_ub,reduced_subprob_y_ub);
				bound_reduction = 1;
				
				if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
				{
					fathomed_by_dominated_local_ideal_segments++;
					fathoming = 1;
					if(printing_in_setbranch) printf("fathoming because left side already dominated\n");
					*useraction_p = CPX_CALLBACK_SET;
					nodecnt = 0;
				  	goto TERMINATE;
				}
				
				if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
				else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
				else goto BRANCHING;
			}
		}
		else
		{
			if(printing_in_setbranch) printf("the right partial ideal is dominated. Move to consideration of ob2.\n");
			right_seg_dom = 1;
			reduced_subprob_x_ub = x_ws[obj1_index];
			reduced_subprob_y_lb = x_ws[obj2_index];
			if(printing_in_setbranch) printf("after changing:\n");
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
					reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_ub,reduced_subprob_y_ub);
			bound_reduction = 1;
			
			if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
			{
				fathomed_by_dominated_local_ideal_pts++;
				fathoming = 1;
				if(printing_in_setbranch) printf("fathoming because left side already dominated\n");
				*useraction_p = CPX_CALLBACK_SET;
				nodecnt = 0;
			  	goto TERMINATE;
			}
			
			if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
			else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
			else goto BRANCHING;
		}
	}
	
	/************** Here we solve the single objective MIP associated with f_1 if we determine that its needed *********/	
	
	SOLVE_OB1_MIP:
	
	if(right_side_dom) goto SOLVE_OB2_LP;
	
	ob1_mip_ran = 1;
/*	if(exact_mips) goto SKIP_THIS5;*/
	if(userhandle_current && userhandle_current->ob1_still_feas && seqnum != 0)
	{
		if(printing_in_setbranch) printf("obj1 soln from parent node is still feasible\n");
		for(j=0;j<cur_numcols;j++) x1[j] = userhandle_current->x1[j];
		ob1_mip_objvals[0] = x1[obj1_index];
		ob1_mip_objvals[1] = x1[obj2_index];
		
		if(!userhandle_up->x1) userhandle_up->x1 = calloc (cur_numcols,sizeof(double));
		if(!userhandle_down->x1) userhandle_down->x1 = calloc (cur_numcols,sizeof(double));
	
		for(i=0;i<cur_numcols;i++) 
		{
			userhandle_up->x1[i] = x1[i];
			userhandle_down->x1[i] = x1[i];
		}
			
		if(printing_in_setbranch) printf("plotting\n");
		if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",x1[obj1_index],x1[obj2_index]);
		
		if(printing_in_setbranch) printf("therefore ob1 mip solution still optimal\n");
		ob1_mip_opt = 1;
		mip_solved = 1;
		
		if(going_back) goto BRANCHING;
		
		add_check = mock_insert(1,x1[obj1_index],x_ws[obj2_index],0,0,0,&tree);
		if(add_check)
		{
			if(printing_in_setbranch) printf("the right partial ideal is not dominated\n");
/*			projection = (-1.)*(x1[obj1_index]-x_ws[obj1_index])+x_ws[obj2_index];*/
			projection = (slope)*(x1[obj1_index]-x_ws[obj1_index])+x_ws[obj2_index];
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],x1[obj1_index],x_ws[obj2_index],projection);
/*			printf("slope used: %lf\n",slope);*/
/*			add_check = mock_insert(2,x1[obj1_index],projection,x_ws[obj1_index],x_ws[obj2_index],-1.,&tree);*/
			add_check = mock_insert(2,x1[obj1_index],projection,x_ws[obj1_index],x_ws[obj2_index],slope,&tree);
			if(add_check)
			{
				if(printing_in_setbranch) printf("the right partial ideal segment is not dominated\n");
				if(ws_mip_ran || ws_lp_int_feas || left_side_dom || right_side_dom) goto BRANCHING;
				else
				{
					if(printing_in_setbranch) printf("since here we haven't actually solved the ws mip, go back and solve ws mip\n");
					goto SOLVE_WS_MIP;
				}
			}
			else
			{
				if(printing_in_setbranch) printf("the right partial ideal segment is dominated. Move to consideration of ob2.\n");
				right_pt_dom = 1;
				reduced_subprob_x_ub = x_ws[obj1_index];
				reduced_subprob_y_lb = x_ws[obj2_index];
				if(printing_in_setbranch) printf("after changing:\n");
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_ub,reduced_subprob_y_ub);
				bound_reduction = 1;
				
				if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
				{
					fathomed_by_dominated_local_ideal_segments++;
					if(printing_in_setbranch) printf("fathoming because left side already dominated\n");
					fathoming = 1;
					*useraction_p = CPX_CALLBACK_SET;
					nodecnt = 0;
				  	goto TERMINATE;
				}
				
				if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
				else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
				else goto BRANCHING;
			}
		}
		else
		{
			if(printing_in_setbranch) printf("the right partial ideal is dominated. Move to consideration of ob2.\n");
			right_seg_dom = 1;
			reduced_subprob_x_ub = x_ws[obj1_index];
			reduced_subprob_y_lb = x_ws[obj2_index];
			if(printing_in_setbranch) printf("after changing:\n");
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
					reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_ub,reduced_subprob_y_ub);
			bound_reduction = 1;
			
			if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
			{
				fathomed_by_dominated_local_ideal_pts++;
				fathoming = 1;
				if(printing_in_setbranch) printf("fathoming because left side already dominated\n");
				*useraction_p = CPX_CALLBACK_SET;
				nodecnt = 0;
			  	goto TERMINATE;
			}
			
			if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
			else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
			else goto BRANCHING;
		}
	}
		
	if( ws_mip_ran && (!userhandle_current || !(userhandle_current->ws_still_feas)))
	{
		if(bound_reduction)
		{
			int ind4[4] = {obj1_index,obj1_index,obj2_index,obj2_index};
			char lu4[4] = {'L','U','L','U'};
			double bds4[4] = {reduced_subprob_x_lb,reduced_subprob_x_ub,reduced_subprob_y_lb,reduced_subprob_y_ub};
			status = CPXchgbds (env_just_solve_mips, nodelp_copy2, 4, ind4, lu4, bds4);
			if ( status ) {
		   		printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
		  	}
		}
		mip_ob1 = CPXcloneprob (env_just_solve_mips, nodelp_copy2, &status);
		if ( status ) {
    			if(printing_in_setbranch) printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
		    	goto TERMINATE;
  		}
  	}
  	else
  	{
  		if(bound_reduction)
		{
			int ind4[4] = {obj1_index,obj1_index,obj2_index,obj2_index};
			char lu4[4] = {'L','U','L','U'};
			double bds4[4] = {reduced_subprob_x_lb,reduced_subprob_x_ub,reduced_subprob_y_lb,reduced_subprob_y_ub};
			status = CPXchgbds (env, nodelp2, 4, ind4, lu4, bds4);
			if ( status ) {
		   		printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
		  	}
		}
		mip_ob1 = CPXcloneprob (env_just_solve_mips, nodelp2, &status);
		if ( status ) {
    			if(printing_in_setbranch) printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
		    	goto TERMINATE;
  		}
  		CPXchgprobtype(env_just_solve_mips, mip_ob1, CPXPROB_MILP);
  		status = CPXchgctype(env_just_solve_mips, mip_ob1, cur_numcols, indexes, xctype);
  	}
  	
  	chg_coefs(env_just_solve_mips, mip_ob1, indexes, -10000000.);
	
	if(!mip_solved)
	{
		if(printing_in_setbranch) printf("here: %d\n",__LINE__);
		mip_solved = 1;
		num_nodes_with_mips_solved++;
	}

	if(printing_in_setbranch) printf("optimizing mip\n");
	start_mipsolve = clock();
	nzcnt = 0;
	prev_numsols = 0;
	
	num_starts = CPXgetnummipstarts (env_just_solve_mips, global_mip);
			    	
    	if(num_starts > global_num_starts)
    	{
    		global_num_starts = num_starts;
    		global_startspace = cur_numcols*global_num_starts;
    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
		global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
		global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
		global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
    	}

	status = CPXgetmipstarts (env_just_solve_mips, global_mip, &nzcnt, global_beg, global_varindices, 
			   global_values, global_effortlevel, global_startspace,
			   &surplus, 0, num_starts-1);
			   
	status = CPXaddmipstarts (env_just_solve_mips, mip_ob1, num_starts, nzcnt, global_beg, global_varindices,
			   global_values, global_effortlevel, NULL);

	CPXmipopt (env_just_solve_mips, mip_ob1);

	numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, mip_ob1);
	numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, mip_ob1);

	num_starts = numsolns - prev_numsols + numrep;

	prev_numsols = numsolns;

    	if(num_starts > global_num_starts)
    	{
    		global_num_starts = num_starts;
    		global_startspace = cur_numcols*global_num_starts;
    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
		global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
		global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
		global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
    	}

	status = CPXgetmipstarts (env_just_solve_mips, mip_ob1, &nzcnt, global_beg, global_varindices, 
			   global_values, global_effortlevel, global_startspace,
			   &surplus, 0, num_starts-1);
			   
	status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,
			   global_values, global_effortlevel, NULL);
			  	
	finish_mipsolve = clock();
	duration_mipsolve = (double)(finish_mipsolve- start_mipsolve) / CLOCKS_PER_SEC;
	if(duration_mipsolve > max_time_to_solve_a_mip) max_time_to_solve_a_mip = duration_mipsolve;
	time_solving_mips += duration_mipsolve;
	if(printing_in_setbranch) printf("time to solve ob1 mip at seqnum %d: %lf\n",seqnum,duration_mipsolve);
	
	cumulative_time = (double)(finish_mipsolve - start_BB) / CLOCKS_PER_SEC;
	if(cumulative_time > max_time) goto BRANCHING;
	
	num_solns = CPXgetsolnpoolnumsolns (env_just_solve_mips, mip_ob1);
  	
  	insert_check = 0;
  	
  	if(num_solns >= prev_numsolns) times_to_run = num_solns - prev_numsolns;
  	else times_to_run = num_solns;
  	prev_numsolns = num_solns;
  	
  	for(j=0;j<times_to_run;j++)
  	{
  		status = CPXgetsolnpoolx (env_just_solve_mips, mip_ob1, j, x1, 0, cur_numcols-1);
	      	insert_check = mock_insert(1,x1[obj1_index],x1[obj2_index],0,0,0,&tree);
	/*      printf("inserting (%lf,%lf)\n",SE_extreme_x,SE_extreme_y);*/
	/*      printf("insert check: %d\n",insert_check);*/
	      	if(insert_check)
	      	{
	/*      	printf("insert was succesful 1\n");*/
	      		//stored_x[x_rotation] = x;
	/*      	printf("x_rotation: %d\n",x_rotation);*/
	      		if(branch_iterations < 5) for(i=0;i<cur_numcols;i++)
	      		{
	      			stored_x[x_rotation][i] = x1[i];
	/*      		printf("stored val: %lf\n",stored_x[x_rotation][i]);*/
	      		}
	      		x_rotation = (x_rotation + 1) % num_x_to_store;
	      		insert_check = 0;
	      		if(check_for_stopping_PSA_full)
			{
				check_for_stopping_PSA_full = 0;
				PSA_full(env,NULL,x1,NULL,NULL);
				check_for_stopping_PSA_full = 1;
			}
			else PSA_full(env,NULL,x1,NULL,NULL);
	      	}
      	}
			  	
  	lpstat = CPXgetstat (env_just_solve_mips, mip_ob1);
	if(printing_in_setbranch) printf("solve status: %d\n",lpstat);
	
	if(keep_solving_infeasible_MIPs) while(lpstat == 108)
	{
		CPXmipopt (env_just_solve_mips, mip_ob1);
	 	lpstat = CPXgetstat (env_just_solve_mips, mip_ob1);
	}
	
	if(lpstat == 103 || lpstat == 119)
	{
		if(printing_in_setbranch) printf("infeasible mip. Fathoming\n");
		fathoming = 1;
		*useraction_p = CPX_CALLBACK_SET;
		nodecnt = 0;
	  	goto TERMINATE;
	}
	else if(lpstat == 108)
	{
		goto BRANCHING;
	}
		  	
  	if(printing_in_setbranch) printf("getting x\n");
	status = CPXgetx (env_just_solve_mips, mip_ob1, x1, 0, cur_numcols-1);
	if(status) 
	{
		printf ("(%d) CPXgetx, Failed to get x values, error code %d\n", __LINE__,status);
    		goto TERMINATE;
	}
	ob1_mip_objvals[0] = x1[obj1_index];
	ob1_mip_objvals[1] = x1[obj2_index];

	if(printing_in_setbranch) printf("plotting\n");
	if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",x1[obj1_index],x1[obj2_index]);
	
	if(!userhandle_up->x1) userhandle_up->x1 = calloc (cur_numcols,sizeof(double));
	if(!userhandle_down->x1) userhandle_down->x1 = calloc (cur_numcols,sizeof(double));
	
	for(i=0;i<cur_numcols;i++) 
	{
		userhandle_up->x1[i] = x1[i];
		userhandle_down->x1[i] = x1[i];
	}
	
	if(lpstat == 101 || lpstat == 102)
	{
		if(printing_in_setbranch) printf("the ob1 mip solution was optimal\n");
		ob1_mip_opt = 1;
		
		if(ob1_mip_objvals[0] < reduced_subprob_x_ub || ob1_mip_objvals[1] > reduced_subprob_y_lb)
		{
			reduced_subprob_x_ub = ob1_mip_objvals[0];
			if(there_will_only_be_points && integer_bb && integer_objective == 1) reduced_subprob_x_ub -= smallest_coef;
			reduced_subprob_y_lb = ob1_mip_objvals[1];
			if(there_will_only_be_points && integer_bb && integer_objective == 2) reduced_subprob_y_lb += smallest_coef;
			if(printing_in_setbranch) printf("after changing:\n");
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
					reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_ub,reduced_subprob_y_ub);
			bound_reduction = 1;
		}
		
		add_check = mock_insert(1,x1[obj1_index],x1[obj2_index],0,0,0,&tree);
		if(add_check)
		{
			if(printing_in_setbranch) printf("adding solution\n");
			for(i=0;i<cur_numcols;i++)
	      		{
		      		stored_x[x_rotation][i] = x1[i];
		      	}
		      	x_rotation = (x_rotation + 1) % num_x_to_store;
	 		add_check = 0;
	     		if(check_for_stopping_PSA_full)
			{
				check_for_stopping_PSA_full = 0;
				PSA_full(env,NULL,x1,NULL,NULL);
				check_for_stopping_PSA_full = 1;
			}
			else PSA_full(env,NULL,x1,NULL,NULL);
		}
		
		if(going_back) goto BRANCHING;
		
		RECHECK_IDEALS:
		
		add_check = 1;
		if(x1[obj1_index] >= x_ws[obj1_index] && x1[obj2_index] <= x_ws[obj2_index]) add_check = mock_insert(1,x1[obj1_index],
														x_ws[obj2_index],0,0,0,&tree);
		if(add_check)
		{
			if(printing_in_setbranch) printf("the right partial ideal is not dominated\n");
			projection = (slope)*(x1[obj1_index]-x_ws[obj1_index])+x_ws[obj2_index];
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],x1[obj1_index],x_ws[obj2_index],projection);
/*			printf("slope used: %lf\n",slope);*/
			if(x1[obj1_index] >= x_ws[obj1_index] && x1[obj2_index] <= x_ws[obj2_index]) add_check = mock_insert(2,
									x1[obj1_index],projection,x_ws[obj1_index],x_ws[obj2_index],slope,&tree);
			if(add_check)
			{
				if(printing_in_setbranch) printf("the right partial ideal segment is not dominated\n");
				if(ws_mip_ran || ws_lp_int_feas || left_side_dom || right_side_dom) goto BRANCHING;
				else
				{
					if(printing_in_setbranch) printf("since here we haven't actually solved the ws mip, go back and solve ws mip\n");
					goto SOLVE_WS_MIP;
				}
			}
			else
			{
				if(printing_in_setbranch) printf("the right partial ideal segment is dominated. Move to consideration of ob2.\n");
				right_seg_dom = 1;
				reduced_subprob_x_ub = x_ws[obj1_index];
				reduced_subprob_y_lb = x_ws[obj2_index];
				if(printing_in_setbranch) printf("after changing:\n");
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_ub,reduced_subprob_y_ub);
				bound_reduction = 1;
				
				if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
				{
					fathomed_by_dominated_local_ideal_segments++;
					fathoming = 1;
					if(printing_in_setbranch) printf("fathoming because left side already dominated\n");
					*useraction_p = CPX_CALLBACK_SET;
					nodecnt = 0;
				  	goto TERMINATE;
				}
				
				if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
				else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
				else goto BRANCHING;
			}
		}
		else
		{
			if(printing_in_setbranch) printf("the right partial ideal is dominated. Move to consideration of ob2.\n");
			right_pt_dom = 1;
			reduced_subprob_x_ub = x_ws[obj1_index];
			reduced_subprob_y_lb = x_ws[obj2_index];
			if(printing_in_setbranch) printf("after changing:\n");
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
					reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_ub,reduced_subprob_y_ub);
			bound_reduction = 1;
			
			if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
			{
				fathomed_by_dominated_local_ideal_pts++;
				fathoming = 1;
				*useraction_p = CPX_CALLBACK_SET;
				nodecnt = 0;
			  	goto TERMINATE;
			}
			
			if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
			else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
			else goto BRANCHING;
		}
	}
	else if(lpstat == 105 || lpstat == 107 || lpstat == 109 || lpstat == 111 || lpstat == 113 || lpstat == 116 || lpstat == 106)
	{
		if(printing_in_setbranch) printf("the mip solution was less than optimal\n");
		ob1_mip_opt = 0;
		add_check = mock_insert(1,x1[obj1_index],x1[obj2_index],0,0,0,&tree);
		if(add_check)
		{
			if(printing_in_setbranch) printf("adding solution\n");
			for(i=0;i<cur_numcols;i++)
			{
	      			stored_x[x_rotation][i] = x_ws[i];
	      		}
			x_rotation = (x_rotation + 1) % num_x_to_store;
	      		add_check = 0;
	      		if(check_for_stopping_PSA_full)
			{
				check_for_stopping_PSA_full = 0;
				PSA_full(env,NULL,x1,NULL,NULL);
				check_for_stopping_PSA_full = 1;
			}
			else PSA_full(env,NULL,x1,NULL,NULL);
		}
		
		if(going_back) goto BRANCHING;
		
		status = CPXgetbestobjval (env_just_solve_mips, mip_ob1, &best_bound);
		double obj_off;
	    	status = CPXgetobjoffset(env_just_solve_mips, mip_ob1, &obj_off );
		best_bound -= obj_off;
		reduced_subprob_x_ub = fmin(reduced_subprob_x_ub, best_bound);
		
		goto BRANCHING;
				
		projection = (-1/1000.)*ob1_lp_objvals[1]+best_bound;
		if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",projection,ob1_lp_objvals[1]);
		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",-65.,best_bound-5./(-1000.),(-1000.)*(-best_bound-65),-5.);

		add_check = mock_insert(1,projection,ob1_lp_objvals[1],0,0,0,&tree);
		if(add_check)
		{
			if(printing_in_setbranch) printf("the projection onto the best objval level curve is not dominated\n");
			goto BRANCHING;
		}
		else
		{
			if(printing_in_setbranch) printf("the projection onto the best ob1 level curve is dominated\n");
			projection = (-1/1000.)*x_ws[obj2_index]+best_bound;
			if(projection <= x_ws[obj1_index])
			{
				printf("projection from ws soln dominated\n");
				exit(0);
			}
			else
			{
				if(printing_in_setbranch) printf("projection from ws soln not dominated\n");
				goto BRANCHING;
			}
		}
	}
	else
	{
		printf("the status of the FAILED mipsolve: %d\n",lpstat);
		exit(0);
	}
	
	/**** Here we solve and process the LP solution associated with objective f_2 ************************/
	
	SOLVE_OB2_LP:
	
/*	if(left_side_dom || (!ws_mip_ran && ws_lp_sol_on_top_bd) ) goto BRANCHING;*/
	
	if(userhandle_current && userhandle_current->ob2_still_feas) goto SOLVE_OB2_MIP;
	
	lp_ob2 = CPXcloneprob (env, nodelp2, &status);
  	if ( status ) {
    		printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
		goto TERMINATE;
  	}
  	chg_coefs(env, lp_ob2, indexes, -.000000001);
  		
  	ob2_lp_been_solved = 1;
  	status = CPXlpopt (env, lp_ob2);
 	if ( status ) {
   		printf ("%s(%d): CPXlpopt, Failed to solve relaxation,  error code %d\n",
   			__FILE__, __LINE__, status);
			goto TERMINATE;
  	}
  	
  	status = CPXgetx (env, lp_ob2, x2, 0, cur_numcols-1);
	if(status) 
	{
		printf ("(%d) CPXgetx, Failed to get x values, error code %d\n", __LINE__,status);
   			goto TERMINATE;
	}
	ob2_lp_objvals[0] = x2[obj1_index];
	ob2_lp_objvals[1] = x2[obj2_index];
	
	if(ob2_lp_objvals[1] < reduced_subprob_y_ub)
	{
/*		reduced_subprob_x_lb = ob2_lp_objvals[0];*/
		reduced_subprob_y_ub = ob2_lp_objvals[1];
		if(printing_in_setbranch) printf("after changing:\n");
		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
				reduced_subprob_y_lb,reduced_subprob_y_ub);
		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
							reduced_subprob_y_lb,reduced_subprob_y_lb);					
		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
							reduced_subprob_y_ub,reduced_subprob_y_ub);
		bound_reduction = 1;
	}

	if(printing_in_setbranch) printf("plotting ob2 lp soln\n");
	if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",x2[obj1_index],x2[obj2_index]);
	
	all_feas = 1;
	for(i=0;i<total_num_integer;i++)
	{
		k = integer_indices[i];
	  	diff = x2[k] - floor(x2[k]);
/*		if(printing_in_setbranch) printf("diff%d: %lf\n",i,diff);*/
  		if( diff >= .00001 && diff <= .99999)
  		{
  			if(printing_in_setbranch) printf("changing frac index to %d (%d)\n",integer_indices[i],__LINE__);
/*  			printf("plot(%lf,%lf,'ko');\n",x2[obj1_index],x2[obj2_index]);*/
		  	frac_index = k;
		  	frac_val = x2[k];
		  	frac_values[k] = frac_val;
	 		if(frac_scores[i] > .0001)
	  		{
		  		frac_scores[i] += 1.;
/*		  		num_frac++;*/
		  	}
		  	else
		  	{
		  		frac_scores[i] += multiplier*k;
	  			num_frac++;
		  	}
/*		  	printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
		  	all_feas = 0;
/*		  	printf("num_frac: %d\n",num_frac);*/
/*			printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
		}
	}
	if(all_feas == 1)
	{
		if(printing_in_setbranch) printf("ob2 lp solution is integer feasible\n");
		all_feas = 0;
		
		ob2_mip_opt = 1;
		ob2_mip_objvals[0] = x2[obj1_index];
		ob2_mip_objvals[1] = x2[obj2_index];
		
		if(!userhandle_up->x2) userhandle_up->x2 = calloc ((cur_numcols),sizeof(double));
		if(!userhandle_down->x2) userhandle_down->x2 = calloc ((cur_numcols),sizeof(double));
	
		for(i=0;i<cur_numcols;i++) 
		{
			userhandle_up->x2[i] = x2[i];
			userhandle_down->x2[i] = x2[i];
		}
		
		ob2_lp_int_feas = 1;
	}

	add_check = mock_insert(1,x2[obj1_index],x2[obj2_index],0,0,0,&tree);
	if(add_check)
	{
		if(printing_in_setbranch) printf("the ob2 lp solution is not dominated\n");
		if(ob2_lp_int_feas == 1)
		{
			add_check = mock_insert(1,x2[obj1_index],x2[obj2_index],0,0,0,&tree);
			if(add_check)
			{
				if(printing_in_setbranch) printf("adding solution\n");
				for(i=0;i<cur_numcols;i++)
		      		{
			      		stored_x[x_rotation][i] = x2[i];
			      	}
			      	x_rotation = (x_rotation + 1) % num_x_to_store;
		 		add_check = 0;
		     		if(check_for_stopping_PSA_full)
				{
					check_for_stopping_PSA_full = 0;
					PSA_full(env,NULL,x2,NULL,NULL);
					check_for_stopping_PSA_full = 1;
				}
				else PSA_full(env,NULL,x2,NULL,NULL);
			}
			
			add_check = mock_insert(1,x_ws[obj1_index],x2[obj2_index],0,0,0,&tree);
			if(add_check)
			{
				if(printing_in_setbranch) printf("the left partial ideal is not dominated\n");
				
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],
								x2[obj2_index]);
				projection = (1./slope)*(x2[obj2_index]-x_ws[obj2_index])+x_ws[obj1_index];
				add_check = mock_insert(2,x_ws[obj1_index],x_ws[obj2_index],projection,x2[obj2_index],slope,&tree);
				if(add_check)
				{
					endpoint2_x = projection;
					endpoint2_y = x2[obj2_index];
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],
									x2[obj2_index]);
					if(printing_in_setbranch) printf("the left partial ideal segment is not dominated\n");
					
					status = CPXgetbase (env, lp_ob2, cstat, rstat);
					if(printing_in_setbranch) printf("calling PSA left\n");
					int reduce_val = PSA_left(env, lp_ob2, x2, cstat, rstat, seqnum, 0, lp_ob2);
	  				if(reduce_val == 2)
	  				{
	  					if(printing_in_setbranch) printf("fathoming node for completed PSA\n");
					  	fathomed_by_PSA_completion++;
						fathoming = 1;
						*useraction_p = CPX_CALLBACK_SET;
						nodecnt = 0;
					  	goto TERMINATE;
	  				}
	  				else if(reduce_val == 1)
	  				{
	  					reduced_subprob_x_lb = sub_pr2_x_lb;
		  				reduced_subprob_y_ub = sub_pr2_y_ub;
		  				if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
		  				bound_reduction = 1;
	  				}
					if(!ws_mip_ran && !ws_lp_int_feas && !right_side_dom) goto SOLVE_WS_MIP;
					else goto BRANCHING;
				}
				else
				{
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],
									x2[obj2_index]);
					if(printing_in_setbranch) printf("Left partial ideal segment is dominated. Fathom (%d)\n",__LINE__);
					if(right_pt_dom || right_side_dom) fathomed_by_dominated_local_ideal_pts++;
					else if(right_seg_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
					else
					{
						reduced_subprob_x_lb = x_ws[obj1_index];
		  				reduced_subprob_y_ub = x_ws[obj2_index];
		  				if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
		  				bound_reduction = 1;
		  				goto BRANCHING;
					}
					fathoming = 1;
					*useraction_p = CPX_CALLBACK_SET;
				  	nodecnt = 0;
					goto TERMINATE;
				}
			}
			else
			{
				left_side_dom = 1;
				if(printing_in_setbranch) printf("Left partial ideal is dominated. Fathom\n");
				if(right_pt_dom || right_side_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
				else if(right_seg_dom) fathomed_by_dominated_local_ideal_segments++;
				else
				{
					reduced_subprob_x_lb = x_ws[obj1_index];
	  				reduced_subprob_y_ub = x_ws[obj2_index];
	  				if(printing_in_setbranch) printf("after changing:\n");
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_ub,reduced_subprob_y_ub);
	  				bound_reduction = 1;
	  				goto BRANCHING;
				}
				fathoming = 1;
				*useraction_p = CPX_CALLBACK_SET;
			  	nodecnt = 0;
				goto TERMINATE;
			}
		}
		else
		{
			if(printing_in_setbranch) printf("ob2 lp solution is not integer feasible\n");
			
			if(!right_side_dom && ws_lp_dom && !ob1_lp_dom && remove_dominated_middle)
			{
/*				printf("running this 2\n");*/
				if(!ws_sol_interior) 
				{
					status = CPXgetbase (env, nodelp_copy, cstat_ws, rstat_ws);
					if(status) 
					{
						printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
						goto SOLVE_OB2_MIP;
					}
				}
				else 
				{
					CPXgetbase (env, nodelp, cstat_ws, rstat_ws);
					if(status) 
					{
						printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
						goto SOLVE_OB2_MIP;
					}
				}
				if(printing_in_setbranch) printf("status: %d\n",status);
		
				sub_pr1_x_ub = x_ws[obj1_index];
				sub_pr1_y_lb = x_ws[obj2_index];
				sub_pr2_x_lb = x_ws[obj1_index];
				sub_pr2_y_ub = x_ws[obj2_index];
			  	PSA_right_check = PSA_reduce_right(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, indexes, seqnum);
			  	PSA_left_check = PSA_reduce_left(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, indexes);
			  	
				if(PSA_right_check == 2 && PSA_left_check == 2)
		  		{
			  		if(printing_in_setbranch) printf("fathoming node %d for PSA completion (%d)\n",seqnum,__LINE__);
				  	fathomed_by_dominated_lb++;
					fathoming = 1;
					*useraction_p = CPX_CALLBACK_SET;
					nodecnt = 0;
			  		goto TERMINATE;
		  		}
			  	else if(PSA_right_check == 2)
			  	{
					if(printing_in_setbranch) printf("the left subproblem is now empty by PSA completion\n");
			  		reduced_subprob_x_lb = x_ws[obj1_index];
			  		reduced_subprob_y_ub = x_ws[obj2_index];
			  		reduced_subprob_x_ub = ubs[0];
					reduced_subprob_y_lb = lbs[1];
					if(printing_in_setbranch) printf("after changing:\n");
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_ub,reduced_subprob_y_ub);
					bound_reduction = 1;
					left_side_dom = 1;
			  		goto SOLVE_OB1_LP;
			  	}
			  	else if(PSA_left_check == 2)
				{
		  			if(printing_in_setbranch) printf("the right subproblem is now empty by PSA completion\n");
			  		reduced_subprob_x_lb = lbs[0];
			  		reduced_subprob_y_ub = ubs[1];
					reduced_subprob_x_ub = x_ws[obj1_index];
		  			reduced_subprob_y_lb = x_ws[obj2_index];
		  			if(printing_in_setbranch) printf("after changing:\n");
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_ub,reduced_subprob_y_ub);
		  			bound_reduction = 1;
		  			right_side_dom = 1;
			  		goto SOLVE_OB2_LP;
			  	}
			  	else if(pareto_branching && sub_pr2_x_lb > sub_pr1_x_ub)
				{
			  		goto PARETO_BRANCH;
			  	}
			  	else if(!ob1_lp_been_solved)
			  	{
			  		if(printing_in_setbranch) printf("none of the above\n");
			  		goto SOLVE_OB1_LP;
			  	}
			  }
			  else if(!right_side_dom && ws_lp_dom && remove_dominated_middle)
			  {
/*			  	printf("running this 3\n");*/
				if(!ws_sol_interior) 
				{
					status = CPXgetbase (env, nodelp_copy, cstat_ws, rstat_ws);
					if(status) 
					{
						printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
						goto SOLVE_OB2_MIP;
					}
				}
				else 
				{
					CPXgetbase (env, nodelp, cstat_ws, rstat_ws);
					if(status) 
					{
						printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
						goto SOLVE_OB2_MIP;
					}
				}
				if(printing_in_setbranch) printf("status: %d\n",status);
		
				sub_pr1_x_ub = x_ws[obj1_index];
				sub_pr1_y_lb = x_ws[obj2_index];
				sub_pr2_x_lb = x_ws[obj1_index];
				sub_pr2_y_ub = x_ws[obj2_index];
			  	PSA_right_check = PSA_reduce_right(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, indexes, seqnum);
			  	
			  	if(PSA_right_check == 2)
			  	{
					if(printing_in_setbranch) printf("the left subproblem is now empty by PSA completion\n");
			  		reduced_subprob_x_lb = x_ws[obj1_index];
			  		reduced_subprob_y_ub = x_ws[obj2_index];
					if(printing_in_setbranch) printf("after changing:\n");
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_ub,reduced_subprob_y_ub);
					bound_reduction = 1;
					left_side_dom = 1;
			  		goto SOLVE_OB1_LP;
			  	}
			  	else if(pareto_branching && sub_pr2_x_lb > sub_pr1_x_ub)
				{
			  		goto PARETO_BRANCH;
			  	}
			}
			
			if((ws_mip_ran || ob1_mip_ran) && !ob2_mip_ran && !going_back_for_lp2) goto SOLVE_OB2_MIP;
			else goto BRANCHING;
		}
	}
	else
	{
		if(printing_in_setbranch) printf("the ob2 lp solution is dominated\n");
		ob2_lp_dom = 1;
		add_check = mock_insert(1,x_ws[obj1_index], x2[obj2_index],0,0,0,&tree);
		if(add_check)
		{
			if(printing_in_setbranch) printf("the left partial ideal is not dominated\n");
			projection = (1./slope)*(x2[obj2_index]-x_ws[obj2_index])+x_ws[obj1_index];
			add_check = mock_insert(2,x_ws[obj1_index],x_ws[obj2_index],projection,x2[obj2_index],slope,&tree);
			if(add_check)
			{
				endpoint2_x = projection;
				endpoint2_y = x2[obj2_index];
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],x2[obj2_index]);
				if(printing_in_setbranch) printf("the left partial ideal segment is not dominated\n");
				
				status = CPXgetbase (env, lp_ob2, cstat, rstat);
				
				if(ob2_lp_int_feas && ob1_lp_int_feas && ws_lp_int_feas)
				{
					if(printing_in_setbranch) printf("calling PSA left\n");
					int reduce_val = PSA_left(env, lp_ob2, x2, cstat, rstat, seqnum, 0, lp_ob2);
	  				if(reduce_val == 2)
	  				{
	  					if(printing_in_setbranch) printf("fathoming node for completed PSA\n");
					  	fathomed_by_PSA_completion++;
						fathoming = 1;
						*useraction_p = CPX_CALLBACK_SET;
						nodecnt = 0;
					  	goto TERMINATE;
	  				}
	  				else if(reduce_val == 1)
	  				{
	  					reduced_subprob_x_lb = sub_pr2_x_lb;
		  				reduced_subprob_y_ub = sub_pr2_y_ub;
		  				if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
		  				bound_reduction = 1;
	  				}
				}
				else //if(ws_lp_dom && (!ob1_lp_been_solved || ob1_lp_dom))
				{
					if(printing_in_setbranch) printf("calling PSA reduce left\n");
					int PSA_reduce_left_val = PSA_reduce_left(env, lp_ob2, x2, cstat, rstat, indexes);
	  				if(PSA_reduce_left_val == 2)
	  				{
	  					if(printing_in_setbranch) printf("fathoming node for dominated lower bound\n");
						fathomed_by_dominated_lb++;
						fathoming = 1;
						*useraction_p = CPX_CALLBACK_SET;
						nodecnt = 0;
					  	goto TERMINATE;
	  				}
	  				else if(PSA_reduce_left_val == 1)
	  				{
	  					reduced_subprob_x_lb = sub_pr2_x_lb;
		  				reduced_subprob_y_ub = sub_pr2_y_ub;
		  				if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
		  				bound_reduction = 1;
		  				bound_reduced_from_PSA_reduce_left = 1;
		  				goto BRANCHING;
	  				}
  				}
	  				
				if((!ob2_lp_int_feas && (right_side_dom || right_seg_dom || right_pt_dom))
					 || (ws_mip_ran && ob1_mip_ran && !ob2_mip_ran && !ob2_lp_int_feas && !going_back_for_lp2)) goto SOLVE_OB2_MIP;
				else if(!right_side_dom && !right_seg_dom && !right_pt_dom && ws_mip_ran && !ob1_mip_ran && !ob1_lp_int_feas) goto SOLVE_OB1_MIP;
				else if(!ws_mip_ran && !ws_lp_int_feas) goto SOLVE_WS_MIP;
				else goto BRANCHING;
			}
			else
			{
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],x2[obj2_index]);
				if(printing_in_setbranch) printf("Left partial ideal segment is dominated. Fathom (%d)\n",__LINE__);
				if(right_pt_dom || right_side_dom) fathomed_by_dominated_local_ideal_pts++;
				else if(right_seg_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
				else
				{
					reduced_subprob_x_lb = x_ws[obj1_index];
	  				reduced_subprob_y_ub = x_ws[obj2_index];
	  				if(printing_in_setbranch) printf("after changing:\n");
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_ub,reduced_subprob_y_ub);
	  				bound_reduction = 1;
	  				goto BRANCHING;
				}
				fathoming = 1;
				*useraction_p = CPX_CALLBACK_SET;
			  	nodecnt = 0;
				goto TERMINATE;
			}
		}
		else
		{
			left_side_dom = 1;
			if(printing_in_setbranch) printf("the left partial ideal is dominated. Fathom\n");
			if(right_pt_dom || right_side_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
			else if(right_seg_dom) fathomed_by_dominated_local_ideal_segments++;
			else
			{
				reduced_subprob_x_lb = x_ws[obj1_index];
  				reduced_subprob_y_ub = x_ws[obj2_index];
  				if(printing_in_setbranch) printf("after changing:\n");
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_ub,reduced_subprob_y_ub);
  				bound_reduction = 1;
  				goto BRANCHING;
			}
			fathoming = 1;
			*useraction_p = CPX_CALLBACK_SET;
		  	nodecnt = 0;
		  	goto TERMINATE;
		}
	}
	
	/************* Here we solve the single objective MIP associated with objective f_2 if we determine its needed *************/
	
	SOLVE_OB2_MIP:
			
	ob2_mip_ran = 1;

/*	if(exact_mips) goto SKIP_THIS7;*/
	if(userhandle_current && userhandle_current->ob2_still_feas)
	{
		if(printing_in_setbranch) printf("obj2 soln from parent node is still feasible\n");
		for(j=0;j<cur_numcols;j++) x2[j] = userhandle_current->x2[j];
        	ob2_mip_objvals[0] = x2[obj1_index];
		ob2_mip_objvals[1] = x2[obj2_index];
		
		if(!userhandle_up->x2) userhandle_up->x2 = calloc (cur_numcols,sizeof(double));
		if(!userhandle_down->x2) userhandle_down->x2 = calloc (cur_numcols,sizeof(double));
	
		for(i=0;i<cur_numcols;i++) 
		{
			userhandle_up->x2[i] = x2[i];
			userhandle_down->x2[i] = x2[i];
		}
				
		ob2_mip_opt = 1;
		mip_solved = 1;
		
		if(going_back) goto BRANCHING;
		
		add_check = mock_insert(1,x_ws[obj1_index],x2[obj2_index],0,0,0,&tree);
		if(add_check)
		{
			if(printing_in_setbranch) printf("the left partial ideal is not dominated\n");
			
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],x2[obj2_index]);
			projection = (1./slope)*(x2[obj2_index]-x_ws[obj2_index])+x_ws[obj1_index];
			add_check = mock_insert(2,x_ws[obj1_index],x_ws[obj2_index],projection,x2[obj2_index],slope,&tree);
			if(add_check)
			{
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],
								x2[obj2_index]);
				if(printing_in_setbranch) printf("the left partial ideal segment is not dominated\n");
				goto BRANCHING;
			}
			else
			{
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],
								x2[obj2_index]);
				if(printing_in_setbranch) printf("Left partial ideal segment is dominated. Fathom (%d)\n",__LINE__);
				if(right_pt_dom) fathomed_by_dominated_local_ideal_pts++;
				else if(right_seg_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
				else 
				{
					if(printing_in_setbranch) printf("somehow fathoming even though neither right pt nor seg has been labelled as dominated\n");
					
						reduced_subprob_x_lb = x_ws[obj1_index];
		  				reduced_subprob_y_ub = x_ws[obj2_index];
		  				if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
		  				bound_reduction = 1;
		  				goto BRANCHING;
					
				}
				fathoming = 1;
				*useraction_p = CPX_CALLBACK_SET;
			  	nodecnt = 0;
				goto TERMINATE;
			}
		}
		else
		{
			left_side_dom = 1;
			if(printing_in_setbranch) printf("the left partial ideal is dominated. Fathom\n");
			if(right_pt_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
			else if(right_seg_dom) fathomed_by_dominated_local_ideal_segments++;
			else 
			{
				if(printing_in_setbranch) printf("somehow fathoming even though neither right pt nor seg has been labelled as dominated\n");
				
						reduced_subprob_x_lb = x_ws[obj1_index];
		  				reduced_subprob_y_ub = x_ws[obj2_index];
		  				if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
		  				bound_reduction = 1;
		  				goto BRANCHING;
					
			}
			fathoming = 1;
			*useraction_p = CPX_CALLBACK_SET;
		  	nodecnt = 0;
		  	goto TERMINATE;
		}
	}
		
	if(!mip_ob2)
	{
		if(ob1_mip_ran && (!userhandle_current || !(userhandle_current->ob1_still_feas)))
		{
			if(bound_reduction)
			{
				int ind4[4] = {obj1_index,obj1_index,obj2_index,obj2_index};
				char lu4[4] = {'L','U','L','U'};
				double bds4[4] = {reduced_subprob_x_lb,reduced_subprob_x_ub,reduced_subprob_y_lb,reduced_subprob_y_ub};
				status = CPXchgbds (env_just_solve_mips, mip_ob1, 4, ind4, lu4, bds4);
				if ( status ) {
			   		printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
			  	}
			}
			mip_ob2 = CPXcloneprob (env_just_solve_mips, mip_ob1, &status);
		  	if ( status ) {
		    		printf ("(%d) CPXcloneprob, Failed to clone problem, error code %d\n", __LINE__, status);
				goto TERMINATE;
			}
		}
		else if(ws_mip_ran  && (!userhandle_current || !(userhandle_current->ws_still_feas)))
		{
			if(bound_reduction)
			{
				int ind4[4] = {obj1_index,obj1_index,obj2_index,obj2_index};
				char lu4[4] = {'L','U','L','U'};
				double bds4[4] = {reduced_subprob_x_lb,reduced_subprob_x_ub,reduced_subprob_y_lb,reduced_subprob_y_ub};
				status = CPXchgbds (env_just_solve_mips, nodelp_copy2, 4, ind4, lu4, bds4);
				if ( status ) {
			   		printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
			  	}
			}
			mip_ob2 = CPXcloneprob (env_just_solve_mips, nodelp_copy2, &status);
		  	if ( status ) {
		    		printf ("(%d) CPXcloneprob, Failed to clone problem, error code %d\n", __LINE__, status);
				goto TERMINATE;
			}
		}
		else
		{
			WHERE_YOU_SHOULD:
			if(bound_reduction)
			{
				int ind4[4] = {obj1_index,obj1_index,obj2_index,obj2_index};
				char lu4[4] = {'L','U','L','U'};
				double bds4[4] = {reduced_subprob_x_lb,reduced_subprob_x_ub,reduced_subprob_y_lb,reduced_subprob_y_ub};
				status = CPXchgbds (env, nodelp2, 4, ind4, lu4, bds4);
				if ( status ) {
			   		printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
			  	}
			}
			mip_ob2 = CPXcloneprob (env_just_solve_mips, nodelp2, &status);
			if ( status ) {
	    			printf ("(%d) CPXcloneprob, Failed to clone problem, error code %d\n", __LINE__, status);
			    	goto TERMINATE;
	  		}
	  		CPXchgprobtype(env_just_solve_mips, mip_ob2, CPXPROB_MILP);
	  		status = CPXchgctype(env_just_solve_mips, mip_ob2, cur_numcols, indexes, xctype);
		}
	}
	
	if(bound_reduced_from_PSA_reduce_left)
	{
		int chg_indices[2] = {obj1_index,obj2_index};
		char chg_lu[2] = {'L','U'};
		double chg_bds[2] = {sub_pr2_x_lb,sub_pr2_y_ub};
		status = CPXchgbds (env_just_solve_mips, mip_ob2, 2, chg_indices, chg_lu, chg_bds);
		if ( status ){
		 	printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
		 	goto TERMINATE;}
	}
	
	chg_coefs(env_just_solve_mips, mip_ob2, indexes, -.000000001);
	
	if(!mip_solved) 
	{
		if(printing_in_setbranch) printf("here: %d\n",__LINE__);
		mip_solved = 1;
		num_nodes_with_mips_solved++;
	}

	if(printing_in_setbranch) printf("optimizing mip\n");
	
	RESOLVE_MIP:
	
	start_mipsolve = clock();
	nzcnt = 0;
	prev_numsols = 0;
	
	num_starts = CPXgetnummipstarts (env_just_solve_mips, global_mip);
	  	
    	if(num_starts > global_num_starts)
    	{
    		global_num_starts = num_starts;
    		global_startspace = cur_numcols*global_num_starts;
    		int *temp = (int *) realloc (global_beg, (global_num_starts*2)*sizeof(int));
    		if(temp != NULL) global_beg = temp;
    		else 
    		{
    			if(global_beg) free(global_beg);
    			global_beg = NULL;
    			global_beg = (int *) malloc ((global_num_starts*2)*sizeof(int));
    			printf("failed on assignment to global beg\n");
    		}
		int *temp2 = (int *) realloc (global_varindices, (global_startspace*2)*sizeof(int));
		if(temp2 != NULL) global_varindices = temp2;
    		else 
    		{
    			if(global_varindices) free(global_varindices);
    			global_varindices = NULL;
    			global_varindices = (int *) malloc ((global_startspace*2)*sizeof(int));
    			printf("failed on assignment to global varind\n");
    		}
		double *temp3 = (double *) realloc (global_values, (global_startspace*2)*sizeof(double));
		if(temp3 != NULL) global_values = temp3;
    		else 
    		{
    			if(global_values) free(global_values);
    			global_values = NULL;
    			global_values = (double *) malloc ((global_startspace*2)*sizeof(double));
    			printf("failed on assignment to global val\n");
    		}
		int *temp4 = (int *) realloc (global_effortlevel, (global_startspace*2)*sizeof(int));
		if(temp4 != NULL) global_effortlevel = temp4;
    		else 
    		{	
    			if(global_effortlevel) free(global_effortlevel);
    			global_effortlevel = NULL;
    			global_effortlevel = (int *) malloc ((global_startspace*2)*sizeof(int));
    			printf("failed on assignment to global eff\n");
    		}
    	}
    	
	int mip_starts_should_be_added = 1;
	if(num_starts)
	{
		status = CPXgetmipstarts (env_just_solve_mips, global_mip, &nzcnt, global_beg, global_varindices, 
				   global_values, global_effortlevel, global_startspace,
				   &surplus, 0, num_starts-1);
		if(seqnum == 14) for(i=0;i<nzcnt;i++) 
		{
			if(global_varindices[i] > cur_numcols || global_varindices[i] < 0)
			{
				printf("for some reason a mip start has become corrupted. Delete it!\n");
				int k = 0;
				while(i > global_beg[k]) k++;
				printf("value of k: %d\n",k);
				status = CPXdelmipstarts (env_just_solve_mips, global_mip, k-1,k-1);
				status = CPXgetmipstarts (env_just_solve_mips, global_mip, &nzcnt, global_beg, global_varindices, 
				   	global_values, global_effortlevel, global_startspace,
				   	&surplus, 0, num_starts-2);
				break;
			}
		}
				   
		if(mip_starts_should_be_added) status = CPXaddmipstarts (env_just_solve_mips, mip_ob2, num_starts, nzcnt, global_beg, global_varindices,
				   global_values, NULL, NULL);
	}
	
	CPXmipopt (env_just_solve_mips, mip_ob2);

	numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, mip_ob2);
	numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, mip_ob2);

	num_starts = numsolns - prev_numsols + numrep;

	prev_numsols = numsolns;

    	if(num_starts > global_num_starts)
    	{
    		global_num_starts = num_starts;
    		global_startspace = cur_numcols*global_num_starts;
    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
		global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
		global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
		global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
    	}

	status = CPXgetmipstarts (env_just_solve_mips, mip_ob2, &nzcnt, global_beg, global_varindices, 
			   global_values, global_effortlevel, global_startspace,
			   &surplus, 0, num_starts-1);
			   
	status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,
			   global_values, global_effortlevel, NULL);
  	
  	finish_mipsolve = clock();
	duration_mipsolve = 
		(double)(finish_mipsolve- start_mipsolve) / CLOCKS_PER_SEC;
	if(duration_mipsolve > max_time_to_solve_a_mip) max_time_to_solve_a_mip = duration_mipsolve;
	time_solving_mips += duration_mipsolve;
	if(printing_in_setbranch) printf("time to solve ob2 mip at seqnum %d: %lf\n",seqnum,duration_mipsolve);
	
	cumulative_time = (double)(finish_mipsolve - start_BB) / CLOCKS_PER_SEC;
	if(cumulative_time > max_time) goto BRANCHING;
	
	num_solns = CPXgetsolnpoolnumsolns (env_just_solve_mips, mip_ob2);
  	
  	insert_check = 0;
  	
  	if(num_solns >= prev_numsolns) times_to_run = num_solns - prev_numsolns;
  	else times_to_run = num_solns;
  	prev_numsolns = num_solns;
  	
  	for(j=0;j<times_to_run;j++)
  	{
  		status = CPXgetsolnpoolx (env_just_solve_mips, mip_ob2, j, x2, 0, cur_numcols-1);
	      	insert_check = mock_insert(1,x2[obj1_index],x2[obj2_index],0,0,0,&tree);
	      	if(insert_check)
	      	{
	      		if(branch_iterations < 5) for(i=0;i<cur_numcols;i++)
	      		{
	      			stored_x[x_rotation][i] = x2[i];
	      		}
	      		x_rotation = (x_rotation + 1) % num_x_to_store;
	      		insert_check = 0;
	      		if(check_for_stopping_PSA_full)
			{
				check_for_stopping_PSA_full = 0;
				PSA_full(env,NULL,x2,NULL,NULL);
				check_for_stopping_PSA_full = 1;
			}
			else PSA_full(env,NULL,x2,NULL,NULL);
	      	}
      	}
  	
  	lpstat = CPXgetstat (env_just_solve_mips, mip_ob2);
	if(printing_in_setbranch) printf("solve status: %d\n",lpstat);
	
	if(keep_solving_infeasible_MIPs) while(lpstat == 108)
	{
		CPXmipopt (env_just_solve_mips, mip_ob2);
	 	lpstat = CPXgetstat (env_just_solve_mips, mip_ob2);
	}

	if(lpstat == 103 || lpstat == 119 || lpstat == 115)
	{
		if(printing_in_setbranch) printf("infeasible mip. Fathoming\n");
		fathoming = 1;
		*useraction_p = CPX_CALLBACK_SET;
		nodecnt = 0;
	  	goto TERMINATE;
	}
	else if(lpstat == 108)
	{
		goto BRANCHING;
		printf("Warning: Fathoming a node as infeasible because no feasible solution was found within time limit %lf s. This may cause incorrect solutions\n",time_limit);
		if(printing_in_setbranch) printf("infeasible mip. Fathoming\n");
		fathoming = 1;
		*useraction_p = CPX_CALLBACK_SET;
		nodecnt = 0;
	  	goto TERMINATE;
	}
	  	
  	if(printing_in_setbranch) printf("getting x\n");
  	status = CPXgetx (env_just_solve_mips, mip_ob2, x2, 0, cur_numcols-1);
	if(status) 
	{
		printf ("(%d) CPXgetx, Failed to get x values, error code %d\n", __LINE__,status);
		printf("the solve status was: %d\n",lpstat);
		goto TERMINATE;
	}
	ob2_mip_objvals[0] = x2[obj1_index];
	ob2_mip_objvals[1] = x2[obj2_index];

	if(printing_in_setbranch) printf("plotting\n");
	if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",x2[obj1_index],x2[obj2_index]);
	
	if(!userhandle_up->x2) userhandle_up->x2 = calloc (cur_numcols,sizeof(double));
	if(!userhandle_down->x2) userhandle_down->x2 = calloc (cur_numcols,sizeof(double));
	
	for(i=0;i<cur_numcols;i++) 
	{
		userhandle_up->x2[i] = x2[i];
		userhandle_down->x2[i] = x2[i];
	}
	
/*	for(i=0;i<total_num_integer;i++) printf("x%d: %lf\n",integer_indices[i],x2[integer_indices[i]]);*/

	if(lpstat == 101 || lpstat == 102)
	{
		if(printing_in_setbranch) printf("the ob2 mip solution was optimal\n");
		ob2_mip_opt = 1;
		
		if(ob2_mip_objvals[0] > reduced_subprob_x_lb || ob2_mip_objvals[1] < reduced_subprob_y_ub)
		{
			reduced_subprob_x_lb = ob2_mip_objvals[0];
			if(there_will_only_be_points && integer_bb && integer_objective == 1) reduced_subprob_x_lb += smallest_coef;
			reduced_subprob_y_ub = ob2_mip_objvals[1];
			if(there_will_only_be_points && integer_bb && integer_objective == 2) reduced_subprob_y_ub -= smallest_coef;
			if(printing_in_setbranch) printf("after changing:\n");
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
					reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_ub,reduced_subprob_y_ub);
			bound_reduction = 1;
		}
		
		add_check = mock_insert(1,x2[obj1_index],x2[obj2_index],0,0,0,&tree);
		if(add_check)
		{
			if(printing_in_setbranch) printf("adding solution\n");
			for(i=0;i<cur_numcols;i++)
	      		{
		      		stored_x[x_rotation][i] = x2[i];
		      	}
		      	x_rotation = (x_rotation + 1) % num_x_to_store;
	 		add_check = 0;
	     		if(check_for_stopping_PSA_full)
			{
				check_for_stopping_PSA_full = 0;
				PSA_full(env,NULL,x2,NULL,NULL);
				check_for_stopping_PSA_full = 1;
			}
			else PSA_full(env,NULL,x2,NULL,NULL);
		}
		
		if(going_back) goto BRANCHING;
		
/*		printf("points_only: %d, ob2_m_o[1]: %lf, red_s_y_lb: %lf, difference: %lf\n",points_only,ob2_mip_objvals[1],reduced_subprob_y_lb,ob2_mip_objvals[1] - reduced_subprob_y_lb);*/
/*		exit(0);*/
		
		if(ob2_mip_objvals[1] - reduced_subprob_y_lb <= .0000001 || (points_only && ob2_mip_objvals[1] - reduced_subprob_y_lb <= .01001) || 
			fabs(reduced_subprob_x_lb - reduced_subprob_x_ub) < .00001 || fabs(reduced_subprob_y_lb - reduced_subprob_y_ub) <= .00001 )
		{
			fathoming = 1;
			*useraction_p = CPX_CALLBACK_SET;
	  		nodecnt = 0;
	  		goto TERMINATE;
		}
		else if(points_only)
		{
/*			printf("here :) \n");*/
/*			if(points_only > 1) exit(0);*/
			if(ob2_mip_objvals[0] - reduced_subprob_x_lb < -.00001)
			{
				ob2_mip_ran = 0;
				ob2_mip_opt = 0;
				goto BEGINNING;
			}
			int chg_indices[2] = {obj1_index,obj2_index};
			char chg_lu[2] = {'L','U'};
			double chg_bds[2] = {x2[obj1_index] + .01, x2[obj2_index]};
			status = CPXchgbds (env_just_solve_mips, mip_ob2, 2, chg_indices, chg_lu, chg_bds);
    			if ( status ){
     			 	printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
     			 	goto TERMINATE;}
     			status = CPXchgbds (env, nodelp2, 2, chg_indices, chg_lu, chg_bds);
    			if ( status ){
     			 	printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
     			 	goto TERMINATE;}
/*     			points_only = 0;*/
/*			points_only++;*/
     			ob2_mip_opt = 0;
     			ob2_mip_ran = 0;
			goto BEGINNING;
		}
		
		add_check = mock_insert(1,x_ws[obj1_index],x2[obj2_index],0,0,0,&tree);
		if(add_check)
		{
			if(printing_in_setbranch) printf("the left partial ideal is not dominated\n");
			projection = (1./slope)*(x2[obj2_index]-x_ws[obj2_index])+x_ws[obj1_index];
			add_check = mock_insert(2,x_ws[obj1_index],x_ws[obj2_index],projection,x2[obj2_index],slope,&tree);
			if(add_check)
			{
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],x2[obj2_index]);
				if(printing_in_setbranch) printf("the left partial ideal segment is not dominated\n");
				goto BRANCHING;
			}
			else
			{
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],x2[obj2_index]);
				if(printing_in_setbranch) printf("Left partial ideal segment is dominated. Fathom (%d)\n",__LINE__);
				if(right_pt_dom || right_side_dom) fathomed_by_dominated_local_ideal_pts++;
				else if(right_seg_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
				else
				{
					reduced_subprob_x_lb = x_ws[obj1_index];
	  				reduced_subprob_y_ub = x_ws[obj2_index];
	  				if(printing_in_setbranch) printf("after changing:\n");
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
							reduced_subprob_y_lb,reduced_subprob_y_lb);					
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_ub,reduced_subprob_y_ub);
	  				bound_reduction = 1;
	  				goto BRANCHING;
				}
				fathoming = 1;
				*useraction_p = CPX_CALLBACK_SET;
		  		nodecnt = 0;
		  		goto TERMINATE;
			}
		}
		else
		{
			left_side_dom = 1;
			if(printing_in_setbranch) printf("the left partial ideal is dominated. Fathom\n");
			if(right_pt_dom || right_side_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
			else if(right_seg_dom) fathomed_by_dominated_local_ideal_segments++;
			else
			{
				reduced_subprob_x_lb = x_ws[obj1_index];
  				reduced_subprob_y_ub = x_ws[obj2_index];
  				if(printing_in_setbranch) printf("after changing:\n");
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_ub,reduced_subprob_y_ub);
  				bound_reduction = 1;
  				goto BRANCHING;
			}
			fathoming = 1;
			*useraction_p = CPX_CALLBACK_SET;
		  	nodecnt = 0;
		  	goto TERMINATE;
		}

	}
	else if(lpstat == 105 || lpstat == 107 || lpstat == 109 || lpstat == 111 || lpstat == 113 || lpstat == 106)
	{
		if(printing_in_setbranch) printf("the mip solution was less than optimal\n");
		ob2_mip_opt = 0;
		add_check = mock_insert(1,x2[obj1_index],x2[obj2_index],0,0,0,&tree);
		if(add_check)
		{
			if(printing_in_setbranch) printf("adding solution\n");
			for(i=0;i<cur_numcols;i++)
	      		{
	      			stored_x[x_rotation][i] = x2[i];
	      		}
		      	x_rotation = (x_rotation + 1) % num_x_to_store;
		      	add_check = 0;
			if(check_for_stopping_PSA_full)
			{
				check_for_stopping_PSA_full = 0;
				PSA_full(env,NULL,x2,NULL,NULL);
				check_for_stopping_PSA_full = 1;
			}
			else PSA_full(env,NULL,x2,NULL,NULL);
		}
		
		if(going_back) goto BRANCHING;
		
		status = CPXgetbestobjval (env_just_solve_mips, mip_ob2, &best_bound);
		
		double obj_off;
	    	status = CPXgetobjoffset(env_just_solve_mips, mip_ob2, &obj_off );
		best_bound -= obj_off;
		reduced_subprob_y_ub = fmin(reduced_subprob_y_ub, best_bound);
		
		goto BRANCHING;
		
		reduced_subprob_y_ub = best_bound;
		projection = (-.001)*(ob2_lp_objvals[0]-best_bound);
		if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",ob2_lp_objvals[0],projection);
		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",-65.,best_bound-5./(-.001), (-.001)*(-best_bound-65),-5.);

		add_check = mock_insert(1,ob1_lp_objvals[1],projection,0,0,0,&tree);
		if(add_check)
		{
			if(printing_in_setbranch) printf("projection onto the level curve is not dominated\n");
			goto BRANCHING;
		}
		else
		{
/*			if(printing_in_setbranch) */
			printf("projection onto the ob2 level curve is dominated\n");
			exit(0);
		}
	}
	else
	{
		printf("the status of the FAILED ob2 mipsolve: %d\n",lpstat);
		exit(0);
	}
	
	/************** Here we start the branching process ********************************************/	
	
	BRANCHING:
	
	finish_time = clock();
	
	time_processing_nodes += (double)(finish_time - start_time) / CLOCKS_PER_SEC;
	
	start_time = clock();
	
	if(printing_in_setbranch) printf("at branching\n");	
	
	if(exact_mips) goto SKIP_THIS9;
	if(depth < 3)
	{
		if(printing_in_setbranch) printf("depth: %d\n",depth);
		if(!ws_mip_ran && !ws_lp_int_feas)
		{
			if(printing_in_setbranch) printf("haven't solved ws mip, going back\n");
			going_back = 1;
			goto SOLVE_WS_MIP;
		}
		else if(!ob1_mip_ran && !ob1_lp_int_feas && !right_side_dom && !right_pt_dom && !right_seg_dom)
		{
			if(printing_in_setbranch) printf("haven't solved ob1 mip, going back\n");
			going_back = 1;
			goto SOLVE_OB1_MIP;
		}
		else if(!ob2_mip_ran && !ob2_lp_int_feas && !left_side_dom)
		{
			if(printing_in_setbranch) printf("haven't solved ob2 mip, going back\n");
			going_back = 1;
			goto SOLVE_OB2_MIP;
		}
	}
	
	SKIP_THIS9:
	
	if(exploit_objective_gaps && seqnum == 0 && fathoming != 1)
	{
		if(printing_in_setbranch) printf("depth: %d\n",depth);
		if(!ws_mip_ran && !ws_lp_int_feas)
		{
			if(printing_in_setbranch) printf("haven't solved ws mip, going back\n");
			going_back = 1;
			goto SOLVE_WS_MIP;
		}
		else if(!ob1_mip_ran && !ob1_lp_int_feas && !right_side_dom && !right_pt_dom && !right_seg_dom)
		{
			if(printing_in_setbranch) printf("haven't solved ob1 mip, going back\n");
			going_back = 1;
			goto SOLVE_OB1_MIP;
		}
		else if(!ob2_mip_ran && !ob2_lp_int_feas && !left_side_dom)
		{
			if(printing_in_setbranch) printf("haven't solved ob2 mip, going back\n");
			going_back = 1;
			goto SOLVE_OB2_MIP;
		}
	}
	
	if(printing_in_setbranch) printf("skip this 9\n");
	
	if(!ob1_lp_been_solved && !ob1_mip_ran && !right_side_dom)
	{
		if(exact_mips || ws_mip_opt ) 
		{
			if(printing_in_setbranch) printf("going back to solve ob1 lp\n");
			going_back_for_lp1 = 1;
			goto SOLVE_OB1_LP;
		}
	}
	if(!ob2_lp_been_solved && !ob2_mip_ran && !left_side_dom)
	{
		if(exact_mips || (ws_mip_opt && ob1_mip_opt && (right_pt_dom || right_side_dom || right_seg_dom))) 
		{
			if(printing_in_setbranch) printf("going back to solve ob2 lp\n");
			going_back_for_lp2 = 1;
			goto SOLVE_OB2_LP;
		}
	}
	
	if(ws_mip_ran && ob1_mip_ran && !ob2_mip_ran && !left_side_dom) goto SOLVE_OB2_MIP;
	
	COPY_MIPSTARTS:

	if(from_pareto) goto PARETO_BRANCH;
	
	double x_vals[2] = {reduced_subprob_x_lb,reduced_subprob_x_ub};
	double y_vals[2] = {reduced_subprob_y_lb,reduced_subprob_y_ub};
	if(bound_reduction)
	{
		/******** Here we attempt to reduce portions of the search region which are dominated by the primal solutions **********/
		if(printing_in_setbranch) printf("attempting to reduce search region dimensions, originally: %lf,%lf,%lf,%lf\n",x_vals[0],x_vals[1],
														y_vals[0],y_vals[1]);
		reduce_box(x_vals,y_vals, tree);
		if(x_vals[0] > reduced_subprob_x_lb || y_vals[0] > reduced_subprob_y_lb) 
		{
			if(printing_in_setbranch) printf("reduced box bound before branching\n");
			if(printing_in_setbranch) printf("new vals: %lf,%lf,%lf,%lf\n",x_vals[0],x_vals[1],y_vals[0],y_vals[1]);
			reduced_subprob_x_lb = x_vals[0];
			reduced_subprob_x_ub = x_vals[1];
			reduced_subprob_y_lb = y_vals[0];
			reduced_subprob_y_ub = y_vals[1];
/*			bound_reduction = 1;*/
			if(printing_in_setbranch) printf("after changing:\n");
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
					reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_ub,reduced_subprob_y_ub);
		}
		int ind4[4] = {obj1_index,obj1_index,obj2_index,obj2_index};
		char lu4[4] = {'L','U','L','U'};
		double bds4[4] = {reduced_subprob_x_lb,reduced_subprob_x_ub,reduced_subprob_y_lb,reduced_subprob_y_ub};
		status = CPXchgbds (env, nodelp2, 4, ind4, lu4, bds4);
		if ( status ) {
	   		printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
	  	}
	}
	
	SKIP_THIS10:
	if(printing_in_setbranch) printf("skip this 10\n");
	if(ws_mip_opt && ob1_mip_opt && ob2_mip_opt)
	{
		int all_same = 1;
		for(i=0;i<total_num_integer;i++)
		{
			k = integer_indices[i];
			if( fabs(x_ws[k] - x1[k]) > .00001 || fabs(x_ws[k] - x2[k]) > .00001 )
			{
				all_same = 0;
				break;
			}
		}
		if(all_same)
		{
			if(printing_in_setbranch) printf("all 3 mips optimal and have same values of integer variables. Run PSA_full, then Fathom.\n");
			if(check_for_stopping_PSA_full)
			{
				check_for_stopping_PSA_full = 0;
				PSA_full(env,NULL,x1,NULL,NULL);
				check_for_stopping_PSA_full = 1;
			}
			else PSA_full(env,NULL,x1,NULL,NULL);
			*useraction_p = CPX_CALLBACK_SET;
			fathoming = 1;
	  		nodecnt = 0;
	  		goto TERMINATE;
		}
	}
	else if(ws_mip_opt && ob1_mip_opt && left_side_dom)
	{
		int all_same = 1;
		for(i=0;i<total_num_integer;i++)
		{
			k = integer_indices[i];
			if( fabs(x_ws[k] - x1[k]) > .00001)
			{
				all_same = 0;
				break;
			}
		}
		if(all_same)
		{
			if(printing_in_setbranch) printf("Left side was dominated and ws and ob1 mips optimal and have same values of integer variables. Run PSA_full, then Fathom.\n");
			if(check_for_stopping_PSA_full)
			{
				check_for_stopping_PSA_full = 0;
				PSA_full(env,NULL,x1,NULL,NULL);
				check_for_stopping_PSA_full = 1;
			}
			else PSA_full(env,NULL,x1,NULL,NULL);
			*useraction_p = CPX_CALLBACK_SET;
			fathoming = 1;
	  		nodecnt = 0;
	  		goto TERMINATE;
		}
	}
	else if(ws_mip_opt && ob2_mip_opt && right_side_dom)
	{
		int all_same = 1;
		for(i=0;i<total_num_integer;i++)
		{
			k = integer_indices[i];
			if(fabs(x_ws[k] - x2[k]) > .00001 )
			{
				all_same = 0;
				break;
			}
		}
		if(all_same)
		{
			if(printing_in_setbranch) printf("right side was dominated and ws and ob2 mips optimal and have same values of integer variables. Run PSA_full, then Fathom.\n");
			PSA_full(env,NULL,x_ws,NULL,NULL);
			*useraction_p = CPX_CALLBACK_SET;
			fathoming = 1;
	  		nodecnt = 0;
	  		goto TERMINATE;
		}
	}
	
	int temp_pareto = 0; // This used to be implemented differently, but now we should never enter the pareto branching portion of the code unless SENT here.
	if(temp_pareto)
	{
		PARETO_BRANCH:
		number_pareto_branches++;
		if(!exact_mips && !mipstarts_copied)
		{
			from_pareto = 1;
			goto COPY_MIPSTARTS;
		}
		*useraction_p = CPX_CALLBACK_SET;
		if(printing_in_setbranch) printf("pareto branching on seqnum: %d\n",seqnum);
		
		double x_vals[2] = {reduced_subprob_x_lb,sub_pr1_x_ub};
		double y_vals[2] = {sub_pr1_y_lb,reduced_subprob_y_ub};
		int time2 = 0;
		
		if(!bound_reduction)			
		{
			if(1) //(!exact_mips)
			{
				reduce_box(x_vals,y_vals, tree);
				if(x_vals[0] > reduced_subprob_x_lb || y_vals[0] > sub_pr1_y_lb) 
				{
					reduced_subprob_x_lb = x_vals[0];
					sub_pr1_x_ub = x_vals[1];
					sub_pr1_y_lb = y_vals[0];
					reduced_subprob_y_ub = y_vals[1];
					bound_reduction = 1;
					goto DO_BOUND_REDUCTION;
				}
/*				else*/
/*				{*/
/*					sub_pr1_x_ub = reduced_subprob_x_ub;*/
/*					sub_pr1_y_lb = reduced_subprob_y_lb;*/
/*				}*/
			
				x_vals[0] = sub_pr2_x_lb;
				x_vals[1] = reduced_subprob_x_ub;
				y_vals[0] = reduced_subprob_y_lb;
				y_vals[1] = sub_pr2_y_ub;
				reduce_box(x_vals,y_vals, tree);
				if(x_vals[0] > sub_pr2_x_lb || y_vals[0] > reduced_subprob_y_lb) 
				{
					sub_pr2_x_lb = x_vals[0];
					reduced_subprob_x_ub = x_vals[1];
					reduced_subprob_y_lb = y_vals[0];
					sub_pr2_y_ub = y_vals[1];
					bound_reduction = 1;
					time2 = 1;
					goto DO_BOUND_REDUCTION;
				}
			}
			
			*useraction_p = CPX_CALLBACK_SET;
		
			vars2[0] = obj1_index;
			vars2[1] = obj2_index;
			varlu2[0] = 'U';
			varlu2[1] = 'L';
			varbd2[0] = sub_pr1_x_ub;
			varbd2[1] = sub_pr1_y_lb;
/*			if(!exact_mips) */
			status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 2, vars2, 
								varlu2, varbd2, objval, userhandle_up, &seqnum1);
			if(status)
			{
				printf("Failed to create up branch (line %d). Error code %d\n",__LINE__,status);
				goto TERMINATE;
			}
								
			x_vals[0] = sub_pr2_x_lb;
			x_vals[1] = reduced_subprob_x_ub;
			y_vals[0] = reduced_subprob_y_lb;
			y_vals[1] = sub_pr2_y_ub;
								
			varlu2[0] = 'L';
			varlu2[1] = 'U';
			varbd2[0] = sub_pr2_x_lb;
			varbd2[1] = sub_pr2_y_ub;
/*			if(!exact_mips) */
			status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 2, vars2, 
								varlu2, varbd2, objval, userhandle_down, &seqnum2);
			if(status)
			{
				printf("Failed to create down branch (line %d). Error code %d\n",__LINE__,status);
				goto TERMINATE;
			}
		}
		else
		{
			if(printing_in_setbranch) printf("also reducing search region\n");
			
			x_vals[0] = reduced_subprob_x_lb;
			x_vals[1] = sub_pr1_x_ub;
			y_vals[0] = sub_pr1_y_lb;
			y_vals[1] = reduced_subprob_y_ub;
			reduce_box(x_vals,y_vals, tree);
			if(x_vals[0] > reduced_subprob_x_lb || y_vals[0] > sub_pr1_y_lb) 
			{
				reduced_subprob_x_lb = x_vals[0];
				sub_pr1_x_ub = x_vals[1];
				sub_pr1_y_lb = y_vals[0];
				reduced_subprob_y_ub = y_vals[1];
				bound_reduction = 1;
			}
			
			DO_BOUND_REDUCTION:
			
			*useraction_p = CPX_CALLBACK_SET;
			
			vars3[0] = obj1_index;
			vars3[1] = obj1_index;
			vars3[2] = obj2_index;
			vars3[3] = obj2_index;
			varlu3[0] = 'U';
			varlu3[1] = 'L';
			varlu3[2] = 'L';
			varlu3[3] = 'U';
			varbd3[0] = sub_pr1_x_ub;
			varbd3[1] = reduced_subprob_x_lb;
			varbd3[2] = sub_pr1_y_lb;
			varbd3[3] = reduced_subprob_y_ub;
			status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 4, vars3, 
								varlu3, varbd3, objval, userhandle_up, &seqnum1);
			if(status)
			{
				printf("Failed to create up branch (line %d). Error code %d\n",__LINE__,status);
				goto TERMINATE;
			}
				
			if(!time2)
			{				
				x_vals[0] = sub_pr2_x_lb;
				x_vals[1] = reduced_subprob_x_ub;
				y_vals[0] = reduced_subprob_y_lb;
				y_vals[1] = sub_pr2_y_ub;
				reduce_box(x_vals,y_vals, tree);
				if(x_vals[0] > sub_pr2_x_lb || y_vals[0] > reduced_subprob_y_lb) 
				{
					sub_pr2_x_lb = x_vals[0];
					reduced_subprob_x_ub = x_vals[1];
					reduced_subprob_y_lb = y_vals[0];
					sub_pr2_y_ub = y_vals[1];
	/*				printf("after changing: %lf,%lf,%lf,%lf\n",sub_pr2_x_lb,reduced_subprob_x_ub,reduced_subprob_y_lb,sub_pr2_y_ub);*/
					bound_reduction = 1;
				}
			}
												
			varlu3[0] = 'L';
			varlu3[1] = 'U';
			varlu3[2] = 'U';
			varlu3[3] = 'L';
			varbd3[0] = sub_pr2_x_lb;
			varbd3[1] = reduced_subprob_x_ub;
			varbd3[2] = sub_pr2_y_ub;
			varbd3[3] = reduced_subprob_y_lb;
			status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 4, vars3, 
								varlu3, varbd3, objval, userhandle_down, &seqnum2);
			if(status)
			{
				printf("Failed to create down branch (line %d). Error code %d\n",__LINE__,status);
				goto TERMINATE;
			}
		}
		goto TERMINATE;
		
	}
	
	int desired_index_exists = 0;

	double add_val = 2.;
	int num_found_ws_1 = 0, num_found_ws_2 = 0, num_found_1_2 = 0, doing_it_again = 0;
	
	if(!ws_mip_opt || !ob1_mip_opt || !ob2_mip_opt) add_val = 3.;

	if(nodecnt ==2 && bdcnt == 2)
	{
		//if(num_x_frac_at_both > 0 || psa_int_pivot > 0)//;frac_index = -1;
		int found_one = 0;
		if(printing_in_setbranch) printf("remember to fix this so branching occurs correctly for mip situations\n");
		if(mip_solved) // || frac_index == -1) //(seqnum == branch_seqnum)
		{
			if(printing_in_setbranch) printf("performing user branching\n");
/*			printf("d_i_e: %d, w_m_o: %d, w_l_i_f: %d, o1_m_o: %d, o1_l_i_f: %d, o2_m_o: %d, o2_l_i_f: %d\n",desired_index_exists,ws_mip_opt,*/
/*				ws_lp_int_feas,ob1_mip_opt,ob1_lp_int_feas,ob2_mip_opt,ob2_lp_int_feas);*/
			if(!desired_index_exists && !right_side_dom && (ws_mip_opt || ws_lp_int_feas) && (ob1_mip_opt || ob1_lp_int_feas))
			{
				if(printing_in_setbranch) printf("ws and ob1 mips optimal. Looking for a branching var that will cause at most 2 of the solns to be optimal at each child\n");
				AGAIN1:
				for(i=0;i<total_num_integer;i++)
				{
/*					if(x1) printf("x1_%d: %lf\t",i,x1[i]);*/
/*					if(x_ws) printf("xws_%d: %lf\t",i,x_ws[i]);*/
/*					if(x2) printf("x2_%d: %lf\t",i,x2[i]);*/
					k = integer_indices[i];
					if(fabs(x_ws[k] - x1[k]) >= .00001)
					{
						if(printing_in_setbranch) printf("x_ws%d: %lf, x1_%d: %lf\n",k,x_ws[k],k,x1[k]);
						frac_index = k;
						frac_val = fmax(x_ws[k],x1[k])-.5;
						frac_values[k] = frac_val;
						if(frac_scores[i] > .0001)
				  		{
					  		frac_scores[i] += add_val;
			/*		  		num_frac++;*/
					  	}
					  	else
					  	{
					  		frac_scores[i] += add_val + multiplier*k;
				  			num_frac++;
					  	}
					  	num_found_ws_1++;
/*						found_one = 1;*/
/*						break;*/
					}
				}
				if(doing_it_again) goto OUTTA_HERE;
			}
			if(!desired_index_exists && !found_one && !left_side_dom && (ws_mip_opt || ws_lp_int_feas) && (ob2_mip_opt || ob2_lp_int_feas))
			{
				if(printing_in_setbranch) printf("ws and ob2 mips optimal. Looking for a branching var that will cause at most 2 of the solns to be optimal at each child\n");
				AGAIN2:
				for(i=0;i<total_num_integer;i++)
				{
					k = integer_indices[i];
					if(fabs(x_ws[k] - x2[k]) >= .00001)
					{
						if(printing_in_setbranch) printf("x_ws%d: %lf, x2_%d: %lf\n",k,x_ws[k],k,x2[k]);
						frac_index = k;
						frac_val = fmax(x_ws[k],x2[k])-.5;
						frac_values[k] = frac_val;
						if(frac_scores[i] > .0001)
				  		{
					  		frac_scores[i] += add_val;
			/*		  		num_frac++;*/
					  	}
					  	else
					  	{
					  		frac_scores[i] += add_val + multiplier*k;
				  			num_frac++;
					  	}
					  	num_found_ws_2++;
/*						found_one = 1;*/
/*						break;*/
					}
				}
				if(doing_it_again) goto OUTTA_HERE;
			}
			if(!desired_index_exists && !found_one && !left_side_dom && !right_side_dom && (ob1_mip_opt || ob1_lp_int_feas) && 
				(ob2_mip_opt || ob2_lp_int_feas))
			{
				if(printing_in_setbranch) printf("ob1 and ob2 mips optimal. Looking for a branching var that will cause at most 2 of the solns to be optimal at each child\n");
				AGAIN3:
				for(i=0;i<total_num_integer;i++)
				{
					k = integer_indices[i];
					if(fabs(x1[k] - x2[k]) >= .00001)
					{
						if(printing_in_setbranch) printf("x1%d: %lf, x1_%d: %lf\n",k,x1[k],k,x2[k]);
						frac_index = k;
						frac_val = fmax(x1[k],x2[k])-.5;
						frac_values[k] = frac_val;
						if(frac_scores[i] > .0001)
				  		{
					  		frac_scores[i] += add_val;
			/*		  		num_frac++;*/
					  	}
					  	else
					  	{
					  		frac_scores[i] += add_val + multiplier*k;
				  			num_frac++;
					  	}
					  	num_found_1_2++;
/*						found_one = 1;*/
/*						break;*/
					}
				}
				if(doing_it_again) goto OUTTA_HERE;
			}
		}
		if(num_found_ws_1 < num_found_ws_2 && num_found_ws_1 < num_found_1_2)
		{
			add_val++;
			doing_it_again = 1;
			goto AGAIN1;
		}
		if(num_found_ws_2 < num_found_ws_1 && num_found_ws_2 < num_found_1_2)
		{
			add_val++;
			doing_it_again = 1;
			goto AGAIN2;
		}
		if(num_found_1_2 < num_found_ws_1 && num_found_1_2 < num_found_ws_2)
		{
			add_val++;
			doing_it_again = 1;
			goto AGAIN3;
		}
		
		OUTTA_HERE:
		
		if(!desired_index_exists && num_frac)
		{
			if(printing_in_setbranch) printf("choosing branch variable based on lps\n");
			qsort(frac_scores, total_num_integer, sizeof(double), comparison_function3);
/*			for(i=0;i<total_num_integer;i++) printf("%lf\n",frac_scores[i]);*/
			int index = (int) (1./multiplier * (frac_scores[0] - (int) frac_scores[0])+.00001);
			frac_index = index;
			frac_val = frac_values[index];
			if(printing_in_setbranch) printf("index: %d, val: %lf\n",index,frac_val);
			found_one = 1;
		}
		else if(!desired_index_exists && !found_one) goto BRANCHASCPLEX;
		
		if(frac_index != -1)
		{		
			if(printing_in_setbranch) printf("there is an index we want to use for branching: %d\n",frac_index);	
			double ub2[1] = {0.};
			double lb2[1] = {0.};
			status = CPXgetub (env, nodelp, ub2, frac_index, frac_index);
			if ( status ) {
				printf ("(%d) Failed to get ub's. Error code %d\n", __LINE__,status);
				printf("cur_numcols: %d\tfrac_index: %d\n",cur_numcols,frac_index);
				exit(0);
			}
			status = CPXgetlb (env, nodelp, lb2, frac_index, frac_index);
			if ( status ) {
				printf ("(%d) Failed to get lb's. Error code %d\n", __LINE__,status);
				printf("cur_numcols: %d\tfrac_index: %d\n",cur_numcols,frac_index);
				exit(0);
			}
			if(ub2[0] == lb2[0])
			{
				if(printing_in_setbranch) printf("(%d) using CPLEX default branch\n", __LINE__);
/*				printf("indices[0]: %d\n",indices[0]);*/
				goto BRANCHASCPLEX;
			}
		
			if(!bound_reduction)
			{
				if(printing_in_setbranch) printf("search region size hasn't been reduced\n");
				if(bd_reduction_during_branching)
				{
					if(printing_in_setbranch) printf("attempting to reduce bounds\n");
					status = CPXgetub (env, nodelp2, ub_, 0, cur_numcols-1);
					status = CPXgetlb (env, nodelp2, lb_, 0, cur_numcols-1);
				
					k = 0;
					for(i=0;i<cur_numcols;i++)
					{
						if(i == integer_indices[k])
						{
							low_up[2*k] = 'L';
							low_up[2*k+1] = 'U';
							bds_br1[2*k] = lb_[i];
							bds_br2[2*k] = lb_[i];
							bds_br1[2*k+1] = ub_[i];
							bds_br2[2*k+1] = ub_[i];
							br_ind[2*k] = integer_indices[k];
							br_ind[2*k+1] = integer_indices[k];
/*							printf("x_%d: %lf to %lf \t %lf to %lf\n",i,bds_br1[2*k],bds_br1[2*k+1],bds_br2[2*k],bds_br2[2*k+1]);*/
							k++;
						}
						if(k >= total_num_integer) break;
					}
				
					bd_red_start_time = clock();
					if(printing_in_setbranch) printf("doing bound reduction\n");
					int reduced_worked = reduce_bds_before_branching( (CPXENVptr) env, nodelp2, bds_br1, bds_br2, 0, frac_index, frac_val, 
											(ub_[obj2_index]-lb_[obj2_index])/(lb_[obj1_index]-ub_[obj1_index]));
					bd_red_finish_time = clock();
					time_tightening_variable_bounds += (double)(bd_red_finish_time - bd_red_start_time) / CLOCKS_PER_SEC;
					
/*					k = 0;*/
/*					for(i=0;i<cur_numcols;i++)*/
/*					{*/
/*						if(i == integer_indices[k])*/
/*						{*/
/*							printf("x_%d: %lf to %lf \t %lf to %lf\n",i,bds_br1[2*k],bds_br1[2*k+1],bds_br2[2*k],bds_br2[2*k+1]);*/
/*							k++;*/
/*						}*/
/*						if(k >= total_num_integer) break;*/
/*					}*/
									
					*useraction_p = CPX_CALLBACK_SET;
				
					if(reduced_worked)
					{
						if(printing_in_setbranch) printf("branching at line %d\n",__LINE__);
/*						printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
						if(printing_in_setbranch) printf("reduction worked\n");
						status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 2*total_num_integer, br_ind, 
										low_up, bds_br1, objval, userhandle_up, &seqnum1);
			   			if ( status )
			   			{  
/*			   				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
			   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
			   				goto TERMINATE;
			   			}
			   			status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 2*total_num_integer, br_ind, 
											low_up, bds_br2, objval, userhandle_down, &seqnum2);
			   			if ( status )
			   			{  
/*			   				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
			   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
			   				goto TERMINATE;
			   			}
			   			
			   			if(printing_in_setbranch) printf("successfully created branches: %d and %d\n",seqnum1,seqnum2);
				   		
			   			goto TERMINATE;
		   			}
		   			else
		   			{
		   				if(printing_in_setbranch) printf("it didn't work\n");
		   				goto BRANCHASCPLEX;
		   			}
	   			}
				
/*				printf("changing branch variable to %d\n",frac_index);*/

				*useraction_p = CPX_CALLBACK_SET;

				vars[0] = frac_index;
				varlu[0] = 'U';
				varbd[0] = (int) frac_val;
	
				if(printing_in_setbranch) printf("branching at line %d\n",__LINE__);
/*				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
				status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 1, vars, 
									varlu, varbd, objval, userhandle_up, &seqnum1);
	   			if ( status )
	   			{  
/*	   				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
	   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
	   				goto TERMINATE;
	   			}
	   			
				varlu[0] = 'L';
				varbd[0] = varbd[0]+1;
				status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 1, vars, 
									varlu, varbd, objval, userhandle_down, &seqnum2);
	   			if ( status )
	   			{
	   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
	   				goto TERMINATE;
	   			}
		   		
				frac_index = -1;
				num_x_frac_at_both = 0;
				psa_int_pivot = 0;
			}
			else
			{
				if(bd_reduction_during_branching)
				{
					status = CPXgetub (env, nodelp2, ub_, 0, cur_numcols-1);
					status = CPXgetlb (env, nodelp2, lb_, 0, cur_numcols-1);
				
					if(printing_in_setbranch) printf("changing branch variable to %d\n",frac_index);
					int i,k = 0;
					for(i=0;i<cur_numcols;i++)
					{
						if(i == integer_indices[k])
						{
							low_up[2*k] = 'L';
							low_up[2*k+1] = 'U';
							bds_br1[2*k] = lb_[i];
							bds_br2[2*k] = lb_[i];
							bds_br1[2*k+1] = ub_[i];
							bds_br2[2*k+1] = ub_[i];
							br_ind[2*k] = integer_indices[k];
							br_ind[2*k+1] = integer_indices[k];
							k++;
						}
						if(k >= total_num_integer) break;
					}
				
					low_up[2*total_num_integer] = 'L';
					low_up[2*total_num_integer+1] = 'U';
					low_up[2*total_num_integer+2] = 'L';
					low_up[2*total_num_integer+3] = 'U';
					bds_br1[2*total_num_integer] = reduced_subprob_x_lb;
					bds_br1[2*total_num_integer+1] = reduced_subprob_x_ub;
					bds_br1[2*total_num_integer+2] = reduced_subprob_y_lb;
					bds_br1[2*total_num_integer+3] = reduced_subprob_y_ub;
					bds_br2[2*total_num_integer] = reduced_subprob_x_lb;
					bds_br2[2*total_num_integer+1] = reduced_subprob_x_ub;
					bds_br2[2*total_num_integer+2] = reduced_subprob_y_lb;
					bds_br2[2*total_num_integer+3] = reduced_subprob_y_ub;
					br_ind[2*total_num_integer] = obj1_index;
					br_ind[2*total_num_integer+1] = obj1_index;
					br_ind[2*total_num_integer+2] = obj2_index;
					br_ind[2*total_num_integer+3] = obj2_index;
				
					bd_red_start_time = clock();
					if(printing_in_setbranch) printf("doing bound reduction\n");
					
/*					printf("$$$$$$$$$$$$$$$$$$$\n");*/
/*					for(i=0;i<total_num_integer;i++) printf("x%d, bd1: %lf to %lf\t bd2: %lf to %lf\n",integer_indices[i],bds_br1[2*i],bds_br1[2*i+1],bds_br2[2*i],bds_br2[2*i+1]);*/
/*					printf("$$$$$$$$$$$$$$$$$$$\n");*/
					
					int reduced_worked = reduce_bds_before_branching( (CPXENVptr) env, nodelp2, bds_br1, bds_br2, 0, frac_index, frac_val,
											(ub_[obj2_index]-lb_[obj2_index])/(lb_[obj1_index]-ub_[obj1_index]));
					bd_red_finish_time = clock();
					time_tightening_variable_bounds += (double)(bd_red_finish_time - bd_red_start_time) / CLOCKS_PER_SEC;
										
					*useraction_p = CPX_CALLBACK_SET;
				
					if(reduced_worked)
					{
/*						printf("$$$$$$$$$$$$$$$$$$$\n");*/
/*						for(i=0;i<total_num_integer;i++) printf("x%d, bd1: %lf to %lf\t bd2: %lf to %lf\n",integer_indices[i],bds_br1[2*i],bds_br1[2*i+1],bds_br2[2*i],bds_br2[2*i+1]);*/
/*						printf("$$$$$$$$$$$$$$$$$$$\n");*/
						if(printing_in_setbranch) printf("branching at line %d\n",__LINE__);
/*						printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
						status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 2*total_num_integer+4, br_ind, 
											low_up, bds_br1, objval, userhandle_up, &seqnum1);
			   			if ( status )
			   			{  
/*			   				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
			   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
			   				goto TERMINATE;
			   			}
			   			status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 2*total_num_integer+4, br_ind, 
											low_up, bds_br2, objval, userhandle_down, &seqnum2);
			   			if ( status )
			   			{  
/*			   				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
			   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
			   				goto TERMINATE;
			   			}
				   		goto TERMINATE;
		   			}
	   			}
	   			
/*				printf("changing branch variable to %d and also reducing subproblem size\n",frac_index);*/

				*useraction_p = CPX_CALLBACK_SET;

				if(printing_in_setbranch) printf("branching at line %d\n",__LINE__);
/*				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
				vars3[0] = frac_index;
				varlu3[0] = 'U';
				varbd3[0] = (int) frac_val;
/*				printf("frac val: %lf, low: %lf, up: %lf\n",frac_val,varbd3[0],varbd3[0]+1);*/
				vars3[1] = obj1_index;
				vars3[2] = obj1_index;
				vars3[3] = obj2_index;
				vars3[4] = obj2_index;
				varlu3[1] = 'L';
				varlu3[2] = 'U';
				varlu3[3] = 'L';
				varlu3[4] = 'U';
				varbd3[1] = reduced_subprob_x_lb;
				varbd3[2] = reduced_subprob_x_ub;
				varbd3[3] = reduced_subprob_y_lb;
				varbd3[4] = reduced_subprob_y_ub;
				
				status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 5, vars3, 
									varlu3, varbd3, objval, userhandle_up, &seqnum1);
	   			if ( status )
	   			{  
/*	   				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
	   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
	   				goto TERMINATE;
	   			}
	   			
				varlu3[0] = 'L';
				varbd3[0] = varbd3[0]+1;

				status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 5, vars3, 
									varlu3, varbd3, objval,userhandle_down, &seqnum2);
	   			if ( status )
	   			{
	   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
	   				goto TERMINATE;
	   			}
		   		
				frac_index = -1;
				num_x_frac_at_both = 0;
				psa_int_pivot = 0;
			}
		}
		else
		{
			BRANCHASCPLEX:
			if(printing_in_setbranch) printf("(%d) using CPLEX default branch\n", __LINE__);
			*useraction_p = CPX_CALLBACK_SET;
			
			if(printing_in_setbranch) printf("branching at line %d\n",__LINE__);
			status = CPXbranchcallbackbranchasCPLEX(env, cbdata, wherefrom, 0, userhandle_up, &seqnum1);
			if ( status )
	   		{
	   			printf("%s(%d): CPXbranchcallbackbranchasCPLEX, Failed to branch %d\n", __FILE__, __LINE__, status);
	   			goto TERMINATE;
			}
			
			status = CPXbranchcallbackbranchasCPLEX(env, cbdata, wherefrom, 1, userhandle_down, &seqnum2);
			if ( status )
	   		{
	   			printf("%s(%d): CPXbranchcallbackbranchasCPLEX, Failed to branch %d\n", __FILE__, __LINE__, status);
	   			goto TERMINATE;
			}
		}
	}
	
	
	/************** Here we start the branching process ********************************************/	
	
	BRANCHING2:
	
	if(last_cutcallback_seqnum != seqnum) goto TERMINATE;
	
	low_up = (char *) malloc ((2*total_num_integer+16)*sizeof(char));
	lb_ = (double *) malloc ((cur_numcols)*sizeof(double));
	ub_ = (double *) malloc ((cur_numcols)*sizeof(double));
	bds_br1 = (double *) malloc ((2*total_num_integer+16)*sizeof(double));
	bds_br2 = (double *) malloc ((2*total_num_integer+16)*sizeof(double));
	br_ind = (int *) malloc ((2*total_num_integer+16)*sizeof(int));
	
	if(printing_in_setbranch) printf("(%d) at branching\n",__LINE__);
	
	if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
	if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
	if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
						reduced_subprob_y_lb,reduced_subprob_y_lb);					
	if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
						reduced_subprob_y_ub,reduced_subprob_y_ub);
	
	double x_vals2[2] = {reduced_subprob_x_lb,reduced_subprob_x_ub};
	double y_vals2[2] = {reduced_subprob_y_lb,reduced_subprob_y_ub};
	if(bound_reduction)
	{
		/******** Here we attempt to reduce portions of the search region which are dominated by the primal solutions **********/
		if(printing_in_setbranch) printf("attempting to reduce search region dimensions, originally: %lf,%lf,%lf,%lf\n",x_vals2[0],x_vals2[1],
														y_vals2[0],y_vals2[1]);
		reduce_box(x_vals2,y_vals2, tree);
		if(x_vals2[0] > reduced_subprob_x_lb || y_vals[0] > reduced_subprob_y_lb) 
		{
			if(printing_in_setbranch) printf("reduced box bound before branching\n");
			if(printing_in_setbranch) printf("new vals: %lf,%lf,%lf,%lf\n",x_vals2[0],x_vals2[1],y_vals2[0],y_vals2[1]);
			reduced_subprob_x_lb = x_vals2[0];
			reduced_subprob_x_ub = x_vals2[1];
			reduced_subprob_y_lb = y_vals2[0];
			reduced_subprob_y_ub = y_vals2[1];
/*			bound_reduction = 1;*/
			if(printing_in_setbranch) printf("after changing:\n");
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
					reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_ub,reduced_subprob_y_ub);
		}
	}
	
/*	if(ws_mip_opt && ob1_mip_opt && ob2_mip_opt)*/
/*	{*/
/*		int all_same = 1;*/
/*		for(i=0;i<total_num_integer;i++)*/
/*		{*/
/*			k = integer_indices[i];*/
/*			if( fabs(x_ws[k] - x1[k]) > .00001 || fabs(x_ws[k] - x2[k]) > .00001 )*/
/*			{*/
/*				all_same = 0;*/
/*				break;*/
/*			}*/
/*		}*/
/*		if(all_same)*/
/*		{*/
/*			if(printing_in_setbranch) printf("all 3 mips optimal and have same values of integer variables. Run PSA_full, then Fathom.\n");*/
/*			if(check_for_stopping_PSA_full)*/
/*			{*/
/*				check_for_stopping_PSA_full = 0;*/
/*				PSA_full(env,NULL,x1,NULL,NULL);*/
/*				check_for_stopping_PSA_full = 1;*/
/*			}*/
/*			else PSA_full(env,NULL,x1,NULL,NULL);*/
/*			*useraction_p = CPX_CALLBACK_SET;*/
/*			fathoming = 1;*/
/*	  		nodecnt = 0;*/
/*	  		goto TERMINATE;*/
/*		}*/
/*	}*/
/*	else if(ws_mip_opt && ob1_mip_opt && left_side_dom)*/
/*	{*/
/*		int all_same = 1;*/
/*		for(i=0;i<total_num_integer;i++)*/
/*		{*/
/*			k = integer_indices[i];*/
/*			if( fabs(x_ws[k] - x1[k]) > .00001)*/
/*			{*/
/*				all_same = 0;*/
/*				break;*/
/*			}*/
/*		}*/
/*		if(all_same)*/
/*		{*/
/*			if(printing_in_setbranch) printf("Left side was dominated and ws and ob1 mips optimal and have same values of integer variables. Run PSA_full, then Fathom.\n");*/
/*			if(check_for_stopping_PSA_full)*/
/*			{*/
/*				check_for_stopping_PSA_full = 0;*/
/*				PSA_full(env,NULL,x1,NULL,NULL);*/
/*				check_for_stopping_PSA_full = 1;*/
/*			}*/
/*			else PSA_full(env,NULL,x1,NULL,NULL);*/
/*			*useraction_p = CPX_CALLBACK_SET;*/
/*			fathoming = 1;*/
/*	  		nodecnt = 0;*/
/*	  		goto TERMINATE;*/
/*		}*/
/*	}*/
/*	else if(ws_mip_opt && ob2_mip_opt && right_side_dom)*/
/*	{*/
/*		int all_same = 1;*/
/*		for(i=0;i<total_num_integer;i++)*/
/*		{*/
/*			k = integer_indices[i];*/
/*			if(fabs(x_ws[k] - x2[k]) > .00001 )*/
/*			{*/
/*				all_same = 0;*/
/*				break;*/
/*			}*/
/*		}*/
/*		if(all_same)*/
/*		{*/
/*			if(printing_in_setbranch) printf("right side was dominated and ws and ob2 mips optimal and have same values of integer variables. Run PSA_full, then Fathom.\n");*/
/*			PSA_full(env,NULL,x_ws,NULL,NULL);*/
/*			*useraction_p = CPX_CALLBACK_SET;*/
/*			fathoming = 1;*/
/*	  		nodecnt = 0;*/
/*	  		goto TERMINATE;*/
/*		}*/
/*	}*/
	
	int temp_pareto2 = 0; // This used to be implemented differently, but now we should never enter the pareto branching portion of the code unless SENT here.
	if(temp_pareto2)
	{
		PARETO_BRANCH2:
		number_pareto_branches++;
		*useraction_p = CPX_CALLBACK_SET;
		if(printing_in_setbranch) printf("pareto branching on seqnum: %d\n",seqnum);
		
		double x_vals2[2] = {reduced_subprob_x_lb,sub_pr1_x_ub};
		double y_vals2[2] = {sub_pr1_y_lb,reduced_subprob_y_ub};
		int time2 = 0;
		
/*		printf("node est's: %lf, %lf\n",nodeest[0],nodeest[1]);*/
		
		if(!bound_reduction)			
		{
			if(1) //(!exact_mips)
			{
				reduce_box(x_vals2,y_vals2, tree);
				if(x_vals2[0] > reduced_subprob_x_lb || y_vals2[0] > sub_pr1_y_lb) 
				{
					reduced_subprob_x_lb = x_vals2[0];
					sub_pr1_x_ub = x_vals2[1];
					sub_pr1_y_lb = y_vals2[0];
					reduced_subprob_y_ub = y_vals2[1];
					bound_reduction = 1;
					goto DO_BOUND_REDUCTION2;
				}
			
				x_vals2[0] = sub_pr2_x_lb;
				x_vals2[1] = reduced_subprob_x_ub;
				y_vals2[0] = reduced_subprob_y_lb;
				y_vals2[1] = sub_pr2_y_ub;
				reduce_box(x_vals2,y_vals2, tree);
				if(x_vals2[0] > sub_pr2_x_lb || y_vals2[0] > reduced_subprob_y_lb) 
				{
					sub_pr2_x_lb = x_vals2[0];
					reduced_subprob_x_ub = x_vals2[1];
					reduced_subprob_y_lb = y_vals2[0];
					sub_pr2_y_ub = y_vals2[1];
					bound_reduction = 1;
					time2 = 1;
					goto DO_BOUND_REDUCTION2;
				}
			}
			
			*useraction_p = CPX_CALLBACK_SET;
		
			vars2[0] = obj1_index;
			vars2[1] = obj2_index;
			varlu2[0] = 'U';
			varlu2[1] = 'L';
			varbd2[0] = sub_pr1_x_ub;
			varbd2[1] = sub_pr1_y_lb;
/*			if(!exact_mips) */
			status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 2, vars2, 
								varlu2, varbd2, nodeest[1], userhandle_up, &seqnum1);
			if(status)
			{
				printf("Failed to create up branch (line %d). Error code %d\n",__LINE__,status);
				goto TERMINATE;
			}
								
			x_vals2[0] = sub_pr2_x_lb;
			x_vals2[1] = reduced_subprob_x_ub;
			y_vals2[0] = reduced_subprob_y_lb;
			y_vals2[1] = sub_pr2_y_ub;
								
			varlu2[0] = 'L';
			varlu2[1] = 'U';
			varbd2[0] = sub_pr2_x_lb;
			varbd2[1] = sub_pr2_y_ub;
/*			if(!exact_mips) */
			status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 2, vars2, 
								varlu2, varbd2, nodeest[0], userhandle_down, &seqnum2);
			if(status)
			{
				printf("Failed to create down branch (line %d). Error code %d\n",__LINE__,status);
				goto TERMINATE;
			}
		}
		else
		{
			if(printing_in_setbranch) printf("also reducing search region\n");
			
			x_vals2[0] = reduced_subprob_x_lb;
			x_vals2[1] = sub_pr1_x_ub;
			y_vals2[0] = sub_pr1_y_lb;
			y_vals2[1] = reduced_subprob_y_ub;
			reduce_box(x_vals2,y_vals2, tree);
			if(x_vals2[0] > reduced_subprob_x_lb || y_vals2[0] > sub_pr1_y_lb) 
			{
				reduced_subprob_x_lb = x_vals2[0];
				sub_pr1_x_ub = x_vals2[1];
				sub_pr1_y_lb = y_vals2[0];
				reduced_subprob_y_ub = y_vals2[1];
				bound_reduction = 1;
			}
			
			DO_BOUND_REDUCTION2:
			
			*useraction_p = CPX_CALLBACK_SET;
			
			vars3[0] = obj1_index;
			vars3[1] = obj1_index;
			vars3[2] = obj2_index;
			vars3[3] = obj2_index;
			varlu3[0] = 'U';
			varlu3[1] = 'L';
			varlu3[2] = 'L';
			varlu3[3] = 'U';
			varbd3[0] = sub_pr1_x_ub;
			varbd3[1] = reduced_subprob_x_lb;
			varbd3[2] = sub_pr1_y_lb;
			varbd3[3] = reduced_subprob_y_ub;
			status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 4, vars3, 
								varlu3, varbd3, nodeest[1], userhandle_up, &seqnum1);
			if(status)
			{
				printf("Failed to create up branch (line %d). Error code %d\n",__LINE__,status);
				goto TERMINATE;
			}
				
			if(!time2)
			{				
				x_vals2[0] = sub_pr2_x_lb;
				x_vals2[1] = reduced_subprob_x_ub;
				y_vals2[0] = reduced_subprob_y_lb;
				y_vals2[1] = sub_pr2_y_ub;
				reduce_box(x_vals2,y_vals2, tree);
				if(x_vals2[0] > sub_pr2_x_lb || y_vals2[0] > reduced_subprob_y_lb) 
				{
					sub_pr2_x_lb = x_vals2[0];
					reduced_subprob_x_ub = x_vals2[1];
					reduced_subprob_y_lb = y_vals2[0];
					sub_pr2_y_ub = y_vals2[1];
	/*				printf("after changing: %lf,%lf,%lf,%lf\n",sub_pr2_x_lb,reduced_subprob_x_ub,reduced_subprob_y_lb,sub_pr2_y_ub);*/
					bound_reduction = 1;
				}
			}
						
			varlu3[0] = 'L';
			varlu3[1] = 'U';
			varlu3[2] = 'U';
			varlu3[3] = 'L';
			varbd3[0] = sub_pr2_x_lb;
			varbd3[1] = reduced_subprob_x_ub;
			varbd3[2] = sub_pr2_y_ub;
			varbd3[3] = reduced_subprob_y_lb;
			status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 4, vars3, 
								varlu3, varbd3, nodeest[0], userhandle_down, &seqnum2);
			if(status)
			{
				printf("Failed to create down branch (line %d). Error code %d\n",__LINE__,status);
				goto TERMINATE;
			}
		}
		goto TERMINATE;
		
	}
	
	int desired_index_exists2 = 0;

	double add_val2 = 2.;
	int num_found_ws_1_ = 0, num_found_ws_2_ = 0, num_found_1_2_ = 0, doing_it_again_ = 0;
	
	if(!ws_mip_opt || !ob1_mip_opt || !ob2_mip_opt) add_val2 = 3.;

	if(nodecnt ==2 && bdcnt == 2)
	{
		//if(num_x_frac_at_both > 0 || psa_int_pivot > 0)//;frac_index = -1;
		int found_one = 0;
		if(printing_in_setbranch) printf("remember to fix this so branching occurs correctly for mip situations\n");
		if(mip_solved) // || frac_index == -1) //(seqnum == branch_seqnum)
		{
			if(printing_in_setbranch) printf("performing user branching\n");
/*			printf("d_i_e: %d, w_m_o: %d, w_l_i_f: %d, o1_m_o: %d, o1_l_i_f: %d, o2_m_o: %d, o2_l_i_f: %d\n",desired_index_exists,ws_mip_opt,*/
/*				ws_lp_int_feas,ob1_mip_opt,ob1_lp_int_feas,ob2_mip_opt,ob2_lp_int_feas);*/
			if(!desired_index_exists2 && !right_side_dom && (ws_mip_opt || ws_lp_int_feas) && (ob1_mip_opt || ob1_lp_int_feas))
			{
				if(printing_in_setbranch) printf("ws and ob1 mips optimal. Looking for a branching var that will cause at most 2 of the solns to be optimal at each child\n");
				AGAIN1_:
				for(i=0;i<total_num_integer;i++)
				{
					k = integer_indices[i];
/*					if(x1) printf("x1_%d: %lf\t",k,x1[k]);*/
/*					if(x_ws) printf("xws_%d: %lf\t",k,x_ws[k]);*/
/*					if(x2) printf("x2_%d: %lf\n",k,x2[k]);*/
/*					if(userhandle_current->x1) printf("uh_x1_%d: %lf\t",k,userhandle_current->x1[k]);*/
/*					if(userhandle_current->x_ws) printf("uh_xws_%d: %lf\t",k,userhandle_current->x_ws[k]);*/
/*					if(userhandle_current->x2) printf("uh_x2_%d: %lf\n",k,userhandle_current->x2[k]);*/
					if(fabs(x_ws[k] - x1[k]) >= .00001)
					{
						if(printing_in_setbranch) printf("x_ws%d: %lf, x1_%d: %lf\n",k,x_ws[k],k,x1[k]);
						frac_index = k;
						frac_val = fmax(x_ws[k],x1[k])-.5;
						frac_values[k] = frac_val;
						if(frac_scores[i] > .0001)
				  		{
					  		frac_scores[i] += add_val2;
			/*		  		num_frac++;*/
					  	}
					  	else
					  	{
					  		frac_scores[i] += add_val2 + multiplier*k;
				  			num_frac++;
					  	}
					  	num_found_ws_1_++;
/*						found_one = 1;*/
/*						break;*/
					}
				}
				if(doing_it_again_) goto OUTTA_HERE2;
			}
			if(!desired_index_exists2 && !found_one && !left_side_dom && (ws_mip_opt || ws_lp_int_feas) && (ob2_mip_opt || ob2_lp_int_feas))
			{
				if(printing_in_setbranch) printf("ws and ob2 mips optimal. Looking for a branching var that will cause at most 2 of the solns to be optimal at each child\n");
				AGAIN2_:
				for(i=0;i<total_num_integer;i++)
				{
					k = integer_indices[i];
					if(fabs(x_ws[k] - x2[k]) >= .00001)
					{
						if(printing_in_setbranch) printf("x_ws%d: %lf, x2_%d: %lf\n",k,x_ws[k],k,x2[k]);
						frac_index = k;
						frac_val = fmax(x_ws[k],x2[k])-.5;
						frac_values[k] = frac_val;
						if(frac_scores[i] > .0001)
				  		{
					  		frac_scores[i] += add_val2;
			/*		  		num_frac++;*/
					  	}
					  	else
					  	{
					  		frac_scores[i] += add_val2 + multiplier*k;
				  			num_frac++;
					  	}
					  	num_found_ws_2_++;
/*						found_one = 1;*/
/*						break;*/
					}
				}
				if(doing_it_again_) goto OUTTA_HERE2;
			}
			if(!desired_index_exists2 && !found_one && !right_side_dom && !left_side_dom && (ob1_mip_opt || ob1_lp_int_feas) && 
				(ob2_mip_opt || ob2_lp_int_feas))
			{
				if(printing_in_setbranch) printf("ob1 and ob2 mips optimal. Looking for a branching var that will cause at most 2 of the solns to be optimal at each child\n");
				AGAIN3_:
				for(i=0;i<total_num_integer;i++)
				{
					k = integer_indices[i];
					if(fabs(x1[k] - x2[k]) >= .00001)
					{
						if(printing_in_setbranch) printf("x1_%d: %lf, x2_%d: %lf\n",k,x1[k],k,x2[k]);
						frac_index = k;
						frac_val = fmax(x1[k],x2[k])-.5;
						frac_values[k] = frac_val;
						if(frac_scores[i] > .0001)
				  		{
					  		frac_scores[i] += add_val2;
			/*		  		num_frac++;*/
					  	}
					  	else
					  	{
					  		frac_scores[i] += add_val2 + multiplier*k;
				  			num_frac++;
					  	}
					  	num_found_1_2_++;
/*						found_one = 1;*/
/*						break;*/
					}
				}
				if(doing_it_again_) goto OUTTA_HERE2;
			}
		}
		if(num_found_ws_1_ < num_found_ws_2_ && num_found_ws_1_ < num_found_1_2_)
		{
			add_val2++;
			doing_it_again_ = 1;
			goto AGAIN1_;
		}
		if(num_found_ws_2_ < num_found_ws_1_ && num_found_ws_2_ < num_found_1_2_)
		{
			add_val2++;
			doing_it_again_ = 1;
			goto AGAIN2_;
		}
		if(num_found_1_2_ < num_found_ws_1_ && num_found_1_2_ < num_found_ws_2_)
		{
			add_val2++;
			doing_it_again_ = 1;
			goto AGAIN3_;
		}
		
		OUTTA_HERE2:
		
		if(!desired_index_exists2 && num_frac)
		{
			if(printing_in_setbranch) printf("choosing branch variable based on lps\n");
			qsort(frac_scores, total_num_integer, sizeof(double), comparison_function3);
/*			for(i=0;i<total_num_integer;i++) printf("%lf\n",frac_scores[i]);*/
			int index = (int) (1./multiplier * (frac_scores[0] - (int) frac_scores[0])+.00001);
/*			if(seqnum == 4) index = 36;*/
			frac_index = index;
			frac_val = frac_values[frac_index];
			if(printing_in_setbranch) printf("index: %d, val: %lf\n",index,frac_val);
			found_one = 1;
		}
		else if(!desired_index_exists2 && !found_one) goto BRANCHASCPLEX2;
		
		if(frac_index != -1)
		{		
			if(printing_in_setbranch) printf("there is an index we want to use for branching: %d\n",frac_index);	
			double ub2[1] = {0.};
			double lb2[1] = {0.};
			status = CPXgetub (env, nodelp, ub2, frac_index, frac_index);
			if ( status ) {
				printf ("(%d) Failed to get ub's. Error code %d\n", __LINE__,status);
				printf("cur_numcols: %d\tfrac_index: %d\n",cur_numcols,frac_index);
				exit(0);
			}
			status = CPXgetlb (env, nodelp, lb2, frac_index, frac_index);
			if ( status ) {
				printf ("(%d) Failed to get lb's. Error code %d\n", __LINE__,status);
				printf("cur_numcols: %d\tfrac_index: %d\n",cur_numcols,frac_index);
				exit(0);
			}
			if(ub2[0] == lb2[0])
			{
				if(printing_in_setbranch) printf("upper and lower bounds matched for desired index, CPLEX default branch\n");
/*				printf("indices[0]: %d\n",indices[0]);*/
				goto BRANCHASCPLEX2;
			}
		
			if(!bound_reduction)
			{
				if(printing_in_setbranch) printf("search region size hasn't been reduced\n");
				if(bd_reduction_during_branching)
				{
					if(printing_in_setbranch) printf("attempting to reduce bounds\n");
					status = CPXgetub (env, nodelp, ub_, 0, cur_numcols-1);
					status = CPXgetlb (env, nodelp, lb_, 0, cur_numcols-1);
					ub_[obj1_index] = reduced_subprob_x_ub;
					ub_[obj2_index] = reduced_subprob_y_ub;
					lb_[obj1_index] = reduced_subprob_x_lb;
					lb_[obj2_index] = reduced_subprob_y_lb;
				
					k = 0;
					for(i=0;i<cur_numcols;i++)
					{
						if(i == integer_indices[k])
						{
							low_up[2*k] = 'L';
							low_up[2*k+1] = 'U';
							bds_br1[2*k] = lb_[i];
							bds_br2[2*k] = lb_[i];
							bds_br1[2*k+1] = ub_[i];
							bds_br2[2*k+1] = ub_[i];
							br_ind[2*k] = integer_indices[k];
							br_ind[2*k+1] = integer_indices[k];
/*							printf("x_%d: %lf to %lf \t %lf to %lf\n",i,bds_br1[2*k],bds_br1[2*k+1],bds_br2[2*k],bds_br2[2*k+1]);*/
							k++;
						}
						if(k >= total_num_integer) break;
					}
				
					bd_red_start_time = clock();
					if(printing_in_setbranch) printf("doing bound reduction\n");
					int reduced_worked = reduce_bds_before_branching( (CPXENVptr) env, nodelp, bds_br1, bds_br2, 0, frac_index, frac_val, 
											(ub_[obj2_index]-lb_[obj2_index])/(lb_[obj1_index]-ub_[obj1_index]));
					bd_red_finish_time = clock();
					time_tightening_variable_bounds += (double)(bd_red_finish_time - bd_red_start_time) / CLOCKS_PER_SEC;
					
/*					k = 0;*/
/*					for(i=0;i<cur_numcols;i++)*/
/*					{*/
/*						if(i == integer_indices[k])*/
/*						{*/
/*							printf("x_%d: %lf to %lf \t %lf to %lf\n",i,bds_br1[2*k],bds_br1[2*k+1],bds_br2[2*k],bds_br2[2*k+1]);*/
/*							k++;*/
/*						}*/
/*						if(k >= total_num_integer) break;*/
/*					}*/
									
					*useraction_p = CPX_CALLBACK_SET;
				
					if(reduced_worked)
					{
						if(printing_in_setbranch) printf("branching at line %d\n",__LINE__);
/*						printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
						if(printing_in_setbranch) printf("reduction worked\n");
						status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 2*total_num_integer, br_ind, 
										low_up, bds_br1, nodeest[1], userhandle_up, &seqnum1);
			   			if ( status )
			   			{  
/*			   				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
			   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
			   				goto TERMINATE;
			   			}
			   			status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 2*total_num_integer, br_ind, 
											low_up, bds_br2, nodeest[0], userhandle_down, &seqnum2);
			   			if ( status )
			   			{  
/*			   				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
			   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
			   				goto TERMINATE;
			   			}
			   			
			   			if(printing_in_setbranch) printf("successfully created branches: %d and %d\n",seqnum1,seqnum2);
				   		
			   			goto TERMINATE;
		   			}
		   			else
		   			{
		   				if(printing_in_setbranch) printf("it didn't work\n");
		   				goto BRANCHASCPLEX2;
		   			}
	   			}
				
/*				printf("changing branch variable to %d\n",frac_index);*/

				*useraction_p = CPX_CALLBACK_SET;

				vars[0] = frac_index;
				varlu[0] = 'U';
				varbd[0] = (int) frac_val;
	
				if(printing_in_setbranch) printf("branching at line %d\n",__LINE__);
/*				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
				status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 1, vars, 
									varlu, varbd, nodeest[1], userhandle_up, &seqnum1);
	   			if ( status )
	   			{  
/*	   				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
	   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
	   				goto TERMINATE;
	   			}
	   			
				varlu[0] = 'L';
				varbd[0] = varbd[0]+1;
				status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 1, vars, 
									varlu, varbd, nodeest[0], userhandle_down, &seqnum2);
	   			if ( status )
	   			{
	   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
	   				goto TERMINATE;
	   			}
		   		
				frac_index = -1;
				num_x_frac_at_both = 0;
				psa_int_pivot = 0;
			}
			else
			{
				if(bd_reduction_during_branching)
				{
					status = CPXgetub (env, nodelp, ub_, 0, cur_numcols-1);
					status = CPXgetlb (env, nodelp, lb_, 0, cur_numcols-1);
				
					if(printing_in_setbranch) printf("changing branch variable to %d\n",frac_index);
					int i,k = 0;
					for(i=0;i<cur_numcols;i++)
					{
						if(i == integer_indices[k])
						{
							low_up[2*k] = 'L';
							low_up[2*k+1] = 'U';
							bds_br1[2*k] = lb_[i];
							bds_br2[2*k] = lb_[i];
							bds_br1[2*k+1] = ub_[i];
							bds_br2[2*k+1] = ub_[i];
							br_ind[2*k] = integer_indices[k];
							br_ind[2*k+1] = integer_indices[k];
							k++;
						}
						if(k >= total_num_integer) break;
					}
				
					low_up[2*total_num_integer] = 'L';
					low_up[2*total_num_integer+1] = 'U';
					low_up[2*total_num_integer+2] = 'L';
					low_up[2*total_num_integer+3] = 'U';
					bds_br1[2*total_num_integer] = reduced_subprob_x_lb;
					bds_br1[2*total_num_integer+1] = reduced_subprob_x_ub;
					bds_br1[2*total_num_integer+2] = reduced_subprob_y_lb;
					bds_br1[2*total_num_integer+3] = reduced_subprob_y_ub;
					bds_br2[2*total_num_integer] = reduced_subprob_x_lb;
					bds_br2[2*total_num_integer+1] = reduced_subprob_x_ub;
					bds_br2[2*total_num_integer+2] = reduced_subprob_y_lb;
					bds_br2[2*total_num_integer+3] = reduced_subprob_y_ub;
					br_ind[2*total_num_integer] = obj1_index;
					br_ind[2*total_num_integer+1] = obj1_index;
					br_ind[2*total_num_integer+2] = obj2_index;
					br_ind[2*total_num_integer+3] = obj2_index;
				
					bd_red_start_time = clock();
					if(printing_in_setbranch) printf("doing bound reduction\n");
					
/*					printf("$$$$$$$$$$$$$$$$$$$\n");*/
/*					for(i=0;i<total_num_integer;i++) printf("x%d, bd1: %lf to %lf\t bd2: %lf to %lf\n",integer_indices[i],bds_br1[2*i],bds_br1[2*i+1],bds_br2[2*i],bds_br2[2*i+1]);*/
/*					printf("$$$$$$$$$$$$$$$$$$$\n");*/
					
					int reduced_worked = reduce_bds_before_branching( (CPXENVptr) env, nodelp, bds_br1, bds_br2, 0, frac_index, frac_val,
											(ub_[obj2_index]-lb_[obj2_index])/(lb_[obj1_index]-ub_[obj1_index]));
					bd_red_finish_time = clock();
					time_tightening_variable_bounds += (double)(bd_red_finish_time - bd_red_start_time) / CLOCKS_PER_SEC;
										
					*useraction_p = CPX_CALLBACK_SET;
				
					if(reduced_worked)
					{
/*						printf("$$$$$$$$$$$$$$$$$$$\n");*/
/*						for(i=0;i<total_num_integer;i++) printf("x%d, bd1: %lf to %lf\t bd2: %lf to %lf\n",integer_indices[i],bds_br1[2*i],bds_br1[2*i+1],bds_br2[2*i],bds_br2[2*i+1]);*/
/*						printf("$$$$$$$$$$$$$$$$$$$\n");*/
						if(printing_in_setbranch) printf("branching at line %d\n",__LINE__);
/*						printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
						status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 2*total_num_integer+4, br_ind, 
											low_up, bds_br1, nodeest[1], userhandle_up, &seqnum1);
			   			if ( status )
			   			{  
/*			   				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
			   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
			   				goto TERMINATE;
			   			}
			   			status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 2*total_num_integer+4, br_ind, 
											low_up, bds_br2, nodeest[0], userhandle_down, &seqnum2);
			   			if ( status )
			   			{  
/*			   				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
			   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
			   				goto TERMINATE;
			   			}
				   		goto TERMINATE;
		   			}
	   			}
	   			
/*				printf("changing branch variable to %d and also reducing subproblem size\n",frac_index);*/

				*useraction_p = CPX_CALLBACK_SET;

				if(printing_in_setbranch) printf("branching at line %d\n",__LINE__);
/*				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
				vars3[0] = frac_index;
				varlu3[0] = 'U';
				varbd3[0] = (int) frac_val;
/*				printf("frac val: %lf, low: %lf, up: %lf\n",frac_val,varbd3[0],varbd3[0]+1);*/
				vars3[1] = obj1_index;
				vars3[2] = obj1_index;
				vars3[3] = obj2_index;
				vars3[4] = obj2_index;
				varlu3[1] = 'L';
				varlu3[2] = 'U';
				varlu3[3] = 'L';
				varlu3[4] = 'U';
				varbd3[1] = reduced_subprob_x_lb;
				varbd3[2] = reduced_subprob_x_ub;
				varbd3[3] = reduced_subprob_y_lb;
				varbd3[4] = reduced_subprob_y_ub;
				
				status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 5, vars3, 
									varlu3, varbd3, nodeest[1], userhandle_up, &seqnum1);
	   			if ( status )
	   			{  
/*	   				printf("frac_index: %d, frac_Val: %lf\n",frac_index,frac_val);*/
	   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
	   				goto TERMINATE;
	   			}
	   			
				varlu3[0] = 'L';
				varbd3[0] = varbd3[0]+1;

				status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom, 5, vars3, 
									varlu3, varbd3, nodeest[0],userhandle_down, &seqnum2);
	   			if ( status )
	   			{
	   				printf("%s(%d): CPXbranchcallbackbranchbds, Failed to branch %d\n", __FILE__, __LINE__, status);
	   				goto TERMINATE;
	   			}
		   		
				frac_index = -1;
				num_x_frac_at_both = 0;
				psa_int_pivot = 0;
			}
		}
		else
		{
			BRANCHASCPLEX2:
			if(printing_in_setbranch) printf("(%d) using CPLEX default branch\n", __LINE__);
			*useraction_p = CPX_CALLBACK_SET;
			
			if(printing_in_setbranch) printf("branching at line %d\n",__LINE__);
			status = CPXbranchcallbackbranchasCPLEX(env, cbdata, wherefrom, 0, userhandle_up, &seqnum1);
			if ( status )
	   		{
	   			printf("%s(%d): CPXbranchcallbackbranchasCPLEX, Failed to branch %d\n", __FILE__, __LINE__, status);
	   			goto TERMINATE;
			}
			
			status = CPXbranchcallbackbranchasCPLEX(env, cbdata, wherefrom, 1, userhandle_down, &seqnum2);
			if ( status )
	   		{
	   			printf("%s(%d): CPXbranchcallbackbranchasCPLEX, Failed to branch %d\n", __FILE__, __LINE__, status);
	   			goto TERMINATE;
			}
		}
	}

 	TERMINATE:
 	
/* 	if(x_ws && (x_ws[obj1_index] != x_ws[obj1_index] || x_ws[obj2_index] != x_ws[obj2_index]))*/
/* 	{*/
/* 		printf("x_ws solns got f'ed up. Exitting\n");*/
/* 		exit(0);*/
/* 	}*/
 	
 	if(printing_in_setbranch) printf("Created branches: %d and %d\n",seqnum1,seqnum2);
 	
 	finish_time = clock();
 	
 	time_branching += (double)(finish_time - start_time) / CLOCKS_PER_SEC;
 	
/* 	printf("arrived at terminate\n");*/
 	
 	if(branch_iterations < 5 && fathoming != 1)
 	{
 		if(printing_in_setbranch) printf("starting user heuristic\n");
	 	i = 0;
	 	if(pure_binary == 1) i += binary_heuristic(env,lp_1,lp_2);
		else i += mixed_heuristic(env);
/*		printf("number of solns added from user heuristic: %d\n",i);*/
	}
	
	if(printing_in_setbranch) printf("fathoming value: %d\n",fathoming);
	
	if(fathoming == 1)
 	{
/* 		printf("going to fathom\n");*/
 		*useraction_p = CPX_CALLBACK_SET;
		nodecnt = 0;
	}
 
 	if(userhandle_current)
 	{
/* 		printf("checking whether or not x solns are stored\n");*/
 		if(userhandle_current->x_ws) free_and_null((char **) &userhandle_current->x_ws);
/* 		printf("freeing current x1\n");*/
 		if(userhandle_current->x1) free_and_null((char **) &userhandle_current->x1);
/* 		printf("x2?\n");*/
 		if(userhandle_current->x2) free_and_null((char **) &userhandle_current->x2);
 		if(userhandle_current->prob) CPXfreeprob(env, &(userhandle_current->prob));
 		free(userhandle_current);
/* 		printf("trying to set userhandle_current to NULL\n");*/
 		userhandle_current = NULL;
/* 		printf("the end\n");*/
 	}
 	if(fathoming == 1)
 	{
/* 		printf("fathoming, must free local userhandles\n");*/
 		if(userhandle_up)
	 	{
	 		if(userhandle_up->x_ws) free_and_null((char **) &userhandle_up->x_ws);
/*	 		printf("freeing up x1\n");*/
	 		if(userhandle_up->x1) free_and_null((char **) &userhandle_up->x1);
/*	 		printf("freeing up x2\n");*/
	 		if(userhandle_up->x2) free_and_null((char **) &userhandle_up->x2);
	 		if(userhandle_up->prob) CPXfreeprob(env, &(userhandle_up->prob));
	 		free_and_null((char **) &userhandle_up);
/*	 		printf("trying to set userhandle_up to NULL\n");*/
	 		userhandle_up = NULL;
	 	}
	 	if(userhandle_down)
	 	{
	 		if(userhandle_down->x_ws) free_and_null((char **) &userhandle_down->x_ws);
/*	 		printf("freeing down x1\n");*/
	 		if(userhandle_down->x1) free_and_null((char **) &userhandle_down->x1);
/*	 		printf("freeing down x2\n");*/
	 		if(userhandle_down->x2) free_and_null((char **) &userhandle_down->x2);
	 		if(userhandle_down->prob) CPXfreeprob(env, &(userhandle_down->prob));
	 		free_and_null((char **) &userhandle_down);
/*	 		printf("trying to set userhandle_down to NULL\n");*/
	 		userhandle_down = NULL;
	 	}
	 	fathoming = 0;
 	}
 	
/* 	SKIP_THIS11:*/
 	
/* 	printf("(%d) freeing\n",__LINE__);*/
 	if(x_ws != NULL) free_and_null((char **) &x_ws);
/* 	printf("(%d) freeing\n",__LINE__);*/
 	if(x1 != NULL) free_and_null((char **) &x1);
/* 	printf("(%d) freeing\n",__LINE__);*/
 	if(x2 != NULL) free_and_null((char **) &x2);
/* 	printf("(%d) freeing\n",__LINE__);*/
 	if(feas != NULL) free_and_null((char **) &feas);
/* 	printf("(%d) freeing\n",__LINE__);*/
 	if(low_up != NULL) free_and_null((char **) &low_up);
/* 	printf("(%d) freeing\n",__LINE__);*/
	if(lb_ != NULL) free_and_null((char **) &lb_);
/* 	printf("(%d) freeing\n",__LINE__);*/
	if(ub_ != NULL) free_and_null((char **) &ub_);
/* 	printf("(%d) freeing\n",__LINE__);*/
	if(bds_br1 != NULL) free_and_null((char **) &bds_br1);
/* 	printf("(%d) freeing\n",__LINE__);*/
	if(bds_br2 != NULL) free_and_null((char **) &bds_br2);
/* 	printf("(%d) freeing\n",__LINE__);*/
	if(br_ind != NULL) free_and_null((char **) &br_ind);
/* 	printf("(%d) freeing\n",__LINE__);*/
	if(cstat != NULL) free_and_null((char **) &cstat);
/* 	printf("(%d) freeing\n",__LINE__);*/
	if(rstat != NULL) free_and_null((char **) &rstat);
/* 	printf("(%d) freeing\n",__LINE__);*/
	if(cstat_ws != NULL) free_and_null((char **) &cstat_ws);
/* 	printf("(%d) freeing\n",__LINE__);*/
	if(rstat_ws != NULL) free_and_null((char **) &rstat_ws);
/*	printf("(%d) freeing\n",__LINE__);*/
	if(var_ubs != NULL) free_and_null((char **) &var_ubs);
/* 	printf("(%d) freeing\n",__LINE__);*/
	if(var_lbs != NULL) free_and_null((char **) &var_lbs);
/*	if(ub != NULL) free_and_null((char **) &ub);*/
/* 	printf("freeing\n");*/
/*	if(lb != NULL) free_and_null((char **) &lb);*/
/* 	printf("(%d) freeing\n",__LINE__);*/
	if(indexes != NULL) free_and_null((char **) &indexes);
/*	printf("done freeing\n");*/
/*	printf("freeing\n");*/
/*	if(frac_values != NULL) free_and_null((char **) &frac_values);*/
	
/* 	printf("freeing\n");*/
	CPXfreeprob(env_just_solve_mips, &nodelp_copy);
	CPXfreeprob(env_just_solve_mips, &nodelp_copy2);
/* 	printf("freeing\n");*/
	CPXfreeprob(env, &lp_ob1);
/* 	printf("freeing\n");*/
	CPXfreeprob(env, &lp_ob2);
/* 	printf("freeing\n");*/
	CPXfreeprob(env_just_solve_mips, &mip_ob1);
/* 	printf("freeing\n");*/
	CPXfreeprob(env_just_solve_mips, &mip_ob2);
/* 	printf("freeing last\n");*/
	CPXfreeprob(env, &nodelp2);

/*	remove_dual_bd(seqnum);*/

	//  status = CPXsetintparam ((CPXENVptr) env, CPX_PARAM_SCRIND, CPX_ON);
/*	printf("returning %d\n",status);*/
	return (status);

} /* END usersetbranch */

int num_new_cuts = 0;
int *rmatbeg = NULL;
int *rmatind = NULL;
double *rmatval = NULL;
double *rhs_s = NULL;
double *local_lb = NULL;
double *local_ub = NULL;

	/*************************************************************************************
	
		This branching callback is used while generating local cuts at a node
		it terminates the MIP after the root node is processed so as to not take
		too much time
		
	*************************************************************************************/

int CPXPUBLIC
usersetbranch2 (CPXCENVptr   env,
               void         *cbdata,
               int          wherefrom,
               void         *cbhandle,
               int          brtype,
               int          sos,
               int          nodecnt,
               int          bdcnt,
               const int    *nodebeg,
               const int    *indices,
               const char   *lu,
               const double *bd,
               const double *nodeest,
               int          *useraction_p)
{
/*  	printf("inside setbranch2\n");*/
  	int i;
  	int nzcnt;
  	int surplus,status;
  	num_new_cuts = 0;
  	
  	CPXLPptr nodelp2 = NULL;
	status = CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp2);
	if ( status ) {
		printf ("CPXgetcallbacknodelp, Failed to catch nodelp, error code %d\n", status);
		goto TERMINATE;
	}
	
	int numrows = CPXgetnumrows (env, nodelp2);
	if(numrows > sub_prob_numrows) 
	{
		num_new_cuts = numrows - sub_prob_numrows;
/*		printf("%d cuts have been added\n",num_new_cuts);*/
		rmatbeg = (int *) malloc ((num_new_cuts+1) * sizeof (int));
		rmatind = (int *) malloc (num_new_cuts*cur_numcols * sizeof (int));
		rmatval = (double *) malloc (num_new_cuts*cur_numcols * sizeof (double));
		rhs_s = (double *) malloc (num_new_cuts * sizeof (double));
		status = CPXgetrows (env, nodelp2, &nzcnt, rmatbeg, rmatind, rmatval, num_new_cuts*cur_numcols, &surplus, sub_prob_numrows, numrows-1);
		if(status) printf("(%d) there was an error\n",__LINE__);
		status = CPXgetrhs  (env, nodelp2, rhs_s, sub_prob_numrows, numrows-1);
		if(status) printf("(%d) there was an error\n",__LINE__);
		rmatbeg[num_new_cuts] = nzcnt;
	}
	if(!bd_reduction_during_branching)
	{
		local_lb = malloc (cur_numcols * sizeof (double));
		local_ub = malloc (cur_numcols * sizeof (double));
		int one = 1;
		status = CPXgetlb (env, nodelp2, local_lb, 0, cur_numcols-1);
		if ( status ) {
			printf ("(%d) Failed to get lb's. Error code %d\n", __LINE__,status);
			goto TERMINATE;
		}
		status = CPXgetub (env, nodelp2, local_ub, 0,  cur_numcols-1);
		if ( status ) {
			printf ("(%d) Failed to get ub's. Error code %d\n", __LINE__,status);
			goto TERMINATE;
		}
	}

	*useraction_p = CPX_CALLBACK_SET;
	nodecnt = 0;
	
 	TERMINATE:

	return (status);

} /* END usersetbranch2 */

int CPXPUBLIC
usersetbranch3 (CPXCENVptr   env,
               void         *cbdata,
               int          wherefrom,
               void         *cbhandle,
               int          brtype,
               int          sos,
               int          nodecnt,
               int          bdcnt,
               const int    *nodebeg,
               const int    *indices,
               const char   *lu,
               const double *bd,
               const double *nodeest,
               int          *useraction_p)
{
/*  	printf("inside setbranch3\n");*/
  	int i;
  	int nzcnt;
  	int surplus,status;
  	
  	CPXLPptr nodelp3 = NULL;
	status = CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp3);
	if ( status ) {
		printf ("CPXgetcallbacknodelp, Failed to catch nodelp, error code %d\n", status);
		goto TERMINATE;
	}
	status = CPXwriteprob (env, nodelp3, "myprob2.lp", "LP");
	
	lp1_get_cuts_copy = CPXcloneprob (env, nodelp3, &status);

	*useraction_p = CPX_CALLBACK_SET;
	nodecnt = 0;
	
 	TERMINATE:

	return (status);

} /* END usersetbranch3 */
  
  
  
	/*************************************************************************************
	
		This callback is used for processing integer feasible solutions
		found during BB. Basically, we utilize anything found by the 
		heuristic, but tell CPLEX to reject all incumbents so that BB does
		not terminate early due to single objective bounds meet.
		
	*************************************************************************************/

     /**************************************************************************************/
     /*	 _  _  ____  ____  ____    __  __ _   ___  _  _  _  _  ____  ____  __ _  ____      */
     /*	/ )( \/ ___)(  __)(  _ \  (  )(  ( \ / __)/ )( \( \/ )(  _ \(  __)(  ( \(_  _)     */
     /*	) \/ (\___ \ ) _)  )   /   )( /    /( (__ ) \/ (/ \/ \ ) _ ( ) _) /    /  )(       */
     /*	\____/(____/(____)(__\_)  (__)\_)__) \___)\____/\_)(_/(____/(____)\_)__) (__)      */
     /**************************************************************************************/
int seqnum_of_prev_pareto_branch = -1;

int CPXPUBLIC
userincumbent (CPXCENVptr env,
               void       *cbdata,
               int        wherefrom,
               void       *cbhandle,
               double     objval,
               double     *x,
               int        *isfeas_p,
               int        *useraction_p)
{
	int status = 0;
	int nodesleft;
	int seqnum = -1;
	
	CPXLPptr nodelp = NULL;
	
	int      *cstat = NULL;
  	int      *rstat = NULL;
  	int 	 *indices = NULL;
  	double	 *red_cost1 = NULL;
  	double	 *red_cost2 = NULL;
  	double	 *red_cost3 = NULL;
  	double	 *red_cost4 = NULL;
	
/*	printf("entering userincumbent. wherefrom is %d\n", wherefrom);*/
	
	/* Initialize useraction to indicate no user node selection */

	*useraction_p = CPX_CALLBACK_DEFAULT;
	
	/* Get sequence number of current node. */
	
	status = CPXgetcallbacknodeinfo(env,
                                 	cbdata,
                                 	wherefrom,
                                 	0,
                                 	CPX_CALLBACK_INFO_NODE_SEQNUM_LONG,
                                 	&seqnum);
/*	printf("The sequence number of a node with an integer feasible solution: %d\n", seqnum);*/
	
	if ((status = CPXgetcallbackinfo (env, cbdata, wherefrom, 
		CPX_CALLBACK_INFO_NODES_LEFT, &nodesleft)))
		goto TERMINATE;
		
	if(nodesleft <=0) *isfeas_p = 1; 	/* accept the incumbent */
	else *isfeas_p = 0; 			/*reject the incumbent and branch */

	// although we don't want Cplex to use this incumbent, we'll just
	// use its value and insert it in our data structure

	*useraction_p = CPX_CALLBACK_SET;
	
	int add_check = mock_insert(1,x[obj1_index],x[obj2_index],0,0,0,&tree);
/*	int add_check;*/

	if(add_check) 
	{
		PSA_full(env,NULL,x,NULL,NULL);
	}

	TERMINATE:
	
/*	printf("arrived at terminate\n");*/
	
	if(cstat)  free_and_null((char **) &cstat);
  	if(rstat)  free_and_null((char **) &rstat);
  	if(indices)    free_and_null((char **) &indices);
  	if(red_cost1)  free_and_null((char **) &red_cost1);
  	if(red_cost2)  free_and_null((char **) &red_cost2);
  	if(red_cost3)  free_and_null((char **) &red_cost3);
  	if(red_cost4)  free_and_null((char **) &red_cost4);


	return (status);

} /* END userincumbent */

int CPXPUBLIC
userincumbent2 (CPXCENVptr env3,
               void       *cbdata,
               int        wherefrom,
               void       *cbhandle,
               double     objval,
               double     *x,
               int        *isfeas_p,
               int        *useraction_p)
{

	*useraction_p = CPX_CALLBACK_DEFAULT;
	int status = 0;
	
	CPXLPptr nodelp = NULL;
	status = CPXgetcallbacknodelp (env3, cbdata, wherefrom, &nodelp);
	if ( status ) {
		printf ("CPXgetcallbacknodelp, Failed to catch nodelp, error code %d\n", status);
	}
	
	int add_check = mock_insert(1,x[obj1_index],x[obj2_index],0,0,0,&tree);
	int i;

	if(add_check) 
	{
    		for(i=0;i<cur_numcols;i++)
		{
		      	stored_x[x_rotation][i] = x[i];
		/*      printf("stored val: %lf\n",stored_x[x_rotation][i]);*/
		}
		x_rotation = (x_rotation + 1) % num_x_to_store;
		PSA_full(env_global,NULL,x,NULL,NULL);
	}
	
	return (status);

} /* END userincumbent2 */

int CPXPUBLIC
userincumbent_presolve (CPXCENVptr env,
               void       *cbdata,
               int        wherefrom,
               void       *cbhandle,
               double     objval,
               double     *x,
               int        *isfeas_p,
               int        *useraction_p)
{

	*useraction_p = CPX_CALLBACK_DEFAULT;
	int status = 0;
	
	CPXLPptr nodelp = NULL;
	status = CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp);
	if ( status ) {
		printf ("CPXgetcallbacknodelp, Failed to catch nodelp, error code %d\n", status);
	}
	
	
	int add_check = mock_insert(1,x[obj1_index],x[obj2_index],0,0,0,&tree);
	int i;

	if(add_check) 
	{
    		for(i=0;i<cur_numcols-2;i++)
		{
		      	stored_x[x_rotation][i] = x[i];
		}
		x_rotation = (x_rotation + 1) % num_x_to_store;
		PSA_full(env,NULL,x,NULL,NULL);
	}
	
	return (status);

} /* END userincumbent_presolve */

int prev_snum = -1;
int found_local_cuts = 0;
/*int cut_callback_count = 0;*/
int cutcallback_iteration = 0;
/*int solve_extra_mips_during_cut_callback = 0;*/

	/*************************************************************************************
	
		This callback is used for two purposes. First, it is where locally
		valid cutting planes are generated. Second, a weighted sum single
		objective MIP is solved here which allows us to generate another
		local cut along the level curve associated with the best known dual
		bound at termination. Essentially, here we attempt to tighten the
		dual bound at a node.
		
	*************************************************************************************/

     /**************************************************************************************/
     /*		  ___  _  _  ____     ___   __   __    __    ____   __    ___  __ _ 	   */
     /*		 / __)/ )( \(_  _)   / __) / _\ (  )  (  )  (  _ \ / _\  / __)(  / )	   */
     /*		( (__ ) \/ (  )(    ( (__ /    \/ (_/\/ (_/\ ) _ (/    \( (__  )  ( 	   */	
     /*		 \___)\____/ (__)    \___)\_/\_/\____/\____/(____/\_/\_/ \___)(__\_)       */
     /**************************************************************************************/

int CPXPUBLIC
 cutcallback (CPXCENVptr env,
           void *cbdata,
           int wherefrom,
           void *cbhandle,
           int *useraction_p)
{
/*	if(printing_in_setbranch) */
/*	printf("in the cut callback");*/

	start_time = clock();
	cumulative_time = (double)(start_time - start_BB) / CLOCKS_PER_SEC;

	int status = 0;
	int seqnum = -1, surplus = 0;
	num_new_cuts = 0;
	int i;
	int *indices = NULL;
	int add_check = 0, iteration = 0;
	double objvals[2] = {0.,0.};
	double original_slope = 0.;
	
	CPXLPptr nodelp2 = NULL, nodelp_copy = NULL, nodelp_copy2 = NULL, lp_ob1 = NULL, lp_ob2 = NULL, mip_ob1 = NULL, mip_ob2 = NULL, cut_prob = NULL;
/*	double *x_ws = NULL;*/
/*   	double *x1 = NULL;*/
/*   	double *x2 = NULL;*/
	int *feas = NULL;
	char *low_up = NULL;
	double *lb_ = NULL;
	double *ub_ = NULL;
	double *bds_br1 = NULL;
	double *bds_br2 = NULL;
	int *br_ind = NULL;
	int *indexes = NULL;
	int *cstat = NULL;
	int *rstat = NULL;
	int *cstat_ws = NULL;
	int *rstat_ws = NULL;
/*	user_data *userhandle_up = NULL;*/
/*	user_data *userhandle_down = NULL;*/
	double *var_ubs = NULL;
	double *var_lbs = NULL;
	
	cutcallback_iteration = 0;
	
	CPXLPptr nodelp = NULL, nodelp_mip = NULL;
	CPXENVptr env4;
	    
	env4 = CPXopenCPLEX(&status);
	if(status)
	{
		printf("failed to open another copy of CPLEX\n");
		exit(0);
	}
	
	if(its_been_only_points)
	{
/*		printf("%d\n",__LINE__);*/
		*useraction_p = CPX_CALLBACK_DEFAULT;
		goto TERMINATE;
	}
	
	status = CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp);
	if ( status ) {
		printf ("CPXgetcallbacknodelp, Failed to catch nodelp, error code %d\n", status);
		goto TERMINATE;
	}
	
	status = CPXgetcallbacknodeinfo(env,
                                 	cbdata,
                                 	wherefrom,
                                 	0,
                                 	CPX_CALLBACK_INFO_NODE_SEQNUM_LONG,
                                 	&seqnum);
        if ( status ) {
		printf ("CPXgetcallbacknodeinfo, Failed to get seqnum, error code %d\n", status);
		goto TERMINATE;
	}
	
	status = CPXgetcallbacknodex (env, cbdata, wherefrom, objvals, obj1_index, obj2_index);
                               
/*	if(printing_in_setbranch) */
/*	printf(" and the seqnum is: %d\n", seqnum);*/
	
/*	printf("__________________________________________________\n");*/
/*	PSA_all(env,nodelp);*/
/*	printf("__________________________________________________\n");*/
	
	sub_prob_numrows = CPXgetnumrows (env, nodelp);
	double ran = (double) rand() / ( (double) RAND_MAX);
	
	double lbs[2] = {0.,0.};
	double ubs[2] = {0.,0.};
	
	if(seqnum != prev_snum) // && ran < .05)
	{
		prev_snum = seqnum;
		last_cutcallback_seqnum = seqnum;
		
/*		user_data *userhandle_current = NULL;*/
		status = CPXgetcallbacknodeinfo(	env,
				         		cbdata,
				         		wherefrom,
				         		0,
				        		CPX_CALLBACK_INFO_NODE_USERHANDLE,
				         		&userhandle_current);
		double *x   = (double *) malloc (cur_numcols * sizeof (double));
		
		int num_starts = CPXgetnummipstarts (env_just_solve_mips, global_mip);
		int prev_numsols = 0, numsolns = 0, numrep = 0;
				         		
		var_ubs   = (double *) malloc (cur_numcols * sizeof (double));
		var_lbs   = (double *) malloc (cur_numcols * sizeof (double));
		
		status = CPXgetlb (env, nodelp, var_lbs, 0, cur_numcols-1);
		if ( status ) {
			printf ("(%d) Failed to get lb's for objectives. Error code %d\n",__LINE__, status);
			exit(0);
			goto TERMINATE;
		}

		status = CPXgetub (env, nodelp, var_ubs, 0, cur_numcols-1);
		if ( status ) {
			printf ("(%d) Failed to get ub's for objectives. Error code %d\n",__LINE__, status);
			goto TERMINATE;
		}		
				         		
		int ws_feas = 1, ob1_feas = 1, ob2_feas = 1;
		if(seqnum != 0 && userhandle_current)
		{
			int k = 0;
			if(!(userhandle_current->x_ws)) 
			{	
/*				if(printing_in_setbranch) printf("uh_cur x_ws didn't exist\n");*/
				ws_feas = 0;
			}
			if(!(userhandle_current->x1))
			{
/*				if(printing_in_setbranch) printf("uh_cur x1 didn't exist\n");*/
				ob1_feas = 0;
			}
			if(!(userhandle_current->x2))
			{
/*				if(printing_in_setbranch) printf("uh_cur x2 didn't exist\n");*/
				ob2_feas = 0;
			}
			while( (ws_feas || ob1_feas || ob2_feas) && k < cur_numcols )
			{
				
				if(ws_feas)
				{
/*					printf("k: %d, lb: %lf, x_ws: %lf, ub: %lf\n",k,var_lbs[k],userhandle_current->x_ws[k],var_ubs[k]);*/
					if(userhandle_current->x_ws[k] != userhandle_current->x_ws[k] || userhandle_current->x_ws[k] - var_lbs[k] < -.00001 
						|| userhandle_current->x_ws[k] - var_ubs[k] > .00001) ws_feas = 0;
				}
				if(ob1_feas)
				{
/*					printf("k: %d, lb: %lf, x1: %lf, ub: %lf\n",k,var_lbs[k],userhandle_current->x1[k],var_ubs[k]);*/
					if(userhandle_current->x1[k] != userhandle_current->x1[k] || userhandle_current->x1[k] - var_lbs[k] < -.00001 
						|| userhandle_current->x1[k] - var_ubs[k] > .00001) ob1_feas = 0;
				}
				if(ob2_feas)
				{
/*					printf("k: %d, lb: %lf, x2: %lf, ub: %lf\n",k,var_lbs[k],userhandle_current->x2[k],var_ubs[k]);*/
					if(userhandle_current->x2[k] != userhandle_current->x2[k] || userhandle_current->x2[k] - var_lbs[k] < -.00001 
						|| userhandle_current->x2[k] - var_ubs[k] > .00001) ob2_feas = 0;
				}
				k++;
			}
			if(ws_feas)
			{
				if(printing_in_setbranch) printf("previous ws mip soln still feasible here\n");
				userhandle_current->ws_still_feas = 1;
			}
			if(ob1_feas)
			{
				if(printing_in_setbranch) printf("previous ob1 mip soln still feasible here\n");
				userhandle_current->ob1_still_feas = 1;
			}
			if(ob2_feas)
			{
				if(printing_in_setbranch) printf("previous ob2 mip soln still feasible here\n");
				userhandle_current->ob2_still_feas = 1;
			}
		}
				         		
				         		
		if(!generate_local_cuts) goto SOLVE_THE_MIP;
		
		add_check = mock_insert(1,objvals[0],objvals[1],0,0,0,&tree);
		if((!userhandle_current || !(userhandle_current->ws_still_feas)) && !add_check)
		{
/*			printf("original add check unsuccesful, checking with modified slope\n");*/
			indices   = (int *) malloc (cur_numcols * sizeof (int));
			for(i=0;i<cur_numcols;i++) indices[i] = i;
			
			CPXLPptr nodelp3 = CPXcloneprob (env, nodelp, &status);
	
			status = CPXgetlb (env, nodelp3, lbs, obj1_index, obj2_index);
			if ( status ) {
				printf ("(%d) Failed to get lb's for objectives. Error code %d\n",__LINE__, status);
				exit(0);
				goto TERMINATE;
			}
	
			status = CPXgetub (env, nodelp3, ubs, obj1_index, obj2_index);
			if ( status ) {
				printf ("(%d) Failed to get ub's for objectives. Error code %d\n",__LINE__, status);
				goto TERMINATE;
			}
	
			double slope = (ubs[1]-lbs[1])/(lbs[0]-ubs[0]);
	
			if(slope != slope) //slope > -.001 || slope < -1000. || slope != slope) 
			{
				if(printing_in_setbranch) printf("setting slope to -1\n");
				slope = -1.;
			}
	
			chg_coefs(env,nodelp3,indices,slope);
			status = CPXlpopt (env, nodelp3);
			
			status = CPXgetx (env, nodelp3,  objvals, obj1_index, obj2_index);
  			if(status)
  			{
				printf("Failed to get x-values from CPLEX. Status: %d Line: %d\n",status,__LINE__);
				exit(0);
  			}
			add_check = mock_insert(1,objvals[0],objvals[1],0,0,0,&tree);
			CPXfreeprob(env, &nodelp3);
/*			printf("(%d) freeing indices\n",__LINE__);*/
			if(indices)    free_and_null((char **) &indices);
		}
		if( (userhandle_current && userhandle_current->ws_still_feas) || !add_check) goto SOLVE_THE_MIP;
		if((!userhandle_current || !(userhandle_current->ws_still_feas)) && add_check)
		{
			indices   = (int *) malloc (cur_numcols * sizeof (int));
			for(i=0;i<cur_numcols;i++) indices[i] = i;
			
			nodelp_mip = CPXcloneprob (env4, nodelp, &status);
	
			status = CPXgetlb (env4, nodelp_mip, lbs, obj1_index, obj2_index);
			if ( status ) {
				printf ("(%d) Failed to get lb's for objectives. Error code %d\n",__LINE__, status);
				exit(0);
				goto TERMINATE;
			}
	
			status = CPXgetub (env4, nodelp_mip, ubs, obj1_index, obj2_index);
			if ( status ) {
				printf ("(%d) Failed to get ub's for objectives. Error code %d\n",__LINE__, status);
				goto TERMINATE;
			}			
	
			double slope = (ubs[1]-lbs[1])/(lbs[0]-ubs[0]);
	
			if(slope > -.001 || slope < -1000. || slope != slope) slope = -1.;
			original_slope = slope;
	
			CPXchgprobtype(env4, nodelp_mip, CPXPROB_MILP);
			status = CPXchgctype(env4, nodelp_mip, cur_numcols, indices, xctype);
			if(status) printf("(%d) there was an error\n",__LINE__);
			status = CPXsetintparam (env4, CPX_PARAM_REPEATPRESOLVE, 0);
			if(status) printf("(%d) there was an error\n",__LINE__);
		 	status = CPXsetintparam (env4, CPX_PARAM_PREIND, CPX_OFF);
		 	if(status) printf("(%d) there was an error\n",__LINE__);
/*		  	status = CPXsetintparam (env4, CPX_PARAM_ADVIND, 1);*/
		  	status = CPXsetintparam (env4, CPX_PARAM_ADVIND, 0);
		  	if(status) printf("(%d) there was an error\n",__LINE__);
		  	status = CPXsetintparam (env4, CPX_PARAM_REDUCE, 0);
		  	if(status) printf("(%d) there was an error\n",__LINE__);
		  	status = CPXsetintparam (env4, CPX_PARAM_HEURFREQ, -1);
		  	if(status) printf("(%d) there was an error\n",__LINE__);
			
			status = CPXsetintparam ((CPXENVptr) env4, CPX_PARAM_SCRIND, CPX_OFF);
			status = CPXsetintparam ((CPXENVptr) env4, CPX_PARAM_THREADS, 1);
	    		if ( status ) {
	    			printf ("Failure to set threads to 1, error %d.\n",status);
	    			goto TERMINATE;
	    		}
	    		
	    		status = CPXsetbranchcallbackfunc(env4, usersetbranch2,  NULL) || CPXsetincumbentcallbackfunc(env4, userincumbent2, NULL); 
			if(status) printf("(%d) there was an error\n",__LINE__);
			
			status = CPXsetdblparam (env4, CPX_PARAM_TILIM, time_limit);
			if ( status ) {
		    		printf ("Failed to set solution time limit to time limit, error %d.\n",status);
			   	goto TERMINATE;
			}
			
			START_OVER:
			
			chg_coefs(env4,nodelp_mip,indices,slope);
			int nzcnt = 0;
			prev_numsols = 0;
			
			num_starts = CPXgetnummipstarts (env_just_solve_mips, global_mip);
				    	
		    	if(num_starts > global_num_starts)
		    	{
/*		    		printf("reallocation: (%d)\n",__LINE__);*/
		    		global_num_starts = num_starts;
		    		global_startspace = cur_numcols*global_num_starts;
		    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
				global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
				global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
				global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
		    	}

			status = CPXgetmipstarts (env_just_solve_mips, global_mip, &nzcnt, global_beg, global_varindices, 
					   global_values, global_effortlevel, global_startspace,
					   &surplus, 0, num_starts-1);
					   
			status = CPXaddmipstarts (env4, nodelp_mip, num_starts, nzcnt, global_beg, global_varindices,
					   global_values, global_effortlevel, NULL);

			CPXmipopt (env4, nodelp_mip);

			numsolns = CPXgetsolnpoolnumsolns (env4, nodelp_mip);
			numrep = CPXgetsolnpoolnumreplaced (env4, nodelp_mip);

			num_starts = numsolns - prev_numsols + numrep;
/*			printf("number of mip starts to use: %d\n",num_starts);*/
	
			prev_numsols = numsolns;
	
		    	if(num_starts > global_num_starts)
		    	{
/*		    		printf("reallocation: (%d)\n",__LINE__);*/
		    		global_num_starts = num_starts;
		    		global_startspace = cur_numcols*global_num_starts;
		    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
				global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
				global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
				global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
		    	}

			status = CPXgetmipstarts (env4, nodelp_mip, &nzcnt, global_beg, global_varindices, 
					   global_values, global_effortlevel, global_startspace,
					   &surplus, 0, num_starts-1);
					   
			status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,
					   global_values, global_effortlevel, NULL);
					   
		  	if(num_new_cuts) 
		  	{
/*			  	printf("%d cuts were added when examining seqnum %d (iteration %d)\n",num_new_cuts, seqnum, iteration);*/
				if(iteration == 0) found_local_cuts++;
			  	for(i=0;i<num_new_cuts;i++) 
			  	{
/*			  		printf("index: %d, rmatbeg of index: %d\n",i,rmatbeg[i]);*/
			  		status = CPXcutcallbackaddlocal   (env,
									   cbdata,
						       			   wherefrom,
									   rmatbeg[i+1]-rmatbeg[i],
									   rhs_s[i],
									   'L',
									   &rmatind[rmatbeg[i]],
									   &rmatval[rmatbeg[i]]);
					if(status) printf("(%d) there was an error\n",__LINE__);
				}
			}
			if(!bd_reduction_during_branching && local_lb)
			{
				double *lowbd = malloc (cur_numcols * sizeof (double));
				double *upbd = malloc (cur_numcols * sizeof (double));
				double one = 1.;
				status = CPXgetlb (env, nodelp, lowbd, 0, cur_numcols-1);
				if ( status ) {
					printf ("(%d) Failed to get lb's. Error code %d\n", __LINE__,status);
					goto TERMINATE;
				}
				status = CPXgetub (env, nodelp, upbd, 0, cur_numcols-1);
				if ( status ) {
					printf ("(%d) Failed to get ub's. Error code %d\n", __LINE__,status);
					goto TERMINATE;
				}
				
				for(i=0;i<cur_numcols;i++)
				{
					if(local_lb[i] > lowbd[i]) 
					{
						status = CPXcutcallbackaddlocal   (env,
									   cbdata,
						       			   wherefrom,
									   1,
									   local_lb[i],
									   'G',
									   &i,
									   &one);
						if(status) printf("(%d) there was an error\n",__LINE__);
					}
					if(local_ub[i] < upbd[i]) 
					{
						status = CPXcutcallbackaddlocal   (env,
									   cbdata,
						       			   wherefrom,
									   1,
									   local_ub[i],
									   'L',
									   &i,
									   &one);
						if(status) printf("(%d) there was an error\n",__LINE__);
					}
				}
				free(lowbd);
				free(upbd);
			}
			
			finish_time = clock();
			double tim = (double)(finish_time - start_time) / CLOCKS_PER_SEC;
			time_generating_cuts += tim;

			SOLVE_THE_MIP:
			
			if(printing_in_setbranch) printf("at solve_the_mip\n");
			
			start_time = clock();
			
			if(!indices) 
			{
				indices   = (int *) malloc (cur_numcols * sizeof (int));
				for(i=0;i<cur_numcols;i++) indices[i] = i;
			}
			
			if(!x_ws) x_ws = (double *) malloc (cur_numcols*sizeof(double));
   			if(!x1) x1 = (double *) malloc (cur_numcols*sizeof(double));
   			if(!x2) x2 = (double *) malloc (cur_numcols*sizeof(double));
			
/*			CPXLPptr nodelp2 = NULL, nodelp_copy = NULL, lp_ob1 = NULL, lp_ob2 = NULL, mip_ob1 = NULL, mip_ob2 = NULL;*/
			int lpstat = 0, depth = 0;
			int      i = -1, j = -1, k = -1, bestj = -1;
		   	double   maxinf = -CPX_INFBOUND;
		   	double xj_inf = 0., best_bound = 0., projection = 0., diff = 0., objval = 0.;
		   	clock_t start_mipsolve, finish_mipsolve;
		   	
/*		   	status = CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp);*/
/*			if ( status ) {*/
/*				printf ("CPXgetcallbacknodelp, Failed to catch nodelp, error code %d\n", status);*/
/*				goto TERMINATE;*/
/*			}*/
		   	
		   	cur_numcols = CPXgetnumcols (env, nodelp);
			cur_numrows = CPXgetnumrows (env, nodelp); 
		   	
		   	frac_index = -1;
		   	frac_val = 0.;
		   	
		   	double *x = NULL;
/*		   	int *feas = NULL;*/
		 
		   	char     varlu[1];
		   	double   varbd[1];
		   	int 	 vars[1];
		   	char     varlu2[2];
		   	double   varbd2[2];
		   	int	 vars2[2];
		   	char     varlu3[5];
		   	double   varbd3[5];
		   	int 	 vars3[5];
		   	seqnum1 += 2;
		   	seqnum2 += 2;
		   	
/*		   	int seqnum = -1;*/
		   	int add_check = 0, all_feas = 0, surplus = 0, infeas = 0, ws_lp_dom = 0, ob1_lp_dom = 0, ob2_lp_dom = 0;
		   	ws_mip_opt = 0;
		   	ob1_mip_opt = 0;
		   	ob2_mip_opt = 0;
		   	mip_solved = 0;
		   	bound_reduction = 0;
		   	left_side_dom = 0;
		   	right_side_dom = 0;
		   	cut_problem_created = 0;
		   	int ws_mip_ran = 0, ob1_mip_ran = 0, ob2_mip_ran = 0, ob1_lp_been_solved = 0, ob2_lp_been_solved = 0, ws_lp_int_feas = 0;
		   	int ob1_lp_int_feas = 0, ob2_lp_int_feas = 0, right_pt_dom = 0, right_seg_dom = 0;
		   	int ws_still_feas = 1, ob1_still_feas = 1, ob2_still_feas = 1, going_back = 0, ws_sol_interior = 0;
		   	int mipstarts_copied = 0;
		   	int from_pareto = 0, bound_reduced_from_PSA_reduce_left = 0, going_back_for_lp1 = 0, going_back_for_lp2 = 0, ws_lp_sol_on_top_bd = 0;
		   	int ws_lp_sol_on_right_bd = 0;
		   	
		   	double endpoint1_x = 0., endpoint2_x = 0., endpoint1_y = 0., endpoint2_y = 0.;
		   	
		   	int PSA_right_check = 0, PSA_left_check = 0;
		   	
/*		   	CPXLPptr cut_prob = NULL;*/
		   	
/*		   	double *x_ws = NULL;*/
/*		   	double *x1 = NULL;*/
/*		   	double *x2 = NULL;*/
/*			feas = NULL;*/
/*			char *low_up = NULL;*/
/*			double *lb_ = NULL;*/
/*			double *ub_ = NULL;*/
/*			double *bds_br1 = NULL;*/
/*			double *bds_br2 = NULL;*/
/*			int *br_ind = NULL;*/
/*			int *indexes = NULL;*/
/*			int *cstat = NULL;*/
/*			int *rstat = NULL;*/
/*			int *cstat_ws = NULL;*/
/*			int *rstat_ws = NULL;*/
/*			user_data *userhandle_up = NULL;*/
/*			user_data *userhandle_down = NULL;*/
	
			/*************** Exit early if too much time has been used *************************/
	
/*			if(cumulative_time > max_time || branch_iterations > max_nodes )*/
/*			{*/
/*				if(!printed_yet)*/
/*				{*/
/*					printed_yet = 1;*/
/*					int nodesleft = 0;*/
/*			*/
/*					duration_BB = cumulative_time;*/
/*			*/
/*					if ((status = CPXgetcallbackinfo (env, cbdata, wherefrom, */
/*					CPX_CALLBACK_INFO_NODES_LEFT, &nodesleft))) goto TERMINATE;*/
/*			*/
/*					if(cumulative_time > max_time) printf("Branch and Bound time has exceeded limit of %lf! Terminating after %lf seconds.\n",*/
/*						max_time,cumulative_time); */
/*					else printf("Branch and Bound has processed the maximum number of nodes, %d! Terminating after %lf seconds.\n",max_nodes,cumulative_time);*/
/*					printf("There are currently %d open nodes left, closing them and calculating global dual bound.\n",nodesleft); */
/*				}*/
/*				CPXLPptr nodelp = NULL;*/
/*				status = CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp);*/
/*				if ( status ) {*/
/*					printf ("CPXgetcallbacknodelp, Failed to catch nodelp, error code %d\n", status);*/
/*					goto TERMINATE;*/
/*				}*/
/*				build_dual_bd(env,nodelp);*/
/*		*/
/*				fathoming = 1;*/
/*				*useraction_p = CPX_CALLBACK_SET;*/
/*				nodecnt = 0;*/
/*		  		goto TERMINATE;	*/
/*			}*/
		   	
/*		   	x_ws = (double *) malloc (cur_numcols*sizeof(double));*/
/*		   	x1 = (double *) malloc (cur_numcols*sizeof(double));*/
/*		   	x2 = (double *) malloc (cur_numcols*sizeof(double));*/
			feas = (int *) malloc (cur_numcols*sizeof(int));
			low_up = (char *) malloc ((2*total_num_integer+16)*sizeof(char));
			lb_ = (double *) malloc ((cur_numcols)*sizeof(double));
			ub_ = (double *) malloc ((cur_numcols)*sizeof(double));
			bds_br1 = (double *) malloc ((2*total_num_integer+16)*sizeof(double));
			bds_br2 = (double *) malloc ((2*total_num_integer+16)*sizeof(double));
			br_ind = (int *) malloc ((2*total_num_integer+16)*sizeof(int));
			indexes = (int *) malloc (cur_numcols*sizeof(int));
			cstat = (int *) malloc (cur_numcols*sizeof(int));
			rstat = (int *) malloc (cur_numrows*sizeof(int));
			cstat_ws = (int *) malloc (cur_numcols*sizeof(int));
			rstat_ws = (int *) malloc (cur_numrows*sizeof(int));
			for(i=0;i<cur_numcols;i++) indexes[i] = i;
	
/*			printf("(%d) mallocing uh's up & down\n",__LINE__);*/
			userhandle_up = (user_data*) malloc( sizeof( user_data ) );
			userhandle_up->x_ws = NULL;
			userhandle_up->x1 = NULL;
			userhandle_up->x2 = NULL;
			userhandle_up->ws_still_feas = 0;
			userhandle_up->ob1_still_feas = 0;
			userhandle_up->ob2_still_feas = 0;
			userhandle_down = (user_data*) malloc( sizeof( user_data ) );
			userhandle_down->x_ws = NULL;
			userhandle_down->x1 = NULL;
			userhandle_down->x2 = NULL;
			userhandle_down->ws_still_feas = 0;
			userhandle_down->ob1_still_feas = 0;
			userhandle_down->ob2_still_feas = 0;
			
			userhandle_up->prob = NULL;
			userhandle_down->prob = NULL;
			
			if(show_progress)
			{
				userhandle_up->prob = CPXcloneprob (env, nodelp, &status);
				if ( status ) {
					printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
					goto TERMINATE;
				}
				userhandle_down->prob = CPXcloneprob (env, nodelp, &status);
				if ( status ) {
					printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
					goto TERMINATE;
				}
			}
	
			double ws_lp_objvals[2] = {0.,0.};
			double ob1_lp_objvals[2] = {0.,0.};
			double ob2_lp_objvals[2] = {0.,0.};
			double ws_mip_objvals[2] = {0.,0.};
			double ob1_mip_objvals[2] = {0.,0.};
			double ob2_mip_objvals[2] = {0.,0.};
	
		   	for(i=0;i<total_num_integer;i++) frac_scores[i] = 0.;
		   	num_frac = 0;
		   	multiplier = convert_it(cur_numcols);
		   	
		   	sub_pr1_x_ub = 0.;
			sub_pr1_y_lb = 0.;
			sub_pr1_x_lb = 0.;
			sub_pr1_y_ub = 0.;
			sub_pr2_x_lb = 0.;
			sub_pr2_y_ub = 0.;
			sub_pr2_x_ub = 0.;
			sub_pr2_y_lb = 0.;
	
			within_PSA_score = 3.;
	
			double lbs[2] = {0.,0.};
			double ubs[2] = {0.,0.};
			
			if(printing_in_setbranch) printf("About to process node %d\n",seqnum);
			nodelp2 = CPXcloneprob (env, nodelp, &status);
			if ( status ) {
				printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
				goto TERMINATE;
			}
	
			BEGINNING:
	
/*			if(infeasible_in_cutcallback)*/
/*			{*/
/*				fathoming = 1;*/
/*				infeasible_in_cutcallback = 0;*/
/*				*useraction_p = CPX_CALLBACK_SET;*/
/*				nodecnt = 0;*/
/*		  		goto TERMINATE;*/
/*			}*/

			if(cumulative_time > max_time || branch_iterations > max_nodes || (show_progress && break_early))
			{
				show_progress = 0;
				if(cumulative_time > max_time + max_time_build_dual_bound && !approximate_dual_bd)
				{
					approximate_dual_bd = 1;
					printf("Maximimum time allowed for generating dual bound has been exceeded. Approximating bound from here on.\n");
				}
				if(cumulative_time > max_time + 2*max_time_build_dual_bound)
				{
					if(!quit_generating_dual_bd)
					{
						quit_generating_dual_bd = 1;
						printf("Double the maximimum time allowed for generating dual bound has been exceeded. Exitting!!\nWARNING: This indicates that the dual bound reported will most likely not be accurate and thus the quality of the primal solutions is UNKNOWN!!!");
					}
					fathoming = 1;
			  		goto TERMINATE;	
				}
				if(!printed_yet)
				{
					printed_yet = 1;
					int nodesleft = 0;
			
					duration_BB = cumulative_time;
			
					if ((status = CPXgetcallbackinfo (env, cbdata, wherefrom, 
					CPX_CALLBACK_INFO_NODES_LEFT, &nodesleft))) goto TERMINATE;
			
					if(cumulative_time > max_time) printf("Branch and Bound time has exceeded limit of %lf! Terminating after %lf seconds.\n",
						max_time,cumulative_time); 
					else if(branch_iterations > max_nodes) printf("Branch and Bound has processed the maximum number of nodes, %d! Terminating after %lf seconds.\n",max_nodes,cumulative_time);
					else if(show_progress && break_early) printf("Duality gap is below desired threshold of %lf! Terminating after %lf seconds.\n",
						duality_gap_limit,cumulative_time);
					printf("There are currently %d open nodes left, closing them and calculating global dual bound.\n",nodesleft);
					fprintf(bb_results,"%lf\t",cumulative_time);
					fprintf(bb_results,"%d\t",nodesleft);
					fclose(bb_results);
					bb_results = fopen ("bb_results.txt", "a+");
					
					time_per = max_time_build_dual_bound/nodesleft; 
				}
				CPXLPptr nodelp = NULL;
				status = CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp);
				if ( status ) {
					printf ("CPXgetcallbacknodelp, Failed to catch nodelp, error code %d\n", status);
					goto TERMINATE;
				}
				printf("(%d) about to call build_dual_bd\n",__LINE__);
				build_dual_bd(env,nodelp);
		
				fathoming = 1;
		  		goto TERMINATE;	
			}
	
			status = CPXgetlb (env, nodelp2, lbs, obj1_index, obj2_index);
			if ( status ) {
				printf ("Failed to get lb's for objectives. Error code %d\n", status);
				goto TERMINATE;
			}
	
			status = CPXgetub (env, nodelp2, ubs, obj1_index, obj2_index);
			if ( status ) {
				printf ("Failed to get ub's for objectives. Error code %d\n", status);
				goto TERMINATE;
			}
	
			reduced_subprob_x_lb = lbs[0];
			reduced_subprob_y_ub = ubs[1];
			reduced_subprob_x_ub = ubs[0];
			reduced_subprob_y_lb = lbs[1];
	
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_lb,reduced_subprob_y_lb);					
			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
								reduced_subprob_y_ub,reduced_subprob_y_ub);
	
			if(! mock_insert(1,reduced_subprob_x_ub,reduced_subprob_y_ub,0,0,0,&tree) )
			{
				if(printing_in_setbranch) printf("the entire search region is dominated. Fathom.\n");
				fathoming = 1;
		  		goto TERMINATE;
			}
	
			slope = (ubs[1]-lbs[1])/(lbs[0]-ubs[0]);
			if(slope != slope) //slope < -1000. || slope > -.001) 
			{
				if(printing_in_setbranch) printf("setting slope to -1\n");
				slope = -1.;
			}
			if(printing_in_setbranch) printf("slope: %lf\n",slope);
		   	
/*		   	if(objective_space_fathoming)*/
/*		   	{*/
/*			   	if (branch_iterations > -1) pareto_branching = 1;*/
/*				if (branch_iterations > 1000)*/
/*				{*/
/*					int nodesleft = 0;*/
/*					if ((status = CPXgetcallbackinfo (env, cbdata, wherefrom, */
/*					CPX_CALLBACK_INFO_NODES_LEFT, &nodesleft))) goto TERMINATE;*/
/*		*/
/*					if(nodesleft > 1000) */
/*					{*/
/*						pareto_branching = 0;*/
/*						remove_dominated_middle = 0;*/
/*					}*/
/*					else pareto_branching = 1;*/
/*				}*/
/*			}*/
	
		/*	printf("cplex index: %d, value: %lf\n",indices[0],bd[0]);*/
	
		/*	if (branch_iterations > 100000000)*/
		/*	{*/
		/*		printf("**************reached iteration limit*****************\n");*/
		/*		exit(0);*/
		/*	}*/
		   	
		   	/* Get the objective value associated with this node. This is used for providing CPLEX with an estimated value of the solution an the next node. */
		
			status = CPXgetcallbacknodeobjval (env, cbdata, wherefrom, &objval);
		  	if ( status ) {
		      		fprintf (stdout, "Can't get node objective value.");
				goto TERMINATE;
			}
	
			status = CPXgetcallbacknodex (env, cbdata, wherefrom, x_ws, 0, cur_numcols-1);
			if ( status ) {
				printf ("CPXgetcallbacknodelp, Failed to get node x, error code %d\n", status);
				goto TERMINATE;
			}
			ws_lp_objvals[0] = x_ws[obj1_index];
			ws_lp_objvals[1] = x_ws[obj2_index];
	
			double proj1 = (slope)*(ubs[0]-x_ws[obj1_index])+x_ws[obj2_index];
			double proj2 = (1./slope)*(ubs[1]-x_ws[obj2_index])+x_ws[obj1_index];
			endpoint1_x = ubs[0];
			endpoint2_x = proj2;
			endpoint1_y = proj1;
			endpoint2_y = ubs[1];
			if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",x_ws[obj1_index],x_ws[obj2_index]);
			
			if(seqnum == 0 && generate_disjunctive_cuts_from_dominated_columns && num_dom_cols)
			{
				int j = 0;
/*				if(printing_in_setbranch) */
				printf("Attempting to generate cuts based on dominated columns\n");
				for(j=0;j<num_dom_cols;j++)
				{
/*					printf("j: %d, dominator: %d, ub: %lf, x_ws: %lf\t dominated: %d, lb: %lf, x_ws: %lf\n",j,dominating_indices[j],var_ubs[dominating_indices[j]],x_ws[dominating_indices[j]],dominated_indices[j],var_lbs[dominated_indices[j]],x_ws[dominated_indices[j]]);*/
					
					if(var_ubs[dominating_indices[j]] > x_ws[dominating_indices[j]] && var_lbs[dominated_indices[j]] < x_ws[dominated_indices[j]])
					{
						int z = 0;
						int numrows = cur_numrows;
						if(!cut_problem_created)
						{
							clock_t s = clock();
							cut_prob = create_cut_prob(env, nodelp);
							cut_problem_created = 1;
							clock_t f = clock();
							double ti = (double)(f - s) / CLOCKS_PER_SEC;
							if(printing_in_setbranch) printf("Time to create cut problem: %lf\n",ti);
						}
						clock_t s = clock();

						CPXLPptr cut_prob_copy = CPXcloneprob (env, cut_prob, &status);
	
						int rmatbeg[2] = {0,2};
						char sense[2] = {'G','G'};
						int rmatind[4] = {dominating_indices[j],cur_numcols,dominated_indices[j]+cur_numcols+1,2*cur_numcols+1};
						double rmatval[4] = {-1.,-var_ubs[dominating_indices[j]],1.,var_lbs[dominated_indices[j]]};
						double *pi = (double *) malloc ( (cur_numcols+2)*sizeof(double) );
						int *chg_ind = (int *) malloc ( cur_numcols*sizeof(int) );
						double *new_obj = (double *) malloc ( (cur_numcols+2)*sizeof(double) );
						int *obj_ind = (int *) malloc ( (cur_numcols+2)*sizeof(int) );
			/*							double *ones = (double *) malloc ( (2*numrows+4+5*cur_numcols)*sizeof(double));*/
						double *ones = (double *) malloc ( (2*numrows+4+4*cur_numcols)*sizeof(double));
			/*							int *rowlist = (int *) malloc ( (2*numrows+4+5*cur_numcols)*sizeof(int));*/
						int *rowlist = (int *) malloc ( (2*numrows+4+4*cur_numcols)*sizeof(int));
			/*							int *collist = (int *) malloc ( (2*numrows+4+5*cur_numcols)*sizeof(int));*/

						status = CPXaddrows (env, cut_prob_copy, 0, 2, 4, NULL, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
				
						for(i=0;i<cur_numcols;i++) chg_ind[i] = 2*numrows + 4*cur_numcols + i;
				
						status = CPXchgrhs (env, cut_prob_copy, cur_numcols, chg_ind, x_ws);
	
						for(i=0;i<cur_numcols;i++)
						{
							new_obj[i] = 0.;
							obj_ind[i] = i;
						}
			
						status = CPXchgobj (env, cut_prob_copy, cur_numcols, obj_ind, new_obj);
	
						double o = -1.;
						int ind[2] = {cur_numcols, 2*cur_numcols+1};
						double val[2] = {-1.,-1.};
						for(i=0;i<2*numrows+4+4*cur_numcols;i++) 
						{
							ones[i] = 1.;
							rowlist[i] = i;
						}
						status = CPXaddrows (env, cut_prob_copy, 0, 1, 2, &o, NULL, &z, ind, val, NULL, NULL);
			
						o = 1.;
						status = CPXaddcols (env, cut_prob_copy, 1, 2*numrows+4*cur_numcols, &o, &z, rowlist,ones,NULL,NULL,NULL);
	
						status = CPXlpopt (env, cut_prob_copy);
					 	if ( status ) {
					   		printf ("%s(%d): CPXlpopt, Failed to solve prob, error code %d\n", __FILE__, __LINE__, status);
							goto TERMINATE;
						}
						int lpstat_cut_prob = CPXgetstat (env, cut_prob_copy);
						if(printing_in_setbranch) printf("lpstat for cut prob: %d\n",lpstat_cut_prob);
	
						if(lpstat_cut_prob == 1 || lpstat_cut_prob == 5)
						{
							number_disj_cuts++;
							int temp_numrows = CPXgetnumrows (env, cut_prob_copy);
							double last_pi[1] = {0.};
				
							status = CPXgetpi(env, cut_prob_copy, pi, 2*numrows + 4*cur_numcols, 2*numrows + 5*cur_numcols+1);
							status = CPXgetpi(env, cut_prob_copy, last_pi, temp_numrows-1, temp_numrows-1);
				
/*							if(printing_in_setbranch) printf("last pi: %le\n",last_pi[0]);*/
				
/*							if(printing_in_setbranch) for(i=0;i<cur_numcols+1;i++) printf("pi_%d: %le\n",i,pi[i]);*/
				
							status = CPXcutcallbackaddlocal   (env,
											   cbdata,
								       			   wherefrom,
											   cur_numcols,
											   last_pi[0],
											   'L',
											   indices,
											   pi);
/*							char sen[1] = {'L'};			   */
/*							status = CPXaddusercuts (env, lp1, 1, cur_numcols, last_pi, sen, &z, indices, */
/*				    							pi, NULL);*/
										
						
/*							printf("*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*\n");*/
/*							print_inorder(tree,1);*/
/*							printf("*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*\n");*/
/*						*/
/*						*/
/*							printf("_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-\n");*/
/*							PSA_all(env,nodelp);*/
/*							printf("_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-\n");*/
/*						*/
/*							status = CPXwriteprob (env, nodelp, "myprob1.lp", "LP");*/
											   
/*							status = CPXaddrows (env, nodelp, 0, 1, cur_numcols, last_pi, sen, rmatbeg, indices,*/
/*											pi, NULL, NULL);*/
										
/*							status = CPXwriteprob (env, nodelp, "myprob2.lp", "LP");*/
/*										*/
/*							printf("_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-\n");*/
/*							PSA_all(env,nodelp);*/
/*							printf("_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-\n");*/
						}
						clock_t f = clock();
						double ti = (double)(f - s) / CLOCKS_PER_SEC;
						if(printing_in_setbranch) printf("Time to generate cut: %lf\n",ti);
						time_generating_disjunction_cuts += ti;
			
						free(pi);
						free(chg_ind);
						free(new_obj);
						free(obj_ind);
						free(ones);
						free(rowlist);
/*						status = CPXwriteprob (env, cut_prob_copy, "myprob.lp", "LP");*/
						CPXfreeprob(env, &cut_prob_copy);
/*						exit(0);*/
					}
				}
			}
/*			exit(0);*/
	
/*			status = CPXgetcallbacknodeinfo(	env,*/
/*				                 		cbdata,*/
/*				                 		wherefrom,*/
/*				                 		0,*/
/*				               			CPX_CALLBACK_INFO_NODE_SEQNUM_LONG,*/
/*				               			&seqnum);*/
/*			if ( status ) {*/
/*				printf ("CPXgetcallbacknodeinfo, Failed to get seqnum, error code %d\n", status);*/
/*				goto TERMINATE;*/
/*			}*/
/*	*/
/*			if(printing_in_setbranch) printf("at branching seqnum is: %d\n",seqnum);*/
	
		/*	if(seqnum == 0) printf("iterations: %d\n",branch_iterations);*/

			if(seqnum == 0 && objective_space_fathoming == 0) 
			{
				if(printing_in_setbranch) printf("(%d) changing dom mid to 0\n",__LINE__);
				remove_dominated_middle = 0;
			}
	
			if(printing_in_setbranch){
			printf("_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-\n");
			PSA_all(env,nodelp);
			printf("_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-_*-\n");}
	
			status = CPXgetcallbacknodeinfo(	env,
				                 		cbdata,
				                 		wherefrom,
				                 		0,
				                 		CPX_CALLBACK_INFO_NODE_DEPTH,
				                 		&depth);
			if ( status ) {
				printf ("CPXgetcallbacknodeinfo, Failed to get node depth, error code %d\n", status);
				goto TERMINATE;
			}
	
/*			status = CPXgetcallbacknodeinfo(	env,*/
/*					         		cbdata,*/
/*					         		wherefrom,*/
/*					         		0,*/
/*					        		CPX_CALLBACK_INFO_NODE_USERHANDLE,*/
/*					         		&userhandle_current);*/
	
			if(printing_in_setbranch){
			printf("*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*\n");
			print_inorder(tree,1);
			printf("*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*\n");
			}
			
/*			if( seqnum == 3)*/
/*			{*/
/*				status = CPXwriteprob (env, nodelp, "myprob2.lp", "LP");*/
/*				exit(0);*/
/*			}*/
	

			if(there_will_only_be_points) 
			{
				points_only = 0;
			}
			if(points_only && seqnum == 0)
			{
				clock_t current_time = clock();	
				cumulative_time = (double)(current_time - start_BB) / CLOCKS_PER_SEC;
				status = CPXsetdblparam (env_just_solve_mips, CPX_PARAM_TILIM, max_time - cumulative_time);
				if ( status ) {
				    	printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);
				   	exit(0);
				}
				goto SOLVE_OB2_MIP;
			}
			else
			{
				if(time_vs_node_lim)
				{
					status = CPXsetintparam (env_just_solve_mips, CPX_PARAM_NODELIM, 0);
				    	if ( status ) {
				    		printf ("Failure to set MIP node limit to 0, error %d.\n",status);
				    		exit(0);
				    	}
				}
				else
				{
					status = CPXsetdblparam (env_just_solve_mips, CPX_PARAM_TILIM, time_limit);
					if ( status ) {
					    	printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);
					   	exit(0);
					}
				}
			}

			nodelp_copy = CPXcloneprob (env, nodelp2, &status);
			if ( status ) {
				printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
				goto TERMINATE;
			}
	
			if(userhandle_current && userhandle_current->ws_still_feas) goto SOLVE_WS_MIP;
	
			chg_coefs(env, nodelp_copy, indexes, slope);
	
		/*	if(x_ws[obj1_index] > reduced_subprob_x_lb && x_ws[obj2_index] < reduced_subprob_y_ub && */
		/*				x_ws[obj1_index] < reduced_subprob_x_ub && x_ws[obj2_index] > reduced_subprob_y_lb) */
		/*	{*/
		/*		ws_sol_interior = 1;*/
		/*		slope = initial_slope;*/
		/*		goto CHECK_FEASIBILITY;*/
		/*	}*/
	
			status = CPXlpopt (env, nodelp_copy);
		 	if ( status ) {
		   		printf ("%s(%d): CPXlpopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
				goto TERMINATE;
			}
	
			lpstat = CPXgetstat (env, nodelp_copy);
			if(lpstat == 3 || lpstat == 4)
			{
/*				printf("the ws lp is infeasible, fathom\n");*/
				fathoming = 1;
		  		goto TERMINATE;
			}
	
		  	status = CPXgetx (env, nodelp_copy, x_ws, 0, cur_numcols-1);
			if(status) 
			{
				printf ("(%d) CPXgetx, Failed to get x values, error code %d\n", __LINE__,status);
				goto TERMINATE;
			}
			ws_lp_objvals[0] = x_ws[obj1_index];
			ws_lp_objvals[1] = x_ws[obj2_index];
			endpoint1_x = x_ws[obj1_index];
			endpoint2_x = x_ws[obj1_index];
			endpoint1_y = x_ws[obj2_index];
			endpoint2_y = x_ws[obj2_index];
	
			if(x_ws[obj1_index] > reduced_subprob_x_lb && x_ws[obj2_index] < reduced_subprob_y_ub && 
						x_ws[obj1_index] < reduced_subprob_x_ub && x_ws[obj2_index] > reduced_subprob_y_lb) ws_sol_interior = 1;
			else
			{
				if(x_ws[obj2_index] >= reduced_subprob_y_ub) 
				{
		/*			ob2_lp_been_solved = 1;*/
					ws_lp_sol_on_top_bd  = 1; 
				}
				if(x_ws[obj1_index] >= reduced_subprob_x_ub) 
				{
		/*			ob1_lp_been_solved = 1;*/
					ws_lp_sol_on_right_bd = 1;
				}
			}
	
			/*************** Here we check to see if the solutions found 
					at the parent node are still feasible at this one ***************************/
	
			CHECK_FEASIBILITY:
	
/*			if(printing_in_setbranch) printf("cplex index: %d\n",indices[0]);*/

			if(control_node_selection)
			{
				userhandle_up->f1 = x_ws[obj1_index];
				userhandle_down->f1 = x_ws[obj1_index];
				userhandle_up->f2 = x_ws[obj2_index];
				userhandle_down->f2 = x_ws[obj2_index];
			}
	
			for(i=0;i<total_num_integer;i++)
			{
				k = integer_indices[i];
			  	diff = x_ws[k] - floor(x_ws[k]);
		/*		printf("diff%d: %lf\n",i,diff);*/
		  		if( diff >= .00001 && diff <= .99999)
		  		{
		  			if(printing_in_setbranch) printf("changing frac index to %d (%d)\n",integer_indices[i],__LINE__);
		/*  			printf("plot(%lf,%lf,'go');\n",x_ws[obj1_index],x_ws[obj2_index]);*/
			  		frac_index = k;
			  		frac_val = x_ws[k];
					if(frac_scores[i] > 0.0001) frac_scores[i] += 1.;
					else 
					{
						frac_scores[i] += 1. + multiplier*k;
						num_frac++;
					}
					frac_values[k] = frac_val;
				}
			}
			if(num_frac == 0)
			{
				if(printing_in_setbranch) printf("ws lp solution is integer feasible\n");
		/*		for(i=0;i<total_num_integer;i++) printf("x%d: %lf\n",integer_indices[i],x_ws[integer_indices[i]]);*/
				all_feas = 0;
		
				ws_mip_opt = 1;
				ws_mip_objvals[0] = x_ws[obj1_index];
				ws_mip_objvals[1] = x_ws[obj2_index];
				
				if(!userhandle_up->x_ws) userhandle_up->x_ws = calloc ((cur_numcols+1),sizeof(double));
				if(!userhandle_down->x_ws) userhandle_down->x_ws = calloc ((cur_numcols+1),sizeof(double));
	
				for(i=0;i<cur_numcols;i++) 
				{
					userhandle_up->x_ws[i] = x_ws[i];
					userhandle_down->x_ws[i] = x_ws[i];
				}
				userhandle_up->x_ws[cur_numcols] = slope;
				userhandle_down->x_ws[cur_numcols] = slope;
		
				ws_lp_int_feas = 1;
			}
	
			if(indices[0] > integer_indices[0] && indices[0] < cur_numcols)
			{
				k = 0;
				for(i=0;i<total_num_integer;i++)
				{
					k = integer_indices[i];
					if(k == indices[0])
					{
						k = i;
						break;
					}
				}
	
		/*		printf("cplex index: %d, cols: %d, k: %d, i: %d\n",indices[0],cur_numcols,k,i);*/
				if(frac_scores[k] > .0001) frac_scores[k] += 1.;
				else frac_scores[k] += 1. + multiplier*indices[0];
				frac_values[indices[0]] = x_ws[indices[0]];
				num_frac++;
			}
	
			/*************** We now start to process the weighted sum LP solution ****************************/

			if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",x_ws[obj1_index],x_ws[obj2_index]);
			add_check = mock_insert(1,x_ws[obj1_index],x_ws[obj2_index],0,0,0,&tree);
			if(add_check)
			{
				if(printing_in_setbranch) printf("weighted sum solution not dominated\n");
				if(ws_lp_int_feas)
				{
					if(printing_in_setbranch) printf("weighted sum solution is integer feasible\n");
			
					ws_mip_opt = 1;
					ws_mip_objvals[0] = x_ws[obj1_index];
					ws_mip_objvals[1] = x_ws[obj2_index];
			
					add_check = mock_insert(1,x_ws[obj1_index],x_ws[obj2_index],0,0,0,&tree);
					if(add_check)
					{
						if(printing_in_setbranch) printf("adding solution\n");
						for(i=0;i<cur_numcols;i++)
				      		{
				      			stored_x[x_rotation][i] = x_ws[i];
				      		}
				      		x_rotation = (x_rotation + 1) % num_x_to_store;
				      		add_check = 0;
				      		PSA_full(env,NULL,x_ws,NULL,NULL);
					}
			
					if(!ws_sol_interior) 
					{
						status = CPXgetbase (env, nodelp_copy, cstat_ws, rstat_ws);
						if(status) 
						{
							printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
							goto TERMINATE;
						}
					}
					else 
					{
						CPXgetbase (env, nodelp, cstat_ws, rstat_ws);
						if(status) 
						{
							printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
							goto TERMINATE;
						}
					}
					if(printing_in_setbranch) printf("status: %d\n",status);
			
			  		sub_pr1_x_ub = x_ws[obj1_index];
					sub_pr1_y_lb = x_ws[obj2_index];
					sub_pr2_x_lb = x_ws[obj1_index];
					sub_pr2_y_ub = x_ws[obj2_index];
				  	PSA_right_check = PSA(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, seqnum, 0, nodelp_copy);
				  	PSA_left_check = PSA_left(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, seqnum, 0, nodelp_copy);
				  	if(within_PSA_score > 1.) within_PSA_score -= 1.;
				  	
					if(PSA_right_check == 2 && PSA_left_check == 2)
			  		{
				  		if(printing_in_setbranch) printf("fathoming node %d for PSA completion (%d)\n",seqnum,__LINE__);
					  	fathomed_by_PSA_completion++;
						fathoming = 1;
				  		goto TERMINATE;
			  		}
				  	else if(PSA_right_check == 2)
				  	{
						if(printing_in_setbranch) printf("the left subproblem is now empty by PSA completion\n");
				  		reduced_subprob_x_lb = sub_pr2_x_lb;
				  		reduced_subprob_y_ub = sub_pr2_y_ub;
				  		reduced_subprob_x_ub = ubs[0];
						reduced_subprob_y_lb = lbs[1];
						if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
						bound_reduction = 1;
						left_side_dom = 1;
				  		goto SOLVE_OB1_LP;
				  	}
				  	else if(PSA_left_check == 2)
					{
			  			if(printing_in_setbranch) printf("the right subproblem is now empty by PSA completion\n");
				  		reduced_subprob_x_lb = lbs[0];
				  		reduced_subprob_y_ub = ubs[1];
						reduced_subprob_x_ub = sub_pr1_x_ub;
			  			reduced_subprob_y_lb = sub_pr1_y_lb;
			  			if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
			  			bound_reduction = 1;
			  			right_side_dom = 1;
				  		goto SOLVE_OB2_LP;
				  	}
				  	else if(pareto_branching && sub_pr2_x_lb > sub_pr1_x_ub)
					{
						pareto = 1;
						goto TERMINATE;
/*				  		goto PARETO_BRANCH;*/
				  	}
				  	else
				  	{
				  		if(printing_in_setbranch) printf("none of the above\n");
				  	}
			
					if(!ob1_lp_been_solved) goto SOLVE_OB1_LP;
					else if(!ob1_lp_int_feas && !ob1_mip_ran) goto SOLVE_OB1_MIP;
					else goto BRANCHING;
				}
				else
				{
					if(printing_in_setbranch) printf("weighted sum solution is not integer feasible\n");
			
					/****************************************************************************
						 This is a scheme for detecting if solutions are all singletons in the 
						 objective space. We then branch in the objective space 
					****************************************************************************/
			
		/*			DO_IT_AGAIN:*/
		/*			;*/
		/*			closest_nodes *two_nodes = find_two_nodes_right_of_val(reduced_subprob_x_lb, reduced_subprob_y_ub, tree);*/
		/*			if(!two_nodes) */
		/*			{*/
		/*				if(points_only || its_been_only_points) goto SOLVE_OB2_MIP;*/
		/*				else goto AFTER_THIS1;*/
		/*			}*/
		/*			if(printing_in_setbranch) printf("the two nodes:\n");*/
		/*			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->closest->nw_x,*/
		/*				x_ideal-two_nodes->closest->se_x,y_ideal-two_nodes->closest->nw_y,y_ideal-two_nodes->closest->se_y);*/
		/*			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->next->nw_x,x_ideal-two_nodes->next->se_x,*/
		/*				y_ideal-two_nodes->next->nw_y,y_ideal-two_nodes->next->se_y);*/
		/*			if(two_nodes->closest->type == 2 && fabs(two_nodes->closest->nw_x - two_nodes->closest->se_x) < .00000001 && */
		/*				fabs(two_nodes->closest->nw_y - two_nodes->closest->se_y) < .00000001) two_nodes->closest->type = 1;*/
		/*			if(two_nodes->next->type == 2 && fabs(two_nodes->next->nw_x - two_nodes->next->se_x) < .00000001 && */
		/*				fabs(two_nodes->next->nw_y - two_nodes->next->se_y) < .00000001) two_nodes->next->type = 1;*/
		/*			*/
		/*			if(two_nodes->closest->type == 1 && ((fabs(two_nodes->closest->nw_x - two_nodes->next->nw_x) < .0000001 && */
		/*				two_nodes->closest->nw_y - two_nodes->next->nw_y >= -.0000001) || */
		/*				(fabs(two_nodes->closest->nw_x - two_nodes->next->se_x) < .0000001 && */
		/*				two_nodes->closest->nw_y - two_nodes->next->se_y >= -.0000001)))*/
		/*			{*/
		/*				if(printing_in_setbranch) */
		/*					printf("the closest node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
		/*				delete_node(two_nodes->closest);*/
		/*				free(two_nodes);*/
		/*				goto DO_IT_AGAIN;*/
		/*			}*/
		/*			*/
		/*			if(x_ideal - two_nodes->closest->nw_x - reduced_subprob_x_ub > .0001 || */
		/*				y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_lb < -.0001 || */
		/*				x_ideal - two_nodes->closest->nw_x - reduced_subprob_x_lb < -.0001 || */
		/*				y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_ub > .0001)*/
		/*			{*/
		/*				if(printing_in_setbranch) printf("closest node is outside region\n");*/
		/*				free(two_nodes);*/
		/*				goto SOLVE_WS_MIP;*/
		/*			}*/
		/*			*/
		/*			else if(two_nodes->next->type == 1 && ((fabs(two_nodes->next->nw_x - two_nodes->closest->nw_x) < .0000001 && */
		/*				two_nodes->next->nw_y - two_nodes->closest->nw_y >= -.0000001) || */
		/*				(fabs(two_nodes->next->nw_x - two_nodes->closest->se_x) < .0000001 && */
		/*				two_nodes->next->nw_y - two_nodes->closest->se_y >= -.0000001)))*/
		/*			{*/
		/*				if(printing_in_setbranch) */
		/*					printf("the next node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
		/*				delete_node(two_nodes->next);*/
		/*				free(two_nodes);*/
		/*				goto DO_IT_AGAIN;*/
		/*			}*/
		/*			*/
		/*			if(!there_will_only_be_points && two_nodes->closest->type == 1 && two_nodes->next->type == 1)*/
		/*			{*/
		/*				points_only = 1;*/
		/*				its_been_only_points = 1;*/
		/*				status = CPXsetdblparam (env_just_solve_mips, CPX_PARAM_TILIM, pow(10.,75.));*/
		/*				if ( status ) {*/
		/*				    	printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);*/
		/*				   	exit(0);*/
		/*				}*/
		/*				if(x_ideal - two_nodes->next->se_x < x_ideal - two_nodes->closest->nw_x )*/
		/*				{*/
		/*					if(y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_lb > .0001)*/
		/*					{*/
		/*						printf("there was only one node right of the val. Pareto branch based on its location\n");*/
		/*						sub_pr1_x_ub = reduced_subprob_x_ub;*/
		/*						sub_pr1_y_lb = y_ideal - two_nodes->closest->nw_y;*/
		/*					 	sub_pr2_x_lb = x_ideal - two_nodes->closest->nw_x + .001;*/
		/*						sub_pr2_y_ub = y_ideal - two_nodes->closest->nw_y;*/
		/*						free(two_nodes);*/
		/*						goto PARETO_BRANCH;*/
		/*					}*/
		/*				}*/
		/*				else if(x_ideal - two_nodes->next->se_x - reduced_subprob_x_ub <= .0000001 && */
		/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001)*/
		/*				{*/
		/*					double ran = (double) rand() / ( (double) RAND_MAX);*/
		/*					if(printing_in_setbranch) printf("left pt of next node is also inside the search region.\n");*/
		/*					if(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y) < -.0001 && */
		/*						y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001 && */
		/*						y_ideal - two_nodes->next->se_y - reduced_subprob_y_ub <= .0000001)*/
		/*					{*/
		/*						if(printing_in_setbranch) printf("separation between y-values, split\n");*/
		/*						printf("y separation is %lf percent of y_range\n", */
		/*							100.*(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y))/y_range);*/
		/*						sub_pr1_x_ub = reduced_subprob_x_ub;*/
		/*						sub_pr1_y_lb = y_ideal - two_nodes->closest->nw_y;*/
		/*					 	sub_pr2_x_lb = x_ideal - two_nodes->closest->nw_x + .001;*/
		/*						sub_pr2_y_ub = y_ideal - two_nodes->closest->nw_y;*/
		/*						free(two_nodes);*/
		/*						goto PARETO_BRANCH;*/
		/*					}*/
		/*					else if(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x) > .0001)*/
		/*					{*/
		/*						if(printing_in_setbranch) printf("separation between x-values, split\n");*/
		/*						printf("x separation is %lf percent of x_range\n", */
		/*							100.*(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x))/x_range);*/
		/*						exit(0);*/
		/*						sub_pr1_x_ub = ((x_ideal - two_nodes->next->se_x) + (x_ideal - two_nodes->closest->nw_x))/2.;*/
		/*						sub_pr1_y_lb = y_ideal - two_nodes->next->se_y + .001;//reduced_subprob_y_lb;*/
		/*					 	sub_pr2_x_lb = ((x_ideal - two_nodes->next->se_x) + (x_ideal - two_nodes->closest->nw_x))/2.;*/
		/*						sub_pr2_y_ub = reduced_subprob_y_ub;*/
		/*						free(two_nodes);*/
		/*						goto PARETO_BRANCH;*/
		/*					}*/
		/*				}*/
		/*				else*/
		/*				{*/
		/*					printf("the 2nd node was out of the region. Go back to beginning\n");*/
		/*					ob2_mip_ran = 0;*/
		/*					ob2_mip_opt = 0;*/
		/*					goto BEGINNING;*/
		/*				}*/
		/*			}*/
		/*			else */
		/*			{*/
		/*				points_only = 0;*/
		/*				its_been_only_points = 0;*/
		/*				status = CPXsetdblparam (env_just_solve_mips, CPX_PARAM_TILIM, time_limit);*/
		/*				if ( status ) {*/
		/*				    	printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);*/
		/*				   	exit(0);*/
		/*				}*/
		/*				*/
		/*				goto AFTER_THIS1;*/
		/*				*/
		/*				if(x_ideal - two_nodes->next->se_x < x_ideal - two_nodes->closest->nw_x )*/
		/*				{*/
		/*					if(y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_lb > .0001)*/
		/*					{*/
		/*						goto AFTER_THIS1;*/
		/*					}*/
		/*				}*/
		/*				else if(x_ideal - two_nodes->next->se_x - reduced_subprob_x_ub <= .0000001 && */
		/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001)*/
		/*				{*/
		/*					double ran = (double) rand() / ( (double) RAND_MAX);*/
		/*					if(printing_in_setbranch) printf("left pt of next node is also inside the search region.\n");*/
		/*					if( (y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y) < -.0001 && */
		/*						y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001 && */
		/*						y_ideal - two_nodes->next->se_y - reduced_subprob_y_ub <= .0000001) && */
		/*						(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x) > .0001) )*/
		/*					{*/
		/*						if(printing_in_setbranch) printf("separation between x and y-values, split\n");*/
		/*						double y_sep = fabs(100.*(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y))/y_range);*/
		/*						double x_sep = fabs(100.*(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x))/x_range);*/
		/*						printf("y separation is %lf percent of y_range\n", */
		/*							100.*(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y))/y_range);*/
		/*						printf("x separation is %lf percent of x_range\n", */
		/*							100.*(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x))/x_range);*/
		/*						if(x_sep > 5. || y_sep > 5.)*/
		/*						{*/
		/*							printf("(%d) exploiting separation\n",__LINE__);*/
		/*							sub_pr1_x_ub = reduced_subprob_x_ub;*/
		/*							sub_pr1_y_lb = y_ideal - two_nodes->closest->nw_y;*/
		/*						 	sub_pr2_x_lb = x_ideal - two_nodes->closest->nw_x + .001;*/
		/*							sub_pr2_y_ub = y_ideal - two_nodes->closest->nw_y;*/
		/*							free(two_nodes);*/
		/*							goto PARETO_BRANCH;*/
		/*						}*/
		/*						else goto AFTER_THIS1;*/
		/*					}*/
		/*				}*/
		/*			}*/
		/*			AFTER_THIS1:*/
		/*			*/
		/*			free(two_nodes);*/
					goto SOLVE_WS_MIP;
/*					goto BRANCHING;*/
				}
			}
			else
			{
				if(printing_in_setbranch) printf("weighted sum solution dominated\n");

				ws_lp_dom = 1;
		
				if(ws_lp_int_feas)
				{
					if(printing_in_setbranch) printf("ws soln was int feas\n");
					if(!ws_sol_interior) 
					{
						status = CPXgetbase (env, nodelp_copy, cstat_ws, rstat_ws);
						if(status) 
						{
							printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
							goto TERMINATE;
						}
					}
					else 
					{
						CPXgetbase (env, nodelp, cstat_ws, rstat_ws);
						if(status) 
						{
							printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
							goto TERMINATE;
						}
					}
					if(printing_in_setbranch) printf("status: %d\n",status);
			
					sub_pr1_x_ub = x_ws[obj1_index];
					sub_pr1_y_lb = x_ws[obj2_index];
					sub_pr2_x_lb = x_ws[obj1_index];
					sub_pr2_y_ub = x_ws[obj2_index];
				  	PSA_right_check = PSA(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, seqnum, 0, nodelp_copy);
				  	PSA_left_check = PSA_left(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, seqnum, 0, nodelp_copy);
				  	
					if(PSA_right_check == 2 && PSA_left_check == 2)
			  		{
				  		if(printing_in_setbranch) printf("fathoming node %d for PSA completion (%d)\n",seqnum,__LINE__);
					  	fathomed_by_PSA_completion++;
						fathoming = 1;
/*						*useraction_p = CPX_CALLBACK_SET;*/
/*						nodecnt = 0;*/
				  		goto TERMINATE;
			  		}
				  	else if(PSA_right_check == 2)
				  	{
						if(printing_in_setbranch) printf("the left subproblem is now empty by PSA completion\n");
				  		reduced_subprob_x_lb = x_ws[obj1_index];
				  		if(there_will_only_be_points && integer_bb && integer_objective == 1) reduced_subprob_x_lb += smallest_coef;
				  		reduced_subprob_y_ub = x_ws[obj2_index];
				  		if(there_will_only_be_points && integer_bb && integer_objective == 2) reduced_subprob_y_ub -= smallest_coef;
				  		reduced_subprob_x_ub = ubs[0];
						reduced_subprob_y_lb = lbs[1];
						if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
						bound_reduction = 1;
						left_side_dom = 1;
				  		goto SOLVE_OB1_LP;
				  	}
				  	else if(PSA_left_check == 2)
					{
			  			if(printing_in_setbranch) printf("the right subproblem is now empty by PSA completion\n");
				  		reduced_subprob_x_lb = lbs[0];
				  		reduced_subprob_y_ub = ubs[1];
						reduced_subprob_x_ub = x_ws[obj1_index];
						if(there_will_only_be_points && integer_bb && integer_objective == 1) reduced_subprob_x_ub -= smallest_coef;
			  			reduced_subprob_y_lb = x_ws[obj2_index];
			  			if(there_will_only_be_points && integer_bb && integer_objective == 2) reduced_subprob_y_lb += smallest_coef;
			  			if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
			  			bound_reduction = 1;
			  			right_side_dom = 1;
				  		goto SOLVE_OB2_LP;
				  	}
				  	else if(sub_pr2_x_lb > sub_pr1_x_ub)
					{
						if(generate_disjunctive_cuts_from_obj_space_disjunctions && (PSA_iter_cnt > 2 || PSA_left_iter_cnt > 2))
						{
							if(printing_in_setbranch) printf("Attempting to generate cuts based on the disjunction created by removing dominated portion of objective space\n");
					
							int z = 0;
							int numrows = cur_numrows;
							if(!cut_problem_created)
							{
								clock_t s = clock();
								cut_prob = create_cut_prob(env, nodelp2);
								cut_problem_created = 1;
								clock_t f = clock();
								double ti = (double)(f - s) / CLOCKS_PER_SEC;
								if(printing_in_setbranch) printf("Time to create cut problem: %lf\n",ti);
							}
							clock_t s = clock();
			
							CPXLPptr cut_prob_copy = CPXcloneprob (env, cut_prob, &status);
					
							int rmatbeg[4] = {0,2,4,6};
							char sense[4] = {'G','G','G','G'};
							int rmatind[8] = {obj1_index,cur_numcols,obj2_index,cur_numcols,obj1_index+cur_numcols+1,2*cur_numcols+1,
										obj2_index+cur_numcols+1,2*cur_numcols+1};
							double rmatval[8] = {1.,sub_pr1_x_ub,-1.,-sub_pr1_y_lb,-1.,-sub_pr2_x_lb,1.,sub_pr2_y_ub};
							double *pi = (double *) malloc ( (cur_numcols+2)*sizeof(double) );
							int *chg_ind = (int *) malloc ( cur_numcols*sizeof(int) );
							double *new_obj = (double *) malloc ( (cur_numcols+2)*sizeof(double) );
							int *obj_ind = (int *) malloc ( (cur_numcols+2)*sizeof(int) );
/*							double *ones = (double *) malloc ( (2*numrows+4+5*cur_numcols)*sizeof(double));*/
							double *ones = (double *) malloc ( (2*numrows+4+4*cur_numcols)*sizeof(double));
/*							int *rowlist = (int *) malloc ( (2*numrows+4+5*cur_numcols)*sizeof(int));*/
							int *rowlist = (int *) malloc ( (2*numrows+4+4*cur_numcols)*sizeof(int));
/*							int *collist = (int *) malloc ( (2*numrows+4+5*cur_numcols)*sizeof(int));*/
			
							status = CPXaddrows (env, cut_prob_copy, 0, 4, 8, NULL, sense, rmatbeg, rmatind,
										rmatval, NULL, NULL);
								
							for(i=0;i<cur_numcols;i++) chg_ind[i] = 2*numrows + 4*cur_numcols + i;
								
							if(printing_in_setbranch) printf("PSA iter cnt: %d\t PSA left iter cnt: %d\n",PSA_iter_cnt,PSA_left_iter_cnt);
							if(PSA_iter_cnt <= 2 || PSA_left_iter_cnt <= 2)
							{
								if(PSA_iter_cnt > 2) 
								{
									if(printing_in_setbranch) printf("Using an intermediate solution found during rr in place of x_ws\n");
									if(printing_in_setbranch) printf("plot(%lf,%lf,'ko');\n",temp_x_r[obj1_index],temp_x_r[obj2_index]);
									status = CPXchgrhs (env, cut_prob_copy, cur_numcols, chg_ind, temp_x_r);
								}
								else if(PSA_left_iter_cnt > 2) 
								{
									if(printing_in_setbranch) printf("Using an intermediate solution found during rl in place of x_ws\n");
									if(printing_in_setbranch) printf("plot(%lf,%lf,'ko');\n",temp_x_l[obj1_index],temp_x_l[obj2_index]);
									status = CPXchgrhs (env, cut_prob_copy, cur_numcols, chg_ind, temp_x_l);
								}
							}
							else status = CPXchgrhs (env, cut_prob_copy, cur_numcols, chg_ind, x_ws);
					
							for(i=0;i<cur_numcols;i++)
							{
								new_obj[i] = 0.;
								obj_ind[i] = i;
							}
							
							status = CPXchgobj (env, cut_prob_copy, cur_numcols, obj_ind, new_obj);
					
							double o = -1.;
							int ind[2] = {cur_numcols, 2*cur_numcols+1};
							double val[2] = {-1.,-1.};
							for(i=0;i<2*numrows+4+4*cur_numcols;i++) 
							{
								ones[i] = 1.;
								rowlist[i] = i;
							}
							status = CPXaddrows (env, cut_prob_copy, 0, 1, 2, &o, NULL, &z, ind, val, NULL, NULL);
							
							o = 1.;
							status = CPXaddcols (env, cut_prob_copy, 1, 2*numrows+4*cur_numcols, &o, &z, rowlist,ones,NULL,NULL,NULL);
					
							status = CPXlpopt (env, cut_prob_copy);
						 	if ( status ) {
						   		printf ("%s(%d): CPXlpopt, Failed to solve prob, error code %d\n", __FILE__, __LINE__, status);
								goto TERMINATE;
							}
							int lpstat_cut_prob = CPXgetstat (env, cut_prob_copy);
							if(printing_in_setbranch) printf("lpstat for cut prob: %d\n",lpstat_cut_prob);
					
							if(lpstat_cut_prob == 1 || lpstat_cut_prob == 5)
							{
								number_disj_cuts++;
								int temp_numrows = CPXgetnumrows (env, cut_prob_copy);
								double last_pi[1] = {0.};
								
								status = CPXgetpi(env, cut_prob_copy, pi, 2*numrows + 4*cur_numcols, 2*numrows + 5*cur_numcols+1);
								status = CPXgetpi(env, cut_prob_copy, last_pi, temp_numrows-1, temp_numrows-1);
								
								if(printing_in_setbranch) printf("last pi: %le\n",last_pi[0]);
								
								if(printing_in_setbranch) for(i=0;i<cur_numcols+1;i++) printf("pi_%d: %le\n",i,pi[i]);
								
								status = CPXcutcallbackaddlocal   (env,
												   cbdata,
									       			   wherefrom,
												   cur_numcols,
												   last_pi[0],
												   'L',
												   indices,
												   pi);
							}
							clock_t f = clock();
							double ti = (double)(f - s) / CLOCKS_PER_SEC;
							if(printing_in_setbranch) printf("Time to generate cut: %lf\n",ti);
							time_generating_disjunction_cuts += ti;
							
							free(pi);
							free(chg_ind);
							free(new_obj);
							free(obj_ind);
							free(ones);
							free(rowlist);
							CPXfreeprob(env, &cut_prob_copy);
						}
						if(pareto_branching)
						{
							pareto = 1;
/*							printf("going to terminate and also pareto branching\n");*/
							goto TERMINATE;
						}
				  	}
				  	else
				  	{
				  		if(printing_in_setbranch) printf("none of the above\n");
		/*		  		for(i=0;i<total_num_integer;i++)*/
		/*				{*/
		/*					k = integer_indices[i];*/
		/*					if(k == frac_index) */
		/*					{*/
		/*						k = i;*/
		/*						break;*/
		/*					}*/
		/*				}*/
		/*				frac_scores[k] += 3.;*/
		/*				frac_values[frac_index] = frac_val;*/
				  	}
				}
				else if(generate_disjunctive_cuts_from_obj_space_disjunctions || remove_dominated_middle)
				{
					if(!ws_sol_interior) 
					{
						status = CPXgetbase (env, nodelp_copy, cstat_ws, rstat_ws);
						if(status) 
						{
							printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
							goto TERMINATE;
						}
					}
					else 
					{
						CPXgetbase (env, nodelp, cstat_ws, rstat_ws);
						if(status) 
						{
							printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
							goto TERMINATE;
						}
					}
					if(printing_in_setbranch) printf("status: %d\n",status);

					sub_pr1_x_ub = x_ws[obj1_index];
					sub_pr1_y_lb = x_ws[obj2_index];
					sub_pr2_x_lb = x_ws[obj1_index];
					sub_pr2_y_ub = x_ws[obj2_index];
				  	PSA_left_check = PSA_reduce_left(env, nodelp2, x_ws, cstat_ws, rstat_ws, indexes);
				  	PSA_right_check = PSA_reduce_right(env, nodelp2, x_ws, cstat_ws, rstat_ws, indexes, seqnum);
				  	
				  	if(PSA_left_check == 2 && PSA_right_check == 2)
				  	{
				  		fathomed_by_dominated_lb++;
				  		fathoming = 1;
				  		goto TERMINATE;
				  	}
					if(PSA_left_check == 2)
					{
			  			if(printing_in_setbranch) printf("the right subproblem is now empty by PSA completion\n");
						reduced_subprob_x_ub = x_ws[obj1_index];
			  			reduced_subprob_y_lb = x_ws[obj2_index];
			  			if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
			  			bound_reduction = 1;
			  			right_side_dom = 1;
				  		goto SOLVE_OB2_LP;
				  	}
				  	if(PSA_right_check == 2)
					{
			  			if(printing_in_setbranch) printf("the left subproblem is now empty by PSA completion\n");
						reduced_subprob_x_lb = x_ws[obj1_index];
			  			reduced_subprob_y_ub = x_ws[obj2_index];
			  			if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
			  			bound_reduction = 1;
			  			left_side_dom = 1;
				  		goto SOLVE_OB1_LP;
				  	}
				  	if(sub_pr2_x_lb > sub_pr1_x_ub)
					{
						if(generate_disjunctive_cuts_from_obj_space_disjunctions && (PSA_rr_iter_cnt > 2 || PSA_rl_iter_cnt > 2))
						{
							if(printing_in_setbranch) printf("Attempting to generate cuts based on the disjunction created by removing dominated portion of objective space\n");
					
							int z = 0;
							int numrows = cur_numrows;
							if(!cut_problem_created)
							{
								clock_t s = clock();
								cut_prob = create_cut_prob(env, nodelp2);
								cut_problem_created = 1;
								clock_t f = clock();
								double ti = (double)(f - s) / CLOCKS_PER_SEC;
								if(printing_in_setbranch) printf("Time to create cut problem: %lf\n",ti);
							}
							clock_t s = clock();
			
							CPXLPptr cut_prob_copy = CPXcloneprob (env, cut_prob, &status);
					
							int rmatbeg[4] = {0,2,4,6};
							char sense[4] = {'G','G','G','G'};
							int rmatind[8] = {obj1_index,cur_numcols,obj2_index,cur_numcols,obj1_index+cur_numcols+1,2*cur_numcols+1,
										obj2_index+cur_numcols+1,2*cur_numcols+1};
							double rmatval[8] = {1.,sub_pr1_x_ub,-1.,-sub_pr1_y_lb,-1.,-sub_pr2_x_lb,1.,sub_pr2_y_ub};
							double *pi = (double *) malloc ( (cur_numcols+2)*sizeof(double) );
							int *chg_ind = (int *) malloc ( cur_numcols*sizeof(int) );
							double *new_obj = (double *) malloc ( (cur_numcols+2)*sizeof(double) );
							int *obj_ind = (int *) malloc ( (cur_numcols+2)*sizeof(int) );
/*							double *ones = (double *) malloc ( (2*numrows+4+5*cur_numcols)*sizeof(double));*/
							double *ones = (double *) malloc ( (2*numrows+4+4*cur_numcols)*sizeof(double));
/*							int *rowlist = (int *) malloc ( (2*numrows+4+5*cur_numcols)*sizeof(int));*/
							int *rowlist = (int *) malloc ( (2*numrows+4+4*cur_numcols)*sizeof(int));
/*							int *collist = (int *) malloc ( (2*numrows+4+5*cur_numcols)*sizeof(int));*/
			
							status = CPXaddrows (env, cut_prob_copy, 0, 4, 8, NULL, sense, rmatbeg, rmatind,
										rmatval, NULL, NULL);
								
							for(i=0;i<cur_numcols;i++) chg_ind[i] = 2*numrows + 4*cur_numcols + i;
								
							if(printing_in_setbranch) printf("PSA rr iter cnt: %d\t PSA rl iter cnt: %d\n",PSA_rr_iter_cnt,PSA_rl_iter_cnt);
							if(PSA_rr_iter_cnt <= 2 || PSA_rl_iter_cnt <= 2)
							{
								if(PSA_rr_iter_cnt > 2) 
								{
									if(printing_in_setbranch) printf("Using an intermediate solution found during rr in place of x_ws\n");
									if(printing_in_setbranch) printf("plot(%lf,%lf,'ko');\n",temp_x_r[obj1_index],temp_x_r[obj2_index]);
									status = CPXchgrhs (env, cut_prob_copy, cur_numcols, chg_ind, temp_x_r);
								}
								else if(PSA_rl_iter_cnt > 2) 
								{
									if(printing_in_setbranch) printf("Using an intermediate solution found during rl in place of x_ws\n");
									if(printing_in_setbranch) printf("plot(%lf,%lf,'ko');\n",temp_x_l[obj1_index],temp_x_l[obj2_index]);
									status = CPXchgrhs (env, cut_prob_copy, cur_numcols, chg_ind, temp_x_l);
								}
							}
							else status = CPXchgrhs (env, cut_prob_copy, cur_numcols, chg_ind, x_ws);
					
							for(i=0;i<cur_numcols;i++)
							{
								new_obj[i] = 0.;
								obj_ind[i] = i;
							}
							
							status = CPXchgobj (env, cut_prob_copy, cur_numcols, obj_ind, new_obj);
					
							double o = -1.;
							int ind[2] = {cur_numcols, 2*cur_numcols+1};
							double val[2] = {-1.,-1.};
							for(i=0;i<2*numrows+4+4*cur_numcols;i++) 
							{
								ones[i] = 1.;
								rowlist[i] = i;
							}
							status = CPXaddrows (env, cut_prob_copy, 0, 1, 2, &o, NULL, &z, ind, val, NULL, NULL);
							
							o = 1.;
							status = CPXaddcols (env, cut_prob_copy, 1, 2*numrows+4*cur_numcols, &o, &z, rowlist,ones,NULL,NULL,NULL);
					
							status = CPXlpopt (env, cut_prob_copy);
						 	if ( status ) {
						   		printf ("%s(%d): CPXlpopt, Failed to solve prob, error code %d\n", __FILE__, __LINE__, status);
								goto TERMINATE;
							}
							int lpstat_cut_prob = CPXgetstat (env, cut_prob_copy);
							if(printing_in_setbranch) printf("lpstat for cut prob: %d\n",lpstat_cut_prob);
					
							if(lpstat_cut_prob == 1 || lpstat_cut_prob == 5)
							{
								number_disj_cuts++;
								int temp_numrows = CPXgetnumrows (env, cut_prob_copy);
								double last_pi[1] = {0.};
								
								status = CPXgetpi(env, cut_prob_copy, pi, 2*numrows + 4*cur_numcols, 2*numrows + 5*cur_numcols+1);
								status = CPXgetpi(env, cut_prob_copy, last_pi, temp_numrows-1, temp_numrows-1);
								
/*								printf("last pi: %le\n",last_pi[0]);*/
/*								*/
/*								for(i=0;i<cur_numcols+1;i++) printf("pi_%d: %le\n",i,pi[i]);*/
								
								status = CPXcutcallbackaddlocal   (env,
												   cbdata,
									       			   wherefrom,
												   cur_numcols,
												   last_pi[0],
												   'L',
												   indices,
												   pi);
							}
							clock_t f = clock();
							double ti = (double)(f - s) / CLOCKS_PER_SEC;
							if(printing_in_setbranch) printf("Time to generate cut: %lf\n",ti);
							time_generating_disjunction_cuts += ti;
							
							free(pi);
							free(chg_ind);
							free(new_obj);
							free(obj_ind);
							free(ones);
							free(rowlist);
							CPXfreeprob(env, &cut_prob_copy);
						}
						if(pareto_branching)
						{
							pareto = 1;
							goto TERMINATE;
						}
				  	}
				}
/*				goto SOLVE_WS_MIP;*/
		
		/*		DO_IT_AGAIN3:*/
		/*		;*/
		/*		closest_nodes *two_nodes = find_two_nodes_right_of_val(reduced_subprob_x_lb, reduced_subprob_y_ub, tree);*/
		/*		if(!two_nodes) */
		/*		{*/
		/*			if(points_only || its_been_only_points) goto SOLVE_OB2_MIP;*/
		/*			else goto AFTER_THIS2;*/
		/*		}*/
		/*		if(printing_in_setbranch) printf("the two nodes:\n");*/
		/*		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->closest->nw_x,*/
		/*			x_ideal-two_nodes->closest->se_x,y_ideal-two_nodes->closest->nw_y,y_ideal-two_nodes->closest->se_y);*/
		/*		if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->next->nw_x,x_ideal-two_nodes->next->se_x,*/
		/*			y_ideal-two_nodes->next->nw_y,y_ideal-two_nodes->next->se_y);*/
		/*		if(two_nodes->closest->type == 2 && fabs(two_nodes->closest->nw_x - two_nodes->closest->se_x) < .00000001 && */
		/*			fabs(two_nodes->closest->nw_y - two_nodes->closest->se_y) < .00000001) two_nodes->closest->type = 1;*/
		/*		if(two_nodes->next->type == 2 && fabs(two_nodes->next->nw_x - two_nodes->next->se_x) < .00000001 && */
		/*			fabs(two_nodes->next->nw_y - two_nodes->next->se_y) < .00000001) two_nodes->next->type = 1;*/
		/*		if(two_nodes->closest->type == 1 && ((fabs(two_nodes->closest->nw_x - two_nodes->next->nw_x) < .0000001 && */
		/*			two_nodes->closest->nw_y - two_nodes->next->nw_y >= -.0000001) || */
		/*			(fabs(two_nodes->closest->nw_x - two_nodes->next->se_x) < .0000001 && */
		/*			two_nodes->closest->nw_y - two_nodes->next->se_y >= -.0000001)))*/
		/*		{*/
		/*			if(printing_in_setbranch) */
		/*				printf("the closest node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
		/*			delete_node(two_nodes->closest);*/
		/*			free(two_nodes);*/
		/*			goto DO_IT_AGAIN3;*/
		/*		}*/
		/*		*/
		/*		if(x_ideal - two_nodes->closest->nw_x - reduced_subprob_x_ub > .0001 || */
		/*			y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_lb < -.0001 || */
		/*			x_ideal - two_nodes->closest->nw_x - reduced_subprob_x_lb < -.0001 || */
		/*			y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_ub > .0001)*/
		/*		{*/
		/*			if(printing_in_setbranch) printf("closest node is outside region\n");*/
		/*			free(two_nodes);*/
		/*			goto SOLVE_WS_MIP;*/
		/*		}*/
		/*		*/
		/*		else if(two_nodes->next->type == 1 && ((fabs(two_nodes->next->nw_x - two_nodes->closest->nw_x) < .0000001 && */
		/*			two_nodes->next->nw_y - two_nodes->closest->nw_y >= -.0000001) || */
		/*			(fabs(two_nodes->next->nw_x - two_nodes->closest->se_x) < .0000001 && */
		/*			two_nodes->next->nw_y - two_nodes->closest->se_y >= -.0000001)))*/
		/*		{*/
		/*			if(printing_in_setbranch) */
		/*				printf("the next node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
		/*			delete_node(two_nodes->next);*/
		/*			free(two_nodes);*/
		/*			goto DO_IT_AGAIN3;*/
		/*		}*/
		/*		*/
		/*		if(!there_will_only_be_points && two_nodes->closest->type == 1 && two_nodes->next->type == 1)*/
		/*		{*/
		/*			points_only = 1;*/
		/*			its_been_only_points = 1;*/
		/*			status = CPXsetdblparam (env_just_solve_mips, CPX_PARAM_TILIM, pow(10.,75.));*/
		/*			if ( status ) {*/
		/*			    	printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);*/
		/*			   	exit(0);*/
		/*			}*/
		/*			if(x_ideal - two_nodes->next->se_x < x_ideal - two_nodes->closest->nw_x)*/
		/*			{*/
		/*				if(y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_lb > .0001)*/
		/*				{*/
		/*					printf("there was only one node right of the val. Pareto branch based on its location\n");*/
		/*					sub_pr1_x_ub = reduced_subprob_x_ub;*/
		/*					sub_pr1_y_lb = y_ideal - two_nodes->closest->nw_y;*/
		/*				 	sub_pr2_x_lb = x_ideal - two_nodes->closest->nw_x + .001;*/
		/*					sub_pr2_y_ub = y_ideal - two_nodes->closest->nw_y;*/
		/*					free(two_nodes);*/
		/*					goto PARETO_BRANCH;*/
		/*				}*/
		/*			}*/
		/*			else if(x_ideal - two_nodes->next->se_x - reduced_subprob_x_ub <= .0000001 && */
		/*				y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001)*/
		/*			{*/
		/*				if(printing_in_setbranch) printf("left pt of next node is also inside the search region.\n");*/
		/*				if(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y) < -.0001 && */
		/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001 && */
		/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_ub <= .0000001)*/
		/*				{*/
		/*					if(printing_in_setbranch) printf("separation between y-values, split\n");*/
		/*					printf("y separation is %lf percent of y_range\n", */
		/*							100.*(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y))/y_range);*/
		/*					sub_pr1_x_ub = reduced_subprob_x_ub;*/
		/*					sub_pr1_y_lb = y_ideal - two_nodes->closest->nw_y;*/
		/*				 	sub_pr2_x_lb = x_ideal - two_nodes->closest->nw_x + .001;*/
		/*					sub_pr2_y_ub = y_ideal - two_nodes->closest->nw_y;*/
		/*					free(two_nodes);*/
		/*					goto PARETO_BRANCH;*/
		/*				}*/
		/*				else if(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x) > .0001)*/
		/*				{*/
		/*					if(printing_in_setbranch) printf("separation between x-values, split\n");*/
		/*					printf("y separation is %lf percent of y_range\n", */
		/*							100.*(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x))/x_range);*/
		/*					sub_pr1_x_ub = ((x_ideal - two_nodes->next->se_x) + (x_ideal - two_nodes->closest->nw_x))/2.;*/
		/*					sub_pr1_y_lb = y_ideal - two_nodes->next->se_y + .001;//reduced_subprob_y_lb;*/
		/*				 	sub_pr2_x_lb = ((x_ideal - two_nodes->next->se_x) + (x_ideal - two_nodes->closest->nw_x))/2.;*/
		/*					sub_pr2_y_ub = reduced_subprob_y_ub;*/
		/*					free(two_nodes);*/
		/*					goto PARETO_BRANCH;*/
		/*				}*/
		/*			}*/
		/*		}*/
		/*		else */
		/*		{*/
		/*			points_only = 0;*/
		/*			its_been_only_points = 0;*/
		/*			status = CPXsetdblparam (env_just_solve_mips, CPX_PARAM_TILIM, time_limit);*/
		/*			if ( status ) {*/
		/*			    	printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);*/
		/*			   	exit(0);*/
		/*			}*/
		/*			*/
		/*			goto AFTER_THIS2;*/
		/*			*/
		/*			if(x_ideal - two_nodes->next->se_x < x_ideal - two_nodes->closest->nw_x )*/
		/*			{*/
		/*				if(y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_lb > .0001)*/
		/*				{*/
		/*					goto AFTER_THIS2;*/
		/*				}*/
		/*			}*/
		/*			else if(x_ideal - two_nodes->next->se_x - reduced_subprob_x_ub <= .0000001 && */
		/*				y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001)*/
		/*			{*/
		/*					double ran = (double) rand() / ( (double) RAND_MAX);*/
		/*				if(printing_in_setbranch) printf("left pt of next node is also inside the search region.\n");*/
		/*				if( (y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y) < -.0001 && */
		/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001 && */
		/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_ub <= .0000001) && */
		/*					(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x) > .0001) )*/
		/*				{*/
		/*					if(printing_in_setbranch) printf("separation between x and y-values, split\n");*/
		/*					double y_sep = fabs(100.*(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y))/y_range);*/
		/*					double x_sep = fabs(100.*(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x))/x_range);*/
		/*					printf("y separation is %lf percent of y_range\n", */
		/*						100.*(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y))/y_range);*/
		/*					printf("x separation is %lf percent of x_range\n", */
		/*						100.*(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x))/x_range);*/
		/*					if(x_sep > 5. || y_sep > 5.)*/
		/*					{*/
		/*						printf("(%d) exploiting separation\n",__LINE__);*/
		/*						sub_pr1_x_ub = reduced_subprob_x_ub;*/
		/*						sub_pr1_y_lb = y_ideal - two_nodes->closest->nw_y;*/
		/*					 	sub_pr2_x_lb = x_ideal - two_nodes->closest->nw_x + .001;*/
		/*						sub_pr2_y_ub = y_ideal - two_nodes->closest->nw_y;*/
		/*						free(two_nodes);*/
		/*						goto PARETO_BRANCH;*/
		/*					}*/
		/*					else goto AFTER_THIS2;*/
		/*				}*/
		/*			}*/
		/*		}*/
		/*		AFTER_THIS2:*/
		/*		*/
		/*		free(two_nodes);*/
				goto SOLVE_OB1_LP;
		
			}
	
			goto SOLVE_OB1_LP;
	
			/*************** Here we solve the weighted sum single objective MIP if its determined that we should ******/

			SOLVE_WS_MIP:

		/*	right_side_dom = 0;*/
		/*	left_side_dom = 0;*/
			ws_mip_ran = 1;
	
			if(bound_reduction)
			{
				int ind4[4] = {obj1_index,obj1_index,obj2_index,obj2_index};
				char lu4[4] = {'L','U','L','U'};
				double bds4[4] = {reduced_subprob_x_lb,reduced_subprob_x_ub,reduced_subprob_y_lb,reduced_subprob_y_ub};
				status = CPXchgbds (env, nodelp2, 4, ind4, lu4, bds4);
				if ( status ) {
			   		printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
			  	}
			}
	
			nodelp_copy2 = CPXcloneprob (env_just_solve_mips, nodelp2, &status);
			if ( status ) {
				printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
				goto TERMINATE;
			}
			CPXchgprobtype(env_just_solve_mips, nodelp_copy2, CPXPROB_MILP);
			status = CPXchgctype(env_just_solve_mips, nodelp_copy2, cur_numcols, indexes, xctype);
	
		/*	if(exact_mips) goto SKIP_THIS4;*/
			if(userhandle_current && userhandle_current->ws_still_feas)
			{
				if(printing_in_setbranch) printf("ws soln from parent is still feasible here\n");
				for(j=0;j<cur_numcols;j++) x_ws[j] = userhandle_current->x_ws[j];
				ws_mip_objvals[0] = x_ws[obj1_index];
				ws_mip_objvals[1] = x_ws[obj2_index];
				slope = userhandle_current->x_ws[cur_numcols];
				
				if(!userhandle_up->x_ws) userhandle_up->x_ws = calloc ((cur_numcols+1),sizeof(double));
				if(!userhandle_down->x_ws) userhandle_down->x_ws = calloc ((cur_numcols+1),sizeof(double));
	
				for(i=0;i<cur_numcols;i++) 
				{
					userhandle_up->x_ws[i] = x_ws[i];
					userhandle_down->x_ws[i] = x_ws[i];
				}
				userhandle_up->x_ws[cur_numcols] = slope;
				userhandle_down->x_ws[cur_numcols] = slope;
				
				if(control_node_selection)
				{
					userhandle_up->f1 = x_ws[obj1_index];
					userhandle_down->f1 = x_ws[obj1_index];
					userhandle_up->f2 = x_ws[obj2_index];
					userhandle_down->f2 = x_ws[obj2_index];
				}
	
				if(printing_in_setbranch) printf("plotting\n");
/*				if(printing_in_setbranch) */
/*				printf("plot(%lf,%lf,'mo');\n",x_ws[obj1_index],x_ws[obj2_index]);*/
	
				if(printing_in_setbranch) printf("the mip solution therefore stays optimal\n");
				ws_mip_opt = 1;
				mip_solved = 1;
		
				if(going_back) goto BRANCHING;
		
				if(!ob1_lp_been_solved) goto SOLVE_OB1_LP;
				else if(!ob1_lp_int_feas) goto SOLVE_OB1_MIP;
				else goto BRANCHING;
			}
	
		/*	printf("slope before any changes: %lf\n",slope);*/
		
			if(bound_reduction) slope = (reduced_subprob_y_ub - reduced_subprob_y_lb)/(reduced_subprob_x_lb - reduced_subprob_x_ub);
		
			if(slope != prev_slope || same_seq_flag != 1)
			{
				if(printing_in_setbranch) printf("slope: %lf\n",slope);
				prev_slope = slope;
				chg_coefs(env_just_solve_mips, nodelp_copy2, indexes, slope);
				if(ubs[0]-lbs[0] < prob_width/width_divisor) pareto_branching = 0;
			}
	
			if(!mip_solved) 
			{
				if(printing_in_setbranch) printf("here: %d\n",__LINE__);
				mip_solved = 1;
				num_nodes_with_mips_solved++;
			}
	
			if(printing_in_setbranch) printf("optimizing mip\n");
			start_mipsolve = clock();
	
			int num_starts = CPXgetnummipstarts (env_just_solve_mips, global_mip);
			nzcnt = 0;
			prev_numsols = 0;
		
		    	if(num_starts > global_num_starts)
		    	{
		    		global_num_starts = num_starts;
		    		global_startspace = cur_numcols*global_num_starts;
		    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
				global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
				global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
				global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
		    	}
	
			status = CPXgetmipstarts (env_just_solve_mips, global_mip, &nzcnt, global_beg, global_varindices, 
					   global_values, global_effortlevel, global_startspace,
					   &surplus, 0, num_starts-1);
					   
			status = CPXaddmipstarts (env_just_solve_mips, nodelp_copy2, num_starts, nzcnt, global_beg, global_varindices,
					   global_values, global_effortlevel, NULL);

			CPXmipopt (env_just_solve_mips, nodelp_copy2);
			
			lpstat = CPXgetstat (env_just_solve_mips, nodelp_copy2);
/*			printf("lpstat: %d\n",lpstat);*/

			if(keep_solving_infeasible_MIPs) while(lpstat == 108)
			{
				CPXmipopt (env_just_solve_mips, nodelp_copy2);
			 	lpstat = CPXgetstat (env_just_solve_mips, nodelp_copy2);
			}
			
			if(lpstat == 103)
			{
/*				printf("fathoming since the mip is infeasible\n");*/
				fathoming = 1;
				goto TERMINATE;
			}				
			
			int z = 0;
		    	double usercut_coefs[2] = {1.,-1./slope};
		    	int usercut_indices[2] = {obj1_index,obj2_index};
		    	char usercut_sense[1] = {'L'};
		    	double usercut_rhs = 0.;
/*		    	printf("slope: %lf\n",slope);*/
		    	
		    	if(lpstat == 101 || lpstat == 102) 
		    	{
		    		status = CPXgetobjval (env_just_solve_mips, nodelp_copy2, &usercut_rhs);
		    	}
		    	else 
		    	{
		    		status = CPXgetbestobjval (env_just_solve_mips, nodelp_copy2, &usercut_rhs);
		    	}
		    	
		    	double obj_off;
		    	status = CPXgetobjoffset(env_just_solve_mips, nodelp_copy2, &obj_off );
			usercut_rhs -= obj_off;
			
/*						    		double tx = -6295.;*/
/*					double tx2 = -754.;*/
/*				    	double ty = (tx - usercut_rhs)*slope;*/
/*				    	double ty2 = (tx2 - usercut_rhs)*slope;*/
/*				    	*/
/*				    	printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",tx,tx2,ty,ty2);*/
/*				    	exit(0);*/
			
			if(printing_in_setbranch) printf("adding user cut: x%d + %lfx%d <= %lf\n",obj1_index,-1./slope,obj2_index,usercut_rhs);
		    	
		    	status = CPXcutcallbackaddlocal   (env,
							   cbdata,
				       			   wherefrom,
							   2,
							   usercut_rhs,
							   'L',
							   usercut_indices,
							   usercut_coefs);

			int numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, nodelp_copy2);
			int numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, nodelp_copy2);

			num_starts = numsolns - prev_numsols + numrep;

			prev_numsols = numsolns;

		    	if(num_starts > global_num_starts)
		    	{
		    		global_num_starts = num_starts;
		    		global_startspace = cur_numcols*global_num_starts;
		    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
				global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
				global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
				global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
		    	}

			status = CPXgetmipstarts (env_just_solve_mips, nodelp_copy2, &nzcnt, global_beg, global_varindices, 
					   global_values, global_effortlevel, global_startspace,
					   &surplus, 0, num_starts-1);
					   
			status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,
					   global_values, global_effortlevel, NULL);
		  	
		  	finish_mipsolve = clock();
			double duration_mipsolve = (double)(finish_mipsolve- start_mipsolve) / CLOCKS_PER_SEC;
			if(duration_mipsolve > max_time_to_solve_a_mip) max_time_to_solve_a_mip = duration_mipsolve;
			time_solving_mips += duration_mipsolve;
			if(printing_in_setbranch) printf("time to solve mip at seqnum %d: %lf\n",seqnum,duration_mipsolve);
			
			cumulative_time = (double)(finish_mipsolve - start_BB) / CLOCKS_PER_SEC;
			if(cumulative_time > max_time) goto BRANCHING;
	
			int num_solns = CPXgetsolnpoolnumsolns (env_just_solve_mips, nodelp_copy2);
		  	
		  	int insert_check = 0;

			if(num_solns >= prev_numsolns) times_to_run = num_solns - prev_numsolns;
		  	else times_to_run = num_solns;
		  	prev_numsolns = num_solns;
		  	
		  	for(j=0;j<times_to_run;j++)
		  	{
		  		status = CPXgetsolnpoolx (env_just_solve_mips, nodelp_copy2, j, x_ws, 0, cur_numcols-1);
			      	insert_check = mock_insert(1,x_ws[obj1_index],x_ws[obj2_index],0,0,0,&tree);
			      	if(insert_check)
			      	{
			      		if(branch_iterations < 5) for(i=0;i<cur_numcols;i++)
			      		{
			      			stored_x[x_rotation][i] = x_ws[i];
			      		}
			      		x_rotation = (x_rotation + 1) % num_x_to_store;
			      		insert_check = 0;
			      		PSA_full(env_global,NULL,x_ws,NULL,NULL);
			      	}
		      	}
		  	
		  	lpstat = CPXgetstat (env_just_solve_mips, nodelp_copy2);
		  	if(printing_in_setbranch) printf("solve status: %d\n",lpstat);

			if(lpstat == 103 || lpstat == 119)
			{
				if(printing_in_setbranch) printf("infeasible mip. Fathoming\n");
				fathoming = 1;
/*				*useraction_p = CPX_CALLBACK_SET;*/
/*				nodecnt = 0;*/
			  	goto TERMINATE;
			}
			else if(lpstat == 108)
			{
				goto BRANCHING;
				printf("Warning: Fathoming a node as infeasible because no feasible solution was found within time limit %lf s. This may cause incorrect solutions\n",time_limit);
				if(printing_in_setbranch) printf("infeasible mip. Fathoming\n");
				fathoming = 1;
/*				*useraction_p = CPX_CALLBACK_SET;*/
/*				nodecnt = 0;*/
			  	goto TERMINATE;
			}
		  	
		  	if(printing_in_setbranch) printf("getting x\n");
		  	status = CPXgetx (env_just_solve_mips, nodelp_copy2, x_ws, 0, cur_numcols-1);
			if(status) 
			{
				printf ("(%d) CPXgetx, Failed to get x values, error code %d\n", __LINE__,status);
				printf("status of the solve: %d\n",lpstat);
				goto TERMINATE;
			}
			ws_mip_objvals[0] = x_ws[obj1_index];
			ws_mip_objvals[1] = x_ws[obj2_index];
	
			if(!userhandle_up->x_ws) userhandle_up->x_ws = calloc ( (cur_numcols+1),sizeof(double));
			if(!userhandle_down->x_ws) userhandle_down->x_ws = calloc ( (cur_numcols+1),sizeof(double));
	
			for(i=0;i<cur_numcols;i++) 
			{
/*				printf("xws_%d: %lf\n",i,x_ws[i]);*/
				userhandle_up->x_ws[i] = x_ws[i];
				userhandle_down->x_ws[i] = x_ws[i];
			}
			userhandle_up->x_ws[cur_numcols] = slope;
			userhandle_down->x_ws[cur_numcols] = slope;
	
			if(x_ws[obj1_index] - reduced_subprob_x_ub >= -.0000001) 
			{
				right_side_dom = 1;
				if(reduced_subprob_y_lb < x_ws[obj2_index])
				{
					reduced_subprob_y_lb = x_ws[obj2_index];
					bound_reduction = 1;
				}
			}
			if(x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001) 
			{
				left_side_dom = 1;
				if(reduced_subprob_y_ub > x_ws[obj2_index])
				{
					reduced_subprob_y_ub = x_ws[obj2_index];
					bound_reduction = 1;
				}
			}

			if(printing_in_setbranch) printf("plotting2\n");
/*			if(printing_in_setbranch) */
			if(printing_in_setbranch) printf("plot(%lf,%lf,'mo');\n",x_ws[obj1_index],x_ws[obj2_index]);

			if(lpstat == 101 || lpstat == 102)
			{
				if(printing_in_setbranch) printf("the mip solution was optimal\n");
				ws_mip_opt = 1;
				
				if(control_node_selection)
				{
					userhandle_up->f1 = x_ws[obj1_index];
					userhandle_down->f1 = x_ws[obj1_index];
					userhandle_up->f2 = x_ws[obj2_index];
					userhandle_down->f2 = x_ws[obj2_index];
				}
				
				add_check = mock_insert(1,x_ws[obj1_index],x_ws[obj2_index],0,0,0,&tree);
				if(add_check)
				{
					if(printing_in_setbranch) printf("adding solution\n");
					for(i=0;i<cur_numcols;i++)
			      		{
			      			stored_x[x_rotation][i] = x_ws[i];
			      		}
			      		x_rotation = (x_rotation + 1) % num_x_to_store;
			      		add_check = 0;
			      		PSA_full(env,NULL,x_ws,NULL,NULL);
				}
				else if(objective_space_fathoming)
				{
					int add_check_1 = mock_insert(1,x_ws[obj1_index]+.0001*x_range,x_ws[obj2_index],0,0,0,&tree);
					int add_check_2 = mock_insert(1,x_ws[obj1_index],x_ws[obj2_index]+.0001*y_range,0,0,0,&tree);
					if(!add_check_1 && !add_check_2)
					{
						sub_pr1_x_ub = x_ws[obj1_index];
						sub_pr1_y_lb = x_ws[obj2_index]+.0001*y_range;
					 	sub_pr2_x_lb = x_ws[obj1_index]+.0001*x_range;
						sub_pr2_y_ub = x_ws[obj2_index];
						pareto = 1;
						goto TERMINATE;
					}
					else if(!add_check_1)
					{
						sub_pr1_x_ub = x_ws[obj1_index];
						sub_pr1_y_lb = x_ws[obj2_index];
					 	sub_pr2_x_lb = x_ws[obj1_index]+.0001*x_range;
						sub_pr2_y_ub = x_ws[obj2_index];
						pareto = 1;
						goto TERMINATE;
					}
					else if(!add_check_2)
					{
						sub_pr1_x_ub = x_ws[obj1_index];
						sub_pr1_y_lb = x_ws[obj2_index]+.0001*y_range;
					 	sub_pr2_x_lb = x_ws[obj1_index];
						sub_pr2_y_ub = x_ws[obj2_index];
						pareto = 1;
						goto TERMINATE;
					}
				}
		
				if(its_been_only_points)
				{
					sub_pr1_x_ub = x_ws[obj1_index];
					sub_pr1_y_lb = x_ws[obj2_index];
				 	sub_pr2_x_lb = x_ws[obj2_index] + .001;
					sub_pr2_y_ub = x_ws[obj2_index];
					pareto = 1;
					goto TERMINATE;
/*					goto PARETO_BRANCH;*/
				}
		
				if(x_ws[obj2_index] - reduced_subprob_y_ub >= -.0000001) reduced_subprob_x_lb = x_ws[obj1_index];
				if(x_ws[obj2_index] - reduced_subprob_y_lb <= .0000001) reduced_subprob_x_ub = x_ws[obj1_index];
		
				if(!exact_mips && x_ws[obj1_index] > reduced_subprob_x_lb && x_ws[obj2_index] < reduced_subprob_y_ub && 
					x_ws[obj1_index] < reduced_subprob_x_ub && x_ws[obj2_index] > reduced_subprob_y_lb)
				{
					if(printing_in_setbranch) printf("weighted sum solution in interior of search region\n");
			
					if(!its_been_only_points && !points_only)
					{
						if(going_back) goto BRANCHING;
		
						if(!ob1_lp_int_feas && !ob1_mip_ran) goto SOLVE_OB1_MIP;
						else if(ob1_mip_opt) goto RECHECK_IDEALS;
						else goto BRANCHING;
					}
			
					sub_pr1_x_ub = x_ws[obj1_index];
					sub_pr1_y_lb = x_ws[obj2_index];
				 	sub_pr2_x_lb = x_ws[obj1_index] + .001;
					sub_pr2_y_ub = x_ws[obj2_index];
					add_check = mock_insert(1,x_ws[obj1_index],reduced_subprob_y_ub,0,0,0,&tree);
					if(!add_check)
					{
						if(printing_in_setbranch) printf("left side is dominated\n");
						left_side_dom = 1;
					}
					int add_check2 = mock_insert(1,reduced_subprob_x_ub,x_ws[obj2_index],0,0,0,&tree);
					if(!add_check2)
					{
						if(printing_in_setbranch) printf("right side is dominated\n");
						right_side_dom = 1;
					}
					if(left_side_dom && right_side_dom)
					{
						fathomed_by_dominated_local_ideal_pts++;
						fathoming = 1;
/*						*useraction_p = CPX_CALLBACK_SET;*/
/*						nodecnt = 0;*/
					  	goto TERMINATE;
					}
					else if(left_side_dom)
					{
						reduced_subprob_x_lb = x_ws[obj1_index];
						reduced_subprob_y_ub = x_ws[obj2_index];
						bound_reduction = 1;
					  	goto SOLVE_OB1_LP;
					}
					else if(right_side_dom)
					{
						reduced_subprob_x_ub = x_ws[obj1_index];
						reduced_subprob_y_lb = x_ws[obj2_index];
						bound_reduction = 1;
					  	goto SOLVE_OB2_LP;
					}
					else if(!ws_lp_dom && pareto_branching) 
					{
						pareto = 1;
						goto TERMINATE;
/*						goto PARETO_BRANCH;*/
					}
				}
				else if(!exact_mips && !ws_lp_dom && pareto_branching && x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
				{
					if(printing_in_setbranch) printf("weighted sum solution on left boundary of search region\n");
		/*			DO_IT_AGAIN1:*/
		/*			;*/
		/*			closest_nodes *two_nodes = find_two_nodes_right_of_val(reduced_subprob_x_lb, reduced_subprob_y_ub, tree);*/
		/*			if(!two_nodes) */
		/*			{*/
		/*				if(points_only || its_been_only_points) goto SOLVE_OB2_MIP;*/
		/*				else goto AFTER_THIS3;*/
		/*			}*/
		/*			if(printing_in_setbranch) printf("the two nodes:\n");*/
		/*			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->closest->nw_x,*/
		/*				x_ideal-two_nodes->closest->se_x,y_ideal-two_nodes->closest->nw_y,y_ideal-two_nodes->closest->se_y);*/
		/*			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->next->nw_x,x_ideal-two_nodes->next->se_x,*/
		/*				y_ideal-two_nodes->next->nw_y,y_ideal-two_nodes->next->se_y);*/
		/*			if(two_nodes->closest->type == 2 && fabs(two_nodes->closest->nw_x - two_nodes->closest->se_x) < .00000001 && */
		/*				fabs(two_nodes->closest->nw_y - two_nodes->closest->se_y) < .00000001) two_nodes->closest->type = 1;*/
		/*			if(two_nodes->next->type == 2 && fabs(two_nodes->next->nw_x - two_nodes->next->se_x) < .00000001 && */
		/*				fabs(two_nodes->next->nw_y - two_nodes->next->se_y) < .00000001) two_nodes->next->type = 1;*/
		/*			*/
		/*			if(two_nodes->closest->type == 1 && ((fabs(two_nodes->closest->nw_x - two_nodes->next->nw_x) < .0000001 && */
		/*				two_nodes->closest->nw_y - two_nodes->next->nw_y >= -.0000001) || */
		/*				(fabs(two_nodes->closest->nw_x - two_nodes->next->se_x) < .0000001 && */
		/*				two_nodes->closest->nw_y - two_nodes->next->se_y >= -.0000001)))*/
		/*			{*/
		/*				if(printing_in_setbranch) */
		/*					printf("the closest node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
		/*				delete_node(two_nodes->closest);*/
		/*				free(two_nodes);*/
		/*				goto DO_IT_AGAIN1;*/
		/*			}*/
		/*			*/
		/*			if(x_ideal - two_nodes->closest->nw_x - reduced_subprob_x_ub > .0001 || */
		/*				y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_lb < -.0001 || */
		/*				x_ideal - two_nodes->closest->nw_x - reduced_subprob_x_lb < -.0001 || */
		/*				y_ideal - two_nodes->closest->nw_y - reduced_subprob_y_ub > .0001)*/
		/*			{*/
		/*				if(printing_in_setbranch) printf("closest node is outside region\n");*/
		/*				free(two_nodes);*/
		/*				goto SOLVE_OB1_LP;*/
		/*			}*/
		/*			*/
		/*			else if(two_nodes->next->type == 1 && ((fabs(two_nodes->next->nw_x - two_nodes->closest->nw_x) < .0000001 && */
		/*				two_nodes->next->nw_y - two_nodes->closest->nw_y >= -.0000001) || */
		/*				(fabs(two_nodes->next->nw_x - two_nodes->closest->se_x) < .0000001 && */
		/*				two_nodes->next->nw_y - two_nodes->closest->se_y >= -.0000001)))*/
		/*			{*/
		/*				if(printing_in_setbranch) */
		/*					printf("the next node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
		/*				delete_node(two_nodes->next);*/
		/*				free(two_nodes);*/
		/*				goto DO_IT_AGAIN1;*/
		/*			}*/
		/*			if(x_ideal - two_nodes->next->se_x - reduced_subprob_x_ub <= .0000001)*/
		/*			{*/
		/*				if(printing_in_setbranch) printf("left pt of next node is also inside the search region.\n");*/
		/*				if(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y) < -.0001 && */
		/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb >= -.0000001 && */
		/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_ub <= .0000001)*/
		/*				{*/
		/*					if(printing_in_setbranch) printf("separation between y-values, split\n");*/
		/*					printf("y separation is %lf percent of y_range\n", */
		/*							100.*(y_ideal - two_nodes->next->se_y - (y_ideal - two_nodes->closest->nw_y))/y_range);*/
		/*					sub_pr1_x_ub = reduced_subprob_x_ub;*/
		/*					sub_pr1_y_lb = ((y_ideal - two_nodes->next->se_y) + (y_ideal - two_nodes->closest->nw_y))/2.;*/
		/*				 	sub_pr2_x_lb = x_ideal - two_nodes->closest->nw_x + .001;//reduced_subprob_x_lb;*/
		/*					sub_pr2_y_ub = ((y_ideal - two_nodes->next->se_y) + (y_ideal - two_nodes->closest->nw_y))/2.;*/
		/*					free(two_nodes);*/
		/*					goto PARETO_BRANCH;*/
		/*				}*/
		/*				else if(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x) > .0001)*/
		/*				{*/
		/*					if(printing_in_setbranch) printf("separation between x-values, split\n");*/
		/*					printf("y separation is %lf percent of y_range\n", */
		/*							100.*(x_ideal - two_nodes->next->se_x - (x_ideal - two_nodes->closest->nw_x))/x_range);*/
		/*					sub_pr1_x_ub = ((x_ideal - two_nodes->next->se_x) + (x_ideal - two_nodes->closest->nw_x))/2.;*/
		/*					sub_pr1_y_lb = y_ideal - two_nodes->next->se_y + .001;//reduced_subprob_y_lb;*/
		/*				 	sub_pr2_x_lb = ((x_ideal - two_nodes->next->se_x) + (x_ideal - two_nodes->closest->nw_x))/2.;*/
		/*					sub_pr2_y_ub = reduced_subprob_y_ub;*/
		/*					free(two_nodes);*/
		/*					goto PARETO_BRANCH;*/
		/*				}*/
		/*				else if(x_ideal - two_nodes->next->se_x - reduced_subprob_x_ub < -.001 && */
		/*					y_ideal - two_nodes->next->se_y - reduced_subprob_y_lb > .001)*/
		/*				{*/
		/*					if(printing_in_setbranch) printf("left pt of next node is strictly inside the search region.\n");*/
		/*					double ran = (double) rand() / ( (double) RAND_MAX);*/
		/*					if(ran < .5)*/
		/*					{*/
		/*						sub_pr1_x_ub = (x_ideal - two_nodes->closest->nw_x);*/
		/*						sub_pr1_y_lb = (y_ideal - two_nodes->closest->nw_y);//reduced_subprob_y_lb;*/
		/*					 	sub_pr2_x_lb = (x_ideal - two_nodes->closest->nw_x);*/
		/*						sub_pr2_y_ub = reduced_subprob_y_ub;*/
		/*						free(two_nodes);*/
		/*						goto PARETO_BRANCH;*/
		/*					}*/
		/*					else*/
		/*					{*/
		/*						sub_pr1_x_ub = reduced_subprob_x_ub;*/
		/*						sub_pr1_y_lb = (y_ideal - two_nodes->closest->nw_y);*/
		/*					 	sub_pr2_x_lb = (x_ideal - two_nodes->closest->nw_x);//reduced_subprob_x_lb;*/
		/*						sub_pr2_y_ub = (y_ideal - two_nodes->closest->nw_y);*/
		/*						free(two_nodes);*/
		/*						goto PARETO_BRANCH;*/
		/*					}*/
		/*				}*/
		/*			}*/
		/*			AFTER_THIS3:*/
		/*			*/
		/*			free(two_nodes);*/
					goto SOLVE_OB1_LP;
			
				}
				else if(!exact_mips && !ws_lp_dom && pareto_branching && x_ws[obj1_index] - reduced_subprob_x_ub >= -.0000001)
				{
					if(printing_in_setbranch) printf("weighted sum solution on right boundary of search region\n");
		/*			DO_IT_AGAIN2:*/
		/*			;*/
		/*			closest_nodes *two_nodes = find_two_nodes_left_of_val(reduced_subprob_x_ub, reduced_subprob_y_lb, tree);*/
		/*			if(!two_nodes) */
		/*			{*/
		/*				if(points_only || its_been_only_points) goto SOLVE_OB2_MIP;*/
		/*				else goto AFTER_THIS4;*/
		/*			}*/
		/*			if(printing_in_setbranch) printf("the two nodes:\n");*/
		/*			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->closest->nw_x,*/
		/*				x_ideal-two_nodes->closest->se_x,y_ideal-two_nodes->closest->nw_y,y_ideal-two_nodes->closest->se_y);*/
		/*			if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-two_nodes->next->nw_x,x_ideal-two_nodes->next->se_x,*/
		/*				y_ideal-two_nodes->next->nw_y,y_ideal-two_nodes->next->se_y);*/
		/*			*/
		/*			if(two_nodes->closest->type == 2 && fabs(two_nodes->closest->nw_x - two_nodes->closest->se_x) < .00000001 && */
		/*				fabs(two_nodes->closest->nw_y - two_nodes->closest->se_y) < .00000001) two_nodes->closest->type = 1;*/
		/*			if(two_nodes->next->type == 2 && fabs(two_nodes->next->nw_x - two_nodes->next->se_x) < .00000001 && */
		/*				fabs(two_nodes->next->nw_y - two_nodes->next->se_y) < .00000001) two_nodes->next->type = 1;*/
		/*			if(two_nodes->closest->type == 1 && ((fabs(two_nodes->closest->nw_x - two_nodes->next->nw_x) < .0000001 && */
		/*				two_nodes->closest->nw_y - two_nodes->next->nw_y >= -.0000001) || */
		/*				(fabs(two_nodes->closest->nw_x - two_nodes->next->se_x) < .0000001 && */
		/*				two_nodes->closest->nw_y - two_nodes->next->se_y >= -.0000001)))*/
		/*			{*/
		/*				if(printing_in_setbranch) */
		/*					printf("the closest node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
		/*				delete_node(two_nodes->closest);*/
		/*				free(two_nodes);*/
		/*				goto DO_IT_AGAIN2;*/
		/*			}*/
		/*			else if(two_nodes->next->type == 1 && ((fabs(two_nodes->next->nw_x - two_nodes->closest->nw_x) < .0000001 && */
		/*				two_nodes->next->nw_y - two_nodes->closest->nw_y >= -.0000001) || */
		/*				(fabs(two_nodes->next->nw_x - two_nodes->closest->se_x) < .0000001 && */
		/*				two_nodes->next->nw_y - two_nodes->closest->se_y >= -.0000001)))*/
		/*			{*/
		/*				if(printing_in_setbranch) */
		/*					printf("the next node is dominated by, or is a repeat of, one of the endpoints of the next node\n");*/
		/*				delete_node(two_nodes->next);*/
		/*				free(two_nodes);*/
		/*				goto DO_IT_AGAIN2;*/
		/*			}*/
		/*			*/
		/*			if(x_ideal - two_nodes->closest->se_x - reduced_subprob_x_ub > .0001 || */
		/*				y_ideal - two_nodes->closest->se_y - reduced_subprob_y_lb < -.0001 || */
		/*				x_ideal - two_nodes->closest->se_x - reduced_subprob_x_lb < -.0001 || */
		/*				y_ideal - two_nodes->closest->se_y - reduced_subprob_y_ub > .0001)*/
		/*			{*/
		/*				if(printing_in_setbranch) printf("closest node is outside region\n");*/
		/*				free(two_nodes);*/
		/*				goto SOLVE_OB2_LP;*/
		/*			}*/
		/*			*/
		/*			if(x_ideal - two_nodes->next->nw_x - reduced_subprob_x_lb >= -.0000001 && */
		/*				y_ideal - two_nodes->next->nw_y - reduced_subprob_y_ub <= .0000001 )*/
		/*			{*/
		/*				if(printing_in_setbranch) printf("right pt of next node is also inside the search region.\n");*/
		/*				if(x_ideal - two_nodes->closest->se_x - (x_ideal - two_nodes->next->nw_x) > .0001)*/
		/*				{*/
		/*					sub_pr1_x_ub = ((x_ideal - two_nodes->closest->se_x) + (x_ideal - two_nodes->next->nw_x))/2.;*/
		/*					sub_pr1_y_lb = y_ideal - two_nodes->closest->se_y;//reduced_subprob_y_lb;*/
		/*				 	sub_pr2_x_lb = ((x_ideal - two_nodes->closest->se_x) + (x_ideal - two_nodes->next->nw_x))/2.;*/
		/*					sub_pr2_y_ub = reduced_subprob_y_ub;*/
		/*					free(two_nodes);*/
		/*					goto PARETO_BRANCH;*/
		/*				}*/
		/*				else if(y_ideal - two_nodes->closest->se_y - (y_ideal - two_nodes->next->nw_y) < -.0001 && */
		/*					y_ideal - two_nodes->closest->se_y - reduced_subprob_y_ub <= .0000001 && */
		/*					y_ideal - two_nodes->closest->se_y - reduced_subprob_y_lb >= -.0000001)*/
		/*				{*/
		/*					sub_pr1_x_ub = reduced_subprob_x_ub;*/
		/*					sub_pr1_y_lb = ((y_ideal - two_nodes->closest->se_y) + (y_ideal - two_nodes->next->nw_y))/2.;*/
		/*				 	sub_pr2_x_lb = x_ideal - two_nodes->next->nw_x;//reduced_subprob_x_lb;*/
		/*					sub_pr2_y_ub = ((y_ideal - two_nodes->closest->se_y) + (y_ideal - two_nodes->next->nw_y))/2.;*/
		/*					free(two_nodes);*/
		/*					goto PARETO_BRANCH;*/
		/*				}*/
		/*				else if(x_ideal - two_nodes->next->nw_x - reduced_subprob_x_lb > .001 && */
		/*					y_ideal - two_nodes->next->nw_y - reduced_subprob_y_ub < -.001)*/
		/*				{*/
		/*					double ran = (double) rand() / ( (double) RAND_MAX);*/
		/*					if(ran < .5)*/
		/*					{*/
		/*						sub_pr1_x_ub = (x_ideal - two_nodes->next->se_x);*/
		/*						sub_pr1_y_lb = (y_ideal - two_nodes->next->se_y);//reduced_subprob_y_lb;*/
		/*					 	sub_pr2_x_lb = (x_ideal - two_nodes->next->se_x);*/
		/*						sub_pr2_y_ub = reduced_subprob_y_ub;*/
		/*						free(two_nodes);*/
		/*						goto PARETO_BRANCH;*/
		/*					}*/
		/*					else*/
		/*					{*/
		/*						sub_pr1_x_ub = reduced_subprob_x_ub;*/
		/*						sub_pr1_y_lb = (y_ideal - two_nodes->next->se_y);*/
		/*					 	sub_pr2_x_lb = (x_ideal - two_nodes->next->se_x);//reduced_subprob_x_lb;*/
		/*						sub_pr2_y_ub = (y_ideal - two_nodes->next->se_y);*/
		/*						free(two_nodes);*/
		/*						goto PARETO_BRANCH;*/
		/*					}*/
		/*				}*/
		/*			}*/
		/*			AFTER_THIS4:*/
		/*			*/
		/*			free(two_nodes);*/
					goto SOLVE_OB2_LP;
				}
		
				if(going_back) goto BRANCHING;
		
				if(!ob1_lp_been_solved) goto SOLVE_OB1_LP;
				else if(!ob1_lp_int_feas && !ob1_mip_ran) goto SOLVE_OB1_MIP;
				else if(ob1_mip_opt) goto RECHECK_IDEALS;
				else goto BRANCHING;
			}
			else if(lpstat == 105 || lpstat == 107 || lpstat == 109 || lpstat == 111 || lpstat == 113 || lpstat == 116 || lpstat == 106)
			{
				if(printing_in_setbranch) printf("the mip solution was less than optimal\n");
				ws_mip_opt = 0;
				add_check = mock_insert(1,x_ws[obj1_index],x_ws[obj2_index],0,0,0,&tree);
				if(add_check)
				{
					if(printing_in_setbranch) printf("adding solution\n");
					for(i=0;i<cur_numcols;i++)
			      		{
			      			stored_x[x_rotation][i] = x_ws[i];
			      		}
			      		x_rotation = (x_rotation + 1) % num_x_to_store;
			      		add_check = 0;
			      		PSA_full(env,NULL,x_ws,NULL,NULL);
				}
		
				if(going_back) goto BRANCHING;
		
				status = CPXgetbestobjval (env_just_solve_mips, nodelp_copy2, &best_bound);
				projection = -ws_lp_objvals[1]+best_bound;
				if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",projection,ws_lp_objvals[1]);
		/*		printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",-65.,best_bound-5./(-1.),(-1.)*(-best_bound-65),-5.);*/
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",-65.,best_bound-5./(slope),(slope)*(-best_bound-65),-5.);
	
				add_check = mock_insert(1,projection,ws_lp_objvals[1],0,0,0,&tree);
				if(add_check)
				{
					if(printing_in_setbranch) printf("the projection onto the best objval level curve is not dominated, go to branching\n");
					goto BRANCHING;
				}
				else
				{
					if(printing_in_setbranch) printf("the projection onto the best objval level curve is dominated\n");
					if(printing_in_setbranch) printf("trying this: setting obj vals to be this projected point and checking domination\n");
					x_ws[obj1_index] = projection;
					x_ws[obj2_index] = ws_lp_objvals[1];
					if(!ob1_lp_been_solved) goto SOLVE_OB1_LP;
					else if(!ob1_lp_int_feas) goto SOLVE_OB1_MIP;
					else goto BRANCHING;
				}
			}
			else
			{
				printf("the status of the FAILED mipsolve: %d\n",lpstat);
				exit(0);
			}
	
			/*********** Here we solve and process the single objective LP solution associated with f_1 ************/
	
			SOLVE_OB1_LP:
	
		/*	printf("arrived at solve of ob1 lp\n");*/

			if(userhandle_current && userhandle_current->ob1_still_feas) goto SOLVE_OB1_MIP;
	
			if(ws_lp_sol_on_right_bd) 
			{
				ob1_lp_been_solved = 1;
				goto SOLVE_OB2_LP;
			}
	
			lp_ob1 = CPXcloneprob (env, nodelp2, &status);
		  	if ( status ) {
		    		printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
				goto TERMINATE;
			}
	
			chg_coefs(env, lp_ob1, indexes,-10000000.);
	
		  	ob1_lp_been_solved = 1;
		  	status = CPXlpopt (env, lp_ob1);
		 	if ( status ) {
		   		printf ("%s(%d): CPXlpopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
				goto TERMINATE;
			}
	
		  	status = CPXgetx (env, lp_ob1, x1, 0, cur_numcols-1);
			if(status) 
			{
				printf ("(%d) CPXgetx, Failed to get x values, error code %d\n", __LINE__,status);
				goto TERMINATE;
			}
			ob1_lp_objvals[0] = x1[obj1_index];
			ob1_lp_objvals[1] = x1[obj2_index];
	
		/*	if(x1[obj2_index] > x_ws[obj2_index]) goto SOLVE_WS_MIP;*/
	
			if(ob1_lp_objvals[0] < reduced_subprob_x_ub)
			{
				reduced_subprob_x_ub = ob1_lp_objvals[0];
		/*		reduced_subprob_y_lb = ob1_lp_objvals[1];*/
				if(printing_in_setbranch) printf("after changing:\n");
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_ub,reduced_subprob_y_ub);
				bound_reduction = 1;
			}

			if(printing_in_setbranch) printf("plotting ob1 lp soln\n");
			if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",x1[obj1_index],x1[obj2_index]);
	
			all_feas = 1;
			for(i=0;i<total_num_integer;i++)
			{
				k = integer_indices[i];
			  	diff = x1[k] - floor(x1[k]);
		/*		printf("diff%d: %lf\n",i,diff);*/
		  		if( diff >= .00001 && diff <= .99999)
		  		{
		  			if(printing_in_setbranch) printf("changing frac index to %d (%d)\n",integer_indices[i],__LINE__);
		/*  			printf("plot(%lf,%lf,'o');\n",x1[obj1_index],x1[obj2_index]);*/
			  		frac_index = k;
			  		frac_val = x1[k];
			  		frac_values[k] = frac_val;
			  		if(frac_scores[i] > 0.0001)
			  		{
				  		frac_scores[i] += 1.;
		/*		  		num_frac++;*/
				  	}
				  	else
				  	{
				  		frac_scores[i] += multiplier*k;
			  			num_frac++;
				  	}
		/*		  	printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
			  		all_feas = 0;
		/*	  		printf("num_frac: %d\n",num_frac);*/
		/*			printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
				}
			}
			if(all_feas == 1)
			{
				if(printing_in_setbranch) printf("ob1 lp solution is integer feasible\n");
				all_feas = 0;
		
				ob1_mip_opt = 1;
				ob1_mip_objvals[0] = x1[obj1_index];
				ob1_mip_objvals[1] = x1[obj2_index];
				
				if(!userhandle_up->x1) userhandle_up->x1 = calloc ((cur_numcols),sizeof(double));
				if(!userhandle_down->x1) userhandle_down->x1 = calloc ((cur_numcols),sizeof(double));
	
				for(i=0;i<cur_numcols;i++) 
				{
					userhandle_up->x1[i] = x1[i];
					userhandle_down->x1[i] = x1[i];
				}
		
				ob1_lp_int_feas = 1;
			}
		
			add_check = mock_insert(1,x1[obj1_index],x1[obj2_index],0,0,0,&tree);
			if(add_check)
			{
				if(printing_in_setbranch) printf("the ob1 lp solution is not dominated\n");
				if(ob1_lp_int_feas == 1)
				{
					add_check = mock_insert(1,x1[obj1_index],x1[obj2_index],0,0,0,&tree);
					if(add_check)
					{
						if(printing_in_setbranch) printf("adding solution\n");
						for(i=0;i<cur_numcols;i++)
				      		{
					      		stored_x[x_rotation][i] = x1[i];
					      	}
					      	x_rotation = (x_rotation + 1) % num_x_to_store;
				 		add_check = 0;
				     		if(check_for_stopping_PSA_full)
						{
							check_for_stopping_PSA_full = 0;
							PSA_full(env,NULL,x1,NULL,NULL);
							check_for_stopping_PSA_full = 1;
						}
						else PSA_full(env,NULL,x1,NULL,NULL);
					}
			
					add_check = mock_insert(1,x1[obj1_index],x_ws[obj2_index],0,0,0,&tree);
					if(add_check)
					{
						if(printing_in_setbranch) printf("the right partial ideal is not dominated\n");
		/*				projection = (-1.)*(x1[obj1_index]-x_ws[obj1_index])+x_ws[obj2_index];*/
						projection = (slope)*(x1[obj1_index]-x_ws[obj1_index])+x_ws[obj2_index];
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],x1[obj1_index],x_ws[obj2_index],projection);
		/*				add_check = mock_insert(2,x1[obj1_index],projection,x_ws[obj1_index],x_ws[obj2_index],-1.,&tree);*/
		/*				printf("slope used: %lf\n",slope);*/
						add_check = mock_insert(2,x1[obj1_index],projection,x_ws[obj1_index],x_ws[obj2_index],slope,&tree);
						if(add_check)
						{
							endpoint1_x = x1[obj1_index];
							endpoint1_y = projection;
							if(printing_in_setbranch) printf("the right partial ideal segment is not dominated (%d)\n",__LINE__);
							if(ws_mip_ran)
							{
								if(printing_in_setbranch) printf("the ws lp soln was not dominated, so PSA will not produce a dominated dual bound.\n");
								goto BRANCHING;
							}
							else
							{
								if(printing_in_setbranch) printf("ws lp soln was dominated. Consider using PSA here to check for a dominated dual bound.\n");
								status = CPXgetbase (env, lp_ob1, cstat, rstat);
								if(printing_in_setbranch) printf("calling PSA\n");
								int reduce_val = PSA(env, lp_ob1, x1, cstat, rstat, seqnum, 0, lp_ob1);
				  				if(reduce_val == 2)
				  				{
				  					if(printing_in_setbranch) printf("fathoming node for completed PSA\n");
								  	fathomed_by_PSA_completion++;
									fathoming = 1;
/*									*useraction_p = CPX_CALLBACK_SET;*/
/*									nodecnt = 0;*/
								  	goto TERMINATE;
				  				}
								else 
								{
					  				if(reduce_val == 1)
							  		{
							  			reduced_subprob_x_ub = sub_pr1_x_ub;
						  				reduced_subprob_y_lb = sub_pr1_y_lb;
						  				if(printing_in_setbranch) printf("after changing:\n");
										if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
											reduced_subprob_x_lb,reduced_subprob_x_lb,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
										if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
											reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
										if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
											reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_lb);					
										if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
											reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
						  				bound_reduction = 1;
							  		}
									if(printing_in_setbranch) printf("PSA did not complete. Try solving MIPs.\n");
									if(!ws_lp_int_feas) goto SOLVE_WS_MIP;
									else goto BRANCHING;
								}
							}
						}
						else
						{
							if(printing_in_setbranch) printf("the right partial ideal segment is dominated. We can begin exploring ob2\n");
							right_pt_dom = 1;
							reduced_subprob_x_ub = x_ws[obj1_index];
							reduced_subprob_y_lb = x_ws[obj2_index];
							if(printing_in_setbranch) printf("after changing:\n");
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
												reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
												reduced_subprob_y_lb,reduced_subprob_y_lb);					
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
												reduced_subprob_y_ub,reduced_subprob_y_ub);
							bound_reduction = 1;
					
							if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
							{
								fathomed_by_dominated_local_ideal_segments++;
								if(printing_in_setbranch) printf("fathoming because left side already dominated\n");
								fathoming = 1;
/*								*useraction_p = CPX_CALLBACK_SET;*/
/*								nodecnt = 0;*/
							  	goto TERMINATE;
							}
					
							if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
							else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
							else goto BRANCHING;
						}
					}
					else
					{
						if(printing_in_setbranch) printf("the right partial ideal is dominated. We can begin exploring ob2\n");
						right_seg_dom = 1;
						reduced_subprob_x_ub = x_ws[obj1_index];
						reduced_subprob_y_lb = x_ws[obj2_index];
						if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
						bound_reduction = 1;
				
						if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
						{
							fathomed_by_dominated_local_ideal_pts++;
							fathoming = 1;
							if(printing_in_setbranch) printf("fathoming because left side already dominated\n");
/*							*useraction_p = CPX_CALLBACK_SET;*/
/*							nodecnt = 0;*/
						  	goto TERMINATE;
						}
				
						if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
						else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
						else goto BRANCHING;
					}
				}
				else
				{
/*					printf("left_side_dom: %d, ws_lp_dom: %d, remove_dominated_middle: %d\n",left_side_dom,ws_lp_dom,remove_dominated_middle);*/
					if(!left_side_dom && ws_lp_dom && remove_dominated_middle)
					{
		/*				printf("running this 1\n");*/
						if(!ws_sol_interior) 
						{
							status = CPXgetbase (env, nodelp_copy, cstat_ws, rstat_ws);
							if(status) 
							{
								printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
								goto TERMINATE;
							}
						}
						else 
						{
							CPXgetbase (env, nodelp, cstat_ws, rstat_ws);
							if(status) 
							{
								printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
								goto TERMINATE;
							}
						}
						if(printing_in_setbranch) printf("status: %d\n",status);
		
						sub_pr1_x_ub = x_ws[obj1_index];
						sub_pr1_y_lb = x_ws[obj2_index];
						sub_pr2_x_lb = x_ws[obj1_index];
						sub_pr2_y_ub = x_ws[obj2_index];
					  	PSA_left_check = PSA_reduce_left(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, indexes);
					  	
						if(PSA_left_check == 2)
						{
				  			if(printing_in_setbranch) printf("the right subproblem is now empty by PSA completion\n");
							reduced_subprob_x_ub = x_ws[obj1_index];
				  			reduced_subprob_y_lb = x_ws[obj2_index];
				  			if(printing_in_setbranch) printf("after changing:\n");
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
												reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
												reduced_subprob_y_lb,reduced_subprob_y_lb);					
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
												reduced_subprob_y_ub,reduced_subprob_y_ub);
				  			bound_reduction = 1;
				  			right_side_dom = 1;
					  		goto SOLVE_OB2_LP;
					  	}
					  	else if(pareto_branching && sub_pr2_x_lb > sub_pr1_x_ub)
						{
							pareto = 1;
							goto TERMINATE;
/*					  		goto PARETO_BRANCH;*/
					  	}
					}
			
					if(!ob1_mip_ran && !going_back_for_lp1) goto SOLVE_OB1_MIP;
					else if(!ob2_lp_been_solved) 
					{
						if(!ws_lp_sol_on_top_bd) goto SOLVE_OB2_LP;
						else
						{
							ob2_lp_been_solved = 1;
							goto BRANCHING;
						}
					}
					else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
					else goto BRANCHING;
				}
			}
			else
			{
				if(printing_in_setbranch) printf("the ob1 lp solution is dominated\n");
				ob1_lp_dom = 1;
		
				add_check = mock_insert(1,x1[obj1_index],x_ws[obj2_index],0,0,0,&tree);
				if(add_check)
				{
					if(printing_in_setbranch) printf("the right partial ideal is not dominated\n");
		/*			projection = (-1.)*(x1[obj1_index]-x_ws[obj1_index])+x_ws[obj2_index];*/
					projection = (slope)*(x1[obj1_index]-x_ws[obj1_index])+x_ws[obj2_index];
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],x1[obj1_index],x_ws[obj2_index],projection);
		/*			add_check = mock_insert(2,x1[obj1_index],projection,x_ws[obj1_index],x_ws[obj2_index],-1.,&tree);*/
		/*			printf("slope used: %lf\n",slope);*/
					add_check = mock_insert(2,x1[obj1_index],projection,x_ws[obj1_index],x_ws[obj2_index],slope,&tree);
					if(add_check)
					{
						endpoint1_x = x1[obj1_index];
						endpoint1_y = projection;
						if(printing_in_setbranch) printf("the right partial ideal segment is not dominated\n");
						if(printing_in_setbranch) printf("consider solving mip(s) to reduce bounds\n"); 
/*						if(ws_mip_ran) */
/*						{*/
/*							if(printing_in_setbranch) printf("the ws lp soln was not dominated, PSA will not produce dominated dual bound.\n");*/
/*							if(!ob1_mip_ran && !ob1_lp_int_feas && !going_back_for_lp1) goto SOLVE_OB1_MIP;*/
/*							else goto BRANCHING;*/
/*						}*/
/*						else*/
						{
							status = CPXgetbase (env, lp_ob1, cstat, rstat);
							if(ws_lp_int_feas && ob1_lp_int_feas)
							{
								int reduce_val = PSA(env, lp_ob1, x1, cstat, rstat, seqnum, 0, lp_ob1);
				  				if(reduce_val == 2)
				  				{
				  					if(printing_in_setbranch) printf("fathoming node for completed PSA\n");
								  	fathomed_by_PSA_completion++;
									fathoming = 1;
/*									*useraction_p = CPX_CALLBACK_SET;*/
/*									nodecnt = 0;*/
								  	goto TERMINATE;
				  				}
				  				else if(reduce_val == 1)
						  		{
						  			reduced_subprob_x_ub = sub_pr1_x_ub;
					  				reduced_subprob_y_lb = sub_pr1_y_lb;
					  				if(printing_in_setbranch) printf("after changing:\n");
									if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
										reduced_subprob_x_lb,reduced_subprob_x_lb,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
									if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
										reduced_subprob_x_ub,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
									if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
										reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
									if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
										reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_ub,reduced_subprob_y_ub);
					  				bound_reduction = 1;
						  		}
						  		goto BRANCHING;
							}
							else //if(ws_lp_dom)
							{
								if(printing_in_setbranch) printf("calling PSA reduce right\n");
								int PSA_reduce_right_val = PSA_reduce_right(env, lp_ob1, x1, cstat, rstat, indexes, seqnum);
				  				if(PSA_reduce_right_val == 2)
				  				{
				  					if(printing_in_setbranch) printf("fathoming node for dominated lower bound\n");
								  	fathomed_by_dominated_lb++;
									fathoming = 1;
/*									*useraction_p = CPX_CALLBACK_SET;*/
/*									nodecnt = 0;*/
								  	goto TERMINATE;
				  				}
							
								if(PSA_reduce_right_val == 1)
				  				{
				  					reduced_subprob_x_ub = sub_pr1_x_ub;
					  				reduced_subprob_y_lb = sub_pr1_y_lb;
					  				if(printing_in_setbranch) printf("after changing:\n");
									if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
										reduced_subprob_x_lb,reduced_subprob_x_lb,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
									if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
										reduced_subprob_x_ub,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
									if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
										reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
									if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",
										reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_ub,reduced_subprob_y_ub);
					  				bound_reduction = 1;
					  				if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
				  				}
				  				
								if(printing_in_setbranch) printf("lp dual bound not dominated. Try solving MIPs.\n");
								if(!ws_mip_ran) goto SOLVE_WS_MIP;
								else goto SOLVE_OB1_MIP;
		/*						goto BRANCHING;*/
							}
							goto BRANCHING;
						}
					}
					else
					{
						if(printing_in_setbranch) printf("the right partial ideal segment is dominated. Move to consideration of ob2.\n");
						right_pt_dom = 1;
						reduced_subprob_x_ub = x_ws[obj1_index];
						reduced_subprob_y_lb = x_ws[obj2_index];
						if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
						bound_reduction = 1;
				
						if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
						{
							fathomed_by_dominated_local_ideal_segments++;
							fathoming = 1;
							if(printing_in_setbranch) printf("fathoming because left side already dominated\n");
/*							*useraction_p = CPX_CALLBACK_SET;*/
/*							nodecnt = 0;*/
						  	goto TERMINATE;
						}
				
						if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
						else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
						else goto BRANCHING;
					}
				}
				else
				{
					if(printing_in_setbranch) printf("the right partial ideal is dominated. Move to consideration of ob2.\n");
					right_seg_dom = 1;
					reduced_subprob_x_ub = x_ws[obj1_index];
					reduced_subprob_y_lb = x_ws[obj2_index];
					if(printing_in_setbranch) printf("after changing:\n");
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_ub,reduced_subprob_y_ub);
					bound_reduction = 1;
			
					if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
					{
						fathomed_by_dominated_local_ideal_pts++;
						fathoming = 1;
						if(printing_in_setbranch) printf("fathoming because left side already dominated\n");
/*						*useraction_p = CPX_CALLBACK_SET;*/
/*						nodecnt = 0;*/
					  	goto TERMINATE;
					}
			
					if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
					else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
					else goto BRANCHING;
				}
			}
	
			/************** Here we solve the single objective MIP associated with f_1 if we determine that its needed *********/	
	
			SOLVE_OB1_MIP:
	
			if(right_side_dom) goto SOLVE_OB2_LP;
	
			ob1_mip_ran = 1;
		/*	if(exact_mips) goto SKIP_THIS5;*/
			if(userhandle_current && userhandle_current->ob1_still_feas && seqnum != 0)
			{
				if(printing_in_setbranch) printf("obj1 soln from parent node is still feasible\n");
				for(j=0;j<cur_numcols;j++) x1[j] = userhandle_current->x1[j];
				ob1_mip_objvals[0] = x1[obj1_index];
				ob1_mip_objvals[1] = x1[obj2_index];
				
				if(!userhandle_up->x1) userhandle_up->x1 = calloc (cur_numcols,sizeof(double));
				if(!userhandle_down->x1) userhandle_down->x1 = calloc (cur_numcols,sizeof(double));
	
				for(i=0;i<cur_numcols;i++) 
				{
					userhandle_up->x1[i] = x1[i];
					userhandle_down->x1[i] = x1[i];
				}
			
				if(printing_in_setbranch) printf("plotting\n");
				if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",x1[obj1_index],x1[obj2_index]);
		
				if(printing_in_setbranch) printf("therefore ob1 mip solution still optimal\n");
				ob1_mip_opt = 1;
				mip_solved = 1;
		
				if(going_back) goto BRANCHING;
		
				add_check = mock_insert(1,x1[obj1_index],x_ws[obj2_index],0,0,0,&tree);
				if(add_check)
				{
					if(printing_in_setbranch) printf("the right partial ideal is not dominated\n");
		/*			projection = (-1.)*(x1[obj1_index]-x_ws[obj1_index])+x_ws[obj2_index];*/
					projection = (slope)*(x1[obj1_index]-x_ws[obj1_index])+x_ws[obj2_index];
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],x1[obj1_index],x_ws[obj2_index],projection);
		/*			printf("slope used: %lf\n",slope);*/
		/*			add_check = mock_insert(2,x1[obj1_index],projection,x_ws[obj1_index],x_ws[obj2_index],-1.,&tree);*/
					add_check = mock_insert(2,x1[obj1_index],projection,x_ws[obj1_index],x_ws[obj2_index],slope,&tree);
					if(add_check)
					{
						if(printing_in_setbranch) printf("the right partial ideal segment is not dominated\n");
						if(ws_mip_ran || ws_lp_int_feas || left_side_dom || right_side_dom) goto BRANCHING;
						else
						{
							if(printing_in_setbranch) printf("since here we haven't actually solved the ws mip, go back and solve ws mip\n");
							goto SOLVE_WS_MIP;
						}
					}
					else
					{
						if(printing_in_setbranch) printf("the right partial ideal segment is dominated. Move to consideration of ob2.\n");
						right_pt_dom = 1;
						reduced_subprob_x_ub = x_ws[obj1_index];
						reduced_subprob_y_lb = x_ws[obj2_index];
						if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
						bound_reduction = 1;
				
						if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
						{
							fathomed_by_dominated_local_ideal_segments++;
							if(printing_in_setbranch) printf("fathoming because left side already dominated\n");
							fathoming = 1;
/*							*useraction_p = CPX_CALLBACK_SET;*/
/*							nodecnt = 0;*/
						  	goto TERMINATE;
						}
				
						if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
						else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
						else goto BRANCHING;
					}
				}
				else
				{
					if(printing_in_setbranch) printf("the right partial ideal is dominated. Move to consideration of ob2.\n");
					right_seg_dom = 1;
					reduced_subprob_x_ub = x_ws[obj1_index];
					reduced_subprob_y_lb = x_ws[obj2_index];
					if(printing_in_setbranch) printf("after changing:\n");
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_ub,reduced_subprob_y_ub);
					bound_reduction = 1;
			
					if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
					{
						fathomed_by_dominated_local_ideal_pts++;
						fathoming = 1;
						if(printing_in_setbranch) printf("fathoming because left side already dominated\n");
/*						*useraction_p = CPX_CALLBACK_SET;*/
/*						nodecnt = 0;*/
					  	goto TERMINATE;
					}
			
					if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
					else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
					else goto BRANCHING;
				}
			}
		
			if( ws_mip_ran && (!userhandle_current || !(userhandle_current->ws_still_feas)))
			{
				if(printing_in_setbranch) printf("copying mip1 from ws mip\n");
				if(bound_reduction)
				{
					if(printing_in_setbranch) printf("changing bounds\n");
					int ind4[4] = {obj1_index,obj1_index,obj2_index,obj2_index};
					char lu4[4] = {'L','U','L','U'};
					double bds4[4] = {reduced_subprob_x_lb,reduced_subprob_x_ub,reduced_subprob_y_lb,reduced_subprob_y_ub};
					status = CPXchgbds (env_just_solve_mips, nodelp_copy2, 4, ind4, lu4, bds4);
					if ( status ) {
				   		printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
				  	}
				}
				mip_ob1 = CPXcloneprob (env_just_solve_mips, nodelp_copy2, &status);
				if ( status ) {
		    			if(printing_in_setbranch) printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
				    	goto TERMINATE;
		  		}
		  	}
		  	else
		  	{
		  		if(printing_in_setbranch) printf("copying mip1 from nodelp2\n");
		  		if(bound_reduction)
				{
					if(printing_in_setbranch) printf("changing bounds\n");
					int ind4[4] = {obj1_index,obj1_index,obj2_index,obj2_index};
					char lu4[4] = {'L','U','L','U'};
					double bds4[4] = {reduced_subprob_x_lb,reduced_subprob_x_ub,reduced_subprob_y_lb,reduced_subprob_y_ub};
					status = CPXchgbds (env, nodelp2, 4, ind4, lu4, bds4);
					if ( status ) {
				   		printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
				  	}
				}
				mip_ob1 = CPXcloneprob (env_just_solve_mips, nodelp2, &status);
				if ( status ) {
		    			if(printing_in_setbranch) printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
				    	goto TERMINATE;
		  		}
		  		CPXchgprobtype(env_just_solve_mips, mip_ob1, CPXPROB_MILP);
		  		status = CPXchgctype(env_just_solve_mips, mip_ob1, cur_numcols, indexes, xctype);
		  	}
		  	
		  	chg_coefs(env_just_solve_mips, mip_ob1, indexes, -10000000.);
	
			if(!mip_solved)
			{
				if(printing_in_setbranch) printf("here: %d\n",__LINE__);
				mip_solved = 1;
				num_nodes_with_mips_solved++;
			}

			if(printing_in_setbranch) printf("optimizing mip\n");
			start_mipsolve = clock();
			nzcnt = 0;
			prev_numsols = 0;
	
			num_starts = CPXgetnummipstarts (env_just_solve_mips, global_mip);
					    	
		    	if(num_starts > global_num_starts)
		    	{
		    		global_num_starts = num_starts;
		    		global_startspace = cur_numcols*global_num_starts;
		    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
				global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
				global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
				global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
		    	}

			status = CPXgetmipstarts (env_just_solve_mips, global_mip, &nzcnt, global_beg, global_varindices, 
					   global_values, global_effortlevel, global_startspace,
					   &surplus, 0, num_starts-1);
					   
			status = CPXaddmipstarts (env_just_solve_mips, mip_ob1, num_starts, nzcnt, global_beg, global_varindices,
					   global_values, global_effortlevel, NULL);

			CPXmipopt (env_just_solve_mips, mip_ob1);

			numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, mip_ob1);
			numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, mip_ob1);

			num_starts = numsolns - prev_numsols + numrep;

			prev_numsols = numsolns;

		    	if(num_starts > global_num_starts)
		    	{
		    		global_num_starts = num_starts;
		    		global_startspace = cur_numcols*global_num_starts;
		    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
				global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
				global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
				global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
		    	}

			status = CPXgetmipstarts (env_just_solve_mips, mip_ob1, &nzcnt, global_beg, global_varindices, 
					   global_values, global_effortlevel, global_startspace,
					   &surplus, 0, num_starts-1);
					   
			status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,
					   global_values, global_effortlevel, NULL);
					  	
			finish_mipsolve = clock();
			duration_mipsolve = (double)(finish_mipsolve- start_mipsolve) / CLOCKS_PER_SEC;
			if(duration_mipsolve > max_time_to_solve_a_mip) max_time_to_solve_a_mip = duration_mipsolve;
			time_solving_mips += duration_mipsolve;
			if(printing_in_setbranch) printf("time to solve ob1 mip at seqnum %d: %lf\n",seqnum,duration_mipsolve);
			
			cumulative_time = (double)(finish_mipsolve - start_BB) / CLOCKS_PER_SEC;
			if(cumulative_time > max_time) goto BRANCHING;
	
			num_solns = CPXgetsolnpoolnumsolns (env_just_solve_mips, mip_ob1);
		  	
		  	insert_check = 0;
		  	
		  	if(num_solns >= prev_numsolns) times_to_run = num_solns - prev_numsolns;
		  	else times_to_run = num_solns;
		  	prev_numsolns = num_solns;
		  	
		  	for(j=0;j<times_to_run;j++)
		  	{
		  		status = CPXgetsolnpoolx (env_just_solve_mips, mip_ob1, j, x1, 0, cur_numcols-1);
			      	insert_check = mock_insert(1,x1[obj1_index],x1[obj2_index],0,0,0,&tree);
			/*      printf("inserting (%lf,%lf)\n",SE_extreme_x,SE_extreme_y);*/
			/*      printf("insert check: %d\n",insert_check);*/
			      	if(insert_check)
			      	{
			/*      	printf("insert was succesful 1\n");*/
			      		//stored_x[x_rotation] = x;
			/*      	printf("x_rotation: %d\n",x_rotation);*/
			      		if(branch_iterations < 5) for(i=0;i<cur_numcols;i++)
			      		{
			      			stored_x[x_rotation][i] = x1[i];
			/*      		printf("stored val: %lf\n",stored_x[x_rotation][i]);*/
			      		}
			      		x_rotation = (x_rotation + 1) % num_x_to_store;
			      		insert_check = 0;
			      		if(check_for_stopping_PSA_full)
					{
						check_for_stopping_PSA_full = 0;
						PSA_full(env,NULL,x1,NULL,NULL);
						check_for_stopping_PSA_full = 1;
					}
					else PSA_full(env,NULL,x1,NULL,NULL);
			      	}
		      	}
					  	
		  	lpstat = CPXgetstat (env_just_solve_mips, mip_ob1);
			if(printing_in_setbranch) printf("solve status: %d\n",lpstat);
			
			if(keep_solving_infeasible_MIPs) while(lpstat == 108)
			{
				CPXmipopt (env_just_solve_mips, mip_ob1);
			 	lpstat = CPXgetstat (env_just_solve_mips, mip_ob1);
			}
	
			if(lpstat == 103 || lpstat == 119)
			{
				if(printing_in_setbranch) printf("infeasible mip. Fathoming\n");
				fathoming = 1;
/*				*useraction_p = CPX_CALLBACK_SET;*/
/*				nodecnt = 0;*/
			  	goto TERMINATE;
			}
			else if(lpstat == 108)
			{
				goto BRANCHING;
			}
				  	
		  	if(printing_in_setbranch) printf("getting x\n");
			status = CPXgetx (env_just_solve_mips, mip_ob1, x1, 0, cur_numcols-1);
			if(status) 
			{
				printf ("(%d) CPXgetx, Failed to get x values, error code %d\n", __LINE__,status);
		    		goto TERMINATE;
			}
			ob1_mip_objvals[0] = x1[obj1_index];
			ob1_mip_objvals[1] = x1[obj2_index];

			if(printing_in_setbranch) printf("plotting\n");
			if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",x1[obj1_index],x1[obj2_index]);
	
			if(!userhandle_up->x1) userhandle_up->x1 = calloc (cur_numcols,sizeof(double));
			if(!userhandle_down->x1) userhandle_down->x1 = calloc (cur_numcols,sizeof(double));
	
			for(i=0;i<cur_numcols;i++) 
			{
				userhandle_up->x1[i] = x1[i];
				userhandle_down->x1[i] = x1[i];
			}
	
			if(lpstat == 101 || lpstat == 102)
			{
				if(printing_in_setbranch) printf("the ob1 mip solution was optimal\n");
				ob1_mip_opt = 1;
		
				if(ob1_mip_objvals[0] < reduced_subprob_x_ub || ob1_mip_objvals[1] > reduced_subprob_y_lb)
				{
					reduced_subprob_x_ub = ob1_mip_objvals[0];
					if(there_will_only_be_points && integer_bb && integer_objective == 1) reduced_subprob_x_ub -= smallest_coef;
					reduced_subprob_y_lb = ob1_mip_objvals[1];
					if(there_will_only_be_points && integer_bb && integer_objective == 2) reduced_subprob_y_lb += smallest_coef;
					if(printing_in_setbranch) printf("after changing:\n");
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_ub,reduced_subprob_y_ub);
					bound_reduction = 1;
				}
		
				add_check = mock_insert(1,x1[obj1_index],x1[obj2_index],0,0,0,&tree);
				if(add_check)
				{
					if(printing_in_setbranch) printf("adding solution\n");
					for(i=0;i<cur_numcols;i++)
			      		{
				      		stored_x[x_rotation][i] = x1[i];
				      	}
				      	x_rotation = (x_rotation + 1) % num_x_to_store;
			 		add_check = 0;
			     		if(check_for_stopping_PSA_full)
					{
						check_for_stopping_PSA_full = 0;
						PSA_full(env,NULL,x1,NULL,NULL);
						check_for_stopping_PSA_full = 1;
					}
					else PSA_full(env,NULL,x1,NULL,NULL);
				}
		
				if(going_back) goto BRANCHING;
		
				RECHECK_IDEALS:
		
				add_check = 1;
				if(x1[obj1_index] >= x_ws[obj1_index] && x1[obj2_index] <= x_ws[obj2_index]) add_check = mock_insert(1,x1[obj1_index],
																x_ws[obj2_index],0,0,0,&tree);
				if(add_check)
				{
					if(printing_in_setbranch) printf("the right partial ideal is not dominated\n");
					projection = (slope)*(x1[obj1_index]-x_ws[obj1_index])+x_ws[obj2_index];
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],x1[obj1_index],x_ws[obj2_index],projection);
		/*			printf("slope used: %lf\n",slope);*/
					if(x1[obj1_index] >= x_ws[obj1_index] && x1[obj2_index] <= x_ws[obj2_index]) add_check = mock_insert(2,
											x1[obj1_index],projection,x_ws[obj1_index],x_ws[obj2_index],slope,&tree);
					if(add_check)
					{
						if(printing_in_setbranch) printf("the right partial ideal segment is not dominated\n");
						if(ws_mip_ran || ws_lp_int_feas || right_side_dom || left_side_dom) goto BRANCHING;
						else
						{
							if(printing_in_setbranch) printf("since here we haven't actually solved the ws mip, go back and solve ws mip\n");
							goto SOLVE_WS_MIP;
						}
					}
					else
					{
						if(printing_in_setbranch) printf("the right partial ideal segment is dominated. Move to consideration of ob2.\n");
						right_seg_dom = 1;
						reduced_subprob_x_ub = x_ws[obj1_index];
						reduced_subprob_y_lb = x_ws[obj2_index];
						if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
						bound_reduction = 1;
				
						if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
						{
							fathomed_by_dominated_local_ideal_segments++;
							fathoming = 1;
							if(printing_in_setbranch) printf("fathoming because left side already dominated\n");
/*							*useraction_p = CPX_CALLBACK_SET;*/
/*							nodecnt = 0;*/
						  	goto TERMINATE;
						}
				
						if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
						else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
						else goto BRANCHING;
					}
				}
				else
				{
					if(printing_in_setbranch) printf("the right partial ideal is dominated. Move to consideration of ob2.\n");
					right_pt_dom = 1;
					reduced_subprob_x_ub = x_ws[obj1_index];
					reduced_subprob_y_lb = x_ws[obj2_index];
					if(printing_in_setbranch) printf("after changing:\n");
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_ub,reduced_subprob_y_ub);
					bound_reduction = 1;
			
					if(left_side_dom || x_ws[obj1_index] - reduced_subprob_x_lb <= .0000001)
					{
						fathomed_by_dominated_local_ideal_pts++;
						fathoming = 1;
/*						*useraction_p = CPX_CALLBACK_SET;*/
/*						nodecnt = 0;*/
					  	goto TERMINATE;
					}
			
					if(!ob2_lp_been_solved) goto SOLVE_OB2_LP;
					else if(!ob2_mip_ran) goto SOLVE_OB2_MIP;
					else goto BRANCHING;
				}
			}
			else if(lpstat == 105 || lpstat == 107 || lpstat == 109 || lpstat == 111 || lpstat == 113 || lpstat == 116 || lpstat == 106)
			{
				if(printing_in_setbranch) printf("the mip solution was less than optimal\n");
				ob1_mip_opt = 0;
				add_check = mock_insert(1,x1[obj1_index],x1[obj2_index],0,0,0,&tree);
				if(add_check)
				{
					if(printing_in_setbranch) printf("adding solution\n");
					for(i=0;i<cur_numcols;i++)
					{
			      			stored_x[x_rotation][i] = x_ws[i];
			      		}
					x_rotation = (x_rotation + 1) % num_x_to_store;
			      		add_check = 0;
			      		if(check_for_stopping_PSA_full)
					{
						check_for_stopping_PSA_full = 0;
						PSA_full(env,NULL,x1,NULL,NULL);
						check_for_stopping_PSA_full = 1;
					}
					else PSA_full(env,NULL,x1,NULL,NULL);
				}
		
				if(going_back) goto BRANCHING;
		
				status = CPXgetbestobjval (env_just_solve_mips, mip_ob1, &best_bound);
				double obj_off;
			    	status = CPXgetobjoffset(env_just_solve_mips, mip_ob1, &obj_off );
				best_bound -= obj_off;
				reduced_subprob_x_ub = best_bound;
				
				goto BRANCHING;
				
				projection = (-1/1000.)*ob1_lp_objvals[1]+best_bound;
				if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",projection,ob1_lp_objvals[1]);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",-65.,best_bound-5./(-1000.),(-1000.)*(-best_bound-65),-5.);

				add_check = mock_insert(1,projection,ob1_lp_objvals[1],0,0,0,&tree);
				if(add_check)
				{
					if(printing_in_setbranch) printf("the projection onto the best objval level curve is not dominated\n");
					goto BRANCHING;
				}
				else
				{
					if(printing_in_setbranch) printf("the projection onto the best ob1 level curve is dominated\n");
					projection = (-1/1000.)*x_ws[obj2_index]+best_bound;
					if(projection <= x_ws[obj1_index])
					{
						printf("projection from ws soln dominated\n");
						exit(0);
					}
					else
					{
						if(printing_in_setbranch) printf("projection from ws soln not dominated\n");
						goto BRANCHING;
					}
				}
			}
			else
			{
				printf("the status of the FAILED mipsolve: %d\n",lpstat);
				exit(0);
			}
	
			/**** Here we solve and process the LP solution associated with objective f_2 ************************/
	
			SOLVE_OB2_LP:
	
		/*	if(left_side_dom || (!ws_mip_ran && ws_lp_sol_on_top_bd) ) goto BRANCHING;*/
	
			if(userhandle_current && userhandle_current->ob2_still_feas) goto SOLVE_OB2_MIP;
	
			lp_ob2 = CPXcloneprob (env, nodelp2, &status);
		  	if ( status ) {
		    		printf ("CPXcloneprob, Failed to clone problem, error code %d\n", status);
				goto TERMINATE;
		  	}
		  	chg_coefs(env, lp_ob2, indexes, -.000000001);
		  		
		  	ob2_lp_been_solved = 1;
		  	status = CPXlpopt (env, lp_ob2);
		 	if ( status ) {
		   		printf ("%s(%d): CPXlpopt, Failed to solve relaxation,  error code %d\n",
		   			__FILE__, __LINE__, status);
					goto TERMINATE;
		  	}
		  	
		  	status = CPXgetx (env, lp_ob2, x2, 0, cur_numcols-1);
			if(status) 
			{
				printf ("(%d) CPXgetx, Failed to get x values, error code %d\n", __LINE__,status);
		   			goto TERMINATE;
			}
			ob2_lp_objvals[0] = x2[obj1_index];
			ob2_lp_objvals[1] = x2[obj2_index];
	
			if(ob2_lp_objvals[1] < reduced_subprob_y_ub)
			{
		/*		reduced_subprob_x_lb = ob2_lp_objvals[0];*/
				reduced_subprob_y_ub = ob2_lp_objvals[1];
				if(printing_in_setbranch) printf("after changing:\n");
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
						reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_ub,reduced_subprob_y_ub);
				bound_reduction = 1;
			}

			if(printing_in_setbranch) printf("plotting ob2 lp soln\n");
			if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",x2[obj1_index],x2[obj2_index]);
	
			all_feas = 1;
			for(i=0;i<total_num_integer;i++)
			{
				k = integer_indices[i];
			  	diff = x2[k] - floor(x2[k]);
		/*		if(printing_in_setbranch) printf("diff%d: %lf\n",i,diff);*/
		  		if( diff >= .00001 && diff <= .99999)
		  		{
		  			if(printing_in_setbranch) printf("changing frac index to %d (%d)\n",integer_indices[i],__LINE__);
		/*  			printf("plot(%lf,%lf,'ko');\n",x2[obj1_index],x2[obj2_index]);*/
				  	frac_index = k;
				  	frac_val = x2[k];
				  	frac_values[k] = frac_val;
			 		if(frac_scores[i] > .0001)
			  		{
				  		frac_scores[i] += 1.;
		/*		  		num_frac++;*/
				  	}
				  	else
				  	{
				  		frac_scores[i] += multiplier*k;
			  			num_frac++;
				  	}
		/*		  	printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
				  	all_feas = 0;
		/*		  	printf("num_frac: %d\n",num_frac);*/
		/*			printf("(%d) frac_score_%d: %lf\n",__LINE__,i,frac_scores[i]);*/
				}
			}
			if(all_feas == 1)
			{
				if(printing_in_setbranch) printf("ob2 lp solution is integer feasible\n");
				all_feas = 0;
		
				ob2_mip_opt = 1;
				ob2_mip_objvals[0] = x2[obj1_index];
				ob2_mip_objvals[1] = x2[obj2_index];
				
				if(!userhandle_up->x2) userhandle_up->x2 = calloc ((cur_numcols),sizeof(double));
				if(!userhandle_down->x2) userhandle_down->x2 = calloc ((cur_numcols),sizeof(double));
	
				for(i=0;i<cur_numcols;i++) 
				{
					userhandle_up->x2[i] = x2[i];
					userhandle_down->x2[i] = x2[i];
				}
		
				ob2_lp_int_feas = 1;
			}

			add_check = mock_insert(1,x2[obj1_index],x2[obj2_index],0,0,0,&tree);
			if(add_check)
			{
				if(printing_in_setbranch) printf("the ob2 lp solution is not dominated\n");
				if(ob2_lp_int_feas == 1)
				{
					add_check = mock_insert(1,x2[obj1_index],x2[obj2_index],0,0,0,&tree);
					if(add_check)
					{
						if(printing_in_setbranch) printf("adding solution\n");
						for(i=0;i<cur_numcols;i++)
				      		{
					      		stored_x[x_rotation][i] = x2[i];
					      	}
					      	x_rotation = (x_rotation + 1) % num_x_to_store;
				 		add_check = 0;
				     		if(check_for_stopping_PSA_full)
						{
							check_for_stopping_PSA_full = 0;
							PSA_full(env,NULL,x2,NULL,NULL);
							check_for_stopping_PSA_full = 1;
						}
						else PSA_full(env,NULL,x2,NULL,NULL);
					}
			
					add_check = mock_insert(1,x_ws[obj1_index],x2[obj2_index],0,0,0,&tree);
					if(add_check)
					{
						if(printing_in_setbranch) printf("the left partial ideal is not dominated\n");
				
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],
										x2[obj2_index]);
						projection = (1./slope)*(x2[obj2_index]-x_ws[obj2_index])+x_ws[obj1_index];
						add_check = mock_insert(2,x_ws[obj1_index],x_ws[obj2_index],projection,x2[obj2_index],slope,&tree);
						if(add_check)
						{
							endpoint2_x = projection;
							endpoint2_y = x2[obj2_index];
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],
											x2[obj2_index]);
							if(printing_in_setbranch) printf("the left partial ideal segment is not dominated\n");
					
							status = CPXgetbase (env, lp_ob2, cstat, rstat);
							if(printing_in_setbranch) printf("calling PSA left\n");
							int reduce_val = PSA_left(env, lp_ob2, x2, cstat, rstat, seqnum, 0, lp_ob2);
			  				if(reduce_val == 2)
			  				{
			  					if(printing_in_setbranch) printf("fathoming node for completed PSA\n");
							  	fathomed_by_PSA_completion++;
								fathoming = 1;
/*								*useraction_p = CPX_CALLBACK_SET;*/
/*								nodecnt = 0;*/
							  	goto TERMINATE;
			  				}
			  				else if(reduce_val == 1)
			  				{
			  					reduced_subprob_x_lb = sub_pr2_x_lb;
				  				reduced_subprob_y_ub = sub_pr2_y_ub;
				  				if(printing_in_setbranch) printf("after changing:\n");
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
													reduced_subprob_y_lb,reduced_subprob_y_ub);
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
													reduced_subprob_y_ub,reduced_subprob_y_ub);
				  				bound_reduction = 1;
			  				}
							if(!ws_mip_ran && !ws_lp_int_feas && !right_side_dom) goto SOLVE_WS_MIP;
							else goto BRANCHING;
						}
						else
						{
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],
											x2[obj2_index]);
							if(printing_in_setbranch) printf("Left partial ideal segment is dominated. Fathom (%d)\n",__LINE__);
							if(right_pt_dom || right_side_dom) fathomed_by_dominated_local_ideal_pts++;
							else if(right_seg_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
							else
							{
								reduced_subprob_x_lb = x_ws[obj1_index];
				  				reduced_subprob_y_ub = x_ws[obj2_index];
				  				if(printing_in_setbranch) printf("after changing:\n");
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
													reduced_subprob_y_lb,reduced_subprob_y_ub);
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
													reduced_subprob_y_ub,reduced_subprob_y_ub);
				  				bound_reduction = 1;
				  				goto BRANCHING;
							}
							fathoming = 1;
/*							*useraction_p = CPX_CALLBACK_SET;*/
/*						  	nodecnt = 0;*/
							goto TERMINATE;
						}
					}
					else
					{
						left_side_dom = 1;
						if(printing_in_setbranch) printf("Left partial ideal is dominated. Fathom\n");
						if(right_pt_dom || right_side_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
						else if(right_seg_dom) fathomed_by_dominated_local_ideal_segments++;
						else
						{
							reduced_subprob_x_lb = x_ws[obj1_index];
			  				reduced_subprob_y_ub = x_ws[obj2_index];
			  				if(printing_in_setbranch) printf("after changing:\n");
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
												reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
												reduced_subprob_y_lb,reduced_subprob_y_lb);					
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
												reduced_subprob_y_ub,reduced_subprob_y_ub);
			  				bound_reduction = 1;
			  				goto BRANCHING;
						}
						fathoming = 1;
/*						*useraction_p = CPX_CALLBACK_SET;*/
/*					  	nodecnt = 0;*/
						goto TERMINATE;
					}
				}
				else
				{
					if(printing_in_setbranch) printf("ob2 lp solution is not integer feasible\n");
			
					if(!right_side_dom && ws_lp_dom && !ob1_lp_dom && remove_dominated_middle)
					{
		/*				printf("running this 2\n");*/
						if(!ws_sol_interior) 
						{
							status = CPXgetbase (env, nodelp_copy, cstat_ws, rstat_ws);
							if(status) 
							{
								printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
								goto SOLVE_OB2_MIP;
							}
						}
						else 
						{
							CPXgetbase (env, nodelp, cstat_ws, rstat_ws);
							if(status) 
							{
								printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
								goto SOLVE_OB2_MIP;
							}
						}
						if(printing_in_setbranch) printf("status: %d\n",status);
		
						sub_pr1_x_ub = x_ws[obj1_index];
						sub_pr1_y_lb = x_ws[obj2_index];
						sub_pr2_x_lb = x_ws[obj1_index];
						sub_pr2_y_ub = x_ws[obj2_index];
					  	PSA_right_check = PSA_reduce_right(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, indexes, seqnum);
					  	PSA_left_check = PSA_reduce_left(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, indexes);
					  	
						if(PSA_right_check == 2 && PSA_left_check == 2)
				  		{
					  		if(printing_in_setbranch) printf("fathoming node %d for PSA completion (%d)\n",seqnum,__LINE__);
						  	fathomed_by_dominated_lb++;
							fathoming = 1;
/*							*useraction_p = CPX_CALLBACK_SET;*/
/*							nodecnt = 0;*/
					  		goto TERMINATE;
				  		}
					  	else if(PSA_right_check == 2)
					  	{
							if(printing_in_setbranch) printf("the left subproblem is now empty by PSA completion\n");
					  		reduced_subprob_x_lb = x_ws[obj1_index];
					  		reduced_subprob_y_ub = x_ws[obj2_index];
					  		reduced_subprob_x_ub = ubs[0];
							reduced_subprob_y_lb = lbs[1];
							if(printing_in_setbranch) printf("after changing:\n");
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
												reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
												reduced_subprob_y_lb,reduced_subprob_y_lb);					
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
												reduced_subprob_y_ub,reduced_subprob_y_ub);
							bound_reduction = 1;
							left_side_dom = 1;
					  		goto SOLVE_OB1_LP;
					  	}
					  	else if(PSA_left_check == 2)
						{
				  			if(printing_in_setbranch) printf("the right subproblem is now empty by PSA completion\n");
					  		reduced_subprob_x_lb = lbs[0];
					  		reduced_subprob_y_ub = ubs[1];
							reduced_subprob_x_ub = x_ws[obj1_index];
				  			reduced_subprob_y_lb = x_ws[obj2_index];
				  			if(printing_in_setbranch) printf("after changing:\n");
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
												reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
												reduced_subprob_y_lb,reduced_subprob_y_lb);					
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
												reduced_subprob_y_ub,reduced_subprob_y_ub);
				  			bound_reduction = 1;
				  			right_side_dom = 1;
					  		goto SOLVE_OB2_LP;
					  	}
					  	else if(pareto_branching && sub_pr2_x_lb > sub_pr1_x_ub)
						{
							pareto = 1;
							goto TERMINATE;
/*					  		goto PARETO_BRANCH;*/
					  	}
					  	else if(!ob1_lp_been_solved)
					  	{
					  		if(printing_in_setbranch) printf("none of the above\n");
					  		goto SOLVE_OB1_LP;
					  	}
					  }
					  else if(!right_side_dom && ws_lp_dom && remove_dominated_middle)
					  {
		/*			  	printf("running this 3\n");*/
						if(!ws_sol_interior) 
						{
							status = CPXgetbase (env, nodelp_copy, cstat_ws, rstat_ws);
							if(status) 
							{
								printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
								goto SOLVE_OB2_MIP;
							}
						}
						else 
						{
							CPXgetbase (env, nodelp, cstat_ws, rstat_ws);
							if(status) 
							{
								printf("(%d) Failed to get basis. Error code: %d\n",__LINE__,status);
								goto SOLVE_OB2_MIP;
							}
						}
						if(printing_in_setbranch) printf("status: %d\n",status);
		
						sub_pr1_x_ub = x_ws[obj1_index];
						sub_pr1_y_lb = x_ws[obj2_index];
						sub_pr2_x_lb = x_ws[obj1_index];
						sub_pr2_y_ub = x_ws[obj2_index];
					  	PSA_right_check = PSA_reduce_right(env, nodelp_copy, x_ws, cstat_ws, rstat_ws, indexes, seqnum);
					  	
					  	if(PSA_right_check == 2)
					  	{
							if(printing_in_setbranch) printf("the left subproblem is now empty by PSA completion\n");
					  		reduced_subprob_x_lb = x_ws[obj1_index];
					  		reduced_subprob_y_ub = x_ws[obj2_index];
							if(printing_in_setbranch) printf("after changing:\n");
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
												reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
												reduced_subprob_y_lb,reduced_subprob_y_lb);					
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
												reduced_subprob_y_ub,reduced_subprob_y_ub);
							bound_reduction = 1;
							left_side_dom = 1;
					  		goto SOLVE_OB1_LP;
					  	}
					  	else if(pareto_branching && sub_pr2_x_lb > sub_pr1_x_ub)
						{
							pareto = 1;
							goto TERMINATE;
/*					  		goto PARETO_BRANCH;*/
					  	}
					}
			
					if((ws_mip_ran || ob1_mip_ran) && !ob2_mip_ran && !going_back_for_lp2) goto SOLVE_OB2_MIP;
					else goto BRANCHING;
				}
			}
			else
			{
				if(printing_in_setbranch) printf("the ob2 lp solution is dominated\n");
				ob2_lp_dom = 1;
				add_check = mock_insert(1,x_ws[obj1_index], x2[obj2_index],0,0,0,&tree);
				if(add_check)
				{
					if(printing_in_setbranch) printf("the left partial ideal is not dominated\n");
					projection = (1./slope)*(x2[obj2_index]-x_ws[obj2_index])+x_ws[obj1_index];
					add_check = mock_insert(2,x_ws[obj1_index],x_ws[obj2_index],projection,x2[obj2_index],slope,&tree);
					if(add_check)
					{
						endpoint2_x = projection;
						endpoint2_y = x2[obj2_index];
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],x2[obj2_index]);
						if(printing_in_setbranch) printf("the left partial ideal segment is not dominated\n");
				
						status = CPXgetbase (env, lp_ob2, cstat, rstat);
				
						if(ob2_lp_int_feas && ob1_lp_int_feas && ws_lp_int_feas)
						{
							if(printing_in_setbranch) printf("calling PSA left\n");
							int reduce_val = PSA_left(env, lp_ob2, x2, cstat, rstat, seqnum, 0, lp_ob2);
			  				if(reduce_val == 2)
			  				{
			  					if(printing_in_setbranch) printf("fathoming node for completed PSA\n");
							  	fathomed_by_PSA_completion++;
								fathoming = 1;
/*								*useraction_p = CPX_CALLBACK_SET;*/
/*								nodecnt = 0;*/
							  	goto TERMINATE;
			  				}
			  				else if(reduce_val == 1)
			  				{
			  					reduced_subprob_x_lb = sub_pr2_x_lb;
				  				reduced_subprob_y_ub = sub_pr2_y_ub;
				  				if(printing_in_setbranch) printf("after changing:\n");
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
													reduced_subprob_y_lb,reduced_subprob_y_ub);
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
													reduced_subprob_y_ub,reduced_subprob_y_ub);
				  				bound_reduction = 1;
			  				}
						}
						else //if(ws_lp_dom && (!ob1_lp_been_solved || ob1_lp_dom))
						{
							if(printing_in_setbranch) printf("calling PSA reduce left\n");
							int PSA_reduce_left_val = PSA_reduce_left(env, lp_ob2, x2, cstat, rstat, indexes);
			  				if(PSA_reduce_left_val == 2)
			  				{
			  					if(printing_in_setbranch) printf("fathoming node for dominated lower bound\n");
								fathomed_by_dominated_lb++;
								fathoming = 1;
/*								*useraction_p = CPX_CALLBACK_SET;*/
/*								nodecnt = 0;*/
							  	goto TERMINATE;
			  				}
			  				else if(PSA_reduce_left_val == 1)
			  				{
			  					reduced_subprob_x_lb = sub_pr2_x_lb;
				  				reduced_subprob_y_ub = sub_pr2_y_ub;
				  				if(printing_in_setbranch) printf("after changing:\n");
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
													reduced_subprob_y_lb,reduced_subprob_y_ub);
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
													reduced_subprob_y_ub,reduced_subprob_y_ub);
				  				bound_reduction = 1;
				  				bound_reduced_from_PSA_reduce_left = 1;
				  				goto BRANCHING;
			  				}
		  				}
			  				
						if((!ob2_lp_int_feas && (right_side_dom || right_seg_dom || right_pt_dom))
							 || (ws_mip_ran && ob1_mip_ran && !ob2_mip_ran && !ob2_lp_int_feas && !going_back_for_lp2)) goto SOLVE_OB2_MIP;
						else if(!right_side_dom && !right_seg_dom && !right_pt_dom && ws_mip_ran && !ob1_mip_ran && !ob1_lp_int_feas) goto SOLVE_OB1_MIP;
						else if(!ws_mip_ran && !ws_lp_int_feas) goto SOLVE_WS_MIP;
						else goto BRANCHING;
					}
					else
					{
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],x2[obj2_index]);
						if(printing_in_setbranch) printf("Left partial ideal segment is dominated. Fathom (%d)\n",__LINE__);
						if(right_pt_dom || right_side_dom) fathomed_by_dominated_local_ideal_pts++;
						else if(right_seg_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
						else
						{
							reduced_subprob_x_lb = x_ws[obj1_index];
			  				reduced_subprob_y_ub = x_ws[obj2_index];
			  				if(printing_in_setbranch) printf("after changing:\n");
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
												reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
												reduced_subprob_y_lb,reduced_subprob_y_lb);					
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
												reduced_subprob_y_ub,reduced_subprob_y_ub);
			  				bound_reduction = 1;
			  				goto BRANCHING;
						}
						fathoming = 1;
/*						*useraction_p = CPX_CALLBACK_SET;*/
/*					  	nodecnt = 0;*/
						goto TERMINATE;
					}
				}
				else
				{
					left_side_dom = 1;
					if(printing_in_setbranch) printf("the left partial ideal is dominated. Fathom\n");
					if(right_pt_dom || right_side_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
					else if(right_seg_dom) fathomed_by_dominated_local_ideal_segments++;
					else
					{
						reduced_subprob_x_lb = x_ws[obj1_index];
		  				reduced_subprob_y_ub = x_ws[obj2_index];
		  				if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
		  				bound_reduction = 1;
		  				goto BRANCHING;
					}
					fathoming = 1;
/*					*useraction_p = CPX_CALLBACK_SET;*/
/*				  	nodecnt = 0;*/
				  	goto TERMINATE;
				}
			}
	
			/************* Here we solve the single objective MIP associated with objective f_2 if we determine its needed *************/
	
			SOLVE_OB2_MIP:
			
			ob2_mip_ran = 1;

		/*	if(exact_mips) goto SKIP_THIS7;*/
			if(userhandle_current && userhandle_current->ob2_still_feas)
			{
				if(printing_in_setbranch) printf("obj2 soln from parent node is still feasible\n");
				for(j=0;j<cur_numcols;j++) x2[j] = userhandle_current->x2[j];
				ob2_mip_objvals[0] = x2[obj1_index];
				ob2_mip_objvals[1] = x2[obj2_index];
				
				userhandle_up->x2 = calloc (cur_numcols,sizeof(double));
				userhandle_down->x2 = calloc (cur_numcols,sizeof(double));
	
				for(i=0;i<cur_numcols;i++) 
				{
					userhandle_up->x2[i] = x2[i];
					userhandle_down->x2[i] = x2[i];
				}
				
				ob2_mip_opt = 1;
				mip_solved = 1;
		
				if(going_back) goto BRANCHING;
		
				add_check = mock_insert(1,x_ws[obj1_index],x2[obj2_index],0,0,0,&tree);
				if(add_check)
				{
					if(printing_in_setbranch) printf("the left partial ideal is not dominated\n");
			
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],x2[obj2_index]);
					projection = (1./slope)*(x2[obj2_index]-x_ws[obj2_index])+x_ws[obj1_index];
					add_check = mock_insert(2,x_ws[obj1_index],x_ws[obj2_index],projection,x2[obj2_index],slope,&tree);
					if(add_check)
					{
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],
										x2[obj2_index]);
						if(printing_in_setbranch) printf("the left partial ideal segment is not dominated\n");
						goto BRANCHING;
					}
					else
					{
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],
										x2[obj2_index]);
						if(printing_in_setbranch) printf("Left partial ideal segment is dominated. Fathom (%d)\n",__LINE__);
						if(right_pt_dom) fathomed_by_dominated_local_ideal_pts++;
						else if(right_seg_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
						else 
						{
							if(printing_in_setbranch) printf("somehow fathoming even though neither right pt nor seg has been labelled as dominated\n");
					
								reduced_subprob_x_lb = x_ws[obj1_index];
				  				reduced_subprob_y_ub = x_ws[obj2_index];
				  				if(printing_in_setbranch) printf("after changing:\n");
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
													reduced_subprob_y_lb,reduced_subprob_y_ub);
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
													reduced_subprob_y_ub,reduced_subprob_y_ub);
				  				bound_reduction = 1;
				  				goto BRANCHING;
					
						}
						fathoming = 1;
/*						*useraction_p = CPX_CALLBACK_SET;*/
/*					  	nodecnt = 0;*/
						goto TERMINATE;
					}
				}
				else
				{
					left_side_dom = 1;
					if(printing_in_setbranch) printf("the left partial ideal is dominated. Fathom\n");
					if(right_pt_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
					else if(right_seg_dom) fathomed_by_dominated_local_ideal_segments++;
					else 
					{
						if(printing_in_setbranch) printf("somehow fathoming even though neither right pt nor seg has been labelled as dominated\n");
				
								reduced_subprob_x_lb = x_ws[obj1_index];
				  				reduced_subprob_y_ub = x_ws[obj2_index];
				  				if(printing_in_setbranch) printf("after changing:\n");
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
													reduced_subprob_y_lb,reduced_subprob_y_ub);
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
								if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
													reduced_subprob_y_ub,reduced_subprob_y_ub);
				  				bound_reduction = 1;
				  				goto BRANCHING;
					
					}
					fathoming = 1;
/*					*useraction_p = CPX_CALLBACK_SET;*/
/*				  	nodecnt = 0;*/
				  	goto TERMINATE;
				}
			}
		
			if(!mip_ob2)
			{
				if(ob1_mip_ran && (!userhandle_current || !(userhandle_current->ob1_still_feas)))
				{
					if(bound_reduction)
					{
						int ind4[4] = {obj1_index,obj1_index,obj2_index,obj2_index};
						char lu4[4] = {'L','U','L','U'};
						double bds4[4] = {reduced_subprob_x_lb,reduced_subprob_x_ub,reduced_subprob_y_lb,reduced_subprob_y_ub};
						status = CPXchgbds (env_just_solve_mips, mip_ob1, 4, ind4, lu4, bds4);
						if ( status ) {
					   		printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
					  	}
					}
					mip_ob2 = CPXcloneprob (env_just_solve_mips, mip_ob1, &status);
				  	if ( status ) {
				    		printf ("(%d) CPXcloneprob, Failed to clone problem, error code %d\n", __LINE__, status);
						goto TERMINATE;
					}
				}
				else if(ws_mip_ran  && (!userhandle_current || !(userhandle_current->ws_still_feas)))
				{
					if(bound_reduction)
					{
						int ind4[4] = {obj1_index,obj1_index,obj2_index,obj2_index};
						char lu4[4] = {'L','U','L','U'};
						double bds4[4] = {reduced_subprob_x_lb,reduced_subprob_x_ub,reduced_subprob_y_lb,reduced_subprob_y_ub};
						status = CPXchgbds (env_just_solve_mips, nodelp_copy2, 4, ind4, lu4, bds4);
						if ( status ) {
					   		printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
					  	}
					}
					mip_ob2 = CPXcloneprob (env_just_solve_mips, nodelp_copy2, &status);
				  	if ( status ) {
				    		printf ("(%d) CPXcloneprob, Failed to clone problem, error code %d\n", __LINE__, status);
						goto TERMINATE;
					}
				}
				else
				{
					if(bound_reduction)
					{
						int ind4[4] = {obj1_index,obj1_index,obj2_index,obj2_index};
						char lu4[4] = {'L','U','L','U'};
						double bds4[4] = {reduced_subprob_x_lb,reduced_subprob_x_ub,reduced_subprob_y_lb,reduced_subprob_y_ub};
						status = CPXchgbds (env, nodelp2, 4, ind4, lu4, bds4);
						if ( status ) {
					   		printf ("%s(%d): CPXchgbds, Failed to change bounds, error code %d\n", __FILE__, __LINE__, status);
					  	}
					}
					mip_ob2 = CPXcloneprob (env_just_solve_mips, nodelp2, &status);
					if ( status ) {
			    			printf ("(%d) CPXcloneprob, Failed to clone problem, error code %d\n", __LINE__, status);
					    	goto TERMINATE;
			  		}
			  		CPXchgprobtype(env_just_solve_mips, mip_ob2, CPXPROB_MILP);
			  		status = CPXchgctype(env_just_solve_mips, mip_ob2, cur_numcols, indexes, xctype);
				}
			}
	
			if(bound_reduced_from_PSA_reduce_left)
			{
				int chg_indices[2] = {obj1_index,obj2_index};
				char chg_lu[2] = {'L','U'};
				double chg_bds[2] = {sub_pr2_x_lb,sub_pr2_y_ub};
				status = CPXchgbds (env_just_solve_mips, mip_ob2, 2, chg_indices, chg_lu, chg_bds);
				if ( status ){
				 	printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
				 	goto TERMINATE;}
			}
	
			chg_coefs(env_just_solve_mips, mip_ob2, indexes, -.000000001);
	
			if(!mip_solved) 
			{
				if(printing_in_setbranch) printf("here: %d\n",__LINE__);
				mip_solved = 1;
				num_nodes_with_mips_solved++;
			}

			if(printing_in_setbranch) printf("optimizing mip\n");
	
			RESOLVE_MIP:
	
			start_mipsolve = clock();
			nzcnt = 0;
			prev_numsols = 0;
	
			num_starts = CPXgetnummipstarts (env_just_solve_mips, global_mip);
			  	
		    	if(num_starts > global_num_starts)
		    	{
		    		global_num_starts = num_starts;
		    		global_startspace = cur_numcols*global_num_starts;
		    		int *temp = (int *) realloc (global_beg, (global_num_starts*2)*sizeof(int));
		    		if(temp != NULL) global_beg = temp;
		    		else printf("failed on assignment to global beg\n");
				int *temp2 = (int *) realloc (global_varindices, (global_startspace*2)*sizeof(int));
				if(temp2 != NULL) global_varindices = temp2;
		    		else printf("failed on assignment to global varind\n");
				double *temp3 = (double *) realloc (global_values, (global_startspace*2)*sizeof(double));
				if(temp3 != NULL) global_values = temp3;
		    		else printf("failed on assignment to global val\n");
				int *temp4 = (int *) realloc (global_effortlevel, (global_startspace*2)*sizeof(int));
				if(temp4 != NULL) global_effortlevel = temp4;
		    		else printf("failed on assignment to global eff\n");
		    	}
		    	
			int mip_starts_should_be_added = 1;
			if(num_starts)
			{
				status = CPXgetmipstarts (env_just_solve_mips, global_mip, &nzcnt, global_beg, global_varindices, 
						   global_values, global_effortlevel, global_startspace,
						   &surplus, 0, num_starts-1);
				if(seqnum == 14) for(i=0;i<nzcnt;i++) 
				{
					if(global_varindices[i] > cur_numcols || global_varindices[i] < 0)
					{
						printf("for some reason a mip start has become corrupted. Delete it!\n");
						int k = 0;
						while(i > global_beg[k]) k++;
						printf("value of k: %d\n",k);
						status = CPXdelmipstarts (env_just_solve_mips, global_mip, k-1,k-1);
						status = CPXgetmipstarts (env_just_solve_mips, global_mip, &nzcnt, global_beg, global_varindices, 
						   	global_values, global_effortlevel, global_startspace,
						   	&surplus, 0, num_starts-2);
						break;
					}
				}
						   
				if(mip_starts_should_be_added) status = CPXaddmipstarts (env_just_solve_mips, mip_ob2, num_starts, nzcnt, global_beg, global_varindices,
						   global_values, NULL, NULL);
			}
	
			CPXmipopt (env_just_solve_mips, mip_ob2);

			numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, mip_ob2);
			numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, mip_ob2);

			num_starts = numsolns - prev_numsols + numrep;

			prev_numsols = numsolns;

		    	if(num_starts > global_num_starts)
		    	{
		    		global_num_starts = num_starts;
		    		global_startspace = cur_numcols*global_num_starts;
		    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
				global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
				global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
				global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
		    	}

			status = CPXgetmipstarts (env_just_solve_mips, mip_ob2, &nzcnt, global_beg, global_varindices, 
					   global_values, global_effortlevel, global_startspace,
					   &surplus, 0, num_starts-1);
					   
			status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,
					   global_values, global_effortlevel, NULL);
		  	
		  	finish_mipsolve = clock();
			duration_mipsolve = 
				(double)(finish_mipsolve- start_mipsolve) / CLOCKS_PER_SEC;
			if(duration_mipsolve > max_time_to_solve_a_mip) max_time_to_solve_a_mip = duration_mipsolve;
			time_solving_mips += duration_mipsolve;
			if(printing_in_setbranch) printf("time to solve ob2 mip at seqnum %d: %lf\n",seqnum,duration_mipsolve);
			
			cumulative_time = (double)(finish_mipsolve - start_BB) / CLOCKS_PER_SEC;
			if(cumulative_time > max_time) goto BRANCHING;
	
			num_solns = CPXgetsolnpoolnumsolns (env_just_solve_mips, mip_ob2);
		  	
		  	insert_check = 0;
		  	
		  	if(num_solns >= prev_numsolns) times_to_run = num_solns - prev_numsolns;
		  	else times_to_run = num_solns;
		  	prev_numsolns = num_solns;
		  	
		  	for(j=0;j<times_to_run;j++)
		  	{
		  		status = CPXgetsolnpoolx (env_just_solve_mips, mip_ob2, j, x2, 0, cur_numcols-1);
			      	insert_check = mock_insert(1,x2[obj1_index],x2[obj2_index],0,0,0,&tree);
			      	if(insert_check)
			      	{
			      		if(branch_iterations < 5) for(i=0;i<cur_numcols;i++)
			      		{
			      			stored_x[x_rotation][i] = x2[i];
			      		}
			      		x_rotation = (x_rotation + 1) % num_x_to_store;
			      		insert_check = 0;
			      		if(check_for_stopping_PSA_full)
					{
						check_for_stopping_PSA_full = 0;
						PSA_full(env,NULL,x2,NULL,NULL);
						check_for_stopping_PSA_full = 1;
					}
					else PSA_full(env,NULL,x2,NULL,NULL);
			      	}
		      	}
		  	
		  	lpstat = CPXgetstat (env_just_solve_mips, mip_ob2);
			if(printing_in_setbranch) printf("solve status: %d\n",lpstat);
			
			if(keep_solving_infeasible_MIPs) while(lpstat == 108)
			{
				CPXmipopt (env_just_solve_mips, mip_ob2);
			 	lpstat = CPXgetstat (env_just_solve_mips, mip_ob2);
			}

			if(lpstat == 103 || lpstat == 119 || lpstat == 115)
			{
				if(printing_in_setbranch) printf("infeasible mip. Fathoming\n");
				fathoming = 1;
/*				*useraction_p = CPX_CALLBACK_SET;*/
/*				nodecnt = 0;*/
			  	goto TERMINATE;
			}
			else if(lpstat == 108)
			{
				goto BRANCHING;
				printf("Warning: Fathoming a node as infeasible because no feasible solution was found within time limit %lf s. This may cause incorrect solutions\n",time_limit);
				if(printing_in_setbranch) printf("infeasible mip. Fathoming\n");
				fathoming = 1;
/*				*useraction_p = CPX_CALLBACK_SET;*/
/*				nodecnt = 0;*/
			  	goto TERMINATE;
			}
			  	
		  	if(printing_in_setbranch) printf("getting x\n");
		  	status = CPXgetx (env_just_solve_mips, mip_ob2, x2, 0, cur_numcols-1);
			if(status) 
			{
				printf ("(%d) CPXgetx, Failed to get x values, error code %d\n", __LINE__,status);
				goto TERMINATE;
			}
			ob2_mip_objvals[0] = x2[obj1_index];
			ob2_mip_objvals[1] = x2[obj2_index];

			if(printing_in_setbranch) printf("plotting\n");
			if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",x2[obj1_index],x2[obj2_index]);
	
			if(!userhandle_up->x2) userhandle_up->x2 = calloc (cur_numcols,sizeof(double));
			if(!userhandle_down->x2) userhandle_down->x2 = calloc (cur_numcols,sizeof(double));
	
			for(i=0;i<cur_numcols;i++) 
			{
				userhandle_up->x2[i] = x2[i];
				userhandle_down->x2[i] = x2[i];
			}
	
		/*	for(i=0;i<total_num_integer;i++) printf("x%d: %lf\n",integer_indices[i],x2[integer_indices[i]]);*/

			if(lpstat == 101 || lpstat == 102)
			{
				if(printing_in_setbranch) printf("the ob2 mip solution was optimal\n");
				ob2_mip_opt = 1;
		
				if(ob2_mip_objvals[0] > reduced_subprob_x_lb || ob2_mip_objvals[1] < reduced_subprob_y_ub)
				{
					reduced_subprob_x_lb = ob2_mip_objvals[0];
					if(there_will_only_be_points && integer_bb && integer_objective == 1) reduced_subprob_x_lb += smallest_coef;
					reduced_subprob_y_ub = ob2_mip_objvals[1];
					if(there_will_only_be_points && integer_bb && integer_objective == 2) reduced_subprob_y_ub -= smallest_coef;
					if(printing_in_setbranch) printf("after changing:\n");
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
							reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_ub);
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_lb,reduced_subprob_y_lb);					
					if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
										reduced_subprob_y_ub,reduced_subprob_y_ub);
					bound_reduction = 1;
				}
		
				add_check = mock_insert(1,x2[obj1_index],x2[obj2_index],0,0,0,&tree);
				if(add_check)
				{
					if(printing_in_setbranch) printf("adding solution\n");
					for(i=0;i<cur_numcols;i++)
			      		{
				      		stored_x[x_rotation][i] = x2[i];
				      	}
				      	x_rotation = (x_rotation + 1) % num_x_to_store;
			 		add_check = 0;
			     		if(check_for_stopping_PSA_full)
					{
						check_for_stopping_PSA_full = 0;
						PSA_full(env,NULL,x2,NULL,NULL);
						check_for_stopping_PSA_full = 1;
					}
					else PSA_full(env,NULL,x2,NULL,NULL);
				}
		
				if(going_back) goto BRANCHING;
		
				if(ob2_mip_objvals[1] - reduced_subprob_y_lb <= .0000001 || (points_only && ob2_mip_objvals[1] - reduced_subprob_y_lb <= .01001)
					|| fabs(reduced_subprob_x_lb - reduced_subprob_x_ub) < .00001 || 
					fabs(reduced_subprob_y_lb - reduced_subprob_y_ub) <= .00001 )
				{
					fathoming = 1;
/*					*useraction_p = CPX_CALLBACK_SET;*/
/*			  		nodecnt = 0;*/
			  		goto TERMINATE;
				}
				else if(points_only)
				{
					if(ob2_mip_objvals[0] - reduced_subprob_x_lb < -.00001)
					{
						ob2_mip_ran = 0;
						ob2_mip_opt = 0;
						goto BEGINNING;
					}
					int chg_indices[2] = {obj1_index,obj2_index};
					char chg_lu[2] = {'L','U'};
					double chg_bds[2] = {x2[obj1_index] + .001, x2[obj2_index]};
					status = CPXchgbds (env_just_solve_mips, mip_ob2, 2, chg_indices, chg_lu, chg_bds);
		    			if ( status ){
		     			 	printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
		     			 	goto TERMINATE;}
		     			status = CPXchgbds (env, nodelp2, 2, chg_indices, chg_lu, chg_bds);
		    			if ( status ){
		     			 	printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
		     			 	goto TERMINATE;}
/*		     			points_only = 0;*/
		     			ob2_mip_opt = 0;
		     			ob2_mip_ran = 0;
					goto BEGINNING;
				}
		
				add_check = mock_insert(1,x_ws[obj1_index],x2[obj2_index],0,0,0,&tree);
				if(add_check)
				{
					if(printing_in_setbranch) printf("the left partial ideal is not dominated\n");
					projection = (1./slope)*(x2[obj2_index]-x_ws[obj2_index])+x_ws[obj1_index];
					add_check = mock_insert(2,x_ws[obj1_index],x_ws[obj2_index],projection,x2[obj2_index],slope,&tree);
					if(add_check)
					{
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],x2[obj2_index]);
						if(printing_in_setbranch) printf("the left partial ideal segment is not dominated\n");
						goto BRANCHING;
					}
					else
					{
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ws[obj1_index],projection,x_ws[obj2_index],x2[obj2_index]);
						if(printing_in_setbranch) printf("Left partial ideal segment is dominated. Fathom (%d)\n",__LINE__);
						if(right_pt_dom || right_side_dom) fathomed_by_dominated_local_ideal_pts++;
						else if(right_seg_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
						else
						{
							reduced_subprob_x_lb = x_ws[obj1_index];
			  				reduced_subprob_y_ub = x_ws[obj2_index];
			  				if(printing_in_setbranch) printf("after changing:\n");
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
									reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
												reduced_subprob_y_lb,reduced_subprob_y_ub);
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
									reduced_subprob_y_lb,reduced_subprob_y_lb);					
							if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
												reduced_subprob_y_ub,reduced_subprob_y_ub);
			  				bound_reduction = 1;
			  				goto BRANCHING;
						}
						fathoming = 1;
/*						*useraction_p = CPX_CALLBACK_SET;*/
/*				  		nodecnt = 0;*/
				  		goto TERMINATE;
					}
				}
				else
				{
					left_side_dom = 1;
					if(printing_in_setbranch) printf("the left partial ideal is dominated. Fathom\n");
					if(right_pt_dom || right_side_dom) fathomed_by_1_dominated_pt_1_dominated_segment++;
					else if(right_seg_dom) fathomed_by_dominated_local_ideal_segments++;
					else
					{
						reduced_subprob_x_lb = x_ws[obj1_index];
		  				reduced_subprob_y_ub = x_ws[obj2_index];
		  				if(printing_in_setbranch) printf("after changing:\n");
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_lb,
								reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_ub,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_ub);
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_lb,reduced_subprob_y_lb);					
						if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",reduced_subprob_x_lb,reduced_subprob_x_ub,
											reduced_subprob_y_ub,reduced_subprob_y_ub);
		  				bound_reduction = 1;
		  				goto BRANCHING;
					}
					fathoming = 1;
/*					*useraction_p = CPX_CALLBACK_SET;*/
/*				  	nodecnt = 0;*/
				  	goto TERMINATE;
				}

			}
			else if(lpstat == 105 || lpstat == 107 || lpstat == 109 || lpstat == 111 || lpstat == 113  || lpstat == 106)
			{
				if(printing_in_setbranch) printf("the mip solution was less than optimal\n");
				ob2_mip_opt = 0;
				add_check = mock_insert(1,x2[obj1_index],x2[obj2_index],0,0,0,&tree);
				if(add_check)
				{
					if(printing_in_setbranch) printf("adding solution\n");
					for(i=0;i<cur_numcols;i++)
			      		{
			      			stored_x[x_rotation][i] = x2[i];
			      		}
				      	x_rotation = (x_rotation + 1) % num_x_to_store;
				      	add_check = 0;
					if(check_for_stopping_PSA_full)
					{
						check_for_stopping_PSA_full = 0;
						PSA_full(env,NULL,x2,NULL,NULL);
						check_for_stopping_PSA_full = 1;
					}
					else PSA_full(env,NULL,x2,NULL,NULL);
				}
		
				if(going_back) goto BRANCHING;
		
				status = CPXgetbestobjval (env_just_solve_mips, mip_ob2, &best_bound);
				
				double obj_off;
			    	status = CPXgetobjoffset(env_just_solve_mips, mip_ob2, &obj_off );
				best_bound -= obj_off;
				reduced_subprob_y_ub = best_bound;
				
				goto BRANCHING;
				
				projection = (-.001)*(ob2_lp_objvals[0]-best_bound);
				if(printing_in_setbranch) printf("plot(%lf,%lf,'go');\n",ob2_lp_objvals[0],projection);
				if(printing_in_setbranch) printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",-65.,best_bound-5./(-.001), (-.001)*(-best_bound-65),-5.);

				add_check = mock_insert(1,ob1_lp_objvals[1],projection,0,0,0,&tree);
				if(add_check)
				{
					if(printing_in_setbranch) printf("projection onto the level curve is not dominated\n");
					goto BRANCHING;
				}
				else
				{
/*					if(printing_in_setbranch) */
					printf("projection onto the ob2 level curve is dominated\n");
					exit(0);
				}
			}
			else
			{
				printf("the status of the FAILED ob2 mipsolve: %d\n",lpstat);
				exit(0);
			}
	
			/************** Here we start the branching process ********************************************/	
	
			BRANCHING:
	
			finish_time = clock();
	
			time_processing_nodes += (double)(finish_time - start_time) / CLOCKS_PER_SEC;
	
			start_time = clock();
	
			if(printing_in_setbranch) printf("at branching\n");	
	
			if(exact_mips) goto SKIP_THIS9;
			if(depth < 3)
			{
				if(printing_in_setbranch) printf("depth: %d\n",depth);
				if(!ws_mip_ran && !ws_lp_int_feas)
				{
					if(printing_in_setbranch) printf("haven't solved ws mip, going back\n");
					going_back = 1;
					goto SOLVE_WS_MIP;
				}
				else if(!ob1_mip_ran && !ob1_lp_int_feas && !right_side_dom && !right_pt_dom && !right_seg_dom)
				{
					if(printing_in_setbranch) printf("haven't solved ob1 mip, going back\n");
					going_back = 1;
					goto SOLVE_OB1_MIP;
				}
				else if(!ob2_mip_ran && !ob2_lp_int_feas && !left_side_dom)
				{
					if(printing_in_setbranch) printf("haven't solved ob2 mip, going back\n");
					going_back = 1;
					goto SOLVE_OB2_MIP;
				}
			}
	
			SKIP_THIS9:
	
			if(exploit_objective_gaps && seqnum == 0 && fathoming != 1)
			{
				if(printing_in_setbranch) printf("depth: %d\n",depth);
				if(!ws_mip_ran && !ws_lp_int_feas)
				{
					if(printing_in_setbranch) printf("haven't solved ws mip, going back\n");
					going_back = 1;
					goto SOLVE_WS_MIP;
				}
				else if(!ob1_mip_ran && !ob1_lp_int_feas && !right_side_dom && !right_pt_dom && !right_seg_dom)
				{
					if(printing_in_setbranch) printf("haven't solved ob1 mip, going back\n");
					going_back = 1;
					goto SOLVE_OB1_MIP;
				}
				else if(!ob2_mip_ran && !ob2_lp_int_feas && !left_side_dom)
				{
					if(printing_in_setbranch) printf("haven't solved ob2 mip, going back\n");
					going_back = 1;
					goto SOLVE_OB2_MIP;
				}
			}
	
			if(printing_in_setbranch) printf("skip this 9\n");
	
			if(!ob1_lp_been_solved && !ob1_mip_ran && !right_side_dom)
			{
				if(exact_mips || ws_mip_opt ) 
				{
					if(printing_in_setbranch) printf("going back to solve ob1 lp\n");
					going_back_for_lp1 = 1;
					goto SOLVE_OB1_LP;
				}
			}
			if(!ob2_lp_been_solved && !ob2_mip_ran && !left_side_dom)
			{
				if(exact_mips || (ws_mip_opt && ob1_mip_opt && (right_pt_dom || right_side_dom || right_seg_dom))) 
				{
					if(printing_in_setbranch) printf("going back to solve ob2 lp\n");
					going_back_for_lp2 = 1;
					goto SOLVE_OB2_LP;
				}
			}
			
			if(ws_mip_ran && ob1_mip_ran && !ob2_mip_ran && !left_side_dom) goto SOLVE_OB2_MIP;
	
			COPY_MIPSTARTS:

			if(from_pareto) 
			{
				pareto = 1;
				goto TERMINATE;
/*				goto PARETO_BRANCH;*/
			}
	
			SKIP_THIS10:
			if(printing_in_setbranch) printf("skip this 10\n");
			if(ws_mip_opt && ob1_mip_opt && ob2_mip_opt)
			{
				int all_same = 1;
				for(i=0;i<total_num_integer;i++)
				{
					k = integer_indices[i];
					if( fabs(x_ws[k] - x1[k]) > .00001 || fabs(x_ws[k] - x2[k]) > .00001 )
					{
						all_same = 0;
						break;
					}
				}
				if(all_same)
				{
					if(printing_in_setbranch) printf("all 3 mips optimal and have same values of integer variables. Run PSA_full, then Fathom.\n");
					if(check_for_stopping_PSA_full)
					{
						check_for_stopping_PSA_full = 0;
						PSA_full(env,NULL,x1,NULL,NULL);
						check_for_stopping_PSA_full = 1;
					}
					else PSA_full(env,NULL,x1,NULL,NULL);
/*					*useraction_p = CPX_CALLBACK_SET;*/
					fathoming = 1;
/*			  		nodecnt = 0;*/
			  		goto TERMINATE;
				}
			}
			else if(ws_mip_opt && ob1_mip_opt && left_side_dom)
			{
				int all_same = 1;
				for(i=0;i<total_num_integer;i++)
				{
					k = integer_indices[i];
					if( fabs(x_ws[k] - x1[k]) > .00001)
					{
						all_same = 0;
						break;
					}
				}
				if(all_same)
				{
					if(printing_in_setbranch) printf("Left side was dominated and ws and ob1 mips optimal and have same values of integer variables. Run PSA_full, then Fathom.\n");
					if(check_for_stopping_PSA_full)
					{
						check_for_stopping_PSA_full = 0;
						PSA_full(env,NULL,x1,NULL,NULL);
						check_for_stopping_PSA_full = 1;
					}
					else PSA_full(env,NULL,x1,NULL,NULL);
/*					*useraction_p = CPX_CALLBACK_SET;*/
					fathoming = 1;
/*			  		nodecnt = 0;*/
			  		goto TERMINATE;
				}
			}
			else if(ws_mip_opt && ob2_mip_opt && right_side_dom)
			{
				int all_same = 1;
				for(i=0;i<total_num_integer;i++)
				{
					k = integer_indices[i];
					if(fabs(x_ws[k] - x2[k]) > .00001 )
					{
						all_same = 0;
						break;
					}
				}
				if(all_same)
				{
					if(printing_in_setbranch) printf("right side was dominated and ws and ob2 mips optimal and have same values of integer variables. Run PSA_full, then Fathom.\n");
					PSA_full(env,NULL,x_ws,NULL,NULL);
/*					*useraction_p = CPX_CALLBACK_SET;*/
					fathoming = 1;
/*			  		nodecnt = 0;*/
			  		goto TERMINATE;
				}
			}
			
			
			/********************* This is where the old mip solve started **********************/
			
/*			if(!indices)*/
/*			{*/
/*				indices   = (int *) malloc (cur_numcols * sizeof (int));*/
/*				for(i=0;i<cur_numcols;i++) indices[i] = i;*/
/*			}*/
/*			*/
/*			CPXLPptr nodelp_mip2 = CPXcloneprob (env_just_solve_mips, nodelp, &status);*/
/*			CPXchgprobtype(env_just_solve_mips, nodelp_mip2, CPXPROB_MILP);*/
/*			status = CPXchgctype(env_just_solve_mips, nodelp_mip2, cur_numcols, indices, xctype);*/
/*			*/
/*			status = CPXgetlb (env_just_solve_mips, nodelp_mip2, lbs, obj1_index, obj2_index);*/
/*			if ( status ) {*/
/*				printf ("(%d) Failed to get lb's for objectives. Error code %d\n",__LINE__, status);*/
/*				exit(0);*/
/*				goto TERMINATE;*/
/*			}*/
/*	*/
/*			status = CPXgetub (env_just_solve_mips, nodelp_mip2, ubs, obj1_index, obj2_index);*/
/*			if ( status ) {*/
/*				printf ("(%d) Failed to get ub's for objectives. Error code %d\n",__LINE__, status);*/
/*				goto TERMINATE;*/
/*			}*/
/*	*/
/*			slope = (ubs[1]-lbs[1])/(lbs[0]-ubs[0]);*/
/*	*/
/*			if(slope > -.001 || slope < -1000. || slope != slope) slope = -1.;*/
/*			*/
/*			double ran = (double) rand() / ( (double) RAND_MAX);*/
/*			if(ran < .5)*/
/*			{*/
/*				ran = (((double) rand() / ( (double) RAND_MAX)) - .5)/2.;*/
/*				slope += ran;*/
/*			}*/
/*			*/
/*			original_slope = slope;*/
/*			chg_coefs(env_just_solve_mips,nodelp_mip2,indices,slope);*/
/*			*/
/*			num_starts = CPXgetnummipstarts (env_just_solve_mips, global_mip);*/
/*			nzcnt = 0; //, prev_numsolns = 0;    */
/*			*/
/*			printf("number of mip starts to use: %d\n",num_starts);*/
/*				*/
/*		    	if(num_starts > global_num_starts)*/
/*		    	{*/
/*		    		printf("reallocation: (%d)\n",__LINE__);*/
/*		    		global_num_starts = num_starts;*/
/*		    		global_startspace = cur_numcols*global_num_starts;*/
/*		    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));*/
/*				global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));*/
/*				global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));*/
/*				global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));*/
/*		    	}*/

/*			status = CPXgetmipstarts (env_just_solve_mips, global_mip, &nzcnt, global_beg, global_varindices, */
/*					   global_values, global_effortlevel, global_startspace,*/
/*					   &surplus, 0, num_starts-1);*/
/*					   */
/*			status = CPXaddmipstarts (env_just_solve_mips, nodelp_mip2, num_starts, nzcnt, global_beg, global_varindices,*/
/*					   global_values, global_effortlevel, NULL);*/

/*			mip_start_time = clock();*/
/*			CPXmipopt (env_just_solve_mips, nodelp_mip2);*/
/*			mip_finish_time = clock();*/
/*			*/
/*			time_solving_mips += (double)(mip_finish_time - mip_start_time) / CLOCKS_PER_SEC;*/
/*			*/
/*			lpstat = CPXgetstat (env_just_solve_mips, nodelp_mip2);*/
/*			printf("lpstat: %d\n",lpstat);*/
/*			*/
/*			if(lpstat == 103)*/
/*			{*/
/*				*useraction_p = CPX_CALLBACK_DEFAULT;*/
/*				infeasible_in_cutcallback = 1;*/
/*				printf("(%d) freeing var ubs\n",__LINE__);*/
/*				if(var_ubs) free_and_null((char **) &var_ubs);*/
/*				printf("(%d) freeing var lbs\n",__LINE__);*/
/*				if(var_lbs) free_and_null((char **) &var_lbs);*/
/*				printf("(%d) freeing x\n",__LINE__);*/
/*				if(x) free_and_null((char **) &x);*/
/*				goto TERMINATE;*/
/*			}*/

/*			numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, nodelp_mip2);*/
/*			numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, nodelp_mip2);*/

/*			num_starts = numsolns - prev_numsols + numrep;*/
/*			printf("number of mip starts to use: %d\n",num_starts);*/
/*	*/
/*			prev_numsols = numsolns;*/
/*			if(num_starts > 0)*/
/*			{*/
/*				printf("number of mip starts to use: %d\n",num_starts);*/
/*				printf("numsols: %d, prev_numsols: %d, numrep: %d\n",numsolns,prev_numsols,numrep);*/
/*			    	if(num_starts > global_num_starts)*/
/*			    	{*/
/*			    		printf("reallocation: (%d)\n",__LINE__);*/
/*			    		global_num_starts = num_starts;*/
/*			    		global_startspace = cur_numcols*global_num_starts;*/
/*			    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));*/
/*					global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));*/
/*					global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));*/
/*					global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));*/
/*			    	}*/

/*				status = CPXgetmipstarts (env_just_solve_mips, nodelp_mip2, &nzcnt, global_beg, global_varindices, */
/*						   global_values, global_effortlevel, global_startspace,*/
/*						   &surplus, 0, num_starts-1);*/
/*						   */
/*				status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,*/
/*						   global_values, global_effortlevel, NULL);*/
/*				*/
/*				printf("about to check userhandle\n");*/
/*				if(!userhandle_current)*/
/*				{	   */
/*					userhandle_current = (user_data*) malloc( sizeof( user_data ) );*/
/*					userhandle_current->x_ws = (double *) malloc ( (cur_numcols+1)*sizeof(double));*/
/*					printf("setting current x1 to NULL\n");*/
/*					userhandle_current->x1 = NULL;*/
/*					userhandle_current->x2 = NULL;*/
/*					userhandle_current->ws_still_feas = 1;*/
/*					userhandle_current->ob1_still_feas = 0;*/
/*					userhandle_current->ob2_still_feas = 0;*/
/*				}*/
/*				else if(!(userhandle_current->x_ws))*/
/*				{*/
/*					userhandle_current->x_ws = (double *) malloc ( (cur_numcols+1)*sizeof(double));*/
/*				}*/
/*		*/
/*				status = CPXgetx (env_just_solve_mips, nodelp_mip2, userhandle_current->x_ws, 0, cur_numcols-1);*/
/*	  			if(status)*/
/*	  			{*/
/*					printf("Failed to get x-values from CPLEX. Status: %d Line: %d\n",status,__LINE__);*/
/*					exit(0);*/
/*	  			}*/
/*	  			userhandle_current->x_ws[cur_numcols] = slope;*/
/*  			}*/
/*			*/
/*			z = 0;*/
/*		    	double usercut_coefs[2] = {1.,-1./slope};*/
/*		    	int usercut_indices[2] = {obj1_index,obj2_index};*/
/*		    	char usercut_sense[1] = {'L'};*/
/*		    	double usercut_rhs = 0.;*/
/*		    	*/
/*		    	if(lpstat == 101 || lpstat == 102) status = CPXgetobjval (env_just_solve_mips, nodelp_mip2, &usercut_rhs);*/
/*		    	else status = CPXgetbestobjval (env_just_solve_mips, nodelp_mip2, &usercut_rhs);*/
/*		    	*/
/*		    	double obj_off;*/
/*		    	status = CPXgetobjoffset(env_just_solve_mips, nodelp_mip2, &obj_off );*/

/*			usercut_rhs -= obj_off;*/
/*		    	*/
/*		    	status = CPXcutcallbackaddlocal   (env,*/
/*							   cbdata,*/
/*				       			   wherefrom,*/
/*							   2,*/
/*							   usercut_rhs,*/
/*							   'L',*/
/*							   usercut_indices,*/
/*							   usercut_coefs);*/

/*			int num_solns = CPXgetsolnpoolnumsolns (env_just_solve_mips, nodelp_mip2);*/
/*			printf("number of soln pool solns: %d\n",num_solns);*/
/*		  	*/
/*		  	int insert_check = 0, j = 0;*/
/*		  	*/
/*		  	if(num_solns >= prev_numsolns) times_to_run = num_solns - prev_numsolns;*/
/*		  	else times_to_run = num_solns;*/
/*		  	prev_numsolns = num_solns;*/
/*		  	for(j=0;j<times_to_run;j++)*/
/*		  	{*/
/*		  		status = CPXgetsolnpoolx (env_just_solve_mips, nodelp_mip2, j, x, 0, cur_numcols-1);*/
/*		  		printf("plot(%lf,%lf,'go');\n",x[obj1_index],x[obj2_index]);*/
/*			      	insert_check = mock_insert(1,x[obj1_index],x[obj2_index],0,0,0,&tree);*/
/*			      	if(insert_check)*/
/*			      	{*/
/*			      		if(seqnum < 10) for(i=0;i<cur_numcols;i++)*/
/*			      		{*/
/*			      			stored_x[x_rotation][i] = x[i];*/
/*			      			printf("stored val: %lf\n",stored_x[x_rotation][i]);*/
/*			      		}*/
/*			      		x_rotation = (x_rotation + 1) % num_x_to_store;*/
/*			      		insert_check = 0;*/
/*			      		PSA_full(env_global,NULL,x,NULL,NULL);*/
/*			      	}*/
/*		      	}*/
				
/*			if(solve_extra_mips_during_cut_callback)*/
/*			{	*/
/*			*/
/*				double obvals[2] = {0.,0.};*/
/*			*/
/*				status = CPXgetx (env_just_solve_mips, nodelp_mip2, obvals, obj1_index, obj2_index);*/
/*	  			if(status)*/
/*	  			{*/
/*					printf("Failed to get x-values from CPLEX. Status: %d Line: %d\n",status,__LINE__);*/
/*					exit(0);*/
/*	  			}*/
/*						   */
/*				double slope2 = (obvals[1]-lbs[1])/(obvals[0]-ubs[0]);*/
/*				double slope3 = (ubs[1]-obvals[1])/(lbs[0]-obvals[0]);*/
/*			*/
/*				chg_coefs(env_just_solve_mips,nodelp_mip2,indices,slope2);*/
/*			*/
/*				CPXmipopt (env_just_solve_mips, nodelp_mip2);*/
/*			*/
/*				lpstat = CPXgetstat (env_just_solve_mips, nodelp_mip2);*/
/*			*/
/*				if(lpstat == 103)*/
/*				{*/
/*					*useraction_p = CPX_CALLBACK_DEFAULT;*/
/*					printf("(%d) freeing var ubs\n",__LINE__);*/
/*					if(var_ubs) free_and_null((char **) &var_ubs);*/
/*					printf("(%d) freeing var lbs\n",__LINE__);*/
/*					if(var_lbs) free_and_null((char **) &var_lbs);*/
/*					printf("(%d) freeing x\n",__LINE__);*/
/*					if(x) free_and_null((char **) &x);*/
/*					goto TERMINATE;*/
/*				}*/

/*				numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, nodelp_mip2);*/
/*				numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, nodelp_mip2);*/

/*				num_starts = numsolns - prev_numsols + numrep;*/
/*				printf("number of mip starts to use: %d\n",num_starts);*/
/*	*/
/*				prev_numsols = numsolns;*/
/*	*/
/*			    	if(num_starts > global_num_starts)*/
/*			    	{*/
/*			    		printf("reallocation: (%d)\n",__LINE__);*/
/*			    		global_num_starts = num_starts;*/
/*			    		global_startspace = cur_numcols*global_num_starts;*/
/*			    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));*/
/*					global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));*/
/*					global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));*/
/*					global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));*/
/*			    	}*/

/*				status = CPXgetmipstarts (env_just_solve_mips, nodelp_mip2, &nzcnt, global_beg, global_varindices, */
/*						   global_values, global_effortlevel, global_startspace,*/
/*						   &surplus, 0, num_starts-1);*/
/*						   */
/*				status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,*/
/*						   global_values, global_effortlevel, NULL);*/
/*						   */
/*				usercut_coefs[1] = -1./slope2;*/
/*			    	*/
/*			    	if(lpstat == 101 || lpstat == 102) status = CPXgetobjval (env_just_solve_mips, nodelp_mip2, &usercut_rhs);*/
/*			    	else status = CPXgetbestobjval (env_just_solve_mips, nodelp_mip2, &usercut_rhs);*/
/*			    	*/
/*			    	printf("adding user cut: x%d + %lfx%d <= %lf\n",obj1_index,-1./slope,obj2_index,usercut_rhs);*/
/*			    	*/
/*			    	status = CPXcutcallbackaddlocal   (env,*/
/*								   cbdata,*/
/*					       			   wherefrom,*/
/*								   2,*/
/*								   usercut_rhs,*/
/*								   'L',*/
/*								   usercut_indices,*/
/*								   usercut_coefs);*/
/*								 */
/*				chg_coefs(env_just_solve_mips,nodelp_mip2,indices,slope3);*/
/*			*/
/*				CPXmipopt (env_just_solve_mips, nodelp_mip2);*/
/*			*/
/*				lpstat = CPXgetstat (env_just_solve_mips, nodelp_mip2);*/
/*			*/
/*				if(lpstat == 103)*/
/*				{*/
/*					*useraction_p = CPX_CALLBACK_DEFAULT;*/
/*					goto TERMINATE;*/
/*				}*/

/*				numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, nodelp_mip2);*/
/*				numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, nodelp_mip2);*/

/*				num_starts = numsolns - prev_numsols + numrep;*/
/*				printf("number of mip starts to use: %d\n",num_starts);*/
/*	*/
/*				prev_numsols = numsolns;*/
/*	*/
/*			    	if(num_starts > global_num_starts)*/
/*			    	{*/
/*			    		printf("reallocation: (%d)\n",__LINE__);*/
/*			    		global_num_starts = num_starts;*/
/*			    		global_startspace = cur_numcols*global_num_starts;*/
/*			    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));*/
/*					global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));*/
/*					global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));*/
/*					global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));*/
/*			    	}*/

/*				status = CPXgetmipstarts (env_just_solve_mips, nodelp_mip2, &nzcnt, global_beg, global_varindices, */
/*						   global_values, global_effortlevel, global_startspace,*/
/*						   &surplus, 0, num_starts-1);*/
/*						   */
/*				status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,*/
/*						   global_values, global_effortlevel, NULL);*/
/*						   */
/*				usercut_coefs[1] = -1./slope3;*/
/*			    	*/
/*			    	if(lpstat == 101 || lpstat == 102) status = CPXgetobjval (env_just_solve_mips, nodelp_mip2, &usercut_rhs);*/
/*			    	else status = CPXgetbestobjval (env_just_solve_mips, nodelp_mip2, &usercut_rhs);*/
/*			    	*/
/*			    	printf("adding user cut: x%d + %lfx%d <= %lf\n",obj1_index,-1./slope,obj2_index,usercut_rhs);*/
/*			    	*/
/*			    	status = CPXcutcallbackaddlocal   (env,*/
/*								   cbdata,*/
/*					       			   wherefrom,*/
/*								   2,*/
/*								   usercut_rhs,*/
/*								   'L',*/
/*								   usercut_indices,*/
/*								   usercut_coefs);*/
/*			*/
/*			}*/

/*			CPXfreeprob(env_just_solve_mips, &nodelp_mip2);*/

			*useraction_p = CPX_CALLBACK_DEFAULT;
/*			printf("(%d) freeing var_ubs\n",__LINE__);*/
			if(var_ubs) free_and_null((char **) &var_ubs);
/*			printf("(%d) freeing var lbs\n",__LINE__);*/
			if(var_lbs) free_and_null((char **) &var_lbs);
/*			printf("(%d) freeing x\n",__LINE__);*/
			if(x) free_and_null((char **) &x);
			goto TERMINATE;
        	}
        	else
		{
			*useraction_p = CPX_CALLBACK_DEFAULT;
			goto TERMINATE;
		}
/*		printf("(%d) freeing var_ubs\n",__LINE__);*/
		if(var_ubs) free_and_null((char **) &var_ubs);
/*		printf("(%d) freeing var lbs\n",__LINE__);*/
		if(var_lbs) free_and_null((char **) &var_lbs);
/*		printf("(%d) freeing x\n",__LINE__);*/
		if(x) free_and_null((char **) &x);
	}
	else
	{
		*useraction_p = CPX_CALLBACK_DEFAULT;
		goto TERMINATE;
	}
	
	TERMINATE:
	
/*	if(x_ws && (x_ws[obj1_index] != x_ws[obj1_index] || x_ws[obj2_index] != x_ws[obj2_index]))*/
/* 	{*/
/* 		printf("x_ws solns got f'ed up. Exitting\n");*/
/* 		exit(0);*/
/* 	}*/
	
	finish_time = clock();
		 	
 	time_processing_nodes += (double)(finish_time - start_time) / CLOCKS_PER_SEC;
/* 	else time_branching += (double)(finish_time - start_time) / CLOCKS_PER_SEC;*/
	
/*	printf("arrived at terminate of cut callback\n");*/
	
	if(indices)    free_and_null((char **) &indices);
	if(rmatbeg)    free_and_null((char **) &rmatbeg);
	if(rmatind)    free_and_null((char **) &rmatind);
	if(rmatval)    free_and_null((char **) &rmatval);
	if(rhs_s)      free_and_null((char **) &rhs_s);
	if(local_lb)   free_and_null((char **) &local_lb);
	if(local_ub)   free_and_null((char **) &local_ub);
	
	if(nodelp2)
	{
		double l[2] = {0.,0.};
		double u[2] = {0.,0.};
	
		status = CPXgetlb (env, nodelp2, l, obj1_index, obj2_index);
		if ( status ) {
			printf ("Failed to get lb's for objectives. Error code %d\n", status);
			exit(0);
		}

		status = CPXgetub (env, nodelp2, u, obj1_index, obj2_index);
		if ( status ) {
			printf ("Failed to get ub's for objectives. Error code %d\n", status);
			exit(0);
		}
	
		reduced_subprob_x_lb = fmax(l[0],reduced_subprob_x_lb);
		reduced_subprob_y_ub = fmin(u[1],reduced_subprob_y_ub);
		reduced_subprob_x_ub = fmin(u[0],reduced_subprob_x_ub);
		reduced_subprob_y_lb = fmax(l[1],reduced_subprob_y_lb);
	}
	
	CPXfreeprob(env4, &nodelp_mip);
	status=CPXcloseCPLEX(&env4); 	
 	
 	if(fathoming == 1)
 	{
 		if(userhandle_current)
 		{
/*	 		printf("checking whether or not x solns are stored\n");*/
	 		if(userhandle_current->x_ws) free_and_null((char **) &userhandle_current->x_ws);
	/* 		printf("freeing current x1\n");*/
	 		if(userhandle_current->x1) free_and_null((char **) &userhandle_current->x1);
	/* 		printf("x2?\n");*/
	 		if(userhandle_current->x2) free_and_null((char **) &userhandle_current->x2);
	 		if(userhandle_current->prob) CPXfreeprob(env, &(userhandle_current->prob));
	 		free(userhandle_current);
	 		userhandle_current = NULL;
 		}
/* 		printf("fathoming, must free local userhandles\n");*/
 		if(userhandle_up)
	 	{
	 		if(userhandle_up->x_ws) free_and_null((char **) &userhandle_up->x_ws);
/*	 		printf("freeing up x1\n");*/
	 		if(userhandle_up->x1) free_and_null((char **) &userhandle_up->x1);
	 		if(userhandle_up->x2) free_and_null((char **) &userhandle_up->x2);
	 		if(userhandle_up->prob) CPXfreeprob(env, &(userhandle_up->prob));
	 		free_and_null((char **) &userhandle_up);
	 		userhandle_up = NULL;
	 	}
	 	if(userhandle_down)
	 	{
	 		if(userhandle_down->x_ws) free_and_null((char **) &userhandle_down->x_ws);
/*	 		printf("freeing down x1\n");*/
	 		if(userhandle_down->x1) free_and_null((char **) &userhandle_down->x1);
	 		if(userhandle_down->x2) free_and_null((char **) &userhandle_down->x2);
	 		if(userhandle_down->prob) CPXfreeprob(env, &(userhandle_down->prob));
	 		free_and_null((char **) &userhandle_down);
	 		userhandle_down = NULL;
	 	}
 	}
 	
/* 	SKIP_THIS11:*/
 	
/* 	printf("freeing\n");*/
/* 	if(x_ws != NULL) free_and_null((char **) &x_ws);*/
/* 	printf("freeing\n");*/
/* 	if(x1 != NULL) free_and_null((char **) &x1);*/
/* 	printf("freeing\n");*/
/* 	if(x2 != NULL) free_and_null((char **) &x2);*/
/* 	printf("freeing\n");*/
 	if(feas != NULL) free_and_null((char **) &feas);
/* 	printf("freeing\n");*/
 	if(low_up != NULL) free_and_null((char **) &low_up);
/* 	printf("freeing\n");*/
	if(lb_ != NULL) free_and_null((char **) &lb_);
/* 	printf("freeing\n");*/
	if(ub_ != NULL) free_and_null((char **) &ub_);
/* 	printf("freeing\n");*/
	if(bds_br1 != NULL) free_and_null((char **) &bds_br1);
/* 	printf("freeing\n");*/
	if(bds_br2 != NULL) free_and_null((char **) &bds_br2);
/* 	printf("freeing\n");*/
	if(br_ind != NULL) free_and_null((char **) &br_ind);
/* 	printf("freeing\n");*/
	if(cstat != NULL) free_and_null((char **) &cstat);
/* 	printf("freeing\n");*/
	if(rstat != NULL) free_and_null((char **) &rstat);
/* 	printf("freeing\n");*/
	if(cstat_ws != NULL) free_and_null((char **) &cstat_ws);
/* 	printf("freeing\n");*/
	if(rstat_ws != NULL) free_and_null((char **) &rstat_ws);
/* 	printf("freeing\n");*/
	if(var_ubs != NULL) free_and_null((char **) &var_ubs);
/* 	printf("freeing\n");*/
	if(var_lbs != NULL) free_and_null((char **) &var_lbs);
/* 	printf("freeing\n");*/
	if(indexes != NULL) free_and_null((char **) &indexes);
/*	printf("freeing\n");*/
/*	if(frac_values != NULL) free_and_null((char **) &frac_values);*/



/* 	printf("freeing\n");*/
	CPXfreeprob(env_just_solve_mips, &nodelp_copy);
	CPXfreeprob(env_just_solve_mips, &nodelp_copy2);
/* 	printf("freeing\n");*/
	CPXfreeprob(env, &lp_ob1);
/* 	printf("freeing\n");*/
	CPXfreeprob(env, &lp_ob2);
/* 	printf("freeing\n");*/
	CPXfreeprob(env_just_solve_mips, &mip_ob1);
/* 	printf("freeing\n");*/
	CPXfreeprob(env_just_solve_mips, &mip_ob2);
/* 	printf("freeing last\n");*/
	CPXfreeprob(env, &nodelp2);
	CPXfreeprob(env, &cut_prob);

	return(status);
}

double percentage = 0.;
int declared = 0;
double f_quart = -100000000000000., t_quart = -10000000000000.;

int CPXPUBLIC 
userselectnode (CPXCENVptr env,
                void       *cbdata,
                int        wherefrom,
                void       *cbhandle,
                int        *nodenum_p,
                int        *useraction_p)
{
  	int status = 0, i = 0, nodesleft = 0, bestnode = 0, second_bestnode = 0, third_bestnode;

   	*useraction_p = CPX_CALLBACK_DEFAULT;
   	if(!show_progress) goto TERMINATE;

   	status = CPXgetcallbackinfo (env, cbdata, wherefrom,
                                CPX_CALLBACK_INFO_NODES_LEFT,
                                &nodesleft);
   	if ( status )
   	{
   		printf("Trouble getting nodes left from node callback. Error: %d\n",status);
   		goto TERMINATE;
   	}

/*	printf("branch iterations: %d\n",branch_iterations);*/
	if(branch_iterations >= last_time_shown + show_frequency)// || branch_iterations <= 1)
	{
		if(first_time_showing_progress)
		{
			printf("\nNodes processed\t Nodes left\t Duality Gap\n");
			first_time_showing_progress = 0;
		}
		int nodesleft = 0;
		if(branch_iterations!= 1) last_time_shown = branch_iterations;
		
		node *tree_copy = NULL, *tree_copy2 = NULL;
		tree_copy = copy_tree(tree,empty_node);
		tree2 = copy_tree(tree,empty_node);
		
/*		if(branch_iterations == 8 || branch_iterations == 9) */
/*		print_inorder(tree_copy,1);*/
/*		print_inorder(tree2,2);*/
		
		user_data *temp_uh = NULL;
		
		if ((status = CPXgetcallbackinfo (env, cbdata, wherefrom, 
			CPX_CALLBACK_INFO_NODES_LEFT, &nodesleft))) goto TERMINATE;
			
/*		printf("Number of open nodes: %d\n",nodesleft);*/
			
/*		printf("building dual bound\n");*/
/*		printf("**************************\n");*/
		for (i=0;i<nodesleft;i++) 
		{
/*			temp_uh = NULL;*/
	      		status = CPXgetcallbacknodeinfo (env, cbdata, wherefrom,
		                               i,
		                               CPX_CALLBACK_INFO_NODE_USERHANDLE,
		                               &temp_uh);
/*		        if(status) printf("didn't get uh, error: %d\n",status);*/
/*		        if(temp_uh) */
/*		        {*/
/*		        	printf("i: %d\t user_handle existed\n",i);*/
/*		        	if(temp_uh->prob) printf("\t problem also existed\n");*/
/*		        }*/
/*		        else printf("i: %d\t user_handle didnt exist\n",i);*/
			if(temp_uh) 
			{
				build_dual_bd(env,temp_uh->prob);
/*				if(control_node_selection) printf("plot(%lf,%lf,'go');\n",temp_uh->f1,temp_uh->f2);*/
			}
/*			if(branch_iterations == 3 || branch_iterations == 4){*/
/*			printf("**************************\n");*/
/*			print_inorder(tree2,2);*/
/*			printf("**************************\n");}*/
		}
/*		printf("**************************\n");*/
/*		printf("calculating max prox dist\n");*/

		//if(branch_iterations == 8 || branch_iterations == 9){
/*		printf("**************************\n");*/
/*		print_inorder(tree2,2);*/
/*		printf("**************************\n");*/
/*		print_preorder(tree2,NULL);*/
/*		printf("**************************\n");//}*/

		if(hausdorff_or_hypervolume == 2)
		{
			double hy_v2 = get_hypervolume(tree2, 1, 0.);
/*			printf("dual hypervolume: %lf\n",hy_v2);*/
			double hy_v = get_hypervolume(tree_copy, 1, 0.);
/*			printf("primal hypervolume: %lf\n",hy_v);*/
			double pd = 100.*(hy_v2-hy_v)/hy_v2;
/*			printf("percent difference: %lf\n",pd);*/
	/*		exit(0);*/
	
			printf("%d\t\t %d\t\t %lf\n",branch_iterations,nodesleft,pd);//,length_percent_of_max_range);
			
			if(pd < duality_gap_limit) 
			{
				break_early = 1;
				clock_t now_time = clock();
	
				cumulative_time = (double)(now_time - start_BB) / CLOCKS_PER_SEC;
				printf("Time spent until < %lf%% hypervolume gap: %lf seconds\n",duality_gap_limit,cumulative_time);
				exit(0);
			}
		}
		else
		{
			double hd_dist = calculate_max_proximal_hd_dist(tree_copy, tree2, 0.);
			double new_val = get_nadirs(tree_copy, 1, 0.);
			hd_dist = fmax(hd_dist,new_val);
			double percent_of_max_range = 100.*hd_dist/max_range;	
		
	/*		double length = get_length(tree2,0.);	*/
	/*		double length_percent_of_max_range = 100.*length/max_range;*/
		
	/*		if(percent_of_max_range < duality_gap_limit) break_early = 1;*/

			printf("%d\t\t %d\t\t %lf\n",branch_iterations,nodesleft,percent_of_max_range);//,length_percent_of_max_range);

			if(percent_of_max_range < duality_gap_limit) 
			{
				break_early = 1;
				clock_t now_time = clock();
	
				cumulative_time = (double)(now_time - start_BB) / CLOCKS_PER_SEC;
				printf("Time spent until < %lf%% hausdorff gap: %lf seconds\n",duality_gap_limit,cumulative_time);
				exit(0);
			}
		}
		
/*		printf("***********************************\n");*/
/*		print_inorder(tree,1);*/
/*		printf("***********************************\n");*/
/*		print_inorder(tree_copy,1);*/
/*		printf("***********************************\n");*/
/*		print_inorder(tree2,2);*/
/*		printf("***********************************\n");*/
/*		exit(0);*/
		destroy_tree(tree_copy);
		destroy_tree(tree2);
/*		printf("destroying dual bound\n");*/
		tree2 = NULL;
	}
	
/*	double percentage = 0.;*/
	if(allow_changing_control_percentage) percentage = fmax(percentage, starting_percentage);
	else percentage = control_percentage;
	
	if(control_node_selection) 
	{	
		double ran = (double) rand() / ( (double) RAND_MAX);
		if(ran < percentage)
		{
			if(!nodesleft) if ((status = CPXgetcallbackinfo (env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODES_LEFT, &nodesleft))) goto TERMINATE;
			user_data *temp_uh = NULL;
			double best_val = -10000000000000.;
			ran = (double) rand() / ( (double) RAND_MAX);
			if(ran < .5)
			{
				for (i=0;i<nodesleft;i++) 
				{
			      		status = CPXgetcallbacknodeinfo (env, cbdata, wherefrom,
						               i,
						               CPX_CALLBACK_INFO_NODE_USERHANDLE,
						               &temp_uh);
					if(temp_uh->f1 > best_val)
					{
						best_val = temp_uh->f1;
						third_bestnode = second_bestnode;
						second_bestnode = bestnode;
						bestnode = i;
					}
				}
	   		}
	   		else
	   		{
	   			for (i=0;i<nodesleft;i++) 
				{
			      		status = CPXgetcallbacknodeinfo (env, cbdata, wherefrom,
						               i,
						               CPX_CALLBACK_INFO_NODE_USERHANDLE,
						               &temp_uh);
					if(temp_uh->f2 > best_val)
					{
						best_val = temp_uh->f2;
						third_bestnode = second_bestnode;
						second_bestnode = bestnode;
						bestnode = i;
					}
				}
	   		}
	   		
/*	   		ran = (double) rand() / ( (double) RAND_MAX);*/
/*	   		if(ran < .15) bestnode = second_bestnode;*/
/*	   		else if(ran < .4) bestnode = third_bestnode;*/
	   		
	   		*nodenum_p = bestnode;  
	   		*useraction_p = CPX_CALLBACK_SET;
   		}
   		if(allow_changing_control_percentage && percentage < stopping_percentage) percentage += .005;
/*	   	printf("percentage: %lf\n",percentage);*/
   	}
	
	TERMINATE:

   	return (status);

} /* END userselectnode */

