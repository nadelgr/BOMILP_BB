#ifndef BB_BICRITERIA_H

/* Input problems should be available at testP1.lp and testP2.lp files */
/* Only difference between two testP.lp files are the objective function coefficients */
/* Note that these problems should be maximization type and equality constraints */
/* 
   /* Max c1x
   /* subject to Ax = b */
/*
  /* All feasible points found are written to output.txt */
/* Before running this code, FindNadir2.c should be run. Nadir points are then written to init_nadir.txt */ 

/* This version does not update Pareto set or nadir set. Initial points are used */


/*The diffence of this version from version ComputeFeasSol4.c is that nadir points are computed by another code
  This program only reads them */

/* The difference of this version from BB-bicriteriaMILP.c  is that 
   y_1^SE = max{c1x: x \in \tilde X}
   y_2^SE = min{c2x: x \in \tilde X}
   y_1^NW = min{c1x: x \in \tilde X}
   y_2^NW = max{c2x: x \in \tilde X}    */ 

/* Preprocessing is not avaialble */
/* Fathoming Rule 2b is closed. It does not work. */
/* The difference of this code from BB-bicriteriaMILP3.c is that I am leaving the CPXptr dlp free now. */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <time.h>
#include "cplex.h"

#ifdef SOL_tree
#include "max_tree.h"
#else
#ifdef SOL_list
#include "max_list.h"
#else
#error "Please define either SOL_tree or SOL_list"
#endif
#endif

#define REALLOC_STEP 1024
#define BILLION  1000000000L;
extern struct timespec start1, stop1;
extern double accum, struct_time2;
extern struct timespec start2, stop2;
extern double insert_time2;

extern node *tree;
extern node *tree2;
extern node *potential_branch_tree;

extern int printing;
extern int bd_reduction_during_branching;

extern int cnt;
extern double x_ideal;
extern double y_ideal;
extern int obj1_index;
extern int obj2_index;
extern int cut_row_index;
extern int insert_to_potential_branch_tree;
extern int fathomed_by_dominated_search_region;
extern int fathomed_by_dominated_local_ideal_pts;
extern int fathomed_by_dominated_local_ideal_segments;
extern int fathomed_by_1_dominated_pt_1_dominated_segment;
extern int num_nodes_with_mips_solved;
extern int fathomed_by_PSA_completion;
extern int fathomed_by_dominated_lb;
extern int num_x_to_store;
extern int heur_limit;
extern int pure_binary;
extern int num_integer;
extern int branch_iterations;
extern int found_local_cuts;
extern int total_num_integer;

extern double **stored_x;
extern int x_rotation;

extern int *nodes_to_fathom;
extern int *integer_indices;
extern int nodes_to_fathom_size;
extern int num_nodes_to_fathom;

extern double *frac_scores;
extern double *frac_values;

extern int insert_counter;
extern int insert_counter2;
extern int rebuild_count;
extern int another_counter;

extern struct nadir *theta;

extern int insert_at_beginning;

extern double prob_width;

extern int cols_have_been_deleted;
extern double obj1_extra_val,obj2_extra_val;

extern int bds_reduced_by_dominated_ideal_pt;
extern int bds_reduced_by_dominated_ideal_segment;
extern int bds_reduced_by_dominated_lb;
extern int bds_reduced_by_infeasibility;
extern int bds_reduced;

extern int constraints_added_after_new_rows;

//extern double *start_x;
//extern double *start_slack;

extern int *global_beg;
extern int *global_varindices; 
extern double *global_values; 
extern int *global_effortlevel;
extern int global_num_starts;
extern int global_startspace;

extern CPXLPptr global_mip;
extern CPXENVptr env_just_solve_mips;

extern double initial_slope;
extern int check_for_stopping_PSA_full;
extern int prev_numsolns;
extern int times_to_run;
extern int there_will_only_be_points, points_only, its_been_only_points;
extern int integer_objective;
extern double smallest_coef;
extern double integer_objective_lb;
extern double multiplying_factor;
extern int integer_bb;
extern double max_range;
extern double x_range, y_range;

extern clock_t start_BB;
extern clock_t finish_BB;
extern double duration_BB;
extern double max_time;
extern int max_nodes;

extern int limit_bd_reduction;
extern double time_processing_nodes;
extern double time_branching;
extern double time_solving_mips;
extern double time_generating_cuts;
extern double time_tightening_variable_bounds;
extern double time_generating_disjunction_cuts;
extern int number_disj_cuts;
extern int number_pareto_branches;

extern double max_time_phase_2;
extern double max_time_bd_red_b4_branching;

extern int objective_space_fathoming;
extern int generate_local_cuts;

extern int only_infeasibility;
extern int check_bound;
extern double time_limit;
extern int check_for_early_PSA_reduce;
extern int check_for_early_PSA_reduce_after_n_iter;
extern int n_iter;
extern int presolve_phase_1;
extern int preprocessing;
extern int presolve_phase_2;
extern int solve_extra_mips_during_cut_callback;
extern int exploit_objective_gaps;
extern int generate_disjunctive_cuts_from_obj_space_disjunctions;

extern int cut_prob_first_row_to_change_index;
extern int cut_prob_second_row_to_change_index;
extern int *cut_prob_bd_indices_disj1;
extern int *cut_prob_bd_indices_disj2;
extern int cut_prob_first_pi_index;
extern int *cut_prob_pi_indices;
extern CPXLPptr cut_problem, lp1_get_cuts_copy;

extern double ob1_offset, ob2_offset;

extern double NW_extreme_x;
extern double NW_extreme_y;
extern double SE_extreme_x;
extern double SE_extreme_y;

extern double *x_separators;
extern double *y_separators;
extern int *separator_dir;
extern int num_separators;

extern double *temp_x_r;
extern double *temp_x_l;

extern user_data *userhandle_current;
extern user_data *userhandle_up;
extern user_data *userhandle_down;
extern double *x_ws;
extern double *x1;
extern double *x2;

extern int show_progress, show_frequency, break_early;
extern double duality_gap_limit;

extern int control_node_selection, allow_changing_control_percentage;
extern double control_percentage, starting_percentage, stopping_percentage;

extern int generate_disjunctive_cuts_from_dominated_columns;
extern int *dominated_indices;
extern int *dominating_indices;
extern int num_dom_cols;

extern int last_cutcallback_seqnum, keep_solving_infeasible_MIPs;
extern double max_preprocessing_time;
extern double max_phase1_time;
extern double max_time_build_dual_bound;

extern double max_time_to_solve_a_mip;

extern int duality_fixing, singleton_columns, dominating_columns, generate_CPLEX_global_cuts_for_several_weighted_sum_objectives, num_iterations_weights;

extern int time_vs_node_lim, node_limit, preprocessing_parameter, times_to_try_heuristic_per_preprocessing_iteration, CPLEX_cuts, hausdorff_or_hypervolume;
extern int preprocessing_method, search_for_extra_solutions_during_preprocessing, stop_PSA_full_early_after_preprocessing, add_local_cuts_during_preprocessing;

extern double percentage_of_integer_variables_to_check_for_probing_during_presolve, extreme_x, extreme_y;

extern FILE *bb_results;

extern clock_t start_struct_timer;
extern clock_t finish_struct_timer;
extern double struct_time;

extern clock_t start_insert_timer;
extern clock_t finish_insert_timer;
extern double insert_time;

//static int CPXPUBLIC 
//nodeoperations (CPXCENVptr env, void *cbdata, int wherefrom,
//		void *cbhandle, int *useraction_p);

//static int CPXPUBLIC
//usersetbranch  (CPXCENVptr env, void *cbdata, int wherefrom,
//		void *cbhandle, int brtype, int sos, int nodes,
//		int bdcnt, const double *nodeest, const int *nodebeg,
//		const int *indices, const char *lu, const int *bd,
//		int *useraction_p);

//static int CPXPUBLIC 
//userselectnode (CPXCENVptr env, void *cbdata, int wherefrom,
//		void *cbhandle, int *nodeid_p,
//		int *useraction_p);

//static int CPXPUBLIC 
//userincumbent (CPXCENVptr env, void *cbdata, int wherefrom,
//	       void *cbhandle, double objval,
//	       double *x, int *isfeas_p,
//	       int *useraction_p);


int computeextremes (CPXENVptr  env,
		     CPXLPptr   lp1,
		     CPXLPptr   lp2,
		     double *obj_coef1,
		     double *obj_coef2);

int computefeassol (CPXENVptr  env,
		    CPXLPptr   lp1,
		    CPXLPptr   lp2,
		    double *obj_coef1,
		    double *obj_coef2,
		    int nPoints,
		    int *nNadir);

int parametricsimplex(CPXLPptr lp);
double mini(double, double);

static void
free_and_null (char **ptr);


struct pointSeg
{double end1_z1; /* end1_z1 < end2_z1   */
  double end1_z2;
  double end2_z1;
  double end2_z2;
  enum {ISOLATED = 1, SEGMENT} type;       /* 1-isolated, 2-segment  */
};

int lp2stdform (CPXENVptr env, CPXLPptr lp);

#endif
