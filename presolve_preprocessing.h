#include "biobjective_bb.h"
#include "callbacks.h"

int run_presolve_phase1	(int cur_numcols, FILE *bb_results, int num_fixed_phase1, int num_singleton_columns1, int num_singleton_columns2,
				double *obj_coef1, double *obj_coef2, CPXCLPptr redlp1, CPXCLPptr redlp2, CPXENVptr env, int *cmatbeg, int *row_indices, 
				double *column_entries, int cur_numrows, CPXLPptr lp1, char *sense2, int fix_it, char *xctype, int *singleton_column_indices2,
				int *their_rows1, int *their_rows2, char *sym, double l, double u, CPXLPptr lp2, int *singleton_column_indices1, int *col_indices,
				double *row_entries, double *lower_bound, double *upper_bound, int suppress_file_output);

int epsilon_constraint_preprocessing(int yet_another_preprocessing_algorithm, clock_t start_presolve, CPXCENVptr env, CPXLPptr lp1, int *indices, 
					double *obj_coef2, double *ub, double *lb, int cur_numcols, CPXLPptr lp2, double split_pt_denom2, clock_t finish_presolve,
					double preprocessing_time, double duration_presolve, FILE *bb_results, int num_mips_to_solve, int numsols, int suppress_file_output);

int weighted_sum_preprocessing(clock_t start_presolve, CPXLPptr lp1, CPXLPptr lp2, CPXENVptr env, int cur_numcols, int cur_numrows, int *indices, 
				double *obj_coef2, int intervals, int run_it, double split_pt_denom, clock_t finish_presolve, double preprocessing_time, 
				double duration_presolve, FILE *bb_results, double *obj_coef1, int suppress_file_output);

int probing2(clock_t phase2_intermediate_time, double time_so_far, clock_t phase2_start_time, int fix_it2, double *lb_pre, double *ub_pre, char *sym, CPXENVptr env,
		CPXLPptr presolve2_lp, CPXLPptr presolve2_mip, int *indices, double *SE_objectives, int insert_check, double *NW_objectives, double *WS_objectives,
		double temp_x, double slope, CPXLPptr lp1, CPXLPptr lp2, int num_bounds_reduced, int cur_numcols, char *xctype, int cur_numrows,
		double *obj_coef1, double *obj_coef2, double *weighted_coefs_, double temp_y);
