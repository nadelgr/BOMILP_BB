/* File created by Nathan Adelgren, Graduate Assisistant at Clemson University.
Started: 9/28/2014 
Finished: N/A
This work serves as part of my doctoral research under the advisement of Dr. Akshay Gupte. 
The goal of this program is to solve Biobjective Mixed-Integer Linear Programs using a BB
based technique employing the Parametric Simplex Algorithm (PSA) and objective-space branching.
(Small pieces of this code may be copied from a BB based code authored by Drs. Banu Soylu and 
Pietro Belotti. I collaborated on this work.)*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <time.h>
#include <pthread.h>
#include <ctype.h>

#include "cplex.h"
#include "max_tree.h"
#include "bb-bicriteria.h"
#include "callbacks.h"
/*#include "presolve_phase1.h"*/
#include "biobjective_bb.h"
#include "minor_functions.h"

/*******************************************************************/
/*              CPLEX VARIABLES                                    */
/*******************************************************************/

CPXENVptr  env=NULL;
CPXLPptr   lp1 = NULL, lp1_get_cuts_copy = NULL;
CPXLPptr   lp2=NULL;
CPXCLPptr   redlp1=NULL;
CPXCLPptr   redlp2=NULL;

double *obj_coef1=NULL;
double *obj_coef2=NULL;
double *weighted_coefs_=NULL;
int cur_numcols ;
int cur_numrows ;
int status = 0;
int status2 = 0;
double objval = 0.;
double *x = NULL;
int numsols = 0;
int real_numsols = 0;

/******************************************************************/

/*******************************************************************/
/*    Other Global Variables - (Don't Define Too Many)             */
/*******************************************************************/

FILE *bb_results;
FILE *pareto_results;
FILE *output2;

split_pt *first_split_pt = NULL;
split_pt *last_split_pt = NULL;

clock_t intermediate_time_start, intermediate_time_finish;
double  duration_presolve;
clock_t start_presolve, finish_presolve;

int presolve = 0;
int intervals = 1;

int constraints_added_after_new_rows = 0;

double *frac_scores;
double *frac_values;

int num_x_to_store = 20;
int heur_limit = 0;
/*double split_pt_denom = 5.;*/
double split_pt_denom = 5.;
double split_pt_denom2 = 12.5;
int yet_another_preprocessing_algorithm = 1;
int num_mips_to_solve = 20;

int obj1_index;
int obj2_index;
int cut_row_index;
int local_cuts = 1;

double obj1_extra_val = 0.;
double obj2_extra_val = 0.;
double initial_slope = 0.;

int cols_have_been_deleted = 0;

/*struct store_it*/
/*{*/
/*    double ratio;*/
/*    int index;*/
/*};*/

/*******************************************************************/

/*void chg_coefs_(CPXCENVptr env, CPXLPptr prob, int *indices, double slope);*/
void chg_coefs2(CPXCENVptr env, CPXLPptr prob, int *indices, double lambda);



/*********************************************************************************************************************** 

	This function is used to remove from the model variables that have been fixed during presolve.
	
***********************************************************************************************************************/

void reduce_prob(CPXCENVptr env, CPXLPptr prob1, CPXLPptr prob2, int *indices, int reduce_global_mip_too)
{
/*	printf("reducing problem\n");*/
	cur_numrows = CPXgetnumrows (env, prob1);
	int cmatbeg[1] = {0};
	double *up_bd = (double *) malloc (cur_numcols*sizeof(double));
	double *lo_bd = (double *) malloc (cur_numcols*sizeof(double));
	double *rhs = (double *) malloc (cur_numrows*sizeof(double));
	int *row_indices = (int *) malloc (cur_numrows*sizeof(double));
	double *column_entries = (double *) malloc (cur_numrows*sizeof(double));
	int nzcnt,surplus;
	
	status = CPXgetlb (env, prob1, lo_bd, 0, cur_numcols-1);
  	status = CPXgetub (env, prob1, up_bd, 0, cur_numcols-1);
  	status = CPXgetrhs (env, prob1, rhs, 0, cur_numrows-1);
  	double val;
  	int i,j,k,m,n,l;
  	int num_changed = 0;
  	
  	for(i=cur_numcols-1;i>=0;i--)
  	{
  		if(lo_bd[i] == up_bd[i])
  		{
  			for(n=num_dom_cols-1;n>=0;n--)
  			{
  				m = 0;
  				if(dominated_indices[n] == i || dominating_indices[n] == i) m = 1;
				if(dominated_indices[n] > i) dominated_indices[n]--;
				if(dominating_indices[n] > i) dominating_indices[n]--;
				if(m)
				{
					for(l=n;l<num_dom_cols-1;l++)
					{
						dominated_indices[l] = dominated_indices[l+1];
						dominating_indices[l] = dominating_indices[l+1];
					}
					num_dom_cols--;
				}
  			}
/*  			printf("upper and lower bounds match for x%d\n",i);*/
  			if(lo_bd[i] != 0.)
  			{
/*  				printf("modifying objective and rhs vals\n");*/
/*  				printf("bd: %lf\t obj1_extra: %lf\t obj2_extra: %lf\t obj_coef: %lf\t obj_coef: %lf\n",lo_bd[i],obj1_extra_val,obj2_extra_val,obj_coef1[i],obj_coef2[i]);*/
	  			obj1_extra_val += obj_coef1[i]*lo_bd[i];
				obj2_extra_val += obj_coef2[i]*lo_bd[i];
/*				printf("bd: %lf\t obj1_extra: %lf\t obj2_extra: %lf\t obj_coef: %lf\t obj_coef: %lf\n",lo_bd[i],obj1_extra_val,obj2_extra_val,obj_coef1[i],obj_coef2[i]);*/
				status = CPXgetcols (env, prob1, &nzcnt, cmatbeg, row_indices, column_entries, cur_numrows, &surplus, i, i);
/*				printf("nonzero count: %d\n",nzcnt);*/
				for(j=0;j<nzcnt;j++)
				{
					k = row_indices[j];
/*					printf("old rhs%d: %lf\n",k,rhs[k]);*/
					rhs[k] -= column_entries[j]*lo_bd[i]; 
/*					printf("new rhs%d: %lf\n",k,rhs[k]);*/
				}
				num_changed++;
			}
/*			printf("deleting column %d\n",i);*/
			status = CPXdelcols (env, prob1, i, i);
/*			printf("status: %d\n",status);*/
			status = CPXdelcols (env, prob2, i, i);
/*			printf("status: %d\n",status);*/
			if(reduce_global_mip_too) status = CPXdelcols (env_just_solve_mips, global_mip, i, i);
			
/*	  		printf("________***************________\n");*/
/*	  		for(m=0;m<num_dom_cols;m++) printf("Column %d dominates column %d\n",dominating_indices[m],dominated_indices[m]);*/
/*	  		printf("________***************________\n");*/
  		}
  	}
  	if(num_changed != 0)
  	{
/*  		status = CPXwriteprob (env, prob1, "myprob2.lp", "LP");*/
/*  		exit(0);*/
  		printf("%d variables were removed that were fixed to non-zero values. Fixing rhs\n", num_changed);
  		printf("extra obj vals: %lf\t %lf\n",obj1_extra_val,obj2_extra_val);
  		status = CPXchgrhs (env, prob1, cur_numrows, indices, rhs);
/*  		printf("status: %d\n",status);*/
  		status = CPXchgrhs (env, prob2, cur_numrows, indices, rhs);
/*  		printf("status: %d\n",status);*/
		if(reduce_global_mip_too) status = CPXchgrhs (env_just_solve_mips, global_mip, cur_numrows, indices, rhs);
  		cols_have_been_deleted = 1;
  	}
  	free(up_bd);
  	free(lo_bd);
  	free(rhs);
  	free(row_indices);
  	free(column_entries);
}

/*********************************************************************************************************************** 

	This function is used for computing the NW-most and SE-most extreme points of the Pareto set 
	in the objective space.
	
***********************************************************************************************************************/

int compute_feas_extremes	(CPXENVptr  env,
				CPXENVptr  env_global,
		    		CPXLPptr   lp1,
		    		CPXLPptr   lp2,
		    		double *obj_coef1,
				double *obj_coef2,
				int *indices)
{
	CPXLPptr lpclone1 = NULL;
	
	/****************** Copy the Problem *****************************/

  	lpclone1 = CPXcloneprob (env, lp1, &status);
  	if ( status ) 
  	{
    		fprintf (stderr, "Failed to clone problem 1.\n");
    		goto TERMINATE;
  	}
  	
  	/******************************************************************/
  	
  	/****************** Set parameters ****************************/

	status = CPXsetintparam (env, CPX_PARAM_REPEATPRESOLVE, 3);
	status = CPXsetintparam (env, CPXPARAM_Emphasis_MIP,CPX_MIPEMPHASIS_BALANCED);
	status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
/*	status = CPXsetdblparam (env, CPX_PARAM_EPGAP, 1.0e-5);*/
	status = CPXsetdblparam (env, CPX_PARAM_TILIM, 180.);
	if ( status ) {
		printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);
		goto TERMINATE;
	}
	
	/******************************************************************/
  	
  	/****************** Compute the extremes ****************************/
	
	x = (double *) malloc (cur_numcols*sizeof(double));
	int i=0;
	int j=0;
	
/*	chg_coefs_(env, lpclone1, indices, -10000000.);*/
	
  	CPXmipopt (env, lpclone1);
/*	CPXpopulate(env,lpclone1);*/
    	status2 = CPXgetstat(env, lpclone1);
/*    	printf("status2: %d\n",status2);*/
    	
    	if(status2 == 107 || status2 == 108)
    	{
    		status = CPXsetintparam (env, CPXPARAM_Emphasis_MIP,CPX_MIPEMPHASIS_OPTIMALITY);
    		status = CPXsetdblparam (env, CPX_PARAM_TILIM, pow(10.,75.));
		if ( status ) {
			printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);
			goto TERMINATE;
		}
		CPXmipopt (env, lpclone1);
    	}
    	
/*    	status2 = CPXgetstat(env, lpclone1);*/
    	
    	
    	if(status2 != 101 && status2 != 102 && status2 != 129 && status2 != 130)
    	{
		if(status2 == 118)
    		{
    			printf("The first test problem is unbounded. Exiting!\n");
    			exit(0);
    			goto TERMINATE;
    		}
    		else if(status2 == 104 || status2 == 105 || status2 == 107 || status2 == 109 || status2 == 111 || status2 == 113)
 			printf("Warning: CPLEX failed to reach an optimal solution. A suboptimal integer-feasible solution is being used.\n");
    	}

	double obvals[2] = {0.,0.};
  	
  	status = CPXgetx (env, lpclone1, obvals, obj1_index, obj2_index);
  	if(status)
  	{
		printf("Failed to get x-values from CPLEX. Line: %d\n",__LINE__);
  		goto TERMINATE;
  	}
  	
  		
  	SE_extreme_x = obvals[0];
  	SE_extreme_y = obvals[1];
  	
  	int num_solns = CPXgetsolnpoolnumsolns (env, lpclone1);
  	
  	prev_numsolns = num_solns;
  	int insert_check = 0;
  	
  	for(j=0;j<num_solns-1;j++)
  	{
  		status = CPXgetsolnpoolx (env, lpclone1, j, x, 0, cur_numcols-1);
/*  		for(i=0;i<total_num_integer;i++)*/
/*		{*/
/*			printf("x%d: %lf\n",integer_indices[i],x[integer_indices[i]]);*/
/*		}*/
	
	      	insert_check = mock_insert(1,x[obj1_index],x[obj2_index],0,0,0,&tree);
	      	if(insert_check)
	      	{
	      		for(i=0;i<cur_numcols;i++)
	      		{
	      			stored_x[x_rotation][i] = x[i];
	/*      		printf("stored val: %lf\n",stored_x[x_rotation][i]);*/
	      		}
	      		x_rotation = (x_rotation + 1) % num_x_to_store;
	      		insert_check = 0;
	      		PSA_full(env_global,NULL,x,NULL,NULL);
	      	}
      	}
      	
 
	chg_coefs_(env, lpclone1, indices, -.000000001);
	CPXmipopt (env, lpclone1);
    	status2 = CPXgetstat(env, lpclone1);
/*    	printf("status2: %d",status2);*/
    	
    	global_mip = CPXcloneprob (env, lpclone1, &status);
  	if ( status ) 
  	{
    		fprintf (stderr, "Failed to clone problem 1.\n");
    		goto TERMINATE;
  	}
    	
    	int num_starts = CPXgetnummipstarts (env, lpclone1);
  	global_num_starts = num_starts;
  	global_startspace = cur_numcols*num_starts;
/*  	printf("number of mip starts: %d\n",num_starts);*/
  	
/*  	printf("first allocation\n");*/
  	global_beg = (int *) malloc ((global_num_starts)*sizeof(int));
	global_varindices = (int *) malloc ((global_startspace)*sizeof(int));
	global_values = (double *) malloc ((global_startspace)*sizeof(double));
	global_effortlevel = (int *) malloc ((global_startspace)*sizeof(int));
	
	int nzcnt = 0, surplus = 0;
	
	status = CPXgetmipstarts (env, lpclone1, &nzcnt, global_beg, global_varindices, 
                           global_values, global_effortlevel, global_startspace,
                           &surplus, 0, global_num_starts-1);
                           
        status = CPXaddmipstarts (env, global_mip, global_num_starts, nzcnt, global_beg, global_varindices,
                           global_values, global_effortlevel, NULL);
                           
        num_starts = CPXgetnummipstarts (env, global_mip);
    	
    	if(status2 != 101 && status2 != 102 && status2 != 129 && status2 != 130)
    	{
		if(status2 == 118)
    		{
    			printf("The second test problem is unbounded. Exiting!\n");
    			exit(0);
    			goto TERMINATE;
    		}
    		else if(status2 == 104 || status2 == 105 || status2 == 107 || status2 == 109 || status2 == 111 || status2 == 113)
 			printf("Warning: CPLEX failed to reach an optimal solution. A suboptimal integer-feasible solution is being used.\n");
    	}

	status = CPXgetx (env, lpclone1, obvals, obj1_index, obj2_index);
  	if(status)
  	{
		printf("Failed to get x-values from CPLEX. Line: %d\n",__LINE__);
		printf("status of second mip: %d\n",status2);
  		goto TERMINATE;
  	}
  	
  	NW_extreme_x = obvals[0];
  	NW_extreme_y = obvals[1];
  	
	num_solns = CPXgetsolnpoolnumsolns (env, lpclone1);
  	
  	if(num_solns >= prev_numsolns) times_to_run = num_solns - prev_numsolns;
  	else times_to_run = num_solns;
  	prev_numsolns = num_solns;
  	
  	
  	for(j=0;j<times_to_run;j++)
  	{
  		status = CPXgetsolnpoolx (env, lpclone1, j, x, 0, cur_numcols-1);
/*  		for(i=0;i<total_num_integer;i++)*/
/*		{*/
/*			printf("x%d: %lf\n",integer_indices[i],x[integer_indices[i]]);*/
/*		}*/
		
	      	insert_check = mock_insert(1,x[obj1_index],x[obj2_index],0,0,0,&tree);
	      	if(insert_check)
	      	{
	      		for(i=0;i<cur_numcols;i++)
	      		{
	      			stored_x[x_rotation][i] = x[i];
	/*      		printf("stored val: %lf\n",stored_x[x_rotation][i]);*/
	      		}
	      		x_rotation = (x_rotation + 1) % num_x_to_store;
	      		insert_check = 0;
	      		PSA_full(env_global,NULL,x,NULL,NULL);
	      	}
      	}
      	
      	prob_width = SE_extreme_x - NW_extreme_x;
      	if(prob_width < .0001) 
      	{
      		printf("Non-conflicting objectives. Only Pareto point is ideal point: (%lf,%lf)\n",NW_extreme_x,NW_extreme_y);
      		finish_BB = clock();
	   	duration_BB = (double)(finish_BB - start_BB) / CLOCKS_PER_SEC;
	   	printf("Total time: %lf\n",duration_BB);
      		exit(0);
      	}

  	/******************************************************************/
  	
  	/****************** Code for TERMINATE *****************************/

  	TERMINATE:

  	free_and_null ((char **) &x);
  	CPXfreeprob(env, &lpclone1);

	status = CPXsetdblparam (env, CPX_PARAM_TILIM, time_limit);
	if ( status ) {
		printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);
		goto TERMINATE;
	}
  	
  	return (status);
  	/******************************************************************/
  	
} /* End of compute_feas_extremes */

/*********************************************************************************************************************** 

	This function is used for changing the weights on the objective function coefficients so that
	they correspond to a particular slope in the objective space.
	
***********************************************************************************************************************/

void chg_coefs_(CPXCENVptr env, CPXLPptr prob, int *indices, double slope)
{
	int i;
/*	printf("slope: %lf\n",slope);*/
/*	printf("negative reciprocal: %lf\n",(-1./slope));*/
	if(slope == 0.) 
	{
/*		printf("%d\n",__LINE__);*/
		for(i=0;i<obj1_index;i++)
		{
			weighted_coefs_[i] = obj_coef2[i];
		}
	}
	else if(!slope)
	{
/*		printf("%d\n",__LINE__);*/
		for(i=0;i<obj1_index;i++)
		{
			weighted_coefs_[i] = obj_coef1[i];
		}
	}
	else
	{
/*		printf("%d\n",__LINE__);*/
		for(i=0;i<obj1_index;i++)
		{
			weighted_coefs_[i] = (-1./slope)*obj_coef2[i] + obj_coef1[i];
/*			printf("negative reciprocal times obj2_%d: %lf\n",i,(-1./slope)*obj_coef2[i]);*/
/*			printf("index %d, obj1 val: %lf, obj2 val: %lf, new val: %lf\n",i,obj_coef1[i],obj_coef2[i],weighted_coefs_[i]);*/
		}	
	}
	weighted_coefs_[obj1_index] = 0.;
	weighted_coefs_[obj2_index] = 0.;
	if(obj1_index < cur_numcols) 
	{
/*		printf("%d\n",__LINE__);*/
		int status = CPXchgobj (env, prob, obj1_index+2, indices, weighted_coefs_);
		if ( status ) {
			printf ("Failed to get change obj coef. Error code %d\n", status);
		}
	}
	else
	{
/*		printf("%d\n",__LINE__);*/
		int status = CPXchgobj (env, prob, obj1_index, indices, weighted_coefs_);
		if ( status ) {
			printf ("Failed to get change obj coef. Error code %d\n", status);
		}
	}
/*	status = CPXwriteprob (env, prob, "myprob2.lp", "LP");*/
}

void chg_coefs2(CPXCENVptr env, CPXLPptr prob, int *indices, double lambda)
{
	int i;
	if(lambda == 0.) 
	{
/*		printf("%d\n",__LINE__);*/
		for(i=0;i<obj1_index;i++)
		{
			weighted_coefs_[i] = obj_coef2[i];
		}
	}
	else if(lambda == 1.)
	{
/*		printf("%d\n",__LINE__);*/
		for(i=0;i<obj1_index;i++)
		{
			weighted_coefs_[i] = obj_coef1[i];
		}
	}
	else
	{
/*		printf("%d\n",__LINE__);*/
		for(i=0;i<obj1_index;i++)
		{
			weighted_coefs_[i] = (1.-lambda)*obj_coef2[i] + lambda*obj_coef1[i];
		}	
	}
	weighted_coefs_[obj1_index] = 0.;
	weighted_coefs_[obj2_index] = 0.;
	if(obj1_index < cur_numcols) 
	{
/*		printf("%d\n",__LINE__);*/
		int status = CPXchgobj (env, prob, obj1_index+2, indices, weighted_coefs_);
		if ( status ) {
			printf ("Failed to get change obj coef. Error code %d\n", status);
		}
	}
	else
	{
/*		printf("%d\n",__LINE__);*/
		int status = CPXchgobj (env, prob, obj1_index, indices, weighted_coefs_);
		if ( status ) {
			printf ("Failed to get change obj coef. Error code %d\n", status);
		}
	}
/*	status = CPXwriteprob (env, prob, "myprob2.lp", "LP");*/
}

int add_interval(double val, interval *intrvl, int count)
{
/*	printf("adding interval\n");*/
/*	printf("current interval properties: %lf,%lf,%d,%d\n",intrvl->lb,intrvl->ub,intrvl->dir,intrvl->miss_count);*/
/*	if(intrvl->next == NULL) printf("next is null\n");*/
/*	if(!intrvl->next) printf("no next\n");*/
/*	printf("next interval properties: %lf,%lf,%d,%d\n",intrvl->next->lb,intrvl->next->ub,intrvl->next->dir,intrvl->next->miss_count);*/
	interval *temp = intrvl;
	while(temp && temp->next && val > temp->ub) temp = temp->next;
	interval *new = (struct interval*) malloc( sizeof( struct interval ) );
	new->lb = val;
	new->ub = temp->ub;
	new->next = temp->next;
	new->dir = temp->dir;
	new->miss_count = temp->miss_count;
	temp->ub = val;
	temp->next = new;
	return(count+1);
}

/*********************************************************************************************************************** 

	Here we perform the BB.
	
***********************************************************************************************************************/

int main(int argc, char **argv)
{
	struct_time = 0.;
	struct_time2 = 0.;
	insert_time = 0.;
	insert_time2 = 0.;
	accum = 0.;
	srand(1);
	if(argc < 3) 
	{
		printf("Error! Not enough input files.\n\nPlease include 2 input LP files having the same feasible set and different objective functions. Also ensure that in both files variables appear for the first time in the same order.\n\nAdditionally, please include a binary flag: 0 if you do NOT wish to have the final Pareto set displayed in MATLAB readable format, or 1 if you DO wish to display this information.\n\n");
		exit(0);
	}
	int print_Pareto = 0;
	
	/*************************************************************/
	/*          Process command line flags                       */
	/*************************************************************/
	
	int i = 3;
	int opened_bb_results = 0;
	int suppress_file_output = 0;
	int write_tree_results = 0;
	int matlab = 1;
	while(i < argc)
    {
        if(argv[i][0] != '-')
        {
            printf("Invalid Command Line Argument (missing '-'). Exiting.\n");
	        printf("Got %s\n", argv[i]);
            exit(0);
        }
        else
        {
            if(!strcmp(argv[i],"-presolve1"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    presolve_phase_1 = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    presolve_phase_1 = 0;
                }
                else
                {
                    printf("Invalid value for flag '-presolve1', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-preprocessing"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    preprocessing = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    preprocessing = 0;
                }
                else
                {
                    printf("Invalid value for flag '-preprocessing', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-ppType"))
            {
                i++;
                if(atoi(argv[i]) < 0 || atoi(argv[i]) > 3) printf("Invalid value for flag '-ppType', ignoring. Valid values are 0 through 3.\n");
                else preprocessing_method = atoi(argv[i]);
		        i++;
            }
            else if(!strcmp(argv[i],"-ppParam"))
            {
                i++;
                if(atoi(argv[i]) < 0 ) printf("Invalid value for flag '-ppParam', ignoring. Valid values are nonnegative integers.\n");
                else preprocessing_parameter = atoi(argv[i]);
		        i++;
            }
            else if(!strcmp(argv[i],"-ppExtra"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    search_for_extra_solutions_during_preprocessing = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    search_for_extra_solutions_during_preprocessing = 0;
                }
                else
                {
                    printf("Invalid value for flag '-ppExtra', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-ppCuts"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    add_local_cuts_during_preprocessing = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    add_local_cuts_during_preprocessing = 0;
                }
                else
                {
                    printf("Invalid value for flag '-ppCuts', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-ppGenetic"))
            {
                i++;
                if(atoi(argv[i]) < 0 ) printf("Invalid value for flag '-ppGenetic', ignoring. Valid values are nonnegative integers.\n");
                else times_to_try_heuristic_per_preprocessing_iteration = atoi(argv[i]);
		        i++;
            }
            else if(!strcmp(argv[i],"-presolve2"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    presolve_phase_2 = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    presolve_phase_2 = 0;
                }
                else
                {
                    printf("Invalid value for flag '-presolve2', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-osFathom"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    objective_space_fathoming = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    objective_space_fathoming = 0;
                }
                else
                {
                    printf("Invalid value for flag '-osFathom', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-cutStrength"))
            {
                i++;
                if(atoi(argv[i]) < 0 || atoi(argv[i]) > 2) printf("Invalid value for flag '-cutStrength', ignoring. Valid values are 0 through 2.\n");
                else CPLEX_cuts = atoi(argv[i]);
		        i++;
            }
            else if(!strcmp(argv[i],"-cutMulti"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    generate_CPLEX_global_cuts_for_several_weighted_sum_objectives = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    generate_CPLEX_global_cuts_for_several_weighted_sum_objectives = 0;
                }
                else
                {
                    printf("Invalid value for flag '-cutMulti', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-cutNum"))
            {
                i++;
                if(atoi(argv[i]) < 1) printf("Invalid value for flag '-cutNum', ignoring. Valid values are positive integers.\n");
                else num_iterations_weights = atoi(argv[i]);
		        i++;
            }
            else if(!strcmp(argv[i],"-cutLocal"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    generate_local_cuts = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    generate_local_cuts = 0;
                }
                else
                {
                    printf("Invalid value for flag '-cutLocal', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-intObj"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    integer_bb = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    integer_bb = 0;
                }
                else
                {
                    printf("Invalid value for flag '-intObj', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-bndRed"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    bd_reduction_during_branching = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    bd_reduction_during_branching = 0;
                }
                else
                {
                    printf("Invalid value for flag '-bndRed', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-bndInf"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    only_infeasibility = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    only_infeasibility = 0;
                }
                else
                {
                    printf("Invalid value for flag '-bndInf', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-bndFull"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    check_bound = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    check_bound = 0;
                }
                else
                {
                    printf("Invalid value for flag '-bndFull', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-bndLimit"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    limit_bd_reduction = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    limit_bd_reduction = 0;
                }
                else
                {
                    printf("Invalid value for flag '-bndLimit', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-psaEarly"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    check_for_early_PSA_reduce = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    check_for_early_PSA_reduce = 0;
                }
                else
                {
                    printf("Invalid value for flag '-psaEarly', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-psaEarlyN"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    check_for_early_PSA_reduce_after_n_iter = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    check_for_early_PSA_reduce_after_n_iter = 0;
                }
                else
                {
                    printf("Invalid value for flag '-psaEarlyN', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-psaIter"))
            {
                i++;
                if(atoi(argv[i]) < 1) printf("Invalid value for flag '-psaIter', ignoring. Valid values are positive integers.\n");
                else n_iter = atoi(argv[i]);
		        i++;
            }
            else if(!strcmp(argv[i],"-psaStop"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    stop_PSA_full_early_after_preprocessing = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    stop_PSA_full_early_after_preprocessing = 0;
                }
                else
                {
                    printf("Invalid value for flag '-psaStop', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-extraMIPs"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    solve_extra_mips_during_cut_callback = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    solve_extra_mips_during_cut_callback = 0;
                }
                else
                {
                    printf("Invalid value for flag '-extraMIPs', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-gaps"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    exploit_objective_gaps = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    exploit_objective_gaps = 0;
                }
                else
                {
                    printf("Invalid value for flag '-gaps', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-osDisj"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    generate_disjunctive_cuts_from_obj_space_disjunctions = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    generate_disjunctive_cuts_from_obj_space_disjunctions = 0;
                }
                else
                {
                    printf("Invalid value for flag '-osDisj', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-showProgress"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    show_progress = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    show_progress = 0;
                }
                else
                {
                    printf("Invalid value for flag '-showProgress', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-showFreq"))
            {
                i++;
                if(atoi(argv[i]) < 1) printf("Invalid value for flag '-showFreq', ignoring. Valid values are positive integers.\n");
                else show_frequency = atoi(argv[i]);
		        i++;
            }
            else if(!strcmp(argv[i],"-showProgress"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    show_progress = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    show_progress = 0;
                }
                else
                {
                    printf("Invalid value for flag '-showProgress', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-hausHyper"))
            {
                i++;
                if(toupper(argv[i][0]) == 'A')
                {
                    hausdorff_or_hypervolume = 1;
                }
                else if(toupper(argv[i][0]) == 'Y')
                {
                    hausdorff_or_hypervolume = 2;
                }
                else
                {
                    printf("Invalid value for flag '-hausHyper', ignoring. Valid values are 'A' (for h(A)usdorff) and 'Y' (for h(Y)pervolume).\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-nodeSel"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    control_node_selection = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    control_node_selection = 0;
                }
                else
                {
                    printf("Invalid value for flag '-nodeSel', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-nodeSelPerc")) 
            {
                i++;
		        if (argv[i] != NULL && atof(argv[i]) > 0 && atof(argv[i]) < 1) 
		        {
		            control_percentage = atof(argv[i]);
		            i++;
		        } 
		        else 
		        {
                    printf("Invalid value for flag '-nodeSelPerc', ignoring. Valid values are reals in (0,1).\n");
		        }
            }
            else if(!strcmp(argv[i],"-nodeSelVar"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    allow_changing_control_percentage = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    allow_changing_control_percentage = 0;
                }
                else
                {
                    printf("Invalid value for flag '-nodeSelVar', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-nodeSelPerc1")) 
            {
                i++;
		        if (argv[i] != NULL && atof(argv[i]) > 0 && atof(argv[i]) < 1) 
		        {
		            starting_percentage = atof(argv[i]);
		            i++;
		        } 
		        else 
		        {
                    printf("Invalid value for flag '-nodeSelPerc1', ignoring. Valid values are reals in (0,1).\n");
		        }
            }
            else if(!strcmp(argv[i],"-nodeSelPerc2")) 
            {
                i++;
		        if (argv[i] != NULL && atof(argv[i]) > 0 && atof(argv[i]) < 1) 
		        {
		            stopping_percentage = atof(argv[i]);
		            i++;
		        } 
		        else 
		        {
                    printf("Invalid value for flag '-nodeSelPerc2', ignoring. Valid values are reals in (0,1).\n");
		        }
            }
            else if(!strcmp(argv[i],"-domCol"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    dominating_columns = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    dominating_columns = 0;
                }
                else
                {
                    printf("Invalid value for flag '-domCol', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-domColCuts"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    generate_disjunctive_cuts_from_dominated_columns = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    generate_disjunctive_cuts_from_dominated_columns = 0;
                }
                else
                {
                    printf("Invalid value for flag '-domColCuts', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-dualFix"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    duality_fixing = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    duality_fixing = 0;
                }
                else
                {
                    printf("Invalid value for flag '-dualFix', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-singCol"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    singleton_columns = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    singleton_columns = 0;
                }
                else
                {
                    printf("Invalid value for flag '-singCol', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-limitType"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    time_vs_node_lim = 0;
                }
                else if(toupper(argv[i][0]) == 'N')
                {
                    time_vs_node_lim = 1;
                }
                else
                {
                    printf("Invalid value for flag '-limitType', ignoring. Valid values are 'T' (for (T)ime) and 'N' (for (N)ode).\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-proveInfeas"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    keep_solving_infeasible_MIPs = 0;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    keep_solving_infeasible_MIPs = 1;
                }
                else
                {
                    printf("Invalid value for flag '-proveInfeas', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else if(!strcmp(argv[i],"-bbTime")) 
            {
                i++;
		        if (argv[i] != NULL && atof(argv[i]) > 0) 
		        {
		            max_time = atof(argv[i]);
		            i++;
		        } 
		        else 
		        {
                    printf("Invalid value for flag '-bbTime', ignoring. Valid values are positive reals.\n");
		        }
            }
            else if(!strcmp(argv[i],"-bbNodes"))
            {
                i++;
                if(atoi(argv[i]) < 1) printf("Invalid value for flag '-bbNodes', ignoring. Valid values are positive integers.\n");
                else max_nodes = atoi(argv[i]);
		        i++;
            }
            else if(!strcmp(argv[i],"-mipTime")) 
            {
                i++;
		        if (argv[i] != NULL && atof(argv[i]) > 0) 
		        {
		            time_limit = atof(argv[i]);
		            i++;
		        } 
		        else 
		        {
                    printf("Invalid value for flag '-mipTime', ignoring. Valid values are positive reals.\n");
		        }
            }
            else if(!strcmp(argv[i],"-mipNode"))
            {
                i++;
                if(atoi(argv[i]) < 1) printf("Invalid value for flag '-mipNode', ignoring. Valid values are positive integers.\n");
                else node_limit = atoi(argv[i]);
		        i++;
            }
            else if(!strcmp(argv[i],"-presolve1time")) 
            {
                i++;
		        if (argv[i] != NULL && atof(argv[i]) > 0) 
		        {
		            max_phase1_time = atof(argv[i]);
		            i++;
		        } 
		        else 
		        {
                    printf("Invalid value for flag '-presolve1time', ignoring. Valid values are positive reals.\n");
		        }
            }
            else if(!strcmp(argv[i],"-presolve2time")) 
            {
                i++;
		        if (argv[i] != NULL && atof(argv[i]) > 0) 
		        {
		            max_time_phase_2 = atof(argv[i]);
		            i++;
		        } 
		        else 
		        {
                    printf("Invalid value for flag '-presolve2time', ignoring. Valid values are positive reals.\n");
		        }
            }
            else if(!strcmp(argv[i],"-bndRedTime")) 
            {
                i++;
		        if (argv[i] != NULL && atof(argv[i]) > 0) 
		        {
		            max_time_bd_red_b4_branching = atof(argv[i]);
		            i++;
		        } 
		        else 
		        {
                    printf("Invalid value for flag '-bndRedTime', ignoring. Valid values are positive reals.\n");
		        }
            }
            else if(!strcmp(argv[i],"-ppTime")) 
            {
                i++;
		        if (argv[i] != NULL && atof(argv[i]) > 0) 
		        {
		            max_preprocessing_time = atof(argv[i]);
		            i++;
		        } 
		        else 
		        {
                    printf("Invalid value for flag '-ppTime', ignoring. Valid values are positive reals.\n");
		        }
            }
            else if(!strcmp(argv[i],"-dualGapLimit")) 
            {
                i++;
		        if (argv[i] != NULL && atof(argv[i]) > 0) 
		        {
		            duality_gap_limit = atof(argv[i]);
		            i++;
		        } 
		        else 
		        {
                    printf("Invalid value for flag '-dualGapLimit', ignoring. Valid values are positive reals.\n");
		        }
            }
            else if(!strcmp(argv[i],"-dualTime")) 
            {
                i++;
		        if (argv[i] != NULL && atof(argv[i]) > 0) 
		        {
		            max_time_build_dual_bound = atof(argv[i]);
		            i++;
		        } 
		        else 
		        {
                    printf("Invalid value for flag '-dualTime', ignoring. Valid values are positive reals.\n");
		        }
            }
            else if(!strcmp(argv[i],"-probePerc")) 
            {
                i++;
		        if (argv[i] != NULL && atof(argv[i]) >= 0. && atof(argv[i]) <= 100.) 
		        {
		            percentage_of_integer_variables_to_check_for_probing_during_presolve = atof(argv[i]);
		            i++;
		        } 
		        else 
		        {
                    printf("Invalid value for flag '-probePerc', ignoring. Valid values are positive reals.\n");
		        }
            }
            else if(!strcmp(argv[i],"-out"))
            {
                i++;
                if(strcmp(argv[i],"none"))
                {
                    if ((bb_results = fopen (argv[i], "a+"))==NULL)
	            	{
	            		printf("could not open output file, exiting\n");
	              		fprintf (stderr, "could not open output file, exiting\n");
	              		exit(1);
	            	}
	        	}
	        	else suppress_file_output = 1;
	        	opened_bb_results = 1;
                i++;
            }
            else if(!strcmp(argv[i],"-print"))
            {
                i++;
                if(!strcmp(argv[i],"T"))
                {
                    print_Pareto = 1;
                }
                else if(!strcmp(argv[i],"F"))
                {
                    print_Pareto = 0;
                }
                else 
                {
                    print_Pareto = 2;
                    if ((pareto_results = fopen (argv[i], "w"))==NULL)
                	{
                		printf("could not open output file, exiting\n");
                  		fprintf (stderr, "could not open output file, exiting\n");
                  		exit(1);
                	}
                }
                i++;
            }
            else if(!strcmp(argv[i],"-matlab"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    matlab = 1;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    matlab = 0;
                }
                else
                {
                    printf("Invalid value for flag '-matlab', ignoring. Valid values are 'T' and 'F'.\n");
                }
                i++;
            }
            else
            {
                printf("'%s' is an invalid command line flag. Exiting!\nFor a list of valid flags, see 'README.txt'\n", argv[i]);
                exit(0);
            }
        }
     }

	/* initialize Cplex environment *********************************/
	
  	env = CPXopenCPLEX(&status);
  	env_just_solve_mips = CPXopenCPLEX(&status);

  	if(env==NULL)
  	{
    		char errmsg[1024];
    		printf("CPXopenCPLEX, Couldn't open the CPLEX environment, error code %d\n", status);
    		CPXgeterrorstring (env, status, errmsg);
    		printf ("%s", errmsg);
    		goto TERMINATE;
  	}

    	/******************************************************************/
    	
    	/********* Open file for storing results ************************/
    	
    	if(!suppress_file_output)
    	{
        	if(!opened_bb_results)
        	{
	        	if ((bb_results = fopen ("bb_results.txt", "a+"))==NULL)
	        	{
	        		printf("could not open output file, exiting\n");
	          		fprintf (stderr, "could not open output file, exiting\n");
	          		exit(1);
	        	}
        	}
        	
        	if(write_tree_results)
        	{
        	    if ((output2 = fopen ("tree_results.txt", "w"))==NULL)
	        	{
	        		printf("could not open output file, exiting\n");
	          		fprintf (stderr, "could not open output file, exiting\n");
	          		exit(1);
	        	}
	        }
	    }
    	
    	/******************************************************************/

  	/************* Set to 1 thread **********************************/
  	
	status = CPXsetintparam (env, CPX_PARAM_THREADS, 1);
	if ( status ) {
		printf ("Failure to set threads to 1, error %d.\n",status);
		goto TERMINATE;
	}
	
	status = CPXsetintparam (env_just_solve_mips, CPX_PARAM_THREADS, 1);
	if ( status ) {
		printf ("Failure to set threads to 1, error %d.\n",status);
		goto TERMINATE;
	}
  	
  	/******************************************************************/
  	
  	/************* Set any desired CPLEX parameters here **************/
  	
	status = CPXsetdblparam (env, CPXPARAM_MIP_Pool_AbsGap, 0.0);
	if ( status ) {
		printf ("Failed to set solution pool gap to 0, error %d.\n",status);
		goto TERMINATE;
	}
	status = CPXsetdblparam (env, CPXPARAM_MIP_Pool_RelGap, 0.0);
	if ( status ) {
		printf ("Failed to set solution pool gap to 0, error %d.\n",status);
		goto TERMINATE;
	}
	
	status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
	if ( status ) {
		printf ("Failed to set solution pool gap to 0, error %d.\n",status);
		goto TERMINATE;
	}
	
	status = CPXsetintparam (env_just_solve_mips, CPX_PARAM_SCRIND, CPX_OFF);
	if ( status ) {
		printf ("Failed to set solution pool gap to 0, error %d.\n",status);
		goto TERMINATE;
	}
	
	status = CPXsetintparam (env_just_solve_mips, CPX_PARAM_POPULATELIM, 200);
	if ( status ) {
		printf ("Failed to set populate limit to 200, error %d.\n",status);
		goto TERMINATE;
	}
	
	status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	if ( status )  goto TERMINATE;
	
	status = CPXsetintparam (env, CPX_PARAM_PRELINEAR, 0);
	if ( status )  goto TERMINATE;
	
	status = CPXsetintparam (env, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);
	
	status = CPXsetintparam (env, CPXPARAM_Emphasis_MIP,CPX_MIPEMPHASIS_OPTIMALITY); //<-- options are FEASIBILITY, BALANCED, OPTIMALITY, ...
	
/*	status = CPXsetintparam (env, CPX_PARAM_ADVIND, 0);*/
/*		  	if(status) printf("(%d) there was an error\n",__LINE__);*/
	
/* 	status = CPXsetintparam (env, CPX_PARAM_BNDSTRENIND, 0);*/
/*	if(status) printf("(%d) there was an error\n",__LINE__);*/

  	/******************************************************************/
  	
  	/** Create test Problems *******************/

  	lp1 = CPXcreateprob (env,&status,argv[1]);
  	if(lp1==NULL) 
  	{
    		printf("CPXcreateprob, Failed to create LP1, error code %d\n", status);
    		goto TERMINATE;
    	}

  	lp2 = CPXcreateprob (env,&status,argv[2]);
  	if(lp2==NULL) 
  	{
    		printf("CPXcreateprob, Failed to create LP2, error code %d\n", status);
    		goto TERMINATE;
    	}
    	
	/******************************************************************/

  	/** read the testP1 and testP2 from the model file **********/

  	status = CPXreadcopyprob(env,lp1,argv[1],NULL);

  	if ( status ) 
  	{
    		printf ("Could not read input 1, exiting. Error code: %d\n", status);
    		goto TERMINATE;
    	}

  	status = CPXreadcopyprob(env,lp2,argv[2],NULL);

  	if ( status ) 
  	{
    		printf ("Could not read input 2, exiting. Error code: %d\n", status);
    		goto TERMINATE;
    	}
	/******************************************************************/
	
  	/******************************************************************
  		Convert greater than or equal to constraints
  		to less than or equal to.
  	******************************************************************/

	start_BB = clock();
	
	cur_numcols = CPXgetnumcols (env, lp1);
	cur_numrows = CPXgetnumrows (env, lp1);
	
	char *row_sense;
    	row_sense = (char *) malloc ((cur_numrows)*sizeof(char));
    	status = CPXgetsense (env, lp1, row_sense, 0, cur_numrows-1);
    	
    	int j = 0, equalities = 0;
    	for(i = 0;i<cur_numrows;i++)
    	{
    		if(row_sense[i] == 'G')
    		{
    			double coef = 0.;
    			double rhs[1] = {0.};
    			char new_sense[1] = {'L'};
    			for(j=0;j<cur_numcols;j++)
    			{
    				status = CPXgetcoef (env, lp1, i, j, &coef);
    				if(coef != 0.)
    				{
    					status = CPXchgcoef (env, lp1, i, j, -coef);
    					status = CPXchgcoef (env, lp2, i, j, -coef);
    				}
    			}
    			status = CPXgetrhs (env, lp1, rhs, i, i);
    			rhs[0] = -rhs[0];
    			status = CPXchgrhs (env, lp1, 1, &i, rhs);
    			status = CPXchgsense (env, lp1, 1, &i, new_sense);
    			status = CPXchgrhs (env, lp2, 1, &i, rhs);
    			status = CPXchgsense (env, lp2, 1, &i, new_sense);
    		}
    		if(!equalities && row_sense[i] == 'E')
    		{
    			equalities = 1;
/*    			generate_disjunctive_cuts_from_dominated_columns = 0;*/
    		}
    	}

	/******************************************************************
  		Perform CPLEX primal presolve.
  	******************************************************************/
	
	status = CPXpresolve (env, lp1, CPX_ALG_NONE) || CPXpresolve (env, lp2, CPX_ALG_NONE);
/*	printf("status: %d\n",status);*/
	
	status = CPXgetredlp (env, lp1, &redlp1) || CPXgetredlp (env, lp2, &redlp2);
/*	printf("status: %d\n",status);*/

/*	status = CPXwriteprob (env, redlp1, "myprob.lp", "LP");*/
/*	status = CPXwriteprob (env, redlp2, "myprob1.lp", "LP");*/
/*	exit(0);*/
	
	if(redlp1)
	{
		cur_numcols = CPXgetnumcols (env, redlp1);
	  	cur_numrows = CPXgetnumrows (env, redlp1);
/*	  	double ob1_offset = 0., ob2_offset = 0.;*/
		status = CPXgetobjoffset(env, redlp1, &ob1_offset );
		status = CPXgetobjoffset(env, redlp2, &ob2_offset );
		printf("ob1 offset: %lf, ob2 offset: %lf\n",ob1_offset,ob2_offset);

		obj1_extra_val += ob1_offset;
		obj2_extra_val += ob2_offset;
		
/*		lp1 = CPXcloneprob (env, redlp1, &status);*/
/*		lp2 = CPXcloneprob (env, redlp2, &status);*/
		lp1 = (CPXLPptr) redlp1;
		lp2 = (CPXLPptr) redlp2;
  	}
  	else
  	{
  		bd_reduction_during_branching = 0;
  		cur_numcols = CPXgetnumcols (env, lp1);
	  	cur_numrows = CPXgetnumrows (env, lp1);
  	}
  	
  	pass_the_lps(env, lp1, lp2);
  	
  	/******************************************************************/
  	
  	/** get the objective coefficients ********************************/
  	
	int status2, insert_check;
	double ob_v1 = obj1_extra_val;
	double ob_v2 = obj2_extra_val;
	int run_it = 0;

  	obj_coef1 = (double *) malloc ((cur_numcols+2)*sizeof(double));
  	weighted_coefs_ = (double *) malloc ((cur_numcols+2)*sizeof(double));
  	weighted_coefs_[cur_numcols] = 0.;
  	weighted_coefs_[cur_numcols+1] = 0.;
  	
/*  	printf("size of obj_coef arrays: %d\n",cur_numcols+2);*/

 	status = CPXgetobj(env, lp1, obj_coef1, 0, cur_numcols -1);

  	if ( status ) 
  	{
    		printf ("CPXgetobj, Failed to read the objective coef. of reduced myBOMIP1, error code %d\n", status);
    		goto TERMINATE;
    	}

  	obj_coef2 = (double *) malloc ((cur_numcols+2)*sizeof(double));

  	status = CPXgetobj(env,lp2, obj_coef2, 0, cur_numcols -1);

  	if ( status ) 
  	{
    		printf ("CPXgetobj, Failed to read the objective coef. of reduced myBOMIP2, error code %d\n", status);
    		goto TERMINATE;
    	}
    	
    	provide_coefs(obj_coef1,obj_coef2,cur_numcols);
    	
    	/******************************************************************/
    	
    	/******************************************************************
  		Prepare for presolve phase 1.
  	******************************************************************/
    	   	
    	int *indices = NULL;
    	indices = (int *) malloc (max(5*cur_numcols+2,5*cur_numrows+2)*sizeof(int));
	for(i=0;i<max(5*cur_numcols+2,5*cur_numrows+2);i++) indices[i] = i;
	
	double *first_x = (double *) malloc ((cur_numcols)*sizeof(double));
    	
    	char *xctype;
    	xctype = (char *) malloc ((cur_numcols+2)*sizeof(char));
    	
    	status = CPXgetctype (env, lp1, xctype, 0, cur_numcols-1);
    	if ( status ) 
  	{
    		printf ("(%d) CPXgetctype, Failed to get variable types, error code %d\n", __LINE__,status);
    		goto TERMINATE;
    	}

	finish_presolve = clock();
	duration_presolve = (double)(finish_presolve - start_BB) / CLOCKS_PER_SEC;
  	printf("*******************\nTime spent prior to presolve phase 1: %lf\n************************\n",
									duration_presolve);
    	
    	/*************************************************************************************
    	 ____  ____  ____      ____   __   __    _  _  ____     __   __ _  ____ 
	(  _ \(  _ \(  __)___ / ___) /  \ (  )  / )( \(  __)   /  \ (  ( \(  __)
	 ) __/ )   / ) _)(___)\___ \(  O )/ (_/\\ \/ / ) _)   (  O )/    / ) _) 
	(__)  (__\_)(____)    (____/ \__/ \____/ \__/ (____)   \__/ \_)__)(____)
	Variable Fixing Based on Problem Structure
	
	*************************************************************************************/
	
	if(redlp1)
	{
		status = CPXgetobj(env, redlp1, obj_coef1, 0, cur_numcols -1);
	  	status = CPXgetobj(env, redlp2, obj_coef2, 0, cur_numcols -1);
  	}
  	else
  	{
  		status = CPXgetobj(env, lp1, obj_coef1, 0, cur_numcols -1);
	  	status = CPXgetobj(env, lp2, obj_coef2, 0, cur_numcols -1);
  	}
	
	int cmatbeg[1] = {0};
	char sym[1] = {'U'};
	double l,u;
	int *row_indices = (int *) malloc (cur_numrows*sizeof(double));
	int *col_indices = (int *) malloc (cur_numcols*sizeof(double));
	double *column_entries = (double *) malloc (cur_numrows*sizeof(double));
	double *row_entries = (double *) malloc (cur_numcols*sizeof(double));
	char *sense2 = (char *) malloc (cur_numrows*sizeof(char));
	if(redlp1) status = CPXgetsense (env, redlp1, sense2, 0, cur_numrows-1);
	else status = CPXgetsense (env, lp1, sense2, 0, cur_numrows-1);
	int nzcnt,surplus;
	int fix_it = 0;
	int num_fixed_phase1 = 0;
	int *singleton_column_indices1 = (int *) malloc (cur_numcols*sizeof(double));
	int *their_rows1 = (int *) malloc (cur_numcols*sizeof(double));
	int *singleton_column_indices2 = (int *) malloc (cur_numcols*sizeof(double));
	int *their_rows2 = (int *) malloc (cur_numcols*sizeof(double));
	int num_singleton_columns1 = 0;
	int num_singleton_columns2 = 0;
	
	double *lower_bound = (double *) malloc (cur_numcols*sizeof(double));
	double *upper_bound = (double *) malloc (cur_numcols*sizeof(double));

	if(redlp1)
	{
		status = CPXgetlb (env, redlp1, lower_bound, 0, cur_numcols-1);
	  	status = CPXgetub (env, redlp1, upper_bound, 0, cur_numcols-1);
  	}
  	else
	{
		status = CPXgetlb (env, lp1, lower_bound, 0, cur_numcols-1);
	  	status = CPXgetub (env, lp1, upper_bound, 0, cur_numcols-1);
  	}
	
/*	printf("variable fixing phase 1:\n");*/
	if(presolve_phase_1)
	{
		status = run_presolve_phase1	(cur_numcols, bb_results, num_fixed_phase1, num_singleton_columns1, num_singleton_columns2, obj_coef1,		
						 obj_coef2, redlp1, redlp2, env, cmatbeg, row_indices, column_entries, cur_numrows, lp1, sense2, fix_it,	
						 xctype, singleton_column_indices2, their_rows1, their_rows2, sym, l, u, lp2, singleton_column_indices1,
						 col_indices, row_entries, lower_bound, upper_bound, suppress_file_output);
		if(status) goto TERMINATE;
	}
	else if(!suppress_file_output)
	{
		fprintf(bb_results,"%d\t",0);
		fprintf(bb_results,"%d\t",0);
		fprintf(bb_results,"%d\t",0);
		fprintf(bb_results,"%d\t",0);
		fprintf(bb_results,"%d\t",0);
		fprintf(bb_results,"%d\t",0);
	}
	
	free(row_indices);
	free(col_indices);
	free(column_entries);
	free(row_entries);
	free(sense2);
	free(singleton_column_indices1);
	free(their_rows1);
	free(singleton_column_indices2);
	free(their_rows2);
	free(lower_bound);
	free(upper_bound);
	
	if(redlp1)
	{
/*		lp1 = CPXcloneprob (env, redlp1, &status);*/
/*		lp2 = CPXcloneprob (env, redlp2, &status);*/
		lp1 = (CPXLPptr) redlp1;
		lp2 = (CPXLPptr) redlp2;
	}
	
/*	for(i=0;i<num_dom_cols;i++) printf("Column %d dominates column %d\n",dominating_indices[i],dominated_indices[i]);*/
	
	reduce_prob(env,lp1,lp2,indices,0);
	
/*	for(i=0;i<num_dom_cols;i++) printf("Column %d dominates column %d\n",dominating_indices[i],dominated_indices[i]);*/
	
/*	status = CPXwriteprob (env, lp1, "myprob3.lp", "LP");*/
/*	status = CPXwriteprob (env, lp2, "myprob4.lp", "LP");*/
/*	exit(0);*/

  	cur_numcols = CPXgetnumcols (env, lp1);
  	cur_numrows = CPXgetnumrows (env, lp1);

	pass_the_lps(env, lp1, lp2);
	
	status = CPXgetobj(env, lp1, obj_coef1, 0, cur_numcols -1);
/*	printf("status: %d\n",status);*/
  	status = CPXgetobj(env,lp2, obj_coef2, 0, cur_numcols -1);
/*  	printf("status: %d\n",status);*/
  	free_coefs();
    	provide_coefs(obj_coef1,obj_coef2,cur_numcols);
    	weighted_coefs_[cur_numcols] = 0.;
  	weighted_coefs_[cur_numcols+1] = 0.;
  	
  	finish_presolve = clock();
	duration_presolve = (double)(finish_presolve - start_BB) / CLOCKS_PER_SEC;
  	printf("*******************\nTotal time spent until end of presolve phase 1: %lf\n************************\n",
		  									duration_presolve);
	if(!suppress_file_output) fprintf(bb_results,"%lf\t",duration_presolve);
  	
/*  	if(cur_numcols < 40) yet_another_preprocessing_algorithm = 0;*/
    	
    	/********************* Determine whether or not either objective depends only on integer variables *************************/
    	
    	status = CPXgetctype (env, lp1, xctype, 0, cur_numcols-1);
    	if ( status ) 
  	{
    		printf ("(%d) CPXgetctype, Failed to get variable types, error code %d\n", __LINE__,status);
    		goto TERMINATE;
    	}
    	
/*    	provide_xctype(env, xctype,cur_numcols);*/
	
	int obj1_contains_only_ints = 1, obj2_contains_only_ints = 1;
	double smallest_obj1_coef = 10000000., smallest_obj2_coef = 10000000.;
	if(integer_bb)
	{  									
		for(i=0;i<cur_numcols;i++)
		{
			if(xctype[i] == 'C' && fabs(obj_coef1[i]) > .0000001) obj1_contains_only_ints = 0;
			if(xctype[i] == 'C' && fabs(obj_coef2[i]) > .0000001) obj2_contains_only_ints = 0;
			if(xctype[i] == 'B' || xctype[i] == 'I')
			{
				if(fabs(obj_coef1[i]) > .0000001 && fabs(obj_coef1[i]) < smallest_obj1_coef) smallest_obj1_coef = fabs(obj_coef1[i]);
				if(fabs(obj_coef2[i]) > .0000001 && fabs(obj_coef2[i]) < smallest_obj2_coef) smallest_obj2_coef = fabs(obj_coef2[i]);
			}
			if(!obj1_contains_only_ints && !obj2_contains_only_ints) break;
		}
		if(obj1_contains_only_ints || obj2_contains_only_ints) there_will_only_be_points = 1;
		if(obj1_contains_only_ints)
		{
			integer_objective = 1;
			smallest_coef = smallest_obj1_coef;
			printf("Objective 1 only depends on integer variables. The smallest nonzero abs val of a coefficient is: %lf\n",smallest_coef);	
	/*		exit(0);	*/
		}
		else if(obj2_contains_only_ints)
		{
			integer_objective = 2;
			smallest_coef = smallest_obj2_coef;
			printf("Objective 2 only depends on integer variables. The smallest nonzero abs val of a coefficient is: %lf\n",smallest_coef);	
	/*		exit(0);*/
		}
	}
	
	provide_xctype(env, xctype,cur_numcols);
	
/*	status = CPXwriteprob (env, lp1, "myprob2.lp", "LP");*/
/*	exit(0);*/
	
	if(multiplying_factor > 10000.) integer_bb = 0;
	
	if(!integer_bb) there_will_only_be_points = 0;
    	
    	if(integer_bb)
    	{
	    	int did_it = 0;
	    	if(obj1_contains_only_ints)
	    	{
		    	for(i=1;i<total_num_integer;i++)
		    	{
		    		if(!did_it) 
		    		{
		    			smallest_coef = gcd(multiplying_factor*obj_coef1[integer_indices[i-1]],multiplying_factor*obj_coef1[integer_indices[i]]);
		    			if(smallest_coef == 0.) smallest_coef = 2.*3.*5.*7.*11.*13.*17.*1000.;
		    			did_it = 1;
		    		}
		    		else if( fabs(multiplying_factor*obj_coef1[integer_indices[i]]) > .0000001) 
		    			smallest_coef = gcd(smallest_coef,multiplying_factor*obj_coef1[integer_indices[i]]);
		    		if(smallest_coef <= 1.) break;
		    	}
		    	smallest_coef = smallest_coef/((double)multiplying_factor);
	/*	    	printf("smallest coef: %lf, mult fact: %lf\n",smallest_coef,multiplying_factor);*/
		    	printf("All solutions in the objective space will lie on x-values which are multiples of %.9f\n",smallest_coef);
	    	}
	    	else if(obj2_contains_only_ints)
	    	{
		    	for(i=1;i<total_num_integer;i++)
		    	{
		    		if(!did_it) 
		    		{
		    			smallest_coef = gcd(obj_coef2[integer_indices[i-1]],obj_coef2[integer_indices[i]]);
		    			if(smallest_coef == 0.) smallest_coef = 2.*3.*5.*7.*11.*13.;
		    			did_it = 1;
		    		}
		    		else if( fabs(multiplying_factor*obj_coef2[integer_indices[i]]) > .0000001) 
		    			smallest_coef = gcd(smallest_coef,obj_coef2[integer_indices[i]]);
		    		if(smallest_coef <= 1.) break;
		    	}
		    	smallest_coef = smallest_coef/multiplying_factor;
		    	printf("All solutions in the objective space will lie on y-values which are multiples of %lf\n",smallest_coef);
	    	}
    	}
    	
    	/******************************************************************/
    	
    	/********************** Prepare for initial MIP timing test ***********************/
  	
  	stored_x = (double**) malloc(num_x_to_store * sizeof(double*) );
/*  	int i;*/
  	if (stored_x)
	{
	  for(i = 0; i < num_x_to_store; i++)
	  {
	    stored_x[i] = (double*) malloc((cur_numcols+2)* sizeof(double) );
	  }
	}
  	x_rotation = 0;
  	
/*  	pass_the_lps(env, lp1, lp2);*/
  	
  	status = CPXgetlb (env, lp1, stored_x[0], 0, cur_numcols-1);
  	x_rotation++;
  	status = CPXgetub (env, lp1, stored_x[1], 0, cur_numcols-1);
  	x_rotation++;
  	
  	if(pure_binary == 1)
  	{
  		printf("The problem contains only binary variables, no general integers.\n");
  		for(i=0;i<cur_numcols;i++)
	  	{
	  		stored_x[0][i] = 0;
	  		stored_x[1][i] = 1;
	  	}
  	}
  	
  	for(i = x_rotation; i < num_x_to_store; i++)
	{
		int j;
		for(j=0;j<cur_numcols;j++)
		{
	    		stored_x[i][j] = stored_x[0][j];
	    	}
	}
  	
  	obj1_index = cur_numcols;
  	double best_bound = 0.;
  	double obval = 0.;

	status = CPXsetdblparam (env, CPX_PARAM_TILIM, time_limit);
	if ( status ) {
    		printf ("Failed to set solution time limit to 1 seconds, error %d.\n",status);
	   	goto TERMINATE;
	}
  	
  	CPXLPptr pre_lp = CPXcloneprob (env, lp1, &status);
	if ( status ) 
  	{
		fprintf (stderr, "Failed to clone problem 1.\n");
	    	goto TERMINATE;
	}
	
	/******************************************************************/
    	
    	/****************************************************
    		 Solve an initial MIP for 10 seconds 
    		 Whether or not it solves to optimality
    		 determines which type of preprocessing
    		 is used and has other impacts on the
    		 remainder of BB.
    	*****************************************************/

	if(preprocessing_method == 0 && !exact_mips)
	{
		double sec = 10.;
		status = CPXsetdblparam (env, CPX_PARAM_TILIM, sec);
		if ( status ) {
	    		printf ("Failed to set solution time limit to 1 seconds, error %d.\n",status);
		   	goto TERMINATE;
		}
		chg_coefs_(env,pre_lp,indices,-1.0);
	  	CPXmipopt(env, pre_lp);
	  	status2 = CPXgetstat(env, pre_lp);
/*		  	printf("status2: %d\n",status2);*/
	  	if(status2 == 108 || status2 == 105 || status2 == 107 || status2 == 109 || status2 == 111 || status2 == 113) exact_mips = 0;
	  	else exact_mips = 1;
	}
	
	/******************************************************************/
    	
    	/****************************************************
    		 Add two new variables to keep track of
    		 the objective function values.
    	*****************************************************/
	
	double lb[2] = {-CPX_INFBOUND,-CPX_INFBOUND};
	double ub[2] = {CPX_INFBOUND,CPX_INFBOUND};

	double x1_[2] = {0.,0.};
	double x2_[2] = {0.,0.};
	
	status = (CPXnewcols (env, lp1, 2, NULL, lb, ub, NULL, NULL) || CPXnewcols (env, lp2, 2, NULL, lb, ub, NULL, NULL)); 
	if ( status ) 
  	{
    		printf ("CPXnewcols, Failed to add additional variables to the model, error code %d\n", status);
    		goto TERMINATE;
    	}
	
	cur_numcols += 2;
	obj1_index = cur_numcols - 2;
	obj2_index = cur_numcols - 1;
	
	char sense[4] = {'E','E'};
	
	status = (CPXnewrows (env, lp1, 2, NULL, sense, NULL, NULL) || CPXnewrows (env, lp2, 2, NULL, sense, NULL, NULL) ); 
	if ( status ) 
  	{
    		printf ("CPXgetobj, Failed to add additional constraints to the model, error code %d\n", status);
    		goto TERMINATE;
    	}
    	
    	int obj1_row_index = cur_numrows;
    	int obj2_row_index = cur_numrows+1;
    	
    	cur_numrows += 2;
    	
    	int *rowlist;
    	int *collist;
    	double *vallist;
    	
    	rowlist = (int *) malloc (cur_numcols*sizeof(int));
    	collist = (int *) malloc (cur_numcols*sizeof(int));
    	vallist = (double *) malloc (cur_numcols*sizeof(double));
    	
    	/******************************* Set coefficients for first new row *****************/
    	for(i=0;i<cur_numcols;i++)
    	{
    		rowlist[i] = cur_numrows - 2;
    		collist[i] = i;
    		if(i < cur_numcols - 2) vallist[i] = obj_coef1[i];
    		else if(i == cur_numcols - 2) vallist[i] = -1;
    		else vallist[i] = 0;
    	}
    	
    	status = (CPXchgcoeflist (env, lp1, cur_numcols, rowlist, collist, vallist) || CPXchgcoeflist (env, lp2, cur_numcols, rowlist, collist, vallist));
    	if ( status ) 
  	{
    		printf ("CPXgchgcoeflist, Failed to set coefficients for first new row, error code %d\n", status);
    		goto TERMINATE;
    	}
    	
    	double rh1 = -obj1_extra_val;
    	double rh2 = -obj2_extra_val;
    	
    	status = (CPXchgrhs (env, lp1, 1, &rowlist[0], &rh1) || CPXchgrhs (env, lp2, 1, &rowlist[0], &rh1));
    	/*************************************************************************************/
    	/******************************* Set coefficients for second new row *****************/
    	for(i=0;i<cur_numcols;i++)
    	{
    		rowlist[i] = cur_numrows - 1;
    		collist[i] = i;
    		if(i < cur_numcols - 2) vallist[i] = obj_coef2[i];
    		else if(i == cur_numcols - 1) vallist[i] = -1;
    		else vallist[i] = 0;
    	}
    	
    	status = ( CPXchgcoeflist (env, lp1, cur_numcols, rowlist, collist, vallist) || CPXchgcoeflist (env, lp2, cur_numcols, rowlist, collist, vallist));
    	if ( status ) 
  	{
    		printf ("CPXgchgcoeflist, Failed to set coefficients for second new row, error code %d\n", status);
    		goto TERMINATE;
    	}
    	
    	status = (CPXchgrhs (env, lp1, 1, &rowlist[0], &rh2) || CPXchgrhs (env, lp2, 1, &rowlist[0], &rh2));
    	
    	pass_the_lps(env, lp1, lp2);
    	
    	/******************************************************************/
    	
    	/****************************************************
    		 Compute extreme Pareto solns.
    	*****************************************************/
  	
  	compute_feas_extremes   (env_just_solve_mips,env,
  				lp1,
  				lp2,
  				obj_coef1,
  				obj_coef2,
  				indices);

	status = CPXsetdblparam (env, CPX_PARAM_TILIM, time_limit);
	if ( status ) {
    		printf ("Failed to set solution time limit to 1 seconds, error %d.\n",status);
	   	goto TERMINATE;
	}
  				
  	/******************************************************************/
  	
	
	/**************** Set up the problem to be solved prior to use of prepopulate **************************************/
	
	cur_numcols = CPXgetnumcols (env, lp1);
  	cur_numrows = CPXgetnumrows (env, lp1);
  	
  	x_range = SE_extreme_x-NW_extreme_x;
  	y_range = NW_extreme_y-SE_extreme_y;
  	max_range = fmax(x_range,y_range);
  	printf("Objective space f1 range: %lf\nf2 range: %lf\nmax range: %lf ... times 0.0001: %lf\n",x_range,y_range,max_range,0.0001*max_range);
/*  	exit(0);*/
  	
  	double changed_vals[4] = {NW_extreme_x,SE_extreme_y,SE_extreme_x,NW_extreme_y};
  	char changed_lu[4] = {'L','L','U','U'};
  	int changed_indices[4] = {obj1_index,obj2_index,obj1_index,obj2_index};
  	
  	if(obj1_contains_only_ints)
  	{
/*  		printf("this\n");*/
  		integer_objective_lb = NW_extreme_x;
  		changed_vals[0] += smallest_coef;
  		changed_vals[2] -= smallest_coef;
  	}
  	else if(obj2_contains_only_ints)
  	{
/*  		printf("this 2\n");*/
  		integer_objective_lb = SE_extreme_y;
  		changed_vals[1] += smallest_coef;
  		changed_vals[3] -= smallest_coef;
  	}
	
	printf("NWx,NWy,SEx,SEy: %lf,%lf,%lf,%lf\n",NW_extreme_x,NW_extreme_y,SE_extreme_x,SE_extreme_y);
	extreme_x = NW_extreme_x;
	extreme_y = SE_extreme_y;
/*	exit(0);*/
	status = CPXchgbds (env, lp1, 4, changed_indices, changed_lu, changed_vals); 
	if ( status ) 
  	{
    		printf ("CPXchgbds, Failed to change bounds, error code %d\n", status);
    		goto TERMINATE;
    	}
    	status = CPXchgbds (env, lp2, 4, changed_indices, changed_lu, changed_vals); 
	if ( status ) 
  	{
    		printf ("CPXchgbds, Failed to change bounds, error code %d\n", status);
    		goto TERMINATE;
    	}
    	
    	int objective_indices[2] = {obj1_index,obj2_index};
    	
    	status = CPXcopyprotected (env, lp1, 2, objective_indices);
    	status = CPXcopyprotected (env, lp2, 2, objective_indices);
    	
/*    	status = CPXwriteprob (env, lp1, "myprob4.lp", "LP");*/
/*	status = CPXwriteprob (env, lp2, "myprob5.lp", "LP");*/
/*	exit(0);*/

	CPXLPptr pre_lp2 = NULL, presolve2_lp = NULL, presolve2_mip = NULL;

/*	if(there_will_only_be_points) goto START_SOLVING;*/
 
	status = CPXsetintparam (env, CPX_PARAM_REPEATPRESOLVE, 0);
	status = CPXsetintparam (env, CPX_PARAM_PREIND, CPX_OFF);
	
	finish_presolve = clock();
	duration_presolve = (double)(finish_presolve - start_BB) / CLOCKS_PER_SEC;
  	printf("*******************\nTotal time spent until start of prepopulate: %lf\n************************\n",
		  									duration_presolve);
		  									
    	/*************************************************************************************
    	 ____  ____  ____     ____   __  ____  _  _  __     __  ____  ____ 
	(  _ \(  _ \(  __)___(  _ \ /  \(  _ \/ )( \(  )   / _\(_  _)(  __)
	 ) __/ )   / ) _)(___)) __/(  O )) __/) \/ (/ (_/\/    \ )(   ) _) 
	(__)  (__\_)(____)   (__)   \__/(__)  \____/\____/\_/\_/(__) (____)
	
	*************************************************************************************/
	
	double preprocessing_time = 0.;
	
	presolve = preprocessing;
	if(presolve)
	{
		status = CPXsetdblparam (env, CPX_PARAM_TILIM, time_limit);
		heur_limit = times_to_try_heuristic_per_preprocessing_iteration;
		if(preprocessing_method != 0 && search_for_extra_solutions_during_preprocessing == 0) exact_mips = 1;
		else if(preprocessing_method != 0 && search_for_extra_solutions_during_preprocessing == 1) exact_mips = 0;
		if(preprocessing_method == 2 || (preprocessing_method == 0 && exact_mips))
	  	{
			status = epsilon_constraint_preprocessing(yet_another_preprocessing_algorithm, start_presolve, env, lp1, indices, obj_coef2, ub, lb,
								  cur_numcols, lp2, split_pt_denom2, finish_presolve, preprocessing_time, duration_presolve, 
								  bb_results, num_mips_to_solve, numsols, suppress_file_output);
			if(status) goto TERMINATE;
	  	}
	  	else if(preprocessing_method == 1 || preprocessing_method == 3 || (preprocessing_method == 0 && !exact_mips))
	  	{
			if(preprocessing_method == 3 || preprocessing_method == 0) intervals = 1;
			else intervals = 0;
			status = weighted_sum_preprocessing    (start_presolve, lp1, lp2, env, cur_numcols, cur_numrows, indices, obj_coef2, intervals, run_it, 
								split_pt_denom, finish_presolve, preprocessing_time, duration_presolve, bb_results, obj_coef1, suppress_file_output);
			if(status) goto TERMINATE;
	  	}
	  	status = CPXsetdblparam (env, CPX_PARAM_TILIM, max_time);
	}
	
/*	exit(0);*/

	if(stop_PSA_full_early_after_preprocessing) check_for_stopping_PSA_full = 1;
	
	/***************************************************************************************/
	
	/******* Uncomment the following lines to print primal solutions after prepopulate *****/
	
/*	printf("****************************\n Solutions after prepopulate\n**********************************\n");*/
/*	print_inorder(tree,1);*/
/*	printf("****************************\n**********************************\n");*/
	
	/***************************************************************************************/
	
	/*************************************************************************************
     ____  ____  ____      ____   __   __    _  _  ____    ____  _  _   __  
	(  _ \(  _ \(  __)___ / ___) /  \ (  )  / )( \(  __)  (_  _)/ )( \ /  \ 
 	 ) __/ )   / ) _)(___)\___ \(  O )/ (_/\\ \/ / ) _)     )(  \ /\ /(  O )
	(__)  (__\_)(____)    (____/ \__/ \____/ \__/ (____)   (__) (_/\_) \__/ 
	Variable Fixing Based on Pareto Bounds
	
	*************************************************************************************/
	
	printf("variable fixing phase 2:\n");
	double *ub_pre = (double *) malloc (cur_numcols*sizeof(double));
	double *lb_pre = (double *) malloc (cur_numcols*sizeof(double));
	double *x = (double *) malloc (cur_numcols*sizeof(double));
	status = CPXgetlb (env, lp1, lb_pre, 0, cur_numcols-1);
	status = CPXgetub (env, lp1, ub_pre, 0, cur_numcols-1);
	
	int k;
	presolve2_mip = CPXcloneprob (env, lp1, &status);
	presolve2_lp = CPXcloneprob (env, lp1, &status);
	CPXchgprobtype(env, presolve2_lp, CPXPROB_LP);
	
/*	status = CPXwriteprob (env, lp1, "myprob.lp", "LP");*/
/*	printf("-*_*_*_*_*_*_*_*_*_*_*_*-\n");*/
/*	PSA_all(env,presolve2_lp);*/
/*	printf("-*_*_*_*_*_*_*_*_*_*_*_*-\n");*/
	
	int *basis_col = (int *) malloc (cur_numcols*sizeof(int));
	int *basis_row = (int *) malloc (cur_numrows*sizeof(int));
	
	double NW_objectives[2] = {0.,0.};
	double SE_objectives[2] = {0.,0.};
	double WS_objectives[2] = {0.,0.};
	
/*	int insert_check = 0;*/
	int PSA_reduce_val = 0;
	double slope = -1.;
	double temp_x,temp_y;
	int fix_it2 = 0;
	int num_bounds_reduced = 0;
	
	clock_t phase2_start_time, phase2_intermediate_time;
	double time_so_far = 0.;

	phase2_start_time = clock();
	
	if(presolve_phase_2 && presolve)
	{
		status = probing2(phase2_intermediate_time, time_so_far, phase2_start_time, fix_it2, lb_pre, ub_pre, sym, env, presolve2_lp, presolve2_mip, 
				 indices, SE_objectives, insert_check, NW_objectives, WS_objectives, temp_x, slope, lp1, lp2, num_bounds_reduced, cur_numcols, 
				 xctype, cur_numrows, obj_coef1, obj_coef2, weighted_coefs_, temp_y);
		if(status) goto TERMINATE;

		reduce_prob(env,lp1,lp2,indices,1);

	  	cur_numcols = CPXgetnumcols (env, lp1);
	  	cur_numrows = CPXgetnumrows (env, lp1);
	  	
	  	obj1_index = cur_numcols-2;
	  	obj2_index = cur_numcols-1;

		pass_the_lps(env, lp1, lp2);

	  	free_xctype();
	  	status = CPXgetctype (env, lp1, xctype, 0, cur_numcols-1);
	    	if ( status ) 
	  	{
	    		printf ("(%d) CPXgetctype, Failed to get variable types, error code %d\n", __LINE__,status);
	    		return status;
	    	}
	    	
	    	provide_xctype(env, xctype,cur_numcols);
	  	free_coefs();
	    	provide_coefs(obj_coef1,obj_coef2,cur_numcols);
	    	weighted_coefs_[obj1_index] = 0.;
	  	weighted_coefs_[obj2_index] = 0.;

	}
	free(ub_pre);
	free(lb_pre);
	free(basis_col);
	free(basis_row);

	finish_presolve = clock();
	duration_presolve = (double)(finish_presolve - start_BB) / CLOCKS_PER_SEC;
  	printf("*******************\nTotal time spent until end of presolve phase 2: %lf\n************************\n",
		  									duration_presolve);
	if(!suppress_file_output) fprintf(bb_results,"%lf\t",duration_presolve);
	if(!suppress_file_output) fprintf(bb_results,"%d\t",num_bounds_reduced);
	
	/************************************************************************************/
	
	/***********************************************************
	
		Prepare for actually solving the BOMILP
		
	***********************************************************/

	START_SOLVING:
	
	initial_slope = (ub[1]-lb[1])/(lb[0]-ub[0]);
	if(initial_slope > 1000. || initial_slope < .001) initial_slope = -1.;
	
	chg_coefs_(env,lp1,indices,initial_slope);

    	pass_a_mip(lp1);

	status = CPXsetintparam (env, CPX_PARAM_REPEATPRESOLVE, 3);
 	status = CPXsetintparam (env, CPX_PARAM_PREIND, CPX_OFF);  //<-- Turn CPLEX presolve on or off

	if((preprocessing_method == 0 && !exact_mips) || CPLEX_cuts == 2)
	{
	  	printf("performing aggressive cut generation\n");
	  
	/*	printf("setting L and P\n");*/
		  status = CPXsetintparam (env, CPX_PARAM_LANDPCUTS, 	3);
	/*	printf("setting clique\n");*/
	  	  status = CPXsetintparam (env, CPX_PARAM_CLIQUES, 	3);
	/*  	printf("setting covers\n");*/
		  status = CPXsetintparam (env, CPX_PARAM_COVERS, 	3);
	/*	printf("setting disjunctive\n");*/
		  status = CPXsetintparam (env, CPX_PARAM_DISJCUTS, 	3);
	/*	printf("setting frac cuts\n");*/
		  status = CPXsetintparam (env, CPX_PARAM_FRACCUTS, 	2);
	/*	printf("setting gub covers\n");*/
		  status = CPXsetintparam (env, CPX_PARAM_GUBCOVERS, 	2);
	/*	printf("setting implied bounds\n");*/
		  status = CPXsetintparam (env, CPX_PARAM_IMPLBD, 	2);
	/*	printf("setting mir cuts\n");*/
		  status = CPXsetintparam (env, CPX_PARAM_MIRCUTS, 	2);
	/*	printf("setting zero half\n");*/
		  status = CPXsetintparam (env, CPX_PARAM_ZEROHALFCUTS, 2);
	/*	printf("setting flow cover\n");*/
		  status = CPXsetintparam (env, CPX_PARAM_FLOWCOVERS, 	2);
	/*	printf("setting flow path\n");*/
		  status = CPXsetintparam (env, CPX_PARAM_FLOWPATHS, 	2);
	/*	printf("setting mcf\n");*/
		  status = CPXsetintparam (env, CPX_PARAM_MCFCUTS, 2);
	}
	if(!CPLEX_cuts)
	{
		  status = CPXsetintparam (env, CPX_PARAM_LANDPCUTS, 	-1);
	  	  status = CPXsetintparam (env, CPX_PARAM_CLIQUES, 	-1);
		  status = CPXsetintparam (env, CPX_PARAM_COVERS, 	-1);
		  status = CPXsetintparam (env, CPX_PARAM_DISJCUTS, 	-1);
		  status = CPXsetintparam (env, CPX_PARAM_FRACCUTS, 	-1);
		  status = CPXsetintparam (env, CPX_PARAM_GUBCOVERS, 	-1);
		  status = CPXsetintparam (env, CPX_PARAM_IMPLBD, 	-1);
		  status = CPXsetintparam (env, CPX_PARAM_MIRCUTS, 	-1);
		  status = CPXsetintparam (env, CPX_PARAM_ZEROHALFCUTS, -1);
		  status = CPXsetintparam (env, CPX_PARAM_FLOWCOVERS, 	-1);
		  status = CPXsetintparam (env, CPX_PARAM_FLOWPATHS, 	-1);
		  status = CPXsetintparam (env, CPX_PARAM_MCFCUTS, -1);
	}
	
	/************* Optionally add constraints generated as CPLEX cuts **********************/
	/******************* for various weighted sum objectives *******************************/
	
	if(generate_CPLEX_global_cuts_for_several_weighted_sum_objectives)
	{
		CPXLPptr temp_lp = NULL;
		CPXENVptr env_new;
		double lambda = 0.;
		int num_objs_tried = 0, num_added_cuts = 0;
		    
		env_new = CPXopenCPLEX(&status);
		if(status)
		{
			printf("failed to open another copy of CPLEX\n");
			exit(0);
		}
		
		lp1_get_cuts_copy = CPXcloneprob (env_new, lp1, &status);
		
		status = CPXsetintparam (env_new, CPX_PARAM_REPEATPRESOLVE, 0);
	 	status = CPXsetintparam (env_new, CPX_PARAM_PREIND, CPX_OFF);
	  	status = CPXsetintparam (env_new, CPX_PARAM_ADVIND, 0);
	  	status = CPXsetintparam (env_new, CPX_PARAM_REDUCE, 0);
	  	status = CPXsetintparam (env_new, CPX_PARAM_HEURFREQ, -1);
		status = CPXsetintparam ((CPXENVptr) env_new, CPX_PARAM_SCRIND, CPX_OFF);
		status = CPXsetintparam ((CPXENVptr) env_new, CPX_PARAM_THREADS, 1);
		
		status = CPXsetbranchcallbackfunc(env_new, usersetbranch3,  NULL) || CPXsetincumbentcallbackfunc(env_new, userincumbent2, NULL);
/*		status = CPXwriteprob (env_new, lp1_get_cuts_copy, "myprob.lp", "LP");*/
		
		chg_coefs2(env_new, lp1_get_cuts_copy, indices, .0);
/*		printf("about to optimize\n");*/
		status = CPXmipopt (env_new, lp1_get_cuts_copy);
/*		printf("exitting optimization\n");*/
		CPXchgprobtype(env_new, lp1_get_cuts_copy, CPXPROB_MILP);
		status = CPXchgctype(env_new, lp1_get_cuts_copy, cur_numcols, indices, xctype);
		num_objs_tried++;
		
		chg_coefs2(env_new, lp1_get_cuts_copy, indices, 1.);
/*		printf("about to optimize\n");*/
		status = CPXmipopt (env_new, lp1_get_cuts_copy);
/*		printf("exitting optimization\n");*/
		CPXchgprobtype(env_new, lp1_get_cuts_copy, CPXPROB_MILP);
		status = CPXchgctype(env_new, lp1_get_cuts_copy, cur_numcols, indices, xctype);
		num_objs_tried++;
		
/*		status = CPXwriteprob (env_new, lp1_get_cuts_copy, "myprob1.lp", "LP");*/
		
		for(i=0;i<num_iterations_weights;i++)
		{
/*			printf("iteration: %d\n",i+1);*/
			int divisor = pow(2,i+1);
			double size = 1./((double) divisor);
			for(j=0;j<divisor-1;j++)
			{
				lambda = 0. + (j+1)*size;
/*				printf("lambda: %lf\n",lambda);*/
				chg_coefs2(env_new, lp1_get_cuts_copy, indices, lambda);
				status = CPXmipopt (env_new, lp1_get_cuts_copy);
				CPXchgprobtype(env_new, lp1_get_cuts_copy, CPXPROB_MILP);
				status = CPXchgctype(env_new, lp1_get_cuts_copy, cur_numcols, indices, xctype);
				j++;
				num_objs_tried++;
			}
		}
		int new_numrows = CPXgetnumrows (env_new, lp1_get_cuts_copy);
		if(new_numrows > cur_numrows)
		{
			int nzcnt_new = 0, surplus_new = 0;
			num_added_cuts = new_numrows - cur_numrows;
			int *rmatbeg_new = (int *) malloc ((num_added_cuts+1) * sizeof (int));
			int *rmatind_new = (int *) malloc (num_added_cuts*cur_numcols * sizeof (int));
			double *rmatval_new = (double *) malloc (num_added_cuts*cur_numcols * sizeof (double));
			double *rhs_new = (double *) malloc (num_added_cuts * sizeof (double));
			char *sense_new = (char *) malloc (num_added_cuts * sizeof (char));
			status = CPXgetrows (env_new, lp1_get_cuts_copy, &nzcnt_new, rmatbeg_new, rmatind_new, rmatval_new, num_added_cuts*cur_numcols,
				 	     &surplus_new, cur_numrows, new_numrows-1);
			if(status) printf("(%d) there was an error\n",__LINE__);
			status = CPXgetrhs  (env_new, lp1_get_cuts_copy, rhs_new, cur_numrows, new_numrows-1);
			status = CPXgetsense (env_new, lp1_get_cuts_copy, sense_new, cur_numrows, new_numrows-1);
			status = CPXaddrows (env, lp1, 0, num_added_cuts, nzcnt_new, rhs_new, sense_new, rmatbeg_new, rmatind_new, rmatval_new, NULL, NULL);
			status = CPXaddrows (env, lp2, 0, num_added_cuts, nzcnt_new, rhs_new, sense_new, rmatbeg_new, rmatind_new, rmatval_new, NULL, NULL);
			free(rmatbeg_new);
			free(rmatind_new);
			free(rmatval_new);
			free(rhs_new);
			free(sense_new);
			cur_numrows = new_numrows;
		}
/*		status = CPXwriteprob (env, lp1, "myprob1.lp", "LP");*/
		printf("*******************\nAdded %d cuts from %d weighted sum objectives!!\n************************\n",num_added_cuts,num_objs_tried);
/*		exit(0);*/
	}  	
	
	/************* Set up user defined MIP callbacks ***************************************/

 	status =
    	CPXsetbranchcallbackfunc     	(env, usersetbranch,  NULL) ||
    	CPXsetincumbentcallbackfunc 	(env, userincumbent,  NULL); 
    	
    	if(show_progress || control_node_selection) CPXsetnodecallbackfunc (env, userselectnode, NULL);

	if(local_cuts)	CPXsetusercutcallbackfunc       (env, cutcallback,    NULL) ;
  	if ( status ) {
    		printf ("CPXsetsolvecallback failed, error code %d.\n", status);
  	}

  	/* Let MIP callbacks work on the original model, otherwise CPLEX works on the presolved model*/

  	status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
 	if ( status )  goto TERMINATE;
 	
 	status = CPXsetintparam (env_just_solve_mips, CPX_PARAM_MIPCBREDLP, CPX_OFF);
 	if ( status )  goto TERMINATE;

 	status = CPXsetintparam (env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL); //<-- other options are AUTO, DYNAMIC
   	if ( status )  goto TERMINATE;
   	
   	status = CPXsetintparam (env, CPX_PARAM_NODESEL, CPX_NODESEL_BESTEST); //<-- other options are AUTO, DYNAMIC
   	if ( status )  goto TERMINATE;
   	
/*   	status = CPXsetintparam (env, CPXPARAM_Emphasis_MIP,CPX_MIPEMPHASIS_FEASIBILITY); //<-- other options are BALANCED, OPTIMALITY, ...*/
/*   	if ( status )  goto TERMINATE;*/
   	
   	status = CPXsetdblparam (env, CPX_PARAM_TILIM, max_time);
   	status = CPXsetdblparam (env_just_solve_mips, CPX_PARAM_TILIM, time_limit);
 	 
 	/***************************************************************************************/
 	
 	/***************************************************************************************/
 	/*  Calculate starting gap, if desired.		                                       */
 	/***************************************************************************************/
 	
 	if(show_progress)
 	{
 		node *tree_copy = NULL, *tree_copy2 = NULL;
		tree_copy = copy_tree(tree,NULL);
		tree2 = copy_tree(tree,NULL);
		
		CPXLPptr temp_lp1 = CPXcloneprob (env, lp1, &status);
		CPXchgprobtype(env, temp_lp1, CPXPROB_LP);
		
		build_dual_bd(env,temp_lp1);
		
/*		printf("*************************************\n");*/
/*		print_inorder(tree_copy,1);*/
/*		printf("*************************************\n");*/
/*		print_inorder(tree2,2);*/
/*		printf("*************************************\n");*/

		if(hausdorff_or_hypervolume == 2)
		{
			double hy_v2 = get_hypervolume(tree2, 1, 0.);
	/*		printf("dual hypervolume: %lf\n",hy_v2);*/
			double hy_v = get_hypervolume(tree_copy, 1, 0.);
	/*		printf("primal hypervolume: %lf\n",hy_v);*/
			double pd = 100.*(hy_v2-hy_v)/hy_v2;
			
			printf("Nodes processed: %d\t\t Gap measure: %lf\n",branch_iterations,pd);//,length_percent_of_max_range);
			
			if(pd < duality_gap_limit) 
			{
				break_early = 1;
				clock_t now_time = clock();
	
				double cumulative_time = (double)(now_time - start_BB) / CLOCKS_PER_SEC;
				printf("Hypervolume gap is less than %lf%% after preprocessing!!\nTime spent until < %lf%% hypervolume gap: %lf seconds\n",duality_gap_limit,duality_gap_limit,cumulative_time);
				exit(0);
			}
		}
		else
		{
	/*		printf("percent difference: %lf\n",pd);*/
	/*		*/
			double hd_dist = calculate_max_proximal_hd_dist(tree_copy, tree2, 0.);
			double new_val = get_nadirs(tree_copy, 1, 0.);
			hd_dist = fmax(hd_dist,new_val);
			double percent_of_max_range = 100.*hd_dist/max_range;	
	/*		*/
			printf("Nodes processed: %d\t\t Gap measure: %lf\n",branch_iterations,percent_of_max_range);//,length_percent_of_max_range);

			if(percent_of_max_range < duality_gap_limit) 
			{
				break_early = 1;
				clock_t now_time = clock();
	
				double cumulative_time = (double)(now_time - start_BB) / CLOCKS_PER_SEC;
				printf("Hausdorff gap is less than %lf%% after preprocessing!!\nTime spent until < %lf%% hypervolume gap: %lf seconds\n",duality_gap_limit,duality_gap_limit,cumulative_time);
				exit(0);
			}
		}
		destroy_tree(tree_copy);
		destroy_tree(tree2);
		tree2 = NULL;
 	}
 	 
 	/***************************************************************************************/
 	/*  Solve the problem!   			                                       */
 	/***************************************************************************************/
 	
 	status = CPXsetintparam ((CPXENVptr) env, CPX_PARAM_SCRIND, CPX_OFF);
 	
 	frac_scores = (double *) malloc ((cur_numcols+3)*sizeof(double));
   	frac_values = (double *) malloc ((cur_numcols+3)*sizeof(double));
 	
	printf("Beginning optimization ...\n");
/*	if(exact_mips) printf("Single objective MIPs will be solved exactly.\n");*/
/*	else printf("Single objective MIPs will NOT be solved exactly.\n");*/
  	
/*  	double total_time_processing_nodes = 0., total_time_branching = 0., total_time_tightening_variable_bounds = 0., total_time_generating_cuts = 0.;*/
/*  	double total_time_solving_mips = 0.;*/
	double prior_time = 0., this_time = 0., longest_time = 0.;
	
/*	status = CPXwriteprob (env, lp1, "myprob2.lp", "LP");*/
/*	exit(0);*/
	
	if(exploit_objective_gaps)
	{
		changed_indices[0] = obj1_index;
	  	changed_indices[1] = obj2_index;
	  	changed_indices[2] = obj1_index;
	  	changed_indices[3] = obj2_index;
	  	
		printf("*****\nAttempting to exploit gaps in objective space in order to reduce complexity of search\n*****\n");
		double remaining_x_lb = NW_extreme_x, remaining_x_ub = SE_extreme_x, remaining_y_lb = SE_extreme_y, remaining_y_ub = NW_extreme_y;
		find_separations(tree, 1);
		if(num_separators)
		{
			printf("Splitting the Objective Space into subregions based on the following points:\n");
			for(i=num_separators-1;i>-1;i--) //num_separators;i++)
		   	{
/*		   		printf("i: %d\n",i);*/
		   		printf("plot([%lf],[%lf],'-go');\n",x_separators[i],y_separators[i]);
		   	}
		   	if(points_only) printf("Also note that after preprocessing, all known solutions are singletons. Solving MIPs at all nodes until segments are found.\n");
		}
		else printf("There were no gaps to exploit\n");
/*	   	exit(0);*/
		int printed_yet = 0;
	   	for(i=num_separators-1;i>-2;i--) //num_separators;i++)
/*	   	for(i=num_separators-2;i>num_separators-3;i--)*/
	   	{
	   		if(show_progress) printf("\nProcessing subproblem %d\n",num_separators-i);
/*	   		printf("i: %d\n",i);*/
	   		if(duration_BB > max_time && !printed_yet) 
	   		{
	   			printf("Time limit has been exceeded, but not all subproblems have been explored.\n");
	   			printf("Calculating dual bound at open subproblems and closing them.\n");
	   			printed_yet = 1;
	   		}
/*	   		if(i==num_separators - 2) do_it = 1;*/
/*	   		printf("plot([%lf],[%lf],'-go');\n",x_separators[i],y_separators[i]);*/
/*			printf("(%d) separator_dir: %d\n",i,separator_dir[i]);*/
/*			exit(0);*/
  			if(i == -1)
  			{
				changed_vals[0] = remaining_x_lb;
  				changed_vals[1] = remaining_y_lb;
  				changed_vals[2] = remaining_x_ub;
  				changed_vals[3] = remaining_y_ub;
  				status = CPXchgbds (env, lp1, 4, changed_indices, changed_lu, changed_vals);
  			}
  			else 
  			{
/*  				printf("checking direction\n");*/
	  			if(separator_dir[i] == 0)
	  			{
/*	  				printf("direction was 0\n");*/
	  				changed_vals[0] = remaining_x_lb;
/*	  				printf("changed first val\n");*/
	  				if(points_only) changed_vals[1] = y_separators[i] - .01;
	  				else changed_vals[1] = y_separators[i];
/*	  				printf("changed second val\n");*/
	  				if(points_only) changed_vals[2] = x_separators[i] + .01;
	  				else changed_vals[2] = x_separators[i];
/*	  				printf("changed third val\n");*/
	  				changed_vals[3] = remaining_y_ub;
/*	  				printf("changed fourth val\n");*/
	  				status = CPXchgbds (env, lp1, 4, changed_indices, changed_lu, changed_vals);
	  				if(status)
	  				{
	  					printf("(%d) Trouble changing bounds, error code: %d\n",__LINE__,status);
	  				}
/*	  				printf("changed bounds\n");*/
	  				remaining_x_lb = x_separators[i] + .01;
/*	  				printf("changed rem x lb\n");*/
	  			}
	  			else
	  			{
/*	  				printf("direction wasn't 0\n");*/
	  				changed_vals[0] = remaining_x_lb;
	  				changed_vals[1] = y_separators[i];
	  				changed_vals[2] = remaining_x_ub;
	  				changed_vals[3] = remaining_y_ub;
	  				status = CPXchgbds (env, lp1, 4, changed_indices, changed_lu, changed_vals);
	  				remaining_x_lb = x_separators[i];
	  				remaining_y_ub = y_separators[i] - .01;
	  			}
  			}
/*  			printf("going to solve submip %d\n",i);*/

			intermediate_time_start = clock();
			if(i == num_separators-1)
			{
				prior_time = (double)(intermediate_time_start - start_BB) / CLOCKS_PER_SEC;
				printf("===========================\nThe time until starting BB: %lf\n----------------------\n",prior_time);
			}
			
			last_cutcallback_seqnum = -1;
  			status = CPXmipopt (env, lp1);
		  	if ( status ) {
		    		printf ("Failed to optimize MIP, error code %d\n", status);
		    		int substatus = CPXgetsubstat (env, lp1);
		  		printf("substat: %d\n",substatus);
		    		goto TERMINATE;
		  	}
		  	
		  	intermediate_time_finish = clock();
		  	this_time = (double)(intermediate_time_finish - intermediate_time_start) / CLOCKS_PER_SEC;
		  	printf("Time to run BB on subregion %d: %lf\n",num_separators-i,this_time);
		  	if(this_time > longest_time) longest_time = this_time;
		  	
		  	if(i == -1) 
		  	{
		  		printf("----------------------\nMaximum of times to run BB on an individual subregion: %lf\n",longest_time);
		  		printf("Hence, total time solving in parallel could be: %lf\n=========================\n",prior_time + longest_time);
			}  	
		  	
/*  			total_time_processing_nodes += time_processing_nodes;*/
/*  			total_time_branching += time_branching;*/
/*  			total_time_tightening_variable_bounds += time_tightening_variable_bounds;*/
/*  			total_time_generating_cuts += time_generating_cuts;*/
/*  			total_time_solving_mips += time_solving_mips;*/
	   	}
/*		time_processing_nodes = total_time_processing_nodes;*/
/*		time_branching = total_time_branching;*/
/*		time_tightening_variable_bounds = total_time_tightening_variable_bounds;*/
/*		time_generating_cuts = total_time_generating_cuts;*/
/*		time_solving_mips = total_time_solving_mips;*/
	}
	else
	{
		points_only = 0;
		status = CPXmipopt (env, lp1);
	  	if ( status ) {
	    		printf ("Failed to optimize MIP, error code %d\n", status);
	    		int substatus = CPXgetsubstat (env, lp1);
	  		printf("substat: %d\n",substatus);
	    		goto TERMINATE;
	  	}
  	}

	finish_BB = clock();
	 
	int solstat = CPXgetstat (env, lp1);
  	printf ("Solution status %d.\n", solstat); // Note that due to the nature of our callbacks, solution status will likely always be 103.

	int terminated_early = 0;
	if(duration_BB > max_time || branch_iterations > max_nodes || (show_progress && break_early)) terminated_early = 1;
   
  	duration_BB = (double)(finish_BB - start_BB) / CLOCKS_PER_SEC;
  	
  	
  	/**** Print results of interest ************************************/
  	
  	printf("*******************\nTotal time: %lf\n************************\n",duration_BB);
  	printf("\nStructure time: %lf\t Which is %lf percent of BB time.\n************************\n",struct_time, struct_time/duration_BB*100.);
  	if(!suppress_file_output) fprintf(bb_results,"%lf\t",duration_BB);
  	printf("Number of BB nodes explored: %d\n", branch_iterations);
  	if(!suppress_file_output) fprintf(bb_results,"%d\t",branch_iterations);
  	printf("Number of BB nodes at which local cuts were discovered: %d\n",found_local_cuts);
  	if(!suppress_file_output) fprintf(bb_results,"%d\t",found_local_cuts);
  	printf("Number of BB nodes at which branching occured on an objective space disjunction: %d\n",number_pareto_branches);
  	if(!suppress_file_output) fprintf(bb_results,"%d\t",number_pareto_branches);
  	printf("Number of BB nodes at which a disjuntive cut was generated: %d\n",number_disj_cuts);
  	if(!suppress_file_output) fprintf(bb_results,"%d\t",number_disj_cuts);
  	printf("Number of BB nodes at which mips were solved: %d\n",num_nodes_with_mips_solved);
  	if(!suppress_file_output) fprintf(bb_results,"%d\t",num_nodes_with_mips_solved);
  	printf("Number of BB nodes fathomed due to PSA completion: %d\n", fathomed_by_PSA_completion);
  	if(!suppress_file_output) fprintf(bb_results,"%d\t",fathomed_by_PSA_completion);
  	printf("Number of BB nodes fathomed due to dominated local ideal points: %d\n", fathomed_by_dominated_local_ideal_pts);
  	if(!suppress_file_output) fprintf(bb_results,"%d\t",fathomed_by_dominated_local_ideal_pts);
  	printf("Number of BB nodes fathomed due to dominated local ideal segments: %d\n", fathomed_by_dominated_local_ideal_segments);
  	if(!suppress_file_output) fprintf(bb_results,"%d\t",fathomed_by_dominated_local_ideal_segments);
  	printf("Number of BB nodes fathomed due to one dominated local ideal point and one dominated local ideal segment: %d\n", 
  											fathomed_by_1_dominated_pt_1_dominated_segment);
  	if(!suppress_file_output) fprintf(bb_results,"%d\t",fathomed_by_1_dominated_pt_1_dominated_segment);
  	printf("Number of BB nodes fathomed due to dominated lower bound set: %d\n", fathomed_by_dominated_lb);
  	if(!suppress_file_output) fprintf(bb_results,"%d\t",fathomed_by_dominated_lb);
  	printf("\n");
  	printf("Total number of times a bound was reduced during branching: %d\n", bds_reduced);
  	if(!suppress_file_output) fprintf(bb_results,"%d\t",bds_reduced);
  	printf("Number of times a bound was reduced during branching due to infeasibility: %d\n", bds_reduced_by_infeasibility);
  	if(!suppress_file_output) fprintf(bb_results,"%d\t",bds_reduced_by_infeasibility);
  	printf("Number of times a bound was reduced during branching due to dominated local ideal point: %d\n", bds_reduced_by_dominated_ideal_pt);
  	if(!suppress_file_output) fprintf(bb_results,"%d\t",bds_reduced_by_dominated_ideal_pt);
  	printf("Number of times a bound was reduced during branching due to dominated local ideal segment: %d\n", bds_reduced_by_dominated_ideal_segment);
  	if(!suppress_file_output) fprintf(bb_results,"%d\t",bds_reduced_by_dominated_ideal_segment);
  	printf("Number of times a bound was reduced during branching due to dominated lower bound set: %d\n", bds_reduced_by_dominated_lb);
  	if(!suppress_file_output) fprintf(bb_results,"%d\t",bds_reduced_by_dominated_lb);
  	
  	printf("\n\nSome Timing Results (Note that these are not mutually exclusive, for example, mips are solved DURING node processing)\n");
  	printf("Total time processing nodes: %lf ... Which makes for an average of %lf per node\n", time_processing_nodes, 
  		time_processing_nodes/branch_iterations);
  	if(!suppress_file_output) fprintf(bb_results,"%lf\t",time_processing_nodes);
  	printf("Total time branching: %lf ... Which makes for an average of %lf per node\n", time_branching, 
  		time_branching/branch_iterations);
  	if(!suppress_file_output) fprintf(bb_results,"%lf\t",time_branching);
  	printf("Total time tightening variable bounds: %lf ... Which makes for an average of %lf per node\n", time_tightening_variable_bounds, 
  		time_tightening_variable_bounds/branch_iterations);
  	if(!suppress_file_output) fprintf(bb_results,"%lf\t",time_tightening_variable_bounds);
  	printf("Total time generating local cuts: %lf ... Which makes for an average of %lf per node at which cuts were found\n", time_generating_cuts, 
  		time_generating_cuts/found_local_cuts);
  	if(!suppress_file_output) fprintf(bb_results,"%lf\t",time_generating_cuts);
  	printf("Total time generating disjunctive cuts: %lf ... Which makes for an average of %lf per node at which cuts were found\n",
  		time_generating_disjunction_cuts, time_generating_disjunction_cuts/number_disj_cuts);
  	if(!suppress_file_output) fprintf(bb_results,"%lf\t",time_generating_cuts);
  	printf("Total time solving MIPs: %lf\n", time_solving_mips);
  	if(!suppress_file_output) fprintf(bb_results,"%lf\t",time_solving_mips);
  	
  	printf("\n\n_=_=_=_=_=_=_=_=\nMax Time to solve a MIP: %lf\n_=_=_=_=_=_=_=_=\n\n",max_time_to_solve_a_mip);
  	
  	if(!suppress_file_output) fprintf(bb_results,"%lf\t",longest_time);
  	if(!suppress_file_output) fprintf(bb_results,"%lf\t",prior_time+longest_time);
  	
  	if(!suppress_file_output && write_tree_results) fprintf(output2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\n",duration_BB,insert_time,insert_time/duration_BB*100.,struct_time, struct_time/duration_BB*100.,insert_time2,insert_time2/duration_BB*100.,struct_time2, struct_time2/duration_BB*100.,branch_iterations,get_num_inserts(),get_num_nodes(tree),get_tree_depth(tree));
  	
  	/**********************Code for TERMINATE*************************/
  	
  	TERMINATE:

  	/* Free the problem as allocated by CPXcreateprob and
     	CPXreadcopyprob, if necessary *************************************/

  	if ( lp1 != NULL || lp2 != NULL ) 
  	{
    		status = CPXfreeprob (env, &lp1) || CPXfreeprob (env, &lp2);
    		if ( status ) 
    		{
      			printf ("CPXfreeprob failed, error code %d.\n",status);
    		}
  	}
  	/******************************************************************/

  	/* Free the CPLEX environment, if necessary ***********************/
	free_probs();
	close_all_env();
/*  	if(env!=NULL) status=CPXcloseCPLEX(&env);*/
  	/******************************************************************/

	/******************* Print Solutions **********************/
	double hd_dist = 0., percent_of_max_range = 0., length = 0., length_percent_of_max_range = 0.;
	if(terminated_early && tree2)
	{
		if(!suppress_file_output) fprintf(bb_results,"0\t");
		printf("*******************************************\n\n\n");
/*		printf("Attempting to calculate Hausdorff distance between primal and dual bounds\n");*/
		
/*		double hd_dist = 0., percent_of_max_range = 0.;*/
/*		hd_dist = calculate_hd_dist(tree, tree2, 0.);*/
/*		percent_of_max_range = 100.*hd_dist/max_range;*/
/*		*/
/*		printf("Hausdorff distance between primal and dual bounds is %lf which is %lf percent of the max range %lf\n",hd_dist,*/
/*			percent_of_max_range,max_range);*/
			
/*		printf("Calculating maximum proximal Hausdorff distance between each element of primal bound and each element of dual bound (including local nadir pts) ...\n");*/

		if(hausdorff_or_hypervolume == 2)
		{
			double hy_v2 = get_hypervolume(tree2, 1, 0.);
/*			printf("dual hypervolume: %lf\n",hy_v2);*/
			double hy_v = get_hypervolume(tree, 1, 0.);
/*			printf("primal hypervolume: %lf\n",hy_v);*/
			double pd = 100.*(hy_v2-hy_v)/hy_v2;
/*			printf("percent difference: %lf\n",pd);*/
	/*		exit(0);*/
		
			length = get_length(tree2,0.);	
			length_percent_of_max_range = 100.*length/max_range;
		
			printf("Hypervolume of primal bound: %lf\t Hypervolume of dual bound: %lf\t Gap: %lf  (%lf as a percent difference)\n\n",hy_v,
				hy_v2,hy_v2-hy_v,pd);
			if(!suppress_file_output) fprintf(bb_results,"%lf\t%lf\t",hy_v2-hy_v,pd);
			
			printf("Overall length of dual bound is %lf which is %lf percent of the max range %lf\n",length,
				length_percent_of_max_range,max_range);
			if(!suppress_file_output) fprintf(bb_results,"%lf\t%lf",length,length_percent_of_max_range);
	
			if(print_Pareto) 
			{
				if(print_Pareto == 1) printf("\n*********************\n Global Dual bound: \n**********************\n\n\n");
				else fprintf(pareto_results, "\n*********************\n Global Dual bound: \n**********************\n\n\n");
				print_inorder(tree2,2,print_Pareto-1,pareto_results,matlab);
			}
			destroy_tree(tree2);
		}
		else
		{
			hd_dist = calculate_max_proximal_hd_dist(tree, tree2, 0.);
			double new_val = get_nadirs(tree, 1, 0.);
			hd_dist = fmax(hd_dist,new_val);
			percent_of_max_range = 100.*hd_dist/max_range;	
		
			length = get_length(tree2,0.);	
			length_percent_of_max_range = 100.*length/max_range;
		
			printf("Maximum proximal Hausdorff distance is %lf which is %lf percent of the max range %lf\n\n",hd_dist,
				percent_of_max_range,max_range);
			if(!suppress_file_output) fprintf(bb_results,"%lf\t%lf\t",hd_dist,percent_of_max_range);
			
			printf("Overall length of dual bound is %lf which is %lf percent of the max range %lf\n",length,
				length_percent_of_max_range,max_range);
			if(!suppress_file_output) fprintf(bb_results,"%lf\t%lf",length,length_percent_of_max_range);
	
			if(print_Pareto) 
			{
				if(print_Pareto == 1) printf("\n*********************\n Global Dual bound: \n**********************\n\n\n");
				else fprintf(pareto_results, "\n*********************\n Global Dual bound: \n**********************\n\n\n");
				print_inorder(tree2,2,print_Pareto-1,pareto_results,matlab);
			}
			destroy_tree(tree2);
		}
		
	}
	else if(!suppress_file_output) fprintf(bb_results,"1\t0.\t0.\t0.\t0.");
	
	if(!terminated_early) 
	{
		if(print_Pareto == 1) printf("*******************\n Terminated Normally! Primal Solutions: \n************************\n");
		else printf("*******************\n Terminated Normally! \n************************\n");
	}
  	else 
  	{
  		if(print_Pareto == 1) printf("*******************\n Terminated Early! Primal Solutions: \n************************\n");
  		else printf("*******************\n Terminated Early! \n************************\n");
  	}
  	if(print_Pareto == 2) fprintf(pareto_results, "*******************\n Primal Solutions: \n************************\n");
  	
  /*  	print_preorder (tree, output);*/
	prev_node = NULL;
/*	clean_it(tree,0);*/

/*	printf("Total Number of Inserts: %d\n",get_num_inserts());*/
/*	printf("Total Number of Mock Inserts: %d\n",get_num_mock_inserts());*/
/*	printf("Total Number of Nodes: %d\n",get_num_nodes(tree));*/
/*	printf("Final Depth of Tree: %d\n\n\n*******************\n",get_tree_depth(tree));*/
	
	if(print_Pareto) print_inorder(tree,1,print_Pareto-1,pareto_results,matlab);
	
	destroy_tree (tree);
	tree = NULL;
	
	if(terminated_early)
	{
		printf("Maximum proximal Hausdorff distance is %lf which is %lf percent of the max range %lf\n\n",hd_dist,
				percent_of_max_range,max_range);
			
		printf("Overall length of dual bound is %lf which is %lf percent of the max range %lf\n",length,
				length_percent_of_max_range,max_range);
	}
  	/******************************************************************/

	/***************** Close Files ************************************/

  	if (bb_results && !suppress_file_output) 
  	{
  		fprintf(bb_results,"\n");
  		fclose(bb_results);
  	}
  	if(!suppress_file_output && write_tree_results) fclose(output2);
/*  	if (inserted_data) fclose(inserted_data);*/
/*  	if (init_nadir) fclose(init_nadir);*/
  	/*  fclose(all_inserted);*/
  	/*  fclose(testing_struct);*/
/*  	if (init_sol) fclose(init_sol);*/
  	/*  close_files();*/
  	/******************************************************************/
  	
  	/******************* Free Variables *******************************/

  	if(obj_coef1) free_and_null ((char **) &obj_coef1);
	if(obj_coef2) free_and_null ((char **) &obj_coef2);
	if(weighted_coefs_) free_and_null ((char **) &weighted_coefs_);
 	if(indices) free_and_null ((char **) &indices);
	if(rowlist) free_and_null ((char **) &rowlist);
	if(collist) free_and_null ((char **) &collist);
	if(vallist) free_and_null ((char **) &vallist);
  	if(xctype) free_and_null ((char **) &xctype);
  	if(x) free_and_null ((char **) &x);
  	if(first_x) free_and_null ((char **) &first_x);
  	if(frac_scores) free_and_null ((char **) &frac_scores);
  	if(frac_values) free_and_null ((char **) &frac_values);
  	if(row_sense) free_and_null ((char **) &row_sense);
  	if(global_beg) free_and_null ((char **) &global_beg);
  	if(global_varindices) free_and_null ((char **) &global_varindices);
  	if(global_values) free_and_null ((char **) &global_values);
  	if(global_effortlevel) free_and_null ((char **) &global_effortlevel);
  	if(cut_prob_bd_indices_disj1) free_and_null ((char **) &cut_prob_bd_indices_disj1);
  	if(cut_prob_bd_indices_disj2) free_and_null ((char **) &cut_prob_bd_indices_disj2);
  	if(cut_prob_pi_indices) free_and_null ((char **) &cut_prob_pi_indices);
  	if(x_separators) free_and_null ((char **) &x_separators);
  	if(y_separators) free_and_null ((char **) &y_separators);
  	if(separator_dir) free_and_null ((char **) &separator_dir);
  	if(temp_x_r) free_and_null ((char **) &temp_x_r);
  	if(temp_x_l) free_and_null ((char **) &temp_x_l);
  	if(userhandle_current)
	{
/* 		printf("checking whether or not x solns are stored\n");*/
 		if(userhandle_current->x_ws) free_and_null((char **) &userhandle_current->x_ws);
/* 		printf("freeing current x1\n");*/
 		if(userhandle_current->x1) free_and_null((char **) &userhandle_current->x1);
/* 		printf("x2?\n");*/
 		if(userhandle_current->x2) free_and_null((char **) &userhandle_current->x2);
 		free(userhandle_current);
 		userhandle_current = NULL;
	}
	if(userhandle_up)
 	{
 		if(userhandle_up->x_ws) free_and_null((char **) &userhandle_up->x_ws);
/* 		printf("freeing up x1\n");*/
 		if(userhandle_up->x1) free_and_null((char **) &userhandle_up->x1);
 		if(userhandle_up->x2) free_and_null((char **) &userhandle_up->x2);
 		free_and_null((char **) &userhandle_up);
 		userhandle_up = NULL;
 	}
 	if(userhandle_down)
 	{
 		if(userhandle_down->x_ws) free_and_null((char **) &userhandle_down->x_ws);
/* 		printf("freeing down x1\n");*/
 		if(userhandle_down->x1) free_and_null((char **) &userhandle_down->x1);
 		if(userhandle_down->x2) free_and_null((char **) &userhandle_down->x2);
 		free_and_null((char **) &userhandle_down);
 		userhandle_down = NULL;
 	}
	if(x_ws) free_and_null ((char **) &x_ws);
	if(x1) free_and_null ((char **) &x1);
	if(x2) free_and_null ((char **) &x2);
  	free_coefs();
  	free_xctype();
  	if(stored_x)
  	{
	  	for(i=0;i<num_x_to_store;i++)
	  	{
	  		if(stored_x[i]) free(stored_x[i]);
	  	}
	  	free(stored_x);
  	}
  	if(pre_lp) CPXfreeprob(env, &pre_lp);
  	if(pre_lp2) CPXfreeprob(env, &pre_lp2);
  	if(lp1) CPXfreeprob(env, &lp1);
  	if(lp2) CPXfreeprob(env, &lp2);
  	if(presolve2_lp) CPXfreeprob(env, &presolve2_lp);
  	if(presolve2_mip) CPXfreeprob(env, &presolve2_mip);
  	if(global_mip) CPXfreeprob(env_just_solve_mips, &global_mip);
  	if(cut_problem) CPXfreeprob(env, &cut_problem);
  	/******************************************************************/

  	return 0;
  	
  	/************** End of TERMINATE **********************************/
}

/**************************************************************/
/** free and nulll                                           **/
/**************************************************************/

static void free_and_null (char **ptr)
{
  if ( *ptr != NULL ) {
    free (*ptr);
    *ptr = NULL;
  }
} /* END free_and_null */
