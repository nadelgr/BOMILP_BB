#include<stdlib.h>
#include<stdio.h>
#include <math.h>
#include<time.h>
#include "cplex.h"
#include "heuristics.h"

/*********************************************************************************************************************** 

	This function is the genetic algorithm used for generating integer feasible solutions during preprocessing
	when only binary variables are involved in the problem.
	
***********************************************************************************************************************/

int binary_heuristic(CPXCENVptr env, CPXLPptr lp1, CPXLPptr lp2, const int *indices)
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
	
/*	int *indices = NULL;*/
/*	*/
/*	indices = (int *) malloc (numcols*sizeof(int));*/
/*	for(i=0;i<numcols;i++) indices[i] = i;*/
	
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
/*    	if(indices)  free_and_null ((char **) &indices);*/
    	CPXfreeprob(env,&temp_lp);
    	
    	return num_added;
}

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

static void free_and_null (char **ptr)
{
  if ( *ptr != NULL ) {
    free (*ptr);
    *ptr = NULL;
  }
}
