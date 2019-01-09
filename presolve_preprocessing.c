#include "cplex.h"
#include "presolve_preprocessing.h"
#include "bb-bicriteria.h"
#include "max_tree.h"
#include "heuristics.h"

int run_presolve_phase1	(int cur_numcols, FILE *bb_results, int num_fixed_phase1, int num_singleton_columns1, int num_singleton_columns2, 
				double *obj_coef1, double *obj_coef2, CPXCLPptr redlp1, CPXCLPptr redlp2, CPXENVptr env, int *cmatbeg, int *row_indices, 
				double *column_entries, int cur_numrows, CPXLPptr lp1, char *sense2, int fix_it, char *xctype, int *singleton_column_indices2,
				int *their_rows1, int *their_rows2, char *sym, double l, double u, CPXLPptr lp2, int *singleton_column_indices1, int *col_indices,
				double *row_entries, double *lower_bound, double *upper_bound, int suppress_file_output)
{
	clock_t phase1_st, phase1_fi;
	double ti = 0.;
	int i=0, status=0, nzcnt=0, surplus=0, j=0, num_singletons_fixed = 0, num_dom_cols_tightened = 0, num_dom_cols_fixed = 0, cur_size = 10;
	phase1_st = clock();
	if(duality_fixing)
	{
		for(i=0;i<cur_numcols;i++)
		{
			phase1_fi = clock();
			ti = (double)(phase1_fi - phase1_st) / CLOCKS_PER_SEC;
			if(ti > max_phase1_time) 
			{
				if(!suppress_file_output) 
				{
				    fprintf(bb_results,"%d\t",num_fixed_phase1);
				    fprintf(bb_results,"%d\t",num_singleton_columns1+num_singleton_columns2);
				    fprintf(bb_results,"%d\t",num_singletons_fixed);
				    fprintf(bb_results,"%d\t",num_dom_cols);
				    fprintf(bb_results,"%d\t",num_dom_cols_tightened);
				    fprintf(bb_results,"%d\t",num_dom_cols_fixed);
				}
			
				printf("Phase 1 exceeded time limit of %lf. Moving on.\n",max_phase1_time);
				goto END_OF_PHASE1;
			}
	/*		printf("considering x%d\n",i);*/
			if(obj_coef1[i] <= 0. && obj_coef2[i] <= 0.)
			{
	/*			printf("x%d is a candidate for fixing\n",i);*/
				if(redlp1)
				{
					status = CPXgetcols (env, redlp1, &nzcnt, cmatbeg, row_indices, column_entries, cur_numrows, &surplus, i, i);
					if ( status ) 
		  			{
		    				printf ("Unable to get column entries for column %d, error code %d\n", i,status);
		    				return status;
		    			}
	    			}
	    			else
	    			{
	    				status = CPXgetcols (env, lp1, &nzcnt, cmatbeg, row_indices, column_entries, cur_numrows, &surplus, i, i);
					if ( status ) 
		  			{
		    				printf ("Unable to get column entries for column %d, error code %d\n", i,status);
		    				return status;
		    			}
	    			}
	    			for(j=0;j<nzcnt;j++)
	    			{
	/*    				printf("\t%dth col entry: %lf\n",j,column_entries[j]);*/
	    				if((sense2[row_indices[j]] == 'L' && column_entries[j] < 0.) || 
	    					(sense2[row_indices[j]] == 'G' && column_entries[j] > 0.) ||
	    				   	sense2[row_indices[j]] == 'E') break;
	    				if(j == nzcnt-1) fix_it = 1;
	    			}
	    			if(cur_numrows == 0) fix_it = 1;
	    			if(nzcnt == 1 && !fix_it && xctype[i] == 'C' && column_entries[0] <= 0.)
	    			{
	    				singleton_column_indices2[num_singleton_columns2] = i;
	    				their_rows2[num_singleton_columns2] = row_indices[0];
	    				num_singleton_columns2++;	
	    			}
	    			if(fix_it == 1)
	    			{
	/*    				printf("fixing x%d to its lb\n",i);*/
	    				sym[0] = 'U';
	    				if(redlp1) status = CPXgetlb (env, redlp1, &l, i, i);
	    				else status = CPXgetlb (env, lp1, &l, i, i);
	/*    				printf("the lower bound is: %lf\n",l);*/
					if(redlp1)
					{
		    				status = CPXchgbds (env, (CPXLPptr) redlp1, 1, &i, sym, &l);
		    				status = CPXchgbds (env, (CPXLPptr) redlp2, 1, &i, sym, &l);
	    				}
	    				else
					{
		    				status = CPXchgbds (env, (CPXLPptr) lp1, 1, &i, sym, &l);
		    				status = CPXchgbds (env, (CPXLPptr) lp2, 1, &i, sym, &l);
	    				}
	    				fix_it = 0;
	    				num_fixed_phase1++;
	    			}
			}
			if(obj_coef1[i] >= 0. && obj_coef2[i] >= 0.)
			{
	/*			printf("x%d is a candidate for fixing\n",i);*/
				if(redlp1)
				{
					status = CPXgetcols (env, redlp1, &nzcnt, cmatbeg, row_indices, column_entries, cur_numrows, &surplus, i, i);
					if ( status ) 
		  			{
		    				printf ("Unable to get column entries for column %d, error code %d\n", i,status);
		    				return status;
		    			}
	    			}
	    			else
	    			{
	    				status = CPXgetcols (env, lp1, &nzcnt, cmatbeg, row_indices, column_entries, cur_numrows, &surplus, i, i);
					if ( status ) 
		  			{
		    				printf ("Unable to get column entries for column %d, error code %d\n", i,status);
		    				return status;
		    			}
	    			}
	    			for(j=0;j<nzcnt;j++)
	    			{
	    				if((sense2[row_indices[j]] == 'L' && column_entries[j] > 0.) || 
	    					(sense2[row_indices[j]] == 'G' && column_entries[j] < 0.) || sense2[row_indices[j]] == 'E') break;
	    				if(j == nzcnt-1) fix_it = 1;
	    			}
	    			if(cur_numrows == 0) fix_it = 1;
	    			if(nzcnt == 1 && !fix_it && xctype[i] == 'C' && column_entries[0] >= 0.)
	    			{
	    				singleton_column_indices1[num_singleton_columns1] = i;
	    				their_rows1[num_singleton_columns1] = row_indices[0];
	    				num_singleton_columns1++;	
	    			}
	    			if(fix_it == 1)
	    			{
	/*    				printf("fixing x%d to its ub\n",i);*/
	    				sym[0] = 'L';
	    				if(redlp1) status = CPXgetub (env, redlp1, &u, i, i);
	    				else status = CPXgetub (env, lp1, &u, i, i);
	/*    				printf("the upper bound is: %lf\n",u);*/
	    				if(redlp1) 
	    				{
		    				status = CPXchgbds (env, (CPXLPptr) redlp1, 1, &i, sym, &u);
		    				status = CPXchgbds (env, (CPXLPptr) redlp2, 1, &i, sym, &u);
	    				}
	    				else
	    				{
		    				status = CPXchgbds (env, (CPXLPptr) lp1, 1, &i, sym, &u);
		    				status = CPXchgbds (env, (CPXLPptr) lp2, 1, &i, sym, &u);
	    				}
	    				fix_it = 0;
	    				num_fixed_phase1++;
	    			}
			}
		}
		printf("variables fixed due to duality fixing: %d\n",num_fixed_phase1);
	}
	if(!suppress_file_output) fprintf(bb_results,"%d\t",num_fixed_phase1);
	
	/************************************************
		Attempting to do singleton stuffing
	************************************************/
	
	printf("Number of singleton columns: %d\n",num_singleton_columns1+num_singleton_columns2);
	if(!suppress_file_output) fprintf(bb_results,"%d\t",num_singleton_columns1+num_singleton_columns2);
	if(singleton_columns && num_singleton_columns1)
	{		
/*			printf("Their indices:\n");*/
/*			for(i=0;i<num_singleton_columns1;i++) printf("%d found in row %d\n",singleton_column_indices1[i],their_rows1[i]);*/
/*			status = CPXwriteprob (env, redlp1, "myprob.lp", "LP");*/
		
		for(i=0;i<num_singleton_columns1;i++)
		{
			phase1_fi = clock();
			ti = (double)(phase1_fi - phase1_st) / CLOCKS_PER_SEC;
			if(ti > max_phase1_time) 
			{
			    if(!suppress_file_output) 
			    {
				    fprintf(bb_results,"%d\t",num_singletons_fixed);
				    fprintf(bb_results,"%d\t",num_dom_cols);
				    fprintf(bb_results,"%d\t",num_dom_cols_tightened);
				    fprintf(bb_results,"%d\t",num_dom_cols_fixed);
				}
			
				printf("Phase 1 exceeded time limit of %lf. Moving on.\n",max_phase1_time);
				goto END_OF_PHASE1;
			}
			if(their_rows1[i] >= 0)
			{
/*					printf("examining row %d\n",their_rows1[i]);*/
				int k = 0, num_singletons_in_this_row = 0, current_row = their_rows1[i];
				double u_j = 0., l_j = 0., current_entry = 0.;
				if(redlp1) 
				{
					status = CPXgetrows (env, redlp1, &nzcnt, cmatbeg, col_indices, row_entries, cur_numcols, &surplus, their_rows1[i],
							their_rows1[i]);
				}
				else
				{
					status = CPXgetrows (env, lp1, &nzcnt, cmatbeg, col_indices, row_entries, cur_numcols, &surplus, their_rows1[i],
							their_rows1[i]);
				}
				their_rows1[i] -= (cur_numrows+1);
				for(j=0;j<nzcnt;j++)
				{
/*					printf("looking at column %d. The entry: %lf, lb: %lf, ub: %lf\n",col_indices[j],row_entries[j],*/
/*						lower_bound[col_indices[j]],upper_bound[col_indices[j]]);*/
					while(singleton_column_indices1[k] < col_indices[j]) k++;
					if(col_indices[j] == singleton_column_indices1[i]) current_entry = row_entries[j];
					if(singleton_column_indices1[k] == col_indices[j]) 
					{
						their_rows1[k] -= (cur_numrows+1);
						u_j += row_entries[j]*lower_bound[col_indices[j]];
/*						l_j += row_entries[j]*lower_bound[col_indices[j]];*/
						num_singletons_in_this_row++;
					}
					else if(row_entries[j] > 0.) 
					{
						u_j += row_entries[j]*upper_bound[col_indices[j]];
/*						l_j += row_entries[j]*lower_bound[col_indices[j]];*/
					}
					else if(row_entries[j] < 0.) 
					{
						u_j += row_entries[j]*lower_bound[col_indices[j]];
/*						l_j += row_entries[j]*upper_bound[col_indices[j]];*/
					}
/*					printf("for this row u_: %lf\t l_: %lf\n",u_j,l_j);*/
				}
				if(num_singletons_in_this_row == 1)
				{
					double rhs[1] = {0.};
					double delta = current_entry*(upper_bound[singleton_column_indices1[i]]-lower_bound[singleton_column_indices1[i]]);
					if(redlp1) status = CPXgetrhs (env, redlp1, rhs, current_row, current_row);
					else status = CPXgetrhs (env, lp1, rhs, current_row, current_row);
				
/*						printf("for this row:\t u_%d: %lf\t l_%d: %lf\t rhs: %lf\t delta: %lf\t coef: %lf\n",*/
/*							current_row,u_j,current_row,l_j, rhs[0],delta, current_entry);*/
				
					if(delta <= rhs[0] - u_j)
					{
						sym[0] = 'L';
		    				if(redlp1) 
		    				{
			    				status = CPXchgbds (env, (CPXLPptr) redlp1, 1, &singleton_column_indices1[i], sym, 
			    					&upper_bound[singleton_column_indices1[i]]);
			    				status = CPXchgbds (env, (CPXLPptr) redlp2, 1, &singleton_column_indices1[i], sym, 
			    					&upper_bound[singleton_column_indices1[i]]);
		    				}
		    				else
		    				{
			    				status = CPXchgbds (env, (CPXLPptr) lp1, 1, &singleton_column_indices1[i], sym, 
			    					&upper_bound[singleton_column_indices1[i]]);
			    				status = CPXchgbds (env, (CPXLPptr) lp2, 1, &singleton_column_indices1[i], sym, 
			    					&upper_bound[singleton_column_indices1[i]]);
	    					}
	    					num_singletons_fixed++;
					}
				}
				else if(num_singletons_in_this_row > 1)
				{
/*						printf("row %d contains more than 1 singleton, figure out how to deal with it\n",current_row);*/
				
/*						printf("defining arrays\n");*/
					struct store_it *array1 = malloc( num_singletons_in_this_row*sizeof( struct store_it ) );
					struct store_it *array2 = malloc( num_singletons_in_this_row*sizeof( struct store_it ) );
				
/*						printf("first loop\n");*/
					int h = 0;
					for(j=0;j<num_singleton_columns1;j++)
					{
/*							printf("j: %d\n",j);*/
						if( (their_rows1[j] >= 0 && their_rows1[j] == current_row) || 
							(their_rows1[j] < 0 && their_rows1[j] + cur_numrows + 1 == current_row) )
						{
							double val = 0.;
						
/*							array1[h] = malloc( sizeof (struct store_it));*/
/*							array2[h] = malloc( sizeof (struct store_it));*/
							if(obj_coef1[singleton_column_indices1[j]] != 0. || obj_coef2[singleton_column_indices1[j]] != 0.)
							{
								if(redlp1) status = CPXgetcoef (env, redlp1, current_row, singleton_column_indices1[j],
									&val);
								else status = CPXgetcoef (env, lp1, current_row, singleton_column_indices1[j], &val);
							}
							if(obj_coef1[singleton_column_indices1[j]] != 0.) array1[h].ratio = 
								obj_coef1[singleton_column_indices1[j]]/val;
							else array1[h].ratio = 0;
							if(obj_coef2[singleton_column_indices1[j]] != 0.) array2[h].ratio = 
								obj_coef2[singleton_column_indices1[j]]/val;
							else array2[h].ratio = 0;
							array1[h].index = singleton_column_indices1[j];
							array2[h].index = singleton_column_indices1[j];
							h++;
						}
					} 
/*						printf("sorting\n");*/
					qsort(array1, h, sizeof(struct store_it), comparison_function4);
					qsort(array2, h, sizeof(struct store_it), comparison_function4);
/*						printf("done sorting\n");*/
				
					double rhs[1] = {0.};
					if(redlp1) status = CPXgetrhs (env, redlp1, rhs, current_row, current_row);
					else status = CPXgetrhs (env, lp1, rhs, current_row, current_row);
				
/*						printf("second loop\n");*/
					for(j=0;j<h;j++)
					{
/*							printf("j: %d\n",j);*/
						if(array1[j].index != array2[j].index) break;
					
						double val = 0.;
/*							printf("getting coef\n");*/
						if(redlp1) status = CPXgetcoef (env, redlp1, current_row, array1[j].index, &val);
						else status = CPXgetcoef (env, lp1, current_row, array1[j].index, &val);
/*							printf("making delta\n");*/
						double delta = val*(upper_bound[array1[j].index]-lower_bound[array1[j].index]);
				
/*							printf("for this row:\t u_%d: %lf\t l_%d: %lf\t rhs: %lf\t delta: %lf\t coef: %lf\n",*/
/*								current_row,u_j,current_row,l_j, rhs[0],delta, val);*/
				
						if(delta <= rhs[0] - u_j)
						{
							sym[0] = 'L';
			    				if(redlp1) 
			    				{
				    				status = CPXchgbds (env, (CPXLPptr) redlp1, 1, &array1[j].index, sym, 
				    					&upper_bound[array1[j].index]);
				    				status = CPXchgbds (env, (CPXLPptr) redlp2, 1, &array1[j].index, sym, 
				    					&upper_bound[array1[j].index]);
			    				}
			    				else
			    				{
				    				status = CPXchgbds (env, (CPXLPptr) lp1, 1, &array1[j].index, sym, 
				    					&upper_bound[array1[j].index]);
				    				status = CPXchgbds (env, (CPXLPptr) lp2, 1, &array1[j].index, sym, 
				    					&upper_bound[array1[j].index]);
		    					}
		    					num_singletons_fixed++;
		    					u_j -= delta;
						}
					}
					free(array1);
					free(array2);
				}
			}
		}
/*			printf("Number of singleton columns which were fixed: %d\n",num_singletons_fixed);*/
	}
	if(singleton_columns && num_singleton_columns2)
	{		
/*			printf("Their indices:\n");*/
/*			for(i=0;i<num_singleton_columns2;i++) printf("%d found in row %d\n",singleton_column_indices2[i],their_rows2[i]);*/
/*			status = CPXwriteprob (env, redlp1, "myprob.lp", "LP");*/
		
		for(i=0;i<num_singleton_columns2;i++)
		{
			phase1_fi = clock();
			ti = (double)(phase1_fi - phase1_st) / CLOCKS_PER_SEC;
			if(ti > max_phase1_time) 
			{
			    if(!suppress_file_output) 
			    {
				    fprintf(bb_results,"%d\t",num_singletons_fixed);
				    fprintf(bb_results,"%d\t",num_dom_cols);
				    fprintf(bb_results,"%d\t",num_dom_cols_tightened);
				    fprintf(bb_results,"%d\t",num_dom_cols_fixed);
				}
			
				printf("Phase 1 exceeded time limit of %lf. Moving on.\n",max_phase1_time);
				goto END_OF_PHASE1;
			}
			if(their_rows2[i] >= 0)
			{
/*					printf("examining row %d\n",their_rows2[i]);*/
				int k = 0, num_singletons_in_this_row = 0, current_row = their_rows2[i];
				double u_j = 0., l_j = 0., current_entry = 0.;
				if(redlp1) 
				{
					status = CPXgetrows (env, redlp1, &nzcnt, cmatbeg, col_indices, row_entries, cur_numcols, &surplus, their_rows2[i],
							their_rows2[i]);
				}
				else
				{
					status = CPXgetrows (env, lp1, &nzcnt, cmatbeg, col_indices, row_entries, cur_numcols, &surplus, their_rows2[i],
							their_rows2[i]);
				}
				their_rows2[i] -= (cur_numrows+1);
				for(j=0;j<nzcnt;j++)
				{
/*					printf("looking at column %d. The entry: %lf, lb: %lf, ub: %lf\n",col_indices[j],row_entries[j],*/
/*						lower_bound[col_indices[j]],upper_bound[col_indices[j]]);*/
					while(singleton_column_indices2[k] < col_indices[j]) k++;
					if(col_indices[j] == singleton_column_indices2[i]) current_entry = row_entries[j];
					if(singleton_column_indices2[k] == col_indices[j]) 
					{
						their_rows2[k] -= (cur_numrows+1);
/*						u_j += row_entries[j]*lower_bound[col_indices[j]];*/
						l_j += row_entries[j]*lower_bound[col_indices[j]];
						num_singletons_in_this_row++;
					}
					else if(row_entries[j] > 0.) 
					{
/*						u_j += row_entries[j]*upper_bound[col_indices[j]];*/
						l_j += row_entries[j]*lower_bound[col_indices[j]];
					}
					else if(row_entries[j] < 0.) 
					{
/*						u_j += row_entries[j]*lower_bound[col_indices[j]];*/
						l_j += row_entries[j]*upper_bound[col_indices[j]];
					}
/*					printf("for this row u_: %lf\t l_: %lf\n",u_j,l_j);*/
				}
				if(num_singletons_in_this_row == 1)
				{
					double rhs[1] = {0.};
					double delta = current_entry*(upper_bound[singleton_column_indices2[i]]-lower_bound[singleton_column_indices2[i]]);
					if(redlp1) status = CPXgetrhs (env, redlp1, rhs, current_row, current_row);
					else status = CPXgetrhs (env, lp1, rhs, current_row, current_row);
				
/*						printf("for this row:\t u_%d: %lf\t l_%d: %lf\t rhs: %lf\t delta: %lf\t coef: %lf\n",*/
/*							current_row,u_j,current_row,l_j, rhs[0],delta, current_entry);*/
				
					if(delta >= rhs[0] - l_j)
					{
						sym[0] = 'L';
		    				if(redlp1) 
		    				{
			    				status = CPXchgbds (env, (CPXLPptr) redlp1, 1, &singleton_column_indices2[i], sym, 
			    					&upper_bound[singleton_column_indices2[i]]);
			    				status = CPXchgbds (env, (CPXLPptr) redlp2, 1, &singleton_column_indices2[i], sym, 
			    					&upper_bound[singleton_column_indices2[i]]);
		    				}
		    				else
		    				{
			    				status = CPXchgbds (env, (CPXLPptr) lp1, 1, &singleton_column_indices2[i], sym, 
			    					&upper_bound[singleton_column_indices2[i]]);
			    				status = CPXchgbds (env, (CPXLPptr) lp2, 1, &singleton_column_indices2[i], sym, 
			    					&upper_bound[singleton_column_indices2[i]]);
	    					}
	    					num_singletons_fixed++;
					}
				}
				else if(num_singletons_in_this_row > 1)
				{
/*						printf("row %d contains more than 1 singleton, figure out how to deal with it\n",current_row);*/
/*					exit(0);*/
/*					*/ 

/*						printf("defining arrays\n");*/
					struct store_it *array1 = malloc( num_singletons_in_this_row*sizeof( struct store_it ) );
					struct store_it *array2 = malloc( num_singletons_in_this_row*sizeof( struct store_it ) );
				
/*						printf("first loop\n");*/
					int h = 0;
					for(j=0;j<num_singleton_columns2;j++)
					{
/*							printf("j: %d\n",j);*/
						if( (their_rows2[j] >= 0 && their_rows2[j] == current_row) || 
							(their_rows2[j] < 0 && their_rows2[j] + cur_numrows + 1 == current_row) )
						{
							double val = 0.;
						
/*							array1[h] = malloc( sizeof (struct store_it));*/
/*							array2[h] = malloc( sizeof (struct store_it));*/
							if(obj_coef1[singleton_column_indices2[j]] != 0. || obj_coef2[singleton_column_indices2[j]] != 0.)
							{
								if(redlp1) status = CPXgetcoef (env, redlp1, current_row, singleton_column_indices2[j],
									&val);
								else status = CPXgetcoef (env, lp1, current_row, singleton_column_indices2[j], &val);
							}
							if(obj_coef1[singleton_column_indices2[j]] != 0.) array1[h].ratio = 
								obj_coef1[singleton_column_indices2[j]]/val;
							else array1[h].ratio = 0;
							if(obj_coef2[singleton_column_indices2[j]] != 0.) array2[h].ratio = 
								obj_coef2[singleton_column_indices2[j]]/val;
							else array2[h].ratio = 0;
							array1[h].index = singleton_column_indices1[j];
							array2[h].index = singleton_column_indices1[j];
							h++;
						}
					} 
/*						printf("sorting\n");*/
					qsort(array1, h, sizeof(struct store_it), comparison_function4);
					qsort(array2, h, sizeof(struct store_it), comparison_function4);
/*						printf("done sorting\n");*/
				
					double rhs[1] = {0.};
					if(redlp1) status = CPXgetrhs (env, redlp1, rhs, current_row, current_row);
					else status = CPXgetrhs (env, lp1, rhs, current_row, current_row);
				
/*						printf("second loop\n");*/
					for(j=0;j<h;j++)
					{
/*							printf("j: %d\n",j);*/
						if(array1[j].index != array2[j].index) break;
					
						double val = 0.;
/*							printf("getting coef\n");*/
						if(redlp1) status = CPXgetcoef (env, redlp1, current_row, array1[j].index, &val);
						else status = CPXgetcoef (env, lp1, current_row, array1[j].index, &val);
/*							printf("making delta\n");*/
						double delta = val*(upper_bound[array1[j].index]-lower_bound[array1[j].index]);
				
/*							printf("for this row:\t u_%d: %lf\t l_%d: %lf\t rhs: %lf\t delta: %lf\t coef: %lf\n",*/
/*								current_row,u_j,current_row,l_j, rhs[0],delta, val);*/
				
						if(delta >= rhs[0] - l_j)
						{
							sym[0] = 'L';
			    				if(redlp1) 
			    				{
				    				status = CPXchgbds (env, (CPXLPptr) redlp1, 1, &array1[j].index, sym, 
				    					&upper_bound[array1[j].index]);
				    				status = CPXchgbds (env, (CPXLPptr) redlp2, 1, &array1[j].index, sym, 
				    					&upper_bound[array1[j].index]);
			    				}
			    				else
			    				{
				    				status = CPXchgbds (env, (CPXLPptr) lp1, 1, &array1[j].index, sym, 
				    					&upper_bound[array1[j].index]);
				    				status = CPXchgbds (env, (CPXLPptr) lp2, 1, &array1[j].index, sym, 
				    					&upper_bound[array1[j].index]);
		    					}
		    					num_singletons_fixed++;
		    					l_j -= delta;
						}
					}
					free(array1);
					free(array2);
				}
			}
		}
/*			printf("Number of singleton columns which were fixed: %d\n",num_singletons_fixed);*/
	}
	printf("Number of singleton columns which were fixed: %d\n",num_singletons_fixed);
	if(!suppress_file_output) fprintf(bb_results,"%d\t",num_singletons_fixed);
/*		exit(0);*/
	
	/************************************************
	   Attempting to check for dominating columns
	************************************************/
	
	num_dom_cols = 0;
	if(dominating_columns) //generate_disjunctive_cuts_from_dominated_columns) 
	{
/*			CPXdelnames (env, lp1);*/
/*			CPXdelnames (env, lp2);*/
/*			status = CPXwriteprob (env, lp1, "myprob3.lp", "LP");*/
/*			status = CPXwriteprob (env, lp2, "myprob4.lp", "LP");*/
		dominated_indices = (int *) malloc (cur_size*sizeof(int));
		dominating_indices = (int *) malloc (cur_size*sizeof(int));
		for(i=0;i<cur_numcols;i++)
		{
			phase1_fi = clock();
			ti = (double)(phase1_fi - phase1_st) / CLOCKS_PER_SEC;
			if(ti > max_phase1_time) 
			{
			    if(!suppress_file_output) 
			    {
				    fprintf(bb_results,"%d\t",num_dom_cols);
				    fprintf(bb_results,"%d\t",num_dom_cols_tightened);
				    fprintf(bb_results,"%d\t",num_dom_cols_fixed);
				}
			
				printf("Phase 1 exceeded time limit of %lf. Moving on.\n",max_phase1_time);
				goto END_OF_PHASE1;
			}
			if(redlp1)
			{
				status = CPXgetcols (env, redlp1, &nzcnt, cmatbeg, row_indices, column_entries, cur_numrows, &surplus, i, i);
				if ( status ) 
	  			{
	    				printf ("Unable to get column entries for column %d, error code %d\n", i,status);
	    				return status;
	    			}
    			}
    			else
    			{
    				status = CPXgetcols (env, lp1, &nzcnt, cmatbeg, row_indices, column_entries, cur_numrows, &surplus, i, i);
				if ( status ) 
	  			{
	    				printf ("Unable to get column entries for column %d, error code %d\n", i,status);
	    				return status;
	    			}
    			}
    			int nzcnt2 = 0;
    			int *row_indices2 = (int *) malloc (cur_numrows*sizeof(double));
			double *column_entries2 = (double *) calloc (cur_numrows, sizeof(double));
			for(j=0;j<cur_numcols;j++)
			{
				double pr = 0;
/*					if(j == 1545 && i == 2) */
/*					{*/
/*						int surplus = 0;*/
/*						char *cur_colname[100];*/
/*						char cur_colnamestore[100];*/
/*						status = CPXgetcolname (env, lp1, cur_colname, cur_colnamestore,100, &surplus, j,j);*/
/*						printf("surplus: %d, status: %d\n",surplus, status);*/
/*						printf("name of dominating column: %s\n",cur_colname[0]);*/
/*						status = CPXgetcolname (env, lp1, cur_colname, cur_colnamestore,100, &surplus, i,i);*/
/*						printf("surplus: %d, status: %d\n",surplus, status);*/
/*						printf("name of dominated column: %s\n",cur_colname[0]);*/
/*						pr = 1;*/
/*					}*/
/*					printf("checking to see if column %d dominates column %d\n",j,i);*/
				double lb_i = 0., lb_j = 0., ub_i = 0., ub_j = 0.;
				if(redlp1)
				{
					status = CPXgetlb (env, redlp1, &lb_i, i, i);
					status = CPXgetlb (env, redlp1, &lb_j, j, j);
					status = CPXgetub (env, redlp1, &ub_i, i, i);
					status = CPXgetub (env, redlp1, &ub_j, j, j);
				}
				else
				{
					status = CPXgetlb (env, lp1, &lb_i, i, i);
					status = CPXgetlb (env, lp1, &lb_j, j, j);
					status = CPXgetub (env, lp1, &ub_i, i, i);
					status = CPXgetub (env, lp1, &ub_j, j, j);
				}
			
				if( i != j && lb_i > -CPX_INFBOUND && ub_i < CPX_INFBOUND && lb_j > -CPX_INFBOUND && ub_j < CPX_INFBOUND && 
					(lb_i != ub_i || lb_j != ub_j) && ( obj_coef1[j] >= obj_coef1[i] && obj_coef2[j] >= obj_coef2[i] ) &&
					( (xctype[i] == 'C' && xctype[j] == xctype[i] ) || (xctype[i] != 'C'  && xctype[j] != 'C') ) )
				{
					
/*						if(pr) printf("col %d obj coefs: %lf, %lf, col %d obj coef: %lf, %lf\n",j,obj_coef1[j], obj_coef2[j],i,obj_coef1[i], obj_coef2[i] );*/
					int h = 0, k = 0, dominates = 1;
					if(redlp1)
					{
						status = CPXgetcols (env, redlp1, &nzcnt2, cmatbeg, row_indices2, column_entries2, cur_numrows, &surplus, j, j);
						if ( status ) 
			  			{
			    				printf ("Unable to get column entries for column %d, error code %d\n", i,status);
			    				return status;
			    			}
		    			}
		    			else
		    			{
		    				status = CPXgetcols (env, lp1, &nzcnt2, cmatbeg, row_indices2, column_entries2, cur_numrows, &surplus, j, j);
						if ( status ) 
			  			{
			    				printf ("Unable to get column entries for column %d, error code %d\n", i,status);
			    				return status;
			    			}
		    			}
/*			    			printf("nzcnt: %d, nzcnt2: %d\n",nzcnt,nzcnt2);*/
					while(h < nzcnt || k < nzcnt2)
					{
						int h_ = min(h,nzcnt-1);
						int k_ = min(k,nzcnt2-1);
/*							printf("for i: (%d) row index: %d, column entry: %lf\n",h,row_indices[h],column_entries[h]);*/
/*							printf("for j: (%d) row index: %d, column entry: %lf\n",k,row_indices2[k],column_entries2[k]);*/
/*							printf("h_,k_: %d,%d\n",h_,k_);*/
						if((row_indices[h_] < row_indices2[k_] || k >= nzcnt2) && h < nzcnt)
						{
/*								printf("checking 1\n");*/
							if(sense2[row_indices[h_]] == 'E' || column_entries[h_] < 0)
							{
								dominates = 0;
								break;
							}
							else h++;
						}
						else if((row_indices[h_] > row_indices2[k_] || h >= nzcnt) && k < nzcnt2)
						{
/*								printf("checking 2\n");*/
							if(sense2[row_indices[k_]] == 'E' || column_entries2[k_] > 0)
							{
								dominates = 0;
								break;
							}
							else k++;
						}
						else if(row_indices[h_] == row_indices2[k_] && 
							((sense2[row_indices[h_]] == 'E' && column_entries[h_] != column_entries2[k_]) || 
							column_entries[h_] < column_entries2[k_]))
						{
/*								printf("cant work\n");*/
							dominates = 0;
							break;
						}
						else
						{
/*								printf("increasing\n");*/
							h++;
							k++;
						}
					}
					if(dominates) 
					{
						num_dom_cols++;
/*							printf("Column %d dominates column %d!\n",j,i);*/
/*							printf("Column %d dominates column %d!\n",j+1,i+1);*/
/*							int surplus = 0;*/
/*							char *cur_colname[100];*/
/*							char cur_colnamestore[100];*/
/*							status = CPXgetcolname (env, lp1, cur_colname, cur_colnamestore,100, &surplus, j,j);*/
/*							printf("surplus: %d, status: %d\n",surplus, status);*/
/*							printf("name of dominating column: %s\n",cur_colname[0]);*/
/*							int r = 0;*/
/*							printf("stuff for x%d column, %s\n",i+1,cur_colname[0]);*/
/*							for(r=0;r<nzcnt;r++)*/
/*							{*/
/*								printf("row: %d, coef: %lf\n",row_indices[r]+1,column_entries[r]);*/
/*							}*/
/*							status = CPXgetcolname (env, lp1, cur_colname, cur_colnamestore,100, &surplus, i,i);*/
/*							printf("surplus: %d, status: %d\n",surplus, status);*/
/*							printf("name of dominated column: %s\n",cur_colname[0]);*/
/*							printf("stuff for x%d column, %s\n",j+1,cur_colname[0]);*/
/*							for(r=0;r<nzcnt2;r++)*/
/*							{*/
/*								printf("row: %d, coef: %lf\n",row_indices2[r]+1,column_entries2[r]);*/
/*							}*/
/*							pr = 1;*/
						

						if(num_dom_cols >= cur_size)
						{
							cur_size = cur_size*2;
							dominated_indices = (int *) realloc(dominated_indices, cur_size*sizeof(int));
							dominating_indices = (int *) realloc(dominating_indices, cur_size*sizeof(int));
						}
						dominated_indices[num_dom_cols-1] = i;
						dominating_indices[num_dom_cols-1] = j;
					
						int m = 0, nzcnt3 = 0, already_fixed = 0;
						double max_val = -100000., min_val = 100000.;
						if(lb_j != ub_j)
						{
/*							printf("Attempting to fix x%d to u%d via maxl\n",j,j);*/
							for(m=0;m<cur_numrows;m++)
							{
								double rhs = 0., a_m_i = 0., a_m_j = 0., val = 0., l_m_i = 0.;
								if(redlp1)
								{
									status = CPXgetcoef (env, redlp1, m, i, &a_m_i);
									status = CPXgetcoef (env, redlp1, m, j, &a_m_j);
								}
								else
								{
									status = CPXgetcoef (env, lp1, m, i, &a_m_i);
									status = CPXgetcoef (env, lp1, m, j, &a_m_j);
								}
								if(a_m_i < 0. && a_m_j < 0.)
								{
									if(redlp1) status = CPXgetrhs (env, redlp1, &rhs, m, m);
									else status = CPXgetrhs (env, lp1, &rhs, m, m);
									val = rhs;
									l_m_i = a_m_i*lb_i;
									if(redlp1) 
									{
										status = CPXgetrows (env, redlp1, &nzcnt3, cmatbeg, col_indices, row_entries,
											cur_numcols,&surplus, m, m);
									}
									else
									{
										status = CPXgetrows (env, lp1, &nzcnt3, cmatbeg, col_indices, row_entries,
											cur_numcols,&surplus, m, m);
									}
									int n = 0;
									for(n=0;n<nzcnt3;n++)
									{
										if(col_indices[n] != i)
										{
											if(row_entries[n] > 0.) l_m_i += 
												row_entries[n]*lower_bound[col_indices[n]];
											else l_m_i += row_entries[n]*upper_bound[col_indices[n]];
										}
									}
/*									printf("L^%d_%d: %lf\n",m,i,l_m_i);*/
									val -= l_m_i;
									val += a_m_j*ub_j;
									val = val/a_m_j;
									if(val > max_val) max_val = val;
								}
							}
/*							printf("maxl^i_j: %lf, u_j: %lf\n",max_val,ub_j);*/
							if(max_val >= ub_j)
							{
								sym[0] = 'L';
				    				if(redlp1) 
				    				{
					    				status = CPXchgbds (env, (CPXLPptr) redlp1, 1, &j, sym, &ub_j);
					    				status = CPXchgbds (env, (CPXLPptr) redlp2, 1, &j, sym, &ub_j);
				    				}
				    				else
				    				{
					    				status = CPXchgbds (env, (CPXLPptr) lp1, 1, &j, sym, &ub_j);
					    				status = CPXchgbds (env, (CPXLPptr) lp2, 1, &j, sym, &ub_j);
			    					}
/*			    					printf("fixing it\n");*/
	/*		    					exit(0);*/
								already_fixed = 1;
			    					num_dom_cols_fixed++;
							}
							else if(max_val >= lb_j)
							{
								sym[0] = 'L';
				    				if(redlp1) 
				    				{
					    				status = CPXchgbds (env, (CPXLPptr) redlp1, 1, &j, sym, &max_val);
					    				status = CPXchgbds (env, (CPXLPptr) redlp2, 1, &j, sym, &max_val);
				    				}
				    				else
				    				{
					    				status = CPXchgbds (env, (CPXLPptr) lp1, 1, &j, sym, &max_val);
					    				status = CPXchgbds (env, (CPXLPptr) lp2, 1, &j, sym, &max_val);
			    					}
/*			    					printf("fixing it\n");*/
	/*		    					exit(0);*/
								lb_j = max_val;
			    					num_dom_cols_tightened++;
							}
/*							else printf("couldn't fix it\n");*/
				
/*							printf("Attempting to fix x%d to u%d via minu\n",j,j);*/
							if(!already_fixed && obj_coef1[j] >= 0. && obj_coef2[j] >= 0.)
							{
								max_val = -100000.;
								min_val = 100000.;
								for(m=0;m<cur_numrows;m++)
								{
	/*								printf("examining row %d, cur max val: %lf\n",m,max_val);*/
									double rhs = 0., a_m_i = 0., a_m_j = 0., val = 0., u_m_i = 0.;
									if(redlp1)
									{
										status = CPXgetcoef (env, redlp1, m, i, &a_m_i);
										status = CPXgetcoef (env, redlp1, m, j, &a_m_j);
									}
									else
									{
										status = CPXgetcoef (env, lp1, m, i, &a_m_i);
										status = CPXgetcoef (env, lp1, m, j, &a_m_j);
									}
									if(a_m_i > 0. && a_m_j > 0.)
									{
										if(redlp1) status = CPXgetrhs (env, redlp1, &rhs, m, m);
										else status = CPXgetrhs (env, lp1, &rhs, m, m);
										val = rhs;
										u_m_i = a_m_i*lb_i;
										if(redlp1) 
										{
											status = CPXgetrows (env, redlp1, &nzcnt3, cmatbeg, col_indices,
												row_entries,cur_numcols, &surplus, m, m);
										}
										else
										{
											status = CPXgetrows (env, lp1, &nzcnt3, cmatbeg, col_indices, row_entries,
												cur_numcols, &surplus, m, m);
										}
										int n = 0;
										for(n=0;n<nzcnt3;n++)
										{
											if(col_indices[n] != i)
											{
												if(row_entries[n] < 0.) u_m_i += 
													row_entries[n]*lower_bound[col_indices[n]];
												else u_m_i += row_entries[n]*upper_bound[col_indices[n]];
											}
										}
/*										printf("U^%d_%d: %lf\n",m,i,u_m_i);*/
										val -= u_m_i;
										val += a_m_j*lb_j;
										val = val/a_m_j;
										if(val > max_val) max_val = val;
									}
								}
/*								printf("minu^i_j: %lf, u_j: %lf\n",max_val,ub_j);*/
								if(max_val >= ub_j)
								{
									sym[0] = 'L';
					    				if(redlp1) 
					    				{
						    				status = CPXchgbds (env, (CPXLPptr) redlp1, 1, &j, sym, &ub_j);
						    				status = CPXchgbds (env, (CPXLPptr) redlp2, 1, &j, sym, &ub_j);
					    				}
					    				else
					    				{
						    				status = CPXchgbds (env, (CPXLPptr) lp1, 1, &j, sym, &ub_j);
						    				status = CPXchgbds (env, (CPXLPptr) lp2, 1, &j, sym, &ub_j);
				    					}
/*				    					printf("fixing it\n");*/
	/*			    					exit(0);*/
				    					num_dom_cols_fixed++;
								}
								else if(max_val >= lb_j)
								{
									sym[0] = 'L';
					    				if(redlp1) 
					    				{
						    				status = CPXchgbds (env, (CPXLPptr) redlp1, 1, &j, sym, &max_val);
						    				status = CPXchgbds (env, (CPXLPptr) redlp2, 1, &j, sym, &max_val);
					    				}
					    				else
					    				{
						    				status = CPXchgbds (env, (CPXLPptr) lp1, 1, &j, sym, &max_val);
						    				status = CPXchgbds (env, (CPXLPptr) lp2, 1, &j, sym, &max_val);
				    					}
/*				    					printf("fixing it\n");*/
	/*			    					exit(0);*/
									lb_j = max_val;
				    					num_dom_cols_tightened++;
								}
/*								else printf("couldn't fix it\n");			*/
							}
						}
				
						if(lb_i != ub_i)
						{
/*							printf("Attempting to fix x%d to l%d via minl\n",i,i);*/
							already_fixed = 0;
							max_val = -100000.;
							min_val = 100000.;
							for(m=0;m<cur_numrows;m++)
							{
								double rhs = 0., a_m_i = 0., a_m_j = 0., val = 0., l_m_j = 0.;
								if(redlp1)
								{
									status = CPXgetcoef (env, redlp1, m, i, &a_m_i);
									status = CPXgetcoef (env, redlp1, m, j, &a_m_j);
								}
								else
								{
									status = CPXgetcoef (env, lp1, m, i, &a_m_i);
									status = CPXgetcoef (env, lp1, m, j, &a_m_j);
								}
								if(a_m_i > 0. && a_m_j > 0.)
								{
									if(redlp1) status = CPXgetrhs (env, redlp1, &rhs, m, m);
									else status = CPXgetrhs (env, lp1, &rhs, m, m);
									val = rhs;
									l_m_j = a_m_j*ub_j;
									if(redlp1) 
									{
										status = CPXgetrows (env, redlp1, &nzcnt3, cmatbeg, col_indices, row_entries,
											cur_numcols,&surplus, m, m);
									}
									else
									{
										status = CPXgetrows (env, lp1, &nzcnt3, cmatbeg, col_indices, row_entries,
											cur_numcols,&surplus, m, m);
									}
									int n = 0;
									for(n=0;n<nzcnt3;n++)
									{
										if(col_indices[n] != j)
										{
											if(row_entries[n] > 0.) l_m_j += 
												row_entries[n]*lower_bound[col_indices[n]];
											else l_m_j += row_entries[n]*upper_bound[col_indices[n]];
										}
									}
/*									printf("L^%d_%d: %lf\n",m,j,l_m_j);*/
									val -= l_m_j;
									val += a_m_i*lb_i;
									val = val/a_m_i;
									if(val < min_val) min_val = val;
								}
							}
/*							printf("minl^j_i: %lf, l_i: %lf\n",min_val,lb_i);*/
							if(min_val <= lb_i)
							{
								sym[0] = 'U';
				    				if(redlp1) 
				    				{
					    				status = CPXchgbds (env, (CPXLPptr) redlp1, 1, &i, sym, &lb_i);
					    				status = CPXchgbds (env, (CPXLPptr) redlp2, 1, &i, sym, &lb_i);
				    				}
				    				else
				    				{
					    				status = CPXchgbds (env, (CPXLPptr) lp1, 1, &i, sym, &lb_i);
					    				status = CPXchgbds (env, (CPXLPptr) lp2, 1, &i, sym, &lb_i);
			    					}
/*			    					printf("fixing it\n");*/
	/*		    					exit(0);*/
								already_fixed = 1;
			    					num_dom_cols_fixed++;
							}
							else if(min_val <= ub_i)
							{
								sym[0] = 'U';
				    				if(redlp1) 
				    				{
					    				status = CPXchgbds (env, (CPXLPptr) redlp1, 1, &i, sym, &min_val);
					    				status = CPXchgbds (env, (CPXLPptr) redlp2, 1, &i, sym, &min_val);
				    				}
				    				else
				    				{
					    				status = CPXchgbds (env, (CPXLPptr) lp1, 1, &i, sym, &min_val);
					    				status = CPXchgbds (env, (CPXLPptr) lp2, 1, &i, sym, &min_val);
			    					}
/*			    					printf("fixing it\n");*/
	/*		    					exit(0);*/
								ub_i = min_val;
			    					num_dom_cols_tightened++;
							}
/*							else printf("couldn't fix it\n");*/
				
/*							printf("Attempting to fix x%d to l%d via maxu\n",i,i);*/
							if(!already_fixed && obj_coef1[j] <= 0. && obj_coef2[j] <= 0.)
							{
								max_val = -100000.;
								min_val = 100000.;
								for(m=0;m<cur_numrows;m++)
								{
									double rhs = 0., a_m_i = 0., a_m_j = 0., val = 0., u_m_j = 0.;
									if(redlp1)
									{
										status = CPXgetcoef (env, redlp1, m, i, &a_m_i);
										status = CPXgetcoef (env, redlp1, m, j, &a_m_j);
									}
									else
									{
										status = CPXgetcoef (env, lp1, m, i, &a_m_i);
										status = CPXgetcoef (env, lp1, m, j, &a_m_j);
									}
									if(a_m_i < 0. && a_m_j < 0.)
									{
										if(redlp1) status = CPXgetrhs (env, redlp1, &rhs, m, m);
										else status = CPXgetrhs (env, lp1, &rhs, m, m);
										val = rhs;
										u_m_j = a_m_j*ub_i;
										if(redlp1) 
										{
											status = CPXgetrows (env, redlp1, &nzcnt3, cmatbeg, col_indices,
												row_entries,cur_numcols, &surplus, m, m);
										}
										else
										{
											status = CPXgetrows (env, lp1, &nzcnt3, cmatbeg, col_indices, row_entries,
												cur_numcols, &surplus, m, m);
										}
										int n = 0;
										for(n=0;n<nzcnt3;n++)
										{
											if(col_indices[n] != j)
											{
												if(row_entries[n] < 0.) u_m_j += 
													row_entries[n]*lower_bound[col_indices[n]];
												else u_m_j += row_entries[n]*upper_bound[col_indices[n]];
											}
										}
/*										printf("U^%d_%d: %lf\n",m,j,u_m_j);*/
										val -= u_m_j;
										val += a_m_i*ub_i;
										val = val/a_m_i;
										if(val < min_val) min_val = val;
									}
								}
/*								printf("maxu^j_i: %lf, l_i: %lf\n",min_val,lb_i);*/
								if(min_val <= lb_i)
								{
									sym[0] = 'U';
					    				if(redlp1) 
					    				{
						    				status = CPXchgbds (env, (CPXLPptr) redlp1, 1, &i, sym, &lb_i);
						    				status = CPXchgbds (env, (CPXLPptr) redlp2, 1, &i, sym, &lb_i);
					    				}
					    				else
					    				{
						    				status = CPXchgbds (env, (CPXLPptr) lp1, 1, &i, sym, &lb_i);
						    				status = CPXchgbds (env, (CPXLPptr) lp2, 1, &i, sym, &lb_i);
				    					}
/*					    					printf("fixing it\n");*/
	/*			    					exit(0);*/
				    					num_dom_cols_fixed++;
								}
								else if(min_val <= ub_i)
								{
									sym[0] = 'U';
					    				if(redlp1) 
					    				{
						    				status = CPXchgbds (env, (CPXLPptr) redlp1, 1, &i, sym, &min_val);
						    				status = CPXchgbds (env, (CPXLPptr) redlp2, 1, &i, sym, &min_val);
					    				}
					    				else
					    				{
						    				status = CPXchgbds (env, (CPXLPptr) lp1, 1, &i, sym, &min_val);
						    				status = CPXchgbds (env, (CPXLPptr) lp2, 1, &i, sym, &min_val);
				    					}
				    					ub_i = min_val;
/*				    					printf("fixing it\n");*/
	/*			    					exit(0);*/
				    					num_dom_cols_tightened++;
								}
/*								else printf("couldn't fix it\n");			*/
							}
						}
					}
				}
			}
			free(row_indices2);
			free(column_entries2);
		}
	}
	printf("Number of times one column dominated another: %d\n",num_dom_cols);
	if(!suppress_file_output) fprintf(bb_results,"%d\t",num_dom_cols);
	printf("Number of times a variable bound was tightened due to dominance: %d\n",num_dom_cols_tightened);
	if(!suppress_file_output) fprintf(bb_results,"%d\t",num_dom_cols_tightened);
	printf("Number of variables fixed due to dominated columns: %d\n",num_dom_cols_fixed);
	if(!suppress_file_output) fprintf(bb_results,"%d\t",num_dom_cols_fixed);
	
/*		for(i=0;i<num_dom_cols;i++) printf("col %d dominates col %d\n",dominating_indices[i],dominated_indices[i]);*/
/*		*/
/*		exit(0);*/
	END_OF_PHASE1:
	
	return 0;

}

int epsilon_constraint_preprocessing(int yet_another_preprocessing_algorithm, clock_t start_presolve, CPXCENVptr env, CPXLPptr lp1, int *indices, 
					double *obj_coef2, double *ub, double *lb, int cur_numcols, CPXLPptr lp2, double split_pt_denom2, clock_t finish_presolve,
					double preprocessing_time, double duration_presolve, FILE *bb_results, int num_mips_to_solve, int numsols, int suppress_file_output)
{
	int status = 0, i = 0, j = 0;
	exact_mips = 1;
	
	if(preprocessing_parameter < 0) 
		preprocessing_parameter = ceil( 4.*(.185*cur_numcols)*(.185*cur_numcols)/((.185*cur_numcols)*(.185*cur_numcols)+1.5*cur_numcols+10000.) + 2.);
	
	if(!yet_another_preprocessing_algorithm)
	{
  		start_presolve = clock();
  		int num_boxes = 0;
/*				printf("beginning presolve n_ya\n");*/
		CPXLPptr presolve_lp = CPXcloneprob (env, lp1, &status);
	  	if ( status ) 
	  	{
	    		fprintf (stderr, "Failed to clone problem 1.\n");
	    		return status;
	  	}
	  	CPXLPptr presolve_lp2 = CPXcloneprob (env, presolve_lp, &status);
	  	if ( status ) 
	  	{
	    		fprintf (stderr, "Failed to clone problem 1.\n");
	    		return status;
	  	}
	  	status = CPXchgobj (env, presolve_lp2, obj1_index+2, indices, obj_coef2);
		if ( status ) {
			printf ("Failed to get change obj coef. Error code %d\n", status);
			return status;
		}
	
		double ubs[2] = {0.,0.};
		double lbs[2] = {0.,0.};
	
		status = CPXgetlb (env, presolve_lp, lbs, obj1_index, obj2_index);
		if ( status ) {
			printf ("Failed to get lb's for objectives. Error code %d\n", status);
			return status;
		}
		status = CPXgetub (env, presolve_lp, ubs, obj1_index, obj2_index);
		if ( status ) {
			printf ("Failed to get ub's for objectives. Error code %d\n", status);
			return status;
		}
	
		box *first_box = (struct box*) malloc( sizeof( struct box ) );
		first_box->x_ub = ub[0];
		first_box->x_lb = lb[0];
		first_box->y_ub = ub[1];
		first_box->y_lb = lb[1];
		first_box->next = NULL;
		num_boxes++;
	
	  	char lu[4] = {'L','L','U','U'};
	  	double bds[4] = {0.,0.,0.,0.};
	  	int vars[4] = {obj1_index,obj2_index,obj1_index,obj2_index};
	    	double *x = (double *) malloc (cur_numcols*sizeof(double));
	    	int presolve_iteration = 1;
	    	if(heur_limit) heur_limit = 12;
		int num_split_pts = 1;
		double x_size = (first_box->x_ub - first_box->x_lb);
		double y_size = (first_box->y_ub - first_box->y_lb);
		double bound;
		int split_x_y;
		int status2;
		if(x_size >= y_size)
		{
		 	bound = (first_box->x_ub + first_box->x_lb)/2.;
		 	bds[0] = bound;
		 	bds[1] = first_box->y_lb;
		 	bds[2] = first_box->x_ub;
		 	bds[3] = first_box->y_ub;
		 	status = CPXchgbds (env, presolve_lp2, 4, vars, lu, bds);
		 	CPXmipopt (env, presolve_lp2);
	    		status2 = CPXgetstat(env, presolve_lp2);
		 	split_x_y = 0;
		}
		else 
		{
			bound = (first_box->y_ub + first_box->y_lb)/2.;
			bds[0] = first_box->x_lb;
		 	bds[1] = bound;
		 	bds[2] = first_box->x_ub;
		 	bds[3] = first_box->y_ub;
		 	status = CPXchgbds (env, presolve_lp, 4, vars, lu, bds);
		 	CPXmipopt (env, presolve_lp);
	    		status2 = CPXgetstat(env, presolve_lp);
			split_x_y = 1;
		} 
	    	int num_inserted = num_split_pts;
	    	int mips_solved = 1;
	    	
/*	    	printf("the status of the solve: %d\n",status2);*/
	    	
	    	if(status2 != 101 && status2 != 102 && status2 != 129 && status2 != 130)
	    	{
			if(status2 == 118)
	    		{
	/*    			printf("The first test problem is unbounded. Exiting!\n");*/
	    			return status;
	    		}
/*			    		else if(status2 == 104 || status2 == 105 || status2 == 107 || status2 == 109 || status2 == 111 || status2 == 113 )*/
/*			 		   printf("Warning: CPLEX failed to reach an optimal solution. A suboptimal integer-feasible solution is being used.\n");*/
	    	}
	    	
	    	if(status2 == 103) goto TRY_ANOTHER_PREPROCESSING;

		if(split_x_y)
		{
			status = CPXgetx (env, presolve_lp, x, 0, cur_numcols-1);
		  	if(status)
		  	{
				printf("(%d) Failed to get x-values from CPLEX.\n",__LINE__);
				printf("(%d) the status of the mipsolve: %d\n",__LINE__,status2);
				printf("status: %d\n",status);
		  		return status;
		  	}
	  	}
	  	else
	  	{
	  		status = CPXgetx (env, presolve_lp2, x, 0, cur_numcols-1);
		  	if(status)
		  	{
				printf("(%d) Failed to get x-values from CPLEX.\n",__LINE__);
		  		return status;
		  	}
	  	}
	  	
	      	int insert_check = 0;
	      	insert_check = mock_insert(1,x[obj1_index],x[obj2_index],0,0,0,&tree);
	      	if(insert_check)
	      	{
	      		for(i=0;i<cur_numcols-2;i++)
	      		{
	      			stored_x[x_rotation][i] = x[i];
/*			      			printf("stored val: %lf\n",stored_x[x_rotation][i]);*/
	      		}
	      		x_rotation = (x_rotation + 1) % num_x_to_store;
	      		insert_check = 0;
	      		PSA_full(env,NULL,x,NULL,NULL);
	      	}
	      	
	      	box *next_box = (struct box*) malloc( sizeof( struct box ) );
	      	if(split_x_y)
	      	{
			next_box->x_ub = first_box->x_ub;
			next_box->x_lb = x[obj1_index];
			next_box->y_ub = bound;
			next_box->y_lb = first_box->y_lb;
			first_box->x_ub = x[obj1_index];
			first_box->y_lb = x[obj2_index];
		}
		else
		{
			next_box->x_ub = first_box->x_ub;
			next_box->x_lb = x[obj1_index];
			next_box->y_ub = x[obj2_index];
			next_box->y_lb = first_box->y_lb;
			first_box->x_ub = bound;
			first_box->y_lb = x[obj2_index];
		}
		next_box->next = NULL;
		first_box->next = next_box;
		num_boxes++;
	      	
		int run_it = 0;
		if(pure_binary == 1) run_it += binary_heuristic(env,lp1,lp2,indices);
		else run_it += mixed_heuristic(env);
		int j = 0;
	
		int cur_num_boxes = num_boxes;

		while(presolve_iteration < 3 || (presolve_iteration < 40 && num_inserted >= cur_num_boxes/split_pt_denom2))
		{
			presolve_iteration++;
			
			finish_presolve = clock();
			preprocessing_time = (double)(finish_presolve - start_presolve) / CLOCKS_PER_SEC;
			if(preprocessing_time > max_preprocessing_time) 
			{
				printf("Preprocessing exceeded time limit of %lf. Moving to presolve phase 2.\n",max_preprocessing_time);
				goto END_PREPROCESSING1;
			}
			
			cur_num_boxes = num_boxes;
			num_inserted = 0;
			box *temp_box = first_box;
			for(j=0;j<cur_num_boxes;j++)
			{	
				x_size = (temp_box->x_ub - temp_box->x_lb);
				y_size = (temp_box->y_ub - temp_box->y_lb);
				if(x_size >= y_size)
				{
				 	bound = (temp_box->x_ub + temp_box->x_lb)/2.;
				 	bds[0] = bound;
				 	bds[1] = temp_box->y_lb;
				 	bds[2] = temp_box->x_ub;
				 	bds[3] = temp_box->y_ub;
				 	status = CPXchgbds (env, presolve_lp2, 4, vars, lu, bds);
				 	CPXmipopt (env, presolve_lp2);
			    		status2 = CPXgetstat(env, presolve_lp2);
				 	split_x_y = 0;
				}
				else 
				{
					bound = (temp_box->y_ub + temp_box->y_lb)/2.;
					bds[0] = temp_box->x_lb;
				 	bds[1] = bound;
				 	bds[2] = temp_box->x_ub;
				 	bds[3] = temp_box->y_ub;
				 	status = CPXchgbds (env, presolve_lp, 4, vars, lu, bds);
				 	CPXmipopt (env, presolve_lp);
			    		status2 = CPXgetstat(env, presolve_lp);
					split_x_y = 1;
				} 
				mips_solved++;
			
				if(status2 == 103 || status2 == 119)
				{
					if(split_x_y) temp_box->y_ub = bound;
					else temp_box->x_ub = bound;
				}
				else
				{
					if(split_x_y)
					{
						status = CPXgetx (env, presolve_lp, x, 0, cur_numcols-1);
					  	if(status)
					  	{
							printf("(%d) Failed to get x-values from CPLEX.\n",__LINE__);
					  		return status;
					  	}
				  	}
				  	else
				  	{
				  		status = CPXgetx (env, presolve_lp2, x, 0, cur_numcols-1);
					  	if(status)
					  	{
							printf("(%d) Failed to get x-values from CPLEX.\n",__LINE__);
/*									printf("status2: %d\n",status2);*/
					  		return status;
					  	}
				  	}
				  	
				  	if((split_x_y == 0 && fabs(x[obj1_index] - temp_box->x_ub) < .0001) ||
				  	  (split_x_y == 1 && fabs(x[obj2_index] - temp_box->y_ub) < .0001))
				  	{
				  		if(split_x_y) temp_box->y_ub = bound;
						else temp_box->x_ub = bound;
				  	}
					else
					{
						box *next_box = (struct box*) malloc( sizeof( struct box ) );
					      	if(split_x_y)
					      	{
							next_box->x_ub = temp_box->x_ub;
							next_box->x_lb = x[obj1_index];
							next_box->y_ub = bound;
							next_box->y_lb = temp_box->y_lb;
							temp_box->x_ub = x[obj1_index];
							temp_box->y_lb = x[obj2_index];
						}
						else
						{
							next_box->x_ub = temp_box->x_ub;
							next_box->x_lb = x[obj1_index];
							next_box->y_ub = x[obj2_index];
							next_box->y_lb = temp_box->y_lb;
							temp_box->x_ub = bound;
							temp_box->y_lb = x[obj2_index];
						}
						if(temp_box->next) next_box->next = temp_box->next;
						else next_box->next = NULL;
						temp_box->next = next_box;
						num_boxes++;
						temp_box = temp_box->next;
					}
				}
				temp_box = temp_box->next;
			  	
			      	int insert_check = 0;
			      	insert_check = mock_insert(1,x[obj1_index],x[obj2_index],0,0,0,&tree);
			      	if(insert_check)
			      	{
			      		num_inserted++;
			      		for(i=0;i<cur_numcols-2;i++)
			      		{
			      			stored_x[x_rotation][i] = x[i];
			      		}
			      		x_rotation = (x_rotation + 1) % num_x_to_store;
			      		insert_check = 0;
			      		PSA_full(env,NULL,x,NULL,NULL);
			      	}
		      	}
		      	
			if(pure_binary == 1) run_it += binary_heuristic(env,lp1,lp2,indices);
			else run_it += mixed_heuristic(env);
		}
	
		box *temp = first_box;
		box *next = first_box->next;
		for(j=0;j<num_boxes;j++)
		{
			free(temp);
			temp = next;
			if(j != num_boxes-1) next = next->next;
		}
		
		END_PREPROCESSING1:
		
		finish_presolve = clock();
		duration_presolve = (double)(finish_presolve - start_presolve) / CLOCKS_PER_SEC;
	  	duration_BB = (double)(finish_presolve - start_BB) / CLOCKS_PER_SEC;
	  	
  		printf("*******************\n");
  		printf("Prepopulate time: (Alone) %lf\t (Cumulative) %lf\nIterations: %d, MIPs solved: %d Num sols added from heur: %d\n",
  									duration_presolve,duration_BB,presolve_iteration,mips_solved,run_it);
  		if(!suppress_file_output) fprintf(bb_results,"%lf\t",duration_BB);
  		printf("*******************\n");
  		
	  	if(x) free_and_null ((char **) &x);
	  	status = CPXfreeprob (env, &presolve_lp) || CPXfreeprob (env, &presolve_lp2);
    		if ( status ) 
    		{
      			printf ("CPXfreeprob failed, error code %d.\n",status);
    		}
    	}
    	else
    	{
    		TRY_ANOTHER_PREPROCESSING:
    		
    		start_presolve = clock();
    		double epsilon = 0.;
  		double step = 0.;
	/*	printf("beginning presolve\n");*/
		CPXLPptr presolve_lp = CPXcloneprob (env_just_solve_mips, lp1, &status);
	  	if ( status ) 
	  	{
	    		fprintf (stderr, "Failed to clone problem 1.\n");
	    		return status;
	  	}
	  	CPXLPptr presolve_lp2 = CPXcloneprob (env_just_solve_mips, presolve_lp, &status);
	  	if ( status ) 
	  	{
	    		fprintf (stderr, "Failed to clone problem 1.\n");
	    		return status;
	  	}
	  	CPXLPptr presolve_lp_ws = CPXcloneprob (env_just_solve_mips, presolve_lp, &status);
	  	if ( status ) 
	  	{
	    		fprintf (stderr, "Failed to clone problem 1.\n");
	    		return status;
	  	}
		chg_coefs_(env_just_solve_mips,presolve_lp2,indices,-.01);
		chg_coefs_(env_just_solve_mips,presolve_lp_ws,indices,-1.);
		chg_coefs_(env_just_solve_mips,presolve_lp,indices,-100.);
	
		double ubs[2] = {SE_extreme_x,NW_extreme_y}; 
		double lbs[2] = {NW_extreme_x,SE_extreme_y};
	
		status = CPXgetlb (env, lp1, lbs, obj1_index, obj2_index);
		if ( status ) {
			printf ("Failed to get lb's for objectives. Error code %d\n", status);
			return status;
		}
		status = CPXgetub (env, lp1, ubs, obj1_index, obj2_index);
		if ( status ) {
			printf ("Failed to get ub's for objectives. Error code %d\n", status);
			return status;
		}
		
		step = (ubs[1] - lbs[1])/(num_mips_to_solve*3);
	
	  	char lu[1] = {'L'};
	    	double *x = (double *) malloc (cur_numcols*sizeof(double));
	    	int presolve_iteration = 1;
	    	if(heur_limit) heur_limit = 12;
	    	
	    	int num_starts = CPXgetnummipstarts (env_just_solve_mips, global_mip);
	    	
	    	if(num_starts > global_num_starts)
	    	{
/*			    		printf("reallocation: (%d)\n",__LINE__);*/
	    		global_num_starts = num_starts;
	    		global_startspace = cur_numcols*global_num_starts;
	    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
			global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
			global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
			global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
	    	}
	    	
	    	int nzcnt = 0, surplus = 0;

		status = CPXgetmipstarts (env_just_solve_mips, global_mip, &nzcnt, global_beg, global_varindices, 
				   global_values, global_effortlevel, global_startspace,
				   &surplus, 0, num_starts-1);
				   
		status = CPXaddmipstarts (env_just_solve_mips, presolve_lp_ws, num_starts, nzcnt, global_beg, global_varindices,
				   global_values, global_effortlevel, NULL);
	    	 
		CPXmipopt (env_just_solve_mips, presolve_lp_ws);
		
		num_starts = CPXgetnummipstarts (env_just_solve_mips, presolve_lp_ws);
/*			  	printf("number of mip starts: %d\n",num_starts);*/
	    	
	    	if(num_starts > global_num_starts)
	    	{
/*			    		printf("reallocation: (%d)\n",__LINE__);*/
	    		global_num_starts = num_starts;
	    		global_startspace = cur_numcols*global_num_starts;
	    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
			global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
			global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
			global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
	    	}

		status = CPXgetmipstarts (env_just_solve_mips, presolve_lp_ws, &nzcnt, global_beg, global_varindices, 
				   global_values, global_effortlevel, global_startspace,
				   &surplus, 0, num_starts-1);
				   
		status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,
				   global_values, global_effortlevel, NULL);
		
	    	int status2 = CPXgetstat(env_just_solve_mips, presolve_lp_ws);
	    	int num_inserted = 0;
	    	int mips_solved = 1;
	    	double usercut_rhs = 0.;
	    	
/*			    	printf("status2: %d\n",status2);*/
	    	
	    	if(status2 != 101 && status2 != 102 && status2 != 129 && status2 != 130)
	    	{
			if(status2 == 118)
	    		{
	/*    			printf("The first test problem is unbounded. Exiting!\n");*/
	    			return status;
	    		}
	    		else if(status2 == 104 || status2 == 105 || status2 == 107 || status2 == 109 || status2 == 111 || status2 == 113 ||
					status2 == 128)
	    		{
	 			status = CPXgetbestobjval (env_just_solve_mips, presolve_lp_ws, &usercut_rhs);
	 			if(status)
	 			{
	 				printf("(%d) Failed to get best objval. Error code %d\n",__LINE__,status);
	 				return status;
	 			}
	 		}
	    	}
	    	else status = CPXgetobjval (env_just_solve_mips, presolve_lp_ws, &usercut_rhs);
	    	
	    	int z = 0;
	    	double usercut_coefs[2] = {1.,1.};
	    	int usercut_indices[2] = {obj1_index,obj2_index};
	    	char usercut_sense[1] = {'L'};
	    	
	    	double obj_off;
    		status = CPXgetobjoffset(env_just_solve_mips, presolve_lp_ws, &obj_off );
    		
    		usercut_rhs = usercut_rhs - obj_off + ob1_offset + ob2_offset;
    		
/*		    		double tx = -6295.;*/
/*					double tx2 = -754.;*/
/*			    	double ty = (tx - usercut_rhs)*(-1.);*/
/*			    	double ty2 = (tx2 - usercut_rhs)*(-1.);*/
/*	    	*/
/*			    	printf("obj offset: %lf, another possibility: %lf\n",obj_off, obj1_extra_val+ obj2_extra_val);*/
/*	    	*/
/*			    	printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",tx,tx2,ty,ty2);*/
    		
/*    		usercut_rhs -= obj_off;*/
    		
/*		    		ty = (tx - usercut_rhs)*(-1.);*/
/*			    	ty2 = (tx2 - usercut_rhs)*(-1.);*/
/*			    	*/
/*			    	printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",tx,tx2,ty,ty2);*/
    		
/*    		usercut_rhs += obj1_extra_val+ obj2_extra_val;*/
    		
/*			    	ty = (tx - usercut_rhs)*(-1.);*/
/*			    	ty2 = (tx2 - usercut_rhs)*(-1.);*/
	    	
/*			    	printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",tx,tx2,ty,ty2);*/

	    	if(add_local_cuts_during_preprocessing)
		{
		    	status = CPXaddusercuts (env, lp1, 1, 2, &usercut_rhs, usercut_sense, &z, usercut_indices, usercut_coefs, NULL);
		    	if(status)
		  	{
				printf("(%d) Failed to user cut. Error code %d\n",__LINE__,status);
		  		return status;
		  	}
	  	}
	  	
/*			  	double tx = -4475.;*/
/*			    	double tx2 = -3925.;*/
/*			    	double ty = (tx - usercut_rhs)*(-1.);*/
/*			    	double ty2 = (tx2 - usercut_rhs)*(-1.);*/
/*			    	*/
/*			    	printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",tx,tx2,ty,ty2);*/
	  	
/*			  	printf("(%d) adding user cut: x%d + %lfx%d <= %lf\n",__LINE__,obj1_index,1.,obj2_index,usercut_rhs);*/

		int num_solns = CPXgetsolnpoolnumsolns (env_just_solve_mips, presolve_lp_ws);
	  	
	  	if(num_solns >= prev_numsolns) times_to_run = num_solns - prev_numsolns;
	  	else times_to_run = num_solns;
	  	prev_numsolns = num_solns;
	  	
	  	int insert_check = 0;
	  	
	  	for(j=0;j<times_to_run;j++)
	  	{
	  		status = CPXgetsolnpoolx (env_just_solve_mips, presolve_lp_ws, j, x, 0, cur_numcols-1);
/*			  		printf("plot(%lf,%lf,'go');\n",x[obj1_index],x[obj2_index]);*/
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
		      		PSA_full(env,NULL,x,NULL,NULL);
		      	}
	      	}
	      	
	      	epsilon = x[obj2_index] - step;
	      	double limit_x = lbs[0];
	      	double limit_y = lbs[1];
	      	double future_epsilon = x[obj1_index] - step;
	      	double future_step;
	      	future_step = (ubs[0] - lbs[0])/(num_mips_to_solve*3);
	      	
		int run_it = 0;
		if(pure_binary == 1) run_it += binary_heuristic(env,lp1,lp2,indices);
		else run_it += mixed_heuristic(env);

		int num_inserted2 = 0;
		int j = 0;
		double prev_val = 0.;
		int prev_numsols = 0, numsolns = 0, numrep = 0;

		while(presolve_iteration < 3 || (presolve_iteration < 100 && epsilon > limit_y))
		{
/*					printf("epsilon: %lf\n",epsilon);*/
			presolve_iteration++;
			
			finish_presolve = clock();
			preprocessing_time = (double)(finish_presolve - start_presolve) / CLOCKS_PER_SEC;
			if(preprocessing_time > max_preprocessing_time) 
			{
				printf("Preprocessing exceeded time limit of %lf. Moving to presolve phase 2.\n",max_preprocessing_time);
				goto END_PREPROCESSING2;
			}
			
			status = CPXchgbds (env_just_solve_mips, presolve_lp, 1, &obj2_index, lu, &epsilon);

			num_starts = CPXgetnummipstarts (env_just_solve_mips, global_mip);
		    	
		    	if(num_starts > global_num_starts)
		    	{
/*				    		printf("reallocation: (%d)\n",__LINE__);*/
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
					   
			status = CPXaddmipstarts (env_just_solve_mips, presolve_lp, num_starts, nzcnt, global_beg, global_varindices,
					   global_values, global_effortlevel, NULL);

            status = CPXwriteprob (env_just_solve_mips, presolve_lp, "myprob.lp", "LP");
/*            exit(0);*/
			CPXmipopt (env_just_solve_mips, presolve_lp);
			
			numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, presolve_lp);
			numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, presolve_lp);
			
			num_starts = numsolns - prev_numsols + numrep;
/*			  		printf("number of mip starts to use: %d\n",num_starts);*/
	  		
	  		prev_numsols = numsolns;
	  		
		    	if(num_starts > global_num_starts)
		    	{
/*				    		printf("reallocation: (%d)\n",__LINE__);*/
		    		global_num_starts = num_starts;
		    		global_startspace = cur_numcols*global_num_starts;
		    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
				global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
				global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
				global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
		    	}

			status = CPXgetmipstarts (env_just_solve_mips, presolve_lp, &nzcnt, global_beg, global_varindices, 
					   global_values, global_effortlevel, global_startspace,
					   &surplus, 0, num_starts-1);
					   
			status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,
					   global_values, global_effortlevel, NULL);
			
			if(!exact_mips) CPXpopulate(env_just_solve_mips,presolve_lp);
    			status2 = CPXgetstat(env_just_solve_mips, presolve_lp);
    			int numsols = CPXgetsolnpoolnumsolns (env_just_solve_mips, presolve_lp);

			double val = 0.;
			
			if(numsols >= prev_numsolns) times_to_run = numsols - prev_numsolns;
		  	else times_to_run = numsols;
		  	prev_numsolns = numsols;
			
			num_inserted2 = 0;
			for(j=0;j<times_to_run;j++)
			{
				status = CPXgetsolnpoolx (env_just_solve_mips, presolve_lp, j, x, 0, cur_numcols-1);
				insert_check = 0;
				insert_check = mock_insert(1,x[obj1_index],x[obj2_index],0,0,0,&tree);
				if(insert_check)
				{
				      	num_inserted++;
				      	num_inserted2++;
				      	for(i=0;i<cur_numcols-2;i++)
				      	{
				      		stored_x[x_rotation][i] = x[i];
				/*     		printf("stored val: %lf\n",stored_x[x_rotation][i]);*/
			  		}
				      	x_rotation = (x_rotation + 1) % num_x_to_store;
				      	insert_check = 0;
					PSA_full(env,NULL,x,NULL,NULL);
			      	}
			}
			
			prev_val = x[obj2_index];
			if(num_inserted2)
			{
			      	if(step > (ub[1] - lb[1])/100.) step = step/(1.+ (double) preprocessing_parameter);
			}
		      	else
		      	{
		      		if(step < (ub[1] - lb[1])/10.) step = step*fmax(5.-(double) preprocessing_parameter, 1.);
		      	}
					
			mips_solved++;				      	
		      	epsilon -= step;
		      	
			if(pure_binary == 1) run_it += binary_heuristic(env,lp1,lp2,indices);
			else run_it += mixed_heuristic(env);
/*			printf("presolve_iteration: %d\n",presolve_iteration);*/
		}
		
		epsilon = future_epsilon;
		step = future_step;
		prev_numsols = 0;
		
/*				printf("starting 2nd half of prepopulate\n");*/
		while(presolve_iteration < 3 || (presolve_iteration < 200 && epsilon > limit_x))
		{
/*					printf("epsilon: %lf\n",epsilon);*/
			presolve_iteration++;
			
			finish_presolve = clock();
			preprocessing_time = (double)(finish_presolve - start_presolve) / CLOCKS_PER_SEC;
			if(preprocessing_time > max_preprocessing_time) 
			{
				printf("Preprocessing exceeded time limit of %lf. Moving to presolve phase 2.\n",max_preprocessing_time);
				goto END_PREPROCESSING2;
			}
			
			status = CPXchgbds (env_just_solve_mips, presolve_lp2, 1, &obj1_index, lu, &epsilon);

			num_starts = CPXgetnummipstarts (env_just_solve_mips, global_mip);
		    	
		    	if(num_starts > global_num_starts)
		    	{
/*				    		printf("reallocation: (%d)\n",__LINE__);*/
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
					   
			status = CPXaddmipstarts (env_just_solve_mips, presolve_lp2, num_starts, nzcnt, global_beg, global_varindices,
					   global_values, global_effortlevel, NULL);

			CPXmipopt (env_just_solve_mips, presolve_lp2);
			
			numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, presolve_lp2);
			numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, presolve_lp2);
			
			num_starts = numsolns - prev_numsols + numrep;
/*			  		printf("number of mip starts to use: %d\n",num_starts);*/
	  		
	  		prev_numsols = numsolns;
	  		
		    	if(num_starts > global_num_starts)
		    	{
/*				    		printf("reallocation: (%d)\n",__LINE__);*/
		    		global_num_starts = num_starts;
		    		global_startspace = cur_numcols*global_num_starts;
		    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
				global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
				global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
				global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
		    	}

			status = CPXgetmipstarts (env_just_solve_mips, presolve_lp2, &nzcnt, global_beg, global_varindices, 
					   global_values, global_effortlevel, global_startspace,
					   &surplus, 0, num_starts-1);
					   
			status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,
					   global_values, global_effortlevel, NULL);
			
			if(!exact_mips) CPXpopulate(env_just_solve_mips,presolve_lp2);
    			status2 = CPXgetstat(env_just_solve_mips, presolve_lp2);
    			numsols = CPXgetsolnpoolnumsolns (env_just_solve_mips, presolve_lp2);

			if(numsols >= prev_numsolns) times_to_run = numsols - prev_numsolns;
		  	else times_to_run = numsols;
		  	prev_numsolns = numsols;

			num_inserted2 = 0;
			for(j=0;j<times_to_run;j++)
			{
/*						printf("(%d) j: %d\n",__LINE__,j);*/
				status = CPXgetsolnpoolx (env_just_solve_mips, presolve_lp2, j, x, 0, cur_numcols-1);
				insert_check = 0;
				insert_check = mock_insert(1,x[obj1_index],x[obj2_index],0,0,0,&tree);
				if(insert_check)
				{
				      	num_inserted++;
				      	num_inserted2++;
				      	for(i=0;i<cur_numcols-2;i++)
				      	{
				      		stored_x[x_rotation][i] = x[i];
			  		}
				      	x_rotation = (x_rotation + 1) % num_x_to_store;
				      	insert_check = 0;
					PSA_full(env,NULL,x,NULL,NULL);
			      	}
			}
			
			prev_val = x[obj2_index];
			
/*					printf("num inserted 2: %d\n",num_inserted2);*/
			if(num_inserted2)
			{
			      	if(step > (ub[1] - lb[1])/100.) step = step/(1.+ (double) preprocessing_parameter);
			}
		      	else
		      	{
		      		if(step < (ub[1] - lb[1])/10.) step = step*fmax(5.-(double) preprocessing_parameter, 1.);
		      	}
					
			mips_solved++;
		      	epsilon -= step;
		      	
			if(pure_binary == 1) run_it += binary_heuristic(env,lp1,lp2,indices);
			else run_it += mixed_heuristic(env);
/*			printf("presolve_iteration2: %d\n",presolve_iteration);*/
		}
		
  		END_PREPROCESSING2:
	
		finish_presolve = clock();
		duration_presolve = (double)(finish_presolve - start_presolve) / CLOCKS_PER_SEC;
	  	duration_BB = (double)(finish_presolve - start_BB) / CLOCKS_PER_SEC;
  		printf("*******************\n");
  		printf("Prepopulate time: (Alone) %lf\t (Cumulative) %lf\nIterations: %d, MIPs solved: %d Num sols added from heur: %d\n",
  									duration_presolve,duration_BB,presolve_iteration,mips_solved,run_it);
  		if(!suppress_file_output) fprintf(bb_results,"%lf\t",duration_BB);
  		printf("*******************\n");
  		
	  	if(x) free_and_null ((char **) &x);
	  	status = CPXfreeprob (env, &presolve_lp) || CPXfreeprob (env, &presolve_lp2) || CPXfreeprob (env, &presolve_lp_ws);
    		if ( status ) 
    		{
      			printf ("CPXfreeprob failed, error code %d.\n",status);
    		}
    	}
	return 0;
}

int weighted_sum_preprocessing(clock_t start_presolve, CPXLPptr lp1, CPXLPptr lp2, CPXENVptr env, int cur_numcols, int cur_numrows, int *indices, 
				double *obj_coef2, int intervals, int run_it, double split_pt_denom, clock_t finish_presolve, double preprocessing_time, 
				double duration_presolve, FILE *bb_results, double *obj_coef1, int suppress_file_output)
{
	start_presolve = clock();
	int num_slopes = 0, status = 0, nzcnt = 0, surplus = 0, status2 = 0, i = 0;
	int j;
/*	printf("beginning presolve alt\n");*/

	if(preprocessing_parameter < 0) 
		preprocessing_parameter = (int) (3.*(.185*cur_numcols)*(.185*cur_numcols)/((.185*cur_numcols)*(.185*cur_numcols)+1.5*cur_numcols+10000.));

	CPXLPptr presolve_lp = CPXcloneprob (env_just_solve_mips, lp1, &status);
  	if ( status ) 
  	{
    		fprintf (stderr, "Failed to clone problem 1.\n");
    		return status;
  	}
  	
/*		  	status = CPXwriteprob (env_just_solve_mips, presolve_lp, "myprob.lp", "LP");*/
/*		  	int numc = CPXgetnumcols (env_just_solve_mips, presolve_lp);*/
/*		  	printf("curnumcols: %d, obj2index: %d\n",numc,obj2_index);*/
/*		  	exit(0);*/
	
	double ubs[2] = {0.,0.};
	double lbs[2] = {0.,0.};
	
	status = CPXgetlb (env, lp1, lbs, obj1_index, obj2_index);
	if ( status ) 
  	{
  		printf("curnumcols: %d, obj2index: %d\n",cur_numcols,obj2_index);
    		printf ("(%d) Failed to get lbs. Error code %d\n",__LINE__,status);
    		return status;
  	}
	status = CPXgetub (env, lp1, ubs, obj1_index, obj2_index);
	if ( status ) 
  	{
  		printf("curnumcols: %d, obj2index: %d\n",cur_numcols,obj2_index);
    		printf ("(%d) Failed to get ubs. Error code %d\n",__LINE__,status);
    		return status;
  	}

	int num_starts = CPXgetnummipstarts (env_just_solve_mips, global_mip);
	int prev_numsols = 0;
		    	
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
	if ( status ) 
  	{
    		printf ("(%d) Failed to get mip starts. Error code %d\n",__LINE__,status);
/*		    		return status;*/
  	}
			   
	status = CPXaddmipstarts (env_just_solve_mips, presolve_lp, num_starts, nzcnt, global_beg, global_varindices,
			   global_values, global_effortlevel, NULL);
	if ( status ) 
  	{
    		printf ("(%d) Failed to add mip starts. Error code %d\n",__LINE__,status);
/*		    		return status;*/
  	}

	if(exact_mips) CPXmipopt (env_just_solve_mips, presolve_lp);
	else CPXpopulate (env_just_solve_mips, presolve_lp);
	
	int numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, presolve_lp);
	int numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, presolve_lp);
	
	num_starts = numsolns - prev_numsols + numrep;
/*	  		printf("number of mip starts to use: %d\n",num_starts);*/
	
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

	status = CPXgetmipstarts (env_just_solve_mips, presolve_lp, &nzcnt, global_beg, global_varindices, 
			   global_values, global_effortlevel, global_startspace,
			   &surplus, 0, num_starts-1);
	if ( status ) 
  	{
    		printf ("(%d) Failed to get mip starts. Error code %d\n",__LINE__,status);
/*		    		return status;*/
  	}
			   
	status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,
			   global_values, global_effortlevel, NULL);
	if ( status ) 
  	{
    		printf ("(%d) Failed to add mip starts. Error code %d\n",__LINE__,status);
/*		    		return status;*/
  	}
	
/*			CPXpopulate(env,presolve_lp);*/
    	status2 = CPXgetstat(env_just_solve_mips, presolve_lp);
/*		    	printf("status2: %d\n",status2);*/
    	
    	int numsols = CPXgetsolnpoolnumsolns (env_just_solve_mips, presolve_lp);
    	double objvals[2] = {0.,0.};

	int insert_check = 0;
	double *x = (double *) malloc (cur_numcols*sizeof(double));
	
	if(numsols >= prev_numsolns) times_to_run = numsols - prev_numsolns;
  	else times_to_run = numsols;
  	prev_numsolns = numsols;

/*			printf("solutions found for objective 1: \n");*/
    	if(status2 != 103) for(i=0;i<times_to_run;i++)
    	{
    		status = CPXgetsolnpoolx (env_just_solve_mips, presolve_lp, i, objvals, obj1_index, obj2_index);
    		if ( status ) 
	  	{
	    		printf ("(%d) Failed to get x's. Error code %d\n",__LINE__,status);
	    		return status;
	  	}
/*		    		printf("plot(%lf,%lf,'go');\n",objvals[0],objvals[1]);*/
    		if(i == 0)
    		{ 
    			ubs[0] = objvals[0];
    			lbs[1] = objvals[1];
    		}
    		insert_check = mock_insert(1,objvals[0],objvals[1],0,0,0,&tree); 
    		if(insert_check)
	      	{
	      		status = CPXgetsolnpoolx (env_just_solve_mips, presolve_lp, i, x, 0, cur_numcols-2);
	      		if ( status ) 
		  	{
		    		printf ("(%d) Failed to get x's. Error code %d\n",__LINE__,status);
		    		return status;
		  	}
/*			      		printf("status %d\n",status);*/
	      		for(j=0;j<cur_numcols-2;j++)
	      		{
	      			stored_x[x_rotation][j] = x[j];
	      		}
	      		x_rotation = (x_rotation + 1) % num_x_to_store;
	      		insert_check = 0;
	      		PSA_full(env,NULL,x,NULL,NULL);
	      	}
    	}
    	
    	status = CPXchgobj (env_just_solve_mips, presolve_lp, obj1_index+2, indices, obj_coef2);
	if ( status ) {
		printf ("Failed to get change obj coef. Error code %d\n", status);
		return status;
	}
	
	if(exact_mips) CPXmipopt (env_just_solve_mips, presolve_lp);
	else CPXpopulate (env_just_solve_mips, presolve_lp);
	
	status2 = CPXgetstat(env_just_solve_mips, presolve_lp);
/*		    	printf("status2: %d\n",status2);*/
			
	numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, presolve_lp);
	numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, presolve_lp);
	
	num_starts = numsolns - prev_numsols + numrep;
/*	  		printf("number of mip starts to use: %d\n",num_starts);*/
	
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

	if(num_starts > 0)
	{
		status = CPXgetmipstarts (env_just_solve_mips, presolve_lp, &nzcnt, global_beg, global_varindices, 
				   global_values, global_effortlevel, global_startspace,
				   &surplus, 0, num_starts-1);
		if ( status ) 
	  	{
	    		printf ("(%d) Failed to get mip starts. Error code %d\n",__LINE__,status);
/*			    		return status;*/
	  	}
				   
		status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,
				   global_values, global_effortlevel, NULL);
		if ( status ) 
	  	{
	    		printf ("(%d) Failed to add mip starts. Error code %d\n",__LINE__,status);
/*			    		return status;*/
	  	}
	
	    	status2 = CPXgetstat(env_just_solve_mips, presolve_lp);
	    	numsols = CPXgetsolnpoolnumsolns (env_just_solve_mips, presolve_lp);
	    	
	    	if(numsols >= prev_numsolns) times_to_run = numsols - prev_numsolns;
	  	else times_to_run = numsols;
	  	prev_numsolns = numsols;
  	}
  	else times_to_run = 0;
	
/*			printf("solutions found for objective 2: \n");*/
    	for(i=0;i<times_to_run;i++)
    	{
    		status = CPXgetsolnpoolx (env_just_solve_mips, presolve_lp, i, objvals, obj1_index, obj2_index);
    		if ( status ) 
	  	{
	  		printf("curnumcols: %d, obj2index: %d\n",cur_numcols,obj2_index);
	    		printf ("(%d) Failed to get objvals. Error code %d\n",__LINE__,status);
	    		return status;
	  	}
/*		    		printf("plot(%lf,%lf,'go');\n",objvals[0],objvals[1]);*/
    		if(i == 0)
    		{ 
    			ubs[1] = objvals[1];
    			lbs[0] = objvals[0];
    		}
    		insert_check = mock_insert(1,objvals[0],objvals[1],0,0,0,&tree); 
    		if(insert_check)
	      	{
	      		status = CPXgetsolnpoolx (env_just_solve_mips, presolve_lp, i, x, 0, obj1_index);
	      		if ( status ) 
		  	{
		    		printf ("(%d) Failed to get x's. Error code %d\n",__LINE__,status);
		    		return status;
		  	}
	      		for(j=0;j<cur_numcols-2;j++)
	      		{
	      			stored_x[x_rotation][j] = x[j];
	      		}
	      		x_rotation = (x_rotation + 1) % num_x_to_store;
	      		insert_check = 0;
	      		PSA_full(env,NULL,x,NULL,NULL);
	      	}
    	}
    	
    	double prob_width,prob_height;
    	
/*		    	if(exact_mips)*/
/*		    	{*/
    		prob_width = ubs[0] - lbs[0];
    		prob_height = ubs[1] - lbs[1];
/*		    	}*/
/*		    	else*/
/*		    	{*/
/*		    		prob_width = x1_[0] - x2_[0];*/
/*		    		prob_height = x2_[1] - x1_[1];*/
/*		    	}*/
    	
/*		    	printf("prob height: %lf, width: %lf\n",prob_height,prob_width);*/
    	
    	interval *first_h_interval = NULL;
    	interval *first_v_interval = NULL;
    	int num_h_intervals = 0;
    	int num_v_intervals = 0;
    	
    	if(intervals)
    	{
    		first_h_interval = (struct interval*) malloc( sizeof( struct interval ) );
/*		    		if(exact_mips)*/
/*		    		{*/
    			first_h_interval->lb = lbs[0];
    			first_h_interval->ub = ubs[0];
/*		    		}*/
/*		    		else*/
/*		    		{*/
/*		    			first_h_interval->lb = x2_[0];//-prob_width;*/
/*		    			first_h_interval->ub = x1_[0];//+prob_width;*/
/*		    		}*/
    		first_h_interval->dir = 0;
    		first_h_interval->miss_count = 0;
    		first_h_interval->next = NULL;
    		num_h_intervals++;
    		first_v_interval = (struct interval*) malloc( sizeof( struct interval ) );
/*		    		if(exact_mips)*/
/*		    		{*/
	    		first_v_interval->lb = lbs[1];
	    		first_v_interval->ub = ubs[1];
/*		    		}*/
/*		    		else*/
/*		    		{*/
/*		    			first_v_interval->lb = x1_[1]-5*prob_height;*/
/*			    		first_v_interval->ub = x2_[1];//+prob_height;*/
/*		    		}*/
    		first_v_interval->dir = 1;
    		first_v_interval->miss_count = 0;
    		first_v_interval->next = NULL;
    		num_v_intervals++;
    	}
    	
    	slope_val *first_slope = (struct slope_val*) malloc( sizeof( struct slope_val ) );
	first_slope->slope = -(ubs[1]-lbs[1])/(ubs[0]-lbs[0]);
	first_slope->ub = 0.;
	first_slope->lb = -12.5;
	first_slope->next = NULL;
	first_slope->prev = NULL;
	num_slopes++;
	
    	int presolve_iteration = 1;
    	if(heur_limit) heur_limit = 12;
      	
	run_it = 0;
	if(pure_binary == 1) run_it += binary_heuristic(env_just_solve_mips,lp1,lp2,indices);
	else run_it += mixed_heuristic(env);
	
	int cur_num_slopes = num_slopes;
	
	int num_inserted = cur_num_slopes;
	int mips_solved = 2;
	
	int pop_lim = 200, times_num_inserted_below = 0;

/*			printf("solutions found during weighted sum search:\n");*/
	while(presolve_iteration < 40 && times_num_inserted_below <= preprocessing_parameter) //num_inserted >= cur_num_slopes/split_pt_denom)
	{
		cur_num_slopes = num_slopes;
		presolve_iteration++;
/*		printf("presolve_iteration: %d, num_inserted: %d, cur_num_slopes/split_pt_denom: %lf\n",presolve_iteration, num_inserted, cur_num_slopes/split_pt_denom);*/
		
		finish_presolve = clock();
		preprocessing_time = (double)(finish_presolve - start_presolve) / CLOCKS_PER_SEC;
		if(preprocessing_time > max_preprocessing_time) 
		{
			printf("Preprocessing exceeded time limit of %lf. Moving to presolve phase 2.\n",max_preprocessing_time);
			goto END_PREPROCESSING3;
		}
		
		num_inserted = 0;
		slope_val *temp_slope = first_slope;
		for(j=0;j<cur_num_slopes;j++)
		{	
/*					printf("slope: %lf\n",temp_slope->slope);*/
			slope_val *next = temp_slope->next;
			chg_coefs_(env_just_solve_mips,presolve_lp,indices,temp_slope->slope);
/*					exit(0);*/
			
			num_starts = CPXgetnummipstarts (env_just_solve_mips, global_mip);
			prev_numsols = 0;
				    	
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
			if ( status ) 
		  	{
		    		printf ("(%d) Failed to get mip starts. Error code %d\n",__LINE__,status);
/*		    		return status;*/
		  	}
					   
			status = CPXaddmipstarts (env_just_solve_mips, presolve_lp, num_starts, nzcnt, global_beg, global_varindices,
					   global_values, global_effortlevel, NULL);
			if ( status ) 
		  	{
		    		printf ("(%d) Failed to add mip starts. Error code %d\n",__LINE__,status);
/*		    		return status;*/
		  	}
			
			if(exact_mips) CPXmipopt (env_just_solve_mips, presolve_lp);
			else CPXpopulate (env_just_solve_mips, presolve_lp);
			
			numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, presolve_lp);
			numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, presolve_lp);
			
			num_starts = numsolns - prev_numsols + numrep;
/*			  		printf("number of mip starts to use: %d\n",num_starts);*/
	  		
	  		prev_numsols = numsolns;
	  		
		    	if(num_starts > global_num_starts)
		    	{
/*				    		printf("reallocation: (%d)\n",__LINE__);*/
		    		global_num_starts = num_starts;
		    		global_startspace = cur_numcols*global_num_starts;
		    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
				global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
				global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
				global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
		    	}

			status = CPXgetmipstarts (env_just_solve_mips, presolve_lp, &nzcnt, global_beg, global_varindices, 
					   global_values, global_effortlevel, global_startspace,
					   &surplus, 0, num_starts-1);
			if ( status ) 
		  	{
		    		printf ("(%d) Failed to get mip starts. Error code %d\n",__LINE__,status);
/*				    		return status;*/
		  	}
					   
			status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, global_varindices,
					   global_values, global_effortlevel, NULL);
			if ( status ) 
		  	{
		    		printf ("(%d) Failed to add mip starts. Error code %d\n",__LINE__,status);
/*				    		return status;*/
		  	}
			
		    	status2 = CPXgetstat(env_just_solve_mips, presolve_lp);
/*				    	printf("status2: %d\n", status2);*/
		    	
			if(status2 == 128)
			{
				pop_lim = 2*pop_lim;
				status = CPXsetintparam (env_just_solve_mips, CPX_PARAM_POPULATELIM, pop_lim);
				if ( status ) {
					printf ("Failed to set populate limit to 200, error %d.\n",status);
					return status;
				}
			}

			double usercut_rhs = 0.;
	    	
		    	if(status2 != 101 && status2 != 102 && status2 != 129 && status2 != 130)
		    	{
/*		    		printf("approximating\n");*/
				if(status2 == 104 || status2 == 105 || status2 == 107 || status2 == 109 || status2 == 111 || status2 == 113 ||
					status2 == 128)
		    		{
		 			status = CPXgetbestobjval (env_just_solve_mips, presolve_lp, &usercut_rhs);
		 			if(status)
		 			{
		 				printf("(%d) Failed to get best objval. Error code %d\n",__LINE__,status);
		 				return status;
		 			}
		 		}
		    	}
		    	else 
			{
/*				printf("true value\n");*/
				status = CPXgetobjval (env_just_solve_mips, presolve_lp, &usercut_rhs);
			    	if ( status ) 
			  	{
			    		printf ("(%d) Failed to get obj val. Error code %d\n",__LINE__,status);
			    		return status;
			  	}
/*			  	printf("objval: %lf\n",usercut_rhs);*/
/*			  	status = CPXgetx (env_just_solve_mips, presolve_lp, X, 0, obj2_index);*/
/*	  			if(status)*/
/*	  			{*/
/*					printf("Failed to get x-values from CPLEX. Status: %d Line: %d\n",status,__LINE__);*/
/*					exit(0);*/
/*	  				return status;*/
/*	  			}*/
		    	}
		    	
		    	int z = 0;
		    	double usercut_coefs[2] = {1.,-1./temp_slope->slope};
		    	int usercut_indices[2] = {obj1_index,obj2_index};
		    	char usercut_sense[1] = {'L'};
		    	
		    	double obj_off;
	    		status = CPXgetobjoffset(env_just_solve_mips, presolve_lp, &obj_off );
	    		if ( status ) 
		  	{
		    		printf ("(%d) Failed to get obj offset. Error code %d\n",__LINE__,status);
		    		return status;
		  	}
				  	
			usercut_rhs = usercut_rhs - obj_off + ob1_offset - 1./temp_slope->slope*ob2_offset;
	    		
/*			    		double tx = -6295.;*/
/*					double tx2 = -754.;*/
/*				    	double ty = (tx - usercut_rhs)*temp_slope->slope;*/
/*				    	double ty2 = (tx2 - usercut_rhs)*temp_slope->slope;*/
/*				    	*/
/*				    	printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",tx,tx2,ty,ty2);*/
/*	    		*/
/*	    		usercut_rhs -= obj_off;*/
/*	    		*/
/*			    		ty = (tx - usercut_rhs)*temp_slope->slope;*/
/*				    	ty2 = (tx2 - usercut_rhs)*temp_slope->slope;*/
/*				    	*/
/*				    	printf("plot([%lf,%lf],[%lf,%lf],'-co');\n",tx,tx2,ty,ty2);*/
/*				    	*/
/*			usercut_rhs += 2*obj_off;*/
/*	    		*/
/*			    		ty = (tx - usercut_rhs)*temp_slope->slope;*/
/*				    	ty2 = (tx2 - usercut_rhs)*temp_slope->slope;*/
/*				    	*/
/*				    	printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",tx,tx2,ty,ty2);*/
/*	    		*/
/*	    		printf("another number: %lf\n",obj1_extra_val - temp_slope->slope*obj2_extra_val);*/
/*	    		*/
/*			usercut_rhs += obj1_extra_val - temp_slope->slope*obj2_extra_val;*/
/*			    		*/
/*				    	ty = (tx - usercut_rhs)*temp_slope->slope;*/
/*				    	ty2 = (tx2 - usercut_rhs)*temp_slope->slope;*/
/*				    	*/
/*				    	printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",tx,tx2,ty,ty2);*/
		    	
		    	if(add_local_cuts_during_preprocessing)
		    	{
			    	status = CPXaddusercuts (env, lp1, 1, 2, &usercut_rhs, usercut_sense, &z, usercut_indices, 
			    							usercut_coefs, NULL);
			  	if ( status ) 
			  	{
			    		printf ("(%d) Failed to add cut. Error code %d\n",__LINE__,status);
			    		return status;
			  	}
		  	}
/*		  	exit(0);*/
		  	
/*				  	double tx = -7600.;*/
/*				    	double tx2 = -6700.;*/
/*				    	double ty = (tx - usercut_rhs)*temp_slope->slope;*/
/*				    	double ty2 = (tx2 - usercut_rhs)*temp_slope->slope;*/
/*				    	*/
/*				    	printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",tx,tx2,ty,ty2);*/
/*				  	*/
/*				  	printf("(%d) adding user cut: x%d + %lfx%d <= %lf\n",__LINE__,obj1_index,-1./temp_slope->slope,obj2_index,usercut_rhs);*/
/*				  	exit(0);*/

		    	numsols = CPXgetsolnpoolnumsolns (env_just_solve_mips, presolve_lp);

			if(numsols >= prev_numsolns) times_to_run = numsols - prev_numsolns;
		  	else times_to_run = numsols;
		  	prev_numsolns = numsols;

			int num_added = 0;
			double best_obj_val = 0.;
			double best_obj2_val = 0.;
		    	for(i=0;i<times_to_run;i++)
		    	{
/*				    		printf("(%d) i: %d\n",__LINE__,i);*/
		    		status = CPXgetsolnpoolx (env_just_solve_mips, presolve_lp, i, objvals, obj1_index, obj2_index);
		    		if ( status ) 
			  	{
			    		printf ("(%d) Failed to get x's. Error code %d\n",__LINE__,status);
			    		return status;
			  	}
/*				    		printf("plot(%lf,%lf,'go');\n",objvals[0],objvals[1]);*/
		    		best_obj_val = objvals[0];
		    		best_obj2_val = objvals[1];
/*				    		printf("plot(%lf,%lf,'o');\n",objvals[0],objvals[1]);*/
		    		insert_check = mock_insert(1,objvals[0],objvals[1],0,0,0,&tree); 
		    		if(insert_check)
			      	{
/*					      		printf("inserted: (%lf,%lf)\n",objvals[0],objvals[1]);*/
			      		num_added++;
			      		status = CPXgetsolnpoolx (env_just_solve_mips, presolve_lp, i, x, 0, cur_numcols-2);
			      		if ( status ) 
				  	{
				    		printf ("(%d) Failed to get x's. Error code %d\n",__LINE__,status);
				    		return status;
				  	}
					int k;
			      		for(k=0;k<cur_numcols-2;k++)
			      		{
			      			stored_x[x_rotation][k] = x[k];
			      		}
			      		x_rotation = (x_rotation + 1) % num_x_to_store;
			      		insert_check = 0;
			      		PSA_full(env,NULL,x,NULL,NULL);
			      	}
		    	}
		    	
		    	if(num_added) num_inserted++;
		    	
			if(intervals) 
			{
				num_h_intervals = add_interval(best_obj_val,first_h_interval,num_h_intervals);
				num_v_intervals = add_interval(best_obj2_val,first_v_interval,num_v_intervals);
			}
/*				    	num_inserted++;*/
	    		slope_val *next_slope = (struct slope_val*) malloc( sizeof( struct slope_val ) );
	    		next_slope->slope = (temp_slope->slope + temp_slope->lb)/2.;
	    		next_slope->lb = temp_slope->lb;
	    		next_slope->ub = temp_slope->slope;
	    		next_slope->next = temp_slope->next;
	    		next_slope->prev = temp_slope;
/*				    	printf("adding slope: %lf\n",next_slope->slope);*/
	    		num_slopes++;
	    		temp_slope->lb = temp_slope->slope;
	    		temp_slope->slope = (temp_slope->ub + temp_slope->slope)/2.;
/*				    	printf("modifying current slope: %lf\n",temp_slope->slope);*/
	    		if(temp_slope->next) temp_slope->next->prev = next_slope;
	    		temp_slope->next = next_slope;
			
			mips_solved++;
			temp_slope = next;
	      	}
	      	
		if(pure_binary == 1) run_it += binary_heuristic(env,lp1,lp2,indices);
		else run_it += mixed_heuristic(env);
/*		printf("presolve_iteration: %d, num_inserted: %d, cur_num_slopes/split_pt_denom: %lf\n",presolve_iteration, num_inserted, cur_num_slopes/split_pt_denom);*/
		if(num_inserted < cur_num_slopes/split_pt_denom) times_num_inserted_below++;
	}
	
	int interval_iterations = 0;
	int first_flip = 0;
	if(intervals)
	{
		int cur_num_h_intervals = num_h_intervals;
		int cur_num_v_intervals = num_v_intervals;
		
		CPXLPptr presolve_lp1 = CPXcloneprob (env_just_solve_mips, presolve_lp, &status);
	  	if ( status ) 
	  	{
	    		fprintf (stderr, "Failed to clone problem 1.\n");
	    		return status;
	  	}
	  	status = CPXchgobj (env_just_solve_mips, presolve_lp, obj1_index+2, indices, obj_coef2);
		if ( status ) {
			printf ("Failed to get change obj coef. Error code %d\n", status);
			return status;
		}
		if ( status ) 
	  	{
	    		printf ("(%d) Failed to change obj. Error code %d\n",__LINE__,status);
	    		return status;
	  	}
	
		num_inserted = cur_num_h_intervals + cur_num_v_intervals;
		int vars[2] = {obj1_index,obj1_index};
		int vars2[2] = {obj2_index,obj2_index};
		char low[2] = {'L','U'};
		double bound[2] = {0.,0.};

		while(presolve_iteration < 100 && (cur_num_h_intervals + cur_num_v_intervals) > 0 && 
		      num_inserted >= (cur_num_h_intervals+cur_num_v_intervals-2)/(split_pt_denom/1.7))
		{
			cur_num_h_intervals = num_h_intervals;
			cur_num_v_intervals = num_v_intervals;
			presolve_iteration++;
			interval_iterations++;
			
			finish_presolve = clock();
			preprocessing_time = (double)(finish_presolve - start_presolve) / CLOCKS_PER_SEC;
			if(preprocessing_time > max_preprocessing_time) 
			{
				printf("Preprocessing exceeded time limit of %lf. Moving to presolve phase 2.\n",max_preprocessing_time);
				goto END_PREPROCESSING3;
			}
			
			while(first_h_interval && (first_h_interval->miss_count >= 1 + preprocessing_parameter || 
			     (first_h_interval->ub - first_h_interval->lb) < prob_width/100.))
			{
				finish_presolve = clock();
				preprocessing_time = (double)(finish_presolve - start_presolve) / CLOCKS_PER_SEC;
				if(preprocessing_time > max_preprocessing_time) 
				{
					printf("Preprocessing exceeded time limit of %lf. Moving to presolve phase 2.\n",max_preprocessing_time);
					goto END_PREPROCESSING3;
				}
			
				interval *itvl = first_h_interval->next;
				free(first_h_interval);
				first_h_interval = itvl;
				num_h_intervals--;
/*						printf("removing first interval\n");*/
			}
			while(first_v_interval && (first_v_interval->miss_count >= 1 + preprocessing_parameter || 
			     (first_v_interval->ub - first_v_interval->lb) < prob_height/100.))
			{
				finish_presolve = clock();
				preprocessing_time = (double)(finish_presolve - start_presolve) / CLOCKS_PER_SEC;
				if(preprocessing_time > max_preprocessing_time) 
				{
					printf("Preprocessing exceeded time limit of %lf. Moving to presolve phase 2.\n",max_preprocessing_time);
					goto END_PREPROCESSING3;
				}
			
				interval *itvl = first_v_interval->next;
				free(first_v_interval);
				first_v_interval = itvl;
				num_v_intervals--;
/*						printf("removing first interval\n");*/
			}
			cur_num_h_intervals = num_h_intervals;
			cur_num_v_intervals = num_v_intervals;
			interval *temp = first_h_interval;
			for(i=0;i<cur_num_h_intervals-1;i++)
			{
/*						printf("considering interval: %lf to %lf\n", temp->next->lb, temp->next->ub);*/
				int removed = 0;
				if(temp->next->miss_count >= 2 || (temp->next->ub - temp->next->lb) < prob_width/20.)
				{
					interval *itvl = temp->next->next;
					free(temp->next);
					temp->next = itvl;
					num_h_intervals--;
					removed = 1;
				}
				if(!removed) temp = temp->next;
			}
			cur_num_h_intervals = num_h_intervals;
			temp = first_v_interval;
			for(i=0;i<cur_num_v_intervals-1;i++)
			{
/*						printf("considering interval: %lf to %lf\n", temp->next->lb, temp->next->ub);*/
				int removed = 0;
				if(temp->next->miss_count >= 2 || (temp->next->ub - temp->next->lb) < prob_height/20.)
				{
					interval *itvl = temp->next->next;
					free(temp->next);
					temp->next = itvl;
					num_v_intervals--;
					removed = 1;
				}
				if(!removed) temp = temp->next;
			}
			cur_num_v_intervals = num_v_intervals;
			temp = first_h_interval;
			for(i=0;i<num_h_intervals;i++)
			{
				temp = temp->next;
			}
			temp = first_v_interval;
			for(i=0;i<num_v_intervals;i++)
			{
				temp = temp->next;
			}
			num_inserted = 0;
			interval *temp_interval = first_h_interval;
			for(j=0;j<cur_num_h_intervals;j++)
			{	
/*						printf("interval: %lf to %lf\n",temp_interval->lb,temp_interval->ub);*/
				interval *next = NULL;
				next = temp_interval->next;
				bound[0] = (temp_interval->lb + temp_interval->ub)/2.;
				bound[1] = temp_interval->ub;
				int num_added = 0;
				double best_obj_val = 0.;
				status = CPXchgbds (env_just_solve_mips, presolve_lp, 2, vars, low, bound);
				if ( status ) 
			  	{
			    		printf ("(%d) Failed to change bound. Error code %d\n",__LINE__,status);
			    		return status;
			  	}
				if(exact_mips) CPXmipopt (env_just_solve_mips, presolve_lp);
				else CPXpopulate (env_just_solve_mips, presolve_lp);
			
				numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, presolve_lp);
				numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, presolve_lp);
			
				num_starts = numsolns - prev_numsols + numrep;
/*				  		printf("number of mip starts to use: %d\n",num_starts);*/
		  		
		  		prev_numsols = numsolns;
		  		
			    	if(num_starts > global_num_starts)
			    	{
/*					    		printf("reallocation: (%d)\n",__LINE__);*/
			    		global_num_starts = num_starts;
			    		global_startspace = cur_numcols*global_num_starts;
			    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
					global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
					global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
					global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
			    	}

				if(num_starts > 0)
				{
					status = CPXgetmipstarts (env_just_solve_mips, presolve_lp, &nzcnt, global_beg, global_varindices, 
							   global_values, global_effortlevel, global_startspace,
							   &surplus, 0, num_starts-1);
					if ( status ) 
				  	{
				    		printf ("(%d) Failed to get mip starts. Error code %d\n",__LINE__,status);
/*					    		return status;*/
				  	}
							   
					status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, 
						global_varindices,global_values, global_effortlevel, NULL);
					if ( status ) 
				  	{
				    		printf ("(%d) Failed to add mip starts. Error code %d\n",__LINE__,status);
/*					    		return status;*/
				  	}
			  	}
				
			    	status2 = CPXgetstat(env_just_solve_mips, presolve_lp);
			    	numsols = CPXgetsolnpoolnumsolns (env_just_solve_mips, presolve_lp);

				if(numsols >= prev_numsolns) times_to_run = numsols - prev_numsolns;
			  	else times_to_run = numsols;
			  	prev_numsolns = numsols;

				for(i=0;i<times_to_run;i++)
			    	{
			    		status = CPXgetsolnpoolx (env_just_solve_mips, presolve_lp, i, objvals, obj1_index, obj2_index);
			    		if ( status ) 
				  	{
				    		printf ("(%d) Failed to get x's. Error code %d\n",__LINE__,status);
				    		return status;
				  	}
			    		best_obj_val = objvals[0];
/*					    		printf("plot(%lf,%lf,'ko');\n",objvals[0],objvals[1]);*/
			    		insert_check = mock_insert(1,objvals[0],objvals[1],0,0,0,&tree); 
			    		if(insert_check)
					{
				      		num_added++;
				      		status = CPXgetsolnpoolx (env_just_solve_mips, presolve_lp, i, x, 0, cur_numcols-2);
				      		if ( status ) 
					  	{
					    		printf ("(%d) Failed to get x's. Error code %d\n",__LINE__,status);
					    		return status;
					  	}
		/*		      		printf("status %d\n",status);*/
						int k;
				      		for(k=0;k<cur_numcols-2;k++)
				      		{
				      			stored_x[x_rotation][k] = x[k];
					      	}
					      	x_rotation = (x_rotation + 1) % num_x_to_store;
						insert_check = 0;
				      		PSA_full(env,NULL,x,NULL,NULL);
				      	}
			    	}
			    	
			    	if(num_added)
			    	{
			    		num_inserted++;
			    		temp_interval->miss_count = 0;
			    	} 
			    	else temp_interval->miss_count++;
/*					    	{*/
/*						    	num_inserted++;*/
				    	num_h_intervals++;
				    	interval *next_interval = (struct interval*) malloc( sizeof( struct interval ) );
				    	next_interval->ub = temp_interval->ub;
				    	next_interval->lb = bound[0];
				    	next_interval->dir = temp_interval->dir;
				    	next_interval->miss_count = temp_interval->miss_count;
				    	next_interval->next = temp_interval->next;
				    	temp_interval->next = next_interval;
/*					    	}*/
			    	temp_interval->ub = bound[0];
			
				mips_solved++;
				temp_interval = next;
/*					exit(0);*/
		      	}
		      	temp_interval = first_v_interval;
		      	for(j=0;j<cur_num_v_intervals;j++)
			{	
/*						printf("interval: %lf to %lf\n",temp_interval->lb,temp_interval->ub);*/
				interval *next = NULL;
				next = temp_interval->next;
				bound[0] = (temp_interval->lb + temp_interval->ub)/2.;
				bound[1] = temp_interval->ub;
				int num_added = 0;
				double best_obj_val = 0.;
				
				status = CPXchgbds (env_just_solve_mips, presolve_lp1, 2, vars2, low, bound);
				if ( status ) {
					printf ("(%d) Failed to chg bds. Error code %d\n", __LINE__,status);
					return status;
				}
				
				num_starts = CPXgetnummipstarts (env_just_solve_mips, global_mip);
		    	
			    	if(num_starts > global_num_starts)
			    	{
/*					    		printf("reallocation: (%d)\n",__LINE__);*/
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
				if ( status ) 
			  	{
			    		printf ("(%d) Failed to get x's. Error code %d\n",__LINE__,status);
			    		return status;
			  	}
						   
				status = CPXaddmipstarts (env_just_solve_mips, presolve_lp1, num_starts, nzcnt, global_beg, global_varindices,
						   global_values, global_effortlevel, NULL);
				if ( status ) 
			  	{
			    		printf ("(%d) Failed to add mip starts. Error code %d\n",__LINE__,status);
			    		return status;
			  	}

				if(exact_mips) CPXmipopt (env_just_solve_mips, presolve_lp);
				else CPXpopulate (env_just_solve_mips, presolve_lp);
			
				numsolns = CPXgetsolnpoolnumsolns (env_just_solve_mips, presolve_lp1);
				numrep = CPXgetsolnpoolnumreplaced (env_just_solve_mips, presolve_lp1);
			
				num_starts = numsolns - prev_numsols + numrep;
/*				  		printf("number of mip starts to use: %d\n",num_starts);*/
		  		
		  		prev_numsols = numsolns;
		  		
			    	if(num_starts > global_num_starts)
			    	{
/*					    		printf("reallocation: (%d)\n",__LINE__);*/
			    		global_num_starts = num_starts;
			    		global_startspace = cur_numcols*global_num_starts;
			    		global_beg = (int *) realloc (global_beg, (global_num_starts)*sizeof(int));
					global_varindices = (int *) realloc (global_varindices, (global_startspace)*sizeof(int));
					global_values = (double *) realloc (global_values, (global_startspace)*sizeof(double));
					global_effortlevel = (int *) realloc (global_effortlevel, (global_startspace)*sizeof(int));
			    	}

				if(num_starts > 0)
				{
					status = CPXgetmipstarts (env_just_solve_mips, presolve_lp1, &nzcnt, global_beg, global_varindices, 
							   global_values, global_effortlevel, global_startspace,
							   &surplus, 0, num_starts-1);
					if ( status ) 
				  	{
				    		printf ("(%d) Failed to get mip starts. Error code %d\n",__LINE__,status);
/*					    		return status;*/
				  	}
							   
					status = CPXaddmipstarts (env_just_solve_mips, global_mip, num_starts, nzcnt, global_beg, 
						global_varindices,global_values, global_effortlevel, NULL);
					if ( status ) 
				  	{
				    		printf ("(%d) Failed to add mip starts. Error code %d\n",__LINE__,status);
/*					    		return status;*/
				  	}
			  	}
				
/*						CPXpopulate(env,presolve_lp1);*/
			    	status2 = CPXgetstat(env_just_solve_mips, presolve_lp1);
			    	numsols = CPXgetsolnpoolnumsolns (env_just_solve_mips, presolve_lp1);

				if(numsols >= prev_numsolns) times_to_run = numsols - prev_numsolns;
			  	else times_to_run = numsols;
			  	prev_numsolns = numsols;

			    	for(i=0;i<times_to_run;i++)
				{
			    		status = CPXgetsolnpoolx (env_just_solve_mips, presolve_lp1, i, objvals, obj1_index, obj2_index);
			    		if ( status ) 
				  	{
				    		printf ("(%d) Failed to get x's. Error code %d\n",__LINE__,status);
				    		return status;
				  	}
			    		best_obj_val = objvals[0];
/*					    		printf("plot(%lf,%lf,'ko');\n",objvals[0],objvals[1]);*/
			    		insert_check = mock_insert(1,objvals[0],objvals[1],0,0,0,&tree); 
			    		if(insert_check)
				      	{
					      	num_added++;
					      	status = CPXgetsolnpoolx (env_just_solve_mips, presolve_lp1, i, x, 0, cur_numcols-2);
					      	if ( status ) 
					  	{
					    		printf ("(%d) Failed to get x's. Error code %d\n",__LINE__,status);
					    		return status;
					  	}
		/*		     		printf("status %d\n",status);*/
						int k;
				      		for(k=0;k<cur_numcols-2;k++)
				      		{
				      			stored_x[x_rotation][k] = x[k];
				      		}
				      		x_rotation = (x_rotation + 1) % num_x_to_store;
				      		insert_check = 0;
					      	PSA_full(env,NULL,x,NULL,NULL);
					}
				}
			    	
			    	if(num_added)
			    	{
			    		num_inserted++;
			    		temp_interval->miss_count = 0;
			    	} 
			    	else temp_interval->miss_count++;
			    	num_v_intervals++;
			    	interval *next_interval = (struct interval*) malloc( sizeof( struct interval ) );
			    	next_interval->ub = temp_interval->ub;
			    	next_interval->lb = bound[0];
			    	next_interval->dir = temp_interval->dir;
			    	next_interval->miss_count = temp_interval->miss_count;
			    	next_interval->next = temp_interval->next;
			    	temp_interval->next = next_interval;
			    	temp_interval->ub = bound[0];
			
				mips_solved++;
				temp_interval = next;
		      	}
			if(pure_binary == 1) run_it += binary_heuristic(env,lp1,lp2,indices);
			else run_it += mixed_heuristic(env);
		}
		status = CPXfreeprob (env_just_solve_mips, &presolve_lp1);
    		if ( status ) 
    		{
      			printf ("CPXfreeprob failed, error code %d.\n",status);
    		}
	}
	
	slope_val *temp = first_slope;
	slope_val *next = NULL;
	if(first_slope && first_slope->next) next = first_slope->next;
/*			printf("num slopes: %d\n", num_slopes);*/
	if(first_slope != NULL) for(j=0;j<num_slopes;j++)
	{
/*				printf("j: %d\n",j);*/
		if(temp != NULL) free(temp);
		temp = next;
		if(j != num_slopes-1) next = next->next;
	}
	
	if(intervals) 
	{
/*				printf("num intervals: %d\n",num_intervals);*/
		interval *temp = first_h_interval;
		interval *next = NULL;
		if(first_h_interval && first_h_interval->next) next = first_h_interval->next;
/*			printf("num slopes: %d\n", num_slopes);*/
		if(first_h_interval != NULL) for(j=0;j<num_h_intervals;j++)
		{
			if(temp != NULL) free(temp);
			temp = next;
			if(j != num_h_intervals-1) next = next->next;
		}
		if(first_v_interval && first_v_interval->next) next = first_v_interval->next;
/*			printf("num slopes: %d\n", num_slopes);*/
		if(first_v_interval != NULL) for(j=0;j<num_v_intervals;j++)
		{
			if(temp != NULL) free(temp);
			temp = next;
			if(j != num_v_intervals-1) next = next->next;
		}
	}
	
  	END_PREPROCESSING3:
	
	finish_presolve = clock();
	duration_presolve = (double)(finish_presolve - start_presolve) / CLOCKS_PER_SEC;
  	duration_BB = (double)(finish_presolve - start_BB) / CLOCKS_PER_SEC;
  	printf("*******************\n");
  	printf("Prepopulate time: (Alone) %lf\t (Cumulative) %lf\nIterations: %d, MIPs solved: %d Num sols added from heur: %d\n",
  									duration_presolve,duration_BB,presolve_iteration,mips_solved,run_it);
  	if(!suppress_file_output) fprintf(bb_results,"%lf\t",duration_BB);
  	printf("*******************\n");
  	
  	if(x) free_and_null ((char **) &x);
  	status = CPXfreeprob (env_just_solve_mips, &presolve_lp);
	if ( status ) 
	{
		printf ("CPXfreeprob failed, error code %d.\n",status);
	}
	return 0;
}

int probing2(clock_t phase2_intermediate_time, double time_so_far, clock_t phase2_start_time, int fix_it2, double *lb_pre, double *ub_pre, char *sym, CPXENVptr env,
		CPXLPptr presolve2_lp, CPXLPptr presolve2_mip, int *indices, double *SE_objectives, int insert_check, double *NW_objectives, double *WS_objectives,
		double temp_x, double slope, CPXLPptr lp1, CPXLPptr lp2, int num_bounds_reduced, int cur_numcols, char *xctype, int cur_numrows,
		double *obj_coef1, double *obj_coef2, double *weighted_coefs_, double temp_y)
{
/*	printf("slope: %lf\n",slope);*/
	int i = 0, k = 0, status = 0;
	double ran = (double) rand() / ( (double) RAND_MAX);
	int ten_percent = (int) (total_num_integer*(percentage_of_integer_variables_to_check_for_probing_during_presolve/100.));
	int starting_index = (int) (ran*((double) total_num_integer - (double) ten_percent));
	if(starting_index+ten_percent >= total_num_integer) starting_index = 0;

/*	printf("total num int: %d\n",total_num_integer);*/

	if(!limit_bd_reduction)
	{
		starting_index = 0;
		ten_percent = total_num_integer;
	}

	for(i=starting_index;i<starting_index+ten_percent;i++)
	{
/*		printf("i: %d\n",i);*/
		phase2_intermediate_time = clock();
		time_so_far = (double)(phase2_intermediate_time - phase2_start_time) / CLOCKS_PER_SEC;
	
		if(time_so_far > max_time_phase_2) 
		{
			printf("Exceeded maximum time of %lf in phase 2. Exitting!\n",max_time_phase_2);
			goto END_OF_PHASE2;
		}
	
		k = integer_indices[i];
		fix_it2 = 0;
		if(lb_pre[k] != ub_pre[k])
		{
			REPEAT_IT:
			fix_it2 = 0;
			/**************fix to lower bounds*************/
/*				printf("fixing x%d to %d:\n",k,(int)lb_pre[k]);*/
			sym[0] = 'U';
    			status = CPXchgbds (env, presolve2_lp, 1, &k, sym, &lb_pre[k]);
    			status = CPXchgbds (env, presolve2_mip, 1, &k, sym, &lb_pre[k]);
    			
/*	    			printf("-*_*_*_*_*_*_*_*_*_*_*_*-\n");*/
/*				PSA_all(env,presolve2_lp);*/
/*				printf("-*_*_*_*_*_*_*_*_*_*_*_*-\n");*/
    			
			chg_coefs_(env,presolve2_lp,indices,-.01);
			status = CPXlpopt (env, presolve2_lp);
			int lpstat = CPXgetstat (env, presolve2_lp);
/*	  			printf("the status of the solve: %d\n",lpstat);*/
  			if(lpstat == 3)
  			{
  				fix_it2 = 1;
  				goto FIXING1;
  			}
			status = CPXgetx (env, presolve2_lp, SE_objectives, obj1_index, obj2_index);
  			if(status)
  			{
				printf("Failed to get x-values from CPLEX. Status: %d Line: %d\n",status,__LINE__);
				exit(0);
  				return status;
  			}
  			insert_check = mock_insert(1,SE_objectives[0],SE_objectives[1],0,0,0,&tree);
  			if(!insert_check)
  			{
/*	  				printf("SE pt dominated when fixing x%d to its lb\n",k);*/
/*		  			printf("plot(%lf,%lf,'o');\n",SE_objectives[0],SE_objectives[1]);*/
	  			chg_coefs_(env,presolve2_lp,indices,-100.);
				status = CPXlpopt (env, presolve2_lp);
				status = CPXgetx (env, presolve2_lp, NW_objectives, obj1_index, obj2_index);
	  			if(status)
	  			{
					printf("Failed to get x-values from CPLEX. Line: %d\n",__LINE__);
	  				return status;
	  			}
	  			insert_check = mock_insert(1,NW_objectives[0],NW_objectives[1],0,0,0,&tree);
  				if(!insert_check)
  				{
/*	  					printf("NW pt also dominated when fixing x%d to its lb\n",k);*/
/*			  			printf("plot(%lf,%lf,'o');\n",NW_objectives[0],NW_objectives[1]);*/
		  			insert_check = mock_insert(1,NW_objectives[0],SE_objectives[1],0,0,0,&tree);
		  			if(!insert_check)
		  			{
/*			  				printf("local ideal pt dominated when fixing x%d to its lb\n",k);*/
/*			  				printf("plot(%lf,%lf,'go');\n",NW_objectives[0],SE_objectives[1]);*/
		  				fix_it2 = 1;
		  			}
		  			else
		  			{
		  				chg_coefs_(env,presolve2_lp,indices,-1.);
		  				status = CPXlpopt (env, presolve2_lp);
/*				  			printf("the status of the WS mip solve: %d\n",lpstat);*/
						status = CPXgetx (env, presolve2_lp, WS_objectives, obj1_index, obj2_index);
/*						printf("plot(%lf,%lf,'o');\n",WS_objectives[0],WS_objectives[1]);*/
						insert_check = mock_insert(1,WS_objectives[0],WS_objectives[1],0,0,0,&tree);
						if(!insert_check)
						{
							temp_x = 1/slope*(SE_objectives[1]-WS_objectives[1])+WS_objectives[0];
							temp_y = slope*(NW_objectives[0]-WS_objectives[0])+WS_objectives[1];
/*								printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",temp_x,NW_objectives[0],SE_objectives[1],temp_y);*/
							insert_check = mock_insert(2,NW_objectives[0],temp_y,temp_x,SE_objectives[1],slope,&tree);
							if(!insert_check)
							{
/*									printf("Changing lb of x%d to %d!\n",k,(int)(lb_pre[k]+1));*/
								fix_it2 = 1;
							}
						}
		  			}
		  		}
  			}
  			FIXING1:
  			sym[0] = 'U';
    			status = CPXchgbds (env, presolve2_lp, 1, &k, sym, &ub_pre[k]);
    			status = CPXchgbds (env, presolve2_mip, 1, &k, sym, &ub_pre[k]);
    			if(fix_it2 == 1)
    			{
/*	    				printf("Changing lb of x%d to %d!\n",k,(int)(lb_pre[k]+1));*/
    				sym[0] = 'L';
    				double bound[1] = {lb_pre[k]+1.};
    				status = CPXchgbds (env, lp1, 1, &k, sym, bound);
    				status = CPXchgbds (env, lp2, 1, &k, sym, bound);
    				status = CPXchgbds (env_just_solve_mips, global_mip, 1, &k, sym, bound);
    				status = CPXchgbds (env, presolve2_lp, 1, &k, sym, bound);
    				status = CPXchgbds (env, presolve2_mip, 1, &k, sym, bound);
    				num_bounds_reduced++;
    				status = CPXgetlb (env, lp1, lb_pre, 0, cur_numcols-1);
    				if(lb_pre[k] + 1 < ub_pre[k])
    				{
/*	    					printf("x%d has had its lower bound increased by 1. Attempting to increase further\n",k);*/
    					fix_it2 = 0;
    					goto REPEAT_IT;
    				}
    			}
    			
    			/**************fix to upper bounds*************/
    			if(lb_pre[k] != ub_pre[k] && (!fix_it2 || xctype[k] != 'B'))
    			{
    				REPEAT_IT2:
    				fix_it2 = 0;
/*	    				printf("fixing x%d to %d:\n",k,(int)ub_pre[k]);*/
				sym[0] = 'L';
	    			status = CPXchgbds (env, presolve2_lp, 1, &k, sym, &ub_pre[k]);
	    			status = CPXchgbds (env, presolve2_mip, 1, &k, sym, &ub_pre[k]);
	    			
/*		    			printf("-*_*_*_*_*_*_*_*_*_*_*_*-\n");*/
/*					PSA_all(env,presolve2_lp);*/
/*					printf("-*_*_*_*_*_*_*_*_*_*_*_*-\n");*/
	    			
				chg_coefs_(env,presolve2_lp,indices,-.01);
				status = CPXlpopt (env, presolve2_lp);
				lpstat = CPXgetstat (env, presolve2_lp);
	  			if(lpstat == 3) 
	  			{
	  				fix_it2 = 1;
	  				goto FIXING2;
	  			}
				status = CPXgetx (env, presolve2_lp, SE_objectives, obj1_index, obj2_index);
	  			if(status)
	  			{
					printf("Failed to get x-values from CPLEX. Line: %d\n",__LINE__);
					exit(0);
	  				return status;
	  			}
	  			insert_check = mock_insert(1,SE_objectives[0],SE_objectives[1],0,0,0,&tree);
	  			if(!insert_check)
	  			{
/*		  				printf("SE pt dominated when fixing x%d to its lb\n",k);*/
/*			  			printf("plot(%lf,%lf,'o');\n",SE_objectives[0],SE_objectives[1]);*/
		  			chg_coefs_(env,presolve2_lp,indices,-100.);
					status = CPXlpopt (env, presolve2_lp);
					status = CPXgetx (env, presolve2_lp, NW_objectives, obj1_index, obj2_index);
		  			if(status)
		  			{
						printf("Failed to get x-values from CPLEX. Line: %d\n",__LINE__);
		  				return status;
		  			}
		  			insert_check = mock_insert(1,NW_objectives[0],NW_objectives[1],0,0,0,&tree);
	  				if(!insert_check)
	  				{
/*		  					printf("NW pt also dominated when fixing x%d to its lb\n",k);*/
/*				  			printf("plot(%lf,%lf,'o');\n",NW_objectives[0],NW_objectives[1]);*/
			  			insert_check = mock_insert(1,NW_objectives[0],SE_objectives[1],0,0,0,&tree);
			  			if(!insert_check)
			  			{
/*				  				printf("local ideal pt dominated when fixing x%d to its lb\n",k);*/
/*				  				printf("plot(%lf,%lf,'go');\n",NW_objectives[0],SE_objectives[1]);*/
			  				fix_it2 = 1;
			  			}
			  			else
			  			{
			  				chg_coefs_(env,presolve2_lp,indices,-1.);
			  				status = CPXlpopt (env, presolve2_lp);
/*					  			printf("the status of the WS mip solve: %d\n",lpstat);*/
							status = CPXgetx (env, presolve2_lp, WS_objectives, obj1_index, obj2_index);
/*							printf("plot(%lf,%lf,'o');\n",WS_objectives[0],WS_objectives[1]);*/
							insert_check = mock_insert(1,WS_objectives[0],WS_objectives[1],0,0,0,&tree);
							if(!insert_check)
							{
								temp_x = 1/slope*(SE_objectives[1]-WS_objectives[1])+WS_objectives[0];
								temp_y = slope*(NW_objectives[0]-WS_objectives[0])+WS_objectives[1];
/*									printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",temp_x,NW_objectives[0],*/
/*										SE_objectives[1],temp_y);*/
								insert_check = mock_insert(2,NW_objectives[0],temp_y,temp_x,SE_objectives[1],slope,&tree);
								if(!insert_check)
								{
/*										printf("Changing ub of x%d to %d!\n",k,(int)(ub_pre[k]-1));*/
									fix_it2 = 1;
								}
							}
			  			}
			  		}
	  			}
	  			FIXING2:
	  			sym[0] = 'L';
	    			status = CPXchgbds (env, presolve2_lp, 1, &k, sym, &lb_pre[k]);
	    			status = CPXchgbds (env, presolve2_mip, 1, &k, sym, &lb_pre[k]);
	    			if(fix_it2 == 1)
	    			{
/*		    				printf("Changing ub of x%d to %d!\n",k,(int)(ub_pre[k]-1));*/
	    				sym[0] = 'U';
	    				double bound[1] = {ub_pre[k]-1.};
	    				status = CPXchgbds (env, lp1, 1, &k, sym, bound);
	    				status = CPXchgbds (env, lp2, 1, &k, sym, bound);
	    				status = CPXchgbds (env_just_solve_mips, global_mip, 1, &k, sym, bound);
	    				status = CPXchgbds (env, presolve2_lp, 1, &k, sym, bound);
    					status = CPXchgbds (env, presolve2_mip, 1, &k, sym, bound);
    					num_bounds_reduced++;
    					status = CPXgetub (env, lp1, ub_pre, 0, cur_numcols-1);
    					if(ub_pre[k] - 1 > lb_pre[k])
	    				{
/*		    					printf("x%d has had its upper bound reduced by 1. Attempting to reduce further\n",k);*/
	    					fix_it2 = 0;
	    					goto REPEAT_IT2;
	    				}
	    			}
    			}
		}
	}
	END_OF_PHASE2:

	printf("number of times a bound was reduced during phase 2: %d\n",num_bounds_reduced);

	return 0;
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
