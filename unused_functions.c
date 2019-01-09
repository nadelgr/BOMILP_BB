/*********************************************************************************************************************** 

	This function is currently unused, but was designed for adding a globally valid cut along
	the level curve in the objective space associated with the best known dual bound after solving
	any weighted sum MIPs during preprocessing.
	
***********************************************************************************************************************/

void add_mip_cut(CPXCENVptr env, CPXLPptr prob1, CPXLPptr prob2, CPXLPptr prob3, double slp, double bound, 
			double *ob1coef, double *ob2coef, double *wtd_coefs, int *indices)
{
	int status,i;
	double rhs = 0.;
	double val = -1.;
	double lb = -CPX_INFBOUND;
	double ub = bound + .001;
	char sense[1] = {'E'};
	int rmatbeg[1] = {0};
	
/*	printf("slope: %lf\n",slp);*/
/*	*/
/*	for(i=0;i<cur_numcols;i++) */
/*	{*/
/*		printf("wtd%d: %lf, ob1_%d: %lf, obj2_%d: %lf\n",i,wtd_coefs[i],i,ob1coef[i],i,ob2coef[i]);*/
/*		wtd_coefs[i] = (-1./slp)*ob2coef[i] + ob1coef[i];*/
/*		printf("new wtd%d: %lf\n",i,wtd_coefs[i]);*/
/*	}*/
	
	status = CPXaddrows (env, prob1, 0, 1, cur_numcols, &rhs,
                      sense, rmatbeg, indices, wtd_coefs, NULL, NULL);
        status = CPXaddrows (env, prob2, 0, 1, cur_numcols, &rhs,
                      sense, rmatbeg, indices, wtd_coefs, NULL, NULL);
        status = CPXaddrows (env, prob3, 0, 1, cur_numcols, &rhs,
                      sense, rmatbeg, indices, wtd_coefs, NULL, NULL);
        
        status = CPXaddcols (env, prob1, 1, 1, NULL, rmatbeg,
                      &cur_numrows, &val, &lb, &ub, NULL);
        status = CPXaddcols (env, prob2, 1, 1, NULL, rmatbeg,
                      &cur_numrows, &val, &lb, &ub, NULL);
        status = CPXaddcols (env, prob3, 1, 1, NULL, rmatbeg,
                      &cur_numrows, &val, &lb, &ub, NULL);
                      
/*        status = CPXwriteprob (env, prob2, "myprob2.lp", "LP");*/
/*        exit(0);*/
        
        cur_numrows++;
        cur_numcols++;
        obj1_index = cur_numcols;
}

/*********************************************************************************************************************** 

	This function is currently unused, but was designed for adding no good constraints associated 
	with discovered solutions for binary problems.
	
***********************************************************************************************************************/

void add_binary_cut(CPXCENVptr env, CPXLPptr prob1, CPXLPptr prob2, double *x, int *indices, double *coefs)
{
	int status,i,k;
	double rhs = -1.;
	char sense[1] = {'L'};
	int rmatbeg[1] = {0};
	
	for(i=0;i<total_num_integer;i++)
	{
		k = integer_indices[i];
		if(x[k] == 1)
		{
			coefs[i] = 1.;
			rhs++;
		}
		else coefs[i] = -1.;
	}
	
	status = CPXaddrows (env, prob1, 0, 1, total_num_integer, &rhs,
                      sense, rmatbeg, integer_indices, coefs, NULL, NULL);
        if(prob2) status = CPXaddrows (env, prob2, 0, 1, total_num_integer, &rhs,
                      sense, rmatbeg, integer_indices, coefs, NULL, NULL);
        cur_numrows++;
        
/*        for(i=0;i<cur_numcols;i++) printf("x%d: %lf\n",i,x[i]);*/
/*        status = CPXwriteprob (env, prob1, "myprob1.lp", "LP");*/
/*        status = CPXwriteprob (env, prob2, "myprob2.lp", "LP");*/
/*        exit(0);*/
}


/*********************************************************************************************************************** 

	The following two functions are unused, but were originally designed for using equality 
	constraints to substitute variables out of the model. After thinking about it, though, I
	realized that it isn't helpful when variables have bounds.
	
***********************************************************************************************************************/

double** make_two_d_array(int arraySizeX, int arraySizeY) 
{
	int i;
	double** theArray;
	theArray = (double**) malloc(arraySizeX*sizeof(double*));
	for (i = 0; i < arraySizeX; i++) theArray[i] = (double*) malloc(arraySizeY*sizeof(double));
   	return theArray;
} 

void remove_equalities(CPXCENVptr env, CPXLPptr prob1, CPXLPptr prob2, char *xctype, double *obj1_coef, double *obj2_coef, int *indices)
{
	double **A = make_two_d_array(cur_numrows,cur_numcols);
	
	int nzcnt,surplus;
	int *cmatbeg = (int*) malloc(cur_numcols*sizeof(int));
	int *cmatind = (int*) malloc(cur_numcols*cur_numrows*sizeof(int));
	double *cmatval = (double*) malloc(cur_numcols*cur_numrows*sizeof(double));
	double *rhs = (double *) malloc (cur_numrows*sizeof(double));
	char *sense = (char *) malloc (cur_numrows*sizeof(char));
	int *rows_to_delete = (int*) malloc((cur_numrows+1)*sizeof(int));
	int *cols_to_delete = (int*) malloc((cur_numcols+1)*sizeof(int));
	double *new_vals = (double*) malloc(cur_numcols*sizeof(double));
	int num_rows_to_delete = 0;
	int num_cols_to_delete = 0;
	rows_to_delete[num_rows_to_delete] = -1;
	cols_to_delete[num_cols_to_delete] = -1;
	
	status = CPXgetcols (env, prob1, &nzcnt, cmatbeg, cmatind,
                      cmatval, cur_numcols*cur_numrows, &surplus, 0,
                      cur_numcols-1);
                      
        status = CPXgetrhs (env, prob1, rhs, 0, cur_numrows-1);
        
        status = CPXgetsense (env, prob1, sense, 0, cur_numrows-1);
        
        int i,j;
        int k = 0;
        
        for(i=0;i<cur_numrows;i++) for(j=0;j<cur_numcols;j++) A[i][j] = 0;
        
        for(i=0;i<nzcnt;i++)
        {
        	if(i == cmatbeg[k+1]) k++;
        	A[cmatind[i]][k] = cmatval[i];
        }
        k = 0;
        
/*        for(i=0;i<cur_numrows;i++) */
/*        {*/
/*        	for(j=0;j<cur_numcols;j++) printf("%lf ",A[i][j]);*/
/*        	printf("\n");*/
/*        }*/
        
        int l;
        for(i=0;i<cur_numrows;i++)
        {
/*        	printf("considering row %d\n",i);*/
        	k = 0;
        	if(sense[i] == 'E')
        	{
/*        		printf("the row is an equality constraint\n");*/
        		while( ( (A[i][k] < 0.00001 && A[i][k] > -0.000001) || xctype[k] == 'I' || xctype[k] == 'B') && k < cur_numcols) k++;
        		if(k < cur_numcols)
        		{
/*        			printf("solving for x%d\n",k);*/
        			if(obj1_coef[k] != 0) 
        			{
        				for(j=0;j<cur_numcols;j++)
        				{
        					new_vals[j] = obj1_coef[j] - obj1_coef[k]*A[i][j]/A[i][k];
        				}
        				obj1_extra_val += obj1_coef[k]*rhs[i]/A[i][k];
        				for(j=0;j<cur_numcols;j++)
        				{
        					obj1_coef[j] = new_vals[j];
/*        					printf("obj1_coef%d: %lf\n",j,obj1_coef[j]);*/
        				}
        			}
        			if(obj2_coef[k] != 0) 
        			{
        				for(j=0;j<cur_numcols;j++)
        				{
        					new_vals[j] = obj2_coef[j] - obj2_coef[k]*A[i][j]/A[i][k];
        				}
        				obj2_extra_val += obj2_coef[k]*rhs[i]/A[i][k];
        				for(j=0;j<cur_numcols;j++)
        				{
        					obj2_coef[j] = new_vals[j];
/*        					printf("obj2_coef%d: %lf\n",j,obj2_coef[j]);*/
        				}
        			}
        			
        			for(l=0;l<cur_numrows;l++)
        			{
/*        				printf("changing row coefs\n");*/
        				if(l != i)
					{
						for(j=0;j<cur_numcols;j++)
						{
							new_vals[j] = A[l][j] - A[l][k]*A[i][j]/A[i][k];
						}
						rhs[l] -= A[l][k]*rhs[i]/A[i][k];
						A[l][k] = 0;
						for(j=0;j<cur_numcols;j++)
						{
							A[l][j] = new_vals[j];
/*							printf("A%d,%d: %lf\n",l,j,A[l][j]);*/
						}
					}
					else
					{
						for(j=0;j<cur_numcols;j++)
						{
							new_vals[j] = A[i][j]/A[i][k];
						}
						rhs[l] = rhs[i]/A[i][k];
						A[l][k] = 0;
						for(j=0;j<cur_numcols;j++)
						{
							A[l][j] = new_vals[j];
/*							printf("A%d,%d: %lf\n",l,j,A[l][j]);*/
						}
					}
				}
				sense[i] = 'L';
        			rows_to_delete[num_rows_to_delete] = i;
        			rows_to_delete[num_rows_to_delete+1] = -1;
				cols_to_delete[num_cols_to_delete] = k;
				cols_to_delete[num_cols_to_delete+1] = -1;
				num_rows_to_delete++;
				num_cols_to_delete++;
        		}
        	}
        }
        
        status = CPXchgobj (env, prob1, cur_numcols, indices, obj1_coef);
        status = CPXchgobj (env, prob2, cur_numcols, indices, obj2_coef);
        status = CPXchgrhs (env, prob1, cur_numrows, indices, rhs);
        status = CPXchgrhs (env, prob2, cur_numrows, indices, rhs);
        
/*        qsort(rows_to_delete, num_rows_to_delete, sizeof(int), comparison_function2);*/
        qsort(cols_to_delete, num_cols_to_delete, sizeof(int), comparison_function2);
        
        k = 0;
        for(i=0;i<cur_numrows;i++)
        {
/*        	if(i != rows_to_delete[k])*/
/*        	{*/
        		l = 0;
			for(j=0;j<cur_numcols;j++)
			{
				if(j != cols_to_delete[l])
				{
					status = CPXchgcoef (env, prob1, i, j, A[i][j]);
					status = CPXchgcoef (env, prob2, i, j, A[i][j]);
				}
				else l++;
			}
/*		}*/
/*		else k++;*/
        }
        
        for(i=0;i<num_cols_to_delete;i++) 
        {
/*        	status = CPXdelrows (env, prob1, rows_to_delete[i], rows_to_delete[i]);*/
/*        	status = CPXdelrows (env, prob2, rows_to_delete[i], rows_to_delete[i]);*/
        	status = CPXdelcols (env, prob1, cols_to_delete[i], cols_to_delete[i]);
        	status = CPXdelcols (env, prob2, cols_to_delete[i], cols_to_delete[i]);
        }
        
        status = CPXchgsense (env, prob1, num_rows_to_delete, rows_to_delete, sense);
        status = CPXchgsense (env, prob2, num_rows_to_delete, rows_to_delete, sense);
        
        free(cmatbeg);
	free(cmatind);
	free(cmatval);
	free(rhs);
	free(sense);
	free(rows_to_delete);
	free(cols_to_delete);
	free(new_vals);
	for(i=0;i<cur_numrows;i++) free(A[i]);
	free(A);	
}
